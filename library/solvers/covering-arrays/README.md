# covering-arrays — generic CA(N;t,k,v) search (simulated annealing + restarts)

**Class:** *search-then-certify.* Find an object — an `N x k` array over `0..v-1` such that
every `t`-subset of columns shows all `v^t` tuples — and certify it with a short independent
checker. Hard to find, trivial to verify, so a local search plus a frozen checker is the whole
toolchain. Same shape as code/design searches: the move engine is incremental, the gate is
somebody else's program.

**Solver:** `ca_sa.c` (one file, `-lm` only, 263 lines).
`sha256 56078a8966a2fa824803251673b09ef61c1825e8b441063740b3a94ac8b754be`.
Build: `make` (⇒ `cc -O2 -Wall -Wextra -std=c11 -march=native -o ca_sa ca_sa.c -lm`).

```
./ca_sa t k v N [--iters M] [--restarts R] [--seed S] [--init FILE] [--out FILE]
                 [--row0zero 0|1] [--T0 X] [--T1 X] [--quiet]
```
Exit 0 iff the *solver's* uncovered-tuple count reached 0. **The exit status is not a
verdict** — certify the output with the frozen field checker:
`python3 ../../../software-testing/covering-arrays/field/checker.py FILE t k v N`.

**Checker:** not copied here on purpose — the field has exactly one frozen copy,
`software-testing/covering-arrays/field/checker.py`,
`sha256 97745376c1fea9efe8cd2b78e9d447e13609e78c58a8bf8d27f80f12b635551c`. A second copy
would be a second thing to drift.

**Structure used** (and what each one buys):
- Cost = number of uncovered `(t-subset, tuple)` pairs; incremental bookkeeping touches only
  `C(k-1,t-1)` counters per move, so one move is O(`C(k-1,t-1)·t`) regardless of `N`.
- Row 0 pinned to all zeros: per-column symbol permutations are automorphisms of the covering
  property, so this is WLOG and removes a factor `v^k` (for the calibrated target, 3^9 = 19,683).
- Balanced init (each symbol `floor(N/v)` or `ceil(N/v)` times per column) satisfies the
  necessary condition `r_{j,s} >= v` and starts with far fewer uncovered tuples (~32-55 instead
  of ~200 for a uniform-random start at 13x9).
- `S_k` and `S_N` are automorphisms; no canonical form imposed on them, so the output is not
  normalised — no loss, just no dedup.
- Move = pick a random cell (rows 1..N-1), a random different symbol; accept if `d <= 0` or
  `exp(-d/T) > u`; geometric cooling from `T0=1.6` to `T1=0.02` over `--iters` moves.
- `--init FILE` re-scores an existing array (`--iters 0`) — this is how the cost function was
  cross-checked against the checker, not just assumed sign-consistent.

**Speed measured on this machine** (nothing parallel, `-O2 -march=native`):

| instance | density | budget | wall | result |
|---|---|---|---|---|
| CA(9;2,4,3) perfect OA | 1.000 | 300k x 20, seed 7 | 5 ms | cost 0, checker PASS |
| CA(11;2,5,3) | 0.818 | 300k x 20, seed 7 | 8 ms | cost 0, checker PASS |
| CA(12;2,7,3) | 0.750 | 300k x 20, seed 7 | 109 ms | cost 0, checker PASS |
| CA(13;2,9,3) | 0.692 | 200k x 20, seed 13 | 303 ms | cost 2 — **stall** |
| CA(13;2,9,3) | 0.692 | 200k x 20, seed 1 | 287 ms | cost 0, checker PASS |
| CA(12;2,9,3) negative control | 0.750 | 300k x 10, seed 11 | 217 ms | cost 11, checker FAIL |

≈1.3e7 moves/s at k=9,t=2. Cost table capped at `C(k,t)·v^t <= 8e7` int32 counters, which is
the practical limit: t=2, v=3 tolerates k up to ~7000, t=3,v=3 k up to ~40, t=4 v=2 k up to ~15.
For larger `t·k·v` the counter table is the wall (the program refuses rather than thrashing).

**Protocol that works** (the useful lesson): at a fixed move budget the outcome is strongly
seed-dependent — CA(13;2,9,3) ends at cost 2 for seed 13 and at cost 0 for seed 1 in the same
4.0e6 moves. So *sweep seeds at a fixed small budget* rather than raising budget on one seed;
`selftest.sh` shows the shape of that loop.

**Validation on record:** `selftest.sh` (below) rebuilds from source and re-certifies the three
tighter known instances with the frozen checker. Cost-function cross-check: `--init` on the
checker's good OA ⇒ 0; on the doctored OA (one row replaced) ⇒ exactly 6, one lost tuple per
column pair, matching the checker's FAIL message. Negative control above shows the machine does
not manufacture a solution below a known size.

**Used by:** [[cal-ca-2-3-k9-n13]] (calibration, PASS by the frozen checker — artifact in
`library/seeds/covering-arrays/covering-array_N13_t2_k9_v3.txt`).

## row_sa.py — Python row-level SA (added 2026-11-22, verify 1 of [[cal-ca-3-3-k6-n33]])

`row_sa.py`, `sha256 afa0e462f7af4e679c277636ca91c85850f40f3468012cdeff1b11748fbc1fb2`
(verbatim copy of the problem-folder `solver.py`; stdlib only). Two engines: cell-level `CA`
and row-level `RowSearch`.

```
python3 row_sa.py t k v N [--method cell|rowsa] [--seed S] [--iters M] [--time SEC]
                    [--restarts R] [--t0 X] [--t1 X] [--out FILE]
```

When to prefer it over `ca_sa.c`: **t=3 with small v^k** (row moves replace a whole row, i.e.
all C(k,t) of its elements at once, which breaks the strength-3 coupling that makes cell moves
plateau). Measured: CA(33;3,6,3) cost 136 → 0 in 1.2 s at seed 1 (cell-level SA plateaued at
10–25/540 over 4 restarts / 25 s on the same instance). Same certification rule as ca_sa.c:
solver output is not a verdict — run the frozen field checker. Provenance: builds its v^k
candidates from `itertools.product`, opens files for write only, reads nothing.

**Used by:** [[cal-ca-3-3-k6-n33]] (calibration, PASS by the frozen checker — artifact in the
problem folder, `array_seed1.txt`).

## rowcover.py — row-level max-coverage + AC-preserving tabu (added 2026-11-22, attempt 1 of [[cal-ca-3-2-k12-n15]])

`rowcover.py`, `sha256 3d0604680dde0a3644490cd7e10586d20da3c7b127d1afcce09070e2d3c58a27`
(verbatim copy of the problem-folder `solver/rowcover.py`; stdlib only, single process).

```
python3 rowcover.py --t T --k K --v V --N N [--seed S] [--restarts R] [--iters I]
                    [--greedy-runs G] [--tenure L] [--slots S] [--full-u U]
                    [--lns-stall K] [--lns-r R] [--time SEC] [--out FILE] [--trace FILE]
```

Third engine for the same class, and the strongest of the three at *row* granularity:

- **Element convention:** universe element = `(triple_index, tuple_code)`, `tuple_code` is the
  `t`-digit value with the first column of the triple most significant, element id =
  `triple_index * v**t + tuple_code`. A candidate row `code` has digit `j = (code // v**j) % v`;
  it covers `C(k,t)` elements, one per triple. Cost = `popcount(U)` = uncovered elements.
- **Best-of-G max-coverage greedy** (randomised tie-break among the top 4 rows) for the start.
- **Moves are whole-row replacements on a *set* of N rows.** Exact incremental `U` / `C1`
  (covered-exactly-once) bookkeeping:
  `delta(b) = popcount(a & C1) - popcount((U | (a & C1)) & b)`, all integer; a full recount
  every 200 moves aborts on drift.
- **The trick that carries it:** for the victim row `a`, `a & C1` are the constraints a swap can
  destroy. `Subcubes.codes_for` enumerates the subcube of rows satisfying all of them, so any
  candidate from it gives `delta <= 0` — the neighbourhood is monotone. Plus deficit-targeted
  subcubes (rows fixing 1–2 uncovered elements), an exhaustive single-swap scan when
  `|U| <= --full-u`, and ruin-and-recreate LNS (`2..--lns-r` rows removed, greedily re-inserted).
- `Subcubes.codes_for` is brute-force verified in the problem folder
  (`probes/test_subcubes.py`, `TOTAL BAD: 0` over (t,k,v) = (3,6,2), (3,5,3), (2,7,2), (3,9,2), (4,7,2)).

**Measured (this machine, single process):**

| instance | density | budget | wall | result |
|---|---|---|---|---|
| CA(12;3,11,2) | 0.667 | seed 7, 6 restarts x 4000 | 0.19 s | cost 0 restart 0, checker PASS |
| CA(15;3,12,2) | 0.533 | seed 1, 20 restarts x 3000 | 12.1 s | cost 0 restart 0, checker PASS |
| CA(16;3,14,2) | 0.500 | seed 21, 40 restarts x 400 | 123 s | **best 6/2912, NOT solved** |

**Limits, measured not guessed:** the k=14/N=16 instance is maximally tight (16·C(14,3) = 5824 =
2 x 2912, so the mean multiplicity per triple-pattern is exactly 2); single row swaps plateau at
6–10 uncovered there over 40 restarts and 9 parameter settings. For t=3 binary with k >= 14 use
`ca_sa.c` (cell SA) and keep this for k <= 13. Solver output is not a verdict — run the frozen
field checker.

**Used by:** [[cal-ca-3-2-k12-n15]] (calibration; array in the problem folder, frozen checker PASS).

## Sequence covering arrays: the antipodal-pair quotient (`pairq`, `pairdfs`, `pairq_model.py`)

**Class:** *sequence covering arrays* SCA(N;t,v) — N permutations of v symbols so that every
ordered t-tuple of distinct symbols is a subsequence of some row. From [[scan-4-7-26]].

- `pairq_model.py` builds, independently of the C code, the reversal classes of S_v (pairs
  {π, rev π}), the 24 orders of each t-subset, and the 12 reversal classes of those orders.
  `python3 pairq_model.py build` runs its own checks (spread ≥ C(v,t)(N−t!), U ≤ spread −
  floor, class counts 210 each) and writes `table.bin`.
- `pairq.c` (gcc -O2 -static, ~10 KB) searches 13 pairs = 26 rows for t=4,v=7 with the exact
  cost **F = Σ_{i<j} A(i,j) = Σ_T c_T**, where A(i,j) = #{T : the two classes agree on T} and
  c_T is that T's collision count. **F ≥ C(v,t) always** (N/2 pairs into t!/2 classes) and
  **F = C(v,t) ⟺ PASS inside the reversal-closed space**; the tabu move is a whole-pair swap,
  so it is not a longer run of the single-row family. Keeps a cost vector so one move is three
  contiguous scans: **29,868 moves/s measured** (vs 305/s for the python single-row version).
- `pairdfs.c` is the exact side: c_T is monotone under adding a class and c_T ≥ 1, so every
  prefix of a solution satisfies c_T ≤ 1 — a very tight DFS pruning condition.
- `cross_check.py` verifies the C tables against the python model (up to per-t relabelling) and
  that F = floor ⟺ all t! orders present; `emit_best.py` turns a class list into rows and
  re-derives spread from the definition.
- **Measured limit (t=4,v=7,N=26):** tabu best F = 69 (floor 35); CP-SAT on the 420-requirement
  covering model UNKNOWN in 240 s; monotone-prefix DFS max depth 9 of 13 over 82 restarts. All
  three agree the reversal-closed space is rigid; the free space is at least as good (U=65 vs 68
  at the same spread 138). Solver output is not a verdict — run `field/sca_check.py`.

## coleq — column extension as a panchromatic 3-colouring (t=3, v=3), 2026-11-23

`coleq.c` (static, std-only) adds one column to a verified `CA(N;3,K,3)` by 3-colouring the N rows so
every pair-class `(a,b,u,v)` shows all three colours (T011). Filed from
`problems/ca-3-3-k45-n91/attempts/01.md`.

```
coleq base out secs seed noise_pm weights workers [tenure] [kbase]
```
`secs` is WALL clock. `kbase` = base column count; if `base` has `kbase+1` columns the last one is used
as a warm-start colouring (so rounds can continue from the best object). Prints `cost=<empty
(class,colour) incidences>` — equal to the number of uncovered triples of the extended array.

Why it is fast and why it aims: with `b[r]` = #classes of r where r is the sole holder of its colour and
`g[r][z]` = #classes of r empty in z, `delta(r,z) = b[r] - g[r][z]` is exact and O(1); both are kept
incrementally (only 0<->1 count transitions update them), so a move costs O(members) instead of
O(members*C(K,2)). Breakout weights, min-conflicts + noise, and kicks from the best colouring on top.

**Control before any number counts:** drop the parent's last column and re-find it. It must pass in ~0 s.
**Both bugs this solver shipped with were caught by that rule and by counting the written file's defects
rather than reading its log** (a sign error made it maximise; a writeback bug printed cost=87 while
writing an all-zero column).
