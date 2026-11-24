#!/usr/bin/env python3
"""Build targets.tsv for ca_sweep_t.py from the ca-tools tracker copy, for the given strengths.
usage: make_targets_t.py <tracker-table-dir> <problem-dir> <t> [<t> ...]
Same filters as the strength-3 sweep: the seed must be the tracker's current best for (t,v,k), not a Torres-Jimenez
array or one of unknown origin (tracker Rule 6), above the v^t lower bound, plain (not gzipped), shape-checked, and
small enough for ca_focus_t (C(k,t) * v^t <= 2e7)."""
import collections, csv, math, os, shutil, sys

S, P, ts = sys.argv[1], sys.argv[2], {x for x in sys.argv[3:]}
os.makedirs(os.path.join(P, "seeds"), exist_ok=True)
best = collections.defaultdict(lambda: 10 ** 9)
for r in csv.DictReader(open(os.path.join(S, "ca_table.csv"))):
    key = (r["t"], r["v"], r["k"]); best[key] = min(best[key], int(r["N"]))
out, skip = [], collections.Counter()
for r in csv.DictReader(open(os.path.join(S, "MANIFEST.csv"))):
    if r["t"] not in ts:
        continue
    t, v, k, N, f, d = int(r["t"]), int(r["v"]), int(r["k"]), int(r["N"]), r["file"], r["source_detail"]
    if "TJ" in d or "Torres" in d or "CA_Database" in d: skip["tj-or-unknown-origin"] += 1; continue
    if best[(r["t"], r["v"], r["k"])] != N: skip["not the best"] += 1; continue
    if N <= v ** t: skip["at lower bound"] += 1; continue
    if math.comb(k, t) * v ** t > 2e7: skip["too big"] += 1; continue
    p = os.path.join(S, "arrays", f)
    if f.endswith(".gz") or not os.path.exists(p): skip["no plain file"] += 1; continue
    L = [l.split() for l in open(p) if l.strip()]
    if len(L) != N or any(len(x) != k for x in L): skip["shape mismatch"] += 1; continue
    out.append((t, v, k, N, f, d.replace("\t", " ")[:80]))
out.sort(key=lambda o: (math.comb(o[2], o[0]) * o[1] ** o[0], o[0], o[1]))
with open(os.path.join(P, "targets.tsv"), "w") as fh:
    fh.write("# seed-file\tt\tv\tk\tN\tsource (seed = current tracker entry; goal N-1)\n")
    for t, v, k, N, f, d in out:
        shutil.copy(os.path.join(S, "arrays", f), os.path.join(P, "seeds", f))
        fh.write(f"seeds/{f}\t{t}\t{v}\t{k}\t{N}\t{d}\n")
print(len(out), "targets;", dict(skip), "; by t:", dict(collections.Counter(o[0] for o in out)))
