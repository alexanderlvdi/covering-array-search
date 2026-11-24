# covering-array-search

Search programs behind Alexander LaVerdi's covering-array submissions to chrispbsmolen/ca-tools (November 2026).

- `ca_focus.c` (strength 3) and `ca_focus_t.c` (any strength 2-6): row-deletion simulated annealing with tabu repair. Given an array with N rows, delete a row and repair coverage by annealing; on success write an N-1 array.
- `ca_sweep.py` / `ca_sweep_t.py`: successive-halving sweep over many table entries at once (half the entries dropped each round, time per entry doubled).
- `ca_chain.c`: column extension at fixed N.
- `check_t.py`: independent exhaustive checker (numpy).

Build: `cc -O3 -march=native -o ca_focus ca_focus.c -lm` (same for the others). The searches are stochastic (seeded from the clock); the submitted arrays are the certificates. The searches were run by an AI-agent research system (Claude) that the author built.
License: MIT.
