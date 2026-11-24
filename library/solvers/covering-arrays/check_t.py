#!/usr/bin/env python3
"""Independent covering-array check for any strength: check_t.py <array> t v  -> prints missing tuple count."""
import itertools, sys
import numpy as np
A = np.loadtxt(sys.argv[1], dtype=np.int64, ndmin=2); t = int(sys.argv[2]); v = int(sys.argv[3])
N, k = A.shape
assert A.min() >= 0 and A.max() < v
miss = 0
for cs in itertools.combinations(range(k), t):
    code = np.zeros(N, dtype=np.int64)
    for c in cs:
        code = code * v + A[:, c]
    miss += v ** t - len(np.unique(code))
print(f"N={N} t={t} k={k} v={v} missing={miss} -> {'PASS' if miss == 0 else 'FAIL'}")
