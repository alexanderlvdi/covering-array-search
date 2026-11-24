#!/usr/bin/env python3
"""Aimed sweep over many strength-3 covering-array table entries at once.
# search-kind: scan
# why: nobody knows which table entries are soft, so every entry with a stored seed array gets a short aimed
#      probe one row below the table (ca_focus), then the time is concentrated on the entries that came closest
#      (successive halving: half the entries dropped each round, time per entry doubled). Not a blind search:
#      each run is the aimed solver, and the allocation follows the measured uncovered counts.

usage: ca_sweep.py <targets.tsv> <outdir> <worker> <workers> <hours> [probe_seconds=60]
targets.tsv lines: seed-file <tab> v <tab> k <tab> N      (the seed passes with N rows; the goal is N-1, then N-2 ...)
Writes <outdir>/PASS_CA_<M>_3_<k>_<v>.txt for every level reached, <outdir>/sweep_w<i>.tsv (every run) and
<outdir>/sweep_w<i>.json (state after each round). Every PASS file still has to go through the frozen checker.
"""
import json, math, os, re, shutil, subprocess, sys, time

SA = os.environ.get("SA") or os.path.expanduser("~/research-workspace/library/solvers/covering-arrays/ca_focus")


def focus(M, k, v, init, out, seed, secs, tight):
    args = [SA, str(M), str(k), str(v), init, out, str(seed), str(int(secs))]
    args += ["500", "0.17", "5", "2", "6"] if tight else ["300", "0.20", "6"]
    try:
        p = subprocess.run(args, capture_output=True, text=True, timeout=secs + 180)
    except subprocess.TimeoutExpired:
        return None
    m = re.search(r"uncovered_tuples=(\d+)", p.stdout)
    return int(m.group(1)) if m else None


def main(a):
    tsv, out, w, W, hours = a[1], a[2], int(a[3]), int(a[4]), float(a[5])
    probe = float(a[6]) if len(a) > 6 else 60.0
    os.makedirs(out, exist_ok=True)
    end = time.time() + hours * 3600 - 120
    T = []
    for i, line in enumerate(l for l in open(tsv) if l.strip() and not l.startswith("#")):
        if i % W != w:
            continue
        f, v, k, N = line.rstrip("\n").split("\t")[:4]
        if not os.path.isabs(f):
            f = os.path.join(os.path.dirname(os.path.abspath(tsv)), f)      # seed paths are relative to targets.tsv
        tag = "v%s_k%s" % (v, k)
        T.append({"tag": tag, "v": int(v), "k": int(k), "N0": int(N), "M": int(N) - 1, "cur": f, "unc": None,
                  "passes": 0, "alive": True, "runs": 0})
    log = open(os.path.join(out, "sweep_w%d.tsv" % w), "a")
    seed, rnd = 1000 * (w + 1), 0

    def save():
        json.dump({"round": rnd, "targets": T}, open(os.path.join(out, "sweep_w%d.json" % w), "w"), indent=0)

    while time.time() + 150 < end and any(t["alive"] for t in T):
        secs = min(probe * 2 ** min(rnd, 10), 1500)
        live = sorted((t for t in T if t["alive"]), key=lambda t: (-t["passes"], t["unc"] if t["unc"] is not None else 10 ** 9))
        for t in live:
            if not t["alive"]:
                continue
            again = 0
            while again < 4:                       # after a pass, go straight for the next row down
                if time.time() + min(secs, 120) > end:
                    break
                budget = min(secs, max(30, end - time.time()))
                dst = os.path.join(out, "w%d_%s_M%d_%s.txt" % (w, t["tag"], t["M"], "b" if t["cur"].endswith("_a.txt") else "a"))
                seed += 7
                u = focus(t["M"], t["k"], t["v"], t["cur"], dst, seed, budget, t["unc"] is not None)
                t["runs"] += 1
                log.write("%s\t%s\tM=%d\tround=%d\tsecs=%d\tuncovered=%s\n" % (time.strftime("%H:%M:%S"), t["tag"], t["M"], rnd, budget, u))
                log.flush()
                if u == 0:
                    shutil.copy(dst, os.path.join(out, "PASS_CA_%d_3_%d_%d.txt" % (t["M"], t["k"], t["v"])))
                    t.update(cur=dst, M=t["M"] - 1, unc=None, passes=t["passes"] + 1)
                    again += 1
                    continue
                if u is None:                      # the solver failed (bad file, crash): count it, never spin on it
                    t["fails"] = t.get("fails", 0) + 1
                    if t["fails"] >= 2:
                        t["alive"] = False
                        log.write("%s\t%s\tdropped: solver failed twice\n" % (time.strftime("%H:%M:%S"), t["tag"]))
                if u is not None and (t["unc"] is None or u < t["unc"]):
                    t.update(unc=u, cur=dst)
                break
            save()
            if time.time() > end:
                break
        live = [t for t in T if t["alive"]]
        ranked = sorted(live, key=lambda t: (-t["passes"], t["unc"] if t["unc"] is not None else 10 ** 9))
        keep = max(3, math.ceil(len(ranked) / 2))
        for t in ranked[keep:]:
            t["alive"] = False
            log.write("%s\t%s\tM=%d\tdropped after round %d at uncovered=%s\n" % (time.strftime("%H:%M:%S"), t["tag"], t["M"], rnd, t["unc"]))
        log.flush()
        rnd += 1
        save()
    done = [(t["tag"], t["N0"], t["M"] + 1) for t in T if t["passes"]]
    print("SWEEP worker %d of %d: %d targets, rounds=%d, entries improved=%d %s" % (w, W, len(T), rnd, len(done), done))


if __name__ == "__main__":
    main(sys.argv)
