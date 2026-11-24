#!/usr/bin/env python3
"""Aimed sweep over many covering-array table entries of any strength (2..6), one row below each.
# search-kind: scan
# why: same as ca_sweep.py (strength 3, 6 records on 2026-11-23), extended to strengths 2 and 4+ with ca_focus_t:
#      a short aimed probe one row below every entry with a stored seed, then successive halving toward the closest.

usage: ca_sweep_t.py <targets.tsv> <outdir> <worker> <workers> <hours> [probe_seconds=60]
targets.tsv lines: seed-file <tab> t <tab> v <tab> k <tab> N   (seed passes with N rows; goal N-1, then lower)
Writes <outdir>/PASS_CA_<M>_<t>_<k>_<v>.txt per level reached, sweep_w<i>.tsv (every run), sweep_w<i>.json (state).
Every PASS file still has to go through the frozen checker.
"""
import json, math, os, re, shutil, subprocess, sys, time

SA = os.environ.get("SA_T") or os.path.expanduser("~/research-workspace/library/solvers/covering-arrays/ca_focus_t")


def focus(t, M, k, v, init, out, seed, secs, tight):
    args = [SA, str(t), str(M), str(k), str(v), init, out, str(seed), str(int(secs))]
    args += ["500", "0.17", "5"] if tight else ["300", "0.20", "6"]
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
        f, t, v, k, N = line.rstrip("\n").split("\t")[:5]
        if not os.path.isabs(f):
            f = os.path.join(os.path.dirname(os.path.abspath(tsv)), f)
        T.append({"tag": f"t{t}_v{v}_k{k}", "t": int(t), "v": int(v), "k": int(k), "N0": int(N), "M": int(N) - 1,
                  "cur": f, "unc": None, "passes": 0, "alive": True, "runs": 0, "fails": 0})
    log = open(os.path.join(out, "sweep_w%d.tsv" % w), "a")
    seed, rnd = 1000 * (w + 1), 0

    def save():
        json.dump({"round": rnd, "targets": T}, open(os.path.join(out, "sweep_w%d.json" % w), "w"), indent=0)

    while time.time() + 150 < end and any(x["alive"] for x in T):
        secs = min(probe * 2 ** min(rnd, 10), 1500)
        live = sorted((x for x in T if x["alive"]), key=lambda x: (-x["passes"], x["unc"] if x["unc"] is not None else 10 ** 9))
        for x in live:
            if not x["alive"]:
                continue
            again = 0
            while again < 4:
                if time.time() + min(secs, 120) > end:
                    break
                budget = min(secs, max(30, end - time.time()))
                dst = os.path.join(out, "w%d_%s_M%d_%s.txt" % (w, x["tag"], x["M"], "b" if x["cur"].endswith("_a.txt") else "a"))
                seed += 7
                u = focus(x["t"], x["M"], x["k"], x["v"], x["cur"], dst, seed, budget, x["unc"] is not None)
                x["runs"] += 1
                log.write("%s\t%s\tM=%d\tround=%d\tsecs=%d\tuncovered=%s\n" % (time.strftime("%H:%M:%S"), x["tag"], x["M"], rnd, budget, u))
                log.flush()
                if u == 0:
                    shutil.copy(dst, os.path.join(out, "PASS_CA_%d_%d_%d_%d.txt" % (x["M"], x["t"], x["k"], x["v"])))
                    x.update(cur=dst, M=x["M"] - 1, unc=None, passes=x["passes"] + 1)
                    again += 1
                    continue
                if u is None:
                    x["fails"] += 1
                    if x["fails"] >= 2:
                        x["alive"] = False
                        log.write("%s\t%s\tdropped: solver failed twice\n" % (time.strftime("%H:%M:%S"), x["tag"]))
                elif x["unc"] is None or u < x["unc"]:
                    x.update(unc=u, cur=dst)
                break
            save()
            if time.time() > end:
                break
        live = [x for x in T if x["alive"]]
        ranked = sorted(live, key=lambda x: (-x["passes"], x["unc"] if x["unc"] is not None else 10 ** 9))
        for x in ranked[max(3, math.ceil(len(ranked) / 2)):]:
            x["alive"] = False
            log.write("%s\t%s\tM=%d\tdropped after round %d at uncovered=%s\n" % (time.strftime("%H:%M:%S"), x["tag"], x["M"], rnd, x["unc"]))
        log.flush()
        rnd += 1
        save()
    done = [(x["tag"], x["N0"], x["M"] + 1) for x in T if x["passes"]]
    print("SWEEP worker %d of %d: %d targets, rounds=%d, entries improved=%d %s" % (w, W, len(T), rnd, len(done), done))


if __name__ == "__main__":
    main(sys.argv)
