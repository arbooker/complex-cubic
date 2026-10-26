#!/usr/bin/env python3
"""Prove the 58 conjecturally norm-Euclidean fields ARE norm-Euclidean.

This is the converse half of the classification.  Cavallar and Lemmermeyer
established it in the 1990s with a floating-point covering carrying a 1% margin
"as a protection against rounding errors"; src/euclid redoes it in interval
arithmetic, discarding a cube only when the upper end of a rigorous enclosure
of sup|N| over the cube falls below 1, so insufficient precision costs a
subdivision rather than a false discard.

src/euclid is single threaded by design, so parallelism lives here: a pool of
one process per core, each field an independent job.  Fields differ in cost by
three orders of magnitude, which is why this is the right level for it -- the
work inside one field is not what needs spreading over cores.

The jobs are dispatched in decreasing order of difficulty.  With independent
jobs of very unequal length on a fixed number of cores, longest-first is the
standard scheduling rule -- it is what keeps the tail from being a single long
job started last while eleven cores stand idle.  For these fields the spread is
wide (seconds to minutes), so the ordering is worth a good deal of makespan.
It also puts the interesting fields first: -876, -744 and -696 are the ones
that can fail or expose a disagreement between the backends, and they are
finished before the small discriminants have begun.

Difficulty is not guessed but measured, by a first pass whose cube count is a
deterministic proxy for the work -- deterministic being the point, since a
wall-clock pre-pass on a busy machine would order the jobs by noise.

    python3 src/run_euclid.py --jobs 12
    python3 src/run_euclid.py --binary src/euclid      --out results/euclid_arb.jsonl
    python3 src/run_euclid.py --binary src/euclid_mpfi --out results/euclid_mpfi.jsonl
"""
import argparse, json, os, re, subprocess, sys, tempfile, time
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from test_known import EUCLIDEAN

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NFDATA = os.path.join(ROOT, "data", "nf", "pari.dat")
FAST = os.path.join(ROOT, "src", "euclid_platt")
ARB = os.path.join(ROOT, "src", "euclid")
ENV1 = dict(os.environ, OMP_NUM_THREADS="1")


def field_records(discs):
    want, out = set(str(d) for d in discs), {}
    with open(NFDATA) as f:
        for line in f:
            d = line.split("|", 1)[0]
            if d in want:
                if int(d) in out:
                    raise ValueError("ambiguous field discriminant: %s" % d)
                out[int(d)] = line
    missing = want - {str(d) for d in out}
    if missing:
        raise ValueError("missing field records: %s" % ", ".join(sorted(missing)))
    return out


@contextmanager
def atomic_results(path):
    """Publish only a complete run; interruption leaves the old file intact."""
    path = os.path.abspath(path)
    pending = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", delete=False,
                                         dir=os.path.dirname(path),
                                         prefix="." + os.path.basename(path) + ".") as out:
            pending = out.name
            yield out
            out.flush()
            os.fsync(out.fileno())
        os.replace(pending, path)
        pending = None
    finally:
        if pending is not None:
            os.unlink(pending)


def run(binary, path, depth, gbox):
    t0 = time.time()
    p = subprocess.run([binary, "--nf", path, "--k", "1.0",
                        "--depth", str(depth), "--gbox", str(gbox)],
                       capture_output=True, text=True, env=ENV1)
    return time.time() - t0, p.stdout


def parse(txt):
    def grab(pat):
        m = re.search(pat, txt)
        return int(m.group(1)) if m else -1
    return ("PROVED" in txt), grab(r"deepest=(\d+)"), grab(r"cubes=(\d+)")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 12)
    ap.add_argument("--out", default=os.path.join(ROOT, "work", "euclidean.jsonl"))
    ap.add_argument("--depth", type=int, default=20)
    ap.add_argument("--gbox", type=int, default=8)
    ap.add_argument("--tmp", default="/tmp")
    ap.add_argument("--binary", default=ARB,
                    help="prover to drive; src/euclid (Arb) by default. "
                         "src/euclid is the Arb build, src/euclid_mpfi the "
                         "MPFI one (see src/Makefile)")
    args = ap.parse_args()
    if args.jobs < 1:
        ap.error("--jobs must be positive")
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)

    # The published list is also used by test_known.py.  Results are output,
    # never the work list: a partial previous run must not remove any fields.
    discs = list(EUCLIDEAN)
    recs = field_records(discs)
    with tempfile.TemporaryDirectory(prefix="euclid-", dir=args.tmp) as tmp:
        paths = {}
        for d in discs:
            paths[d] = os.path.join(tmp, "ref%d.dat" % d)
            with open(paths[d], "w") as f:
                f.write(recs[d])
        return run_fields(args, discs, paths)


def run_fields(args, discs, paths):
    print("ordering %d fields by cube count..." % len(discs))
    cost = {}
    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        for d, (_, txt) in zip(discs, ex.map(
                lambda d: run(args.binary, paths[d], args.depth, args.gbox), discs)):
            cost[d] = parse(txt)[2]
    order = sorted(discs, key=lambda d: -cost[d])
    print("hardest first: %s ... easiest last: %s" %
          (", ".join("%d (%d cubes)" % (d, cost[d]) for d in order[:3]),
           ", ".join(str(d) for d in order[-3:])))
    print("running on %d cores\n" % args.jobs)

    t_start, done = time.time(), 0
    results = []
    with atomic_results(args.out) as out, ThreadPoolExecutor(max_workers=args.jobs) as ex:
        futures = [(d, ex.submit(run, args.binary, paths[d], args.depth, args.gbox))
                   for d in order]
        for d, fut in futures:
            secs, txt = fut.result()
            proved, depth, cubes = parse(txt)
            done += 1
            rec = {"disc": d, "status": "proved" if proved else "OUTSTANDING",
                   "depth": depth, "cubes": cubes, "secs": round(secs, 2)}
            results.append(rec)
            out.write(json.dumps(rec) + "\n")
            out.flush()
            print("[%2d/%2d] %7d  %-11s depth %2d  cubes %7d  %7.1fs  (elapsed %.0fs)"
                  % (done, len(order), d, rec["status"], depth, cubes, secs,
                     time.time() - t_start), flush=True)
            if not proved:
                print("  !! not proved -- investigate before trusting the rest",
                      flush=True)

    wall = time.time() - t_start
    summed_elapsed = sum(r["secs"] for r in results)
    nproved = sum(1 for r in results if r["status"] == "proved")
    print("\n%d/%d proved: %.1fs proof-pass wall time; "
          "%.1fs summed subprocess elapsed time (%d workers). "
          "Scheduling pre-pass is excluded from these timings."
          % (nproved, len(results), wall, summed_elapsed, args.jobs))
    return 0 if nproved == len(results) else 1


if __name__ == "__main__":
    sys.exit(main())
