#!/usr/bin/env python3
"""Verify self-contained certificates; optionally prove exact coset minima.

The proof is performed by certify, not by this driver. Records are assigned
round-robin to single-threaded processes. Output is replaced atomically only
after every process succeeds and every record has the required verdict.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from run_euclid import atomic_results

ROOT = Path(__file__).resolve().parent.parent


def parse_output(output, records, argmin):
    lines = [s for s in output.splitlines() if s.strip() and not s.startswith("#")]
    width = 2 if argmin else 1
    if len(lines) != width * len(records):
        raise ValueError("wrong number of verifier output lines")
    parsed = []
    for i, record in enumerate(records):
        group = lines[width*i:width*(i+1)]
        if argmin:
            parts = group[0].split()
            if len(parts) != 4 or parts[0] != "MINIMAL":
                raise ValueError("missing MINIMAL representative")
            coordinates = [int(x) for x in parts[1:]]
        verdict = group[-1]
        if not verdict.startswith("CERTIFIED "):
            raise ValueError("certificate not accepted: " + verdict)
        if argmin and not verdict.startswith("CERTIFIED m = "):
            raise ValueError("exact minimum not certified: " + verdict)
        columns = record.split("|")
        row = dict(disc=int(columns[0]),
                   poly=[int(x) for x in columns[1].split(",")],
                   status="certified", msg=verdict)
        if argmin:
            row["minimum_coordinates"] = coordinates
        parsed.append(row)
    return parsed


def run_chunk(binary, records, argmin):
    command = [str(binary)] + (["--argmin"] if argmin else [])
    result = subprocess.run(command, input="\n".join(records) + "\n",
                            text=True, capture_output=True,
                            env=dict(os.environ, OMP_NUM_THREADS="1"))
    if result.returncode != 0:
        raise RuntimeError("verifier process failed (%d): %s\n%s" %
                           (result.returncode, result.stderr[-1000:],
                            result.stdout[-1000:]))
    return parse_output(result.stdout, records, argmin)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--src", type=Path, default=ROOT / "results/certificates.txt")
    parser.add_argument("--out", type=Path, default=ROOT / "work/certified.jsonl")
    parser.add_argument("--argmin", action="store_true")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    if args.src.resolve() == args.out.resolve():
        parser.error("input and output must be different files")
    records = [s for s in args.src.read_text().splitlines()
               if s.strip() and not s.startswith("#")]
    if not records:
        parser.error("no certificate records")
    jobs = min(args.jobs, len(records))
    chunks = [records[i::jobs] for i in range(jobs)]
    rows = [None] * len(records)
    start = time.perf_counter()
    print("%d certificates, %d workers%s" %
          (len(records), jobs, ", exact minima" if args.argmin else ""), flush=True)
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        pending = {pool.submit(run_chunk, ROOT / "src/certify", chunk, args.argmin): i
                   for i, chunk in enumerate(chunks)}
        for future in as_completed(pending):
            i = pending[future]
            for j, row in enumerate(future.result()):
                rows[i + jobs*j] = row
            print("chunk %d/%d verified" % (i+1, jobs), flush=True)
    if any(row is None for row in rows):
        raise RuntimeError("incomplete verifier output")
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with atomic_results(args.out) as out:
        for row in rows:
            out.write(json.dumps(row) + "\n")
    print("%d CERTIFIED%s; %.3f seconds wall time" %
          (len(rows), ", %d MINIMAL" % len(rows) if args.argmin else "",
           time.perf_counter()-start))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, RuntimeError) as error:
        sys.exit(str(error))
