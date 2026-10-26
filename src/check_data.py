#!/usr/bin/env python3
"""Check release integrity and coverage (not the mathematical proofs).

Use --write-manifest only when intentionally preparing a new release.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import sys

from test_known import EUCLIDEAN

ROOT = Path(__file__).resolve().parent.parent
CERT_HASH = "67e9f2cc4c9c83bff7bc66f97770fc11b5c2a0b878f70bb820a6c64a679a328d"
BINARIES = {"certify", "euclid", "euclid_mpfi", "euclid_platt",
            "ne_sweep", "ne_chain", "test_rounding"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def release_files(root):
    files = []
    for path in sorted(root.rglob("*")):
        rel = path.relative_to(root)
        if any(p in {"work", ".git", "__pycache__"} for p in rel.parts):
            continue
        if rel.as_posix() == "SHA256SUMS" or path.suffix == ".pyc":
            continue
        if rel.parent == Path("src") and path.name in BINARIES:
            continue
        require(not path.is_symlink(), "release contains symlink: " + str(rel))
        if path.is_file():
            files.append(rel.as_posix())
    return files


def manifest(root, write=False):
    names = release_files(root)
    actual = {name: hashlib.sha256((root / name).read_bytes()).hexdigest()
              for name in names}
    target = root / "SHA256SUMS"
    if write:
        target.write_text("".join(actual[name] + "  " + name + "\n" for name in names))
    else:
        expected = {}
        for line in target.read_text().splitlines():
            digest, name = line.split("  ", 1)
            require(name not in expected, "duplicate manifest path")
            expected[name] = digest
        require(actual == expected, "release checksums or file list differ from SHA256SUMS")
    return len(names)


def nf_records(path):
    records = {}
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("|")
        require(len(parts) == 6, "malformed algebraic-input record")
        key = (int(parts[0]), tuple(map(int, parts[1].split(","))))
        require(key not in records, "duplicate algebraic-input field")
        require(int(parts[5]) == key[0], "inconsistent field discriminant")
        require(int(parts[2]) > 0, "invalid basis denominator")
        require(len(list(map(int, parts[3].split(",")))) == 9, "invalid basis")
        require(len(list(map(int, parts[4].split(",")))) == 3, "invalid unit")
        records[key] = parts
    return records


def check(root):
    with (root / "data/fields.csv").open(newline="") as inp:
        rows = list(csv.DictReader(inp))
    keys = [(int(r["disc"]), tuple(int(r["a"+str(i)]) for i in range(4))) for r in rows]
    fields = set(keys)
    require(len(keys) == len(fields) == 17283, "incorrect field count or duplicate field")
    require(len({r["label"] for r in rows}) == len(rows), "duplicate field label")
    require(all(-170520 < d < 0 and p[3] == 1 for d, p in fields), "invalid field range")
    pari = nf_records(root / "data/nf/pari.dat")
    magma = nf_records(root / "data/nf/magma.dat")
    require(set(pari) == set(magma) == fields, "algebraic input does not match field list")
    payload = (root / "results/certificates.txt").read_bytes()
    require(hashlib.sha256(payload).hexdigest() == CERT_HASH, "wrong final certificate set")
    certificates = set()
    types = Counter()
    for line in payload.decode().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("|")
        require(len(parts) == 6, "malformed certificate record")
        key = (int(parts[0]), tuple(map(int, parts[1].split(","))))
        require(key in pari and parts[:5] == pari[key][:5], "certificate field data mismatch")
        require(key not in certificates, "duplicate certificate field")
        certificates.add(key)
        cert = parts[5].split()
        require((cert[0] == "b" and len(cert) == 7) or
                (cert[0] == "n" and len(cert) == 5), "invalid certificate format")
        list(map(int, cert[1:]))
        types[cert[0]] += 1
    require(types == {"b": 16825, "n": 400}, "incorrect certificate type counts")
    covering = [json.loads(s) for s in (root / "results/euclidean.jsonl").read_text().splitlines()]
    require(len(covering) == 58 and {r["disc"] for r in covering} == set(EUCLIDEAN),
            "incorrect covering list")
    require(all(r["status"] == "proved" for r in covering), "incomplete covering results")
    covered = {key for key in fields if key[0] in EUCLIDEAN}
    require(len(covered) == 58, "ambiguous covering discriminant")
    require(not certificates & covered and certificates | covered == fields,
            "certificates and covering fields do not partition the field list")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-manifest", action="store_true")
    args = parser.parse_args()
    check(ROOT)
    count = manifest(ROOT, args.write_manifest)
    print("PASS: %d release files; 17283 fields = 17225 certificates + 58 covering fields" % count)


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError, IndexError) as error:
        sys.exit(str(error))
