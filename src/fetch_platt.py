#!/usr/bin/env python3
"""Download and checksum the unmodified, optional Platt header dependency."""
import argparse
import hashlib
from pathlib import Path
import sys
import urllib.request

REVISION = "42b21426718e542daa2b006dc05ea2d7f26426e6"
URL = ("https://raw.githubusercontent.com/djplatt/code/" + REVISION +
       "/includes/int_double14.2.h")
SHA256 = "37cc444d50998799e4783a3fb8b7ad038d379c6ddebf5827d41f320278ea308e"
ROOT = Path(__file__).resolve().parent.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--destination", type=Path, default=ROOT / "work/deps/platt")
    args = parser.parse_args()
    target = args.destination / "int_double14.2.h"
    if target.exists():
        payload = target.read_bytes()
    else:
        print("Downloading " + URL, flush=True)
        with urllib.request.urlopen(URL, timeout=60) as response:
            payload = response.read()
    if hashlib.sha256(payload).hexdigest() != SHA256:
        sys.exit("Platt header checksum mismatch; no file was replaced")
    if not target.exists():
        target.parent.mkdir(parents=True, exist_ok=True)
        # Exclusive creation protects an existing local copy from replacement.
        with target.open("xb") as out:
            out.write(payload)
    print("Verified unmodified Platt header: " + str(target))


if __name__ == "__main__":
    main()
