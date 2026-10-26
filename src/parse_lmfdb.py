#!/usr/bin/env python3
"""Normalise an LMFDB number-field download into data/fields.csv.

Input lines look like

    "3.1.23.1"\t[1, 0, -1, 1]\t-23\t"3T2"\t[]

The polynomial is given by its coefficient list in ASCENDING order,
i.e. [a0,a1,a2,a3] means a0 + a1*x + a2*x^2 + a3*x^3.

Output CSV columns:  label,disc,a0,a1,a2,a3   (ascending, a3 == 1)
"""
import sys, csv, re, ast

LINE = re.compile(r'^"([^"]+)"\s+(\[[^\]]*\])\s+(-?\d+)\s+"([^"]*)"\s+(\[[^\]]*\])\s*$')

def main(src, dst):
    rows = []
    with open(src) as fh:
        for line in fh:
            if line.startswith('#') or not line.strip():
                continue
            m = LINE.match(line.rstrip('\n'))
            if not m:
                continue
            label, poly, disc, gal, clg = m.groups()
            coeffs = ast.literal_eval(poly)
            disc = int(disc)
            assert len(coeffs) == 4 and coeffs[3] == 1, (label, coeffs)
            assert disc < 0, label
            assert ast.literal_eval(clg) == [], (label, clg)   # h = 1
            rows.append((label, disc, *coeffs))
    rows.sort(key=lambda r: (-r[1], r[0]))
    with open(dst, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['label', 'disc', 'a0', 'a1', 'a2', 'a3'])
        w.writerows(rows)
    print(f'{len(rows)} fields -> {dst}', file=sys.stderr)

if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
