#!/usr/bin/env python3
"""Compare the integral bases produced by PARI and by Magma.

The certifier checks everything it uses except one thing: that the basis it is
given spans the *maximal* order O_K.  That is the only assumed input, so it is
the one worth getting from two independent systems.

A basis of O_K is not unique, so the test is not equality of the records but
equality of the lattices they span: with A and B the two bases written over the
power basis of K, the change of basis A^{-1}B must be an integer matrix of
determinant +-1.  The field discriminants are compared too, since a
non-maximal order would show up there.

    python3 src/cmp_nfdata.py pari.txt magma.txt
"""
import sys
from fractions import Fraction


def parse(path):
    out = {}
    with open(path) as source:
        for line in source:
            line = line.strip()
            if not line or '|' not in line:
                continue
            f = line.split('|')
            if len(f) < 6:
                continue
            disc = int(f[0])
            poly = tuple(int(x) for x in f[1].split(','))
            den = int(f[2])
            ent = [int(x) for x in f[3].split(',')]
            eps = tuple(int(x) for x in f[4].split(','))
            dK = int(f[5])
            rows = [[Fraction(ent[3*i + j], den) for j in range(3)] for i in range(3)]
            out[(disc, poly)] = (rows, eps, dK)
    return out


def inv3(m):
    d = (m[0][0]*(m[1][1]*m[2][2] - m[1][2]*m[2][1])
         - m[0][1]*(m[1][0]*m[2][2] - m[1][2]*m[2][0])
         + m[0][2]*(m[1][0]*m[2][1] - m[1][1]*m[2][0]))
    if d == 0:
        return None, 0
    c = [[0]*3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            a = [[m[r][s] for s in range(3) if s != i] for r in range(3) if r != j]
            mi = a[0][0]*a[1][1] - a[0][1]*a[1][0]
            c[i][j] = mi * (-1)**(i + j) / d
    return c, d


def main(p1, p2):
    A, B = parse(p1), parse(p2)
    keys = sorted(set(A) & set(B))
    only = (set(A) ^ set(B))
    n_ok = n_bad = 0
    for k in keys:
        ra, ea, da = A[k]
        rb, eb, db = B[k]
        if da != db:
            print('  d=%d discriminants differ: %d vs %d' % (k[0], da, db))
            n_bad += 1
            continue
        ia, det = inv3(ra)
        if ia is None:
            print('  d=%d singular basis' % k[0]); n_bad += 1; continue
        # Rows are basis vectors, so B = T A and the change of basis is
        # T = B A^{-1}, not A^{-1} B.  With the factors the other way round the
        # test still passes whenever both bases happen to be triangular in the
        # same order, which is why getting it backwards showed up on only 3 of
        # the 17283 fields.
        T = [[sum(rb[i][l]*ia[l][j] for l in range(3)) for j in range(3)]
             for i in range(3)]
        integral = all(x.denominator == 1 for row in T for x in row)
        dT = (T[0][0]*(T[1][1]*T[2][2] - T[1][2]*T[2][1])
              - T[0][1]*(T[1][0]*T[2][2] - T[1][2]*T[2][0])
              + T[0][2]*(T[1][0]*T[2][1] - T[1][1]*T[2][0]))
        if integral and abs(dT) == 1:
            n_ok += 1
        else:
            n_bad += 1
            print('  d=%d bases span different lattices (integral=%s det=%s)'
                  % (k[0], integral, dT))
    print('  %d fields compared: %d agree, %d disagree%s'
          % (len(keys), n_ok, n_bad,
             ('; %d present in only one file' % len(only)) if only else ''))
    if not keys:
        print('  no common fields: comparison failed')
    return 1 if n_bad or only or not keys else 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1], sys.argv[2]))
