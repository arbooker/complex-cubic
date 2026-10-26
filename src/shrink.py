#!/usr/bin/env python3
"""Replace each certificate's theta by the smallest representative of its coset.

Which representative of theta + beta O_K is stored is free -- the certifier
enumerates the coset either way -- so it may as well be the one that prints
shortest.  That is a closest-vector problem: reduce theta modulo the lattice
beta O_K, whose coordinate form is the multiplication matrix of beta.  LLL the
basis, then Babai nearest-plane.

This is a different question from the minimal-*norm* representative that
certify --argmin returns, and the two answers are far apart. The norm
is a product of embeddings, so |N| can be tiny while the embeddings are split
wide and the coordinates are enormous; minimising it is not minimising size and
on high-regulator fields it is actively terrible.  Here the objective is the
coordinate vector itself.

The coset is unchanged, so
m and |N beta| are unchanged, and theta not in beta O_K is a property of the
coset rather than the representative.  The output is re-certified from scratch.

    python3 src/shrink.py --jobs 12
"""
import argparse, os, sys, time
from fractions import Fraction


def polmul_mod(a, b, f):
    """(a*b) mod f, for f monic of degree 3 given as [a0,a1,a2,1]."""
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    while len(r) > 3:
        c = r.pop()
        d = len(r) - 3
        for k in range(3):
            r[d + k] -= c * f[k]
    while len(r) < 3:
        r.append(0)
    return r


def solve3(A, rhs):
    """Exact solve A y = rhs, A a 3x3 of Fractions."""
    M = [list(map(Fraction, A[i])) + [Fraction(rhs[i])] for i in range(3)]
    for c in range(3):
        p = next(r for r in range(c, 3) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(3):
            if r != c and M[r][c] != 0:
                fac = M[r][c]
                M[r] = [x - fac * y for x, y in zip(M[r], M[c])]
    return [M[i][3] for i in range(3)]


class Field(object):
    def __init__(self, f, den, basis):
        self.f, self.den, self.basis = f, den, basis
        # columns of this are omega_j as polynomials over den
        self.B = [[Fraction(basis[j][k], den) for j in range(3)] for k in range(3)]

    def to_poly(self, coords):
        p = [0, 0, 0]
        for j in range(3):
            for k in range(3):
                p[k] += coords[j] * self.basis[j][k]
        return p                      # element = p/den

    def from_poly(self, p, denom):
        """basis coordinates of the element p/denom"""
        y = solve3(self.B, [Fraction(c, denom) for c in p])
        assert all(v.denominator == 1 for v in y), 'not integral'
        return [int(v) for v in y]

    def mul(self, u, v):
        pr = polmul_mod(self.to_poly(u), self.to_poly(v), self.f)
        return self.from_poly(pr, self.den * self.den)

    def mulmat(self, u):
        """columns = coordinates of u*omega_j"""
        cols = []
        for j in range(3):
            e = [0, 0, 0]
            e[j] = 1
            cols.append(self.mul(u, e))
        return [[cols[j][i] for j in range(3)] for i in range(3)]


def lll(basis, delta=Fraction(99, 100)):
    """LLL on a list of integer column vectors (as rows here), exact."""
    b = [list(map(Fraction, v)) for v in basis]
    n = len(b)

    def gso(b):
        bs, mu = [], [[Fraction(0)] * n for _ in range(n)]
        for i in range(n):
            v = list(b[i])
            for j in range(i):
                d = sum(x * y for x, y in zip(bs[j], bs[j]))
                mu[i][j] = (sum(x * y for x, y in zip(b[i], bs[j])) / d) if d else Fraction(0)
                v = [x - mu[i][j] * y for x, y in zip(v, bs[j])]
            bs.append(v)
        return bs, mu

    bs, mu = gso(b)
    k = 1
    while k < n:
        for j in range(k - 1, -1, -1):
            q = mu[k][j]
            r = int(q + Fraction(1, 2)) if q >= 0 else -int(-q + Fraction(1, 2))
            if r:
                b[k] = [x - r * y for x, y in zip(b[k], b[j])]
                bs, mu = gso(b)
        nk = sum(x * x for x in bs[k])
        nk1 = sum(x * x for x in bs[k - 1])
        if nk >= (delta - mu[k][k - 1] ** 2) * nk1:
            k += 1
        else:
            b[k], b[k - 1] = b[k - 1], b[k]
            bs, mu = gso(b)
            k = max(k - 1, 1)
    return [[int(x) for x in v] for v in b]


def babai(basis, target):
    """nearest-plane: return target reduced modulo the lattice"""
    b = [list(map(Fraction, v)) for v in basis]
    n = len(b)
    bs, mu = [], []
    for i in range(n):
        v = list(b[i])
        row = []
        for j in range(i):
            d = sum(x * y for x, y in zip(bs[j], bs[j]))
            c = (sum(x * y for x, y in zip(b[i], bs[j])) / d) if d else Fraction(0)
            row.append(c)
            v = [x - c * y for x, y in zip(v, bs[j])]
        bs.append(v)
        mu.append(row)
    t = list(map(Fraction, target))
    for i in range(n - 1, -1, -1):
        d = sum(x * x for x in bs[i])
        if d == 0:
            continue
        c = sum(x * y for x, y in zip(t, bs[i])) / d
        r = int(c + Fraction(1, 2)) if c >= 0 else -int(-c + Fraction(1, 2))
        if r:
            t = [x - r * y for x, y in zip(t, b[i])]
    return [int(x) for x in t]


def shrink_line(line):
    fl = line.rstrip('\n').split('|')
    f = [int(x) for x in fl[1].split(',')]
    den = int(fl[2])
    bb = [int(x) for x in fl[3].split(',')]
    basis = [bb[0:3], bb[3:6], bb[6:9]]
    eps = [int(x) for x in fl[4].split(',')]
    K = Field(f, den, basis)
    tok = fl[5].split()
    if tok[0] == 'b':
        beta = [int(x) for x in tok[1:4]]
        th = [int(x) for x in tok[4:7]]
        lead = tok[:4]
    else:
        n = int(tok[1])
        th = [int(x) for x in tok[2:5]]
        acc = [1, 0, 0]
        for _ in range(n):
            acc = K.mul(acc, eps)
        beta = [acc[0] - 1, acc[1], acc[2]]     # eps^n - 1
        lead = tok[:2]
    L = K.mulmat(beta)                          # columns generate beta*O_K
    cols = [[L[i][j] for i in range(3)] for j in range(3)]
    red = lll(cols)
    small = babai(red, th)
    if sum(len(str(abs(x))) for x in small) < sum(len(str(abs(x))) for x in th):
        th = small
    fl[5] = ' '.join(lead + [str(x) for x in th])
    return '|'.join(fl) + '\n'


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=os.cpu_count())
    ap.add_argument('--src', default='results/certificates.txt')
    ap.add_argument('--out', required=True,
                    help='new output file; the input is never overwritten')
    a = ap.parse_args()
    if os.path.realpath(a.out) == os.path.realpath(a.src):
        ap.error('--out must differ from --src')
    if a.jobs < 1:
        ap.error('--jobs must be positive')
    head, recs = [], []
    for l in open(a.src):
        if l.startswith('#'):
            head.append(l)
        elif l.strip():
            recs.append(l)
    sys.stderr.write('%d certificates, %d workers\n' % (len(recs), a.jobs))
    t0 = time.time()
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=a.jobs) as ex:
        out = list(ex.map(shrink_line, recs, chunksize=64))
    with open(a.out, 'w') as fh:
        fh.writelines(head)
        fh.writelines(out)
    sys.stderr.write('done in %.1f min\n' % ((time.time() - t0) / 60))
    return 0


if __name__ == '__main__':
    sys.exit(main())
