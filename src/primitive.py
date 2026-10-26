#!/usr/bin/env python3
"""Make each Method I certificate primitive: divide out gcd(theta, beta).

A certificate bounds M(K, xi) for xi = theta/beta, and M depends only on xi.
So if delta = gcd(theta, beta) is not a unit the pair is not in lowest terms:
(beta/delta, theta/delta) presents the same xi, and since
m(theta' + beta' O) = m/|N delta| and |N beta'| = |N beta|/|N delta|, the
assertion m >= |N beta| survives unchanged.  The primitive pair is the
canonical form of a certificate.

Only Method I is worth doing.  Method II stores beta implicitly as
`n <n>` -- two tokens, because beta = eps^n - 1 is reconstructed -- so dividing
it out would force beta to be written explicitly as three large integers and
the certificate would grow.  Method I already writes beta out in full, so
dividing costs nothing in encoding.

It gains nothing in size either, and the reason is worth recording.  A
principal ideal of norm N does *not* generally have a generator of size
N^(1/3): it has one of size N^(1/3)*exp(Delta), where Delta is the offset of
log|sigma_1| from log(N)/3 along the unit orbit, and the orbit is a grid of
step R_K -- about 50 nats, or 22 digits, for these fields.  Generator size is
therefore governed by the regulator, not by the norm.  Dividing out delta
shrinks the norm by |N delta| ~ 3, worth (1/3)^(1/3) = 0.69, i.e. -0.16
digits, which is invisible against a +-12 digit spread from where beta and
beta/delta happen to land on that grid.  Measured over the 4406 affected
certificates: 1836 got smaller, 258 unchanged, 2314 larger, mean +0.55 digits.

Two details are needed to avoid making it actively worse:

  * bnfisprincipal returns a generator that is arbitrary within its unit
    orbit, so beta/delta must be unit-reduced.  Scaling theta and beta by the
    same unit leaves xi = theta/beta literally unchanged, so this is free.
  * the reduction must minimise T_2, not |Delta|.  T_2 is proportional to
    exp(2*Delta) + 2*exp(-Delta), which is asymmetric -- overshooting upward
    costs exp(Delta), downward only exp(-Delta/2) -- so on a coarse grid the
    optimum sits below zero.  Rounding Delta to nearest instead costs about
    0.3 digits per certificate.

Needs PARI (bnfinit), so unlike the other reductions it is GRH-conditional at
generation time.  That does not touch the result: src/certify re-derives
everything unconditionally and takes nothing here on trust.

    python3 src/primitive.py --src results/certificates.txt --jobs 12
"""
import argparse, glob, os, subprocess, sys, tempfile, time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
GP = os.path.join(HERE, 'primitive.gp')


def solve3(A, rhs):
    M = [list(map(Fraction, A[i])) + [Fraction(rhs[i])] for i in range(3)]
    for c in range(3):
        p = next(r for r in range(c, 3) if M[r][c] != 0)
        M[c], M[p] = M[p], M[c]
        pv = M[c][c]
        M[c] = [x / pv for x in M[c]]
        for r in range(3):
            if r != c and M[r][c] != 0:
                fa = M[r][c]
                M[r] = [x - fa * y for x, y in zip(M[r], M[c])]
    return [M[i][3] for i in range(3)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', default='results/certificates.txt')
    ap.add_argument('--out', required=True,
                    help='new output file; the input is never overwritten')
    ap.add_argument('--jobs', type=int, default=os.cpu_count())
    a = ap.parse_args()
    out = a.out
    if os.path.realpath(out) == os.path.realpath(a.src):
        ap.error('--out must differ from --src')
    if a.jobs < 1:
        ap.error('--jobs must be positive')

    head, recs = [], []
    for l in open(a.src):
        if l.startswith('#'):
            head.append(l)
        elif l.strip():
            recs.append(l.rstrip('\n').split('|'))

    tmp = tempfile.mkdtemp(prefix='prim')
    inp = os.path.join(tmp, 'in.txt')
    with open(inp, 'w') as fh:
        for i, fl in enumerate(recs):
            tok = fl[5].split()
            if tok[0] != 'b':
                continue
            den = int(fl[2])
            bb = [int(x) for x in fl[3].split(',')]
            basis = [bb[0:3], bb[3:6], bb[6:9]]

            def topoly(c):
                p = [0, 0, 0]
                for j in range(3):
                    for k in range(3):
                        p[k] += c[j] * basis[j][k]
                return p
            beta = [int(x) for x in tok[1:4]]
            th = [int(x) for x in tok[4:7]]
            fh.write('%d;%s;%d;%s;%s\n' % (i, fl[1], den,
                     ','.join(map(str, topoly(th))),
                     ','.join(map(str, topoly(beta)))))
    n = sum(1 for _ in open(inp))
    sys.stderr.write('%d Method I certificates, %d workers\n' % (n, a.jobs))

    subprocess.run(['split', '-n', 'l/%d' % a.jobs, '-d', inp,
                    os.path.join(tmp, 'c_')], check=True)
    t0 = time.time()
    procs = []
    for p in sorted(glob.glob(os.path.join(tmp, 'c_*'))):
        script = 'FIN="%s"; FOUT="%s";\n%s' % (p, p + '.out', open(GP).read())
        procs.append(subprocess.Popen(['gp', '-q', '-s', '400000000'],
                                      stdin=subprocess.PIPE,
                                      stdout=subprocess.DEVNULL,
                                      stderr=subprocess.DEVNULL,
                                      text=True).communicate(script))
    red = {}
    for fn in glob.glob(os.path.join(tmp, 'c_*.out')):
        for l in open(fn):
            i, tp, bp = l.strip().split(';')
            red[int(i)] = ([Fraction(x) for x in tp.split(',')],
                           [Fraction(x) for x in bp.split(',')])

    napp = 0
    for i, fl in enumerate(recs):
        if i not in red:
            continue
        den = int(fl[2])
        bb = [int(x) for x in fl[3].split(',')]
        basis = [bb[0:3], bb[3:6], bb[6:9]]
        B = [[Fraction(basis[j][k], den) for j in range(3)] for k in range(3)]
        tp, bp = red[i]
        tc, bc = solve3(B, tp), solve3(B, bp)
        if not all(v.denominator == 1 for v in tc + bc):
            continue
        fl[5] = ' '.join(['b'] + [str(int(v)) for v in bc]
                                + [str(int(v)) for v in tc])
        napp += 1
    with open(out, 'w') as fh:
        fh.writelines(head)
        for r in recs:
            fh.write('|'.join(r) + '\n')
    sys.stderr.write('made %d primitive, %.1f min\n' % (napp, (time.time() - t0) / 60))
    return 0


if __name__ == '__main__':
    sys.exit(main())
