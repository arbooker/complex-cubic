#!/usr/bin/env python3
"""Phase 3: settle a field by sampling residues of eps^n - 1.

Why this works where the others did not
---------------------------------------
Everything reduces to finding one theta with

    m(theta + beta O_K) >= |N(beta)|,                                   (*)

since then M(K, theta/beta) >= 1 and K is not norm-Euclidean.  Phase 1 hunts
for a *small* modulus beta whose unit group mod beta is small, and enumerates
all |N(beta)| residues.  That fails for the fields left here for a concrete
reason: the good moduli are divisors of eps^t - 1, and finding them means
factoring a number of size e^{t Reg} -- e^{394} for d = -131183.

Take beta = eps^n - 1 itself instead.  Then eps^n = 1 (mod beta), so the unit
group mod beta has order exactly 2n, the smallest it can be, and no
factorisation is needed -- beta is handed to us.  The price is that |N(beta)| =
e^{n Reg} residues cannot be enumerated.  So sample them.

The success rate is not a guess.  The residues excluded by (*) are the w*gamma_a
for ideals a of norm below |N(beta)|, at most 2n per ideal, and there are
rho_K |N(beta)| such ideals, so the expected number of times a residue is
excluded is

    lambda = 2 n rho_K.

The same number comes out geometrically, which is the better check: the region
{|N(x)| < 1} has volume 2 pi R in one fundamental domain of the eps-action
(substitute s = log|N(x)|, t = log|x^(1)|; the Jacobian is 1/2 and the s
integral is 1), and O_K has covolume sqrt|d_K|/2, so a random xi has on average
4 pi R / sqrt|d_K| = 2 rho_K neighbours with |N(xi - eta)| < 1.  Hence

    P(a uniform theta satisfies (*)) ~ e^{-2 n rho_K},

which makes n = 1 optimal.  This is the right order but an over-estimate of the
rate, and knowingly so: lambda is a mean, and P(count = 0) is only e^{-lambda}
for a Poisson count, whereas lattice points in a region are far more regular
than Poisson, so the count concentrates near its mean and vanishing is rarer.
Measured: d = -156703 has 2 rho_K = 8.22, predicting 1 in 3700, and gives 1 in
4000; d = -10200 has 2 rho_K = 5.58, predicting 1 in 265, and gives 1 in ~2000.
Hence no budget is computed anywhere below -- see the note on scheduling.

The cost of a trial splits cleanly.  A doomed theta is abandoned at the first
coset point below |N(beta)| -- that is what --screen does -- so it inspects
nsl/lambda = (R/log 2)/(2 rho_K) = sqrt|d_K| / (4 pi log 2) slices, about 45
here, independent of the regulator.  What does grow with R is the cost of a
slice, since the working precision must carry the range e^{n R} of |x^(1)|; a
trial is O(R_K sqrt|d_K|), measured at 0.025 ms for d = -10200 (R = 45) and
0.97 ms for d = -131183 (R = 394).

Screening decides |N| by exact integer arithmetic, so it never discards a good
theta; the floating-point geometry only chooses where to look.  Each survivor
is nonetheless re-evaluated, and the enumeration re-run with a wider box to
confirm it was exhaustive, before anything is recorded.

    python3 src/sample.py --all --jobs 12
    python3 src/sample.py --disc -156703
"""
import argparse, json, math, os, subprocess, sys, threading, time

HERE  = os.path.dirname(os.path.abspath(__file__))
ROOT  = os.path.dirname(HERE)
CHAIN = os.path.join(HERE, 'ne_chain')
MAIN  = os.path.join(ROOT, 'results', 'results.jsonl')
OUT   = os.path.join(ROOT, 'results', 'sample')

KNOWN_EUCLIDEAN = {
    -23, -31, -44, -59, -76, -83, -87, -104, -107, -108, -116, -135, -139,
    -140, -152, -172, -175, -200, -204, -211, -212, -216, -231, -239, -243,
    -244, -247, -255, -268, -300, -324, -356, -379, -411, -419, -424, -431,
    -440, -451, -460, -472, -484, -492, -499, -503, -515, -516, -519, -543,
    -628, -652, -687, -696, -728, -744, -771, -815, -876,
}


def orbit_lengths(reg, floor=12.0, count=6):
    """The orbit lengths to try, in order of cost.

    A single n will not do, and d_K = -680 is the proof: its maximiser is fixed
    by eps^3, so it is *invisible* at every n not divisible by 3 -- not merely
    expensive, unreachable.  Picking one n from a size heuristic (n = 5 there)
    means the field never settles however long it runs.  d_K = -327 is the same
    story at eps^2.

    So walk n upward from the least one whose residue ring can hold an excluded
    residue at all (e^{n Reg} >~ e^{12}).  lambda = 2 n rho_K makes each
    increment cost a factor e^{2 rho_K}, so small n first; but every n has to be
    tried, because divisibility, not size, decides whether a given maximiser is
    reachable at all.
    """
    n0 = max(1, int(math.ceil(floor / max(reg, 1.0))))
    return [n0 + k for k in range(count)]


# There is deliberately no budget model here.  e^{-2 n rho_K} is the right
# order for the success rate -- it is E[count] under the volume computation --
# but it is an over-estimate of P(count = 0), because lattice points in a
# region are far more regular than a Poisson process, so the count concentrates
# near its mean and vanishing is rarer than e^{-mean}.  Measured: d = -10200
# has 2 rho_K = 5.58, predicting 1 success in 265, and delivers 1 in ~2000.
# Rather than calibrate that gap, allocate trials round-robin: every unsettled
# field gets the next chunk before any field gets the one after.  A field then
# receives exactly as much work as it needs, the total is the sum over fields of
# what each needs, and nothing depends on predicting the rate.


def _eval(ps, n, theta, widen=0, timeout=900):
    cmd = [CHAIN, '--poly', ps, '--n', str(n)]
    if widen:
        cmd += ['--widen', str(widen)]
    q = subprocess.run(cmd, input=' '.join(theta) + '\n',
                       capture_output=True, text=True, timeout=timeout)
    w = q.stdout.split()
    if len(w) < 2 or not int(w[1]):
        return None
    return int(w[0]), int(w[1])


def confirm(ps, n, theta, timeout=900):
    """Re-evaluate a screened survivor, and check the enumeration was complete.

    A slice contributes coset points with |y0| in [1/2,1) and y1^2+y2^2 <= 1.05,
    so ||y|| <= sqrt(2.05) and each coefficient of u = tgt - W^{-1}y lies within
    sqrt(2.05)*||row of W^{-1}|| + 1/2 of round(tgt).  The default box is that
    bound plus margin, so it is exhaustive and its minimum is the true
    m(theta + beta O_K); widening it must therefore change nothing, and that is
    what is checked here.

    Fincke-Pohst (--exact) is not used for this.  It agrees wherever it runs,
    but the target has e^{n Reg}-sized entries -- 48 digits for d = -98488 --
    and qfcvp0's Cholesky loses them, returning an empty list rather than an
    error: it silently reported "no point below the cap" for two certificates
    that are perfectly good.
    """
    a = _eval(ps, n, theta, 0, timeout)
    if a is None:
        return None
    b = _eval(ps, n, theta, 3, timeout)
    if b != a:
        sys.stderr.write('  UNSTABLE enumeration for %s n=%d: %r vs %r\n'
                         % (ps, n, a, b))
        return None
    return a


def main():
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--disc', type=int)
    g.add_argument('--all', action='store_true')
    ap.add_argument('--jobs', type=int, default=os.cpu_count())
    ap.add_argument('--chunk', type=int, default=200000,
                    help='trials per task; the unit of both parallelism and '
                         'round-robin fairness, so equal trials -- not equal '
                         'time -- is the right measure')
    ap.add_argument('--rounds', type=int, default=100000)
    ap.add_argument('--stubborn', type=int, default=40,
                    help='rounds at the cheapest orbit length before trying '
                         'others; only a field whose maximiser is unreachable '
                         'at that n needs the others, and those are rare')
    ap.add_argument('--limit', type=int, default=0)
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--results', default=None,
                    help='results.jsonl to read (default results/results.jsonl)')
    ap.add_argument('--outdir', default=None,
                    help='where to write part.jsonl (default results/sample)')
    ap.add_argument('--keep-euclidean', action='store_true',
                    help='do not skip the 58 conjectured discriminants')
    a = ap.parse_args()

    main_path = a.results or MAIN
    outdir = a.outdir or OUT
    rs = [json.loads(l) for l in open(main_path)]
    if a.disc is not None:
        todo = [r for r in rs if r['disc'] == a.disc]
    else:
        todo = [r for r in rs if r['verdict'] != 'NOT_EUCLIDEAN'
                and (a.keep_euclidean or r['disc'] not in KNOWN_EUCLIDEAN)]
    todo.sort(key=lambda r: r['rho'])          # cheapest first
    if a.limit:
        todo = todo[:a.limit]

    os.makedirs(outdir, exist_ok=True)
    part = os.path.join(outdir, 'part.jsonl')
    done = set()
    if os.path.exists(part):
        for line in open(part):
            try:
                d = json.loads(line)
            except json.JSONDecodeError:
                continue
            if d.get('result') == 'NOT_EUCLIDEAN':
                done.add((d['disc'], tuple(d['poly'])))
    todo = [r for r in todo if (r['disc'], tuple(r['poly'])) not in done]

    sys.stderr.write('%d fields, %d workers, %d trials/chunk, %d rounds max\n'
                     % (len(todo), a.jobs, a.chunk, a.rounds))
    if a.dry_run:
        for r in todo[:20]:
            sys.stderr.write('  d=%-9d rho %5.2f reg %6.1f  n in %s\n'
                             % (r['disc'], r['rho'], r['reg'],
                                orbit_lengths(r['reg'])))
        return

    lock = threading.Lock()
    settled = set()                      # fields already proven; skip their tasks
    nset = [0]
    fh = open(part, 'a')
    t0 = time.time()

    def gen():
        # The orbit length is *not* cycled from the start.  Each increment of n
        # costs a factor e^{2 rho_K} in the success rate -- at rho_K = 7 that is
        # 1e-6 per step -- so for the fields that are actually slow, every chunk
        # spent at n > n0 is thrown away, and rotating over six values of n made
        # them six times slower.  Only a field that has resisted n0 for a long
        # time can be one whose maximiser is unreachable there (d_K = -680 is
        # fixed by eps^3 and invisible at any n not divisible by 3), and only
        # then is it worth spreading the effort.
        for c in range(a.rounds):
            for r in todo:               # already sorted by rho: cheapest first
                if (r['disc'], tuple(r['poly'])) in settled:
                    continue
                ns = orbit_lengths(r['reg'])
                n = ns[0] if c < a.stubborn else ns[(c - a.stubborn) % len(ns)]
                yield r, n, a.chunk, c
    tasks = gen()

    def work():
        while True:
            with lock:
                try:
                    r, n, chunk, c = next(tasks)
                except StopIteration:
                    return
                key = (r['disc'], tuple(r['poly']))
                if key in settled:
                    continue
            ps = ','.join(map(str, r['poly']))
            seed = abs(hash((r['disc'], tuple(r['poly']), c))) % (2**31 - 1) + 1
            try:
                p = subprocess.run([CHAIN, '--poly', ps, '--n', str(n),
                                    '--screen', '--random', str(chunk),
                                    '--seed', str(seed)],
                                   capture_output=True, text=True)
            except Exception as e:                       # noqa: BLE001
                sys.stderr.write('  d=%d chunk %d failed: %r\n'
                                 % (r['disc'], c, e))
                continue
            for line in p.stdout.splitlines():
                if not line.startswith('SURVIVOR'):
                    continue
                w = line.split()
                theta = w[2:5]
                try:
                    got = confirm(ps, n, theta)
                except subprocess.TimeoutExpired:
                    got = None
                if not got:
                    continue
                m, N = got
                if m < N:
                    continue                  # screening artefact, keep going
                gg = math.gcd(m, N)
                rec = {'disc': r['disc'], 'poly': r['poly'], 'label': r['label'],
                       'rho': r['rho'], 'reg': r['reg'], 'n': n, 'seed': seed,
                       'theta': theta, 'm': str(m), 'Nbeta': str(N),
                       'M_lower_bound': '%d/%d' % (m // gg, N // gg),
                       'ratio': m / N, 'method': 'sample',
                       'result': 'NOT_EUCLIDEAN'}
                with lock:
                    key = (r['disc'], tuple(r['poly']))
                    if key in settled:
                        break
                    settled.add(key)
                    nset[0] += 1
                    fh.write(json.dumps(rec) + '\n')
                    fh.flush()
                    sys.stderr.write(
                        '  d=%-9d rho %5.2f  M(K) >= %.6f   settled %d/%d'
                        '  (%.0f min)\n'
                        % (r['disc'], r['rho'], m / N, nset[0], len(todo),
                           (time.time() - t0) / 60))
                break

    ths = [threading.Thread(target=work, daemon=True) for _ in range(a.jobs)]
    for t in ths:
        t.start()
    for t in ths:
        t.join()
    sys.stderr.write('settled %d of %d in %.1f min\n'
                     % (nset[0], len(todo), (time.time() - t0) / 60))


if __name__ == '__main__':
    main()
