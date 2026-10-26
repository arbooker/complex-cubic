#!/usr/bin/env python3
"""Regression tests: the published values that any correct implementation of
Method I has to reproduce. See the accompanying paper and README.

    python3 src/test_known.py
"""
import json, os, subprocess, sys

HERE   = os.path.dirname(os.path.abspath(__file__))
ROOT   = os.path.dirname(HERE)
SWEEP  = os.path.join(HERE, 'ne_sweep')
FIELDS = os.path.join(ROOT, 'data', 'fields.csv')

# the 58 conjecturally norm-Euclidean discriminants (Cavallar-Lemmermeyer)
EUCLIDEAN = [
    -23, -31, -44, -59, -76, -83, -87, -104, -107, -108, -116, -135, -139,
    -140, -152, -172, -175, -200, -204, -211, -212, -216, -231, -239, -243,
    -244, -247, -255, -268, -300, -324, -356, -379, -411, -419, -424, -431,
    -440, -451, -460, -472, -484, -492, -499, -503, -515, -516, -519, -543,
    -628, -652, -687, -696, -728, -744, -771, -815, -876,
]

# d_K -> exact lower bound that must be met or beaten
BOUNDS = {
    -199:  (1, 1),            # M_1(K) = 1 exactly: not norm-Euclidean
    -351:  (1, 1),
    -367:  (1, 1),
    -327:  (101, 99),         # M_1(K) = 101/99
    -1055: (1483, 1370),      # published lower bound, must be matched exactly
    -999:  (1, 1),            # published 294557/272112; any bound >= 1 settles it
}


def run(discs, args):
    rows = {}
    with open(FIELDS) as fh:
        for line in fh.readlines()[1:]:
            p = line.split(',')
            if int(p[1]) in discs:
                rows[int(p[1])] = line
    missing = set(discs) - set(rows)
    assert not missing, f'not in fields.csv: {sorted(missing)}'
    inp = 'label,disc,a0,a1,a2,a3\n' + ''.join(rows[d] for d in sorted(rows))
    out = subprocess.run([SWEEP] + args, input=inp,
                         capture_output=True, text=True)
    if out.returncode != 0:
        sys.exit(f'ne_sweep failed: {out.stderr[-500:]}')
    return {r['disc']: r for r in map(json.loads, out.stdout.splitlines())}


def main():
    fails = []

    res = run(BOUNDS, ['--refine', '2.0', '--pbound', '500000',
                       '--nbcap', '400000', '--budget', '4000000'])
    for d, (num, den) in BOUNDS.items():
        r = res[d]
        if r['verdict'] != 'NOT_EUCLIDEAN':
            fails.append(f'd={d}: expected NOT_EUCLIDEAN, got {r["verdict"]}')
            continue
        gm, gn = (int(x) for x in r['M_lower_bound'].split('/'))
        if gm * den < num * gn:
            fails.append(f'd={d}: bound {r["M_lower_bound"]} weaker than '
                         f'{num}/{den}')
        else:
            print(f'  ok  d={d:6d}  M(K) >= {r["M_lower_bound"]:>12s} '
                  f'(published {num}/{den})')

    res = run(EUCLIDEAN, ['--refine', '1.0', '--pbound', '500000',
                          '--nbcap', '400000', '--budget', '4000000'])
    bad = [d for d in EUCLIDEAN if res[d]['verdict'] != 'SURVIVOR']
    if bad:
        fails.append('conjecturally norm-Euclidean fields wrongly settled: '
                     f'{bad}  (a certificate for one of these would be a '
                     'counterexample to Cavallar-Lemmermeyer, or a bug)')
    else:
        print(f'  ok  all {len(EUCLIDEAN)} conjecturally norm-Euclidean '
              'discriminants survive, as they must')
        worst = max(res[d].get('best_ratio', 0) for d in EUCLIDEAN)
        print(f'      (largest m/|N beta| reached among them: {worst:.4f})')

    if fails:
        print('\nFAILURES:')
        for f in fails:
            print('  ' + f)
        sys.exit(1)
    print('\nall regression tests passed')


if __name__ == '__main__':
    main()
