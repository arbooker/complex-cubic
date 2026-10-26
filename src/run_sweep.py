#!/usr/bin/env python3
"""Parallel, resumable driver for ne_sweep, with the tiered escalation.

Tier 1 (cheap) runs over every field with a small modulus cap; almost all
fields are settled there.  Tier 2 re-runs whatever survived with a much larger
search, and tier 3 larger still.  Only fields that survive every tier are
reported as SURVIVOR, i.e. handed on to Phase 2 (proving M(K) < 1).

    python3 src/run_sweep.py                     # full sweep, all tiers
    python3 src/run_sweep.py --jobs 8
    python3 src/run_sweep.py --tier 1            # just the first tier
"""
import argparse, json, os, subprocess, sys, tempfile, time
from concurrent.futures import ThreadPoolExecutor

HERE    = os.path.dirname(os.path.abspath(__file__))
ROOT    = os.path.dirname(HERE)
SWEEP   = os.path.join(HERE, 'ne_sweep')
FIELDS  = os.path.join(ROOT, 'data', 'fields.csv')
OUTDIR  = os.path.join(ROOT, 'results')

# (name, extra ne_sweep arguments).  Each tier is applied only to the fields
# the previous tier left open.  |N(beta)| has to reach about 2*exp(2*rho_K),
# so the caps are chosen to cover successively larger rho_K: tier 1 settles
# rho_K <~ 4.3 (about 94% of the fields), tier 2 <~ 5.7, tier 3 <~ 6.9, and
# tier 4 the remaining handful (rho_K goes up to 8.2 in this range).
# All tiers parallelise over fields, one worker per chunk of the input.  At the
# tier-4 setting a worker peaks at about 620 MB (mostly the prime-generator
# cache), so a dozen of them fit comfortably; if that ever stops being true,
# ne_sweep's own `--jobs` forks workers over the modulus pool of a single field
# *after* the cache is built, so they share it copy-on-write -- pass
# {'per_field': k} in the tier options to use it.
TIERS = [
    ('tier1', ['--refine', '1.3', '--pbound',    '50000',
               '--nbcap',    '25000', '--budget',    '150000'], {}),
    ('tier2', ['--refine', '1.3', '--pbound',   '500000',
               '--nbcap',   '400000', '--budget',   '4000000'], {}),
    ('tier3', ['--refine', '1.0', '--pbound',  '3000000',
               '--nbcap',  '3000000', '--budget',  '30000000'], {}),
    ('tier4', ['--refine', '1.0', '--pbound', '20000000',
               '--nbcap', '20000000', '--budget', '150000000'],
              {'drop_known_euclidean': True}),
]

# Cavallar-Lemmermeyer's 58 conjecturally norm-Euclidean discriminants.  Phase 1
# cannot settle them by construction, so spending the expensive tiers on them is
# pure waste; they stay in the output as survivors, which is the correct verdict.
KNOWN_EUCLIDEAN = {
    -23, -31, -44, -59, -76, -83, -87, -104, -107, -108, -116, -135, -139,
    -140, -152, -172, -175, -200, -204, -211, -212, -216, -231, -239, -243,
    -244, -247, -255, -268, -300, -324, -356, -379, -411, -419, -424, -431,
    -440, -451, -460, -472, -484, -492, -499, -503, -515, -516, -519, -543,
    -628, -652, -687, -696, -728, -744, -771, -815, -876,
}


def run_chunk(rows, args, tag, shard=0):
    """Run one worker.  Its stdout goes straight to a per-shard .jsonl file
    rather than into a pipe: ne_sweep flushes after every field, so the file
    shows partial results while the sweep is still running, and a tier that is
    interrupted can be resumed from what is already on disk."""
    with tempfile.NamedTemporaryFile('w', suffix='.csv', delete=False) as fh:
        fh.write('label,disc,a0,a1,a2,a3\n')
        fh.writelines(rows)
        name = fh.name
    part = os.path.join(OUTDIR, f'{tag}.part{shard:02d}.jsonl')
    try:
        with open(name) as fin, open(part, 'a') as fout:
            p = subprocess.run([SWEEP] + args, stdin=fin, stdout=fout,
                               stderr=subprocess.PIPE, text=True)
        if p.returncode != 0:
            sys.stderr.write(f'[{tag}] ne_sweep exited {p.returncode}: '
                             f'{(p.stderr or "").strip()[-200:]}\n')
        return part
    finally:
        os.unlink(name)


def read_parts(tag):
    """Everything a tier has produced so far, across all its shards."""
    res = {}
    import glob
    for f in sorted(glob.glob(os.path.join(OUTDIR, f'{tag}.part*.jsonl'))):
        for line in open(f):
            line = line.strip()
            if not line:
                continue
            try:
                r = json.loads(line)
            except json.JSONDecodeError:
                continue                       # a truncated final line
            res[r['label']] = r
    return res


def sweep(rows, args, jobs, tag, per_field=False):
    """Split rows across `jobs` workers, one chunk each, interleaved so that
    the (rare, slow) survivors spread evenly over the workers.  With
    per_field=True the split is inverted: fields are done one at a time and the
    parallelism is handed to ne_sweep's own --jobs, which forks over the
    modulus pool and shares the prime cache."""
    if per_field:
        # `per_field` fields at a time, each getting jobs//per_field forked
        # modulus workers.  Filling the prime cache is single-threaded and
        # dominates, so running a few fields side by side is what keeps the
        # machine busy; the forked workers then share that field's cache.
        nf = max(1, min(per_field, jobs))
        args = args + ['--jobs', str(max(1, jobs // nf))]
        chunks = [[r] for r in rows]
    else:
        chunks = [rows[i::jobs] for i in range(jobs)]
        nf = jobs

    done = read_parts(tag)                     # resume: skip what is on disk
    todo = [[r for r in c if r.split(',')[0] not in done] for c in chunks]
    nskip = sum(len(c) for c in chunks) - sum(len(c) for c in todo)
    if nskip:
        sys.stderr.write(f'[{tag}] resuming: {nskip} fields already on disk\n')

    t0 = time.time()
    with ThreadPoolExecutor(max_workers=nf) as ex:
        list(ex.map(lambda ic: run_chunk(ic[1], args, tag, ic[0]),
                    [(i, c) for i, c in enumerate(todo) if c]))
    res = list(read_parts(tag).values())
    sys.stderr.write(f'[{tag}] {len(res)}/{len(rows)} fields in '
                     f'{time.time()-t0:.1f}s\n')
    return res


def main():
    global OUTDIR
    ap = argparse.ArgumentParser()
    ap.add_argument('--jobs', type=int, default=os.cpu_count())
    ap.add_argument('--fields', default=FIELDS)
    ap.add_argument('--refine', default=None,
                    help='override the per-tier --refine setting')
    ap.add_argument('--tier', type=int, default=0, help='run only this tier')
    ap.add_argument('--no-certify', action='store_true')
    ap.add_argument('--outdir', default=OUTDIR,
                    help='where to write results.jsonl / certs.txt / summary.csv')
    a = ap.parse_args()

    OUTDIR = a.outdir
    os.makedirs(OUTDIR, exist_ok=True)
    with open(a.fields) as fh:
        rows = [l for l in fh.readlines()[1:] if l.strip()]

    common = (['--no-certify'] if a.no_certify else [])
    if a.refine:
        TIERS[:] = [(t, [x for i, x in enumerate(e)
                         if e[max(i-1, 0)] != '--refine' and x != '--refine']
                        + ['--refine', a.refine], o) for t, e, o in TIERS]
    settled, pending = {}, rows
    tiers = TIERS if not a.tier else [TIERS[a.tier - 1]]

    for tag, extra, topt in tiers:
        if not pending:
            break
        todo = pending
        if topt.get('drop_known_euclidean'):
            todo = [r for r in pending
                    if int(r.split(',')[1]) not in KNOWN_EUCLIDEAN]
            sys.stderr.write(f'[{tag}] skipping {len(pending)-len(todo)} '
                             'conjecturally norm-Euclidean fields\n')
        res = sweep(todo, common + extra, a.jobs, tag,
                    per_field=topt.get('per_field', 0))
        by_label = {r['label']: r for r in res}
        nxt = []
        for row in pending:
            lab = row.split(',')[0]
            r = by_label.get(lab)
            if r is None:                        # crashed worker: retry later
                nxt.append(row)
            elif r['verdict'] == 'NOT_EUCLIDEAN':
                r['tier'] = tag
                settled[lab] = r
            else:
                r['tier'] = tag
                settled[lab] = r                 # keep the latest verdict
                nxt.append(row)
        pending = nxt
        sys.stderr.write(f'[{tag}] {len(pending)} fields still unsettled\n')
        write_output(settled)      # checkpoint after each tier


def write_output(settled):
    out = sorted(settled.values(), key=lambda r: -r['disc'])
    with open(os.path.join(OUTDIR, 'results.jsonl'), 'w') as fh:
        for r in out:
            fh.write(json.dumps(r) + '\n')

    ne   = [r for r in out if r['verdict'] == 'NOT_EUCLIDEAN']
    surv = [r for r in out if r['verdict'] == 'SURVIVOR']
    part = [r for r in out if r['verdict'] == 'PARTIAL']
    print(f'{len(out)} fields: {len(ne)} NOT_EUCLIDEAN, {len(surv)} SURVIVOR, '
          f'{len(part)} PARTIAL/ERROR')
    ds = sorted(r['disc'] for r in surv) if surv else []
    print('survivors:', ds if len(ds) <= 80 else
          f'{len(ds)} fields, |d| from {-max(ds)} to {-min(ds)} '
          f'(full list in results.jsonl)')

    # certificates, in the format understood by src/archive/verify.gp (archived)
    with open(os.path.join(OUTDIR, 'certs.txt'), 'w') as fh:
        for r in ne:
            if 'beta_pow' not in r:
                continue                          # settled by h_K > 1
            b, c = r['beta_pow'], r['alpha_pow']
            fh.write('[%d,%d,%d,%d,%d,%s,%s,%s,%s,%s,%s,%s,%s]\n' % (
                r['disc'], *r['poly'], b[0], b[1], b[2], b[3],
                c[0], c[1], c[2], c[3]))

    # summary table, one row per 1000 in |d_K|, in the shape of
    # Cavallar-Lemmermeyer, Euclidean windows, section 7
    with open(os.path.join(OUTDIR, 'summary.csv'), 'w') as fh:
        fh.write('d_lo,d_hi,not_euclidean,survivor,partial,total\n')
        step = 1000
        hi = max(-r['disc'] for r in out)
        for lo in range(0, hi + step, step):
            g = [r for r in out if lo < -r['disc'] <= lo + step]
            if not g:
                continue
            fh.write('%d,%d,%d,%d,%d,%d\n' % (
                lo, lo + step,
                sum(r['verdict'] == 'NOT_EUCLIDEAN' for r in g),
                sum(r['verdict'] == 'SURVIVOR' for r in g),
                sum(r['verdict'] == 'PARTIAL' for r in g), len(g)))


if __name__ == '__main__':
    main()
