/* ne_sweep.c -- Phase 1 of the norm-Euclidean complex cubic sweep.
 *
 * For a complex cubic field K = Q[x]/(f) this program searches for a
 * *non-Euclideanity certificate*: a pair (beta, alpha) of algebraic integers,
 * beta a non-zero non-unit, with
 *
 *      |N(alpha - beta*eta)| >= |N(beta)|      for every eta in O_K.        (*)
 *
 * Since M(K, alpha/beta) = min_eta |N(alpha-beta*eta)| / |N beta|, (*) gives
 * M(K) >= 1 and hence K is not norm-Euclidean.
 *
 * How (*) is decided exactly.  Put b = (beta) and R = O_K/b, a finite ring of
 * order |N beta|.  A gamma with 0 < |N gamma| < N(b) generates an ideal a of
 * norm < N(b); conversely the generators of a are exactly {w*gamma_a : w in
 * O_K^*}, and O_K^* = {+- eps^j} for a fundamental unit eps.  Hence the set of
 * residue classes represented by such gamma is the finite, exactly computable
 * set
 *
 *      S = { w * (gamma_a mod b) : a principal, 0 < N(a) < N(b), w in U },
 *      U = <-1, eps> mod b.
 *
 * Any class outside S other than the class of 0 witnesses (*) -- alpha need
 * not be invertible mod b, a zero divisor is just as good, and admitting those
 * classes is what makes small primes useful.  We record, for every class, the
 * least |N gamma| realising it, which also gives the quantitative bound
 * M(K) >= m / |N beta|.
 *
 * Choice of b.  We take b = prod p_i^{k_i} over distinct prime ideals with the
 * order of eps mod b small: with u = |U| there are ~ N(b)/u orbits of classes
 * and ~ rho*N(b) ideals of norm < N(b), where
 * rho = Res_{s=1} zeta_K = 2*pi*h*Reg/sqrt|d_K|.  So the expected number of
 * unmarked orbits is (N(b)/u) * exp(-rho*u): small u is what matters, and the
 * natural source of small u is the factorisation of eps^t - 1, whose small
 * prime power divisors we find by trial-dividing the integer |N(eps^t - 1)|.
 *
 * Build: see src/Makefile. Docs: README and the accompanying paper.
 */

#include <pari/pari.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <signal.h>
#include <sys/wait.h>
#include <sys/mman.h>

#define DEG      3      /* cubic fields only */
#define MAXCAND  1024   /* candidate prime powers for the modulus search */
#define MAXSEL   12     /* prime powers in the chosen modulus */
#define UMAX     64     /* reject a modulus where |<-1,eps>| exceeds this */

static int g_debug = 0;
#define DBG(...) do { if (g_debug) fprintf(stderr, __VA_ARGS__); } while (0)

/* ================================================================== */
/* the residue ring R = O_K/b, with b given by its HNF basis           */
/* ================================================================== */

typedef struct {
  long H[DEG][DEG];              /* HNF basis of b, columns; upper triangular */
  long diag[DEG];                /* H[i][i] */
  long long stride[DEG];         /* 1, diag[0], diag[0]*diag[1] */
  long long NB;                  /* prod diag = N(b) */
  long long mt[DEG][DEG][DEG];   /* w_i*w_j = sum_k mt[i][j][k]*w_k */
} ring_t;

typedef long res_t[DEG];         /* reduced residue: 0 <= v[i] < diag[i] */

static void ring_reduce(const ring_t *R, __int128 *v, res_t out)
{
  int i, k;
  for (i = DEG - 1; i >= 0; i--) {
    __int128 d = R->diag[i], q = v[i] / d;
    if (v[i] % d < 0) q--;                        /* floor division */
    if (q) for (k = 0; k <= i; k++) v[k] -= q * R->H[k][i];
  }
  for (i = 0; i < DEG; i++) out[i] = (long)v[i];
}

static void ring_mul(const ring_t *R, const res_t a, const res_t b, res_t out)
{
  __int128 v[DEG] = {0, 0, 0};
  int i, j, k;
  for (i = 0; i < DEG; i++) {
    if (!a[i]) continue;
    for (j = 0; j < DEG; j++) {
      long long c;
      if (!b[j]) continue;
      c = (long long)a[i] * b[j];
      for (k = 0; k < DEG; k++)
        if (R->mt[i][j][k]) v[k] += (__int128)c * R->mt[i][j][k];
    }
  }
  ring_reduce(R, v, out);
}

static long long ring_index(const ring_t *R, const res_t v)
{ return v[0] * R->stride[0] + v[1] * R->stride[1] + v[2] * R->stride[2]; }

static void ring_unindex(const ring_t *R, long long idx, res_t v)
{
  v[2] = (long)(idx / R->stride[2]); idx %= R->stride[2];
  v[1] = (long)(idx / R->stride[1]); idx %= R->stride[1];
  v[0] = (long)idx;
}

static void ring_scalar(const ring_t *R, long s, res_t v)
{ __int128 w[DEG] = {0, 0, 0}; w[0] = s; ring_reduce(R, w, v); }

/* ---- PARI glue ---- */

static int ring_set_multable(ring_t *R, GEN nf)
{
  int i, j, k;
  for (j = 0; j < DEG; j++) {
    GEN M = zk_multable(nf, col_ei(DEG, j + 1));   /* column i = w_i*w_j */
    for (i = 0; i < DEG; i++)
      for (k = 0; k < DEG; k++) {
        GEN c = gcoeff(M, k + 1, i + 1);
        if (is_bigint(c)) return 0;
        R->mt[i][j][k] = itos(c);
      }
  }
  return 1;
}

/* fill the lattice part of R from an ideal in HNF; 0 if it does not fit */
static int ring_set_lattice(ring_t *R, GEN hnf)
{
  int i, j;
  GEN N = gen_1;
  for (j = 0; j < DEG; j++)
    for (i = 0; i < DEG; i++) {
      GEN c = gcoeff(hnf, i + 1, j + 1);
      if (is_bigint(c)) return 0;
      R->H[i][j] = itos(c);
    }
  for (i = 0; i < DEG; i++) {
    R->diag[i] = R->H[i][i];
    if (R->diag[i] <= 0) return 0;
    N = mulis(N, R->diag[i]);
  }
  if (is_bigint(N)) return 0;
  R->NB = itos(N);
  R->stride[0] = 1;
  R->stride[1] = R->diag[0];
  R->stride[2] = (long long)R->diag[0] * R->diag[1];
  return 1;
}

/* reduce an algebraic integer (t_COL over Z in the integral basis, possibly
 * with huge entries -- eps has coefficients of size e^{Reg/2}) into R */
static void ring_from_nfelt(const ring_t *R, GEN hnf, GEN x, res_t out)
{
  pari_sp av = avma;
  GEN r = ZC_hnfrem(x, hnf);
  __int128 v[DEG];
  int i;
  for (i = 0; i < DEG; i++) v[i] = (__int128)itos(gel(r, i + 1));
  ring_reduce(R, v, out);
  set_avma(av);
}

/* U = <-1, eps> inside R^*, listed; returns |U|, or 0 if |U| > cap */
static int unit_orbit(const ring_t *R, GEN hnf, GEN eps, res_t *U, int cap)
{
  res_t e, x, one, mone, t;
  int n = 0, i, k, w;
  ring_scalar(R, 1, one);
  ring_scalar(R, -1, mone);
  ring_from_nfelt(R, hnf, eps, e);
  memcpy(x, one, sizeof(res_t));
  for (;;) {
    if (n + 2 > cap) return 0;
    memcpy(U[n++], x, sizeof(res_t));
    ring_mul(R, x, mone, t);
    memcpy(U[n++], t, sizeof(res_t));
    ring_mul(R, x, e, x);
    if (!memcmp(x, one, sizeof(res_t))) break;
  }
  for (i = 0, w = 0; i < n; i++) {            /* dedupe: -1 may be a power */
    for (k = 0; k < w; k++) if (!memcmp(U[i], U[k], sizeof(res_t))) break;
    if (k == w) memcpy(U[w++], U[i], sizeof(res_t));
  }
  return w;
}

/* ================================================================== */
/* candidate moduli                                                   */
/* ================================================================== */

typedef struct {
  GEN pr;              /* the prime ideal p */
  long k;              /* exponent: the candidate is p^k */
  int grp;             /* index of p; at most one power of each p is used */
  long long Q;         /* N(p^k) */
  long t;              /* multiplicative order of eps mod p^k */
} cand_t;

/* order of eps in (O_K/a)^*, capped at UMAX; 0 if larger */
static long unit_order_mod(GEN nf, GEN a, GEN eps)
{
  pari_sp av = avma;
  ring_t R;
  res_t e, x, one;
  long t;
  if (!ring_set_lattice(&R, a)) { set_avma(av); return 0; }
  if (!ring_set_multable(&R, nf)) { set_avma(av); return 0; }
  ring_scalar(&R, 1, one);
  ring_from_nfelt(&R, a, eps, e);
  memcpy(x, e, sizeof(res_t));
  for (t = 1; t <= UMAX; t++) {
    if (!memcmp(x, one, sizeof(res_t))) { set_avma(av); return t; }
    ring_mul(&R, x, e, x);
  }
  set_avma(av);
  return 0;
}

typedef struct { int n; int idx[MAXSEL]; double logNB; } subset_t;

typedef struct {
  const cand_t *c;
  int nc;
  double target;       /* need log(prod N(p^k)) >= target */
  double best;         /* best log(prod Q) so far */
  subset_t cur, bestset;
  long nodes;
} dfs_t;

/* Smallest N(b) = prod N(p^k) that reaches exp(target), using at most one
 * power of each prime.  Branch and bound over subsets. */
static void dfs(dfs_t *D, int start, double lq, double lphi)
{
  int j;
  if (D->nodes++ > 100000) return;
  for (j = start; j < D->nc; j++) {
    double nlq, nlphi;
    int i, clash = 0;
    for (i = 0; i < D->cur.n; i++)
      if (D->c[D->cur.idx[i]].grp == D->c[j].grp) { clash = 1; break; }
    if (clash) continue;
    nlq   = lq + log((double)D->c[j].Q);
    nlphi = nlq;
    if (nlq >= D->best) continue;
    D->cur.idx[D->cur.n++] = j;
    if (nlphi >= D->target) {
      D->best = nlq;
      D->bestset = D->cur;
      D->bestset.logNB = nlq;
    } else if (D->cur.n < MAXSEL) {
      dfs(D, j + 1, nlq, nlphi);
    }
    D->cur.n--;
  }
}

/* ================================================================== */
/* the marking kernel                                                 */
/* ================================================================== */

typedef struct { long long norm; res_t r; } pgen_t;

typedef struct {
  const ring_t *R;
  const pgen_t *P;
  int np;
  long long bound;      /* strict norm bound */
  const res_t *U;
  int u;
  unsigned *firsthit;
  long long nideals;
} mark_t;

static void mark_res(mark_t *M, long long n, const res_t r)
{
  int w;
  for (w = 0; w < M->u; w++) {
    res_t z;
    long long idx;
    ring_mul(M->R, r, M->U[w], z);
    idx = ring_index(M->R, z);
    if ((unsigned)n < M->firsthit[idx]) M->firsthit[idx] = (unsigned)n;
  }
}

/* every ideal coprime to b of norm < bound, as a product of the P[j] */
static void enum_ideals(mark_t *M, int start, long long norm, const res_t res)
{
  int j;
  for (j = start; j < M->np; j++) {
    long long np = M->P[j].norm, n2 = norm, lim = (M->bound - 1) / np;
    res_t r2;
    if (norm > lim) break;                    /* P is sorted by norm */
    memcpy(r2, res, sizeof(res_t));
    while (n2 <= lim) {
      n2 *= np;
      ring_mul(M->R, r2, M->P[j].r, r2);
      M->nideals++;
      mark_res(M, n2, r2);
      enum_ideals(M, j + 1, n2, r2);
    }
  }
}

/* ================================================================== */
/* cache of (prime ideal, verified generator), shared by all the moduli   */
/* tried for one field                                                    */
/* ================================================================== */

typedef struct { GEN pr, g; long long norm; } pcache_t;

static pcache_t *g_pc = NULL;
static int g_npc = 0, g_pccap = 0;
static long long g_pc_upto = 0;     /* cache is complete for norms < g_pc_upto */
static int g_pc_bad = 0;            /* a generator failed its exact check */

static int pcache_cmp(const void *a, const void *b)
{
  long long x = ((const pcache_t *)a)->norm, y = ((const pcache_t *)b)->norm;
  return (x > y) - (x < y);
}

/* Keep the cache sorted by norm.  The fill produces it in order of the
 * rational prime, which is badly out of order -- a prime of residue degree 3
 * has norm p^3 but sits beside the norm-p entries -- and run_modulus needs it
 * sorted.  Sorting here, once per fill, replaces an O(n^2) insertion sort per
 * modulus that dominated everything once the cache was pre-filled. */
static void pcache_sort(void)
{ qsort(g_pc, (size_t)g_npc, sizeof(pcache_t), pcache_cmp); }

static void pcache_reset(void)
{
  int i;
  for (i = 0; i < g_npc; i++) { gunclone(g_pc[i].pr); gunclone(g_pc[i].g); }
  g_npc = 0; g_pc_upto = 0; g_pc_bad = 0;
}

/* make the cache complete for all prime ideals of norm < bound */
/* ------------------------------------------------------------------ *
 * Parallel fill.
 *
 * Filling the table is the sequential half of the per-field cost (about 49%
 * at |N beta| = 3e6, and it is the larger share the smaller rho_K is), so for
 * the few fields that get a tier to themselves it is worth spreading over the
 * cores.  The work is embarrassingly parallel: each rational prime needs an
 * independent idealprimedec + bnfisprincipal against a read-only bnf.
 *
 * We fork rather than thread.  This PARI is built PARI_MT_ENGINE "single", and
 * even a threaded build would be delicate here: bnfisprincipal0 calls
 * bnf_build_cycgen, one of several lazy builders that obj_insert their result
 * back into the bnf, so a nominally read-only shared object is written on
 * first use.  Under fork that is harmless -- children get private
 * copy-on-write pages.
 *
 * The usual cost of fork, serialising results back, is avoided with a shared
 * anonymous mapping created before the fork: each child writes fixed-format
 * records into its own disjoint slice, so there is no IPC, no locking and no
 * PARI parsing.  A t_INT is written as (sign, word count, words) and rebuilt
 * verbatim; parent and children are the same binary, so the word order of the
 * kernel does not matter.
 *
 * The children do the expensive part (bnfisprincipal); the parent still
 * performs the exact verification (g) = p on every generator it accepts, so
 * the rigour contract is unchanged by running this in parallel.
 * ------------------------------------------------------------------ */

typedef struct {
  int      k;             /* number of workers */
  size_t   slice;         /* bytes per worker */
  size_t  *used;          /* [k], in shared memory */
  unsigned long *sum;     /* [k]: checksum of the bytes each worker wrote */
  int     *err;           /* [k]: 1 = overflow, 2 = generator check failed */
  unsigned char *buf;     /* k slices of `slice` bytes */
  void    *map;
  size_t   maplen;
} arena_t;

static int arena_init(arena_t *A, int k, size_t slice)
{
  size_t head = ((size_t)k * (sizeof(size_t) + sizeof(unsigned long) + sizeof(int))
                 + 63) & ~(size_t)63;
  A->maplen = head + (size_t)k * slice;
  A->map = mmap(NULL, A->maplen, PROT_READ | PROT_WRITE,
                MAP_SHARED | MAP_ANONYMOUS, -1, 0);
  if (A->map == MAP_FAILED) return 0;
  A->k = k; A->slice = slice;
  A->used = (size_t *)A->map;
  A->sum  = (unsigned long *)((unsigned char *)A->map + (size_t)k * sizeof(size_t));
  A->err  = (int *)((unsigned char *)A->map
                    + (size_t)k * (sizeof(size_t) + sizeof(unsigned long)));
  A->buf  = (unsigned char *)A->map + head;
  memset(A->map, 0, head);
  return 1;
}

static void arena_free(arena_t *A) { if (A->map) munmap(A->map, A->maplen); }

/* room check, then raw appends into worker i's slice */
static int arena_room(arena_t *A, int i, size_t n)
{
  if (A->used[i] + n > A->slice) { A->err[i] = 1; return 0; }
  return 1;
}
static void arena_put(arena_t *A, int i, const void *src, size_t n)
{
  memcpy(A->buf + (size_t)i * A->slice + A->used[i], src, n);
  A->used[i] += n;
}

/* FNV-1a over a worker's records.  The exact check (g) = p is done by the
 * worker that computed g; this is what lets the parent accept the result
 * without redoing the algebra, while still detecting any corruption of the
 * shared mapping between writing and reading. */
static unsigned long fnv1a(const unsigned char *p, size_t n)
{
  unsigned long h = 1469598103934665603UL;
  size_t i;
  for (i = 0; i < n; i++) { h ^= p[i]; h *= 1099511628211UL; }
  return h;
}

/* a t_INT as sign / word count / words */
static size_t int_wire_size(GEN x)
{ return 2 * sizeof(long) + (size_t)(lgefint(x) - 2) * sizeof(long); }

static void arena_put_int(arena_t *A, int i, GEN x)
{
  long sgn = (long)signe(x), n = lgefint(x) - 2, j;
  arena_put(A, i, &sgn, sizeof(long));
  arena_put(A, i, &n, sizeof(long));
  for (j = 0; j < n; j++) arena_put(A, i, x + 2 + j, sizeof(long));
}

static GEN wire_get_int(unsigned char **p)
{
  long sgn, n, j;
  GEN z;
  memcpy(&sgn, *p, sizeof(long)); *p += sizeof(long);
  memcpy(&n,   *p, sizeof(long)); *p += sizeof(long);
  if (!sgn) { *p += (size_t)n * sizeof(long); return gen_0; }
  z = (sgn > 0) ? cgetipos(n + 2) : cgetineg(n + 2);
  for (j = 0; j < n; j++) { memcpy(z + 2 + j, *p, sizeof(long)); *p += sizeof(long); }
  return z;
}

#define ARENA_MAX 2.0e9        /* bytes of shared arena for one fill pass */
#define MEM_BUDGET 8.0e9       /* bytes the search workers may use between them */

static int g_fill_jobs = 1;

static double wallclock(void)
{ struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec + 1e-9 * ts.tv_nsec; }

static int pcache_fill_parallel(GEN bnf, GEN nf, long long bound, int k);

/* returns 1 if the parallel path completed, 0 to fall back to sequential */
static int pcache_fill_parallel(GEN bnf, GEN nf, long long bound, int k)
{
  arena_t A;
  double pi_est = (double)bound / log((double)bound);
  size_t est = (size_t)(1.6 * pi_est) + 1024;
  size_t slice = (est / (size_t)k + 4096) * 208;
  long long lo = g_pc_upto;
  int i, ok = 1;
  pid_t pids[64];
  double t0 = wallclock(), t1;

  /* If the whole range would need too large an arena, split it: the fill is
   * incremental, so doing it in a few sub-ranges is exactly equivalent.  The
   * earlier version gave up and fell back to the sequential fill instead,
   * which disabled parallelism precisely at the bounds where it matters. */
  if ((double)k * (double)slice > ARENA_MAX) {
    long long mid = lo + (bound - lo) / 2;
    if (mid <= lo || mid >= bound) return 0;
    if (!pcache_fill_parallel(bnf, nf, mid, k)) return 0;
    return pcache_fill_parallel(bnf, nf, bound, k);
  }
  if (!arena_init(&A, k, slice)) return 0;

  for (i = 0; i < k; i++) {
    pid_t pid = fork();
    if (pid < 0) { for (--i; i >= 0; i--) kill(pids[i], SIGKILL); arena_free(&A); return 0; }
    if (pid == 0) {
      forprime_t T;
      long p;
      long idx = 0;
      u_forprime_init(&T, 2, bound - 1);
      while ((p = u_forprime_next(&T))) {
        pari_sp av2 = avma;
        GEN dec;
        long s;
        if ((idx++ % k) != i) continue;                 /* stride: balances */
        if ((long long)p < lo && (long long)p * p * p < lo) continue;
        dec = idealprimedec(nf, utoi(p));
        for (s = 1; s < lg(dec); s++) {
          GEN q = gel(dec, s), g;
          long long Q = 1;
          int j;
          size_t need;
          for (j = 0; j < pr_get_f(q); j++) Q *= p;
          if (Q >= bound || Q < lo) continue;
          g = gel(bnfisprincipal0(bnf, q, nf_GEN | nf_FORCE), 2);
          /* exact check that (g) = q, done here rather than in the parent:
           * |N g| = N q and g lies in q */
          if (typ(g) != t_COL ||
              !equalii(absi(nfnorm(nf, g)), utoipos((ulong)Q)) ||
              !gequal0(nf_to_Fq(nf, g, nfmodprinit(nf, q))))
            { A.err[i] = 2; goto child_done; }
          need = 2 * sizeof(long);
          for (j = 0; j < DEG; j++) need += int_wire_size(gel(g, j + 1));
          if (!arena_room(&A, i, need)) goto child_done;
          { long pp = p, ss = s;
            arena_put(&A, i, &pp, sizeof(long));
            arena_put(&A, i, &ss, sizeof(long)); }
          for (j = 0; j < DEG; j++) arena_put_int(&A, i, gel(g, j + 1));
        }
        set_avma(av2);
      }
    child_done:
      A.sum[i] = fnv1a(A.buf + (size_t)i * A.slice, A.used[i]);
      _exit(0);
    }
    pids[i] = pid;
  }
  for (i = 0; i < k; i++) waitpid(pids[i], NULL, 0);
  t1 = wallclock();

  for (i = 0; i < k; i++)
    if (A.err[i] ||
        A.sum[i] != fnv1a(A.buf + (size_t)i * A.slice, A.used[i])) ok = 0;
  if (ok) {
    GEN lastdec = NULL;
    long lastp = 0;
    for (i = 0; i < k && ok; i++) {
      unsigned char *p = A.buf + (size_t)i * A.slice, *end = p + A.used[i];
      while (p < end) {
        pari_sp av2 = avma;
        long pp, ss, j;
        GEN g, q;
        long long Q = 1;
        memcpy(&pp, p, sizeof(long)); p += sizeof(long);
        memcpy(&ss, p, sizeof(long)); p += sizeof(long);
        g = cgetg(DEG + 1, t_COL);
        for (j = 0; j < DEG; j++) gel(g, j + 1) = wire_get_int(&p);
        if (pp != lastp || !lastdec) {          /* clone: it must outlive av2 */
          if (lastdec) gunclone(lastdec);
          lastdec = gclone(idealprimedec(nf, utoi(pp)));
          lastp = pp;
        }
        if (ss < 1 || ss >= lg(lastdec)) { ok = 0; set_avma(av2); break; }
        q = gel(lastdec, ss);
        for (j = 0; j < pr_get_f(q); j++) Q *= pp;
        if (g_npc >= g_pccap) {
          g_pccap = g_pccap ? 2 * g_pccap : 4096;
          g_pc = realloc(g_pc, sizeof(pcache_t) * g_pccap);
          if (!g_pc) { fprintf(stderr, "out of memory (prime cache)\n"); exit(2); }
        }
        g_pc[g_npc].pr = gclone(q); g_pc[g_npc].g = gclone(g); g_pc[g_npc].norm = Q;
        g_npc++;
        set_avma(av2);
      }
    }
    if (lastdec) gunclone(lastdec);
    if (ok) { g_pc_upto = bound; pcache_sort(); }
  }
  DBG("      fill(%lld) x%d: workers %.2fs, merge %.2fs, %d records\n",
      bound, k, t1 - t0, wallclock() - t1, g_npc);
  arena_free(&A);
  return ok;
}

static void pcache_fill(GEN bnf, GEN nf, long long bound)
{
  forprime_t T;
  long p;
  pari_sp av = avma;
  if (bound <= g_pc_upto) return;
  if (g_fill_jobs > 1 && bound > 200000 &&
      pcache_fill_parallel(bnf, nf, bound, g_fill_jobs)) { set_avma(av); return; }
  u_forprime_init(&T, 2, bound - 1);
  while ((p = u_forprime_next(&T))) {
    pari_sp av2 = avma;      /* the clones live off-stack, so this is safe and
                              * necessary: without it the stack grows by one
                              * idealprimedec + bnfisprincipal per prime */
    GEN dec;
    long s;
    if ((long long)p < g_pc_upto && (long long)p * p * p < g_pc_upto) continue;
    dec = idealprimedec(nf, utoi(p));
    for (s = 1; s < lg(dec); s++) {
      GEN q = gel(dec, s), g;
      long long Q = 1;
      int j;
      for (j = 0; j < pr_get_f(q); j++) Q *= p;
      if (Q >= bound || Q < g_pc_upto) continue;
      g = gel(bnfisprincipal0(bnf, q, nf_GEN | nf_FORCE), 2);
      /* exact check that (g) = q: |N g| = N q and g lies in q */
      if (typ(g) != t_COL ||
          !equalii(absi(nfnorm(nf, g)), utoipos((ulong)Q)) ||
          !gequal0(nf_to_Fq(nf, g, nfmodprinit(nf, q)))) { g_pc_bad = 1; set_avma(av); return; }
      if (g_npc >= g_pccap) {
        g_pccap = g_pccap ? 2 * g_pccap : 4096;
        g_pc = realloc(g_pc, sizeof(pcache_t) * g_pccap);
        if (!g_pc) { fprintf(stderr, "out of memory (prime cache)\n"); exit(2); }
      }
      /* clone: run_modulus resets the PARI stack under our feet */
      g_pc[g_npc].pr = gclone(q); g_pc[g_npc].g = gclone(g); g_pc[g_npc].norm = Q;
      g_npc++;
    }
    set_avma(av2);
  }
  g_pc_upto = bound;
  pcache_sort();
  set_avma(av);
}

/* ================================================================== */

typedef struct {
  double safety;        /* steers the modulus search (expected unmarked orbits) */
  long pbound;          /* rational primes scanned for candidate moduli */
  long long nbcap;      /* largest allowed |N(beta)| */
  double refine;        /* enumerate ideals of norm < refine*|N beta| */
  int certify;          /* run bnfcertify (removes the GRH hypothesis) */
  double budget;        /* total |N(beta)| the field may spend */
  int jobs;             /* >1: fork this many workers over the modulus pool */
} opts_t;

typedef struct {
  const char *verdict;
  long long NB, m;
  int m_exact;
  char beta[512];      /* coordinates in the integral basis */
  long alpha[DEG];
  char beta_pow[640];  /* (b0 + b1 x + b2 x^2)/D, basis-independent */
  char alpha_pow[640];
  int u, ntried;
  double best_ratio;   /* max over all moduli tried of m/|N beta|; >= 1 iff settled */
  long long br_m, br_NB;
  long long nideals;
  double rho, reg;
  long h;
  int certified;
  double ms;
  char note[128];
} result_t;

static unsigned *g_firsthit = NULL;
static long long g_firsthit_size = 0;

static unsigned *get_firsthit(long long n)
{
  if (n > g_firsthit_size) {
    free(g_firsthit);
    g_firsthit = malloc((size_t)n * sizeof(unsigned));
    if (!g_firsthit) { fprintf(stderr, "out of memory (%lld cells)\n", n); exit(2); }
    g_firsthit_size = n;
  }
  return g_firsthit;
}

/* write x in O_K as (c0 + c1*t + c2*t^2)/D with integers c_i, D > 0, where t
 * is the root of the defining polynomial: a representation independent of the
 * integral basis, so that other systems can read the certificate */
static void write_pow(GEN nf, GEN x, char *buf, size_t n)
{
  pari_sp av = avma;
  GEN a = nf_to_scalar_or_alg(nf, x), den, num;
  char *c[DEG], *d;
  int i;
  num = Q_remove_denom(a, &den);
  if (!den) den = gen_1;
  if (typ(num) != t_POL) num = scalarpol(num, 0);
  for (i = 0; i < DEG; i++) {
    GEN t = (degpol(num) >= i) ? gel(num, i + 2) : gen_0;
    c[i] = GENtostr(t);
  }
  d = GENtostr(den);
  snprintf(buf, n, "[%s,%s,%s,%s]", c[0], c[1], c[2], d);
  for (i = 0; i < DEG; i++) pari_free(c[i]);
  pari_free(d);
  set_avma(av);
}

/* Replace beta by eps^j*beta with T2(eps^j*beta) minimal.  bnfisprincipal can
 * return a generator with coordinates of size e^{Reg}; multiplying by units
 * costs nothing and makes the certificate small enough to print and to
 * re-check by hand. */
static GEN unit_reduce(GEN nf, GEN beta, GEN eps)
{
  pari_sp av = avma;
  GEN G = nf_get_G(nf), inv = nfinv(nf, eps), cur = beta;
  GEN t = gnorml2(RgM_RgC_mul(G, cur));
  int dir, it;
  for (dir = 0; dir < 2; dir++) {
    GEN mul = dir ? inv : eps;
    for (it = 0; it < 4096; it++) {
      GEN nxt = nfmul(nf, cur, mul), tn = gnorml2(RgM_RgC_mul(G, nxt));
      if (gcmp(tn, t) >= 0) break;
      cur = nxt; t = tn;
    }
  }
  return gc_GEN(av, cur);
}

/* N(prod pr[i]^kk[i]), or 0 if it does not fit in a long long */
static long long modulus_norm(GEN nf, GEN *pr, const long *kk, int m)
{
  pari_sp av = avma;
  GEN N = gen_1;
  int i;
  for (i = 0; i < m; i++)
    N = mulii(N, powiu(pr_norm(pr[i]), (ulong)kk[i]));
  if (is_bigint(N)) { set_avma(av); return 0; }
  { long long r = itos(N); set_avma(av); return r; }
}

/* Try one modulus b = prod pr[i]^kk[i].  Returns 1 iff a witness was found. */
static int run_modulus(GEN bnf, GEN nf, GEN eps, GEN *pr, const long *kk, int m,
                       const opts_t *O, result_t *out)
{
  pari_sp av = avma;   /* NB: pcache_fill allocates above this mark, so the
                        * cache must be filled before anything is discarded */
  ring_t R;
  GEN b, beta;
  res_t U[2 * UMAX + 4];
  pgen_t *P;
  int u, i, np = 0, cap = 1024, ok = 0;
  long long bound, bestidx = -1;
  unsigned *fh, bestval = 0;
  mark_t M;

  b = idealpow(nf, pr[0], stoi(kk[0]));
  for (i = 1; i < m; i++) b = idealmul(nf, b, idealpow(nf, pr[i], stoi(kk[i])));
  b = idealhnf(nf, b);
  if (!ring_set_lattice(&R, b) || R.NB > O->nbcap) { set_avma(av); return 0; }
  if (!ring_set_multable(&R, nf)) { set_avma(av); return 0; }

  {                                            /* beta with (beta) = b */
    GEN v = bnfisprincipal0(bnf, b, nf_GEN | nf_FORCE);
    beta = gel(v, 2);
    if (typ(beta) != t_COL) { set_avma(av); return 0; }
    if (!gequal(idealhnf(nf, beta), b)) { set_avma(av); return 0; }   /* exact */
    beta = unit_reduce(nf, beta, eps);
  }

  u = unit_orbit(&R, b, eps, U, 2 * UMAX + 4);
  if (!u) { set_avma(av); return 0; }

  bound = (long long)(O->refine * (double)R.NB);
  if (bound < R.NB) bound = R.NB;
  if (bound > 2000000000LL) { set_avma(av); return 0; }  /* first-hit is u32 */

  /* residues of the generators of every prime ideal of norm < bound; the
   * primes dividing b are included, because a residue class that is a zero
   * divisor is just as good a witness -- it is only the class of 0 that has
   * to be excluded, and no gamma with |N gamma| < N(b) lies in it. */
  pcache_fill(bnf, nf, bound);
  if (g_pc_bad) {
    snprintf(out->note, sizeof(out->note), "generator check failed");
    set_avma(av); return 0;
  }
  cap = g_npc + 1;
  P = malloc(sizeof(pgen_t) * cap);
  for (i = 0; i < g_npc; i++) {
    if (g_pc[i].norm >= bound) continue;
    P[np].norm = g_pc[i].norm;
    ring_from_nfelt(&R, b, g_pc[i].g, P[np].r);
    np++;
  }
  /* P is already sorted by norm: it is a filtered copy of the cache, which
   * pcache_sort keeps in order. */

  fh = get_firsthit(R.NB);
  { long long k; for (k = 0; k < R.NB; k++) fh[k] = 0xffffffffu; }

  M.R = &R; M.P = P; M.np = np; M.bound = bound;
  M.U = (const res_t *)U; M.u = u; M.firsthit = fh; M.nideals = 1;
  { res_t one; ring_scalar(&R, 1, one);
    mark_res(&M, 1, one);                      /* the unit ideal: gamma a unit */
    enum_ideals(&M, 0, 1, one); }

  { long long k;                               /* latest-hit non-zero class */
    for (k = 1; k < R.NB; k++)                 /* k = 0 is the class of 0 */
      if (fh[k] > bestval) { bestval = fh[k]; bestidx = k; }
  }
  DBG("      NB=%lld u=%d np=%d ideals=%lld best=%u\n",
      R.NB, u, np, M.nideals, bestval);
  { long long bm = (bestval == 0xffffffffu) ? bound : (long long)bestval;
    double ratio = (double)bm / (double)R.NB;
    if (ratio > out->best_ratio) {
      out->best_ratio = ratio; out->br_m = bm; out->br_NB = R.NB;
    } }

  if (bestidx >= 0 && (long long)bestval >= R.NB) {
    res_t v;
    ring_unindex(&R, bestidx, v);
    out->verdict = "NOT_EUCLIDEAN";
    out->NB      = R.NB;
    out->m       = (bestval == 0xffffffffu) ? bound : (long long)bestval;
    out->m_exact = (bestval != 0xffffffffu);
    out->u       = u;
    out->nideals = M.nideals;
    for (i = 0; i < DEG; i++) out->alpha[i] = v[i];
    { char *s0 = GENtostr(gel(beta, 1)), *s1 = GENtostr(gel(beta, 2)),
           *s2 = GENtostr(gel(beta, 3));
      snprintf(out->beta, sizeof(out->beta), "[%s,%s,%s]", s0, s1, s2);
      pari_free(s0); pari_free(s1); pari_free(s2); }
    /* also in the power basis 1, x, x^2, so that the certificate does not
     * depend on which integral basis a given CAS happens to choose */
    write_pow(nf, beta, out->beta_pow, sizeof(out->beta_pow));
    { GEN al = mkcol3(stoi(v[0]), stoi(v[1]), stoi(v[2]));
      write_pow(nf, al, out->alpha_pow, sizeof(out->alpha_pow)); }
    ok = 1;
  }
  free(P);
  set_avma(av);
  return ok;
}

/* One shard-worker record sent back over the pipe.  result_t holds a
 * `const char *verdict`, which is meaningless in another address space, so it
 * is zeroed on the way out and restored on the way in. */
typedef struct { int idx; int found; result_t r; } shard_msg_t;

/* Shared control block for the modulus workers.  `next` is a lock-free work
 * queue: a worker takes the next index with an atomic increment rather than
 * being handed a fixed shard up front.  That matters here because the moduli
 * span six orders of magnitude in cost (|N beta| from a handful to 10^7), so
 * any static assignment leaves most workers idle behind the one that drew the
 * expensive end of the list.  Results go in the shared slots too, so there is
 * no pipe and no signalling. */
typedef struct {
  volatile long next;        /* next unclaimed modulus */
  volatile long found_idx;   /* smallest index that produced a certificate */
  volatile int  found;
} shardctl_t;

static void run_pool_forked(GEN bnf, GEN nf, GEN eps, const cand_t *cand,
                            const subset_t *pool, const int *trylist, int ntry,
                            const opts_t *O, result_t *out)
{
  int nj = O->jobs, i;
  pid_t pids[64];
  size_t hdr = (sizeof(shardctl_t) + 63) & ~(size_t)63;
  size_t maplen = hdr + (size_t)(nj > 64 ? 64 : nj) * sizeof(shard_msg_t);
  void *map;
  shardctl_t *ctl;
  shard_msg_t *slot;
  int bestidx = -1;
  result_t best;

  if (nj > ntry) nj = ntry;
  if (nj > 64) nj = 64;
  {  /* each worker allocates 4 bytes per residue class of the largest modulus
      * it may draw, so the worker count has a memory ceiling of its own */
    long long biggest = 0;
    int t;
    for (t = 0; t < ntry; t++) {
      long long nb = pool[trylist[t]].logNB > 0 ? (long long)exp(pool[trylist[t]].logNB) : 0;
      if (nb > biggest) biggest = nb;
    }
    if (biggest > 0) {
      int lim = (int)(MEM_BUDGET / (4.0 * (double)biggest));
      if (lim < 1) lim = 1;
      if (lim < nj) {
        DBG("   pool: %d workers -> %d (memory: %.1f GB of first-hit array each)\n",
            nj, lim, 4.0 * biggest / 1e9);
        nj = lim;
      }
    }
  }
  maplen = hdr + (size_t)nj * sizeof(shard_msg_t);
  map = mmap(NULL, maplen, PROT_READ | PROT_WRITE,
             MAP_SHARED | MAP_ANONYMOUS, -1, 0);
  if (map == MAP_FAILED) nj = 1;
  ctl  = (shardctl_t *)map;
  slot = (shard_msg_t *)((unsigned char *)map + hdr);
  memset(map, 0, maplen);
  ctl->found_idx = ntry;

  for (i = 0; i < nj; i++) {
    pid_t p = fork();
    if (p < 0) { nj = i; break; }
    if (p == 0) {
      result_t r = *out;
      slot[i].found = 0; slot[i].idx = -1;
      for (;;) {
        long k = __sync_fetch_and_add(&ctl->next, 1);
        GEN pr[MAXSEL];
        long ex[MAXSEL];
        int pi, kk;
        if (k >= ntry) break;
        /* someone already has a certificate from an earlier modulus */
        if (ctl->found && ctl->found_idx <= k) break;
        pi = trylist[k];
        for (kk = 0; kk < pool[pi].n; kk++) {
          pr[kk] = cand[pool[pi].idx[kk]].pr;
          ex[kk] = cand[pool[pi].idx[kk]].k;
        }
        if (run_modulus(bnf, nf, eps, pr, ex, pool[pi].n, O, &r)) {
          long cur;
          slot[i].found = 1; slot[i].idx = (int)k;
          do {                                   /* keep the least index */
            cur = ctl->found_idx;
            if (cur <= k) break;
          } while (!__sync_bool_compare_and_swap(&ctl->found_idx, cur, k));
          ctl->found = 1;
          break;
        }
      }
      r.verdict = NULL;
      slot[i].r = r;
      _exit(0);
    }
    pids[i] = p;
  }
  for (i = 0; i < nj; i++) waitpid(pids[i], NULL, 0);

  for (i = 0; i < nj; i++) {
    if (slot[i].r.best_ratio > out->best_ratio) {
      out->best_ratio = slot[i].r.best_ratio;
      out->br_m = slot[i].r.br_m; out->br_NB = slot[i].r.br_NB;
    }
    if (slot[i].found && (bestidx < 0 || slot[i].idx < bestidx)) {
      bestidx = slot[i].idx; best = slot[i].r;
    }
  }
  if (bestidx >= 0) {
    double br = out->best_ratio;
    long long bm = out->br_m, bn = out->br_NB;
    *out = best;
    out->verdict = "NOT_EUCLIDEAN";
    if (br > out->best_ratio) { out->best_ratio = br; out->br_m = bm; out->br_NB = bn; }
  }
  DBG("   pool: %d workers drew %ld of %d moduli%s\n",
      nj, ctl->next < ntry ? ctl->next : (long)ntry, ntry,
      bestidx >= 0 ? " (stopped early on a certificate)" : "");
  munmap(map, maplen);
}

static void solve_field(GEN pol, const opts_t *O, result_t *out)
{
  pari_sp av = avma;
  GEN bnf, nf, eps;
  double rho;
  cand_t cand[MAXCAND];
  int nc = 0, ngrp = 0;
  clock_t t0 = clock();
  double t_start = wallclock();

  memset(out, 0, sizeof(*out));
  out->verdict = "SURVIVOR";
  pcache_reset();

  bnf = Buchall(pol, nf_FORCE, DEFAULTPREC);
  nf  = bnf_get_nf(bnf);
  out->h   = itos(bnf_get_no(bnf));
  out->reg = rtodbl(bnf_get_reg(bnf));
  if (O->certify) out->certified = (bnfcertify(bnf) == 1);

  if (out->h != 1) {          /* norm-Euclidean => PID, so h > 1 already settles it */
    out->verdict = "NOT_EUCLIDEAN";
    snprintf(out->note, sizeof(out->note), "class number %ld > 1", out->h);
    out->ms = 1000.0 * (clock() - t0) / CLOCKS_PER_SEC;
    set_avma(av);
    return;
  }

  eps = algtobasis(nf, gel(bnf_get_fu(bnf), 1));
  rho = 2.0 * M_PI * out->reg / sqrt(fabs(gtodouble(nf_get_disc(nf))));
  out->rho = rho;

  /* ---- candidate prime powers p^k with small ord(eps mod p^k) ----
   *
   * ord(eps mod p) divides t exactly when p divides eps^t - 1, so instead of
   * scanning the rational primes we factor the integers |N(eps^t - 1)| for a
   * few small t.  Trial division up to pbound costs one pass over a big
   * integer rather than pbound calls to idealprimedec, and it reaches primes
   * far larger than a scan could -- which is what fields with a large
   * regulator (large rho) need. */
  {
    static const long TVALS[] = {1, 2, 3, 4, 6, 8, 12};
    ulong *qs = malloc(sizeof(ulong) * MAXCAND);
    int nq = 0, iq, it;
    forprime_t T;
    long p;
    /* the smallest primes are always worth having, factorisation or not */
    u_forprime_init(&T, 2, O->pbound < 200 ? O->pbound : 200);
    while ((p = u_forprime_next(&T)) && nq < MAXCAND) qs[nq++] = (ulong)p;
    for (it = 0; it < (int)(sizeof(TVALS) / sizeof(long)) && nq < MAXCAND; it++) {
      pari_sp av2 = avma;
      GEN e1 = nfadd(nf, nfpow(nf, eps, stoi(TVALS[it])), gen_m1);
      GEN N, F;
      long r;
      if (gequal0(e1)) { set_avma(av2); continue; }
      N = absi(nfnorm(nf, e1));
      if (!signe(N)) { set_avma(av2); continue; }
      F = Z_factor_limit(N, (ulong)O->pbound);
      for (r = 1; r <= nbrows(F) && nq < MAXCAND; r++) {
        GEN q = gcoeff(F, r, 1);
        ulong qq;
        int j, seen = 0;
        if (is_bigint(q) || cmpiu(q, (ulong)O->pbound) > 0) continue;
        qq = itou(q);
        for (j = 0; j < nq; j++) if (qs[j] == qq) { seen = 1; break; }
        if (!seen) qs[nq++] = qq;
      }
      set_avma(av2);
    }
    for (iq = 0; iq < nq && nc < MAXCAND; iq++) {
      pari_sp av2 = avma;
      GEN dec = idealprimedec(nf, utoi(qs[iq]));
      long s;
      for (s = 1; s < lg(dec) && nc < MAXCAND; s++) {
        GEN q = gel(dec, s);
        long long Q = 1, Qk;
        long k, j;
        int used = 0;
        for (j = 0; j < pr_get_f(q); j++) Q *= (long long)qs[iq];
        if (Q > O->nbcap) continue;
        for (k = 1, Qk = Q; Qk <= O->nbcap && nc < MAXCAND; k++, Qk *= Q) {
          long t = unit_order_mod(nf, idealhnf(nf, idealpow(nf, q, stoi(k))), eps);
          if (!t) break;
          cand[nc].pr  = gclone(q);
          cand[nc].k   = k;
          cand[nc].grp = ngrp;
          cand[nc].Q   = Qk;
          cand[nc].t   = t;
          nc++; used = 1;
        }
        if (used) ngrp++;
      }
      set_avma(av2);
    }
    free(qs);
  }
  if (g_debug) {
    int i;
    fprintf(stderr, "  rho=%.4f  candidates:", rho);
    for (i = 0; i < nc; i++)
      fprintf(stderr, " (N=%lld,k=%ld,t=%ld)", cand[i].Q, cand[i].k, cand[i].t);
    fprintf(stderr, "\n");
  }

  /* ---- assemble a pool of moduli and try them, cheapest first ---- */
  {
    static const int TARGETS[] = {1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 24};
    const int nT = (int)(sizeof(TARGETS) / sizeof(int));
    subset_t pool[4096];
    int npool = 0, si, ti, k, kk;
    double budget = 0.0;

    for (si = 0; si < 18; si++) {
      double safety = O->safety * pow(2.0, si - 6);
      for (ti = 0; ti < 2 * nT; ti++) {
        int Tt = TARGETS[ti / 2], i, ne = 0, map[MAXCAND], forced;
        cand_t elig[MAXCAND];
        /* u = |<-1,eps> mod b| is lcm(t_i) or twice it; try both guesses */
        double U = (ti & 1) ? 2.0 * Tt : (double)Tt;
        double target = log(safety) + log(U) + rho * U;
        if (rho * U > 40.0) continue;
        for (i = 0; i < nc; i++)
          if (Tt % cand[i].t == 0) { map[ne] = i; elig[ne++] = cand[i]; }
        if (!ne) continue;
        for (i = 1; i < ne; i++) {             /* sort eligible by norm */
          cand_t ct = elig[i]; int mt = map[i], j = i - 1;
          while (j >= 0 && elig[j].Q > ct.Q) { elig[j+1] = elig[j]; map[j+1] = map[j]; j--; }
          elig[j+1] = ct; map[j+1] = mt;
        }
        /* Minimising N(b) alone never picks up a prime with phi(p^k) small
         * (e.g. a degree-1 prime over 2), yet such a factor reshuffles the
         * residue classes at almost no statistical cost and is exactly what
         * the classical moduli beta | eps^t - 1 contain.  So we also try
         * moduli in which the few smallest eligible prime powers are forced. */
        for (forced = 0; forced <= 3 * 3 + 2; forced++) {
          /* `forced/3` smallest eligible prime powers are pre-selected, with
           * `forced%3` controlling which of several equal-norm candidates is
           * taken -- distinct prime ideals over the same rational prime give
           * genuinely different moduli of the same cost. */
          dfs_t D;
          double lq0 = 0.0;
          int nf = 0, want = forced / 3, off = forced % 3;
          if (want == 0 && off > 0) continue;
          D.c = elig; D.nc = ne; D.target = target;
          D.best = log((double)O->nbcap);
          D.cur.n = 0; D.bestset.n = 0; D.nodes = 0;
          for (i = off; i < ne && nf < want; i++) {
            int j, clash = 0;
            for (j = 0; j < D.cur.n; j++)
              if (elig[D.cur.idx[j]].grp == elig[i].grp) { clash = 1; break; }
            if (clash) continue;
            D.cur.idx[D.cur.n++] = i;
            lq0 += log((double)elig[i].Q);
            nf++;
          }
          if (nf < want) continue;
          if (lq0 >= D.best) continue;
          dfs(&D, 0, lq0, lq0);
          if (!D.bestset.n) continue;
          for (k = 0; k < D.bestset.n; k++) D.bestset.idx[k] = map[D.bestset.idx[k]];
          for (kk = 0; kk < npool; kk++)
            if (pool[kk].n == D.bestset.n &&
                !memcmp(pool[kk].idx, D.bestset.idx, sizeof(int) * D.bestset.n)) break;
          if (kk == npool && npool < (int)(sizeof(pool) / sizeof(pool[0])))
            pool[npool++] = D.bestset;
        }
      }
    }
    for (k = 1; k < npool; k++) {              /* cheapest first */
      subset_t t = pool[k];
      int j = k - 1;
      while (j >= 0 && pool[j].logNB > t.logNB) { pool[j + 1] = pool[j]; j--; }
      pool[j + 1] = t;
    }

    /* the moduli we are actually going to try, with their exact norms */
    { static int trylist[512];
      static long long tryNB[512];
      int ntry = 0;
      long long maxbound = 0;
      for (k = 0; k < npool && ntry < 512; k++) {
        GEN pr[MAXSEL];
        long ex[MAXSEL];
        long long NB;
        for (kk = 0; kk < pool[k].n; kk++) {
          pr[kk] = cand[pool[k].idx[kk]].pr;
          ex[kk] = cand[pool[k].idx[kk]].k;
        }
        NB = modulus_norm(nf, pr, ex, pool[k].n);
        if (!NB || NB > O->nbcap) continue;
        if (budget + (double)NB > O->budget) break;
        budget += (double)NB;
        trylist[ntry] = k; tryNB[ntry] = NB; ntry++;
        { long long bd = (long long)(O->refine * (double)NB);
          if (bd > maxbound) maxbound = bd; }
      }
      out->ntried = ntry;

      if (O->jobs > 1 && ntry > 1) {
        /* Staged fill.  trylist is sorted by |N beta|, so instead of filling
         * the generator table to the largest bound any modulus might need and
         * only then trying them, we advance in geometric stages: fill to the
         * bound the cheap moduli need, try those, extend, try the next band.
         * A field that settles on an early modulus never pays for the big
         * fill, and one that exhausts everything pays only
         * 1 + 1/F + 1/F^2 + ... = F/(F-1) times a single full fill.  The
         * search is bounded by --nbcap and --budget, both intrinsic: the same
         * invocation does the same work and returns the same certificate on
         * any machine.  Within a
         * stage the workers draw moduli from a shared queue (run_pool_forked);
         * fork means the extended cache cannot be shared backwards, so the
         * stages are separated by a join. */
        const long long STAGE_F = 8, STAGE_MIN = 50000;
        int start = 0;
        while (start < ntry) {
          long long limit = (long long)(O->refine * (double)tryNB[start]) * STAGE_F;
          long long sb = 0;
          int end = start;
          if (limit < STAGE_MIN) limit = STAGE_MIN;
          /* Once few enough moduli remain to schedule in one go, stop staging
           * and hand the whole tail to the queue.  Splitting the expensive end
           * into bands makes the makespan the *sum* of each band's slowest
           * modulus instead of a single balanced schedule, which costs more
           * than the fill it saves. */
          if (ntry - start <= 2 * O->jobs) limit = (long long)1 << 62;
          while (end < ntry) {
            long long b = (long long)(O->refine * (double)tryNB[end]);
            if (b < tryNB[end]) b = tryNB[end];
            if (b > limit && end > start) break;
            if (b > sb) sb = b;
            end++;
          }
          DBG("   stage: moduli [%d,%d) of %d, fill to %lld  (%.1fs elapsed)\n",
              start, end, ntry, sb, wallclock() - t_start);
          pcache_fill(bnf, nf, sb);
          run_pool_forked(bnf, nf, eps, cand, pool, trylist + start,
                          end - start, O, out);
          if (!strcmp(out->verdict, "NOT_EUCLIDEAN")) break;
          start = end;
        }
        (void)maxbound;
      } else {
        for (k = 0; k < ntry; k++) {
          GEN pr[MAXSEL];
          long ex[MAXSEL];
          int pi = trylist[k];
          DBG("   try #%d: m=%d NB=%lld\n", pi, pool[pi].n, tryNB[k]);
          for (kk = 0; kk < pool[pi].n; kk++) {
            pr[kk] = cand[pool[pi].idx[kk]].pr;
            ex[kk] = cand[pool[pi].idx[kk]].k;
          }
          if (run_modulus(bnf, nf, eps, pr, ex, pool[pi].n, O, out)) break;
        }
      }
    }
  }

  { int i; for (i = 0; i < nc; i++) gunclone(cand[i].pr); }
  out->ms = (O->jobs > 1) ? 1000.0 * (wallclock() - t_start)
                          : 1000.0 * (clock() - t0) / CLOCKS_PER_SEC;
  set_avma(av);
}

/* ================================================================== */

static long long gcdll(long long a, long long b)
{ while (b) { long long t = a % b; a = b; b = t; } return a; }

static void usage(void)
{
  fprintf(stderr,
    "usage: ne_sweep [options] < fields.csv > results.jsonl\n"
    "  --safety S    steers the modulus search (default 3)\n"
    "  --pbound P    trial-division limit when looking for moduli (default 100000)\n"
    "  --nbcap N     largest allowed |N(beta)| (default 200000)\n"
    "  --budget B    total |N(beta)| spent per field (default 6*nbcap)\n"
    "  --jobs N      parallelise ONE field N ways: forked workers both build\n"
    "                the prime-generator table and try the modulus pool (default 1)\n"
    "  --refine R    enumerate ideals of norm < R*|N(beta)| (default 1)\n"
    "  --no-certify  skip bnfcertify: results become GRH-conditional\n"
    "  --skip K --take N   process only input rows K..K+N-1\n"
    "  --debug       trace the modulus search on stderr\n");
  exit(1);
}

int main(int argc, char **argv)
{
  opts_t O = { 3.0, 100000, 200000LL, 1.0, 1, -1.0, 1 };
  volatile long skip = 0, take = -1, row = 0, done = 0;
  char line[512];
  int i;

  for (i = 1; i < argc; i++) {
    if      (!strcmp(argv[i], "--safety") && i + 1 < argc) O.safety = atof(argv[++i]);
    else if (!strcmp(argv[i], "--pbound") && i + 1 < argc) O.pbound = atol(argv[++i]);
    else if (!strcmp(argv[i], "--nbcap")  && i + 1 < argc) O.nbcap  = atoll(argv[++i]);
    else if (!strcmp(argv[i], "--budget") && i + 1 < argc) O.budget = atof(argv[++i]);
    else if (!strcmp(argv[i], "--refine") && i + 1 < argc) O.refine = atof(argv[++i]);
    else if (!strcmp(argv[i], "--jobs")   && i + 1 < argc) O.jobs   = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--no-certify")) O.certify = 0;
    else if (!strcmp(argv[i], "--debug")) g_debug = 1;
    else if (!strcmp(argv[i], "--skip") && i + 1 < argc) skip = atol(argv[++i]);
    else if (!strcmp(argv[i], "--take") && i + 1 < argc) take = atol(argv[++i]);
    else usage();
  }
  if (O.budget < 0) O.budget = 6.0 * (double)O.nbcap;
  g_fill_jobs = O.jobs;   /* --jobs parallelises the fill as well as the search */

  pari_init(1 << 26, 500000);            /* 64 MB, grows on demand */
  paristack_setsize(1 << 26, 1L << 31);  /* hard ceiling 2 GB per worker */

  while (fgets(line, sizeof(line), stdin)) {
    char label[64];
    long long disc;
    long a[4];
    result_t r;
    GEN pol;
    if (line[0] == 'l') continue;                         /* csv header */
    if (sscanf(line, "%63[^,],%lld,%ld,%ld,%ld,%ld",
               label, &disc, &a[0], &a[1], &a[2], &a[3]) != 6) continue;
    if (row++ < skip) continue;
    if (take >= 0 && done >= take) break;
    done++;

    pol = mkpoln(4, stoi(a[3]), stoi(a[2]), stoi(a[1]), stoi(a[0]));
    setvarn(pol, 0);
    pari_CATCH(CATCH_ALL) {                    /* never let one field kill the run */
      GEN err = pari_err_last();
      char *msg = pari_err2str(err);
      memset(&r, 0, sizeof(r));
      r.verdict = "PARTIAL";
      snprintf(r.note, sizeof(r.note), "pari error: %.100s", msg ? msg : "?");
      pari_free(msg);
      pcache_reset();
    } pari_TRY {
      solve_field(pol, &O, &r);
    } pari_ENDCATCH;

    printf("{\"label\":\"%s\",\"disc\":%lld,\"poly\":[%ld,%ld,%ld,%ld],"
           "\"verdict\":\"%s\",\"h\":%ld,\"reg\":%.10g,\"rho\":%.6g,"
           "\"certified\":%s,\"ntried\":%d,\"runtime_ms\":%.3f",
           label, disc, a[0], a[1], a[2], a[3], r.verdict, r.h,
           r.reg, r.rho, r.certified ? "true" : "false", r.ntried, r.ms);
    if (r.br_NB) {
      long long g2 = gcdll(r.br_m, r.br_NB);
      printf(",\"best_ratio\":%.6f,\"best_bound\":\"%lld/%lld\"",
             r.best_ratio, r.br_m / g2, r.br_NB / g2);
    }
    if (r.NB) {
      long long g = gcdll(r.m, r.NB);
      printf(",\"Nbeta\":%lld,\"m\":%lld,\"m_exact\":%s,"
             "\"beta\":%s,\"alpha\":[%ld,%ld,%ld],"
             "\"beta_pow\":%s,\"alpha_pow\":%s,"
             "\"orbit\":%d,\"nideals\":%lld,\"M_lower_bound\":\"%lld/%lld\"",
             r.NB, r.m, r.m_exact ? "true" : "false",
             r.beta,
             r.alpha[0], r.alpha[1], r.alpha[2],
             r.beta_pow, r.alpha_pow,
             r.u, r.nideals, r.m / g, r.NB / g);
    }
    if (r.note[0]) printf(",\"note\":\"%s\"", r.note);
    printf("}\n");
    fflush(stdout);
  }
  return 0;
}
