/* euclid.cpp -- rigorous proof that a complex cubic field IS norm-Euclidean.
 *
 * The converse direction to src/certify.  A certificate there exhibits one
 * point with M(K,xi) >= 1; here we must show no such point exists, i.e. that
 * every xi has some gamma in O_K with |N(xi - gamma)| < 1.
 *
 * The method is Cavallar and Lemmermeyer's: cover a fundamental domain of
 * (K tensor R)/O_K by cubes, discard any cube on which a single gamma keeps
 * |N| below k, and interleave that with an elimination step using the unit
 * action.  Their implementation evaluated the cube bound in floating point
 * with a 1% margin "as a protection against rounding errors", and says of the
 * result only that "experiments have led us to trust our results".  A discard
 * decided by a rounded comparison is not a proof, and a wrongly discarded cube
 * would hide an exceptional point and so declare a field norm-Euclidean that
 * is not.  Here every bound is an enclosure and a cube is discarded only when
 * the upper end of that enclosure is below k, so insufficient precision costs
 * a subdivision, never a false discard.
 *
 * ---- interval arithmetic: three builds from this one file --------------
 *
 * The covering bound is the whole inner loop.  -DIALIB selects how it is
 * evaluated; everything else is shared.
 *
 *   -DIALIB=ARB    (default) Arb balls throughout.  Slower, but needs no
 *                  unusual build flags and never touches the SSE control word.
 *                  A separately written implementation of the same algorithm,
 *                  so it checks the *algorithm*.
 *   -DIALIB=PLATT  Platt's int_double -- hardware doubles with SSE directed
 *                  rounding.  Several times faster.  Needs -msse2
 *                  -mfpmath=sse and src/int_double14.2.h.
 *   -DIALIB=MPFI   True inf/sup interval arithmetic at IA_PREC bits (53 by
 *                  default).  Shares the *same source* as PLATT -- only the
 *                  interval type differs -- so the sequence of operations is
 *                  identical by construction.  Interval +, * and abs each have
 *                  a unique tightest representable enclosure, so at 53 bits a
 *                  correct MPFI and a correct Platt must agree bit for bit:
 *                  any difference is a defect in the SSE path, which is the
 *                  one that produced the shipped results.  Raising IA_PREC
 *                  does NOT tighten the bound: the enclosure widths here are
 *                  set by the width of the cube being bounded, not by
 *                  rounding, and 53, 113 and 256 bits give identical output.
 *
 * Everything else -- the algebra, the offset set, the elimination, Theorem 1,
 * the order in which gamma are tried -- is shared, so the three builds differ
 * only in how three affine forms are multiplied together.  They are expected
 * to return the same verdict for every field, and do.  PLATT and MPFI agree
 * cube for cube as well; ARB subdivides slightly further on 28 of the 58 and
 * never less, its balls being symmetric about a midpoint and so unable to
 * track an asymmetric interval as tightly as an inf/sup pair.
 *
 * The program is single threaded.  Fields are independent and vastly unequal
 * in cost, so parallelism belongs in the driver, which runs one process per
 * core longest-job-first; see src/run_euclid.py.  Threading inside a field
 * bought nothing once the offset set and the mean-value bound had removed the
 * work it was hiding.
 *
 * Two things here are not optimisations but requirements -- without either the
 * computation does not terminate at all, rather than merely running slowly:
 *
 *   - the covering bound is a mean-value form (see norm_bound).  The naive
 *     product of the three embeddings' enclosures has a first-order dependency
 *     error that never shrinks relative to the cube, so refinement does not
 *     converge;
 *   - the elimination walks a precomputed offset set (see build_offsets)
 *     rather than the image bounding box, and that set is found by flood fill
 *     rather than by scanning the box.
 *
 * The embeddings are not taken on trust from a package: the roots of f are
 * enclosed by arb_fmpz_poly_complex_roots, as in src/certify.
 *
 * Usage:  euclid --nf FIELD.dat [--k 1.0] [--depth 20] [--gbox 8] [-v]
 *                [--selftest]
 *
 * Exit status 0 if the field is proved norm-Euclidean, 1 otherwise.
 */
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpz_poly.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_vec.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/arb_fmpz_poly.h>
/* ---- interval backend ------------------------------------------------- */
#define ARB   1
#define MPFI  2
#define PLATT 3
#ifndef IALIB
#define IALIB ARB
#endif
#if IALIB != ARB && IALIB != MPFI && IALIB != PLATT
#error "IALIB must be one of ARB, MPFI, PLATT"
#endif

#if IALIB == PLATT
#include <int_double14.2.h>
typedef int_double ival;

/* Keep control-register safeguards in our wrapper, not in Platt's header. */
static inline void read_mxcsr(unsigned int &control)
{
  asm volatile ("stmxcsr %0" : "=m" (control) : : "memory");
}

/* Round-down for the duration of norm_bound only.  Platt sets round-down once
 * and globally, which is right for a program that is nothing but interval
 * arithmetic; here it is a small part of a process that spends most of its
 * time in FLINT, Arb and the C++ library, and leaving those to run under a
 * rounding mode they never agreed to is not safe. */
/* The opaque read/write pointer operands make the kernel's input loads depend
 * on setting MXCSR.  The read/write result operand makes restoration depend
 * on the completed bound, and the caller's use of it depend on restoration.
 * Memory clobbers also order the interval operators' memory operands.  Merely
 * marking ldmxcsr volatile would not stop arithmetic moving across it. */
#define IA_RND_SAVE(cube, gamma) unsigned int _cw_saved, _cw_down;             \
                        read_mxcsr(_cw_saved);                              \
                        _cw_down = (_cw_saved & 0x1FBF) | 0x2000;            \
                        asm volatile ("ldmxcsr %2"                           \
                            : "+r" (cube), "+r" (gamma) : "m" (_cw_down)    \
                            : "memory")
#define IA_RND_RESTORE(bound) asm volatile ("ldmxcsr %1"                      \
                            : "+m" (bound) : "m" (_cw_saved) : "memory")

#elif IALIB == MPFI
#include "mpfi.h"
#ifndef IA_PREC
#define IA_PREC 53                 /* 53 reproduces the double-precision path */
#endif
#define IA_RND_SAVE(cube, gamma) ((void)0)
#define IA_RND_RESTORE(bound)   ((void)0)

/* An inf/sup interval with int_double's interface, so that norm_bound below
 * compiles unchanged against either.  `right` holds the *negated* upper
 * endpoint, which is Platt's convention (it lets both endpoints round the same
 * way under a single SSE mode); reproducing it here is what keeps the shared
 * source free of backend conditionals.
 *
 * The mpfi_t carries the real value at IA_PREC bits; left/right are outward-
 * rounded double snapshots refreshed after each operation.  They are read only
 * to extract a final bound or to test the sign, never fed back into the
 * arithmetic, so no precision is lost along a chain of operations. */
class ival {
public:
  mpfi_t v;
  double left, right;
  void sync() {
    mpfr_t t; mpfr_init2(t, IA_PREC);
    mpfi_get_left(t, v);   left  =  mpfr_get_d(t, MPFR_RNDD);
    mpfi_get_right(t, v);  right = -mpfr_get_d(t, MPFR_RNDU);
    mpfr_clear(t);
  }
  ival()                    { mpfi_init2(v, IA_PREC); mpfi_set_ui(v, 0); sync(); }
  ival(double l)            { mpfi_init2(v, IA_PREC); mpfi_set_d(v, l);  sync(); }
  ival(double l, double r)  { mpfi_init2(v, IA_PREC); mpfi_interv_d(v, l, r); sync(); }
  ival(const ival &o)       { mpfi_init2(v, IA_PREC); mpfi_set(v, o.v);
                              left = o.left; right = o.right; }
  ival &operator=(const ival &o) {
    if (this != &o) { mpfi_set(v, o.v); left = o.left; right = o.right; }
    return *this;
  }
  ~ival() { mpfi_clear(v); }
};
static inline ival operator+(const ival &a, const ival &b)
  { ival r; mpfi_add(r.v, a.v, b.v); r.sync(); return r; }
static inline ival operator*(const ival &a, const ival &b)
  { ival r; mpfi_mul(r.v, a.v, b.v); r.sync(); return r; }
static inline ival abs(const ival &a)
  { ival r; mpfi_abs(r.v, a.v); r.sync(); return r; }

#endif  /* IALIB backend type */

#include <array>
#include <algorithm>
#include <set>
#include <string>
#include <vector>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cmath>

static const int D = 3;
/* 64 bits is ample: every quantity here is O(1) and the cubes reach 2^-20.
 * Precision carries no rigour burden in any case -- an enclosure too wide to
 * decide a cube costs a subdivision, never a false discard. */
static const slong PREC = 64;

struct Cube { long x, y, z; int d; };      /* the box [x/2^d,(x+1)/2^d] etc. */
typedef std::array<long,3> Key;

static std::vector<std::string> split(const std::string &s, char c)
{
  std::vector<std::string> v; size_t p = 0;
  while (true) {
    size_t q = s.find(c, p);
    v.push_back(s.substr(p, q == std::string::npos ? q : q - p));
    if (q == std::string::npos) break;
    p = q + 1;
  }
  return v;
}

/* ---- the order ------------------------------------------------------- */

/* multiplication table of the order, built exactly from the integral basis;
 * needed only to get the matrix by which eps acts on coordinates. */
static fmpz *TAB;
static fmpz *cc(int i, int j, int l) { return TAB + ((i*D + j)*D + l); }

static int build_table(const fmpz_poly_t f, const fmpz_mat_t B, const fmpz_t den)
{
  int ok = 1;
  fmpq_mat_t Bq, Binv;
  fmpq_mat_init(Bq, D, D); fmpq_mat_init(Binv, D, D);
  for (int j = 0; j < D; j++)
    for (int l = 0; l < D; l++)
      fmpq_set_fmpz_frac(fmpq_mat_entry(Bq, l, j), fmpz_mat_entry(B, j, l), den);
  if (!fmpq_mat_inv(Binv, Bq)) { fmpq_mat_clear(Bq); fmpq_mat_clear(Binv); return 0; }
  for (int i = 0; i < D && ok; i++)
    for (int j = 0; j < D && ok; j++) {
      fmpq_poly_t pi, pj, pr, fq;
      fmpq_poly_init(pi); fmpq_poly_init(pj); fmpq_poly_init(pr); fmpq_poly_init(fq);
      fmpq_poly_set_fmpz_poly(fq, f);
      for (int l = 0; l < D; l++) {
        fmpq_t t; fmpq_init(t);
        fmpq_set_fmpz_frac(t, fmpz_mat_entry(B, i, l), den);
        fmpq_poly_set_coeff_fmpq(pi, l, t);
        fmpq_set_fmpz_frac(t, fmpz_mat_entry(B, j, l), den);
        fmpq_poly_set_coeff_fmpq(pj, l, t);
        fmpq_clear(t);
      }
      fmpq_poly_mul(pr, pi, pj);
      fmpq_poly_rem(pr, pr, fq);
      fmpq_mat_t v, rr; fmpq_mat_init(v, D, 1); fmpq_mat_init(rr, D, 1);
      for (int l = 0; l < D; l++) {
        fmpq_t t; fmpq_init(t); fmpq_poly_get_coeff_fmpq(t, pr, l);
        fmpq_set(fmpq_mat_entry(v, l, 0), t); fmpq_clear(t);
      }
      fmpq_mat_mul(rr, Binv, v);
      for (int l = 0; l < D; l++) {
        if (!fmpz_is_one(fmpq_denref(fmpq_mat_entry(rr, l, 0)))) { ok = 0; break; }
        fmpz_set(cc(i, j, l), fmpq_numref(fmpq_mat_entry(rr, l, 0)));
      }
      fmpq_mat_clear(v); fmpq_mat_clear(rr);
      fmpq_poly_clear(pi); fmpq_poly_clear(pj); fmpq_poly_clear(pr); fmpq_poly_clear(fq);
    }
  fmpq_mat_clear(Bq); fmpq_mat_clear(Binv);
  return ok;
}

/* ---- the offset set --------------------------------------------------- */

/* The image of the grid cell c + [0,1]^3 under an integer matrix M is
 *
 *     M*c + M*[0,1]^3,
 *
 * and M*c is an *integer* vector.  So every cube's image is a translate, by a
 * lattice vector, of the one fixed parallelepiped P = M*[0,1]^3: the set of
 * cells the image meets is the same for every cube.  Computing it once turns
 * the elimination's per-cube test into a walk over OFF instead of a scan of
 * the image bounding box.  For M the matrix of a unit, P has volume
 * |det M| = 1 -- a single cell's worth -- but the expanding embedding shears
 * it into a long thin sliver, so it meets about the *sum* of its extents in
 * cells while its bounding box holds their *product*.
 *
 * OFF must be a superset of the cells P really meets.  An extra cell can only
 * make a test report "meets" when it need not have, which leaves a cube in
 * play, and that is the safe direction; a missing cell would permit a wrong
 * elimination, i.e. a false proof.  So the test is exact integer arithmetic
 * and closed (touching counts as meeting), with no floating point at all.
 *
 * Exactness is the separating-axis theorem, which is complete for two convex
 * polytopes on the face normals of each together with the pairwise cross
 * products of their edge directions.  The cell has normals and edges e_l; P
 * has edges the columns of M and, since |det M| = 1 makes the inverse
 * integral, normals the rows of Minv.  Every axis and projection is an
 * integer.
 */
static void build_offsets(long M[D][D], long Minv[D][D], std::vector<Key> &OFF)
{
  long lo0[D], hi0[D];
  for (int l = 0; l < D; l++) {
    lo0[l] = hi0[l] = 0;
    for (int j = 0; j < D; j++)
      { if (M[l][j] < 0) lo0[l] += M[l][j]; else hi0[l] += M[l][j]; }
  }
  std::vector<Key> axes;
  for (int l = 0; l < D; l++) {
    Key n = {0,0,0}; n[l] = 1; axes.push_back(n);
    Key m = { Minv[l][0], Minv[l][1], Minv[l][2] }; axes.push_back(m);
  }
  for (int l = 0; l < D; l++)
    for (int j = 0; j < D; j++) {
      long v[D] = { M[0][j], M[1][j], M[2][j] };          /* column j */
      Key n;                                              /* e_l x v  */
      n[0] = (l==1 ? -v[2] : (l==2 ?  v[1] : 0));
      n[1] = (l==2 ? -v[0] : (l==0 ?  v[2] : 0));
      n[2] = (l==0 ? -v[1] : (l==1 ?  v[0] : 0));
      if (n[0] || n[1] || n[2]) axes.push_back(n);
    }
  std::vector<long> plo(axes.size()), phi(axes.size()),
                    blo(axes.size()), bhi(axes.size());
  for (size_t a = 0; a < axes.size(); a++) {
    long pl = 0, ph = 0, bl = 0, bh = 0;
    for (int j = 0; j < D; j++) {
      long w = 0;                                      /* w = (n^T M)_j */
      for (int l = 0; l < D; l++) w += axes[a][l] * M[l][j];
      if (w < 0) pl += w; else ph += w;
      if (axes[a][j] < 0) bl += axes[a][j]; else bh += axes[a][j];
    }
    plo[a] = pl; phi[a] = ph; blo[a] = bl; bhi[a] = bh;
  }
  /* Walk the cells the parallelepiped meets, rather than testing every cell of
   * its bounding box.  The box is the very thing the offset set exists to
   * avoid: for eps it is small, but for eps^{-1} the shear runs the other way
   * and it reaches 10^9 cells, which cost 39 seconds of pure startup before
   * this was flood-filled instead.
   *
   * A flood fill is complete here because the cells meeting a convex body are
   * connected under the 26-neighbourhood: given two of them, pick a point of
   * the body in each; the segment between lies in the body, consecutive cells
   * along it share at least a point, and every cell it passes through meets
   * the body.  So starting from the cell holding the centre of P and expanding
   * over neighbours reaches every cell and no others.  The test applied to
   * each candidate is the same exact one; only the enumeration changed.
   */
  OFF.clear();
  long seed[D];
  for (int l = 0; l < D; l++) {
    long ssum = 0;                       /* 2 * centre of P, in coordinate l */
    for (int j = 0; j < D; j++) ssum += M[l][j];
    seed[l] = (ssum >= 0) ? ssum/2 : -(((-ssum)+1)/2);   /* floor(ssum/2) */
  }
  std::set<Key> seen;
  std::vector<Key> stack;
  { Key z = {seed[0],seed[1],seed[2]}; seen.insert(z); stack.push_back(z); }
  while (!stack.empty()) {
    Key z = stack.back(); stack.pop_back();
    OFF.push_back(z);
    for (int dx = -1; dx <= 1; dx++)
    for (int dy = -1; dy <= 1; dy++)
    for (int dz = -1; dz <= 1; dz++) {
      if (!dx && !dy && !dz) continue;
      Key w = { z[0]+dx, z[1]+dy, z[2]+dz };
      if (seen.count(w)) continue;
      bool sep = false;
      for (size_t a = 0; a < axes.size() && !sep; a++) {
        long d = axes[a][0]*w[0] + axes[a][1]*w[1] + axes[a][2]*w[2];
        /* the cell projects to [d+blo, d+bhi], P to [plo, phi]; strict, so
         * that a touching pair counts as meeting */
        if (phi[a] < d + blo[a] || d + bhi[a] < plo[a]) sep = true;
      }
      seen.insert(w);
      if (!sep) stack.push_back(w);
    }
  }
  (void)lo0; (void)hi0;
}

/* ---- the covering bound ----------------------------------------------- */

/* Paper notation: W1[j] = (omega_j)_R, W2[j] = (omega_j)_C.
 * For xi - gamma, the local names L,R,I (or s1,sr,si) denote z,x,y,
 * respectively.  The cube coordinates t_j are coefficients in the integral
 * basis; they are not the embedding coordinates x,y,z. */
static arb_t W1[D];     /* (omega_j)_R, the real embedding    */
static acb_t W2[D];     /* (omega_j)_C, one complex embedding */

/* Plain-double copies of the embeddings, and the Lipschitz constants of the
 * two places.  These drive the prefilter below and the self-test; neither
 * carries any rigour burden. */
static double w1d[D], w2rd[D], w2id[D], lip1, lip2;

/* A cheap rejection test for a candidate gamma, in ordinary floating point.
 *
 * Most gamma are hopeless for most cubes, and an uncoverable cube has to try
 * every one of them; deciding each with an enclosure is the bulk of the run.
 * This settles the hopeless ones with a handful of doubles instead: it lower
 * bounds the absolute values of the real and complex embeddings over the cube by their values at the centre
 * less the Lipschitz movement across it, and rejects when even that lower
 * bound is far above k.
 *
 * It is deliberately generous, and it is the one place where floating point is
 * free: a wrong rejection here can only skip a gamma that might have covered
 * the cube, which costs a subdivision.  It can never discard a cube, so it
 * cannot make a false proof.  A wrong acceptance costs nothing at all -- the
 * enclosure is then evaluated and has the final word. */
static bool plausible(const Cube &C, const long *g, double k)
{
  double s1 = 0, sr = 0, si = 0;
  double h = ldexp(1.0, -(C.d + 1));
  for (int j = 0; j < D; j++) {
    double cj = ldexp((double)(2*((j==0)?C.x:((j==1)?C.y:C.z)) + 1), -(C.d+1)) - g[j];
    s1 += w1d[j]*cj; sr += w2rd[j]*cj; si += w2id[j]*cj;
  }
  double b1 = fabs(s1) - h*lip1;                 if (b1 < 0) b1 = 0;
  double m2 = sqrt(sr*sr + si*si) - h*lip2;      if (m2 < 0) m2 = 0;
  return b1 * m2 * m2 < 4.0 * k;                 /* only prunes the hopeless */
}

#if IALIB != ARB
/* Enclosures of the embeddings as machine intervals, derived from the Arb
 * enclosures (themselves from certified roots) by rounding outward.  Both
 * non-Arb backends start from exactly these doubles, so PLATT and MPFI are
 * given identical input as well as an identical operation sequence. */
static ival W1I[D], W2RI[D], W2II[D];
#endif

#if IALIB == ARB

/* A rigorous upper bound for sup |N(xi - gamma)| over the cube, as an Arb
 * ball reduced to a double rounded up.  Returns -1 if unusable.
 *
 * With xi - gamma = sum_j t_j omega_j and L,R,I = z,x,y the three
 * embedding coordinates, |N| = (R^2 + I^2)|L| = (x^2+y^2)|z|, and these
 * coordinates are affine in t.  Two enclosures are computed and the smaller
 * returned; both are valid.
 *
 * (1) The naive one, |L| evaluated on the box times (R^2 + I^2) evaluated on
 *     the box.  Each of L, R, I is an *exact* enclosure, since each coordinate
 *     occurs once in each; the loss is entirely in multiplying them as though
 *     they were independent when they are three functions of the same point.
 *     That error is first order in the width of the cube, so it shrinks no
 *     faster than the cube does and defeats the same proportion of cubes at
 *     every depth -- subdivision then never converges.  It is worst exactly
 *     where the problem is hardest: near a maximiser of M(K,.) the gradient of
 *     |N| nearly vanishes, so the true variation over a cube is second order
 *     while this bound still charges first order for it.
 *
 * (2) The mean-value form.  For t = m + delta over the cube, |delta_j| <= r,
 *
 *       f(t) = f(m) + sum_j delta_j (df/dt_j)(chi)   for some chi in the cube,
 *
 *     so f <= f(m) + r sum_j sup |df/dt_j|, every sup taken as an enclosure.
 *     The overestimate is now second order.  The derivative
 *
 *       df/dt_j = sgn(L) W1_j (R^2+I^2) + |L| 2 (R W2R_j + I W2I_j)
 *
 *     must be enclosed in one piece: bounding its two terms separately and
 *     adding magnitudes reproduces (1) exactly, because the entire gain is the
 *     cancellation between the variation of |L| and that of R^2 + I^2.
 *
 *     |L| is not differentiable where L vanishes, so this needs L to keep one
 *     sign on the cube.  Where it does not, (1) stands -- and there |N| is
 *     near zero anyway, which is the easy case to cover.
 *
 * The cube's endpoints and centre are exact dyadics, so they enter the
 * arithmetic without error; all widening is Arb's own.
 */
static double norm_bound(const Cube &C, const long *g)
{
  arb_t L, R, I, Lm, Rm, Im, q, t, acc, dfj, fm, r, one;
  arb_t cbox[D], cmid[D];
  acb_t z, u;
  arb_init(L); arb_init(R); arb_init(I);
  arb_init(Lm); arb_init(Rm); arb_init(Im);
  arb_init(q); arb_init(t); arb_init(acc); arb_init(dfj);
  arb_init(fm); arb_init(r); arb_init(one);
  acb_init(z); acb_init(u);
  for (int j = 0; j < D; j++) { arb_init(cbox[j]); arb_init(cmid[j]); }

  long xs[D] = { C.x, C.y, C.z };
  for (int j = 0; j < D; j++) {
    arb_t lo, hi;
    arb_init(lo); arb_init(hi);
    arb_set_si(lo, xs[j]);     arb_mul_2exp_si(lo, lo, -C.d);
    arb_set_si(hi, xs[j] + 1); arb_mul_2exp_si(hi, hi, -C.d);
    arb_sub_si(lo, lo, g[j], PREC);
    arb_sub_si(hi, hi, g[j], PREC);
    arb_union(cbox[j], lo, hi, PREC);            /* the whole cube  */
    /* the centre (2x+1)/2^(d+1) - g, again an exact dyadic */
    arb_set_si(cmid[j], 2*xs[j] + 1);
    arb_mul_2exp_si(cmid[j], cmid[j], -(C.d + 1));
    arb_sub_si(cmid[j], cmid[j], g[j], PREC);
    arb_clear(lo); arb_clear(hi);
  }

  /* the three embeddings, over the cube and at its centre */
  arb_zero(L); arb_zero(R); arb_zero(I);
  arb_zero(Lm); arb_zero(Rm); arb_zero(Im);
  for (int j = 0; j < D; j++) {
    arb_mul(t, W1[j], cbox[j], PREC); arb_add(L, L, t, PREC);
    arb_mul(t, acb_realref(W2[j]), cbox[j], PREC); arb_add(R, R, t, PREC);
    arb_mul(t, acb_imagref(W2[j]), cbox[j], PREC); arb_add(I, I, t, PREC);
    arb_mul(t, W1[j], cmid[j], PREC); arb_add(Lm, Lm, t, PREC);
    arb_mul(t, acb_realref(W2[j]), cmid[j], PREC); arb_add(Rm, Rm, t, PREC);
    arb_mul(t, acb_imagref(W2[j]), cmid[j], PREC); arb_add(Im, Im, t, PREC);
  }

  /* (1) the naive product */
  arb_sqr(q, R, PREC); arb_sqr(t, I, PREC); arb_add(q, q, t, PREC);
  arb_abs(t, L);
  arb_mul(t, t, q, PREC);
  arf_t f; arf_init(f);
  arb_get_ubound_arf(f, t, PREC);
  double ub = arf_get_d(f, ARF_RND_UP);

  /* (2) the mean-value form, when L keeps one sign on the cube */
  if (!arb_contains_zero(L)) {
    int neg = arb_is_negative(L);
    arb_set_si(r, 1); arb_mul_2exp_si(r, r, -(C.d + 1));   /* half-width */
    arb_abs(t, L);                                          /* |L| on the box */
    arb_zero(acc);
    for (int j = 0; j < D; j++) {
      /* dfj = sgn(L) W1_j q + |L| 2 (R W2R_j + I W2I_j), in one piece */
      arb_mul(dfj, W1[j], q, PREC);
      if (neg) arb_neg(dfj, dfj);
      arb_t s; arb_init(s);
      arb_mul(s, R, acb_realref(W2[j]), PREC);
      arb_mul(one, I, acb_imagref(W2[j]), PREC);
      arb_add(s, s, one, PREC);
      arb_mul_2exp_si(s, s, 1);
      arb_mul(s, s, t, PREC);
      arb_add(dfj, dfj, s, PREC);
      arb_abs(dfj, dfj);
      arb_add(acc, acc, dfj, PREC);
      arb_clear(s);
    }
    /* f(m) + r * acc */
    arb_sqr(fm, Rm, PREC); arb_sqr(one, Im, PREC); arb_add(fm, fm, one, PREC);
    arb_abs(one, Lm);
    arb_mul(fm, fm, one, PREC);
    arb_mul(acc, acc, r, PREC);
    arb_add(fm, fm, acc, PREC);
    arb_get_ubound_arf(f, fm, PREC);
    double ubc = arf_get_d(f, ARF_RND_UP);
    if (ubc == ubc && ubc < ub) ub = ubc;
  }
  arf_clear(f);

  arb_clear(L); arb_clear(R); arb_clear(I);
  arb_clear(Lm); arb_clear(Rm); arb_clear(Im);
  arb_clear(q); arb_clear(t); arb_clear(acc); arb_clear(dfj);
  arb_clear(fm); arb_clear(r); arb_clear(one);
  acb_clear(z); acb_clear(u);
  for (int j = 0; j < D; j++) { arb_clear(cbox[j]); arb_clear(cmid[j]); }
  if (!(ub == ub)) return -1.0;
  return ub;
}

#else  /* PLATT or MPFI -- one source, two interval types */

#if IALIB == PLATT
/* Initialise the sign masks Platt's multiply/divide asm uses.  This is his
 * _fpu_rndd() minus the crlibm-derived log constants, which nothing here
 * needs; omitting the mask initialisation silently corrupts every product, so
 * it is done explicitly rather than relying on the header.
 *
 * It deliberately does NOT touch the SSE control word.  Platt sets round-down
 * once and globally, which is right for a program that is nothing but interval
 * arithmetic; here the interval arithmetic is a small part of a process that
 * spends most of its time in FLINT, Arb and the C++ library, and leaving those
 * to run with the rounding mode in a state they never agreed to is not safe --
 * any of them may use floating point internally.  The mode is instead set and
 * restored inside bound_double, the only place it is needed.  The ival
 * constructors do not need it: they only negate, which is exact. */
static void platt_init()
{
  d_zero_zero.left = 0.0;     d_zero_zero.right = 0.0;
  d_zero.left = 0.0;          d_zero.right = _nze;
  d_neg_zero.left = _nze;     d_neg_zero.right = _nze;
  d_neg_neg_zero.left = _nze; d_neg_neg_zero.right = 0.0;
}
#else
static void platt_init() { }        /* nothing to initialise for MPFI */
#endif

/* A rigorous upper bound for sup|N| over the cube, in hardware interval
 * arithmetic.
 *
 * Arb is doing multiprecision work on quantities that are all O(1) here, and
 * once the elimination step was made fast that arithmetic became the dominant
 * cost.  Platt's ival does the same job with SSE directed rounding, so
 * the result is still an enclosure, not an estimate -- no error analysis of
 * mine stands behind it.  The coordinate endpoints are exact dyadics, so the
 * only widening is the outward rounding of each operation.
 */
static double norm_bound(const Cube &input_cube, const long *input_gamma)
{
  /* Round-down for the duration of this function only; see platt_init.  The
   * pair of stmxcsr/ldmxcsr is not free, but it is the price of not leaving
   * every other library in the process running under a rounding mode it does
   * not expect. */
  const Cube *cube = &input_cube;
  const long *g = input_gamma;
  IA_RND_SAVE(cube, g);
  const Cube &C = *cube;

  ival s1(0.0), sr(0.0), si(0.0);     /* over the whole cube */
  ival m1(0.0), mr(0.0), mi(0.0);     /* at its centre        */
  long xs[D] = { C.x, C.y, C.z };
  double r = ldexp(1.0, -(C.d + 1));        /* half-width, exact    */
  for (int j = 0; j < D; j++) {
    double lo = ldexp((double)xs[j], -C.d) - g[j];
    double hi = ldexp((double)(xs[j] + 1), -C.d) - g[j];
    ival c(lo, hi);                    /* exact: dyadic endpoints */
    s1 = s1 + c * W1I[j];
    sr = sr + c * W2RI[j];
    si = si + c * W2II[j];
    /* the midpoint (2x+1)/2^(d+1) - g is again an exact dyadic */
    ival cm(ldexp((double)(2 * xs[j] + 1), -(C.d + 1)) - g[j]);
    m1 = m1 + cm * W1I[j];
    mr = mr + cm * W2RI[j];
    mi = mi + cm * W2II[j];
  }
  ival n = abs(s1) * (sr * sr + si * si);
  double ub = -n.right;                       /* the upper endpoint */
  if (!(ub == ub)) { IA_RND_RESTORE(ub); return -1.0; }

  /* The naive product above is the weak step.  s1, sr and si are each an
   * *exact* interval -- every coordinate occurs once in each -- but combining
   * them as independent intervals throws away the fact that they are the three
   * embeddings of one and the same point.  That dependency error is first
   * order in the cube's width, so it shrinks no faster than the cube does and
   * the fraction of cubes it defeats stays constant however deep we go.  It is
   * worst exactly where it hurts: near a maximiser of M(K,.) the gradient of
   * |N| very nearly vanishes, so the true variation over a cube is second
   * order while this bound still charges first order for it.
   *
   * The mean-value form charges the true first-order term instead.  Writing
   * f = (R^2+I^2)|L| with L,R,I = z,x,y affine in the integral-basis
   * coordinates t, for t = m + delta in the cube
   *
   *     f(t) = f(m) + sum_j delta_j (d f / d t_j)(chi),  |delta_j| <= r,
   *
   * for some chi in the cube, so f <= f(m) + r * sum_j sup |d f/d t_j|,
   * with every sup taken in interval arithmetic over the cube.  Cancellation
   * between the variation of |L| and that of R^2+I^2 is now kept rather than
   * discarded, and the overestimate drops to second order.
   *
   * f(m) is evaluated at the exact dyadic centre, so it carries only rounding
   * width.  Both bounds are valid enclosures; return whichever is smaller.
   *
   * |L| is not differentiable where L vanishes, so the mean-value step needs
   * L to keep one sign on the cube.  When it does not, the naive bound stands
   * -- and there |N| is near zero anyway, which is the easy case to cover.
   */
  double slo = s1.left, shi = -s1.right;
  if (slo > 0.0 || shi < 0.0) {
    ival q = sr * sr + si * si;
    ival L = abs(s1);
    ival sgn = (slo > 0.0) ? ival(1.0) : ival(-1.0);
    ival acc(0.0);
    for (int j = 0; j < D; j++) {
      /* Keep the derivative in one piece.  Bounding its two terms separately
       * and adding magnitudes would reproduce the naive product exactly: the
       * whole point is that d|L|/dt_j and d(R^2+I^2)/dt_j tend to oppose each
       * other near a maximiser, and that cancellation must survive. */
      ival dfj = sgn * W1I[j] * q
                     + L * (sr * W2RI[j] + si * W2II[j]) * 2.0;
      acc = acc + abs(dfj);
    }
    ival nc = abs(m1) * (mr * mr + mi * mi) + acc * r;
    double ubc = -nc.right;
    if (ubc == ubc && ubc < ub) ub = ubc;
  }
  IA_RND_RESTORE(ub);
  return ub;
}

#endif  /* IALIB */

/* ---- main ------------------------------------------------------------- */

int main(int argc, char **argv)
{
  const char *nfpath = NULL;
  double k = 1.0;
  int maxdepth = 20, gbox = 8;
  bool verbose = false, selftest = false;
  for (int i = 1; i < argc; i++) {
    if      (!strcmp(argv[i], "--nf")    && i+1 < argc) nfpath   = argv[++i];
    else if (!strcmp(argv[i], "--k")     && i+1 < argc) k        = atof(argv[++i]);
    else if (!strcmp(argv[i], "--depth") && i+1 < argc) maxdepth = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--gbox")  && i+1 < argc) gbox     = atoi(argv[++i]);
    else if (!strcmp(argv[i], "-v"))                    verbose  = true;
    else if (!strcmp(argv[i], "--selftest"))            selftest = true;
    else { fprintf(stderr, "unknown option: %s\n", argv[i]); return 2; }
  }
  if (!nfpath) { fprintf(stderr, "usage: euclid_ref --nf FIELD.dat\n"); return 2; }
  /* A proof program must not turn a bad option into a mathematical verdict.
   * --depth -1 skipped the loop with `residual` still empty and fell straight
   * through to PROVED, and --k > 1 proves only M(K) < k, which is not a proof
   * of norm-Euclideanity, yet returned the same status and exit code. */
  if (maxdepth < 0 || maxdepth > 40) {
    fprintf(stderr, "--depth must be in [0, 40]\n"); return 2; }
  if (gbox < 0 || gbox > 64) {
    fprintf(stderr, "--gbox must be in [0, 64]\n"); return 2; }
  if (!(k > 0.0) || !(k <= 1.0)) {
    fprintf(stderr, "--k must satisfy 0 < k <= 1: success is meant to prove "
                    "norm-Euclideanity, and k > 1 would not\n"); return 2; }

  fmpz_poly_t f; fmpz_mat_t B; fmpz_t den; fmpz *EPS = _fmpz_vec_init(D);
  fmpz_poly_init(f); fmpz_mat_init(B, D, D); fmpz_init(den);
  long disc = 0; bool have_eps = false;
  {
    FILE *fp = fopen(nfpath, "r");
    if (!fp) { fprintf(stderr, "cannot open %s\n", nfpath); return 2; }
    static char buf[1 << 16];
    if (!fgets(buf, sizeof(buf), fp)) { fclose(fp); return 2; }
    fclose(fp);
    std::vector<std::string> fld = split(std::string(buf), '|');
    if (fld.size() < 5) { fprintf(stderr, "malformed record\n"); return 2; }
    disc = atol(fld[0].c_str());
    /* R4: validate before indexing.  Same class of defect as C4 in
     * certify.cpp -- split() returns whatever the line held, so a short record
     * was undefined behaviour, and an unchecked fmpz_set_str would let the
     * proof run on unintended algebraic input. */
    std::vector<std::string> co = split(fld[1], ',');
    std::vector<std::string> bb = split(fld[3], ',');
    if (co.size() < 4 || bb.size() < (size_t)(D * D)) {
      fprintf(stderr, "field record has too few entries\n"); return 2; }
    for (int i = 0; i < 4; i++) {
      fmpz_t t; fmpz_init(t);
      if (fmpz_set_str(t, co[i].c_str(), 10) != 0)
        { fprintf(stderr, "bad integer in f\n"); fmpz_clear(t); return 2; }
      fmpz_poly_set_coeff_fmpz(f, i, t); fmpz_clear(t);
    }
    { fmpz_t lc; fmpz_init(lc); fmpz_poly_get_coeff_fmpz(lc, f, D);
      if (!fmpz_is_one(lc))
        { fprintf(stderr, "f is not monic\n"); fmpz_clear(lc); return 2; }
      fmpz_clear(lc); }
    if (fmpz_set_str(den, fld[2].c_str(), 10) != 0 || fmpz_is_zero(den))
      { fprintf(stderr, "bad or zero denominator\n"); return 2; }
    for (int j = 0; j < D; j++)
      for (int l = 0; l < D; l++)
        if (fmpz_set_str(fmpz_mat_entry(B, j, l), bb[3*j + l].c_str(), 10) != 0)
          { fprintf(stderr, "bad integer in basis\n"); return 2; }
    std::vector<std::string> ee = split(fld[4], ',');
    if (ee.size() >= 3) {
      for (int i = 0; i < D; i++)
        if (fmpz_set_str(EPS + i, ee[i].c_str(), 10) != 0)
          { fprintf(stderr, "bad integer in eps\n"); return 2; }
      have_eps = true;
    }
  }
  TAB = _fmpz_vec_init(D*D*D);
  if (!build_table(f, B, den)) { fprintf(stderr, "basis is not an order\n"); return 2; }

  /* Certified root enclosures, then the embeddings of the basis by Horner.
   * Nothing here is an approximation: arb_fmpz_poly_complex_roots returns
   * enclosures it has proved isolate the roots. */
  arb_t rr; acb_t ww;
  arb_init(rr); acb_init(ww);
  {
    acb_ptr rts = _acb_vec_init(D);
    int nr = 0;
    arb_fmpz_poly_complex_roots(rts, f, 0, PREC);
    for (int i = 0; i < D; i++)
      if (arb_contains_zero(acb_imagref(rts + i)))
        { arb_set(rr, acb_realref(rts + i)); nr++; }
    if (nr != 1) { fprintf(stderr, "cannot isolate the real root\n"); return 2; }
    for (int i = 0; i < D; i++)
      if (arb_is_positive(acb_imagref(rts + i))) acb_set(ww, rts + i);
    _acb_vec_clear(rts, D);
  }
  for (int j = 0; j < D; j++) {
    arb_init(W1[j]); acb_init(W2[j]);
    arb_t p, c, dd; acb_t qq, cq;
    arb_init(p); arb_init(c); arb_init(dd); acb_init(qq); acb_init(cq);
    arb_set_fmpz(dd, den);
    arb_zero(p); acb_zero(qq);
    for (int l = D-1; l >= 0; l--) {
      arb_set_fmpz(c, fmpz_mat_entry(B, j, l));
      acb_set_arb(cq, c);
      arb_mul(p, p, rr, PREC); arb_add(p, p, c, PREC);
      acb_mul(qq, qq, ww, PREC); acb_add(qq, qq, cq, PREC);
    }
    arb_div(W1[j], p, dd, PREC);
    acb_div_arb(W2[j], qq, dd, PREC);
    arb_clear(p); arb_clear(c); arb_clear(dd); acb_clear(qq); acb_clear(cq);
  }
  for (int j = 0; j < D; j++) {
    w1d[j]  = arf_get_d(arb_midref(W1[j]), ARF_RND_NEAR);
    w2rd[j] = arf_get_d(arb_midref(acb_realref(W2[j])), ARF_RND_NEAR);
    w2id[j] = arf_get_d(arb_midref(acb_imagref(W2[j])), ARF_RND_NEAR);
  }
  lip1 = lip2 = 0;
  for (int j = 0; j < D; j++) {
    lip1 += fabs(w1d[j]);
    lip2 += sqrt(w2rd[j]*w2rd[j] + w2id[j]*w2id[j]);
  }
#if IALIB != ARB
  platt_init();
  for (int j = 0; j < D; j++) {
    arf_t lo, hi; arf_init(lo); arf_init(hi);
    /* FLOOR/CEIL, not DOWN/UP: Arb's DOWN and UP are toward and away from
     * zero, which inverts the endpoints of a negative interval. */
    arb_get_interval_arf(lo, hi, W1[j], PREC);
    W1I[j]  = ival(arf_get_d(lo, ARF_RND_FLOOR), arf_get_d(hi, ARF_RND_CEIL));
    arb_get_interval_arf(lo, hi, acb_realref(W2[j]), PREC);
    W2RI[j] = ival(arf_get_d(lo, ARF_RND_FLOOR), arf_get_d(hi, ARF_RND_CEIL));
    arb_get_interval_arf(lo, hi, acb_imagref(W2[j]), PREC);
    W2II[j] = ival(arf_get_d(lo, ARF_RND_FLOOR), arf_get_d(hi, ARF_RND_CEIL));
    arf_clear(lo); arf_clear(hi);
  }
#endif

  /* the matrix of eps on coordinates, and its inverse (integral, |det| = 1) */
  long ME[D][D], MEI[D][D];
  std::vector<Key> OFF, OFFI;
  /* fmpz_get_si is not a checked narrowing: an out-of-range value is truncated
   * silently.  A wrong matrix entry can drop a cell from the offset set or send
   * a cube to the wrong image cell, and unit elimination would then discard a
   * cube still containing an exceptional point -- a false proof.  Every entry
   * must fit, and be small enough that the products and sums formed later
   * cannot overflow: coordinates run to 2^maxdepth, each image coordinate is a
   * sum of D such products, and the offsets add at most one more cell. */
  long entry_cap;
  {
    /* Two constraints bound the entries, and the tighter one wins.
     *
     * (a) Depth arithmetic: an image coordinate is a sum of D products of an
     *     entry with a grid coordinate below 2^maxdepth, plus an offset.
     * (b) build_offsets is *quadratic* in the entries, which the first version
     *     of this cap missed (Codex R3).  Its separating axes are rows of
     *     M^{-1} and cross products of columns of M, so an axis entry is O(E);
     *     w = (n^T M)_j then sums D products of two entry-sized values, the
     *     projections sum D of those, and d = n . cell pairs an axis entry with
     *     a cell coordinate of size O(D*E).  All of it stays below 10 E^2. */
    const long lim = (long)1 << 62;
    long a = lim;
    for (int sh = 0; sh <= maxdepth; sh++) a /= 2;
    a /= (D + 1);
    long lo = 1, hi = 3037000499L;          /* largest E with 10E^2 <= 2^62 */
    while (lo < hi) {
      long mid = lo + (hi - lo + 1) / 2;
      if (mid <= lim / 10 / mid) lo = mid; else hi = mid - 1;
    }
    entry_cap = (a < lo) ? a : lo;
    if (entry_cap < 1) entry_cap = 1;
  }
  /* Compare on the fmpz, never on a narrowed value: fmpz_fits_si admits
   * LONG_MIN, for which labs() is itself undefined (Codex R3). */
  fmpz_t ENTRY_CAP; fmpz_init(ENTRY_CAP); fmpz_set_si(ENTRY_CAP, entry_cap);
  if (have_eps) {
    fmpz_mat_t M; fmpz_mat_init(M, D, D);
    for (int l = 0; l < D; l++)
      for (int j = 0; j < D; j++) {
        fmpz_zero(fmpz_mat_entry(M, l, j));
        for (int i = 0; i < D; i++)
          fmpz_addmul(fmpz_mat_entry(M, l, j), EPS + i, cc(i, j, l));
      }
    for (int l = 0; l < D; l++)
      for (int j = 0; j < D; j++) {
        const fmpz *e = fmpz_mat_entry(M, l, j);
        if (fmpz_cmpabs(e, ENTRY_CAP) > 0) {
          fprintf(stderr, "eps matrix entry exceeds the exact machine range "
                          "for depth %d\n", maxdepth);
          return 2;
        }
        ME[l][j] = fmpz_get_si(e);
      }
    fmpq_mat_t A, Ai; fmpq_mat_init(A, D, D); fmpq_mat_init(Ai, D, D);
    fmpq_mat_set_fmpz_mat(A, M);
    if (!fmpq_mat_inv(Ai, A)) { fprintf(stderr, "eps not invertible\n"); return 2; }
    for (int l = 0; l < D; l++)
      for (int j = 0; j < D; j++) {
        if (!fmpz_is_one(fmpq_denref(fmpq_mat_entry(Ai, l, j))))
          { fprintf(stderr, "eps^{-1} not integral\n"); return 2; }
        const fmpz *e = fmpq_numref(fmpq_mat_entry(Ai, l, j));
        if (fmpz_cmpabs(e, ENTRY_CAP) > 0) {
          fprintf(stderr, "eps^{-1} matrix entry exceeds the exact machine "
                          "range for depth %d\n", maxdepth);
          return 2;
        }
        MEI[l][j] = fmpz_get_si(e);
      }
    fmpq_mat_clear(A); fmpq_mat_clear(Ai); fmpz_mat_clear(M);
    build_offsets(ME, MEI, OFF);
    build_offsets(MEI, ME, OFFI);
    if (verbose)
      printf("   offsets: eps image meets %ld cells, eps^{-1} image %ld\n",
             (long)OFF.size(), (long)OFFI.size());
  }

  /* ---- cover and eliminate, interleaved, depth by depth --------------- */

  if (selftest) {
    /* bound_double is now the mean-value form, and every proof rests on it
     * being a genuine upper bound for sup|N| over the cube.  Check it the only
     * way that does not just re-run the same algebra: evaluate |N| directly at
     * many points of the cube in plain double and confirm the bound dominates
     * every one of them.  Sampling can only ever find a violation, never prove
     * its absence, but a bound that is wrong by a first-order amount -- the
     * failure mode that matters -- is caught immediately. */
    long tested = 0, viol = 0; double worst = 1e300;
    unsigned long seed = 12345;
    for (int dpt = 1; dpt <= 16; dpt++)
      for (int trial = 0; trial < 400; trial++) {
        seed = seed*6364136223846793005ULL + 1442695040888963407ULL;
        long N2 = 1L << dpt;
        Cube C; C.d = dpt;
        C.x = (long)((seed>>17) % (unsigned long)N2);
        C.y = (long)((seed>>29) % (unsigned long)N2);
        C.z = (long)((seed>>41) % (unsigned long)N2);
        long g[D];
        for (int j = 0; j < D; j++)
          { seed = seed*6364136223846793005ULL + 1442695040888963407ULL;
            g[j] = (long)((seed>>33) % 17UL) - 8; }
        double u = norm_bound(C, g);
        if (!(u == u) || u < 0) continue;
        double mx = -1e300;
        for (int sm = 0; sm < 3000; sm++) {
          double s1 = 0, sr = 0, si = 0;
          for (int j = 0; j < D; j++) {
            seed = seed*6364136223846793005ULL + 1442695040888963407ULL;
            double f = (double)((seed>>11) & 0xFFFFF) / (double)0x100000;
            if (sm < 8) f = (sm >> j) & 1;            /* the corners too */
            long xj = (j==0)?C.x:((j==1)?C.y:C.z);
            double cj = ldexp((double)xj + f, -C.d) - g[j];
            s1 += w1d[j]*cj; sr += w2rd[j]*cj; si += w2id[j]*cj;
          }
          double v2 = fabs(s1)*(sr*sr+si*si);
          if (v2 > mx) mx = v2;
        }
        tested++;
        if (mx > u) viol++;
        double slack = (u - mx) / (u > 1e-12 ? u : 1e-12);
        if (slack < worst) worst = slack;
      }
    printf("SELFTEST  %ld cube/gamma pairs x 3000 sampled points:"
           " %ld violations, tightest slack %.3e\n", tested, viol, worst);
    if (viol) { printf("SELFTEST FAILED\n"); return 3; }
    return 0;
  }

  /* The gamma to try, ordered by size.  The order is a pure efficiency
   * device -- it changes which gamma is found first, never which cubes are
   * coverable -- but without it each cube pays for all (2*gbox+1)^3 of them,
   * and the natural covering gamma for a cube of the unit fundamental domain
   * is small.  Paired with remembering the last gamma that worked, this is
   * what keeps the reference version merely slow rather than unusable. */
  std::vector<Key> gammas;
  for (long a = -gbox; a <= gbox; a++)
  for (long b = -gbox; b <= gbox; b++)
  for (long c2 = -gbox; c2 <= gbox; c2++)
    { Key g = {a,b,c2}; gammas.push_back(g); }
  struct BySize {
    bool operator()(const Key &p, const Key &q) const {
      long dp = p[0]*p[0]+p[1]*p[1]+p[2]*p[2];
      long dq = q[0]*q[0]+q[1]*q[1]+q[2]*q[2];
      return dp != dq ? dp < dq : p < q;
    }
  };
  std::sort(gammas.begin(), gammas.end(), BySize());
  Key lastgood = {0,0,0};

  std::vector<Cube> live, residual;
  live.push_back((Cube){0,0,0,0});
  long ncubes = 0, ncovered = 0;
  int dep = 0, deepest = 0;
  clock_t t0 = clock();

  while (!live.empty() && dep <= maxdepth) {
    std::vector<Cube> keep;
    for (size_t ci = 0; ci < live.size(); ci++) {
      const Cube &C = live[ci];
      ncubes++;
      bool done = false;
      { long g[D] = { lastgood[0], lastgood[1], lastgood[2] };
        double u = norm_bound(C, g);
        if (u >= 0.0 && u < k) done = true; }
      for (size_t gi = 0; gi < gammas.size() && !done; gi++) {
        long g[D] = { gammas[gi][0], gammas[gi][1], gammas[gi][2] };
        if (!plausible(C, g, k)) continue;
        double u = norm_bound(C, g);
        if (u >= 0.0 && u < k) { done = true; lastgood = gammas[gi]; }
      }
      if (done) ncovered++; else keep.push_back(C);
    }
    long nafter_cover = (long)keep.size();

    /* Eliminate by the unit action, to a fixed point.  |N| is constant on
     * eps-orbits, so a cube whose image eps*T meets nothing still in play
     * cannot hold an exceptional point and can be dropped.  Dropping cubes
     * lets others become droppable, hence the iteration. */
    if (!keep.empty() && have_eps) {
      long N = 1L << dep;
      std::set<Key> res;
      for (size_t i = 0; i < keep.size(); i++)
        { Key q = { keep[i].x, keep[i].y, keep[i].z }; res.insert(q); }
      bool changed = true;
      while (changed) {
        changed = false;
        std::vector<Key> gone;
        for (std::set<Key>::iterator it = res.begin(); it != res.end(); ++it) {
          long cv[D] = { (*it)[0], (*it)[1], (*it)[2] };
          long base[D];
          for (int l = 0; l < D; l++) {
            long t = 0;
            for (int j = 0; j < D; j++) t += ME[l][j] * cv[j];
            base[l] = t;
          }
          bool meets = false;
          for (size_t oi = 0; oi < OFF.size() && !meets; oi++) {
            Key p = { ((base[0]+OFF[oi][0]) % N + N) % N,
                      ((base[1]+OFF[oi][1]) % N + N) % N,
                      ((base[2]+OFF[oi][2]) % N + N) % N };
            if (res.count(p)) meets = true;
          }
          if (!meets) gone.push_back(*it);
        }
        for (size_t i = 0; i < gone.size(); i++) { res.erase(gone[i]); changed = true; }
      }
      keep.clear();
      for (std::set<Key>::iterator it = res.begin(); it != res.end(); ++it)
        keep.push_back((Cube){ (*it)[0], (*it)[1], (*it)[2], dep });
    }

    if (verbose)
      printf("   depth %2d: %7ld live -> %7ld covered -> %7ld after elimination\n",
             dep, (long)live.size(), nafter_cover, (long)keep.size());
    deepest = dep;
    if (keep.empty()) { live.clear(); break; }
    if (dep == maxdepth) { residual = keep; break; }
    live.clear();
    for (size_t i = 0; i < keep.size(); i++)
      for (int t = 0; t < 8; t++)
        live.push_back((Cube){ 2*keep[i].x + (t&1), 2*keep[i].y + ((t>>1)&1),
                               2*keep[i].z + ((t>>2)&1), dep + 1 });
    dep++;
  }

  /* Theorem 1 of Cavallar-Lemmermeyer: if the eps^{-1} image of T meets
   * nothing outside T, and a single beta puts the eps image inside the covered
   * region together with T itself, then the only point of T that can be
   * exceptional is beta/(eps-1).  At most one image cell may land in T, so
   * there is such a beta exactly when precisely one image cell is still in the
   * residual set and it is congruent to T mod N. */
  if (!residual.empty() && have_eps) {
    long N = 1L << deepest;
    std::set<Key> res;
    for (size_t i = 0; i < residual.size(); i++)
      { Key q = { residual[i].x, residual[i].y, residual[i].z }; res.insert(q); }
    std::set<Key> cands;
    long unresolved = 0;
    for (size_t i = 0; i < residual.size(); i++) {
      long cv[D] = { residual[i].x, residual[i].y, residual[i].z };
      long b1[D], b2[D];
      for (int l = 0; l < D; l++) {
        long t1 = 0, t2 = 0;
        for (int j = 0; j < D; j++)
          { t1 += ME[l][j] * cv[j]; t2 += MEI[l][j] * cv[j]; }
        b1[l] = t1; b2[l] = t2;
      }
      bool cond2 = true;
      for (size_t oi = 0; oi < OFFI.size() && cond2; oi++) {
        long A = b2[0]+OFFI[oi][0], B2 = b2[1]+OFFI[oi][1], C2 = b2[2]+OFFI[oi][2];
        Key p = { (A%N+N)%N, (B2%N+N)%N, (C2%N+N)%N };
        if (res.count(p) && (((A-cv[0])%N) || ((B2-cv[1])%N) || ((C2-cv[2])%N)))
          cond2 = false;
      }
      bool found = false; long bb[D] = {0,0,0};
      if (cond2) {
        long nhit = 0, p0[D] = {0,0,0};
        for (size_t oi = 0; oi < OFF.size() && nhit < 2; oi++) {
          long A = b1[0]+OFF[oi][0], B2 = b1[1]+OFF[oi][1], C2 = b1[2]+OFF[oi][2];
          Key p = { (A%N+N)%N, (B2%N+N)%N, (C2%N+N)%N };
          if (!res.count(p)) continue;
          nhit++; p0[0]=A; p0[1]=B2; p0[2]=C2;
        }
        if (nhit == 1 && !((p0[0]-cv[0])%N) && !((p0[1]-cv[1])%N) && !((p0[2]-cv[2])%N)) {
          found = true;
          for (int l = 0; l < D; l++) bb[l] = (p0[l]-cv[l])/N;
        }
      }
      if (found) { Key q = {bb[0],bb[1],bb[2]}; cands.insert(q); }
      else unresolved++;
    }
    printf("   Thm 1: %ld residual cubes -> %ld candidate point(s), %ld unresolved"
           "  [diagnostic only -- candidates are NOT evaluated]\n",
           (long)residual.size(), (long)cands.size(), unresolved);
    /* The reduction narrows the possible exceptional points of a residual cube
     * to beta/(eps-1).  It does NOT show that those points are covered.
     * Clearing `residual` here reached PROVED without ever evaluating
     * M(K, beta/(eps-1)) -- and a field that is not norm-Euclidean can have an
     * exceptional point fixed by the unit action, which is exactly the case
     * this branch would wave through.  Until each candidate is evaluated
     * rigorously the honest outcome is RESIDUAL, so the reduction now reports
     * and proves nothing.  (No shipped result depended on it: `residual` is
     * filled only on reaching maxdepth, and all 58 finish by depth 12.) */
  }

  double secs = (double)(clock() - t0) / CLOCKS_PER_SEC;
  printf("%ld | k=%.6f | cubes=%ld covered=%ld residual=%ld deepest=%d | %.2fs\n",
         disc, k, ncubes, ncovered, (long)residual.size(), deepest, secs);
  if (residual.empty()) {
    printf("PROVED  every point of the fundamental domain is %.4f-covered:"
           " M(K) < %.4f\n", k, k);
    return 0;
  }
  printf("RESIDUAL %ld cubes at depth %d\n", (long)residual.size(), deepest);
  return 1;
}
