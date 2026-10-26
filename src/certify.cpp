/* certify.cpp -- rigorous certification of a norm-Euclidean disproof.
 *
 * Decides, for a monic integer cubic f with one real place, a modulus beta and
 * a residue theta given as exact integer coordinates in an integral basis of
 * O_K, whether
 *
 *     m(theta + beta O_K) = min { |N(theta - beta eta)| : eta in O_K } >= |N beta|,
 *
 * By the certificate criterion (Proposition 1.3 of the paper), this certifies
 * M(K) >= 1: K is not norm-Euclidean.
 *
 * No computer algebra system is involved.  The algebraic input -- an integral
 * basis of O_K and a unit eps -- is read from a file prepared once per field by
 * src/nfdata.gp (PARI) or src/nfdata.m (Magma); src/cmp_nfdata.py checks that
 * the two span the same lattice.  Everything else is done here:
 *
 *   exact integer/rational arithmetic   FLINT   (fmpz, fmpq_mat)
 *   certified roots and geometry        Arb     (arb, acb, arb_mat)
 *   lattice reduction                   fplll
 *
 * Checked here, exactly, so that they are not assumptions:
 *   - the basis contains 1 and its multiplication table is integral, so it
 *     spans an order;
 *   - |N eps| = 1, so eps is a unit;
 *   - beta | eps^n - 1, so the coset is eps^n-stable
 *     (the unit-invariance lemma, Lemma 4.1 of the paper);
 *   - the reduction transform is unimodular.
 * Assumed (this is the whole of the rigour contract's algebraic input):
 *   - that the basis spans the *maximal* order.
 *
 * Notation relative to the paper:
 *   theta (stored as th) denotes alpha; eta is an algebraic integer (the paper uses gamma);
 *   s1 = a_i, s2 = b_i; W is the reduced basis A_i in real-first order;
 *   TG (called tgt in comments) is the coset coordinate centre c_i;
 *   nsl is the implementation's slice count, with a certified upper bound.
 * Embedding vectors and matrix rows use (z,x,y) = (real,Re complex,Im complex),
 * not the paper's (x,y,z).  If P maps (x,y,z) to (z,x,y), the code's scaling
 * is P D_i P^{-1} and its W corresponds to P A_i; coefficient vectors such
 * as c_i do not change under this row permutation.  The exhaustive-box
 * lemma cited below is Lemma 4.2 of the paper.  Integral-basis
 * coordinates in the input/output records are unaffected.
 *
 * eps is not assumed to be fundamental and the class number is never used, so
 * no step of this verification is GRH-conditional.
 *
 *   certify < records
 *
 * Each input line is self-contained: the five field columns of
 * results/certificates.txt followed by the certificate itself, so one process
 * verifies any number of certificates for any number of fields.
 * where each certificate line is
 *   n <n> <t0> <t1> <t2>              (beta = eps^n - 1)
 *   b <b0> <b1> <b2> <t0> <t1> <t2>   (beta given explicitly)
 */
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/fmpq.h>
#include <flint/fmpq_mat.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_vec.h>
#include <flint/fmpz_poly.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/arb_mat.h>
#include <flint/arb_fmpz_poly.h>
#include <fplll.h>

#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <string>
#include <vector>

using namespace fplll;

static const int D = 3;
static slong PREC = 512;

/* ---------------------------------------------------------------- */
/* the order: basis, multiplication table, exact norms               */
/* ---------------------------------------------------------------- */
struct Order {
  fmpz_poly_t f;                 /* the defining polynomial, monic */
  fmpz_mat_t B;                  /* omega_j = (sum_k B[j][k] x^k) / den */
  fmpz_t den;
  fmpz *tab;                     /* structure constants c[i][j][l] */
  fmpz_t eps[D];                 /* a unit, in basis coordinates */

  fmpz *c(int i, int j, int l) { return tab + ((i * D + j) * D + l); }
};

/* omega_i * omega_j expressed in the basis; returns 0 if not integral */
static int build_table(Order &O)
{
  int ok = 1;
  fmpq_mat_t Bq, Binv;
  fmpq_mat_init(Bq, D, D); fmpq_mat_init(Binv, D, D);
  /* Bq[k][j] = coefficient of x^k in omega_j  (columns are basis elements) */
  for (int j = 0; j < D; j++)
    for (int k = 0; k < D; k++) {
      fmpq_set_fmpz_frac(fmpq_mat_entry(Bq, k, j),
                         fmpz_mat_entry(O.B, j, k), O.den);
    }
  if (!fmpq_mat_inv(Binv, Bq)) { fmpq_mat_clear(Bq); fmpq_mat_clear(Binv); return 0; }

  for (int i = 0; i < D && ok; i++)
    for (int j = 0; j < D && ok; j++) {
      /* product of the two basis elements as a polynomial, reduced mod f */
      fmpq_poly_t pi, pj, pr, fq;
      fmpq_poly_init(pi); fmpq_poly_init(pj);
      fmpq_poly_init(pr); fmpq_poly_init(fq);
      fmpq_poly_set_fmpz_poly(fq, O.f);
      for (int k = 0; k < D; k++) {
        fmpq_t t; fmpq_init(t);
        fmpq_set_fmpz_frac(t, fmpz_mat_entry(O.B, i, k), O.den);
        fmpq_poly_set_coeff_fmpq(pi, k, t);
        fmpq_set_fmpz_frac(t, fmpz_mat_entry(O.B, j, k), O.den);
        fmpq_poly_set_coeff_fmpq(pj, k, t);
        fmpq_clear(t);
      }
      fmpq_poly_mul(pr, pi, pj);
      fmpq_poly_rem(pr, pr, fq);
      /* express in the basis: coords = Binv * (coefficients of pr) */
      fmpq_mat_t v, r;
      fmpq_mat_init(v, D, 1); fmpq_mat_init(r, D, 1);
      for (int k = 0; k < D; k++) {
        fmpq_t t; fmpq_init(t);
        fmpq_poly_get_coeff_fmpq(t, pr, k);
        fmpq_set(fmpq_mat_entry(v, k, 0), t);
        fmpq_clear(t);
      }
      fmpq_mat_mul(r, Binv, v);
      for (int l = 0; l < D; l++) {
        if (!fmpz_is_one(fmpq_denref(fmpq_mat_entry(r, l, 0)))) { ok = 0; break; }
        fmpz_set(O.c(i, j, l), fmpq_numref(fmpq_mat_entry(r, l, 0)));
      }
      fmpq_mat_clear(v); fmpq_mat_clear(r);
      fmpq_poly_clear(pi); fmpq_poly_clear(pj);
      fmpq_poly_clear(pr); fmpq_poly_clear(fq);
    }
  fmpq_mat_clear(Bq); fmpq_mat_clear(Binv);
  return ok;
}

/* multiplication matrix of x: (M_x)[l][j] = sum_i x_i c[i][j][l] */
static void mulmat(fmpz_mat_t M, Order &O, fmpz *x)
{
  for (int l = 0; l < D; l++)
    for (int j = 0; j < D; j++) {
      fmpz_zero(fmpz_mat_entry(M, l, j));
      for (int i = 0; i < D; i++)
        fmpz_addmul(fmpz_mat_entry(M, l, j), x + i, O.c(i, j, l));
    }
}

static void nfmul(fmpz *z, Order &O, fmpz *x, fmpz *y)
{
  fmpz_t acc[D];
  for (int l = 0; l < D; l++) { fmpz_init(acc[l]); fmpz_zero(acc[l]); }
  for (int i = 0; i < D; i++)
    for (int j = 0; j < D; j++) {
      fmpz_t t; fmpz_init(t);
      fmpz_mul(t, x + i, y + j);
      for (int l = 0; l < D; l++) fmpz_addmul(acc[l], t, O.c(i, j, l));
      fmpz_clear(t);
    }
  for (int l = 0; l < D; l++) { fmpz_set(z + l, acc[l]); fmpz_clear(acc[l]); }
}

static void nfnorm(fmpz_t nrm, Order &O, fmpz *x)
{
  fmpz_mat_t M; fmpz_mat_init(M, D, D);
  mulmat(M, O, x);
  fmpz_mat_det(nrm, M);
  fmpz_abs(nrm, nrm);
  fmpz_mat_clear(M);
}

/* is x divisible by y in O_K?  solve M_y z = x over Q and test integrality */
static int nfdivides(Order &O, fmpz *y, fmpz *x)
{
  fmpz_mat_t My; fmpz_mat_init(My, D, D);
  mulmat(My, O, y);
  fmpq_mat_t A, b, z;
  fmpq_mat_init(A, D, D); fmpq_mat_init(b, D, 1); fmpq_mat_init(z, D, 1);
  fmpq_mat_set_fmpz_mat(A, My);
  for (int l = 0; l < D; l++) fmpq_set_fmpz(fmpq_mat_entry(b, l, 0), x + l);
  int ok = fmpq_mat_solve_fraction_free(z, A, b);
  if (ok)
    for (int l = 0; l < D; l++)
      if (!fmpz_is_one(fmpq_denref(fmpq_mat_entry(z, l, 0)))) { ok = 0; break; }
  fmpq_mat_clear(A); fmpq_mat_clear(b); fmpq_mat_clear(z);
  fmpz_mat_clear(My);
  return ok;
}

/* ---------------------------------------------------------------- */
/* embeddings: certified enclosures                                  */
/* ---------------------------------------------------------------- */
static void embed_real(arb_t out, Order &O, fmpz *x, const arb_t r)
{
  arb_t acc, p, c, dd;
  arb_init(acc); arb_init(p); arb_init(c); arb_init(dd);
  arb_zero(acc);
  arb_set_fmpz(dd, O.den);
  for (int j = 0; j < D; j++) {
    arb_zero(p);
    for (int k = D - 1; k >= 0; k--) {           /* Horner */
      arb_set_fmpz(c, fmpz_mat_entry(O.B, j, k));
      arb_mul(p, p, r, PREC); arb_add(p, p, c, PREC);
    }
    arb_div(p, p, dd, PREC);
    arb_set_fmpz(c, x + j);
    arb_mul(p, p, c, PREC);
    arb_add(acc, acc, p, PREC);
  }
  arb_set(out, acc);
  arb_clear(acc); arb_clear(p); arb_clear(c); arb_clear(dd);
}

static void embed_cx(acb_t out, Order &O, fmpz *x, const acb_t w)
{
  acb_t acc, p, c;
  arb_t dd;
  acb_init(acc); acb_init(p); acb_init(c); arb_init(dd);
  acb_zero(acc);
  arb_set_fmpz(dd, O.den);
  for (int j = 0; j < D; j++) {
    acb_zero(p);
    for (int k = D - 1; k >= 0; k--) {
      acb_set_fmpz(c, fmpz_mat_entry(O.B, j, k));
      acb_mul(p, p, w, PREC); acb_add(p, p, c, PREC);
    }
    acb_div_arb(p, p, dd, PREC);
    acb_set_fmpz(c, x + j);
    acb_mul(p, p, c, PREC);
    acb_add(acc, acc, p, PREC);
  }
  acb_set(out, acc);
  acb_clear(acc); acb_clear(p); acb_clear(c); arb_clear(dd);
}

/* ---------------------------------------------------------------- */
/* fplll: a unimodular transform for the scaled lattice              */
/* ---------------------------------------------------------------- */
/* V is a real 3x3 whose COLUMNS are the basis vectors.  Returns T (columns
 * operate on columns: W = V*T) as an integer matrix, verified unimodular. */
static int reduce_basis(fmpz_mat_t T, arb_mat_t V)
{
  /* scale so that the smallest entry is comfortably above 1, then round:
   * fplll wants integers, and a uniform scaling changes neither the lattice's
   * shape nor the unimodularity of the transform */
  slong emin = WORD_MAX, emax = WORD_MIN;
  {
    mag_t mx, mn, t;
    mag_init(mx); mag_init(mn); mag_init(t);
    mag_zero(mx); mag_inf(mn);
    for (int i = 0; i < D; i++)
      for (int j = 0; j < D; j++) {
        arb_get_mag(t, arb_mat_entry(V, i, j));
        if (mag_cmp(t, mx) > 0) mag_set(mx, t);
        if (!mag_is_zero(t) && mag_cmp(t, mn) < 0) mag_set(mn, t);
      }
    emin = mag_is_inf(mn) ? 0 : (slong)MAG_EXP(mn);
    emax = (slong)MAG_EXP(mx);
    mag_clear(mx); mag_clear(mn); mag_clear(t);
  }
  slong shift = 80 - emin;                    /* smallest entry -> ~2^80 */
  if (emax + shift > 100000) return 0;        /* absurdly skew; refuse */

  ZZ_mat<mpz_t> b(D, D), u(D, D);
  for (int j = 0; j < D; j++)                 /* rows of b = basis vectors */
    for (int i = 0; i < D; i++) {
      arb_t s; arb_init(s);
      arb_mul_2exp_si(s, arb_mat_entry(V, i, j), shift);
      fmpz_t z; fmpz_init(z);
      arf_get_fmpz(z, arb_midref(s), ARF_RND_NEAR);
      mpz_t m; mpz_init(m);
      fmpz_get_mpz(m, z);
      b[j][i] = m;
      mpz_clear(m); fmpz_clear(z); arb_clear(s);
    }
  u.gen_identity(D);
  int st = lll_reduction(b, u, LLL_DEF_DELTA, LLL_DEF_ETA, LM_WRAPPER);
  if (st != RED_SUCCESS) return 0;
  /* new_j = sum_k u[j][k] old_k, so in column form T = u^T */
  for (int j = 0; j < D; j++)
    for (int k = 0; k < D; k++) {
      mpz_t m; mpz_init(m);
      u[j][k].get_mpz(m);
      fmpz_set_mpz(fmpz_mat_entry(T, k, j), m);
      mpz_clear(m);
    }
  fmpz_t d; fmpz_init(d);
  fmpz_mat_det(d, T);
  int uni = fmpz_is_pm1(d);
  fmpz_clear(d);
  return uni;
}

/* ---------------------------------------------------------------- */
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

/* Set up the order, the certified roots and log|eps| from one record.
 *
 * This used to run once per process, against a field named by --nf.  But every
 * certificate is for a different field, so that pinned the design to one
 * process -- and one temporary file -- per certificate, with the field and the
 * certificate paired only by convention.  Nothing enforced the pairing, and a
 * mismatch is not detectable downstream: giving certify the wrong field for a
 * certificate yields a confident REFUTED, not an error.  The record travels on
 * the certificate's own line instead, and is read here, once per line.
 *
 * Returns 0 if the record cannot be used, having said why on stderr; the
 * caller reports REFUSED for that line and continues with the next.
 */
static int setup_field(Order &O, arb_t r, acb_t w, fmpz *ep, double *LOGEPS,
                       const std::vector<std::string> &fld)
{
  /* ---- the field, from the certificate's own line --------------- */
  {
    auto ints = [](const std::string &t) { return split(t, ','); };
    auto co = ints(fld[1]);
    auto bb = ints(fld[3]);
    auto ee = ints(fld[4]);
    /* Check the shape before indexing: split() returns whatever the line
     * happened to hold, and operator[] past the end is undefined behaviour.
     * Every string->integer conversion is checked too -- fmpz_set_str reports
     * failure rather than throwing, and an unchecked failure would leave a
     * garbage entry that the rest of this function would treat as data. */
    if (co.size() < 4 || bb.size() < (size_t)(D * D) || ee.size() < (size_t)D) {
      fprintf(stderr, "field record has too few entries\n"); return 0;
    }
    for (int i = 0; i < 4; i++) {
      fmpz_t t; fmpz_init(t);
      if (fmpz_set_str(t, co[i].c_str(), 10) != 0)
        { fprintf(stderr, "bad integer in f\n"); fmpz_clear(t); return 0; }
      fmpz_poly_set_coeff_fmpz(O.f, i, t); fmpz_clear(t);
    }
    if (fmpz_set_str(O.den, fld[2].c_str(), 10) != 0 || fmpz_is_zero(O.den))
      { fprintf(stderr, "bad or zero denominator\n"); return 0; }
    for (int j = 0; j < D; j++)
      for (int k = 0; k < D; k++)
        if (fmpz_set_str(fmpz_mat_entry(O.B, j, k), bb[3 * j + k].c_str(), 10) != 0)
          { fprintf(stderr, "bad integer in basis\n"); return 0; }
    for (int i = 0; i < D; i++)
      if (fmpz_set_str(O.eps[i], ee[i].c_str(), 10) != 0)
        { fprintf(stderr, "bad integer in eps\n"); return 0; }
  }
  { fmpz_t lc; fmpz_init(lc); fmpz_poly_get_coeff_fmpz(lc, O.f, D);
    if (!fmpz_is_one(lc)) { fprintf(stderr, "f is not monic\n"); return 0; }
    fmpz_clear(lc); }
  if (!build_table(O)) {
    fprintf(stderr, "the given basis is not closed under multiplication: "
                    "not an order\n"); return 0;
  }
  { /* omega_1 = 1 exactly.  Checking |N omega_1| = 1 is not enough: that says
     * only that omega_1 is a unit, and eps itself would pass.  Coordinate 0 is
     * used throughout as the rational integer 1 -- the power accumulator starts
     * at (1,0,0), eps^n - 1 is formed by subtracting one from it, and the
     * divisibility and stability tests use the same convention -- so a first
     * basis element u != 1 would silently substitute u*eps^n for eps^n and
     * check the wrong condition for the slicing argument.  Test the identity
     * itself: multiplication by omega_1 must be the identity map, i.e. the
     * structure constants c[0][j][l] must be delta[j][l]. */
    for (int j = 0; j < D; j++)
      for (int l = 0; l < D; l++) {
        int want = (j == l);
        if ((want && !fmpz_is_one(O.c(0, j, l))) ||
            (!want && !fmpz_is_zero(O.c(0, j, l)))) {
          fprintf(stderr, "the first basis element is not 1\n");
          return 0;
        }
      }
  }
  { fmpz *ep = _fmpz_vec_init(D);
    for (int i = 0; i < D; i++) fmpz_set(ep + i, O.eps[i]);
    fmpz_t nn; fmpz_init(nn); nfnorm(nn, O, ep);
    if (!fmpz_is_one(nn)) { fprintf(stderr, "eps is not a unit\n"); return 0; }
    fmpz_clear(nn); _fmpz_vec_clear(ep, D); }

  /* ---- certified roots of f ------------------------------------ */
  arb_zero(r); acb_zero(w);
  {
    acb_ptr rts = _acb_vec_init(D);
    slong nr = 0;
    arb_fmpz_poly_complex_roots(rts, O.f, 0, PREC);
    for (int i = 0; i < D; i++)
      if (arb_contains_zero(acb_imagref(rts + i)))
        { arb_set(r, acb_realref(rts + i)); nr++; }
    if (nr != 1) { fprintf(stderr, "cannot isolate the real root\n"); return 0; }
    for (int i = 0; i < D; i++)
      if (arb_is_positive(acb_imagref(rts + i))) acb_set(w, rts + i);
    if (acb_is_zero(w)) { fprintf(stderr, "no complex root\n"); return 0; }
    _acb_vec_clear(rts, D);
  }
  /* eps must expand the real place; replace it by its inverse if not */
  for (int i = 0; i < D; i++) fmpz_set(ep + i, O.eps[i]);
  {
    arb_t e1; arb_init(e1);
    embed_real(e1, O, ep, r);
    arb_abs(e1, e1);
    arb_t one; arb_init(one); arb_one(one);
    if (arb_lt(e1, one)) {
      /* eps^{-1}: solve M_eps y = 1 exactly */
      fmpz_mat_t M; fmpz_mat_init(M, D, D); mulmat(M, O, ep);
      fmpq_mat_t A, b, z;
      fmpq_mat_init(A, D, D); fmpq_mat_init(b, D, 1); fmpq_mat_init(z, D, 1);
      fmpq_mat_set_fmpz_mat(A, M);
      fmpq_one(fmpq_mat_entry(b, 0, 0));
      if (!fmpq_mat_solve_fraction_free(z, A, b))
        { fprintf(stderr, "cannot invert eps\n"); return 0; }
      for (int l = 0; l < D; l++) {
        if (!fmpz_is_one(fmpq_denref(fmpq_mat_entry(z, l, 0))))
          { fprintf(stderr, "eps^{-1} is not integral\n"); return 0; }
        fmpz_set(ep + l, fmpq_numref(fmpq_mat_entry(z, l, 0)));
      }
      fmpq_mat_clear(A); fmpq_mat_clear(b); fmpq_mat_clear(z); fmpz_mat_clear(M);
    } else if (!arb_gt(e1, one)) {
      fprintf(stderr, "cannot separate |eps^(1)| from 1\n"); return 0;
    }
    arb_clear(e1); arb_clear(one);
  }

  /* The working precision has to carry the dynamic range e^{n ell} of
   * the slice decomposition.  Estimate ell = log|eps_R| from the supplied
   * unit; ell equals R_K only when eps is fundamental.  The certificate's n
   * then fixes the precision.  This is only a precision heuristic, not a
   * trusted regulator computation.  Insufficient precision can prevent
   * certification of the inverse and cause slices to be refused. */
  { arb_t e1; arb_init(e1);
    embed_real(e1, O, ep, r);
    arb_abs(e1, e1);
    arb_log(e1, e1, PREC);
    *LOGEPS = arf_get_d(arb_midref(e1), ARF_RND_UP);
    if (!(*LOGEPS > 0)) { fprintf(stderr, "log|eps| not positive\n"); return 0; }
    arb_clear(e1);
  }
  return 1;
}

int main(int argc, char **argv)
{
  /* --argmin also prints the coset element attaining the minimum.  The
   * enumeration has it in hand anyway: the search is exhaustive over a
   * fundamental domain of the unit action and |N| is constant on orbits, so
   * the element attaining mymin is a minimal-norm representative of
   * theta + beta O_K, not merely the smallest one happened upon.  Printing it
   * lets the certificates be rewritten with a shorter, near-canonical theta;
   * see src/canonicalise.py. */
  int argmin = 0;
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--argmin")) argmin = 1;
    else { fprintf(stderr, "usage: certify [--argmin] < records\n"); return 1; }
  }

  Order O;
  fmpz_poly_init(O.f);
  fmpz_mat_init(O.B, D, D);
  fmpz_init(O.den);
  O.tab = _fmpz_vec_init(D * D * D);
  for (int i = 0; i < D; i++) fmpz_init(O.eps[i]);


  /* Allocated once, re-filled for each input line. */
  arb_t r; acb_t w; arb_init(r); acb_init(w);
  fmpz *ep = _fmpz_vec_init(D);
  double LOGEPS;
  arb_t sq2; arb_init(sq2);
  const slong PREC0 = PREC;

  /* ---- one field and one certificate per input line -------------- */
  static char line[1 << 16];
  long nbad = 0;
  while (fgets(line, sizeof(line), stdin)) {
    std::vector<std::string> fld = split(std::string(line), '|');
    if (fld.size() < 6) {
      printf("REFUSED malformed input line\n"); fflush(stdout); nbad++; continue;
    }
    PREC = PREC0;             /* the previous line raised it; start clean */
    if (!setup_field(O, r, w, ep, &LOGEPS, fld)) {
      printf("REFUSED the field record could not be verified\n");
      fflush(stdout); nbad++; continue;
    }
    printf("# order verified: omega_1 = 1, table integral, |N eps| = 1; "
           "log|eps| = %.3f\n", LOGEPS);
    const char *cert = fld[5].c_str();
    fmpz *bta = _fmpz_vec_init(D), *th = _fmpz_vec_init(D);
    fmpz *epn = _fmpz_vec_init(D);
    long n = 1;
    char kind = 0;
    { char b0[8192], b1[8192], b2[8192], t0[8192], t1[8192], t2[8192];
      /* A parse failure here used to fall through a bare continue: no verdict
       * printed, nbad untouched, allocations leaked.  The driver noticed only
       * because it compares verdict and input counts; run standalone, certify
       * would omit the record and still exit 0.  Every nonempty line must now
       * produce exactly one verdict, and any malformed line makes the process
       * fail.  Integer conversions are checked, not assumed. */
      int ok = 1;
      if (sscanf(cert, " %c", &kind) != 1) ok = 0;
      else if (kind == 'n') {
        if (sscanf(cert, " n %ld %8191s %8191s %8191s", &n, t0, t1, t2) != 4) ok = 0;
        else if (n < 0 || n > 4096) ok = 0;
      } else if (kind == 'b') {
        if (sscanf(cert, " b %8191s %8191s %8191s %8191s %8191s %8191s",
                   b0, b1, b2, t0, t1, t2) != 6) ok = 0;
        else if (fmpz_set_str(bta + 0, b0, 10) != 0 ||
                 fmpz_set_str(bta + 1, b1, 10) != 0 ||
                 fmpz_set_str(bta + 2, b2, 10) != 0) ok = 0;
      } else ok = 0;
      if (ok && (fmpz_set_str(th + 0, t0, 10) != 0 ||
                 fmpz_set_str(th + 1, t1, 10) != 0 ||
                 fmpz_set_str(th + 2, t2, 10) != 0)) ok = 0;
      if (!ok) {
        printf("REFUSED malformed certificate record\n"); fflush(stdout); nbad++;
        _fmpz_vec_clear(bta, D); _fmpz_vec_clear(th, D); _fmpz_vec_clear(epn, D);
        continue;
      }
    }
    /* eps^n, and beta = eps^n - 1 when not given */
    { fmpz *acc = _fmpz_vec_init(D);
      fmpz_one(acc + 0);
      for (long i = 0; i < n; i++) nfmul(acc, O, acc, ep);
      for (int l = 0; l < D; l++) fmpz_set(epn + l, acc + l);
      _fmpz_vec_clear(acc, D);
    }
    if (kind == 'n') {
      for (int l = 0; l < D; l++) fmpz_set(bta + l, epn + l);
      fmpz_sub_ui(bta + 0, bta + 0, 1);
    } else {
      /* find the least t with beta | eps^t - 1 */
      long t; fmpz *acc = _fmpz_vec_init(D), *dd = _fmpz_vec_init(D);
      fmpz_one(acc + 0);
      for (t = 1; t <= 64; t++) {
        nfmul(acc, O, acc, ep);
        for (int l = 0; l < D; l++) fmpz_set(dd + l, acc + l);
        fmpz_sub_ui(dd + 0, dd + 0, 1);
        if (nfdivides(O, bta, dd)) break;
      }
      if (t > 64) {
        printf("REFUSED no t <= 64 with beta | eps^t - 1\n"); fflush(stdout); nbad++;
        _fmpz_vec_clear(acc, D); _fmpz_vec_clear(dd, D);
        _fmpz_vec_clear(bta, D); _fmpz_vec_clear(th, D); _fmpz_vec_clear(epn, D);
        continue;
      }
      n = t;
      for (int l = 0; l < D; l++) fmpz_set(epn + l, acc + l);
      _fmpz_vec_clear(acc, D); _fmpz_vec_clear(dd, D);
    }
    /* exact: the coset must be eps^n-stable */
    { fmpz *dd = _fmpz_vec_init(D);
      for (int l = 0; l < D; l++) fmpz_set(dd + l, epn + l);
      fmpz_sub_ui(dd + 0, dd + 0, 1);
      if (!nfdivides(O, bta, dd)) {
        printf("REFUSED beta does not divide eps^n - 1\n"); fflush(stdout); nbad++;
        _fmpz_vec_clear(dd, D);
        _fmpz_vec_clear(bta, D); _fmpz_vec_clear(th, D); _fmpz_vec_clear(epn, D);
        continue;
      }
      _fmpz_vec_clear(dd, D);
    }
    /* exact: alpha (stored as th) must not lie in beta O_K;
     * see the certificate criterion, Proposition 1.3 of the paper. */
    if (nfdivides(O, bta, th)) {
      printf("REFUSED theta in beta O_K\n"); fflush(stdout); nbad++;
        _fmpz_vec_clear(bta, D); _fmpz_vec_clear(th, D); _fmpz_vec_clear(epn, D);
      continue;
    }

    /* set the precision from this certificate's n, and redo the enclosures */
    PREC = (slong)(4.0 * (n + 2) * LOGEPS / M_LN2) + 1024;
    { acb_ptr rts = _acb_vec_init(D);
      arb_fmpz_poly_complex_roots(rts, O.f, 0, PREC);
      for (int i = 0; i < D; i++)
        if (arb_contains_zero(acb_imagref(rts + i))) arb_set(r, acb_realref(rts + i));
      for (int i = 0; i < D; i++)
        if (arb_is_positive(acb_imagref(rts + i))) acb_set(w, rts + i);
      _acb_vec_clear(rts, D);
    }
    arb_set_ui(sq2, 2); arb_sqrt(sq2, sq2, PREC);

    fmpz_t NB; fmpz_init(NB); nfnorm(NB, O, bta);

    /* embeddings of the lattice basis beta*omega_j, and of theta */
    arb_t Br1[D]; acb_t Bc[D]; arb_t ct1; acb_t ct2;
    for (int j = 0; j < D; j++) { arb_init(Br1[j]); acb_init(Bc[j]); }
    arb_init(ct1); acb_init(ct2);
    for (int j = 0; j < D; j++) {
      fmpz *e = _fmpz_vec_init(D), *bw = _fmpz_vec_init(D);
      fmpz_one(e + j);
      nfmul(bw, O, bta, e);
      embed_real(Br1[j], O, bw, r);
      embed_cx(Bc[j], O, bw, w);
      _fmpz_vec_clear(e, D); _fmpz_vec_clear(bw, D);
    }
    embed_real(ct1, O, th, r);
    embed_cx(ct2, O, th, w);

    /* number of slices: 2^nsl must cover |eps^n| at the real place */
    long nsl;
    { arb_t E; arb_init(E); embed_real(E, O, epn, r); arb_abs(E, E);
      mag_t hi; mag_init(hi); arb_get_mag(hi, E);
      nsl = (long)MAG_EXP(hi) + 1;
      if (nsl < 1) nsl = 1;
      mag_clear(hi); arb_clear(E);
    }

    fmpz_t mymin; fmpz_init(mymin); fmpz_zero(mymin);
    fmpz *mymin_at = _fmpz_vec_init(D);      /* the element attaining it */
    int refused = 0;
    /* The cap the box is built for.  The exhaustive-box lemma covers
     * norms < CAP and no further, so the least norm found is m only when it
     * lands at or below CAP.  The verdict needs CAP = |N beta| and nothing more
     * -- excluding everything below |N beta| is the whole certificate -- but
     * --argmin claims a *minimal* representative, and that claim needs the value
     * to be inside the proved range.  So under --argmin the sweep is repeated
     * with the cap raised to the value found.
     *
     * This terminates after one extra sweep, but the cap must be set one above
     * the value found, not equal to it (Codex R2).  The exhaustive-box lemma
     * covers norms *strictly* below CAP, so CAP := mymin would leave the point just
     * found on the boundary, where the lemma says nothing -- and that point was
     * only ever known to lie in the enumerated rectangle, which is rebuilt from
     * a different LLL basis at each cap, so membership is not inherited.  At
     * CAP := mymin + 1 the point is strictly inside the region the lemma covers,
     * so the next sweep provably enumerates it and the value is settled.
     * Nothing can appear below |N beta| that the first sweep missed, either, so
     * the verdict never changes -- only what can be said about the value.
     *
     * (min_proved also holds when mymin == CAP exactly: no point below CAP
     * exists, and a point of norm CAP is in hand, so m = CAP.) */
    /* |N theta| is a cap in its own right: theta lies in the coset, so
     * m <= |N theta|, and at any cap strictly above |N theta| the exhaustive-box
     * lemma puts theta itself inside the box -- so the value found is at or below the cap
     * and m is proved in a single sweep.  It is not used as the first cap
     * because the cost is linear in it and theta is size-reduced, not
     * norm-reduced: |N theta|/|N beta| has median 13 and mean 153 over the
     * shipped file.  It serves instead as the fallback for the second sweep,
     * which makes the escalation terminate in two sweeps unconditionally and
     * bounds the second sweep's cost. */
    fmpz_t Nth; fmpz_init(Nth); nfnorm(Nth, O, th); fmpz_abs(Nth, Nth);
    fmpz_add_ui(Nth, Nth, 1);
    fmpz_t CAP; fmpz_init(CAP); fmpz_set(CAP, NB);
    int min_proved = 0;
    for (int round = 0; round < 16 && !refused; round++) {
    fmpz_zero(mymin);
    for (long i = 0; i < nsl && !refused; i++) {
      arb_mat_t V, W, Wi, C0, TG;
      arb_t s1, s2, t1;
      arb_mat_init(V, D, D); arb_mat_init(W, D, D); arb_mat_init(Wi, D, D);
      arb_mat_init(C0, D, 1); arb_mat_init(TG, D, 1);
      arb_init(s1); arb_init(s2); arb_init(t1);
      arb_one(s1); arb_mul_2exp_si(s1, s1, i + 1);
      arb_set_fmpz(s2, CAP); arb_mul_2exp_si(s2, s2, -i); arb_sqrt(s2, s2, PREC);
      for (int j = 0; j < D; j++) {
        arb_div(arb_mat_entry(V, 0, j), Br1[j], s1, PREC);
        arb_div(arb_mat_entry(V, 1, j), acb_realref(Bc[j]), s2, PREC);
        arb_div(arb_mat_entry(V, 2, j), acb_imagref(Bc[j]), s2, PREC);
      }
      arb_div(arb_mat_entry(C0, 0, 0), ct1, s1, PREC);
      arb_div(arb_mat_entry(C0, 1, 0), acb_realref(ct2), s2, PREC);
      arb_div(arb_mat_entry(C0, 2, 0), acb_imagref(ct2), s2, PREC);

      fmpz_mat_t T; fmpz_mat_init(T, D, D);
      if (!reduce_basis(T, V)) { refused = 1; goto slice_done; }
      for (int a1 = 0; a1 < D; a1++)
        for (int b1 = 0; b1 < D; b1++) {
          arb_zero(arb_mat_entry(W, a1, b1));
          for (int k = 0; k < D; k++) {
            arb_t tt; arb_init(tt);
            arb_set_fmpz(tt, fmpz_mat_entry(T, k, b1));
            arb_mul(tt, tt, arb_mat_entry(V, a1, k), PREC);
            arb_add(arb_mat_entry(W, a1, b1), arb_mat_entry(W, a1, b1), tt, PREC);
            arb_clear(tt);
          }
        }
      if (!arb_mat_inv(Wi, W, PREC)) { refused = 1; goto slice_done; }
      arb_mat_mul(TG, Wi, C0, PREC);

      { long h[D], hup[D]; fmpz *ctr = _fmpz_vec_init(D);
        for (int j = 0; j < D; j++) {
          arb_t nr2, tt; arb_init(nr2); arb_init(tt); arb_zero(nr2);
          for (int k = 0; k < D; k++) {
            arb_sqr(tt, arb_mat_entry(Wi, j, k), PREC);
            arb_add(nr2, nr2, tt, PREC);
          }
          arb_sqrt(nr2, nr2, PREC);
          arb_mul(nr2, nr2, sq2, PREC);          /* the exhaustive-box lemma: sqrt2*||W^-1_j|| */
          arb_set_d(tt, 0.5); arb_add(nr2, nr2, tt, PREC);
          arf_t hi2; arf_init(hi2); arb_get_ubound_arf(hi2, nr2, PREC);
          /* R1: keep this exact.  Routing the certified bound through a double
           * and casting to long is undefined when the bound is infinite or out
           * of range, so the intended h > 80 refusal was not guaranteed to run
           * and a bad conversion could land on a small or negative h -- an
           * undersized box, in the direction that manufactures certificates.
           * Take the floor into an fmpz, range-check exactly, and narrow only
           * after the check has passed. */
          h[j] = 0;
          if (!arf_is_finite(hi2)) refused = 1;
          else {
            fmpz_t hz; fmpz_init(hz);
            arf_get_fmpz(hz, hi2, ARF_RND_FLOOR);
            if (fmpz_sgn(hz) < 0) fmpz_zero(hz);
            if (fmpz_cmp_ui(hz, 80) > 0) refused = 1;
            else h[j] = fmpz_get_si(hz);
            fmpz_clear(hz);
          }
          arf_clear(hi2); arb_clear(nr2); arb_clear(tt);
          /* Centre c_i (TG): round to a nearest integer, widening if needed.
           * The paper writes floor(c_i + 1/2).  Arb's tie convention can
           * differ, but either nearest integer has error at most 1/2, which
           * is all the exhaustive-box bound requires. */
          arf_t lo3, hi3; arf_init(lo3); arf_init(hi3);
          fmpz_t zl, zh; fmpz_init(zl); fmpz_init(zh);
          arb_get_lbound_arf(lo3, arb_mat_entry(TG, j, 0), PREC);
          arb_get_ubound_arf(hi3, arb_mat_entry(TG, j, 0), PREC);
          arf_get_fmpz(zl, lo3, ARF_RND_NEAR);
          arf_get_fmpz(zh, hi3, ARF_RND_NEAR);
          /* The target is an enclosure, not a point.  The true centre may be
           * anywhere in [lo3, hi3], so the integers that can satisfy the exhaustive-box lemma
           * are the union of the boxes about every possible centre, which is
           * [zl - h, zh + h].  Widening by one about zl is correct only when
           * zh - zl <= 1; a wider enclosure would drop coset points, and
           * dropping points inflates m -- the one direction that manufactures
           * a false certificate.  A span too large to enumerate is a precision
           * failure, so refuse the slice rather than guess. */
          { fmpz_t sp; fmpz_init(sp); fmpz_sub(sp, zh, zl);
            long span = 0;
            if (fmpz_sgn(sp) < 0 || !fmpz_fits_si(sp)) refused = 1;
            else { span = fmpz_get_si(sp); if (span > 80) { refused = 1; span = 0; } }
            hup[j] = h[j] + span;
            fmpz_clear(sp); }
          fmpz_set(ctr + j, zl);
          arf_clear(lo3); arf_clear(hi3); fmpz_clear(zl); fmpz_clear(zh);
        }
        if (!refused)
        for (long e0 = -h[0]; e0 <= hup[0]; e0++)
        for (long e1b = -h[1]; e1b <= hup[1]; e1b++)
        for (long e2 = -h[2]; e2 <= hup[2]; e2++) {
          fmpz *u = _fmpz_vec_init(D), *eta = _fmpz_vec_init(D),
               *g = _fmpz_vec_init(D);
          fmpz_add_si(u + 0, ctr + 0, e0);
          fmpz_add_si(u + 1, ctr + 1, e1b);
          fmpz_add_si(u + 2, ctr + 2, e2);
          for (int l = 0; l < D; l++) {
            fmpz_zero(eta + l);
            for (int k = 0; k < D; k++)
              fmpz_addmul(eta + l, fmpz_mat_entry(T, l, k), u + k);
          }
          nfmul(g, O, bta, eta);
          for (int l = 0; l < D; l++) fmpz_sub(g + l, th + l, g + l);
          if (!_fmpz_vec_is_zero(g, D)) {
            fmpz_t v; fmpz_init(v); nfnorm(v, O, g);
            if (!fmpz_is_zero(v) && (fmpz_is_zero(mymin) || fmpz_cmp(v, mymin) < 0))
              { fmpz_set(mymin, v); _fmpz_vec_set(mymin_at, g, D); }
            fmpz_clear(v);
          }
          _fmpz_vec_clear(u, D); _fmpz_vec_clear(eta, D); _fmpz_vec_clear(g, D);
        }
        _fmpz_vec_clear(ctr, D);
      }
    slice_done:
      fmpz_mat_clear(T);
      arb_mat_clear(V); arb_mat_clear(W); arb_mat_clear(Wi);
      arb_mat_clear(C0); arb_mat_clear(TG);
      arb_clear(s1); arb_clear(s2); arb_clear(t1);
    }
    if (refused) break;
    min_proved = (!fmpz_is_zero(mymin) && fmpz_cmp(mymin, CAP) <= 0);
    if (min_proved || !argmin) break;
    /* Raise the cap to min(mymin + 1, |N theta| + 1) -- and to |N theta| + 1
     * outright if the box came up empty, where there is no found value to aim
     * at.  Whichever is taken, a known coset point (the minimiser, or theta
     * itself) has norm strictly below the new cap, so the next sweep is
     * guaranteed to enumerate it and land at or below that cap.  At most two
     * sweeps. */
    { fmpz_t nxt; fmpz_init(nxt);
      if (fmpz_is_zero(mymin)) fmpz_set(nxt, Nth);
      else { fmpz_add_ui(nxt, mymin, 1);
             if (fmpz_cmp(Nth, nxt) < 0) fmpz_set(nxt, Nth); }
      fmpz_set(CAP, nxt); fmpz_clear(nxt); }
    }

    if (argmin && !refused && !fmpz_is_zero(mymin)) {
      /* Balance the minimiser within its unit orbit.
       *
       * |N| is constant on the orbit {eps^(nk) gamma}, which stays inside the
       * coset because beta divides eps^n - 1, so every member is an equally
       * valid minimal-norm representative -- but they differ enormously in
       * size, each step scaling gamma_R by |eps^n| and gamma_C by its inverse
       * square root.  The enumeration walks slices in increasing |gamma_R| and
       * keeps the first element attaining the minimum, so what it reports is
       * the member with the *smallest* real embedding, hence the largest
       * complex one and the largest coordinates: for Method I certificates,
       * theta of 188 digits for an m of a few dozen.
       *
       * Minimising T_2 = |gamma_R|^2 + 2|gamma_C|^2 subject to
       * |gamma_C|^2|gamma_R| = m gives |gamma_R| = |gamma_C| = m^(1/3) and
       * T_2 = 3m^(2/3).  The orbit only offers |gamma_R| in steps of
       * n*log|eps|, so take the step nearest that. */
      { fmpz *inv = _fmpz_vec_init(D), *cur = _fmpz_vec_init(D);
        fmpq_mat_t A, Ai; fmpq_mat_init(A, D, D); fmpq_mat_init(Ai, D, D);
        /* matrix of multiplication by eps^n, then its (integral) inverse */
        for (int l = 0; l < D; l++)
          for (int j2 = 0; j2 < D; j2++) {
            fmpz_t acc2; fmpz_init(acc2); fmpz_zero(acc2);
            for (int i2 = 0; i2 < D; i2++)
              fmpz_addmul(acc2, epn + i2, O.tab + ((i2*D + j2)*D + l));
            fmpq_set_fmpz(fmpq_mat_entry(A, l, j2), acc2);
            fmpz_clear(acc2);
          }
        int okinv = fmpq_mat_inv(Ai, A);
        /* omega_1 = 1 is checked above, so eps^-n is the first column */
        if (okinv)
          for (int l = 0; l < D; l++)
            if (fmpz_is_one(fmpq_denref(fmpq_mat_entry(Ai, l, 0))))
              fmpz_set(inv + l, fmpq_numref(fmpq_mat_entry(Ai, l, 0)));
            else okinv = 0;
        if (okinv) {
          arb_t s1b; arb_init(s1b);
          embed_real(s1b, O, mymin_at, r); arb_abs(s1b, s1b);
          arb_log(s1b, s1b, PREC);
          double logs1 = arf_get_d(arb_midref(s1b), ARF_RND_NEAR);
          double logm = fmpz_dlog(mymin);
          double step = (double)n * LOGEPS;
          /* Minimise T_2, not |Delta|.  With log|gamma_R| = log(m)/3 + Delta
           * the T_2 form is proportional to exp(2*Delta) + 2*exp(-Delta),
           * which is asymmetric: overshooting upward costs exp(Delta),
           * downward only exp(-Delta/2).  Rounding Delta to nearest is
           * therefore the wrong objective whenever the orbit step n*log|eps|
           * is coarse, as it is here -- some 50 nats, or 22 digits.  Compared
           * in the log domain because Delta reaches several hundred and the
           * exponentials would otherwise overflow. */
          long k0 = (step > 0) ? (long)llround((logm/3.0 - logs1)/step) : 0;
          long kk = k0; double best = 0; int have = 0;
          for (long jj = -3; jj <= 3; jj++) {
            double dl = (logs1 + (double)(k0 + jj) * step) - logm/3.0;
            double a = 2.0*dl, b = -dl + log(2.0);
            double M = (a > b) ? a : b;
            double sc = M + log(exp(a - M) + exp(b - M));
            if (!have || sc < best) { best = sc; kk = k0 + jj; have = 1; }
          }
          for (long t2 = 0; t2 < labs(kk); t2++)
            nfmul(mymin_at, O, mymin_at, kk > 0 ? epn : inv);
          /* the norm is an invariant of the orbit; check we kept it */
          { fmpz_t chk; fmpz_init(chk); nfnorm(chk, O, mymin_at);
            if (!fmpz_equal(chk, mymin))
              printf("BALANCE FAILED norm changed\n");
            fmpz_clear(chk); }
          arb_clear(s1b);
        }
        fmpq_mat_clear(A); fmpq_mat_clear(Ai);
        _fmpz_vec_clear(inv, D); _fmpz_vec_clear(cur, D);
      }
      char *a0 = fmpz_get_str(NULL, 10, mymin_at + 0),
           *a1 = fmpz_get_str(NULL, 10, mymin_at + 1),
           *a2 = fmpz_get_str(NULL, 10, mymin_at + 2);
      /* Minimality is proved only when the attaining point lies at or below
       * the cap the box was built for: the exhaustive-box lemma makes the enumeration
       * exhaustive for norms < |N beta|, so a larger value is an upper bound
       * for m, not m itself.  Label the two cases apart -- canonicalise.py
       * substitutes only what is marked MINIMAL. */
      printf("%s %s %s %s\n", min_proved ? "MINIMAL" : "CANDIDATE", a0, a1, a2);
      flint_free(a0); flint_free(a1); flint_free(a2);
    }
    if (refused) { printf("REFUSED a slice could not be certified\n"); nbad++; }
    else if (fmpz_is_zero(mymin)) {
      char *sn = fmpz_get_str(NULL, 10, NB);
      printf("CERTIFIED no coset point below |N beta| = %s\n", sn);
      flint_free(sn);
    } else if (fmpz_cmp(mymin, NB) >= 0) {
      char *sm = fmpz_get_str(NULL, 10, mymin), *sn = fmpz_get_str(NULL, 10, NB);
      /* The box is exhaustive for norms < CAP only.  A point found in it above
       * CAP bounds m from above but does not determine it -- another coset
       * point with norm between the two need not lie in the box.  With
       * --argmin the cap was raised until the value fell inside the proved
       * range, so there m really is determined. */
      if (min_proved)
        printf("CERTIFIED m = %s >= |N beta| = %s\n", sm, sn);
      else
        printf("CERTIFIED |N beta| = %s <= m <= %s\n", sn, sm);
      flint_free(sm); flint_free(sn);
    } else {
      char *sm = fmpz_get_str(NULL, 10, mymin), *sn = fmpz_get_str(NULL, 10, NB);
      printf("REFUTED m = %s < |N beta| = %s\n", sm, sn);
      flint_free(sm); flint_free(sn); nbad++;
    }
    fflush(stdout);
    for (int j = 0; j < D; j++) { arb_clear(Br1[j]); acb_clear(Bc[j]); }
    arb_clear(ct1); acb_clear(ct2); fmpz_clear(NB); fmpz_clear(mymin);
    fmpz_clear(CAP); fmpz_clear(Nth);
    _fmpz_vec_clear(mymin_at, D);
    _fmpz_vec_clear(bta, D); _fmpz_vec_clear(th, D); _fmpz_vec_clear(epn, D);
  }
  return nbad ? 4 : 0;
}
