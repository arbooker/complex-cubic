/* ne_chain.c -- exact Euclidean minimum at points with a finite unit orbit.
 *
 * For xi = theta/beta with beta | eps^n - 1, the coset S = theta + beta*O_K is
 * stable under multiplication by eps^n and |N| is constant on its orbits, so
 *
 *     M(K, xi) = min_{gamma in S} |N gamma| / |N beta|
 *
 * and the minimum may be sought in a single fundamental domain of the eps^n
 * action, |x^(1)| in [1, E) with E = |(eps^n)^(1)|.  That interval spans
 * e^{n*Reg}, so it is walked in binary slices; within a slice the region
 * {|x^(1)| ~ s, |N| < cap} is a box, and the lattice is LLL-reduced against
 * that box's metric before enumerating.
 *
 * Notation relative to the paper: theta is alpha, xi = alpha/beta, and
 * gamma is an algebraic integer.  ct denotes the embedding vector stored as
 * (z,x,y) = (real,Re complex,Im complex), rather than the paper's (x,y,z).
 * Likewise y[0],y[1],y[2] correspond to v_3,v_1,v_2; s1,s2 are a_i,b_i.
 * W is the reduced basis A_i with its rows permuted to this order, and tgt
 * is c_i.  The member named A below is W^{-1} times the real-first scaling
 * matrix, NOT the paper's A_i.  nsl remains a useful internal slice count.
 *
 * Two things make this fast enough to use on thousands of candidate theta:
 *
 *  - the slice geometry does not depend on theta, so the per-slice reduced
 *    bases are computed once per (field, n) and reused;
 *  - only candidates whose approximate norm is below the cap are re-evaluated
 *    exactly in PARI.
 *
 * The per-slice search is a complete Fincke-Pohst enumeration (qfcvp) of the
 * coset inside the slice's box, not Babai rounding with a fudge box: LLL exists
 * precisely to make the complete enumeration cheap, and a box around the Babai
 * point provably misses vectors -- which is what kept d = -680 at 0.93.
 *
 * Floating point therefore only decides *where* to look; every reported value
 * is an exact integer norm.
 *
 * Precision.  The search region spans e^{n*Reg} in the real direction and Reg
 * reaches ~400 here, so the scaled bases and their LLL transforms carry
 * entries of that size: both are kept as PARI GENs at a precision set from
 * n*Reg, never squeezed into machine words.  What *is* O(1), and so done in
 * doubles, is the reduced basis and the enumeration inside a slice -- the
 * reduced lattice has covolume O(1) by construction.  An earlier version tried
 * to hold the transform in a long and lost 17 slices of 238 at Reg = 403.
 *
 *   ne_chain --poly a0,a1,a2,a3 --n N [--cap C] < thetas.txt
 *     one "t0 t1 t2" per line
 *   ne_chain --poly ... --centres [--auto-n] < centres.txt
 *     one "m0 m1 m2 D" per line, a cube centre (m_i / 2^D) from ne_cover;
 *     theta = round((eps^n - 1) * xi) is formed here, in exact arithmetic,
 *     because eps^n - 1 has coordinates of size e^{n*Reg/2} and cannot go
 *     through a double.  The 27 nearest lattice points are tried, since the
 *     cube's image is a parallelepiped of volume N(eps^n-1)*delta^3 and holds
 *     only a handful.
 *
 * Prints "m Nbeta t0 t1 t2" for the best theta seen.
 */

#include <pari/pari.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define DEG 3
typedef unsigned __int128 u128;

/* frac(x) as a fixed-point integer in units of 2^-128.  x is O(e^{n Reg}) and
 * only its fractional part matters, so this is where the multiprecision ends
 * and machine arithmetic begins. */
static u128 frac_u128(GEN x)
{
  pari_sp av = avma;
  GEN f = gsub(x, gfloor(x));                 /* in [0,1) */
  GEN y = gfloor(gmul(f, int2n(128)));
  u128 hi = (u128)itou(shifti(y, -64));
  u128 lo = (u128)itou(remi2n(y, 64));
  set_avma(av);
  return (hi << 64) | lo;
}

/* (a * t) mod 2^128, for a a 128-bit fraction and t a trial index */
static inline u128 mul128(u128 a, unsigned long t)
{
  u128 lo = (u128)(unsigned long)a * t;
  u128 hi = (u128)(unsigned long)(a >> 64) * t;
  return lo + (hi << 64);
}

/* Is theta in beta O_K?  Then xi = theta/beta lies in O_K, M(K,xi) = 0, and the
 * coset is beta O_K itself, whose least non-zero norm is |N(beta)| exactly --
 * so the search reports m = |N(beta)| and the ratio comes out at exactly 1,
 * looking like a certificate while proving nothing.  The certificate criterion (Proposition 1.3 of the paper) needs
 * alpha not in beta O_K, and this is where that hypothesis is enforced. */
static int in_beta_OK(GEN nf, GEN bta, GEN th)
{
  pari_sp av = avma;
  GEN q = nfdiv(nf, th, bta);
  long j; int z = 1;
  for (j = 1; j <= DEG; j++)
    if (typ(gel(q, j)) != t_INT) { z = 0; break; }
  set_avma(av); return z;
}

/* a fixed-point fraction as a double in [-1/2, 1/2) */
static inline double frac_to_double(u128 f)
{
  double d = (double)(unsigned long)(f >> 64) / 18446744073709551616.0;
  return (d >= 0.5) ? d - 1.0 : d;
}

#define KBOX 1            /* (unused: the slice search is a full qfcvp) */
#define NTH  512          /* hard ceiling on theta candidates per cube */

typedef struct {
  double W[DEG][DEG];     /* LLL-reduced basis, scaled metric: entries are O(1) */
  double Wi[DEG][DEG];    /* its inverse */
  GEN    T;               /* unimodular transform; entries reach e^{n*Reg} */
  GEN    Wg, Wi_g, G;     /* reduced basis, its inverse, and the Gram matrix */
  double c1, c2;
  long h[DEG];        /* half-widths of a complete search box */
  GEN is1, is2;       /* 1/s1, 1/s2: the slice scaling, precomputed */
  GEN A;              /* W^{-1} D_i: takes theta's embeddings straight to tgt */
  GEN qg;             /* A * ct(1): the step of tgt per unit increment of theta */
  GEN pg;             /* A * ct(theta0): the base of the current line */
  GEN qgb;            /* A * ct(delta): its step, delta being beta-sized */
  double (*we)[DEG];  /* W*e for every e of the box, nearest first */
  long *wee[DEG];     /* the e that produced it, for the exact branch */
  long nwe;
  u128 fq[DEG], fp[DEG];   /* their fractional parts, in units of 2^-128 */          /* scaled coordinates of theta are formed per theta */
} slice_t;

static void inv3(const double a[DEG][DEG], double o[DEG][DEG])
{
  double d = a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
           - a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
           + a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
  if (fabs(d) < 1e-300) d = 1e-300;
  o[0][0] =  (a[1][1]*a[2][2]-a[1][2]*a[2][1])/d;
  o[0][1] = -(a[0][1]*a[2][2]-a[0][2]*a[2][1])/d;
  o[0][2] =  (a[0][1]*a[1][2]-a[0][2]*a[1][1])/d;
  o[1][0] = -(a[1][0]*a[2][2]-a[1][2]*a[2][0])/d;
  o[1][1] =  (a[0][0]*a[2][2]-a[0][2]*a[2][0])/d;
  o[1][2] = -(a[0][0]*a[1][2]-a[0][2]*a[1][0])/d;
  o[2][0] =  (a[1][0]*a[2][1]-a[1][1]*a[2][0])/d;
  o[2][1] = -(a[0][0]*a[2][1]-a[0][1]*a[2][0])/d;
  o[2][2] =  (a[0][0]*a[1][1]-a[0][1]*a[1][0])/d;
}

int main(int argc, char **argv)
{
  long a[4] = {0,0,0,1};
  long n = 1;
  double cap = 0;
  int i, j;
  long prec;
  GEN pol, nf, bnf, eps, epsn, bta, NBg, M, Bg[DEG], capg;
  GEN bta_given = NULL;
  double E, NB;
  slice_t *sl;
  int nsl;
  /* theta has ~n*Reg/log10 digits per coordinate -- 182 for R_K = 419, so an
   * input line runs to ~550 characters and a 512-byte buffer truncated it
   * silently, splitting one theta across two fgets and evaluating neither.
   * That is exactly what made confirmation reject good survivors. */
  static char line[65536];
  long bestm_set = 0;
  int centres = 0, exact = 0, screen = 0;
  long rnd = 0, seed = 1, nsurv = 0, widen = 0;
  char *betastr = NULL;   /* --beta: verify against a modulus we are given */
  double hitsum = 0; long hitcnt = 0;
  long nskip = 0;
  long maxth = 8;   /* theta kept per cube: any point of an exceptional cube
                     * will do, so there is no gain in enumerating them all */
  GEN bestm = NULL, bestth = NULL;
  long bestc[4] = {0,0,0,0};   /* the cube that produced the best theta */

  for (i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--poly") && i+1 < argc)
      sscanf(argv[++i], "%ld,%ld,%ld,%ld", &a[0],&a[1],&a[2],&a[3]);
    else if (!strcmp(argv[i], "--n") && i+1 < argc) n = atol(argv[++i]);
    else if (!strcmp(argv[i], "--cap") && i+1 < argc) cap = atof(argv[++i]);
    else if (!strcmp(argv[i], "--centres")) centres = 1;
    else if (!strcmp(argv[i], "--exact")) exact = 1;
    else if (!strcmp(argv[i], "--screen")) screen = 1;
    else if (!strcmp(argv[i], "--widen") && i+1 < argc) widen = atol(argv[++i]);
    else if (!strcmp(argv[i], "--beta") && i+1 < argc) betastr = argv[++i];
    else if (!strcmp(argv[i], "--random") && i+1 < argc) rnd = atol(argv[++i]);
    else if (!strcmp(argv[i], "--seed") && i+1 < argc) seed = atol(argv[++i]);
    else if (!strcmp(argv[i], "--maxtheta") && i+1 < argc) maxth = atol(argv[++i]);
    else { fprintf(stderr, "see the header of %s\n", __FILE__); return 1; }
  }

  pari_init(1L<<28, 500000);
  paristack_setsize(1L<<26, 3L<<28);   /* 64 MB, ceiling 768 MB: a dozen
                                        * of these run side by side */
  pol = mkpoln(4, stoi(a[3]), stoi(a[2]), stoi(a[1]), stoi(a[0]));
  setvarn(pol, 0);
  bnf = Buchall(pol, nf_FORCE, DEFAULTPREC);
  nf  = bnf_get_nf(bnf);
  eps = algtobasis(nf, gel(bnf_get_fu(bnf), 1));

  /* An explicit modulus (--beta) verifies a Method-I certificate.  Nothing in
   * the unit-invariance and exhaustive-box lemmas (Lemmas 4.1-4.2 of the
   * paper) requires beta = eps^n - 1: all that is used is that the
   * coset alpha + beta O_K is stable under multiplication by some unit, and
   * for eps^t that says exactly beta | eps^t - 1.  Method I chooses its moduli
   * to be divisors of eps^t - 1 for small t, so every one of its certificates
   * has the invariance already.  The least such t is found here by exact
   * arithmetic -- it is a checked precondition, not an assumption -- and
   * verification then costs O(t*R_K) slices however large |N beta| is. */
  { GEN bgiven = NULL;
    if (betastr) {
      char s0[4096], s1[4096], s2[4096];
      long t;
      if (sscanf(betastr, "%4095s %4095s %4095s", s0, s1, s2) != 3) {
        fprintf(stderr, "--beta wants three integers\n"); return 2;
      }
      bgiven = mkcol3(gp_read_str(s0), gp_read_str(s1), gp_read_str(s2));
      { /* normalise eps at the default precision: only |eps^(1)| vs 1 matters */
        GEN Mt = nf_get_M(nf), e1 = gen_0;
        for (j = 0; j < DEG; j++)
          e1 = gadd(e1, gmul(gel(eps, j+1), gmael(Mt, j+1, 1)));
        if (gcmp(gabs(e1, DEFAULTPREC), gen_1) < 0) eps = nfinv(nf, eps);
      }
      for (t = 1; t <= 64; t++) {
        pari_sp av = avma;
        GEN d = nfadd(nf, nfpow(nf, eps, stoi(t)), gen_m1);
        GEN q = nfdiv(nf, d, bgiven);
        int ok = 1;
        for (j = 0; j < DEG; j++) if (typ(gel(q, j+1)) != t_INT) { ok = 0; break; }
        set_avma(av);
        if (ok) break;
      }
      if (t > 64) {
        fprintf(stderr, "no t <= 64 with beta | eps^t - 1: the coset is not "
                        "unit-stable and this method does not apply\n");
        return 2;
      }
      n = t;
      bta_given = gclone(bgiven);
      if (getenv("NE_VERB")) fprintf(stderr, "# beta | eps^%ld - 1\n", t);
    }
  }

  /* Precision is dictated by the dynamic range e^{n*Reg} of the region, plus
   * slack to separate lattice points inside a slice. */
  { double reg = rtodbl(bnf_get_reg(bnf));
    /* the region spans e^{n*Reg}, but the theta enumeration inverts the
     * multiplication matrix of beta, whose entries are e^{n*Reg/2}, so
     * intermediate quantities reach e^{3n*Reg/2}; budget for that */
    /* The 6x is for --centres, which inverts the multiplication matrix of
     * beta: its entries are e^{n*Reg/2} and the inverse cancels some
     * 3*n*Reg/2 worth of digits.  The slice geometry alone only has to
     * represent the range e^{n*Reg} of |x^(1)| with working slack, and
     * precision is the dominant cost of a trial, so do not pay for what this
     * path does not use. */
    prec = nbits2prec((long)((centres ? 6.0 : 2.0) * n * reg / M_LN2) + 512);
    if (bta_given) {
      /* beta is no longer eps^n - 1: its own coordinates can reach
       * |N beta|^{1/3} e^{R_K/2}, and the slice geometry has to represent them
       * against a range of e^{t R_K}. */
      long b = 0, j2;
      for (j2 = 0; j2 < DEG; j2++) {
        long e = expi(gel(bta_given, j2+1));
        if (e > b) b = e;
      }
      prec = nbits2prec((long)(2.0 * (n + 1) * reg / M_LN2) + 4*b + 1024);
    }
    if (getenv("NE_PREC")) prec = nbits2prec(atol(getenv("NE_PREC")));
    if (getenv("NE_VERB")) fprintf(stderr, "# prec = %ld bits\n", prec2nbits(prec));
    if (prec < DEFAULTPREC) prec = DEFAULTPREC;
    nf = nfnewprec_shallow(nf, prec);
  }
  M = nf_get_M(nf);

  /* use the unit that expands the real place */
  { GEN e1 = gen_0;
    for (j = 0; j < DEG; j++)
      e1 = gadd(e1, gmul(gel(eps, j+1), gmael(M, j+1, 1)));
    if (gcmp(gabs(e1, prec), gen_1) < 0) eps = nfinv(nf, eps);
  }
  epsn = nfpow(nf, eps, stoi(n));
  bta  = bta_given ? bta_given : nfadd(nf, epsn, gen_m1);
  NBg  = absi(nfnorm(nf, bta));
  NB   = gtodouble(NBg);
  if (cap <= 0) cap = screen ? NB : 2.0 * NB;
  /* keep the cap exact: dbltor() of a C double carries only 53 bits, and every
   * slice scaling derives from it, so the whole geometry silently drops to
   * 53-bit precision -- fatal once |x^(1)| reaches e^{n*Reg}. */
  capg = (cap == 2.0 * NB) ? shifti(NBg, 1)
       : (cap == NB)       ? NBg
       : gfloor(dbltor(cap));

  /* embeddings of the lattice basis beta*w_j, as PARI reals: these span
   * e^{n*Reg} and must not go through a double */
  { GEN b1 = gen_0, b2 = gen_0, e1 = gen_0;
    for (j = 0; j < DEG; j++) {
      b1 = gadd(b1, gmul(gel(bta, j+1), gmael(M, j+1, 1)));
      b2 = gadd(b2, gmul(gel(bta, j+1), gmael(M, j+1, 2)));
      e1 = gadd(e1, gmul(gel(epsn, j+1), gmael(M, j+1, 1)));
    }
    E = gtodouble(glog(gabs(e1, prec), prec));      /* log E, to size the walk */
    for (j = 0; j < DEG; j++) {
      GEN c = cgetg(DEG+1, t_COL);
      GEN p2 = gmul(b2, gmael(M, j+1, 2));
      gel(c,1) = gmul(b1, gmael(M, j+1, 1));
      gel(c,2) = greal(p2);
      gel(c,3) = gimag(p2);
      Bg[j] = c;
    }
  }
  if (!(E > 0)) { fprintf(stderr, "|eps^n| <= 1 at the real place\n"); return 2; }

  nsl = (int)ceil(E / M_LN2);
  if (nsl < 1) nsl = 1;
  sl = malloc(sizeof(slice_t) * nsl);

  for (i = 0; i < nsl; i++) {
    pari_sp av = avma;
    GEN s1 = gmul2n(gen_1, i+1);                    /* 2^(i+1) */
    GEN s2 = gsqrt(gdiv(capg, gmul2n(gen_1, i)), prec);
    GEN V = cgetg(DEG+1, t_MAT), W, T;
    for (j = 0; j < DEG; j++) {
      GEN c = cgetg(DEG+1, t_COL);
      gel(c,1) = gdiv(gel(Bg[j],1), s1);
      gel(c,2) = gdiv(gel(Bg[j],2), s2);
      gel(c,3) = gdiv(gel(Bg[j],3), s2);
      gel(V,j+1) = c;
    }
    T = qflll0(V, 0);                               /* transform: stays a GEN */
    W = RgM_mul(V, T);                              /* reduced: entries O(1) */
    { int k;
      for (j = 0; j < DEG; j++) for (k = 0; k < DEG; k++)
        sl[i].W[k][j] = gtodouble(gcoeff(W, k+1, j+1));
    }
    inv3(sl[i].W, sl[i].Wi);
    sl[i].T    = gclone(T);
    sl[i].Wg   = gclone(W);
    sl[i].Wi_g = gclone(RgM_inv(W));
    /* qfcvp wants an integral form: scale and round.  The reduced basis has
     * entries O(1), so 2^40 leaves ~12 digits of headroom, and the bound is
     * scaled to match. */
    { GEN Gr = RgM_mul(shallowtrans(W), W);
      sl[i].G = gclone(ground(RgM_Rg_mul(Gr, int2n(40)))); }
    /* A complete box, not a Babai guess.  A coset point of this slice has
     * |y0| in [1/2,1) and 2|y0|(y1^2+y2^2) < 1.05, hence y1^2+y2^2 <= 1.05 and
     * ||y|| <= sqrt(2.05); with u = tgt - W^{-1}y that bounds each coefficient
     * of u by ||row j of W^{-1}|| * sqrt(2.05), plus 1/2 because the box is
     * centred on round(tgt).  Enumerating that box is exhaustive, so the
     * minimum returned is the true one -- and a +-1 box is not: it silently
     * overlooks coset points, which inflates the reported m and can manufacture
     * a certificate that is not one.
     *
     * This bound is the one place where a floating-point error could cost
     * soundness rather than time, so it is not taken from the double-precision
     * inverse.  ||W^{-1}|| is read off the multiprecision inverse, and the
     * inverse is validated a posteriori: with R = I - (computed W^{-1}) W,
     * the true inverse satisfies ||W^{-1}|| <= ||computed|| / (1 - ||R||), and
     * the bound below is inflated by that factor.  If ||R|| is not comfortably
     * below 1 the slice is refused outright rather than trusted. */
    { GEN Res = RgM_sub(matid(DEG), RgM_mul(sl[i].Wi_g, W));
      double resn = 0, fac;
      int k;
      for (j = 0; j < DEG; j++) {
        double rs = 0;
        for (k = 0; k < DEG; k++) rs += fabs(gtodouble(gcoeff(Res, j+1, k+1)));
        if (rs > resn) resn = rs;
      }
      if (!(resn < 0.5)) {
        fprintf(stderr, "slice %d: inverse residual %.3g too large; refusing "
                        "to bound the search box\n", i, resn);
        return 3;
      }
      fac = 1.0 / (1.0 - resn);
      for (j = 0; j < DEG; j++) {
        double r2 = 0;
        for (k = 0; k < DEG; k++) {
          double v = gtodouble(gcoeff(sl[i].Wi_g, j+1, k+1));
          r2 += v * v;
        }
        /* sqrt(2.05)*||W^{-1}_j||*fac + 1/2, then a unit of margin for the
         * double evaluation of the row norm itself */
        r2 = sqrt(2.05 * r2) * fac * (1.0 + 1e-9) + 0.5;
        sl[i].h[j] = (long)ceil(r2) + 1 + widen;
        if (sl[i].h[j] < 1) sl[i].h[j] = 1;
        if (sl[i].h[j] > 60) {
          fprintf(stderr, "slice %d: box bound %ld exceeds 60; refusing\n",
                  i, sl[i].h[j]);
          return 3;
        }
      }
    }
    /* s2 = sqrt(cap/2^i) is a square root at full precision, and it depends
     * only on the slice -- recomputing it per trial cost more than everything
     * else in the inner loop put together. */
    sl[i].is1 = gclone(ginv(s1));
    sl[i].is2 = gclone(ginv(s2));
    /* c0 = D_i ct and tgt = W^{-1} c0, so tgt = (W^{-1} D_i) ct and c0 need
     * never be formed at all. */
    { GEN D = mkmat3(mkcol3(sl[i].is1, gen_0, gen_0),
                     mkcol3(gen_0, sl[i].is2, gen_0),
                     mkcol3(gen_0, gen_0, sl[i].is2));
      sl[i].A = gclone(RgM_mul(RgM_inv(W), D)); }
    /* W*e for the whole box depends only on the slice, so tabulate it once and
     * let a trial cost three subtractions per candidate instead of nine
     * multiplications.  Nearest-first, because the scan stops at the first
     * coset point below |N(beta)| and the near ones are where it is. */
    { long e0, e1_, e2, k = 0, h0 = sl[i].h[0], h1 = sl[i].h[1], h2 = sl[i].h[2];
      long nb = (2*h0+1) * (2*h1+1) * (2*h2+1);
      sl[i].nwe = nb;
      sl[i].we  = malloc(sizeof(double) * DEG * nb);
      for (j = 0; j < DEG; j++) sl[i].wee[j] = malloc(sizeof(long) * nb);
      for (e0 = -h0; e0 <= h0; e0++)
      for (e1_ = -h1; e1_ <= h1; e1_++)
      for (e2 = -h2; e2 <= h2; e2++, k++) {
        sl[i].we[k][0] = sl[i].W[0][0]*e0 + sl[i].W[0][1]*e1_ + sl[i].W[0][2]*e2;
        sl[i].we[k][1] = sl[i].W[1][0]*e0 + sl[i].W[1][1]*e1_ + sl[i].W[1][2]*e2;
        sl[i].we[k][2] = sl[i].W[2][0]*e0 + sl[i].W[2][1]*e1_ + sl[i].W[2][2]*e2;
        sl[i].wee[0][k] = e0; sl[i].wee[1][k] = e1_; sl[i].wee[2][k] = e2;
      }
      { long a, b2;                        /* insertion sort by |W e| */
        for (a = 1; a < nb; a++) {
          double ka[DEG]; long ke[DEG];
          double na;
          for (j = 0; j < DEG; j++) { ka[j] = sl[i].we[a][j]; ke[j] = sl[i].wee[j][a]; }
          na = ka[0]*ka[0] + ka[1]*ka[1] + ka[2]*ka[2];
          for (b2 = a - 1; b2 >= 0; b2--) {
            double nb2 = sl[i].we[b2][0]*sl[i].we[b2][0]
                       + sl[i].we[b2][1]*sl[i].we[b2][1]
                       + sl[i].we[b2][2]*sl[i].we[b2][2];
            if (nb2 <= na) break;
            for (j = 0; j < DEG; j++) {
              sl[i].we[b2+1][j] = sl[i].we[b2][j];
              sl[i].wee[j][b2+1] = sl[i].wee[j][b2];
            }
          }
          for (j = 0; j < DEG; j++) { sl[i].we[b2+1][j] = ka[j]; sl[i].wee[j][b2+1] = ke[j]; }
        }
      }
    }
    /* theta -> theta + 1 shifts tgt by A*ct(1), and ct(1) = (1,1,0) */
    { GEN q = RgM_RgC_mul(sl[i].A, mkcol3(gen_1, gen_1, gen_0));
      sl[i].qg = gclone(q);
      for (j = 0; j < DEG; j++) sl[i].fq[j] = frac_u128(gel(q, j+1)); }
    sl[i].c1 = gtodouble(gdiv(gen_1, s1));          /* unused directly; kept */
    sl[i].c2 = gtodouble(gdiv(gen_1, s2));
    set_avma(av);
  }

  if (screen && rnd) {
    /* Fast screening.
     *
     * Along theta = theta0 + t the target of every slice is linear in t,
     * tgt_i = p_i + t q_i, and the geometry only ever uses tgt modulo 1: the
     * integer part is discarded by the Babai rounding.  So the e^{n Reg}-sized
     * multiprecision work is done once per line, and each trial costs three
     * 128-bit multiply-adds per slice instead of nine 1800-bit multiplications.
     *
     * Exact arithmetic is not abandoned, it is confined.  Rejecting theta needs
     * |N(x)| < |N(beta)| to be *true*, and the scaled coordinates give
     * |N(x)|/|N(beta)| to a relative 1e-15, so a rejection below 1 - 1e-9 is
     * certain and is taken in doubles; only the vanishingly rare candidate
     * inside that band is decided by an exact integer norm.  A survivor is
     * still only a candidate, and is re-evaluated in full before it is
     * recorded.
     *
     * The line is short (LINE trials) and restarted from a fresh random
     * theta0, because theta -> theta + 1 moves xi by 1/beta, which is O(1) at
     * the complex place but e^{-n Reg} at the real one: a single long line
     * would explore a slice of the torus rather than the whole of it.
     */
    const long LINE = 4096;
    long done = 0, t, i, j;
    setrand(stoi(seed));
    while (done < rnd) {
      pari_sp avblk = avma;
      long L = (rnd - done < LINE) ? rnd - done : LINE;
      GEN th0, dlt, ct0, ctd, mb = zk_multable(nf, bta);
      GEN v = mkcol3(randomi(int2n(64)), randomi(int2n(64)), randomi(int2n(64)));
      GEN w = mkcol3(randomi(int2n(64)), randomi(int2n(64)), randomi(int2n(64)));
      th0 = ground(gdiv(RgM_RgC_mul(mb, v), int2n(64)));
      dlt = ground(gdiv(RgM_RgC_mul(mb, w), int2n(64)));
      { GEN c1 = gen_0, c2 = gen_0, e1 = gen_0, e2 = gen_0;
        for (j = 0; j < DEG; j++) {
          c1 = gadd(c1, gmul(gel(th0, j+1), gmael(M, j+1, 1)));
          c2 = gadd(c2, gmul(gel(th0, j+1), gmael(M, j+1, 2)));
          e1 = gadd(e1, gmul(gel(dlt, j+1), gmael(M, j+1, 1)));
          e2 = gadd(e2, gmul(gel(dlt, j+1), gmael(M, j+1, 2)));
        }
        ct0 = mkcol3(c1, greal(c2), gimag(c2));
        ctd = mkcol3(e1, greal(e2), gimag(e2));
      }
      for (i = 0; i < nsl; i++) {
        sl[i].pg = RgM_RgC_mul(sl[i].A, ct0);
        sl[i].qgb = RgM_RgC_mul(sl[i].A, ctd);
        for (j = 0; j < DEG; j++) {
          sl[i].fp[j] = frac_u128(gel(sl[i].pg, j+1));
          sl[i].fq[j] = frac_u128(gel(sl[i].qgb, j+1));
        }
      }
      for (t = 0; t < L; t++) {
        int rejected = 0;
        for (i = 0; i < nsl && !rejected; i++) {
          double d0 = frac_to_double(sl[i].fp[0] + mul128(sl[i].fq[0], t));
          double d1 = frac_to_double(sl[i].fp[1] + mul128(sl[i].fq[1], t));
          double d2 = frac_to_double(sl[i].fp[2] + mul128(sl[i].fq[2], t));
          double b0 = sl[i].W[0][0]*d0 + sl[i].W[0][1]*d1 + sl[i].W[0][2]*d2;
          double b1 = sl[i].W[1][0]*d0 + sl[i].W[1][1]*d1 + sl[i].W[1][2]*d2;
          double b2 = sl[i].W[2][0]*d0 + sl[i].W[2][1]*d1 + sl[i].W[2][2]*d2;
          long k, nwe = sl[i].nwe;
          for (k = 0; k < nwe && !rejected; k++) {
            double y0 = b0 - sl[i].we[k][0];
            double y1 = b1 - sl[i].we[k][1];
            double y2 = b2 - sl[i].we[k][2];
            long e0 = sl[i].wee[0][k], e1_ = sl[i].wee[1][k], e2 = sl[i].wee[2][k];
            double q = 2.0 * fabs(y0) * (y1*y1 + y2*y2);
            if (q >= 1.05) continue;              /* outside the search region */
            if (q < 1.0 - 1e-9) {                 /* certain: |N| < |N(beta)| */
              rejected = 1; hitsum += i + 1; hitcnt++; break;
            }
            { /* the boundary band: decide it exactly */
              pari_sp av3 = avma;
              GEN r = ground(RgC_add(sl[i].pg, RgC_Rg_mul(sl[i].qgb, stoi(t))));
              GEN u = mkcol3(addis(gel(r,1), e0), addis(gel(r,2), e1_),
                             addis(gel(r,3), e2));
              GEN th = RgC_add(th0, RgC_Rg_mul(dlt, stoi(t)));
              GEN g = gsub(th, nfmul(nf, bta, RgM_RgC_mul(sl[i].T, u)));
              if (!gequal0(g) && cmpii(absi(nfnorm(nf, g)), NBg) < 0) {
                rejected = 1; hitsum += i + 1; hitcnt++;
              }
              set_avma(av3);
            }
          }
        }
        if (!rejected) {
          GEN th = RgC_add(th0, RgC_Rg_mul(dlt, stoi(t)));
          if (in_beta_OK(nf, bta, th)) continue;
          { char *sn = GENtostr(NBg), *t0 = GENtostr(gel(th,1)),
               *t1 = GENtostr(gel(th,2)), *t2 = GENtostr(gel(th,3));
          printf("SURVIVOR %s %s %s %s\n", sn, t0, t1, t2);
          fflush(stdout);
          pari_free(sn); pari_free(t0); pari_free(t1); pari_free(t2);
          nsurv++; }
        }
      }
      done += L;
      set_avma(avblk);
    }
    { double lam = (hitcnt && hitsum > 0) ? nsl * hitcnt / hitsum : 0.0;
      printf("# screened %ld theta, %ld survivors, nsl %d, lambda %.4f\n",
             rnd, nsurv, nsl, lam); }
    return nsurv ? 0 : 3;
  }

  if (rnd) setrand(stoi(seed));
  { long iter = 0;
  while (rnd ? (iter++ < rnd) : (fgets(line, sizeof(line), stdin) != NULL)) {
    pari_sp av = avma;
    volatile long skipped_local = 0;
    long t[DEG], dep;
    pari_CATCH(CATCH_ALL) {
      if (!nskip) {
        GEN e = pari_err_last();
        char *m = pari_err2str(e);
        fprintf(stderr, "# first error: %s\n", m ? m : "?");
        if (m) pari_free(m);
      }
      nskip++; set_avma(av); continue;
    }
    pari_TRY {
    GEN th, ctg[DEG], thlist[NTH];
    int nth = 1, ith;
    GEN mymin = NULL;
    volatile int rejected = 0;
    if (rnd) {
      /* xi uniform in the torus (K tensor R)/O_K, theta = round(beta*xi).
       * Sampling xi rather than theta is what makes it uniform: the residues
       * of O_K/beta correspond to a fundamental cell of beta^{-1}O_K, and the
       * cell is extremely skew (beta has embeddings of size e^{n*Reg/2}), so
       * picking theta coordinates from a box would sample it very unevenly. */
      GEN v = mkcol3(randomi(int2n(64)), randomi(int2n(64)), randomi(int2n(64)));
      thlist[0] = ground(gdiv(RgM_RgC_mul(zk_multable(nf, bta), v), int2n(64)));
      nth = 1;
    } else if (centres) {
      /* The integer theta with theta/beta inside the cube are the lattice
       * points of (1/delta)*Mbeta^{-1} Z^3 in a unit box about xi/delta; there
       * are about N(beta)*delta^3 of them, but spread over a coordinate range
       * of size |beta|*delta ~ e^{n*Reg/2}*delta, so rounding beta*xi and
       * nudging by +-1 misses them entirely.  Reduce and enumerate instead. */
      GEN xi, A = NULL, T = NULL, W, Wi, tgt = NULL, m0 = NULL;
      int simple;
      int d0, d1, d2;
      if (sscanf(line, "%ld %ld %ld %ld", &t[0], &t[1], &t[2], &dep) != 4) continue;
      xi = mkcol3(gdiv(stoi(t[0]), int2n(dep)),
                  gdiv(stoi(t[1]), int2n(dep)),
                  gdiv(stoi(t[2]), int2n(dep)));
      /* How many theta have theta/beta inside the cube?  About
       * N(beta)*delta^3.  When that is large the lattice beta^{-1}O_K is far
       * finer than the cube and plain rounding lands well inside it; the LLL
       * apparatus is only needed when the count is O(1), which happens for
       * small regulators.  Taking the general path always is not merely
       * wasteful: for Reg ~ 200, Mbeta^{-1} has 92-digit entries and its LLL
       * transform 184-digit ones, so forming the reduced basis cancels some
       * 300 digits and the rounding loses all meaning. */
      { GEN volg = gdiv(NBg, powiu(gen_2, 3*(dep-1)));
        simple = (gcmp(volg, powuu(10, 6)) > 0);
      }
      if (!simple) {
        GEN den, Ainv = RgM_inv(zk_multable(nf, bta));
        GEN Aint = Q_remove_denom(Ainv, &den);
        T = qflll0(Aint, 0);
        A = RgM_gtofp(RgM_Rg_mul(Ainv, int2n(dep-1)), prec);
      }
      if (!simple) { W = RgM_mul(A, T); Wi = RgM_inv(W); }
      nth = 0;
      if (simple) {                       /* fine lattice: round and nudge */
        GEN t0 = ground(nfmul(nf, bta, xi));
        int e0, e1_, e2;
        for (e0 = -1; e0 <= 1; e0++)
        for (e1_ = -1; e1_ <= 1; e1_++)
        for (e2 = -1; e2 <= 1; e2++)
          if (nth < maxth)
            thlist[nth++] = gadd(t0, mkcol3(stoi(e0), stoi(e1_), stoi(e2)));
        goto have_theta;
      }
      tgt = RgC_Rg_mul(xi, int2n(dep-1));
      m0 = RgM_RgC_mul(Wi, tgt);
      for (d0 = -3; d0 <= 3 && nth < maxth && nth < NTH; d0++)
      for (d1 = -3; d1 <= 3 && nth < maxth && nth < NTH; d1++)
      for (d2 = -3; d2 <= 3 && nth < maxth && nth < NTH; d2++) {
        GEN u = mkcol3(gadd(ground(gel(m0,1)), stoi(d0)),
                       gadd(ground(gel(m0,2)), stoi(d1)),
                       gadd(ground(gel(m0,3)), stoi(d2)));
        GEN cand = RgM_RgC_mul(T, u);
        GEN chk = RgM_RgC_mul(A, cand);        /* must land in the unit box */
        int ok = 1, q;
        for (q = 0; q < DEG; q++)
          if (gcmp(gabs(gsub(gel(chk,q+1), gel(tgt,q+1)), prec),
                   gdivgs(gen_1, 2)) > 0) { ok = 0; break; }
        if (ok) thlist[nth++] = cand;
      }
    have_theta:
      if (!nth) { set_avma(av); continue; }
    } else {
      { static char w0[16384], w1[16384], w2[16384];
        if (sscanf(line, "%16383s %16383s %16383s", w0, w1, w2) != 3) continue;
        thlist[0] = mkcol3(gp_read_str(w0), gp_read_str(w1), gp_read_str(w2));
      }
    }
    for (ith = 0; ith < nth; ith++) {
    th = thlist[ith];
    if (in_beta_OK(nf, bta, th)) continue;
    mymin = NULL;
    { GEN c1 = gen_0, c2 = gen_0;
      for (j = 0; j < DEG; j++) {
        c1 = gadd(c1, gmul(gel(th, j+1), gmael(M, j+1, 1)));
        c2 = gadd(c2, gmul(gel(th, j+1), gmael(M, j+1, 2)));
      }
      ctg[0] = c1; ctg[1] = greal(c2); ctg[2] = gimag(c2);
    }
    for (i = 0; i < nsl; i++) {
      pari_sp av2 = avma;
      GEN c0 = NULL, tgt, vecs;
      long nv, iv;
      /* The admissible region of the slice lies in the ball of squared radius
       * 2 about the target.  --exact enumerates it completely (Fincke-Pohst);
       * the default takes the 27 points of a +-1 box about the Babai point.
       * On every case tested the two return the same minimum, and qfcvp costs
       * ~2 ms against ~2 us -- decisive when it is called once per slice per
       * theta and a large regulator means 300+ slices. */
      tgt = RgM_RgC_mul(sl[i].A, mkcol3(ctg[0], ctg[1], ctg[2]));
      if (exact) {
        c0 = mkcol3(gmul(ctg[0], sl[i].is1), gmul(ctg[1], sl[i].is2),
                    gmul(ctg[2], sl[i].is2));
        vecs = gel(qfcvp0(sl[i].G, tgt, mulii(stoi(3), int2n(40)), NULL, 0), 3);
        nv = lg(vecs) - 1;
        for (iv = 1; iv <= nv; iv++) {
          pari_sp av3 = avma;
          GEN u = gel(vecs, iv);
          GEN y = RgC_sub(c0, RgM_RgC_mul(sl[i].Wg, u));
          double y0 = gtodouble(gel(y,1)), y1 = gtodouble(gel(y,2)),
                 y2 = gtodouble(gel(y,3));
          if (!(2.0 * fabs(y0) * (y1*y1 + y2*y2) < 1.05))
            { set_avma(av3); continue; }
          { GEN eta = RgM_RgC_mul(sl[i].T, u);
            GEN g = gsub(th, nfmul(nf, bta, eta));
            if (!gequal0(g)) {
              GEN v = absi(nfnorm(nf, g));
              if (signe(v) && (!mymin || cmpii(v, mymin) < 0)) {
                if (mymin) gunclone(mymin);
                mymin = gclone(v); set_avma(av3);
                if (screen && cmpii(mymin, NBg) < 0)
                  { rejected = 1; hitsum += i + 1; hitcnt++; break; }
                continue;
              }
            }
          }
          set_avma(av3);
        }
      } else {
        /* c0 is of size e^{n*Reg} while y = c0 - W*u is O(1), so y can only be
         * formed at full precision -- but only once.  Every other candidate of
         * the box is y_base - W*e with e in {-h..h}^3, which is O(1) minus
         * O(1): pure double arithmetic, and the prefilter is all it feeds.
         * Anything that passes is still decided by an exact integer norm, so
         * the doubles never decide the answer, only where to look. */
        long e0, e1_, e2;
        long h0 = sl[i].h[0], h1 = sl[i].h[1], h2 = sl[i].h[2];
        /* y_base = c0 - W r = W (tgt - r), and tgt - r is O(1) by
         * construction: the e^{n Reg}-sized cancellation is three
         * subtractions, and the reduced basis W is O(1), so the product runs
         * in machine doubles like the rest of the box. */
        GEN r = ground(tgt);
        GEN dl = RgC_sub(tgt, r);
        double d0 = gtodouble(gel(dl,1)), d1 = gtodouble(gel(dl,2)),
               d2 = gtodouble(gel(dl,3));
        double b0 = sl[i].W[0][0]*d0 + sl[i].W[0][1]*d1 + sl[i].W[0][2]*d2;
        double b1 = sl[i].W[1][0]*d0 + sl[i].W[1][1]*d1 + sl[i].W[1][2]*d2;
        double b2 = sl[i].W[2][0]*d0 + sl[i].W[2][1]*d1 + sl[i].W[2][2]*d2;
        pari_sp av4 = avma;
        for (e0 = -h0; e0 <= h0 && !rejected; e0++)
        for (e1_ = -h1; e1_ <= h1 && !rejected; e1_++)
        for (e2 = -h2; e2 <= h2 && !rejected; e2++) {
          double y0 = b0 - (sl[i].W[0][0]*e0 + sl[i].W[0][1]*e1_ + sl[i].W[0][2]*e2);
          double y1 = b1 - (sl[i].W[1][0]*e0 + sl[i].W[1][1]*e1_ + sl[i].W[1][2]*e2);
          double y2 = b2 - (sl[i].W[2][0]*e0 + sl[i].W[2][1]*e1_ + sl[i].W[2][2]*e2);
          if (!(2.0 * fabs(y0) * (y1*y1 + y2*y2) < 1.05)) continue;
          { GEN u = mkcol3(addis(gel(r,1), e0), addis(gel(r,2), e1_),
                           addis(gel(r,3), e2));
            GEN eta = RgM_RgC_mul(sl[i].T, u);
            GEN g = gsub(th, nfmul(nf, bta, eta));
            if (!gequal0(g)) {
              GEN v = absi(nfnorm(nf, g));
              if (signe(v) && (!mymin || cmpii(v, mymin) < 0)) {
                if (mymin) gunclone(mymin);
                mymin = gclone(v);
                if (screen && cmpii(mymin, NBg) < 0)
                  { rejected = 1; hitsum += i + 1; hitcnt++; }
              }
            }
          }
          set_avma(av4);
        }
      }
      set_avma(av2);
      if (rejected) break;
    }
    if (screen && !rejected) {
      /* survivor: no coset point below N(beta) was found, so m/N(beta) >= 1.
       * The default box enumeration can miss points, so this is a candidate,
       * not yet a proof -- re-run it with --exact, and then independently. */
      char *sn = GENtostr(NBg), *t0 = GENtostr(gel(th,1)),
           *t1 = GENtostr(gel(th,2)), *t2 = GENtostr(gel(th,3));
      printf("SURVIVOR %s %s %s %s\n", sn, t0, t1, t2);
      fflush(stdout);
      pari_free(sn); pari_free(t0); pari_free(t1); pari_free(t2);
      nsurv++;
    }
    if (!screen && mymin && (!bestm_set || cmpii(mymin, bestm) > 0)) {
      if (bestm_set) { gunclone(bestm); gunclone(bestth); }
      bestm = gclone(mymin); bestth = gclone(th); bestm_set = 1;
      if (centres) { bestc[0]=t[0]; bestc[1]=t[1]; bestc[2]=t[2]; bestc[3]=dep; }
    }
    if (mymin) { gunclone(mymin); mymin = NULL; }
    }
    } pari_ENDCATCH;
    (void)skipped_local;
    set_avma(av);
  } }
  if (nskip) fprintf(stderr, "# %ld input lines skipped (PARI error)\n", nskip);
  if (screen) {
    { double lam = (hitcnt && hitsum > 0) ? nsl * hitcnt / hitsum : 0.0;
      printf("# screened %ld theta, %ld survivors, nsl %d, lambda %.4f\n",
             rnd, nsurv, nsl, lam); }
    return nsurv ? 0 : 3;
  }

  { char *sn = GENtostr(NBg);
    if (bestm_set) {
      char *sm = GENtostr(bestm);
      char *t0 = GENtostr(gel(bestth,1)), *t1 = GENtostr(gel(bestth,2)),
           *t2 = GENtostr(gel(bestth,3));
      printf("%s %s %s %s %s", sm, sn, t0, t1, t2);
      if (centres) printf("  # cube %ld %ld %ld %ld",
                          bestc[0], bestc[1], bestc[2], bestc[3]);
      printf("\n");
      pari_free(sm); pari_free(t0); pari_free(t1); pari_free(t2);
    } else printf("0 %s 0 0 0\n", sn);
    pari_free(sn);
  }
  return 0;
}
