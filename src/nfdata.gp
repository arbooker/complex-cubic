/* nfdata.gp -- emit the algebraic input the certifier needs, once per field.
 *
 * The certifier (src/certify) does not compute algebraic number theory.  It is
 * handed the little that cannot be derived from f by exact integer arithmetic:
 *
 *     the integral basis of O_K, and a unit eps.
 *
 * Everything else it checks itself.  Splitting it out this way has two points.
 * It removes the last dependency of the certifier on any computer algebra
 * system, so the program that decides the result is FLINT + Arb + fplll and
 * nothing else; and it makes the one input that is *assumed* rather than
 * checked -- that the basis spans the maximal order -- a small explicit file
 * that a second, independent system can be asked to produce (src/nfdata.m does
 * the same job in Magma).  Disagreement between the two would be visible
 * immediately; agreement is worth having precisely because the two compute it
 * by different means.
 *
 * The unit need not be fundamental: the certifier checks |N eps| = 1 and
 * beta | eps^n - 1 exactly, and needs nothing else about it.
 *
 * Output, one line per field:
 *   disc | a0,a1,a2,a3 | den | b00,b01,b02,b10,... | e0,e1,e2 | d_K
 * meaning omega_j = (sum_k b_jk x^k)/den and eps = sum_j e_j omega_j.
 *
 *   echo 'read("src/nfdata.gp"); ne_nfdata("data/fields.csv")' | gp -q
 */
{
ne_nfdata(csv) =
  my(fh, line, first = 1);
  fh = fileopen(csv, "r");
  while(line = filereadstr(fh),
    if(first, first = 0; next);
    my(p = strsplit(line, ","));
    if(#p < 6, next);
    my(d = eval(p[2]), co = vector(4, i, eval(p[2+i])));
    my(f = Pol(Vecrev(co), 'x));
    my(nf = nfinit(f), bnf = bnfinit(f, 1));
    my(zk = nf.zk, den = 1);
    /* the lcm of the coefficient denominators, taken coefficient by
       coefficient: denominator() of a t_POL does not give it */
    for(j = 1, 3, for(l = 0, 2, den = lcm(den, denominator(polcoef(zk[j], l)))));
    my(B = vector(9), k = 0);
    for(j = 1, 3,
      my(q = zk[j] * den);
      for(l = 0, 2, k++; B[k] = polcoef(q, l)));
    my(e = nfalgtobasis(bnf, bnf.fu[1]));
    print(d, "|", co[1], ",", co[2], ",", co[3], ",", co[4],
          "|", den, "|", B[1], ",", B[2], ",", B[3], ",", B[4], ",", B[5],
          ",", B[6], ",", B[7], ",", B[8], ",", B[9],
          "|", e[1], ",", e[2], ",", e[3], "|", nf.disc);
  );
  fileclose(fh);
}
