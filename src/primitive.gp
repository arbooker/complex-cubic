default(realprecision, 500);
{
inf = fileopen(FIN); ouf = fileopen(FOUT, "w"); nred = 0; nfail = 0;
while(1,
  L = filereadstr(inf); if(L == 0, break);
  v = strsplit(L, ";");
  pol = Pol(Vecrev(apply(x->eval(x), strsplit(v[2],","))));
  den = eval(v[3]);
  th  = Pol(Vecrev(apply(x->eval(x), strsplit(v[4],","))))/den;
  be  = Pol(Vecrev(apply(x->eval(x), strsplit(v[5],","))))/den;
  nf  = nfinit(pol);
  ta  = nfalgtobasis(nf, th); ba = nfalgtobasis(nf, be);
  g   = idealadd(nf, ta, ba);
  if(idealnorm(nf, g) != 1,
    bnf = bnfinit(pol, 1);
    pr  = iferr(bnfisprincipal(bnf, g, 1), E, 0);
    if(type(pr) == "t_VEC" && #pr >= 2 && type(pr[2]) == "t_COL" && #pr[2] == 3,
      dd = pr[2];
      qt = iferr(nfbasistoalg(nf, ta)/nfbasistoalg(nf, dd), E, 0);
      qb = iferr(nfbasistoalg(nf, ba)/nfbasistoalg(nf, dd), E, 0);
      t2 = iferr(nfalgtobasis(nf, qt), E, 0);
      b2 = iferr(nfalgtobasis(nf, qb), E, 0);
      if(type(t2) == "t_COL" && type(b2) == "t_COL"
         && denominator(t2) == 1 && denominator(b2) == 1,
        uu = bnf.fu[1];
        lu = log(abs(nfeltembed(nf, uu)[1]));
        lb = log(abs(nfeltembed(nf, b2)[1]));
        nn = abs(nfeltnorm(nf, b2));
        k0 = round((log(nn)/3 - lb)/lu);
        kk = k0; best = -1;
        for(j = -3, 3,
          dl = (lb + (k0+j)*lu) - log(nn)/3;
          tv = exp(2*dl) + 2*exp(-dl);
          if(best < 0 || tv < best, best = tv; kk = k0+j);
        );
        if(kk != 0,
          up = nfeltpow(nf, uu, kk);
          b2 = nfeltmul(nf, b2, up);
          t2 = nfeltmul(nf, t2, up);
        );
        nred++;
        pt = Vecrev(lift(nfbasistoalg(nf, t2)));
        pb = Vecrev(lift(nfbasistoalg(nf, b2)));
        pt = concat(pt, vector(3 - #pt));
        pb = concat(pb, vector(3 - #pb));
        filewrite(ouf, Str(v[1], ";", pt[1], ",", pt[2], ",", pt[3],
                           ";", pb[1], ",", pb[2], ",", pb[3]));
      , nfail++);
    , nfail++);
  );
);
fileclose(inf); fileclose(ouf);
print("reduced ", nred, " failed ", nfail);
}
quit
