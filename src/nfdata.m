/* nfdata.m -- the Magma counterpart of src/nfdata.gp.
 *
 * Emits the same per-field record: the integral basis of O_K and a unit, in
 * the same format, so that the two can be compared.  The certifier assumes
 * exactly one thing it cannot check -- that the basis spans the *maximal*
 * order -- and that assumption is the reason for having two sources for it.
 *
 * A basis of O_K is not unique, so agreement does not mean identical output:
 * the test is that the change of basis between the two is an integer matrix of
 * determinant +-1, i.e. that they span the same lattice.  src/cmp_nfdata.py
 * does that comparison.
 *
 *   magma -b csv:=data/fields.csv src/nfdata.m
 */
/* Magma wraps output at the terminal width by default, and eps runs to
   hundreds of digits here, so records would be split across lines. */
SetAutoColumns(false);
SetColumns(0);

csvname := csv;
fh := Open(csvname, "r");
_ := Gets(fh);                      // header
while true do
    line := Gets(fh);
    if IsEof(line) then break; end if;
    p := Split(line, ",");
    if #p lt 6 then continue; end if;
    d  := StringToInteger(p[2]);
    co := [StringToInteger(p[3]), StringToInteger(p[4]),
           StringToInteger(p[5]), StringToInteger(p[6])];
    /* not P<x> := ... : the angle-bracket form inside a loop hangs Magma */
    P := PolynomialRing(Rationals());
    x := P.1;
    f := co[1] + co[2]*x + co[3]*x^2 + co[4]*x^3;
    K := NumberField(f);
    O := MaximalOrder(K);
    B := Basis(O);
    /* common denominator of the basis in the power basis of K */
    den := 1;
    for b in B do
        c := Eltseq(K ! b);
        for e in c do den := Lcm(den, Denominator(e)); end for;
    end for;
    ent := [];
    for b in B do
        c := Eltseq(K ! b);
        for l in [1..3] do Append(~ent, Numerator(c[l] * den)); end for;
    end for;
    U, mU := UnitGroup(O);
    eps := O ! mU(U.2);                       // a unit of infinite order
    ec := Eltseq(O ! eps);                    // coordinates in the basis B
    printf "%o|%o,%o,%o,%o|%o|", d, co[1], co[2], co[3], co[4], den;
    for i in [1..9] do
        printf "%o", ent[i];
        if i lt 9 then printf ","; end if;
    end for;
    printf "|%o,%o,%o|%o\n", ec[1], ec[2], ec[3], Discriminant(O);
end while;

quit;
