// sing.m : singularities of a plane sextic given by its 28 coefficients (linesearch order) in a file "[c1, ..., c28]".
// usage: magma -b file:=out/model_A_red_coeffs.txt [out:=out/model_A_sing.gp] sing.m
// The singular subscheme is decomposed into its primary components; for each: the reduced point (rational, or of degree
// d > 1), the local length tau (Tjurina number: 1 for an ordinary node, n for an A_n/D_n/E_n simple singularity), the
// multiplicity of the curve and the tangent cone.  A double point with tau = n is an A_n singularity of the sextic, i.e. an
// A_n singularity of the double plane.  If out is given, writes for finalize.gp:
//   singpts = [[a,b,c], ...] (the rational singular points, primitive integer coordinates); singtau = [tau, ...];
//   singmult = [multiplicity, ...]; singdeg = <degree of the singular subscheme>; singnonrat = [[degree, tau], ...].
c := eval Read(file);
P2<x0,x1,x2> := ProjectiveSpace(Rationals(), 2);
mons := [];
for i in [6..0 by -1] do for j in [6-i..0 by -1] do k := 6-i-j; Append(~mons, x0^i*x1^j*x2^k); end for; end for;
S := &+[ c[n]*mons[n] : n in [1..28] ];
C := Curve(P2, S);
printf "file %o: irreducible %o, reduced %o\n", file, IsIrreducible(C), IsReduced(C);
Sg := SingularSubscheme(C);
deg := Degree(Sg);
printf "singular subscheme: degree %o, dimension %o\n", deg, Dimension(Sg);
Q, P := PrimaryDecomposition(DefiningIdeal(Sg));
ptlist := []; taulist := []; multlist := []; nonrat := [];
for k in [1..#Q] do
  Zp := Scheme(P2, P[k]); Zq := Scheme(P2, Q[k]);
  if IsEmpty(Zp) then continue; end if;      // the irrelevant component (x0, x1, x2) of the homogeneous ideal
  dp := Degree(Zp); tau := Degree(Zq);
  if dp eq 1 then
    p := RationalPoints(Zp)[1];
    m := Multiplicity(C, p);
    cs := Coordinates(p); den := LCM([Denominator(x) : x in cs]); cs := [Integers() | x*den : x in cs]; g := GCD(cs); cs := [x div g : x in cs];
    if cs[1] lt 0 or (cs[1] eq 0 and cs[2] lt 0) or (cs[1] eq 0 and cs[2] eq 0 and cs[3] lt 0) then cs := [-x : x in cs]; end if;
    typ := m eq 2 select Sprintf("A%o", tau) else (m eq 3 select Sprintf("D/E%o (multiplicity 3)", tau) else Sprintf("multiplicity %o, tau %o", m, tau));
    printf "  rational singular point %o : multiplicity %o, tau %o, ordinary node? %o, type %o, tangent cone %o\n", cs, m, tau, IsNode(C, p), typ, DefiningPolynomial(TangentCone(C, p));
    Append(~ptlist, cs); Append(~taulist, tau); Append(~multlist, m);
  else
    printf "  non-rational singular point of degree %o (tau %o per conjugate: local length %o), ideal %o\n", dp, tau div dp, tau, Basis(P[k]);
    Append(~nonrat, [dp, tau]);
  end if;
end for;
printf "summary: %o rational singular points with tau %o (multiplicities %o), non-rational components %o\n", #ptlist, taulist, multlist, nonrat;
if assigned out then
  // one statement per line (gp's read cannot parse a sequence printed by Magma over several lines)
  ptstr := "[" cat Join([Sprintf("[%o,%o,%o]", p[1], p[2], p[3]) : p in ptlist], ",") cat "]";
  nrstr := "[" cat Join([Sprintf("[%o,%o]", q[1], q[2]) : q in nonrat], ",") cat "]";
  lst := func<L | "[" cat Join([Sprint(x) : x in L], ",") cat "]">;
  PrintFile(out, Sprintf("singpts = %o;\nsingtau = %o;\nsingmult = %o;\nsingdeg = %o;\nsingnonrat = %o;", ptstr, lst(taulist), lst(multlist), deg, nrstr) : Overwrite := true);
  printf "written %o\n", out;
end if;
