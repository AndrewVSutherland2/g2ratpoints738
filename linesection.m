// linesection.m -- among the sextics of data/models.txt, C_738 is a line section of exactly one, its own model S_3
// (the row t0453), by exactly one line over an algebraic closure of Q; on the others no line, rational or not,
// cuts out a curve isomorphic to C_738.
//
//   magma -b linesection.m                    all models of data/models.txt (a few hours)
//   magma -b only:=p100 linesection.m         one model (a few seconds);  only:=t0453 is the model S_3
//   magma -b lo:=1 hi:=40 linesection.m       models number lo..hi of the file (for running several processes)
// linesection.log is the output of 16 such pieces of 149 models each, one after the other.
//
// Statement proved for a sextic S: the number of lines l of the plane, over Qbar, such that the binary sextic S|l has nonzero
// discriminant and the same Igusa-Clebsch invariants as the sextic f of C_738, as a point of the weighted projective space
// P(2,4,6,10).  Two binary sextics with nonzero discriminant have the same point there if and only if the genus 2 curves
// y^2 = (sextic) are isomorphic over Qbar (Igusa), so that number is the number of lines whose section is Qbar-isomorphic to
// C_738.  A section over Q that is isomorphic to C_738 over Q is in particular counted.
//
// The invariants I2, I4, I6, I10 of a section are polynomials in the parameters of the line, and I10 is its discriminant up
// to a constant.  With J = (J2, J4, J6, J10) the invariants of f, and J2 != 0, the condition is
//     I2^2 J4 = I4 J2^2,   I2^3 J6 = I6 J2^3,   I2^5 J10 = I10 J2^5,   I10 != 0
// (if I2 = 0 the equations force I4 = I6 = I10 = 0, which is excluded; otherwise mu = I2/J2 satisfies I4 = mu^2 J4,
// I6 = mu^3 J6, I10 = mu^5 J10, and conversely).  All lines are covered by three charts:
//   (1) x2 = A x0 + B x1: the ideal (E1, E2, E3) of Q[A, B].  If it is the unit ideal there is no solution over Qbar.
//       Otherwise it must have dimension 0, and the solutions with I10 != 0 are the points of the components of its radical
//       decomposition that do not contain I10; their number over Qbar is the sum of the degrees of those components.
//       When the Groebner basis says "unit ideal", this is proved a second time without Groebner bases: a common zero
//       (a, b) of E1, E2, E3 would make a a common root of the resultants Res_B(E1, E2) and Res_B(E1, E3), and the script
//       checks that these two polynomials in A are coprime.
//   (2) a x0 + x1 = 0, the lines through (0:0:1) other than x0 = 0: the gcd of E1, E2, E3 in Q[a], with every factor it has in
//       common with I10 removed; the number of solutions over Qbar is the degree of its squarefree part.
//   (3) the line x0 = 0.
// The invariants are computed over Q(A, B) and Q(a) from the universal formulas, and the script checks that they are
// polynomials, so that they specialise correctly, also where the section has a root at infinity.
SetColumns(0); SetQuitOnError(true);
Q := Rationals(); PX<x> := PolynomialRing(Q);
f := PX![6838225, -46160230, 40133975, 172162750, 60026500, -16144220, 743425];
J := IgusaClebschInvariants(f);
assert J[1] ne 0 and J[4] ne 0;
P3<x0, x1, x2> := PolynomialRing(Q, 3);
mons := [x0^i*x1^j*x2^(6-i-j) : j in [6-i..0 by -1], i in [6..0 by -1]];
R2<A, B> := PolynomialRing(Q, 2); F2 := FieldOfFractions(R2); PU<u> := PolynomialRing(F2);
PA<a> := PolynomialRing(Q); FA := FieldOfFractions(PA); PV<v> := PolynomialRing(FA);
Eqs := func< I | [I[1]^2*J[2] - I[2]*J[1]^2, I[1]^3*J[3] - I[3]*J[1]^3, I[1]^5*J[4] - I[4]*J[1]^5] >;

// returns n1, n2, n3 (numbers of lines over Qbar in the three charts) and the list of rational lines among them
function Lines(S)
  rat := [];
  // chart 1
  I := IgusaClebschInvariants(Evaluate(S, [PU | u, 1, A*u + B]));
  assert forall{c : c in I | Denominator(c) eq 1};
  I := [R2 | Numerator(c) : c in I];
  n1 := 0;
  if I[4] ne 0 then
    Id := ideal< R2 | Eqs(I) >;
    d := Dimension(Id);
    error if d gt 0, "chart 1: the solution set is not finite";
    if d eq -1 then
      E := Eqs(I);
      r12 := UnivariatePolynomial(Resultant(E[1], E[2], B)); r13 := UnivariatePolynomial(Resultant(E[1], E[3], B));
      error if r12 eq 0 or r13 eq 0 or Degree(GCD(r12, r13)) ne 0, "chart 1: the resultants do not confirm that there is no solution";
    end if;
    if d eq 0 then
      for P in RadicalDecomposition(Id) do
        if I[4] in P then continue; end if;
        n1 +:= QuotientDimension(P);
        if QuotientDimension(P) eq 1 then
          pt := Variety(P)[1]; den := LCM(Denominator(pt[1]), Denominator(pt[2]));
          Append(~rat, [Integers() | pt[1]*den, pt[2]*den, -den]);
        end if;
      end for;
    end if;
  end if;
  // chart 2
  I := IgusaClebschInvariants(Evaluate(S, [PV | v, -a*v, 1]));
  assert forall{c : c in I | Denominator(c) eq 1};
  I := [PA | Numerator(c) : c in I];
  n2 := 0;
  if I[4] ne 0 then
    g := GCD(Eqs(I));
    error if g eq 0, "chart 2: the equations vanish identically";
    while Degree(GCD(g, I[4])) gt 0 do g := g div GCD(g, I[4]); end while;
    g := g div GCD(g, Derivative(g));
    n2 := Degree(g);
    for r in Roots(g) do Append(~rat, [Integers() | Numerator(r[1]), Denominator(r[1]), 0]); end for;
  end if;
  // chart 3
  n3 := 0;
  g0 := Evaluate(S, [PX | 0, x, 1]);
  if Degree(g0) ge 5 and IsSquarefree(g0) then
    I := IgusaClebschInvariants(g0);
    if I[4] ne 0 and forall{e : e in Eqs(I) | e eq 0} then n3 := 1; Append(~rat, [1, 0, 0]); end if;
  end if;
  return n1, n2, n3, rat;
end function;

models := [];
for ln in Split(Read("data/models.txt"), "\n") do
  if #ln eq 0 or ln[1] eq "#" then continue; end if;
  F := Split(ln, "|"); c := [StringToInteger(t) : t in Split(F[3], " ") | #t gt 0];
  assert #c eq 28;
  Append(~models, <&cat Split(F[1], " "), &cat Split(F[2], " ") cat " nodes", c>);
end for;
NM := #models; OWN := "t0453";
assert #{m[1] : m in models} eq NM;
s3 := [StringToInteger(t) : t in Split(Read("data/S3_coefficients.txt"), " \n") | #t gt 0 and t[1] ne "#"];
assert [m[3] : m in models | m[1] eq OWN] eq [s3];           // the row t0453 is the sextic S_3 of the note
todo := models;
if not assigned only then only := ""; end if;
if not assigned lo then lo := "1"; end if;
if not assigned hi then hi := IntegerToString(NM); end if;
full := only eq "" and lo eq "1" and hi eq IntegerToString(NM);
if only ne "" then todo := [m : m in todo | m[1] eq only]; error if #todo eq 0, "no such model";
else todo := todo[StringToInteger(lo)..Min(StringToInteger(hi), #todo)]; end if;
bad := 0; nother := 0; nown := 0;
for m in todo do
  t := Cputime();
  S := &+[m[3][k]*mons[k] : k in [1..28]];
  n1, n2, n3, rat := Lines(S);
  printf "%o %o: lines over Qbar with a section isomorphic to C_738: %o + %o + %o in the three charts%o  (%o s)\n", m[1], m[2], n1, n2, n3,
    #rat gt 0 select "; rational:" cat &cat[Sprintf(" (%o, %o, %o)", l[1], l[2], l[3]) : l in rat] else "", Cputime(t);
  if m[1] eq OWN then
    nown +:= 1;
    // its own model: exactly one line, the line (66, 31, -52) of the note, and its section is isomorphic to C_738 over Q
    ok := [n1, n2, n3] eq [1, 0, 0] and #rat eq 1 and rat[1] in {[66, 31, -52], [-66, -31, 52]};
    if ok then
      g := Evaluate(S, [PX | 6*x - 2, 4*x + 16, 10*x + 7]);
      ok := g eq 64*f and IsIsomorphic(HyperellipticCurve(g), HyperellipticCurve(f));
    end if;
    printf "the model S_3: %o\n", ok select "exactly one line over Qbar, the line (66, 31, -52), and S_3 on it is 64 f" else "FAILED";
    if not ok then bad +:= 1; end if;
  else
    nother +:= 1;
    if n1 + n2 + n3 ne 0 then bad +:= 1; end if;
  end if;
end for;
if bad eq 0 and full then
  assert nother eq NM - 1 and nown eq 1;
  printf "LINE SECTION CHECK PASSED: C_738 is the section of S_3 by exactly one line over Qbar, and it is not a line section of any of the other %o sextics\n", NM - 1;
elif bad eq 0 then
  printf "partial run: %o models, no failure\n", #todo;
else
  printf "LINE SECTION CHECK FAILED on %o models\n", bad;
end if;
quit;
