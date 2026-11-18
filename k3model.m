// k3model.m -- a Magma certificate for the double plane w^2 = S_3(z0,z1,z2) that carries the curve C_738.
//   magma -b k3model.m                                  (a few seconds)
//   magma -b conics:=data/S3_conics_810.txt k3model.m   (checks the table of conics as well)
// It proves, from the sextic S_3 and the list of 39 lines below and nothing else:
//  1. S_3 is absolutely irreducible and its only singular points are three rational ordinary nodes, so the minimal
//     resolution X of the double plane is a K3 surface with three contracted curves E_1, E_2, E_3.
//  2. S_3 restricts to a perfect square on each of the 39 lines: 21 through no node, 15 through one, 3 through two.
//     Each such line has two lifts L+ and L- to X, smooth rational curves defined over Q.
//  3. The 78 lifts and the E_n span a lattice of rank 20 and discriminant -163.  The Picard number of a K3 surface in
//     characteristic 0 is at most 20 and 163 is squarefree, so this lattice is all of NS(X), every class is defined
//     over Q, and the transcendental lattice is the unique even positive definite binary form of determinant 163,
//     [2,1;1,82].  A K3 surface of Picard number 20 is determined by its transcendental lattice (Shioda-Inose), so
//     X is X_163.
//  4. From the lattice alone: the classes x with x^2 = -2, x.D = 0 form the root system A_1^3; the irreducible
//     classes with x^2 = -2, x.D = 1 are exactly the 78 lifts, so S_3 is a square on exactly 39 lines (78 forced points);
//     the irreducible classes with x^2 = -2, x.D = 2 number 1620, so there are exactly 810 conics on which S_3 is a
//     square.  (An irreducible class of the last kind is a lift of such a conic, or the pullback of a line through
//     two nodes on which S_3 is not a square; by step 2 all three lines through two nodes are among the 39.)
//     (D is the pullback of a line.  A class x with x^2 = -2 and x.D > 0 is effective, and it is
//     irreducible if and only if x.C >= 0 for every smooth rational curve C with C.D < x.D.)
//  5. #X(F_p) = p^2 + 20p + 1 + a_p for the good primes p < 200, where a_p is the coefficient of the CM newform of
//     weight 3 and level 163: a_p = (u^2 - 163 v^2)/2 if 4p = u^2 + 163 v^2, and a_p = 0 if p is inert in
//     Q(sqrt(-163)).
//  6. The line 66 z0 + 31 z1 - 52 z2 = 0, parametrised by (6x - 2 : 4x + 16 : 10x + 7), cuts out y^2 = f_738(x);
//     the 39 lines meet it at 39 distinct values of x, each giving two rational points of C_738.
SetColumns(0); SetQuitOnError(true);
// coefficients of S_3: for d = 0..6 and j = 0..d, the coefficient of z0^(6-d) z1^(d-j) z2^j
S3coef := [25, -1930, 3240, 969, -1734, 1996, 5752, -6326, -996, -864, 2748, -1328, -3588, 960, -176,
           -144, 216, -548, 1376, 16, 0, 36, 12, -311, -52, 676, 0, 0];
// the 39 lines: [a, b, c] is a z0 + b z1 + c z2 = 0
lines := [[0,2,-3], [1,-2,-2], [1,3,-2], [1,42,-62], [5,-11,-4], [5,0,-2], [5,1,4], [5,4,4],
          [6,-4,-3], [10,-3,-8], [10,0,3], [11,1,2], [12,-2,-5], [17,6,10], [19,12,4], [20,26,5],
          [25,-21,4], [25,-3,10], [35,-60,6], [40,30,31], [59,56,-60], [0,1,-2], [1,-2,0], [1,0,0],
          [1,1,-2], [1,1,0], [1,3,-4], [2,3,0], [3,-2,0], [3,-1,2], [4,3,2], [4,5,-2],
          [5,-7,0], [5,1,-2], [6,-11,0], [7,4,6], [8,-1,2], [13,8,10], [23,-2,4]];
f738coef := [6838225, -46160230, 40133975, 172162750, 60026500, -16144220, 743425];       // constant term first
Q := Rationals(); P<z0,z1,z2> := PolynomialRing(Q, 3); R<x> := PolynomialRing(Q);
mons := []; for d in [0..6] do for j in [0..d] do Append(~mons, z0^(6-d)*z1^(d-j)*z2^j); end for; end for;
S3 := &+[S3coef[i]*mons[i] : i in [1..28]];
n := #lines; assert n eq 39;

// 1. singular points
P2 := ProjectiveSpace(P); B := Curve(P2, S3);
assert IsAbsolutelyIrreducible(B);
Sing := SingularSubscheme(B); assert Dimension(Sing) eq 0 and Degree(Sing) eq 3;
nodes := [Eltseq(p) : p in RationalPoints(Sing)]; assert #nodes eq 3;
printf "ok 1: S_3 is absolutely irreducible; its singular subscheme has degree 3 and consists of the\n";
printf "      rational points %o, %o, %o: three ordinary nodes\n", nodes[1], nodes[2], nodes[3];

// 2. the 39 lines.  For a line l with basis p1, p2 of its points, S_3(s p1 + p2) = g(s)^2; the lift L+ is
//    w = G(a, b) at the point a p1 + b p2, where G is g made homogeneous of degree 3, and L- is w = -G.
LD := [];
for l in lines do
  K := Basis(Kernel(Matrix(Q, 3, 1, l))); p1 := Eltseq(K[1]); p2 := Eltseq(K[2]);
  ok, g := IsSquare(Evaluate(S3, [p1[i]*x + p2[i] : i in [1..3]])); assert ok;
  Append(~LD, <p1, p2, g>);
end for;
function wplus(ld, v)
  c := Solution(Matrix(Q, 2, 3, [ld[1], ld[2]]), Vector(Q, v));
  return &+[Coefficient(ld[3], k)*c[1]^k*c[2]^(3-k) : k in [0..3]];
end function;
onnode := [[&+[lines[i][k]*nd[k] : k in [1..3]] eq 0 : nd in nodes] : i in [1..n]];
nn := [#[b : b in onnode[i] | b] : i in [1..n]];
for i in [1..n] do for t in [1..3] do if onnode[i][t] then      // the line meets the sextic exactly twice at the node
  g := LD[i][3]; c := Solution(Matrix(Q, 2, 3, [LD[i][1], LD[i][2]]), Vector(Q, nodes[t]));
  if c[2] ne 0 then assert Evaluate(g, c[1]/c[2]) eq 0 and Evaluate(Derivative(g), c[1]/c[2]) ne 0;
  else assert Degree(g) eq 2; end if;
end if; end for; end for;
assert [#[i : i in [1..n] | nn[i] eq k] : k in [0..3]] eq [21, 15, 3, 0];
printf "ok 2: S_3 is a square on all 39 lines; 21 pass through no node, 15 through one, 3 through two\n";

// 3. intersection numbers on X.  Index 2i-1 = L_i+, 2i = L_i-, 2n+t = E_t.
//    L.L = E.E = -2;  L.E = 1 if the line passes through the node;  L_i+.L_i- = 3 - (number of nodes on the line),
//    because L_i+ + L_i- + (the E_t on the line) = D and D.L = 1;  two lines through a common node have
//    disjoint lifts;  two lines meeting at a point v off the sextic have L_i^a.L_j^b = 1 exactly when the two
//    lifts take the same value w at v.
N := 2*n + 3; M := ZeroMatrix(Integers(), N, N);
for i in [1..n] do
  M[2*i-1,2*i-1] := -2; M[2*i,2*i] := -2; M[2*i-1,2*i] := 3 - nn[i]; M[2*i,2*i-1] := 3 - nn[i];
  for t in [1..3] do if onnode[i][t] then
    for e in [0,1] do M[2*i-e,2*n+t] := 1; M[2*n+t,2*i-e] := 1; end for;
  end if; end for;
end for;
for t in [1..3] do M[2*n+t,2*n+t] := -2; end for;
for i in [1..n] do for j in [i+1..n] do
  v := Eltseq(Basis(Kernel(Transpose(Matrix(Q, 2, 3, [lines[i], lines[j]]))))[1]);
  if exists{nd : nd in nodes | Rank(Matrix(Q, 2, 3, [nd, v])) eq 1} then continue; end if;
  sv := Evaluate(S3, v); assert sv ne 0;
  wi := wplus(LD[i], v); wj := wplus(LD[j], v); assert wi^2 eq sv and wj^2 eq sv;
  for e in [0,1] do a := 2*i - e; b := wi eq wj select 2*j - e else 2*j - 1 + e; M[a,b] := 1; M[b,a] := 1; end for;
end for; end for;
assert Rank(M) eq 20;
Ker := KernelMatrix(M); S, U, V := SmithForm(Ker); assert Set(Diagonal(S)) eq {1};
// 20 rows completing a basis of the kernel to a basis of Z^81; G is the Gram matrix of the lattice of the 81 curves
T := V^-1; Bas := RowSubmatrix(T, Nrows(Ker) + 1, 20);
G := Bas*M*Transpose(Bas);
assert Determinant(G) eq -163;
// the class of an integer combination c of the 81 curves, in the basis Bas
coord := func< c | Vector(Integers(), [(c*V)[k] : k in [Nrows(Ker)+1..N]]) >;
e := func< k | Vector(Integers(), [i eq k select 1 else 0 : i in [1..N]]) >;
ip := func< a, b | (Matrix(a)*G*Transpose(Matrix(b)))[1,1] >;
pull := func< i | coord(e(2*i-1) + e(2*i) + &+[Parent(e(1)) | e(2*n+t) : t in [1..3] | onnode[i][t]]) >;
Dcl := pull(1); assert forall{i : i in [1..n] | pull(i) eq Dcl};       // every line pulls back to the same class D
assert ip(Dcl, Dcl) eq 2;
printf "ok 3: the 78 lifts and the 3 contracted curves span a lattice of rank 20 and discriminant %o,\n", Determinant(G);
printf "      which is therefore NS(X); D^2 = %o\n", ip(Dcl, Dcl);

// 4. counting curves in the lattice.  Pd(x) = (x.D)^2 - x^2 is a positive definite integral form on NS(X).
dv := Matrix(Dcl)*G; Pd := Transpose(dv)*dv - G; assert IsPositiveDefinite(Pd) and Determinant(Pd) eq 163;
LP := LatticeWithGram(Pd);
sv := [Vector(Integers(), Eltseq(Coordinates(LP, t[1]))) : t in ShortVectors(LP, 6)]; sv := sv cat [-v : v in sv];
deg := func< v | ip(v, Dcl) >;
roots := [v : v in sv | deg(v) eq 0 and ip(v, v) eq -2];
Ecl := [coord(e(2*n+t)) : t in [1..3]];
assert Set(roots) eq Set(Ecl) join {-v : v in Ecl};
Lcl := [coord(e(k)) : k in [1..2*n]]; assert #Set(Lcl) eq 78;
d1 := [v : v in sv | deg(v) eq 1 and ip(v, v) eq -2];
d1irr := [v : v in d1 | forall{E : E in Ecl | ip(v, E) ge 0}];
assert Set(d1irr) eq Set(Lcl);
d2 := [v : v in sv | deg(v) eq 2 and ip(v, v) eq -2];
d2irr := [v : v in d2 | forall{C : C in Ecl cat Lcl | ip(v, C) ge 0}];
assert #d2irr eq 1620;
perp := Matrix(Basis(Kernel(Transpose(dv))));  Lperp := perp*G*Transpose(perp);
assert Determinant(Lperp) eq -326 and IsPositiveDefinite(-Lperp);
printf "ok 4: D-perp is negative definite of rank 19 and determinant 326, with roots +-E_1, +-E_2, +-E_3 (A_1^3);\n";
printf "      of the %o classes with x^2 = -2, x.D = 1 the irreducible ones are exactly the 78 lifts:\n", #d1;
printf "      exactly 39 lines, 78 forced points;\n";
printf "      of the %o classes with x^2 = -2, x.D = 2, %o are irreducible: exactly 810 conics\n", #d2, #d2irr;

// 5. point counts.  Each node is replaced by a conic, so #X(F_p) = sum over the plane of (1 + chi(S_3)) + 3p.
function ap(p)
  if KroneckerSymbol(-163, p) ne 1 then return 0; end if;
  ok, u, v := NormEquation(163, 4*p); assert ok and u^2 + 163*v^2 eq 4*p;
  return (u^2 - 163*v^2) div 2;
end function;
good := []; 
for p in PrimesInInterval(5, 200) do
  F := GF(p); Pp<a0,a1,a2> := PolynomialRing(F, 3); Sp := Pp ! S3;
  // skip the primes at which the sextic has more than its three nodes (bad reduction of the model)
  if Degree(SingularSubscheme(Curve(ProjectiveSpace(Pp), Sp))) ne 3 then continue; end if;
  cnt := 0;
  for pt in [[F|1,a,b] : a, b in F] cat [[F|0,1,b] : b in F] cat [[F|0,0,1]] do
    s := Evaluate(Sp, pt); cnt +:= s eq 0 select 1 else (IsSquare(s) select 2 else 0);
  end for;
  assert cnt + 3*p eq p^2 + 20*p + 1 + ap(p);
  Append(~good, p);
end for;
printf "ok 5: #X(F_p) = p^2 + 20p + 1 + a_p for the %o primes %o\n", #good, good;

// 6. the line of C_738
f738 := R ! f738coef; par := [6*x - 2, 4*x + 16, 10*x + 7];
assert 66*par[1] + 31*par[2] - 52*par[3] eq 0 and Evaluate(S3, par) eq 64*f738;
onl := [l[1]*par[1] + l[2]*par[2] + l[3]*par[3] : l in lines]; assert forall{h : h in onl | Degree(h) eq 1};
xs := [Roots(h)[1][1] : h in onl];
assert #Set(xs) eq 39 and forall{t : t in xs | Evaluate(f738, t) ne 0 and IsSquare(Evaluate(f738, t))};
printf "ok 6: S_3(6x - 2, 4x + 16, 10x + 7) = 64 f_738(x); the 39 lines meet this line at 39 distinct x,\n";
printf "      each with f_738(x) a nonzero square\n";

// 7. (optional) the table of conics: one conic per row, coefficients of z0^2, z0 z1, z0 z2, z1^2, z1 z2, z2^2
if assigned conics then
  rows := [[StringToInteger(t) : t in Split(r, " ")] : r in Split(Read(conics), "\n") | #r gt 0 and r[1] ne "#"];
  qm := [z0^2, z0*z1, z0*z2, z1^2, z1*z2, z2^2];
  qs := {}; U2<uu,vv> := PolynomialRing(Q, 2);
  for r in rows do
    q := &+[r[k]*qm[k] : k in [1..6]]; C := Conic(P2, q); assert IsNonsingular(C);
    phi := Parametrization(C); pe := DefiningPolynomials(phi);
    h := Evaluate(Evaluate(S3, [Evaluate(pp, [x, 1]) : pp in pe]), x);
    assert IsSquare(R ! h) and (Degree(R ! h) mod 2 eq 0);
    Include(~qs, q / LeadingCoefficient(q));
  end for;
  assert #rows eq 810 and #qs eq 810;
  printf "ok 7: S_3 is a square on each of the 810 distinct conics of %o, which by step 4 are all of them\n", conics;
end if;
print "ALL CHECKS PASSED";
quit;
