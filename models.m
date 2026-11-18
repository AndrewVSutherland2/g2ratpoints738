// models.m -- every sextic of data/models.txt is a double-plane model of the K3 surface X_163.
//   magma -b models.m                    (all rows and the classification)
//   magma -b lo:=1 hi:=100 models.m      (a range of rows; the classification is skipped)
// For each row (name | number of nodes | 28 coefficients of S | singular points | lines) it proves, from the row alone:
//  1. the singular subscheme of the sextic S = 0 has degree k (the number of nodes stated) and consists of the k
//     rational points listed.  A singular point at which that subscheme is reduced is an ordinary node, so the
//     double plane w^2 = S has k ordinary double points and its minimal resolution X is a K3 surface with k
//     contracted curves E_1, .., E_k.  (A plane sextic with at most three nodes and no other singular point is
//     irreducible: two components of degrees a and 6 - a meet in a(6 - a) >= 5 points, counted with multiplicity.)
//  2. S restricts to the square of a polynomial with rational coefficients on each line listed.  So each line has
//     two lifts L+ and L- to X, smooth rational curves defined over Q.  Their intersection numbers: L.L = E.E = -2;
//     L.E = 1 if the line passes through the node (the pullback D of the line contains E once, D.E = 0, and the
//     covering involution exchanges L+ and L-); L+.L- = 3 - (number of nodes on the line), because D.L+ = 1; two
//     lines through a common node have disjoint lifts; two lines that meet at a point v off the sextic have
//     L_i^a.L_j^b = 1 exactly when the two lifts take the same value w at v.  (Two of these lines cannot meet at a
//     smooth point of the sextic, where both would be the tangent.)  The count of the lines through a node that
//     touch a branch of the sextic there is printed; the intersection numbers are the same for them.
//  3. the lifts and the E_t span a lattice of rank 20 and discriminant -163.  As the Picard number of a K3 surface
//     in characteristic 0 is at most 20 and 163 is squarefree, this lattice is NS(X), every class is defined over Q,
//     and the transcendental lattice is [2,1;1,82]; so X is X_163 (Shioda-Inose), over Q by the descent argument
//     of the note.  The pullback D of a line has D^2 = 2.
//  4. in that lattice the classes x with x^2 = -2, x.D = 0 are exactly the +-E_t, and the irreducible classes with
//     x^2 = -2, x.D = 1 are exactly the lifts: the list of lines of the row is complete.
// After the last row it sorts the models into isomorphism classes.  Two models are isomorphic (as double planes of
// X_163, over Q) exactly when their lattices D-perp are isometric: the discriminant group of NS(X) is cyclic of
// order 163, so the orbit of D under O(NS(X)) is determined by the isometry class of D-perp, and an isometry that
// takes one nef class to another can be changed by reflections into one that preserves the nef cone.
// Output: one line per row (name, nodes, lines through 0 / 1 / 2 / 3 nodes, forced points), then the classes.
SetColumns(0); SetQuitOnError(true);
if not assigned file then file := "data/models.txt"; end if;
rows := [r : r in Split(Read(file), "\n") | #r gt 0 and r[1] ne "#"];
ints := func< s | [StringToInteger(t) : t in Split(s, " ,") | #t gt 0] >;
full := not assigned lo and not assigned hi;
a := assigned lo select StringToInteger(lo) else 1; b := assigned hi select StringToInteger(hi) else #rows;
Q := Rationals(); P<z0,z1,z2> := PolynomialRing(Q, 3); R<x> := PolynomialRing(Q); P2 := ProjectiveSpace(P);
mons := []; for d in [0..6] do for j in [0..d] do Append(~mons, z0^(6-d)*z1^(d-j)*z2^j); end for; end for;
prim := function(v) v := [Q | t : t in v]; m := LCM([Denominator(t) : t in v]); v := [Integers() | t*m : t in v]; g := GCD(v); v := [t div g : t in v];
  if [t : t in v | t ne 0][1] lt 0 then v := [-t : t in v]; end if; return v; end function;
names := []; grams := []; keys := [* *];
for r in [a..b] do
  F := Split(rows[r], "|"); name := &cat Split(F[1], " "); k := StringToInteger(&cat Split(F[2], " "));
  c := ints(F[3]); assert #c eq 28; S := &+[c[i]*mons[i] : i in [1..28]];
  nodes := [prim(ints(t)) : t in Split(F[4], ";")]; lines := [ints(t) : t in Split(F[5], ";")]; n := #lines;
  assert #nodes eq k and #Set(nodes) eq k and #{prim(l) : l in lines} eq n;
  // 1. singular points
  Sing := SingularSubscheme(Curve(P2, S)); assert Dimension(Sing) eq 0 and Degree(Sing) eq k;
  assert {prim(Eltseq(p)) : p in RationalPoints(Sing)} eq Set(nodes);
  // 2. the lines
  LD := [];
  for l in lines do
    K := Basis(Kernel(Matrix(Q, 3, 1, l))); p1 := Eltseq(K[1]); p2 := Eltseq(K[2]);
    ok, g := IsSquare(Evaluate(S, [p1[i]*x + p2[i] : i in [1..3]])); assert ok;
    Append(~LD, <p1, p2, g>);
  end for;
  wplus := function(ld, v)
    cc := Solution(Matrix(Q, 2, 3, [ld[1], ld[2]]), Vector(Q, v));
    return &+[Coefficient(ld[3], e)*cc[1]^e*cc[2]^(3-e) : e in [0..3]];
  end function;
  onnode := [[&+[lines[i][e]*nd[e] : e in [1..3]] eq 0 : nd in nodes] : i in [1..n]];
  nn := [#[t : t in onnode[i] | t] : i in [1..n]];
  ntan := 0;       // lines through a node that touch a branch of the sextic there
  for i in [1..n] do for t in [1..k] do if onnode[i][t] then
    g := LD[i][3]; cc := Solution(Matrix(Q, 2, 3, [LD[i][1], LD[i][2]]), Vector(Q, nodes[t]));
    if cc[2] ne 0 then assert Evaluate(g, cc[1]/cc[2]) eq 0; if Evaluate(Derivative(g), cc[1]/cc[2]) eq 0 then ntan +:= 1; end if;
    else assert Degree(g) le 2; if Degree(g) lt 2 then ntan +:= 1; end if; end if;
  end if; end for; end for;
  // 3. the intersection matrix.  Index 2i-1 = L_i+, 2i = L_i-, 2n+t = E_t.
  N := 2*n + k; M := ZeroMatrix(Integers(), N, N);
  for i in [1..n] do
    M[2*i-1,2*i-1] := -2; M[2*i,2*i] := -2; M[2*i-1,2*i] := 3 - nn[i]; M[2*i,2*i-1] := 3 - nn[i];
    for t in [1..k] do if onnode[i][t] then
      for e in [0,1] do M[2*i-e,2*n+t] := 1; M[2*n+t,2*i-e] := 1; end for;
    end if; end for;
  end for;
  for t in [1..k] do M[2*n+t,2*n+t] := -2; end for;
  for i in [1..n] do for j in [i+1..n] do
    v := Eltseq(Basis(Kernel(Transpose(Matrix(Q, 2, 3, [lines[i], lines[j]]))))[1]);
    if exists{nd : nd in nodes | Rank(Matrix(Q, 2, 3, [[Q | t : t in nd], v])) eq 1} then continue; end if;
    sv := Evaluate(S, v); assert sv ne 0;
    wi := wplus(LD[i], v); wj := wplus(LD[j], v); assert wi^2 eq sv and wj^2 eq sv;
    for e in [0,1] do aa := 2*i - e; bb := wi eq wj select 2*j - e else 2*j - 1 + e; M[aa,bb] := 1; M[bb,aa] := 1; end for;
  end for; end for;
  assert Rank(M) eq 20;
  Ker := KernelMatrix(M); Sm, U, V := SmithForm(Ker); assert Set(Diagonal(Sm)) eq {1};
  T := V^-1; Bas := RowSubmatrix(T, Nrows(Ker) + 1, 20); G := Bas*M*Transpose(Bas);
  assert Determinant(G) eq -163;
  coord := func< cc | Vector(Integers(), [(cc*V)[e] : e in [Nrows(Ker)+1..N]]) >;
  ev := func< e | Vector(Integers(), [i eq e select 1 else 0 : i in [1..N]]) >;
  ip := func< u, v | (Matrix(u)*G*Transpose(Matrix(v)))[1,1] >;
  pull := func< i | coord(ev(2*i-1) + ev(2*i) + &+[Parent(ev(1)) | ev(2*n+t) : t in [1..k] | onnode[i][t]]) >;
  Dcl := pull(1); assert forall{i : i in [1..n] | pull(i) eq Dcl} and ip(Dcl, Dcl) eq 2;
  // 4. roots and lines in the lattice.  Pd(x) = (x.D)^2 - x^2 is positive definite on NS(X).
  dv := Matrix(Dcl)*G; Pd := Transpose(dv)*dv - G; assert IsPositiveDefinite(Pd) and Determinant(Pd) eq 163;
  LP := LatticeWithGram(Pd);
  sh := [Vector(Integers(), Eltseq(Coordinates(LP, t[1]))) : t in ShortVectors(LP, 3)]; sh := sh cat [-v : v in sh];
  Ecl := [coord(ev(2*n+t)) : t in [1..k]]; Lcl := [coord(ev(e)) : e in [1..2*n]]; assert #Set(Lcl) eq 2*n;
  assert {v : v in sh | ip(v, Dcl) eq 0 and ip(v, v) eq -2} eq Set(Ecl) join {-v : v in Ecl};
  assert {v : v in sh | ip(v, Dcl) eq 1 and ip(v, v) eq -2 and forall{E : E in Ecl | ip(v, E) ge 0}} eq Set(Lcl);
  perp := Matrix(Basis(Kernel(Transpose(dv)))); Lperp := -perp*G*Transpose(perp);
  assert Determinant(Lperp) eq 326 and IsPositiveDefinite(Lperp);
  cnt := [#[i : i in [1..n] | nn[i] eq e] : e in [0..3]]; assert &+cnt eq n;
  printf "%o %o: %o lines (%o, %o, %o, %o through 0, 1, 2, 3 nodes%o), %o forced points; NS of rank 20, discriminant -163\n", name, k, n, cnt[1], cnt[2], cnt[3], cnt[4],
    ntan eq 0 select "" else Sprintf("; %o touching a branch at a node", ntan), 2*n;
  if full or assigned gramout then
    LL := LatticeWithGram(Lperp);
    key := <k, cnt, [#ShortVectors(LL, m, m) : m in [2, 4, 6]]>;
    Append(~names, name); Append(~grams, Lperp); Append(~keys, key);
    if assigned gramout then PrintFile(gramout, Sprintf("%o|%o|%o|%o|%o", name, k, cnt, key[3], Eltseq(Lperp))); end if;
  end if;
end for;
printf "ROWS %o..%o PASSED: each is a model of X_163 with exactly the singular points and lines listed\n", a, b;
if full then
  // isomorphism classes: models with equal invariants are compared with IsIsometric
  cls := []; rep := [];             // cls[i] = index of the class of model i; rep = first member of each class
  for i in [1..#names] do
    found := 0;
    for j in [1..#rep] do
      if keys[rep[j]] eq keys[i] and IsIsometric(LatticeWithGram(grams[rep[j]]), LatticeWithGram(grams[i])) then found := j; break; end if;
    end for;
    if found eq 0 then Append(~rep, i); found := #rep; end if;
    Append(~cls, found);
  end for;
  for kk in [1..3] do
    members := [i : i in [1..#names] | keys[i][1] eq kk]; cl := {cls[i] : i in members};
    printf "%o nodes: %o sextics in %o isomorphism classes\n", kk, #members, #cl;
    for j in cl do mm := [names[i] : i in members | cls[i] eq j]; if #mm gt 1 then printf "  the same model: %o\n", mm; end if; end for;
  end for;
  printf "MODELS CHECK PASSED: %o sextics, %o isomorphism classes\n", #names, #rep;
end if;
quit;
