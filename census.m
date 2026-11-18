// census.m -- the models of data/models.txt against the lattice census of the note (data/census.txt).
//   magma -b gramout:=grams.txt models.m      (writes the lattice D-perp of every model; grams.txt must not exist)
//   magma -b grams:=grams.txt census.m
// data/census.txt lists the classes of polarizations of X_163 of degree 2 with one, two or three nodes that the random
// walks of the census met (Section 3 of the note): 7, 186 and 2192 lattices L = D-perp(-1), positive definite of rank 19.
// The script proves:
//  1. every row of data/census.txt is an even positive definite lattice of rank 19 and determinant 326 with exactly as
//     many pairs of roots as its type says (A1, A1^2, A1^3), and no two rows are isometric;
//  2. the lattice D-perp of every model of data/models.txt is isometric to exactly one row, of the same number of nodes,
//     and different models go to different rows.
// So the models constructed are classes of the census: all seven with one node, all 186 with two, and 2188 of the
// 2192 with three.  The four rows that are not models are printed.
// (Two polarizations give isomorphic models exactly when their lattices D-perp are isometric; see models.m.)
SetColumns(0); SetQuitOnError(true);
if not assigned grams then error "usage: magma -b grams:=FILE census.m, with FILE written by magma -b gramout:=FILE models.m"; end if;
nodesof := AssociativeArray(); nodesof["A1"] := 1; nodesof["A1^2"] := 2; nodesof["A1^3"] := 3;
cen := [* *];
for ln in Split(Read("data/census.txt"), "\n") do
  if #ln eq 0 or ln[1] eq "#" then continue; end if;
  F := [&cat Split(t, " ") : t in Split(ln, "|")[1..3]] cat [Split(ln, "|")[4]];
  G := Matrix(Integers(), 19, 19, [StringToInteger(t) : t in Split(F[4], " ,[]") | #t gt 0]);
  assert IsSymmetric(G) and IsPositiveDefinite(G) and Determinant(G) eq 326 and forall{i : i in [1..19] | IsEven(G[i,i])};
  L := LatticeWithGram(G); key := [#ShortVectors(L, m, m) : m in [2, 4, 6]];
  assert key[1] eq nodesof[F[2]];                  // pairs of roots = contracted curves, for A1^k
  Append(~cen, <F[1], F[2], L, key, StringToInteger(F[3])>);
end for;
assert #{c[1] : c in cen} eq #cen;
buck := AssociativeArray();
for i in [1..#cen] do
  k := cen[i][4];
  if IsDefined(buck, k) then
    assert forall{j : j in buck[k] | not IsIsometric(cen[j][3], cen[i][3])};        // 1: pairwise non-isometric
    Append(~buck[k], i);
  else buck[k] := [i]; end if;
end for;
printf "%o classes of the census: %o, %o, %o with one, two, three nodes; pairwise non-isometric\n", #cen, #[c : c in cen | c[2] eq "A1"], #[c : c in cen | c[2] eq "A1^2"], #[c : c in cen | c[2] eq "A1^3"];
mrows := [r : r in Split(Read("data/models.txt"), "\n") | #r gt 0 and r[1] ne "#"];
names := [&cat Split(Split(r, "|")[1], " ") : r in mrows]; forced := [2*#Split(Split(r, "|")[5], ";") : r in mrows];
rows := [r : r in Split(Read(grams), "\n") | #r gt 0];
assert [Split(r, "|")[1] : r in rows] eq names;      // one lattice for every model, in the order of models.txt
used := AssociativeArray();
for n in [1..#rows] do
  r := rows[n]; F := Split(r, "|"); k := StringToInteger(F[2]); key := eval F[4]; L := LatticeWithGram(Matrix(Integers(), 19, 19, eval F[5]));
  error if not IsDefined(buck, key), F[1], "is not a class of the census";
  hits := [i : i in buck[key] | IsIsometric(cen[i][3], L)];
  error if #hits ne 1, F[1], "is not a class of the census";
  i := hits[1];
  error if nodesof[cen[i][2]] ne k or IsDefined(used, i), F[1], ": wrong number of nodes, or a class met twice";
  error if forced[n] ne cen[i][5], F[1], ": the number of forced points differs from the census";
  used[i] := F[1]; printf "%o %o\n", F[1], cen[i][1];
end for;
for t in ["A1", "A1^2", "A1^3"] do
  ids := [i : i in [1..#cen] | cen[i][2] eq t];
  printf "%o: %o classes in the census, %o of them are models of the list%o\n", t, #ids, #[i : i in ids | IsDefined(used, i)],
    #[i : i in ids | not IsDefined(used, i)] eq 0 select "" else Sprintf("; not: %o", [cen[i][1] : i in ids | not IsDefined(used, i)]);
end for;
printf "CENSUS CHECK PASSED: the %o models are %o classes of the census\n", #rows, #Keys(used);
quit;
