// minred.m -- minimisation and reduction of a plane sextic by Magma's MinRedTernaryForm (Elsenhans-Stoll), for the model
// pipeline (pipeline.sh, between reduce.gp and sing.m) and for existing model files (minred_file.sh).
//
//   magma -b file:=COEFFS.txt out:=OUT.gp minred.m
//
// COEFFS.txt holds the 28 coefficients of S in the order of src/linesearch.c (for i = 6..0, j = 6-i..0: the coefficient
// of x0^i x1^j x2^(6-i-j); brackets and commas are tolerated).  OUT.gp receives, as gp assignments,
//   MRN  = the 3x3 rational matrix N of the change of coordinates y = N x (the convention of reduce.gp: the new sextic
//          is S(N^-1 y) up to a square factor, points go to N p, lines to l N^-1),
//   MRG  = the 28 coefficients of the form g returned by MinRedTernaryForm, and
//   MRC  = the constant with S(N^-1 y) = MRC g(y).
// The double plane w^2 = S is only determined up to squares, so the caller must NOT replace S by g when MRC is not a
// square: it transforms S by N and removes square factors (normvec in reducelib.gp), which keeps the quadratic twist.
// The identity S(N^-1 y) = MRC g(y) is checked here before anything is written.
SetColumns(0); SetQuitOnError(true);
Q := Rationals(); P<x0,x1,x2> := PolynomialRing(Q, 3);
mons := []; for d in [0..6] do for j in [0..d] do Append(~mons, x0^(6-d)*x1^(d-j)*x2^j); end for; end for;
c := [StringToInteger(t) : t in Split(Read(file), " \n,[]") | #t gt 0 and t[1] ne "#"];
error if #c ne 28, "minred.m: expected 28 coefficients in " cat file;
f := &+[c[i]*mons[i] : i in [1..28]];
g, T := MinRedTernaryForm(f);
// MinRedTernaryForm does not document which of T, T^-1 and their transposes realises g: find the matrix A with
// f(x A) = const * g(x) (x a row vector) by trying them; then x_old = x_new A, so N = (A^T)^-1 in the convention above.
found := false;
for A in [T, T^-1, Transpose(T), Transpose(T)^-1] do
  h := Evaluate(f, [&+[A[i][j]*P.i : i in [1..3]] : j in [1..3]]);
  if IsDivisibleBy(h, g) and TotalDegree(h div g) eq 0 then N := Transpose(Matrix(Q, A))^-1; cst := LeadingCoefficient(h div g); found := true; break; end if;
end for;
error if not found, "minred.m: could not match the form returned by MinRedTernaryForm with a transform of the input";
Ni := N^-1;
assert Evaluate(f, [&+[Ni[j][i]*P.i : i in [1..3]] : j in [1..3]]) eq cst*g;      // S(N^-1 y) = cst g(y)
gc := [MonomialCoefficient(g, m) : m in mons]; assert &+[gc[i]*mons[i] : i in [1..28]] eq g;
row := func< v | "[" cat &cat[Sprint(v[i]) cat (i lt #v select ", " else "") : i in [1..#v]] cat "]" >;
PrintFile(out, Sprintf("\\\\ minred.m: MinRedTernaryForm on %o\nMRN = [%o; %o; %o];\nMRG = %o;\nMRC = %o;", file,
  row(Eltseq(N[1]))[2..#row(Eltseq(N[1]))-1], row(Eltseq(N[2]))[2..#row(Eltseq(N[2]))-1], row(Eltseq(N[3]))[2..#row(Eltseq(N[3]))-1], row(gc), cst) : Overwrite := true);
printf "minred.m: max |coeff| %o -> %o, det N = %o, constant %o\nwritten %o\n", Max([Abs(t) : t in c]), Max([Abs(t) : t in gc]), Determinant(N), cst, out;
quit;
