/* minred.gp : applies the change of coordinates found by minred.m (Magma MinRedTernaryForm) to the reduced model
   out/model_<NAME>_red.gp and rewrites that file and out/model_<NAME>_red_coeffs.txt, so that sing.m, sqlines.sh and
   finalize.gp verify the minimised model.  The model before this step is kept as out/model_<NAME>_red0.gp.
   The sextic is transformed by N and normalised up to squares (normvec), which keeps the quadratic twist of the double
   plane; it must be proportional to the form returned by MinRedTernaryForm, and every known line must still be a
   square line.  usage: ./run.sh minred.gp logs/minred_<NAME>.log 'NAME="<NAME>"' */
\r cplib.gp
\r reducelib.gp
if(type(NAME) != "t_STR", error("set NAME"));
read(concat(["out/model_", NAME, "_red.gp"]));
read(concat(["out/model_", NAME, "_minred.gp"]));
{
my(c0 = c, S0 = S, M0 = M, c1, S1, k);
c1 = normvec(transformvec(S0, MRN)); S1 = polfromvec(c1);
k = 0; for(n = 1, 28, if(MRG[n] != 0, k = c1[n]/MRG[n]; break));
if(k == 0 || c1 != k*MRG, error("minred.gp: the transformed sextic is not proportional to the form of MinRedTernaryForm"));
if(!issquare(MRC*k), error("minred.gp: the normalised sextic is not S(N^-1 y) up to a square"));
writefile(concat(["out/model_", NAME, "_red0.gp"]), Str("\\\\ the reduced model before minred.gp\nc = ", c0, ";\nM = ", M0, ";"));
c = c1; S = S1; M = MRN*M0;
trit = translines(trit, MRN); nodelines = translines(nodelines, MRN);
sing = vector(#sing, n, prim(MRN*sing[n]~)~);
if(nodeimg != 0, nodeimg = prim(MRN*nodeimg~)~);
if(vector(#trit, n, issqline(S, trit[n])) != vector(#trit, n, 1) || vector(#nodelines, n, issqline(S, nodelines[n])) != vector(#nodelines, n, 1), error("minred.gp: a known line is not a square line of the minimised model"));
print("minimised model: max |coeff| ", vecmax(abs(c0)), " -> ", vecmax(abs(c)), "; det N = ", matdet(MRN), "; twist factor ", k, " (MinRedTernaryForm form times this); singular points ", sing);
print("max |coefficient| of the known lines: ", vecmax(vector(#trit + #nodelines, n, vecmax(abs(concat(trit, nodelines)[n])))));
}
writefile(concat(["out/model_", NAME, "_red.gp"]), Str("\\\\ reduced and minimised model (models/reduce.gp, minred.gp) of ", NAME, ": y = M x, S(y) = S_raw(M^-1 y) up to a square\nNAME = \"", NAME, "\";\nPAIR = ", PAIR, ";\nc = ", c, ";\nS = ", S, ";\nM = ", M, ";\nsing = ", sing, ";\nnodeimg = ", nodeimg, ";\ntrit = ", trit, ";\nnodelines = ", nodelines, ";\n"));
writefile(concat(["out/model_", NAME, "_red_coeffs.txt"]), Str(c));
print("written out/model_", NAME, "_red.gp (minimised)");
\q
