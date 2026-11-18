load "pencil.m";
// Z_t = conj(F_b) + E : pencil |F| inside H^0(32 D)
Zt := [<"L",15,-1,3>, <"L",17,-1,5>, <"L",29,-1,3>, <"L",32,-1,2>, <"L",33,-1,6>, <"L",37,-1,4>, <"L",39,-1,4>, <"L",52,1,1>, <"C",199,1,2>];
t0 := Cputime(); rows, MA, MB := PencilRows(32, Zt, 1);
M := Transpose(Matrix(Q, rows)); print Nrows(M), Ncols(M), Cputime(t0);
Mp := ChangeRing(M*LCM([Denominator(c) : c in Eltseq(M)]), GF(NextPrime(2^60)));
print "kernel dimension mod p:", Dimension(Kernel(Mp)), Cputime(t0);
ker := Kernel(M); print "kernel dimension over Q:", Dimension(ker), Cputime(t0);
B := Basis(ker);
B := [ v*LCM([Denominator(c) : c in Eltseq(v)]) : v in B ];
B := [ v/GCD([Integers() ! c : c in Eltseq(v)]) : v in B ];
print "digits of the largest coefficient:", [ #IntegerToString(Max([Abs(Integers() ! c) : c in Eltseq(v)])) : v in B ];
L := LLL(Matrix(Integers(), [[Integers() ! c : c in Eltseq(v)] : v in B]));
print "after LLL:", [ #IntegerToString(Max([Abs(c) : c in Eltseq(L[i])])) : i in [1..Nrows(L)] ];
PrintFile("t_pencil.txt", <MA, MB, [Eltseq(L[i]) : i in [1..Nrows(L)]]> : Overwrite := true);
quit;
