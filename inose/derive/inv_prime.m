// magma -b I:=<n> inv_prime.m : echelon basis of the solutions modulo the n-th prime below 2^61
load "invlib.m";
n := StringToInteger(I); p := 2^61; for i in [1..n] do p := PreviousPrime(p); end for;
KK := 34; nn := #Exps2(KK) + #Exps2(KK-3);
ker, EA, EB := InvKernel(p, KK, (3*nn) div 2 + 60 : Seed := 1000 + n);
B := BasisMatrix(ker);
fn := "k34_p" cat I cat ".txt";
PrintFile(fn, p : Overwrite := true); PrintFile(fn, Nrows(B));
SetColumns(0);
for i in [1..Nrows(B)] do PrintFile(fn, [Integers() ! B[i,j] : j in [1..Ncols(B)]]); end for;
PrintFile(fn, "DONE");
quit;
