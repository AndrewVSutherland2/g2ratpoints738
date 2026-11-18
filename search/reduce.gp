\\ reduce.gp -- reduced models of the curves cut out by a list of lines.  Called by search_model.sh as
\\   gp -q reduce.gp < /dev/null   with the environment variables RED_IN and RED_OUT.
\\ Input rows "A B C f0 .. f6": a line and the sextic f = f0 + f1 x + .. + f6 x^6 it cuts out (from linesearch -exactbatch).
\\ Output rows "A B C g0 .. g6" with g = hyperellred(f / s^2), s^2 the largest square dividing the content of f.
\\ Rows whose sextic has degree less than 5, a repeated root, or is a square are dropped.  The last row is "DONE n".
default(parisizemax, 4000000000);
{
  my(R = readstr(getenv("RED_IN")), out = getenv("RED_OUT"), n = 0, buf = List());
  for(i = 1, #R,
    my(w = select(t -> #t > 0, strsplit(R[i], " ")));
    if(#w < 10, next);
    my(f = Pol(Vecrev(vector(7, j, eval(w[3 + j])))), c = content(f), s = 1, fa = factor(c), g, v);
    for(j = 1, #fa~, s *= fa[j, 1]^(fa[j, 2] \ 2));
    f /= s^2;
    if(poldegree(f) < 5 || poldisc(f) == 0 || issquare(f), next);
    g = hyperellred(f)[1]; v = Vecrev(g); v = concat(v, vector(7 - #v));
    listput(buf, Str(w[1], " ", w[2], " ", w[3], concat(vector(7, j, Str(" ", v[j]))), "\n")); n++;
    if(#buf == 2000, write1(out, concat(Vec(buf))); buf = List()));
  if(#buf, write1(out, concat(Vec(buf))));
  write1(out, Str("DONE ", n, "\n"));
}
\q
