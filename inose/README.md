# The surface: Inose's equation and Elkies' sextic

The K3 surface X_163 has an equation that needs no table of coefficients: Inose's elliptic fibration attached to the
elliptic curve with complex multiplication by the ring of integers of Q(sqrt(-163)), whose j-invariant is -640320^3.
With a = 53360 = 640320/12 it is

    y^2 = x^3 - 3 a^2 x - (t + 1/t + 2 a^3 + 2)

(fibres of type II* at t = 0 and t = infinity, of type I2 at t = -1).  The search starts from Elkies' sextic
`../data/elkies163.txt`.  This directory holds an explicit birational map between the two surfaces, in both directions,
with a proof of each, so that the input of the search is tied to the equation above by a computation that takes seconds.

## The surfaces

    Inose:  y^2 = x^3 - 3 a^2 x - (t + 1/t + 2 a^3 + 2),  a = 53360
    Y:      w'^2 = t (u^3 (t+1) - 3 a u^2 t - t^2)       the same surface as a double plane:
                                                         x = u (t+1)/t - a,  y = w' (t+1)/t^2
    X:      w^2 = S(x0,x1,x2)                            Elkies' sextic, node at (0:0:1)

## The maps

* `map_data.m` — X -> Y:  t = (At + Bt w)/Gb,  u = (A0 + B0 w)/(Ai + Bi w),  w' = w J/8, where J is the Jacobian
  determinant of (t, u) with respect to (x1, x2) on the chart x0 = 1.  At, Bt are forms of degree 32 and 29 (coefficients
  of up to 58 digits), A0, B0, Ai, Bi of degree 24 and 21 (46 digits); Gb is a product of eight lines and a conic
  (given in `verify.m`).
* `inverse_data.m` — Y -> X:  (x0 : x1 : x2) = (X0 : X1 : X2) with X_j = XjA(t,u) + XjB(t,u) w', of degree 34 and 31
  (71 digits).
* `verify.m` — proves the first map: the equation of Y holds identically on X (both sides are polynomials of degree at
  most 244 in x1 and they agree on 260 lines x1 = c); the fibre t = 0 is eight lines and a conic in the configuration of a
  fibre of type II*; the map is birational (the fibres of t are connected and u has degree 2 on them).
  `magma -b verify.m` (4 seconds); output in `verify.log`, ending with INOSE MAP CHECK PASSED.
* `verify_inverse.m` — proves that the second map is the inverse of the first: phi(psi(q)) = q exactly at 102 points on
  each of 100 fibres, which is more than the degree of the difference on a fibre; the argument is in the file.
  `magma -b LO:=1 HI:=100 verify_inverse.m` (the log is the output of 20 pieces of 5 fibres, about a minute on 16 cores).

Run both from this directory.

## What the fibration looks like on Elkies' sextic

    zero section        the exceptional curve of the node (0:0:1)
    fibre t = 0         a0 + 2 a1 + 3 a2 + 4 a3 + 5 a4 + 6 a5 + 4 a6 + 2 a7 + 3 a8, one lift each of
                        a0: 9 x0 - x1 (through the node)   a1: 80 x0 - 11 x1 + x2     a2: 72 x0 - 27 x1 + x2
                        a3: 112 x0 - 49 x1 + 3 x2          a4: 33 x0 - 22 x1 + x2     a5: 80 x0 - 57 x1 + 3 x2
                        a6: 108 x0 - 99 x1 + 5 x2          a8: 24 x0 - 19 x1 + x2
                        a7: the conic 768 x0^2 - 1960 x0 x1 + 100 x0 x2 + 1152 x1^2 - 117 x1 x2 + 3 x2^2
    fibre t = infinity  the image of that under (x1,x2) -> (-x1,-x2); this symmetry of the sextic is t -> 1/t
    fibres              curves of degree 32;  At^2 - Bt^2 S = Ga Gb, the products of the components of the two fibres
    fibre t = -1        two components, of degrees 8 and 24

In the basis F (fibre), O (zero section), a1..a8, b1..b8 (the components of the two II* fibres that do not meet O),
r (the component of the I2 fibre that does not meet O) and P (the section of height 163/2 that generates the
Mordell-Weil group), the class of a line of Elkies' plane is

    D_E = -18 F + 30 O - (31,61,90,118,145,98,50,73).a - (31,61,90,118,145,98,50,73).b - 11 r + 2 P.

Its intersection number with the class of a line of the plane of Y is 42.

## How the two data files are computed

`derive/derive.sh WORKDIR` recomputes both files (about two minutes on 16 cores) and reproduces them byte for byte.
The function t is the quotient of two sections of 32 D_E that vanish to prescribed orders on lifted lines, the node curve
and a conic (`pencil.m`, `t_q.m`; exact linear algebra), and u is obtained in the same way from sections of 24 D_E
(`u_q.m`); `fit.m` normalises them and `mkdata.m` writes `map_data.m`.  The inverse is found by interpolation modulo 32
primes (`invlib.m`, `inv_prime.m`), rational reconstruction (`crt.py`) and lattice reduction (`lll.m`, `mkinv.m`).
Neither proof depends on how the files were found: `verify.m` and `verify_inverse.m` use only the two data files.
