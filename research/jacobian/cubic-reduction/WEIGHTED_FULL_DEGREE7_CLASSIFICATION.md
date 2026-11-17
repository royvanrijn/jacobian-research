# Full marked degree-seven weighted classification

**Theorem (KCUB13).** Let F be a polynomial map in three variables over a
characteristic-zero field. Assume identity linear part, det JF=1,
source/output weights (1,-1,-2), ordinary degree at most seven,
coefficients -3 on a^2b and -1 on a^3c in F1, and a collision between
(0,1,c0) and (-1,1,4). Allow every compatible nonlinear monomial.
Then F is the established 13-term map

\[
\begin{aligned}
F_1&=a-3a^2b-a^3c,\\
F_2&=3(1+2ab)^2ac+
 \bigl(9(1+2ab)^2-6(1+2ab)-2\bigr)b,\\
F_3&=(1+2ab)^3c+
 \bigl(6(1+2ab)^2+2(1+2ab)\bigr)b^2,
\end{aligned}
\]

and c0=-8. There are 23 allowed nonlinear atoms before solving;
all ten additions to this 13-term support vanish. Hence 13 is the
minimum raw nonlinear count **in this marked weighted class**, and
the coefficient vector itself is unique. This is not a global sparsity
minimum or an unrestricted three-variable degree-seven classification.

The conclusion follows over the algebraic closure, so it also holds
over the original characteristic-zero field.

## Universal quotient equation

Set u=ab,w=a^2c and write F=(af,a^{-1}g,a^{-2}h).
Weight-compatible source atoms make F polynomial. Direct differentiation
on the a!=0 chart gives the polynomial identity

\[
J=g(f_u h_w-f_w h_u)+2h(g_u f_w-g_w f_u)
    +f(g_u h_w-g_w h_u)=\det JF. \tag{1}
\]

The independent checker reconstructs the full universal three-variable
determinant and verifies (1); the quotient formula is not trusted as input.
The first collision equation forces

\[
f=1-3u-w+R u^3+P u^2+Q u w+T w^2,
\qquad R=P-4Q+16T. \tag{2}
\]

P,Q,T are called f2,f3,f4 in the saved input. The second collision
eliminates g1 linearly. The model has 21 coefficient parameters and
19 free parameters after these two complete linear eliminations.
The third collision equation only determines c0: the coefficient of c
at a=0 is one. The eight g-atoms and nine h-atoms are all retained.

## A smooth cubic boundary cannot occur

We first prove RT=0, without a coefficient-ring elimination.
If f=f_u=f_w=0 anywhere, (1) gives J=0. Thus the affine curve C:f=0
must be smooth for a Keller completion.

Suppose RT!=0 and put y=2Tw+Qu-1=f_w. Completing the square gives

\[
y^2=(Qu-1)^2-4T(Ru^3+Pu^2-3u+1)=:c(u). \tag{3}
\]

Here c has degree three, with leading coefficient -4TR. A multiple
root of c would give a singular point of C, already impossible.
Consequently the smooth projective completion E has exactly one
point O outside C. These assertions can also be seen directly from
Y^2Z=c3 U^3+c2 U^2Z+c1 UZ^2+c0 Z^3: its sole point at Z=0 is
[0:1:0], smooth since the derivative with respect to Z is one there.
In its Y=1 chart, Z=c3 U^3+terms of higher order. Therefore

\[
\operatorname{ord}_O u=-2,\qquad
\operatorname{ord}_O y=-3,\qquad
\operatorname{ord}_O(du/y)=0.
\]

The differential omega=du/y is regular and nowhere zero everywhere
on E. At finite points with y!=0 this is immediate. At a simple root
r of c, use y as local parameter: u-r=y^2/c'(r)+higher terms, so
omega has a nonzero constant leading coefficient. The preceding
orders give the same fact at O. This supplies the local facts needed
below without assuming them from a genus label.

Restrict g,h to C and put rho=h/g^2. The tangent derivation
f_w partial_u-f_u partial_w and (1), restricted to f=0, give

\[
g\,dh-2h\,dg=-\omega,
\qquad d\rho=-\omega/g^3. \tag{4}
\]

In particular g is not identically zero. If g is a nonzero constant,
(4) says dh is regular everywhere on E. Any pole of h would give a
pole of dh in characteristic zero. Thus h is regular on the complete
curve, hence constant, contradicting (4).

Otherwise g has its only pole at O, of some order m>0. Its finite
zeros have total multiplicity m. At a zero of order k, the right
side of (4) has pole order 3k. The derivative of a rational function
has pole order one more than that function's pole order. Therefore
rho has a pole of order 3k-1. But h is regular there, so h/g^2 has
pole order at most 2k. Hence 3k-1<=2k forces k=1. All m zeros of g
are simple, and rho has exactly m double poles, total degree 2m.

At O the differential d rho has zero order 3m. Hence rho is regular
there and rho-rho(O) has zero order 3m+1. A nonconstant rational
function on a smooth complete curve has equal total zero and pole
degrees; here this would require 3m+1<=2m, impossible. This proves
RT=0. The argument is algebraic and applies after any algebraically
closed characteristic-zero base extension.

## All residual parameter branches

Every identity below is stored with explicit rational polynomial
multipliers and independently multiplied, not inferred from a solver
status. Write I for the determinant and collision coefficient equations.

| Branch | Replayed polynomial in the necessary ideal | Consequence |
|---|---|---|
| T=0 | Q^16 | Q=0 |
| T=Q=0 | P^3, after the boundary constraints below | P=0 |
| R=0 | T^2(Q-8T)^14 | T=0 or Q=8T |
| R=0, Q=8T, P=16T | T^8(T-7/4), after boundary constraints | T=0 or T=7/4 |
| Previous line, T=7/4 | 1 | This exceptional value is impossible |

For the two rows requiring boundary constraints, f=0 is a polynomial
line with the following parameterizations:

\[
\begin{array}{c|cc}
 &u(t)&w(t)\\ \hline
T=Q=0&t&1-3t+Pt^2+Pt^3\\
R=0,\ Q=8T&-1+t-Tt^2&t-4(-1+t-Tt^2).
\end{array}
\]

They extend to polynomial coordinates (t,v) with v=f and Jacobian
partial(u,w)/partial(t,v)=-1. Thus G(t)=g(u(t),w(t)) and
H(t)=h(u(t),w(t)) satisfy

\[
G H'-2G'H=1. \tag{5}
\]

For completeness, (5) forces deg G<=1, deg H<=2. G cannot be zero.
If deg H differs from twice deg G, the leading coefficient in the
left side has degree deg G+deg H-1 and nonzero coefficient; a
constant result forces the degree sum to be one. If deg H=2deg G,
subtract a suitable constant multiple of G^2 from H; (5) is unchanged
and H's degree decreases below 2deg G, reducing to the previous
case. Constant G gives linear H directly. Consequently all
coefficients of t^k for k>1 in G and k>2 in H vanish. These are
necessary linear conditions on the tail coefficients. The checker
reconstructs them, their coordinate identities and (5), before
using any certificate containing them.

The branch table is exhaustive: RT!=0 was excluded geometrically.
When T=0 its first two rows force Q=P=0. When R=0 and T!=0, the
last three rows give a contradiction, including the sole nonzero
exceptional value. Therefore P=Q=T=R=0.

## Reduction, verification and retained failures

With these four first-output additions zero, the map lies in the
already classified [19-atom tail class](WEIGHTED_TAIL_CLASSIFICATION.md)
(KCUB11), with beta=1. That classification forces the displayed map,
its 13 nonzero nonlinear coefficients and c0=-8.

Replay only the independent audit with:

```sh
.venv/bin/python research/jacobian/cubic-reduction/verify_weighted_newfirst7.py
```

[validation.json](weighted-newfirst7/validation.json) records four
cleared parameter identities and one exceptional unit identity:
1,834 multiplier monomials, maximum multiplier degree 15. Its
standard-library Fraction engine reconstructs all 23 atoms, the
universal determinant and linear eliminations, the quotient formula,
the completed square, both boundary charts and Wronskians, every
necessary linear constraint, and all five polynomial identities.
The smooth projective curve and pole-order argument above is a
written mathematical proof, **not a formalized machine proof**.

The original sparse real query remains UNKNOWN at 30 seconds.
Unaugmented cubicgraph and parabolic parameter-basis queries timed
out at 20 seconds. Boundary constraints resolve those branches;
no timeout is relabeled UNSAT. A polynomial-renderer failure after
rational linear elimination is preserved as error-v0.json, alongside
the corrected source and input snapshots. A generic full-field unit
certificate was also exported but has an additional large exceptional
factor; it is exploratory and unnecessary for this proof. The boundary
argument avoids any inference from that generic result.

This closes the 23-atom marked search and supersedes its old open
handoff. Changing weights, fixed coefficients or collision marks,
and all global dimension and sparsity minima, remain outside the theorem.

[KCUB14](WEIGHTED_UNMARKED_DEGREE7_CLASSIFICATION.md) subsequently removes
all collision assumptions and extends the result to any nonzero pair of
the two first-output coefficients by diagonal conjugacy. Changing those
marks within the nonzero sector is now closed; zero coefficients and
other weights remain outside the broader theorem.
