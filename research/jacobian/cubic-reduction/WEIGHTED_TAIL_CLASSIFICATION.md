# Complete19-atom weighted tail classification

Fix F1=a-3a^2b-a^3c, identity linear part and source/output weights
(1,-1,-2). Allow **every** compatible nonlinear monomial of ordinary
degree at most seven in F2,F3: there are eight and nine respectively,
so19 possible nonlinear atoms including the two fixed first-output terms.
Require a collision between (0,beta,c0) and(-1,beta,1+3beta), beta nonzero.

Over every characteristic-zero field, the only such coefficient vector
is the established13-term map in [KCUB10](RAW13_PATTERN_RIGIDITY.md),
with beta=1 and c0=-8. **All six additional tail coefficients vanish.**
Thus these additions cannot give a12-term map. This is a complete
classification of the stated slice, not a global raw-support bound.
Changing the first output, weight pattern or collision configuration
remains outside it.

## Two-variable determinant reduction

On a nonzero a chart put u=ab,w=a^2c. Write

\[
F=(a f(u,w),a^{-1}g(u,w),a^{-2}h(u,w)),\qquad f=1-3u-w.
\]

The allowed source monomials make F polynomial. The coordinate
transformation (a,b,c) to(a,u,w) has determinant a^3. Applying the same
transformation on the target shows

\[
\det JF=\det D(fg,f^2h)/f^2.
\]

The quotient is a polynomial and the identity extends to all source
points. In v=f coordinates, set G(u,v)=g(u,1-3u-v) and similarly H.
The Keller equation is

\[
(G+vG_v)H_u-G_u(2H+vH_v)=1. \tag{1}
\]

The first collision equation is automatic. The second is
G(-beta,0)=-beta; the third uniquely recovers c0, because
F3(0,beta,c0)=c0+I beta^2. The model includes an inverse variable for
beta, so it covers exactly the stated nonzero-beta family.

## Necessary boundary degrees

At v=0, equation(1) reads

\[
G_0H_0'-2G_0'H_0=1. \tag{2}
\]

For polynomials G0,H0, this forces deg G0<=1 and deg H0<=2.
Indeed, if deg H0=2 deg G0, subtract its leading multiple of G0^2;
this preserves(2) and lowers that degree. Otherwise the leading degree
of the left side is deg G0+deg H0-1. In characteristic zero its leading
coefficient does not vanish, forcing the sum of degrees to be one.
The constant-G0 case gives linear H0 directly.

The permitted monomials consequently give

\[
G=G_0+vG_1+v^2G_2,\quad
H=H_0+vH_1+v^2H_2+v^3H_3,
\]

with degrees in u bounded respectively by(1,2,1) and(2,3,1,0).
In particular H3 is constant. The coefficients of v^5,v^4,v^3 in(1)
are

\[
\begin{aligned}
0&=-5G_2'H_3,\\
0&=3G_2H_2'-4G_2'H_2-5G_1'H_3,\\
0&=3G_2H_1'-3G_2'H_1+2G_1H_2'-4G_1'H_2-5G_0'H_3.
\end{aligned} \tag{3}
\]

## A linear-g exclusion

No polynomial completion is possible when g=u+lambda w. Put
kappa=1-3lambda. Equation(1) becomes

\[
(\kappa u+\lambda-2\lambda v)H_u-\kappa vH_v-2\kappa H=1.
\]

If kappa=0, the nonconstant factor lambda(1-2v) divides the left
side, impossible. Otherwise change to q=kappa u+lambda(1-v).
The vector field sends q to kappa q and v to-kappa v, so the equation
is kappa(q H_q-v H_v-2H)=1. Its polynomial solutions have the form

\[
H=-1/(2\kappa)+\sum_{j\ge0} C_j q^{j+2}v^j.
\]

At source(u,w)=(0,0), equivalently(q,v)=(0,1), this contradicts
h(0,0)=0. This argument has no ordinary-degree cutoff.

## Exhaustive coefficient cases

The following three reduced systems have explicit rational unit
identities. All degree and linear eliminations are independently checked.

1. If G2' is nonzero, equations(3) first force H3=0, then H2=0,
   and then deg H1<=1. Returning to u,w leaves
   h=w+I u^2+J u w. The **small-h** certificate excludes this entire
   branch, including all compatible g tails and the collision equations.
   Hence G2'=0; write G2=K.
2. If H3 is nonzero, equation(3) forces deg G1<=1. If K=0,
   it further makes G1 constant, giving the excluded linear-g case.
   If K is nonzero, the u^2 coefficient of the third equation forces
   deg H1<=2. In source coordinates this sets the u^2w term of g,
   and the u^4,u^3w terms of h, to zero. The **high-w** unit identity,
   with an inverse for the two nonzero leading coefficients, excludes
   this case. Thus H3=0.
3. If K is nonzero, equation(3) now gives H2'=0 and deg H1<=2.
   The **mid-w** unit identity excludes this case. Hence K=0.

We have G2=H3=0. The remaining v^3 equation is
G1 H2'-2G1'H2=0. If H2 is nonzero and G1 is nonzero, their leading
degrees must satisfy deg H2=2 deg G1. The bounds force both constant,
again giving linear g. If G1=0, g=u, also excluded. Thus H2=0.

These conclusions remove all six new source monomials:
a^3b^4,a^4bc^2,a^3c^2 in F2 and a^3bc^2,a^4c^3,a^2c^2 in F3.
The remaining normalized coefficient/collision equations are exactly
the already independently certified raw13 equations. They force the
known coefficient vector and beta=1.

## Exact verification and retained failures

[Independent verification](verify_weighted_tail7.py) reconstructs the
full22-variable universal determinant (three source and nineteen
parameter variables), the necessary boundary equations, all branch
parameterizations and three unit identities. It also checks the
abstract top-v identities, the linear-g coordinate change and the
reduction to KCUB10. [Validation](weighted-tail7/validation.json).
The three certificates have15,665 multiplier monomials and maximum
degree seven. [Certificate manifest](weighted-unit-certificates/manifest.json).

The direct19-atom solve and its linearly reduced generic solve timed
out, and the real sparse query returned UNKNOWN. These statuses are
preserved. The classification follows from the exhaustive structural
cases and replayed identities, not from those caps. A later direct
unit lift for a separate theta-zero model timed out after the three
branch certificates were already complete; that unfinished lift is
not used here.

[KCUB13](WEIGHTED_FULL_DEGREE7_CLASSIFICATION.md) subsequently allows all four
new first-output atoms in the beta1 marked class and still proves uniqueness
of the13-term representative. [KCUB14](WEIGHTED_UNMARKED_DEGREE7_CLASSIFICATION.md) subsequently removes
all collision marks in the nonzero-first-coefficient sector. Zero first
coefficients and other weights remain open.
