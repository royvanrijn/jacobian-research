# Relative critical loci close every y=a(x)p+b(x) reparameterization

> The later [critical-component theorem](CRITICAL_COMPONENTS.md),
> `HC4MYGR0/1/2/3`, closes nonlinear y=R(x,p) with a simple boundary
> zero, and excludes every opposite x=a(y)p+b(y) with a(0) nonzero.
> In particular it excludes the matching-critical-determinant control
> below by its reducible critical geometry, without mistaking that
> determinant test for sufficiency.

14 November 2026. **HC4MYGK2** proves that every graph
y=a(x)p+b(x), with nonzero polynomial a and arbitrary polynomial b,
has no unit-Hessian pullback after any polynomial automorphism of
the three retained dual variables over K[x]. Degrees and word lengths
are unrestricted, and a need not be a unit or avoid zeros.

This is stronger than the [native-chart invertibility theorem](BASE_ELIMINATED_GRAPHS.md)
on this class. The proof is elementary and uses neither Tang nor a
plane Jacobian conjecture. The key is a Hessian invariant on the
relative critical locus, where nonlinear chain-rule corrections vanish.

## 1. Preserve the critical locus, not the whole Hessian

For f(s,z1,z2,z3), write

\[
 I_f=(f_{z1},f_{z2},f_{z3}),\qquad
 \mathcal C_f=V(I_f).
\]

Let Phi(s,Z)=(s,phi(s,Z)), where phi is a polynomial automorphism
over K[s], and put g=f composed with Phi. Its relative Jacobian
determinant is a field unit d. The relative gradient ideals satisfy

\[
 I_g=\Phi^*I_f,                                 \tag{1}
\]

since the relative coordinate Jacobian and its inverse are polynomial.
The full Hessian chain rule is

\[
 \operatorname{Hess}g=J_\Phi^T\operatorname{Hess}f\,J_\Phi
       +\sum_i (f_{zi}\circ\Phi)\operatorname{Hess}\phi_i.
\]

There is no correction from f_s, since the base coordinate is unchanged.
Thus, as a polynomial identity modulo I_g,

\[
 \det\operatorname{Hess}g
   \equiv d^2\,\Phi^*(\det\operatorname{Hess}f)\pmod {I_g}. \tag{2}
\]

**HC4MYGK0.** The relative critical schemes are isomorphic under Phi,
and their Hessian determinant classes differ only by the nonzero
constant factor d^2. In particular a nonconstant determinant on any
critical component, a singular critical point, or different constant
values on two components prevents unit Hessian after every such
automorphism. An actual polynomial inverse is needed to transport
the entire critical locus; a merely Keller substitution is not used.
Constant affine changes of the free base are harmless as well.

Equation (2) is not nonlinear invariance away from the critical locus.
For example f=s*z1+z2*z3+z1^3 has Hessian determinant one. The relative
automorphism phi=(P+Q^2+sZ,Q+Z^2,Z) preserves that determinant on the
critical curve, but the pulled-back determinant at (s,P,Q,Z)=(0,1,0,0)
is -323. This is an exact control on the scope of (2).

## 2. The complete critical curve for y=R(x,p), away from x=0

Return to the original Meng--Yang family Psi=LA^2+MA+NB, LN nonzero,
and take any polynomial graph y=R(x,p), independent of q,r. Use x
as the free base and (p,q,r) as the retained duals. On x nonzero use
the rational computation coordinates (p,T=A,q). Eliminating r gives
the literal identity, with u=1+xR,

\[
 f=LT^2+\left(M-\frac{N(2-3xR)}{x^2}\right)T
       +\frac{Nu(u+1)p}{x^2}+\frac{2N(2u+1)}xq. \tag{3}
\]

These coordinates are used only to compute a critical branch. They
are not asserted to be a global polynomial parameterization.

The q derivative forces 2u+1=0. The T derivative and then the p
derivative give exactly

\[
 2xR(x,p)+3=0,\qquad R_p\ne0,
\]
\[
 T_c=\frac{13N/(2x^2)-M}{2L},\qquad
 q_c=-\frac{3T_c}{4x}+\frac1{16x^2R_p}.          \tag{4}
\]

r is recovered polynomially in p,T,q over K[x,x^(-1)]. If R_p=0
on the displayed curve, the p derivative of (3) equals -N/(4x^2),
so there is no critical point there. Thus (4) classifies every
relative critical point with x nonzero.

On this curve the critical value is

\[
 D=-\frac{Np}{4x^2}
       -\frac{(M-13N/(2x^2))^2}{4L}.             \tag{5}
\]

The derivation along the curve is

\[
 \delta=\partial_x+
       \frac{3/(2x^2)-R_x}{R_p}\partial_p.        \tag{6}
\]

It is well defined on every component after inverting x*R_p, since
delta(2xR+3)=(2xR+3)/x. In the rational coordinates (p,T,q), the
dual Hessian determinant on (4) is -32*L*N^2*R_p^2. The relative
coordinate Jacobian from (p,q,r) is x^3. The Schur complement for
the remaining x row is exactly delta^2 D. Consequently

\[
 \left.\det\operatorname{Hess}f\right|_{\mathcal C_f}
      =-32LN^2x^6R_p^2\,\delta^2D.              \tag{7}
\]

This proves **HC4MYGK1**, an exact all-degree critical-curve and
determinant formula. HC4MYGK0 requires (7) to be the same nonzero
constant on every such component if any retained-dual automorphism
is to give unit Hessian. It is a necessary condition, not a proof
that a graph passing it has any unit-Hessian parameterization.

## 3. Exclude every y=a(x)p+b(x), with arbitrary polynomial coefficients

Let a be a nonzero polynomial. The curve in (4) is the rational branch

\[
 p_c=\frac{-3/(2x)-b(x)}{a(x)},\qquad x a(x)\ne0.
\]

On it (5) becomes

\[
 D(x)=\frac{N(3+2xb(x))}{8a(x)x^3}
            -\frac{(M-13N/(2x^2))^2}{4L}.       \tag{8}
\]

In particular its denominator divides a(x)*x^4: explicitly

\[
 D=\frac{2LNx(3+2xb)-a(2Mx^2-13N)^2}{16La x^4}. \tag{9}
\]

If an actual relative dual automorphism made the Hessian determinant
a nonzero constant, (2) and (7) would force

\[
 D''=\frac{K}{a(x)^2x^6},\qquad K\ne0.          \tag{10}
\]

Write d=deg a. At infinity the right side has exact order 2d+6.
The polynomial part of the rational function D must therefore have
degree at most one. Subtract that affine part and let n be the
vanishing order of the remaining nonzero proper rational function
at infinity. Characteristic zero gives n+2=2d+6, hence n=2d+4.
But (9) bounds its denominator degree by d+4, so n<=d+4. Thus d=0.
This proves the contradiction for every nonconstant a without a
coefficient bound, root enumeration or assumption on deg b.

For constant a nonzero, the x^(-3) coefficient of (8) is 3N/(8a).
The polynomial b term has powers at least x^(-2) and cannot cancel
it. Therefore D'' has nonzero x^(-5) coefficient 9N/(2a), whereas
(10) has only the x^(-6) term. This final case is also impossible.

This proves **HC4MYGK2** over every characteristic-zero field:
no graph y=a(x)p+b(x), a nonzero polynomial, admits a nonzero
constant Hessian after any polynomial retained-dual automorphism
over K[x]. Zeros of a cause no gap: the contradiction is an identity
on the nonempty rational open branch x*a(x) nonzero. No inverse of
p->a(x)p+b(x) is asserted globally.

If a is identically zero this argument does not apply. That is the
base-only family already treated by
[HC4MYGU0/3](RELATIVE_COORDINATE_ENDPOINT.md): b(0)=0 forbids unit
Hessian, while b(0) nonzero gives the unresolved unit-ideal plane
Keller endpoint. The exception is mathematically essential.

## 4. Controls and the remaining boundary

The argument is asymmetric in the two base coordinates. As an exact
control take L=N=1, M=0 and the opposite graph

\[
 x=p+\frac{1013}{54}y^2.                        \tag{11}
\]

Its generic relative critical branch has x=-3/(2y) and critical
value D(y)=y/6-y^4/486. The native dual Hessian determinant there is
-162/y^2, so its full critical determinant is 4. The boundary
critical line y=p=r=0 also has determinant 4. Nevertheless the full
Hessian determinant at (y,p,q,r)=(0,1,0,0) is 1556. These computed
critical components do not certify a unit Hessian or a polynomial
reparameterization achieving one. They prevent treating a passed
branch condition, or an exchange of x and y, as the theorem in
Section 3. No closure of all opposite-chart affine-p graphs is claimed.
The later HC4MYGR3 nevertheless excludes this control and the full
a(0) nonzero subclass by their critical components.

For a general graph y=R(x,p), the next exact gate is the entire
localized curve 2xR+3=0 and the function (7), including every
component and the remaining x=0 critical locus. Graphs involving
q,r in the omitted base coordinate, opposite-chart reparameterizations
outside the proved cases, and nonlinear changes mixing the free base
with retained duals remain. The unit-ideal omitted-point endpoint and
the distinct nonlinear-p singleton-zero properness endpoint also remain.

The [checker](verify_relative_critical.py) verifies the actual Hessian
chain-rule remainder, relative inverse controls, the critical curve,
all original-family rational identities, both Schur determinants,
the denominator bound, and the opposite-chart controls. The infinity
order proof and its unrestricted degree quantifiers are written
mathematics, not a finite degree test. [Actual receipt](receipt-relative-critical.json).

```sh
timeout 60s .venv/bin/python research/hessian/hc4-triangular-graphs/verify_relative_critical.py
```

There is no new external classification input. No unit-ideal Keller
pair, plane surjectivity or singleton-zero properness theorem,
formal verification, independent external review, or unrestricted
HC4/JC2 solution is asserted.
