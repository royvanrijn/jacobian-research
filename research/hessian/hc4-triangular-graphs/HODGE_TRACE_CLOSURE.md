# One Hodge polynomial closes arbitrary critical-value patterns

14 November 2026. A single trace's Hodge--Deligne polynomial gives a
stronger obstruction than its Euler number. It closes every trace
`H(U)+VW` with arbitrary deg H>=2, rather than only pure powers. It also
closes nonlinear old-r-dependent y boundaries and polynomial pullbacks
whose trace map is branched rather than a coordinate triple.

All degrees, boundary multiplicities and retained-dual word lengths are
unrestricted. The argument uses retained Tang and the retained plane
affine-component theorem, plus the explicitly imported basic
Hodge--Deligne invariant. It does not resolve unrestricted HC4, JC2 or
either distinct plane Keller endpoint.

## 1. The additional invariant and its exact input

Write `E(X;u,v)` for the compactly supported Hodge--Deligne polynomial
and put t=uv, a formal variable shorthand. We use additivity for
constructible stratifications and multiplicativity for products. Thus

\[
 E(A^1)=t,\quad E(A^2)=t^2,\quad
 E(C)=t-g(u+v)+1-r                                         \tag{1}
\]

for a smooth connected affine curve whose completion has genus g and
r omitted points. Every such affine curve has r>=1.

These standard statements are explicitly imported from
[Ricolfi, Sections 1.1 and 3.2, equation (3.6)](https://arxiv.org/html/1907.08123v1),
where the E-polynomial is the ring homomorphism from the Grothendieck
ring of complex varieties and the smooth projective curve formula is
given. Removing r points proves (1). The
[import record](sources/hodge-trace-import.json) pins the inspected
[primary v1 HTML](sources/ricolfi-hodge-1907.08123v1.html).

We use no multiplicativity for arbitrary finite covers. Every product
below comes from an actual trace-stratum isomorphism. We also do not
identify a general variety's whole mixed Hodge structure from its
E-polynomial. Only the constrained curve formula (1) is used for that
direction of the argument.

## 2. A single-trace criterion with no critical-point hypothesis (HC4MYHC0)

Let f(s,z1,z2,z3) be a complex four-variable polynomial. Suppose for
some scalar s0 and generic lambda its entire trace value fibre satisfies

\[
                 E(f(s_0,-)=\lambda)=t^2+\kappa t,
                             \qquad \kappa>0.               \tag{2}
\]

Then **nonzero constant Hessian determinant forces a polynomial gradient
inverse**. No generic-base Euler, isolated trace point, finite projection
or rationality hypothesis is imposed. Actual polynomial equivalences
of the three-variable trace, and actual retained-dual automorphisms over
C[s], preserve (2).

Assume a unit-Hessian nonautomorphic gradient. Retained
[HCTANG2](../tang-hc4-audit/STRUCTURAL_REVIEW.md) supplies a constant
jointly affine source two-plane. Adapt it to the free base using constant
affine changes, preserving the trace being compared.

### Rank zero: the t coefficient cannot be positive

The form is

\[
 f=D(s,U)+P F(s,U)+Q G(s,U),\qquad J(F,G)\in C^\times.
\]

At s=s0, F0,G0 cannot both be the zero polynomial: otherwise both
are divisible by s-s0, and their Jacobian is divisible by that factor.
Let g(U) be their nonzero gcd. A Bezout row completion over C[U] puts
the **trace** in the polynomial coordinates

\[
                       D_0(U)+g(U)V,
                 \qquad W\text{ free}.                     \tag{3}
\]

This nonlinear coordinate change identifies trace varieties only; it
is not used as Hessian invariance. If g has m distinct roots, choose
lambda avoiding the finitely many D0 values at those roots. The value
fibre of (3) is `(A1 minus m points) times A1`. Hence

\[
                         E=t^2-mt,\qquad m\ge0,             \tag{4}
\]

contradicting (2).

### Rank one: genus and missing points force actual affine lines

The form is

\[
 f=D(U,Z)+sA(U,Z)+VB(U,Z),\qquad J(A,B)\in C^\times.        \tag{5}
\]

The **whole** relative critical curve is isomorphic to the smooth,
nonempty plane curve C=B=0: its other equations solve for s,V by
the unit Jacobian. For generic lambda, let d0 be the number of points
`D+s0*A=lambda` on C, summed over all nonconstant components. Values
on constant components are avoided.

Off C the trace value equation solves V, giving A2 minus C. At those
d0 points it contributes one affine V line each. Therefore

\[
                E(f(s_0,-)=\lambda)=t^2-E(C)+d_0t.           \tag{6}
\]

Let c be the number of components of C, G the sum of their completion
genera and R the total number of omitted points. Formula (1) gives

\[
 E(C)=ct-G(u+v)+c-R,
\]
\[
 E(f(s_0,-)=\lambda)=t^2+(d_0-c)t+G(u+v)+(R-c).              \tag{7}
\]

Comparing with (2) forces G=0 and R=c. Each component omits at least
one point, so every component has genus zero and exactly one omitted
point: it is an **actual closed A1**, not just a rational normalization
of a singular curve. The retained
[HC4MYGR0 theorem](CRITICAL_COMPONENTS.md#1-a-plane-keller-fibre-with-an-affine-line-component)
makes (A,B) a polynomial automorphism. Its two outputs recover U,Z;
the remaining gradient outputs recover s,V by the inverse transpose
of its Jacobian after subtracting grad D. Thus grad f is an automorphism,
contradicting the hypothesis. Both ranks are closed.

More generally, if a unit-Hessian member with nonautomorphic gradient
has a trace E=t^2+kappa*t for any integer kappa, it must be in rank zero
relative to that base, and kappa=-m<=0 as in (4). The positive condition
in (2) is essential to this criterion. The retained relative-coordinate
endpoint has trace fibres A2 and kappa=0; it is not closed here.

## 3. All traces linear in one variable with a positive root budget (HC4MYHC1)

Suppose a trace is polynomially equivalent, with actual inverse, to

\[
                         H(U,V)+a(V)W,
                        \qquad a\ne0.                       \tag{8}
\]

Let rho1,...,rhom be the **distinct** roots of a. Set
`d_i=deg_U H(U,rho_i)`, using d_i=0 if that polynomial is constant.
For generic lambda, outside these m V values the equation solves W
with U free. At V=rho_i it has d_i distinct U roots with W free;
if d_i=0, avoid its constant value and it contributes nothing. Hence

\[
       E((8)=\lambda)=t^2+\left(\sum_i(d_i-1)\right)t.        \tag{9}
\]

If the sum is positive, HC4MYHC0 proves gradient invertibility under
unit Hessian, for arbitrary off-hyperplane terms and trace equivalences.
Repeated roots of a count once, not with multiplicity. Degree drops
of H at individual rho_i and constant specializations are retained.

In particular every trace `H(U)+bVW`, b nonzero and deg H=d>=2,
has E=t^2+(d-1)t. Thus **all univariate H are closed, with arbitrary
critical points and critical values**. This strengthens the earlier
[pure-power trace corollary HC4MYHB1](ISOLATED_TRACE_MIXING.md).
The example `f=sU+H(U)+bVW` retains unit Hessian and an explicit
triangular gradient inverse. No unit-Hessian nonexistence is asserted.

## 4. Nonlinear old-r-dependent y boundaries (HC4MYHC2)

Take any genuine graph y=R(x,p,q,r), followed by any actual retained
three-dual automorphism over K[x]. All off-x-zero terms are arbitrary.
Write R0=R(0,p,q,r). Either of the following boundary classes is closed
for nonautomorphic unit-Hessian gradients.

### A. Any polynomial of a unit-affine r coordinate

Suppose

\[
 R_0=P(ar+b(p,q)),\qquad a\in K^\times,\quad \deg P\ge1.    \tag{10}
\]

The trace coordinate change w=ar+b(p,q) fixes p,q. The native trace is

\[
 h=Lp^2+Mp+4NpP(w)^2+NqP(w).                               \tag{11}
\]

This is (8), with U=p,V=w,W=q, coefficient a(V)=N*P(V) and
degree d_i=2 at every distinct root of P. Therefore E=t^2+m*t,
where m>=1 counts those distinct roots. Arbitrary nonlinear P and all
its multiplicities are permitted. The former unit-affine r boundary
is just the degree-one case. An actual retained-dual automorphism
preserves this entire trace fibre even when the displayed shape changes.

### B. Any polynomial of a nonvertical plane coordinate

Suppose instead

\[
 R_0=P(C(p,r)),\qquad (C,D)\text{ a plane polynomial
 automorphism},\quad \deg P\ge1,                            \tag{12}
\]

with C not an affine polynomial in p alone. For the m distinct roots
alpha_j of P, the zero curves Gamma_j=C=alpha_j are disjoint closed
copies of A1. On each, p is a nonconstant polynomial in its parameter
D. Indeed, if p were constant rho, that curve would be the line p=rho;
its prime equation C-alpha_j would then be a unit multiple of p-rho,
making C affine in p, contrary to the assumption.

Put `e_j=deg_D p|Gamma_j>=1`. Off the union of Gamma_j in the (p,r)
plane, the native trace equation solves q, giving A2 minus those m
copies of A1. On Gamma_j its value is Lp^2+Mp, so a generic value
has2e_j points, each supplying a free q line. Therefore

\[
             E(h=\lambda)=t^2+
                        \left(2\sum_j e_j-m\right)t,        \tag{13}
\]

whose coefficient is at least m>0. HC4MYHC0 applies. The degrees e_j
may differ between components; equality is not assumed. This also includes
plane coordinates nonlinear in r outside (10).
The excluded affine-p coordinate has e_j=0 instead, so this sign
argument does not cover it.

## 5. Branched trace triples and arbitrary mixing (HC4MYHC3)

In actual affine source coordinates (s,z1,z2,z3), substitute arbitrary
polynomials (X,Y,P,Q,R) into the scaled Meng--Yang potential. Assume
X(0,z)=0. Suppose a trace coordinate triple (u,v,w) is an actual
polynomial automorphism of A3, and on the hyperplane s=0 one has

\[
                    P=u,\quad Q=v,\quad Y=J(w),
                           \deg J\ge1.                     \tag{14}
\]

Then the trace is (11), and every actual unit-Hessian pullback has a
polynomial gradient inverse. All off-hyperplane corrections, including
X=s times any polynomial, and all values of R are unrestricted.

For deg J>1 the displayed physical triple (Y,P,Q) is a finite branched
map and has no polynomial inverse. This is intentionally a larger
trace class than the coordinate-triple hypothesis of HC4MYHB3.
It does not replace that hypothesis by a merely Keller assumption:
the branched trace is admitted by its explicit positive Hodge budget.
No global inverse or embedding is asserted for the broader substitutions.

The source hyperplane must remain affine in the actual Hessian variables.
A nonlinear source divisor is not straightened using nonexistent Hessian
invariance. Trace automorphisms identify whole value varieties only.

## Descent, checks and remaining scope

The arguments are first over C. Characteristic-zero candidates and all
required coordinate inverses use finitely many coefficients; embed that
field into C and descend any polynomial gradient inverse by uniqueness.
Neither graph degrees nor word lengths acquire a bound in this passage.

Remaining old-r-dependent y boundaries must avoid both (10) and (12),
as well as the earlier closures. Other opposite q/r-dependent bases,
mixing outside the trace criteria, uncovered embeddings and both distinct
plane Keller endpoints remain open. A zero or negative Hodge t coefficient
is not relabeled a positive-budget theorem. No unrestricted HC4/JC2
solution, formal verification, independent external review or novelty
priority is claimed.

[verify_hodge_trace.py](verify_hodge_trace.py) and its
[receipt](receipt-hodge-trace.json) check the polynomial strata, both
Tang-rank formulas, genus/puncture coefficient identities, nonlinear
boundary coordinate controls and retained input hashes. The all-degree
curve and rank arguments are written proofs. The basic Hodge invariant
is an explicit external mathematical input, not proved by finite controls.
