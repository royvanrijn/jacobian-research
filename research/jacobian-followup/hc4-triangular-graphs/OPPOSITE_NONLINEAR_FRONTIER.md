# Opposite nonlinear graphs: global closure without finite-base escapes

Put `x=R(s,p), y=s` in the Meng--Yang potential, with `L*N != 0`.
Allow every actual polynomial automorphism of the retained three duals
over K[s]. Degrees and word lengths are unrestricted. Work first over C;
characteristic-zero nonexistence and polynomial invertibility descend
through the finitely generated coefficient field as in the preceding notes.

**HC4MYGX0. A simple zero of R(0,p) forces gradient invertibility under
actual unit Hessian, for arbitrary nonlinear R.** With an active critical
curve as well, unit Hessian is impossible. This extends the earlier
opposite affine-p boundary theorem without assuming linear p dependence.

**HC4MYGX2. Every case with no finite-base escapes of the active critical
curve is closed for nonautomorphic unit-Hessian gradients.** The proof
handles all cyclic covering degrees and inactive multiplicities at once.

**HC4MYGX3. Every remaining nonautomorphic unit-Hessian opposite graph
must have a finite-base escape and at least six distinct generic roots of
`3+2sR` in p.** Thus at most five distinct generic roots already force
gradient invertibility, even if repeated factors make deg_p R arbitrarily
large. Six is a necessary distinct-root threshold, not a bounded search
proposal or an assertion of existence at that threshold.

Base-only graphs retain their exact unit-ideal plane-surjectivity endpoint.
The general opposite nonlinear chart with escapes is not fully closed here.

## 1. Critical geometry, the boundary and every value fibre (HC4MYGX0)

Write `z=R(s,p)` for the old x coordinate and put `F=3+2sR`. On z nonzero,
eliminating the old r through T=A gives the literal potential

\[
 LT^2+\left(M-\frac{N(2-3zs)}{z^2}\right)T
       +\frac{N(1+zs)(2+zs)p}{z^2}+\frac{2N(3+2zs)}zq.
\]

Its complete relative critical set away from s=0 is

\[
 F=0,\quad R_p\ne0,\quad z=-3/(2s),\quad
 T=\frac{13N/(2z^2)-M}{2L},\quad
 q=\frac{17T}{12z}+\frac{p}{12z}-\frac1{24R_p},                 \tag{1}
\]

with r recovered from A. If R_p=0 the remaining p derivative is
`-N/(4z^2)`, so such points are not critical. At z=0 with s nonzero
the q derivative is Ns, so no additional critical points are missed.

The critical value and full native Hessian class are

\[
 D=-\frac{Ns^2p}{9}-\frac{(M-26Ns^2/9)^2}{4L},\qquad
 K=-162LN^2\frac{R_p^2}{s^2}D''.                               \tag{2}
\]

Prime means differentiation along the critical curve with respect to s.
The relative (p,T,q) Hessian determinant is `-72LN^2 R_p^2/z^4`;
the native relative-coordinate Jacobian squared is z^6, and the Schur
complement is D''. This proves (2), including the asymmetric s^-2
weight. It is not obtained by interchanging the old base variables.
The [relative critical-scheme chain rule](RELATIVE_CRITICAL_LOCUS.md)
requires K to be one nonzero constant after any allowed unit pullback.

At s=0, let rho be a simple zero of R0=R(0,p), v=R0'(rho) nonzero.
There is a literal closed critical line

\[
 p=\rho,\qquad
 r=-\frac{(2L\rho+M)(1+3vq)}{2Nv},\qquad q\text{ free}.
\]

Its full native Hessian is `4N^4 v^2`. One way to see independence
from all higher R derivatives is to use the local coordinate z=R(s,p).
At z=s=0 the bottom (q,r) Hessian block vanishes and the cross block
against (s,z) is `[[N,3H'(rho)],[0,2N]]`, up to transpose. Its squared
determinant is4N^4, and the native coordinate Jacobian contributes v^2.
The critical condition removes every nonlinear chain-rule correction.
Thus the line is an actual component, not a rational normalization.

For generic fixed s nonzero, let t count distinct roots of R(s,p) and
let r count distinct roots of F(s,p). Off both sets the value equation
solves q after choosing T, giving Euler1-t-r. At a root of R, the
original potential is `H(p)+4Ns^2*p+Ns*q`, so every value fibre has one affine
r line. At a root of F it has two affine lines, except at its critical
value D, when it has one. The roots of R and F are disjoint. Hence

\[
 \chi(f_s^{-1}(\lambda))=1+r-\#\{p_i:F(s,p_i)=0,D(s,p_i)=\lambda\}. \tag{3}
\]

In particular the generic Euler value is1+r>=2. The general
[HC4MYGR1 affine-component theorem](CRITICAL_COMPONENTS.md) now proves
the simple-boundary-root invertibility conclusion. If any active
p-dependent multiplicity-one F factor exists, it gives a separate
component away from s=0, so reducibility forbids unit Hessian.
An identically zero R0 is the earlier composite-trace obstruction;
a multiple boundary zero at p=-M/(2L) gives a critical two-plane.

## 2. The same global constraints, with the opposite Hessian weight (HC4MYGX1)

It remains to consider positive p degree and

\[
 R_0\ne0,\quad R_0\text{ has no simple zero},\quad
                         R_0(-M/(2L))\ne0.                     \tag{4}
\]

There are no boundary critical points. The finite-normalization proof
of [HC4MYGM0/1](CRITICAL_NORMALIZATION.md) applies with s replacing its
base: at every finite nonzero s, D is affine in p with nonzero coefficient
and the prefactor in (2) is a unit. Its conductor valuation excludes
all finite R_p zeros on active components, including singularities and
intersections. Thus the whole critical curve C=S=0 is closed, smooth,
etale over C*, and satisfies `(S,F_p)=(1)`. It is nonempty by the
empty-critical-scheme/Euler obstruction HC4MYGN0. A gradient inverse
would make C an affine line, incompatible with its nonconstant unit s.
Tang therefore supplies the entire rank-one plane Keller-fibre realization.

Set tau=D', E=D-s*tau, and let h be the generic degree of E+sigma*tau
on all C components. The scalar sigma is independent of the curve
coordinate s. Formula (3) gives the same two global identities as before:

\[
 h+\epsilon=r,\qquad
 \prod_{i\ {\rm inactive}}(Y-D_i(s))
      =\prod_q(Y-\alpha_q-\beta_qs)^{k_q}.                     \tag{5}
\]

Here epsilon sums the s-projection indices at all finite-base escapes;
the product on the right is over boundary points with E,tau finite.
The proofs are the whole-curve Euler comparison and the every-value
defect comparison of [HC4MYGW0/1](GLOBAL_CRITICAL_BUDGET.md) and
[HC4MYGY0](Y_GRAPH_CLOSURE.md), using (3). All inactive multiplicities
are counted once in r. The product's left side is squarefree, so k_q=1
and there is one distinct affine value for each inactive root. In particular

\[
 p_i=p_*(s)-\frac9{Ns^2}(\alpha_i+\beta_i s),\qquad
 p_*=-\frac{169N}{9L}s^2+\frac{13M}{L}-\frac{9M^2}{4LN}s^{-2}.  \tag{6}
\]

Every inactive root is rational with no finite nonzero pole.
The identity `dE=-s*d tau` still holds. A finite-value puncture at
s=infinity would have, for s=t^-e, `ord(tau-beta)=e+1`, because k_q=1.
Equation (2) would then have odd valuation `2*ord R_p+4e+1`, impossible.
Finite-base escapes have poles of E,tau, so all finite-value punctures
in (5) lie over s=0.

At any point over s=infinity, E has a pole of positive **even** order.
Indeed if tau has pole k then E has pole h=e+k; if tau is regular,
write `ord(tau-beta)=d`. Regular E would be the just-excluded finite-value
puncture; d=e would produce a logarithm, so d<e and E has pole h=e-d.
In both cases the valuation of (2) is `2*ord R_p+4e-h`, forcing h even.
It follows that every infinity point costs at least two in the global
pencil degree. This parity is specific to the s^-2 Hessian weight.

At a finite-base escape, p has pole m and s-s0 has order e. Formula (2)
forces m even and `ord R_p=e+m/2`, exactly as in HC4MYGM2 since its
prefactor is a unit there. The full pole budget is therefore

\[
 r=\sum_{\rm escapes}(m+2e)
          +\sum_{s=0,\infty}\max(0,-\operatorname{ord}E,-\operatorname{ord}\tau).
                                                                    \tag{7}
\]

No finite projection or sufficiency conclusion is inferred from these
necessary identities.

## 3. Close all finite critical projections (HC4MYGX2)

Suppose there are no finite-base escapes. Each connected component of C
is a finite etale cover of C*. Its Euler characteristic is zero. Its
projective normalization has genus g and at least two omitted points
(a zero and a pole of s); `2-2g-number_of_punctures=0` forces g=0 and
exactly two punctures. Thus it is C*, with a parameter t for which s=t^e.
This is a statement about this entire normalized curve, not about the
generic fibre of its residual plane polynomial.

Both p and R_p are Laurent polynomials in t, and R_p is a unit, hence
a monomial. Equation (2), applied term by term, gives

\[
 D=\alpha+\beta t^e+\gamma t^j,\quad \gamma\ne0,\quad
 j\notin\{0,e\},\quad j\text{ even}.
\]

Because C(s,p)=C(s,D)=C(t), `gcd(e,j)=1`; thus e is odd. A negative j
would make E,tau both finite at infinity, impossible by Section 2. If
j>e, the finite-value puncture at zero has degree j-e, so j=e+1.
The other possibility is0<j<e. The generic pencil degree contributed
by this component is respectively e+1 or e.

Use a common cyclic cover for all C components. Every active root p_i
is Laurent on it. Its difference from any other distinct F root has
no zero on C*: otherwise a finite active intersection or ramification
would give F_p=0. That difference is therefore a Laurent monomial.
This elementary unit fact controls all covering degrees simultaneously.

Let b count distinct inactive roots. The difference between an active
D=alpha+beta*s+gamma*s^(j/e) and an inactive affine D forces the two
affine parts equal: a sum with three distinct powers can be a Laurent
unit only if its constant and s coefficients vanish. Consequently b<=1.
Also all active exponents theta=j/e agree: two different nonaffine powers
cannot have monomial difference.

Write a for the total active degree, v for total inactive multiplicity
and n=a+v. The leading p coefficient of F is `c*s^d`. It has no root
on C*: every root p_i is finite there and a vertical content factor
would intersect the surjective active cover. If the leading coefficient
vanished while all roots stayed finite, every coefficient would vanish,
producing precisely such a content factor. At zero, F's constant p
coefficient is a unit, so d is the sum of all p-pole orders.
When all root differences have power theta-2, the root product and (2)
give `R_p proportional to s^(d-1+(n-1)*(theta-2))`, whereas (2)
requires exponent `2-theta/2`. Therefore

\[
             d=2n+1-(n-\tfrac12)\theta.                         \tag{8}
\]

If b=0, there are no finite-value punctures. Thus0<j<e on every
component, so0<theta<1 and e>=3. Comparing conjugates forces the same
affine parts on all components; all distinct root differences are indeed
constant multiples of s^(theta-2). Put `delta=-M^2/(4L)`.
If alpha differs from delta, every p root has pole two at zero, so d=2n
and (8) forces `theta=2/(2n-1)`. Its reduced denominator2n-1 exceeds
the component covering degree e<=n, impossible. If alpha=delta the
nonaffine term dominates at zero, so `d=n*(2-theta)`; (8) would force
theta=-2, also impossible.

If b=1, all affine parts equal the inactive alpha,beta. Exactly one
component has a finite-value puncture, hence j=e+1. Since theta agrees
on all components, there can be no other component with0<j<e, so there
is exactly one component and a=e. If alpha differs from delta, d=2n
again contradicts theta>1. If alpha=delta, beta must be nonzero since
the inactive root cannot be finite at s=0. Every root then has pole one,
so d=n. Equation (8) gives `e=2v-1`.

This last numerical profile must still pass the boundary. Put

\[
 z_0=9\beta/N,\quad
 Z=s(p-13M/L)+169N s^3/(9L)+z_0.
\]

The inactive factor is Z, and the one active factor is `Z^e-c^e*s`,
with c nonzero. No vertical content is allowed. Therefore

\[
 F=\frac3{z_0^{e+v}}(Z^e-c^e s)Z^v.
\]

Its boundary trace is

\[
 R(0,p)=\frac{3(e+v)}{2z_0}(p-13M/L)-\frac{3c^e}{2z_0^e},
\]

a nonconstant linear polynomial with a simple zero. This contradicts
the remaining boundary assumption (4). The exceptional cyclic profile
is thus excluded by the actual boundary, rather than incorrectly discarded
by exponent counting. All no-escape cases under (4) are closed.

## 4. The remaining distinct-root threshold (HC4MYGX3)

A nonautomorphic unit-Hessian positive-p opposite graph cannot have a
simple boundary root, zero trace or critical multiple root, by Section 1.
It must satisfy (4), and Section 3 requires at least one finite-base escape.
That escape costs m+2e>=4 in (7). Every active component has a point
over infinity, where the cost is even and at least two. Hence

\[
                                  r\ge6.
\]

Thus every positive-p opposite graph with at most five distinct generic
F roots has invertible gradient under actual unit Hessian, regardless of
the total degree caused by inactive multiplicities. This is a distinct-root
bound, not an import of a plane geometric-degree theorem. General cases
with r>=6 and finite-base escapes remain open. Their whole critical curve,
every-value factorization, even infinity poles and local conductor budgets
must all be retained before attempting another parameterization.

[verify_opposite_nonlinear.py](verify_opposite_nonlinear.py) and its
[receipt](receipt-opposite-nonlinear.json) check the asymmetric critical
formula, the extra zero-R fibre stratum, boundary cross block, Laurent
classification identities, and the exceptional cyclic profile's actual
simple boundary. The all-degree covering and Laurent-unit arguments are
written proofs. Inputs are retained Tang, local normalization and the
classical line statements already used in HC4MYGR1. No new generic-rational
fibre theorem, formal verification, independent external review or novelty
priority is claimed. The two distinct plane endpoints, q/r-dependent bases,
actual free-base/dual mixing and unrestricted HC4/JC2 remain open.
