# A quadratic critical-face obstruction with negative Laurent tails

**QCL1.** A weighted quantum rescaling extends the finite-mate obstruction
to negative Laurent coefficients, including the first inverse-square
boundary potential. Higher generic momentum orders are allowed above
the critical face. The argument uses a polynomial centralizer at that
face to remove the finitely many parameter poles of a hypothetical mate.

## The Laurent support criterion

Let K be a characteristic-zero field, R=K[[t]], L=K((t)), and
eta=t^e*u(t) with e>=1 and u in R a unit. In the finite scalar
Ore algebra
\[
             L(v)[a;\eta\delta],\qquad
                    \delta=v^{-1}\partial_v,
\]
suppose S is a finite Laurent polynomial in v and a, with all scalar
coefficients in R and classical reduction
\[
                         S\bmod t=g_0a^2+\mu_0v,
                         \qquad g_0\mu_0\ne0.               \tag{1}
\]
The exponents of a are nonnegative; the exponents of v may be any
integers. For a nonzero normal coefficient c_(j,m)(t) of v^m*a^j,
put r_(j,m)=ord_t c_(j,m) and
\[
                  F(j,m)=5r_{j,m}+e(j+2m-2).                \tag{2}
\]
The two anchor coefficients of a^2 and v have r=0 and F=0.
Assume:

1. Every nonzero coefficient has F>=0.
2. Among the terms with F=0, the maximum a order is two.

Then no finite scalar T over L(v), or any finite algebraic coefficient
cover, satisfies [S,T]=kappa for a nonzero constant kappa in L.
No generic upper bound on ord_a S or ord_a T is required. The second
condition concerns only the critical face, not the complete operator.

In particular, if every correction to the two anchors has F>0, the
conclusion holds. More strongly, the terms with F=0 may include
t^2*v^-4 in an ordinary quantum parameter eta=t, and arbitrary other
critical terms of a order at most two. The inverse-square coefficient
is unrestricted. Neither a one-coefficient p-power cancellation nor a
larger generic momentum order above this face supplies a finite mate.

## A finite mate cannot hide behind parameter poles

First reduce an algebraic-cover mate to L(v) by the normalized
coefficient-field trace of [QPC1](UNIVERSAL_P_POWER_OBSTRUCTION.md).
Scale it by eta/kappa, so that [S,T]=eta. Introduce s and perform the
finite changes
\[
 t=s^5,\qquad v=s^{2e}V,\qquad a=s^e A,\qquad
 \mathcal S=s^{-2e}S,\qquad \mathcal T=s^{-3e}T.             \tag{3}
\]
These are identities over a generic scalar field, not a specialization
of an unbounded formal operator series. They give
\[
 [A,V]=u(s^5)/V,\qquad [\mathcal S,\mathcal T]=u(s^5).       \tag{4}
\]
A coefficient of S with order r produces in mathcal S the s order
5r+e(j+2m-2), exactly F. Hence mathcal S is integral and has a
finite reduction P at s=0. The order-two face hypothesis gives
\[
 P=G(V)A^2+H(V)A+J(V),\qquad
 G(V)\longrightarrow g_0,\quad H(V)=O(1),\quad
                  J(V)=\mu_0V+O(1)                         \tag{5}
\]
at V=infinity. Indeed every new critical term has r>=1. For j=2,
equation F=0 forces m<0; for j=1 or j=0 it forces m<=0.
Thus no new coefficient changes the displayed leading infinity data.

Here is the exact coefficient ring in which pole removal is valid.
On K((s))(V) use the s-Gauss valuation: the valuation of a rational
function is the minimum scalar s valuation of its numerator
coefficients minus that of its denominator coefficients. Its valuation
ring O has residue field K(V), and delta_V=V^-1*partial_V preserves O.
For example a denominator whose Gauss valuation is zero has nonzero
reduction, so its inverse is integral and its derivative remains
integral. The finite Ore algebra O[A;u(s^5)delta_V] is therefore
closed under multiplication and reduction. Every coefficient of
mathcal T has a finite s pole; there are only finitely many such
coefficients, even for a rational-position mate.

We prove below that
\[
             \operatorname{Cent}_{K(V)[A;u(0)\delta_V]}(P)
                                      =K[P],                \tag{6}
\]
and that P has no finite nonzero constant mate. Suppose mathcal T
has pole order N>0. Equation (4) shows that its leading coefficient
\(\overline{s^N\mathcal T}\) commutes with P. By (6) it equals q(P)
for a finite polynomial q over K. Replace
\[
              \mathcal T\longmapsto
                       \mathcal T-s^{-N}q(\mathcal S).       \tag{7}
\]
This is an actual finite operator, retains (4), and strictly lowers
the pole order. After at most N repetitions mathcal T is integral.
Its reduction is then a finite rational-coefficient mate of P with
commutator u(0) nonzero, a contradiction. No infinite polynomial in
mathcal S, generic order bound, or bounded formal gauge is used.

## Polynomial centralizer and mate exclusion at the quadratic face

Extend constants algebraically and adjoin a square root rho of G.
Set E=rho*A, whose derivation is
epsilon=u(0)*rho*delta_V. A coefficient-preserving shift of E removes
the linear coefficient in P, giving
\[
                 P=E_1^2+U,
                 \epsilon\sim u(0)\sqrt{g_0}\,V^{-1}\partial_V,
                 U\sim\mu_0 V.                            \tag{8}
\]
One may see the last asymptotic directly: with
ell=H/rho-u(0)*delta_V(rho), write E_1=E+ell/2 and
\[
                  U=J-ell^2/4-\epsilon(ell)/2.
\]
Here rho tends to a nonzero constant and ell is bounded at infinity,
whereas J has the nonzero linear growth in (5). Consequently
epsilon(U) is nonzero and tends to zero: its leading term is a
nonzero multiple of V^-1.

The universal quadratic scalar-commutator normal form in
[QSO2, equation (8)](QUADRATIC_OPERATOR_MATE_CLASSIFICATION.md)
also applies to a commuting operator: after even-power subtraction,
its remaining odd part is a finite sum of the operators B_(2j+1)
with scalar constant coefficients. Commutation forces
\[
                   \sum_{j=0}^{q}c_j R_{j+1}(U)
                                  =\hbox{constant}.         \tag{9}
\]
The same universal Lenard polynomials have leading pure coefficient
beta_l=binom(2l,l)/(2*4^l), nonzero in characteristic zero.
Since U has V degree one and epsilon lowers V degree by two,
\[
                    R_l(U)\sim\beta_l\mu_0^l V^l.           \tag{10}
\]
To check dominance in every l, a weight-2l differential monomial
with q factors of U and total derivative order d has 2q+d=2l
and V degree at most q-2d. This equals l only for q=l,d=0;
every derivative monomial has smaller degree. The largest j with
c_j nonzero in (9) therefore cannot cancel. No odd part survives,
so the centralizer over the enlarged coefficient curve is polynomial
in P. Descend to K(V)[A] by comparing highest A coefficients in a
polynomial in P: its scalar coefficients lie in K(V) and are
constants, hence lie in K. This proves (6).

For a nonzero constant mate, QSO2's all-order theorem says that in
the normalized quadratic chart epsilon(U) must be a nonzero
constant. Its nonzero decaying asymptotic from (8) contradicts this.
This proves the second assertion needed in the pole-removal argument.
No classification of arbitrary higher-order centralizers is imported.

## A canonical critical normalization for Long's actual finite images

An actual finite exact-R image has a finite Laurent normal form by
[QEP1's equality of localized quantum algebras](EXACT_POSITION_CENTRALIZER_AND_COMPLETION.md).
Take the constant field K to include the rational exact-position
parameter r, and write its v^0 part as a finite polynomial G_t(a).
For Long,
\[
                   G_0(a)=4a-3ra^2,
\]
whose derivative has the simple root alpha_0=2/(3r).
There is a unique alpha(t) in K[[t]], reducing to alpha_0, with
G_t'(alpha(t))=0. This follows by coefficient induction: at every
new t order the unknown coefficient is multiplied by
G_0''(alpha_0)=-6r, a unit in K. Replace a by a+alpha(t) and
subtract the scalar G_t(alpha(t)). This is an actual finite scalar
Ore change; delta annihilates alpha. It preserves existence of a
mate, Laurent v support, and finiteness of momentum order.

The normalized image has precisely the classical form (1), with
g_0=-3r and mu_0=2. Its v^0 linear and constant coefficients vanish
identically, which avoids a spurious obstruction from a removable
scalar shift. The shift may live in the enlarged scalar field
K((t)), with K already containing the rational parameter r; it is
being used for a necessary nonexistence test, not
as a bounded target Weyl automorphism preserving global admission.

It follows that every actual finite scalar mate for an exact-R Long
image must have, in this critical normalization, at least one of:

- a coefficient with F<0;
- a critical coefficient F=0 of momentum order at least three.

There is a broad immediate corollary: **no regular finite polynomial
deformation in v,a of Long's classical S has a finite scalar mate**,
even allowing arbitrary rational-position or finite algebraic-cover
coefficients in that mate. There is no restriction on the positive
v support, mixed coefficient shape, or generic momentum order in
this corollary. Indeed the critical scalar shift preserves nonnegative
v powers and removes all v^0 constant and linear terms. Every
remaining correction has r>=1 and either m=0,j>=2 or m>=1,j>=0.
In both cases F>0. Its critical face is just g_0*A^2+mu_0*V,
so QCL1 applies. A bounded exact-R scalar proposal therefore needs
an actual negative-v tail, rather than only a different polynomial
mixed coefficient face. Global admission may permit such a tail;
this corollary does not exclude it.

These are checkable necessary conditions on an individual finite
candidate, not a claim that such a candidate exists. With eta=t,
a term t^r*v^m*a^j can escape the quadratic-face test only on or
below
\[
                         5r+j+2m=2.                        \tag{11}
\]
All possible escapes have negative v powers. For example the
first r=2 scalar equality term is t^2*v^-4; QCL1 still excludes
it if the full critical face stays quadratic. An even t branch
has only even r, so (11) forces an even j. Its first possible
nonquadratic critical order is therefore at least four, including
t^2*v^-6*a^4. A more polar coefficient can instead have F<0.
Generic higher momentum terms strictly above the face do not help.

[Exact controls](check_critical_laurent_face.py), with
[the retained receipt](critical_laurent_face_controls.json), verify the weighted
rescaling, a higher-order term above a quadratic Laurent face, the
finite quadratic normalization and infinity asymptotics, several
Lenard leading terms, and finite pole removal in a positive Airy
control. They also retain positive first-order and Airy pairs to
mark the quadratic and infinity hypotheses. Universal centralizer
exhaustion and mate exclusion depend on QSO2's written proof, and
the pole-removal and critical-root inductions are written arguments.

Negative-weight faces, higher-order critical faces, nonpositional
first images and genuinely nonscalar carriers remain open. No
bounded full Long lift, DC1/DC2 solution, formal verification,
independent whole-proof replay, external review or novelty
priority is claimed.
