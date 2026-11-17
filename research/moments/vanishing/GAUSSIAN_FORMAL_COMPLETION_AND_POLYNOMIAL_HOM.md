# Exact Gaussian completion and polynomial Hom

This algebraic proof works over the stated coefficient rings. It is an internal written argument; independent review and formal verification remain outstanding.

## 1. The exact source has a rank-one formal completion

Let A be any commutative Q-algebra and P in A[X_1,...,X_n]. Write

    V=A[t,X],  nabla_i=partial_i-X_i+t*P_i,
    H_P=V/sum_i nabla_i V,  partial_t[A]=[partial_t A+P A].

Here H_P is the **exact** Gaussian integration quotient, not an arbitrary
annihilating presentation for one of its periods. The nabla_i commute with
the displayed time connection, so H_P is a left Weyl module over A. The
standard formal Gaussian functional is

    Phi_P([B])=E[B(t,X) exp(tP(X))] in A[[t]].                (1)

Integration by parts makes (1) well-defined and D_t-linear. E denotes the
normalized product Gaussian Wick functional, defined over every Q-algebra.

**Truncated completion theorem.** For every k>=1, the map induced by (1) is
an isomorphism of A[t]/(t^k)-modules

    H_P/t^k H_P  ->  A[t]/(t^k).                            (2)

To prove it, work in R=(A[t]/t^k)[X]. Multiplication by the finite exponential

    U=exp(tP) mod t^k

is invertible, with inverse exp(-tP) mod t^k. The identity

    (partial_i-X_i)(U B)=U*nabla_i B                       (3)

identifies the quotient in (2) with the unperturbed Gaussian quotient.
The latter is A[t]/t^k via E. In one variable,

    X^j = (j-1)X^(j-2) - (partial_X-X)X^(j-1),

and the case j=1 starts the reduction. Iterating over variables gives both
surjectivity and the kernel equality with sum_i(partial_i-X_i)R; E kills
that sum and E(1)=1. This works over A without reducedness or domain hypotheses.
The image of [1] in (2) is E exp(tP), a unit modulo t^k. Thus [1] is also a
basis of the quotient, although its image under the normalized period map
need not be the constant1.

Consequently the ordinary t-adic completion is

    lim_k H_P/t^k H_P = A[[t]].                            (4)

These quotients in (2) are structure-sheaf modules, **not D_t quotients**:
partial_t does not preserve t^k H_P. The time action is continuous with a
one-step filtration loss and extends to the completed module; under (4) it
is ordinary differentiation. No time-infinity regularity assumption occurs.

There is a second useful form of the theorem:

    Hom_Dt(H_P,A[[t]]) = A*Phi_P.                           (5)

Indeed any such map psi gives a functional C(Q)=psi([Q])(0) on A[X].
Evaluating the integration-by-parts relations at t=0 shows
C(partial_i Q-X_i Q)=0, so C(Q)=a E(Q), where a=C(1).
For t-independent Q, partial_t^m[Q]=[Q P^m]. Therefore

    psi([Q])=a sum_m E(Q P^m)t^m/m!.

This proves uniqueness, including when A has zero divisors. Equivalently,
a D_t-linear map is t-adically continuous and (2) identifies its continuous
extension with a differential endomorphism of A[[t]], necessarily a scalar.

## 2. The actual period kernel is the infinite-divisibility part

Equations (1)--(2) identify the entire actual kernel intrinsically:

    K_P=ker Phi_P=intersection_(k>=0) t^k H_P.              (6)

It is a D_t-submodule, since differentiation lowers divisibility by at most
one. Moreover t*K_P=K_P. If h in K_P, write h=t*h_1 using (6). Then
t*Phi_P(h_1)=Phi_P(h)=0, and multiplication by t is injective on A[[t]],
so h_1 also belongs to K_P. Conversely every submodule N with tN=N maps
to zero: its image lies in every t^k A[[t]]. Thus K_P is the largest
t-surjective submodule of H_P. Multiplication by t need not be injective
on K_P; punctual modules at the origin must not be discarded.

The actual image is H_P/K_P and is t-adically separated. This is a structural
description, not an implemented algorithm for the infinite intersection.
Computing this infinite intersection is a separate algorithmic problem.

Rank one in (4) does **not** mean that H_P or its actual generic connection
has rank one. The source is finite over D_t over a field, not necessarily
finite over A[t], so a Krull-intersection argument cannot force (6) to vanish.
For an abstract control, O_t direct_sum K[t,t^-1] has completion O_t^,
while its second summand is nonzero and t-invertible. The completion statement alone places no bound on actual generic rank.

## 3. Safety is existence of an abstract polynomial quotient

Now take A=K, a characteristic-zero field, and O=K[t]. The exact source H_P
is holonomic by the retained algebraic integration construction. Since
O embeds in K[[t]], (5) gives

    dim_K Hom_Dt(H_P,O) <= 1.                              (7)

If a nonzero map exists, its formal extension is a*Phi_P with a!=0.
Dividing by a proves that **every** original multiplier period is polynomial.
The converse is immediate. Hence

    P is GMC-safe  iff  Hom_Dt(H_P,O) != 0.                 (8)

No pure-null assumption is needed for this eventual-mixed statement. At a
pure point the polynomial quotient is normalized by [1] -> 1. A nonzero
submodule of O is all of O: differentiate a nonzero polynomial to a nonzero
constant. Thus the maximal quotient of H_P that is a finite direct sum of
copies of O is zero in the unsafe case and a **single copy of O** in the safe
case. It equals the actual image precisely in the safe case.

This removes the need to compute the actual kernel or compare a selected
period with a polynomial solution basis for the pointwise safety decision.
Compute an exact finite presentation of H_P, then polynomial Hom. Over an
exactly presented effective field this is a terminating theoretical method.
[Tsai--Walther,math/0007139v1](https://arxiv.org/html/math/0007139v1),
Sections1--3, especially Algorithm3.3, supplies explicit polynomial Hom bases
for finitely presented holonomic modules. The [whole-source construction](WHOLE_GAUSSIAN_SOURCE_HOM_AND_EFFECTIVE_CUTOFFS.md)
records the current finite-input implementation and its limits.

Exactness of the integration presentation is essential in (7)--(8). For
the Gaussian constant phase P=1, the true source is D_t/D_t(partial_t-1),
and the period exp(t) is not polynomial. The larger scalar presentation
D_t/D_t[partial_t(partial_t-1)] also annihilates exp(t) but admits the
polynomial solution1. Its extra polynomial quotient does not certify safety.
Criterion (8) requires the exact integration source and uniqueness statement (5).

## 4. Coefficient rings retain a normalization ideal

Over a general A, evaluation psi -> psi([1])(0) embeds
Hom_Dt(H_P,A[t]) into A. Its image is the ideal

    J_P={a in A : a*Phi_P(H_P) is contained in A[t]}.        (9)

Safety is equivalent to 1 in J_P, not merely J_P!=0. An arbitrary nonzero
scalar need not be invertible. There is no assertion that J_P is the ordinary
pure-moment ideal or one of the prescribed mixed cutoff ideals.

The complete one-variable linear Gaussian family P=aX has exact presentation

    H_(aX)=A[t]e,  partial_t e=a^2*t*e.                     (10)

The relation X=at in the Gaussian quotient and its monomial recurrence prove
(10), with no localization in a. Its normalization ideal is

    J_(aX)=union_(N>=0) Ann_A(a^N).                         (11)

The formal image of e is exp(a^2t^2/2); polynomiality of a scalar multiple is
equivalent to eventual annihilation of its even powers of a, which is
equivalent to annihilation by some power. Conversely this also makes all
marked moments polynomial, since every mixed m-th moment has the factor a^m.
For A=Q[a,b]/(ab), the ideal in (11) is (b), which is nonzero while the
normalized period is unsafe. This is an exact coefficient-ring control,
not a counterexample over a field.
