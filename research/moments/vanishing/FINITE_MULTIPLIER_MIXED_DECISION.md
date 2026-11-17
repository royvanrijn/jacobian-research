# A finite multiplier criterion and pointwise decision for eventual Gaussian moments

15 November 2026. Internal all-degree, all-dimension deduction from classical
algebraic holonomic integration. No novelty, independent replay, external
review or formal verification is claimed. The general module-generation
algorithm is not implemented or executed. The
[checker](verify_finite_multiplier_decision.py) verifies the scalar
polynomiality test on supplied operators, a complete elementary safe
module, and exact source/low-order identities for the retained quartic
failure. It does not certify the general integration theorem.

## 1. What becomes decidable

Let K be a specified effective finitely generated characteristic-zero field,
with an exact presentation over Q. Examples include number fields and
fraction fields of specified prime quotients of rational polynomial rings.
Let P in K[X_1,...,X_n], in arbitrary dimension and degree. All Gaussian
moments below are formal Wick moments. Write

    F_Q(t)=E(Q(X) exp(tP(X))).

**Finite multiplier theorem.** A computable integer b(P)>=0 has the property

    every Q has E(QP^m)=0 for all sufficiently large m
      iff F_(X^alpha)(t) is a polynomial for every |alpha|<=b(P).   (1)

Consequently there is a terminating exact pointwise decision algorithm for
the eventual mixed property. If it fails, the algorithm produces a monomial
multiplier of degree at most b(P), a scalar differential certificate and a
finite coefficient test proving its series is nonpolynomial, hence has
infinitely many nonzero mixed moments. It can also find the least possible
degree of a failure multiplier.

All-order pure nullity is separately decidable by the retained
[holonomic zero oracle](GENERIC_POINT_FINITE_MOMENT_ALGORITHM.md).
Thus one can decide whether a supplied exactly represented polynomial is a
pure-null GMC counterexample. This is an algorithmic criterion, not a
geometric classification of all safe or failing coefficient families.
For arbitrary complex multipliers, the same finite test suffices after
extending K to C: each multiplier is a finite linear combination of monomials.

No degree-independent bound, efficient complexity estimate, numerical value
of b(P) in a new large example, uniform linear mixed cutoff, or general
solver implementation is asserted. Arbitrary black-box complex coefficients
without decidable exact relations are outside the algorithm's input model.

## 2. Use the whole integration module over the polynomial Weyl algebra

Put D_t=K[t]<partial_t>, and define

    V=K[t,X],
    nabla_i A=partial_(X_i) A+(t P_i-X_i) A,
    nabla_t A=partial_t A+P A,
    H=V / sum_i nabla_i V.                                    (2)

The operators nabla_i commute with nabla_t. Hence H is a left D_t-module.
The underlying exponential connection is holonomic; algebraic integration
makes H holonomic and, in particular, finitely generated over D_t.
The classical integration algorithm computes polynomial representatives
A_1(t,X),...,A_r(t,X) of a generating set and their relations.

More explicitly, the exponential connection is cyclic over the full
polynomial Weyl algebra, with relations partial_t-P and
partial_(X_i)-(tP_i-X_i). Normal ordering identifies its quotient with V.
In [Oaku, section 3, the algebraic integration-module construction preceding
definition 4](https://www.lab.twcu.ac.jp/oaku/Oaku2013.pdf), a computable
integration b-function supplies an integer bound b such that the classes
of X^alpha, |alpha|<=b, generate H over D_t. The computation also gives the
whole finite module presentation. This is the classical imported input,
not a new Gaussian theorem of that paper. A negative generator bound would
make H zero, which is impossible here because the functional below takes
[1] to a series with constant coefficient one.

Retain this full presentation. Eliminating down to the integration ideal
annihilating [1] loses the other amplitude classes. Likewise a basis after
localizing to K(t) is insufficient: division by a polynomial in t can turn
a polynomial into a rational function with poles.

Define the formal Gaussian map

    Lambda:H -> K[[t]],       [A] -> E(A(t,X) exp(tP)).         (3)

It is well-defined by coefficientwise Gaussian integration by parts, and
is D_t-linear:

    Lambda(nabla_t A)=partial_t Lambda(A),
    Lambda(t A)=t Lambda(A).

No exponential convergence or analytic boundary argument is required.
The subspace K[t] is stable under both t and partial_t. Therefore all of
Lambda(H) lies in K[t] exactly when the finitely many generator images do.
Every Q is an amplitude in H, and every polynomial t-dependent amplitude
is a finite sum of t-powers times fixed Q. This proves (1), including both
directions of the quantifier over all multipliers.

The argument uses a polynomial Weyl-algebra generating set, not a
finite-dimensional Jacobi ring and not a rational-function basis.

## 3. Deciding whether one specified holonomic series is polynomial

Suppose a series f in K[[t]] has computable coefficients and a nonzero
scalar annihilator

    L=sum_(j=0)^r p_j(t) partial_t^j,       p_j in K[t].        (4)

All F_Q have such an operator with a polynomial Gaussian total-derivative
certificate. They can be computed using holonomic module operations, or
by the finite-degree-box search from the earlier zero oracle, now with
seed Q instead of 1.

For nonzero coefficients put

    h=max_j(deg_t p_j-j),
    C(z)=sum_(deg p_j-j=h) lc_infinity(p_j) (z)_j.              (5)

C is nonzero by independence of falling factorials. If a nonzero polynomial
f has degree N, the coefficient of t^(N+h) in Lf is lc(f) C(N).
Thus its degree must be a nonnegative integer root of C. Compute these
finitely many roots in K. Integer-root extraction over a presented finitely
generated field reduces, after clearing denominators and prime normal
forms, to a gcd of finitely many polynomials in Q[z].

Let B be the largest such root; if there are none put B=-1. Form the
candidate T=sum_(k=0)^B [t^k]f*t^k, with T=0 for B=-1.
It remains essential to decide f-T=0, rather than merely inspect this
truncation. Put S=L(T). A nonzero scalar annihilator of f-T is

    M=L                         if S=0,
    M=(S partial_t-S') composed L otherwise.                  (6)

Indeed L(f-T)=-S, and the left factor kills S. M is nonzero because the
scalar differential-operator ring is a domain. Apply the origin indicial
zero test to f-T: for M=sum q_j partial_t^j set

    mu=min_j(ord_0 q_j-j),
    B_0(z)=sum_(ord_0 q_j-j=mu) lc_0(q_j) (z)_j.               (7)

Every possible first nonzero coefficient index is a nonnegative integer
root of the nonzero B_0. A finite exact jet therefore decides f=T.
If equality holds, f is polynomial and its exact degree is available.
Otherwise (5) excludes every polynomial solution of higher degree, so f
is nonpolynomial. A detected coefficient above B, together with (4)--(7),
is a finite certificate for infinitely many nonzero coefficients.

This reasoning allows singular or irregular points; it is not an
ordinary-point initial-value argument. For example f=t^5/(1-t) has
infinity degree bound four for the valid operator

    t(1-t) partial_t-5+4t.

Its truncation through degree four is zero, but its origin indicial root
five detects the missing term. Inspecting only the infinity-bound jet
would falsely declare f polynomial.

The two endpoint arguments are classical in scalar differential algebra.
[Abramov--Ryabenko, Programming 2025, no.1, pp.51--58,
DOI 10.31857/S0132347425010066](https://journals.rcsi.science/0132-3474/article/download/287079/264547),
sections 2 and 4, discusses degree/valuation indicial bounds. Its general
operator-composition theorem and difference/q-difference cases are not
needed here. The direct differential proof above is included in full.

## 4. The terminating pointwise algorithm and exact outputs

1. Compute the full polynomial integration-module generators from section 2.
2. For each monomial generator, apply the scalar polynomiality decision.
3. Accept the eventual mixed property exactly when every test accepts.
4. If desired, separately test F_1=1 to determine pure nullity.

Every computation terminates by the cited holonomic algorithms, finite
integer-root computations and finite scalar jet tests. If a general
integration presentation initially uses t-dependent amplitude generators,
expand their finitely many t coefficients into fixed polynomial amplitudes;
testing all of those is sufficient as well. The monomial presentation gives
the explicit degree bound in (1).

If failure is found, test monomials in increasing degree up to b(P) to
find the least witness degree. Linear combinations cannot introduce a
nonpolynomial series when all basis monomials give polynomials.

If P passes, for any requested d one can compute the exact last possible
nonzero mixed index among all multipliers of degree at most d: compute
the polynomial F_(X^alpha) for |alpha|<=d and take their largest degree.
A monomial attaining that maximum proves sharpness at that d. This is a
computable cutoff function of d; no linear-growth theorem is used or
claimed. Generator representations also give certified, possibly coarser,
bounds for an individual Q.

There is a module-theoretic formulation. The nonzero image Lambda(H)
is a D_t-submodule of K[[t]]. The eventual mixed property holds exactly
when

    Lambda(H)=K[t].                                           (8)

If the image lies in K[t], it equals K[t], because any nonzero polynomial
generates K[t] under derivatives and multiplication by t. This is the
actual Gaussian image, not the mere existence of some abstract trivial
quotient of H. Generic rank one or rationality of all image series does
not suffice, as the next example shows.

## 5. The retained quartic failure has a rank-one rational image with a pole defect

Use circular coordinates E(ZW)=E(T^2)=1 and the retained polynomial

    B=W(1+Z)-(1+3Z/2+Z^2/2) T^2.                              (9)

The [affine-W classification](CUBIC_AFFINE_W_CLASSIFICATION.md) already
proves F_1=1 and F_Z=t/(1-t). The following describes its entire Gaussian
image without a new coefficient search.

Introduce sources alpha Z+beta W+gamma T and set

    q=1-t,
    g=(t+beta)/q,
    R=1+t(3-t)beta+t beta^2.

The retained Gaussian--Lagrange identity and quadratic T integration give

    E exp(tB+alpha Z+beta W+gamma T)
      =R^(-1/2) exp(alpha*(t+beta)/q+gamma^2*q^2/(2R)).        (10)

To check the prefactor, the fixed-point Jacobian is q and the quadratic
T factor is K=1+t(2+3g+g^2)=R/q^2; their normalized square roots cancel.
All identities are formal at zero, with the square root of constant term
one. At any fixed source degree, R^(-1) and R^(-1/2) contribute polynomials
in t; the only possible denominator is a power of q from the alpha factor.
Consequently every polynomial multiplier has

    F_Q in K[t,(1-t)^(-1)].

The image contains 1 and t/(1-t), hence also 1/(1-t). Derivatives of the
last function give every higher pole power, and multiplication by t gives
the polynomial part. Therefore

    Lambda(H_B)=K[t,(1-t)^(-1)].                              (11)

This image has generic rank one over K(t). Nonetheless it contains
nonpolynomial series, and B fails GMC. More precisely,

    0 -> K[t] -> Lambda(H_B) -> delta_(t=1) -> 0              (12)

is exact as a sequence of polynomial Weyl-algebra modules. The quotient is
generated by the class of (1-t)^(-1), killed by t-1; its derivatives give
all principal parts at 1, so it is the usual simple point module delta_1.
The sequence does not split: the middle term lies in the field K(t), so
has no nonzero element killed by t-1.

This identifies a structural obstruction invisible to generic rank:
a point-supported quotient measuring poles at a finite parameter value.
It does not say that every failure has this form. In particular, rational
trivial quotients and generic seed ranks in the remaining normalized
constant-b classification still require their own geometric analysis.

For comparison P=W has full integration module K[t] generated by [1].
The total-derivative relations reduce W A to partial_Z A and Z A to
t A+partial_W A. They reduce every monomial class to a multiple of [1],
and partial_t[1]=[W]=0. The Gaussian map sends [1] to 1, proving this full
module is exactly K[t]. Explicitly,

    F_(Z^a W^b)=a!/(a-b)! t^(a-b) if a>=b, and zero otherwise.

Thus a one-element module test controls all multipliers in this safe case.

## 6. Specialization, literature and remaining work

At each exactly represented P the full generators and scalar operators
must be recomputed or specialized with a proved valid module certificate.
Neither generic generator bounds nor generic indicial data may be reused
blindly at exceptional points. The ordinary-ideal algorithm from the
[previous turn](PRIMARY_THICKENING_MOMENT_IDEAL_ALGORITHM.md) concerns
all-order zero identities and does not itself decide eventual tails.

General holonomic parameter families illustrate the danger: the equation

    (1+t) f'-a f=0,       f(0)=1

has f=(1+t)^a. Over Q(a) there is no nonnegative integral infinity root;
at a=N>=0 the series is polynomial of degree N. Such degrees are unbounded.
This is a holonomic example, not a claim that this family arises as
Gaussian moments or that Gaussian safety loci are nonconstructible.

One implication does descend correctly: if the generic field of a rational
prime coefficient component passes the full eventual mixed test, then the
whole component passes. For each fixed monomial Q its generic polynomial
series says all sufficiently high polynomial moment coefficients lie in
that prime. Those identities vanish at every complex point, including
certificate denominator exceptions. A generic failure does not permit
cutting away every point where one finite mixed moment is nonzero: such a
point might have a later finite cutoff. No finite stratification or
coefficient-locus constructibility theorem is claimed.

Selected recent moment searches were reconciled with the retained Gaussian
counterexample and binary theorem; no new unrestricted positive GMC theorem
is imported. The current [Long source](sources/long_gaussian_de4082c.tex),
retrieved through the GitHub connector at blob
de4082caa41f537c9712a38b934a83172a6658d9, was reread for the introduction,
ternary theorem/proof and the half-pair generating identity. Its older
sentence leaving dimension two open is superseded by the retained binary
theorem; it is not the current local open list. The full source identity
(10) and image extension (11)--(12) above are internal deductions from
the already retained Gaussian--Lagrange mechanism.
The new reusable foundation is the full classical integration
module, rather than an additional low-degree calculation. The source audit
records the 2025 indicial paper, its actual sections read, and the full
integration-module passage. No producer package was executed.

This closes theoretical pointwise decidability of the eventual mixed
property, a computable failure-multiplier degree bound, and the rank-one
rational-image obstruction (11)--(12). It does not implement a general
solver, classify the exceptional coefficient loci, compute a new cubic
number, prove uniform linear cutoffs or remove independent-review
obligations. Structural normal forms, higher-rank quotients and valid
operator-order transfers remain open.
[Source audit](finite_multiplier_source_audit.json).

