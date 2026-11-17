# Gaussian logarithm exclusion and rational-primitive rigidity

17 November 2026. The user proposed the pure unit together with the
original marked period -log(1-t), and required ALL original markers to
stay regular. The candidate is excluded, even without purity or any
regularity premise on the other markers. A stronger consequence excludes
rational primitives and embedded unipotent extensions in the whole actual
generic image. General regular transcendental realization remains OPEN.
These are internal written proofs; independent replay, formal verification
and external review remain false; priority unassessed.

## 1. A Gaussian Frobenius congruence with all original moments retained

First take diagonal covariance diag(c_1,...,c_n), with c_i in Z_(p),
where p is an ODD prime. Work over any commutative Z_(p)-algebra A.
Gaussian expectation is A-linear, with

    E[X^alpha]=product_i c_i^(alpha_i/2) (alpha_i-1)!!

if every alpha_i is even, and zero otherwise. If alpha_i>=p for any i,
this moment is zero modulo p: an odd exponent gives zero, while an even
one is at least p+1 and its odd double factorial contains p. Consequently

    E[(X_1^p,...,X_n^p) A[X]] subset pA.                (1)

For ANY polynomials F,Q in A[X], the freshman's-dream identity modulo p
and (1) give

    E[Q F^p] = F(0)^p E[Q] modulo pA.                 (2)

There is no degree restriction on F or Q. In particular, for any phase P
and EVERY integer m>=0, the actual original moments obey

    E[Q P^(m+p)] = P(0)^p E[Q P^m] modulo pA.         (3)

These are RAW moments; no factorial denominator is reduced modulo p.
For standard covariance these statements hold at EVERY odd p whenever
the finitely many input coefficients are p-integral. At p=2 they fail:
E[X^2]=1, whereas X(0)^2 E[1]=0.

For a fixed nonsingular rational symmetric covariance C, choose a rational
congruence diagonalization C=U diag(c_i) U^T. Apply (2)-(3) to P(UY),
Q(UY); the origin is unchanged. Exclude the FINITE set of primes in
the denominators of this fixed change and of c_i. Thus the same assertion
holds after one finite arithmetic localization for C. No real positivity,
convergence, purity, actual-rank or image-kernel hypothesis is needed.
The diagonalization is over Q, so no constant-field extension is required.

## 2. The proposed original logarithmic marker is impossible

Use the user's exponential normalization

    Phi_Q(t)=sum_(m>=0) E[Q P^m] t^m/m!.

The proposed identity Phi_Q=-log(1-t) requires

    E[Q]=0, E[QP^m]=(m-1)! for every m>=1.            (4)

Over rational or number-field coefficients, choose any sufficiently large
good odd prime and a residue field above it. At m=p, (2) and E[Q]=0
give E[QP^p]=0 modulo that field. But (p-1)! is a UNIT there (indeed
-1 in the prime field). This contradicts (4).

Here is the same proof over ANY characteristic-zero field K, including
transcendental coefficients. All coefficients of P,Q and of the fixed
covariance change lie in ONE finitely generated Z-domain A inside K,
after a finite localization. All moment identities (4) hold in A, since
raw Gaussian moments have integral polynomial expressions in those
coefficients. Spec A has positive-characteristic fibers at arbitrarily
large primes: its dominant finite-type image in Spec Z contains an open
set. Choose any such fiber of odd characteristic p outside the fixed bad
set. Equation (2) gives zero, whereas (p-1)! remains nonzero. This is a
single exact fiber contradiction, not specialization of an inferred rank.

Therefore there is NO polynomial Gaussian phase P and ORIGINAL polynomial
marker Q over ANY characteristic-zero field with Phi_Q=-log(1-t).
This holds whether P is pure or impure and whether the other original
markers are regular, irregular, convergent or divergent.

It also excludes the literal identity over EVERY NONZERO Q-algebra:
reduce through a prime ideal and pass to its fraction field. The scalar
coefficient 1 of t in the target survives. This last argument does not
assert general rational-primitive rigidity over rings, or address a
nilpotent multiple of a logarithm that vanishes in every residue field.

## 3. An elementary arithmetic rational-primitive lemma

Let K be any characteristic-zero field. Suppose h in K((t)) has ALL
coefficients in one finitely generated Z-domain A inside K, and

    h'=r(t) in K(t).                                  (5)

Then h is RATIONAL. This statement needs no convergence or D-finiteness.
It is a finite-prime-denominator criterion, stronger than merely being
a G-function. The logarithm is a G-function but fails this coefficient-
ring hypothesis, since its coefficients include 1/p for every prime p.

Here is a full proof. Split the finite poles of r over a finite constant
extension, write

    r=S(t)+sum_(a,j>=1) c_(a,j)/(t-a)^j,

and enlarge the finite arithmetic ring by those constants, coefficients,
the needed denominators, and inverses of distinct-pole differences.
If a residue c_(a,1) is nonzero, invert it as well. Choose a positive-
characteristic fiber of characteristic p larger than the polynomial
degree plus one and all pole orders. Preserve the nonzero rational
denominator; a formal Laurent expansion of (5) still holds in the fiber.
On its scalar Laurent-series field, D^p=0. Therefore

    0=D^p h=D^(p-1) r
     =-sum_a c_(a,1)/(t-a)^p.                         (6)

For j>=2, the derivative multiplier j(j+1)...(j+p-2) contains p.
For j=1, the multiplier is (p-1)!=-1 modulo p; the polynomial part
has vanished by the bound on p. Distinct poles cannot cancel, and the
chosen nonzero residue remains a unit. This contradiction proves that
ALL characteristic-zero residues are zero. The needed fiber exists
after every finite localization, by the same finite-type argument used
in Section 2. No algebraic function-field fiber is being identified with
an unrelated formal branch here: r is rational and only constants were
adjoined.

A rational differential with zero residues has a rational primitive:
integrate its polynomial part and its higher-order polar terms. The
primitive descends to K by taking the normalized trace over the finite
splitting extension. Its difference from h has derivative zero in
K((t)), hence lies in K. This proves the assertion.

The rational residue criterion is classical; it is explicitly recalled
in [Fuernsinn--Pannier2510.00892v3, Section 2](https://arxiv.org/html/2510.00892v3#S2).
Their multiplicative first-order/Kronecker algorithm is not required
for this additive proof. The finite-type fiber existence is
[Chevalley's constructibility theorem, Stacks054K](https://stacks.math.columbia.edu/tag/054K).
No conjecture about higher-rank p-curvature is used.

## 4. Every actual Gaussian rational primitive is rational

Let P be ANY polynomial Gaussian phase over ANY characteristic-zero
field, with fixed nonsingular rational covariance, and let

    V=K(t) Phi(K[X]) subset K((t))

be the FULL actual generic image. For ANY h in V,

    h' in K(t) => h in K(t).                           (7)

For zero or constant h this is immediate. Otherwise, if r=h' is nonzero
rational, h satisfies

    h''-(r'/r)h'=0.                                   (8)

This is globally Fuchsian: r'/r has only simple finite poles and is
O(1/t) at infinity. In the coordinate s=1/t, the coefficient of the
first derivative is still O(1/s) and the zero-order coefficient is zero.
Thus the ACTUAL minimal cyclic module of h is globally regular, even
if (8) is nonminimal or the other original markers are irregular.

The retained [regular-core theorem](GAUSSIAN_G_FUNCTIONS_AND_GLOBAL_REGULAR_CORE.md),
Section 4, now forces the selected critical action psi to be ZERO.
The [critical arithmetic theorem](CRITICAL_NORMALIZATION_AND_RANK_ONE_ARITHMETIC.md),
Section 2, supplies one finite arithmetic coefficient ring for EVERY
original period after removing exp(psi). Since psi=0, it supplies that
ring for the unnormalized periods themselves. For finitely many rational
time multipliers defining h, enlarge that ring by their coefficients
and inverses of their lowest nonzero denominator coefficients. Their
Laurent expansions then have coefficients in this same finite ring.
The lemma in Section 3 applies and proves (7). If the covariance's
arithmetic model uses a finite constant extension K'/K, apply the lemma
there and take the normalized trace: it fixes every coefficient of h
and descends the resulting rational function to K(t).

For a PURE phase, psi=0 already follows from the nonzero regular unit
period 1. So the primitive question in the user's proposal closes
directly: if Phi_(PQ)=1/(1-t), then Phi_Q' equals that rational function,
and (7) would make Phi_Q rational. But dt/(1-t) has residue -1 at t=1
and no rational primitive. Such an original marker Q cannot exist.
Having an original rational period alone does NOT supply an original
polynomial-marker primitive.

Because V is stable under differentiation, induction also proves

    h in V, D^N h in K(t), N>=1 => h in K(t).           (9)

No all-marker regularity assumption is required for (7) or (9). The
field argument does not promote them to arbitrary coefficient rings.
The direct proof of the particular logarithm in Section 2 is independent
of both imported critical-action theorems.

## 5. Embedded unipotent generic modules are excluded

Call a differential module over K(t) unipotent here if it has a finite
flag with every successive quotient the TRIVIAL rank-one module.
Every such differential SUBMODULE U of the actual V above has rank at
most ONE. Indeed, a first trivial submodule has a nonzero horizontal
scalar series, hence a nonzero constant in K; rescale it to 1. If a
second step existed, a lift h of its trivial quotient would satisfy
h' in K(t)1. Equation (7) makes h rational, contradicting independence
of 1 and h over K(t). The same argument applies to the first two steps
of any longer unipotent flag.

Consequently the user's nonsplit rank-two system cannot occur even as
an embedded actual Gaussian submodule, let alone as the WHOLE image.
More generally, the nonrational part of any finite iterated-integral
system with trivial successive quotients cannot be embedded this way.
For example classical polylogarithms Li_s(t), s>=1, generate unipotent
systems over K(t) by the equations D Li_s=Li_(s-1)/t and
Li_1=-log(1-t); none can be an actual Gaussian generic period.

In particular, for the target logarithm L=-log(1-t),

    V intersect (K(t)+K(t)L)=V intersect K(t).          (10)

One can also see (10) directly. If a+bL lies in V with b!=0, its
derivative minus (b'/b) times itself is the rational function
a'-(b'/b)a+b/(1-t). This is nonzero, since otherwise a/b would be a
rational primitive of -1/(1-t). Thus 1 belongs to V, then L belongs
to V, contradicting (7).

These conclusions concern EMBEDDED actual scalar-function submodules.
They do NOT exclude an arbitrary unipotent QUOTIENT of a larger source,
assert semisimplicity of every Gaussian connection, or prove algebraicity
of all regular Gaussian images. Regular irreducible transcendental
connections, such as an elliptic G-function connection, lie outside
this trivial-successive-quotient argument. Algebraic covers and nontrivial
algebraic composition factors require additional arguments; none are
silently trivialized by a tensor product that may lose actual rank.

## 6. Where the old nonsplit extensions keep their irregular term

The retained [Euler-unit classification](GLOBAL_EULER_AND_GAUSSIAN_PULLBACK_LATTICES.md)
has

    A_b=1-t b(theta), theta=tD,
    A_b h_b=1, h_b(0)=1,
    a_(m+1)=b(m)a_m for h_b=sum_m a_m t^m.

If deg b=r>=1 and b has no nonnegative integer root, the quotient has
positive slope 1/r at zero, and h_b has exact Gevrey order r. For
b(z)=a z+c with a!=0, the homogeneous solution is explicitly

    y=t^(-c/a) exp(-1/(a t)).                          (11)

The exponential factor is the irregular term. It cannot be removed
by a rational gauge while a!=0. More generally, finite coefficient
specialization that lowers deg b leaves a positive slope while the
degree is positive. If b acquires a nonnegative integer root, the
selected h_b terminates and its actual cyclic image collapses to K[t], although
an ambient differential source may retain an irregular kernel.

At deg b=0, say b=c, the relation is only multiplication:

    (1-ct)h_b=1, h_b=(1-ct)^-1.                       (12)

The extra GENERIC differential rank disappears; no logarithmic extension
survives. At c=0 the result is h_b=1. In the Kloosterman/squared-spectator
constructions these are the same Euler mechanisms after their retained
transfers. Their [full all-marker proofs](GAUSSIAN_KLOOSTERMAN_QUOTIENTS_AND_FAITHFUL_SYMMETRIC_IMAGES.md)
remain sealed. No new limit of their actual image is inferred from a
nonminimal specialized source operator. Singular rescaling outside this
normalized family would require a new formal period and marker analysis;
the universal logarithm exclusion still applies to any valid Gaussian
specialization.

Thus the requested cancellation does not leave the logarithmic residue
in these constructions. This is an actual-image statement, not a claim
that every nonsplit differential extension is irregular in general.
The abstract logarithmic system is a perfectly valid regular nonsplit
extension; it fails the Gaussian arithmetic constraint.

## 7. Remaining frontier and retained evidence

The user's minimal candidate and the wider original rational-primitive/
unipotent-submodule route are now excluded within the written theorem
scopes above. A regular transcendental ORIGINAL Gaussian image remains
OPEN in other module types. At psi=0, number-field candidates must have
all original markers conjugate-convergent and the stronger finite-prime
arithmetic coefficient property, not merely the general G-function
denominator bound. A single transcendental period with other irregular
markers would not close that full-image gap.

The separate rank-one safety=>purity/primitive-vacuum question remains
OPEN; see its [current reduction](SAFE_SEED_DIVISIBILITY_AND_PHASE_MARKER_FILTRATION.md).
All Oaku operator identities, torsion witnesses, denominators and
specialization thresholds remain unchanged. Sabbah/Maisonobe retain
convergent rational-slope credit; arbitrary nonconvergent formal and
all-ring extensions remain separately scoped, priority unassessed.
No terminal UNKNOWN restart, deep finite campaign, delegation, outreach,
general image-kernel computation or assurance promotion was used.
