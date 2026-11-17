# Nonnegative weighted Gaussian phases: finite Morse reduction and exact safety

16 November2026. Internal all-degree/all-ring proof. No independent,
formal or external assurance, and no priority claim. The
[finite producer/checker](verify_semipositive_gaussian_reduction.py),
[certificate](semipositive_gaussian_reduction_certificate.json) and
[source audit](semipositive_gaussian_reduction_source_audit.json) retain
exact finite identities; universal statements use the argument below.
All earlier sealed packets remain unchanged.

## 1. The entire architecture and its zero-weight boundary

Fix a nonsingular rational Gaussian quadratic form Q0 and nonnegative
rational coordinate weights w_i for which Q0 has weight2. Let P be ANY
polynomial, in any dimension and degree, homogeneous of weight2 for
these SAME weights. Coefficients have weight0. This is a hypothesis on
a supplied weighted Gaussian coordinate system, not a claim every phase
admits one. It is not the ordinary homogeneous-degree-two hypothesis.

Nondegeneracy pairs weights w and 2-w. The zero-weight variables Z_i
pair with equally many weight2 variables W_i. Constant linear changes
within these blocks put the Gaussian form and phase into

    Q0=Z dot W+(1/2)T^T B0 T,
    P=W dot H(Z)+F(Z,T).                                 (1)

B0 is nonsingular rational, all residual weights satisfy 0<w(T_j)<2,
and F has weight2 in T. In particular F has no T-constant or T-linear
term. Every W occurs linearly with a polynomial coefficient H_i(Z),
since any additional positive-weight factor would exceed weight2.
Arbitrary degree in the zero-weight variables is allowed. The covariance
is the inverse precision; this fixed split form is equivalent to standard
Gaussian coordinates over Q(i), as in the earlier algebraic lift.

Write

    F=(1/2)T^T A(Z)T+F_(>=3)(Z,T),
    B(t,Z)=B0-t A(Z), sigma(t,Z)=det(B)/det(B0).

All ring statements below hold over EVERY commutative Q-algebra R,
including non-Noetherian rings and rings with zero divisors.

## 2. A finite polynomial, volume-preserving relative Morse lemma

Over R[t,Z,sigma^-1] there is a POLYNOMIAL triangular automorphism
T=Psi(t,Z,S), identity at t=0 and with identity linear part, such that

    Q0_T(Psi)-t F(Z,Psi)=(1/2)S^T B(t,Z)S,
    det(partial Psi/partial S)=1.                        (2)

Both Psi and its inverse preserve the residual weights and are polynomial
in S. Every coefficient uses only rational constants, t,Z and sigma^-1.
This is stronger than merely asserting a formal or analytic Morse map;
the Jacobian and its effect on the ORIGINAL measure are retained.

Here is a finite construction. Suppose the current potential is q_B plus
terms of ordinary S-degree>=k, k>=3. Let R_k be its degree-k part and set

    delta=-(1/k) B^-1 grad R_k, S -> S+delta(S).           (3)

The cross term of q_B is -S dot grad(R_k)/k=-R_k by ordinary Euler
homogeneity; the remaining error has degree>=k+1. Every delta_i has
weight w_i and ordinary degree k-1>=2. Each variable in one of its
monomials therefore has weight STRICTLY BELOW w_i. Ordering coordinates
by weight makes(3) triangular, with Jacobian1 and polynomial inverse.
The inverse matrix B^-1 respects the same complementary-weight blocks
as B0, so no wrong-weight term is introduced.

There are no monomials of weight2 and ordinary degree greater than
floor(2/min_j w_j). Therefore after finitely many degree eliminations
the error is exactly zero. Composition preserves triangularity, weights,
Jacobian1 and the stated coefficient ring. Every correction is O(t),
so the full map is identity at t=0. This proves(2) over arbitrary R:
only rational units and the explicitly inverted sigma occur.

Classical formal splitting/Morse theory is credited to its usual sources.
We read the full existence argument of Greuel--Pfister,
[2507.17078v1,Theorem2.1(2)](https://arxiv.org/html/2507.17078v1),
and the statement/proof boundary of Theorem2.2. The published record is
[Journal of Algebra689(2026),610-628](https://doi.org/10.1016/j.jalgebra.2025.10.023).
Their general field/algebraic-series result is not silently promoted to
our arbitrary coefficient rings or polynomial Jacobian1 assertion:
the positive-weight finite construction(3) proves those here. No
singularity-isolation, analytic contour or Artin-approximation theorem
is needed for this construction.

## 3. Exact original periods and exact generic rank

Let g in tR[[t]]^n be the UNIQUE vector fixed point and set

    g=t H(g), J=det(I-t H_Z(g)),
    sigma_g=sigma(t,g), f=J^-1 sigma_g^(-1/2).           (4)

The empty determinant is1. All roots and inverses use their constant-one
formal branches. Formal circular Gaussian contraction in the Z,W pairs
and the volume-preserving map(2) give

    Phi[1]=f,
    Phi[Q(Z)]=f Q(g).                                    (5)

The multivariable residue formula is coefficientwise: expanding each W_i
exponential uses E[Z^a W^b]=delta_ab a!, and the change of normal
coordinates Z-tH(Z) has determinant J, a unit modulo t. It holds over
EVERY R. The conditional residual Gaussian uses precision B(t,Z),
contributes sigma^-1/2, and applies the polynomial amplitude substitution
Psi. Its Jacobian is exactly1; no volume factor has been dropped.

For full W-momentum amplitudes introduce sources eta and solve
g_eta=eta+tH(g_eta). Differentiating(5) and its residual-Gaussian version
recovers EVERY original W^alpha marker. The implicit derivatives require
only (I-tH_Z(g))^-1. Wick contraction after(2) requires only B(t,g)^-1.
Thus all specified source-marker images are obtained by FINITELY many
operations involving g, J^-1, sigma_g^-1, Psi and its finite derivatives,
with the common factor f. This retains all amplitudes and is not a
full actual-image-kernel computation.

Over a characteristic-zero FIELD put L=K(t,g_1,...,g_n). The fixed-point
branch is algebraic: its implicit Jacobian at t=0 is identity, and its
dominating component has relative dimension zero. Equivalently the
classical algebraic implicit-function theorem applies. It follows that
L/K(t) is finite. The EXACT FULL generic image is

    K(t) M = f L,       rank_actual=[L:K(t)].             (6)

The upper inclusion follows from the full-amplitude calculation: all
matrix inverses and derivatives are rational in t,g. The lower inclusion
uses original position monomials in(5); their K(t)-span is the field L.
Multiplication by the nonzero algebraic f is injective.

In the NONPURE case f may lie outside L. One must not substitute
[K(t,g,f):K(t)] for the rank in(6). Already P=T^2/2 has f=(1-t)^-1/2,
algebraic function-field degree2 but actual differential rank1. Nor does
the degree of a nonminimal input polynomial give the rank. In the PURE
case f=1, so(6) is exactly the actual algebraic field L.

Every original period is algebraic. The full generic connection has
FINITE geometric monodromy and is regular at EVERY place including
infinity: after a finite cover it is an algebraic gauge of a trivial
connection. Local Puiseux/Kummer bases give only rational Euler exponents
and no positive slopes. The global polynomial lattice is not asserted
semisimple or minimal merely from this generic conclusion.

Consequently regular TRANSCENDENTAL realizability is IMPOSSIBLE throughout
this entire nonnegative weighted-degree-two architecture. This is an
exact obstruction to that proposed route, not a proof of impossibility
for arbitrary Gaussian phases. The
[universal algebraic lift](ALGEBRAIC_FIXED_POINT_LIFTS_AND_UNIVERSAL_GAUSSIAN_IMAGES.md)
is included in(1): all its normal positions have weight0, normal momenta
weight2, and T weight1. Thus this architecture already realizes every
algebraic formal branch, while its transcendental boundary is precise.

## 4. One determinant identity and finitely many fixed linear safety markers

Over EVERY R, purity is exactly the single formal algebraic identity

    P pure-null iff J^2 sigma_g=1.                       (7)

Indeed f(0)=1 and f^-2=J^2 sigma_g. If its square is1, the constant-one
formal branch is1: the other factor f+1 is a unit. This reasoning is
valid with arbitrary nilpotents and zero divisors.

Under purity, ALL-original-marker safety is equivalent to

    g_i(t) in R[t], EVERY zero-weight coordinate i.       (8)

Necessity is(5), since the ORIGINAL LINEAR marker Z_i has period g_i.
For sufficiency suppose g=p(t) is a polynomial vector. Both J and sigma_g
are now polynomials. Identity(7) makes them UNITS in R[t], with explicit
inverses J^-1=J sigma_g and sigma_g^-1=J^2. Their adjugate matrices,
every finite source jet and every coefficient in the polynomial Morse
map evaluated at p are therefore polynomial in t. The common seed
factor is1. The full-amplitude construction in Section3 proves every
marked period polynomial. There is no reduction to closed field points
and no uncontrolled nilradical bound.

This gives a complete safe normal-form criterion for the architecture:
the pure safe members are precisely those admitting p in tR[t]^n with

    p=t H(p),
    det(I-t H_Z(p))^2 det(B0-t A(p))/det(B0)=1.           (9)

Higher residual terms F_(>=3) are COMPLETELY UNRESTRICTED within their
stated positive weights. They affect polynomial marker formulas and may
affect global lattices or cutoffs, but they change NEITHER the original
pure period nor the exact generic image fL nor pure safety. No equality
of ordinary MIXED cutoff ideals is inferred. The full ordinary PURE
coefficient ideals do agree with the quadratic-residual model, because
the pure period is LITERALLY the same series.

Over fields, if p is polynomial, (9) forces J=1 and sigma_g=1 individually
by polynomial units and constant coefficient1. Over general R their
nonconstant coefficients may instead be nilpotent; they must be retained.
This extends the earlier compensator family's two-marker normal form to
ANY number of zero-weight positions and ANY positive residual weights.
It does not give a universal linear-marker test for all Gaussian phases.

The all-ring polynomial-Hom normalization ideal also becomes exact:

    N_P={a in R: a g_i(t) in R[t] for EVERY i}.           (10)

To prove it, multiplication by a identifies vanishing over R/Ann(a)
with scalar-annihilated original coefficients. Purity survives this
quotient, and(8) applies there. Thus(10) detects a times EVERY marked
tail, not merely the coordinate tails. Nonzero normalization alone
does not certify safety; safety requires1 in this ideal.

## 5. No zero weights: complete all-ring eventual safety without assuming purity

If all Gaussian weights are STRICTLY POSITIVE, there are no Z,W pairs.
Then J=1, g is empty and sigma=det(B0-tA)/det(B0) is a polynomial of
degree at most r, the number of residual variables. The nonlinear
part is removed by the finite Jacobian1 map. Over EVERY R,

    all-marker safety
      iff Phi[1] is polynomial
      iff sigma is a UNIT of R[t]
      iff its finitely many nonconstant coefficients generate
          a NILPOTENT ideal I.                          (11)

If f=Phi[1] is polynomial, f^2 sigma=1 is an identity of polynomials,
so sigma is a polynomial unit. Conversely sigma=1+nilpotent polynomial
has a finite binomial square-root inverse. Both B^-1 and the finite
Morse map now have polynomial t-coefficients, proving all-marker safety.
No purity premise is needed here. Over fields(11) reduces to sigma=1;
eventual purity, all-positive purity and safety are equivalent.

The exact normalization ideal in this case is

    N_P=Gamma_I(R)=union_N Ann_R(I^N).                   (12)

This follows by applying(11) over R/Ann(a). The full ordinary all-positive
pure ideal is generated by the r coefficients of sigma-1, equivalently
by the FIRST r original pure moments. The reciprocal-square relation
f^-2=sigma and the binomial series give both literal ordinary ideal
inclusions; no radical replacement occurs. Higher weighted nonlinear
coefficients are free parameters in this ordinary pure scheme. The
field and arbitrary-ring conclusions use finite determinants, not an
unproved Gaussian positivity conjecture.

There is also an elementary independent charge proof and a numerical
bound. Assign charge c_i=w_i-1; the Gaussian is charge-invariant. Every
quadratic phase term has charge0, whereas a residual monomial of ordinary
degree k>=3 has charge2-k<=-1. For a marker of degree d only at most d
nonlinear factors can survive. Thus its period is the FINITE sum

    sum_(j=0..d) t^j E[Q F_(>=3)^j exp(t P_2)]/j!.

For r>=1, if I^L=0, then deg_t f<=r(L-1) and every entry of B^-1 has degree at most
rL-1. With D=max(2,deg P), Wick contraction in this finite sum gives the
explicit sufficient affine cutoff

    m>r(L-1)+ceil(1+(D+1)(rL-1)/2) deg Q.               (13)

For pure inputs sigma=1, take L=1: the intercept vanishes and the same
formula is an all-ring linear cutoff. This conservative bound is not
asserted optimal. The charge argument also explains directly why higher
weighted nonlinear coefficients disappear from the entire pure ideal.
For r=0 the phase is zero and slope0 suffices.

## 6. Boundary, provenance and remaining general work

The original
[weak one-sided neutral-face theorem](NEUTRAL_FACE_ALL_DEGREE_CLASSIFICATION.md)
uses circular CHARGE w(Z)=-w(W), makes the Gaussian quadratic invariant,
and imports a curve exclusion. Here Gaussian weights add to2 and the
Gaussian quadratic has weight2; the zero-weight algebraic fixed points
can be unsafe. These are different architectures. Neither theorem is
being used to justify the other's hypotheses or analytic assurance.

Ordinary Morse theory alone would not give a finite polynomial map or
Jacobian1. Negative weights or a zero-weight residual coordinate destroy
the finite triangular argument; a zero-weight coordinate must instead
be kept among Z and its fixed-point markers tested. The univariate cubic
Gaussian phase has an infinite ordinary Morse expansion and lies outside
the supplied-weight hypothesis. Known irregular squared-spectator and
pullback families remain valid and are not folded into this algebraic class.

This closes an entire all-degree covariant-weight architecture, with
an exact full generic module and ordinary all-ring safety/normalization.
Regular transcendental realization OUTSIDE this architecture, exhaustive
Gaussian safe normal forms, arbitrary actual generic connections and
remaining sharp real-cutoff geometry stay open. The published splitting
lemma is credited at its actual field/formal scope; the finite weighted
application above has no established novelty/priority claim.
Sabbah/Maisonobe retain convergent rational-slope credit; nonconvergent
formal extension and all-ring descent priority stay separately unassessed.
Oaku identities, torsion witnesses, specialization thresholds, protected
endpoints and both terminal UNKNOWN jobs are unchanged. Broad goal active.
