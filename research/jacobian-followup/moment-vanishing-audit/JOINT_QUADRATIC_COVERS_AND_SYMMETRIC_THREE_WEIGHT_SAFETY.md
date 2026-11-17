# Joint quadratic covers and symmetric three-weight safety

16 November 2026. Internal written proofs using retained geometric inputs.
Independent replay, formal verification and external review are not claimed.
Both COMPLETE exceptional MZ cases remain **OPEN**; review is deferred.

On BOTH exceptional quotients the ENTIRE sector with left weights contained
in {-n,0,n}, n!=0, is now closed, with arbitrary neutral terms, mixed
Peter--Weyl degrees and right multiplicities. Over C an eventually pure
phase in this sector has only ONE nonzero weight. Over ANY commutative
C-algebra, the exact criterion is nilpotence of the neutral part AND of
the product of the two charged parts. Every such phase is safe against
ALL original markers, with an explicit nilpotence-dependent cutoff.

Every polynomial recomposition is classified over C as well. This extends
the Long construction to arbitrary independent neutral functions; no
special relation between the neutral part and the balanced product remains.
The support condition is specifically {-n,0,n}. Unequal opposite magnitudes,
other three-weight supports and larger supports are NOT closed here.

See the [checker](verify_haar_three_weight.py),
[receipt](haar_three_weight_verification.json) and
[source audit](haar_three_weight_source_audit.json).

## 1. Two neutral functions determine all unweighted moments

Let a commutative complex algebra have a rational Cstar action and a
normalized invariant linear functional h. Write

    P=p+U+V, wt(p)=0, wt(U)=n, wt(V)=-n, n!=0,
    q=UV, d=p^2-4q.                                    (1)

For EVERY m, weight projection gives

    h(P^m)=sum_(0<=j<=m/2) m!/[j!j!(m-2j)!]
                                  *h(p^(m-2j)q^j).     (2)

No relation between p and q is assumed. On the exceptional algebras the
full joint moment functional of p,q is exactly their normalized compact
sphere period in C[Q2]^O or C[Q2]^I.

More generally, take a smooth complex surface X with regular functions
p,q, an algebraic 2-form Omega and a compact real 2-cycle M of mass1.
Suppose q is not identically zero. Form

    W: ab=q(x), S=p(x)+a+b,
    Theta=Omega wedge da/a=Omega wedge dS/(a-b).         (3)

The conic threefold is integral; its normalization is unnecessary on its
smooth locus. The compact cycle parametrized by x in the base cycle M
and |z|=1, a=z, b=q(x)/z, lies in a!=0 and is smooth. The actual pullback
of Theta is Omega wedge dz/z. Normalize by 1/(2pi i), so its mass is1.
The constant term of (p+z+q/z)^m is exactly the sum in (2). Whenever
the base joint p,q period represents h, as it does on the exceptional
auxiliary sphere, this proves

    h(R(P))=integral_C R(S) Theta/(2pi i)
        for EVERY polynomial R.                       (4)

This transfers the ENTIRE unweighted polynomial functional. It does not
insert a Jacobian or marker, take a scalar branch norm, or reflect purity
from a composite to a primitive phase.

## 2. Compactify the joint map, rather than choosing one height

Let Z be a smooth projective compactification of the base surface X on which

    [1:p:d]: Z -> P2

is a morphism. Resolve its graph as an isomorphism on X and arrange that
the pullback of the line at infinity and the volume divisor are SNC.
Locally write the three homogeneous sections as [F:A:B], with

    p=A/F, d=B/F.

They have NO common zero. At each point of F=0, at least one of A,B is
a unit. This is the essential joint-map step; resolving just one height
would not establish this assertion for two independent neutral functions.

The conic function field is the quadratic cover

    y^2=S^2-2Sp+d, y=a-b.

Its normalization over Z times A1_S is finite, hence projective over
A1_S. With Y=F y, the local equation and ORIGINAL form are

    Y^2=F J, J=B-2SA+S^2F,
    Theta=F Omega wedge dS/Y.                          (5)

At F=0 and A a unit, partial_S J=-2A is a unit. Thus J=0 is smooth
and transverse to every local component of F=0. At a point with A=0,
B is a unit and J is a unit for every finite S. Therefore the branch
divisor of (5), after removing square factors of F, is SNC EVERYWHERE
above F=0 and finite S. This includes all moving intersections and
corners; there is no exceptional parameter that was omitted.

If F=unit*product r_i^b_i and Omega has orders k_i there, removal of
the square part gives exactly the volume prefactor

    product r_i^(k_i+ceil(b_i/2)).                      (6)

Thus the sufficient infinity condition is k_i+ceil(b_i/2)>=0. The
explicit SNC-double-cover blowups in Section2 of the
[quadratic suspension theorem](QUADRATIC_SUSPENSION_AND_ALL_LONG_RECOMPOSITIONS.md)
resolve the cover while preserving the actual residue form. On F!=0,
p,q are regular. If Omega is holomorphic there, the finite conic
calculation ab=q in the same proof gives a resolution with the original
Theta holomorphic, including multiple zeros and critical points of q.
No classification of canonical threefold singularities is imported.

We have proved a general JOINT quadratic-cover criterion: if Omega is
holomorphic where F!=0 and (6) has nonnegative exponents, a smooth
projective model over A1_S carries the ORIGINAL Theta. The compact
cycle is retained in its smooth open locus. The retained
[projective original-volume period theorem](SEMISTABLE_LOG_PERIODS_AND_SINGLE_ISOTYPIC_SAFETY.md)
then excludes eventual vanishing of the moments of EVERY nonzero
polynomial H(S). Nonzero constants use mass1; S is nonconstant along a
generic conic. This is a direct all-H conclusion.

## 3. Exceptional neutral functions supply the joint criterion automatically

For EVERY nonconstant f in C[Q2]^O or C[Q2]^I, the earlier valuation
proof supplies the following TWO properties for ALL boundary divisorial
valuations E on every smooth model:

    if f is finite at E, ord_E(Omega)>=0;
    if f has pole b>0 at E,
             ord_E(Omega)+ceil(b/2)>=0.                (7)

For clarity these are more than the no-pure-height conclusion. The
first is the original finite-volume theorem; the second is Section3
of the quadratic suspension proof. Write f=H/u^D, a=ord_E u,
h=ord_E H, and A the log discrepancy. Then k=A-1-2a. At infinity
b=Da-h, and 2(k+1)+b=2A+(D-4)a-h>0. Its proof uses
A>=max(a,h/m), the simple-root SNC inequality, and D>2m. Orbit sizes
6/12 give D>=3m/6m and D>=3 for EVERY nonconstant invariant f,
including all lower terms and complex coefficients. No new estimate
depending on the degree of p relative to d is required.

Assume at least one of p,d is nonconstant. Over F!=0 they are both
finite, so use the first part of (7) for any nonconstant one. This
makes Omega holomorphic there, by the divisorial criterion on smooth Z.
At an infinity component E, b=ord_E F is the MAXIMUM of the pole
orders of p,d: since [F:A:B] has no base point, at least one numerator
is a unit at its generic point. That function is nonconstant and has
pole order b. The second part of (7) gives exactly (6). Hence the
joint criterion in Section2 applies to EVERY such p,q with q!=0.

## 4. Complete field classification and every recomposition

The original algebra A_Gamma=C[SL2]^Gamma, Gamma=2O,2I, is a domain
with only constant units. If U,V are nonzero, q=UV cannot be constant:
a nonzero constant product makes both functions units, contrary to
their nonzero weights. Thus q is nonconstant, and at least one of
p,d=p^2-4q is nonconstant. Sections1--3 give

    U,V both nonzero => EVERY nonzero polynomial H
                       has H(P) NOT EVEN EVENTUALLY PURE.       (8)

If UV=0, one charged function is zero. Weight projection then gives

    h(H(P)^m)=h(H(p)^m) for EVERY m.                    (9)

If p is nonconstant, a nonzero H gives a nonzero polynomial in that
neutral height, which is not eventually pure by torus closure. If p
is a constant c, the right side of (9) is H(c)^m. Its eventual
vanishing is equivalent to H(c)=0. In that case H(c+U) or H(c+V) has
only strictly one-sided weights, or is zero, and is safe against ALL
original markers. Therefore the exact ALL-RECOMPOSITION law is

    H(P) eventually pure iff
       H=0, or [UV=0, p=c in C, H(c)=0].               (10)

In particular

    P eventually pure iff p=0 and UV=0,
    and every such field phase is one-weight or zero and SAFE.  (11)

This closes the ENTIRE symmetric three-weight sector with arbitrary
neutral terms and mixed isotypic components, in ALL degrees. All
nonlinear recompositions are safe whenever eventually pure. The Long
phase is a special case with neutral part 1+2uv and charged terms
u(1+uv),v; its special discriminant is no longer needed.

The proof uses the exceptional neutral volume properties (7), not
merely that the neutral kernel has zero eventual radical. The abstract
three-weight control in the
[weight-lifting theorem](NONCLOSED_ORBITS_AND_EXACT_WEIGHT_LIFTING.md)
therefore remains valid: a zero neutral radical alone does not justify
lifting from affinely dependent weights. No abstract contradiction or
unproved strengthening of that lifting theorem is asserted.

## 5. Exact criterion over EVERY commutative C-algebra

Let R be any commutative C-algebra and work in A_Gamma tensor R. Let
P=p+U+V have the weights (1). Then

    P eventually pure iff p AND q=UV are nilpotent elements.     (12)

We prove the necessary direction before giving a direct marked bound.
The algebra, its torus weight spaces and the normalized Haar functional
have a Qbar form. Express the three functions in finitely many elements
of a Qbar-defined filtered weight basis. At any prime of R, the finitely
many residue coefficients generate a finitely generated field over Qbar.
Embed that field into C fixing Qbar. Every moment equation and weight
is preserved. The field classification (11) makes p and q zero in that
residue algebra. Injectivity of the field embedding makes their residue
coefficients zero already before embedding.

Thus the finitely many coefficients of p and q belong to every prime
of R, hence are nilpotent. They generate a nilpotent ideal (no
Noetherian hypothesis is needed for a finitely generated ideal of
nilpotent elements). Consequently p and q are nilpotent in A_Gamma
tensor R. This argument does not claim that R itself embeds into C or
that its nilpotents disappear in the original moment equations.

Conversely suppose p^nu=0 and q^kappa=0, nu,kappa>=1. Formula (2)
has no surviving term once m>nu+2kappa-3, so P is eventually pure.
More strongly, let Q be ANY fixed original marker and B_Q the largest
absolute torus weight occurring in it. In a term of QP^m write c for
the number of neutral factors, r and s for the charged counts. Weight
zero requires |r-s|<=B_Q/|n|. A nonzero term has c<=nu-1 and
min(r,s)<=kappa-1, since its product contains p^c q^min(r,s). Therefore

    m=c+2min(r,s)+|r-s|
      <=nu+2kappa-3+B_Q/|n|.

Hence the explicit FULL original-marker cutoff is

    m>nu+2kappa-3+B_Q/|n|
       => h_R(QP^m)=0.                                (13)

For the standard left torus one may take B_Q<=deg_matrix(Q). The
argument proves literal zero of every weight-zero term; it uses no
positivity, domain assumption, reducedness, cancellation estimate or
uniform numerical nilpotence exponent. It establishes MZ safety over
EVERY R. Nilpotence characterizes EVENTUAL purity; it need not kill
the first moments. The field recomposition classification (10) is not
copied verbatim to arbitrary rings.

## 6. Full target and the precise remaining weight gap

The conclusion transports to every maximal left torus and under the
retained algebraic automorphism invariance of Haar. Ring safety also
transports after scalar extension; an unchanged matrix-degree bound
is not asserted for a nonlinear automorphism.

Every full failure must avoid support contained in {-n,0,n} after
EVERY algebraic automorphism and every maximal left torus. This is
stronger than the prior two-weight exclusion, but it does NOT say
every failure has at least four weights: an UNEQUALLY SPACED
three-weight support is not covered. A shifted arithmetic progression
whose middle weight is nonzero is not covered either. Such supports,
larger supports, and the general mixed-weight component obstruction
remain in the full problem. No new coefficient/support/prefix campaign
is proposed.

The distinct Hopf phase has a further charged weight and is not
automatically of form (1). Its general nonlinear composites are not
excluded by this theorem. Nor is every regular non-equivariant Bezout
pullback reduced to this sector. The FULL target stays actual
unweighted exceptional purity forcing all corrected local torus orders
zero, or an actual invariant pure failure. Both COMPLETE exceptional
MZ cases remain OPEN. Review remains DEFERRED, not the next gate.
