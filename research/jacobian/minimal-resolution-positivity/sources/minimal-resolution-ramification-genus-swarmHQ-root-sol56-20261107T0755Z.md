# Ramification arithmetic genus on the actual minimal resolution

Claim ID: MINIMAL-RESOLUTION-RAMIFICATION-GENUS-1. Producer: swarmHQ (ROOT/Astra configuration). Independent reviewer: gpt-5.6-sol / ultra configuration; actual hosted identity UNVERIFIED. Evidence tier MANUAL, relative to named foundations and the retained smooth-finite interface. Lifecycle PROVISIONAL; independent CONFIRMED_AT_NAMED_FOUNDATIONS. Novelty unassessed. JC2 remains unresolved.

## Reviewed statements

Work over C. Let g:Y -> P^2 be a generically finite morphism of degree N>1 from a smooth irreducible projective surface. Suppose U is an open subset isomorphic to the whole A^2, g is etale on U, and g(U) misses a target line L. Assume g contracts no smooth rational (-1)-curve. Put H=g^*L and let R_g be the effective Cartier divisor of det(dg), including its multiplicities on contracted curves. Then, conditional on the named foundations and the supplied PROJECTIVE-RAMIFICATION-GENUS-1 interface,

    p_a(R_g) = 1 - chi(O_{R_g}) >= 1.

The statement is proposed in particular for g=h o q, where h:X -> P^2 is the finite projective normalization in the source function field of a hypothetical degree N>1 plane Keller map, and q:Y -> X is its minimal resolution. It does not assume X smooth, Q-Gorenstein, or rational-singular; it does not assume g finite, H ample, or the original affine Keller map finite/surjective.

## Named imports and accepted interface

Imports: finite normalization of a projective variety in a finite function-field extension; Zariski's Main Theorem for quasi-finite birational maps to a normal target; projective minimal resolution of a normal complex surface, with no exceptional smooth rational (-1)-curves; smooth surface cone/contraction theorem; Hodge index for a nef divisor of positive square; adjunction/Riemann--Roch; classification and Picard/canonical formulas of Hirzebruch surfaces; divisor-class localization; proper quasi-finite implies finite; and characteristic-zero DVR ramification/different and e*f<=N. No fresh source-proof audit is claimed.

PROJECTIVE-RAMIFICATION-GENUS-1 is used only as the supplied accepted statement: if h:X -> P^2 is finite of degree N>1 with X smooth projective and a whole U=A^2 open on which h is etale and misses a target line, then p_a(R_h)>=1. Its MANUAL/PROVISIONAL, independently confirmed named-foundation scope and its retained CUBIC-SCROLL-DONOR-1 dependency are unchanged. Exact interface: AUDIT.md26885-26894 at publicc9f508a77890db2a70597826a51bc9a34ea20da7. No fresh re-audit of its proof.

## Proof of the reviewed statement

Y is rational because it contains the dense open U=A^2. The actual pullback H is nef with H^2=N>0. The determinant of the differential is a nonzero section of omega_Y tensor g^*omega_{P^2}^{-1}, so its effective divisor R=R_g satisfies R linearly equivalent to K_Y+3H. This constructs R directly on the smooth Y: no pullback of a possibly non-Q-Cartier canonical divisor on X is used. Since g is etale on U, every prime of R lies outside U.

Suppose p_a(R)<=0 and put Q=2K_Y+3H=K_Y+R. Adjunction gives Q.R=2p_a(R)-2<0, so Q is not nef. On the K_Y-nonnegative part of the curve cone, Q is nonnegative because H is nef. Thus a K_Y-negative extremal ray has Q negative.

A birational contraction has a smooth rational (-1)-curve E. If H.E=0, then g contracts E, forbidden by the hypothesis. Otherwise H.E is a positive integer and Q.E=-2+3H.E>=1, a contradiction. A contraction to a point forces Y=P^2; H=O(m) with m>0, and Q negative forces m=1 and N=1, also impossible.

The remaining extremal contraction pi:Y -> B is to a smooth curve with connected fibres. Rationality gives B=P^1. For a general fibre F, F^2=0 and K_Y.F=-2. H.F cannot be zero: Hodge index and H^2>0 would give F^2<0 for the nonzero numerical class F. Hence H.F is a positive integer. Q.F<0 gives -4+3H.F<0, so H.F=1 and R.F=1.

Every curve contracted by pi belongs to its single extremal numerical ray. Thus each prime component of a scheme fibre is numerically a positive multiple of F and has self-intersection zero. A connected reducible fibre has a component meeting another component; from C.F_b=0 its self-intersection is then negative. This contradiction proves every fibre is irreducible. If F_b=mC, then 1=H.F_b=m(H.C), hence m=1. Flatness over the smooth base curve preserves arithmetic genus zero, so each reduced irreducible fibre is smooth P^1. Therefore Y is a P^1-bundle over P^1, hence a Hirzebruch surface F_n, n>=0. This step uses the extremal contraction, not ampleness of H.

Write Pic(F_n)=Z[S] plus Z[F], S^2=-n, S.F=1, F^2=0. Since H.F=1, H=S+bF. Nefness gives b>=n and N=H^2=2b-n>1. If b>n, H is ample. No positive-dimensional fibre of g can then exist, since a curve in such a fibre would have H-degree zero. Proper quasi-finiteness makes g finite, and PROJECTIVE-RAMIFICATION-GENUS-1 contradicts p_a(R)<=0.

It remains to consider b=n. Then N=n>=2, H.S=0, and g contracts S. The differential has rank at most one at the generic point of S, so S occurs in R with positive coefficient. But R.F=1 and S.F=1; therefore its coefficient is exactly one and R has no other horizontal prime. All remaining primes are fibres of pi. The canonical formula gives

    R = S + sum_j r_j F_j,    sum_j r_j = 2n-2 = 2N-2,

where the F_j are distinct fibres and each r_j is a positive integer. The equality of the coefficient sum follows from R linearly equivalent to K_Y+3H=S+(2n-2)F.

Every F_j has H.F_j=1, so g does not contract it. At its generic point, the coefficient in R equals e_j-1, where e_j is its DVR ramification index over the target image divisor. One may justify this by factoring g through the finite normalization: the birational map to that normal surface is an isomorphism at a noncontracted prime's generic point, so the usual finite DVR different formula applies there. The function-field extension has degree N, giving e_j<=N and hence r_j<=N-1. Since the positive coefficient sum is 2(N-1), at least two distinct F_j must occur. Thus support(R) contains at least three distinct curves: S and two fibres.

All those primes lie in Y minus U. However divisor-class localization gives

    O(Y)^* -> O(U)^* -> direct_sum_boundary Z[Z] -> Cl(Y) -> Cl(U) -> 0.

The unit groups are both C^* and Cl(U)=0, so the boundary-prime group maps isomorphically to Pic(Y)=Z^2. There are exactly two boundary curves. The three distinct ramification curves contradict this. Both ruled cases are impossible, proving the proposed p_a(R_g)>=1.

## Application to the actual finite projective normalization

Let F:A^2 -> A^2 be a hypothetical Keller map of generic degree N>1. It is etale and quasi-finite. Normalize target P^2 in C(x,y), obtaining a normal projective X and finite h:X -> P^2 of degree N. Zariski's Main Theorem identifies the whole source A^2 with an open U of X over the target affine chart. Equivalently, on that chart the integral closure of C[F_1,F_2] in C(x,y) lies in the integrally closed source ring C[x,y]; the resulting quasi-finite birational map to the normal normalization is an open immersion. No equality with the whole finite inverse image of A^2 is asserted.

The minimal resolution q:Y -> X is an isomorphism on the smooth U. Thus Y is smooth projective, contains the same whole U, and g=h o q is etale there and misses the target line at infinity. Since h is finite, a curve is g-contracted exactly when it is q-exceptional. Minimality of q excludes exceptional smooth rational (-1)-curves. All hypotheses of the proposed general statement hold, without making X smooth or g finite. Hence this would be a necessary condition on the actual minimal resolution of any counterexample.

## Checks, limits, and review targets

The negative degree sign is used only on a K-negative extremal ray. H is nef, never treated as ample until b>n. The irreducibility-of-fibres argument uses the one-ray contraction and connectedness; it must not be applied to an arbitrary ruling on a blown-up surface. The different exponent e-1 is used only for noncontracted F_j, never for the contracted S. Residue degree need not be one for the e_j<=N inequality. No equality of open complements, rationality of singularities, or Cartier property on singular X is assumed.

Resolution sensitivity check: if beta:Y' -> Y blows up a boundary point and g'=g o beta, then H'=beta^*H and R_{g'}=beta^*R_g+E. Adjunction gives p_a(R_{g'})=p_a(R_g)-1. The new exceptional (-1)-curve is g'-contracted, exactly violating the minimality hypothesis. Thus arbitrary boundary blowups/SNC refinements cannot be silently substituted for the minimal model. This is a divisor transformation check, not a constructed degree N>1 Keller example.

Arithmetic genus can be positive even when every component is rational, because of incidence, singularities or multiplicity. This candidate does not force p_a(R_g)<=0, reduced ramification, properness, surjectivity, an actual Keller pair, or JC2. It removes the smooth-finite restriction only from this necessary arithmetic-genus obstruction, if the proof survives review. No novelty claim.

Hostile review targets: the cone argument with merely nef H; exclusion of birational rays from minimality; irreducible/reduced fibres; the nonample F_n endpoint and coefficient sum; the DVR bound on each noncontracted fibre; the whole-A^2 boundary basis; and the actual normalization/ZMT/minimal-resolution application. Check the boundary-blowup genus formula. Verdicts: CONFIRMED_AT_NAMED_FOUNDATIONS / GAP_FOUND / FALSE_WITH_CONTROL / INCOMPLETE. No promotion before independent review.

## Independent review and integration

The [independent review](minimal-resolution-ramification-review-sol56-20261107T0740Z.md), full SHA256 f5140956ebc907d7147f54ce7bb1e7d526496e5138350acfaa88a7f66c04634a, confirms the general statement, actual-normalization application, and boundary-blowup formula at the named foundations. The proof above preserves the frozen producer argument; its proposal/review-target wording records what was submitted, not unresolved review debt. No mathematical repair was needed.

The [evidence bundle](../box/minimal-resolution-ramification-genus-swarmHQ-20261107T0755Z/README.md) contains exact producer bytes, the supplied finite-interface excerpt, and compact review custody. The producer was frozen at HQ commitf0ddb20c0ee4c63faf40c6c1cc9d95c9de68d7a3 before review; its canonical artifact basis ise2fca4ba314e0abfac3768fcc3f13241871dace9. ROOT verified the seal, all18 immutable inputs/tool/prompt pre-post hashes, all37 reading windows covering1743 nonSTATE+83 STATE LF, and read the whole125-LF report. A minor accounting correction counts30 help LF conservatively against the review's29, within the100-LF cap; the sealed review is unchanged.

This removes the smoothness/finiteness restriction only from this necessary arithmetic-genus obstruction. It supplies no opposite genus bound, no integrality or target descent, and no JC2 conclusion. It is not a theorem about arbitrary further resolutions: each additional boundary blowup lowers the arithmetic genus by one. [The earlier smooth-finite theorem](projective-ramification-genus-swarmHQ-root-sol56-20261107T0715Z.md) remains an unchanged conditional dependency.

## OPEN(S) RAISED

None. The existing whole-polynomial/no-fold contradiction gap remains open; no automatic successor or new OPEN identifier is introduced.

## COLLISIONS

status: EMPTY

- NONE — the report contains no explicitly raised `OPEN[...]` entries.

Corpus basis: `c9f508a77890db2a70597826a51bc9a34ea20da7` (Git blobs only).

<!-- BODY-END -->

## Seal

- Body definition: every byte through the unique standalone `<!-- BODY-END -->` line,
  including its terminating newline; this seal is outside the body.
- Body bytes: `11613`.
- Body SHA-256:
  `ea3a23cc70bf1ff782c9f782979804cfe1b44c125aa6aeb9fee5ae0acc368c25`.
- Frozen basis: `c9f508a77890db2a70597826a51bc9a34ea20da7`.
