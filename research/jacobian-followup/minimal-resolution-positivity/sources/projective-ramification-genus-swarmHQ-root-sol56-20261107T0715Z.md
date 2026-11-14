# Arithmetic genus of ramification in a smooth finite projective completion

Claim ID: PROJECTIVE-RAMIFICATION-GENUS-1.
Producer: swarmHQ (ROOT/Astra configuration). Independent reviewer: gpt-5.6-sol / ultra configuration; actual hosted model identity is unverified.
Evidence tier: MANUAL, relative to the named standard foundations and retained cubic-scroll interface. Lifecycle: PROVISIONAL; independent verdict CONFIRMED_AT_NAMED_FOUNDATIONS. Novelty unassessed. JC2 remains unresolved.

## Statement

Work over C. Let X be a smooth irreducible projective surface, h:X -> P^2 a finite morphism of degree N>1, and L a target line. Assume there is an open U isomorphic to the whole affine plane A^2, h is etale on U, and h(U) misses L. Put H=h^*L and let R be the full ramification divisor, including different multiplicities. Then

    p_a(R) = 1 - chi(O_R) >= 1.

In particular R cannot be a nonempty reduced SNC forest of rational curves. The arithmetic-genus convention includes disconnected and nonreduced divisors.

This is a conditional obstruction to smooth finite projective completions. It supplies neither smoothness of an actual projective normalization nor finiteness after resolving it.

## Dependencies and exact inherited interface

The named imports are: adjunction/Riemann--Roch for Cartier divisors on a smooth projective surface; the smooth surface cone and contraction theorem; the classification of smooth rational ruled surfaces; finite-map Riemann--Hurwitz and tame DVR different exponent e-1 in characteristic zero; divisor-class localization; and the canonical/Picard formulas for Hirzebruch surfaces. These are mathematical imports, not newly audited source proofs.

The retained CUBIC-SCROLL-DONOR-1 interface excludes a dominant everywhere-defined morphism A^2 -> F_1 minus (support Ram(h) union support h^*L) for every finite cubic h:F_1 -> P^2 and every target line L. It imposes no bound on the degree of the first morphism. Its MANUAL/BOOK-relative reviewed scope is used unchanged; its proof was not re-audited here. See [the retained integration](cubic-scroll-donor-exclusion-root-20260912.md), full SHA-256 f5abd9086559b1cfdd6900baafcf339c1a8ad394917ab3141d23bc999403545a, and AUDIT.md lines20885-20920 at basis bc9a5972e4215486b0bb3fc6ae419f07781dbb9f. The exact interface excerpt is included in the evidence bundle.

## Proof

X is rational because it contains U=A^2 as a dense open. Finiteness makes H ample and H^2=N. Riemann--Hurwitz gives R linearly equivalent to K_X+3H; R is effective and supported outside U. It is nonzero: otherwise h is finite etale, and the simply connected projective plane admits no nontrivial connected finite etale cover.

Suppose p_a(R)<=0. Let Q=2K_X+3H=K_X+R. Adjunction gives

    Q.R = 2p_a(R)-2 < 0,

so Q is not nef. On the K_X-nonnegative part of the curve cone, Q is nonnegative. The cone theorem therefore supplies a K_X-negative extremal ray with Q negative.

A birational contraction has a (-1)-curve E; then Q.E=-2+3H.E>=1, impossible. A contraction to a point forces X=P^2. Writing H=O(m), negativity gives -6+3m<0, hence m=1 and N=1, also impossible. The remaining contraction pi:X -> B is to a smooth curve. Its general fibre F has K_X.F=-2. Since H.F is a positive integer, -4+3H.F<0 implies H.F=1. Therefore R.F=1.

Here X really is a Hirzebruch surface, without an extra relative-minimality assumption. Rationality gives B=P^1. Any scheme fibre F_b=sum m_i E_i satisfies 1=H.F_b=sum m_i(H.E_i). Ampleness makes every summand a positive integer, so the fibre is reduced and irreducible. Flatness over the smooth curve gives arithmetic genus zero for every fibre, which makes each smooth and rational. Thus pi is a smooth P^1-bundle over P^1 and X=F_n for n>=0.

Since R.F=1, R has one horizontal prime C, of coefficient one and C.F=1. The effective divisor T=h^*L likewise has a unique horizontal prime B_0 of coefficient one. C and B_0 are distinct. If C occurred in T, finiteness would give h(C)=L. Its coefficient one in R equals e_C-1, so e_C=2. Its coefficient in T would then be two, contradicting T.F=1. This argument does not assume residue degree one.

Every prime in R or T lies in X minus U. Divisor-class localization gives the exact sequence

    O(X)^* -> O(U)^* -> direct_sum_boundary Z[Z] -> Cl(X) -> Cl(U) -> 0.

The unit groups are both C^* and Cl(U)=0. Thus boundary-prime classes form an integral basis of Pic(X)=Z^2. There are exactly two boundary curves, already the distinct C and B_0. Consequently no vertical component of R or T exists: R=C and T=B_0. Possible isolated missing points do not change the codimension-one sequence. Only the inclusion U into X minus (support R union support T) is used.

In the standard basis S,F of Pic(F_n), write S^2=-n, S.F=1, F^2=0. Then

    H=S+bF, b>n, N=2b-n,
    K_X=-2S-(n+2)F, R=S+(3b-n-2)F.

The integral boundary basis R,T=H has determinant of absolute value one, giving

    |b-(3b-n-2)| = |2-N| = 1.

Since N>1, N=3. The conditions 2b-n=3, b>n, n>=0 force n=1 and b=2. Thus the same h is a finite cubic map F_1 -> P^2. The original open immersion U=A^2 is a dominant everywhere-defined morphism into F_1 minus (support R union support h^*L), precisely what CUBIC-SCROLL-DONOR-1 excludes. Therefore p_a(R)>=1.

A nonempty reduced rational SNC forest with c connected components has chi(O_R)=c and p_a(R)=1-c<=0. The asserted forest consequence follows.

## Checks, review, and limits

The independent review attacked the contraction cases, singular fibres, different/pullback multiplicities, the boundary class sequence in the presence of possible isolated missing points, the determinant, and exact compatibility with the cubic-scroll interface. It confirmed each at the named-foundation scope. [Original sealed review](projective-ramification-genus-review-sol56-20261107T0655Z.md), full SHA-256 b8df842acdfe27a6cca42c582ec321fc1077c0c97bf34f83fbb17af6ccb4b5f7. Some display escapes in that original are malformed; its bytes are preserved. This integration supplies readable notation without changing its claim.

The frozen producer proposal, exact inherited interface, and compact review-custody evidence are in [the bundle](../box/projective-ramification-genus-swarmHQ-20261107T0715Z/README.md). The producer proposal was frozen at HQ commit5eeea7ec5472ba0eb4145eeddba064bdbc556af7 before review; its canonical artifact basis is8b11a9d642cb0ae05ea57ceb1edbddf03a6d0083. Independent input-window auditing verified all required 1785 non-STATE and82 STATE LF, including103 LF of mathematical proposal/interface, against the frozen packet. The configured reviewer differs from the producer configuration. No computational theorem or formal proof is claimed.

The degree-one contraction endpoint, possible disconnected R, and vertical-component alternatives were checked rather than silently excluded. A positive arithmetic genus does not imply a component of positive geometric genus: incidence cycles, singularities and multiplicities remain possible. Thus the earlier positive-geometric-genus donor obstruction does not close this gap. No counterexample satisfying all hypotheses was found; this is not an exhaustive search.

No singular normalization, arbitrary resolution, nonlinear target affine chart, generic-cover condition, or higher-degree donor exclusion is asserted. Resolving a singular finite normalization can destroy finiteness and ampleness, exactly where this proof uses them. Esterov and the recent generic-cover manuscript are not premises of this proof. Selected historical comparisons do not establish novelty.

## OPEN(S) RAISED

None. The existing singular-normalization/whole-source descent gap is unchanged; no new OPEN identifier or automatic successor is introduced.

## COLLISIONS

status: EMPTY

- NONE — the report contains no explicitly raised `OPEN[...]` entries.

Corpus basis: `bc9a5972e4215486b0bb3fc6ae419f07781dbb9f` (Git blobs only).

<!-- BODY-END -->

## Seal

- Body definition: every byte through the unique standalone `<!-- BODY-END -->` line,
  including its terminating newline; this seal is outside the body.
- Body bytes: `7951`.
- Body SHA-256:
  `781a4dde72e01a80342d9a1486930b83d4855b34d2cce6ec51ee312ffd4c1c30`.
- Frozen basis: `bc9a5972e4215486b0bb3fc6ae419f07781dbb9f`.
