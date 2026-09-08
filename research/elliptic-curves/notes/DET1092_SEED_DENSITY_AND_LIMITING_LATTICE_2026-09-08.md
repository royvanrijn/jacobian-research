# Why the conic seeds look alike, and what a seed does not explain

## Result and current evidence

**New deduction from verified generic data and established specialization
theory.** The orbit8044 factory has two different asymptotic properties:

1. Its split parent parameters form a thin set: there are `Theta(B)` of height
   at most `B`, among `Theta(B^2)` rational parent parameters. Nevertheless,
   all but finitely many rational factory addresses supply 18 independent
   sections after specialization.
2. The normalized height lattice of those **displayed 18 sections** tends to
   one fixed positive definite lattice. The seed's squared orthogonal distance
   from the inherited span, normalized by the inherited squared minimum,
   tends to **3/8**. The same limit holds with the fixed-basis median used in
   the completed landscape comparison. All parity-CVP minima have corresponding
   uniform limits. Neither assertion requires a point search or an exceptional
   point as input.

This explains a seed-incidence bottleneck and the similar geometry of the
factory's seeds. It does **not** explain the arithmetic existence of all the
additional points on302, or prove that the bounded null fibres lack them.

**Verified applications, completed snapshot.** The separately produced
[M18 comparison](M18_LANDSCAPE_COMPARISON_2026-09-08.md) has a sealed package
with 15 states on 12 distinct `j`-invariants. Its lightweight hash, archive and
outcome-binding checker passes again here; no search or CVP is repeated.

| Completed experiment | Certified subgroup outcome | Boundary |
|---|---|---|
| Six new orbit8044 factory fibres | 18 to18, each114 charts | Bounded no gain, not exact rank18 |
| Conic fibre `s=1926/2699` | 18 to21,568 charts | Terminal independent V3 replay, then bounded no gain |
| Three different known seeds on302 | 18 to31 in855,637,385 charts | One repeated fibre; calibrated, known-seed runs |

The remaining302 panel case has no sealed outcome in that snapshot. The
[earlier seed-pair audit](../../artifacts/generated-results/elliptic-curves/curve302_seed_pair_arithmetic_v1.json)
already proves rank19 for M17 plus the first two successful seeds. They are
independent over the generic rational span, not two presentations of a single
exceptional direction. These summaries are retrospective evaluation only;
none enters the construction/checker below.

**Verified numerical diagnostic, not exact canonical heights.** The six
factory Schur ratios in the sealed rounded metrics are

```
0.3836779983, 0.3725537371, 0.3831436200,
0.3706199230, 0.3878195334, 0.3712697802.
```

They surround the predicted `3/8 = 0.375`. This is consistency with an
asymptotic theorem, not a certified finite-height error estimate. The rank21
success has ratio0.3898957, while the bounded null `2980/1967` has0.3899692.
The three completed302 seeds have much larger ratios2.254–5.709. Neither
large seed distance nor a count of cheap cosets is a demonstrated general
amplification predictor. The conic does not split at302, so its specialized
18-section theorem does not describe those302 seed packets.

## 1. Audit: what is reused

**Verified applications already completed.** We use, rather than re-prove:

- The determinant1092 parent's full geometric and rational MW rank17, height
  Gram `G`, and24 irreducible nodal fibres.
- The [orbit8044 base-change theorem](DET1092_RANK18_BASE_CHANGE_AND_INITIAL_UNLOCK_2026-09-08.md):
  a rational conic, a degree-two map, an explicit section `Q`,
  `Q + sigma(Q) = P_w`, `w^t G w = 10`, and independent18-section height Gram.
- The [factory parametrization](ORBIT8044_SEED_FACTORY_2026-09-08.md), which
  is a different coordinate on that same conic, not a new generic construction.
- The [half-lattice propositions](ADAPTIVE_HALF_LATTICE_VISIBILITY_2026-09-07.md#3-intrinsic-half-lattice-distance-and-covering):
  nesting decreases intrinsic CVP distance, but does not by itself control
  chart-coordinate height or a finite scheduling policy.

The old conic progression already gives an explicit infinite set with
point-by-point independence certificates. The new all-but-finite statement
below concerns *all rational addresses* on the conic; it supplies no effective
exceptional list or guarantee that a bounded certificate backend succeeds.

## 2. Seed incidence: a thin image with abundant constructed seeds

**Verified application.** Write the current factory map as

\[
S(u)=\frac{5193-125501628u-193042839045u^2}
 {35630-130270680u-1324582766150u^2}.
\]

For its homogeneous numerator `F(a,b)`, denominator `Y(a,b)` and retained
ordinate `Z(a,b)`, the checker verifies

\[
Z^2=4940136120823281Y^2-21931265870123820FY
       +51437496474840100F^2.
\]

The conic is smooth and the degree-two forms `F,Y` have no common projective
zero. Thus the map from `P1_u` to the conic is an isomorphism: its composition
with the conic's degree-two projection already has degree two. Away from the
finite exceptional chart/branch parameters, its rational image is precisely
the rational splitting locus of this bisection.

**New deduction, elementary height count.** For primitive integers `[a:b]`,
put `H(u)=max(|a|,|b|)`. The certificate's4-by4 Sylvester matrix and adjugate
give integer linear-form identities for `R a^3` and `R b^3` in the ideal
`(F,Y)`, where

\[
R=-535452514534253135498797198080000\ne0.
\]

Consequently `gcd(F(a,b),Y(a,b))` divides `R`: take valuations of both
identities and use `gcd(a,b)=1`. Bounding those same identities at the
archimedean place and bounding the forms by their coefficient sums gives

\[
\frac{H(u)^2}{K}\le H(S(u))\le U H(u)^2,
\quad
K=214157082266752155892137085276800,
\quad U=1324713072460.
\]

These deliberately coarse constants are exact, not a useful small-height
estimate. Since there are `Theta(T^2)` rational projective points of height
at most `T`, and a degree-two map has at most two preimages of any point,

\[
\#\{s\in S(\mathbf P^1(\mathbf Q)):H(s)\le B\}=\Theta(B).
\]

The proportion inside all rational parameters of height at most `B` is
therefore `Theta(1/B)`. This explains why scanning generic parameters and
parametrizing the conic are very different incidence experiments. It does
not predict the exact three splits in the earlier ten-million-address
intake, whose finite selection is not a uniform rational-height sample.

The exact deck transformation is also retained:

\[
\iota(u)=\frac{-368291203u+3162606137}
 {-117574559717537993u+368291203}.
\]

The checker verifies `S(iota(u))=S(u)`, `iota^2=id` projectively, and reversal
of the conic ordinate. The two sheets cannot supply two independent extras:
the already proved relation is `Q + sigma(Q) = P_w`.

## 3. The limiting seed lattice

**Established literature.** Silverman's height-limit and specialization
theorems give convergence of specialized height pairings divided by the
parameter's logarithmic height, with consistent height conventions. Positive
generic regulator then implies independence at sufficiently large parameter
height. See [Silverman, *Canonical Heights in Families*, pp.24–26](https://swc-math.github.io/aws/2024/2024SilvermanLecture4Slides.pdf).

**Verified application.** The parent's rational `j` map has degree24, and its
pullback has degree48, so the family is nonisotrivial. The existing18-section
Gram is

\[
\mathcal G=\begin{pmatrix}2G&Gw\\w^tG&8\end{pmatrix},
\qquad\det\mathcal G=429391872.
\]

The checker rederives it from generic data and verifies the exact identity

\[
\boxed{\|(x,k)\|_{\mathcal G}^2
       =2\|x+kw/2\|_G^2+3k^2.}
\tag{1}
\]

In particular, the seed's squared distance from the inherited real span is3.
Every nonzero generic section on the original24I1 K3 has height
`4+2(P.O)>=4`, and displayed sections attain4. Thus the inherited lattice on
the cover has squared minimum8. The median of the first17 displayed diagonal
entries is also8; the intrinsic minimum and this coordinate-dependent median
are different statistics which happen to agree in the generic limit.

**New deduction.** Let `A(u)` be the actual canonical-height Gram of the
specialized displayed18 sections. For a positive constant `c` accounting only
for a chosen height convention,

\[
A(u)=c\,h(u)\mathcal G+o(h(u)),\qquad h(u)=\log H(u).
\tag{2}
\]

There are only finitely many matrix entries, so the height-limit theorem
gives this in matrix norm. Positive definiteness of `mathcal G` gives
independence for all sufficiently large `h(u)`. Northcott finiteness over
`Q` then proves that only finitely many rational addresses can fail
independence (also remove the finite singular/undefined set).

Schur complement is continuous near a positive definite matrix. Applying
it to (2), and likewise to the inherited squared minimum, proves

\[
\boxed{
\frac{\operatorname{dist}(Q_u,\operatorname{span}_{\mathbf R}M_{17,u})^2}
 {\lambda_1(M_{17,u})^2}\longrightarrow\frac38.}
\tag{3}
\]

The same limit holds on replacing the denominator by the fixed-basis median.
The numerator is unchanged by translating the seed by an inherited point,
negating it, or choosing the conjugate conic sheet. Formula (3), using the
lattice minimum, is intrinsic under unimodular changes of the inherited basis.
The limiting lattice in (2) transforms by congruence, not by a new geometry.

**New deduction: uniform parity-CVP limits, without enumeration.** If
`mu_A(r)=min_{z congruent r mod2} z^t A z`, matrix convergence implies, for
every `epsilon>0` and sufficiently large `h(u)`,

\[
(1-\epsilon)\mathcal G\preceq A(u)/(c h(u))
 \preceq(1+\epsilon)\mathcal G.
\]

Taking minima gives the same two multiplicative bounds on every `mu`,
simultaneously for all `2^18` parities. No subset census is needed. Sorted
minima converge; counts below a threshold stabilize if that threshold avoids
all limiting minima. Exact minimizing multiplicities can change at ties,
and no assertion about retained reducer output or chart heights follows.

For constructive use, (1) also reduces the generic template to17-dimensional
coset problems. Define `mu_G(r)` modulo2 and `nu_G(v)` modulo4. Then

\[
\begin{aligned}
\mu_{\mathcal G}(r,0)&=\min\{2\mu_G(r),\ 12+2\mu_G(r+w)\},\\
\mu_{\mathcal G}(r,1)&=3+\tfrac12\nu_G(2r+w).
\end{aligned}
\tag{4}
\]

Proof: replacing `(x,k)` by `(x-2w,k+4)` preserves its parity and the first
term in (1). A minimum therefore has `k=0,2,-2` in the even case or `k=1,-1`
in the odd case. The two nonzero signs agree by negation. Substitution gives
(4). This is an exact identity for the generic template, not an instruction
to retune the running V3 policy.

The rounded numerical V3 metrics are not exact height certificates. The
analytic theorem applies to them only if their normalized numerical errors
tend to zero; the present code's finite precision is not such a proof.
Likewise, convergence does not prohibit useful finite-height deviations or
an arithmetic predictor using information beyond the displayed subgroup.

## 4. The next arithmetic theorem, rather than another score

**Established algebra, applied to the existing cover.** Let the reduced
parent short equation be `E_s: y^2=x^3+A(s)x+B(s)` and put

\[
q(s)=4940136120823281-21931265870123820s+51437496474840100s^2.
\]

Its explicit quadratic twist is

\[
T_s:\quad Y^2=X^3+q(s)^2 A(s)X+q(s)^3B(s).
\]

The polynomials `A,B` are the existing reduced-parent equations, not fitted
to any successful fibre. Over `K=Q(s)`, the two eigenspaces of the quadratic
involution, after tensoring with `Q`, give

\[
\operatorname{rank}E(K(\sqrt q))
 =17+\operatorname{rank}T(K).
\tag{5}
\]

Indeed `1+sigma` and `1-sigma` project to the invariant and anti-invariant
spaces; the latter is the twist under `(X,Y)=(qx,q^{3/2}y)`. The known
`Q-sigma(Q)` supplies one nontorsion twist section. These facts were already
used by the trace/twist and base-change notes. **The full twist rank is not
proved by the18-section Gram.**

This is now a useful exact fork for theorem hunting:

- Prove `rank T(Q(s))=1`. Then the conic pullback's full rational generic
  rank is exactly18; every subsequent direction is a genuine specialization
  jump relative to this construction.
- Or construct a second independent twist section. It immediately supplies
  a generic rank19 construction and, outside finitely many rational factory
  addresses, two independent seeds without quartic point search.

Neither branch is settled here. Different twist squareclasses or two sheets
alone do not establish additional independence. A new Picard/Frobenius upper
bound or a second exact section would address (5); another scan of M18
landscape scores would not.

**Conjectural interpretation, not a theorem about rank incidence.** The302
runs are consistent with several entry seeds opening an atlas that contains
successively cheap representatives of already-existing additional points.
The sealed non302 rank21 success shows that such amplification is not unique
to302. The similar cheap geometry of null seeds shows why seed construction
and later rank incidence must remain separate. We still lack a prospective
arithmetic condition explaining why302 has at least31 independent points,
and a uniform chart-height theorem turning subgroup distances into the
observed finite-search closure.

## Reproduction and scope

The [arithmetic certificate](../../artifacts/generated-results/elliptic-curves/det1092_seed_geometry_v1.json)
contains the exact Gram decomposition, Sylvester/adjugate height bounds,
conic/deck identities and nonconstant `j` degrees. Its checker imports only
the old generic parent, conic construction and its replay, and the factory
parameter map. The specialization and counting arguments are proofs above,
not conclusions inferred from numerical sampling. This is an ordinary
mathematical proof with an exact CAS check, not formal verification.

```sh
timeout 25s sage -python research/elliptic-curves/cas/verify_det1092_seed_geometry.sage
timeout 25s python3 research/elliptic-curves/cas/check_m18_landscape_comparison.py
```

The new arithmetic check completes in under one second. No point search,
new parameter, class group, descent census, policy change or background job
was started. No later point coordinates or V3 artifacts were construction
inputs. The M18 metadata and terminal outcomes are a separately labelled
retrospective comparison only.
