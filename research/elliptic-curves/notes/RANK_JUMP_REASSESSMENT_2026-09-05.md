# Rank jumps: what the search is actually testing

Reviewed against `18ab9abf830923d9aee35ea0e2d2a32d0628fb85`.
This is a structural audit with a written generic-family calculation and cheap
exact diagnostics. It does not change the existing specialization, Selmer or
rank-32 entries in `MATH_STATUS.json`, rerun their expensive certificates, or
claim a new general rank-jump theorem. The underlying surface/descent theory
is established; the application to our fixed-field pencil is the deduction
below. Historical protocols and artifacts remain unchanged.

## What the failed experiments do and do not say

The [MW18 comparison](MW18_DEEP_CENTRE_CALIBRATION_2026-09-05.md) found
22 presentation-weighted extra directions with deep centres versus zero with
nearest-first. All 600 charts completed. Its predeclared gate failed, but the
experiment had no adaptive wave. Two of its five presentations are the same
curve. This is useful detection evidence, not a failed existence theorem.

The [height-population pilot](FIBRE_HEIGHT_POPULATION_2026-09-05.md) found no
new points in either arm. Its improved MW16 models still have roughly 940-bit
largest minimal coefficients, whereas the calibration examples have 184--289
bits. MW18 used singleton centres, not the new deep-centre policy. Eight curves
per arm and zeros in both arms do not distinguish rare jumps from inadequate
exposure. They do not establish that height is irrelevant.

The [fixed-field comparison](FIXED_FIELD_COMPARISON_2026-09-05.md) is different:
it proves that many transported classes are insoluble. All five new odd-dimensional
spaces have maximal restricted alternating rank. The remaining one-dimensional
radicals are not rational points and have no whole-curve rank interpretation.

The strongest mistake was to treat these three outcomes as measurements of the
same missing phenomenon. They measure **candidate incidence**, **coordinate
visibility**, and **global solubility**, respectively.

## The fixed-field family has generic rank zero

Let A,B be rational, B nonzero, and delta=-4A^3-27B^2 nonzero. Our family is

\[
 E_u:y^2=F_u(x),\quad
 F_u(x)=x^3+2Au x^2+(A+3Bu+A^2u^2)x+B+ABu^2-B^2u^3.
\]

Put D(u)=1+Au^2+Bu^3. Then

\[
 \boxed{\operatorname{rank}E/\overline{\mathbf Q}(u)=1,\qquad
 \operatorname{rank}E/\mathbf Q(u)=
 \begin{cases}1&B\in\mathbf Q^{\times2},\\0&B\notin\mathbf Q^{\times2}.
 \end{cases}}
\]

This is a statement about the function-field curve, **not** an upper bound on
any specialized E_u(Q). In particular it does not contradict the twenty
independent points at u=0 or impose rank one on a tested deformation.

### Proof and exact checks

Direct calculation gives

\[
 \Delta(E_u)=16\delta D(u)^2,\qquad
 c_4(E_u)=16(A^2u^2-9Bu-3A),
\]

Also \(\operatorname{disc}(D)=\delta\). For A nonzero,

\[
 \operatorname{Res}_u(c_4,D)=4096(4A^3+27B^2)^2.
\]

When A=0 the degree of c4 drops: c4=-144Bu is still coprime to D,
since D(0)=1. The code checks coprimality directly in both cases.
Consequently there are precisely three geometric finite singular fibres, all
I2. At infinity set v=1/u, X=v^2x, Y=v^3y. The integral equation is

\[
 Y^2=X^3+2AvX^2+(A^2v^2+3Bv^3+Av^4)X-B^2v^3+ABv^4+Bv^6.
\]

Its discriminant has valuation six and its c4 valuation at least two. The
coefficient of v^3 in a6 is nonzero. In characteristic zero this is a minimal
I0* fibre. Equivalently the three distinct leading root slopes are the squares
of the roots of x^3+Ax+B; they remain distinct since B and delta are nonzero.

Thus the minimal surface is geometrically rational, with chi=1, geometric
Picard rank ten and reducible-fibre root lattice 3A1+D4. Shioda--Tate gives
10-2-(3+4)=1. The rational-surface criterion, rank formula and intersection
height formula used here are in [Schuett--Shioda, sections 8 and 11][SS].

Over Q(sqrt(B))(u) there is an explicit section

\[
 S=(Bu^2,\sqrt B D(u)),\qquad F_u(Bu^2)=B D(u)^2.
\]

It is disjoint from O. At each finite singular fibre,

\[
 \partial_xF_u(Bu^2)=(A+3Bu)D(u),
\]

so S goes through the node of the Weierstrass fibre and meets the nonidentity
component of the resolved I2. At infinity its reduction is
(X,Y)=(B,B sqrt(B)), a smooth point of Y^2=X^3, so it meets the identity
component there. Its Shioda height is therefore

\[
 \langle S,S\rangle=2\chi-3(1/2)=1/2>0.
\]

S spans the one-dimensional geometric Mordell--Weil vector space. Constant-field
Galois sends it to itself or its negative according to the character of B.
The invariant subspace is consequently one-dimensional for square B and zero
otherwise. No full class group, point search or specialization heuristic is
used in this argument. We do not need to claim S generates the integral group.

The pinned anchor has

```
B = 167347710468055045100164888198438918505621536951206
floor(sqrt(B)) = 12936294309733952752041790
```

and the two adjacent squares strictly bracket B. Hence **our particular
fixed-field pencil has arithmetic generic rank zero**. Its rank-at-least-20
u=0 fibre is already an exceptional specialization of this rank-zero pencil.
The prior language suggesting transport of a generic rank-20 structure was
not justified: we transported cohomology classes, not sections.

The Sage-free `fixed_cubic_geometry.py` recomputes the identities as polynomial
coefficient identities over Q[u], checks separability and coprimality, and emits
the proof's precise scope. B=0 or delta=0 is rejected rather than extrapolated.

### A useful second identity, not a high-rank claim

\[
 F_u(-Au-1/u)=-D(u)^2/u^3.
\]

Thus the base change u=-v^2, v nonzero, has the rational point

\[
 (x,y)=(Av^2+v^{-2},(1+Av^4-Bv^6)/v^3).
\]

At v=1 this is exactly the already known point (A+1,A-B+1). This explains a
point by an identity, not by hoping that a large Selmer space is soluble.
It does not show that the point belongs to the inherited twenty-class span,
or that several independent inherited classes have been transported.

There is also a coordinate issue. Under x=c^2x', y=c^3y' on the anchor,
the same generator construction has u'=c^2u. An integer box |u|<=3 is not an
isomorphism-invariant notion of being close to the high-rank fibre u=0.
This is checked by an exact scaling regression. Re-scaling the box is a new
population experiment, not evidence of restored rank transport.

## Jump size is not an intrinsic quality of a curve

The calculation above embeds any separable short cubic with nonsquare B into
a generic-rank-zero family through that same curve. One can manufacture a very
large *displayed jump* retrospectively without explaining how to find the curve.
Absolute rank is intrinsic; the decomposition into generic rank plus jump is
relative to a chosen family and subgroup.

Likewise, if a degree-d base change makes k formerly exceptional directions
generic at a rational lift of the same fibre, it changes the accounting

\[
 r+j=(r+k)+(j-k),
\]

not the fibre's rank. The image of a nontrivial rational base map is a thin set,
and height growth changes the affordable population. Higher generic rank is
not automatically a better record search. The MW18 anchors have ten remaining
demonstrated directions; fourteen is the target increment, not an observed
transfer. This must stay explicit in schedulers and comparison reports.

A rational point on one fibre is not automatically a new divisor class of the
fixed K3 surface. A bounded atlas of low-degree multisections need not account
for all exceptional points. Conversely, extending one point to a multisection
only supplies that certified section after base change; it does not extend an
entire exceptional subgroup.

## The overlooked link between generic sections and cubic class pressure

[Gillibert--Levin, Theorem 2.1][GL] injects a suitable mod-p section subgroup
into Pic(C)[p], where C compactifies the nonzero p-torsion cover. Component-group
corrections are part of the theorem, not optional bookkeeping. For p=2 this is
a trigonal curve. In the geometrically connected 24-I1 K3 setting, its 24 simple
branch points give, by Riemann--Hurwitz, 2g-2=-6+24, hence g=10. With no component
corrections the generic MW17 already contributes seventeen dimensions to the
geometric construction in that theorem.

This is a reason to expect large inherited cubic class information before any
exceptional points are found, not a reverse implication producing new points.
Applying the paper's number-field specialization theorem additionally requires
its hypotheses, including the relevant universal-bad-reduction check. We do
not assert those extra hypotheses for a new family from a discriminant alone.

The existing full/localized class-image distinctions and saturation corrections
in [the soluble-cover theory](../../elkies-k3/RATIONAL_SOLUBILITY_AND_RESIDUAL_SELMER_THEOREMS.md)
remain necessary. No finite genus bound here bounds specialized cubic class
ranks or the ranks of E_t(Q).

## Why preserving E[2] did not preserve the points

An E[2] identification identifies the ambient H1, but not the local Kummer
images, the global rational-point image, or the higher descent data. The
[Selmer-companion theorem of Mazur--Rubin][MR] has additional torsion and local
compatibility hypotheses. Even an isomorphism of Selmer groups would not by
itself identify their rationally soluble subspaces. We must not call a family
rank-preserving just because its cubic algebra is fixed.

The new CT profiles are consistent with this loss of structure. As an
**illustrative algebraic null**, a uniformly random alternating F2 matrix of
odd dimension 2m+1 has maximal-rank probability

\[
 (2-2^{-2m})\prod_{i=1}^m(1-2^{-(2i-1)}),
\]

about 0.839 in dimensions 13, 15 and 17. A one-dimensional radical is forced
by parity and is common in this model; it is not a special positive signal.
`alternating_rank_distribution` checks the full distribution by an exact
bordering recurrence, with exhaustive matrix tests through dimension five.
The five elliptic curves are NOT independent uniform matrices: no p-value,
actual frequency law or arithmetic conclusion is deduced from this analogy.
The broader random alternating-matrix connection is developed in [BKLPR][BKLPR].

The soluble classes form a subspace. Two insoluble basis classes can have a
soluble sum. Rational classes lie in the full CT radical, not merely in a
maximal isotropic subspace; a restricted radical may still be obstructed by
external classes or higher-divisible Sha. These are different losses of
information and need separate labels.

## Replace blind budget escalation with an observability audit

For a known point R, its raw pointed coordinate is

\[
 t_Q(R)=(y_R+y_Q)/(x_R-x_Q).
\]

For the retained horizontal map t=(as+b)/(cs+d), compute exactly

\[
 s=[dT-bU:-cT+aU]\quad\text{when }t=[T:U].
\]

Primitive normalization then gives the minimum affine search height for this
point. Handle the tangent value at R=-Q, the two known endpoints and parameter
infinity explicitly. There is no need to enumerate a single rational slope.

`search_observability.py` independently recomputes the final quartic identity
from the full retained transformation chain. It distinguishes:

- a point outside the chosen box;
- a point in an unsearched/timed-out interval;
- a point in completed coverage but absent from the output;
- a returned point, whose independence is a separate campaign certificate.

This is a retrospective diagnostic. It does not know the height of an unseen
point, and a missed basis representative does not prove its quotient direction
was missed: another representative or combination might have been recovered.
A transcript's completeness still trusts its pinned enumeration worker.

For **prospective exposure controls**, mask one or more already certified
independent generic points on each candidate and try to recover them from the
remaining subgroup. This places positive controls on the actual large-coefficient
fibres, instead of extrapolating from small record curves. `masked_control`
extracts the principal Gram block and writes the withheld oracle separately.
The original independence certificate remains required. Recovery measures
WITHHELD_KNOWN_DIRECTIONS, not a newly increased rank. Its geometry also changes
when the reference subgroup is reduced, so it is a useful diagnostic, not a
complete calibration theorem for exceptional directions.

Only after the search policy is frozen may the oracle be opened to explain
misses. Keep coefficient/coordinate selection and public points separate.
Chart-counts and trillions of numerator/denominator pairs are not independent
rank-jump trials; the population unit is a distinct curve with recorded exposure.

## The research order changes

**First, diagnose rather than rescan.** Use retained MW18 charts and known
points to locate the missing coverage. Separately perform masked generic-point
controls on the already prepared ordinary fibres. Do not move failed thresholds
retroactively or confuse the initial-only trial with an adaptive protocol.

**Second, validate selection on the right population.** Compare modest,
frozen prime-prefix Nagao/Mestre policies with disjoint validation primes,
reporting actual curves and exposure. Elkies's original account explicitly
uses thousands of primes, many millions of specializations and shared residue
trace tables; it also explains the half-lattice/fake-descent point search [E07].
Our bounded short-prime screens are not automatically that selector at a smaller
cost. There is no guarantee that extending the prime bound will help; test it
on a held-out population rather than naming another record-sized search.

**Third, require point witnesses from a construction hypothesis.** Search a
well-defined incidence problem (a candidate fibre together with a certified
new point or independent pullback section), not just a large ambient Selmer
space. Small-point-conditioned construction and staged descent have precedents
in [Watkins et al.][W]. Keep the gain and height cost of each additional section
explicit. The existing seventeen rank-zero product twists stay excluded; their
exclusion does not close every V4 base change or the existing rank-19 routes.

The non-thin rank-jump theorems [LS] and [GS] are important existence results under
specific surface/multisection hypotheses. They are not positive-density
predictions for our finite boxes, and do not promise jumps of fourteen on
these Picard-rank-nineteen K3 models. A new family should declare which
hypotheses and which point-producing map it uses before a large campaign.

## Reproduction and changes

No existing historical certificate, frozen selection rule, or mathematical
status entry is rewritten. These independent diagnostics can be run before
allocating new arithmetic work:

```sh
python3 -m unittest discover -s elliptic-curves/tests -p test_fixed_cubic_geometry.py
python3 -m unittest discover -s elliptic-curves/tests -p test_search_observability.py
python3 -m unittest discover -s elliptic-curves/tests -p test_rank_jump_audit_cli.py
python3 elliptic-curves/cas/audit_rank_jump.py family \
  --input artifacts/generated-results/elliptic-curves/fixed_cubic_field_fermigier_rank20_local_kummer_u2_v1.json
python3 elliptic-curves/cas/audit_rank_jump.py ct-null --dimension 17
```

`visibility --record CHART.json --oracle POINTS.json` consumes a shared search
record (or its checkpoint envelope) and a list of known points. `mask --input
INPUT.json --withhold 0,2 --search-input BLIND.json --oracle ORACLE.json` splits
an input with `curve`, `points`, `metric_gram`; it does not launch a search.
Output files must be new. These commands neither require nor certify a full
Selmer group. The written generic-rank proof and its hypotheses must accompany
any downstream use; specialized ranks remain unchanged.

The new tests are ordinary Python tests. Full Sage, PARI, remote-conic and
historical certificate campaigns are separate checks; this patch does not
claim to have rerun them. The prior report's ten status/checker-source hash
mismatches are not repaired by silently repinning old proofs.

## Primary sources and provenance

[SS]: https://arxiv.org/abs/0907.0298
[GL]: https://arxiv.org/abs/1811.08166
[MR]: https://arxiv.org/abs/1203.0620
[BKLPR]: https://arxiv.org/abs/1304.3971
[E07]: https://arxiv.org/abs/0709.2908
[W]: https://doi.org/10.5802/pmb.9
[LS]: https://arxiv.org/abs/1907.01987
[GS]: https://arxiv.org/abs/2505.15159

SS: Schuett and Shioda, *Elliptic surfaces*, sections 8, 11.
GL: Gillibert and Levin, *Elliptic surfaces over P1 and large class groups of
number fields*, Theorems 2.1, 2.7 (their hypotheses matter).
MR: Mazur and Rubin, *Selmer companion curves*, Theorem 3.1.
BKLPR: Bhargava, Kane, Lenstra, Poonen and Rains, *Modeling the distribution
of ranks, Selmer groups, and Shafarevich--Tate groups of elliptic curves*.
E07: Elkies, *Three lectures on elliptic surfaces and curves of high rank*,
Lecture III; the half-lattice search is established prior art.
W: Watkins, Donnelly, Elkies, Fisher, Granville and Rogers, *Ranks of quadratic
twists of elliptic curves*, especially the small-point and descent strategy.
LS: Loughran and Salgado, *Rank jumps on elliptic surfaces and the Hilbert property*.
GS: Garbagnati and Salgado, *Rank jumps and multisections of elliptic fibrations
on K3 surfaces*; the headline small-Picard-number theorem is not our hypothesis.
