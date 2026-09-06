# Parameter covers have room for large blocks, but capacity is not incidence

A quadratic cover of the original R17 parameter line branched at two
smooth fibres can have at most **20 new geometric generic directions**;
four smooth branch points raise the bound to **22**. Thus the
[one-auxiliary-point multiplicity bound](ONE_AUXILIARY_POINT_HAS_A_MULTIPLICITY_BOUND.md)
does not rule out a large block arising from a finite parameter cover.

The new coefficient calculation proves more: the published R17 surface
discriminant is irreducible of degree **24** over \(\mathbb Q\). Any
rational branch divisor of degree below 24 therefore avoids every
singular fibre. These low-degree covers necessarily lie in the
smooth-branch cases above.

The 37 retained quadratic covers and 11 retrospectively fitted quartic
covers all pass the exact geometry checks. Their capacity is uniform;
it does not distinguish high-gain controls from the lower-gain control.
The canonical discriminant quadratic cover, which is defined before
any points, has no rational lift at any of seven retained R17 parameters.

## Generic rank accounting on a quadratic parameter cover

Let \(E/\mathbb Q(u)\) be a semistable elliptic K3 surface with geometric
reducible-fibre root rank \(R\). Let \(B\to\mathbb P^1_u\) be a connected
quadratic cover with \(b>0\) geometric branch points, of which \(r\)
lie at multiplicative fibres. Thus \(b\) is even and
\[
g(B)=b/2-1.
\]
Write the quadratic extension as adjoining \(\sqrt{d(u)}\). Decomposition
into invariant and anti-invariant sections gives
\[
\operatorname{rank}E(\overline{\mathbb Q}(B))
=\operatorname{rank}E(\overline{\mathbb Q}(u))
 +\operatorname{rank}E^{(d)}(\overline{\mathbb Q}(u)).
\]
The same decomposition holds arithmetically over \(\mathbb Q\) for a
cover defined over \(\mathbb Q\).

Twisting a smooth fibre produces \(I_0^*\): Euler number increases by six
and root rank by four. Twisting \(I_n\) produces \(I_n^*\): Euler number
again increases by six, while root rank increases from \(n-1\) to \(n+4\),
a gain of five. Consequently
\[
\chi(E^{(d)})=2+b/2,\qquad
R(E^{(d)})=R+4b+r.
\]
For an elliptic surface over \(\mathbb P^1\) in characteristic zero,
\(h^{1,1}=10\chi\). The Hodge bound and Shioda--Tate yield
\[
\boxed{
\operatorname{rank}E^{(d)}(\overline{\mathbb Q}(u))
\le 10(2+b/2)-2-(R+4b+r)=18-R+b-r.}
\]
The local twist rules and the elliptic-surface formulas used here are
standard; see [Schütt--Shioda, *Elliptic surfaces*, §§5–6](https://arxiv.org/pdf/0907.0298).

For the retained Picard-rank-19 families this gives:

| Original fibration | Original geometric generic rank | Smooth branch points | Cover genus | Twist rank upper bound | Pullback generic rank upper bound |
|---|---:|---:|---:|---:|---:|
| R17, \(24I_1\) | 17 | 2 | 0 | 20 | 37 |
| R17, \(24I_1\) | 17 | 4 | 1 | 22 | 39 |
| A1/MW16, \(I_2+22I_1\) | 16 | 2 | 0 | 19 | 35 |
| A1/MW16, \(I_2+22I_1\) | 16 | 4 | 1 | 21 | 37 |

The A1 rows use the existing geometric fibration data in
[the A1/MW16 atlas](../../elkies-k3/ICARM_A1_MW16_ATLAS_2026-09-04.md);
no new A1 cover is constructed here. A branch point on a multiplicative
fibre reduces the corresponding twist bound by one. The artifact's
general table includes such conditional geometric configurations; it
does not assert that every configuration occurs over \(\mathbb Q\).

These are **generic rank bounds for the pullback and twist**, not upper
bounds on ranks of individual rational fibres. Even a fibre admitting a
lift to \(B\) can gain further directions beyond the generic pullback.
Likewise, a bound of 20 is not evidence that 20 new sections exist.

## The original R17 singular fibres form one degree-24 orbit

The input retains the original published coefficient polynomials \(A,B\).
Exact arithmetic verifies
\[
\Delta=-16(4A^3+27B^2),\qquad
(\deg A,\deg B,\deg\Delta)=(8,12,24),
\]
with \(\Delta\) squarefree and coprime to \(A\). Thus all 24 finite
singular fibres are \(I_1\); the model is smooth at infinity.

Rational polynomial factorization returns one degree-24 factor. An
independent finite-polynomial certificate proves its irreducibility
without trusting that factorization:

| Prime | Irreducible factor degrees modulo the prime | Remaining possible proper rational factor degrees |
|---|---|---|
| 167 | \(2,2,6,14\) | \(2,4,6,8,14,16,18,20,22\) |
| 181 | \(1,6,17\) | \(6,18\) |
| 191 | \(8,16\) | none |

After removing integer content, the reductions have full degree.
Any proper rational factor would, by Gauss's lemma, reduce to a product
of some modular irreducible factors at each prime. Its degree would
belong to all three subset-sum sets, which have empty intersection.
The pure-Python verifier checks the factor products and every modular
factor's irreducibility with Rabin's criterion.

Therefore the singular divisor is a single closed point of degree 24
on the rational parameter line. A branch divisor defined over
\(\mathbb Q\) containing one of its geometric points must contain its
entire orbit. In particular, for \(b<24\), necessarily \(r=0\).
This includes branching at infinity, which is a smooth fibre.

The first modular dictionary through 97 had only bad displayed
reductions, and its independent irreducibility result remains
**UNKNOWN** in its original artifact. A separate bounded completion
through 503 closed at 191. Bad reduction in that dictionary was not
treated as evidence of reducibility.

## Exact retrospective cover checks

The [capacity protocol](PARAMETER_COVER_CAPACITY_PROTOCOL.json) freezes
all 37 quadratic branch polynomials already retained in the five-fibre
census, and the eleven quartics already fitted to the rank-28 control.
No new parameter or point is searched.

For each polynomial \(q\), the calculation verifies separability,
\(\gcd(q,\Delta)=1\), and the actual twisted invariants
\[
A_q=q^2A,\qquad B_q=q^3B,\qquad
\Delta_q=q^6\Delta,\qquad c_{4,q}=q^2c_4.
\]
Degree bounds give a smooth fibre at infinity and the finite configurations
\[
24I_1+2I_0^*\quad\text{or}\quad24I_1+4I_0^*.
\]
Their root ranks are 8 and 16; their \(\chi\)'s are 3 and 4.
The resulting bounds are 20 and 22, respectively.

The high/low distinction is absent from this capacity:

| Retained R17 parameter | Observed quotient over MW17 | Retained quadratic hits | Per-cover twist upper bound |
|---|---:|---:|---:|
| \(-2/377\) | 8 | 6 | 20 |
| \(-308/251\) | 9 | 3 | 20 |
| \(2456/135\) | 10 | 2 | 20 |
| \(-9529/5471\) | 11 | 1 | 20 |
| \(3/8\) | 4 | 25 | 20 |

The existing
[branch and specialization analysis](BRANCH_BLOCKS_AND_SPECIALIZATION.md)
proves that the 25 lifts in the last row span only four exceptional
directions. This computation neither changes that result nor turns the
upper bound 20 into a lower bound. Each known anti-invariant section
retains only its already proved contribution; hidden additional sections
on these twists remain uncomputed.

All eleven quartics at \(-9529/5471\) have capacity 22 individually.
Their separate fitted covers are not one common 11-dimensional section
block, and they use exceptional points retrospectively. The same
family-level bounds also apply at the recent R17 high/observed-zero
parameters; no cover construction for that pair is asserted here.

## The canonical discriminant cover misses the retained fibres

A natural coefficient-defined candidate is
\[
B_\Delta:\quad v^2=\Delta(u).
\]
It branches at all 24 singular fibres, so \(b=r=24\), has genus 11,
and its quadratic twist has
\[
\chi=14,\quad R=120,\quad
\operatorname{rank}E^{(\Delta)}(\overline{\mathbb Q}(u))\le18.
\]
Its large branch divisor spends most of its Hodge capacity on reducible
fibres. Large cover genus is not a guarantee of many new sections.

A rational lift at \(u_0\) requires \(\Delta(u_0)\) to be a rational
square. Exact rational-square tests fail at all seven retained parameters:
the five in the table above and the recent pair
\[
-2300/843,\qquad -1561/3133.
\]
Thus this canonical cover cannot explain their exceptional directions by
generic sections evaluated at a rational lift. The arithmetic exclusion
is specific to this cover; it does not exclude other quadratic covers
or high rank on those fibres.

For an irreducible specialized 2-division cubic, its discriminant square
condition would reduce the Galois group from \(S_3\) to \(A_3\).
The failure agrees with the retained \(S_3\) controls. This confirms that
that particular Galois-sign event is not needed for their observed jumps.

## What the missing block would have to prove

For a two-smooth-branch-point R17 twist,
\[
\operatorname{rank}E^{(d)}(\overline{\mathbb Q}(u))=\rho(E^{(d)})-10.
\]
A block of fourteen generic geometric directions therefore requires
Picard rank at least 24, within the Hodge ceiling of 30. Those section
directions must additionally be rational over \(\mathbb Q(u)\) to
explain rational points by this construction.

The desired chain is now:

1. **Incidence:** an original-family choice of \(d(u)\), made without
   exceptional points, forces a sufficiently large rational section
   space of the twist.
2. **Solubility:** \(d(u_0)\) is square and the explicit sections evaluate
   at a rational lift of the specialization.
3. **Independence:** those evaluated directions remain independent
   modulo the original generic subgroup.

The capacity theorem supplies only a necessary upper bound for step 1.
The retained examples supply isolated sections and retrospective lifts,
not a large common section space. Neither a high Picard rank without
Galois control nor a square branch value alone closes the chain.

## Ranked implications

1. **Viable mechanism:** a finite quadratic parameter cover is not
   geometrically too small to carry a \(+8,+10\), or \(+14\) block.
   This remains a substantive opening distinct from the one-point
   auxiliary-curve mechanism.
2. **Excluded candidate:** the canonical discriminant cover has no
   rational lift on the seven tested R17 fibres. Its Galois-sign
   condition is not their explanation.
3. **Weak feature:** low-degree smooth branching is automatic for
   rational R17 branch divisors of degree below 24. Its rank capacity
   cannot distinguish the matched controls.
4. **Missing theorem/computation:** certify a large rational twist
   section space for a point-independent \(d(u)\), then certify its
   specialized quotient image. No enlarged bisection search is justified
   merely by unused capacity.
5. **For Agent 1:** these are structural **incidence** upper bounds and
   a failed **solubility** gate, not a new rank score. **Visibility**
   policies and candidate populations stay unchanged.

## Evidence and replay

- [Frozen model and branch inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_parameter_cover_capacity_inputs_v1.json)
- [Twist geometry and capacity table](../../artifacts/generated-results/elliptic-curves/rank_jump_parameter_cover_capacity_v1.json)
- [Initial verification and seven discriminant-cover tests](../../artifacts/generated-results/elliptic-curves/rank_jump_parameter_cover_capacity_verification_v1.json)
- [Bounded modular irreducibility completion](../../artifacts/generated-results/elliptic-curves/rank_jump_surface_discriminant_irreducibility_v1.json)
- [CAS-free three-prime certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_surface_discriminant_modular_verification_v1.json)

An initial capture reached JSON serialization with a Sage integer and
failed; its [failure record](../../artifacts/generated-results/elliptic-curves/rank_jump_parameter_cover_capture_failure_v1.json)
and partial local checkpoint are preserved. The successful capture
uses explicit integer conversion and a separate checkpoint.

    sage -python elliptic-curves/rank-jump/parameter_cover_capacity.py check
    sage -python elliptic-curves/rank-jump/verify_parameter_cover_capacity.py check
    sage -python elliptic-curves/rank-jump/surface_discriminant_irreducibility.py check
    python3 elliptic-curves/rank-jump/verify_surface_discriminant_modular.py check

Only new rank-jump files and the analysis index are changed. No active
search, model, mathematical-status entry, or K3 proof source is modified.
