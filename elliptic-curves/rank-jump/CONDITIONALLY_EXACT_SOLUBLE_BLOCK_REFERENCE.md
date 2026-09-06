# A complete +6 reference for the missing additional strict block

The supplemental fibre **a1-fibration-05 at t=3/17** now has an exact
local-boundary calculation using only its equation and sixteen generic
sections. Its generic strict dimension is **6**, its inherited minus-one
CT-switch rank is **6**, and its generic first-pair governing field again
has degree **192**. The generic subgroup fills the entire Selmer boundary.

Joining the repository's completed rank/Selmer certificate afterward proves,
**under that certificate's class-character GRH assumption**, that the full
strict Selmer group has dimension **12** and is entirely rational. Exactly
six dimensions remain after quotienting its generic strict subspace. This
is a fully accounted-for +6 reference, not an exact zero-gain control or a
new member of the frozen matched panel.

The [comparison artifact](../../artifacts/generated-results/elliptic-curves/rank_jump_bounded_gain_reference_comparison_v1.json)
separates the masked measurements from that conditional inference. Its
additional class representatives and additional CT matrix were **not**
computed independently of exceptional points.

## Inputs and authority

The arithmetic [input](../../artifacts/generated-results/elliptic-curves/rank_jump_bounded_gain_reference_inputs_v1.json)
contains one token, its model, sixteen marked generic sections, and prime
hints for the equation discriminant. The generic projection comes from
the retained retrospective input. The hints come only from the equation's
factorization field in the completed rank-22 proof; every prime and the
full product are verified anew. No exceptional coordinate, point-derived
class, character, anchor, relation or rank label enters a worker.

The separate outcome comes from
[the canonical completion proof](../notes/SMALL_CONDUCTOR_CLASS_COMPLETION_PROOF_2026-09-06.md)
and [certificate](../../artifacts/generated-results/elliptic-curves/small_conductor_class_completion_v1.json),
under `EC-SMALL-CONDUCTOR-CLASS-TARGET-20260906` in MATH_STATUS.json.
That proof uses all twenty-two known rational directions, including
exceptional points, and proves

\[
 \dim\operatorname{Cl}(K)/2=16,\qquad
 \operatorname{rank}E(\mathbb Q)=\dim\operatorname{Sel}_2(E)=22
\]

conditional on GRH for the nontrivial quadratic ordinary ideal-class
characters specified there. The rank lower bound 22 is unconditional.
The current work checks the source binding and joins its stated theorem;
it does not rerun the extensive class-group completion proof.

The field discriminant agrees exactly between the masked computation and
that certificate:

```text
128900477062442043600727490102612931938219670661531295245188752203875468
```

The arithmetic processes read no outcome manifest. Selection remains
retrospective, and this compact MW16-05 reference is not a same-family or
same-scale match to the MW16-01 +11 example. No frozen panel row is replaced.

## Exact boundary and conditional block accounting

Use the [strict-block definitions and exact sequence](FRESH_STRICT_BLOCK_NECESSITIES.md):
G is the generic Kummer subgroup, U the strict Selmer group, I the full
localized Selmer image, and

\[
 m=\dim G,\quad k=\dim(G\cap U),\quad
 c_S=\dim U=\dim\operatorname{Cl}(\mathcal O_{K,S_K})/2.
\]

Here S contains infinity, 2 and all bad primes. The masked arithmetic gives

\[
 m=16,\quad k=6,\quad \dim\operatorname{loc}(G)=10,\quad
 \ell=\dim\prod_{v\in S}L_v=11.
\]

The equation-defined derivative class −disc(f)f′(theta) has norm disc(f)^4
and real sign pattern (−,+,−). Its pairing with the real point-image
generator (+,−,−) is nonzero. Reciprocity therefore gives dim I≤ell−1=10.
The generic image already has dimension ten, so **dim I=10 exactly** and

\[
 I=\operatorname{loc}(G),\qquad
 \boxed{\operatorname{Sel}_2(E)/G\simeq U/(G\cap U),\quad
        \dim\operatorname{Sel}_2(E)=c_S+10.}
\tag{1}
\]

Equation (1) is unconditional and uses no exceptional point. It means that
every additional class admits a strict representative after adding a generic
class; it does not assert that each originally supplied representative is
itself strict.

The rank lower bound 22, joined afterward, forces **at least six additional
strict rational directions unconditionally**. The conditional upper bound
closes the remaining intervals:

\[
 c_S=12,\qquad \dim U/(G\cap U)=6,\qquad
 \Sha(E/\mathbb Q)[2]=0
 \quad\text{under the stated GRH assumption.}
\tag{2}
\]

The irreducible cubic gives E(Q)[2]=0. Equality of rank and 2-Selmer dimension
then makes every 2-Selmer class rational. Thus all twelve strict directions
are soluble, and exactly six contribute beyond G. Equation (2) also implies
that the bad-prime ideal classes span dimension four in the sixteen-dimensional
ordinary class quotient. Both class-dimension conclusions are **label-derived
conditional deductions**, not independently measured incidence features.

This determines the dimension of the additional arithmetic quotient. It
does not identify its minimal simultaneous-solubility variety, its genus or
degree, or a condition on varying t that produces six rational lifts.

## The inherited switch and comparison

Four entries of the initial Artin matrix had nonunit residues. Their gcds
were 269,269,269,41. The separate bounded completion evaluates those good-prime
Frobenius contributions with local-power tests and independent generic point
characters. All four agree. The resulting inherited form is

\[
 \mathrm{CT}_{E^{(-1)}}|_{G\cap U}=A+A^\mathsf T
 \simeq H\oplus H\oplus H,
\]

of rank six. Original CT vanishes on this inherited subspace because its
classes come from generic rational sections. None of the nonzero classes
in this transported six-dimensional subspace can be rational on the minus-one
twist. This is a simultaneous obstruction switch on a **generic** block;
it does not give a rank bound on the whole twist.

All rows below retain the frozen panel's labels; a reported zero is censored.
Each has its own cubic field. The last row is the new supplemental reference.

| Fibre | Reported gain | Generic strict k | Inherited CT-switch rank | Additional boundary cap a | Additional strict rational dimension |
|---|---:|---:|---:|---:|---:|
| 074d9, 2824/885 | 0 | 9 | 8 | 0 | unknown |
| 103b2, −1049/2296 | 0 | 0 | 0 | 2 | unknown |
| 103b2, 3726/881 | +10 | 0 | 0 | 1 | ≥9 |
| 11952, −1171/1683 | 0 | 10 | 10 | 0 | unknown |
| 11952, −2448/11 | +10 | 0 | 0 | 5 | ≥5 |
| 11952, 130349/28916 | 0 | 5 | 4 | 0 | unknown |
| 11952, 110314/102227 | +10 | 0 | 0 | 4 | ≥6 |
| ICARM356, generic17 | +12 | 1 | 0 | 1 | ≥11 |
| ICARM385, generic17 | +12 | 0 | 0 | 8 | ≥4 |
| ICARM398, generic16 | +14 | 0 | 0 | 9 | ≥5 |
| MW16-05, 3/17 | exactly +6 under GRH | 6 | 6 | 0 | exactly 6 under GRH; ≥6 unconditional |

Every fixed generic first-pair governing field has degree192, including the
reference. It supplies no separation. Nor does a large inherited switch:
larger switches occur on observed-zero controls than on this exact +6 fibre,
while several +10 fibres have no inherited strict block. These observations
reject those quantities as sufficient explanations of the retained jumps;
they do not establish population-level independence or a predictive law.

The fresh +11 boundary calculation remains incomplete because its frozen
equation factorization is incomplete. This supplemental result does not fill
that entry or identify any additional-class CT matrix in the fresh panel.

## Lessons and next falsifiable target

1. **Incidence, highest priority:** obtain an unramified additional class
   supply independently of exceptional points. This reference now has a
   precise calibration endpoint: a complete masked calculation must recover
   c_S=12, and an additional quotient of dimension six, under the same GRH
   assumptions if it uses a conditional upper bound. Returning only the six
   generic strict classes measures an incomplete extractor, not a small jump.
   The earlier norm-projection dictionary's ramified survivors illustrate
   why square norm, independence and bad-place tests do not suffice.
2. **Solubility, still open prospectively:** explain why an independently
   constructed strict excess is rational. Here the completed outcome theorem
   settles solubility after the fact. It does not produce a point-independent
   sufficient condition in t. A nonzero CT form would obstruct solubility;
   zero CT alone still does not prove it on another fibre.
3. **Weak explanations:** generic pair field degree, inherited strict size,
   and inherited CT-switch rank have now been calibrated against both
   censored controls and a conditionally exact positive-gain fibre. No
   candidate-scoring change follows. Local boundary saturation locates the
   additional block but does not predict that it exists or is rational.
4. **Useful alternative:** the completed class proof supplies a reusable
   character-exclusion method which can certify a proposed ideal-class span
   before a full relation presentation is found. Its retained point-derived
   anchors and adaptively selected relation waves must not be imported as a
   masked feature. Removing point coordinates from an adaptively selected
   dataset would not erase its selection dependence.

A small falsifiable next experiment is an **input-provenance and capacity
test**, before any new relation campaign: use this single frozen equation
and its generic sections to propose an unramified-class construction; freeze
its finite generator set and bounds before evaluating it. Measure its exact
additional span modulo G, recording all outside-S valuation coordinates.
The success endpoint is six independent strict classes beyond G; any smaller
span disproves completeness of that construction on this known-soluble
reference. No new elliptic points, family sweep, or use of the six exceptional
classes is required. Only a method which supplies such independent classes
is ready for the original fixed103b2 high/low CT discrimination test.

There is no new visibility measurement. For Agent1 the immediately usable
information is the necessary incidence accounting (1), the failure of the
inherited proxies, and a reference on which class-extraction coverage can be
tested. No prospective selector is claimed.

## Reproduction

The [initial protocol](BOUNDED_GAIN_REFERENCE_PROTOCOL.json) caps factor-hint
verification at15 seconds, local/Artin work at60, and boundary/octic stages
at30 each. The [completion protocol](BOUNDED_GAIN_REFERENCE_COMPLETION_PROTOCOL.json)
allows30 seconds for the four small repairs. There is one worker, no full
class-group computation, no general factorization and no point search.

The [verification certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_bounded_gain_reference_verification_v1.json)
passes sixteen generic independence checks,160 local-character checks,
114 separate local-power checks, six ideal-square identities,32 direct
Jacobi entries, four repaired Artin entries and nine independent modular
octic replays. It also checks the derivative norm and exact boundary, and
audits the conditional label join. It explicitly does not replay the full
external class-group proof.

```sh
timeout 60 sage -python elliptic-curves/rank-jump/verify_bounded_gain_reference.py check
```

The two captures preserve the partial and completed results separately.
All artifacts bind their inputs and code. This work changes no active search
file, worker protocol, candidate population or mathematical-status entry.
