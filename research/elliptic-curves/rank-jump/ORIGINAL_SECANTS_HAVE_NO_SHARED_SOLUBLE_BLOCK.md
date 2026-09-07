# Original generic secants do not expose the observed jump blocks

The complete fixed dictionary of **1,504 generic-basis secants** gives no
rational intercept points and no repeated quadratic class on any of the
six retained high/low controls. Their quadratic characters span at least
19--21 dimensions. This rules out a small common quadratic cover for this
dictionary at these fibres, including the three known high-gain fibres.

Unlike [the preceding fitted-pencil experiment](SECANT_PENCILS_DO_NOT_EXPLAIN_THE_ORIGINAL_JUMP.md),
this calculation keeps the original elliptic curve. It tests a concrete
possible solubility event using only marked generic points. Its negative
result does not exclude other multisections, other choices of generic
basis, or exceptional points outside this dictionary.

## Mechanism and falsifiable test

For the original short model \(E:y^2=F(x)=x^3+Ax+B\), take two marked
generic points \(P_i=(a,p)\), \(P_j=(b,q)\). For each relative sign
\(\epsilon\in\{1,-1\}\), their secant has x-axis intercept

\[
x_{ij,\epsilon}=\frac{a\epsilon q-bp}{\epsilon q-p},
\qquad C_{ij,\epsilon}=F(x_{ij,\epsilon}).
\]

The potential point

\[
R_{ij,\epsilon}=(x_{ij,\epsilon},\sqrt{C_{ij,\epsilon}})
\]

lies on the original curve. A square \(C\) is an exact rational-solubility
event. Two nonsquare values whose ratio is square give points over the
same quadratic extension, equivalently rational points on one quadratic
twist. This is a possible shared arithmetic object; independence and
non-torsion would still require proof. A rational hit would likewise need
independence modulo the original generic subgroup before counting toward
the jump.

The [protocol](ORIGINAL_SECANT_CLASSES_PROTOCOL.json) freezes every
unordered pair of marked generic basis points with both relative signs:
\(m(m-1)\) secants on a rank-\(m\) marked family. It introduces no new
parameters and uses no exceptional points. The computation input contains
only six indexed short models and their generic point lists. Outcome and
family labels are joined afterwards by the verifier. This is a
retrospective diagnostic, not a prospective validation cohort.

All denominators and slopes are nonzero, all intercepts are finite, and
all \(C\)'s are nonzero in these six cases. Thus there are no degenerate
rows hiding missing outcomes.

## Three paired comparisons

Here \(m\) denotes the original generic rank and \(w\) the retained
certified independent specialized subgroup rank. The observed quotient
rank is \(q=w-m\). The zero controls have observed quotient zero; their
full curve rank remains unknown. No value of \(w\) below is used as a
full curve-rank upper bound.

| Original family | Parameter | \(m,w,q\) | Secants | Rational hits | Repeated classes | Character rank, lower bound |
|---|---:|---:|---:|---:|---:|---:|
| A1/MW16-05 | 307/206 | 16,25,9 | 240 | 0 | 0 | 19 |
| A1/MW16-05 | -3158/1291 | 16,16,0 observed | 240 | 0 | 0 | 20 |
| A1/MW16-04 | -1647/91 | 16,25,9 | 240 | 0 | 0 | 21 |
| A1/MW16-04 | -2177/2397 | 16,16,0 observed | 240 | 0 | 0 | 21 |
| published R17 | -2300/843 | 17,24,7 | 272 | 0 | 0 | 21 |
| published R17 | -1561/3133 | 17,17,0 observed | 272 | 0 | 0 | 21 |

Thus no pair supports the hypothesis that high gain comes from these
generic secant characters becoming rational or merging into one shared
quadratic class. The rank bounds use a dictionary with only 21 bits and
are not exact character ranks. In particular, the 19 versus 20 comparison
is not evidence for a lower true character rank on the high-gain fibre.

## Exact certificates without factorization

For nonzero \(C\in\mathbb Q\), the producer records its sign and the
valuation parity and unit-square bit at each of

\[
3,5,7,11,13,17,19,23,29,31.
\]

These give a homomorphism from rational squareclasses to
\(\mathbb F_2^{21}\). Different signatures certify different global
classes. Equal signatures are only a filter collision: the producer
checks their exact rational ratio with integer square roots. There are
110 such ratio checks in total, and none gives a square.

An independent verifier uses Sage rational arithmetic, recomputes every
intercept by the determinant formula, and tests **every pair of values**
without the local filter: 188,432 exact nonsquare ratio checks. It also
tests all 1,504 values for rational solubility independently of the
producer's integer-square-root implementation.

Six \(C\)'s have zero signature in the limited local dictionary but are
still globally nonsquare: respectively 1,1,0,3,1,0 across the six rows.
This is a concrete example of why a finite local filter cannot stand in
for global solubility. These are equations \(z^2=C\), not everywhere
locally soluble torsors or examples of Sha; no full local-solubility
claim is made for those six values.

For each curve the verifier retains indices of 19, 20, or 21 independent
signature rows and checks their rank over \(\mathbb F_2\). The corresponding
rational squareclasses are therefore independent. Their multiquadratic
compositum has degree at least \(2^{19}\), \(2^{20}\), or \(2^{21}\),
respectively. This is a field-degree lower bound, not a Mordell--Weil
rank statement. These quadratic fields are fields of the constructed
intercept points, not the halving fields of known exceptional points.

## A specialization certificate on the original family

There is a useful consequence beyond a null count. Let \(C_1(u),\ldots,C_k(u)\)
be the same constructions from the marked generic sections of the
original family over \(\mathbb Q(u)\). At any retained fibre where they
are regular and nonzero, specialization defines a homomorphism on their
squareclass span. Indeed, if

\[
\prod_i C_i(u)^{e_i}=h(u)^2,
\]

then the left side is a unit at the specialization place, so \(h\) is a
unit there too. Evaluating produces the same square relation over
\(\mathbb Q\). Consequently:

- A specialized nonsquare proves that the corresponding function is
  not a square over the original parameter field.
- Two distinct specialized classes prove that the corresponding
  functions do not share a squareclass generically.
- Independent specialized classes prove independence of the
  corresponding generic characters.

Applied to the retained marked sections, the experiment certifies that
the 240 or 272 functions in each such dictionary are generically
nonsquare and pairwise distinct. It certifies a generic character rank
of at least 20 for A1/MW16-05 and at least 21 for A1/MW16-04 and published
R17, taking the stronger bound in each pair. This deduction uses the
panel's existing identification of the inputs as original generic
sections; it does not reconstruct their formulas or assert exact
generic twist ranks.

This rules out a generic identity equating two of these squareclasses.
It does **not** rule out collisions or square values at other special
parameters. The stronger practical finding is that neither event occurs
at the three tested high-gain parameters where such an explanation was
needed.

The statements depend on the marked basis and the prescribed secants.
They are invariant under rational rescaling between short Weierstrass
models: \(x=c^2x'\), \(y=c^3y'\) sends an intercept to \(c^2x'_0\)
and multiplies \(C'\) by the square \(c^6\). A basis change can produce
different secants and is not covered by this negative result.

## Consequences for the mechanism ranking

1. **Still the strongest constructive model — simultaneous solubility:**
   one quadratic condition can produce multiple independent directions
   in the previously proved split-cubic example. Its production
   analogue must use a more specific original-family construction than
   this entire generic-basis secant dictionary.
2. **Still the strongest production evidence — Selmer incidence:** the
   existing strict ideal-class and CT structures describe arithmetic
   blocks that actually contain observed exceptional classes. Their
   remaining gap is global rational solubility, as demonstrated by the
   [rational-versus-Sha controls](NORM_LIFTS_CAN_BE_ENTIRELY_SHA.md).
3. **Downgraded:** fitted secant pencils and original generic-basis
   intercept square conditions both fail their paired production tests.
   Enlarging this secant dictionary without a new structural identity
   is not justified by the present evidence.
4. **Missing implication:** a point-free original-family condition must
   select a block of the relevant exceptional Selmer classes and prove
   simultaneous rational solubility, with independent images modulo the
   original generic subgroup. The current calculation supplies neither
   implication and prevents crediting this particular construction with
   doing so.
5. **For Agent 1:** no new candidate score follows. These square tests
   are solubility tests for prescribed points, and character sharing is
   only a candidate common-cover structure. Neither is a rank predictor;
   chart recovery remains a separate visibility endpoint.

The next useful experiment should return to the observed strict blocks:
test whether their explicit 2-cover equations admit a common low-degree
geometric construction before solving for rational points. An exact
geometric factorization or quotient relation would distinguish that work
from increasing a visibility/search budget. No such new campaign is
launched in this change.

## Artifacts and replay

- [Generic-only inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_original_secant_class_inputs_v1.json)
- [All 1,504 constructions and exact classes](../../artifacts/generated-results/elliptic-curves/rank_jump_original_secant_classes_v1.json)
- [Independent exhaustive verification and rank witnesses](../../artifacts/generated-results/elliptic-curves/rank_jump_original_secant_class_verification_v1.json)

```sh
python3 elliptic-curves/rank-jump/original_secant_classes.py check
sage -python elliptic-curves/rank-jump/verify_original_secant_classes.py check
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_original_secant_classes.py
```

The positive regression recovers the rational intercept of the previous
oblique example; separate tests cover squareclass multiplication, a
finite-local-filter false positive, and degenerate secants. Capture is
checkpointed per curve with a 30-second worker limit. These new files and
certificates are independent of Agent 1's outputs; no live search or
mathematical-status entry is changed.
