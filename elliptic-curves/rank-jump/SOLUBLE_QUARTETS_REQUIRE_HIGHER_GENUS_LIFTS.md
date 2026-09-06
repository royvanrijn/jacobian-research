# The observed four-cover blocks require a genus-17 rational lift

**Follow-up:** [exact Frobenius moments](NATIVE_SINGLE_COVER_CANNOT_EXPLAIN_THE_WHOLE_PLUS8.md)
now bound the full arithmetic generic rank of the native 1795d twist by 7.
Even all hidden sections on that one cover cannot supply the whole observed
+8 quotient. Other common-cover and product-character possibilities remain
open.

**Solubility:** neither successful four-cover block can come from a
nonconstant rational or elliptic parametrization of those fixed cover
equations. Both simultaneous-lift curves have **genus 17**, exact rational
and geometric **gonality 8**, and minimum geometric degree **4** to a
genus-one curve. The previously obstructed quartet has the same genus and
map-degree bounds. These are proved restrictions on the proposed construction, not
rank estimates or explanations of why the successful parameters occur.

For the observed +8 fibre, the missing implication can now be located
precisely: a positive-rank, globally soluble pair carrier must lift through
a **ramified degree-four cover of genus 17** to make the other two covers
soluble. There is no rational section of this map. The genus-one carrier
has infinitely many rational points, while only finitely many can lift.

A second bounded audit checks the alternative of several displayed
sections sharing one quadratic field. All **39,119** frozen finite-chart
atlas equations define different quadratic extensions of Q(t). The atlas
therefore supplies no repeated-field group to test for multiple independent
sections. Additional, undisplayed sections on a fixed twist remain an open
possibility.

## Fixed systems and exact branch checks

The [compression protocol](SOLUBLE_QUARTET_COMPRESSION_PROTOCOL.json) fixes
the two previously successful quartets and the existing obstructed control:

| System | Atlas labels, with orbit- prefix omitted | Observed quotient over the marked MW17 subgroup | Carrier genus | Gonality |
|---|---|---:|---:|---:|
| 08234-003 | 01333, 0b2d0, 13109, 19e45 | +7 | 17 | 8 |
| 08234-009 | 0911e, 0a037, 1795d, 18f5d | +8 | 17 | 8 |
| Obstructed A,B,C,D | 030cb, 03da0, 07086, 11278 | Not a fibre observation | 17 | 8 |

These are comparisons between fixed cover systems on the same R17 family.
The obstructed quartet is not a newly chosen low-rank curve. The two high
fibres' four constructed points have a certified quotient rank **at least
3 and at most 4**; neither this calculation nor their four square values
prove rank four. Full ranks of the original fibres remain UNKNOWN.

For each quartet, exact polynomial gcds verify that the four quadratics
are squarefree with pairwise disjoint geometric branch divisors. All have
degree two. Consequently they are independent in
\(\overline{\mathbf Q}(t)^*/\overline{\mathbf Q}(t)^{*2}\): a branch
point belonging to any included factor has odd valuation in a nonempty
product. No branch at infinity occurs. The same assertions hold over Q.

This proves that
\[
C_n:\quad u_i^2=f_i(t)\quad(1\le i\le n)
\]
is geometrically connected of degree \(2^n\) over the parameter line.
At each of its \(2n\) branch points the inertia group has order two.
Riemann–Hurwitz gives
\[
2g(C_n)-2=-2\cdot2^n+2n\cdot2^{n-1},\qquad
g(C_n)=1+2^{n-1}(n-2).
\]
Thus one, two, three and four independent native radicands give genera
0, 1, 5 and 17. Their character quotients
\(H_S:y^2=\prod_{i\in S}f_i(t)\) have genus \(|S|-1\).
For a quartet these comprise four genus-zero, six genus-one, four
genus-two and one genus-three quotient. The positive-genus dimensions sum
to 17. A small quotient's genus does not equal the simultaneous carrier's
genus.

## Exact lower bounds on compression maps

Work first over an algebraic closure. Forgetting any one root gives a
degree-two map \(C_n\to C_{n-1}\). Suppose there is also a degree-d map
\(C_n\to B\), with B of genus h equal to zero or one.
If this map does not factor through that root-forgetting map, the two
function subfields generate Qbar(C_n), since the root extension has prime
degree two. The
[Castelnuovo–Severi inequality, Proposition 3.2](https://math.mit.edu/~poonen/papers/dyn_gonality.pdf)
then yields
\[
g(C_n)\le2g(C_{n-1})+dh+(2-1)(d-1),
\quad\text{hence}\quad d(h+1)\ge2^{n-1}.
\]

If d is below that threshold, the map must factor through **every**
root-forgetting map. Its functions are then invariant under all n
independent sign involutions, so lie in Qbar(t). For h=0 its degree is a
multiple of \(2^n\), contradicting the assumed bound. For h=1 this would
give a nonconstant map from P1 to a genus-one curve, which is impossible.

The bounds are attained: projection to a single genus-zero conic has
degree \(2^{n-1}\), and projection to any genus-one pair carrier has
degree \(2^{n-2}\). Each tested single conic also has an explicitly
verified retained rational point, so its identification with P1 is defined
over Q and the gonality upper bound holds over Q too.

| Number of radicands | Genus | Exact geometric gonality | Minimum geometric degree to genus one |
|---:|---:|---:|---:|
| 2 | 1 | 2 | 1 |
| 3 | 5 | 4 | 2 |
| 4 | 17 | 8 | 4 |

These maps go **from** the simultaneous carrier to a simpler curve.
They do not parametrize its rational points. In the other direction,
Riemann–Hurwitz implies that a degree-e map from a curve B onto C4
requires \(g(B)\ge16e+1\). In particular no rational or elliptic source
can supply all four roots as rational functions of a varying parameter.
A constant map to the already observed point is not a construction of
new specializations.

## The positive +8 comparison leaves a ramified lifting problem

Use the [already certified soluble pair](GLOBAL_CARRIER_SOLUBILITY_AND_SPECIALIZATION.md)
\(C_{1795d,0911e}\) inside the +8 quartet. It has a rational point and its
Jacobian has exact rank 2. The remaining roots have labels 0a037 and
18f5d. On the pair carrier, each remaining quadratic has eight simple
geometric zeros: its two roots in t each have four unramified preimages.
Its poles have even orders. The two zero divisors are disjoint.

Adjoining one remaining root therefore gives a double cover branched at
eight points and genus 5. Adjoining both gives a degree-four cover of the
genus-one pair carrier, with 16 branch points of inertia two, and genus
17. This is the exact additional simultaneous-solubility problem.
Positive auxiliary rank and vanishing Sha[2] of the *pair's Jacobian*
do not solve these ramified lifts.

The pair carrier has infinitely many rational points. The genus-17
carrier has finitely many by
[Faltings' finiteness theorem](https://math.uchicago.edu/~drinfeld/Deligne%27s_conjecture_Manin_conf/Faltings_argument/Faltings.pdf).
Since the parameter maps are finite, infinitely many rational parameters
admit the two pair lifts but fail to admit all four. This proves that even
the **two individual square conditions** cannot globally replace the
quartet. No new parameter or rational point was generated to establish
that conclusion. The number or full list of quartet points is not known.

## A product-square shortcut looks perfect in the small sample

The calculation evaluates all 15 product characters on all 32 frozen
parameters, for each of the three systems: 1,440 exact tests. At a fixed
t away from the branch locus, let
\[
\rho_t:\mathbf F_2^4\longrightarrow
\mathbf Q^*/\mathbf Q^{*2},\qquad
e_i\longmapsto[f_i(t)].
\]
This is the squareclass map of **parameter-cover radicands**, not a
Selmer map for the original elliptic fibre. The kernel is computed by
testing every product for being a square.

| System | rank rho = 0 | rank rho = 2 | rank rho = 3 | rank rho = 4 |
|---|---:|---:|---:|---:|
| +7 quartet | 1 | 0 | 0 | 31 |
| +8 quartet | 1 | 0 | 0 | 31 |
| Obstructed quartet | 0 | 1 | 2 | 29 |

For each positive quartet, every nonempty product is a square at the
one successful parameter, and none is a square at any other sampled
parameter. Thus *every one-character shortcut*, including the product of
all four quadratics, has perfect agreement with full lifting on this
sample. This does not establish equivalence: the covers were selected
retrospectively, and the sample contains no intermediate splitting for
those quartets.

Generically a one-character quotient forgets three of the four sign
directions. Its function field has degree two over Q(t), whereas the
simultaneous carrier has degree sixteen; the remaining lift has degree
eight. Whether the product-of-four genus-three quotient has additional
rational parameters outside the full quartet image remains UNKNOWN.
The infinite-counterexample argument above applies to the positive-rank
pair carrier, not automatically to this genus-three quotient.

## Does one existing cover already carry several displayed sections?

The [separate atlas protocol](ATLAS_COMMON_COVER_PROTOCOL.json) tests a
different way around the genus growth. If several independent directions
are sections over the **same** quadratic function field, one square
condition can make them rational together, without forming C4.

For squarefree quadratic polynomials, equality of squareclasses in
Q(t) forces identical branch divisors, hence proportional polynomials;
the proportionality constant must itself be a rational square. The audit
groups by primitive signed integer coefficients and tests scalar ratios
only within matching groups. An independent check normalizes by dividing
by the leading coefficient instead.

There are 39,119 equations, 39,119 different branch polynomials and
39,119 different quadratic fields. There are **zero** repeated-field
groups. This uses the full frozen finite-chart atlas, with no exceptional
points, ranks or parameter selection. The sole previously excluded inverted
chart remains outside the scope. The grouping does not prove each twist
has rank one: it only proves that this displayed section roster supplies
no second candidate over the same quadratic field.

## Ranked conclusions and next proof gates

1. **Solubility, established:** the observed quartet is a rational point
   on an intrinsically genus-17 carrier. The pair-to-quartet implication
   requires a ramified degree-four lift and fails for infinitely many
   rational pair parameters. Locating its exceptional soluble lifts is
   the remaining arithmetic problem for these fixed labels.
2. **Incidence leading to solubility, still viable:** several independent
   sections on one common quadratic cover would avoid this genus cost.
   The current atlas supplies no such group. A bounded section computation
   on one fixed native twist, with a generic independence certificate,
   would test for undisplayed directions. That is a different computation
   from enlarging parameter or chart searches.
3. **Solubility, possible but unproved:** vary the cover labels in a
   mathematically defined construction, so several lifts arise from one
   auxiliary event without insisting on these four fixed radicands.
   No such positive construction is supplied here.
4. **Disproved for these systems:** a rational or elliptic parametrization
   supplying all four fixed roots; or replacing the quartet everywhere
   by the known positive-rank pair condition.
5. **Weak evidence:** perfect product-square agreement in this 32-fibre
   retrospective sample. It cannot justify a prospective selector or
   resolve the uncovered original quotient directions.
6. **Visibility:** no conclusion about chart exposure, heights or recovery
   budgets follows. Those are separate from these existence obstructions.

Agent 1 could eventually use a proved common-cover section block as a
sufficient construction rule. The present result supplies constraints on
that rule, not a new scoring feature or a reason to exclude high-rank
fibres. The active search and mathematical-status registry are unchanged.

## Certificates and replay

```sh
python3 elliptic-curves/rank-jump/soluble_quartet_compression.py check
python3 elliptic-curves/rank-jump/atlas_common_cover.py check
sage -python elliptic-curves/rank-jump/verify_quartet_compression.py check
```

The producer uses rational polynomial arithmetic and exhaustive square
characters on retained parameters. The independent verifier uses Sage
polynomials, homogeneous integer square tests and F2 linear algebra.
The geometric theorems and their hypotheses are proved above; they are
not inferred from a green arithmetic check.

- [Frozen quartet inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_soluble_quartet_compression_inputs_v1.json)
- [Geometry and 1,440 character tests](../../artifacts/generated-results/elliptic-curves/rank_jump_soluble_quartet_compression_v1.json)
- [Full finite-atlas common-field audit](../../artifacts/generated-results/elliptic-curves/rank_jump_atlas_common_cover_v1.json)
- [Independent arithmetic verification](../../artifacts/generated-results/elliptic-curves/rank_jump_quartet_compression_verification_v1.json)
