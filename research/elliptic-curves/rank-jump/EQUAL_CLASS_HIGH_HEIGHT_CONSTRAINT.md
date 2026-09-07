# The remaining high-control section must meet the zero section

The single-prime follow-up **does not close the high control's rank**.
It does yield a useful exact constraint when combined with its rational
component data: every nonzero rational mixed section must meet the zero
section. If such a section is not divisible by two, its height is at least
`9/2`, and its intersection multiplicity with zero belongs to an explicit
square sequence.

The low equal-class control remains at geometric mixed rank zero. Both
claims concern the new fixed-cubic construction, not the rank of either
production curve or its original generic quotient.

## One additional reduction, with a fixed stopping point

The [protocol](EQUAL_CLASS_HIGH_FOLLOWUP_PROTOCOL.json) selects the first
coefficient-eligible prime above 17 and at most 31. Primes 19 and 23 are
ineligible; 29 is eligible. No weighted traces were inspected before this
choice. The one new count gives

| High equal-class surface | Value |
|---|---|
| complementary traces over degrees 1, 2, 3 | `(-25,335,7865)` |
| finite-base contributions | `(-16,358,7811)` |
| smooth-infinity contributions | `(-9,-23,54)` |
| geometric reduction Picard rank | 18 |
| NS discriminant squareclass | `-1` |

The new discriminant agrees with the earlier primes 13 and 17. Agreement
does not prove characteristic-zero Picard rank 18. The high mixed rank
remains `0–1`, and the full new base's arithmetic rank remains `3–4`.
No further prime is added under this protocol.

The [counts](../../artifacts/generated-results/elliptic-curves/rank_jump_equal_class_high_followup_counts_v1.json),
[report](../../artifacts/generated-results/elliptic-curves/rank_jump_equal_class_high_followup_v1.json),
and [independent replay](../../artifacts/generated-results/elliptic-curves/rank_jump_equal_class_high_followup_verification_v1.json)
retain all five reductions and bind the prior certificates. Independent
NumPy sums on 8,613 Frobenius orbits verify every one of the 25,259 new
base parameters. The earlier four arrays and their independent certificate
are reused unchanged. All five Frobenius polynomials and discriminants
are checked by companion matrices. Counting and independent replay have
respective 40- and 60-second caps.

## A conditional height theorem

Use the conventional twist coordinates `(X,Y)` on `E_u^d`, and let
`m=P.O` on its smooth minimal elliptic K3 surface. This is intersection
with the zero section over the entire projective base, including infinity.
For a nonzero rational section, `m` is a nonnegative integer.

The [additive-component theorem](ADDITIVE_COMPONENT_KUMMER_COMPATIBILITY.md)
puts every rational section on the identity component at each of the three
`I0*` fibres. Their height corrections therefore vanish. At the three
`I2` fibres, the Kummer bit `e` gives a common component flag: `e=1`
meets all three nonidentity components, and `e=0` meets all three identity
components. Over a splitting field this is the valuation pattern
`(0,1,1)` and its conjugates: the zero entry labels the single root of the
nodal cubic, while the odd entries label its two colliding roots.
In the local `I2` resolution their odd valuations mean the nonidentity
component. Moreover, `e=0` is equivalent to divisibility by two over
`Q(u)`.

The height formula, with `chi=2`, now reads
\[
h(P)=4+2m-\frac32 e.
\]
The formula and component corrections are standard
[Schütt–Shioda, sections 11.8–11.10](https://arxiv.org/pdf/0907.0298).
In particular there is no nonzero rational torsion: the displayed height
is strictly positive for either bit.

Suppose a nonzero rational section exists. The established upper bound
then forces geometric mixed rank one and Picard rank 18. Equal-rank
specialization to prime 13 (or 17 or 29) forces characteristic-zero
NS discriminant squareclass `-1`.

The trivial lattice has determinant
\[
\det\bigl(U\oplus D_4(-1)^3\oplus A_1(-1)^3\bigr)
=(-1)4^3(-2)^3=512.
\]
The Mordell–Weil discriminant formula, modulo rational squares, gives
\[
[\operatorname{disc}\mathrm{NS}]=[-512h(P)]=[-2h(P)].
\]
This squareclass calculation does not require `P` to be a geometric
primitive generator: its index and the geometric torsion order contribute
only squares. Consequently:

* If `P` is not divisible by two, then
  \[
  5+4m=n^2,\qquad n\text{ odd},\quad n\ge3.
  \]
  Thus `m=1,5,11,19,29,...` and `h(P)=n^2/2 >= 9/2`.
* If `P` is divisible by two, then `m+2` is a square. Repeated halving in
  the finitely generated torsion-free rational group ends in a section
  with `e=1`. Hence `h(P)>=18` and `m>=7`; the apparent value `m=2`
  from the squareclass condition alone cannot occur.

In particular **there is no nonzero rational section disjoint from `O`**.
This does not prove rank zero. Nor does it assert a geometric height bound
over arbitrary constant extensions, where additive component invariants
can change.

## The smallest remaining construction has an exact form

At the first allowed value `m=1`, the intersection with `O` is supported
at a rational base point, possibly infinity. Write it as the zero of a
homogeneous linear form `q` in the base coordinates. Since the Weierstrass
line bundle is `O(2)`, a section at this level has the form
\[
X=A_6/q^2,\qquad Y=B_9/q^3,
\]
where the subscripts denote homogeneous degree, and neither numerator
vanishes at the zero of `q`. For the homogeneous Weierstrass coefficients
`a2,a4,a6` of degrees `4,8,12`, respectively, the exact condition is
\[
\boxed{B_9^2=A_6^3+a_2A_6^2q^2+a_4A_6q^4+a_6q^6.}
\]
Both sides have degree 18. The noncancellation conditions are essential:
they exclude lower-contact or spurious solutions created by clearing
denominators. The component pattern must also be checked on a returned
solution.

This is a finite-degree **rational-solubility target**, available from the
equation and generic points. We have not solved it or run a coefficient
search. Failure at this first allowed height would still leave higher
heights open. A polynomial section with degree bounds `(4,6)` would be
disjoint from zero and is already excluded; an unrestricted polynomial
expression can meet zero at infinity and is not excluded by that wording.

## Mechanism ranking and next decision

1. **Incidence:** an extra geometric divisor remains possible on the high
   construction and is rigorously absent on the low one. Three matching
   reductions do not certify its existence. Further small-prime agreement
   is not the next endpoint.
2. **Solubility:** the remaining rational section, if any, must satisfy
   both additive compatibility and the new height/contact condition.
   The first exact construction target is the degree-18 identity above.
3. **Weak explanation:** common Kummer classes and matching reduction
   discriminants do not force a rational section. No auxiliary-Jacobian
   or Kummer-surface identification has been established here; it cannot
   be inferred from the fibre labels or discriminant squareclass alone.
4. **Visibility consequence:** excluding the disjoint-zero ansatz explains
   why that specific coordinate construction cannot succeed. It is not a
   specialization-rank predictor or a reason to change Agent 1's charts.

Even a successful rational mixed section would supply at most one direction
in this new family, and its image modulo the original rank-16 subgroup
would still need proof. The central multi-direction question therefore
remains the production Selmer block's simultaneous global solubility.
Before investing in a section computation, the next decision must compare
the exact construction above with an independent computation of the retained
cup-obstruction bits. Neither route should be expanded merely because its
current numerical evidence is compatible with success.

Replay from the repository root:

```sh
sage -python elliptic-curves/rank-jump/equal_class_high_followup.py check
sage -python elliptic-curves/rank-jump/equal_class_high_followup.py replay
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_equal_class_high_followup.py
```
