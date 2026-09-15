# Q80 genus-zero branch pairs must collide at131

On the literal direct11952 Q80 MW17 parent, every surviving genus-zero
`k=0` common-quartic solution must have **a double branch reduction at one
of sixteen finite base values modulo131**. This is a restriction on the
actual coefficient system, with no coefficient-height or denominator bound.
It excludes every branch divisor with distinct geometric reductions.

More generally, the conclusion holds over Q_131 for any two rational
anti-invariant sections on a quadratic cover with two smooth geometric
branch fibres, if their specializations are distinct nonzero 2-torsion
points at both branch fibres. No polynomial or height hypothesis is needed
for this implication.

In the fixed Q80 coordinate, choose a primitive integral binary quadratic
`q(T,Z)` for the residual branch divisor. Necessarily

```
q(T,Z) = c*(T-b*Z)^2 mod131,       c !=0,
b in {2,9,12,27,32,35,38,46,48,63,74,75,103,105,107,115}.     (1)
```

In particular its leading coefficient is a131-unit and its discriminant
is divisible by131. The literal cover scalar is still part of the
construction. Equation (1) concerns the branch map in this fixed model;
it does not assert bad reduction of the abstract genus-zero curve.

**This norm-character gate alone leaves all five coefficient allocations
inside the collision residues.** The subsequent
[coefficient proof](Q80_GENUS_ZERO_COLLISION_CLOSURE_2026-09-14.md) closes them.
Genus-one allocations and other parents are outside this theorem. The later
[nodal orientation](Q80_GENUS_ONE_K2_NODAL_ORIENTATION_2026-09-14.md) and
[cubic reciprocity](Q80_GENUS_ONE_K1_CUBIC_RECIPROCITY_2026-09-14.md) proofs
leave five genus-one allocations open. No positive [MW17 two-gain cover](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
is constructed. The earlier same-branch proofs keep their own assurance;
this component does not use their norm-eight singular-pencil replay.

## 1. Coupled norm conditions, including zero ordinates and poles

Use the smooth cubic root curve `S: theta^3+A(t)*theta+B(t)=0` and the
actual seventeen inherited sections `T=(X_T,Y_T)`. The
[single-branch proof](Q80_SINGLE_BRANCH_RECIPROCITY_2026-09-14.md) establishes
that `alpha_T=X_T-theta` has even divisor and that the twist representative
`D*(x-theta)` has odd divisor precisely at the complementary cubic roots
of each nonzero branch specialization.

For an even-divisor function, evaluation modulo squares also makes sense
at a zero or pole: write `alpha=u*z^(2m)` and take the class of u at the
point. Changing the local parameter changes this value by a square.
The norm tame-symbol form of
[Weil reciprocity, Corollary2.5](https://arxiv.org/pdf/1410.5391), then says
that the product of these regularized evaluations over the odd divisor
of the twist representative is a square. This follows directly because
the other tame-symbol factors have even exponents; common support is allowed.

At a smooth cubic fibre with a selected root e, the complementary-root
evaluation has the same squareclass as

```
chi_e(T) = 1                    if T=O,
           X_T-e                if X_T !=e,
           3*e^2+A              if T=(e,0).                         (2)
```

The last case uses `X_T-theta=Y_T^2/unit` at the simple root. If T meets
a different 2-torsion root, the same calculation cancels that zero in
the complementary norm. Poles have square leading coefficient in the
elliptic formal parameter, giving the first case. The weighted infinity
change multiplies abscissas by a fourth power.

Let b run through the closed branch points and e_i(b) be the root selected
by the i-th anti-invariant section. Each branch cubic splits over its
residue field, since two distinct roots are selected. For every inherited T,

```
product_b Norm_{k(b)/Q_131}(chi_{e_i(b)}(T)) is square,    i=1,2.    (3)
```

These are simultaneous conditions for both sections. A single character
at one branch point need not vanish separately.

## 2. Complete residue calculation

The parent has good reduction at131, with squarefree degree24 discriminant
and no singular fibre over P1(F_131). The root curve has good reduction.
The even-divisor classes extend there: their horizontal valuations are even,
and their vertical valuation is zero after square scaling. Equivalently,
the local square/unit calculations above extend across the model. Thus their
squareclasses specialize, even when their original numerator or denominator
vanishes at the reduced point.

Write `F_131^2=F_131[i]`, `i^2=-1`. Encode an element `a+b*i` as `a+131*b`.
For all132 rational projective base values and all8515 conjugate pairs
of nonrational quadratic base values, the certificate retains every cubic
root over this field, with multiplicity, and all seventeen character bits.
Quadratic characters are computed after the norm to F_131. For rational
roots the separate F_131 characters are also retained.

| Residue calculation | Complete result |
|---|---:|
| Rational fibres with three rational roots |16|
| Rational fibres with one rational root and two quadratic roots |62|
| Rational fibres with no quadratic-field root |54|
| Ordered pairs of distinct rational roots in fully split fibres |96|
| Collisions between their ordered pairs of17-bit character vectors |0|
| Nonrational quadratic base orbits with no cubic root |2820|
| Nonrational quadratic base orbits with exactly one cubic root |4238|
| Smooth nonrational quadratic base orbits with three roots |1456|
| Nodal nonrational quadratic base orbits |1|
| Quadratic root triples with two zero character vectors, outside the16 rational split fibres |0|

The nodal orbit is retained. At `t=100+31*i` its root multiset is

```
(122+38*i, 122+38*i, 18+55*i),
character vectors (13412,13412,0).                              (4)
```

An inherited section cannot meet the node, by the local total-space
equation `xy=t`. Its evaluation at the double root is consequently a
nonzero unit, or the square pole case. At the simple root the character
is zero: the other two roots coincide and their product evaluation is a
square. This handles a smooth characteristic-zero branch fibre that
reduces to the nodal fibre. Omitting this orbit would leave a gap.

## 3. Why every local branch algebra forces (1)

There are only two geometric branch points. Over Q_131 their algebra is
either split or quadratic.

If it is split, both branch cubics split over Q_131. Their smooth rational
reductions must lie in the sixteen fully split fibres. For each section,
(3) says its character vector at the first point equals its vector at
the second. Injectivity of the96 ordered vector pairs forces the same
base reduction and the same ordered root choices. Infinity is included
in the roster and is not fully split.

If the branch field is unramified quadratic and the base parameter has
nonrational residue, one of the8515 rows applies. Each of the two selected
roots would need zero norm-character vector. None of the smooth split
rows allows this. The nodal row (4) does not either: the only zero vector
belongs to the simple root, which has a unique nearby lift, so two
distinct selected roots cannot both reduce there.

An unramified quadratic branch field may instead have rational base
residue. A rational one-root fibre splits over the quadratic residue
field, but its root-character vectors have only one zero vector; it is
again the unique rational simple root. Such a fibre cannot receive two
distinct selected roots. A root-free rational fibre has no root over
the quadratic residue field. Thus only the sixteen fully split fibres
remain, and the two conjugate branch parameters have the same reduction.

For a ramified quadratic branch field, every residue lies in F_131.
All rational parent fibres are smooth, so distinct roots over the branch
field have distinct residues. Its cubic must therefore split completely
over F_131, at one of the same sixteen values. Both conjugate branch
parameters reduce there. These cases exhaust quadratic local algebras and
prove (1), including reduction collisions and arbitrary denominators.

## 4. The collision residues really survive these norm tests

The conclusion is not a local nonexistence theorem for the remaining
coefficient systems. For example, consider the rational branch pair

```
b_1=2,   b_2=133,
q(t)=(t-2)*(t-133)=t^2-135*t+266.
```

Both parent cubics have simple roots reducing to19,54,58. Hensel lifting
pairs these roots over Q_131. Each inherited character is locally constant
on the corresponding root residue disc, so its two values have square
product. Choosing two distinct root labels satisfies all conditions (3).
The same argument works for arbitrarily close distinct branch values in
each of the sixteen discs. Raising the modulus of these character tests
alone cannot remove this boundary.

This supplies only local branch data satisfying necessary conditions.
It supplies no two section functions, no global cubic splitting, and no
solution of the fixed-parent coefficient equations. Although this example's
conic has rational points, that alone does not make it a two-gain source.
The [subsequent arithmetic gate](Q80_GENUS_ZERO_COLLISION_CLOSURE_2026-09-14.md)
uses the coefficient equations near the collision and closes the polynomial
genus-zero cases. This example remains valid necessary branch data.
Genus-one degree-two/three/four norm systems remain separate problems.

## 5. Reproducible evidence

The [input](../artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1/input.json)
binds the literal parent, the previous exact17-section checker, and both
new implementations. The [Sage producer](scripts/certify_q80_degree_two_reciprocity.sage)
uses finite-field factorization and square testing. The
[independent checker](scripts/verify_q80_degree_two_reciprocity.py) implements
integer arithmetic in F_131[i]. It proves roster completeness by full
three-root factorization, a nonsquare quadratic discriminant, or
`gcd(x^(131^2)-x,f)=1`; it checks every character and the96-pair injection.
Derivative evaluation cancels numerator/denominator zeros independently
of the producer's rational-function simplification.

The previous exact rational identities and good-reduction certificate are
also replayed. Six new tests cover extension norms, omitted roots, nodal
multiplicity, 2-torsion regularization, poles/cancellation and false sections.
The [independent receipt](../artifacts/generated-results/elkies-k3-q80-degree-two-reciprocity-v1/independent-replay.json)
records1.594 seconds; producer CPU was9.570 seconds. Each process has a
60-second CPU/4GiB cap. The discovery preflight is retained separately.
Written reciprocity and specialization arguments are not formally verified.

```
python3 research/elkies-k3/scripts/verify_q80_degree_two_reciprocity.py
python3 -m unittest discover -s research/tests -p 'test_q80_degree_two_reciprocity.py' -q
```
