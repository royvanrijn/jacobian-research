# Q80 rational-bisection branch fields exclude all same-branch genus-zero cases

For every smooth rational bisection of the direct11952 alternate-Q80 MW17
parent, the smooth parent fibre over its quadratic branch field has
**no nonzero rational 2-torsion**. The finite certificate covers all
39,147 translation classes. Consequently the anti-invariant height lattice
of Q(C)-rational sections on each cover has minimum **exactly 12**, and none
of these covers supports a polynomial common-quartic section of height eight.

Combined with the [singular-pencil gate](Q80_GENUS_ZERO_SINGULAR_GATE_2026-09-14.md),
this excludes **all three same-branch genus-zero cases `k=2`**, including
the formerly open cross heights `+/-2`. The subsequent
[reciprocity gate](Q80_SINGLE_BRANCH_RECIPROCITY_2026-09-14.md) excludes
all four `k=1` cases. The later
[coefficient collision proof](Q80_GENUS_ZERO_COLLISION_CLOSURE_2026-09-14.md)
closes all five genus-zero `k=0` cases. The remaining
genus-one strata, other parents and the positive
[infinite rank-at-least-19 construction](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
remain **OPEN**. No exact twist-rank or whole-family genus bound follows.

## 1. The finite arithmetic theorem

Write the fixed parent as `E_t:y^2=x^3+A(t)x+B(t)`. The
[complete retained rational-bisection atlas](R17_NORM12_BISECTION_CHARACTER_EXHAUSTION_2026-09-03.md#3-complete-alternate-frame)
gives 39,147 branch quadratics `q_i(t)`, all irreducible over Q, with branch
divisors disjoint from the parent's 24 nodal fibres. These are all smooth
rational bisections modulo translation by inherited rational sections.
Translation preserves the map to the original t-line and its branch field.

For each `K_i=Q[t]/(q_i)`, with `alpha_i=t mod q_i`, the new certificate proves

```
x^3+A(alpha_i)*x+B(alpha_i) has no root in K_i,
E_alpha_i(K_i)[2] = {O}.                              (1)
```

Here is the exact witness criterion. Normalize q_i to a primitive integer
quadratic. At a prime p require:

1. all coefficients of A and B are p-integral;
2. the leading coefficient and discriminant of q_i are p-units;
3. a simple root `t0 in F_p` of q_i;
4. the smooth cubic `x^3+A(t0)x+B(t0)` has no root in F_p.

The root t0 defines a degree-one prime of K_i over p. Its alpha_i is
integral there. A hypothetical K_i root x of the monic cubic is also
integral there: a negative valuation would make x^3 the unique term of
least valuation. Reduction would give an F_p root, contradicting item4.
Thus a single witness excludes every algebraic root in that fixed branch
field, without a rational-height bound. Different fields may use different
primes. A branch polynomial with irreducible reduction, a repeated root,
leading-coefficient loss or a nonintegral parent is deferred at that prime.

The [frozen input](../artifacts/generated-results/elkies-k3-q80-rational-branch-torsion-v1/input.json)
uses the existing compact atlas, the fixed parent coefficients and the first
64 primes at least101. The first40 primes, through313, close every field;
the other24 are unused. The final survivor is excluded at313. Each prime
has its own checkpoint, including the three stages that add no exclusion.
No trace compilation, neighbour enumeration, bisection search, point search
or exceptional-specialization selection was performed.

The [result](../artifacts/generated-results/elkies-k3-q80-rational-branch-torsion-v1/result.json)
records 39,147 distinct index/residue witnesses and complete coverage. The
producer enumerates F_p roots of the cubic. The independent checker instead
computes `gcd(x^p-x,x^3+a*x+b)` by polynomial modular exponentiation. It
checks every used fibre roster, every simple branch root, coefficient
integrality, exact quadratic irreducibility over Q and attachment to the
retained orbit masks and frame words. The constant squareclass is immaterial
to (1), because only the branch field is used; it remains mandatory when
accepting a positive quadratic-cover construction.

## 2. The anti-invariant minimum is twelve

Let C be the normalization of any retained smooth rational bisection and
`L=Q(C)`. Its degree-two map to the t-line has precisely the two conjugate
branch values of q_i. A ramification point of C has residue field K_i,
independently of the nonzero scalar in its equation. For a Q(C)-rational
anti-invariant section P, specialization commutes with the involution and
gives `P=-P`. Equation (1) forces its value to be O.

Over an algebraic closure there are two ramification points. Any nonzero
section P therefore has `P.O>=2`. The pullback has `chi=4` and only
irreducible fibres, so the standard Shioda formula gives

```
height(P)=8+2*(P.O)>=12.                              (2)
```

This uses the intersection normalization of the height pairing from
[Schutt--Shioda, section11.8](https://arxiv.org/pdf/0907.0298).
The bisection's existing anti-trace has height12 by the retained atlas,
so the minimum is exactly12. This does not bound the rank or assert that
all height12 vectors are multiples of that anti-trace. In particular,
additional independent sections of height at least12 remain possible.

A section `(x(t),s(t)*sqrt(q))` in the full polynomial genus-zero chart
misses O everywhere, including infinity in the weighted coordinates.
It has height8 and is incompatible with (2). Thus a successful Q80
common-quartic genus-zero construction must use a branch field outside
this complete smooth rational-bisection atlas.

## 3. Closing the two remaining same-branch cases

Use the complete twelve-stratum chart and height formulas from the
[genus-zero gate](Q80_GENUS_ZERO_SINGULAR_GATE_2026-09-14.md#1-the-full-genus-zero-coefficient-chart).
Suppose `k=2`, so at both branch values P1 and P2 choose the same nonzero
2-torsion point. The inherited Q80 unramified Kummer theorem supplies one
rational parent section T and actual points

```
R_plus=(P1+P2+T)/2,       R_minus=(P1-P2+T)/2,
R_plus+sigma(R_plus)=R_minus+sigma(R_minus)=T.           (3)
```

This is actual divisibility, not formal division of a lattice word. The
[previous descent proof](COMMON_QUARTIC_BRANCH_STRATA_2026-09-14.md#3-the-unramified-arithmetic-input-on-q80)
checks the good, nodal and infinity places and the required arithmetic
Brauer input. Its proof applies to the actual quadratic branch divisor.

For `j=0` or2 the heights of the anti-traces in (3) are12 and20 in
some order. The bisection images have

```
B^2=(height(R-sigma(R))-16)/2,
p_a(B)=1+B^2/2.                                       (4)
```

Thus the image corresponding to height12 has square -2 and arithmetic
genus zero. Its anti-trace is nonzero, so it is a genuine integral
bisection with normalization C. An integral arithmetic-genus-zero curve
on a smooth K3 is smooth and geometrically rational. After inherited
translation it belongs to the complete atlas of section1, so its quadratic
cover is one of those excluded for polynomial height8 sections by (2).
This contradicts the original P1 (and also P2).

Equivalently, at a ramification point equation (3) gives
`R_plus-R_minus=P2`, a nonzero K_i-rational 2-torsion point, directly
contradicting (1). The equality follows from doubling and the torsion-free
pullback; all points extend over the smooth branch fibre. This formulation
does not rely on an affine chord denominator being nonzero there.

The middle case `j=1` has two anti-traces of height16 and would give an
integral bisection of arithmetic genus1 with rational normalization. The
previous complete singular-pencil gate excludes it. Therefore every `k=2`
case is empty, including higher contacts, degree drops and infinity. No
rational point on the covering conic was needed for these exclusions.

## 4. Replay boundary and remaining construction gate

The [producer](scripts/certify_q80_rational_branch_torsion.py) used1.338
component CPU seconds within the120-second/4GiB cap. The
[independent checker](scripts/verify_q80_rational_branch_torsion.py) took
1.714 wall seconds; its [receipt](../artifacts/generated-results/elkies-k3-q80-rational-branch-torsion-v1/independent-replay.json)
records all39,147 exclusions. These timings exclude authoring and are
not an end-to-end research-speed comparison.

```
python3 research/elkies-k3/scripts/verify_q80_rational_branch_torsion.py
python3 -m unittest discover -s research/tests -p 'test_q80_rational_branch_torsion.py' -q
```

Five regressions retain a valid split-prime witness and reject actual
residue torsion, the wrong branch root, ramification, leading-coefficient
loss and a singular parent fibre. The new finite-field arithmetic has an
independent implementation. The atlas' original section construction and
complete lattice enumeration remain inherited evidence. The k=2,j=1
conclusion also inherits the explicitly recorded norm-eight singular-pencil
replay gap; the new torsion test does not repair that gap. Written geometry
and descent are not proof-assistant verification.

The present branch-field theorem left `k=0,j=0,...,4` and `k=1,j=0,...,3`.
The later [reciprocity proof](Q80_SINGLE_BRANCH_RECIPROCITY_2026-09-14.md)
closes k=1 using the additional inherited-section squareclasses. The five
k=0 cases require degree-two rational effective divisors on the genus31
splitting curve. Their absence, the remaining full coefficient-system
exclusion and the positive MW17 construction remain unknown. Enlarging
either completed finite bank will not resolve this gate.
