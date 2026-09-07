# Kihara: full first-parent geometry and a restricted parent intake

The first new Kihara parent, at path parameter **3/2**, now has an explicit
full geometric Mordell–Weil basis and Néron–Severi lattice. The minimal
constant field defining them is **Q(sqrt(-3))**. The full geometric NS
rank is eighteen, determinant **-756**, and discriminant group
Z/3 plus Z/6 plus Z/42. Its rational/geometric MW ranks remain **12/13**:
the recovered section is not rational and produces no rational rank gain.

The existing parent intake is restricted in another, directly checkable
way. Writing v=p/q for the Kihara parent coordinate, every nonzero rational
parameter on the retained rank14 path gives **-1/2<v<0**. Positive v is
absent from this intake. This is a coordinate coverage statement; a
different value or sign is not by itself proof of a different unmarked
surface. New parents still require arithmetic and isomorphism checks.

Authority: `EC-KIHARA-FIRST-PARENT-FULL-GEOMETRIC-NS-20260907` and
`EC-KIHARA-UNRESTRICTED-SECTION-AND-PATH-COVERAGE-20260907` in
[`MATH_STATUS.json`](../../MATH_STATUS.json).
The [accounting and source bindings](../../artifacts/generated-results/elliptic-curves/kihara_parent_coverage_followup_v1.json)
retain every completed and failed stage. No new parameter population or
point-search boxes were used here.

## Recovering the missing section

The [preceding exact Picard proof](KIHARA_FIRST_PARENT_PICARD_AND_RANK_2026-09-07.md)
provides the full rational MW basis of rank twelve, determinant189, exact
geometric rank thirteen, trivial geometric torsion, and fibre configuration
20I1+split I4. The quotient by T->-T is a rational elliptic surface with
10I1+split I2 and geometric MW rank seven. Six independent invariant
rational sections descend to it. Its remaining geometric direction has a
nontrivial quadratic character.

The first eight good-prime quotient counts on each of the five retained
parents give forty signs, all consistent with the character of Q(sqrt(-3)).
Those finite signs do not identify a number field. On each quotient the
negative sign, together with the six rational directions, does prove
rational MW rank six. It does not determine the full K3 rank on the other
parents.

[Scholten, section4](https://arxiv.org/abs/math/9709235) supplies the
minimal-section construction for the E7* quotient. Normalize its infinity
node to12 and translate the base to remove the cubic coefficient of A.
Seeking X=12s^2+a*s+b and Y=c*s^2+d*s+e gives five exact equations. Their
saved Groebner basis has quotient dimension56. A conversion to lexicographic
order timed out; that failed run and its saved basis remain intact.

Multiplication by c+a on the saved quotient algebra gives a characteristic
polynomial with32 linear and12 simple quadratic factors. A quadratic factor
with discriminant -3 times a rational square supplies coordinates in
Q(sqrt(-3)). All five original equations and the full section identity pass.
The [portable bundle](../../artifacts/generated-results/elliptic-curves/kihara_quadratic_section_bundle_v1.json)
contains the section and exact transport to the original K3.

## Full geometric basis, without an unproved saturation assumption

Let P be this section, bar its quadratic conjugation, and D=P-bar(P).
Independent direct group-law and height checks give

```
h(P)=3, P.O=0, infinity component(P)=2;
h(D)=4, D.O=0, infinity component(D)=0.
```

Every nonzero geometric section on20I1+I4 has height at least3. Hence D is
primitive: a proper division would have height at most1. Exact geometric
rank thirteen and rational rank twelve imply that the anti-invariant
subspace has dimension one. P exhibits its nontrivial character, so all
geometric sections are defined over Q(sqrt(-3)): any action trivial on their
rational span is trivial on the torsion-free group.

For any geometric section R, R-bar(R)=mD for some integer m. Therefore
R-mP is rational and lies in the already certified full rational basis.
Thus that basis together with P is the **full geometric MW group**.
The exact trace word in the old full basis is

```
P+bar(P) = [3,0,-2,0,0,-1,-1,-1,1,-1,-1,0].
```

The [standalone result](../../artifacts/generated-results/elliptic-curves/kihara_quadratic_section_full_NS_v1.json)
contains all91 pairings, the height determinant189, full rank18 divisor
Gram matrix of determinant-756 and its Galois involution. The latter
preserves the Gram matrix and has fixed rank17. This closes the previously
unknown full geometric lattice and field on this first parent only. Its
all-Q-fibration generic rank bound15 remains unchanged.

## The unrestricted parent line

Set q=1 and v=p/q. Reconstruct the quartic Jacobian from the six centres

```
0, (2v^2+v+2)^2, 2(v+1)^2(2v^2+v+1), 4v^2-v+4,
v(2v-1)(2v^2+4v+5), 4v^4+8v^3+9v^2-2v+2.
```

After normalizing the quotient's leading coefficients, seek
X=-36s^2+k*s+l and Y=sqrt(-3)*(96s^3-6k*s^2+n*s+o).
The coefficient equations determine n,o and leave an elimination cubic in l.
An initial assumption that this ideal would give two linear equations was
false. The saved cubic instead has one linear and one quadratic factor over
Q(v). Extracting its unique rational component gives an explicit section
over **Q(v)(sqrt(-3))(T)**, with X and Y/sqrt(-3) polynomial of degrees4 and6
in T. This is not restricted to the rank14 path.

The [independent verifier](../cas/verify_kihara_universal_anti.sage)
reconstructs the original product, square approximation and Jacobian
directly from those centres. It verifies the full identity using only the
[section bundle](../../artifacts/generated-results/elliptic-curves/kihara_universal_anti_bundle_v1.json),
without the producer's normalized model, equations or Groebner basis.
The [result](../../artifacts/generated-results/elliptic-curves/kihara_universal_anti_replay_v1.json)
records polynomial coordinates with no parameter denominators. Their
Y/sqrt(-3) leading coefficient is nonzero for rational v outside
{0,1/2,-1}. Singular specializations still require a separate check. This
section is nonrational wherever its specialized Y is nonzero; it supplies
no new rational direction. No universal full lattice or exact rank claim
is made.

For the retained path, put z=t^2. Direct polynomial identities give

```
v=-z(8+3z)/(6(z+2)(z+4)),
v+1/2=(5z+12)/(3(z+2)(z+4)).
```

Both signs are strict when z>0. Conversely the inverse equation in z is
(6v+3)z^2+(36v+8)z+48v=0, of discriminant16(9v^2+4).
Rational inverse roots must also be rational squares to come from a
rational t. Enlarging a t-height box therefore does not sample the full
rational v-line. A broader parent intake should operate on that line
directly, deduplicate surfaces, and certify the retained rational sections
before comparing equal point-search exposure. This proof does not assert
that the omitted coordinates contain better specializations.

## Replay and exposure

Run the two portable verifiers with their bundles in an empty directory:

```sh
sage -python verify_kihara_quadratic_section.sage --input kihara_quadratic_section_bundle_v1.json
sage -python verify_kihara_universal_anti.sage --input kihara_universal_anti_bundle_v1.json
python3 elliptic-curves/cas/record_parent_coverage_followup.py kihara --check
```

The first standalone proof takes19.077485987 supervised seconds; the
unrestricted identity replay takes0.689139260. All seven new stages total
**151.735045815 seconds**, including the121.032160787-second timeout and
the1.713220510-second failed linear-extraction assumption. The timeout
includes termination grace. All stages used one worker and2GiB; separate
60- or120-second limits were frozen in their protocols. No failed arithmetic
was silently discarded or rerun with a larger population.

The [subsequent fixed positive-ratio intake](KIHARA_THREE_POSITIVE_PARENT_EXPANSION_2026-09-07.md)
now constructs and independently checks three such parents, distinct over Q
from all previously retained ones, with section spans7,9,12. The [subsequent point pilot](KIHARA_POSITIVE_PARENT_POINT_PILOT_2026-09-07.md)
now records294 completed boxes and three certified directions.
