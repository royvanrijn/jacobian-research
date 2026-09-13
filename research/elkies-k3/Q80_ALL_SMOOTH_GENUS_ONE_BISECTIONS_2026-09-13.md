# Smooth genus-one bisection injectivity on alternate Q80

On the [direct11952 alternate-Q80 MW17 parent](R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md),
**two independent new directions on one quadratic cover cannot both have
smooth genus-one bisection images**, even after inherited section translation.
This now covers every such image on this parent, rather than only those
disjoint from the zero section.

The new calculation constructs the full **49 minimum-norm-twelve pencils**
and excludes all 1,176 pairs between them and all 3,131,933 pairs with the
[63,917 minimum-norm-eight pencils](Q80_COMPLETE_GENUS_ONE_BRANCH_INJECTIVITY_2026-09-13.md).
Every individual pencil is injective. Together the two certificates cover
**63,966 pencils and 2,045,792,595 unordered pairs**.

This is an exclusion for smooth images of arithmetic genus one. A genus-one
cover can map birationally to a singular curve of larger arithmetic genus on
the K3; such images remain outside the theorem. Other parents and rational
normalizations are also outside this new assertion. The
[two-gain construction goal](CORRELATED_QUADRATIC_GAINS_2026-09-12.md) remains open.

## The missing full pencil

Write the parent as `y^2=x^3+A(t)x+B(t)`. For a retained height-twelve trace
`tau=(Nx/h^2,Ny/h^3)`, choose a finite chart in which

```
h is monic, degree(h)=4, gcd(h,Nx)=1,
degree(Nx)<=12, degree(Ny)<=18,
degree(M0)<8, M0*Nx+Ny=0 mod h^2.
```

Set `c=coefficient(M0,t^7)`. It is nonzero for all 49 classes, as certified
by their good reductions. For `[u:v]` in the complete parameter line put

```
g=v*t-u,                 M=g*M0-v*c*h^2,
Q=(M^4-6*M^2*Nx*g^2-8*M*Ny*g^3-3*Nx^2*g^4-4*A*h^4*g^4)/h^6.
```

Then **Q is polynomial of bidegree at most (4,4)** in the old base and the
new homogeneous parameter. The congruence and the Weierstrass identity imply
divisibility by `h^6`, by the same local expansion as the
[regular chord theorem](RANK_MUTATION_AND_LIFT_THEOREMS.md#proposition-f12-the-height-twelve-regular-quartic).
Equivalently, over the parameter function field apply that theorem to
`M/g`, which has the same residue modulo `h^2` as `M0`, and then clear `g^4`.
The monic normalization makes the degree-eight terms of `M` cancel, so
`degree_t(M)<=7`. The numerator has degree at most28; division by the monic
degree24 polynomial `h^6` gives the stated degree bound. Homogeneity gives
parameter degree four.

For `v=1`, the chord of slope `m=M/(g*h)` through `-tau` has residual
discriminant `h^2*Q/g^4`. On `W^2=Q(t,u,1)` its points are given by

```
x=(M^2-Nx*g^2+h^3*W)/(2*g^2*h^2),
y=m*(x-Nx/h^2)-Ny/h^3.
```

Changing the sign of W gives the conjugate point; the two points sum to tau.
These are identities on the generic curve, with the usual smooth extension
at apparent poles. At the moving contact `t=u`,

```
Q(u,u,1)=c^4*h(u)^2.
```

Thus the two rational points are `W=+/-c^2*h(u)`. They come from the old zero
and tau. The parameter `[1:0]` recovers the unique regular quartic from the
old norm-twelve construction. A fixed family `M0+lambda*h^2` misses the moving
finite pole and does not substitute for this full pencil.

## Divisor and completeness

Use `NS=U+M(-1)` with coordinates `(a,d,w)` and pairing
`(a,d,w)^2=2*a*d-w^2`. The old fibre is `(1,0,0)` and the old zero is
`(-1,1,0)`. For a minimum-norm-twelve word w,

```
D=(3,2,w)=O+tau-F,        D^2=0, D.F=2, D.O=1.
```

This class is primitive. It is nef and has no orthogonal roots, as follows
directly from the minimum-norm condition. An effective root C has
`d=C.F>0`, since the old frame is rootless. Write `C=(a,d,z)` with
`2*a*d-z^2=-2`. Then

```
D.C=(|z-d*w/2|^2-2)/d.
```

For odd d the numerator is at least `12/4-2=1`, using the minimum of
`w+2M`. For even d, `z-d*w/2` lies in M. If nonzero its norm is at least4,
so again the numerator is positive. If it is zero, the root equation would
give `2*a*d=3*d^2-2`, impossible for even d by reduction modulo4. Hence every
effective root meets D positively. Riemann--Roch makes D effective because
`D.F>0` rules out `-D`; a negative intersection with an irreducible component
of an effective D would require an effective root. This proves nefness.
The primitive nef isotropic pencil lemma on a K3 now gives a base-point-free
genus-one pencil `|D|` with `h^0(D)=2`. Its fibres have no reducible components.

For completeness, the explicit functions producing this linear system can
be checked without a new neighbour enumeration. Put

```
ell=(y+Ny/h^3)/(x-Nx/h^2),
f0=(-M0+h*ell)/h^2,       f1=t*f0+c.
```

The congruence cancels vertical poles at the zeros of h. The only horizontal
poles are at most simple along O and tau. The degree bounds make both
functions vanish along the old infinity fibre. Thus `f0,f1` belong to
`H^0(O+tau-F)`. They are independent because `c!=0` and `ell` is not a
function of t. They form a basis, and `u=f1/f0` gives precisely the moving
slope above. This is also the linear Riemann--Roch construction used by the
[retained direct two-neighbour compiler](scripts/compile_r17_norm12_orbit11952_qq.sage).

Every smooth genus-one bisection has square zero. In these coordinates its
class is `(n,2,w)` with `w^2=4n`. Translation by an inherited section changes
w by twice a lattice vector and preserves smoothness. Minimize the norm in
its parity class. The
[complete retained parity certificate](R17_PRODUCT_TATE_COHOMOLOGY_REDUCTION_2026-09-04.md)
has nonzero isotropic minima `4:1313, 8:63917, 12:49`, and the zero class has
minimum zero. Minima zero and four would give intersection `n-2<0` with O,
which is impossible for an irreducible bisection distinct from O. Therefore
every smooth genus-one bisection translates into one of the 63,917 norm-eight
or 49 norm-twelve pencils used here. The complete parity enumeration is an
inherited theorem, not recomputed in this calculation.

## Rational points on the covering family

In the norm-twelve fibration, O is a zero section and tau is another section:
both have intersection one with D. Their mutual intersection is the old
`tau.O=4`. Since the new frame is rootless, the Shioda height of tau in the
new fibration is `4+2*4=12`, so it is nontorsion. In particular these covering
curves have a generic nontorsion rational point, and all but finitely many
rational parameters with smooth fibre have positive rank by specialization.

The obstruction found here is the lack of a second shared branch image;
the norm-twelve covering family itself has a source of infinitely many
rational points. No arithmetic conclusion is inferred for a hypothetical
isolated branch collision without checking its fibre.

## Exact projective comparison and replay

The [frozen packet](../artifacts/generated-results/elkies-k3-q80-norm12-moving-pencils-v1/input.json.gz)
contains only the generic model, basis, height Gram and all49 retained trace
words. It fixes primes `521,523,541`, a limit of60 CPU seconds per stage,
4 GiB of address space, and blocks of512 norm-eight images. No exceptional
points, old product targets or specialized ranks enter the construction.
The first two primes finish the comparison, so541 is unused.

At each used prime the old discriminant remains squarefree of degree24 and
the basis height Gram is unchanged. The recorded trace frames have four
finite poles, with no cancellation. Thus the characteristic-zero pole
polynomial has degree four and unit leading coefficient. Joint primitive
reduction of numerator and denominator, the unit resultant of `h,Nx`, and
division by monic powers of h make the entire moving branch matrix integral
at the prime. Its reduction is the recorded matrix. All98 new frames use
the original finite base chart.

The replay checks that no rational projective parameter over either finite
field gives a zero branch vector. Any rational equality of branch images
therefore reduces to a recorded collision, including parameters at infinity
and rational parameters of unbounded height. Equal quadratic fields from
smooth genus-one bisections require proportional binary branch quartics in
the fixed original base coordinate. Proportionality is used only as a
necessary condition; no scalar squareclass is discarded in a positive claim.

| Prime | Norm-twelve pencils | Norm-eight pencils needed | Projective parameters | Cross pairs left | Norm-twelve pairs left |
|---:|---:|---:|---:|---:|---:|
|521|49|63,917|33,390,252|40|0|
|523|49|40|46,636|0|0|

The [constructor](scripts/compare_q80_norm12_moving_pencils.sage) used3.557
and0.734 CPU seconds. The
[independent replay](scripts/verify_q80_norm12_moving_pencils.py) took8.166
seconds. It verifies98 norm-twelve trace frames and one newly required
norm-eight frame; the remaining norm-eight frames are inherited from the
previous theorem and checked against its hash manifest. All **33,436,888**
projective images in this extension are evaluated again.

The new trace Weierstrass identities have degree at most36 and are checked
at37 distinct field elements. The moving-branch identities have bidegree at
most `(28,4)` before division and are checked on a complete29-by-5 grid.
These are exact interpolation certificates with checked degree bounds.
The rational-point diagonal is also checked at9 values, matching its
degree-eight bound.

To identify each reported trace with its prescribed word, the checker uses
27 smooth fibres and direct finite elliptic-curve group arithmetic. At most
four are trace poles, leaving at least23 affine comparisons. Two distinct
height-twelve sections agreeing there would have a difference of height at
least `4+2*23=50`, whereas Cauchy--Schwarz bounds it by48. Thus identification
is exact. The branch-image replay uses Horner evaluation and normalization
by the last nonzero coordinate; the producer uses matrix multiplication
and the first nonzero coordinate.

At521,47 new branch matrices are invertible; at523, all49 are invertible.
Consequently every characteristic-zero parameter map is a rational normal
quartic followed by an invertible linear map, and is injective. Together
with complete cross-pair exclusions and the inherited norm-eight theorem,
this proves the stated injectivity result. Two lifts with the same bisection
image differ only by deck conjugation; after undoing inherited translations,
they represent at most one new direction modulo the inherited subgroup.

```
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 research/elkies-k3/scripts/verify_q80_norm12_moving_pencils.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -m unittest discover -s research/tests -p 'test_q80_norm12_moving_pencils.py' -q
```

The [replay receipt](../artifacts/generated-results/elkies-k3-q80-norm12-moving-pencils-v1/independent-replay.json)
and [execution receipt](../artifacts/generated-results/elkies-k3-q80-norm12-moving-pencils-v1/execution.json)
retain the precise checks, timings and one startup failure before any frame
was constructed. Seven regressions cover trace relabelling, interpolation
bounds, monic cancellation, replacing a pencil by a fixed member, the
23-fibre identification bound, shared buckets, and a missing final stage.

## Scope correction for earlier regular-slope exclusions

The moving term also corrects the carrier scope of the earlier
[rational-V4 deep-trace screen](R17_RATIONAL_V4_DEEP_TRACE_EXHAUSTION_2026-09-04.md).
Its formula `M0+lambda*h^2` allows only a constant lambda. The full norm-twelve
genus-one pencil instead has `lambda=-c/(t-u)`. For generic u its lift meets
zero once and has height `2*4+2=10` on the quadratic pullback. The trace has
pullback height24, so the anti-invariant difference `2R-tau` has height
`4*10-24=16`, or twist height8. These are precisely moving-contact carriers
omitted by the assertion that every height-ten half-point has the earlier
fixed regular slope.

The historical screen remains a complete comparison against its 4,358,409
targets for the 49 regular-slope families actually tested. It is not an
exhaustion of all norm-twelve integral coboundaries. The canonical claim and
its foundry routing restriction are narrowed accordingly, with the previous
source and certificates retained. No V4 comparison is rerun. Separate
arithmetic rank-zero proofs for the seventeen selected product twists remain
valid and those targets are not reopened. The parity census used in this
note is independent of the faulty carrier-exhaustion implication.

## Remaining construction gate

A solution on this parent with a genus-one quadratic base must use a singular
bisection image of arithmetic genus at least two in at least one independent
direction. This does not exclude such a solution or bound arbitrary twist
rank. Repeating either complete smooth pencil layer cannot produce it.
The next gate must concern singular images with certified normalization,
or a distinct parent with a new arithmetic mechanism; it must still construct
the actual squareclass and independent sections and certify infinitely many
rational points on the base.
