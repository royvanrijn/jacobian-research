# Complete norm-eight genus-one branch injectivity on alternate Q80

The **63,917 complete genus-one pencils** in the minimum-norm-eight,
pole-order-zero layer of the direct11952 alternate-Q80 parent have pairwise
disjoint rational projective branch images. Every individual pencil is also
injective. Thus no two smooth genus-one members in this entire layer provide
the requested two-direction quadratic cover.

All **2,042,659,486 pairs** are excluded, with independent arithmetic replay.
This replaces a small sample by a complete layer on this parent. It does not
exclude other carrier layers, singular normalizations or arbitrary quadratic
twist sections. The [positive construction goal](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
remains open.

## Frozen generic source and exact scope

The parent is the certified
[direct11952 MW17 fibration](R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md).
The retained
[complete norm-eight inventory](R17_ALTERNATE_Q80_PRODUCT_BISECTION_INVERSION_2026-09-03.md#2-complete-norm-eight-layer)
contains exactly 63,917 minimum-norm-eight classes modulo twice its saturated
Mordell--Weil lattice. Its old comparison with seventeen product quartics
does not compare these pencils with each other.

The new [frozen packet](../artifacts/generated-results/elkies-k3-q80-complete-genus-one-pairs-v1/input.json.gz)
projects only the generic equation, generic section coordinates, height Gram,
complete generic trace words and their provenance. The complete table was
already present and was not reconstructed. Its hash is checked before freezing.
No previous product target, exceptional point, specialization parameter or
specialized rank enters the worker.

The input fixes primes `509,521,523,541`, 240 CPU seconds per prime, a 4 GiB
address-space limit, and checkpoints every 512 pencils. Later prime stages
contain exactly the vertices still needed by the preceding pair comparison.
A separately frozen five-class supplement at the already listed prime541
finishes within-pencil injectivity.

The complete class enumeration and the generic parent height lattice are
inherited mathematical dependencies. The new independent replay verifies the
arithmetic comparison, not a fresh Fincke--Pohst enumeration or parent-rank proof.

## Projective branch maps

For a trace `T=(Nx/h^2,Ny/h^3)` of height eight, choose a chart with
`degree(h)=2` and set `M=M0+lambda*h^2`, where

```
degree(M0)<4,          M0*Nx+Ny=0 mod h^2.
```

The [regular chord theorem](RANK_MUTATION_AND_LIFT_THEOREMS.md#proposition-f11-the-height-eight-genus-one-bisection-pencil)
gives five coefficient polynomials of degree at most four:

```
q0=(M0^4-6*M0^2*Nx-8*M0*Ny-3*Nx^2-4*A*h^4)/h^6,
q1=(4*M0^3-12*M0*Nx-8*Ny)/h^4,
q2=6*(M0^2-Nx)/h^2,    q3=4*M0,    q4=h^2.
```

Thus `q_lambda=sum(qj*lambda^j)` has coefficient vector
`B_T*(1,lambda,lambda^2,lambda^3,lambda^4)^t`. Homogenization includes
`lambda=infinity`. For a smooth genus-one member, equal quadratic extensions
over the fixed original parameter require proportional binary branch quartics.
Constant proportionality alone would not prove equality of rational
squareclasses; this calculation only uses it as a necessary condition.

If a trace has a pole at infinity, use `t=c+1/z`, trying `c=0,1,2` after
the ordinary chart. A degree-two pole divisor cannot contain all three values,
so one chart makes both poles finite. Transform the branch polynomial back
by `q_t(t)=(t-c)^4*q_z(1/(t-c))`. All comparisons use the original `t`.
The checkpoints use 124,876 ordinary charts, 211 charts with `c=0`, and one
with `c=1`; these counts precede the five-class embedding supplement.

## Why the reductions give exclusions over Q

At each used prime the short model has a squarefree degree-24 discriminant,
so all fibres are irreducible and the infinity model is minimal. Rational
functions are reduced after scaling numerator and denominator together to
primitive integral polynomials. Reducing their coefficients separately can
incorrectly reject a perfectly valid rational-function reduction.

The checker verifies that the reduced 17-section height Gram equals the
generic Gram. Therefore every selected word still has height eight. The
verified reduced trace has two finite poles in its chosen chart.

These facts justify integral lifting of the branch map. Normalize the
characteristic-zero trace pole polynomial primitively at the prime. Its
degree is at most two by the generic height bound. The reduced trace has
two poles, so its leading coefficient is a unit and no pole cancels on
reduction. Normalize that coefficient to one. Then `h,Nx,Ny` are integral,
and `gcd(Nx,h)=1` on reduction makes the inverse of `Nx mod h^2` integral.
The uniquely defined `M0` and the quotients defining all five `qj` are
integral and reduce to the recorded polynomials. The coordinate changes used
above have unit determinant, so the argument also applies to those charts.

Every rational pencil parameter has a reduction in `P1(Fp)`. The replay
checks that **none** of these projective parameters has a zero branch vector.
Consequently any rational equality of branch images reduces to one of the
computed finite collisions. Disjoint images exclude that pair over Q, without
a bound on parameter height. The matrices need not be invertible at every
comparison prime; nonzero branch vectors at all rational reduction parameters
are the condition used here.

## Independent trace and polynomial certification

The replay uses Python fractions, integer polynomial arithmetic and NumPy
integer arrays. It imports no Sage and does not call the constructor.

First it verifies the generic basis equations over each finite function field
and computes its height Gram independently by rational-function addition.
On a rootless elliptic K3 a nonzero section with reduced abscissa `n/d` has

```
height=max(4+degree(d),degree(n)).
```

For each reported trace, coprimality of `h,Nx` and the degree bounds prove
height eight. Its Weierstrass identity has degree at most24 and is checked
at25 distinct field elements. Each branch-coefficient identity has degree
at most16 and is checked at17 distinct elements. These are exact polynomial
identity tests by interpolation; the degree bounds are checked explicitly.

To identify the trace with its prescribed generic word, choose17 distinct
smooth parameter fibres. At most two are trace poles, leaving at least15
affine comparisons with direct finite elliptic-curve group arithmetic.
If two distinct height-eight sections agreed at15 smooth fibres, their
difference would meet zero at least15 times. Its height would then be at
least `4+2*15=34`. Cauchy--Schwarz instead bounds that height by
`(sqrt(8)+sqrt(8))^2=32`. Thus the reported trace is exactly the prescribed
word over the finite function field. This is not a probabilistic point test.

Finally, the replay computes all projective images by Horner evaluation and
normalizes using the last nonzero coordinate. The constructor uses matrix
multiplication and the first nonzero coordinate. The independently recomputed
collision buckets agree, and bitsets verify complete pair coverage while
retaining every pair incident to an untested vertex.

## Complete result

| Prime | Pencils compared | Projective parameters | Pairs remaining |
|---:|---:|---:|---:|
|521|63,917|33,364,674|118,251|
|523|61,153|32,044,172|9|
|541|18|9,756|0|

The total is **65,418,602 projective parameter images**. Prime509 is unusable:
the discriminant and its derivative have a degree-three gcd there. It supplies
no exclusion. All used primes come from the original frozen panel.

For each of the 63,917 pencils, an invertible branch matrix is also certified
at one good prime. The first stage certifies63,799; the second leaves five
uncertified, and the
[five-class supplement](../artifacts/generated-results/elkies-k3-q80-complete-genus-one-pairs-v1/embedding-result.json.gz)
certifies those at541. Thus each characteristic-zero matrix is invertible
and its rational-normal-quartic parameter map is injective.

The [independent replay](../artifacts/generated-results/elkies-k3-q80-complete-genus-one-pairs-v1/independent-replay-final.json)
verifies125,093 trace frames including the supplement. Eight regressions
cover a relabelled valid trace, interpolation degree bounds, joint rational
normalization, projective scaling and infinity, zero branch vectors,
adaptive pair coverage, a missing final certificate and a bad model prime.

```
OPENBLAS_NUM_THREADS=1 python3 research/elkies-k3/scripts/verify_q80_complete_genus_one_pairs.py
OPENBLAS_NUM_THREADS=1 python3 -m unittest discover -s research/tests -p 'test_q80_complete_genus_one_pairs.py' -q
```

The [execution receipt](../artifacts/generated-results/elkies-k3-q80-complete-genus-one-pairs-v1/execution.json)
retains component timings and both startup failures. The bad-prime rejection
and the corrected coefficientwise rational-function reduction happened before
any pencil was constructed; the frozen input was unchanged.

## Boundary for the construction goal

There is no shared smooth genus-one quadratic cover from two different members
in this complete minimum-norm-eight, zero-intersection layer, up to inherited
section translation. This does not exhaust singular higher-arithmetic-genus
carriers or genus-one carriers with a larger minimum intersection with zero.
It does not prove that an arbitrary twist has rank at most one: not every
integral twist direction necessarily lies in this bisection image. No
specialized rank upper bound follows.

The next geometric gate is to determine which retained higher-intersection
classes supply full genus-one pencils beyond their fixed regular members.
Their divisor, moving intersection with zero and complete branch map must be
proved before another comparison. Repeating this 63,917-pencil layer is
unnecessary.
