# The order-four mechanism requires an isotropic transcendental lattice

**Written structural obstruction, with exact local lattice checks.** Let X be
a projective K3 over a subfield of R, with geometric Picard rank19 and every
geometric NS class Galois invariant. If its actual transcendental lattice T
is anisotropic over Q, every real symplectic involution of X has a real fixed
point. Consequently no pointed elliptic involution can act on the rational
base without real fixed points.

In particular, the [order-four conic mechanism](ORDER_FOUR_LIFT_CORRELATED_GAIN_MECHANISM_2026-09-15.md)
for two correlated gains requires **T isotropic over Q** under these marking
and Picard-rank hypotheses. This is necessary, not sufficient. The
[two-gain objective](CORRELATED_QUADRATIC_GAINS_2026-09-12.md) remains open;
no obstruction to nonsymmetry constructions is claimed.

## 1. The invariant real lattice

Put L=H²(X(C),Z), let c denote complex conjugation, and set K=L^c. As in the
[real marking argument](REAL_MARKING_CONTENT_OBSTRUCTION_2026-09-14.md), c acts
as -1 on invariant divisor classes in untwisted Betti cohomology. Therefore
K lies in T. On the real period plane, c has one eigenvalue+1 and one-1.
Since rank T=3, rank K is either1 or2.

The lattice K is primitive, even and 2-elementary. Here is the elementary
integral argument for the last assertion. Every integral functional on K
extends to L because K is primitive; unimodularity of L represents the
extension by some x in L. Its orthogonal projection to K tensor Q is
(x+c(x))/2. Thus every element of K* has twice its value in K, so2A_K=0.
Its discriminant group needs at most rank K generators.

If rank K=2, it has signature(1,1). Indeed an ample class is in the negative
c-eigenspace, and the positive period plane supplies one positive direction
to each eigenspace; L has only three positive directions. Hence K has negative
determinant, with absolute value1,2 or4. An even binary Gram

```
[2a b]
[b 2d]
```

has determinant4ad-b² congruent to0 or3 modulo4. This excludes determinant-2.
The remaining determinants-1 and-4 make the associated rational binary
quadratic form isotropic, since its discriminant b²-4ad is a square. This
would give a nonzero rational isotropic vector in T. Under the anisotropy
hypothesis, rank K must therefore be1.

For background on real K3 invariant lattices and their topology, see
[Heckman–Rieken, section5.3 and Theorem5.11](https://link.springer.com/article/10.1007/s00208-017-1587-2).
The rank-two determinant argument above does not require a lattice census.

## 2. Real topology forces a fixed point

The topological Lefschetz formula for c gives

```
chi(X(R)) = 2 + trace(c|H²) = 2 + (1-21) = -18.
```

The real locus is a compact orientable surface, possibly disconnected.
Choose a holomorphic2-form omega with c*omega=conjugate(omega). Its restriction
to the real locus is a nowhere-zero real2-form and supplies the orientation.
A real symplectic involution eta preserves omega and therefore this orientation.

If eta had no real fixed point, X(R) -> X(R)/eta would be an unramified double
cover of compact orientable surfaces. Its quotient would have Euler
characteristic-9. Every compact orientable surface, and every disjoint union
of them, has even Euler characteristic. This contradiction proves the fixed
point assertion. This proof uses only the Euler characteristic; no particular
connected-component classification is needed.

A pointed nonsymplectic elliptic involution can be composed with fibre
inversion to give a symplectic involution with the same action on the base.
The two commute because the involution preserves the origin. If that base
action had no real fixed point, neither involution could fix any real point
of X. The preceding argument excludes this as well.

## 3. Consequence for correlated gain by order four

The retained rational-conic construction requires the parent's base involution
to have fixed-point field Q(i), hence no real fixed point. Section2 rules this
out for every elliptic fibration and every such involution on X when T is
anisotropic. Changing frames, the rational zero section or the equation
cannot avoid this obstruction. The genus-one order-four construction over Q
is separately impossible by its action on regular differentials.

This restricts the parent search for this mechanism to isotropic T, unless
the full rational NS or Picard-rank19 hypothesis changes. It does not claim
that an isotropic T admits an appropriate rational marking, rootless frame,
fixed-point field or new section.

## 4. Retained parent applications

The actual integral lattices are

```
948:  [[-2,0,1],[0,4,0],[1,0,118]],
1092: [[-2,1,0],[1,2,2],[0,2,220]],
1020: [[-2,1,0],[1,2,0],[0,0,204]].
```

The [checker](scripts/verify_real_c4_anisotropy_gate.py) verifies each literal
matrix against its retained source, checks its determinant, and tests all702
vectors modulo9 that are primitive at3. None has norm0 modulo9. A rational
isotropic vector could be scaled to an integral vector primitive at3 and
would give such a residue, so all three lattices are anisotropic over Q.
The [certificate](../artifacts/generated-results/elkies-k3-real-c4-anisotropy-gate-v1/result.json)
retains the matrices, source hashes, counts and elementary determinant/parity
checks. The full rational marking and Picard rank are inherited from the
actual parent theorems; the topology is written mathematics, not machine proof.

For determinant1020 this excludes the mechanism on **all** its elliptic
fibrations, strengthening the earlier one-frame computation without a frame
census. For948 and1092 it explains the arithmetic failure independently of
the existing complete frame gates. Those gates remain stronger for other
symmetry questions because they also exclude the specified pointed involutions
with real fixed base points.

The isotropic Inose-type lattice U+<622> is retained as a boundary control:
(1,0,0) is a nonzero null vector. This theorem does not exclude it. Its separate
root obstruction and fixed-point-field gate retain their own scopes.

```
python3 research/elkies-k3/scripts/verify_real_c4_anisotropy_gate.py
```

This is a new deduction in the repository, not a claim of literature novelty,
formal verification or external review. No new two-gain cover is constructed.
