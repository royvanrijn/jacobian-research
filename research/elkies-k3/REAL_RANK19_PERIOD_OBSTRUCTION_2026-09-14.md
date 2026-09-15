# Rank-19 real marking: an isotropic vector or a vector of square two

Let X be a projective K3 over a subfield of R, with geometric Picard rank
19 and its full geometric NS Galois invariant. Its actual integral
transcendental lattice T must satisfy at least one of:

1. T contains a nonzero rational isotropic vector;
2. T contains an integral vector of square 2.

This necessary condition strengthens the
[content obstruction](REAL_MARKING_CONTENT_OBSTRUCTION_2026-09-14.md).
It is not a sufficient condition for a rational marking, a non-CM period,
or an MW17 fibration.

## Proof

Let L=H²(X(C),Z), and s=-c where c is Betti pullback by complex
conjugation. As in the content proof, s is an integral isometric involution
fixing NS pointwise. On the positive real Hodge period plane it has one
positive and one negative eigenvalue. Write M=L intersect ker(s+I).
Then M is primitive in L and lies in T. Since T has signature (2,1), M
has either rank one and signature (1,0), or rank two and signature (1,1).

Both eigensublattices of an involution of a unimodular integral lattice
have discriminant groups killed by two. Here is the needed argument.
Given y in M*, primitivity lets its integer-valued functional on M extend
to L, and unimodularity represents that extension by x in L. The
orthogonal projection of x onto M tensor Q is y, so 2y=x-s(x) belongs to
M. Thus 2M* is contained in M. This standard integral-involution fact
is also used in the real-K3 discussion in
[Heckman–Rieken, Hyperbolic geometry and moduli of real curves of genus three](https://doi.org/10.1007/s00208-017-1587-2).
The proof here does not require a classification of real K3 surfaces.

If M has rank one, it is positive, even, and has discriminant group of
exponent at most two. Its generator therefore has square 2.

If M has rank two, its negative determinant has absolute value 1, 2, or
4. An even binary Gram matrix has determinant 4ac-b², hence determinant
congruent to 0 or 3 modulo four. This rules out -2, leaving -1 or -4.
In either case minus the determinant is a rational square, so its binary
quadratic space of signature (1,1) is rationally isotropic. Since M is
contained in T, so is T. This proves the dichotomy.

The full geometric NS and actual primitive T are essential. A partially
rational marking, an unsaturated subgroup, a Picard rank other than 19,
or a field without a real embedding is outside this theorem's scope.

## Cheap exact corollary and catalogue application

If T is rationally anisotropic and all its integral squares are divisible
by four, neither alternative holds. The latter condition is checked on a
Gram matrix G by requiring every diagonal entry to be divisible by four
and every off-diagonal entry to be even. No vector search is needed.

The [certificate](../artifacts/generated-results/elkies-k3-real-rank19-period-v1.json)
applies this corollary to the retained 827-row catalogue after the content
gate. It finds 31 additional exclusions, including determinants 736 and
480. For each, rational diagonalization is recomputed from the literal
Gram matrix, and the resulting ternary Clifford Hilbert symbol is -1 at
the retained obstruction prime. The latter is an exact local anisotropy
certificate, not a bounded rational-point miss.

The combined filters leave 722 unresolved rows from the historical
820-row queue. Determinant 800 is the only remaining member of the
retained 21-row coarse-low-genus diagnostic shortlist. Its passing result
is not a full marked-curve or point certificate. No equation work is
authorized by these filters alone.

## Reproduction

```sh
sage -python research/elkies-k3/scripts/audit_real_rank19_period.py --check
```

The audit is limited to the retained 827 inputs, with 31 exact local
recomputations. The Hodge and integral-involution proof is written, not
machine formalized. Independent implementation and external review are
unclaimed. The earlier content packet and historical planner remain
unchanged.
