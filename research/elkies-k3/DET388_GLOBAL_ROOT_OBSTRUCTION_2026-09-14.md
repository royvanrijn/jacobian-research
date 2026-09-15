# Determinant 388: the admitted NS has no MW17 fibration

Every positive even rank-17 lattice with the discriminant form of the
admitted determinant388 frame contains a root. Hence **the actual NS of
K3-f5168aef816b1b25 admits no rootless elliptic frame and no MW17 fibration**.

The [full rational rank19 marking](DET388_RATIONAL_MARKING_SOURCE_2026-09-14.md)
remains valid. This is a geometric obstruction after successful arithmetic
admission, not a retraction of that marking. It supersedes the UNKNOWN
rootless endpoint and the need to enlarge the bounded neighbour search.
It is not a determinant-only theorem about every possible genus at388.

## A smaller auxiliary certificate

Let F be the positive source frame and use the positive even auxiliary

```text
 2  1  1  1 -1  1  0
 1  2  0  1  0  0  0
 1  0  2  1  0  0  0
 1  1  1  2  0  0  0
-1  0  0  0  2 -1  0
 1  0  0  0 -1  8 -1
 0  0  0  0  0 -1 14
```

Call it K. Its determinant is388. The first five vectors generate D5;
an exact unimodular change to its Cartan Gram is in the certificate.
The remaining two vectors, after orthogonal projection away from D5,
have Gram

```text
 7 -1
-1 14
```

Their total squared norm is therefore21. Exact cyclic discriminant-form
arithmetic gives q_K opposite to q_F. The checker constructs the index388
primitive graph glue of K+F and checks that its rank24 Gram is positive,
even and unimodular.

For **any** frame in the same discriminant-form genus, the same gluing
argument embeds K primitively in a Niemeier lattice with that frame as its
orthogonal complement. The Leech lattice is impossible because K contains
D5 roots. Thus a global frame exclusion reduces to all D5 embeddings into
the rooted Niemeier lattices.

## Two-vector root-avoidance bound

For an irreducible ADE root system R with Cartan matrix C, suppose two
vectors u,v pair integrally with every root, and no root is orthogonal to
both. Define, for each subset S of its simple nodes,

    b_R(S) = 1_S^t C^-1 1_S + 1_Z^t C_Z^-1 1_Z,
    Z = complement of S.

The second term is zero for empty Z. Then

    ||u||² + ||v||² >= min_S b_R(S).

Norms here may be taken after projection into the root span, so the bound
also applies to vectors in a larger Euclidean space.

To prove it, move u into a closed Weyl chamber. Its integral Dynkin labels
are positive integers on a support S and zero elsewhere. Positivity of
entries of inverse finite Cartan matrices gives
||u||² at least 1_S^t C^-1 1_S. The roots orthogonal to u form the parabolic
root system on Z. Since none may also be orthogonal to v, its projection
is regular there. Move that projection into a chamber; its integral simple
labels are all at least one, giving squared norm at least
1_Z^t C_Z^-1 1_Z. Adding proves the bound.

Bounds add across orthogonal root systems. For large components the checker
uses a weaker bound from an explicit orthogonal subsystem: A_n contains
A_a+A_b with a+b=n-1, and D_n contains D_a+D_b with a+b=n. It computes the
remaining support minima by exhausting at most512 subsets per component.
No assumption about optimality of these weaker bounds is needed.

## All Niemeier cases

The retained [complete D5 orbit calculation](ROOTLESS_J2_COMPLETENESS_TRACK_2026-08-31.md)
proves that13 rooted Niemeier classes have16 D5 anchor orbits. In E6, E7
and E8 there is one internal Weyl orbit each; all other rooted classes and
the Leech lattice are excluded from containing D5. This existing
classification is a named theorem input, not inferred from the new search.
The checker re-enumerates the roots orthogonal to each retained anchor and
verifies the following exact lower bounds for two integral-pairing vectors
to avoid them all.

| Niemeier root system | Anchor | Required norm sum, at least |
|---|---:|---:|
| D24 | 1 | 69/2 |
| D16+E8 | 1 | 51 |
| D16+E8 | 2 | 38 |
| 3E8 | 1 | 72 |
| 2D12 | 1 | 34 |
| A17+E7 | 1 | 49/2 |
| D10+2E7 | 1 | 32 |
| D10+2E7 | 2 | 45 |
| A15+D9 | 1 | 89/4 |
| 3D8 | 1 | 38 |
| A11+D7+E6 | 1 | 74/3 |
| A11+D7+E6 | 2 | 133/6 |
| 4E6 | 1 | 63/2 |
| 2A9+D6 | 1 | 30 |
| 4D6 | 1 | 30 |
| 2A7+2D5 | 1 | 99/4 |

The smallest bound is133/6, strictly greater than21. Project the last two
K generators away from D5. They still pair integrally with its orthogonal
roots, and their norm sum is21. Hence in every case a root is orthogonal to
both, and thus to all of K. It belongs to F. This contradicts rootlessness
and proves the global exclusion.

## Discovery, replay and exact boundary

A bounded rank-seven auxiliary search found63 classes whose weighted
isometry mass sum equals3169/5760; class61 supplied the displayed K.
The search was bounded by20000 neighbours,256 classes and120 seconds;
it reached the mass after4689 neighbours in about seven seconds. The
[discovery packet](../artifacts/generated-results/elkies-k3-det388-global-root-gate-v1/aux-classification.json)
and failed candidate auxiliaries are retained. Neither that classification
nor its completeness is needed for the new theorem: the explicit K,
discriminant anti-isometry and all16 root-energy inequalities suffice.

```sh
sage -python research/elkies-k3/scripts/certify_det388_global_root_obstruction.py --check
```

The [certificate](../artifacts/generated-results/elkies-k3-det388-global-root-gate-v1/certificate.json)
contains primitive glue, the projected Gram, every support-bound value,
all16 residual root decompositions and pinned input hashes. Written gluing,
Niemeier classification, the existing D5 orbit theorem and the energy lemma
are explicit proof inputs. Independent implementation, formal verification,
external review and literature novelty are unclaimed.

This closes the determinant388 MW17 branch. Other NS genera, including any
other genus of the same determinant, and the overall search for a different
arithmetic MW17 K3 remain open. No maximum MW rank below17 is asserted.
