# Determinant 622: the Inose NS admits no MW17 fibration

Every positive even rank-17 lattice with the discriminant form of
E8+E8+<622> contains a root. Consequently the NS
U+E8(-1)^2+<-622> of the [level311 Inose source](DET622_INOSE_RATIONAL_MARKING_SOURCE_2026-09-15.md)
admits no MW17 elliptic fibration. Existence of a constant twist with full
rational saturated marking remains proved. Computing that twist or its
primitive section is no longer needed for the different-NS MW17 objective.
This is an exclusion for this NS genus, not for every lattice of determinant622.

## Explicit auxiliary and reduction

Use the positive even rank-seven lattice K with Gram

```text
 2  1 -1  1  1 -1  1
 1  2  0  1  0  0  0
-1  0  2 -1 -1  0 -1
 1  1 -1  2  0  0  1
 1  0 -1  0  2 -1  1
-1  0  0  0 -1 12  1
 1  0 -1  1  1  1 16
```

The first five generators span a primitive D5. The projections u,v of the
last two generators away from D5 have Gram

```text
43/4  7/4
 7/4 59/4
```

Thus their norm sum is51/2. The checker proves that K has determinant622,
and that its cyclic discriminant quadratic form is opposite to that of
F=E8+E8+<622>. It constructs the primitive index622 graph glue of K+F and
checks the resulting positive even unimodular rank24 Gram. The same graph
gluing works for every F with this discriminant form.

Such a glue is a Niemeier lattice. It cannot be Leech, since K contains D5.
The [retained complete D5 orbit theorem](ROOTLESS_J2_COMPLETENESS_TRACK_2026-08-31.md)
reduces to16 anchors in13 rooted Niemeier classes. This completeness theorem
is an explicit input. The checker enumerates the orthogonal roots at every
anchor. If F were rootless, no such root could be orthogonal to both u and v.
Their projections pair integrally with these roots.

## Fifteen strict energy exclusions

Use the [two-vector Weyl-support lemma](DET388_GLOBAL_ROOT_OBSTRUCTION_2026-09-14.md#two-vector-root-avoidance-bound):
for Cartan matrix C the norm sum is at least the minimum, over simple-node
supports S with complement Z, of

    1_S^t C^-1 1_S + 1_Z^t C_Z^-1 1_Z.

For A_n the checker exhausts every support, using
(C^-1)ij=min(i,j)(n+1-max(i,j))/(n+1), with one-based indices.
The complement is a union of A_k zero runs, each contributing
k(k+1)(k+2)/12. This includes all131072 supports for A17 and gives the
stronger bounds A11>=22, A15>=159/4, A17>=905/18. Larger D types use the
same weaker orthogonal-subsystem bounds as the determinant388 proof.

Fifteen of the16 anchor bounds are strictly greater than51/2; their exact
root decompositions and bounds are in the certificate. The sole remaining
case is the D5 anchor in 2A7+2D5, whose orthogonal roots are 2A7+D5.
The coarse support bound there is99/4 and needs refinement.

## Sharp A7 equality excludes the remaining case

Represent A7 in the sum-zero subspace of R^8. Because u and v pair
integrally with its roots e_i-e_j, their coordinates can be written as
integer pairs p_i=(a_i,b_i) centered at their mean. A root orthogonal to
both is precisely a repeated pair. Root avoidance therefore requires eight
distinct points of Z^2. Their norm sum is

    sum_i ||p_i - mean(p)||^2 >= 39/4.

Here is a finite exact proof with its equality information. Integer
translation puts the mean modulo Z^2 at (a/8,b/8), 0<=a,b<8. For each of
these64 possibilities, sum the eight smallest distances to Z^2. The checker
enumerates [-4,4]^2 and verifies that its eighth distance is less than
(5-7/8)^2, a lower bound for every omitted point. The minimum is39/4,
attained only at (a,b)=(1,1),(1,7),(7,1),(7,7). In each case the eighth and
ninth distances are distinct, so the minimizing set is unique. Its centered
coordinate Gram has both diagonal entries39/8. These sets also have the
specified mean. This certifies both the bound and all equality cases.

The support lemma gives6 for D5. Thus 2A7+D5 needs at least
2*(39/4)+6=51/2, exactly the available trace. Rootlessness would force
equality everywhere. In particular the two A7 projections consume39/4 of
each of u's and v's squared norms. The remaining D5 projections therefore
have squared norms at most1 and5, respectively.

In the standard D5 model, its dual vectors have either all integer or all
half-integer coordinates. Norm at most1 leaves just0 and the ten signed
coordinate vectors: the checker exhausts their doubled coordinates in
[-2,2]^5. For any of these vectors, the orthogonal roots contain D4.
Avoiding all those roots forces v to be regular on D4. Its integral Dynkin
labels, after a Weyl move, are all at least one, giving squared norm at
least rho(D4)^2=14. This contradicts the available bound5. Thus the final
anchor also has a root in K-perp=F, completing the global exclusion.

## Reproducibility and boundary

```sh
sage -python research/elkies-k3/scripts/certify_det622_global_root_obstruction.py --check
```

The [certificate](../artifacts/generated-results/elkies-k3-det622-root-gate-v1/certificate.json)
contains the auxiliary, glue, projected Gram,16 anchor checks,64 centroid
bounds, sharp equality Grams, and all small D5 dual vectors. The checker
imports Niemeier classification, the retained complete D5 anchor theorem,
and the written energy and gluing arguments. Independent implementation,
formal verification, external review and literature novelty are unclaimed.

Discovery first constructed a primitive norm622 vector in E8 and its
rank-seven complement, then a bounded p=3 neighbour enumeration found182
classes whose weighted mass equals59683/11520 after14358 neighbours in
about20 seconds. Limits were20000 neighbours,256 classes and120 seconds,
with a checkpoint for each new class. Class158 supplied K. Discovery and
all inputs are retained in the packet; its genus completeness is not
needed for the theorem, which uses only the displayed K.

This closes the determinant622 branch for MW17. It neither retracts the
arithmetic source nor proves a maximum MW rank below17. The overall search
for a different arithmetic MW17 K3 remains open.
