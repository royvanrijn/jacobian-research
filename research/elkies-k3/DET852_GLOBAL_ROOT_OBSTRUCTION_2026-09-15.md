# Determinant 852: the admitted NS has no MW17 fibration

Every positive rank17 frame in the admitted NS genus of
K3-4ff75fec54d01662 contains a root. Hence this NS has no MW17 elliptic
fibration. Its [full saturated rational marking](DET852_RATIONAL_MARKING_SOURCE_2026-09-15.md)
remains proved. This closes the geometric852 branch and supersedes the
remaining-case UNKNOWN, without retracting its arithmetic admission.

## Exhausting the final anchor

The [previous reduction](DET852_ONE_ANCHOR_REDUCTION_2026-09-15.md)
uses a primitive opposite-discriminant rank7 auxiliary K containing D5.
Primitive gluing and all16 Niemeier D5 cases exclude15 cases. The only
remaining ambient has root system2A7+2D5, with residual roots2A7+D5.
The last two generators of K have squared norms14 and18. After projection
away from the anchored D5 their Gram is

```text
51/4  3/4
 3/4 67/4
```

The final checker uses the actual integral Niemeier Gram and retained
anchor. It puts K's D5 block into the Cartan basis by the certified integral
isometry. Isometries of D5 consist of its Weyl group and an order-two
outer diagram symmetry. Ambient root reflections realize the Weyl group.
We explicitly check both remaining orientations, swapping the two terminal
D5 nodes. No extension of that outer symmetry to the ambient glue is assumed.

For each orientation let A be its five anchored simple roots. Enumerate
all ambient roots orthogonal to A, choose simple systems for their three
components, and call their stacked rank19 basis B. The combined root basis
D=(A;B) has rank24. If a vector has prescribed root pairings p, its unique
ambient row coordinates are

    p * (G D^t)^-1.

Integrality of these coordinates is exactly membership in the actual
Niemeier lattice. Root-lattice membership alone is not substituted for
this glue test.

## All sixth vectors modulo the residual Weyl group

Move the projected sixth vector u to the closed dominant chamber of each
residual component. The relevant Weyl reflections preserve the Niemeier
lattice and fix anchored D5 pointwise. Its Dynkin labels are therefore
nonnegative integers. For a component Cartan matrix C, its squared norm is
l^t C^-1 l. The checker exhausts every such label vector with norm at most
51/4 by an exact positive-coefficient recursion.

If Z is the set of zero labels, the roots orthogonal to u form the parabolic
system on Z. Any root-avoiding seventh vector v must be regular there;
its projected norm is at least 1^t C_Z^-1 1. These necessary bounds sum
across components. Requiring total u norm51/4 and total required v norm at
most67/4 leaves52 label combinations. Imposing the fixed D5 pairings and
actual ambient integrality leaves13 sixth vectors, for each orientation.
All52 possibilities are exhausted before the glue test; no heuristic
selection or random neighbour walk enters this step.

## All root-avoiding seventh vectors

For each of these13 vectors, use the Weyl stabilizer of u to make v
strictly dominant on each zero-node subsystem. This fixes both A and u.
Thus v has positive integer labels on the zero nodes and arbitrary integer
labels elsewhere. It has projected norm67/4 and projected pairing3/4
with u.

If the regular lower bounds on the three components are b1,b2,b3, the
individual norm of component j is at most

    67/4 - sum_i b_i + b_j.

The checker enumerates every dual-root-lattice vector up to this bound:
it clears denominators in C^-1 and calls exact positive-definite quadratic
vector enumeration. Stored-vector counts are checked against the full
reported counts, so truncation cannot silently discard candidates. Both
signs are retained; zero is irrelevant because every zero-node subsystem
here is nonempty and strict positivity is required. Filter the vectors by
that positivity condition and record their exact norms and pairings with u.

Across the three components, require norms to sum to67/4 and pairings to
sum to3/4. **There are no matches for any of the26 oriented sixth cases.**
The largest candidate lists in a case have sizes39,39,18; most are much
smaller. The contradiction arises even before imposing the seventh vector's
ambient glue and the rank7 embedding's primitivity. Relaxing these latter
constraints cannot hide a valid embedding.

Therefore the remaining anchor cannot contain K with rootless complement.
Combined with the preceding fifteen exclusions, every frame of this
actual NS contains a root. The NS cannot support MW17.

## Reproducibility and boundary

```sh
sage -python research/elkies-k3/scripts/certify_det852_global_root_obstruction.py --check
```

The [certificate](../artifacts/generated-results/elkies-k3-det852-final-anchor-v1/certificate.json)
retains both D5 orientations, all26 integral sixth candidates and every
allowed component label vector for the seventh, with exact norms and
pairings. The earlier auxiliary, glue and complete15-case reduction are
pinned inputs. The final calculation takes seconds and has finite bounds
derived from the exact projected Gram. Discovery scripts and checkpoints
are preserved in the same directory.

The proof uses the prior primitive-glue and complete-anchor theorems,
standard Weyl-chamber facts and exact quadratic-vector enumeration.
Independent implementation, formal verification, external review and
literature novelty are unclaimed. This is an actual NS-genus exclusion,
not a determinant-wide theorem or an optimal bound below MW17.
The different arithmetic MW17 objective remains open;600 other retained
rows remain unadmitted. Computing852 period coordinates or a new equation
on this NS is no longer a prerequisite for that objective.
