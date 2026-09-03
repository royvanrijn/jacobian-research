# Small-genus defect-graph dynamics — 2026-09-03

## Result

Complete good-prime graphs in three mass-closed positive even ternary genera
show that unrestricted neighbour connectivity and defect-directed
reachability are genuinely different.

For the calibration used here, the physical defect of a lattice `K` is its
complete signed norm-two set `Phi(K)`.  A state has zero support precisely
when it is rootless.  A good-`p` line is directed when it is nonorthogonal
modulo `p` to every current root.  The survival theorem H0i then says that
every old root dies; the nonzero affine layers of the birth--death theorem
H0i.1 may still create replacement roots.

The exact census gives:

| genus | classes | zero states | prime set | unrestricted graph | directed distances from defective states | closed defective SCCs |
| --- | ---: | ---: | --- | --- | --- | --- |
| rank 3, determinant 112 | 4 | 1 | `{3}` | strongly connected | `infinity, 1, infinity` | one two-state cycle |
| rank 3, determinant 112 | 4 | 1 | `{5}` | strongly connected | `1, 1, 1` | none |
| rank 3, determinant 126 | 3 | 1 | `{5}` | strongly connected | `1, 2` | none |
| rank 3, determinant 316 | 9 | 6 | `{3}` | strongly connected | `1, infinity, infinity` | two singletons |
| rank 3, determinant 316 | 9 | 6 | `{5}` | strongly connected | `1, 1, 1` | none |

Every row is a complete quotient graph on integral isometry classes, not a
sample.  The determinant-112, -126, and -316 reciprocal-automorphism sums are
respectively

```text
3/4, 3/4, 39/16,
```

and equal the exact Minkowski--Siegel masses of their genera.

The generated certificate is
[`../artifacts/generated-results/elkies-k3-small-genus-defect-graphs-v1.json`](../artifacts/generated-results/elkies-k3-small-genus-defect-graphs-v1.json).
It retains every state Gram matrix, automorphism order, theta shells through
norm twenty, full signed physical root coordinates and Gram matrix, root
complement, every projective isotropic line, its child class, and its exact
death/birth counts.

## What the traps say

In determinant 112, all three defective states have exactly two signed roots.
At `p=3`, two of them form the closed directed cycle

```text
D0 <-> D2,
```

while `D1` has a directed edge to the unique zero state.  The full
3-neighbour graph is nevertheless strongly connected.  Thus mass closure or
ordinary Kneser connectivity does not imply connectivity after imposing the
current-witness-killing rule.

The determinant-316 control is sharper.  At `p=3`, one defective class is a
directed self-trap and a second, with four signed roots, has no directed
isotropic line at all.  Both are escaped in one step at `p=5`.  These are
therefore exact **fixed-prime traps**, not traps for the graph containing all
good primes.

This answers one interpretation of the trap question positively: a closed
defective region can exist after the prime set is fixed, even when the
underlying good-`p` graph is connected and a zero state exists in the genus.
It does not exhibit a region closed under every good prime.

## Distance and the failure of defect count

The determinant-126 directed 5-graph contains the exact shortest path

```text
D1 -> D0 -> Z0
```

with signed defect counts

```text
2 -> 2 -> 0.
```

The first edge kills both old roots and replaces them with two new roots.
Thus equal defect count can occur at distances two and one, even in a
three-state mass-complete genus.  This is the smallest control here showing
that defect cardinality is not a distance function.

The useful local quantity is instead the **directed transition profile**:
for each good prime, retain the multiset of child physical signatures over
all current-witness-killing isotropic lines.  In determinant 126 the two
equal-defect states have respectively four and zero directed lines landing
immediately at zero.  This vector-valued profile detects the distance
asymmetry that the defect count misses.

Other exact separators occur, but should not be promoted prematurely.  In
both 3-primary trap controls the automorphism-group order separates the
reachable and unreachable defective states:

```text
determinant 112: reachable {4}, unreachable {8};
determinant 316: reachable {4}, unreachable {8,16}.
```

The theta shells and root-complement classes also distinguish the states.
None of these data is monotone in the examples, and no universal obstruction
theorem is claimed.

## Literature consequence for the all-prime question

Chenevier's [`p`-neighbour statistics](https://arxiv.org/abs/2104.06846)
separate two facts that matter here.  In a genus of rank greater than two
which is one spinor genus, every good-`p` Kneser graph is connected; moreover,
for a fixed source and target class, sufficiently large compatible primes
have a neighbour in the target class.  The paper also formulates the
asymptotic with level structures, which is the relevant setting when the
defect mask depends on a discriminant marking rather than only on the
unmarked lattice.

Consequently, if a zero-defect state lies in the same spinor/level component
and the physical defect is preserved by that marking, the large-prime theorem
predicts something stronger than a descending path: for all sufficiently
large compatible primes there is a one-click neighbour isometric to the
zero-defect target.  Such an edge necessarily kills every current physical
witness by H0i.  Since the state set is finite, choosing one such prime for
each source also gives a finite, but non-effective, prime set sufficient for
directed reachability.

The remaining global obstruction is therefore not ordinary class-set
connectivity inside a one-spinor-genus component.  It is the marked
spinor/level decomposition, together with making the large-prime statement
effective enough to produce a usable prime set.  This qualification is
essential for the rank-15 reverse masks: an unmarked isometry class need not
preserve the distinguished discriminant summand and graph multiplier.

## Reproduce

Generate the certificate or byte-check it with:

```bash
sage -python elkies-k3/scripts/analyze_small_genus_defect_graphs.sage
sage -python elkies-k3/scripts/analyze_small_genus_defect_graphs.sage --check
```

The run is dependency-free beyond Sage/PARI and takes only a few seconds on
the recorded workstation.  It performs no random sampling.

## Boundary and next experiment

What is proved is an exact finite computation for root defects in three
ternary genera.  It proves fixed-prime traps and nontrivial directed distance.
It does not prove an all-good-prime trap, a universal finite-prime bound, a
scalar Lyapunov function, or the corresponding statement for the rank-15
Q80 completion mask.

The next rank-15 experiment should therefore be more targeted than a large
unfiltered genus walk:

1. attach the Q80 discriminant marking and bridge multiplier to the state;
2. determine its spinor/level component, not only its unmarked genus;
3. enumerate complete neighbour orbits for the smallest feasible good prime;
4. store the full physical-witness Gram/incidence signature and the complete
   directed destination profile;
5. close a component by mass or a level-aware neighbour theorem before
   interpreting an SCC as a trap.

The ternary controls say exactly what must be kept: the prime-labelled
transition profile is already more informative than defect count, while the
underlying unrestricted graph can hide directed traps completely.
