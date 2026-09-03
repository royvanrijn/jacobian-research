# Small-genus defect-graph dynamics — 2026-09-03

<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-DEFECT-GRAPH-REACHABILITY e02f950eba79b32a -->

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

| genus | classes / zero | analyzed good primes | singleton directed distances | minimum sufficient prime sets |
| --- | ---: | --- | --- | --- |
| rank 3, determinant 112 | `4 / 1` | `{3,5,11}` | `p=3: infinity,1,infinity`; `p=5,11: 1,1,1` | `{5}`, `{11}` |
| rank 3, determinant 126 | `3 / 1` | `{5,11,13}` | `p=5: 1,2`; `p=11,13: 1,1` | `{5}`, `{11}`, `{13}` |
| rank 3, determinant 316 | `9 / 6` | `{3,5,7}` | `p=3: 1,infinity,infinity`; `p=5,7: 1,1,1` | `{5}`, `{7}` |

Every two-prime and three-prime union in the table has universal zero
reachability, maximum distance one, and one directed SCC.  At every singleton,
pair, and triple the certificate labels every SCC, records the condensation
edges, and stores an exact shortest path with a prime and isotropic-line
witness for each step.

Every row is a complete quotient graph on integral isometry classes, not a
sample.  The determinant-112, -126, and -316 reciprocal-automorphism sums are
respectively

```text
3/4, 3/4, 39/16,
```

and equal the exact Minkowski--Siegel masses of their genera.

The generated certificate is
[`../artifacts/generated-results/elkies-k3-small-genus-defect-graphs-v2.json`](../artifacts/generated-results/elkies-k3-small-genus-defect-graphs-v2.json).
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

The unions make this statement completely explicit.  In determinant 112,
both `{3,5}` and `{3,11}` merge the trapped SCC with the zero basin.  In
determinant 316, both `{3,5}` and `{3,7}` do the same.  The corresponding
three-prime unions are strongly connected as directed graphs.  Thus the trap
is destroyed by adjoining one suitable good prime; it is not merely bypassed
while some other defective SCC remains closed.

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

## Minimum-prime-set experiment

For a declared list of good primes, the checker exhausts every nonempty subset
and solves the finite optimization problem

```text
minimize |S| subject to every state reaching a zero state in G_S.
```

The minimum cardinality is one in all three controls, with the complete lists
of minimizers displayed in the first table.  This is still useful negative
routing information: the least analyzed prime can fail even when another
single prime succeeds, so “use the smallest good prime” is not an invariant
strategy.  No present control requires a genuinely mixed two-prime route.

## The all-prime question is a finite-level theorem

<!-- status-consumer: EC-K3-INTEGRAL-RANK-TRANSFER-MARKED-ROOTLESS-REACHABILITY 354cc7a9fc81f33e -->

Chenevier's
[`p`-neighbour statistics](https://doi.org/10.24033/bsmf.2852), Theorem 5.9,
already treat arbitrary compact-open level structure.  The discriminant/glue
markings used by the reverse-mask calculus fit that setting: after putting
the finitely many discriminant, bridge, and marking primes in `S`, the
stabilizer of a distinguished finite quadratic summand, graph subgroup, or
graph anti-isometry is a compact-open subgroup `K^lev` of the finite adelic
orthogonal group.  The marked states are therefore Chenevier's finite class
set

```text
X(K^lev)=O(V) backslash O_V(A_f)/K^lev.
```

The intrinsic marked component is the fibre of the spinor-genus map `s`
after quotienting by `S_1(K^lev)`, the subgroup of good-prime spinor
displacements.  This definition uses the adelic level structure, not graph
reachability.

Theorem H0i.3 now proves the sharp equivalence suggested by the controls:

```text
zero-support reachability under some finite compatible good-prime set
<=> a zero-support state exists in the same marked spinor/level component.
```

Indeed, every good-prime edge multiplies `s` by an element `delta_p` of
`S_1(K^lev)`, which proves the forward implication.  Conversely, for a
rootless target `z` in the same component, put
`a=s(z)*s(x)^(-1)`.  Chenevier's Remark 5.11 gives infinitely many primes in
arithmetic progressions with `delta_p=a`, and Theorem 5.9 gives a direct
`p`-neighbour marked-isomorphic to `z` for every sufficiently large such
prime.  If any parent physical witness survived that line, H0i.1 would place
the same norm-two vector in the zero affine layer of the child completion,
contradicting rootlessness.  The one-click edge is therefore automatically
directed and has no replacement birth.

Since the marked class set is finite, one sufficiently large prime for each
required element of `S_1(K^lev)` gives a finite sufficient prime set for the
whole component.  This is an existence theorem, not an effective algorithm:
the proof gives neither the threshold nor a small-prime bound.

Sage's exact proper-spinor-kernel quotient has order one in each of the three
genera above.  Hence the fixed-`3` traps survive after the ordinary spinor
obstruction has already disappeared.  The experiments and the theorem
together separate two scales of state data:

1. the marked spinor/level component controls eventual all-prime access;
2. at a specified small prime, the prime-labelled incidence of the individual
   physical witnesses with all isotropic lines controls the actual directed
   SCC and distance.

The finite-level hypothesis is essential.  An unmarked isometry need not
preserve the distinguished discriminant summand and graph multiplier, while
an infinite marking that fixes exact rational vectors, a full embedded core,
a nef chamber, or an equation need not have open stabilizer.  Those stronger
markings, effectiveness, and equation-level lifting are outside the theorem.

## Reproduce

Generate the certificate or byte-check it with:

```bash
sage -python elkies-k3/scripts/analyze_small_genus_defect_graphs.sage
sage -python elkies-k3/scripts/analyze_small_genus_defect_graphs.sage --check
```

The run is dependency-free beyond Sage/PARI and takes only a few seconds on
the recorded workstation.  It performs no random sampling.

## Boundary and next experiment

The exact computation proves fixed-prime traps, nontrivial directed distance,
minimum sufficient sets inside the three declared prime lists, and
one-proper-spinor-genus status for all three ternary controls.  Theorem H0i.3
separately proves that no all-good-prime directed trap can persist in a finite
marked component containing a rootless state.  It gives no effective prime
threshold, universal small-prime bound, scalar Lyapunov function, or complete
rank-15 small-prime graph.

The next rank-15 experiment should therefore be more targeted than a large
unfiltered genus walk:

1. attach the Q80 discriminant marking and bridge multiplier to the state;
2. determine its finite-level spinor component explicitly, not only its
   unmarked genus;
3. enumerate complete neighbour orbits for the smallest feasible good prime;
4. store the full physical-witness Gram/incidence signature and the complete
   directed destination profile;
5. use the finite-level reachability theorem to distinguish a genuine
   component obstruction from a merely fixed-prime SCC.

The ternary controls say exactly what must be kept: the prime-labelled
transition profile is already more informative than defect count, while the
underlying unrestricted graph can hide directed traps completely.
