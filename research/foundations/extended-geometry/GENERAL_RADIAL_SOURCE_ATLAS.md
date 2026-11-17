# The general radial source atlas

## Result

The source-cover comparison is now uniform on every nonresonant radial
stratum of a simultaneous polynomial root-cluster degeneration.

Let the labelled clusters have multiplicities

\[
 \mu_1,\ldots,\mu_r,
 \qquad N=\sum_i\mu_i,
\]

and local forms

\[
 H(c_i+x_iX)
 =
 x_i^{\mu_i}u_i
 \prod_{a=1}^{\mu_i}(X-\alpha_{ia})
 +O(x_i^{\mu_i+1}),\qquad u_i\ne0.                  \tag{0.1}
\]

The branch scale of cluster \(i\) is
\(\mu_i v(x_i)\).  An ordered set partition of the clusters records which
branch scales are equal and their order.  On the target bubble where cluster
\(i\) becomes active, its complete source contribution is forced:

\[
\begin{array}{c|c|c}
\text{target level} & \text{source contribution} & \text{degree}\\ \hline
\text{before the active level} & S=V^{\mu_i}
  \text{ power connector} & \mu_i\\
\text{at the active level} &
  S=u_i\prod_a(X-\alpha_{ia})
  \text{ local polynomial tail} & \mu_i\\
\text{after the active level} &
  \mu_i\text{ identity strands }S=V & \mu_i.
\end{array}                                         \tag{0.2}
\]

Thus every target bubble has total source degree \(N\), every component
satisfies Riemann--Hurwitz, and every target node has an explicit index
partition of \(N\).  The degree-six thirteen-type atlas is the special case
\((\mu_1,\mu_2,\mu_3)=(2,2,2)\).

The companion
[polynomial monodromy-forest theorem](POLYNOMIAL_MONODROMY_FORESTS.md)
determines the source dual graphs and node indices on the nonradial
simple-branch resonance refinements as well.  The
[recursive resonance atlas](RECURSIVE_RESONANCE_ATLAS.md) supplies the
algebraic placement and contraction descent of their target-flag divisors.
Abstract existence and the coarse comparison also follow independently by
[finite normalization](SOURCE_GRAPH_FINITE_NORMALIZATION.md).

## 1. Derivation of the component rule

The factors outside cluster \(i\) are units at \(c_i\).  Rescaling
\(W=c_i+x_iX\) therefore gives (0.1).  The finite critical values of the
local degree-\(\mu_i\) polynomial all have scale \(x_i^{\mu_i}\), proving
the weight

\[
 w_i=\mu_i v(x_i).                                  \tag{1.1}
\]

At its active target scale, (0.1) is the local polynomial tail in (0.2).
Its zero divisor is the labelled cluster-root divisor

\[
 D_{i,0}=\sum_a[\alpha_{ia}],
\]

with multiplicities if the local roots collide, and its pole divisor is
\(\mu_i[\infty]\).  The
[source-vertex rigidity theorem](SOURCE_VERTEX_RIGIDITY.md) reconstructs the
map from those two fibers; the chosen target residue scale fixes \(u_i\).

Before this scale, only the totally ramified attaching behavior is visible.
Both node fibers are \(\mu_i\)-fold, so two-fiber rigidity gives the power
connector \(S=V^{\mu_i}\).  After the scale has passed, the \(\mu_i\)
branches near the roots are separated and contribute \(\mu_i\) unramified
identity strands.

There is no additional vertexwise Hurwitz choice in any of the three cases.

## 2. Degree and Riemann--Hurwitz

A power connector of degree \(\mu\) is totally ramified over both ends, so
its ramification contribution is

\[
 (\mu-1)+(\mu-1)=2\mu-2.
\]

A local polynomial tail is totally ramified at infinity and has total finite
critical ramification \(\mu-1\).  Its contribution is again
\(2\mu-2\).  An identity strand has contribution zero.  Hence every row of
(0.2) satisfies

\[
 R_E=2\deg(E\to P)-2.                                \tag{2.1}
\]

Each cluster contributes total degree \(\mu_i\) on every target bubble,
regardless of its state in (0.2).  Summing over clusters gives \(N\).

## 3. Node partitions and saturation

Let \(\ell_i\) be the active level of cluster \(i\).  At the target node
after level \(q\), cluster \(i\) contributes

\[
 \begin{cases}
  \mu_i\text{ indices equal to }1,&\ell_i\le q,\\
  \text{one index equal to }\mu_i,&\ell_i>q.
 \end{cases}                                        \tag{3.1}
\]

The indices always sum to \(N\).  If \(A_q\) is the set of clusters not yet
active, the number of normalized branches at that node is

\[
 \frac{\prod_{i\in A_q}\mu_i}
      {\operatorname{lcm}_{i\in A_q}(\mu_i)},        \tag{3.2}
\]

with value one when \(A_q\) is empty.  This is exactly the general
[labelled-node saturation formula](LABELLED_NODE_SATURATION.md).  Actual
inertia is still the inverse image of the diagonal phase subgroup under the
label-preserving cover-automorphism character; (3.2) alone is not an inertia
order.

For example, strict scales with multiplicities \((2,3,4)\) have successive
node partitions

\[
 (2,3,4),\qquad
 (1,1,3,4),\qquad
 (1,1,1,1,1,4),                                    \tag{3.3}
\]

and normalization branch counts \(2,1,1\).  Unequal multiplicities therefore
show why the lcm, rather than a common cluster multiplicity, is the natural
saturation invariant.

### Full-chain radial inertia

Normalization branch counts at separate nodes cannot be multiplied to
obtain inertia.  Let

\[
 B_0\mid\cdots\mid B_k
\]

be the ordered scale partition in the level order of (3.1), so \(B_0\) is
first active, and define

\[
 M_j=\operatorname{lcm}_{i\in B_j}\mu_i,\qquad
 L_j=\operatorname{lcm}_{i\in B_j\cup\cdots\cup B_k}\mu_i.
\]

Matching connector rotations from the labelled zero screen to the rigid
stationary root vertex gives

\[
 |I_{\rm radial}|
   =\prod_{j=0}^{k}\frac{L_j}{M_j}.                 \tag{3.4}
\]

The labelled zero fibers fix the local polynomial-tail actions, while the
stationary transposition tree has trivial centralizer; hence only the cyclic
power-connector rotations remain.  If
\(t_j\in\mathbf Z/L_j\) is the diagonal phase at the \(j\)-th node, a
cluster \(i\in B_\ell\) imposes
\(t_0+\cdots+t_\ell=0\pmod{\mu_i}\).  Given the earlier phases, there are
exactly \(L_j/M_j\) choices at level \(j\).  Equal multiplicities make all
factors one, but unequal multiplicities can leave inertia.  For example,
the strict orders \((2)\mid(3)\) and \((3)\mid(2)\) have orders \(3\) and
\(2\).

This formula is proved and checked by the
[recursive resonance atlas](RECURSIVE_RESONANCE_ATLAS.md); it replaces the
incorrect componentwise claim that a root label automatically kills every
connector rotation.

## 4. Radial face compatibility

Ordered set partitions are the faces of the first weighted braid
subdivision.  Refining a block inserts one target node.  Formula (3.1) gives
the same degree \(N\) on its two sides and equal indices on paired source
branches, so it is precisely the characteristic-monoid datum required for
the admissible smoothing

\[
 \tau=s_j^{e_j}.
\]

Coarsening deletes that target node and smooths its source preimages.
Consequently the harmonic source trees, component degrees, and saturated
node branches are compatible on the whole radial face poset.

This statement is deliberately limited to the radial poset.  On an
equal-scale face, leading critical-value residues may collide.  The
wonderful target theorem inserts further resonance bubbles, and their
source flag divisors are constructed and checked against the third-fiber
equations by the recursive resonance atlas.  Their dual graphs and node
indices are supplied by the monodromy-forest theorem.

## 5. Interface with the completed atlas

Combining the four current results gives:

1. the weighted braid walls \(w_i=w_j\) are the first target layer;
2. the complete target refinement is the normalized wonderful pullback;
3. every radial source tree and component map is given by (0.2);
4. every radial node normalization is given by (3.2); and
5. full-chain radial inertia is given by (3.4); and
6. the recursive resonance theorem constructs all later flag divisors and
   contractions and extracts their simultaneous characters; and
7. the H2 marked/unmarked quotient is formal over the resulting labelled
   graph.

In particular, neither radial source enumeration nor a fresh Hurwitz search
is needed in higher degrees.

## Reproduction

Run:

```bash
.venv/bin/python scripts/verify_general_radial_source_atlas.py
```

The checker covers multiplicities \(2\) through \(6\) for one through four
labelled clusters: 780 multiplicity profiles, 48,580 ordered scale types,
149,630 target bubbles, 149,630 target nodes, and 15,805 label
permutations.  It verifies degree, componentwise Riemann--Hurwitz, every node
partition, lcm saturation, and the unequal \((2,3,4)\) example.  It also
checks (3.4) independently by dynamic phase counting on all \(48{,}580\)
ordered types: \(42{,}158\) have nontrivial inertia, with 81 distinct orders
and maximum order \(3000\) in the bounded range.
