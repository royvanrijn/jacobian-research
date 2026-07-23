# The stable-target graph as a wonderful pullback

## Result

The target half of the recursive branch-scale fan conjecture is a theorem.

Pass to an etale splitting atlas on which the finite critical points, hence
their critical values

\[
 \lambda_1,\ldots,\lambda_b,
\]

are labelled.  Retain the fixed target marks \(0\) and \(\infty\).  The
smooth locus gives a rational map to the configuration space
\(M_{0,b+2}\).  Its stable-target graph is

\[
 B^{\mathrm{tgt}}
 =
 \operatorname{Norm}
 \overline{\Gamma\bigl(X^\circ\dashrightarrow
   \overline M_{0,b+2}\bigr)}.                      \tag{0.1}
\]

Because \(\overline M_{0,b+2}\) is the wonderful compactification of the
labelled genus-zero configuration, (0.1) is equivalently the normalization
of the principal component of the iterated blowups of the pulled-back
diagonal ideals.  Its cone complex is the saturated pullback of the tropical
\(\overline M_{0,b+2}\) tree fan.

Thus:

- monomial branch scales determine the first radial walls;
- equal leading residues are simply later diagonal centers in the same
  building set;
- nested resonances are indexed by nested collision subsets;
- permutation descent is canonical.

The **coarse source-cover** comparison is now formal: the
[finite-normalization theorem](SOURCE_GRAPH_FINITE_NORMALIZATION.md) shows
that the normalized polynomial admissible-cover closure is finite birational
over (0.1), hence has coarse space \(B^{\mathrm{tgt}}\).  No further coarse
source blowup occurs.  For an explicit stack presentation, the
[source-vertex rigidity theorem](SOURCE_VERTEX_RIGIDITY.md) now shows that
there is no further Hurwitz-class choice once the source tree and the
weighted inverse-image divisors of the target flags are fixed: two fibers
determine each rational component map up to scale, and a third fiber fixes
the scale.  Thus the wonderful pullback does not by itself construct the
polynomial source enhancement, but the missing datum is precise.
The [general radial source atlas](GENERAL_RADIAL_SOURCE_ATLAS.md) constructs
that enhancement for every ordered first-scale stratum and arbitrary cluster
multiplicities.  The
[polynomial monodromy-forest theorem](POLYNOMIAL_MONODROMY_FORESTS.md)
determines the source dual trees, degrees, and node indices on every nested
simple-branch resonance stratum.  Only the algebraic positions, descent, and
simultaneous higher-codimension inertia characters of their target-flag
divisors remain; generic node inertia is given by
[monodromy centralizers](MONODROMY_INERTIA_CHARACTERS.md).

## 1. Boundary building set

Let \(n=b+2\), label the target marks by

\[
 0,1,\ldots,n-2,\infty,
\]

and distinguish \(\infty\).  A boundary divisor of
\(\overline M_{0,n}\) is represented uniquely by a subset

\[
 S\subset\{0,\ldots,n-2\},
 \qquad 2\le |S|\le n-2.                            \tag{1.1}
\]

It records the split \(S\mid S^c\), with \(\infty\in S^c\).  Two boundary
divisors meet precisely when

\[
 S\subset T,\qquad T\subset S,\qquad\text{or}\qquad S\cap T=\varnothing.
                                                               \tag{1.2}
\]

Hence boundary strata are nested sets, equivalently stable target trees.
The number of boundary divisors and zero-dimensional strata is

\[
 2^{n-1}-n-1,
 \qquad
 (2n-5)!!,                                          \tag{1.3}
\]

respectively.

The wonderful blowup description is standard; see Li,
[Wonderful compactification of an arrangement of
subvarieties](https://arxiv.org/abs/math/0611412), and the original
Fulton--MacPherson
[configuration-space compactification](https://annals.math.princeton.edu/1994/139-1/p06).
Kapranov's Chow-quotient description identifies the genus-zero stable-curve
model
([arXiv:alg-geom/9210002](https://arxiv.org/abs/alg-geom/9210002)).

## 2. Pullback ideals

On an affine target chart containing the fixed value \(0\), let
\(\lambda_i\) denote the moving branch values.  For a collision subset
\(S\), the diagonal ideal is

\[
 I_S=
 \begin{cases}
  (\lambda_i:i\in S\setminus\{0\}),&0\in S,\\[2mm]
  (\lambda_i-\lambda_{i_0}:i\in S\setminus\{i_0\}),
    &0\notin S,
 \end{cases}                                        \tag{2.1}
\]

where \(i_0\in S\) is arbitrary.  Changing \(i_0\) does not change the
ideal.

Construct the wonderful model by blowing up the diagonals in any
building-set order, for example increasing dimension after dominant
transforms.  At each stage, the graph of the induced rational lift from the
current principal component is the blowup of the pulled-back center ideal.
Iterating this graph construction, discarding components that do not
dominate \(X\), and normalizing at the end gives (0.1):

\[
 B^{\mathrm{tgt}}
 =
 \operatorname{Norm}\operatorname{Prin}
 \left(
   \operatorname{Bl}_{I_{S_k}}\cdots
   \operatorname{Bl}_{I_{S_1}}X
 \right).                                           \tag{2.2}
\]

Formula (2.2) is the algebraic version of the recursive initial-form
algorithm.  If the first nonzero terms of \(\lambda_i\) are monomials, the
normalization saturates their exponent lattice and produces the weighted
braid fan.  On an equality face, the dominant transforms of the larger
ideals in (2.1) are generated by differences of leading residues.  Their
blowups are exactly the later resonance modifications.

## 3. Tropical formulation

Let

\[
 \operatorname{trop}(\lambda):
 \Sigma_X\longrightarrow
 \mathbb R^b/\mathbb R(1,\ldots,1)
\]

be the tropical critical-value map on a toroidal root-cluster chart.  The
cone complex of (2.2) is

\[
 \Sigma_{B^{\mathrm{tgt}}}
 =
 \operatorname{Sat}
 \left(
   \Sigma_X
     \times_{\mathbb R^b/\mathbb R\mathbf1}
   \mathcal M_{0,b+2}^{\mathrm{trop}}
 \right).                                           \tag{3.1}
\]

This is a saturated fiber product, not merely the fan cut out by the first
hyperplanes \(v(\lambda_i)=v(\lambda_j)\).  The tropical tree fan retains
successive residue collisions automatically.

Existing tropical admissible-cover theory proves compatibility of the
classical and tropical branch maps; see
Cavalieri--Markwig--Ranganathan,
[Tropicalizing the space of admissible
covers](https://arxiv.org/abs/1401.4626).

## 4. Recovery of the two experiments

### Degree five

There are two relevant moving branch scales.  The target is
\(\overline M_{0,4}\cong\mathbf P^1\), and the target map is locally

\[
 [\lambda_1:\lambda_2]
 =
 [u_1x^3:u_2y^2].
\]

The graph of this rational map is the blowup of
\((u_1x^3,u_2y^2)\); units do not change its normalized blowup.  Hence

\[
 B^{\mathrm{tgt}}_5
 =\operatorname{Norm}\operatorname{Bl}_{(x^3,y^2)}X,
\]

recovering the degree-five obstruction.

### Degree six

With target marks \(0,\infty,A,B,C\), formula (1.3) gives ten boundary
divisors.  The six size-two subsets give

\[
 A=0,\ B=0,\ C=0,\ A=B,\ A=C,\ B=C.
\]

The four size-three subsets give the four triple centers

\[
 [1:0:0],\quad[0:1:0],\quad[0:0:1],\quad[1:1:1].
\]

Thus (2.2) is exactly the four-point blowup of \(\mathbf P^2\) found in the
sextic experiment.  The first three centers form the radial permutohedral
layer; the diagonal center is the triple-resonance refinement.  They were
not two unrelated constructions.

## 5. Equivariance and descent

Permuting the critical-value labels permutes the subsets \(S\), their ideals
\(I_S\), and every nested set.  Graph closure, principal component, and
normalization are functorial for these automorphisms.  Consequently the
unlabelled target graph is the finite quotient of (0.1); no separate
chartwise descent choices are required.

Combined with
[Labelled node saturation](LABELLED_NODE_SATURATION.md), the corrected
pipeline is now

\[
\boxed{
\begin{array}{c}
\text{normalized wonderful pullback of branch diagonals}
\\ \Downarrow\\
\text{construct the flag-complete polynomial source enhancement}
\\ \Downarrow\\
\text{reconstruct the unique source-vertex maps}
\\ \Downarrow\\
\text{apply node saturation and the required root stacks}
\\ \Downarrow\\
\text{take the formal marked/unmarked subgroup quotient}.
\end{array}}
\]

Finite normalization constructs the first arrow abstractly and identifies
its coarse space with the target graph.  Its combinatorial part is explicit
on both the radial and nested simple-branch resonance strata.  Algebraic
target-flag coordinates and all contraction maps are constructed by the
[recursive resonance atlas](RECURSIVE_RESONANCE_ATLAS.md).
Source-vertex reconstruction, node saturation, and the subgroup quotient
remain independent general theorems.

## 6. Explicit source-stack completion

Let \(B_N^{\mathrm{tgt}}\) be (0.1) and let
\(\Gamma_N^{\mathrm{ACV}}\) be the selected normalized graph in the
admissible-cover stack.  Forgetting the source cover gives

\[
 \Gamma_N^{\mathrm{ACV}}\longrightarrow B_N^{\mathrm{tgt}}. \tag{6.1}
\]

After selecting the generic polynomial component, normalizing its closure,
and passing to coarse spaces, (6.1) is an isomorphism by finite
normalization.  The recursive resonance theorem gives an explicit
presentation of the stack above that coarse graph.

For a source component over a stable target vertex, let \(D_h\) be the
weighted inverse-image divisor of an incident target flag \(h\).  The
source-vertex theorem proves that any two \(D_h\)'s determine the component
map up to one scalar and that a third flag fixes that scalar.  The explicit
**flag-complete source enhancement** atlas over the nonradial wonderful
boundary strata therefore consists of:

1. construction of the actual divisor positions \(D_h\) on the
   monodromy-forest source components in labelled root-residue coordinates;
2. the third-fiber residue equations on every source vertex;
3. compatibility of the divisor data under nested-set edge contraction;
4. componentwise Riemann--Hurwitz and node-index matching; and
5. comparison of the surviving cover automorphisms with the character test
   in the labelled-node theorem.

The recursive resonance theorem proves items 1--5: framed screens give the
root-residue coordinates, normalized initial forms give every \(D_h\),
affine composition proves edge contraction, and the full-centralizer fiber
product gives the actual simultaneous characters.  The degree-six central
square-cubic invariant and the radial/Maxwell atlases verify this package on
the `(2,2,2)` chart.  For general multiplicities, the radial components and
all simple-branch resonance dual graphs and node indices agree with the same
general construction.

## Reproduction

Run:

```bash
.venv/bin/python scripts/verify_branch_wonderful_pullback.py
.venv/bin/python scripts/verify_source_vertex_rigidity.py
.venv/bin/python scripts/verify_recursive_resonance_atlas.py
.venv/bin/python scripts/verify_branch_scale_fan.py
.venv/bin/python scripts/verify_degree_six_branch_target_graph.py
```

The first command enumerates all boundary divisors and maximal nested sets
for four through seven target marks, checks the full permutation action, and
recovers the degree-five weighted blowup and the degree-six six-line/four-
center construction from the same building set.
The second command checks the componentwise divisor reconstruction in
degrees one through seven and recovers the central, cluster-tail, and cyclic-
tail maps in the sextic atlas.
