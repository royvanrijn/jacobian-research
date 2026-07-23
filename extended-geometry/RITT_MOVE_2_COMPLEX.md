# The Ritt move 2-complex and the degree-thirty braid

Work over a field of characteristic zero, with monic-original polynomial
factors.  This note upgrades the relation graph used in
[the general Hessian--Ritt intersection programme](GENERAL_HESSIAN_RITT_INTERSECTIONS.md)
to a coefficient-decorated 2-complex.  The upgrade is forced already in
degree thirty: the two reduced paths through the `S_3` braid have the same
reduced image and the same normalization, but they do not define the same
scheme on the standard complete-decomposition chart.

The calculation is exact.  It distinguishes four notions which a graph alone
conflates:

1. equality of set-theoretic images;
2. equality of reduced image schemes;
3. equality after normalization;
4. equality of the unreduced path incidence schemes.

Only the first three agree in the degree-thirty braid.

## 1. The combinatorial 2-complex

Fix a multiset of indecomposable degrees.  The **Ritt move 2-complex** has:

* a vertex for every normalized complete decomposition;
* an edge for every elementary adjacent power or Chebyshev Ritt move;
* a commuting square for
  \(s_i s_j=s_j s_i\), \(\lvert i-j\rvert>1\);
* a braid hexagon for
  \(s_i s_{i+1}s_i=s_{i+1}s_is_{i+1}\).

The 2-cells identify move paths, not vertices.  This distinction matters:
identifying all vertices in a connected graph remembers only Ritt
connectivity, whereas filling its Coxeter relations records which
identifications between paths are coherent.

For pairwise-coprime distinct degrees the implementation in
[`jcsearch/ritt_complex.py`](../jcsearch/ritt_complex.py) constructs this
Coxeter 2-skeleton.  As a structural regression, the degrees `(2,3,5,7)`
give the boundary of the three-dimensional permutohedron:

\[
  24\ {\rm vertices},\qquad 36\ {\rm edges},\qquad
  6\ {\rm commuting\ squares}+8\ {\rm braid\ hexagons},
\]

and hence Euler characteristic \(24-36+14=2\).

## 2. Coefficient systems, not bare cells

A vertex with degree word \(\mathbf d=(d_1,\ldots,d_r)\) carries its affine
space of normalized factors, of dimension

\[
 \sum_i(d_i-1).
\]

An edge carries the power or Chebyshev coefficient correspondence between
the two endpoint charts.  A path therefore means an iterated fiber product
of correspondences.  A 2-cell should supply a comparison between the two
path correspondences around its boundary.

This is naturally a correspondence-valued weak 2-functor rather than an
ordinary graph-valued invariant.  There are three possible strengths for the
comparison attached to a 2-cell:

* equality of reduced dense images;
* isomorphism of normalizations;
* equality, or derived equivalence, of the full path schemes.

The degree-thirty calculation below proves that the third cannot be imposed
automatically.

On the Dickson cell there is a common coefficient space
\(\mathbb A^2_{t,z}\).  For a word \(\mathbf d\), put
\(s_i=\prod_{j>i}d_j\), \(c_i=D_{s_i}(t,z)\), and

\[
 F_i(X)=D_{d_i}(X+c_i,z^{s_i})-D_{d_i}(c_i,z^{s_i}).
                                                               \tag{2.1}
\]

Every \(F_i\) is monic and original, and the Dickson composition identity
gives

\[
 F_1\circ\cdots\circ F_r
 =D_{\prod d_i}(W+t,z)-D_{\prod d_i}(t,z).                     \tag{2.2}
\]

Thus (2.1) is an explicit coefficient map from the same affine plane to
every vertex and every Chebyshev edge correspondence in the cell.

## 3. The degree-thirty braid

For `(2,3,5)` the one-skeleton is a hexagon with six vertices and six edges.
Its unique 2-cell compares

\[
\begin{aligned}
\gamma_L:\quad&
235\longrightarrow325\longrightarrow352\longrightarrow532,\\
\gamma_R:\quad&
235\longrightarrow253\longrightarrow523\longrightarrow532.
\end{aligned}                                                  \tag{3.1}
\]

Here a word is ordered outer to inner.  Both paths are carried by the
Dickson polynomial

\[
 D_{30}(W+t,z)-D_{30}(t,z).                                   \tag{3.2}
\]

The filled hexagon has Euler characteristic \(6-6+1=1\).  At the level of
normalizations, (2.1) makes (3.1) a genuinely coherent braid.

## 4. Path ideals on the `2 o 3 o 5` chart

Let

\[
\begin{aligned}
P(W)&=W^2+p_1W,\\
Q(W)&=W^3+q_2W^2+q_1W,\\
R(W)&=W^5+r_4W^4+r_3W^3+r_2W^2+r_1W,
\end{aligned}
\]

and work in

\[
 S=\mathbb Q[p_1,q_1,q_2,r_1,r_2,r_3,r_4]
\]

with \(f=P\circ Q\circ R\).  The base chart automatically carries the cuts
of outer degree `2` and `6`.  For \(d\in\{3,5,10,15\}\), let \(J_d\) be the
pullback of the canonical residual ideal for the cut
\(\mathcal C_{d,30/d}\).

The endpoint, two path, and full-boundary ideals are

\[
\begin{aligned}
 I_E&=J_5+J_{15},\\
 I_L&=J_3+J_5+J_{15},\\
 I_R&=J_5+J_{10}+J_{15},\\
 I_{\partial}&=J_3+J_5+J_{10}+J_{15}.                         \tag{4.1}
\end{aligned}
\]

Let \(K\) be the graph ideal of the parametrization

\[
\begin{aligned}
r_4&=5t,&r_3&=10t^2-5z,&r_2&=10t^3-15zt,\\
r_1&=5t^4-15zt^2+5z^2,&
q_2&=3D_5(t,z),\\
q_1&=3(D_5(t,z)^2-z^5),&
p_1&=2D_{15}(t,z).
\end{aligned}                                                  \tag{4.2}
\]

The inverse coordinates are already polynomial:

\[
 t={r_4\over5},\qquad
 z={2r_4^2-5r_3\over25}.                                      \tag{4.3}
\]

Consequently \(S/K\cong\mathbb Q[t,z]\); the common reduced image is smooth,
normal, and has no ordinary normalization conductor.

## 5. Exact scheme comparison

Groebner reduction in `S` gives

\[
 \boxed{I_E=I_L\subsetneq K=I_R=I_{\partial}.}                 \tag{5.1}
\]

Moreover,

\[
 K^3\not\subset I_L,\qquad K^4\subset I_L.                    \tag{5.2}
\]

It follows at once that

\[
 \sqrt{I_L}=K.
\]

Thus the two paths have the same reduced image and the same normalization,
but different nilpotent structures.  The first path has nilpotence index
four; the second path and the full 2-cell boundary are already reduced.
Adding the missing degree-`10` incidence kills the entire thickening, while
adding the degree-`3` incidence to the endpoints does nothing
scheme-theoretically.

The nilpotent discrepancy is not generic on the Dickson surface.  If
\(\mathcal N=K/I_L\), exact ideal quotient gives

\[
 \operatorname{Ann}_{S/K}(\mathcal N)
 =((2r_4^2-5r_3)^2)=(z^2),\qquad
 \sqrt{\operatorname{Ann}_{S/K}(\mathcal N)}=(z).             \tag{5.3}
\]

More precisely, the checker verifies

\[
 (I_L:K)+K=K+((2r_4^2-5r_3)^2).
\]

For a nonreduced scheme the map to the normalization factors through its
reduction and is not an inclusion of coordinate rings.  Hence (5.3) is best
called the **normalization-defect annihilator**, or conductor analogue, not
the ordinary conductor of a birational ring inclusion.  It is the doubled
monomial divisor scheme-theoretically, with reduced support exactly the
degeneration \(z=0\) where the Dickson cell meets the power cell.

The tangent calculation sees the same asymmetry:

| point on \(\mathbb A^2_{t,z}\) | \(\dim T(S/I_L)\) | \(\dim T(S/I_R)\) |
|---|---:|---:|
| generic test point `(1,2)` | 2 | 2 |
| generic point of `z=0`, tested at `(1,0)` | 3 | 2 |

The reduced Dickson surface has tangent dimension two everywhere.  Therefore
the left path has exactly one excess first-order direction along the
monomial divisor, despite its nilpotence persisting to order four.

Two zero-dimensional slices show that this is not a product thickening.
Write \(h=-2r_4^2+5r_3=-25z\), and fix \(t=1\) by \(r_4=5\).  Exact vector
space dimensions and tangent dimensions give

\[
\begin{array}{c|c|c|c}
\text{slice ideal}&\text{length}&\text{embedding dimension}&\text{algebra}\\
\hline
I_L+(h,r_4-5)&2&1&\mathbb Q[\varepsilon]/(\varepsilon^2),\\
I_L+(I_L:K,r_4-5)&5&1&\mathbb Q[\varepsilon]/(\varepsilon^5).
\end{array}                                                    \tag{5.4}
\]

There is no conflict between (5.2) and the length-five conductor slice.  The
nilradical \(K/I_L\) has index four, while in the conductor slice the reduced
transverse parameter also becomes nilpotent.  In the length-five algebra the
image of \(K\) is \((\varepsilon^2)\): the successive quotient lengths by
\(K,K^2,K^3,K^4\) are \(2,4,5,5\).  Equation (5.4) therefore shows that the
order-four defect is coupled to the doubled base divisor; it is not the
product of the reduced divisor with
\(\mathbb Q[\varepsilon]/(\varepsilon^4)\).

The same exact calculation on the opposite `5 o 3 o 2` coefficient chart
gives identical results.  The path

\[
 532\longrightarrow352\longrightarrow325\longrightarrow235
\]

which omits the composite cut `10` is nonreduced with nilpotence index four,
whereas the path through `523` and `253`, which omits the prime cut `3`, and
the full boundary are reduced.  If `n_1` is the inner quadratic coefficient
and `m_2` the middle cubic quadratic coefficient, then

\[
 t={n_1\over2},\qquad
 z={3n_1^2-4m_2\over24},
\]

and the defect annihilator restricts to
\(((3n_1^2-4m_2)^2)=(z^2)\).  Thus the ideal-theoretic asymmetry is unchanged
when either endpoint of the braid is used as the normalization chart.

Repeating the calculation at the other four vertices gives the complete
scheme table:

| opposite endpoint charts | omitted composite cut | omitted prime cut on the other path | nilpotence index | defect annihilator |
|---|---:|---:|---:|---|
| `235`, `532` | 10 | 3 | 4 | `(z^2)` |
| `325`, `523` | 15 | 2 | 3 | `(z^2)` |
| `352`, `253` | 6 | 5 | 4 | `(z^4)` |

In every row the path omitting the composite cut is the nonreduced path and
the path omitting the prime cut is the reduced Dickson surface.  Both
endpoint charts give the same exponent and annihilator.  Consequently the
reduced support \(z=0\), the reduced/thick orientation, and exact
polynomial--Hessian synchronization are invariant under changing endpoint
chart.  The higher scheme structure is nevertheless **sector-dependent**:
it is attached to which composite cut the half-braid omits, not merely to
the underlying unlabelled hexagonal 2-cell.

The annihilator slices distinguish the three labelled sectors even more
sharply:

| omitted cut | exact slice algebra | Hilbert vector | \(K\)-adic quotient lengths |
|---:|---|---|---|
| 10 | \(\mathbb Q[u]/(u^5)\) | `(1,1,1,1,1)` | `(2,4,5,5)` |
| 15 | \(\mathbb Q[u,v]/(u^2,v^2)\) | `(1,2,1)` | `(2,4,4,4)` |
| 6 | \(\mathbb Q[u,v]/(u^4,v^2)\) | `(1,2,2,2,1)` | `(4,7,8,8)` |

These are exact presentations, not merely algebras with the displayed
Hilbert vectors.  In the cut-`15` chart the eliminated ideal is
\((a^2,b^2-6ab)\), and the change \(v=b-3a\) gives the table.  In the
cut-`6` chart it is

\[
 \left(a^4,\,
 500b^2-1050a^2b-53a^3b\right);
\]

completing the square modulo \(a^4\) gives
\(\mathbb Q[a,v]/(a^4,v^2)\).  All three slices are Artin Gorenstein with
one-dimensional socle.

This also supplies the first derived decoration.  Their conormal modules
over the smooth transverse ambient spaces are free of ranks `1,2,2`, and
their residue-field Koszul Tor ranks are respectively `(1,1)`, `(1,2,1)`,
and `(1,2,1)`.  The internal generator orders `(5)`, `(2,2)`, and `(4,2)`
distinguish the sectors even when the ungraded Tor ranks agree.  Thus the
transverse local algebra and its filtered conormal data are labelled 2-cell
data: neither can be recovered from the unlabelled Coxeter hexagon or the
common normalization alone.

More intrinsically, the comparison with the reduced path is the augmentation
\(A\to\mathbb Q\).  Its kernels have lengths `4`, `3`, and `7`.  At the
augmentation, the two-term cotangent complex of \(A/\mathbb Q\) has
\((\dim H_0,\dim H_1)=(1,1),(2,2),(2,2)\): every presentation relation has
order at least two, so its Jacobian differential vanishes at the closed
point.  The intrinsic residue-field Poincare series are

\[
 P_{\mathbb Q}^{A}(T)=\frac1{1-T},\qquad
 \frac1{(1-T)^2},\qquad \frac1{(1-T)^2}.                       \tag{5.4}
\]

Thus the intrinsic Tor ranks begin `(1,1,1,...)` in the cut-`10` sector and
`(1,2,3,4,...)` in the other two.  Again the ungraded series does not
distinguish cuts `15` and `6`; their relation orders and conductor filtrations
do.  This completes the transverse derived comparison.  What remains is to
glue it along the non-transverse \(z=0\) divisor and across adjacent cells.

## 6. Ambient polynomial and Hessian tangent spaces

The tangent asymmetry is not an artifact of pulling the incidence ideals
back to the `2 o 3 o 5` chart.  For every cut \(d\mid30\), differentiate the
normalized composition map

\[
 \mathcal P_d\times\mathcal P_{30/d}\longrightarrow\mathcal P_{30}
\]

at the Dickson factorization.  Intersect the resulting six tangent images
inside the 29-dimensional monic-original coefficient space.  Repeating the
calculation after deleting the linear coefficient gives the corresponding
28-dimensional Hessian coefficient space.  The two calculations agree:

| point | endpoint cuts | \(\gamma_L\) cuts | \(\gamma_R\) cuts | all six cuts |
|---|---:|---:|---:|---:|
| `(t,z)=(1,2)` | 2 | 2 | 2 | 2 |
| `(t,z)=(1,0)` | 3 | 3 | 2 | 2 |

Every individual normalized composition map has its full expected tangent
rank in both ambient spaces.  Hence the excess direction on the left path
survives in the actual polynomial and Hessian incidence intersections.  The
degree-`10` cut removes it, while the degree-`3` cut does not.

There is also no first-order synchronization defect at either test point.
Deleting the linear coefficient is injective on each composition tangent
image, and the polynomial and Hessian common intersections have the same
dimension.  Therefore projection identifies the common polynomial tangent
space with the intersection of the projected Hessian tangent spaces.  Any
unsynchronized Hessian component must be nonlinear or supported away from
these certified Dickson opens.

On all six complete-decomposition charts the comparison is stronger than
first order.  Restore
the degree-one coefficient residual in the canonical reconstruction for
every requested cut.  For the endpoint, both paths, and the full boundary,
Groebner reduction gives

\[
 I^{\rm polynomial}_{\bullet}
 =I^{\rm Hessian}_{\bullet}.                                  \tag{6.1}
\]

Thus the order-four braid thickening is already present in the synchronized
ordinary-composition intersection.  Forgetting the linear coefficient
neither creates nor enlarges it on these charts.  Equation (6.1) is an exact
scheme-theoretic Hessian-transfer certificate for this braid component,
though not a proof of transfer for components outside the
complete-decomposition atlas.

Cycling the braid gives a uniform prime/composite rule.  At a generic
Dickson point, any five of the six cuts have two-dimensional tangent
intersection.  At `z=0`, omit one cut at a time:

| omitted outer degree | 2 | 3 | 5 | 6 | 10 | 15 |
|---|---:|---:|---:|---:|---:|---:|
| tangent dimension | 2 | 2 | 2 | 3 | 3 | 3 |

Thus every prime outer-degree cut is first-order redundant at the monomial
point, whereas each composite outer-degree cut kills one braid-excess
direction.  This holds identically before and after Hessian projection and
describes all three rotations of the `S_3` braid, not only the displayed
choice of opposite vertices.

## 7. Consequences for the conjecture

The reduced component statement suggested by the relation-graph theorem
survives this test:

> **Reduced Ritt-complex conjecture.**  Irreducible components of multiple
> Hessian-composition intersections are dense reduced images of connected
> coefficient-decorated sub-2-complexes of the Ritt move complex, after paths
> related by its 2-cells are identified at the level of normalization.

The unreduced enhancement needs more structure.  A braid 2-cell cannot mean
literal equality of its two path fiber products.  It must retain at least
the nilpotent comparison module, and a derived formulation should retain the
corresponding excess conormal or Tor class.  In degree thirty the coherence
defect is supported on the power--Chebyshev overlap \(z=0\).

The exact ideal and Hessian-transfer comparisons are now verified on all six
complete-decomposition charts.  Independently, five factor-chart
ideal-membership certificates form the cut spanning tree
`2-6-3-15-5-10`, so the all-six degree-thirty Hessian intersection is
globally synchronized even on components missing this atlas.  This is not
yet a global component classification: the calculation does not exclude
additional synchronized off-chart components and does not prove the general
synchronization conjecture for arbitrary subintersections.  What it does
prove is that any global formulation which remembers scheme structure needs labelled 2-cells
together with nontrivial coherence data.

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_degree30_ritt_2_complex.py
```

The checker constructs the Coxeter 2-complex, verifies all six Dickson
vertex coefficient maps and edge identities, builds the four ideals in
(4.1) on all six vertex charts, certifies (5.1)--(5.4) and the complete
sector table, and computes both path intersections in the ambient polynomial
and Hessian tangent spaces exactly.
