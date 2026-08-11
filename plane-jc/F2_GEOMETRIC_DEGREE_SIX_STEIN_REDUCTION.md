# F2 geometric-degree-six Stein and cubic-germ reduction

> **Status.**  This note proves an unconditional localization theorem at
> geometric degree six and a complete classification under the additional
> hypothesis that the terminal Cartier slice in the local degree-three Stein
> factor is normal (equivalently here, conductor-free).
> Every affine nonproperness component punctures the target boundary at the
> unique terminal value `125/729`, and all of its local branching is carried
> by the single three-sheet packet there.  For a normal terminal slice the germ
> has the one-parameter cubic normal form
> `z=w^3+a(pi)w`.  It has exactly three possible reduced ramification
> patterns: one `k=1` three-cycle row, two `k=1` transposition rows, or one
> `k=2` transposition row.  In particular a mixed transposition/three-cycle
> pair cannot occur.  Every one of these patterns reproduces the already
> known terminal logarithmic cokernel `R/(w^3)`, so that cokernel does not
> distinguish them.  Nonnormal terminal slices remain open, even when the
> ambient Stein surface is smooth; their positive normalization conductor is
> precisely the
> point information missed by the terminal determinant packet.  Thus this is
> a finite reduction, not an exclusion of degree six or of `(75,125)`.

The residue arithmetic, cubic discriminants, parity classification,
endpoint blowup models, and three basic conductor orders are replayed by
[`verify_f2_geometric_degree_six_stein_reduction.py`](../scripts/verify_f2_geometric_degree_six_stein_reduction.py).

## 1. The terminal fiber exhausts degree six

Let

\[
 F:\mathbb A^2\longrightarrow\mathbb A^2
\]

be a hypothetical F2 `(75,125)` Keller map of geometric degree `d=6`, and
let

\[
 X\longrightarrow Z\mathop{\longrightarrow}^{\bar F}Y          \tag{1.1}
\]

be the Stein factorization on a proper resolved compactification.  The
finite map `bar F` has degree six.  The certified terminal divisor maps to
the target boundary component with residue cover

\[
 h(s)=\frac{125s(s+1)^5}{(9s^2+15s+5)^3}.                       \tag{1.2}
\]

Its finite nonzero branch value and its fiber partition are

\[
 \lambda_0=\frac{125}{729},\qquad
 h^{-1}(\lambda_0)=(3,1,1,1).                                  \tag{1.3}
\]

At every other finite nonzero value the partition is `(1^6)`.

Fix a target-boundary point `y_lambda` with residue `lambda`.  For every
point `z_i` of the finite Stein fiber above it, let `N_i` be the local degree
of `bar F`.  A terminal residue point of local index `q_i` contributes at
least `q_i` to `N_i`.  Therefore

\[
 6=\sum_iN_i\ge\sum_iq_i=6.                                   \tag{1.4}
\]

Every inequality in (1.4) is an equality.  There are no additional local
sheets over the point, and each terminal residue point has local surface
degree exactly its residue index.

There is one separation point hidden in this statement.  Two distinct
points of the terminal normalization cannot be identified by `X->Z`.  Such
an identification would put both points in one connected Stein fiber, hence
would require a connected chain of contracted curves joining them.  Every
curve over the target boundary lies in the source boundary.  But this
boundary is obtained from the line at infinity by point blowups, so its dual
graph is a tree; a chain joining two distinct points of the same horizontal
terminal vertex would create a cycle.  Therefore the terminal points remain
distinct in `Z`.

If `lambda` is finite, nonzero, and different from `lambda_0`, the six local
degrees are all one.  A finite flat local algebra of rank one is the base
ring, so all six germs are etale.  No affine branch divisor can meet the
target boundary there.  Consequently:

\[
 \boxed{\text{every affine nonproperness component punctures at }
        \lambda_0=125/729.}                                   \tag{1.5}
\]

At `lambda_0`, completion of the finite algebra splits it into three rank-one
etale factors and one local rank-three factor.  Hence every affine branch
germ lies in the branch locus of that same cubic factor.  The three simple
sheets cannot carry it.

This argument uses only finite-degree conservation and the terminal
passport.  It does not assume that the cubic Stein point is smooth.

## 2. Immediate global consequences

A meridian of an irreducible branch divisor in a three-sheet local cover is
either a transposition or a three-cycle.  Its moved-sheet cost is therefore
`2` or `3`.  Orevkov's residue-degree identity gives, at `d=6`,

\[
 \sum_j A_j\le d-1=5.                                          \tag{2.1}
\]

Thus there are at most two affine nonproperness components, and two
three-cycle components are impossible.  Before using the local algebra, the
only numerical multisets are

\[
 (2),\quad(3),\quad(2,2),\quad(2,3).                           \tag{2.2}

Section 3 removes the last multiset whenever the terminal Cartier slice of
the cubic Stein factor is normal.

## 3. Normal terminal Stein slice: complete local classification

Assume now that the terminal Cartier slice of the local degree-three Stein
factor is normal.  Its completed local ring is `k[[w]]`, with target residue
parameter `z=w^3`.  Since this regular Cartier section lies in a
two-dimensional Cohen--Macaulay local ring, the Stein point itself is
regular and `pi` is a source parameter.  The characteristic-zero
one-parameter `A_2` normal form then gives formal or analytic coordinates in
which

\[
 \boxed{(\pi,w)\longmapsto
        (\pi,z=w^3+a(\pi)w),\qquad a(0)=0.}                    \tag{3.1}
\]

A target translation has removed the constant term.  The ramification and
branch equations are

\[
 3w^2+a(\pi)=0,qquad
 \Delta(\pi,z)=-4a(\pi)^3-27z^2=0.                            \tag{3.2}

There are three and only three cases.

### 3.1 Cyclic specialization

If `a=0`, the reduced ramification curve is `w=0`, with generic index three.
Its image is the smooth branch `z=0`, transverse to `pi=0`.  It is one
`k=1` component with cubic inertia:

\[
 (A,k)=(3,1).                                                  \tag{3.3}

### 3.2 Even order

Suppose `r=ord_pi(a)<infinity` is even.  After adjoining a harmless unit
square root, (3.2) has two smooth ramification branches

\[
 w=\mathord\pm c\pi^{r/2}.                                    \tag{3.4}

Each is simply ramified and each has `pi` itself as normalization parameter.
Their images are two distinct `k=1` target branches.  Thus the pattern is

\[
 (A_1,k_1)+(A_2,k_2)=(2,1)+(2,1).                             \tag{3.5}

### 3.3 Odd order

If `r=ord_pi(a)` is odd, the reduced ramification curve is irreducible.  Its
normalization has

\[
 \pi=t^2,\qquad w=t^r\cdot\text{unit},\qquad
 z=t^{3r}\cdot\text{unit}.                                   \tag{3.6}

It is generically simply ramified, while its image has boundary contact
`k=2`.  The pattern is

\[
 (A,k)=(2,2).                                                  \tag{3.7}

For `r=1` the image is the ordinary cusp `(pi,z)=(t^2,t^3)`.
For larger odd `r` it is the higher unibranch germ `(2,3r)`.

Combining the three cases gives the promised exhaustive smooth-point list:

\[
 \boxed{
 (3;k=1),\qquad(2,2;k=1,1),\qquad(2;k=2).}                    \tag{3.8}

In particular the numerical possibility `(2,3)` in (2.2) is not the branch
divisor of one smooth cubic germ.

## 4. Why the terminal `R/(w^3)` packet cannot choose a case

The exact terminal calculation gives a cyclic logarithmic cokernel
`R/(w^3)` at the endpoint over `lambda_0`.  All three cases in (3.8) have
that same resolved endpoint.

For `a=0`, the terminal and ramification curves are already transverse.  In
logarithmic bases the map is diagonal with last entry `3w^3`.

For `ord(a)=1`, the ramification curve is tangent to the terminal curve.
After the two point blowups that separate them, a chart at the
terminal--exceptional node is

\[
 \pi=vw^2,qquad z=(1+v)w^3.                                  \tag{4.1}

The logarithmic determinant is

\[
 w^3(3+v)\cdot\text{unit}.                                   \tag{4.2}

For `r=ord(a)>=2`, one blowup separates the terminal tangent direction from
the ramification direction.  At the terminal--exceptional node,

\[
 \pi=vw,qquad
 z=w^3+v^rw^{r+1},                                            \tag{4.3}

and the logarithmic determinant is

\[
 w^3\left(3+v^rw^{r-2}\right).                               \tag{4.4}

Each parenthesis is a unit.  Since the first logarithmic row is unimodular,
all cases reduce to

\[
 \boxed{\operatorname{coker}(d\bar F^{\log})\simeq R/(w^3).} \tag{4.5}

This explains structurally why the terminal determinant and its generic
Smith data could not reveal the missing affine branch.

On the current common fan, the immediate component on the carrier side of
the terminal divisor is a `(-1)` component of valency two.  A smooth cubic
cusp attachment would meet it at a second point and make it trivalent.  This
is exactly the sole local shape not removed by the existing contracted
divisor gate.  This incidence match is a necessary compatibility, not a
construction of the global cover.

## 5. Nonnormal terminal slices are a conductor problem

The finite local cubic algebra is flat: a normal two-dimensional local
domain finite over the regular target is Cohen--Macaulay, hence free of rank
three.  Its closed length-three local fiber has one of the two algebra
types

\[
 k[\epsilon]/(\epsilon^3),
 \qquad
 k[u,v]/(u,v)^2.                                               \tag{5.1}

That dichotomy does **not** decide normality of the terminal slice.  The first
two rows below both have curvilinear closed fiber, although only the first is
normal.

The restriction to the target boundary makes the normalization parameter
visible.  Put `z=w^3`.  Three basic rank-three orders illustrate the missing
information:

\[
\begin{array}{c|c|c|c}
\text{boundary order}&\text{basis over }k[[z]]&
 \delta=\ell(\widetilde A/A)&\operatorname{ord}_z\operatorname{disc}\\ \hline
k[[w]]&(1,w,w^2)&0&2\\
k[[w^2,w^3]]&(1,w^2,w^4)&1&4\\
k[[w^3,w^4,w^5]]&(1,w^4,w^5)&2&6.
\end{array}                                                     \tag{5.2}

The last closed fiber is square-zero.  The discriminant order rises by
`2delta`, while the normalized residue map remains `z=w^3` in every row.
There are higher-conductor orders as well; (5.2) is diagnostic, not an
exhaustive classification of singular surface sections.

### 5.1 Exact conductor--contact identity

The same discriminant calculation gives an exact classification without
classifying the order itself.  Let `A_0` be the terminal Cartier slice of
the local cubic algebra and let

\[
 A_0\subset \widetilde A_0=k[[w]],\qquad z=w^3,
 \qquad \delta_T=\ell_k(\widetilde A_0/A_0).                  \tag{5.3}
\]

Both rings are free rank-three lattices over the DVR `k[[z]]`.  The
normalized lattice has discriminant `-27z^2`, while passage to a sublattice
of colength `delta_T` multiplies its discriminant by the square of the
index.  Therefore

\[
 \operatorname{ord}_z\operatorname{disc}(A_0/k[[z]])
 =2+2\delta_T.                                                \tag{5.4}
\]

The boundary itself is not a branch component: its generic normal index is
one.  Hence intersecting the discriminant divisor of the finite cubic map
with the target boundary counts precisely the affine branch germs through
the terminal value.  If `C_j` has boundary contact `k_j` and generic tame
inertia index `e_j` inside the cubic packet, its discriminant coefficient is
`e_j-1`.  Consequently

\[
 \boxed{\sum_j(e_j-1)k_j=2+2\delta_T.}                        \tag{5.5}
\]

Combining (5.5) with the four inertia multisets in (2.2) gives the complete
degree-six conductor--contact atlas:

\[
\begin{array}{c|c}
\text{generic inertia}&\text{exact contact equation}\\ \hline
(3)&k=1+\delta_T\\
(2)&k=2+2\delta_T\\
(2,2)&k_1+k_2=2+2\delta_T\\
(3,2)&2k_3+k_2=2+2\delta_T.
\end{array}                                                    \tag{5.6}
\]

Here every `k_j` lies in `1,...,24` by affine purity.  In particular the
first family has `0<=delta_T<=23`, the second has
`0<=delta_T<=11`, and the transposition contact `k_2` in every mixed row is
even.  The mixed row also spends the entire Orevkov budget
`3+2=d-1`, so every residue-defect term in that identity vanishes.

At `delta_T=0`, (5.6) gives exactly the three normal-slice patterns (3.8):
`(3;k=1)`, `(2;k=2)`, and `(2,2;k_1=k_2=1)`; the mixed family is empty.
Thus positive conductor is not merely one possible explanation of an extra
contact.  Equation (5.5) measures *all* of the missing boundary contact
exactly.

Ambient smoothness does not remove this regime.  The exact cyclic cubic map

\[
 (x,y)\longmapsto(u=x^3-y^2,\ v=y)                            \tag{5.7}
\]

has smooth source and target.  The target boundary `u=0` pulls back to the
cusp `x^3=y^2`, whose normalization is `(x,y)=(t^2,t^3)`.  Hence its
boundary order is `k[[t^2,t^3]]`, its terminal residue map is still
`v=t^3`, and its normalization conductor has length one.  The ramification
curve `x=0` has generic index three and maps to

\[
 u=-v^2,                                                       \tag{5.8}
\]

a smooth branch of contact `k=2` with the target boundary.  Thus even a
smooth Stein point can carry the extra pattern `(A,k)=(3,2)` once its
terminal slice is nonnormal.  More generally `x^3=u+v^m` produces higher
contact orders when the boundary slice is irreducible.  These models are
local compatibility examples, not Keller maps.

This is exactly the normalization mismatch anticipated by the logarithmic
perfect-complex programme.  The terminal `R/(w^3)` packet sees the normalized
degree-three branch but not the conductor of the contracted Stein slice.
Any exclusion of the conductor regime must either:

1. recover and bound this conductor through the localized `ch_2`/`Fitt_1`
   filtration; or
2. classify the normal rank-three surface algebras compatible with the
   compiled carrier chain and exclude their resolution graphs.

## 6. Exact remaining degree-six ledger

The degree-six problem is now divided as follows.

| local cubic regime | target pattern | present disposition |
| --- | --- | --- |
| normal slice, `a=0` | one `k=1` cubic-inertia component | covered by the existing one-component degree-six cubic-row exclusion |
| normal slice, `ord(a)` even | two `k=1` simple-inertia components | unresolved multi-component gluing case |
| normal slice, `ord(a)` odd | one `k=2` simple-inertia component | unresolved cusp-at-infinity case |
| nonnormal terminal slice | one of the four finite conductor--contact families (5.6) | unresolved normalization/contracted-fiber case |

Thus a degree-six counterexample, if one exists, must lie in one of the last
three rows.  Conversely, excluding those three rows closes geometric degree
six completely.  The last row includes smooth ambient Stein points such as
(5.7); it must not be shortened to “singular Stein point.”

## Sources

- Rick Miranda,
  [*Triple Covers in Algebraic Geometry*](https://www.math.colostate.edu/~miranda/preprints/TripleCoversInAG.pdf),
  Sections 2--5, for finite flat cubic algebras, the discriminant, and the
  smooth total-ramification cusp model.
- S. Yu. Orevkov,
  [*On three-sheeted polynomial mappings of C2*](https://www.math.univ-toulouse.fr/~orevkov/jc86.pdf),
  Lemma 4.2, for the global local-degree budget.

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_geometric_degree_six_stein_reduction.py
```
