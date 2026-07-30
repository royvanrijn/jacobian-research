# Jelonek's coefficient-space corollary and the deformation project

Work over `C`.  For integers `n,d>=1`, let

\[
 {\cal V}(n,d)=
 \bigl(\mathbb C[x_1,\ldots,x_n]_{\leq d}\bigr)^n
\]

be the affine coefficient space of polynomial self-maps of degree at most
`d`, and let

\[
 X(n,d)=\{F\in{\cal V}(n,d):\det DF=1\}                 \tag{1}
\]

be the coefficient scheme cut out by the coefficients of `det DF-1`.
Statements about irreducible components below concern the irreducible
components of the underlying reduced variety.  The scheme structure becomes
essential in the tangent-space project in Section 3.

## 1. Jelonek's component theorem

Jelonek proves that the subset

\[
 {\cal A}(n,d)=
 \{F\in X(n,d):F\text{ is a polynomial automorphism}\}
\]

is Zariski closed in `X(n,d)`.  Consequently every irreducible component
`Omega` of `X(n,d)` satisfies exactly one of the following:

1. `Omega` is contained in `A(n,d)`; or
2. `Omega\setminus A(n,d)` is a nonempty dense open subset, so a general
   point of `Omega` is a counterexample to the Jacobian conjecture.

He also proves

\[
 \dim\Omega\geq n^2-1.                                 \tag{2}
\]

These statements are Jelonek's
[Lemma 2.1 and Theorem 2.2](https://arxiv.org/abs/2607.20597),
not new results of this repository.

## 2. Immediate corollary for the explicit maps

Let

\[
 F_N:\mathbb A^3_{\mathbb Q}\longrightarrow\mathbb A^3_{\mathbb Q},
 \qquad N\geq3,
\]

be any of the explicit determinant-one weighted maps constructed in
[the all-degree rational-fiber theorem](../verified/ALL_DEGREE_RATIONAL_FIBERS.md).
It has a complete fiber of `N` distinct points, hence its complexification is
not injective and is not a polynomial automorphism.

Fix any coefficient bound `d>=deg(F_N)`, and regard `F_N` as a point of
`X(3,d)`.  If `Omega` is any irreducible component of `X(3,d)` containing
that point, then `Omega` cannot be contained in the closed automorphism
locus.  Jelonek's theorem therefore gives:

> **Coefficient-space corollary (Jelonek).** Every irreducible component of
> `X(3,d)` containing `F_N` is generically noninvertible and has dimension at
> least
> \[
> 3^2-1=8.
> \]

More precisely, the nonautomorphisms form a dense open subset of every such
component.  The corollary supplies neither the number of components through
`F_N` nor their dimensions beyond the lower bound.  It also does not imply
that `F_N` is a smooth or reduced point.

## 3. The ordinary tangent space

The incremental coefficient-space project starts with the full scheme (1),
not with a weighted-support ansatz.  Put `F=F_N` and write a first-order
deformation as

\[
 F_\epsilon=F+\epsilon G,\qquad
 G\in{\cal V}(3,d),\qquad \epsilon^2=0.
\]

Linearizing the determinant gives

\[
 \det D(F+\epsilon G)
 =1+\epsilon\,\operatorname{tr}
   \bigl(\operatorname{adj}(DF)\,DG\bigr).              \tag{3}
\]

Thus the linearized Jacobian operator is

\[
 L_F:{\cal V}(3,d)\longrightarrow
 \mathbb C[x,y,z]_{\leq3(d-1)},\qquad
 G\longmapsto
 \operatorname{tr}\bigl(\operatorname{adj}(DF)\,DG\bigr), \tag{4}
\]

and

\[
 \boxed{T_FX(3,d)=\ker L_F,\qquad
 \dim T_FX(3,d)=3\binom{d+3}{3}-\operatorname{rank}L_F.} \tag{5}
\]

Equation (4), expanded monomial by monomial, is the requested system of
linearized Jacobian equations.  Since `det DF=1`, one may equivalently replace
`adj(DF)` by `(DF)^{-1}`; the adjugate form is preferable for exact sparse
coefficient arithmetic.

Jelonek's bound already implies

\[
 \dim T_FX(3,d)\geq8,                                  \tag{6}
\]

but equality in (6) would not by itself prove that there is a unique smooth
eight-dimensional component through `F`.

## 4. The unrestricted infinitesimal quotient collapses

Let `U` and `V` be polynomial vector fields on the target and source,
respectively.  The dual-number target and source automorphisms
`id+epsilon U` and `id+epsilon V` induce the first-order left--right
variation

\[
 \delta_{U,V}F=U\circ F+DF\cdot V.                     \tag{7}
\]

It belongs to the fixed-Jacobian tangent space precisely when

\[
 (\operatorname{div}U)\circ F+\operatorname{div}V=0,   \tag{8}
\]

and it belongs to the chosen coefficient box when its degree is at most
`d`.  Define the unrestricted infinitesimal left--right subspace by

\[
 {\cal O}^{LR}_{F,d}
 =
 \left\{
 U\circ F+DF\cdot V:
 \begin{array}{l}
 \deg(U\circ F+DF\cdot V)\leq d,\\
 (\operatorname{div}U)\circ F+\operatorname{div}V=0
 \end{array}
 \right\}
 \subseteq T_FX(3,d),                                  \tag{9}
\]

where `U,V` range over all polynomial vector fields.  No rank computation is
needed: the source part alone exhausts the tangent space.

Indeed, for any `G in T_FX(3,d)`, put

\[
 V=(DF)^{-1}G=\operatorname{adj}(DF)G.                 \tag{10}
\]

This is a polynomial vector field, and

\[
 F\circ(\operatorname{id}+\epsilon V)
 =F+\epsilon\,DF\cdot V
 =F+\epsilon G.                                        \tag{11}
\]

The chain rule, or the linearized determinant identity, gives
`div(V)=L_F(G)=0`.  Thus (11) is a determinant-one source automorphism over
the dual numbers.  Consequently

\[
 \boxed{{\cal O}^{LR}_{F,d}=T_FX(3,d),\qquad
 {\cal N}^{LR}_{F,d}=0.}                               \tag{12}
\]

This is the first-order form of the repository's
[formal orbit-triviality theorem](FORMAL_ORBIT_TRIVIALITY.md).  In
particular, the proposed unrestricted quotient

\[
 {\cal N}^{LR}_{F,d}
 =T_FX(3,d)/{\cal O}^{LR}_{F,d}
\]

cannot detect any of the stable parameters.

## 5. What ordinary coefficient deformation theory can detect

The collapse (12) sharply separates two questions.

The raw tangent space (5) remains useful.  Together with local equations it
can show that `F_N` is singular or nonreduced in `X(3,d)`, and primary
decomposition of the completed local ring can show that several components
pass through it.  Jelonek then says that every one of those components is
generically counterexample-producing and has dimension at least eight.

By contrast, the normalized weighted seed space has dimension `N-3`, and
decorated normalization detects an `(N-3)`-dimensional family of stable
left--right classes on a nonempty clean open.  Differentiating the family at
the chosen seed gives

\[
 T_{H_N}{\cal S}_N\longrightarrow T_{F_N}X(3,d),
 \qquad \dim T_{H_N}{\cal S}_N=N-3.                    \tag{13}
\]

Composing (13) with the unrestricted LR quotient is identically zero by
(12).  This is not a contradiction.  The dual-number source automorphism in
(11) need not algebraize to one polynomial automorphism family that
trivializes a positive-dimensional seed family, and its required complexity
can grow without a uniform bound.  Stable moduli live in precisely that
formal-to-algebraic and bounded-complexity gap.

Thus the meaningful comparisons are:

- the rank of (13) inside the raw coefficient tangent space;
- `dim T_{F_N}X(3,d)` versus the dimensions and tangent spaces of the local
  irreducible components;
- quadratic and higher local obstructions to raw coefficient directions;
- bounded-degree or Rees-filtered contact spaces before the cutoff is
  allowed to grow, followed by the algebraization/no-escape problem for an
  entire parameter family.

For a vector-field cutoff `b`, one may define

\[
 {\cal O}^{LR,\leq b}_{F,d}\subseteq T_FX(3,d)          \tag{14}
\]

by imposing `deg U,deg V<=b` in (9).  These spaces measure the complexity
needed to kill a tangent direction, but they are not stable-moduli tangent
spaces.  In fact (10) gives a finite cutoff that kills every fixed
first-order direction.  Detecting stable moduli therefore requires uniform
higher-order or family-level algebraization, not merely a larger
first-order matrix.

For the foundational case `N=3`, the stable-parameter count is zero.  The
existing
[sixteen-monomial coefficient calculation](FOUNDATIONAL_WEIGHTED_COEFFICIENT_SCHEME.md)
finds a nonreduced tangent direction in a much smaller normalized ansatz and
then proves that direction is polynomial source-orbit tangent.  It does not
compute (5) for the full bounded-degree scheme.

## 6. Exact full-box tangent calculation in degrees four through six

Write a general point of the minimal full coefficient box as

\[
 \Phi_i=\sum_{|\alpha|\leq d}c_{i,\alpha}x^\alpha
 \qquad(1\leq i\leq3).
\]

The bounded-degree determinant-one coefficient scheme used here is the
explicit affine scheme

\[
 X(3,d)=\operatorname{Spec}
 \frac{\mathbb Q[c_{i,\alpha}]}
 {\left([x^\beta](\det D\Phi-1):
          |\beta|\leq3(d-1)\right)}.                 \tag{15}
\]

Zero coefficient equations are harmless in (15).  There are
`3 binomial(d+3,3)` coefficient variables and `binomial(3d,3)` displayed
coefficient slots.

For the integer-root maps in
[the all-degree rational-fiber theorem](../verified/ALL_DEGREE_RATIONAL_FIBERS.md),
exact sparse elimination over `Q` gives:

| `N` | coordinate degrees of `F_N` | `d` | variables | determinant slots | nonzero rows of `L_F` | `rank_Q L_F` | `dim T_F X(3,d)` |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | `(12,11,4)` | 12 | 1365 | 7140 | 3135 | 1307 | 58 |
| 5 | `(17,16,4)` | 17 | 3420 | 20825 | 9065 | 3332 | 88 |
| 6 | `(22,21,4)` | 22 | 6900 | 45760 | 19808 | 6777 | 123 |

These are ranks over `Q`, not finite-field estimates.  The sparse echelon
calculation uses the graded-lexicographic leading monomial of each column and
retains exact rational pivots.  The largest numerator or denominator bit
sizes in the normalized pivots are respectively `33`, `101`, and `204`.

The tangent space of the normalized seed family has the concrete basis

\[
 h_j(W)=W^{j+2}(W-1)^2,\qquad 0\leq j\leq N-4.       \tag{16}
\]

Differentiating the weighted suspension in these directions gives rank
exactly `N-3` inside the raw tangent spaces in the table.  For each resulting
coefficient tangent `G_j`, the checker also constructs

\[
 V_j=\operatorname{adj}(DF_N)G_j,\qquad
 DF_N\,V_j=G_j,\qquad \operatorname{div}V_j=0.       \tag{17}
\]

The maximum coordinate degrees of these canonical source trivializers are:

| `N` | seed-tangent rank | maximum degrees of `V_0,V_1,...` |
|---:|---:|---:|
| 4 | 1 | `25` |
| 5 | 2 | `30,35` |
| 6 | 3 | `35,40,45` |

Thus the visible stable parameters are honest reduced family directions in
the full coefficient scheme and are linearly independent before quotienting.
Their canonical first-order source gauges already leave the original
coefficient-degree boxes.  This is a concrete filtered-contact datum, not a
stable tangent invariant: a target gauge may lower a particular filtered
representative, and higher orders require the full two-sided resource
spectrum.

The raw bounded schemes already have nonzero quadratic Kuranishi maps.  Put

\[
 \gamma_N=1+a_Nxy+x^2z,\qquad
 (a_4,a_5,a_6)=\left(-\frac53,-\frac{11}7,-\frac{57}{34}\right),
\]

and define

\[
 P_{N,m}=x^m\gamma_N^{m-1},\qquad
 G_{N,m}=(yP_{N,m},P_{N,m},0).
\]

For `(N,m)=(4,3),(5,3),(6,4)`, exact differentiation gives

\[
 L_{F_N}(G_{N,m})=0,\qquad
 Q_{F_N}(G_{N,m})
 =-x^{2m+2}\gamma_N^{2m-2}.                         \tag{18}
\]

Here `Q_F(G)` is the coefficient of `t^2` in
`det D(F+tG)`.  Exact reduction by the rational echelon basis of
`im(L_F)` gives nonzero canonical remainders with respectively `10`, `3`,
and `9` terms.  Thus (18) is nonzero in `coker(L_F)`, so no correction `K`
of degree at most `d` can make

\[
 F_N+tG_{N,m}+t^2K
\]

Keller modulo `t^3`.  Thus the first nonzero raw bounded-box Kuranishi map is
quadratic for each of `N=4,5,6`, and all three coefficient-scheme points fail
formal smoothness.  This does not distinguish a reduced singularity from
nilpotent structure.

### 6.1 The complete quartic quadratic map

For `N=4`, relation-tracked exact elimination gives a 58-vector basis of
`ker(L_F)`.  Reducing all

\[
 \binom{58+1}{2}=1711
\]

polarized quadratic determinant classes by a **full** canonical echelon
normal form in `coker(L_F)` gives

\[
 \boxed{\operatorname{rank}_{\mathbb Q}K_{2,F_4}=53.} \tag{19}
\]

Exactly `727` basis pairs are nonzero before taking their span, and their
remainders use `338` canonical cokernel monomials.  The largest rational
pivot numerator or denominator in the final quadratic image has `214` bits.

The determinant-preserving affine source--target parameter space has
dimension `23`, but its orbit at `F_4` has dimension `22`: the missing
dimension is the weighted `G_m` stabilizer.  Four further exact target-shear
directions are

\[
 (C^2,0,0),\quad(C^3,0,0),\quad
 (0,C^2,0),\quad(0,C^3,0),                           \tag{20}
\]

where `C=F_{4,3}`.  They raise the reduced orbit-family tangent rank from
`22` to `26`.  Adding the normalized quartic seed direction raises it to
`27`.  Hence there is an explicit reduced 27-dimensional family through
`F_4`; this is a subvariety of the coefficient scheme, not yet a proof that
its closure is an irreducible component.

Modulo each of `32003`, `32009`, and `32027`, quotienting only the affine
orbit gives a 36-variable normal slice.  The seed is one distinguished
coordinate, has zero quadratic pairing with every normal coordinate, and
the remaining quadratic ideal again has rank `53`.  Thirteen other displayed
coordinate axes survive the quadratic equations.  A complete order-three
lookahead kills eight of those thirteen axes and allows five to lift through
order three.  Two of the five are the first-coordinate target shears in
(20).  These are modular discovery calculations.  Beyond order three, a
greedy echelon lift is not an obstruction certificate because kernel choices
made at earlier orders can change later equations; the seed family itself
is the control example.

The exact characteristic-zero checker is
[`verify_quartic_full_box_kuranishi.py`](../scripts/verify_quartic_full_box_kuranishi.py).
The modular slice compiler and bounded-jet research screen is
[`research_quartic_coefficient_kuranishi.py`](../scripts/research_quartic_coefficient_kuranishi.py).

Most importantly, there is **no first nonzero Kuranishi obstruction after
the unrestricted source--target jet quotient**.  The formal
source-triviality theorem applies at every Artin order and makes that quotient
the one-point functor; allowing target jets cannot restore an obstruction.
Consequently the proposed order-by-order quotient and the request for a
later nonzero unfiltered Kuranishi class are incompatible.  The `N-3` stable
parameters first appear only at the reduced algebraization/global-boundary
level, or in a specified bounded-degree/Rees filtration.  They do not appear
at any finite unfiltered Kuranishi order.

The raw bounded scheme (15) has the quadratic classes (18) precisely because
the correcting source automorphisms may leave its degree box.  The full
quartic quadratic map is now known by (19), but determining the radical and
associated primes of its 36-variable affine-normal slice, its higher
restrictions, and the reduced and embedded completed local components
remains open.  Neither the large tangent dimensions nor (18)--(19) proves
nonreducedness or a component count.

Run the exact audit with

```bash
.venv/bin/python scripts/verify_all_degree_coefficient_tangents.py
```

## 7. Revised exact computational deliverables

For each selected `F_N`, beginning with the foundational map and the minimal
bound `d=deg(F_N)`:

1. retain a certified sparse basis of `T_{F_N}X(3,d)`, beyond the exact ranks
   now known for `N=4,5,6` (complete for `N=4`);
2. compute the full raw bounded-box quadratic Kuranishi map for `N=5,6`,
   extending the quartic calculation (19);
3. compute minimal primes and embedded associated primes of the completed
   local coefficient ideal where feasible;
4. for the separate complexity question, compute the filtered spaces
   (14) and higher-order contact growth, then test whether the formal
   trivializers algebraize uniformly along the seed family.

Sparse ranks should first be audited modulo several good primes and then
certified over `Q`.  The output should distinguish ambient degree `d`,
geometric degree `N`, raw tangent dimension, local component dimensions,
filtered orbit dimensions, and higher-order obstruction ranks.  The
unrestricted first-order LR-normal dimension is always zero and should not
be reported as a candidate stable-moduli dimension.
