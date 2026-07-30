# Local coefficient components versus decorated stable moduli

This note records Jelonek's coefficient-space corollary and organizes the
coefficient, filtered-contact, and boundary-decoration calculations into the
single programme `OP-CCDM`.

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

### 6.2 Reduced families and their generic tangent excess

The four shears in (20) are the first case of a uniform reduced family.
Put

\[
 d_N=5N-8,\qquad k_N=\left\lfloor\frac{d_N}{4}\right\rfloor .
\]

Starting from the normalized seed family, add `C^k` independently to the
first and second outputs for every `2<=k<=k_N`, and then take the affine
left--right orbit.  Its tangent rank is

\[
 r_N=22+2(k_N-1)+(N-3).                                \tag{21}
\]

The three summands are respectively the 22-dimensional affine orbit, the
target shears, and the visible boundary parameters.  At selected rational
points of these families, exact characteristic-zero elimination gives:

| `N` | `d_N` | `k_N` | reduced-family rank `r_N` | full generic tangent dimension | transverse tangent excess |
|---:|---:|---:|---:|---:|---:|
| 4 | 12 | 3 | 27 | 49 | 22 |
| 5 | 17 | 4 | 30 | 80 | 50 |
| 6 | 22 | 5 | 33 | 109 | 76 |

Thus the special integer-root points have tangent dimensions `58,88,123`,
while these explicit rational points of the reduced families have the smaller
dimensions `49,80,109`.  Neither row proves that the displayed reduced family
is a component: the transverse tangent excess may contain nilpotent
thickening, additional reduced branches, or both.  The exact checker is
[`verify_generic_coefficient_family_tangents.py`](../scripts/verify_generic_coefficient_family_tangents.py).

For the quartic generic point, quotienting the 27 displayed reduced
directions leaves a 22-dimensional normal tangent space.  The normal
quadratic Kuranishi image has rank `22` modulo each of `32003`, `32009`, and
`32027`; at `32003`, its 22 independent quadrics use 76 cokernel monomials.
Five coordinate axes survive quadratically.  Three fail the complete cubic
axis test, while two lift through order three.  A noncanonical greedy choice
of corrections makes the latter two fail at order four, but that is not an
order-four obstruction certificate because other lower-order kernel choices
can change the result.  The reproducible modular compiler is
[`research_quartic_generic_component.py`](../scripts/research_quartic_generic_component.py).
With the optional full cubic compiler, `2021` of the `2024` normal cubic
monomial classes are nonzero and their cokernel coefficients span `305`
independent cubic equations modulo `32003`.  The five pure-cube coefficients
agree term-for-term with the separate coordinate-axis implementation.
The resulting 327-generator quadratic-plus-cubic homogeneous-layer input is
[`quartic_generic_component_order3_mod32003.sing`](../artifacts/generated-results/quartic_generic_component_order3_mod32003.sing).
It is not the completed local ideal: lower-order tangent choices couple
successive arc equations.  Both the quadratic-only primary decomposition and
this stronger homogeneous-layer standard-basis run exhausted practical memory
before producing a basis (the latter at a 12 GB cap).  Hence no radical,
minimal-prime, or embedded prime is claimed.

### 6.3 Filtered source quotients

The source-only trivializer `V=adj(DF_N)G` is unique.  This makes its degree
an unambiguous first-order filtration, even though it is not the optimal
two-sided left--right filtration.  Modulo `32003`, the full tangent spaces are
exhausted at source cutoffs

\[
 b_{\rm all}(4)=33,\qquad b_{\rm all}(5)=48,\qquad
 b_{\rm all}(6)=63.                                    \tag{22}
\]

The visible seed subspaces disappear one direction at a time at:

| `N` | seed quotient dimension at `b=0` | source cutoffs at which its directions become trivial |
|---:|---:|---:|
| 4 | 1 | `25` |
| 5 | 2 | `30,35` |
| 6 | 3 | `35,40,45` |

The full breakpoint profiles are recorded in
[`filtered_source_tangent_profiles_mod32003.json`](../artifacts/generated-results/filtered_source_tangent_profiles_mod32003.json)
and generated by
[`research_filtered_source_tangent_profile.py`](../scripts/research_filtered_source_tangent_profile.py).
These are good-prime filtered computations, not characteristic-zero rank
certificates.  They show concretely that every visible direction eventually
vanishes in a finite first-order source quotient.  A stable parameter
therefore cannot be defined as a direction that survives every fixed
first-order degree cutoff.

### 6.4 All-degree location of the stable parameters

Combining formal source triviality with the decorated-normalization theorem
gives the requested all-degree statement, but its conclusion is a separation
of levels rather than a late finite-order obstruction:

> **Formal-versus-algebraic location theorem.** For every `N>=4`, the
> fixed-Jacobian deformation functor of `F_N` modulo unrestricted polynomial
> source jets is the one-point functor on local Artin algebras.  Hence every
> finite unfiltered Kuranishi obstruction group and every finite-order
> unfiltered contact invariant vanishes.  On a nonempty reduced
> boundary-clean seed open, however, decorated normalization has image
> dimension `N-3` and is invariant under stable polynomial left--right
> equivalence.  Thus the `N-3` stable parameters first appear at the exact
> reduced/global-boundary level, not at any infinitesimal or finite Artin
> order.

The first assertion is
[formal source triviality](FORMAL_ORBIT_TRIVIALITY.md#2-formal-source-triviality-theorem);
the second is the
[degreewise decorated-normalization theorem](DECORATED_NORMALIZATION_INVARIANT.md).
What remains open is a sharper **filtered algebraization theorem**: determine
the optimal two-sided complexity growth needed to approximate an entire seed
family and prove a uniform no-algebraization statement from that growth.

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

## 7. Unified programme: local component versus decorated moduli (`OP-CCDM`)

There is not yet a distinguished object called *the* nonproper component.
Let

\[
 \widehat{\mathcal O}_{N,d}
 =\widehat{\mathcal O}_{X(3,d),F_N}
\]

and let `p_{N,d,i}` run through its minimal primes.  Jelonek's theorem says
that every corresponding reduced branch is generically nonproper, but it
does not select one branch or identify their scheme structures.  Thus
“the tangent space to the nonproper component” must mean the tangent spaces
and tangent cones of **each** `p_{N,d,i}`, together with their intersections
and embedded structure.

The five calculations now form one comparison:

| layer | present status | unified deliverable |
|---|---|---|
| Keller coefficient tangent | Formula (5) holds in all degrees; the exact ranks at `F_4,F_5,F_6` are `58,88,123`. | Retain certified tangent bases and compute the tangent space and cone of every completed reduced branch, not only the raw scheme tangent. |
| Generically nonproper local branches | Jelonek makes every reduced branch through `F_N` generically nonproper; the quartic raw quadratic map has rank `53`, but no minimal prime is known. | Determine the minimal and embedded primes, branch dimensions, intersections, and bounded-box obstruction maps; extend the full quadratic calculation to `N=5,6`. |
| Formal source trivialization | For every Artin deformation the source trivializer exists uniquely; on tangent vectors it is `adj(DF_N)G`. | Measure its two-sided degree, parameter-height, pole-order, and stabilization growth along reduced families rather than forming the zero unrestricted tangent quotient. |
| First non-algebraizable filtered class | The source-degree profiles (22) are known modulo one good prime. On the `F_2` ordinary-degree face, linear Rees strictness fails and the `p=2` quadratic symbol contributes a new length-one class. | Find the first intrinsic two-sided filtered class that survives all lower gauges and prove that it cannot algebraize uniformly along the marked-root family. |
| Coefficient-to-decoration comparison | Decorated normalization is defined on the reduced boundary-clean seed locus and its image there has dimension `N-3`. | Construct the family-level decoration map on every boundary-clean reduced coefficient branch, compute its rank/kernel on reduced arcs or tangent cones, and compare it with the first surviving filtered class. |

The last row is deliberately a **reduced family-level** comparison.  The
current decorated-normalization theorem is not asserted as a natural
transformation on arbitrary nonreduced Artin coefficient schemes.  Treating
it as one and differentiating would conflict with formal source triviality.
The required construction must record the normalization, Fitting divisor,
conductor, node pairing, and intrinsic boundary marks in families before any
tangent or associated-graded map is claimed.

The foundational nilpotent gluing formerly tracked as `OP-FC` is the `N=3`
calibration of the second row.  The linear and quadratic faces formerly
tracked as `OP-LR-REES` and `OP-LR-II` are calibrations of the fourth row.
They are no longer separate programme endpoints.  The adjacent
`OP-LR-NE` remains the valuative problem of extending an LR equivalence
through a marked-cover compactification; it can supply a no-escape theorem,
but it does not replace the local coefficient calculation.

For each selected `F_N`, use the minimal coefficient bound first.  Sparse
ranks should be audited modulo several good primes and then certified over
`Q`.  Every output must distinguish ambient degree `d`, geometric degree
`N`, raw scheme tangent dimension, reduced-branch tangent and component
dimensions, filtered orbit dimensions, decoration rank, and higher
obstruction ranks.  The unrestricted LR-normal dimension is always zero and
must not be reported as a candidate stable-moduli dimension.
