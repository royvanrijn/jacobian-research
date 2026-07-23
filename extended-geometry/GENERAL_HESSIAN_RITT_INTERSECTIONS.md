# General Hessian--Ritt intersections: relation graphs and the transfer problem

Work over an algebraically closed field of characteristic zero.  This note
extracts the general theorem suggested by the degree-six, degree-eight,
degree-twelve, and degree-eighteen computations, and records the decisive
degree-twenty-four and degree-thirty tests.

The main conclusion has two logically separate parts.

1. For ordinary polynomial composition, the desired multi-intersection
   theorem is essentially already the tame multi-collision theorem of
   [Ziegler](https://arxiv.org/abs/1402.5945), built on the description of
   complete decompositions by
   [Mueller--Zieve](https://arxiv.org/abs/0807.3578).  Its combinatorial object
   is a normalized **relation graph**, not a single refinement tree.
2. For the projective Hessian incidences of this project, one additional
   theorem is needed: projection away from the constant and linear terms must
   not create extra, unsynchronized intersection components.  This is the
   genuinely new Hessian-transfer problem.

Thus the degree-eighteen nested-power threefold is not a counterexample to a
general composition theorem.  It is the exponential block predicted by the
relation-graph normal form.  What it disproves is only the narrower attempt to
index higher intersections by one common refinement word.

## 1. The synchronized polynomial problem

Let

\[
 \mathcal P_N=\{f=W^N+c_{N-1}W^{N-1}+\cdots+c_1W:f(0)=0\}
\]

be the space of monic original degree-`N` polynomials.  For an ordered
factorization

\[
 \mathbf d=(d_1,\ldots,d_r),\qquad \prod_i d_i=N,
\]

let `D_d` be the image of the normalized composition map

\[
 \mathcal P_{d_1}\times\cdots\times\mathcal P_{d_r}
 \longrightarrow \mathcal P_N.
\]

The map is injective on its image in characteristic zero after the monic
original normalization.  Hence

\[
 \dim \mathcal D_{\mathbf d}=\sum_i(d_i-1).       \tag{1.1}
\]

For a set `D` of ordered factorizations, put

\[
 \mathcal D_{N,D}=\bigcap_{\mathbf d\in D}\mathcal D_{\mathbf d}.
                                                               \tag{1.2}
\]

This is the synchronized problem: every requested decomposition belongs to
the same polynomial `f`, including its linear coefficient.

## 2. Gcd refinement and the relation graph

Two ordered factorizations can first be replaced, without changing (1.2), by
their gcd refinements.  Iterating the operation gives a normalized collection
`D*` with the following properties.

* Every word in `D*` is a permutation of the same factor multiset `B`.
* Factors with non-coprime degrees retain their relative order.
* Only adjacent coprime factors can change order.

This is the multiword form of Engstroem reduction.  It can also be seen from
the lattice of intermediate fields between `k(f)` and `k(W)`: meet and join
have the gcd and lcm degrees forced by the totally ramified place at infinity.

The **relation graph** `G_D` has the occurrences in the multiset `B` as its
vertices.  For every normalized word, orient every pair in the order in which
it occurs in that word, and take the union over all words.  Thus:

* a one-way edge records an order that never changes;
* a two-way edge records a required Ritt interchange;
* strongly connected components record groups of factors coupled by several
  compatible Ritt moves;
* the condensation of the strongly connected components is acyclic and gives
  the unique outer-to-inner composition order of the blocks.

This graph is the precise replacement for the proposed common-refinement
tree.  A tree or a single word sees associativity, but it cannot see two
compatible power moves as in degree eighteen or a braid as in degree thirty.

## 3. The ordinary composition-intersection theorem

The characteristic-zero specialization of Ziegler's tame multi-collision
normal form gives the following result.

> **Relation-graph composition theorem.**  Normalize the requested ordered
> factorizations and form their relation graph `G`.  The intersection
> `D_(N,D)` depends only on `G`.  If
> `G_1 -> ... -> G_s` are its strongly connected components in condensation
> order, then every polynomial in the intersection has a unique composition
> `f=f_1 o ... o f_s` with `f_j` belonging to the collision family of `G_j`.
> A nontrivial strongly connected block has exactly two normal-form types:
>
> 1. an **exponential** or power type, assembled from polynomials
>    `W^k R(W)^e` prescribed by the two-way neighborhoods; and
> 2. a **trigonometric** type
>    `D_m(W+t,z)-D_m(t,z)`, where `D_m` is the monic Dickson polynomial.
>
> The two types may overlap or one may lie in the closure of the other.
> Irreducible components are obtained only after normal-form families with the
> same dense image, or with contained image, are identified.

For a two-way constraint between coprime degrees `d>e`, write
`d=se+k`, `1<=k<e`.  Its exponential core is

\[
 \left(W^{ke}R(W^e)^e\right)^{[t]},qquad \deg R=s,             \tag{3.1}
\]

where `[t]` denotes original source translation.  The trigonometric core is

\[
 D_{de}(W+t,z)-D_{de}(t,z).                                    \tag{3.2}
\]

These are precisely the power and Chebyshev cases of Ritt's second theorem.
When `e=2`, the two-factor Chebyshev family is contained in the power normal
form.  In larger diagrams the containment question must be decided for the
whole block, not edge by edge.

This theorem answers four of the motivating questions for the ordinary
composition problem.

* Literal nested components are the blocks with no directed cycle.
* Multiple refinement paths define the same component exactly when their
  normal forms have the same dense image.
* Extra components arise only from exponential and trigonometric strongly
  connected blocks.
* Braid compatibility is a block condition, not an independent condition for
  every edge.

## 4. Pair intersections and their dimensions

The pair case has a particularly compact closed form.  Request outer degrees
`a` and `c`, and set

\[
 g=\gcd(a,c),\qquad
 \ell={a\over g},\qquad m={c\over g},\qquad
 r={N\over\operatorname{lcm}(a,c)}.                            \tag{4.1}
\]

The residual degrees `ell,m` are coprime.

If `a|c` or `c|a`, there is no Ritt core.  The intersection is the literal
common-refinement locus with word

\[
 \left(\min(a,c),{\max(a,c)\over\min(a,c)},{N\over\max(a,c)}\right),
                                                               \tag{4.2}
\]

after deleting entries equal to one.  Its marked Hessian dimension is

\[
 \sum_i(d_i-1)-1.                                             \tag{4.3}
\]

If neither degree divides the other, strip the common outer degree `g` and
common inner degree `r`.  With `L=min(ell,m)`, `M=max(ell,m)`, and
`s=floor(M/L)`, the two synchronized normal-form families have marked
dimensions

\[
 \boxed{\dim_{\rm pow}=g+r+s-2,\qquad
        \dim_{\rm Dickson}=g+r-1.}                            \tag{4.4}
\]

Formula (4.4) counts arbitrary common outer and inner decorations, the core
normal-form parameters, and then removes the one-dimensional marked source
dilation.  It recovers every degree-six and degree-twelve pair dimension in
the existing atlas.

## 5. The general dimension rule

For a vertex `v` in a nontrivial strongly connected block, let `e_v` be the
product of the degrees of all vertices joined to `v` by a two-way edge.  The
undirected two-way subgraph is connected, so `e_v>1`.  Define

\[
 \epsilon(d,e)=
 \begin{cases}
 d-1,&e=1,\\
 \lfloor d/e\rfloor,&e>1,
 \end{cases}                                                  \tag{5.1}
\]

where `e` is this two-way-neighborhood degree.  An
exponential block has affine monic-original dimension

\[
 1+\sum_{v\in G}\epsilon(d_v,e_v),                            \tag{5.2}
\]

the first parameter being the common original shift.  A trigonometric block
has affine dimension two, from `(t,z)`.  A singleton block of degree `d` has
dimension `d-1`.

For a choice of exponential or trigonometric normal form on every strongly
connected block, add the block dimensions.  The projective marked-Hessian
dimension is

\[
 \boxed{\dim_{\rm marked}=\sum_j\dim_{\rm affine}(G_j)-1.}     \tag{5.3}
\]

Finally discard choices whose image is contained in another choice.  This is
the desired dimension formula from the refinement poset, with one important
correction: the poset must carry its two-way Ritt graph and the power versus
Dickson labels.  The bare divisor poset does not determine (5.2).

## 6. The Hessian-transfer conjecture

Let

\[
 \pi_N:\mathcal P_N\longrightarrow\mathcal K_N
\]

forget the constant and linear terms, equivalently retain the monic second
derivative.  The Hessian-composition locus is the projectivized image

\[
 \mathcal C_{a,b}=\overline{\pi_N(\mathcal D_{(a,b)})}.
\]

For several requested cuts there is always an inclusion

\[
 \overline{\pi_N\left(\bigcap_i\mathcal D_i\right)}
 \subseteq \bigcap_i\overline{\pi_N(\mathcal D_i)}.           \tag{6.1}
\]

The reverse inclusion is not formal: image does not commute with
intersection.  A point on the right can a priori have decompositions

\[
 F+\lambda_iW=A_i\circ B_i
\]

with different linear coefficients `lambda_i`.  Ritt's theorems apply only
after these vertical parameters synchronize.

The correct new conjecture is therefore:

> **Hessian synchronization and transfer conjecture.**  On the exact-degree
> leading-coefficient chart in characteristic zero, every irreducible
> component of a multiple Hessian-composition intersection is generically
> synchronized.  Equivalently, (6.1) is equality on reduced supports after
> projectivization.  The normalization of every component is consequently a
> relation-graph power/Dickson family, and its dimension is (5.3).

This formulation is stronger and more precise than the earlier
common-refinement-tree conjecture.  It also isolates the only step not already
covered by the tame multi-collision literature.

The canonical atlas makes synchronization testable.  For each requested word
it reconstructs a unique linear coefficient `lambda_d` from the factor
parameters.  The conjecture says that on every minimal prime of the Hessian
intersection,

\[
 \lambda_{\mathbf d}=\lambda_{\mathbf e}                     \tag{6.2}
\]

for all requested words.  Thus (6.2), rather than another full coefficient
elimination, is the next general symbolic target.

## 7. The low-degree diagrams in relation-graph form

| degree | normalized basis / graph | dominant all-order family | affine dimension | marked dimension |
|---:|---|---|---:|---:|
| 6 | two-way edge `2--3` | power, containing Dickson | 2 | 1 |
| 8 | literal word `2 o 2 o 2` | common refinement | 3 | 2 |
| 12 | path through basis `2,2,3` | Dickson `D_12` | 2 | 1 |
| 18 | path `3--2--3` with ordered repeated cubics | exponential nested power | 3 | 2 |
| 24 | star with center `3` and three ordered quadratic leaves | Dickson `D_24` | 2 | 1 |
| 30 | full braid on `2,3,5` | Dickson `D_30` | 2 | 1 |

For degree eighteen, (5.2) gives one translation plus one power parameter for
each cubic, exactly

\[
 \bigl((W-s)^3+u(W-s)\bigr)
 \stackrel{A}{\longmapsto}
 \bigl(A(B(W-s))\bigr)^2,\qquad A(Z)=Z^3+vZ.                \tag{7.1}
\]

This explains both the unique affine threefold and why it is larger than the
Dickson curve.

## 8. New degree-twenty-four path certificate

Start on the generic prime-word chart `3 o 2 o 2 o 2`.  It already carries
the cuts `(3,8)`, `(6,4)`, and `(12,2)`.  Pulling back the canonical residuals
for `(2,12)`, `(4,6)`, and `(8,3)` gives one minimal prime of affine dimension
two.  In chart parameters it is

\[
\begin{aligned}
 b_1&=-{b_3^4\over16}+{b_2b_3^2\over4}+{b_2^2\over4},\\
 a_2&=-{3b_2^4\over32}+{3b_1b_2^2\over8}+{3b_1^2\over8},\\
 a_1&=-{3b_1^4\over64}+{a_2b_1^2\over8}+{a_2^2\over4}.
\end{aligned}                                                \tag{8.1}
\]

Set

\[
 t={b_3\over2},\qquad z=-{2b_2-b_3^2\over8}.
\]

Then the polynomial on (8.1) is identically

\[
 \boxed{D_{24}(W+t,z)-D_{24}(t,z).}                           \tag{8.2}
\]

Thus the entire compatible four-vertex path on this chart has one dense
Dickson image.  The checker also fixes `z=1` and proves exactly that the
endpoint polynomial in `t` is squarefree and not exhausted by the derivative
gap, marked-Hessian, primitive-root discriminant, or `H''(1)+2` deletions.
Hence this component meets the marked clean open.

## 9. New degree-thirty braid certificate

Start on the generic prime-word chart `2 o 3 o 5`.  Impose the four remaining
cuts `(3,10)`, `(5,6)`, `(10,3)`, and `(15,2)`.  Exact minimal-prime
computation again gives one prime of affine dimension two.

Let `t,z` be parameters and put

\[
 c_5=D_5(t,z),\qquad c_{15}=D_{15}(t,z).
\]

The prime has the parametrization

\[
\begin{aligned}
 r_4&=5t,& r_3&=10t^2-5z,& r_2&=10t^3-15zt,\\
 r_1&=5t^4-15zt^2+5z^2,&
 q_2&=3c_5,& q_1&=3(c_5^2-z^5),& p_1&=2c_{15}.
\end{aligned}                                                \tag{9.1}
\]

Substitution gives the exact identity

\[
 \boxed{D_{30}(W+t,z)-D_{30}(t,z).}                           \tag{9.2}
\]

Therefore the two reduced paths around the `S_3` braid have the same dense
image; there is no path-labelled second component.  The same exact clean-open
gcd test as in degree twenty-four proves that (9.2) meets the marked clean
open.

## 10. What remains to prove

The experiments now point to a short, theorem-driven program.

1. **Prove synchronization.**  Establish (6.2) generically on every minimal
   prime, preferably by a differential or ramification argument rather than
   degree-by-degree elimination.
2. **Transfer the normal form.**  Once synchronization holds, apply the
   relation-graph theorem blockwise and the finite marked-Hessian
   normalization to obtain the component classification and (5.3).
3. **Control scheme structure.**  Ziegler's theorem classifies reduced
   points/families; nilpotents, embedded primes, and excess conormal modules
   remain separate intersection-theoretic questions.
4. **Globalize the degree-24/30 certificates.**  The new computations are on
   one prime-word chart.  A full Hessian-space calculation, or the general
   synchronization theorem, is needed to exclude components not meeting that
   chart.
5. **Produce clean witnesses blockwise.**  The Dickson blocks now have exact
   algebraic clean-open certificates in degrees 24 and 30.  Exponential
   blocks should admit a uniform specialization argument extending the
   rational degree-eighteen witness.

The best theorem target is therefore not “all components are common
refinement trees.”  It is:

\[
 \boxed{\text{Hessian intersections}
 =\text{ synchronized relation-graph multi-collisions},}
\]

with power/Dickson block normal forms, dense-image identification, and the
additive block dimension formula (5.3).

## Reproduction

```bash
.venv/bin/python scripts/explore_hessian_ritt_component_census.py
.venv/bin/python scripts/explore_degree18_hessian_ritt_power_core.py
.venv/bin/python scripts/explore_degree24_hessian_ritt_all_orders.py
.venv/bin/python scripts/explore_degree30_hessian_ritt_braid.py
```
