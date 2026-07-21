# C08--C11 rebuilt from coincident-root geometry

This note gives one proof chain for the contact strata, their closures,
components, normalizations, and set-theoretic intersections.  Fix an
algebraically closed field `k` of characteristic zero and an integer `n>=3`.

## 1. Coincident-root normalizations and the quotient slice

For a partition

\[
 \lambda=(\lambda_1,\ldots,\lambda_\ell),\qquad
 \lambda_i\ge2,\qquad \sum_i\lambda_i=n,
\]

let `Delta_lambda` be the monic coincident-root locus and let

\[
 M_\lambda(W)=\prod_{i=1}^{\ell}(W-r_i)^{\lambda_i}.               \tag{1}
\]

Permutations of roots carrying the same multiplicity form
`G_lambda=product_mu S_(a_mu)`.  Put

\[
 X_\lambda=\mathbb A^\ell/G_\lambda.
\]

Equivalently, if `Q_mu` is the monic polynomial whose roots are the roots
carrying multiplicity `mu`, then `X_lambda` has polynomial coordinates given
by the coefficients of the `Q_mu`, and

\[
 (Q_2,Q_3,\ldots)\longmapsto\prod_{\mu\ge2}Q_\mu^\mu.              \tag{2}
\]

The exact-root locus `X_lambda^ex` is obtained by removing the discriminants
of the `Q_mu` and the pairwise resultants `Res(Q_mu,Q_nu)`.

Define retained-coefficient functions

\[
 \Phi(M)=M(1)-M(0)-M'(0),\quad
 D(M)=M'(1)-M'(0),\quad
 A(M)=M''(1)-2D(M).                                                \tag{3}
\]

The root model for the exceptional stratum is

\[
 Y_\lambda^{\rm adm}
 =X_\lambda^{\rm ex}\cap V(\Phi)\cap D(D A).                      \tag{4}
\]

Its map to normalized seeds is

\[
 \rho(M)={-M+M'(0)W+M(0)\over D(M)}.                              \tag{5}
\]

The image is the exact full-contact stratum `E_lambda`.

## 2. Finiteness of the coincident-root map

Put

\[
 p_j=\sum_{i=1}^{\ell}\lambda_i r_i^j,\qquad1\le j\le\ell.       \tag{6}
\]

The forms `p_1,...,p_ell` have no common projective zero.  Indeed, group the
nonzero `r_i` by distinct values `u_1,...,u_e`, and let `Lambda_a` be the sum
of the positive integer weights in group `a`.  The first `e` equations are

\[
 \sum_{a=1}^e\Lambda_a u_a^j=0,\qquad1\le j\le e.
\]

The determinant of `(u_a^j)` is a nonzero product of the `u_a` and a
Vandermonde determinant.  Hence every `Lambda_a=0`, impossible in
characteristic zero.  Thus all `r_i=0`.

The homogeneous base-locus criterion now shows that

\[
 (r_i)\longmapsto(p_1,\ldots,p_\ell)                               \tag{7}
\]

is finite.  Weighted Newton identities are triangular and identify the
`p_j` with the highest `ell` nonleading coefficients
`c_1,...,c_ell` of `M_lambda`.  Passing to the finite equal-weight quotient
therefore gives a finite map

\[
 X_\lambda\longrightarrow\mathbb A^\ell_{c_1,\ldots,c_\ell}.      \tag{8}
\]

On `X_lambda^ex`, unique factorization recovers every `Q_mu`, so (2) is
generically one-to-one.  The invariant ring defining `X_lambda` is normal;
therefore (2) is the finite birational normalization of `Delta_lambda`.

Since `ell<=n/2`, the coefficient `c_j` has index `n-j>=2`.  It survives
projection modulo `<1,W>`.  On the normalized slice it is recovered from a
seed by

\[
 c_j={h_{n-j}\over h_n},\qquad1\le j\le\ell.                     \tag{9}
\]

Consequently (5), both before and after restriction by (4), is finite onto
its image.  This is the finiteness input used in every closure and
normalization argument below.

## 3. Nonemptiness of every full-contact stratum

It is not enough to observe that `Phi` is a nonzero polynomial; its zero set
must meet the exact and admissible opens.  We prove this uniformly.

First consider the maximally collided polynomial

\[
 M_r(W)=(W-r)^n,qquad
 \phi_n(r)=(1-r)^n-(-r)^n-n(-r)^{n-1}.                             \tag{10}
\]

The polynomial `phi_n` has degree `n-2` and is squarefree.  To see this, put
`y=-r`, `x=1-r=y+1`, and `z=x/y`.  A zero cannot have `y=0`, and (10) becomes

\[
 \sum_{i=0}^{n-1}z^i-n=0.
\]

After removing the extraneous factor at `z=1`, its roots are those of

\[
 g_n(z)={z^n-nz+n-1\over(z-1)^2}.                                 \tag{11}
\]

If `z^n-nz+n-1` and its derivative have a common root, then
`z^(n-1)=1`, and substitution gives `(n-1)(1-z)=0`.  Thus `z=1`; its
multiplicity is exactly two because the second derivative there is
`n(n-1)`.  Hence `g_n`, and therefore `phi_n`, is squarefree.

At a zero of `phi_n`, one always has `D(M_r)!=0`.  Otherwise
`x^(n-1)=y^(n-1)`, and (10) becomes `(1-n)y^(n-1)=0`, which is incompatible
with `x-y=1` and the derivative equality.

At least one zero also satisfies `A(M_r)!=0`.  Both `A(M_r)` and `phi_n`
have degree `n-2`; their leading-coefficient ratio is `-2`.  But at `r=0`,

\[
 \phi_n(0)=1,qquad A(M_0)=n(n-3)\ne-2.
\]

Since `phi_n` is squarefree, `A` cannot vanish at all of its roots.  Choose a
root `r_0` with `D A!=0`.

Now split this root according to an arbitrary `lambda`.  Choose distinct
`c_i` and put

\[
 r_i(\epsilon)=r_0+\epsilon c_i,qquad
 M_\epsilon=\prod_i(W-r_i(\epsilon))^{\lambda_i}.
\]

Let

\[
 F_\epsilon(z)=M_\epsilon(z)-M_\epsilon(0)-zM_\epsilon'(0).
\]

Then `F_0(1)=0` and `partial F_0/partial z(1)=D(M_(r_0))!=0`.  The formal
implicit-function theorem gives a unit `beta_epsilon`, with constant term
one, satisfying `F_epsilon(beta_epsilon)=0`.  The monic rescaling

\[
 \widetilde M_\epsilon(W)
 =\beta_\epsilon^{-n}M_\epsilon(\beta_\epsilon W)                  \tag{12}
\]

has `Phi=0`.  Over `k((epsilon))` its roots are distinct with exact
multiplicities `lambda`, while `D` and `A` remain units.  Thus (4) is
nonempty after a field extension, hence is a nonempty finite-type scheme over
the algebraically closed field `k`.  Therefore:

\[
 \boxed{\mathcal E_\lambda\ne\varnothing
        \quad\text{for every full-contact partition }\lambda.}    \tag{13}
\]

## 4. Dimension and codimension

The irreducible normal variety `X_lambda` has dimension `ell`.  The function
`Phi` is nonzero—at the all-zero root tuple it equals one—and (13) shows that
its zero locus is nonempty.  Every irreducible component of the principal
divisor `V(Phi)` therefore has dimension `ell-1`.  The exact-root and
admissibility conditions remove proper closed subsets and, by (13), leave a
nonempty open.

The finite seed map (5) preserves dimension.  Since the normalized
admissible seed space has dimension `n-3`,

\[
 \boxed{\dim\mathcal E_\lambda=\ell(\lambda)-1},\qquad
 \boxed{\operatorname{codim}_{\mathcal A_n}\mathcal E_\lambda
        =n-\ell(\lambda)-2}.                                      \tag{14}
\]

## 5. Exact closure order

Write `lambda preceq mu` if the parts of `lambda` can be partitioned into
blocks whose sums are the parts of `mu`; equivalently, `mu` is obtained by
colliding roots of type `lambda`.

If `lambda preceq mu`, repeat construction (12), now splitting each root of
an arbitrary admissible polynomial of exact type `mu` according to the
corresponding block.  Since `D!=0`, the same implicit rescaling keeps
`Phi=0`; admissibility remains open.  Hence

\[
 \mathcal E_\mu\subseteq\overline{\mathcal E_\lambda}.             \tag{15}
\]

Conversely, the finite map from the retained-collision model
`Z_lambda=V(Phi) intersect D(DA) subset X_lambda` has closed image in the
normalized admissible seed space, equal to `closure(E_lambda)`.  Indeed, no
component retained by `D(DA)` is trapped
in a root diagonal: restriction of `Phi` to any diagonal is nonzero because
the all-zero root tuple lies on every diagonal and has `Phi=1`.  Thus the
exact-root open is dense in the retained model.  A point over a limiting seed
still gives a genuine
factorization (2), now with some roots collided.  Its exact multiplicities
are sums of blocks of parts of `lambda`.  The monic omitted polynomial is
unique by C07, so if the seed has exact type `mu`, this limiting polynomial
must be that type-`mu` polynomial.  Therefore `lambda preceq mu`.

Thus the order is exact:

\[
 \boxed{
 \mathcal E_\mu\subseteq\overline{\mathcal E_\lambda}
 \quad\Longleftrightarrow\quad \lambda\preceq\mu.}                \tag{16}
\]

Finiteness is essential in the converse: it prevents a root factorization
from escaping while its retained coefficients converge.

## 6. Maximal `2/3` strata are irreducible

The minimal partitions for `preceq` are exactly

\[
 \lambda=2^a3^b,\qquad2a+3b=n,
\]

because two and three are the indecomposable elements of
`{2,3,4,...}`.  Write

\[
 M=Q^2R^3,\qquad \deg Q=a,\quad\deg R=b.
\]

In quotient coordinates set

\[
 x=Q(0),\ u=Q'(0),\ X=Q(1),\qquad
 y=R(0),\ v=R'(0),\ Y=R(1).
\]

Then

\[
 \Phi=X^2Y^3-x^2y^3-2xuy^3-3x^2y^2v.                             \tag{17}
\]

If `b>=3`, the endpoint data `(y,v,Y)` are independent and (17) is
primitive linear in `v`; its coefficient `-3x^2y^2` is coprime to the
constant term, since reduction modulo `x` or `y` retains `X^2Y^3`.
Gauss's lemma proves irreducibility.

If `a>=3` and `b<=2`, the data `(x,u,X)` are independent and (17) is
quadratic in `X`.  Its constant term

\[
 -xy^2(xy+2uy+3xv)
\]

has `x`-adic valuation one, while the leading coefficient `Y^3` has
valuation zero.  The corresponding quotient in the coefficient fraction
field is not a square, so the quadratic is irreducible.

The only remaining pairs are

\[
 (0,1),(0,2),(1,1),(1,2),(2,0),(2,1),(2,2).
\]

After using the degree-zero, one, and two endpoint relations, each is linear
or quadratic with the following certificate:

| `(a,b)` | Variable | Certificate |
|---|---|---|
| `(0,1)` | `y` | `Phi=3y+1` |
| `(0,2)` | `y` | `Disc=3(v+1)^3(3v-1)` |
| `(1,1)` | `x` | `Disc=4y^2(6y^2+8y+3)` |
| `(1,2)` | `x` | `Disc=-4y^2G`, with `G(0,y)=y(3y^2+3y+1)` |
| `(2,0)` | `x` | `Phi=u^2+2u+2x+1` |
| `(2,1)` | `u` | the discriminant has `x`-adic valuation one |
| `(2,2)` | `u` | the discriminant has `x`-adic valuation one |

The odd valuations exclude squares in the relevant fraction fields; in the
`(1,2)` case, specialization at `v=0` gives the displayed odd `y`-valuation.
These seven exact identities are verified independently by
`verify_maximal_phi_irreducibility.py`; they are the complete finite endpoint
list, not a bounded-degree extrapolation.

Hence `V(Phi) subset X_(2^a3^b)` is irreducible.  Its exact admissible open is
nonempty by (13), so `E_(2^a3^b)` is irreducible.  By (16), every other
stratum lies in one of these closures, and no minimal type lies in another.
Thus the irreducible components are exactly

\[
 \boxed{\mathcal C_{a,b}=\overline{\mathcal E_{2^a3^b}}.}          \tag{18}
\]

## 7. The proposed smooth models are the normalizations

Retain root collisions and define

\[
 \widetilde{\mathcal C}_{a,b}
 =V(\Phi)\cap D(D A)\subset\mathbb A^a_Q\times\mathbb A^b_R.      \tag{19}
\]

On this open, `x` and `y` are nonzero.  If, for example, `x=0`, then
`Phi=0` forces `X=0` or `Y=0`; hence `M'(0)=M'(1)=0`, contradicting
`D!=0`.  The argument for `y` is identical.

If `a>=3`, the independent coordinate `u` satisfies

\[
 {\partial\Phi\over\partial u}=-2xy^3\ne0;
\]

if `b>=3`, the independent coordinate `v` satisfies

\[
 {\partial\Phi\over\partial v}=-3x^2y^2\ne0.
\]

The same seven endpoint pairs remain.  For each, the ideal generated by
`Phi` and all first partials saturates to the unit ideal by `D A`; the exact
identities are checked in `verify_component_normalization.py`.  Thus (19) is
smooth in every case.  It is integral by the irreducibility just proved, so
it is normal.

The normalized coefficient map

\[
 \pi_{a,b}:\widetilde{\mathcal C}_{a,b}\longrightarrow\mathcal A_n,
 \qquad(Q,R)\longmapsto\rho(Q^2R^3)                                \tag{20}
\]

is finite by (8)--(9).  Its image is closed, contains the dense exact stratum,
and therefore equals `C_(a,b)`.  Over the exact `2/3` stratum, unique
factorization recovers `Q` and `R`, so (20) is generically one-to-one and
birational.  A finite birational map from a normal integral source is the
normalization.  Therefore

\[
 \boxed{\pi_{a,b}:\widetilde{\mathcal C}_{a,b}
        \longrightarrow\mathcal C_{a,b}
        \text{ is the normalization morphism}.}                   \tag{21}
\]

This proves C11 rather than merely proposing a smooth cover.

## 8. Intersections are common coarsenings

Let `lambda,mu` be arbitrary full-contact partitions.  If a seed belongs to
both `closure(E_lambda)` and `closure(E_mu)`, finiteness of the two retained
root models supplies limiting factorizations of both types.  They cannot
represent distinct omitted pencil values by C07.  Their monic omitted
polynomials are therefore equal.  Its exact multiplicity partition `nu` is
obtained by colliding parts of both `lambda` and `mu`.

Conversely, if `nu` is a common coarsening, (16) puts `E_nu` in both
components.  Hence, as sets inside the admissible seed space,

\[
 \boxed{
 \overline{\mathcal E_\lambda}\cap\overline{\mathcal E_\mu}
 =\bigcup_{\substack{\lambda\preceq\nu\\\mu\preceq\nu}}
   \mathcal E_\nu.}                                                \tag{22}
\]

For maximal `2/3` types these closures are the components `C_lambda` and
`C_mu`, giving the proposed component-intersection formula.  Equation (22)
is deliberately set-theoretic.  The scheme intersection can be nonreduced
and can have nontrivial conductor multiplicity; C12 and C17--C20 study those
stronger local structures.

## 9. Consequences and proof status

For `lambda=2^a3^b`,

\[
 \dim\mathcal C_{a,b}=a+b-1,
 \qquad\operatorname{codim}_{\mathcal A_n}\mathcal C_{a,b}
 =a+2b-2.
\]

The uniform arguments in Sections 2--5 are prose proofs; computations only
handle the seven endpoint-rank irreducibility and smoothness identities.
Low-degree coincident-root eliminations are regressions, not evidence needed
for the all-degree quantifiers.

For classical context, see Chipalkatti,
[*On equations defining coincident root loci*](https://arxiv.org/abs/math/0110224),
and Feher--Nemethi--Rimanyi,
[*Coincident root loci of binary forms*](https://arxiv.org/abs/math/0311312).
The finite normalization, slice nonemptiness, and exact admissible closure
statements used here are proved directly above.
