# Uniform exceptional-seed theorem

## Normalized seed space

Fix an algebraically closed field `k` of characteristic zero and `n>=3`.  The
normalized admissible seed space is the open subset

\[
\mathcal A_n=\left\{H:\deg H=n,\quad
H(0)=H'(0)=H(1)=0,\quad H'(1)=-1,\quad H''(1)\ne-2\right\}.
\]

Writing

\[
H(W)=\sum_{j=3}^{n}h_j(W^j-W^2),
\qquad
h_n={-1-\sum_{j=3}^{n-1}(j-2)h_j\over n-2},
\]

exhibits `A_n` as an open subset of affine `(n-3)`-space.

Let `N_n` be the locus of normalized seeds whose weighted map omits a target
value over the algebraic closure.  The triangular `C=0` chart covers the
entire target plane `C=0`, so omission can occur only on `C!=0`.

## Full-contact strata

For a partition `lambda=(lambda_1,...,lambda_l)` of `n` with every part at
least two, set

\[
M_\lambda(W)=\prod_{i=1}^{\ell}(W-r_i)^{\lambda_i},\qquad
\Phi_\lambda=M_\lambda(1)-M_\lambda(0)-M_\lambda'(0),
\]

and

\[
D_\lambda=M_\lambda'(1)-M_\lambda'(0).
\]

On `Phi_lambda=0`, away from root diagonals, `D_lambda!=0`, and the weighted
forbidden divisor, define

\[
a=-D_\lambda^{-1},\quad
s={M_\lambda'(0)\over D_\lambda},\quad
t=-{M_\lambda(0)\over D_\lambda},\quad
H={-M_\lambda+M_\lambda'(0)W+M_\lambda(0)\over D_\lambda}.
\]

Roots carrying equal multiplicities are quotiented by their finite
permutation group.  Denote the image of this map by `E_lambda`.

## Theorem

For every `n>=3`,

\[
\boxed{\mathcal N_n=
\bigsqcup_{\substack{\lambda\vdash n\\\lambda_i\ge2}}\mathcal E_\lambda.}
\]

The union is disjoint because every seed has at most one omitted inverse-
pencil value.  Mason--Stothers excludes two distinct full-contact polynomials
except for the even-degree all-double equality case; there the difference of
two monic squares factors as `(A-B)(A+B)` and has degree at least two.  See
`UNIQUE_OMITTED_VALUE.md`.

Every nonempty full-contact stratum satisfies

\[
\boxed{\dim\mathcal E_\lambda=\ell(\lambda)-1},
\qquad
\boxed{\operatorname{codim}_{\mathcal A_n}\mathcal E_\lambda
=n-\ell(\lambda)-2}.
\]

### Proof

For `C!=0`, finite source points are reconstructed exactly from simple roots
of `H(W)-sW+t`.  A target is omitted precisely when that polynomial has no
simple root.  Its complete multiplicity partition therefore has all parts at
least two and gives one of the displayed factorizations.  Conversely, every
such full-contact factorization has no reconstructible root.  Together with
the direct `C=0` calculation, this proves the covering.  The unique omitted-
value theorem and unique factorization make it the displayed disjoint union.

The root parameter space has dimension `ell(lambda)`.  The polynomial
`Phi_lambda` is nonzero—for example its specialization with all roots zero is
`1`.  Nonemptiness of its exact admissible zero set requires a separate
argument: [COINCIDENT_ROOT_REBUILD.md](../../extended-geometry/COINCIDENT_ROOT_REBUILD.md#3-nonemptiness-of-every-full-contact-stratum)
first finds an admissible zero on the maximally collided `(n)` locus, then
splits it according to `lambda` while preserving `Phi=0` by the formal
implicit-function theorem.  Thus the retained principal divisor has
dimension `ell(lambda)-1`.  Quotienting equal parts is finite and does not
change dimension.

It remains to prove that projection to seed coefficients does not lower the
dimension.  Put

\[
p_j=\sum_i\lambda_i r_i^j
\]

and write `c_j` for the coefficient of `W^(n-j)` in the monic polynomial
`M_lambda`.  Newton identities are triangular from `(p_1,...,p_ell)` to
`(c_1,...,c_ell)`, and direct differentiation gives

\[
\det\left({\partial(c_1,\ldots,c_\ell)\over
\partial(r_1,\ldots,r_\ell)}\right)
=(-1)^\ell\left(\prod_i\lambda_i\right)
\prod_{i<j}(r_j-r_i).
\]

This weighted Vandermonde is nonzero off the root diagonals.  Moreover,
`ell<=floor(n/2)`, so all these top monic coefficients are recovered from the
normalized seed by

\[
c_j={h_{n-j}\over h_n}\qquad(1\le j\le\ell).
\]

The coefficient map is therefore immersive on the root space and has rank
`ell-1` after restriction to `Phi_lambda=0`.  This proves the dimension
formula.  Subtraction from `dim A_n=n-3` gives the codimension formula.

## Collision partial order and closures

For permitted partitions of the same integer, define

\[
\lambda\preceq\mu
\]

when the parts of `lambda` can be divided into blocks whose sums are the parts
of `mu`.  Equivalently, `mu` is obtained from `lambda` by merging parts.  This
is a partial order: singleton blocks give reflexivity, composing block
partitions gives transitivity, and mutual refinement forces equal lengths and
therefore identical sorted partitions.

### Collision-closure theorem

If `lambda <= mu`, then

\[
\boxed{\mathcal E_\mu\subseteq\overline{\mathcal E_\lambda}}.
\]

Indeed, choose blocks witnessing the relation and replace every root in one
block by the corresponding root of `M_mu`.  Then

\[
M_\lambda\longrightarrow M_\mu,\qquad
\Phi_\lambda\longrightarrow\Phi_\mu.
\]

It remains to ensure that nearby split roots stay on `Phi_lambda=0`, rather
than merely specializing to it.  Choose distinct generic velocities inside
each block, set `r_i(epsilon)=u_j+epsilon*c_i`, and call the resulting
polynomial `M_epsilon`.  Put

\[
F_\epsilon(z)=M_\epsilon(z)-M_\epsilon(0)-zM_\epsilon'(0).
\]

At the collision point, `F_0(1)=0` and

\[
F_0'(1)=M_\mu'(1)-M_\mu'(0)=D_\mu\ne0.
\]

The formal implicit-function theorem therefore gives
`beta_epsilon in k[[epsilon]]`, with constant term one, such that
`F_epsilon(beta_epsilon)=0`.  The monic polynomial

\[
\widetilde M_\epsilon(W)=
\beta_\epsilon^{-n}M_\epsilon(\beta_\epsilon W)
\]

then satisfies `Phi_lambda=0` exactly.  This formal arc is enough for Zariski
closure.  Its normalization converges to the original seed, while `D!=0` and
the weighted admissibility factor remain units because their constant terms
are nonzero.  This proves the inclusion for every admissible point, including
points where the original root hypersurface is singular.

## Multiple omitted values

`multiple_omission_incidence(degree, partitions)` imposes several
factorizations of one normalized seed.  For each pair it returns the cleared
numerators of `s_i-s_j` and `t_i-t_j`.  Their common vanishing describes one
omitted value represented by two colliding factorizations; separate
Rabinowitsch gates force genuinely distinct omitted values.

For two partitions, `two_omission_incidence` also returns the compact
root-space presentation

\[
aM_\lambda-bM_\mu=\alpha W+\beta,
\]

its coefficient equations, collision discriminant, and off-collision
saturation factor.

### Uniform two-omission theorem

Let `lambda!=mu` be any two full-contact partitions of `n`, with all parts at
least two. Away from their collision discriminants, the affine-difference
system is empty:

\[
\boxed{
aM_\lambda-bM_\mu=\alpha W+\beta,\quad ab\ne0
\quad\Longrightarrow\quad\text{no solutions}.}
\]

Both contact polynomials are monic, so the leading coefficient gives `a=b`.
After division, write `P-Q=L` with `L` affine-linear.  If `L=0`, unique
factorization and exact distinct roots force the same multiplicity partition,
contrary to `lambda!=mu`.  Thus `L!=0`.

The polynomials `P,Q,L` are pairwise coprime.  Indeed, a common root of `P`
and `Q` would occur to multiplicity at least two in `P-Q=L`, impossible for a
nonzero affine polynomial; the other pairwise gcds reduce to the same
observation.  Polynomial Mason--Stothers now gives

\[
n\le \deg\operatorname{rad}(PQL)-1
\le\ell(\lambda)+\ell(\mu).
\]

Define the excess above double contact by

\[
\epsilon(\lambda)=\sum_i(\lambda_i-2)=n-2\ell(\lambda).
\]

Consequently

\[
\ell(\lambda)+\ell(\mu)
=n-{\epsilon(\lambda)+\epsilon(\mu)\over2}.
\]

The excess is nonnegative, and its unique zero is the all-double partition.
Since `lambda!=mu`, the total excess is positive, so the right side is at most
`n-1`, contradicting Mason--Stothers.  Thus abc rigidity separates every pair
of distinct full-contact types; restricting to maximal 2/3 types is not
necessary.

The reason the maximal types themselves use only twos and threes is instead
the contact-atom principle: `{2,3,...}` is the additive semigroup generated by
its indecomposable elements two and three.  See
`CONTACT_ATOM_PRINCIPLE.md`.

In degree six, the `(2,2,2)` and `(3,3)` exact strata are disjoint.  Their
closures meet exactly on the `(6)` common-value collision boundary.  Pulling
the quartic equation of `E_(2,2,2)` back to a rational parameter on
`E_(3,3)` gives, up to invertible factors, the square of the two-root
discriminant.

## Degree eight as a theorem test

The two partitions using only parts two and three are

\[
(2,2,2,2),\qquad(3,3,2).
\]

The theorem gives dimensions `3,2` and codimensions `2,3`, respectively.
Both closures contain `E_(6,2)`, and that common boundary contains `E_(8)`.

The possible off-diagonal intersection can also be eliminated exactly.  Write

\[
Q(W)^2-A(W)^3(W-b)^2=L(W),\qquad \deg L\le1,
\]

with `Q` monic quartic and `A=W^2-a_1W+a_2`.  Comparing coefficients from
degrees seven through four determines `Q`.  The next two coefficients become

\[
-{3\over128}(a_1-2b)(a_1^2-4a_2)^2,
\]

and

\[
{1\over512}(a_1^2-4a_2)^2
(11a_1^2-12a_1b+4a_2-24b^2).
\]

Away from the collision divisor `a_1^2-4a_2=0`, these force `a_1=2b` and
then `a_2=b^2`, a contradiction.  Hence there are no isolated exact
off-diagonal intersections in degree eight.

## Component-poset consequence

The minimal elements of the collision order are precisely the partitions
using only twos and threes.  They index the maximal stratum closures: a part
two or three cannot be split into permitted parts, while every `m>=4` can be
split as `2+(m-2)` and iterated.  The closure theorem therefore places every
exceptional stratum in the closure of a maximal 2/3 stratum.  The
two-omission theorem shows that distinct maximal types are disjoint away from
collision boundaries.

## Uniform irreducibility of maximal strata

Let `a` be the number of double roots and `b` the number of triple roots.  In
equal-part quotient coordinates write

\[
M=Q(W)^2R(W)^3,\qquad \deg Q=a,\quad\deg R=b.
\]

Set

\[
x=Q(0),\quad u=Q'(0),\quad X=Q(1),
\qquad
y=R(0),\quad v=R'(0),\quad Y=R(1).
\]

Then

\[
\Phi_{2^a3^b}
=X^2Y^3-x^2y^3-2xuy^3-3x^2y^2v.                 \tag{1}
\]

For a monic polynomial of degree at least three, its values at zero, its
derivative at zero, and its value at one are independent affine-linear
coordinates.  In degrees two, one, and zero the only relations are

\[
X=1+x+u,\qquad (u,X)=(1,1+x),\qquad (x,u,X)=(1,0,1),
\]

and similarly for `(y,v,Y)`.

### Theorem

For every maximal 2/3 partition of inverse degree at least three,

\[
\boxed{\Phi_{2^a3^b}\text{ is irreducible over }k.}
\]

If `b>=3`, equation (1) is primitive linear in the independent coordinate
`v`.  Its coefficient `-3x^2y^2` is coprime to its constant term: reduction
modulo `x` or `y` leaves the nonzero term `X^2Y^3`, with the degree-zero,
one, and two relations for `Q` giving the same conclusion.  Gauss's lemma
proves irreducibility.

If `a>=3` and `b<=2`, equation (1) is primitive quadratic in the independent
coordinate `X`.  Its constant term is

\[
-xy^2(xy+2uy+3xv),
\]

which has `x`-adic valuation exactly one, while the leading coefficient is
`Y^3`.  The corresponding element of the coefficient fraction field is not
a square, so the quadratic is irreducible.

Only seven endpoint-rank cases remain.  They have the following elementary
certificates; square factors in the displayed discriminants are harmless.

| `(a,b)` | Variable | Certificate |
|---|---|---|
| `(0,1)` | `y` | `Phi=3y+1` |
| `(0,2)` | `y` | `Disc=3(v+1)^3(3v-1)` |
| `(1,1)` | `x` | `Disc=4y^2(6y^2+8y+3)` |
| `(1,2)` | `x` | `Disc=-4y^2G`; `G(0,y)=y(3y^2+3y+1)` |
| `(2,0)` | `x` | `Phi=u^2+2u+2x+1` |
| `(2,1)` | `u` | the discriminant has `x`-adic valuation one |
| `(2,2)` | `u` | the discriminant has `x`-adic valuation one |

For `(1,2)`, a square in the coefficient fraction field would, by unique
factorization, be a constant times a polynomial square; specializing `v=0`
would then have even `y`-valuation, contrary to the displayed formula.
The other odd factors and valuations are immediate from the table.  Thus
every remaining quadratic has nonsquare discriminant.  This completes the
proof in all degrees.

### Irreducible-component theorem

The hypersurface `Phi_lambda=0` is irreducible for every maximal 2/3 type;
removing collision and admissibility divisors leaves a dense irreducible open
set.  Its seed-image closure is therefore irreducible.  Every other stratum
lies in one of these closures by the collision theorem.  Finally, Mason's
inequality also excludes containment of one maximal closure in the boundary
of another.

Retaining collisions while imposing `D!=0` and weighted admissibility gives a
smooth quotient-coordinate hypersurface.  Its finite generically degree-one
map to the seed component is the normalization morphism; see
`COMPONENT_NORMALIZATION.md`.

For the last assertion, the map from `ell` marked roots to the top `ell`
monic coefficients is finite, not merely generically finite.  Its leading
homogeneous forms are equivalent by Newton identities to the weighted power
sums

\[
p_j=\sum_i\lambda_i r_i^j,\qquad 1\le j\le\ell.
\]

They have no common projective zero: after grouping equal nonzero root values,
the first `m` equations form an invertible Vandermonde matrix times nonzero
positive integral weights.  Thus roots cannot escape while their top
coefficients remain bounded.  Any boundary point inside `A_n` therefore comes
from an actual collision, hence from a proper coarsening of the second
partition.  Replacing that maximal partition by a proper coarsening decreases
its length, so the Mason defect only increases.  No maximal closure is
contained in another.

Consequently,

\[
\boxed{
\text{the irreducible components of }\overline{\mathcal N_n}
\text{ are indexed by the partitions of }n\text{ using only }2\text{ and }3.}
\]

Writing `C_lambda=closure(E_lambda)` for a maximal 2/3 partition, the closure
order is exact: for every full-contact partition `nu`,

\[
\boxed{
\mathcal E_\nu\subseteq\mathcal C_\lambda
\quad\Longleftrightarrow\quad
\lambda\preceq\nu.}
\]

The forward implication is the collision theorem.  For the converse,
finiteness of the top-coefficient map turns any limiting factorization into a
genuine coarsening, while Mason excludes a second omitted value of a
noncoarsening type.

For `lambda=2^a3^b`, write `C_lambda=closure(E_lambda)`.  The component
dimensions are

\[
\dim\mathcal C_\lambda=a+b-1,
\qquad
\operatorname{codim}_{\mathcal A_n}\mathcal C_\lambda
=n-a-b-2=a+2b-2.
\]

The unique top-dimensional component minimizes `b`, so

\[
\operatorname{codim}_{\mathcal A_n}\overline{\mathcal N_n}
=\left\lceil{n\over2}\right\rceil-2.
\]

The same finiteness and Mason argument gives the set-theoretic intersection
formula inside `A_n`:

\[
\mathcal C_\lambda\cap\mathcal C_\mu
=\bigcup_{\substack{\lambda\preceq\nu\\\mu\preceq\nu}}
\mathcal E_\nu.
\]

Indeed, the two limiting omitted polynomials cannot have distinct target
values by Mason.  They are therefore equal, and their exact multiplicity
partition is a common coarsening.  The reverse inclusion is the collision-
closure theorem.

These are set-theoretic statements.  The scheme-theoretic intersections may
carry multiplicities or embedded components and are not identified here.

In particular, the number of components is

\[
\#\left\{b:0\le b\le\left\lfloor{n\over3}\right\rfloor,
\quad b\equiv n\pmod 2\right\},
\]

where `b` is the number of triple parts and `(n-3b)/2` is the number of
double parts.

If `m=floor(n/3)`, this count is

\[
c_n=
\begin{cases}
\lfloor m/2\rfloor+1,&n\text{ even},\\
\lfloor(m+1)/2\rfloor,&n\text{ odd},
\end{cases}
\qquad
\sum_{n\ge0}c_nz^n={1\over(1-z^2)(1-z^3)}.
\]

## Executable certificate

Run:

```bash
python scripts/verify_uniform_exceptional_seed_theorem.py
python scripts/verify_maximal_phi_irreducibility.py
python scripts/verify_contact_atom_principle.py
```

The script checks the weighted-Vandermonde identity symbolically, audits the
collision order and specialization identities through degree ten, checks the
Mason defect for maximal pairs through degree fifty, and retains the exact
degree-eight intersection calculation as a regression.  The irreducibility
script verifies the two stable endpoint-coordinate arguments, all seven small
cases, and exact quotient-coordinate factorizations through degree fourteen.
The atom script checks all distinct full-contact pairs through degree
twenty-four and the threshold-`r` generalization.
