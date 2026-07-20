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
\bigcup_{\substack{\lambda\vdash n\\\lambda_i\ge2}}\mathcal E_\lambda.}
\]

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
the direct `C=0` calculation, this proves the union formula.

The root parameter space has dimension `ell(lambda)`.  The polynomial
`Phi_lambda` is nonzero—for example its specialization with all roots zero is
`1`—so its hypersurface has dimension `ell(lambda)-1`.  The excluded factors
are proper divisors, and tangent-chord normalization supplies points in the
admissible open set.  Quotienting equal parts is finite and does not change
dimension.

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

Let `lambda!=mu` be partitions of `n` whose parts are only two and three.
Away from their collision discriminants, the affine-difference system is
empty:

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

If `b_lambda` denotes the number of parts equal to three, then

\[
\ell(\lambda)={n-b_\lambda\over2}.
\]

Two distinct feasible 2/3 partitions have different numbers of threes of the
same parity, so `b_lambda+b_mu>=2`.  Consequently

\[
\ell(\lambda)+\ell(\mu)
=n-{b_\lambda+b_\mu\over2}\le n-1,
\]

contradicting Mason--Stothers.

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

This proves the proposed component poset at the level of full-contact strata
and their closures.  To identify each maximal stratum closure with one
irreducible component in the strict scheme-theoretic sense, one additional
issue remains: prove irreducibility of `Phi_lambda` (or determine its
irreducible factors) uniformly for every maximal 2/3 partition.

## Executable certificate

Run:

```bash
python scripts/verify_uniform_exceptional_seed_theorem.py
```

The script checks the weighted-Vandermonde identity symbolically, audits the
collision order and specialization identities through degree ten, checks the
Mason defect for all maximal pairs through degree fifty, and retains the exact
degree-eight intersection calculation as a regression.
