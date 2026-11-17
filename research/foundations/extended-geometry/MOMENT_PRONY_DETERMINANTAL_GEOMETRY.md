# Moment--Prony determinantal geometry

The Gaussian Moments Conjecture is already false in three real variables.
The useful role of Gaussian moments here is therefore not another
counterexample search.  It is that the optimal mixed moments are polynomial
coordinates on the normalized seed space, so they can carry the entire
Keller-geometric package.

This note replaces equation-by-equation triangular transport by three
intrinsic matrix constructions:

1. a log-Prony Toeplitz--Hessenberg matrix for equal-multiplicity omitted
   partitions; and
2. saturated Christoffel--Hankel/Sylvester matrices for mixed omitted
   partitions; and
3. a Ritt--Krylov coefficient matrix for vertical composition strata.

Both are constructed from the optimal moment vector itself.  Degree six and
degree eight are computed exactly.  The scheme structures are compared, not
only the reduced zero sets.

There is also a precise limitation.  Ordinary Hankel rank remembers the
number of Prony nodes but forgets their fixed multiplicity weights.  For a
mixed partition such as `3+3+2`, the natural marked fixed-weight incidence has
a compact subresultant presentation, but its naive Fitting minors thicken the
total-collision fiber.  Thus the strongest version of the proposed
determinantal conjecture is false without a reduced-image or saturation
clause.

## 1. Optimal moments and the reversed omitted polynomial

Fix characteristic zero, degree `N>=4`, and bridge scale one.  Write

\[
 g(U)=U+\sum_{j=3}^{N-1}\mu_jU^j,\qquad
 \eta=g^{\langle-1\rangle},\qquad h(Z)=\frac Z{\eta(Z)}.
 \tag{1}
\]

The coefficients of `h=1+H` through degree `N-2` follow from the reciprocal
Toeplitz recurrence.  The endpoint equations

\[
 H(1)=0,\qquad H'(1)=-1
 \tag{2}
\]

then recover the last two coefficients.  Hence (1)--(2) construct

\[
 H_{\boldsymbol\mu}(W)=\sum_{j=2}^N c_j(\boldsymbol\mu)W^j
 \tag{3}
\]

directly from
`\boldsymbol\mu=(\mu_3,\ldots,\mu_{N-1})`.

Work on the exact-degree open `c_N!=0`.  The coefficients relevant to an
omitted polynomial are those of degree at least two: its linear and constant
terms can be changed by the pencil parameters.  Reverse and monicize these
known coefficients:

\[
 A_{\boldsymbol\mu}(T)
 =1+\sum_{j=1}^{N-2}
   \frac{c_{N-j}(\boldsymbol\mu)}{c_N(\boldsymbol\mu)}T^j.
 \tag{4}
\]

No seed-locus equation has entered (4).  It is computed from the mixed
moments by the recurrence (1).

## 2. The log-Prony Hessenberg matrix

Let `m>=2` divide `N`, and put `r=N/m`.  Define

\[
 \frac{A_{\boldsymbol\mu}'(T)}{A_{\boldsymbol\mu}(T)}
 =\sum_{j\ge1}\ell_jT^{j-1},\qquad \rho_j=\frac{\ell_j}{m}.
 \tag{5}
\]

For every `n`, form the Toeplitz--Hessenberg matrix

\[
 K_n^{(m)}=
 \begin{pmatrix}
 \rho _1&-1&0&\cdots&0\\
 \rho _2&\rho _1&-2&\ddots&\vdots\\
 \rho _3&\rho _2&\rho _1&\ddots&0\\
 \vdots&\ddots&\ddots&\ddots&-(n-1)\\
 \rho _n&\rho _{n-1}&\cdots&\rho _2&\rho _1
 \end{pmatrix}.
 \tag{6}
\]

If

\[
 A_{\boldsymbol\mu}(T)^{1/m}
 =1+\sum_{n\ge1}b_nT^n,
 \tag{7}
\]

then logarithmic differentiation gives

\[
 nb_n=\sum_{j=1}^n\rho_jb_{n-j}.
 \tag{8}
\]

The continuant identity for (8) is

\[
 \boxed{\det K_n^{(m)}=n!\,b_n.}
 \tag{9}
\]

This is the appropriate moment matrix.  A standard Hankel matrix of root
power sums is too coarse; (6) includes the fixed multiplicity through the
factor `1/m`.

### Theorem 2.1 -- pure-partition Moment--Prony equations

On `c_N!=0`, the scheme of normalized seeds for which an affine translate
of `H` is an `m`-th power of a degree-`r` polynomial is

\[
 \boxed{
 \det K_{r+1}^{(m)}
 =\det K_{r+2}^{(m)}
 =\cdots
 =\det K_{N-2}^{(m)}=0.}
 \tag{10}
\]

Equivalently, its ideal is `(b_(r+1),...,b_(N-2))`.

This is scheme-theoretic.  Indeed, the formal-root recurrence has

\[
 b_n=\frac1m a_n+F_n(a_1,\ldots,a_{n-1}),
 \tag{11}
\]

where `A=1+\sum a_nT^n`.  Thus
`(a_1,\ldots,a_(N-2)) -> (b_1,\ldots,b_(N-2))`
is a triangular polynomial automorphism after inverting `m`.  In the
`b`-coordinates, (10) is simply the reduced coordinate subspace

\[
 b_{r+1}=\cdots=b_{N-2}=0.
 \]

This also explains the exact-degree qualification.  Clearing powers of
`c_N` gives equations on the affine closure, but their unsaturated boundary
may contain spurious components.

## 3. Degree six from consecutive leading minors

Use

\[
 (u,v,\omega)=(\mu_3,\mu_4,\mu_5).
 \]

The recurrence gives

\[
 c_2=u,\qquad c_3=v,\qquad c_4=\omega-2u^2,
 \tag{12}
\]

followed by the two endpoint coefficients.

For the all-double partition `2+2+2`, Theorem 2.1 gives one equation:

\[
 \det K_4^{(2)}=0.
 \tag{13}
\]

After primitive denominator clearing, (13) is the negative of the previously
displayed sextic `\mathscr F_6`.  It was obtained here without substituting
moments into that old equation.

For the all-triple partition `3+3`, the two consecutive equations are

\[
 \boxed{\det K_3^{(3)}=\det K_4^{(3)}=0.}
 \tag{14}
\]

The first primitive determinant is the negative of `\mathscr G_6`, the
degree-six `3 o 2` Ritt hypersurface.  The second cuts that surface down to
the actual all-triple curve.  Thus the earlier equations acquire an intrinsic
interpretation:

\[
 \begin{aligned}
 \mathcal C_{2+2+2}&=V(\det K_4^{(2)}),\\
 \mathcal C_{3+3}&=V(\det K_3^{(3)},\det K_4^{(3)}).
 \end{aligned}
 \tag{15}
\]

The nilpotent intersection is also visible in the determinant presentation.
On the rational normalization parameter `q` of the all-triple curve,

\[
 \det K_4^{(2)}
 =\text{unit}\cdot
 (45q^4-30q^3+15q^2-6q+1)^2.
 \tag{16}
\]

The quartic factor is squarefree.  Hence each of the four all-six collision
points has transverse algebra `k[\epsilon]/(\epsilon^2)`.  The new ideal has
not reduced away the primitive dual number.

## 4. The Ritt--Krylov matrix

Let `N=ab` and let `f_{\boldsymbol\mu}` be the monic polynomial obtained from
`H_{\boldsymbol\mu}` after discarding its affine terms.  There is a unique
monic approximate inner polynomial

\[
 B_{\boldsymbol\mu}(W)
 =W^b+\beta_{b-1}W^{b-1}+\cdots+\beta_1W
 \tag{17}
\]

whose successive coefficients make the coefficients of
`f-B^a` in degrees `N-1,...,N-b+1` vanish.  This is again a triangular
recurrence in the optimal moments, not elimination.

Let `\operatorname{coeff}(P)` denote the column vector of coefficients in
degrees `0,...,N`.  Define the `(N+1) x (a+3)` matrix

\[
 \boxed{
 \mathcal K_{a,b}(\boldsymbol\mu)=
 \left[
 \operatorname{coeff}(1)\ \ 
 \operatorname{coeff}(W)\ \ 
 \operatorname{coeff}(B)\ \cdots\
 \operatorname{coeff}(B^a)\ \ 
 \operatorname{coeff}(f)
 \right].}
 \tag{18}
\]

### Theorem 4.1 -- Ritt--Krylov determinantal equations

On `c_N!=0`, the vertical `(a,b)` Ritt-composition locus is the maximal-minor
scheme of `\mathcal K_(a,b)`.

Indeed, rank drop in (18) says exactly

\[
 f-\alpha W-\beta\in k[B],                         \tag{19}
\]

which is the desired affine vertical composition.  This statement is also
scheme-theoretic.  The columns before the last have a unit pivot minor on
rows

\[
 0,1,b,2b,\ldots,ab.                               \tag{20}
\]

Row reduction by this unit minor shows that all maximal minors generate the
same ideal as the fixed-pivot minors obtained by adjoining one nonpivot row.
The top `b-1` such minors vanish from the approximate-root recurrence; the
remaining

\[
 N-a-b
\]

fixed-pivot minors are the canonical residual equations.  No radical or
set-theoretic replacement occurs.

The construction is uniform in `(a,b,N)` and direct in the optimal moment
coordinates.

## 5. Degree eight without seed-equation elimination

Degree eight has five optimal coordinates

\[
 (\mu_3,\mu_4,\mu_5,\mu_6,\mu_7).
 \]

The all-double locus `2+2+2+2` is

\[
 \boxed{\det K_5^{(2)}=\det K_6^{(2)}=0.}          \tag{21}
\]

After primitive expansion, these equations contain respectively `686` and
`1317` monomials.  Their compact representation is therefore the two leading
minors, not a printed elimination polynomial.

The two Ritt factor orders are the fixed-pivot maximal-minor ideals of

\[
 \mathcal K_{2,4}(\boldsymbol\mu),\qquad
 \mathcal K_{4,2}(\boldsymbol\mu).                 \tag{22}
\]

They each have two generators.  In universal monic coefficient coordinates,
if `R_5,R_6` are the two square-root truncation equations and `D_2,D_3`
are the two fixed-pivot minors of `\mathcal K_(2,4)`, exact calculation gives

\[
 D_3=R_5,\qquad R_6=D_2-2c_7R_5.                  \tag{23}
\]

Thus

\[
 (R_5,R_6)=(D_2,D_3)                              \tag{24}
\]

as ideals.  This proves scheme-theoretically, before any coordinate
specialization, that the degree-eight all-double locus equals the
`2 o 4` Ritt locus.  Applying the moment-coordinate isomorphism preserves
(24), including every nonreduced base change and intersection.

For the intersection of the two Ritt orders, use the `2 o 4` chart

\[
 B=W^4+u_3W^3+u_2W^2+u_1W,\qquad f=B^2+a_1B.
 \tag{25}
\]

The two pulled-back `4 o 2` minors are

\[
 \frac14(u_2-u_3^2)E,\qquad \frac14E,
 \quad
 E=8u_1-4u_2u_3+u_3^3.                            \tag{26}
\]

Consequently their scheme ideal is exactly `(E)`, not `(E^2)` or merely
`\sqrt{(E)}`.  Since `\partial E/\partial u_1=8`, it is reduced.  This is the
degree-eight `2 o 2 o 2` refinement locus, recovered from the moment matrices
without eliminating back to a seed equation.

## 6. Why ordinary Hankel minors do not prove the full conjecture

For a polynomial with distinct nodes `x_i` and multiplicities `m_i`, its
root power sums are

\[
 p_k=\sum_i m_ix_i^k.                              \tag{27}
\]

The ordinary Hankel matrix `(p_(i+j))` detects the number of distinct nodes,
but not the prescribed weights.  For example, the weight vectors

\[
 (2,2,2)\quad\text{and}\quad(4,1,1)
\]

on the same three nodes both give Hankel rank three.  The second polynomial
is not a square.  Thus standard Hankel minors alone cannot define the
all-double component.

Equal weights are repaired by the factor `1/m` in (5), which converts the
weighted logarithm into a unit-weight root sequence.  Unequal weights require
more data.

### The marked `3+3+2` incidence

Distinguish the weight-two node `x`.  For power sums `p_1,...,p_6`, set

\[
 q_k=\frac{p_k-2x^k}{3}.                           \tag{28}
\]

Let `e_n(q)` be obtained from the Newton recurrence

\[
 ne_n=\sum_{k=1}^n(-1)^{k-1}e_{n-k}q_k.
 \tag{29}
\]

Then the degree-eight `3+3+2` marked Prony incidence is

\[
 \boxed{e_3=e_4=e_5=e_6=0.}                       \tag{30}
\]

These equations say that `(q_k)` comes from two unit-weight nodes.  The
optimal moments enter directly: compute `p_1,...,p_6` from (4) by Newton's
identities and substitute them in (30).  The result is four equations over
the five optimal degree-eight moments and the one marked node `x`; no seed
equation is translated or eliminated.

Because `e_3` is cubic in `x` with constant leading coefficient, the incidence
is finite over moment space.  Reducing `e_4,e_5,e_6` modulo `e_3` and taking
their three multiplication matrices in the free cubic algebra gives one
`3 x 9` generalized Sylvester/Fitting matrix.  Its `3 x 3` minors describe the
support of the finite image.

They do **not** automatically give the reduced component scheme.  At total
collision

\[
 p_k=8t^k,
\]

all four relations in (30) are divisible by

\[
 (x-t)^3.                                          \tag{31}
\]

The marked incidence has a cubic fiber there.  Fitting minors retain this
thickness, whereas the irreducible omitted-partition component itself is
reduced.  Replacing the scheme-theoretic image (the annihilator/elimination
ideal) by the Fitting ideal would therefore change nilpotents precisely where
the partition geometry is most delicate.

This is not a defect of moment coordinates.  It is the familiar distinction
between:

* the reduced or scheme-theoretic image of a finite Prony incidence; and
* a Fitting scheme carrying its fiber lengths.

## 7. The unmarked Christoffel--Hankel construction

The marked Fitting obstruction does not prevent a uniform unmarked answer.
It says that the collision closure must be taken **after** reconstructing
support and weights on the distinct-node open.

Let

\[
 \lambda=(m_1,\ldots,m_r),\qquad m_i\ge2,\qquad
 N=\sum_i m_i                                      \tag{32}
\]

be any omitted partition that is not all-double.  Then `N>=2r+1`, so the
optimal reversed polynomial supplies root power sums at least through
`p_(2r-1)`.  Put `p_0=N` and form

\[
 \mathsf H=(p_{i+j})_{0\le i,j<r},\qquad
 \mathsf H^+=(p_{i+j+1})_{0\le i,j<r}.             \tag{33}
\]

The generalized-eigenvalue or matrix-pencil form of Prony recovery is
standard; see Harmouch--Khalil--Mourrain,
[*Structured low rank decomposition of multivariate Hankel
matrices*](https://arxiv.org/abs/1701.05805).  In the present univariate exact
setting, define

\[
 \begin{aligned}
 d&=\det\mathsf H,\\
 D_\lambda(T)&=\det(T\mathsf H-\mathsf H^+),\\
 v(T)&=(1,T,\ldots,T^{r-1})^{\mathsf t},\\
 \Xi_\lambda(T)&=v(T)^{\mathsf t}
                  \operatorname{adj}(\mathsf H)v(T).
 \end{aligned}                                     \tag{34}
\]

On the open `d!=0`, the moment sequence has a unique `r`-node Prony
decomposition

\[
 p_k=\sum_{i=1}^r w_ix_i^k.                        \tag{35}
\]

The Vandermonde factorization gives the two identities

\[
 \boxed{
 D_\lambda(T)=d\prod_i(T-x_i),\qquad
 \Xi_\lambda(x_i)=\frac d{w_i}.}                  \tag{36}
\]

For every distinct multiplicity value `m` occurring `n_m` times in
`\lambda`, set

\[
 C_m(T)=m\Xi_\lambda(T)-d.                         \tag{37}
\]

By (36), the support polynomial and `C_m` have exactly the nodes of weight
`m` in common.  Therefore the fixed-weight condition is the Sylvester rank
drop

\[
 \boxed{
 \operatorname{rank}\operatorname{Sylv}
 (D_\lambda,C_m)
 \le \deg D_\lambda+\deg C_m-n_m
 \quad\text{for every }m.}                        \tag{38}
\]

Equivalently, all minors one size larger than the bound in (38) vanish.

If the optimal moments supply later power sums `p_k`, for `k>=2r`, impose
the support recurrence.  Writing

\[
 D_\lambda(T)=\sum_{j=0}^r d_jT^j,
\]

these flat-extension equations are

\[
 \boxed{
 \sum_{j=0}^r d_jp_{k-r+j}=0,
 \qquad 2r\le k\le N-2.}                           \tag{39}
\]

Let `J_\lambda` be generated by the minors in (38) and the equations (39).
The required collision closure is

\[
 \boxed{
 I_\lambda^{\mathrm{mom}}
 =J_\lambda:(\det\mathsf H)^\infty.}              \tag{40}
\]

Equation (40), not the marked Fitting ideal, is the uniform unmarked
subresultant presentation.  In the polynomial affine moment ring, first clear
the powers of `c_N` introduced by the reversed monic polynomial and then
saturate by `c_N\det\mathsf H`; (40) is the equivalent formulation inside
the exact-degree localization.

### Theorem 7.1 -- all omitted partitions

On the exact-degree optimal-moment chart, every omitted-partition component
has a uniform intrinsic determinantal ideal:

1. the all-double partition is given by the log-Prony minors (10);
2. every other partition is given by the saturated
   Christoffel--Hankel/subresultant ideal (40).

These are equalities of schemes.

To prove the scheme assertion, pass to `d!=0`.  The ordered distinct-node
cover is finite etale, and (36) gives regular node and weight coordinates.
The Sylvester rank conditions say that at least `n_m` of the weight
coordinates equal `m`.  Since `\sum_m n_m=r`, the conditions assign every
node exactly one prescribed weight.  Locally the first nonzero
subresultant minors contain the corresponding linear weight differences, so
the locus is reduced.  Permuting the ordered nodes acts transitively on all
assignments with multiplicity multiset `\lambda`; the etale quotient is
therefore the irreducible reduced fixed-weight Prony locus.  Contracting its
ideal from `d!=0` is precisely the saturation (40), hence gives its prime
collision closure.

This proof also explains why saturation is essential.  It discards the
fiber-length scheme of the marked map, not the genuine collision points.
Those points return in the closure with the reduced component structure.
Nilpotents in intersections are then obtained correctly by summing the
prime component ideals.

## 8. Degree five and eight mixed tests

The shortest mixed case is `3+2` in degree five.  Here `r=2` and the optimal
data end at `p_3=p_(2r-1)`.  No flat-extension equation is needed: the two
Sylvester resultant conditions already give the full component.  At the
opposite extreme, a one-node partition `(N)` needs no weight resultant at
all: `p_0=N` fixes its sole weight, and (39) is the complete support
recurrence.

For `3+3+2` in degree eight, the unmarked matrices are

\[
 \mathsf H,\mathsf H^+\in\operatorname{Mat}_{3\times3},
 \tag{41}
\]

and the two Sylvester conditions are

\[
 \begin{aligned}
 \operatorname{rank}\operatorname{Sylv}(D,C_2)&\le6,\\
 \operatorname{rank}\operatorname{Sylv}(D,C_3)&\le5.
 \end{aligned}                                    \tag{42}
\]

Both Sylvester matrices are `7 x 7`.  There is one flat-extension equation,
for `p_6`.  The component ideal is the saturation of these minors and that
extension equation by `det(H)`.

For the parameterization

\[
 p_k=2x^k+3y^k+3z^k,
\]

the checker obtains

\[
 \det\mathsf H
 =18(x-y)^2(x-z)^2(y-z)^2,\qquad
 D(T)=\det\mathsf H\,(T-x)(T-y)(T-z),              \tag{43}
\]

and the two Sylvester matrices have exact ranks six and five at a distinct
rational test point.

The mixed/all-double scheme intersection retains the expected primitive
dual number.  Pulling the two degree-eight all-double equations back to the
`3+3+2` normalization gives

\[
 \begin{aligned}
 R_5&=-3(y-z)^4(2x-y-z),\\
 R_6&=-(y-z)^4
 (12xy+12xz-7y^2-10yz-7z^2).
 \end{aligned}                                    \tag{44}
\]

The unordered transverse normalization coordinate is
`\delta=(y-z)^2`.  Thus (44) contains `\delta^2`, proving generic
intersection length two along the `6+2` common coarsening.  The saturation in
(40) has removed the spurious cubic marked fiber while preserving the
correct quadratic intersection nilpotent.

## 9. Completed Moment--Prony theorem

Combining Theorems 2.1, 4.1, and 7.1 gives the proposed result in corrected
scheme-theoretic form.

> **Moment--Prony determinantal theorem.**
>
> On every fixed-scale exact-degree optimal Gaussian moment chart:
>
> 1. every omitted-value partition stratum has a uniform intrinsic
>    log-Hessenberg or saturated Christoffel--Hankel/subresultant ideal;
> 2. every vertical Ritt-composition stratum has the uniform Ritt--Krylov
>    maximal-minor ideal (18); and
> 3. sums of these exact component ideals compute their intersections with
>    the intended nilpotents.

The saturation clause is mathematically necessary and is the only correction
to the original “one moment matrix” formulation.  Unsaturated ordinary
Hankel or marked Fitting minors can remember fiber length instead of the
reduced component closure.

This is the direct connection to algebraic statistics and sparse moment
recovery: Keller partitions are fixed-weight Prony models, Ritt loci are
Krylov rank-drop models, and the collision algebra distinguishes component
closures from fiber-length Fitting schemes.

## 10. Reproduction

The reusable implementation is
[`jcsearch/moment_prony.py`](../jcsearch/moment_prony.py).  Run

```bash
.venv/bin/python scripts/verify_moment_prony_determinantal_geometry.py
```

The checker verifies:

1. ordinary Hankel-rank failure for fixed weights;
2. the degree-six all-double and `3+3` leading-minor ideals;
3. exact agreement with the independent degree-six Ritt equations;
4. the four length-two degree-six collision intersections;
5. both degree-eight all-double leading minors in five optimal moments;
6. both degree-eight Ritt--Krylov ideals;
7. scheme equality of the all-double and `2 o 4` ideals;
8. reducedness of the degree-eight Ritt intersection ideal;
9. the marked `3+3+2` subresultant incidence; and
10. the cubic total-collision thickness obstructing a naive global Fitting
    claim;
11. the saturated unmarked Christoffel--Hankel presentation for `3+3+2`;
12. the minimal-length degree-five `3+2` presentation;
13. the one-node support recurrence; and
14. the length-two degree-eight mixed/all-double intersection.
