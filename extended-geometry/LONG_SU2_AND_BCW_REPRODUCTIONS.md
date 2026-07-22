# Complete local reproductions: `SU(2)`, the BCW 79-variable route, and a 22-variable optimization

This note supplies the two proofs left deliberately conditional in the first
external-consequences audit.  It has two distinct outcomes.

1. The `SU(2)` integration formula and Long's lifted witness are now proved
   locally from normalized surface measure on the unit three-sphere.
2. The complete `3 -> 39 -> 79` Bass--Connell--Wright route is constructed and
   checked exactly.  The fixed-dimensional implication from that map to
   `not GMC(158)` is proved locally in a companion note, following and
   crediting Derksen--van den Essen--Zhao and Zhao.
3. A repository-derived common-factor optimization replaces the conservative
   degree-lowering stage by `3 -> 16`.  The cubic component vector has exact
   rational rank seven, so rank-compressed homogenization needs only 24
   variables.  Its two-dimensional constant-kernel quotient gives a
   22-variable cubic-homogeneous collision and improves the route-based
   consequence to `not GMC(44)`.
   This is not a formula or dimension claim attributed to Long.

Neither local proof changes the provenance of Christopher D. Long's external
results or constitutes external review.

## 1. A self-contained `SU(2)` integration proof

### 1.1 `SU(2)` as the unit three-sphere

Every element of `SU(2)` has a unique presentation

\[
 g(\alpha,\beta)=
 \begin{pmatrix}
  \alpha&-\overline\beta\\
  \beta&\overline\alpha
 \end{pmatrix},
 \qquad |\alpha|^2+|\beta|^2=1.
\]

Thus `SU(2)` is the unit sphere `S^3` in `C^2=R^4`.  Under the standard
identification with the unit quaternions, left multiplication by
`q=(a,b,c,d)` is represented on `R^4` by

\[
 L_q=
 \begin{pmatrix}
 a&-b&-c&-d\\
 b&a&-d&c\\
 c&d&a&-b\\
 d&-c&b&a
 \end{pmatrix}.
\]

Direct multiplication gives

\[
 L_q^T L_q=(a^2+b^2+c^2+d^2)I_4.
\]

For a unit quaternion, left multiplication is therefore orthogonal.  The
normalized surface measure on `S^3` is left invariant, has total mass one,
and hence is the normalized Haar measure on `SU(2)`.  This uses only the
uniqueness of normalized Haar probability on a compact group.

### 1.2 Exact Hopf-coordinate density

Away from the two measure-zero coordinate circles, put

\[
 \alpha=\sqrt{1-x}\,e^{i\theta_2},
 \qquad
 \beta=\sqrt{x}\,e^{i\theta_1},
 \qquad 0<x<1.
\]

As a map to `R^4`, this is

\[
 \Phi(x,\theta_1,\theta_2)=
 \left(
 \sqrt{1-x}\cos\theta_2,
 \sqrt{1-x}\sin\theta_2,
 \sqrt{x}\cos\theta_1,
 \sqrt{x}\sin\theta_1
 \right).
\]

The three coordinate derivatives are orthogonal and their squared norms are

\[
 \frac{1}{4x(1-x)},\qquad x,\qquad 1-x.
\]

Consequently the Gram determinant is `1/4`, so the surface-volume element is

\[
 \frac12\,dx\,d\theta_1\,d\theta_2.
\]

Integrating this density over `0<=x<=1` and both phases gives total area
`2 pi^2`.  Dividing by that total, normalized Haar measure is exactly

\[
 dx\,\frac{d\theta_1}{2\pi}\,\frac{d\theta_2}{2\pi}.
\]

In particular, `x` is uniform on `[0,1]` and the two phases are independent
normalized circle variables.  The orthogonality and Gram calculations are
checked symbolically by
[`verify_long_su2_haar.py`](../scripts/verify_long_su2_haar.py).

### 1.3 Monomial integration

Use Long's coordinate order

\[
 g=\begin{pmatrix}a&c\\b&d\end{pmatrix},
 \qquad
 a=\alpha,\quad b=\beta,\quad
 c=-\overline\beta,\quad d=\overline\alpha.
\]

For nonnegative integers `r,s,t,u`, Hopf coordinates give

\[
 a^r b^s c^t d^u
 =(-1)^t(1-x)^{(r+u)/2}x^{(s+t)/2}
   e^{i(s-t)\theta_1}e^{i(r-u)\theta_2}.
\]

The two circle integrals vanish unless `r=u` and `s=t`.  In the surviving
case, the beta integral gives

\[
 \int_{SU(2)}a^r b^s c^t d^u\,dg
 =(-1)^s\delta_{r,u}\delta_{s,t}
   \int_0^1(1-x)^r x^s\,dx
 =(-1)^s\delta_{r,u}\delta_{s,t}
   \frac{r!s!}{(r+s+1)!}.                         \tag{1.1}
\]

Now define the polynomial torus substitution

\[
 \beta(z_1,z_2,x)=
 ((1-x)z_2,xz_1,-z_1^{-1},z_2^{-1}).              \tag{1.2}
\]

The monomial `a^r b^s c^t d^u` maps to

\[
 (-1)^t(1-x)^r x^s z_1^{s-t}z_2^{r-u}.
\]

Its torus/Beta integral is again exactly (1.1).  By linearity, for every
polynomial `P(a,b,c,d)`,

\[
 \int_{SU(2)}P\,dg
 =\int_0^1\int_{T^2}P(\beta(z_1,z_2,x))
   \frac{dz_1}{2\pi i z_1}\frac{dz_2}{2\pi i z_2}\,dx.             \tag{1.3}
\]

This is a complete local proof of the Müger--Tuset formula used by Long, not
merely a check conditional on that formula.  Müger and Tuset remain the cited
external source for the formula.

### 1.4 Long's witness

For

\[
 F=(1+c)(ad+b),\qquad G=-c,
\]

the substitution (1.2) gives

\[
 F\longmapsto(1-z_1^{-1})((1-x)+xz_1),
 \qquad
 G\longmapsto z_1^{-1}.
\]

The general beta/binomial computation already proved in the
[external-consequences note](EXTERNAL_CONSEQUENCES_AND_PROVENANCE.md) now
applies through the locally proved formula (1.3), yielding for every `n>=1`

\[
 \int_{SU(2)}F^n\,dg=0,
 \qquad
 \int_{SU(2)}GF^n\,dg=\frac{(-1)^{n-1}}{n+1}\ne0.
\]

Thus the displayed `SU(2)` identities are independently reproduced in full.

## 2. The exact `3 -> 39 -> 79` BCW route

### 2.1 Starting normalization and support

Long uses the determinant-one presentation `L` which satisfies

\[
 L=\operatorname{diag}(1/2,1/2,-1/2)
   \circ F\circ\operatorname{diag}(1,2,2)
\]

for the repository's foundational determinant `-2` map `F`.  Composing the
target with `(u,v,w) -> (-w,v,u)` gives a map `G` with identity linear part:

\[
 \begin{aligned}
 G_1={}&x-3x^2y-x^3z,\\
 G_2={}&y+3xz+24xy^2+12x^2yz+36x^2y^3+12x^3y^2z,\\
 G_3={}&z+8y^2+6xyz+28xy^3+12x^2y^2z
          +24x^2y^4+8x^3y^3z.
 \end{aligned}
\]

Its terms of degrees `4,5,6,7` occur with exact counts

\[
 3,\quad2,\quad2,\quad1.                            \tag{2.1}
\]

The determinant, normalization, support counts, and three-point collision are
all checked exactly by the route script.

### 2.2 One stable degree-lowering step

Suppose one coordinate `f_i` contains a monomial `ab`, with `deg a=p`,
`deg b=q`, `p,q>=2`.  Adjoin variables `Y,Z` and replace the stabilized map by

\[
 \widetilde f_i=f_i-(Y+a)(Z+b),
 \qquad
 \widetilde f_{n+1}=Y+a,
 \qquad
 \widetilde f_{n+2}=Z+b,                            \tag{2.2}
\]

leaving the other coordinates unchanged.  This is an exact polynomial
left--right equivalence.  Indeed, first apply the source automorphism

\[
 (x,Y,Z)\longmapsto(x,Y+a(x),Z+b(x))
\]

to `f x id`, then apply the target automorphism which subtracts the product of
the last two outputs from output `i`.

The term `ab` cancels.  The new nonlinear terms have degrees at most

\[
 p+1,\quad q+1,\quad p,\quad q,\quad2.              \tag{2.3}
\]

The determinant, generic degree, noninvertibility, and collision schemes are
preserved.  On a collision point `x`, adjoining `Y=-a(x), Z=-b(x)` transports
the collision explicitly.

### 2.3 Eighteen balanced steps

Let `c(d)=0` for `d<=3`.  Factoring a degree-`d` monomial into balanced degrees
`p+q=d` and using (2.3) gives

\[
 c(d)\leq1+c(p+1)+c(q+1)+c(p)+c(q).
\]

The choices

\[
 4=2+2,\quad5=2+3,\quad6=3+3,\quad7=3+4
\]

give the upper bounds

\[
 c(4)\leq1,\qquad c(5)\leq2,\qquad
 c(6)\leq3,\qquad c(7)\leq5.
\]

Together with (2.1), the total is

\[
 3c(4)+2c(5)+2c(6)+c(7)\leq18.
\]

The exact implementation performs the degree sequence

\[
 7,6,6,5,5,5,\underbrace{4,\ldots,4}_{12\text{ times}},
\]

and terminates with a noninjective determinant-one map

\[
 K(X)=X+Q(X)+C(X),\qquad X\in\mathbb A^{39},        \tag{2.4}
\]

where `Q` is quadratic homogeneous and `C` is cubic homogeneous.  Since each
step adds two variables, the dimension is exactly

\[
 3+2\cdot18=39.
\]

### 2.4 Nilpotent form in 78 variables

On variables `(X,Y)` in `A^39 x A^39`, define

\[
 U(X,Y)=(X+Q(X)+Y,\;Y-C(X))=(X,Y)+N(X,Y).           \tag{2.5}
\]

This map is stably equivalent to `K`.  Starting with `K x id`, apply

\[
 A(X,Y)=(X,Y-C(X))
\]

on the source and then

\[
 P(R,S)=(R+S,S)
\]

on the target.  The result is exactly (2.5).

It remains to prove that `JN` is nilpotent, not merely that `det(I+JN)=1`.
For an indeterminate `t`, put

\[
 E_t(X)=X+tQ(X)+t^2C(X)=t^{-1}K(tX).
\]

Hence `det DE_t=1`.  More explicitly,

\[
 (X,Y)+tN(X,Y)
 =P_t\circ(E_t\times\operatorname{id})\circ A_t(X,Y),
\]

where

\[
 A_t(X,Y)=(X,Y-tC(X)),\qquad P_t(R,S)=(R+tS,S).
\]

Therefore

\[
 \det(I+tJN)=1                                      \tag{2.6}
\]

as a polynomial identity in `t`.  All nonconstant coefficients of the
characteristic polynomial of `JN` vanish, so `JN` is nilpotent.

### 2.5 Cubic homogenization in 79 variables

Adjoin `T` and homogenize the degree-one, degree-two, and degree-three parts of
`N`:

\[
 H(X,Y,T)=(YT^2+Q(X)T,\;-C(X),\;0).
\]

Every component of `H` is cubic homogeneous.  The map

\[
 V(X,Y,T)=(X,Y,T)+H(X,Y,T)                          \tag{2.7}
\]

lives in

\[
 2\cdot39+1=79
\]

variables.  To verify its Jacobian, write `Z=(X,Y)`.  On `T!=0`,

\[
 H(Z,T)=T^3N(Z/T),
 \qquad
 \partial_ZH=T^2JN(Z/T).
\]

Equation (2.6), with `t=T^2`, gives `det DV=1` on this dense open, hence
everywhere by polynomial identity.

At `T=1`, map (2.7) restricts to `(U,1)`.  If `X_i` are the three transported
collision points of `K`, then

\[
 (X_i,C(X_i),1)
\]

are three distinct points with one common image under `V`.  Thus `V` is an
explicit noninvertible cubic-homogeneous Keller map in 79 variables.

The complete construction and collision are checked by
[`verify_long_bcw_79_route.py`](../scripts/verify_long_bcw_79_route.py).
That script writes the exact
[79-variable sparse artifact](../artifacts/generated-results/long_bcw_79_counterexample.json),
which is replayed without SymPy by
[`audit_long_bcw_79_independent.py`](../scripts/audit_long_bcw_79_independent.py).

### 2.6 The final Gaussian implication

For each fixed `r`, the needed implication is

\[
 \mathrm{GMC}(2r)\Longrightarrow\mathrm{SIC}(r).
\]

and `SIC(r)` forces a cubic-homogeneous Keller map in `r` variables to be
invertible.  The companion
[fixed-dimensional proof](FIXED_GMC_SIC_PROOF.md) derives both statements from
Gaussian contraction, an uncountable-field countable-union lemma, and a
coefficient proof of the formal inversion identity.  Applying its
contrapositive to (2.7), with `r=79`, gives

\[
 \neg\mathrm{GMC}(158).
\]

This is now locally proved, but it is still nonconstructive at the Gaussian
witness step because the passage from pointwise thresholds to one uniform
threshold uses a countable-union argument.  The proof reproduces only the
fixed-dimensional implication needed here, not all results of DVEZ or Zhao.

## 3. Shared-factor and rank-compressed optimization: `3 -> 16 -> 24`

The 79-variable construction above remains the exact reproduction of Long's
conservative route.  It is not dimension-minimal: its 18 steps expose two
fresh factors independently even when a factor has already appeared as an
output coordinate.

### 3.1 Reusable-factor elementary equivalence

Suppose a previous source shear has exposed a polynomial `a(x)` as the output

\[
 A=Y+a(x).                                           \tag{3.1}
\]

If coordinate `i` contains `c a(x)b(x)` and `i` is not the `A` coordinate,
adjoin only one new variable `Z`, apply the source shear

\[
 Z\longmapsto Z+b(x),
\]

and then the elementary target shear

\[
 T_i\longmapsto T_i-cA T_Z.                         \tag{3.2}
\]

After the source shear, `T_Z=Z+b(x)`, so (3.2) cancels `cab`.  Both
automorphisms have determinant one.  If `b` was exposed previously too, no
new variable is needed: subtract `cAB` directly.  If `a=b` was not exposed,
one variable suffices—expose `A=Y+a` and subtract `cA^2`.  These are exact
stable left--right equivalences, not arithmetic-circuit substitutions made
outside the map category.

At a transported collision point, set every new source variable to the
negative value of its exposed factor.  All exposed-factor outputs are then
zero, so each target shear preserves the common image exactly.

### 3.2 Deterministic factor search

For each current top-degree monomial, enumerate every unordered factor split
`ab` with both degrees between two and `d-2`.  Retain a registry of outputs of
the form (3.1), and rank candidates lexicographically by

\[
 (\text{new maximum degree},
   \sum_{\deg m>3}(\deg m-3)^2,\
   \text{number of high terms},\
   \text{new variables},\
   -\text{reusable factors}).                        \tag{3.3}
\]

Freezing the best trace from a deterministic width-24 beam search produces 17
target cancellations with degree
sequence

\[
 7,6,6,5,5,5,5,\underbrace{4,\ldots,4}_{10\text{ times}},
\]

but the corresponding new-variable counts are

\[
 2,0,1,0,1,1,1,0,0,1,1,1,1,1,1,0,1.              \tag{3.4}
\]

Their sum is 13.  The resulting determinant-one map therefore has degree at
most three in

\[
 3+13=16                                             \tag{3.5}
\]

variables, identity linear part, and the transported rational three-point
collision.  Five zero-cost cancellations in (3.4) reuse both exposed factors;
the square cancellation at the second step is the first of them.

This is a certified upper bound, not a minimality theorem.  The greedy score
is deliberately recorded so SAT, MILP, dynamic-programming, or beam searches
can seek a still smaller exposure registry without changing the certificate
format.

### 3.3 Rank-compressed cubic homogenization

Write the 16-variable map as `K=X+Q+C`, with `Q,C` homogeneous of degrees two
and three.  Let `k` be the row rank over `Q` of the coefficient matrix of the
component vector `C`.  Choose independent component polynomials
`c=(c_1,...,c_k)` and the unique constant matrix `B` such that

\[
 C(X)=B c(X).                                       \tag{3.6}
\]

Only `Y in A^k` is needed.  Put

\[
 U(X,Y)=(X+Q(X)+BY,Y-c(X)).                         \tag{3.7}
\]

This is stably left--right equivalent to `K`: precompose `K times id` with
`A(X,Y)=(X,Y-c(X))`, then postcompose with `P(R,S)=(R+BS,S)`.  More strongly,
for

\[
 E_t(X)=X+tQ(X)+t^2C(X)=t^{-1}K(tX),
\]

define `A_t(X,Y)=(X,Y-tc(X))` and
`P_t(R,S)=(R+tBS,S)`.  Direct substitution gives

\[
 \operatorname{id}+tN
 =P_t\circ(E_t\times\operatorname{id})\circ A_t,
 \qquad N=(Q+BY,-c).                                \tag{3.8}
\]

All source and target shears have determinant one, and
`det DE_t=det DK(tX)=1`; hence `det(I+tJN)=1`.  Equivalently, the Schur
complement in the homogenized Jacobian is

\[
 I+tJQ+t^2B Jc=DE_t.
\]

Therefore

\[
 V_{n+k+1}(X,Y,T)
 =(X,Y,T)+(TQ(X)+T^2BY,-c(X),0)                     \tag{3.9}
\]

is cubic homogeneous and Keller.  If `K(p)=q`, then
`V(p,c(p),1)=(q,0,1)`, so every collision transports.

For the frozen 16-variable trace, precisely the components numbered
`0,1,2,3,4,6,8` are nonzero.  Exact rational row reduction shows that all
seven are independent.  Thus `k=7`, not merely `k<=7`, and (3.9) has
`16+7+1=24` variables.  Its transported rational three-point collision makes
it noninvertible.  The locally proved fixed-dimensional implication now gives

\[
 \boxed{\neg\mathrm{GMC}(48)}.                     \tag{3.10}
\]

This improves only the nonexplicit route-based dimension bound.  Long's
direct three-real-Gaussian witness remains far stronger and independently
authored.

The SymPy generator
[`verify_shared_bcw_33_route.py`](../scripts/verify_shared_bcw_33_route.py)
writes the exact
[33-variable sparse artifact](../artifacts/generated-results/shared_bcw_33_counterexample.json).
The dependency-free
[`audit_shared_bcw_33_independent.py`](../scripts/audit_shared_bcw_33_independent.py)
replays all factor exposures and target shears from the original map and
reconstructs the cubic collision without importing the generator.
The general rank factorization is implemented in
[`rank_compressed_bcw_homogenization.py`](../scripts/rank_compressed_bcw_homogenization.py).
The generator
[`verify_rank_compressed_bcw_24_route.py`](../scripts/verify_rank_compressed_bcw_24_route.py)
writes the exact
[24-variable sparse artifact](../artifacts/generated-results/rank_compressed_bcw_24_counterexample.json),
and
[`audit_rank_compressed_bcw_24_independent.py`](../scripts/audit_rank_compressed_bcw_24_independent.py)
independently recomputes the rational rank, replays (3.8) with sparse
polynomials, reconstructs the map, and checks the collision using only the
standard library.

### 3.4 What remains after 48

The bound 48 is still an upper bound, not a minimality theorem.  The executable
[`search_rank_aware_bcw.py`](../scripts/search_rank_aware_bcw.py) reconstructs
the monomial-factor beam search, deduplicates exact polynomial states, and
scores the actual rational objective `s+rank(C)` after every candidate.  Both
the legacy degree-first and genuinely rank-first orderings, at width 128,
finish with `s+rank(C)=20`; neither finds a 23-variable map.  Across the 232
completed traces retained by the degree-first run, no trace passes eight exact
necessary samples of `det(I+sJQ+tJC)=1`.  These are finite search results, not
lower-bound proofs.

The present 16-variable map itself definitely cannot avoid doubling.  At
`X=(1,...,1)` its known scaling family still has

\[
 \det(I+JQ+JC)=\det(I+2JQ+4JC)=1,
\]

but the independent specializations give

\[
 \det(I+JC)=-4160,\qquad \det(I+JQ)=-78.            \tag{3.11}
\]

Thus the desired two-parameter identity fails before homogenization.  The
exact checker
[`verify_two_parameter_bcw_obstruction.py`](../scripts/verify_two_parameter_bcw_obstruction.py)
records this obstruction.  Further improvement now requires a broader
equivalence class—most plausibly polynomial-factor reuse and multi-term
cancellation—or a larger search that escapes the width-128 monomial beam.

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_long_su2_haar.py
.venv/bin/python scripts/verify_long_xz_mathieu.py
.venv/bin/python scripts/verify_long_bcw_79_route.py
python3 scripts/audit_long_bcw_79_independent.py
.venv/bin/python scripts/verify_shared_bcw_33_route.py
python3 scripts/audit_shared_bcw_33_independent.py
.venv/bin/python scripts/verify_rank_compressed_bcw_24_route.py
python3 scripts/audit_rank_compressed_bcw_24_independent.py
.venv/bin/python scripts/verify_constant_kernel_bcw_22_route.py
python3 scripts/audit_constant_kernel_bcw_22_independent.py
.venv/bin/python scripts/verify_two_parameter_bcw_obstruction.py
python3 scripts/verify_fixed_gmc_sic_bridge.py
```

The first two scripts jointly certify the complete `SU(2)` proof.  The next
two construct and independently replay Long's conservative 79-variable
route.  The next pair record and replay the repository's shared-factor
baseline, and the following pair construct and independently replay its
rank-compressed 24-variable homogenization.  The next pair construct and
independently replay its 22-variable constant-kernel quotient.  The next script checks the
two-parameter shortcut obstruction; the final script checks the coefficient
skeleton of the fixed-dimensional implication.  None of these
replaces the repository's separate 95-variable cubic-homogeneous artifact.

The external theorem inputs and construction sources are:

- Christopher D. Long, [*Small Counterexamples to the Gaussian Moments
  Conjecture*](https://arxiv.org/abs/2607.18186), arXiv:2607.18186v1;
- Christopher D. Long, [*Counterexamples to the (xz)-Conjecture and the
  Mathieu Conjecture for (SU(2))*](https://arxiv.org/abs/2607.19012),
  arXiv:2607.19012v1;
- Hyman Bass, Edwin H. Connell, and David Wright,
  [*The Jacobian Conjecture: Reduction of Degree and Formal Expansion of the
  Inverse*](https://doi.org/10.1090/S0273-0979-1982-15032-7), Bulletin of the
  AMS 7 (1982), 287--330;
- Harm Derksen, Arno van den Essen, and Wenhua Zhao,
  [*The Gaussian Moments Conjecture and the Jacobian
  Conjecture*](https://arxiv.org/abs/1506.05192), Israel Journal of Mathematics
  219 (2017), 917--928; and
- Wenhua Zhao, [*Images of commuting differential operators of order one with
  constant leading coefficients*](https://arxiv.org/abs/0902.0210), Journal
  of Algebra 324 (2010), 231--247.
