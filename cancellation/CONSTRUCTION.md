# The cancellation construction

This document contains the uniform map, its constant Jacobian, the exact polynomiality
criterion, generic reconstruction, and a symbolic full collision.  Arithmetic
of the parameter polynomial and boundary comparison are separate results in
[ARITHMETIC.md](ARITHMETIC.md) and
[BOUNDARY_GEOMETRY.md](BOUNDARY_GEOMETRY.md).

Throughout, `k` has characteristic zero, `m,r>=1`, `C in k^*`, and
`h in K[A]` for a finite separable extension `K/k`.  Put

\[
 A=1+xy^m,\qquad
 B=A^{r+1}z+y^{m+1}h(A),\qquad
 P=AB,\qquad Q=y+xB,                                      \tag{1}
\]

and initially in `K[x,y,z,A^{-1}]` put

\[
 R=C\int_0^{x/A}\bigl(1-t(Q-Pt)^m\bigr)^r\,dt.            \tag{2}
\]

The construction has two logically independent parts: the Jacobian identity
holds for every `h`, while polynomiality is a finite cancellation condition
on the jet of `h` at `A=0`.

## 1. Universal Jacobian

Set

\[
 s=x/A,\qquad y=Q-sP,\qquad D=1-s(Q-sP)^m.
\]

On the source,

\[
 D=A^{-1},\qquad x=s/D,qquad B=PD,
\]

and

\[
 z=PD^{r+2}-(Q-sP)^{m+1}D^{r+1}h(D^{-1}).                 \tag{3}
\]

These are precisely the cancellation inputs to the
[boundary-cancelled incidence lemma](CONTROLLED_BOUNDARY_SUSPENSIONS.md#1-boundary-cancelled-incidence-lemma).
The controlled divisor and the two chart factors are

\[
 D=1-s(Q-Ps)^m=A^{-1},\qquad
 R_s=CD^r,
\]

\[
 J_\alpha=
 \det\frac{\partial(s,P,Q)}{\partial(x,y,z)}
 =-A^r=-D^{-r},\qquad J_\beta=1.
\]

Thus the lemma has `u=C` and `kappa=-1`, and gives

\[
 \boxed{\det\frac{\partial(P,Q,R)}{\partial(x,y,z)}=-C}.   \tag{4}
\]

Neither `h` nor its derivatives enter this determinant ledger.

### Controlled-boundary plane core

The determinant proof has a plane core analogous to, but different from, the
weighted tangent-map core.  Put

\[
 D(s,P,Q)=1-s(Q-Ps)^m,
 \qquad
 R(s,P,Q)=C\int_0^s D(t,P,Q)^r\,dt.
\]

For each fixed `P`, define

\[
 \chi_{m,r;P}:\mathbb A^2_{s,Q}\longrightarrow\mathbb A^2_{Q,R},
 \qquad
 (s,Q)\longmapsto(Q,R(s,P,Q)).
\]

This is the coordinate-preserving plane core in the lemma. Its Jacobian is
the inverse-equation derivative, up to orientation:

\[
 \boxed{\det D\chi_{m,r;P}=-C D(s,P,Q)^r.}               \tag{4a}
\]

The common lemma now supplies the three-dimensional determinant; no second
chain-rule calculation is needed here. The family-specific issue is instead
the polynomiality gate: `R` is automatically polynomial in `(s,P,Q)`, while
the finite cancellation operator below decides whether it remains
polynomial after returning through the rational chart to `(x,y,z)`.

This should not be conflated with the weighted suspension.  There the
vertical maps are polynomial and the plane core has simple ramification; here
the compensating source chart is birational and cancels the arbitrary power
`D^r`.

The weighted and quadratic-gauge dictionaries, including their distinct
boundary normalizations and collision-fibre clauses, are recorded alongside
this one in the
[common incidence lemma](CONTROLLED_BOUNDARY_SUSPENSIONS.md#3-the-three-established-families-are-instances).

The same `(m,r)=(1,1)` divisor also admits a second polynomial plane gauge.
Its inverse equation has a quadratic target tilt, it contains the
foundational map as the seed `G(S)=S^3+S`, and it engineers a prescribed
complete root fiber in every degree.  The determinant, polynomial lift,
inverse reconstruction, `G_m` discriminant normalization, and symmetric
monodromy are proved in the
[root-engineered quadratic-gauge theorem](ROOT_ENGINEERED_QUADRATIC_GAUGE.md).
Its complete canonical boundary ledger and the resulting unit-rank
obstruction prove that, in every degree at least four, it is not stably
polynomially left--right equivalent to a boundary-clean weighted map; see the
[quadratic-versus-weighted stable-separation theorem](../verified/QUADRATIC_WEIGHTED_STABLE_SEPARATION.md).
On the coefficient-torus locus, its internal stable orbits are classified
exactly by a two-dimensional scaling action, giving moduli dimension `N-4`;
see the
[quadratic-gauge stable-moduli theorem](../verified/QUADRATIC_GAUGE_STABLE_MODULI.md).
Its stable intersection with the full cancellation family is now exact:
the foundational cubic is the only common class.  In every degree at least
four, the ramified-stratum Fitting support and the thick boundary contact
give independent obstructions; see the
[quadratic-gauge/cancellation stable-intersection theorem](../verified/QUADRATIC_CANCELLATION_STABLE_INTERSECTION.md).
The corresponding marked-line criterion and the complete
root-preserving, `P`-fibration-preserving affine rechart search through
`deg X<=4` are in the
[degree-four incidence-suspension classification](../verified/INCIDENCE_SUSPENSION_DEGREE_FOUR_CLASSIFICATION.md);
the rational `X=S^3,S^4` charts have unavoidable source poles.

## 2. Finite cancellation operator

Every positive `z`-degree term in (2) is automatically polynomial.  The only
possible pole is the `z`-free part

\[
 R|_{z=0}=\frac{Cx}{A^{r+1}}\Phi_{m,r}(A,h(A)),             \tag{5}
\]

where

\[
 \Phi_{m,r}(A,H)=
 \int_0^1\left[A+(1-A)u\{1-(1-A)H(1-u)\}^{m}\right]^rdu.  \tag{6}
\]

Define

\[
 \mathcal L_{m,r}(h)=\Phi_{m,r}(A,h(A))\bmod A^{r+1}.      \tag{7}
\]

The integral only divides the coefficient of `u^d` by `d+1`, so this is a
finite algebraic operator over `Q`.

### Polynomiality theorem

For every `h in K[A]`,

\[
 R\in K[x,y,z]
 \quad\Longleftrightarrow\quad
 \mathcal L_{m,r}(h)=0.                                   \tag{8}
\]

Indeed, `A=1+xy^m` is prime and does not divide `x`, so (5) is regular exactly
when `A^{r+1}` divides `Phi`.  An independent formula useful for exact
implementations is

\[
 \frac{R|_{z=0}}{Cx}=
 \sum_{j=0}^r\sum_{k=0}^{mj}
 (-1)^j\binom rj\binom{mj}k
 \frac{j!k!}{(j+k+1)!}
 \frac{(A-1)^{j+k}h(A)^k}{A^{j+1}}.                       \tag{9}
\]

## 3. Parameter roots and the unique cancellation jet

Let `q=h(0)` and `n=mr`.  The constant term of (7), after monic
normalization, is

\[
 \boxed{M_{m,r}(q)=
 \sum_{j=0}^{n}(-1)^j\binom{n+r+1}{j}q^{n-j}}.             \tag{10}
\]

This polynomial is squarefree in characteristic zero.  For every root `q`
there is a unique normalized jet

\[
 h_q(A)=q+h_1A+\cdots+h_rA^r                               \tag{11}
\]

annihilated by (7), obtained recursively from

\[
 h_d=-\frac{[A^d]\Phi_{m,r}
  (A,q+h_1A+\cdots+h_{d-1}A^{d-1})}
 {\partial_q\Phi_{m,r}(0,q)},\qquad1\le d\le r.           \tag{12}
\]

Adding `A^{r+1}g(A)` does not change the cancellation operator and is removed
by a polynomial source shift in `z`.  Thus (10)--(12), rather than a table of
small maps, define all normalized branches.  Their irreducibility,
discriminants, and Galois groups are recorded only in
[ARITHMETIC.md](ARITHMETIC.md).

## 4. Exact generic reconstruction

For target coordinates `(P,Q,R)`, define

\[
 \Psi_{P,Q,R}(T)=
 C\int_0^T(1-t(Q-Pt)^m)^r\,dt-R.                           \tag{13}
\]

It is irreducible over `K(P,Q,R)`, separable, and has degree

\[
 N=r(m+1)+1.                                               \tag{14}
\]

If `s` is a root and `D=1-s(Q-Ps)^m` is nonzero, the unique point over it is

\[
 y=Q-sP,\qquad A=D^{-1},\qquad x=sD^{-1},
\]

\[
 z=PD^{r+2}-y^{m+1}D^{r+1}h(D^{-1}).                      \tag{15}
\]

The resultant of `Psi` and `D` is nonzero (specialize `P=Q=0`), so all `N`
generic roots satisfy `D!=0`.  Formula (15) therefore proves exact generic
degree `N`, not just an upper bound.

## 5. Uniform full collision

At the target `(P,Q,R)=(1,0,0)`, put `epsilon=(-1)^m`.  Then

\[
 \Psi(T)/C=T E_{m,r}(\epsilon T^{m+1}),
\qquad
 E_{m,r}(U)=\sum_{j=0}^r
 \frac{(-1)^j\binom rj}{j(m+1)+1}U^j.                    \tag{16}
\]

This degree-`N` polynomial has no repeated root.  A common root with its
derivative would require `epsilon T^{m+1}=1`, but

\[
 E_{m,r}(1)=\int_0^1(1-v^{m+1})^r\,dv\ne0.
\]

For every root `s`, let

\[
 D_s=1-\epsilon s^{m+1},\qquad y_s=-s,\qquad x_s=s/D_s,
\]

\[
 z_s=D_s^{r+2}-(-s)^{m+1}D_s^{r+1}h(D_s^{-1}).            \tag{17}
\]

The `N` distinct points (17) all map to `(1,0,0)`.  This is a symbolic
all-parameter collision, independent of case tables.

## 6. Relation to the cubic seed

For `(m,r,q)=(1,1,3)`, one has `h=3+9A`; the resulting map is a linear
reparametrization of the foundational Keller map.  Every noncubic member is stably distinguished from
a generic weighted seed by [the boundary distinction
theorem](BOUNDARY_GEOMETRY.md).  Rigidity under enlargements of the present
ansatz and target-fixed parameter rigidity are intentionally collected in
[RIGIDITY.md](RIGIDITY.md), not claimed as part of this construction theorem.

## 7. Scope and verification

The argument establishes the determinant (4), polynomiality criterion (8),
parameter formula and recurrence (10)--(12), exact reconstruction
(13)--(15), and collision (16)--(17) for every `m,r>=1`.

The public reproduction target is `make verify-master`; see
[REPRODUCE.md](../REPRODUCE.md).  Optional generators and exploratory
regressions are retained under `archive/tooling/`.

The detailed former master note and its component derivations are retained in
[the cancellation archive](../archive/cancellation-components/).
