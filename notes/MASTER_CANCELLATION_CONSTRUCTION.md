# The master cancellation construction

This note treats the maps below as one construction.  The small cases are
specializations and regressions, not independent primitive discoveries.

Throughout, `k` is a characteristic-zero field, `m,r>=1`, `C in k^*`, and
`h in K[A]` for a finite separable extension `K/k`.  Put

\[
 A=1+xy^m,\qquad
 B=A^{r+1}z+y^{m+1}h(A),\qquad
 P=AB,\qquad Q=y+xB,                                      \tag{1}
\]

and, initially in `K[x,y,z,A^{-1}]`, put

\[
 R=C\int_0^{x/A}\bigl(1-t(Q-Pt)^m\bigr)^r\,dt.             \tag{2}
\]

The integral is the formal polynomial antiderivative with zero constant term.
The two logically separate questions are:

1. why `(P,Q,R)` always has constant Jacobian in the localization; and
2. for which `h` the rational expression (2) belongs to `K[x,y,z]`.

The first question has a uniform answer independent of `h`.  The second is a
finite cancellation problem.

## 1. The localized Jacobian theorem

Set

\[
 s={x\over A},\qquad y=Q-sP,\qquad
 D=1-s(Q-sP)^m.
\]

The identities in (1) give

\[
 D=A^{-1},\qquad x={s\over D},\qquad B=PD,
\]

and hence

\[
 z=PD^{r+2}-(Q-sP)^{m+1}D^{r+1}h(D^{-1}).                 \tag{3}
\]

These are exact reconstruction formulas wherever `D!=0`.

### Theorem 1

In `K[x,y,z,A^{-1}]`,

\[
 \det {\partial(P,Q,R)\over\partial(x,y,z)}=-C.            \tag{4}
\]

### Proof

First replace `z` by `B`.  The determinant of
`(x,y,z) -> (x,y,B)` is `A^(r+1)`.  At fixed `y`, the identity
`x=s/(1-sy^m)` gives

\[
 \det {\partial(s,y,B)\over\partial(x,y,B)}=A^{-2}.
\]

At fixed `s`, use `P=AB` and `Q=y+sP`.  The two-by-two determinant in
the variables `(y,B)` is `-A`; the terms containing `partial P/partial y`
cancel.  Therefore

\[
 \det {\partial(s,P,Q)\over\partial(x,y,z)}=-A^r.          \tag{5}
\]

The cyclic reordering `(s,P,Q) -> (P,Q,s)` has positive sign.  Holding
`P,Q` fixed in (2) gives

\[
 {\partial R\over\partial s}
 =C\bigl(1-s(Q-Ps)^m\bigr)^r=CD^r=CA^{-r}.
\]

Multiplying this with (5) proves (4).  Notice that neither `h` nor its
derivatives survive.  This is the universal Jacobian argument; the generated
symbolic checks in `verify_master_universal.py` are regressions for it, not its
proof.  QED

## 2. The finite cancellation operator

After setting `t=xu/A`, expansion of (2) by the powers `j` of the outer
binomial and `k` of `B` gives terms proportional to

\[
 {x^{j+k+1}y^{mj-k}B^k\over A^{j+1}},
 \qquad 0\le j\le r,\quad 0\le k\le mj.                   \tag{6}
\]

If a term in `B^k` uses `ell>=1` copies of `A^(r+1)z`, its net power of
`A` is at least

\[
 (r+1)\ell-(j+1)\ge0.
\]

Thus every positive `z`-degree term is automatically polynomial.  Only the
`z`-free part can have a pole.  It is

\[
 R\big|_{z=0}={Cx\over A^{r+1}}\Phi_{m,r}(A,h(A)),          \tag{7}
\]

where the completely explicit polynomial is

\[
 \Phi_{m,r}(A,H)=
 \int_0^1\left[
 A+(1-A)u\{1-(1-A)H(1-u)\}^{m}
 \right]^rdu.                                             \tag{8}
\]

Define the finite cancellation operator

\[
 \boxed{\quad
 \mathcal L_{m,r}(h)
   =\Phi_{m,r}(A,h(A))\bmod A^{r+1}.
 \quad}                                                   \tag{9}
\]

The integral in (8) merely divides the coefficient of `u^d` by `d+1`, so
(9) is an explicit finite algebraic operator over `Q`.

### Theorem 2

For every `h in K[A]`,

\[
 R\in K[x,y,z]
 \quad\Longleftrightarrow\quad
 \mathcal L_{m,r}(h)=0.                                  \tag{10}
\]

### Proof

The preceding expansion proves that the positive `z`-degree part is regular
and gives (7) for the remaining part.  The polynomial `A=1+xy^m` is prime in
`K[x,y]` and does not divide `x`.  Hence the right side of (7) is regular if
and only if `A^(r+1)` divides `Phi`.  This is exactly (9).  QED

An equivalent double-sum formula, useful for independent implementations, is

\[
 {R|_{z=0}\over Cx}=
 \sum_{j=0}^r\sum_{k=0}^{mj}
 (-1)^j{r\choose j}{mj\choose k}
 {j!k!\over(j+k+1)!}
 {(A-1)^{j+k}h(A)^k\over A^{j+1}}.                        \tag{11}
\]

## 3. Uniform algebraic branches

Let `q=h(0)`.  The constant term of (9) is

\[
 \Phi_{m,r}(0,q)=
 \int_0^1u^r(1-q+qu)^{mr}\,du.                            \tag{12}
\]

Writing `n=mr`, its monic normalization is

\[
 \boxed{\quad
 \mathcal M_{m,r}(q)=
 \sum_{j=0}^{n}(-1)^j{n+r+1\choose j}q^{n-j}.
 \quad}                                                   \tag{13}
\]

Thus the apparent parameter tables are sections of one binomial expansion.
Equivalently,

\[
 \Phi_{m,r}(0,q)
 ={1\over r+1}\,{}_2F_1(-mr,1;r+2;q)                     \tag{14}
\]

and, after normalization, it is the Jacobi polynomial
`P_(mr)^(r+1,-mr-r-1)(1-2q)`.  It is therefore a terminating
hypergeometric/truncated-binomial object.  It is not generically a Laguerre
polynomial, and no Pade interpretation is needed for the cancellation.

The hypergeometric differential equation also shows that (13) is squarefree
in characteristic zero.  Indeed a common zero of the function and its
derivative away from `q=0,1` would force the solution of the second-order ODE
to vanish identically, while (12) is nonzero at both `0` and `1`.

For every root `q` of (13), there is a unique jet

\[
 h_q(A)=q+h_1A+\cdots+h_rA^r                              \tag{15}
\]

with `L_(m,r)(h_q)=0`.  It is obtained without a table by the recurrence

\[
 h_d=-{[A^d]\Phi_{m,r}
   (A,q+h_1A+\cdots+h_{d-1}A^{d-1})
   \over \partial_q\Phi_{m,r}(0,q)},
 \qquad 1\le d\le r.                                    \tag{16}
\]

All entries are evaluated in `Q(q)`.  Adding `A^(r+1)g(A)` for arbitrary
`g` does not change the cancellation operator.  Formula (16), rather than a
finite coefficient table, is the uniform coefficient specification.

The previously separated labels are now specializations:

- `(N_m)` is `r=1`, with parameter polynomial
  `N_m(q)=sum_(j=0)^m (-1)^j binom(m+2,j)q^(m-j)` and the one-step
  recurrence (16).
- `(H_r)` is `m=1`, with
  `H_r(q)=sum_(j=0)^r (-1)^j binom(2r+1,j)q^(r-j)` and the same recurrence.
- `(M_(m,r))` is the full formula (13)--(16); the first two branches are its
  coordinate axes in parameter space, not different constructions.

Four exact displayed regressions are

| label | parameter polynomial | cancellation polynomial `h(A)` |
|---|---|---|
| `N_1` | `q-3` | `q+9A` |
| `N_2` | `q^2-4q+6` | `q+(4q-6)A` |
| `H_2` | `q^2-5q+10` | `q+(5q+10)A/3+(25q+50)A^2/9` |
| `M_(2,2)` | `q^4-7q^3+21q^2-35q+35` | `q+(q^3/3-11q^2/6+35q/6-35/6)A+(7q^3/18-7q^2/4+245q/36-70/9)A^2` |

These rows are generated by `generate_master_regression.py`; they are not the
definition of the branches.

### Coefficient fields and conjugacy

For a chosen root `q`, the smallest coefficient field of the normalized map
is `K_q=Q(q)`: formula (16) puts every coefficient in that field, and the
constant coefficient of `h` recovers `q`.  Roots lying in the same irreducible
factor of (13) give Galois-conjugate maps.  The displayed fields are:

- `N_1`: `Q`;
- `N_2`: `Q(sqrt(-2))`, with its two maps conjugate;
- `H_2`: `Q(sqrt(-15))`, with its two maps conjugate;
- `M_(2,2)`: `Q(q)` for a root of the displayed irreducible quartic; the four
  maps are conjugate and its splitting-field group is `D_4`.

Uniform irreducibility and Galois group of (13) are **not** claimed.  The
reciprocal is a truncated binomial polynomial, whose general irreducibility is
a genuine arithmetic problem; see
[Filaseta--Kumchev--Pasechnik](https://arxiv.org/abs/math/0409523),
[Khanduja--Khassa--Laishram](https://arxiv.org/abs/1306.0758), and the more
recent [Laishram--Yadav](https://digitalcommons.isical.ac.in/journal-articles/4864/).
The exact scripts factor the displayed cases only.

## 4. Generic fiber degree and exact reconstruction

For target coordinates `(P,Q,R)`, define

\[
 \Psi_{P,Q,R}(T)=
 C\int_0^T(1-t(Q-Pt)^m)^r\,dt-R.                           \tag{17}
\]

Its degree in `T` is

\[
 N=r(m+1)+1,                                               \tag{18}
\]

because the leading coefficient is a nonzero rational multiple of `P^(mr)`.
Over `K(P,Q,R)`, (17) is irreducible: over `K(P,Q)` the polynomial
`R-C integral` is linear in the indeterminate `R`, defines the prime graph of
the nonconstant polynomial in `T`, and is primitive; Gauss's lemma then gives
irreducibility over the fraction field.  It is separable in characteristic
zero, also directly from

\[
 \Psi'(T)=C(1-T(Q-PT)^m)^r.                               \tag{19}
\]

If `s` is a root and `D=1-s(Q-Ps)^m!=0`, the unique source point above it is

\[
 y=Q-sP,\quad A=D^{-1},\quad x=sD^{-1},
\]

\[
 z=PD^{r+2}-y^{m+1}D^{r+1}h(D^{-1}).                     \tag{20}
\]

Conversely, substitution in (1)--(2) recovers the target.  The resultant of
`Psi` and `D` is not the zero polynomial (specialize `P=Q=0`, when `D=1`).
Therefore every one of the `N` generic roots satisfies `D!=0`, and (20)
reconstructs exactly `N` distinct source points.  This proves the generic
degree, not merely an upper bound.

## 5. A symbolic collision for every `(m,r)`

Take the target `(P,Q,R)=(1,0,0)` and put `epsilon=(-1)^m`.  Then

\[
 \Psi(T)/C=T E_{m,r}(\epsilon T^{m+1}),                   \tag{21}
\]

where

\[
 E_{m,r}(U)=\sum_{j=0}^r
 {(-1)^j\binom rj\over j(m+1)+1}U^j.                     \tag{22}
\]

Equation (21) has degree `N`.  Its only possible common roots with its
derivative would satisfy `epsilon T^(m+1)=1`, but

\[
 E_{m,r}(1)=\int_0^1(1-v^{m+1})^r\,dv\ne0.
\]

Thus it has exactly `N` distinct roots over an algebraic closure.  For every
root `s`, set

\[
 D_s=1-\epsilon s^{m+1},\quad y_s=-s,\quad x_s=s/D_s,
\]

\[
 z_s=D_s^{r+2}-(-s)^{m+1}D_s^{r+1}h(D_s^{-1}).            \tag{23}
\]

The nonvanishing just proved gives `D_s!=0`, and (20) proves symbolically that
all `N` distinct points (23) map to `(1,0,0)`.  This is the uniform collision
proof for arbitrary `m,r`; no case table is involved.

For a fixed map over `K_q`, the collision points are defined over the
compositum of `K_q` with the splitting field of (21).  The Galois group of
(21) permutes the points while fixing the map; conjugating `q` conjugates the
map and its entire marked collision.  No smaller collision field is asserted
without factoring (21) for the chosen pair.

## 6. Relation with the weighted-seed construction

The member `(m,r,q)=(1,1,3)` has `h=3+9A`.  If `F=(F_1,F_2,F_3)` is the
original cubic map in `FOUNDATIONAL_GEOMETRY.md`, direct substitution gives

\[
 (P,Q,R)(x,y,z)=(3F_1,F_2,F_3/2)(x,y,z/3).                \tag{24}
\]

So the first member is a linear reparametrization of the existing primitive
cubic map.

This does not make the whole master construction a disguised weighted-seed
pencil.  The primitive polynomial (17) has derivative (19).  For `r>=2` it is
generically a perfect `r`-th power and every critical point has ramification
index `r+1`.  A generic admissible weighted seed has inverse derivative
`H'(T)-P`, with simple critical points and ramification index two.  Ramification
indices survive source and target reparametrization, so the `r>=2` maps are
not reparametrizations of generic weighted seeds.

The precise conclusion is therefore:

- the master family strictly extends the original cubic construction;
- both constructions are pullbacks of the universal marked-root incidence;
- neither family is presently proved to contain the entire other weighted
  family; and
- classification up to polynomial source/target equivalence remains open.

## 7. Formal closure operations

If `F` and `G` are polynomial Keller maps, their Cartesian product, extension
by identity coordinates, and composition (when dimensions match) again have
constant nonzero Jacobian.  Their fiber degrees multiply in the generically
finite case, and a collision transports through these operations whenever the
other factor is held fixed.

These are formal closure corollaries only.  They do not define new primitive
`N`, `H`, or `M` branches and are not counted as equal-status discoveries in
the claim ledger.

## 8. Exact status and open classification target

Proved uniformly in this note:

- the localized determinant (4);
- equivalence of polynomiality and the finite operator (10);
- the parameter formula (13), the coefficient recurrence (16), and the
  resulting normalized simple-root jets;
- irreducibility, separability, exact degree, and reconstruction for the
  generic fiber polynomial (17);
- the `N`-point collision (21)--(23) for every `m,r`; and
- the comparison (24) and the ramification obstruction for `r>=2`.

Verified by finite exact computation:

- the four displayed maps and coefficient fields;
- generated Jacobian identities for `1<=m,r<=3`;
- the Hensel recurrence for the four displayed instances; and
- squarefreeness/factorization data in the bounded ranges named by the scripts.

Not proved uniformly:

- irreducibility or the Galois group of `M_(m,r)` over `Q`;
- the arithmetic minimal collision field for every pair;
- equivalence or inequivalence among arbitrary tail deformations
  `h_q+A^(r+1)g`; and
- a global classification, up to polynomial coordinate equivalence and field
  descent, of all polynomial cancellation branches represented by
  `L_(m,r)`.

That last item is the central open classification target.  The finite-jet
recurrence is the local input to such a theorem, not a substitute for its
arithmetic and equivalence statements.

## 9. Verification entry points

- `scripts/master_cancellation.py`: exact operator and recurrence library;
- `scripts/verify_master_universal.py`: universal formulas plus bounded
  structural regressions;
- `scripts/verify_master_instances.py`: exact polynomial maps, Jacobians, and
  collision polynomials for every displayed row; and
- `scripts/generate_master_regression.py m r`: generator for a new pair.

