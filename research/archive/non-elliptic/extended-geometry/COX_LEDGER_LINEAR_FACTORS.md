# Multi-boundary Cox ledgers for ordered linear factors

This note completes the Cox-ledger calculation for every number of ordered
linear binary factors.  It produces constant-Jacobian morphisms on explicit
Cox--Kummer boundary complements and separates two notions which coincide
for three factors but diverge afterwards:

1. cancelling the total Jacobian unit; and
2. retaining one primitive reconstruction coordinate for every independent
   boundary character.

The second construction uses exactly the unit rank.  The first uses only one
coordinate.  Consequently the proposed lower bound by unit rank is false for
unrestricted torsor-Keller morphisms; it becomes valid only after imposing a
separated-character reconstruction condition.

Work over a characteristic-zero field `k`.

## 1. Boundary lattice in every arity

Fix `s>=3`.  Let

\[
 L_i=u_iX+v_iY,\qquad
 P=\prod_{i=1}^sL_i=\sum_{j=0}^sp_jX^{s-j}Y^j,
\]

and write

\[
 r_{ij}=u_iv_j-v_iu_j.
\]

On `X_s=(P^1)^s`, remove every pairwise collision divisor `R_(ij)` and the
irreducible divisor

\[
 E=(m=0),\qquad
 m=p_0+p_s=\prod_i u_i+\prod_i v_i.                  \tag{1}
\]

Irreducibility follows from Gauss's lemma, exactly as in the three-factor
case: as a polynomial in `(u_1,v_1)`, the two coefficients are coprime
monomials.

The full boundary-class matrix has columns

\[
 [R_{ij}]=e_i+e_j,\qquad [E]=(1,\ldots,1).            \tag{2}
\]

The edge columns generate the even-coordinate-sum lattice in `Z^s`.
Therefore

\[
 \operatorname{coker}B_s\simeq
 \begin{cases}
 0,&s\text{ odd},\\
 \mathbb Z/2\mathbb Z,&s\text{ even},
 \end{cases}                                         \tag{3}
\]

because the `E` column has odd coordinate sum exactly when `s` is odd.
Divisor localization gives

\[
 \boxed{
 \operatorname{rank}\mathcal O(U_s)^*/k^*
 =\binom{s}{2}+1-s
 =\frac{(s-1)(s-2)}2=:r_s.
 }                                                    \tag{4}
\]

Thus the free unit rank is `r_s`, while the finite Smith obstruction has
order one or two according to the parity of `s`.

## 2. A minimal-index normalization tree

Choose a bipartite spanning tree `G` on `{1,...,s}`.  Let its bipartition
have sizes `alpha,beta`.  Form the square scaling matrix

\[
 M_G=\left(
 e_i+e_j\ (ij\in E(G))\ \middle|\ (1,\ldots,1)
 \right).                                            \tag{5}
\]

An elementary incidence-matrix calculation gives

\[
 |\det M_G|=|\alpha-\beta|.                           \tag{6}
\]

Choose `G` with

\[
 (\alpha,\beta)=
 \begin{cases}
 ((s-1)/2,(s+1)/2),&s\text{ odd},\\
 (s/2-1,s/2+1),&s\text{ even}.
 \end{cases}
\]

Such a tree can be taken to be the double star obtained by joining one
vertex in the first part to every vertex in the second and one vertex in the
second part to every remaining vertex in the first.  Then

\[
 d_s:=|\det M_G|=
 \begin{cases}
 1,&s\text{ odd},\\
 2,&s\text{ even}.
 \end{cases}                                         \tag{7}
\]

This realizes the minimum allowed by (3).

Define the normalized Cox--Kummer slice

\[
 Y_{s,G}=
 \left\{
 r_e=1\ (e\in E(G)),\quad m=1,\quad
 r_f\ne0\ (f\notin E(G))
 \right\}\subset\mathbb A^{2s}.                      \tag{8}
\]

The target squarefree chart is

\[
 T_s=
 \left\{
 m(P)=1,\quad\operatorname{Disc}(P)\ne0
 \right\}\simeq
 \mathbb A^s\setminus V(\operatorname{Disc}).         \tag{9}
\]

For odd `s`, the normalization is unimodular and `Y_(s,G)` is the
root-free Cox slice of the projective complement.  For even `s`, it is the
unavoidable `mu_2` Kummer slice.  In characteristic two the latter is not
etale, which is why the theorem is stated in characteristic zero.

## 3. Universal ambient determinant

Let

\[
 \Theta_G:\mathbb A^{2s}\longrightarrow\mathbb A^{2s}
\]

send the factor coefficients to the `s+1` coefficients of `P` followed by
the `s-1` resultants `r_e`, `e in E(G)`.

### Proposition 3.1

Up to the orientation determined by the order of the edges,

\[
 \boxed{
 \det D\Theta_G
 =\pm\det(M_G)
 \left(\prod_{1\le i<j\le s}r_{ij}\right)
 \left(\prod_{e\in E(G)}r_e\right).
 }                                                    \tag{10}
\]

### Proof

Work on the dense chart

\[
 L_i=a_i(X-x_iY).
\]

Then

\[
 r_{ij}=a_ia_j(x_i-x_j),\qquad
 P=A\prod_i(X-x_iY),\qquad A=\prod_i a_i.
\]

The coefficient differential contributes the Vandermonde
`prod_(i<j)(x_i-x_j)`.  Modulo the root differentials, each selected
resultant contributes

\[
 d\log r_{ij}=d\log a_i+d\log a_j.
\]

Together with `d log A=sum_i d log a_i`, these logarithmic scaling rows are
exactly `M_G`.  Taking their wedge gives `det(M_G)` and one extra factor
`r_e` for every selected edge.  Converting from `(a_i,x_i)` back to
`(u_i,v_i)` supplies the remaining powers of the `a_i`.  They combine with
the Vandermonde to give `prod_(i<j)r_(ij)`, proving (10) on a dense chart and
hence polynomially everywhere.  QED

Replace one coefficient output by `m`; this is a determinant-one linear
target change.  Taking the complete-intersection residue along (8), all
selected resultants and `m` become one.  Multiplication

\[
 \mu_{s,G}:Y_{s,G}\longrightarrow T_s,\qquad
 (L_1,\ldots,L_s)\longmapsto P                       \tag{11}
\]

therefore has residue Jacobian

\[
 \boxed{
 J_{\mu_{s,G}}=\pm d_s
 \prod_{f\notin E(G)}r_f.
 }                                                    \tag{12}
\]

There are exactly

\[
 \binom{s}{2}-(s-1)=r_s
\]

unselected edges.  Thus every free boundary direction appears once in the
Jacobian ledger.

Every geometric squarefree binary form has `s!` orderings of its projective
linear factors.  The scaling equations have kernel `mu_(d_s)`, so

\[
 \deg(\mu_{s,G})=d_s\,s!.                             \tag{13}
\]

Formula (12) shows directly that the map is etale in characteristic zero.

## 4. Separated-character suspension

Introduce one affine coordinate `z_f` for every nonedge of `G` and put

\[
 Z_f=\frac{z_f}{r_f}.
\]

Define

\[
 \widehat\mu_{s,G}:
 Y_{s,G}\times\mathbb A^{r_s}
 \longrightarrow
 T_s\times\mathbb A^{r_s},
\qquad
 (L_\bullet,(z_f))\longmapsto(P,(Z_f)).              \tag{14}
\]

Each quotient is regular because every `r_f` is a unit on `Y_(s,G)`.
The vertical derivative block is diagonal with determinant

\[
 \prod_{f\notin E(G)}r_f^{-1}.
\]

Combining this with (12) gives

\[
 \boxed{
 J_{\widehat\mu_{s,G}}=\pm d_s.
 }                                                    \tag{15}
\]

Scaling one target coordinate by `d_s^{-1}` gives determinant one.  The
geometric degree remains `d_s s!`, since each `z_f` is uniquely reconstructed
as `r_fZ_f`.

This is the rank-sharp separated Cox ledger:

\[
 \text{extra dimension}
 =r_s
 =\operatorname{rank}\mathcal O(U_s)^*/k^*.           \tag{16}
\]

For `s=3`, `r_s=1`, `d_s=1`, and (14) is the construction in
[the first three-factor ledger](COX_LEDGER_THREE_FACTOR.md).

## 5. The unrestricted lower bound is false

Let

\[
 \Delta_G=\prod_{f\notin E(G)}r_f.
\]

Instead of separating the boundary characters, introduce just one coordinate
and define

\[
 \mu^{\mathrm{compressed}}_{s,G}:
 Y_{s,G}\times\mathbb A^1_z
 \longrightarrow
 T_s\times\mathbb A^1_Z,\qquad
 (L_\bullet,z)\longmapsto
 \left(P,\frac z{\Delta_G}\right).                   \tag{17}
\]

Its Jacobian is again `+-d_s`.  For every `s>=4`,

\[
 1<r_s,
\]

so (17) is an explicit counterexample to

\[
 \text{extra suspension dimension}\ge
 \operatorname{rank}\mathcal O(U_s)^*/k^*
\]

if “suspension” merely means cancellation of the total determinant unit.

There is an even more basic dependence on conventions.  On a variety with
units, replacing the chosen source volume form `omega` by
`Delta_G omega` makes the unsuspended map (11) have constant Jacobian.
Unlike the affine-space condition, “constant” is not intrinsic until a Cox
or residue volume form has been fixed.

The quantitative conjecture therefore has a correct formulation only after
adding both requirements:

1. use the canonical residue/Cox volume forms; and
2. require the reconstruction coordinates to detect the independent
   boundary characters separately.

Under those requirements, (14) attains the predicted number exactly.

## 6. Why the construction does not extend across collisions

At the generic point of an unselected collision divisor `R_f`, let `t=r_f`
be its DVR parameter.  The primitive target function is `z_f/t`.  In the
product local ring

\[
 A[z_f],\qquad A\text{ a DVR with uniformizer }t,
\]

it has valuation `-1` at the prime `(t)`.  Hence it cannot extend regularly
across the generic collision point.

The only elementary repair is to impose `z_f=tZ_f`, namely to pass to an
affine-modification chart centered at `(t,z_f)`.  That changes the source
chart and its canonical form; it is precisely the new mixed-jet problem
which the direct multiplicative ansatz does not solve.

Thus (14) has `r_s` independently visible boundary valuations, but it is not
a polynomial Keller map on affine space and does not yet turn those
valuations into dicritical divisors of such a map.

## 7. Finite fields

The weighted family's universal `C=0` excess is absent: the target chart is
`m=1`, and there is no distinguished `C=0` plane.

The symmetric ordered-factor core supplies a different obstruction.  Over a
good finite field containing a split squarefree target:

- for odd `s`, that target has `s!` normalized ordered lifts;
- for even `s` and odd characteristic, it has `2s!` lifts whenever the
  `mu_2` normalization is rational.

The primitive affine coordinates are uniquely reconstructed and do not
change these counts.  Consequently neither (14) nor (17) is a permutation
on such a field.  Avoiding the weighted `C=0` excess is therefore necessary
but not sufficient; a permutation-oriented construction must also break or
twist the symmetric ordered-factor cover.

## 8. Remaining affine-space problem

The Cox/Kummer part of the linear-factor programme is now complete:

1. the full free and torsion boundary lattice is known for every `s`;
2. a minimal-index normalization is explicit;
3. the ambient determinant is factored for every normalization tree;
4. both separated and compressed constant-Jacobian suspensions are explicit;
5. the separated dimension equals the unit rank; and
6. the direct extension and symmetric finite-field obstructions are exact.

What remains is materially different: construct an affine modification which
extends the primitive quotients across the collision divisors while retaining
a constant polynomial Jacobian, or prove that a specified class of such
modifications cannot do so.  Nothing here claims that this final
affine-space problem is solved.

The same frontier appears independently for the exceptional Davenport cover.
Its [complete Cox-boundary audit](DAVENPORT_COX_BOUNDARY_OBSTRUCTION.md)
proves that stable straightening and every coordinate-preserving polynomial
suspension fail, leaving a nontrivial source/target affine modification as
the common unresolved step.

Allowing the non-tree collision units to vanish and adjoining their product
as an oriented discriminant coordinate gives the
[all-arity polynomial Cox maps](ORIENTED_LINEAR_FACTOR_COX_MAPS.md).  They
have constant residue Jacobian and at least `s-1` distinct dicritical
boundary primes.

## 9. Reproduction

Run

```bash
.venv/bin/python scripts/verify_linear_factor_cox_ledgers.py
```

The checker verifies the Smith parity, minimal-index trees, symbolic
determinant factorizations for the first arities, exact random-integer
specializations through `s=8`, suspension cancellation, and the degree
formula.
