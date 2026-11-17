# Ordinary-degree-six boundary audit

## 0. Result and status

This note records an exact audit of several proposed low-degree extensions of
the one-root compiler.  It supports the following deliberately scoped
statement over a characteristic-zero field.

> **Scoped degree-six exclusion.** There is no polynomial Keller map of
> ordinary degree at most six in the union of the following declared
> templates:
>
> 1. the normalized linear--quadratic factorization slice, with the linear
>    coefficient and product coordinates retained;
> 2. the balanced `2+2` quotient, its elementary Cox multiplication cover,
>    either standard resolution chart, or their plain two-chart gluing;
> 3. an affine-linear two-reconstruction-coordinate core with monomial
>    `kappa=S^2` coefficients, pulled back through the minimal reciprocal
>    chart of Section 3;
> 4. the `D^3` plane core lifted by the standard weighted source chart with
>    a target Jacobian supported only on `C=0`, or by the standard reciprocal
>    chart with any `z`-linear denominator clearing.

In the first template the maximum coordinate degree has exact minimum seven
in the retained gauge and under affine-linear left--right changes.  In the
third template the determinant forces a degree-seven term before the final
polynomiality test, and that test now excludes both monomial Wronskian
profiles for every polynomial boundary jet: `(1,2)` is impossible already
on the plane core, while `(0,3)` fails separately for zero, constant
nonzero, and nonconstant boundary jets.  Every `z`-linear standard
reciprocal `D^3` lift has degree floor eleven.

This is not a classification of all Keller maps of geometric degree three or
four.  It is also not a theorem that every presentation with at most two
boundary relations belongs to one of the four templates.  Ordinary degree is
not invariant under arbitrary nonlinear source--target equivalence, so the
degree floors below always retain the displayed coordinates and charts.

The repository-level search conclusion as of **2026-08-02** is therefore:

> No explicit noninjective Keller map
> `A^3 -> A^3` of ordinary degree at most six has been found in the
> repository or in the literature sources audited in Section 7.

This last sentence is a dated negative search report, not an exhaustive
theorem.

## 1. Normalized asymmetric linear--quadratic slice

Put

\[
 L=aT+bS,\qquad Q=cT^2+dTS+eS^2
\]

and impose

\[
 a^2e-abd+b^2c=1,\qquad ad+bc=1.                    \tag{1.1}
\]

The normalized factorization variety and its affine-three-space coordinates
are proved in
[the canonical factorization note](../verified/NORMALIZED_FACTORIZATION_MODEL.md).
Here we allow an arbitrary polynomial reparametrization that retains `a` and

\[
 b=1+ay.                                               \tag{1.2}
\]

The second equation in (1.1) gives

\[
 d={1-bc\over a},
\]

so polynomiality first forces `c=1+ah`.  Substitution in the formula

\[
 e={1+b-2b^2c\over a^2}
\]

then makes the coefficient of `a` in its numerator equal to `-2h-3y`
modulo `a`.  Therefore every polynomial parametrization in this gauge has

\[
 \boxed{c=1-\frac32ay+a^2g(a,y,z).}                  \tag{1.3}
\]

With `d,e` recovered from (1.1), let

\[
 p=ac,\qquad r=ae+bd,\qquad s=be.                    \tag{1.4}
\]

Direct exact differentiation gives

\[
 \boxed{\det {\partial(p,r,s)\over\partial(a,y,z)}=-g_z}       \tag{1.5}
\]

and

\[
 p_z=a^3g_z,\qquad
 r_z=-3a(1+ay)^2g_z,\qquad
 s_z=-2(1+ay)^3g_z.                                  \tag{1.6}
\]

The Keller condition is consequently

\[
 g=\lambda z+g_0(a,y),\qquad \lambda\in k^*.          \tag{1.7}
\]

The triangular source change `z -> lambda*z+g_0(a,y)` removes `g_0` and
normalizes the vertical slice.  More importantly for the degree bound, even
before that change the last coordinate contains

\[
 -2\lambda a^3y^3z,                                  \tag{1.8}
\]

which has ordinary degree seven and cannot be canceled by the
`z`-independent term `g_0`.  Taking `g=lambda*z` gives coordinate-degree
profile

\[
 (4,6,7).                                             \tag{1.9}
\]

Thus seven is the exact floor for the normalized product map in this
fiber-preserving gauge.  Maximum ordinary coordinate degree is preserved by
invertible affine source substitutions, and an invertible linear target
matrix cannot kill every nonzero top homogeneous component.  Hence the same
floor holds throughout its affine-linear left--right orbit.  It is not
asserted to be a minimum under arbitrary nonlinear polynomial left--right
changes.

## 2. Balanced `2+2` factorization

Write two ordered monic quadratics as

\[
 Q_\pm=T^2+(h\mathbin\pm u)T+(b\mathbin\pm v).        \tag{2.1}
\]

Factor swap is `(u,v) -> (-u,-v)`.  Its invariant ring is

\[
 k[h,b,A,B,C]/(B^2-AC),\qquad
 A=u^2,\quad B=uv,\quad C=v^2.                        \tag{2.2}
\]

After fixing `h`, as required by the normalized cubic coefficient, the
source is

\[
 \mathbb A^1_b\times V(B^2-AC).                       \tag{2.3}
\]

The cone is the affine toric surface with character lattice

\[
 M=\{(i,j)\in\mathbb Z^2:i+j\equiv0\pmod2\}.
\]

The toric divisor map embeds `M` in `Z^2` with index two, so

\[
 \operatorname{Cl}V(B^2-AC)\simeq\mathbb Z/2.         \tag{2.4}
\]

Polynomial extension by `b` preserves this class group.  The quotient
source is therefore not `A^3`.

The Cox/index-two cover uses the affine coordinates `(b,u,v)`.  The three
nonleading coefficients of `Q_+Q_-` are

\[
 X=2b+h^2-u^2,\qquad
 Y=2hb-2uv,\qquad
 Z=b^2-v^2.                                           \tag{2.5}
\]

Their exact determinant and the factor resultant are

\[
 \begin{aligned}
 \det {\partial(X,Y,Z)\over\partial(b,u,v)}
   &=8(bu^2-huv+v^2),\\
 \operatorname{Res}(Q_+,Q_-)
   &=4(bu^2-huv+v^2).                                \tag{2.6}
 \end{aligned}
\]

Thus the elementary Cox multiplication map is not Keller.  A generic
quartic has three unordered partitions of its four roots into two pairs;
ordering the two factors gives six.  The Cox cover therefore also changes
the generic factorization degree from three to six.  The degree change is
recorded as geometry, not used by itself as a Keller obstruction.

On the standard `A`-chart

\[
 A=x,\qquad B=xy,\qquad C=xy^2,                       \tag{2.7}
\]

the same product coefficients have determinant

\[
 \boxed{4x(b-hy+y^2).}                                \tag{2.8}
\]

It has two unavoidable affine divisors.  The complementary `C`-chart is

\[
 C=x',\qquad B=x'y',\qquad A=x'y'^2,                  \tag{2.9}
\]

and its determinant is

\[
 \boxed{-4x'(by'^2-hy'+1).}                           \tag{2.10}
\]

On the overlap,

\[
 x'=xy^2,\qquad y'=y^{-1}.                            \tag{2.11}
\]

Both determinants are, up to sign, four times the pullback of the same
resultant `bA-hB+C`.  Gluing the two charts gives `A^1_b` times the minimal
resolution of the quadratic cone.  Over every fixed value of `b`, its zero
section contains a closed complete exceptional `P^1`.  The threefold is
therefore not affine and hence is not `A^3`; the product map also retains the
resultant divisor on both charts.  Thus a plain two-chart resolution does
not solve the source or Keller conditions.  This excludes neither a
different affine modification nor conductor gluing that changes the source
and introduces a compensating target ledger.

## 3. Two affine-linear reconstruction coordinates

Let

\[
 \begin{aligned}
 B&=b_0(P,S)+b_1(P,S)Q,\\
 C&=c_0(P,S)+c_1(P,S)Q,
 \end{aligned}                                       \tag{3.1}
\]

and suppose the controlled divisor is

\[
 D=D_0(P,S)-\kappa(P,S)Q.                             \tag{3.2}
\]

The condition

\[
 \det {\partial(P,B,C)\over\partial(P,S,Q)}=\lambda D
\]

is equivalent coefficient by coefficient in `Q` to

\[
 \boxed{
 \begin{aligned}
 b_{1,S}c_1-b_1c_{1,S}&=-\lambda\kappa,\\
 b_{0,S}c_1-b_1c_{0,S}&= \lambda D_0.
 \end{aligned}}                                      \tag{3.3}
\]

These Wronskian identities explain the linear and quadratic reconstruction
coordinates already used in the repository.  For the first new monomial
case, take `kappa=S^2` and

\[
 b_1=\beta S^m,\qquad c_1=\gamma S^n,qquad
 \beta\gamma\ne0.                                    \tag{3.4}
\]

The first equation of (3.3) forces

\[
 m+n=3,\qquad m\ne n.                                \tag{3.5}
\]

Up to swapping `B,C`, the only nonnegative exponent types are

\[
 \{m,n\}=\{0,3\}\quad\hbox{or}\quad\{1,2\}.          \tag{3.6}
\]

This is an exhaustion only for the monomial coefficient hypothesis (3.4),
not for arbitrary polynomial solutions of the Wronskian equations.

The two profiles behave differently.  For `(m,n)=(1,2)`, the left side of
the second equation in (3.3) is divisible by `S`, whereas in the minimal
chart

\[
 D_0=1+PS^3.                                         \tag{3.7}
\]

Its constant term is nonzero, so this profile has no polynomial `b_0,c_0`
at the plane-core level.  For `(m,n)=(0,3)`, constant rescaling normalizes
`beta=gamma=1` and `lambda=3`; the second equation integrates to

\[
 c_0=S^3b_0-3\int S^2b_0\,dS-3S-\frac34PS^4+f(P).    \tag{3.8}
\]

The minimal reciprocal chart for this case is

\[
 t=1+x^2y,\qquad
 P=tq,\qquad S={x\over t},\qquad Q=ty+xq.             \tag{3.9}
\]

It satisfies

\[
 D=1-S^2Q+PS^3={1\over t},\qquad
 \det {\partial(P,S,Q)\over\partial(x,y,q)}=1.        \tag{3.10}
\]

Allow the complete polynomial `z`-linear clearing

\[
 q=q_0(x,y)+\eta(x,y)z.                               \tag{3.11}
\]

Thus \(q_0,\eta\in k[x,y]\).  This is also forced if the displayed `P` and
`Q` are required to be polynomial: then both `t*q_0` and `x*q_0` are
polynomial, and `gcd(t,x)=1`.

The source-chart determinant is exactly `eta`.  Combining it with the
relative determinant `lambda*D` shows that the Keller condition is

\[
 \eta=\mu t,\qquad \mu\in k^*.                        \tag{3.12}
\]

The retained target coordinate therefore is

\[
 P=\mu t^2z+tq_0(x,y),                                \tag{3.13}
\]

and contains `x^4y^2z`, of ordinary degree seven.  Consequently every
polynomial Keller realization in this precise chart has degree at least
seven.

There is a further exact polynomiality gate for the remaining `(0,3)`
profile.  If `t` divides `q_0`, a polynomial triangular change of `z`
removes `q_0`; set `q=t z`.  Since `Q` is then polynomial, polynomiality of
`B=b_0+Q` forces every monomial `P^iS^j` of `b_0` to satisfy `j<=2i`.
For such a monomial the `b_0`-dependent part of (3.8) is

\[
 S^3P^iS^j-3\int S^2P^iS^j\,dS
 ={j\over j+3}P^iS^{j+3}.                            \tag{3.14}
\]

Terms with `j=2i>0` create an uncancellable pole of order three, so they
must vanish.  The remaining order-two corrections have `j=2i-1`, hence
positive `z`-degree.  But the fixed order-two residue is

\[
 \left.t^2\left(-3S-\frac34PS^4+S^3Q\right)
 \right|_{t=0}
 =-x+\frac14x^4z.                                    \tag{3.15}
\]

Its nonzero `z^0` term cannot be canceled.  Thus the removable-jet `(0,3)`
subchart does not polynomialize.

It remains to inspect `q_0 mod t != 0`; this case also has an exact
elimination.  In the boundary quotient

\[
 k[x,y]/(t)\simeq k[x,x^{-1}],\qquad y=-x^{-2},
\]

put

\[
 \rho(x)=q_0(x,-x^{-2}),\qquad U=x\rho(x).            \tag{3.16}
\]

Write `H=b_0+PS` and introduce the coefficient operator

\[
 \mathcal T(H)=S^3H-3\int S^2H\,dS,                  \tag{3.17}
\]

where the integral has zero `S`-constant term.  Since `Q=PS+ty`, equations
(3.8) give the exact pullback identities

\[
 B=H+ty,
 \qquad
 C=\mathcal T(H)+{x^3y\over t^2}-{3x\over t}+f(P).   \tag{3.18}
\]

First suppose that `U` is nonconstant.  Every polynomial `H` in `k[P,S]`
has a finite diagonal expansion

\[
 H=\sum_{r\in\mathbb Z}P^r h_r(PS).
\]

Because `P` has boundary order one and `PS` has the nonconstant residue
`U(x)`, polynomiality of `B` forces every term with `r<0` to vanish: the
lowest such term would have nonzero residue `h_r(U(x))`.  Thus `r>=0`.
For a polynomial `h`, exact monomial integration gives

\[
 \mathcal T(P^r h(U))=P^{r-3}\mathcal L(h)(U),
 \qquad
 \mathcal L(U^j)={j\over j+3}U^{j+3}.                \tag{3.19}
\]

The triple-pole coefficient forces `h_0` to be constant, since constants
are exactly the kernel of \(\mathcal L\).  The double-pole coefficient would
then require

\[
 U(x)R(U(x))=x^{-1},                                 \tag{3.20}
\]

where `R` is a polynomial with zero constant term.  If it is nonzero, the
polynomial `uR(u)` has degree at least two.  Degrees of nonconstant rational
maps of `P^1` multiply under composition, whereas `x^{-1}` has degree one.
Equation (3.20) is impossible.

It remains only that `U=b` is a nonzero constant.  Then
`q_0+bxy` is divisible by `t`, so a polynomial triangular change normalizes

\[
 q=-bxy+\mu tz,qquad
 PS=b+t(\mu xz-b),qquad
 {P\over t}\bigg|_{t=0}={b\over x}.                 \tag{3.21}
\]

Put `delta=(PS-b)/P`; its boundary value is

\[
 \delta_0={\mu x^2z\over b}-x.
\]

The Laurent polynomial \(\mathcal T(H)\) lies in `k[P,P^{-1},PS]`.
Expanding at \(PS=b+P\delta\), and using that (3.18) permits no pole worse
than two, shows that its double-pole coefficient is
\(P^{-2}K(\delta)\) for a polynomial `K`.  Consequently

\[
 \left.t^2\mathcal T(H)\right|_{t=0}
 ={x^2\over b^2}K\left({\mu x^2z\over b}-x\right).   \tag{3.22}
\]

The fixed double-pole residue in (3.18) is `-x`.  If `K` is nonconstant,
(3.22) depends on `z`; if `K` is constant, it is a scalar multiple of
`x^2`.  Neither can cancel `-x`.  This excludes the constant nonzero jet.
Together with (3.15) and the nonconstant argument, it exhausts every
polynomial `q_0`.  Hence neither monomial profile in (3.6) produces a
polynomial Keller map in the standard chart.  The forced degree-seven term
in (3.13) remains useful as a lower bound for broader constructions that
introduce an additional pole-canceling target ledger.

## 4. The total-ramification `D^3` core

The linear-section quartic core can be normalized to

\[
 D=q-w,\qquad
 \Phi(w,q)=\left(q,{D^4-q^4\over4}\right).            \tag{4.1}
\]

Then

\[
 \det D\Phi=D^3,                                     \tag{4.2}
\]

and the inverse equation is quartic in `D`.  This verifies the proposed
plane skeleton; the issue is its threefold lift.

### 4.1 Standard weighted pure-`C` ledger

Retain the weighted source chart used by the known family.  Its Jacobian is,
up to a constant,

\[
 J_\alpha=x^3D^2,\qquad C=xD.                         \tag{4.3}
\]

Suppose more generally that the target-chart Jacobian is supported only on
`C=0`.  Since the target coordinate ring is a UFD, it is a unit times
`C^M` for some `M>=0`.  The elementary chart

\[
 \beta_{m,n}(A,B,C)=(BC^m,AC^n,C),\qquad m,n\ge0,
\]

realizes `M=m+n`, but the following ledger depends only on the divisor.
Combining (4.2) and (4.3), the composite determinant is a unit times

\[
 {x^3D^5\over(xD)^M}=x^{3-M}D^{5-M}.                 \tag{4.4}
\]

Constancy would require both `M=3` and `M=5`.  Hence no target chart whose
Jacobian divisor is supported only on `C=0` can lift this core through the
standard weighted source chart.  A target ledger with another prime remains
outside this argument.

### 4.2 Standard reciprocal lift

Put

\[
 t=1+xy,\qquad
 P=tq,\qquad S={x\over t},\qquad Q=y+xq.              \tag{4.5}
\]

Then

\[
 D=1-SQ+PS^2={1\over t},\qquad
 \det {\partial(P,S,Q)\over\partial(x,y,q)}={1\over t}.
                                                               \tag{4.6}
\]

For the complete `z`-linear clearing `q=q_0(x,y)+eta(x,y)z`, the source
determinant is `eta/t`.  Multiplication by the core determinant
`D^3=t^-3` is constant exactly when

\[
 \eta=\mu t^4,\qquad \mu\in k^*.
\]

Thus

\[
 P=\mu t^5z+tq_0(x,y)                                \tag{4.7}
\]

contains `mu*x^5y^5z`, of ordinary degree eleven.  This is an exact floor
for every `z`-linear clearing in the standard reciprocal chart, not an
obstruction to a
nonmultiplicative lift or a distributed target ledger.

## 5. Nodal conductor gluing

For the elementary node algebra

\[
 A=k+t(t-1)k[t]\subset B=k[t],\qquad I=t(t-1)B,
\]

the conductor quotients are

\[
 A/I\simeq k,\qquad B/I\simeq k\times k,              \tag{5.1}
\]

with the first unit group embedded diagonally in the second.  The
normalization/conductor exact sequence therefore contains the quotient

\[
 {(k^*)^2\over k^*_{\rm diag}}\simeq\mathbb G_m.       \tag{5.2}
\]

This is the exact extra gluing character.  The existing
[one-chart obstruction](CONDUCTOR_FIRST_ONE_CHART_OBSTRUCTION.md) excludes
the separated nodal chart, and the
[symmetric three-boundary Cox-fill obstruction](CONDUCTOR_THREE_BOUNDARY_COX_FILL_OBSTRUCTION.md)
shows that the minimal symmetric fill is still not affine three-space.

Equation (5.2) does **not** prove that every nodal construction needs
exactly three primitive boundary relations.  A different asymmetric or
nonprincipal ledger could mix the unit and divisor classes.  “A third
boundary relation is needed” is therefore retained only as the design
heuristic suggested by the known two-relation compiler, not as a theorem.

## 6. Boundary-signature census

| signature | degree-at-most-six status | logical strength |
|---|---|---|
| geometric degree 3, one selected root, normalized `(1,2)` product slice | excluded; exact retained-gauge floor 7 | proved in Section 1 |
| geometric degree 3, two selected roots, raw `2+2` quotient | fixed-`h` source has class group `Z/2` | proved in Section 2 |
| same, elementary Cox cover | affine source, but multiplication Jacobian is twice the resultant; ordered degree is 6 | proved in Section 2 |
| same, standard two-chart resolution | both charts retain the resultant divisor; their gluing contains an exceptional `P^1` and is not affine | proved in Section 2; no claim about other modifications |
| geometric degree 4, weighted `A^1` core | known coordinates already exceed degree 6 | existing explicit-family calculation |
| geometric degree 4, cancellation/`G_m` core | known coordinates already exceed degree 6 | existing explicit-family calculation |
| geometric degree 4, monomial cubic reconstruction in (3.9) | `(1,2)` impossible; `(0,3)` fails for every polynomial boundary jet | exact elimination, Section 3; broader target ledgers open |
| geometric degree 4, total ramification `D^3` | every standard weighted pure-`C` ledger fails; every standard reciprocal `z`-linear lift has floor 11 | proved in Section 4 |
| nodal conductor gluing | separated one-chart and minimal symmetric fill excluded | existing theorems; general asymmetric gluing open |
| direct two-center three-puncture reciprocal core | excluded in every degree | [`TPR1`](../cancellation/PUNCTURE_RANK_FRONTIER.md) |
| nonlinear double-incidence three-puncture completion | bounded affine/quadratic screens only | [still open](THREE_PUNCTURE_NONLINEAR_COMPLETION_FRONTIER.md) |
| asymmetric Cox filling with a nontrivial target divisor ledger | no construction or exhaustion theorem | open |

The table does not prove that all rational two-boundary signatures have been
enumerated.  The phrase “the search has collapsed to three boundaries” is a
useful prioritization statement, not a completed small-boundary
classification.

## 7. Literature audit

The search through 2026-08-02 did not locate a displayed noninjective Keller
map `A^3 -> A^3` of maximum ordinary coordinate degree four, five, or six.
This is not an exhaustive bibliographic theorem.

Two recent sources need their degree notions separated carefully.

1. Jelonek's
   [*On mappings with Jacobian one*](https://arxiv.org/abs/2607.20597)
   defines `X(n,d)` using maximum ordinary coordinate degree at most `d`.
   Corollary 2.3 states that, if `X(n,d)` is irreducible and `n>=3,d>=6`,
   its generic point is a counterexample.  Its proof only cites an external
   X post; it displays no degree-six map.  The repository's foundational
   witness has profile `(7,6,4)`.  The corollary should therefore not be
   quoted, by itself, as an explicit degree-six existence result.
2. Migus's
   [*Generic degrees of real polynomial Keller maps with non-dense image*](https://arxiv.org/abs/2607.21572v2)
   uses `d` for **geometric generic degree** in the family `G_d`.  The paper
   explicitly gives ordinary coordinate-degree profile
   `(7,6d+2,6d)`.  Its realization of generic degree four is consequently
   not an ordinary-degree-four construction.  The same paper records
   ordinary component degrees `(12,11,4)` for Gallagher's separate
   generic-degree-four map, so that example also lies outside the
   ordinary-degree-at-most-six search.

The graded and homogeneous results reconciled in
[the mixed-sign literature note](MIXED_SIGN_GRADINGS_LITERATURE_RECONCILIATION.md)
are consistent with the census, but they do not supply an exhaustion of
mixed-homogeneous maps without a chosen grading.

## 8. Remaining high-value searches

The exact calculations prioritize three related programs without claiming
that they are exhaustive.

1. **Nodal `2+2` conductor gluing.** The plain complementary-chart
   resolution is now excluded.  A survivor must change that gluing and use
   an asymmetric source--target ledger capable of absorbing both the
   index-two class and the conductor gluing character.  The first task is to
   prove the complete local degree sums, not to assume a three-prime count.
2. **Cubic reconstruction with a boundary ledger.** Both monomial profiles
   are now excluded when `q_0` is polynomial and no extra ledger is present.
   The next cases are nonmonomial solutions of the Wronskian equations, or
   a controlled reciprocal pole canceled by a separate target prime.
3. **Three-puncture Laurent cores.** Enumerate primitive valuation matrices
   with determinant `+/-1`, enforce all local degree sums `sum(e*f)=3` or
   `4`, and only then solve the determinant and polynomiality equations.

Thus coefficient reduction of the known degree-seven representative is no
longer the best-supported route inside these templates.  The live geometric
frontier is an asymmetric multi-chart or multi-prime boundary architecture.

## 9. Exact reproduction

Run

```bash
.venv/bin/python scripts/verify_ordinary_degree_six_boundary_audit.py
```

The checker verifies (1.1)--(1.9), the cone lattice index, both resolution
charts and their transition, the Wronskian determinant and exponent census,
the `(1,2)` divisibility obstruction, all three boundary-jet branches of
the `(0,3)` polynomiality obstruction, both general `z`-linear reciprocal
ledgers and degree floors, the pure-`C` weighted incompatibility, and the
rank-one nodal conductor character.  It performs no unbounded coefficient
search and no literature search.
