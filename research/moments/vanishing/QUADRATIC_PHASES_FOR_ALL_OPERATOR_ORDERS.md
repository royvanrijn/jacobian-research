# Quadratic phases for all operator orders

15 November 2026. Internal written proof, independent review outstanding;
no novelty claim. This complements the [order-two operator theorem](ORDER_TWO_DRIFT_REDUCTION_AND_UNIFORM_SAFETY.md):
the operator order is now arbitrary, while the **polynomial phase** has degree
at most two. The [registry](../../MATH_STATUS.json) is authoritative.

## 1. Statement and source boundary

Let K have characteristic zero, let n>=1, and let
\[
 P(x)=\tfrac12x^THx+\ell^Tx+p_0,\qquad
 \lambda(\xi)=a_0+b^T\xi+\xi^TA\xi+\lambda_{\ge3}(\xi),
\]
with H,A symmetric and \(\Lambda=\lambda(\partial)\). The symbol may be a
polynomial of arbitrary order, or a formal power series: every such series
acts locally finitely on polynomials. P and all multipliers remain polynomials.

For P!=0, the following are equivalent:
\[
 \begin{split}
 &\Lambda^mP^m=0\quad(1\le m\le n+1);\\
 &a_0=0,\quad Hb=0,\quad b^T\ell=0,\quad(AH)^n=0;\\
 &\Lambda^mP^m=0\quad(m\ge1).
 \end{split}                                                    \tag{1}
\]
For P=0, purity and safety are automatic, including a_0!=0. Every pure pair
in (1) satisfies full GVC, with the common, nonoptimal bound
\[
 \boxed{\ \Lambda^m(QP^m)=0\quad\text{if }m>3(n+1)\deg Q.\ }    \tag{2}
\]
The bound is independent of operator order and every higher symbol coefficient.
For affine phases, slope2 suffices and is sharp in general; for constant phases,
slope1 suffices. The first n+1 test is sufficient, not claimed optimal.

For context, [van den Essen--Willems--Zhao0903.1478v2](https://arxiv.org/html/0903.1478v2)
proves the one-variable theorem and the separated binary drift theorem
(Theorems2.3,3.1). The retained version is 16 January 2011; its HTML body date
of August 2026 is not a new arXiv version. Selected complete proofs and the
formal-phase counterexample3.5 were read. This note supplies an all-dimensional
quadratic-phase proof; it does not attribute that theorem to the source.
Formal **operator symbols** here must not be confused with formal nonpolynomial
phases, for which the source gives failures. The only imported operator lemma
is [de Bondt's commutator bound](https://arxiv.org/html/1206.2836v2).

## 2. Pure nullity is determined by the first two symbol orders

If a_0!=0, Lambda is a unit on bounded-degree polynomials, so its kernel is zero.
Assume a_0=0. If b=0, every term of \(\Lambda^m\) involving a symbol piece
of order>=3 has derivative order>2m and kills P^m. Thus all pure moments equal
those for \(\Delta_A\), and the preceding theorem proves (1).

If b!=0, choose linear coordinates (z,y) with \(D_b=\partial_y\).
The first pure equation only uses the two-jet of lambda:
\[
 \partial_yP+\Delta_AP=0.
\]
Since P is quadratic, \(\Delta_AP\) is constant. Hence
\[
 P=\alpha y+p(z),\qquad\deg p\le2.                             \tag{3}
\]
Give each z coordinate weight1 and y weight2. The phase has weight at most2.
Write
\[
 \Lambda=B+R,\qquad B=\partial_y+\Delta_{A_0},                  \tag{4}
\]
where A_0 is A's z,z block. B has derivative weight2; every term of R has
derivative weight at least3. This includes the quadratic cross derivatives
involving y, not just the higher ordinary symbol orders. In \(\Lambda^mP^m\),
every term involving R therefore has derivative weight>2m and is zero.
The same argument applies to the ordinary order-two truncation of Lambda.

Consequently, **after the first pure equation holds**, all pure moments of
Lambda equal those of its two-jet. Apply the preceding order-two theorem to
the first n+1 tests. It forces \(\alpha=0\) and nilpotence of \(A_0H_z\).
In the original coordinates this is exactly Hb=0, b^Tell=0 and nilpotence of
AH: with P independent of y, AH has diagonal blocks A_0H_z and zero, so its
nilpotence is equivalent to that of the active block. Conversely these
conditions give purity for the order-two truncation, and the same weighted
argument gives every pure identity for Lambda. This proves (1).

The first-equation qualification matters. For
\(\Lambda=\partial_y+\partial_x^2\partial_y\) and P=xy,
\(\Lambda P=x\ne0\), and
\[
 \Lambda^2P^2=2x^2+8\ne2x^2=\partial_y^2P^2.
\]
Thus unrestricted two-jet equality of pure expressions would be false.

## 3. Bounding the higher symbol pieces on every multiplier

Use B from (4) if b!=0; if b=0 use B=Delta_A and ordinary weights1.
The pure criterion ensures B is pure on P. The scalar inverse potential of
this leading order-two operator is polynomial in t. More explicitly, for
its quadratic matrix A_B put J=A_BH and
\[
 S_t=(I-tJ)^{-1}=\sum_{j=0}^{n-1}t^jJ^j.
\]
The potential is
\[
 U_t=\tfrac12x^THS_tx+\ell^TS_tx+p_0
        +\tfrac t2\ell^TS_tA_B\ell.                            \tag{5}
\]
This is the full scalar potential, including its time-dependent constant.
It satisfies the Hamilton--Jacobi equation and the pure leading-operator
identity. Its time degree N is at most n. For B with drift, D_bP=0 and the
same potential applies. Thus
\[
 B^jP^{j+k}=0\quad\text{if }j>Nk,\quad k\ge1.                 \tag{6}
\]
No inverse-degree theorem or general integration producer is needed here.

Let d=deg Q>=1. In the binomial expansion
\[
 \Lambda^m(QP^m)=\sum_{r=0}^m\binom mr R^rB^{m-r}(QP^m),        \tag{7}
\]
terms with r>2d vanish by weight: Q has weight<=2d, P^m has weight<=2m,
and the derivative weight is at least2m+r. When b=0, r>d already vanishes,
but the common weaker bound suffices. For r<=2d put j=m-r-d. If
\(m>3(N+1)d\), then j>N(r+d). Equation (6) yields
\[
 B^{m-r-d}P^m=B^jP^{j+(r+d)}=0.
\]
The commutator lemma, with multiplier Q of degree d, gives
\(B^{m-r}(QP^m)=0\). Every remaining term of (7) is therefore zero.
Since N<=n, equation (2) follows. Constant Q is covered directly by purity.
For a formal symbol, each application to a fixed polynomial still has only
finitely many contributing derivatives, so the same argument applies.

This does **not** replace the whole higher-order heat logarithm by a scalar
inverse. The potential in (5) belongs to B, and the remainder is retained in
the finite multiplier estimate (7).

## 4. Sharp low-degree examples and the global phase minimum

For a nonconstant affine phase choose x=P up to a nonzero affine rescaling.
The first pure equation removes the derivative in its direction from the
linear symbol. Give x weight1 and all transverse coordinates weight2. Every
remaining symbol monomial has derivative weight>=2. Since Q has weight<=2d,
the marked expression vanishes for m>2d. The previously retained
\(\partial_x^2+\partial_y,P=x,Q=y^d\) attains every boundary m=2d.

Higher symbol terms really do change the quadratic mixed cutoff. For
\[
 \Lambda=\partial_y+\partial_x^3,\quad P=x^2,\quad Q=y^d,
\]
all pure identities hold, but
\[
 \Lambda^{3d}(y^dx^{6d})
    =\binom{3d}{d}d!(6d)!
    =\frac{(3d)!(6d)!}{(2d)!}\ne0.                             \tag{8}
\]
Weights x=1,y=3 prove that slope3 is sharp for this pair. Its two-jet is only
partial_y, which would kill the marked sequence after d. This is also a sharp
instance of the 3(N+1) bound with N=0 in (7). It is not a new counterexample:
the entire tail still vanishes eventually.

Because all phases of degree<=2 are safe for arbitrary operators, while DVG1
retains an independently replayed homogeneous cubic GVC failure, the least
phase degree of a GVC counterexample, over all dimensions and finite operator
orders, is **exactly three**. The prior result separately fixes minimum
operator order at three **among cubic failures**. Neither statement determines
the least operator order of arbitrary-degree failures or changes the ternary
Delta-power minimum.

## 5. Closed loci and an implemented finite certificate

For every finite coefficient family with deg P<=2, the pure field locus is
the closed union of P=0 and the equations in (1). Every point on it is safe,
with the common bound (2). Its dependence on the operator is a cylinder over
the two-jet: every higher symbol coefficient is free. This includes families
of unboundedly high declared operator order; no order enters the bound.
The assertion concerns field points and radicals, not full ordinary ideals
or arbitrary nonreduced algebras.

The [implemented classifier](verify_quadratic_phase_all_orders.py) extracts
the two-jet from a supplied finite constant-coefficient polynomial symbol and
emits these matrix equations. Its implemented coefficient domain is rational
data and polynomial coefficient parameters, with exact polynomial zero tests;
symbolic coefficient parameters retain the defining equations instead of
claiming a solved parameter decomposition. It also accepts a supplied
two-jet, certifying all formal higher-order extensions at once. This is a
finite quadratic-phase classifier, not a general higher-degree GVC solver.
The [receipt](quadratic_phase_all_orders_verification.json) and
[source audit](quadratic_phase_all_orders_source_audit.json) retain exact
versions, adversarial controls and the finite-check boundary.

Independent review of the weighted pure reduction and multiplier estimate is
outstanding. Phases of degree>=3 with higher-order operators, general Gaussian
uniformity, ordinary nilpotents, higher-degree quadraticization and the protected
numeric frontiers remain separate. Formal nonpolynomial phases and multipliers
are not covered by this theorem.
