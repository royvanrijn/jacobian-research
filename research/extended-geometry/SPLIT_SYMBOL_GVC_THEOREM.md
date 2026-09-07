# The split-symbol Generalized Vanishing theorem

## 1. Statement

Let \(k\) be a characteristic-zero field, let \(\Lambda\) be a homogeneous
constant-coefficient differential operator of order \(d\) in \(n\)
variables, and let \(P\) be an arbitrary polynomial.

> **Theorem 1.1 — split-symbol GVC.**  Suppose that, after scalar
> extension to an algebraic closure, the symbol of \(\Lambda\) is a
> product of linear forms.  If
> \[
>  \Lambda^m(P^m)=0\qquad(m\geq1),
>  \tag{1.1}
> \]
> then, for every fixed polynomial \(Q\),
> \[
>  \Lambda^m(QP^m)=0\qquad(m\gg0).
>  \tag{1.2}
> \]

There is no restriction on \(\deg P\).  Every homogeneous binary operator
symbol splits over an algebraic closure, so:

> **Corollary 1.2.**  Every homogeneous constant-coefficient operator in
> two variables satisfies GVC for every polynomial \(P\).

Thus the degree-raising class
\(\deg P>\operatorname{ord}\Lambda\) is closed for homogeneous binary
operators.  The remaining two-variable operator frontier is genuinely
nonhomogeneous.

The same argument gives arbitrary SIC multipliers on the balanced binary
rank-one cone. If \(F=A(\xi)P(z)\) and
\(g=\sum_\nu B_\nu(\xi)Q_\nu(z)\), then
\[
 \mathcal E_2(gF^m)
 =\sum_\nu B_\nu(\partial_z)
   \bigl(A(\partial_z)^m(Q_\nu P^m)\bigr).
\]
Each expression in parentheses is eventually zero by Corollary 1.2.
Hence the complete rank-one Segre cone is SIC-safe; in bidegree
\((4,4)\), every counterexample has coefficient-matrix rank at least two.

## 2. Translated complete polarization

The zero operator is trivial, so assume its symbol is nonzero.  After
scalar extension, write
\[
 \sigma_\Lambda(\xi)
 =c\prod_{i=1}^d\ell_i(\xi),
 \qquad
 \ell_i(\partial)=D_{v_i}.
 \tag{2.1}
\]
The nonzero scalar \(c\) is irrelevant to vanishing, so suppress it.
Introduce \(t_1,\ldots,t_d\) and put
\[
 Vt=\sum_{i=1}^d t_iv_i.
\]
Taylor expansion gives, for every polynomial \(R\),
\[
 [t_1^m\cdots t_d^m]R(z+Vt)
 =\frac{1}{(m!)^d}
   D_{v_1}^m\cdots D_{v_d}^mR(z).
 \tag{2.2}
\]

Define the Laurent polynomial with polynomial \(z\)-coefficients
\[
 H_z(t)=\frac{P(z+Vt)}{t_1\cdots t_d}.
 \tag{2.3}
\]
Unlike the top-degree specialization used in the first proof, retaining
the translation variable gives the exact identity
\[
 \boxed{
 \Lambda^m(P^m)(z)
 =(m!)^d\operatorname{CT}_t(H_z(t)^m).
 }
 \tag{2.4}
\]
This holds for every degree of \(P\); the right side may be a nonconstant
polynomial in \(z\).

Likewise,
\[
 \boxed{
 \Lambda^m(QP^m)(z)
 =(m!)^d\operatorname{CT}_t
 \left(Q(z+Vt)H_z(t)^m\right).
 }
 \tag{2.5}
\]

## 3. One generic support gives a uniform separator

Let \(S\subset\mathbb Z^d\) be the set of \(t\)-exponents whose
coefficients in \(H_z(t)\) are nonzero polynomials in \(z\).  This is a
finite set.  Because the ground field is infinite, there is a point
\(z^\circ\) over the algebraic closure at which none of these finitely
many coefficient polynomials vanishes.  Hence
\[
 \operatorname{Supp}(H_{z^\circ})=S.
 \tag{3.1}
\]

The premise (1.1) and (2.4) imply
\[
 \operatorname{CT}(H_{z^\circ}{}^m)=0
 \qquad(m\geq1).
 \tag{3.2}
\]
The Duistermaat--van der Kallen constant-term theorem therefore gives
\[
 0\notin\operatorname{Newt}(H_{z^\circ})
 =\operatorname{conv}(S).
 \tag{3.3}
\]
Strict rational separation supplies an integral weight
\(w\in\mathbb Z^d\), scaled so that
\[
 w\mathbin{\cdot}\beta\geq1
 \qquad(\beta\in S).
 \tag{3.4}
\]

This one weight works uniformly in \(z\).  Indeed,
\(\operatorname{Supp}(H_z)\subseteq S\) for every specialization \(z\),
so every monomial of \(H_z^m\) has \(w\)-weight at least \(m\).

## 4. A fixed multiplier cannot cross the gap

Let \(T\) be the finite set of \(t\)-exponents whose coefficients in
\(Q(z+Vt)\) are nonzero polynomials in \(z\), and put
\[
 K=-\min\bigl(0,\min_{\eta\in T}w\mathbin{\cdot}\eta\bigr).
 \tag{4.1}
\]
Every monomial of
\[
 Q(z+Vt)H_z(t)^m
\]
has \(w\)-weight at least \(m-K\).  For \(m>K\), this weight is positive,
so the product has no constant term.  Equation (2.5) proves (1.2).

The separator even gives the explicit eventual bound
\[
 m>K.
 \tag{4.2}
\]
It depends on the exposed support of \(P(z+Vt)\), the chosen integral
separator, and the translated support of \(Q\), but not on coefficients
after specialization.

All data lie in a finitely generated characteristic-zero subfield.  After
embedding it in \(\mathbb C\), applying the constant-term theorem, and
proving the polynomial identity, the conclusion descends to the original
field.

## 5. A structured nonhomogeneous extension

The same proof gives a useful factor-unit continuation.  Suppose
\[
 \Lambda=\Lambda_0\Gamma,
 \tag{5.1}
\]
where \(\Lambda_0\) is homogeneous with split symbol and \(\Gamma\) is a
constant-coefficient operator with nonzero order-zero term.

> **Corollary 5.1 — split factor times a differential unit.**  Every
> operator (5.1) satisfies GVC for arbitrary \(P\).

To prove this, note that \(\Gamma\) acts injectively on the polynomial
ring.  If \(R_e\) is the top homogeneous part of a nonzero polynomial
\(R\), then the top homogeneous part of \(\Gamma R\) is \(cR_e\), where
\(c\ne0\) is the constant term of \(\Gamma\).  Hence every
\(\Gamma^m\) is injective.  Since constant-coefficient operators commute,
\[
 0=\Lambda^m(P^m)
 =\Gamma^m\Lambda_0^m(P^m)
\]
implies \(\Lambda_0^m(P^m)=0\).  Theorem 1.1 gives
\(\Lambda_0^m(QP^m)=0\) eventually, and applying \(\Gamma^m\) proves the
claim.

For a nonhomogeneous binary operator, this covers exactly the natural
class in which every homogeneous piece is divisible by one fixed
homogeneous binary factor and the residual quotient has nonzero constant
term.  It allows arbitrary \(\deg P\).

## 6. Consequences and limits

The earlier exact results remain stronger in a different direction:

- the [cubic theorem](TWO_VARIABLE_CUBIC_GVC_THEOREM.md) replaces the
  infinite premise by the first four moments, identifies the Segre
  nullcone exactly, and gives the sharp simple cutoff \(m>\deg Q\);
- the [low-root calculations](TWO_VARIABLE_LOW_ROOT_GVC_THEOREMS.md) give
  the same explicit cutoff for all two-root balanced symbols and finite
  radical certificates on quartic strata.

The translated-support proof supplies an eventual bound but not a uniform
finite pure-moment certificate or a nullcone classification.  It does not
apply to a genuinely nonhomogeneous operator whose higher homogeneous
pieces are not all divisible by one split homogeneous factor.

The full-rank two-pair SIC witness is not contradicted: it is a
nonseparable form, not a product \(A(\zeta)P(z)\).  The
Dvorsky counterexample is also not contradicted: its five-variable cubic
symbol
\[
 \xi_t(\xi_a\xi_d-\xi_b\xi_c)
\]
contains a nonsplit rank-four quadratic factor.  Likewise, an
ordinary-Laplacian symbol in more than two variables is not a product of
linear forms.

No computer algebra is required for Theorem 1.1.  Its external input is
the Duistermaat--van der Kallen theorem:

- J. J. Duistermaat and W. van der Kallen,
  [*Constant terms in powers of a Laurent polynomial*](https://webspace.science.uu.nl/~kalle101/powers.pdf),
  Indagationes Mathematicae 9 (1998), 221--231.

For the surrounding Vanishing Conjecture context, see also
A. van den Essen, R. Willems, and W. Zhao,
[*Some Results on the Vanishing Conjecture of Differential Operators with Constant Coefficients*](https://arxiv.org/abs/0903.1478).
