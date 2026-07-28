# The split-symbol Generalized Vanishing theorem

## 1. Statement

Let \(k\) be a characteristic-zero field, let
\(\Lambda\) be a homogeneous constant-coefficient differential operator of
order \(d\) in \(n\) variables, and let \(P\) be any polynomial of degree
at most \(d\).

> **Theorem 1.1 — split-symbol GVC.** Suppose that, after scalar extension
> to an algebraic closure, the symbol of \(\Lambda\) is a product of linear
> forms.  If
> \[
>  \Lambda^m(P^m)=0\qquad(m\geq1),
>  \tag{1.1}
> \]
> then, for every fixed polynomial \(Q\),
> \[
>  \Lambda^m(QP^m)=0\qquad(m\gg0).
>  \tag{1.2}
> \]

Every homogeneous binary operator symbol splits over an algebraic closure.
Therefore:

> **Corollary 1.2.** If \(\Lambda\) is homogeneous of order \(d\) in two
> variables, then GVC holds for every \(P\) of degree at most \(d\).

This includes nonhomogeneous \(P\), and proves the GVC conclusion on the
complete rank-one Segre cone in every balanced two-pair bidegree.  It does
not assert that moment-zero Segre points lie in the
\(\mathrm{SL}_2\)-nullcone, and it does not prove unrestricted SIC(2),
nonhomogeneous operators, the case \(\deg P>d\), or ordinary-Laplacian GVC
in higher rank.

## 2. Complete polarization of the operator factors

The zero operator is trivial, so assume its symbol is nonzero.  After
scalar extension, write

\[
 \sigma_\Lambda(\xi)
 =c\prod_{i=1}^d\ell_i(\xi),
 \qquad
 \ell_i(\partial)=D_{v_i},
 \tag{2.1}
\]

where \(D_{v_i}\) is the directional derivative along \(v_i\).  The
nonzero scalar \(c\) is irrelevant to vanishing, so suppress it.

Introduce polarization variables \(t_1,\ldots,t_d\) and put

\[
 Vt=\sum_{i=1}^d t_i v_i.
\]

Taylor expansion gives, for every polynomial \(R\),

\[
 [t_1^m\cdots t_d^m]R(z+Vt)
 =\frac{1}{(m!)^d}
   D_{v_1}^m\cdots D_{v_d}^mR(z).
 \tag{2.2}
\]

Let \(P_d\) be the degree-\(d\) homogeneous part of \(P\), with \(P_d=0\)
allowed.  Since \(\Lambda^m\) has order \(dm\), every lower-degree term of
\(P^m\) is killed and

\[
 \Lambda^m(P^m)=\Lambda^m(P_d^m).
\tag{2.3}
\]

Apply (2.2) to \(R=P_d^m\).  The coefficient of total \(t\)-degree \(dm\)
cannot contain \(z\).  Hence

\[
\begin{aligned}
 \Lambda^m(P^m)
 &=(m!)^d[t_1^m\cdots t_d^m]P_d(Vt)^m\\
 &=(m!)^d\operatorname{CT}
   \left(\frac{P_d(Vt)}{t_1\cdots t_d}\right)^m.
\end{aligned}
\tag{2.4}
\]

This identity is the bridge from a split higher-order operator to one
Laurent polynomial.

## 3. The degenerate restriction

Let \(W=\operatorname{span}\{v_1,\ldots,v_d\}\).  If

\[
 B_0(t):=P_d(Vt)
\]

is zero, then \(P_d\) vanishes identically on \(W\).  In coordinates
adapted to \(W\), every monomial of \(P_d\) contains at least one transverse
variable, while every lower homogeneous part already has degree at most
\(d-1\).  Thus \(P\) has \(W\)-degree at most \(d-1\).

All \(dm\) derivatives in \(\Lambda^m\) lie in \(W\), whereas

\[
 \deg_W(QP^m)\leq\deg Q+(d-1)m<dm
\]

for \(m>\deg Q\).  Thus (1.2) holds in this case.  We may henceforth assume
\(B_0\ne0\).

## 4. The Newton-weight separator

Define the Laurent polynomial

\[
 H(t)=\frac{B_0(t)}{t_1\cdots t_d}.
\tag{4.1}
\]

By (1.1) and (2.4), every positive power of \(H\) has zero constant term.
The Duistermaat--van der Kallen constant-term theorem implies

\[
 0\notin\operatorname{Newt}(H).
\tag{4.2}
\]

Strict rational separation of the finite lattice polytope in (4.2)
supplies an integral weight

\[
 w=(w_1,\ldots,w_d)\in\mathbb Z^d
\]

such that

\[
 w\mathbin{\cdot}\beta\geq1
 \qquad
 \text{for every }\beta\in\operatorname{Supp}(H).
\tag{4.3}
\]

Put \(W_0=w_1+\cdots+w_d\).  If \(t^\gamma\) occurs in \(B_0\), then
\(\gamma-(1,\ldots,1)\) occurs in \(H\), so

\[
 w\mathbin{\cdot}\gamma\geq W_0+1.
\tag{4.4}
\]

The unit gap in (4.4) grows linearly under powers.

## 5. A fixed multiplier cannot cross the gap

Separate the unique full-\(t\)-degree part of the translated polynomial:

\[
 P(z+Vt)=B_0(t)+B_{<}(z,t),
\tag{5.1}
\]

where every \(t\)-monomial in \(B_<\) has total degree at most \(d-1\).
This includes both the \(z\)-dependent pieces of \(P_d(z+Vt)\) and all
lower homogeneous parts of \(P\).

Let \(q=\deg Q\).  A term of

\[
 Q(z+Vt)P(z+Vt)^m
\tag{5.2}
\]

that contributes to \(t_1^m\cdots t_d^m\) receives at most \(q\) units of
total \(t\)-degree from the first factor.  Every choice from \(B_<\) loses
at least one of the \(dm\) possible \(t\)-degrees in the second factor.
Consequently at most \(q\) of its \(m\) factors come from \(B_<\), and at
least \(m-q\) come from \(B_0\).

Let

\[
 M=\min\{w\mathbin{\cdot}\gamma:
 t^\gamma\text{ occurs in }P(z+Vt)\},
\qquad
 L=\min(M,W_0+1).
\]

By (4.4), the \(t\)-weight of the selected \(P\)-term is at least

\[
 (m-q)(W_0+1)+qL
 =m(W_0+1)-C,
\qquad
 C=q(W_0+1-L).
\tag{5.3}
\]

There are only finitely many \(t\)-exponents \(\eta\) in \(Q(z+Vt)\).
Choose \(K\) with

\[
 |w\mathbin{\cdot}\eta|\leq K
\]

for all of them.  If a product term in (5.2) had total exponent
\((m,\ldots,m)\), the \(P\)-term would have weight at most

\[
 mW_0+K.
\tag{5.4}
\]

For \(m>C+K\), the lower bound (5.3) is strictly larger than (5.4), a
contradiction.  Thus the coefficient in (2.2) is zero for all sufficiently
large \(m\), proving (1.2).

All data used in this argument lie in a finitely generated
characteristic-zero subfield.  After embedding that field in
\(\mathbb C\), applying the constant-term theorem, and proving the
polynomial identity, the identity descends to the original field.

## 6. Consequences and limits

For \(n=2\), every degree-\(d\) operator symbol is a binary form and hence
splits.  Corollary 1.2 therefore closes every pair with homogeneous
\(\Lambda\) and \(\deg P\leq\operatorname{ord}\Lambda\).  The earlier
exact results remain stronger in a different direction:

- the [cubic theorem](TWO_VARIABLE_CUBIC_GVC_THEOREM.md) replaces the
  infinite premise by the first four moments, identifies the Segre
  nullcone exactly, and gives the sharp simple cutoff \(m>\deg Q\);
- the [low-root calculations](TWO_VARIABLE_LOW_ROOT_GVC_THEOREMS.md) give
  the same explicit cutoff for all two-root symbols and a first-five-moment
  radical certificate on the quartic \((2,1,1)\) orbit.

The split-symbol proof gives an eventual bound depending on the separating
weight; it does not supply a uniform finite pure-moment cutoff.  For the
squarefree quartic orbit, the remaining question is therefore a finite
moment--nullcone certificate, not the GVC conclusion itself.  For GVC(2),
the separable frontier has moved to nonhomogeneous operators or to
polynomials whose degree exceeds the homogeneous operator order.

The Dvorsky counterexample is not contradicted: its five-variable cubic
symbol

\[
 \xi_t(\xi_a\xi_d-\xi_b\xi_c)
\]

contains a nonsplit rank-four quadratic factor.  Likewise, an
ordinary-Laplacian symbol in more than two variables is a nondegenerate
quadratic form of rank greater than two and is not a product of two linear
forms.

No computer algebra is required for Theorem 1.1.  Its external input is
the Duistermaat--van der Kallen theorem for constant terms of powers of a
Laurent polynomial:

- J. J. Duistermaat and W. van der Kallen,
  [*Constant terms in powers of a Laurent polynomial*](https://webspace.science.uu.nl/~kalle101/powers.pdf),
  Indagationes Mathematicae 9 (1998), 221--231.

For the surrounding Vanishing Conjecture context, see also
A. van den Essen, R. Willems, and W. Zhao,
[*Some Results on the Vanishing Conjecture of Differential Operators with Constant Coefficients*](https://arxiv.org/abs/0903.1478).
