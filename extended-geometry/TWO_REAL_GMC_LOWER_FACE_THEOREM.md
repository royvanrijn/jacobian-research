# The lower-face prime theorem for GMC(2)

## 1. Statement

Write every polynomial in two circular Gaussian variables as
\[
P=\sum_{k\in S}T^kB_k(U),
\qquad U=ZW,\qquad
\mathbb E(P^m)=\mathcal L(\operatorname{CT}_T P^m),               \tag{1.1}
\]
where \(S\subset\mathbb Z\) is finite and
\(\mathcal L(U^j)=j!\).

> **Lower-face theorem.** If
> \[
> \mathcal L(\operatorname{CT}_T P^m)=0
> \qquad(m\geq1),                                 \tag{1.2}
> \]
> then \(0\notin\operatorname{conv}(S)\).  Equivalently, every nonzero
> rotational weight of \(P\) has the same strict sign.

> **Corollary.** The Gaussian Moments Conjecture holds in two real
> variables.

This contains the three-level, unit-star, arbitrary-star, and circuit
results without a graph hypothesis.

## 2. The exposed lower face

Discard zero coefficient polynomials and put
\[
\nu_k=\operatorname{ord}_U B_k,\qquad
b_k=[U^{\nu_k}]B_k\ne0.                          \tag{2.1}
\]
Assume for contradiction that \(0\in\operatorname{conv}(S)\).  Minimize
\[
\sum_k\lambda_k\nu_k                             \tag{2.2}
\]
subject to
\[
\lambda_k\geq0,\qquad
\sum_k\lambda_k=1,\qquad
\sum_k k\lambda_k=0.
\]
Let the minimum be \(\rho\).  Linear-programming duality supplies
\(\theta\in\mathbb Q\) such that
\[
\nu_k\geq\rho+\theta k\qquad(k\in S).             \tag{2.3}
\]
Let
\[
S_0=\{k:\nu_k=\rho+\theta k\},\qquad
P_0(T)=\sum_{k\in S_0}b_kT^k.                    \tag{2.4}
\]
An optimal solution of (2.2) is supported on \(S_0\), so
\[
0\in\operatorname{conv}(S_0).                    \tag{2.5}
\]

By the
[Duistermaat--van der Kallen constant-term theorem](https://doi.org/10.1016/S0019-3577(98)80020-7),
(2.5) implies that for some \(r\geq1\),
\[
c:=\operatorname{CT}_T(P_0^r)\ne0.               \tag{2.6}
\]
Put
\[
F_m(U)=\operatorname{CT}_T(P^m).                 \tag{2.7}
\]
For every weight-zero monomial of length \(m\), (2.3) gives radial degree
at least
\[
\sum\nu_k\geq \rho m+\theta\sum k=\rho m.         \tag{2.8}
\]
At \(m=r\), equality occurs exactly on the exposed face.  Therefore
\[
n:=\rho r\in\mathbb Z_{\geq0},\qquad
F_r(U)=cU^n+\text{terms of degree greater than }n.               \tag{2.9}
\]

## 3. Prime isolation after the factorial functional

Use the literal finite-type domain generated over \(\mathbb Z\) by all
coefficients of \(P\), with \(c^{-1}\) adjoined.  For all but finitely many
rational primes \(p\), this ring has a reduction of characteristic \(p\)
in which \(c\) remains nonzero.

Fix such a prime.  Equation (2.8) shows that every term of \(F_{rp}\) has
degree at least \(np\).  Hence \((np)!\) can be cancelled from
\[
\mathcal L(F_{rp})=0.                             \tag{3.1}
\]
In characteristic \(p\), Frobenius and constant-term extraction give
\[
F_{rp}
=\operatorname{CT}_T(P^{rp})
=\operatorname{CT}_T((P^r)^p)
=F_r^p.                                          \tag{3.2}
\]
Using (2.9),
\[
F_r^p\equiv
c^pU^{np}+\sum_{j>n}c_j^pU^{jp}\pmod p.          \tag{3.3}
\]
After applying \(\mathcal L/(np)!\), every higher term in (3.3) vanishes
modulo \(p\), because
\[
\frac{(jp)!}{(np)!}
\]
crosses the multiple \((n+1)p\).  Equation (3.1) therefore reduces to
\[
c^p=0,
\]
contradicting the choice of reduction.  This proves the lower-face theorem.

## 4. The GMC consequence

The conclusion \(0\notin\operatorname{conv}(S)\) says that either every
weight of \(P\) is positive or every weight is negative.  Let \(Q\) be any
fixed polynomial.  The weights of \(Q\) are bounded, whereas every weight
of \(P^m\) tends uniformly to the same strict side as \(m\to\infty\).
Consequently
\[
\operatorname{CT}_T(QP^m)=0
\]
for all sufficiently large \(m\), and therefore
\[
\mathbb E(QP^m)=0.
\]
This is precisely \(\operatorname{GMC}(2)\).

## 5. Interpretation

The decisive object is not a circuit graph or even the Hilbert basis of its
relation semigroup.  It is the exposed lower face of the radial-order
Newton polygon
\[
\operatorname{conv}\{(k,\nu_k):k\in S\}.
\]
The constant-term theorem finds a nonzero angular moment on that face;
prime dilation transports it to a Gaussian moment whose lowest factorial
term cannot cancel.  Leaf removal, star rigidity, and toric tie analysis
are finite presentations of this single lower-face obstruction.
