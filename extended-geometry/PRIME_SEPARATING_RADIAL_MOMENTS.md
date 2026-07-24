# Weighted constant-term functionals and prime-separating radial moments

## Research assessment (24 July 2026)

The lower-face proof is fundamentally a theorem about the weighted
constant-term functional
\[
 \Phi_a(F)=\sum_{d\geq0}a_d[U^dT^0]F
 \quad\text{on}\quad \mathbb C[U,T,T^{-1}],
\]
not about probability.  Its only inputs are the
Duistermaat--van der Kallen constant-term theorem, Frobenius compatibility,
and a normalized divisibility condition on the weights \(a_d\).

Rotationally invariant probability measures give an important application.
If \(\sigma\) is such a measure on \(\mathbb C\), \(Z\sim\sigma\), and all
polynomial moments exist, write
\[
 \mathbb E[Z^a\bar Z^b]=\delta_{ab}\mu_a,\qquad \mu_0=1.
\]
Equivalently, \(X=|Z|^2\) has the Stieltjes moment sequence
\(\mu_n=\mathbb E[X^n]\).

The lower-face proof of \(\operatorname{GMC}(2)\) does extend to a useful
arithmetic class of such measures.  In particular, it proves the polynomial
moment property for
\[
 d\sigma_s(z)=\frac{1}{\pi s!}|z|^{2s}e^{-|z|^2}\,dA(z),
 \qquad s\in\mathbb Z_{\geq0},
\]
and, more strongly, for every positive **rational** Gamma shape.

The natural publishable statement is first an abstract sufficient criterion
for prime-separating weight sequences, followed by Stieltjes, Gamma,
Gamma-product, factorial-ratio, and \(q\)-factorial corollaries.  It should
not be advertised as a classification: prime separation is not necessary,
and even its restriction to \(T\)-independent polynomials contains the
general problem of classifying moment kernels in \(\mathbb C[U]\).

There is now also a sharp negative theorem.  If \(w(U)\) is any nonzero
polynomial, then the kernel of
\[
 F\longmapsto\int_0^1\operatorname{CT}_T F(U,T)\,w(U)\,dU
\]
is not Mathieu--Zhao.  Long's \((xz)\)-counterexample is the case \(w=1\).
The integer-Beta family gives the moment sequences
\[
 a_d=\frac{(s)_d}{(s+t)_d}
     =\frac{(s)_t}{(d+s)_t},
 \qquad s,t\in\mathbb Z_{>0}.
\]
Every one of these kernels is therefore non-Mathieu.  Thus failure of prime
separation is not merely a defect of the proof: compact polynomial radial
weights and reciprocal Gamma ratios support exact two-sided cancellations.
On the other hand, geometric weights show that prime separation is not
necessary.  The correct status is therefore a strong sufficient theorem
with examples on both sides of its boundary, not an if-and-only-if
criterion.

## 1. Abstract weighted constant-term problem

Fix a nonzero complex sequence \(a=(a_d)_{d\geq0}\) and define
\[
 \Phi_a:\mathbb C[U,T,T^{-1}]\longrightarrow\mathbb C,\qquad
 \Phi_a(F)=\sum_{d\geq0}a_d[U^dT^0]F.                         \tag{1.1}
\]
The sum is finite for every polynomial \(F\).  The classification problem is
to determine those \(a\) for which \(\ker\Phi_a\) is a Mathieu--Zhao
subspace.

There are two immediate necessary cautions.

1. Unless \(\Phi_a=0\), one must have \(a_0\ne0\).  Otherwise
   \(1\in\ker\Phi_a\), whereas a proper Mathieu--Zhao subspace cannot contain
   \(1\).
2. If \(a_{dm}=0\) for every \(m\geq1\), then the Mathieu--Zhao property
   forces \(a_{dm+e}=0\) for all sufficiently large \(m\), for every
   \(e\geq0\).

The theorem below gives a broad sufficient class.  Its conclusion is
stronger than the Mathieu--Zhao property:

> **Abstract lower-face conclusion.** If \(a\) is prime-separating and
> \(\Phi_a(F^m)=0\) for all sufficiently large \(m\) (in particular, for
> every \(m\geq1\)), then the \(T\)-weights of \(F\) all have the same
> strict sign.  Consequently
> \(\Phi_a(GF^m)=0\) for every fixed \(G\) and all \(m\gg0\).

Thus the arithmetic task is to recognize prime-separating sequences.  The
definition and proof are given in Section 3.

## 2. Polynomial moment property

Say that \(\sigma\) has the polynomial moment property (PMP) if, for
\(P,Q\in\mathbb C[z,\bar z]\),
\[
 \mathbb E[P(Z,\bar Z)^m]=0\quad(m\geq1)
 \quad\Longrightarrow\quad
 \mathbb E[Q(Z,\bar Z)P(Z,\bar Z)^m]=0\quad(m\gg0).
\]
This is the statement that the kernel of integration is a
Mathieu--Zhao subspace.  It is the one-complex-variable/two-real-variable
analogue of the Gaussian Moments Conjecture as formulated by
[Derksen--van den Essen--Zhao](https://arxiv.org/abs/1506.05192).

Every polynomial has a rotational-weight decomposition
\[
 P=\sum_{k\in S}T^kB_k(U),\qquad U=z\bar z,
\]
for which
\[
 \mathbb E[P^m]=\mathcal L_\mu(\operatorname{CT}_T P^m),
 \qquad \mathcal L_\mu(U^d)=\mu_d.
\]
Thus only the radial moment functional changes when the Gaussian law is
replaced by another rotationally invariant law.

## 3. The correct arithmetic hypothesis

It is not quite enough to require only
\(a_{jp}/a_{np}\equiv0\pmod p\).  In the Frobenius step, coefficients in
degrees not divisible by \(p\) vanish modulo \(p\).  A negative \(p\)-adic
valuation in their moment ratios could undo that vanishing.  The proof needs
integrality of all intermediate ratios as well.

### Definition (normalized rational prime separation)

A sequence \(a=(a_d)_{d\geq0}\) of nonzero rational numbers has
**normalized prime separation** (or is **prime-separating**) if, for every
pair \(0\leq n<M\), there are infinitely many rational primes \(p\) such
that
\[
 \frac{a_d}{a_{np}}\in\mathbb Z_{(p)}
 \quad(np\leq d\leq Mp),                                      \tag{PS1}
\]
and
\[
 \frac{a_{jp}}{a_{np}}\in p\mathbb Z_{(p)}
 \quad(n<j\leq M).                                             \tag{PS2}
\]
It is enough to have arbitrarily large such primes.  One may replace
\(\mathbb Q\) and \(\mathbb Z_{(p)}\) by a fixed number field and suitable
prime ideals, but the rational version contains the cleanest
applications.

More precisely, if \(K\) is a number field and \(a_d\in K^\times\), require
infinitely many prime ideals \(\mathfrak p\) of residue characteristic \(p\)
such that
\[
 \frac{a_d}{a_{np}}\in\mathcal O_{K,\mathfrak p}
 \quad(np\leq d\leq Mp),                                      \tag{PS1K}
\]
and
\[
 \frac{a_{jp}}{a_{np}}\in\mathfrak p\mathcal O_{K,\mathfrak p}
 \quad(n<j\leq M).                                             \tag{PS2K}
\]
The proof below works verbatim over this integral model.  Degree-one primes
are needed only in applications, such as \(q\)-factorials, where the size of
the residue field controls a multiplicative order.

### Theorem (normalized prime-separation forces one-sided support)

If \(a=(a_d)_{d\geq0}\) is a prime-separating rational sequence, then
\(\ker\Phi_a\) is a Mathieu--Zhao subspace.  More precisely, if
\[
 \Phi_a(P^m)=0\qquad(m\gg0),
\]
then \(0\notin\operatorname{conv}(S)\), where
\(S=\{k:B_k\ne0\}\) is the rotational-weight support of \(P\).

In particular, every rotationally invariant measure whose radial moment
sequence is prime-separating has PMP.

The same conclusion holds for a number-field sequence satisfying
(PS1K)--(PS2K).  This is the promised reusable form of the GMC(2) method:

> **Any one-angular-variable constant-term functional with normalized
> prime separation forces one-sided angular support.**

### Proof

Put
\[
 \nu_k=\operatorname{ord}_U B_k,\qquad
 b_k=[U^{\nu_k}]B_k.
\]
If \(0\in\operatorname{conv}(S)\), minimize the radial order over convex
weight-zero combinations.  A rational supporting line produces a face
\(S_0\), and the theorem of
[Duistermaat--van der Kallen](https://webspace.science.uu.nl/~kalle101/powers.pdf)
gives an \(r\geq1\) such that
\[
 F_r(U):=\operatorname{CT}_T(P^r)
        =cU^n+\sum_{j>n}c_jU^j,\qquad c\ne0.                   \tag{3.1}
\]
The same supporting-line inequality says that every term of
\[
 F_{rp}(U)=\operatorname{CT}_T(P^{rp})
\]
has degree at least \(np\).  Choose an integer \(M>n\), independent of
\(p\), such that \(Mp\) bounds its degree.

Work in the finite-type characteristic-zero domain generated by the
coefficients of \(P\), with \(c^{-1}\) adjoined.  Its morphism to the
relevant arithmetic base has a nonempty fiber away from finitely many
primes; after passing to an algebraic residue field, \(c\) remains nonzero.
Choose a sufficiently large prime supplied by prime separation, so that
the assumed eventual vanishing includes the exponent \(rp\), and divide
\(\Phi_a(P^{rp})=0\) by \(a_{np}\).  Condition (PS1) makes the
resulting identity \(p\)-integral.  In characteristic \(p\),
\[
 F_{rp}=\operatorname{CT}_T((P^r)^p)=F_r^p.
\]
Consequently, all coefficients outside degrees \(jp\) vanish modulo \(p\).
Condition (PS1) ensures that this vanishing survives application of the
normalized moment functional.  For \(j>n\), condition (PS2) kills the
remaining degree-\(jp\) terms.  The normalized degree-\(np\) term is
\(c^p\), so the reduced moment identity gives \(c^p=0\), a contradiction.

Therefore \(0\notin\operatorname{conv}(S)\).  All weights of \(P\) have the
same strict sign, while \(Q\) has bounded weight support, so
\(\operatorname{CT}_T(QP^m)=0\) for \(m\gg0\).  This proves PMP.

## 4. Gamma laws

For \(\alpha>0\), define the normalized circular Gamma law
\[
 d\sigma_\alpha(z)
 =\frac{1}{\pi\Gamma(\alpha)}
   |z|^{2\alpha-2}e^{-|z|^2}\,dA(z).
\]
Then \(X=|Z|^2\sim\operatorname{Gamma}(\alpha,1)\), and
\[
 \mu_n=(\alpha)_n=\frac{\Gamma(n+\alpha)}{\Gamma(\alpha)}.
\]

### Integer shape

For \(\alpha=s+1\), \(s\in\mathbb Z_{\geq0}\),
\[
 \mu_n=\frac{(n+s)!}{s!}.
\]
If \(d\geq np\), then
\[
 \frac{\mu_d}{\mu_{np}}
 =\frac{(d+s)!}{(np+s)!}\in\mathbb Z.
\]
If \(j>n\) and \(p>s\), the product contains \((n+1)p\), hence
\[
 p\mid\frac{\mu_{jp}}{\mu_{np}}.
\]
Thus \(\sigma_{s+1}\) has PMP.  This verifies the proposed extension.

### Rational shape

Let \(\alpha=A/B>0\) in lowest terms.  For \(p\nmid B\),
\[
 \frac{\mu_d}{\mu_{np}}
 =\prod_{k=np}^{d-1}\frac{A+Bk}{B}
\]
is \(p\)-integral.  If \(j>n\), the interval
\([np,jp-1]\) contains a complete residue system modulo \(p\), so one of its
integers \(k\) satisfies \(A+Bk\equiv0\pmod p\).  Therefore (PS2) holds.

> **Gamma corollary.** Every circular Gamma law of positive rational shape
> has PMP.

The scale parameter is irrelevant: a radial dilation conjugates the
polynomial moment problem.

## 5. A broader recurrence criterion

The Gamma calculation is a case of a simple reusable criterion.

### Proposition (polynomial-ratio moments)

Suppose that \(\mu=(\mu_n)_{n\geq0}\) is a sequence of nonzero rational
numbers and
\[
 \frac{\mu_{n+1}}{\mu_n}=H(n)
\]
for a nonconstant \(H\in\mathbb Q[x]\).  Then \(\mu\) is prime-separating
and, if it is a Stieltjes moment sequence, every rotationally invariant
representing measure has PMP.

Indeed,
\[
 \frac{\mu_d}{\mu_{np}}=\prod_{k=np}^{d-1}H(k).
\]
Away from the fixed denominators of \(H\), this is \(p\)-integral.  By
Schur's theorem, a nonconstant integer polynomial has a root modulo
infinitely many primes.  For each such prime, every interval of \(p\)
consecutive integers contains a representative of that root, proving
(PS2).

No positivity is needed for the abstract weighted-functional theorem.
The Stieltjes hypothesis is needed only to interpret the sequence as the
moments of a positive measure; positivity of \(H(n)\) alone is not asserted
to produce such a measure.

### Concrete families

1. **Products of Gamma variables.**  If
   \(X=\prod_{i=1}^qY_i\), with independent
   \(Y_i\sim\operatorname{Gamma}(\alpha_i,1)\) and
   \(\alpha_i\in\mathbb Q_{>0}\), then
   \[
   \mu_n=\prod_i(\alpha_i)_n,\qquad
   \frac{\mu_{n+1}}{\mu_n}=\prod_i(n+\alpha_i).
   \]
   Hence the circular law obtained from \(|Z|^2=X\) and an independent
   uniform angle has PMP.

2. **Integer powers of Gamma variables.**  If
   \(X=Y^a\), \(a\in\mathbb Z_{>0}\), and
   \(Y\sim\operatorname{Gamma}(b,1)\) with
   \(b\in\mathbb Q_{>0}\), then
   \[
   \mu_n=\frac{\Gamma(an+b)}{\Gamma(b)},\qquad
   \frac{\mu_{n+1}}{\mu_n}=(an+b)_a.
   \]
   Products of finitely many such variables work as well.

The moment sequences \(\Gamma(an+b)/\Gamma(b)\), their positive powers, and
their representing densities belong to the established theory of
infinitely divisible Stieltjes moment sequences; see
[Berg, *A two-parameter extension of the Urbanik semigroup*](https://arxiv.org/abs/1802.00993)
and [Berg, *On powers of Stieltjes moment sequences, II*](https://arxiv.org/abs/math/0412340).
The prime-separation argument above covers integer \(a\) and integer powers
of these sequences directly.  It does not automatically cover arbitrary
real \(a\) or arbitrary positive convolution powers.

## 6. Factorial ratios, hypergeometric sequences, and \(q\)-factorials

### Positive-excess integral factorial ratios

Consider
\[
 a_n=\frac{\prod_{i=1}^r(\alpha_i n)!}
           {\prod_{j=1}^s(\beta_j n)!},
 \qquad
 \Delta=\sum_i\alpha_i-\sum_j\beta_j,
                                                               \tag{6.1}
\]
with positive integral \(\alpha_i,\beta_j\), and put
\[
 L(x)=\sum_i\lfloor\alpha_i x\rfloor
      -\sum_j\lfloor\beta_jx\rfloor .                          \tag{6.2}
\]
Landau's criterion says that \(a_n\) is integral for every \(n\) exactly
when \(L(x)\geq0\) for every \(x\in[0,1]\).

Assume integrality and \(\Delta>0\).  For fixed \(n<M\) and all sufficiently
large \(p\), Legendre's formula gives, uniformly for \(np\leq d\leq Mp\),
\[
 v_p(a_d)=L(d/p),\qquad v_p(a_{np})=n\Delta.                   \tag{6.3}
\]
Since \(L(n+x)=n\Delta+L(x)\), the first quantity is at least the
second.  At the endpoints,
\[
 v_p(a_{jp})-v_p(a_{np})=(j-n)\Delta>0.                       \tag{6.4}
\]
Hence every integral factorial ratio of positive factorial excess is
prime-separating.

Balanced ratios with \(\Delta=0\), such as
\(\binom{2n}{n}\), fail this endpoint-isolation test.  This does not show
that their kernels are not Mathieu--Zhao; it shows only that the present
prime-isolation proof cannot decide them.  General integer-valued
hypergeometric sequences must likewise be separated according to their
\(p\)-adic excess rather than integrality alone.

Factorial excess therefore gives the following useful, but not complete,
frontier:

* positive integral excess plus Landau integrality implies prime separation;
* slope-balanced integral ratios such as \(\binom{2n}{n}\) are not decided
  by endpoint isolation;
* balanced Gamma ratios with a positive shift deficit can already fail.
  The integer-Beta sequence
  \[
   \frac{(s)_n}{(s+t)_n}=\frac{(s)_t}{(n+s)_t}
  \]
  gives an infinite counterexample family by the theorem below.

### \(q\)-factorials

Let
\[
 a_n=(q;q)_n=\prod_{k=1}^n(1-q^k)                             \tag{6.5}
\]
with nonzero algebraic \(q\) that is not a root of unity.  Use the
number-field version (PS1K)--(PS2K) and choose rational primes that split
completely in \(\mathbb Q(q)\).  Chebotarev supplies infinitely many such
primes.  Away from finitely many of them,
\(\bar q\in\mathbb F_p^\times\), and \(q\) and \(q^{-1}\) are integral
units.  The normalized ratio
\[
 \frac{a_d}{a_{np}}=\prod_{k=np+1}^{d}(1-q^k)                 \tag{6.6}
\]
is integral.  Every interval of \(p\) consecutive exponents contains a
multiple of \(\operatorname{ord}(\bar q)\leq p-1\), so the endpoint ratios
with \(d=jp\), \(j>n\), vanish modulo \(p\).  Thus these algebraic
\(q\)-factorials satisfy the number-field prime-separation criterion.

Pure cyclotomic specialization in characteristic zero is not a substitute
for this reduction: the proof still needs
\((X+Y)^p=X^p+Y^p\).  A genuinely characteristic-zero \(q\)-proof would
require an additional quantum-Frobenius structure.

### Candidate ledger

| Sequence | Result of the lower-face method |
|---|---|
| \(n!\), \((n+s)!\) | Prime-separating |
| \((\alpha)_n\), positive rational \(\alpha\) | Prime-separating |
| Products of rational rising factorials | Prime-separating |
| \(\prod_{k<n}H(k)\), nonconstant \(H\in\mathbb Q[k]\) | Prime-separating when all terms are nonzero |
| Integral factorial ratios with \(\Delta>0\) | Prime-separating |
| Balanced integral factorial ratios such as \(\binom{2n}{n}\) | Not decided by prime isolation |
| Integer-Beta ratios \((s)_n/(s+t)_n\) | Not Mathieu--Zhao; explicit two-sided circuit |
| Moments of any nonzero polynomial density on \([0,1]\) | Not Mathieu--Zhao |
| Algebraic non-torsion \(q\)-factorials | Prime-separating over a number-field integral model |
| Geometric sequences \(c\lambda^n\) | Mathieu--Zhao by evaluation, not by prime separation |

## 7. Why this is not yet a classification theorem

### Prime separation is sufficient, not necessary

Let \(\sigma\) be uniform measure on the circle \(|z|=R\).  Then
\(\mu_n=R^{2n}\), whose ratios are \(p\)-adic units rather than
prime-separating.  Nevertheless PMP follows directly from the
Duistermaat--van der Kallen theorem applied to the Laurent polynomial
\(P(RT,R/T)\).  Thus no theorem of the form
\[
 \text{PMP}\quad\Longleftrightarrow\quad\text{prime separation}
\]
can hold.

Algebraically, if \(a_n=c\lambda^n\), then
\[
 \Phi_a(F)=c\,\operatorname{CT}_T F(\lambda,T).                \tag{7.1}
\]
This is the pullback of the ordinary constant-term functional under the
algebra homomorphism \(U\mapsto\lambda\), so its kernel is
Mathieu--Zhao by Duistermaat--van der Kallen.  This proves directly that
prime separation is sufficient rather than necessary.

### Measures versus moment sequences

PMP for polynomials depends only on \((\mu_n)\).  A Stieltjes moment sequence
can have several representing measures.  Berg's results, for example, show
moment indeterminacy in part of the generalized-Gamma product-convolution
family.  Consequently, the arithmetic object to classify is the radial
moment functional, not the measure itself; all representing circular
measures then share the same polynomial integrals.

### No polynomial interval density has a Mathieu kernel

Let \(0\ne w\in\mathbb C[U]\), and define
\[
 \Phi_w(F)=\int_0^1\operatorname{CT}_T F(U,T)\,w(U)\,dU.       \tag{7.2}
\]
Write
\[
 w(U)=(1-U)^{t-1}h(U),\qquad h(1)\ne0,                       \tag{7.3}
\]
where \(t\geq1\), and put \(D=\deg h\).  For the same three-weight circuit
\[
 f(U,T)=(1-T^{-1})((1-U)+UT),                                \tag{7.4}
\]
one has, for every \(m\geq D+t\),
\[
 \Phi_w(f^m)=0,\qquad
 \Phi_w(T^{-t}f^m)
 =(-1)^{m+t}\frac{h(1)}{t\binom{m+t}{t}}\ne0.                \tag{7.5}
\]

Indeed, expand \(h(U)=\sum_{r=0}^Dc_rU^r\).  The coefficient of \(T^k\)
in \(((1-U)+UT)^m\), integrated against
\(U^r(1-U)^{t-1}\), is
\[
 \binom mk B(k+r+1,m-k+t).                                   \tag{7.6}
\]
As a function of \(k\), after removing a nonzero factor depending only on
\(m,r,t\), this is a polynomial of degree \(r+t-1\).  The pure constant
term is its \(m\)-th alternating finite difference and vanishes for
\(m\geq r+t\).

The mixed calculation is the unnormalized \(s=r+1\) case of the
integer-Beta identity below.  Its value is
\[
 (-1)^{m+t}\frac{1}{t\binom{m+t}{t}},                         \tag{7.7}
\]
independent of \(r\).  Summing over \(r\) replaces \(1\) by
\(\sum_rc_r=h(1)\), proving (7.5).

Set \(N=D+t\), \(P=f^N\), and \(Q=T^{-t}\).  Then
\[
 \Phi_w(P^\ell)=0,\qquad \Phi_w(QP^\ell)\ne0
 \quad(\ell\geq1).                                            \tag{7.8}
\]

> **Polynomial-density counterexample theorem.** The kernel of
> \(\Phi_w\) is not a Mathieu--Zhao subspace for any nonzero polynomial
> interval weight \(w\).  The witness always comes from a power of the
> same three-weight circuit (7.4).

The result includes every positive polynomial probability density on the
interval and, after \(U=|z|^2\), every circular disk law with nonzero
polynomial radial density.  It gives a broad negative answer for
one-angular-weight Mathieu problems with compact polynomial radial weights.

### Integer-Beta Gamma ratios as an exact subfamily

Let \(s,t\in\mathbb Z_{>0}\), and define the normalized Beta functional
\[
 \mathcal B_{s,t}(h)
 =\frac1{B(s,t)}\int_0^1h(U)U^{s-1}(1-U)^{t-1}\,dU.           \tag{7.9}
\]
Its moments are the balanced Gamma ratios
\[
 a_d=\mathcal B_{s,t}(U^d)
 =\frac{(s)_d}{(s+t)_d}
 =\frac{(s)_t}{(d+s)_t}.                                     \tag{7.10}
\]
For every sufficiently large prime \(p\),
\[
 \frac{a_{jp}}{a_{np}}
 =\frac{(np+s)_t}{(jp+s)_t}\equiv1\pmod p,                   \tag{7.11}
\]
so endpoint prime separation fails maximally.

Put
\[
 N=s+t-1.                                                     \tag{7.12}
\]
Then for every \(m\geq N\),
\[
 \Phi_{s,t}(f^m)=0,\qquad
 \Phi_{s,t}(T^{-t}f^m)
 =(-1)^{m+t}
   \frac{\binom{s+t-1}{t}}{\binom{m+t}{t}}\ne0,               \tag{7.13}
\]
where \(\Phi_{s,t}=\mathcal B_{s,t}\circ\operatorname{CT}_T\).
Consequently, with \(P=f^N\) and \(Q=T^{-t}\),
\[
 \Phi_{s,t}(P^\ell)=0,\qquad
 \Phi_{s,t}(QP^\ell)\ne0
 \quad(\ell\geq1).                                            \tag{7.14}
\]

> **Integer-Beta counterexample theorem.** For every
> \(s,t\in\mathbb Z_{>0}\), the kernel of \(\Phi_{s,t}\) is not a
> Mathieu--Zhao subspace.  Equivalently, every circular disk law for which
> \(U=|Z|^2\) has a \(\operatorname{Beta}(s,t)\) distribution fails PMP.

To prove (7.13), the coefficient of \(T^k\) in
\(((1-U)+UT)^m\), after applying \(\mathcal B_{s,t}\), is
\[
 \frac{m!}{(s+t)_m}
 \binom{k+s-1}{s-1}
 \binom{m-k+t-1}{t-1}.                                      \tag{7.15}
\]
Apart from the first factor, this is a polynomial in \(k\) of degree
\(s+t-2=N-1\).  Taking the constant term after multiplication by
\((1-T^{-1})^m\) is its \(m\)-th alternating finite difference, hence is
zero for \(m\geq N\).

For the mixed moment, the angular equation is \(k=j+t\).  Apply the same
finite-difference identity to
\[
 \binom{j+t+s-1}{s-1}\binom{m-j-1}{t-1}.                    \tag{7.16}
\]
The omitted values \(j=m-t+1,\ldots,m-1\) vanish, while the polynomial
value at \(j=m\) uses
\(\binom{-1}{t-1}=(-1)^{t-1}\).  The remaining boundary term simplifies
to the nonzero expression in (7.13).

For \(s=t=1\), \(N=1\), this is exactly Long's counterexample:
\[
 a_d=\frac1{d+1},\qquad
 \Phi_{1,1}(T^{-1}f^m)=\frac{(-1)^{m-1}}{m+1}.
\]
Under \(U=|z|^2\), this first member is normalized area measure on the
complex unit disk.  The general theorem shows that the failure persists
for every integer Jacobi/Beta radial density, not only uniform area.

Consequently, every mixed one-interval/one-or-more-angular-variable
\((xz)\) problem with uniform interval weight also fails by padding Long's
first example.  More generally, the theorem supplies failures for all
integer-Beta interval weights.  This sharply
limits any attempted extension based only on the existence of one angular
weight: arithmetic separation of the radial functional is essential.

### A general Bernstein finite-difference obstruction

The preceding counterexamples admit an abstract formulation dual to prime
separation.  For a radial functional \(\mathcal L\), put
\[
 b_{m,k}
 =\binom mk\mathcal L\!\left(U^k(1-U)^{m-k}\right)
 \qquad(0\leq k\leq m).                                      \tag{7.17}
\]
For the circuit \(f\) in (7.4), direct angular extraction gives
\[
 \Phi(f^m)=\sum_{k=0}^m(-1)^k\binom mk b_{m,k},                \tag{7.18}
\]
and, for a fixed \(q\geq1\),
\[
 \Phi(T^{-q}f^m)
 =\sum_{j=0}^{m-q}(-1)^j\binom mj b_{m,j+q}.                  \tag{7.19}
\]

> **Bernstein finite-difference criterion.** Suppose that, for all large
> \(m\), the values \(b_{m,k}\) are the restrictions of polynomials
> \(R_m(k)\) of uniformly bounded degree.  Then \(f\) lies in the radical
> of \(\ker\Phi\).  If, for some fixed \(q\), the truncated difference
> (7.19) is nonzero for infinitely many \(m\), then \(\ker\Phi\) is not
> Mathieu--Zhao.

The first assertion is simply the vanishing of a sufficiently high finite
difference.  The second is the definition of the Mathieu obstruction.  For
polynomial interval densities, (7.6) makes the Bernstein degree uniformly
bounded and (7.7) computes a surviving boundary term.

This criterion identifies a second structural axis for the classification
problem:

* prime separation prevents cancellation by isolating one radial layer;
* bounded Bernstein complexity creates cancellation by turning the pure
  moment into a high finite difference.

Nonintegral Beta laws, arcsine/central-binomial moments, and more general
balanced hypergeometric sequences usually have unbounded Bernstein
complexity and no prime gap.  They lie between the two theorems and form the
next genuine frontier.

### The univariate obstruction to a complete classification

Restricting \(\Phi_a\) to \(T\)-independent polynomials gives the moment
functional
\[
 \mathcal L_a\!\left(\sum f_nU^n\right)=\sum f_na_n.           \tag{7.20}
\]
Therefore any full classification of the weighted constant-term kernels
contains, as a necessary subproblem, the classification of codimension-one
moment kernels in \(\mathbb C[U]\).  The radical characterization of
univariate Mathieu subspaces by van den Essen and Zhao gives an abstract
criterion, but not a closed arithmetic classification of arbitrary
sequences.

## 8. Recommended paper framing

A concise companion section or note could be organized as follows:

1. define the abstract functional \(\Phi_a\) and prime separation using
   (PS1)--(PS2);
2. prove the weighted prime-separating lower-face theorem;
3. state the positive-rational-shape circular Gamma corollary;
4. give the polynomial-ratio, positive-excess factorial-ratio,
   Gamma-product, powered-Gamma, and algebraic \(q\)-factorial families;
5. contrast the geometric positive example (no prime separation needed)
   with the polynomial-density and integer-Beta counterexample theorems;
6. pose the genuine classification problem:
   characterize sequences \(a\) for which \(\ker\Phi_a\) is a
   Mathieu--Zhao subspace.

The integer-shape Gamma theorem is essentially immediate once the GMC(2)
lower-face proof is available.  The rational-shape and polynomial-ratio
extensions make the result substantially cleaner and broader.  As a
standalone publication, the strongest version should include either this
broader family or a necessity/alternative-mechanisms analysis; the
integer-\(s\) corollary alone is better presented alongside the GMC(2)
paper.

## References

- J. J. Duistermaat and W. van der Kallen,
  [*Constant terms in powers of a Laurent polynomial*](https://webspace.science.uu.nl/~kalle101/powers.pdf),
  Indagationes Mathematicae 9 (1998), 221--231.
- H. Derksen, A. van den Essen, and W. Zhao,
  [*The Gaussian Moments Conjecture and the Jacobian Conjecture*](https://arxiv.org/abs/1506.05192).
- C. Berg,
  [*A two-parameter extension of the Urbanik semigroup*](https://arxiv.org/abs/1802.00993).
- C. Berg,
  [*On powers of Stieltjes moment sequences, II*](https://arxiv.org/abs/math/0412340).
- A. van den Essen and W. Zhao,
  [*Mathieu Subspaces of Univariate Polynomial Algebras*](https://arxiv.org/abs/1012.2017).
- J. W. Bober,
  [*Factorial ratios, hypergeometric series, and a family of step functions*](https://arxiv.org/abs/0709.1977).
- C. D. Long,
  [*Counterexamples to the \(xz\)-Conjecture and the Mathieu Conjecture for
  \(SU(2)\)*](https://arxiv.org/abs/2607.19012).
