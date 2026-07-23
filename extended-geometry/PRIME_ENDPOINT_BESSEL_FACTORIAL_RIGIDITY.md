# Prime-endpoint rigidity for the Bessel--factorial transform

## 1. Statement

Let \(R\) be a characteristic-zero domain and define the factorial
functional
\[
\mathcal L_R\!\left(\sum_{j\geq0}f_jU^j\right)
   =\sum_{j\geq0}j!f_j.
\]
For \(C,D\in R[U]\), put
\[
H_m(C,D)=
\sum_{r=0}^{\lfloor m/2\rfloor}
\frac{m!}{(r!)^2(m-2r)!}C^{m-2r}D^r.             \tag{1.1}
\]

The purpose of this note is to prove the following result.

> **Bessel--factorial rigidity theorem.** Let \(C,D\in\mathbb C[U]\).  If
> \[
> \mathcal L_{\mathbb C}(H_m(C,D))=0
> \qquad\text{for every }m\geq1,
> \]
> then \(C=D=0\).

The proof uses only reduction modulo rational primes and the \(U\)-adic
orders of \(C\) and \(D\).

## 2. The prime-endpoint lemma

> **Prime-endpoint lemma.** Let \(p\) be an odd prime and let \(C,D\) belong
> to a polynomial ring over an \(\mathbb F_p\)-algebra.  Then
> \[
> H_p(C,D)=C^p,\qquad
> H_{2p}(C,D)=C^{2p}+2D^p.                        \tag{2.1}
> \]

**Proof.** For \(1\leq r\leq(p-1)/2\), the numerator of
\[
\frac{p!}{(r!)^2(p-2r)!}
\]
contains one factor \(p\), while its denominator contains none.  Thus every
nonzero-\(r\) coefficient of \(H_p\) vanishes modulo \(p\).

For \(1\leq r<p\), Legendre's formula gives
\[
v_p((2p)!)=2,\quad v_p(r!)=0,\quad
v_p((2p-2r)!)\leq1.
\]
Hence every nonendpoint coefficient of \(H_{2p}\) vanishes modulo \(p\).
At the second endpoint,
\[
\frac{(2p)!}{(p!)^2}=\binom{2p}{p}\equiv2\pmod p.
\]
This proves (2.1). \(\square\)

The lemma is an identity before the factorial functional is applied.  The
remaining argument shows that one of its endpoints also survives
\(\mathcal L\).

## 3. Direct reduction of finitely generated coefficients

We use the following standard spreading fact.

> **Finite-type reduction lemma.** Let \(R\) be a finitely generated
> integral domain over \(\mathbb Z\) of characteristic zero.  Then, for all
> but finitely many rational primes \(p\), there is a prime ideal
> \(\mathfrak p\subset R\) with \(\mathfrak p\cap\mathbb Z=(p)\).

**Proof.** The finite-type morphism
\(\operatorname{Spec}R\to\operatorname{Spec}\mathbb Z\) contains the generic
point in its image because \(R\otimes_{\mathbb Z}\mathbb Q\ne0\).  By
Chevalley's theorem, its constructible image therefore contains a nonempty
open subset of \(\operatorname{Spec}\mathbb Z\).  Every prime in that open
subset has a nonempty fiber. \(\square\)

This replaces algebraic-coefficient reduction.  Given finitely many complex
coefficients and selected nonzero coefficients \(c,d\), use the literal
subring
\[
R=\mathbb Z[\text{coefficients},c^{-1},d^{-1}]
   \subset\mathbb C.                              \tag{3.1}
\]
It is a finite-type domain.  For almost every \(p\), reduction at a prime
\(\mathfrak p\) above \(p\) is available, and the images of \(c,d\) remain
nonzero because they are units in \(R\).  No Nullstellensatz specialization,
number field, Dedekind extension, or discrete valuation is required.

If literal valuation language is desired, choose a local ring
\(R_{\mathfrak p}\) at one of these characteristic-\(p\) places and a
valuation ring of \(\operatorname{Frac}(R)\) dominating it.  The selected
coefficients have value zero, while \(p\) and every element of
\(\mathfrak pR_{\mathfrak p}\) have positive value.  The proof below then
says that a sum with a unique value-zero endpoint cannot vanish.  Reduction
modulo \(\mathfrak p\) is the shorter version of exactly this argument.

## 4. Endpoint isolation after the factorial functional

Assume first that \(C,D\ne0\), and write
\[
a=\operatorname{ord}_U C,\qquad
e=\operatorname{ord}_U D.                         \tag{4.1}
\]
Adjoin the inverses of the coefficients of \(U^a\) and \(U^e\) to the ring
\(R\) in (3.1), and choose a sufficiently large odd prime and reduction
\(R\to\kappa(\mathfrak p)\) as in Section 3.

If \(e\geq2a\), use the equation of index \(p\).  Every monomial in
\(C^{p-2r}D^r\) has degree at least
\[
a(p-2r)+er\geq ap.                                \tag{4.2}
\]
Consequently \((ap)!\) divides every factorial produced by
\(\mathcal L_R(H_p)\), so it may be cancelled in the characteristic-zero
domain \(R\).  After cancellation and reduction modulo \(\mathfrak p\), the
terms with Bessel index \(r\geq1\) vanish by the coefficient divisibility in
the prime-endpoint lemma.  Within \(C^p\), non-Frobenius multinomial terms
vanish modulo \(p\), and the only remaining term is
\[
c_a^p.                                            \tag{4.3}
\]
Indeed, the higher Frobenius monomials of \(C^p\) contain a factorial ratio
\((jp)!/(ap)!\) with \(j>a\), which crosses the multiple \((a+1)p\).
Thus (4.3) is zero, contradicting that \(c_a\) is a unit.

If \(e<2a\), use the equation of index \(2p\).  Every mixed monomial has
degree at least
\[
2a(p-r)+er=ep+(2a-e)(p-r)\geq ep.                 \tag{4.4}
\]
Cancel \((ep)!\).  The terms with Bessel index \(1\leq r<p\) vanish by the
coefficient divisibility in the prime-endpoint lemma.  Every monomial of
\(C^{2p}\) has degree at least
\(2ap>ep\), and its factorial quotient crosses \((e+1)p\).  Frobenius on
the remaining endpoint \(2D^p\) kills its non-Frobenius multinomial terms
and leaves only
\[
2d_e^p,                                           \tag{4.5}
\]
again contradicting that the odd-prime reduction of \(d_e\) is nonzero.

If \(D=0\ne C\), the first argument applies directly to \(H_p=C^p\).  If
\(C=0\ne D\), the second applies directly to
\(H_{2p}=\binom{2p}{p}D^p\).  Hence the only pair with all transformed
moments zero is \(C=D=0\), proving the theorem.

## 5. Three-level Gaussian consequence

Let
\[
P=WA(U)+C(U)+ZB(U),\qquad U=WZ.
\]
Angular constant-term extraction gives
\[
\mathbb E(P^m)=
\mathcal L_{\mathbb C}\!\left(
H_m(C,UAB)
\right).                                         \tag{5.1}
\]
The theorem forces \(C=0\) and \(UAB=0\).  Since \(\mathbb C[U]\) is a
domain,
\[
C=0,\qquad A=0\ \text{or}\ B=0.                  \tag{5.2}
\]
Thus every pure-moment-zero three-level polynomial has one-sided rotational
support.  For every fixed multiplier \(Q\), rotational weight then gives
\(\mathbb E(QP^m)=0\) for all sufficiently large \(m\).  Therefore the
entire three-level family satisfies \(\operatorname{GMC}(2)\) in every
degree.

The prime-endpoint lemma is the reusable part of the proof.  For larger
support graphs, the corresponding problem is to determine whether one
primitive zero-weight endpoint remains isolated after both prime reduction
and the factorial functional.

## 6. Exact audits

Two implementations independently check the arithmetic skeleton:

- [`audit_prime_endpoint_rigidity_independent.py`](../scripts/audit_prime_endpoint_rigidity_independent.py)
  uses a pure-Python sparse-polynomial implementation to verify (2.1) and
  both normalized endpoint cases;
- [`verify_two_real_gmc_three_level_rigidity.py`](../scripts/verify_two_real_gmc_three_level_rigidity.py)
  repeats the calculation through independent SymPy expansions.

These are finite regressions for the coefficient and factorial congruences;
the all-degree coverage comes from the proof above.
