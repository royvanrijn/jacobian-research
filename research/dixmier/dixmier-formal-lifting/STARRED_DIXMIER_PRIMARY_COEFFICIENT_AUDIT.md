# A coefficient correction in the starred Dixmier preprint

This note audits one displayed step in Moskowicz--Valqui,
*The starred Dixmier conjecture for \(A_1\)*,
arXiv:1401.5141v1, 21 January 2014, Proposition 2.6, printed pages 4--5.
The exact inspected version is
[retained here](sources/moskowicz-valqui-1401.5141v1.pdf).

**Finding.** The claimed odd part of the next coefficient in its
commutator calculation is incorrect. The correct odd part is a
derivative of the preceding coefficient and vanishes automatically
when that preceding coefficient vanishes. Thus the displayed argument
does not establish the deduction \(h''=0\).

This is an obstruction to importing that argument as checked proof.
It is not a counterexample to the proposition or to the starred
Dixmier theorem. The later publisher edition was not inspected in
this audit; no assertion is made about its proof.

## 1. Conventions and the precise identity

Work over a characteristic-zero field \(K\) in
\[
 A_1(K)=K\langle X,Y\mid YX-XY=1\rangle.
\]
Use the normal order \(X^iY^j\), exactly as in the inspected source.
Let \(\beta\) be the coefficient-fixing antiinvolution
\[
 \beta(X)=X,\qquad \beta(Y)=-Y.
\]
For an operator \(A\), its \(X\)-degree is its degree when expanded as
\(\sum_i X^i a_i(Y)\).

Let \(\beta(P_0)=P_0\), \(\beta(P_1)=-P_1\), and suppose their
positive \(X\)-degrees are \(r,s\). Write
\[
 \begin{split}
 P_0&=X^r g(Y)+X^{r-1}g_1(Y)+\cdots,\\
 P_1&=X^s h(Y)+X^{s-1}h_1(Y)+\cdots .
 \end{split}                                                   \tag{C1}
\]
Primes below differentiate polynomials in \(Y\). Adjoint parity gives
\[
 \begin{split}
 &g\text{ even},\qquad h\text{ odd},\\
 &g_1=g_{10}+\frac r2g',\qquad
 h_1=\frac s2h'+h_{11},
 \end{split}                                                   \tag{C2}
\]
where \(g_{10}\) is even and \(h_{11}\) is odd. These are the source's
equations (2.2)--(2.3).

**Proposition CorrectNextOddCoefficient.** Put
\[
 T=s g'h-r h'g.
\]
The coefficient of \(X^{r+s-1}\) in \([P_0,P_1]\) is \(T\). If \(C\)
is the coefficient of \(X^{r+s-2}\), then
\[
 \boxed{\quad
 C_{\rm odd}=\frac{r+s-1}{2}\,T'
 =\frac{r+s-1}{2}
       \bigl(sg''h+(s-r)g'h'-rh''g\bigr).
 \quad}                                                       \tag{C3}
\]
In particular \(T=0\) already implies \(C_{\rm odd}=0\); this condition
does not yield a remaining term \(2rgh''\).

**Proof.** The normal-order identity
\[
 a(Y)X^j=\sum_{k=0}^j\binom jk X^{j-k}a^{(k)}(Y)
\]
gives
\[
 \begin{split}
 C={}&\binom s2g''h-\binom r2h''g
 +(s-1)g'h_1-rh_1'g\\
 &\hspace{24mm}+s g_1'h-(r-1)h'g_1 .
 \end{split}                                                   \tag{C4}
\]
Terms omitted from (C1) contribute only to lower \(X\)-degrees.
Substituting (C2), the coefficients of the three odd products
\(g''h,h''g,g'h'\) are respectively
\[
 \frac{s(r+s-1)}2,\qquad
 -\frac{r(r+s-1)}2,\qquad
 \frac{(s-r)(r+s-1)}2 .
\]
The remaining expression
\[
 (s-1)g'h_{11}-r h_{11}'g
       +s g_{10}'h-(r-1)h'g_{10}
\]
is even. This proves (C3). \(\square\)

The source assumes \([P_0,P_1]=1/2\). In its branch \(r,s>0\),
the positive-degree coefficient \(T\) must vanish. Equation (C3)
shows that the next odd-part equation is then automatic, including
the case \(r+s-2=0\), where the constant commutator is even.
The source instead retains the term \(2rgh''\) and uses it to conclude
\(h''=0\). The displayed substitution does not supply that conclusion.

## 2. A small exact witness to the failed identity

Let
\[
 H=XY+\frac12,\qquad P_0=H^2,\qquad P_1=H^3.
\]
Then \(\beta(H)=-H\), so \(P_0\) is symmetric and \(P_1\) is
skew-symmetric. Direct normal ordering gives
\[
 \begin{split}
 P_0&=X^2Y^2+2XY+\frac14,\\
 P_1&=X^3Y^3+\frac92X^2Y^2+\frac{13}4XY+\frac18 .
 \end{split}                                                   \tag{C5}
\]
Here
\[
 r=2,\quad s=3,\quad g=Y^2,\quad h=Y^3,\quad T=0.
\]
The operators commute, so their actual next odd coefficient is zero.
The alleged remaining term is
\[
                         2rgh''=24Y^3\ne0.                    \tag{C6}
\]

This witness has commutator zero, not \(1/2\). It therefore does not
contradict the proposition's conclusion under its full hypotheses.
Its role is narrower: the disputed substitution is an algebraic
coefficient identity, and (C5) satisfies its adjoint parities, leading
coefficient equation and all positive-\(X\)-degree commutator
vanishings. The identity already fails on that input.

For a genuine nonzero-commutator control in the source's conventions,
take
\[
 P_0=-\frac X6+Y^6+2Y^2+3,\qquad P_1=3Y .
                                                               \tag{C7}
\]
Then \(\beta(P_0)=P_0\), \(\beta(P_1)=-P_1\), and
\([P_0,P_1]=1/2\). This confirms the sign convention and supplies an
actual instance of the displayed target normal form.

## 3. What this audit establishes

The correction (C3) is a symbolic identity for arbitrary positive
integers \(r,s\) in characteristic zero. It identifies a missing
deduction in the inspected proof. Another argument could still prove
the same proposition or theorem; neither a theorem counterexample nor
an impossibility of repair is claimed.

Consequently the inspected preprint does not, by this displayed
calculation alone, provide a checked route from a selfadjoint
polynomial Weyl image to the proposed unrestricted normal form.
No stronger result is imported from it here.

The existing
[selfadjoint power-of-two classification](SELFADJOINT_POWER_OF_TWO_STRING_CLASSIFICATION.md)
uses the retained
[coprime algebraic string obstruction](ALL_COPRIME_ALGEBRAIC_STRING_OBSTRUCTION.md).
The starred preprint is not one of its premises. A fresh internal
read-through of that separate analytic and algebraic route is recorded
in the source ledger with its exact inspected steps. This audit does
not upgrade its assurance or extend its scope.

## 4. Sources and exact controls

Primary version: V. Moskowicz and C. Valqui,
*The starred Dixmier conjecture for \(A_1\)*,
https://arxiv.org/abs/1401.5141v1 ,
Proposition 2.6 and its proof, printed pages 4--5.
The rendered page containing the disputed display was inspected.
This excludes a text-extraction or normal-order ambiguity.

The later bibliographic record is *Communications in Algebra*
43(8) (2015), 3073--3082,
https://doi.org/10.1080/00927872.2014.907418 .
That record identifies the publication; its inaccessible full proof
is outside this audit.

The [checker](verify_starred_dixmier_primary_coefficient_audit.py)
derives (C3) symbolically, uses exact Weyl multiplication to verify
(C5)--(C7), and rejects the alleged coefficient identity.
The [source ledger](starred_dixmier_primary_coefficient_audit_source_ledger.json)
pins the proof, checker, inspected primary PDF and separate retained
dependencies reviewed here.

No external human review, formal verification, historical priority,
or counterexample to the starred Dixmier theorem is asserted.

~~~sh
python3 research/dixmier/dixmier-formal-lifting/verify_starred_dixmier_primary_coefficient_audit.py \
  --repo-root . --output /tmp/starred_dixmier_coefficient_review.json
~~~
