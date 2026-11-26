# Strong Factorial windows for homogeneous binary quadratics

Let \(K\) be a field of characteristic zero. Define the factorial functional on \(K[x,y]\) by
\[
 \mathcal L(x^i y^j)=i!j!.
\]
For a polynomial \(f\), write \(\mathcal N(f)\) for its number of nonzero monomials.

**Theorem 1.** For every nonzero homogeneous quadratic
\[
 f=A x^2+Bxy+C y^2\in K[x,y]
\]
and every integer \(n\ge1\), at least one of
\[
 \mathcal L(f^n),\ldots,
 \mathcal L(f^{n+\mathcal N(f)-1})
 \tag{1}
\]
is nonzero. The length three is sharp when all three coefficients are nonzero.

This proves the precise Strong Factorial assertion for homogeneous binary quadratics. In the standard notation \(SFC(m)\), \(m\) denotes the number of variables; the window length is \(\mathcal N(f)\). Thus full \(SFC(2)\) implies Theorem 1. The restricted theorem is not a proof of full \(SFC(2)\), nor does it cover arbitrary inhomogeneous binary quadratics.

## 1. The interval identity

Work first over \(\mathbb C\). Put
\[
 q(u)=A u^2+B u(1-u)+C(1-u)^2,
 \qquad I_m=\int_0^1 q(u)^m\,du .
 \tag{2}
\]
On expanding each monomial and using the beta integral,
\[
 \mathcal L(f^m)=(2m+1)!\,I_m .
 \tag{3}
\]
Indeed the term indexed by \(i+j+k=m\) has integral
\(\int_0^1 u^{2i+j}(1-u)^{j+2k}\,du
 =(2i+j)!(j+2k)!/(2m+1)!\).

Set \(\alpha=A-B+C\). If \(\alpha\ne0\), define
\[
 h=\frac{2C-B}{2\alpha},\qquad
 k=C-\alpha h^2=\frac{4AC-B^2}{4\alpha}.
\]
Then \(q(u)=\alpha(u-h)^2+k\). Integrating the derivative of \((u-h)q(u)^m\) gives, for every \(m\ge1\),
\[
 (2m+1)I_m-2mk I_{m-1}
   =(1-h)A^m+hC^m .
 \tag{4}
\]
This follows because
\[
 \frac{d}{du}\bigl((u-h)q(u)^m\bigr)
  =(2m+1)q(u)^m-2mkq(u)^{m-1}.
\]
The identity is valid for complex coefficients; only integration of polynomials is involved.

If \(\alpha=0\), then \(q(u)=C+(A-C)u\). If \(A=C\), the phase is constant, and its positive moments are nonzero unless \(q=0\). If \(A\ne C\),
\[
 I_m=\frac{A^{m+1}-C^{m+1}}{(m+1)(A-C)}.
 \tag{5}
\]
Two consecutive vanishing moments would give
\(A^{n+1}=C^{n+1}\) and \(A^{n+2}=C^{n+2}\).
If either endpoint is zero both are zero, contradicting \(A\ne C\). Otherwise division gives \(A=C\), again a contradiction. Hence the linear branch has no two consecutive zero moments.

## 2. No three consecutive zeros

Suppose
\(I_n=I_{n+1}=I_{n+2}=0\).
Section 1 excludes \(\alpha=0\). Apply (4) at \(m=n+1,n+2\):
\[
 \begin{split}
 (1-h)A^{n+1}+hC^{n+1}&=0,\\
 (1-h)A^{n+2}+hC^{n+2}&=0.
 \end{split}
 \tag{6}
\]
If \(A,C\) are distinct and nonzero, their two endpoint sequences are linearly independent, so (6) forces \(1-h=h=0\), impossible. If \(A=C\ne0\), its first line is \(A^{n+1}=0\), also impossible.

If \(A=0,C\ne0\), the first line gives \(h=0\). Thus \(B=2C\) and \(q(u)=C(1-u^2)\). But
\(I_n/C^n=\int_0^1(1-u^2)^n\,du>0\).
If \(C=0,A\ne0\), it gives \(h=1\); then \(B=2A\), \(q(u)=Au(2-u)\), and
\(I_n/A^n=\int_0^1[u(2-u)]^n\,du>0\).
Finally if \(A=C=0\), then \(q(u)=Bu(1-u)\), whose positive moments are nonzero whenever \(B\ne0\). If \(B=0\), the phase is the excluded zero polynomial.

These cases exhaust all endpoints. Together with (3), they prove the length-three assertion for every nonzero homogeneous quadratic. The smaller support sizes require the sharper arguments that follow.

## 3. Exact windows for one or two monomials

A monomial \(f=c x^i y^j\) has
\(\mathcal L(f^n)=c^n(ni)!(nj)!\ne0\), so its window has length one.

For adjacent support \(f=A x^2+Bxy\), assume \(A,B\ne0\) and \(I_n=I_{n+1}=0\). The linear branch is excluded by (5). Otherwise \(C=0\), and (4) at \(m=n+1\) yields
\((1-h)A^{n+1}=0\).
Consequently \(h=1\) and \(B=2A\). The positive integral
\(I_n/A^n=\int_0^1[u(2-u)]^n\,du\)
is a contradiction. Swapping \(x,y\), an exact symmetry of \(\mathcal L\), treats the other adjacent support.

For the opposite support \(f=A x^2+C y^2\), define
\[
 P_n(A,C)=\frac{\mathcal L((A x^2+C y^2)^n)}{2^n n!}
  =\sum_{j=0}^n (2j-1)!!(2n-2j-1)!!\,A^jC^{n-j},
 \tag{7}
\]
where \((-1)!!=1\). All coefficients are odd integers. Thus modulo 2,
\[
 \overline P_n(A,C)=h_n(A,C)
       :=\sum_{j=0}^n A^jC^{n-j}.
 \tag{8}
\]
Consecutive \(h_n,h_{n+1}\) have no nonzero projective common zero over \(\overline{\mathbb F}_2\):
\[
 h_{n+1}-A h_n=C^{n+1},
 \tag{9}
\]
and \(C=0\) then forces \(A=0\). Equivalently, after \(C=1\), (9) equals 1, and the point at infinity is excluded by the nonzero leading coefficients.

The integer univariate resultant of \(P_n(T,1)\) and \(P_{n+1}(T,1)\) reduces modulo 2 to a nonzero resultant, because both degrees are retained. It is therefore an odd, nonzero integer. Thus these two polynomials have no common complex root. The leading coefficients also exclude a common projective zero with \(C=0\). This proves the length-two window for the opposite support. No theorem about a general reduction of complex coefficients modulo a prime is being assumed: the integer resultant is reduced.

All nonempty supports are covered, proving Theorem 1 over \(\mathbb C\).

## 4. Characteristic-zero transfer and a sharp witness

The finitely many coefficients of a hypothetical counterexample over \(K\) generate a finitely generated field over \(\mathbb Q\), which embeds into \(\mathbb C\). Its coefficient nonvanishing and finite moment equations survive the embedding, contradicting the complex result. Coefficients need not be algebraic. This proves Theorem 1 over every characteristic-zero field.

Let \(s^2=-15\), and put
\[
 f=(5+s)x^2-20xy+(5-s)y^2.
 \tag{10}
\]
Direct factorial expansion gives
\[
 \mathcal L(f)=\mathcal L(f^2)=0,
 \qquad \mathcal L(f^3)=-115200\ne0.
 \tag{11}
\]
All three coefficients are nonzero. Thus the universal length-two assertion is false for three-term homogeneous quadratics; the support-length-three theorem cannot be uniformly shortened. This witness is not a counterexample to Strong Factorial.

The accompanying standard-library checker verifies (3), the denominator-cleared version of (4), the coefficient identity in (7), and (11) by native polynomial expansion in exact rational arithmetic. It also checks explicit corruptions. These finite controls validate the transcription and witness only; Sections 1--3 are the all-order proof.

## 5. Source comparison and remaining scope

The standard definition is [Edo--van den Essen, *The Strong Factorial Conjecture*, arXiv:1304.3956v2, Definitions 2.7 and Conjecture 2.8](https://arxiv.org/html/1304.3956v2). In that notation the two in \(SFC(2)\) counts variables, while \(\mathcal N(f)\) determines the required window. One must distinguish it from a local abbreviation meaning two monomials.

[Rocks, author preprint, Theorem 3.10](https://people.math.sc.edu/rocks/BradyRocksSFC2015.pdf) proves the Strong Factorial assertion for powers of sums of variable monomials whose exponents are powers of the same prime. Its specialization \(p=2,r=1,m=2\) includes the opposite sparse support in Section 3. Equations (7)--(9) give that specialization directly. The full general theorem is not a proof premise here.

[Liu--Sun, *The factorial conjecture and images of locally nilpotent derivations*, Theorem 2.6, Bull. Aust. Math. Soc. 101 (2020), 71--79](https://doi.org/10.1017/S0004972719000546) proves ordinary Factorial for all homogeneous binary polynomials: vanishing at every positive power forces the phase to be zero. That assertion has a different quantifier from the finite support-size windows proved here.

All external papers in this section provide definitions or comparisons. The proof above uses no deep external nonvanishing theorem. No priority claim is made.

An arbitrary binary polynomial of total degree at most two can have six monomials, including constant and linear terms. For such a polynomial the radial integral mixes its homogeneous layers through different factorial factors; the interval identity (4) does not transfer unchanged. Higher-degree arbitrary binary polynomials and full \(SFC(3)\) are also outside Theorem 1. The witness (10) refutes only the overstrong universal two-moment window for three-term support. Neither this restricted theorem nor the witness establishes a converse to the implication from full \(SFC(2)\).
## 6. Implications and the remaining targets

Here full \(SFC(d)\) denotes the standard assertion for every polynomial in
\(d\) variables over \(\mathbb C\), with its actual nonzero support size as
window length. It is equivalent to its assertion over every characteristic-zero
field by the finite coefficient-field embedding in Section 4.

The proved specialization implications are
\[
 \text{unrestricted SFC}\ \Longrightarrow\ SFC(3)\
 \Longrightarrow\ SFC(2)\
 \Longrightarrow\ \text{all binary polynomials of total degree at most two}\
 \Longrightarrow\ \text{Theorem 1}.
\]
To pass from three variables to two, take a polynomial independent of the third
variable: its support size is unchanged and \(0!=1\) preserves every moment.
Restricting total degree, and then restricting to homogeneous degree two, proves
the other specialization implications. They are valid implications even when
their premises remain open. No reverse specialization is proved.

Equation (3) is a reversible reduction: for each fixed coefficient triple and
finite set of indices, all the factorial moments vanish if and only if all the
corresponding interval moments vanish, since every multiplier \((2m+1)!\)
is nonzero. This exact equivalence concerns the same coefficient triple and
the same indices. It supplies no reduction of inhomogeneous quadratics to
homogeneous ones.

## 7. Exact replay and review scope

The [native checker](verify_binary_homogeneous_quadratic_strong_factorial.py)
has no third-party dependencies. From the repository root run

    python3 research/factorial/image-rigidity/verify_binary_homogeneous_quadratic_strong_factorial.py

The retained [finite-control receipt](binary_homogeneous_quadratic_strong_factorial_verification.json)
records 66 exact transcription controls and 12 rejected corruptions. A fresh
separate execution from outside the repository reproduced that result. The
entire written proof, including the sparse cases, field extension and witness,
was reconstructed by internal reviewers. These checks do not constitute
external specialist review or formal verification. The universal theorem
is established by the written argument, rather than by a finite prefix.
