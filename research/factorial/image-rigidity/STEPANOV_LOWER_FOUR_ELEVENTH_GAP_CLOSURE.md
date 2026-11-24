# Lower exponent four: closing the exceptional prime gap eleven

Internal written proof, 21 November 2026. Inspected repository base:
`19eb9912fb60142f8b620dac24896931d4ff3e46`.

## 1. Result

**Theorem.** Over every characteristic-zero field, every nonzero polynomial

\[
f(X)=X^4(A+BX^{11})
\]

has a nonzero ordinary factorial moment in every pair of consecutive positive
orders. Consequently, combining the retained prime-gap families, **lower
exponent four satisfies ordinary SFC at every prime gap**.

This closes the exceptional prime in the retained lower-four theorem. It does
not close all composite gaps or unrestricted SFC. The proof uses an integral
auxiliary scalar and an exact finite-field obstruction. It assumes no ordering
or positivity on the arbitrary coefficient field.

## 2. The actual pair and its remaining class

Put \(N=n+1\), \(\gamma=4n+1\), and

\[
T_{n,\gamma,11}(t)=\sum_{k=0}^n\binom nk(\gamma)_{11k}t^k,
\qquad P_n(u)=u^nT_{n,\gamma,11}(1/u).
\]

The upper moment uses \(T_{n+1,\gamma+4,11}\), with its actual shifted
shape. The ordinary phase has moments
\(B^m(4m)!P_m(A/B)\) when \(B\ne0\). Monomial degenerates have nonzero
moments automatically.

Let

\[
C=\gamma+11n+11=15N-3,
\qquad
\kappa=\frac{11N}{(C)_4}
=\frac{11}{15C(C+1)(C+2)}.
\tag{1}
\]

The retained native-prime criterion and the extra ordinary digit proved in
[STEPANOV_NATIVE_PRIME_OFFSET_COROLLARIES.md](STEPANOV_NATIVE_PRIME_OFFSET_COROLLARIES.md)
show that a common root can occur only if

\[
11^3\mid C,
\qquad N\equiv1065\pmod{1331}.
\tag{2}
\]

Indeed the possible denominator residues correspond to
\(\sigma_4(0)=11\), \(\sigma_4(-1)=-1\), and
\(\sigma_4(-2)=-1\). Only the first is zero modulo eleven.
In (2), \(n\ge2\), \(\gamma\equiv0\pmod{11}\), and the factors
\(15,C+1,C+2\) are units.

Every lower reciprocal root satisfies \(v_{11}(u)\ge1\), because each
coefficient of degree \(k\) in \(T\) contains \(k\) complete eleven-factor
blocks. This is an integer-shape assertion, not a claim about arbitrary
ramified integral Gamma parameters.

## 3. A fixed unit minor

Use the actual auxiliaries at the same lower shape:

\[
B_0=1,\quad B_{j+1}=(z+11j+\gamma+10)B_j-zB_j',
\qquad L_j=z^j(D_z-1)^jB_{n+1}.
\]

Let \(Z(u)\) have cyclic columns \(-L_{10-m}\) for
\(1\le m\le10\), and last column \(B_n\), in the basis
\(1,z,\ldots,z^{10}\) modulo \(z^{11}-u\). The retained full-jet
concomitant, with its nonzero-jet argument, gives a kernel vector at every
nonzero lower root. The precise sources are
[DIAGONAL_PRIME_GAP_SFC.tex](DIAGONAL_PRIME_GAP_SFC.tex), `thm:offset-rank`,
[LOWER_TWO_PRIME_GAP_SFC.tex](LOWER_TWO_PRIME_GAP_SFC.tex), `thm:mixed`,
and [HIGHER_GAP_ADJOINT_RANK.tex](HIGHER_GAP_ADJOINT_RANK.tex), `thm:rank`.

Delete row zero and column eleven to form a ten-by-ten minor. Its constant
determinant is an eleven-adic unit. Here is the proof specialized to the
present residue, rather than an appeal to generic rank.

At \(u=0\), write \(A_0=z-zD_z\),
\(\Pi_n(X)=\prod_{r=0}^{n-1}(X+\gamma+10+11r)\). The cyclic basis
map \(g(X)\mapsto g(A_0)1\) is unimodular and identifies the quotient
with \(\mathbb Z_{11}[X]/\prod_{r=0}^{10}(X+r)\). Also

\[
B_n=\Pi_n(A_0)1,\qquad
L_j=(-1)^j(A_0)_j(A_0+C-1)B_n.
\]

For a prospective homogeneous null vector with last component zero, put

\[
F(X)=w_{11}-(X+C-1)\sum_{j=0}^{9}(-1)^jw_{10-j}(X)_j.
\tag{3}
\]

The ten retained rows say \(F\Pi_n=\epsilon\), modulo the characteristic
polynomial, for a constant \(\epsilon\). Modulo eleven,
\(\Pi_n=(X-1)^n\). Evaluation at one gives \(\epsilon=0\);
evaluation at the other ten residues gives zeros of \(F\) there. Since
\(F(1)=w_{11}=0\) and \(\deg F\le10\), \(F=0\) and the whole
homogeneous vector vanishes. This proves invertibility of the fixed minor.

It follows that the normalized vector
\(w(u)=(w_1(u),\ldots,w_{10}(u),1)\) belongs to
\(\mathbb Z_{11}[[u]]^{11}\), converging for \(v_{11}(u)>0\).
At an actual lower root it is the actual jet divided by its nonzero last
component. No lower-root simplicity is assumed here.

As in [the lower-three prime proof](STEPANOV_LOWER_THREE_PRIME_CLOSURE.md),
the characteristic-zero Euler operator is used before cyclic projection,
or at \(u=0\). It does not act on the nonzero characteristic-zero cyclic
quotient. Variable-\(u\) cyclic operator calculations below are modulo eleven.

## 4. The two uniform scalar coefficients

Define

\[
S(u)=(C+2)(C+3)w_1(u)+(C+3)w_2(u)+w_3(u).
\tag{4}
\]

The exact offset-four upper relation at a lower root is

\[
1=\kappa uS(u).
\tag{5}
\]

This is the rising-basis expansion of
\(((E)_4-(C)_4)/(E-C)\), with its constant term dropped at a lower
zero. We claim

\[
S_0\equiv55\pmod{121},\qquad S_1\equiv5\pmod{11}.
\tag{6}
\]

At \(u=0\), the exceptional eigenvalue is \(-10\); its value under
\(\Pi_n\) is divisible by \(11^n\), while all other eigenvalues give
units. The same unit Vandermonde interpolation used above yields

\[
F(X)\equiv\frac{(X)_{10}}{(1-C)_{10}}\pmod{11^n}.
\tag{7}
\]

The denominator is a unit. If \(f_{10},f_9,f_8\) denote the top three
coefficients of \(F\), expansion of (3) gives exactly

\[
S=(3C^2+69C+546)f_{10}-(2C+30)f_9+f_8.
\tag{8}
\]

The top three coefficients of \((X)_{10}\) are \(1,45,870\).
Therefore

\[
S_0\equiv\frac{3(C^2-7C+22)}{(1-C)_{10}}\pmod{11^n}.
\tag{9}
\]

When \(121\mid C\), the numerator is \(66\) modulo \(121\),
and the denominator is \(-1\) modulo eleven. Its constant factor gives
\(S_0/11\equiv-6=5\pmod{11}\), proving the first part of (6).

For the second part, work in \(\mathbb F_{11}[u]/(u^2)\). The cyclic
characteristic polynomial is \(X^{11}-X-u\), with exceptional root
\(\rho=1-u\). Since \(n\ge2\), \((\rho-1)^n=0\). Thus

\[
F(X)=-H_u(X),\qquad H_u(X)=\frac{X^{11}-X-u}{X-\rho},
\]

where the monic polynomial quotient has \(H_u(1)=-1\).
Its top three coefficients are \(1,\rho,\rho^2\). Formula (8) modulo
eleven becomes

\[
S(u)=-7+8\rho-\rho^2=5u\pmod{u^2}.
\]

This proves (6) and, for every integral \(x\),

\[
S(11x)/11\equiv5(1+x)\pmod{\mathfrak m}.
\tag{10}
\]

## 5. The lower residue polynomial

Under (2), \(n\equiv8\pmod{11}\). Let
\(q=\gamma/11=C/11-N\), so \(q\equiv-N\pmod{11}\).
Wilson's product in each complete block gives
\(11^{-k}(\gamma)_{11k}\equiv(-1)^k(q)_k\pmod{11}\).
Hence the monic integral polynomial
\(11^{-n}P_n(11x)\) reduces to \(x^{n-8}\ell(x)\), where

\[
\begin{aligned}
\ell(x)
&=\sum_{k=0}^{8}\binom8k9^{\underline k}x^{8-k}\\
&=x^8+6x^7+3x^6+9x^5+7x^4+6x^3+x^2+4x+1\quad\text{in }\mathbb F_{11}[x].
\end{aligned}
\tag{11}
\]

The polynomial is squarefree and \(\ell(-1)=10\ne0\). Therefore every
actual lower root, including roots with \(v_{11}(u)>1\), has
\(v_{11}(S(u))=1\) by (10).

Put \(t=v_{11}(C)\ge3\). Equation (5) forces

\[
v_{11}(u)=t-2.
\tag{12}
\]

## 6. Excluding both valuation ranges

If \(t\ge4\), every common root lies in \(v_{11}(u)>1\). The divided
difference of \(uS(u)\) at two distinct points of this disk is

\[
S_0+(u+v)S_1+\sum_{j\ge2}S_j\sum_{r=0}^{j}u^{j-r}v^r.
\]

Its first term has valuation one and every other term has larger valuation.
It cannot vanish. Thus (5) permits at most one distinct common root. Rational
conjugacy and characteristic-zero separability force any common irreducible
factor over \(\mathbb Q\) to be linear. Such a root is impossible by the
exact factorial sum-of-rational-squares identity in Section 7 of
[the lower-three prime proof](STEPANOV_LOWER_THREE_PRIME_CLOSURE.md): one of
the consecutive powers of \(X^4(u+X^{11})\) is an even power of a nonzero
rational polynomial.

If \(t=3\), write \(C=11^3d\) with \(d\) a unit. Every common root
has \(u=11x\) with \(x\) a unit. Reduction of (5) gives

\[
x^2+x-6\bar d=0.
\tag{13}
\]

The complete finite-field check is

\[
\gcd_{\mathbb F_{11}[x]}\bigl(\ell(x),x^2+x-6d\bigr)=1
\qquad(d=1,2,\ldots,10).
\tag{14}
\]

These ten exact gcd identities, retained with polynomial Bezout witnesses in
the accompanying receipt, rule out the remaining range. Equations (2) and
(12)--(14) prove the theorem for every terminal order. They are not a
bounded terminal scan.

Rational coprimality supplies a Bezout identity over \(\mathbb Q\), which
remains valid over every characteristic-zero coefficient field.

## 7. All prime gaps at lower exponent four

The retained `cor:lower-four-prime` in
[PRIME_PLUS_ONE_EXPONENT_GAP_SFC.tex](PRIME_PLUS_ONE_EXPONENT_GAP_SFC.tex)
already covers every prime gap greater than four except possibly eleven.
The present theorem supplies eleven. Gap three is that note's prime-plus-one
family. Gap two is covered by
[LOWER_FOUR_QUADRATIC_SECTOR_SFC.tex](LOWER_FOUR_QUADRATIC_SECTOR_SFC.tex).
Thus every prime gap is now covered, without asserting coverage of every
composite gap.

## 8. Verification and limits

The universal input is the fixed-minor interpolation, the first two scalar
coefficients, the root-radius argument, and the injective divided difference.
The finite arithmetic obligation is precisely (14).

[The exact checker](verify_stepanov_lower_four_eleventh_gap.py) verifies
the coefficient identity (8), (9), the dual-number calculation, all ten
Bezout witnesses, the lower reduction, and actual normalized-jet expansions
at representative indices. The latter are transcription controls, not a proof
of the all-order assertions. Results are retained in
[the receipt](stepanov_lower_four_eleventh_gap_verification.json).

The argument retains the actual shifted shapes and seed. It does not assume
generic rank, does not import an unrestricted positivity claim, and has no
Riemann-hypothesis premise. Assurance is internal written review and exact
replay, not external refereeing, formal verification, or a priority claim.
