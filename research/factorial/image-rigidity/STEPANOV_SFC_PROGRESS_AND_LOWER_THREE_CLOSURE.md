# SFC progress: complete lower exponent three and further transfer families

Internal research synthesis, 21 November 2026. Latest inspected GitHub main:
`19eb9912fb60142f8b620dac24896931d4ff3e46`.

The active research direction is now
[structural closure](SFC_STRUCTURAL_CLOSURE.md): bounds uniform in all
exponents, exact determinant stability, and prime-power scalar transfer.
The completed two-support families below remain available as checks.

## 1. What is now proved

The new proofs close the entire ordinary univariate family

\[
\boxed{f(X)=X^3(A+BX^d),\qquad d\ge1.}
\tag{1}
\]

For every characteristic-zero coefficient field, every nonzero phase (1), and
every \(n\ge1\), the moments \(\mathcal L(f^n)\) and
\(\mathcal L(f^{n+1})\) cannot both vanish. Here
\(\mathcal L(X^j)=j!\). Together with the retained lower exponents zero,
one, and two, this covers **every ordinary univariate binomial whose least
exponent is at most three**.

There are three additional advances beyond this complete layer:

1. Every support \(\{3^r,3^{r+1}\}\), \(r\ge1\), is covered at every
   starting order. This includes \(\{9,27\}\), \(\{27,81\}\), and
   arbitrarily large exponents.
2. Lower exponent four is now covered at **every prime gap**, including the
   previously exceptional gap eleven. Gap eight is also closed.
3. Explicit multivariate families with arbitrarily many monomials satisfy
   all-order SFC, using polynomial inverse rigidity and an exact Gamma-to-
   ordinary block-sum identity.

**Unrestricted SFC remains open in this investigation.** These are complete
theorems on specified families, not a proof for every support. No full ordinary
binomial theorem or general multivariate two-support theorem is claimed.

## 2. The complete lower-three proof

For ordinary binomials define

\[
H_{m;a,d}(u)=\sum_{k=0}^{m}\binom mk(am+1)_{dk}u^{m-k}.
\tag{2}
\]

These are monic integer polynomials. If \(B\ne0\),

\[
\mathcal L([X^a(A+BX^d)]^m)=B^m(am)!H_{m;a,d}(A/B).
\tag{3}
\]

The two actual shapes are \(an+1\) and \(a(n+1)+1\). All proofs below
retain that change. A rational Bezout identity for the two polynomials in
(2) remains valid over every characteristic-zero field. If a coefficient is
zero, the phase is a monomial and all its positive moments are nonzero.

For \(a=3\), the following disjoint cases cover every positive gap:

| Gap \(d\) | Proof supplying every starting order | Status in this investigation |
|---:|---|---|
| 1 | [UNIFORM_PADE_PARABOLA_SFC.tex](UNIFORM_PADE_PARABOLA_SFC.tex), `thm:main` | Retained theorem |
| 2 | [PRIME_PLUS_ONE_EXPONENT_GAP_SFC.tex](PRIME_PLUS_ONE_EXPONENT_GAP_SFC.tex), `thm:family`, including its explicit prime-two input | Retained theorem |
| 3 | [DIAGONAL_PRIME_GAP_SFC.tex](DIAGONAL_PRIME_GAP_SFC.tex), `thm:diagonal-prime` | Retained theorem |
| 4 | [LOWER_THREE_FOURTH_POWER_FROBENIUS_SFC.tex](LOWER_THREE_FOURTH_POWER_FROBENIUS_SFC.tex), `thm:coprime` | Retained theorem |
| 5 | [STEPANOV_LOWER_THREE_PRIME_CLOSURE.md](STEPANOV_LOWER_THREE_PRIME_CLOSURE.md), Theorem A | New all-order closure |
| 6 | [STEPANOV_ADJACENT_POWERS_OF_THREE_SFC.md](STEPANOV_ADJACENT_POWERS_OF_THREE_SFC.md), at \(r=1\) | New all-order closure |
| 7 | [STEPANOV_LOWER_THREE_PRIME_CLOSURE.md](STEPANOV_LOWER_THREE_PRIME_CLOSURE.md), Theorem B, and [the complete Hensel certificate](STEPANOV_SEVENTH_GAP_HENSEL_CERTIFICATE.md) | New all-order closure |
| 8 | [STEPANOV_LOWER_THREE_EIGHTH_GAP_CLOSURE.md](STEPANOV_LOWER_THREE_EIGHTH_GAP_CLOSURE.md) | New all-order closure |
| \(d\ge9\) | [RELATIVE_GAP_THREE_STRONG_FACTORIAL.md](RELATIVE_GAP_THREE_STRONG_FACTORIAL.md), with \(a=3\) | Retained theorem |

This proves (1). The relative-gap input includes its retained analytic tail
and complete finite certificate; the synthesis does not replace that input
with a new proof. The new gaps five through eight are independent of the
support-uniform eventual theorem and its unevaluated cutoff. Their closure
does not depend on increasing the finite terminal prefix beyond twelve.

## 3. The mechanisms that closed the four gaps

### 3.1. Prime gaps five and seven

A fixed minor of the **actual** concomitant normalizes the actual rising jet.
For \(p=5,7\), its scalar \(S(u)\) is an integral convergent power series
satisfying

\[
S_0\equiv p\pmod{p^2},\qquad S_1\equiv-1\pmod p.
\]

After the preceding \(p^3\) reduction, every common root satisfies
\(1=\kappa uS(u)\). On the disk \(v_p(u)>1\), the divided difference
of \(uS(u)\) has valuation exactly one, so there is at most one common
root. Rational conjugacy and an exact factorial sum-of-rational-squares
identity exclude that possibility.

The boundary disk has one residual candidate at gap five, so the same
one-root argument closes it. At gap seven, only a simple residual quadratic
remains, at \(N\equiv446\pmod{2401}\). The normalized moments modulo
\(7^4\) truncate at degree 27 and have period \(7^5\) in \(N\).
Seven exact Hensel witnesses therefore cover that entire infinite class.
Every upper evaluation is nonzero. Two separately written verifiers replay
those identities, including a two-by-two matrix representation of the
unramified quadratic ring.

### 3.2. Gap six and all adjacent powers of three

For \(a=3^r\), \(d=2a\), define

\[
P_n(x)=3^{-an}H_{n;a,2a}(3^a x).
\]

Exact factorial splitting makes \(P_n\) an integer polynomial with unit
constant term at three, and gives the coefficientwise identity

\[
P_n(x)\equiv Q_n(x):=
\sum_{k=0}^{n}(-1)^k\binom{n+2k}{3k}x^{n-k}\pmod9.
\tag{4}
\]

The mod-three reductions of \(Q_n,Q_{n+1}\) can have large common
factors. They are not claimed coprime. Instead their rational generating
function proves that **every common residual root is simple in the lower
polynomial**. Hensel uniqueness puts any hypothetical actual common root
in a finite unramified extension of \(\mathbb Q_3\).

For roots \(\tau_i\) of \(\tau^3+\tau^2-x\), let \(S_j\) and
\(V_j\) be the power sums of the roots and of their pairwise products.
The exact identities

\[
Q_N+2xQ_{N-1}=S_{3N}
=S_N^3-3S_NV_N+3x^N
\tag{5}
\]

give the contradiction. A common root and (4) force
\(S_{3N}\in9\mathcal O\). In the unramified ring, reduction of (5)
forces \(S_N\in3\mathcal O\); then (5) gives
\(S_{3N}\equiv3x^N\not\equiv0\pmod9\), since \(x\) is a unit.
This proof has no positivity input.

### 3.3. Gap eight

Three exact Artin–Hasse splittings, with their actual descended seeds,
give seven nonresonant scalar equations. Interpolation at seven integer
nodes has an invertible matrix \(I+2vU(v)\), and determines the monic
relation polynomial with a controlled \(8v\) error. In each remaining
valuation range, the upper mixed condition becomes a restricted power
series with unit linear coefficient and all higher coefficients even.
It is injective on the entire relevant disk. Rational conjugacy and the
same explicit rational Gram identity exclude its last possible common root.

The proof checks the last split to scalar equations directly. It never
assumes that a characteristic-zero Euler operator acts on a quotient by a
nonzero cyclic parameter.

## 4. Lower exponent four

[The gap-eleven theorem](STEPANOV_LOWER_FOUR_ELEVENTH_GAP_CLOSURE.md)
uses the same fixed-minor mechanism, now with

\[
S_0/11\equiv5,\qquad S_1\equiv5\pmod{11}.
\]

The deep disk is again injective. In the boundary disk the lower residue
polynomial is coprime to \(x^2+x-6d\) for all ten nonzero
\(d\in\mathbb F_{11}\). These are ten exact finite-field Bezout
certificates, not ten sampled terminal orders. This supplies the one missing
prime in the retained lower-four prime-gap theorem.

[The gap-eight theorem](STEPANOV_LOWER_FOUR_EIGHTH_GAP_CLOSURE.md)
uses the generic coefficient valuations from the three-stage argument.
The monic relation has constant valuation four; the actual upper correction
has valuation at least seven. The contradiction holds also at the boundary
root radius, so no residue progression is left.

After the retained degree, dyadic, prime-gap and relative-gap results, a
useful containing frontier for completing lower exponent four is

\[
d\in\{6,9,10\},\qquad N\ge13.
\tag{6}
\]

Here gap four is covered by the retained
[dyadic theorem](DYADIC_BINOMIAL_SFC.tex), gap two by the retained
[lower-four sector theorem](LOWER_FOUR_QUADRATIC_SECTOR_SFC.tex), and
\(d\ge12\) by relative-gap-three. Other retained arithmetic filters can
further reduce (6); it is a containing target set, not a claim that every
listed pair defeats every available criterion.

## 5. Progress toward unrestricted support

The [common-power transfer](POWER_DIAGONAL_LINEAR_SFC_TRANSFER.md) proves
all-order SFC for

\[
f(Y)=(Y_1\cdots Y_r)^q\sum_i c_iY_i,\qquad q\ge1.
\tag{7}
\]

It identifies the moments exactly with the inverse coefficients of
\(z\prod_i(1-c_i z^q)\). Structural zero coefficient classes turn a
support-size moment window into a long consecutive inverse-coefficient
window, excluded by the retained polynomial inverse-rigidity theorem.

The [weighted radial extension](WEIGHTED_RADIAL_LINEAR_SFC_TRANSFER.md)
strengthens the critical-value argument by excluding the value zero.
For arbitrary positive integers \(p_i\), take independent blocks
\(S_i=\sum_{j=1}^{p_i}X_{ij}\). Then

\[
G(X)=\prod_i S_i^{q p_i}\sum_i c_iS_i
\tag{8}
\]

has a nonzero ordinary factorial moment in every window of length equal
to the number of **distinct nonzero coefficient values** among the \(c_i\).
This is at most its actual number of monomials. The exact block-sum identity
\(\mathcal L(S_i^k)=(p_i)_k\) supplies the ordinary transfer.

The input from the credited critical-value proof has been internally
reconstructed. The extension does not assume arbitrary complex moments are
positive. It does not delete Gamma weights or identify independent variables
without proving the corresponding moment identity.

## 6. Failed routes and remaining mathematical gap

The successful arguments do not validate the generic determinant proposals
previously rejected. Pure Hankel determinants can vanish at nonzero
binomials; full shifted-jet determinants reproduce the unknown resultant;
product-Wronskians alone give no improved common-factor bound. Those exact
countercontrols remain in
[the original investigation](STEPANOV_AUXILIARY_DETERMINANTS_SFC2.md).

The new gap-six proof also records a useful failure: adjacent residual gcds
modulo three can have degree ten already at lower order 39, and degree
thirteen at 52. A uniformly small residual-root count was false. Simple
residual roots and the next congruence digit supplied the missing invariant.

For the ordinary unequal-exponent phase
\(F=(\prod_iY_i^{r_i})\sum_i c_iY_i\), there is now a precise further
target. Set

\[
a(z)=z\prod_i(1-c_i z)^{r_i},\qquad
H'(z)=\prod_i(1-c_i z)^{r_i-1},\qquad H(0)=0.
\]

Direct coefficient extraction gives

\[
\mathcal L(F^m)=m!\prod_i(r_i m)!(m+1)
 [t^{m+1}]H(a^{-1}(t)).
\tag{9}
\]

The desired contact bound is
\(\operatorname{ord}_0(H-b\circ a)\le\deg b+s\), for the matching
special zero-fiber profile and active support count \(s\). It is not
ordinary inverse rigidity: arbitrary polynomial observables admit arbitrarily
long gaps. This bound is a genuine remaining theorem to prove, not a corollary
already established here. Even this special family is narrower than arbitrary
integer support matrices, so it would still not by itself prove unrestricted
SFC.

Two exact countercontrols explain this distinction. With
\(a=z(1-z)^2\) and \(H_M=a+z^M\), one has
\(H_M(a^{-1}(t))=t+t^M+O(t^{M+1})\), allowing arbitrarily long gaps
for an unrestricted observable. Even requiring \(a'=H'D\) with a small
\(\deg D\) is insufficient: \(a=z-z^2-z^3\),
\(H=z-3z^2/2\), \(D=1+z\) satisfies that identity, but the cubic
coefficient of \(H(a^{-1})\) is zero. The special zero-fiber profile in
(9) must supply additional information.

## 7. Reproducibility and assurance

Every new theorem has a written proof, exact algebraic controls, and a retained
verification receipt. The cubic family, prime-gap calculations, and scalar-leaf
argument received separate internal review. The seventh-gap branch has two
algorithmically distinct certificate replays.

The [package manifest](stepanov_investigation_manifest.json) lists all new files,
their hashes, proof scope, and the retained baseline source-guard failures.
The [original investigation](STEPANOV_AUXILIARY_DETERMINANTS_SFC2.md), Section 9,
records those failures explicitly. No old source pin has been changed to
manufacture a successful replay. This synthesis does not depend on the failing
support-uniform eventual-source chain.

Assurance remains internal written mathematics and exact reproducible
arithmetic. No claim of external peer review, formal verification, publication,
or literature priority is made. The global registry and generated navigation
have not been promoted by this investigation.
