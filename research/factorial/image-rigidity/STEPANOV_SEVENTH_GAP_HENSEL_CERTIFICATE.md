# An exact certificate for the last lower-three seventh-gap branch

## Scope

This note proves a uniform exclusion on the last residue class left by the
prime-local argument in
[the companion prime-closure note](STEPANOV_LOWER_THREE_PRIME_CLOSURE.md).
It is a written mathematical argument with a small exact certificate, not
an extrapolation from sampled terminal orders. The companion note supplies
the confinement of a common root to this branch; the certificate here
does not assume that all roots of a moment polynomial are units after
scaling.

Let

\[
T_{m,\gamma,7}(t)=\sum_{k=0}^{m}\binom mk(\gamma)_{7k}t^k,
\qquad
H_m(u)=u^mT_{m,3m+1,7}(1/u).
\]

All coefficients are rational integers. Put \(N=n+1\). The remaining
branch has

\[
N\equiv446\pmod{2401},\qquad v_7(u)=1,
\qquad \bar x^2-\bar x+3=0,
\quad x=u/7.
\tag{1}
\]

**Theorem.** No common root of \(H_n,H_{n+1}\) satisfies (1).

The proof uses reduction in an unramified quadratic extension, not an
ordering or positivity assumption on the coefficient field. Coprimality
over \(\mathbb Q\), once combined with the companion confinement theorem,
extends by a rational Bezout identity to every characteristic-zero field.

## 1. Integral normalized moments and a finite truncation

For any nonnegative integer \(m\), define

\[
F_m(z)=T_{m,3m+1,7}(z/7)
=\sum_{k=0}^{m} b_{m,k}z^k,
\qquad
b_{m,k}=\binom mk\frac{(3m+1)_{7k}}{7^k}.
\tag{2}
\]

Each coefficient is an integer. Split the rising product into \(k\)
blocks of length seven and divide each block's unique multiple of seven
by seven. If \(q=\lceil(3m+1)/7\rceil\), the resulting product contains
\((q)_k\), and its other factors are units in \(\mathbb Z_7\).
Consequently

\[
v_7(b_{m,k})\ge\lfloor k/7\rfloor.
\tag{3}
\]

In particular, modulo \(7^4=2401\), every term with \(k\ge28\)
vanishes. Set

\[
f_m(z)=\sum_{k=0}^{27} b_{m,k}z^k\pmod{2401},
\tag{4}
\]

with \(b_{m,k}=0\) if \(k>m\). Evaluation at any integral point in
any extension of \(\mathbb Q_7\) preserves \(F_m\equiv f_m\pmod{7^4}\).

## 2. Exact period in the terminal order

For \(0\le k\le27\),

\[
b_{m+7^5,k}\equiv b_{m,k}\pmod{7^4}.
\tag{5}
\]

To check the binomial factor, use Vandermonde's identity and, for
\(1\le j\le27\),

\[
v_7\binom{7^5}{j}=5-v_7(j)\ge4.
\tag{6}
\]

For completeness, (6) follows from
\(\binom{7^5}{j}=(7^5/j)\binom{7^5-1}{j-1}\): the second factor
is a unit, since \(7^5-r\) and \(r\) have the same valuation for
\(1\le r<7^5\). Thus the binomial factor is periodic modulo \(7^4\).

In each rising-product block used above, the position of the multiple of
seven is unchanged by \(m\mapsto m+7^5\). After division by seven,
that factor changes by \(3\cdot7^4\); every other factor changes by
\(3\cdot7^5\). Hence the normalized product is also unchanged modulo
\(7^4\). This proves (5), including at the ends of the binomial sum.

It follows that every order in (1) is covered by the following seven
representatives modulo \(7^5=16807\):

\[
N_i=446+2401i,\qquad 0\le i\le6.
\tag{7}
\]

This period theorem, rather than a numerical range estimate, is what
makes the finite certificate uniform in \(N\).

## 3. The quadratic extension and root uniqueness

Let \(\theta\) satisfy

\[
\theta^2-\theta+3=0.
\]

Its discriminant \(-11\) is \(3\) modulo seven, a nonsquare. Thus
\(E=\mathbb Q_7(\theta)\) is the unramified quadratic extension and
its integer ring is \(\mathbb Z_7[\theta]\). Work modulo \(2401\)
in

\[
R=(\mathbb Z/2401\mathbb Z)[\theta]/(\theta^2-\theta+3).
\tag{8}
\]

Write elements as pairs \((a,b)\), meaning \(a+b\theta\). Their
norm is \(a^2+ab+3b^2\), so an element is a unit precisely when its
reduction modulo seven is nonzero. In (1), choose the residue embedding
with \(\bar x=\bar\theta\). The other embedding is its conjugate.
The variable in (2) is

\[
z=7/u=1/x,\qquad
\bar z=\bar\theta^{-1}=(1-\bar\theta)/3=5+2\bar\theta.
\tag{9}
\]

The lower normalized moment has a simple zero at (9) modulo seven.
This is verified directly by its degree-at-most-27 reduction in the
certificate; equivalently the lower unit-root polynomial is

\[
\ell_7(x)=x^4+6x^3+x^2+2x+1
=(x-3)(x+3)(x^2-x+3)\quad\text{in }\mathbb F_7[x],
\]

and these three factors are pairwise coprime and squarefree.

Hensel's uniqueness theorem therefore applies in the residue disk (9).
One can also see the needed finite uniqueness directly: for two points
in that disk, the divided difference of the lower polynomial is a unit,
since its reduction is the nonzero derivative at (9). Consequently a
point where \(F_{N-1}\) vanishes must agree modulo \(7^4\) with the
displayed lift below. This argument applies to a point in an arbitrarily
ramified extension as well; it does not require that a hypothetical
common root was initially known to belong to \(E\).

## 4. The complete certificate

The table gives \(z_i=a_i+b_i\theta\) modulo \(2401\) and the pair
\((A_i,B_i)\) such that

\[
f_{N_i-1}(z_i)=0,\qquad
f_{N_i}(z_i)=343(A_i+B_i\theta)\pmod{2401}.
\tag{10}
\]

Every lower derivative is a unit, and every \(z_i\) reduces to
\(5+2\theta\) modulo seven.

| \(i\) | \(N_i\) | \(a_i\) | \(b_i\) | \(A_i\) | \(B_i\) |
|---:|---:|---:|---:|---:|---:|
| 0 | 446 | 1517 | 2053 | 4 | 6 |
| 1 | 2847 | 1860 | 1710 | 0 | 4 |
| 2 | 5248 | 2203 | 1367 | 3 | 2 |
| 3 | 7649 | 145 | 1024 | 6 | 0 |
| 4 | 10050 | 488 | 681 | 2 | 5 |
| 5 | 12451 | 831 | 338 | 5 | 3 |
| 6 | 14852 | 1174 | 2396 | 1 | 1 |

None of the pairs \((A_i,B_i)\) is zero in \(\mathbb F_7^2\).
The upper moment therefore cannot vanish at the unique lower lift.
The conjugate residue disk gives the conjugate nonzero value. Equations
(3)--(7) cover every terminal order in (1), proving the theorem.

All 28 coefficients of each of the fourteen polynomials, along with the
seven lifts, are retained in
[the exact witness](stepanov_seventh_gap_hensel_witness.json).
The [standalone verifier](verify_stepanov_seventh_gap_hensel.py)
reconstructs every coefficient by direct integer multiplication, checks
all displayed identities, derivative units, residue embeddings, and
nonzero upper values. It also rejects deliberately corrupted witnesses.
The [verification receipt](stepanov_seventh_gap_hensel_verification.json)
records the replay and source hashes.

The uniform proof consists of the valuation truncation, the period
identity, and the simple-root argument. The finite arithmetic is only
the seven exact identities (10).

## 5. Scope and failed shortcuts

- The test is not an assertion that all roots are unramified. The
  companion note first confines common roots to the simple quadratic
  unit branch; uniqueness then identifies its lift.
- A common residual quadratic modulo seven does not itself exclude a
  characteristic-zero common factor. The nonzero value modulo \(7^4\)
  is essential.
- Checking only \(N=446\) would be insufficient: coefficients modulo
  \(7^4\) can change as \(N\) increases by \(2401\). All seven
  classes modulo \(16807\) are included.
- Integer truncation at degree 27 is justified by (3), not a numerical
  approximation of a convergent infinite series.
- The calculation proves the stated two-support branch, not SFC for
  arbitrary support or even all ordinary binomials.

No claims about Alpöge's RH argument are used. This certificate has no
analytic input, and does not assert publication, external refereeing,
formal verification, or a literature priority claim.
