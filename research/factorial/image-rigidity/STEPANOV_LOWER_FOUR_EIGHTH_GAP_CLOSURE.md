# Ordinary support {4,12} from the three-stage coefficient obstruction

## Statement

Let `K` be any field of characteristic zero. For every nonzero
`f(X)=X^4(A+B X^8)` in `K[X]` and every integer `n>=1`, the two ordinary
factorial moments `L_1(f^n)` and `L_1(f^(n+1))` cannot both vanish.
Equivalently,

\[
 P_n(t)=\sum_{k=0}^{n}\binom nk(4n+1)_{8k}t^k
\]

satisfies `gcd(P_n,P_(n+1))=1` over `Q[t]` for every `n>=1`.
The coefficient degenerates are included. This is an ordinary factorial
statement and does not assert arbitrary positive Gamma shapes.

The new input is the generic three-stage coefficient lemma in
[STEPANOV_LOWER_THREE_EIGHTH_GAP_CLOSURE.md](STEPANOV_LOWER_THREE_EIGHTH_GAP_CLOSURE.md),
Section 2 (the necessary three-stage coefficient pattern). Its stage class and proof do not depend on lower exponent
three. That lemma has already been accepted after internal root review.
The mixed-offset identity is Theorem `thm:offset-rank` of the retained
[DIAGONAL_PRIME_GAP_SFC.tex](DIAGONAL_PRIME_GAP_SFC.tex).

## 1. An integer-shape offset criterion at gap eight

**Proposition.** Let `n>=1`, let `gamma` be a positive integer, and let
`1<=h<=8`. Write

\[
 T_{n,\gamma,8}(t)=\sum_{k=0}^{n}\binom nk(\gamma)_{8k}t^k,
 \quad C=\gamma+8n+8,
 \quad \kappa_h=\frac{8(n+1)}{(C)_h}.
\]

If `v_2(kappa_h)>=-2`, then `T_(n,gamma,8)` and
`T_(n+1,gamma+h,8)` are coprime over `Q[t]` and over every
characteristic-zero field.

**Proof.** Work in a complete algebraically closed extension of `Q_2`,
with `v_2(2)=1`, valuation ring `O`, and maximal ideal `m`. A common root
is nonzero; set `u=-1/t`. Each block of `8k` consecutive integer factors
contributes at least `7k` to its factorial valuation. The lower monic
reciprocal therefore gives

\[
                              v_2(u)\ge7.
 \tag{1}
\]

Put `c=gamma+7`, `A=c-zD_z+z`, and
`Pi_(n,8)(X)=prod_(ell=0)^(n-1)(X+8ell)`. Keep the shape `gamma` fixed
in every auxiliary. Before cyclic evaluation, define

\[
 \ell_j(X)=(X+8n)\prod_{r=0}^{j-1}(c-X-r),
 \qquad 0\le j\le6.
\]

Then `L_j=Pi_(n,8)(A) ell_j(A)1`. The exact mixed-rank theorem has
columns `L_j+kappa_h*u*v_(7-j)*B_n`, where

\[
 v_k=(C+k+1)_{h-1-k}\quad(0\le k<h),\qquad v_k=0\quad(k\ge h).
\]

A dependence can be normalized to coefficients `b_j in O`, one a
unit. Put

\[
 T=\sum_{j=0}^{6}b_jv_{7-j}\in\mathcal O,
 \qquad Q(X)=\sum_{j=0}^{6}b_j\ell_j(X)+\kappa_huT.
 \tag{2}
\]

It gives a zero cyclic column of `Pi_(n,8)(A)Q(A)1`. The polynomials
`ell_j` have successive degrees `1,...,7` and leading coefficients
`1,-1,...,1`. Thus their nonconstant coefficient matrix is integral
triangular with unit diagonal, and at least one nonconstant coefficient
of `Q` is a unit. By (1) and `v_2(kappa_h)>=-2`, its correction has
valuation at least five. Hence `Q` is integral and primitive.

The operator and seed are in the exact stage class at `q=8`, with
integer `c`, original seed 1, and radius at least seven. The retained
three-stage coefficient lemma applies. Writing `Q=sum_(j=0)^7 c_jX^j`,
it gives

\[
                    v_2(c_0)=4,\qquad v_2(c_1)=2.
 \tag{3}
\]

All coefficients are integral. Consequently `Q(-8n)` has valuation
exactly four: its constant term has valuation four, its linear term
has valuation at least five, and every term of degree at least two has
valuation at least six. On the other hand, each `ell_j(-8n)=0`, so
(2) gives

\[
 Q(-8n)=\kappa_huT,\qquad v_2(\kappa_huT)\ge5.
\]

This contradiction proves coprimality over `Q_2` and hence over `Q`.
A rational Bezout identity transfers the result to every
characteristic-zero field. All coefficients of `Q`, including the
outer scalar `u`, are held constant by the inner Euler operators. ∎

The proof uses the closed radius `v_2(u)>=7`; it does not evaluate the
last scalar-leaf series at a unit argument.

## 2. Application to lower exponent four

Take `gamma=4n+1` and `h=4`. Then

\[
 C=12n+9,\qquad
 \kappa_4=\frac{8(n+1)}{(12n+9)(12n+10)(12n+11)(12n+12)}
          =\frac2{3(12n+9)(12n+10)(12n+11)}.
 \tag{4}
\]

The first and third displayed linear factors are odd, while
`12n+10=2(6n+5)` has valuation exactly one. Thus `v_2(kappa_4)=0`
for every integer `n>=1`, and the proposition applies at every order.

For `AB!=0`, the actual moment normalization is

\[
 \mathcal L_1(f^j)=A^j(4j)!P_j(B/A).
\]

If either coefficient is zero, the phase is a nonzero monomial and
every positive factorial moment is nonzero. This proves the theorem.

## 3. Explicit actual upper row

For an independently checkable offset-four specialization, write
`Q=sum c_jX^j` in (2). Only the last three mixed columns contribute to
`T`, and

\[
 T=b_4+(C+3)b_5+(C+2)(C+3)b_6.
\]

The top three coefficients of the exact `ell_j` basis give

\[
 b_6=c_7,\qquad b_5=-c_6-(16n+33)c_7,
\]
\[
 b_4=c_5+(12n+30)c_6+(144n^2+480n+545)c_7.
\]

Substitution yields the integral functional

\[
                 T=c_5+18c_6+(96n^2+168n+281)c_7.
 \tag{5}
\]

This is the actual mixed row. Its integrality, rather than a presumed
generic rank or an arbitrary cyclic seed, is what makes the valuation
contradiction work.

There are actual lower roots on the radius-seven boundary. If
`n=1 (mod4)`, the coefficient of `t` in `P_n` has valuation exactly
seven; all coefficients of degree `k>=2` have valuation strictly greater
than `7k`. Thus there is exactly one lower root with `v_2(u)=7`.
For example, `n=1` has a single lower root on that boundary. This is
why the closed-radius coefficient obstruction is useful here.

## 4. Exact replay and assurance

Run:

```sh
python3 research/factorial/image-rigidity/verify_stepanov_lower_four_eighth_gap.py
```

The replay reconstructs the complete offset-four upper functional and
checks its identity by polynomial interpolation with an explicit degree
bound: its coefficient functions in `n` have degree at most two, so
the three exact values `n=0,1,2` suffice for that identity. It also
checks the literal kappa normalization, exact integer coefficient
valuation controls, and the radius-seven boundary control. Those
finite coefficient checks supplement the block-counting proof; they
do not substitute for the unrestricted proposition.

The [receipt](stepanov_lower_four_eighth_gap_verification.json) records
the exact replay scope. The root reviewer independently reconstructed
the closed-radius valuation contradiction, the unit value of `kappa_4`,
and the more general `v_2(kappa_h)>=-2` criterion. These are internal
mathematical checks, not external review or formal verification. The
all-order conclusion comes from the proof and its reviewed stage lemma,
not from a finite common-root search.
