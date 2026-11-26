# Fixed weighted shifts of the least target multipliers

Let \(n=X^L\), with \(L\) fixed before the height limit,
and let \(M=\lceil K/L\rceil\). For a target \(c\le X\)
of a genuine positive primitive row, write
\[
 Z_M(c)=\prod_{p\mid c}p^{(M-v_p(c))_+},\qquad
 H_M(c)=cZ_M(c),\qquad \gamma_\pm=1/3\pm\xi .
\]
Assume the original target cap \(\operatorname{rad}(c)\le
n^{\gamma_+}\). The original caps on both summands and
any additional least-repair eligibility can be retained.
All multipliers below are the actual least multipliers.
[Least target repair (MPR)](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

## A strict half bound at a prescribed weighted multiplier shift

**LM1.** There is an absolute \(\Delta>0\) such that for
sufficiently small fixed \(\xi>0\), sufficiently large
fixed \(K\), and sufficiently small fixed \(L>0\), the
number of actual target pairs \(c_1,c_2\le X\) with
\[
                   A Z_M(c_1)-B Z_M(c_2)=d\ne0
\]
is \(O_{K,L,\varepsilon}(n^{1/6-\Delta}X^\varepsilon)\).
Here \(A,B\) are prescribed positive integers at most
\(X^{C_0}\) for any fixed \(C_0\). The bound is uniform
in their values and the nonzero signed shift, with the
fixed \(C_0\) suppressed in the constant and in the
choice of \(K\).

Put \(Y_0=X^{M L\gamma_+}\), so \(Z_M(c_i)\le Y_0\).
Fix \(0<\theta<1\) and let \(k_0=M-\lfloor\theta M\rfloor\).
On each original target the primes with \(v_p(c_i)>\theta M\)
have radical at most \(X^{1/(\theta M)}\le n^{1/(\theta K)}\).
Fix the whole corresponding supported parts \(E_i\) of
\(Z_M(c_i)\) at this price, up to a height-subpower factor.
The free integers \(z_i=Z_M(c_i)/E_i\) are \(k_0\)-full
and at most \(Y_0\). The original equation becomes
\((AE_1)z_1-(BE_2)z_2=d\).

To make this a primitive fixed-shift equation, choose
the gcd of its two complete terms among the divisors of
\(|d|\). Fix the full \(z_i\)-parts on the primes of \(|d|\)
at uniform subpower cost. The chosen gcd divides both
supported coefficient parts, and dividing through preserves
the free floors on both sides. The resulting nonzero shift
and coefficients have bounded polynomial height in \(Y_0\).
Indeed one can choose \(K\) large enough in terms of \(C_0\)
that \(X^{C_0}\le Y_0\). All supported coefficients before
division are at most the corresponding complete term,
which is at most \(Y_0^2\). SHF1 in its primitive
arbitrary-coefficient form therefore supplies an absolute
\(\delta_0>0\) and an upper exponent
\((1/2-\delta_0)/k_0+\varepsilon'\) in base \(Y_0\),
at all sufficiently large fixed \(k_0\).
[Primitive coefficient form of SHF1](ABC_FIXED_SHIFT_STRICT_HALF_SAVING.md#why-equality-of-the-raw-and-full-height-budgets-is-impossible)

After fixing \(Z_M(c_i)\), its prime powers determine all the
original valuations below \(M\). The remaining original primes
have valuations at least \(M\), so their support costs at most
\(n^{1/K}\) per target, with exact exponents at height-subpower
cost. Using \(M/k_0\le1/(1-\theta)\) gives a total count exponent,
before height-subpower factors, at most
\[
                 \frac{2}{\theta K}+\frac2K+
                 \frac{(1/2-\delta_0)\gamma_+}{1-\theta}.
\]
Choose \(\theta,\xi\) small first in terms of \(\delta_0\),
then \(K\) large and then \(L\) sufficiently small that
\(k_0\) is above the fixed-shift threshold. This makes the
display at most \(1/6-\Delta\). All constants are fixed
during the height limit. This proves the claimed count on
primitive targets without any choice of common padding.

## A mixed original and repaired pair record

**LM2.** Under the conditions of LM1, prescribe the exact
reduced positive ratio \(a_1/a_2=s/t\) and the exact nonzero
signed difference
\[
                       F=Z_M(c_1)a_1-Z_M(c_2)a_2 .
\]
Then the number of genuine original ordered pairs with this
record is \(O(n^{1/6-\Delta}X^\varepsilon)\).
Indeed \(a_1=g s,\ a_2=g t\) for a common integer \(g\mid |F|\).
Choose \(g\) at divisor cost, apply LM1 to
\(sZ_M(c_1)-tZ_M(c_2)=F/g\), and recover both rows from
their targets and the now known \(a_i\). The result is
unchanged when \(a,b\) exchange roles.

A joint catalogue of price \(f\) bounds a pair sector by
\(n^{f+1/6-\Delta}X^\varepsilon\). A sector of full pair
exponent \(2/3-o(1)\) therefore needs \(f\ge1/2+\Delta-o(1)\).
A catalogue of price strictly below \(1/6+\Delta\) also
gives a strict row saving on a subfamily in which each row
occurs as a first component of at least one retained pair.
For \(F=0\), fix the first row and thus the repaired summand.
The second multiplier divides it, so there are at most
\(n^{1/K}X^\varepsilon\) second rows by target reconstruction
and the exact value of the second original summand. Such
pairs are power-small among all pairs of a near-full family.

There is a completion using only actual shifts. Prescribe
the exact nonzero weighted multiplier shift of LM1 and
the exact original normalized shift
\(h=(c_2/g_c)a_1-(c_1/g_c)a_2\), \(g_c=\gcd(c_1,c_2)\).
First recover both complete targets by LM1. In the original
shell \(X/2<c_i\le X\) with balanced caps on both summands,
FL1 counts the first original row at exponent
\(1/6+\tau/2+\eta\) when each summand cap is
\(n^{1/3+\tau}\). The exact normalization and \(h\) then
give \(a_2\) and complete the second row. Hence this
two-shift record has pair exponent at most
\(1/3-\Delta+\tau/2+\eta\), retaining all three original caps.
[Fixed-target original count (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

## Multiplier expansion and the repaired scale

**LM3.** Write \({\cal Z}\) for the dictionary of actual
least multipliers of any subset in the original shell with
these caps. For any \(\varepsilon>0\), its differences obey
\[
 \bigl|{\cal Z}-{\cal Z}\bigr|-1
      \ \gg_{\varepsilon}\ |{\cal Z}|(|{\cal Z}|-1)
                           n^{-1/6+\Delta}X^{-\varepsilon}.
\]
Indeed each nonzero exact difference has at most
\(n^{1/6-\Delta}X^\varepsilon\) representations by LM1.
Equivalently its additive energy is at most
\(|{\cal Z}|^2+O(|{\cal Z}|^2 n^{1/6-\Delta}X^\varepsilon)\).
On a balanced near-full family of exponent \(1/3-\tau\),
FL1 and target reconstruction give
\(\log_n|{\cal Z}|\ge1/6-3\tau/2-1/K-o(1)\), with an
arbitrarily small fixed count allowance. Thus the multiplier
dictionary has a strict additive expansion power if \(K\)
is large and the tolerances are small relative to \(\Delta\).
The same nonzero difference bound holds for the distinct
repaired \(a\)-values inside each genuine fixed-\(a\) fibre,
since those values are a fixed dilation of a multiplier
dictionary.

In a balanced shell also satisfying
\(\operatorname{rad}(c)\ge c^{\gamma_-L}\), one has
\(H_M(c)\ge c^{M\gamma_-L}\) and
\(Z_M(c)\ge 2^{-M\gamma_-L}X^{M\gamma_-L-1}\).
In terms of an upper repaired height \(H_M(c)\le Y\),
take \(X=Y^{1/(M\gamma_-L)}\) for the LM1–LM2 upper count.
Their exact-record exponent becomes
\((1/6-\Delta)/(M\gamma_-)\) before \(Y^\varepsilon\).
It is strictly below \(1/(2M)\) when \(\xi<2\Delta\).
An original incidence count of exponent below
\(\gamma_-\) in base \(n\) likewise yields a strict
least-repaired quotient count below \(1/M\).

The normalized cross products in FP1–FP3 are invariant
under repairing either original row. For the common-sum
construction the repaired common lcm \(W'\) is a multiple
of the original common lcm \(W\). Each scaled summand
acquires the same factor \(W'/W\), canceled by primitive
normalization of the products. For the direct cross
products both sides acquire the factor
\(Z_M(c_1)Z_M(c_2)\), which also cancels. The multiplier
shift is a nonprojective record and retains the new fullness
input. At a deficient prime \(p\mid c\), for every
\(k\mid M\), the valuation \(M-v_p(c)\) of the multiplier
is divisible by \(k\) exactly when \(v_p(c)\) is.
Thus least repair raises these numerical floors without
providing new exact grid divisibility. The open incidence
step is to extract a cheap catalogue of genuine mixed
records from arbitrary balanced families.
[Primitive product normalization (FP1–FP3)](ABC_FIELD_FREE_CROSS_PRODUCT_TRANSFERS.md)
