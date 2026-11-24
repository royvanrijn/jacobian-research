# The residual ordinary binomial Strong Factorial problem

The [original polynomial triangle](POLYNOMIAL_TRIANGLE_AND_TERMINAL_SIX.md)
defines the actual moment polynomials in (T1) and the gap cutoff in (T2).
Work throughout with ordinary factorial moments and their rational polynomials.

## Residual region

Write \(N_*\) for the unevaluated common cutoff in
[SupportUniformEventualOrdinaryBinomialSFC](SUPPORT_UNIFORM_EVENTUAL_ORDINARY_BINOMIAL_SYNTHESIS.md#3-a-finite-covering-of-every-positive-lower-exponent-support).
For each \(N\), keep \(D_N\) from (T2). The remaining ordinary binomial
problem is equivalent to rational coprimality at all integer triples
\((a,d,N)\) satisfying
\[
\begin{gathered}
25\le N<N_*,\qquad 2\le d<\min\{D_N,128(N-1)\},\qquad
5\le a,\qquad d<3a,\\
d^2<128(a+d),\qquad
a<\min\{d(N-2),\,10^{24}d^2,\,10^{20}d^4\},\\
(a,d)\notin\{(p^w,p^s):p\text{ prime},\ 0\le w\le s,\ s\ge1\},\\
(a,d)\notin\{(2p^w,p^s):p\text{ prime},\ 1\le w<s\}.
\end{gathered}
\tag{T5}
\]
At prime gap \(d=p\), impose the two further necessary conditions
\[
 \gamma=a(N-1)+1,\quad C=\gamma+pN,\qquad
 v_p\!\left(\frac{pN}{(C)_a}\right)\le-2,\qquad
 [Y^2](\gamma+Y)_a\equiv0\pmod p.
 \tag{T6}
\]
The two polynomials are the actual pair \(P_{N-1,a,d},P_{N,a,d}\)
from (T1). Equivalence holds for the conjecture over \(\mathbb C\) and
for the conjecture over all characteristic-zero fields.

The [terminal-through-twenty-four theorem](ORDINARY_BINOMIAL_TERMINALS_THROUGH_TWENTY_FOUR.md)
covers \(2\le N\le24\). The common eventual cutoff covers \(N\ge N_*\).
The triangle note's Section 1 gives \(d<D_N\), and
[sharp shape separation](GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md#3-shape-separation-and-a-sharper-soft-complement)
excludes \(a\ge d(N-2)\). The
[quadratic](UNIFORM_QUADRATIC_STRONG_FACTORIAL_TAIL.md) and
[quartic](UNIFORM_QUARTIC_STRONG_FACTORIAL_TAIL.md) all-terminal tails
give the other two upper bounds on \(a\).
The [effective small-ratio theorem](EFFECTIVE_SMALL_RATIO_ORDINARY_BINOMIAL_SFC.md)
excludes \(d^2\ge128(a+d)\) at every order. Its remaining
inequality gives \(a>d^2/128-d\). Together with \(a<d(N-2)\)
it gives \(d<128(N-1)\). Thus for fixed \(N\) there are
\(O(N^3)\) candidate support pairs, using at most \(d(N-2)\)
lower exponents per gap.

The [constant-plus-power](ALL_UNIVARIATE_TWO_TERM_SFC.tex),
[adjacent](UNIFORM_PADE_PARABOLA_SFC.tex),
[lower-one](ORDINARY_LOWER_ONE_STRONG_FACTORIAL.md),
[lower-two](ORDINARY_LOWER_TWO_STRONG_FACTORIAL.md),
[lower-three](STEPANOV_SFC_PROGRESS_AND_LOWER_THREE_CLOSURE.md) and
[lower-four](#lower-exponent-four) results give \(d\ge2\) and \(a\ge5\).
The [relative-gap theorem](RELATIVE_GAP_THREE_STRONG_FACTORIAL.md)
covers \(d\ge3a\) for \(a\ge2\).
The [same-prime support theorem](SFC_PRIME_POWER_DIVIDED_DIFFERENCE_CLOSURE.md#4-every-lower-power-of-the-same-prime-is-closed)
covers \(a=p^w,d=p^s\) for \(0\le w\le s\).
The [upper-response argument at offset \(2p^w\)](PRIME_POWER_HASSE_COPRIMALITY.md#offset-twice-a-prime-power)
covers \(a=2p^w,d=p^s\) for \(1\le w<s\).
The [arbitrary-offset native-prime criterion](NATIVE_PRIME_ARBITRARY_OFFSET_COPRIMALITY.md)
gives (T6), including the exclusion of \(a=2p,d=p\) at every order.

Each restriction excludes a proved domain. Coprimality on (T5)–(T6)
and on the complement gives the full result by (T1) and rational Bézout.
Conversely, a nonconstant rational gcd on this region gives a common
complex root, hence a complex counterexample. The ranges are finite
because \(N_*\) is finite. A sufficient numerical \(N_*\) and the
coprimality of the whole remaining region are open.

## Lower exponent four

Every nonzero ordinary binomial \(X^4(A+BX^d)\), for every \(d\ge1\),
has a nonzero ordinary factorial moment in each consecutive positive
pair over every characteristic-zero field.

Adjacency covers \(d=1\) and the relative-gap theorem covers \(d\ge12\).
The [all-prime-gap theorem](STEPANOV_LOWER_FOUR_ELEVENTH_GAP_CLOSURE.md#7-all-prime-gaps-at-lower-exponent-four)
covers all prime gaps. The [dyadic theorem](DYADIC_BINOMIAL_SFC.tex)
covers \(d=4\), and the [gap-eight theorem](STEPANOV_LOWER_FOUR_EIGHTH_GAP_CLOSURE.md)
covers \(d=8\). These complements leave \(6,9,10\).
At all three gaps the [terminal-through-twelve theorem](STEPANOV_AUXILIARY_DETERMINANTS_SFC2.md)
and the eventual theorem already covered the terminals outside \(13\le N<N_*\).
Within that prefix the [prime-power Hasse criterion](PRIME_POWER_HASSE_COPRIMALITY.md#lower-four-gap-nine-and-the-offset-twice-the-prime)
closes gap nine. The [gap-six proof](ORDINARY_LOWER_FOUR_GAP_SIX.md)
and [gap-ten proof](LOWER_FOUR_GAP_TEN_STRONG_FACTORIAL.md) now cover
the two other gaps at every order. Hence the whole lower-four
statement follows without an evaluation of \(N_*\).
