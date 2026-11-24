# The residual ordinary binomial Strong Factorial problem

The [original polynomial triangle](POLYNOMIAL_TRIANGLE_AND_TERMINAL_SIX.md)
defines the actual moment polynomials in (T1) and the gap cutoff in (T2).
Work throughout with ordinary factorial moments and their rational polynomials.

## Residual region

Write \(N_*\) for the unevaluated common cutoff in
[SupportUniformEventualOrdinaryBinomialSFC](SUPPORT_UNIFORM_EVENTUAL_ORDINARY_BINOMIAL_SYNTHESIS.md#3-a-finite-covering-of-every-positive-lower-exponent-support).
For each \(N\), keep \(D_N\) from the triangle note's (T2). The remaining ordinary binomial
problem is equivalent to rational coprimality at the triples of integers
\((a,d,N)\) satisfying all the following conditions.

\[
\begin{gathered}
13\le N<N_*,\qquad 2\le d<D_N,\qquad 4\le a,\qquad d<3a,\\
a<\min\{d(N-2),\,10^{24}d^2,\,10^{20}d^4\},\\
a=4\ \Longrightarrow\ d\in\{6,9,10\}.
\end{gathered}
\tag{T5}
\]

The two polynomials are the actual pair \(P_{N-1,a,d},P_{N,a,d}\)
in the triangle note's (T1). Equivalence holds for the ordinary conjecture over \(\mathbb C\)
and for the conjecture over all characteristic-zero fields.

The [terminal-through-twelve theorem](STEPANOV_AUXILIARY_DETERMINANTS_SFC2.md)
covers \(2\le N\le12\). The common cutoff covers \(N\ge N_*\), and the triangle note's Section 1 gives the bound \(d<D_N\). The
[sharper shape separation](GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md#3-shape-separation-and-a-sharper-soft-complement)
gives coprimality already at \(a=d(N-2)\) and above. The independent
[quadratic](UNIFORM_QUADRATIC_STRONG_FACTORIAL_TAIL.md) and
[quartic](UNIFORM_QUARTIC_STRONG_FACTORIAL_TAIL.md) all-terminal tails
supply the other two upper bounds on the candidate lower exponent.
The [constant-plus-power](ALL_UNIVARIATE_TWO_TERM_SFC.tex),
[adjacent](UNIFORM_PADE_PARABOLA_SFC.tex),
[lower-one](ORDINARY_LOWER_ONE_STRONG_FACTORIAL.md),
[lower-two](ORDINARY_LOWER_TWO_STRONG_FACTORIAL.md) and
[lower-three](STEPANOV_SFC_PROGRESS_AND_LOWER_THREE_CLOSURE.md) results
give \(d\ge2\) and \(a\ge4\). The
[relative-gap theorem](RELATIVE_GAP_THREE_STRONG_FACTORIAL.md) covers
\(d\ge3a\) for \(a\ge2\).

For \(a=4\) the remaining gaps before the final restrictions are
\(2,\ldots,11\). All primes are covered by the
[lower-four prime-gap theorem](STEPANOV_LOWER_FOUR_ELEVENTH_GAP_CLOSURE.md#7-all-prime-gaps-at-lower-exponent-four).
The [dyadic theorem](DYADIC_BINOMIAL_SFC.tex) covers \(d=4\);
the [gap-eight theorem](STEPANOV_LOWER_FOUR_EIGHTH_GAP_CLOSURE.md) covers
\(d=8\). Thus lower-four SFC alone is equivalent to coprimality at
\(d\in\{6,9,10\}\) and \(13\le N<N_*\). All three gaps satisfy the
other inequalities of (T5) for \(a=4\) in that order range.

Each restriction excludes only a domain already proved. Conversely,
coprimality on (T5) and on its complement gives the full result by the triangle note's (T1)
and rational Bézout. A failure on (T5) gives a common complex root and
hence a genuine complex counterexample. The ranges are finite since
\(N_*\) is finite. Evaluating a sufficient \(N_*\) and deciding the
remaining coprimalities are open obligations.
