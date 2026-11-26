# Reciprocal moments and missing denominator support

A short sequence of consecutive reciprocal moments records the full
prime powers of its active integer labels, apart from a numerical
divisor controlled by their pairwise differences. For primitive
additive triples near the small-budget count exponent, compact
positive-order bands make this an exact cancellation test. The
argument also applies to the complete values in least target repairs.

## A local annihilator and a numerical divisor bound

Let \(\mathcal W\) be a set of \(m\ge1\) distinct nonzero signed
integers of absolute value at most \(Y\ge2\). Assign nonzero rational
coefficients \(\nu_W\) to the values \(W\in\mathcal W\), and let a
prime pool \(P\) contain the prime supports of all coefficient
numerators and denominators. For a fixed integer \(q\ge1\), put
\[
 t_s=\sum_{W\in\mathcal W}\nu_W W^{-s}.
\]
Let \(J\ge1\), and let \(D\) be the least common multiple of the
reduced denominators of \(t_q,\ldots,t_Q\), where \(Q=q+J-1\).
A zero moment has denominator one. At a prime \(p\) write
\[
          m_p=\#\{W\in\mathcal W:p\mid W\}.
\]
Whenever \(p\notin P\), \(p\mid W\), and \(m_p\le J\), one has
\[
 Qv_p(W)
     \le v_p(D)+\sum_{\substack{V\in\mathcal W\\V\ne W}}v_p(V-W).
                                                               \tag{RV1}
\]
In particular, \(J=m\) or any upper bound for the \(m_p\) suffices.
For the full part \(W_{\rm miss}\) supported on primes outside both
\(P\) and \(D\), assume \(m_p\le J\) at its primes. Then
\[
 W_{\rm miss}^{\,Q}
   \ \mid\ \prod_{\substack{V\in\mathcal W\\V\ne W}} |V-W|,
 \qquad
 W_{\rm miss}\le (2Y)^{(m-1)/Q}.                               \tag{RV2}
\]
There is also an estimate with \(D\) retained. If \(W_{\rm out}\)
is the full part outside \(P\) and \(m_p\le J\) at its primes,
then
\[
 W_{\rm out}^{\,Q}
   \ \mid\ D\prod_{\substack{V\in\mathcal W\\V\ne W}} |V-W| .
                                                               \tag{RV3}
\]
These statements require neither a bound on the number of distinct
prime factors of the labels nor comparability of their orders at a
prime.

For the proof at \(p\mid W\), use the integral polynomial
\[
       \prod_{\substack{V\in\mathcal W\\ V\ne W,\ p\mid V}}
                      (VT-1)=\sum_{r=0}^{m_p-1}d_r T^r .
\]
Take the last \(m_p\) moments in the list. Then
\(\sum_r d_r t_{Q-m_p+1+r}\) has valuation at least \(-v_p(D)\).
In its expansion every other \(p\)-divisible label vanishes, while
the terms with \(p\)-unit labels have integral coefficients and
unit denominators at \(p\). The remaining term is
\[
     \nu_W W^{-Q}
         \prod_{\substack{V\in\mathcal W\\V\ne W,\ p\mid V}}(V-W).
\]
Its valuation is at least \(-v_p(D)\). The factors omitted from the
product are units at \(p\), which proves (RV1). Multiplying the
prime-power divisibilities gives (RV2) and (RV3), since all
differences in their products are nonzero integers.
This uses the elementary annihilator underlying Vandermonde
interpolation. The related \(p\)-adic power-sum literature also
uses local interpolation to isolate summands.
[Laohakosol and Pitman, Applications of \(p\)-adic interpolation to exponential polynomials and sums of powers](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/23393C3B19454D4ADCD32B3AB60C164E/S1446788700017900a.pdf/applications_of_padic_interpolation_to_exponential_polynomials_and_sums_of_powers.pdf)

## Label cancellation in capped additive triples

Use \(n=X^L\), with \(L\) fixed before the height limit. Consider
at most \(k\) primitive rows \((a_i,b_i,c_i)\) with \(a_i+b_i=c_i\)
and \(X/2<c_i\le X\), where \(k\) is fixed. Choose polynomial-height
nonzero rational weights \(w_i\) and place the numerator and
denominator primes of all nonzero subsums of the weights in a
polynomial-height pool \(P\). Additional prescribed primes may
be in \(P\). For each of the three signed complete role values
\(W=a_i,b_i,-c_i\), suppose
\[
 \sum_{\substack{p\mid W,\ p\notin P\\ Lv_p(W)\ge r_0}}
                   \log_n p\ \ge\ \gamma                       \tag{RV4}
\]
for fixed \(r_0,\gamma>0\). Any one fixed compact band of this
mass suffices. The proof only needs the lower order cutoff.

Use a moment list of length \(3k\), starting at a fixed order \(q\).
At each order aggregate either the three signed reciprocal powers
per row, or use one reciprocal power for each role with separate
output lists. In the first case the output is
\[
        t_s=\sum_{i=1}^k w_i\bigl(a_i^{-s}+b_i^{-s}+(-c_i)^{-s}\bigr),
        \qquad q\le s<q+3k .
\]
In the second case use the three output lists \(t_{s,W}\) with
a separate sum over the rows for each role. Let \(D\) be the
least common multiple of all reduced output denominators in the
chosen variant. For any fixed \(0\le\mu<\gamma\), take \(q\) with
\[
                       q r_0(\gamma-\mu)>3k-1.
\]
If the outside radical of \(D\) is at most \(n^\mu\), every
signed complete value cancels after collecting its coefficients.
In the separate-list variant the cancellation holds in each role
separately. Consequently all output moments in the chosen variant
are zero.

Indeed, merge equal signed values in the pooled list, or merge
each role list separately. Every surviving coefficient is a
nonzero subset sum of the \(w_i\), with one possible repeated
summand per row if \(a_i=b_i\). For large \(X\) the latter case
is impossible in a primitive height shell, since \(a_i=b_i\)
forces \((a_i,b_i,c_i)=(1,1,2)\). Alternatively enlarge \(P\) by
the finitely many coefficient subsums with multiplicities at most
three, which also handles small shells. Any surviving label
has at most \(\mu\) normalized radical mass in the support of
\(D\) outside \(P\). Its remaining primes from (RV4) contribute
a full factor at least \(X^{r_0(\gamma-\mu)}\). With
\(m\le3k\) this contradicts (RV2) for all sufficiently large \(X\).
With the role lists and weights \(+1,-1\), signed cancellation
is precisely a trade on the complete role labels.

There is a quantitative form without assuming the output support
is small. For each surviving label \(W\), the primes in (RV4) that
are absent from \(D\) carry normalized radical mass at most
\[
                    \frac{m-1}{q r_0}+o_X(1).
\]
Thus the output denominator support must contain at least
\(\gamma-(m-1)/(q r_0)-o_X(1)\) mass from that band of each
active label. These lower bounds add over pairwise disjoint
bands, for example over the three roles of any one primitive row
when all three labels of that row survive. The test controls the
union of output denominator supports, so the output numerators
may vary without restriction.

At the global balanced endpoint the finite-band consequence of
(SH1) gives an absolute uniform lower bound on some compact band
of each role for all but a power-small class of rows. For a
fixed polynomial-height pool, imposing (RV4) with a smaller
\(\gamma\) loses only a power-small class by (EB5), which discounts
positive raw radical mass carried by that pool. For fibres with bounded caps and any fixed positive normalized count
exponent, (FT1) and moment tightness give a corresponding uniform
compact-band test on the two free roles. A sequence of violations
at that exponent would have a free-letter limit supported at zero,
contrary to its positive nonatomic mass. If the weights cancel the fixed role, it gives no active
label. If all its primes lie in \(P\), it can be retained as an
additional untested label in (RV1), while the free-role conclusion
is unchanged.
[Uniform compact bands at the unrestricted endpoint](ABC_PRIMITIVE_SECANT_RADICAL_AMPLIFICATION.md)
· [Prescribed-pool radical discount](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Fixed-fibre compactness input](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#nonatomic-mass-outside-prescribed-prime-pools)

## Least target repair and separation of moment sums

Let \(Z_M(c)\) be the least multiplier used to make the repaired
target \(M\)-full, and suppose \(ML\) is bounded by a fixed constant.
If the original raw cap is \(n^B\), the complete repaired values
have absolute size at most
\[
                   X^{1+BML}.
\]
At each original prime of any specified role the corresponding
repaired role retains at least its original order. The original
band in (RV4) therefore still supplies full powers of numerical
weight at least \(r_0(\gamma-\mu)\) in base \(X\) after removing
the output support. The preceding proof applies to the actual
repaired complete labels on taking
\[
                    q r_0(\gamma-\mu)>(3k-1)(1+BML)
\]
with a fixed margin and a fixed upper bound on the right.
In particular, a moment vector of the repaired triples with
small outside denominator radical forces exact repaired label
cancellation under the separate-role test. The rows counted
are the primitive quotients with their prescribed least repairs.
No arbitrary common multiplier enters this assertion.
[Least target repair](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

As a consequence, given a good primitive row set \(\mathcal E\)
satisfying (RV4) outside a fixed pool, record the separate-role
moment vectors under the addition law in
\(\mathbb Q/\mathbb Z[P^{-1}]\) in each coordinate. Here
\(\mathbb Z[P^{-1}]\) consists of rationals supported on \(P\) in
their reduced denominators. Fix \(h\ge1\). Use moments for at
least \(2h\) consecutive orders per role, since there are at most
\(2h\) active labels in each role of the difference. Take \(q\)
as above with \(k=2h\) and \(\mu=0\). Two ordered \(h\)-sums of the recorded
vectors agree in this quotient exactly when their three role
multisets agree. The implication from equality follows from the
separate-list cancellation result, using the two input lists with
weights \(+1,-1\). The reverse holds as an equality of rational
vectors. Each ordered left list has at most \((h!)^2\) ordered
lists sharing the same multisets, since choosing the ordered
summands determines every target. Hence the \(h\)-fold sumset in
the quotient has at least \(|\mathcal E|^h/(h!)^2\) elements.
The same conclusion holds for the actual least-repaired labels
with the adjusted choice of \(q\). Here the use of a whole moment
list gives a support test. A single low-order moment does not
separate all label lists this way.
