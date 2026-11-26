# Consecutive powers of the scalar projective coordinate

Put \(F(z)=z-z^{-1}+(1-z)^{-1}\) on \(0<z<1\). Two adjacent positive
powers \(F^q,F^{q+1}\) give a power saving for sums whose moving
denominators have small radical. The proof uses one compact
positive-valuation band of the original letters. The gain is linear
in the number of distinct active rows at every fixed order, with
one \(q\) independent of that order.

## Statement

Write \(n=X^L\), with \(X\to\infty\) first at each fixed small \(L>0\).
Fix \(B>0\). There are two settings.

- For arbitrary primitive positive \(a_i+b_i=c_i\) with
  \(X/2<c_i\le X\) and \(\operatorname{rad}(a_ib_ic_i)\le n^B\),
  set \(\kappa=3\).
- For a fixed complete-letter fibre with two moving-role normalized
  outside radical caps of sum at most \(B\), set \(\kappa=4\).
  More generally, take a fixed complete letter for each slot,
  with its support in one prescribed prime pool \(P_*\) of
  polynomial height and fixed height exponent. Use a shell
  \(X_i/2<c_i\le X_i\) with \(X^\theta\le X_i\le X\) for each
  slot, where \(\theta>0\) is fixed. Cap each original letter radical by \(n^C\) for fixed \(C\). The additional moving-role caps are on primes
  outside \(P_*\).

In the first setting \(P_*=\varnothing\). Here and below
a pool of polynomial height means a subset of primes dividing
a nonzero integer of polynomial height. Fixed settings and cap constants
are allowed in constants of the estimates. Let \(k\ge2\) be fixed, set
\(f_i=F(a_i/c_i)\), and give the slots fixed nonzero rational
weights \(w_i\) of height polynomial in \(X\). The rows are
distinct, in their original orientation. Let \(P\supseteq P_*\) be
a prescribed polynomial-height pool containing the numerator and
denominator primes of the weights. Impose
\[
 \operatorname{rad}_{P^c}\!\left(\operatorname{den}t_q\,
                             \operatorname{den}t_{q+1}\right)
       \le n^\mu,\qquad
 t_s=\sum_{i=1}^k w_i f_i^s,
 \tag{SPA1}
\]
where \(\operatorname{den}\) denotes reduced positive denominator and
\(\operatorname{rad}_{P^c}\) retains only primes outside \(P\).
The outputs \(t_s\) need not be fixed.

There exist a positive integer \(q\) and constants \(\delta,\mu_0>0\),
all independent of \(k\), such that for \(0\le\mu\le\mu_0\)
the count in (SPA1) is
\[
                 \ll X^{\{k(B/(2\kappa)-\delta)+\xi\}L+\varepsilon}.
 \tag{SPA2}
\]
Here \(\xi,\varepsilon>0\) are arbitrary, and the threshold for small
\(L\) may depend on \(k,\xi\) and the fixed height exponents of the
weights and of the additional pool \(P\). The constants \(q,\delta,\mu_0\)
depend on the cap and height parameters in the settings and, in the
fibre case, the fixed height exponent of \(P_*\). No claim is made
for an order \(k\) increasing in the height limit.

Two known row estimates are used. First, after fixing previous rows,
a row with new outside radical mass at most \(v\) costs at most
\(n^{v/\kappa+\xi}X^\varepsilon\). This is the effective radical count
in the first setting. In the second, with respective new masses
\(v_U,v_V\), the complete-letter bound costs at most
\(n^{\min(v_U,v_V)/2+\xi}X^\varepsilon\), and
\(\min(v_U,v_V)/2\le(v_U+v_V)/4\). The conversion uses budget
\(L/\theta\) in the shell \(X_i\) and the caps divided by
\(\log_n(X_i^{L/\theta})\), a quantity between \(1\) and \(1/\theta\).
All prior-row primes are included in the pricing pool.
[Effective radical counting with fixed supports](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Complete-letter bound](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

Second, for uniform constants \(0<r_0<R_0\), \(\gamma,\sigma>0\),
at most \(n^{B/\kappa-\sigma}X^\varepsilon\) rows in each slot
fail to have raw mass at least \(\gamma\) on a compact band
\(r_0\le Lv_p\le R_0\) of one tracked role outside \(P_*\).
For the global setting this is the compact-band assertion at the
balanced endpoint, applied with the raw cap \(n^B\). The usual
conversion from \(X^{BL}\) to a relative radical budget in
\(X/2<c\le X\) uses an arbitrarily small fixed budget padding.
For the complete-letter fibre it follows from positive nonatomic mass
on the free letters outside a prescribed pool. Indeed a failure for
all compact bands and every positive exponent loss gives a sequence
of fibres at exponent approaching \(B/4\) in base \(n\), with
vanishing mass on each compact positive-order interval of the
offending role. In the converted shell budget \(L/\theta\) its count
exponent is bounded away from zero. Tightness of the bounded-mass
degree measures follows from bounded first moment, so the
normalized limit of that role is supported at zero, contradicting
positive nonatomic mass. One makes the exact-layer refinements at
each fixed budget before selecting the diagonal. We may decrease
\(\gamma\) to at most \(B/4\).
[Uniform positive bands at the global endpoint](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#positive-scales-in-every-endpoint-letter)
· [Positive nonatomic mass outside fixed pools in a fibre](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#nonatomic-mass-outside-prescribed-prime-pools)

## Precision from two scalar powers

The reduced denominator of \(f_i\) is \(a_ib_ic_i\). Indeed
its unreduced numerator on that denominator is
\(a_i^2b_i+(a_i-b_i)c_i^2\), a unit at every prime of the
denominator by primitivity. Consequently its height is
\(O(X^3)\), and it has a pole at precisely the original
tracked primes outside \(P\). It is injective in the ratio,
since
\[
              F'(z)=\frac{(z^2-z+1)^2}{z^2(1-z)^2}>0.
\]
Thus it determines a primitive ordered row.

Let \(E\) be the support of the output denominators outside \(P\).
A prime outside \(P\) occurring on exactly one row survives in
\(E\) by the unique negative valuation test. At a prime \(p\)
outside \(P\cup E\) occurring on exactly two rows \(i,j\),
their pole orders \(e=-v_p(f_i)=-v_p(f_j)\) agree, again
by that test for the \(q\)-th powers. Subtract \(f_j\) times
the local \(q\)-sum from the local \((q+1)\)-sum to obtain
\[
                       v_p(f_i-f_j)\ \ge\ (q-1)e .
 \tag{SPA3}
\]
In fact the local sums \(w_i f_i^s+w_j f_j^s\) are
\(p\)-integral, and \(f_j\) has valuation \(-e\); the
difference is \(w_i f_i^q(f_i-f_j)\). After clearing the
two reduced denominators, (SPA3) is divisibility by
\(p^{(q+1)e}\). For a distinct pair the cross numerator
is a nonzero integer of size \(O(X^6)\). Thus the
numerical packet shared at such exactly-two poles by
any one pair has weight at most
\[
                 H=\frac7{q+1}
 \tag{SPA4}
\]
in base \(X\), for all sufficiently large \(X\). This is
also an upper bound on any subpacket from a compact band.

Now expose rows in an order. For a current row \(i\), let
\(S_s=-\sum_{\ell\ {\rm earlier}} w_\ell f_\ell^s\).
At a prime \(p\) outside \(P\cup E\) whose last occurrence
is on the current row, the later values are integral and
the output is integral at \(p\). Hence
\[
       S_s\equiv w_i f_i^s\pmod{\mathbb Z_p},\qquad s=q,q+1.
\]
In particular \(v_p(S_q)=-qe\), and, on subtracting
\(f_i S_q\) from \(S_{q+1}\),
\[
                    f_i\equiv S_{q+1}/S_q
                                   \pmod{p^{(q-1)e}} .
 \tag{SPA5}
\]
The prime and its pole order \(e\) come from the known denominator
of \(S_q\), of polynomial height at every fixed \(k,q\).
If \(S_q=0\) there are no such primes. Choose a subset of these
primes as a record, at height-subpower cost. If its total
full numerical weight exceeds \(7/(q+1)\), two current
candidates agreeing with the record must have zero cross
numerator by (SPA5). Thus there is at most one current
row per record. This comparison uses the forced common
orders \(e\) in both candidate denominators, yielding
the exponent \(q+1\) after clearing the denominators.

## Pricing, including rows without a compact band

In the chosen setting count the raw mass \(b_i\) on row \(i\) outside
\(P_*\), so \(b_i\le B\). Let \(m_p\) count its incidences over the
rows. Write \(D=kB-\sum b_i\), and let \(S\) be the total incidence
of the additional pool \(P\setminus P_*\) plus
\(\sum_{p\notin P,\ m_p\ge3}(m_p-2)\log_n p\).
For any order put \(v_i\) equal to the raw mass first seen at
row \(i\) outside \(P\). The singleton test gives
\[
          \sum_i v_i \le (kB-D-S+\mu)/2 .
 \tag{SPA6}
\]
These identities follow by counting prime incidences at each
multiplicity. We will reserve the term eligible for primes outside
\(P\cup E\) of multiplicity exactly two. If \(m_i\) is the
raw mass on row \(i\) which is not eligible, then
\[
                       \sum_i m_i\le 3S+2\mu .
 \tag{SPA7}
\]
Indeed \(m\le3(m-2)\) for \(m\ge3\), all pool incidences
are in \(S\), and the primes of \(E\) of multiplicity at
most two contribute at most \(2\mu\).

Choose
\[
 d_0<\min(B,\gamma,\kappa\sigma)/100,\qquad
 p_0=\min(1/4,\kappa\sigma/(8B)).
\]
If \(D+S\ge d_0 k\), (SPA6) and the ordinary pricing
save \(d_0 k/(2\kappa)\), apart from the target allowance
\(\mu/(2\kappa)\). Otherwise suppose at least \(k/4\)
rows are bad, meaning they lack the band specified above.
Thin the bad indices independently with probability \(p_0\),
and put all selected indices \(I\) first.
Their total new mass is their total raw mass less the additional
pool incidence and the repetition losses inside \(I\). The expected
pool mass and losses from primes of multiplicity at least three
are at most \(3p_0 S\). The loss from each outside double
prime can occur only if both locations are selected, so
the total expected loss there is at most \(p_0^2 kB/2\).
The selected deficit has expectation at most \(p_0 D\).
Replacing the costs of these initial rows by the compact-band
exceptional bound gives expected saving, relative to
\(\sum_i v_i/\kappa\), at least
\[
       p_0\sigma k/4 - p_0(D+3S)/\kappa-p_0^2 kB/(2\kappa)
                   \ \ge\ p_0\sigma k/8 .
\]
Some selection realizes this saving; summing over selections
and orders costs a bounded factor at fixed \(k\).

## Many recovered rows in one exposure order

In the remaining case \(D+S<d_0 k\) and fewer than \(k/4\)
rows are bad. If \(\mu\le\gamma/100\), (SPA7) and the bound
on \(D+S\) show that at least \(k/2\) rows are good and have
\(b_i\ge3B/4,\ m_i\le\gamma/2\). Indeed the two latter
inequalities fail at most \(4D/B\) and \(2(3S+2\mu)/\gamma\)
times, respectively. On each such row the eligible raw
mass is at least \(B/2\). Inside its recorded compact band
the eligible full numerical weight \(w_i^*\) is at least
\(w_0=r_0\gamma/2>0\).

Choose each exposure time independently and uniformly
on \([0,1]\). For one of these rows \(i\), condition on its
time \(u\in[1/4,3/4]\). For each other row \(j\), let \(x_{ij}\)
be the numerical weight of their eligible exactly-two
packet inside the selected compact band on \(i\). Let
\(Y_i=\sum_{j\ {\rm earlier}}x_{ij}\), and let \(Z_i\) be the raw
mass of all eligible primes on \(i\) whose partner is later.
The earlier-part indicators are independent after conditioning
on \(u\). By (SPA4),
\[
 \mathbb E(Y_i\mid u)=u w_i^*\ge w_i^*/4,\qquad
 \operatorname{Var}(Y_i\mid u)\le H w_i^* .
\]
Chebyshev therefore bounds \(\mathbb P(Y_i<w_i^*/8\mid u)\)
by \(64H/w_0\). On the other hand \(0\le Z_i\le B\) and
\(\mathbb E(Z_i\mid u)\ge B/8\), so
\(\mathbb P(Z_i\ge B/16\mid u)\ge1/15\).
Choose \(q\) large enough in terms of \(w_0\) to have
\(H\le w_0/4096\). Both events
\(Y_i\ge w_i^*/8,\ Z_i\ge B/16\) then occur jointly with
conditional probability at least \(1/32\), and with
unconditional probability at least \(1/64\) when the time
also lies in \([1/4,3/4]\). Some order realizes both
conditions on at least \(k/128\) rows.

On all these rows the prime powers counted by \(Y_i\) are
last occurrences. Their numerical weight exceeds
\(7/(q+1)\), so (SPA5) gives unique recovery per subpower
record. The prime powers counted by \(Z_i\) are first
occurrences and belong to \(v_i\). Thus each recovered
row removes a charge of at least \(B/(16\kappa)\).
There is a saving of at least \(kB/(2048\kappa)\) on this
class, apart from the target allowance.

In each class count in the selected order after binning
the raw prices on a fixed finite mesh. At recovered slots
use the bound from (SPA5) for every prefix, recording
only primes of the known denominator of \(S_q\). At chosen
bad slots use the compact-band exceptional estimate.
Dropping subsequent restrictions only enlarges the count.
Sum over the finitely many possible orders, slot choices
and price classes. The three classes give a fixed linear
saving before the target allowance and rounding losses.
Choosing \(\mu_0\) sufficiently small absorbs the target
allowance. This proves (SPA2).

## Prescribed powerful target with arbitrary gcd

Let \(A+B=N\) with \(A\) \(p\)-full, \(B\) \(q_1\)-full,
\(p,q_1\ge2\), and with no coprimality requirement. Put
\(L=1/p+1/q_1,\ n=N^L\) and assume
\(\operatorname{rad}(N)\le n^C\). For any fixed \(k\ge2\),
count distinct representations satisfying (SPA1) for
\(f_i=F(A_i/N)\), with the additional requirement that \(P\)
contains \(N\)'s prime support. At sufficiently small fixed
\(L\) the conclusion is
\[
                        \ll N^{\{k(1/8-\delta)+\xi\}L+\varepsilon}
 \tag{SPA8}
\]
with constants \(q,\mu_0,\delta\) depending on \(C\) and
independent of the fixed order. No ordering of \(A,B\) is
required in the relation.

Indeed, fix all \(g_i=\gcd(A_i,B_i)\) at divisor cost since
\(g_i\mid N\). On the primitive quotient the fixed letter is
\(W_i=N/g_i\). Its outside free caps are at most
\(N^{1/p}\) and \(N^{1/q_1}\), with combined normalized
budget one in base \(n\); the prime powers outside \(N\)
retain their original fullness. Expose the short quotients
\(W_i<N^{1/8}\) first. After fixing the supported parts
on \(N\), the smaller free dictionary costs at most
\(N^{\min(1/p,1/q_1)/8+\varepsilon}\le n^{1/16}N^\varepsilon\)
per such slot. The other letter follows from the fixed target.
Put the denominator supports of these exposed rows in the
additional pool. Subtracting their contributions preserves
(SPA1) outside the enlarged pool.

On the tall quotient slots apply (SPA2) with \(\theta=1/8\),
\(\kappa=4\), \(B=1\), \(P_*=\operatorname{supp}(N)\).
For one residual row the singleton test gives cost at most
\(n^{\mu/4+\xi}N^\varepsilon\), and no residual row requires
no cost. Since short slots cost at most \(n^{1/16}N^\varepsilon\),
decreasing the uniform gap if necessary gives (SPA8) in
every case. The coordinate \(F(A/N)\) determines the actual
fixed-target representation. The short quotient split is the
same one used in scalar exact pair recovery.
[Scalar exact relation bound on powerful targets](ABC_PROJECTIVE_RECIPROCAL_GENUS_TWO_RECOVERY.md#scalar-relations-on-powerful-representations-of-a-fixed-target)

On a least-repaired family the two scalar powers in (SPA1)
can be evaluated on the ratio of the actual repaired letters:
the ratio is unchanged by the multiplier. This retains both
outputs and the same weights and does not count arbitrary
additional multipliers. The saving controls tuples producing
two sums with small outside denominator radical. An application
to the unrestricted count would require a lower bound for the
abundance of such tuples, which is still unknown.
[Least target repair and primitive quotients](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

## Relation to existing moment methods

The two-order calculation is a short local interpolation identity.
Many consecutive moments allow elimination of all but one of the
\(p\)-divisible labels by a polynomial annihilator. The two-order
argument instead uses exactly-two poles and last occurrences,
then orders the rows so that many recoveries replace positive
raw prices. Bloom and Kuperberg count odd linear relations of
bounded numerators and denominators. Their Theorem 3 does not
impose the second power sum or track the radical of a moving
denominator. The closely related consecutive signed-reciprocal
argument also applies at every fixed order. Using powers of the
projective coordinate here reduces each recovery to one rational
value, making a single positive-order band and a linear-in-order
saving available.
[Laohakosol and Pitman, Applications of p-adic interpolation to exponential polynomials and sums of powers](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/23393C3B19454D4ADCD32B3AB60C164E/S1446788700017900a.pdf/applications_of_padic_interpolation_to_exponential_polynomials_and_sums_of_powers.pdf)
· [Bloom and Kuperberg, Theorem 3](https://arxiv.org/abs/2312.09021)
· [Exact local annihilator](ABC_RECIPROCAL_VANDERMONDE_SUPPORT.md#a-local-annihilator-and-a-numerical-divisor-bound)
· [Consecutive signed reciprocal powers](ABC_RECIPROCAL_POWERS_VARIABLE_TARGET.md)
