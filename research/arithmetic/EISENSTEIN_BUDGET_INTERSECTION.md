# Where the Eisenstein sieve and divisor budgets meet

The one-sided composite mark and the current shell sieve do
not overlap, even with an exact one-sided mark. A two-sided
correction gives a positive limiting margin at
\(\eta=0.116\), with marked prime products bounded by
\(N^{0.495}\). The roughness and product budgets fit
separately just below \(N^{1/2}\). The current joint
absolute-error estimate charges their product, which explains
why the two separate fits do not furnish a pointwise proof.

## The margin and the errors it must absorb

Work with the value-normalized, squarefree, coprime Eisenstein
pairs \(a+b=N\), \(N\equiv2\bmod6\), at \(N^\eta\) roughness.
Let \(R_N\) be the total, \(C_N\) the number with a composite
first summand, and \(Q_N\) any count of a subset with both
summands composite. Then inclusion–exclusion gives the exact
one-sided inequality
\[
                   G_N\ge R_N-2C_N+Q_N.                       \tag{1}
\]
Here \(C_N\) can be replaced by an upper bound and \(Q_N\)
by a lower bound. With a common main scale \(R_*\), put
\(c=C_*/R_*\), \(q=Q_*/R_*\), and \(m=1-2c+q\).
Errors measured in units of \(R_*\) must meet
\[
                  e_R^-+2 e_C^++e_Q^-<m.                     \tag{2}
\]
For three equal relative errors measured against the separate
main terms, the sufficient bound is \(e<m/(1+2c+q)\).
If the two mark estimates can be proved as **ratios to the
actual \(R_N\)**, the needed budget is \(2e_C^++e_Q^-<m\)
in units of \(R_N\). A separate lower bound for \(R_N\)
then sets the scale of the resulting prime count.
This version allows a common multiplicative error in all
three counts.

For comparison, the previous proper-divisor inequality
\(G_N\ge R_N-D_N\) has a budget \(e_R^-+e_D^+<m\) in
units of \(R_*\). Its approximately \(5.46\%\) margin at
\(1/5\) allowed about \(2.8\%\) relative error on each
of the two main terms. The two-sided correction below
has a smaller budget and shorter individual moduli.

## The one-sided and two-sided shape limits

Use the constants \(c_j(u)\) of the
[fixed-factor distribution](EISENSTEIN_ROUGH_FACTOR_DISTRIBUTION.md#the-factor-count-distribution)
and write
\[
                 T(u)=\sum_{j\ge1}2^{1-j}c_j(u).
\]
Only \(j<1/u\) contribute, with \(c_1=1\) for the prime
profile. In the unnormalized first-summand shape measure,
\(T(u)\) is the mass of profiles all of whose exponents are
above \(u\). This uses the split-prime density \(1/2\)
and counts distinct values. On almost all targets
\(C_N/R_N\) tends to \(1-1/T(\eta)\).
Consequently an exact one-sided composite mark has limiting
margin \(2/T(\eta)-1\). It changes sign at a unique
\(\eta\) in \((0.164,0.165)\). The sign assertions at the
endpoints are certified below. Uniqueness follows from
nested simplices and the strictly decreasing \(c_2(u)\).
Even the semiprime piece alone is enough to give a negative
margin at \(\eta\le 1/(1+e^2)\simeq0.11920\).

For the new lower mark, choose finitely many bin endpoints
\(\eta=u_0<u_1<\cdots<u_m<1/2\). For \(a>1\) write
\(x(a)=\log P^-(a)/\log a\), where \(P^-(a)\) is its least
prime factor. Count a pair in \(Q_N\) if
\[
 x(a)\in(u_{i-1},u_i],\qquad x(b)\in(u_{j-1},u_j],
 \qquad u_i+u_j\le\tau .
\]
A prime has \(x=1\), so every selected pair has two composite
summands. Each marked product satisfies the exact bound
\[
           P^-(a)P^-(b)\le a^{u_i}b^{u_j}\le N^\tau .
\]
The bins are disjoint, so the count uses each pair only once.
The unnormalized bin masses are
\(w_i=T(u_{i-1})-T(u_i)\). With
\(W=\sum_{u_i+u_j\le\tau}w_i w_j\),
the limiting \(Q_N/R_N\) is \(W/T(\eta)^2\).

This shape comparison follows by applying the
[fixed-box refinement](EISENSTEIN_DIVISOR_MODULUS.md#factor-shapes-and-the-constants)
to both summands at once. That proof uses bounded-degree
prime products with a fixed interval specification for each
labelled factor. Each choice on either summand has the same
major-arc multiplier. The minor-arc bound for the product
of the two box-restricted sums follows from the same
supremum and Parseval estimates. The endpoint conditions
on the least factor are finitely many hyperplane cuts in
the factor exponents and have measure-zero boundaries.
Inner and outer finite box approximations then give the
claimed bin limits. For any fixed tolerance and \(A>0\),
the exceptional set in \((X,2X]\) is \(O_A(X/\log^A X)\),
with the implied constant depending on the tolerance and
parameters.

**Theorem (separate budgets fit).** Take
\(\eta=29/250=0.116\), \(\tau=99/200=0.495\), and the
bin endpoints \(\eta\) and the multiples of \(1/200\)
from \(0.120\) through \(0.375\). The limiting value of
\((R_N-2C_N+Q_N)/R_N\) on the almost-all sets is in
\([0.004498,0.004734]\). In detail,
\[
 C_*/R_*\in[0.585153,0.585198],\qquad
 Q_*/R_*\in[0.174894,0.175040].
\]
The mark recovers enough two-composite pairs to exceed
\(2C_*/R_*-1\). The common relative error \(0.19\%\) in
(2) fits the certified margin. In ratio form, an upper
error of \(0.1\) percentage point in \(C_N/R_N\) and
a lower error of \(0.2\) percentage point in \(Q_N/R_N\)
fit, provided \(R_N\) is positive.
A sufficient pointwise input is therefore \(R_N>0\),
\(C_N/R_N\le0.586198\), and \(Q_N/R_N\ge0.172894\)
for all sufficiently large \(N\equiv2\bmod6\), on this
same value-normalized shell. Those bounds together
give \(G_N\ge0.000498 R_N\).

The bound \(N^{0.495}\) concerns the product of the two
marked primes. Checking that those primes are **least**
and estimating \(C_N\) are additional arithmetic inputs.
For example the composite band with least first factor
\(>a^{0.495}\) consists only of semiprimes and carries
\(0.4148\%\) to \(0.4149\%\) of \(R_*\). Leaving that
entire band uncertain would cost more than (2) allows
in the upper estimate for \(C_N\).

## The product budget

Write \(\beta=4.267\) for the sufficient lower-sieve
ratio used in the
[dimension-two argument](RANDOM_SPLIT_PRIME_SUMSETS.md#rough-decompositions-with-controlled-supports).
A direct application asks for a roughness level
\(D>N^{\beta\eta}\). On the macroscopic shell, the
[Hecke absolute joint remainder](EISENSTEIN_BOUNDED_PRIME_FACTORS.md#all-target-valuations-and-joint-cofactor-levels)
allows a summed cofactor-times-sieve budget below \(N^{1/2}\)
for the stated coordinate tests. For the two-sided mark,
the two separate exponents at \(\eta=0.116\) are
\[
                  \beta\eta=0.494972,\qquad \tau=0.495.
\]
They almost coincide, and each separately lies below \(1/2\).
Multiplying the supports in the absolute joint estimate
asks for a level beyond \(N^{0.989972}\) even before budgeting
leastness and the value normalization. The marked primes
exceed the roughness cutoff while the sieve modulus is
formed from primes below it, so the two moduli are coprime.
There is no saving from taking their least common multiple.
Orienting the split norm factors yields
\(pQ(z')+qQ(w')=N\), still at total index \(N\), since
these are separate factors of the components.

The exponents here are bookkeeping for this direct use
of the sieve and marks. They do not give the estimates
at those levels. For example lower-sieve positivity
alone does not supply the \(0.19\%\) error on \(R_N\)
used in the separate-main-term interpretation of (2).

**Proposition (best common ceiling for the separate levels).**
With an exact marginal and the full least-prime product mark,
the infimum common exponent for both separate levels to
allow a positive shape margin is strictly between \(0.487\)
and \(0.489\).

Indeed at a shared ceiling \(s\), the feasibility condition
is \(\eta\le s/\beta,\ \tau\le s\). The sign of the full-cap
margin is the sign of \(1-M(\eta,\tau)\), with \(M\)
defined below. The most favorable sign under this ceiling
is attained at \((s/\beta,s)\) by monotonicity. Using bin
width at most \(1/2000\), the inner sum uses upper endpoints
whose sum is at most \(s\). For the upper bound, a pair
in the full mark must have upper endpoints whose sum
is at most \(s+2/2000\). The two rational enclosures in the
checker give a negative full-cap upper bound at \(s=0.487\)
and a positive inner-bin lower bound at \(s=0.489\).
The full-cap mass varies continuously, since the composite
profile densities have measure-zero boundaries at the
cutoffs. Hence these signs give strict bounds on the
infimum, also when the separate levels must fit strictly
below their ceiling.

For the earlier direct proper-divisor rule the mark
has unshortened power \(1-\eta\). The two separate
requirements \(\beta\eta\) and \(1-\eta\) meet at
\(\eta=1000/5267\), with level exponent \(4267/5267\)
(about \(0.81014\)). The limiting proper-divisor margin
there is \(1.3352\%\) to \(1.3374\%\).
Using a joint product instead asks for about \(1.62028\).
This comparison excludes any saving from switching the
large divisor into a smaller cofactor.

**Proposition (no meeting under the direct joint cap).**
For the capped two-composite correction in (1), with
the exact marginal \(C_N\) granted, there is no positive
factor-shape margin under the simultaneous restrictions
\(\beta\eta+\tau\le1/2\) and \(\eta>0,\tau\ge0\).

To prove this, let \(M(\eta,\tau)\) be the unnormalized
mass of two composite profiles whose two least exponents
sum to more than \(\tau\). Even the full cap, allowing
every product that meets the limit, has normalized
margin \((1-M(\eta,\tau))/T(\eta)^2\). This follows from
\(T^2-2T(T-1)+(T-1)^2=1\). In the shell contribution
to the shape limit both \(a\) and \(b\) exceed
\(X/(\log X)^4\), so the exponents relative to \(\log N\)
differ from those relative to \(\log a\) or \(\log b\) by
\(O(\log\log X/\log X)\). The endpoints outside this shell
are negligible by the fixed-box proof. Since the sum of the
two composite minima has no atom at \(\tau\), the actual
product bound and the sum bound have the same limiting mass.
The missed mass increases as either \(\eta\) or \(\tau\)
decreases. A lower rectangle sum gives
\(M(0.08,0.5)>1.06722\), so all \(\eta\le0.08\) fail.
For \(\eta>0.08\), the cap has
\(\tau\le1/2-\beta\eta<2\eta\), so there is no two-sided
correction. The constraint \(\tau\ge0\) forces
\(\eta\le1/(2\beta)<0.164\), where the one-sided sign is
negative. The argument concerns this cap and this joint
budget. A signed distribution theorem that avoids
multiplying the two levels would be a different input.

## Does normalization lower the sieve dimension?

At a split prime \(p\nmid N\), the normalized local **value
model** gives each residue \(a\bmod p^2\) equal weight
before the squarefree restriction. The admissible classes
exclude \(0,N\bmod p^2\), and \(2(p-1)\) of the remaining
\(p^2-2\) classes have \(p\mid a(N-a)\). The local sifting
density is therefore \(2/p+O(p^{-2})\) at split primes,
which have density \(1/2\). This is a dimension-one
model. A raw coordinate shell instead has leading density
\(4/p\) at split primes and dimension two.
Even if the **normalized** remainders were available in
a dimension-one linear sieve at level \(D>N^{2\eta}\),
the square-root mark at \(\eta=0.220\) would need about
\(0.440+0.500=0.940\) in the unsplit joint budget.
The individual roughness bound \(0.440\) would fit.

This dimension change is a local calculation, rather than
a transfer of the raw Hecke estimate. On squarefree values,
expanding the normalization gives
\(2^{-\Omega(a)}=\sum_{d\mid a}(-1/2)^{\Omega(d)}\).
Truncating that expansion at \(d<N^{1/2}\) gives \(1\)
for prime values \(a>N^{1/2}\), whereas their weight is
\(1/2\). Thus the simplest truncation fails on the very
portion whose coefficient we want to detect.

Even under that favorable dimension-one model, the direct
joint cap has no positive limiting margin **with a
nonempty two-sided correction**. For a nonempty correction
one needs \(\tau>2\eta\), so \(2\eta+\tau\le1/2\) forces
\(\eta<1/8\). The range \(\eta\le0.08\) was already excluded.
In \(0.08<\eta<1/8\) the cap forces \(\tau<0.34\),
and the same rational rectangle argument gives
\(M(1/8,0.34)>1.37250\). With larger \(\eta\), the cap leaves
only the one-sided criterion and does not estimate its
marginal. Thus lowering the formal dimension by itself
does not establish a joint distribution result.

## Exact enclosures

Run `python3 research/arithmetic/verify_eisenstein_budget_intersection.py`
from the repository root. At rank \(j\ge3\) the checker
divides the simplex of (1) in the factor-count note into
a grid of width \(h=(1-j\eta)/192\).
For full cells convexity bounds each reciprocal integral
by midpoint and trapezoid values. A boundary cell with
integer slack \(b\) occupies the fraction
\[
 \frac{1}{d!}\sum_{i=0}^b(-1)^i\binom di(b-i)^d,\qquad d=j-1,
\]
of the cube, by inclusion–exclusion on the unit-cube faces.
Endpoint bounds enclose its integrand.
The cell products at each integer slack are coefficients
of a truncated polynomial power, so the enclosures remain
practical for eight or more factors. Rounding is directed
with exact integers. Rank two uses the positive atanh series.
The checker also compares the new full-factor enclosures
against the previous box bounds at \(1/5\).

For the selected-bin correction, form
\(A_i=T(\eta)-T(u_i)\) and let \(j(i)\) be the largest bin
endpoint index allowed with \(u_i\).
Abel summation rewrites the count as
\(\sum_i A_i(A_{j(i)}-A_{j(i+1)})\), with a terminal zero.
Replacing both cumulative masses by monotone lower or upper
envelopes bounds the bin sum in the corresponding direction.
For a missed-mass lower bound, use disjoint rectangles
whose lower corners already meet the strict-sum condition
away from a boundary of measure zero. These two calculations
check both the positive separate-level model and the
negative joint-budget model.
