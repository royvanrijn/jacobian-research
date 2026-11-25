# Valuation profiles surviving the external arithmetic rows

Let \(M\to\infty\) and take an interval of length at
most \(H=\lceil\sqrt{M\log M}\rceil\) in \([M,2M]\)
containing \(n\gg H/\log M\) primes. Such an interval
exists by partitioning and the prime number theorem.
Write \(p_{\min}\) for the least of its primes and
\(\bar p\) for their mean. Set \(u=1/\bar p\) and
\(\lambda=L=3u\). On each of the three letters put
radical mass \(u/n\) on each of those prime degrees,
and zero mass elsewhere. Define
\[
                 E=uH/M+u/n,\qquad D=2u-2E .
\]
Every letter has numerical height moment one. At any
visible truncation \(m\), replacing degrees above \(m\)
by \(m+1\) only reduces the moment, so every tail
height constraint is satisfied. The profile is
microscopically spread over many degree positions,
while its scaled degree law converges to an atom.
On these data the inequalities below permit
\(D=2\lambda/3-O(\lambda^{3/2}\sqrt{\log(1/\lambda)})\).

For an actually nonatomic formal limit, take instead
all primes in \([M,(1+\gamma)M]\) for a fixed
\(0<\gamma<1\), and set \(H=\gamma M\) in the same
construction. The empirical scaled degree law is
nonatomic by the prime number theorem. The following
comparisons then permit
\(D\ge (2/3-O(\gamma))\lambda-o(\lambda)\).
One may choose any sufficiently small fixed \(\gamma\)
before taking \(M\to\infty\).

## Zhang's moving selectors and the inherited bounds

In Zhang's formulation each pair of letters costs
\(2u\). Any divisor at least two meets at most one
prime degree on each letter. Thus the Fourier row
(2.11) has right-hand side at least \(2u-u/(2n)\),
and the determinant and Thue rows (2.13)–(2.15)
have right-hand sides at least \(3u-2u/n\).
These comparisons allow moving indices at arbitrary
depth, beyond the indices in the checked forest.

For any released set of coordinates of total mass \(s\),
the excess height is at least \((p_{\min}-1)s\).
The geometry alternative (2.12) consequently satisfies
\[
       \max(L-s,L-1+(p_{\min}-1)s)
                     \ \ge\ L-1/p_{\min}\ \ge\ 2u-uH/M .
\]
Split the first comparison at \(s=1/p_{\min}\).
For the second use \(\bar p-p_{\min}\le H\) and
\(p_{\min}\ge M\). For every exposed pair \(i,j\)
and locked set \(K\), the first arm of the
primitive-slope estimate (3.2) is at least
\(2u-2u/n\). The cube arm of the square-full bound
(3.6) is \(L\) since all support degrees are odd.
The fifth layer is empty, so the complementary
arm (4.8) is available in every orientation.
[Zhang, §§2–5](https://doi.org/10.5281/zenodo.22151481)

Each selector in the two- and three-order divisibility
rows recovers at most \(u/n\) per letter. Thus their
right-hand sides are at least \(3u-3u/n\).
For the grouped fourth-moment selectors, the saved term
\(t_j\) is at most \(3u/n\) on each letter, so the best
fourth-moment bounds remain the inherited ones. The
permuted Hölder rows also allow \(2u-u/(2n)\).
[Divisibility rows](ABC_ALL_DEPTH_VALUATION_TRANSFER.md)
· [Grouped fourth moments](ABC_GROUPED_VALUATION_FOURTH_MOMENTS.md)

## Li and Zeng's moving same-target collision rows


Use the notation of the same-target rows in Zeng
(4.1)–(4.6). The original early-six free sets have
no mass here. The diagonal arm
\(f_A+f_B+\alpha_C\) is then \(3u\), so the original
rows allow the profile for all their modulus and
energy selectors.

The following comparisons also verify the formal
moving-index extension of these formulas. Allow any
fixed exposed pair and any fixed collections of
released and frozen indices. For a subset of mass
\(m\) and numerical weight \(h\) on one letter,
\(|m-uh|\le uH/M\). Indeed each degree differs from
\(\bar p\) by at most \(H\), and its mass is \(u/n\).
In the collision estimate,
\(K=2(s_A-p-q)+\kappa\ge2s_A-2u/n\), since
\(\kappa\ge p+q\) and each exposed mass is at most \(u/n\).
Each grouped difference selector of order at least
two removes at most one prime level. The direct
cap is \(s_B\). The one-coordinate residue cap
is at least \(s_B-u/n\), and the older two-coordinate
refinements can reduce it by at most \(2u/n\).

For the whole-value interval cap,
\[
 f_C+(h_B-w_C)_+
       \ \ge\ u w_C+(h_B-w_C)_+-uH/M
       \ \ge\ s_B-2uH/M .
\]
The first \(w_C\) bound is for the exact frozen
weighted height. A smaller certified lower bound on
that height only increases the interval cap. The
second inequality uses \(u\le1\) and the same
height-to-mass comparison. Thus for the minimum
\(\delta\) of the available caps,
\(f_C+\delta\ge s_B-2E\). The off-diagonal arm alone
has right-hand side
\[
 f_A+f_B+(u+f_C+K+\delta)/2
                 \ \ge\ 2u+f_B/2-u/n-E
                 \ \ge\ D .
\]
The named arithmetic rows, including these moving
versions, therefore permit the same profile. The
claim concerns these rows of the exponent program
and makes no arithmetic existence assertion.
[Zeng, §§3–5](https://hexagonmath.org/pdf/2611.00621v1)
· [Li, radical budget, version 4](https://rstation.ai/problems/math/papers/17/v4)
· [Li, squares and twisted energy, version 4](https://rstation.ai/problems/math/papers/28/v4)

## The diagonal floor and an actual arithmetic separation

Suppose \(N\) triples survive all frozen records.
The same-target collision step gives
\(N^2\le |W|(N+E_{\rm off})\), where \(W\) is the
remaining third-letter dictionary. The diagonal term
with \(|W|\le X^{\alpha_C-f_C+o(1)}\) costs
\(f_A+f_B+\alpha_C\) in the original exponent after the
three freezes are restored. On a balanced profile this
is at least \(L/3\). Improving the estimate only for
\(E_{\rm off}\) leaves this diagonal arm of the
same-target calculation in place.

There is an arithmetic separation on the narrowing-band
version of the profile. The narrow actual degree range
is a bounded-rank valuation grid of width \(o(1/L)\)
on the two free letters. The complementary microscopic
cap makes its high-valuation primes individually
affordable for a common full modulus. Thus the
two-letter congruence estimate yields normalized
limiting count exponent zero there for actual boxes
meeting these hypotheses. This uses information
absent from the audited external rows.
[Small-prime modulus selection and release](ABC_VALUATION_ENTROPY_GRID.md#choosing-a-full-divisor-inside-exact-valuation-layers)
