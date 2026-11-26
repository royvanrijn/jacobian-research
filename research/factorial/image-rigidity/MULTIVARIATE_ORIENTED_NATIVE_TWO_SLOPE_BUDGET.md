# Oriented native budgets for two shared Newton slopes

Fix an original ordinary two-monomial profile with distinct
exponent columns. Omit equal rows. For each remaining row let
\(a_i\) be its smaller exponent, \(d_i\) its positive gap,
and \(\sigma_i\in\{+,-\}\) its orientation. Choose either
order \(m\) of an adjacent pair \(n,N=n+1\) with \(n\ge2\),
and write \(m'\) for the other order.

For each \(j\in\{m,m'\}\) define literal prefix products
\[
 A_{j,k}^\sigma=\prod_{i:\sigma_i=\sigma}
                       \prod_{t=1}^{k d_i}(j a_i+t),
 \qquad B_j^\sigma=A_{j,j}^\sigma .
\]
Use empty products equal to one. These quantities require
no moment polynomial to be expanded. The direction of the
profile can be reversed in the tests below.

**Theorem.** Any physical common root forces the integer
divisibility
\[
                (A_{m,1}^{+})^2 \ \mid\
                         (B_m^-)^2 (B_{m'}^+)^{m}.            \tag{1}
\]
Thus at any prime the marked-order increasing first
blocks can be paid for only by opposite-direction factors
at the marked order or same-direction factors in the
other order. The powers in (1) use the fact that a
rational common factor has even degree. No bound on
the number of opposite or same-direction rows, their
heights, gaps, or prime depths is needed.

A sharper factorization-free test uses the rational
\[
 D_{k,\ell}=
    \frac{(A_{m,k}^+ A_{m,m-k}^-/B_m^-)^\ell}
         {(B_{m'}^+/(A_{m',m'-\ell}^+ A_{m',\ell}^-))^k},
                  \quad 1\le k\le m,\quad 2\le\ell\le m'.     \tag{2}
\]
A common root forces the gcd of the positive reduced
numerators of all \(D_{k,\ell}\) to be one. One may
multiply \(D_{k,\ell}\) by
\(\binom mk^\ell\binom{m'}{\ell}^{\,k}\), obtaining a
test that includes the binomial contributions.
These are necessary conditions on the actual original
profiles, at all primes simultaneously.

## Proof of the native inequalities

At each order strip the smaller-end factorial from
every row, which changes only a scalar. Write
\(\xi_{j,k}\) for the valuation of the resulting
coefficient at index \(k\) at a fixed prime \(p\).
Then
\[
 \begin{aligned}
  \xi_{m,k}-\xi_{m,0}
    &=
       v_p(A_{m,k}^+ A_{m,m-k}^-/B_m^-)
                         +v_p\binom mk,\\
  \xi_{m',m'}-\xi_{m',m'-\ell}
    &=
       v_p(B_{m'}^+/(A_{m',m'-\ell}^+ A_{m',\ell}^-))
                         -v_p\binom{m'}{\ell}.
 \end{aligned}                                                 \tag{3}
\]
The smallest unit slope at order \(m\) is the minimum
of the secants \((\xi_{m,k}-\xi_{m,0})/k\).
The second-largest unit slope at order \(m'\) is at
most the largest suffix secant
\((\xi_{m',m'}-\xi_{m',m'-\ell})/\ell\) with \(\ell\ge2\).
Indeed start at the left endpoint of the lower
polygon edge containing the second-to-last unit
step and take the secant to the last endpoint.
This secant has length at least two and averages
slopes all at least the second-largest one.
The slope lists count each edge by its width.

The even moment is strictly positive on the real line,
so any rational common factor has degree at least two
if nonconstant. By the Newton polygon product rule,
the smallest slope of order \(m\) cannot exceed the
second-largest slope of order \(m'\) at a collision.
Consequently some \(k,\ell\) in (3) satisfy the
corresponding secant comparison. This proves (2),
with or without the binomial contributions: the
first binomial term is nonnegative and the second
nonpositive. In particular, their removal can only
weaken the necessary comparison.

If \(v_p(A_{m,1}^+)-v_p(B_m^-)>0\), the smallest
marked slope is at least this difference divided by
\(m\), since all increasing prefixes contain their
first block and every decreasing valuation drop is
at most that of the full decreasing ramp.
Otherwise the claim (1) at this prime is automatic.
Every suffix secant of length at least two in the
other order is at most \(v_p(B_{m'}^+)/2\),
since the decreasing ramps can only lower the
suffix difference. This proves (1) prime by prime.

## Exact corrections after common low slopes are removed

For a prime \(p\) and a cutoff \(b\ge0\) define
\(a^\sigma_{j,k}\) to be the sum of \((v_p(f)-b)_+\)
over the factors \(f\) in \(A_{j,k}^\sigma\),
counting multiplicity. These are precisely the high
layers of the valuation. The lower layers through
\(p^b\) have a common linear slope in the two orders,
by the factorial floors. After subtracting it let
\(E_j\) be a bound on the oscillation of the remainder
of the lower layers at order \(j\), without the
binomial term. One may use \(E_j=b\) times the number
of unequal rows, or account for cancellations by the
[cyclic discrepancy theorem](MULTIVARIATE_RETAINED_CYCLIC_DISCREPANCY.md).

Put
\[
 d^+_m(k)=a^+_{m,k},\quad
 d^-_m(k)=a^-_{m,m}-a^-_{m,m-k},\quad
 s^+_{m'}(\ell)=a^+_{m',m'}-a^+_{m',m'-\ell},\quad
 s^-_{m'}(\ell)=a^-_{m',\ell}.
\]
After common slope removal the proof of (3) gives the
following necessary condition at a collision:
\[
  \min_{1\le k\le m}\frac{d^+_m(k)-d^-_m(k)-E_m}{k}
   \ \le\
  \max_{2\le\ell\le m'}
               \frac{s^+_{m'}(\ell)-s^-_{m'}(\ell)+E_{m'}}{\ell}.
                                                               \tag{4}
\]
The binomial contributions can again be dropped in
these directions. With \(k=m\) and \(\ell=2\) used
only as bounds on the lengths, (4) implies
\[
 2a^+_{m,1}
   \ \le\ 2a^-_{m,m}+m a^+_{m',m'}+2E_m+m E_{m'}.             \tag{5}
\]
In deriving (5) it suffices to consider
\(a^+_{m,1}-a^-_{m,m}-E_m>0\), since the other
case gives the inequality directly. The common
affine subtraction is justified by the same floor
identity as in the
[first-block depth theorem](MULTIVARIATE_FIRST_BLOCK_DEPTH_AND_MOVING_RESIDUES.md).
Here (4) retains the actual orientation of both
arrays, and no designated row is required.
The Newton polygon comparison uses
[Filaseta and Rocks, Theorem 2.1](https://doi.org/10.4064/cm6761-3-2016).

## An adaptive all-degree prime mass criterion

Take an integer \(J\ge4\) at least the
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
At each pair \(2\le n<N<J\), one may choose any marked
order \(m\), orientation, and row occurring \(r\ge1\)
times, with smaller exponent \(L\) and gap \(d\).
Reverse if necessary to make the row increasing.
Let \(\mathcal P_m\) be any set of primes, and let
\[
 F_m=\prod_{1\le t\le d}(mL+t),\qquad
 X_{\mathcal P_m}=\prod_{p\in\mathcal P_m}p^{v_p(X)}
\]
for a positive integer \(X\).
From (1) a collision forces
\[
       (F_m)_{\mathcal P_m}^{\,2r}
              \mid ((B_m^-)_{\mathcal P_m})^2
                    ((B_{m'}^+)_{\mathcal P_m})^m.             \tag{6}
\]
Thus only the indicated hostile ramps at the actual
mark primes need to be charged. Failure of (6) at
one such choice for every low pair proves SFC(2)
in every degree. The first pair separates by even
positivity and the high pairs by the cutoff.
There is no requirement to use the same selected
row, orientation, mark order, or prime set across
the terminals.

For a progression consequence, suppose \(d\) grows
and take for \(\mathcal P_m\) the primes exceeding
\(nNd\) in the product of those first-block factors
with \(t\equiv1\bmod m\). Write \(Q\) for the number
of terms of this progression, so \(Q\ge d/m\).
At a possible collision the
[absolute curvature bound](MULTIVARIATE_ABSOLUTE_CURVATURE_12.md)
gives \(L+d>d^2/12\), hence \(L>d^2/24\) for
large \(d\). For each prime at most \(nNd\), delete
one term with its largest prime-power part. The
remaining contributions of this prime are bounded
by its valuation in \((Q-1)!\): at every power its
divisible indices are in one residue class since
the progression step \(m\) is a unit at the prime.
Primes dividing \(m\) do not divide any term.
Writing \(P\) for the product of the selected terms,
\[
 \log P_{\mathcal P_m}\ \ge\
       Q\log(mL)-\pi(nNd)\log(mL+d)-\log((Q-1)!).
\]
The elementary bound \(\pi(t)=O(t/\log t)\) gives
\(\pi(nNd)\le Q/8\) for sufficiently large \(d\) in
terms of \(J\). Also \(\log(mL+d)\le(9/8)\log(mL)\)
and \(\log Q\le(3/5)\log(mL)\) for large \(d\).
The right side is then at least
\(Q\log(mL)/4\). These thresholds are computable
using the elementary estimates in the
[short-window progression argument](MULTIVARIATE_INDEPENDENT_HEIGHT_WINDOW_CLOSURE.md#unbounded-windows-and-unbounded-gaps).
Consequently a necessary condition at such an
original profile is
\[
          2\log (B_m^-)_{\mathcal P_m}
             +m\log (B_{m'}^+)_{\mathcal P_m}
                         \ \ge\ \frac{rd}{2J}\log L.           \tag{7}
\]
This yields an effective all-degree family with
arbitrary retained heights and gaps whenever the
directed hostile budgets are smaller at every low
pair. The same-direction marked ramps and opposite
unmarked ramps have no cost in (7).

## Strict ordinary cases

At orders \(3,4\), use two increasing unit rows with
smaller exponents \(L=(17^7-1)/4\) and \(A=17^3-1\).
At \(17\) their stripped valuation arrays are
\((0,0,0,3)\) and \((0,7,7,7,10)\), respectively.
The two polygons share slope \(3\), with width one.
In (1), the marked first-block depth is \(7\),
the hostile same-direction other-order depth is
\(3\), and the opposite depth is zero.
Thus (1) would require \(14\le12\).
The maximal other slope on its own would not give
that contradiction.

For an illustration of benign large depths at the
same pair, use one increasing unit row at smaller
exponent \(4\), one increasing unit row at
\(A=(17^e-1)/2\), and one decreasing unit row at
\(B=(17^e-1)/3\), with \(e\ge2\) even.
At \(17\) the first row supplies the upper first
mark of depth one. The second row supplies upper
depth \(e\) at the second block and no lower hit.
The third supplies lower first-block depth \(e\)
in the decreasing orientation and no upper hit.
Both hostile budgets at this prime vanish.
Equation (1) excludes a common factor for all \(e\),
although the total retained prime depth diverges.
