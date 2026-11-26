# Terminal-fourteen shared profiles at low outer means

Let \(u\) be the normalized order-thirteen reflected outer mean of a
fully shared ordinary exponent profile. The order-thirteen and
order-fourteen factorial moment polynomials have disjoint root
arguments under independent positive rescalings on the following
bands.

| \(u\) | Reciprocal cap |
| --- | --- |
| \([0.04,0.10]\) | \(1/8\) |
| \([0.10,0.11]\) | \(1/9\) |
| \([0.11,0.12]\) | \(1/8\) |
| \([0.12,0.15]\) | \(1/10\) |
| \([0.15,0.16]\) | \(1/11\) |
| \([0.16,0.21]\) | \(3/40\) |
| \([0.21,0.25]\) | \(1/14\) |

The cap in the table is on every \(d/(13\ell+1)\) in the profile.
A rational Bézout identity gives coprimality of the original
polynomials over every characteristic-zero field.

## Exact row envelopes

Use the [matched grouped logarithms](MULTIVARIATE_TERMINAL_THIRTEEN_SHARED_LIGHT.md#matched-logarithmic-enclosure)
with degrees thirteen and fourteen and \(R=169/168\). The checker
fixes these constants before the order-dependent caches are populated.
The upper skew and outer-loss ratio bounds use the
[matched excess-loss estimate](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md)
and the [linked lower estimate](MULTIVARIATE_MATCHED_OUTER_RATIO_MONOTONICITY.md#coupling-to-the-lower-order-tanh-moments).
Column reversal allows one nonnegative signed coordinate to cover
both orientations.

Write \(c\) for the declared reciprocal cap and \(u_-\) for the
lower mean endpoint of a slab or rectangle. Use a reciprocal
floor \(f\le c\). Every row with largest reciprocal above \(f\)
has gap \(d\) bounded by
\[
               (1+12 f/(1+f))^d u_-^2\le1.
\]
In fact each reciprocal on the row is at least its largest
reciprocal divided by one plus that reciprocal. The row mean
is at least the profile mean, since all the other row means
are at most one. For each possible \(d\), the requirements
\(f<d/(13\ell+1)\le c\) give a finite integer base interval.
The squared row mean
\[
         \prod_{h=1}^d\frac{13\ell+h}{13\ell+12d+h}
\]
increases with \(\ell\), so binary search finds the smallest
allowed base with squared row mean at least \(u_-^2\).
For fixed \(d\) the reciprocal and the ratio of the row outer
losses both decrease with this base. For the loss ratio this is
the [real-base monotonicity theorem](MULTIVARIATE_FACTORIAL_GAUSSIAN_ANNULUS_HOMOTOPY.md#the-row-ratio-decreases-under-inflation).
It therefore suffices to bound that one base per gap.
For rows at or below the floor, the termwise bound with cap
\(f\) applies. The profile loss ratio is a positive weighted
mean of its row ratios.

The replay uses \(f=\min(c,\max(1/80,c/2))\) for this search.
For the four preliminary slabs starting at .04, .08, .13
and .18 it instead replays the original enumeration of every
admissible base with floor \(1/80\). Both enumerations use the
same row-product condition and outward rational log bounds.

## Angular bounds from grouped root disks

The model gives a sum of grouped interval polynomials at each
order. It proposes centers from a midpoint polynomial and shifts
the groups separately in outward complex ball arithmetic. A
Rouché inequality against the linear coefficient certifies one
root in each disk. The argument check uses the standard disk
margin when adequate and a projected Newton enclosure otherwise.

At a disk center \(c\), write the shifted polynomial as
\(\sum_m b_m w^m\), where \(|w|\le r\) for its certified root.
Choose a nonzero midpoint \(s\) of the linear coefficient
enclosure. Then the root displacement satisfies
\[
 w=-b_0/s+e,\qquad
 |e|\le
   \frac{r\,|b_1-s|_{\max}+\sum_{m\ge2}r^m |b_m|_{\max}}
        {|s|_{\min}}.
\]
The linear part \(-b_0/s\) retains the individual grouped scalar
parameters. For consecutive proposed centers \(c,d\), bound
\[
             \Im(\overline c d)+\Im(\overline c\,(-b_{0,d}/s_d))
                              -\Im(\overline d\,(-b_{0,c}/s_c))
\]
from below by projecting the groups before summing. Subtract
\(|c|E_d+|d|E_c+r_c r_d\), where \(E_c,E_d\) are the displacement
error bounds. A strict positive result orders the two root
arguments in the upper half-plane. Conjugation accounts for
the lower roots. The odd-order real disk certifies one negative
real root by uniqueness and conjugation. These disks exhaust
both degrees.

## Exact replay

The [checker](verify_multivariate_terminal_fourteen_shared_low_u.py)
covers twenty-one closed slabs of width \(1/100\) with
14,353 accepted parameter rectangles.

~~~sh
.venv/bin/python research/factorial/image-rigidity/verify_multivariate_terminal_fourteen_shared_low_u.py
~~~
