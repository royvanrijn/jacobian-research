# Quadratic signatures of private complete letters near the test scale

Private odd valuations can supply a growing set of actual complete
letters whose quadratic signatures balance at almost every test
prime. A second moment over all integers and Burgess's inequality
put the test scale just above the square root of the individual
radical cap. The sign balance does not require a prime-by-prime
choice of the letters.

## A uniform interval second moment

Let \(n=X^L\) at fixed \(L>0\). Fix \(t>0,\beta>t/2\) and let
\(Q=n^\beta\). Take \(m\) positive integers \(V_i\le X\) with
canonical squarefree kernels \(D_i\le n^t\), where \(V_i/D_i\)
is a square. Assume the \(D_i\) are distinct. Write \(\chi_i(j)\)
for the primitive quadratic Dirichlet character of
\(\mathbf Q(\sqrt{D_i})\), using the trivial character for \(D_i=1\).
For any complex numbers \(w_i\) of absolute value at most one,
\[
   \sum_{Q<j\le2Q}\left|\sum_{i=1}^m w_i\chi_i(j)\right|^2
          \ll Qm+Qm^2 n^{-\delta},                         \tag{QS1}
\]
for some \(\delta>0\) depending on \(\beta,t\). Uniformity also
holds for \(\beta\) in a fixed compact subset of \((t/2,\infty)\).

Indeed, the diagonal terms cost \(O(Qm)\). For \(i\ne k\) the
product character is nonprincipal, because a distinct pair of
squarefree kernels has nonsquare product. Its modulus is at most
\(64n^{2t}\), by multiplying the two fundamental discriminants.
Write it as a primitive quadratic character of conductor \(D\)
times the exclusion of the extra prime factors of that modulus.
Separate the possible component of conductor dividing eight by
progressions modulo eight. The remaining primitive conductor is
odd and squarefree. In each progression, Möbius inversion on the
extra odd prime factors reduces the sum to \(n^{o(1)}\) interval
sums of the primitive character, with lengths at most \(Q+O(1)\).
Indeed each divisor with nonempty contribution is at most \(2Q\).
If a divisor is not coprime to the primitive conductor the
contribution is zero. Otherwise multiply its character value
out and re-index the resulting progression by the inverse of
eight modulo the odd conductor. If the primitive odd conductor
is one, the original primitive component modulo eight is
nonprincipal. In that case the same inclusion-exclusion and its
bounded complete-period sums give \(n^{o(1)}\) directly.

For any fixed integer \(h\ge2\), Burgess bounds each of the
remaining interval sums of length at most \(Q+O(1)\) by
\[
                  O_{h,\epsilon}
           \left(Q^{1-1/h}n^{\,2t(h+1)/(4h^2)+\epsilon}\right).
\]
Choose \(h\) large enough to have
\(\beta>t(h+1)/(2h)\). Choose the exponent allowance
\(\epsilon\) smaller than the resulting fixed saving. The
factorization count of the extra divisors is subpower and fits
inside that saving. This proves (QS1).
[Burgess, On character sums and L-series II](https://doi.org/10.1112/plms/s3-13.1.524)
· [Trevino, explicit cubefree form](https://arxiv.org/abs/2511.17778)

At an odd prime \(p\), the character \(\chi_i(p)\) equals the
Legendre symbol \((V_i/p)\) unless \(p\mid V_i\) but \(p\nmid D_i\).
A single \(V_i\le X\) has \(O_{L,\beta}(1)\) primes in the test
interval. Hence for any fixed \(\xi>0\), at most
\(O_{L,\beta,\xi}(1)\) primes there have more than \(\xi m\)
nonunit values \(V_i\). For the remaining comparison, (QS1)
and Chebyshev bound the number of primes where
\(\left|\sum_i w_i\chi_i(p)\right|>\xi m\) by
\(O_\xi(Q/m+Qn^{-\delta})\). The prime number theorem therefore
gives a relative exceptional fraction
\[
 O_{L,\beta,t,\xi}\left(\frac{\log Q}{m}+n^{-\delta}\log Q\right).
                                                            \tag{QS2}
\]
If \(m/\log X\) tends to infinity, almost every prime has
\(m/2+O(\xi m)\) square and \(m/2+O(\xi m)\) nonsquare
members when \(w_i=1\). The signs here concern actual unit
values. Assigned fixed weights satisfy the corresponding
small-bias estimate.

## Actual private selections at saturation

In a central abc shell let \(n=X^L\), and assume a family with
\(\operatorname{rad}(abc)\le c^L\) lies at the balanced \(1/3\)
limiting count scale, with each complete radical at most
\(n^{1/3+\kappa}\) for small fixed \(\kappa>0\).
Fix \(0<u<2/15\). The weighted square-ray bound
discards at power saving the triples whose odd-valuation
radical on any pair has mass at most \(u\). Use the shell count
exponent margin and then take sufficiently small fixed budgets
before the height limit. This is uniform over the three pair
choices.
[Strict odd-mass bound on each pair](ABC_WEIGHTED_SQUARE_RAY_CLASS_MOMENTS.md#exclusion-at-the-balanced-endpoint)

On the remaining near-full family select private packets of
\(m\gg(\log X)^C\) rows for some fixed \(C>1\), using (GP2)
with per-role loss less than \(u/4\). Each row offers a letter
with odd-valuation mass at least \(u/2\). Its retained packet
contains an odd-valuation prime other than two, since the
normalized loss is less than \(u/4\). Choose that complete
letter as \(V_i\). In the sampling proof of (GP2), primes
shared with any other original sampled row are removed, so
this prime is absent from all other chosen values even
before their stripping. Thus the selected squarefree kernels
are distinct, and each is at most \(n^{1/3+\kappa}\).
Equations (QS1)–(QS2) apply above
\(\beta>(1/3+\kappa)/2\).
[Private packets and shared-prime removal](ABC_GROWING_PRIVATE_PRIME_PACKETS.md#drawing-the-private-packets)

There is a fixed-letter analogue. Assume the fixed-letter
near-half hypotheses of (GP2), with free radical caps
\(n^{1/3+\kappa}\) and a fibre of size at least \(n^r\)
where \(r\) is close to \(1/6\). Apply (AF1) to the two
free affine forms as functions of one variable. Their supports
are disjoint by primitivity. Outside the fixed-letter prime
pool their joint odd-valuation mass exceeds \(r-\eta\)
on all but a power-smaller fraction. Use this lower bound
in place of \(u\) in the preceding selection. The same
character bounds hold inside the original fixed-letter fibre.
[Uniform pair power-defect bounds](ABC_GROWING_AFFINE_DEFECT_FIBRES.md#growing-families-of-affine-forms)

If the global family also has a fixed positive mass of odd
valuations on primes \(\ge n^b\) on every retained row, one can
use the polynomial-size private selection in (GP2) on that
range. Choose \(\kappa_1>0\) below \(b\) and the family count
exponent, and keep \(m\gg n^{\kappa_1}\) rows, losing less than
a fixed fraction of the stated odd mass. Select one complete
letter of each row using a private odd prime in that range.
Then (QS2) gives a power-saving exceptional fraction. The
polylogarithmic version does not assume mass at polynomially
large primes.

In the central least-repair regime without a large escaped
target tail, \(N=n^{1/3}\) up to the balance tolerances.
The limiting prime threshold on that scale
for this quadratic test is \(Q>N^{1/2+\varepsilon}\).
This includes scales of interest for high target popularity.
It controls the distribution of one quadratic mode across a
sample of actual complete values. The capped signed Jacobi
test depends on the simultaneous third-letter constraint and
arbitrary eligible pairs. Equation (QS1) does not bound that
test or the weighted residue moment of the full odd
dictionary.
[The capped signed Jacobi comparison (CPH1)](ABC_CAPPED_SIGNED_POPULARITY.md#capped-comparisons-at-the-filtered-prime-scale)
