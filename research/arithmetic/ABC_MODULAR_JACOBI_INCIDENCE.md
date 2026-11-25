# Multiplicative Fourier reduction for three-support incidences

External prime tests with all three radical caps reduce the
unrestricted count to a trilinear correlation of Jacobi sums.
Bounded-order characters give a power-smaller contribution at the
required prime size. For a factorized anatomic box, the remaining
single-letter transforms have useful mean-square bounds, but the
trilinear correlation is still an additional estimate.

Use \(N=X^\nu\) for a fixed \(\nu>0\) and \(Q=N^\beta\).
Average over the primes \(q\in[Q,2Q]\) as \(X\to\infty\).
Exponents and the number of factors in an anatomic box
are fixed before this height limit. The three dictionaries
consist of positive integers at most a fixed power of \(X\);
weighted counts below can also count all presentations,
with subpower multiplicity. In a balanced box take their
sizes or total weights at most \(N^{1+o(1)}\).

## Jacobi sums and the removal of bounded-order modes

For residue weights on the units \(w_A,w_B,w_C\) modulo
\(q\), let \(D_i=\sum_r w_i(r)\) and set
\[
 \widehat w_i(\chi)=\sum_{r\ne0}w_i(r)\overline{\chi(r)},\qquad
 E_i(q)=\sum_{r\ne0}w_i(r)^2.
\]
The weights are nonnegative. Count just unit modular triples
\(a+b=c\). By substituting \(a=uc,\ b=(1-u)c\) and
multiplicative Fourier inversion,
\[
 I_q=\frac1{(q-1)^2}\sum_{\chi,\psi}
 J_q(\chi,\psi)\widehat w_A(\chi)\widehat w_B(\psi)
                     \overline{\widehat w_C(\chi\psi)},\qquad
 J_q(\chi,\psi)=\sum_{u\ne0,1}\chi(u)\psi(1-u).
 \tag{MF1}
\]
The three principal modes give
\((q-2)D_A D_B D_C/(q-1)^2\). All other Jacobi sums
have absolute value at most \(q^{1/2}\), by the Gauss
sum identity when \(\chi,\psi,\chi\psi\) are nontrivial.
If one is principal the bound follows directly by
orthogonality.

For each letter designate at most \(O(1)\) characters
exceptional, including its principal character. The set
may depend on the prime. For example, in the factorized
setting below use the characters annihilated by at least
one of that letter's fixed powers. For a term with
\(\chi\) belonging to the exceptional set for \(A\),
fix \(\chi\), bound its transform by \(D_A\), and use
multiplicative Parseval and Cauchy on \(\psi\):
\[
 \sum_\psi |\widehat w_B(\psi)\widehat w_C(\chi\psi)|
                         \le(q-1)\sqrt{E_B(q)E_C(q)}.
\]
The modes with \(\psi\) or \(\chi\psi\) exceptional have
the same bound in their respective orientation. Thus
the total nonprincipal contribution from at least one
exceptional coordinate is at most a constant times
\[
             Q^{-1/2}\sum_{\rm cyc}
                         D_A\sqrt{E_B(q)E_C(q)}.              \tag{MF2}
\]

For arbitrary fixed integer dictionaries of size at
most \(N^{1+o(1)}\) with subpower weights,
\[
        \mathbb E_q E_i(q)
           \ll N^{1+o(1)}+N^{2+o(1)}/Q.                      \tag{MF3}
\]
The integer diagonal contributes its subpower-weighted
size, and a nonzero difference of polynomial height
is divisible by \(O_\nu(1)\) test primes in the dyadic
interval. The number of those primes is \(\gg Q/\log Q\).
Unit filtering may only decrease this bound.
In particular if \(Q=N^\beta\) with \(\beta>2\), (MF2)
and Cauchy on the prime average give \(O(N^{2-\beta/2+o(1)})\),
which is power-below \(N\). The principal term is
\(O(N^{3-\beta+o(1)})\), also power-below \(N\).

Suppose the candidate domain is a Cartesian product of
three dictionaries. Each dictionary includes its own
height and radical cap. A coupled cap can be covered
by subpower many boxes of individual caps, with any
fixed exponent slack charged in the bound. The identity
does not represent an arbitrary coupled tuple filter as
three separate weights.

Now assume each dictionary satisfies a radical cap \(N\),
so its integers are units at these test primes.
The averaged support-test identity then identifies the
average of the remaining sum in (MF1), with all three
character coordinates outside the exceptional sets,
with the genuine triple count up to a power-small error.
For presentation weights this follows either by testing
each exact tuple or by bounding all presentation
multiplicities by \(X^{o(1)}=N^{o(1)}\).
For the unweighted support test, tuples with one
support have subpower many potential representatives.
After summing the subpower differences in formulating
upper bounds, a uniform bound \(N^{1-\eta}\) for those
averaged trilinear sums supplies a power saving in the
corresponding triple incidence. This requires a
bound on the signed correlation in (MF1).
[Prime-test reduction retaining the three caps](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md#an-average-finite-field-test-equivalent-to-true-occurrence)

## Higher moments on an actual factorized letter

Write a letter as a fixed coefficient times
\(\prod_j x_j^{e_j}\), with the variables \(x_j\) ranging
independently in arbitrary sets \(A_j\subset[1,N^{p_j}]\).
Here \(e_j\ge1\), \(p_j>0\), and the lists are fixed
before the height limit. Bound a variable of zero
exponent cost by fixing it at subpower cost or by using
an arbitrarily small fixed positive cap.
Set \(D=\prod_j |A_j|\), \(p_*=\max_j p_j\), and
\(s=\sum_j p_j\). Assume
\[
                               2\beta>s+p_* .
 \tag{MF4}
\]
In particular \(Q>N^{p_*}\).
Discard the primes dividing the fixed coefficient,
at most \(O_{\nu,\beta}(1)\) for coefficients of
polynomial height. For \(\mathcal G_q\) the characters
with \(\chi^{e_j}\ne1\) for every \(j\), one has
\[
      \sum_{q\in[Q,2Q]}\sum_{\chi\in\mathcal G_q}
          \prod_j
       \left|\sum_{x\in A_j}\chi^{e_j}(x)\right|^2
               \ll Q^2 D\,X^{o(1)} .                        \tag{MF5}
\]
The outer sum is over the retained primes.

To prove this take \(m_j=\lfloor 2\beta/p_j\rfloor\), so
\(N^{m_j p_j}\le Q^2\) and
\(\sum_j 1/m_j\le s/(2\beta-p_*)<1\).
The \(m_j\)-th power of the character sum is an ordinary
Dirichlet polynomial in the products of its \(x\)'s.
Its coefficient square sum is at most
\(|A_j|^{m_j} X^{o(1)}\) by the divisor bound and
its length is at most \(Q^2\). The map
\(\chi\mapsto\chi^{e_j}\) has fibres at most \(e_j\).
For nonprincipal characters at the test primes,
the character large sieve thus bounds the summed
\(2m_j\)-th moment by
\(Q^2 |A_j|^{m_j} X^{o(1)}\).
Hölder on the pairs \((q,\chi)\), with the unused
reciprocal exponent assigned to the constant function
one, proves (MF5).
[Montgomery–Vaughan, Theorem 1 and its prime-character consequence by Gauss sums](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf)

For the large sieve application one can use Theorem 1
in its additive form. For a prime \(q\) and
\(S(t/q)=\sum_m a_m e^{2\pi i mt/q}\), the Gauss transform
and Parseval in \(\mathbb F_q^\times\) give
\[
   \sum_{\chi\ne 1}\left|\sum_m a_m\chi(m)\right|^2
          \le \frac{q-1}{q}\sum_{t=1}^{q-1}|S(t/q)|^2 .
\]
The prime-denominator fractions on the right across
\([Q,2Q]\) are \(1/(4Q^2)\)-spaced.

In a prime-resolved anatomic box the three factorized
sets can be made primitive by assigning primes to
individual slots, with all coefficient supports
treated separately. Within each slot restrict its
allowed integers to its own prime pool. Among actual
primitive triples, averaging over the assignments
keeps a subpower fraction consistent with some fixed
assignment, because each triple uses \(o(\log X)\)
primes and the number of slots is fixed. Every tuple
from the resulting Cartesian product has the
disjoint support property, and its radical budgets
follow from its slot caps. The moment proof holds
for arbitrary fixed allowed root sets. Its
most-primes consequence is uniform in the choice of
a fixed box and assignment, with the set of good
primes allowed to depend on that choice.

## An asymptotic at smaller prime scales for arbitrary cap profiles

Let \(0<\beta<1\) be fixed. For any one factorized
letter with \(s=\sum_j p_j\le1\) and
\(D=N^{1+o(1)}\), define
\[
 \kappa_\beta=\min\left(\frac{\beta}{3},\
                \frac{1-\beta}{2},\ \frac{2\beta(1-\beta)}{3}\right)>0 .
\]
For every fixed \(\epsilon>0\), all but \(o(1)\) of the
test primes satisfy
\[
 \sum_{\chi\in\mathcal G_q}|\widehat w(\chi)|^2
          \ll q D N^{1-\beta-\kappa_\beta+\epsilon}.
 \tag{MF6}
\]
All characters are now extended to vanish at zero, including
the principal character. Discard primes dividing the fixed
coefficient as before. A root can be larger than the prime
in this range.

For any positive integers \(m_j\) with \(\sum 1/m_j\le1\),
the previous proof allows Dirichlet polynomial lengths
beyond \(Q^2\). The large sieve then has factor
\(Q^2+N^{m_jp_j}\), so the bound in (MF5) becomes
\(Q^2 D N^{E+o(1)}\) with
\[
                 E=\sum_j (p_j-2\beta/m_j)_+ .
\]
Take \(m_j=\lceil 1/p_j\rceil\) for every slot of positive
cap. If all \(p_j\le1/2\), then \(m_jp_j\le3/2\) and
\(E\le(1-4\beta/3)_+\). Otherwise let \(p>1/2\) be the
unique larger cap and \(r=1-p\). With these moment orders
the cost is at most
\[
                    (p-\beta)_+
                       +r(1-2\beta/(1+r))_+ .
\]
There is also a bound with \(E\le\max(r,1-2\beta)\), by using just the
second moment on the large slot and bounding the other
transforms trivially. Use this when \(r\le(1-\beta)/2\).
If \((1-\beta)/2<r<1-\beta\), the moment-order bound
saves \(r\min(1,2\beta/(1+r))\) under \(1-\beta\).
If \(r\ge1-\beta\), it gives at most
\(\frac12(1-4\beta/3)_+\). These bounds imply
\(E\le1-\beta-\kappa_\beta\) in all cases. Dividing
by the number of test primes and applying Markov proves
(MF6).

For three dictionaries with these hypotheses, the full
modular incidence, including zero residue classes, is
\[
                I_q^{\mathrm{all}}\sim D_A D_B D_C/q
                   \quad\text{for all but \(o(1)\) of the primes.}
 \tag{MF7}
\]
Use exact tuple multiplicities for the residue weights.
Each slot loses a proportion \(O(Q^{-1}N^\epsilon)\)
to divisibility by \(q\) for all but \(o(1)\) primes,
by counting the bounded number of test prime factors of
each root and applying Markov. Thus the unit populations
in (MF1) are asymptotic to the full populations \(D_i\).
The principal term has the required size. For the
all-good sum use \(|J_q|\le q^{1/2}\) and
\[
 \sum_{\chi,\psi}|F(\chi)G(\psi)H(\chi\psi)|
      \le(q-1)^{1/2}\|F\|_2\|G\|_2\|H\|_2 .
\]
By (MF6), relative to the principal term this is
\(O(N^{-3\kappa_\beta/2+o(1)+3\epsilon/2})\).
For exactly one exceptional coordinate, Cauchy in
the two good transforms gives relative error
\(O(N^{-\beta/2-\kappa_\beta+o(1)+\epsilon})\).
At least two exceptional coordinates determine the
third among boundedly many pairs, and give relative
error \(O(q^{-1/2})\) after excluding the principal
term. Take \(\epsilon\) sufficiently small.

To include nonunits, (MF6) and the \(O(1)\) exceptional
transforms show by Parseval that
\(\sum_{r\ne0} w_i(r)^2\ll D_i^2/q\).
When exactly one dictionary contributes a zero residue,
the other two residues must agree or be negatives.
Cauchy bounds those pair counts by \(O(D_jD_k/q)\),
and the zero mass is \(o(D_i)\). When two residues
are zero, the third is zero as well. The product of
the three zero masses is at most
\(D_A D_B D_C\,O(Q^{-3}N^{3\epsilon})\).
This proves (MF7).

There is also a version with fixed slack. If
\(\beta\in[h,1-h]\) for fixed \(0<h\le1/2\),
the cap sums may be at most \(1+\tau\) and the
populations at least \(N^{1-\tau}\), whenever
\(0<\tau<h^2/24\). Use \(N'=N^{1+\tau}\) and
\(\beta'=\beta/(1+\tau)\) in the moment proof.
The population loss relative to \(N'\) is at most
\(2\tau/(1+\tau)\), which is smaller than
\(\kappa_{\beta'}\ge h^2/3\). Substituting these
quantities in the relative bounds above still gives
power savings in every error term.

This conclusion allows every atomic and diffuse cap
profile in its population range. For exact valuation
slots the integer presentations of each letter have
subpower multiplicity. In a primitive box obtained by
assigning disjoint prime slots, the same is true
with all three radical caps. The true integer count
is at most \(N^{2+o(1)}\) by choosing two letters in
the balanced setting. With fixed slack as above the
same comparison follows by comparing the principal
term with the product of the two largest populations.
Consequently the modular incidences are predominantly
false throughout this range. The size hypotheses apply
to each fixed dictionary used for this theorem.

## The noise exponent for repeated independent prime tests

In the balanced setting of (MF7), draw \(k\) primes
independently with replacement from \([Q,2Q]\), with
fixed \(k\). For a fixed candidate tuple \(t\) with
nonzero integer additive error, write \(p_t\) for its
probability of passing a single test. The divisor bound
gives \(p_t\ll(\log Q)/Q\). The mean single-prime
incidence and the negligible true contribution imply
\(\sum_{\rm false}p_t\ge N^{3-\beta-o(1)}\).
There are at most \(N^{3+o(1)}\) tuples. Convexity gives
\[
 \sum_{\rm false} p_t^k
    \ge \frac{(\sum_{\rm false}p_t)^k}{(N^{3+o(1)})^{k-1}}
    \ge N^{3-k\beta-o(1)} .
\]
The divisor bound gives the matching upper exponent
\(N^{3-k\beta+o(1)}\). This holds with the same false
integer representatives passing all tests. When each
union support has subpower many representatives, the
same exponent holds for support records. Removing
every support that has a true occurrence costs
\(N^{2+o(1)}\) candidate tuples before the single-prime
estimate, since the true count is at most \(N^{2+o(1)}\).
For \(k\beta<2\), the false exponent exceeds one even
when the product of the drawn primes exceeds \(N\).

The lower bound uses independent draws with replacement.
It permits repetition of a prime that divides the
fixed error. There is a separate construction for
distinct primes in a range of total sizes.
[Distinct-prime false incidences](ABC_MODULAR_FALSE_INCIDENCES.md)

## Scope at the prime sizes needed for true occurrence

The averaging reduction needs a test of size
\(Q>N^{2+\delta}\), or several tests of comparable
total modulus size. At that scale (MF2) already removes
bounded-order modes. The mean-square bound (MF5) is
compatible with diagonal-size fluctuations, but in
(MF1) applying absolute values and the convolution
estimate still gives a contribution of order
\(Q^{1/2}(D_A D_B D_C)^{1/2}\) before cancellation.
In particular an independent bound for each of the
three Fourier norms at diagonal strength does not
give the required power saving. The specific missing
estimate is cancellation among the Jacobi phases
weighted by the three high-order multiplicative
Fourier factors of the actual dictionaries. For one
fixed prime exponent \(\beta>2\), an upper bound
\(O_{\lambda,\varepsilon}(N^{1-\eta}X^\varepsilon)\)
for the average signed all-good sum, with absolute
\(\eta>0\) on every balanced endpoint box, would suffice.
The averaged genuine count would then be at most
\(N^{1-\min(\eta,\beta/2-1)+o(1)}X^\varepsilon\).
One can take \(N=X^{(1/3+\tau)\lambda}\) on the
balanced sector and choose its slack \(\tau>0\) small
against this absolute gap, using the existing endpoint
reduction for the other profiles.
The candidate height cutoff is fixed throughout the
prime average. In particular the small-budget limit
is taken only after bounds at each fixed budget.

The Fourier estimates for smooth-number convolutions
in congruence classes address a different population.
Jain's construction of smooth summands bounds the
largest prime factors and uses multiplicative
convolutions across the full smooth families. Here
the product of distinct prime factors must be
bounded on all three letters, and the needed estimate
is an upper bound for a sparse additive incidence.
[Jain, Sums of Integers with Small Prime Factors](https://arxiv.org/html/2611.21136v1)
