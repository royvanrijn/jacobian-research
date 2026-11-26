# Native prime certificates in growing mixed profile boxes

Let \(H\) be a large positive integer. Choose a nonempty list of
positive integer gaps \(d_i\) and arbitrary orientations for the
corresponding unequal rows. Independently choose each smaller
exponent \(A_i\) uniformly from \(\{H,\ldots,2H\}\). All the
gaps and orientations are held fixed within one such box, and the
list and its length can depend on \(H\). Set
\[
                     D=\sum_i d_i,\qquad d=\max_i d_i .
\]
One may also append a retained ordinary profile fixed within that
box. Its rows may vary with \(H\), with any signed gaps and
orientations. Write \(d_\rho^{\rm abs}\) and \(M_\rho\) for their
absolute gaps and larger exponents, and put
\[
             B_R=\sum_{\rho:\,d_\rho^{\rm abs}>0}
                         d_\rho^{\rm abs}\log(2+JM_\rho),
\]
where \(J\ge4\) is a fixed
[absolute angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Equal retained rows can be dropped.

**Theorem.** For any fixed \(C>0\) and \(0<\varepsilon<1\),
suppose that
\[
 D\le CH,\qquad d\le H^{1-\varepsilon},\qquad
                      B_R\le \exp((\log H)^{9/16}).             \tag{1}
\]
Then a proportion \(1-O((\log H)^{-1/8})\) of the profiles have
disjoint consecutive moment root arguments in every degree.
The implicit constant is uniform in the entire varying gap,
orientation and retained profile data, depending only on
\(C,\varepsilon,J\). In particular the conclusion, with the same
uniform error, holds whenever \(\sum_i d_i^2\le CH\) and the
retained bound in (1) holds. This includes variable row counts
of order \(H\) with genuinely mixed orientations and positive
diffuse curvature.

The [fixed-data density theorem](MULTIVARIATE_NATIVE_PRIME_CHORD_ANGULAR_CERTIFICATE.md)
uses a single variable row when the number of rows and their
gaps are fixed. Here the total number of ramp factors can be of
order \(H\), so one first-block candidate can have many possible
competitors. We count successful primes across the whole box
in that regime, and use a deterministic progression product to
obtain a large prime when the maximum gap dominates.

## Native witness and prime inputs

At each terminal \(3\le N<J\) we seek \(p>N+1\) such that
exactly one variable row contributes a single factor divisible
by \(p\) in its full upper strip
\(\{N A_i+1,\ldots,N A_i+Nd_i\}\). The factor lies among its
first \(d_i\) elements and has valuation one. All other variable
and retained upper strips are \(p\)-units. This certifies disjoint
arguments by the
[native upper chord theorem](MULTIVARIATE_NATIVE_PRIME_CHORD_ANGULAR_CERTIFICATE.md#checking-the-polygon-from-the-factorial-blocks).
Indeed after orienting the marked row increasingly, its sums
in that theorem are \(E_+=1\) and \(E_-=U=V=0\).
One only needs such a prime for each of the finitely many
upper orders. The lowest pair follows by even positivity and
the higher ones follow by the absolute cutoff.

For this terminal let \(\mathcal V\) be the set of primes
dividing the full upper strips of retained rows. Then
\[
                         \sum_{p\in\mathcal V}\log p\le J B_R. \tag{2}
\]
We use the fixed-modulus prime number theorem in arithmetic
progressions, as follows also from
[Thorner and Zaman, Corollary 1.4](https://arxiv.org/pdf/2108.10878)
with fixed modulus (first sum \(\log p\), then use partial
summation). The interval upper bound for primes used below is
[Montgomery and Vaughan, Theorem 2, equation (1.12)](https://personal.science.psu.edu/rcv4/personal/Publications/large_sieve.pdf).
That equation bounds the number of primes in any positive
interval of length \(y>1\) by \(2y/\log y\). It also gives
\(\pi(x)=O(x/\log(2x))\), with a constant adjusted for the
first few integers.

Write \(L=\log H\) and \(W=H+1\), the number of choices for each
base. For any prime \(p>N+1\), conditional on the marked row
and \(p\), another row's chance of one or more upper-strip hits is
\[
                         \le Nd_i(1/p+1/W).
\]
Indeed each prescribed factor is affine in its independently
uniform base, with coefficient \(N\) invertible modulo \(p\).
In the uses below \(p\le O_J(H)\), so the union bound over all
other variable rows is \(O_J(D/p)\).

## Many possible isolated primes

Consider the primes \(p\equiv1\pmod N\) in
\([3NH/2,7NH/4]\) outside \(\mathcal V\). The prime number theorem
gives \(\gg_N H/L\) candidates before the deletions, and (2)
loses only \(O_J(B_R/L)=o(H/L)\) under (1).
Greedily choose primes separated by distances greater than
\(Nd\). One choice discards at most \(O_N(d/\log(2d))\)
additional candidates by the interval upper bound. Thus the
selected set \(\mathcal S\) satisfies
\[
                    |\mathcal S|\gg_N
                    \frac{H\log(2d)}{dL}.                       \tag{3}
\]

For each row and each selected prime \(p\), a strip can contain
no multiple of \(p\) except \(p\) itself. In fact its factors
are positive and at most \(N(2H+d)<2p\) for large \(H\).
Every admissible offset \(h\equiv1\pmod N\) gives a permissible
base \(A_i=(p-h)/N\) in the interior of \([H,2H]\), since
\(d=o(H)\). Thus the full-strip and first-block probabilities are
respectively
\[
                     \beta_i=\frac{d_i}{W},\qquad
                     \phi_i=\frac{\lceil d_i/N\rceil}{W}.
\]
Events that a row's full strip hits two different selected
primes are disjoint, since two factors of that row differ
by less than \(Nd\).

Let \(X_p\) be the indicator that exactly one row has a
first-block hit and all other rows have no full-strip hit.
The expectation of \(X_p\) is the same for all selected primes,
and equals
\[
       P=Z\sum_i\frac{\phi_i}{1-\beta_i},\qquad
                       Z=\prod_i(1-\beta_i).
\]
For large \(H\), \(\beta_*=\max_i\beta_i<1/2\) and (1)
implies \(Z\ge \exp(-2D/W)\), so \(P\gg_{C,N}D/H\).
For distinct \(p,q\) in \(\mathcal S\), the two first-block
hit rows must be distinct, and independence of the bases gives
\[
 \begin{split}
 {\mathbb E}(X_pX_q)
   &=\sum_{i\ne k}\phi_i\phi_k\prod_{j\ne i,k}(1-2\beta_j)\\
   &\le Z^2\left(\sum_i\frac{\phi_i}{(1-\beta_i)^2}\right)^2
      \le (1-\beta_*)^{-2}P^2.
 \end{split}
\]
It follows from the variance bound and Chebyshev's inequality
that the probability of no witness among \(\mathcal S\) is at
most \(1/(|\mathcal S|P)+O(\beta_*)\). By (3) this is
\[
                O_{C,N}\left(\frac{dL}{D\log(2d)}+\frac dH\right).
                                                                  \tag{4}
\]
This proves the required rate if
\(d\ge e^{L^{5/8}}\) and \(D/d>L^{1/2}\), or if
\(d<e^{L^{5/8}}\) and \(D>e^{L^{11/16}}\).

## A wide first block with a small total competing width

Suppose \(d\ge e^{L^{5/8}}\) and \(D/d\le L^{1/2}\).
In a row of maximum gap \(d\), use all first-block factors
of the form \(x_h=NA+h\), \(1\le h\le d,\ h\equiv1\bmod N\).
There are \(k=\lceil d/N\rceil\) such factors.
Put \(P_0=c_{N,\varepsilon}d\log d\), with a fixed sufficiently
small positive constant. Then for large \(H\),
\(\pi(P_0)\le(\varepsilon/4)k\). No prime dividing \(N\)
divides any \(x_h\). For every prime \(p\le P_0\) remove
its contribution from one prime-power-maximal factor.
At each power \(p^e\), the remaining divisible indices
number at most \(\lfloor(k-1)/p^e\rfloor\), because the
original indices form a residue class modulo \(p^e\).
The product of all small-prime powers after the removals
divides \((k-1)!\). Each removed contribution is at most
\(O_N(H)\), as are the \(x_h\). Every \(x_h\ge H\), so if
\(M_{>P_0}\) is the large-prime part of their product,
\[
 \log M_{>P_0}
  \ge k\log H-\pi(P_0)\log(O_N(H))-\log((k-1)!)
  \ge(\varepsilon/2)k\log H.                                  \tag{5}
\]
Here \(k\le d\le H^{1-\varepsilon}\) gives the last estimate.
In particular a prime divisor above \(P_0\) is always present.

The chance that any of these factors is divisible by \(p^2\)
for any \(p>P_0\) is \(O_N(d/P_0)\). For a given \(p\),
the number of possible factor values divisible by \(p^2\)
is at most \(\lfloor C_N H/p^2\rfloor\); each value has at
most \(k\) preimages as pairs \((A,h)\). Sum and divide
by \(W\). On the complement of this exceptional event,
every prime above \(P_0\) occurs at most once in the product
of the marks, since \(P_0>d\) for large \(H\). The expected
logarithmic mass of the primes of \(\mathcal V\) dividing
the marks and exceeding \(P_0\), counted once, is
\[
              O_N\left(\frac{k}{P_0}
                        \sum_{p\in\mathcal V}\log p\right)
                            =O_J(k B_R/P_0).
\]
For this bound a given prime meeting the marks is \(O_N(H)\),
and each of the \(k\) candidate congruences has probability
\(\le1/p+1/W=O_N(1/p)\).
Using (5) and Markov's inequality, the probability of no
prime of valuation one outside \(\mathcal V\) in the chosen
first block is therefore
\[
                   O_{J,\varepsilon}\left(
                        \frac d{P_0}+\frac{B_R}{P_0 L}\right).
\]
Choose a successful prime as a function of only \(A\) and
the fixed retained data. Since \(P_0>Nd\) for large \(H\),
every other factor of the marked row is a unit. Conditional
on this choice, the probability of interference from other
variable rows is \(O_N(D/P_0)\). Under (1) and the two size
conditions starting this subsection, the resulting total
failure bound is \(O_{C,J,\varepsilon}(L^{-1/8})\).

## A narrow first block with a small total width

In the remaining regime \(d<e^{L^{5/8}}\) and
\(D\le e^{L^{11/16}}\), use the first factor \(NA+1\) of a
single variable row. Set \(P_1=e^{L^{3/4}}\) and
\(\sigma=1-1/\log P_1\). By the elementary smooth-number
sum in the fixed-data density proof, the fraction of
possible bases for which this factor has no prime above
\(P_1\) is at most
\[
     \frac{(C_NH)^\sigma}{W}
                 \prod_{p\le P_1}(1-p^{-\sigma})^{-1}.
\]
The logarithm of the product is \(O(\log\log P_1)\).
Indeed \(p^{-\sigma}\le e/p\) for the first-power terms,
and the bound \(\sum_{p\le P_1}1/p=O(\log\log P_1)\)
follows by partial summation of \(\pi(x)=O(x/\log(2x))\).
The higher-power terms are bounded uniformly when
\(\sigma\ge3/4\). The displayed fraction is consequently
\(O_N(e^{-L^{1/4}/2})\). The same multiples-of-\(p^2\)
count as above, with one offset, bounds the square-hit
probability at primes above \(P_1\) by \(O_N(1/P_1)\).
Retained contamination has probability
\(O_J(B_R/(P_1\log P_1))\) by (2), and independent variable
contamination is \(O_N(D/P_1)\). All are \(o(L^{-1/8})\)
under the stipulated bounds. These are all possible
regimes, which proves the theorem.

The proof uses the retained rows solely through (2). In
particular larger retained budgets are allowed within
individual regimes. For (4) one can use \(B_R=o(H)\).
For the wide first-block comparison one can use
\(B_R=O(d L^{3/2})\), and in the narrow first-block
comparison one can use \(B_R\le e^{L^{11/16}}\).
