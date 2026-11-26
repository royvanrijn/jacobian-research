# Affine collision tests for least-repaired target popularity

A random affine map of two distinct integers into the circle gives
two independent uniform values. Applying it to \(a\) and \(c-b\)
thins actual representations of \(c\) at a prescribed rate and
rejects false candidates at the square of that rate. Extra tests
can require the two values to land in the same bucket without
prescribing the bucket. They reject more false candidates without
further thinning true ones. This gives a bounded occurrence test
for the support multiplicity of a minimally repaired target.

## Actual complete candidates and one randomized test

Write \(N=Y^{1/M}\) at a fixed repair order \(M\). More generally,
allow positive integer representatives \(a,b,c\le N^D\) for fixed
\(D>0\), and at most \(N^{1+o(1)}\) target labels. All eligible
triples are fixed before testing. They can retain the original
complete representatives, pairwise coprimality, the three radical
caps, the exact cutoff \(H_M(c)\le Y\), and the short real window.
Let \(B_c\) count the eligible candidate pairs for a complete
target \(c\), and let \(m(c)\) count its true pairs \(a+b=c\).
For a target radical \(R\) define
\[
 B_R=\sum_{\operatorname{rad}(c)=R} B_c,\qquad
 m_{\rm all}(R)=\sum_{\operatorname{rad}(c)=R}m(c),\qquad
 m_*(R)=\max_{\operatorname{rad}(c)=R}m(c).
\]
An empty class has all three counts zero. A fixed radical permits
\(N^{o(1)}\) complete targets by uniform fixed-support enumeration,
so \(m_*(R)\le m_{\rm all}(R)\le N^{o(1)}m_*(R)\).
For each fixed \(c\), distinct true pairs have distinct complete
\(a\) values, because \(b=c-a\). With fixed squarefree supports
\(R_a,R_b,R\), the number of complete witnesses is \(N^{o(1)}\)
uniformly. Thus the number of occurring pairs of summand supports
at \(R\) and \(m_*(R)\) agree up to \(N^{o(1)}\) factors.
[Weighted support occurrence and the exact least-repair cutoff](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md#exact-weighted-form-on-the-repaired-target-radical-alone)

Take \(H\ge1\). Choose independent Haar \(\alpha,t\) in
\(\mathbb T=\mathbb R/\mathbb Z\), and put
\(h(x)=\alpha x+t\) on integer inputs. Require
\[
             h(a),\ h(c-b)\in[0,1/H).
\]
Independently, take finitely many further affine Haar maps.
For the \(i\)-th such map, split \(\mathbb T\) into \(K_i\) equal
half-open intervals, with integer \(K_i\ge1\), and require the
images of \(a,c-b\) to lie in the same interval. This second kind
of test allows every bucket, with \(V=\prod_i K_i\) and \(V=1\)
when there are no further maps. A test round may also use a
block of distinct random primes from fixed dyadic scales
\([N^{\beta_i},2N^{\beta_i}]\), \(\beta_i>0\), with sum \(\beta\).
Require equality and \(abc\) a unit at each such prime. For
an empty prime block take \(\beta=0\) and omit these conditions.

For \(x\ne y\) the map
\((\alpha,t)\mapsto(\alpha x+t,\alpha y+t)\) is an onto
endomorphism of \(\mathbb T^2\), since its integer determinant
is \(x-y\ne0\). It therefore sends Haar measure to Haar
measure on \(\mathbb T^2\). Distinct inputs pass the first
test with probability \(H^{-2}\), and each further test
with probability \(K_i^{-1}\). Equal inputs pass with
probability \(H^{-1}\) overall. These are the exact affine
sampling and collision identities. [Carter–Wegman, universal hash families](https://research.ibm.com/publications/universal-classes-of-hash-functions--1)

A false complete candidate passes the full round with probability
\[
                   \ll H^{-2}V^{-1}N^{-\beta+o(1)} .
 \tag{HT1}
\]
It has unequal inputs \(a,c-b\). At each optional prime it must
also have \(p\mid a+b-c\), a nonzero integer of polynomial height.
There are only boundedly many such prime divisors at each fixed
positive scale. The prime number theorem bounds the probability
by \(N^{-\beta_i+o(1)}\). Conditioning the fixed number of independent
primes to be distinct has a bounded cost. For the empty block,
the false probability is exactly \(H^{-2}V^{-1}\).

In a fixed complete fibre, the number \(X_c\) of true hits before
optional prime exclusions has
\[
 \mathbb E X_c=m(c)/H,\qquad
 \mathbb E X_c^2=m(c)/H+m(c)(m(c)-1)/H^2 .
\]
This uses pairwise independence on the distinct integer inputs
\(a\) of true witnesses. By Cauchy–Schwarz the probability of a
true hit is at least \(m(c)/(H+m(c)-1)\) for \(m(c)>0\).
Conditional on the affine choices and on a hit, select one
witness by a fixed rule. The optional prime block retains
that witness with probability \(1-o(1)\), uniformly, because the
primes excluding it divide \(abc\), again of polynomial height.
Taking one target of maximal fibre size gives the same lower
bound for the radical event with \(m_*(R)\).

## The bounded target test and its comparison

Admit \(R\) in a round if any eligible complete representative at
that \(R\) passes. Complete representatives and summand pairs
may vary across rounds. Define \(S_j(H,V,\beta)\) to count the
radicals admitted in all \(j\) independent rounds, with fixed \(j\).
For \(h\ge1\) put \(T_*(h)=\#\{R:m_*(R)\ge h\}\). If
\(B_R\le B\) uniformly, then, with \(\kappa_N=o(1)\) uniform,
\[
 \begin{split}
 \sum_R\left((1-\kappa_N)
          \frac{m_*(R)}{H+m_*(R)-1}\right)^j
                &\le \mathbb E S_j(H,V,\beta)\
 &\le T_*(h)+
     O\!\left(N^{1+o(1)}
       \min\{1,\ N^{o(1)}h/H+B H^{-2}V^{-1}N^{-\beta+o(1)}\}^{j}\right).
 \end{split}
 \tag{HT2}
\]
The fraction on the left is defined as zero for \(m_*(R)=0\).
When \(h=1\) the term \(N^{o(1)}h/H\) on the right is
omitted. The same statement for complete targets uses \(m(c)\)
on the left and in the tail, and gives the smaller upper error
with \(h/H\) in place of \(N^{o(1)}h/H\).

Indeed, a union bound over true witnesses gives per-round
probability at most \(m_{\rm all}(R)/H\), and (HT1) and the
candidate union bound give the other term. For \(h=1\) all
targets below threshold have no true witnesses. Use independence
of the \(j\) rounds, not of candidate collisions within a round.
The lower bound uses the same complete target for every round.

In the central least-repair sector refine to individual balanced
summand radical caps \(N^{a}\) with \(a=1+O(\xi)\). On the relevant
original height shells, the real filter in (SI3) supplies
\(B\le N^{a+\eta+o(1)}\) for arbitrary fixed \(\eta>0\).
This counts eligible complete representatives after summing the
subpower choices for the original target on each support.
The real window and every support, height and repair restriction
are deterministic and identical in all rounds.
[Real-filtered candidate bounds (SI3)](ABC_SHORT_INTERVAL_RADICAL_FILTER.md#consequences-for-real-filtered-popularity-tests)

For a fixed \(u>0\), take \(H=N^u\), \(V=N^{\nu+o(1)}\), and
\(h=N^{u-\gamma}\) with \(0<\gamma<u\). If
\[
                   2u+\nu+\beta>a+\eta,
 \tag{HT3}
\]
the lower bound in (HT2) is at least a fixed positive multiple
of \(T_*(N^u)\). The upper error is
\(N^{1-j\mu+o(1)}\), with
\(\mu=\min\{\gamma,2u+\nu+\beta-a-\eta\}>0\).
Take \(j\mu>u+\sigma\) to make the error smaller than
\(N^{1-u-\sigma}\). The same comparison for the target support
pair count follows with arbitrarily small fixed exponent slack
in \(H\) and the tail thresholds. All parameters are fixed
before \(N\) tends to infinity at the given budget.

For occurrence \(u=0\), take \(H=h=1\), so the first bucket
is the entire circle. True pairs always pass the affine
tests, and the error exponent is \(1-j(\nu+\beta-a-\eta)+o(1)\).
In particular one can take an empty auxiliary prime block
and one free bucket test with \(\nu>a+\eta\). Near the
balanced upper threshold \(u=1/2\), with \(\nu=0\), the
auxiliary prime block can have arbitrarily small fixed
\(\beta>0\) as the balance and filter allowances decrease.
At the exact border \(2u+\nu+\beta=a+\eta\) the first-moment
noise bound has no power margin.

## A global noise bound in a shorter real window

If the real filter is \(\lvert a+b-c\rvert\le W\), with \(W\ge1\)
fixed before the affine choices, the probability that *any* false
candidate passes one round is at most
\[
     O_g\!\left(\frac{W}{HV}N^{-\beta+o(1)}\right),
 \tag{HT4}
\]
where \(g\) is the fixed number of free affine gates. For a
nonzero error \(\delta=a+b-c\), the first gate requires
\(\|\alpha\delta\|_{\mathbb T}\le1/H\), and each free gate
requires its corresponding \(\|\alpha_i\delta\|_{\mathbb T}
\le1/K_i\). These are necessary conditions independent of
the translations and the target. Each has probability at
most \(2/H\) or \(2/K_i\), respectively, by uniformity of a
nonzero integer multiple of Haar measure. The optional
prime test has rejection scale \(N^{-\beta+o(1)}\) for this
fixed error. Sum over the at most \(2W\) nonzero integers
in the window. One can therefore replace the false term in
the per-round upper bound (HT2) by its minimum with (HT4).
This does not give the same substitution for the signed
kernel, whose removed terms need not vanish when there are
no gated false hits.

For instance, keep \(W=N^\gamma\) with fixed \(0<\gamma<u\)
and use only the first affine gate with \(H=N^u\). Then almost
every affine round has no false candidate anywhere, and the
gate is one separated product. At \(u=0\), one free bucket
gate with \(K=N^{\nu+o(1)}\), \(\nu>\gamma\), gives the same
global absence of false hits. The short real filter retains
every true pair. An ordinary congruence modulo a prime
\(q>W\) also tests exact equality in this window. The affine
versions show how the real width, thinning and separation
cost interact. They do not yield an upper bound for the
occupied true supports in the reduced window.
[Separation ranks of randomized tests](ABC_RANDOMIZED_ADDITIVE_TEST_SEPARATION_RANK.md)

## A signed test at a fixed occurrence threshold

In a nonempty optional prime block replace each modular
indicator on a unit candidate by its Jacobi kernel with a
prescribed bounded number of coordinate characters removed,
including the principal characters. The kernel differs
pointwise from the exact indicator by \(O(1/p)\),
uniformly in the candidate. Take zero for both on nonunits.
At each complete target \(c\), let \(G(c)\) be the sum over
affine-gated eligible pairs of the product of those kernels.
Use the bounded ramp \(f(x)=\min(1,\max(0,4x-1))\) and give
each radical the maximum of \(f(\operatorname{Re}G(c))\)
over its exact target representatives. Multiply these
quantities across \(j\) independent rounds and sum over radicals.
[Exact kernel resummations and removed modes](ABC_CAPPED_SIGNED_POPULARITY.md#capped-comparisons-at-the-filtered-prime-scale)

On a true candidate unit at every prime, the product kernel
has real part \(1-o(1)\) uniformly. For a false candidate,
at one prime its absolute value is bounded by the exact
indicator plus \(O(1/p)\). Taking the product of these
bounds over the fixed block and expanding shows that its
expected absolute value is \(O(N^{-\beta+o(1)})\): each
exact-indicator factor or error factor at prime scale
\(\beta_i\) costs \(N^{-\beta_i+o(1)}\) by the same
bounded-divisor estimate. The affine gates are independent
and give the factor \(H^{-2}V^{-1}\). Thus the expected
total absolute false signed mass for a radical is
\(O(B_R H^{-2}V^{-1} N^{-\beta+o(1)})\).

If a true unit witness at one representative was gated and
the total absolute false mass is less than \(1/4\), its
ramp is one for sufficiently large \(N\). For a radical
below the lower multiplicity threshold, either a true pair
was gated or the weighted test is zero unless the false
absolute mass reaches \(1/4\). Markov's inequality shows
that under the same strict margin (HT3) the upper comparison
and error exponent of (HT2) still hold. At the high threshold
\(m_*(R)\ge N^u\), the lower bound is a fixed positive
multiple of the number of high targets. The ramp is
at the fixed threshold \(1/2\), independently of \(u\).
There is no estimation of an uncapped high moment of
the original target count. When the prime block is empty
one can use the exact affine-gated occurrence as before.

## Alternative hashes and the remaining analytic count

A uniform residue hash modulo a random prime, or a fixed
number of independent distinct primes of product scale
\(N^\alpha\), gives a modular variant. Require
\(a\equiv r,\ b\equiv c-r\) on that modulus. A false pair
must have all hash primes divide the discrepancy and also
match the uniform hash residue. This costs
\(N^{-2\alpha+o(1)}\). For distinct true \(a\)'s, colliding
modulo the hash product costs \(N^{-\alpha+o(1)}\). The
same second moment gives hit probability
\(N^{\min(0,u-\alpha)-o(1)}\) at true multiplicity \(N^u\).
Using one hash on each summand at distinct primes gives
the same single uniform CRT residue on \(a\) when \(c\)
is fixed. Prescribing independent residues for both
summands at the same modulus instead discards most
complete targets by the compatibility condition. It
does not give an additional independent sample inside
one fibre.

Each free \(K\)-bucket gate is a sum of \(K\) separated
products on the two arguments, so using all the free gates
gives \(V\) terms by direct expansion. The exact rank
comparison for randomized tests describes when that cost
can be reduced.
[Separation rank bounds](ABC_RANDOMIZED_ADDITIVE_TEST_SEPARATION_RANK.md)

For comparison, the mean of the uncapped exact candidate
count with no prime tests is precisely
\(m_{\rm all}(R)/H+(B_R-m_{\rm all}(R))/(H^2V)\).
This follows by summing the single-candidate probabilities.
The needed new upper bound is on the *bounded* test of
actual occupied target supports. These identities and the
noise estimates alone do not give it. In particular the
target-support multiplicity problem remains equivalent
to a power saving in the range of bounded tests above.
[Remaining least-repair multiplicity estimate (FL8)](ABC_FIXED_LETTER_RADICAL_CAPS.md#the-remaining-multiplicity-estimate)
