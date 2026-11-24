# An indistinguishability bound for rough Eisenstein divisor tests

The pointwise rough-shell theorem supplies the scale of the
signed least-prime inequality. A different approach would avoid
testing leastness by replacing that inequality with a finite
combination of ordinary divisibility tests on the already-rough
shell. The factor-shape model rules out this replacement even
if every such test is evaluated exactly up to modulus \(N^{0.554}\).
At roughness \(N^{0.116}\), even the proposed square-root saving
in the joint error would only reach marked products below
\(N^{0.505028}\).

## The observations

Fix the roughness exponent \(0<\eta\le500/4267\).
A factor profile \(A\) is an unordered list of exponents
\(x_i>\eta\) with sum one. The prime profile is \(P=(1)\).
The [fixed-factor comparison](EISENSTEIN_ROUGH_FACTOR_DISTRIBUTION.md)
and its [fixed-box refinement](EISENSTEIN_DIVISOR_MODULUS.md#factor-shapes-and-the-constants)
give the distinct-value measure \(\mu_\eta\) on profiles.
It has prime mass one and rank-\(j\) density
\[
  \frac{2^{1-j}}{j!\,x_1\cdots x_j}\,dx_1\cdots dx_{j-1},
  \qquad x_j=1-x_1-\cdots-x_{j-1},
\]
in labelled simplex coordinates before passing to unordered
profiles. Its total mass is \(T(\eta)\). On pairs the measure
is \(\nu=\mu_\eta\otimes\mu_\eta\). The normalized shape
comparison divides by \(T(\eta)^2\).

A marked divisor test at exponent \(\theta<1\) sums a bounded
kernel over selected subsets of the two profiles, with the
combined sum of the selected exponents at most \(\theta\).
The kernel may depend on the sizes and assignment of all
selected exponents. The empty selection gives the total
mass. This describes finite sums of divisor marks with
combined modulus at most \(N^\theta\) in the macroscopic
shape limit. Arbitrary fixed-position windows can be
retained by keeping the shape comparison's common
archimedean variable fixed in the construction below.

**Theorem (prime-free factor-shape model with the same divisor moments).**
For all \(0<\eta\le500/4267\) and all \(\theta\le0.554\), there
is a nonnegative measure on the same rough factor-profile pairs
with the same observations as \(\nu\) at exponent \(\theta\)
and zero prime-pair mass. Consequently a pointwise lower
weight for the prime-pair indicator built from these tests
has limiting mean at most zero. The construction also
preserves the first-summand and second-summand counts of
composites with least exponent above any fixed \(b>0.495\).

## Two exact exchanges

It suffices to build the alternative measure at a cap
\(\theta\ge1/2\). Write \(E_\theta\) for the set of
both-composite pairs with minima summing to more than
\(\theta\), and let \(M=\nu(E_\theta)\).
A pair in \(E_\theta\) has no eligible cross marks.
Thus all observations agree between the following two
unit-mass pairs of profiles taken together:
\[
                 (P,P),\ (A,B)
                \quad\longleftrightarrow\quad
                 (A,P),\ (P,B).                              \tag{1}
\]
The same-side observations are unchanged by this exchange.

For the second exchange take a triple \(A=(x_1,x_2,x_3)\)
with every \(x_i<1-\theta\). Its mass among first-summand
profiles is denoted by \(U\). Each two-factor sum in \(A\)
exceeds \(\theta\), and for every \(i\) the companion
semiprime \(S_i=(x_i,1-x_i)\) has exactly one visible factor,
\(x_i\). Consequently
\[
                 2(P,P)+(A,P)
                 \quad\longleftrightarrow\quad
                 (S_1,P)+(S_2,P)+(S_3,P).                    \tag{2}
\]
This preserves the empty mark and every singleton mark.
There are no larger marks on either side. The exchanged
profiles all meet the same roughness cutoff. Interchanging
the two summands gives another copy of (2).

Integrate (1) against \(\nu\) on \(E_\theta\), and integrate
both versions of (2) against the triple-profile measure.
Using source-minus-destination signs, the resulting
signed measure has prime-pair mass
\[
                             L=M+4U.                         \tag{3}
\]
Outside the prime-pair cell the three source layers
\(E_\theta\), \((A,P)\), and \((P,A)\) are disjoint and each
is dominated by \(\nu\). If \(L\ge1\), execute all three
exchanges with factor \(1/L\). The resulting measure is
nonnegative, removes the full prime-pair atom, and leaves
every observation unchanged. The first exchange preserves
the full marginals. The second replaces visible composites
by visible composites or primes by visible composites.
At the cap \(0.554\) all the modified least exponents
in (2) lie below \(0.446\), so the two large-least-factor
counts in the theorem are preserved.

Both \(M\) and \(U\) are nondecreasing as \(\eta\) decreases
or the modulus cap decreases. At \(\eta=500/4267\) and
\(\theta=0.554\), the exact enclosure proves \(M+4U>1.00186\).
Using this construction at \(0.554\) preserves all tests at
smaller caps as well. This proves the theorem.

The same construction applies to the raw coordinate
shape measure by itself. After renormalizing its prime
mass to one, the rank-\(j\) weight is \(2^{j-1}\) times
the distinct-value weight. The two source masses can
only increase. The theorem concerns either weight bank
in isolation. In fact the joint banks can see this particular
exchange. Put \(w(A)=2^{|A|-1}\). In (1) the source-minus-
destination difference in raw mass, in prime units, is
\((w(A)-1)(w(B)-1)\), positive on \(E_\theta\).
In (2) it is \(2+4-3\cdot2=0\). Hence the replacement
of the distinct-value model reduces its raw total by
\[
 \frac1L\int_{E_\theta}(w(A)-1)(w(B)-1)\,d\nu\ \ge\ \frac{M}{L}.
\]
An independently evaluated exact composite marginal likewise
retains information changed in (2).

## Joint raw and value observations

**Proposition (prime-pair mass can vary with both sets of observations fixed).**
At \(\eta=0.116\) and modulus cap \(N^{0.506}\), a nonnegative
replacement of \(\nu\) preserves both the raw and distinct-value
marked moments and decreases the prime-pair mass by more than
\(15\%\). This construction leaves some prime-pair mass.

For \(A\) in the balanced triple layer of (2), put
\(\ell_A=\delta_A+2\delta_P-\sum_i\delta_{S_i}\).
Given any fixed profile on the second side,
\(\ell_A\) annihilates the distinct-value marks on the first
side, including cross selections. Let \(\rho\) be the
probability measure obtained by restricting the semiprime
profile measure to minima in \((\theta-\eta,1/2)\) and
normalizing. This restriction has mass
\[
 W=\tfrac12\log\frac{1-\theta+\eta}{\theta-\eta}.
\]
Every composite in \(\ell_A\) has minimum \(>\eta\), so
there are no cross marks against the support of \(\rho\).
The signed tensor
\[
                  \ell_A\otimes(\delta_P-\tfrac12\rho)        \tag{4}
\]
annihilates the raw tests as well. Indeed the \(w\)-mass
of \(\ell_A\) is \(4+2-3\cdot2=0\), and the \(w\)-mass
of the second factor is \(1-\tfrac12\cdot2=0\).
Same-side raw tests therefore vanish on (4).

The source layers of the integrated tensor are
\((A,P)\), \((P,P)\), and the positive companion-pair
measure \(\tfrac12K\otimes\rho\), where \(K\) is the
first-side intensity of all companions \(S_i\).
Write \(h=1-\theta\).
At a companion minimum \(x\), the allowed interval of
the second coordinate of the triple has endpoints
\(\ell=\max(\eta,1-x-h)\) and \(u=1-x-\ell\).
If nonempty, the ratio of companion intensity to
semiprime density is \(\tfrac12\log(u/\ell)\), by integrating
\(1/(24xy(1-x-y))\) and summing over the three choices
of marked coordinate. Hence the ratio is bounded above
by \(\tfrac12\log(h/\eta)\).
The replay proves
\((3/10)\tfrac12\log(h/\eta)<W\),
so executing the integrated exchange with coefficient
\(3/5\) stays within the available composite-pair mass.
It also stays within the available prime-pair mass,
since \(U\le \log^2(1/\eta)/8\). It removes
\(2(3/5)U>0.15\) of the unit prime-pair mass.
Thus even the combined observations leave a quantified
ambiguity, though the proof of complete deletion above
concerns one set of observations at a time.

## Consequence for the signed pointwise inequality

The comparison from the
[pointwise leastness note](EISENSTEIN_POINTWISE_LEASTNESS.md#signed-reduction-and-the-leastness-budget)
asks for \(R_N-2C_N(B)+Q_N\ge0.001R_N\).
Here \(C_N(B)\) counts first-summand composites visible at
the fixed cutoff \(B\) near \(\sqrt N\), and \(Q_N\) is the
selected cross count with both summands composite.
In the alternative shape model the total remains \(R\),
both marginal invisible bands still have mass at most
\(0.0002R\) if the cutoff is chosen sufficiently near
the half exponent, and the prime-pair atom is zero.
Inclusion–exclusion and \(Q\le\#\{\text{both composite}\}\)
then give \(R-2C(B)+Q\le0.0004R\).
Thus the identical tested moments coexist with different
answers to the required signed inequality. The construction
also leaves the fixed-position distribution intact if the
same exchange is applied at every position.

In the joint-remainder model
\(N^{1/2+o(1)}(B_{\rm mark}D)^{1/2}\), the rough sieve uses an exponent
above \(4.267\cdot0.116=0.494972\). To fit inside a power
saving on \(N\) the mark exponent must be below \(0.505028\).
Even a cap of \(0.554\) permits the exchange. At that cap
the combined exponent would be \(1.048972\), and the
hypothesized square-root dependence gives error exponent
\(1.024486\). A power-law bound \(N^{1/2+o(1)}(B_{\rm mark}D)^\alpha\)
would need \(\alpha<0.5/1.048972\) even to accommodate
that cap. This bound on \(\alpha\) is a necessary
bookkeeping condition for this route, with no sufficiency
claim at longer moduli.

The model concerns factor-shape information. A distribution
estimate that retains leastness after switching to \(pr+qs=N\),
uses arithmetic structure beyond the marked exponent profiles,
or combines different weight banks is outside its scope.
Friedlander and Iwaniec's
[asymptotic sieve](https://arxiv.org/abs/math/9811186) is a
classical example of getting extra information through a
bilinear hypothesis. Here the two exchanges identify what
is missing for these specific Eisenstein marks.

## Exact lower certificates

The [standard-library replay](verify_eisenstein_divisor_exchange.py)
bounds \(M\) by subtracting an upper rectangle sum for the
retained cross pair mass from \((T(\eta)-1)^2\).
The rectangles use lower endpoints to decide membership in
the outer covering. Prefix masses and the Abel formula
give an upper value from the monotone upper CDF envelope.
For \(U\), retain square cells wholly inside the triple
domain in the two free coordinates. The density
\(1/(24xy(1-x-y))\) is log-convex there, so its value at the
midpoint gives a lower integral on each full square by Jensen.
The replay uses exact rational bounds throughout.