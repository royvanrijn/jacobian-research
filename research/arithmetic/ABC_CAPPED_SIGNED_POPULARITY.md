# A capped signed test for target multiplicity

The real near-equality filter gives a signed test for the target
multiplicity distribution. Clipping each target contribution keeps
larger multiplicities from dominating a lower popularity threshold.
Fixed polynomial moments face a degree obstruction to implementing
this cap. A separate interpolation example shows that pointwise fibre
bounds allow critical target multiplicity even when every summand
support is used exactly once.

## Capped comparisons at the filtered prime scale

Use \(N\to\infty\) and let the complete integer labels lie in
\([1,N^D]\), with \(D>0\) fixed. The set \(C\) has at most
\(N^{1+o(1)}\) distinct complete targets. For each \(c\in C\),
fix the eligible pairs \((a,b)\) and their nonnegative weights
\(\omega_c(a,b)\) independently of the test primes. Write \(B_c\)
for their total weight, assume \(B_c\le B\), and let
\(m(c)=\sum_{a+b=c}\omega_c(a,b)\). On central balanced shells,
the real filter \(|a+b-c|\le X^\theta\) allows
\(B=N^{1+\eta+o(1)}\) for prescribed \(\eta>0\) and a sufficiently
small fixed \(\theta>0\).
[Real-filtered candidate counts (SI3)](ABC_SHORT_INTERVAL_RADICAL_FILTER.md#consequences-for-real-filtered-popularity-tests)

Sample a prime \(q\) uniformly in \([Q,2Q]\), with \(Q=N^\beta\)
and fixed \(\beta>0\). Count modular pairs in \(I_q(c)\) only
with \(q\nmid abc\). Define \(G_q(c)\) by removing prescribed
sets of at most \(K\) characters from each coordinate of the
Jacobi expansion, including each principal character, with
\(K\) fixed. For arbitrary eligibility conditions, sum the
individual Jacobi kernels over eligible unit pairs. Put both
counts equal to zero if \(q\mid c\).

For a threshold \(H>0\), use the bounded ramp
\[
       f_H(t)=\min(1,\max(0,\,2t/H-1)).
\]
It vanishes for \(t\le H/2\) and equals one for \(t\ge H\).
Let \(T_C(h)=\#\{c\in C:m(c)\ge h\}\). Sample \(j\) independent
primes \(q_i\), with \(j\) fixed, and put
\[
  S_j(H)=\sum_{c\in C}\prod_{i=1}^j
                          f_H(\operatorname{Re}G_{q_i}(c)).
\]
There is \(C_0=C_0(D,\beta,K)>0\) such that, writing
\(\zeta=C_0B\log(2Q)/(QH)\), one has
\[
 (1-j\zeta)T_C(2H)\ \le\ \mathbb E S_j(H)\
       \le T_C(H/4)+|C|\min\{1,4\zeta\}^{j}.                    \tag{CPH1}
\]
These estimates are uniform in the eligibility sets and weights.

Indeed \(I_q(c)=m(c)-U_q(c)+F_q(c)\), where \(U_q(c)\) is
omitted true weight on nonunit triples and \(F_q(c)\) is
false weight from nonzero discrepancies divisible by \(q\).
Every nonzero discrepancy and every product \(abc\) has
\(O_{D,\beta}(1)\) prime divisors in the test range, so the
prime number theorem gives
\[
                 \mathbb E_q(U_q(c)+F_q(c))
                       \ll_{D,\beta} B_c\log(2Q)/Q.
\]
The exact resummations (PM2) and (PM3) give a kernel of
magnitude at most \(1/(q-1)\) on each pair when one character
coordinate is specified. An intersection fixes both indices
and also costs \(O(1/q)\), using the Jacobi bound with the
all-principal term kept separately. Inclusion-exclusion gives
\[
 |I_q(c)-G_q(c)|\ll_K B_c/q,\qquad
 \mathbb E_q|\operatorname{Re}G_q(c)-m(c)|
                    \le C_0 B_c\log(2Q)/Q.
\]
[Exact character resummations (PM2, PM3)](ABC_MODULAR_POPULARITY_NONPRINCIPAL_MODES.md#fixed-target-expansion)

If \(m(c)\ge2H\), the ramp equals one unless the absolute
difference from \(m(c)\) exceeds \(H\). If \(m(c)<H/4\),
the ramp is zero unless the difference exceeds \(H/4\).
Apply Markov's inequality in each case, multiply the
expectations over independent primes, and sum over targets.
For \(0\le\zeta\le1\) use \((1-\zeta)^j\ge1-j\zeta\).
For \(\zeta>1\) the lower bound is automatic.

In particular, with \(H=N^u\), \(B\le N^{1+\eta+o(1)}\)
and \(\beta>1-u+\eta\), one can choose
\(j(\beta+u-1-\eta)>u\) by a fixed margin. The additive
error is then smaller than \(N^{1-u}\) by a fixed power.
Both high and intermediate multiplicities contribute at
most one per target to \(S_j(H)\). In the unweighted count at
\(u=0\), take \(H=1/2\) so that the lower bound includes every
occurring target. For \(u>0\) taking \(H=N^u/2\) puts the true
tail at \(N^u\) in the lower bound, with the same prime-scale
condition. Conversely, the true tail bound at a fixed smaller
exponent than \(u\) bounds \(T_C(H/4)\) in the upper bound.
The exponent slack can be taken smaller than the saving.
[Independent prime blocks (MP1)](ABC_MODULAR_TARGET_POPULARITY.md#support-multiplicities-and-independent-blocks)

### Counting radicals

In the unweighted setting, write \(m^*(R)\) for the largest
\(m(c)\) among the complete targets with radical \(R\), and
use zero for an empty class. One may replace the summand in
\(S_j(H)\) by
\[
 \max_{\operatorname{rad}(c)=R}
             \prod_{i=1}^j f_H(\operatorname{Re}G_{q_i}(c)).
\]
For the lower bound use one representative achieving
\(m^*(R)\). For the upper bound when \(m^*(R)<H/4\), a union
bound over the uniformly subpower many representatives
costs \(N^{o(1)}(4\zeta)^j\). The comparison also holds
if the representative can vary by prime, using the union
bound before independence. Neither version requires a
common summand pair across the blocks.

The number of true summand support pairs at \(R\) differs
from \(m^*(R)\) by at most subpower factors. Each occurring
pair has subpower many exact representatives and the number
of complete target representatives is subpower. Thus (CPH1)
also tests true support multiplicity with arbitrarily small
fixed exponent slack in the threshold. Eligibility may
retain the original primitive quotients and the exact
least-repair cutoff.
[Support occurrence and fixed-support enumeration](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md#an-average-finite-field-test-equivalent-to-true-occurrence)

## A degree obstruction to fixed polynomial moments

A polynomial approximating \(f_H\) within \(1/4\) on
\([0,M]\), \(M\ge2H\), must have degree
\[
                             \gg \sqrt{M/H}.                   \tag{CPH2}
\]
Its supremum norm is at most \(5/4\). Its values at \(H/2\)
and \(H\) differ by at least \(1/2\), so the derivative
has norm at least \(1/H\). Markov's polynomial inequality
bounds the derivative by \(2d^2/M\) times the polynomial norm.

The same lower bound holds when \(H,M\) are integers and the
approximation is required only at integers. The derivative
has norm at least \(1/(2H)\) from the values at
\(\lfloor H/2\rfloor\) and \(H\). If \(d^2>M/2\) the lower
bound already holds. Otherwise compare a point of maximum
absolute value with its nearest integer. Their values differ
by at most half the derivative norm. Therefore
\[
             \|p\|_{[0,M]}
              \le 5/4+(d^2/M)\|p\|_{[0,M]}
              \le 5/4+\tfrac12\|p\|_{[0,M]} .
\]
Apply Markov's inequality again with the resulting norm bound.
[Borwein and Erdélyi, Markov and Bernstein type inequalities](https://people.tamu.edu/~terdelyi/papers-online/subsets.pdf)

Thus at fixed \(u<1/2\), testing \(H=N^u\) uniformly over
the multiplicity range up to \(M\asymp N^{1/2}\) by direct
polynomial approximation takes degree growing as a fixed
power of \(N\). The obstruction concerns this specific
implementation of a uniform cap. A power saving in the
capped signed sum on the actual third dictionary remains
open.

## A pointwise model with full summand dictionaries

Take integers \(s\ge2\) and \(0\le b\le s/2\) and an odd
prime \(q>s+1\). Decompose the polynomials of degree less
than \(s\) over \(\mathbb F_q\) as \(G\oplus H_0\), where
\(G\) consists of polynomials of degree less than \(s-b\),
and \(H_0=t^{s-b}\mathbb F_q[t]_{<b}\). Take the abstract
additive triples
\[
             (x+h,\ x-h,\ 2x),\qquad x\in G,\ h\in H_0 .
                                                               \tag{CPH3}
\]
Each summand label occurs once, since its two components
determine \(x,h\). The two summand parts each have \(q^s\)
labels. The target part has \(q^{s-b}\) labels, each of
degree \(q^b\).

Choose \(s\) distinct nonzero evaluation points \(t_i\).
At fixed \(x\), prescribing the evaluations of \(x+h\)
or \(x-h\) at \(k\) distinct points leaves at most
\(q^{(b-k)_+}\) choices for \(h\). The rank is \(\min(k,b)\),
since \(t_i^{s-b}\ne0\) and the remaining interpolation
matrix is Vandermonde. On the full summand parts the
evaluation vectors range freely over \(\mathbb F_q^s\).

For the integer encodings put \(e_i=s+i-1\),
\(E=s(3s-1)/2\), and \(L=3s/E\). Let \(Q_0=D_s q\log q\)
with \(D_s\) sufficiently large, and fix
\(0<\vartheta<(\log2)/(4E)\). For every slot \(i\) of each
letter \(W\), assign \(q\) distinct primes \(p_{W,i,z}\)
to its \(q\) symbols. All lists are disjoint. Use the
interval \([Q_0,(1+\vartheta)Q_0]\) for the target and
\([2^{-1/E}Q_0,(1+\vartheta)2^{-1/E}Q_0]\) for the
summands. The prime number theorem supplies the lists for
large \(q\) at fixed \(s\). Encode a polynomial \(f\) by
\[
               W(f)=\prod_i p_{W,i,f(t_i)}^{e_i}.
\]
The integer values and radicals are injective on each part
by interpolation and unique factorization. The full
summand parts become full product dictionaries. Encoded
triples are pairwise coprime. For \(X=(1+\vartheta)^E Q_0^E\),
the targets lie in \([Q_0^E,X]\subset (X/2,X]\), and the
summands in \([Q_0^E/2,X/2]\). Their combined radical is
bounded by \(2^{-2s/E}(1+\vartheta)^{3s}Q_0^{3s}
<Q_0^{3s}\le C^L\). The normalized individual radical
masses tend to \(1/3\) at fixed \(s\). After this height
limit and then \(s\to\infty\), each letter has limiting
exact valuation measure \((1/3)\operatorname{Unif}[2,4]\).

Fix a target and a prime pool supported on an integer at
most \(X^{D_0}\), \(D_0\) fixed. Only \(O_{s,D_0}(1)\)
encoding primes can lie in the pool, since each exceeds
\(2^{-1/E}Q_0\). There are boundedly many capture records.
If one free letter has \(k\) captured primes, they fix
\(k\) evaluation positions and leave at most
\(q^{(b-k)_+}\) rows. Since \((b-k)_+\le(s-k)/2\), an
outside radical at most \(X^{L\ell}\) leaves
\(O_{s,D_0}(X^{L\ell/2})\) rows: its \(s-k\) remaining
primes all exceed \(q\). Applying the bound to either
free letter gives the minimum of the two exponents.
Fibres with a summand fixed have one row. This verifies
the full pointwise prime-pool bounds of (FL1).
[Fixed-letter half-radical bound (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

In base \(N=((1+\vartheta)Q_0)^s\), the two summand support
exponents are \(1,1\), and the target support and
multiplicity exponents are \(1-u,u\), with \(u=b/s\).
A fixed repair order at least \(\max_i e_i\) sends the
target to its radical raised to that order, preserving
this base. Any limiting \(u\in[0,1/2]\) is obtained by
choosing \(b/s\) to tend to it. No repeated summand
support can be extracted from this model.

The integer encodings need not satisfy the abstract
addition or the real window in (CPH1). Independently
permuting the target assignments at each slot makes the
expected number of true equalities at most one. For each
edge its integer sum fixes at most one ordered target prime
tuple, assigned to its target label with probability at
most \(q^{-s}\), and there are \(q^s\) edges. The pointwise
model leaves the real-filtered arithmetic problem open.
[Earlier intermediate multiplicity model](ABC_TARGET_MULTIPLICITY_POINTWISE_OBSTRUCTION.md)

## Identical finite moments with different critical tails

Fix \(K\ge1\) and \(0\le u<1/2\). For each
\[
        0<\rho<\min\bigl(1/2-u,\,(1-u)/K\bigr)
\]
there are two pointwise graphs of the type (CPH3) on the same
complete integer dictionaries with all of the following properties
in the height-first limit. Their lists of target multiplicities
have identical power moments of orders \(0,\ldots,K\), with common
total edge count \(N^{1+o(1)}\) in repaired base \(N\).
Each summand support is used at most once in either graph, and both
graphs obey the original pointwise prime-pool bounds.
At a threshold of exponent \(u\), their tail counts have exponents
\(1-u\) and \(1-u-\rho\), respectively. The examples impose
the same integer height and radical caps as in (CPH3).
They do not impose integer equality or the real filter.
\[
 \text{Exact moments through any fixed order can leave
 a power-sized ambiguity in an intermediate target tail.}
                                                               \tag{CPH4}
\]

To construct the multiplicity lists, fix a sufficiently large
integer \(s\) and take \(b=\lfloor us\rfloor\),
\(d=\lfloor\rho s\rfloor\ge1\). Then \(b+d\le s/2\) and
\(Kd<s-b\). Let \(q\) tend to infinity through odd primes,
take \(H=q^b\) and \(R=\lfloor q^d/(K+1)\rfloor\).
Use the \(K+2\) nodes
\[
                     t_0=0,\quad t_1=1,\quad
                     t_{i+1}=iR\quad (1\le i\le K).
\]
Set \(D_0=K!R^K\prod_{i=1}^K(iR-1)\) and
\[
                    w_j=\frac{D_0}{\prod_{\ell\ne j}(t_j-t_\ell)}.
\]
Each \(w_j\) is an integer. The two low-node denominators
have moduli \(K!R^K\) and \(\prod_{i=1}^K(iR-1)\).
The denominator at \(iR\) has modulus
\(R^K i!(K-i)!(iR-1)\), which also divides \(D_0\).
In particular \(|w_1|=K!R^K\); for \(j\ge2\),
\(|w_j|\asymp_K |w_1|/R\).
For every polynomial \(P\) of degree at most \(K\),
Lagrange interpolation on the \(K+2\) nodes and comparison
of coefficients in degree \(K+1\) give
\[
                        \sum_j w_j P(Ht_j)=0 .                 \tag{CPH5}
\]
Thus the positive and negative parts of these integer
weights give identical moments at multiplicities \(m_j=Ht_j\).

To realize the lists, for each \(j\ge1\) allow a separate copy of the
polynomial model (CPH3) with \(b\) high coordinates if \(j=1\),
and \(b+d\) if \(j\ge2\). View the copies in one additive group
by appending tags \(j,j,2j\) on the three roles in a fixed cyclic
group of odd order greater than \(2(K+1)\).
Each target of the first copy allows
at least \(H\) incident edges, and each target of the other
copies allows \(q^{b+d}\ge m_j\). The number of available targets
in the first copy is \(q^{s-b}\), and in each other it is
\(q^{s-b-d}\). Choose
\[
                       A_0=\left\lfloor
                  c_K q^{s-b}/|w_1|\right\rfloor
\]
with a sufficiently small positive \(c_K\).
Since \(Kd<s-b\), one has \(A_0|w_1|\asymp_K q^{s-b}\).
The \(j\)-th copy can supply \(A_0|w_j|\) distinct targets,
by the capacity bounds and \(|w_j/w_1|\ll_K 1/R\) for \(j\ge2\).
Choose precisely that many targets for the graph of the sign
of \(w_j\), and keep an arbitrary \(m_j\) of each one's edges.
There are no other edges. Thinning edges in a copy cannot
increase any fixed-target prime-pool fibre. Summand roles retain
their one-use property by disjointness of the copies.

Embed all copies in the same three integer dictionaries using
disjoint prime palettes for each role, slot and copy.
The palettes fit the prime intervals in (CPH3) after increasing
its constant \(D_s\) by a factor depending on \(K\).
Both graphs have identical target dictionaries, including unused
labels. One can pad this dictionary to exponent one using
unused polynomials of degree less than \(s\) in its existing
palettes. Unused labels have multiplicity zero. Thus moments of
positive order match by (CPH5) and the zeroth moment matches
because the target set is the same. At threshold \(H\) the sign
of \(w_1\) supplies \(\asymp_K q^{s-b}\) targets. The other sign
has \(\asymp_K q^{s-b-d}\) targets by the weight ratios and
the fact that some node of the other sign is positive.
The first graph has \(\asymp_K q^s\) edges, using the
\(H\)-node for the lower bound and all nodes for the upper.
Both edge counts are identical by the first moment.
The numerical height, radical, coprimality and pool arguments
of (CPH3) apply within every copy. At fixed \(s\) the repaired base
is \(q^{s+o_q(1)}\). Taking that limit first and then
\(s\to\infty\) gives the asserted exponents, including occurrence
at \(u=0\) where \(H=1\). Equal power moments also give equal
factorial moments of the same orders.

This example separates tail counts while the total count remains
at the critical exponent in both graphs. A proof of first-moment
saving on the actual integer graph would directly prove the desired
count. The ambiguity concerns inferring a saving at one prescribed
tail level from finitely many moment values and the pointwise
conditions. The moment matching by interpolation here is classical.
The additional property is simultaneous realization under the
balanced pointwise arithmetic restrictions.
[Prékopa, The discrete moment problem and linear programming](https://doi.org/10.1016/0166-218X(90)90068-N)
· [Wu–Yang, Chebyshev polynomials, moment matching, and optimal estimation of the unseen](https://arxiv.org/abs/1504.01227)
