# Aggregate deep packet residue blocks

A deep target divisor can be distributed over primes of subpower size.
The residue-block inverse still applies to complete composite divisors.
Its controlling exponent is the empirical price of the divisor, which
can be smaller than the price of its prime support. A list of
almost disjoint divisors can then be compared with the height of one
original summand or a rational root that represents it.

## Composite packet inverse

Write \(c_0=1/6\), \(n=X^L\). The edges are primitive triples
\(A+B=C\) with \(X/2<C\le X\). All three complete dictionaries have
size at most \(n^{c_0+\kappa}\), and each original letter has
radical at most \(n^{2c_0+\kappa}\). Fix
\(0<\beta\le\beta+\gamma<c_0\), \(S>0\) and \(w>0\).
Assign one divisor \(M=M(C)\) to each complete target value.
Call this divisor its packet. Assume
\[
   n^\beta\le\operatorname{rad}(M)\le n^{\beta+\gamma},
   \qquad X^w\le M\le X,\qquad
   v_p(M)=v_p(C)\le S/L\quad(p\mid M).
   \tag{AP1}
\]
There is no lower bound on the size of its individual primes.
For example, suppose \(C\) has normalized radical mass at
least \(\beta\) on primes \(p\le n^\gamma\) with
\(a\le Lv_p(C)\le S\). Greedily take their full powers until the
normalized radical mass reaches \(\beta\). This gives (AP1)
with \(w=a\beta\). The choice depends only on the complete target.

**AP2.** Fix a sufficiently small \(\eta>0\) in terms of the
packet parameters. Take \(0\le\kappa,\tau\le\eta/100\).
For all sufficiently small fixed \(L>0\) and sufficiently large
\(X\) at this \(L\), a family with at least \(n^{2c_0-\tau}\)
edges has a subgraph with at least \(n^{2c_0-\tau-\eta/100}\)
edges and a number \(s\) satisfying
\[
                   \beta/2-\eta\le s\le \beta+\gamma+\eta .
   \tag{AP2}
\]
The subgraph uses at most \(n^{s+\kappa}\) different packets,
each assigned to fewer than \(2n^{c_0-s}\) complete targets.
For each surviving packet \(M\), group its retained edges by
\(r=A\bmod M\), so \(B\equiv-r\bmod M\). The number of
its occupied blocks is between
\(\tfrac12 n^{s-2\eta}\) and \(n^{s+\kappa+\eta}\).
Each occupied block has between
\(n^{2(c_0-s)-3\eta}\) and \(2n^{2(c_0-s)+\eta}\)
edges and between
\(\tfrac12 n^{c_0-s-4\eta}\) and \(2n^{c_0-s+\eta}\)
different values in each of its three roles. A block is an actual
set of original triples with this packet and these residues.

For the proof, first discard targets of degree less than
\(n^{c_0-\eta/10}\). This costs at most
\(n^{2c_0+\kappa-\eta/10}\) edges. Bin packets by the number
\(\nu_M\) of remaining complete targets assigned to them, using
dyadic intervals. Some bin \([H,2H)\) retains at least
\(n^{2c_0-\tau-\eta/200}\) edges for sufficiently large height.
Set \(s=c_0-\log_n H\). Each packet in that bin has at least
\(H\) targets and the complete target dictionary has size at most
\(n^{c_0+\kappa}\). Hence the bin has at most \(n^{s+\kappa}\)
packets. It also has at least
\(\tfrac12 n^{s-\tau-2\kappa-\eta/200}\) packets, by its edge
lower bound, the upper bound on \(\nu_M\), and the fact that
each target lies on at most \(n^{c_0+\kappa}\) edges.

The prime support of \(M\) can be enumerated in at most
\(n^{\beta+\gamma}\) ways. At each fixed support, enumerating
all height-\(X\) values costs \(X^\epsilon\) uniformly for
every \(\epsilon>0\) by fixed-support counting (EB3). This
bounds \(s\) above as in (AP2). To bound \(s\) below fix one
packet and one \(A\)-value. The target radical outside the
packet is at most \(n^{2c_0+\kappa-\beta}\). Apply (FL1)
with fixed letter \(A\), the empty pool on \(B\), and pool \(M\)
on \(C\). With total normalized allowance \(\eta/10\), this
bounds the degree at fixed \(A,M\) by
\(O(n^{c_0-\beta/2+\kappa/2+\eta/10})\).
Summing over \(A\) and all packets in the bin proves the lower
bound in (AP2). This step uses the original radical caps and
uniformity in the pool.
[Fixed-support counting (EB3)](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)
· [Prescribed-pool fixed-letter counting (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

For each residue \(r\) and packet \(M\) let \(a,b\) count
the ambient \(A,B\) values in the classes \(r,-r\), and let
\(e\) count its edges in the bin. Every two complete letters
determine an edge, giving \(e\le\min(ab,\nu_M a,\nu_M b)\).
For \(a>n^{c_0-s+\eta}\) this implies
\(e\le 8 n^{-\eta}\binom a2\) for large \(X\).
Each unordered pair of distinct \(A\)'s enters the sum of such
binomial coefficients at most \(\tau_{\rm div}(|A-A'|)\) times,
where \(\tau_{\rm div}\) denotes the divisor function. This
holds for composite and overlapping packets because the packet
labels are distinct positive integers. The resulting exceptional
edge count is \(O(n^{2c_0+2\kappa-\eta+\eta/100})\) by the
ordinary divisor bound at fixed \(L\). Remove also the
corresponding oversized \(B\)-classes.

If \(a<n^{c_0-s-\eta}\), then at fixed \(M\) the total of \(ab\)
over these classes is at most
\(n^{c_0-s-\eta}|\mathcal B|\), since the \(B\)-classes are
disjoint. Sum over the packet dictionary and similarly on
the other side. This costs \(O(n^{2c_0+2\kappa-\eta})\).
There are now at most \(n^{s+\kappa+\eta}\) remaining blocks
per packet, by the lower bound on each ambient \(A\)-class.
Discard blocks with \(e<n^{2(c_0-s)-3\eta}\).
This costs at most \(n^{2c_0+2\kappa-2\eta}\) edges.
Finally discard packets supporting fewer than
\(n^{2c_0-s-\eta/2}\) surviving edges, at cost at most
\(n^{2c_0+\kappa-\eta/2}\). These losses fit the stated
edge lower bound. The block bounds follow from \(e\le
2n^{2(c_0-s)+\eta}\) and from the fact that the degree at
each vertex inside a block is at most
\(2n^{c_0-s+\eta}\), by pair determination and the
ambient class and target count bounds.
[Single-prime residue-block inverse (RB1)](ABC_DEEP_RESIDUE_BLOCK_INVERSE.md#the-marked-deep-residue-block-theorem)

## Almost disjoint packets at a common summand

**AP3.** Fix an integer \(k\ge1\), \(B_*>0\) and
\(0<\xi<\beta\). Choose \(\eta\) sufficiently small in terms
of \(\xi,\beta,c_0-\beta-\gamma\), and apply (AP2). At each
\(A\)-value prescribe any positive integer \(Q_A\le X^{B_*}\).
For all sufficiently small fixed \(L\), there is a fixed
ordered list of \(k\) packets \(M_1,\ldots,M_k\) which is
admissible for at least \(n^{c_0-4k\eta}\) \(A\)-values.
Admissible means each packet occurs on an edge at \(A\) in
the retained subgraph and
\[
 \prod_{\substack{p\mid M_i\\p\mid Q_A M_1\cdots M_{i-1}}}p
                       \le n^\xi\qquad (1\le i\le k).
   \tag{AP3}
\]
In particular the numerical least common multiple of the
packets after stripping all prime powers on \(Q_A\) is at
least \(X^{k(w-S\xi)}\). The proof also works with fixed
summand role \(B\) in place of \(A\).

Indeed, remove \(A\)-values whose retained degree is less than
\(n^{c_0-\eta/4}\), costing at most
\(n^{2c_0+\kappa-\eta/4}\) edges. At least
\(n^{c_0-\eta/10}\) different \(A\)-values remain for large
\(X\), because each degree is at most \(n^{c_0+\kappa}\).
For any fixed \(A\) and any previous list, put
\(Q=Q_A M_1\cdots M_{i-1}\le X^{B_*+k}\).
An edge at \(A\) violating the next constraint has target
radical outside \(Q\) at most \(n^{2c_0+\kappa-\xi}\).
By (FL1) with that pool on \(C\), their number is
\(O(n^{c_0-\xi/2+\kappa/2+\eta/10})\), uniformly in the
choice of \(Q\). With \(\eta\) small in terms of \(\xi\),
at least half the incident edges remain eligible. Each new
packet accounts for fewer than \(2n^{c_0-s}\) of those edges.
Thus there are at least \(\tfrac14 n^{s-\eta/4}\) choices at
each step. They are distinct from earlier choices because
their whole radical is at least \(n^\beta>n^\xi\).
Multiply these lower bounds and average over the at most
\(n^{k(s+\kappa)}\) lists in the subgraph to get the claim.
For the height assertion, keep from \(M_i\) only its prime
powers outside \(Q_A M_1\cdots M_{i-1}\). The remaining parts
are pairwise coprime, and their product divides the stated
least common multiple. Each deletion costs at most \(X^{S\xi}\)
by the upper order bound in (AP1).

The bounded length \(k\) is fixed before \(L\). The use of
(FL1) is uniform in the actual pools and compatible with
composite packets built from subpower primes. The estimate
gives a polynomial modulus of arbitrarily large fixed height
exponent by taking larger fixed lists. The superpolynomial load
in (RB3) uses separate positive-size primes.
[Single-prime load (RB3)](ABC_DEEP_RESIDUE_BLOCK_INVERSE.md#each-typical-summand-value-obeys-many-independent-residue-constraints)

## Height exclusion and power classes

**AP4.** Fix \(k,w,\beta,\gamma,S\) as above and an effective
price upper bound \(\bar s\) for the bin in (AP2).
Suppose the original \(A\)-values satisfy
\[
 A=\rho t^d,\qquad
 |\mathcal D|\le n^\sigma,\quad \rho\in\mathcal D\subset\mathbb Q_{>0},
 \qquad H(\rho)\le X^{B},\quad H(t)\le X^u,
\]
where \(H\) is reduced rational height and \(d\) is chosen
from a finite list of positive integers at fixed \(L\), before
height tends to infinity. If
\[
                        kw>2u,\qquad \sigma+k\bar s<c_0,
   \tag{AP4}
\]
then no bin with \(s\le\bar s\) can support the (AP2) edge
lower bound. For integer \(t\), use \(kw>u\).
If all bins meeting the (AP2) edge threshold have
\(s\le\bar s\), the family has a uniform power saving
\(O(n^{2c_0-\delta})\) for some \(\delta>0\). For a direct
height exclusion on the original \(A\) take
\(\rho=1,d=1,t=A,\sigma=0,u=1\) with integer roots.

To prove this choose \(\xi>0\) small enough to preserve the
strict height margin after replacing \(w\) by \(w-S\xi\),
then take \(\eta>0\) small enough to absorb all displayed
\(O(k\eta)\) losses. Choose one representation for each
\(A\) before applying (AP3), with \(Q_A\) the product of
the reduced numerator and denominator of \(\rho\) and the
integer \(d\). For sufficiently large \(X\) at fixed \(L\)
this pool has height at most \(X^{2B+1}\). On the resulting
fixed list, each \(A\) uses one of at most
\(n^{k(s+\kappa+\eta)}\) possible residue lists by (AP2).

Fix \(\rho,d\) and one such residue list. Strip the powers
supported on \(Q_A\) from each \(M_i\), obtaining \(M_i'\).
These moduli depend only on the fixed data \(\rho,d,M_i\).
Each \(A\) is a unit modulo \(M_i\) by primitivity, so any
\(t\) is a unit modulo \(M_i'\). At any prime \(p\mid M_i'\)
the equation \(t^d=A/\rho\) has at most \(d\) roots on the
units modulo \(p^v\), since \(p\nmid d\). The total number
of root classes modulo the least common multiple of the
\(M_i'\) is at most \(d^{\sum_i\omega(M_i')}\), hence
\(X^\epsilon\) for each \(\epsilon>0\) at fixed \(L\) and
fixed order list. This follows also from
\(\omega(M)=O(\log X/\log\log X)\) uniformly for \(M\le X\).
Compatible classes combine by CRT even if moduli overlap.
Two distinct rational roots with the same class have a
nonzero cross determinant of absolute value at most
\(2X^{2u}\), divisible by the least common multiple, which by
(AP3) exceeds \(2X^{2u}\) for large \(X\). For integers use
a nonzero difference of size at most \(2X^u\).
Thus all fixed-list \(A\)-values number at most
\(O(n^{\sigma+k(s+\kappa+\eta)+\eta})\), contradicting
(AP3) by the price margin and the choice of \(\eta\).
The bounded fixed order list is absorbed into the inner
height threshold and implied constant. It does not enter
\(\delta\) or the small-budget threshold.

The default price in (AP4) is \(\bar s=\beta+\gamma+\eta\).
A stronger empirical bound can be used if it is uniform
over the bins that carry near-full exponent. For example,
for arbitrarily small fixed packets with compact mean
depth \(w/\beta>3\) bounded away from \(3\), an empirical
price bound \(\bar s\le(1/2+o(1))\beta\) as packet size
decreases would give direct height exclusion. In this
example fix a sufficiently small packet size first, and
then choose an integer \(k\) in the nonempty interval
\(1/w<k<c_0/\bar s\) before the budget and height limits.
For the other profiles the missing input is a stronger empirical
price bound, short root height, or a reconstruction inequality
with weaker hypotheses.

The pair-divisibility count used above is related to
Gallagher's larger sieve. Here the marked divisors have
edge-dependent multiplicities, are not required to be
coprime, and the full packet dictionary is controlled by
the third letter. The step producing lists uses the
separate outside radical cap in (FL1).
[P. X. Gallagher, A larger sieve](https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/18/0/97714/a-larger-sieve)


## A one-letter obstruction to compulsory packet compression

**AP5.** The complete target count, radical cap, a diffuse compact
degree profile, and uniform avoidance of polynomial-height prime
pools are compatible with full packet price \(\beta\).
For each sufficiently small fixed \(L>0\) there exist arbitrarily
large height shells and dictionaries \(\mathcal C\subset(X/2,X]\),
with \(n=X^L\) and \(|\mathcal C|=n^{c_0+o_X(1)}\), such that
\[
 \operatorname{rad}(C)\le n^{2c_0+o_X(1)},\qquad
 \max_{p\mid C}\log_n p=o_X(1).
\]
All valuations \(Lv_p(C)\) lie in a common compact interval in
\((0,\infty)\). Their normalized radical measures have a nonatomic
limit when \(L\to0\) after the height limit. For any fixed
\(0<\beta<c_0\), each packet \(M\) of radical at least \(n^\beta\)
divides at most \(n^{c_0-\beta+o_X(1)}\) members. In particular,
assigning packets adaptively or retaining several choices per
target cannot produce a packet dictionary of price less than
\(\beta-o_X(1)\) while covering \(n^{c_0-o_X(1)}\) targets.

Take an integer \(m\to\infty\) at fixed \(L\), let \(q=2^m\),
and choose \(m\) distinct evaluation points \(\alpha_i\) in
\(\mathbb F_q\). Set
\[
 e_0=\lceil 3/L\rceil,\qquad
 K=\lfloor c_0 Le_0m\rfloor .
\]
Choose positive integers \(e_i\) with exact mean \(e_0\), with
their frequencies tending to uniform on the integer interval
\([e_0-\lfloor1/L\rfloor,e_0+\lfloor1/L\rfloor]\).
Indeed, assign the same number of copies of every integer in that
interval, and place any remaining copies at \(e_0\).
For all small fixed \(L\) this also gives \(1\le K<m\) for
large \(m\).

Use pairwise disjoint lists of \(q\) primes
\(\{p_{i,z}:z\in\mathbb F_q\}\) for \(1\le i\le m\) in
\([qm^4,2qm^4]\). The prime number theorem supplies more than
\(mq\) primes there for large \(m\). For each polynomial
\(f\in\mathbb F_q[x]\) with \(\deg f<K\), form the integer
\[
                  C(f)=\prod_{i=1}^m p_{i,f(\alpha_i)}^{e_i}.
\]
The \(q^K\) integers are distinct, since two polynomials agreeing
on all \(m\) evaluation points are equal. Their logarithmic heights
lie in an interval of length \(O_L(m)\), with
\[
 \log C(f)=e_0m\bigl(\log q+O(\log m)\bigr).
\]
Choose a dyadic shell containing at least \(q^K/O_L(m)\) members.
Its upper endpoint \(X\) satisfies the same logarithmic estimate.
The dictionary exponents in the first display follow, since
\(Le_0\ge3\) and
\[
 \log_n p_{i,z}=\frac{1+o_X(1)}{Le_0m}.
\]
At a fixed \(L\) the normalized radical measures converge
uniformly over the shell members to a uniform grid of total mass
\(1/(Le_0)\). As \(L\to0\) the grid converges weakly to the
measure of total mass \(1/3\) with constant density on \([2,4]\).

If \(M\) divides a code integer, each of its distinct prime
factors specifies an evaluation of \(f\) at a different coordinate.
If \(M\) has \(\ell\) distinct primes, interpolation gives at most
\(q^{\max(K-\ell,0)}\) such polynomials. For
\(\operatorname{rad}(M)\ge n^\beta\) the prime-size estimate gives
\(\ell\ge (\beta-o_X(1))Le_0m\), proving the bound on packet
multiplicity. It holds for all divisors, including the exact full
packets in (AP1). Greedy selection realizes radical mass
\(\beta+o_X(1)\) for any fixed \(0<\beta<c_0\), so the lower
price is attained to the stated precision by fixed-support
enumeration.

The same calculation yields a uniform prescribed-pool bound.
For each \(B_*>0\), \(Q\le X^{B_*}\) includes at most
\(O_{B_*,L}(m)\) of the designated primes. The possible sets from
\(Q\) cost \(2^{O_{B_*,L}(m)}=n^{o_X(1)}\). For
\(0<\xi\le2c_0\) the number of members whose radical on \(Q\) is
at least \(n^\xi\) is therefore at most
\[
                 n^{\max(c_0-\xi,0)+o_X(1)}.
\]
The error is uniform over all pools of the stated height.

The construction concerns the target dictionary alone. It leaves
open whether such packet multiplicities can coexist with near-full
additive incidence and the two other radical caps. Compression,
if forced in that setting, has to use additional hypotheses
beyond the one-letter data above. The interpolation input is the
standard polynomial evaluation code. The earlier construction
fixes the number of slots before the height limit and places a
positive normalized prime mass in each slot. Here that number
grows with height at every fixed budget, so the individual
normalized prime masses tend to zero in the inner limit. The
uniform prescribed-pool bound holds at that same fixed budget.
[I. S. Reed and G. Solomon, Polynomial Codes Over Certain Finite Fields](https://doi.org/10.1137/0108018)
· [Earlier diffuse target construction](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-diffuse-target-family-with-no-cheap-heavy-lock-record)
