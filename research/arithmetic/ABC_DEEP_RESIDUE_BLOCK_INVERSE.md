# Deep residue blocks with all three radical caps

A marked prime of logarithmic size \(\beta\) on the target of a
near-minimum complete-dictionary graph creates blocks of dictionary
exponent \(1/6-\beta\). This applies to the actual residue classes at
a prime power of fixed positive numerical height. Almost all marked
edges occur in blocks that are nearly quadratic in their three local
dictionary sizes. Apart from the depth requirement, the argument
works with arbitrary valuations of the three letters.

## The marked deep residue block theorem

Put \(c_0=1/6\), let \(n=X^L\), and fix
\[
  0<\beta\le c_0,\qquad 0<w\le1,\qquad
  0<\eta<\beta/30,\qquad
  0\le\kappa,\gamma,\tau\le \eta/1000.
\]
All these parameters are fixed before the height limit. Let \(\mathcal E\)
be a family of primitive positive triples \(A+B=C\), \(X/2<C\le X\),
with each of the three complete dictionaries of size at most
\(n^{c_0+\kappa}\) and each individual radical at most \(n^{2c_0+\kappa}\).
The argument applies also to signed permutations of the three roles,
with the shell restriction on the largest original letter. In that
case use the signed values of the two summands in the residue equations.

On every edge choose one prime \(p\) with
\[
 n^\beta\le p\le n^{\beta+\gamma},\qquad
 P_p:=p^{\lceil w\log X/\log p\rceil}\mid C .
\]
The marks may vary arbitrarily between edges, including edges with the
same target. For fixed \(p\), partition the marked edges by the unit
residue \(r=A\bmod P_p\). The \(B\)-residue is \(-r\), and the target
residue is zero. Call the edges with mark \(p\) and residue \(r\) a
block. Write \(m=n^{c_0-\beta}\).

**RB1.** Uniformly for all sufficiently small fixed \(L\), one can
discard \(O(n^{2c_0-\eta/40})\) edges, after which the following hold.

- Every nonempty block has at least \(m^2 n^{-3\eta}\) edges.
  Each of its three actual original dictionaries has size between
  \(m n^{-4\eta}\) and \(m n^\eta\). In fact the upper bounds on the
  two summand dictionaries hold for their entire ambient residue
  classes, and the full target dictionary has at most \(m n^{\eta/4}\)
  entries divisible by \(p\).
- Each surviving prime has at least \(n^{2c_0-\beta-\eta}\) surviving
  edges in at least \(n^{\beta-3\eta}\) distinct residue blocks. In
  each summand dictionary at least \(n^{c_0-7\eta}\) different
  original values occur on those blocks, all in at most
  \(n^{\beta+\kappa+\eta}\) unit residue classes modulo \(P_p\).

If the initial count is at least \(n^{2c_0-\tau}\), then there are
\(\gg n^{\beta-\eta/4-\kappa-\tau}\) surviving primes. The implied
constant in the exceptional count is uniform in the graph, prime
marks and dictionaries, at fixed parameters and \(L\).

To prove the target bound, let \(f_p\) be the fraction of the
distinct ambient target values divisible by \(p\), before any
thinning. Adopt the prime prices
\[
   q_p=\max(f_p,1/p),\qquad
   w_X(p)=-\log_X q_p .
\]
The prices satisfy the uniform partition sum condition (EB3).
Indeed \(\sum_p f_p\) is the average number of prime divisors of one
target value and is \(O(\log X/\log\log X)\).
If more than \(m n^{\eta/4}\) ambient targets are divisible by a
marked \(p\ge n^\beta\), then \(p f_p>n^{\eta/4-\kappa}\).
An edge with such a mark has its total effective radical price
discounted by more than \((\eta/4-\kappa)L\).
The raw budget on the whole triple is \((1+3\kappa)L\) in base \(X\).
Apply (EB4) at the padded budget \((1+4\kappa+\eta/1000)L\), which
is a valid relative budget with respect to the largest original
letter in the shell for sufficiently large \(X\).
The effective price cap of these edges is at most
\((1+3\kappa-\eta/5)L\). With the allowance in (EB4) chosen to be
\(\eta/200\), and the height allowance chosen to be
\(\varepsilon=\eta L/200\), the total count of these edges is
\(O(n^{2c_0-\eta/40})\). The same prices account for all the marked
primes at once. No sum over primes is needed in this step.
[Effective radical counting with arbitrary prime prices (EB3–EB4)](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)

Now let \(a_{p,r},b_{p,r}\) be the ambient counts in the two
summand residue classes, and write \(e_{p,r}\) for the marked
edge count in the block. The number \(C_0(p)\) of ambient targets
divisible by \(p\) is at most \(m n^{\eta/4}\).
For each residue \(r\) we have
\[
 e_{p,r}\le
 \min\{a_{p,r}b_{p,r},\ C_0(p)a_{p,r},\ C_0(p)b_{p,r}\}.
 \tag{RB2}
\]
Each inequality uses uniqueness from the corresponding pair of
original complete letters. In particular, a class with
\(a_{p,r}>m n^\eta\) accounts for at most
\(4 n^{-3\eta/4}\binom{a_{p,r}}2\) marked edges,
and the symmetric statement holds for \(b_{p,r}\).
A pair of distinct ambient summand values belongs to such a common
residue class for at most \(1/w\) different marked primes. The
product of the corresponding pairwise coprime \(P_p\) divides
their nonzero difference, which has absolute value at most \(X\).
Summing the pair counts, all oversized blocks together account
for \(O_w(n^{2c_0+2\kappa-3\eta/4})\) edges.

For a block with \(a_{p,r}<m n^{-\eta}\), bound its edge count by
\(m n^{-\eta}b_{p,r}\). At fixed \(p\), the classes on the right
are disjoint, so sum this bound over \(r\) and use the analogous
bound for small \(b_{p,r}\). There are at most \(n^{\beta+\gamma}\)
primes in the window even with the integer count in place of a
prime count. Altogether the small blocks contain at most
\(2n^{2c_0+\kappa+\gamma-\eta}\) edges. For each \(p\) there are
at most \(n^{\beta+\kappa+\eta}\) remaining blocks, by the
ambient lower bounds on their side sizes. The blocks having
fewer than \(m^2 n^{-3\eta}\) edges now account for at most
\(n^{2c_0+\kappa+\gamma-2\eta}\) edges across all \(p\).
In each retained block all three actual dictionary sizes are
at most \(m n^\eta\), and pair uniqueness bounds each from below
by its edge count divided by \(m n^\eta\).

Finally discard primes with fewer than \(n^{2c_0-\beta-\eta}\)
edges left, losing at most \(n^{2c_0+\gamma-\eta}\) edges.
Each remaining block has at most \(m^2 n^{2\eta}\) edges, so
at least \(n^{\beta-3\eta}\) blocks survive for such a prime.
The two residue class lists are injectively paired, giving the
stated summand occupancies. By (RB2) each prime can carry at most
\(C_0(p)n^{c_0+\kappa}\) edges. Near the full scale the
exceptional count is negligible and this upper bound gives the
number of surviving primes. This proves (RB1).
[Graph path congruences and the same-side transfer](ABC_GRAPH_PATH_DIVISOR_LOCKS.md#a-difference-moment-inequality-for-the-actual-edge-layers)

### The upper end of the prime-size range

The estimates in (RB1) include \(\beta=c_0\). In that case
\(m=1\), and the lower bounds on the sizes and edge counts
within a nonempty block are automatic. The bound on oversized
classes and the lower bounds on the number of occupied classes
and primes still apply.

More generally, if an original role \(T\) has complete
dictionary size at most \(n^{a_T}\), then any prime of \(T\) has
empirical frequency at least \(n^{-a_T}\) in that dictionary.
With prices formed from the maximum of the three frequencies
and \(1/p\), an original row gains discount at least
\[
           L\sum_T\sum_{p\mid T}(\log_n p-a_T)_+ .
\]
The partition sum stays sublogarithmic by adding the three
frequency sums. Thus, under the raw cap \(n^{1+3\kappa}\),
a sum bounded below by a fixed \(\xi>3\kappa\) gives a
strict saving by (EB4), for sufficiently small fixed \(L\).
In particular, primes of size at least \(n^{c_0+\xi}\) are
power-small on minimum dictionaries with sufficiently small
cap and dictionary tolerances, whether or not the primes
are deep. The argument bounds the upper end of the prime-size
range and gives no fixed margin at \(p=n^{c_0+o(1)}\).

## Each typical summand value obeys many independent residue constraints

For a surviving summand value \(A\), put
\[
 Q(A)=\prod_{\substack{p:\ \text{a surviving edge at }A\\
                                 \text{is marked }p}}P_p ,
\]
and define \(Q(B)\) in the other role. For a given marked prime
all edges at a fixed \(A\) lie in one of its residue blocks,
and there are at most \(C_0(p)\le m n^{\eta/4}\) such edges.
The number of prime factors of \(Q(A)\) is at most
\(\log Q(A)/(w\log X)\). The same holds at \(B\).

**RB3.** At most \(2w^{-1}n^{2c_0+\kappa-7\eta/4}\)
surviving edges have an endpoint with
\(Q(A)<X^{n^{\beta-2\eta}}\) or \(Q(B)<X^{n^{\beta-2\eta}}\).
In fact an endpoint with this bound is incident to at most
\[
 w^{-1} n^{\beta-2\eta}\,m n^{\eta/4}
           =w^{-1}n^{c_0-7\eta/4}
\]
marked edges, and each whole summand dictionary has at most
\(n^{c_0+\kappa}\) values. Thus a near-full marked family
has almost all of its edges at both endpoints subject to a
superpolynomial product of distinct moduli. More generally, the
same bound works with any list of residue records containing all
the incident marked constraints at the vertex. An actual integer
can represent those residues, and no lower bound on the least
representative follows merely from the size of the product.

## A uniform saving from a power class on one summand

The residue inverse can use a power grid on just one original
summand. Both other letters may have arbitrary valuations.
It is enough that the coefficients of the powers form an affordable
dictionary. They need not be shared divisors of all the values.

**RB4.** Fix \(\beta,w\) as in (RB1), an integer \(k\ge1\), and
constants \(u,\sigma,B\ge0\) such that
\[
                 kw>2u,\qquad \sigma+k\beta<c_0 .
\]
On one summand, say \(A\), suppose every value in the marked graph
can be represented as
\[
              A=\rho t^d,\qquad
              \rho\in\mathcal D\subset\mathbb Q_{>0},\qquad
              H(\rho)\le X^B,\quad H(t)\le X^u,\quad
              |\mathcal D|\le n^\sigma ,
\]
where \(H\) is the maximum of the absolute reduced numerator and
the positive reduced denominator. The order \(d\) can be chosen
from any finite list of positive integers determined at fixed \(L\) before the height
limit. If \(t\) is always an integer, the condition \(kw>u\)
suffices. Under these hypotheses, for sufficiently small fixed
\(L\) the marked graph has \(O(n^{2c_0-\delta})\) edges for some
\(\delta>0\), with \(\kappa,\gamma\) sufficiently small in terms of
the fixed margins. The small-budget threshold and \(\delta\) are
independent of the finite list of orders. The inner height
threshold and implied constants may depend on that list.

To prove it, choose \(\eta\) small in terms of the fixed margins
and \(\tau=\eta/1000\), with \(\kappa,\gamma\le\eta/1000\).
Suppose the original count is at least \(n^{2c_0-\tau}\).
Apply (RB1), then remove the exceptional edges of (RB3).
At least \(\gg n^{c_0-\tau-\kappa}\) different \(A\)-values remain,
since any one has degree at most \(n^{c_0+\kappa}\).
Each such value had, in the (RB1) graph, at least
\(n^{\beta-2\eta}\) incident marked primes. This follows since
\(\log P_p\le\log X\) at every mark, while
\(\log Q(A)\ge n^{\beta-2\eta}\log X\).
Use the prime incidence before the (RB3) removals, where the
residue occupation bounds are valid. For each \(A\) choose its
representation. Remove from its incident list the
\(O_{B,L,\beta}(1)\) band primes dividing the numerator or
denominator of its chosen coefficient.
Each value still gives \(\gg_k n^{k(\beta-2\eta)}\) lists of
\(k\) distinct admissible primes. There are at most
\(n^{k(\beta+\gamma)}\) ordered lists in total, so one fixed
list is shared by at least
\[
         \gg_k n^{c_0-\tau-\kappa-2k\eta-k\gamma}
\]
different values.

By (RB1), the reductions of these values at the \(k\) moduli
lie among at most \(n^{k(\beta+\kappa+\eta)}\) residue lists.
At fixed \(\rho,d\), all the roots \(t\) are units at these
primes. For each prescribed unit value of \(A\) modulo \(P_p\),
the congruence \(t^d\equiv A/\rho\bmod P_p\) has at most \(d\)
unit residue solutions once \(p>d\), by reduction modulo \(p\)
and unique Hensel lifting at each simple root. Since the list
of orders is fixed at the height limit, this holds for all
marked primes for large \(X\). Combine the residues by CRT.
Two different rational \(t\)'s in the same class have a
nonzero cross determinant of absolute value at most \(2X^{2u}\)
divisible by \(\prod P_p\ge X^{kw}\), a contradiction for large
\(X\). For integer roots use their difference at most \(X^u\).
Allowing all coefficients and orders gives at most
\[
                    O\bigl(n^{\sigma+k(\beta+\kappa+\eta)}\bigr)
\]
values on the fixed list. The implied constant includes the sum
of \(d^k\) over the finite order list. This contradicts the
lower bound by the strict margin in \(\sigma+k\beta<c_0\),
so one can take \(\delta=\tau\).

For example, suppose \(A=D a^d\) with integer \(D,a>0\),
\(d\ge d_0>1/w\) from a finite list at fixed \(L\), and
\(\operatorname{rad}(D)\le n^{\sigma_0}\) for
\(\sigma_0+\beta<1/6\). The coefficient need only satisfy this
radical bound. Uniform support enumeration gives a dictionary
of size \(n^{\sigma_0+\epsilon}\) for any fixed
\(\epsilon>0\), and \(a\le X^{1/d_0}\). Taking \(k=1\) gives
a strict count saving in this sector. This permits a power
grid that varies between rows, without a grid hypothesis on
the other original summand or the target.
[Uniform fixed-support enumeration (EB3)](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)

### Rational carriers in minimum dictionaries

The same (RB4) conclusion holds for rational carriers on the
minimum-dictionary graph. Replace \(A=\rho t^d\) by
\[
                       A=\rho\,\Phi(t),\qquad H(t)\le X^u,
                       \quad t\in\mathbb P^1(\mathbb Q).
\]
Use the same coefficient dictionary and allow \(\Phi\) from a
finite list of nonconstant rational maps at each fixed \(L\).
Their degrees and list size are fixed before the height limit.
Represent each map as \([P:Q]\) by coprime homogeneous integer
forms of the same degree, with coefficients bounded by a fixed
power of \(X\) at that \(L\). The maps themselves may vary with
\(X\) within these bounds. The condition \(kw>2u\) applies to
projective parameters. For integer affine parameters \(kw>u\)
suffices.

Fix a representation for each \(A\). For fixed \(\rho,\Phi\), only
boundedly many actual \(A\)-values have \(A/\rho\) a critical value
of the map in characteristic zero. The finite list and the
coefficient cap therefore bound the number of these \(A\)-values
by \(O(n^\sigma)\). Each has degree at most \(n^{c_0+\kappa}\),
so their incident edges have a fixed power loss because
\(\sigma+k\beta<c_0\). For every remaining \(A\), write
\(A/\rho=a_0/b_0\) in lowest terms, and exclude the band primes
dividing the numerators or denominators of \(\rho\) and \(A/\rho\),
the homogeneous resultant of \(P,Q\), or the binary discriminant
of \(b_0P-a_0Q\). This discriminant is nonzero by the exclusion of
exact critical values. In degree one the resultant condition alone
ensures a nonzero linear fibre at the retained primes. These
excluded primes divide a polynomial-height integer at fixed \(L\),
so only \(O(1)\) of them lie in the band for each \(A\), uniformly
in the maps and values.

The prime-list pigeonhole proof of (RB4) now goes through after
this bounded removal. At a fixed retained prime, coefficient, map
and marked unit residue \(R\) of \(A\), any corresponding parameter
solves the projective congruence
\[
                         P-R\rho^{-1}Q\equiv0\pmod {P_p}.
\]
The fibre polynomial reduces modulo \(p\) to a nonzero scalar
multiple of \(b_0P-a_0Q\) for any eligible lift. Its projective
roots in reduction are simple. The resultant condition prevents
both map coordinates from vanishing, so a root above a unit value
is not a pole. In either chart with unit denominator Hensel's
lemma gives a unique lift modulo \(P_p\) above each projective
root. There are at most \(\deg\Phi\) such classes. A list of \(k\)
classes therefore fixes a rational parameter by the same cross
determinant argument as before. For affine integer parameters,
use a difference bounded by \(2X^u\). The finite list and the
degree factors cost constants at fixed \(L\), and the exponent
comparison is unchanged.

### Bounded multi-grids without minimum-dictionary assumptions

**RB5.** A single deep prime already works when \(A\) is a
monomial in boundedly many roots of height \(n^{O(1)}\).
Fix \(r\ge1\) and \(C_1,B\ge0\). Suppose the values on one
summand have representations
\[
 A=\rho\prod_{i=1}^r x_i^{d_i},
 \qquad x_i\le n^{C_1},\qquad \gcd(x_i,x_j)=1\ (i\ne j),
 \qquad H(\rho)\le X^B,
\]
with the same coefficient dictionary cap \(n^\sigma\) as in
(RB4). The positive degree vectors may range over any finite
list at each fixed \(L\). Suppose that the opposite complete
summand dictionary has at most \(n^b\) entries. There is no
minimum-dictionary assumption on the other roles and no
restriction on valuations of the opposite summand or the
marked letter.
For any fixed \(a>0\), at sufficiently small fixed \(L\), the
marked edges for which the chosen \(\rho\) is a \(p\)-unit number at
most
\[
                       O(n^{b+\sigma+\beta+\gamma+a}) .
 \tag{RB5}
\]
The small-budget threshold depends on \(r,C_1,\beta,\gamma,w,a\)
and is independent of the degree list at the height limit.
If \(\rho\) is always integral, every mark has the unit property
by primitivity. For rational coefficients one may treat the
remaining marks separately, for instance by (RB3) in its
minimum-dictionary regime. A marked edge at which its coefficient
is not a unit uses one of only \(O_{B,L,\beta}(1)\) possible
band primes per value. This does not give (RB5) for these edges
at arbitrary dictionary sizes.

Indeed, at fixed opposite summand, \(p\), coefficient and degree
vector, the remaining product satisfies a fixed unit congruence
modulo \(P_p\). Put \(j_p=v_p(P_p)\).
Each free root is at most
\(P_p^{C_1/(\beta j_p)}\), and
\(j_p\ge w/((\beta+\gamma)L)\).
Apply (FC2) with a fixed amplification order \(\ell\) large
enough that \((\beta+\gamma)r(r-1)/(2\ell)<a/2\).
Its cost is at most \(p^{r(r-1)/(2\ell)} X^\varepsilon\),
where the small height allowance can be chosen as a fixed
fraction of \(aL\). The hypotheses of (FC2) hold because the
original roots are pairwise coprime and are units at \(p\).
The degree list contributes a constant at the height limit.
Sum over the opposite dictionary, the coefficient records,
and the at most \(n^{\beta+\gamma}\) possible marked primes.
This proves (RB5); it does not require the block inverse
for the edges on which the coefficients are units.

For example, take an arbitrary bounded number of exact valuation
grids on \(A\), using disjoint prime supports for their roots. If
the complement has radical at most \(n^{\sigma_0}\), it gives
integer coefficients with dictionary exponent
\(\sigma_0+\epsilon\) for any fixed \(\epsilon>0\).
The roots have height at most \(n^{2c_0+\kappa}\) under the
original cap when they are products of the primes at specified
exact degrees. Thus the whole marked sector saves below \(1/3\)
if \(b+\sigma_0+\beta+\gamma<1/3\) with a fixed margin.
This condition applies also away from minimum dictionaries.
[Amplified monomial release over full moduli (FC2)](ABC_OMITTED_COORDINATE_GRAPH.md#signed-monomial-congruences-over-full-moduli)

## The remaining lift problem

For fixed positive \(\beta,w\), a sufficient new estimate for this
sector is a power saving in the number of primes \(p\in
[n^\beta,n^{\beta+\gamma}]\) for which two actual low-radical dictionaries
support the simultaneous occupancies from (RB1) at \(P_p\).
The relevant residue classes must match additively and each
block must contain its own equally small dictionary of actual
target values with the original third cap. An upper bound of
\(n^{\beta-\sigma}\) such primes, with \(\sigma>0\) uniform
at sufficiently small fixed \(L\), contradicts (RB1) if the
tolerances are chosen small compared with \(\sigma\).

The pair-divisor argument is the collision step of the larger sieve.
At these occupancies its saving is limited: order \(n^\beta\) primes
with \(n^{\beta+o(1)}\) occupied classes at each prime use order one
in the normalized height allowance of each pair. Consequently
a further constraint on the distribution of high residue lifts
or on the simultaneous matched blocks is needed. The result applies
to marked positive-size prime windows in minimum complete dictionaries.
For primes of normalized logarithmic size tending to zero, the
effective price discount and the occupancy margins above also tend
to zero. General balanced profiles need a separate extraction or count.
[Gallagher's larger sieve](https://www.impan.pl/en/publishing-house/journals-and-series/acta-arithmetica/all/18/0/97714/a-larger-sieve)
· [Compact-prime packing in a minimum dictionary](ABC_MINIMUM_DICTIONARY_COMPACT_PRIME_PACKING.md)
