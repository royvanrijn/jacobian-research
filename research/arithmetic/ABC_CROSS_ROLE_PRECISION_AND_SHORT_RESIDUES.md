# Many-prime precision and short representatives in additive graphs

The same-side difference moments in the three orientations consume one
joint numerical budget when primes are assigned by original letter.
A separate short-term abc estimate applies to high-precision classes
whose residue representatives have small height and small radical.
Actual primitive graphs show that the difference moments can remain
large under the target radical cap alone, even with many fresh primes
at diffuse scaled orders and all three dictionaries at the critical size.

## One difference budget for three orientations

Let \(\mathcal E\) be a finite graph of positive primitive triples
\(a+b=c\le X\), with \(N_W\) distinct complete values on role
\(W\in\{A,B,C\}\). Each pair of labels from different roles gives
at most one row. Write \(n=X^\lambda\). Let \(\mathcal J_A,\mathcal J_B,\mathcal J_C\)
be three sets of prime-depth pairs, with mutually disjoint prime
supports and all primes at most \(n^\beta\).
For \(W\) fixed put \(\mathcal J_{\ne W}=\bigcup_{T\ne W}\mathcal J_T\).
As in (GL5), \(\mathcal W_{\mathcal J}(z)\) sums \(\log p\) over
\((p,j)\in\mathcal J\) with \(p^j\mid z\).
If the set of complete values on each role has width \(B_W\), then
for \(h>2\),
\[
 \begin{split}
 (1-2/h)\sum_{e\in\mathcal E}\sum_{W\in\{A,B,C\}}
                    \mathcal W_{\mathcal J_W}(W_e)
 &\le \sum_W\sum_{x<x'\ {\rm on}\ W}
                       \mathcal W_{\mathcal J_{\ne W}}(x-x')\\
 &\quad+O\!\left(hn^\beta\log X\sum_WN_W\right)\\
 &\le \sum_W\binom{N_W}{2}\log\max(1,B_W)
             +O\!\left(hn^\beta\log X\sum_WN_W\right).
 \end{split}                                                    \tag{CPR1}
\]
Here each pair \(x,x'\) has fixed role \(W\).
The first inequality follows by applying (GL6) with each of the
three roles as target. The signed permutations of \(a+b=c\)
have the same path argument. Use the complete values themselves
as cofactors with fixed factors 1. The disjointness of the prime
supports ensures that a prime-depth contribution from \(x-x'\)
is used at most once across the two orientations with \(W\) on
the summand side. The second inequality is the product bound
on the retained powers in each nonzero difference. The \(O\)-term
is absolute. The component imbalance estimate (GL7) likewise
sums over the same three disjoint difference budgets.
[Path layers and difference moments (GL5–GL7)](ABC_GRAPH_PATH_DIVISOR_LOCKS.md#a-difference-moment-inequality-for-the-actual-edge-layers)

For instance, let \(N_W\le N\), \(B_W\le X\), and
\(N=n^{1/6+o_X(1)}\). For \(0\le\eta<3\), if the tested weight is at least
\((3-\eta)|\mathcal E|\log X\), fix \(0<\beta<1/6\)
and take \(h=n^\theta\) with \(0<\theta<1/6-\beta\).
For every fixed \(\lambda\) and sufficiently large \(X\),
\[
             |\mathcal E|\le
                    \left(\frac{3}{2(3-\eta)}+o_X(1)\right)N^2.
\]
The hypothesis is about the tested weight and the actual colored
graph. The usual fixed-truncation prime coloring preserves count
exponents, with subpower losses. The displayed density bound has no
power-saving implication at relative density \(n^{-o(1)}\).

## Rank tests for a fixed pattern with mixed-role congruences

Let \(x\in\mathbb Z^D\) parametrize an integer additive pattern, and
let a record specify \(\ell_i(x)\equiv b_i\bmod M_i\) for finitely
many integer linear forms. In an actual graph these may be complete
letters, signed differences or forms obtained by substituting its
exact additive equations. For each prime \(p\) and depth \(j\),
let \(r_{p,j}\) be the rank over \(\mathbb F_p\) of the forms
with \(p^j\mid M_i\). For \(1\le k\le D\), put
\[
      P_k=\prod_p p^{\sum_{j\ge1}(r_{p,j}-(D-k))_+}.
\]
Every \(k\)-minor of a matrix of differences of points satisfying
one fixed record is divisible by \(P_k\). In particular, if all
points lie in a box of coordinate widths \(B_1,\ldots,B_D\) and
\[
                     P_k>k!\max_{|I|=k}\prod_{i\in I} B_i ,
 \tag{CPR3}
\]
they lie in an affine subspace of dimension less than \(k\).

Indeed, at \(p\) order the forms by decreasing modulus depth
and greedily keep the forms independent modulo \(p\).
Complete the retained forms to a basis over \(\mathbb Z_p\)
with supplementary forms of depth zero. The count of basis
forms of depth at least \(j\) is \(r_{p,j}\).
After this unimodular coordinate change, each difference has
coordinates divisible by the respective powers of \(p\).
The \(p\)-adic order of each \(k\)-minor is at least the sum of
the smallest \(k\) depths, which is
\(\sum_j(r_{p,j}-(D-k))_+\). Invert the unimodular coordinate
change and use the expansion of minors to obtain the same
divisibility in the original coordinates. The real bound in
(CPR3) is the permutation expansion of each determinant.
If for each active set at \(p\) its rank over \(\mathbb F_p\)
equals its rational rank, the local bound is exact on the
full homogeneous lattice. Each discarded form then belongs
to the integral \(p\)-adic span of chosen forms of at least its own depth,
as follows by solving at a unit minor. The diagonal lattice
in the chosen basis attains the sum of the smallest depths.
The local lattices at finitely many primes can be combined
by the Chinese remainder theorem.

For one four-cycle, with two complete coordinates from each
summand side, the four target forms have rank three together
at every prime. A prime power on the whole cycle supplies
its depth to every two-minor of the differences of valid
parameter vectors. A private prime power on a single edge
supplies no such two-minor precision in dimension four.
More generally, if every tested prime touches just one form,
then \(P_k=1\) for \(k<D\). Deeper private powers of this
kind only strengthen the full-rank determinant bound.
This identifies which overlap is needed to force more than
one rank drop through this divisibility argument.

## A uniform count from a short signed residue

Fix the constants \(h_*,\sigma>0\) in (SH1), decreasing them
if necessary so that \(h_*<1\) and \(\sigma<1/3\).
Fix \(B,q,\xi,\rho\ge0\) and a constant \(K\ge1\).
Consider the distinct integers \(W\in[X/K,X]\) satisfying
\(\operatorname{rad}(W)\le n^B\), for each of which there exist
positive integers \(Q,t\) and an integer \(r\) such that
\[
 \begin{gathered}
      W=Qt+r,\qquad \gcd(Q,r)=1,\qquad 0<|r|\le X^{h_*/2},
      \qquad \operatorname{rad}(|r|)\le n^\rho,\\
      \operatorname{rad}(Q)\le n^q,\qquad \operatorname{rad}(t)\le K n^\xi.
 \end{gathered}
\]
For every \(\eta>0\), all sufficiently small fixed \(\lambda>0\)
and every \(\varepsilon>0\), the dictionary size is
\[
 O_{B,q,\xi,\rho,K,\eta,\lambda,\varepsilon}
 \left(X^{\{\min(\xi,\rho)+(1/3-\sigma)(B+q+\xi+\rho)+\eta\}\lambda
                +\varepsilon}\right).
                                                               \tag{CPR2}
\]
The estimate is uniform over the witnesses, which may vary on
every integer \(W\). When \(r\) is fixed, or is supported on
the primes of one prescribed integer of polynomial height, the
term \(\min(\xi,\rho)\) can be omitted.

Put \(g=\gcd(W,|r|)=\gcd(Qt,|r|)\). Since \(r\) is a unit
modulo \(Q\), \(g\mid t\). The equation between \(W,Qt,|r|\),
after division by \(g\), is a positive primitive additive equation
with largest letter \(C'\) between \(X^{1-h_*/2}/K\) and \(2X\).
The short summand satisfies
\[
   \frac{|r|}{g}\le
       \left(\frac{\max(W,Qt)}{g}\right)^{h_*}
\]
for large \(X\), since \(X^{h_*/2}\le \max(W,Qt)^{h_*}\)
and \(g^{h_*-1}\le1\). Put \(P=B+q+\xi+\rho\).
The radical of the primitive equation is at most \(K n^P\).
For any small \(\eta_1>0\), choose finitely many exponent
intervals \([h_-,h_+]\) covering this range of \(C'\), with
\(h_->0\), \(h_+/h_-\le1+\eta_1\); only sufficiently large
heights need to be covered. On the interval
\(X^{h_-}\le C'\le X^{h_+}\), use the fixed height-relative
budget \(\lambda'=(P+\eta_1)\lambda/h_-\), which is admissible
for large \(X\). Taking \(\lambda\) sufficiently small for
this finite mesh, (SH1) up to height \(X^{h_+}\) has exponent
at most \((1/3-\sigma)(P+\eta_1)(1+\eta_1)\lambda\), up to the
height-exponent allowance. The number of possible \(g\) is at most
\(O_\varepsilon(n^{\min(\xi,\rho)}X^\varepsilon)\), since its radical
is bounded by both \(\operatorname{rad}(t)\) and
\(\operatorname{rad}(|r|)\). For a prescribed support of polynomial
height, supported-power enumeration gives \(X^{o_X(1)}\) choices.
A primitive row and \(g\) recover \(W\) up to the choice of one of
its three letters. Taking \(\eta_1\) sufficiently small gives (CPR2).
[Short-term abc saving (SH1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-strict-saving-when-the-short-term-is-small)
[Supported powers and radical enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

For example, set \(B=1/3\) and \(\xi=1/6\).
For \(q,\rho>0\) sufficiently small in terms of \(\sigma\),
the coefficient in (CPR2) is strictly below \(1/6\) with an
absolute margin after fixing \(q,\rho,\eta\).
This counts original values in a union of classes whose moduli
\(Q\ge X^{1-\xi\lambda}\) have radical at most \(n^q\) and have
short radical-capped unit representatives \(r\), since
\(t=(W-r)/Q\le 2n^\xi\).
The modulus may combine any number of primes and may depend
on the vertex. If all opposite dictionaries have size at most
\(n^{1/6+\eta}\), then edges using these values are power-smaller
than \(n^{1/3}\). The estimate needs a short signed representative
with a small radical, in addition to the aggregate congruence.
Absolute bounds on \(\operatorname{rad}(Q)\) and \(\log Q\)
alone do not impose that short-residue hypothesis.

## A primitive example with the target cap and abundant precision

Set \(c_0=1/6\). Fix an integer \(d\) sufficiently large that
\[
 1/d<c_0/2,\qquad
 (1/3-\sigma)(1/2+2/d)<1/6.
\]
Put \(z=1/d\). For each integer \(R\ge 3\) fix an integer
\(b_R\ge100dR\) divisible by \(R-1\), with \(b_R\to\infty\).
Use the fixed budget \(\lambda_R=1/b_R\), and put
\[
 \mu_R=\frac{z}{4R},\qquad
 \beta_R=\frac{2z}{R},\qquad
 J=\lfloor n^{\mu_R}\rfloor,\qquad
 N=\lfloor n^{c_0-\mu_R}\rfloor,\qquad n=X^{\lambda_R}.
\]
At each fixed \(R\), for every sufficiently large \(X\), there
is a primitive additive graph with all letters in \([X/5,X]\),
all \(c>X/2\), and with the following properties.

- The three complete dictionaries all have size \(\asymp JN\),
  and \(|\mathcal E|\asymp JN^2=n^{1/3-\mu_R+o_X(1)}\).
  The target has genuine cap
  \(\operatorname{rad}(C)\le n^{1/6+2z}<n^{1/3}\).
- At primes at most \(n^{\beta_R}\), using depths at least
  \(d/(2\lambda_R)\), every target has tested prime-power
  contribution at least \((1/2-o_R(1))\log X\) in (GL6).
  The sum of deep difference moments of the two summand
  dictionaries is \(\gg|\mathcal E|\log X\) as well.
  All but a vanishing proportion of this tested edge weight
  occurs in components of size at least \(n^\theta\)
  with both sides within a factor \(n^\delta\), for any fixed
  \(0<\theta<c_0-\mu_R\) and \(\delta>0\).
- Uniformly in a prescribed integer \(P\le X^B\) with \(B\)
  fixed, all but \(O_{R,B}(1/J)\) of the edges have no tested
  prime dividing \(P\). The scaled degree measure of the
  displayed target packets tends to \(z\,\operatorname{Unif}[d-1,d+1]\)
  after the height limit and then \(R\to\infty\).
- For \(R\) sufficiently large, all but a power-small fraction
  of the edges violate both summand caps
  \(\operatorname{rad}(A),\operatorname{rad}(B)\le n^{1/3}\).
  In particular the example concerns the target-capped
  form of the difference-moment hypothesis.

To construct the graph, let
\[
 e_\ell=b_R\left(d-1+\frac{2(\ell-1)}{R-1}\right)\quad(1\le\ell\le R),
 \quad E_0=Rb_Rd,\quad U=(X^{1-c_0\lambda_R}/10^4)^{1/E_0}.
\]
Take \(J\) disjoint \(R\)-element prime sets
\(\{p_{i,1},\ldots,p_{i,R}\}\) in
\([U,(1+1/(10E_0))U]\). At fixed \(R\) this is a fixed proportional
prime interval, so the prime number theorem supplies
\(n^{(1-c_0\lambda_R)z/R+o_X(1)}\) primes.
The exponent is greater than \(\mu_R\).
Write \(M_i=\prod_\ell p_{i,\ell}^{e_\ell}\). These moduli are pairwise
coprime and
\[
  X^{1-c_0\lambda_R}/10^4\le M_i
                \le e^{1/10}X^{1-c_0\lambda_R}/10^4 .
\]
In each block \(0\le i<J\), choose \(A_i\equiv1\bmod M_i\)
and \(B_i\equiv-1\bmod M_i\) to be the least such integers
respectively above \(X/4+iX/(16J)\) and \(X/3+iX/(16J)\).
Take the two lists \(A_i+M_i u\) and \(B_i+M_i v\)
for \(0\le u,v<N\). Since \(M_i(N+1)<X/(32J)\) for large
\(X\), the blocks of either list occupy disjoint intervals.
Their positives sum to \(C=M_i(s_i+u+v)\) with
\(s_i=(A_i+B_i)/M_i\asymp n^{c_0}\). The targets lie in
\((X/2,X]\). Different blocks have no common target:
otherwise their coprime moduli would divide the other
block's quotient, which is \(O(n^{c_0})<M_i\).
The radical of any target is at most
\[
  \operatorname{rad}(M_i)(s_i+u+v)
                       \ll_R n^{c_0+(1-c_0\lambda_R)z}.
\]
This is at most \(n^{c_0+2z}\) for large \(X\).

Retain precisely the primitive edges within each block.
A prime \(p\) dividing both summands cannot divide \(M_i\),
and must divide \(s_i+u+v\). For \(p\nmid M_i\) the two
congruences on \(u,v\) exclude at most \((N/p+1)^2\)
pairs. Sum over \(p\ll n^{c_0}\) to get an upper bound
\[
    N^2\sum_{p}p^{-2}+O_R(N\log n+n^{c_0})
\]
for the nonprimitive pairs in a block. The leading constant is
less than one and the error is \(o_R(N^2)\) in the height limit
since \(2(c_0-\mu_R)>c_0\). Thus there are \(\gg N^2\)
primitive edges per block. At least \(\gg N\) vertices on
each summand side are used, and each target has at most \(N\)
representations per block, so at least \(\gg N\) targets are used.

Use the tests \(p_{i,\ell}\), at every integer depth from
\(\lceil d b_R/2\rceil\) to \(e_\ell\).
For large \(X\) no such power can divide a target in a different
block, since it is coprime to the other block's modulus and
\[
  p_{i,\ell}^{\lceil db_R/2\rceil}
       \ge X^{(1-c_0\lambda_R)/(2R)+o_X(1)}
       \gg_R n^{c_0}.
\]
Every tested layer therefore has precisely the retained edges
from one block. The tested weight on each target in that block
is at least \(((1-c_0\lambda_R)/2+o_X(1))\log X\),
by subtracting \((db_R/2)\log\operatorname{rad}(M_i)\) from
\(\log M_i\). Every distinct pair of used labels on either
summand side in one block has difference divisible by \(M_i\),
so summing over these pairs gives the deep moment lower bound.
In one tested layer the components of size \(<n^\theta\)
have at most \(O(Nn^\theta)\) edges, by summing size times
\(n^\theta\). The components with larger-to-smaller ratio
at least \(n^\delta\) have at most \(O(N^2 n^{-\delta})\) edges.
These errors are \(o_X(N^2)\), uniformly over the tested blocks
and layers.
Since the blocks use disjoint tested primes and each is
of size \(X^{(1-c_0\lambda_R)/E_0+o_X(1)}\), a prescribed
\(P\le X^B\) can meet \(O_{R,B}(1)\) blocks.
The prime weights are asymptotically \((1-c_0\lambda_R)z/R\)
in base \(n\), and the scaled orders form an equally spaced
mesh on \([d-1,d+1]\). The quotient of the target is \(O(n^{c_0})\),
so it contributes no additional positive-scale moment in the
small-budget limit.

Finally every used \(A\) satisfies \(A-1=M_it\) with
\(t\asymp n^{c_0}\), and every used \(B\) satisfies \(B+1=M_it'\)
with \(t'\asymp n^{c_0}\). Thus both lists admit (CPR2) with
\(q=2z,\ \rho=0,\ \xi=c_0\) and \(|r|=1\), for large \(X\).
The portion of either list satisfying its \(n^{1/3}\) cap has
exponent at most \(c_0-\delta_0\) for some fixed \(\delta_0>0\)
by the choice of \(d\). Every vertex is incident to at most
\(N\) edges, so the fraction of qualifying edges is at most
\(n^{-\delta_0+o_X(1)}\) on each side. This verifies the final
property and identifies the omitted radical tests on the example.
In addition the summand blocks are short ordinary progressions
modulo \(M_i\), so the example does not satisfy the regularized
all-scale progression bound (AC1).

At fixed \(R\) all count and component estimates take the height
limit first. Passing to \(R\to\infty\) removes the exponent deficit
\(\mu_R\) while making each tested prime's normalized individual
cost tend to zero. An application to the unrestricted endpoint
still needs to use both summand radicals. One sufficient input
for the high-precision case is a modulus and short residue of
the kind counted in (CPR2), but the target weights and component
moments do not imply that input for arbitrary profiles.
[Regularized complete-letter progression bounds (AC1–AC2)](ABC_DIFFUSE_FACTOR_GRAPH.md#an-all-scale-progression-bound-for-complete-letters)

### Freedom of the residue representatives

In the block construction replace the two starting residues \(1,-1\)
by \(r_i,-r_i\) for an arbitrary unit \(r_i\bmod M_i\).
The height, target cap, primitive-edge lower bound and all tested
moment estimates hold uniformly in these choices. The proof used
only that each residue is a unit and their sum is zero.

Fix \(\alpha,\gamma\ge0\) with \(\alpha+\gamma<1-c_0\lambda_R\).
The number of classes \(r_i\bmod M_i\) with some
\(1\le t\le X^\alpha\) and \(|v|\le X^\gamma\) satisfying
\(t r_i\equiv v\bmod M_i\) is at most
\[
 \begin{split}
 (2\lfloor X^\gamma\rfloor+1)
                  \sum_{1\le t\le X^\alpha}\gcd(t,M_i)
 &\le (2\lfloor X^\gamma\rfloor+1)X^\alpha
                            \sum_{d\mid M_i}1\\
 &=X^{\alpha+\gamma+o_X(1)} .
 \end{split}
\]
The bound uses that one linear congruence has at most
\(\gcd(t,M_i)\) solutions. To bound the sum of gcds,
majorize each gcd by \(\sum_{d\mid t,\ d\mid M_i}d\)
and sum first over \(d\). The modulus has
\(X^{1-c_0\lambda_R+o_X(1)}\) unit classes, since it has
\(R\) prime divisors at fixed \(R\). Thus almost all unit
classes avoid every such short rational representation,
uniformly over the tested blocks.

Also at most \(n^{1/3+o_X(1)}\) classes contain any \(W\le X\)
with \(\operatorname{rad}(W)\le n^{1/3}\), by the supported-power
count. This also bounds the number of classes whose negatives
contain such a value. One can therefore choose every \(r_i\)
so that neither summand list contains a value satisfying
its radical cap and no block class has a short rational
representative of the indicated heights. This shows why the
target cap alone cannot supply the short-residue hypothesis
in (CPR2).

## Relation to inverse-sieve hypotheses

Walsh's Theorem 2.4 assumes that the same set is supported on few
residue classes modulo every prime in a set of positive logarithmic
prime weight \(\sum (\log p)/p\). Croot, Mao and Yip's Theorem 1.11
assumes that the reductions of the same set lie densely in bounded-rank
small progressions for a positive proportion of a specified dyadic
prime interval. A layer in (GL5) only supplies agreement of certain
pairs on the two summand sides at a prime power. The participating
vertices can vary with the prime. In the construction the tested
primes even have \(\sum (\log p)/p=o(1)\) at fixed \(R\), because their
number is \(R J\) and their common size has a strictly larger exponent
than \(\mu_R\). Thus the inverse-sieve hypotheses require a separate
argument even when most tested precision is deep.
[Walsh, Theorem 2.4](https://arxiv.org/pdf/1105.1551)
· [Croot, Mao and Yip, Theorem 1.11](https://arxiv.org/html/2510.08862v2#S1)
