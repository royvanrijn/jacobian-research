# Mean size of the full affine secant gcd

In a balanced abc chart with small complete-letter multiplicity, the
common divisor of two complete-coordinate differences has negligible
mean logarithmic size at the occurrence endpoint. The estimate counts
all its prime powers, including shallow layers. A one-residue bound
controls overlap with the primes of earlier secants, and random
interpolation controls the least common multiple of the secants.

## A one-residue count at low combined depth

Put \(n=X^L\), with \(L\) fixed before \(X\) tends to infinity, and
write \(d=1/3\). For every \(\eta>0\), a prescribed unit residue
modulo \(D\le X\) of depth at least \(s\), with \(0<s\le d/4\),
contains at most
\[
                  O_{L,\eta,\varepsilon}
                      \bigl(n^{d+\tau-s+\eta}X^\varepsilon\bigr)
                                                               \tag{SGM1}
\]
positive integers \(T\le X\) with \(\operatorname{rad}(T)\le
n^{d+\tau}\). The bound is uniform in the residue and \(D\).
Here \(\tau\ge0\) has a fixed bounded range.

To prove it, fix the prime powers of orders above a fixed truncation
\(J\) at normalized cost at most \(1/(L(J+1))+\eta/4\).
Indeed their radical is at most \(X^{1/(J+1)}\), and their values
on any given support have height-subpower many possibilities.
Choose \(J\) large at the given \(L,\eta\). Box the remaining exact
layer roots, with total normalized cap \(r_0\le d+\tau+\eta/4\)
and arbitrarily small further losses. One may enlarge a cap to
make \(r_0\ge d\). This is an upper envelope, so enlargement
does not affect the represented values. With no varying roots the
count after fixing the coefficient is at most one.

If some cap \(a_i\ge s\), fix every other root. The remaining
powered residue equation has height-subpower many torsion classes,
each containing \(O(1+n^{a_i}/D)\) choices. The cost is at most
\(r_0-s+o_X(1)\). Otherwise each \(a_i<s\). For positive caps take
\(h_i=\lceil s/a_i\rceil\). Then
\[
             h_i a_i\le 2s,\qquad
              \sum_i h_i^{-1}\ \ge\ r_0/(2s)\ \ge\ 2 .
\]
In this step use the depth-\(s\) interpolants of (JRM1).
For \(w=1/\lceil s/a_i\rceil\) their value is \(sw\), since
\(w\le a_i/s\). They are linear of slope \(s\) from zero to
this weight. Assign reciprocal weights up to these bounds
with total exactly two. Thus \(G_2(s;\mathbf a)\ge2s\).
The interpolants and their optimized gains are nondecreasing
in depth, so the same gain applies at the actual depth of
\(D\). Apply (JRM1) with \(\vartheta=1,\xi=0\), using the
\(O(\log X)\) unit-width log bins meeting the complete range.
This gives the same count \(r_0-s+o_X(1)\). Add the initial
coefficient and binning costs within the stated allowance.
[Joint real and residue list moments (JRM1)](ABC_JOINT_REAL_RESIDUE_MOMENTS.md#an-arbitrary-list-of-congruence-and-logarithmic-bins)
· [Uniform full-factor enumeration (EB3–EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## Conditional secant overlap on a prescribed pool

Let \(S\) be a family of \(N\ge n^{d-\tau}\) primitive rows
\(A+B=C\), \(X/2<C\le X\), with
\[
 \operatorname{rad}(ABC)<C^L,\qquad
 \operatorname{rad}(T)\le n^{d+\tau}\quad(T=A,B,C).
\]
Assume the maximum degree of each complete letter is at most
\(n^{u+\tau}\). Identify each row with its complete \((A,B)\)
pair. For distinct rows \(x,y\), put
\[
                g(x,y)=\gcd(|A_x-A_y|,|B_x-B_y|).
\]
This gcd is positive, at most \(X\), and divides the difference
on all three roles. Set \(g(x,x)=1\). For a prime pool \(\mathcal P\)
write \(g_{\mathcal P}\) for the full part supported there.

**SGM2.** Fix \(0<\beta<1/8\), \(0<a<1/3\), and \(b>0\) with
\[
                  (1+b)(1-a)<1,\qquad
                    u<(2/3-2a)\beta .
\]
Take \(\tau\) sufficiently small in terms of the strict margin.
For each fixed anchor \(x\) and each prime pool with at most
\(O((\log X)^{1+b})\) primes, the fraction of \(y\in S\) with
\[
                        g(x,y)_{\mathcal P}\ge n^\beta
\]
is \(O_L(n^{-\delta})\) for a fixed \(\delta>0\) depending on
the margins and sufficiently small fixed \(L\). This holds
uniformly in the anchor and pool. The pool may be chosen from
the anchor and any previously sampled rows. Its choice must
precede the choice of \(y\).

For a selected pair, if a prime \(p\ge n^\beta\) divides
\(g(x,y)_{\mathcal P}\), choose the single modulus \(p\).
The pool gives only height-subpower many possibilities.
At least one original letter of \(x\) is a unit at \(p\), and
its value on \(y\) has the same residue. Apply (SGM1) with
\(s=2\beta/3<d/4\), then multiply by the complete-letter degree cap.

Otherwise, use the prime factors with multiplicity to choose a
divisor \(q\mid g(x,y)_{\mathcal P}\) in \([n^\beta,n^{2\beta}]\).
The first product that reaches \(n^\beta\) suffices. At every
prime of \(q\), at least two of the three original letters
of \(x\) are units. For a role \(T\), let \(q_T\) be the full
part of \(q\) on primes at which \(T_x\) is a unit. Then
\(\prod_T q_T\ge q^2\), so at least one \(q_T\ge q^{2/3}\).
The selected complete letter of \(y\) lies in the prescribed
unit class of \(T_x\) modulo \(q_T\). The choice of this role
can be made from \(x,q\). Apply (SGM1) with \(s=2\beta/3<d/4\),
and multiply by the degree cap. The count for a fixed \(q\)
is at most \(n^{d+u-2\beta/3+2\tau}\), up to arbitrary inner
allowances.

The number of possible \(q\) supported on the pool is at most
\[
 n^{2a\beta}\prod_{p\in\mathcal P}(1-p^{-a})^{-1}
                     = n^{2a\beta+o_X(1)}.
\]
In fact the logarithm of the Euler product is
\(O_a((1+|\mathcal P|)^{1-a})=o(\log X)\), using the order
of the distinct primes. This version of the estimate allows
the pool to have more than sublogarithmically many primes.
Sum over \(q\) and divide by \(N\). The margin
\((2/3-2a)\beta-u\) absorbs \(3\tau\) and the counting
allowances, proving (SGM2). In the single-prime branch the
margin \(2\beta/3-u\) is at least as large.

For a pool supported on a prescribed integer of fixed
polynomial height, the choices of \(q\le n^{2\beta}\) have
height-subpower cost directly. In this fixed-pool case the
same tail estimate holds whenever \(u<2\beta/3\), with \(\tau\)
and counting allowances smaller than the margin.

## Interpolation and the entire gcd

**SGM3.** Under the hypotheses on \(S\) above, for every
\(0<\beta<1/8\) with \(u<2\beta/3\), sufficiently small
\(\tau,L>0\) give
\[
                   \mathbb E_{x,y\in S}\log_n g(x,y)
                               \le \beta+o_X(1).               \tag{SGM3}
\]
The expectation uses independent uniform rows, and conditioning
on distinct rows gives the same bound. The error tends to zero
at each fixed budget. The threshold on \(\tau,L\) may depend on
\(\beta\) and the margin \(2\beta/3-u\). In particular, in the
height-first balanced limit with \(u=0\) the mean tends to zero.
More generally the outer upper bound is \(3u/2\) if \(u<1/12\)
and the balance and size allowances tend to zero.

Choose \(a>0\) sufficiently small that
\(u<(2/3-2a)\beta\). Use a degree
\(r=\lfloor(\log X)^c\rfloor\), where \(0<c<c_*\) in (JSI1)
and \((1+2c)(1-a)<1\). Let
\(m=\lfloor\binom{r+2}{2}/4\rfloor\) and choose independent
uniform anchors \(y_1,\ldots,y_m\), independent of \(x\).
Write \(g_i=g(x,y_i)\), \(P_0=1\), and
\(P_j=\operatorname{lcm}(g_1,\ldots,g_j)\).
Each \(P_j\le X^m\) has \(O(m\log X)=O((\log X)^{1+2c})\)
distinct primes, even though its numerical height cap is
no longer a fixed power of \(X\).

Apply (SGM2) to the next independent anchor, conditionally on
\(x,y_1,\ldots,y_{i-1}\), with pool \(\operatorname{supp}(P_{i-1})\).
Since \(\gcd(g_i,P_{i-1})\le g_{i,\operatorname{supp}(P_{i-1})}\le X\),
its conditional mean logarithm is at most
\(\beta\log n+O_L(n^{-\delta}\log X)\). Use the exact identity
\[
      \log P_m
           =\sum_{i=1}^m\log g_i
                         -\sum_{i=1}^m\log\gcd(g_i,P_{i-1}).
\]
The expectation of \(\log P_m\) is the left side of (JSI1)
with all prime powers at most \(X\). Indeed a proper hit at
\(p^j\) means \(p^j\mid g_i\) for at least one \(i\). Summing
its weight \(\log p\) over the depths counts the log of that
least common multiple exactly. Take \(\tau\) below \(\sigma_*/4\)
as required by the size hypothesis of (JSI1). This gives
\[
            \mathbb E\log_n g(x,y)
               \le \beta+O_L(r/m+n^{-\sigma_*/2}+n^{-\delta}),
\]
because \(\log X/\log n=1/L\). All terms in the error tend to
zero at each fixed \(L\). The diagonal contributes zero, and
conditioning on distinct rows divides by \(1-1/N\).
The final limiting assertions allow \(\beta\) to decrease
after the height limits. No parameter \(\beta\) has to shrink
with height at a fixed budget.
[Random interpolation and weighted proper hits (JSI1)](ABC_SHALLOW_JOINT_INTERPOLATION.md#interpolation-against-proper-congruence-hits)

## Extra gcd in a fixed shallow residue cell

**SGM4.** On the same original radical and complete-degree caps,
suppose every row has a prescribed divisor \(M\) on one marked
complete role \(W\). Fix an absolute residue on each of the other
two complete roles modulo \(M\), as in an actual \(W\)-locked cell.
Write \(\log_n M=\sigma+o_X(1)\). Assume this cell has at least
\(n^{d-\sigma-\tau}\) rows and that
\[
       0\le \sigma<\min\{d/4,\sigma_*/4\},\qquad
                  2u<\beta< d/2-2\sigma .
\]
Here \(\sigma_*\) is the absolute curve gap in (JSI1).
At sufficiently small fixed \(\tau,L>0\), the mean on
distinct independent rows of the cell satisfies
\[
                \mathbb E_{\rm distinct}\log_n(g(x,y)/M)
                                      \le \beta+o_X(1).       \tag{SGM4}
\]
Consequently the limiting upper bound for the mean extra depth is
\(2u\) under the displayed strict range for \(u\). At zero complete
multiplicity it is zero. A cell with a substantial additional
cardinality loss beyond \(\sigma\) is outside this assertion.

For distinct rows set \(g'=g/M\), an integer because all the
coordinate differences are divisible by \(M\). On the diagonal
set \(g'=1\). If the full pool-supported part of \(g'\) reaches
\(n^\beta\), select a divisor \(q'\) just as in (SGM2), taking
a single pool prime if it is at least \(n^\beta\), and otherwise
using \(q'\in[n^\beta,n^{2\beta}]\). Of the two opposite
complete roles, at least one is a unit at each prime of \(q'\).
One of these roles \(T\) therefore has a unit part
\(q'_T\ge (q')^{1/2}\). Both opposite roles are units at every
prime of \(M\), by primitivity and \(M\mid W_x\).
The congruence on \(T\) holds modulo the product \(M q'_T\),
including when the two factors have overlapping support:
the whole product divides \(g\), because \(q'_T\mid g'\).
Its depth is at least \(\sigma+\beta/2-o_X(1)\), and the
modulus is at most \(X\). Apply (SGM1) at this depth, with
\(\sigma+\beta/2<d/4\). After summing over the same pool
dictionary and dividing by the cell size, the strict
margin is \((1/2-2a)\beta-u\), up to chart and count allowances.
The pre-existing loss \(\sigma\) cancels against the contribution
of \(M\) to the one-residue saving.

Choose \(a>0\) small relative to \(\beta/2-u\), and then
the interpolation exponent \(c\) as in (SGM3). The pool overlap
estimate now applies to the successive \(g'_i\). Require
\(\sigma+\tau<\sigma_*/4\) so that the cell itself meets the size
hypothesis of (JSI1), with its original primitive values and
original radical cap. The least common multiple of the modified
\(g'_i\) is at most the least common multiple of the original
\(g_i\). The interpolation bound therefore applies to its
logarithm, and the identical overlap identity proves (SGM4).
The diagonal convention affects only a factor \(1+O(1/N)\)
on passing to distinct rows.
[Absolute residue cells (RB4)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#the-actual-residue-cell-sizes)

## Scope

The previous secant peeling bounds (BC3) and (BC5) control the
part of a secant gcd outside an anchor-dependent pool. They
allow a large part supported on that pool. The deep-collision
estimate (JSI5) treats the average contribution of high layers
of individual primes without a pool. The bound (SGM3) treats
the entire average, including an accumulation of shallow
layers across many primes. It uses the small complete-letter
multiplicity through the uniform conditional pool estimate.

This concerns unconditioned pairs in a near-full family.
A residue cell or a selected relation of power-small density
can still contain much larger gcds. The small mean also
does not bound the radical of the normalized difference
vector. In particular it does not force a saving in the
unrestricted small-budget count.
[Earlier secant peeling and local pools (BC3–BC5)](ABC_BOUNDED_DEGREE_CURVES_AND_SECANT_GCD.md#local-prime-pools-for-common-secant-gcds)
· [Deep secant gcds (JSI5)](ABC_SHALLOW_JOINT_INTERPOLATION.md#deep-secant-gcds-near-zero-complete-multiplicity)
