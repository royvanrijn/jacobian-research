# Growing affine fibres and the radical of weighted differences

A growing list of distinct shifted values of a query has disjoint
fresh prime support once the primes of its fixed pair determinants
are removed. The determinant product may have quasi-polynomial
height. Its number of distinct primes, rather than that height,
is the relevant input. This gives simultaneous bounds for every
power order and for arbitrary weighted products of the shifts.

For a prime pool \(\mathcal P\), a nonzero rational \(z\), and
\(q\ge2\), write
\[
 {\cal R}_{\mathcal P}(z)=\prod_{\substack{p\notin\mathcal P\\v_p(z)\ne0}}p,
 \qquad
 {\cal D}_{q,\mathcal P}(z)=
       \prod_{\substack{p\notin\mathcal P\\q\nmid v_p(z)}}p .
\]

## Growing families of affine forms

Fix \(B\ge1\), \(r>0\), and \(0<\eta<\min(1,r)\).
Set \(\tau=\eta/(64B)\), fix \(0<L\le1\) before the height
limit and put \(n=X^L\). Let
\[
                   2\le m\le 2(\log X)^{\tau L}.
\]
Take pairwise nonproportional affine forms
\(F_i(x)=a_i x+b_i\), \(1\le i\le m\), with integer coefficients
of magnitude at most \(X^B\) and nonzero slopes. The forms can vary
with height. Assume \(\mathcal Q\) is a set of integer queries on
which \(0<|F_i(x)|\le X^B\) for all \(i\). Let \(\mathcal P\)
be the union of the prime factors of the nonzero determinants
\(a_i b_j-a_j b_i\), together with any specified pool
of cardinality at most \((\log X)^{1+2\tau L}\).

Except for \(O_{B,r,L,\eta}(n^{r-\eta/2})\) queries, all the
following inequalities hold simultaneously:
\[
 \begin{split}
       \log_n{\cal R}_{\mathcal P}(F_i(x)) &> r-\eta
                 &&(1\le i\le m),\\
       \log_n{\cal D}_{q,\mathcal P}(F_i(x)/F_j(x)) &> r-\eta
                 &&(i\ne j,\ q\ge2),\\
       \log_n{\cal D}_{d,\mathcal P}(F_i(x)) &> r-\eta
                 &&(1\le i\le m,\ d\ge 8B/(\eta L)).
 \end{split}                                                   \tag{AF1}
\]
The bounds are uniform in the forms, the prescribed pool and the
query set. In particular, if \(|\mathcal Q|\ge n^r\), they hold
on all but a power-smaller fraction. The simultaneous quantifier
permits choosing \(q\) and the pair separately on each query.

### Counting on a pool by its cardinality

The determinant pool has cardinality
\[
                    O_B(m^2\log X/\log\log X).
\]
Indeed every nonzero determinant has magnitude at most \(2X^{2B}\),
and the product of the first \(k\) primes exceeds \(k!\).
Consequently for all sufficiently large \(X\) the full pool
has cardinality at most \((\log X)^{1+3\tau L}\).
Adding any prime set of squarefree product at most \(X^{2B}\)
preserves this bound up to an absolute factor.

Here is the uniform enumeration used below. For a prime set of
cardinality \(O((\log X)^{1+\alpha})\), and fixed
\(\sigma>\alpha/(1+\alpha)\) with \(\sigma<1\), the sum of
\(p^{-\sigma}\) over that set is at most
\(O_\sigma((\log X)^{(1+\alpha)(1-\sigma)})=o(\log X)\).
Use the sum over the smallest possible integers for this upper
bound. Rankin's inequality and its Euler product then bound the
number of supported positive integers up to \(X^B\) by
\(X^{B\sigma+o_X(1)}\), uniformly in the pool. For positive
rationals of reduced height at most \(X^B\), apply the same bound to
both numerator and denominator, giving
\(X^{2B\sigma+o_X(1)}\). Neither count depends on a height
bound for the product of the pool.
[Supported-integer enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#counting-by-a-sublogarithmic-prime-partition-sum)

For (AF1) take \(\alpha=3\tau L\) and \(\sigma=\eta L/(8B)\).
They satisfy the strict Euler-product margin. If the first
inequality fails, enumerate the outside squarefree radical up to
\(n^{r-\eta}\). It has product at most \(X^B\). For each choice,
the supported integer count bounds the possible values of
\(|F_i(x)|\), each giving at most two queries.

For the second inequality, put \(z=|F_i(x)/F_j(x)|\).
For fixed \(q\), represent its class modulo rational \(q\)-th
powers by the signed valuation remainders of magnitude at most
the original valuation magnitudes. This preserves the height cap
\(X^B\). If the second inequality fails, the representative is
supported on the pool and one outside squarefree product up to
\(n^{r-\eta}\). Any occurring outside product is at most \(X^{2B}\).
The rational enumeration just proved bounds the possible
representatives for each choice by \(X^{2B\sigma+o_X(1)}\).
Each class has \(O_{B,\epsilon}(X^\epsilon)\) queries for any
\(\epsilon>0\) by the fixed affine ratio-class estimate (RP1),
using either sign pattern as needed. The fixed-class part of that
proof uses a polynomial-height determinant and a uniform
fixed-value binary form count. It places no restriction on the pool.

For the last inequality in (AF1), enumerate the \(d\)-free
part of \(|F_i(x)|\) on \(\mathcal P\) and a squarefree outside
support of product up to \(n^{r-\eta}\). The supported part
costs \(X^{B\sigma+o_X(1)}\). The complementary \(d\)-th root
has size at most \(X^{B/d}\le n^{\eta/8}\).
The affine form is determined by its absolute value up to at
most two queries. This unary count fits the same exceptional bound.

There are only \(O_B(1+\log X)\) orders to test. Above
\(B\log X/\log 2\) every nonzero valuation of \(z\) is a defect,
so one order represents the whole remaining range. Taking the
union over orders and pairs costs a height-subpower factor.
Since \(2B\sigma=\eta L/4\), choosing the inner
\(\epsilon<\eta L/8\) and taking the height limit at fixed \(L\)
gives (AF1).
[Uniform fixed affine ratio-class count (RP1)](ABC_RECTANGLE_POWER_CLASS_DEFECT.md#a-uniform-affine-power-class-count)

## Adaptive weights and disjoint prime powers

Outside the determinant pool a prime divides at most one
\(F_i(x)\). To see this, eliminate \(x\) from two affine forms.
Thus for a query satisfying (AF1), the first inequality gives
\[
   \log_n{\cal R}_{\mathcal P}\!\left(\prod_{i=1}^m F_i(x)\right)
                             >m(r-\eta).
                                                               \tag{AF2}
\]
This does not require a polynomial height bound for the displayed
product.

There is a simultaneous weighted version. For any \(q\ge2\)
and any integers \(w_1,\ldots,w_m\), put
\(d_i=q/\gcd(q,w_i)\) with \(\gcd(q,0)=q\).
Make a graph on the indices, joining \(i,j\) if
\(\gcd(d_i,d_j)>1\). Also allow a singleton edge at each index satisfying
\(d_i\ge8B/(\eta L)\). Let \(\nu^*\) be the maximum total
weight of nonnegative edge weights with incident sum at most
one at each vertex. A singleton charges just its own vertex. Then
\[
 \log_n{\cal D}_{q,\mathcal P}
             \left(\prod_{i=1}^m F_i(x)^{w_i}\right)
                              \ \ge\ \nu^*(r-\eta).
                                                               \tag{AF3}
\]
For each edge choose a prime \(\ell\mid\gcd(d_i,d_j)\).
On primes outside \(\mathcal P\), the contribution of \(F_i\)
is tested by \(d_i\mid v_p(F_i)\), irrespective of the sign of
\(w_i\). Its defect mass is therefore at least the mass tested
by \(\ell\), and likewise at \(j\). By (AF1) their sum exceeds
\(r-\eta\), since the ratio defect is the disjoint union
of the two single-form defects. At an allowed singleton, the last inequality of (AF1) gives
the same lower bound with a single form. Weight these inequalities
by any admissible matching and sum. Each form is charged with total
weight at most one, proving (AF3).
For prime \(q\) and at least two active weights, the graph is
complete on the active indices. Its ordinary edges give half
the number of active indices, and at prime orders at least
\(8B/(\eta L)\) the singleton edges give the full number. Equations (AF2)–(AF3) hold for all weights and orders
on the same query set. The weights can depend on the query.

One can subsequently remove all primes of an integer
\(J\le X^H\), even if \(J\) is chosen after observing the
query and the list of forms. The right sides of (AF2)–(AF3)
lose at most \(H/L\), since at most \(\log_n J\) distinct-prime
mass is removed. In (AF1), (AF2) and (AF3), the prime packets
can consist of the full powers in the individual form values.

## Degree of the simultaneous radical extension

Put \(u=r-\eta>0\), and take a query satisfying (AF1).
Fix a number field \(K\), an integer \(q\ge2\), and
\(J\le X^H\) divisible by \(q\,|\operatorname{Disc}(K)|\).
One may choose \(K,q,J\) after observing the query. Then
\[
 \left[K\bigl(\mu_q,F_1(x)^{1/q},\ldots,F_m(x)^{1/q}\bigr)
                :K(\mu_q)\right]
   \ \ge\ q^{\,\max(0,m-C)},\qquad
   C=\max\{1,\lceil 2H/(Lu)\rceil\}.
 \tag{AF5}
\]
To prove it, fix a prime \(\ell\mid q\). An index is lost if
there is no prime outside \(\mathcal P\) and \(J\) with
valuation on its form prime to \(\ell\).
If there are \(k\ge2\) lost indices, every pair of them has
combined \(\ell\)-defect mass exceeding \(u\) by (AF1), and
all that mass is supported on \(J\). Summing over pairs and
dividing by \(k-1\) shows that these indices use more than
\(ku/2\) of its normalized radical mass, since different
forms share no primes outside \(\mathcal P\). This proves
\(k\le C\), including \(k=0,1\).

The primes retained are unramified in \(K(\mu_q)\), since
they divide neither the field discriminant nor \(q\).
Take normalized valuations at one prime above each retained
rational prime and reduce modulo \(q\). Values from different
forms have disjoint coordinate supports, so the order of the
subgroup generated by their valuation vectors is the product
of the individual vector orders. For each \(\ell^b\parallel q\)
at least \(m-C\) vectors have \(\ell\)-primary order \(\ell^b\).
Kummer theory identifies the degree in (AF5) with the order
of the subgroup of \(K(\mu_q)^\times/(K(\mu_q)^\times)^q\)
generated by the signed form values. Its valuation image
gives the stated lower bound.
[Milne, Fields and Galois Theory, Theorem 5.30](https://www.jmilne.org/math/CourseNotes/FT.pdf)

The same conclusion applies to the least \(R\)-full
multipliers of \(|F_i(x)|\) if \(q\mid R\) and
\(2B/(RL)<u\). Replace \(u\) by \(u-2B/(RL)\) in \(C\).
For one form, primes with \(v_p(F_i)\ge R\) have normalized
radical mass at most \(B/(RL)\), directly from the height cap.
At all other primes dividing the form, the least multiplier
has valuation \(R-v_p(F_i)\), congruent to \(-v_p(F_i)\)
modulo \(\ell\) for each \(\ell\mid q\).
Deleting the high valuations therefore leaves the needed
pair inequalities with the reduced margin. The supports
are still disjoint outside \(\mathcal P\).
[Least-multiplier formula](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

## An actual fixed-letter fibre

Consider any set of at least \(n^r\) distinct primitive positive
triples with target \(T\le X\) fixed. Choose
\(m\le(\log X)^{\tau L}\) distinct actual representations
\(a_i+b_i=T\) from it, and take the query \(a+b=T\) from the same
fibre. For this application choose \(m\le(\log X)^{\tau L}\), so two additional forms also fit the bound. Omit the selected representations from the query set, or
retain the \(O(m)\) omitted rows as exceptions. Apply the preceding
result to \(F_i(a)=a-a_i\). In addition to their pair determinants,
put all primes of \(T\) and of the anchor values \(a_i b_i\)
into \(\mathcal P\). This adds \(O(m\log X/\log\log X)\) primes.
A prime outside \(\mathcal P\) dividing \(a-a_i\) cannot divide
\(a\), since it would divide \(a_i\), and cannot divide \(T-a\),
since it would divide \(b_i\). Hence all radical in (AF2) is
fresh relative to the entire query triple and the entire list
of anchors. The defect in (AF3) is fresh in the same sense.

One can include the additional forms \(a\) and \(T-a\) as
tests in (AF1), increasing \(m\) by two. The corresponding
determinants have supports on the anchor terms and \(T\).
Their two-form inequalities give, for instance,
\[
 \log_n{\cal D}_{q,\mathcal P}(a-a_i)
    \ >\ r-\eta
       -\min\{\log_n{\cal D}_{q,\mathcal P}(a),
              \log_n{\cal D}_{q,\mathcal P}(T-a)\}.
 \tag{AF4}
\]
Thus when either actual free letter has a small power defect,
each difference must carry the complementary defect. The proof
also works for a fixed positive summand after changing the sign
of the other affine form. All the original height, primitivity
and radical caps can be imposed on the query set throughout.

In a near-half fibre at the balanced endpoint, \(r\) can approach
\(1/6\). These conclusions rule out a bounded fresh budget for
a growing weighted product with a growing fractional matching
of effective orders. They do not bound the radical of each
difference from above. To turn the amplification into an
unrestricted count saving one still needs to derive a sufficiently
cheap weighted product from the actual additive constraints.
[Growing private packets in fixed-letter stars](ABC_GROWING_PRIVATE_PRIME_PACKETS.md#drawing-the-private-packets)
