# Power-class defects of affine target rectangles

On two independent affine target columns, a prescribed power class
of the complete target ratio determines the two reduced coefficients
of a fixed nonzero norm equation. This gives a subpower codegree in
each class. A uniform version counts rectangles whose ratio has a
small radical outside the chosen power grid, even if the power order
and the allowed prime pool are chosen separately for each reference
star.

## A uniform affine power-class count

For a positive rational \(z\) and an integer \(q\ge2\), set
\[
 {\cal D}_{q,S}(z)=\prod_{\substack{p\notin S\\q\nmid v_p(z)}}p .
\]
Here \(S\) is a set of primes and \(v_p(z)\) is the signed exponent
in the reduced fraction. The product is \(1\) when there is no such
prime. The height of a reduced rational is the maximum of the
absolute values of its numerator and denominator.

Fix \(K>0\), \(X\ge2\), and two integer affine forms
\[
       T_j(x)=a_j x+b_j,\qquad T_{j'}(x)=a_{j'}x+b_{j'}
\]
whose coefficients have size at most \(X^K\), with \(a_j,a_{j'}\ne0\)
and \(D=a_{j'}b_j-a_jb_{j'}\ne0\). Consider integer \(x\) with
\(1\le T_j(x),T_{j'}(x)\le X^K\). For any \(q\ge2\) and any
prescribed class of \(T_j(x)/T_{j'}(x)\) in
\(\mathbb Q_{>0}^\times/(\mathbb Q_{>0}^\times)^q\), the number
of such \(x\) is \(O_{K,\varepsilon}(X^\varepsilon)\), uniformly
in \(q\), the class and the forms.

Indeed, \(g=\gcd(T_j(x),T_{j'}(x))\) divides \(|D|\), giving
a divisor count for \(g\). The reduced positive numerator and
denominator are coprime, so their exponents modulo \(q\) are
determined by the prescribed class. At a prime with class residue
\(r\in\{1,\dots,q-1\}\), either the numerator has remainder
\(r\), or the denominator has remainder \(q-r\). A prime with
residue zero contributes to neither coefficient. The class leaves
the choice of which term carries each nonzero residue. The two
choices on each nonzero residue give at most
\(2^{\omega(AB)}=X^{o_X(1)}\) assignments, where \(A/B\) is
any realization of the class with height at most \(X^K\).
Fix the assignments. There are then fixed positive coprime
integers \(c,d\le X^K\) with
\[
      T_j(x)=g c u^q,\qquad T_{j'}(x)=g d v^q,\qquad
      a_{j'}c u^q-a_jd v^q=D/g .
\]
The last expression on the right is nonzero. If \(q\) exceeds a
constant depending on \(K,\varepsilon\), the inequality
\(u\le X^{K/q}\) suffices, since \(u\) determines \(x\). For the
remaining finitely many \(q\), use the coefficient-uniform
fixed-value binary form bound proved below. This gives the
claimed estimate after distributing the inner exponent allowance
among the divisor count, residue assignment count and fibre
count.

### The fixed-value bound with moving coefficients

For bounded degree, consider a squarefree homogeneous integer
form \(F(U,V)\) of degree at least two and polynomial coefficient
height. For nonzero \(k\) of polynomial height,
the number of integer pairs of polynomial height with
\(F(U,V)=k\) is subpower, uniformly in \(F,k\), as long as
all degree and height exponents are fixed. In the binomial
application the absence of repeated factors follows from
\(q\ne0\) in characteristic zero and the nonzero leading and
trailing coefficients.

Factor \(F\) into primitive integral irreducible factors, after
removing content. Their coefficient heights are polynomial, by
a bounded integral change removing any zero leading coefficients
and the elementary height bound from the roots and leading
coefficient of the resulting univariate polynomial. If all
factors are linear, fix the values of two nonproportional factors
at divisor cost and solve the resulting linear system. Otherwise
fix the nonzero value of an irreducible factor \(P\) of degree
\(d\ge2\). With \(a=P(1,0)\ne0\) and a root \(\theta\) of
\(P(T,1)\), the element \(a(U-\theta V)\) is integral in
\(\mathbb Q(\theta)\) and has norm \(a^{d-1}P(U,V)\).
An ideal of this norm can be assigned to at most
\(\tau(|a^{d-1}P(U,V)|)^d\) possibilities. In each ideal the
generators differ by units. All embedding moduli of the
generators are polynomially bounded above, and inverse
polynomially bounded below by the nonzero integral norm.
Two distinct unit log vectors in degree at most \(d\) have a
uniform positive minimum distance. Indeed units whose embedding
moduli lie in \([e^{-1},e]\) have bounded integral minimal
polynomial coefficients and hence belong to a fixed finite
set. This also bounds the size of the log map's kernel.
The box of side \(O(\log X)\) therefore has only
\(O((1+\log X)^d)\) possible generators per ideal. The generator
determines \(U,V\) because \(\theta\) is irrational. This is the
uniform fixed-value argument used for primitive norm factors.
[Uniform weighted-power factors](ABC_MINIMAL_REPAIR_CYCLOTOMIC_SPLIT.md#direct-target-bounds-for-reducible-weighted-powers)
· [Linear and nonlinear factors on a primitive carrier](ABC_RATIONAL_PROJECTIVE_CARRIERS.md#linear-and-nonlinear-factor-counts-on-a-primitive-carrier)

For \(q=2\) the binomial fibre also follows from the
coefficient-uniform diagonal quadratic bound of Vaughan and
Wooley, restated in Browning and Verzobio's powerful-triple
paper. The argument above uses the same fixed nonzero
norm mechanism at every bounded order, then the height
bound at unbounded orders.
[Browning and Verzobio, Lemma 2.3](https://arxiv.org/html/2608.24512v1#S2.SS1)

## A bound uniform over defects and power orders

Keep the two affine forms fixed. Prescribe a positive rational
\(\chi\) of height at most \(X^K\) and a prime pool \(S\) with
prime product at most \(X^K\). For \(Z\ge1\), the count of
eligible rows \(x\) for which
\[
              {\cal D}_{q,S}\bigl(\chi T_j(x)/T_{j'}(x)\bigr)
                    \le Z\quad\text{for some integer }q\ge2
\]
is
\[
                         O_{K,\varepsilon}(Z X^\varepsilon).
                                                               \tag{RP1}
\]
The implied constant is independent of the pool, rational anchor
and affine forms. The same bound applies to any fixed subfamily.

To see this, put \(z=\chi T_j(x)/T_{j'}(x)\). Its reduced height
is at most \(X^{2K}\), and every nonzero valuation has absolute
value at most \(2K\log_2 X\). Above this bound the defect is
independent of \(q\) and equals the radical of \(z\) outside
\(S\), so one representative order suffices for the entire range.
For any remaining \(q\), enumerate the outside defect support
by its squarefree product at most \(Z\). For a fixed such
support, represent the class of \(z\) by reducing each valuation
modulo \(q\) to the remainder of the same sign and no larger
absolute value. The representative has height at most \(X^{2K}\)
and lies on that support together with \(S\). Any occurring
outside defect has product at most \(X^{4K}\). The fixed-support
count at polynomial height gives \(X^{o_X(1)}\) possible class
representatives, uniformly in \(q\) and the support. This use
of signed remainders preserves the height when a negative
exponent has very small absolute value and \(q\) is large.
Multiplying the resulting class by \(\chi^{-1}\) fixes the
class of \(T_j/T_{j'}\). Apply the affine count, and absorb the
\(O_K(1+\log X)\) orders into the exponent allowance.
[Uniform fixed-support count](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

### Several simultaneous affine targets

There is a fixed-class extension to a given list of \(s\)
pairwise nonproportional affine forms as above, with \(s\) fixed.
Prescribe a class of \(\prod_{i=1}^{s}T_i(x)^{e_i}\) modulo
rational \(q\)-th powers for integral \(e_i\). If
\[
 \gcd\!\left(\frac{q}{\gcd(e_i,q)},\frac{q}{\gcd(e_j,q)}\right)>1
\]
for at least one pair \(i\ne j\), the count is
\(O_{K,s,\varepsilon}(X^\varepsilon)\), uniformly in the
weights, class, forms and \(q\).

Outside the prime pool of all pair determinants, each prime
occurs on at most one of the target values. On this part of
target \(i\), its valuation \(v\) satisfies \(e_i v\equiv r\bmod q\),
with \(r\) the prescribed product class residue. Whenever
soluble, the congruence specifies \(v\) modulo \(q/\gcd(e_i,q)\).
Choose a prime \(\ell\) dividing both displayed moduli. For
\(T_i/T_j\) this fixes the residue at every prime outside the
pool up to three choices at each nonzero product class prime,
according as it belongs to \(i\), \(j\), or neither. The number of
such choices is subpower. Indeed the nonzero class support lies
in the support of the product of the \(s\) polynomial-height
targets of any actual realization. Any prime where the product
class has zero residue has zero residue modulo \(\ell\) in
both selected targets outside the pool. Fix the full parts of
the two targets on the pool by fixed-support enumeration at
subpower cost. The \(\ell\)-class of the ratio is now fixed, and
the two-column affine bound applies. In particular every
prime-order class relation involving two or more active columns
has subpower common rows.

## Multiplicative rectangles in the arithmetic graph

Let the full target on an edge of a bipartite graph be
\[
              T_{xj}=Fxy_j+\sigma Gk_j>0 ,
             \qquad \sigma\in\{-1,1\},\quad F,G>0 ,
\]
with polynomial-height labels and fixed \(F,G\). Suppose the
row labels \(x\) are distinct and the column labels
\(j=(y_j,k_j)\) are distinct positive primitive pairs. Then
the affine determinant is nonzero for every two columns.
For an ordered rectangle with distinct rows and columns
write
\[
             Q(x,x';j,j')=
                \frac{T_{xj}T_{x'j'}}{T_{xj'}T_{x'j}} .
\]
For each reference star \((x,j,j')\) choose a prime pool
\(S_{x,j,j'}\) with polynomial prime product. With \(d(x)\)
the degree of \(x\), the number of ordered rectangles
satisfying \({\cal D}_{q,S_{x,j,j'}}(Q)\le Z\) for some
\(q\ge2\) is
\[
          O_\varepsilon\!\left(
                 Z X^\varepsilon\sum_x d(x)(d(x)-1)\right).
                                                               \tag{RP2}
\]
Fix the reference star and use (RP1) on \(x'\), with
\(\chi=T_{xj'}/T_{xj}\) and invert \(Q\), which leaves the
defect the same. Then sum over the stars. This conclusion
holds without a compressed-core hypothesis or a bound on
the ordinary radical of a projective secant.

For comparison with the complete-label endpoint, put \(n=X^L\)
at fixed \(L\), take both part sizes at most \(n^{c+o(1)}\),
and suppose there are \(n^{2c-o(1)}\) edges. The total
number of ordered rectangles is at least \(n^{4c-o(1)}\),
by two applications of Cauchy–Schwarz, first to the row
degrees and then to the ordered column-pair codegrees.
There are at most \(n^{3c+o(1)}\) reference stars. For
\(Z=n^\beta\) and \(\beta<c\) separated by a fixed positive
margin, (RP2) shows that the rectangles satisfying the
adaptive low-defect condition form a power-smaller set.
The inner exponent allowance is first taken sufficiently
small at fixed \(L\). At the complete-label minimum-mass
endpoint \(c=1/6\). The same proof applies to unequal part
caps and without taking any small-budget limit.
[Complete-label arithmetic graph](ABC_DIFFUSE_FACTOR_GRAPH.md#the-complete-label-reduction-close-to-minimum-mass)

More generally write \(A,B,M\) for the two part sizes and
the edge count. If at least a fraction \(\vartheta>0\) of
all ordered rectangles is low-defect in (RP2), then
\[
 M\ll_\varepsilon A+
       B\sqrt{A(1+Z X^\varepsilon/\vartheta)}.                 \tag{RP3}
\]
Indeed, let \(S_2=\sum_x d(x)(d(x)-1)\). The total rectangle
count is at least \(S_2^2/[B(B-1)]-S_2\) when \(B\ge2\).
Apply (RP2) and use \(M^2/A-M\le S_2\). If \(B\le1\), the
bound holds from \(M\le A\). This gives an edge saving even
when the low-defect share of rectangles is allowed a small
power loss and the power orders differ across rectangles.

### An edge bound for adaptive power defects

Prescribe a positive rational \(C_j\) and a prime pool \(S_j\)
for each column, and allow a positive rational \(D_x\) depending
only on the row. All rational heights and prime products are
polynomial in \(X\). Suppose every edge has an integer \(q\ge2\)
for which
\({\cal D}_{q,S_j}(T_{xj}/(C_jD_x))\le Z\). The order may vary
with the column or with each edge. Then
\[
                  M\ll_\varepsilon X^\varepsilon
                          \bigl(A+B\sqrt A\,Z\bigr).            \tag{RP4}
\]
To prove this, colour each edge by one admissible order, using
one common order to represent all orders above a fixed multiple
of \(\log X\), as in (RP1). In one colour, a common row of two
columns has
\[
 {\cal D}_{q,S_j\cup S_{j'}}\!\left(
       \frac{C_{j'}T_{xj}}{C_jT_{xj'}}\right)\le Z^2 ,
\]
since the row factor cancels. Consequently its colour codegree is
\(O(Z^2 X^\varepsilon)\) by (RP1). Cauchy–Schwarz on the coloured
row degrees and the sum of the column codegrees bounds the edge
count of this colour by \(O(X^\varepsilon(A+B\sqrt A\,Z))\).
Sum over \(O(1+\log X)\) colours, using a smaller inner exponent
allowance for each application. If one order works throughout
each row, the same proof applies by (RP1) directly to each
column pair, without colours or a height restriction on \(D_x\).

For example, if each column consists of shifted targets
\(T_{xj}=C_j u_{xj}^{q_j}\) with a fixed column coefficient
and arbitrary column orders \(q_j\ge2\), (RP4) gives
\(M\ll_\varepsilon X^\varepsilon(A+B\sqrt A)\).
At the complete-label endpoint this lowers the edge exponent
from \(2c\) to \(3c/2\). More generally a normalized defect cap
\(Z=n^\beta\) on each edge gives a strict saving if
\(\beta<c/2\) by a fixed margin. This provides a transfer from
shifted exact-power graphs even when the exact exponent varies
between columns, and keeps a quantitative allowance for the
off-grid part.

### The original supports and cofactor ratios

In the primitive signed triples
\[
        (Fxy_j,\ Gk_j,\ T_{xj}),
\]
the reference pool in (RP2) can contain all primes on all
three letters of the two reference edges. It then contains
the primes of \(F,G,y_j,y_{j'},k_j,k_{j'}\). Outside that pool,
a prime on either of the two non-target letters of either
query edge divides \(x'\). It cannot divide either query
target by primitivity. Thus the defect in (RP2), for this
pool, is unchanged if one also removes the primes on the
non-target letters of the two query edges. This allows all
three reference letters and both non-target query letters to
be discounted without paying for a new dictionary.

For a factorization \(T_{xj}=W_{xj}R_{xj}\), (RP2) concerns
the full targets. The defect of the full ratio is at most
that of the cofactor ratio times the defect of the core ratio,
by submultiplicativity of \({\cal D}_{q,S}\). In particular the
same count applies to cofactor ratios on rectangles on which
the core ratio is a rational \(q\)-th power. More generally,
both defect costs may be charged in \(Z\). When the reference
pool contains the reference targets it is enough to charge
the two query core defects, since the reference core factors
are supported on the pool. This is an arithmetic complement
to the compressed-core weighted bound, which uses physical
height and a fixed cofactor ratio instead.
[Compressed-core weighted codegrees](ABC_COMPRESSED_CORE_WEIGHTED_CODEGREES.md)

### Least target repairs in a near-full graph

At a fixed integer order \(r\ge2\) write the least target
repair and the repaired letter as
\[
       \mu_r(T)=\prod_{0<v_p(T)<r}p^{r-v_p(T)},
       \qquad H_r(T)=\mu_r(T)T .
\]
This is the actual least common multiplier on the triple when
the requirement is to make the designated target \(r\)-full.
A prime can contribute to \({\cal D}_{r,\varnothing}(H_r(T))\)
only if its original target valuation exceeded \(r\). For
\(T\le X^K\) it follows that
\[
       {\cal D}_{r,\varnothing}(H_r(T))\le X^{K/r}.
\]
Let \(Q_\mu\) be the cross ratio of the four \(\mu_r(T_{xj})\)
in the same orientation as \(Q\). Take the reference prime
pool to include the two reference targets. The corresponding
cross ratio \(Q_H\) of the repaired letters has
\({\cal D}_{r,S}(Q_H)\le X^{2K/r}\): only the two query
letters can contribute outside the pool. Since
\(Q_H=Q Q_\mu\), (RP2) therefore gives
\[
  \#\bigl\{\text{ordered rectangles}:
           {\cal D}_{r,S}(Q_\mu)\le Z\bigr\}
      \ll_{\varepsilon}
          ZX^{2K/r+\varepsilon}\sum_x d(x)(d(x)-1).            \tag{RP5}
\]
In a dense graph with part exponent \(c\) in base \(n=X^L\),
and at orders with \(rL>2K/c\) by a fixed margin, the repairs
on almost all rectangles must have defect exceeding
\(n^{c-2K/(rL)-\eta}\) for any fixed positive \(\eta\) smaller
than the exponent on the right before subtracting \(\eta\).
The same conclusion holds with a larger reference pool under
the polynomial height requirement. Thus a repair model whose
rectangle cross ratios are all rational \(r\)-th powers forces
the edge bound \(O_\varepsilon(A+B\sqrt A\,X^{K/r+\varepsilon})\)
by (RP3). In particular multiplicative separation of the
repairs modulo \(r\)-th powers forces this bound, including
when the separating row and column factors are arbitrary
common rational multipliers. No height bound on those factors
enters (RP5).

If the three original letters are primitive and the requirement
is simultaneous \(r\)-fullness, the least common multiplier
inserts \(p^r\) at each prime occurring at an original order
between \(1\) and \(r-1\). At such a prime, a positive
multiplier exponent introduces the prime on the other two
letters and must be at least \(r\). The multiplier has
trivial \(r\)-power class. The arithmetic obstruction in (RP5)
is a statement about the least repair of the designated
target, as used in the primitive quotient count.
[Least target repair and true support occurrence](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md)

An independent affine scale is essential to the two-column
power-class count. For homogeneous forms in two varying
primitive coordinates, the forms \(T_1(U,V)=U\) and
\(T_2(U,V)=V\) have \(\gg_q X^{2/q}\) positive coprime pairs
of size at most \(X\) with \(T_1/T_2\) a rational \(q\)-th
power, obtained by taking coprime \(q\)-th power bases below
\(X^{1/q}\). The nonzero fixed term after elimination is
precisely what supplies the norm bound in the affine case.

The reference supports can be discounted at any radical caps.
It remains to show that a large fraction of the remaining query
target ratios has small defect under the three original caps.
Weak limits alone cannot supply this step. For any fixed \(q\),
equal masses on the degrees in \(q\mathbb Z\) with \(Le\in[2,4]\)
and on the corresponding degrees in \(q\mathbb Z+1\) both tend
to the same uniform measure, while their defect masses are zero
and the full mass respectively. The bound therefore supplies
a discrete arithmetic restriction on near-full graphs,
including at nonminimum nonatomic mass, while the step
forcing low defects or another power saving from the
unrestricted three-letter budgets remains open.
