# Adaptive power defects in stars of actual targets

Prescribed relative odds make the complete targets in a common-column
star into binary linear forms in a primitive pair. A small defect in
a weighted product of those targets then bounds the number of stars.
The order can vary at each height, and the bound does not require
an affine realization by complete labels. A separate application
counts stars of least target repairs whose classes separate into
row and column factors.

## Uniform counts for weighted target products

Take \(m\) ordered primitive positive triples with common signs
\[
        T_i=\sigma U_i+\tau V_i>0,\qquad
        U_i,V_i,T_i\le X.
\]
Prescribe the relative odds
\[
          r_i=\frac{U_iV_1}{V_iU_1},\qquad r_1=1,
\]
with all \(r_i\) distinct. Let \(S\) be a fixed prime pool with
\(\prod_{p\in S}p\le X^A\) for some fixed \(A\). For positive rational
\(f\), write
\[
       D_{q,S}(f)=\prod_{\substack{p\notin S\\q\nmid v_p(f)}}p .
\]
The valuation is the signed valuation of the reduced rational. For
\(B\ge1\), there are the following uniform bounds.

- For two targets, the count with
  \(D_{q,S}(T_1/T_2)\le B\) at any fixed \(q\ge2\) is
  \(O_{A,\varepsilon}(B X^{2/q+\varepsilon})\). \(\tag{AH1}\)
- For three targets, the count with
  \(D_{q,S}(T_1T_2/T_3^2)\le B\) for some \(q\ge3\) is
  \(O_{A,\varepsilon}(B X^\varepsilon)\). \(\tag{AH2}\)
- For four targets, the same subpower bound holds for
  \(D_{q,S}(T_1T_2/(T_3T_4))\le B\) at some \(q\ge2\).
  \(\tag{AH3}\)

In particular, the orders allowed in (AH2) need not be prescribed
or bounded. The two-target estimate is useful when the order is
large in terms of the radical budget. The three-target assertion
cannot include order two with the displayed weights.

### Reduction to a primitive pair

Write \(x=\gcd_i U_i,\ y=\gcd_i V_i\) and
\(U_i=d_i x,\ V_i=e_i y\). The pairs \((U_i,V_i)\) are primitive,
so \((x,y)=1\). At a prime outside the numerators and denominators
of the \(r_i\), the signed valuations of \(U_i/V_i\) are all the
same. The two original nonnegative valuations have minimum zero
and are determined by their difference. Thus all the quotients
\(d_i,e_i\) are supported on the fixed ratio pool, with only
\(X^{o_X(1)}\) possibilities for their complete values. Fix them
and use
\[
                    L_i(x,y)=\sigma d_i x+\tau e_i y .
\]
These forms are pairwise nonproportional by the distinct odds
hypothesis. Add the primes of their pair determinants to the
pool. Its prime product still has polynomial height, uniformly
in the fixed data. Outside this enlarged pool no prime can divide
two \(L_i(x,y)\), by primitivity of \((x,y)\). The reconstruction
from \(x,y\) to the ordered witness triples is exact for the fixed
coefficient record.
[Pair determinants and coefficient records](ABC_MULTI_TARGET_MIXED_KERNEL_COUNTS.md#fixed-cross-ratios-and-the-common-factors)

For the two-target assertion, outside the enlarged pool each
prime appears in at most one of \(L_1,L_2\). Their \(q\)-free
coefficient parts therefore share a total support radical at most
\(B\) there. Enumerate that radical and then the two coefficient
parts at cost \(B X^{o_X(1)}\). The latter are positive integers
of height at most \(X\), supported on the prescribed pool and
the enumerated radical. The fixed-support count at polynomial
height is uniform in \(q\), including when the possible residue
exponents are large. Writing \(L_i=D_i u_i^q\), there are at most
\(X^{2/q}\) choices for the two positive integral roots. Two
distinct forms determine \((x,y)\) from their values. This proves
(AH1).
[Uniform fixed-support count](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-radical-count-contradiction-to-the-claimed-smooth-number-lower-bound)

### Three targets at bounded orders

If \(q\) has an odd prime divisor \(p\), reduction to the \(p\)-class
of the same rational product only decreases the support cap. The
class equals that of \(T_1T_2T_3^{p-2}\), since the quotient is a
\(p\)-th power. For fixed \(p\), use the three-target prime-order
bound (PT1), including the separately proved cubic case. The
constant can depend on \(p\).
[Three actual targets and cubic variant](ABC_THREE_TARGET_PRIME_ORDER_STARS.md#count-for-three-actual-targets-at-prescribed-odds-quotients)
· [The cubic variant](ABC_THREE_TARGET_PRIME_ORDER_STARS.md#the-cubic-variant-from-the-linear-relation-of-the-targets)

The remaining bounded orders \(q\ge3\) are divisible by four.
Reduce to \(q=4\) and fix the class representative
\(\Delta\in\mathbb Q_{>0}\) of \(T_1T_2/T_3^2\), taken with signed
residue exponents of absolute value at most three. Its outside
support has radical at most \(B\), and its height is polynomial
in \(X\). There are \(B X^{o_X(1)}\) possible representatives
by fixed-support enumeration. The equation is now
\[
               L_1 L_2/L_3^2=\Delta z^4,\qquad z\in\mathbb Q_{>0}.
\]
Express \(L_3=aL_1+bL_2\), where \(a,b\in\mathbb Q^\times\) have
polynomial height, and set \(s=aL_1/L_3\). Thus
\(s(1-s)=ab\Delta z^4\). Set \(u=1-2s\) and \(h=-4ab\Delta\ne0\).
Then \(u^2=1+hz^4\). The map
\[
                    (\xi,\eta)=(h z^2,h z u)
\]
takes these points to \(E:\eta^2=\xi^3+h\xi\), with rational
two-torsion \((0,0)\). The coefficients and the displayed point
coordinates have polynomial height in \(X\). Since \(z\ne0\),
the map has bounded fibres, and \(s\) determines the primitive
pair \(x,y\). Dujella's uniform bounded-height estimate gives
\(X^{o_X(1)}\) completions for each record. This proves the
three-target estimate at order four and hence at every fixed
multiple of four.
[Dujella, Theorem 1.1 and Corollary 1.3](https://arxiv.org/html/2312.03655v1)

### Large adaptive orders and the parity exception

For every \(q\), outside the determinant pool a prime appearing
on \(L_1\) or \(L_2\) has exponent in the weighted product equal
to its own exponent there. Thus its \(q\)-free contribution on
these two forms is detected without any information about the
third form. The proof of (AH1) gives \(O(BX^{2/q+\varepsilon})\)
uniformly, as in that proof. Choose \(q\ge Q_\varepsilon\) with
\(2/Q_\varepsilon<\varepsilon/4\) and take the inner error
allowance smaller than the remaining margin. For \(3\le q<
Q_\varepsilon\) the finite list of bounded-order estimates just
proved is uniform after enlarging the constant. The valuations
of the complete rational product have absolute value at most
\(4\log_2 X\), so above that order all defects have stabilized.
Summing over \(O(1+\log X)\) orders, with smaller initial exponent
allowances, proves (AH2). At order two the exponents of
\(T_1T_2/(T_3T_4)\) agree modulo two with those of
\(T_1T_2T_3T_4\). The four-target mixed-kernel count therefore
gives (AH3) at parity. At \(q\ge3\), outside the full
pair-determinant pool each prime occurs on at most one form,
so the cap for the four-target quotient implies that for
\(T_1T_2/T_3^2\) on this complement. Apply (AH2) to the
first three forms, retaining the coefficient record for the
fourth. This completes (AH3) uniformly over the orders.
[Four actual targets and the mixed kernel](ABC_MULTI_TARGET_MIXED_KERNEL_COUNTS.md#the-mixed-kernel-count)

For a counterexample at order two with three targets, take
\(x=v^2-u^2,\ y=2u^2-v^2\), with coprime positive \(u,v\equiv1\bmod6\)
in fixed disjoint intervals of length \(\asymp X^{1/2}\) where
\(u<v<\sqrt2\,u\). Scale the intervals by a small enough absolute
constant that all \(i x+y,\ i x,\ y\) for \(i=1,2,3\) are
at most \(X\). The change from \((u^2,v^2)\) to \((x,y)\)
is unimodular, so \(x,y\) are primitive. Also \((y,6)=1\), hence
\((i x,y)=1\) for \(1\le i\le3\). The triples \((i x,y,i x+y)\)
have fixed distinct relative odds \(1,2,3\), and their first two
targets are \(u^2,v^2\). Thus \(T_1T_2/T_3^2\) is a square
for every pair, with at least a constant times \(X\) distinct
tuples by the elementary coprime-pair density in the chosen
residue classes. That density follows by Möbius inversion: the
sum on common divisors coprime to six has positive main
coefficient \(\prod_{p>3}(1-p^{-2})\), and the endpoint errors
are \(O(X^{1/2}\log X)\) on these fixed intervals. This is a
height example for the uniform star count and does not assert
a small radical budget on all the letters.

## An adaptive graph consequence

Use a simple bipartite graph of actual witness triples, with
\(R\) rows, \(C\) columns and \(E\) edges. Assume the full triple
determines its edge labels. For each ordered triple of distinct
rows, suppose the three relative odds arising at common columns
have at most \(M\) possible vectors. In a common column the
directions are distinct, since equal primitive odds determine
the entire witness. Suppose the condition in (AH2) holds for
each such star, with the prime pool fixed after the row tuple
and odds vector have been chosen. Counting ordered stars gives
\[
                 E\ll C+R C^{2/3}(M B X^\varepsilon)^{1/3}.
                                                               \tag{AH4}
\]
Indeed, if \(E\ge6C\), the columns of degree at least three
carry enough edges to give \(\gg E^3/C^2\) ordered stars by
the power mean, whereas (AH2) gives
\(O(R^3 M B X^\varepsilon)\) stars. For (AH1), use ordered
pairs and obtain
\[
                 E\ll C+R C^{1/2}(M B X^{2/q+\varepsilon})^{1/2},
                                                               \tag{AH5}
\]
where \(q\) is fixed for the count. For four rows use (AH3)
and replace \(1/3,2/3\) in the second term of (AH4) by
\(1/4,3/4\) respectively.

For example, suppose each edge has
\(D_{q,S}(T)\le n^{f}\), where \(n=X^L\) and the order \(q\ge3\)
can vary by edge. This assumption gives (AH4) with \(B=n^{3f}\)
up to an inner subpower factor. Color edges by one eligible
order, with all \(q>\log_2 X\) represented by a single such
order, and apply (AH4) in each of the \(O(1+\log X)\) colors.
For a mixed product, a prime can contribute to the defect
only if it contributes on an edge in that star. The odds
bounds and all dyadic part caps are inherited by the subgraphs.
The same edge conclusion for graphs that allow \(q=2\) must
keep the separate four-row estimate on the parity subgraph.

## Separable least target repairs without affine realization

For fixed order \(r\ge2\), put
\[
          \mu_r(T)=\prod_{0<v_p(T)<r} p^{r-v_p(T)},\qquad
          H_r(T)=T\mu_r(T).
\]
The off-grid prime support of \(H_r(T)\) is contained in
the support of the original prime powers with orders at
least \(r\). Consequently \(D_{r,\varnothing}(H_r(T))\le X^{1/r}\)
when \(T\le X\). Suppose the edge classes of \(\mu_r(T)\) modulo
rational \(r\)-th powers separate into row and column factors.
No height restriction on the separating factors is needed.

For a fixed ordered pair or triple of rows with common columns,
the product of the repair factors with respective weights
\((1,-1)\) or \((1,1,-2)\) has the same class at every common
column, since these weights sum to zero. Choose one reference
common column, and put all its target prime supports in the
pool \(S\). Its prime product is at most \(X^m\). At every other
common column the weighted repair class contributes no defect
outside \(S\). The identity \(H_r(T)=T\mu_r(T)\) and the
preceding upper bound for the repaired defects show that the
weighted product of the original targets has outside defect
at most \(X^{m/r}\), for \(m=2,3\). The choice of the pool can
be made for each row tuple before fixing the odds vector.
Applying (AH5) and, for \(r\ge3\), (AH4), gives respectively
\[
 \begin{split}
       E&\ll C+R C^{1/2} M^{1/2}X^{2/r+\varepsilon},\\
       E&\ll C+R C^{2/3} M^{1/3}X^{1/r+\varepsilon}.
 \end{split}                                                   \tag{AH6}
\]
Here \(M\) bounds the odds vector count for the corresponding
star size. More generally, if the weighted repair class
takes at most \(J\) different values on each such row tuple,
choose a separate reference column for each class. The same
argument applies with \(MJ\) in place of \(M\).

In the actual outside projection with nonatomic masses
\(c_U,c_V,c_T\), suppose the target letter supplies a cheap
full block, for instance if \(c_T>0\). The usual outside
pair-fibre and odds-quotient arguments give simple graphs
with part exponents at most \(c_U,c_V\), with \(M\) of
vanishing normalized exponent. Freezing the common block
costs as little as desired. The least-repair hypothesis is
inherited by the retained edges. Therefore, if \(r\ge3\)
is fixed at the height limit and may vary with \(L\), the
normalized eligible count exponent satisfies
\[
 d\le \max\left(c_V,\ c_U+
       \min\left\{\frac12 c_V+\limsup\frac{2}{rL},
                   \frac23 c_V+\limsup\frac1{rL}\right\}\right).
                                                               \tag{AH7}
\]
The roles of \(U,V\) can also be reversed. In the symmetric
minimum nonatomic sector, the three-row term is strictly less
than \(1/3\) when \(\liminf rL>18\). The two-row term tends
to \(1/4\) when \(rL\to\infty\). This uses odds at actual
stars and does not pay the exponent cost of an affine or
complete-label conversion. The hypotheses concern least
repair of a designated target in an original primitive
triple; they do not apply to trivial classes of multipliers
which make all three letters full simultaneously.
[Outside odds bounds and anchor construction](ABC_DIFFUSE_FACTOR_GRAPH.md#anchor-sets-for-arbitrary-nonatomic-masses)
· [Affine least-target repair comparison](ABC_RECTANGLE_POWER_CLASS_DEFECT.md#least-target-repairs-in-a-near-full-graph)

The unrestricted count still requires a restriction on the
repair class dictionaries or a different use of the original
additive relation. A small radical budget alone has not been
used to force the separation hypothesis.
