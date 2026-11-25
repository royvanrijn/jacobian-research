# Valuation defects of quotient determinants

Cross determinants of a sparse family may have large radicals.
A different restriction concerns the primes whose valuations are
not divisible by a given order. A bounded number of prescribed
anchors converts that restriction into a count on a curve. The
anchor vertices need no repeated factor labels or fixed common
target. The proof uses the actual residual cofactor pairs of the
primitive triples.

## Secants as linear forms on a primitive pair

Fix \(X\ge2\), coprime \(F,G\ge1\), and signs \(\sigma,\tau\).
Consider a family of positive primitive triples with complete
values at most \(X\),
\[
              U_i=F u_i,\qquad V_i=G v_i,\qquad
              T_i=\sigma F u_i+\tau G v_i>0.
\]
Keep one representative per distinct target. The positive pairs
\((u_i,v_i)\) are primitive and have distinct projective directions.
For \(i\ne j\) use the quotient determinant
\[
     H_{ij}=\frac{|u_i v_j-u_jv_i|}
                     {\gcd(u_i,u_j)\gcd(v_i,v_j)\gcd(T_i,T_j)}.
                                                               \tag{SD1}
\]
This is a nonzero integer. The target gcd divides the raw
determinant by taking differences in both complete expressions
for \(T\), using \(\gcd(F,G)=1\). It is coprime to the two
summand gcds by primitivity. All three denominator factors
divide the product \(u_iv_iT_i\) of the anchor data.
[Three orientations of the quotient determinant](ABC_AFFORDABLE_FIBRE_DETERMINANT_COMPLEXITY.md#three-orientations-of-the-same-quotient-determinant)

For a set of primes \(S\), write
\[
          \mathfrak d_q(b;S)
            =\prod_{\substack{p\notin S\\q\nmid v_p(|b|)}}p
               \qquad (b\ne0).
\]
Thus \(\mathfrak d_q\) records only the prime support of the
\(q\)-th-power-free part outside \(S\). Take any fixed \(q\ge2\)
and put \(m=4\) when \(q=2\), and \(m=3\) otherwise.

Choose \(s\) distinct anchors with \(1\le s\le m\), and choose
\(m-s\) of the actual complete query letters \(U,V,T\).
Fix a cap \(B\ge1\) and an extra prime pool \(S_0\) whose product is at most
\(X^C\) for a fixed \(C\). For any \(s\) anchors, the number of
query vertices \(j\) distinct from the anchors satisfying
\[
 \prod_{\substack{p\notin S_0\\
       p\ \text{in the }q\text{-defect support of some }
       H_{ij}\ \text{or selected query letter}}}p\le B
                                                               \tag{SD2}
\]
is \(O_{q,C,\varepsilon}(B X^\varepsilon)\), uniformly in the
anchors and the fixed factors. The same count holds if the
pool in (SD2) is enlarged by the prime supports of the anchors
and of the nonzero pair determinants of the selected linear
forms. Those enlarged pools may depend on the anchors.

To see this, the complete query letters are the forms
\(Fu,Gv,\sigma Fu+\tau Gv\) at the primitive pair \((u,v)\).
For an anchor \(i\) use the form \(L_i(u,v)=u_i v-v_i u\).
These give \(m\) pairwise nonproportional forms. The zero of
an anchor form is its own positive direction, where all three
complete letters are nonzero. Distinct anchors have different
directions. Each form value is nonzero at every eligible query.
Outside the support of \(u_iv_iT_i\) the valuations of
\(|L_i(u,v)|\) and \(H_{ij}\) agree exactly by (SD1).

Put the fixed original coefficient supports, the anchor
supports, and the pair determinants of the \(m\) forms into
the pool. Their product has polynomial height with fixed
\(m,C\). At a prime outside this pool, at most one form value
can have positive valuation at a primitive pair. Otherwise
the two linear congruences at that pair would make their nonzero
coefficient determinant vanish modulo the prime. The outside
defect support of the
product of the form values is consequently exactly the
union of their outside defect supports. Coefficients and
actual values have polynomial height uniformly.

One can choose the joint outside defect support at cost \(B\).
Its allocations and the separate power-free parts on the
fixed pool cost \(X^{o_X(1)}\), by the usual bound on
\(c^{\omega(a)}\) for \(a\) of polynomial height and \(c\)
fixed. The fixed number of possible sign patterns is harmless.
For parity the normalized value of the four-form product
gives a point on a split elliptic double cover. Moving one
rational branch to infinity yields a monic cubic model
with rational \(2\)-torsion and polynomial model and
point heights. The fixed-field bounded-height elliptic
estimate gives \(X^{o_X(1)}\) query directions for each
defect record.

For \(q\ge3\), use a fixed odd prime divisor \(p\mid q\),
or use four when \(q\) is a power of two. The primes of
the resulting lower-order defects off the pool belong to
the same joint record. Three fixed cube-free value
coefficients give a signed primitive diagonal cubic, after
accounting for the common gcd of its roots. That gcd
cubed divides a fixed form determinant, so restoring it
is a divisor cost. Four as an order gives the diagonal
quartic with its bounded-height rational-\(2\)-torsion
elliptic map. For \(p\ge5\), the weights \(1,1,p-2\)
on the three form values give the cyclic curve
\(s(1-s)=h z^p\) after a fractional linear change of
parameter, with \(h\ne0\) of polynomial height. Its rational
point count is subpower by the fixed cyclotomic isogeny
descent and uniform Mordell–Lang. In each case the
projective direction determines the positive primitive
pair, and hence the actual query triple. This proves (SD2).
[Four linear forms and the elliptic completion bound](ABC_MULTI_TARGET_MIXED_KERNEL_COUNTS.md#four-forms-and-elliptic-points-of-bounded-height)
· [The cubic, quartic and cyclic completion bounds](ABC_THREE_TARGET_PRIME_ORDER_STARS.md#the-unweighted-product-for-any-fixed-order-at-least-three)

## Sparse pair bounds and adaptive orders

Take two families from the same fixed divisor record, of
sizes \(M\) on the anchor side and \(N\) on the query side.
Join distinct vertices \(i,j\) when
\(\mathfrak d_q(H_{ij};S_0)\le B_1\).
On the query side require also that the joint outside
\(q\)-defect support of the \(m-s\) selected complete
letters have product at most \(B_0\). Empty selection
means \(B_0=1\). Then the directed incidence count satisfies
\[
       E\ll_{q,C,\varepsilon}
       N+ M N^{1-1/s}(B_0 B_1^s)^{1/s}X^\varepsilon .
                                                               \tag{SD3}
\]
For \(s=1\) the term \(N\) may be dropped.
In fact \(s\) anchors with a common neighbor obey the cap
in (SD2) with \(B=B_0 B_1^s\). If \(E>4sN\), the degree-at-least-\(2s\)
queries contain at least half the edges. Counting ordered
\(s\)-tuples of neighbors at such queries gives
\(\gg_s E^s/N^{s-1}\) configurations by the power mean.
Counting instead by the anchor tuple bounds the same
quantity by \(M^s B_0 B_1^s X^\varepsilon\). The other
case contributes \(O_s(N)\). For \(s=1\) simply apply
the one-anchor bound directly.

In particular, write \(n=X^L\) and suppose one family
has \(N\ge n^r\), with fixed \(r,L>0\). For each fixed
\(q\ge3\), all but a power-smaller fraction of its ordered
pairs have outside \(q\)-defect mass at least
\(r/3-\eta\) in base \(n\), for every \(\eta>0\).
For \(q=2\) the analogous threshold is \(r/4-\eta\).
These use \(M=N\) and \(s=m\) in (SD3), with no power
conditions on any original letter. If the query side
has joint original defect at most \(n^f\) on \(m-s\)
specified letters, the threshold improves to
\((r-f)/s-\eta\) when positive. In the asymmetric statement
this fraction is relative to \(MN\), provided \(M\ge n^r\)
also, so that the first term of (SD3) is negligible.

For the unrestricted-letter statement the two thresholds
hold simultaneously for every possible order \(q\), even
when the order is chosen separately on each pair.
First remove a power-smaller fraction of pairs satisfying
\[
           \operatorname{rad}_{\notin S_0}(H_{ij})
                                         \le n^{r-\eta_0}
\]
for some fixed \(0<\eta_0<r/4\).
For each fixed anchor the three gcds from (SD1) have only
divisor-many choices. The number of \(H_{ij}\le X^2\)
with the indicated radical cap is at most
\(n^{r-\eta_0}X^\varepsilon\), by enumeration of the
part supported on \(S_0\) and then the one-variable
radical count. Orient each pair so that \(u_i\ge u_j\)
at the anchor \(i\). For each resulting raw determinant
the integral query points lie on an affine line with
primitive direction \((u_i,v_i)\), so there is at most
one with \(1\le u_j\le u_i\). Counting each unordered
pair once and restoring both orientations costs at most
a factor two. This fixed-anchor count does not require
dyadic boxes.
[Fixed-anchor divisor count](ABC_AFFORDABLE_FIBRE_DETERMINANT_COMPLEXITY.md#a-fixed-anchor-bound)

If a remaining pair violates the threshold for its
chosen \(q\ge3\), the primes outside \(S_0\) with
positive valuation divisible by \(q\) on \(H_{ij}\)
have product at least \(n^{r-\eta_0-r/3}\).
Since \(H_{ij}\le X^2\), this forces
\(q\le 2/(L(r-\eta_0-r/3))\), a finite bound at the
fixed budget. Apply (SD3) to these finitely many orders
and to parity. All constants may depend on \(C,L,r,\eta\)
before the height limit.

For the estimates using only secants, the extra prime pool
can vary with the anchor, with a fixed polynomial-height
bound. For \(s\) anchors put the union of their pools into
the single fixed pool of (SD2). The common-neighbor bound
is unchanged. For the adaptive-order screening with such
pools, work in a common dyadic box for \(u,v\): the length
of the primitive direction \((u_i,v_i)\) is comparable
coordinatewise to the query caps, so there are \(O(1)\)
query points per signed raw determinant for every fixed
anchor. Thus the fixed-anchor screening works with its
own pool, without changing the orientation of a pair.

These estimates make a small exact-power defect a
sparse relation on actual determinant values. They
require no graph density at the start. They exclude
saturation in an affordable record if a subfamily at
its putative exponent has a density of such pairs
exceeding the bound (SD3). In the larger-defect regime
the argument supplies only a necessary condition
on a surviving family.

## Power defects at bounded scaled orders on fresh primes

Combine (SD3) and its adaptive version with simultaneous
real and projective peeling. Fix a positive lower bound
\(r\) for the fibre exponent and a positive affordable
margin. For any fixed sufficiently small \(\epsilon>0\)
there is a thinning with loss at most \(\epsilon\) in the
normalized exponent, and the following conclusions for all
but a power-smaller fraction of its ordered pairs in late
balanced endpoint boxes. Let \(M\) be the retained
projective modulus, \(h\) the retained angular height
exponent after dividing by \(M\), and \(S_i\) the primes
of \(M U_i V_i T_i\) for the anchor \(i\). There is a fixed
\(\kappa_0\), depending on \(r\) and the fixed margins but
independent of \(\epsilon\), such that with
\(k=\lceil\kappa_0/L\rceil\) the prime powers
\[
          p^{v_p(H_{ij})},\qquad
                p\notin S_i,\qquad 0<v_p(H_{ij})<k
\]
carry at least \((1-\epsilon)h\) numerical weight.
Their primes also avoid the query's complete letters.
Among these bounded-order prime powers,
the normalized radical mass where \(q\nmid v_p(H_{ij})\)
is at least \(r/3-\epsilon\) for every \(q\ge3\)
simultaneously, and at least \(r/4-\epsilon\) for parity.

In fact (PS3) gives a uniform positive lower bound on
\(h\), and the no-cutoff bound in (PS5) holds with
\(\kappa_0\) sufficiently large in terms of \(r\), since
\(h\le2+o_X(1)\). Choose an arbitrarily small fixed loss
\(\gamma\) in the fibre exponent and a small relative
high-order weight tolerance \(\vartheta\). After
peeling, (PS6) bounds the numerical weight at orders
\(\ge k\) outside \(S_i\) by \(\vartheta h\) on almost
all pairs. This part has normalized radical mass at
most \(\vartheta h/(kL)\). The positive gap for \(h\)
allows the angular and support tolerance in (PS2)
to be chosen so that the lower numerical weight at
orders below \(k\) exceeds \((1-\epsilon)h\).
Apply the adaptive-order bound to the retained dyadic
fibre, using the anchor-dependent pools \(S_i\) and
\(R=r-\gamma\) as the exponent lower bound. Removing
the high-order part reduces a \(q\)-defect mass by
at most \(\vartheta h/(kL)\), uniformly in every
choice of \(q\). Taking the fixed \(\gamma,\vartheta\)
and the inner allowances small in terms of \(\epsilon\)
gives both defect thresholds. The exclusion of the
query support follows because a secant prime outside
the fixed factors that divides a query letter also
divides the corresponding anchor letter. The
height limit is always taken before varying the
radical budget.
[Real and projective peeling and the full-packet count](ABC_AFFORDABLE_FIBRE_PROJECTIVE_SECANT_PEELING.md#a-consequence-without-any-secant-power-hypothesis)
