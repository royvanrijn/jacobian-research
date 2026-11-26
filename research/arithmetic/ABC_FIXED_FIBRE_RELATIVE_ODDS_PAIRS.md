# Close coefficient pairs in genuine fixed-letter fibres

Cross-ratio recovery turns a short relative-height condition
into a pointwise pair count. The fixed original letter
provides either a positive-depth full modulus or a short
real window for the free ratio. The real proof uses the
product of the two rational height exponents and handles
orders that grow arbitrarily slowly with the budget.

## A product-of-heights count in a real monomial strip

**FRP1.** For a reduced rational number \(u\), let \(H(u)\) be its multiplicative
height. Fix \(b,\tau>0\), an integer \(k\ge1\), and \(R,E\ge0\).
There is an absolute \(c_*>0\) with the following property. Suppose
\[
       \max(R,E)\le \tau/10,\qquad R E\le c_* b\tau .
\]
For every \(K>0\), the number of distinct outputs \(r^k z\) from positive
reduced rational pairs with
\[
       H(r)\le X^R,\qquad H(z)\le X^E,\qquad
       |r^kz/K-1|\le C X^{-\tau}
\]
is \(O(X^{b+\varepsilon})\) for each fixed \(\varepsilon>0\),
uniformly in \(K\). The implied constant can depend on the displayed
parameters and \(C,\varepsilon\). In particular, \(k,R,E,b\) may all
depend on a fixed budget \(L\) before taking the height limit. The
numerical conditions have absolute constants independent of \(k\).
This is the two-coordinate Taylor argument adapted from
[the unequal-height strip count (RUW1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights).
For method context, Heath-Brown's unequal-height approximate
determinant uses a related product of logarithmic height
scales. The proof here completes along an explicit monomial
curve by a one-variable zero estimate.
[Heath-Brown, Square-free values of \(n^2+1\), §3 and Lemma 1](https://www.impan.pl/shop/en/publication/transaction/download/product/82579)

Normalize each variable by the lower endpoint of its dyadic box.
This raises its reduced rational height by at most a factor two. The
new coordinates \(x,y\) lie in \([1,2]\). The normalized target
\(K'\) varies over a fixed positive compact interval, depending
only on \(k\), in any nonempty box. The strip condition gives
\[
                    y=f(x)+O_k(X^{-\tau}),\qquad
                    f(t)=K't^{-k}.
\]
Derivatives of any fixed order of \(f\) are bounded on \([1/2,3]\),
uniformly in \(K'\). There are \(O((1+\log X)^2)\) boxes.
If \(R\le b/10\) or \(E\le b/10\), enumerate the coordinate of
the smaller height at cost \(O(X^{b/5})\). Once that coordinate
is fixed, the other lies in an interval of length \(O_k(X^{-\tau})\)
by the strip and the lower derivative bound for \(f\) on
\([1,2]\). This is shorter than the spacing of two distinct
rationals of the allowed height.

Otherwise put
\[
  D_x=\lceil A E/b\rceil,\quad D_y=\lceil A R/b\rceil,\quad
  N=(D_x+1)(D_y+1),\quad m=\lfloor N/2\rfloor,\quad S=D_xR+D_yE
\]
with \(A\) a sufficiently large absolute integer. Since
\(R,E>b/10\), one has \(S\le(2A+20)RE/b\) and
\(A^2RE/b^2\le N\le(A+20)^2RE/b^2\).
Choose \(A\) large enough that \(b(m-1)/2>2S\).
Then choose \(c_*\) small enough that \(b(m-1)<\tau\) and
\(8S<\tau\).

Partition \([1,2]\) into intervals of length at most \(X^{-b}\).
In such an interval take any \(m\)-row minor of the evaluations of
the rectangular monomials in \(x,y\) of bidegrees at most
\((D_x,D_y)\). Taylor-expand the reference evaluations
\(x^i f(x)^j\) through order \(m-2\) around the interval midpoint.
Each evaluation row is a combination of the same Taylor vectors
with coefficients \(O(X^{-b\ell})\) in order \(\ell\).
The remainder row, including the \(y-f(x)\) error, is bounded by
\(O(X^{-b(m-1)})\). Terms in the determinant using the same
Taylor vector twice vanish. The remaining determinant is
\(O(X^{-bm(m-1)/2})\). Multiplying each evaluation row by the
product of its two reduced denominators to powers \(D_x,D_y\)
makes the minor integral, at multiplication cost \(O(X^{mS})\).
All these integral minors therefore vanish for sufficiently
large \(X\).

The rank on each interval is below \(m\). There is a nonzero
integer polynomial \(F\) of the stated bidegree and coefficient
height \(O(X^{2S})\) vanishing on the samples of the interval.
Indeed, evaluate coefficient vectors from an integer cube of
side \(\lceil X^{2S}\rceil\) on a maximal independent set of
cleared rows. Each row gives \(O(X^{3S})\) possible values,
whereas there are \(\gg X^{2SN}\) coefficient vectors and fewer
than \(N/2\) independent rows. Subtract two coefficient vectors
with identical evaluations.

Factor \(F\) over \(\mathbb Q\) and take primitive integer
irreducible factors \(G\). Each has bidegree at most
\((D_x,D_y)\) and coefficient height \(O(X^{2S})\),
with constants at fixed degrees. For the height assertion,
multiplication on the product of the unit coefficient spheres
in two fixed finite degree spaces has a positive minimum norm,
since the polynomial ring is a domain. The integer cofactor
of a primitive integral factor has norm at least one.

For each \(G\) put \(J=xG_x-kyG_y\). If \(G\mid J\),
then \(J=\gamma G\) for a constant \(\gamma\), allowing
\(J=0\). Thus every monomial of \(G\) has \(i-kj=\gamma\)
and \(G(x,y)=x^\gamma P(x^ky)\) on the positive chart
for a nonzero polynomial \(P\) of degree at most \(D_y\).
This case contributes at most \(D_y\) exact outputs.
If \(G\nmid J\), their common zero set is finite by Bézout.
In fact if they shared a geometric irreducible component,
all conjugate components of \(G\) would also divide \(J\),
and hence \(G\) would divide \(J\).

At any remaining sample point, \(J(x,y)\ne0\). Clearing
only its evaluation denominator gives
\(|J(x,y)|\gg X^{-S}\). Put \(g(t)=G(t,f(t))\).
On \([1/2,3]\) its second derivative is \(O(X^{2S})\).
Since \(G(x,y)=0\) and \(y=f(x)+O_k(X^{-\tau})\),
\[
         |g(x)|\ll X^{2S-\tau},\qquad
         \left|g'(x)-J(x,y)/x\right|\ll X^{2S-\tau}.
\]
The bound \(8S<\tau\) makes \(|g'(x)|\gg X^{-S}\).
On the interval of radius \(X^{-\tau+4S}\) around \(x\)
the derivative has constant sign and remains of this
order or greater. Its variation there is \(O(X^{6S-\tau})\).
The change in \(g\) from the sample point to either endpoint
is at least a constant times \(X^{-\tau+3S}\), so \(g\)
has a zero in this interval. Multiplying \(g(t)\) by
\(t^{kD_y}\) makes a polynomial of degree at most
\(D_x+kD_y\). It is nonzero in this case because its
Laurent polynomial has nonzero derivative at the sample.
There are thus \(O(1)\) possible nearby roots. Each root
can be assigned to at most one distinct \(x\), since
\(\tau-4S>2R\) and the normalized \(x\)-height gives
spacing \(\gg X^{-2R}\). For fixed \(x\), the same
spacing argument gives at most one \(y\) in the strip
because \(\tau>2E\).

Each factor and each interval consequently contributes
\(O(1)\) exact outputs. Undo the normalization and sum
the intervals and dyadic boxes to obtain the stated count.

## Relative coefficient covers with a real gap

**FRP2.** In the setting of APC4, replace the third divisor
condition by the requirement that the third original
letter \(W\) is at most \(X^{1-\tau}\), where \(\tau>0\)
is fixed. Then it is the short summand for sufficiently
large \(X\) in the shell \(X/2<c\le X\). On the tall
roles \(U,V\) one has \(U/V=1+O(X^{-\tau})\).
Fix the signed projected coefficient record and its
relative-height subset and choose an actual representative
\(Q_0\). For \(z=Q/Q_0\), the exact quotient \(r^kz\)
is within relative \(O(X^{-\tau})\) of \(1/Q_0\).
Use \(R=C_r/k\), \(E=c\min(1,kL)\). Then
\(RE\le C_r cL\). Fix \(b=\eta L/4\), choose \(c\)
small, and choose the fixed lower cutoff \(k_0\)
large in terms of the strip gap. The preceding
count gives \(O(n^{\eta/4}X^\varepsilon)\) distinct
outputs at each order and in each subset.

Orders up to \(J/L\) are fixed in number before the
height limit. Above \(J/L\), enumerate \(r\) at
cost \(O(n^{2C_r/J})\) and use the real strip and
rational spacing to determine \(z\), with \(2c<\tau\).
Nontrivial \(r\) gives only \(O(\log X)\) possible
orders and the root-one case is collapsed as in (APC4).
Primitivity recovers the original row from each exact
output. The light branch now costs \(at+b_1+\eta\),
since the cover is required only after the projected
coefficient record and the selected pool. The heavy
branch is unchanged and uses (EB5) on the full original
radical budget. The resulting exponent in (APC4b)
has \(u=b_1\) in place of \(\mu+b_1\), including
when the complete short summand varies.

Combining this with the congruence proof gives one
statement for any fixed complete letter \(W\), under
the original three caps. If \(W\le X^{1/2}\), use
the real strip with \(\tau=1/2\). If \(W>X^{1/2}\)
and \(\operatorname{rad}(W)\le X^{CL}\), choose
a fixed \(\kappa>0\) with \(C\kappa<1/8\) and
take the full part of \(W\) having valuations at
least \(\lfloor\kappa/L\rfloor\). The omitted
part has height at most \(X^{1/8}\); the fixed
full part has height at least \(X^{3/8}\) and
scaled depth at least \(\kappa/2\) at small
budgets. Apply the light-branch congruence
argument of (APC4) with this one fixed modulus
and pay no divisor record. Thus after fixing the
complete letter, the per-center relative
completion cost is arbitrarily small uniformly
over its full range of numerical heights.
[Conditional covers and output congruences (APC4)](ABC_ADAPTIVE_POWER_COEFFICIENT_PRICES.md#price-of-a-cover-by-small-relative-coefficient-heights)



## All close pairs in one fixed complete-letter fibre

**FRP3.** Fix the value and role of a complete letter \(W\) in the shell
\(X/2<c\le X\) of genuine positive primitive triples. Fix the
orientation of the other roles, so that
\[
                       W=\sigma U_i+\tau' V_i,\qquad
                       \sigma,\tau'\in\{1,-1\}.
\]
Suppose \(\operatorname{rad}(W)\le X^{CL}\), with \(C\) fixed.
Allow the representations and common order to be
chosen for each pair. For each represented row use
\[
           U_i/V_i=Q_i r_i^k,\qquad H(r_i)\le X^{C_r/k},
\]
where all coefficient and root primes lie on the two original
free roles. The height cap for \(Q_i\) can be \(X^{C_Q}\).
For a fixed allowance \(\eta>0\), choose a fixed \(k_0\) large
and \(c_0>0\) small in terms of \(C,C_r,\eta\).
For all sufficiently small fixed \(L\), the number of ordered
pairs of distinct rows admitting a common order \(k\ge k_0\)
and satisfy
\[
                  H(Q_1/Q_2)\le X^{c_0\min(1,kL)}
\]
is \(O_{C,C_r,C_Q,L,\eta,\varepsilon}
(n^\eta X^\varepsilon)\), uniformly in the fixed value \(W\).
All original caps on the free roles may also be kept.
The exact cross-ratio count handles all eligible
pairs together, with no per-center catalogue.

Indeed, put \(r=r_1/r_2\), \(z=Q_1/Q_2\), and let \(w=r^kz\)
be the exact cross ratio of their complete odds. If
\(W>X^{1/2}\), its fixed high-depth part \(D\), as above,
gives \(w\equiv1\bmod D\). Apply the distinct-output
congruence argument of (APC4) with root height constant
\(2C_r\). If \(W\le X^{1/2}\), both free letters of
each row are tall, with \(U_i/V_i=1+O(X^{-1/2})\).
Thus \(w=1+O(X^{-1/2})\), and the real-strip lemma
gives the distinct-output bound with the same new
height constant. In both cases choose an allowance
strictly smaller than \(\eta\). There are only finitely
many orders up to \(J/L\) before the height limit.
For orders above \(J/L\), count \(r\) at cost
\(O(n^{4C_r/J})\) and determine \(z\) by rational
spacing in the fixed congruence or real strip.
The orders with either original \(r_i\ne1\) number
\(O(\log X)\). When both roots are one, count their
outputs in a single order state using the same
relative height bound. Summing gives
\(O(n^\eta X^\varepsilon)\) distinct outputs with
allowance left for height-subpower factors.

One now recovers the original pair at height-subpower
cost from each \(w\), using the fixed-target cross-ratio
identity. Write \(w=A/B\) in lowest terms, and
\(s_i=\sigma U_i,\ t_i=\tau' V_i\), so \(s_i+t_i=W\).
The signed fractions \(s_i/t_i\) have cross ratio \(w\).
Since \(w\ne1\) for distinct primitive rows, the identity
\[
           ((A-B)s_1-AW)((A-B)s_2+BW)=-AB W^2
\]
determines \(s_1,s_2\) up to divisor-many choices,
and hence both rows. The right side has nonzero
polynomial height. This is the fixed-target
specialization of
[complete cross-ratio recovery (TR1)](ABC_DIFFUSE_FACTOR_GRAPH.md#recovery-from-a-matched-target-ratio).
The identity applies to signed summands by expansion.
This last recovery step is what avoids paying a
separate reference catalogue for the close pairs.


### Removing close pairs over a target dictionary

Suppose a family has at most \(n^{\alpha+o_X(1)}\) different
complete \(W\)'s in that role. Summing FRP3 over this actual
dictionary gives \(O(n^{\alpha+\eta}X^\varepsilon)\)
eligible ordered close pairs overall. Deleting one row
from each nonempty close-pair edge therefore removes at
most this many original rows. If the original family has
exponent \(d>\alpha\) by a fixed margin, choose the inner
allowances and \(\eta\) smaller than that margin. A family
of the same exponent remains, with all distinct pairs
in an original \(W\)-fibre separated at the stated width.
This applies to full centered coefficients when the
roots are formed by nearest rounding, since \(kq_e\le2e\)
on every retained root prime. After fixing actual full
parts, it applies to the remaining coefficients whenever
their quotient gives the original complete ratio through
the one fixed prefactor of that record. The pair bound
at occurrence uses the actual dictionary price and need
not leave a large family if \(\alpha=d\).

## Moving targets with a common full divisor

**FRP5.** Allow the fixed-letter value \(W\) of FRP3
to vary between pairs, retaining within each pair
the same complete letter \(W\) and its original role.
Assume either \(W\le X^{1-\tau}\) on the selected
pairs for fixed \(\tau>0\), or all target values
in the count have a prescribed common divisor
\(D\ge X^w\) with scaled orders of \(D\) bounded
below by \(\sigma>0\). Under the corresponding
FRP1 or FRP3 width and order hypotheses, the
count of distinct ordered close pairs is
\(O(n^\eta X^\varepsilon)\) *over all such \(W\)*.

Indeed the proof of FRP3 bounds the exact cross
ratio \(A/B\) using one fixed real strip about
one or one fixed unit congruence at \(D\).
This step is uniform in the moving complete
letter. In the real case the short letter is
a summand and both free terms are tall.
In the modulus case the signed ratios
\(s_i/t_i\) are units congruent to \(-1\)
modulo \(D\), where \(s_i+t_i=W\).
For either case, primitivity gives
\(\gcd(W,s_1t_2)=\gcd(W,s_2t_1)=1\) and
\[
        s_1t_2-s_2t_1=W(s_1-s_2).
\]
Thus the gcd removed in reducing the ratio
is prime to \(W\), and \(W\mid A-B\).
For distinct rows \(A-B\ne0\), and the divisor
bound gives height-subpower possible targets
\(W\) from one exact ratio. Apply the fixed-target
factorization in FRP3 after fixing one such \(W\).
If the two summands in a pair on the same target
also agree modulo its prescribed \(D\), the
cross difference is divisible by \(D^2\).
The denominator factors remain units and the
count can use this larger full modulus.

Thus a catalogue of \(n^{\mu+o_X(1)}\) prescribed
deep divisors costs \(n^{\mu+\eta}X^\varepsilon\)
for all close distinct pairs, whenever an eligible
target can be assigned one of these divisors.
This can be combined with the actual complete-target
price from FRP3 by taking the smaller count.
With a fixed positive compact nonatomic band on
the target role, greedy packing supplies arbitrarily
cheap full divisors of fixed positive numerical
weight on the primitive target. The width and
fullness thresholds in the pair estimate are then
chosen using that weight before the budget limit.
This permits deletion of a power-small set of
original rows to remove all such close pairs
in those charts, even when the complete-target
dictionary has full exponent. It provides
a sparsity test for same-target pairs and does
not infer a positive number of such pairs at
the occurrence endpoint.
[Cheap full divisors in a positive nonatomic band](ABC_ADAPTIVE_POWER_COEFFICIENT_PRICES.md#a-one-grid-estimate-in-a-genuine-fixed-target-fibre)
· [Complete cross-ratio recovery](ABC_DIFFUSE_FACTOR_GRAPH.md#recovery-from-a-matched-target-ratio)

## A doubled full modulus for two moving queries

Use the RT5 orientation of a residue cell. A complete common
query \(V\) is divisible by its prescribed chart divisor \(F\),
and the fixed distinct anchors \(b_1,b_2\) are congruent modulo
\(F\). The computed values \(W_i(V)=\epsilon V+\epsilon' b_i\)
are positive and are units at \(F\) by original primitivity.
Put \(\Delta=W_2(V)-W_1(V)=\epsilon'(b_2-b_1)\ne0\),
and write \(|\Delta|_F\) for the full part of \(|\Delta|\)
supported on the primes dividing \(F\). Set
\[
                   D=F|\Delta|_F,\qquad
                   R(V)=W_1(V)/W_2(V).
\]
For any two queries \(V,V'\) the exact cross difference
is
\[
  W_1(V)W_2(V')-W_2(V)W_1(V')
                                =\epsilon\Delta(V-V').
\]
All its denominator factors are units at \(F\).
Consequently \(R(V)/R(V')\equiv1\pmod D\).
The modulus contains \(F^2\) and has the same radical
as the original chart divisor. Its numerical size can
use the actual extra orders of \(\Delta\) at those primes.
This calculation uses the original query role, before
any normalization of the computed weighted difference.
[Role switch preserving the chart lock (RT5)](ABC_BALANCED_FIXED_SHIFT_REVERSE_TRANSFER.md#a-role-switch-that-preserves-the-full-shift-divisor)

**FRP4.** The number of distinct-query pairs with small
relative coefficient height has a uniform subpower bound
in this configuration.
Assume all original complete values are at most \(X^B\)
for fixed \(B>0\). For fixed \(w,\sigma>0\) require
\[
                  F\ge X^w,\qquad
               p^a\parallel F\ \Longrightarrow\ La\ge\sigma .
\]
For each ordered pair of distinct queries suppose its
exact output has a representation
\[
         R(V)/R(V')=z r^k,\qquad H(r)\le X^{A/k},\qquad
                  H(z)\le X^{c\min(1,kL)},
\]
with positive rational units \(r,z\) at \(D\), a fixed
\(A>0\), and an order \(k\ge k_0\) that may vary by pair.
For example, if two computed-operand quotients on this
pair have full representations \(R(V)=Q_V r_V^k\)
and \(R(V')=Q_{V'}r_{V'}^k\), with
\(H(r_V),H(r_{V'})\le X^{A/(2k)}\), the condition is
\(H(Q_V/Q_{V'})\le X^{c\min(1,kL)}\). Fixed
coefficient factors common to both queries cancel.
The congruence argument uses a representation of the
entire quotient. A bound on an observable omitting
some valuation slots suffices only if their remaining
factor can also be absorbed in the stated numerical
height and unit bounds.

For every fixed \(\eta>0\) there are constants
\(c>0,k_0,L_0\), depending on \(B,w,\sigma,A,\eta\),
such that the number of these ordered distinct pairs
at sufficiently small fixed \(0<L<L_0\) is
\(O(n^\eta X^\varepsilon)\). The bound is uniform
in the fixed anchors and chart divisor, retaining both
original genuine rows of each query and all their caps.
It also works with \(D\) replaced by any prescribed
divisor of \(F|\Delta|_F\) satisfying the corresponding
numerical size and full-order lower bounds.

First count distinct exact outputs. We have
\(X^{2w}\le D\le 2X^{2B}\), and \(D\) is
\(K\)-full for \(K=\lfloor2\sigma/L\rfloor\).
For a fixed amplification order \(\ell\) the two
relative coordinates of amplified outputs have
height constants
\[
 C_1=\max(1,\ell K A/(2wk)),\qquad
 C_2=\max(1,\ell K c\min(1,kL)/(2w)).
\]
The ratio \(C_1C_2/K\) is at most
\[
                \frac1K+\frac{\ell A}{2wk}
                  +\frac{\ell c}{2w}
                  +\frac{\ell^2\sigma A c}{2w^2}.
\]
Choose \(\ell\) large in terms of the desired
count allowance and \(B/\sigma\). Then choose
\(k_0\) large and \(c,L_0\) small so the displayed
ratio is below the absolute threshold of the
two-base unit congruence count. The outputs have
polynomial height, and their fixed residue is
one modulo \(D\), so the amplified output count
is \(O(D^{1/(\ell K)+\varepsilon'})\) per fixed
order. Since \(KL\ge\sigma\) at small budgets,
this is \(O(n^{\eta/4}X^\varepsilon)\) with the
allowances chosen suitably. For \(k\le J/L\)
there are finitely many orders at each fixed
budget, and the constants can depend on this
finite list.
[Rational monomial output amplification](ABC_SIGNED_VALUATION_PHASES.md#a-rational-congruence-bound-for-distinct-exact-outputs)

For \(k>J/L\), choose the rational \(r\) of height
at most \(n^{A/J}\) at cost \(O(n^{2A/J})\).
The residue \(z\equiv r^{-k}\pmod D\) gives at
most one \(z\), since two such reduced fractions
of height at most \(X^c\) have cross determinant
of absolute value at most \(2X^{2c}<D\) for
\(c<w\) and large \(X\). Put \(2A/J<\eta/4\).
If \(r\ne1\) there are only \(O_A(\log X)\)
orders, since the rational height \(H(r)\) is
at least two. If \(r=1\) then \(z\equiv1\bmod D\)
and the same uniqueness gives the exact output
one. For distinct queries this cannot occur:
\(W_1/(W_1+\Delta)\) is injective because
\(\Delta\ne0\).

It remains to count pairs after fixing an exact
output \(a/b\ne1\) in lowest terms. Put
\(x=W_1(V)\) and \(y=W_1(V')\). Then
\[
     \frac{a}{b}=\frac{x(y+\Delta)}{y(x+\Delta)},\qquad
      ((a-b)x+a\Delta)((a-b)y-b\Delta)=-ab\Delta^2 .
\]
The second equality is the split factorization of
TR1 with the signed shift fixed in both queried
pairs. The nonzero product has polynomial height.
Its divisor-many signed factorizations give \(x,y\),
and then the original \(V,V'\), since the two
computed functions have slopes of absolute value
one. No unknown target product or normalization
enters the count. Sum the divisor allowance and
the finite-order and tail counts to get the stated
pair bound. In a codegree of exponent \(q>0\)
this rules out an exponent-full population of
such small relative coefficient pairs, whenever
the doubled lock and complete ratio presentation
are available.
[Split recovery for a complete cross ratio (TR1)](ABC_DIFFUSE_FACTOR_GRAPH.md#recovery-from-a-matched-target-ratio)


### Exact query gcds retain the full primitive shift

There is a stronger version when the monomial
coordinates are units at the retained factors of
the following larger modulus. For the two queries
write \(g=\gcd(W_1(V),W_2(V))\) and
\(g'=\gcd(W_1(V'),W_2(V'))\). Each divides
\(|\Delta|\) and is coprime to \(F\), so the pair
\(g,g'\) can be fixed at height-subpower cost.
Put \(J=|\Delta|/\operatorname{lcm}(g,g')\).
Each query's normalized computed pair is coprime
and its difference is \(\Delta/g\) or
\(\Delta/g'\). It is therefore unit at \(J\),
and also at \(F\). Moreover
\(F\gcd(g,g')\mid V-V'\), using the fixed
chart divisor and \(W_1(V)-W_1(V')=\epsilon(V-V')\).
After dividing the determinant by \(g g'\)
one obtains the congruence
\[
                         R(V)/R(V')\equiv1\pmod{FJ}.
\]
The modulus still contains \(F^2\), since
\(\operatorname{lcm}(g,g')\) is prime to \(F\).
It also retains the common primitive shift
primes outside \(F\). Use its \(K\)-full part,
with \(K=\lfloor2\sigma/L\rfloor\), to retain
all of \(F^2\) and any additional powers of
sufficient order. This gives a prescribed
polynomial-height modulus in each divisor
state, with no additional radical dictionary.
The same pair bound applies if \(r,z\) are
units on this full part, as they are for
canonical decompositions of the two normalized
computed pairs using the same order \(k\).
The signed factorization still uses the
original fixed shift \(\Delta\) and recovers
the two original queries.


### A short real window for moving queries

There is also a real version of the distinct-query
count. Keep the fixed distinct anchors and computed
affine values, and suppose the pairs under consideration
satisfy
\[
        |\Delta|\,|V-V'|\le
                    X^{-\tau}W_2(V)W_1(V')
\]
for a fixed \(\tau>0\). For example this holds up to a
fixed multiplicative constant when both factors on
the right are of order \(X\) and the queries belong
to an interval of length \(X^{1-\tau}\) in the
original height shell. The displayed cross difference
then gives \(R(V)/R(V')=1+O(X^{-\tau})\). Under the
same two-base height assumptions as FRP4, use FRP1
with \(R=A/k\), \(E=c\min(1,kL)\), and \(b\) a small
multiple of \(L\). For \(k>J/L\), enumerate \(r\) and
use rational spacing in the real strip for \(z\).
The exact nontrivial cross ratio recovers the two
queries by the same fixed-shift factorization.
This gives \(O(n^\eta X^\varepsilon)\) such
pairs for sufficiently small \(c,L\) and large
fixed \(k_0\), with no chart modulus needed.

### Consequence for actual switched residue cells

After fixing the two anchors, suppose a genuine
codegree has a conditional cover of price \(b\)
by classes of full coefficient ratios, each of
sufficiently small relative diameter at an
admissible common order. Assign each query to
one class. If there are \(S\) queries, the
number of distinct ordered close pairs is at
least \(S^2/n^{b+o_X(1)}-S\). Under the original
unit-lock conditions, FRP4 bounds it by
\(O(n^\eta X^\varepsilon)\) for arbitrarily small
fixed height-budget allowance, with width
chosen accordingly. Hence a fixed margin
between the codegree exponent and \(b\) is
impossible. In the RT5 switched cell, a
conditional price strictly below \(2u-c_W\)
therefore gives the codegree saving. The local
target-divisor data and anchors are prescribed
before applying the bound. In particular a
cover price on a larger projection must first
be converted to its actual conditional price.
An argument forcing such a cover or a
power-large population of close pairs in
actual balanced cells is still missing.
