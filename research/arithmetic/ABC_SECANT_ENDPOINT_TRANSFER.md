# Primitive secants and products at the radical threshold

A product of ratios from a balanced exceptional family can produce a
second exceptional family. At the radical threshold the output is
itself a balanced endpoint. Its two products retain explicit
information about the original degrees. For primitive secants this
transfer gives a lower bound on the nonatomic radical of the secant,
a real-height restriction, and a bound on its high-valuation tail. A
separate factorization argument counts fresh prime support in
polynomials of the same products, including exact-power secants.

The statements below describe consequences of a putative endpoint
family. The existence of a large set of small-radical secants remains
an additional hypothesis in the transfer.

## The height-first setting

Use a sequence of balanced endpoint boxes, with \(n=X^L\), \(L\to0\)
and \(L/\lambda\to1\), and \(n^{1/3-o(1)}\) distinct primitive rows of
height at most \(X\). Exponents in base \(n\) are taken first as
\(X\to\infty\) with a fixed box budget \(L\), then as \(L\to0\).
An inner subpower loss means \(X^\varepsilon\) for every
\(\varepsilon>0\), with \(L\) fixed. Each row has three complete
positive letters. Their radical exponents tend to \(1/3\) in base
\(n\), their numerical exponents in base \(X\) tend to \(h_T\), and
their scaled-degree measures tend weakly to \(\nu_T\). In exact
valuation notation the measure of a positive integer \(Q\) at budget
\(L\) is
\[
       \mu_{Q,L}=\sum_{p^v\parallel Q}\frac{\log p}{L\log X}\,
                                    \delta_{Lv}.
\]
Each \(\nu_T\) has mass \(1/3\). Write \(c_T\) for its nonatomic mass.

One may use exact valuation layers in taking these limits. Indeed use
the anatomic truncation in BBLT, Lemma 2.2, with parameter \(o(L)\);
outside the initial coefficient and the one packed coordinate its
roots are the actual squarefree valuation layers. For cutoff \(d\)
the radical of actual layers of degree \(>d\) has exponent at most
\(1/d\) in base \(X\). The discarded coefficient and the root of the
packed coordinate have exponent \(o(L)\). Thus the normalized weak
measures differ by \(o(1)\) in total variation for dyadically refined
caps, with analogous statements when the height is bounded by a
fixed power of \(X\).
[Anatomic reduction, BBLT Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)
· [Endpoint valuation measures and the pair inequalities](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)

Fix \(m\ge2\) and nonzero integers \(k_1,\ldots,k_m\). In slot \(i\)
choose an ordered pair of distinct letters. Write \(U_i,V_i\) for the letters
as they occur respectively in the numerator and denominator, after
exchanging their order if needed when \(k_i<0\). For an ordered
\(m\)-tuple of rows set \(P=\prod_i U_i^{|k_i|}\),
\(Q=\prod_i V_i^{|k_i|}\), \(g=\gcd(P,Q)\), and
\[
                          A=P/g,\qquad B=Q/g.
\]
The fresh part of an integer on such a tuple consists of its full
prime powers outside the supports of all three letters of all the
input rows.

### Small intersections and output multiplicities

For a prescribed polynomial-height prime pool, (DRR4) implies that
rows with shared radical more than \(n^\eta\) on the union of their
three letters have exponent at most \(1/3-\eta/3+o(1)\).
Use one of the three pairs of letters whose shared radical is at least
\(n^{2\eta/3}\). The pool may depend on the complementary tuple, since
the bound is uniform after the complement is fixed. In particular,
except on \(n^{m/3-\eta/3+o(1)}\) tuples, every two input rows have
shared radical at most \(n^\eta\). In subsequent applications a
family of exponent \(m/3-o(1)\) may therefore take the total shared
radical to have exponent \(o(1)\). Formally one first discards the
bad tuples for fixed \(\eta>0\) and decreases \(\eta\) along
a diagonal sequence of budgets. The same argument works for a fixed
finite collection of such tests.
[Uniform anisotropic overlap inequality (DRR4)](ABC_DIVISOR_RECORD_RECIPROCITY.md#anisotropic-intersections-for-an-arbitrary-record)

On these tuples each of \(A,B\) has radical exponent \(m/3+o(1)\).
Every prime exclusive to an input row and dividing either selected
letter survives with its full exponent multiplied by \(|k_i|\) in
the corresponding primitive product. The total radical mass lost or
merged is at most the shared radical mass, times a constant depending
on \(m\). Moreover one output pair has at most \(n^{o(1)}\)
preimages, including inner subpower allowances. To verify this,
first enumerate the joint shared radical as a squarefree integer at
most \(n^{o(1)}\). Then record the full part of each selected letter on
these primes. For each given radical the number of such records is
inner subpower by fixed-support enumeration. For the primes outside
this record allocate each power in \(A,B\) to its slot in at most
\(m^{\omega(AB)}=X^{o_X(1)}\) ways, and divide its exponent by
\(|k_i|\). This determines both selected letters in each slot and
hence the original triples. If a component is \(1\), it needs no
allocation.

The actual maximum output height exponent \(H\) can be binned and
has a positive lower bound. In fact each endpoint letter has positive
radical mass on a fixed compact positive degree interval, by (SH1).
Removing \(o(1)\) radical mass from that interval removes only
\(o(1)\) numerical weight. The remaining exclusive prime powers give a
positive numerical height in each of \(A,B\). The upper bound for \(H\)
is \(\sum_i |k_i|\). Binning \(\log_X\max(A,B)\) with any fixed mesh
costs no exponent; taking the mesh to zero along budgets gives a
limiting \(H>0\).
[Positive scales in every endpoint letter (SH1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#positive-scales-in-every-endpoint-letter)

### Control of nonfresh primes for fixed binary polynomials

Let \(F(Z,W)\in\mathbb Z[Z,W]\) be a fixed homogeneous polynomial of
positive degree with \(F(1,0)F(0,1)\ne0\). Except on a
power-smaller tuple set, the radical of the nonfresh part of
\(F(A,B)\) has exponent \(o(1)\). Here and below tuples with
\(F(A,B)=0\) can be discarded. With all but one input row fixed,
a nonzero fixed polynomial has boundedly many roots in the remaining
ratio. The map from a primitive positive triple to the ratio of its
two selected letters is injective. Thus zero evaluations account for
\(O_F(\prod_{j\ne i}N_j)\) tuples, where \(N_j\) is the number of
available rows in slot \(j\).

For the nonfresh claim first remove bad shared-prime tuples as above.
A prime exclusive to one slot and dividing a chosen letter there
cannot divide \(F(A,B)\), except possibly if it divides the fixed
number \(F(1,0)F(0,1)\). A prime exclusive to slot \(i\) on its omitted
letter makes \(U_i/V_i\) equal to a fixed sign modulo that prime.
Thus if it divides \(F(A,B)\), it divides the numerator of
\[
                   F\left(\chi_i^{|k_i|}
                            \prod_{j\ne i}(U_j/V_j)^{|k_j|},\,1\right),
                 \qquad \chi_i\in\{-1,1\}.
\]
All denominators here are units at the exclusive prime. The numerator
has polynomial height and depends on the complement alone. It is
nonzero except on a power-smaller set of complements: fix all but one
of the \(m-1\) remaining ratios and use the same bounded root argument.
This is where \(m\ge2\) is used. Apply the uniform overlap test with
this numerator as the prescribed pool. Summing over slots and adding
back the negligible joint shared radical gives the claim. The
complementary congruence is the one used in the fresh monomial
threshold (FR3).
[Fresh radicals in products of ratios (FR3)](ABC_POLYNOMIAL_TRANSFER_RADICAL_THRESHOLD.md#fresh-radicals-in-products-of-ratios)

## The limiting output profile

Take \(F=aZ+bW\) for fixed nonzero integers \(a,b\). Suppose a tuple
family of exponent \(m/3-o(1)\) satisfies
\[
         \operatorname{rad}\bigl(F(A,B)_{\rm fresh}\bigr)
                                      \le n^{m/3+o(1)}.            \tag{ST1}
\]
Dividing \((|aA|,|bB|,|aA+bB|)\) by \(\gcd(aA,bB)\) and ordering
the signs gives a positive primitive abc triple. The divisor is
supported on \(ab\) and is at most \(|ab|\). Bounded coefficients
and this division change measures by \(o(1)\) on all three letters.
The output maximum height is comparable to \(\max(A,B)\). The map
from \((A,B)\) to the primitive triple has bounded multiplicity,
because a labelled signed ratio recovers \(A,B\). The previous
section therefore gives \(n^{m/3-o(1)}\) distinct primitive outputs
of maximum height \(X^{H+o(1)}\) with total radical at most
\(n^{m+o(1)}\).

Apply the two-cheapest-letter bound (GD3) and the limiting one-third
count at the output. The output budget can be taken to be
\((m+o(1))L/H\), which tends to zero. The two pure letters already
have raw mass \(m/3+o(1)\) in base \(n\). The third must also have
raw mass \(m/3+o(1)\) on a family of the original exponent. Otherwise
a radical bin bounded below \(m/3\) by a fixed constant, taken
together with either pure letter, gives a smaller exponent by (GD3).
Consequently these outputs give a genuine normalized \(1/3\)
endpoint. For this step, one first uses finite bins in \(\log_X\) of
the height, then dyadic output shells and bins in the radical exponents, and in the anatomic root
coordinates at each fixed output budget. Constants and the auxiliary
height exponent allowances are fixed before sending \(X\) to infinity.
Then the bin errors are decreased along the budget sequence.
[Two-cheapest-letter bound (GD3)](ABC_GROUPED_ROOT_DICTIONARY_ENTROPY.md)

Write \(D_s(r)=sr\) for dilation by \(s\).
In the normalized output measures the two pure laws are
\[
       \frac1m\sum_{i=1}^m(D_{|k_i|m/H})_*\nu_{U_i},
       \qquad
       \frac1m\sum_{i=1}^m(D_{|k_i|m/H})_*\nu_{V_i}.                \tag{ST2}
\]
This follows for bounded tests from the exact exclusive valuation
multiplication, with \(o(1)\) mass on the discarded supports.
The third law is \(m^{-1}(D_{m/H})_*\zeta\), where \(\zeta\) is a
weak limit of the fresh degree measures of \(F(A,B)\) in base \(n\).
Tightness of these measures follows from the bounded numerical height.
Let \(c(\zeta)\) denote its nonatomic mass. The two output pair
inequalities with the third letter give
\[
 c(\zeta)\ \ge\
       \max\left\{\frac m3-\sum_i c_{U_i},\
                   \frac m3-\sum_i c_{V_i}\right\}.              \tag{ST3}
\]
Both sums refer to the actual letter choices in the input slots.
The output description uses weak degree measures and the actual
reduced height. It permits an escaped first moment in an input or an
output letter.

### Primitive cross secants

Take \(m=2\) and the same pair \(U,V\) in opposite orientations, so
\(P=U_1V_2\), \(Q=U_2V_1\) and \(F=A-B\).
Let \(h_D\) be the limiting numerical exponent of \(|A-B|\) in
base \(X\), after binning. The two pure output letters have the same
limiting numerical exponent \(H\), since \(P,Q\) have the same
limiting numerical exponent and the same divisor is removed. Both
carry the law \(\tfrac12(D_{2/H})_*(\nu_U+\nu_V)\).
For \(\Delta=c_U+c_V-1/3\ge0\), (VM1)–(VM2) and (ST3) yield
\[
\begin{split}
 c(\zeta)&\ge 2/3-c_U-c_V,\\
 H-h_D&\le CH\Delta,\\
 \left(h_D-\int s\,d\zeta(s)\right)
             +\int_{[HR_*/2,\infty)} s\,d\zeta(s)&\le CH\Delta .
\end{split}                                                     \tag{ST4}
\]
For the third line the output weak measure is
\(\tfrac12(D_{2/H})_*\zeta\), and its first moment, in height units
\(H\), is \(H^{-1}\int s\,d\zeta\). The escaped moment
includes any numerical mass lost by omitting the nonfresh part.
[Height and tail bounds from the pair excess (VM1)–(VM2)](ABC_VARIABLE_MESH_DEGREE_GROUPING.md#bounds-from-the-pair-excess)

The actual reduced height also satisfies
\[
          0\le h_U+h_V-H\le e_U+e_V,\qquad
          e_T=h_T-\int s\,d\nu_T(s).                              \tag{ST4a}
\]
Indeed \(\log_X g\) is bounded by the numerical contribution of the
shared primes on the two factors of \(P=U_1V_2\). For any fixed
scaled cutoff \(R\), the contribution at degrees at most \(R\) is
\(o(1)\) by the intersection bound. The limiting upper bound for the
contribution above \(R\), taken at a continuity point, is the sum of
the input numerical heights minus the two moments below \(R\).
Let \(R\to\infty\) to obtain (ST4a).

At the minimum input nonatomic total \(1/2\), every \(c_T=1/6\)
and the input letters have numerical height \(1\) with no escaped
moment. Thus (ST4a) gives \(H=2\). Any near-full relation with fresh
radical at the threshold \(n^{2/3+o(1)}\) has
\[
 c(\zeta)\ge1/3,\qquad h_D=2,\qquad
       \int s\,d\zeta=2,\qquad
       \operatorname{supp}\zeta\subseteq[0,R_*].
\]
In particular at that endpoint there can be at most \(o(1)\)
normalized numerical weight on nonfresh primes of the primitive
difference. The nonatomic lower bound is in radical units. It is
\(1/3\) here, compared with the \(r/2=1/6\) consequence of (SW3)
for a fibre of exponent \(r=1/3\) in its positive-scale setting.
[The earlier numerical window relation (SW3)](ABC_FINITE_WINDOW_SECANT_GRAPHS.md#a-lower-bound-for-diffuse-secant-cost)

The three pairs of original letters give the same raw signed
determinant up to sign, since the original letters satisfy an
additive relation. Dividing by the primitive gcd in each chart
affects only shared primes, hence the fresh part is identical
in the three charts. As a result the first bound in (ST4) also
holds with \(\min_{U\ne V}(c_U+c_V)\) in place of \(c_U+c_V\).
The other bounds apply in each chart with its actual binned
numerical exponents. The full-height and tail claims are statements
in the iterated limit; for example, an imposed fixed positive
deficit \(H-h_D\) cannot coexist with \(\Delta\to0\) in a
near-full threshold sequence.

## Fresh radical growth in polynomial values

Suppose \(F\) as above has \(s\) distinct irreducible homogeneous
factors over \(\mathbb Q\). For every fixed \(0<\eta< s m/3\),
the tuples satisfying
\[
                  \operatorname{rad} F(A,B)_{\rm fresh}
                              \le n^{s m/3-\eta}                 \tag{ST5}
\]
have power-smaller count: there is \(\delta>0\) such that at
sufficiently late input budgets the count is
\(O(n^{m/3-\delta}X^\varepsilon)\), for every \(\varepsilon>0\).
The implied constants may depend on the parameters held fixed at
the height limit. There is no upper bound on the raw radical of
\(F(A,B)\) in this statement.

First consider one factor \(G\). For any fixed \(\eta_1>0\),
tuples for which its fresh radical has exponent at most
\(m/3-\eta_1\) have power-smaller count. If they had exponent
\(m/3-o(1)\), the overlap and multiplicity argument would give
\(n^{m/3-o(1)}\) distinct output pairs with radical of \(G(A,B)\)
at most \(n^{m/3-\eta_1+o(1)}\). For a linear \(G\) the
two-cheapest-letter bound on the derived abc triple contradicts
that exponent. Its height is bounded below by a fixed power of
\(X\), and its full raw radical is \(n^{O(1)}\), so the
relative budget required by (GD3) tends to zero.

For an irreducible \(G\) of degree \(d\ge2\), fix a nonzero
integer value \(M=G(A,B)\). If \(c\) is the leading coefficient
of \(G(Z,1)\) and \(\theta\) one of its roots, then \(c\theta\)
is an algebraic integer and
\[
  \operatorname{Norm}_{\mathbb Q(\theta)/\mathbb Q}
           (cA-c\theta B)=c^{d-1}M.
\]
The ideal of the generator has at most
\(\tau(|c^{d-1}M|)^d\) choices. Every conjugate has size
\(O_G(X^{\sum|k_i|})\) and is at least an inverse fixed power of
\(X\), by the nonzero integral norm. In a fixed ideal all possible
generators differ by units. The logarithmic vectors of the units
form a discrete group with finite kernel, since there are only
finitely many algebraic integers with all conjugates bounded.
Consequently a box of side \(O_G(\log X)\) has
\(O_G((1+\log X)^d)\) such vectors. Each generator recovers
\(A,B\) since \(\theta\notin\mathbb Q\). The count per \(M\) is
\(X^{o_X(1)}\). Counting the values \(M\) by their radicals gives
a contradiction to the near-full pair exponent. This argument is
the fixed-field analogue of Pell recovery, and extends the
cyclotomic norm count used for minimally repaired quotients.
[Cyclotomic norm recovery](ABC_MINIMAL_REPAIR_CYCLOTOMIC_SPLIT.md#primitive-cyclotomic-factors-and-their-fixed-value-counts)

Thus outside a power-smaller set every factor has fresh radical
larger than \(n^{m/3-\eta_1}\). The supports for different factors
are disjoint off a fixed finite set of primes, by their
resultants. To see this, for \(p\nmid B\) a common prime divisor
gives a common root of the two dehomogenized polynomials.
For \(p\mid B\), primitivity gives \(p\nmid A\), so such a prime
must divide a fixed leading coefficient. Add the distinct
normalized radical costs and take \(\eta_1<\eta/s\). This
proves (ST5). In this proof and for the fixed profile restrictions
after (ST4), failure of a uniform positive exponent gap would
give a sequence with losses tending to zero at late fixed
budgets, contradicting the associated near-full exclusion.

### A power-class restriction on small-radical secants

The exponent in this consequence is uniform over exact power orders
chosen separately at each budget and even each tuple. For each fixed
\(0<\eta<m/3\), apart from a power-smaller tuple set, simultaneously for
all \(q\ge2\) for which \(A/B\in(\mathbb Q^\times)^q\) and \(A\ne B\),
\[
        \operatorname{rad}\bigl((A-B)_{\rm fresh}\bigr)
                             > n^{(m/3-\eta)\tau(q)}.             \tag{ST6}
\]
In particular at \(m=2\), a relation at the threshold \(2/3\)
cannot concentrate in power-class equality, including classes and
orders that vary with the pair.

Suppose to the contrary there is a near-full exponent of tuples
failing (ST6). Discard the bad intersections. Exact \(q\)-th powers
\(A=A_0^q\), \(B=B_0^q\) satisfy
\(\log_X A\ge q\log_X\operatorname{rad}(A)\), and the pure radical
lower bounds show \(q=O(1/L)\). At each fixed budget we may select
one such \(q\) by pigeonholing, without losing exponent. On the
selected sequence one has \(qL\to0\). Otherwise pass to a subsequence
with \(qL\to b>0\). The exact degree measures of \(A,B\) are supported
on \(qL\mathbb Z_{\ge0}\). Their weak limits must be supported on
\(b\mathbb Z_{\ge0}\), by testing compact sets separated from that
lattice. But the unrescaled versions of (ST2), obtained solely from
exclusive valuations, are finite sums of fixed positive dilations of
input measures. At least one pure limit has positive nonatomic mass,
since \(c_{U_i}+c_{V_i}\ge1/3\) for every slot. This is a contradiction.

By the cyclotomic factorization of \(A_0^q-B_0^q\), distinct factors
\(\Phi_d(A_0,B_0)\), \(d\mid q\), have disjoint prime support outside
\(q\). After removing these primes, some factor has fresh radical at
most \(n^{m/3-\eta}\). They have negligible radical in the height
limit at fixed \(q\). The nonfresh control for \(F=A-B\) also applies
to every factor, so the chosen factor has total radical at most
\(n^{m/3-\eta+o(1)}\). Select its index at each budget. For nonlinear
factors the fixed-field norm argument in (ST5) counts the possible
pairs at that radical cost. The constants can depend on \(q\)
since the height limit comes first. For either linear factor use the
ordinary two-cheapest-letter bound. The input pair \(A_0,B_0\) has
the same radicals as \(A,B\) and height exponent comparable to
\(1/q\) in base \(X\). The relative raw budget of the triple is
\(O(qL)=o(1)\), while the absolute saving in base \(n\) remains a
fixed fraction of \(\eta\). The same finite output height and radical
bins justify this application. Either case contradicts near-full
multiplicity-adjusted exponent.

There is also a degree-residue consequence without assuming exact
powers. Fix \(b>0\). On a near-full tuple sequence with any adaptive
\(q\ge b/L\), the sum of the radical exponents on primes of \(A,B\)
with valuations not divisible by \(q\) is at least
\[
                    \sum_i(c_{U_i}+c_{V_i})\ \ge\ m/3             \tag{ST7}
\]
in the limit inferior. To see this, extract a subsequence realizing a
putative lower limit and take weak limits of the submeasures on
divisible valuations. The portion with \(qL>C\) contributes at most
\(O(1/C)\) divisible mass by the first moment bound. For any
fixed \(C\), the orders with \(b\le qL\le C\) have finitely many
choices at each budget, so a near-full relation in this range
can fix \(q\) before taking the height limit. Along a subsequence
with \(qL\to a\), the divisible limits are supported on
\(a\mathbb Z_{\ge0}\). They are dominated by the
respective full weak measures, and can use only those measures'
atomic mass. Tightness of all four submeasures makes the
complementary limiting masses at least the two nonatomic masses.
The equality of a pure limit's nonatomic mass with the sum of those
of its slot contributions follows from positivity of the measures.
For a near-full tuple family the small-overlap restriction may be
imposed first, so this argument also gives a power-smaller set of
exceptions to any fixed positive gap below the bound.

## Relation to exceptional-set and powerful-triple counts

The anatomic reduction is from Bernert, Browning, Lichtman and
Teräväinen. The whole-value modulus certificates in Zeng's
*Complete modulus locking* bound a single exceptional-set count
through a finite exponent program. Browning and Verzobio's
*Sums of three powerful numbers* counts primitive triples with
prescribed integer fullness thresholds by a coefficient-uniform
generalized Fermat surface estimate. The results above instead
derive constraints on products of independently chosen primitive
input ratios, by counting their primitive outputs and charging all
three original supports. The fixed-field proof reuses the norm
mechanism in the repository's minimally repaired cyclotomic
sector. It does not estimate the whole minimally repaired quotient
count.
[Anatomic reduction](https://arxiv.org/html/2410.12234v2#S2)
· [Complete modulus locking, Theorem 1.1](https://hexagonmath.org/pdf/2611.00621v1)
· [Powerful triples, Theorems 1.1–1.2 and 4.1](https://arxiv.org/html/2608.24512v1)
· [Minimal repair and the primitive quotient](POWERFUL_TRIPLES_ABC_ENDPOINT_TRANSFER.md)

## Remaining question

There is no bound here forcing enough unrestricted input pairs
to have fresh radical at most \(n^{2/3+o(1)}\). Even given that
cap, (ST4) allows secants whose height is full, whose scaled
degrees are bounded but spread out, and whose exact power classes
usually differ. Weak convergence alone cannot supply (ST6)'s
power-class hypothesis: for any fixed \(q\), the grids
\(Lq\mathbb Z\) and \(L(q\mathbb Z+1)\) approximate the same
positive-scale weak laws. The transfer and factorization results
specify two arithmetic constraints on any proposed argument
using such a secant relation.
