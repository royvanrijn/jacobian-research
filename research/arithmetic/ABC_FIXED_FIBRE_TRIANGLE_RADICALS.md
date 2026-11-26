# Primitive difference triangles inside fixed abc fibres

Three distinct complete values on a fixed-letter fibre give three
distances on a line. Divide their distances by the common gcd to
obtain a positive primitive additive triple. A positive-exponent
fibre cannot have many triangles of small normalized height.
A full-exponent relation whose normalized triangle radicals meet the
critical cap gives a second family at the one-third coefficient,
for every positive fixed-fibre exponent.

## Nonatomic mass outside prescribed prime pools

Fix \(C>0\) and use primitive positive triples \(a+b=c\) with
\(X/2<c\le X\), individual radicals at most \(n^C\) and
\(n=X^L\), where \(L\) is fixed before \(X\to\infty\).
On one fixed complete \(T\)-fibre write \(U,V\) for the other
roles. A prime pool is supported on an integer at most \(X^K\)
for some fixed \(K\). One may prescribe different pools for \(U,V\).

Consider a height-first sequence of such fibres with limiting
normalized count exponent \(d>0\). Refine to anatomic boxes for
the parts of \(U,V\) outside their pools and take weak scaled
valuation limits there. The nonatomic masses \(c_{U,\mathrm{out}}\)
and \(c_{V,\mathrm{out}}\) of those limits satisfy
\[
               d\le \min(c_{U,\mathrm{out}},c_{V,\mathrm{out}}).
                                                               \tag{FT1}
\]
The fixed \(T\) and the prescribed pools can vary with height
and budget. In particular this assertion applies to pools chosen
from witnesses in the fibre before counting the remaining rows.

To prove this, fix the entire supported parts \(U_0,V_0\) in
the prescribed pools at height-subpower cost, and write
\(U=U_0 U_1,\ V=V_0 V_1\). The stated cost is uniform in the
chosen pools, so it is compatible with the inner height limit.
Do not put \(U_0,V_0\) in any of the moving monomials below.
Their numerical heights need not tend to zero. The lower
nonatomic inequality of (FL2) applies to \(U_1,V_1\) by the
following verification of its congruences and real strips.
[Uniform supported enumeration and prime prices](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Fixed-letter valuation bound (FL2)](ABC_FIXED_LETTER_RADICAL_CAPS.md#strict-positive-scale-inequalities-on-the-two-free-letters)

Suppose some residual letter, say \(U_1\), has positive
nonatomic mass. A small compact positive band supplies an
arbitrarily cheap full divisor \(F\mid U_1\) of fixed positive
numerical weight. Fix it. The signed original equation makes
\(V_1\) a prescribed unit modulo \(F\), since the constant
\(V_0\) and the fixed letter \(T\) are units there by original
primitivity. Freeze the residual parts of \(V_1\) outside its
selected atomic windows and zero window, at normalized cost at
most \(c_{V,\mathrm{out}}\) plus an arbitrarily small allowance.
For this grouping one can first freeze the tail above a large
fixed scaled order at arbitrarily small radical cost, then
select finitely many positive atoms accounting for all but
a small part of the remaining atomic mass. Shrinking their
windows gives a bounded list of powered roots of height
\(X^{O(L)}\) and a rational correction of numerical height
exponent tending to zero. The zero window contributes only
vanishing numerical weight. The monomial congruence bound
(GW) counts the remaining part with arbitrarily small
normalized loss. Its exact value determines \(V_1\), hence
the original row. This gives \(d\le c_{V,\mathrm{out}}\).
Thus \(V_1\) also has positive nonatomic mass, and the
symmetric argument gives \(d\le c_{U,\mathrm{out}}\).
[Grouped roots with a vanishing numerical correction (GW)](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)

If both residual nonatomic masses vanish, use atom grouping
on both residual letters at arbitrarily small freezing cost.
When \(T\ge X^{1/2}\), its radical cap supplies a full divisor
of fixed positive numerical weight and fullness of order
\(1/L\), as in (FL1). The signed equation prescribes a unit
residue for \(U_1/V_1\) at that modulus, including the fixed
unit ratio \(V_0/U_0\). Apply (GW) to the grouped residual
ratio. If \(T<X^{1/2}\), the two complete free letters are
the tall terms, differing by \(T\), so
\[
         \frac{U_1/V_1}{V_0/U_0}=1+O(X^{-1/2}).
\]
Use the relative real strip count (RUW1), which is uniform
in the positive centre \(V_0/U_0\). In both cases the exact
residual ratio determines the original row, since the
complete ratio and the fixed signed relation with \(T\)
do. This excludes a positive normalized count in the
both-atomic case.
[Fixed-letter full divisor or real gap (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
· [Real strip count uniform in its centre (RUW1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights)

## A small-feature bound on a single fibre

Fix \(r>0,\ C>0\) and a rank bound \(R\ge1\). There are
\(\theta,\delta,L_0>0\) such that, for every fixed
\(0<L<L_0\), the following fibre has
\(O_{r,C,R,L}(n^{r-\delta})\) rows. It is a fixed complete
\(T\)-fibre as above in which one free complete value \(W\)
is an affine real function of an assigned integer feature
vector \(z\in[-X^\theta,X^\theta]^R\). The affine function
can vary with the fibre and the height. In particular, the
steps of an arithmetic progression presenting \(W\) can
be arbitrarily large.

Suppose this fails through fixed budgets \(L_j\to0\),
with feature exponents and putative gaps tending to zero.
At every fixed budget there are arbitrarily large violating
heights. Take the inner height subsequences and refinements
before using the limiting valuation result (FT1). Choose a
small normalized lock cost \(p_0>0\) with \(Rp_0<r/8\).
We keep a class of exponent at least \(r-Rp_0-o(1)\)
during at most \(R\) rank descents.

If the current affine feature span has positive rank and
\(W\) moves, choose two witness rows with different
\(W\)-values and let \(D\) be their nonzero difference.
Fix the entire \(D\)-supported part of \(W\) at height-subpower
cost. Apply (FT1) to the outside parts at the current stage.
After anatomic refinements and a subsequence these have
nonatomic mass at least \(r-Rp_0\) in \(W\). Choose a
compact positive band in which a positive part survives
after removing neighborhoods of the finitely many atoms
of mass at least \(p_0/4\). These neighborhoods can be chosen
to lose arbitrarily little nonatomic mass. The remaining
band can be chosen as a finite union of continuity intervals.
For sufficiently late boxes the available individual prime
powers have normalized radical costs at most \(p_0/2\),
and their total cost exceeds \(p_0\). Otherwise a single
large mass in the selected compact band would give a
limiting atom of mass at least \(p_0/2\) there. A greedy
selection gives a full divisor \(F\mid W\) outside the
witness pool, with radical cost at most \(p_0\) and
\(F\ge X^{w_*}\) for a fixed \(w_*>0\) on this subsequence.
The positive lower band endpoint gives the numerical
weight. Fix one such value of \(F\) at cost
\(n^{p_0+o_X(1)}\) by radical enumeration.

If the new feature rows still span the previous affine
space of rank \(j\le R\), choose \(j\) independent
differences among retained rows and a nonzero minor
\(\Delta\) on \(j\) coordinate columns. These differences
are a linear basis for all differences of the previous
span. By integer Cramer's rule, \(\Delta\) times the
original witness difference vector is an integer combination
of the retained difference vectors. Applying the linear
part of the real affine function for \(W\) gives
\(F\mid \Delta D\). Since \((F,D)=1\), this forces
\(F\mid\Delta\), which contradicts
\[
             0<|\Delta|\le j!\,(2X^\theta)^j<X^{w_*}.
\]
The last inequality holds for sufficiently late fixed
budgets in the contradiction sequence and large enough
inner heights. Thus the affine rank drops.

At each rank one may re-extract the outside profiles after
the new witness pool is prescribed. The number of descents
is fixed and their total count loss stays below \(r/8\).
No uniform lower bound for a lock weight over unspecified
limiting profiles is being assumed. If \(W\) is constant
at any stage, including rank zero, it and \(T\) determine
a single row, contradicting the retained positive exponent.
This proves the uniform bound by the failure diagonal.
[Global complete-letter rank descent (CF1)](ABC_COMPLETE_LETTER_SMALL_FEATURES.md)

For each fixed doubling constant \(K\), a free complete-value
dictionary with doubling at most \(K\) on such a fibre obeys
the same bound, with \(\delta\) also depending on \(K\).
Green and Ruzsa embed it in a generalized arithmetic
progression of bounded rank and volume \(O_K(n^{C+o_X(1)})\).
Here the dictionary itself has size at most \(n^{C+o_X(1)}\)
by direct radical enumeration. The progression coordinates
fit the feature box for sufficiently small fixed \(L\).
[Green and Ruzsa, Theorem 1.1](https://arxiv.org/pdf/math/0505198)

## Primitive triangles of three values on the fibre

Let \(\mathcal A\) be the set of complete \(W\)-values on
one fixed \(T\)-fibre. Put \(N=|\mathcal A|\). For an increasing
triple \(x_1<x_2<x_3\) of values in \(\mathcal A\), set
\[
 d_{ij}=x_j-x_i,\qquad g=(d_{12},d_{23}),\qquad
 D_{ij}=d_{ij}/g\quad (i<j).
\]
Then \(D_{12}+D_{23}=D_{13}\) is a positive primitive equation.
Its height is \(D_{13}\), and its raw radical is
\(\operatorname{rad}(D_{12}D_{23}D_{13})\).

Suppose \(N\ge n^r\) for fixed \(r>0\). There are
\(\theta,\sigma,L_0>0\) such that for every fixed
\(0<L<L_0\) the triangles of normalized height at most
\(X^\theta\) number \(O_{r,C,L}(N^3 n^{-\sigma})\).
For fixed distinct \(x_1,x_3\), choose \(g\mid(x_3-x_1)\)
at height-subpower cost. Then \(x_2=x_1+g u\) with
\(u\) an integer in \([1,X^\theta]\). Apply the rank-one
fibre bound to this arithmetic progression and sum over
the fixed anchors. Choose the divisor allowance smaller
than its exponent gap. The promise \(N\ge n^r\) converts
that gap into the stated fraction of triangles.

Now fix \(B,K>0\) and \(0<\eta<2r\). Restrict to triangles
with raw radical at most \(n^B\). On all but a power-small
fraction of the full \(N^3\) triangle scale, each of the
three normalized sides has outside radical greater than
\(n^{2r-\eta}\). For a side tested, the pool may be any
prescribed integer at most \(X^K\) chosen using the two
endpoints of a different side, and independent of the
remaining query. Different side tests can use different
fixed choices of the complementary anchor side and pool
rule.

Indeed omit the short triangles first. Fix the endpoints
of that other side and record \(g\), a divisor of their
nonzero difference. The normalized value on the anchor
side is fixed, with raw radical at most \(n^B\). The other
two true quotient sides have raw cap at most \(n^B\).
On the tested side discount the given pool. Apply (FL1)
in dyadic shells of the normalized height \(Y\ge X^\theta\),
with the fixed complete anchor-side value. Use the fixed
budget unit \(L\) at each shell. A cap \(X^{eL}\) there
becomes \(Y^{\ell L}\) with \(\ell=e\log X/\log Y\),
bounded uniformly in \(B,\theta\). The pool has polynomial
height in \(Y\). This gives at most \(n^{r-\eta/2+\varepsilon}\)
choices of the free original value up to bounded orientations
and inner divisor losses. Indeed the tested directed
difference times \(g\) recovers the third value from its
fixed endpoint. Once the allowed normalized losses have
been chosen below \(\eta/4\), summation over anchors gives
the power saving. In particular for a common prescribed
pool, the total normalized outside radical exceeds
\(n^{6r-3\eta}\) on this remaining bounded-raw relation:
the three normalized sides are pairwise coprime.
[Fixed complete letters with pools (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
· [Height conversion for fixed-endpoint counts](ABC_FIXED_ENDPOINT_SHIFT_RADICAL_COUNTS.md#a-physical-difference-on-one-role)

## An incidence transfer at the critical triangle cap

For any set \(\mathcal Q\) of increasing triples from a
real set of size \(N\), the ratio dictionary
\((x_2-x_1)/(x_3-x_1)\) has cardinality
\[
                   \gg \frac{|\mathcal Q|^2}{N^4\log(2N)}.
                                                               \tag{FT2}
\]
Indeed two triples with equal ratio give three collinear
points of the Cartesian square, by pairing their first,
second and third entries. Each line in that square has
at most \(N\) points. The rich-line bound
\(O(N^4/k^3+N^2/k)\) from Szemerédi–Trotter, summed
dyadically with weight \(k^3\), gives
\(O(N^4\log(2N))\) collinear triples. Cauchy–Schwarz proves
(FT2). This is the ratio second-moment argument of
Murphy, Roche-Newton and Shkredov.
[Variations on the sum-product problem II, Lemmas 2.2 and 3.1](https://doi.org/10.1137/17M112316X)

For integer inputs the exact ratio determines the full
positive primitive difference triple. In the height-first
endpoint convention, suppose \(N=n^{r+o(1)}\) for fixed
\(r>0\), and a relation of \(n^{3r-o(1)}\) triangles has
raw radical at most \(n^{6r+o(1)}\). Apply the short-height
and side bounds at fixed exponents strictly below \(r\).
Discard short triangles and retain tall quotient triples.
By (FT2) their dictionary has exponent at least \(2r\).
The one-third bound in
height shells gives the matching upper exponent \(2r\).
Each individual output side has radical at most
\(n^{2r+o(1)}\) after removing the power-smaller cases in the
preceding section and using the total cap. To formalize the
height normalization, choose a dyadic output shell of upper endpoint \(Y\) retaining
the same exponent. At each fixed original budget \(L\), along
an inner subsequence \(\log Y/\log X\to H_L\in[\theta,1]\).
If the cap at that budget is \(n^{6r+\beta_L+o_X(1)}\),
choose a fixed output budget just above
\((6r+\beta_L)L/H_L\), with normalized rounding tolerance
tending to zero along the budgets. Here \(\beta_L\to0\)
and the choice is independent of the inner heights on the
selected subsequence. The ratio of the normalized count
exponent to the raw radical budget tends to \(1/3\) in the outer
limit. Thus the output is a balanced endpoint at the
new budget scale.
[Limiting one-third bound](ABC_OMITTED_COORDINATE_GRAPH.md#a-limiting-small-budget-coefficient-of-13)

This result applies also to the intermediate repaired
multiplicity range. In the central normalization
\(Z=Y_{\rm rep}^{1/M}=n^{1/3+o(1)}\), a target-radical fibre
of size \(Z^{u+o(1)}\) contains a complete-target subfibre
of size \(n^{u/3+o(1)}\). This follows by enumerating the
height-subpower many complete values on the fixed radical.
For every fixed \(u>0\) the short-triangle exclusion then
holds. A full triangle relation at raw cap
\(n^{2u+o(1)}\) transfers to the absolute radical scale described above.
This scale is smaller than the original one for \(u<1/2\).
The minimally repaired representation selects primitive quotient rows.
[Least-repair radical fibres](ABC_FIXED_LETTER_RADICAL_CAPS.md#size-of-a-minimally-repaired-target-radical-fibre)

The construction concerns normalized triangles in a single
complete fibre, where the primitive affine secant of two
original rows would have a zero coordinate. It supplies no
upper bound on the raw radicals of these triangles. Forcing
a sufficiently large relation at a cap below the critical
exponent \(6r\) remains a separate arithmetic question.

### Transfer of the high-powerful saving

On the critical radical scale of the triangles, any
primitive-output sector with a strict small-budget
coefficient loss has a strict count loss on the input
triangles. More precisely, a bound \(n^{2r-\beta+o(1)}\)
on the number of output triples within that sector, with
fixed \(\beta>0\), gives a triangle count of at most
\(n^{3r-\beta/2+o(1)}\) by (FT2). Discard the power-smaller
short normalized heights first so that the same output
height and budget conversion is available.

In particular, at the critical cap \(n^{6r+o(1)}\), the
triangle class in which two output letters meet the
high-fullness condition of the two-full-letter theorem
is power smaller than \(n^{3r}\). Apply that theorem in
the derived height shell at the rounded relative budget
from the preceding proof. The normalized output
exponent is at most \(2r-6r\delta+o(1)\) for the fixed
saving \(\delta>0\) in that theorem. Its positive gap
is transferred by (FT2). The full orders and heights
may be grouped by fixed meshes at each budget; using
fullness with the original theorem's fixed margin
absorbs rounding. The same argument applies to the
numerically short-letter sector. These are tests on
the normalized output differences of actual fixed-fibre
values, even when the three original letters of each
input have unrestricted valuation profiles.
[Two full letters at the small-budget scale](POWERFUL_TRIPLES_ABC_ENDPOINT_TRANSFER.md#corollary-a-small-budget-saving-with-two-full-letters)
· [Short-letter sector (SH1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-strict-saving-when-the-short-term-is-small)

The exact-cell refinement strengthens the pinned
transfer. It fixes the true scaling gcd on a
full-exponent subrelation, giving an injective
output with two minimum complete dictionaries
even at intermediate repaired multiplicity.
For a global tagged relation on minimum input
dictionaries, popularity of physical gaps gives
all three minimum output dictionaries.
[Exact scaling of selected triangle relations](ABC_HEREDITARY_MINIMUM_TRIANGLE_TRANSFER.md)
