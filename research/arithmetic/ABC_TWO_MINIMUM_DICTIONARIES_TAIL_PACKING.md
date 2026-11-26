# Two minimum dictionaries and normalized triangle profiles

Two complete dictionaries at the minimum size control the valuation
profile of the third letter. The third dictionary can be larger.
This applies to critical difference triangles from a fixed complete
fibre of any positive exponent. Pinning the opposite vertex makes
a specified side the letter to be tested. A second estimate controls
relative short intervals and high-fullness residue classes in an
original fibre. It applies to the affine images of the triangle
sides before any count solely on the primitive outputs.

## Prime-power packing with two dictionary caps

Put \(n=X^L\), with \(L>0\) fixed before the height limit. Fix
\(r,\omega>0\). Consider positive primitive triples, with maximum
height in a dyadic shell \(Y/2<\max W\le Y\), where
\(X^\omega\le Y\le X\). Assume the three individual radicals
are at most \(n^{2r+\kappa}\). Fix three distinct roles \(A,B,C\)
and impose complete dictionary caps \(n^{r+\kappa}\) on \(A\)
and \(B\). There is no dictionary cap on \(C\).

Let the graph have \(E\ge n^{2r-\tau}\) edges. Select \(m\ge1\)
distinct rows in role \(A\) with common columns of role \(B\).
The count of ordered stars with distinct rows gives at least
\[
                \gg_m n^{\,r-m\tau-(2m-1)\kappa}                 \tag{TM1}
\]
common columns for some row list. To see this, write \(N_A,N_B\)
for the sizes of these two dictionaries. The mean column degree
is at least \(E/N_B\to\infty\) if \(\tau+\kappa<r\). Columns of
degree at least \(2m\) carry a fixed positive fraction of the
edges for large \(X\). Their count of ordered choices of \(m\)
distinct rows is \(\gg_m E^m/N_B^{m-1}\) by the power mean.
Divide by at most \(N_A^m\) row lists.

Write the complete signed \(C\)-value in row \(i\) as
\(C_i(t)=\pm t+a_i\), where the signed column value is \(t\).
The slopes are equal and the offsets are distinct. Suppose the
column also lies in a fixed residue class modulo \(G\le X^{K}\)
and has absolute height at most \(X^{h_B+\zeta}\). Here \(K\)
is fixed and \(\zeta>0\) may be made arbitrarily small.
If \(G\ge X^\beta\), the integer parameter in
\(t=t_0+Gv\) lies in an interval of length
\(O(X^{h_B-\beta+\zeta})\). Put the primes of \(G\) and of
\(\prod_{i<j}|a_i-a_j|\) in one prescribed pool \(P\). It has
polynomial height, uniformly for a fixed \(m,K\). Discard the
columns having normalized radical overlap at least \(\xi\) with
\(P\) on any of the \(C_i(t)\). For that fixed row, (FL1) with
outside \(C\)-cap \(2r+\kappa-\xi\) counts at most
\[
                   O(n^{\,r+\kappa/2-\xi/2+\epsilon_1})         \tag{TM2}
\]
such columns, with any fixed \(\epsilon_1>0\).
This use of (FL1) is on the shell of height \(Y\), so replace
each cap \(u\) in base \(n\) by \(u\log X/\log Y\) in base
\(Y^L\). These factors are bounded by \(1/\omega\).
The inner exponent losses in (FL1) can be made smaller than
\(\epsilon_1 L\) in base \(X\), uniformly in the pools of
bounded polynomial height in \(Y\). Choose \(\tau,\kappa\)
and \(\epsilon_1\) after \(\xi,m\), so that the bound in
(TM2) is power smaller than (TM1).
[Pointwise prescribed-pool count (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

Outside the pool, chosen prime powers from different rows are
coprime, and their congruences in \(v\) have unit slopes. If
records on the surviving \(C_i(t)\) have total normalized raw
cost at most \(p_*<r\) and combined numerical weight greater than
\(h_B-\beta\) by a fixed margin, then the reconstruction argument
in (KP3) gives at most \(n^{p_*+o(1)}\) surviving columns.
Indeed their product has height at most \(X^m\); enumerate
its radical and then the powers and ordered factors at subpower
cost. The fixed congruences put \(v\) in one residue class
modulo a product larger than the interval length. This contradicts
(TM1) at small enough fixed tolerances. The same argument proves
the marked and filler packet tests with cost ceiling \(r\)
and numerical threshold \(h_B-\beta\). The loss in (TM1) records
the entire dependence on the two dictionary deficits.
[Prime packing and uniform record enumeration](ABC_MINIMUM_DICTIONARY_COMPACT_PRIME_PACKING.md#several-large-primes-with-smaller-prime-contributions)

With arbitrary row and column caps \(n^{a}\) and \(n^{b}\),
and edge exponent at least \(d\), the same star count gives
common-column exponent at least
\[
                      s_m=m d-ma-(m-1)b
\]
when \(d>b\). For individual raw caps \(\ell_A,\ell_B,\ell_C\),
(FL1) strips target pool overlap of mass \(\xi\) when
\(\min\{\ell_B,\ell_C-\xi\}/2<s_m\) by a fixed margin.
When the outside cap is negative the corresponding class
is empty. Any pack that leaves enough primes after this strip,
has total record cost less than \(s_m\), and exceeds the
parameter width gives the same contradiction. This extension
is a sector test; it does not assert that a suitable pack
survives a large loss in the common-column exponent.

Here are support consequences for a limiting profile on \(C\).
Suppose \(h_C\) is its limiting numerical height, \(H>0\)
is the available exponent of the parameter interval, and the
strict inequalities just described can be tested with
threshold \(H+\zeta\) for arbitrary fixed \(\zeta>0\).
Then
\[
 \begin{split}
 \operatorname{supp}\nu_{C,\rm na}
     &\subseteq [0,H/r],\\
 \operatorname{supp}\nu_C
     &\subseteq
       [0,\min\{(h_C+H)/r,\ \max\{h_C/r,2H/r\}\}].
 \end{split}                                                    \tag{TM3}
\]
For the nonatomic bound, on a compact band of positive
nonatomic mass above \(s>H/r\), remove neighborhoods of its
finitely many atoms of mass above a sufficiently small
threshold. On the remaining compact continuity intervals
the individual prime sizes are uniformly small in normalized
cost, while positive total mass remains. This follows from
weak convergence, since a sequence of individual primes of
cost at least the threshold would give an atom in the chosen
compact set. Take enough rows to accumulate raw mass a little
above \(H/s\). Greedy selection overshoots this cost by at
most the individual prime threshold. By (TM2) the pool can
be stripped with arbitrarily small cost loss per row,
preserving the positive mass and the strict CRT margin.

For the first finite cutoff, use the same greedy selection
on any compact band with lower endpoint \(s\). Each selected
prime has cost at most \((h_C+o(1))/s\), because the full
prime power divides \(C\). Thus a positive-mass band with
\((H+h_C)/s<r\) is excluded. To obtain the second cutoff,
take \(s>h_C/r\), so each individual prime in the band
has raw cost less than \(r\) by a fixed margin. If any
of these powers has numerical weight larger than
\(H+\eta\) for fixed \(\eta>0\), test that power alone
with one row. Its cost is less than \(r\), and its mass
is bounded away from zero by the finite upper endpoint
of the band. Hence it survives a sufficiently small
pool omission in (TM2) and reconstructs the column.
This sector is power smaller. The remaining prime
costs in the band are at most \((H+\eta)/s\).
Use these instead for the greedy overshoot and let \(\eta\)
be small. The band is excluded when \(s>2H/r\) as well.
Continuity windows yield (TM3). Only compact bands are
tested in this argument.

## Arbitrary positive fixed-fibre exponent and escaping moments

Start with a fixed complete-letter fibre having
\(n^{r+o(1)}\) distinct free complete values, with a fixed
bound on its individual normalized radical caps.
Take an arbitrary relation of \(n^{3r-o(1)}\) increasing
triangles in that fibre whose primitive difference triples
have total radical at most \(n^{6r+o(1)}\).
Discard short normalized heights and small side radicals
using the fixed-fibre triangle counts. At fixed margins,
the remaining output triangles have height at least
\(X^\omega\) for a fixed \(\omega>0\), and the three
individual radicals all have exponent \(2r+o(1)\) in base \(n\).
The bounds on the sides follow from testing each at a fibre
threshold strictly smaller than \(r\) and then using the total
radical cap. Perform these discards before the exact-gcd and
pinning selections, so that the later count losses are smaller
than the fixed discard gaps.
[Short normalized heights and side radicals](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#primitive-triangles-of-three-values-on-the-fibre)

Refine the original triangle relation into height and valuation
boxes before pinning. At fixed \(L\), dyadic choices for a
fixed exact degree truncation cost height subpower. High-order
tails have radical mass at most the reciprocal of their scaled
cutoff. The usual height-first diagonal refinement gives
limiting numerical heights \(h_i\), weak degree laws \(\nu_i\)
and escaping moments
\[
                       e_i=h_i-\int s\,d\nu_i(s)
\]
for the three output side positions. A single such box
retains the full exponent. The following deductions concern
any selected limiting profiles in these boxes.

Fix a side \(i\) to be tested. In the exact-gcd selection,
pin the vertex opposite side \(i\) rather than either endpoint
of that side. The two incident side dictionaries each have
size at most \(n^{r+o(1)}\): the fixed vertex and exact scaling
gcd invert the map from an original endpoint to its incident
quotient side. Once the vertex and scaling have been fixed,
the two incident quotient values also determine both other
original endpoints. Thus the pinned primitive output has
\(n^{2r-o(1)}\) distinct triples. A middle vertex is an
allowed pin. In that case its two incident sides enter the
reconstruction with opposite signs and the target is the
unpinned side. A different choice of pin can be made for
each tested side, preserving the profiles on the original
selected box.
[Exact scaling and pinning in a single fixed fibre](ABC_HEREDITARY_MINIMUM_TRIANGLE_TRANSFER.md#a-single-fibre-at-an-arbitrary-positive-exponent)

Fix a large common continuity cutoff \(R\) for the output
measures. For each role \(j\), let \(F_j\) be the full factor
on scaled orders at least \(R\). Freeze all three \(F_j\)
at total normalized cost at most \(3/R\) plus inner allowance:
the original output heights are at most \(X\), so each tail
has radical at most \(n^{1/R}\). The fixed factors are pairwise
coprime. Since side \(i\) is zero modulo \(F_i\), the
residue classes on the two pinned roles are bijectively paired
by their signed additive equation. If one of such a pair
has fewer than \(n^{r-\eta}\) labels, all corresponding
edges together number at most \(2n^{2r-\eta+o(1)}\):
an edge is determined by the two pinned labels, and summing
the sizes of the opposite paired classes costs at most its
whole dictionary. Choose \(\eta\) larger than \(3/R\) plus the preexisting
count and dictionary deficits and inner losses. The other
class pairs number at most \(n^{\eta+o(1)}\), so one can
fix a pair at this additional cost. Each pinned role \(j\)
now lies in a fixed residue class modulo \(F_i F_j\),
since it is already divisible by \(F_j\).

For a fixed cutoff, the weight of \(F_j\) tends to
\(h_j-\int_{[0,R)}s\,d\nu_j(s)\). As \(R\to\infty\) this
decreases to \(e_j\). Take \(R\) large and the initial deficits
small relative to all desired fixed counting gaps; in the
outer diagonal take \(R\to\infty\) and then shrink the residue
losses along the budgets. Each cutoff and every residue
threshold is fixed before its budget's height limit.
Using pinned role \(j\) as the column therefore gives
parameter exponent arbitrarily close from above to
\(h_j-e_i-e_j\). Put
\[
                H_i=\min_{j\ne i}(h_j-e_j)-e_i .
\]
There is a constant \(\theta_{r,\omega}>0\) such that
\[
                             H_i\ge\theta_{r,\omega}.           \tag{TM4}
\]
Otherwise after the tail and residue selections, fix the
other pinned role by averaging. The resulting primitive
output fibre still has exponent close to \(r\), while the
column lies in an affine progression with integer parameters
of size at most \(X^{H_i+o(1)}\). If this exponent is
negative, the fibre has at most one value. If it is sufficiently
small positive, the fixed-fibre rank-one count gives a
contradiction. For that count use the output height shell
\(Y\) and a fixed fibre threshold \(r/2\) in base \(Y^L\);
the conversion between log heights is bounded by \(1/\omega\).
The feature threshold can thus be chosen absolute for the
fixed parameters \(r,\omega\).
[Small features in a fixed complete fibre](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#a-small-feature-bound-on-a-single-fibre)

Apply (TM3) with \(H=H_i\), using the better of the two
pinned columns. This yields, for each of the three side
positions of the original selected triangle profile,
\[
 \begin{split}
 \operatorname{supp}\nu_{i,\rm na}
     &\subseteq[0,H_i/r],\\
 \operatorname{supp}\nu_i
     &\subseteq
        [0,\min\{(h_i+H_i)/r,\ \max\{h_i/r,2H_i/r\}\}].
 \end{split}                                                    \tag{TM5}
\]
In particular for a wholly nonatomic side its total mass
is \(2r\). Its finite first moment \(h_i-e_i\) is then at
most \(2H_i\), so
\[
                  h_i+e_i\le2\min_{j\ne i}(h_j-e_j).              \tag{TM6}
\]
These statements do not assume that all three output
dictionaries are minimum in a single pinning. They apply
also at any positive intermediate multiplicity on the
original fibre. Neither compact support nor (TM4) bounds
the full raw radical of the difference triangles.

## Single-letter numerical strips and full residues

The variable-mesh grouping argument also applies to one
set of distinct complete values, rather than a ratio of
two letters. Here is a formulation in base \(n=X^L\).
Suppose the values are positive, at most \(X\), with
radical at most \(n^C\) for fixed \(C\). In a height-first
profile let \(\mu\) be the degree measure, and let \(c_\mu\)
be its nonatomic mass. Assume \(c_\mu\ge c_*>0\) and let
\(d\) be the normalized count exponent of the set.

There is a constant \(k_*>0\), depending only on \(c_*,C\)
and fixed height bounds for any modulus, with the following
properties.

- If all values \(W\) are in a relative real strip
  \(\lvert W/K-1\rvert\ll X^{-\tau}\) around one fixed
  positive \(K\), with fixed \(0<\tau\le1\), then
  \(d\le c_\mu-k_*\tau\).
- Suppose the values are all in a prescribed unit residue
  modulo \(M\le X^B\) with limiting numerical weight
  \(w=\lim\log_X M>0\). There is a fixed
  \(R_*=R_*(c_*,C,B)>0\) such that if \(M\) is
  \(\lfloor R_*/L\rfloor\)-full, then
  \(d\le c_\mu-k_*\min(w,1)\).

The strip centre, unit residue, and modulus may move with
height. The assertions concern the specified limiting
profile, with the same height-first convention as (VM1).
[Variable mesh grouping (VM3)](ABC_VARIABLE_MESH_DEGREE_GROUPING.md#one-grouped-root-on-the-free-pair)
· [Real strip count (RUW1)](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-real-strip-with-unequal-root-heights)
· [Unequal-height congruence count (UW)](ABC_QUANTITATIVE_DEGREE_WINDOWS.md#a-product-bound-for-the-admissible-wider-height)

Here are the details of the grouping. Choose
\(0<\vartheta\le1/8\), a compact positive continuity
interval \(J\) with nonatomic mass at least \(c_*/2\),
and \(m=\lceil 4/\vartheta\rceil\). Put \(t_0\) equal
to the nonatomic mass in \(J\) divided by \(m\).
Then \(t_0\ge c_*\vartheta/10\). Choose finitely many
positive atoms leaving omitted atomic mass at most
\(t_0/16\). Simultaneously approximate their scales
\(s_a\) on a mesh of spacing \(a\), with
\(\|s_a/a\|\le\vartheta t_0/(16m(C+1))\).
A sufficiently small such \(a>0\) exists by simultaneous
Dirichlet approximation, as in (VM3). Take \(a\) below one
and below half the lower endpoint of \(J\) and half every
selected positive atomic scale. For every scale
\(s\in J\) some integer \(1\le j\le m\) satisfies
\(\|js/a\|\le 1/m\), by comparing the \(m+1\) consecutive
multiples of \(s/a\) on the unit circle. Thus one mesh
of spacing \(b=a/j\) captures at least \(t_0\) nonatomic
mass within normalized distance \(\vartheta/2\) of a
mesh point. Enlarge the window radius to a continuity
radius at most \(\vartheta\), and denote its captured
nonatomic mass by \(t\ge t_0\).

Keep those windows, all the selected atomic windows
shrinking to their atoms, and a shrinking zero window.
Freeze all other full prime powers at normalized cost
at most \(c_\mu-t+t_0/16+o(1)\). This includes any high
tail not retained, since the frozen block is counted by
its radical and fixed-support enumeration. In the kept
positions assign each positive scaled order to the
appropriate integer multiple of \(b\), and use
\(q=\lfloor b/L\rfloor\) as its representative degree.
After removing the fixed frozen factor, the exact
value is \(W_{\mathrm w}Z^q\) with
\[
 H(Z)\le X^{(O_C(1/b)+o(1))L},\qquad
 H(W_{\mathrm w})\le X^{2\vartheta b t+o(1)} .
\]
The first estimate uses the first-moment bound of one
and the total raw cap \(C\). For the second, the captured
nonatomic positions have discrepancy at most
\(\vartheta b\) per unit of mass, the finitely many
selected atoms contribute total at most
\(\vartheta b t_0/16+o(1)\), and the other kept atoms
contribute at most \(\vartheta b t_0/16+o(1)\).
The zero window and the rounding of \(q\) have vanishing
height cost. Once frozen factors are fixed one can take
one presentation for each exact value, so the exact-value
fibre bound in the unamplified numerical estimates is one.
No coprimality between \(Z\) and \(W_{\mathrm w}\) is needed.

For the real strip take \(\vartheta\) a sufficiently
small fixed multiple of \(\tau\). In (RUW1) choose the
one-variable chart parameter \(t/8\). The narrow height
product there is \(O_C(1/(bt))\), so its wider-height
hypothesis allows \(2\vartheta bt\). The chart cost is
\(t/4\) in base \(n\). For the congruence use
\(\vartheta\) a sufficiently small multiple of
\(\min(w,1)\), and take the fullness parameter with
\(kL\to8w/t\). This is available by the stated fixed
threshold \(R_*\), since \(t\ge c_*\vartheta/10\)
and \(w\le B\). In (UW) the narrow height constant
relative to \(M\) is \(O_C(1/(bt))\), and the wider
relative numerical exponent is at most
\(2\vartheta bt/w+o(1)\). Their product is sufficiently
small by the choice of \(\vartheta\). The chart cost
is again \(2w/(8w/t)=t/4\). Dividing by the frozen
factor changes only the fixed centre or unit residue.
In both cases the total cost is at most
\(c_\mu-11t/16+o(1)\).
This proves the two assertions with fixed choices
at the limiting profile before the inner height limits.

## Consequences for the affine images and the scaling gcd

A fixed original complete-letter fibre
of exponent \(d\ge d_0>0\) has nonatomic mass at least
\(d\) on each moving letter by (FT1). For a
fixed first endpoint, let the positive complete original
values on its two moving roles be \(a,b\). Its other
endpoint is represented on these roles by
\[
                        a+Gz,\qquad b\pm Gz
\]
in a suitable signed orientation. If the numerical height of
\(Gz\) is smaller than that of the corresponding anchor
by \(\tau>0\), the image lies in a relative strip of
precision \(X^{-\tau+o(1)}\). Each image determines the
original row on the fixed fibre. The real-strip conclusion
therefore charges that positive height gap against the
nonatomic excess of the image letter above \(d\).

In particular take an original minimum complete-dictionary
endpoint of minimum nonatomic total, so each letter has
nonatomic mass \(c_0=1/6\). In any fixed-role fibre pair
one of the two moving letters is at least \(X/4\) on the
first endpoint. If the gap in either moving coordinate is
at most \(X^{1-\tau}\), the corresponding other endpoint
on that tall role is in a relative strip of precision
\(O(X^{-\tau})\). Fix the fibre and first endpoint, and
apply the single-letter strip count to this role. If its
count retained exponent \(c_0\) in the original profiled
boxes, the estimate would contradict \(c_\mu=c_0\).
The estimate is uniform in the strip centre. Summing over
the at most \(n^{2c_0+o(1)}\) fixed choices gives a
power-small class of short fibre pairs at the
\(n^{3c_0}\) scale in these profiles. Each pair extends
to at most \(n^{c_0+o(1)}\) tagged triangles.
Consequently all three physical sides of any full tagged
triangle relation have limiting numerical height one.
After an exact common scaling gcd \(G\) is fixed, all
three corresponding normalized side heights are
\(1-\lim\log_XG\). This conclusion uses the original
arithmetic caps even when no radical cap on the
difference triangles is assumed.
[Nonatomic mass on a fixed fibre (FT1)](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#nonatomic-mass-outside-prescribed-prime-pools)
· [Exact tagged scaling selection](ABC_HEREDITARY_MINIMUM_TRIANGLE_TRANSFER.md#fixing-the-exact-scaling-and-the-three-residue-classes)

There is a parallel constraint on the exact scaling gcd.
On a fixed fibre of exponent \(d\ge d_0\), suppose
both free roles lie in prescribed residues modulo the
same integer \(G\le X\). For each prime in the
\(\lfloor R_*/L\rfloor\)-full tail of \(G\), at least
one of the two free roles is a unit, by original
primitivity. Assign the prime power to one such role.
This is a choice from the two fixed residues, independent
of the member of the fibre. The two resulting full
moduli are coprime, their product is the full tail, and
each gives a unit residue for its assigned role.
Write \(c_U,c_V\) for the two nonatomic masses.
If one of the resulting moduli has limiting weight zero,
use the unconditional \(d\le c_W\) on its role.
Otherwise use the unit-residue conclusion on that role.
Adding gives
\[
  \limsup\ \log_X
       \prod_{\substack{p^e\parallel G\\eL\ge R_*}}p^e
                 \ \ll_{d_0,C}\ c_U+c_V-2d,                    \tag{TM7}
\]
with \(R_*\) enlarged by a fixed margin if needed for
integer endpoints. Every limit can be taken on a
subsequence attaining the displayed upper limit.
At minimum nonatomic original saturation and full fibre
exponent \(d=c_0\), this high-order numerical tail of
the exact scaling gcd vanishes. The assertion imposes
no bound on \(\operatorname{rad}(G)\). It leaves the
numerical weight of moderate orders in \(G\) unrestricted.
