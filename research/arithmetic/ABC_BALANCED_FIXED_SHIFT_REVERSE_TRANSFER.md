# One-sided powerful-number bounds on balanced abc residue cells

A fixed shift can save when just one variable term nearly
exhausts its numerical height at an active fullness floor.
The factors below the floor are charged by their actual
joint dictionary in the fixed-shift context. This applies
to ordinary complete-letter fibres and to the off-diagonal
codegrees forced in balanced residue cells.

## Fixed-shift contexts and conditional charges

Put \(n=X^L\) and use the height-first conventions of
(RB1)–(RB8). The original rows are positive and primitive,
\(A+B=C\), with \(X/2<C\le X\). Keep all three radical caps,
and any eligibility restriction on the original rows such
as a least-repair cutoff. The discussion uses fixed finite
meshes for normalized exponents and fixed anatomic lists
at every fixed budget. All small-budget margins are chosen
before the height limit. Normalized count allowances can
be made arbitrarily small at sufficiently small fixed \(L\).
[Balanced cohorts and residue cells](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md)

The argument uses two types of context. In the first,
fix one original complete letter \(T\) and count its two
free letters \(U,V\) satisfying a nonzero signed shift.
In the second, take an actual residue cell from a fixed
lock \(D\) on \(W\). Fix two different complete values
\(V_1,V_2\) on one free role and count common complete values
\(U\) on the other. The two computed letters satisfy
\[
 W_i=sU+tV_i,\quad W_1-W_2=t(V_1-V_2),\qquad s,t\in\{-1,1\},
\]
with both original rows retained. Indeed, all of these
shifts and labels come from the same original additive law.
In the codegree context the common \(D\mid W_i\) divides
their difference. Remove the parts of both computed
letters supported on that difference, after fixing those
parts at height-subpower cost. The remaining factors
have the original exact prime valuations and the shift
has fixed positive integer coefficients and coprime
complete terms. For general positive polynomial-height
column weights the same argument applies to the
weighted difference (FS2): include the column weights
among the fixed supported data, enumerate the gcd at
divisor cost, and divide it out. Condition on that
record too. Fixing supported parts does not discard
any original radical restriction or increase a free
factor's height.
[Weighted shifts and fixed-support stripping](ABC_FULL_SHIFT_DICTIONARIES_AND_CODEGREES.md#weighted-differences-in-a-varying-shift-graph)

A **unit shift modulus** below means a full divisor of the
*normalized* shift with a fixed positive numerical weight
and minimum scaled valuation bounded below by a positive
constant. It is coprime to the two normalized complete
terms. The \(W\)-lock on both codegree targets is removed
in the normalization and supplies no such modulus.
For a fixed original \(T\) of positive numerical height,
primitivity makes its full high-degree part a unit shift
modulus. To see that this part has positive weight,
the degrees below \(a/L\) have weight at most
\(a\log_n\operatorname{rad}(T)\); take fixed \(a>0\) small
compared to the lower height bound. A fixed relative real
height gap, as in (FS1), can be used in place of a unit
shift modulus when it exists for the normalized terms.
[Modular and real fixed-shift hypotheses](ABC_FULL_SHIFT_DICTIONARIES_AND_CODEGREES.md#a-fixed-shift-count-without-a-radical-cap-on-the-shift)

In a fixed context prescribe one or both deficient-factor
records. Each is an exact integer product on selected
original prime-support positions, with the slot list
available at subpower cost at fixed truncation. The
records may be joint and may include any additional
full factors that are worth fixing. Let \(b\) be a bound
on their conditional joint dictionary exponent. For the
remaining low positions below a fixed positive scaled
cutoff \(a\), let \(s_a\) be the corresponding joint
dictionary exponent, after the first record is fixed.
One may always use the raw sum of remaining low radical
exponents for \(s_a\). All root or bounded-group dictionary
caps in the cited estimates are the actual conditional
caps in the resulting fibres. The arbitrary-shift
(FS5) uses singleton or fixed bounded-size groups.
The full-modulus and real-gap (FS1) also permits
the specified shrinking-window groups. A bound on one complete
dictionary is not a bound on the sum of its singleton
position caps.

There is also an average-cost formulation. If \(J\) is
the selected row or common-neighbour variable, \(K\)
the fixed-shift context and \(F\) its deficient record,
then
\[
 H(J\mid K)=H(F\mid K)+H(J\mid K,F).
\]
Use the fibre size bound on the last term and take its
conditional average. This charges the actual joint
entropy without asserting that either marginal is
independent. If the residual hypotheses vary, the
formula uses their bounds on the stipulated bins.
The cost of supported-part stripping is subpower and
can be absorbed separately.

## An optimized one-sided cutoff bound

Choose one restrictive variable part \(V\). Fix its
positive deficient factors below scaled order \(\kappa>0\),
with the option to leave a shrinking window at zero
of radical exponent at most \(z\) and numerical height
\(X^{o(1)}\). The remaining high part \(V_{\rm hi}\) has
floor \(\kappa/L+O(1)\), a numerical height cap
\(H=\limsup\log_X V_{\rm hi}\), and weight excess
\[
 E=\limsup\left(\log_X V_{\rm hi}
              -\frac{\lfloor\kappa/L\rfloor}{\log X}
                         \log\operatorname{rad}(V_{\rm hi})\right)\ge0.
\]
Keep both the remaining raw mass and numerical height
after the fixed records and support stripping. Write
\(W_U,W_V\) for any allowed conditional summed dictionary
caps in base \(n\) and set
\[
       m=\frac{\min\{W_U,W_V,H/\kappa+z\}}2,\qquad
       m_0=\frac{H/\kappa+z}2 .
\]
Use \(a<\kappa\) for the low-position cost \(s_a\). The
opposite side has no prescribed valuation floor.

**RT1.** The normalized fixed-context count exponent is
at most
\[
 b+\max\{m,s_a\}+o(1)
\]
for arbitrary shifts, or \(b+m+o(1)\) with a unit shift
modulus or the stated real gap. If \(m_0\) is bounded
away from zero, there is a fixed positive saving from
\(b+m_0\) when \(s_a<m_0\) by a fixed margin, or with the
unit modulus or real gap and no condition on \(s_a\).
The strict saving is uniform as the finite parameters
stay in compact admissible ranges, including the lower
bounds on \(\kappa\) and, where used, the terminal margin.
In particular, one can use the minimum over the two
orientations and over any fixed finite collection of
cutoffs and record strategies. The cap and cost for each
strategy must be measured together in its actual
conditional fibre.

For the first bounds apply (FS5) or (FS1) after the
record is fixed. In (FS5), the cost at the terminal step
can be the joint dictionary of the low parts. Its proof
uses those parts only when fixing the low remainder
after numerical concentration. Fixing all low parts,
instead of just the remainder, costs no more than
the stated joint bound. The high rearrangement and
concentration steps use the allowed group costs.
The upper \(H/\kappa+z\) follows from fullness and height.
[Fixed-shift estimates (FS1 and FS5)](ABC_FULL_SHIFT_DICTIONARIES_AND_CODEGREES.md#arbitrary-shifts-with-high-free-degrees)

For the strict assertion, suppose there is a sequence
where the residual count tends to \(m_0>0\). The applicable
half bound forces the restrictive letter's remaining
raw mass to tend to \(H/\kappa+z\). In particular its high
mass tends to \(H/\kappa\) and its high numerical height
tends to \(H\). Write the high part as
\[
 V_{\rm hi}=(\operatorname{rad}V_{\rm hi})^{\lfloor\kappa/L\rfloor}
                       \prod_{p^e\parallel V_{\rm hi}}
                                  p^{\,e-\lfloor\kappa/L\rfloor}.
\]
The correction has numerical height \(X^{o(1)}\);
the retained zero window has that same height property.
This gives the compression in (SHF2). Its proof counts
the low terminal by at most \(s_a\) by the same joint
fixing argument as in (FS5). This contradicts a fixed
margin \(m_0-s_a>0\). A restrictive side with zero
limiting high height is already a wide cofactor and
obeys the same compression argument. With a unit
shift modulus or real gap, the first switch in (SHF2)
makes the other side compressible too, by the cheap
full-block dichotomy. Then release the two compressed
values by the terminal step of (FS1). This costs
\(o(1)\) in base \(n\), because their exact coprime ratio,
or the recovered single value in the short-free-letter
case, determines the pair. This proves strictness in
that case as well. If no uniform positive margin and
small budget threshold existed, take a sequence of
violations with the lost margin and budget tending to
zero. Take the height limit first at each budget and
apply the preceding argument to the limiting parameters.
[Compressed fixed shifts (SHF2)](ABC_FIXED_SHIFT_STRICT_HALF_SAVING.md#a-compressed-term-in-an-arbitrary-fixed-shift-fibre)
· [Unit-modulus and real-strip terminals (FS1)](ABC_FULL_SHIFT_DICTIONARIES_AND_CODEGREES.md#proof-by-concentration-of-the-free-roots)

**RT2.** More generally, suppose after the records the
restrictive letter compresses into a fixed number of
positive scaled-order windows and a vanishing-height
correction, as in (SHF2). The normalized bound becomes
\(b+s_a+o(1)\) for arbitrary shifts, and \(b+o(1)\)
under the unit modulus or real gap. In particular
\(E\to0\) in the one-sided cutoff is sufficient.
For each fixed normalized allowance \(\eta>0\),
fixed number of windows and fixed lower scaled
cutoff, these bounds hold with \(2\eta\) in place
of \(o(1)\) when the numerical correction is at
most \(X^\omega\) for some \(\omega>0\) and
at all sufficiently small fixed budgets. This follows
by taking a violating sequence with both width and
budget tending to zero, and applying the vanishing-width
argument. It gives a criterion without using a numeric
value of the SHF1 gap. In an original fixed-\(T\) fibre
with positive numerical height, a one-sided compressed
record of cost \(b<u\) with fixed margin therefore
excludes that sector. In weighted codegrees one must
include the low-position cost unless an independent
terminal is supplied.
[Fixed-width compressed-core estimate](ABC_COMPRESSED_CORE_WEIGHTED_CODEGREES.md#setting-and-finite-width-hypothesis)

## From local savings to balanced graph savings

**RT3.** In a regular balanced near-saturating cohort
let \(\alpha\) be the common complete dictionary
exponent, \(u=1/3-\alpha\), and \(s=1/3-2u\). In an
ordinary fixed-complete fibre a bound strictly below
\(u\) with fixed margin contradicts regularity. On a
cheap fixed-\(W\) chart, each typical cell has free
dictionaries of size at most \(n^{u+o(1)}\), and there
are at most \(n^{s+o(1)}\) typical cells. Suppose in
each retained cell one orientation has every distinct
column-pair codegree bounded by \(n^{q+o(1)}\), with
the same \(q<u\). If the number of rows in a cell is
\(e\) and the row and column dictionary sizes are
\(x,y\), respectively, Cauchy–Schwarz gives
\[
                       e^2/x\le e+y^2n^{q+o(1)} .
\]
Thus the sum of the rows on the chart has exponent at
most \(s+\max(u,(3u+q)/2)\), before the chart, pruning
and count allowances. This is below \(1/3\) by
\(\min(u,(u-q)/2)\). For near-saturating bins choose the chart loss and
RB4 cell tolerances smaller than the proposed power
margin. For an endpoint contradiction it suffices
to apply the bound to one typical exponent-full cell,
whose existence follows from (RB4). The average-cost
version uses the uniform law on oriented two-edge
paths in a typical cell, conditioned on its two
distinct column anchors. At \(u>0\) this law has
conditional entropy at least \((u-o(1))\log n\),
because the number of paths is at least
\(e^2/x-e=n^{3u-o(1)}\) and the anchor dictionary
has size at most \(y^2\). This proves the same
contradiction if the average of the RT1 or RT2 prices,
including the conditional record cost, is below \(u\)
by a fixed margin. No rectangle is required to leave
the original arithmetic graph.
[Actual residue cell sizes (RB4)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#the-actual-residue-cell-sizes)

For \(u=0\) the diagonal bound already reaches the
local exponent; there may be one row per cell.
Even a subpower bound for every nondegenerate codegree
would then give only the existing global exponent.
Uniform closure of the full range requires an
additional reduction in the count of occupied cells,
or in their record prices, as \(u\) approaches zero.
An existing primitive matching construction has
two original radical caps and simultaneous affordable
record saturation, with the third radical unrestricted.
It supplies no power-many local codegrees.
[Two-cap matching construction](ABC_OVERLAPPING_FULL_RECORD_COUNTEREXAMPLES.md#fixed-parameters-and-complete-integers)

## Codegree thresholds at the actual chart price

**RT9.** Let a fixed finite-price \(W\)-lock chart of an
original regular balanced cohort retain row exponent at
least \(d-f-o(1)\), where \(d=1/3\).
Here \(f\) is the row exponent lost to the chart and
the original complete dictionary exponent is
\(\alpha=d-u\). Take a cell bin with \(n^{z+o(1)}\)
occurring unit residues, \(n^{m+o(1)}\) rows per cell
and free dictionary exponents \(v\ge w\), after swapping
the orientation if necessary. A bin retaining the chart
exponent satisfies
\[
 z+v\le\alpha+o(1),\qquad z+2v\le d+o(1),\qquad
 m\ge\max(u+v,2v)-f-o(1).
\]
The first inequality uses disjoint complete dictionaries
in distinct cells. For the second apply RB2 to the
larger dictionaries, on the same free role throughout
the bin, with the list of its \(n^{z+o(1)}\) occurring
unit classes. The last follows from \(z+m\ge d-f-o(1)\).
An arbitrarily small fixed exponent mesh costs only a
fixed number of bins before the height limit. The errors
in these displays include that mesh and all chart
regularity and residue-count allowances.

If \(f<u\) by a fixed margin, then \(m>v\) by a fixed
margin. Put the distinct anchors on the smaller side.
The uniform law on the resulting two-edge paths,
conditioned on their anchors and the cell, has entropy
in base \(n\) at least
\[
 2m-v-2w-o(1)\ \ge\
 \max(2u-v,v)-2f-o(1)\ \ge\ u-2f-o(1).
\]
Indeed, each cell supplies at least
\(e^2/x-e\) paths, with \(e\) its number of rows and
\(x\) the larger dictionary size. Divide the path
count by the number of possible ordered anchor pairs
in the entropy bound. Any average fixed-shift price
that falls below \(u-2f\) by a fixed margin therefore
excludes a full-exponent bin. Apply the original-fibre
lift on these paths when the computed shift has lost
its fixed full divisor. All shifted counts are taken
on the actual original lifts.

For a uniform codegree upper exponent \(q\), there is
an equivalent direct count in every bin. The path bound
gives \(m\le\max(v,(v+2w+q)/2)+o(1)\), including the
diagonal term. From the first two inequalities,
\(z+3v/2\le(\alpha+d)/2+o(1)\). The chart row exponent
is consequently at most
\[
              d-\min\{u,(u-q)/2\}+o(1).
\]
Thus \(q+2f<u\) by a fixed margin suffices after
accounting for chart price. The entropy threshold, capped below at zero, is
sharp for the abstract size constraints. For \(0\le f<u\)
use \(n^{d-2u}\) copies of one bipartite graph on two
sets of size \(n^u\) with independent edge probability
\(n^{-f}\). Rounding the sizes affects no exponent.
The degrees concentrate around \(n^{u-f}\).
A codegree is binomial with mean \(n^{u-2f}\), so
for any fixed \(\eta>0\), exponential binomial tails
and the union bound give maximum codegree at most
\(n^{\max(0,u-2f)+\eta}\) with probability tending to
one. This attains the entropy threshold as a graph
model. A realization with all three radical caps
would require additional arithmetic.
The optimized threshold still tends to zero as
\(\alpha\) tends to \(1/3\).
[Original-fibre codegree bounds](ABC_BALANCED_ORIGINAL_FIBRE_LIFT.md#counting-a-common-neighbour-in-its-original-fibre)
· [Optimized occupied-cell count](ABC_JOINT_REAL_RESIDUE_MOMENTS.md#optimized-list-prices-after-full-factor-records)

## Collision-weighted codegrees without a cell size mesh

**RT10.** On one actual chart, let the shared cell label
be a function of either complete free letter, with the
sign of a residue absorbed on one side. Write \(E\) for
the number of original rows, \(x_j,y_j\) for the two
dictionary sizes in cell \(j\), and \(D_U=\sum_jx_j\),
\(D_V=\sum_jy_j\). Keep the indicated original signed
additive relation and all three radical caps. With distinct
anchors on the \(V\) side and uniform codegree cap \(Q\),
\[
 E\le D_U+\left(Q D_U\sum_j y_j^2\right)^{1/2}.
\]
Indeed \(e_j^2/x_j\le e_j+Q y_j^2\) by counting ordered
paths through the \(U\)-vertices in cell \(j\).
Sum over \(j\) and apply Cauchy–Schwarz once more to
obtain \(E^2/D_U\le E+Q\sum_j y_j^2\), which implies
the bound. A variable cap \(Q_j\) uses
\(\sum_j Q_jy_j^2\) in place of \(Q\sum_jy_j^2\).
The \(U,V\) dictionaries are counted once per complete
value because the shared label is a function of either
value.

For a deep unit lock, let \(r_V\) be the actual remaining
raw cap **after** fixing full factors of \(V\). Multiplying
the residue labels by the inverse of the fixed part gives
unit classes for the quotient. Apply RB2 to the \(j\)
most populated classes. The \(j\)-th largest dictionary
then has size at most \(n^{r_V/2+o(1)}j^{-1/2}\), so
\(\sum_j y_j^2\le n^{r_V+o(1)}\) with a height-subpower
harmonic factor. Summing chart states of exponent \(b\)
therefore gives the row exponent bound
\[
         b+\max\{\alpha_U,(\alpha_U+r_V+q)/2\}+o(1),
\]
where \(\alpha_U\) caps \(\log_n D_U\) and \(q\) caps
the actual fixed-shift codegree exponent. One can minimize
over either orientation and over available record strategies,
using the residual caps and the fixed numerical heights
from the same record. For one chart retaining exponent
\(d-b-o(1)\), with \(d-b>\alpha_U\) by a fixed margin,
the path comparison has threshold \(2d-2b-\alpha_U-r_V\).
In particular a record capturing mass \(f_V\) on the
anchor role can raise the balanced threshold from
\(u-2b\) to \(u+f_V-2b\); its actual price \(b\) remains
in the comparison.

For any shared label, let \(h_2\) be its normalized
collision entropy for the uniform law on the
\(D_V\) distinct complete values. Then
\(\sum_j y_j^2=D_V^2n^{-h_2}\) exactly. This permits
the optimized JRM4 bound on \(h_2\) when the shared
label is an actual shallow residue or an actual
shared combination of a congruence and a real bin.
For other mixed labels the assumption that either
complete value determines the joint cell must be
checked before using this estimate. Here all heights
and valuation caps in JRM4 pertain to the residual
root product in its own full-factor state. A fixed
coefficient only translates the logarithmic bins and
permutes the residue classes.

There is an average fixed-shift version. The actual
ordered distinct-anchor path population has size
at least \(E^2/D_U-E\). Its conditioning alphabet
consisting of the cell and the ordered anchors has
size at most \(\sum_j y_j^2\). When \(E/D_U\) exceeds
a fixed power of \(n\), the uniform path law has
conditional entropy in base \(n\) at least
\[
 2\log_n E-\log_n D_U-2\log_n D_V+h_2-o(1).
\]
Compare this to the conditional record entropy
plus the average certified residual fixed-shift count.
Using the original-fibre lift for each path retains
the third radical cap and the original numerical
height, including when a common full divisor is removed
in the normalized difference. The weighted estimate
does not posit independence of the residual factors.
[Optimized conditional collision entropy (JRM4)](ABC_JOINT_REAL_RESIDUE_MOMENTS.md#collision-and-pointwise-versions)
· [Original-fibre fixed-shift counting (LFR1)](ABC_BALANCED_ORIGINAL_FIBRE_LIFT.md#counting-a-common-neighbour-in-its-original-fibre)

### A necessary extra discount beyond raw capture

**RT11.** Suppose an exponent-full balanced cohort is covered
by joint row-level records of exponent \(b\) capturing
prescribed full radical masses \(f_V,f_Z\) on the anchor
role \(V\) and a distinct role \(Z\) used for the codegree
count. These masses are in arbitrarily narrow bins.
Assume the third role supplies a cheap lock, and use the
regular cohort on that lock before charging these records.
At positive \(u\) a cheap lock is available on each original
role. The two selected factor records are functions of the
joint record, with any fixed support choices and positional
data priced separately. LFR3 gives
\[
                         2b\ge f_V+f_Z-o(1).
\]
It applies to the two letters opposite the cheap lock
on the third role, even when both records depend on the
whole row. Its joint entropy lower bound compares with
the upper bound \(H(R)\le b+o(1)\) on the actual record
alphabet. The preliminary lock is chosen at arbitrarily
small price before the height and box allowances.

In the deep unit-lock application of RT10, the gain on
the anchor-cap side is \(f_V\). If the only upper bound
on the residual codegree is the raw half price
\((d-f_Z)/2\), the proposed comparison against
\(u+f_V-2b\) cannot produce a fixed saving. More
precisely, a certificate \(q=(d-f_Z)/2-\Delta\)
must have
\[
               \Delta> d/2-u+2b-f_V-f_Z/2
                      \ \ge\ d/2-u+f_Z/2-o(1)
\]
to win by that comparison. This obstruction concerns
full-factor captures paid as row-level chart records.
A factor removed using the particular two-anchor
difference, a better conditional dictionary on the
selected paths, and a list deficit in the actual
occupied cells are separate possible inputs. The
formula isolates the additional discount required
from those inputs or from numerical height and
compression on the original fibres.
[Two-record entropy charge (LFR3)](ABC_BALANCED_ORIGINAL_FIBRE_LIFT.md#full-factor-prices-in-the-shared-residue-register)

## A role switch that preserves the full shift divisor

**RT5.** On a typical \(W\)-locked residue cell, instead
use \(W\) as the common-row label, the same \(V\) as
the column label, and \(U\) as the computed target.
The anchors \(V_1,V_2\) still agree modulo \(D\) and
their nonzero signed difference is the shift. The two
computed \(U_i\) are units modulo \(D\), by original
primitivity. Hence their common gcd is coprime to \(D\)
and the normalized shift retains \(D\). The unit-modulus
branch of RT1 and RT2 applies with arbitrary low-degree
cost on the opposite side. If the \(W\)-dictionary in
the cell has exponent at most \(c_W\), this orientation
needs a codegree gap below \(2u-c_W\), since its edge
exponent is \(2u\) and the \(V\)-dictionary exponent is
\(u\). In particular the default bound \(c_W\le\alpha\)
makes this route nonvacuous when \(\alpha<2/9\).
A better \(W\)-dictionary bound in the actual cell
improves the threshold. With this orientation the
computed target's own fixed factors are stripped before
its variable height and fullness are charged.

There is an all-role version. Suppose all three cheap
locks have been fixed and write \(J=(J_A,J_B,J_C)\)
for the three residue registers in (RB6)–(RB7). Use the
uniform law on that chart, with entropies in units of
\(\log n\). Write \(s=2\alpha-1/3\) and \(z_*\) for
the common pairwise register mutual information.
Put \(v=\alpha-2s+z_*\), so
\(0\le v\le u\) and \(v\ge 2/3-3\alpha\). One has
\[
 H(\text{row}\mid J)=\frac13-3s+2z_*+o(1)=2v+o(1).
\]
For each \(T\), both registers opposite the other
two locks are functions of \(T\), so
\(H(T\mid J)\le\alpha-2s+z_*+o(1)=v+o(1)\).
Any pair of complete values determines the row.
Subadditivity gives equality for every \(T\) and
\(I(U;V\mid J)=o(1)\) for each complete pair.

**RT6.** If \(v>0\), one can extract an actual common
three-register cell at exponent \(2v'\) with \(v'\ge v\)
up to arbitrarily small losses and with all three
complete dictionary exponents at most \(v'\), where
\(v'\le u\).
First prune the atypical one-lock cells for all three
locks by (RB4), taking the combined chart price much
smaller than the fixed cell tolerance. This is a
power-small deletion and leaves all three complete
marginals at most \(n^{u+o(1)}\) in a common-register
cell, in the iterated limit. Then choose a cell with row exponent at least
\(2v-\tau\) and each conditional pair mutual information
tending to zero. Such cells exist for each fixed
\(\tau>0\): the mean row exponent is \(2v+o(1)\),
the exponent is bounded above, and the averaged
nonnegative mutual informations tend to zero. Within
the chosen cell use the uniform row law on \(e\) rows.
Let \(I_T=-\log_n P(T)\) denote the information of
the sampled complete value. For each complete pair,
\(Z=I_U+I_V-\log_n e\) has mean the pair mutual
information in the cell and satisfies
\(\Pr\{Z<-\xi\}\le n^{-\xi}\) for \(\xi>0\).
Indeed the pair uniquely specifies the row and the
product marginal probabilities sum to at most one
over the observed pairs. Integrating this lower-tail
bound gives \(\mathbb E(Z_-) \le 1/\log n\), and hence
the upper tail above every fixed \(\xi\) tends to zero.
The three pair sums show that on a positive fraction
of the rows each \(I_T\) is within a vanishing allowance
of \(\frac12\log_n e\), with the tolerances decreased
diagonally. The retained values in each dictionary
thus number at most \(n^{\frac12\log_n e+o(1)}\).
All rows retain the original caps. Passing to the
limiting cell size gives the claimed \(v'\).

In such a cell the RT5 role switch can use any lock,
with the codegree threshold now \(v'\). The three
full divisors are pairwise coprime by original
primitivity, so a fixed divisor belonging to a computed
target does not destroy the unit shift modulus from
the switched role. Positivity of \(v\) is automatic
on \(\alpha<2/9\) when all three roles supply cheap
locks. It is not a consequence of the register bounds
on the rest of the range, and an extraction of one
cell places no upper bound on the number of cells.
[Conditional entropy of the registers (RB5–RB7)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#conditional-entropy-identities-and-three-registers)

## Arbitrary weighted shifts at radical-scale height

**RT7.** In an RT1 context put \(n=X^L\) and use only
singleton or fixed bounded-size same-role groups. After
normalizing the two weighted values to a coprime pair
\(P,Q\) and fixing the deficient record, let
\(W_P,W_Q\) be the sums of the actual conditional group
caps and \(b\) the joint record price. Include coefficients
and discarded full factors in the fixed data. For every
\(\eta>0\), uniformly for sufficiently small fixed \(L\)
and \(\max(P,Q)>n^{C}\) with a sufficiently large fixed
\(C=C(\eta)\), the bound is
\[
             O\!\left(n^{b+\min(W_P,W_Q)/2+\eta}
                              X^\varepsilon\right).
\]
The threshold can also depend on the aggregate raw cap,
group size and polynomial height bounds. The caps can vary
with the record if the estimate is applied in the actual
record bins. An independent full modulus, numerical real
gap, and a low-position charge are unnecessary. In the
cutoff setting one can include \(H/\kappa+z\) in the
minimum, with the prelimit height and zero-window bounds.
Apply (MCB2) to \(P+S=Q\), or \(P+Q=S\), with the
normalized nonzero shift \(S\) fixed and of group cost
zero. The largest complete value is at least \(\max(P,Q)\).
Sum over records. All restrictions of the original rows
remain in place. This supplies the basic half bound
also when the common lock has disappeared in the quotient.
The strict SHF1–SHF3 branch still requires a hypothesis
such as active-floor equality or discrete compression.
If the fixed record puts either variable term into
(MCB3), with a bounded number of powered roots and
an exact integral cofactor of height \(n^{O(1)}\),
the price improves further to \(b+\eta\) at any
nonzero weighted shift with the stated radical-scale
height floor. The integral cofactor condition is stronger
than an \(X^{o(1)}\) numerical height bound.
[Small integral corrections (MCB3)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#compressed-fixed-shifts-without-a-unit-shift-modulus)
[Arbitrary-coefficient half bound (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)

In an ordinary fixed-complete fibre the positive real
height floor is automatic on the original shell. In an
unswitched codegree it must be checked after normalization.
The following criterion permits the normalization for
almost all oriented paths without any assumption on the
radical of the anchor differences.

**RT8.** In a typical \(W\)-locked cell with \(u>0\),
for every fixed \(C>0\) the oriented paths with
distinct anchors for which
\[
        \max(W_1,W_2)/\gcd(W_1,W_2)\le n^C
\]
have count \(n^{2u+\eta}X^\varepsilon\) for arbitrary
\(\eta>0\), after making the fixed chart and cell
allowances sufficiently small. They are power-smaller
than the \(n^{3u-o(1)}\) total paths at every fixed
positive \(u\). No hypothesis on the overlap mass or
the valuation profile of \(W\) is needed.

Fix the first original row and enumerate the
divisor \(g=\gcd(W_1,W_2)\) of \(W_1\) at
height-subpower cost. In the second original row
the common complete value \(U\) is fixed and
\(W_2=g j_2\) with \(j_2\le n^C\). Apply (MCB3)
with no high powered roots on \(W_2\), coefficient
\(g\), and integral variable \(j_2\). The other
variable is the original \(V_2\), with its
radical cap. Its anatomic coefficient fixing can
be paid at an arbitrarily small normalized cost.
The original fixed shift \(U\) and the two free
letters are pairwise coprime, and their maximum
is at least \(X/2\). Thus the completion costs
\(n^\eta X^\varepsilon\) with arbitrarily small
\(\eta\), uniformly in \(g\) and the first row.
There are \(n^{2u+o(1)}\) first rows. Unlike the
weighted normalization used for the other paths,
this argument counts the two free terms in their
original equation with the third cap present.
[Integral-cofactor completion (MCB3)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#compressed-fixed-shifts-without-a-unit-shift-modulus)
· [Oriented paths in a typical cell (RT3)](#from-local-savings-to-balanced-graph-savings)

The negligible-path cutoff can be pushed to
\(n^{L^{-\xi}}\) for some fixed \(\xi>0\) depending
on \(u\) and the count allowance. On the complementary
paths the arbitrary-shift half bound then permits
groups in a bounded number of sufficiently thin
valuation windows, even when a group contains a
growing number of exact orders. The proof chooses
the width after the path and counting thresholds.
[Thin windows and growing short index blocks (TW3–TW4)](ABC_THIN_WINDOW_ARBITRARY_SHIFT_RELEASE.md#consequences-for-progressions-and-balanced-paths)

For genuine residue-cell codegrees one can also lift
each candidate to either original fixed-anchor fibre.
This gives a half-cap and nonatomic bound at the
original numerical height with all three radical caps.
The supported-part cutoff for these paths can then be
taken at a fixed positive numerical height exponent.
[Original-fibre bounds (LFR1–LFR2)](ABC_BALANCED_ORIGINAL_FIBRE_LIFT.md#counting-a-common-neighbour-in-its-original-fibre)

## Residual inequalities and a diffuse test profile

For a fixed orientation, context, and cutoff, evaluate
the actual joint deficiency cost \(b\), the conditional
remaining low price \(s_a\), the dictionary caps, and
the numerical heights after stripping. At saturation
the bounds in RT1 and RT2 cannot be strictly below \(u\),
either uniformly or in the average of RT3. In particular
when \(m_0\) is active and the strictness hypothesis of
RT1 holds, a profile with \(b+m_0\le u\) is excluded.
If a remaining factor compresses with vanishing
correction, it is necessary that \(b+s_a\ge u\) for
an arbitrary shift, or \(b\ge u\) with the independent
terminal. These inequalities describe the cutoff and
compression tests with their conditional data. They
also allow further discounts in the joint deficient
record. A weak limiting degree measure by itself
does not specify those discounts.

To see the limitation of using raw mass instead, let the
restrictive residual measure before choosing a cutoff
have mass \(\ell\), atom \(z\) at zero and complete
variable numerical height \(h\). At a continuity cutoff
\(\kappa\), write \(b_{\rm raw}=\nu((0,\kappa))\) and
\(H(\kappa)=h-\int_{(0,\kappa)}x\,d\nu(x)\). Thus \(H\)
includes the escaped moment. Keeping the zero window
has the raw height-price identity
\[
 b_{\rm raw}+\frac{H(\kappa)/\kappa+z}{2}
     =\frac{\ell}{2}+\frac{b_{\rm raw}}2
       +\frac{H(\kappa)-\kappa(\ell-z-b_{\rm raw})}{2\kappa}.
                                                               \tag{RT4}
\]
The last numerator is nonnegative. Already locked
prime powers are absent from both sides of this
identity. In an unstripped balanced tall profile with
\(\ell=1/3\), raw pricing therefore requires a strict
saving against the deficit and excess on top of
\(1/6\), whereas \(u\) can be any value down to zero.
For a fixed cut at a positive atom one can make the
same calculation with the actual prelimit portion of
its mass retained. One can also price the opposite
letter by its smaller half-bound cap, as in RT1.

For example, take on each letter the formal measure
\((1/3)\operatorname{Unif}[2,4]\) and height one. For
\(2\le\kappa\le4\), the raw active-floor score on the left
of (RT4) is
\[
                        \frac{\kappa}{8}
                          -\frac13+\frac{2}{3\kappa}.
\]
Its global minimum, allowing also \(\kappa<2\) and
\(\kappa>4\), is \((\sqrt3-1)/3>1/6\), attained at
\(\kappa=4/\sqrt3\). Using a raw dictionary half bound
instead gives \(1/6\) with no fixed margin. Neither
retained high distribution can have zero correction
using a fixed number of shrinking windows while
retaining positive mass. Even allowing arbitrary
removal of parts of the measure, any positive-mass
remainder is nonatomic. Accordingly the raw
compression and active-floor equality tests do not
exclude this model. Taking away chart masses tending
to zero from this bounded spectrum leaves the same
limiting example. This is a statement about the
profile tests and does not construct arithmetic
triples or rule out a stronger quantitative gap for
a broad discrete profile. At \(u=0\) the local
obstruction is independent of the cutoff profile.

The one-sided argument uses the coefficient-uniform
fixed-shift estimates for multiple exact valuation
layers. The distinct-coefficient case of the generalized
Fermat result of Browning and Verzobio gives counts
with one prescribed power on each variable. That
estimate can be used independently when actual
divisibility grids have been extracted. A lower
valuation floor and a small joint deficient dictionary
alone do not make those powers exact.
[Browning–Verzobio, Theorem 1.2](https://arxiv.org/html/2608.24512v1#S1)
