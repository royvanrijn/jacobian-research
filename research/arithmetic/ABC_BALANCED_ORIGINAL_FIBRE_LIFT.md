# Original-fibre bounds for balanced residue-cell codegrees

The common neighbours in a balanced residue cell come from two
original fixed-letter fibres. These fibres retain the original
height and all three radical caps. The full-divisor records can
also be tested against the entropy of the shared residue itself.

## Counting a common neighbour in its original fibre

Use the height-first conventions of (RB1)–(RB8), with
\(n=X^L\), \(X/2<C\le X\), and all three original radical caps.
Take an actual residue cell. Fix two distinct complete anchors
\(V_1,V_2\) in one role, let \(U\) be the common neighbour, and
write \(W_i=sU+tV_i\) for the other complete value, with the signs
prescribed by the original roles. The same formulation applies
after switching the common-neighbour role. For each \(i\),
\[
                  U\longmapsto (U,V_i,W_i)
\]
is an injection into an original fixed-\(V_i\) fibre. The
injection also applies to one original fixed-letter fibre
viewed directly. In a weighted-difference application use this
argument when the candidate lists have these original lifts.

Prescribe a joint full-factor record \(F\) of price \(b\) in base
\(n\). Fix any parts of the free values supported on prescribed
pools, possibly different on the two original lifts. The pools
may depend on the anchors and on \(F\). In particular they may
include the support of the nonzero anchor difference, all
weights in the derived equation, and recorded primes. Each pool
is supported on an integer of fixed polynomial height. Its
complete supported part has height-subpower many choices,
uniformly in the pool. Group and bin **after** these fixings.
On each original lift, the residual factors have their original
exact valuations. The supported parts and every fixed
coefficient are absorbed into that lift's coefficients.

Write \(w_Z\) for an allowed conditional sum of bounded-group
dictionary caps on a free letter \(Z\) of either lift, and
\(c_Z\) for its residual nonatomic mass when its weak profile has
been extracted. These are caps on the conditioned sample,
including the requirement that the other row exists. They need
not be available on all rows of the larger fixed-letter fibre.

**LFR1.** In the iterated limit the per-record count exponent is
at most
\[
       \min_{Z}\{w_Z/2,\ c_Z\}.
\]
For the \(c_Z\) branch use a box with the stated residual weak
profile. For the \(w_Z/2\) branch the estimate is uniform on
fixed finite cost bins at sufficiently small fixed budgets, with
an arbitrary positive count allowance and \(X^\varepsilon\).
In a genuine codegree either lift can be used, so the minimum
may include both computed values and the common value if its
record and remaining factors are priced as well.

Indeed, on the original lift the complete anchor has group
cost zero, the triple is primitive, and its largest value
exceeds \(X/2\). Apply (MCB2) to the actual subset of that
fibre to get the half-dictionary bound. Apply (FT1) with the
original anchor and the fixed support pools to get the
nonatomic bound. A residual exact valuation box and its group
caps describe the same subset under the injection. The
support-fixing choices have subpower cost before the outer
limit. Even when the normalized weighted difference absorbs
the marked full modulus, the full-modulus or real-gap terminal
of the original fixed-letter argument is available. In
particular the remaining low positions on the opposite
letter carry no separate terminal charge in LFR1.
[Balanced codegrees and their weighted shifts](ABC_BALANCED_FIXED_SHIFT_REVERSE_TRANSFER.md#from-local-savings-to-balanced-graph-savings)
· [Arbitrary fixed coefficients (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)
· [Fixed-letter fibre outside prescribed pools (FT1)](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#nonatomic-mass-outside-prescribed-prime-pools)

Fix deficient factors on a restrictive free role \(Z\), leaving
a high part of scaled floor \(\kappa+o(1)>0\) and optionally
a zero window of mass at most \(z\) and vanishing numerical
height. Measure \(H\) as a numerical height cap on the residual
high part alone. The raw half cap is then at most
\(m_0=(H/\kappa+z)/2\); the actual residual grouped cap may
give a smaller bound. Add \(b\) after taking the minimum over
the applicable branches. If \(m_0\) has a fixed positive lower
bound and \(\kappa\) stays in a compact positive interval,
there is a uniform strict gap from \(b+m_0\) at all sufficiently
small fixed budgets, for the specified height, cap and zero
window bounds. To prove this, take a violating sequence with
per-record exponent tending to \(m_0\), taking the height limit
first. FL1 with the actual residual raw caps forces the high mass to tend to
\(H/\kappa\) and the remaining zero mass to tend to \(z\).
The excess numerical weight above the floor then tends to
zero. The high part is its radical to the floor order times
an integral correction of vanishing numerical weight. The
residual measure is supported at the limiting floor and at
zero, including when \(H\) tends to zero. This contradicts
FT1. Thus the active-floor branch of SHF1–SHF3 lifts without
a terminal condition on the normalized difference.

More generally, if the restrictive part groups into a fixed
number of shrinking positive scaled-degree windows and a
cofactor of vanishing numerical weight as in (SHF2), its
remaining measure is atomic. A fixed truncation at a large
scaled degree may first be paid at arbitrarily small
normalized radical cost. Equivalently one may use integral
compression with boundedly many powered short roots:
the correction has negligible mass above any positive
scaled degree and the sum of valuations contributed by the
short roots has bounded radical-weighted mean. Truncating
their multiplicities at any fixed bound gives finitely many
scaled sums. After a subsequence, the mass on sums tending
to infinity vanishes by the numerical height cap, and one
can then increase the multiplicity bound to prove atomicity.
In either case the bound becomes \(b+o(1)\).
For fixed root-count and root-height bounds and each fixed
count allowance \(\eta>0\), some fixed numerical cofactor
width \(\omega>0\) suffices for count \(b+\eta\), at all
sufficiently small fixed budgets. A violating sequence with
\(\omega\to0\) gives the atomic obstruction just proved.
These uniform assertions use finite cutoff and width tests.
The nonatomic inequality itself is a limiting statement on
refined profiles. Averaging it over changing contexts requires
a uniform bound on the bins used in the average.
[The active-floor and compression proofs (SHF1–SHF3)](ABC_FIXED_SHIFT_STRICT_HALF_SAVING.md)
· [Atomicity of compressed remainders (CRP1)](ABC_COMPRESSED_RECORD_ATOMIC_PRICE.md#crp1-a-compressed-remainder-is-atomic)

## A fixed numerical height floor after support stripping

**LFR2.** Suppose the distinct-anchor path count is at least
\(n^{2y+t-o(1)}\), with \(t>0\), and the anchor dictionary
has size at most \(n^{y+o(1)}\). Given arbitrarily small fixed
\(\eta>0\), there is a fixed \(\theta>0\) such that the paths
on which either free computed letter has an outside part at
most \(X^\theta\), relative to any pool prescribed by the
anchor pair, number at most \(n^{2y+\eta}X^\varepsilon\) with
arbitrarily small inner allowance and chart losses. At
\(t>0\) choose these allowances smaller than \(t\).

Fix the anchors and the supported part \(G\) of the indicated
computed value. On the appropriate original fixed-letter
fibre that value is \(Gz\), with a single integer feature
\(z\le X^\theta\). Apply (FT2) with rank one and a count
threshold less than \(\eta\). Sum over the height-subpower
many choices of \(G\) and both orientations. This also
bounds the paths where either \(W_i/\gcd(W_1,W_2)\) is
at most \(X^\theta\), by using the prime pool of the anchor
difference, which contains the gcd support. For an
unweighted residue-cell path this strengthens the
radical-scale short-quotient exclusion. The same proof
bounds the short outside part after a joint record of
price \(b\), at the cost of multiplying the per-context
bound by the number of those records.
[Rank-one fixed-fibre bound (FT2)](ABC_FIXED_FIBRE_TRIANGLE_RADICALS.md#a-small-feature-bound-on-a-single-fibre)

## Full-factor prices in the shared residue register

Use normalized entropy in base \(n\) on one fixed \(W\)-lock
chart. Let \(\Theta\) be the shared residue of \(U,V\) modulo
the lock \(D\), with the sign absorbed on one side. Suppose a
record \(R\), which can depend on the whole row, fixes a full
part \(D_U(R)\) on \(U\). Write
\(f_U=\mathbb E\log_n\operatorname{rad}D_U(R)\).
One can also use a prescribed prime pool, by fixing
the complete supported part at subpower cost. If the
raw cap on \(U\) is \(r_U\), then
\[
          2H(U\mid R)\le r_U-f_U+H(\Theta\mid R)+o(1).
\]
The outer error includes arbitrary count allowances and
the inner height-subpower error. This inequality holds
before balanced saturation, with the numerical weight
of the lock bounded below as in RB2.

To prove it, in each \(R\)-fibre order the residues by
nonincreasing conditional probability, and bin their ranks
in dyadic intervals. For one bin with ranks at most \(2^j\),
the residual \(U/D_U(R)\) has radical cap
\(n^{r_U-\log_n\operatorname{rad}D_U(R)}\). Its allowed
residue classes are obtained by multiplying the bin's at
most \(2^j\) unit classes by the inverse of \(D_U(R)\).
The inverse exists by primitivity. Apply (RB2) to this
alphabet. The entropy price of the bin index is
\(O(\log\log X)/\log n\), since \(D\) has polynomial
height. For a residue of rank \(k\), its conditional
probability is at most \(1/k\). Thus the average normalized
logarithm of the size bound \(2^j\) is at most the
conditional residue entropy plus \(\log_n2\). Average
the alphabet bounds, using a common uniform counting
allowance in (RB2).

**LFR3.** Put \(\Delta_U=r_U+H(\Theta)-2H(U)\).
The resulting necessary information cost is
\[
 f_U\le 2I(U;R\mid\Theta)+I(\Theta;R)+\Delta_U+o(1).
\]
In a balanced chart \(\Delta_U=o(1)\) in the chart and
count allowances, by (RB3)–(RB5), and likewise for \(V\).
If \(R_U\) and \(R_V\) fix full factors of \(U,V\)
with respective captured masses \(f_U,f_V\), then
\[
 \begin{split}
 H(R_U,R_V)&\ge \tfrac12(f_U+f_V)
      +\tfrac12\bigl(I(R_U;\Theta\mid R_V)
                          +I(R_V;\Theta\mid R_U)\bigr)-o(1),\\
 H(R_U\mid\Theta)&\ge \tfrac12(f_U-s)-o(1),
 \qquad H(R_V\mid\Theta)\ge \tfrac12(f_V-s)-o(1).
 \end{split}
\]
These records may be selected jointly from the original
row. Let \(\ell_i=H(R_i\mid i,\Theta)\) for \(i=U,V\).
The first inequality gives
\(f_i\le 2H(R_i\mid\Theta)+I(R_i;\Theta)-2\ell_i+o(1)\).
By conditional data processing and the chain rule,
\[
 I(R_U;R_V\mid\Theta)
     \le I(U;V\mid\Theta)+\ell_U+\ell_V.
\]
The two leakage terms cancel when summing the
inequalities and passing to the joint entropy.
Use \(I(U;V\mid\Theta)=o(1)\) and
\[
 2I(R_U,R_V;\Theta)-I(R_U;\Theta)-I(R_V;\Theta)
   = I(R_U;\Theta\mid R_V)+I(R_V;\Theta\mid R_U).
\]
For the separate conditional bounds use \(H(\Theta)=s+o(1)\).
At \(u=0\), the first inequality gives
\(f_U\le I(\Theta;R)+o(1)\) since \(H(U\mid\Theta)=o(1)\).
In particular a deterministic single-role factor record
with fixed captured mass and raw alphabet bound
\(n^{f_U+o(1)}\) must have entropy \(f_U+o(1)\), even
though the per-cell codegrees have no positive target
exponent. A record that fixes a full
part with a larger raw mass than the indicated information
cost contradicts balanced saturation by a fixed margin.
The inequalities use the actual joint data and give
conditions at every \(\alpha\). They do not assert a cheap
record having the required capture.
[Absolute residue alphabet bound and entropy identities (RB2–RB5)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#regular-balanced-cohorts-and-absolute-residue-pricing)

### Conditional costs of exact valuation layers

In a chart with fixed exact valuation lists and dyadic
raw caps, let \(x_i\) be the residual roots on \(U\) after
record \(R\) is fixed, with orders \(e_i\) and heights at
most \(n^{a_i}\), where \(\sum_i a_i\le A\). A tail
omitted from these lists can be fixed with arbitrarily
small normalized radical cost before choosing the
height limit. For the roots kept one has
\[
        2H(U\mid R)\le \sum_i H(x_i\mid R)
                                + H(\Theta\mid R)+o(1).
\]
This refinement prices the actual conditional layer
alphabets. It requires exact powers on each coordinate.

For finite root alphabets of sizes \(M_i\) in one
\(R\)-fibre, apply the proof of the absolute residue
bound to the product of \(x_i^{e_i}\) and to a list
of at most \(K\) residue classes. For the character sum
at position \(i\), expand its moment of order
\(2h_i\), with \(h_i=\lceil A/a_i\rceil\) for positive
\(a_i\). The collapsed products are at most \(n^{2A}\),
smaller than the fixed lock for sufficiently small
\(L\). The number of \(e_i\)-th roots of unity modulo
the lock is height-subpower at these fixed orders.
Fix one collapsed product and the root of unity to
recover at most one other product, with divisor-many
ordered preimages. The norm is at most
\(M_i^{1/2}X^{o_X(1)}\). Cauchy–Schwarz for the residue
list and Hölder, using \(\sum_i h_i^{-1}\le1\),
bound the count by \(K^{1/2}\prod_i M_i^{1/2}X^{o_X(1)}\).
A known coefficient multiplies the residue list by a
unit. This estimate counts a Cartesian envelope of
the actual tuples, so it requires no probabilistic
independence. Now rank-bin each individual root and
\(\Theta\) in each \(R\)-fibre as above, and use this
alphabet estimate in each joint bin. The number of
coordinates is fixed in the inner limit. The entropy
of the bin choices is height-subpower, and their
averaged logarithmic alphabet bounds are the displayed
conditional entropies. The orders and moment indices here are fixed at the
fixed-budget height limit before taking the outer limit.
[Degree-one residue list moment (AT1)](ABC_ABSOLUTE_TARGET_RAY_MOMENTS.md#algebraic-integers-with-a-ray-residue-dictionary)

### Modular overlap of complementary factor records

Let \(U=\prod_{i=1}^k F_i\) be a deterministic,
disjoint partition into full prime factors on the chart,
with \(k\) fixed. Write \(\Phi_i=F_i\bmod D\), so the
complete shared residue determines the product of the
\(\Phi_i\), and conversely. For a tuple \(Y\) write
\(T(Y)=\sum_i H(Y_i)-H(Y)\). The same RB2 alphabet
estimate, now on the integer \(F_i\) itself, gives
\(2H(F_i)\le \mathbb E\log_n\operatorname{rad}F_i+
H(\Phi_i)+o(1)\). A small finite mass mesh handles
variable mass caps: charge its bin at subpower cost
and use the upper bin endpoint. Sum over \(i\)
and use \(H((F_i)_i)=H(U)\) and
\(2H(U)=r_U+H(\Theta)+o(1)\). This yields
\[
 2T((F_i)_i)-T((\Phi_i)_i)
                     \le H((\Phi_i)_i\mid\Theta)+o(1).
\]
Since \(T((\Phi_i)_i)\le T((F_i)_i)\) by data
processing, the part of the joint factor entropy lost
to dependence is bounded by the split-residue entropy
at a fixed complete residue. For example, at \(u=1/6\)
take two complementary records of mass and entropy
\(1/6\) each. Their total dependence is \(1/6+o(1)\).
Given \(\Theta\), either record's residue determines
the other residue, so \(H((\Phi_i)_i\mid\Theta)\)
is at most \(1/6+o(1)\) and must attain that exponent.
Each record is therefore determined by its residue
and the complete register up to vanishing normalized
conditional entropy. This identifies a concrete
modular matching requirement on such high-entropy
complementary records.

## Optimized residual tests and robust diffuse profiles

Let \(t>0\) be the target exponent in a comparison:
\(t=u\) in an ordinary fibre or ordinary residue-cell
codegree, \(t=2u-c_W\) after switching roles, or \(t=v'\)
in a dense three-register component. On a full-exponent
bin the lift requires, for every available record of
price \(b<t\), residual caps \(w_Z\ge2(t-b)\) and
nonatomic masses \(c_Z\ge t-b\) on every free lift role
to which the estimates apply. For each chosen cutoff
one must also have \(H/\kappa+z\ge2(t-b)\), with positive
active-floor equality excluded. Compression at the
stated width, or a numerically short outside part of a
fixed-pool record, supplies the stronger estimates of
LFR1–LFR2. In the path average formulation use the
average of certified per-context bounds, with the
joint record entropy paid once by the chain rule.
These conditions describe the residual region of
these comparisons using actual conditional data.

There is a uniform upper bound on free stripping in a
full-exponent path comparison. For a fixed pair of
anchors and either computed letter, (FL1) with a pool
prescribed by that pair bounds the candidates with
pool radical mass at least \(1/3-2t+\xi\) by
\(n^{t-\xi/2+o(1)}X^\varepsilon\), for fixed \(\xi>0\).
Summing over at most \(n^{2y+o(1)}\) anchor pairs shows
that this event is power-small at path exponent
\(2y+t\). Thus the fully stripped residual mass is at
least \(2t\) in exponent-full bins in the limit.
For a compact band made almost entirely nonatomic
by trimming sufficiently many atom neighbourhoods,
(FT1) supplies an additional test. Capturing
nonatomic mass greater than \(c_Z-t\) on a prescribed
pool in a fixed-anchor box contradicts full conditional
exponent \(t\). With record price \(b\), the two
capture thresholds are \(1/3-2(t-b)\) for raw mass
and \(c_Z-(t-b)\) for nonatomic mass, provided the
corresponding box is full at exponent \(t-b\).
[Outside radical half bound (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

In particular the nonatomic tests of (CQ3) and (CQ4)
on the ordinary positive-\(u\) cell can use
\(\min(1/3-2u,c_W-u)\) in place of the overlap cap
\(1/3-2u\). For the packet forming either test, approximate
its compact nonatomic bands by continuity regions in
which the remaining atomic mass is less than any prescribed
small allowance. If a pool removes more packet mass than
\(c_W-u\) by a fixed margin, pass to weak limits of the
submeasures in and out of the pool. Domination by the
original profile and continuity of the chosen regions
show that the removed nonatomic mass exceeds
\(c_W-u\) by a fixed margin. FT1 on the appropriate
fixed-anchor lift then bounds the number of such queries
strictly below \(n^u\). This argument is uniform in the
pools for the chosen profile and margins, by a violating
sequence with exponents tending to \(u\). The star
packing and two-copy least-common-multiple proofs use
the remaining packet and the same fixed numerical
width. Shrinking the allowances gives the stated
overlap cap for the nonatomic tests.
[Compact bands and two-copy quantiles (CQ3–CQ4)](ABC_BALANCED_CELL_FULL_FACTOR_QUANTILES.md#compact-high-tail-constraints-at-any-positive-multiplicity)

For a formal test of these inequalities use full
numerical heights and initial measures
\(\mu=\frac13\operatorname{Unif}[2,4]\) on every
role. Allow arbitrary submeasure excision up to
\(1/3-2t\) at zero price, and charge a further joint
raw record of mass \(F\) at \(\min(t,F)\). Where
\(F<t\) and the remaining dictionary prices equal
their raw caps, every restrictive residual measure
has mass \(m\ge2t-F\), first moment at least \(2m\),
and nonatomic mass \(m\). Its raw half-charge plus
the record costs at least \(t+F/2\). For a
retained floor \(\kappa>0\), the density bound \(1/6\)
gives
\[
         \int(x-\kappa)\,d\mu_{\rm rem}(x)\ \ge\ 3m^2
\]
on the retained measure supported in \([\kappa,4]\).
Indeed its mass within distance \(a\) to the right of
\(\kappa\) is at most \(a/6\); integrate the remaining
mass until \(a=6m\). Consequently the raw active-floor
price is at least
\[
       F+\frac m2+\frac{3m^2}{8}
       \ \ge\ t+\frac32t^2
                     \qquad(0<t\le1/6).
\]
The last expression is minimized at \(F=0\) after
putting \(m=2t-F\). This separates it from the
active-floor equality. Similarly, the mass within
distance \(a\) of \(J\) representative degrees is at
most \(Ja/3\). Their minimum total absolute
correction is therefore at least \(3m^2/(2J)\).
An omitted part carried by a numerical cofactor can
be counted with degree centre zero. Thus the model
stays outside every bounded-window, vanishing-width
test even after arbitrary allowed stripping.
Its positive residual numerical height also
survives sufficiently small height cutoffs at each
fixed \(t\).

The full initial measures and the saturated
projection costs are individually realized by the
evaluation-code integer dictionaries. The additional
zero-price stripping calculation is a bound for
submeasures with those formal residual prices.
Realizing that combined data inside a three-cap
additive graph would require further arithmetic.
A stronger quantitative fixed-shift bound away from
the active floor, or a discount in the actual
conditional costs, could change the comparison.
The same weak spectrum can arise on fine divisibility
grids, where the separate \(7/24\) estimate already
gives a saving. Hence a weak profile cannot classify
the arithmetic survivors without exact valuation
information.
[Fine valuation grids and diffuse weak limits](ABC_FINE_VALUATION_DIVISIBILITY.md#consequence-for-a-genuinely-fine-valuation-grid)
[Projection-cost test dictionaries and the register model](ABC_BALANCED_FIXED_SHIFT_REVERSE_TESTS.md#surviving-conditional-profiles-across-the-full-range)

## The remaining three-letter test

For the tall diffuse formal profile, (OE1) allows
individual energy up to exponent \(2p+1/6\) on
a dictionary of local exponent \(p\). The dense
three-register lower bound \(3p\), where \(p\le1/6\),
is compatible, with equality in the powers at
\(p=1/6\). The signed compressed-shift estimate
requires a genuine supply of shifts with effective
price \(D+3v_{\rm eff}<1/3+3u\) in the regular
cohort. A common full divisor alone leaves the
outside radical of its shift quotient uncontrolled.
For reciprocal sums, typical positive-fibre tuples
have private denominator radical of exponent
\(4ku-o(1)\) in a \(k\)-slot word, so ordinary
additive energy cannot be used as a small-denominator
reciprocal-image bound. Full-factor numerical
packing gives separate high-tail and width
constraints, which the tall \([2,4]\) test profile
satisfies at the default query height.
[Power upper bound for energy (OE1)](ABC_OCCURRENCE_DICTIONARY_ADDITIVE_ENERGY.md#additive-energy-with-positive-scale-radical-mass)
· [Energy, shift and reciprocal comparisons](ABC_BALANCED_FIXED_SHIFT_REVERSE_TESTS.md)
· [Numerical tail and two-copy quantile checks](ABC_BALANCED_CELL_FULL_FACTOR_QUANTILES.md#coverage-of-the-numerical-test)

At positive \(u\) bounded away from zero, a uniform
discount of the average certified lifted path price
below \(u\) by a fixed margin would give the (RT3)
saving. The discount may come from an actual joint
record, smaller residual dictionaries, or a stronger
fixed-shift count on broad spectra. As \(u\to0\), even a subpower off-diagonal
codegree bound gives only the RT3 gap \(u/2\).
A sufficient complementary estimate there is an
upper bound \(n^{s-\delta+o(1)}\) on occupied
typical cells on a chart of price much less than
\(\delta\), with the RB4 local dictionary bounds.
LFR3 gives a conditional full-factor information
test at this end of the range. The distribution
estimate for actual occupied cells remains open.
