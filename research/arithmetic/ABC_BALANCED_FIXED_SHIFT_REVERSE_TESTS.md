# Energy and compressed shifts after balanced fixed-shift tests

The regular residue-cell reduction forces additive energy in
large multiplicity cells and gives a precise threshold for
compressed differences. Reciprocal sums behave differently.
On an actual positive-exponent fibre, generic signed reciprocal
words retain a substantial radical in their denominators.

## Additive energy and the scale of inverse results

Use the RB1–RB8 conventions. In a regular balanced cohort the
three complete dictionaries have exponent
\(\alpha=1/3-u\), where \(u\in[0,1/6]\) is the common
complete-letter degree exponent. Write \(n=X^L\). In a typical
cell of (RB4), the two complete free dictionaries
\(\mathcal U,\mathcal V\) have size \(n^{u+o(1)}\) and carry
\(n^{2u-o(1)}\) actual rows as the chosen allowances tend to
zero. Each row uses an original marked value \(W\) from a
dictionary of size at most \(n^{\alpha+o(1)}\).

Cauchy–Schwarz on the signed sums giving \(W\) gives
\[
 E_+(\mathcal U,\mathcal V)\ge n^{4u-\alpha-o(1)}
                           =n^{5u-1/3-o(1)}.
\]
Here the mixed energy can be written as the inner product
of the two difference representation functions and is at
most \(\sqrt{E_+(\mathcal U)E_+(\mathcal V)}\).
Including the diagonal, the mixed lower exponent is
\(\max(2u,5u-1/3)\). Using the other dictionary's cubic
upper bound gives an individual lower exponent
\(\max(2u,7u-2/3)\) on each side. The mixed lower bound
first exceeds the pair-count exponent when \(u>1/9\),
and the individual lower bound does so when \(u>2/15\).
At \(u=1/6\), each individual energy has
exponent \(3u\), using \(E_+(\mathcal W)\le|\mathcal W|^3\)
for the other dictionary. This forces energy at the cubic
exponent. It permits a power-small density at every fixed
budget when the outer allowance is fixed, and also a
subpower factor tending to zero.

Consequently the existing \(o(|\mathcal W|^3)\) energy
exclusion at the minimum nonatomic endpoint is compatible
with this lower bound. The bounded-rank affine sparsity
test also does not contradict it. The digit example in
the affine-rank note gives actual integers of full numerical
height with cubic energy exponent and only polylogarithmic
intersection with every fixed-rank feature system tested
there. This example has no imposed three-letter radical
condition. The construction with two free radical caps
and unrestricted target radical supplies a different
arithmetic obstruction, at complete-label degree exponent
zero. It does not realize the positive-multiplicity cells.
[Balanced cell sizes (RB4)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#the-actual-residue-cell-sizes)
· [Qualitative additive-energy exclusion](ABC_COMPLETE_LETTER_ADDITIVE_ENERGY.md#theorem-equivalent-doubling-and-energy-exclusions)
· [Digit example](ABC_AFFORDABLE_FIBRE_AFFINE_RANK.md#high-energy-despite-polylogarithmic-intersections)
· [Two free radical caps](ABC_OVERLAPPING_FULL_RECORD_COUNTEREXAMPLES.md#fixed-parameters-and-complete-integers)

In the three-register extraction (RT6), a positive
dense component has all three dictionary exponents
\(v'\) and row exponent \(2v'\). On this component the
same Cauchy–Schwarz argument gives mixed energy with
exponent \(3v'\). This strengthens the energy signal
when the global marked-letter dictionary is too large
for the one-lock lower bound. It still leaves subpower
losses at the local scale.
[Dense components from three registers (RT6)](ABC_BALANCED_FIXED_SHIFT_REVERSE_TRANSFER.md#a-role-switch-that-preserves-the-full-shift-divisor)

The same-target coloured moment also gives a power
saving on certain nontrivial collisions, including
threefold collisions at \(u>5/36\). Those bounds act on
fibres of the word map and give a lower bound for its
image. The remaining scalar target moment needs an
upper bound for that image over the genuine target
multiplicity tail. An ordinary additive energy lower
bound gives neither this upper bound nor a transfer
from ordinary sums to exact reciprocal sums.
[Coloured collisions (CE1–CE2)](ABC_SAME_TARGET_COLOURED_PRODUCT_ENERGY.md#a-free-support-discount-for-common-target-threefold-energy)
· [Scalar-image obligation (SR3)](ABC_SCALAR_MASK_TARGET_MOMENTS.md#a-scalar-image-for-weighted-target-occurrence)

An SHF2 reduction obtained by paying all of one original
free letter's nonatomic mass does not improve the existing
fixed-letter test. For an original fixed-letter fibre of
exponent \(u>0\), (FL2) already requires each free letter's
nonatomic mass to be at least \(u\). The potential stronger
uses of SHF2 and SHF3 concern **reduced difference fibres**.
They need a cheaper compression or a high original valuation
floor on the residual terms, together with control of the
cost of any free low positions. Cubic energy supplies
neither valuation hypothesis.
[Fixed-letter test (FL2)](ABC_FIXED_LETTER_RADICAL_CAPS.md#strict-positive-scale-inequalities-on-the-two-free-letters)
· [Compressed and high-floor shift tests (SHF2–SHF3)](ABC_FIXED_SHIFT_STRICT_HALF_SAVING.md#a-compressed-term-in-an-arbitrary-fixed-shift-fibre)

## The effective shift threshold in an actual regular cohort

Take two rows of one complete target \(C\), and on a tall
summand role \(W\) write their values as \(U_1,U_2=U_1+h\).
Use the tall normalized threshold of (SCM2), or any stated
fixed positive height threshold from (SCM9). Put
\[
 g=\gcd(U_1,h),\qquad (u_0,H,v_0)=(U_1/g,h/g,U_2/g).
\]
Let \(\mathcal P\) be a prescribed polynomial-height prime
pool. The compressed cost \(d_W\) adds the outside radical
of \(H\) and the support of the multiplier missing from
\(u_0,H,v_0,\mathcal P\). In units of \(\log n\), put
\[
 \beta=\log_n\operatorname{rad}_{\mathcal P^c}(\gcd(g,H)),
 \quad
 \gamma=\log_n\operatorname{rad}_{\mathcal P^c}(\gcd(g,u_0v_0)),
 \quad
 \delta=\log_n\operatorname{rad}_{\mathcal P}(g),\qquad
 D_W=d_W-2\beta-\gamma-2\delta .
\]
In the third term \(\operatorname{rad}_{\mathcal P}(g)\)
means the product of the primes of \(g\) **in** the pool.
This is the signed effective price (SCM8). It prices the
source pair before its target is counted.

Use the raw shift cap of (SCM9). If the pool itself has
bounded normalized radical and \(D_W\) is bounded above,
the raw cap may be deduced. In fact
\(d_W\le D_W+2r_W\), since the three recorded multiplier
supports are disjoint and have combined mass at most the
raw source cap \(r_W\). The outside radical of \(h\) has
cost \(d_W+\gamma\le D_W+3r_W\). Merely prescribing a pool
of polynomial height does not bound the raw shift cost in
this way.

For prescribed lower bounds \(z_C,z_V\) on the actual
overlaps from (SCM10), the completion price is at most
\[
       v_{\rm eff}=
       \min\{u,(r_C-z_C)/2,(r_V-z_V)/2\},
\]
up to the stated cohort and height allowances. The bound
\(u\) comes from the degree of either fixed complete
source value in the regular cohort. Therefore the number
of qualifying same-target ordered source pairs with
\(D_W\le D\) is
\[
           O\!\left(n^{(2r_W+D)/3+v_{\rm eff}+\eta}
                        X^\varepsilon\right).
\]
For \(u>0\) the full ordered-pair population in that cohort
has exponent \(1/3+u\). In balanced caps the qualifying
population is power-smaller when
\[
                      D+3v_{\rm eff}<1/3+3u
\]
by a fixed margin. In particular \(D_W<1/3\) by a fixed
margin on a tall side is already power-small, without
an extra codegree estimate. This conclusion holds for a
prescribed pool and compares with the original population
of pairs. It does not construct qualifying shifts. The raw cost arm of (SCM9) also gives a saving if
\(d_W+v_{\rm eff}<1/6+u\) by a fixed margin. A
common locked divisor may divide a shift, but numerical
height of the quotient and the other residue-cell
identities do not give the needed outside radical cap.
[Signed compressed prices and target discounts (SCM8–SCM10)](ABC_COMPRESSED_SHIFT_TARGET_TAILS.md#endpoint-support-and-original-target-discounts)

## Generic reciprocal sums have large outside denominators

Fix one complete target \(C\) in the original shell with
\(D_C=n^{u+o(1)}\) primitive rows, where \(u>0\) is fixed.
Fix \(k\ge2\) and a polynomial-height prime pool
\(\mathcal P\). For any \(0<\sigma<2u\), all but a
power-small proportion of \(k\)-tuples have, on both
free roles of every slot, private radical mass at least
\(2u-\sigma\) in units of \(\log n\). Here a prime of
one slot is private if it lies outside \(\mathcal P\)
and the full prime supports of the other \(k-1\) rows.

Indeed, condition on all other rows and apply (FL1) to
the chosen slot with their full supports adjoined to the
prescribed pool. On either chosen free role, candidates
with private mass less than \(2u-\sigma\) have exponent
at most \(u-\sigma/2+\eta\). Uniformity of (FL1) in the
pool and fixed target makes this bound usable for every
choice of the other rows. Summing over roles and slots
gives the claimed power-small exceptional proportion
with the allowances smaller than the fixed margin.

Write \(F(x)=x-x^{-1}+(1-x)^{-1}\). For weights
\(w_i\ne0\) whose numerator and denominator primes lie in
\(\mathcal P\), any such private prime is a unique pole
in the sum \(\sum_i w_i F(A_i/C)^{q_i}\) for positive
integers \(q_i\). On a prime of \(A_i\) the inverse of
\(A_i/C\) is the only pole of \(F\), and on a prime of
\(B_i\) it is the inverse of \(1-A_i/C\). The weight
is a unit and all other slots are integral at that prime.
A pole survives summation and these private sets are
disjoint. The reduced denominator therefore has outside
radical at least \(n^{4ku-2k\sigma}\). This statement
uses the simultaneous free radical hypotheses in (FL1).
It does not need an auxiliary upper bound on the exponents
\(q_i\). It restricts attempts to force small-denominator
reciprocal words on a power-dense set of tuples in a
fixed positive-exponent fibre.

There is also a simple obstruction to identifying ordinary
additive energy with reciprocal energy on a fixed target.
On the ordered summand branch \(0<x<1/2\),
\(F''(x)=2(1-x)^{-3}-2x^{-3}<0\). On a prescribed
ordinary sum \(s=x+y\), the sum \(F(x)+F(s-x)\) is
strictly increasing in \(x\) for \(x<s/2\) in the
admissible region. Thus a nontrivial collision of the
ordinary sums of two unordered pairs does not give a
collision of their reciprocal sums. Ordinary integer
rectangles from a residue cell do not supply the
small-denominator hypothesis of the Vandermonde test.
[Fixed-letter radical caps (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
· [Variable-denominator reciprocal sums (RS9)](ABC_THREE_CAP_RATIO_SUMSETS.md#odd-relations-with-variable-output-denominator)
· [Reciprocal Vandermonde test](ABC_RECIPROCAL_VANDERMONDE_SUPPORT.md#label-cancellation-in-capped-additive-triples)

## Surviving conditional profiles across the full range

Write \(t\) for the local comparison threshold. It is
\(u=1/3-\alpha\) in ordinary fibres and unswitched
codegrees, \(2u-c_W\) in the RT5 switch, and \(v'\le u\)
in a positive common-register cell. With three locks
the register constraints give
\[
 \max(0,2/3-3\alpha)\le v\le u.
\]
In particular \(v=0\) is allowed by these constraints
throughout \(\alpha\ge2/9\). For \(t>0\), on a sector
where (RT7) applies, the residual price tests at
saturation are
\[
 b+\tfrac12\min(W_P,W_Q,H/\kappa+z)\ge t
\]
for every available record and cutoff with these
hypotheses. An active height equality is excluded when
the strict terminal hypothesis of (RT1) holds. If a
chosen record leaves a compression, the sharper (RT2)
test is \(b+s_a\ge t\), or \(b\ge t\) when it has an
independent terminal. For varying context prices these
conditions concern the comparison averaged over actual
paths as in (RT3). They do not say each individual path
has a high price. If both free caps are initially
\(1/3\), write \(\Delta_i=1/3-W_i\) including support
stripping and conditional dictionary savings. The
ordinary threshold \(t=u\) requires
\[
                 \Delta_i-2b\le s=1/3-2u
\]
for either free role. A violation by a fixed margin
is enough for a local gain even at a diffuse profile.
The domain where updated caps and heights achieve
this inequality is a condition on the actual arithmetic
records.
The original-fibre lift strengthens these conditions for
genuine residue-cell paths, removing the separate low-price
terminal and adding the remaining nonatomic-mass cap. The
same argument gives a fixed numerical short-quotient exclusion
and a diffuse test stable under the full allowed support
stripping.
[Original-fibre codegrees and residual conditions](ABC_BALANCED_ORIGINAL_FIBRE_LIFT.md#optimized-residual-tests-and-robust-diffuse-profiles)

A formal family survives these comparisons for every
\(0<t\le u\). Give each free role the diffuse measure
\(\tfrac13\mathrm{Unif}[2,4]\) at full height. Let \(F\)
denote the total mass of selected raw positions, with
common recorded positions charged once in the raw
union. Take the conditional record price to be
\(\min(t,F)\). When \(F<t\) by a fixed margin,
assume the remaining bounded groups and finitely
many shrinking windows have their residual raw caps.
In this case the half-cap price is at least
\[
                     \tfrac16+\tfrac12 F\ge t .
\]
Indeed no free role loses more than \(F\).
The height-price branch is bounded below by the
same expression, because residual height divided
by its floor is at least residual mass. If \(F\ge t\),
the record alone reaches the comparison threshold.
No positive-mass residual of the uniform measure has
vanishing correction in a fixed number of shrinking
windows. At \(u=1/6\) the no-record half-cap score
is exactly \(u\), while every active height-floor
score is strictly larger than \(u\) by (RT4). At
\(u<1/6\) the no-record half-cap score already
exceeds \(u\) by a fixed margin for fixed \(u\).

The price law has actual root-dictionary models.
Use the polynomial evaluation encoding in (RO2) on
the union of the slot lists, with interpolation
dimension chosen to give exponent \(t\). Each
projected full-factor record of raw exponent \(F\)
has exponent \(\min(t,F)+o(1)\). If \(F<t\) by a
fixed margin, conditioning on that record leaves
all disjoint fixed-size groups, or groups of
vanishing raw mass, at their raw cardinalities.
The primes form integer power products of the
assigned heights and radicals. This realizes the
projection data of fixed contexts without requiring
a prescribed additive relation or prescribed
residue cells. Separately, the register entropies
have an abstract product model for every allowed
\(\alpha,z_*\): give the registers a shared field
block of exponent \(z_*\) and independent private
blocks of exponent \(s-z_*\). Give each complete
value its two opposite registers and one inner
coordinate of exponent \(v\), with the three inner
coordinates related by addition over the field.
These two models verify the respective price and
entropy obstructions. Neither realizes a saturating
family of primitive integer sums.

At \(t=0\) the codegree comparison gives no
strict loss. In fact even \(q=0\) in (RT3) gives
only a gap \(u/2\), tending to zero as
\(\alpha\to1/3\). For (RT5), \(c_W\ge u\) on a typical
cell since the marked and one free complete value
determine the row. That switch has the same ceiling.
A uniform gain across the range must therefore also
reduce the number of occupied cells or the price
of their records near zero multiplicity.
[Exact projection charges (RO2)](ABC_RESIDUE_OCCURRENCE_ROOT_PROJECTION_OBSTRUCTION.md#an-arithmetic-dictionary-from-finite-field-evaluations)
· [Register bounds and local thresholds (RB5–RB7, RT3, RT5–RT8)](ABC_BALANCED_FIXED_SHIFT_REVERSE_TRANSFER.md#from-local-savings-to-balanced-graph-savings)

## Formal profiles and external techniques

For example, three formal measures
\(\frac13\mathrm{Unif}[0,6]\) have the balanced masses,
first moments one, and no atomic part or escaped weight.
They illustrate why numerical height and the small-budget
mass balance alone give neither a near-three lower degree
floor nor a bounded list of compressed roots. Such data
are a check on the applicability conditions. They are
not an arithmetic construction or a verification of every
other profile inequality in the repository. The arithmetic
counterexamples above have their own stated variations or
missing radical restrictions.

The few-prime-factor result of Agrawal, Bloom and
Petridis has an alternative between growth of a fixed-fold
product and additive-energy control. For complete integer
dictionaries of polynomial height, divisor counting already
gives the near-maximal growth of the fixed-fold product.
That alternative leaves the additive branch optional.
Fang, He and Zhu's products of translates give expansion
in a fixed translate, where the same divisor-count barrier
already applies to the integer dictionaries here. Jain's
Theorem 2 gives averaged distribution for smooth numbers
in progressions through level \(153/256-\varepsilon\).
It tests largest prime factors, a one-variable smooth
population and an averaged modulus range, with logarithmic
relative error. It does not give a uniform power estimate
for two shifted outside radicals in the same cell.
[Agrawal–Bloom–Petridis, Theorems 8–9](https://arxiv.org/html/2512.04931v2#S6)
· [Fang–He–Zhu, Theorems 2–3](https://arxiv.org/html/2611.30174v1)
· [Jain, Theorem 2](https://arxiv.org/html/2611.21136v1#S2)
