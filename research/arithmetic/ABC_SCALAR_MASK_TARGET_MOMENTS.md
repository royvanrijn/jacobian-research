# Scalar products of positive valuation masks

A single unit residue per probe suffices to recover bounded tuples
of actual abc rows at arbitrarily small normalized cost on a
full-divisor chart. The scalar record uses the original integer
letters. For a lock on the target it encodes the free pair \(a,b\)
and makes a product-image formulation of target multiplicities
available at arbitrarily small individual moduli.

## One scalar per probe on two free letters

Use a height-first sequence of balanced boxes as in (PP1). Put
\(n=X^L\), let every letter have height at most \(X^B\) for a
fixed \(B\), and let \(\nu_A,\nu_B,\nu_C\) be the limiting scaled
valuation measures. The total radical cap is \(n^{1+o(1)}\),
and each pair of letters has total nonatomic mass at least
\(1/3\). The notation \(o(1)\) for profile allowances refers
to the outer small-budget limit. Constants in \(X^\varepsilon\)
may depend on the fixed budget. Use a prescribed polynomial-height
prime pool containing \(2\) and the exceptional coefficient
supports, with negligible normalized radical contribution.

Fix an orientation \(U,V,T\) with positive nonatomic mass on \(T\).
Use one affordable complete divisor \(F\mid T\) of positive
numerical weight and of fullness order comparable to \(1/L\).
Its values can be charted at an arbitrarily small fixed normalized
exponent cost. In each chart \(F\) is fixed and the given signed
addition imposes \(U/V\equiv 1\) or \(-1\bmod F\). Both free
letters are units modulo \(F\). The block defining \(F\) comes
from compact positive nonatomic slots as in (MD3).
[Affordable full divisors and atomic completion](ABC_MODULAR_DIGIT_SUPPORT_AND_ATOMIC_COMPLETION.md#completing-all-balanced-profiles-after-small-modulus-support-recovery)

For any fixed mask height allowance \(b>0\), choose a bounded
list of positive quotient probes as in (PP1), with arbitrarily
small allowed missed nonatomic mass. At each probe \(j\) put
\[
       G_j(z)=M_j(U)M_j(V),\qquad
       R_j(z)=G_j(z)\pmod{D_j}.
\]
Here \(z\) is the original row, \(G_j(z)\) is an odd positive
integer, and \(D_j\) is a power of two. The masks are applied
to the two free letters only. This is a scalar record on the
actual integer values. The probe selection in (PP1) gives
\(G_j(z)\le n^{b+o(1)}\). One can also use individual caps
\(G_j(z)\le n^{b_j}\), with any required box slack included in
\(b_j\). Denote the resulting record vector by \(R(z)\).
[Positive quotient masks and support coverage (PP1)](ABC_POSITIVE_PROJECTION_TRANSCRIPT_PRODUCTS.md#bounded-lists-of-short-positive-quotients)

For a fixed word \(w=(w_1,\ldots,w_m)\) of nonzero integer
weights, define \(W=\sum_i|w_i|\). Choose \(D_j\) to be the
first power of two exceeding \(2n^{Wb_j}\). For a common
chosen upper bound \(b_j\le b\) this gives arbitrarily small
positive normalized log modulus per probe, by taking \(b\)
small after fixing the word.

**SR1.** Write \(\tau(w)\) for the minimum size of a set of
indices meeting every nonempty subset of indices with weight
sum zero. For fixed \(w,\eta>0\), the probes and chart can be
chosen so that on any subset of \(N\) rows of the chart, a
prescribed scalar word
\[
                       \prod_{i=1}^{m}R(z_i)^{w_i}
\]
has at most \(O(N^{\tau(w)} n^\eta X^\varepsilon)\) ordered
completions. This includes repeated rows. In particular the
bound is \(O(n^\eta X^\varepsilon)\) for every fixed positive
word and, more generally, for every fixed word with all its
nonempty subset sums nonzero. With positive weights there is
no need for the high-valuation support record in the proof.

For the lift, write \(W_+=\sum_{w_i>0}w_i\) and
\(W_-=-\sum_{w_i<0}w_i\). The integral numerator and denominator
products in the rational word at probe \(j\) are bounded by
\(n^{W_+b_j}\) and \(n^{W_-b_j}\), respectively. They are odd.
If two such rationals agree modulo \(D_j\), their cross
determinant is divisible by \(D_j\) and has absolute value
at most \(2n^{Wb_j}\). They coincide. Thus the exact rational
word at each probe and the union of its numerator and
denominator supports are fixed by the scalar word. The prime
product of that support is of polynomial height for every
nonempty fibre, since it divides a product of the original
letters of one realizing tuple.

Fix the rows in a minimum set of indices meeting the
zero-sum subsets and divide out their contributions. This
costs at most \(N^{\tau(w)}\) and leaves a word with no
nonempty zero-sum subset. Use the same rational lift for the
remaining word, with the same moduli. Choose a large fixed
\(H>0\) and record the union of primes of valuation exceeding
\(H/L\) anywhere in the remaining tuple. Its normalized
radical is at most \(3mB/H\) by the height bound. There are at
most \(n^{3mB/H}\) possible squarefree unions, and the exact
parts of any one realizing row on a fixed such union have
subpower multiplicity by fixed-support enumeration.

In (PP1) take the phase parameter \(\xi\) small enough that
\(\xi H<1\). The order gaps in that proof give
\(eh/(q(q+h))\le \xi Le\), so the mask exponent on a slot
with \(Le\le H\) is either zero or one. Each prime occurs
on at most one role of a given row by primitivity. At any
probe, a prime outside the recorded high-valuation union
that is detected on some row therefore has word valuation
equal to a nonempty subset sum of the weights. It is
visible in the rational word. For an all-positive or
all-negative word, positivity of the masks gives visibility
without a bound on the mask exponent.

The resulting fixed pool contains every prime from the
two free letters detected by the probes on each remaining
row. To count each row, apply the atomic completion argument
of (MD3) with \(F\) as the modulus. Fix the exact parts of
\(U,V\) on the supplied pool at subpower cost. Outside the
pool and outside finitely many positive atomic windows and
the zero window, the remaining radical has arbitrarily small
normalized mass, by the support guarantee in (PP1) and the
choice of essential atoms. Enumerate its support and fix its
exact powers. In each positive window pool the remaining
prime roots at one representative exponent. The width of
the windows tends to zero in the scaled degrees, so the
exponent discrepancies, together with the zero window,
give a rational correction of numerical height \(X^{o(1)}\).
Apply amplified (GW) to the unit congruence for \(U/V\)
modulo \(F\). This counts the ratio at arbitrarily small
normalized cost. Its coprime numerator and denominator
determine the complete letters \(U,V\), and addition gives
\(T\). First choose \(H\) large and then take all support
and completion allowances small enough in terms of \(w,\eta\).
Multiplying the remaining row counts proves (SR1). All
profile choices precede the height limit at each budget.
[Bounded roots and wide-cofactor congruences (GW)](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)

With weights \(1\) in \(k\) slots, (SR1) gives essentially
diagonal multiplicative energy \(O(N^k n^\eta X^\varepsilon)\)
at every fixed order \(k\) on the chart. This follows by
bounding the maximum positive-product multiplicity and
summing the multiplicities, which total \(N^k\). It does
not require the three-row elliptic completion. The earlier
two-coordinate record has the additional feature of being
unchanged by common projective scaling. Here the scalar
uses the masks of the actual free integer letters.
[Two-coordinate positive products (PP2–PP4)](ABC_POSITIVE_PROJECTION_TRANSCRIPT_PRODUCTS.md#two-row-products-have-small-fibres)
· [Signed projective completion (SC3)](ABC_SMALL_MODULUS_TRIPLE_PRODUCT_COMPLETION.md#signed-words-and-separation-of-weighted-slot-sums)

## Signed products on tuples with small support intersections

**SR2.** For any fixed word, consider only tuples for which
the radical of the intersection of the full prime supports
of any two rows is at most \(n^\rho\). Then the scalar word
has \(O(n^{\binom m2\rho+\eta}X^\varepsilon)\) completions
among these tuples. The same bound holds after excluding a
fixed polynomial-height pool in the pair intersection
condition, with that pool adjoined to the completion record.

Indeed the radical of the union of all the pair intersections
outside the prescribed pool is at most \(n^{\binom m2\rho}\).
Enumerate that squarefree union. Every detected mask prime that cancels at its probe is either
in the pool or occurs on at least two original rows. A detected prime in
only one row has a nonzero valuation in the word, whatever
its mask exponent. With the pair union adjoined, the same
atomic completion counts all the rows. This proof needs
neither binary masks nor restrictions on zero-sum subsets.

For arbitrary fixed \(\rho>0\), the pair-intersection condition
can be imposed on almost all ordered tuples across any chart
that retains at least \(n^{1/3-\kappa}\) rows with
\(\kappa<\rho/6\), once the budgets and the auxiliary losses
are sufficiently small. For each fixed preceding row its
complete support is a prescribed polynomial-height pool.
Sharing radical at least \(n^\rho\) discounts the total
radical budget by \(\rho\), and (EB4–EB5) bound the number of
such partners by \(n^{(1-\rho)/3+o(1)}X^\varepsilon\).
Sum over positions of the pair. At fixed \(L\) the inner
\(\varepsilon\) may be taken small relative to \(\rho L\).
This gives a power-small proportion of bad tuples, uniformly
over the prescribed charts. It is a bound on completions
restricted to the retained tuples. An unrestricted signed
energy estimate must also account for the discarded tuples.
[Effective radical prices on polynomial-height supports](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

## A scalar image for weighted target occurrence

**SR3.** Choose \(T=C\) when the target admits the affordable
lock. Keep the original eligibility conditions and the
actual cutoff \(H_M(c)\le Y\). Work in one chart, and let
\(m(c)\) count admissible original rows with complete target
\(c\), with \(m_*(R)=\max_{\operatorname{rad}(c)=R}m(c)\).
For each radical \(R\) with \(m_*(R)\ge D_0\), choose one
maximal representative and retain exactly \(D_0\) of its
rows, where \(D_0\ge1\) is an integer. Use the same scalar
record for all these choices, with a positive word of
fixed length \(k\). Let \(\mathcal W_{k,D_0}\) be the set
of scalar word values obtained by taking the \(k\) rows
over a single retained representative, allowing repetitions.
Then with \(Z(D_0)=\#\{R:m_*(R)\ge D_0\}\),
\[
 \frac{Z(D_0)D_0^k}{O(n^\eta X^\varepsilon)}
            \le |\mathcal W_{k,D_0}|
            \le Z(D_0)D_0^k .
\]
Indeed the total number of ordered tuples is \(Z(D_0)D_0^k\).
The uniform multiplicity bound from (SR1) applies across
targets, since the scalar completion recovers the actual
rows. The upper bound just counts tuples. One may also
use all \(m_*(R)\) rows of each representative, obtaining
the same comparison with \(\sum_R m_*(R)^k\). Fixed-support
enumeration identifies the complete and radical multiplicity
counts to within subpower factors, as in the least-repair
reduction.
[Target radical multiplicity and the exact repair cutoff](ABC_MINIMAL_REPAIR_SUPPORT_OCCURRENCE.md#exact-weighted-form-on-the-repaired-target-radical-alone)

Thus, on balanced primitive height shells where the benchmark
count is \(n^{1/3+o(1)}\), an upper estimate of
\(n^{1/3+(k-1)u-\sigma}\) for this image at
\(D_0=\lceil n^u\rceil\) would give a power saving for that
target multiplicity bin. It would have to hold uniformly
over the original eligible witnesses and prescribed charts,
including \(u=0\). When \(C\) has no affordable lock, one
can orient the scalar on the two free letters of another
positive nonatomic role. The same image comparison remains
valid with tuples still sharing \(c\), although their
records then include one factor from the target. At least
one orientation is always available under the balanced
profile hypotheses.

The bare size of the residue space cannot give this saving
by the present lift and support covering alone. On a nonempty
box covering all but \(\zeta\) of the free pair's nonatomic
mass, the sum of the mask height caps is at least
\(c_U+c_V-\zeta-o(1)\), where \(c_U,c_V\) denote these masses.
Indeed the product of the masks on a row is divisible by the
radical of the union of detected free primes. At the endpoint
\(c_U+c_V\ge1/3\). For all-positive weights \(1\), using the
modulus recipe of (SR1) only gives the ambient-size estimate
\(n^{k\sum_j b_j+o(1)}\). For \(k>1\) and \(u\) in the range
\(0\le u\le 1/6+o(1)\) this cannot reach even the unsaved
target \(1/3+(k-1)u\) as \(\zeta\) tends to zero. The
remaining estimate concerns the subset of the scalar
residue space reached by genuine same-target rows with
their three original caps. The support argument controls
representation multiplicity of this image, while the
fixed-target and minimal-repair arithmetic still have to
bound the image itself.

For the full implication, a sufficient image loss is one absolute
\(\sigma>0\) across all limiting balanced profiles and on a
fixed finite mesh of \(u\), with spacing \(\Delta<\sigma/8\).
Work inside charts whose total pricing cost \(f\) is less
than \(\sigma/8\), take \(\eta<\sigma/8\), and choose the
inner height losses smaller still. The fibre bound (FL1)
allows the mesh to end just above \(1/6\), up to the balance
allowance. For the multiplicity bin starting at \(n^u\),
(SR3) bounds its weight by the number of radicals times
\(n^{u+\Delta}\), with subpower cost for complete
representatives. Summing charts and bins gives the exponent
\(1/3-\sigma+\eta+\Delta+f+o(1)\). These bounds must hold
at each fixed budget in the height limit, with a uniform
power margin as the budget tends to zero. The one-letter
imbalance estimate deals with the complementary unbalanced
families. The same primitive height shells can retain the
least-repair condition since \(c\le H_M(c)\le Y\).
[Fixed-letter bound (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
· [One-letter imbalance estimate](ABC_WIDE_COFACTOR_CONGRUENCES.md#small-radical-on-any-one-letter)
· [Exact repair transfer](POWERFUL_TRIPLES_COMMON_DIVISOR_BOUND.md#a-reverse-transfer-by-repairing-only-the-target)

## Comparison with p-adic projection estimates

Shen's Theorem 1.7 gives an expansion estimate for the
same set under addition and multiplication in the \(p\)-adic
unit group, with a prescribed multiscale nonconcentration
condition. The restricted-graph projection estimate in
Theorem 1.8 likewise needs nonconcentration of the second
input and of the unit coefficient set. In the scalar
record, the product operation acts on the masks, whereas
the exact additive relation acts on the original integers.
An application would need a bridge between these objects
and proof of the requisite scale bounds. Taking the
complete projective values as inputs retains their exact
additive relation but does not follow from the scalar
product moment estimate.
[Shen, Theorems 1.7–1.8](https://arxiv.org/pdf/2607.12251#page=4)
