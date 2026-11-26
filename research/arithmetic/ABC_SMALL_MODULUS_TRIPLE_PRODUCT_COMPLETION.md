# Small-modulus completion for three positive projection records

Three positive quotient records at a small radical-scale
modulus suffice to recover the three original rows at
arbitrarily small normalized cost inside a full-divisor chart.
The probes distinguish diffuse valuation slots with a fixed
positive scaled separation and avoid a prescribed finite list
of atomic windows. If a prime still disappears from the
product of the three records, its three detected valuations
agree to the chosen scaled precision. The original additive
relation at two higher moduli then recovers the two exact
complementary products.

## Setting and statement

Use a height-first sequence of balanced original valuation boxes
as in (PP1), with \(n=X^L\), all three numerical heights at most
\(X^B\) for fixed \(B\), and fixed limiting scaled measures
\(\nu_U\) for \(U=A,B_0,C\). Here \(A+B_0=C\) and \(B_0\)
denotes the second summand. Each measure has mass \(1/3\),
and each pair has total nonatomic mass at least \(1/3\).
All rows are positive and primitive. Exponents of \(X\) at
fixed \(L\) are taken first, and \(L\) then tends to zero.
Use a prescribed pool \(P\) of polynomial prime product
containing 2 and the fixed exceptional factors of the boxes,
with negligible normalized radical contribution. The
projected roots and positive masks are
\[
 x_q(W)=\prod_{p\notin P}p^{\lfloor v_p(W)/q\rfloor},\qquad
 M_q(W)=x_q(W)/x_{q+h}(W),\qquad
 \phi_q(e)=\lfloor e/q\rfloor-\lfloor e/(q+h)\rfloor.
\]
The integer \(h>0\) will be fixed for the list of probes at a
given budget.

**SC1.** Fix \(\delta,\eta>0\). At sufficiently late budgets,
one can choose \(O_{\delta,\eta,\mathrm{profile}}(1+\log(1/L))\)
probes \(q\) and a chart of the same type as in (PP2), with two
fixed affordable full divisors on two distinct letters of
positive nonatomic mass. Take \(D\) to be the first power of
two exceeding \(n^\delta\). At each probe give a row the two
unit residues
\[
                 \left(\frac{M_q(A)}{M_q(C)},\
                       \frac{M_q(B_0)}{M_q(C)}\right)\pmod D .
\]
For every fixed product of three such vectors the number of
ordered actual row triples in the chart is
\[
                               O(n^\eta X^\varepsilon).
                                                               \tag{SC1}
\]
The full-divisor chart can have dictionary cost smaller than
any prescribed positive normalized margin. All fixed parameters
of the construction are chosen before the sufficiently late
budgets in the sequence. The statement includes repeated rows
and primes shared between distinct rows.

In particular the third multiplicative energy on any subset of
\(T\) actual rows in the chart is \(O(T^3 n^\eta X^\varepsilon)\).
This counts ordered tuples of original rows. The two full
divisors, used only for completion, need not divide the small
modulus \(D\).
[Balanced profiles and affordable divisors](ABC_MODULAR_DIGIT_SUPPORT_AND_ATOMIC_COMPLETION.md#completing-all-balanced-profiles-after-small-modulus-support-recovery)
· [Positive quotient records](ABC_POSITIVE_PROJECTION_TRANSCRIPT_PRODUCTS.md#two-row-products-have-small-fibres)

## Mask height and separation of scaled degrees

The combined limiting measure has an atomless part supported
on the positive axis. Choose a compact carrier \(K\subset[r,R]\)
for all but an arbitrarily small mass of this nonatomic part,
disjoint from all atoms, with \(r>0\). Use slowly shrinking
continuity neighborhoods of \(K\) to select its original exact
valuation slots. Also choose a finite list of distinct positive
atoms containing all but a small mass of the atomic part on the
positive axis. Let the essential windows be slowly shrinking
intervals around these atoms, and a shrinking interval at zero.
The mass outside the carrier slots and these essential windows
is arbitrarily small at late budgets. The widths may be chosen
so that the windows avoid the carrier slots. Let \(s_L\to0\)
satisfy \(s_L\ge L^{1/4}\), with their maximum width \(o(s_L)\).
Put
\[
             Q=\lceil s_L/L\rceil,\qquad
             h=\lfloor \xi LQ^2\rfloor,\qquad Q\le q<2Q
\]
for a fixed sufficiently small \(\xi>0\). Then \(LQ\to0\),
\(LQ^2\to\infty\), \(h=o(Q)\). All choices depending on \(L\)
are fixed before the height limit at that budget.

For compact scaled degrees in \([r,R]\), take \(\xi R<1\)
and let \(q\) be uniform in the indicated range. The mask
values \(\phi_q(e)\) are zero or one. A value is one exactly
when the fractional part of \(e/q\) is strictly less than
\[
            \frac{eh}{q(q+h)}
                      =\xi Le\,(q/Q)^{-2}+o(1).
\]
The strict inequality handles the case of a floor endpoint.
For the mean, the one-phase estimate used in (MR4) holds also
on any fixed positive-length subinterval of \(q/Q\in[1,2]\).
Indeed the same Riemann sum error and integration by parts give
\(O(|te|/Q^2+Q/|te|)\) for each fixed nonzero integral Fourier
frequency \(t\), with constants allowed to depend on the
subinterval. This tends to zero uniformly when \(Le\) is in
the compact positive band. Partition \([1,2]\) into fixed
small intervals, use upper and lower step bounds on the phase
threshold on each, and then let the interval width decrease.
It follows uniformly on the compact band that
\[
                         \mathbb E\phi_q(e)=\xi Le/2+o(1).
                                                               \tag{SC2}
\]
Here \(1/2=\int_1^2 u^{-2}\,du\).

We discard probes at which the three mask product heights are
not certifiably \(n^{\delta/12}\) or less, with a strict margin.
The proof of (PP1) applies also with the present \(Q\): its
uniform phase estimates use only \(LQ\to0\) and \(LQ^2\to\infty\),
and its tail bound uses \(h\le \xi LQ^2\).
In more detail, with \(w_e\) the combined slot masses and
\(u_e\) the carrier slot masses, set
\[
 F_q=\sum_e w_e\phi_q(e),\qquad G_q=\sum_e u_e\phi_q(e).
\]
An upper estimate for \(F_q\) may use the unboxed tail bound
\(O(\xi+1/R')\) above a fixed compact cutoff \(R'\) on \(Le\).
The mixed estimate
\(\mathbb E F_qG_q\ll \xi^2RR'+\xi^2+\xi r'+\xi/R'+o(1)\)
from the proof of (PP1) holds for fixed cutoffs \(r',R'\).
By choosing \(r'\) small and \(R'\) large, and then \(\xi\)
small, the total carrier hit weight on discarded heavy
probes is as small a fixed multiple of \(\xi\) as desired.
The cutoff and bounds are in terms of the fixed boxes,
not of any particular row.
[Positive quotient mask mass (PP1)](ABC_POSITIVE_PROJECTION_TRANSCRIPT_PRODUCTS.md#bounded-lists-of-short-positive-quotients)
· [Phase estimates (MR4)](ABC_MULTI_SCALE_MIXED_RESIDUE_RECOVERY.md#a-bounded-number-of-small-modulus-probes)

Also discard probes at which some essential window can have
nonzero mask. The zero window has \(e<Q\) and is invisible.
For a positive atom \(a_0\), use the nearest integer \(e_0\)
to \(a_0/L\). At both orders \(q,q+h\), the variation in the
real quotient of an order \(e\) inside the window from the
central quotient is \(o(1)\), by the width condition. Changes
in the floors are possible only within \(o(1)\) of an integer
in either central phase, an event of frequency \(o(1)\) by
univariate equidistribution. This estimate is uniform for
the finite atom list. Outside these boundary probes, discard
when the central order has \(\phi_q(e_0)>0\). The carrier
hit weight lost by this last discard is \(O(\xi^2)+o(1)\),
with constants allowed to depend on the chosen finite atoms
and carrier. In fact for \(e\) on the carrier and fixed
\(Le_0\to a_0>0\), the two hits lie in phase intervals of
length \(O(\xi)\). The same bivariate Fourier estimate as
in (MR4) applies off finitely many near resonances
\(|Lte+Lse_0|\le\theta\), with fixed integral \(t,s\).
Their weighted carrier mass tends to zero as \(\theta\)
decreases because the carrier is atomless. Summing over
the finite list proves this claim.

Fix a small positive scaled tolerance \(\omega<r/2\).
After these discards, remove from the carrier the slots
whose hit probability lost on discarded probes exceeds
\(\xi\omega/10\). Their normalized mass is as small as
desired, by the preceding estimates and Markov's inequality,
choosing \(\xi\) after the fixed tolerances. Call the
remaining exact slots good. Every good slot has hit
probability at least \(\xi r/4\) on the admitted probes,
measured in the original uniform sample space. For two
distinct good slots with \(L(f-e)\ge\omega\), (SC2) gives
\[
 \mathbb E\!\left[(\phi_q(f)-\phi_q(e))
                         \mathbf 1_{\mathrm{admitted}}\right]
                 \ge \xi\omega/3
\]
at late budgets, after reducing the numerical constant if
needed. Indeed the full mean gap is \(\xi\omega/2+o(1)\),
and at most the lost hits of \(f\) can decrease it. Thus the
probability of an admitted probe separating them has
the same positive lower bound. At fixed tolerances there are
\(O(L^{-1})\) compact slots and \(O(L^{-2})\) eligible pairs.
A constant times \(1+\log(1/L)\) independent draws, dropping
discarded probes, hit each good slot and separate all the
eligible pairs with positive probability, by the union bound.
Fix such a list before the height limit. The constants in
its size may depend on the tolerances and profile. It has
the following properties.

- The full mask product on each row at each listed probe
  has height at most \(n^{\delta/12}\), with the inner
  margins chosen small enough for rational reconstruction.
- Every essential window is invisible to every listed probe.
- Every good slot is detected, and identical vectors of hits
  on two good slots imply scaled degree distance less than
  \(\omega\).

The slots that belong to neither the essential windows nor
the good set have arbitrarily small total normalized mass.
The size of this exceptional mass is chosen before the
height limit, with the finite list of atoms enlarged first
if necessary. The measures used in the estimates are those
outside the prescribed prime pool \(P\).

## The visible pool and the hidden exponents

For a given modular product, lift both rational components at
each probe exactly. The product of numerator and denominator
in each component is at most \(n^{\delta/4+o(1)}\), so the cross
determinant test modulo \(D\) gives uniqueness. Include the
supports of these lifted rationals in a known prime pool.
Their joint prime product is polynomial in \(X\), because
it divides the product of the original letters of any
realization. Enumerate the union of primes in the exceptional
slots across the three rows at its small normalized radical
cost, and include it and \(P\) in the pool. On the resulting
known pool fix the actual factors of the nine original
letters at subpower cost, uniformly. These factors include
the full valuations, not only the observed masks.
[Fixed-support enumeration](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

Let \(p\) outside that pool occur at a good slot in some row.
It is hit at a listed probe. At that probe the contributions
of the three rows at \(p\) to the two lifted product
valuations are nonnegative multiples of
\((1,0),(0,1),(-1,-1)\), with the direction determined by
the letter of the row carrying \(p\). Their sum is zero
outside the known pool. Thus \(p\) occurs on all three rows,
once in each of the three roles, with equal contributions
at the hit. The other two orders lie on good slots, since
the essential windows are invisible and the exceptional
slots have been recorded. At every probe the three
contributions must be equal, so the three orders have
identical hit vectors. Their pairwise scaled differences
are less than \(\omega\).

Write the exact complementary products as
\[
                U_*=\prod_{i=1}^3 A_i/C_i,\qquad
                V_*=\prod_{i=1}^3 B_{0,i}/C_i .
\]
On the hidden good-slot primes, the numerator and
denominator of either product restricted to these primes
have heights at most \(X^{O(\omega)}\). This follows by
summing the exponent differences, at most
\(\omega/L\) per prime, against the logarithmic radical
cap of the three rows. The same conclusion holds for
the third ratio \(U_*/V_*\).

It remains to handle hidden primes occurring only in essential
windows. For each prime put its memberships in the nine
original letters into one vector. Each position in the
vector is absent, in the window at zero, or in one of
the finitely many positive atom windows. Make one integer
root per nonempty vector, the product of its primes.
These roots have disjoint supports and height \(X^{O(L)}\).
In any of the three aggregate ratios choose a representative
order \(\lfloor a_0/L\rfloor\) for each positive atom and
order zero in the window at zero. Summing these orders
with the signs of the particular ratio gives the exponent
of a vector root. Its scaled limiting signed sum is either
nonzero, in which case the order is comparable in absolute
value to \(1/L\), or zero, in which case the exponent is
\(O(1)\) and can be absorbed into a short correction.
For all vectors the deviations from the representative
orders have aggregate numerical weight \(o(1)\) by the
shrinking windows. Consequently every aggregate ratio
equals a fixed known-pool factor times a monomial in a
bounded number of disjoint short integer roots, multiplied
by one rational factor of height \(X^{O(\omega)+o(1)}\).
No price for knowing the membership vectors of individual
primes is incurred here: the existence of these disjoint
group roots is what enters the bounded-root count.

## Recovering the two exact product ratios

Choose two distinct letters with positive nonatomic mass;
at least two qualify by the pairwise nonatomic condition.
On each, freeze a full divisor from compact positive
nonatomic slots with positive limiting numerical log weight
and fullness parameter comparable to \(1/L\), as in (MD3).
They can be recorded at arbitrarily small exponent cost by
choosing the raw mass of those slots small. On a chart
their values are fixed and each divides its own letter in
every row. A full divisor on \(B_0\) gives \(U_*\equiv1\);
a divisor on \(A\) gives \(V_*\equiv1\); a divisor on \(C\)
gives \(U_*/V_*\equiv-1\), in each case modulo that divisor.
These are unit congruences by original primitivity and
exact addition. Whichever two divisors were chosen, they
give residues of two independent aggregate ratios. The
known-pool factors in those ratios are units at the
corresponding moduli, so may be absorbed in the prescribed
residues.

Apply the amplified wide-cofactor count (GW) at each
modulus. Its hypotheses are the disjoint short integer
roots, signed narrow exponents comparable to \(1/L\),
fullness and positive numerical weight of the modulus,
and a sufficiently short rational correction. These follow
from the preceding representation by taking \(\omega\)
sufficiently small. To make the last choice with fixed
constants, first fix the atom lists, full-divisor scale
bounds and the required normalized counting allowance.
The vanishing-width formulation of amplified (GW) then
implies a positive fixed upper width for the correction
in this allowance at sufficiently late budgets. Otherwise
use a violating subsequence with width tending to zero
and \(L\to0\), contradicting that formulation with the
fixed number of roots, scale bounds and signed exponent
bounds. Thus \(\omega>0\) may be fixed before the probe
construction, with a constant in the logarithmic list
size depending on it. If there are no nonconstant narrow
roots, the short correction is fixed by rational
reconstruction at its full modulus. The two residue
counts supply the two exact ratios at arbitrarily small
normalized cost. In either case the resulting pair
\(U_*,V_*\) is determined, since any two of the three
aggregate ratios are independent.

Finally apply the exact complementary product completion
(EC1) to the two products of \(A_i/C_i,B_{0,i}/C_i\).
It gives subpower multiplicity for the original tuples,
uniformly in the recovered products. The small exceptional
record price, two ratio counts and all inner allowances
can be chosen within the exponent in (SC1). This
completes the proof.
[Amplified wide-cofactor count (GW)](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)
· [Completion for exact complementary products (EC1)](ABC_EXACT_COMPLEMENTARY_PRODUCT_TRIPLE_COMPLETION.md#completion-from-two-exact-products)

## Signed words and separation of weighted slot sums

**SC2.** Fix a word of \(m\) nonzero integer weights
\(w_1,\ldots,w_m\), with \(m\) fixed before the limits.
The probes can be chosen with the same logarithmic bound
so that one prescribed modular word, the product of the
row records raised to these weights, permits at most
\(O(n^\eta X^\varepsilon)\) distinct pairs of exact word values
\[
             \prod_i (A_i/C_i)^{w_i},\qquad
             \prod_i (B_{0,i}/C_i)^{w_i}.
                                                               \tag{SC3}
\]
If at least three of the weights are equal, on a chart of
\(T\) original rows a fixed modular word has at most
\(O(T^{m-3}n^\eta X^\varepsilon)\) ordered original tuple
representations. In particular the signed six-slot energy
test with weights \(1,1,1,-1,-1,-1\) gives
\(O(T^3n^\eta X^\varepsilon)\) tuples for the identity record.
The same energy bound follows by summing squared three-row
product multiplicities from (SC1).

In the probe construction let \(W_*=\sum_i|w_i|\) and use
a mask-height cap \(n^{\delta/(4W_*)}\) with inner margins.
Replace the separation test on pairs of good slots by
separation of weighted sums. Consider at most \(m\) orders
\(e_j\) on good slots, counted with repetitions, with
coefficients \(d_j\) chosen from the weights and their
negatives using at most one occurrence of each weight index.
If
\[
                       \left|L\sum_j d_j e_j\right|\ge\omega,
\]
we require a listed probe with \(\sum_j d_j\phi_q(e_j)\ne0\).
In selecting the good slots lower the lost-hit threshold
to \(\xi\omega/(10W_*)\), keeping an arbitrarily small
exceptional radical mass by the same mixed estimates.
By (SC2) the signed expectation of the mask sum has
absolute value at least \(\xi\omega/2+o(1)\) before the
discard. The change caused by discarding probes is at
most \(\xi\omega/10\) since \(\sum_j|d_j|\le W_*\).
The admitted mask sum is bounded by \(W_*\) in absolute
value, so an admitted probe separates each specified sum
with probability bounded below by a positive constant
depending on the fixed parameters. The number of tests is
\(O_{\mathrm{profile},m}(L^{-m})\). The same union bound gives
a fixed list of \(O(1+\log(1/L))\) probes. Essential windows
are still invisible and high-cost probes are still excluded.

The two rational word components at each probe lift
exactly, since their height products are at most
\(n^{\delta/4+o(1)}\). Fix the pool from their supports
and record the exceptional factors as in (SC1), now over
\(m\) rows. On any prime outside that pool, both signed
mask differences between roles vanish at all listed
probes. Only good slots can appear in a nonzero
contribution, since essential windows are invisible.
Hence in each of the two role differences the signed
valuation sum from good slots has absolute value less
than \(\omega/L\), by the weighted separation property.
The same bound with \(2\omega/L\) holds between the two
non-target roles. All such good-slot terms go into
the short rational corrections. The remaining atomic
memberships now range over the \(3m\) original letters,
still a bounded set of positions. The grouped roots are
disjoint as before, with exponents given by the signed
weighted representative orders. Those with zero leading
scaled exponent go into the correction, and the other
exponents are comparable to \(1/L\).

The same two actual full divisors give two independent
word ratio congruences. The congruences on summand locks
have residue \(1\), and the ratio of the two summands on
a target lock has residue \((-1)^{\sum_i w_i}\). Apply
amplified (GW) as before, choosing the fixed tolerance
\(\omega\) sufficiently small after the finite groups
and full-modulus scales. This yields the stated small
number of pairs in (SC3) after summing over the cheap
exceptional factors. With three equal nonzero weights,
fix the other rows, divide out their contributions and
take the unique positive rational roots for those three
weights, when they exist. The exact completion (EC1)
now applies to the remaining three rows and proves the
tuple bound.

A list bounded independently of \(L\) suffices when the
probe steps are allowed to vary. The finite construction
uses exact Möbius inversion on the scaled valuation
intervals.
[Finite interval probes (FP1–FP2)](ABC_FINITE_PROBE_AND_SIGNED_COMPLEMENTARY_COMPLETION.md#finite-probes-by-möbius-inversion)
