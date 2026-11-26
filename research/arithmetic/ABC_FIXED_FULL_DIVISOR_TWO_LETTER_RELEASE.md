# Two-letter compression with a fixed large radical-capped divisor

The larger-cofactor terminal can be replaced by a second
modular release when one prescribed coefficient contains
a full divisor of positive numerical height. The two
compressed corrections can then have small fixed numerical
height. No rate of convergence of a weak valuation profile
is needed for the corresponding limiting argument.

## A second modulus from a fixed coefficient

**FD1.** Fix \(A,B,H,w,\kappa>0\) and an integer \(m\ge0\),
and put \(n=X^L\). Use the two compressed roles and third
anatomic role of (DTC1), with the same total root bound
\(m\), narrow height caps \(n^B\), and third-role normalized
raw cap \(A\). Count positive primitive sums of complete
height at most \(X^H\). Suppose one of the two prescribed
fixed coefficients, say \(F_U\), contains a divisor \(M\)
with
\[
            M\ge X^w,\qquad M\text{ is }
                           \lfloor\kappa/L\rfloor\text{-full}.
\]
For every fixed \(\eta>0\) there are \(\omega,L_0>0\)
such that for every fixed \(0<L<L_0\) and every
\(\varepsilon>0\), the count with both integral cofactors
\(u_0,v_0\le X^\omega\) is
\[
                              O(n^\eta X^\varepsilon).
                                                               \tag{FD1}
\]
Both fixed coefficients may otherwise be arbitrary.
The count is uniform in the coefficients and in \(M\)
with these bounds, at each fixed budget and finite
fixed root list. Within-letter overlaps in the compressed
representations are permitted. All original radical and
repair constraints can be retained.

Suppose this fails, choosing a violating sequence with
\(L\to0\) and allowed numerical cofactor exponent tending
to zero. At each budget take the height limit with
the original root lists fixed. The normalized row count
along this sequence is bounded below by \(\eta\).
Colour variable primes by their original complete role at
height-subpower cost so cross-role product ratios do not
cancel upon amplification.

Any block of original roots of \(W\) whose powered
product has positive limiting numerical weight can be
fixed as a full modulus, after dropping its low orders.
Indeed if its weight is at least \(w_1>0\), the orders
below \(\lfloor w_1/(2(A+1)L)\rfloor\) account for at most
\(w_1/2+o(1)\) of the numerical weight. The remaining
power product is a modulus of positive weight and
fullness \(\gg 1/L\). It is a unit modulus for the
complete ratio \(U/V\) by primitivity. After absorbing
the fixed coefficients, the two compressed descriptions
give boundedly many roots of height \(n^{O(1)}\) and a
rational correction of height \(X^{o(1)}\).
Apply the amplified (GW) estimate. Each amplified
coordinate divides the numerator or denominator of the
exact signed monomial because there is no cross-role
cancellation, so the exact fibres and the ordered
decompositions are divisor-bounded at fixed depth.
An exact complete ratio determines the row by
primitivity and the additive relation. Thus this
release has arbitrarily small normalized cost.
[Full-modulus wide cofactor estimate (GW)](ABC_WIDE_COFACTOR_CONGRUENCES.md#a-bounded-number-of-small-roots)
· [Integral cofactors and role colouring (DTC1)](ABC_TWO_COMPRESSED_LETTER_RELEASE.md#releasing-the-compressed-ratio)

It follows that any block of \(W\) with positive
limiting numerical weight must cost at least
\(\eta/2\) to fix. Refine its raw root sizes dyadically,
set aside all positions of raw exponent at least
\(\eta/8\) in base \(n\), and greedily pack the others
into blocks of cost at most \(\eta/4\). The number of
positions set aside and the number of blocks are
bounded in terms of \(A,\eta\). Uniformity of the
first release shows that the powered product of each
cheap block has vanishing numerical weight. Indeed,
a subsequence with positive weight would contradict
the cost bound just proved. Hence \(W\) now has a
description by boundedly many original powered roots
of height \(n^{O(1)}\) and one cofactor of height
\(X^{o(1)}\).

Use the fixed divisor \(M\mid U\). The same amplified
(GW) bound on the compressed ratio \(V/W\) gives
vanishing normalized cost, with no charge for fixing
this modulus within the fibre. The exact ratio
determines the row. This contradicts the positive
count exponent and proves (FD1) with some small
fixed \(\omega\). The argument also applies if the
third variable part starts at vanishing numerical
height, since it may then be its own cofactor.
[Concentration by cheap full blocks (RSG1)](ABC_RADICAL_SCALE_REAL_GAP.md#concentration-and-release-at-a-variable-gap)

## A two-divisor sector at fixed numerical width

**FD2.** Fix \(H,h_0,A_0,\eta>0\). In actual primitive
sums with maximum at most \(X^H\) and individual radical
caps at most \(n^{A_0}\), suppose two prescribed roles
\(U,V\) admit a joint dictionary of divisors
\((F_U,F_V)\) of size at most \(n^f\). Assume
\(\max(U,V)\ge X^{h_0}\) and that the quotients satisfy
\(U/F_U,V/F_V\le X^\omega\) for a sufficiently small fixed
\(\omega>0\) chosen after the caps and the allowance.
Then the count is
\[
                            O(n^{f+\eta}X^\varepsilon).
                                                               \tag{FD2}
\]
The same assertion permits the two quotients to
have boundedly many powered roots of height \(n^B\)
and integral cofactors at most \(X^\omega\) if
one of the two fixed divisors itself is at least
\(X^{h_0/2}\). In either case there is a strict balanced
saving at \(f<1/3\) with a fixed margin, including
zero complete-letter multiplicity.

For the first assertion take \(\omega<h_0/2\) so
one of the two record divisors has height at least
\(X^{h_0/2}\). The radical cap on the corresponding
original letter makes its large divisor have a
full sub-divisor of positive numerical weight.
Indeed the prime powers in that fixed divisor of
orders below \(\lfloor h_0/(4(A_0+1)L)\rfloor\)
have combined weight at most \(h_0/4\) in base
\(X\). The other prime powers are the required
modulus. Fix the two divisors, anatomic-factorize
the third role at the usual arbitrarily small
allowance, and apply (FD1) with zero narrow roots
on the quotient-only input, or with the stipulated
bounded root lists for the second assertion.
Sum over the joint records. Separate the two
possible roles supplying the fixed large divisor
at constant cost.
[Stretched short quotients without a fixed positive-weight modulus (DTC2)](ABC_TWO_COMPRESSED_LETTER_RELEASE.md#a-joint-two-divisor-record-across-the-balanced-range)
· [Anatomic factorization, Lemma 2.2](https://arxiv.org/html/2410.12234v2#S2)

## Conditional nonatomic mass after a fixed record

**FD3.** Consider a fixed-record fibre of actual primitive
sums, in the height-first limiting convention at \(L\to0\).
The record fixes positive coefficients on any of the
three roles, has bounded individual radical caps in base
\(n\), and contains at least one fixed divisor of positive
limiting numerical height. For each \(T\), let \(\nu'_T\)
be its limiting measure of the remaining exact prime
orders after fixing the full supported parts on the
record primes, and let \(c'_T\) be the nonatomic mass
of that measure. If the normalized fibre exponent is
\(e\), then
\[
                             e\le c'_U+c'_V
                       \quad\text{for each pair }U,V .
                                                               \tag{FD3}
\]
Equivalently a joint record of price \(f\) can support
an original count exponent \(d\) only if
\(d\le f+c'_U+c'_V\) in any common residual profile
sector chosen with vanishing exponent loss. For
nonuniform records the statement applies to the
conditional fibres with their actual profiles.

To prove this, fix the full parts on the supports of
the supplied divisors at height-subpower cost.
For either chosen free role select finitely many
positive atoms capturing all but an arbitrarily small
fixed margin of its positive atomic mass, plus a
window at zero. Concentrate the selected positive
atoms in windows of widths tending to zero in the
outer budget limit. The remaining raw mass is at
most \(c'_T\) plus the margin: the mass of exact
orders tending to infinity can be made arbitrarily
small by the numerical first-moment bound. Fix the
remaining full factors at that raw price, using
height-subpower enumeration once the radical is
chosen. A finite exact-order boxing costs only an
inner subpower. In a selected window put the
product of original roots to the minimum order, and
include the positive order differences in its
integral correction. The zero window has vanishing
numerical height because of the original raw cap.
This gives boundedly many powered roots and a
correction at \(X^{o(1)}\) height on each of the two
roles, after paying at most \(c'_U+c'_V\) plus the
chosen margin.
[Atomic residuals and exact window grouping (CRP3)](ABC_COMPRESSED_RECORD_ATOMIC_PRICE.md#crp3-the-converse-at-vanishing-numerical-width)

The original fixed divisor of positive numerical
weight has a positive-weight full sub-divisor by
the radical cap, as in (FD2). If it belongs to
\(U\) or \(V\), use (FD1). If it belongs to the
third role, it directly supplies the modulus to
release \(U/V\) by the amplified (GW) estimate.
In either case the extra cost tends to zero in
normalized units. Take the chosen margin in the
atomic captures and count allowances arbitrarily
small after fixing the limiting profile, which
proves (FD3).

In particular, in a balanced sequence at
\(d=1/3\), a joint record on two prescribed roles
which leaves purely atomic residual measures
must have price at least \(1/3\) in the limiting
normalization. In fact a pair of original roles
has total nonatomic mass at least \(1/3\).
Purely atomic residual measures force the record
to capture that mass on the two roles in the
limiting sense of (CRP1). At least one of the
captured measures has positive mass in a compact
positive-degree band, so full refixing supplies
the fixed divisor of positive numerical weight.
Thus (FD3) applies even when the residual atomic
windows have no stretched cofactor rate. This
conclusion is a constraint on the actual joint
record price, with all original caps retained.

The same inequality gives a reciprocity condition when the
two original nonatomic masses are critical. Write \(c_U,c_V\)
for those masses and \(\beta_U=c_U-c'_U,\ \beta_V=c_V-c'_V\)
for the amount captured by the record. Then
\[
                     f\ \ge\
                     d-(c_U+c_V)+\beta_U+\beta_V .
\]
In particular if \(c_U+c_V=d\), a two-role record capturing
positive nonatomic mass has price at least the raw captured
nonatomic mass. The necessary modulus weight follows by
taking a positive compact band of the captured measure,
as above. This tests an arbitrary two-role support record
with common residual profiles, even when it spans many
original valuation groups.
[Nonatomic mass on every balanced pair](ABC_WIDE_COFACTOR_CONGRUENCES.md#nonatomic-mass-in-the-scaled-degrees)
