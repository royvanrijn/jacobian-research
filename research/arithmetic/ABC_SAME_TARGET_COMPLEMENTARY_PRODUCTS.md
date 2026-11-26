# Signed complementary products on a common target

When three rows \(A_i+B_i=C\) have the same complete target,
the two exact complementary products with weights \(1,1,-1\)
determine all three rows at subpower cost, away from repeated
rows. Other mixed-sign weights give a power saving under the
original radical caps. At the same time, the positive complete
products occupy almost every allowed residue pair at the level
of exponents in moduli supported on the target primes.

## Subpower completion for the basic mixed-sign word

Fix \(X\ge2\). Let \((A_i,B_i,C)\), \(1\le i\le3\),
be pairwise distinct positive primitive rows with common
\(C\le X\), and put \(x=A_1/C,\ y=A_2/C,\ z=A_3/C\).

**CT1.** For any prescribed positive rationals \(u,v\), the
number of such ordered triples satisfying
\[
  xy/z=u,\qquad (1-x)(1-y)/(1-z)=v
\]
is \(O_\varepsilon(X^\varepsilon)\), uniformly in the
two rationals.

For every \(p\mid C\), both word valuations are
\(-v_p(C)\) because all the summands are units at \(p\).
For any other prime, a negative valuation can arise only
from the third row's summand in the corresponding role,
and its two summands are coprime. Thus the gcd of the two
reduced denominators of \(u,v\) is exactly \(C\).

Put \(t=u+v\), \(K=1-v\), and \(\Delta=x-y\). The equations
give \(xy=uz\) and \(x+y=K+tz\), so
\[
 (t^2z+Kt-2u)^2-(t\Delta)^2=4uv(t-1).                \tag{CT2}
\]
If \(t=1\), the original equations give
\((1-x/u)(y-u)=0\). If \(x=u\) then \(z=y\), and if
\(y=u\) then \(z=x\). Hence \(t=1\) has no pairwise
distinct solution.

If \(t\ne1\), the right side in (CT2) is nonzero.
Let \(d_u,d_v\) be the reduced denominators of \(u,v\)
and set \(J=C(d_ud_v)^2\). Then
\(J(t^2z+Kt-2u)\) and \(Jt\Delta\) are integers.
The difference of their squares is the prescribed nonzero
integer \(4J^2uv(t-1)\) of height \(X^{O(1)}\), since
any nonempty fibre has \(H(u),H(v)\le X^3\).
Factor that integer as the difference-of-squares product.
Each ordered pair of signed factors determines \(z,\Delta\),
since \(t>0\). The sum and difference recover \(x,y\)
and \(C\) is already known. The divisor bound proves (CT1).

Taking positive rational roots gives the same assertion
for weights \((h,h,-h)\), for each fixed nonzero integer
\(h\), whenever the prescribed values admit the roots. In
the full-divisor chart of (FP1) or (SC3), a prescribed signed
modular word of that pattern therefore has
\(O(n^\eta X^\varepsilon)\) distinct three-row
representations on a common target. The uniform exact
completion is what allows summation over the affordable
dictionary of exact product pairs.
[Finite signed probe recovery (FP1)](ABC_FINITE_PROBE_AND_SIGNED_COMPLEMENTARY_COMPLETION.md#finite-probes-by-möbius-inversion)
· [Small-modulus signed word recovery (SC3)](ABC_SMALL_MODULUS_TRIPLE_PRODUCT_COMPLETION.md#signed-words-and-separation-of-weighted-slot-sums)

## Target recovery for arbitrary mixed three-slot weights

Let \(w_i\) be any fixed nonzero integer weights, and consider
one common target \(C\) and complete primitive rows of height
at most \(X\), with exact word values
\[
 u=\prod_i (A_i/C)^{w_i},\qquad
 v=\prod_i (B_i/C)^{w_i},\qquad s=\sum_i w_i.
\]
If \(s\ne0\), at every prime dividing \(C\) the two word
valuations both equal \(-s v_p(C)\). Thus \(C^{|s|}\)
divides the gcd of the two reduced numerators if \(s<0\),
or the gcd of the two reduced denominators if \(s>0\).
Either prescribed gcd gives at most \(X^\varepsilon\)
possible \(C\)'s. If there is only one negative weight
and \(s>0\), the denominator gcd is exactly \(C^s\),
by the same role argument as in (CT1).

If \(s=0\) and \(u\ne v\), both rational word values are
units at every prime of \(C\), and \(u\equiv v\bmod C\)
follows from \(B_i\equiv-A_i\bmod C\). Hence \(C\)
divides the numerator of \(u-v\) in lowest terms. There
are again at most \(X^\varepsilon\) possibilities.

For a word with only one negative weight, the case
\(s=0,\ u=v\) forces all rows equal. Indeed let \(d\)
be the sum of the positive weights, so the negative
weight is \(-d\). With \(C\) canceled, the reduced
denominator of the common value divides both \(A_-^d\)
and \(B_-^d\) of the negative row. It is therefore one.
Weighted Hölder for the positive rows yields
\[
 u^{1/d}
   =\left(\prod_{i:w_i>0}(A_i/C)^{w_i}\right)^{1/d}
       +\left(\prod_{i:w_i>0}(B_i/C)^{w_i}\right)^{1/d}
   \le1 .
\]
A positive integer \(u\) is thus one. Equality makes
the positive rows identical, and the equation for \(u\)
then makes the negative row identical to them. With one
positive row the latter equation alone suffices. Reversing
all signs reciprocates \(u,v\) and covers one positive
weight as well.

**CT3.** Fix positive integer weights \(a,b,d\) and
\(\rho,\eta>0\). On pairwise distinct triples of rows
sharing a target \(X/2<C\le X\), require each primitive
letter's radical to be at most \(n^\rho\), where
\(n=X^L\) and \(L\) is sufficiently small and fixed
before the height limit. A prescribed exact word with
weights \((a,b,-d)\) has
\(O(n^{\rho/4+\eta}X^\varepsilon)\) representations.
The constant is uniform in the prescribed pair and the
same conclusion holds with all signs reversed.

The target-recovery argument leaves only subpower many
candidate targets. Fix one. Use the known polynomial
prime pool of the proof of (SP1), enlarged by the primes
of \(C\). Off the pool, the positive weights imply that
any prime of either positive row also occurs in the
negative row in the same role. An overlap of the
two positive rows in their \(A\)-role requires \(v=1\)
and an overlap in their \(B\)-role requires \(u=1\).
The two equalities cannot both hold for pairwise distinct
rows. If \(s\ne0\) the valuations at \(C>1\) already
exclude either equality. If \(s=0\), the preceding Hölder
argument excludes their conjunction. Choose a
nonexceptional free role. The two positive rows have
disjoint outside supports there, and their combined
radical divides the corresponding radical of the negative
row. At least one of them has an outside cap \(n^{\rho/2}\)
on that free role. The fixed-letter bound (FL1) with the
known pools counts that row at exponent \(\rho/4+\eta\).
After fixing either positive row, there are at most two
distinct positive real completions by the derivative
argument in (SP1). This proves (CT3), after adjusting the
small prescribed exponent allowances.

If \(u,v\ne1\) and \(d\) fails to divide one positive
weight, the exact count improves to \(O(X^\varepsilon)\)
without the radical caps. In fact choose a positive
weight \(a\) with \(q=d/\gcd(a,d)>1\). Every off-pool
prime on either free role of that positive row occurs
only on that row and the negative row. Its exponents
satisfy \(a e_+=d e_-\), hence \(q\mid e_+\).
The two free letters in the chosen row have the forms
\(\alpha r^q,\beta s^q\), with \(\alpha,\beta\) supported
on the known pool. Fix its two pool factors at subpower
cost. The fixed target then gives
\(\alpha r^q+\beta s^q=C\), with polynomial-height
coefficients and coprime roots. Use the coefficient-uniform
fixed-level binary quadratic or Thue–Mahler bound, as
appropriate. Once the chosen row is known, the same
derivative bound completes the other two. This argument
applies for any fixed \(q\ge2\).

The same estimates transfer to a prescribed modular word
after (FP1) recovers the exact pair in its chart. For
a longer fixed word, one may fix the remaining slots
and apply the corresponding three-slot assertion when
their three rows are pairwise distinct. The corresponding
statements count least-repaired primitive quotients in an
original height shell because the common multiplier leaves
both projective coordinates unchanged.
[Fixed-letter bound (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)
· [Mixed-sign support and derivative argument (SP1)](ABC_FINITE_PROBE_AND_SIGNED_COMPLEMENTARY_COMPLETION.md#mixed-sign-completion-using-the-three-radical-caps)
· [Fixed-level binary recovery](ABC_REPAIRED_WEIGHTED_TARGET_DIFFERENCES.md#proof-of-recovery-and-enumeration)

## Local range of positive complete products

Let \(C\ge2\), and allow all free residues
\(a_i\in(\mathbb Z/Q\mathbb Z)^\times\), where every
prime of \(Q\) divides \(C\). Then \(C-a_i\) is also a
unit. At any fixed \(k\ge2\), consider the pair
\[
       U=\prod_{i=1}^k a_i,\qquad
       V=\prod_{i=1}^k(C-a_i)\pmod Q .
\]

**CT4.** Every image pair has \(V\equiv(-1)^k U\) modulo
\(\gcd(Q,C)\). The number of image pairs is at least
\[
      \frac{\varphi(Q)}{32\cdot6^{\omega(Q)}}
                         \frac{Q}{\gcd(Q,C)}.                \tag{CT4}
\]
Thus at polynomial modulus height the loss from the full
set of unit pairs satisfying the unavoidable congruence is
subpower. This holds for arbitrary depths at the target
primes.

It suffices to prove the lower bound for \(k=2\). Fix
each additional \(a_i=1\) when \(k>2\), so the two
product coordinates are rescaled by prescribed units.
For \(k=2\) put \(S=a_1+a_2\), so that
\(V=U+C^2-CS\). At \(p^e\parallel Q\), write
\(r=\max(e-v_p(C),0)\). Given \(U\bmod p^e\), the
residue \(V\bmod p^e\) varies injectively with \(S\bmod p^r\)
over the allowed congruence class. If \(r=0\), every
unit \(U\) occurs by setting \(a_1=1\).

For odd \(p\) and \(r>0\), choose \((U,S)\bmod p\) from
an unordered pair of distinct nonzero roots. There are
\((p-1)(p-2)/2\) choices of those pairs, with distinct
outputs. The proportion of the \((p-1)p\) initial unit
choices of \((U,S)\) is \((p-2)/(2p)\ge1/6\).
Each simultaneous lift \(U\bmod p^e,\ S\bmod p^r\) of
a chosen output is admissible. Extend \(S\) arbitrarily
to modulus \(p^e\) and lift the two simple roots of
\(X^2-SX+U\) by Hensel's lemma.

For \(p=2\) and \(r>0\), prescribe \(U\equiv3\bmod8\)
and \(S\equiv4\bmod8\), projecting either condition
to its indicated modulus when that modulus is below
\(8\). In that case choose representatives extending
the prescribed residue to the stated class modulo \(8\).
The even representative \(S\) has
\((S/2)^2-U\equiv1\bmod8\), so this radicand is a
square unit modulo \(2^e\). This square-root criterion
follows by lifting odd roots from depth three, since
changing an odd root by \(2^{j-1}\) changes its square
by \(2^j\) modulo \(2^{j+1}\) for \(j\ge3\).
If \(q_0\) is a square root, the odd residues
\(S/2+q_0,S/2-q_0\) give the required product and sum.
The restrictions retain at least one quarter of the
unit \(U\)'s and one eighth of the \(S\)'s, even at
depths below three. Apply the Chinese remainder theorem
to the independent prime-power outputs. This proves
(CT4).

The one-coordinate complete scalar is also locally large.
At an odd prime of \(C\) the image of \(a(C-a)\) modulo
\(p^e\) is exactly the negative square-unit coset.
Complete the square as \(C^2/4-(a-C/2)^2\). Translation
by \(C/2\) preserves the units. Subtracting \(C^2/4\)
from a square unit \(s\) preserves the square-unit class,
because \(1-C^2/(4s)\) has a square root by Hensel
lifting. Translation is injective, so this gives the
whole coset. At \(2^e\), with \(C/2\) even, completing
the square gives a translate of the odd squares. With
\(C/2\) odd, choose \(a-C/2=2t\) with \(t\) odd.
The values are a translate of \(4t^2\bmod 2^e\) and
are at least a fixed fraction of the units. Multiplying
independent local choices and fixing \(k-1\) inputs
gives
\[
 \#\left\{\prod_{i=1}^k a_i(C-a_i)\bmod Q\right\}
              \gg 2^{-\omega(Q)}\varphi(Q).
\]
This remains a subpower loss at polynomial modulus height.

In particular, at \(Q=C\) all the input unit residues
lift to positive primitive rows of the actual common
height, simply by taking representatives in \([1,C-1]\).
At larger depths (CT4) concerns local residues. It
does not impose the two original free radical caps
or an Archimedean interval on their lifts. The exact
factorization for two positive rows is
\(U-V=C(A_1+A_2-C)\). The additional integer factor
is a cross difference of complete free letters. The
local range theorem leaves a global estimate for the
actual free-letter radical restrictions, or another
restriction on the occupied scalar-mask image, as the
remaining input.
[Same-target scalar mask image (SR3)](ABC_SCALAR_MASK_TARGET_MOMENTS.md#a-scalar-image-for-weighted-target-occurrence)
