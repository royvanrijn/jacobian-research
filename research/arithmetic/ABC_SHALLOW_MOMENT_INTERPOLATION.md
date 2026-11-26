# Interpolated residue moments for arbitrary valuation profiles

The absolute residue bound can use products longer than the
modulus. Keeping their precise multiplicity in each congruence
class gives a list bound at shallow depth even when some exact
valuation layers have large radical mass.

## A residue-list bound from interpolated moments

At a fixed budget \(L>0\) put \(n=X^L\). Consider representations
\[
              W=c\prod_{i=1}^m x_i^{e_i}
\]
of polynomial height, with \(1\le x_i\le n^{a_i}\),
\(\sum_i a_i=r\), and all coordinates and the coefficient units
modulo \(D\). Here \(D\) has polynomial height and
\(\sigma=\log_n D>0\). The coordinate lists and orders are fixed
before the height limit. The tuples can come from an arbitrary
subset of the root box. For entropy bounds choose one presentation
per complete integer, or state the sample on tuples. Zero cap
coordinates can be removed and harmless dyadic errors absorbed in
the height-subpower allowance.

For each \(i\), let \(g_i\) be the continuous piecewise-linear
interpolant on \([0,1]\) of
\[
         g_i(0)=0,\qquad
         g_i(1/h)=\min\{a_i,\sigma/h\},\quad h=1,2,\ldots.
\]
The interpolant is concave, including at zero. Define the loss
in the full raw count by
\[
       G(\sigma;\mathbf a)=
       \max_{\substack{0\le w_i\le1\\\sum_iw_i\le1}}\,
                                  \sum_i g_i(w_i).
\]
In particular \(G\le\min(r,\sigma)\).
For fixed caps, \(G\) is \(1\)-Lipschitz in \(\sigma\):
each interpolated value at weight \(w_i\) changes by at
most \(w_i|\sigma-\sigma'|\), and the weights sum to at
most one. The same construction handles dyadic size
bounds and limiting depths without a cost proportional
to the number of layers.

**SMI1.** For every list of \(K\) unit residues modulo \(D\), the
number of tuples landing in the list is at most
\[
                K^{1/2} n^{r-G/2}X^\varepsilon
\]
up to a constant for the fixed budget, orders and \(\varepsilon\).
The bound is uniform in the modulus at fixed polynomial height,
the coefficient, the root subsets and the residue list.

For the proof enlarge to the Cartesian product of the root
alphabets. Write \(S_i(\chi)=\sum_x\chi^{e_i}(x)\) over the
\(i\)-th alphabet. Expand the \(2h\)-th moment under normalized
character measure. Fix the first \(h\)-fold product of roots at
cost at most \(n^{ha_i}\) in ordered choices. The other collapsed
product belongs to one of the residues obtained by multiplying
the first product by an \(e_i\)-th root of unity. There are at
most \((2e_i)^{\omega(D)}=X^{o_X(1)}\) such residues. In each,
there are at most \(1+n^{ha_i}/D\) integers in the product
range, each with height-subpower many ordered decompositions.
Thus at the reciprocal weight \(1/h\) the norm is at most
\[
       \|S_i\|_{2h}\ \ll\
              n^{a_i-\min(a_i,\sigma/h)/2}X^{o_X(1)} .
\]
The term \(n^{ha_i}/D\) keeps track of products longer than
the modulus. Interpolation between successive moment indices
gives the \(g_i(w_i)\) saving, with the sup norm at \(w_i=0\).
Apply Cauchy–Schwarz to the character sum of the residue list,
then Hölder on the \(\prod_i S_i\) factor using the weights
\(w_i\). The unused weight is assigned to the constant one.
This proves the claim. The maximizing reciprocal weights can
depend on the modulus and caps. There is no problem with
uniform constants from moments of unbounded index: approximate
the weights from below by a fixed positive cutoff, replacing
smaller weights by the sup norm. The lost total gain is at
most \(m\sigma\) times the cutoff, since \(g_i(w)\le\sigma w\).
Here \(m\) and a polynomial height bound for \(D\) are fixed
at each budget, so the loss can be absorbed in any
\(X^\varepsilon\). All moments used after this choice are
bounded before the height limit.
[Degree-one absolute residue moment at deep moduli](ABC_ABSOLUTE_RESIDUE_TARGET_TAILS.md#one-full-divisor-and-the-multiplicity-threshold)
· [Root blocks at varying depths (RO3)](ABC_RESIDUE_OCCURRENCE_ROOT_PROJECTION_OBSTRUCTION.md#a-residue-bound-at-several-depths)

## Full half-depth occupancy and a deeper bound

**SMI2.** If \(0<\sigma\le r/2\), then \(G(\sigma;\mathbf a)=\sigma\)
for every cap profile. More generally this equality holds if
\[
               \sum_{i:a_i>0}\frac1{\lceil\sigma/a_i\rceil}\ge1 .
\]
Indeed \(g_i(w)=\sigma w\) for
\(w\le1/\lceil\sigma/a_i\rceil\). If an \(a_i\) is at least
\(\sigma\), one such slot suffices. Otherwise
\(1/\lceil\sigma/a_i\rceil\ge a_i/(2\sigma)\), which proves the
uniform assertion.

At joint depth \(\sigma=r\) put \(w_i=a_i/r\).
If \(w_i\) lies between \(p=1/(k+1)\) and \(q=1/k\),
then
\[
 \frac{a_i-g_i(w_i)}{a_i}
      =\frac{(w_i-p)(q-w_i)}{(q-p)w_i}
      \le\frac{\sqrt q-\sqrt p}{\sqrt q+\sqrt p}
      \le 3-2\sqrt2 .
\]
The maximum over \(w_i\) is at \(\sqrt{pq}\), and the largest
bound among consecutive reciprocal integers is at \(k=1\).
Thus \(G(r;\mathbf a)\ge(2\sqrt2-2)r\), including
macroscopic individual layers. The exact profile formula often
gives a higher gain. For example, if the two caps at depth
\(r\) are \(3r/5,2r/5\), the two first segments have total
reciprocal weight \(5/6\) and slope \(r\). The remaining
\(1/6\) should be assigned to the second slot, which has
interpolation slope \(2r/5\) versus \(r/5\) on the first.
Thus this method gives \(G=9r/10\). This calculation identifies
a loss of the moment interpolation at full joint depth and
makes no assertion about the actual residue distribution for
those exponents.

If \(\max_i a_i\le a_*\), the deeper estimate can instead
use
\[
               G(\sigma;\mathbf a)\ge
                            \min(r,(\sigma-a_*)_+).
\]
For a positive \(t=\min(r,(\sigma-a_*)_+)\) take
\(h_i=\lceil t/a_i\rceil\). Then \(h_i a_i\le\sigma\) and
\(g_i(1/h_i)=a_i\ge t/h_i\). If the reciprocal weights
sum to at most one, all can be used and the gain is \(r\).
Otherwise scale them down to sum one. Concavity and
\(g_i(0)=0\) give gain at least \(t\). This also explains
why fine splitting can reach the full-depth gain at
\(\sigma=r\) in the iterated limit.

## Application to occupied residue registers

**SMI3.** For the uniform law on an actual subset of size
\(n^{d+o(1)}\) in the above box, write \(p_\theta\) for the
probability of one residue. The bound
\[
                 \sum_\theta p_\theta^2
                 \ \ll\ n^{2(r-d)-G}X^\varepsilon
\]
is uniform with the same conventions. In fact the \(j\)-th
largest residue count is at most
\(n^{r-G/2}j^{-1/2}X^\varepsilon\), by applying SMI1 to the
first \(j\) residues. The harmonic sum is height-subpower.
The order-two Rényi entropy and therefore the Shannon entropy
are at least \(G-2(r-d)-o(1)\) in base \(n\). For a uniform
law on triples with a free letter of degree at most \(n^{v+o(1)}\),
one can instead use \(d=d_{\rm rows}-v\) in the
collision estimate: multiply the residue list count by the
degree cap and divide by the row count. For a fixed
full-factor record use the actual exponent \(d\) and
remaining raw cap \(r\) in its own fibre. These bounds apply
to every prescribed list and every bin with the required
unit and height condition.

In particular on a balanced occurrence chart with complete
degree exponent \(u\), one has \(r=1/3+o(1)\) and
\(d\ge1/3-u-o(1)\) for a free role, with the small chart cost
included in the outer allowance. Each shallow probe of depth
\(\sigma\le1/6\) on a marked opposite role thus has collision
and Shannon entropy at least \(\sigma-2u-o(1)\), without a
granularity assumption. At \(u=0\) it saturates the alphabet
size at every such depth. For any set of coprime probes of
total depth at most \(1/6\) whose residues are all functions
of the same complete value, the same holds jointly by CRT.
In particular the sum of their marginal entropies exceeds
the joint entropy by \(o(1)\). The statement covers
coprime probes on either or both opposite roles, after
fixing genuine full divisors.

For two coprime shallow registers from different marked roles
of depth \(1/6+o(1)\) each, the pair is encoded by their third
complete role. Its collision entropy is at least
\((2\sqrt2-2)/3-2u-o(1)\) with arbitrary caps, and it reaches
\(1/3-2u-o(1)\) when that role's maximum cap tends to zero.
Thus the argument applies to actual occupied cells even where
complete-letter codegrees give no positive saving. It gives
lower bounds on occupied joint alphabets. A strict
one-third improvement would need an incompatible arithmetic
upper bound for their simultaneous lifts using all three caps,
or a saving in their full-factor record prices.
[Balanced chart degrees and registers (RB3–RB7)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#the-actual-residue-cell-sizes)
· [Conditional exact-layer entropies (LFR3)](ABC_BALANCED_ORIGINAL_FIBRE_LIFT.md#conditional-costs-of-exact-valuation-layers)


## Conditional depth after recording a full factor

**SMI4.** At zero complete multiplicity on a fine-layer
chart, let \(F\mid U\) be any full-factor record on the
free role \(U\), of captured raw mass \(f+o(1)\). Its
choice may depend on the original row. For a prescribed
unit probe \(D\) with \(\sigma=\log_n D\), at balanced
saturation
\[
  H(U\bmod D\mid F)=\min(\sigma,d-f)+o(1),\qquad
  I(F;U\bmod D)=\min(\sigma,d)-\min(\sigma,d-f)+o(1).
\]
These are entropy statements for the conditional law.
For a fixed finite collection of probes on \(U\),
use the depth of their least common multiple in both
formulas.

For the proof, \(H(F\mid U)=o(1)\) since a polynomial-height
integer has height-subpower many divisors, even if the
record choice is adaptive. The record alphabet exponent
is at most \(f+o(1)\) by radical enumeration. In a fixed
record fibre, the remaining raw mass is at most \(d-f+o(1)\),
hence so is the remaining complete-value entropy.
Since \(H(U)=d+o(1)\), both inequalities are equalities
up to \(o(1)\). The same conditioning may further bin the
remaining root sizes at height-subpower cost at the fixed
budget. This gives residual caps with sum \(d-f+o(1)\)
and maximum tending to zero, since each retains the
original exact valuation positions. In each bin apply
SMI1 to \(U/F\). Dyadic conditional residue ranks, as
in LFR3, give
\[
          2H(U\mid F)
            \le 2(d-f)-\min(\sigma,d-f)
                         +H(U\bmod D\mid F)+o(1),
\]
using the fine-layer gain in SMI2 and averaging the bins.
The upper entropy bound is the smaller of the alphabet
depth and the conditional complete-value entropy.
This proves the first formula, and the unconditional
entropy \(H(U\bmod D)=\min(d,\sigma)+o(1)\) gives the second.

For example, when the two opposite shallow registers have
depth \(d/2\) and a full record has \(f\le d/2\), it tells
vanishing normalized information about either register
viewed alone. The pair still determines the record up to
vanishing normalized conditional entropy. The exact upper
bound on the count of F choices and the retained
full-exponent law are essential. This property of one
record does not make several records from different roles
conditionally independent.
[Information pricing in a fixed lock chart (LFR3)](ABC_BALANCED_ORIGINAL_FIBRE_LIFT.md#full-factor-prices-in-the-shared-residue-register)
