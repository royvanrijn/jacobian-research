# Paired Gamma factors and radial separation of factorial moments

Increasing and decreasing exponent rows contribute rising
factorials in opposite directions. Pairing these factors can
reduce the root-modulus uncertainty independently of the total
curvature. In particular, exact mirror pairs can be added without
weakening any annular separation already established for the
remaining rows.

## The beta-binomial annulus

For \(A,B>0\) and \(j\ge1\) define
\[
 P_{j;A,B}(z)=\sum_{k=0}^j\binom jk(A)_k(B)_{j-k}z^k .
\]

**Theorem.** Every root of \(P_{j;A,B}\) satisfies
\[
                 \min(1,B/A)\le |z|\le\max(1,B/A).
                 \tag{1}
\]
Both inequalities are strict when \(j\ge2\) and \(A\ne B\).
For \(A=B\) all roots lie on the unit circle.

The consecutive coefficient ratio gives
\[
 z(1-z)P''+[1-B-j-(A-j+1)z]P'+jAP=0.
\]
Indeed the coefficient of \(z^k\) on the left is
\((k+1)(k+1-B-j)[z^{k+1}]P+(j-k)(A+k)[z^k]P=0\).
There are no roots at zero or one by coefficient positivity.
A repeated root elsewhere would leave a nonzero lowest power
from the second derivative, so all roots are simple.

For a root \(z_i\), subtract \(j-1\) from both sides of its
logarithmic derivative identity to obtain
\[
  2\sum_{l\ne i}\Re\frac{z_i}{z_i-z_l}-(j-1)
       =\Re\frac{B+A z_i}{1-z_i}
       =\frac{B-A|z_i|^2+(A-B)\Re z_i}{|1-z_i|^2}.          \tag{2}
\]
The left side is nonpositive at a root of least modulus
and nonnegative at a root of greatest modulus, since
\(2\Re(1/(1-w))-1=(1-|w|^2)/|1-w|^2\).
By reversal it suffices to consider \(A\ge B\). At least
modulus \(r\), the numerator on the right is at least
\((1+r)(B-Ar)\), giving \(r\ge B/A\). At greatest
modulus it is at most \((1-r)(B+Ar)\), giving \(r\le1\).
If \(A>B\) and the greatest modulus equals one, the
numerator there is \((A-B)(\Re z_i-1)<0\), a contradiction.
If \(j\ge2\), \(A>B\), and the least modulus equals \(B/A\),
then all summands in the formula for the left side of (2)
must vanish, forcing every root to have modulus \(B/A\).
But their product of moduli is
\((B)_j/(A)_j>(B/A)^j\). This proves strictness.

## Matched-factor bounds for the ordinary profiles

For order \(j\), write \(\ell_i=\min(a_i,b_i)\),
\(d_i=|b_i-a_i|\) and \(\sigma_i=\operatorname{sgn}(b_i-a_i)\)
after removing equal rows. Form a list of shapes
\[
        A_{i,h}(j)=\frac{j\ell_i+h}{d_i},\qquad 1\le h\le d_i,
\]
with signs \(\sigma_i\). Put \(C=\prod_i d_i^{\sigma_i d_i}\),
independent of \(j\). Grouping the factors of each factorial
by their residues modulo \(d_i\) shows that, up to a scalar,
the polynomial in \(z=Ct\) is the binomial-basis composition
of the factors with multipliers \((A_{i,h})_k\) for sign
plus and \((A_{i,h})_{j-k}\) for sign minus.

**Theorem.** At each order make any partial matching of the
positive shapes to the negative shapes. In the variable \(z\)
each matched pair \((A,B)\) contributes the closed annulus of
radii
\[
 \left[\max\!\left(\min(1,B/A),\frac{B}{A+j-1}\right),\
       \min\!\left(\max(1,B/A),\frac{B+j-1}{A}\right)\right].
                                                               \tag{3}
\]
An unmatched positive \(A\) contributes
\([1/(A+j-1),1/A]\); an unmatched negative \(B\) contributes
\([B,B+j-1]\). Let \(L_j,U_j\) be the separate products of the
lower and upper radii. Then \(L_j\le C|t|\le U_j\) at every root.
At \(j\ge2\) each bound is strict if there is an unmatched
shape or an unequal matched pair. Thus disjoint intervals
at adjacent orders give consecutive coprimality. If intervals
just touch, a strict bound at the contact also suffices.

For (3) combine (1) with the individual annuli of the
[Gamma factor bounds](GAMMA_PRODUCT_SHARP_LOWER_ANNULUS_AND_SHAPE_SEPARATION.md#2-every-positive-shape-a-sharper-product-annulus).
The same individual bounds give the unmatched radii.
To multiply radii, the
[closed-disc composition lemma](GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md#1-a-self-contained-closed-disc-convolution-lemma)
multiplies upper radii. Reversing the factors and the result
gives the lower product bounds. Strictness follows by shrinking
the radius of the appropriate disc for one factor. For the
lower bound shrink after reversal. When \(A=B\), (3) reduces
exactly to \([1,1]\). Because the scale factor \(C\) is the same
at both orders, their radii are compared in the same variable.

One can re-center the result around the geometric mean of
the endpoints. For a paired factor this multiplies its root
moduli by \(((A)_j/(B)_j)^{1/j}\). This multiplier and \(A/B\)
are between the same two endpoints \(1\) and \(A/B\).
By (1), paired contributions to the absolute log-radius
bound are therefore at most \(|\log(A/B)|\). Unmatched shapes
give the usual one-factor contributions. In particular, if
every shape is matched, the upper bound on the absolute log
radius given by the sum of these costs is minimized by
matching the sorted lists. The inequality
\[
 |x-u|+|y-v|\le |x-v|+|y-u|\qquad(x\le y,\ u\le v)
\]
applied to the logarithms removes inversions in any matching.
The intersections in (3) may improve this centered estimate.

## Consequences for ordinary profiles

At orders five and six, consider the three fully shared rows
\((5,6),(7,6),(29,30)\). Pair the first two and keep the third
unmatched. Here \(C=1\). The new lower order-five bound is
\(31/4500\), and the new upper order-six bound is \(37/5425\).
The latter is smaller by exact cross multiplication. Taking
all factors separately gives overlapping intervals, with
lower order-five radius \(31/4500\) and upper order-six radius
\(6/775\). This example has no exact mirror pair
and has total curvature \(12/35\).

More generally, adjoining any collection of exact mirrored pairs
\((b,b+d),(b+d,b)\) preserves any prior annular separation of
the remaining support. The new shapes match exactly at each order
and their scale factors cancel. In particular, a support with one
row \((a,a+d)\), \(a,d>0\), and any such mirror pairs has
consecutive coprimality at every lower order \(j\) for which
\(a\ge(j-1)d\). Indeed its reciprocal root moduli are trapped in
the single row's interval
\[
            \big[(ja+1)_d,\ (ja+(j-1)d+1)_d\big].
\]
The next lower endpoint is \(((j+1)a+1)_d\), which is no
smaller than the old upper endpoint. The higher order has
strict bounds. Either orientation of the unpaired row is
allowed by reversal.

**Corollary.** There is one integer \(K\), independent of the
number and sizes of the mirror pairs, such that the last family
satisfies SFC(2) at every order whenever \(a\ge Kd\).

In fact the profile is fully shared or has private variables on
both sides, according to whether a mirrored pair is private.
The [uniform eventual theorem for that union](MULTIVARIATE_ALL_DEGREE_REDUCTION.md)
provides one absolute lower-order cutoff \(K\). For smaller
orders the displayed radial comparison gives the conclusion.
This yields an unevaluated absolute ratio threshold. The explicit
finite-order comparison is independent of the eventual theorem.
For ordinary exponent rows all these coprimality statements extend
to characteristic-zero coefficient fields by the rational Bézout
identity.
