# Rational projective parametrizations at the balanced endpoint

Suppose that a family of primitive triples in a fixed
small-budget box is carried by a common map
\[
 [s:t]\longmapsto [F_A(s,t):F_B(s,t):F_C(s,t)],
 \qquad F_A+F_B=F_C ,
\]
of fixed degree \(k\), given by integral homogeneous
binary forms with no common geometric zero. Assume
the coefficient heights are \(X^{o_X(1)}\) in the
inner height limit, uniformly in the possible map.
The coefficients are fixed within the family being
counted. Height losses \(O_k(\epsilon)\) instead
apply if the coefficient bound is \(X^\epsilon\);
strict margins can absorb sufficiently small such
losses. The infinitesimal order in the statement is
the usual height limit at fixed small budget \(L\),
followed by \(L\to0\).

If a limiting \(1/3\) profile has
positive nonatomic mass on all three letters, it
cannot be realized on a fixed-degree carrier as above
of degree \(k\ge3\). More generally, a positive
nonatomic letter with at least three distinct geometric
boundary contacts suffices, including a short letter
with a positive real-height gap. The secondary degree-two case is discussed below.

## Height, local branch weights, and the real correction

Choose coprime integers \(s,t\) and set
\(H=\max(|s|,|t|)\). Since the three forms have no common
zero and one is the sum of the others, they are pairwise
coprime. The homogeneous resultant identities imply
that the gcd of their evaluations divides a nonzero
integer of subpower height, and that
\[
  H^k X^{o_X(1)}
   = \max(|F_A(s,t)|,|F_B(s,t)|,|F_C(s,t)|)
   = X^{1+o_X(1)}
\]
where the first equality allows both signs of the
subpower exponent. Indeed the resultant times
\(s^{2k-1}\), or times \(t^{2k-1}\), is an integral
degree-\((k-1)\) linear combination of the two forms
with coefficients of subpower height. Their values
are the primitive terms up to their common gcd and
sign. Hence \(H=X^{1/k+o_X(1)}\).

Write \(F_T=c_T\prod_j f_{Tj}^{m_{Tj}}\), where \(f_{Tj}\)
are distinct primitive irreducible homogeneous binary
forms of degrees \(d_{Tj}\), and put \(n_T=\sum_j d_{Tj}\).
Each factor and scalar have subpower coefficient height
by the fixed-degree factor bound. For a tall letter, all
its factor values attain their upper height bounds
\(X^{d_{Tj}/k+o_X(1)}\), because their weighted product
attains the sum of the caps. There are at least two such
letters in the signed primitive equation. If there is
a short letter \(S\) of height exponent \(h<1\), use a
bounded real affine chart of the parameter. Distinct
geometric roots of the forms are separated by
\(X^{-o_X(1)}\). This follows by applying the squarefree
discriminant and root bounds to the product of the
distinct dehomogenized factors. Roots outside the chart
do not create small factors in it. Since the short
form loses a positive exponent, exactly one root gives
an exponentially small distance. This root is real
because conjugate nonreal roots are separated. Record it
among at most \(k\) choices and write \(m_*\) for its
multiplicity. Factors other than its own attain their
caps, and the chart interval for the parameter has width
\[
 X^{-(1-h)/m_*+o_X(1)}.
\]
This follows by factoring the dehomogenized polynomial,
since all other root distances, the leading coefficients,
and the height of the projective representatives are
known up to subpower factors.

Define \(D=k(1-h)/m_*\), taking \(D=0\) if all heights
tend to one. The scaled total branch weight
\[
 W_T=k\sum_j\log_X|f_{Tj}(s,t)|
\]
is \(n_T+o_X(1)\) on a tall term and
\(n_S-D+o_X(1)\) on the short term. The determinant of
two primitive parameters in the same real record is
bounded by \(X^{(2-D)/k+o_X(1)}\): in the real chart it
is the product of the two homogeneous denominators times
their real coordinate difference. Without a real gap
this is the usual bound \(2H^2\).

## Projective residue records, including exceptional primes

Choose a fixed product \(J\ne0\) containing the scalar
coefficients, the homogeneous pairwise resultants, and
the derivative resultants of each nonconstant
dehomogenized factor in both affine charts. For primes
outside \(J\), a prime dividing a complete term belongs
to exactly one boundary factor, with original exponent
\(e=m_j u\) when its exponent in \(f_j(s,t)\) is \(u\).
On recording a parameter chart with denominator a
unit and its root branch, this locks the reduced
projective parameter to precision \(p^u\).

On the exceptional support the precision loses at most \(O_k(\log_X|J|)=o_X(1)\)
weight on the entire exceptional support. Indeed
on a chart with denominator a \(p\)-adic unit the
integer resultant identities imply that two distinct
factors cannot simultaneously have valuations greater
than \(v_p(J)\). For a chosen factor with valuation
\(u>2v_p(J)\), its derivative has valuation at most
\(v_p(J)\), and Hensel gives a root in \(\mathbb Z_p\)
to which the chart coordinate is congruent modulo
\(p^{u-v_p(J)}\). Contributions with \(u\le2v_p(J)\)
may be discarded in the precision. A factor constant
on a chart can contribute only via its coefficient,
which is included in \(J\). These statements also
account for the evaluated gcd, since it divides one
of the resultants. The full supported part of each
integer factor, for any factor and prime of \(J\),
has \(X^{o_X(1)}\) possible values: use the fixed-support
divisor estimate for integers of polynomial height
supported on a fixed polynomial-height integer.
Thus one may record it outright, keep its precision
as above, and charge no normalized root cost. There
are at most \(2k\) chart-branch choices per retained
prime and hence divisor-many combined records.

For any selected original prime-power blocks on the
letters, let \(W\) be their total scaled branch weight:
each good prime with exponent \(e\) assigned to a factor
of multiplicity \(m_j\) contributes
\(k(e/m_j)\log_X p\). Also allow the exceptional support
for free, with its aggregate precision loss \(o_X(1)\).
If \(W>2-D\) with a strict margin, two parameters in
one record have cross determinant divisible by
\(X^{W/k+o_X(1)}\). Indeed the chosen parameter chart at
each prime determines their projective cross
determinant modulo \(p^{u-o(u)}\), interpreted with the
explicit precision above. The archimedean bound forces
that determinant to be zero. Thus each lock record
determines at most one projective parameter and one
complete primitive triple.

The standard root-coordinate record cost applies after
splitting each squarefree degree coordinate among the
factors \(f_{Tj}\) and refining the associated dyadic
heights. The number of size refinements is subpower
at fixed box truncation and degree. Prelimit positions
in bounded-scale nonatomic pieces can be selected to
any smaller prescribed mass with arbitrarily small
rounding error. For example, remove windows covering
the finitely many atomic masses larger than a chosen
tolerance, then greedily select positions from the
remaining piece. The tolerance and the windows may
then shrink, using weak convergence and the vanishing
mass of the omitted atoms. Beyond a sufficiently
large fixed scale cutoff, lock the entire tail at
arbitrarily small normalized root cost if desired.
Its radical has the corresponding small bound by the
numerical height cap. For each fixed tail radical the
values supported on it are divisor-many in the
subpower sense, so this costing is uniform.

## Consequences at the endpoint

At normalized count \(1/3\), all letter masses are
\(1/3\). If \(n_T\ge3\) and \(T\) has positive
nonatomic mass, take all its factors except a small
positive-mass diffuse piece at bounded scale.
The cost is strictly below \(1/3\) and the numerical
branch loss can be made arbitrarily small. For a
short \(T\) this compares \(n_T-D\) to \(2-D\), so
the real deficit cancels. For a tall \(T\) the
inequality is at least as strong. The projective
residue count contradicts saturation.

If two letters with positive nonatomic mass have
\(n_T=n_U=2\), first take the full lock on \(T\), of
scaled numerical weight at least \(2-D+o_X(1)\).
Within a factor class the scaled rate per normalized
root mass is \(k r/m_{Tj}\), where \(r=eL\) is the
scaled original degree. If a bounded-scale nonatomic
piece of \(U\) has strictly higher projected rate
than a bounded-scale nonatomic piece of \(T\), replace
a small positive mass from the latter by equal mass
from the former. This gains strictly above \(2-D\).
Omit a smaller piece also to get strict cost below
\(1/3\). Reversing \(T,U\) proves that all projected
rates in their nonatomic parts must be equal to a
common number. This is impossible, since the finite
factor splitting and the nonzero scalar maps
\(r\mapsto kr/m_{Tj}\) preserve nonatomicity on each
part.

Finally, Riemann–Hurwitz for the degree-\(k\) rational
function \(F_A/F_C\) gives
\(\sum_T n_T\ge k+2\), by using total ramification over
\(0,1,\infty\). Thus for \(k\ge3\), either some \(n_T\ge3\)
or at least two \(n_T=2\). When all three nonatomic parts
are positive, each is excluded by the argument above.
[The Stacks Project, Riemann–Hurwitz](https://stacks.math.columbia.edu/tag/0C1B)

For the leftover degree-two boundary count \((2,1,1)\)
with all three nonatomic parts positive, if there is
a short gap and the two-contact letter is tall, its
full weight \(2\) strictly exceeds \(2-D\). If the
two-contact letter is short, its contacts are both
simple and \(D=2(1-h)\). Each tall one-contact letter
has weight one. This again suffices when \(h<1/2\),
while at \(h=1/2\) both tall letters are at equality
and the same exchange argument applies. These
additional cases use the real archimedean improvement.

### Fixed-budget degrees growing in the outer limit

The degree can tend arbitrarily to infinity as the budget
shrinks, while remaining fixed in each inner limit. For
that case use the total unscaled branch weight. The loss
of weight on the short term relative to \(n_S/k\) is at
most \(1-h\), since the original height loses \(1-h\)
and each factor multiplicity is at least one.
The total unscaled branch weight is therefore at least
\[
 \frac{\sum_T n_T}{k}-(1-h)\ge h+\frac2k .
\]
At an actual endpoint \(h\) is bounded below by the
[short-height exclusion](ABC_SHORT_HEIGHT_AND_REAL_GAP.md#a-strict-saving-when-the-short-term-is-small).
Thus some letter has unscaled branch weight at least a
fixed positive number. If all three weak measures have
positive nonatomic mass, omit a small fixed positive mass
from a bounded diffuse band of the chosen letter, whose
choice may be fixed along a subsequence. The unscaled
branch weight loss is at most the omitted numerical
weight on that letter, independent of \(k\). The remaining
weight exceeds the ordinary projective threshold \(2/k\)
when \(k\to\infty\), and the record cost is below \(1/3\).
All coefficient errors are taken arbitrarily small in
the inner limit at the chosen fixed degree. Together with
the bounded-degree cases, this excludes every sequence of
such maps of degree at least three with three positive
nonatomic masses.

### The split quadratic case

For \(k=2\) with boundary counts \((2,1,1)\), each
one-contact form is a scalar times the square of a rational
linear form, because its unique geometric root is rational.
If the two-contact form is split over \(\mathbb Q\), the
three complete terms therefore involve four distinct
rational linear forms. Refine the roots and their size
parameters among these forms. Outside the fixed resultant
support, their radical root costs sum to at most \(1+o(1)\)
in normalized units, since no prime can belong to two
distinct forms. The exceptional support has subpower product
and costs \(o(1)\) in the inner limit. Choose the three
cheapest forms, whose combined mass is at most \(3/4+o(1)\).

The scaled values of these three forms satisfy a fixed
signed linear relation with nonzero coefficients of subpower
height, given by their two-by-two minors. The gcd of any
two scaled values divides the corresponding nonzero
determinant, by the primitivity of \((s,t)\). Hence the
primitive additive triple obtained from their absolute
values has largest term \(X^{1/2+o_X(1)}\), since two
independent linear forms recover the parameters with
subpower coefficients and \(H=X^{1/2+o_X(1)}\).
Its radical is at most \(X^{(3/4+o(1))L}\), including the
subpower coefficient support. Its ordered values and their
signs determine the original projective parameter because
the linear projection to two of the scaled forms is
invertible. The unrestricted limiting one-third count
applied to the derived triple thus counts these parameters
at exponent at most \(1/4+o(1)\) in units of \(L\log X\).
The relative budget in that application is
\((3/2+o(1))L\), fixed with an arbitrarily small allowance
before the inner height limit. This excludes the split
quadratic case at a balanced endpoint, with no diffusion
assumption. In contrast the irreducible quadratic case
with two ramified boundary fibres and the two-contact
letter of height greater than \(1/2\) as the short letter, or
with all three heights one, still contains configurations
outside these lock and linear-form criteria.
