# Growth of projective shapes in an affordable fibre

A prescribed projective shape confines every private secant prime
to the determinants of its fixed marks. In an affordable
small-budget fibre, typical secants have nearly their entire
numerical weight on private primes. A fixed shape then leaves
essentially only one free row in a typical tuple. This constrains
attempts to use a small projective-shape dictionary together with
a power-class count.

## Secant valuations and the content of a map

Fix \(s\ge4\) distinct rational marks with primitive integer
representatives \(w_i\in\mathbb Z^2\). Let \(v_i\in\mathbb Z^2\)
be primitive representatives of a tuple with this projective
shape. Thus \(v_i=\pm M w_i/c_i\) for a nonsingular integral
matrix \(M\) with relatively prime entries and
\(c_i=\gcd((M w_i)_1,(M w_i)_2)>0\).
Write \(d_{ij}=|\det(w_i,w_j)|\), \(S_{ij}=|\det(v_i,v_j)|\),
and let \(P\) consist of the primes dividing \(\prod_{i<j} d_{ij}\).

**PG1.** At each prime \(p\notin P\), at most one \(c_i\)
is divisible by \(p\), and
\[
      v_p(S_{ij})=v_p(\det M)-v_p(c_i)-v_p(c_j).
\]
If \(p\mid S_{ij}\) for at least one edge, then it divides
at least \(\binom{s-1}{2}\) secants. In particular any prime
exclusive to one secant lies in \(P\).

The valuation identity is the determinant of the two columns.
For \(p\notin P\) the vectors \(w_i,w_j\) are independent
modulo \(p\). Reduction of \(M\) is nonzero by primitivity
of its entries, so it cannot kill both vectors. All edges
avoiding the possible exceptional slot have valuation
\(v_p(\det M)\), the largest valuation of any edge. This
proves the assertions. Also the difference between the largest
and smallest secant valuation is \(\sum_i v_p(c_i)\), so
\[
    v_p\!\left(\frac{\operatorname{lcm}_{i<j}S_{ij}}
                         {\gcd_{i<j}S_{ij}}\right)
                         =\sum_i v_p(c_i) \qquad (p\notin P).
\]
The weighted secant correction in (PM1) therefore has residue
\(-\sum_i e_i v_p(c_i)\bmod q\), and its defect support
outside \(P\) is that of the displayed least-common-multiple
to gcd quotient, since every \(e_i\) is a unit modulo \(q\)
and at most one \(c_i\) is involved. This locates the moving
correction in the row contents when a shape is fixed.
[Secant-corrected classes (PM1–PM3)](ABC_PROJECTIVE_SHAPE_SECANT_CLASS_COUNTS.md#the-secant-correction)

## A fixed shape with one original anchor

Use the balanced affordable-fibre setting of (CR1), with
the original fixed factors \(F,G,D\), distinct targets and
at least one of \(F,G\) of positive numerical weight bounded
below by the affordable margin. Use \(n=X^L\), taking height
first, and write \(N\) for the number of rows in the fibre
with \(N\ge n^r\), \(r>0\) fixed. The shape can be read on the
primitive cofactor pairs \((u_i,v_i)\) of the two varying
terms. This gives the same shape as the ratios of any two
complete original letters, because \(F,G\) are fixed and
the signed equation relates the three ratios by a common
fractional linear map.

**PG2.** Fix \(s\ge4\) and \(\eta>0\). At all sufficiently
late budgets the fibre has a thinning of size \(T\ge Nn^{-\eta}\)
and a family of ordered distinct \(s\)-tuples containing all but
a power-smaller fraction of the tuples in that thinning at each
fixed budget, such that each prescribed projective shape occurs
at most
\[
                            O(T n^\eta X^\varepsilon)        \tag{PG2}
\]
times among those tuples. The implied constant is uniform in
shapes that have a representative of fixed polynomial height.
One can use a single height exponent for every occurring
shape, since normalizing three distinct primitive rational
ratios to \(0,1,\infty\) has polynomial coefficient height.

Choose \(\rho>0\) so that \(2\rho<\eta/2\) and let \(\nu>0\)
be the projective packing gap at class exponent \(\rho\)
and the given affordable margin. Apply (CR1) with a loss
less than \(\eta\) and a packet tolerance \(\epsilon_0\)
so small that \(\epsilon_0 h<\nu/2\), where \(h\) is the
effective numerical scale. This is achieved by fixing the
tolerance before the budget, since \(h\) has a fixed upper
bound from the original height. Let \(M_0\) be the common
projective modulus after that thinning, with \(D\mid M_0\).
The raw secants in the thinning satisfy
\[
                 |\det((u_i,v_i),(u_j,v_j))|
                       \ll M_0 X^h,
\]
and each private packet \(Q_{ij}\) has numerical weight at
least \((1-\epsilon_0)h\) and is coprime to \(M_0\).

Fix a prescribed shape and one actual anchor, say row 1.
The shape pool \(P\) in (PG1) has polynomial prime product
uniformly. Every prime of \(Q_{12}\) lies in that pool, since
it is exclusive to one raw secant on a retained tuple. List
the possible full \(Q_{12}\) supported on \(P\) at subpower
cost. Once a packet is fixed, row 2 is in a single projective
class modulo \(M_0Q_{12}\). Indeed both moduli impose the
class of the fixed anchor, and the primitive coordinates
identify that class prime by prime. The gap between the
numerical cap for the raw secants and the modulus is at most
\(\epsilon_0 h+o_X(1)<\nu\). The projective packing gap
therefore bounds this class by \(n^\rho X^\varepsilon\)
rows, since it may be applied to the quotient modulus
\(M_0Q_{12}/D\). The same independent bound applies to row 3
using \(Q_{13}\). Three distinct rows determine the moving
projective map and every remaining row. Sum over the \(T\)
anchors. The subpower factors and \(2\rho<\eta\) prove (PG2).
[Typical private secant packets (CR1)](ABC_CROSS_RATIO_FRESH_PRIME_OBSTRUCTION.md#exclusive-prime-packets-on-bounded-tuples)
· [Uniform projective packing gap](ABC_AFFORDABLE_FIBRE_SECANT_VALUATION_WEIGHTS.md#setup-and-the-projective-packing-gap)

## Consequence for a shape dictionary

**PG3.** If the retained family in (PG2) has exponent \(r_*\)
in base \(n\), any shape dictionary covering a relation of
\(n^{s r_*-\tau}\) of its typical tuples must have exponent
at least \((s-1)r_*-\tau-\eta\), up to inner allowances.
In particular the cross ratios of typical ordered quadruples
have exponent at least \(3r_*-\eta\). The estimate is
simply the tuple count divided by the uniform bound (PG2).

The two-anchor uncorrected power-class count (PSO2) would
require the shape and selected defect at joint cost below
\((s-2)r_*-\tau\) for a saving on such a relation.
Here the shape alone costs more than that threshold.
The same restriction clarifies the use of the secant
correction. Its prime-pool price on a prescribed shape
can be low because the private primes are in the
shape pool, while the number of shapes still carries
a power cost. These conclusions concern the typical
tuples of the retained affordable fibre. The exceptional
fraction in (CR1) has a power loss at each fixed budget;
that loss is not asserted to be uniform as \(L\to0\).
A thin, specially selected tuple relation or a family
outside the affordable-fibre setup needs a separate count.
[Two-anchor uncorrected count (PSO2)](ABC_PROJECTIVE_SHAPE_TWO_CLASS_OBSTRUCTION.md#an-uncorrected-class-after-fixing-two-original-rows)
