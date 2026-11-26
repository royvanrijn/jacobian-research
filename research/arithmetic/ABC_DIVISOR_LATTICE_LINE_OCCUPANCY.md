# Small-radical points on divisor lattices

Fix full divisors on the three roles of an actual primitive sum.
The lattice of the complete values can be sliced into lines with
arbitrary progression steps. When the determinant is close enough
to the area of the complete-value box, each line has a short index
range. The radical caps then give a subpower count per line, even
when the direction or its coordinate gcds have large numerical height.
This replaces a volume bound by its square root in the resulting
divisor-record count.

## DL1. A short lattice line has subpower occupancy

Fix \(h_0,H,a_0,\eta>0\), and put \(n=X^L\). There is a
\(\xi>0\) such that for every sufficiently small fixed \(L>0\)
the following bound is uniform. Let
\[
         T(k)=b_T+a_T k\qquad(T=A,B,C),\qquad a_T,b_T\in\mathbb Z,
\]
with a nonzero direction vector \((a_A,a_B,a_C)\), and let the
allowed integer parameters \(k\) lie in an interval of length at
most \(n^{L^{-\xi}}\). Count positive primitive triples of these
values with \(A+B=C\), maximum between \(X^{h_0}\) and \(X^H\),
and \(\operatorname{rad}(T)\le n^{a_0}\) for all three roles.
For every \(\varepsilon>0\) their count is
\[
                             O(n^\eta X^\varepsilon).
                                                               \tag{DL1}
\]
The bound is uniform in the step of a coordinate progression
and its prime support.

If there are two points, the entire affine line is in the plane
\(A+B=C\). A line through the origin has at most one positive
primitive point, because positive collinear integer triples with
coordinate gcd one are equal. Work off the origin. Define
\(g_T=\gcd(a_T,b_T)>0\) for each role. If a role is constant
and an eligible point exists, its normalized value is one.
The normalized affine forms \(T(k)/g_T\) have coprime
coefficients. The three \(g_T\) are pairwise coprime, since
they divide the corresponding terms of any one eligible row.

Take \(\theta=1/2\) in (TW2), with zero powered roots on the
compressed side, opposite raw cap slightly greater than \(a_0\),
and the fixed polynomial height bound. Choose \(\xi>0\) smaller
than both \(\theta\) and the cofactor-width exponent supplied
by (TW2). Allow a fixed margin in the exponent allowances to
cover fixed-support decompositions and dyadic refinements. Set
\(Z=n^{2L^{-\theta}}\).

For each varying role \(T\), restrict first to its normalized
values \(T/g_T\ge Z\). They are in one unit progression of step
\(M_T=|a_T|/g_T\). If this subset is nonempty choose its least
value \(U\). Every other value is \(V=U+M_T j\) with
\(1\le j\le n^{L^{-\xi}}\). Put \(g=\gcd(U,V)\). The unit
condition implies \(g\mid j\). Enumerate \(g\mid U\) at
height-subpower cost. The primitive equation
\[
                       V/g-M_T(j/g)=U/g
\]
has fixed nonzero shift and an integral cofactor \(j/g\) within
the (TW2) cap. The other variable has radical at most \(n^{a_0}\)
and can be put into an anatomic box with arbitrarily small
normalized cost for fixing the exceptional coefficient. The maximum
is at least \(Z/g\ge n^{L^{-\theta}}\) for small \(L\), and each
complete value is at most \(X^H\). The step is also at most \(X^H\)
if a second value exists. Thus (TW2) gives the required subpower
count of \(V\). A varying coordinate is injective on the line, so
this bounds the number of original rows in this branch.
[Integral cofactors beyond fixed normalized width (TW2)](ABC_THIN_WINDOW_ARBITRARY_SHIFT_RELEASE.md#fixed-shifts-with-a-wider-integral-correction)

It remains to count rows with \(T/g_T<Z\) in every role. If two
of these rows exist, shift the parameter so that the first is
at zero. Let \(\mathbf b'\) and \(\mathbf a'\) be the vectors
of normalized intercepts and normalized slopes. The intercepts
are at most \(Z\) in magnitude. Each slope is the difference of
the two positive normalized values in its role divided by a
nonzero integer parameter difference, so it too has magnitude
at most \(Z\). The two vectors are independent because the
original line avoids the origin. Their cross product is a nonzero
integer multiple of the primitive vector
\((g_A,g_B,-g_C)\), since both normalized vectors are orthogonal
to it. Consequently every \(g_T\le 2Z^2\), and the first row has
maximum at most \(2Z^3\). For small enough \(L\) this is less than
\(X^{h_0}\), since \(\log_X(2Z^3)=6L^{1-\theta}+o_X(1)\).
This contradiction completes the line count. A fixed-width
version with parameter length \(n^{O(1)}\) follows using (MCB3)
and a sufficiently large fixed power of \(n\) in place of \(Z\).
[Fixed-shift completion with integral cofactors (MCB3)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#compressed-fixed-shifts-without-a-unit-shift-modulus)

## DL2. Slicing the actual complete-value lattice

Use the height and radical caps of (DL1). Fix positive divisors
\(D_A,D_B,D_C\) to be imposed on the original complete values,
and put \(P=D_A D_B D_C\). A nonempty record makes these
divisors pairwise coprime. Choose any two complete roles \(U,V\)
and an invertible real linear map \(\mathcal S\) on the pair
of complete coordinates. Suppose the record also stipulates
\(\|\mathcal S(U,V)\|_\infty\le1\). Write
\[
                             \Delta=P|\det\mathcal S|.
\]
For some fixed \(\xi_1>0\), depending on the allowances of (DL1),
if \(\Delta\ge n^{-L^{-\xi_1}}\) the count in the record is
\[
                    O\!\left((1+\Delta^{-1/2})
                               n^\eta X^\varepsilon\right).
                                                               \tag{DL2}
\]
In particular, if \(U\le X_U\) and \(V\le X_V\) are actual
upper caps, take the diagonal map and
\(\Delta=P/(X_U X_V)\). The two roles with the smallest caps
give the strongest coordinate-box result. The shared maximum-height floor suffices
for these two roles.

In the chosen complete coordinates the lattice consists of pairs
with \(D_U\mid U\), \(D_V\mid V\), and with the third
divisor dividing the associated signed sum or difference.
After dividing out the first two divisors, the third congruence
has index equal to the third divisor because its coefficients are
units. Thus the determinant before applying \(\mathcal S\) is
\(P\), and the determinant after applying it is \(\Delta\).

Choose a shortest nonzero vector \(v_1\) in this real lattice.
It is primitive in the lattice and extends to a basis. Reduce a
companion \(v_2\) by a multiple of \(v_1\) so that its projection
along \(v_1\) has length at most \(|v_1|/2\). With
\(\ell_i=|v_i|\), one has \(\ell_2\ge\ell_1\) and the sine
of the angle is bounded below, hence
\(\ell_1\ell_2\asymp\Delta\). The projection perpendicular
to \(v_1\) shows that the integer coefficient of \(v_2\) in the
unit box has absolute value \(O(1/\ell_2)\).
If \(\ell_2\) is larger than the relevant absolute constant,
that coefficient is zero. In this case all points lie on the
origin line, and primitivity leaves at most one row.

Otherwise \(\ell_2=O(1)\). There are
\(O(1+1/\ell_2)=O(1+\Delta^{-1/2})\) possible coefficients
of \(v_2\) and therefore at most this many parallel lines.
The coefficient of \(v_1\) on each line is in an interval of
length \(O(1/\ell_1)=O(\Delta^{-1})\), since
\(\ell_1\ell_2\asymp\Delta\). Pull back the lines to the
integer complete coordinates. The coefficients give integer
parameters for affine linear forms in all three original roles.
Taking \(\xi_1<\xi\) absorbs the absolute interval constants
into the bound \(n^{L^{-\xi}}\). Apply (DL1) on every line.
This proves (DL2), with the stated uniformity. Both the
linear change of coordinates and the lattice lines are fixed
within the chosen record. The same argument works for any
rank-two sublattice of the complete integer coordinates, with
its determinant in place of the divisor product. A translated
lattice requires additional control on the remaining line in
the large second-vector case.

For example a ratio interval \(|A/C-r|\le d\) with
\(C\le X\) allows
\[
                \mathcal S(A,C)=
                   ((A-rC)/(Xd),\,C/X),\qquad
                       \Delta=P/(d X^2).
\]
Thus the area in (DL2) can also incorporate a prescribed ratio
width. This use is homogeneous in the complete coordinates,
so the origin-line case is unchanged.

The classical linear lattice bound in this setting is
\(O(1+\Delta^{-1})\) for primitive rows. It follows from the
same origin-line split and counting all lattice parameters in
the other case. The square root in (DL2) instead counts only
the parallel lines and uses the actual radical caps on them.
[Classical geometry bound in Bernert–Browning–Lichtman–Teräväinen, Proposition 4.1](https://arxiv.org/html/2410.12234v2#S4)

## DL3. Summing actual records

Let an actual shell family be covered by records consisting of
the three divisors, along with any chosen coordinate caps or
ratio windows. If \(\Delta_R\) is the scaled determinant for
record \(R\), and the (DL2) threshold holds for all counted records,
the total count is at most
\[
              O\!\left(n^\eta X^\varepsilon
                       \sum_R(1+\Delta_R^{-1/2})\right).
                                                               \tag{DL3}
\]
In particular if the actual joint dictionary has size at most
\(n^f\), and every record has \(\Delta_R\ge n^{-v}\) for
some \(v\ge0\) in the (DL2) range, its count exponent is
at most \(f+v/2+\eta\). In a balanced family this gives
a fixed saving whenever \(f+v/2<1/3\) by a fixed margin.
The condition applies down to zero complete-letter multiplicity.
All original eligibility conditions, including a least-repair
cutoff, may remain imposed.

For a coordinate box an equivalent sufficient condition with
a fixed exponent margin \(\sigma>0\) is
\[
     P\ge X_U X_V\,n^{-2(1/3-f-\sigma)},\qquad f<1/3-\sigma .
\]
If ratio windows are used, replace \(X_U X_V\) by \(dX^2\)
and charge the windows in the joint record price. The range
\(v\le L^{-\xi_1}\) describes the per-record estimate.
A one-third saving from the displayed simple tally uses a
bounded \(v\) satisfying the stated inequality.

The progression step on a sliced line is a component of a
lattice basis vector, and its radical can be uncontrolled.
This allows direction steps outside the small-radical range in
(TC4). Likewise (PCP4) prices the total remaining quotients in
the original record. The short index in (DL1) is the parameter
along a line of that record. Large actual box area relative to
the divisor product can leave the line parameters too long
for this proof. For example, on the formal balanced
\(\tfrac13\mathrm{Unif}[2,4]\) spectra and at raw total
record price \(f<1/3\), the divisor weight is at most
\(4f<4/3\), before any discount or ratio narrowing. In
height-one coordinate boxes this is short of the numerical
weight two within the (DL2) margin. The combined degree density in this model
is one half on the interval from two to four. Taking its top
mass f sharpens the numerical bound to 4f minus f squared,
which is at most 11/9 for f at most 1/3. At that raw price a
ratio interval would need numerical width at most X to the
power minus 7/9, up to a vanishing height exponent, just to
enter the lattice range. Extracting an inexpensive record or ratio window on arbitrary
balanced profiles remains open.
[Unit blocks at a large small-radical step (TC4)](ABC_TARGET_PACKET_CRT_MOMENTS.md#unit-progression-cells-with-a-large-small-radical-step)
· [Fixed divisor records and short quotients (PCP4)](ABC_SHORT_PROGRESSION_PAIR_COUNTS.md#arbitrary-span-pair-counts-and-fixed-full-divisors)

## DL4. The gcd of the actual projective minors

Keep the height caps, the three radical caps and the homogeneous
box of (DL2), with area parameter
\(\mathcal B=|\det\mathcal S|^{-1}\). Let \(\mathcal F\) be any
specified family of \(M\ge2\) distinct eligible primitive rows
in this box. Define
\[
        G_{\mathcal F}
          =\gcd_{i<j}|A_iB_j-A_jB_i|.
\]
Replacing \((A,B)\) by either other pair of roles changes only
the signs of these minors. At every sufficiently small fixed
\(L\), for every fixed \(\eta,\varepsilon>0\) with the budget
cutoff chosen in terms of \(\eta\) and the fixed caps,
\[
            M\ll n^\eta X^\varepsilon
                         \sqrt{\mathcal B/G_{\mathcal F}} .
                                                               \tag{DL4}
\]
Equivalently \(G_{\mathcal F}\ll
\mathcal B M^{-2} n^\eta X^\varepsilon\), on resizing the
allowances. The determinant may be computed inside a fixed
divisor record or on the whole specified family.

Two distinct positive primitive rows are noncollinear.
Thus the lattice generated by the actual complete coordinate
pairs has rank two. The determinant of a lattice basis is a
linear combination of the pair minors, since both basis vectors
are integral combinations of these pairs. Conversely that
determinant divides every pair minor. The determinant is
therefore \(G_{\mathcal F}\). In a divisor record this is
a multiple of \(D_A D_B D_C\). Using the generated lattice
for one specified family only uses the original record. Its
basis is fixed in that family, so the price stays that
of the original record.

Put \(\Delta=G_{\mathcal F}/\mathcal B\). This is at most two
by taking any pair of rows in the transformed unit box.
If \(\Delta\ge n^{-L^{-\xi_1}}\), apply (DL2) to this
generated lattice and absorb the term one into
\(O(\Delta^{-1/2})\). If \(\Delta<n^{-L^{-\xi_1}}\), use
the fixed-support count on two complete roles. Their bounded
radical caps give \(M\ll n^{2a_0}X^\varepsilon\).
Since \(\Delta^{-1/2}>n^{L^{-\xi_1}/2}\), this gives
(DL4) at sufficiently small budgets as well. This last
step makes the determinant threshold unnecessary for the
inequality on the gcd itself.
[Uniform fixed-support counting](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

In a fixed coordinate box in which the smaller two complete
upper caps are \(X_U,X_V\), (DL4) shows a record of exponent
\(r>0\) forces \(G_{\mathcal F}\le X_U X_V n^{-2r+o(1)}\)
in the height-first small-budget convention. In an interval of
ratios of width \(d\), the corresponding area is \(dX^2\).
Classical projective spacing in that interval gives
\(M\le 1+dX^2/G_{\mathcal F}\), with one power of \(M\) in the
resulting bound for the gcd. At the central scale
\(X_A,X_B,X_C\asymp X\), the geometry term in BBLT Proposition
4.1, with all three divisor records fixed, is
\[
  1+\frac{(X_A/D_A)(X_B/D_B)(X_C/D_C)}
           {\max(X_A,X_B,X_C)}
       =1+\frac{X_U X_V}{D_A D_B D_C}.
\]
For the matched record with \(D_A D_B D_C\asymp X^2 n^{-v}\)
and \(v>0\) fixed, this gives exponent \(v\). The bound with
the three radical caps has exponent \(v/2+\eta\). The
comparison concerns the per-record geometry term. Applying
it to an unrestricted family still requires inexpensive
records of the indicated determinant sizes.
[Projective spacing for common-lock lattices (SL5)](ABC_SPLIT_CUBIC_LOCK_LATTICES.md#the-exact-span-of-the-full-lock-primitive-locus)
· [BBLT, Proposition 4.1](https://arxiv.org/html/2410.12234v2#S4)

## DL5. A charged divisor-record inequality

Fix a box of actual complete heights and let
\(h_1,h_2\) be the two smallest upper height exponents.
For every prime \(p\le X^H\) assign \(q_p\in(0,1]\),
and write
\[
   B=\sum_{p\le X^H}q_p,\qquad
   C(P)=-\log_X\prod_{p\mid P}q_p,\qquad w=\log_X P,
   \qquad P=D_A D_B D_C .
\]
Keep the original three radical caps. The number of rows
admitting a choice of divisors with
\[
                         C(P)+\tfrac12(h_1+h_2-w)_+\le b
                                                               \tag{DL5}
\]
is at most \(O(X^{b+\eta L+\varepsilon}e^B)\), for fixed
\(\eta,\varepsilon>0\) at sufficiently small fixed budgets.
One may also use a prescribed homogeneous ratio window
and its area exponent in place of \(h_1+h_2\). A collection
of window records must be charged separately.

Indeed, for charge cap \(b_0\) each allowed prime support
contributes at least \(X^{-b_0}\) to
\(\prod_{p\le X^H}(1+q_p)\le e^B\).
Allocation of the primes to three roles and enumeration
of their powers up to \(X^H\) cost \(X^{o_X(1)}\)
uniformly at fixed support. Thus there are at most
\(X^{b_0+o_X(1)}e^B\) divisor records of that charge.
Bin the nonnegative numerical deficit
\((h_1+h_2-w)_+\) in intervals of width \(1/\log X\).
For a bin of deficit \(v\) the charge is at most
\(b-v/2+O(1/\log X)\). The bound in (DL2), or the
trivial two-letter cap when outside its determinant range
as in (DL4), contributes \(X^{v/2+\eta L+\varepsilon}\)
after the usual allowance refinements. The product gives
(DL5) after summing the logarithmic number of bins.
This improves the first-power deficit of (SP6) to a half
deficit under the three original radical caps.
[Charged angular packing (SP6)](ABC_SMOOTH_SUPPORT_DIVISOR_PACKING.md#charged-radicals-and-angular-packing)

For an actual family of size at least \(n^d\) let \(f_p\)
be its frequency of prime \(p\) in the radical of a row.
Taking \(q_p=\max(f_p,p^{-1})\) gives \(B=o(\log X)\),
as in (SP6). Consequently for every fixed \(\sigma>0\),
outside a power-smaller set of rows, every choice of the
three divisors in a fixed complete-height box has
\[
               \frac{C(P)}{L}
                +\frac{(h_1+h_2-w)_+}{2L}>d-\sigma .
\]
Here the allowances are chosen smaller than \(\sigma\)
before the fixed budget and its height limit. Dyadic
complete-height boxes cost subpower. The inequality
quantifies the determinant or record price required at
a saturated count. Finding an affordable divisor
selection that violates it over the unrestricted
balanced profiles remains open.
