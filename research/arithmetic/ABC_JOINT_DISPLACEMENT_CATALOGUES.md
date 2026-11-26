# Joint displacement bounds on genuine three-cap charts

An exact pair of complete displacements restricts all
three original dictionaries at once. Full divisors turn
these restrictions into a stronger row-count bound. In the
high-depth sector, even the pairs that share all three
registers can already exceed this upper bound. The
argument uses the additive equation and the actual
radical caps on both rows.

## Exact displacement fibres

**JD1.** Let \(\mathcal F\) be a chart of distinct positive
primitive triples \(A+B=C\le X\), with fixed
\(D_T\mid T\) in every role. Missing locks are set to
one. By primitivity the three \(D_T\) are pairwise coprime.
Write \(D=D_A D_B D_C\), \(n=X^L\), and
\[
                   \theta=\log_n(1+X/D).
\]
Take actual pairs with fixed exact signed displacement
\((h_A,h_B,h_C)\), where \(h_A+h_B=h_C\) and all three
entries are nonzero. Let the first and second rows
satisfy the original caps imposed on the chart.
Suppose the pairs of complete values in role \(T\)
with difference \(h_T\) have count
\(O(n^{d_T+\eta}X^\varepsilon)\) for arbitrary fixed
allowances. Let \(n^{u_T}\) bound the maximum degree
of a complete \(T\) in the original chart. Then the
number of ordered rows realizing the displacement
costs at most
\[
  O\left(n^{F+\eta}X^\varepsilon\right),\qquad
  F=\min\left\{\min_T(d_T+u_T),\
      \min_{T\ne T'}(d_T+d_{T'}),\
      \frac{d_A+d_B+d_C+\theta}{2}\right\}.
                                                               \tag{JD1}
\]
The \(d_T\) may be actual uniform upper caps or
any of the fixed-shift estimates below. The statement
is uniform in the exact displacement.

The first row determines the second one, so
the \(T\) dictionary in this fibre has size at
most its shifted value-pair count. Fixing any
complete \(T\), or any two letters, gives
the first two bounds. For the third use the
following finite argument. Suppose the first-row
set has size \(E\) and complete dictionaries
of sizes at most \(K_A,K_B,K_C\).
For each role delete rows at vertices of degree
less than \(E/(12K_T)\) in the current set.
Iterate until no such vertex remains. Each
deleted vertex is charged once and the total
deletion is at most \(E/4\).
On the surviving set, every active residue
\(A=-B\bmod D_C\) has at least
\(E/(12K_A)\) distinct \(B\)-values, since
any pair of complete \(A,B\) fixes a row.
The classes are disjoint, so their count
is at most \(12K_AK_B/E\).
Apply the same argument on the other two
locks, using signed residue compatibility at
each. A triple of those registers fixes both
\(A,B\bmod D\) by CRT and thus has at most
\((1+X/D)^2\) actual lifts. Multiplication
gives \(E^4\ll(K_A K_B K_C)^2(1+X/D)^2\).
Together with the fixed-shift projection caps
this proves JD1, without requiring regularity
of the translated subfamily beforehand.

### Fixed-shift inputs and the typical half prices

Outside any prescribed polynomial-height pools,
let \(r_T\) cap the radical mass of \(T/D_T\)
in both rows. To apply the arbitrary-shift
count NS1, divide out \(D_T\), and use the exact
valuation layers with fixed singleton dictionaries
after binning and fixing any exceptional tails.
If \(b_{T,1},b_{T,2}\) bound the corresponding
low-degree masses at a fixed positive cutoff in
scaled degree, one may take
\[
 d_T\le
   \min\{r_T,\ \max(2r_T/3,\min(b_{T,1},b_{T,2})),\
              \max(r_T/2,b_{T,1}+b_{T,2})\}
\]
up to the fixed small allowances. When there
is a **different-role** deep full lock \(D_{T'}\)
and \(D_{T'}\mid h_T\), FS1 gives \(d_T=r_T/2\)
regardless of the low-degree mass. Indeed
\(D_{T'}\) is prime to both \(T\)-values. Divide
by \(D_T\) and then by the gcd of the two
cofactors. The gcd divides \(h_T/D_T\), has
divisor-many choices, and cannot remove any part
of \(D_{T'}\). The fixed shift consequently
retains the full high-degree divisor.
[Nonprimitive shift estimate (NS1)](ABC_NONPRIMITIVE_FIXED_SHIFT_ENERGY.md#two-dictionary-bounds-for-any-nonzero-shift)
· [Fixed shifts with a deep divisor (FS1)](ABC_FULL_SHIFT_DICTIONARIES_AND_CODEGREES.md#a-fixed-shift-count-without-a-radical-cap-on-the-shift)

**JD2.** On original balanced shells satisfying
CF1's individual cap tolerances, the same
half prices hold on all but a power-small
fraction of all original pairs. More precisely,
for a chart with size \(N\), assume cap
\(n^{1/3+\kappa_1}\) per original letter
with the fixed rank-one CF1 tolerance
\(\kappa_1\). There are absolute
\(\sigma,h_0>0\) such that at most
\(O(N n^{1/3-\sigma}X^\varepsilon)\)
ordered pairs have any reduced role quotient
of height at most \(X^{h_0}\).
For a fixed first row and role \(T\) enumerate
\(g=\gcd(T,T')\mid T\) at divisor cost.
If the reduced height is at most \(X^{h_0}\),
the second value is \(g z\) with \(1\le
z\le X^{h_0}\). This complete value has an
actual affine feature of rank one, so CF1
applies uniformly after taking \(h_0\) below
half its determinant threshold. Sum over
the first rows, the divisor choices and the
three roles. For \(N\ge n^{1/3-\kappa}\)
with small fixed \(\kappa<\sigma\) this
discards only a power-small fraction.
When the chosen pair sector is sparse,
compare its lower exponent against this
exceptional bound before discarding.

For each fixed nonzero role shift, divide the
two complete \(T\)'s by their gcd, choosing
that gcd amongst the divisors of the shift.
The resulting signed primitive equation is
a positive primitive additive triple after
ordering its three absolute terms. Its fixed
shift is a prescribed complete coefficient,
and its height exceeds \(X^{h_0}\) on the
remaining pairs. Strip prescribed supports
from both moving terms at height-subpower
cost, including the support of \(D_T\).
Apply MCB2 with a fixed third coefficient
and exact-order singleton root slots on
the moving terms. This supplies \(d_T=r_T/2\)
uniformly, with no cutoff on low orders.
The same height floor holds because only
the normalization gcd alters the complete
quotient values. This argument is for
pairs in the original balanced chart, and
does not assert the tall-quotient hypothesis
for all pairs in an arbitrary sparse sector.
[Small affine feature count (CF1)](ABC_COMPLETE_LETTER_SMALL_FEATURES.md#a-strict-saving-for-small-affine-minors)
· [Arbitrary coefficients with a height floor (MCB2)](ABC_MULTIPLE_SHORT_COFACTOR_BOUNDARY.md#bounded-group-dictionaries-with-arbitrary-fixed-coefficients)

## Actual shift catalogues from the signed registers

**JD3.** If selected nonzero-displacement pairs in
a near-full chart have count \(n^{p-o(1)}\) and
their actual exact vector catalogue has price \(q_*\),
then \(p\le q_*+F+o(1)\). Here \(F\) is from JD1
on the selected set, and \(o(1)\) represents the
stated height-first fixed allowances. The following
different input prices are available.

For \(I\) any subset of the roles, let \(q_I\) price
the **joint** difference vectors of the selected
signed registers \(B\bmod D_A,A\bmod D_B,A\bmod D_C\)
on \(I\). The shift in role \(T\) is then determined
modulo \(M_T=D_T\prod_{W\in I\setminus\{T\}}D_W\).
Indeed the own-lock residue is zero, and at any
other lock the original additive equation gives
the displacement from its signed register difference.
Any two shifts determine the third, so
\[
 q_*\le q_I+
     \min_{T\ne T'}\left[
       \log_n(1+X/M_T)+\log_n(1+X/M_{T'})\right]+o(1).
                                                               \tag{JD3}
\]
Thus full register differences give
\(q_*\le q_{\{A,B,C\}}+2\theta+o(1)\).
For one register \(C\) the two usable moduli
include \(D_A D_C,D_B D_C\).
For any chosen single nonzero \(h_T\) one
can instead enumerate only this shift,
then use \(d_T\) for the value-pair count
and \(2u_T\) for the choices of original rows.
That comparison can be stronger when one
complete dictionary has low multiplicities.

Without a register record one has the
alternative
\[
                    q_*\le\log_n(1+X^2/D)+o(1).
\]
The complete positions lie in the integer lattice
\(D_A\mid A,\ D_B\mid B,\ D_C\mid A+B\),
of index \(D\) by CRT. If there are two original
rows their position vectors are independent,
because a primitive row is determined by its
positive rational proportion \(A/B\).
Take a shortest nonzero lattice vector \(v\)
in sup norm, and one independent original row.
Their determinant has modulus at least \(D\)
and at most \(2X\|v\|_\infty\).
The lines parallel to \(v\) that meet
\([-X,X]^2\) are indexed by their
determinants with \(v\), in \(D\mathbb Z\).
There are \(O(1+X\|v\|_\infty/D)\) such
lines. Each contains \(O(1+X/\|v\|_\infty)\)
lattice points because \(v\) is shortest.
Using \(D/(2X)\le\|v\|_\infty\le X\)
gives \(O(1+X^2/D)\) positions in the
difference box, with no choice of a residue
catalogue. The case of at most one row is
immediate.

One may also price \(q_*\) from a catalogue
of primitive shift directions. Divide the
three signed displacements by their common
positive gcd \(g\). For all nonzero coordinates
the absolute values form an actual primitive
additive triple. If its height has a fixed
positive numerical floor and its outside
radical cap is \(n^\beta\), EB4 yields
a direction catalogue of exponent
\(\beta/3+\eta\). Fix the height bins
first, as in NGT3, to obtain the height-first
count. For a prescribed direction \(v\),
the scale is a multiple of
\[
             g_0(v)=\prod_T D_T/\gcd(D_T,v_T).
\]
Thus if \(g/g_0(v)\) has radical at most
\(n^\xi\) outside its prescribed pool, or
has numeric size at most \(n^\mu\),
the scale adds at most
\(\min(\xi,\mu)\) when both caps apply.
This prices the actual scale rather than
summing unrestricted common multipliers.
[Effective radical budgets (EB4)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)
· [Height bins for primitive directions (NGT3)](ABC_SECANT_NORMALIZED_GCD_TAIL.md#counting-complete-primitive-directions-below-a-fixed-radical-cap)

At balanced caps, on full-scale good pairs
with three tight locks \(D\ge X\), JD1
and JD2 give
\[
         F\le\min(1/6+u,1/4)+o(1).
\]
Here \(u\) is the common complete multiplicity
exponent. Thus the joint difference catalogue
of the three registers must have price at
least \(\max(1/2-u,5/12)-o(1)\).
In particular the pointwise bound using
an original complete degree is stronger
near occurrence, while the simultaneous
shifted fibre is stronger at large
multiplicity. A fixed discount in the
catalogue below the displayed threshold
contradicts near saturation on a sufficiently
cheap chart. For one tight pair of locks and
one actual additive difference list the
single-coordinate comparison can instead
be used, in particular at occurrence.

## A deep same-register sector without additive compression

**JD4.** Assume the three locks on the chart
have positive numerical weight and a fixed
positive scaled fullness floor. Let the actual
row-count exponent be \(d\), and let \(t\)
cap the price of the joint signed register
support. If \(d>t+\theta\) by a fixed
margin, then
\[
                 2d-t
                 \ \le\ \frac{r_A+r_B+r_C}{4}
                             +\frac52\theta+o(1).
                                                               \tag{JD4}
\]
Indeed Cauchy–Schwarz gives at least
\(n^{2d-t-o(1)}\) pairs with all three same
registers. For a fixed first row in such
a cell, a zero complete shift leaves one
other original value determined modulo \(D\),
so there are at most \(3(1+X/D)\) such
second rows. These pairs are negligible
under the displayed margin. All other
displacements are multiples of \(D\) in
each coordinate. Their actual vector
catalogue therefore costs at most
\(2\theta+o(1)\). The shift in any role
retains the deep divisor on either different
role after dividing out the gcd of the two
moving terms, as in the FS1 arm. Thus
\(d_T=r_T/2\) in JD1 even when some of the
shifted roles contain many low exponents.
Use the last bound in JD1 and sum the
exact displacements to obtain JD4.

For instance with near-full balanced exponent
\(1/3\), default individual raw caps \(1/3\),
and \(\alpha\in(1/6,2/9)\), the actual support
is at most \(n^{3s+o(1)}\) by the single
register caps, where \(s=2\alpha-1/3\).
Write \(h=2/3-3\alpha>0\). Ordinary height
packing requires \(\theta\ge h-o(1)\),
since the \(n^{3s+o(1)}\) cells lift
at most \(O(n^{2\theta})\) rows each.
JD4 excludes the wider open band
\[
                 \theta < \min\{2h,\ 4h/5+1/30\}
\]
by fixed margins. In particular the
elementary maximal-depth boundary
\(\theta=h+o(1)\) is excluded throughout
this open range of \(\alpha\). The conclusion
also applies at more general dictionary
exponents whenever the **actual joint** support
has price \(t=d-2\theta+o(1)\), with
\(d=1/3-o(1)\) and \(0<\theta<1/6\).
Use the raw caps after the fixed support
choices and include every actual chart
price in \(d\). The zero-surplus endpoint
\(\alpha=1/6,\theta=1/6\) gives equality
at the default cap. An occurrence chart
can have \(t=d\) and therefore no
off-diagonal same-register pairs.

### One shared register and conditional difference prices

If only the \(C\)-register is required equal
for the selected pairs, and its lock is
deep, then \(d_A,d_B\le1/6+o(1)\) at
default caps by FS1 since both shifts
retain \(D_C\). One may use the full
complete dictionary exponent \(\alpha\)
for \(d_C\), with no further hypothesis.
If the actual *shared* difference
catalogue of the other two registers,
over all selected same-\(C\)-register cells,
has price \(q_{\mathrm{cond}}\), then the
pair exponent is at most
\[
 q_{\mathrm{cond}}+\frac16+\frac{\alpha}{2}
                  +\frac52\theta+o(1).
\]
The two conditioned register differences
and the zero \(C\)-register difference
supply the hypotheses of JD3. In an RB4
chart of fixed \(u>0\), its actual same-cell
pairs with all shifts nonzero have exponent
\(s+4u=1/3+2u\), since fixing one complete
letter inside the cell leaves at most
\(n^{u+o(1)}\) choices there.
Thus \(q_{\mathrm{cond}}\le2u+o(1)\)
already gives a strict saving when
\(D\ge X\) and \(u>0\), with the
actual chart and regularization allowances
inside the gap \(u/2\). It is not enough
to have small difference lists separately
in each cell, since those lists may
vary with the cell. This is a precise
remaining catalogue obligation.
[Regular genuine residue cells (RB4)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#the-actual-residue-cell-sizes)
· [Tight locks at fixed shifts (NJ5)](ABC_NORMALIZED_JOINT_MODULAR_MOMENTS.md#tight-full-locks-and-the-third-register)
