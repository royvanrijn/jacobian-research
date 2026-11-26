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
radical cap is \(n^\beta\), fix the actual
powers on the prescribed pool at subpower
cost. MCB2 on the other prime powers costs
at most half the two cheapest direction
masses, hence yields a catalogue of
exponent \(\beta/3+\eta\). The fixed
positive height floor suffices for MCB2. For a prescribed direction \(v\),
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

### CF1 on projected deep registers

**JD6.** Work in the original shell \(X/2<C\le X\), and impose
\(\operatorname{rad}(T)\le n^{1/3+\kappa_1}\) on all three letters,
where \(\theta_1,\delta_1,\kappa_1,L_1>0\) are the rank-one CF1
constants. On a chart choose a complete role \(T\) and a subset \(I\)
of the other roles' signed registers. Fixing these registers determines
one residue of \(T\) modulo \(M_T=D_T\prod_{W\in I}D_W\), by JD3 and CRT.
Let \(H_T=\max T-\min T\) on the chart, and suppose \(M_T\ge H_T X^{-\theta_1/2}\). On the
cell \(T=T_0+M_T z\), all differences of the integer feature \(z\)
have magnitude at most \(H_T/M_T\le X^{\theta_1/2}\). Thus CF1 counts
\(O_{L,\varepsilon}(n^{1/3-\delta_1}X^\varepsilon)\) rows there, uniformly
in the modulus and register choices. This proof uses only divisibility
of the fixed pairwise coprime \(D_T\). It has no scaled fullness
requirement.

If the charts cost \(b_0\) and the joint support of the selected
projection on each chart costs \(t_I\), the exponent is at most
\(1/3-\delta_1+b_0+t_I\), with the usual height allowances. Use the best
available \(T,I\) independently on each chart and include the bounded
choices of these roles and subsets. In particular \(D\ge X/n^K\) for
fixed \(K\) suffices for both foreign registers, for
\(0<L<L_1\) with \(KL<\theta_1/2\). There is no cost for the own register
since \(T=0\bmod D_T\).

In the balanced regular dictionary setup with foreign single-register
price \(s=2\alpha-1/3\), this prices the two registers at \(2s\), or at
less if their actual joint projection is cheaper. The minimum
\(\alpha=1/6,\theta=1/6\) left at equality in JD4 is therefore excluded
on charts of negligible price. JD5 also excludes that limiting case
through the joint register information bound. JD6 applies directly
to the actual finite charts and can use the diameter of a shorter
role. With a fixed margin the two-register
version also excludes \(\alpha<1/6+(\delta_1-b_0)/4\) in the same
depth range. If just one foreign register yields sufficient numerical
depth, replace the divisor \(4\) in this condition by \(2\).
[Rank-one small-feature count (CF1)](ABC_COMPLETE_LETTER_SMALL_FEATURES.md#a-strict-saving-for-small-affine-minors)

### Two short CRT quotient coordinates give arbitrarily small count

**JD12.** Fix \(B,C,\nu,\eta>0\). Use original primitive shell rows
with all three caps \(\operatorname{rad}(T)\le n^C\) and fixed
pairwise coprime role divisors \(D_R\mid R\). Suppose at least one
\(D_R\ge X^\nu\). A projection of the signed register tuple fixes
the residues of any two distinct complete roles \(T,W\) modulo
\[
 M_Q=D_Q\prod_{R\in I\setminus\{Q\}} m_R
 \quad(Q=T,W),
\]
where \(I\) is the set of fixed native registers and
\(m_R\mid D_R\) is the modulus of that register projection.
Let \(H_Q\) be the complete-value diameter of that role on the chart,
or an upper bound on the diameter in the relevant projected cell.
If \(H_T/M_T,H_W/M_W\le n^B\), then for sufficiently small **fixed**
\(L\) every projected cell has
\[
 O_{B,C,\nu,\eta,L,\varepsilon}(n^\eta X^\varepsilon)
\]
rows. More generally the same conclusion holds with the two quotients
bounded by \(X^\theta\) for sufficiently small \(\theta>0\) chosen
in terms of \(C,\nu,\eta\) before \(L\).

Choose one actual row in the cell as anchor. Each of the two
anchored complete differences is divisible by its respective \(M_Q\).
The two integer quotients are exact affine features for the two
complete roles. The original additive relation gives the third as
an affine function of these same features. Apply the exact
small-feature count with one prescribed tall original divisor.
The coefficients and both CRT moduli may vary with the cell.
This argument includes rows in an affine line or a point.

If the chart catalogue costs \(b_0\) and the actual joint projection
support costs \(t\) uniformly per chart, the exponent is at most
\(b_0+t+\eta\), with ordinary height allowances. In particular, for
\(D_AD_BD_C\ge X/n^K\) with \(K\) fixed, use the full joint tuple
and any two roles. Both quotients are bounded by \(n^K\), and for
sufficiently small fixed \(L\) some divisor is at least \(X^{1/4}\).
This improves the per-cell rank-one CF1 bound and gives a strict
saving whenever \(b_0+t<1/3\) by a fixed margin. The depth
condition allows a mixed projection with a smaller modulus on a
shorter role. For example, in a balanced chart with \(\alpha<1/3\),
a single native register of price \(s=2\alpha-1/3\) suffices if the
two opposite roles both meet the short-quotient test and the chart
price is below \(1/3-s\) by a fixed margin. Neither this test nor
the full joint depth is an extraction result for arbitrary rows.

For a fixed \(K>0\), this finite-depth result subsumes the JD4
trigger at sufficiently small fixed budgets whenever the chart
has exponent \(d>t+\theta\) by a fixed margin and
\(D\ge X/n^\theta\), \(\theta\le K\). In fact it rules out
\(d>t\) by any fixed margin throughout that depth range. Under the
three caps \(n^C\), the elementary two-letter support count gives
\(d\le 2C+o(1)\). Thus every JD4 trigger with a fixed gap can use
one fixed \(K>2C\) and JD12, taking the cell allowance \(\eta\)
below the gap. At the default balanced register prices
\(s=2\alpha-1/3\), the full tuple costs at most \(3s\), so the
total count exponent is at most \(b_0+6\alpha-1+\eta\).
In particular any fixed \(\alpha<2/9-b_0/6\) is excluded in this
finite-depth range. Actual joint support can improve this
criterion. JD6 can use relative depth in just one role or a
cheaper projection for which the second quotient remains large.

The same arbitrarily small bound holds on the **whole chart**
without charging the projected residue support when the two
least-residue lists, each translated by its own chart minimum,
have fixed-dimensional exact affine integer models with feature
coordinates bounded by \(n^{B'}\) for fixed \(B'\). For each role,
adjoin the integer quotient by \(M_Q\) between the complete value
translated by the chart minimum and its least residue. It is
bounded by \(1+H_Q/M_Q\). Concatenate both sets of features, so
the three original letters are exact affine functions in a
dimension fixed independently of height and budget. Apply the
same tall-divisor result. Bounded doubling on both projected
residue lists is one sufficient condition for these models, by
the lift in JD7.

[Exact small features with one tall original divisor](ABC_AFFORDABLE_FIBRE_AFFINE_RANK.md#exact-small-features-with-one-fixed-tall-divisor)

### Lifting one affine complete letter using a partial register

**JD13.** Work on an original primitive chart with all three
radical caps \(n^C\) and one fixed role divisor at least \(X^\nu\).
Suppose the complete role \(T\) is affine in integer features of
fixed ambient dimension, with actual affine rank \(r\) and all
maximal anchored feature minors of size at most \(Q_0\). The other
roles are \(U,V\). Fix their role divisors \(D_U,D_V\) and the
actual native signed register on \(T\) modulo \(m_T\mid D_T\).
Let \(H_V\) be the complete \(V\)-diameter of this cell. If
\[
         D_U D_V m_T>(r+1)H_V Q_0 ,
\]
then all three complete letters are affine functions on the
same feature space. If \(Q_0\le X^\theta\) for sufficiently small
\(\theta>0\), this cell has \(O(n^\eta X^\varepsilon)\) rows for
every fixed \(\eta>0\) and sufficiently small fixed \(L\).

Project onto \(r\) feature coordinates giving an injective affine
map on the span. For any anchor, augment the \(r\) coordinate
difference columns by the \(V\)-difference column. Each
\((r+1)\)-minor is divisible by \(D_V m_T\), since all the
\(V\)-differences vanish modulo those coprime moduli. It is also
divisible by \(D_U\): modulo \(D_U\) its last column equals the
signed \(T\)-difference column, and the corresponding integer
determinant vanishes exactly because \(T\) is affine. Expansion
along the last column bounds the absolute determinant by
\((r+1)H_VQ_0\). All augmented minors vanish and therefore \(V\)
is affine, including at rank zero. The additive relation gives
\(U\). The integral Cramer coordinates of CF2 have absolute
values at most \(Q_0\) and preserve all these identities. The
fixed-tall-divisor theorem supplies the asserted count.

Sum the actual incoming chart and register costs, as in JD12.
The numerical test involves the diameter of one cell, permitting
a short displacement range even when its absolute \(V\)-height
is large. It works with unknown real affine coefficients for \(T\).

[Integer Cramer coordinates (CF2)](ABC_COMPLETE_LETTER_SMALL_FEATURES.md#normalization-and-a-rank-drop-criterion)
· [All-letter exact features with a tall divisor](ABC_AFFORDABLE_FIBRE_AFFINE_RANK.md#exact-small-features-with-one-fixed-tall-divisor)

### Additive structure of projected registers

**JD7.** Fix \(B,K\ge1\). On the same kind of chart let \(\mathcal R_T\)
be the actual residue set of \(T\bmod M_T\) in the selected projection.
Suppose \(H_T/M_T\le n^B\), \(|\mathcal R_T|\le n^B\), and
\[
             |\mathcal R_T+\mathcal R_T|
                       \le K|\mathcal R_T|\quad\text{in }\mathbb Z/M_T\mathbb Z.
\]
Then there exist \(\delta,\kappa,L_0>0\), depending only on \(B,K\), such
that at every fixed \(L<L_0\) a chart satisfying the three original
caps \(n^{1/3+\kappa}\) has count
\(O_{B,K,L,\varepsilon}(n^{1/3-\delta}X^\varepsilon)\). No residue-support
price is charged in this bound. Summing charts of price \(b_0\)
gives a strict saving if \(b_0<\delta\) by a fixed margin.

Translate \(T\) by its chart minimum; this translates the residue set and preserves doubling. Take the least nonnegative integer representatives of the translated set.
Their ordinary sumset has size at most twice the modular sumset size,
because every residue has at most two lifts in \([0,2M_T-2]\).
The integer Freiman theorem places the representatives in a proper
generalized arithmetic progression of rank at most \(r(2K)\) and size
at most \(C(2K)|\mathcal R_T|\). Its integer coordinates can all be
taken of size at most \(C(2K)n^B\), including in the singleton case.
Append the quotient of the translated value minus its representative, divided by \(M_T\), as an extra
coordinate, of size at most \(n^B\). Then \(T\) is one fixed integer
affine function of these coordinates. Use CF1 with dimension
\(R=r(2K)+1\); each maximal difference minor is
\(O_{K,R}(n^{BR})\) by its determinant expansion. Choose \(L\) with
\(BRL<\theta_R\) and \(L<L_R\). This proves the claim with \(\kappa_R\)
and \(\delta_R\), even when many triples repeat the same residue.

It also suffices that each of the foreign registers used has bounded
doubling in its own cyclic group and size at most \(n^B\). Model
the integer representatives of each by its own proper progression
as above. Compose their affine models with the integer linear CRT
formula for the projected residue. Reduce the constant and each
progression direction modulo \(M_T\). The representative produced
as an affine function of the progression coordinates differs from
\(T-\min T\) by a multiple of \(M_T\), whose coefficient has
magnitude \(O_{K,R}(n^B)\), by the bounds on the reduced directions,
the coordinates, and \(H_T/M_T\). This extra integer coordinate
puts the complete \(T\) in the required fixed-rank model. No
independence assumption on the foreign registers is needed.

More generally the same proof works whenever the least representatives
admit a bounded-rank affine feature model with coordinates \(n^{O(1)}\).
A fixed doubling bound is one sufficient structural condition. A
doubling ratio of \(n^{o(1)}\) can still grow with height, and the
argument does not give a rank bound uniform in height in that
case. Extraction of enough affordable charts with either JD6 or JD7
is open across the complete balanced dictionary interval.
[Green–Ruzsa, Theorem 1.1](https://arxiv.org/pdf/math/0505198)

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
inside the gap \(u/2\). The estimate requires a common catalogue across
cells, since separate lists may vary with the
cell. Extracting that common catalogue remains
open.
[Regular genuine residue cells (RB4)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#the-actual-residue-cell-sizes)
· [Tight locks at fixed shifts (NJ5)](ABC_NORMALIZED_JOINT_MODULAR_MOMENTS.md#tight-full-locks-and-the-third-register)

### A rank-two affine index bound

**JD8.** In a genuine chart write its positions as \(P=(A,B)\).
Let its width in affine determinants be
\[
 {\cal A}=\max_{P_0,P_1,P_2}
             |\det(P_1-P_0,P_2-P_0)|.
\]
Suppose after fixing a chosen tuple of actual partial
signed registers \(J_R\bmod m_R\) with \(m_R\mid D_R\),
the remaining chart has rank two. Its position
differences lie in the exact congruence lattice of index
\[
                       I=D_AD_BD_C\prod_R m_R.
\]
Indeed a given role lock first imposes one independent
divisibility condition at its primes. The other signed
coordinate fixed modulo \(m_R\) gives its additional
index \(m_R\). This is immediate in the \(A,B\) basis
for the \(A\) or \(B\) locks and in the \(A,A+B\) basis
for the \(C\) lock. Coprimality across roles gives CRT.
The actual index of the lattice spanned by all affine
differences can only be larger. In rank two it is the
gcd of the nonzero affine determinants, by the integer
Smith normal form.

Choose one occupied position \(P_0\) and an integer
basis of any such containing lattice. Write \(P-P_0\)
in that basis as integer features \(z\in\mathbb Z^2\).
The complete letters are fixed affine functions of
\(z\). Every maximal rank-two affine minor of these
features equals its original affine determinant divided
by \(I\). Thus if \({\cal A}/I\le X^{\theta_2}\), where
\(\theta_2\) is the CF1 rank-two constant, CF1 supplies
its strict row saving. The feature minors suffice even with large individual coordinate sizes. If the affine
rank is one, a coordinate line is handled by the fixed
complete-letter bound FL1 and a noncoordinate line
with nonzero intercept by AS1. A line through the origin
contains at most one positive primitive row. Thus
degenerate rank does not obstruct the strict chart count.

Uniformity includes all three original caps, with
tolerances at most the CF1 rank-two tolerance and the
AS1/FL1 gap. Sum over the actual incoming chart and
partial-register records; their total exponent must
be below the applicable strict saving. No fixed-order
or fullness hypothesis on the locks is used. In
particular, if \(D_AD_BD_C\ge X^{2-\theta_2/2}\), the
bound is unconditional within that actual lock chart
for all balanced dictionary exponents, since
\({\cal A}=O(X^2)\) and the unused numerical margin
absorbs constants. With smaller \(D\), the affine
determinant condition also covers a short transverse
strip with both complete coordinate diameters large,
where the rank-one JD6 diameter test may give nothing.
Near-saturation at occurrence must retain
\({\cal A}/I>X^{\theta_2}\) on every chart having
sufficiently small actual record price.

[Small feature minors (CF1)](ABC_COMPLETE_LETTER_SMALL_FEATURES.md#a-strict-saving-for-small-affine-minors)
· [Noncoordinate line count (AS1)](ABC_AFFINE_SLICE_RADICAL_COUNTS.md#noncoordinate-lines-and-projection-sizes)
· [Fixed-letter count (FL1)](ABC_FIXED_LETTER_RADICAL_CAPS.md#a-pointwise-count-with-prescribed-prime-pools)

### An arbitrarily small exponent from a heavy lock and affine area

Suppose the actual affine rank in JD8 is two and one fixed
original role divisor is at least \(X^\nu\) with \(\nu>0\)
fixed. For each \(\eta>0\), the criterion
\(\mathcal A/I\le X^\theta\), with sufficiently small
\(\theta,L\) chosen in terms of \(\eta,\nu\) and the raw cap,
gives count \(O(n^\eta X^\varepsilon)\) on the chart.

Take one actual maximum-area triangle with one vertex \(P_0\)
and use its two anchored differences as a real coordinate
basis. Every chart difference from \(P_0\) has coordinates
of magnitude at most one by Cramer's rule and the maximum
area property. The difference lattice now has covolume
\(I/\mathcal A\), and the two unit coordinate vectors are
in it. A shortest vector and a shortest nonzero orthogonal
projection lift give a basis of bounded lengths, after
subtracting an integer multiple of the shortest vector
from the lift. Hence all differences have integer basis
coordinates \(O(\mathcal A/I)\). The original complete
letters are exact affine functions of these two features,
so the fixed heavy divisor allows the stronger
[exact-feature count](ABC_AFFORDABLE_FIBRE_AFFINE_RANK.md#exact-small-features-with-one-fixed-tall-divisor).
This uses rank two; for arbitrary rank one the old JD8
line-count branch still applies.

### Modular energy from the original witness degrees

**JD9.** Let a regular balanced chart have \(n^{1/3-o(1)}\) rows,
common complete dictionary exponent \(\alpha\) and original
degree exponent \(u=1/3-\alpha\). Each complete \(T\)-value
on the pruned chart has at least \(n^{u-o(1)}\) witnesses.
Fix an original role \(T\) and an actual projection modulus
\(M=D_T\prod_{R\in I}D_R\) as in JD6, with \(I\) a subset
of the opposite registers. Write \(S\) for the occupied
\(T\bmod M\) residues and choose one distinct complete
\(T\)-representative per residue. Let \(H_T\) be the actual
diameter of complete \(T\)-values and
\[
                          v=\log_n(1+H_T/M).
\]
Let \(n^{r_T+o(1)}\) cap the radical of the two representatives
in role \(T\) after fixing the primes of any already
prescribed polynomial-height pool. This pool may contain
the support of \(D_T\), and the prescribed parts are fixed
at height-subpower cost when applying shift counts.

There exist fixed positive \(\sigma,h_0\), uniform in
the balanced chart, such that
\[
 E_{+,\mathbb Z/M\mathbb Z}(S)
  \ \ll |S|^2
   \left(n^{\,r_T/2+v+o(1)}
             + n^{\,\alpha-\sigma+o(1)}+1\right),
\]
with ordinary arbitrary small fixed count allowances.
The \(o(1)\) here follows the repository's inner height
and then outer small-budget convention.

For any nonzero exact difference \(t\), choose the gcd
of the two original representatives at divisor cost from
\(|t|\). If the reduced \(T\)-quotient has height at least
\(X^{h_0}\), MCB2 applied to the signed primitive shift
after gcd normalization gives at most
\(n^{r_T/2+o(1)}\) such ordered value pairs, as in JD2.
The third auxiliary term is the fixed exact shift, of
dictionary cost zero. The prescribed supported powers
of the two free values are fixed before applying that
count, and numerical height remains that of the
primitive quotient. Each cyclic difference has at
most \(O(1+H_T/M)=O(n^v)\) actual exact lifts. Summing
the squared good-pair multiplicities therefore costs
at most \(|S|^2n^{r_T/2+v+o(1)}\).

For the remaining pairs fix a first complete representative.
JD2, or its CF1 rank-one proof, gives at most
\(n^{1/3-\sigma+o(1)}\) original second ROWS whose complete
\(T\)-value has reduced quotient height below \(X^{h_0}\)
against it. This condition depends only on the two complete
\(T\)-values. Every eligible second representative therefore
contributes all its \(n^{u-o(1)}\) chart witnesses. The number
of bad representative pairs is at most
\(|S|n^{\alpha-\sigma+o(1)}\). Each cyclic difference has at
most \(|S|\) ordered residue representations, giving a
bad-energy contribution of at most
\(|S|^2n^{\alpha-\sigma+o(1)}\). The diagonal is harmless,
and mixed energy terms are bounded by splitting the squared
multiplicity using \((a+b+c)^2\le3(a^2+b^2+c^2)\).
This counts representative differences rather than all
pairs of original rows with the same value.

Since \(|S|\ge n^{\alpha-v-o(1)}\) by direct lift counting,
the energy is power-small relative to \(|S|^3\) when
\(v<\sigma\) and \(2v<\alpha-r_T/2\) by fixed margins.
In an actual JD5 near-tight chart at absolute depth,
\(|S|=n^{\alpha+o(1)}\) for the two-opposite-register
projection of any \(T\), and \(\alpha\ge2/9\). If the lift
exponent has \(v=0+o(1)\), this gives an absolute power gap
at least \(\min(\sigma,1/18)\) before allowances for every
role, including occurrence. At higher lift exponents use
the explicit two thresholds rather than assuming
\(|S|=n^\alpha\).

At occurrence there is also the direct OE2 transfer from exact
energy by Cauchy–Schwarz over difference lifts. The regular-witness
division above applies throughout the surviving dictionary range
with small lift exponent and it needs no role-specific positive-scale
mass condition.
[Occurrence energy bounds (OE1–OE2)](ABC_OCCURRENCE_DICTIONARY_ADDITIVE_ENERGY.md)

### Native register energy from actual pairs of original rows

**JD10.** Fix a native register \(J_R\) and the locks on the
two other roles \(U,V\), with product \(D=D_A D_B D_C\).
For the selected chart let \(H_U,H_V\) be the two complete
value diameters, set \(v_T=\log_n(1+H_T/(D_TD_R))\)
for \(T=U,V\), and put \(\theta=\log_n(1+X/D)\).
Call a pair good when all three primitive role quotients
satisfy the tall condition of JD2. Let \(F\) be any JD1
exponent for a fixed exact displacement on the good pairs;
in particular one can use the three estimates \(d_T=r_T/2\)
of JD2 with the actual outside caps and original degree bounds.
The good pair count per fixed nonzero difference
\(\eta=J_R-J'_R\bmod D_R\) is at most
\[
                            n^{v_U+v_V+F+o(1)} .
\]
Indeed each opposite exact shift is fixed modulo \(D_TD_R\).
Its own-lock residue vanishes and its foreign-lock residue
is the prescribed signed register difference.
There are \(O(n^{v_T})\) integer lifts. The two shifts
determine the third, so JD1 bounds the number of actual
pairs in each of those vector fibres. Pairs with any
vanishing exact shift are JD2 exceptions.

For a regular RB4 chart with \(n^s\) native cells of
\(n^{2u}\) actual rows each up to exponent tolerances, it follows that
\[
 E_+(J_R)\ \ll\
  n^{2s+o(1)}
  \left(1+n^{v_U+v_V+F-4u}+n^{s-\sigma}\right),
\]
where \(E_+\) is cyclic additive energy on the distinct native
residue labels and \(\sigma>0\) is the JD2 quotient saving.
For the proof, call a label pair good-rich if at least half
its original row-pair lifts are good. Each such pair contributes
at least \(n^{4u-o(1)}\) lifts to its fixed-difference
bound. The other distinct label pairs contribute at least
that scale to the total JD2 exceptional count
\(n^{2/3-\sigma+o(1)}\). As \(s=1/3-2u\), there are at most
\(n^{2s-\sigma+o(1)}\) such pairs. Their cyclic energy
with themselves is bounded by their total number times \(n^{s+o(1)}\).
Similarly sum the good-rich energy, using the pointwise
difference estimate after the lift rebate. The diagonal adds
\(n^{2s+o(1)}\), and Cauchy–Schwarz controls cross terms.
This uses the count of actual original row pairs and all three
original radical caps in JD1.

There is also a single-side estimate. Each occupied native
label has \(n^{u+o(1)}\) distinct opposite complete values.
For either \(T=U,V\), each native label pair has
\(n^{2u+o(1)}\) such complete value pairs. JD2's anchored
bound, divided by the actual original witness degrees, gives
at most \(n^{2\alpha-\sigma+o(1)}\) quotient-exceptional
complete value pairs. Divide by the lift count for a bad-rich
native pair and use \(s=\alpha-u\). For good-rich pairs the
fixed nonzero native difference has \(O(n^{v_T})\) complete
shift lifts, each with at most \(n^{r_T/2+o(1)}\) value
pairs. Thus one can replace the middle exponent in the
displayed energy bound by the minimum of its stated value
and \(r_T/2+v_T-2u\) for the two choices of \(T\).

The numerical lift counts can be replaced by actual catalogues
at each difference. Let \(r_R(\eta)\) count the native residue
pairs of difference \(\eta\), and let \(L_{UV}(\eta)\)
count the joint exact opposite shifts realized by good original
row pairs there. The good-rich term in the energy bound is
\[
       n^{F-4u+o(1)}
                 \sum_{\eta\ne0} r_R(\eta)L_{UV}(\eta).
\]
This follows by summing the pointwise bound for good-rich
pairs against their respective multiplicities. The same argument
with one opposite dictionary gives
\(n^{r_T/2-2u+o(1)}\sum_{\eta\ne0} r_R(\eta)L_T(\eta)\)
for \(L_T(\eta)\) counting shifts of tall-quotient complete
\(T\)-value witness pairs in the two native cells. The other
two role quotients need not be tall for this single-role list.
One can therefore use
a small weighted joint lift average even when the two opposite
numerical lift ranges are large. The remaining energy terms are
unchanged. If the native label is projected further modulo
\(m_R\mid D_R\), use the corresponding exact lift lists or
replace \(D_R\) by \(m_R\) in the numerical estimate. On any
regular retained projection of exponent \(t\) carrying the
full original row exponent \(d\) and complete dictionary
exponent \(\alpha\), the assumed regular projected fibre
sizes are \(n^{d-t+o(1)}\) original rows and
\(n^{\alpha-t+o(1)}\) opposite complete values. The same
proof replaces the two rebates by \(2(d-t)\) and
\(2(\alpha-t)\), respectively. The diagonal and exceptional
parts now have exponents \(2t\) and \(3t-\sigma\), at
\(d=1/3-o(1)\).

Another consequence of the fixed-difference estimate is the
following test on any selected good pairs. If their exponent
is \(p\) and the actual single register difference catalogue
costs \(q_R\), then
\[
                            p\le q_R+v_U+v_V+F+o(1).
\]
At occurrence, if the chart has exponent \(1/3-o(1)\)
and the pair population has exponent \(2/3-o(1)\), JD2
allows us to take only good pairs. With default caps and
\(v_U,v_V=o(1)\), the original degrees are subpower
so \(F\le1/6+o(1)\); thus \(q_R\ge1/2-o(1)\).
The proof retains the full three-role restrictions, including
an arbitrary outside radical profile in the quotient count.
The result is a difference growth condition. At positive
multiplicity, very deep relative two-lock charts also obey
the stricter JD5 entropy obstruction below.
[Regular native cells (RB4)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#the-actual-residue-cell-sizes)

## Tight full locks force almost complete register information

**JD5.** For each fixed allowance \(\eta>0\) there are
\(\theta,L_0>0\) with the following properties at
\(0<L<L_0\) fixed. Work with original primitive shell
rows and all three original individual caps
\(n^{1/3+o(1)}\), or any fixed larger bound on their
normalized exponents. Fix coprime locks \(D_R\mid R\) on
the indicated roles and the actual additive signs.
If two fixed locks satisfy \(D_R D_S\ge X^{1-\theta}\),
the number of rows over each fixed complete third letter
is \(O(n^\eta X^\varepsilon)\). With three fixed locks
and \(D_A D_B D_C\ge X^{1-\theta}\), fixing any one
complete letter and the three exact registers leaves at most
\(O(n^\eta X^\varepsilon)\) rows.

For the first assertion, fix the complete third letter.
It supplies the two cross residues, so one free locked
letter is fixed modulo the product \(D_R D_S\) by CRT.
For the second, fix one complete letter and the three
registers. Each free letter has its own zero residue and
the two opposite residues. Thus a free letter is fixed
modulo \(D_A D_B D_C\). In both cases it has the form
\(W=w_0+D j\) for a fixed \(w_0,D\) and a single integer
\(0\le j\le X/D\). Apply CS1 to the original fixed-letter
fibre, using the rank-one feature \(j\), the affine function
\(w_0+Dj=W\), and the numerical height bound \(W\le X\).
Taking \(\theta\) no larger than the CS1 threshold proves
the uniform bounds.
[Fixed-letter affine feature count (CS1)](ABC_CELL_PRIMITIVE_SECANT_BOUNDS.md#small-features-on-a-fixed-complete-letter-fibre)

In a regular balanced chart under RB5–RB7, write \(\alpha\)
for the complete dictionary exponent, \(u=1/3-\alpha\)
for the complete-letter fibre exponent, \(s=1/3-2u\)
for the register entropy, and \(z\) for the common pairwise
register mutual information, in units of \(\log n\).
All the following comparisons absorb the height-subpower
loss at each fixed budget and then the arbitrarily small
chart losses in the outer limit. The two-lock result
excludes every fixed \(u>0\) along a sequence of allowable
charts in which \((1-\log_X(D_R D_S))_+\) tends to zero.

For three locks, put \(\mathbf J=(J_A,J_B,J_C)\). The
second bound gives \(H(\mathrm{row}\mid T,\mathbf J)
\le(\eta+o(1))\log n\). The two registers opposite \(T\)
are already determined by \(T\). Hence
\[
 (u-\eta-o(1))\log n
 \ \le\ H(J_T\mid T)
 \ \le\ H(J_T\mid J_{T'},J_{T''})
 \ =\ (s-z+o(1))\log n.
\]
It follows that \(z\le1/3-3u+\eta+o(1)\) and
\(H(\mathbf J)\ge(1/3-2\eta-o(1))\log n\). The first
comparison also uses regularity to identify
\(H(\mathrm{row}\mid T)=(u+o(1))\log n\).
For an allowable sequence of charts with
\((1-\log_X(D_A D_B D_C))_+\to0\), choose \(\eta\) first,
then the associated numerical and budget cutoffs. Letting
\(\eta\) decrease, RB7 now gives the necessary conditions
\[
         u\le1/9,\qquad \alpha\ge2/9,\qquad
         z=1/3-3u,\qquad H(\mathbf J)=(1/3+o(1))\log n .
\]
Also \(H(T\mid J_{T'},J_{T''})=o(\log n)\) for each
role, since that register pair is a function of \(T\) and
its entropy \(2s-z+o(1)=\alpha+o(1)\) reaches the full
complete-letter entropy. Thus almost full numerical depth
excludes the whole limiting range \(1/6\le\alpha<2/9\)
and forces equality in the register-information constraint
throughout the rest of the range.
[Balanced register identities (RB5–RB7)](ABC_BALANCED_DICTIONARY_RESIDUE_CELLS.md#conditional-entropy-identities-and-three-registers)

The same entropy conclusion holds if relative depth is available
on just one free role \(U\). In each fixed-register cell choose a
rank-one coordinate anchored at an occupied \(U\)-value. Its
absolute value is at most \(H_U/D\), where \(H_U\) is the chart
diameter of this role. Fix one complete letter \(T\ne U\) and
apply CS1 with this feature. The proof needs
\(H(\mathrm{row}\mid T,\mathbf J)=o(\log n)\) for this one \(T\)
only, because \(s,z,u\) are already common to the roles under
RB5–RB7. Thus it suffices that \(D/H_U\ge X^{-o(1)}\).
For two locks the same reasoning uses either of the two free
roles, and only its diameter relative to the two-lock product.
As before choose the loss allowance before the numerical threshold
and the sufficiently small fixed budget.

## A primitive CRT model with near-diagonal energy

**JD11.** There is a tall three-lock CRT model with every
surviving JD5 exponent \(\alpha\in[2/9,1/3]\), primitive
positive rows, the full limiting registration entropy,
and near-Sidon modular dictionaries. It drops the
original radical caps.

Take comparable fixed prime-power locks
\(D_A=2^k,D_B=3^\ell,D_C=5^j\) of order \(H\), and set
\(X=D_AD_BD_C\). Pick phase numerators \(a,b,c\) in
fixed-relative intervals for
\(x=a/D_A\) around \(0.05\), \(y=b/D_B\) around \(0.75\),
and \(z=c/D_C\) around \(0.35\). Define
\[
        A=X(y-z),\qquad B=X(z-x),\qquad C=X(y-x).
\]
These are positive integers with exact \(A+B=C\),
each prescribed role divisible by its own lock, and
\(X/2<C<X\). The signed native registers are unit
multiples of \(a,b,c\) modulo their respective locks.
By congruence restrictions at 2,3,5 one can require
that all original role locks are exact and that each
of the two noncorresponding complete values is a
unit at every one of these primes. For example at a
prime of \(D_A\), take \(a\) a unit and choose \(b,c\)
so the residual \(D_C b-D_B c\) is a unit.
Rules at the other primes are independent by CRT.

For a large fixed \(P_0\), also impose one fixed
residue class on each phase numerator at each prime
\(p\le P_0\) outside 2,3,5 making common divisibility
impossible (for example \(a,b\equiv0,c\equiv1\)).
Each ambient phase box is now an arithmetic progression
of cardinality \(\asymp H\) and fixed step. A positive
fraction, in fact uniformly most, of ambient phase
triples is primitive. For primes outside the prescribed
small ones and at most the smallest box cardinality,
the two congruences for \(p\mid A,B\) require two
independent phase classes on fixing the remaining
phase, so have probability at most \(4/p^2\).
The sum is arbitrarily small by taking \(P_0\) large.
For larger primes, first fix \(b,c\). The nonzero
\(A\) has \(O(\log X)\) prime divisors. For each
possible remaining common prime the suitable \(a\)
classes have probability \(O(1/H)\), since the phase
box has order \(H\). This contribution tends to zero.

Use \(G=n^{z_0+o(1)}\) groups, each of \(m=n^{u+o(1)}\)
labels per role, where \(z_0=1/3-3u\ge0\). Choose
each role's \(Gm=n^{s+o(1)}\) distinct numerators
uniformly from its ambient phase box and partition
them into the groups. Keep only group-aligned Cartesian
triples. There are \(n^{z_0+3u+o}=n^{1/3+o}\)
registered rows, of which a positive fraction can
be chosen primitive by the preceding average.
Each sampled native set can simultaneously be chosen
strictly Sidon modulo its own lock (up to unordered
symmetry). The expected number of genuine collisions
is \(O(n^{4s}/H)\), which tends to zero at fixed
sufficiently small \(L\), by fixing all but one
sampled label. A two-label doubled collision is
excluded by the short relative phase interval.

A complete \(T\)-label determines the two opposite
phase labels exactly by CRT, so there are at most
\(Gm^2=n^{\alpha+o}\) labels and at most \(m\) rows
per complete label. Pruning to fixed fractions of
the permitted degrees yields a regular primitive
chart of the desired exponent, and similarly each
native dictionary has exponent \(s\). Its opposite
joint dictionaries have exponent \(\alpha\) and
the full native tuple support exponent \(1/3\).
All JD5 entropy equalities therefore hold at the
exponent scale. The three complete role modular
energies are at most four times the corresponding
squared dictionary sizes: a modular equal sum
must preserve the two unordered pairs of
opposite Sidon phase labels independently.
This persists after all deletions.

Most sampled rows have large radicals on all three roles. For example \(A/D_A=D_C b-D_B c\) is an
injective function of the two independent phase
numerators in their ranges, giving \(\asymp H^2\)
possible values. Only \(n^{1/3+o(1)}\) integers of
polynomial height have radical at most \(n^{1/3}\),
by fixed-support enumeration. The probability of
the target cap in any role thus tends to zero,
by the analogous calculation in that role.
One can delete every low-radical edge and retain
the same exponent. This establishes that the original radical caps are necessary
for a contradiction using the stated CRT, entropy and energy
conditions. The model itself meets the positivity, primitivity
and shell conditions.

## JD14. Costs of extracting deep projected registers

Use a regular balanced original cohort of exponent \(d=1/3\),
complete dictionary exponent \(\alpha\) and witness exponent
\(u=d-\alpha\), with all three original radical caps. Errors
below are in the fixed-budget height limit followed by the
small-budget limit.

Suppose an exponent-full selected population is covered by actual
full-divisor records \(D_A,D_B,D_C\). Choose a definite certificate
per row. Ordinary size and finite projection refinements may be
made before counting. Use only the native registers on roles
\(R,S\); let \(T\) be the remaining role. If two of the projected
quotients satisfy the small-height requirement of JD12, then one
of those short quotient roles is \(R\) or \(S\). Its projection
modulus divides \(D_R D_S\), including for partial registers.
For fixed \(T\), the two native registers are already determined
by \(T\) and \(D_R,D_S\), with the prescribed signs. The own
record \(D_T\) has only height-subpower choices from the divisors
of \(T\). Consequently the width in the remaining fixed-\(T\)
fibre is bounded by the already required projected-cell width.
The short role is in one fixed class modulo \(D_R D_S\) by the
signed sum and CRT. Its integer quotient after subtracting an
anchor is bounded by the prescribed short feature width. CS1 then gives arbitrarily small normalized
row multiplicity on this fibre, uniformly in the fixed values
and divisors. This is also the relative-width form of JD5.

It follows that the pair of opposite divisor records must have
conditional entropy at least \(u-o(1)\) given the exact remaining
complete letter \(T\). Indeed these two records and \(T\) determine
the full original row up to arbitrarily small normalized
entropy, while \(H(\mathrm{row}\mid T)=u+o(1)\). In
particular their actual joint catalogue has exponent at least
\(u-o(1)\). Summing that catalogue and charging the two projected
registers by their trivial complete-letter cap \(\alpha\) can
never give a fixed strict saving. One needs an actual projected
support discount, a joint entropy discount between charts and
registers, or a different cover.

If just one native register on \(R\) is used to shorten two
quotients, they are the two opposite roles. Fixing one complete
opposite value and the two locks on the other opposite role and
on \(R\) gives the same fixed-fibre argument. Hence each of the
two corresponding joint lock catalogues needs exponent at least
\(u-o(1)\). If \(R\)'s record has vanishing normalized price,
each opposite-lock marginal catalogue already needs exponent at
least \(u-o(1)\). On a near-full RB4 reference chart the default
native price is \(s=2\alpha-d\), so merely summing the two
opposite marginal prices and the native cap has no strict gain:
\(2u+s=d\). This last conclusion concerns separately summed prices.

These comparisons apply to adaptive, actual divisors on the original three-cap rows.

### The cost of a cover by bounded-doubling sets

On a regular RB4 chart fix a native register \(J_R\) of size
\(n^{s+o(1)}\), with \(s=d-2u>0\). Its occupied native
cells contain at most \(n^{2u+o(1)}\) original rows.
Suppose an opposite role \(V\) satisfies
\[
       v_V=\log_n(1+H_V/(D_V D_R)) < d-r_V/2
\]
by a fixed margin, where \(r_V\) is the actual outside
radical cap used in JD10. The single-side JD10 estimate for the full native
register gives
\[
 E_+(J_R)\ll n^{3s-\gamma+o(1)},\qquad
 \gamma=\min\{s,\sigma,d-r_V/2-v_V\}>0.
\]
Here \(\sigma\) is the genuine three-cap exceptional
quotient gap of JD2. All estimates use exact integer
difference lifts from the appropriate cyclic modulus.

For a subset \(P\subset J_R\) with
\(|P+P|\le n^\beta |P|\), Cauchy–Schwarz gives
\(E_+(P)\ge |P|^3/n^\beta\). It follows that
\[
 |P|\le n^{s-(\gamma-\beta)/3+o(1)},\qquad
 \#\{\mathrm{rows\ above\ }P\}
       \le n^{d-(\gamma-\beta)/3+o(1)}.
\]
In particular bounded doubling, and even subpower doubling,
cannot cover an exponent-full portion of this deep native
chart with subpower many such pieces. Any such cover has
positive catalogue price at least \((\gamma-\beta)/3\).
The argument applies on original three-cap charts.

The corresponding full projected version on a near-tight
JD5 chart of depth \(D\ge X/n^{O(1)}\) can use JD9.
For the projection of a complete role, with \(\alpha\ge
2/9\), the additive energy exponent is \(3\alpha-\gamma\)
with \(\gamma\) any fixed margin smaller than
\(\min(\sigma,1/18)\) at vanishing lift cost. The same
argument bounds the covered row exponent by
\(d-(\gamma-\beta)/3+v+o(1)\), where \(v\) is the
actual complete lift exponent for that projection.

A sparse subset of an affine integer box of
parameter volume \(n^B\) can have large doubling.
The energy argument in that generality only gives
\(|P|^4/n^B\le E_+(P)\), hence for an exponent-full
subset \(B\ge \alpha+\gamma\) in the full-projection
example. JD12 allows any fixed \(B\) at sufficiently small
budget, so sparse fixed-rank representations with larger
boxes remain possible.

### Escape-only locking across the balanced range

Let \(e_A,e_B,e_C\) be the escaped first moments of
the three original scaled-degree measures. For every
fixed \(c>0\), an original balanced endpoint with
\(\alpha\le 2/9-c\) satisfies
\[
                         e_A+e_B+e_C\le 1-\theta(c)
\]
for some fixed \(\theta(c)>0\). Otherwise choose a
height-first contrary sequence with total escape
approaching or exceeding one. Select the entire
high-order tail on each role, at a common large
bounded continuity cutoff. The total raw cost is
arbitrarily small by the bounded numerical first
moments, while the product has limiting numerical
weight at least the total escape. If a role has
zero escaped moment, augment its lock with an
arbitrarily cheap fixed positive-depth nonatomic
block. Such blocks exist on all roles since
\(u=1/3-\alpha>0\) by the fixed-fibre nonatomic
bound. Fixing the three divisors loses only the
arbitrarily small total cost.

On a large surviving chart, RB4 for the three
cheap original locks bounds the joint native
support by \(n^{3s+o(1)}\), where
\(3s=6\alpha-1\le 1/3-6c\). The height deficit
of the full product tends to zero. Choose the
JD12 counting allowance smaller than the margin
before taking its positive depth threshold and
then the smaller fixed budgets. The row count
after paying charts is below \(n^{1/3}\) by a
fixed margin, a contradiction.

For any positive limiting witness exponent \(u\),
FW2 gives all three strict inequalities
\[
                          e_R+e_S<h_{\min(R,S)}
\]
with uniform margin on a fixed positive \(u\) stratum.
Thus an escape-only selection cannot supply depth
relative to the complete heights on two roles after
omitting any native register. If one instead uses
still shorter projected-cell diameters, the earlier
conditional record floor excludes an exponent-full
selection with vanishing chart cost: two selected
native registers require the summed pair on at
least one of their own roles, and one selected
native register requires it on both opposite roles.
If the shortest complete height exponent is
\(h_{\rm short}\le1/2\), the two FW2 inequalities
involving it already give \(e_A+e_B+e_C<1\), so
even the three-register escape selection cannot
be fully deep in that case. At larger heights
and \(\alpha\ge2/9\), total escape at least one
may yield inexpensive three-lock charts. On those
deep charts the count of JD12 at any prescribed small
allowance forces the joint register entropy to reach
\(1/3+o(1)\) at arbitrarily small price for the locks.
A fixed entropy discount still needs arithmetic compression.
[Fixed-complete fibre features (CS1)](ABC_CELL_PRIMITIVE_SECANT_BOUNDS.md#small-features-on-a-fixed-complete-letter-fibre)
· [Fibre valuation and escape tests](ABC_WIDE_COFACTOR_CONGRUENCES.md)
