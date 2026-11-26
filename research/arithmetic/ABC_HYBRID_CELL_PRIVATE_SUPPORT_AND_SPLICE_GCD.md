# Private supports and common divisors of hybrid abc splices

A shallow residue and a short real interval can be imposed
after discounting the primes of an actual anchor. Within a
typical half-depth cell, this forces the opposite capped
supports to be almost private. A separate list argument
controls the common divisor of a splice and its nearby
genuine target. At sub-half depth the imposed modulus
accounts for almost all of this gcd. At any strict
subterminal depth it accounts for the powerful part up
to arbitrarily small normalized loss.
[Joint real and congruence moments (JRM1–JRM2)](ABC_JOINT_REAL_RESIDUE_MOMENTS.md)
[Positive splices of actual matching rows (SRS1–SRS2)](ABC_SHALLOW_REGISTER_SPLICES_AND_CARRIES.md)

## A fixed prime pool inside a hybrid cell

Put \(n=X^L\), and take the height limit at fixed \(L\).
Consider a chart of genuine positive primitive rows
\(A+B=C\), \(X/2<C\le X\), with the three original
radical caps retained. The normalized comparison exponent
is \(d=1/3\). On a role \(U\) assume raw cap \(d+o(1)\)
and a fine exact-layer box as in (SRL4), with an optional
fixed coefficient of polynomial height. The maximum
normalized cap of a single variable exact-layer root
tends to zero. The notation \(o(1)\) allows the stipulated
outer box and chart errors and any prescribed inner
counting allowance. All fixed positive margins are
chosen before taking sufficiently small fixed budgets.

Let a prescribed integer modulus \(M\) of depth
\(\sigma=\log_n M\) be a unit modulus for the values
being counted. A hybrid cell fixes a unit residue
at \(M\) and a log interval of length \(O(n^{-\xi})\)
for \(U\), where \(s=\sigma+\xi\) is the combined
depth. The interval and the residue may depend on
fixed conditioning data. For a pure residue cell
one can take \(\xi=0\) and cover the full log range
by \(O(\log X)\) intervals. Let \(u\) be an upper
exponent for the complete degree of \(U\) in the
original chart. In a complete-label matching we
take \(u=0\).

**HPS1.** For any prescribed prime pool \(P\) supported
on a polynomial-height integer, the number of rows in
one such cell on which the \(P\)-supported radical of
\(U\) is at least \(n^f\) is at most
\[
    O_{L,\varepsilon}\!\left(
      n^{\,u+\max\{(d-f)/2,\ d-f-s\}+\eta}
                     X^\varepsilon\right) .                  \tag{HPS1}
\]
Here \(f>0\) and \(\eta>0\) are fixed, the chart and
fine-layer errors are sufficiently small in terms of
the desired counting allowance, and the bound is
uniform in the pool and prescribed cell. If \(d-f\)
is negative by a fixed margin the event is empty.

Fix the exact \(P\)-supported factor at height-subpower
cost and divide it out. This is legitimate for a
pool selected from an anchor row, since that row is
already fixed. The factor is a unit at \(M\) and its
logarithm translates the window. The remaining exact
layers still have the original orders. One can record
the supported pieces in each actual layer separately,
with the same height-subpower cost, if the fixed
coefficient shares some of their primes. The variable
squarefree radical factors are disjoint, and after tight dyadic
binning their total cap is at most \(d-f+o(1)\).
Their individual caps are bounded by the original fine
caps, up to the vanishing bin allowance. The fixed
coefficient is handled by the same division. At each
fixed budget the additional dyadic refinement costs
height-subpower. The original row and its two other
letters are left unchanged.
[Uniform enumeration on fixed supports](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)

Apply (JRM1) with one residue-window pair and
\(\vartheta=1\). For clarity, with residual total cap
\(r'\) and maximum cap \(a_*\), its gain satisfies
\[
                     G_2(s)\ge \min(r',2(s-a_*)_+).
\]
For positive \(t=\min(r',2(s-a_*)_+)\), put
\(h_i=\lceil t/(2a_i)\rceil\) at positive-cap positions.
Then \(h_i a_i\le s\), so the slot gain at \(1/h_i\)
is \(a_i\). If the reciprocal weights have sum at most
two, use all of them and gain \(r'\). Otherwise scale
them down to sum two. Concavity of the slot gains and
\(a_i\ge t/(2h_i)\) give gain at least \(t\).
Thus the free complete-value exponent is at most
\(\max(r'/2,r'-s)+o(1)\). Restore the true degree
bound \(n^{u+o(1)}\) for each complete value to obtain
(HPS1). This proof also gives the same bound for a
direct outside-pool radical cap \(n^{d-f}\) after
fixing the supported part.

## Almost private factors in an occupied cell

Take a complete-label matching and a hybrid cell on
\(U\) containing \(n^{d-s-o(1)}\) actual rows, with
\(0\le s<d\) fixed. Such row exponents hold on typical
occupied cells in the fine-layer zero-multiplicity
setting of (JRM2), when the total mixed alphabet is
at most \(n^{s+o(1)}\). Use \(P\) equal to the support
of one anchor \(U_0\), with any additional prescribed
polynomial-height pool adjoined if needed. For
ordered pairs of rows in the cell, those with raw
radical overlap on \(U\) at least \(f\) are power-smaller
than the total, provided
\[
                 \min\{f,(d+f)/2-s\}>0.
\]
Indeed choose the first row, then apply (HPS1) to
the second. Comparison with its cell exponent
\(d-s\) gives the displayed loss before allowances.
At \(s=d/2\) this gives a loss \(f/2\) for each fixed
\(f>0\). At depths \(s\le(d-f)/2\) the loss is \(f\).

In particular a pure shallow register of depth \(d/2\)
on a marked role gives this privacy simultaneously on
both opposite roles in a typical \(n^{d/2-o(1)}\) cell,
under the fine-layer hypothesis on both. The sign in a
shared register changes only the prescribed unit residue
on that role. Adding a real interval on just one letter
guarantees the combined-depth assertion on that letter.
A simultaneous assertion on a different role uses its
own available congruence and real conditions. Privacy
inside these small cells is stronger than applying the
global prescribed-pool discount to the second row and
then conditioning on the entire cell.

Here is an application to exact target dictionaries
from a half-depth shared register. Work on a complete
matching with the primes coloured by original roles,
as in (IM1). The opposite roles of the marked letter
are \(U,V\). Fix \(k\ge2\), prescribed positive integer
orders \(a_i,b_i\), and prescribed nonzero rational
weights \(w_i,z_i\) of polynomial height. Use a pool
\(P\) containing all the weight supports. On a typical
\(k\)-tuple in a full-exponent half-depth cell consider
\[
                R=\sum_{i=1}^k
                     \bigl(w_i U_i^{-a_i}+z_i V_i^{-b_i}\bigr).
\]
For any fixed \(\epsilon_0>0\), outside a power-small
proportion of tuples the participating radicals have
outside mass at least \(d-\epsilon_0\) on each opposite
role, and every pairwise shared radical on either role
has mass less than \(\epsilon_0\). These conclusions
follow from (HPS1) for a single value and for an anchored
pair. A full outside cap below \(d-\epsilon_0\) is
handled by the direct outside-cap version of (HPS1).
A union bound suffices for the fixed \(k\).
[Global role colouring (IM1)](ABC_OCCURRENCE_INDUCED_ADDITIVE_MATCHINGS.md#extracting-an-exact-induced-matching-under-all-three-caps)

Every private prime outside \(P\) is then a unique pole
of its summand and survives in the reduced denominator
of \(R\). This denominator has outside radical at least
\(n^{2kd-O_k(\epsilon_0)}\). The same argument gives a
list bound from the denominator alone. Let \(J\) be the
product of the primes outside \(P\) appearing in two or
more same-role slots. Its mass is \(O_k(\epsilon_0)\),
so it can be recorded at that cost. Outside \(P\) and
\(J\) the denominator's primes can be assigned to their
slots at height-subpower cost, with the role known by
the colouring. The valuation on an assigned private
prime is the original valuation times the corresponding
order \(a_i\) or \(b_i\), so it specifies the original
prime power. Enumerate the full parts on \(P\cup J\)
by the fixed-support bound. The resulting complete
\(U_i\)'s specify the original rows in the matching.
Thus at least \(n^{kd/2-O_k(\epsilon_0)-o(1)}\) distinct
exact denominator labels occur in that cell. This applies
even if the marked complete value varies from row to
row. Small reciprocal images on the half-depth cells
would require a different arithmetic restriction.
[Private poles in a fixed complete fibre](ABC_BALANCED_FIXED_SHIFT_REVERSE_TESTS.md#generic-reciprocal-sums-have-large-outside-denominators)

## Gcds of a splice with its compared true target

Use the complete matching and the near-splice setting
of (SRS2). Write \(M=M_W\), where \(W\) is the marked
computed role. Let \(\sigma=\log_n M\), let the real
difference scale be \(Xn^{-\xi}\), and put \(s=\sigma+\xi\)
with \(0\le s<d\). The compared original role is tall
and has the fine-layer hypothesis. In a sum splice
one can use \(W^*=A_i+B_j\) and \(W_j=C_j\), comparing
\(A_i\) with \(A_j\). In a difference splice take
\(W^*=C_i-V_j>0\) and \(W_j=C_j-V_j\), comparing the
tall targets \(C_i,C_j\). Denote the two compared values
by \(Z_i,Z_j\). They are \(\asymp X\), by the selection
in SRS2. Put \(G=\gcd(W^*,W_j)\) and \(g=G/M\). Then
\[
               M\mid G\mid W_j,\qquad G\mid Z_i-Z_j .
\]
Every prime on \(G\) has the marked colour of the original
\(W_j\), so the compared complete values remain units at
\(G\). The full near-splice population has exponent
\(2d-s\) by (SRS2).

**HPS2.** If \(s<d/2\), then for every fixed
\(\zeta>0\) all but a power-small fraction of the near
splices have \(g<n^\zeta\). The loss on the exceptional
population is at least
\(\min\{\zeta,d/2-s\}\) before the chart and counting
allowances.

Fix the second original row and enumerate the divisors
\(G\mid W_j\) at height-subpower cost. For \(g\ge n^\zeta\)
the first compared value lies in a single hybrid cell
with combined depth at least \(s+\zeta\). Its real
interval has logarithmic width \(O(n^{-\xi})\) since
\(Z_j\asymp X\). The \(K=1,\vartheta=1\) list bound
on the original compared role gives at most
\(n^{\max(d/2,d-s-\zeta)+o(1)}\) choices. Sum over the
\(n^{d+o(1)}\) second rows. This proves the assertion.
At \(s=d/2\) this fixed-anchor list bound alone has
the same exponent as the reference population.
[Near splices at two scales (SRS2)](ABC_SHALLOW_REGISTER_SPLICES_AND_CARRIES.md#near-splices-in-both-residue-and-numerical-scales)

**HPS3.** At every fixed \(s<d\), let \(P\) be a fixed
polynomial-height prime pool. If an additional divisor
\(E\mid g\) has depth at least \(\zeta\) and outside
radical mass at most \(b\), the number of such near
splices is at most
\[
               n^{2d-\min(d,s+\zeta)+b+o(1)} .
                                                               \tag{HPS3}
\]
This is a power loss if \(b<\min(d-s,\zeta)\) by a
fixed margin. Consequently, if the marked original
role also has no positive-size individual prime in
its box, then for every fixed \(\eta_0>0\), almost all
near splices by a power margin satisfy
\[
                   \frac{g}{\operatorname{rad}_{P^c}(g)}
                             < n^{\eta_0}.
\]
In particular the powerful part of the extra gcd
cannot carry a fixed positive normalized depth, nor
can any part supported on the prescribed fixed pool.

To prove the first assertion, enumerate \(E\) at
cost \(n^{b+o(1)}\) by its outside radical and uniform
fixed-support enumeration including \(P\). The range
of \(E\) has polynomial height since \(E\mid W_j\).
For a fixed \(E\) from this list, both compared
values are units modulo \(ME\), and their residues
agree there. This assertion uses \(E\mid G/M\), so
it holds even if \(E\) and \(M\) are not coprime.
Apply the joint collision bound (JRM2) to the whole
comparative-role dictionary at modulus \(ME\) and
the same real depth. Real differences up to
\(Xn^{-\xi}\) occupy only boundedly many neighboring
log bins for tall values. Cauchy–Schwarz between
those occupancies costs a bounded factor. Fine
layers give the full mixed gain
\(\min(d,s+\zeta)-o(1)\). Summing over \(E\) proves
(HPS3).

For the last assertion let \(a=d-s>0\).
First consider the case where the outside radical of
all primes of \(g\) with exponent at least two has
mass at most \(a/2\). Take \(E\) to be the full part
of \(g\) supported on those primes together with the
full \(P\)-supported part of \(g\). If the displayed
ratio is at least \(n^{\eta_0}\), \(E\) has depth
at least its outside radical mass plus \(\eta_0\).
Bin that mass in a sufficiently fine fixed mesh.
Then (HPS3) saves by at least
\(\min(a/2,\eta_0)\) before mesh and counting errors.
In the remaining case choose just squares of repeated
outside primes with total radical mass between \(a/4\)
and \(a/3\). This is possible by greedy selection
at late boxes, since the full repeated support has
mass above \(a/2\) and the maximum individual prime
cost tends to zero. The resulting \(E\) has depth
twice its radical mass and again gives a fixed gap
in (HPS3). This proves the assertion.

For a squarefree \(E\) outside the prescribed pool, its
generic record cost \(b\) equals its depth \(\zeta\).
The preceding list estimate then gives no saving from
that record. No argument here forces the large fresh
radical of a splice to fall in the cap of the original
computed letter, or bounds the number of occupied
hybrid cells. The estimates quantify the local support
requirements of an endpoint profile.
[Joint real and congruence collision bound (JRM2)](ABC_JOINT_REAL_RESIDUE_MOMENTS.md#joint-collision-entropy-and-conditional-full-factors)


## A joint radical floor for a splice and its shift

**HPS4.** Work in the near-splice chart of SRS2, keeping
all three original radical caps and a marked computed
role \(W\) of height at least \(X^{h_0}\) for fixed
\(h_0>0\). In a complete-label matching let \(W^*\)
be the positive splice output compared with \(W_j\)
on the original anchor row. Write \(s=\sigma+\xi\)
for the combined depth and assume \(s<d/2\) as in
HPS2. The common fixed pool \(P\) contains the
prescribed full locks and has bounded normalized
radical mass \(\log_n\operatorname{rad}(P)=O(1)\),
as for the fixed union of the original locks in SRS2. For every fixed \(\delta>0\),
outside a power-small fraction of the
\(n^{2d-s+o(1)}\) near splices one has
\[
  \log_n\operatorname{rad}_{P^c}(W^*)
     +\log_n\operatorname{rad}_{P^c}(|W^*-W_j|)
   > 3(2d-s)-r_W+\mu_W-\delta .
\]
The stated threshold is used when positive. At
balanced caps and \(d=1/3\) it tends to \(7/6+\mu_W\)
as \(s\) tends to \(d/2\) from below and the
fixed gap decreases.

In fact HPS2 lets us discard a power-small
fraction and take \(G=\gcd(W^*,W_j)=M g\), with
\(g<n^\zeta\) for an arbitrarily small fixed
\(\zeta>0\). Here \(M\) is the imposed shallow
modulus of depth \(\sigma\). The shift
\(W^*-W_j\) is nonzero because the compared
source values come from distinct rows of a
complete matching. Divide the pair and its shift
by \(G\). These three nonzero integers give a
positive primitive additive triple after choosing
its sign orientation. Its largest term is at
least \(W_j/G\ge X^{h_0/2}\) for sufficiently
small fixed \(L\). If the left side of the desired
bound is at most \(b\), the effective radical
budget of this derived triple, with \(P\) free,
is at most \(r_W-\mu_W+b\), up to the stipulated
box errors. Its raw budget is also bounded in
base \(n\). The effective one-third bound
therefore gives at most
\(n^{(r_W-\mu_W+b)/3+\eta}X^\varepsilon\)
such primitive normalized triples. Use a fixed
mesh in the physical height exponent between
\(h_0/2\) and a fixed polynomial-height upper
bound to obtain the indicated arbitrary small
\(\eta\) after conversion to base \(n\).

Supply \(g\) and the finitely many sign and
role choices at cost \(n^\zeta\) times a constant.
Restoring \(G\) gives \(W_j\) and \(W^*\).
The complete \(W_j\) identifies its original
anchor in the matching. Then the known signed
splice formula recovers the remaining original
source value from \(W^*\) and the fixed
anchor operand. That complete value identifies
the other row. Thus the direction of the
inequality survives lifting back to the actual
pair count. Take \(\zeta\) and the counting
allowance smaller than the desired margin
\(\delta/3\) and compare with the full population.
[Near-splice population at mixed depth (SRS2)](ABC_SHALLOW_REGISTER_SPLICES_AND_CARRIES.md#near-splices-in-both-residue-and-numerical-scales)
· [Effective count on prescribed supports (EB4–EB5)](ABC_EFFECTIVE_RADICAL_BUDGET.md#fixed-supports-and-polynomial-prime-prices)


## Physical shifts with small outside radical are rare in a cell

**HPS5.** Take a complete matching in the same fine-layer
chart and compare a fixed role \(U\) in a hybrid cell of
depth \(s<d/2\) and size \(n^{d-s-o(1)}\). Suppose
each row has a fixed positive radical mass \(\gamma\)
on a compact positive scaled-degree band \([a_0,a_1]\)
of \(U\). This is available in the regular endpoint
boxes with positive mass on the indicated role.
For any fixed \(B,K>0\) and any pool \(P_x\) of
prime support of polynomial height \(X^K\) prescribed
by the anchor \(x\), the pairs with
\[
                \operatorname{rad}_{P_x^c}
                          (|U(y)-U(x)|)\le n^B
\]
are power-smaller than the cell pair exponent.
The same conclusion holds across the typical
hybrid cells producing the SRS2 near splices,
so that below half depth the physical
difference itself generically has outside
radical beyond every fixed normalized cap.
The budget threshold can depend on that cap and
on the fixed cell and profile margins.

Indeed take any fixed \(0<f<\gamma/2\). Apply HPS1
with the prime support of \(U(x)\) to suppress
the rows in the cell whose \(U\) overlaps the
anchor's support at mass at least \(f\). Both
\((d-f)/2\) and \(d-f-s\) are strictly below
the target neighbour exponent \(d-s\).
In the remainder put \(g_0=\gcd(U(x),U(y))\).
The band outside the anchor's support contributes
numerical height at least \(a_0(\gamma-f)\)
in base \(X\) to \(U(y)/g_0\). Fix \(g_0\mid U(x)\)
at divisor cost. In the primitive normalized
difference triple the prescribed complete endpoint
\(U(x)/g_0\) has the original radical cap, the
other complete quotient has the same cap, and
the outside-\(P_x\) radical of the difference
quotient is at most \(n^B\). Apply FL1 as in
FE1, dividing into height shells between the
positive numerical floor and \(X\).
This gives per-anchor exponent at most
\(d/2+\eta\) after the subpower divisor choices
and the arbitrarily small allowances. The complete
values in the matching determine their rows. Taking
\(\eta<d/2-s\) by a fixed margin proves the claim.
Adjoining earlier polynomial-height shift supports
to the pool is allowed. The use of HPS1 localizes
the overlap exception in FE1 to the hybrid cell scale.
[Fixed-endpoint physical difference bound (FE1)](ABC_FIXED_ENDPOINT_SHIFT_RADICAL_COUNTS.md#a-physical-difference-on-one-role)


## An interpolation bound for the mean squarefree excess

**HPS6.** Consider a pure shallow residue cell of depth
\(\sigma<d\) on a marked target \(C\) in a complete-label
matching, with a common tall summand \(Z=A\) on its rows.
Assume the fine layers on the compared role and the small
individual outside primes on the marked target required by
HPS3's gcd-excess conclusion. Take an exponent-full cell
with \(N\ge n^{d-\sigma-\tau}\) rows. Use cells typical
for the HPS3 pair estimate, so the conditional exception
there is power-small. Indeed if the global bad pair
count is \(n^{2d-\sigma-\delta+o(1)}\), the cells of
this size with bad-pair fraction at least \(n^{-\delta/2}\)
have total row mass at most \(n^{d-\delta/2+\tau+o(1)}\).
Choose \(\tau<\delta/4\). If \(M\) is
the imposed shallow modulus and
\[
       g'_{ij}=\frac{\gcd(C_j,A_i+B_j)}{M}
                =\frac{\gcd(C_j,Z_i-Z_j)}{M}
                        \quad(i\ne j),\qquad g'_{ii}=1 ,
\]
then \(\mathbb E_{i\ne j}\log_n g'_{ij}\le\sigma+o(1)\)
in the height-first limit with arbitrarily small chart
allowances. At depth zero this gives vanishing mean
normalized log gcd for the splice and the true target.
If the fixed HPS3 pool contains a common full target lock
of raw mass \(\mu\), the mean bound is
\(\max(0,\sigma-\mu)+o(1)\).
The mean here concerns splice gcds, which differ from
the common gcds of two original coordinate differences.

Fix \(m\) independently sampled anchors in this cell.
Take a nonzero integer polynomial \(F(Z,W)\) of degree
at most \(r\) vanishing at every anchor \((Z_j,C_j)\),
where \(m=\lfloor\binom{r+2}{2}/4\rfloor\). The
coefficient pigeonhole used in JSI1 gives values of
absolute size at most \(X^{3r}\) in the height box.
Divide \(F\) by its largest power of \(W\); it still
vanishes at the anchors, since \(C_j>0\). Its
specialization \(F(Z,0)\) is now a nonzero polynomial
of degree at most \(r\). For each query \(i\) and
distinct anchor \(j\),
\(\gcd(C_j,Z_i-Z_j)\mid F(Z_i,0)\). This follows by
reducing the two arguments in the anchor evaluation.
As complete \(Z_i\)'s are distinct, at most \(r\) queries
have \(F(Z_i,0)=0\). If \(D_m\) is the lcm of the
\(g'_{ij}\) for these anchors, then
\[
       \mathbb E_i\log_n D_m
           \le \frac{3r}{L}+\frac{rm}{LN}.
\]
This interpolation uses the nonzero specialization to a
coordinate axis rather than a curve count for the input
chart, which can be much smaller than the family in JSI1.

To price overlap, let \(P_{\rm prev}\) be the pool of
primes of the earlier marked values. It has polynomial
height with exponent \(m\) fixed at the budget.
For the next independent \(j\), the count with
\(\log_n\operatorname{rad}_{P_{\rm prev}}(C_j)\ge
\sigma+2\eta\) is at most
\(n^{d-\sigma-2\eta+o(1)}\) by outside radical
enumeration and the complete matching. For
\(\tau\ll\eta\) this is a power-small fraction of
the cell. All earlier common factors use primes
in \(P_{\rm prev}\). Thus
\(\log_n\gcd(g'_{ij},D_{\rm prev})\) is bounded by
the raw overlap of \(C_j\) with \(P_{\rm prev}\)
plus the powerful and fixed-pool excess of \(g'_{ij}\)
from HPS3. The two expected contributions are at
most \(\sigma+2\eta+o(1)\) and \(\eta+o(1)\),
respectively, using the HPS3 test with threshold
\(\eta\). The fixed pool in that test need not
contain \(P_{\rm prev}\). This follows by writing
\(g'_{ij}\) as its squarefree outside-pool part times
the excess, and using the upper bound for the first
factor on the earlier support. If the fixed pool contains
the common full target lock, that lock is also in
\(P_{\rm prev}\) after the first anchor. Its radical
mass \(\mu\) can be subtracted when bounding the
squarefree outside-pool part of the overlap. The
first overlap is zero, so this sharpens the mean
overlap to \(\max(0,\sigma-\mu)+3\eta+o(1)\).

Use \(\sum\log g'_j=\log D_m+
\sum\log\gcd(g'_j,D_{\rm prev})\) and average.
Choose the fixed \(r\) at the budget large enough that
\(3r/(Lm)<\eta\). With \(N\) of positive exponent,
\(rm/(LN)\) is negligible in the height limit.
The diagonal convention has negligible effect.
This gives \(\sigma+O(\eta)\) and hence the stated
iterated bound. A bound \(o(1)\) on the conditional
marked-role raw overlap, in this same cell, would
remove the term \(\sigma\). The established half-depth
privacy concerns the opposite roles.
[Interpolation through fixed-height anchors (JSI1)](ABC_SHALLOW_JOINT_INTERPOLATION.md#interpolation-against-proper-congruence-hits)
