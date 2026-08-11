# F2 generic `k=1` affine-row logarithmic Chern packet

> **Status.**  This note computes the divisorial logarithmic `ch_2` term of
> any logarithmically immersed normalized target curve, then specializes it
> to the generic nodal `k=1` affine-purity row.  If its reduced boundary
> component `E` maps with residue degree `f` to the normalized target
> quintic, has transverse index `e`, and has self-intersection `E^2=-n`, its
> cyclic divisorial contribution on the current target model is
> `e*f*(b-7)-e^2*n/2`.  Here `0<=b<=8` records how many smooth centers of the
> special carrier extraction the target quintic follows after the mandatory
> `(5,2)` extraction.  Away from `lambda=125/729` one has `b=0`; at the
> special value `1<=b<=8`.  The resolved-terminal nonattachment theorem
> shows that the terminal passport does not determine `e`.  Conditional on
> a root-plus-one-affine-packet filtration, the squarefree/double degree-floor
> residuals are
> `(e^2*n-2*e*f*(b-7)-20-s_X)/2` and
> `(e^2*n-2*e*f*(b-7)+17-s_X)/2`.  Their parity and positivity give exact
> all-signature compiler tests.  On the ordinary-cusp face, the local
> attachment dichotomy shows that the finite cusp ledger is not universally
> `2f`: for `h` points in the residue-degree-`f` fiber, `c` of them actual
> boundary nodes, the minimal isolated ledger is `2f-h+c`, ranging from `f`
> for unramified smooth-boundary folds to `2f` for node-saturated SNC
> attachments.  These are conditional K-theory sieves, not construction of
> the pullback, proof that `k=1` occurs, or an exclusion of `(75,125)`.

The local module, conormal degrees, target intersection changes, and the two
degree-floor sieves are replayed by
[`verify_f2_affine_k1_log_ch2.py`](../scripts/verify_f2_affine_k1_log_ch2.py).
The reusable exact arithmetic is in
[`log_node_profiles.py`](../jcsearch/log_node_profiles.py).

<!-- status-consumer: LCCT1 2fc6ecea7a7c8b49 -->

<!-- status-consumer: LCHB1 176bf85520516fa6 -->

<!-- status-consumer: PF2K1C1 358a6ba820e8b2f1 -->

<!-- status-consumer: PF2PPA1 b24c4d80c2f8230e -->

## 1. Setup

Let

\[
 \nu:\mathbb P^1\longrightarrow \bar C\subset Y                 \tag{1.1}
\]

be the normalization of a generic `k=1` target from
[`F2_AFFINE_TARGET_K1_COLLISION.md`](F2_AFFINE_TARGET_K1_COLLISION.md).
Its four affine singularities are ordinary nodes, so `nu` is unramified and
is a branchwise immersion at both preimages of every node.  It has one point
over the target boundary.  Let a reduced source boundary component `E` map
to this normalization with residue degree `f`, and let its transverse
ramification index over `C` be `e>=2`.

The exact implicit equation and the eight-point normalization conductor are
given in
[`F2_AFFINE_TARGET_K1_IMPLICIT_CONDUCTOR.md`](F2_AFFINE_TARGET_K1_IMPLICIT_CONDUCTOR.md).
They identify the finite target singularity correction that must accompany
the generic divisorial formula below.

<!-- status-consumer: PF2K1I1 a7582c1e36140840 -->

Write `D_Y` for the target boundary already containing both the `(5,2)` and
special `(5,36)` extraction clusters.  The normalization map is logarithmic
at infinity, and the logarithmic conormal sequence is

\[
 0\longrightarrow K_\nu
 \longrightarrow \nu^*\Omega_Y^1(\log D_Y)
 \longrightarrow \Omega_{\mathbb P^1}^1(\log\infty)
 \longrightarrow0.                                      \tag{1.2}
\]

The ordinary conormal sequence underlying (1.2) is the standard exact
sequence for an immersion; see the Stacks Project,
[*Exact sequences of differentials and conormal sheaves*](https://stacks.math.columbia.edu/tag/06BB).
The logarithmic form follows directly in the SNC coordinates at the unique
point at infinity.

## 2. Generic local Smith profile

At a smooth point of the target normalization, choose an affine transverse
equation `x` and a normalization coordinate `y`.  At the generic point of
`E`, source coordinates may be chosen so that

\[
 x=u^e\cdot\text{unit},\qquad y=\phi(t)+u(\cdots).        \tag{2.1}
\]

In the logarithmic source basis `(du/u,dt)`, after invertible row and column
operations the differential has generic matrix

\[
 \begin{pmatrix}e u^e&0\\0&1\end{pmatrix}.                \tag{2.2}
\]

Consequently

\[
 \operatorname{Fitt}_0=(u^e),\qquad
 \operatorname{Fitt}_1=R,qquad
 \operatorname{coker}\simeq R/(u^e).                     \tag{2.3}
\]

Thus the generic divisorial packet is cyclic and its determinant divisor is

\[
 D=eE.                                                     \tag{2.4}
\]

Ramification of `E->P^1`, a failure of the target normalization to immerse,
or a collision with another boundary packet may make `Fitt_1` nonunit at
finitely many points.  Those corrections must be retained as separate
point-supported quotients.  They do not change the generic calculation
(2.2)--(2.4).

## 3. The conormal kernel degree

Taking degrees in (1.2), using
`deg Omega_P1(log infinity)=-1`, gives

\[
 \deg K_\nu=L_Y\mathbin{\cdot}\bar C+1.                  \tag{3.1}
\]

Pullback along the residue-degree-`f` map from `E` gives

\[
 \deg_E K_{\mathrm{red}}
 =f\left(L_Y\mathbin{\cdot}\bar C+1\right).              \tag{3.2}
\]

On the Cartier thickening `D=eE`, the filtration by powers of the ideal of
`E` has `e` graded pieces.  Tensoring all pieces by the kernel line changes
each Euler characteristic by `deg_E K_red`, so

\[
 \deg_D K=e\,\deg_EK_{\mathrm{red}}
 =ef\left(L_Y\mathbin{\cdot}\bar C+1\right).              \tag{3.3}
\]

The cyclic cokernel formula from
[`LOG_CYCLIC_COKERNEL_TWIST.md`](LOG_CYCLIC_COKERNEL_TWIST.md) now gives

\[
 \boxed{
 \operatorname{ch}_2(T_E^{\log})
 =ef\left(L_Y\mathbin{\cdot}\bar C+1\right)
  +\frac{e^2E^2}{2}.}                                    \tag{3.4}
\]

This removes the kernel Gauss degree from the list of unknowns for every
immersed target normalization.  Only the target log intersection, the
finite-cover data `(e,f)`, the source self-intersection, and point failures
of cyclicity remain.

More generally, let `C_tilde` be a smooth target normalization of genus `g`
with a reduced logarithmic puncture divisor `S` of degree `s`.  Then

\[
 \deg\Omega_{\widetilde C}^1(\log S)=2g-2+s,
\]

and exactly the same conormal and Cartier-filtration argument gives the
degree-independent immersed-curve packet

\[
 \boxed{
 \operatorname{ch}_2(T_E^{\log})
 =ef\left(L_Y\mathbin{\cdot}\bar C-(2g-2+s)\right)
  +\frac{e^2E^2}{2}.}                                  \tag{3.5}
\]

Formula (3.4) is the rational one-puncture specialization `(g,s)=(0,1)`.
As before, ramification of the residue cover and failures of immersion or
cyclicity contribute separate point-supported terms in the exact global
filtration; they are not silently included in (3.5).

## 4. The `k=1` target intersection

On `(P^2,L_infinity)`, the logarithmic canonical class is `-2H`.  Since the
target is a quintic,

\[
 (K_{P^2}+L_\infty)\mathbin{\cdot}\bar C=-10.             \tag{4.1}
\]

In target coordinates `(a,b)=(1/(-Q),P/(-Q))`, the infinity branch has
orders `(5,2)`.  The first blowup extracting `(5,2)` is a smooth-boundary
blowup at a point of curve multiplicity two, so it raises (4.1) by two.  The
other three blowups in that cluster are boundary-node blowups and are log
crepant.  Hence, before the special carrier centers,

\[
 L_Y\mathbin{\cdot}\bar C=-8.                            \tag{4.2}
\]

The carrier target cluster contains eight smooth-boundary blowups followed
by four node blowups.  Let `b` be the number of its smooth centers traversed
by the strict transform of `C`.  If the curve meets a different point of the
`(5,2)` divisor then `b=0`; otherwise `1<=b<=8`.  The strict transform is
smooth at all these infinitely near points, so every traversed smooth center
raises the intersection by one:

\[
 \boxed{L_Y\mathbin{\cdot}\bar C=-8+b,\qquad0\le b\le8.} \tag{4.3}
\]

The integer `b` has an exact fixed-target jet description.  After the seven
certified carrier shears, write

\[
 w=h-\frac{125}{729}-c_1\pi-c_2\pi^2-\cdots-c_7\pi^7.   \tag{4.4}
\]

For `lambda=125/729`, put `r=ord_u(w|_C)`, allowing `r=infinity` if the
restriction vanishes identically.  Since `ord_u(pi)=1`, the curve valuation
is `(1,r)`.  The eight smooth carrier blowups insert the rays

\[
 (1,1),(1,2),\ldots,(1,8),
\]

and a direct chart calculation gives

\[
 \boxed{b=\min(r,8).}                                   \tag{4.5}
\]

Thus `b` is computable from only the first eight fixed-target infinity jets.
The first test can be read directly from the parametrization without
expanding the implicit equation.  If

\[
 \begin{aligned}
 P&=A u^{-3}(1+\alpha u+O(u^2)),\\
 -Q&=C u^{-5}(1+\beta u+O(u^2)),
 \end{aligned}                                         \tag{4.6}
\]

then

\[
 \pi=\frac{A^3}{C^2}u+O(u^2),\qquad
 h=\frac{A^5}{C^3}
 \left(1+(5\alpha-3\beta)u+O(u^2)\right).              \tag{4.7}
\]

On the special-residue locus `A^5/C^3=125/729`, one therefore has `b=1`
unless

\[
 \boxed{
 \frac{125}{729}(5\alpha-3\beta)
 =c_1\frac{A^3}{C^2}.}                                 \tag{4.8}
\]

Successive values `b>=j` are obtained by the successive vanishing of the
coefficients of `u,...,u^(j-1)` in (4.4).  This turns the contact problem
into seven explicit jet equations in the lower Laurent data.

The exact fixed-coordinate compiler is now given in
[`F2_AFFINE_K1_CARRIER_JET_FACTORIZATION.md`](F2_AFFINE_K1_CARRIER_JET_FACTORIZATION.md).
After the raw carrier centers are scaled by `kappa^j/lambda`, its weighted
triangular untransport removes `P0,Q0,Gamma`: the first four normalized jets
recover `a,b,c,d`, and the final three are explicit compatibility residuals
of weights `5,6,7`.  If `P0,Q0,Gamma` are fixed independently, the
seven-step test has four reconstruction rows and three genuine obstruction
rows.  If they remain free, however, the raw seven-jet Jacobian is
`3*Res(p',q')`, so on the immersed locus those parameters generically absorb
the residuals.  No target-normalization parameter may be silently set to
zero, and the carrier jet alone is then not an exclusion.  On the
`E_6+A_1` subfamily the situation reverses: its normalized carrier jets obey
an exact binomial closed formula and, after all three transports, still form
a codimension-three locus in the seven fixed centers.  Its carrier contact
is therefore an exact finite test once those centers are exposed.  At the
`E_8` endpoint this specializes further to a prime codimension-four locus,
given by four explicit scale-free equations in the raw centers.

<!-- status-consumer: PF2K1JF1 7bc57f390f0531b5 -->

The primitive carrier count sharpens this limitation: the `(3,5)` raw
parameter space and the seven matched centers both have dimension seven, so
the existing fan stops at saturation.  Its first normalization-invariant raw
equation occurs only at jet eight.  On the nonimmersion divisor the generic
four-node packet becomes `A_2+3A_1`; the conductor preimage changes from
eight simple points to a `2+6` split and the raw jet Jacobian has corank one.
This marks the target value at which a boundary attachment must be tested;
it is not itself a localized `ch_2` class.  The unibranch SNC theorem shows
that a minimal two-boundary point contributes `2q_p`, where `q_p` is its
local residue ramification index.  The complementary cusp-fold theorem
shows that a smooth-boundary point has lower Fitting exponent `2q_p-1`, and
its exact unramified fold has point correction one.  Both numbers are
distinct from the transverse divisorial index `e` in (4.9).  The source
incidence, not the target cusp alone, selects the applicable packet.

<!-- status-consumer: PCJDP1 d4c16bb71dfc6b80 -->
<!-- status-consumer: LUAF1 b0279670ffbd3fa5 -->
<!-- status-consumer: LCAD1 7b9c15d3dfae0337 -->

Writing `E^2=-n`, formula (3.4) becomes

\[
 \boxed{
 A_{k=1}(e,f,n,b)=ef(b-7)-\frac{e^2n}{2}.}                \tag{4.9}
\]

For the minimal purity signature `(e,f,n)=(2,1,1)`, the arithmetic expression
is

\[
 \boxed{A_{\min}(b)=2b-16.}                             \tag{4.10}
\]

The puncture theorem determines only the target contact.  If
`lambda!=125/729`, then `b=0`.  If `lambda=125/729`, the curve passes the
first carrier center and `1<=b<=8`; its higher infinity jet determines the
exact value.  Because no affine divisor can be extracted over the resolved
terminal neighborhood, the terminal residue index three does **not** force
the affine transverse index.  Formula (4.9), with `e>=2`, is the applicable
statement in both cases.

## 5. The puncture-compatible degree-floor sieve

Let `s_X` be the number of further smooth source-boundary blowups relative
to the present `27/48` graphs.  They change the global budget by `-s_X/2`.
Assume, only for this calculation, an exact filtration whose only
divisorial quotient beyond the known root class `27` is one affine packet
(4.9), with all remaining correction carried by a finite-length sheaf `Z`.

For the squarefree degree floor `d=6`, the root-subtracted virtual residual
was `-10`.  Therefore

\[
 \boxed{
 \ell(Z_{\rm sq})
 =\frac{e^2n-2ef(b-7)-20-s_X}{2}.}                       \tag{5.1}
\]

Effectivity and integrality force

\[
 \boxed{
 s_X\equiv e^2n\pmod2,\qquad
 s_X\le e^2n-2ef(b-7)-20.}                              \tag{5.2}
\]

For the double-row floor `d=12`, the root-subtracted residual was `17/2`.
The same calculation gives

\[
 \boxed{
 \ell(Z_{\rm dbl})
 =\frac{e^2n-2ef(b-7)+17-s_X}{2},}                        \tag{5.3}
\]

and hence

\[
 \boxed{
 s_X\equiv e^2n+1\pmod2,\qquad
 s_X\le e^2n-2ef(b-7)+17.}                               \tag{5.4}
\]

The opposite parities in (5.2) and (5.4) are a concrete source-chain test.
They arise before any node length is evaluated.

For the smallest arithmetic signature `(e,f,n)=(2,1,1)`, one has
`A=2b-16`, and the formulas reduce to

\[
 \boxed{
 \ell(Z_{\rm sq})=\frac{12-4b-s_X}{2},\qquad
 \ell(Z_{\rm dbl})=\frac{49-4b-s_X}{2}.}                \tag{5.5}
\]

At the smallest parity-compatible values `s_X=0` and `s_X=1`, respectively,
the squarefree benchmark over `b=0,...,8` is
`6,4,2,0,-2,...,-10`, while the double benchmark is
`24,22,20,...,8`.  These numbers become actual lengths only after the exact
filtration and the assumed source signature have been proved.  In
particular, they do not assert that `s_X=0`, `E^2=-1`, or `f=1` occurs.

### 5.1 Conditional cusp-incidence refinement

Suppose the target lies on the generic nonimmersion stratum.  Let its
ordinary cusp have `h` boundary preimages on the normalization of `E`, and
let `c` of those points be actual source-boundary nodes.  If every point is
one of the exact minimal fold/SNC packets, the cusp-attachment dichotomy gives

\[
 \boxed{B_{\rm cusp}=2f-h+c,\qquad f\le B_{\rm cusp}\le2f.} \tag{5.6}
\]

More generally this expression is a lower bound for the sum of isolated
`Fitt_1` colengths.  Equality as a finite `ch_2` subledger still requires the
assumed exact filtration.  Once it holds, (5.1) and (5.3) must satisfy

\[
\boxed{\begin{aligned}
e^2n-2ef(b-7)-20-s_X&\ge2B_{\rm cusp}
 &&\text{(squarefree)},\\
e^2n-2ef(b-7)+17-s_X&\ge2B_{\rm cusp}
 &&\text{(double)}.
\end{aligned}}                                  \tag{5.7}
\]

Equivalently, after the forced cusp subquotients are booked, the remaining
unidentified finite quotient satisfies

\[
 \ell(Z_{{\rm rest},\rm sq})
 =\frac{e^2n-2ef(b-7)-20-s_X-2B_{\rm cusp}}{2},\tag{5.7a}
\]

\[
 \ell(Z_{{\rm rest},\rm dbl})
 =\frac{e^2n-2ef(b-7)+17-s_X-2B_{\rm cusp}}{2}.\tag{5.7b}
\]

The parity gates are unchanged; the effectivity upper bounds for `s_X`
drop by `2B_cusp`.  Here `e` is still the transverse divisorial index; the
local residue indices `q_p` enter only through the fiber sum and `h`.

For the smallest signature `(e,f,n)=(2,1,1)`, the smooth-fold endpoint
`B_cusp=1` and the boundary-node endpoint `B_cusp=2` give, respectively,

\[
\begin{array}{c|cc}
&\text{smooth fold}&\text{boundary node}\\ \hline
\text{squarefree}&4b+s_X\le10&4b+s_X\le8\\
\text{double}&4b+s_X\le47&4b+s_X\le45.
\end{array}                                                   \tag{5.8}
\]

At the smallest parity-compatible values `s_X=0,1`, the corresponding
smooth-fold remaining-residual rows over `b=0,...,8` are
`5,3,1,-1,...,-11` and `23,21,19,...,7`; the boundary-node rows are
`4,2,0,-2,...,-12` and `22,20,18,...,6`.
This is conditional on the cusp stratum and attachment hypotheses; it is not
an exclusion of the full F2 row.

The same order gate is branchwise for other degenerations.  For a complete
fiber over a unibranch value of multiplicity `m_C`, it gives the general
isolated-Fitting lower bound

\[
 B_C\ge m_C f-h+c\ge(m_C-1)f.                  \tag{5.9}
\]

The stronger value `m_Cf` is the node-saturated endpoint, not the universal
one.  At the monomial `(3,5)` cusp, `m_C=3`; the smallest universal fiber
bound is therefore `2f`, while an all-node fiber has value `3f`.  For the
smallest signature these two endpoints give squarefree/double inequalities

\[
\begin{array}{c|cc}
&B_C=2f&B_C=3f\\ \hline
\text{squarefree}&4b+s_X\le8&4b+s_X\le6\\
\text{double}&4b+s_X\le45&4b+s_X\le43.
\end{array}                                                   \tag{5.10}
\]

In contrast, nodes and the individual branches of an ordinary triple point
are immersive and force no lower bound through this theorem.  This is why
the conserved affine delta/conductor total cannot replace the branchwise
attachment calculation.

## 6. Consequence and next compiler target

The global logarithmic budget now reaches the generic `k=1` affine row: its
divisorial kernel degree is no longer unknown.  The calculation also shows
that the route does not automatically exclude the row.  Several
nonnegative minimal-signature residuals survive, while other contacts are
conditionally impossible in a one-packet squarefree filtration.

The next source calculation has only three leading geometric tasks:

1. evaluate the now-explicit fixed-target carrier eight-jet and four affine
   node fibers after the complete F2 source pair is supplied;
2. locate the boundary divisor away from the resolved terminal neighborhood
   and compute its actual `(e,f,E^2)` signature; and
3. resolve the source pullback far enough to determine `s_X` and all
   boundary points where `Fitt_1` ceases to be the unit ideal.

The fixed-coordinate and affine-local parts of Step 1 are provided by the
[`Keller-pullback theorem`](F2_AFFINE_K1_KELLER_PULLBACK.md): the pullback is
squarefree, and every affine conductor point is an ordinary node with
normalization defect one in one of four explicit fibers.  Only the fiber
counts and the boundary
valuation calculation require the source polynomials.

<!-- status-consumer: PF2K1PB1 6f837229017243c4 -->

The
[`all-stratum conductor theorem`](F2_AFFINE_TARGET_K1_CONDUCTOR_CONSERVATION.md)
also removes the need for a separate total-conductor computation on the
degenerate `k=1` locus: affine delta remains four and the exact conductor
divisor retains degree eight.  Only its distribution among singular fibers
and boundary attachments remains stratum-dependent.

<!-- status-consumer: PF2K1CC1 f152c82ef2d54c32 -->

The
[`tame-node packet theorem`](F2_AFFINE_K1_TAME_NODE_PACKET.md) removes one
possible source of that finite point length.  Every fs tame Kummer toroidal
packet over a resolved target node is log-étale; this remains true for the
collided model `z^e=x*y` after resolving its entire exceptional chain.
Hence its logarithmic cokernel and localized `ch_2` are zero.  The remaining
`s_X` calculation must isolate the non-toroidal part of the completed source
pullback rather than reuse the target conductor length.  Its first finite
test is exact: compute the exponent-matrix rank and, only in rank one, the
two determinant first jets.  A nonzero first jet gives smooth support and
zero normalization mismatch.

<!-- status-consumer: PF2K1TN1 521fb57f7e6abc1f -->

The
[`affine strict-log-étale resolution theorem`](AFFINE_KELLER_STRICT_LOG_ETALE_RESOLUTION.md)
removes all affine singularity types from this relative point ledger.  Every
embedded resolution of the target conductor pulls back strictly étale, so
nodes, cusps, tacnodes, and higher multiple points have zero relative log
cokernel.  Their weighted finite-fiber counts constrain boundary escape but
cannot contribute to `Z`.  Only the unresolved compactification boundary is
left for the rank/first-jet test.  At a unibranch value the local isolated
Fitting lower bound is `q_p*m_C-1+epsilon_p`; the exact SNC packet has
`epsilon_p=1`, while the ordinary smooth-boundary fold realizes
`epsilon_p=0` and point length one.

<!-- status-consumer: PAER1 60eb24b2232d159e -->

Once those are known, (4.9) and the global budget leave only an exact finite
point-length comparison.  The target conductor atlas for charts
`k=2,...,24` remains outside this theorem, but its affine singularity types
will not require separate relative logarithmic packets.

There is now one necessary correction to the one-packet sieve used in this
note.  Every immersed, distinct-image `k=1` collision partition and the
generic `A_2+3A_1` and `2A_2+2A_1` cusp strata have affine target
complement `Z`.  A component meridian with a positive affine fixed-sheet
remainder cannot act transitively by itself, so a connected Keller
normalization needs a second ramified affine target component.  Formulas
(5.1)--(5.7) remain conditional tests for one selected packet, but the
complete filtration on these seven strata needs two affine divisorial
packets and has source-component floors `29/50`.

The `E_6+A_1` stratum is different: its noncyclic complement admits a
transitive degree-six action of cycle type `2+2+1+1`.  At that cusp the
unibranch lower ledger has `m_C=3`, hence costs at least `2f`; its actual
source incidence and the global Chern filtration are the next obstruction.

<!-- status-consumer: PF2K1M1 fafcbb3c2e6ceb2b -->

For the `E_8` endpoint, exact enumeration proves that the same degree-six
cycle type has one conjugacy class, with image `A_5`.  Hence its Chern
compiler has exactly two distinct `(e,f)=(2,1)` packets: the preferred
longitude preserves each transposition orbit and rules out one `(2,2)`
packet.  If `N` is the sum of their self-intersection magnitudes, their
combined multiplicity-three residual condition is
`28+4N-8b-s_X>=0`.

More generally, the exhaustive `A_5` coset atlas has fixed-sheet degrees
`d=6,10,15,30` and respectively `r=2,4,6,14` distinct `(2,1)` packets.
The squarefree doubled residual is
`7d-62+4N-4r(b-6)-s_X`; on the double row it is
`7d-67+4N-4r(b-6)-s_X`.  For minimal `N=r` and maximal `b=8`, only degrees
six and ten are negative in the squarefree ledger, so Chern effectivity
alone does not close the higher icosahedral actions.

<!-- status-consumer: PF2K1E8M1 bbb282c6bcfa62fc -->

The universal `M^2=1` orbifold atlas also permits nontrivial central gluing.
Its fixed-sheet F2 degrees are
`6,10,12,15,20,24,30,40,60,120`; peripheral rows have `f=1,2,4`.  If `q`
is the number of actual rows and `R` their total residue degree, the formulas
become `7d-62+4N-4R(b-6)-s_X` and
`7d-67+4N-4R(b-6)-s_X`, with component floors `27+q` and `48+q`.
Minimal maximal contact is negative only at squarefree degrees `6,10,12`.

<!-- status-consumer: PF2K1E8O1 4251750ed4e43c89 -->

The later complete-chain theorem shows that `N` is not an independent stable
variable.  Combining the affine cycle with the extraction-root cycle and
all node matching cancels `d,b,s_X`, and the strict-transform squares,
leaving point budget exactly `u-1`.  Every simple-inertia E8 row has cusp
lower `2R>u-1`; any survivor therefore needs a negative noncyclic correction
of the exact deficit listed in that theorem.

<!-- status-consumer: PF2K1CB1 5cc386dba344a867 -->

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_affine_k1_log_ch2.py
```
