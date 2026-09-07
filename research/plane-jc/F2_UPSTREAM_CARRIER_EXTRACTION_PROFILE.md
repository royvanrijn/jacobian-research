# F2 upstream carrier-extraction logarithmic profile

> **Status.**  This note proves the local logarithmic profile at every toric
> node of the six-blowup carrier-extraction skeleton forced by the certified
> F2 `(75,125)` common edge.  The carrier-zero ladder maps unimodularly to the already
> extracted target fan, so all of its boundary nodes are log-etale.  At the
> other end, where the first exceptional divisor meets the strict line at
> infinity, the logarithmic cotangent cokernel is
> `R/(W^3 U^18)`.  Its canonical branchwise splitting has finite quotient
> `R/(W^3,U^18)` of length `54`.  This is a nonzero degree-one matching
> defect, not finite-support torsion in the cokernel.  It does not by itself
> contradict Keller geometry, exclude `(75,125)`, or prove `JC(2)`.

The subsequent
[`cyclic boundary blowup theorem`](LOG_CYCLIC_BOUNDARY_BLOWUP_CONSERVATION.md)
shows that the raw length `54` must be combined with component
self-intersections.  The resulting birationally stable Cartier charge is
`D_root^2/2=27`.  The subsequent
[`kernel-line theorem`](LOG_CYCLIC_COKERNEL_TWIST.md) proves that the actual
cyclic cokernel contribution is `deg(K_root)+27`, so the remaining twist
problem is exactly the degree of the restricted logarithmic kernel line.  The
[`contracted-packet theorem`](LOG_KERNEL_GAUSS_DEGREE.md) makes this
`27-e_root<=27`.  The
[`tangential-coordinate theorem`](LOG_TANGENTIAL_KERNEL_TRIVIALIZATION.md)
then uses the full divisibility `f^*z=W^3U^18*unit` to prove `e_root=0`.
The exact cyclic root contribution is therefore `27`.
<!-- status-consumer: LCBBC1 b3eb4679f781c55f -->
<!-- status-consumer: LCCT1 2fc6ecea7a7c8b49 -->
<!-- status-consumer: LKGD1 8a357250b5005186 -->
<!-- status-consumer: LTKT1 32ac27318f16c20c -->

The exact lattice, support, differential-order, and matching-length checks
are replayed by
[`verify_f2_upstream_carrier_extraction.py`](../scripts/verify_f2_upstream_carrier_extraction.py).
The reusable cyclic-SNC matching calculation is in
[`log_node_profiles.py`](../jcsearch/log_node_profiles.py).

## 1. Certified source support

Write the original affine coordinates as `x,y`.  The common carrier edge is

\[
 P=C_R^3+\text{higher carrier order},\qquad
 -Q=\frac95C_R^5+\text{higher carrier order},
 \tag{1.1}
\]

where

\[
 C_R=x(xy^5-1)^2R(xy^5),\qquad R(0)\ne0.       \tag{1.2}
\]

Consequently every source exponent `(i,j)` satisfies

\[
 -5i+j\ge-15\quad(P),\qquad
 -5i+j\ge-25\quad(Q),                           \tag{1.3}
\]

together with `i+j<=75` or `i+j<=125`.  For the monomial valuations

\[
 \nu_k(x)=-k,\qquad \nu_k(y)=1,qquad 1\le k<5,
\]

these inequalities have the unique minimizers `(3,0)` for `P` and `(5,0)`
for `-Q`.  Hence

\[
 \nu_k(P)=-3k,\qquad \nu_k(-Q)=-5k.             \tag{1.4}
\]

At `k=5` the complete carrier edge appears.  Thus no omitted lower Laurent
coefficient can change the leading orders used below.

## 2. The carrier-zero ladder is unimodular

On the carrier-zero side put

\[
 a=x^{-1},\qquad v=xy^5=\frac{y^5}{a}.
\]

The five boundary rays from the first carrier-side exceptional through the
carrier are

\[
 E_k=(k,1),\qquad 1\le k\le5,                   \tag{2.1}
\]

in `(ord(a),ord(y))` coordinates.  The nonboundary endpoint is `(0,1)`.

Use the target functions

\[
 \pi=\frac{P^3}{(-Q)^2},
 \qquad \omega=h-\text{the seven certified target shears},
 \qquad h=\frac{P^5}{(-Q)^3}.                    \tag{2.2}
\]

At `v=0`, the carrier Wronskian calculation gives

\[
 \pi=a\cdot\text{unit},\qquad
 \omega=a^7y\cdot\text{unit}.                  \tag{2.3}
\]

The exponent map is therefore

\[
 M_0=\begin{pmatrix}1&0\\7&1\end{pmatrix},
 \qquad \det M_0=1.                              \tag{2.4}
\]

It sends the complete source ladder to

\[
 (0,1),(1,8),(2,15),(3,22),(4,29),(5,36).       \tag{2.5}
\]

These are consecutive rays of the certified target fan.  Every adjacent
source and target cone is regular, and the local logarithmic exponent matrix
has determinant one.  Thus every boundary node on this side of the
carrier-extraction packet is log-etale.  No additional source or target
blowup is required.

## 3. Coordinates at the extraction root

The first blowup is centered at `[1:0:0]` on the original line at infinity.
At the node between its exceptional divisor and the strict line at infinity,
use

\[
 W=\frac yx,\qquad U=\frac1y,
 \qquad x=(WU)^{-1},\quad y=U^{-1}.              \tag{3.1}
\]

Here `W=0` is the first exceptional divisor and `U=0` is the strict line at
infinity.  The degree bounds and carrier halfspaces imply

\[
 P=W^{-15}U^{-75}A(W,U),\qquad
 -Q=W^{-25}U^{-125}B(W,U),                       \tag{3.2}
\]

where `A` and `B` are units at `(0,0)`.  Indeed `(15,60)` and `(25,100)`
are the unique carrier-edge endpoints of maximal `x`-degree; every other
allowed monomial acquires positive `W` or `U` order in (3.2).

It follows that

\[
 \pi=W^5U^{25}\cdot\text{unit}.                 \tag{3.3}
\]

The carrier-infinity profile in `PF2CLP1` gives

\[
 \operatorname{ord}_{W=0}(\omega)=36.            \tag{3.4}
\]

## 4. The exact two-form forces the missing exponent

Target shears do not change the wedge with `pi`.  Since the exponent change
from `(P,-Q)` to `(pi,h)` has determinant one and the Keller Jacobian is a
nonzero constant,

\[
 d\pi\wedge d\omega
 =\text{unit}\cdot\frac{P^7}{(-Q)^6}\,dx\wedge dy.
 \tag{4.1}
\]

Equation (3.1) gives

\[
 dx\wedge dy=\text{unit}\cdot W^{-2}U^{-3}
 dW\wedge dU.
\]

Substitution in (4.1) yields the exact order

\[
 d\pi\wedge d\omega
 =\text{unit}\cdot W^{43}U^{222},dW\wedge dU.  \tag{4.2}
\]

Absorb the unit in (3.3) into `W`, so that `pi=W^5U^25`.  On a monomial
`W^rU^s`, the operator `d pi wedge d(-)` is a nonzero scalar times

\[
 W^{r+4}U^{s+24},dW\wedge dU                    \tag{4.3}
\]

unless `s=5r`, in which case the monomial is a function of the primitive
target normal and lies in the shear kernel.  Because `ord_W(omega)=36`,
(4.2)--(4.3) force the leading kernel monomial

\[
 \omega=W^{36}U^{180}\cdot\text{unit},           \tag{4.4}
\]

while the first nonshear term has exponent

\[
 (r,s)=(43-4,222-24)=(39,198).                   \tag{4.5}
\]

The difference between (4.5) and (4.4) is `(3,18)`.

## 5. Local logarithmic normal form

The primitive target normal and carrier residue are

\[
 T=\frac{\pi^{29}}{\omega^4},\qquad
 \zeta=\frac{\omega^5}{\pi^{36}}.                \tag{5.1}
\]

Equations (3.3) and (4.4) give

\[
 T=WU^5\cdot\text{unit}.                         \tag{5.2}
\]

After subtracting a harmless local function of `T` from `zeta-zeta(0,0)`,
equation (4.5) gives a target tangential coordinate

\[
 z=W^3U^{18}\cdot\text{unit}.                    \tag{5.3}
\]

The exponent determinant is

\[
 \det\begin{pmatrix}1&5\\3&18\end{pmatrix}=3.  \tag{5.4}
\]

The target has only the smooth boundary `T=0`, while both `W=0` and `U=0`
are source-boundary components.  Relative to the target basis
`(dlog(T),dz)` and source basis `(dlog(W),dlog(U))`, one matrix entry is a
unit and the determinant is

\[
 3W^3U^{18}\cdot\text{unit}.                    \tag{5.5}
\]

Smith reduction over the completed regular local ring therefore gives

\[
 \boxed{\operatorname{coker}(\theta_f^{\log})
       \simeq R/(W^3U^{18}).}                    \tag{5.6}
\]

In particular `Fitt_1=R` and `Fitt_0=(W^3U^18)`.

## 6. The finite branch-matching quotient

The ideals `(W^3)` and `(U^18)` have intersection `(W^3U^18)`.  Hence there
is a canonical exact sequence

\[
0\longrightarrow R/(W^3U^{18})
 \longrightarrow R/(W^3)\oplus R/(U^{18})
 \longrightarrow R/(W^3,U^{18})\longrightarrow0. \tag{6.1}
\]

The final quotient has monomial basis

\[
 W^iU^j,\qquad 0\le i<3,quad0\le j<18,
\]

and therefore

\[
 \boxed{\operatorname{length}R/(W^3,U^{18})=54.} \tag{6.2}
\]

The cyclic module in (5.6) is a hypersurface Cohen--Macaulay module, so
`H^0_{(W,U)}=0`.  The finite quotient (6.2) is instead the thickened
branch-matching defect that injects into degree-one local cohomology, exactly
as predicted by `LCDS1`.  Calling (6.2) a degree-zero torsion class would be
incorrect.

## 7. Consequence and remaining claim boundary

The upstream carrier-extraction calculation has two different outcomes:

1. the carrier-zero ladder is entirely log-etale and contributes no defect;
2. the extraction-root node contributes a forced nonzero matching quotient
   of length `54`.

This supplies the first explicit nonzero class required by the full nodal
programme, but only in degree-one local cohomology.  A contradiction still
requires an independent global Keller theorem that annihilates, cancels, or
cannot accommodate this exact class.  No such theorem is proved here.

The subsequent affine-purity theorem resolves the generic existence and
counting part of the first task: it forces a new component and raises the
source floors to `28/49`.  It also proves that the coarse ledger is
underdetermined at every remaining geometric degree.  The remaining
geometric and global tasks are therefore:

1. recover the affine nonproperness curve, factor its pullback, and locate
   the proximity chain of the forced ramification component;
2. compile any additional global resolution centers or non-toric special
   points on the extraction components;
3. place the exact cyclic root class `27` in a global
   localized-second-Chern or conductor identity and determine whether other
   nodes can cancel it; and
4. independently complete or exclude the lower Laurent coefficient rows.

The subsequent
[`outgoing terminal-tail theorem`](F2_OUTGOING_TERMINAL_TAIL.md) proves that
the omitted tail maps unimodularly to the existing target fan and contributes
no additional logarithmic cokernel.

<!-- status-consumer: PF2OTT1 af25012e34020e11 -->

<!-- status-consumer: PF2APF1 192055eb737d3140 -->

The source graphs need no additional blowups for this theorem itself.  The
later purity frontier supplies a different component, giving the current
global lower bounds `28/49`.

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_upstream_carrier_extraction.py
```
