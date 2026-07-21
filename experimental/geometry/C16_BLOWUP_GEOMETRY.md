# C16: explicit blow-up geometry

Work over an algebraically closed field of characteristic zero.  Let `H` be
an admissible seed and put

\[
 E=H(W)-BCW+cAC^2,\qquad
 \gamma={BC-H'(W)\over c},
\]

on the marked-root incidence `I_H=V(E)`.  The reconstruction functions are

\[
 x={C\over\gamma},\qquad y={W-\gamma\over C},\qquad
 z={\gamma N\over b_0C^2},\qquad
 N=\gamma^2+(a_0-1)\gamma-a_0W.                 \tag{1}
\]

This note resolves the rational reconstruction map, computes every boundary
valuation and records the precise sense in which the resolution is minimal.
It also corrects a terminology error in the earlier Kummer description:
`mu` and `m_0-1` are ramification indices, not residue degrees of the
dicritical prime divisors.

## 1. The canonical graph blow-up

After embedding the source as `[1:x:y:z]` in `P^3`, (1) is represented on
`I_H` by

\[
 [j_0:j_1:j_2:j_3]
 =[\gamma C^2:C^3:\gamma C(W-\gamma):\gamma^2N/b_0].              \tag{2}
\]

Let `J=(j_0,j_1,j_2,j_3)`.  The closure of the graph is
`Bl_J(I_H)`, and the resolved inverse graph is

\[
 X_H^{\rm gr}=\operatorname{Nor}(\operatorname{Bl}_J(I_H)).       \tag{3}
\]

Thus the blow-up center is the base scheme `V(J)`, not merely its reduced
support.  Its four standard charts are explicit:

| chart | regular ratios |
|---|---|
| `j_0 != 0` | `(j_1/j_0,j_2/j_0,j_3/j_0)=(x,y,z)` |
| `j_1 != 0` | `(j_0/j_1,j_2/j_1,j_3/j_1)=(1/x,y/x,z/x)` |
| `j_2 != 0` | `(j_0/j_2,j_1/j_2,j_3/j_2)=(1/y,x/y,z/y)` |
| `j_3 != 0` | `(j_0/j_3,j_1/j_3,j_2/j_3)=(1/z,x/z,y/z)` |

Equation (3) is canonical and minimal among normal modifications on which
reconstruction extends: any such modification makes the pullback of `J`
invertible, hence factors uniquely through `Bl_J(I_H)`, and normality gives
the factorization through its normalization.  This is the intrinsic global
answer.  The following local blow-ups describe its divisorial valuations and
give smooth transverse models.

## 2. The discriminant divisor

On `C != 0`, the incidence is smooth because `E_A=cC^2` is a unit.  At the
generic point of `E=E_W=0`, use `u=E_W` as a transverse parameter.  Since
`gamma=-u/c`, for a generic repeated root `r != 0` one obtains

\[
\begin{array}{c|rrrr}
 &\gamma&C&x&y\\ \hline
 \operatorname{ord}_{D_\Delta}&1&0&-1&0
\end{array}
\]

and

\[
 \operatorname{ord}_{D_\Delta}(z)=
 \begin{cases}1,&a_0\ne0,\\2,&a_0=0.\end{cases}                 \tag{4}
\]

The projective reconstruction limit is `[0:1:0:0]`.  No blow-up of the
incidence is required here: `D_Delta` is a polar divisor, not an exceptional
divisor.  It maps birationally to the discriminant by

\[
 B={H'(r)\over C},\qquad
 A={rH'(r)-H(r)\over cC^2},                       \tag{5}
\]

and the finite root cover has transverse ramification index two.  Its
discrepancy is zero.

## 3. A nonzero primitive root

Let `rho != 0,1` be a root of multiplicity `mu`.  Put `w=W-rho` and make the
exact coordinate change

\[
 b=B(\rho+w)-cAC.
\]

Then `E=H(rho+w)-Cb`.  After an etale unit extraction
`H(rho+w)=q^mu`, the completed local model is

\[
                         Cb=q^\mu,                              \tag{6}
\]

with `A` a smooth parameter.  It is `A_(mu-1) x A^1`.  For `mu=1` it is
smooth.  For `mu>=2`, its minimal toric resolution has exceptional divisors

\[
 E_{\rho,1},\ldots,E_{\rho,\mu-1},
\]

obtained by inserting the rays `(i,mu-i)`, `1<=i<=mu-1`, in the cone of
(6).  The two affine toric charts adjacent to `E_(rho,i)` use the consecutive
rays `(i-1,mu-i+1),(i,mu-i)` and
`(i,mu-i),(i+1,mu-i-1)`.  Their monomial coordinate maps are determined by

\[
 \operatorname{ord}_{E_{\rho,i}}(C)=i,\quad
 \operatorname{ord}_{E_{\rho,i}}(b)=\mu-i,\quad
 \operatorname{ord}_{E_{\rho,i}}(q)=1.             \tag{7}
\]

Explicitly, the resolution is covered by `mu` smooth charts
`U_i=Spec k[A,u_i,v_i]`, `0<=i<=mu-1`, with

\[
 C=u_i^{i+1}v_i^i,\qquad
 b=u_i^{\mu-i-1}v_i^{\mu-i},\qquad q=u_iv_i.       \tag{7a}
\]

The boundary `v_0=0` is `D_b`, the boundary `u_(mu-1)=0` is `D_rho`, and
the intermediate coordinate divisors glue to the `E_(rho,i)`.  This is a
sequence of normalized toric blow-ups.  Its first center is the singular line
`V(C,b,q)` (with coordinate `A`); each later center is the remaining invariant
singular line.  On a transverse surface slice the centers are points.

On every transverse surface slice the divisor intersections are

```text
D_b -- E_(rho,1) -- ... -- E_(rho,mu-1) -- D_rho,
```

where each exceptional curve has self-intersection `-2`, adjacent curves
meet once transversely and all other pairs are disjoint.  In the threefold
the corresponding divisors meet along the smooth `A`-lines.  This is the
Du Val resolution, so

\[
 a(E_{\rho,i};I_H)=0.                              \tag{8}
\]

Write

\[
 g_i=\min(2i,\mu-1).
\]

Because `B=(b+cAC)/(rho+q)`, direct valuation of (1) gives

\[
\begin{array}{c|rrrrr}
 &C&W&\gamma&x&y\\ \hline
 E_{\rho,i}&i&0&g_i&i-g_i&-i
\end{array}                                                       \tag{9}
\]

and

\[
 \operatorname{ord}_{E_{\rho,i}}(z)=
 \begin{cases}g_i-2i,&a_0\ne0,\\2g_i-2i,&a_0=0.\end{cases}      \tag{10}
\]

Every `E_(rho,i)` maps to the target line `B=C=0` (with `A` free), so none
is dicritical.  The strict divisor `D_rho` has

\[
 (\operatorname{ord}C,\operatorname{ord}\gamma,
   \operatorname{ord}x,\operatorname{ord}y)
 =(\mu,\mu-1,1,-\mu),                              \tag{11}
\]

and `ord(z)=-(mu+1)` if `a_0!=0`, while `ord(z)=-2` if `a_0=0`.
It maps onto `C=0` with residue degree one and ramification index `mu`.
The `mu` nearby escaping sheets are the ramified sheets of this one divisor.

The total transform of `gamma=0` records its intersections with the
discriminant:

\[
 \operatorname{div}(\gamma)=D_\Delta+(\mu-1)D_\rho
   +\sum_{i=1}^{\mu-1}\min(2i,\mu-1)E_{\rho,i}+\cdots .           \tag{12}
\]

## 4. The zero cluster

Write `H(W)=W^mL(W)`, `L(0)=h_0!=0`.  Blow up the smooth ambient center

\[
                         Z_0=V(C,W).                              \tag{13}
\]

The strict transform has the following two charts.

* In the `C`-chart, `W=CR`, and after removing `C^2` its equation is

  \[
  C^{m-2}R^mL(CR)-BR+cA=0.                        \tag{14}
  \]

  For `m=2` the special fiber is
  `h_0R^2-BR+cA=0`, giving two generic finite branches.  For `m>=3` it is
  `R=cA/B`, giving one finite branch.

* In the `W`-chart, `C=WT`, and after removing `W^2` its equation is

  \[
  W^{m-2}L(W)-BT+cAT^2=0.                          \tag{15}
  \]

  Put `b'=B-cAT`.  If `r=m-2`, etale unit extraction changes (15) into

  \[
                           Tb'=q^r.                              \tag{16}
  \]

For `m=2` there is no escaping component.  For `m=3`, (16) is smooth.  For
`m>=4`, resolve the `A_(m-3)` singularity by inserting the rays
`(i,r-i)`, `1<=i<=r-1`.  The exceptional divisors `F_i` satisfy

\[
 \operatorname{ord}_{F_i}(T)=i,\quad
 \operatorname{ord}_{F_i}(b')=r-i,\quad
 \operatorname{ord}_{F_i}(q)=1.                   \tag{17}
\]

Explicitly, the `r` smooth charts `V_i=Spec k[A,u_i,v_i]`,
`0<=i<=r-1`, are

\[
 T=u_i^{i+1}v_i^i,\qquad
 b'=u_i^{r-i-1}v_i^{r-i},\qquad q=u_iv_i.          \tag{17a}
\]

After the initial center (13), the first toric center is `V(T,b',q)` and
each later center is the remaining invariant singular line.  On a transverse
surface these are point centers.  The endpoints in (17a) are `F_fin` and
`D_0`; the intermediate coordinate divisors are the `F_i`.

The transverse dual graph is

```text
F_fin -- F_1 -- ... -- F_(m-3) -- D_0.
```

The `F_i` are `-2` curves on a transverse surface; adjacent components meet
once and nonadjacent components are disjoint.  The first blow-up has ambient
codimension two and the incidence has multiplicity two along `Z_0`, so
adjunction gives coefficient `1-2=-1` on both components of its exceptional
divisor.  Resolution of (16) is crepant over the first transform, while
`div(q)` pulls back with multiplicity one along every `F_i`.  Consequently,
relative to the original Gorenstein incidence,

\[
 a(D_0;I_H)=a(F_{\rm fin};I_H)=a(F_i;I_H)=-1.       \tag{18}
\]

For `1<=i<=m-3`, put `g_i=min(2i+1,m-1)`.  In the original coordinates,

\[
\begin{array}{c|rrrrr}
 &C&W&\gamma&x&y\\ \hline
 F_i&i+1&1&g_i&i+1-g_i&-i
\end{array}                                                       \tag{19}
\]

and

\[
 \operatorname{ord}_{F_i}(z)=
 \begin{cases}g_i-2i-1,&a_0\ne0,\\2g_i-2i-2,&a_0=0.\end{cases}  \tag{20}
\]

Each `F_i` maps to `B=C=0` and is not dicritical.  The escaping divisor
`D_0` has

\[
 (\operatorname{ord}W,\operatorname{ord}C,
   \operatorname{ord}\gamma,\operatorname{ord}x,
   \operatorname{ord}y)
 =(1,m-1,m-1,0,-(m-2)),                            \tag{21}
\]

with `ord(z)=-(m-2)` for `a_0!=0` and `ord(z)=0` for `a_0=0`.  It maps onto
`C=0` with residue degree one and ramification index `m-1`.  The finite
component has `ord(C)=ord(W)=ord(gamma)=1` and regular `x,y,z`; it maps
generically birationally onto `C=0` by `B=cAT`.  For `m=2`, the quadratic
in (14) instead records the two generic finite sheets.  These finite
components are exceptional for the auxiliary blow-up (13), but not polar
and hence not boundary divisors of the normalized reconstruction graph.

## 5. Exhaustion, images and minimality

There is no root at `W=infinity`, since the homogenized pencil restricts
there to the nonzero leading coefficient of `H`.  Over `C!=0`, (1) shows
that only `gamma=0`, hence `D_Delta`, can be polar.  Over `C=0`, every center
lies over a root of `H`; Sections 3 and 4 resolve all such roots, while the
simple root `W=1` has `gamma=1` and is the regular `x=0` chart.  The list is
therefore exhaustive.

The divisorial images are:

| divisor | target image | residue degree | ramification index |
|---|---:|---:|---:|
| `D_Delta` | `V(Q_H)` | `1` | `2` generically |
| `D_rho` | `V(C)` | `1` | `mu` |
| `D_0` (`m>=3`) | `V(C)` | `1` | `m-1` |
| `E_(rho,i)`, `F_i` | `V(B,C)` | not generically finite | -- |
| `F_fin` | `V(C)` | `1` generically | `1` |

There are three distinct minimality statements.

1. `X_H^gr` in (3) is the unique minimal normal graph modification by the
   universal property of blowing up.
2. The `A`-chains in (6) and (16) are the unique minimal smooth resolutions
   on transverse normal surface slices.  Their exceptional curves are all
   `-2`, so none can be contracted while retaining smoothness.
3. The blow-up (13) is necessary to separate the finite zero-cluster branch
   from the escaping valuation.  The remaining `A_(m-3)` chain is then
   forced by surface minimality.

A globally minimal smooth resolution of an arbitrary threefold is not a
canonical notion.  C16 therefore claims the canonical minimal normalized
graph and the minimal transverse toric resolutions, not an unspecified
global MMP-minimal log resolution.  Further blow-ups may be used to make the
entire union of strict discriminant and boundary divisors simple normal
crossings, but they introduce no new dicritical valuations.

## Executable audit

Run

```bash
python scripts/verify_c16_blowup_geometry.py
```

It verifies (2), the two blow-up charts (14)--(15), both toric valuation
tables, the reconstruction valuations, and the residue-degree/ramification
distinction for bounded multiplicities.  These computations audit the
displayed identities; the all-multiplicity toric argument is the proof.
