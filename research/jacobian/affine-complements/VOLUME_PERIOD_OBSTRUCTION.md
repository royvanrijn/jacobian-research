# A global volume-period obstruction leaves one boundary-pair orbit

14 November 2026. PF2PERIOD1 gives an obstruction at every target pole
height. On N=S\B, the prescribed two-form with divisor E-3B must be
exact if it equals dP wedge dQ for regular functions P,Q. In the full
two-chart boundary family its period is b or a. This excludes the entire
marked-contact family and forces a=0 in the other family. The remaining
b!=0 cases are one complex projective-bundle automorphism orbit.

This is a necessary restriction on any actual compatible quartic
normalization. It does not prove that the remaining pair has a finite
map, nor that any actual quartic normalization admits a ruled completion.

## Full boundary family and two affine charts of N

Use the [complete boundary normal forms](SOURCE_COORDINATES.md), with
q the contact point, and a arbitrary, b!=0. Write Bcal and Ecal for
the Cox equations. The marked point on B is B cap C0.
The charts below cover N, not just the smaller source open U.
Their x,y are local coordinates on N; they are distinct from the global
source coordinates of PF2COORD1.

On t!=0, normalize t=Bcal=1. If q is marked, B column is(s,t^2), so
x=s, y=alpha and beta=1-xy. On s!=0, normalize s=Bcal=1 and take
u=t, v=-beta. The overlap is

\[
 x=u^{-1},\qquad y=u^3v+u.\tag{1}
\]

If q is away from the marked point, B column is(t,s^2). On t!=0 use
x=s, y=beta and alpha=1-x^2y. On s!=0 use u=t, v=-alpha. Then

\[
 x=u^{-1},\qquad y=u^3v+u^2.\tag{2}
\]

Both chart domains are A2 and their overlap is Spec C[u,u^-1,v].
They follow directly from the Cox weights(f,f,h,h-f): changing t=1
to s=1 multiplies alpha by u^-1 and beta by u^-2.
The different cocycle terms u and u^2 must not be interchanged.
In either case dx wedge dy=-u du wedge dv.

## The prescribed rational two-form

Define sigma=f0(x,y) dx wedge dy, where

| Family | f0 on t!=0 | f1 on s!=0, with sigma=f1 du wedge dv |
| --- | --- | --- |
| q marked | b+ax+x^3-x^4y | -bu-a+v |
| q away | b+ax-x^2+x^4y | -bu-a-v |

The transition formulas show f1=-u*f0(1/u,u^3v+u^r) exactly. Each
coefficient is the local equation of E up to a nonzero chart constant.
Thus sigma is regular on N and has zero divisor E intersect N.
At a general point of B, the fibre coordinate y has a simple pole;
the coefficient is linear in y and dy has a double pole. Hence sigma
has pole order three along B. There are no other boundary components:

\[
 \operatorname{div}_S(\sigma)=E-3B.\tag{3}
\]

Any proposed dP wedge dQ with this same divisor is a nonzero constant
multiple of sigma. Their ratio has zero divisor on projective S and
is therefore constant. In particular exactness of sigma is necessary.

## Why a residue on the overlap obstructs a global primitive

Every closed polynomial one-form A dx+B dy on an affine-plane chart
is exact in characteristic zero. An explicit primitive is

\[
 H(x,y)=\int_0^x A(t,0)dt+\int_0^y B(x,t)dt;
\]

the closedness identity B_x=A_y verifies both derivatives.
No analytic convergence or topological comparison theorem is needed.

Choose local polynomial one-form primitives eta0,eta1 of sigma on
the two charts. If sigma=d omega globally on N, then eta_i-omega is
closed and hence dH_i on each affine chart. Therefore eta0-eta1 is
exact on the overlap. If that difference is g(u)du, an overlap primitive
has zero derivative in v and lies in C[u,u^-1]. The derivative of a
Laurent polynomial has zero u^-1 coefficient. Consequently

\[
 \operatorname{Res}_{u=0}(\eta0-\eta1)=0\tag{4}
\]

is necessary. Also dP wedge dQ=d(P dQ) for all regular P,Q on N, so
this criterion applies at all pole heights, without requiring finiteness.

## The two exact residue computations

Take F0=integral(f0 dy), with zero y-independent integration term,
and eta0=-F0 dx. In the marked family choose

eta1=[(bu+a)v-v^2/2]du.

Substitution of(1) gives

\[
 \eta0-\eta1=
 \left(\frac b u+\frac a{u^2}+\frac1{2u^4}\right)du.
 \tag{5}
\]

Its residue is b, which is nonzero in every allowed marked family.
Thus no regular P,Q on N have nonzero dP wedge dQ of divisor E-3B.
In particular no finite quartic map with this volume divisor exists,
at any pole height, when q=B cap C0. The original positive boundary
control belongs to this excluded family, although its abstract A2
complement and canonical classes are valid.

In the away family choose

eta1=[(bu+a)v+v^2/2]du.

Substitution of(2) gives

\[
 \eta0-\eta1=
 \left(b+\frac a u-\frac1{2u^2}\right)du.
 \tag{6}
\]

Its residue is a. Hence a=0 is necessary in this family as well.
The identities are over Q[a,b,b^-1] and retain every parameter.

## The remaining zero-period family is exact, and is not excluded here

For a=0 the following local expressions glue to a global regular
one-form eta on N:

\[
 \begin{aligned}
 \eta|_{t\ne0}
 &=-\left(by-x^2y+\tfrac12x^4y^2+\tfrac12\right)dx,\\
 \eta|_{s\ne0}
 &=\left(buv+b+\tfrac12v^2\right)du.
 \end{aligned}\tag{7}
\]

Their differentials are sigma in their respective charts. Thus the
zero-period condition is not merely an inconclusive computation: the
remaining form really has a global algebraic primitive.
A global primitive is weaker than a decomposition dP wedge dQ with
P,Q regular, and weaker still than a finite degree-four map. Neither
of these stronger properties is constructed or excluded by(7).

## A unique surviving boundary-pair orbit

In the away chart Bcal=alpha*t+beta*s^2 and
Ecal_(a,b)=(as+bt)Bcal-alpha*s^2. For lambda nonzero the Cox change

\[
 (s,t,\alpha,\beta)\longmapsto
 (\lambda s,t,\alpha,\lambda^{-2}\beta)
\]

induces an automorphism of F1 preserving B and q. It satisfies

\[
 Ecal_{a,b}\circ\psi_\lambda
 =\lambda^2 Ecal_{a/\lambda,b/\lambda^2}.\tag{8}
\]

With a=0 choose lambda^2=b over C. Thus every remaining b!=0 pair
is isomorphic, as a pair inside S, to the fixed model

\[
 \boxed{B=(t,s^2),\qquad E=(t^2-s^2,s^2t).}\tag{9}
\]

Their determinant is s^4. This normalization preserves the actual open
N and source open U under the compatible-completion hypothesis; it
is not a birational substitution losing omitted divisors.

## Consequence for the active programme

Any compatible completion of the actual one-omitted-prime quartic packet
must have boundary pair in the orbit(9). The two-omitted-prime packet
already has no compatible completion of the stipulated type.
The separate [degree-transfer theorem](POLE_SEVEN_EXCLUSION.md) continues
to exclude every target pole height through seven on the remaining model.

The next reduced pairs(6,8),(7,8) therefore require analysis
only on this fixed zero-period model, not over the whole boundary family.
The original ambitious goal remains active: higher-pole realizability,
actual compatible embeddings and either full quartic packet are open.

The [producer](compute_volume_period.py) and
[exact identities](volume_period_certificates.json) retain the chart
calculation. The [checker](verify_volume_period.py) checks the Cox chart
transitions, both local volume forms and primitives, the residues, the
remaining global primitive and the normalization identity. The affine
cover, divisor interpretation and global exactness obstruction are
written proofs, without formal verification, independent whole-proof
implementation, external review or a novelty claim.

The subsequent [all-height pole-three theorem](POLE_THREE_ALL_HEIGHTS.md)
excludes every(3,n) presentation, independently of the2022 degree bound.
The current first open reduced pairs are (6,8),(7,8).

The adopted frontier(6,8),(7,8) incorporates the subsequent
[PF2RULED12 full(5,8) exclusion](POLE_FIVE_EIGHT_EXCLUSION.md); it is
not an additional exclusion proved by this earlier note.
