# All-degree q-graph obstruction from a square on one base axis

> [HC4MYGL1/2](CRITICAL_COORDINATE_JET.md) subsequently closes all
> polynomial base and relative dual automorphisms in this chart by a
> critical first jet and an imported line-injectivity theorem. The
> elementary square-axis result below retains its smaller hypotheses
> and needs no such additional import.

14 November 2026. **HC4MYGQ2** closes every standard polynomial graph
eliminating q. It permits arbitrary relative affine changes of the retained
p,r variables and every polynomial base automorphism whose old y coordinate
is affine in the parameter base. No degree bound is imposed.

Unlike the [curve-trace method](CURVE_TRACE.md), this proof does not need
the graph to disappear on a nonconstant base curve. The decisive restriction
is a polynomial square plus an affine dual term on y=0. **HC4MYGQ1** is
the elementary affine-plane classification; the global conclusion uses
only the reproduced [Tang alternative](../tang-hc4-audit/STRUCTURAL_REVIEW.md).
No Rentschler theorem or JC2 assumption is needed for this q-chart result.

## 1. A general square-axis criterion

Let f(x,y,P,Q) be polynomial over a characteristic-zero field, with

\[
 f(x,0,P,Q)=L a(x,P,Q)^2+M a(x,P,Q)+2Nx r(x,P,Q),   \tag{1}
\]

where LN is nonzero, a is polynomial, r is affine in P,Q with polynomial
x coefficients, and a(0,P,Q), r(0,P,Q) are an affine coordinate pair.
After a constant affine dual change we may normalize

\[
 a(0,P,Q)=P,\qquad r(0,P,Q)=Q.                    \tag{2}
\]

**HC4MYGQ1.** If f has a jointly affine constant two-plane and its
Hessian determinant is a nonzero constant, its gradient is polynomially
invertible. This is proved below without importing Tang.

**HC4MYGQ2.** Every f satisfying (1)-(2) with nonzero constant Hessian
determinant has a polynomially invertible gradient. A hypothetical
nonautomorphic gradient would have the plane in HC4MYGQ1 by HCTANG2.
The standard characteristic-zero descent applies after proving this
over C. Constant affine parameter changes preserve both conclusions.

## 2. The polynomial-square degree lemma

Suppose S(x,U,V) is polynomial, S(0,U,V)=ell(U,V) is affine, and

\[
 S^2=\ell^2+xH(U,V)+xE(x,U,V),                    \tag{3}
\]

where E is affine in U,V. Then H is affine in U,V.

Indeed, if d=deg_(U,V) H>2, the right side has dual degree d and
d=2m with m=deg_(U,V) S>=2. Its highest homogeneous part satisfies
S_m^2=xH_d. But S_m vanishes at x=0, so S_m is divisible by x and
its square by x^2, a contradiction since H_d is independent of x.

Thus d<=2 and S is affine in U,V. Write its linear part as
s1(x)U+s2(x)V. Comparing quadratic dual parts of (3) shows that each
si(x)^2 has x degree at most one. Each si is therefore constant.
Their values at x=0 already give ell's linear part, so H's quadratic
part is zero. This proves the lemma in all degrees, including zero
or constant H. The affine bound on E is essential to this argument.

## 3. Exhaust the affine planes by base-projection rank

Let W be a jointly affine constant two-plane. Its projection onto the
(x,y) base has rank 0, 1 or 2. All source changes used in this proof
are constant affine changes, so unit Hessian remains unit.

### Rank zero

Then W is the whole pure-dual plane, and f is globally affine in P,Q.
At x=y=0, however, f=LP^2+MP by (2), a contradiction.

### Rank one

W contains a nonzero pure-dual direction. Its second derivative at
the origin base slice LP^2+MP forces that direction to be partial_Q.
Thus f_QQ=0 globally. On y=0, the term 2Nx r in (1) is affine in Q.
Polynomial Q degree then forces a to be independent of Q: otherwise
the square has degree twice a's positive Q degree, which neither
M a nor the affine term can cancel. Consequently

\[
 f_Q(x,0,P,Q)=C_0(x):=2Nx k(x),\qquad k(0)=1,      \tag{4}
\]

where k is the Q coefficient of r. In particular C0(0)=0 and
C0'(0)=2N is nonzero.

If the other plane generator moved x but not y, its mixed derivative
with partial_Q, restricted to y=0, would force C0'=0, a contradiction.
Therefore normalize that generator to (eta,1,tau,sigma). Constant-flow
transport from y=0 gives f_Q=C0(x-eta y), independent of P.

In coordinates S=y, T=x-eta y, U=P-tau y, V=Q-sigma y, joint affineness
now gives

\[
 f=C_0(T)V+S F(T,U)+B(T,U),\qquad
 \det\operatorname{Hess}f=(C_0'(T)F_U)^2.
\]

Unit determinant makes C0' and F_U nonzero constants. Equation (4)
fixes C0(T)=2NT, while F=alpha U+beta(T), alpha nonzero constant.
This triangular cotangent form has the explicit inverse already used
in HC4MYGF1: recover T from f_V, U from f_S, S from f_U, and V from f_T.

### Rank two

Subtract constant base-linear terms from P,Q to make W the base plane.
Write the resulting coordinates as U,V; they agree with P,Q at the
base origin. Then

\[
 f=LU^2+MU+xF(U,V)+yG(U,V),\qquad
 \det\operatorname{Hess}f=\operatorname{Jac}(F,G)^2. \tag{5}
\]

Restrict to y=0. The transformed a is still polynomial and has
a(0)=U; the transformed r is still affine in U,V. Comparing (1) and
(5) gives

\[
 (2La+M)^2=(2LU+M)^2+4LxF-8LNx r.                \tag{6}
\]

Apply the degree lemma to S=2La+M, H=4LF and E=-8LN r. It follows
that F is affine in U,V. The constant nonzero Jacobian in (5) makes
F nonconstant. Choose a complementary affine coordinate to F; then
G has a nonzero constant derivative in that coordinate, so (F,G) is
a triangular polynomial automorphism. Recover U,V from F,G and x,y
from the other two gradient outputs using D(F,G)^T. This proves
HC4MYGQ1 and hence HC4MYGQ2.

## 4. Applying the criterion to the graph eliminating q

For the scaled Meng--Yang family, put q=R(x,y,P,Q), and allow any
relative affine dual automorphism (p,r)=phi(x,y,P,Q). On y=0 the
defining linear forms satisfy exactly

\[
 A=p+3xR-x^3r,\qquad B=2xr.
\]

Thus the pullback has (1), with a=A|y=0. At the base origin its
(a,r) pair is the affine coordinate pair (p,r), giving (2) after
constant affine normalization. All polynomial degrees and all affine
dual coefficient degrees are unrestricted.

The same proof permits a polynomial base automorphism beta(s,t) when
its old y coordinate is affine in s,t. Normalize this affine coordinate
to t by a constant affine base change. Since beta is an automorphism,
its other coordinate is x=c*s+H(t), c a nonzero constant. A further
constant affine change in s makes x=s on t=0. The square-axis hypothesis
then holds on that straight parameter axis. This includes arbitrary
nonlinear triangular base changes preserving the y-axis coordinate;
it does not include every nonlinear base automorphism.

The determinant hypothesis always concerns the actual pulled-back
potential. It is not inferred from a determinant before reparameterization.
The native q-chart conclusion is gradient invertibility, not absence
of every unmarked unit-Hessian graph.

## Scope, controls and reproduction

All three standard dual-elimination graph charts p,q,r are now excluded
as gradient-counterexample routes in their native parameters. Their
reparameterization extensions have different scopes: the r chart permits
every base and relative dual automorphism; the p and q results permit
all relative affine dual changes, with the additional q base restriction
just stated. The later critical-jet result covers arbitrary polynomial
base and retained-dual automorphisms in the q chart. The later
[HC4MYGB2](BASE_ELIMINATED_GRAPHS.md) also closes both native
base-eliminated charts. Nonlinear reparameterizations of those charts
and genuinely nonlinear base/dual mixing remain separate.

The proof cannot simply assume r affine on the axis when it is not.
For example f=LU^2+MU+2Nx(V+U^2)+yU has a unit Hessian and (1) with
a=U, r=V+U^2. Its x coefficient is quadratic, not affine. Its gradient
is nevertheless invertible. This is a control on the degree lemma's
hypothesis, not a counterexample to the new obstruction.

The [checker](verify_square_axis.py) verifies the square identity,
normal forms, representative inverses, affine-coordinate normalization,
and the explicit graph construction in the companion
[p-chart endpoint note](P_GRAPH_SINGLETON_ENDPOINT.md).
The unbounded degree lemma and plane exhaustion are written proofs.
[Actual receipt](receipt-square-axis.json).

```sh
timeout 60s .venv/bin/python research/jacobian-followup/hc4-triangular-graphs/verify_square_axis.py
```

No plane degree bound, JC2 assumption, formal verification, independent
external review or novelty priority is claimed. The only external
classification used for the q-chart conclusion is HCTANG2.
