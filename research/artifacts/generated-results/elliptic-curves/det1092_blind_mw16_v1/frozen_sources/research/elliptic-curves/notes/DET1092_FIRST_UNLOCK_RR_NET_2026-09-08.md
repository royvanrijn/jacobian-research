# First 302 unlock: the full next Riemann--Roch net

## Result, with the unresolved arithmetic stated first

**Verified application.** The rejected historical genus-zero bisection is
only a one-dimensional RR solution. Increasing the bound by one fibre gives
a three-dimensional solution space whose restriction at `t=0` is the **whole
pointed-chart chord pencil**. An exact rational transformation identifies its
specialized residual discriminant with the complete saved winning quartic,
up to a rational-function square. The saved coordinate then reconstructs an
independent eighteenth point without reading its saved point coordinates.

**Unknown.** This does not yet explain why302 has additional independent
rational points, predict their square values from the parent, or distinguish
302 from the null fibres. It replaces the incorrect identification of the
pointed chart with one minimal multisection by a verified larger object.
The full mechanism-for-rank-gain goal remains open.

Construction and independent replay:

- [generic-input certificate](../../artifacts/generated-results/elliptic-curves/det1092_first_centre_rr_net_v1.json);
- [independent chart-bridge and rank replay](../../artifacts/generated-results/elliptic-curves/det1092_first_centre_rr_net_replay_v1.json);
- [constructor](../cas/construct_det1092_first_centre_rr_net.sage) and
  [checker](../cas/verify_det1092_first_centre_rr_net.sage).

## Exact construction from the generic sections

**Verified application.** Retain the historical centre word

```text
w=(1,-1,-2,1,0,-1,1,-1,1,-1,1,0,1,0,0,-1,0),  w.G.w=10.
```

Let `T=P_{-w}`. The prior RR calculation gave one line
`A0+A1*x+A2*y=0` through `T`, with bounds `(9,5,3)`, in
`H0(3O+9F)`. Its residual curve `C_min` is nonsplit at zero.

The next calculation uses bounds `(10,6,4)` in `H0(3O+10F)`.
Its exact matrix has size20 by23 and rank20. Its kernel is

\[
\operatorname{span}_{\mathbf Q}\{A,tA,B\},
\]

where the certificate displays every coefficient of `A` and `B`. The new
vector `B` is the first rational-kernel basis vector outside `span(A,tA)`,
made primitive. No splitting condition or exceptional coordinate enters this
choice. The centre choice itself is retrospective, not a prospective policy.

The affine net is

\[
f_i(t;u,v)=B_i(t)+(u+vt)A_i(t),\qquad i=0,1,2.
\]

For the original Weierstrass coefficients, write `L=f0+f1*x`. Its exact
residual quadratic is defined by the polynomial identity

\[
L^2-a_1xLf_2-a_3Lf_2-f_2^2(x^3+a_2x^2+a_4x+a_6)
=(x-x_T)(a x^2+b x+c).
\]

The trace identity ensures divisibility for every net member. Off the zeros
and poles of the displayed coefficients needed by these maps, put

\[
\mathcal D(t;u,v)=b^2-4ac,\qquad
s^2=\mathcal D(t;u,v),\qquad
x=\frac{-b+s}{2a},\quad y=-\frac{f_0+f_1x}{f_2}.
\]

This supplies an explicit equation-side square condition and maps, without
fitting a multisection to the exceptional point. Degenerate parameters need
separate treatment; no universal irreducibility claim is made for all members.

## What the next system adds on the 302 fibre

**Verified application.** The matrix of coefficient values at zero has
rank2 and kernel `span(tA)`. Its image is exactly the two-dimensional space
of line equations through `T(0)`. For a finite literal-coordinate slope `m`,

\[
u(m)=-\frac{B_1(0)+mB_2(0)}{A_1(0)+mA_2(0)}.
\]

The missing denominator-zero value is present in the homogeneous net.
Parameter `v` disappears from the restriction. Thus a chord on302 determines
one projective restriction class, not a unique multisection on the surface.

Let `z` denote the frozen historical search coordinate, and `ell(z)` its
saved rational short-model slope. The model transformation implies
`m(z)=(ell(z)-3)/6`. The checker supplies exact coefficient records for the
Mobius function `u(m(z))` and a rational function `k(z)` and proves the
**whole-function identity**

\[
\mathcal D(0;u(m(z)),v)=k(z)^2F_{\rm saved}(z).
\]

This is independent of `v` and is not an interpolation at one successful
point. Substituting the completed witness `z=-1714/2373` verifies the square
and reconstructs a rational point on literal302. A fresh finite-group
calculation at the18 fixed certificate primes gives rank18 for that point
and the17 specialized generic sections. The replay never reads the stored
exceptional point coordinates. It does use the saved successful chart
coordinate as a **retrospective verification witness**, not as a discovery.

## Genus, and the remaining construction bottleneck

**New deduction from the established divisor formula.** The residual class
in this next system is

\[
D=C_{\min}+F=2O+5F+\phi(w),\quad D^2=2,
\quad p_a(D)=2,\quad D.O=1,\quad D.C_{\min}=0.
\]

The projective line `span(A,tA)` consists of the old bisection plus a fibre.
It must not be counted as a supply of new irreducible curves.
The separately fixed member `B` has a residual discriminant whose squarefree
part is an explicit sextic with nonzero discriminant, coprime to the parent
discriminant. Its normalized double cover therefore has genus2.

**Established literature.** For six simple branch points, Riemann--Hurwitz
gives `2g-2=2*(-2)+6`, hence `g=2`; see
[Stacks, Riemann--Hurwitz](https://stacks.math.columbia.edu/tag/0C1B).
The surface intersection and height framework is reviewed in
[Schuett--Shioda, Elliptic Surfaces](https://arxiv.org/abs/0907.0298).

**New deduction.** The quotient formula at norm10 requires
`b=4+g/2`, so arithmetic genus1 is impossible in this degree-two parity
class. This does not exclude geometric genus1: a singular arithmetic-genus2
member can normalize to genus1. That distinction makes singular members a
meaningful equation-side construction problem. The subsequent
[singular-member gate](DET1092_RR_NET_SINGULAR_MEMBER_GATE_2026-09-08.md)
identifies the dense-open singular-member locus with a genus-nine halving
curve; exceptional members remain unclassified.

**Unknown, not a conjectural theorem.** We have not found an oracle-free
rule selecting a useful low-genus member that splits at zero or a source
of additional non-inherited rational points on the fixed genus2 curve.
The subsequent [single-seed construction](DET1092_SINGLE_SEED_COVER_2026-09-08.md)
does certify an independent eighteenth elliptic section over each of two
existing genus2 cover function fields, and exact splitting at the frozen
control addresses. This does not supply a generic-input successful member
selector. The square condition at zero is
exactly equivalent to the old chart condition; merely rewriting it does not
explain its solubility or exceptional independence. An arithmetic explanation
must add something beyond this equivalence.

## Limits and replay

One20x23 RR kernel and one fixed-member squarefree decomposition were used.
No point search, orbit census, parameter sweep, broad class-group calculation,
or modification of the running eight-fibre V3 pilot was made. Both computations
finished in seconds; no background job is running for this supplement.

```sh
sage -python research/elliptic-curves/cas/verify_det1092_first_centre_rr_net.sage
```

The preceding two specific obstructions remain valid; see the
[initial-unlock note](DET1092_RANK18_BASE_CHANGE_AND_INITIAL_UNLOCK_2026-09-08.md).
