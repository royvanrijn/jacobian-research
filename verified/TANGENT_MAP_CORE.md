# The tangent-map core and its weighted suspension

This note isolates the two-dimensional mechanism inside every admissible
weighted marked-root Keller map.  It consolidates the denominator-free
Jacobian calculation, inverse pencil, discriminant normalization,
reconstruction pole, and Hessian Fitting divisor into one diagram.

Work over a characteristic-zero field `k`.  Let `c in k^*` and let
`H in k[W]` have degree `n>=2`.  Define the **plane tangent map**

\[
 \Phi_H:\mathbb A^2_{W,\gamma}\longrightarrow\mathbb A^2_{s,t},
 \qquad
 (W,\gamma)\longmapsto
 \bigl(H'(W)+c\gamma,\,
 W(H'(W)+c\gamma)-H(W)\bigr).
 \tag{1}
\]

## Core theorem

The map `Phi_H` has the following properties.

1. It is finite of generic degree `n`.  Its inverse equation is
   \[
   H(W)-sW+t=0,
   \qquad
   \gamma=\frac{s-H'(W)}c.
   \tag{2}
   \]
2. Its Jacobian is
   \[
   \det\frac{\partial(s,t)}{\partial(W,\gamma)}=-c^2\gamma.
   \tag{3}
   \]
   Hence its critical divisor is exactly `gamma=0`, with multiplicity one.
3. The restriction to that divisor is
   \[
   \nu_H:\mathbb A^1_W\longrightarrow\mathbb A^2_{s,t},
   \qquad
   W\longmapsto\bigl(H'(W),WH'(W)-H(W)\bigr),
   \tag{4}
   \]
   the normalization of the reduced repeated-root discriminant of (2).
4. If `D_H` denotes that reduced discriminant, then
   \[
   \Omega_{\mathbb A^1_W/D_H}
   \simeq k[W]/(H''(W))\,dW,
   \qquad
   \operatorname{Fitt}_0\Omega_{\mathbb A^1_W/D_H}=(H''(W)).
   \tag{5}
   \]

For finiteness, eliminate `gamma` from (1).  The result is precisely the
monic-up-to-a-unit degree-`n` equation in (2), and `gamma` is then recovered
linearly.  The two rows of the differential of (1) are

\[
 (H''(W),c),\qquad
 (c\gamma+WH''(W),cW),
\]

which proves (3).  On `gamma=0`, equations (1) give (4).  A repeated root of
(2) satisfies exactly these two equations.  The coordinate degrees `n-1`
and `n` in (4) are coprime, so the parameterization is birational; its
finiteness makes it the normalization.  Finally,

\[
 dH'=H''(W)dW,
 \qquad
 d(WH'-H)=WH''(W)dW,
\]

and the transitivity sequence for differentials gives (5), including all
multiplicities.

## Universal-incidence normal form

There is an even simpler normal form hidden in (1).  The triangular source
automorphism

\[
 \theta_H(W,\gamma)=(W,s)=(W,H'(W)+c\gamma)
\]

has determinant `c` and inverse
`\gamma=(s-H'(W))/c`.  In these coordinates,

\[
 \Phi_H=\Psi_H\circ\theta_H,
 \qquad
 \Psi_H(W,s)=(s,Ws-H(W)).
 \tag{6}
\]

Thus the tangent map is, up to a triangular source coordinate, the projection
of the smooth universal marked-root incidence

\[
 \mathcal J_H=V\bigl(H(W)-sW+t\bigr)
 \subset\mathbb A^2_{s,t}\times\mathbb A^1_W.
\]

Indeed `t=sW-H(W)` identifies `J_H` with `A^2_(W,s)`, and its projection to
`A^2_(s,t)` is exactly `Psi_H`.  In particular

\[
 \det D\Psi_H=H'(W)-s,
\]

so the ramification equation, inverse pencil, and discriminant normalization
are all properties of one smooth incidence projection.

The raw weighted incidence is the base change of this universal projection
along

\[
 \mu_0:\mathbb A^3_{A,B,C}\longrightarrow\mathbb A^2_{s,t},
 \qquad
 (A,B,C)\longmapsto(BC,cAC^2):
\]

```text
 I_H = V(H(W)-BCW+cAC^2)  ----->  J_H
              |                          |
              v                          v
          A^3_(A,B,C)  ------ mu_0 ---> A^2_(s,t).
```

Consequently the normalized marked-root model is not an additional inverse
construction: it is the normalization of a ramified base change of the
universal plane incidence, followed by restriction to the
regular-reconstruction open.

## Weighted-suspension theorem

Now assume that `H` is an admissible weighted seed, with the constants
`a_0,b_0` of the weighted construction.  Put

\[
 v=xy,\qquad S=x^2z,\qquad
 \gamma=1+a_0v+b_0S,\qquad
 W=(1+v)\gamma,\qquad C=x\gamma.
\tag{7}
\]

For the weighted map `G_H=(A,B,C)`, define

\[
 \mu(A,B,C)=(s,t,C)=(BC,cAC^2,C),
 \qquad
 \rho(x,y,z)=(W,\gamma,C).
\]

Then the following square commutes:

```text
 A^3_(x,y,z)  -------- G_H -------->  A^3_(A,B,C)
      |                                      |
      | rho                                  | mu
      v                                      v
 A^3_(W,gamma,C) -- Phi_H x id_C -->  A^3_(s,t,C).
```

Indeed, the two nontrivial bottom coordinates are

\[
 s=BC=H'(W)+c\gamma,
 \qquad
 t=cAC^2=W(H'(W)+c\gamma)-H(W).
\tag{8}
\]

The four relevant determinants are

\[
 \det\frac{\partial(x,v,S)}{\partial(x,y,z)}=x^3,
 \qquad
 \det\frac{\partial(W,\gamma,C)}{\partial(x,v,S)}=b_0\gamma^2,
\]

\[
 \det D(\Phi_H\times\mathrm{id}_C)=-c^2\gamma,
 \qquad
 \det D\mu=-cC^3.
\]

Taking determinants in the square and using `C=x gamma` gives

\[
 (-cC^3)\det DG_H
 =(-c^2\gamma)(b_0x^3\gamma^2)
 =-b_0c^2C^3.
\]

Cancellation on the dense open `C!=0`, followed by polynomial continuation,
yields

\[
 \boxed{\det DG_H=b_0c.}
\tag{9}
\]

Thus `G_H` is a **weighted suspension** of the non-etale plane tangent map.
The plane core loses one factor `gamma` along its ramification divisor.  The
weighted third coordinate supplies the complementary factors in the two
vertical maps, so the three-dimensional Jacobian is constant.  This is not a
direct-product suspension: the ramified coordinate changes in the square are
essential.

### Foundational Poisson-square chart

For the foundational seed `H(W)=W^2(1-W)`, the invariant-plane pair
`P=C^2A` and `Q=CB` is a scaled form of the same tangent core:

\[
 Q/4=H'(W)+\gamma,
 \qquad
 P/4=W(H'(W)+\gamma)-H(W).
\]

In oriented invariant coordinates it satisfies
`[P/2,Q]=C^2`.  Thus the inverse pencil, the weighted suspension, and the
plane Poisson-square equation are three coordinate presentations of one
mechanism.  The exact coordinate change and its three-layer
weighted-Wronskian reduction are proved in the
[weighted tangent-suspension bridge](../extended-geometry/WEIGHTED_TANGENT_SUSPENSION.md).

## Reconstruction as the same mechanism

For a target point put `s=BC` and `t=cAC^2`.  Equation (2) is exactly the
weighted inverse pencil

\[
 H(W)-BCW+cAC^2=0.
\]

Moreover its `W`-derivative is

\[
 E_W=H'(W)-BC=-c\gamma.
\]

Consequently

\[
 \gamma=-\frac{E_W}{c},
 \qquad
 x=\frac C\gamma=-\frac{cC}{E_W}.
\]

Away from `C=0`, simple roots are exactly the regularly reconstructible
points and repeated roots are exactly the reconstruction poles.  On the
critical divisor `gamma=0`, the same formula becomes the discriminant
normalization (4), and differentiating that restriction produces the Fitting
divisor (5).  These are therefore different faces of one plane-core
factorization, rather than separate features of the weighted construction.

## A second suspension type: cancellation maps

The cancellation construction has a parallel controlled-boundary core.  Put

\[
 D(s,P,Q)=1-s(Q-Ps)^m,
\qquad
 R(s,P,Q)=C\int_0^s D(t,P,Q)^r\,dt.
\]

For each fixed `P`, the polynomial plane map

\[
 \chi_{m,r;P}:\mathbb A^2_{s,Q}\longrightarrow\mathbb A^2_{Q,R},
 \qquad (s,Q)\longmapsto(Q,R(s,P,Q))
\]

has

\[
 \det D\chi_{m,r;P}=-C D(s,P,Q)^r.                       \tag{10}
\]

The full three-dimensional map is the family of these plane maps over `P`.
Its source chart satisfies

\[
 \det\frac{\partial(s,P,Q)}{\partial(x,y,z)}
 =-A^r=-D^{-r},
\]

whereas

\[
 \det\frac{\partial(P,Q,R)}{\partial(s,P,Q)}
 =C D^r.
\]

Their product is `-C`.  Thus cancellation maps are a second, birational
suspension type: a core loses the boundary power `D^r`, and a rational source
chart contributes its inverse.  The finite cancellation condition on `h`
is exactly what makes the resulting coordinate `R` polynomial back in
`(x,y,z)`.

This differs essentially from the weighted construction.  The weighted
square uses polynomial vertical maps and a simple core ramification divisor;
the cancellation chart is birational and cancels an arbitrary boundary
power.  A broader classification problem is now concrete: classify
controlled-boundary plane maps and their polynomial or birational
constant-Jacobian suspensions.

The full polynomiality and reconstruction theorem for this second type is in
the [cancellation construction](../cancellation/CONSTRUCTION.md).
