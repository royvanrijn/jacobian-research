# Weighted tangent suspension and the plane Poisson-square equation

This note identifies the plane equation hidden in the foundational weighted
Keller map and compares it, with orientations and supports made explicit, to
the Laurent/Newton cascade used in the plane `(72,108)` exclusion.

Work over a characteristic-zero field.  Put

\[
 u=xy,\qquad v=x^2z,
\]

and consider the foundational weight pattern

\[
 F=(x^{-2}A(u,v),x^{-1}B(u,v),xC(u,v)).             \tag{1}
\]

The apparent negative powers cancel for the polynomial supports considered
below.  Assume first that `C` is a nonconstant affine-linear polynomial.

## 1. The general plane-core formula

Choose affine coordinates `(s,t)` on the invariant plane such that

\[
 C=s+2,
 \qquad
 \delta=\det\frac{\partial(s,t)}{\partial(u,v)}=-1. \tag{2}
\]

The second condition fixes the orientation and scale left unspecified by the
phrase “`t` runs along the `C`-level sets.”  For an arbitrary such coordinate
with determinant `delta`, the formula below acquires the factor `-delta`.

Define

\[
 P=C^2A,\qquad Q=CB.                                \tag{3}
\]

Then

\[
 \boxed{
 \det JF=C^{-2}(P_tQ_s-P_sQ_t).}                   \tag{4}
\]

Although (4) is written in the localization at `C`, after multiplication by
`C^2` it is a polynomial identity, so it holds across `C=0` as well.

### Proof

The invariant-coordinate determinant lemma gives

\[
 \det JF=
 \det\begin{pmatrix}
 -2A&A_u&A_v\\
 -B&B_u&B_v\\
 C&C_u&C_v
 \end{pmatrix}.                                    \tag{5}
\]

Changing the derivative columns from `(u,v)` to `(s,t)` multiplies the
determinant by `delta`.  In `(s,t)` coordinates, substitute

\[
 A=P/C^2,\qquad B=Q/C,\qquad C_s=1,\quad C_t=0.
\]

A direct expansion gives

\[
 \det\begin{pmatrix}
 -2A&A_s&A_t\\
 -B&B_s&B_t\\
 C&1&0
 \end{pmatrix}
 =C^{-2}(P_sQ_t-P_tQ_s).                            \tag{6}
\]

Multiplication by `delta=-1` proves (4).  More generally,

\[
 \det JF=-\delta C^{-2}(P_tQ_s-P_sQ_t).             \tag{7}
\]

Thus the normalization in (2) is essential, not cosmetic.

## 2. Poisson-square reduction

With the oriented bracket

\[
 \{P,Q\}_{t,s}=P_tQ_s-P_sQ_t,
\]

equation (4) says that

\[
 \boxed{\det JF=\lambda
 \quad\Longleftrightarrow\quad
 \{P,Q\}_{t,s}=\lambda C^2.}                       \tag{8}
\]

Equivalently, use plane coordinates `(X,T)=(C,t)` and the standard
orientation

\[
 [R,S]_{X,T}=R_XS_T-R_TS_X.
\]

Then

\[
 \det JF=-C^{-2}[P,Q]_{X,T}.                        \tag{9}
\]

For the foundational determinant `-2`, put `\widehat P=P/2`.  The equation
becomes exactly

\[
 \boxed{[\widehat P,Q]_{X,T}=X^2.}                 \tag{10}
\]

This is a plane polynomial pair with Jacobian equal to a square.  It is not a
plane Keller counterexample: its Jacobian is intentionally nonconstant and
vanishes on `X=0`.

If `C` is constant, coordinates satisfying (2) do not exist.  That separate
boundary is the triangular locus already described in the foundational
coefficient-scheme note; the Poisson-square chart is the nonconstant-`C`
chart.

## 3. The explicit foundational pair

At the foundational point, take

\[
 s=-3u-v,\qquad t=-u.                               \tag{11}
\]

Then `det partial(s,t)/partial(u,v)=-1`, `C=s+2`, and direct substitution in
the invariant polynomials gives

\[
 \begin{aligned}
 A&=(t-1)\bigl(s(t-1)^2+t(2t-3)\bigr),\\
 B&=-3s(t-1)^2-2t(3t-4),\\
 C&=s+2.
 \end{aligned}                                      \tag{12}
\]

Consequently

\[
 P=(s+2)^2A,\qquad Q=(s+2)B                         \tag{13}
\]

satisfy the compact identity

\[
 \boxed{P_tQ_s-P_sQ_t=-2(s+2)^2.}                  \tag{14}
\]

Equations (4) and (14) recover `det JF=-2` without expanding a
three-variable Jacobian.

## 4. Relation with the tangent pencil

The same pair is the plane tangent map underlying the weighted marked-root
construction.  Put

\[
 \gamma=\frac{s+2}{2},\qquad
 W=(1-t)\gamma,\qquad H(W)=W^2(1-W).                \tag{15}
\]

Substitution in (13) gives

\[
 \frac Q4=H'(W)+\gamma,
 \qquad
 \frac P4=W\bigl(H'(W)+\gamma\bigr)-H(W).          \tag{16}
\]

Thus

\[
 (W,\gamma)\longmapsto(Q/4,P/4)
\]

is exactly the plane tangent map `Phi_H` with `c=1`.  Its inverse equation is

\[
 \boxed{4H(W)-QW+P=0,}                              \tag{17}
\]

the scaled form of `H(W)-S W+T=0`.  The source change in (15) has determinant
`(s+2)/4`, while the tangent map has determinant `-gamma`; their product,
together with the two target factors `4`, is precisely (14).  This identifies
the Poisson square, the inverse pencil, and the weighted suspension as one
mechanism.

## 5. The sixteen-monomial scheme as a Wronskian cascade

On the normalized chart `C=2-3u-v`, use (11) and put `X=C`.  For the general
sixteen-monomial ansatz, substitution

\[
 u=-t,\qquad v=3t-X+2
\]

has the band form

\[
 \widehat P=p_3(t)X^3+p_2(t)X^2,
 \qquad
 Q=q_2(t)X^2+q_1(t)X.                              \tag{18}
\]

In the coefficient labels of
[the foundational scheme](FOUNDATIONAL_WEIGHTED_COEFFICIENT_SCHEME.md),

\[
\begin{aligned}
p_3={}&\tfrac12(A_{11}t-A_{21}t^2+A_{31}t^3-1),\\
p_2={}&\tfrac12(A_{40}-3A_{31})t^4
 +\tfrac12(3A_{21}-A_{30}-2A_{31})t^3\\
&+\tfrac12(-3A_{11}+A_{20}+2A_{21})t^2
 +(\tfrac32-A_{11})t+1,\\
q_2={}&-B_{01}+B_{11}t-B_{21}t^2,\\
q_1={}&(3B_{21}-B_{30})t^3
 +(-3B_{11}+B_{20}+2B_{21})t^2\\
&+(3B_{01}-2B_{11}-1)t+2B_{01}.
\end{aligned}                                       \tag{19}
\]

For

\[
 \mathcal W_{i,j}(f,g)=ifg'-jf'g,
\]

coefficient comparison in (10) is only the three-layer cascade

\[
\begin{aligned}
\mathcal W_{3,2}(p_3,q_2)&=0, &(K4)\\
\mathcal W_{3,1}(p_3,q_1)+\mathcal W_{2,2}(p_2,q_2)&=0, &(K3)\\
\mathcal W_{2,1}(p_2,q_1)&=1. &(K2)
\end{aligned}                                       \tag{20}
\]

The exact checker verifies that the coefficients of (20) generate the same
normalized ideal as the invariant determinant and reduce to the same
triangular presentation, including `(B_01-3)^2`.  At the foundational point,

\[
 p_3=\tfrac12(t-1)^3,\quad
 p_2=\tfrac12(t-2)(t-1),\quad
 q_2=-3(t-1)^2,\quad q_1=-2(2t-3),                 \tag{21}
\]

so `(K4)`, `(K3)`, and `(K2)` become `0,0,1` immediately.

There is also a band-level view of the nilpotent direction.  Substitute the
ten linear relations of the triangular presentation and write
`B_01=3+epsilon`.  The three residual layers are

\[
\begin{aligned}
K4={}&-\frac{\varepsilon^2}{4}(t-2)(t-1)^2,\\
K3={}&-\frac{\varepsilon^2}{24}
 (3t^4-12t^3+33t^2-64t+36),\\
K2-1={}&-\frac{\varepsilon^2}{48}
 (t-2)(t+2)(3t^2-8t+12).
\end{aligned}                                       \tag{22}
\]

Thus the tangent direction cancels every layer to first order, while the
leading homogeneous Wronskian `K4` already forces `epsilon^2=0`.  The lower
layers confirm the same obstruction; they are not where nilpotence first
appears.

This replaces the normalized three-variable viewpoint by four univariate
bands and three weighted-Wronskian layers.  It does not by itself decompose
the constant-`C` corner `C_10=C_01=0`.  The plane-core construction itself
extends to either one-sided chart where exactly one coefficient is nonzero.
After diagonal normalization these are

\[
 (C_{10},C_{01})=(0,-1),\qquad(-3,0).
\]

Exact Gröbner certificates give the unit ideal on both charts.  Hence the
punctured boundary `C_10 C_01=0`, `(C_10,C_01)\ne(0,0)`, contains no Keller
point in this ansatz.  Set-theoretically, all remaining boundary geometry is
concentrated at the constant-`C` corner.  Its exact primary decomposition has
two nonreduced primary components; their radicals are explicit affine
three-spaces of triangular automorphisms, as described in the
[coefficient-scheme note](FOUNDATIONAL_WEIGHTED_COEFFICIENT_SCHEME.md#the-two-reduced-components-at-the-constant-c-corner).
The same note also computes the degree-ten toric closure of the open
diagonal-gauge orbit and its two attachment lines.  Thus the entire reduced
global coefficient scheme, not only its boundary, is explicit.

## 6. Precise comparison with the plane `(72,108)` cascade

The published plane reduction produces Laurent pairs with

\[
 [P,Q]_{x,y}=x^2.
\]

After `t=xy^2` and `z=y^{-1}`, their upper bands begin

\[
 P=A(t)z^2+\cdots,\qquad Q=D(t)z^3+\cdots,
\]

and the coefficient of `z^4` is

\[
 \mathcal W_{2,3}(A,D)=2AD'-3A'D=t^2.              \tag{23}
\]

Equations (20) and (23) use exactly the same band-bracket rule

\[
 [f(t)z^i,g(t)z^j]
 =z^{i+j-1}\mathcal W_{i,j}(f,g).                  \tag{24}
\]

The differences are structural and should not be suppressed:

- the weighted pair is polynomial in `(X,t)`, while the plane exclusion uses
  a Laurent ring and Newton-polygon support restrictions;
- the weighted bands have heights `(3,2)` and the forced square occurs at
  layer `X^2`, after two homogeneous compatibility layers;
- the `(72,108)` upper bands have heights `(2,3)` and `t^2z^4` occurs already
  in the leading layer;
- the weighted pair is a completely understood ramified plane model, not a
  candidate plane Keller map.

The genuine bridge is therefore exact at the level of the Poisson-square PDE
and its weighted-Wronskian coefficient cascade.  Newton support and the
location of the forced layer remain programme-specific.

The exact certificate is
[`verify_weighted_tangent_suspension.py`](../scripts/verify_weighted_tangent_suspension.py).
