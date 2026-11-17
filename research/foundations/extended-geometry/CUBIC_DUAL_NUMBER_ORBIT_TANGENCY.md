# The cubic dual number and polynomial orbit tangency

The dual-number tangent in the normalized sixteen-monomial coefficient
scheme is not a tangent to the three-dimensional normalized `PGL_2` action.
It is nevertheless zero in the tangent quotient by the full polynomial
left--right automorphism functor. These are different statements.

## 1. Normalized binary `sl_2` action

Write

\[
L=aT+bS,\qquad Q=cT^2+dTS+eS^2
\]

and use `E=T partial_S`, `F=S partial_T`, and
`H=T partial_T-S partial_S`. Their raw coefficient fields are

\[
\begin{array}{c|ccccc}
 &\dot a&\dot b&\dot c&\dot d&\dot e\\ \hline
E&b&0&d&2e&0\\
F&0&a&0&2c&d\\
H&a&-b&2c&0&-2e.
\end{array}
\]

Put `m=ad+bc` and `R=a^2e-abd+b^2c`. The raw `SL_2` action fixes `R`.
If `kappa=dot(m)`, the unique first-order factor rescaling that restores
`m=R=1` is

\[
\boxed{\dot L=\rho(L)+\kappa L,\qquad
       \dot Q=\rho(Q)-2\kappa Q.}
\]

On the normalized slice,

\[
\kappa_E=2(ae+bd),\qquad \kappa_F=3ac,\qquad \kappa_H=1.
\]

## 2. Pullback through the global chart

For the chart `Phi(a,y,z)=(a,b,c,d,e)`, the corrected fields pull back to

\[
\begin{aligned}
V_E(a)={}&1+2ay-6a^2z+12a^2y^2-12a^3yz+9a^3y^3-6a^4y^2z,\\
V_E(y)={}&-6z+11y^2-12ayz+9ay^3-6a^2y^2z,\\
V_E(z)={}&\tfrac12(48a^3y^2z^2-126a^2y^3z+96a^2yz^2
 +81ay^4-198ay^2z+48az^2+99y^3-56yz),\\[2mm]
V_F(a)={}&-\tfrac32a^2(-2a^2z+3ay-2),\\
V_F(y)={}&\tfrac12(6a^2z-9ay+8),\\
V_F(z)={}&-\tfrac34(16a^3z^2-42a^2yz+27ay^2+18az-21y),\\[2mm]
V_H={}&2a\partial_a-2y\partial_y-4z\partial_z.
\end{aligned}
\]

The last field is twice the already visible residual-torus field.

## 3. Corresponding target cubic action

Write

\[
P=p_0T^3+T^2S+p_2TS^2+p_3S^3.
\]

Normalization acts on the product by `dot(P)=rho(P)-kappa P`. Hence

\[
\begin{aligned}
W_E&=(1-2p_0p_2,\ 3p_3-2p_2^2,\ -2p_2p_3),\\
W_F&=(-3p_0^2,\ 2-3p_0p_2,\ p_2-3p_0p_3),\\
W_H&=(2p_0,\ -2p_2,\ -4p_3).
\end{aligned}
\]

Exact comparison gives

\[
\boxed{DG\,V_i=W_i\circ G}\qquad(i=E,F,H).
\]

In the announced target coordinates `(u,v,w)=(p3,2p2,2p0)`, these become

\[
\begin{aligned}
\widetilde W_E&=(-uv,\ 6u-v^2,\ 2-vw),\\
\widetilde W_F&=((v-3wu)/2,\ 4-3wv/2,\ -3w^2/2),\\
\widetilde W_H&=(-4u,\ -2v,\ 2w).
\end{aligned}
\]

After the source conjugacy `(a,y,z)=(x,y,-z_original/2)`, the same identity
holds for the announced map `F`.

## 4. The dual-number tangent

Let `H_def` be the explicit deformation from the weighted coefficient
scheme. Coefficient comparison proves that `H_def` is not a constant linear
combination of `DF V_E`, `DF V_F`, and `DF V_H`. Thus the normalized
`PGL_2` orbit does not specifically account for the dual number.

For the full polynomial automorphism functor, however, there is a decisive
general observation. Since

\[
\det DF=-2,
\]

the matrix `(DF)^(-1)=-adj(DF)/2` has polynomial entries. Therefore

\[
\boxed{V_{\rm def}=(DF)^{-1}H_{\rm def}\text{ is polynomial},\qquad
H_{\rm def}=DF\,V_{\rm def}.}
\]

Over the dual numbers, `id+epsilon V_def` is a polynomial automorphism with
inverse `id-epsilon V_def`. The checker also verifies `div(V_def)=0`, as
forced by preservation of the fixed Jacobian determinant to first order.
Thus the class of `H_def` in the tangent quotient by polynomial left--right
equivalence is zero, already using the source side and `W=0`.

This does not say that `V_def` integrates to an algebraic one-parameter
subgroup over a reduced base. In fact, it does not.

The vector field commutes with the source weight Euler field

\[
E_w=x\partial_x-y\partial_y-2z\partial_z
\]

and therefore induces a derivation `bar(V)` on the invariant plane
`k[u,v]=k[xy,x^2z]`. Give `u,v` ordinary degree one. The degree-seven
homogeneous part of this induced field is

\[
\boxed{
\bar V_{\rm top}
=\frac32u^4(3u+v)^2
 \left(-u\partial_u+(6u+v)\partial_v\right).}
\]

Put `f=(3/2)u^4(3u+v)^2` and
`K=-u partial_u+(6u+v) partial_v`. Then `K(f)=-2f`, and hence

\[
(fK)^n(u)=(-1)^n(1\cdot3\cdot5\cdots(2n-1))f^n u\ne0.
\]

Thus the top derivation is not locally nilpotent. A locally nilpotent
filtered derivation has locally nilpotent top homogeneous part, so `V_def`
is not locally nilpotent. The same formula gives unbounded degrees and shows
that it is not locally finite either. Consequently `V_def` generates neither
an algebraic additive nor multiplicative one-parameter subgroup.

The absence of a one-parameter subgroup does not prevent an ordinary reduced
curve in the automorphism group. Such a curve exists here by the following
general characteristic-zero lemma.

### Divergence-free shear lemma

Every divergence-free polynomial vector field on `A^3` is a finite sum of
locally nilpotent derivations. Consequently it is the tangent at the identity
to a polynomial `A^1`-curve in the special polynomial automorphism group.

Indeed, for `V=(P,Q,R)` with `P_x+Q_y+R_z=0`, set

\[
A_2=\int R\,dx,\qquad
A_3=-\int Q\,dx+\int P(0,y,z)\,dy.
\]

Then

\[
V=(\partial_yA_3-\partial_zA_2)\partial_x
  -(\partial_xA_3)\partial_y
  +(\partial_xA_2)\partial_z.
\]

Decompose `A_3`, degree by degree in `(x,y)` over `k[z]`, into powers
`c(z)(x+lambda*y)^d`; this is ordinary binary polarization. Its Hamiltonian
field is

\[
d c(z)(x+\lambda y)^{d-1}
(\lambda\partial_x-\partial_y),
\]

which is locally nilpotent. Similarly, polarizing `A_2` in `(x,z)` over
`k[y]` gives locally nilpotent fields proportional to
`-lambda partial_x+partial_z`. This proves the lemma.

For `V_def`, the exact checker performs this construction rather than merely
invoking it. The two potentials polarize into `46+87=133` locally nilpotent
shears. If these fields are `D_1,...,D_133`, then

\[
\alpha_t=\exp(tD_1)\circ\cdots\circ\exp(tD_{133})
\]

is a polynomial automorphism for every `t`, satisfies `alpha_0=id`, and has
`alpha'_0=V_def`. Therefore

\[
\frac d{dt}\bigg|_{t=0}F\circ\alpha_t=H_{\rm def}.
\]

The quadratic obstruction inside the rigid weighted-support slice records
that this orbit curve must leave that slice at second order. The dual number
is therefore precisely an orbit-tangency artifact: it is not a genuine
infinitesimal moduli class, even when “geometric orbit” is required to mean a
reduced algebraic curve rather than only a dual-number-valued point.

The all-order version is the
[formal orbit-triviality theorem](FORMAL_ORBIT_TRIVIALITY.md): modulo
unrestricted polynomial source automorphisms, every deformation over a local
Artin base is trivial.  Consequently further deformation searches must impose
a complexity filtration, address algebraization over a reduced base, or use
global normalization-boundary invariants; increasing the infinitesimal order
inside this support ansatz cannot produce global moduli.

The exact certificate is
[`verify_cubic_orbit_tangency.py`](../scripts/verify_cubic_orbit_tangency.py).
