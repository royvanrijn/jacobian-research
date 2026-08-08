# Developable-gradient-image reduction for the final rank-three `[4]` HC4 stratum

## Status

This note starts **after** `HC4RSD70`.  Hence every rank-three `[4]` packet
whose associated quasi-translation has a linear invariant is already closed by
the Wronskian unit obstruction and the fixed-kernel theorem.

The only remaining relative-nilpotent packet has

\[
S=\operatorname{Hess}\psi,\qquad
T=\operatorname{Hess}A,\qquad
N=S^{-1}T,
\]

with

\[
\det(S+sT)=\delta\in K^*,\qquad
\operatorname{rank}T=3
\]

generically and `N` a single nilpotent `[4]` block generically.  Let

\[
F=\nabla A.
\]

Let \(g\) generate the prime relation of the three-dimensional gradient
image

\[
Y=\overline{F(\mathbb A^4)}\subset\mathbb A^4,
\]

and let \(k\) be the primitive polynomial generator of \(\ker T\) obtained
from the associated singular-Hessian quasi-translation.

> **Theorem HC4RSD71 — developable gradient-image reduction.**  In the final
> rank-three `[4]` packet:
>
> 1. the source kernel field is a quasi-translation,
>    \[
>    D:=k\cdot\nabla,\qquad Dk=0,\qquad DF=0;
>    \]
> 2. the projective Gauss map of the gradient-image hypersurface \(Y\) has
>    rank at most two; hence the projective closure of \(Y\) has degenerate
>    Gauss image;
> 3. along every kernel orbit, \(T\) satisfies an exact matrix Riccati
>    equation;
> 4. the HC4 gradient \(H=\nabla\psi\) sends every kernel orbit to a polynomial
>    curve of degree at most three contained in one affine tangent hyperplane
>    of \(Y\).
>
> Consequently a genuinely unresolved HC4 packet must be supported on the
> focal/singular geometry of a classical developable threefold in projective
> four-space.  It is not an arbitrary four-variable moving nilpotent flag.

This does **not** close the final linearly-independent quasi-translation case.
It replaces it by a substantially smaller global projective problem.

## 1. Associated quasi-translation

On the smooth generic locus of \(Y\),

\[
\nabla g(F)=\mu k                                      \tag{1.1}
\]

for a nonzero rational scalar \(\mu\).  After passing to the primitive
associated kernel field, the standard singular-Hessian construction gives

\[
Tk=0,\qquad Dk=0,\qquad DF=Tk=0.                    \tag{1.2}
\]

Thus

\[
\phi_t(x)=x+t k(x)                                  \tag{1.3}
\]

is a polynomial \(\mathbb G_a\)-action and every generic fiber of \(F\) is
one of its affine lines.

The previous cofactor calculation also gives the HC4 null equation

\[
k^{\mathsf T}Sk=0,                                  \tag{1.4}
\]

or equivalently

\[
D^2\psi=0.                                          \tag{1.5}
\]

## 2. The gradient image is developable

Put

\[
B=(\operatorname{Hess}g)(F).
\]

Differentiating the normal field gives, before primitive rescaling,

\[
J(\nabla g(F))=BT.                                  \tag{2.1}
\]

Projectively the scalar rescaling in (1.1) is irrelevant.  Hence the
differential of the projective Gauss map of \(Y\), pulled back by \(F\), is
the map induced by \(BT\) modulo the normal line.

For the associated quasi-translation its Jacobian is nilpotent.  A nilpotent
endomorphism of a four-dimensional vector space has projective image rank at
most two after the normal/kernel line is quotiented out: if its ordinary rank
is at most two this is immediate; if its rank is three, it is a regular
nilpotent block and its one-dimensional kernel is contained in its image, so
quotienting by that line lowers the rank from three to two.

Therefore

\[
\boxed{\operatorname{rank}\gamma_Y\le2},            \tag{2.2}
\]

where \(\gamma_Y\) denotes the projective Gauss map.  The projective closure
of \(Y\) is a developable threefold.

If the Gauss rank is two, a general Gauss fiber is a projective line and its
focal scheme has degree two.  Classical focal theory puts that focal scheme
inside the singular locus of the projective closure.  Thus any non-cylindrical
survivor must use genuine singular/focal geometry; the smooth-affine
cylinder argument cannot simply be bypassed generically.

For a threefold in \(\mathbb P^4\), the classical Gauss-rank-two
classification leaves the five projective types surviving in \(\mathbb P^4\):
`(1,1)`, `(1,2)`, `(2,1)`, `(2,2)`, and `(3,1)`.  Types `(1,3)` and `(2,3)`
do not occur in \(\mathbb P^4\).  Gauss-rank-one cases are the classical
osculating-scroll/cone cases and are to be treated separately.

## 3. Riccati evolution along a collapsed fiber

Differentiate

\[
T k=0.                                               \tag{3.1}
\]

For a coordinate vector \(e_i\),

\[
(\partial_iT)k+T\,\partial_i k=0.                  \tag{3.2}
\]

Third-derivative symmetry of \(A\) identifies the matrix whose \(i\)-th
column is \((\partial_iT)k\) with \(D T\).  Moreover

\[
Jk=BT.                                               \tag{3.3}
\]

Hence

\[
\boxed{DT=-TBT.}                                    \tag{3.4}
\]

Since \(F\) is constant on a \(D\)-orbit, \(B=(\operatorname{Hess}g)(F)\)
is constant on that orbit.  Thus (3.4) is a constant-coefficient matrix
Riccati equation along each collapsed affine line.  Equivalently, wherever
the rational expression is read,

\[
T(\phi_t(x))
=T(x)\bigl(I+tB(x)T(x)\bigr)^{-1}.                 \tag{3.5}
\]

The inverse is in fact polynomial along the orbit because the associated
quasi-translation Jacobian is nilpotent.

## 4. Cubic HC4-gradient curves inside tangent hyperplanes

Write

\[
H=\nabla\psi,\qquad K=Jk.
\]

Quasi-translation invariance gives

\[
k(\phi_t(x))=k(x)                                  \tag{4.1}
\]

and

\[
J\phi_t=I+tK.                                       \tag{4.2}
\]

Because \(D^2\psi=0\),

\[
\psi(\phi_t(x))=\psi(x)+tD\psi(x).                 \tag{4.3}
\]

Differentiate (4.3) with respect to the starting point.  Put

\[
a=D\psi,\qquad p=Sk.
\]

Since

\[
\nabla a=K^{\mathsf T}H+p,                          \tag{4.4}
\]

we obtain

\[
(I+tK)^{\mathsf T}H(\phi_t(x))
=(I+tK^{\mathsf T})H+t p,
\]

and therefore

\[
\boxed{
H(\phi_t(x))
=H+t(I+tK^{\mathsf T})^{-1}p.
}                                                    \tag{4.5}
\]

The polynomial automorphisms \(\phi_t\) have Jacobian determinant one, so

\[
\det(I+tK)=1
\]

and \(K\) is nilpotent.  Furthermore

\[
k^{\mathsf T}p=k^{\mathsf T}Sk=0.                 \tag{4.6}
\]

Thus \(p\in k^\perp\).  Since \(Kk=0\), the
three-dimensional hyperplane \(k^\perp\) is invariant under \(K^{\mathsf T}\).
The restriction of the nilpotent map \(K^{\mathsf T}\) to this
three-dimensional space has cube zero.  Hence

\[
(K^{\mathsf T})^3p=0,
\]

and (4.5) becomes the exact cubic formula

\[
\boxed{
H(\phi_t(x))
=H+t p-t^2K^{\mathsf T}p+t^3(K^{\mathsf T})^2p.
}                                                    \tag{4.7}
\]

All derivatives in \(t\) lie in \(k^\perp\).  Equivalently,

\[
k^{\mathsf T}H(\phi_t(x))=k^{\mathsf T}H(x)       \tag{4.8}
\]

for all \(t\).  Therefore each collapsed source fiber is sent by \(H\) to a
cubic-or-lower polynomial curve lying in the **single affine tangent
hyperplane** to \(Y\) with normal \(k\).

## 5. Quotient/Krylov form of the four vanished determinant coefficients

Choose rational local quotient coordinates \(u=(u_1,u_2,u_3)\) for the
kernel foliation and a fiber parameter \(t\).  Then

\[
F=F(u),\qquad x=x(u,t).
\]

Put

\[
B_0=\frac{\partial F}{\partial u},\qquad
H_u=\frac{\partial H}{\partial u},\qquad
H_t=\frac{\partial H}{\partial t}.
\]

The identity

\[
\det(T+sS)=s^4\delta                               \tag{5.1}
\]

becomes

\[
\det[B_0+sH_u,H_t]
=s^3\delta\det\frac{\partial x}{\partial(u,t)}.    \tag{5.2}
\]

At a generic point choose a target frame \([B_0,n]\), and write

\[
H_u=B_0M+n r^{\mathsf T},\qquad
H_t=B_0c+q n.                                      \tag{5.3}
\]

The constant term of (5.2) gives \(q=0\): this is again tangency of the
fiber curve.  Expanding the remaining determinant gives

\[
-r^{\mathsf T}\operatorname{adj}(I+sM)c
=C s^2,
\qquad C\ne0.                                      \tag{5.4}
\]

Consequently

\[
\boxed{
r^{\mathsf T}c=0,\qquad
r^{\mathsf T}Mc=0,\qquad
r^{\mathsf T}M^2c\ne0.
}                                                    \tag{5.5}
\]

In particular

\[
c,\ Mc,\ M^2c
\]

are linearly independent.  Thus the residual three-dimensional tangent
problem is **cyclic**: the final `[4]` block descends to a length-three Krylov
chain on the tangent space of the developable gradient image.

This is the finite-dimensional quotient object to attack next.

## 6. Research consequence

The final unresolved HC4 mechanism must simultaneously satisfy all of the
following:

1. a rank-three gradient map \(F=\nabla A\) with affine-line fibers;
2. a linearly-independent associated four-variable quasi-translation;
3. a developable gradient-image threefold \(Y\subset\mathbb A^4\) of Gauss
   rank at most two;
4. nontrivial focal/singular geometry of its projective closure;
5. a cubic-or-lower tangent-hyperplane curve supplied by \(H=\nabla\psi\) on
   every source fiber;
6. a cyclic three-dimensional quotient chain (5.5).

The most promising continuation is therefore **not** another polynomial
ansatz.  It is to run the five classical \(\mathbb P^4\) Gauss-rank-two focal
types through conditions (4.7) and (5.5), starting with `(2,1)` and `(1,2)`.
The cone/join types should be tested separately for immediate fixed-direction
or cotangent reductions.

## 7. External geometric input

For the projective part we use the classical facts that general fibers of a
degenerate Gauss map are linear and carry focal hypersurfaces contained in
the singular locus, and Piontkowski's classification of Gauss-rank-two
projective varieties.  In the specialization to hypersurfaces in
\(\mathbb P^4\), only five of the seven general types occur.
