# Smooth-chart obstruction in the final rank-three `[4]` HC4 stratum

## Status

This note continues `HC4RSD71`.  It studies the smooth Gauss-rank-two locus of

\[
Y=\overline{\operatorname{im}\nabla A}\subset \mathbb A^4
\]

for the only remaining relative-nilpotent HC4 packet.

> **Theorem HC4RSD72 — maximal Gauss-kernel motion.**
> On every generic smooth chart of `Y` whose projective Gauss map has rank two,
> the HC4-selected Gauss-kernel line is constant along each collapsed source
> fiber, while its projective differential transverse to that fiber has rank
> exactly two.  In determinant-one frozen normal coordinates its transverse
> derivative is
> \[
> P=\begin{pmatrix}a&b\\c&a\end{pmatrix},
> \qquad \det P=a^2-bc=\delta\ne0.
> \]
> In particular there is no rank-one intermediate motion.
>
> **Corollary HC4RSD73 — affine focal singularities are compulsory.**
> A genuinely unresolved `[4]` packet cannot have an affinely smooth
> Gauss-rank-two gradient image.  Hence its projective focal scheme must meet
> the affine chart: any survivor is necessarily supported on finite
> singular/focal geometry and associated nonproperness.

The theorem is local and algebraic.  The corollary uses Piontkowski's global
classification of affinely smooth developable varieties of Gauss rank two.

## 1. Hamilton--Jacobi / partial Legendre chart

At a generic smooth point choose target coordinates so that

\[
Y:\quad y_4=\mathcal H(p_1,p_2,p_3).
\]

Choose source coordinates `(q,z)` so the principal `3 x 3` minor
`A_{qq}` is invertible.  Since `F=grad A` takes values in `Y`,

\[
A_z=\mathcal H(A_q).
\]

Put `p=A_q` and take the partial Legendre transform in `q`.  There is a local
function `L(p)` with

\[
q=\nabla L(p)-z\nabla\mathcal H(p).
\]

Write

\[
B=\operatorname{Hess}\mathcal H,
\qquad
M=L''-zB.
\]

Then `M` is invertible and the Jacobian of the coordinate map `(p,z)->(q,z)`
is

\[
R=\begin{pmatrix}M&-\nabla\mathcal H\\0&1\end{pmatrix}.
\]

A direct calculation gives

\[
R^TTR=\begin{pmatrix}M&0\\0&0\end{pmatrix}.      \tag{1.1}
\]

## 2. The companion metric

Let `Psi(p,z)=psi(q(p,z),z)`.  The top `[4]` coefficient gives

\[
\Psi_{zz}=0,
\]

because the `z`-direction at fixed `p` is exactly the Hessian-kernel direction
of `T`.  Hence

\[
\Psi=U(p)+zV(p).
\]

If `v` denotes the first three components of `grad psi` in the original
`(q,z)` coordinates, then

\[
v=M^{-1}(\nabla U+z\nabla V).
\]

Set

\[
u=v_z,\qquad W=v_p.
\]

Hessian symmetry gives the exact transformed companion matrix

\[
R^TSR=
\begin{pmatrix}
MW&Mu\\
u^TM&0
\end{pmatrix}.                                      \tag{2.1}
\]

Consequently

\[
\det(R^T(S+sT)R)
=-\det(M)\,u^TM\operatorname{adj}(W+sI)u.
\]

Since `det R=det M` and `det(S+sT)=delta`,

\[
\boxed{
-u^TM\operatorname{adj}(W+sI)u=\delta\det M.
}                                                    \tag{2.2}
\]

The `s^2` and `s` coefficients are

\[
u^TMu=0,
\qquad
u^TMWu=0.                                           \tag{2.3}
\]

## 3. No rotation along a collapsed source fiber

Differentiate the explicit formula for `v`.  Since `M_z=-B`,

\[
u_z=2M^{-1}Bu.                                   \tag{3.1}
\]

Differentiating the first equation of (2.3) gives

\[
0=3u^TBu.
\]

Thus

\[
u^TBu=0,
\qquad
u^TL''u=0.                                          \tag{3.2}
\]

At the chosen point shift `z` so `L''=M(0)` is nonsingular.  When the Gauss
rank of `Y` is two, `B` has rank two.  The nonsingular ternary conic defined by
`L''` and the rank-two conic defined by `B` have no common projective component.
Their intersection is finite.  Hence the rational map `z -> [u]` is constant.
Write `u=rho(z)d`.  Equation (3.1) then gives either `Bd=0`, or

\[
\rho(z)=\rho(0)(1-\mu z)^{-2}
\]

for a nonzero scalar `mu`.  But `v=grad_q psi` evaluated on the affine kernel
orbit is polynomial in `z`; therefore `u=v_z` is polynomial in `z`.  The pole is
impossible, so

\[
\boxed{u_z=0,\qquad Bu=0.}                          \tag{3.3}
\]

Geometrically, `u` is precisely the Gauss-kernel direction of the smooth graph
`Y` and is constant on each Gauss fiber.

In the original four-variable coordinates this says, for the source
quasi-translation kernel `k`,

\[
D_k(Sk)=0.                                           \tag{3.4}
\]

Third-derivative symmetry and `D_kk=0` give

\[
D_k(Sk)=-2(Jk)^T(Sk),
\]

hence

\[
(Jk)^T(Sk)=0.                                        \tag{3.5}
\]

The cubic orbit formula of `HC4RSD71` therefore collapses to the affine formula

\[
\boxed{
\nabla\psi(x+t k(x))=\nabla\psi(x)+tS(x)k(x).
}                                                    \tag{3.6}
\]

## 4. The transverse motion is forced to have rank two

Freeze a determinant-one affine frame at a generic point and normalize

\[
u=e_3,
\qquad
B=\begin{pmatrix}0&1&0\\1&0&0\\0&0&0\end{pmatrix}.
\]

Rescale the representative `u` so its derivative has no component along `u`
at that point.  Differentiating `Bu=0` and using the full symmetry of the third
derivatives of `mathcal H` says `B u_p` is symmetric.  Therefore the most
general projective first jet is

\[
U:=u_p=
\begin{pmatrix}
a&b&0\\
c&a&0\\
0&0&0
\end{pmatrix}.                                      \tag{4.1}
\]

Write

\[
L''=
\begin{pmatrix}
A_{11}&A_{12}&\alpha\\
A_{12}&A_{22}&\beta\\
\alpha&\beta&0
\end{pmatrix},                                      \tag{4.2}
\]

where the bottom-right zero is (3.2), and put

\[
W=W_0+zU.
\]

Symmetry of `MW` gives, denoting the entries of `W_0` by `w_ij`,

\[
w_{13}=-\alpha b-\beta a,
\qquad
w_{23}=-\alpha a-\beta c.                           \tag{4.3}
\]

The coefficient of `s` in (2.2) is

\[
\alpha^2b+2\alpha\beta a+\beta^2c=0.               \tag{4.4}
\]

The coefficient of `z` is

\[
2\alpha\beta\bigl(\delta-(a^2-bc)\bigr)=0.         \tag{4.5}
\]

If `alpha beta != 0`, this already gives

\[
a^2-bc=\delta.
\]

If `beta=0`, nonsingularity of `L''` gives `alpha A_22 != 0`; (4.4) gives
`b=0`, and the constant term of (2.2), after the three remaining symmetry
relations are used, is

\[
\alpha^2A_{22}(\delta-a^2)=0.
\]

Thus again `a^2-bc=delta`.  The case `alpha=0` is symmetric and gives
`c=0`, `a^2=delta`.

Therefore in every case

\[
\boxed{\det U_{\mathrm{proj}}=a^2-bc=\delta\ne0.}   \tag{4.6}
\]

Up to the harmless nonzero square introduced by a non-unit frozen frame, the
statement is invariant: the differential of the projective Gauss-kernel line
has rank exactly two.

## 5. Affinely smooth gradient images are impossible

Piontkowski's theorem says that an affinely smooth, non-conical developable
variety of dimension `n` and Gauss rank two is the union of a **one-dimensional
family** of `(n-1)`-planes.  These planes are the Gauss-fiber cones: all Gauss
fibers inside one such plane have the same point at infinity, hence the same
affine ruling direction.

For `n=3`, the ruling-direction map of an affinely smooth noncone therefore
has image of dimension at most one.  This contradicts (4.6).

If the projective closure is a cone with vertex at infinity, the affine part has
a fixed direction; this is already the linearly-dependent/fixed branch closed
by `HC4RSD65` and `HC4RSD70`.

Consequently a genuinely unresolved packet cannot be affinely smooth:

\[
\boxed{
\operatorname{Sing}(\overline Y)\cap\mathbb A^4\ne\varnothing.
}                                                    \tag{5.1}
\]

More precisely, the degree-two focal scheme on a general Gauss line cannot be
pushed entirely to infinity.  The surviving mechanism necessarily uses finite
focal/singular geometry and nonproperness.

## 6. Research implication

The final `[4]` branch has now lost both smooth alternatives:

1. fixed/linearly-dependent rulings are closed;
2. affinely smooth moving rulings are incompatible with HC4.

Any counterexample must therefore combine

- a rank-three gradient map with affine-line fibers;
- a Gauss-rank-two developable image with finite focal singularities;
- a maximally varying Gauss-line direction map;
- and an HC4 gradient which maps every collapsed source line to a parallel
  affine line by (3.6).

This is exactly the sort of nonproper tangent/direction-sweep geometry that
appears in the new 2026 counterexamples to `JC_n` for `n>=3`; the remaining
question is whether Hessian symmetry prevents the pole-absorption mechanism in
four variables.