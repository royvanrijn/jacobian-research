# The dilation response to scalar correlation at a three-cell maximizer

Let \(A\) be an actual balanced planar three-cell Gaussian noise-stability
maximizer at \(0<r<1\). The mixed variation with respect to shape and
scalar correlation can vanish on every balanced shape direction only
if the dilation speed is the unique minimizer of a positive quadratic
form among speeds with its radial mass flux. The minimizer and the
mixed score constants are determined by a resolvent on the two-dimensional
mass space. A second expression quantifies the shape response in terms
of the part of the dilation that is orthogonal to translations.

## The mixed variation and the mass flux

Use the outward normal of label \(i\) on each interface \(e=(i,j)\).
Write \(d_e=e_i-e_j\), \(U_i=T_r1_{A_i}\), \(a_e=|\nabla(U_i-U_j)|\),
and \(d\nu=\phi_2\,d\mathcal H^1\) on the interfaces. Work in
\(\mathcal H=L^2(a\,d\nu)\) and write \(K_a=a^{-1}K\) for the
interaction operator in this space. Set \(Q_0=I-K_a\) as an operator.
The flux is \(Cf=\sum_e d_e\int_{\Sigma_e} f\,d\nu\), taking values
in the centered label plane \(E=\mathbf1^\perp\). It is bounded by
(IB1). Compact regular speeds of zero flux are dense in \(\ker C\)
by compact approximation followed by fixed compact flux corrections.
[Global reciprocal-gradient and volume-flux bounds, IB1](THREE_CELL_GLOBAL_INTERFACE_BOUNDS.md#a-global-matrix-bound-on-volume-flux)

The two-dimensional translation subspace \(V=\{v\cdot N:v\in\mathbb
R^2\}\) is the \(1/r\) eigenspace of \(K_a\). On \(Z=V^\perp\) the
operator \(K_a\) satisfies \(0\preceq K_a|_Z\preceq I\), and \(Q_0\)
is nonnegative on \(\ker C\). Write \(C_V=C|_V\) and \(C_Z=C|_Z\).
The map \(C_V:V\to E\) is invertible by full centroid rank, since
\(C(v\cdot N)_i=-v\cdot\int_{A_i}x\,d\gamma_2(x)\).
[Negative index and translation eigenspace, G17](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#negative-index-and-the-higher-hermite-blocks-at-a-maximizer)
[Active dimension and full centroid rank](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-active-dimension-is-impossible)

Let \(b=x\cdot N\), \(d=Cb\), and
\(g_e=\partial_s T_s(1_{A_i}-1_{A_j})|_{s=r}\) on \(\Sigma_e\).
The dilation identity yields
\[
 B b=-r(1-r^2)\frac{g}{a},\qquad
 B=I-r^2K_a,\qquad
 d_i=\int_{A_i}(2-|x|^2)\,d\gamma_2.
                                                               \tag{MR1}
\]
The weighted polynomial integrability gives \(b\in\mathcal H\).
Consequently \(g/a\in\mathcal H\), and the mixed pairing
\(\int gf\,d\nu\) extends to all of \(\mathcal H\). The operator
\(B\) satisfies \(B|_V=(1-r)I\) and \((1-r^2)I\preceq B|_Z\preceq I\).
[Dilation identity and integrability, G6](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#translation-equations-on-an-unconstrained-space-of-speeds)
[Radial mass flux, G21](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#balanced-dilation-at-an-endpoint-tie)

## Resolvent on the mass plane

The pairing with \(g\) vanishes on every compact regular balanced
speed exactly when \(g/a=C^*\lambda\) for a unique \(\lambda\in E\).
Indeed these speeds are dense in \(\ker C\), and \(C\) is onto by
invertibility of \(C_V\). Equivalently, \(g_e=\lambda_i-\lambda_j\)
on each regular interface, since there \(C^*\lambda\) has value
\((\lambda_i-\lambda_j)/a_e\).

Define the positive definite operator \(F_r=CB^{-1}C^*\) on \(E\).
Under this vanishing mixed pairing,
\[
 \begin{split}
 b&=B^{-1}C^*F_r^{-1}d,\\
 \lambda&=-\frac{F_r^{-1}d}{r(1-r^2)},\\
 (1-r)^{-1}C_VC_V^*+C_ZC_Z^*
    &\preceq F_r\preceq
        (1-r)^{-1}C_VC_V^*+(1-r^2)^{-1}C_ZC_Z^* .
 \end{split}                                                   \tag{MR2}
\]
The first two identities follow by applying \(C\) to (MR1). The
operator bounds follow by inverting the bounds for the two blocks
of \(B\). Positive definiteness follows from \(C\) being onto and
\(B^{-1}\) being positive definite.

In particular \(b\) minimizes \(\langle f,Bf\rangle_a\) uniquely
over \(Cf=d\). On that affine space \(f=b+h\) has \(h\in\ker C\).
The cross term \(\langle Bb,h\rangle_a\) vanishes and the quadratic
increment \(\langle h,Bh\rangle_a\) is positive for \(h\ne0\).
Conversely this minimizing property implies vanishing of that
cross term and thus of the balanced mixed pairing.

For a nonsimplex maximizer, the vanishing pairing requires \(d\ne0\)
and \(\lambda\cdot d<0\). If \(d=0\) the minimizer is zero, so
\(b=0\). The zero distributional radial derivative of each indicator
then makes the partition conical about the origin, and conical
rigidity gives the simplex. For \(d\ne0\), the sign follows from
\(\lambda\cdot d=-\langle d,F_r^{-1}d\rangle/[r(1-r^2)]\).
[Conical rigidity](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#cones-centered-at-an-arbitrary-point)

## The residual on translation-orthogonal speeds

Define \(c_r=(1-r)/r\), \(R=C_V^{-1}C_Z:Z\to V\),
\(w=C_V^{-1}d\), and \(z=P_Zb\). Then \(b=w-Rz+z\), and every
balanced speed is \(h=y-Ry\) for a unique \(y\in Z\).
On \(V\), \(Q_0=-c_r I\). Thus
\[
 L=Q_0|_Z-c_r R^*R,\qquad
 0\preceq L\preceq I,\qquad
 \langle h,Q_0h\rangle_a=\langle y,Ly\rangle_a.
\]
Here the lower bound follows from balanced stability and the upper
bound from positivity of \(K_a\). Put
\(S=(1+r)(I+R^*R)+(r/c_r)L\), so \(S\succeq (1+r)I\).
Using \(B/r=rQ_0+(r^{-1}-r)I\) in (MR1) gives
\[
 -r(1+r)\int gh\,d\nu
     = \langle S z-R^*w,y\rangle_a,\qquad h=y-Ry .                \tag{MR3}
\]
In fact the translation block of \(B/r\) is \(c_rI\) and contributes
\(-c_r\langle w-Rz,Ry\rangle_a\). The other block contributes
\(\langle [r Q_0|_Z+(r^{-1}-r)I]z,y\rangle_a\).
Substituting \(Q_0|_Z=L+c_rR^*R\) gives (MR3).

If the balanced mixed pairing vanishes, then \(z=S^{-1}R^*w\) and
\(\|R^*w\|_a\ge(1+r)\|z\|_a\). For a nonsimplex maximizer
\(z\ne0\), by rigidity of cones about any center. Thus even
a nonzero radial flux cannot produce vanishing mixed pairing
when its balancing translation \(w\) is orthogonal to the flux
directions \(R(Z)\). The positive lower bound on \(\|z\|_a\)
over nonsimplex maximizers on compact correlation intervals
makes this a uniform requirement on \(\|R^*w\|_a\).
[Uniform nonconical dilation residual, G24--G25](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#a-uniform-gap-from-conical-geometry)

There is a curvature bound involving the same coupling. Suppose
\(A\) ties the simplex \(T\) at \(r\), distinct simplex-optimal
correlations accumulate at \(r\), and the fixed deficit
\(D(s)=S_s(T)-S_s(A)\) has \(D'(r)=0\). The joint shape-correlation
Taylor inequality (G10) extends to all balanced weighted speeds
by density and (MR1). With \(h=y-Ry\) in that inequality, (MR3)
and \(L\preceq I\) give
\[
 \|S z-R^*w\|_a \le r(1+r)\sqrt{D''(r)/2},\qquad
 D''(r)\ge \frac{2}{r^2}
        \left(\|z\|_a-\frac{\|R^*w\|_a}{1+r}\right)_+^2 .
                                                               \tag{MR4}
\]
The second inequality uses the triangle inequality and
\(\|S z\|_a\ge(1+r)\|z\|_a\).
[Mixed contact inequality, G10](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#a-first-crossing-obstruction-at-the-level-of-hermite-weights)

## An abstract example with a vanishing balanced pairing

Fix \(0<r<1\) and choose \(0<\kappa<1\) and
\(0<\alpha<\sqrt{(1-\kappa)/c_r}\).
On \(\mathbb R^2\oplus\mathbb R^2\oplus\mathbb R\) with the
Euclidean inner product, let \(V\) be the first block,
\(K_a=\operatorname{diag}(r^{-1}I,\kappa I,1)\), and
\(C(v,y,t)=v+\alpha y\). The translation eigenvalue is \(1/r\),
and \(0\preceq K_a|_Z\preceq I\). On balanced speeds
\((-\alpha y,y,t)\) one has
\[
 \langle (-\alpha y,y,t),Q_0(-\alpha y,y,t)\rangle
          =(1-\kappa-c_r\alpha^2)|y|^2 .
\]
The sole balanced zero direction is the last block, which
can represent a rotation in this linear model. For any
nonzero \(\lambda\in\mathbb R^2\), define
\[
       \frac{g}{a}=C^*\lambda,\qquad
       b=-r(1-r^2)B^{-1}C^*\lambda .
\]
This obeys (MR1) with \(d=Cb\), and the pairing with \(g\) vanishes
on balanced speeds. The second block of \(b\) is nonzero, since
it equals \(-r(1-r^2)\alpha\lambda/(1-r^2\kappa)\).
Thus the operator inequalities and flux identities used above
permit a nonzero dilation component orthogonal to translations
when the balanced mixed pairing vanishes. In the Gaussian problem
\(b\) and \(C\) are determined by the same interfaces, a geometric
requirement outside the finite-dimensional example.
