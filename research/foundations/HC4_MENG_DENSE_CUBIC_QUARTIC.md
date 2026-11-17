# Dense cubic--quartic reduction in the Meng `HC_4` chart

## Status

This note closes the complete mixed-degree potential

\[
 \psi=q_2+h_3+h_4
\]

in the collision-normalized Meng chart.  Here \(q_2\) is any nondegenerate
quadratic form in four variables, while \(h_3,h_4\) are arbitrary
homogeneous cubic and quartic forms.  There is no support restriction.

> **Theorem `HC4CQ1`.**  If
> \(\det\operatorname{Hess}(\psi)\) is a nonzero constant, then
> \(\nabla\psi\) cannot identify a nonzero antipodal pair.

The proof uses the Gordan--Noether classification in dimensions at most
four, the known Hessian conjectures in dimensions two and three, and Moh's
plane Jacobian degree bound.  The accompanying checker verifies the exact
determinant identities between those structural inputs.

This is not a proof of `HC_4`.  Adding a sextic layer, mixing three
homogeneous degrees, or changing the coisotropic embedding lies outside the
theorem.

We may extend the ground field to its algebraic closure throughout.  A
collision and the constant-Hessian identity survive scalar extension, so
exclusion there implies exclusion over the original characteristic-zero
field.

## 1. Homogeneous determinant layers

Put

\[
 H_0=\operatorname{Hess}(q_2),\qquad
 A=\operatorname{Hess}(h_3),\qquad
 B=\operatorname{Hess}(h_4).
\]

Scaling the variables by \(\lambda\) gives

\[
 \det(H_0+\lambda A+\lambda^2B)=\det H_0.             \tag{1.1}
\]

The coefficient of \(\lambda^8\) is \(\det B\), so it vanishes.
The Gordan--Noether theorem implies that \(h_4\) has a nonzero constant
kernel direction.  More precisely, if the generic rank of \(B\) is \(r\),
iterate the theorem after removing one inessential variable.  This gives a
constant kernel \(K\) of dimension \(4-r\), and coordinates in which

\[
 B=\begin{pmatrix}0&0\\0&C_r\end{pmatrix},
 \qquad \det C_r\ne0.
\]

For \(r=3,2,1\), the highest remaining coefficient in (1.1) is respectively

\[
\begin{array}{c|c}
r&\text{coefficient}\\ \hline
3&\lambda^7:\quad \det(C_3)\det(A|_K),\\
2&\lambda^6:\quad \det(C_2)\det(A|_K),\\
1&\lambda^5:\quad \det(C_1)\det(A|_K).
\end{array}                                           \tag{1.2}
\]

For \(r=0\), the coefficient of \(\lambda^4\) is \(\det A\).
Thus the Hessian block of \(h_3\) along \(K\) has zero determinant.  Turning
this into one constant direction requires a short rank analysis:

- If \(r=3\), then \(\dim K=1\), and (1.2) directly says
  \(D_v^2h_3=0\).
- If \(r=2\), the symmetric \(2\times2\) Hessian block has linear-form
  entries \(a,b,c\) satisfying \(ac=b^2\).  Unique factorization forces
  its nonzero rows to be constant-proportional, so it has a constant
  kernel direction.
- If \(r=1\), the Hessian restriction has the form
  \(\operatorname{Hess}(f_3)(u)+z\operatorname{Hess}(f_2)\) on the
  three-dimensional kernel; the omitted \(z^2f_1(u)+z^3f_0\) terms have
  zero \(u\)-Hessian.  Gordan--Noether
  makes \(\operatorname{Hess}(f_3)\) constant-kernel.  If its rank is two,
  one, or zero, the next coefficients of
  \[
    \det\bigl(\operatorname{Hess}(f_3)
       +z\operatorname{Hess}(f_2)\bigr)
  \]
  are respectively
  \[
    \det(P_2)(f_2)_{mm},\qquad
    p_{11}\det(\operatorname{Hess}(f_2)|_{\ker P_1}),\qquad
    \det\operatorname{Hess}(f_2).
  \]
  Their vanishing supplies a constant kernel direction of
  \(\operatorname{Hess}(f_3)\) which is isotropic for
  \(\operatorname{Hess}(f_2)\).
- If \(r=0\), Gordan--Noether applies directly to \(h_3\).

In every case there is a nonzero \(v\in K\) such that

\[
 D_vh_4=0,\qquad D_v^2h_3=0.                          \tag{1.3}
\]

Choose a coordinate \(t\) in this direction and write the other three
coordinates as \(u\).  Then

\[
 h_4=b_4(u),\qquad h_3=t\,a_2(u)+b_3(u).              \tag{1.4}
\]

## 2. The nonisotropic direction

Subtract the common gradient value of the antipodal pair from the
potential.  This does not change its Hessian and makes both points critical.
Write

\[
 \psi=\frac{\kappa}{2}t^2+t\,s(u)+\phi(u).
\]

If \(\kappa\ne0\), the critical equation in \(t\) has the polynomial
solution \(t=-s/\kappa\).  Its Schur complement is

\[
 \bar\psi(u)=\phi(u)-\frac{s(u)^2}{2\kappa},
\]

and

\[
 \det\operatorname{Hess}(\psi)
 =\kappa\det\operatorname{Hess}(\bar\psi).
\]

The two critical points descend to distinct points: if their \(u\)
coordinates agreed, the unique critical value of \(t\) would also agree.
This would give a three-variable constant-Hessian collision, contradicting
\(\mathrm{HC}_3\).

## 3. The isotropic bordered form

It remains to take \(\kappa=0\):

\[
 \psi=t\,s(u)+\phi(u),\qquad
 s=c+\ell(u)+a_2(u).                                  \tag{3.1}
\]

Here \(c\) comes from subtracting the common gradient value,
\(\ell\ne0\) because \(q_2\) is nondegenerate, and
\(C=\operatorname{Hess}(a_2)\) is constant.  Put
\(v=\ell+\nabla a_2\).  The coefficient of \(t^2\) in the bordered Hessian
determinant is

\[
 -v^{\mathsf T}\operatorname{adj}(C)v.                \tag{3.2}
\]

Equation (3.2) first forces \(\det C=0\).  If \(\operatorname{rank}C=2\),
write \(\operatorname{adj}C=\rho\,mm^{\mathsf T}\); then (3.2) gives
\(m^{\mathsf T}\ell=0\).  If the rank is at most one, choose
\(m\in\ker C\cap\ker\ell\).  In either case \(s\) is independent of the
coordinate \(m\).  Write the other two coordinates as \(x,y\).

Let

\[
 R=(\nabla s)^{\mathsf T}
   \operatorname{adj}(\operatorname{Hess}_{x,y}s)\nabla s.
\]

The coefficient of \(t\) factors exactly as

\[
 -\phi_{mm}R.                                         \tag{3.3}
\]

The polynomial ring is a domain, so there are two cases.

## 4. The cotangent-lift case

If \(\phi_{mm}=0\), then

\[
 \phi=m\,g(x,y)+h(x,y)
\]

and

\[
 \psi=t\,s(x,y)+m\,g(x,y)+h(x,y).
\]

Direct block expansion gives

\[
 \det\operatorname{Hess}(\psi)
 =\operatorname{Jac}(s,g)^2.                          \tag{4.1}
\]

Here \(\deg s\le2\) and \(\deg g\le3\).  Moh's degree bound makes the plane
Keller map \((s,g)\) a polynomial automorphism.  Formula (4.1) is a
cotangent lift, so the full four-variable gradient is also a polynomial
automorphism.  It has no collision.

## 5. The binary bordered case

Suppose \(R=0\).  Normalize the nonzero linear form \(\ell\) to \(x\) and
write

\[
 s=c+x+\alpha x^2+\beta xy+\gamma y^2.
\]

The constant coefficient of \(R\) is \(2\gamma\).  After \(\gamma=0\), the
coefficient of \(x\) is \(-2\beta^2\).  Hence

\[
 \beta=\gamma=0,\qquad s=c+x+\alpha x^2.
\]

If \(\alpha\ne0\), then \(\nabla s=0\) at
\(x=-1/(2\alpha)\).  The first row and column of the bordered Hessian vanish
there, contradicting its nonzero constant determinant.  Therefore
\(\alpha=0\) and \(s=c+x\).

Both antipodal critical points satisfy \(s=0\), which forces \(c=x=0\).
Expansion along the \(t,x\) block gives

\[
 \det\operatorname{Hess}(\psi)
 =-\det\operatorname{Hess}_{y,m}(\phi).               \tag{5.1}
\]

At \(x=0\), \(\mathrm{HC}_2\) makes the gradient of
\(\phi(0,y,m)\) injective.  The two antipodal critical points therefore
have \(y=m=0\).  The remaining equation
\(t+\phi_x=0\) fixes \(t\) uniquely.  Since both \(t\) and \(-t\) must
equal that same value, \(t=0\), so the two points coincide.  This is the
final contradiction.

## Reproduction

Run:

```bash
.venv/bin/python \
  scripts/verify_hc4_meng_dense_cubic_quartic_reduction.py
```

The shared leading-layer and bordered-block identities are checked by
`scripts/verify_hc4_meng_dense_rank_three_sextic_reduction.py`, which the
command above replays silently.

The external structural inputs are the
[Gordan--Noether classification](https://arxiv.org/abs/1501.05168),
the known Hessian conjecture in dimensions at most three, and
[Moh's plane degree bound](https://www.math.purdue.edu/~ttm/jacobian.pdf).
