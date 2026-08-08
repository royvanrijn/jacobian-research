# Final Frobenius closure of the rank-three `[4]` HC4 stratum

## Status

This note closes the focal/Frobenius gap left by `HC4RSD71--74`.

> **Theorem HC4RSD75 — final `[4]` Frobenius closure.**
> On the generic rank-three locus of the last relative-nilpotent `HC4` packet,
> let
> \[
> S=\operatorname{Hess}\psi,\qquad
> T=\operatorname{Hess}A,\qquad
> N=S^{-1}T,
> \]
> with `N` a single nilpotent Jordan block of length four. Then the complete
> Jordan flag
> \[
> \ker N\subset \ker N^2\subset \ker N^3
> \]
> is Frobenius-integrable. Hence, on the constant-rank generic locus, `N` is
> locally strictly upper triangularizable.

This is a genuine closure of the **moving-flag differential obstruction**.  It
is not, by itself, a proof that the triangularizing coordinates can be chosen
affine or polynomial.  That affine-straightening question is the next global
step before one may identify the generic `[4]` packet with the already closed
fixed-affine-flag / `HC2` / `JC2` endpoints.

The proof uses only one scalar identity. The missing Frobenius scalar is
exactly one value of the second fundamental form of the gradient-image
hypersurface, and nilpotence of its polynomial normal field forces that value
to vanish.

## 1. Polynomial normal and quasi-translation

Let
\[
F=\nabla A
\]
and let `g` generate the prime relation of
\[
Y=\overline{F(\mathbb A^4)}\subset\mathbb A^4.
\]
On the smooth generic locus put
\[
n=\nabla g(F).
\tag{1.1}
\]
Then `n` spans `ker T`.  The primitive associated kernel field `k` is a
quasi-translation and `n=\mu k`, where `mu` is constant along the kernel
orbits because it factors through `F`. Hence `n` is itself a quasi-translation:
\[
D_n n=0.
\tag{1.2}
\]
Therefore its Jacobian `J_n` is nilpotent and
\[
\operatorname{tr}J_n=0.
\tag{1.3}
\]
Writing
\[
B=(\operatorname{Hess}g)(F),
\]
differentiation gives
\[
J_n=BT.
\tag{1.4}
\]

## 2. Adapted self-adjoint Jordan chain

At a generic point choose an `S`-adapted Jordan chain with
\[
Ne_1=0,\quad Ne_2=e_1,\quad Ne_3=e_2,\quad Ne_4=e_3,
\]
scale `e1=n`, and normalize
\[
S=
\begin{pmatrix}
0&0&0&1\\
0&0&1&0\\
0&1&0&0\\
1&0&0&0
\end{pmatrix}.
\tag{2.1}
\]
Then
\[
T=SN=
\begin{pmatrix}
0&0&0&0\\
0&0&0&1\\
0&0&1&0\\
0&1&0&0
\end{pmatrix}.
\tag{2.2}
\]
Put
\[
\ell=Se_1=e_4^*,\qquad m=Se_2=e_3^*.
\tag{2.3}
\]
The line `ell` is the Gauss line of the developable gradient image.

## 3. The trace is the missing asymptotic scalar

Write `B=(b_ij)`.  Since `ell=e4*` is the radical of the second fundamental
form of `Y`,
\[
b_{24}=b_{34}=b_{44}=0.
\tag{3.1}
\]
Direct multiplication gives
\[
\operatorname{tr}(BT)=2b_{24}+b_{33}.
\tag{3.2}
\]
Therefore, on the developable locus,
\[
\operatorname{tr}J_n=b_{33}=II_Y(m,m).
\tag{3.3}
\]
Nilpotence of `J_n` now yields
\[
\boxed{II_Y(m,m)=0.}
\tag{3.4}
\]
Thus the second Krylov direction is automatically asymptotic.

## 4. The same scalar is the Frobenius obstruction

Let
\[
\lambda=Se_1.
\]
Self-adjointness gives
\[
\ker N^3=(Ke_1)^{\perp_S}=\ker\lambda.
\tag{4.1}
\]
Because `S` is a Hessian, third-derivative symmetry gives
\[
d\lambda\sim -(SJ_n-(SJ_n)^{\mathsf T}).
\tag{4.2}
\]
Restricted to
\[
\ker\lambda=\langle e_1,e_2,e_3\rangle,
\]
the three independent coefficients are
\[
b_{44},\qquad b_{34},\qquad b_{33}-b_{24}.
\tag{4.3}
\]
The Gauss-line equations (3.1) kill the first two and `b24`; (3.4) kills
`b33`. Hence
\[
\boxed{\lambda\wedge d\lambda=0.}
\tag{4.4}
\]
and `ker N^3` is Frobenius-integrable.

The line `ker N` is automatically integrable. On the Gauss-rank-two locus,
`ker N^2` is the kernel of the map from the source to the two-dimensional
Gauss image, hence is also Frobenius-integrable. Thus the complete Jordan flag
is integrable.

For an operator similar to one nilpotent Jordan block, integrability of the
complete flag is equivalent to local strict upper triangularizability; see
Bolsinov--Konyaev--Matveev.

## 5. What remains

`HC4RSD75` eliminates the last **differential** modulus in the `[4]` flag.  The
remaining global issue is sharper:

> Does a polynomial Hessian pencil whose regular nilpotent Jordan flag is
> Frobenius-integrable admit a **constant affine** flag, or otherwise reduce
> canonically to the plane-cotangent (`JC2`) packet?

A nonlinear triangularizing coordinate change cannot simply be substituted
into the existing Hessian normal forms, because the Hessian property is tied
to the original flat affine connection.

## 6. Verification

Run

```bash
.venv/bin/python scripts/verify_hc4_final_frobenius_closure.py
```

The checker verifies the canonical self-adjoint Jordan pair, the trace formula,
and the equality of the remaining Frobenius coefficient with the same scalar
`b33` after imposing the Gauss-line radical equations.

## 7. External reference

A. V. Bolsinov, A. Yu. Konyaev, V. S. Matveev,
*On the Jordan-Chevalley decomposition problem for operator fields in small
dimensions and Tempesta-Tondo conjecture*, Journal of Geometry and Physics 218
(2025), 105656; arXiv:2503.10208.
