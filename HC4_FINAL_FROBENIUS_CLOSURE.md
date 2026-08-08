# Final Frobenius closure of the rank-three `[4]` HC4 stratum

## Status

This note closes the last moving-flag gap left by `HC4RSD71--74`.

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
>
> Consequently the genuinely moving rank-three `[4]` obstruction does not
> exist. Combined with the earlier triangular/fixed-direction analysis, the
> complete relative-nilpotent HC4 branch reduces to the already isolated
> `HC2` / plane-cotangent (`JC2`) endpoints.

The proof is pointwise and uses only one scalar identity.  The missing
Frobenius scalar is exactly one value of the second fundamental form of the
gradient-image hypersurface, and nilpotence of its polynomial normal field
forces that value to vanish.

## 1. Polynomial normal and quasi-translation

Let

\[
F=\nabla A
\]

and let `g` generate the prime relation of the three-dimensional gradient
image

\[
Y=\overline{F(\mathbb A^4)}\subset\mathbb A^4.
\]

On the smooth generic locus put

\[
n=\nabla g(F).
\tag{1.1}
\]

Then `n` spans `ker T`.  The primitive associated kernel field `k` is a
quasi-translation and `n=\mu k` for a rational first integral `mu` of the
kernel flow.  Since `F` is constant on the kernel orbits, so is `mu`; hence

\[
D_n n=0.
\tag{1.2}
\]

Thus `n` itself is a polynomial quasi-translation on the generic locus.  Its
Jacobian

\[
J_n=J(n)
\]

is nilpotent; in particular

\[
\operatorname{tr}J_n=0.
\tag{1.3}
\]

Writing

\[
B=(\operatorname{Hess}g)(F),
\]

differentiation of (1.1) gives the exact identity

\[
J_n=BT.
\tag{1.4}
\]

## 2. An adapted self-adjoint Jordan chain

Fix a generic point.  Since `N` is `S`-self-adjoint and is a single nilpotent
Jordan block, choose an `S`-adapted Jordan chain `e1,e2,e3,e4` such that

\[
Ne_1=0,\qquad Ne_2=e_1,\qquad Ne_3=e_2,\qquad Ne_4=e_3,
\tag{2.1}
\]

and

\[
S=
\begin{pmatrix}
0&0&0&1\\
0&0&1&0\\
0&1&0&0\\
1&0&0&0
\end{pmatrix}.
\tag{2.2}
\]

Scale the chain so that

\[
e_1=n.
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
\tag{2.3}
\]

The tangent space of `Y` is `im T`, i.e. the span of the target vectors
corresponding to `e2,e3,e4`.

Put

\[
\ell=Se_1=e_4^*,\qquad m=Se_2=e_3^*.
\tag{2.4}
\]

The line `ell` is the Gauss line of the developable gradient image, as in
`HC4RSD71--74`.

## 3. The missing scalar is `II_Y(m,m)`

Write the symmetric matrix `B` as

\[
B=(b_{ij}).
\]

Because `ell=e4*` is the radical of the second fundamental form of `Y`,

\[
B(\ell,e_2^*)=B(\ell,e_3^*)=B(\ell,e_4^*)=0.
\]

In coordinates,

\[
b_{24}=b_{34}=b_{44}=0.
\tag{3.1}
\]

A direct multiplication of `BT`, using (2.3), gives

\[
\operatorname{tr}(BT)=2b_{24}+b_{33}.
\tag{3.2}
\]

Thus on the developable locus

\[
\operatorname{tr}J_n=b_{33}=II_Y(m,m).
\tag{3.3}
\]

But `J_n` is nilpotent, so (1.3) gives

\[
\boxed{II_Y(m,m)=0.}
\tag{3.4}
\]

This is exactly the focal/asymptotic identification that remained open after
`HC4RSD74`: the second Krylov direction is forced to be asymptotic.

## 4. The same scalar is the Frobenius obstruction for `ker N^3`

Define the one-form

\[
\lambda=Se_1.
\tag{4.1}
\]

Self-adjointness gives

\[
\ker N^3=(\operatorname{im}N^3)^{\perp_S}
=(Ke_1)^{\perp_S}=\ker\lambda.
\tag{4.2}
\]

Since `S` is a Hessian, its third derivatives are totally symmetric.  Hence,
for `J_n=BT`,

\[
d\lambda
\quad\text{is represented by}\quad
-(SJ_n-(SJ_n)^{\mathsf T}).
\tag{4.3}
\]

Restricting to

\[
\ker\lambda=\langle e_1,e_2,e_3\rangle
\]

in the adapted frame, the three independent entries are

\[
b_{44},\qquad b_{34},\qquad b_{33}-b_{24}.
\tag{4.4}
\]

The Gauss-line radical conditions (3.1) kill the first two and `b24`; the
nilpotent trace identity (3.4) kills `b33`. Therefore

\[
\boxed{\lambda\wedge d\lambda=0.}
\tag{4.5}
\]

Thus `ker N^3` is Frobenius-integrable.

## 5. The complete Jordan flag

The line distribution `ker N` is automatically integrable.

On the Gauss-rank-two locus, `ker N^2` is the pullback of the Gauss-fiber
line together with the source kernel line. Equivalently it is the kernel of
the composite map from the source to the two-dimensional Gauss image.
Therefore

\[
\ker N^2
\]

is Frobenius-integrable.

Section 4 proves the only remaining condition:

\[
\ker N^3
\]

is Frobenius-integrable.

For an operator similar to one nilpotent Jordan block, integrability of the
complete image/kernel flag is equivalent to the existence of local coordinates
in which the operator is strictly upper triangular.  This is the elementary
flag criterion recalled in Bolsinov--Konyaev--Matveev and is compatible with
their dimension-four triangularization theorem.

Hence `N` is locally strictly upper triangularizable on the generic
constant-rank locus. This removes the final genuinely moving `[4]` stratum.

## 6. What was actually needed

The projective focal classification was useful for locating the gap, but is
not needed to close it.  The final argument uses only:

1. the gradient image is developable and `ell=Sn` is its Gauss line;
2. `n=grad g(F)` is a quasi-translation, hence `Jn` is nilpotent;
3. `N` is `S`-self-adjoint;
4. `S` is a Hessian.

The crucial identity is

\[
\boxed{
\operatorname{tr}J_n
=II_Y(m,m)
=\text{the sole remaining Frobenius scalar}
}
\]

after imposing the Gauss-line radical equations.

## 7. Verification

Run

```bash
.venv/bin/python scripts/verify_hc4_final_frobenius_closure.py
```

The checker verifies the canonical self-adjoint Jordan pair, the trace formula,
and the equality of the remaining Frobenius coefficient with the same scalar
`b33` after imposing the Gauss-line radical equations.

## 8. External reference

A. V. Bolsinov, A. Yu. Konyaev, V. S. Matveev,
*On the Jordan-Chevalley decomposition problem for operator fields in small
dimensions and Tempesta-Tondo conjecture*, Journal of Geometry and Physics 218
(2025), 105656; arXiv:2503.10208.  The paper recalls that, for a single
nilpotent Jordan block, triangular coordinates imply integrability of the full
Jordan flag and states the converse; its dimension-four theorem supplies an
intrinsic tensor criterion for the same triangularization problem.
