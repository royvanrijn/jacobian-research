# Kernel-line formula for a cyclic logarithmic cokernel

> **Status.**  For an injective map of rank-two bundles whose cokernel is an
> invertible module on its Cartier determinant divisor `D`, this note proves
> that the cokernel line is `K tensor O_D(D)`, where `K` is the kernel line of
> the restricted bundle map on `D`.  Consequently
> `ch_2(coker)=deg_D(K)+D^2/2`.  For the F2 extraction-root component this
> becomes `deg(K_root)+27`.  The subsequent contracted-packet theorem proves
> `deg(K_root)=-e_root` for a nonnegative logarithmic Gauss degree, so the F2
> contribution is `27-e_root`.  The subsequent tangential-coordinate theorem
> computes `e_root=0`, giving the exact cyclic root contribution `27`.

The determinant-line and Chern-character identities are replayed by
[`verify_log_cyclic_cokernel_twist.py`](../scripts/verify_log_cyclic_cokernel_twist.py).

## 1. Setup

Let `X` be a smooth surface and let

\[
 \theta:E\longrightarrow F                              \tag{1.1}
\]

be an injective map of rank-two vector bundles.  Suppose its determinant
vanishes on an effective Cartier divisor `D`, `Fitt_1(coker(theta))=O_X`, and

\[
 \operatorname{coker}(\theta)=i_*L,qquad i:D\hookrightarrow X, \tag{1.2}
\]

where `L` is invertible on `D`.  Locally this is the Smith model
`diag(1,s)` with `D=(s)`.

Restrict (1.1) to `D`.  It has rank one.  Write

\[
 K=\ker(\theta|_D),\qquad I=\operatorname{im}(\theta|_D). \tag{1.3}
\]

The Tor term in the restricted cokernel sequence identifies `K` with
`L tensor O_D(-D)`; the determinant calculation below gives the same result.

## 2. Determinant-line identity

On `D`, the rank-one kernel, image, and cokernel give

\[
 \det(E)|_D=K\otimes I,qquad
 \det(F)|_D=I\otimes L.                           \tag{2.1}
\]

Therefore

\[
 L\otimes K^{-1}
 =\bigl(\det(F)\otimes\det(E)^{-1}\bigr)|_D.     \tag{2.2}
\]

The determinant of `theta` is a section of
`det(F) tensor det(E)^(-1)` with divisor `D`, so this line bundle is `O_X(D)`.
Hence

\[
 \boxed{L\simeq K\otimes O_D(D).}                \tag{2.3}
\]

Thus the cokernel twist is not an independent datum: it is the restricted
kernel line plus the normal bundle of `D`.

## 3. Codimension-two Chern character

Grothendieck--Riemann--Roch for the Cartier embedding gives

\[
 \operatorname{ch}_2(i_*L)
 =\deg_D(L)-\frac12D^2.                           \tag{3.1}
\]

From (2.3),

\[
 \deg_D(L)=\deg_D(K)+D^2.                         \tag{3.2}
\]

Substitution yields

\[
 \boxed{
 \operatorname{ch}_2(\operatorname{coker}\theta)
 =\deg_D(K)+\frac12D^2.}                         \tag{3.3}
\]

The untwisted module `O_D` has `ch_2=-D^2/2`; the actual cyclic cokernel has
the opposite normal-bundle contribution plus the kernel degree.  Confusing
these two sheaves changes the answer by `D^2`.

## 4. Logarithmic interpretation

For a resolved surface map,

\[
 E=f^*\Omega_Y^1(\log D_Y),\qquad
 F=\Omega_X^1(\log D_X).                          \tag{4.1}
\]

The line `K` consists of target logarithmic covectors whose pullback vanishes
to first order along the determinant divisor.  On a contracted component it
is a subline of the fixed target cotangent fiber.  Its degree records the
variation of this annihilated covector direction: equivalently, a
logarithmic Gauss/kernel-line map to `P^1`.

If a single target covector generates `K` without zeros on all of `D`, then
`K` is trivial and (3.3) reduces to `D^2/2`.  Local constancy near one node is
not sufficient; zeros or changes of direction at other special points alter
`deg_D(K)`.

## 5. F2 extraction-root specialization

`LCBBC1` gives

\[
 D_{\rm root}^2=54,qquad \frac12D_{\rm root}^2=27. \tag{5.1}
\]

Therefore the cyclic logarithmic cokernel on the completed root packet must
satisfy

\[
 \boxed{
 \operatorname{ch}_2(\mathcal T^{\log}_{\rm root})
 =\deg(K_{\rm root})+27.}                        \tag{5.2}
\]

In the extraction-root chart, the target tangential covector `dz` generates
the kernel and has no local zero.  What remains unknown is whether this
generator extends without zeros and without a direction change across the
entire determinant component, including every non-toric special point and
future resolution center.

The global twist problem has therefore become the explicit task

\[
 \text{compute }\deg(K_{\rm root}).               \tag{5.3}
\]

The
[`contracted-packet Gauss theorem`](LOG_KERNEL_GAUSS_DEGREE.md) sharpens this
once more.  The root components map to one target point, so

\[
 K_{\rm root}=\gamma_{\rm root}^*O_{P^1}(-1),
 \qquad \deg(K_{\rm root})=-e_{\rm root}\le0,     \tag{5.4}
\]

and (5.2) is `27-e_root<=27`.

Finally, the
[`tangential-coordinate theorem`](LOG_TANGENTIAL_KERNEL_TRIVIALIZATION.md)
uses `f^*z=W^3U^18*unit` itself.  Its logarithmic differential is divisible by
the full determinant ideal, so the fixed covector `dz` trivializes `K_root`.
Hence `e_root=0`, `deg(K_root)=0`, and (5.2) equals `27`.

<!-- status-consumer: LKGD1 8a357250b5005186 -->

<!-- status-consumer: LTKT1 32ac27318f16c20c -->

## 6. Claim boundary and next step

This theorem applies only where `Fitt_1=O_X` and the cokernel is invertible on
its Cartier support.  Noncyclic nodes require their full presentation and may
carry an additional codimension-two correction.

For F2, the next calculation is to transport the kernel generator around the
first exceptional and strict line-at-infinity components, recording its zero
and pole divisor at:

1. the carrier endpoint;
2. the extraction-root node;
3. every non-toric special point on those components; and
4. any point where the outgoing/purity resolution attaches.

Only after this multidegree and the remaining determinant components are
known can (5.2) be inserted into the global logarithmic Chern identity.

## Reproduction

```bash
.venv/bin/python scripts/verify_log_cyclic_cokernel_twist.py
```
