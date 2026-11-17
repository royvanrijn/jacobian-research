# Tangential-coordinate trivialization of a cyclic logarithmic kernel

> **Status.**  This note proves that a fixed target tangential coordinate
> trivializes the kernel line whenever its pullback divisor contains the full
> cyclic logarithmic determinant divisor.  The statement includes all
> nilpotent thickness: if `f^*z=s*b`, where `s` cuts out `D`, then
> `d(f^*z)` is divisible by `s` in logarithmic differentials.  For the F2
> extraction-root packet, `f^*z=W^3U^18*unit` and
> `D_root=3E+18L`; hence `K_root=O_D*dz`, the Gauss degree is exactly zero,
> and `ch_2(T_root^log)=27`.  This fixes the complete cyclic root contribution
> but does not control noncyclic attachment points or the global Chern ledger.
> The later affine-purity frontier forces a separate new source component and
> raises the global component floors to `28/49`; its Chern contribution remains
> undetermined until the target curve and pullback factorization are known.

<!-- status-consumer: PF2APF1 192055eb737d3140 -->

The logarithmic divisibility, residual determinant, and Chern arithmetic are
replayed by
[`verify_log_tangential_kernel_trivialization.py`](../scripts/verify_log_tangential_kernel_trivialization.py).

## 1. Logarithmic divisibility lemma

Let `X` be smooth with SNC boundary `B`, and let

\[
 D=\sum_i m_iB_i
\]

be an effective boundary-supported Cartier divisor.  In an SNC chart write

\[
 s=\prod_i x_i^{m_i}.                            \tag{1.1}
\]

If a regular function `g` has divisor containing `D`, write `g=sb`.  Then

\[
 dg=s\left(db+b\sum_i m_i\,d\log x_i\right).    \tag{1.2}
\]

Thus

\[
 dg\in s\,\Omega_X^1(\log B).                   \tag{1.3}
\]

This is divisibility by the full nonreduced ideal of `D`, not merely
vanishing after restriction to `D_red`.  It is preserved by every boundary
blowup after replacing `D` by its total transform.

## 2. Kernel trivialization theorem

Let

\[
 \theta_f^{\log}:f^*\Omega_Y^1(\log C)
 \longrightarrow\Omega_X^1(\log B)              \tag{2.1}
\]

have unit `Fitt_1` and cyclic cokernel on `D`.  Suppose the support of `D`
maps into a target chart carrying a regular coordinate `z` such that `dz` is
part of a frame of `Omega_Y^1(log C)`, and suppose

\[
 \operatorname{div}(f^*z)\ge D.                 \tag{2.2}
\]

By (1.3), `theta_f^log(dz)` lies in
`I_D Omega_X^1(log B)`, so `dz` belongs to the kernel of the restricted map.
It is a nowhere-zero member of the fixed target frame.  Since the restricted
kernel is a line bundle,

\[
 \boxed{K\simeq O_D\cdot dz.}                   \tag{2.3}
\]

Consequently its kernel-direction map is constant even on the nonreduced
scheme `D`:

\[
 \deg_D(K)=0,\qquad e(\gamma)=0.                 \tag{2.4}
\]

Combining (2.3) with `LCCT1` gives

\[
 L\simeq O_D(D),\qquad
 \operatorname{ch}_2(\operatorname{coker}\theta_f^{\log})
 =\frac12D^2.                                   \tag{2.5}
\]

## 3. Why the divisor inequality is global

It is enough to compute the order of `f^*z` at one point where each
irreducible component of `D` is present and the normalized coefficient is a
unit.  Those are divisorial valuations and are constant along the component.
If the computed order is at least `m_i` on every component, then
`div(f^*z)-D` has no negative codimension-one coefficient.  Smoothness of
`X`, hence normality, makes the quotient by a local equation of `D` regular;
isolated points cannot create a hidden pole.

This also explains why zeros or poles of a chosen *formula* for a kernel
ratio are irrelevant here.  The covector `dz` itself is a fixed unimodular
element of the target logarithmic frame.

## 4. F2 extraction-root packet

At the node of the first exceptional component `E=(W=0)` and the strict line
at infinity `L=(U=0)`, `PF2UCE1` constructs target coordinates with

\[
 f^*T=WU^5\cdot\text{unit},\qquad
 f^*z=W^3U^{18}\cdot\text{unit}.                \tag{4.1}
\]

The logarithmic determinant has the same component multiplicities:

\[
 D_{\rm root}=3E+18L.                            \tag{4.2}
\]

The unit in (4.1) at `E intersect L` simultaneously proves

\[
 \operatorname{ord}_E(f^*z)=3,qquad
 \operatorname{ord}_L(f^*z)=18.                 \tag{4.3}
\]

Hence (2.2) holds on the complete cyclic root packet.  More explicitly, with
`s=W^3U^18` and `f^*z=sb`,

\[
 d(f^*z)=s\bigl(db+b(3\,d\log W+18\,d\log U)\bigr). \tag{4.4}
\]

Relative to the first target covector `dlog T`, the residual determinant at
the root node is `18-5*3=3`, a unit.  Thus the same calculation both recovers
unit `Fitt_1` and proves that `dz` spans the restricted kernel.

The Gauss-degree alternative from `LKGD1` therefore collapses to

\[
 \boxed{e_{\rm root}=0,\qquad
 \operatorname{ch}_2(\mathcal T^{\log}_{\rm root})=27.} \tag{4.5}
\]

No additional multiplicity-`3/18` jet calculation is needed: equation (4.4)
already works modulo the full ideal `(W^3U^18)`.  Ordinary boundary blowups
also preserve this conclusion by total transform.

## 5. Claim boundary and path forward

The theorem covers the cyclic root determinant packet wherever `Fitt_1` is a
unit and the fixed target coordinate remains in one logarithmic frame.  It
does not extend automatically through a point where the cokernel becomes
noncyclic, through a determinant component mapping onto a target curve, or
through an attachment for which the target image leaves this chart.

Thus the root twist/Gauss-degree problem is closed, but `(75,125)` is not.
The next required geometry is now external to the tame root packet:

1. compile the purity-forced affine ramification row;
2. compute full presentations at every noncyclic attachment and their
   codimension-two corrections; and
3. insert the exact positive root term `27` into the complete logarithmic
   Chern identity and prove whether the remaining terms can cancel it.

The subsequent
[`outgoing terminal-tail theorem`](F2_OUTGOING_TERMINAL_TAIL.md) closes the
formerly first item as a unimodular log-etale fan map with zero correction.

<!-- status-consumer: PF2OTT1 af25012e34020e11 -->

The theorem does not exclude `(75,125)` or prove `JC(2)`.

## Reproduction

```bash
.venv/bin/python scripts/verify_log_tangential_kernel_trivialization.py
```
