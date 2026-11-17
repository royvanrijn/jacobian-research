# Cyclic logarithmic boundary charge under blowups

> **Status.**  This note proves that the raw finite matching length at an SNC
> node is not birationally invariant, but its sum with the component
> self-intersection term is.  For a cyclic logarithmic cokernel supported on
> the Cartier boundary divisor `D=sum m_i C_i`, the conserved quantity is
> `D^2/2`, equivalently the untwisted Cartier contribution
> `-ch_2(O_D)`.  Applied to the F2 extraction-root cycle
> `D_root=3E+18L`, with `(E^2,L^2,E.L)=(-6,0,1)`, it gives the stable charge
> `D_root^2/2=27`.  This stabilizes the earlier raw length `54` under further
> source-boundary blowups.  Subsequent kernel-line, Gauss-degree, and
> tangential-coordinate theorems identify the actual F2 cyclic root
> contribution exactly as `27`; they do not prove the global Chern ledger
> contradictory.

The symbolic blowup identities, exact graph transformations, and F2 charge
are checked by
[`verify_log_cyclic_boundary_blowup_conservation.py`](../scripts/verify_log_cyclic_boundary_blowup_conservation.py).
The reusable graph calculation is in
[`log_node_profiles.py`](../jcsearch/log_node_profiles.py).

## 1. Cyclic matching at one node

Let `R` be a two-dimensional regular local ring with parameters `u,v`.  For
positive integers `a,b`, put

\[
 M_{a,b}=R/(u^av^b).                              \tag{1.1}
\]

Since `(u^a) intersect (v^b)=(u^av^b)`, there is a canonical exact sequence

\[
0\longrightarrow M_{a,b}
 \longrightarrow R/(u^a)\oplus R/(v^b)
 \longrightarrow R/(u^a,v^b)\longrightarrow0.   \tag{1.2}
\]

The final quotient has length

\[
 \ell_{u,v}(a,b)=ab.                              \tag{1.3}
\]

This is the branch-matching length.  The module `M_(a,b)` is a hypersurface
Cohen--Macaulay module, so it has no finite-support submodule.  Thus (1.3)
belongs to the degree-one normalization/matching sequence, not to
`H^0_(u,v)(M_(a,b))`.

## 2. Global SNC decomposition

Let `X` be a smooth surface and let

\[
 D=\sum_i m_iC_i                               \tag{2.1}
\]

be an effective Cartier divisor supported on an SNC curve with no triple
points.  The componentwise version of (1.2) gives

\[
0\longrightarrow O_D
 \longrightarrow\bigoplus_i O_{m_iC_i}
 \longrightarrow
 \bigoplus_{p\in C_i\cap C_j}O_{X,p}/(u^{m_i},v^{m_j})
 \longrightarrow0.                              \tag{2.2}
\]

Consequently the total finite matching length is

\[
 L_{\rm node}(D)=\sum_{C_i\cap C_j\ne\varnothing}m_im_j. \tag{2.3}
\]

This number alone depends on the selected boundary model.

## 3. The corrected Cartier charge

Define

\[
 \mathcal Q(D)=
 \frac12\sum_i m_i^2C_i^2
 +\sum_{C_i\cap C_j\ne\varnothing}m_im_j.       \tag{3.1}
\]

Expanding the self-intersection of (2.1) gives

\[
 \boxed{\mathcal Q(D)=\frac12D^2.}              \tag{3.2}
\]

The exact sequence

\[
0\longrightarrow O_X(-D)\longrightarrow O_X
 \longrightarrow O_D\longrightarrow0
\]

also gives

\[
 \operatorname{ch}(O_D)=1-e^{-D},\qquad
 \operatorname{ch}_2(O_D)=-\frac12D^2=-\mathcal Q(D). \tag{3.3}
\]

Thus (3.1) is precisely the negative codimension-two Chern character of the
untwisted cyclic Cartier cokernel.

## 4. Blowup of a boundary node

Suppose `C_1` and `C_2` meet, with multiplicities `a,b` and
self-intersections `s_1,s_2`.  Blow up their node.  The strict transforms
have self-intersections `s_1-1,s_2-1`; the exceptional curve has
self-intersection `-1` and total-transform multiplicity `a+b`.  The old
local contribution to `2Q` is

\[
 a^2s_1+b^2s_2+2ab.                              \tag{4.1}
\]

The new contribution is

\[
\begin{split}
 &a^2(s_1-1)+b^2(s_2-1)-(a+b)^2\\
 &\qquad+2a(a+b)+2b(a+b),
\end{split}                                      \tag{4.2}
\]

which expands exactly to (4.1).

Notice that the raw matching length changes from

\[
 ab\quad\text{to}\quad a(a+b)+b(a+b)=(a+b)^2.  \tag{4.3}
\]

The component self-intersection correction is therefore essential.

## 5. Blowup of a smooth boundary point

If a smooth point of a multiplicity-`a` component of self-intersection `s`
is blown up, the total transform has multiplicity `a` on both the strict
transform and exceptional curve.  Its doubled contribution is

\[
 a^2(s-1)-a^2+2a^2=a^2s.                         \tag{5.1}
\]

Together with Section 4, this proves invariance under every sequence of
ordinary admissible source-boundary blowups.

Equivalently, for the blowup `g:X'->X` and `D'=g^*D`, one has

\[
 (D')^2=D^2,
 \qquad Rg_*[O_{D'}]=[O_D]                       \tag{5.2}
\]

in numerical intersection theory and coherent `K`-theory.  The second
identity follows from the Cartier resolutions and the projection formula.

## 6. The F2 extraction-root charge

In the F2 carrier-extraction skeleton, let `E` be the first exceptional
component and `L` the strict line at infinity.  The exact source graph gives

\[
 E^2=-6,\qquad L^2=0,qquad E\cdot L=1.          \tag{6.1}
\]

`PF2UCE1` gives logarithmic determinant multiplicities `3` and `18`, so put

\[
 D_{\rm root}=3E+18L.                            \tag{6.2}
\]

The raw node length is `3*18=54`, while

\[
\begin{split}
 \mathcal Q(D_{\rm root})
 &=\frac12\bigl(3^2(-6)+18^2(0)\bigr)+3\cdot18\\
 &=-27+54\\
 &=\boxed{27}.
\end{split}                                      \tag{6.3}
\]

After blowing up the root node, the multiplicities are `(3,18,21)`, the
self-intersections are `(-7,-1,-1)`, and the raw matching length becomes

\[
 3\cdot21+18\cdot21=441.                         \tag{6.4}
\]

The self-intersection term becomes `-414`, so the corrected charge remains
`27`.  Thus the root packet has a stable nonzero Cartier charge even though
its raw finite matching length is model-dependent.

## 7. The cyclic twist and what remains

The actual logarithmic cotangent cokernel need only be locally cyclic.  If it
is globally `i_*L` along its determinant divisor, then the
[`kernel-line theorem`](LOG_CYCLIC_COKERNEL_TWIST.md) proves

\[
 L\simeq K\otimes O_D(D),                        \tag{7.1}
\]

where `K` is the kernel line of the differential restricted to `D`.  Hence

\[
 \operatorname{ch}_2(i_*L)
 =\deg_D(K)+\frac12D^2.                          \tag{7.2}
\]

Thus the former abstract twist ambiguity is the concrete degree of a
logarithmic kernel-direction line.  At the F2 extraction root the exact
contribution is

\[
 \operatorname{ch}_2(\mathcal T^{\log}_{\rm root})
 =\deg(K_{\rm root})+27.                         \tag{7.3}
\]

The local presentation `R/(W^3U^18)` does not determine this global degree.
Nor does either theorem cover nodes where `Fitt_1` is nontrivial and the
cokernel is not cyclic.

For this contracted root packet the
[`Gauss-degree theorem`](LOG_KERNEL_GAUSS_DEGREE.md) gives

\[
 \deg(K_{\rm root})=-e_{\rm root},\qquad
 \operatorname{ch}_2(\mathcal T^{\log}_{\rm root})=27-e_{\rm root}\le27,
 \tag{7.4}
\]

where `e_root` is the degree of a globally generated kernel-direction pencil.

For F2, the
[`tangential-coordinate theorem`](LOG_TANGENTIAL_KERNEL_TRIVIALIZATION.md)
then proves `e_root=0`: the fixed covector `dz` pulls back into the full ideal
of `D_root`.  Thus the cyclic root contribution is exactly `27`.

Therefore the next global tasks are:

1. compile the entire logarithmic determinant cycle and all noncyclic nodal
   corrections; the outgoing terminal tail is now closed, while affine purity
   forces a new component whose target curve, pullback factorization, and
   proximity chain remain to be recovered;
2. compare the resulting `ch_2` with the global logarithmic Chern identity;
3. prove that the exact root term `27` cannot be cancelled or absorbed by the
   remaining components.

Until these are done, neither the charge `27` nor the raw matching length
`54` excludes `(75,125)`.

<!-- status-consumer: LCCT1 2fc6ecea7a7c8b49 -->

<!-- status-consumer: LKGD1 8a357250b5005186 -->

<!-- status-consumer: LTKT1 32ac27318f16c20c -->

<!-- status-consumer: PF2APF1 192055eb737d3140 -->

## Reproduction

```bash
.venv/bin/python scripts/verify_log_cyclic_boundary_blowup_conservation.py
```
