# A cellular cotangent prototype for the Hessian--Ritt braid

> **Status.** This note constructs and verifies an actual cellular chain
> complex of coefficient modules.  It is an exact associated-graded model
> for the conormal data already certified in degrees thirty and forty-two.
> It is not yet a proof that this complex is quasi-isomorphic to the
> completed cotangent complex of every Hessian intersection.  In particular,
> it does not compute the degree-forty-two transitivity extension class or
> prove filtered \(H^2\)-vanishing in all orders.

The purpose of the prototype is to make the proposal in
[the Hessian--Ritt deformation-complex programme](HESSIAN_RITT_DEFORMATION_COMPLEX.md)
literal.  Vertices, Ritt moves, and Coxeter cells now carry modules, the
signed restriction maps are matrices, and \(H^1,H^2\) are computed from
their totalization.

## 1. The two coefficient blocks

Let \(X\) be the filled braid hexagon.  It has six vertices, six move edges,
and one braid two-cell.  Let \(B=\mathbb Q[[\tau,z]]\) be the completed
Dickson component and let \(T_B\) denote either its tangent module or a
finite Artin jet of its structure module.

The reduced coefficient block is the ordinary cellular complex

\[
 C^\bullet(X;T_B):
 T_B^6\xrightarrow{\delta_0}T_B^6
 \xrightarrow{\delta_1}T_B.                                  \tag{1.1}
\]

Orienting each move edge and the braid boundary gives
\(\delta_1\delta_0=0\).  Since \(X\) is a disk,

\[
 H^0(1.1)=T_B,\qquad H^1(1.1)=H^2(1.1)=0.                    \tag{1.2}
\]

Now take one three-edge half-braid
\(\gamma=(v_0,v_1,v_2,v_3)\) and a coefficient module \(D\) supported on
the power boundary.  Relative to its endpoints, its cellular complex is

\[
 C^\bullet(\gamma,\partial\gamma;D):
 D^{\{v_1,v_2\}}\xrightarrow{
 \left(\begin{smallmatrix}1&0\\-1&1\\0&-1\end{smallmatrix}\right)}
 D^{\{e_{01},e_{12},e_{23}\}}\longrightarrow0.               \tag{1.3}
\]

The first map is injective and its cokernel is canonically one copy of
\(D\).  Hence

\[
 H^0(1.3)=0,\qquad H^1(1.3)=D,\qquad H^2(1.3)=0.              \tag{1.4}
\]

The prototype totalization for one labelled sector is

\[
 \mathcal C_D^\bullet
 =
 C^\bullet(X;T_B)
 \oplus C^\bullet(\gamma,\partial\gamma;D).                  \tag{1.5}
\]

Thus the vertex deformation modules, move modules, and braid-cell module
are respectively

\[
 \begin{aligned}
 \mathcal C_D^0&=T_B^6\oplus D^2,\\
 \mathcal C_D^1&=T_B^6\oplus D^3,\\
 \mathcal C_D^2&=T_B,
 \end{aligned}                                               \tag{1.6}
\]

and

\[
 H^0(\mathcal C_D)=T_B,\qquad
 H^1(\mathcal C_D)=D,\qquad
 H^2(\mathcal C_D)=0.                                        \tag{1.7}
\]

Equation (1.5) is a genuine complex for any module \(D\).  The theorem still
to be proved is that the cellwise cotangent totalization of the original
coefficient correspondences reduces to (1.5), with the displayed \(D\), in
the derived completed category.

## 2. Commuting and braid cells

The implementation orients both Coxeter cell types.  For four coprime
factor degrees `(2,3,5,7)`, the two-skeleton has

\[
 24\text{ vertices},\quad36\text{ moves},\quad
 6\text{ commuting squares},\quad8\text{ braid hexagons}.    \tag{2.1}
\]

Its rational cellular matrices have ranks \(23,13\), so

\[
 (\dim H^0,\dim H^1,\dim H^2)=(1,0,1).                       \tag{2.2}
\]

The last class is topological: this two-skeleton is the boundary of the
three-dimensional permutohedron.  Its \(H^2\) is killed by the
permutohedron three-cell.  It must not be misidentified with a Ritt
obstruction class.  In the three-factor degree-thirty and degree-forty-two
braids, the filled hexagon is already a disk and has no such topological
\(H^2\).

## 3. Degree thirty

Reduced component completeness is known here: the all-cut reduced
intersection is the Dickson plane and its power boundary is \(z=0\).
For each of the three rotations, take \(D_j\) to be the
composite-omitting path-to-boundary defect module.  Its conormal fiber at
the monomial point is one-dimensional.  Applying (1.5) with
\(\dim T_B=2\) and \(\dim D_j=1\) gives exact matrix dimensions

\[
 \dim(\mathcal C^0,\mathcal C^1,\mathcal C^2)=(14,15,2),
 \qquad
 \operatorname{rank}(\delta_0,\delta_1)=(12,2),              \tag{3.1}
\]

and therefore

\[
 (\dim H^0,\dim H^1,\dim H^2)=(2,1,0)                        \tag{3.2}
\]

in every labelled sector.

The equality of these linear answers is exactly why completion is needed.
The established nonlinear decorations are:

| omitted composite cut | omitted prime cut | \(\operatorname{Ann}_B D_j\) known on the base | path nilpotence | transverse slice | slice Hilbert vector |
|---:|---:|---|---:|---|---|
| `10` | `3` | \((z^2)\) | `4` | \(\mathbb Q[u]/(u^5)\) | `(1,1,1,1,1)` |
| `15` | `2` | \((z^2)\) | `3` | \(\mathbb Q[u,v]/(u^2,v^2)\) | `(1,2,1)` |
| `6` | `5` | \((z^4)\) | `4` | \(\mathbb Q[u,v]/(u^4,v^2)\) | `(1,2,2,2,1)` |

The point-cotangent homology pairs of the three transverse slices are
\((1,1),(2,2),(2,2)\).  These are local complete-intersection invariants,
not the cellular \(H^1\) in (3.2).  The latter is the single path mismatch;
the former records the number of transverse generators and relations after
fixing the base.

This separation prevents two common errors:

1. the cellular \(H^1\)-rank does not determine the nilpotence index; and
2. \(H^2=0\) for the displayed linear totalization does not reconstruct the
   completed algebra.

## 4. Degree forty-two: sector and spectator

On the `2 o 7 o 3` chart the exact completed ideal flag is

\[
 I_6\subsetneq I_7=I_{\partial}\subsetneq K.                 \tag{4.1}
\]

Put

\[
 D_{\rm sec}=I_{\partial}/I_6,\qquad
 D_{\rm sp}=K/I_{\partial}.                                  \tag{4.2}
\]

These are the two coefficient modules underlying the successive conormal
layers.  The exact annihilation calculations give

\[
 z^8D_{\rm sec}=0,\quad z^7D_{\rm sec}\ne0,\qquad
 zD_{\rm sp}=0,\quad D_{\rm sp}\ne0.                         \tag{4.3}
\]

Each has a one-dimensional first conormal fiber.  The associated-graded
cellular model is

\[
 \operatorname{gr}\mathcal C^\bullet_{42}
 =
 C^\bullet(X;T_B)
 \oplus C^\bullet(\gamma,\partial\gamma;D_{\rm sec})
 \oplus C^\bullet(\gamma,\partial\gamma;D_{\rm sp}).         \tag{4.4}
\]

At first conormal order its matrices have dimensions and ranks

\[
 (16,18,2),\qquad(14,2),                                     \tag{4.5}
\]

so

\[
 (\dim H^0,\dim H^1,\dim H^2)=(2,2,0),                       \tag{4.6}
\]

with the canonical associated-graded identification

\[
 H^1\cong D_{\rm sec}\oplus D_{\rm sp}.                      \tag{4.7}
\]

The known completed jets give a stronger regression.  For
\(\mathfrak m=(n_1,\ldots,n_7,\tau,z)\), the coefficient dimensions through
orders \(q=1,2,3\) are:

| \(q\) | Dickson base | sector layer | spectator layer | cellular \((H^0,H^1,H^2)\) |
|---:|---:|---:|---:|---|
| 1 | 1 | 0 | 0 | `(1,0,0)` |
| 2 | 3 | 1 | 1 | `(3,2,0)` |
| 3 | 6 | 5 | 3 | `(6,8,0)` |

The order-three values \(5\) and \(3\) show that the two rank-one conormal
directions immediately acquire different higher structure.  In particular,
neither layer should be replaced without proof by a bare cyclic module
\(B/(z^m)\).

The direct sum in (4.4) is the associated graded of the exact flag (4.1).
The completed transitivity triangle may carry a nontrivial extension
coupling \(D_{\rm sec}\) and \(D_{\rm sp}\).  Nothing in (4.4)--(4.7)
asserts that it splits.

## 5. What the prototype establishes

The executable calculation proves the following finite statements.

1. The signed vertex-to-move and move-to-cell matrices satisfy
   \(\delta_1\delta_0=0\) for every constructed commuting and braid cell.
2. The degree-thirty sector totalizations have \(H^1\)-dimension one and
   \(H^2=0\).
3. The separated degree-forty-two associated graded has two
   one-dimensional conormal \(H^1\)-summands, with exact base annihilation
   exponents \(8\) and \(1\).
4. The order-two and order-three totalizations reproduce the exact jet
   differences from the completed ideal flag.

This is enough to isolate a reusable theorem target:

> **Cellular reduction target.** After completion along a reduced
> power/Dickson component, the cellwise cotangent homotopy limit is filtered
> quasi-isomorphic to a reduced-component cellular block plus relative path
> blocks of the form (1.3), with the transitivity extensions and higher
> brackets retained.

Proving that target would make (1.7) formal and would reduce
degree-specific synchronization to three universal tasks: identify the
move coefficient modules, calculate the extension between adjacent
filtration layers, and prove vanishing of the genuinely algebraic
\(H^2\)-classes.  The present computation verifies the shape but not that
general comparison theorem.

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_hessian_ritt_cellular_cotangent_prototype.py
```

The command writes
`artifacts/generated-results/hessian_ritt_cellular_cotangent_prototype.json`.
