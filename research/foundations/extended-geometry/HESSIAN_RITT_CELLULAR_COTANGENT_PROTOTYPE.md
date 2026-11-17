# A cellular cotangent prototype for the Hessian--Ritt braid

> **Status.** This note constructs and verifies an actual cellular chain
> complex of coefficient modules.  It is an exact associated-graded model
> for the conormal data already certified in degrees thirty and forty-two.
> It is not yet a proof that this complex is quasi-isomorphic to the
> completed cotangent complex of every Hessian intersection.  In particular,
> the degree-forty-two ideal-module extension is proved non-split by a
> finite tensor obstruction, and its conormal image proves that the
> cotangent transitivity connecting morphism is nonzero.  Its first
> homology sequence is now proved short exact, with the apparent truncated
> kernel separated as base-change Tor.  The individual higher cotangent
> homology modules, the cellwise derived comparison, and filtered
> \(H^2\)-vanishing are not proved.

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

Here the prime-omitting path is already the full Dickson boundary.  Thus
the spectator relative cotangent term is zero in degree thirty: (3.2) is a
sector-only block, not a split sector--spectator extension.  This is the
degenerate predecessor of the nonzero degree-forty-two transitivity class
in (4.22).

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
| 4 | 10 | 13 | 6 | `(10,19,0)` |

The order-three values \(5\) and \(3\) show that the two rank-one conormal
directions immediately acquire different higher structure.  In particular,
neither layer should be replaced without proof by a bare cyclic module
\(B/(z^m)\).

The direct sum in (4.4) is the associated graded of the exact flag (4.1).
It does **not** split at the next computed order.  Put

\[
 D_{\rm tot}=K/I_6.
\]

Modulo \(\mathfrak m^3\), exact standard-monomial reduction gives

\[
 \dim(D_{\rm sec},D_{\rm sp},D_{\rm tot})=(5,3,8)
\]

and the sequence

\[
 0\longrightarrow D_{\rm sec}\longrightarrow D_{\rm tot}
 \longrightarrow D_{\rm sp}\longrightarrow0                 \tag{4.8}
\]

is non-split as a \(\mathbb Q[\tau,z]\)-module.  In adapted rational bases,
the \(\tau\)-coupling is zero and the only nonzero block is

\[
 C_z=
 \begin{pmatrix}
 0&1&10\\
 0&0&0\\
 0&0&0\\
 0&0&0\\
 0&0&0
 \end{pmatrix}.                                               \tag{4.9}
\]

Changing a vector-space splitting by
\(H\in\operatorname{Hom}_{\mathbb Q}(D_{\rm sp},D_{\rm sec})\)
changes the pair of coupling blocks by

\[
 \left(
 T_{\rm sec}H-HT_{\rm sp},
 Z_{\rm sec}H-HZ_{\rm sp}
 \right).                                                     \tag{4.10}
\]

The coboundary map (4.10) has rank seven.  Appending the cocycle (4.9)
raises the rank to eight.  More explicitly, the functional

\[
 \ell(C)=C_z[0,1]+10C_z[1,1]                                 \tag{4.11}
\]

annihilates every coboundary and satisfies \(\ell(C)=1\).  This is an exact
finite-jet extension certificate.  At order two all base actions vanish and
the corresponding \(1+1\) extension splits.

Thus the correct order-three cellular coefficient is the non-split module
\(D_{\rm tot}\), not merely \(D_{\rm sec}\oplus D_{\rm sp}\).  Its
underlying vector-space totalization still has
\((\dim H^0,\dim H^1,\dim H^2)=(6,8,0)\), but its \(H^1\) now retains the
first scheme-theoretic coupling.  The calculation does not yet prove that
this class persists in the untruncated completion or identifies it with the
full derived cotangent transitivity class.

The fourth jet confirms persistence one order further:

\[
 \dim(D_{\rm sec},D_{\rm sp},D_{\rm tot})=(13,6,19).
                                                                  \tag{4.12}
\]

The change-of-splitting map now has rank \(52\), while adjoining the
coupling raises the rank to \(53\).  A primitive obstruction functional
supported on eight \(\tau\)-entries and two \(z\)-entries evaluates to
\(-3\).  The complete matrices and functional are stored in the generated
artifact.  This second non-splitting rules out an order-three anomaly, but
still does not prove formal persistence at every order.

### 4.1 A formal tensor obstruction

The jet calculations above use images inside truncated quotient rings.
For a formal splitting theorem the invariant operation is instead to
present the modules first and then tensor their presentations.

Let

\[
 R=\mathbb Q[[n_0,\ldots,n_6,\tau,z]]
\]

and present \(D_{\rm tot}=K/I_6\) and
\(D_{\rm sp}=K/I_\partial\) using the seven normal generators of \(K\).
Singular's exact `modulo` operation computes the relation modules.  Tensor
both presentations with

\[
 Q=
 R/\bigl((\tau,z)^2+(n_0,\ldots,n_6)^2\bigr).                \tag{4.13}
\]

The resulting surjection has dimensions

\[
 0\longrightarrow\mathbb Q^8
 \longrightarrow\mathbb Q^{12}
 \longrightarrow\mathbb Q^4
 \longrightarrow0.                                          \tag{4.14}
\]

All nine coordinate actions commute and the projection intertwines them.
A section must therefore solve the projection equations together with nine
families of commutator equations.  The change-of-splitting map has rank
\(27\); adjoining the actual coupling raises the rank to \(28\).  A
primitive annihilating functional is supported on

\[
\begin{aligned}
 &25C_{n_0}[4,0]+10C_{n_0}[5,1]
 -20C_{n_1}[4,0]\\
 &\qquad
 +4C_{n_1}[4,1]-10C_{n_1}[5,0],
\end{aligned}                                                \tag{4.15}
\]

and evaluates to \(240\).  Hence the tensor-quotient surjection has no
\(Q\)-linear section compatible with the \(R\)-action.

This proves the untruncated consequence:

\[
 \boxed{
 0\to I_\partial/I_6\to K/I_6\to K/I_\partial\to0
 \text{ is non-split over }\widehat R.}                       \tag{4.16}
\]

Indeed, an \(\widehat R\)-linear splitting would remain a splitting after
tensoring with (4.13), contradicting (4.15).  This is a theorem about the
completed ideal-module flag.  Its cotangent meaning is detected after
passing to the conormal modules below.

### 4.2 The conormal shadow of transitivity

For a surjection \(A\to B=A/J\), the first cotangent homology is the
conormal module \(H_1(L_{B/A})=J/J^2\).  Hence the two right-hand terms of
the degree-forty-two transitivity triangle have

\[
\begin{aligned}
 N_6&=H_1(L_{B/A_6})=K/(I_6+K^2),\\
 N_\partial&=H_1(L_{B/A_\partial})
             =K/(I_\partial+K^2).
\end{aligned}                                                \tag{4.17}
\]

The map in the triangle induces the canonical conormal projection
\(\pi:N_6\to N_\partial\).  Presenting both modules before reduction and
then tensoring over \(B=\mathbb Q[[\tau,z]]\) with

\[
 Q_B=B/(\tau,z)^2                                             \tag{4.18}
\]

gives the exact finite-dimensional sequence

\[
 0\longrightarrow\mathbb Q^4
 \longrightarrow\mathbb Q^6
 \overset{\pi\otimes Q_B}{\longrightarrow}\mathbb Q^2
 \longrightarrow0.                                           \tag{4.19}
\]

All seven normal actions vanish, as they must on a \(B\)-module.  A
\(B\)-linear section of the completed projection would induce a section of
(4.19).  The finite section equations instead have change-of-splitting
rank \(5\), while adjoining the actual coupling raises the rank to \(6\).
The primitive functional

\[
 2C_z[0,0]+5C_z[1,0]                                         \tag{4.20}
\]

annihilates every coboundary and evaluates to \(2\).  Thus the completed
conormal projection has no \(B\)-linear section.

Let

\[
 B\otimes^{\mathbf L}_{A_\partial}L_{A_\partial/A_6}
 \longrightarrow L_{B/A_6}\longrightarrow
 L_{B/A_\partial}\overset{\partial}{\longrightarrow}
 \left(B\otimes^{\mathbf L}_{A_\partial}
 L_{A_\partial/A_6}\right)[1]                                \tag{4.21}
\]

be transitivity.  If \(\partial=0\), the distinguished triangle splits, so
the middle-to-right map has a derived section and \(\pi\) has a section on
\(H_1\).  Equation (4.20) rules this out.  Therefore

\[
 \boxed{\partial\ne0.}                                       \tag{4.22}
\]

This identifies a nonzero first-Postnikov shadow of the actual cotangent
connecting morphism.  It does not compute the higher cotangent homology or
prove that the cellular totalization is the full cotangent homotopy limit.

### 4.3 The quadratic overlap and base-change Tor

Write \(I=I_6\) and \(J=I_\partial\).  Since a surjection has zero
cotangent \(H_0\), the first homology of the left term in (4.21) is

\[
 S=H_1\!\left(
 B\otimes^{\mathbf L}_{A_\partial}L_{A_\partial/A_6}
 \right)
 =J/(I+KJ).                                                   \tag{4.23}
\]

The kernel of its map to \(N_6\) is the quadratic overlap

\[
 Q_{\rm ov}=
 \frac{J\cap(I+K^2)}{I+KJ}.                                  \tag{4.24}
\]

This is a finite \(B\)-module because \(KQ_{\rm ov}=0\).  Exact ideal
intersection and Artin--Rees cutoff calculations give

\[
 Q_{\rm ov}/(\tau,z)^2Q_{\rm ov}=0.                          \tag{4.25}
\]

Concretely, at cutoff \(\mathfrak m^5\) the numerator and denominator
quotient rings both have length \(38\), and the \(1286\)-generator
intersection with \(\mathfrak m^5\) reduces identically to zero modulo the
denominator.  Nakayama's lemma therefore gives \(Q_{\rm ov}=0\) after
completion.  Thus the first transitivity homology sequence is

\[
 \boxed{
 0\longrightarrow S\longrightarrow N_6
 \longrightarrow N_\partial\longrightarrow0.}               \tag{4.26}
\]

There is a useful warning about finite jets.  A separate cutoff-four
certificate computes

\[
 \dim_{\mathbb Q}S/(\tau,z)^2S=22-16=6,                      \tag{4.27}
\]

whereas (4.19) has a four-dimensional kernel.  Tensoring (4.26) with
\(Q_B\) is not left exact:

\[
 \operatorname{Tor}_1^B(Q_B,N_\partial)\longrightarrow
 S\otimes_BQ_B\longrightarrow N_6\otimes_BQ_B.               \tag{4.28}
\]

The two-dimensional kernel in (4.28) is the image of ordinary base-change
Tor.  It is not an image of higher cotangent homology in the completed
\(H_1\) sequence.  Equivalently, the map
\(H_2(L_{B/A_\partial})\to S\) in the transitivity long exact sequence has
zero image.  The individual \(H_i\) for \(i\ge2\) may still be nonzero and
are not computed here.

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
5. The order-two sector--spectator extension splits, while the order-three
   and order-four extensions are non-split with exact obstruction
   functionals.
6. A presentation-first finite tensor quotient proves that the completed
   ideal-module extension itself is non-split.
7. Passing to the actual conormal modules gives a non-split
   \(4\to6\to2\) quotient and proves that the relative cotangent
   transitivity connecting morphism is nonzero.
8. The completed quadratic overlap vanishes, making the first transitivity
   homology sequence short exact; the two missing dimensions after
   base-square reduction are ordinary base-change Tor.

This is enough to isolate a reusable theorem target:

> **Cellular reduction target.** After completion along a reduced
> power/Dickson component, the cellwise cotangent homotopy limit is filtered
> quasi-isomorphic to a reduced-component cellular block plus relative path
> blocks of the form (1.3), with the transitivity extensions and higher
> brackets retained.

The ideal-flag and arbitrary-length Postnikov portions of this target are
now proved in
[cellular Postnikov transitivity](CELLULAR_POSTNIKOV_TRANSITIVITY.md).
They give the universal overlap formula, the multi-layer conormal
filtration, finite-quotient splitting detection, and the separation of
base-change Tor.  The
[cotangent-descent comparison](HESSIAN_RITT_COTANGENT_DESCENT_COMPARISON.md)
now proves descent for the full bar diagram of the actual derived
intersection and shows that the finite two-skeleton controls \(H^0,H^1\).
What remains here is coefficient effectivity of the finite Ritt compression
and the higher-cell calculation needed for algebraic \(H^2\).

Proving this remaining comparison would make (1.7) formal and would reduce
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

The finite-jet extension is replayed separately by

```bash
.venv/bin/python scripts/verify_degree42_ritt_cellular_extension.py
```

Its cached source residuals can be rebuilt, rather than silently refreshed,
with

```bash
.venv/bin/python scripts/research_degree42_cellular_extension.py \
  --order 2 --rebuild-source
```

The formal tensor obstruction is replayed by

```bash
.venv/bin/python scripts/verify_degree42_ritt_tensor_extension.py
```

The conormal transitivity obstruction is replayed by

```bash
.venv/bin/python scripts/verify_degree42_ritt_conormal_transitivity.py
```

The quadratic-overlap and base-change-Tor separation is replayed by

```bash
.venv/bin/python scripts/verify_degree42_ritt_postnikov_overlap.py
```

The arbitrary-length finite-module tower is replayed by

```bash
.venv/bin/python scripts/verify_cellular_postnikov_transitivity.py
```
