# Completion of the degree-forty-two cellular cotangent prototype

> **Status.** The factor/move/labelled-cell diagrams in degrees forty-two
> and thirty, their cellular totalizations, and their filtration
> cohomology are computed exactly.  The split associated-graded prototype
> has \(H^2=0\) in every filtration degree.  It is nevertheless not the
> filtered cotangent complex: in degree forty-two its sector--spectator
> coefficient extension first becomes non-split modulo
> \(\mathfrak m^3\), and the completed conormal projection has a nonzero
> cotangent-transitivity connecting morphism.  This is the requested
> explicit higher obstruction to the naive cellwise reduction.  The
> extension-retaining Postnikov tower survives this obstruction; global
> coefficient effectivity on the other five degree-forty-two charts and
> the individual cotangent \(H_i\), \(i\ge2\), remain open.

This note closes the literal prototype requested in
[the deformation-complex programme](HESSIAN_RITT_DEFORMATION_COMPLEX.md).
It collects the exact HRCELL1--HRCELL5 calculations without promoting the
associated graded to an unproved derived equivalence.

## 1. The labelled factor and move diagrams

Write a factor word outer-to-inner.  In degree forty-two the vertices are

\[
237,\ 273,\ 327,\ 372,\ 723,\ 732.
\]

The six oriented adjacent Ritt moves, in the executable ordering, are

\[
\begin{array}{lll}
237\to327,&237\to273,&273\to723,\\
327\to372,&372\to732,&723\to732.
\end{array}                                                \tag{1.1}
\]

There is one labelled braid two-cell.  Its two half-braids are

\[
\begin{aligned}
237&\to327\to372\to732,\\
237&\to273\to723\to732.                                    \tag{1.2}
\end{aligned}
\]

With the vertex and edge order just displayed, the scalar cellular
coboundaries are

\[
\delta _0=
\begin{pmatrix}
-1&0&1&0&0&0\\
-1&1&0&0&0&0\\
0&-1&0&0&1&0\\
0&0&-1&1&0&0\\
0&0&0&-1&0&1\\
0&0&0&0&-1&1
\end{pmatrix},
\qquad
\delta _1=
\begin{pmatrix}1&-1&-1&1&1&-1\end{pmatrix}.                \tag{1.3}
\]

Direct multiplication gives \(\delta _1\delta _0=0\), and the ranks are
\((5,1)\).  Thus the filled braid is a cellular disk.

The degree-thirty control is obtained by replacing \(7\) by \(5\):

\[
\begin{aligned}
&235,\ 253,\ 325,\ 352,\ 523,\ 532,\\
&235\to325\to352\to532,\qquad
235\to253\to523\to532.                                     \tag{1.4}
\end{aligned}
\]

Its scalar matrices are again (1.3).  The labels, rather than the abstract
incidence, distinguish the two degrees.

## 2. The labelled power-boundary block

Let \(D\) be a coefficient module on one labelled power half-braid.
Relative to the common endpoints, its cellular block is

\[
D^2\xrightarrow{P_D}D^3\longrightarrow0,\qquad
P_D=
\begin{pmatrix}
1&0\\
-1&1\\
0&-1
\end{pmatrix}\otimes1_D.                                  \tag{2.1}
\]

The map is split-injective as an underlying abelian-group map and has
cokernel \(D\).  Consequently

\[
H^0=0,\qquad H^1=D,\qquad H^2=0                           \tag{2.2}
\]

for every \(D\), without a dimension or flatness assumption.

Let \(T_B\) be the completed Dickson-base coefficient.  The degree-forty-two
associated-graded totalization is therefore

\[
\begin{aligned}
C^0&=T_B^6\oplus D_{\rm sec}^2\oplus D_{\rm sp}^2,\\
C^1&=T_B^6\oplus D_{\rm sec}^3\oplus D_{\rm sp}^3,\\
C^2&=T_B,                                                  \tag{2.3}\\
d^0&=\operatorname{diag}
(\delta _0\otimes T_B,P_{D_{\rm sec}},P_{D_{\rm sp}}),\\
d^1&=(\delta _1\otimes T_B,0,0).
\end{aligned}
\]

Hence

\[
H^0=T_B,\qquad
H^1=D_{\rm sec}\oplus D_{\rm sp},\qquad
H^2=0.                                                      \tag{2.4}
\]

The exact ideal flag is

\[
I_6\subset I_\partial\subset K,\qquad
D_{\rm sec}=I_\partial/I_6,\quad D_{\rm sp}=K/I_\partial,
                                                                  \tag{2.5}
\]

with minimal base annihilators \(z^8\) and \(z\), respectively.

In degree thirty \(I_\partial=K\), so \(D_{\rm sp}=0\).  For each of its
three labelled sectors,

\[
C^0=T_B^6\oplus D_j^2,\quad
C^1=T_B^6\oplus D_j^3,\quad
C^2=T_B,                                                    \tag{2.6}
\]

and \((H^0,H^1,H^2)=(T_B,D_j,0)\).  This is the sector-only control,
not a split two-layer example.

## 3. First Postnikov layers

The degree-forty-two maximal-ideal quotients give:

| filtration order \(q\) | \(\dim T_B\) | \(\dim D_{\rm sec}\) | \(\dim D_{\rm sp}\) | \((\dim H^0,\dim H^1,\dim H^2)\) |
|---:|---:|---:|---:|---:|
| 1 | 1 | 0 | 0 | \((1,0,0)\) |
| 2 | 3 | 1 | 1 | \((3,2,0)\) |
| 3 | 6 | 5 | 3 | \((6,8,0)\) |
| 4 | 10 | 13 | 6 | \((10,19,0)\) |

Because (1.3) is a disk complex and (2.1) has no degree-two term, the same
argument proves

\[
\boxed{H^2_{\rm cell}(\operatorname{gr}_q)=0
\text{ for every filtration degree }q.}                    \tag{3.1}
\]

There is therefore no “first nonzero \(H^2\)” inside the split prototype.
Any nonzero \(H^2\) of the actual derived intersection must occur in an
internal higher-cotangent row of the cellular spectral sequence, or in the
mapping cone measuring failure of coefficient effectivity.  It cannot be
manufactured from the incidence topology of the filled three-factor braid.

For degree thirty the three conormal controls all have

\[
(\dim C^0,\dim C^1,\dim C^2)=(14,15,2),\quad
\operatorname{rank}(d^0,d^1)=(12,2),\quad
(\dim H^0,\dim H^1,\dim H^2)=(2,1,0).                      \tag{3.2}
\]

Their different annihilators and transverse algebras remain as recorded in
HRCELL1; none creates cellular \(H^2\).

## 4. Comparison with HRCELL1--HRCELL5

The five calculations now have a precise division of labor.

| result | retained information |
|---|---|
| HRCELL1 | split associated-graded cellular matrices (2.3) |
| HRCELL2 | the first non-split sector--spectator jets |
| HRCELL3 | completed ideal-module non-splitting |
| HRCELL4 | nonzero cotangent-transitivity connecting morphism |
| HRCELL5 | zero first Postnikov overlap; finite-jet excess is ordinary base-change Tor |

At order two the \(1+1\) sector--spectator extension splits.  At order
three,

\[
0\to\mathbb Q^5\to\mathbb Q^8\to\mathbb Q^3\to0           \tag{4.1}
\]

does not split equivariantly.  The change-of-splitting map has rank seven;
adjoining the actual coupling raises the rank to eight.  The primitive
functional

\[
\ell(C)=C_z[0,1]+10C_z[1,1]                               \tag{4.2}
\]

annihilates every coboundary and evaluates to \(1\).  At order four the
corresponding ranks are \(52\) and \(53\), and the stored primitive
functional evaluates to \(-3\).

Presentation-first reduction of the completed conormal projection gives

\[
0\to\mathbb Q^4\to\mathbb Q^6\to\mathbb Q^2\to0.           \tag{4.3}
\]

Its change-of-splitting rank is five and its augmented rank is six.  The
functional

\[
2C_z[0,0]+5C_z[1,0]                                       \tag{4.4}
\]

vanishes on all coboundaries and evaluates to \(2\).  Hence the completed
conormal projection has no module section and the connecting map in the
cotangent transitivity triangle is nonzero.

HRCELL5 separately proves

\[
\frac{I_\partial\cap(I_6+K^2)}{I_6+KI_\partial}=0,          \tag{4.5}
\]

so this obstruction is not a hidden first-Postnikov overlap.  Nor is it the
two-dimensional kernel created after a non-flat finite base change; that
kernel is an ordinary Tor image.

## 5. Success criterion and exact boundary

The split cellwise homotopy limit does **not** reduce to (2.3) as a filtered
object.  The exact failure is:

\[
\boxed{
\text{filtration order }3:
\quad
K/I_6\not\simeq
(I_\partial/I_6)\oplus(K/I_\partial),
\quad
\partial_{\rm cot}\ne0.}                                  \tag{5.1}
\]

This is a higher obstruction in the coefficient/Postnikov direction, not a
topological \(H^2\) of the braid disk.  Consequently the correct possible
reduction target is the extension-retaining Postnikov tower of
HRCELL2--HRCELL5, not the direct-sum HRCELL1 associated graded.

The full bar cotangent complex and its cellular subdivision are formal by
[cotangent descent](HESSIAN_RITT_COTANGENT_DESCENT_COMPARISON.md).  What is
not proved is that the actual completed coefficient diagram on all six
degree-forty-two charts factors coherently through this finite Ritt face
category.  Computing the transported flags on the other five charts is
still required before claiming a global degree-forty-two
quasi-isomorphism.

## Reproduction

Run

```bash
.venv/bin/python scripts/verify_ritt_cellular_prototype_completion.py
```

The command consumes the pinned HRCELL2, HRCELL4, and HRCELL5 artifacts and
writes
`artifacts/generated-results/ritt_cellular_prototype_completion.json`.
