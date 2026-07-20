# Canonical-family image and nonproperness theorem

This note completes the affine image calculation for the canonical primitives

\[
H_d(W)=W^d(1-W),\qquad d\ge2.
\]

It builds on the universal inverse and `S_n` monodromy theorem in
`WEIGHTED_SEED_THEOREM.md`.

## The canonical maps

Put

\[
p_d=H_d'=dW^{d-1}-(d+1)W^d,
\qquad
q_d=Wp_d-H_d=(d-1)W^d-dW^{d+1},
\]

and

\[
a_d=-{2d-1\over2(d-1)},\quad
u=1+xy,\quad
\gamma=1+a_dxy+x^2z,\quad W=u\gamma.
\]

The weighted construction gives

\[
G_d=\left(
{u+q_d(W)/\gamma^2\over x^2},
{1+p_d(W)/\gamma\over x},
x\gamma
\right).
\]

The quotients cancel, `G_d` is polynomial, `det DG_d=1`, and its generic
degree is `d+1`. For a target `(A,B,C)`, its inverse polynomial is

\[
E_d(W)=W^d(1-W)-BCW+AC^2.
\]

On `C!=0`, affine source points correspond exactly to the simple roots of
`E_d`.

## Saturated discriminant and nonproperness

Let

\[
D_d(s,t)=\operatorname{Disc}_W\bigl(W^d(1-W)-sW+t\bigr).
\]

For `d>=2`, substitution has exact `C`-adic order `d`:

\[
D_d(BC,AC^2)=C^dQ_d(A,B,C),
\]

which defines a polynomial `Q_d` without a hidden denominator stratum. At
`C=0`,

\[
Q_2(A,B,0)=B^2-4A,
\]

whereas for `d>=3`,

\[
Q_d(A,B,0)=\epsilon_d(d-1)^{d-1}B^d,
\qquad \epsilon_d\in\{+1,-1\}.
\]

Repeated roots on `C!=0` are reconstruction poles and yield escaping source
branches. Outside the pulled-back discriminant, all roots are finite and
simple and the reconstruction denominators are units. The projective inverse
has leading coefficient `-1`, so there is no unexamined root-at-infinity
chart.

The nonproperness sets are

\[
S_{G_2}=V(Q_2)
\]

and, for every `d>=3`,

\[
S_{G_d}=V(C)\cup V(Q_d).
\]

The difference at `d=2` is genuine. Its generic `C=0` fiber already contains
all three sheets. For `d>=3`, the direct `C=0` fiber has fewer than the generic
`d+1` points; étaleness prevents distinct bounded branches from coalescing, so
the remaining branches escape. Hence the entire plane `C=0` is a real
nonproperness component, not an artifact of clearing `C^d`.

## Direct fibers on `C=0`

The `x=0` chart maps bijectively to the target plane. Explicitly,

\[
B={y\over2(d-1)},\qquad
A={16d^2-19d+6\over8(d-1)}y^2-2(d-1)z.
\]

Thus every `(A,B,0)` has one distinguished finite source point.

On `gamma=0`, the cubic member `d=2` satisfies

\[
(4A-B^2)x^2=-1.
\]

It contributes two points when `4A!=B^2` and none when `4A=B^2`. Therefore
the `C=0` fiber sizes for `d=2` are respectively `3` and `1`.

For every `d>=3`, the same chart simplifies to

\[
A={1+xy\over x^2},\qquad B={1\over x}.
\]

It contributes one point when `B!=0` and none when `B=0`. Consequently

\[
\#G_d^{-1}(A,B,0)=
\begin{cases}
2,&B\ne0,\\
1,&B=0.
\end{cases}
\]

## Exact image

The lacunary polynomial

\[
P(W)=W^{d+1}-W^d+sW-t=-E_d(W)
\]

has no simple root only in the following two cases:

\[
\begin{array}{c|c|c}
d&(s,t)&P(W)\\ \hline
2&(1/3,1/27)&(W-1/3)^3\\
3&(1/8,-1/64)&(W^2-W/2-1/8)^2.
\end{array}
\]

Here is the uniform obstruction. Any root of multiplicity at least three must
also annihilate

\[
H_d''(W)=dW^{d-2}\bigl((d-1)-(d+1)W\bigr).
\]

The root `W=0` forces `(s,t)=(0,0)`, where `W=1` remains simple. The other
possibility is the triple root `(d-1)/(d+1)`. Hence a polynomial with no
simple root must be either a square of even degree or that triple factor times
a square of even residual degree.

Coefficient comparison at infinity finishes the classification. In the square
case, the forced coefficients are those of the truncation of
`sqrt(1-Z)`; for degree at least six its first omitted nonzero coefficient
would produce a forbidden term of degree at least two. In the triple-square
case they are the truncation of

\[
\sqrt{1-Z}\,(1-aZ)^{-3/2},
\qquad a={d-1\over d+1}.
\]

The first omitted coefficient vanishes only for `d=2`; for every even `d>=4`
it is strictly positive. Explicitly, if the last series is
`sum c_j Z^j`, then

\[
2(j+1)c_{j+1}
=\bigl(2j(1+a)+3a-1\bigr)c_j-2ajc_{j-1};
\]

substitution of `a=(2m-1)/(2m+1)` gives `c_m>0` for `m>=2` by induction
(while `c_1=0` for `m=1`). This leaves exactly the two displayed low-degree
polynomials.

Since every `C=0` target is already in the image, the affine images are

\[
G_2(\mathbb C^3)=\mathbb C^3\setminus
V(3BC-1,27AC^2-1),
\]

\[
G_3(\mathbb C^3)=\mathbb C^3\setminus
V(8BC-1,64AC^2+1),
\]

and the higher canonical maps are surjective:

\[
G_d(\mathbb C^3)=\mathbb C^3\qquad(d\ge4).
\]

For `C!=0`, the exact fiber size is the number of multiplicity-one roots of
`E_d`. Together with the direct `C=0` table above, this determines every
fiber without losing any boundary stratum.

Run:

```bash
.venv/bin/python scripts/verify_canonical_family_image.py
```

The exact audit covers `d=2,...,8`, verifies the boundary charts and saturated
discriminants, checks both exceptional omitted polynomials, and audits the
lacunary coefficient obstructions in every parity.
