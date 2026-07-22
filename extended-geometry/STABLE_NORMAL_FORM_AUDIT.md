# Stable normal-form audit and explicit consequences

The 190-variable Image and quartic HN witnesses below remain valid artifacts
of the conservative 95-variable reduction.  They are no longer the best
repository dimensions: the [essential 21-variable construction](IMAGE_VANISHING_COUNTEREXAMPLES.md)
gives a direct `SIC(21)` witness and one explicit 42-variable quartic for both
the generalized Laplacian and classical HN Vanishing Conjectures.

## Corrected theorem

Over any characteristic-zero field, the normalized foundational Keller map is connected to an
explicit 95-dimensional map `K=I+H`, with `H` cubic homogeneous, by 22
polynomial stable equivalences, one Segre equivalence, and one final polynomial
stable equivalence. Thus the precise phrase is **stable--Segre equivalent**;
the Segre step is not, on the cited evidence, a polynomial stable equivalence.
The stored collision is transported directly, so this distinction does not
affect the counterexample.

The 95-dimensional map is GZ-paired, hence polynomially stably equivalent, to
an explicit Druzkowski map

\[
G(X)=X-(AX)^{*3}:\mathbb C^{451}\longrightarrow\mathbb C^{451},
\qquad \operatorname{rank}A=95.
\]

Both maps have determinant one and the artifacts contain three distinct exact
rational points with one common image.

## Independent regeneration

The production scripts use SymPy. `audit_stable_normal_form_independent.py` instead implements
sparse polynomial arithmetic and rational Gaussian elimination using only the
Python standard library. Starting with the three foundational map coordinates, it recomputes
the determinant, executes all 22 stable extensions without reading the stored
step list, reconstructs all 148 terms of `H`, transports the collision,
independently polarizes `H` into 415 cubes, and verifies the GZ matrices.

## Dimension audit

The original construction appended all 95 columns of `I_95` to `B_0`, giving
dimension 510. Since `rank(B_0)=59`, exactly 36 independent columns suffice to
make `B` surjective. The audited dimension is therefore

\[
N=415+(95-59)=451.
\]

This is minimal among extensions of this fixed 415-column `B_0`. It is not
asserted to be minimal among all Waring decompositions or Druzkowski reductions.
No smaller cubic-homogeneous construction than 95 is certified here.

## Primary-source hypothesis audit

- Campbell, Theorem 5, applies to a nondegenerate polynomial map and explicitly
  uses the product-removal formula, normalization, Segre step, and final stable
  step used here. The paper states that the proof works over `C` and for Keller
  maps. The foundational Keller map is polynomial, normalized to fix zero with identity linear part,
  and has determinant one.
- Campbell, Theorem 4, says GZ-paired maps are polynomially stably equivalent
  over any ground field. The artifact verifies `BC=I`, `ker B=ker A`, and
  `K(x)=BG(Cx)`.
- Gorni--Zampieri, Proposition 2.1, permits adding only a few columns to make
  `B_0` full rank; its identity block is merely an example. Here `B` and `D`
  have rank 95, `C` is a right inverse, and `A=DB`, so `ker A=ker B`.
- Gorni--Zampieri, Proposition 2.4 supplies the determinant comparison; the
  artifact also proves it directly with Sylvester's identity.

Primary sources: [Campbell](https://arxiv.org/abs/1303.3853) and
[Gorni--Zampieri](https://arxiv.org/abs/1204.4026).

## Explicit Dixmier witness

Let `F=(F_3/2,F_2,F_1)` be the normalized three-dimensional foundational Keller map and
`J=JF`. In the third Weyl algebra, with `[d_i,x_j]=delta_ij`, define

\[
\Psi(x_i)=F_i(x),\qquad
\Psi(d_i)=\delta_i:=\sum_j(J^{-1})_{ji}d_j.
\]

Since `det J=1`, the coefficients are polynomial. Exact calculation verifies
`delta_i(F_k)=delta_ik` and `[delta_i,delta_k]=0`, so this is an endomorphism.
It is not an automorphism: otherwise `Psi^{-1}(x_j)` would commute with every
`x_i`, hence lie in `C[x]` by PBW. Applying `Psi` would produce a polynomial
left inverse to `F`, contradicting the foundational rational collision. All nine coefficients are
stored in `artifacts/generated-results/stable_normal_form_consequence_witnesses.json`. No external
Dixmier--Jacobian reduction is needed.

## Explicit quartic Hessian-nilpotent witness

Write `K(x)=x+H(x)`. In 190 variables `(u,v)`, set

\[
x=(u+iv)/\sqrt2,\quad y=(u-iv)/\sqrt2,\quad
R(u,v)=-\sum_{j=1}^{95}y_jH_j(x).
\]

The square roots cancel because `R` is quartic; the artifact expands it into
2012 terms over `Q(i)`. A block determinant for
`x.y+t y.H(x)`, followed by the orthogonal change of variables, gives

\[
\det(I+t\operatorname{Hess}R)=\det(I+tJH)^2=1.
\]

Thus `Hess(R)` is nilpotent. The map `Z-grad(R)` has the transported collision
`T^t(p_i,0)`, so it is not invertible. Zhao's Theorem 7.2 and formula (3.8)
then imply that `Delta^m R^(m+1)` is not eventually zero. This is an explicit
190-variable quartic witness. [Zhao's paper](https://arxiv.org/abs/math/0409534)
requires precisely a complex homogeneous quartic HN polynomial.

## Explicit Special Image-Conjecture witness

With commuting variables `(zeta,z)` in dimension 190, define

\[
f=(\zeta_1^2+\cdots+\zeta_{190}^2)R(z),\qquad g=R(z).
\]

For `E(zeta^alpha z^beta)=partial_z^alpha z^beta`, one has
`E(f^m)=Delta^m R^m=0` for every positive `m`, whereas
`E(gf^m)=Delta^m R^(m+1)` is not eventually zero. Theorem 3.3 of
[van den Essen--Wright--Zhao](https://arxiv.org/abs/1008.3962) identifies
`ker E` with the image of the commuting operators `zeta_i-partial_i`.
Therefore this pair explicitly violates the special Image Conjecture in
dimension 190.

The least exceptional exponent and dimension minimality remain open.

## Reproduction

```bash
make verify-derived
```
