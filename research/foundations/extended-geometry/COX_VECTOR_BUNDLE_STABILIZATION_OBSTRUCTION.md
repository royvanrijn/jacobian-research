# Cox vector-bundle stabilization obstruction

Adding primitive reconstruction coordinates as an affine vector bundle does
not turn any normalized three-factor coefficient-boundary source into
affine space.  In degree one, the motivic defect survives multiplication by
an affine fiber.  In higher degree, every equivariant linear Cox extension
retains the free residual `mu_d` action, independently of the weights placed
on the new coordinates.

Thus the remaining Cox construction must be a genuinely nonlinear affine
modification, not a vector-bundle suspension of the existing source.

## 1. Degree-one sources

Let `Y_ell` be a normalized source associated with a nonzero linear
coefficient functional.  Over `C`, the
[linear-hyperplane classification](LINEAR_HYPERPLANE_COX_CLASSIFICATION.md)
gives exactly three compactly supported Hodge--Deligne polynomials, with
`q=uv`:

\[
 E_c(Y_\ell)=
 q^3-2q^2+q,\qquad
 q^3+q-1,\qquad
 q^3-2q-2.                                          \tag{1}
\]

Let

\[
 \pi:E\longrightarrow Y_\ell
\]

be a Zariski-locally trivial affine-space bundle of rank `r`; algebraic
vector bundles and ordinary reconstruction-coordinate suspensions are
included.  Then

\[
 [E]=\mathbb L^r[Y_\ell],
\qquad
 E_c(E)=q^rE_c(Y_\ell).                              \tag{2}
\]

None of the three expressions in (2) equals `q^(r+3)`.  Hence

\[
\boxed{
 E_\mathbb C\not\simeq\mathbb A^{r+3}_\mathbb C.
}                                                    \tag{3}
\]

The conclusion is independent of whether the vector bundle is trivial.

## 2. Higher-degree coefficient sources

Let `F` be homogeneous of degree `d>1`, and let

\[
 Y_F=\{r_{13}=r_{23}=F(P)=1\}.
\]

The [coefficient-hypersurface obstruction](COEFFICIENT_HYPERSURFACE_COX_OBSTRUCTION.md)
constructs the free residual action

\[
 \zeta\cdot(L_1,L_2,L_3)
 =
 (\zeta L_1,\zeta L_2,\zeta^{-1}L_3),
\qquad \zeta\in\mu_d.                                \tag{4}
\]

Now let `E -> Y_F` be a rank-`r` vector bundle with a linearization of this
`mu_d` action.  The fiber representation may contain arbitrary characters:

\[
 (z_1,\ldots,z_r)
 \longmapsto
 (\zeta^{w_1}z_1,\ldots,\zeta^{w_r}z_r).             \tag{5}
\]

The induced action on the total space is still free.  A fixed total-space
point would project to a fixed point of (4), but no such base point exists.
This remains true for every choice of weights `(w_1,...,w_r)`.

If the total space were `A^(r+3)`, its finite etale quotient by `mu_d`
would force

\[
 1=\chi_c(\mathbb A^{r+3})
  =d\,\chi_c(E/\mu_d),                               \tag{6}
\]

which is impossible for `d>1`.  Therefore

\[
\boxed{
 E_\mathbb C\not\simeq\mathbb A^{r+3}_\mathbb C.
}                                                    \tag{7}
\]

## 3. Consequence for primitive Cox coordinates

A separated Cox-ledger suspension, a compressed reconstruction coordinate,
or several coordinates transforming through boundary characters all have
source total space of the form covered above as long as the new coordinates
form an affine vector bundle over the normalized source.

The obstruction is not a shortage of dimensions:

> No finite number of vector-bundle Cox directions, with arbitrary linear
> boundary weights, can straighten these sources to affine space.

To evade (3) and (7), the added coordinates must change the source rather
than merely fiber over it.  Possibilities not excluded here include:

1. an affine modification whose center meets the old boundary;
2. a nonlinear complete intersection mixing old and new coordinates;
3. replacing the coefficient-boundary level by a different primitive
   boundary-class presentation; or
4. a cover with different nonsymmetric monodromy.

The triple-root slice supplies an explicit test of the first alternative.
Its source is `G_m^2 times A^1` and fills to `A^3`, but the
[triple-root completion theorem](TRIPLE_ROOT_AFFINE_COMPLETION_OBSTRUCTION.md)
shows that minimal target pole-clearing introduces three additional affine
Jacobian divisors.

## 4. Reproduction

Run

```bash
.venv/bin/python scripts/verify_cox_vector_bundle_stabilization.py
```

The checker verifies the affine-bundle Hodge products for all three linear
orbit types and the weight-independent freeness and Euler divisibility for
higher-degree residual Cox actions.
