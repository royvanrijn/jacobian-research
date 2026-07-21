# Boundary saturation with repeated primitive roots

This note removes the squarefree-extra-root hypothesis from the boundary
theorem. Combined with `OMITTED_VALUE_CLASSIFICATION.md`, it completes the
image and nonproperness schema for arbitrary admissible characteristic-zero
weighted seeds.

## Root profile and saturation exponent

Write an admissible primitive as

\[
H(W)=h,W^{m_0}(W-1)
\prod_{j=1}^r(W-\rho_j)^{\mu_j},
\]

where `m_0>=2`, every `rho_j` is distinct from `0,1`, and `mu_j>=1`. For the
inverse pencil

\[
E(W)=H(W)-BCW+cAC^2,
\]

put

\[
e=m_0+\sum_{j=1}^r(\mu_j-1).
\]

Then the pulled-back discriminant has exact `C`-adic order `e`:

\[
\operatorname{Disc}_W(E)=C^eQ_H(A,B,C),
\qquad C\nmid Q_H.
\]

This is the required saturation factor; no higher power of `C` may be
discarded.

### Why each cluster contributes this amount

Near `W=0`, the Newton polygon of

\[
h_0W^{m_0}-BCW+cAC^2
\]

shows that the zero cluster contributes `m_0` to the discriminant valuation.
For `m_0=2`, both roots have scale `C`. For `m_0>=3`, one branch has scale
`C`, the other `m_0-1` branches have scale `C^{1/(m_0-1)}`, and the squared
pairwise differences again total `m_0`.

Near a nonzero root `rho_j` of multiplicity `mu_j`, the generic local equation
starts as

\[
h_j(W-\rho_j)^{\mu_j}-B\rho_jC.
\]

Its branches have scale `C^{1/\mu_j}`, so their squared pairwise differences
contribute `mu_j-1`. Differences between distinct root clusters tend to
nonzero constants and contribute no additional `C`-power. Summing the cluster
valuations proves the formula for `e`.

## The saturated boundary trace

Let `h_0` be the coefficient of `W^{m_0}` in `H`, and put

\[
M=\sum_j(\mu_j-1).
\]

Up to a nonzero scalar, the restriction of the saturated equation to `C=0`
is

\[
Q_H(A,B,0)=
\begin{cases}
B^M(B^2-4ch_0A),&m_0=2,\\
B^{m_0+M},&m_0\ge3.
\end{cases}
\]

Thus repeated extra roots add exactly the predicted powers of `B`. They do
not erase the double-zero conic, and clearing `C^e` loses no intersection
stratum.

## Direct fibers and escaping repeated branches

The finite source over `C=0` is unchanged by nonzero roots of `H`:

- `x=0` maps bijectively to the target plane and corresponds to `W=1`;
- `gamma=0` has `W=0` and gives the same `3/1` fiber table when `m_0=2`, or
  the same `2/1` table when `m_0>=3`.

No nonzero root `rho_j` supplies another finite chart. If `mu_j>1`, then along
a generic transverse approach the nearby branches satisfy

\[
W-\rho_j=O(C^{1/\mu_j}),
\qquad
\gamma=O(C^{(\mu_j-1)/\mu_j}),
\qquad
x=O(C^{1/\mu_j}).
\]

Since `W` tends to `rho_j!=1` while every bounded `x=0` source has `W->1`,
these branches escape in the remaining source coordinates. Repetition changes
their rates, not their status as boundary branches.

The direct `C=0` fiber has fewer points than the generic inverse degree.
Étaleness prevents distinct bounded sheets from coalescing. Therefore every
target on `C=0` is a nonproper value whenever an extra root is present or
`m_0>=3`.

## Complete arbitrary-seed schema

Let `Omega_H` be the exact finite set returned by
`classify_omitted_values(H,W)`. Then

\[
G_H(\mathbb C^3)=\mathbb C^3\setminus
\{(A,B,C):C\ne0,\ (BC,cAC^2)\in\Omega_H\}.
\]

Apart from the minimal cubic primitive, the nonproperness set is

\[
S_{G_H}=V(C)\cup V(Q_H),
\]

where `Q_H` is obtained by dividing the pulled-back discriminant by the exact
power `C^e` above. This statement covers simple, multiple, split, and nonsplit
extra factors without losing a boundary stratum.

The implementation `boundary_saturation_profile(H,W)` computes `m_0`, the
extra repetition excess, `e`, and the boundary trace directly over the
coefficient field. An irreducible factor of degree `d` occurring with
multiplicity `mu` contributes `d(mu-1)`, correctly counting its roots over the
algebraic closure.

Run:

```bash
.venv/bin/python scripts/verify_repeated_root_boundary.py
```

The exact audit covers single and multiple repeated extra roots, both zero
multiplicity regimes, inverse degrees through eight, direct affine charts, and
a nonsplit repeated quadratic factor.
