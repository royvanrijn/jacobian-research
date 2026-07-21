# Compactification and dicritical divisors of the weighted maps

Let `H` be an admissible primitive of degree `n` over a characteristic-zero
field, normalized by

\[
H(0)=H(1)=0,\qquad H'(0)=0,\qquad H'(1)=-c\ne0.
\]

For the associated weighted map `G_H`, a target `(A,B,C)` has inverse pencil

\[
E(W)=H(W)-BCW+cAC^2.                              \tag{1}
\]

This note gives a compactification of the inverse correspondence and
classifies every divisor at infinity which dominates a divisor in the target.

The intrinsic finite-cover formulation is proved in the archived
[weighted marked-root derivation](../../archive/core-support/WEIGHTED_MARKED_ROOT_MODEL.md).  The
normalization of the marked-root incidence contains the affine source as its
regular-reconstruction open, and the divisors classified below are exactly
the polar divisors in its complement.  The graph model used here retains the
valuations needed for nonproperness and omission.

The explicit graph blow-up, local toric charts, discrepancies and intersection
graphs are retained in the
[geometry-support archive](../../archive/geometry-support/C16_BLOWUP_GEOMETRY.md).

## Resolved inverse graph

Embed the affine source in `P^3` and take the closure of the graph of `G_H` in
`P^3 x A^3_(A,B,C)`.  Normalize that closure and then normalize the graph
of the rational root function

\[
W=(1+xy)(1+axy+bx^2z).
\]

Denote the result by `bar(X)_H`.  The projection

\[
\bar G_H:\bar X_H\longrightarrow\mathbb A^3_{A,B,C}
\]

is proper and restricts to `G_H` on the affine source.  Equation (1), with
`W` homogenized in `P^1`, is the finite root cover controlling this normalized
graph.  Since the coefficient of `W^n` is the nonzero leading coefficient of
`H`, the homogenized pencil has no root at `W=infinity`; all boundary divisors
come from poles of reconstruction, not from a missing root chart.

A boundary prime of `bar(X)_H` is called dicritical when its image is a
divisor in the affine target.

## Theorem

Write the primitive-root profile as

\[
H(W)=hW^{m_0}(W-1)\prod_{j=1}^q(W-\rho_j)^{\mu_j},
\]

where the `rho_j` are distinct and different from `0,1`.  On the normalized
inverse graph, the dicritical divisors are exactly the following.

1. **The discriminant divisor `D_Delta`.**  It is present for every seed.  Its
   normalization has coordinates `(r,C)` and maps by
   \[
   B={H'(r)\over C},\qquad
   A={rH'(r)-H(r)\over cC^2},\qquad C=C.          \tag{2}
   \]
   Its image is the saturated pulled-back discriminant `V(Q_H)`.

2. **The extra-root divisors `D_(rho_j)`.**  For every nonzero primitive root
   `rho_j` of multiplicity `mu_j`, one divisor dominates the target plane
   `C=0`.  After the ramified chart
   \[
   C=\delta^{\mu_j},\qquad W=\rho_j+\kappa\delta,
   \]
   its generic boundary equation is
   \[
   \kappa^{\mu_j}={B\rho_j\over h_j},             \tag{3}
   \]
   where `h_j` is the first nonzero local coefficient of `H` at `rho_j`.
   Intrinsically `D_(rho_j) -> V(C)` has residue degree one and ramification
   index `mu_j`.  The ramified chart displays the `mu_j` nearby escaping
   sheets; it does not make them distinct prime divisors.

3. **The zero-cluster divisor `D_0`.**  It occurs exactly when `m_0>=3`.
   With
   \[
   C=\delta^{m_0-1},\qquad W=\kappa\delta,
   \]
   its generic boundary equation is
   \[
   \kappa^{m_0-1}={B\over h_0}.                  \tag{4}
   \]
   Hence `D_0 -> V(C)` has residue degree one and ramification index `m_0-1`.
   It accounts for the `m_0-1` zero-cluster sheets which do not extend to the
   finite `gamma=0` chart.

There are no other dicritical divisors.  The distinguished simple root `1`
extends to the finite `x=0` chart.  When `m_0=2`, both generic zero-cluster
sheets extend to the finite `gamma=0` chart, so the zero cluster creates no
dicritical divisor.

## Proof by valuations

On `C!=0`, reconstruction gives

\[
\gamma=-{E'(W)\over c},\qquad x=-{cC\over E'(W)}. \tag{5}
\]

The only possible reconstruction pole is `E'(W)=0`.  Together with `E(W)=0`
this gives

\[
BC=H'(r),\qquad cAC^2=rH'(r)-H(r),
\]

which is (2).  This proves the universal discriminant divisor and also shows
that its image is the normalization image of the saturated discriminant.

At `C=0`, equation (1) specializes to `H(W)=0`, so every remaining boundary
valuation is centered at a primitive root.  Near a nonzero root `rho` of
multiplicity `mu`, write

\[
H(W)=h_\rho(W-\rho)^\mu+O((W-\rho)^{\mu+1}).
\]

Substitution of `C=delta^mu`, `W=rho+kappa delta` into (1) gives at first
order

\[
\delta^\mu(h_\rho\kappa^\mu-B\rho),
\]

proving (3).  Moreover

\[
E'(W)\asymp\delta^{\mu-1},\quad
x\asymp\delta,\quad
\gamma\asymp\delta^{\mu-1}.
\]

For `rho!=1`, the invariant `W=(1+xy)gamma` cannot approach `rho` through the
bounded `x=0` chart, where `W` approaches `1`.  These branches therefore lie
at source infinity.

Near zero, write `H(W)=h_0W^{m_0}+O(W^{m_0+1})`.  The scale
`C=delta^(m_0-1)`, `W=kappa delta` gives

\[
\delta^{m_0}\kappa(h_0\kappa^{m_0-1}-B),
\]

and proves (4) after removing the finite `kappa=0` branch.  For `m_0>=3`,
`gamma` has order `m_0-1` while `W` has order one, so `1+xy=W/gamma`
diverges.  These are boundary branches.  For `m_0=2` the quotient remains
finite and gives the two finite `gamma=0` points already computed directly.

The specialization `H(W)=0` exhausts all centers over `C=0`, while (5)
exhausts all poles over `C!=0`; hence the list is complete.

## Nonproperness and omission

The images of the dicritical divisors give the nonproperness set:

\[
S_{G_H}=V(Q_H)\cup
\begin{cases}
V(C),&m_0\ge3\text{ or }q>0,\\
\varnothing,&m_0=2\text{ and }q=0.
\end{cases}                                      \tag{6}
\]

Here repetitions change the ramification indices along the maps to `V(C)`,
not their residue degrees or images.  Formula (6) recovers the exceptional
minimal cubic and the general boundary theorem without relying only on a
discriminant saturation.

Omission occurs on `C!=0` precisely when every root of (1) is multiple.  In
the compactification this means that the entire length-`n` inverse fiber is
supported on `D_Delta`; no affine point remains.  If the multiplicity
partition is `(m_1,...,m_l)`, the normalization of `D_Delta` has `l` points
over the target, carrying boundary lengths `m_1,...,m_l`.  Thus the
contact-partition stratification is exactly the stratification by how the
dicritical discriminant divisor exhausts the inverse fiber.

## Executable valuation audit

Run

```bash
python scripts/verify_dicritical_divisors.py
```

The script checks (2), the Kummer leading equations (3)--(4) for
multiplicities through five or six, and the finite nature of the distinguished
root-one chart.  The Kummer exponents are checked as ramification indices.
The graph blow-up and all exceptional valuations are audited separately by
`verify_c16_blowup_geometry.py`.  These finite multiplicity checks are
regressions for the uniform valuation argument above.
