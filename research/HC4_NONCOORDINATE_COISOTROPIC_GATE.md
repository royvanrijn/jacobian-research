# Smallest non-coordinate coisotropic scalar gate for `HC_4`

## Status

This note records a bounded exact experiment on the first reciprocal
mixed-line coisotropic graphs beyond the coordinate family in `HC4MCP1`.
It finds no survivor of the graph-specialized scalar Schur gate.  No full
four-variable determinant is needed.

For ordered \(i\ne j\), set

\[
 K=q_i+\rho p_j,\qquad L=q_j+\rho p_i,\qquad
 H=\tau K L^2.                                      \tag{0.1}
\]

The search exhausts

\[
 \rho,\tau\in\{-2,-1,1,2\},\qquad
 \lambda\in\{-1,1\},\qquad
 \mu\in\{-1,0,1\}.                                 \tag{0.2}
\]

There are 96 Hamiltonian charts, 128 affine scalar pivots, and 768
specialized scalar-gate trials.  Every trial has an exact nonconstancy
witness modulo \(1000003\), so none can be constant over
characteristic zero.  A post-gate audit also gives a nonconstant parent
Hessian determinant in all 96 charts.  This is the finite-box result
`HC4MCP8`; it is not a parameter-uniform theorem for arbitrary nonzero
\(\rho,\tau,\lambda\), and it does not classify general coisotropic
embeddings.

The later
[Meng--Yang graph obstruction](HC4_MENG_YANG_GRAPH_OBSTRUCTIONS.md) is
parameter-uniform for every affine four-plane and every single graph
\(y_3=R(x_1,x_2,y_1,y_2)\) of degree at most three in the v2
five-variable potential.  Its two-slope continuation also excludes every
graph through degree four over an arbitrary characteristic-zero field.  It
then makes the degree-five plane equation uniquely solvable for the normal
jet, forces a vertical quintic trace, and excludes one complete sparse trace
family by first-transverse coefficients.  The general degree-five graph and
non-graph generating families remain open; none of these parameter-uniform
graph results turns this finite Hamiltonian-chart calculation into a
classification of nonlinear generating families.

Replay the search with

```bash
.venv/bin/python \
  scripts/search_hc4_noncoordinate_coisotropic_scalar_gate.py \
  --output \
  artifacts/generated-results/hc4_noncoordinate_coisotropic_scalar_gate.json
```

The generated JSON has
`sha256:a6a806d649c1af0c1ea2c26e01937817832de502daacda97a1a015faca472eb2`
under the repository `.python-version` and `requirements.txt` locks.

## 1. Coisotropic parameterization

With the source--dual Poisson convention used in the mixed-pivot search,

\[
 \{K,L\}=\rho-\rho=0.                               \tag{1.1}
\]

Hence the time-one flow of (0.1) is polynomial: \(K,L\), and therefore
the Hamiltonian velocity, are constant along the flow.  Moreover
\(\{K,p_i\}=1\) and \(\{L,p_i\}=0\).  Up to the pullback/pushforward sign
convention, the coordinate coisotropic \(p_i=0\) is sent to

\[
 p_i+\tau L^2=0.                                    \tag{1.2}
\]

Both invariant linear directions in (0.1) are mixed when \(\rho\ne0\).
Thus (1.2) is the smallest reciprocal mixed-line nonlinear graph family;
the unit slopes recover the corresponding two-mixed-line cases already
present, but not separately classified as coisotropic graphs, in
`HC4MCP1`.

## 2. The weaker gate

For every coordinate \(t\) in which the pulled-back six-variable
potential is affine, write

\[
 \Phi=tA(w)+B(w).
\]

The script tests the exact necessary scalar condition from `SDX1`,

\[
 \left.
 \det\operatorname{Hess}_w(B+sA)
 \right|_{s=\mu+\lambda A(w)}
 \in K.                                             \tag{2.1}
\]

This is weaker than demanding the whole pencil be identically singular
and is tested before a full descended Hessian determinant.  A pair of
unequal evaluations in \(\mathbf F_{1000003}\) proves that the specialized
integer polynomial in (2.1) is nonconstant over \(\mathbb Q\).

The affine-pivot census by ordered pair is

\[
\begin{array}{c|rrrrrr}
(i,j)&(0,1)&(1,0)&(0,2)&(2,0)&(1,2)&(2,1)\\ \hline
\text{pivots}&32&32&16&16&16&16.
\end{array}
\]

No modular survivor remains.  Since the family already fails (2.1), the
scalar route has no candidate to promote.  The separate parent audit gives
two unequal modular Hessian-determinant values in every chart, so none is
also a constant-Hessian rechart of the foundational parent.  Collision
transfer and complete descended determinants are therefore not evaluated.

## 3. Relation to the exceptional sextic atlas

The independent symmetric-sextic calculation in
[`HC4_EXCEPTIONAL_SCHUR_LOCUS.md`](HC4_EXCEPTIONAL_SCHUR_LOCUS.md)
classifies the full reduced projective quartic incidence:

\[
\{(0,0)\}\times
\mathbb P\langle x^4,y^4,z^4\rangle
\;\sqcup\;
\left\{\left(
\frac15,\frac1{10},[(x^2+y^2+z^2)^2]
\right)\right\}.                                   \tag{3.1}
\]

Thus there are no exceptional parameter curves on that surface.  The
Fermat module in (3.1) is excluded at lower determinant layers by
`HC4QF1`, and the radial line by `HC4QSE2`; collision equations are never
reached.  The coisotropic screen above produces no new four-variable
candidate in its finite parameter box.
