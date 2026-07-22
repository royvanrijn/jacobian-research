# Factorized BCW circuit exploration and the 22-variable quotient

## Result

The compact foundational circuit is the correct next search representation,
but its first exact stable lift does not by itself beat the monomial trace.
More importantly, the requested structural audit finds a two-dimensional
constant Jacobian kernel in the existing 24-variable homogeneous map.  The
map therefore descends to an exact 22-variable cubic-homogeneous Keller
collision.  The current route bound improves from

\[
  \neg\operatorname{GMC}(48)
  \quad\hbox{to}\quad
  \boxed{\neg\operatorname{GMC}(44)}.
\]

The generator is
[`verify_constant_kernel_bcw_22_route.py`](../scripts/verify_constant_kernel_bcw_22_route.py),
and it writes
[`constant_kernel_bcw_22_counterexample.json`](../artifacts/generated-results/constant_kernel_bcw_22_counterexample.json).

## 1. The compact foundational DAG

Put

\[
 U=1+xy,\qquad R=xz,\qquad H=Uz,\qquad V=yU.
\]

The canonical determinant-minus-two foundational map is exactly

\[
\begin{aligned}
 F_1&=U^2H+yV+3V^2,\\
 F_2&=-2y+3U^2R+9UV-6V,\\
 F_3&=5x-3xU-x^2R.
\end{aligned}
\]

Thus the normalized identity-linear map is
`K=(F_3/2,F_2,F_1)`.  This representation exposes common polynomial nodes
which the current exponent-vector registry cannot express.

There is one important stable-equivalence caveat.  Exposing a factor `a` as
an output `A=w+a` and subtracting a circuit expression in `A` creates mixed
terms containing `w`.  Formal replacement of `a` by an independent gate
variable is therefore not a BCW certificate.

For example, simultaneously expose

\[
 U,R,H,V,x^2,y^2
\]

using six source shears and apply the three multi-term target shears

\[
\begin{aligned}
 T_0&\mapsto T_0+\tfrac12 A_{x^2}A_R,\\
 T_1&\mapsto T_1-3A_U^2A_R-9A_UA_V,\\
 T_2&\mapsto T_2-A_U^2A_H-A_{y^2}A_U-3A_V^2.
\end{aligned}
\]

These cancel all displayed high-degree circuit terms at once, and every
operation is an exact determinant-one stable shear.  Nevertheless, expansion
of the mixed residual has maximum degree six.  It has twelve terms above
degree three and provisional cubic-output rank five in dimension nine.  A
width-24 run of the old single-monomial cleanup from this seed finishes with
twenty introduced variables and cubic-output rank eight, much worse than the
certified objective `13+7=20`.  This is not evidence against the circuit
approach; it shows that polynomial factors must remain first-class throughout
the residual cleanup.

The next search state should therefore contain:

- a hash-consed polynomial expression DAG, not only expanded exponent tuples;
- exposed-output nodes with component ownership and lifetime;
- target rewrites which may cancel a sum of products in one shear;
- exact residual polynomials, simplified after each triangular source or
  target gate;
- both the degree-reduction dimension and the exact cubic component rank.

The terminal objective is `n+rank(C)`, equivalently `4+s+rank(C)` when the
three-dimensional source and homogenizing coordinate are separated.  A useful
beam key is

\[
(\max\deg,\;s+\operatorname{rank}C,\;s,\;\operatorname{rank}C,
 \text{high-degree residual cost},-\text{live reusable nodes}).
\]

## 2. Constant kernel of the 24-variable map

For the homogeneous part `H_24` of the certified rank-compressed map, exact
rational coefficient comparison gives

\[
 \dim_{\mathbb Q}\bigcap_X\ker JH_{24}(X)=2.
\]

A basis is

\[
 k_1=-3e_6+e_{10},\qquad
 k_2=-\frac2{11}e_{12}-\frac3{11}e_{13}+e_{15}.
\]

Although all 24 coordinate variables occur in the monomial support of
`H_24`, its essential input span has dimension only 22: it depends on 22
linear forms and is constant on cosets of `K=span(k_1,k_2)`.

Choose a rational matrix `B` with `ker(B)=K` and a section `C` with

\[
 BC=I_{22}.
\]

The generator verifies the exact Gorni--Zampieri-type identities

\[
 H_{24}=H_{24}CB,
 \qquad
 f_{22}=B F_{24} C
       =q+B H_{24}(Cq).
\]

If `S=(C\mid K)`, then in `S`-coordinates the 24-variable map is

\[
 S^{-1}F_{24}S(q,r)
   =\bigl(f_{22}(q),\ r+g(q)\bigr).
\]

This is the requested triangular identity factor: a two-dimensional identity
fiber translated by a polynomial depending only on the quotient coordinate.
Consequently

\[
 \det Df_{22}=\det DF_{24}=1.
\]

Projecting the three stored collision points by `B` leaves three distinct
points and gives a literal common image.  Hence `f_22` is noninvertible.

## 3. Diagnostics after quotienting

For the new 22-variable homogeneous part `H_22`:

- the constant kernel of `JH_22` is zero;
- all 22 quotient variables occur, so there is no further linear input
  quotient of this type;
- 21 output components of `H_22` are nonzero and their coefficient-row span
  has dimension 21;
- the coordinate-dependency graph has two strongly connected components:
  the block `{0,...,20}` and the singleton `{21}`;
- coordinate 21 is the homogenizing variable: its output is the identity,
  while the other components depend on it.

Thus the two constant directions are genuinely removable.  The remaining
singleton is a triangular homogenizing coordinate, not a removable constant
kernel direction; setting it to one recovers a nonhomogeneous map rather than
a smaller cubic-homogeneous witness.

## 4. Search consequence

The certified baseline for future circuit search is now 22, not 24.  A new
degree-reduced candidate `K=X+Q+C` should be scored twice:

1. form the rank-compressed homogenization of dimension `n+rank(C)+1`;
2. quotient the constant kernel of its homogeneous `JH` to a fixed point;

and compare the final quotient dimension.  Constant-kernel removal can change
the winner even when two traces have the same `s+rank(C)`.

