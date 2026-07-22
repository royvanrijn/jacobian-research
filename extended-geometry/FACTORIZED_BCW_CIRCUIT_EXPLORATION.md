# Factorized BCW circuit exploration and essential quotients

## Result

The compact foundational circuit is the correct next search representation,
but its first exact stable lift does not by itself beat the monomial trace.
More importantly, the requested structural audit finds a two-dimensional
constant Jacobian kernel in the existing 24-variable homogeneous map.  The
map therefore descends to an exact 22-variable cubic-homogeneous Keller
collision.  A subsequent essential-dimension search finds a different
17-dimensional trace whose rank-compressed 24-dimensional homogenization has
a three-dimensional constant kernel.  Its exact 21-dimensional quotient
improves the route bound from

\[
  \neg\operatorname{GMC}(44)
  \quad\hbox{to}\quad
  \boxed{\neg\operatorname{GMC}(42)}.
\]

The generator is
[`verify_constant_kernel_bcw_22_route.py`](../scripts/verify_constant_kernel_bcw_22_route.py),
and it writes
[`constant_kernel_bcw_22_counterexample.json`](../artifacts/generated-results/constant_kernel_bcw_22_counterexample.json).

### Constant-kernel quotient proposition

Over a characteristic-zero field, let `F=id+H:A^n->A^n` be a polynomial map
and let the columns of `K` span a constant subspace contained in `ker DH(x)`
for every `x`.  Choose a quotient
`B:k^n->k^(n-d)` with `ker B=im K` and a section `C` with `BC=I`.  Then

\[
 H=HCB,
 \qquad f=BFC=id+BHC.
\]

Writing `S=(C|K)`, there is a polynomial `g` such that

\[
 S^{-1}FS(q,r)=(f(q),r+g(q)).
\]

Indeed, `DH K=0` makes `H` constant on the affine `K`-cosets, which gives
`H=HCB`; applying `S^-1` gives the displayed triangular form.  Consequently
`det DF=det Df`.  Cubic homogeneity descends, and any collision whose
`B`-images remain distinct descends to `f`.  This linear-section/quotient
construction is closely related to Gorni--Zampieri pairing, without requiring
an identification with every form of that correspondence.

The independent standard-library replay
[`audit_constant_kernel_bcw_22_independent.py`](../scripts/audit_constant_kernel_bcw_22_independent.py)
recomputes the kernel from the 24-dimensional source artifact and verifies
every identity in this proposition for the displayed quotient.

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
certified objective `14+6=20`.  This is not evidence against the circuit
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

There is also no more general collision-preserving **linear** quotient.  The
66 coefficient matrices of `JH_22` have only three common invariant row
spaces: zero, the homogenizing covector `span(e_21^T)`, and the full space.
The middle space is constant on all three collision points.  The exact audit
proves this by showing that the induced coefficient algebra modulo
`span(e_21^T)` is the full matrix algebra `M_21(Q)`, while a common invariant
hyperplane is excluded by the zero common kernel of the 65 nilpotent
coefficient matrices.  See
[`audit_bcw_22_linear_quotients.py`](../scripts/audit_bcw_22_linear_quotients.py).

## 4. Search consequence

The certified baseline for future circuit search is now 21, not 22.  Every new
degree-reduced candidate `K=X+Q+C` is scored by the mandatory pipeline in the
[constant-kernel quotient theorem](../verified/CONSTANT_KERNEL_QUOTIENT.md):

1. form the rank-compressed homogenization of dimension `n+rank(C)+1`;
2. quotient the constant kernel of its homogeneous `JH` to a fixed point;
3. verify separation of the projected collision and compute invariant-row
   module diagnostics;
4. compare the final quotient dimension.

Constant-kernel removal can change
the winner even when two traces have the same `s+rank(C)`.

For the present witness, further gains must therefore occur before or beyond
linear quotienting: a better circuit-level stable trace, a trace satisfying
the two-parameter quadratic--cubic determinant identity and hence avoiding
rank-compressed doubling, or a genuinely nonlinear triangular decomposition.

The first automatic essential-dimension experiment found a
17-dimensional quadratic--cubic
trace, cubic-output rank six, 24-dimensional homogenization, a
three-dimensional constant kernel, and final quotient dimension 21.  See the
[verified theorem and experiment note](../verified/CONSTANT_KERNEL_QUOTIENT.md).
The frozen generator and dependency-free replay now certify `not GMC(42)`.
The exact invariant-row-module audit also proves that the only proper module
is the homogenizing covector, which is constant on the collision, so this new
map has no further collision-preserving linear quotient.
No minimality claim is made.
