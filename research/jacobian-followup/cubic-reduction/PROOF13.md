# A cubic Keller collision in thirteen variables with 42 nonlinear terms

There is an explicit rational polynomial map \(K:\mathbb A^{13}\to\mathbb A^{13}\)
with ordinary degree three,

\[
K(0)=0,\qquad JK(0)=I,\qquad \det JK=1,
\]

and three distinct rational points in one fibre. Its nonlinear monomial
occurrence count is **42**. This is a historical alternate construction, dominated by the already
verified twelve-variable,40-term BCRDB1 map. The initial comparison only
with the17-variable branch missed that existing baseline; see the
[correction](BASELINE_CORRECTION.md). The current contribution is the
[36-term twelve-variable conjugate](PROOF12.md).

No new counterexample mechanism, novelty priority, minimum dimension or
support optimum is asserted. The degree-three collision threshold across
varying dimensions was already established by the classical quadratic
midpoint argument and cubic reduction. The foundational map is the one
retained in [KFRAME1](../sparse-nine/PROOF.md).

## Explicit formula

Order the variables as \(X=(x,y,z,a_1,\ldots,a_{10})\), and set \(B=a_7-x^2\).
The first three outputs are

\[
\begin{aligned}
K_1={}&x-\tfrac32x^2y+\tfrac12a_6a_7-\tfrac12a_6x^2-\tfrac12a_7xz,\\
K_2={}&y+3xz+18xy^2+3a_4(xa_2-a_5)-6xa_2\\
 &+3za_5a_7+9xya_5-6a_6(a_8-xy)+6xza_8-3Ba_{10},\\
K_3={}&z+4y^2+3xyz+xa_1(a_4-2)-a_4(a_3+3a_2-3y^2)\\
 &+3xa_2(a_6+3y)+a_3(a_7z+3xy)-3a_5(a_6-xz)-Ba_9.
\end{aligned}
\]

The remaining outputs are

\[
\begin{aligned}
C_1&=a_1-y^3,&C_2&=a_2-y^2,\\
C_3&=a_3-xa_1,&C_4&=a_4-x^2z-3xy,\\
C_5&=a_5-xa_2,&C_6&=a_6-xz,\\
C_7&=a_7-x^2,&C_8&=a_8-xy,\\
C_9&=a_9-za_3,&C_{10}&=a_{10}-za_5.
\end{aligned}
\]

The row counts are
\((4,12,15,1,1,1,2,1,1,1,1,1,1)\), summing to 42.

Three colliding points, in the displayed source order, are

\[
\begin{aligned}
p_0={}&(0,-2,-16,-8,4,0,0,0,0,0,0,0,0),\\
p_+={}&(1,-2,8,-8,4,-8,2,4,8,1,-2,-64,32),\\
p_-={}&(-1,1,5,1,1,-1,2,-1,-5,1,-1,-5,-5).
\end{aligned}
\]

Every one maps to \((0,-2,0,0,\ldots,0)\).

The 42-term count refers to these tangent-normalized coordinates. The
collision fibre is not the origin fibre. Affine normalization placing a
collision at zero and \(e_1\) preserves degree three but need not preserve
this term count. The separate normalized three-variable degree-seven
support interval remains **[10,48]**.

## Construction by reusing output factors

Start from the tangent-normalized foundational map

\[
R=(F_3/2,F_2,F_1),\qquad \det JR=1.
\]

In an elementary degree-reduction step, an existing output can supply a
factor. To add a missing source monomial \(b(X)\), append an identity
coordinate and precompose its source with \(v\mapsto v-b(X)\).
Its output is now \(v-b(X)\). For a different output \(R_i\), the target
operation

\[
R_i\longmapsto R_i-c\,R_j\cdot\bigl(v-b(X)\bigr),\quad i\ne j,
\]

is elementary: \(R_j\) and the new carrier are unchanged. It preserves
the constant Jacobian and transports the full collision. Two existing
factors can be reused without introducing a variable. A square needs
only one new factor.

This differs from the earlier bounded splitting experiment, which
introduced fresh carriers for both factors. The first reused original
coordinate cancels two high-degree terms at once. The retained discovery
path has ten carrier extensions, thirteen product target steps and four
linear target steps. After multiplying by the inverse constant derivative
at zero, it gives a cubic map \(H\) in thirteen variables with 47 nonlinear
occurrences.

All 27 steps, including every coefficient and every transported source
coordinate, are retained in [the discovery witness](reuse-witness.json).
The [primary replay](construct13.py) and
[independent rational replay](verify13.py) agree on the entire map.

Finally set \(S=I-3E_{7,11}\), with one-based indices in the displayed
thirteen-variable order. Equivalently its source substitution is
\(a_4\mapsto a_4-3a_8\). Then

\[
K(X)=S^{-1}H(SX).
\]

This single conjugation reduces the nonlinear counts from
\((4,14,19,1,1,1,1,1,1,1,1,1,1)\) to the displayed 42-term vector.
It fixes the origin and preserves the identity derivative and determinant.

## Independent determinant and collision verification

The [final witness](cubic13-sparse-witness.json) contains every expanded
coefficient and all three points. The
[primary symbolic construction](construct_sparse13.py) uses SymPy 1.14.0.
The [independent checker](verify_sparse13.py) uses standard-library
rational sparse polynomial arithmetic.

The independent checker starts from the foundational degree-seven
coefficient map, verifies its full three-by-three determinant and three
collisions, and replays every carrier, target product and linear step.
It computes the inverse derivative by rational Gaussian elimination,
checks the resulting 47-term map, then replays the final conjugation and
compares every coefficient of the 42-term map.

It also performs a fresh full polynomial determinant check on the final
thirteen-variable map. With the original three variables first, write

\[
JK=\begin{pmatrix}A&B\\C&D\end{pmatrix},\qquad D=I+N.
\]

The complete symbolic carrier block satisfies \(N^3=0\). The checker
verifies \(D^{-1}=I-N+N^2\) on both sides and computes

\[
\det(A-BD^{-1}C)=1.
\]

Thus the thirteen-variable determinant is identically one. This check
keeps every source and carrier variable symbolic. It independently
checks the identity derivative, the degrees, every nonlinear occurrence
and every coordinate of all three colliding points.

[The final validation receipt](sparse-validation13.json) pins the primary
data, the discovery path, both independent arithmetic helpers and the
checker source.

## Reproduction and limits of the search

~~~sh
.venv/bin/python research/jacobian-followup/cubic-reduction/construct13.py
.venv/bin/python research/jacobian-followup/cubic-reduction/verify13.py
.venv/bin/python research/jacobian-followup/cubic-reduction/construct_sparse13.py
.venv/bin/python research/jacobian-followup/cubic-reduction/verify_sparse13.py
.venv/bin/python research/jacobian-followup/cubic-reduction/check_integration13.py
~~~

The [reused-factor search](search_reuse.py) was run with beam widths 60
and 800. Both found a thirteen-variable cubic. The first run is preserved
under [reuse-initial](reuse-initial/reuse-witness.json).

A separate [order experiment](search_reuse_orders.py) tests six
permutations of the original variables and three term orders, with beam
width 800 and maximum dimension twelve. All eighteen completed bounded
runs fail to find a cubic witness. The
[order manifest](orders-manifest.json) pins their checkpoints. These
runs are not a lower-bound certificate: they use monomial carriers,
specified factor choices, beam pruning and a greedy zero-cost closure.

[Selected quadratic graph slices](quadratic-slice-repair13.json) also
fail to produce a smaller cubic witness by the tested quadratic carrier
target corrections. This is a limited discovery result.

The [linear descent](linear-descent13.json) checks rational cancellation
values on every one of the 156 elementary lines at each of two iterations,
finding the shear above. It stops at a rational elementary-line minimum
of 42; it does not search the entire affine orbit or all algebraic frames.
No result from these searches is needed as a premise of the explicit
existence theorem.

For current optimization use the [twelve-variable baseline](PROOF12.md).
The failed bounded searches through dimension twelve do not exclude cubic
maps there: such maps were already verified independently. The overall
simplicity goal is open.

