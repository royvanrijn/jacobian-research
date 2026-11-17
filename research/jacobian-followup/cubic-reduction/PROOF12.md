# A 36-term integral cubic Keller map in twelve variables

The already verified [de Bondt twelve-variable cubic](../de-bondt-2610.24889/README.md),
claim BCRDB1, has **40 nonlinear monomial occurrences** in its displayed
coordinates. Two rational elementary conjugations reduce that count to
**36**. A diagonal conjugation clears the remaining halves without
changing support. The final map and both colliding points are integral.

This is a sparsity improvement at the already known dimension twelve.
No new dimension record, minimum dimension, support optimum or priority
is claimed. [The baseline correction](BASELINE_CORRECTION.md) explains why
the retained 17- and 13-variable experiments are dominated and must not be
reported as current dimension improvements.

The input construction is from
[de Bondt, arXiv:2610.24889v1, section 4.6](https://arxiv.org/html/2610.24889v1),
whose explicit maps were already independently audited in this repository.
The original counterexample mechanism is Alpöge's.

## Explicit map

Use \(X=(x_1,\ldots,x_{12})\). Define the last nine outputs first:

\[
\begin{aligned}
C_4&=x_4+x_1^2x_3+3x_1x_2,&C_5&=x_5+x_2^3,\\
C_6&=x_6+x_1^2,&C_7&=x_7+x_1x_3,\\
C_8&=x_8+x_2^2,&C_9&=x_9+2x_1x_5-3x_2^2,\\
C_{10}&=x_{10}+x_1x_8,&C_{11}&=x_{11}+x_1x_9,\\
C_{12}&=x_{12}+x_1x_2+x_1x_{10}.
\end{aligned}
\]

Put \(d=3x_2-x_7\), \(e=x_2+x_{10}\). The first three outputs are

\[
\begin{aligned}
K_1={}&x_1-x_1^2d+x_6C_7,\\
K_2={}&x_2+3x_1x_3+12x_1ed+12x_4C_{10}
                 +12x_1x_8-12x_{12}C_7,\\
K_3={}&x_3+8x_2^2+6x_1x_2x_3+4x_1x_9d+8x_1x_5
                 +4x_4C_9-4x_{11}C_7.
\end{aligned}
\]

Set \(K=(K_1,K_2,K_3,C_4,\ldots,C_{12})\).
Expansion gives nonlinear row counts
\((4,10,10,2,1,1,1,1,2,1,1,2)\), totaling **36**.
All coefficients are integers and the maximum total degree is three.
Moreover \(K(0)=0\), \(JK(0)=I\), and \(\det JK=1\).

The distinct integer points

\[
\begin{aligned}
p&=(0,1,-8,0,-1,0,0,-1,3,0,0,0),\\
q&=(-1,1,4,-1,-1,-1,4,-1,1,-1,1,0)
\end{aligned}
\]

both map to \((0,1,0,\ldots,0)\).

The term count uses the displayed tangent-normalized coordinates. The
collision fibre is not the origin fibre. It is separate from the
normalized three-variable degree-seven support problem, whose current
bounds remain **[10,48]**.

## Exact conjugation from the existing twelve-variable map

Let \(F\) be the section 4.6 map produced by the existing
[verified source](../de-bondt-2610.24889/verify.py).
Let

\[
S_1=I+\tfrac32E_{9,8},\qquad
S_2=I+\tfrac32E_{11,10},\qquad
D=\operatorname{diag}(1,\ldots,1,\tfrac12,1,\tfrac12,1),
\]

where only the ninth and eleventh diagonal entries differ from one.
With \(U=S_1S_2D\),

\[
K(X)=U^{-1}F(UX).
\]

Equivalently, the only nonidentity source coordinates of \(U\) are

\[
(UX)_9=(x_9+3x_8)/2,\qquad
(UX)_{11}=(x_{11}+3x_{10})/2.
\]

The first two conjugations change the occurrence count
\(40\to39\to36\). The diagonal change preserves every occurrence and
makes all coefficients and both collision points integral.
The target is transformed by the inverse matrix, so the identity
derivative and determinant one are preserved.

This is a coordinate simplification of the established construction,
not a new family or a claim that these linear combinations were absent
from the source paper.

## Independent verification

- [Pinned twelve-variable input](de-bondt12-witness.json) retains all
  coefficients, points and the source checker hash.
- [Rational primary construction](construct_sparse12.py) reconstructs
  the existing source map, applies both shears and computes the full
  polynomial Schur determinant.
- [Integral primary construction](construct_integer12.py) applies \(D\).
- [Final integral witness](cubic12-integer-witness.json) retains every
  expanded coefficient and both integer points.
- [Independent checker](verify_integer12.py) uses only standard-library
  rational sparse arithmetic. Its
  [rational replay](verify_sparse12.py) verifies both the 40-term input
  determinant and the 36-term determinant, both origin derivatives,
  every conjugation coefficient, and the collision. It then verifies
  the complete integral witness and recomputes its determinant.
- [Final receipt](integer-validation12.json) pins the sources and
  intermediate receipt. [Scoped integration](integration12.json) checks
  metadata, hashes, links and corrected optimization status.

For each full determinant, the checker differentiates every output in
all twelve variables and forms the block decomposition

\[
JK=\begin{pmatrix}A&B\\C&I+N\end{pmatrix}.
\]

It verifies \(N^3=0\), both identities for
\((I+N)^{-1}=I-N+N^2\), and the full symbolic identity
\(\det(A-B(I-N+N^2)C)=1\). No carrier is specialized during this check.

From the repository root:

~~~sh
.venv/bin/python research/jacobian-followup/cubic-reduction/construct_sparse12.py
.venv/bin/python research/jacobian-followup/cubic-reduction/construct_integer12.py
.venv/bin/python research/jacobian-followup/cubic-reduction/verify_integer12.py
.venv/bin/python research/jacobian-followup/cubic-reduction/check_integration12.py
~~~

The independent checker needs no SymPy. The primary construction uses the
pinned SymPy 1.14.0 environment and the existing BCRDB1 source.

## Remaining optimization

The [elementary-line descent](linear-descent12.json) tests all rational
coefficient-cancellation values on 132 elementary conjugation lines in
each of three iterations. No further single rational shear lowers 36.
This is not an affine-orbit optimum. The later
[KCUB19 carrier-realization theorem](CARRIER_REALIZATION_LIMIT.md) proves36
sharp only in its specified47-coefficient carrier class and excludes every
one-axis nonlinear source shear with same-axis target correction. Joint
frames, other carrier supports/orders and coordinated nonlinear stages remain open.
The later [KCUB20 theorem](FULL_CUBIC_CARRIER_CLASSIFICATION.md) closes all135
weighted quadratic/cubic carrier atoms in the original order, with minimum36,
full nonlinear-output rank12 and minimum cubic-output rank5. Its standard
minimum-dimensional companion locus is linearly equivalent to the old tensor.

The existing [nine literal graph-deletion obstructions](../de-bondt-2610.24889/README.md)
remain relevant. A dimension improvement needs a nonliteral coordinate,
a source change, sequential nonlinear stages, or another construction.
The full goal remains open: neither minimum cubic dimension nor minimum
term count is proved.

