# NS0031: valid modular arithmetic, unresolved rational marking

<!-- status-consumer: EC-K3-NS0031-QQ-MARKING-OBSTRUCTION cf8720d3d8c90dac -->

The September 4 nonexistence conclusion is withdrawn. Its finite Clifford,
coset and Frobenius calculations remain valid, but the asserted containment
of the full marked period group in rational norm-one units is false.
**Existence of a K3 over QQ with the full rational NS0031 marking is UNKNOWN.**
The model-157 formal branch and F017 physical corridor retain their local and
geometric meanings; neither supplies a rational source.

The [preserved note and ledgers](../archive/repository-cleanup-2026-09-12/ns0031-period-group-review/REVIEW.json)
record the pre-correction claims. The original arithmetic producer and its
certificate are unchanged. Their old K3-obstruction status label is historical,
not a current arithmetic exclusion.

<!-- status-consumer: EC-K3-NS0031-PERIOD-GROUP-COUNTERWITNESS c318a7d819b28d32 -->

## Exact counter-witness to the group containment

For the literal transcendental lattice and negative root

\[
G=\begin{pmatrix}0&0&4\\0&74&1\\4&1&-2\end{pmatrix},
\qquad v=(0,0,1)^t,\qquad v^tGv=-2,
\]

the root reflection is

\[
R=I+v(Gv)^t=\begin{pmatrix}1&0&0\\0&1&0\\4&1&-1\end{pmatrix}.
\]

Direct multiplication gives `R^t G R=G`, `R^2=I`, and `det(R)=-1`.
Moreover `(R-I)G^-1=vv^t` is integral, so `R` acts trivially on `T^dual/T`.
It fixes the plane spanned by `(1,0,2)` and `(0,2,1)` pointwise. That plane
has Gram `[[8,4],[4,298]]`, with positive first minor and determinant 2368.
Thus `R` preserves the chosen positive-plane orientation and belongs to the
literal stable group `O^+(T)^*`.

This is the group relevant to fixing the polarizing NS lattice: see
[Dolgachev, Proposition 3.3 and its following period-quotient discussion](https://arxiv.org/pdf/alg-geom/9502005v2).
A modern formulation for primitive lattice polarizations, with the small ample
cone hypothesis, is [Bragg–Brakkee–Várilly-Alvarado, Theorem 5.6](https://arxiv.org/pdf/2510.11477v1).
The discriminant-kernel condition extends `R` by the identity on NS across the
unimodular K3 gluing. A full NS marking therefore does not justify discarding
this determinant-minus-one isometry. These are complex period statements;
the arithmetic descent of the corrected curve still requires justification.

In the original split Clifford embedding, multiplication by the central
volume element identifies `T_Q` with the trace-zero matrices

\[
B_0=\begin{pmatrix}0&8\\0&0\end{pmatrix},\quad
B_1=\begin{pmatrix}-74&2\\0&74\end{pmatrix},\quad
B_2=\begin{pmatrix}0&-2\\74&0\end{pmatrix}.
\]

They span the trace-zero algebra and satisfy
`B_i B_j+B_j B_i=148*G_ij*I`. The matrix

\[
A=\begin{pmatrix}0&-1\\37&0\end{pmatrix},\qquad\det(A)=37,
\]

acts by `A B_j A^-1=sum_i (-R)_ij B_i`. Since `R` and `-R` have the same
projective period action, `A` represents this stable period transformation.
Any other rational matrix with the same adjoint action is a scalar multiple
of `A`; its determinant is `37*c^2`, never 1 for rational `c`.
The matrix also normalizes the retained integral order, as verified by an
explicit integral involution on its four basis elements. It is a normalizer
action, not a rational norm-one unit.

Consequently the claimed containment in the projective norm-one group fails
already over the complex period domain. The proposed period-induced map to
`X_ns(4) x_{X(1)} X_0(37)` cannot be obtained from that containment.
This does not construct a rational NS0031 K3 or disprove its nonexistence by
some other argument. It identifies the precise failure in the recorded proof.

## Arithmetic result that survives

The original even Clifford order, after conjugation by `diag(4,1)`, has basis
`I`, `[[0,4],[0,0]]`, `[[4,0],[0,0]]`, `[[1,1],[37,0]]`.
Its norm-one group has congruences
`37|C`, `B=C/37 mod4`, `A-D=C/37 mod4` for the entries of a matrix.
The associated norm-one modular curve has index 304, elliptic counts `(0,4)`,
cusp widths `[4,4,148,148]` and genus 23.

[Vélu, Theorem 7 and the following models, pp.175–176](https://www.numdam.org/article/MSMF_1974__37__169_0.pdf)
gives the two noncuspidal rational `X_0(37)` points, with j-invariants
`-7*11^3` and `-7*137^3*2083^3`. Both displayed elliptic curves have good
reduction at 19 and trace -6. Their mod-4 trace/determinant pair `(2,3)` is
absent from the full nonsplit Cartan; twisting preserves containment because
that Cartan contains `+/-I`. Hence neither point lifts to the stated Cartan
fibre product. The fibre product has no noncuspidal rational point.

That valid modular-curve obstruction is conditional input for any application
that separately proves the required rational lift. It no longer closes the
NS0031 rational-marking question. [Elkies, Section 2](https://arxiv.org/pdf/0802.1301v1)
already distinguishes norm-one groups, normalizer quotients and arithmetic
twists; replacing one with another needs an explicit proof.

## Replay and remaining gate

The [counter-witness](../artifacts/generated-results/elkies-k3-ns0031-period-group-counterwitness-v1.json)
and [standard-library verifier](scripts/verify_ns0031_period_group.py) check
the reflection, discriminant action, fixed positive plane, Clifford
identification, adjoint action and order normalizer. Six corruption controls and a valid-witness check accompany the proof;
five additional checks protect its propagation and preserved history. From `research/`, run
`python3 elkies-k3/scripts/verify_ns0031_period_group.py`.

The separate [finite arithmetic replay](scripts/verify_ns0031_marking_arithmetic.py)
uses the [2,771-byte original-input projection](../artifacts/generated-results/elkies-k3-ns0031-replay-inputs-v1.json)
and the [original arithmetic certificate](../artifacts/generated-results/elkies-k3-ns0031-qq-marking-obstruction-v1.json).
`make verify-ns0031-arithmetic` runs that bounded replay and its controls;
`--check-source-projection` optionally compares the three original catalogues.
No mode reconstructs missing inputs. The original Sage checker had already
passed its finite arithmetic replay; repeating it cannot repair the period map.

Original certificate SHA-256:
`49fc6570bf5a6e9411ae617e5c0aac45d04795af02a64f878d158cfa437818ae`.
Original projection SHA-256:
`f3aaaec35ee8a635d35ea295361ee5fd5439bf809b231cc74206d38ca4aeb80e`.

The next proof gate is to determine the full stable projective period group,
its correct model over QQ and its rational non-CM points, or supply another
exact arithmetic obstruction. Neither an empty norm-one cover nor the known
formal local branch decides this. NS0031 is an unresolved research row;
no equation or foundry campaign is authorized by this correction.
