# Parameterized target completion on \(K_{12}\)

## Status and objective

This note proves an **exact bounded obstruction** over \(\mathbf Q\). It is
not a lower bound for Keller counterexamples.

The preceding
[coordinate-pair frontier](K12_TO_K11_COORDINATE_PAIR_FRONTIER.md)
classifies every linear target coordinate \(g_a=\sum a_i y_i\) whose
pullback \(g_a(K)\) is a polynomial graph coordinate. The present
calculation addresses the next question: can parameters in \(g_a\) make a
target completion possible even though every literal coordinate is
obstructed?

An eleven-variable descendant with cubic-output rank at most six would lower
the compressed score from \(12+6=18\) to at most \(11+6=17\). It would
therefore give an eighteen-variable cubic-homogeneous counterexample and a
dimension-36 homogeneous quartic HN consequence.

## The six quadratic graph families

Write \(K_i=z_i+N_i\). A normalized linear target coordinate

\[
g_a(y)=\sum_{i=1}^{12}a_i y_i,\qquad a_j=1,
\]

has a quadratic graph with pivot \(z_j\) when
\(\sum_i a_iN_i\) is independent of \(z_j\) and has no cubic part. Exact
coefficient elimination shows that such families exist precisely for

\[
j=7,8,9,10,11,12.
\]

They have four parameters for \(j=7\) and five parameters otherwise:

| pivot | normalized target coordinate |
|---:|---|
| 7 | \(y_7+a_0y_8+a_1y_9+a_2y_{11}+a_3y_{12}\) |
| 8 | \(a_0y_7+y_8+a_1y_9+a_2y_{10}+a_3y_{11}+a_4y_{12}\) |
| 9 | \(a_0y_7+a_1y_8+y_9+a_2y_{10}+a_3y_{11}+a_4y_{12}\) |
| 10 | \(a_0y_7+a_1y_8+a_2y_9+y_{10}+a_3y_{11}+a_4y_{12}\) |
| 11 | \(a_0y_7+a_1y_8+a_2y_9+a_3y_{10}+y_{11}+a_4y_{12}\) |
| 12 | \(a_0y_7+a_1y_8+a_2y_9+a_3y_{10}+a_4y_{11}+y_{12}\) |

Restriction to the collision slice \(g_a(K)=g_a(K(p))\) produces bad
original components

\[
\{3\},\ \{2,3,4\},\ \{2\},\ \{3\},\ \{1\},\ \{1\},
\]

respectively. All excess terms are tested, not merely the top homogeneous
piece.

## Fixed-minor completion certificate

Choose one bad retained component \(L_i(a)\). For target-degree bound \(d\),
let \(M_{j,d}(a)\) have as columns the coefficients above source degree three
in every nonconstant target monomial of degree at most \(d\) in the other
ten raw retained outputs. Let \(b_j(a)\) be the high-degree coefficient
vector of \(L_i(a)\).

For each recorded parameter point, exact row reduction selects a square
full-column minor with determinant

\[
\Delta_\nu(a)=\det M_{j,d}(a)[R_\nu,*].
\]

Using one additional row, the augmented determinant satisfies

\[
\det[M_{j,d}(a)\mid b_j(a)][R'_\nu,*]
=c_\nu\Delta_\nu(a),\qquad c_\nu\in\mathbf Q^\times.             \tag{1}
\]

Finally, the checker verifies

\[
(\Delta_\nu(a):\nu)=\mathbf Q[a_0,\ldots,a_s]                    \tag{2}
\]

by a Gröbner basis equal to \(\{1\}\). At every parameter point, at least
one \(\Delta_\nu\) is nonzero. Equation (1) then gives full column rank and
strictly larger augmented rank, so the completion system has no solution
over \(\overline{\mathbf Q}\), hence none over \(\mathbf Q\).

For five families one constant minor suffices at target degree two:

| pivot | matrix | column determinant | augmented/column ratio |
|---:|---:|---:|---:|
| 7 | \(479\times55\) | \(-3^{25}\) | \(1\) |
| 9 | \(396\times55\) | \(-3^7/2\) | \(-9\) |
| 10 | \(323\times55\) | \(3^{25}/2\) | \(-1\) |
| 11 | \(537\times55\) | \(-3^{24}\) | \(-1/2\) |
| 12 | \(514\times55\) | \(-2^4 3^{25}\) | \(-1/2\) |

For pivot \(z_8\), two other bad raw outputs add two nonzero linear columns,
so the selected system has 57 columns. Four determinant opens cover its
five-dimensional parameter space. Their determinant ideal has Gröbner basis
\(\{1\}\), and every augmented determinant is \(-6\) times its column
determinant. This closes the exceptional hypersurfaces of the first minor
rather than discarding them as nongeneric.

## Exact result

> **Parameterized quadratic-completion obstruction.** None of the six
> quadratic graph-coordinate families with pivots
> \(z_7,\ldots,z_{12}\) admits a one-stage target completion of degree at
> most two that lowers every retained component to degree at most three.

The five single-defect families \(z_7,z_9,z_{10},z_{11},z_{12}\) admit a
stronger certificate. Adding all 220 cubic target monomials gives 275
nonzero high-degree columns. In every family a constant \(275\times275\)
minor and a constant \(276\times276\) augmented minor have nonzero ratio.
Therefore:

> **Parameterized cubic-completion obstruction.** None of the five
> single-defect quadratic graph families admits a one-stage target
> completion of degree at most three.

The multi-defect \(z_8\) case is completed separately by the
[sparse minor-first cubic certificate](K12_Z8_CUBIC_COMPLETION.md).
Together, the two results exclude cubic target completion for all six
quadratic graph families. They still do not cover:

- target degree at least four;
- larger linear graph families with cubic graph correction;
- nonlinear target coordinates or ordered multi-stage automorphisms.

The \(z_8\) calculation shows how to continue: keep coefficients as sparse
polynomials in graph parameters, use modular evaluation only to select
stable minors and exceptional strata, and reconstruct the selected
determinants exactly over \(\mathbf Q\).

## Reproduction

Run

```bash
make verify-k12-parameterized-completion
```

The checker reconstructs \(K\), derives every graph family, builds the full
high-degree completion matrices, recomputes all fixed determinants and
unit-ideal covers, and writes
[`k12_parameterized_completion_frontier.json`](../artifacts/generated-results/k12_parameterized_completion_frontier.json).
