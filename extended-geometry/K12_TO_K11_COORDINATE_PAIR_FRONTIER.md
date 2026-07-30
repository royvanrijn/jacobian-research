# The \(K_{12}\)-to-\(K_{11}\) coordinate-pair frontier

## Status

This note records an **exact bounded obstruction**, not a lower bound for
Keller counterexamples.  Starting from the certified twelve-variable
degree-three Keller collision \(K\), it proves:

- no linear target coordinate has a raw degree-three graph restriction;
- source coordinates equal to one literal triangular component \(K_j\);
- restriction to the graph \(K_j=0\);
- one parallel target-shear stage in the other raw retained outputs;
- target polynomials of degree at most three.

The two closest literal graph restrictions are also excluded through target
degree four.  Parameterized target completion over nonliteral linear graph
coordinates, nonlinear source coordinates, and multi-stage target
automorphisms remain open.

The parent map and collision are the ones proved in
[the twelve-variable counterexample](../verified/TWELVE_VARIABLE_DEGREE_THREE_KELLER_COUNTEREXAMPLE.md).
This note does not duplicate that theorem.

The optimization target is the compressed score \(n+r\), where \(n\) is the
dimension of a degree-three map and \(r\) is its cubic-output rank.  The
present map has \(12+6=18\).  An eleven-variable descendant retaining rank at
most six would reach \(n+r\leq17\), hence an eighteen-variable
cubic-homogeneous parent and a homogeneous quartic HN consequence in
dimension 36.  This would improve the current `19/38` endpoints.  It would
not be the first known eleven-variable nonhomogeneous map; that separate
public route and the score comparison are audited in the
[dimension-38 cross-construction frontier](HVC38_CROSS_CONSTRUCTION_FRONTIER.md).

## The bounded completion problem

Write

\[
K=(K_1,\ldots,K_{12}),\qquad K_i=z_i+N_i(z).
\]

A literal component \(K_j\) is a source coordinate for this audit when
\(N_j\) is independent of \(z_j\).  Then \(K_j=0\) is the polynomial graph

\[
z_j=-N_j(z_1,\ldots,\widehat z_j,\ldots,z_{12}).
\]

Substitution produces eleven retained raw outputs \(L_k^{(j)}\), \(k\ne j\).
Let \(\pi_{>3}\) retain only monomials of ordinary degree greater than three.
For a bad component \(i\ne j\) and target-degree bound \(d\), define

\[
V_{j,i,d}=
\operatorname{span}_{\mathbf Q}
\left\{
\pi_{>3}\!\left(\prod_{k\ne i,j}(L_k^{(j)})^{a_k}\right):
1\leq\sum a_k\leq d
\right\}.                                                       \tag{1}
\]

A one-stage triangular target repair

\[
y_i\longmapsto y_i-P_i(y_1,\ldots,\widehat y_i,\ldots,
\widehat y_j,\ldots,y_{12}),\qquad \deg P_i\leq d,
\]

can lower the restricted component to degree at most three only if

\[
\pi_{>3}(L_i^{(j)})\in V_{j,i,d}.                               \tag{2}
\]

Thus failure of (2) for one bad component is an exact obstruction to the
whole parallel repair stage.  Notice that linear mixing with every other raw
retained output is included in (1).

## All linear graph coordinates

The literal choices are only points inside a larger exact family.  Let

\[
g_a(y)=\sum_{i=1}^{12}a_i y_i,\qquad
h_a(z)=g_a(K(z)).
\]

For a chosen pivot \(z_j\), the equation
\(h_a(z)=g_a(K(p))\) is a polynomial graph in \(z_j\) precisely when
\(a_j\ne0\) and the nonlinear part
\(\sum_i a_i(K_i-z_i)\) is independent of \(z_j\).  These conditions are
linear in \(a\).  Exact row reduction gives:

| source pivot \(z_j\) | constraint rank | solution-space dimension | normalized family exists | parameters after \(a_j=1\) |
|---:|---:|---:|:---:|---:|
| 1 | 11 | 1 | no | — |
| 2 | 8 | 4 | no | — |
| 3 | 8 | 4 | no | — |
| 4 | 2 | 10 | yes | 9 |
| 5 | 4 | 8 | yes | 7 |
| 6 | 3 | 9 | yes | 8 |
| 7 | 2 | 10 | yes | 9 |
| 8 | 3 | 9 | yes | 8 |
| 9 | 1 | 11 | yes | 10 |
| 10 | 1 | 11 | yes | 10 |
| 11 | 1 | 11 | yes | 10 |
| 12 | 1 | 11 | yes | 10 |

For each of the nine normalized families, substitute its general graph into
all eleven retained raw outputs and set every coefficient above degree three
to zero.  The resulting parameter ideal has reduced Gröbner basis
\(\{1\}\) over \(\mathbf Q\) in all nine cases.  Consequently:

> **Linear raw-graph obstruction.** No linear target coordinate \(g_a\)
> whose pullback is a polynomial graph coordinate produces a raw
> eleven-variable map of degree at most three on the collision slice.

This statement does not yet allow target completion depending on \(a\).  The
completion audit below is exhaustive at the nine literal points \(a=e_j\).

## Exhaustive literal-coordinate result

Exactly \(K_4,\ldots,K_{12}\) are triangular in their matching source
variables.  Every raw graph restriction has degree four or five:

| deleted \(K_j\) | graph degree | raw maximum degree | bad original components | one certified obstruction, rank \(r\to r+1\) |
|---:|---:|---:|---|---:|
| 4 | 3 | 5 | 3, 10 | \(i=3:\ 276\to277\) |
| 5 | 3 | 5 | 2, 3, 7, 9 | \(i=2:\ 277\to278\) |
| 6 | 3 | 5 | 2, 3, 9 | \(i=2:\ 276\to277\) |
| 7 | 2 | 4 | 3 | \(i=3:\ 275\to276\) |
| 8 | 2 | 4 | 2, 3, 4 | \(i=2:\ 277\to278\) |
| 9 | 2 | 4 | 2 | \(i=2:\ 275\to276\) |
| 10 | 2 | 4 | 3 | \(i=3:\ 275\to276\) |
| 11 | 2 | 4 | 1 | \(i=1:\ 275\to276\) |
| 12 | 2 | 4 | 1 | \(i=1:\ 275\to276\) |

Here the target-degree-three basis has

\[
\binom{10}{1}+\binom{11}{2}+\binom{12}{3}
=10+55+220=285
\]

monomials.  Zero high-degree columns and dependencies are retained in the
calculation; the table reports the actual column-space rank.

For the closest cases \(j=11,12\), the raw excess is a single quartic
monomial in component one.  Enlarging (1) through target degree four gives
1,000 target monomials.  The ten linear columns have zero high-degree part,
the remaining 990 columns are independent modulo \(1{,}000{,}003\), and the
required defect adds a 991st pivot in both cases:

\[
\operatorname{rank}V_{11,1,4}=990,\quad
\operatorname{rank}[V_{11,1,4}\mid\pi_{>3}(L_1^{(11)})]=991,
\]

\[
\operatorname{rank}V_{12,1,4}=990,\quad
\operatorname{rank}[V_{12,1,4}\mid\pi_{>3}(L_1^{(12)})]=991.      \tag{3}
\]

These are rational obstructions, not heuristic finite-field failures.  A
rank increase modulo the good prime exhibits a nonzero augmented minor.
The same integer numerator is therefore nonzero over \(\mathbf Q\), so the
rational linear system (2) is inconsistent.

## What this teaches the next search

The \(13\to12\) reduction succeeded because the deleted source coordinate
was an output component and a quadratic target square owned the entire
quartic defect.  At \(12\to11\), literal graph deletion still produces
sparse defects, but none belongs to the bounded algebra generated by the
other raw outputs in one shear stage.  General linear graph coordinates do
not remove the defect without completion.  The next search should therefore
enlarge the *ownership language*, not merely raise the target degree:

1. solve parameterized target-completion equations over the nine classified
   linear graph-coordinate families;
2. if those unit ideals persist after completion, search low-degree
   nonlinear \(g\) while solving simultaneously for its graph and
   high-degree cancellation;
3. independently, allow two ordered target shears and recompute the second
   stage after the first, rather than using only the raw-output span.

The first branch is the most economical backwards attack: its parameter
spaces are now exact and finite-dimensional, and it keeps the coordinate-pair
identity that made the previous reduction work. The quadratic graph
subfamilies have now been pushed further in the
[parameterized completion frontier](K12_PARAMETERIZED_COMPLETION_FRONTIER.md):
quadratic completion is excluded for all six families, and cubic completion
is excluded for all five single-defect families.

## Reproduction

Run

```bash
make verify-k12-coordinate-pair-frontier
```

The checker reconstructs \(K\) from the pinned thirteen-variable formulas,
classifies all linear graph-coordinate families, verifies their raw
degree-three unit ideals, enumerates all literal triangular deletions, builds
every bounded target monomial, and performs sparse modular column
elimination.  It writes
[`k12_coordinate_pair_frontier.json`](../artifacts/generated-results/k12_coordinate_pair_frontier.json).
