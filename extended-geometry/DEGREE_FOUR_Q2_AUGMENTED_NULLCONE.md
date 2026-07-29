# The quartic \(q_2\)-augmented nullcone frontier

## 1. Status

Let
\[
 V_4=\operatorname{End}(\operatorname{Sym}^4),\qquad
 {\cal B}_4=\mathbb Q[\mu_1,\mu_2,\ldots,q_2].
\]
The global question is whether the common zero set of the moments and
\(q_2\) is exactly the \(\operatorname{SL}_2\)-nullcone.  Equivalently,
is \(R_4\) integral over \({\cal B}_4\)?

This note does not settle that question.  It gives the first exact
normal-jet calculation on the generic nonzero-\(\operatorname{Sym}^2\)
branch.  At one deterministic synchronized point, four moment equations
give independent normal pivots.  After formal elimination, the quadratic
normal jet has dimension six and the cubic normal jet has dimension four.
Thus cubic order does not isolate the synchronized nullcone branch.
Quartic and quintic normal terms, and the boundary where the
\(\operatorname{Sym}^2\) component vanishes, remain open.

## 2. Why \(q_2\) is the correct first branch variable

Use the multiplicity-free decomposition
\[
 V_4=\operatorname{Sym}^0\oplus\operatorname{Sym}^2
 \oplus\operatorname{Sym}^4\oplus\operatorname{Sym}^6
 \oplus\operatorname{Sym}^8.
\tag{2.1}
\]
The invariant \(q_2\) is the discriminant pairing on the binary
quadratic component \(F_2\).  Therefore
\[
 q_2=0,\quad F_2\ne0
 \quad\Longrightarrow\quad
 F_2=L^2
\tag{2.2}
\]
over an algebraic closure.  An \(\operatorname{SL}_2\) change normalizes
this to the highest-weight vector \(E\).  This reduces the first global
branch to synchronization of the higher binary-form components with the
same root \(L\).

Let \(E,F\) be the raising and lowering matrices on
\(\operatorname{Sym}^4\).  In the component
\(\operatorname{Sym}^{2r}\), use the integral weight basis
\[
 B_{r,k}=\operatorname{primitive}
 \bigl(\operatorname{ad}(F)^k(E^r)\bigr),
 \qquad 0\leq k\leq2r.
\tag{2.3}
\]
The common-root nullcone condition is
\[
 L^{r+1}\mid F_{2r},
\tag{2.4}
\]
which in this normalization retains precisely \(k<r\).  Hence the
synchronized slice has
\[
 2+3+4=9
\tag{2.5}
\]
allowed higher-component coordinates.  The twelve coordinates with
\(k\geq r\) are normal to it.

## 3. The exact normal-jet calculation

Work modulo \(32003\) and choose the allowed coordinates
\[
 (2,3,4,5,6,7,8,9,10)
\tag{3.1}
\]
in the ordered \(B_{r,k}\) basis with \(k<r\).  This is a one-sided
point, so every moment vanishes there.

The Taylor expansions of
\(\mu_2,\ldots,\mu_{21}\) are computed in the twelve forbidden
coordinates.  The linear coefficient matrix of
\(\mu_2,\ldots,\mu_5\) has rank four.  One exact pivot choice is
\[
 (r,k)=(2,3),(2,4),(3,6),(4,8).
\tag{3.2}
\]
The formal implicit equations solve those four coordinates through
normal degree three in the remaining eight coordinates
\[
 (2,2),(3,3),(3,4),(3,5),
 (4,4),(4,5),(4,6),(4,7).
\tag{3.3}
\]

Global degree-order standard bases of the resulting finite jets give
\[
\begin{array}{c|c|c}
\text{normal jet}&\text{affine dimension}&
\text{standard-basis size}\\ \hline
2&6&4\\
3&4&85.
\end{array}
\tag{3.4}
\]
These are exact statements over \(\mathbb F_{32003}\).  A nonzero minor
also certifies the four linear pivots in characteristic zero.

Equation (3.4) is not a local-isolation theorem.  The cubic jet is still
positive-dimensional, and a finite jet proves formal isolation only
after a maximal-ideal-power containment suitable for Nakayama has been
certified.  No such containment occurs through degree three.

Using native Singular map composition, the quartic jet reaches the
eight-variable standard-basis stage, but both `std` and `slimgb` with
`dp` ordering hit their declared \(300\)-second timeouts.  These
timeouts are a computational frontier and supply no evidence about the
quartic jet dimension.  The retained artifact records the `slimgb` run.

## 4. Interpretation

The calculation rules out two overly optimistic approaches.

1. Merely adding \(q_2\) and checking a Jacobian cannot prove
   synchronization: only four of the twelve forbidden coordinates are
   visible to first order.
2. Moments through cubic normal degree do not isolate the generic
   synchronized branch.  Six normal dimensions survive the quadratic
   jet and four survive the cubic jet.

This does not show that \(q_2\) is insufficient.  The surviving jet cone
may be killed by quartic or quintic Taylor terms.  It also does not
produce a semistable residual point: a positive-dimensional truncated
cone need not lift to an actual moment-zero branch.

The next exact calculation should use coefficientwise formal
composition, or a native Singular map, to finish the quartic and quintic
jets after the four pivots.  The target is either:

- a maximal-ideal-power containment proving formal isolation at the
  tested synchronized point; or
- an explicit formal arc in the residual cone, which can then be tested
  for semistability and algebraized.

Independently, the boundary \(F_2=0\) must be split by the first nonzero
component \(F_4,F_6,F_8\).  Even a successful generic
\(F_2=L^2\) calculation cannot prove global integrality without those
boundary branches.

## 5. Reproduction

Run

```bash
.venv/bin/python scripts/research_degree_four_q2_augmented_nullcone.py \
  --prime 32003 --max-jet 4 --composition native \
  --ordering dp --timeout 300
```

The generated
[`degree_four_q2_augmented_nullcone_local.json`](../artifacts/generated-results/degree_four_q2_augmented_nullcone_local.json)
has SHA-256
`9363d71f0d35f87b05c2ac370c98a0902f3bc9ce85402b0259261493efb8328b`.
The artifact records every allowed and forbidden weight coordinate, the
four pivots, the formal-series term counts, the reduced moment term
counts when Python composition is selected, both completed
standard-basis dimensions, and the bounded quartic timeout.
