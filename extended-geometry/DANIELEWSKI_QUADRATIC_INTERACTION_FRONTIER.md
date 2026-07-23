# The quadratic Danielewski interaction frontier

The minimal coupled theorem gives two complementary obstructions for

\[
P=x(x^2+1),\qquad x=bc,\qquad w=b(x^2+1).
\]

Diagonal ambient-quadratic terms satisfy a universal coefficient lock,
while the two escape terms \(cx,xw\), by themselves, lead to an
inconsistent normalized coefficient system.  This note classifies every
quadratic interaction chart containing at most four extra directions.

## 1. The full ambient-quadratic space

Let

\[
\mathcal Q=
\operatorname{span}_k
\{1,a,c,x,w,a^2,c^2,x^2,w^2,
ac,ax,aw,cx,cw,xw\}.                               \tag{1}
\]

This is the image in the source ring of all ambient polynomials of degree
at most two in \((a,c,x,w)\).  The relation \(cw=P(x)\) is retained through
the displayed generator \(cw=x+x^3\).

Only \(a,c,w\) have source-linear terms:

\[
w=b+\text{terms of source degree five}.
\]

Therefore every triple in \(\mathcal Q\) with nonzero constant Jacobian has
an invertible coefficient matrix on \(a,c,w\).  A constant target-linear
change normalizes it to

\[
(A,C,W)=(a,c,w)+\text{quadratic terms},             \tag{2}
\]

and the desired determinant is \(-1\).

## 2. Common, diagonal, and escape directions

Separate the quadratic monomials into

\[
\begin{aligned}
M&=\{ac,ax,aw\},\\
D&=\{a^2,c^2,x^2,w^2,cw\},\\
E&=\{cx,xw\}.                                      \tag{3}
\end{aligned}
\]

The terms in \(M\) are allowed in every chart.  The preceding coefficient
locks exclude every chart supported in \(M\cup D\), while the normalized
five-equation certificate excludes the chart supported in \(M\cup E\).

It remains to understand interactions between \(D\) and \(E\).

## 3. Exact chart theorem

For \(S\subset D\cup E\), let

\[
\mathcal Q_S=
\operatorname{span}_k
\bigl(\{1,a,c,x,w\}\cup M\cup S\bigr).              \tag{4}
\]

**Theorem.**  If

\[
1\le |S|\le4,\qquad S\cap E\ne\varnothing,
\]

then no triple in \(\mathcal Q_S\) has nonzero constant source Jacobian.

There are exactly

\[
\sum_{j=1}^4
\left(\binom7j-\binom5j\right)
=2+11+25+30=68                                    \tag{5}
\]

such charts.

For each chart, introduce independent coefficients for the three
normalized targets (2), expand

\[
\det\frac{\partial(A,C,W)}{\partial(a,b,c)}+1,
\]

and set every source-monomial coefficient to zero.  The linear equations
are eliminated exactly over \(\mathbb Q\).  In every one of the 68 charts,
the reduced nonlinear coefficient ideal has Groebner basis

\[
\{1\}.                                             \tag{6}
\]

Thus every chart is empty over every characteristic-zero field.  This is
an exhaustive coefficient-ideal calculation, not a bounded-height search.

## 4. New support lower bound

Combine the theorem with the diagonal and escape-only obstructions.
Any constant-Jacobian triple in the full space \(\mathcal Q\), if one
exists, must use at least five members of \(D\cup E\), and at least one
must lie in \(E\).  Equivalently, it must contain either

\[
\boxed{
\begin{array}{ll}
\text{one escape term }cx\text{ or }xw
&\text{and at least four diagonal terms},\\[2mm]
\text{both escape terms }cx,xw
&\text{and at least three diagonal terms}.
\end{array}
}                                                   \tag{7}
\]

The three common mixed terms \(ac,ax,aw\) remain available throughout;
the lower bound concerns additional quadratic support.

This turns the full quadratic problem from 127 possible subsets of
\(D\cup E\) into the 28 subsets satisfying (7).  Several size-five
charts also reduce to the unit ideal in exploratory calculations, but
larger charts can require substantially heavier elimination and are not
included in the chart-by-chart claim.

The remaining 28 charts are now closed simultaneously by the
[full quadratic obstruction](DANIELEWSKI_FULL_QUADRATIC_OBSTRUCTION.md):
column-major coefficient order gives the exact rational Groebner basis
\(\{1\}\) for the maximal 27-variable, 89-equation ideal.

## 5. Reproduction

Run

```bash
.venv/bin/python scripts/verify_danielewski_quadratic_interaction_frontier.py
```

The checker builds the full generic normalized determinant once, specializes
it to all 68 charts, eliminates their linear equations, and verifies an
exact unit Groebner basis over \(\mathbb Q\) in every case.
