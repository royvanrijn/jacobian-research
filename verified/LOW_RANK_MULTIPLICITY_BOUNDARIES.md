# Low-rank boundaries of universal Keller-fiber multiplicity

The universal multiplicity theorem has two deliberately different low-rank
boundaries:

1. rank four over arbitrary characteristic-zero fields;
2. rank three even over number fields.

This note closes the two immediate construction paths.  It does **not** prove
that universal stable multiplicity is false in either case.  It proves that
the present mechanisms cannot settle those cases by the arguments used in
ranks at least four.

## 1. Rank four: the trace-chord quadric can be anisotropic

The weighted rank-four construction requires an isotropic vector of

\[
 \mathcal Q_A(\eta,e,u)
 =\operatorname{Tr}_{A/K}(\eta^2)-2e^2-4u^2
\]

on

\[
 A_0\oplus K^2,\qquad A_0=\ker(\operatorname{Tr}_{A/K}).
\]

Over number fields, local isotropy and Hasse--Minkowski force such a vector.
This cannot be extended formally to every characteristic-zero field.

Let

\[
 K=\mathbb Q((a))((b)),\qquad
 A=K(\sqrt a,\sqrt b).
\]

The independent Laurent parameters `a,b` give a biquadratic field, so `A`
is connected finite etale of rank four.  Every trace-zero element has the
form

\[
 \eta=x\sqrt a+y\sqrt b+z\sqrt{ab},
\]

and

\[
 \operatorname{Tr}_{A/K}(\eta^2)
 =4(ax^2+by^2+abz^2).
\]

After scaling by `1/2`, the trace-chord form is

\[
 q=\langle 2a,2b,2ab,-1,-2\rangle.                    \tag{1.1}
\]

View `K` first as `\mathbb Q((a))((b))`.  Its `b`-adic Springer
decomposition is

\[
 q=q_0\perp bq_1,
\]

with

\[
 q_0=\langle-1,-2,2a\rangle,\qquad
 q_1=\langle2,2a\rangle.                               \tag{1.2}
\]

Over `\mathbb Q((a))`, the `a`-adic residue pairs are

\[
\begin{array}{c|cc}
 &\text{even residue}&\text{odd residue}\\ \hline
q_0&\langle-1,-2\rangle&\langle2\rangle\\
q_1&\langle2\rangle&\langle2\rangle .
\end{array}                                            \tag{1.3}
\]

Every one-dimensional form is anisotropic, and
`\langle-1,-2\rangle` is anisotropic over `\mathbb Q` by the real ordering.
Springer's theorem first makes both `q_0,q_1` anisotropic over
`\mathbb Q((a))`, then makes `q` anisotropic over `K`.

Therefore:

> **Trace-chord obstruction.**  
> There is a characteristic-zero field `K` and a connected quartic finite
> etale `K`-algebra `A` whose weighted trace-chord quadric has no `K`-point.

This proves that the number-field hypothesis in the present rank-four proof
is substantive.  It does not exclude other quartic Keller mechanisms, so the
unrestricted statement `|\mathcal R_K(A)|=\infty` remains open for this
example.

## 2. Rank three: the three current mechanisms collapse

Let `K` have characteristic zero.  Inside the weighted family, the normalized
cubic conditions

\[
 H(0)=H'(0)=H(1)=0
\]

force

\[
 H(W)=cW^2(1-W),\qquad c\ne0.
\]

Source and target scalings remove `c`, so the weighted cubic locus has one
stable class.

For a degree-three quadratic gauge, shear away the quadratic seed
coefficient and divide by the linear coefficient.  The normalized seed is

\[
 h_P(S)=S+a_3PS^3,\qquad a_3\ne0.
\]

The source--target scaling weights `a_3` by `(-2,-1)`, a transitive action on
`\mathbb G_m`.  Thus every cubic quadratic gauge is in the single
foundational class.

Finally, the cancellation degree formula

\[
 N=r(m+1)+1
\]

has the unique positive solution `(m,r)=(1,1)` when `N=3`.  Its parameter
polynomial is `q-3`, and the resulting map is polynomially left--right
equivalent to the same foundational cubic.

The quadratic-gauge realization theorem realizes every rank-three finite
etale algebra.  Since all its cubic seeds are in the foundational orbit, the
same one stable class contains a representative having any prescribed cubic
algebra as a complete fiber.  Weighted or cancellation presentations add no
second class.

Hence:

> **Three-mechanism cubic collapse.**  
> For every rank-three finite etale `K`-algebra `A`, the union of the
> weighted, cancellation, and root-engineered quadratic-gauge mechanisms
> contributes exactly one element to `\mathcal R_K(A)`: the foundational
> cubic class.

This closes the direct cross-family path to cubic multiplicity.  It does not
classify all geometric-degree-three Keller maps.  The remaining question is
precisely whether a Keller map outside these three controlled mechanisms can
have the same complete cubic fiber without being stably equivalent to the
foundational map.

The active route to that classification is the
[cubic closure protocol](../cancellation/CUBIC_CLOSURE_ATTACKS.md).  Its
remaining load-bearing issue is intrinsic extraction and elimination of the
closed-point normalization defect, not another seed-normalization search.

## 3. Updated frontier

The exact frontier is now:

\[
\begin{array}{c|c}
\text{scope}&\text{status}\\ \hline
K\text{ a number field},\ N\ge4
  &|\mathcal R_K(A)|=\infty\\
\operatorname{char}K=0,\ N\ge5
  &|\mathcal R_K(A)|=\infty\\
\operatorname{char}K=0,\ N=4
  &\text{proved when the trace-chord quadric is isotropic;}\\
  &\text{that quadric can be anisotropic}\\
\operatorname{char}K=0,\ N=3
  &\text{all three current mechanisms give one class;}\\
  &\text{unrestricted multiplicity remains open.}
\end{array}
\]

## 4. Exact regression

Run

```bash
.venv/bin/python scripts/verify_low_rank_multiplicity_boundaries.py
```

The checker verifies the unique normalized cubic weighted seed, transitivity
of the cubic quadratic-gauge coefficient action, uniqueness of the
degree-three cancellation type, and the exact diagonal trace form and
Springer decomposition in the biquadratic quartic example.  The two uses of
Springer's theorem are the written anisotropy proof.
