# Low-rank boundaries of universal Keller-fiber multiplicity

The original diagonal/weighted multiplicity mechanisms have two deliberately
different low-rank boundaries:

1. rank four over arbitrary characteristic-zero fields;
2. rank three even over number fields.

This note closes the two immediate construction paths.  The later
[power-shifted quartic gauge theorem](UNIVERSAL_QUARTIC_GAUGE_MULTIPLICITY.md)
settles the unrestricted rank-four question positively over every
characteristic-zero field.  Thus Section 1 is now an exact obstruction to
the trace-chord mechanism, not an open boundary for quartic multiplicity.
The later
[fiber-invisible cubic gauge theorem](UNIVERSAL_CUBIC_GAUGE_MULTIPLICITY.md)
also settles unrestricted rank-three multiplicity positively.  Section 2
remains the exact collapse theorem for the three minimal mechanisms.

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

This proves that the number-field hypothesis in the weighted rank-four proof
is substantive.  It does not exclude other quartic Keller mechanisms.  In
fact, the power-shifted gauge theorem gives
`|\mathcal R_K(A)|=\infty` for this example.

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

This closes the direct cross-family path to cubic multiplicity, but it does
not classify all geometric-degree-three Keller maps.  The fiber-invisible
gauge theorem constructs infinitely many maps outside these three controlled
mechanisms, disproving unrestricted cubic-class uniqueness.

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
  &|\mathcal R_K(A)|=\infty\text{ by power-shifted gauges;}\\
  &\text{the weighted trace-chord submechanism can be obstructed}\\
\operatorname{char}K=0,\ N=3
  &|\mathcal R_K(A)|=\infty\text{ by fiber-invisible cubic lifts;}\\
  &\text{the three minimal mechanisms still give one class.}
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
