# The parameterized \(z_8\) cubic-completion obstruction

## Status

This note proves an **exact bounded theorem** over \(\mathbf Q\). Together
with the
[parameterized completion frontier](K12_PARAMETERIZED_COMPLETION_FRONTIER.md),
it closes one-stage target completion through degree three for every linear
coordinate of the twelve-variable map whose pullback has a quadratic graph.
It is not a dimension-eleven lower bound.

## The remaining family

The only multi-defect quadratic graph family is

\[
g_a(y)=a_0y_7+y_8+a_1y_9+a_2y_{10}+a_3y_{11}+a_4y_{12}.         \tag{1}
\]

Restricting to \(g_a(K)=g_a(K(p))\) and solving for \(z_8\) leaves original
components \(2,3,4\) above degree three. To obstruct the whole parallel
target stage, it suffices to obstruct completion of component 2 using all
other ten raw retained outputs, including the two other bad outputs.

There are

\[
\binom{10}{1}+\binom{11}{2}+\binom{12}{3}
=10+55+220=285
\]

target monomials through degree three. Eight linear columns have no
high-degree part. The other 277 columns occupy 54,977 source-monomial rows.

## Sparse minor-first construction

Expanding the full matrix over
\(\mathbf Q(a_0,\ldots,a_4)\) obscures its sparse structure. Instead, the
checker represents every coefficient by two exponent blocks:

\[
(\text{source exponent in }\mathbf N^{11},
 \text{parameter exponent in }\mathbf N^5).
\]

Products and cancellations are performed in this sparse representation.
Evaluation modulo \(1{,}000{,}003\) selects candidate pivot rows, but no
modular rank claim is used as the final certificate. Only the resulting
\(277\times277\) column minors and \(278\times278\) augmented minors are
reconstructed exactly over \(\mathbf Q[a_0,\ldots,a_4]\).

Three column determinants suffice. Up to nonzero rational constants they
are

\[
\begin{aligned}
\Delta_0={}&
(a_2a_4+1)(2a_2a_4+1)^2(4a_2a_4+1)\\
&\qquad\cdot
(18a_0a_3a_4^2+7a_2a_4+1),\\
\Delta_1={}&a_2^{22}a_4^{159},\\
\Delta_2={}&a_0^2a_3^{29}a_4^{160}.                            \tag{2}
\end{aligned}
\]

For every selected row set, the exact augmented determinant is

\[
\Delta_\nu^{\mathrm{aug}}=\frac97\Delta_\nu.                    \tag{3}
\]

The opens in (2) cover the entire parameter space. If \(a_4=0\), then
\(\Delta_0\ne0\). If \(a_4\ne0\) and \(a_2\ne0\), then
\(\Delta_1\ne0\). Finally, if \(a_4\ne0\), \(a_2=0\), and
\(\Delta_0=0\), then

\[
18a_0a_3a_4^2+1=0,
\]

so \(a_0a_3\ne0\) and \(\Delta_2\ne0\). Equivalently, exact Gröbner
reduction gives

\[
(\Delta_0,\Delta_1,\Delta_2)=(1).                               \tag{4}
\]

At every parameter point, one column determinant is nonzero and its
augmented determinant is also nonzero by (3). Thus the target defect is
outside the full completion span.

## Theorem and consequence

> **Parameterized \(z_8\) cubic-completion obstruction.** The family (1)
> admits no one-stage target completion of degree at most three that lowers
> every retained component to degree at most three.

Combining this theorem with the five constant-minor families proves that
none of the six quadratic graph-coordinate families admits a one-stage
cubic target completion.

The new insight is methodological as well as negative: a very large
parameterized completion problem can have a tiny exact certificate if rows
are selected modularly and only the selected minors are reconstructed.
This “minor-first” compiler is the appropriate next tool for:

- target degree four on the same graph families;
- the larger linear-coordinate families with cubic graph correction;
- ordered two-stage completions, where full fraction-field expansion is
  even less practical.

## Reproduction

Run

```bash
make verify-k12-z8-cubic-completion
```

The generated exact record is
[`k12_z8_cubic_completion_frontier.json`](../artifacts/generated-results/k12_z8_cubic_completion_frontier.json).
