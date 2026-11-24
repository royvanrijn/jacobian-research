# Strong Factorial for identical monomials on disjoint variables

## Scope and attribution

This internal research note proves an all-order, all-exponent ordinary
Strong Factorial family. Its complex-root separation input is the explicit
coefficient extension in Remark 3.4 of Jianxi Mao, Lijie Wang and Sainan
Zheng, *Interlacing on the Unit Circle via Coefficientwise Reciprocals*,
[arXiv:2609.34551v1](https://arxiv.org/abs/2609.34551v1), dated
28 September 2026. Theorem 1.1, Remark 3.4, and their actual proof in
Sections 2–3 were inspected; this is a recent preprint, not a claim of
external review of this application or of priority for the result.

Let `K` be any characteristic-zero field. Write
`Lambda(X^alpha)=prod alpha_i!` for the ordinary factorial functional.
Choose a nonconstant monomial template

\[
 M(X)=\prod_{i=1}^t X_i^{r_i},\qquad r_i\in\mathbb Z_{\ge0},
 \qquad D=\sum_i r_i\ge1.
\]

Let `M(X)` and `M(Y)` use disjoint variable sets and let `c,d in K`.

**Theorem.** If `F=c M(X)+d M(Y)` is nonzero, then for every `n>=1`,

\[
                  (\Lambda(F^n),\Lambda(F^{n+1}))\ne(0,0).
 \tag{1}
\]

Thus every such two-term polynomial satisfies SFC at every starting
order, with no bound on the exponents, number of variables, or coefficients.
If one coefficient is zero, every moment is nonzero. A nonzero monomial
factor in further disjoint variables can also be included, since it
multiplies each moment by a nonzero factorial scalar.

This is a structured two-term family. The template equality and disjoint
variables are hypotheses. The theorem does not close arbitrary binomials
or unrestricted SFC, and does not assert SFC for arbitrary powers of `F`:
adjacent nonvanishing does not by itself control the sampled orders `mn`.

## 1. The exact moment polynomials

Set

\[
 a_k=\frac{\prod_i(r_i k)!}{k!},\qquad
 f_n(z)=\sum_{k=0}^n a_k a_{n-k}z^k,
 \qquad
 q_n(z)=\sum_{k=0}^n\frac{z^k}{a_k a_{n-k}}.
 \tag{2}
\]

The disjointness of the two variable blocks gives, for `c!=0`,

\[
                  \Lambda(F^n)=n!c^n f_n(d/c).
 \tag{3}
\]

Both `f_n` and `q_n` have strictly positive rational, palindromic
coefficients. The coefficientwise reciprocal of `q_n` is exactly `f_n`.
The normalization required by the imported theorem is

\[
 \frac{[z^k]q_n(z)}{\binom nk}
     =\frac{\delta_k\delta_{n-k}}{n!},
 \qquad
 \delta_k=\frac{k!}{a_k}
          =\frac{(k!)^2}{\prod_i(r_i k)!}.
 \tag{4}
\]

The factorial `k!` in the definition of `delta_k` is essential. Mere
log-convexity of `a_k` is a weaker condition and is insufficient for
invoking this theorem.

## 2. Uniform normalized log-concavity

Suppose `D>=2`. Extend the ratio of consecutive `delta` values to positive
real arguments by

\[
 R(x)=\frac{(x+1)^2}
 {\displaystyle\prod_i\prod_{j=1}^{r_i}(r_i x+j)},\qquad x\ge0.
 \tag{5}
\]

Zero exponents contribute empty products. Every denominator factor is
positive, and

\[
 \frac{R'(x)}{R(x)}
 =\frac{2}{x+1}
   -\sum_i\sum_{j=1}^{r_i}\frac{r_i}{r_i x+j}
 \le\frac{2-D}{x+1}\le0,
 \tag{6}
\]

because `j<=r_i` makes each summand at least `1/(x+1)`.
As `delta_(k+1)/delta_k=R(k)`, the positive sequence `delta_k` is
log-concave for all `k`. Therefore for each `n` the sequence on the
right side of (4) is log-concave: it is the pointwise product of a
log-concave sequence and its reversal, multiplied by a positive constant.
This argument includes equality cases, such as the template `XY`.
No real-rootedness of `q_n` is required or asserted.

For `D=1`, one has `a_k=1` and

\[
              f_n(z)=1+z+\cdots+z^n.
 \tag{7}
\]

The identity `f_(n+1)-z f_n=1` proves consecutive coprimality directly.
In this exceptional case `delta_k=k!` is not log-concave, so it is
handled separately instead of forcing the hypothesis of (6).

## 3. The inspected interlacing input

The version used from Mao–Wang–Zheng, Remark 3.4, says that two positive
palindromic polynomials of consecutive degrees, whose coefficients divided
by binomial coefficients are log-concave, have coefficientwise reciprocals
with simple unit-circle roots that strictly interlace. This is a paraphrase.
The following records the proof mechanism and the conditions checked
during review.
For a normalized reciprocal polynomial of degree `N`, let

\[
 G(z)=1+\sum_{k=1}^{\lfloor N/2\rfloor}\gamma_k z^k
\]

be its lower half, with its middle coefficient halved for even `N`.
The coefficient hypothesis gives
`gamma_1<=1/N` and `gamma_1>=2 gamma_2>=...>0`.
Eneström–Kakeya places the reversed half inside the unit disk.
Fejér's inequality `sum_(k=1)^j sin(kt)/k>0` for `0<t<pi`, followed
by summation by parts, gives `Im G(e^(it))>0`. The estimate
`sum (2k-1) gamma_k<1` gives

\[
                 0<\arg G(e^{it})<t/2\qquad(0<t<\pi).
 \tag{8}
\]

The boundary Blaschke phase is `N t-2 arg G(e^(it))`.
For consecutive degrees its difference lies strictly between `0` and
`2t`; conjugation and the opposite signs at `-1` exclude coincidence
on the rest of the circle. The strictly increasing Blaschke phases
then give strict interlacing. Degree one is checked directly.
This uses normalized log-concavity, not a conjectural interlacing
principle for arbitrary log-convex series.

The inspected source was the actual PDF, whose first page identifies
`arXiv:2609.34551v1 [math.CO] 28 Sep 2026`:

- PDF SHA-256: `bf4731acc6af6348c4c5b306eaffff7d0be660f67945b2234458e742ac816e13`.
- Extracted text SHA-256: `9c4f1c16648a3e193c6a0024669e72ec85a8e3c543adcad0598c7b4ded00b5b3`.
- Extracted Sections 2–3 SHA-256: `de91eb83ce4ec7a2e5bd7ff5be59a3c7023c19f7f42da9eb2dcd78f6d405a445`.

The last hash is computed from the exact UTF-8 text beginning at the
first `2. Basic properties` and ending just before the following
`4. Applications` in the retained `pdftotext -layout` extraction.
Extraction hashes identify the reviewed text; another PDF extraction
program can produce different whitespace without changing the theorem.

## 4. Completion over every characteristic-zero field

For `D>=2`, (4)–(6) verify every hypothesis of Remark 3.4 for `q_n`
and `q_(n+1)`. Thus `f_n` and `f_(n+1)` strictly interlace and have
no common complex root. For `D=1`, (7) proves the same coprimality.
Since these are rational polynomials,

\[
                  \gcd_{\mathbb Q[z]}(f_n,f_{n+1})=1.
 \tag{9}
\]

Their rational polynomial Bézout identity remains valid after scalar
extension to every characteristic-zero field. Applying it to `z=d/c`
in (3) proves (1), whether the coefficients are real, complex, algebraic,
or transcendental. There is no real-coefficient assumption on `F`.

## 5. Boundaries and exact verification

The example `a_k=1` is positive and log-convex, but its degree-two
reciprocal polynomial has normalized coefficients `(1,1/2,1)`, which
are not log-concave. This rejects the inference from arbitrary
log-convexity to the hypothesis of Remark 3.4; it is not a
counterexample to consecutive interlacing, since (7) already treats it.

Likewise, root-of-unity filtered moments of powers or products require
the separate hypotheses in `SFC_CYCLIC_MONOMIAL_BLOCKS.md`. They are
not consequences of (9) alone. Shared variables or unequal templates
change (2), so those cases cannot silently inherit this argument.

Run from the repository root:

```text
python research/factorial/image-rigidity/verify_sfc_identical_monomial_pair.py
```

The checker uses exact rational arithmetic. It reconstructs factorial
moments directly, checks (3)–(6) on retained instances, verifies
consecutive coprimality by exact polynomial Euclidean division, tests
the degree-one identity, and rejects the stated incorrect hypothesis
shortcut and deliberately altered coefficient data. These finite
controls support the identities; the all-order proof is Sections 1–4.
An optional `--source-pdf` verifies the retained external PDF hash.
No floating-point root computation or empirical extrapolation is used.
