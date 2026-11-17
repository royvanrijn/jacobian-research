# Historical GVC proposal from the 4 September review

This is the original proposal, retained as development history. The current
result is `GVC2SC` in [MATH_STATUS.json](../MATH_STATUS.json), with its
[canonical proof](../extended-geometry/BINARY_GVC_FINITE_CERTIFICATE.md).
It was strengthened and incorporated into the GVC manuscript on 5 September.

## A concrete GVC consequence to extract

**Proposed corollary for review, derived from the current binary proof.**
For nonzero binary inputs with positive lowest operator order, all pure
powers vanish if and only if, after a linear coordinate change over a finite
algebraic extension, their actual supports admit a strict positive weight
separator:

\[
g:=\min_{\alpha\in\operatorname{supp}\lambda} w\cdot\alpha
 -\max_{\beta\in\operatorname{supp}P}w\cdot\beta>0,
\qquad w\in\mathbb Z_{>0}^2.
\]

This is an extraction of the claimed binary proof, not an independent
verification of that proof or a newly promoted registry theorem.

Here is the deduction, including the degree-equality case:

1. If `deg P < ord_min Lambda`, ordinary degree already separates.
2. If the degrees are equal, the proof of Hall localization still applies:
   two independent derivative directions see all `d = r` polynomial
   factors, so cannot form a deficient set of size at most `r`. The same
   multiplicity gap `e > t` follows. At weight `(1+epsilon,1)` the gap is
   immediately strict. This extension of the stated Hall lemma should be
   written explicitly when adding the corollary.
3. If `deg P > ord_min Lambda`, the envelope proof reaches a common
   threshold. Shifted-ray separation with zero shift makes its two equality
   faces disjoint. A sufficiently small rational perturbation of the
   positive weight separates those faces strictly while preserving all
   off-face inequalities. Clear denominators.
4. Conversely, a strict separator kills every monomial selection in
   `Lambda^m(P^m)` by weighted degree, without cancellation.

For any fixed nonzero multiplier `Q`, put

\[
K_Q=\max_{q\in\operatorname{supp}Q}w\cdot q.
\]

The same certificate proves the explicit cutoff

\[
\Lambda^m(QP^m)=0\qquad
\left(m\ge\left\lfloor K_Q/g\right\rfloor+1\right).
\]

Both supports and `Q` are measured in the chosen coordinates. This is not
a coefficient-independent uniform cutoff for all inputs.

For exact algebraic coefficients, the proof also suggests a terminating
decision procedure: factor the lowest binary symbol, choose one coordinate
frame for each distinct derivative direction, and test the finitely many
strict linear inequalities in a rational slope `s>1`. The Hall direction
must occur among these finitely many frames. The lower-degree and zero or
constant-symbol cases are handled separately. This could replace unbounded
moment testing with a finite certificate. No implementation or complexity
claim for arbitrary characteristic-zero coefficient fields is asserted.

For example,

\[
\Lambda=\partial_x^2+\partial_y^7,
\qquad P=xy^2+y^3,
\qquad w=(3,1)
\]

has operator minimum `6`, polynomial maximum `5`, and gap `1`. For
`Q=y^6`, the certificate gives `m>=7`. An independent symbolic spot check
verified pure powers through eight and the predicted mixed vanishing at
seven and eight. The support argument, not that spot check, supplies the
all-order implication.

