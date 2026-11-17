# One new cube cannot shorten the saved 89-form identity

A stronger [fixed-span result](SPAN_OPTIMALITY.md) now also closes arbitrary
simultaneous changes whose scalar cubes stay in the saved span.

The exact89-form padded companion identity used by [KCUB5](POWER89.md)
cannot be shortened by deleting at least two of its selected forms,
adding **one arbitrary linear cube**, and re-solving the coefficients
of every retained cube. This holds over every characteristic-zero
extension field. The new form is unrestricted; it need not belong to
any previous finite search pool.

This is a replacement obstruction for one particular identity. It does
not prove89 minimal for the tensor, much less for all cubic-linear
Keller maps. A shorter identity can introduce two or more new forms;
sequences with temporary support growth can also escape this result.
The fixed tensor's certified global shared-cube interval remains
[34,89](FIXED_TENSOR_BOUND.md).

## Reduction to two canceled coefficient columns

Write the retained identity as

\[
P=\sum_{a=1}^{89}\beta_a c_a,\qquad c_a=\ell_a^3.
\]

The89 scalar cubics are linearly independent, and all coefficient
columns \(\beta_a\in\mathbb Q^{18}\) are nonzero. Pure deletion is
therefore impossible. Suppose a shorter representation retains a
subset of these forms and introduces one cube \(m^3\) with nonzero
coefficient column \(b\). Selecting an output where \(b\ne0\) shows
that \(m^3\) lies in \(S=\operatorname{span}(c_1,\ldots,c_{89})\).
Write its unique coordinates as \(m^3=\sum_a a_a c_a\).

Every deleted column must satisfy \(\beta_j=b a_j\). In particular,
at least two deleted columns \(j,k\) must be proportional and both
coordinates nonzero. There are exactly134 proportional column pairs.
For each, write \(\beta_k=\rho\beta_j\). A compression requires
\(a_k=\rho a_j\) with \(a_j\ne0\). Over an algebraic closure we may
scale \(m\) by a cube root to normalize \(a_j=1\). It suffices to
exclude the equations

\[
m^3\in S,\qquad a_j(m)=1,\qquad a_k(m)=\rho.
\]

## Exhausting the possible source supports

Join source coordinates i and j when some selected form has both
coefficients nonzero. If both coefficients of m are nonzero, the
monomial \(z_i^2z_j\) occurs in \(m^3\); its occurrence in S forces
that edge. Thus the nonzero coordinate support of m is a clique.

There are15 maximal cliques, of sizes at most six. These cover all225
nonempty cliques. The independent verifier checks coverage by an
exhaustive traversal of all \(2^{18}\) coordinate subsets, independently
of the producer's maximal-clique routine.

For each support, exact cubic-span elimination produces homogeneous
membership equations I and coordinate cubics \(a_1(m),\ldots,a_{89}(m)\).
Pairs where either coordinate is the zero polynomial are impossible.
Across the15 supports,86 remaining pair systems require exclusion.
For each, the retained certificate is an explicit rational identity

\[
1=\sum_f h_f f+h_j(a_j-1)+h_k(a_k-\rho),\qquad f\in I.
\]

Every multiplier has degree at most three. These identities exclude
solutions over any field extension of characteristic zero; no numerical
root search or bounded coefficient assumption is used.

## Independent replay

- [Primary support and ideal construction](power_single_compress.py)
  uses SymPy/FLINT for exact cube-span coordinates and Singular for the
  unit-ideal tests. [Its checkpoint](power-single-compress/checkpoint.json)
  closes all15 support cases.
- [Certificate export](certify_power_single_compress.py) retains all86
  multiplier identities, their source cases and hashes in
  [the certificate manifest](power-single-compress/certification.json).
- [The independent checker](verify_power_single_compress.py) uses only
  standard-library Fraction arithmetic and the retained sparse
  polynomial engine. It reconstructs the89 independent cubes, every
  membership row space and coordinate map, all134 proportional pairs,
  and the complete clique coverage. It then multiplies out all86
  certificate identities coefficient by coefficient.
- [Independent validation](power-single-compress/validation.json) records
  these checks. A positive control verifies
  \((x+2y)^3=-3x^3+6y^3+3(x+y)^3+(x-y)^3\), so shortening an independent
  cube list is not ruled out by independence alone.

The earlier [pair-line search](power_line_compress.py) and its
[separate Fraction/gcd replay](power-line-compress/validation.json)
cover all3916 lines through pairs of saved forms. They are now a
retained discovery route; the unrestricted new-form certificate above
is stronger. A first generic-parameter sampling gap in that route was
corrected and preserved under its v1 directory. An initial Singular
rendering failure is retained under syntax-v0 and was never treated as
an algebraic exclusion.

## Other completed discovery tests

A new256-form pool obtained by completing squares of quadratic-output
pencils increases cube-span rank from121 to154. Its60-second exact
basis-exchange search makes10,554 exchanges and stays at89 forms.
[Checkpoint](power-basis-pencils/checkpoint.json).

For the ordinary cubic, a width12, depth8 nonmonotone conjugation beam
checks69,747 line values, encounters1,329 distinct coefficient maps, and stays at36
terms. [Checkpoint](linear-beam12/checkpoint.json). A separate primary
calculation of all rational cancellation values on twelve tangent-
normalized coordinate translation lines also stays at36.
[Checkpoint](translation-lines12/checkpoint.json). Neither calculation
is a general affine-orbit optimum or a new dimension bound.

The next cube-compression test must introduce at least two new forms
at once, or change the companion beyond its current linear orbit.
Classifying all cubes in S can also guide those simultaneous changes.
All global simplicity questions remain open.
