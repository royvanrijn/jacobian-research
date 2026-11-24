# The next whole zero piece of the cubic cube-zero algebra

Let \(A=A_{3,3}\) be the free symmetric ternary algebra on one generator
over a characteristic-zero field, subject to the polarized identity
\(L_a^3(b)=0\), where \(L_a(b)=\mu(b,a,a)\). Grade counts internal
ternary vertices, so a grade-\(k\) monomial has \(2k+1\) leaves.

**Theorem CubeZeroGrade29.** The entire grade-\(29\) component of \(A\)
is zero. Together with the retained certificates, this gives
\[
 A_k=0\qquad(19\le k\le29).
\]
Equivalently, in every characteristic-zero symmetric ternary algebra
satisfying the Engel identity, each ternary monomial with exactly
29 internal vertices, evaluated at one repeated element, vanishes.

This statement concerns the whole graded component. The already sharp
inverse-degree bound \(29\) concerns a particular weighted sum of trees
and is unchanged. No assertion about all grades above \(29\), or the
dimension of the whole algebra, follows from this computation alone.

## 1. Exact finite certificate

The [existing quotient proof](CUBIC_CUBE_ZERO_SHUFFLE_EXTENSION.md),
Sections 3–4, constructs exact quotient bases through grade \(18\) and
proves the whole components \(19,\ldots,28\) zero. Its rational basis
dimensions are
\[
1,1,1,1,2,3,4,6,9,13,19,26,34,40,35,16,13,9,2.
\]
The same selected pivot blocks remain invertible over
\(\mathbb F_{1000000009}\), so these bases specialize to the actual
modular quotient. This base-change assertion is rechecked here.

Graft three earlier quotient basis elements with grades summing to
\(28\), using unordered triples. There are exactly **46,928** resulting
grade-\(29\) prequotient columns. Lower-root relations have already
been imposed. Each remaining root relation is specified by six
off-path basis elements and one bottom basis element whose total grade,
plus three, is \(29\).

The [compressed certificate](cubic_cube_zero_grade29_certificate.json.gz)
contains **46,928** such substitutions. The
[verifier](verify_cubic_cube_zero_grade29.py) regenerates their actual
root shuffle sums and verifies, by exact finite-field elimination, that
every selected relation increases rank. Thus they span the complete
46,928-dimensional prequotient. Unused root relations need not be
enumerated for this zero-piece conclusion.

The certificate contains integer substitution identifiers, not
precomputed relation coefficients or a claimed row-reduction matrix.
Its compressed size is 127,596 bytes. The exploratory producer and its
search order are unnecessary for verification.

## 2. Why the finite-field calculation proves the statement

Use the integral polarized presentation over
\(R=\mathbb Z[1/6!]\) from the existing proof. Every grade has a finite
presentation by unordered trees and genuine polarized Engel relations.
The normalized shuffle sums differ from their labelled versions only
by factorial factors that are units in \(R\).

The verifier first rebuilds every rational root-relation space through
grade \(18\), checks its selected pivot block modulo the same prime,
and verifies that all rational coefficients specialize. It then
replays each zero-piece certificate from \(19\) through \(28\) before
using its zero component in any higher calculation. Grade \(29\) uses
only these established lower pieces.

Consequently the new full-rank certificate proves
\[
 (A_R\otimes_R\mathbb F_{1000000009})_{29}=0.
\]
The rank of a finite relation matrix over \(\mathbb Q\) is at least its
rank modulo this prime: a minor nonzero modulo the prime is nonzero over
\(\mathbb Q\). Hence \((A_{\mathbb Q})_{29}=0\). Extension of scalars
gives the assertion over every characteristic-zero field.

This uses whole-prequotient full rank. A modular zero of the inverse
coefficient alone, or an incomplete relation rank, would not establish
this conclusion.

## 3. Reproduction and scope

Run from the repository root:

    python3 research/jacobian/subtree-shuffling/verify_cubic_cube_zero_grade29.py \
      --repo-root . \
      --input research/jacobian/subtree-shuffling/cubic_cube_zero_grade29_certificate.json.gz \
      --output /tmp/cubic_cube_zero_grade29_verification.json

The verifier uses only the Python standard library. It authenticates
the retained prefix sources and reconstructs both the rational prefix
and every required modular zero piece. The
[receipt](cubic_cube_zero_grade29_verification.json) records the actual
replay; the [source ledger](cubic_cube_zero_grade29_source_ledger.json)
identifies the proof, certificate and dependencies.

There is a complete internal exact replay, with the earlier integral
model supplying the characteristic-zero argument. No external human
review, formal verification or historical-priority claim is made.

The known dimensions through grade \(18\) sum to \(235\). To deduce that
the entire free algebra has dimension \(235\) by root grafting, it would
suffice to prove the further whole-piece zeros \(A_{30},\ldots,A_{55}\).
Then induction rules out every higher root: its three nonzero children
would each have grade at most \(18\), giving total grade at most
\(1+3\cdot18=55\). Those further zeros are not part of this theorem.
