# Arithmetic obstruction to cotangent and canonical shear words

14 November 2026. **QFL3.** Assume the arithmetic nonlifting theorem for
Long's map in Kanel-Belov--Fedorov, Theorem 4.2. Then Long's polynomial
symplectic endomorphism cannot be a finite composition of cotangent lifts
of polynomial Keller maps, positional gradient shears, and linear
symplectic changes. There is no bound on the degrees or number of factors.

This is a corollary of the cited October 2026 preprint, with an explicit
arithmetic compatibility proof for the entire generating class. It does
not exclude arbitrary bounded formal lifts from QFL2.

## The arithmetic compatibility condition

Let f be a polynomial symplectic map defined over a characteristic-zero
field. Say that a finite Weyl endomorphism Phi is compatible with f on
arithmetic centers if both are defined over an integral finitely generated
Z-algebra Lambda and, after one localization of Lambda, every sufficiently
large positive-characteristic reduction satisfies
\[
 \Phi_\ell(z_i^\ell)=
              \iota_\ell(f_i^{[\ell]}),\qquad
 \iota_\ell(z_i)=z_i^\ell.                           \tag{1}
\]
Here f_i^[ell] raises its scalar coefficients to the ell-th power.
The polynomial variables on the right are the central variables.

If Phi and Psi satisfy (1), their composite satisfies (1) for the
corresponding composite polynomial map. This follows by substitution on
the center and compatibility of coefficient Frobenius with composition.
Any finite word admits one common arithmetic coefficient ring: adjoin
its finitely many coefficients and inverses needed for its formulas.

Thus it suffices to verify (1) for the stated generators.

## Cotangent lifts of arbitrary Keller maps

Let F=(F_1,...,F_n) in K[q]^n have Jacobian determinant a nonzero constant.
Write J=DF, B=J^{-1}, and put
\[
 \delta_i=\sum_j B_{ji}(q)\partial_{q_j}.
\]
The delta_i commute and delta_i(F_j)=delta_ij. They give the exact
Weyl endomorphism
\[
 \Phi_F(q_i)=F_i(q),\qquad
 \Phi_F(D_i)=\widehat\delta_i
                  =\sum_j B_{ji}(q)D_j.             \tag{2}
\]
Its classical symplectic map is the cotangent lift
(F(q), B(q)^t p), with the components understood as in (2).
All coefficients remain polynomial. The determinant and the finitely
many scalar denominators can be inverted in Lambda.

In characteristic ell at any such good reduction, delta_i^ell is an
ordinary derivation. It kills every F_j. Invertibility of J therefore
forces delta_i^ell=0 on K[q].

The order-ell leading symbol of hat(delta_i)^ell is
sum_j B_ji(q)^ell D_j^ell. Subtract this central expression. The
difference has differential order less than ell. It commutes with every
q_j because delta_i^ell=0. A PBW operator of order less than ell
commuting with all q_j is multiplication by a polynomial: commutation
differentiates each positive momentum exponent, which is less than ell.
That multiplication polynomial is zero by application to 1. Hence
\[
          \widehat\delta_i^\ell
                    =\sum_j B_{ji}(q)^\ell D_j^\ell. \tag{3}
\]
Together with Phi_F(q_i^ell)=F_i(q)^ell, this is exactly (1).

This is the elementary crystalline p-curvature calculation; compare the
central identification in
[Kanel-Belov--Fedorov, Section 2.1](https://arxiv.org/html/2610.33919v1).
It uses no degree cutoff on F and does not assume F is invertible.

## Positional gradient shears

For any polynomial H(q) in characteristic zero, the assignment
\[
                 q_i\mapsto q_i,\qquad
                 D_i\mapsto D_i+\partial_i H         \tag{4}
\]
is an exact Weyl automorphism. Its inverse replaces H by -H.

In prime characteristic the differential-operator identity
\[
       (D_i+g)^\ell=D_i^\ell+g^\ell+
                                      \partial_i^{\ell-1}g
\]
holds when g is a position polynomial. See
[Bavula, Theorem 1.3](https://eprints.whiterose.ac.uk/id/eprint/8651/3/Bavula_automorphisms.pdf).
Other position variables can be put in the coefficient field.
For g=partial_i H the derivative term is partial_i^ell H=0.
Thus (4) induces exactly the Frobenius-twisted classical shear on the
center. This works for arbitrary H and all good primes.

Fourier conjugation gives the analogous statement for momentum shears.
Linear symplectic transformations also lift compatibly at every odd
good prime: the ell-th power of a linear combination of Weyl generators
is the corresponding sum of ell-th powers in odd characteristic.
This follows from the restricted Heisenberg identity, since its
commutators are central. Translations have the same compatibility.

Therefore arbitrary finite words in these operations satisfy (1).

## Consequence for Long's map

[Kanel-Belov--Fedorov, Theorem 4.2](https://arxiv.org/html/2610.33919v1)
excludes any finite characteristic-zero Weyl endomorphism whose reductions
on a Zariski-dense arithmetic set have Long's prescribed Frobenius center
map. The theorem also excludes every matrix size.

If Long's classical map were a word in the generators above, equations
(2)--(4) and composition would give precisely such an endomorphism.
A common arithmetic ring can include Z[1/30], as required in that
theorem. The compatible good reductions contain a dense open set,
contradicting Theorem 4.2.

In particular:

- Long's map is not a cotangent lift in any source and target coordinates
  obtained by these exact arithmetic-compatible canonical automorphisms.
- Adding arbitrarily many canonical shears, raising their degrees, or
  composing more Keller cotangent factors cannot repair that route.
- Invertible canonical changes from the stated group can be moved across
  a proposed factorization because their inverses are in the same group.

This conclusion is conditional on the cited preprint theorem. It makes
no claim that every polynomial symplectic automorphism admits such an
arithmetic-compatible Weyl lift.

## What remains outside the exclusion

[QFL1](PROOF.md) supplies formal polynomial lifts for Long's full classical
map. [QFL2](BOUNDED_REALIZATION.md) shows that one common finite spatial
degree bound would yield a genuine Weyl counterexample with automatic
nonsurjectivity. Neither statement requires (1).

A formal lift can have the same classical h=0 reduction while its
nonzero-parameter Frobenius center map changes. QFL2's exact shear control
exhibits that difference even for a bounded polynomial quantum family.
Thus the present exclusion does not prove that every bounded formal
quantization of Long's map is impossible.

The remaining construction target is a bounded full lift outside the
arithmetic-compatible cotangent-word class, or a structural theorem
excluding every bounded representative. No bound, new Weyl counterexample,
formal verification, independent replay, external review or novelty
priority is claimed here.
