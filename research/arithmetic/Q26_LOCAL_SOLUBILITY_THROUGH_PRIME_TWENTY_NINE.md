# The degree-twenty-six good locus is locally soluble through prime29

The normalized degree-twenty-six affine-trace accessory curve has
infinitely many good points over Qp for every prime p<=29, and over
the reals. At5, where the naive integral reduction is highly singular,
one explicit translated polynomial gives a simple Hensel certificate.
The other new primes have small nonsingular affine residues.

These are necessary-condition tests for rational existence. They do
not give a Q-rational point, solubility at every prime, or a Hasse
principle. They do not change the finiteness theorem for normalized
rational models in this slice.

## The direct5-adic certificate

Let F(s,t) be the frozen primitive integer polynomial of the
[high-genus accessory theorem](Q26_HIGH_GENUS_AFFINE_TRACE_CURVE.md).
Its full [coefficient input](q26_high_genus_affine_trace_curve_inputs.json)
has already been exactly bound to the uniform affine-trace construction.

**Theorem 1.** The degree-twelve polynomial

\[
             P(X)=5^{-40}F(1056+1200X,238)                 \tag{1}
\]

has integer coefficients and satisfies

\[
                    \overline P=3X+4X^2\quad\text{in }\mathbb F_5[X].
 \tag{2}
\]

More generally

\[
 H(X,V)=5^{-40}F(1056+1200X,238+5^{41}V)\in\mathbb Z[X,V]
 \tag{3}
\]

has the same reduction(2), independent of V. Thus for every V in Z5
there are unique roots X0(V) congruent0 and X3(V) congruent3, in their
respective residue classes. Each gives a smooth Q5 accessory-curve
point, and each family contains infinitely many good models.

**Proof.** Exact integer expansion of the numerator in(1) gives
the following coefficient valuations, in ascending order:

\[
                41,40,40,41,42,44,46,48,51,52,54,56,58.
 \tag{4}
\]

Dividing every coefficient by5^40 is therefore exact. The retained
[input](q26_local_solubility_through_29_inputs.json) gives every
resulting coefficient, and the checker verifies the entire identity
and reduction(2), without using the discovery blowup tree.

Since F has integer coefficients, every positive V-degree j term
in the numerator of(3) has a factor5^(41j). After division by5^40
it is still integral and divisible by5. The constant V-term is P.
This proves(3) and its reduction for every V, without finite-parameter
sampling. The roots0 and3 of(2) have derivatives3 and2. The simple-root
Hensel lemma, [Milne ANT3.00, Proposition7.31,p113](sources/milne-ant300.pdf),
gives the two unique lifts. Their X-derivatives remain units.

The affine coordinate change has nonzero rational scale1200, so
these are smooth characteristic-zero points of F=0. Distinct V give
distinct slopes t=238+5^41V, hence infinitely many points in each
branch. The geometric good open is nonempty by the
[uniform theorem's Q29 point](ODD_AFFINE_CRITICAL_TRACE_CURVES.md),
and the accessory curve is geometrically integral by the high-genus
theorem. Its bad complement is therefore finite. Removing that
finite set leaves infinitely many good Q5 models in both branches. ∎

The fixed-slope roots at t238 are not individually asserted to be
good; the varying-slope families and finite bad locus supply that
conclusion. There is no claim that these are every Q5 branch.

## Six additional simple affine residues

**Proposition 2.** The following exact residues satisfy F=0 and
have the indicated nonzero gradient modulo p:

| p | s | t | F_s | F_t |
| --- | ---: | ---: | ---: | ---: |
| 7 | 0 | 6 | 4 | 0 |
| 11 | 2 | 10 | 2 | 5 |
| 13 | 1 | 0 | 9 | 3 |
| 17 | 0 | 1 | 5 | 14 |
| 19 | 0 | 8 | 0 | 14 |
| 23 | 0 | 9 | 13 | 11 |

**Proof.** The checker evaluates the full frozen integer F and both
partial derivatives at each residue. Holding one variable near its
residue and using the variable with unit derivative, Hensel gives
infinitely many smooth local points. The same finite-bad-complement
argument as in Theorem1 leaves infinitely many good points over
each corresponding Qp. The residue itself need not satisfy every
critical-good-reduction condition. ∎

The existing Q29 point is compatible with this primitive equation:
F(22,2)=0 modulo29 and its gradient is(14,24) modulo29. This is also
rechecked, as a consistency control of the imported local result.

## Complete scope of the tested places

The [earlier exact toric and real-branch theorem](Q26_SMALL_PRIME_LOCAL_BRANCHES.md)
supplies good Q2, Q3 and real points. The uniform theorem supplies
the Q29 point. Together with Theorem1 and Proposition2, these cover
exactly all primes at most29:

\[
                    2,3,5,7,11,13,17,19,23,29.
\]

No local obstruction occurs at any of those places. Every listed
place has infinitely many good points. This does not assert local
solubility at every larger prime or rational global solubility.
Real solubility also does not mean all finite critical points are real.

All conclusions concern the normalized affine-trace accessory slice.
Its S12 monodromy over the slope line differs from the A26 monodromy
of a reconstructed Laurent pencil over its own parameter. Rational
twists and coordinate rescalings remain outside the normalized
coefficient-point count.

## Exact evidence

The [checker](verify_q26_local_solubility_through_29.py),
[inputs](q26_local_solubility_through_29_inputs.json),
[receipt](q26_local_solubility_through_29_verification.json) and
[source ledger](q26_local_solubility_through_29_source_ledger.json)
retain the direct integer specialization, all coefficient valuations,
both simple-root digit controls, every residue/gradient and the exact
list of primes covered. The proof of the all-V family is the integer
binomial argument, not sampling of V. No exhaustive search or
absence of other branches is used.

No rational Q point, effective rational-point list, all-prime Hasse
claim, formal verification, independent external replay, external
human review or novelty claim is made. Earlier frozen proof packages
and their metadata remain unchanged by this companion note.
