# Local-conductor parameter recovery

`ecsearch.conductor_engineering` implements the reusable part of the
prime-power fingerprint pipeline for one-parameter elliptic families.  Its
scope is deliberately exact but bounded:

1. derive the family `c4`, `c6`, and discriminant polynomials from polynomial
   Weierstrass coefficients and check the declared discriminant identity;
2. globally minimize the target curve with PARI/GP and compute its selected
   local Tate data;
3. solve the family discriminant congruence modulo each selected prime power
   on both charts of `P^1(Q_p)`;
4. compress complete multiple-root sets into disjoint maximal p-adic balls;
5. combine homogeneous rows `C*a+D*b=0 (mod p^e)` by CRT, exactly
   Gauss-reduce their rank-two kernel, and enumerate a declared coefficient
   box; and
6. verify a recovered fibre by the exact integer identity equating its
   j-invariant to the target j-invariant.

The two-chart partition is important.  The affine chart contains all classes
with `b` a p-adic unit and uses `t=a/b`.  The infinity chart contains the
remaining primitive classes, where `a` is a unit and `b/a` is divisible by
`p`.  Thus denominator-divisible parameters are not silently discarded.

## ICARM curve 282 replay

Run the pinned Fermigier example with:

```sh
python3 elliptic-curves/scripts/recover_conductor_parameter.py \
  elliptic-curves/data/conductor-engineering/icarm_curve282_fermigier.json \
  --output \
  artifacts/generated-results/elliptic-curves/icarm_curve282_conductor_parameter_recovery_v1.json \
  --check
```

The replay proves that the submitted curve is already globally minimal.  At
the selected primes its minimal discriminant valuations are

```text
v_5=2, v_11=4, v_13=3, v_23=2, v_31=2.
```

PARI gives multiplicative types `I2`, `I4`, `I3`, `I2`, `I2`; the reductions
at 11 and 13 are split, while those at 5, 23, and 31 are nonsplit.  The exact
declared branches combine to a single bounded candidate, `u=11671/42`, and
the family and target j-invariants agree exactly.

The input specification records exact branch residues because the valuations
alone do **not** determine them.  For example, `v_11(H)>=4` holds on three
whole residue classes modulo 11, each covering 1,331 leaves modulo `11^4`.
Omitting a `branch` field asks the tool to enumerate every maximal ball for
that prime; the search remains exact within the declared lattice box but can
grow rapidly.

For a new target, the `constraints` array may be omitted and replaced by an
`automatic_constraints` object with `maximum_prime`, `minimum_valuation`, and
an optional `excluded_primes` list.  The tool then selects every minimal-
discriminant valuation meeting that threshold.  Family-specific primes where
denominators, coordinate normalization, or a known fixed factor interfere
should be excluded explicitly and explained in the specification.  When the
target record supplies its bad-prime support, the output also divides the
entire minimal discriminant over that support and records the remaining
cofactor; cofactor one is an exact factorization replay.

## Evidence boundary

A returned exact match proves that a nonsingular family fibre has the same
j-invariant as the target, hence is isomorphic over the algebraic closure.
It does not by itself produce a rational change of variables, prove equality
of quadratic-twist class, or identify the submitter's search algorithm.
Likewise, failure in a bounded coefficient box is not a nonexistence proof.

The local profiles are still useful when no match is found: they distinguish
genuinely expensive high-power residue conditions from cheap repeated-root
balls and fixed divisors.  That distinction should be retained when using
the fingerprints to design a forward search.
