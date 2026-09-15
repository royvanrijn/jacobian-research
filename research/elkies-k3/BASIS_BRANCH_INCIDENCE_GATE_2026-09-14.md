# Testing explicit basis relations as exceptional branch seeds

The [finite exceptional-locus theorem](JOINT_HALVING_EXCEPTIONAL_BRANCH_LOCUS_2026-09-14.md)
requires actual low-degree incidence points. A direct source would be two
independent generic sections specializing to zero at the same fibre:
their halves are then both O. Another source would be one such relation
and an explicit nonzero2-torsion value.

The first fixed generic-basis gate finds neither source on any of the four
retained MW17 parents. On Q80, all twelve individual basis pole divisors
have **exactly one** specialization-kernel direction, with no nonzero
residue-field2-torsion. These particular divisors cannot supply either
of the two simultaneous conditions. This is not an exclusion of all
Mordell--Weil combinations or of the required two-gain construction.

## Exact scope of the construction test

Use only the existing17-section bases and their original parameter:

| Parent | Individual pole loci plus pair-equality loci | Accepted reduction prime | Modular joint pairs replayed over Q | Modular mixed pairs replayed over Q |
|---|---:|---:|---:|---:|
| Direct11952 Q80 | 148 | 157 | 148 | 31 |
| Published R17 | 136 | 193 | 143 | 79 |
| Recovered Curve302 parent | 138 | 139 | 3154 | 1111 |
| X1092 class1 | 137 | 131 | 395 | 216 |

For each section P_i, retain its pole divisor when nonempty. For each
unordered distinct pair, retain the zeros of x(P_i)-x(P_j) after cancelling
its denominator. They represent P_i=O or P_i=+/-P_j respectively. Any
two different labels give independent parity vectors e_i or e_i+e_j.
The actual sign can be chosen over each characteristic-zero residue field.

For the mixed test use the numerator of 2y(P_i)+a1*x(P_i)+a3, after
cancellation; its zeros are nonzero2-torsion values of this basis section.
Using y alone on the long Curve302 Weierstrass equation would be wrong.
There are no exact common factors between any two distinct relation loci,
or between a relation locus and any of these17 torsion loci. The minimal
infinity chart also has17 distinct finite abscissas and no zero ordinates,
so it supplies neither type of simultaneous incidence.

The test excludes these incidences at every finite algebraic degree,
not only degree at most4. Its restriction is the fixed basis-derived
relation list, not a degree or height bound on the entire MW group.

## Twelve Q80 specialization kernels

The nonempty Q80 basis pole loci belong to basis indices

```
1,2,3,4,5,7,8,9,10,11,14,15     (zero-based source order).
```

The first two are irreducible quadratics; the other ten are rational
points. Their reduced abscissa denominators are q_i², ordinate
denominators q_i³, and both numerators are coprime to q_i. Thus P_i
specializes to O, exhibiting the nonzero class e_i in

```
K_b = ker(M/2M -> E_b(k(b))/2E_b(k(b))).
```

For each divisor, the fixed ascending prime pool5..997 supplies simple
degree-one residue places of k(b), with a smooth parent reduction.
Evaluate all17 generic sections, including O when a coordinate has a
pole. Enumerate the full finite elliptic group and its quotient by doubles.
The combined quotient matrices have rank16 in all twelve cases. Their
kernels are precisely span(e_i). Rational-root Kummer characters independently
agree with the quotient row spaces at every retained place. The known
exact pole gives the lower bound; these finite quotients give the upper.

The last primes used are respectively

```
229,181,173,193,227,173,131,191,173,157,113,223.
```

Each divisor also has a smooth degree-one reduction with odd elliptic-group
order, so E_b(k(b))[2] is zero. The witness primes are respectively

```
53,23,23,23,89,29,59,53,71,59,31,29.
```

Hence none of these twelve branch fields supports two independent common
halving parities, and none supports a half together with nonzero2-torsion.
This does not prove zero gain on every cover containing such a branch.
A single common trace and a nonzero branch class at a different branch
remain possible in the relevant global scope.

## Reproducibility and arithmetic safeguards

The [producer/checker](scripts/audit_basis_branch_incidence.py) uses only
retained generic equation and basis fields. It does not use specialized
rank labels, parameters or exceptional point coordinates in the source
packets. It constructs every relation polynomial, not a sample of them.

All relation and torsion polynomials are made monic. The first prime in
131..997 at which all coefficients are integral is selected separately
for each parent. This ensures that a genuine common characteristic-zero
factor cannot disappear on reduction. Every modular common-factor pair
is then checked by exact rational polynomial gcd. No modular common
factor is promoted to an arithmetic incidence.

The initial fixed131 prototype did not certify the leading-coefficient
condition and is **not** accepted as evidence. Its scripts and outputs
are [preserved as unaccepted preflight](../artifacts/generated-results/elkies-k3-basis-branch-incidence-v1/preflight/probe_archive.json).
The final calculation was rerun with the monic integrality gate; notably,
Curve302 then has many more modular coincidences, all rejected exactly.

The [accepted certificate](../artifacts/generated-results/elkies-k3-basis-branch-incidence-v1/result.json)
contains the exact polynomials, residue matrices, full point reductions,
source hashes and limits. The complete calculation took about0.8 seconds
under20 CPU seconds and1GiB, using Python and SymPy1.14.0. Reproduce with:

```
.venv/bin/python research/elkies-k3/scripts/audit_basis_branch_incidence.py --record /tmp/basis-branch-incidence-replay.json
```

The two finite-quotient computations share a harness; no independent
whole-theorem replay or formal verification is claimed. Full parent
rank/saturation and the specialization interpretation remain inherited
or written mathematical dependencies.

The next construction must go beyond these individual basis poles and
pairwise signed equalities, for example by producing an exact additional
half of a nonzero specialized combination. Blindly enlarging a coefficient
box or enumerating all parity pairs is not justified by this finite miss.
