# Independent norm-square classes can fail the unramified incidence gate

A fixed equation-defined dictionary of twelve classes on each of the eleven
completed panel fibres generates **no additional strict class beyond the
marked generic subgroup**. This is an exact statement about the supplied
dictionary, not a bound on the full strict Selmer group or the curve rank.

Two combinations survive all local squareness conditions at S and the
outside-S probes through 1999, and finite characters prove that each is
globally independent of G. They nevertheless fail Selmer incidence: an
exact, factor-free norm argument proves ramification at some good prime.
No new Sha or rational class is asserted.

This directly tests a new independent class source after the existing
[native carrier audit](RESEARCH_PIVOT_AFTER_THE_CARRIER_PANEL.md) and the
bounded full-class-group failures. The latter remain relevant:
[Hecke's class-group documentation](https://docs.hecke.thofma.com/v0.28/orders/ideals/#Class-Group)
and [Magma's descent documentation](https://magma.maths.usyd.edu.au/magma/handbook/text/1621)
do not provide a verified bypass of the missing class relations on these
fields. No unchanged class-group call was restarted.

## The independent construction

For the retained PARI maximal-order basis 1,b1,b2, take all coefficient
triples in {−1,0,1}³ whose first nonzero entry is +1, excluding (1,0,0).
There are twelve. For each nonzero alpha define

\[
 \pi(\alpha)=N_{K/\mathbf Q}(\alpha)\alpha,
 \qquad N\pi(\alpha)=N(\alpha)^4.
\]

Thus each projected class belongs to the global norm-square cohomology
space. The rule uses the equation and its maximal order; it uses neither
exceptional points nor exceptional-derived classes. The chosen small
coordinate dictionary is basis-dependent, so its counts are not intrinsic
field statistics or proposed rank features.

The [protocol](FRESH_NORM_PROJECTION_PROTOCOL.json) tests the span of these
classes **together with G**, so a candidate cannot be overlooked merely
because it needs generic correction. It imposes zero local squareclass at
2, infinity and every bad prime, then even valuations at every tested
outside-S prime through 1999 dividing a candidate norm. The quotient here
is initially in coefficient space; nonzero coefficient combinations need
not represent independent global classes.

The [capture](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_norm_projection_v1.json)
records the basis, all 132 norm identities, local matrices, valuation probes,
coefficient kernels and finite-character tests. It deliberately leaves
the full support of the two surviving combinations unresolved.

| Case | Generic strict dimension k | Coefficient kernel at S | Kernel after finite outside-S probes | Candidate coefficient dimension modulo generic strict kernel |
|---|---:|---:|---:|---:|
| 074d9 low, case-01 | 9 | 13 | 10 | 1 |
| 103b2 +10, case-02 | 0 | 1 | 1 | 1 |
| 103b2 low, case-03 | 0 | 1 | 0 | 0 |
| 11952 +10, case-04 | 0 | 0 | 0 | 0 |
| 11952 compact low, case-05 | 10 | 17 | 10 | 0 |
| 11952 larger +10, case-06 | 0 | 1 | 0 | 0 |
| 11952 larger low, case-07 | 5 | 9 | 5 | 0 |
| 11952 +10, case-09 | 0 | 1 | 0 | 0 |
| ICARM356 +12 | 1 | 2 | 1 | 0 |
| ICARM385 +12 | 0 | 0 | 0 | 0 |
| ICARM398 +14 | 0 | 0 | 0 | 0 |

The two surviving combinations each increase the global finite-character
rank from 17 to 18. This certifies independence from G, but not membership
in the strict Selmer group. It also illustrates why counting a bad-place
kernel as a Selmer lower bound would be wrong.

## A factor-free ramification certificate

Let alpha_i=h_i(theta)/d_i, with h_i an integer polynomial of degree less
than three. Select a product of projected classes pi(alpha_i), allowing
multiplication by any generic point class. Write N_i=N(alpha_i).

Suppose a rational prime p satisfies:

1. p is odd and does not divide the defining cubic's discriminant, any
   coefficient denominator d_i, or the content of any selected h_i;
2. v_p(N_i) is odd for one selected i;
3. p divides none of the other selected norms N_j.

Then this projected product is ramified at a prime of K above p.

**Proof.** At such p, the maximal-order residue algebra is the product of
the distinct factors of the defining cubic modulo p. The nonzero polynomial
h_i has degree less than three, so cannot vanish in every component. In
at least one component alpha_i is a unit. All other selected alpha_j are
units in every component because their norms are p-adic units. The rational
scalar product of the selected N_j has odd p-adic valuation. Since K is
unramified at p, its valuation in the chosen component is also odd.
Consequently the projected product has odd valuation there. At an odd
prime this gives a ramified quadratic extension. Generic point classes
have even valuations at this good prime and cannot cancel the obstruction.

An explicit p is unnecessary. For each selected N_i, remove by repeated
gcd every prime factor shared with the cubic discriminant, denominators,
contents and other selected norms. Call the positive remaining integer R_i.
If R_i is not a square, unique factorization guarantees a prime satisfying
conditions 1–3. A floor-square-root inequality proves nonsquareness even
when R_i is a large unfactored composite.

The [fixed survivor audit](FRESH_NORM_ISOLATION_PROTOCOL.json) applies this
argument only to retained combinations, generating no new candidates.
It finds an isolated nonsquare remainder of **409 bits** on case-01 and
**246 bits** on case-02. The certificate retains each removal gcd, the
remaining integer and its floor square root. No integer factorization runs.

Each candidate coefficient quotient has dimension one. Therefore every
combination involving its nonzero quotient vector is the same ramified
class up to generic strict classes and squares. This proves exclusion of
the whole generated strict excess. In higher dimension, separately
ramified basis vectors could cancel each other's ramification; basis-wise
rejection would not suffice. The verifier explicitly checks this dimension
gate rather than applying that invalid inference.

## What changes in the mechanism search

* **Incidence:** this particular independent dictionary is exhausted on all
  eleven fibres. Two apparently promising global classes are not Selmer
  classes at all. Their failure is not a rational/Sha distinction.
* **Solubility:** no new admissible additional class reaches a CT calculation.
  The necessary large strict blocks established by the panel remain real
  mathematical targets, but these small basis combinations do not find them.
* **Next construction gate:** prioritize an unramified character or a proved
  ideal-parity circuit. A larger small-coefficient dictionary needs a reason
  to improve that gate; square norm and local correction at S alone are
  insufficient. The isolated-norm certificate can reject many candidates
  before expensive support factorization.

No prospective rank or visibility score follows. This is a bounded
incidence exclusion and a reusable exact rejection lemma. It does not
exclude different bases, larger dictionaries, other class constructions,
or a simultaneous-solubility carrier for the genuine additional block.

## Verification

The [verification certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_fresh_norm_projection_verification_v1.json)
passes all 132 rational norm-projection identities, 4216 finite local
signature replays, real signs and outside-S valuations, all coefficient
kernel calculations, and 4224 independent finite-field character
evaluations. Norm and product identities use independent rational cubic
arithmetic; local coordinates reuse the established LocalSquareclasses
backend. Both factor-free exclusions are exact integer proofs.

```sh
timeout 60 sage -python elliptic-curves/rank-jump/verify_fresh_norm_projection.py check
```

All workers finish within the frozen caps. No exceptional point input,
new original parameter, point search, class-group call or active-search
modification occurs. The panel's overall coverage remains eleven of sixteen.
