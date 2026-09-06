# Conductor bounds for the complete187-curve inventory

Every curve in inventoryV14 now has an unconditional conductor upper bound.
The [complete audit](../../artifacts/generated-results/elliptic-curves/inventory187_conductor_bounds_v2.json)
passes187 exact integral translations and1656 local Tate calculations.
[Independent PARI replay](../../artifacts/generated-results/elliptic-curves/inventory187_conductor_bounds_pari_replay_v1.json)
agrees with all local conductor exponents, Kodaira symbols, minimal
discriminant valuations and reconstructed bounds. This closes the coverage
gap beyond the earlier31-curve and next12 audits; it supplies no new rank.

## Exact argument

The recorded short models have integral discriminant and coefficient
denominators supported only on2 and3. Exact translations produce integral
models with the same invariants. At2,3 and each discriminant prime through
10000, compute the local conductor exponent by Tate's algorithm. Write the
remaining positive discriminant cofactor asR. At every unprocessed prime
p>=5, the conductor-discriminant inequality gives
f_p<=v_p(Delta_min)<=v_p(Delta_displayed). Consequently

    N divides (product of computed p^f_p) * R.

This proves the upper bound without assuming thatR is prime or squarefree.
An exact conductor is asserted by this audit only whenR=1. The fixed audit
uses one worker and a120-second cap, with a checkpoint after every curve.
Build and exact rerun took4.462 and2.788 seconds respectively.

The first version failed before completing a row in Sage's generic Tate
`uniformizer` branch on a nonintegral short equation. Its frozen source,
protocol and failure log remain under
`artifacts/local/elliptic-curves/inventory187-conductor-bounds-v1`.
Version2 first performs the exact integral translations. Its corrected
results are independently checked with PARI's separate local algorithm.

## New27 curve

For [ID186 at11952 parameter4286/1881](FULL11952_NEW_RANK27_2026-09-06.md),
the [separate minimal-model certificate](../../artifacts/generated-results/elliptic-curves/full11952_conductor_bound_v1.json)
gives

    N <= 626331835852237837959942160026415058213286411653252011129135977016421858574282246171665388209084791081551835278559332642950.

The invariant gcd is75;3 is good. At5 the curve has typeIV and conductor
exponent2. Every other bad prime is multiplicative, so its exponent is1.
The trial discriminant factors are
2^15,5^4,7^6,13^3,37^2,61^2,89^2,149^2,179^3,4919; the106-digit
remaining cofactor is retained exactly in the certificate. Both independent
local algorithms confirm the bound. The conductor itself remains UNKNOWN.

The pinned593-entry catalogue has33 curves with rank lower bound at least27.
Of these,29 have recorded conductors:24 exceed the displayed bound, and5
(IDs363,400,401,402,585) do not. Four lack conductors. This is a comparison
against recorded values, not a universal ranking or conductor record.

## Complete inventory comparison

None of the187 upper bounds is below the recorded catalogue minimum for
its rank threshold. The smallest bound among our rank-at-least27 curves
belongs to ID40 and has120 digits, while the new ID186 bound has123 digits.
For thresholds22,23,24,25,26,27, the smallest inventory bounds belong to
IDs36,162,129,54,92,40 respectively; the JSON retains exact values.
These are minima of the computed **upper bounds**, not proved minima of the
actual conductors. Bounds above a benchmark leave possible improvements
unresolved. The earlier exact conductor proof for ID36 remains separate.

Reproduce with Sage Python:

    elliptic-curves/cas/audit_inventory187_conductor_bounds_v2.sage check
    elliptic-curves/cas/verify_inventory187_conductor_bounds.sage --check
    elliptic-curves/cas/certify_full11952_conductor_bound.sage --check

The audit freezes the full roster and catalogue snapshot. It performs no
large-cofactor factorization, descent, parameter selection or point search.
No isolated portable replay of this supplement is claimed yet.
