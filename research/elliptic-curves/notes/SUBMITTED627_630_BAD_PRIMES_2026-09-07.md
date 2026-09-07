# Bad primes and conductors for ICARM627–630

The first priority is the four curves submitted by Roy van Rijn, identified
by exact equation matches to the existing research inventory:

| ICARM | Local inventory ID | Family / parameter | Current conductor status |
|---|---|---|---|
| 627 | `new-20260906-90` | `a1-fibration-01`, `-1867/270` | exact, 129 digits |
| 628 | `new-20260906-186` | `11952`, `4286/1881` | exact, 123 digits |
| 629 | `new-20260906-71` | `103b2`, `3726/881` | exact, 124 digits |
| 630 | `new-20260906-40` | `074d9`, `2818/1535` | exact, 120 digits |

The [complete prime lists and exact conductor integers](../../artifacts/generated-results/elliptic-curves/submitted628_conductor_v2/PRIMES_TO_PASTE.md)
are ready to paste into ICARM. Each list belongs to its numbered public entry.
**All four prime lists are complete, including628.** The underlying points and
rank bounds are unchanged. No website update was made by these scripts.

## Exact results

For627–630, complete factorizations reconstruct the entire minimal
discriminant. Every factor has a saved exact primality certificate: ECPP for
large primes and the deterministic small-prime certificate where applicable.
The verifier checks both certificate validity and the root prime it certifies.

At every discriminant prime, Sage's generic Tate algorithm and PARI's
`elllocalred` independently agree on the conductor exponent and Kodaira type.
The local scaling and minimal-discriminant valuation prove the submitted
integral equation minimal at every such prime. All other primes are good by
the complete discriminant factorization. Thus the reported bad-prime supports
are complete and the products `N = product p^f_p` are exact.

- [627 certificate](../../artifacts/generated-results/elliptic-curves/submitted627_630_conductors_v1/curve627.json)
- [628 certificate](../../artifacts/generated-results/elliptic-curves/submitted628_conductor_v2/curve628.json)
- [629 certificate](../../artifacts/generated-results/elliptic-curves/submitted627_630_conductors_v1/curve629.json)
- [630 certificate](../../artifacts/generated-results/elliptic-curves/submitted627_630_conductors_v1/curve630.json)

These close all four of the earlier conductor upper bounds without changing
the equations or making a conductor-record claim. The
[current conductor index](../data/submitted_conductors_current.json)
now selects four exact values. The original three-exact/one-unknown index is retained as historical evidence.

## Completion of628

The former88-digit composite is exactly

```text
1297327695885285412767490764778583 *
7644533052528481157308027511138991115839138041559182477
```

Both factors have saved, verified ECPP certificates. Together with the earlier
factor526689608973452707 and the small factors, they reconstruct the full
minimal discriminant. Only5 has conductor exponent2; the other12 bad primes
have exponent1. The exact conductor is

```text
626331835852237837959942160026415058213286411653252011129135977016421858574282246171665388209084791081551835278559332642950
```

CADO-NFS completed the88-digit split in500.62 elapsed seconds with the configured eight-thread
setting, within the declared3600-second allocation. Its official source is
pinned at commit `f85098631fd047777fbdca94c765d77cb7f12299`; the protocol,
binary hashes, build logs, command, relations and resumable database are retained
under `artifacts/local/elliptic-curves/submitted628-nfs-v1/` and the adjacent
`cado-nfs-source/`. The final certificate replays without CADO or factoring.
The earlier failed attempts below remain historical evidence.

## Bounded discovery and evidence

All four public inputs, the fixed roster and exact existing discriminants are
retained in `artifacts/local/elliptic-curves/submitted627-630-conductors-v1/`.
Cross-GCDs against the201-curve inventory found no useful shared factors.
FactorDB lookups on the four public residuals returned only unresolved composites
and supplied no factorization evidence. External factorization results would
in any case have been treated solely as candidate factors.

The first PARI partial-factorization pass had60 seconds per curve. It completed
627/629 and timed out on628/630. A deterministic GMP-ECM schedule then split
627,628 and630 further. The extended schedule for628/630 remained within the
original600-second ECM allocation per curve;630 completed, while628 exhausted
that allocation. A separate180-second PARI attempt on628's remaining cofactor
also timed out. Its first wrapper invocation emitted no factorization because
the GP stack-reset command shared the input line; the corrected multiline input
and both logs are preserved. A 300-second FLINT quadratic-sieve attempt also
timed out and is recorded in the same evidence directory. No incomplete or
failed attempt proves a prime list or conductor.

Exact factorization of the two relevant generic discriminant polynomials over
Q gives an irreducible degree24 factor in each case; specializing these factors
supplies no further proper GCD with the residual composites. The initial attempt
to coerce a rational factor evaluation directly to an integer failed; the corrected
numerator calculation is retained as a diagnostic, not a prime certificate.

Every ECM attempt retains its input, bound, deterministic sigma, output and
time, with a checkpoint after each attempt and each factor split. The exact
conductor certificates replay without rerunning factor discovery.

```sh
sage -python elliptic-curves/cas/certify_submitted627_630_conductors.sage --check --id 627
sage -python elliptic-curves/cas/certify_submitted627_630_conductors.sage --check --id 629
sage -python elliptic-curves/cas/certify_submitted627_630_conductors.sage --check --id 630
sage -python elliptic-curves/cas/certify_submitted628_nfs.sage --check
python3 elliptic-curves/cas/index_submitted_conductors_v2.py --check
```

After completing628, the renewed [201-curve conductor audit](../../artifacts/generated-results/elliptic-curves/inventory201_conductor_bounds_v1.json) covers the entire inventory against a fresh630-entry catalogue. The [completed priority screen](NEW_CURVE_CONDUCTOR_RECORD_SCREEN_2026-09-07.md) proves24 further exact conductors, including all13 priorities; none beats the recorded benchmark. The earlier187-curve audit remains preserved.
