# Five compact A1/MW16 family inputs

The five existing anonymous A1/MW16 fibrations now have exact equations with
141–181 coefficient bits and all eighty marked generic sections in the same
coordinates. Together with the [six compact R17 inputs](COMPACT_SIX_R17_ATLAS_2026-09-05.md),
this supplies eleven compact family inputs for separately frozen prospective
experiments. It does not construct new fibrations, a different K3 surface, or
increase any generic or specialized rank.

| Anonymous fibration | Retained presentation | Literal coefficient bits | Constant scaling only | Compact coefficient bits |
|---|---|---:|---:|---:|
| `a1-fibration-01` | `a1-presentation-01` | 6159 | 5174 | 141 |
| `a1-fibration-02` | `a1-presentation-03` | 6185 | 5231 | 143 |
| `a1-fibration-03` | `a1-presentation-05` | 5931 | 4957 | 181 |
| `a1-fibration-04` | `a1-presentation-06` | 5967 | 4965 | 147 |
| `a1-fibration-05` | `a1-presentation-07` | 6182 | 5219 | 141 |

The size is the largest numerator or denominator bit length among the
polynomial A and B coefficients. It is not a specialized minimal-model
height, a conductor, or a measured point-search speedup.

## Exact equations and sections

The [portable atlas](../../artifacts/generated-results/elliptic-curves/compact_five_mw16_atlas_v1.json)
contains all coefficients and rational-function section coordinates. Each
family retains its original marked 16-by-16 rational height Gram in exactly
the exported section order; all five are positive definite of determinant
474. Generic independence is inherited from the source's marked MW16 bases.
The source and partition are the sanitized
[`a1_mw16_family_template_v1.json`](../data/a1_mw16_family_template_v1.json),
whose nine coordinate presentations represent five fibration classes.
The first presentation in template order is used once per class.

For the stored matrix `(a,b;c,d)` and nonzero scale `u`, the convention is

```text
lambda_old = (a*t+b)/(c*t+d)
x_old = u^2*x_new/(c*t+d)^4
y_old = u^3*y_new/(c*t+d)^6
```

The checker proves the two polynomial identities

```text
(c*t+d)^8 A_old((a*t+b)/(c*t+d)) = u^4 A_new(t)
(c*t+d)^12 B_old((a*t+b)/(c*t+d)) = u^6 B_new(t).
```

All sixteen sections are reconstructed over `Q(t)` from the generic native
11952 sections and the retained source-marking words. The checker verifies
each old-section Mobius base map, its point on the chord quartic, the pointed
quartic-to-Weierstrass normalization, the old short equation, and the final
compact short equation. It compares the reconstructed source pencil with
the template exactly. This is a full rational-function computation, not
interpolation or a test at finitely many parameter values.

The generation and full marking replay both passed for all eighty sections.
Replay with Sage 10.9:

```sh
sage -python elliptic-curves/cas/export_compact_mw16_atlas.sage --check \
  artifacts/generated-results/elliptic-curves/compact_five_mw16_atlas_v1.json \
  --reconstruct
```

Omitting `--reconstruct` checks the coefficient identities, point membership,
Gram and source bindings; it does not redo the full marking derivation.
Replay uses the same reconstruction implementation. Independent external
verification, global minimality and exact generic-rank upper bounds are not
claimed by this export.

## Bounded coordinate computation

The first protocol, `compact-mw16-base-v1`, allowed a composite perfect-power
root through 160 bits. On the first family, the gcd of the primitive A/B
discriminants and their resultant, after stripping primes below 10000, was
an exact 56th power with a 673-bit composite root. The worker rejected that
root; the other four were not dispatched.

Exact gcds with rational coefficients already in the sanitized template
split that root into pieces of 17, 95 and 188 bits, including multiplicities.
The 95-bit piece is a square of a 48-bit prime. The separately frozen
`compact-mw16-base-v2` permits this coefficient-gcd decomposition, recursive
perfect-power extraction, and proved factorization of composite pieces
through 192 bits. No 673-bit integer is sent to general factorization.

Each worker retains the 120-second/1-GiB cap, with one worker at a time and
five total. Prime-local PARI auxiliary minimization with the proved prime
list, followed by auxiliary reduction, proposes a rational base map. Weighted
constant scaling and exact forward and inverse elliptic identities decide
the result. The first family's reduction passes the frozen 25-percent gate
before the other four run. All five coordinate workers finished in less
than four seconds each. Generic-section workers have separate 300-second/
1.5-GiB limits and the same first-family success gate.

The [evidence manifest](../../artifacts/generated-results/elliptic-curves/compact_five_mw16_evidence_v1.json)
and [portable bundle](../../artifacts/generated-results/elliptic-curves/compact_five_mw16_evidence_v1.zip)
retain the failed gate, continuation protocols, equations, generic-section
generation, source snapshots and worker logs. Full marking replay logs are
under `artifacts/local/elliptic-curves/compact-mw16-replay-v1/`.

## Using the broader base

[`compact_mw16_specialization.py`](../cas/compact_mw16_specialization.py)
offers `specialize(family, parameter)` without Sage. For reduced `t=n/d`, it
returns the short coefficients with weights `d^8,d^12` and the sixteen
points with weights `d^4,d^6`, checking membership exactly. It rejects singular
fibres, section poles, incomplete section rosters and incorrect points.
This interface currently accepts finite rational compact parameters.

Small parameter boxes in these new coordinates are new populations. The old
104-fibre MW16 null experiment is unchanged. Preserve all five family
identities, freeze any new parameter/prime/point budgets separately, use the
exact 16-dimensional basis when selecting generic half-lattice centres, and
certify specialized independence and any extra directions on each candidate.
The rank-17-only oracle is unsuitable here; the shared general-dimension
geometry is required. This task runs no prospective parameter or point search.

An additional fixed `t=1` audit checks practical input transport and finite
independence. Its first all-or-nothing exporter stopped when the second
family's mod-2 test through prime 1000 failed to certify all sixteen points.
The retained continuation keeps exactly the same five parameters and prime
limit, records inconclusive full-subgroup tests, and certifies only the
independent subsets it can prove. A finite quotient rank deficit does not
prove point dependence or a rank upper bound; specialized 2-saturation can
also affect this detector. The symbolic generic-section certificate does
not depend on this auxiliary specialization test.

The [completed fixed-input certificate](../../artifacts/generated-results/elliptic-curves/compact_five_mw16_input_audit_v1.json)
has finite quotient ranks `16,14,16,16,15` in family order, and exactly
certifies independent subsets of those sizes. All eighty points pass
membership. Full sixteen-point independence on these particular `t=1`
fibres remains unproved by this test for families 02 and 05. No search
parameter was changed to replace these outcomes.
The [validation manifest](../../artifacts/generated-results/elliptic-curves/compact_five_mw16_validation_v1.json)
and [bundle](../../artifacts/generated-results/elliptic-curves/compact_five_mw16_validation_v1.zip)
retain the full marking replay, the first fixed-input failure, and the
completed fixed-input audit with their source snapshots and logs.

A subsequent [odd-prime diagnostic](../../artifacts/generated-results/elliptic-curves/compact_five_mw16_odd_independence_v1.json)
tests only those same two deficient `t=1` inputs. Neither has equal or
sign-opposite point pairs. Both mod-3 and mod-5 quotient tests, using
reduction primes through 1000, retain dimensions 14 and 15 respectively;
neither supplies a full sixteen-point independence certificate.
The [exact short-relation audit](../../artifacts/generated-results/elliptic-curves/compact_five_mw16_short_relations_v1.json)
then tests all ternary lifts of the mod-3 kernels up to sign: four nonzero
words for family 02 and one for family 05. None is an exact rational zero
sum. The ranks of the displayed subgroups therefore remain in `[14,16]`
and `[15,16]`; these are not intervals for the full curve ranks. No
additional parameters, point searches or rational-halving campaign ran.
Both diagnostic certificates replay; their
[evidence manifest](../../artifacts/generated-results/elliptic-curves/compact_five_mw16_independence_diagnostics_v1.json)
retains the protocols and failed initial tuple/list comparison in a separate
startup, followed by its corrected replay.

```sh
python3 elliptic-curves/cas/audit_compact_mw16_odd_independence.py --check \
  artifacts/generated-results/elliptic-curves/compact_five_mw16_odd_independence_v1.json
sage -python elliptic-curves/cas/audit_compact_mw16_short_relations.sage --check \
  artifacts/generated-results/elliptic-curves/compact_five_mw16_short_relations_v1.json
```

```sh
python3 -m unittest discover -s elliptic-curves/tests \
  -p test_compact_mw16_specialization.py
python3 elliptic-curves/cas/audit_compact_mw16_inputs.py --check \
  artifacts/generated-results/elliptic-curves/compact_five_mw16_input_audit_v1.json
```

The curve-finding task owns the next separately frozen prospective experiment
on these compact models. Higher generic rank and a wider family base remain
different research objectives; this result achieves the latter at the input
and coordinate level.
