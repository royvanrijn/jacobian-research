# Generated elliptic-family discovery

`ecsearch.family_discovery` turns exact family recognition into a bounded
family-discovery pipeline.  The input is a declared construction space, not a
preselected family.  For each target and each generated family, the engine:

1. forms the exact equation `c4(t)^3-j_target*Delta(t)=0`;
2. rejects a family when the homogenized equation has neither an affine root
   nor a root at infinity modulo one of the declared good primes;
3. factors every modular survivor over `QQ` and retains its rational
   parameters;
4. re-specializes each parameter and checks the exact `j`-identity; and
5. checks the weight-4, weight-6, and weight-12 invariant ratios to distinguish
   a `Q`-isomorphism from a geometric `j`-match.

Polynomial Weierstrass families are read from standardized family metadata;
their declared discriminant is checked against the polynomial model before
the search starts.  Six-root Mestre families are generated from normalized
root censuses and construction adapters.  Their even target-`j` equation is
interpolated in `z=T^2`, checked at held-out values, screened modulo primes,
and interpolated and factored exactly only for survivors.

## ICARM 273/282/302 census

This bounded result is registered as
`EC-ICARM-GENERATED-FAMILY-DISCOVERY`. A maintenance-only provenance check
of the committed inputs and output is available without generating families,
running a modular sieve, or factoring a polynomial:

```sh
python3 elliptic-curves/cas/audit_icarm_construction_recognition_artifacts.py
```

Replay the pinned computation with:

```sh
.venv/bin/python elliptic-curves/scripts/discover_record_families.py \
  elliptic-curves/data/family-discovery/icarm_273_282_302.json \
  --output artifacts/generated-results/elliptic-curves/icarm_273_282_302_family_discovery_v1.json \
  --check
```

The construction space contains one polynomial Weierstrass family and 2,333
distinct six-root Mestre families.  The latter are the 2,329 declared census
tuples through normalized diameter 300 together with the distinct outputs of
Fermigier's two-parameter root generator on `1 <= u,v <= 5`, after exact
normalization and deduplication.  The relevant large tuple is therefore not a
hand-entered positive control: the generator emits

```text
(0,29,658,722,981,1036)
```

at `(u,v)=(3,5)` and, by symmetry, `(5,3)`.

The exact result is:

| target | families tested | exact-factor survivors | `Q`-isomorphic matches |
| --- | ---: | ---: | --- |
| ICARM 273 | 2,334 | 113 | none |
| ICARM 282 | 2,334 | 114 | Fermigier `u=11671/42`; six-root Mestre `T=11671/21` |
| ICARM 302 | 2,334 | 146 | none |

The two curve-282 coordinates describe the same fibre with `T=2u`.  Relative
to the public target model, the invariant scales are respectively `882` and
`147`.  Thus the discriminant fingerprints agree prime by prime after the
model change:

```text
v_p(Delta_target) = v_p(Delta_source) + 12*v_p(scale).
```

For example, at the target's repeated primes `2,3,5,7,11,13,23,31`, the
canonical Fermigier source valuations are
`-2,-11,2,-20,4,3,2,2`; applying scale `882` gives the target valuations
`10,13,2,4,4,3,2,2`.

## Evidence boundary

The positive result is an exact rational-isomorphism recognition.  It does
not identify the submitter's search program or prove that local conductor
engineering was used.

The negative results are complete only for the declared construction space
and affine parameters.  They do not exclude a larger root census, a different
construction template or fibration, a parameter at infinity, an isogenous
family, or an unpublished family.  Repeated prime powers remain useful
fingerprints, but without a generated candidate family they do not determine
a parameter.
