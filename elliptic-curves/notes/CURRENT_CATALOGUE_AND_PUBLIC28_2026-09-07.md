# Current catalogue matches and a public28 reproduction

The separately downloaded [ICARM database](https://elliptic-rank.icarm.cloud/database.json)
at **2026-09-06T22:30:32.628536+00:00** contains620 equations. Its SHA256 is
`b378d458aab4b0d09f9fc8e2e382baa7a20236715a8099c52b4573f00239126e`.
The cached web listing still showed614; the archived JSON download is the
input to this check. This update follows discovery and changes no frozen
candidate selection, point budget or earlier catalogue certificate.

## Exact novelty update

All200 V20 inventory equations were compared with the620 downloaded equations.
The [exact rational comparison](../../artifacts/generated-results/elliptic-curves/inventory200_current_catalogue_comparison_v1.json)
and [independent Sage replay](../../artifacts/generated-results/elliptic-curves/inventory200_current_catalogue_sage_replay_v1.json)
agree:198 have no rational-isomorphism match; two are now catalogue-known.
All593 previously pinned equations are still present as identical equations.

| Local ID | Local family and parameter | Earlier certified bound | Catalogue match | Catalogue-reported bound |
|---|---|---:|---:|---:|
| `new-20260905-12` | `08234`, −506/9 | 23 | 600 | 23 |
| `new-20260906-188` | `11952`, 110314/102227 | 27 | 619 | 28 |

The snapshot records entry600 at2026-09-06 15:45:21 and entry619
at2026-09-06 22:22:59. These are reported submission times, not a proof of
first-discovery priority. The old593-equation absence claims remain true
for that snapshot. Present-day claims must acknowledge these two matches.
The inventory remains a collection of distinct equations and valid point
certificates; it must not be described as200 currently catalogue-absent curves.
Absence from the larger catalogue still does not establish universal novelty.

## Independent reproduction of the public28 points

Entry619 supplies28 points and credits submitter `Mundoamundo`. We verify
them separately from all prospective experiments. Its equation is

```text
[1, 0, 0,
 -138217506563605889872043263854781697377016668744941773513580,
 20242130248207841544219808693189001887421294913221456339098328387159912218373571988160400]
```

The exact map to our retained short equation has scale u=1:
`X=x+b2/12`, `Y=y+(a1*x+a3)/2`. Every published point is checked on the
public equation and after transport. The [public reproduction certificate](../../artifacts/generated-results/elliptic-curves/inventory188_public28_reproduction_v1.json)
replays the earlier27-point proof and then audits the55-point union of those27
local points with the28 transported public points, using fixed good primes
through997. It proves rank at least28 both for the public basis and the union.
Its selected union basis consists of all27 old points and public point index26
(zero-based). The finite quotient rank28 is not an upper bound on the curve
or on the rational span of the full union.

The [independent Sage check](../../artifacts/generated-results/elliptic-curves/inventory188_public28_sage_replay_v1.json)
rechecks point membership and transport, enumerates each complete finite
elliptic-curve group, forms its quotient by doubles and computes the resulting
binary matrix ranks. It separately verifies irreducibility of the2-division
cubic at the torsion-exclusion prime. All checks pass.

**This is a public-data reproduction, not a new rank28 discovery or a blind
recovery.** The points stay outside candidate selection and prospective search
inputs. Our earlier bounded search found this equation and certified27; the
published witness now provides a concrete missed-direction control. The [completed exact visibility audit](../../artifacts/generated-results/elliptic-curves/inventory188_public28_visibility_v1.json)
checks both signs of the certified extra public point26 against all49 original
completed charts. All98 observations lie outside the125000 box. The minimum
affine coordinate height is21,632,242,813,382, at chart13 with positive sign,
with exact coordinate13349368513299/21632242813382. There is no
visible-but-unrecorded discrepancy. This excludes only these representatives;
other translates or representatives may be far cheaper. No further point or
parameter search is launched by this audit.

An initial build compared an in-memory tuple representation of the old
certificate directly with JSON lists and stopped before writing a certificate.
The preserved source and failure log remain under the local download directory.
The V2 audit compares their canonical JSON values; the mathematical inputs,
prime bound and point lists are unchanged.

## Reproduction

Sources are under `elliptic-curves/cas/`:

```sh
python3 elliptic-curves/cas/audit_inventory200_current_catalogue.py --check
python3 elliptic-curves/cas/audit_inventory188_public28_v2.py --check
```

The independent Sage sources are `verify_inventory200_current_catalogue.sage`
and `verify_inventory188_public28.sage`. Their already completed supervisor
records, the original download and metadata are under
`artifacts/local/elliptic-curves/inventory200-current-catalogue-v1/`.
Do not overwrite the frozen outputs when reproducing; use an isolated copy.


## Current inventory with explicit publication status

The completed near-finalist trial adds the catalogue-unmatched23-point curve
ID201. Historical V21 therefore contains201 distinct curves with all point
proofs replayed. The [V22 JSON](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v22.json)
and [CSV](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v22.csv)
preserve every discovery ID, mark ID12 and ID188 as current catalogue matches,
and upgrade only ID188 to the separately reproduced public28-point basis.
Its original local search bound27 and original source remain explicit.
There are199 catalogue-unmatched equations. Their largest local certified
bound remains27; the public28 reproduction is not counted as a new discovery.

[The V22 replay](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v22_replay_v1.json)
binds all200 unchanged point bases to the complete V21 replay, separately
rechecks the changed28-point basis, checks current catalogue comparisons and
verifies every CSV cell. The public basis also has its independent Sage
finite-group proof above. The exporter is
`elliptic-curves/cas/export_publication_aware_inventory_v22.py --check`.

The subsequent [fixed translation and own27 geometry control](INVENTORY188_CHART_COVERAGE_2026-09-07.md) checks coverage without another point search.
