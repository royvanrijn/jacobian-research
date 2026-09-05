# A new rank-at-least-22 curve with an exact 76-digit conductor

The compact MW16 family `a1-fibration-05` at parameter `3/17` gives the
global minimal equation

```text
y^2 + x*y = x^3
  - 182451976602578656424609725499140*x
  + 710003150253794219215652666162794189038512805392.
```

It has **22 independently certified rational points** and exact conductor

```text
5651319165610115296564894984823807468501395910978185811787187868498410205790.
```

The [portable proof](../../artifacts/generated-results/elliptic-curves/small_conductor_rank22_proof_v1.json)
contains all 22 points on both the short and integral models, the exact
independence certificate, recursive primality certificates, and local
conductor arguments. Its stable inventory ID is `new-20260905-36`.
The [Python checker](../cas/certify_small_conductor_curve.py) replays the
mathematical proof without Sage, factorization, numerical heights or search.
A [standalone Sage curve file](../../artifacts/generated-results/elliptic-curves/new_rank22_small_conductor_curve.sage)
loads the integral equation and all 22 points directly.
The [evidence manifest](../../artifacts/generated-results/elliptic-curves/prospective_mw16_next12_evidence_v1.json)
and adjacent ZIP retain the inputs and source dependencies.

## What is new, and what is a record comparison

There is no rational-isomorphism match in the pinned
[ICARM catalogue](https://elliptic-rank.icarm.cloud/curves), fetched on
2026-09-05 at 18:57:46 UTC with 586 equations. Its raw SHA256 is
`1ec915b1d108f906791f5361f8150d328ce96e5f41d95d9e78d9a354e175e53a`.
The comparison is against rational isomorphism classes, not textual equation
equality. Catalogue absence is not a proof that nobody has encountered the
curve elsewhere. ICARM is supported by NSF grant DMS2425401.

The new conductor would be **third among rank-at-least-22 entries with
recorded conductors** in that snapshot. The two smaller entries are:

| ICARM ID | Recorded conductor |
|---|---|
| 376 | `176634054787705380330890095321066714644533156566592717278673628271112803086` |
| 575 | `2229394614546176466599146730229608949237869287138056463724778634634017351560` |

IDs 537, 543, 545 and 581 have rank lower bounds at least 22 but no recorded
conductor. No absolute third-smallest-in-the-world claim follows. The earlier
584-row snapshot lacked the conductor of existing entry 575; using it would
incorrectly describe the current placement as second. The refreshed gate
rejected that stale comparison before launching a follow-up worker.

The new conductor is about 31.99 times the recorded rank-at-least-22 minimum.
The recorded rank-at-least-23 minimum, entry 539, is

```text
95221940916620808088830143485247567246687065418433849627482469977720285740879030062.
```

One additional independent point on the new curve would beat that recorded
minimum. **No 23rd independent point or exact rank is proved here.** This is
a concrete low-conductor near-record result; the new rank-at-least-28/32
targets remain open. The catalogue's rank and conductor fields are reported
as external comparison data, not independently reproved for every entry.

## Exact proof

The short model has coefficients

```text
A = -8757694876923775508381266823958721/48
B = 613442721819278218538866218950317441901175299796769/864.
```

The rational group isomorphism to the integral equation is
`x=X-1/12`, `y=Y-x/2`. The checker verifies the invariant identities and each
transported point exactly. Its 22 finite-quotient columns are independent
modulo 2, and a separate good-prime witness excludes rational 2-torsion.
Any integral relation must therefore have all coefficients even; dividing
the relation and repeating gives infinite descent. Thus the 22 points are
independent over the integers. This proves a lower bound only.

The displayed discriminant is positive, with complete factorization

```text
Delta = 2^12 * 3^9 * 5^6 * 13^4 * 17^4 * 19^2 * 29^3 * 71^2
        * 2465779087453622652131442949
        * 519784438179112504122441050306814600881.
```

All ten prime factors have exact recursive Lucas certificates, using 39
prime nodes in total. For each non-base node p the certificate factors p-1
completely into previously proved primes q and supplies a witness a with
`a^(p-1)=1 mod p` and `gcd(a^((p-1)/q)-1,p)=1` for every q. For any prime
divisor l of p these conditions force a to have order p-1 modulo l, hence
l>=p and p is prime. The base case is 2. The checker uses integer powers,
gcds and products, never a probable-prime test.

At every bad prime other than 17, c4 is a unit. The integral model is then
minimal and has multiplicative reduction, with conductor exponent 1.
At 17, `v17(Delta)=4<12` proves minimality and `v17(c4)=2` implies additive
reduction. Since 17>=5 its conductor exponent is 2. All remaining primes
have good reduction by the complete discriminant factorization. Therefore
the model is globally minimal and the exact conductor is `17*rad(Delta)`.
Here the extra factor 17 is essential; a discriminant radical alone would
give the wrong conductor.

These are the standard local reduction criteria used by
[Sage's local-data implementation](https://doc.sagemath.org/html/en/reference/arithmetic_curves/sage/schemes/elliptic_curves/ell_local_data.html).
An independent Sage/Tate computation of every local exponent and minimal
discriminant valuation also passed build and replay, retained in
[`next12_rank22_exact_conductor_v1.json`](../../artifacts/generated-results/elliptic-curves/next12_rank22_exact_conductor_v1.json).
The complete residual factorization was separately bounded to 120 seconds
and 1.5 GiB for the fixed 220-bit cofactor. Lucas-proof generation was capped
at 180 seconds and 1 GiB, for p-1 values of at most 129 bits.

## The fixed next-twelve experiment

This curve was one of twelve preselected addresses, positions 5–16 of the
already frozen full-prime height-1024 ranking for family 05. Choosing this
family was adaptive: its previous four addresses all gave certified lower
bounds 24–25. There was no new population scan, score refit, catalogue
prefilter or replacement address. Every fixed address received a point
attempt using the 43 generic parity classes and explicit 384-bit height
proposals. Point admission used the quotient-only reduction cache.

| Parameter | Certified rank lower bound |
|---|---:|
| `3/17` | 22 |
| `-109/521` | 24 |
| `226/513` | 16 |
| `-107/999` | 16 |
| `93/172` | 20 |
| `-553/423` | 17 |
| `643/7` | 23 |
| `917/74` | 18 |
| `611/438` | 19 |
| `-279/730` | 22 |
| `875/708` | 17 |
| `-436/997` | 18 |

All 516 retained chart/admission records passed exact replay. All twelve
complete point clouds, totalling 5,940 point occurrences up to sign within
each cloud, were independently rechecked through prime 997; no further
admission gain was found. No chart is claimed to exhaust its entire
height-100,000 box. Workers were limited to four concurrent processes,
300 seconds and 1.5 GiB each. The original `3/17` worker stopped after 35
charts; a separately declared 180-second continuation completed its eight
remaining charts. The original checkpoint and censored outcome are retained.

All twelve equations are absent from the old pinned 584-row catalogue and
263 previously retained measured equations, and are mutually distinct over
Q. The four new lower bounds at least 22 extend the consolidated inventory
to **36 distinct curves: 3 at least 25, 8 at least 24, 11 at least 23 and
14 at least 22**, with each curve counted once at its certified lower bound.
The refreshed 586-row comparison finds no match for any of the 36.
The [index](../../artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v2.json)
preserves the old 32 IDs and adds IDs 33–36. Its own frozen comparison remains
the earlier 584-row snapshot; the later comparison is retained separately.

## Replay

From the repository root, or from the extracted evidence ZIP:

```bash
python3 elliptic-curves/cas/certify_small_conductor_curve.py --check artifacts/generated-results/elliptic-curves/small_conductor_rank22_proof_v1.json
python3 elliptic-curves/cas/certify_prospective_mw16_next12_results.py --check artifacts/generated-results/elliptic-curves/prospective_mw16_next12_results_v1.json
python3 elliptic-curves/cas/export_new_high_rank_curve_index_v2.py --check artifacts/generated-results/elliptic-curves/new_high_rank_curve_index_v2.json
```

The first checker also rejects six deliberate corruptions: a Lucas witness,
a discriminant exponent, the conductor, a point, the rank label and the
catalogue placement. These checks complement, rather than replace, the
mathematical argument above. No external submission has been made.

The subsequent [bounded follow-up](SMALL_CONDUCTOR_FOLLOWUP_2026-09-05.md)
retains 127 further chart attempts and a 7,753-point finite-quotient audit,
still at lower bound 22. A short descent attempt did not reach an upper
bound. It also identifies and measures a separate observation-validation
cost, with an optional exact cache that preserves the original proof records.
