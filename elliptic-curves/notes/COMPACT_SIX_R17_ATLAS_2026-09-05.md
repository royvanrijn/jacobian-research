# Compact coordinates for six existing R17 fibrations

All six compiled R17 presentations now have exact equations with maximum
coefficient numerator/denominator sizes 141–169 bits. All 102 generic sections
are transported and checked as rational-function identities. These are new
coordinates on existing fibrations, not new generic-rank records or newly
discovered specialized curves.

| Native chart | Literal coefficient bits | After constant scaling | After base compactification |
|---|---:|---:|---:|
| `103b2` (representing `0e80b`) | 991 | 786 | 141 |
| `11952` | 1371 | 1141 | 142 |
| `074d9` | 1584 | 1207 | 141 |
| `07ca9` | 1621 | 1289 | 141 |
| `08234` | 1603 | 1237 | 169 |
| `08f72` | 1887 | 1503 | 146 |

The [portable atlas](../../artifacts/generated-results/elliptic-curves/compact_six_r17_atlas_v1.json)
contains every coefficient, section, base matrix, constant scale and source
binding. The native generic height forms are retained in their exact basis;
each has determinant 948. For `074d9`, the recorded unimodular word matrix
transports the stored Gram form back to the literal section basis.

## Exact identity

If the retained matrix is `(a,b;c,d)` and the constant scale is `u`, the old
base parameter and elliptic coordinates are

```text
t_old = (a*t+b)/(c*t+d),
x_old = u^2*x_new/(c*t+d)^4,
y_old = u^3*y_new/(c*t+d)^6.
```

The checker proves `ad-bc != 0`, `u != 0` and the full coefficient identities

```text
(c*t+d)^8 A_old(t_old)  = u^4 A_new(t),
(c*t+d)^12 B_old(t_old) = u^6 B_new(t).
```

It reconstructs each of the seventeen old sections, checks its old equation,
transports it by the displayed formula and checks the new equation. The
degree-one base map and Weierstrass isomorphism preserve generic independence.
No specialized point search, numerical height estimate or analytic rank enters
this certificate. Replay with:

```sh
sage -python elliptic-curves/cas/export_compact_r17_atlas.sage --check \
  artifacts/generated-results/elliptic-curves/compact_six_r17_atlas_v1.json
```

## The missing coordinate reduction

Constant scaling alone left the six models large. A bounded Cremona–Stoll
reduction of the auxiliary curve `y^2=A(t)` reduced the first model only from
786 to 778 bits; a separate greedy small-prime base descent reached 772.
Prime-local auxiliary minimization using only small primes reached 708.
Each failed its declared 25-percent improvement gate, so it was not expanded
to the other five families.

The useful diagnostic was the gcd of the discriminants of primitive `A` and
`B` and their resultant. For `103b2`, after stripping primes below 10000,
the remaining gcd is the exact 56th power of the proved prime
`28935871643137959223363528111`. Passing that prime explicitly to PARI's
prime-local auxiliary minimization exposed the large base-coordinate scale.
The resulting rational base map, followed by Cremona–Stoll reduction and
exact weighted constant scaling, gave 141 bits.

The other families were then admitted by the frozen improvement gate.
Four completed under the original factor rule. The two remaining composite
power roots, of 141 and 157 bits, exceeded the original 128-bit factoring
gate; a separate continuation permitted just those two up to 160 bits.
Every worker retained the 120-second and 1-GiB cap. All six completed.
Generation used auxiliary hyperelliptic models to propose maps; the final
elliptic identities and section checks are the proof. There is no assertion
that the resulting elliptic surfaces are globally minimal at every prime.

The failed attempts, protocols, factorizations and worker logs are retained in
`artifacts/local/elliptic-curves/r17-*-base-v1/`,
`r17-base-bad-primes-v1/`, `r17-compactification-v1/` and
`r17-compactification-160bit-v1/`.
The [portable evidence manifest](../../artifacts/generated-results/elliptic-curves/compact_six_r17_coordinate_evidence_v1.json)
and [bundle](../../artifacts/generated-results/elliptic-curves/compact_six_r17_coordinate_evidence_v1.zip)
retain 49 protocol, diagnostic, source and log files, including the failed gates.

## Consequence for the curve search

Large coefficients in the old literal models were partly a coordinate
artifact. This supplies a practical input to a new prospective experiment;
it does not retroactively change the exposure or conclusions of old runs.
New small-parameter boxes in these coordinates are different populations
and must be frozen separately. The retained 43 generic deep masks may be
reused only after matching their generic basis and height form. Any new
specialized rank still requires exact independent point certificates and
post-selection rational-isomorphism comparison.

The first balanced pilot froze height 1024, all 562 primes from 5 through
4093 before retention, and four finalists per family. Each family scores
1,275,854 signed primitive parameters, excluding zero and infinity. The
generic mask sets match the native height forms and their representatives
check at norm 12; the source's specialization-dependent ordering and public
exceptional points are not selection inputs. New specialized numerical
height forms only schedule representatives.

The executed protocol is `compact-six-r17-h1024-v2`: the preceding v1 froze
generic masks but ran no table or point search before its dependency binding
was tightened. The executed selector and catalogue input hashes are pinned.
All six selections completed. Direct Sage cardinality checks independently
validated 240 projective residue rows at primes 5, 7, 11 and 13.

The fixed 24-address roster has five exact catalogue matches, which are
skipped without refilling; its nineteen fresh addresses receive at most
43 charts at height 100000 and four seconds per chart. Four workers run at
once, each capped at 300 seconds and 1.5 GiB, with durable per-chart witnesses.
The batch stops further dispatch on a certified rank-at-least-28 signal for
independent replay. Completed chart counts do not imply full search-box
coverage. The active ledger is retained locally under
`artifacts/local/elliptic-curves/compact-six-r17-h1024-v2/ledger.json`.
