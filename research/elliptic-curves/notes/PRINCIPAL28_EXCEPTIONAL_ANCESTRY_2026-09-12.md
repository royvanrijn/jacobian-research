# Principal control: 11952 at110314/102227, rank at least28

## Result

The principal retrospective control is now the certified28-point fibre,
not the broad-search rank-at-least25 endpoint. The stronger comparison is
**M17 inside D31 versus M17 inside D28**, with known exceptional dimensions
14 and11. Neither full fibre rank is asserted exact.

| Matched ancestry measurement | Curve302 on X1092 at0 | 11952 on X948 at110314/102227 |
|---|---:|---:|
| Independent extra directions | 14 | 11 |
| Global minimum multisection degree, each point | 2 | 2 |
| Best genus constructed, each point | 1 | 1 |
| Global minimum genus | UNKNOWN | UNKNOWN |
| Exact covers checked | 504 | 396 |
| Shared tested quadratic covers between targets | 0 | 0 |
| Genus-one construction pencils used | 1 | 1 |
| Distinct selected genus-one carriers | 14 | 11 |
| Independent selected quadratic radicals over Qbar(t) | 14 | 11 |

The [complete report and per-point index](../../artifacts/generated-results/elliptic-curves/rank_ancestry_principal28_v1/SUMMARY.md)
link every original point, exact transport, covering equation, branch divisor,
lift and map. All900 maps replay, together with900 altered-map negative
controls and fresh finite-independence certificates for31 and28 points.
The504 original302 cover payloads are unchanged.

## Interpretation: similar carriers, not yet similar rank decompositions

Both jumps have the same available low-degree ancestry construction: fibres
of an existing transverse genus-one pencil. Curve302 uses class1, while11952
uses the direct old-R17 pencil; neither requires a new fibration. There is no
qualitative distinction in the extra three302 directions in this atlas.

But one pencil is **not** one common quadratic cover. The selected carriers
have14/11 different pencil parameters. Their exact trace sections agree
within each pencil, but their squarefree branch quartics have pairwise gcd1:
91 checks for302 and55 for11952. Thus their geometric branch supports are
disjoint. The [branch certificate](../../artifacts/generated-results/elliptic-curves/rank_ancestry_principal28_v1/branch-analysis.json)
records each gcd and exact trace.

Valuation at a branch point unique to each quartic proves that the14/11
quadratic radicals are independent even over Qbar(t). Their composita
therefore have degrees16,384 and2,048. For n selected quartics, the4n branch
points have inertia2, giving genus1+2^n(n-1) by Riemann–Hurwitz:212,993 and20,481.
These are minimum common-cover degrees for **these emitted quadratic
extensions only**, not obstructions to all other ancestries of the points.

Accordingly, outcome A holds only in the weak sense of similar existing
low-degree pencils. It does not establish a general mechanism for extreme
specialization jumps. No rich shared-cover rank decomposition, full
Neron–Severi orbit classification, arbitrary-base cover equivalence, or
exceptional-basis-invariant grouping is proved. This is not outcome C:
simple ancestry does exist. A stronger A/B distinction remains open.

## Minimum degree and unresolved genus

A degree-one horizontal curve is a generic section. Since the arithmetic
generic rank is17, a nonzero integer multiple of that section is an integral
combination of the generic basis up to torsion. Specialization at the smooth
target fibre would put its point in the rational span of M17, contrary to
the independently certified extra direction. Hence degree1 is impossible;
the explicit degree2 maps attain minimum degree.

Genus1 is only the minimum **in the frozen atlas**. No exhaustive genus0
incidence calculation was undertaken. The earlier fifteen-orbit rational
bisection exclusions for302 do not close its entire40,917-orbit atlas and
are not silently transferred to the new11952 control.

The atlas comprises36 constructions per target: vertical x, constant-slope
chords through all34 signed generic sections, and one retained alternate
fibration. Per-target genus histograms are identical within a fibre:
302 has(1,30,5,0), and11952 has(1,10,21,4), at genera(1,3,5,7).
These higher-genus counts reflect the chosen presentation and atlas, not
independent intrinsic rank invariants.

The11952 roster uses the
[public28 reproduction](../../artifacts/generated-results/elliptic-curves/inventory188_public28_reproduction_v1.json):
the first27 old subgroup points plus union column53 (zero-based). The first17
match the specialized generic sections in order, with recorded signs. Thus
E1–E10 are the old extra directions and E11 the public addition. This does
not upgrade the local-search endpoint from27 to28. The control remains
`historical_external` / `retrospective_known_high_rank`, outside all frozen
prospective-population selections and correlations.

## One bounded relation-sieve extension: gate still closed

The [skew-sieve certificate](../../artifacts/generated-results/elliptic-curves/two_class_skew_sieve_v1/verified.json)
keeps the two prior equation inputs, maximal orders and small factor bases.
It chooses a dyadic rectangle minimizing the exact coefficient triangle
bound, subject to A*B=2^18, and uses modular-root progressions to divide out
all supported prime powers. A separate primorial-GCD implementation replays
every primitive candidate; an independent Sage matrix checks mod2 rank.
Both fields completed within their frozen180-second caps.

| Fibre | Primitive candidates | Noncanonical relations | Deficiency | Median norm-numerator bits, old→new |
|---|---:|---:|---:|---:|
| 921/653 | 524,289 | 0 | 1,905 | 140→122 |
| 110314/102227 | 325,173 | 1 | 1,890 | 135→131 |

The sole historical relation is again(a,b)=(0,1), so neither field gains a
new independent relation beyond the first pilot. Neither meets the frozen
engineering gate of100 independent noncanonical rows and removal of at
least5% of the initial deficiency. **Stop here; no wider class-group profile
is launched.** Different primitive candidate counts and a larger/skewed box
mean this is not an equal-candidate causal comparison of sieve performance.

The mod2-only presentation, generation caveat and skew selection are guided
by [Klagsbrun–Sherman–Weigandt, sections4–5](https://arxiv.org/abs/1606.07178).
This small prototype is not their production NFS or certified Julia
reduction. Its deficiency bounds only the chosen factor-base image;
generation and missing relations are unresolved, so no full g estimate or
upper bound follows. The existing consequences g(921/653)>=16 and
g(110314/102227)>=18 still depend on known MW lower bounds.

No BNF, broad point search, new local-arithmetic census or correlation was
run. Original census, historical arithmetic and first-pilot bytes are
unchanged. The921/653 fibre is retained as a secondary control.

## Replay and navigation

```sh
timeout 300 sage -python research/elliptic-curves/cas/verify_principal28_ancestry.sage
timeout 60 sage -python research/elliptic-curves/cas/analyze_principal28_covers.sage
python3 research/elliptic-curves/cas/report_principal28_ancestry.py --check
PYTHONPATH=research/elliptic-curves/cas python3 -m unittest discover \
  -s research/elliptic-curves/tests -p 'test_two_class_*.py' -v
```

The construction freeze is preserved, including its original checker hash.
The final replay separately binds a checker-only correction that normalizes
Python tuples to their serialized JSON arrays before certificate equality;
all construction sources and inputs still match the pre-run freeze.

The earlier
[14-versus8 ancestry/strict-descent study](RANK_TRIANGLE_ANCESTRY_AND_DESCENT_2026-09-12.md),
[wide arithmetic census](WIDE_ARITHMETIC_PROFILE_2026-09-12.md), and
[censoring/first relation pilot](TWO_CLASS_RELATION_PILOT_2026-09-12.md)
remain intact as secondary evidence. The open research target is a
non-tautological common-cover mechanism or an independently measured
rank-arithmetic constraint, not further correlation with our own chart
choices. No additional campaign is scheduled.
