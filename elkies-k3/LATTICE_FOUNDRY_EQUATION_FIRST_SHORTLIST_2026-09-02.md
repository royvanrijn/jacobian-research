# Equation-first shortlist after the complete degree-three census

The strict source-rank ordering is useful, but it is not a sufficient proxy
for equation time.  On the five rootless MW17 targets whose degree-three
translation spectra are now complete, a slightly higher-rank source often has
a dramatically cheaper complete Mordell--Weil basis than the MW1 leader.

The exact aggregation is
[`../artifacts/generated-results/elkies-k3-lattice-foundry-equation-first-degree3-top5-v1.json`](../artifacts/generated-results/elkies-k3-lattice-foundry-equation-first-degree3-top5-v1.json),
with SHA-256
`2b45b5f9678da082f96c1b005f533654d5c5178109a24476083bc6f101c12836`.
It retains the rank-first leader, separate pole leaders, and all 27
nondominated source-metric rows rather than forcing incomparable source
fibrations into one weighted score.

## The ideal-cut leaders

Imposing

```text
source MW <= 2,
semistable root support,
at most three reducible-fibre supports,
exact complete MW-basis pole audit
```

selects the following same-surface sources.  The last two columns are complete
counts over all `3^17 = 129,140,163` translation cosets, not extrapolations
from the former 256-coset sample.

| target | source root type | source MW | basis pole profile | rational trisections | genus-one trisections |
|---|---:|---:|---:|---:|---:|
| `NS0028-F005` | `A2+A6+A7` | 2 | `[0,0]` | 19,645,256 | 34,294,400 |
| `NS0011-F002` | `A1+2A7` | 2 | `[1,1]` | 19,023,996 | 33,978,764 |
| `NS0022-F011` | `A1+2A7` | 2 | `[1,1]` | 18,774,826 | 33,788,528 |
| `NS0005-F008` | `A1+2A7` | 2 | `[0,1]` | 18,446,258 | 33,705,930 |
| `NS0001-F001` | `A2+A6+A7` | 2 | `[1,1]` | 18,024,616 | 33,484,468 |

Thus the complete degree-three experiment and the equation-first lattice cut
point to the same first batch.  In particular, `NS0028` is not merely rich in
minimal target vectors: it is the trisection leader and its abstract source
has two pole-zero MW generators.

## Arithmetic warning

This does **not** make `NS0028` equation-ready.  Its normalized `GF(5)` fibre
chart has 25 squarefree models and neither required section in either twist
class.  The normalized `GF(7)` chart has 112 squarefree models: the square
twist contains ten sections of the second type, while the nonsquare twist
contains ten of the first type and two of the second, but no fibre model
carries both.  These are exact exhaustive statements for the displayed
charts only; they are not characteristic-zero obstructions.

Accordingly the main unresolved discriminator is now the rational source
marking (and then the neighbour corridor), not the lattice source rank or the
degree-three richness.  The aggregation initially had no marking audit for
the other four ideal-cut leaders; the NS0005 test below closes their first
finite-field gate but not the characteristic-zero marking problem.

## NS0005 equation gate

The first new marking test targets `NS0005-S001`, the `A1+2A7/MW2` source
with basis pole profile `[0,1]`.  The generalized exact Hermite scanner
exhausts all `5^8` normalized degree-eight `A` polynomials for the
`I2+I8+I8` profile.  It finds 530 compatible signed branches, 98 branches
with the exact three discriminant orders, and 71 squarefree
`I2+I8+I8+6I1` models.  The fibre artifact has SHA-256
`efa7f439e96d397c54e167b9684bbd8431b079b696829e0c8a7f2751a1f6755e`.

The pole-zero generator has depth one at the infinity `I8` support and lies
on the identity components at the two finite supports.  Exhausting all
`71*5^4=44,375` component-adapted X polynomials in each twist class finds two
signed sections on one square-twist model and none in the nonsquare twist.
The two artifacts have SHA-256
`ec530f756c58a21ebe985646e1cce5849b70879d22e23ad5825e0a50ba18c11b`
and
`da9657569ee557f3b6e2db93e3d5d403ef4ba0674e09d8d7255713e5fdb931c2`.

On that surviving model, write the pole-one generator as
`x=N/C^2,y=M/C^3`.  The required component depths `(1,2,2)` leave three
normalized linear denominators and a two-dimensional numerator-X chart for
each.  The complete 75-case scan finds four correctly component-marked
pole-one sections.  Nevertheless none forms the required MW2 basis with the
pole-zero section: the eight smooth pair intersections have degrees
`0,0,1,1,4,4,5,5`, while the height Gram requires degree three.  This exact
GF(5) pair artifact has SHA-256
`22a0b04377d23099f1f4f3852401a81df6757e27136a198845e9af6f077c83c4`.

This is a negative result only for the displayed normalized GF(5) chart.  It
is still a stronger equation precursor than the NS0028 chart: both NS0005
generator types occur on the same fibre model, though not with the required
mutual height.

The complete `GF(7)` repetition strengthens that diagnosis.  Among all
`7^8` normalized `A` polynomials, 2,190 signed branches are Hermite-compatible,
344 have the exact three fibre orders, and 271 are squarefree.  The pole-zero
scan finds eight signed sections on four square-twist models and four signed
sections on one nonsquare-twist model.  On the four square-twist survivors,
the complete pole-one chart finds 14 correctly component-marked sections, but
the 28 smooth pair intersections again have only degrees `0,1,4,5` (with
multiplicities `6,8,8,6`) and never the required degree three.  The nonsquare
survivor has no pole-one section.  Thus both complete normalized `GF(5)` and
`GF(7)` charts are empty at the full MW2 marking gate, despite carrying the
individual generators.

The five `GF(7)` artifacts—fibre census, square/nonsquare pole-zero scans, and
square/nonsquare pair scans—have SHA-256 respectively
`76039e21549b61ad3823c12e1052f4398c4f0c8e020bf406ef029b290e46010e`,
`e065801a3a4d009567235c22304021220c85e408f368e4a1557195d51cc1acbe`,
`ca5db9faef6a02aa426899a241f7d69eede0a2d47a58f940613933858639141c`,
`50b2ba171b7f4c7c4292b161945f65d5ebb4348241be1f1e2c5ea28830979420`,
and
`bfdbd1b12af1d27287229d5a5d5b49c4aa37788cfdf36d66a1c92607048d61c4`.

## NS0031: the first equation-level positive

The same-surface search should not stop at the five targets selected before
the complete degree-three census.  `NS0031-S001` is another semistable
`A1+2A7/MW2` source, with exact complete-basis pole profile `[0,1]`, attached
to thirteen rootless MW17 frames, fifteen MW16 frames, and one MW15 frame.  Its source height
Gram is `[[2,1],[1,41/8]]`.  Although its normalized fibre orders coincide
with the NS0005 stratum, its minimum generators have different component
depths and require smooth intersection degree two, so the NS0005 negative
pair result does not apply.

The exact source-pinned scan exhausts both generator charts in both twist
classes over `GF(5)` and `GF(7)`.  Both `GF(5)` charts are empty at the full
marking gate.  Over `GF(7)`, the nonsquare chart contains 36 pole-zero
sections on 15 models but no complementary pole-one section.  The square
chart is positive: among 8 models with 20 pole-zero sections, it finds 6
correctly component-marked pole-one sections and **two complete MW2 marked
pairs**.  Both occur on normalized fibre model 157 and have the required
smooth intersection degree two and Shioda pairing one.

The four exact artifacts have SHA-256 respectively
`7e795c8573b8f805cc31b64d7ed15de52d752c79a30362eb961d4e945ee295ab`,
`cbf103e50030683aea893be5eb70fe0326b9e52a7efa4bcba1c49feab084ff35`,
`a1e0debd606102305bdd399f79e3aab3b5e442221eea5de33042991528b6c0e6`,
and
`552a8d0042d86ef5568f1254ea5d986d3730f444768b16d60366dec1953b83eb`
for the `GF(5)` square/nonsquare and `GF(7)` square/nonsquare charts.

This is the strongest current equation precursor, not a characteristic-zero
construction.  A finite-field marked pair does not prove a lift over
`Q`, rational source marking, rational parameterization, or cheap neighbour
corridor.  It does change the work order: lift the marked `GF(7)` model and
enumerate the degree-three spectra of its thirteen MW17 target frames before
spending more equation effort on the negative NS0028/NS0005 normalized
charts.

There is also a positive infinitesimal lifting gate.  Writing the normalized
short-Weierstrass coefficients, both section numerators, and the monic pole
denominator as 52 variables gives 59 fibre, section, and exact component-jet
equations.  At model 157 their Jacobian over `GF(7)` has rank 51, so the full
marked tangent space is one-dimensional; an explicit maximal minor is
`1 mod 7`.  Deterministic Newton corrections are consistent for the full
overdetermined system through seven further digits, producing 52 explicit
integer coordinates on which all 59 residuals vanish modulo
`7^8 = 5,764,801`.  The byte-checked tangent and finite-lift certificate has
SHA-256
`0187b44d46d41572961261be344a1bafd3ee204f98754d173f15438699863f55`.
This finite-precision lift is not a proof of an infinite compatible `Z_7`
point.  Extending the lift and certifying localized equation dependence—or
otherwise proving formal smoothness of the overdetermined ideal—remains the
next arithmetic gate.

The thirteen rootless NS0031 frames were then prescreened together.  Their
exact rational-bisection counts occupy the narrow interval 41,885--41,959,
while the sampled degree-three coordinates vary enough to select five frames
for full enumeration.  The pilot artifact has SHA-256
`4f9a8041cf2f8d56d49d6ede1b5c2529c9edda414517ac93fbc29948a7d9f34b`.
It is used only for selection, not for the final score.

The complete census visits all `3^17 = 129,140,163` translation cosets for
each selected frame.  It changes the pilot ordering—`F006` rises from fourth
to second and `F018` falls from second to fifth—and gives:

| NS0031 target | rational trisections | genus-one trisections | maximum coset minimum norm |
|---|---:|---:|---:|
| `F017` | 20,059,924 | 34,485,512 | 26 |
| `F006` | 20,059,202 | 34,487,054 | 26 |
| `F008` | 20,059,160 | 34,486,376 | 26 |
| `F026` | 20,053,748 | 34,489,024 | 26 |
| `F018` | 20,048,322 | 34,488,086 | 26 |

The complete artifact has SHA-256
`31c6ee0b0ee4669295072754965fdc8b560994fe585167ca7d28540cffc74c0e`.
In particular, `NS0031-F017` exceeds the former audited leader `NS0028-F005`
by 414,668 rational and 191,112 genus-one trisection cosets.  NS0031 is now
the strongest current candidate on both independent axes: a low-pole MW2
source with a positive equation-level marking, and unusually rich complete
degree-three target spectra.

The central source scorer now consumes both the NS0031 pilot and complete
census.  Its declared low-degree tie-break remains bisection-first, so it
selects `F018` as NS0031's best audited multisection target because `F018`
ties for the largest rational-bisection count; `F017` is the distinct
trisection-first leader.  Keeping both labels avoids letting a sampled or
single-degree choice silently stand in for the full multisection objective.

## Determinant 720: a better source, but not a richer target

The source-first search on the independently manufactured Golay-octad
determinant-720 K3 now supplies the cleanest equation precursor in the
foundry.  The prescribed-root enumeration is complete across all 23 rooted
Niemeier ambients in its declared rank-15/16, one-to-three-support window.  It
returns 4,823 distinct reduced-Gram rows: 32 MW1 and 4,791 MW2.  The separate
complete Smith-quotient pole audit finds 1,587 rows with a full physical MW
basis of pole at most two.  Among them are 177 semistable all-A rows; 35 have
fibre type `3A5` and basis pole profile `[0,0]`.

The compact ledger is
[`../artifacts/generated-results/elkies-k3-golay-det720-equation-first-shortlist-v1.json`](../artifacts/generated-results/elkies-k3-golay-det720-equation-first-shortlist-v1.json),
SHA-256
`c4930c78a67c2e144e8fdb09cb6325e76bfb2933927be10444e11511834439d0`.
It retains all 177 semistable complete-basis rows while hash-pinning the
153 MB source inventory and 29 MB pole audit.  The latter artifacts have
SHA-256 respectively
`f203027cae98df3f1cf69286dd149a73ef93f35ee9408130d89e39a856dbb7af`
and
`ddeda8fd0059ec7b1e45f09d87be1a0dd6815152ee12bcd2bb3fb313e6e46bf8`.
Distinct reduced Grams are not asserted in general to be distinct
integral-isometry classes.  On the ideal cut (MW2, semistable, at most three
supports, complete pole profile `[0,0]`), however, all 48 rows have now been
classified exactly.  They form only three marked integral-isometry classes:
35 copies of `3A5` represented by `G720-S0128`, four copies of `A11+A4`
represented by `G720-S0260`, and nine copies of `A3+A4+A8` represented by
`G720-S0052`.  Thus the three finite-field charts already tested exhaust the
entire ideal cut.  The classification artifact has SHA-256
`34f0e94abea6b257f56c161a14d14737ab14cd29855f0f30148249023d9b2674`.

The leading row is `G720-S0128`, a semistable `3A5/MW2` source with height
Gram `diag(5/6,4)`, three reducible supports, trivial torsion, and two
pole-zero physical generators.  Its normalized `GF(7)` nonsquare chart has
237 squarefree `3I6+6I1` models and 24 complete marked MW2 pairs, whereas the
other three tested prime/twist charts are empty at the pair gate.  The
positive pair artifact has SHA-256
`fe31c9ef763a4d1d52feb0eea6cbbecd846efd20f45153a49932f6c767f9a5d1`.
The comparison charts `A11+A4` (`G720-S0260`) and `A3+A4+A8`
(`G720-S0052`) carry individual polynomial sections but no complete marked
pair in either twist class at 5 or 7.

At the first marked `G720-S0128` point, the normalized system has 46
variables and 55 displayed equations.  Its Jacobian has rank 45, tangent
dimension one, and a maximal minor equal to 6 mod 7; all equations lift
simultaneously through `7^8`.  More importantly, the apparent overdetermination
is removed exactly.  For

```text
f(X)=X^3+A*X+B,  D=4*A^3+27*B^2,  N=2*A*X+3*B,
8*A^3*f(X) = D*(N-B) - 9*B*N^2 + N^3.
```

Because `A` is a unit at all three marked supports, the `I6` and component
depths `(1,1,3)` force section-residual orders `(2,2,6)` at
`0,1,infinity`.  Thus the residual is
`t^2*(t-1)^2*(c0+c1*t+c2*t^2)`; the three retained residual coefficients
2, 3, and 4 kill `c0,c1,c2`.  The ten omitted equations are automatic, and
the unit 45-by-45 minor proves a one-parameter formal `Z_7` marked family.
The finite-lift and formal-smoothness certificates have SHA-256 respectively
`4cac494b0ace3c8bedc2c7515f4ca7409b53cf155298493849a5a0dacd863f53`
and
`f878f7b5b453b1b22db1446cfa7f2517544fcfb64616d53b3222c6bb1414d9e1`.

Fixing the free coordinate to `s6=10` changes the arithmetic conclusion.  All
46 coordinates reconstruct over `Q` at precision `7^40`, and direct
substitution makes all 55 equations vanish exactly.  The resulting model has
split fibres `3I6+6I1`; its two rational sections have component depths
`(1,1,3)` and `(0,0,0)` and height Gram `diag(5/6,4)`.  Hence it displays a
rank-19 NS sublattice of determinant `-720`.  The fixed lift and exact-QQ
artifacts have SHA-256 respectively
`5c8c900a8fa6e35a6647316bc5059aad1e9252fbdab58ed01b09142b759a8d9c`
and
`63b8a3e97d6e1e41a722283cbfb670ff32e2a4a4da129fdaf14018cda83dfb78`.

Exact point counts at the good primes 17 and 19 give rank-20 reductions with
Artin--Tate discriminant representatives `-256` and `-1440`.  Their ratio
`45/8` is nonsquare, so the two-prime specialization argument proves that the
rational model has geometric Picard rank exactly 19.  The Picard certificate
has SHA-256
`82a46940e54c9be162371d689dfd518ca1b9835b13a24b867c33696e99e6c4f3`.
Saturation of the displayed section/frame lattice—and hence exact identity
with the abstract determinant-720 target NS class—remains a separate gate.

The target-side result is deliberately less flattering.  The rootless MW17
frame `G720-F001` has 3,064 norm-four vectors, 34,848 rational bisection
orbits, and 64,355 genus-one bisection candidates.  Its complete census over
all `3^17 = 129,140,163` translation cosets gives:

| target | rational trisections | genus-one trisections | maximum coset minimum norm |
|---|---:|---:|---:|
| `G720-F001` | 15,717,830 | 31,988,690 | 24 |

The complete artifact has SHA-256
`a9960fa9fdb375b25a6b9f79d919846a56c7e4aa85235a466748e16448f463ed`.
These counts are substantially below the NS0031 leaders despite the simpler
source.  Source equation cost and target multisection richness are therefore
independent Pareto coordinates: determinant 720 is presently the source-first
leader, while NS0031 remains the complete degree-three leader.  A physical
neighbour corridor and rational algebraization are the next discriminators.

The entire degree-two corridor gate is also complete.  The exact target
coset-minimum histogram is `0:1, 4:1532, 6:30176, 8:64355, 10:34848,
12:160`.  Section nonnegativity and isotropic integrality retain the 64,355
minimum-eight and 160 minimum-twelve classes.  Of these 64,515 classes,
64,512 have fibre divisibility one and three have divisibility two.  Splitting
every primitive class produces `rA1` children for `1 <= r <= 11` at minimum
eight and 160 rootless children at minimum twelve; no child has even the
rank/count/determinant signature of `3A5`.  Hence no one-edge degree-two route
reaches `G720-S0128`.  The exact negative
artifact has SHA-256
`5a83518072cb5f4d010e911a6c99de28c9aaecd6c777b3154bcfacccd1a25ad6`.
This does not exclude higher-degree edges or a short multi-edge route.

## Reproduction and proof boundary

Rebuild and byte-check the aggregation with

```bash
python3 elkies-k3/scripts/build_lattice_foundry_equation_first_shortlist.py
python3 elkies-k3/scripts/build_lattice_foundry_equation_first_shortlist.py --check
```

The NS0005 finite-field gates are reproduced by

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0028_source_ansatz_modp.sage \
  --ns-id NS0005 --source-id NS0005-S001 --prime 5 --examples 0 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-source-ansatz-mod5-v1.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_lattice_foundry_ns0028_pole0_section_pairs_modp.sage \
  --mode infinity-pole0 \
  --input artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-source-ansatz-mod5-v1.json \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-infinity-pole0-sections-mod5-v1.json

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_lattice_foundry_ns0005_pole1_pole0_pairs_modp.sage
```

The NS0031 marking gate is reproduced and byte-checked by

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_lattice_foundry_ns0031_a1_2a7_marking_modp.sage \
  --fibres artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-a1-2a7-source-ansatz-mod7-v1.json \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0031-a1-2a7-marking-mod7-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_lattice_foundry_ns0031_marked_gf7_hensel.sage \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/complete_lattice_foundry_degree3_spectrum.py \
  --frame-id NS0031-F017 --frame-id NS0031-F018 \
  --frame-id NS0031-F008 --frame-id NS0031-F006 \
  --frame-id NS0031-F026 --workers 8 --chunk-size 1000000 \
  --float-type dd --audit-precision 256 --audit-stride 4096 \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-degree3-complete-ns0031-pilot-top5-v1.json \
  --check
```

The determinant-720 source, formal-family, and target-spectrum summaries are
rebuilt and checked by

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_golay_det720_3a5_marked_gf7_lift.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_golay_det720_3a5_formal_smoothness.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_golay_det720_3a5_marked_gf7_lift.sage \
  --free-parameter-integer 10 --lift-precision 40 \
  --output artifacts/generated-results/elkies-k3-golay-det720-3a5-s6-10-lift-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_golay_det720_3a5_source_qq.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_golay_det720_3a5_picard19.sage --check

python3 elkies-k3/scripts/build_golay_det720_foundry_adapter.py --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_golay_det720_degree2_source_corridor.sage \
  --workers 16 --chunk-size 512 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/complete_lattice_foundry_degree3_spectrum.py \
  --database artifacts/generated-results/elkies-k3-golay-det720-foundry-adapter-v1.json \
  --frame-id G720-F001 --workers 16 --chunk-size 1000000 \
  --float-type dd --audit-precision 256 --audit-stride 4096 \
  --output artifacts/generated-results/elkies-k3-golay-det720-degree3-complete-v1.json \
  --check

python3 elkies-k3/scripts/build_golay_det720_equation_first_shortlist.py --check
```

All displayed lattice-source data, section-pole profiles, finite-field marking
counts, formal-local statements, rational equation identities, Picard bound,
and degree-three counts are copied from hash-pinned exact artifacts.  The
report constructs no saturation/target-NS identity, rational parameterization,
neighbour route, effective target multisection, or specialization rank jump.
