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
`e5ab02930f87b162d00ebc18a416b03615bbe8e3febbbb10cd85198b9bf601be`.
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

The subsequent saturation gate **rejects this rational point from the
determinant-720 class**.  It has the rational 3-torsion section

```text
x = 3*t^4 + 36/5*t^3 + 258/25*t^2 + 36/5*t + 3,
y = 1728/5*t^2*(t-1)^2,
```

and the displayed height-four section `Q` is twice the rational height-one
section with

```text
x = 3*t^4 + 36/5*t^3 + 3138/25*t^2 + 36/5*t + 3,
y = 41472/25*t^3.
```

These enlarge the displayed frame by index `3*2=6`.  Exact enumeration of
the discriminant form with Smith invariants `[2,6,60]` finds no isotropic
subgroup of order greater than six, so this is the full NS lattice.  Its
determinant is `-720/6^2=-20`, its torsion is `Z/3`, and its free MW height
Gram is `diag(5/6,1)`.  Thus `s6=10` is a simple Picard-19 K3, but not the K3
carrying `G720-F001`.  The rejection artifact has SHA-256
`df50dd78ace634fb1611307cbac802c37ded715ee952c31d649a8aecd78d534e`.

The bounded rational search now covers all six distinct etale marked residue
disks found in the `GF(7)` chart.  It tests 1,478 reduced rational parameters
`a/b` with `|a|,b <= 40`, `7` not dividing `b`, reconstructs only after all 55
equations vanish exactly, and then applies rational 3-division and section
halving gates.  Exactly three rational parameters reconstruct: `s6=10` in
two disks and `s6=-8` in one.  Every one has rational 3-torsion and a rational
half of `Q`, hence belongs to the same wrong determinant-20 locus.  The six
scan artifacts have SHA-256
`60ca764fd3012ae190aa1409ca143026df329976c45a1efc72c7afa92f7b0ca0`,
`cec7abedc99a65dda82b973b7e2b25d287846bc2bc9d176de72ddcd0e1f991fa`,
`6d223be7bf9d9d85bac9aa5cff7a6349d10a1f507073d1e597efe812d004c835`,
`55de7b4a6e9ca533a0b5444f5d5ffc68cfd0d9475407ed3fe56570cc22c69d61`,
`cd22021ee950cd8116255f121f35ff1cedc141f0d715bb4a8ce519ee50c8c379`,
and
`180464017a83addf7c7995a7d084d1b8a5260c4a6bee52ba0ef36bef2d217699`.
This is an exact bounded negative, not a proof that the remaining formal
parameters are irrational or that no determinant-720 rational point exists.

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

## Determinant 500: MW1 source to a rootless MW17 target

Searching the other fibrations of expanded-catalogue surface
`K3-04b86146cc6b284b` gives the first candidate that matches the revised
objective directly.  Its transcendental lattice is
`U(5)+<20>`.  The pinned target `K3-04b86146cc6b284b-F001` is rootless MW17,
while a complete search in the six large-component ambients `D24`,
`D16+E8`, `3E8`, `A24`, `A17+E7`, and `A15+D9` returns 3,101 exact
MW1--3 source Grams in the declared rank-14--17, one-to-three-support window.
There are 47 MW1 rows, 236 MW2 rows, and 2,818 MW3 rows.

Three reduced Grams have source type `A3+A4+A9/MW1`; exact PARI isometries
merge them into one integral class.  Its invariants are

```text
reducible fibres: I4 + I5 + I10
MW group:         Z, with height Gram [5/2]
torsion:          trivial
support count:    3
minimum basis PO: 1
determinant:      200 * (5/2) = 500
```

Thus this is a one-dimensional Picard-19 lattice family with a semistable
root pattern and only one additional section condition beyond the expected
two-dimensional `I4+I5+I10` fibre stratum.  The lattice does not prove that
the remaining singular fibres are all `I1`, nor does it descend the source
marking to `Q`; those are the next equation gates.  The exact source inventory
and compact promotion certificate have SHA-256
`fd651d8daadbaa2b797e34aa949fce68bd14ea7e13952f73f5ea66699dcf2aba`
and
`59381302670d455e3906b27fc3cf58933f088422e3ae661f74bd6137be1eaed0`.

The target is not multisection-leading.  Its exact low-height degree-two
counts are 29,040 rational and 63,895 genus-one bisection translation orbits.
A deterministic 1,024-coset exact-CVP pilot finds 96 rational and 226
genus-one trisection candidates, plus 349 genus-zero-through-two
quadrisection candidates.  The pilot artifact has SHA-256
`5273b0c4631c7b2b885905f23ecd9b7d6e40b42aab914fc413c102c48f740800`.
The complete `3^17` degree-three census gives 12,095,162 rational and
29,878,240 genus-one trisection translation cosets, with maximum coset
minimum norm 24.  Its SHA-256 is
`4e5f1d7011b9f1732343865953cdd8f802a83ae585a675dc8e9bbed5065a2393`.
These counts confirm rather than overturn the pilot: determinant 500 lies on
a genuine Pareto frontier, with dramatically better source-equation cost
than NS0031 or NS0024 but substantially weaker target multisection richness.

The first equation gate is positive but characteristic-5 singular.  The
complete normalized `GF(5)` fibre scan exhausts all `5^8=390,625` choices of
`A`, finding three exact squarefree `I4+I5+I10+5I1` models.  On those three
models the pole-one section chart is especially small: depths `(2,0,5)` at
`I4,I5,I10` give seven independent X-jets, hence one X numerator for each of
the three allowed monic linear denominators.  The square twist contains four
signed marked generators on two models; the nonsquare twist is empty.  The
fibre, square-marking, and nonsquare-marking artifacts have SHA-256
`8b923fb250f5327809c01ff52e3ab47ddbbd0a68ca98c4ef3f410db52025afcc`,
`68aadb4a04e78fdd33a13742339883f2ee3c79d004f9598c4a0a516372f0090f`,
and
`ccd1a45faea89705bdf8c00435b6a98850d3e311057f1d6498c8f66dddce299e`.

Neither distinct marked model supplies a smooth 5-adic lift in the displayed
integral system.  Each 40-variable, 53-equation Jacobian has rank 38 rather
than the expected 39, and the compatible overdetermined lift stops at
`5^2`.  The two tangent artifacts have SHA-256
`0a76a92119ad0e4d22269d109d7fe7f7af2b80c6cd0b82b1cff3f2c2e2de9be1`
and
`55bb7500f02b9a52f962e342148a2276b07b29570f1c5556ea89b8e6268262f9`.
Because the fibre configuration contains `I5`, this is not promoted to a
characteristic-zero obstruction; the prime-7 chart is the next arithmetic
gate.

That prime-7 gate succeeds.  The complete `7^8=5,764,801` fibre scan finds
six squarefree normalized models.  Their square and nonsquare twists contain
respectively two and four signed pole-one sections with the required depths.
All three geometrically distinct marked seeds have rank-39 Jacobians in 40
variables and lift through `7^8`; the pinned square seed has tangent dimension
one and unit minor `6 mod 7`.  The fibre, square-marking,
nonsquare-marking, and pinned Hensel artifacts have SHA-256
`4067b0d98b7ef96ea8591d20537f3cb2211487a46d24b5dbc060782f488e4f95`,
`6395b2e836a407e4662d5ccac814900d9b8161dbc6bea18db05798804cf45dbc`,
`8ba842cfcc2eee7411c31c15c665a74a5ed358b4a3d766b873994d1eaabb407f`,
and
`d4f35d028c3fa90887e29ac4a16853493ee5a1b4899a654d1bf3dbb282fa5b74`.
The exact identity

```text
D=4*A^3+27*B^2, H=2*A*X+3*B*C^2,
8*A^3*(X^3+A*X*C^4+B*C^6)
  = D*C^4*(H-B*C^2) - 9*B*H^2*C^2 + H^3
```

closes formal dependence.  The fibre and component jets force the section
residual to vanish to order four at zero and order ten at infinity, leaving
only coefficients 4--8; these are exactly the five residual rows in the unit
minor.  Hence the 39 retained equations define a one-dimensional formally
smooth marked branch over `Z_7`.  The formal certificate has SHA-256
`c676312311c238506f2ce0ec177b5fc142e95d49ee0a784e458ff730205a7d56`.
Rational algebraization remains open.

A bounded rational probe across all three smooth residue disks tests 87 small
integer values of the free coordinate at precision `7^40`.  One value,
`m4=-20` in the first nonsquare disk, reconstructs all forty coordinates and
gives a strikingly small `Q` model.  It is nevertheless the wrong lattice:
the displayed height-`5/2` pole-one section is exactly five times the rational
pole-zero section

```text
x = 4/3*t^4 - 4/3*t^3 + 5/3*t^2 - 10/3*t + 3,
y = 16*t^3 - 16*t^2,
```

whose component depths are `(2,1,3)` and height is `1/10`.  The displayed
frame has Smith invariants `[5,5,20]`; exhaustive discriminant-form enumeration
finds six isotropic subgroups of order five and none larger.  Thus the explicit
fifth root realizes the maximal index-five enlargement, and the rank-19
primitive closure has determinant `500/5^2=20`.  The exact rational rejection
artifact has SHA-256
`2548f5a2e62e1a4e551ed2fe0592df56a398160d56e1ffcbfbfe69aca6572b2e`.
The scan itself is now a byte-checkable artifact: 86 parameters have no full
coefficientwise rational reconstruction and the sole exact reconstruction is
the rejected determinant-20 point.  Its SHA-256 is
`140910f1521559972a76396c5a2e0873f8edf6e42e977f9a7f9d9b1c4d63a513`.
This is a bounded integral-parameter result, not a rational-point exclusion.
As with determinant 720, every future rational reconstruction must pass
torsion and section-divisibility gates before it is identified with the
foundry lattice.

The same-surface MW2 alternative has also been tested rather than inferred
from its lower rank.  Nine retained `A3+A4+A8/MW2` rows have pole profile
`[0,0]`; exact integral isometries put them in one unmarked frame class.  Their
chosen physical bases split into two profiles of sizes three and six, but
both impose the same normalized equation conditions:

```text
reducible fibres:       I4 + I5 + I9
generator depths:      (1,0,1), (1,0,3)
required smooth P.Q:   1
MW rank / basis poles: 2 / [0,0]
```

This chart is fibre-rich: exhaustive scans find 30 squarefree normalized
models over `GF(5)` and 114 over `GF(7)`, versus 3 and 6 for the MW1 chart.
But it is marking-poor.  Neither twist class at either prime contains a valid
ordered MW2 basis.  The `GF(7)` nonsquare chart does contain 20 individual
sections in the first generator class and 32 in the second; on five models
both classes occur, giving 28 component-matched pairs, but all 28 pairs meet
at a reducible fibre and fail the required smooth-intersection gate.  Thus
this is a useful fallback if the MW1 rational gate stays empty, not presently
a better promoted source.  The isometry, fibre-5, fibre-7, and four marking
artifacts have SHA-256
`c06c0d03c9c74c3bb2112197665d6336543ab7f1c179bb168503ce068712cfbc`,
`130d5d755abfeb8dab673c02a7ca5854c1b4896b78bd9dc788f49cd3cd3dbf1d`,
`9f24bf1f1e0761ee2bba5b3058d82ffd83b1159fb7057deacd03899f107d47e3`,
`3e288d5ca45bed9545b56d7b85ad7b85581b55e2d7d019ab8ed2260784041a0d`,
`38d00092ee2d4952c83dcd144449c2656c27e63903912e02730fe78bd6cba6d2`,
`5c273eceb1f1506a127bc58de27dafed603e79b3a432e27e1207136b5f95c5d8`,
and
`f22c748b9743526951a479bca1a92006bcc5b27ea1559a3ebc6bc6c3791a7f59`.

The first exact same-NS corridor beam from `S0160` is also complete.  It uses
degree-two old fibres, `q=4`, section pole order at most one, rank-first beam
width eight, and depth twelve; MW enumeration is capped at 2,000 vectors only
from MW rank four onward.  The best retained root rank falls from 16 at the
source to 14, 9, 6, and 4 in the first four edges, hence reaches MW13, but it
never reaches the named rootless MW17 frame.  Root rank four recurs at every
depth except seven.  This is a bounded beam miss, not a complete neighbour-
graph rejection: other `q`, degrees, pole orders, repaired presentations, and
wider beams remain open.  The artifact has SHA-256
`aabc732debe6378eafeaac651f0f3a6c5ba181cbc8f55059a53a976499eb9bb2`.

## Determinant 384: the first explicit MW-rank tradeoff candidate

The source-first expansion to lower target rank immediately finds a much
better abstract source on `K3-6ce16abb9de3c7c5`.  Its imported target
`F001` has root type `A2` and MW rank 15.  The complete six-large-ambient,
semistable MW0--2, one-to-three-support cut contains **144** source Grams, all
of MW rank two:

| source root type | supports | rows |
|---|---:|---:|
| `A10+A5` | 2 | 69 |
| `A15` | 1 | 45 |
| `A2+A5+A8` | 3 | 30 |

Fifty rows have a complete physical MW basis through pole two.  The best cut
consists of 24 `A10+A5` rows with pole profile `[0,1]`; exact isometry testing
puts all 24 in one integral frame class and one marked-basis profile:

```text
reducible fibres:       I6 + I11
generator depths:      (0,4), (0,2)
required smooth P.Q:   1
MW rank / basis poles: 2 / [0,1]
```

This validates the proposed tradeoff sharply: an MW2 source for an MW15
target appears in abundance even though five rational-moduli rootless MW17
surfaces have no source at all in the same ideal cut.  The lattice-only
adapter, source census, and frame classifier have SHA-256
`6ea0566c456d1a431a39aabebad06569237fa5081ec7f18529c2cfc2147267ec`,
`e52d1fda05b7a799848a7b4190c7508473428cf977f22f2cff0fd95d3b9c6fe8`,
and
`5e57127ccc2cb14cf7fc2904a00d9c74cc00ebcb239241b7a839bade284d9492`.

The two-support fibre chart is also rich.  Exhaustive scans find 152
squarefree normalized `I6+I11+7I1` models over `GF(5)` and 1,032 over
`GF(7)`.  But the best chart is marking-poor.  At `GF(5)`, the square twist
has 16/0 individual sections in the two generator classes; the nonsquare
twist has 8/16 and four models with both classes, but no required pair.  At
`GF(7)`, the square twist has 108/60 sections with no common model; the
nonsquare twist has 84/228 sections and 72 component-matched pairs on twelve
models, all with the wrong smooth intersection.  Thus all four normalized
charts are empty at the full MW2 marking gate.  The fibre artifacts have
SHA-256
`cb77aa2773c31e83a57c9408f36b3b59925fdff32ebfbb2521079ed351011ecd`
and
`b3a5cf83075844e781e6ae6d84fbd4c4710c191fe58fdfd3094aaa1c22ed9c19`;
the four marking artifacts have SHA-256
`9aa5aa45f512785150a93160dbaf0fe5090d0ed7f5394a4bc47c6e5280b8a00f`,
`9a1107ed0b54fa5042d7b0dbc2262d237147833376bcf570e5a8703bc607559f`,
`64f87cc5cf4ba54c9b5d358b5cd120633a874fee0d066f778251a611b673b275`,
and
`6633682bf4b45d398176807fbbe91142bebb7c3408520906e0a4f0b315c6580a`.

The first same-NS corridor beam reaches root rank four, hence MW13, in five
degree-two edges but does not reach the exact `A2/MW15` target through depth
eight.  Its SHA-256 is
`0d757c26dc5ad4b6ee80f5ae63a6804a882ec2df6c3e1942074e427901f063b4`.
The T-arithmetic curve is still unidentified because the even-Clifford order
is non-Eichler at a ramified prime; the lattice-only adapter explicitly keeps
equation authorization false.  The aggregate candidate certificate has
SHA-256
`266e9e9e94714f224e3559babb696011dfb78c1c3e64efd5c86067a978231d39`.
Accordingly determinant 384 is the new lattice-source leader but remains
behind the formally smooth determinant-500 MW1 branch at the rational-marking
gate.

The next `A2+A5+A8` pole-`[0,1]` class has also been closed at this gate.  Its
17 rows form one integral frame class and two marking profiles with identical
depths `(0,2,1)` and `(1,2,4)`, requiring smooth intersection one or two.
The fibre scans find 82 squarefree models over `GF(5)` and 306 over `GF(7)`.
Both square twists have no pole-one generator.  The nonsquare charts contain
14/10 and 50/42 individual generators at primes 5 and 7, but all 28 and 92
component-matched pairs meet a singular fibre.  Hence neither intersection
profile survives.  The classifier, two fibre scans, and four marking scans
have SHA-256
`8a4c3ab0d035eaf0f5c7aaca38a4372879ecd6781491103eb810dd84d5f00d1a`,
`575345033900da7581d62b8ff2d5b9c40673e74739587883a127b0ace5b73795`,
`87e32f33d10d434e84a8dc2835267451364e2287681eec054bed2885f211b459`,
`ea9b282db0eb75553abff424e65f817133c5f87d3a99ef1228355b3196bdaceb`,
`fcc27b3b5b723247f50368c52e6d426d8fda974b4b4570c6703a3bf54f2a730a`,
`4750975b23a0b0bb1792231803d8fc18afbfd88d1d0fd4c4d64a221c5100746c`,
and
`31419c18c0dc863198fe66fb015f349f5a2901974afd659d0a61c2c82769b184`.

## Determinant 654: abundant MW2 sources, but the cheap markings fail

The next MW-rank tradeoff screen is `K3-14ad03cd7c1848b2`.  Its selected
target `F001` has root type `A1`, hence MW rank 16.  The same complete
six-large-ambient semistable MW0--2 cut returns **420** source Grams: 9
abstract MW1 rows and 411 MW2 rows.  The MW1 count is not yet an equation
advantage: six rows still require glue analysis and the other three have no
physical basis through pole two.  In contrast, 127 MW2 rows have a complete
basis through pole two.

The first marking audit exhausts the five cheapest distinct profiles:

| representative | fibres | basis poles | required smooth `P.Q` |
|---|---|---:|---:|
| `S0050` | `I3+I6+I9` | `[0,0]` | 1 |
| `S0093` | `I3+I6+I9` | `[0,0]` | 1 |
| `S0360` | `3I6` | `[0,0]` | 1 |
| `S0071` | `I6+I11` | `[0,1]` | 4 |
| `S0197` | `I8+I9` | `[0,1]` | 2 |

All five are empty at the full marked-basis gate in the displayed exhaustive
normalized charts.  The three pole-`[0,0]` profiles were tested in both twist
classes at 5 and 7.  `S0050` comes closest: its nonsquare charts have four
component-matched pairs at 5 and eight at 7, but two of the latter meet a
singular fibre and every remaining pair has the wrong smooth intersection.
`S0093` and `S0360` never put both generator classes on one model.  At prime
5, `S0071` has 40/16 individual sections and 32 pairs in the nonsquare chart,
but every pair misses the lattice-required intersection four.  The new
`I8+I9` census has 168 squarefree models; `S0197` has 8/0 sections in the
square twist and 0/8 in the nonsquare twist, so the two classes never coexist.

The lattice-only adapter, full source census, `I8+I9` fibre census, its two
marking scans, and the compact aggregate certificate have SHA-256
`7ddd99f37ab31c07ba4f60b8298382ad63c05939d4b31f2d77d9a691ebb9f68a`,
`c3010b3256f0d291dcb1a5dbd285e8489eb6ab3600e71a18d9bc7158e40e4d33`,
`b32a4adbd94ceaac47f1bfba582d7521b0d89c03b8546c41101bce0f410a3721`,
`7adb8496a70f9d4109418cabb66c6c063c0447f8baeda22cdecb4af2eacd12fb`,
`6fdb0a308c01328e67e92871acaa5adc99f851ee08749bdbc788c50435d86306`,
and
`c7d76508a27805f4075dc742855b07d86d81f95b0b4dc722ca9cbf58ce9aca99`.

The rank-first same-NS beam from pole-zero source `S0050` reaches root ranks
13, 9, 4, and 3 in its first four degree-two edges, then stalls at root rank
3 through depth eight.  It does not hit the exact `A1/MW16` target.  The
artifact has SHA-256
`36e29920315709081852f0bc49652025065810d5590877c5be27c8964d977cb7`.
This is a capped, beam-pruned miss, not a graph obstruction.

Determinant 654 is therefore a strong demonstration that abstract low-MW
source abundance is not enough: the rational-marking gate must be scored
before corridor cost.  It is demoted behind the formally smooth determinant-
500 MW1 branch.  Higher-pole profiles, other normalizations, the unresolved
T-arithmetic curve, wider or mixed-degree corridors, and the rootful target's
`D.F=2,3,4` spectrum remain open.  The existing rootless translation-coset
formula is not applied to this `A1` target.

## Determinant 714: a formally smooth MW1 source opposite MW16

Continuing the source-first expansion produces a genuine new promotion on
`K3-cf7f6c91a3a40d32`.  Its named target `F001` has root type `A1` and MW
rank 16.  Both partner auxiliaries were searched in the complete six-large-
ambient semistable MW0--2 cut.  Partner one gives 463 MW2 rows and no MW1
row; 186 have complete bases through pole two.  Partner two gives 548 rows:
546 MW2 and **two MW1**, with 216 and 2 complete low-pole bases respectively.

The two MW1 rows `S0223` and `S0430` are integrally isometric.  Their common
source profile is

```text
reducible fibres:       I5 + I7 + I7
MW group:               Z, height 102/35
torsion:                trivial
basis pole profile:     [1]
component depths:       (1,2,1)
target:                 A1 / MW16
```

The normalized fibre census exhausts all `5^8` and `7^8` degree-eight `A`
polynomials.  It finds 6 squarefree models over `GF(5)` and 20 over `GF(7)`.
Both `GF(5)` twist classes and the square `GF(7)` class are empty at the MW1
marking gate.  The nonsquare `GF(7)` chart is positive: four signed pole-one
sections occur on two distinct fibre models, all with the exact depths
`(1,2,1)`.  The fibre artifacts have SHA-256
`6e35bc342ce730d914c1f98c852222a3ea3fc5def1585235cd47f896b20d1d4d`
and
`61f442f5a3d18aafd844b282ac74e83a2641d8531ebb30569297dba4bff6a448`;
the positive marking artifact has SHA-256
`60550e97b4ab3f18b0869b479970efdbd372acad29bb07d6b29b8817c4483cf7`.

Both geometrically distinct marked seeds have Jacobian rank 39 in 40
variables and lift through `7^8`; their unit minors are 2 and 6 modulo 7.
The apparent 47-equation overdetermination is again exact.  The universal
node/discriminant identity and the `I5,I7,I7` orders with depths `(1,2,1)`
force the degree-at-most-18 section residual to be
`t^2*(t-1)^4*q(t)` with `deg(q)<=10`.  The unit minor retains residual
coefficients 2 through 12, which kill all eleven coefficients of `q`
triangularly.  Hence the marked branch is one-dimensional and formally smooth
over `Z_7`.  The two Hensel and formal certificates have SHA-256
`e38e8b98f47b10b322530201fe319fc516e6f6cb36aa26d1def0aff02f0eae14`,
`d7b4e460fb607f37a4395c0abf6e05695faa02e172e06668303bbb95b7466ad5`,
and
`6aea3aee698ab94b4b3042ee2ef546919ab62ad3f96a8315419c83d303da3a53`.

A bounded rational scan fixes the unique free coordinate `m8` at every
allowed integer in `[-40,40]` across both residue disks.  All 23 candidates
lift through `7^40`, but none reconstructs all forty coefficients over `Q`.
Its SHA-256 is
`2e5f32db3decc1774236463908feb23570c4574b213b8c8cba7f96b0af9009a6`.
This is not an irrationality result.

The marked equations nevertheless give an exact affine curve chart over
`QQ`: 40 coefficients, 39 independent equations on the certified
localization, and relative dimension one at the `GF(7)` seed.  The candidate
T-arithmetic datum has quaternion discriminant 51 and local level 7.  An
independent Ogg-formula calculation gives genus 21 for `X_0^51(7)` and genus
two for its full Atkin--Lehner quotient, agreeing with the published
low-genus table.  This strongly explains the failed rational scan, but the
Eichler-order identification and a birational map from the marked chart to a
genus-two equation remain open.  The algebraization artifact has SHA-256
`b54e40a1dda8b5a2fee91145e5229bffa4004a8a09a526bcd24ddb9102b39309`.

The route compiler formerly compared only root-rank-zero children with named
targets, even in its declared `A2/MW15` and `A1/MW16` cases.  This was a
terminal-recognition bug, not a neighbour-enumeration or finite-field bug.
The stored old beams did not conceal hits: determinants 384, 654, and 714
bottomed out at root ranks 4, 3, and 2, while their targets have root ranks 2,
1, and 1.  The corrected recognizer passes direct regression checks at target
root ranks 0, 1, and 2.

The regenerated determinant-714 beam retains every marked state.  At depth
six its eight states have root ranks `(2,2,3,3,3,4,4,4)`; the two MW15 states
have the same compiler cost `(1,4,2,6,0)` and differ in the final q4 orbit.
The frontier artifact has SHA-256
`43b61a24864f4c0c1f26d8cf2984848fe2cc360a09f1374e9de67cc6834bbc72`.
A targeted last-edge search from both MW15 states tests 65,963 capped
candidates: degree two at q4, q6, q8; degree three at q6, q9, q12; and degree
four at q8, q12, with pole at most two.  The q4 shell has 92 fully physical
survivors but no target; every other declared shell is empty at the physical
gate.  The four artifacts have SHA-256
`15273eb9400742c6f50700dd02872030c745c43726647512832aecde7f3aba9e`,
`943588934fcf6e3e9bfd587889ed659d0b13b740f65ac2a7f07951db08065d77`,
`30c3b7fdf341d2b989cff420c2056abb1a51f1b5c9ecfc617d69de20822ff586`,
and
`ed947f2652d27e3209d3cc4a6ba64b6b3ffa87e49bd85d5b15e2d411ca149ecb`.

The rootful A1 target spectrum uses `(M/dM)/W(A1)`, not the rootless formula
verbatim.  Its complete degree-two low-height census has 98,304 Weyl/section-
translation orbits: 26,908 rational bisections of minimum norm ten, 47,943
genus-one bisections of minimum norm eight, and 27 orbits whose minimum is
above ten.  In deterministic 256-orbit exact-CVP samples, degree three has 40
rational and 68 genus-one candidates; degree four has 22, 39, and 42
candidates of genera zero, one, and two.  These are lattice spectra, not
irreducibility, arithmetic-descent, or rank-jump theorems.  The artifact has
SHA-256
`99f96f29b8ad8b887b4ef92b3859e6fb0075d94378292954c42aea074ffbcb10`.

The aggregate promotion certificate has SHA-256
`a3fc378d7c2faf73ef26d4eed5e9dc451693b9f03535e1cd0312814d4eb90cc5`.

Determinant 714 is therefore the strongest new rank-tradeoff candidate: a
formally smooth MW1 source opposite MW16, with an exact MW15 corridor
intermediate.  It remains behind determinant 500 on target rank and current
equation simplicity, and it still lacks a rational point or explicit
hyperelliptic model for the marked moduli curve, primitive-closure and Picard
audits, an exact target route, and a resolved T-arithmetic identification.

## Determinant 750: rational moduli, but no ideal source in the large-ambient cut

The next genus-zero rootless MW17 target in the expanded Pareto order is
`K3-10a14a46c14b3150-F001`, of determinant 750.  Its arithmetic ledger
identifies the moduli curve as `Gamma_0(3)`, so it is attractive before source
search.  The surface has two catalogue auxiliary classes.  Both have now been
searched completely in the same six large Niemeier ambients, restricted to
semistable MW0--2 sources with one to three reducible supports.  Neither
auxiliary produces an admissible source Gram.  The first has 12 terminal
embeddings in each of `A17+E7` and `A15+D9`, all nonprimitive; the second has
36 in `A17+E7` and none elsewhere, again with no retained source.

This is an exact rejection only for the declared ideal cut, not a complete
fibration classification: smaller-component ambients, additive fibres, MW3,
and four-support sources remain open.  It nevertheless demotes determinant
750 behind determinant 500 for the present equation-first objective.  The
two search artifacts and compact rejection certificate have SHA-256
`8e822ec7d557b85f693560316e3ed61f7228df7e5b1755b5da11fd31e3a31e87`,
`d94d0e78c9295a93a1f319b1af91c24e094d15cb5cf312ee5c8cf86f63389bae`,
and
`70b5fef7775640fe5757e0fd24e715efc2c0adbaabfdda349686869df6f2be2b`.

The determinant-864 genus-zero `Gamma_0(2)` rootless MW17 surface
`K3-3425921cd7db891f` has the same outcome.  Its two auxiliary searches find
many terminal embeddings—312 for the first auxiliary and 332 for the
second—but no semistable MW0--2 source in the declared cut.  Their SHA-256
hashes are
`af9067cf73a6ecf700b4c853177119467732fd20124b537fbc92c2ebb7745551`
and
`5143e0425bb1c7a89253b7e2c7dd010a00148e53a089e45a8f5aa4f94a05021c`.

The determinant-1296 genus-zero `Gamma_0(3)` surface
`K3-49b947f9626a0481` is now rejected in the same scoped cut.  Each of its two
auxiliary classes has exactly 402 terminal embeddings across the selected
ambients—168/144 in `A24`, 48/72 in `A17+E7`, and 186/186 in `A15+D9`—but
every one is nonprimitive and no source Gram survives.  The two exact search
artifacts have SHA-256
`de144b84c019cb1c8a962bac68c1a2891d14dca2646b105bc867cf5664d88570`
and
`af4e04ebfa08ba2fa0bbe5900c7f79cd9d38252f3c385d46f6ce0cf5a6ad5735`.

The last two rational-moduli rootless MW17 surfaces also miss the ideal cut.
For determinant 1500, `K3-99a0b9b18de6e19b`, the sole auxiliary has 120
terminal embeddings, all in `A15+D9` and all nonprimitive.  For determinant
1728, `K3-dc0e324e4ac40dbc`, the sole auxiliary has 360 terminals in `A24`,
152 in `A17+E7`, and 256 in `A15+D9`; again every terminal is nonprimitive.
The search artifacts have SHA-256
`75155f554e943833308e1233637d061d512c15120816d66812c758dcdd798b0c`
and
`8722e6abeba3a4ba150b3ae3e75cf664993e4b01612dc5614756cd5828a86d36`.

A live rational-moduli optimizer now records the resulting queue.  Among the
six catalogue surfaces with rational moduli and a rootless MW17 frame,
determinant 500 remains the sole active promotion and the other five are
scoped ideal-source rejections.  The rootless-MW17 rational queue is therefore
exhausted at this cut; the next source-first expansion is to rational-moduli
MW15/MW16 targets, not to equation work on the rejected MW17 rows.  This
ledger has SHA-256
`8747e1f950585ad6ae1c6e2aa03b2a26442a9c47829325b4fd22576dcbe93b77`.

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

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_golay_det720_3a5_saturation_rejection.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_golay_det720_3a5_rational_parameters.sage --check

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

The determinant-500 source search and promotion are reproduced by

```bash
python3 elkies-k3/scripts/extract_rank7_catalogue_source_search_target.py \
  --surface-id K3-04b86146cc6b284b --partner-index 2 \
  --frame-id K3-04b86146cc6b284b-F001 \
  --output artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-source-search-target-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_golay_det720_prescribed_root_sources.sage \
  --target artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-source-search-target-v1.json \
  --output artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-prescribed-root-sources-large-a-v1.json \
  --source-root-rank-min 14 --source-root-rank-max 17 \
  --source-support-min 1 --source-support-max 3 \
  --ambient-label A24 --ambient-label A17_E7 --ambient-label A15_D9 \
  --ambient-label D24 --ambient-label D16_E8 --ambient-label 3E8 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/sample_lattice_foundry_multisection_spectrum.sage \
  --database artifacts/generated-results/elkies-k3-rank7-auxiliary-catalogue-v1.json \
  --frame-id K3-04b86146cc6b284b-F001 --sample-count 1024 \
  --height-slack 4 \
  --output artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-multisection-spectrum-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/complete_lattice_foundry_degree3_spectrum.py \
  --database artifacts/generated-results/elkies-k3-rank7-auxiliary-catalogue-v1.json \
  --frame-id K3-04b86146cc6b284b-F001 --workers 8 \
  --chunk-size 1000000 --float-type dd --audit-precision 256 \
  --audit-stride 4096 \
  --checkpoint artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-degree3-complete-v1.checkpoint.json \
  --output artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-degree3-complete-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_k3_04b_equation_first_promotion.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0028_source_ansatz_modp.sage \
  --source artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-prescribed-root-sources-large-a-v1.json \
  --ns-id K3-04b86146cc6b284b \
  --source-id K3-04b86146cc6b284b-S0160 --prime 5 --examples 20 \
  --output artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-fibre-ansatz-mod5-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_k3_04b_a3_a4_a9_pole1_marking_modp.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_k3_04b_a3_a4_a9_marked_gf5_hensel.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0028_source_ansatz_modp.sage \
  --source artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-prescribed-root-sources-large-a-v1.json \
  --ns-id K3-04b86146cc6b284b \
  --source-id K3-04b86146cc6b284b-S0160 --prime 7 --examples 0 \
  --output artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-fibre-ansatz-mod7-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_k3_04b_a3_a4_a9_pole1_marking_modp.sage \
  --fibres artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-fibre-ansatz-mod7-v1.json \
  --output artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-pole1-marking-mod7-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_k3_04b_a3_a4_a9_marked_gf5_hensel.sage \
  --fibres artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-fibre-ansatz-mod7-v1.json \
  --marking artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-pole1-marking-mod7-v1.json \
  --output artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-a3-a4-a9-marked-gf7-hensel-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_k3_04b_a3_a4_a9_formal_smoothness.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_k3_04b_a3_a4_a9_source_qq_rejection.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_k3_04b_a3_a4_a9_rational_parameters.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/classify_k3_04b_semistable_mw2_sources.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0028_source_ansatz_modp.sage \
  --source artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-prescribed-root-sources-large-a-v1.json \
  --ns-id K3-04b86146cc6b284b \
  --source-id K3-04b86146cc6b284b-S2021 --prime 5 --examples 0 \
  --output artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-a3-a4-a8-mw2-fibre-ansatz-mod5-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0028_source_ansatz_modp.sage \
  --source artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-prescribed-root-sources-large-a-v1.json \
  --ns-id K3-04b86146cc6b284b \
  --source-id K3-04b86146cc6b284b-S2021 --prime 7 --examples 0 \
  --output artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-a3-a4-a8-mw2-fibre-ansatz-mod7-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_k3_04b_a3_a4_a8_mw2_marking_modp.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_k3_04b_a3_a4_a8_mw2_marking_modp.sage \
  --fibres artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-a3-a4-a8-mw2-fibre-ansatz-mod7-v1.json \
  --quadratic-twist 3 \
  --output artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-a3-a4-a8-mw2-marking-mod7-nonsquare-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_lattice_foundry_same_ns_compiler_routes.sage \
  --case k304b --q 4 --degree 2 --max-pole 1 --beam-width 8 \
  --max-depth 12 --mw-vector-cap 2000 --cap-from-mw-rank 4 --rank-first \
  --output artifacts/generated-results/elkies-k3-k3-04b86146cc6b284b-same-ns-compiler-routes-rankfirst-cap2000-v1.json \
  --check

python3 elkies-k3/scripts/certify_k3_10a_semistable_source_rejection.py --check

python3 elkies-k3/scripts/extract_rank7_catalogue_source_search_target.py \
  --surface-id K3-49b947f9626a0481 --partner-index 1 \
  --frame-id K3-49b947f9626a0481-F001 \
  --output artifacts/generated-results/elkies-k3-k3-49b947f9626a0481-source-search-target-partner1-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_golay_det720_prescribed_root_sources.sage \
  --target artifacts/generated-results/elkies-k3-k3-49b947f9626a0481-source-search-target-partner1-v1.json \
  --output artifacts/generated-results/elkies-k3-k3-49b947f9626a0481-semistable-mw0-2-sources-large-a-partner1-v1.json \
  --source-root-rank-min 15 --source-root-rank-max 17 \
  --source-support-min 1 --source-support-max 3 --all-a-only \
  --ambient-label D24 --ambient-label D16_E8 --ambient-label 3E8 \
  --ambient-label A24 --ambient-label A17_E7 --ambient-label A15_D9 --check

# Repeat the preceding two commands with partner1 replaced by partner2.

python3 elkies-k3/scripts/extract_rank7_catalogue_source_search_target.py \
  --surface-id K3-6ce16abb9de3c7c5 --partner-index 1 \
  --frame-id K3-6ce16abb9de3c7c5-F001 --lattice-only \
  --output artifacts/generated-results/elkies-k3-k3-6ce16abb9de3c7c5-source-search-target-partner1-lattice-only-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/classify_k3_6ce_a5_a10_mw2_sources.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_k3_6ce_a5_a10_fibre_ansatz_modp.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_k3_6ce_a5_a10_mw2_marking_modp.sage \
  --output artifacts/generated-results/elkies-k3-k3-6ce16abb9de3c7c5-a5-a10-mw2-marking-mod5-square-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_lattice_foundry_same_ns_compiler_routes.sage \
  --case k36ce --q 4 --degree 2 --max-pole 1 --beam-width 8 \
  --max-depth 8 --mw-vector-cap 2000 --cap-from-mw-rank 4 --rank-first \
  --output artifacts/generated-results/elkies-k3-k3-6ce16abb9de3c7c5-same-ns-compiler-routes-rankfirst-cap2000-v1.json \
  --check

python3 elkies-k3/scripts/certify_k3_6ce_equation_first_candidate.py --check

python3 elkies-k3/scripts/extract_rank7_catalogue_source_search_target.py \
  --surface-id K3-14ad03cd7c1848b2 --partner-index 1 \
  --frame-id K3-14ad03cd7c1848b2-F001 --lattice-only \
  --output artifacts/generated-results/elkies-k3-k3-14ad03cd7c1848b2-source-search-target-partner1-lattice-only-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_k3_6ce_a5_a10_fibre_ansatz_modp.sage \
  --source artifacts/generated-results/elkies-k3-k3-14ad03cd7c1848b2-semistable-mw0-2-sources-large-a-partner1-v1.json \
  --source-id K3-14ad03cd7c1848b2-S0197 \
  --surface-id K3-14ad03cd7c1848b2 --determinant 654 \
  --schema elkies-k3.k3-14ad-a7-a8-mw2-fibre-ansatz-modp.v1 \
  --prime 5 \
  --output artifacts/generated-results/elkies-k3-k3-14ad03cd7c1848b2-a7-a8-mw2-fibre-ansatz-mod5-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_lattice_foundry_same_ns_compiler_routes.sage \
  --case k314ad --q 4 --degree 2 --max-pole 1 --beam-width 8 \
  --max-depth 8 --mw-vector-cap 2000 --cap-from-mw-rank 4 --rank-first \
  --output artifacts/generated-results/elkies-k3-k3-14ad03cd7c1848b2-same-ns-compiler-routes-rankfirst-cap2000-v1.json \
  --check

python3 elkies-k3/scripts/certify_k3_14ad_equation_first_candidate.py --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_lattice_foundry_ns0028_source_ansatz_modp.sage \
  --source artifacts/generated-results/elkies-k3-k3-cf7f6c91a3a40d32-semistable-mw0-2-sources-large-a-partner2-v1.json \
  --ns-id K3-cf7f6c91a3a40d32 \
  --source-id K3-cf7f6c91a3a40d32-S0223 --prime 7 --examples 0 \
  --output artifacts/generated-results/elkies-k3-k3-cf7f6c91a3a40d32-a4-2a6-mw1-fibre-ansatz-mod7-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_k3_04b_a3_a4_a9_pole1_marking_modp.sage \
  --fibres artifacts/generated-results/elkies-k3-k3-cf7f6c91a3a40d32-a4-2a6-mw1-fibre-ansatz-mod7-v1.json \
  --sources artifacts/generated-results/elkies-k3-k3-cf7f6c91a3a40d32-semistable-mw0-2-sources-large-a-partner2-v1.json \
  --source-id K3-cf7f6c91a3a40d32-S0223 \
  --surface-id K3-cf7f6c91a3a40d32 \
  --schema elkies-k3.k3-cf7f-a4-2a6-mw1-pole1-marking-modp.v1 \
  --quadratic-twist 3 \
  --output artifacts/generated-results/elkies-k3-k3-cf7f6c91a3a40d32-a4-2a6-mw1-pole1-marking-mod7-nonsquare-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_k3_cf7f_a4_2a6_mw1_formal_smoothness.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/scan_k3_cf7f_a4_2a6_mw1_rational_parameters.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_k3_cf7f_moduli_algebraization.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/search_lattice_foundry_same_ns_compiler_routes.sage \
  --case k3cf7f --q 4 --degree 2 --max-pole 1 --beam-width 8 \
  --max-depth 6 --mw-vector-cap 2000 --cap-from-mw-rank 4 --rank-first \
  --retain-frontier-witnesses \
  --output artifacts/generated-results/elkies-k3-k3-cf7f6c91a3a40d32-same-ns-compiler-routes-rankfirst-frontiers-depth6-v2.json \
  --check

# Resume both root-rank-two states at depth six.  The four checked artifacts
# use respectively: (degree,q,cap) = (2,4,2000), (2,6+8,5000),
# (3,6+9+12,5000), and (4,8+12,3000).

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/sample_lattice_foundry_multisection_spectrum.sage \
  --target-artifact artifacts/generated-results/elkies-k3-k3-cf7f6c91a3a40d32-source-search-target-partner2-lattice-only-v1.json \
  --sample-count 256 --height-slack 4 --pari-stack-gb 4 \
  --output artifacts/generated-results/elkies-k3-k3-cf7f6c91a3a40d32-rootful-a1-multisection-spectrum-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/certify_k3_cf7f_equation_first_candidate.sage --check

python3 elkies-k3/scripts/build_rank7_rational_moduli_source_optimizer.py --check
```

All displayed lattice-source data, section-pole profiles, finite-field marking
counts, formal-local statements, rational equation identities, Picard bound,
and degree-three counts are copied from hash-pinned exact artifacts.  The
report constructs no replacement rational point in the intended determinant-720
NS class, rational marking or parameterization of the determinant-500 source,
neighbour route, effective target multisection, or specialization rank jump.
