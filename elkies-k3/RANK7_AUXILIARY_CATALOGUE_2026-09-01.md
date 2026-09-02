# Rank-seven auxiliary catalogue — 2026-09-01

## Exact imported catalogue

The first surface-first catalogue layer is now replayable. It imports the
exact `N(2A7+2D5)` one-root foundry shell, the bounded symmetry-first `2C`
fixed-lattice seed shell, the bounded `N(6A4)` double-swap shell, and the
exact determinant-720 `N(24A1)` Golay-octad design, then deduplicates in the
required order:

```text
(transcendental lattice T, Neron--Severi lattice NS)
    -> rank-seven partner auxiliaries K
    -> rank-seventeen frames M=K^perp
    -> ambient primitive embeddings.
```

The current exact imported inventory has

```text
161 (T,NS) surface classes
180 rank-seven partner-auxiliary isometry classes
724 frame isometry classes
1,027 retained primitive embedding-orbit records.
```

All 724 imported frames have generic Mordell--Weil rank in the requested
high-rank window: 20 have MW12, 190 have MW13, one has MW14, 113 have MW15,
261 have MW16, and 139 have MW17. This is an
inventory of exact records, not a completeness statement for any determinant
band.

The machine-readable catalogue is
[`../artifacts/generated-results/elkies-k3-rank7-auxiliary-catalogue-v1.json`](../artifacts/generated-results/elkies-k3-rank7-auxiliary-catalogue-v1.json).
Its content-derived `K3-...` surface identifiers supersede sequential
`NS0001...` labels as cross-backend catalogue keys; every imported legacy ID
is retained as an alias.

## Determinant and backend ledger

The initial determinant range `1 <= |det(T)| <= 5000` is split into four
bands:

```text
D0001-0500
D0501-1000
D1001-2000
D2001-5000.
```

There are 24 separate backends: the 23 rooted Niemeier lattices and one Leech
backend. Their Cartesian product with the four determinant bands gives 96
explicit enumeration shards. All 96 remain open. The four rooted backends
with imported or bounded search records carry narrower statements:

- `ROOTED-2A7_2D5` contains the declared 768-element one-root mutation shell
  around the two determinant-948 controls and the exact 7-of-16 coordinate
  seed shell in a pinned `2C` fixed-lattice basis; neither is a complete
  determinant-band census;
- `ROOTED-24A1` contains one exact Golay-octad design found by a bounded,
  nonexhaustive proposal search;
- `ROOTED-4D6` closes one exact component-transposition coordinate shell but
  finds no MW12--17 seed;
- `ROOTED-6A4` imports one exact double-swap coordinate shell canonicalized
  under its complete 240-element chamber-preserving residual group.

The `24A1` backend now also has an exact canonical-prefix layer and a
contiguous rank-seven completion frontier.  The M24 orbit counts of unordered Golay-octad
subsets of sizes one through five are

```text
1, 3, 16, 206, 10547.
```

At every size the independent orbit-stabilizer mass equals
`binomial(759,k)`.  The size-three count `16` independently agrees with the
published octad-triple classification of Kelsey--Rowley.  Completing prefix
indices `0:2000` in eight contiguous shards by two octads, imposing
coordinate-union size at least 19, and saturating in the integral
Construction-A model gives 1,267 shard-local residual-M24 records with primitive
determinant at most 500.  Their determinant distribution is

```text
384: 5, 448: 1, 480: 809, 486: 8, 500: 444,
```

and their MW-rank distribution is `MW12: 260, MW13: 970, MW14: 37`.  Exact
ternary discriminant-form gates retain 1,253 of the 1,267 local records.  Forty
retained spans have primitive-closure index two; they demonstrate why a raw
determinant cutoff is invalid.

The full Weyl-sign quotient is now exact for this declared input.  Using the
24 doubled physical coordinate covectors, intrinsic auxiliary isometries, and
exact M24 transporters across all eight shards, the 1,267 local records collapse to 23
`2^24 semidirect M24` embedding orbits.  Their determinant distribution is
`384:3, 448:1, 480:13, 486:1, 500:5`, their MW-rank distribution is
`MW12:5, MW13:13, MW14:5`, and 18 pass the ternary-genus gate.  Every member
has an explicit row-isometry, coordinate-permutation, and coordinate-sign
witness; every representative has an exact full stabilizer and orbit size.

These records are not merged into the main surface catalogue yet: only 2,000
of 10,547 five-prefix orbits have been completed, and the generator language
contains only positive octad vectors.  Moreover, the ternary output is a
genus gate rather than a class enumeration of all possible transcendental
lattices `T`.

## Cross-Niemeier mod-2 stabilizer priority

The cross-backend scheduler now implements the proposed umbral experiment.
For `A7^2 D5^2`, the exact `Dih_4` section contains classes
`1A,2A,2B,2C,4A`, while the six existing full ambient stabilizers realize
only `1A,2A`; their complement image is `{+I,-I}` and is trivial on `M/2M`.
Future searches therefore prioritize stabilizers containing `2B`, `2C`, or
`4A` in this backend, and analogous non-scalar component permutations in the
other rooted backends.

Across all 23 rooted Niemeiers the deterministic scheduling tiers contain
`1,9,7,6` backends: the exact `A7^2 D5^2` seed, at-least-four-component
permutation envelopes, swap envelopes, and no-repeated-component controls.
This multiplicity ranking is only a proposal filter.  An embedding enters the
high-priority experiment only after its full ambient stabilizer is computed
and an induced complement action satisfies

```text
rank_GF2(g_M - I) > 0.
```

Only then are fixed-point and orbit distributions on the rational subset of
`M/2M` used to rank source searches.

The first exact symmetry-first family is now closed.  The two `4A` elements
are inverse and have a common primitive rank-eight fixed lattice `F` of
determinant 4,096; in the computed chamber section all eight `Dih_4` elements
fix `F` pointwise.  Every primitive corank-one `K` in `F` is therefore fixed
pointwise by literal `2B`, `2C`, and `4A` elements.  The dual-lattice identity

```text
det(K) = det(F) * a H_F^(-1) a^t
```

makes the determinant-5,000 search finite and complete in this family.  It
gives 336 primitive embeddings, three auxiliary isometry classes, and five
frame classes.  Their MW-rank distribution is
`MW13:16, MW15:16, MW17:304`; every requested outer class acts nontrivially on
the complement modulo two, with `4A` moving nine dimensions.

This entire attractive-looking family nevertheless fails the K3 gate.  Every
auxiliary/frame discriminant group has length seven, while a rank-three
transcendental lattice has discriminant length at most three.  Accordingly,
none of the 42 determinant-2,048 or 50 determinant-4,096 even ternary genera
matches.  This is an exact exclusion of the pointwise-`4A`-fixed corank-one
family, not an exclusion of rank-seven auxiliaries with a nontrivial `4A`
action on `K`; the latter is the next symmetry-first search space.

The rank-sixteen `2C` fixed lattice gives the first positive symmetry-first
catalogue increment.  In a pinned integral LLL basis, all
`binomial(16,7)=11,440` coordinate direct summands are tested.  After the
determinant, discriminant-length, and MW12--17 gates, 97 seeds remain.  Closing
under the exact `Dih_4` section gives 97 embedding orbits of size four; each
has a literal `2C` stabilizer whose action on `M/2M` moves six, seven, or eight
dimensions.  The 73 distinct frame discriminant forms all have exactly one
matching ternary genus.

Required `(T,NS)`-first deduplication produces 73 new surface classes, 76
partner-auxiliary classes, and 86 frame classes.  Their post-frame-dedup MW
distribution is `MW12:16, MW13:66, MW14:1, MW15:3`.  Thirteen surfaces have
two inequivalent frames; several pair MW12 and MW13 presentations of the same
surface.  The computation is exact inside the declared coordinate-summand
language and exact under the residual `Dih_4` section.  It does not enumerate
all primitive rank-seven sublattices of `Fix(2C)` and does not yet quotient by
the full Weyl group.

The first exact cross-backend comparison beyond `A7^2D5^2` separates a useful
positive shell from a negative control. For `N(4D6)`, a lifted `S4` section
contains literal component transpositions with primitive rank-sixteen fixed
lattice of determinant 256. All 11,440 coordinate rank-seven summands of its
pinned LLL basis are tested: 183 fail discriminant length and the remaining
11,257 have MW rank below 12. Thus this declared shell contributes no target,
without excluding other `4D6` auxiliaries.

For `N(6A4)`, exhaustive glue testing shows that the chamber-preserving
component/diagram group has order 240 and component-permutation image order
120; the abstract `S6` envelope is therefore not used as if every permutation
lifted. A cycle-shape `1^2 2^2` involution has primitive rank-sixteen fixed
lattice of determinant 256. Its 11,440 coordinate summands leave 161 MW12--13
seeds after all lattice gates. Every selected involution acts nontrivially on
`M/2M`, moving six or seven dimensions. Exact canonicalization under all 240
residual elements gives 161 orbits of size 120; all 42 distinct discriminant
forms have one matching ternary genus.

Local `(T,NS)`-first deduplication gives 42 surfaces, 55 auxiliary classes,
and 128 frames, with post-dedup distribution `MW12:4, MW13:124`. Global
deduplication adds only 39 surfaces: three already occur in `N(2A7+2D5)` and
gain alternative `6A4` frames. This is direct evidence that surface-first
cross-backend merging is necessary. The shell is not all primitive
rank-seven sublattices of `Fix(g)` or `N(6A4)`, and its exact residual quotient
does not replace the open full Weyl quotient.

The other 19 rooted backends are not yet enumerated. The separate Leech
backend now has an exact ambient Gram matrix and two certified `Co0`
generators, but no primitive rank-seven embedding orbit has yet been
enumerated. In particular, the catalogue does not turn the old mutation shell
into a determinant-5000 theorem by relabeling it.

## Orbit policy

For a rooted Niemeier lattice `N^X`, embeddings must be canonicalized under
the full group generated by the Weyl group and the residual umbral group

```text
W(X) normal in Aut(N^X),
G^X = Aut(N^X)/W(X).
```

For `24A1`, this means Golay-code support combinatorics and the residual
`M24` action. For the Leech lattice it means `Co0` orbits, with no Weyl
layer. The 290 Hoehn--Mason orbits classify fixed-point sublattices of the
Leech lattice, not arbitrary primitive rank-seven sublattices, so they are
useful orbit/stabilizer infrastructure but not a completeness shortcut.

The `(T,NS)` key uses the exact ternary Gram matrix together with the NS
discriminant quadratic module, represented in the source frame's canonical
normal basis with its quadratic Gram negated. In every imported record the
discriminant length is at most three. Nikulin's indefinite uniqueness
criterion therefore applies to signature `(1,18)`: the discriminant form
determines the NS isometry class. Different ternary lattices with the same
discriminant form remain different surface records.

## Typed Pareto discovery ledger

The catalogue-wide discovery ordering is
[`../artifacts/generated-results/elkies-k3-rank7-surface-pareto-v1.json`](../artifacts/generated-results/elkies-k3-rank7-surface-pareto-v1.json).
Its universal frontier minimizes the tuple

```text
(- maximum catalogued generic MW rank,
 easiest known exact source MW rank,
 source reducible-fibre support count,
 determinant).
```

All 161 surfaces have these four exact lattice metrics. Forty-eight use the
existing external low-MW source certificates; the other 113 use only their
easiest currently imported frame, with that weaker inventory scope recorded
in the row. The core Pareto frontier has four surfaces:

```text
K3-ebaf00b3723751ba  (MW17 from MW1, support 1, determinant 950)
K3-8188cdcda8c57b2d  (MW17 from MW1, support 2, determinant 948)
K3-f43753fb154e3406  (MW17 from MW17, support 0, determinant 720)
K3-14ad03cd7c1848b2  (MW16 from MW1, support 1, determinant 654).
```

Missing arithmetic data are never imputed. Separate enriched frontiers use
only rows carrying the required evidence: 39 surfaces have exact minimum-pole
data, 117 have nontrivial exact stabilizer evidence, nine have the existing
degree-two-exact/degree-three-and-four-bounded multisection data, and 34 have
rootless target short-vector data. No surface currently has a certified
physical neighbour route, so that coverage frontier is empty. Equation
construction, field of definition, conductor prospects, and moduli
genus/rationality remain typed `UNKNOWN` unless an input artifact certifies
them. Pareto leadership is a search priority, not an optimality or arithmetic
rank theorem.

## Replay

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_prefix_orbits.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_prefix_orbits.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 0 --prefix-stop 250

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 0 --prefix-stop 250 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 250 --prefix-stop 500

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 250 --prefix-stop 500 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 500 --prefix-stop 750

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 500 --prefix-stop 750 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 750 --prefix-stop 1000

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 750 --prefix-stop 1000 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 1000 --prefix-stop 1250

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 1000 --prefix-stop 1250 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 1250 --prefix-stop 1500

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 1250 --prefix-stop 1500 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 1500 --prefix-stop 1750

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 1500 --prefix-stop 1750 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 1750 --prefix-stop 2000

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_24a1_octad_rank7_completion_shard.sage \
  --prefix-start 1750 --prefix-stop 2000 --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_24a1_octad_completion_manifest.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_24a1_octad_completion_manifest.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_24a1_weyl_m24_shard.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_24a1_weyl_m24_shard.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_cross_niemeier_mod2_priority.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_cross_niemeier_mod2_priority.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_2a7_2d5_4a_fixed_rank7.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_2a7_2d5_4a_fixed_rank7.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_2a7_2d5_2c_fixed_high_mw_seed.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_2a7_2d5_2c_fixed_high_mw_seed.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_4d6_swap_fixed_high_mw_seed.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_4d6_swap_fixed_high_mw_seed.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_6a4_double_swap_fixed_high_mw_seed.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/enumerate_6a4_double_swap_fixed_high_mw_seed.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_leech_co0_backend.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_leech_co0_backend.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_rank7_auxiliary_catalogue.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_rank7_auxiliary_catalogue.sage --check

python3 elkies-k3/scripts/build_rank7_surface_pareto.py

python3 elkies-k3/scripts/build_rank7_surface_pareto.py --check
```

## Next exact frontier

The next useful computation is not another mutation of the determinant-948
H3 auxiliary. It is completion of the first orbit-complete backend/band shard, with the
enumerator required to emit:

1. a proof that every primitive rank-seven embedding in the declared shard
   occurs;
2. canonicalization under `W(X) semidirect G^X` (or `Co0` for Leech);
3. the saturated complement and full root system;
4. the exact ternary discriminant-form gate;
5. a shard-completeness certificate consumed by this catalogue.

The natural first rooted target remains `24A1:D0001-0500`.  Full Weyl-sign
canonicalization is closed for the declared `0:2000` input; the next gates are
sharded completion of prefix indices `2000:10547` and enlargement beyond the
positive-octad generator language.  Neither may be replaced by the present
23-orbit partial count.
The natural Leech target is a minimal-vector-generated rank-seven shard with a
declared reduced-basis bound; it must remain separate because rootlessness is
automatic there and the orbit group is `Co0`.

## Literature

- K. Nishiyama, *The Jacobian fibrations on some K3 surfaces and their
  Mordell--Weil groups*, Japanese Journal of Mathematics **22** (1996),
  293--347, DOI `10.4099/math1924.22.293`.
- V. Kelsey and P. Rowley, *M24-Orbits of Octad Triples*, Graphs and
  Combinatorics **34** (2018), 1429--1443, DOI
  `10.1007/s00373-018-1961-1`.
- V. Nikulin, *Integral symmetric bilinear forms and some of their
  applications*, Mathematics of the USSR-Izvestiya **14** (1980), 103--167,
  DOI `10.1070/IM1980v014n01ABEH001060`.
- M. Cheng, J. Duncan, and J. Harvey, *Umbral Moonshine and the Niemeier
  Lattices*, arXiv:`1307.5793`.
- G. Hoehn and G. Mason, *The 290 fixed-point sublattices of the Leech
  lattice*, Journal of Algebra **448** (2016), 618--637,
  DOI `10.1016/j.jalgebra.2015.08.028`.
