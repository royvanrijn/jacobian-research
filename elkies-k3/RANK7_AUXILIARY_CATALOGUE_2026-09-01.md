# Rank-seven auxiliary catalogue — 2026-09-01

## Exact imported catalogue

The first surface-first catalogue layer is now replayable. It imports the
exact `N(2A7+2D5)` one-root foundry shell and the exact determinant-720
`N(24A1)` Golay-octad design, then deduplicates in the required order:

```text
(transcendental lattice T, Neron--Severi lattice NS)
    -> rank-seven partner auxiliaries K
    -> rank-seventeen frames M=K^perp
    -> ambient primitive embeddings.
```

The current exact imported inventory has

```text
49  (T,NS) surface classes
49  rank-seven partner-auxiliary isometry classes
510 frame isometry classes
769 retained primitive embedding records.
```

All 510 imported frames have generic Mordell--Weil rank in the requested
high-rank window: 110 have MW15, 261 have MW16, and 139 have MW17. This is an
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
explicit enumeration shards. All 96 remain open. The two imported rooted
backends carry narrower statements:

- `ROOTED-2A7_2D5` is complete only in the declared 768-element one-root
  mutation shell around the two determinant-948 controls;
- `ROOTED-24A1` contains one exact Golay-octad design found by a bounded,
  nonexhaustive proposal search.

The `24A1` backend now also has an exact canonical-prefix layer and a first
rank-seven completion shard.  The M24 orbit counts of unordered Golay-octad
subsets of sizes one through five are

```text
1, 3, 16, 206, 10547.
```

At every size the independent orbit-stabilizer mass equals
`binomial(759,k)`.  The size-three count `16` independently agrees with the
published octad-triple classification of Kelsey--Rowley.  Completing prefix
indices `0:500` in two contiguous shards by two octads, imposing
coordinate-union size at least 19, and saturating in the integral
Construction-A model gives 291 shard-local residual-M24 records with primitive
determinant at most 500.  Their determinant distribution is

```text
384: 1, 448: 1, 480: 221, 486: 4, 500: 64,
```

and their MW-rank distribution is `MW12: 71, MW13: 195, MW14: 25`.  Exact
ternary discriminant-form gates retain 285 of the 291 local records.  Eight
retained spans have primitive-closure index two; they demonstrate why a raw
determinant cutoff is invalid.

The full Weyl-sign quotient is now exact for this declared input.  Using the
24 doubled physical coordinate covectors, intrinsic auxiliary isometries, and
exact M24 transporters across both shards, the 291 local records collapse to 16
`2^24 semidirect M24` embedding orbits.  Their determinant distribution is
`384:1, 448:1, 480:10, 486:1, 500:3`, their MW-rank distribution is
`MW12:4, MW13:10, MW14:2`, and 13 pass the ternary-genus gate.  Every member
has an explicit row-isometry, coordinate-permutation, and coordinate-sign
witness; every representative has an exact full stabilizer and orbit size.

These records are not merged into the main surface catalogue yet: only 500
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

The other 21 rooted backends are not yet enumerated. The separate Leech
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
  elkies-k3/scripts/canonicalize_24a1_weyl_m24_shard.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/canonicalize_24a1_weyl_m24_shard.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_cross_niemeier_mod2_priority.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_cross_niemeier_mod2_priority.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_leech_co0_backend.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_leech_co0_backend.sage --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_rank7_auxiliary_catalogue.sage

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_rank7_auxiliary_catalogue.sage --check
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
canonicalization is closed for the declared `0:500` input; the next gates are
sharded completion of prefix indices `500:10547` and enlargement beyond the
positive-octad generator language.  Neither may be replaced by the present
16-orbit partial count.
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
