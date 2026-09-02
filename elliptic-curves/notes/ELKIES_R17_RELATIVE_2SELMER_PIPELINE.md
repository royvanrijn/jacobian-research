# Relative 2-Selmer pipeline for the compact R17 family

## Current outcome

The basis-level inputs cover the held-out rank-21 mechanism control, the
rank-25--28 controls, and the first ten candidates in the frozen height-10000
weakest-block Nagao ranking. The original Magma supervisor records all fifteen
jobs as `backend_unavailable`. An open-source Sage/PARI replacement is now
implemented and validated end to end on a small curve; the large R17 controls
remain resource-bounded computations. No incomplete run supplies a Selmer
dimension, unrealized Selmer class, exact-rank statement, or point-search
success.

A separate class-group-free Kummer audit now completes on all fifteen inputs.
It proves, by exact residue squareclasses, that the five known control
subgroups have mod-2 dimensions `21,25,26,27,28`. Thus the public exceptional
points give respectively `4,8,9,10,11` independent directions modulo the
specialized generic MW17 subgroup, realizing `16,256,512,1024,2048` classes
including zero inside those *known* quotients. This labels the known lower
bound; it does not identify the known quotient with the full Selmer quotient.
The same audit proves mod-2 rank 17 for MW17 on every high-Nagao candidate.

The first pinned open run on the rank-21 control used a 300-second wall limit,
a 4 GB RSS limit, a 2 GB PARI stack, and all twelve proved discriminant-prime
hints. It timed out inside `ellrankinit` at 440,283,136 bytes peak observed
RSS, before BNF certification or `ell2cover`. This is a measured backend
bottleneck, not a Selmer result.

The same frozen method was applied without public-point hints to the top
high-Nagao candidate `t=-5643/6760`. Its 120-second diagnostic run also timed
out inside `ellrankinit`, at 230,608,896 bytes peak observed RSS. No candidate
cover search was reached.

The exact input manifest is
[`elkies_2026_relative_2selmer_suite_inputs_v1.json`](../../artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_suite_inputs_v1.json).
The host/backend audit is
[`elkies_2026_relative_2selmer_suite_run_v1.json`](../../artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_suite_run_v1.json).

| case | `t` | certified known rank | held-out quotient directions |
| --- | ---: | ---: | ---: |
| rank-21 mechanism control | `3/8` | 21 | 4 |
| rank-25 control | `-2/377` | 25 | 8 |
| rank-26 control | `-308/251` | 26 | 9 |
| rank-27 control | `2456/135` | 27 | 10 |
| rank-28 control | `-9529/5471` | 28 | 11 |

These are lower-bound controls.  The displayed ranks are not asserted to be
exact until a matching upper bound exists.

## Construction-derived quotient and stop decision

The certified R17 construction supplies a genuine first quotient before any
global cubic-field descent.  Let `G_t` be the specialization of the generic
MW17 subgroup, let

```text
S_t = Sel_2(E_t/QQ),
H_t = image(G_t/(G_t intersect 2*E_t(QQ))) in S_t,
```

and let `B_t` be the subspace of `S_t` spanned by `H_t` and the Kummer classes
of every rational point obtained by specializing a split member of the
complete 39,120-bisection atlas.  Then `B_t/H_t` is a certified subspace of
the rational-point part of `S_t/H_t`.  Writing `K_t` for the full Kummer image
of `E_t(QQ)/2*E_t(QQ)`, quotienting by it gives the exact sequence

```text
0 -> K_t/B_t
  -> S_t/B_t
  -> Sha(E_t/QQ)[2]
  -> 0.
```

This is a useful *relative* decomposition, but it is not a computation of the
middle term.  In particular, the last column below is residual only inside
the certified known rational-point subgroup, not inside the full Selmer
group:

| control | known quotient over MW17 | split-bisection span | known-point residual |
| --- | ---: | ---: | ---: |
| rank 21 | 4 | 4 | 0 |
| rank 25 | 8 | 5 | 3 |
| rank 26 | 9 | 3 | 6 |
| rank 27 | 10 | 2 | 8 |
| rank 28 | 11 | 1 | 10 |

The construction found the points spanning the middle column without being
given the public exceptional coordinates.  Those coordinates were loaded
only afterwards to label the recovered classes.  Thus the rank-21 control is
a complete blind recovery of its four *known* exceptional directions, while
the thinning `5,3,2,1` spans on ranks 25--28 prove that the same mechanism is
not a complete exceptional-direction detector.

There is a second exact decomposition after adjoining independent bisection
fields.  For `L=QQ(t)(sqrt(q_1),...,sqrt(q_k))`, Theorem F4 of
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](../../elkies-k3/RANK_MUTATION_AND_LIFT_THEOREMS.md)
gives

```text
E(L) tensor QQ = direct_sum_I E^(product_(i in I) q_i)(QQ(t)) tensor QQ.
```

The missing product characters `E^(q_i*q_j)` are the only genuinely new
rank contributions in a paired bisection construction.  This decomposition
uses projectors with denominators `2^k`; it is therefore a rational
Mordell--Weil decomposition, not a direct-sum decomposition of the
2-Selmer group.  It can construct new points, but cannot replace the residual
2-descent.

The construction-native evidence is now sufficient for a stop decision:

1. translated or inverted bisections cannot enlarge the specialization
   quotient, by Corollary F2.1 of the theorem note;
2. the complete split census and exact control classes are recorded in
   [`ELKIES_BISECTION_SPECIALIZATION_CONTROLS.md`](ELKIES_BISECTION_SPECIALIZATION_CONTROLS.md);
3. the product-twist census found no stable extra-character signal, while
   proving no twist-rank upper bound, as recorded in
   [`QUADRATIC_TWIST_RANK_CENSUS_2026-08-31.md`](../../elkies-k3/QUADRATIC_TWIST_RANK_CENSUS_2026-08-31.md);
4. tested degree-three/four equations and exact reverse-q12 transport do not
   explain the ten-class rank-28 packet, as recorded in
   [`R17_EXCEPTIONAL_SPECIALIZATION_RELATIONS_2026-09-02.md`](../../elkies-k3/R17_EXCEPTIONAL_SPECIALIZATION_RELATIONS_2026-09-02.md); and
5. repeated PARI and Hecke relation engines stop at the same uncertified
   cubic `S`-class/`S`-unit gate.

Therefore no further generic backend retry or translated-bisection shell is
authorized by the present evidence.  The problem should be reopened only for
one of two materially new inputs:

- a certified non-torsion section on a predeclared product twist
  `E^(q_i*q_j)`, which would instantiate a new character and can be tested
  blindly on the locked controls; or
- a complete residual class-field computation, equivalently a complete
  enumeration of the locally admissible quartic extensions with the given
  cubic resolvent, after quotienting by the explicit MW17 and bisection
  Kummer images.

Until then, the frozen high-Nagao application stops after the exact MW17 and
complete bisection stages.  A zero bisection gain remains “mechanism not
seen,” never a Selmer or rank upper bound.

### Exact known-subgroup Kummer stage

`audit_elkies_2026_known_kummer_quotients.py` avoids both maximal orders and
class groups. Completing the square gives the two-division cubic with root
`theta`; after setting `zeta=4*theta`, it evaluates the equivalent class
`4*x(P)-zeta` in every squarefree residue factor at selected odd primes. A
global square has square residue at every such place, so full row rank is an
exact certificate of global squareclass independence. The primes are selected
only to certify the supplied points and do not form an injective model for
unknown Selmer classes.

The complete five-control run took 11.03 seconds with 241,100 KiB peak RSS;
the complete fifteen-case run took 11.50 seconds with 241,004 KiB peak RSS.
The largest selected primes for ranks 21, 25, 26, 27, and 28 were respectively
`163,223,281,271,283`. This is the successfully frozen open-source
embedding/labeling stage. It cannot construct unrealized quotient covers
because the full quotient basis remains unavailable.

The known subgroup nevertheless supplies a complete explicit-cover control
corpus. `build_elkies_2026_known_quotient_covers.py` enumerates every nonzero
combination of the exceptional basis on each control, adds the corresponding
rational points exactly, and constructs the two quadrics for
`alpha=X(P)-theta`. All 3,851 covers have the verified rational point
`[1:0:0:1]`: respectively `15,255,511,1023,2047` classes on ranks
21, 25, 26, 27, and 28. The compact generated manifest retains all 42 basis
covers and hashes the 25 MiB all-class ledger under `artifacts/local/`.
These are realized positive controls, not unknown Selmer classes or blind
recoveries.

## Frozen computation

### Open-source path

`run_elkies_2026_relative_2selmer_open.py` uses PARI's `ell2cover`, whose
output is a basis of everywhere locally soluble binary-quartic 2-covers with
maps to the elliptic curve. For each selected case it:

1. starts an isolated Sage worker with the minimal model but no generic or
   exceptional points;
2. calls `ellrankinit` and requires `bnfcertify` to return one on every
   cubic-field BNF in the rank context;
3. calls `ell2cover`, records every basis quartic and its map, and runs a
   bounded `hyperellratpoints` search without public point hints;
4. reloads the exact generic and held-out points only after the worker exits;
5. uses exact good-prime quotients to express recovered cover images in the
   known Mordell--Weil basis; and
6. emits point-to-Selmer rows when the blindly recovered basis classes span
   the known image.

A full `ell2cover` return plus successful BNF certification is an actual
2-Selmer basis. A bounded quartic-search miss does not prove that its class is
unrealized. PARI supplies the basis quartics but does not expose addition of
arbitrary cover classes, so explicit covers for non-basis quotient
combinations still use the repository's cubic-etale intersection-of-quadrics
layer.

The newer
`run_elkies_2026_relative_2selmer_checkpointed.py` makes that last layer
explicit and avoids repeating `ellrankinit`. It first replaces the curve cubic
by the instantaneous `polredbest` field model, stores exact maps in both
directions, and persists only a `bnfcertify`-accepted BNF. A modified Simon
descent then computes the field squareclass group and all local conditions in
the reduced model; its real sign and finite local-image maps evaluate the
original curve root through the stored isomorphism. Each resulting class is
transported back to the original cubic before its two quadratic equations are
written.

MW17 is embedded in a separate exact worker. The following blind worker
enumerates all nonzero quotient classes up to a declared cap, constructs both
quadrics for each representative, and searches a declared coefficient box.
Only after it exits are the public exceptional coordinates mapped and the
known exceptional subgroup labels attached. Thus a search hit can measure
recovery of an exceptional direction without having supplied that point.
Every search miss remains bounded evidence only.

This decomposition is integration-tested on small curves: the transported
Selmer dimension agrees with `ell2cover`, the known point maps to the expected
squareclass, and the constructed cover recovers a rational point. On the
rank-21 control, however, the first original-field checkpoint attempt still
failed in BNF relation collection: `--bnf-tech 0.1,4,20` reached 600 seconds
at 256,798,720 bytes observed RSS. It emitted no Selmer dimension. The two
follow-up global engines fail the same gate: a reduced-field, locally tuned
PARI 2.17.4 run times out at 34,402,304 bytes peak RSS, and Hecke 0.40.2
`method=2` with `GRH=false` times out at 3,954,507,776 bytes peak RSS. Neither
timeout is an upper bound. Consequently the frozen protocol stops at the
rank-21 control before any *full-Selmer* work on rank 25--28 or high-Nagao
inputs. The independent known-subgroup audit above is cheap and has been run
on all of them.

### Optional Magma cross-check

For each specialization the generated Magma job performs the following steps.

1. Reconstruct the global minimal fibre and the seventeen specialized generic
   sections from the pinned compact R17 model.
2. Construct multiplication by two and call
   `SelmerGroup([2] : Bound := -1, Raw := true)`, the shared-map form of
   `TwoSelmerGroup`. `Bound=-1` requests unconditional class-group data; no
   GRH class-group bound is installed.
3. Use `DescentMaps([2])` and the returned `AtoS` map to compute the actual
   Selmer coordinate `AtoS(mu(P_i))` of every generic section.  The job aborts
   unless these rows have rank 17.
4. Extend those rows to a basis of the full 2-Selmer group and thereby obtain
   an explicit basis of

   ```text
   Sel_2(E/Q) / image(<P_1,...,P_17>).
   ```

5. Before declaring any public exceptional points, materialize quotient
   representatives with `TwoCover`, record their quartic equations, and run a
   bounded `Points` search.  The generic x-coordinates are the only Selmer
   hints.  Thus the measured recovered quotient span is genuinely blind to the
   public extra-point coordinates.
6. Declare the held-out control points only after the blind phase, map them
   into the same quotient basis, compute their span, and label every stored
   cover class as inside or outside that known rational span.

If a quotient has at most 255 nonzero classes, every nonzero class is built
and searched.  Above that threshold the job builds and searches a canonical
quotient basis.  This preserves a spanning set of unexplained directions
without pretending that an exponential all-class search is feasible.  The
parser explicitly records whether enumeration was exhaustive.

The parser rejects missing stages, source-hash changes, inconsistent
dimensions, incomplete generic/quotient bases, and class-count mismatches.  A
bounded cover-search miss remains a negative experiment, not evidence for a
Tate--Shafarevich class.

## High-Nagao application

The full-Selmer method has been instantiated on these first ten prospective
candidates; its global class-group backend has been executed on the first
candidate only:

```text
-5643/6760, 1452/7817, 4298/8873, -7634/2859, -841/8544,
461/4420, 6695/1353, 1217/151, 9783/7559, 9446/3605.
```

Their exact minimal models, specialized generic sections, Nagao records, and
program hashes are in the input manifest.  They remain heuristic candidates
until the complete descent returns.  In particular, no candidate is promoted
to point search merely because it was selected by Nagao score.

The frozen class-group-free stage has been applied to all ten and certifies
that MW17 specializes with mod-2 rank 17 in every case, using largest
auxiliary primes between 131 and 277. With no extra points supplied, this
does not recover an exceptional direction and supplies no residual Selmer
dimension.

## Backend calibration

Input reconstruction and generation of all fifteen source-pinned jobs took
about 13.1 seconds in the first local run.  This is an input-generation
benchmark, not a descent benchmark.

The official PARI development branch at commit
`6af5b91cfaeb6939331945f301e65bd775f6cdef` adds a six-parameter threaded
relation engine absent from the installed PARI 2.17 backends.  A local
pthread build passes a certified small-cubic smoke test.  On the reduced
rank-21 field, eight threads are the best tested host setting: the 300-second
run reaches factor-base bound `16348` with 1,996 ideals, executes 6,143
random-relation rounds, and still requests 1,635 relations in 1,630 ideals.
It consumes 2,190.47 user CPU seconds and 136,872 KiB peak RSS.  Increasing
the ideal power, using twelve threads, or raising the trial-factor budget is
neutral or worse.  Thus the newer engine improves factor-base traversal but
does not complete the class/unit group or reach certification.  Its internal
relation request is telemetry, not an S-class or Selmer dimension.

Four further PARI technical settings, two Hecke relation methods, three exact
bounded relation ledgers, and an aggressive PARI restart build are pinned in
`elkies_2026_relative_2selmer_open_bottleneck_benchmarks_v2.json`. None
completed the global class/unit group. The first two bounded ledgers collected
449 and 96 large-prime cycles but gained zero relations after quotienting the
declared S-primes; their residual 34-dimensional bounded model is not a
Selmer bound. A third exact-factor, multi-large-prime hypergraph run widens
the factor base to 1,000. It closes 666 hypergraph dependencies in 26.54
seconds, but all resulting rank lies in the canonical S-span: quotient-rank
gain is again zero. Its remaining 125-dimensional bounded model is likewise
not an S-class or Selmer bound.

A further PARI `bnfinit` flag-0 run with `tech=[0.03,4,20]` reached the same
600-second relation-collection stop at 260,067,328 bytes. Flag 0 is now an
explicit runner option and passes a certified small-field smoke test, but it
does not improve this R17 field. A Hecke monitor also checkpoints the exact
bound-240 relation state modulo its 16 visible S columns: 14 principal rows,
augmented rank 23 of 57, and residual dimension 34 after 120 seconds. Missing
large S ideals and factor-base generation make that a bounded model only.

The specialized follow-up removes the first of those caveats. It augments the
factor base by all 25 prime ideals above the declared rational S-primes,
including the two very large primes, without triggering Hecke's prime-counting
heuristic up to their norms. This gives 66 columns, 15 initial principal rows,
augmented mod-2 rank 32, and the same residual dimension 34. Neither the normal
collector nor a new direct collector that enumerates in products of a target
ideal and an S-ideal adds a quotient row in 120 seconds; their peak RSS values
are 1,561,526,272 and 1,574,998,016 bytes. Three seeded archimedean parameter
variants also have zero gain. All S ideals are visible in these last two
models, but factor-base generation, units, and local conditions are still
missing, so dimension 34 remains neither an S-class nor a Selmer bound.

Earlier open-source entry points did not provide the requested basis-level
result on the controls:

- eclib/mwrank failed on the rank-21 minimal model in about one second with
  `lower bound on c too large`, before returning a Selmer rank;
- PARI `ellrank`, with all 21 certified points supplied, did not return within
  a strict 60-second calibration;
- PARI `ellrankinit` with all twelve proved rank-21 factor hints did not return
  within the pinned 300-second open-suite run and peaked at 440,283,136 bytes
  observed RSS;
- the already pinned rank-28 eclib and PARI attempts both timed out after 300
  seconds, and the factor-supplied PARI attempts timed out after 600 and 1800
  seconds.

None of these outcomes is a Selmer upper bound. PARI 2.17.3 additionally
exposes `ell2cover`, which supplies the previously missing full locally
soluble cover basis. The open runner separately certifies the field BNF and
reconstructs the point embedding from blindly recovered cover images. Magma's
raw interface remains an optional independent cross-check.

The mathematical reduction through a cubic etale algebra and an `S`-class/
`S`-unit computation follows Schaefer--Stoll,
[*How to do a p-descent on an elliptic curve*](https://mathe2.uni-bayreuth.de/stoll/papers/p-descent-long.pdf).
The explicit-cover minimization and reduction layer is consistent with
Cremona--Fisher--Stoll,
[*Minimisation and reduction of 2-, 3- and 4-coverings of elliptic curves*](https://arxiv.org/abs/0908.1741).

## Replay

Generate the five controls and first ten candidates:

```bash
python3 elliptic-curves/cas/build_elkies_2026_relative_2selmer_suite.py \
  --output-dir artifacts/local/elliptic-curves/elkies-2026-relative-2selmer-suite-v1 \
  --manifest artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_suite_inputs_v1.json \
  --candidate-count 10 --search-bound 1000 \
  --enumerate-class-limit 255 --overwrite
```

Certify the known Kummer subgroup on all fifteen inputs without a class group:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/audit_elkies_2026_known_kummer_quotients.py \
  --manifest artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_suite_inputs_v1.json \
  --prime-bound 5000 \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_known_kummer_quotients_suite_v1.json \
  --overwrite
```

This command checkpoints after every case. Its `PASS` status certifies only
the supplied Kummer classes and their known relative quotient coordinates.

Build every explicit cover in the known exceptional control subgroups:

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/build_elkies_2026_known_quotient_covers.py \
  --overwrite
```

The generated result is compact; the command writes the complete hashed
all-class ledger to `artifacts/local/elliptic-curves/`.

Run the rank-21 control with the open-source backend and explicit wall/RSS
limits:

```bash
python3 elliptic-curves/cas/run_elkies_2026_relative_2selmer_open.py \
  --manifest artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_suite_inputs_v1.json \
  --output-dir artifacts/local/elliptic-curves/elkies-2026-relative-2selmer-open-v1 \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_open_rank21_300s_v1.json \
  --case control-r21-t3_8 --timeout-per-case 300 \
  --rss-limit-bytes 4000000000 --pari-stack-bytes 2000000000 \
  --search-bound 1000 --certificate-prime-bound 1000 --overwrite
```

Repeat `--case` for controls or candidates, or use `--controls-only`. The
output is self-classifying: only `COMPLETE_CERTIFIED_PARI_TWO_SELMER_BASIS`
contains a full basis and quotient calculation.

For an optional licensed Magma cross-check, supervise every case with explicit
wall/RSS limits:

```bash
python3 elliptic-curves/cas/run_elkies_2026_relative_2selmer_suite.py \
  --manifest artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_suite_inputs_v1.json \
  --log-dir artifacts/local/elliptic-curves/elkies-2026-relative-2selmer-suite-v1/logs \
  --output artifacts/local/elliptic-curves/elkies_2026_relative_2selmer_suite_run.json \
  --timeout-per-case 86400 --rss-limit-bytes 16000000000 --overwrite
```

Only after every job completes, parse the source-matched transcripts:

```bash
python3 elliptic-curves/cas/parse_elkies_2026_relative_2selmer_suite.py \
  --manifest artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_suite_inputs_v1.json \
  --log-dir artifacts/local/elliptic-curves/elkies-2026-relative-2selmer-suite-v1/logs \
  --output artifacts/generated-results/elliptic-curves/elkies_2026_relative_2selmer_suite_results_v1.json
```

The parser cannot turn a timeout, backend failure, or partial log into a
result.  `MATH_STATUS.json` therefore remains unchanged.
