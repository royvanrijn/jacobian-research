# Source-first optimization for the Picard-19 lattice foundry

## Outcome

The foundry objective is now ordered as

```text
same NS/T class with a catalogued MW15--17 frame
-> MW0--2 source fibration
-> rational source marking and low equation complexity
-> cheap certified elliptic-neighbour corridor
-> low-degree multisection richness on the target.
```

Rootlessness is no longer a source-search gate.  A rootless MW17 fibration is
one possible endpoint, not the object from which equation work must start.
The search starts at any exact MW15, MW16, or MW17 frame, fixes its underlying
Neron--Severi class, and searches other primitive fibrations of that same K3.

The first bounded implementation inspected 20 exact source candidates attached
to catalogued high-rank frames.  It found no MW0--2 source.  The best existing
candidate remains the NS0024 semistable-compatible `A3+A4+A6/MW4` source.  The
best source newly reached by starting directly at a high-rank frame is the
NS0005 semistable-compatible `A1+2A3+A6/MW4` source.  NS0005 has 40 catalogued
MW15--17 target frames, so it is a useful broad endpoint class, but the source
still misses the preferred MW band and has no rational marking or equation
certificate.

These are bounded discovery results.  They do not prove that any of the tested
Neron--Severi classes lacks an MW0--2 fibration.

## Source score and proof boundary

The exact ranking artifact orders candidates lexicographically by:

1. the preferred band `MW<=2`, then MW rank;
2. number of reducible-fibre supports;
3. compatibility with a semistable all-`A` configuration;
4. expected fibre-stratum dimension;
5. minimum nonzero-section pole order;
6. known rational marking before unknown marking, then Galois orbit size;
7. expected number of additional coefficient conditions;
8. certified neighbour cost, with an unknown route ranked last;
9. the five audited low-degree multisection coordinates as a final tie-break.

Root rank, MW rank, support count, all-`A` compatibility, and the displayed
minimum pole order are exact lattice computations.  The deformation count

```text
expected fibre-stratum dimension = 18 - root rank = 1 + MW rank
```

and the resulting estimate of `MW rank` additional section conditions needed
to isolate a Picard-19 locus are heuristics until an equation ansatz is
constructed.  A one-dimensional complex lattice-polarized moduli space does
not imply a rational parameter over `QQ`.  Rational source marking, Galois
orbit size, and rational parametrization remain explicitly unknown unless an
arithmetic certificate supplies them.

The final multisection tie-break maximizes, in order, rational bisections,
genus-one bisection candidates, sampled rational trisection candidates,
sampled genus-one trisection candidates, and sampled low-genus quadrisection
candidates.  Only the degree-two entries are complete.  This last coordinate
cannot outrank source feasibility or a certified corridor.

The current leading rows are:

| NS | source root type | MW | supports | all-`A` | high-rank endpoints | route |
|---|---:|---:|---:|---:|---:|---|
| NS0024 | `A3+A4+A6` | 4 | 3 | yes | 5 | certified 13-edge degree-two route to MW17 |
| NS0005 | `A1+2A3+A6` | 4 | 4 | yes | 40 | unknown |
| NS0022 | `A1+A2+A3+A6` | 5 | 4 | yes | 13 | unknown |
| NS0005 | `A1+2A3+A5` | 5 | 4 | yes | 40 | unknown |
| NS0033 | `2A2+A3+D5` | 5 | 4 | no | 40 | unknown |

The ranking is reproduced by

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/score_lattice_foundry_sources.sage --check
```

from
[`../artifacts/generated-results/elkies-k3-lattice-foundry-source-ranking-v2.json`](../artifacts/generated-results/elkies-k3-lattice-foundry-source-ranking-v2.json).

## High-rank-frame search

The source hunter now accepts any exact foundry frame, rather than requiring a
rootless start.  The first direct trials used one catalogued MW15 or MW16 frame
in each of eight Neron--Severi classes, twelve generations, beam width twelve,
60 sampled admissible Kneser neighbours per parent, and 7,981 reduced keys per
run.  The target root rank was fifteen, equivalently source MW at most two.

| starting frame | starting MW | best exact source | source MW |
|---|---:|---|---:|
| NS0002-F003 | 15 | `2A1+2A2+2A3` | 5 |
| NS0005-F001 | 15 | `A1+2A3+A6` | 4 |
| NS0011-F003 | 16 | `2A1+A2+A3+D5` | 5 |
| NS0022-F003 | 15 | `3A2+2A3` | 5 |
| NS0024-F003 | 15 | `2A1+2A3+A4` | 5 |
| NS0028-F001 | 16 | `A1+A2+A3+A5` | 6 |
| NS0032-F001 | 16 | `2A1+A3+A6` | 6 |
| NS0033-F001 | 15 | `2A1+2A2+D5` | 6 |

Every retained row is an exact primitive root-lattice/MW computation.  The
Kneser walk is discovery provenance, not an elliptic-neighbour corridor, so
its edge count is not used as equation cost.  The negative result is complete
only for the declared deterministic beams and samples.

For example, replay the strongest new high-rank-start row with

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/hunt_lattice_foundry_rootful_source.sage \
  --ns-id NS0005 --target-frame-id NS0005-F001 \
  --generations 12 --beam 12 --samples-per-parent 60 \
  --primes 3,7,11,13,17,23 --seed 20262906 \
  --target-root-rank 15 --allow-below-target \
  --output artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-mw2-source-from-high-mw-scout-v1.json \
  --root-adapted-frame-output artifacts/generated-results/elkies-k3-lattice-foundry-ns0005-mw2-source-from-high-mw-scout-root-adapted.txt
```

The other seven JSON artifacts use the same
`elkies-k3-lattice-foundry-nsNNNN-mw2-source-from-high-mw-scout-v1.json`
naming pattern and record their seed, admissible prime list, generation
accounting, and visited-key count.

The next algorithmic improvement should search primitive auxiliary embeddings
with the root target built into the Niemeier enumeration, especially all-`A`
root types of rank 15--17 with two or three supports.  That is more aligned
with the objective than merely lengthening the random-neighbour beams.  Only
after a candidate passes rational marking and source-equation gates should a
physical neighbour corridor be optimized.

## Low-degree multisection spectrum

Proposition F5 in
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md)
reduces rootless degree-`d`, genus-`g` section-nonnegative classes to coset
minima in `M/dM`, with threshold

```text
2*d^2 - 2*g + 2.
```

The degree-two calculation is complete through norm ten.  It exactly
reproduces the published R17 count of 39,120 geometrically rational bisection
translation orbits and finds several foundry endpoints with more; the largest
in this nine-frame batch is NS0032-F011 with 41,421, about 5.9 percent above
R17.  NS0028-F005 has 41,376 and NS0033-F026 has 40,912.  This confirms that
R17 is not extremal even for the exact minimal rational-bisection coordinate.

Genus-one bisection counts are exact lattice-candidate counts, but global
nefness, irreducibility, and arithmetic descent are not yet certified.
Degree-three and degree-four entries use 256 deterministic residue classes per
target and exact closest-vector minima within that sample.  They are not
complete censuses and must not be promoted to geometric curve counts.

The artifact is reproduced by

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/sample_lattice_foundry_multisection_spectrum.sage \
  --sample-count 256 --height-slack 4 \
  --frame-id NS0001-F001 --frame-id NS0002-F007 \
  --frame-id NS0005-F008 --frame-id NS0011-F002 \
  --frame-id NS0022-F011 --frame-id NS0024-F005 \
  --frame-id NS0028-F005 --frame-id NS0032-F011 \
  --frame-id NS0033-F026 --check
```

and stored at
[`../artifacts/generated-results/elkies-k3-lattice-foundry-multisection-spectrum-v1.json`](../artifacts/generated-results/elkies-k3-lattice-foundry-multisection-spectrum-v1.json).

Multisection richness is a secondary discovery coordinate, not a specialization
rank theorem.  The R17 positive controls already show that 39,120 bisections
can leave an extreme specialization largely invisible.  The geometric
motivation for retaining this coordinate is the relation between multisections
and rank jumps studied by Garbagnati--Salgado, while the use of alternative
elliptic fibrations to obtain rank jumps is consistent with Salgado's earlier
two-fibration method:

- A. Garbagnati and C. Salgado,
  [*Rank jumps and Multisections of elliptic fibrations on K3 surfaces*](https://arxiv.org/abs/2505.15159).
- C. Salgado,
  [*On the rank of the fibers of rational elliptic surfaces*](https://arxiv.org/abs/1307.3994).

The same-K3 fibration search is grounded in the Kneser--Nishiyama framework;
the bounded foundry catalogue is not a replacement for a complete fibration
classification:

- K. Nishiyama,
  [*The Jacobian fibrations on some K3 surfaces and their Mordell--Weil groups*](https://doi.org/10.4099/math1924.22.293).
- I. Shimada,
  [*On elliptic K3 surfaces*](https://arxiv.org/abs/math/0505140).
