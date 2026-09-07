# Simultaneous solubility has obstructions invisible to every triple

Four explicit bisections of the original published R17 family have **no
simultaneous rational lift at any parameter**, because their equations have
no simultaneous solution over \(\mathbf Q_{23}\). Nevertheless, every triple
of those four equations has a nonsingular solution modulo 23 and hence a
\(\mathbf Q_{23}\)-point. In the broader frozen experiment every triple is
soluble at all fourteen tested primes.

This is a statement about **simultaneous solubility**, not point-search
visibility or an upper bound on the rank of any elliptic fibre. It gives a
concrete reason that candidate generic directions cannot always be combined
independently into a rational block. It does not identify the global event
making a successful large block rational.

## Experiment and result

The [protocol](LOCAL_SOLUBILITY_BLOCK_PROTOCOL.json) uses the 14 generic
bisections which split in the previous
[completed-cohort experiment](SOLUBILITY_FIRST_ON_COMPLETED_FIBRES.md).
Selection of these labels is retrospective. Their equations require no
exceptional point coordinates. The grouping is four, four, three, two, and
one covers on five distinct fibres, with no repeated label.

All 91 pairs and all 16,384 subsets of this fixed dictionary are examined.
The primes are 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, and 43.
The computation enumerates \(\mathbf P^1(\mathbf Z/p^2)\) at odd primes and
\(\mathbf P^1(\mathbf Z/32)\) at 2. A missing simultaneous square residue is
an exact local obstruction. Positive certificates give one rational base
coordinate at which all the required values are squares in \(\mathbf Q_p\),
checking their exact valuations and unit squareclasses. Passing the necessary
residue test without such a witness remains UNKNOWN.

This enumerates 8,580 **local residue classes**, not new rational candidates
for Agent 1's search. All maps, inputs, and outputs are separate from that
search. No product-twist section computation was run: the earlier pair-cover
work proves the two singleton directions but supplies no extra product
character. This cheaper compatibility experiment tests a prerequisite for
combining those directions first.

| Cover combinations | Total | Soluble at every tested prime | Locally obstructed | Unknown |
|---|---:|---:|---:|---:|
| Pairs observed to split together | 16 | 16 | 0 | 0 |
| Pairs from different observed groups | 75 | 75 | 0 | 0 |
| All triples | 364 | 364 | 0 | 0 |
| All quartets | 1,001 | 985 | 14 | 2 |

The test also finds 16 minimal excluded sets of size five, three of size
six, and one of size seven, in addition to the fourteen quartets. Here
“minimal” means that every immediate proper subset passes the necessary
finite tests; it does not upgrade an UNKNOWN proper subset to a local point.
All sets of thirteen or fourteen covers are excluded by at least one tested
prime. These are bounds on simultaneous lifts within this **selected cover
dictionary**, not bounds on independent rational directions on a curve.

The pairwise-gate hypothesis fails in this dictionary: it distinguishes
none of the 16 observed pairs from the 75 cross-group pairs. All triples
pass as well. Four-way and larger compatibility is a real extra condition,
but most quartets still pass. Only two quartets were observed to split
rationally in the 32-parameter cohort. The other passing quartets are not
proved globally insoluble, and the finite-place tests do not constitute an
everywhere-local certificate: omitted primes and the real place are not
covered.

## An explicit four-way obstruction

Remove the exactly verified rational square contents from these atlas cover
equations, without changing rational or local square values:

| Name | Atlas label | Primitive quadratic f(t) |
|---|---|---|
| A | orbit-030cb | 4865126421024 + 2514185838528 t + 320914613929 t² |
| B | orbit-03da0 | 32749355881225 + 8983755404520 t + 183939299856 t² |
| C | orbit-07086 | -163392175124 + 42862125644 t + 167098581169 t² |
| D | orbit-11278 | -144492039 - 201200094 t + 18383017 t² |

Modulo 23, they reduce to

\[
f_A=16+2t+6t^2,\quad f_B=3+2t+7t^2,\quad
f_C=21+4t,\quad f_D=10+19t+14t^2.
\]

Treat them as homogeneous quadratics at infinity. The common square-value
support of A, B, and C on \(\mathbf P^1(\mathbf F_{23})\), including zero
as a square, is exactly \(\{1,14,21\}\). D is nonsquare at each of these:
its values are respectively 20, 7, and 5. Therefore all four equations
cannot hold simultaneously even modulo 23. Projective enumeration includes
infinity, so it leaves no omitted denominator chart.

Yet deleting any one equation leaves a unit solution:

| Omitted equation | t mod 23 | Square roots of the other three values, in A/B/C/D order with the omitted entry removed |
|---|---:|---|
| A | 13 | 4, 2, 1 |
| B | 18 | 8, 1, 9 |
| C | 17 | 6, 6, 3 |
| D | 1 | 1, 9, 5 |

All listed roots are nonzero modulo 23. Holding the integer parameter fixed,
Hensel's lemma lifts each square root. Thus every triple has an actual
\(\mathbf Q_{23}\)-point, while the quartet has none. Each of the four covers
also has a rational lift individually, inherited from the completed cohort;
this does **not** assert a rational lift of every triple.

The generic bisection construction still supplies independent singleton
characters after the corresponding multiquadratic base change. This example
locates a failure of the next implication precisely:

\[
\text{generic candidate directions, individually soluble}
\ \not\Longrightarrow\quad
\text{simultaneously locally soluble at a common parameter}.
\]

It therefore blocks a proposed rational specialization construction even
though generic character independence and all pairwise local checks pass.

## What this changes next

1. **Solubility: retain joint local conditions on whole proposed blocks.**
   A list of pairwise-compatible covers discards information. The forbidden
   combinations in this experiment are exact construction obstructions
   available from cover equations before exceptional points are supplied.
   They cannot veto a curve's rank or other point constructions.
2. **Solubility: local compatibility alone is a weak explanation of the
   observed positive grouping.** Its complete failure to separate pairs,
   and the 985 surviving quartets, direct the next computation toward global
   rational points on the auxiliary carriers or their torsor obstructions.
   First close all relevant local places for a chosen survivor before calling
   a remaining obstruction global or Sha. Then seek a norm/descent
   trivialization that can explain multiple surviving lifts together.
3. **Incidence: extra sections on a product twist remain a separate open
   gate.** This experiment neither constructs nor excludes one. Generic
   independence of singleton lifts does not settle this gate.
4. **Visibility: none of these counts measures point-search efficiency.**
   The selected successful labels and observed grouping are retrospective;
   no prospective ranking advantage has been validated.

The strongest positive evidence remains the earlier explicit construction
of subblocks of dimensions at least 3, 3, 3, and 2. The new result explains
why merely checking each cover, pair, or even triple for local compatibility
cannot justify assembling a larger block. The missing positive theorem is
still a shared arithmetic condition that makes the surviving system
**globally** rationally soluble and preserves enough quotient independence
to account for an extreme jump.

## Reproduction

```sh
python3 elliptic-curves/rank-jump/local_solubility_blocks.py check
python3 elliptic-curves/rank-jump/verify_local_solubility_blocks.py check
```

The independent verifier re-enumerates all residue classes, checks the
maximal permitted sets and every retained exact local witness, and recomputes
the subset accounting without the producer's closure algorithm. It also
records the small modulo-23 proof and all four triple witnesses explicitly.

- [Frozen inputs](../../artifacts/generated-results/elliptic-curves/rank_jump_local_solubility_block_inputs_v1.json)
- [Full local compatibility certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_local_solubility_blocks_v1.json)
- [Independent verification and modulo-23 certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_local_solubility_block_verification_v1.json)
