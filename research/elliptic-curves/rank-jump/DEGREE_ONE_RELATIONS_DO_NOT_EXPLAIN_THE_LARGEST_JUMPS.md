# A uniform relation census finds rational pair events on moderate-gain fibres

The first family-wide relation test now has an exact result. Among the
519 co-split native pairs in the fixed published-R17 panel, **131 admit a
generic degree-one translated intersection**, but only **18 such relation
components lie over the tested parameters**. Those incidences occur on six
fibres with retained gains +4 through +8. None occurs on the historic
+9/+10/+11 controls.

This is a positive mechanism for simultaneous pair solubility, with an
explicit limit: each incidence identifies two native quotient directions
with one line. Its prevalence does not explain the largest retained jumps.
The [relation comparison CSV](../../artifacts/generated-results/elliptic-curves/rank_jump_degree_one_relation_discrimination_v1.csv)
contains all 165 published-R17 addresses, joined to the prior rank labels
after the equation-only computation. Full curve ranks remain UNKNOWN.

A separate test checks all **2,853 co-split triples**. None has even trace
sum. Consequently every integer-translate triple relation in this panel
has **total intersection degree at least six**. This does not exclude a
rational degree-one component inside a larger intersection algebra.

## A single generic rule, applied uniformly

Let M be the marked published-R17 generic lattice, with height h and
minimum four. A native bisection with trace w has class

\[
B_w=((h(w)-2)/4,\,2,\,w)\in U\oplus M(-1).
\]

Every native trace in the atlas has height ten. For any unordered pair
of distinct atlas labels i,j, the rule is:

\[
R\in w_i+w_j+2M,\quad h(R)=6,\qquad
S=(R+w_i+w_j)/2.
\]

Then

\[
B_{w_i}\cdot(S-B_{w_j})=h(R)/2-2=1.
\]

The two Q-defined rational curves therefore have a zero-dimensional
intersection of length one, a reduced rational point. This is the
[previous norm-six solubility theorem](NORM_SIX_INTERSECTION_FORCES_NATIVE_PAIR_SOLUBILITY.md),
now used as a uniform rule rather than on a selected successful pair.
The theorem does not locate its parameter at an arbitrarily supplied t.

There is at most one norm-six representative up to sign in a coset. If
R and R' were distinct and not opposite, both R+R' and R−R' would be
nonzero vectors of 2M, hence have height at least 16. Their two heights
sum to 24 by the parallelogram identity, a contradiction. A norm-four
representative likewise excludes a norm-six representative in the same
coset. The two signs of R give conjugate translates S and w_i+w_j−S.
The census counts their paired intersection components once, modulo
simultaneous native conjugation.

Counts refer to labelled relation schemes. Different pairs can share an
elliptic point or a quotient line, so 18 incidences do not mean 18 distinct
points or independent constructions.

This rule is defined for the **whole 39,120-cover family atlas**. At a
smooth nonbranch panel fibre, an incidence requires both covers to split.
The prior complete square census therefore restricts the possible
incidences to its 519 co-split pairs without losing any. This is exact
lazy evaluation of the family-wide rule at the fixed addresses; it does
not select pairs by rank or by an exceptional point relation. There is
no need to enumerate or compute intersection parameters outside the panel.
The underlying hypothesis was motivated by earlier retrospective evidence;
this is not a prospective holdout claim.

The 519 pairs happen to give 519 distinct parity cosets. Each coset is
completely enumerated through norm six. The producer uses PARI after
exact lattice reduction; the independent verifier uses rational LDL
bounds and integer endpoint rounding, visiting 95,532 recursion nodes.
It separately verifies the lattice minimum with 4,180 nodes.

## Rational components must occur at the right parameter

For each eligible pair, both native branch choices and both conjugate
translates are checked by exact elliptic group arithmetic at the already
supplied parameter. The tests receive generic lift equations and generic
sections, not the exceptional witness coordinates. The 118 relevant lift
maps satisfy their elliptic equations coefficient by coefficient. Every
specialized trace identity is checked separately.

| Published parameter | Retained gain | Split covers | Globally forced degree-one pairs among them | Components at this parameter, modulo conjugation | Native quotient upper bound from these relations |
|---|---:|---:|---:|---:|---:|
| 3/8 | +4 | 25 | 69 | 7 | 18 |
| −70/61 | +5 | 13 | 25 | 6 | 8 |
| 8/39 | +5 | 7 | 4 | 2 | 5 |
| 1229/894 | +6 | 2 | 1 | 1 | 1 |
| 28/117 | +6 | 11 | 13 | 1 | 10 |
| −2300/843 | +7 | 3 | 2 | 0 | 3 |
| −288/65 | +7 | 4 | 0 | 0 | 4 |
| −2/377 | +8 | 6 | 6 | 0 | 6 |
| −4112/1937 | +8 | 4 | 2 | 1 | 3 |
| −308/251 | +9 | 3 | 0 | 0 | 3 |
| 2456/135 | +10 | 2 | 0 | 0 | 2 |
| −9529/5471 | +11 | 1 | 0 | 0 | 1 |

The upper bounds concern only the points from the displayed native
covers, not the full exceptional quotient or the whole curve. They need
not be sharp. For example, the +7 quartet at −288/65 has known native
quotient rank three, but this pair dictionary gives only the upper bound
four; its relevant dependence is the previously certified triple relation.

The difference between the two middle columns is essential. At −2/377,
all six eligible pair carriers have a rational point forced by geometry,
but none of those distinguished points lies above −2/377. Global carrier
solubility alone does not answer the specialization question.

For the existing R17 score-matched +7/observed-zero case, −2300/843 versus
−1929/3242, the relation counts are 0 versus 0 despite split-cover counts
3 versus 0. The height-matched control 1576/2331 gives the same outcome.
The observed-zero labels are censored, not exact rank-17 certificates.
All tested observed-zero fibres have at most one split cover, so their
lack of pair incidences is already implied by the cover census; it is not
additional evidence of discrimination.

## What the signed graph proves

Choose the positive rational square root in each native cover map and
write its point as P_i. The other branch is w_i(t)−P_i. A branch-sum
relation therefore becomes

\[
\epsilon_i P_i+\epsilon_j P_j\in G_t,
\qquad \epsilon_i,\epsilon_j\in\{1,-1\},
\]

where G_t is the specialized marked generic subgroup. Each row gives a
linear relation on the native quotient over Q. The independent verifier
computes the exact rational rank of these rows; the resulting bound is
n minus that rank. It does not replace this rank with a raw edge count.
At −70/61, six edges have rank five, so the bound is eight rather than
seven.

The producer's `generic_word` is the translate for the selected *branches*.
The joined report also records the canonical signed relation word: subtract
w_i whenever the negative-root branch of cover i was used. Thus the stored
translate is not incorrectly presented as the right-hand side for signed
positive-root points.

A connected component leaves at most one quotient line, and an inconsistent
signed cycle would kill that line over Q. This is why the degree-one
construction cannot itself explain several independent directions merely
by making many labels soluble. The data also show that it explains only
part of the dependence among larger native cover sets.

## The lowest-degree triple gate

For any three native traces, a zero residual trace would require

\[
w_i+w_j+w_k\in2M,\qquad S=(w_i+w_j+w_k)/2.
\]

Signs do not change this parity condition. It is tested on every co-split
triple at every panel address, with no rank-based selection. All 2,853
fail it. This is a structural gate, not an enrichment statistic; no
random-triple baseline is inferred.

The [triple intersection class calculation](MINIMAL_CARRIER_AND_RATIONAL_SPLITTING_OF_A_TWO_DIRECTION_BLOCK.md)
applies uniformly here. The generic pair-sum map has trivial stabilizer
because its two nonzero anti-invariant characters are different, and hence
is birational to its degree-four image. The checked branch divisors are
disjoint from each other and from all 24 nodal fibres. For the image of a
translated signed pair, write

\[
D_z=(h(z)/2+2,4,2z),\qquad B_{w_i}=(2,2,w_i).
\]

Then

\[
B_{w_i}\cdot D_z=h(z-w_i)+2.
\]

The residual z−w_i lies in the triple trace-sum coset modulo 2M, regardless
of the branch signs and integral generic translate. Since that coset is
nonzero in every tested triple, the lattice minimum gives

\[
\boxed{\deg Z_{i,j,k,S}\ge 4+2=6.}
\]

The [barrier certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_triple_degree_barrier_v1.json)
checks squarefreeness, disjointness from the singular fibres, all 519 pair
branch gcds, and all triple parities. The inherited generic-height argument
then supplies the class identity above.

This bound concerns the **total proper intersection length**. A degree-six,
eight, ten or twelve algebra may still contain rational degree-one factors.
The bound does not exclude rational triple blocks; it excludes automatic
rationality from a total length-one intersection and the proposed
zero-residual length-two shortcut on these triples.

## Revised test priorities

- **Solubility, verified but dimension-limited:** degree-one pair intersections
  force rational lifts. Their incidences concentrate here on moderate-gain
  fibres and identify quotient lines; they do not explain large jumps.
- **Solubility, still the strongest open block candidate:** splitting of a
  higher-degree triple relation scheme into rational components, possibly
  with several components or fourth lifts at the same t. A uniform rule for
  such schemes can be fixed using short residual trace cosets.
- **Weak explanation:** counting globally rational pair carriers without
  locating their distinguished points over the tested fibre. The six misses
  at −2/377 give a direct counterexample to that inference.

The next bounded panel test should use residual norm four or six, giving
triple relation schemes of total degree six or eight. Apply the same
short-coset rule across all co-split triples, test its branch relations at
the frozen parameters, and compare the resulting rational incidences and
quotient relations. This targets the first remaining degrees instead of
further elaborating an individual quartet. Full higher-degree component
counts and A1/MW16 coverage remain open; no rank selector is changed.

## Reproduction

The protocols are [DEGREE_ONE_RELATION_PANEL_PROTOCOL.json](DEGREE_ONE_RELATION_PANEL_PROTOCOL.json)
and [TRACE_ZERO_TRIPLE_PANEL_PROTOCOL.json](TRACE_ZERO_TRIPLE_PANEL_PROTOCOL.json).
The initial generic input, all short-coset certificates, per-fibre relation
records, independent verification, and joined CSV are immutable. Completed
coset and fibre workers are separately checkpointed under the ignored local
rank-jump directory. All workers completed within their declared bounds.

```sh
sage -python elliptic-curves/rank-jump/verify_degree_one_relation_panel.py check
sage -python elliptic-curves/rank-jump/triple_degree_barrier.py check
python3 elliptic-curves/rank-jump/report_degree_one_relation_panel.py check
```

The verifier uses neither PARI short-vector enumeration nor Sage elliptic
group arithmetic for its independent checks. PARI is used only to obtain a
unimodular reduction basis for the minimum check; exact LDL enumeration
proves the bound. Generic trace-class identification is inherited from the
pinned atlas, while generic lift identities and all specialized branch
relations are replayed here. No active search file or rank status is changed.
