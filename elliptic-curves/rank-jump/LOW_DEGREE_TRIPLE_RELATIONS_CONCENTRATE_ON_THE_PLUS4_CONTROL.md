# Degree-six/eight triple relations concentrate on the +4 control

The next uniform discrimination test is negative for the largest tested
jumps. Among **2,853 co-split triples at 165 frozen published-R17 addresses**,
575 admit an integral generic translate with total intersection degree six
or eight. Exactly 11 labelled rational incidences occur at the tested
parameters, all at **t=3/8, retained gain +4**. None occurs on a tested fibre
with retained gain +5 through +11.

The [complete comparison CSV](../../artifacts/generated-results/elliptic-curves/rank_jump_low_degree_triple_discrimination_v1.csv)
includes the censored zero controls. The
[independent certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_low_degree_triple_panel_verification_v1.json)
replays every short-vector enumeration and every eligible branch equation.
The [protocol](LOW_DEGREE_TRIPLE_PANEL_PROTOCOL.json) was written before the
calculation. No exceptional point input, new parameter, point search, or
active-search change was used.

This extends the [degree-one pair test](DEGREE_ONE_RELATIONS_DO_NOT_EXPLAIN_THE_LARGEST_JUMPS.md).
It does not claim that all low-degree rational factors on all carriers have
been enumerated. In particular, rational factors of a degree-twelve scheme
are outside a total-degree-at-most-eight rule.

## The condition tested is defined before exceptional points

Let M be the marked integral generic rank-17 Mordell--Weil lattice. Its
minimum is four. A native bisection has two branch maps P_i and its conjugate,
whose sum is the generic trace tau_i of height ten. For a triple put

\[
w=\tau_i+\tau_j+\tau_k,\qquad R\in w+2M,\qquad
S_R=(w+R)/2\in M.
\]

The uniform condition at t is the existence of rational roots satisfying

\[
u_a^2=q_a(t)\quad(a=i,j,k),\qquad
P_i(t,u_i)+P_j(t,u_j)+P_k(t,u_k)=S_R(t).
\tag{1}
\]

Every coefficient and map in (1) comes from the pinned generic atlas.
Thus (1) is an exact, pre-point simultaneous-solubility condition for this
relation carrier. It is stronger than requiring the three quadratics to
split. It is not a closed classification of all possible rational t.
Here it is evaluated only at already frozen addresses.

The [verified intersection calculation](../../artifacts/generated-results/elliptic-curves/rank_jump_triple_degree_barrier_v1.json)
gives total proper intersection length

\[
\ell(Z_{ijk,R})=h(R)+2.
\tag{2}
\]

This uses disjoint native branch pairs over smooth fibres and the nodal
generic model. Every tested trace sum is nonzero in M/2M. Consequently
h(R)>=4; the first two possible lengths are six and eight. Enumerating
all representatives of norm at most six is exhaustive for these lengths,
including both conjugate translates. It does not impose a bound on the
coordinates of R in the original basis.

There is at most one unoriented representative of norm at most six in a
nonzero coset. Otherwise both nonzero vectors R+R' and R-R' lie in 2M,
so each has norm at least sixteen, contradicting the parallelogram
identity h(R+R')+h(R-R')=2h(R)+2h(R')<=24.

The uncut native three-cover carrier has degree eight over t and genus
five. Equation (1) cuts out a finite scheme on it; the six/eight in (2)
are **total intersection lengths**, not genera or guaranteed irreducible
factor degrees. A recorded rational incidence supplies a rational closed
point. Its scheme multiplicity has not been computed.

The full-atlas rule can be evaluated lazily on the co-split triples:
a rational solution of (1) at a smooth nonbranch parameter requires all
three q_a(t) to be rational squares. The previously completed 39,120-cover
census therefore identifies all triples that could contribute at these
addresses. The computation uses no gain-based choice of triples.

## Counts and same-family comparisons

There are 2,827 distinct nonzero parity cosets among the 2,853 triples.
The independent rational LDL enumerations visit 524,488 nodes in total.
Twenty-seven triples have a norm-four representative and 548 have a
norm-six representative. Of the corresponding labelled relation schemes,
four and seven, respectively, have incidences at the tested parameters.
Counts are modulo simultaneous conjugation of all three branches and the
translate. They are not counts of independent elliptic points.

All gains below are inherited retained subgroup gains relative to the
same rank-17 generic subgroup. Full ranks remain UNKNOWN.

| Published t | Retained gain | Split native covers | Co-split triples | Eligible degree 6 / 8 | Incidences degree 6 / 8 | Extra relation rank beyond pairs |
|---|---:|---:|---:|---:|---:|---:|
| 3/8 | +4 | 25 | 2,300 | 26 / 438 | 4 / 7 | 6 |
| -70/61 | +5 | 13 | 286 | 1 / 65 | 0 / 0 | 0 |
| 8/39 | +5 | 7 | 35 | 0 / 7 | 0 / 0 | 0 |
| -3115/2756 | +5 | 3 | 1 | 0 / 1 | 0 / 0 | 0 |
| 28/117 | +6 | 11 | 165 | 0 / 30 | 0 / 0 | 0 |
| 44/35 | +6 | 3 | 1 | 0 / 0 | 0 / 0 | 0 |
| -2300/843 | +7 | 3 | 1 | 0 / 0 | 0 / 0 | 0 |
| -288/65 | +7 | 4 | 4 | 0 / 0 | 0 / 0 | 0 |
| 33/119 | +7 | 7 | 35 | 0 / 5 | 0 / 0 | 0 |
| -2/377 | +8 | 6 | 20 | 0 / 2 | 0 / 0 | 0 |
| -4112/1937 | +8 | 4 | 4 | 0 / 0 | 0 / 0 | 0 |
| -308/251 | +9 | 3 | 1 | 0 / 0 | 0 / 0 | 0 |
| 2456/135 | +10 | 2 | 0 | 0 / 0 | 0 / 0 | 0 |
| -9529/5471 | +11 | 1 | 0 | 0 / 0 | 0 / 0 | 0 |

Three comparisons clarify what the null result means:

1. The existing score-matched +7 / observed-zero pair, -2300/843 versus
   -1929/3242, gives no eligible relation on either fibre. The height-matched
   observed-zero alternative 1576/2331 also gives none; its parameter height
   differs from 2300 by about 1.35%. Native cover splitting separates these
   retained outcomes, but this low-degree triple rule adds no separation.
   Observed zero is censored, not a proof of rank seventeen. Matching and
   search-exposure limitations are inherited from the broader panel.
2. The +6 / +7 pair 28/117 and 33/119 has nearly identical parameter height,
   117 versus 119. There are thirty versus five eligible degree-eight
   translates, yet neither fibre realizes one at its parameter. Generic
   intersection capacity does not imply incidence at the prescribed t.
3. The +7 / +9 pair -288/65 and -308/251 has parameter heights 288 and 308.
   Neither has a degree-six/eight relation. The former's previously proved
   degree-twelve rational triple relation is a positive regression outside
   this dictionary, preventing the null result from being mistaken for an
   absence of every low-degree rational component.

The contrast with 3/8 is descriptive, not a matched statistical experiment:
its parameter is much smaller and it supplies 2,300 of the 2,853 triples.
The observations are neither independent trials nor a random fibre sample.
In particular, the +10/+11 rows have too few native covers even to form a
triple; their zeros primarily expose a mismatch between this native
dictionary and their retained quotient directions.

## The positive result measures dependence

For the canonical positive-root point Q_i, the other branch is tau_i-Q_i.
Each branch relation therefore becomes

\[
\epsilon_i\overline Q_i+\epsilon_j\overline Q_j+
\epsilon_k\overline Q_k=0
\quad\text{in }(E_t(\mathbb Q)\otimes\mathbb Q)/M_t.
\]

The certificate stores both the original branch-sum translate and the
corrected generic word for this signed equation. Negative branches require
subtracting their traces from the original translate.

At 3/8 the eleven triple rows have rank eleven. The seven pair rows have
rank seven. Together their rank is thirteen, so their spans overlap in
dimension five and the triples contribute six new constraints. The
resulting native quotient upper bound improves from eighteen to twelve:

\[
\dim_{\mathbb Q}\langle\overline Q_1,\ldots,\overline Q_{25}\rangle
\le 25-13=12.
\]

This is not a full-curve upper bound, nor an exact native rank. The inherited
retained gain +4 does not by itself prove that all 25 native points lie in
those four retained directions. Nevertheless, the exact relations show
why counting soluble covers or labelled relation events can overcount
directions. This panel supplies no positive evidence that the degree-six/
eight mechanism creates the largest independent jump blocks.

Nor would counting fourth lifts repair the result: at a co-split triple
on a fibre already known to have n split covers, the n-3 remaining native
covers all lift by definition. At 3/8 that gives 22 formal fourth choices
per successful triple. It is a restatement of the original square census,
not an additional predictor or an independence certificate.

## Mechanisms and remaining implications, in priority order

1. **Solubility, unresolved:** simultaneous vanishing of global descent
   obstructions on classes whose rational lifts retain several independent
   quotient directions remains the target. A useful carrier must address
   both this prescribed-parameter global condition and survival under
   specialization. No tested native relation supplies that chain for the
   largest gains.
2. **Solubility, weak as a large-jump explanation here:** degree-one pair
   events and degree-six/eight triple events construct rational relations,
   but concentrate on moderate retained gains and impose dependencies.
   Higher total-degree schemes and other cover dictionaries remain open;
   this null result does not strengthen either as an explanation.
3. **Solubility, contradicted within the tested native model:** the proposed
   compression to a small collision-defect span fails the broader census:
   all seventeen nontrivial compatible blocks attain the maximum n-1.
   Many split covers, raw collision support size, and genus mostly tracking
   the cover count also fail to order the largest retained gains.
4. **Incidence, still required:** determine which pre-point classes produce
   genuinely additional Mordell--Weil directions, then identify a global
   vanishing condition on their common carrier. Rational relations provide
   upper constraints; they are not such an incidence lower bound. Exact
   native quotient accounting and a carrier covering the unexplained
   directions are missing from this test.
5. **Coverage, required before generalization:** complete comparable native
   dictionaries and exact parameter transports for 11952 and A1/MW16,
   including the historic +12...+14 and recent high-gain controls. The
   1,024-cover partial 11952 census cannot support a full-dictionary null
   claim. No ordinary exact-rank control is created by a search miss.

For Agent 1 the actionable information is presently a restriction:
native splitting is a **solubility witness**, and these signed relations
can identify redundant constructions, but neither event count should be
promoted to a rank score. The current calculation has no visibility
endpoint and changes no selector. Repeatedly enlarging the degree bound
on the same selected quartets has lower priority than closing the family
coverage gap and testing whether the relevant carriers cover the retained
quotient directions at all.

## Replay

The frozen producer is `low_degree_triple_panel.py`; it uses one worker,
checkpointed parity cosets, a 180-second coset-batch cap, and a 60-second
cap per nontrivial fibre. Every worker completed. Its existing immutable
output should not be overwritten; a fresh producer run belongs in a
versioned location. Verification and report replay are:

```bash
/home/royvanrijn/.local/bin/sage -python elliptic-curves/rank-jump/verify_low_degree_triple_panel.py check
python3 elliptic-curves/rank-jump/report_low_degree_triple_panel.py check
```

Primary short-vector enumeration uses PARI; the independent enumeration
uses exact rational LDL bounds. Primary group equations use Sage elliptic
points; the independent replay uses a separate explicit rational group
law, including every negative branch test and canonical signed relation.
The geometric degree formula and generic lift data remain bound to their
earlier certificates. The environment is Sage 10.9 / PARI 2.17.3.
