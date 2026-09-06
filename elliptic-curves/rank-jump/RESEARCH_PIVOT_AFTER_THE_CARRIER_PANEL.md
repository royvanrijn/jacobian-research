# Research pivot: explain the unaccounted directions and the rational/Sha switch

The working goal remains an explanation of large specialization jumps.
The native-cover line should now be a regression and an occasional tool,
not the default source of the next experiment. Its recent negative results
do not justify automatically increasing the degree, number of charts, or
size of the generic atlas.

This decision follows the user's request to reconsider approaches and
goals, the [structural reassessment](../notes/RANK_JUMP_REASSESSMENT_2026-09-05.md),
and a new exact coverage audit. It does not change any active search policy
or mathematical-status entry.

## A stronger reason to change direction than a missing correlation

Let V=(E_t(Q) tensor Q)/M_t, where M_t is the rational span of the marked
generic subgroup. Let L be the certified retained j-dimensional quotient,
and let N be the span of the n displayed native sections at the fibre.
If the verified relation rows among those sections have rank c, then

\[
\dim N\le n-c,\qquad
\boxed{\dim L/(L\cap N)\ge\max(0,j-n+c).}
\]

This does not assume N lies in L. It follows directly from
dim(L intersect N)<=dim N. Thus it gives an exact lower bound on how much
of the retained quotient any explanation confined to these marked sections
must leave out, even without computing their full membership in L.

The [coverage certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_marked_carrier_coverage_gate_v1.json)
joins the frozen square census and independently verified pair/triple
relation ranks. It finds:

| Dictionary | Addresses with positive retained gain | Marked span provably cannot cover the retained quotient |
|---|---:|---:|
| Complete published-R17 native bisection atlas | 69 | **62** |
| Partial 1,024-cover native 11952 atlas | 57 | **57** |

These are family/parameter addresses, not asserted independent curves.
There are another 206 exposure observations with no valid dictionary or
transport in this audit. Deficit zero means NOT_EXCLUDED, not proved
coverage. The full curve ranks remain UNKNOWN.

For the largest published controls the deficits are especially clear:

| t | Retained gain | Number of marked native directions available | Retained quotient dimension outside their span, at least |
|---|---:|---:|---:|
| -308/251 | +9 | 3 | **6** |
| 2456/135 | +10 | 2 | **8** |
| -9529/5471 | +11 | 1 | **10** |

This is not a rank bound on the full generic group after a native base
change: undisplayed sections in native or mixed characters may exist.
It does rule out explaining the entire retained jump by reorganizing the
displayed native points, their collision primes, or their relations.
Expanding another family's atlas is worthwhile only with a specific reason
to expect it to reach the missing quotient directions.

## Lessons that change the next experiment

- **A rational construction can mainly create dependencies.** The
  degree-one pair events occur on moderate gains; all eleven degree-six/
  eight triple incidences occur at the +4 control. Their signed equations
  reduce the dimension available to the marked directions. Carrier
  solubility and survival of independent directions need separate proofs.
- **The simplest compressed-defect hypothesis failed.** All seventeen
  nontrivial compatible native blocks have maximum defect span n-1.
  Collision support identifies the remaining lift test, but a short list
  of primes does not explain why its global zero class occurs.
- **The strongest arithmetic contrast is still the rational/Sha switch.**
  In the fixed-cubic controls, u=0 has twenty witnessed rational classes;
  u=-1 has full Selmer dimension only one smaller, relative to a shared
  unknown anchor excess, yet has a certified restricted CT rank sixteen.
  The [relative full-Selmer theorem](RELATIVE_FULL_SELMER_THEOREM.md)
  separates a small incidence change from a large global obstruction.
  It does not prove the total rank difference or transfer this mechanism
  to the original MW17/MW16 families.
- **Most observed-zero controls are not mathematical negatives.** A
  lower-bound rank panel cannot identify the incidence/solubility split
  by itself. More accurate observed counts do not fix that identification
  problem. At least some certified upper bounds are needed.
- **A common block is a hypothesis, not a required answer.** Independent
  rational classes automatically generate many soluble subspaces and use
  one descent algebra. Those facts do not distinguish a special common
  cause from a rare arithmetic tail. A block claim needs structure fixed
  before exceptional points and an appropriate arithmetic baseline.

## Alternative approaches, ranked by what they could resolve

### 1. Directly separate Selmer capacity and rational solubility

On an actual same-family high/low pair, the desired report is

\[
(\text{generic rank},\ \text{certified Selmer bound},\
\text{CT obstruction rank},\ \text{witnessed quotient rank}).
\]

**Incidence endpoint:** establish whether both fibres have enough residual
Selmer capacity to accommodate a large jump. **Solubility endpoint:**
determine whether a substantial part of that capacity is excluded by CT
or higher descent on the lower-gain fibre. Positive point-built CT-zero
rows are controls, not independent predictors.

The initial matched published-R17 pair -2300/843 and -1561/3133 remains
a relevant target: generic rank, coefficient scale and retained exposure
are comparable in the [original panel](ANALYSIS.md). Neither low exact
rank nor a complete descent is currently supplied. The older 3/8 and
Nagao-0001 equation-only probes failed before producing a full Selmer
basis; historical factor-hinted extreme probes stalled in class-group
work. There is no justification for restarting those unchanged jobs with
a larger timeout merely to claim activity.

The next computation in this approach must first identify a new source
of certified Selmer bounds or independently constructed residual classes
that bypasses that recorded bottleneck. A bounded feasibility result is
useful only if it determines which exact arithmetic stage can now close.

### 2. Make exclusion, rather than prediction, a useful intermediate goal

For E(Q)[2]=0, suppose U is a certified upper bound on dim Sel_2(E),
c is the rank of a certified restricted CT matrix, and g is the known
generic rank. Rational point classes lie in the radical of the full
pairing, whose rank is at least c. Therefore

\[
\operatorname{rank}E\le U-c,\qquad j\le U-c-g.
\]

This needs a valid upper envelope and genuine Selmer classes. A partial
Selmer span supplies a lower bound, not U; a CT radical alone supplies
neither. Full pairing computation is not necessary when a smaller exact
nondegenerate minor already suffices.

**Concrete intermediate outcome:** certify one same-family ordinary or
low-gain fibre despite its high search score, or certify an upper bound
below a specified large-jump threshold. That would provide an uncensored
negative control and information Agent 1 can eventually use. It is a
useful result even if no sufficient high-jump condition has been found.
This remains a subgoal of the existing research objective, not a change
to the active rank search or a claim that such a certificate is available.

### 3. Explain a controlled multi-class rational/Sha switch first

The fixed-cubic comparison supplies a rigorously identified common
cohomology space and an exact relative incidence calculation. It is a
better setting for deriving why obstructions change than unrelated
high/low fields with no class identification.

**Solubility endpoint:** derive a normalization-invariant arithmetic
condition controlling a nontrivial part of the CT change; then explicitly
close the remaining global obstruction on a small independent block.
Reproducing a pairing matrix is not this endpoint. Nor does vanishing CT
exclude higher-divisible Sha. Existing norm/Jacobian constructions and
their Sha counterexamples are the required positive and negative
regressions.

The benefit is a possible theorem about simultaneous solubility. The
limitation is transfer: the fixed-field pencil has generic rank zero and
does not preserve the original MW17/MW16 structure. A successful theorem
must subsequently explain an independently selected block in those
families before being called a large-jump mechanism.

### 4. Test whether a special block is needed at all

Compete a discrete common-construction hypothesis with an arithmetic
rare-tail explanation. Both must be evaluated on quantities defined
without exceptional point input, conditional on family and on the local/
Selmer capacity actually known.

The present observed-rank panel cannot decide this: low ranks are censored,
the high controls were selected for high rank, and most full Selmer groups
are unknown. Soluble subspace counts conditioned on the witnessed rank
would be largely tautological. The first requirement is better arithmetic
labels from approaches 1 and 2, not another fitted correlation. Random
alternating matrices remain an illustrative null, not an asserted law for
these particular fibres.

### 5. Retain visibility as a separate, narrower objective

The initial half-lattice burst is worth understanding as a recovery
phenomenon. Existing exact chart incidence and height-compression work
already explains substantial clustering. New visibility work should
answer a precise unresolved exposure question from a completed transcript.
It should not consume most of the incidence/solubility effort or produce
a rank selector. The current turn starts no new chart computation.

## Gates against another unproductive local extension

For the next substantial mechanism experiment, state before computation:

1. Which independent quotient directions or which point-blind candidate
   classes it concerns, and how the generic subgroup is removed.
2. Whether its endpoint is incidence, solubility, or visibility.
3. What positive and negative outcomes would change the next decision.
4. Why the required arithmetic is newly feasible, given the frozen failed
   descents and negative carrier panels.
5. Whether a positive result proves independence, merely allows it, or
   imposes relations that reduce it.

Until one of these approaches supplies stronger evidence, the main effort
should not return to larger native relation degree bounds. There is no
validated positive rank feature to hand to Agent 1. The immediate handoff
is the exact coverage exclusion and the priority of obtaining a certified
arithmetic negative control or a genuine rational/Sha switching theorem.

## Scope and replay of this turn

The new arithmetic is a finite dimension audit over the already verified
326-address panel. It creates no points, parameters, local witnesses,
Selmer groups, rank claims, or changes to live search outputs.

```bash
python3 elliptic-curves/rank-jump/carrier_coverage_gate.py check
```

The producer binds the square census, its verification, the combined
relation report, and its own source. Existing sources consulted for this
decision include the original analysis, the relative full-Selmer theorem,
the structural reassessment, the frozen exceptional descent feasibility
report, and the two incomplete equation-only descent probes. Their
unresolved conclusions are preserved.
