# X1092: Elkies-style pair search and certified M19-to-V3 handoff

**Implementation, not a claimed X1092 rank result.** The code is committed for
execution with Sage/PARI. Eighteen local Python/SymPy regressions pass, including
an independent nontorsion proof for Elkies's published X948 reference pair.
Those tests do not execute X1092 arithmetic or V3. The bounded GitHub Actions
smoke records their actual endpoints separately; an implementation commit,
passing workflow, generic M19 claim or rational carrier point is not itself an
X1092 specialized M19 certificate. No `MATH_STATUS.json` claim changes here.

This is a parallel, opt-in lane. It does not stop, replace, resume or change any
existing campaign, scoring rule, seed packet, theorem gate or V3 numerical code.

## Mathematical route

Reuse the certified reduced X1092 equation, its full generic M17 basis, the
complete norm-ten translation-orbit table and the Euclidean bisection formula.
No known exceptional point, rank label, historical successful parameter or
record equation chooses a pair. Inventory equations are used only for exclusions.

For `C_w = 2O + 4F + phi(w)`, the intersection is `8 - w^t G v`.
The fixed cheap-conic window is ordered by stored word size, then orbit mask.
After constructing its equations, try same-squareclass pairs first, with a
separate bounded allowance that cannot consume the odd-intersection allowance.
Keep full rational scalars in every quadratic.

* **Same cover:** if `q_v = c^2 q_w`, transport both actual points to that same
  rational double cover. Conjugate sheets or two expressions are not two
  directions. A fixed five-address probe must certify all nineteen images with
  both existing finite-rank implementations before generic M19 is claimed.
  Failure is UNKNOWN, then the ordinary pair route proceeds.
* **Distinct covers:** require disjoint branch pairs and odd positive
  intersection. Parametrize the first conic and form the actual second-cover
  quartic, not just its Jacobian or a product-curve quotient. Construct a
  rational origin from the actual odd intersection divisor; polynomial norm
  reduction lowers its degree to one. A rational point at infinity is handled
  by reversing the chart. Check pointed birational maps in the function field.
* **Positive base rank:** test points supplied by the two deck involutions.
  Good-reduction group orders bound rational torsion; a nonzero multiple by
  their gcd proves nontorsion. No full rank, class-group, BSD or generator
  calculation is required. A pair whose bounded witnesses fail is UNKNOWN.

On a distinct-cover compositum, the nineteen-section Gram has inherited block
`4G`, cross columns `2Gw,2Gv`, new diagonal entries `16,16` and off-diagonal
`w^t G v`. Its Schur complement is `diag(6,6)` and determinant
`4^17 * 1092 * 36`. Each conic's generic map, trace attachment, branch avoidance
and height-eight lift is checked. The two independent quadratic characters
then certify a generic M19 subgroup on the genus-one base, not on the K3's
original rational base.

This follows the double-cover construction in Elkies, *An elliptic K3 surface
X/Q(t) with Mordell-Weil rank 17, I*, arXiv:2608.25406v1, especially section 4,
pp. 8-9: https://arxiv.org/abs/2608.25406 . Positive genus-one rank remains a
separate proof from the odd-degree rational-point argument.

## Specialization, scoring and V3

Enumerate a fixed signed-multiple window on the certified genus-one base, or a
fixed rational-address window for a successful same-cover pair. All maps are
rechecked. Keep exclusions for chart poles, singular fibres, explicit height
caps and exact Q-isomorphism duplicates. Equal j alone does not remove twists.
A fixed generator window is not uniform sampling of rational points or heights.

Choose score-independent hash controls **before scoring**. Compute the existing
Mestre-style good-prime score through 997, then extend the shortlisted candidates
and controls through 32749. Bad-reduction terms stay missing, not zero;
parameter-specific removable `p^4/p^6` factors are removed first. Cache repeated
finite curves and cross-check every distinct small-prime curve by a separate
character sum. Save every trace used. No prime at or beyond 65537 is read.
These are heuristic scheduling sums, not rank estimates or bounds.

Every selected fibre must supply the literal seventeen generic images plus the
two explicit lifted points. `FinitePointAdmission`, `memory_rank_certificate`
and `parent_foundry_worker.verify` must independently close rank nineteen.
A failure remains UNKNOWN rather than rank eighteen or a dependent cover.

The existing `parent_foundry_worker` consumes the resulting packet with
`generic_rank=17`, `rank_lower_bound=19` and its original generic prefix intact.
It reuses the unchanged V3 bank, adaptive search, complete-cloud reconciliation
and two finite implementations. Carrier rank does not bypass the seed gate.
A completed V3 result is accepted only with its model, request, prefix, packet
and independent-replay hashes attached.

## Run and resume

From the repository root, with Linux, Sage 10.9, NumPy and `/usr/bin/gp`:

```sh
sage -python -m unittest discover -s research/elliptic-curves/tests \
  -p test_elkies_pair_search.py -v
sage -python research/elliptic-curves/cas/run_x1092_elkies_pairs.py self-test

sage -python research/elliptic-curves/cas/run_x1092_elkies_pairs.py run \
  --folder research/artifacts/local/elliptic-curves/x1092-elkies-pairs-v1
```

The default frozen experiment constructs at most 64 conics, attempts at most
128 same-cover and 128 odd-intersection pairs, enumerates 128 addresses, extends
32 score-selected addresses plus two independent controls, and selects eight
ranked fibres plus those controls. Each admitted M19 receives at most 100 V3
calls. One worker, 3 GiB process-tree RSS, 60-second individual conic attempts,
90-second individual pair/probe attempts, 1800-second mathematical phases and
1800-second V3 jobs are explicit bounds. A bounded stage failure may stop the
campaign before the full finite window is exhausted.

The `prepare` subcommand freezes sources, generic-only inputs, policy and
executable hashes without running a search. `run` prepares only a new folder;
`resume` uses the sealed plan, not new command-line policy options. Completed
phases replay their output hashes. An unreceipted interrupted phase is UNKNOWN
and is not automatically rerun. Preserve it and use a fresh campaign after
investigation. A `STOP` file in the campaign folder prevents the next phase/job;
it does not terminate an already supervised call.

```sh
sage -python research/elliptic-curves/cas/run_x1092_elkies_pairs.py status \
  --folder research/artifacts/local/elliptic-curves/x1092-elkies-pairs-v1
sage -python research/elliptic-curves/cas/run_x1092_elkies_pairs.py resume \
  --folder research/artifacts/local/elliptic-curves/x1092-elkies-pairs-v1
```

The status command prints the latest saved report, not a live process check.
`plan.json`, `manifest.json`, phase supervisor logs/receipts and
`runtime/research/pair-results/` retain conics, selected pairs, carrier maps,
nontorsion witnesses, trace tables, seed proofs and full V3 outputs. A clean
controller exit does not mean a carrier or M19 was found: inspect these gates.

The separate GitHub Actions smoke uses the official Sage 10.9 image, one worker
and smaller fixed limits. It runs only manually or when its own workflow file
changes on main; it is not a scheduled production search. It preserves evidence
as an artifact, including failures, and has read-only repository permission.
