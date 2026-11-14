# Dixmier experiments connected to our recent research

14 November 2026. The best near-term direction is to make an existing
quantization obstruction survive a precisely larger class of corrections
and chart changes. The strongest conceptual alternative is to replace a
bounded correction calculation by a rational differential-operator theorem.
Both can use our recent boundary methods. Neither starts with a claimed
reduction of arbitrary DC2 to the known classical symbols.

## Two immediate findings

The previous review overlooked **QRF2**, which supersedes the cubic
continuation in OP-QO. The live ledger and the canonical residue note
already record a characteristic-zero obstruction on the complete reduced
affine-27 fifth-order lift component. Its X18 functional is independent
of all lift parameters and has an explicit inverse in the cubic field.
The scalar Bezout identity passes a fresh SymPy calculation. The original
401-equation Singular artifact is missing in both old and current artifact
locations, so full replay is unavailable. This is an assurance gap in a
recorded result, not an uncompleted modular-to-rational discovery task.

A new [exact ordering pilot](test_weyl_ordering.py) also shows why our raw
normal-order defects need further interpretation. Under symmetric Weyl
quantization, R and S commute exactly. The other five pairs retain defects:

| Pair | Nonzero powers of hbar in the symmetric commutator defect |
|---|---|
| R,T | 3 |
| R,D | 3 |
| R,S | none |
| T,D | 3,5,7 |
| T,S | 3,5 |
| D,S | 3,5,7 |

The complete series terminates by the total momentum-degree bound.
For R,T its order-three term is

\[
-\frac32 x^8(3qx-2)(3qx-1).
\]

This was checked both by the sparse Moyal engine and independently as
three quarters of the mixed derivative T_ppz. A nonzero defect can still be
removed by a correction; it is not yet a nonzero obstruction class.

## First priority is obstruction transport and support completeness

The requested first transport is now completed as
[QRT1](../dixmier-residue-transport/PROOF.md): arbitrary polynomial
positional shears preserve the strong residue on the transported category.
The nonlinear example also exposes a native correction of period
-1296 lambda. The remaining gate is a strong cocycle for the combined
transported and native spaces; the experiment below is retained as its
original motivation, rather than an unstarted task.

The existing quantum work gives strong dual cocycles, while our newer
[two-divisor graph result](../hc4-triangular-graphs/TWO_DIVISOR_EXTENSION.md)
shows how boundary traces can replace a word-length search by a theorem
covering an infinite class. The proposed quantum analogue is to transport
the correction complex, its boundary lattice and its dual functional
together under exact quantized canonical changes.

Start with one generic admitted symbol and its known nonzero cocycle.
Choose one nontrivial elementary canonical shear with a known inverse and
exact Weyl lift. Transport both the images and the allowed corrections;
compute the induced map on the full six-relation complex and verify that
the transported functional has the same nonzero value. Then compare the
transported support with the native support. Any discrepancy is an explicit
missing correction direction, rather than a reason to claim an invariant.

The first useful result would cover every polynomial shear in one declared
direction, proved by a symbolic shear parameter and an all-degree support
argument. Survival for a few numerical parameters or support cutoffs is
only exploratory evidence. Eventually the target is an obstruction for a
declared orbit and correction category, stronger than a single-symbol
rank jump.

Direct reduction modulo the classical boundary ideals is insufficient:
for f=x(1+xy), the Weyl relation gives

\[
[p_x,f]=\hbar(1+2xy),
\]

which vanishes on neither proof divisor. Those ideals are not preserved
by the quantum differentiation. Use shifted valuation lattices, as already
constructed in QRF1, and explicitly test differential stability instead.

## Second priority is an all-degree rational operator obstruction

[Kanel-Belov and Fedorov, Lemma 4.8](https://arxiv.org/html/2610.33919v1#S4.SS2)
exclude monic rational differential operators of orders two and three
with nonzero constant commutator in their specified differential field.
Their proof derives an equation W''+3W^2=lambda*t+mu and uses local and
infinity valuations. The preceding reduction from a graph module to
those operators is essential.

Our relative weighted families and boundary calculations suggest a test:
for one genuinely different admitted family, derive its localized scalar
operator orders and its coefficient differential field. Record the field's
normalization and all punctures before using a rational-function argument.
Compute the compatibility equations for the first new order pair, then
ask whether valuation balances force a contradiction independent of
coefficient degree.

This could exclude all rational lower-order coefficients in the identified
operator category, a stronger conclusion than another bounded order-seven
calculation. If an exact rational solution exists, the proposed operator
obstruction fails and that solution becomes a candidate for testing
polynomial descent. Neither outcome settles arbitrary Weyl endomorphisms.
The approach is a research proposal, not a theorem transferred from the
orders-two-and-three case.

## A tractable retained endpoint is the degree-nine family

The [degree-nine note](../../archive/non-elliptic/extended-geometry/DC2_DEGREE_NINE_MARKED_ROOT_SEARCH.md)
has an exact relative-rank calculation but only modular discovery of the
order-five section-zero locus. Degrees six through eight already have
their complete inherited order-seven exclusions; repeating them would not
advance the frontier.

For degree nine, inventory the retained polynomial and modular inputs
first. If available, reconstruct the characteristic-zero locus and verify
ideal containments and all boundary saturations exactly. Only then build
order-seven lift equations for its genuine components, preserving possible
nilpotents. A rational unit certificate or a genuine higher-order survivor
would be a concrete result. Missing artifacts do not authorize rebuilding
the old campaign. This endpoint is lower risk than a new structural
theorem, but remains restricted to the inherited filtration.

## A longer route is bounded arithmetic complexity

The October graph-lifting criterion suggests computing Berest–Wilson
invariants of actual plane graph-splitting modules as characteristic varies.
An initial implementation should recover the invariant from two independent
ideal presentations and normalize away the common factor of the associated
graded ideal. Automorphism examples provide controls, not evidence for
the difficult cases. Arbitrary rank-one ideals need not be graph modules;
higher-rank graph modules do not acquire this invariant by the same formula.

Our [new cusp-filtration rigidity](../jc2-cusp-filtration-rigidity/README.md)
could help control degrees once a compatible filtration exists. First
prove compatibility: its Poisson bracket has degree zero, so one cannot
assume it supplies a filtration in which commutators strictly lower degree.
Our [finite-target regularization result](../keller-solidity/PROOF.md)
also rules out trying to manufacture global module finiteness merely by
a fixed finite-order differential regularizer on a nonfinite source.
It supplies no noncommutative finiteness theorem automatically.

The concrete arithmetic target is a bound on presentations for a declared
graph-module family, uniform in characteristic under the paper's hypotheses.
Finite prime data can guide that proof but cannot replace it.

## Reproduction and priority

```sh
.venv/bin/python research/jacobian-followup/dixmier-audit/verify_long_pbw.py
.venv/bin/python research/jacobian-followup/dixmier-audit/test_weyl_ordering.py
```

[Ordering receipt](weyl_ordering_receipt.json) and
[independent scalar controls](pilot_controls.json) record the small pilots.
The ordering pilot uses complete exact rational coefficient support and
an a priori contraction bound, with no candidate search. No new claim is
promoted in MATH_STATUS.json.

Prioritize transport of a strong cocycle for a structural contribution;
the degree-nine locus is the nearer finite certificate target when its
retained inputs are available. The arithmetic-module route requires more
new infrastructure. Long's fixed arithmetic lift and the completed QRF2
cubic branch serve as negative controls rather than new discovery targets.
