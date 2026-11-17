# Canonical receivers and uniform quotient obstructions

Latest extension: [finite universal cubic and quartic self-maps, torus-quotient rank ceiling and canonical volume lifting](../universal-keller-torus-quotient-selfmaps/README.md).

Earlier extension: [all-rank homogeneous receiver pairs, universal monodromy and exact factorial isotropy costs](../universal-keller-homogeneous-receivers/README.md).

Earlier extension: [canonical cubic and quartic torus receivers, exact torus classification and Galois-lattice stable atomicity](../universal-keller-torus-receivers/README.md).

Earlier extension: [canonical cyclic surface receivers, exact omissions, stable atomicity and stronger parity degree bounds](../universal-keller-acyclic-boundary/README.md).

Earlier extension: [the full-degree canonical obstruction in every dimension and exact modular affine exceptions](../universal-keller-fixed-fibres/README.md).

Further extension: [uniform profile atomicity, exceptional boundary classification and equivariant normal form](../universal-keller-structure/README.md).

The [proof](PROOF.md) constructs the full good-generator receiver in every rank N>=3, with trace-zero versions for N>=5 and two-trace versions for N>=6. It retains every choice as a canonical correspondence and glues the polynomial Keller maps on rank-three vector bundles.

Its geometric automorphisms are completely classified: trivial in even rank, and at most mu2 in odd rank. Hence it admits no natural marked-realization section from BS_N. This also excludes every enhancement that still maps to this receiver over BS_N.

Two separate results explain why taking the quotient does not yield a smaller affine Keller map. The global affine coarse quotient contracts one explicit affine line in every degree; normalizing target c gives a source C=1 with a nonconstant unit.

- UARECV1: full root-line receiver, polynomial bundle gluing and canonical degree-N root cover.
- UAINERT1: exact inertia and no sections for receivers retaining a good generator line.
- UAQFAIL1: global coarse-quotient contraction and singular vertices.
- UACTARGET1: affine source obstructions for target-dependent normalization, including even-rank charts.
- UAPROFILE1: exact rational-profile polynomiality, preserved Jacobian and characteristic-zero monodromy.
- UAEQUIV1: fixed-space restrictions and the necessary JC(2) gate for small natural representation bundles.

[Unified theory map](THEORY_MAP.md) · [Literature and gap map](LITERATURE_AND_GAPS.md) · [checker](verify.py) · [receipt](receipt.json) · [sources](sources.json) · [authority](../../MATH_STATUS.json).

    .venv/bin/python research/jacobian-followup/universal-keller-receiver/verify.py

These are construction-specific theorems with written proofs and exact controls. They do not give global dimension lower bounds or a global moduli stack of stable maps. Independent review remains pending.

The [canonical trace-filter receiver](../universal-keller-canonical-filter/PROOF.md) now constructs a different all-rank marked polynomial-bundle self-map without a primitive generator line. Its exact complete fibre is the input algebra, and its two nontrivial factors are explicit; canonical atomicity remains open.
