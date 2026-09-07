# Parent expansion and a generic-section certification correction

The current201-curve high-rank inventory is concentrated on one underlying
determinant-948 K3 surface. It uses twelve family labels: six compact R17
inputs, five compact A1/MW16 inputs and the published-R17 presentation.
The exact fibration-hop and atlas provenance, rather than equality of lattice
determinants alone, identifies that common surface. New parameter values,
different fibrations, carrier pullbacks and point-search coordinates must
not be counted as unrelated parent surfaces.

This follow-up supplies **six additional Q-isomorphism classes of K3 parents**,
all pairwise distinct and all distinct from X948. They are inputs from an
existing construction, not a claim of unpublished surfaces or record curves.
The intake itself performs no elliptic point search or score sweep. A subsequent
[fixed calibration](MESTRE_PARENT_CALIBRATION_2026-09-07.md) now completes588
point boxes on twelve fibres, with two subgroup gains from11 to12.

The later [exact lattice proof](MESTRE_DETERMINANT468_PARENTS_2026-09-07.md)
establishes geometric NS rank19 and determinant468, rational NS rank18,
and generic MW ranks11/12 over Q(T)/Qbar(T) on all six. This closes the
different-geometric-determinant question left open by this intake. A useful
rank16 fibration remains unconstructed; arithmetic rank17 is impossible
on these parents. Neither statement bounds specialized elliptic ranks.

## The outer parameter changes the parent

The exact Mestre–Fermigier two-section component has outer parameter u,
with v=(u²+u+2)/u; its six roots define a quartic family in the separate
fibre parameter T. Freeze u=11,13,17,19,23,29, beyond the previous H(u)≤10
numerical triage, and use only T=1 for the subgroup intake. No candidate
replacement or parameter search is performed.

For each u, symbolic square approximation reconstructs the quartic from
the six labelled roots and verifies the two extra rational polynomial
ordinates. Its Jacobian has coefficient degrees8,12. The exact discriminant
has degree20 with squarefree multiplicity profile16 simple roots and two
double roots; c4 is coprime to it. Infinity has split multiplicative type I4.
Thus the resolved surface is a K3 with fibre types16 I1 +2 I2 +I4.
This differs from the initial intake's assumed20 I1 +I4 pattern. That
first six-row attempt failed its geometry gate before rank intake; its
sources, protocols and failures remain preserved.

Each subsequent T=1 intake checks the fourteen supplied covariant images.
It exactly halves their differences from the first image to obtain rational
divisor-section specializations. Adding the first image back does not
increase the certified finite rank. The14-point seed clouds each certify
**rank at least11**, independently confirmed by complete finite-group
enumeration. This is a lower bound for each parent and its T=1 fibre;
neither exact generic rank nor absence of further directions is claimed.

## Distinct surfaces, not merely different equations

The first three common good primes in the fixed interval101..251 are
131,239,251. The checker preserves discriminant multiplicities and
semistability at these primes and verifies the split infinity node.
It counts the smooth projective surfaces, including the minimal-resolution
corrections: +p at each rational I2 base value and +3p at the split I4
infinity fibre. Every smooth elliptic fibre count is computed both by a
character sum and independently by Sage cardinality.

| Outer u / parent | #S(F131) | #S(F239) | #S(F251) | Certified section lower bound |
|---|---:|---:|---:|---:|
| 11 | 19266 | 61074 | 67410 | 11 |
| 13 | 19403 | 60956 | 67356 | 11 |
| 17 | 19403 | 60956 | 67410 | 11 |
| 19 | 19403 | 61074 | 67050 | 11 |
| 23 | 19403 | 61524 | 67410 | 11 |
| 29 | 19403 | 60956 | 67403 | 11 |
| Existing X948 | 19308 | 61428 | 68022 | Existing MW17 / MW16 fibrations |

Every one of the21 pairs differs at a common good prime. Isomorphic smooth
projective K3 surfaces over Q have the same good-reduction cohomological
traces, hence the same counts. These witnesses prove distinct Q-isomorphism
classes; for K3 surfaces they also exclude Q-birational equivalence.
They do not distinguish geometric isomorphism classes over Qbar or prove
different Néron–Severi lattices. This is not the different-NS MW17 foundry
milestone. The [parent certificate](../../artifacts/generated-results/elliptic-curves/mestre_parent_portfolio_intake_v1.json)
retains the equations, subgroup proofs, all counts and separating primes.

## The old generic-rank13 argument mixed its columns

`EC-MF2S13` previously asserted thirteen independent sections by stacking
finite quotients at four different (u,T) values. Two implementation choices
made its column labels inconsistent across those values:

- `SixRootMestreConstruction` sorts roots by their rational values at every
  specialization. The sorted positions change as u changes.
- Both extra ordinates use a positive rational square root at each fibre.
  Positive values need not specialize one rational-function square-root branch.

A relation among fixed generic sections specializes with the same integer
coefficients. These changing labels/signs broke that premise. The old
finite calculations themselves still replay, but cannot prove the claimed
generic independence.

The correction fixes the source-root permutation at u=-5, derives exact
global rational ordinate polynomials in Q(u)[T], and verifies their square
identities. At u=-3, the old visible positions6/7 and8/9 swap, and the
second extra ordinate has the opposite sign. At u=-1/2 a larger root
permutation occurs and the second sign again changes. With these corrections,
the **same13 probes have rank11 instead of13**. Their torsion exclusion
remains valid, so coherent independent columns prove generic lower bound11.
This audit does not prove the full generic rank equals11.

The [label certificate](../../artifacts/generated-results/elliptic-curves/mestre_component_label_audit_v1.json)
records both matrices and every permutation/sign change. A
[separate finite-group verifier](../../artifacts/generated-results/elliptic-curves/mestre_parent_and_label_independent_v1.json)
enumerates complete groups and quotients by triples, independently reproducing
the13-to11 change, and verifies all six new seed clouds by quotients by doubles.
`EC-MF2S13` is now partial; its former rank13 proof is rejected. The symbolic
component identities and unrelated single-specialization rank13 certificates
remain valid.

## Search integration

[`mestre_parent_adapter.py`](../cas/mestre_parent_adapter.py) exposes
`specialize(outer_u, fibre_T)` with a fixed generic section order and exact
global ordinate branches. Its portable rational-function coefficients are
[exported separately](../../artifacts/generated-results/elliptic-curves/mestre_component_coherent_sections_v1.json).
Regression checks cover all thirteen old probes, all six new T=1 inputs and
parameter poles. This adapter supplies equations and points; it does not
assume they are independent. Each new fibre must pass its own subgroup check.

The lower generic ranks mean these parents need a larger specialization jump
to compete with the MW16/MW17 portfolio. Their small T=1 displayed equations
and genuinely different parent surfaces justify calibration of point
visibility and arithmetic scoring before substantial searches. No superiority
claim follows from parent count alone, and no automatic large campaign is queued.

The preceding retained-rank27 coordinate audit also completed: all343
factor-free maps are exact, and each has an explicit rational coordinate
inside its new125000 box but outside its old box. Such a coordinate need
not be a quartic square. It launched zero point boxes and adds no parent;
the user-directed parent expansion took priority over a follow-up search.

## Replay

```sh
sage -python elliptic-curves/cas/verify_mestre_parent_and_label_audits.sage --check
sage -python elliptic-curves/cas/export_mestre_coherent_sections.sage --check
python3 -m unittest discover -s elliptic-curves/tests -p test_mestre_parent_adapter.py
python3 elliptic-curves/cas/factor_free_rank27_box_audit.py check
```

The immutable construction scripts write new outputs only. Raw parent
protocols, failed V1 geometry gates, V2 intakes and supervised audits are under
`artifacts/local/elliptic-curves/mestre-parent-portfolio-intake-v1` and `v2`.
The earlier map-only audit is under `factor-free-rank27-box-audit-v1`.
