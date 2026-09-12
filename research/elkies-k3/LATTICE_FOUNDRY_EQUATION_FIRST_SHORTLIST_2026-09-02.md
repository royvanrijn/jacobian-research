# Historical equation-first experiments and the NS0031 local precursor

The different-NS foundry now applies the
[full rational-marking gate](DIFFERENT_NS_ARITHMETIC_MW17_FOUNDRY_OBJECTIVE_2026-09-04.md)
before source or equation searches. This note preserves completed local and
lattice comparisons. It is not a live shortlist or a campaign runbook.
Determinants 500, 720, 750, 950/NS0024 and 1184/NS0031 are arithmetically
excluded by the later exact marking theorems. Their local successes remain
useful regression controls.

The [original 1,227-line report](../archive/repository-cleanup-2026-09-12/research__elkies-k3__LATTICE_FOUNDRY_EQUATION_FIRST_SHORTLIST_2026-09-02.md.txt)
retains every historical table, SHA-256, comparison and command. Its dated
promotions and “next” instructions have been superseded.

## Exact NS0031 finite-field and finite-lift certificate

The source key is the pair
`(elkies-k3-lattice-foundry-prescribed-root-sources-all-ns-3e8-all-a-v1.json, NS0031-S001)`.
The identifier alone is insufficient because it is local to that source shard.
This semistable `A1+2A7/MW2` source has height Gram
`[[2,1],[1,41/8]]` and complete-basis pole profile `[0,1]`.

The source-pinned scan exhausts the declared generator charts in both twist
classes at 5 and 7. Both GF(5) charts are empty at the full marking gate.
At GF(7), the nonsquare chart has 36 pole-zero sections on 15 models and no
complementary pole-one section. The square chart has 20 pole-zero sections
on eight models, six correctly component-marked pole-one sections and two
complete marked pairs, both on normalized model 157. Their smooth mutual
intersection is two and their Shioda pairing is one.

At one pair, all 59 normalized fibre, section and component-jet equations
in 52 variables vanish over GF(7). The Jacobian has rank 51 and an explicit
maximal minor equal to one modulo seven. Deterministic Newton corrections
give 52 integer coordinates solving all 59 equations modulo
`7^8 = 5,764,801`. This is the exact scope of
`EC-K3-NS0031-MARKED-SOURCE-PRECURSOR`.

| Evidence | Retained source |
|---|---|
| Complete finite-field marking | [Scanner](scripts/scan_lattice_foundry_ns0031_a1_2a7_marking_modp.sage) · [GF(7) certificate](../artifacts/generated-results/elkies-k3-lattice-foundry-ns0031-a1-2a7-marking-mod7-v1.json) |
| Jacobian and finite lift | [Checker](scripts/certify_lattice_foundry_ns0031_marked_gf7_hensel.sage) · [Coordinates and residuals](../artifacts/generated-results/elkies-k3-lattice-foundry-ns0031-marked-gf7-hensel-v1.json) |
| Infinite compatible formal branch | [Separate formal proof](NS0031_MARKED_FORMAL_BRANCH_2026-09-04.md) |
| Rational marking excluded | [Exact modular obstruction](NS0031_QQ_MARKING_OBSTRUCTION_2026-09-04.md) |

A unit minor for an overdetermined system alone does not prove that its
omitted equations vanish along an infinite lift. The later formal theorem
closes that local obligation with the exact identity

```text
8 A^3 F = D C^4 (H-B C^2) - 9 B H^2 C^2 + H^3,
D=4 A^3+27 B^2, H=2 A X+3 B C^2, F=X^3+A X C^4+B C^6.
```

The prescribed fibre/component orders make all eight omitted residual rows
consequences of the retained 51 equations on the stated unit chart. Neither
the finite lift nor the formal theorem gives a rational marking. The later
global obstruction excludes an exact NS0031 marking over QQ.

The completed bounded rational-coordinate scan fixed `m9=n/d`, with
`|n|<=40`, `1<=d<=40`, coprime `n,d`, `7` not dividing `d`, and `n/d=1 mod 7`.
All 247 values lifted to `7^40` and returned `NO_FULL_RR`. This is a bounded
miss, not a proof of irrationality; the independent marking obstruction is
the reason not to extend this source route.

## Results to reuse before another calculation

| Historical comparison | Retained lesson and boundary |
|---|---|
| NS0028 and NS0005 normalized charts | Individual generators can occur without the required pair. Check their common model, physical components and mutual height; chart failures do not exclude characteristic-zero sources. |
| NS0031 complete degree-three census | Five selected frames already have all `3^17` cosets enumerated. F017 leads trisections, while the bisection-first scorer selects F018. Reuse the census; neither ordering proves a rational marking or rank jump. |
| Determinant720 ideal source cut | The 48 low-pole rows collapse to three marked isometry classes, all tested. The rational `s6=10` model saturates by index six to determinant20; its displayed determinant720 subgroup is not the full NS lattice. |
| Determinant500 MW1 source | The small rational reconstruction has a section divisible by five and saturates to determinant20. Its formal local branch and rootless lattice are valid, but the later exact marking theorem excludes determinant500 over QQ. |
| Determinants384,654,714 corridor beams | The stored beams missed their named targets at their declared caps. Rootful targets require recognition at root ranks two or one; the repaired terminal recognizer found no hidden old hit. |
| Determinant714 rootful spectrum | Use `(M/dM)/W(A1)` for the A1 target. The rootless translation formula cannot be copied unchanged. Its formal branch still supplies no rational marking. |
| Determinants750,864,1296,1500,1728 ideal source cut | Nonprimitive terminal embeddings left no source in the selected six ambients and MW0–2/support window. This is a scoped source-cut result; coarse genus-zero labels do not establish rational marked moduli. Determinant750 has a separate later arithmetic exclusion. |

Source rank, full-basis pole cost, target multisection counts and rational
marking are different coordinates. A lower source rank or a richer target
census cannot replace the marking gate. Keep exact isometry classes distinct
from reduced-Gram duplicates and preserve the physical marked basis.

## Replay boundary

All scripts and inputs remain at their certificate paths. Historical commands
are in the preserved report. Several `--check` modes repeat discovery:
the rational-parameter checker starts 247 Sage lifts before comparison, and
the degree-three census visits 129,140,163 cosets per frame. They are not
navigation checks and were not rerun for this cleanup.

For the current equation-dispatch decision use the
[arithmetic classifier](RANK19_ARITHMETIC_MARKING_CLASSIFIER_2026-09-04.md)
and [different-NS objective](DIFFERENT_NS_ARITHMETIC_MW17_FOUNDRY_OBJECTIVE_2026-09-04.md).

<!-- status-consumer: EC-K3-NS0024-QQ-MARKING-OBSTRUCTION b7f0cf002c0411fe -->
<!-- status-consumer: EC-K3-NS0031-MARKED-SOURCE-PRECURSOR 2e115b35c30a8cea -->
<!-- status-consumer: EC-K3-NS0031-MARKED-FORMAL-BRANCH b31e99bce4edac0a -->
<!-- status-consumer: EC-K3-NS0031-MARKED-RATIONAL-PARAMETER-SCAN ca678e520745dd3c -->
