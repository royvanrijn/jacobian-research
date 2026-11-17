# Complete omitted-divisor accounting

The ledger separates ramification, omission, target image, residue degree,
pole order, normalization index and exceptional ambient length. Current
additional affine image equations remain unknown; no absent curve is inferred
from a missing input or from identity monodromy.

The [remaining input audit](INPUT_GATE.md) specifies the concrete global
data needed to turn the conditional rows into complete image and fibre
accounting, and the normalization calculations those data enable.

| Result | Canonical proof and exact data |
|---|---|
| PF2OMIT19: central image-fibre Euler/Milnor, normalized missing blocks and retained residue covers | [Proof](CENTRAL_IMAGE_FIBRE.md) · [ledger](central_image_ledger.json) · [replay](central_image_receipt.json) |
| PF2OMIT18: nonlinear image-pencil geometric monodromy, connected source and every generic omitted point | [Proof](IMAGE_PENCIL_MONODROMY.md) · [ledger](image_pencil_ledger.json) · [replay](image_pencil_receipt.json) |
| PF2OMIT17: universal additional-image intersection, connected support and mandatory singleton loci | [Proof](UNIVERSAL_IMAGE_INTERSECTION.md) · [ledger](critical_trace_ledger.json) · [replay](critical_trace_receipt.json) |
| PF2OMIT16: affine pole module, disjoint-image and support-connectedness gates | [Proof](DISJOINT_IMAGE_GATE.md) · [ledger](disjoint_image_ledger.json) · [replay](disjoint_image_receipt.json) |
| PF2OMIT15: odd-degree images force affine intersections and singleton packets | [Proof](ODD_IMAGE_INTERSECTIONS.md) · [ledger](odd_intersection_ledger.json) · [replay](odd_intersection_receipt.json) |
| PF2OMIT14: generic/radial genus, distinct puncture and sheet charges, projection Euler jumps and the full two-family catalogue | [Proof](GENERIC_LINE_ACCOUNTING.md) · [ledger](generic_line_ledger.json) · [replay](generic_line_receipt.json) |
| PF2OMIT13: full first-image corank correction, explicit cusp values and the fourth oriented r7 row | [Proof](NONIMMERSED_LEDGER.md) · [ledger](nonimmersed_ledger.json) · [replay](nonimmersed_receipt.json) |
| PF2OMIT12: no quartic common divisor; every image-parametrization fibre has length at most3 | [Proof](NO_FOUR_BRANCHES.md) · [certificate](four_branch_certificate.json) · [ledger](four_branch_ledger.json) · [replay](four_branch_receipt.json) |
| PF2OMIT11: explicit first-image equations in both r7 charts and complete singular-value matching | [Proof](IMAGE_NORM_LEDGER.md) · [matrices/ledger](image_norm_ledger.json) · [replay](image_norm_receipt.json) |
| PF2OMIT10: every singleton packet's inverse coefficient and normalization-index orders | [Proof](SINGLETON_INVERSE_VALUATIONS.md) · [ledger](singleton_valuation_ledger.json) · [replay](singleton_valuation_receipt.json) |
| PF2OMIT9: marked Wright image jets, quartic cusp value and compulsory inverse-index zero | [Proof](WRIGHT_IMAGE_RESTRICTION.md) · [ledger](wright_image_ledger.json) · [replay](wright_image_receipt.json) |
| PF2OMIT8: residual S6 monodromy, primitive source field and complete generic inertia/deletion rows | [Proof](GLOBAL_MONODROMY.md) · [ledger](global_monodromy_ledger.json) · [replay](global_monodromy_receipt.json) |
| PF2OMIT7: simple-boundary regularity, ordinary cusp and complete singleton radial packets | [Proof](SIMPLE_RAMIFICATION_AND_RADIAL.md) · [ledger](simple_radial_ledger.json) · [replay](simple_radial_receipt.json) |
| PF2OMIT6: all residual omitted f1 and pointwise fibre rules | [Proof](RESIDUE_LEDGER.md) · [ledger](residue_ledger.json) · [replay](residue_receipt.json) |
| PF2OMIT5: entire single-prime gap21 stratum excluded without finality | [Proof](R7_INFINITY_COMPLETION.md) · [ledger](infinity_ledger.json) · [replay](infinity_receipt.json) |
| PF2OMIT4: normalized discriminant versus inverse-order index | [Proof](DISCRIMINANT_LEDGER.md) · [ledger](discriminant_ledger.json) · [replay](discriminant_receipt.json) |
| PF2OMIT3: source-finite inverse and exact specialization lengths | [Proof](ADMISSIBLE_INVERSE.md) · [equations](inverse_equations.json) · [replay](inverse_receipt.json) |
| PF2OMIT2: all-prime Euler identity and global counts | [Proof](GLOBAL_LEDGER.md) · [ledger](global_ledger.json) · [replay](global_receipt.json) |
| PF2OMIT1: reconstruction poles and ambient lengths | [Proof](PROOF.md) · [ledger](ledger.json) · [replay](receipt.json) |

Under the immersed r7 degree-six hypothesis with additional ramification,
the other ramified row has (e,E,t)=(2,0,1),(2,1,0),or(3,0,0). All omitted
residue degrees are one. All three affine normalizations are smooth;
the excess-one point is an ordinary cusp with algebra C[w]/w^3. Every
residual cover has monodromy S6 and no proper intermediate source field.
The sign double cover is only a field in its Galois closure.

The known r7 image now has two explicit six-by-six norm equations,
whose corank gives every normalization branch count. The second image
must avoid its singular values in the index-three row; every intersection
then has a singleton source fibre. The grouping H=C excludes thirteen
immersed polynomial-chart parameters having triple values. Actual B,H,
the chosen coefficient parameter and source realization remain unknown.

The forbidden four-branch Laurent locus is now proved empty: the image
parameter map has fibre scheme length at most3 in both complete charts,
including nonimmersed parameters. This rank bound removes no coefficient
parameter; the full affine collision catalogue is recorded below.

At either nonimmersed r7 class the unique cusp coordinates are explicit.
The matrix correction is c=corank-epsilon and lost length=2corank-epsilon.
The new oriented row has two omitted e2,f1 primes, excess1 on the known
first image and0 on the second, with no unramified omitted prime. Its
only possible singleton value is that cusp, exactly when B passes there.

Both affine first-image catalogues are complete:268 nonnodal coefficient
points,163 triple,99 tacnode,6 cusp. Generic target-line and generic
fixed-point radial pullbacks have genus>=3 (>=5 in the index-three row).
The global line ledger separates removed points, lost sheets and different;
a projection Euler jump3 can include an unramified omitted point.
For source degrees(75,125), the inverse pole budget bounds deg B by110,
by100 when H=C, and by55 when H=B. Every additional image now must
meet C, by universal coprime trace moments on both complete charts.
All residual supports are connected; shared-image A and cyclic C
necessarily have singleton values with packets2+2+1 and2+3. Actual
image equations and the normalized global source remain unknown.
The nonlinear H_C pencil has geometric generic group S6, or A6 in row C,
and an integral source generic fibre. Its genus and puncture formulas
use target genus12 and retain every unramified omission.
Its central fibre requires an infinity Milnor jump at least46. The
normalized retained cover records tangency/cusp deletions separately,
with clean-stratum residue partitions and split-section lift conditions.

On the actual marked Wright model, a hypothetical finite quartic pair has
one omitted e2,f1 prime and exactly one ordinary cusp value. Its source
fibre has length1 there, length2 over one immersed branch, length0 over a
double image, and length4 off the image. Restriction formulas locate the
cusp from any proposed pair. Its admissible inverse has reduced leading
coefficient H_C and a compulsory index zero at the cusp; the global pair
is still unknown.

The complete r7 singleton packets A0,A1,A2,B,C require inverse-index
orders at least10,6,4,4,2. Their coefficient and raw-discriminant bounds
retain the distinction between cusp and cyclic length-three factors.

Over a singleton zero, the complete omitted radial index lists are
1+1+1+1+1, 2+1+1+1, 2+2+1, or 3+2, each with f1 and target image the
exceptional P1. Special-direction points are counted once even where
radial and strict boundaries meet. Their lengths and canonical vectors
are explicit; actual affine image polynomials and exceptional target
coordinates remain unresolved. Earlier nulls fixed by later results
are retained as historical checkpoints.

The degree-eight marked-root control retains both ramified and unramified
omission, with zero fibre 8=1+6+1 and 27 pure inverse-coordinate collisions
on its elliptic line. It is not a plane Keller counterexample. The
[Stein scope audit](STEIN_SCOPE_AUDIT.md) keeps the old whole-even-row
claim PF2D6E1 partial because omitted e1 rows were not exhausted.

```sh
timeout 30s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify_central_image.py
timeout 30s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify_image_pencil.py
timeout 30s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify_critical_trace.py
timeout 30s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify_disjoint_images.py
timeout 30s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify_odd_intersections.py
timeout 30s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify_generic_line.py
timeout 30s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify_nonimmersed_ledger.py
timeout 30s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify_no_four_branches.py
timeout 30s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify_image_norm.py
timeout 30s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify_singleton_valuations.py
timeout 30s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify_wright_image.py
timeout 30s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify_global_monodromy.py
timeout 30s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify_simple_radial.py
timeout 30s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify_residue_ledger.py
timeout 30s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify_infinity_completion.py
timeout 30s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify_discriminant.py
timeout 30s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify_inverse.py
timeout 60s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify_global.py
timeout 60s .venv/bin/python research/jacobian/jc2-omitted-divisors/verify.py
```

[Current scoped validation](central_image_validation.json) ·
[Orevkov source](sources/residue_manifest.json) ·
[trace/order sources](sources/trace_manifest.json) ·
[normalization sources](sources/manifest.json) ·
[singleton-zero reduction](../../hessian/hc4-relative-row-orbit/PROOF.md).
Exact controls support written proofs and named imports. They certify
neither formal verification nor external review nor the full open goal.
