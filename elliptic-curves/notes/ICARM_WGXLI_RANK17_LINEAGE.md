# Public `wgxli` rank-17 lineage candidates

Status: **historical fingerprint; the five-member family reconstruction is
now exact**.

The complete 43-chart norm-twelve sweep has identified curves 351, 356, 376,
377, and 385 as untwisted rational fibres of one eight-chart published-R17
class, reconstructed and saturated all seventeen generic sections, and
computed the displayed exceptional quotients.  See
[`../../elkies-k3/R17_NORM12_RECORD_LINEAGE_SWEEP_2026-09-04.md`](../../elkies-k3/R17_NORM12_RECORD_LINEAGE_SWEEP_2026-09-04.md).

<!-- status-consumer: EC-K3-R17-NORM12-RECORD-LINEAGE-ATLAS 00e39f6b05c2688a -->

The fingerprint and inverse-problem discussion below is retained as the
historical route to that exact result.

This note by itself records experimental construction evidence; the later
atlas certificate supplies the theorem that was still open here.

## Outcome

The 2026-09-01 ICARM snapshot contains three additional public candidates for
the labelled curve-351/356 lineage:

| curve | submitted UTC | rank lower bound | role |
| --- | --- | ---: | --- |
| 351 | 2026-08-25 13:52:22 | 25 | original target |
| 356 | 2026-08-25 14:20:06 | 29 | original target |
| 376 | 2026-08-26 13:40:43 | 22 | new candidate fibre |
| 377 | 2026-08-26 13:43:23 | 23 | new candidate fibre |
| 385 | 2026-08-27 02:47:30 | 29 | new candidate fibre |

All five were submitted by `wgxli`.  The ordered first seventeen public
points retain two nontrivial denominator-root anchors in every fibre:

```text
point position       5    17
sqrt(denominator)   71    41
```

Four of the five also retain denominator root `7` at position 16.  At
position 13, curves 351, 356, and 377 have root `679`, while curves 376 and
385 have root `97`.  The complete ordered vectors are stored in the artifact.

The scalar fits between the ordered seventeen-point canonical-height matrices
are:

| pair | fitted scale | relative Frobenius residual | Pearson correlation |
| --- | ---: | ---: | ---: |
| 351/356 | 1.420878 | 0.112201 | 0.974884 |
| 351/376 | 0.880750 | 0.181111 | 0.922835 |
| 351/377 | 0.876573 | 0.158751 | 0.941499 |
| 351/385 | 1.586907 | 0.147310 | 0.951100 |
| 356/376 | 0.621380 | 0.125181 | 0.957047 |
| 356/377 | 0.615095 | 0.136444 | 0.949747 |
| 356/385 | 1.115838 | 0.104980 | 0.969527 |
| 376/377 | 0.977911 | 0.164220 | 0.927683 |
| 376/385 | 1.773492 | 0.141496 | 0.945399 |
| 377/385 | 1.782051 | 0.165960 | 0.925706 |

Curve 385 is therefore a particularly strong third fibre: its residual
against curve 356 is smaller than the original 351/356 residual.

## Complete same-submitter sweep

The hash-pinned database snapshot has 474 curves, of which thirteen were
submitted by `wgxli`:

```text
351, 356, 363, 364, 376, 377, 378, 385, 389, 390, 391, 393, 395.
```

Join two records when the best scalar fit between their ordered first-17
height matrices has relative Frobenius residual at most `0.2`.  The connected
components are exactly

```text
{351,356,376,377,385}, {363,364,378}, {389,390,391}, {393}, {395}.
```

The largest within-target residual is `0.181111`; the smallest cross-component
residual is `0.348333`.  Thus the five-member component is not an artifact of
a finely tuned cutoff.  The second and third components also look like
separate labelled private lineages, but they are not mixed into the 351/356
inverse-interpolation target.

This graph is a bounded numerical classifier on one public snapshot.  It is
not a proof of family membership and is not claimed to recover deleted or
unpublished search output.

## Rootless-K3 interpolation input

For each of the five target curves, the artifact applies the exact canonical
short-model change

```text
X = 36*x + 3*b2,
Y = 216*y + 108*(a1*x+a3)
```

and stores the resulting equation

```text
Y^2 = X^3 - 27*c4*X - 54*c6
```

together with the first seventeen transported points.  This removes the
general-Weierstrass translations before interpolation.  A proposed rootless
K3 model must then admit parameters `t_k` and nonzero scalings `u_k` such that

```text
A_k      = u_k^4 A(t_k),        B_k      = u_k^6 B(t_k),
X_{i,k}  = u_k^2 x_i(t_k),      Y_{i,k}  = u_k^3 y_i(t_k),
```

with

```text
deg A <= 8,  deg B <= 12,  deg x_i <= 4,  deg y_i <= 6.
```

The five-fibre system has a useful interpolation reduction.  Use the weighted
`PGL_2` action on the base to set three distinct parameters to `0,1,-1`, and
use the global Weierstrass scaling to set one `u_k=1`.  The five values then
determine every quartic `x_i` uniquely.  Interpolation leaves:

```text
2  unknown base parameters,
4  remaining fibre scalings,
4  free coefficients of A after five values,
8  free coefficients of B after five values,
34 free coefficients in the seventeen sextic y_i after five values,
--
52 unknowns.
```

Equating coefficients in

```text
y_i(t)^2 = x_i(t)^3 + A(t)*x_i(t) + B(t)
```

gives thirteen equations for each of the seventeen labelled sections.  This
is a finite, heavily overdetermined modular system of 221 coefficient
equations in the 52 reduced unknowns, away from interpolation diagonals and
zero scalings.  The first admissible small primes at which all five stored
fibres are nonsingular and all 170 point coordinates are defined are

```text
17, 53, 67, 79, 83, 101, 137, 149, 157, 163.
```

This count only defines the next exact experiment.  It does not show that the
system has a solution, that the public point labels are the desired generic
sections, or that a modular solution lifts to characteristic zero.  A modular
solver should saturate the distinct-parameter, nonzero-scaling, and
interpolation-denominator factors before CRT/LLL reconstruction.

The first-jet elimination now removes all 34 free ordinate coefficients and
all twelve free surface coefficients before solving. Its complete projective
mod-17, mod-53, and mod-67 charts, including both infinity orientations, have
no solution for literal labels and signs. See
[`ICARM_WGXLI_RANK17_FIRST_JET_ELIMINATION.md`](ICARM_WGXLI_RANK17_FIRST_JET_ELIMINATION.md).
This is an exact necessary-condition obstruction in those reduction charts,
not a characteristic-zero nonexistence theorem.

<!-- status-consumer: EC-ICARM-WGXLI-R17-FIRST-JET 11b13e24c5e42a14 -->

The subsequent bounded rebasing audit exhausts relative diagonal signs and
fingerprint-indistinguishable permutations. It then exhausts one common,
anchor-preserving elementary shear under the recorded coefficient and height
bounds. The sole proposed shear `P4 -> P4-P1` is constructed exactly and also
fails the complete mod-17 and mod-53 projective eliminations. See
[`ICARM_WGXLI_RANK17_BOUNDED_REBASING.md`](ICARM_WGXLI_RANK17_BOUNDED_REBASING.md).

<!-- status-consumer: EC-ICARM-WGXLI-R17-BOUNDED-REBASING 6e0c7b116b5b25c3 -->

The bounded fixed-root Mestre census is not repeated or enlarged here.  Its
negative result for curve 356 remains exactly the previously declared census
boundary.

## Reproduction

From the repository root:

```bash
.venv/bin/python \
  elliptic-curves/cas/analyze_icarm_wgxli_rank17_lineages.py

.venv/bin/python \
  elliptic-curves/cas/analyze_icarm_wgxli_rank17_lineages.py \
  --write-artifact \
  artifacts/generated-results/elliptic-curves/icarm_wgxli_rank17_lineage_v1.json
```

The command requires network access to the hash-pinned ICARM sources and
PARI/GP for the numerical height matrices.  The generated artifact has SHA-256

```text
f875f3917486c78089c6ba618daaef7dea07cb36152dd192d0f00c5d425c1c03
```

before any intentional regeneration following this note.

## Public sources

- [ICARM database](https://elliptic-rank.icarm.cloud/database.json), snapshot
  SHA-256 `18699517c2969c8c3a250ae612d5caae9fb23c379fe054ba3c7fdf2ec2a83e50`.
- [ICARM curve 351](https://elliptic-rank.icarm.cloud/curve/351).
- [ICARM curve 356](https://elliptic-rank.icarm.cloud/curve/356).
- [ICARM curve 376](https://elliptic-rank.icarm.cloud/curve/376).
- [ICARM curve 377](https://elliptic-rank.icarm.cloud/curve/377).
- [ICARM curve 385](https://elliptic-rank.icarm.cloud/curve/385).
