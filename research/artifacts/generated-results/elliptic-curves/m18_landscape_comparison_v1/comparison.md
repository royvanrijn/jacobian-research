# Frozen M18 landscape comparison

All CVPs use the retained rounded integer Gram. Norms are divided by the median of the first 17 diagonal entries. Outcomes are certified lower bounds after bounded V3; 18 is not an upper bound. Quartic sizes cover all 64 refined parities, irrespective of which charts ran.

| M18 state / parameter | Outcome | CVP min / median / max | ≤2 masks | Within 1.25×min | Anchors | Quartic bits median | Tail max/min |
|---|---:|---:|---:|---:|---:|---:|---:|
| 302 recovered-strict-02 | 18→31 | 1.443 / 2.369 / 4.438 | 13 | 5 | 5 | 594 | 3.075 |
| 302 recovered-strict-03 | 18→31 | 1.443 / 2.369 / 7.720 | 13 | 5 | 5 | 594 | 5.349 |
| 302 recovered-strict-01 | pending | 1.443 / 2.284 / 2.959 | 13 | 5 | 5 | 594 | 2.051 |
| 302 residual-strict-03 | 18→31 | 1.443 / 2.369 / 7.835 | 13 | 5 | 5 | 594 | 5.429 |
| 1926/2699 | 18→21 | 1.724 / 2.156 / 2.702 | 16 | 31 | 23 | 415 | 1.568 |
| 2953/1671 | 18→18 | 1.657 / 2.122 / 2.546 | 22 | 27 | 21 | 460 | 1.537 |
| 5193/35630 | 18→18 | 1.410 / 2.018 / 2.450 | 30 | 6 | 5 | 483 | 1.737 |
| 2980/1967 | 18→18 | 1.653 / 2.229 / 2.699 | 13 | 17 | 15 | 459 | 1.632 |
| 1117/2193 | 18→18 | 1.716 / 2.193 / 2.566 | 15 | 27 | 22 | 429 | 1.496 |
| 1004777772/6898189895 | 18→18 | 1.746 / 2.196 / 2.637 | 16 | 31 | 25 | 855 | 1.511 |
| 1609736129/11039275010 | 18→18 | 1.614 / 2.223 / 2.551 | 19 | 19 | 16 | 864 | 1.580 |
| 257306782577/1766023495870 | 18→18 | 1.893 / 2.231 / 2.689 | 9 | 40 | 27 | 1009 | 1.421 |
| 257474118081/1766197190110 | 18→18 | 1.516 / 2.191 / 2.638 | 17 | 1 | 1 | 992 | 1.740 |
| 2073030269/14240022390 | 18→18 | 1.692 / 2.208 / 2.697 | 14 | 22 | 19 | 871 | 1.593 |
| 7159030427/49068265370 | 18→18 | 1.810 / 2.197 / 2.604 | 16 | 34 | 23 | 906 | 1.439 |

## Strict retrospective separations

| Feature | Gain states range | Bounded null range | Gap / total observed range |
|---|---:|---:|---:|
| shell10_cvp_max | 2.7022–7.8348 | 2.4501–2.6966 | 0.001 |
| cvp_max | 2.7022–7.8348 | 2.4501–2.6986 | 0.001 |
| extension1_cvp_max | 2.7022–7.8348 | 2.4501–2.6986 | 0.001 |

This scan is exploratory and uses many correlated features. Curve 302 repeats one fibre and its seeds were obtained retrospectively. No separating threshold has been prospectively validated. Missing panel outcomes are excluded from the separation scan.
