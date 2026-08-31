# Bisection specialization at the high-rank controls

<!-- status-consumer: EC-K3-ELKIES-2026-BISECTION-SPECIALIZATION-CONTROLS 04f49e48e1c1dd88 -->

## Exact result

All 39,120 branch quadratics in the complete published-(R17) bisection
atlas were evaluated at

```text
-2/377, -308/251, 2456/135, -9529/5471, 3/8.
```

This is a complete (5\mathbin{\times}39{,}120=195{,}600)-test census, not a
bounded sub-atlas.  For every nonzero rational square, the replay substitutes
both signs of the square root in the stored section

```text
(x,y)=(x0(t)+x1(t)u, y0(t)+y1(t)u),  u^2=q(t),
```

checks both points on the specialized projective fibre and its global minimal
model, and checks that their elliptic sum is the stored trace section.

The finite-quotient calculation uses the deterministic public complements
already certified for the five fibres.  It scans every usable good-reduction
prime through 1,000 and computes the direct sum of
(E(\mathbf F_p)/2E(\mathbf F_p)).  The resulting summary is:

| $t$ | known rank lower bound | known quotient dimension over $R17$ | split bisections | span in known quotient | finite-quotient escapes | exact rank of generated subgroup after adjoining |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `-2/377` | 25 | 8 | 6 | 5 | 0 | 25 |
| `-308/251` | 26 | 9 | 3 | 3 | 0 | 26 |
| `2456/135` | 27 | 10 | 2 | 2 | 0 | 27 |
| `-9529/5471` | 28 | 11 | 1 | 1 | 0 | 28 |
| `3/8` | 21 | 4 | 25 | 4 | 0 | 21 |

Thus the rank-28 fibre realizes the third proposed outcome: the sole split
bisection sees a one-dimensional subspace of the eleven known exceptional
directions.  The complementary ten-dimensional space is invisible to this
height-ten trace atlas at that fibre.  This does not support another shell
merely by itself, but it gives the exact visibility deficit the next shell
would have to address.

At `t=3/8`, the 25 split bisections span all four known directions beyond
(R17).  This directly validates the bisection atlas as a practical detector
on the conductor control, ICARM curve 394, whose existing exact certificate
proves rank at least 21 and

```text
log(N)=166.252098527727201665...
```

## Individual finite-quotient classes

For each fibre, `Q1,...,Qd` denotes the ordered deterministic public
complement, not the full public point list.  Its source indices are:

| (t) | `Q1,...,Qd` are public points |
| --- | --- |
| `-2/377` | `1,4,6,10,11,16,20,24` |
| `-308/251` | `1,2,4,5,6,8,9,10,12` |
| `2456/135` | `1,10,11,13,18,22,24,25,26,27` |
| `-9529/5471` | `1,2,3,4,7,8,9,11,15,19,22` |
| `3/8` | `1,2,3,4` |

The bit strings below are the exact classes in that basis, modulo the generic
seventeen, for the rank-increasing-prime ensemble selected by the complete
prime-through-1,000 scan.

| (t) | bisection orbit | class in `Q1,...,Qd` |
| --- | --- | --- |
| `-2/377` | `orbit-1cb25` | `10000000` |
|  | `orbit-0cff7` | `10000010` |
|  | `orbit-1ea09` | `00000110` |
|  | `orbit-051a1` | `00100000` |
|  | `orbit-0d4ca` | `00000100` |
|  | `orbit-1d5bb` | `00001000` |
| `-308/251` | `orbit-0da89` | `111010100` |
|  | `orbit-12c1b` | `100010111` |
|  | `orbit-1ea54` | `110010111` |
| `2456/135` | `orbit-195a4` | `0100000110` |
|  | `orbit-00edf` | `1010000000` |
| `-9529/5471` | `orbit-15a68` | `01011001010` |
| `3/8` | `orbit-05980` | `1100` |
|  | `orbit-04d17` | `0100` |
|  | `orbit-101f2` | `0101` |
|  | `orbit-07843` | `0011` |
|  | `orbit-090a3` | `1000` |
|  | `orbit-05443` | `1000` |
|  | `orbit-18fd5` | `0100` |
|  | `orbit-1f786` | `0111` |
|  | `orbit-02d31` | `0011` |
|  | `orbit-055ad` | `1000` |
|  | `orbit-0be21` | `0000` |
|  | `orbit-00ca6` | `0100` |
|  | `orbit-0976c` | `1000` |
|  | `orbit-045c2` | `1000` |
|  | `orbit-0fb68` | `1000` |
|  | `orbit-08e3a` | `0010` |
|  | `orbit-196a3` | `1000` |
|  | `orbit-01926` | `1100` |
|  | `orbit-0888a` | `0100` |
|  | `orbit-06faa` | `1000` |
|  | `orbit-10aaa` | `0000` |
|  | `orbit-06f04` | `0000` |
|  | `orbit-0eba4` | `0100` |
|  | `orbit-01e36` | `1100` |
|  | `orbit-126e6` | `0000` |

The generated artifact records every exact value `q_i(t_0)`, its canonical
positive rational square root, both minimal-model points, the generic
correction vector, the finite-prime ensemble, and hashes of the stacked
signature rows. PARI's height matrix and LLL reduction discover a relation
lattice; every retained relation is then proved independently by exact
rational group addition. The coefficient block on the new points has full
rank `6,3,2,1,25`, proving that adjoining all split points leaves the ranks of
the displayed generated subgroups at `25,26,27,28,21`.

## Rank boundary

No split point escapes the certified public finite-quotient span. More
decisively, the exact relation blocks prove that all split points lie in the
rational spans of the displayed public bases. Hence the generated subgroup
ranks after adjoining are exactly 25, 26, 27, 28, and 21. This remains a
subgroup statement: it gives no upper bound for the full Mordell--Weil group,
so no exact curve-rank claim is made here.

## Replay

```sh
.venv/bin/python \
  elliptic-curves/scripts/evaluate_elkies_2026_bisections_at_controls.py \
  --check
```

The pinned certificate is
[`elkies_2026_bisection_specialization_controls_v1.json`](../../artifacts/generated-results/elliptic-curves/elkies_2026_bisection_specialization_controls_v1.json).
