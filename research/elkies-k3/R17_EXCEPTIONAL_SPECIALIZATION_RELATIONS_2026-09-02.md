# Exceptional specialization relations on the rank-25--28 R17 fibres

## Result

The 38 displayed exceptional quotient basis vectors have now been compared
with the specialized generic R17 subgroup, the complete rational-bisection
atlas, the tested degree-three/four equations, and the exact q12/orbit5867
map back to the 4A1/MW13 fibration.

The main conclusion is negative but sharper than the previous visibility
rank count.

- Only four individual directions on the rank-25 fibre are direct
  specializations of a rational bisection: `Q1,Q3,Q5,Q6`. Their generic field
  degree is exactly two.
- Rank-25 `Q7` and rank-26 `Q2` are exact sums of points from two independent
  quadratic bisection fields. They deform over biquadratic covers of degree
  four. This is the least known degree, not a proved absolute minimum, because
  higher-genus bisections and the complete degree-three universe have not been
  constructed.
- None of the other 32 unit directions is in the compositum span of the
  rational bisections that split at its fibre. In particular, none of the
  eleven rank-28 basis vectors is individually a known bisection
  specialization. The familiar rank-28 class

  ```text
  Q2-Q4+Q5-Q8+Q10
  ```

  is one collective relation. Calling `Q2` “visible” only records the chosen
  row-reduction pivot; it does not make `Q2` a branch of that bisection.
- The 12 rational-bisection equations that split across the four controls have
  12 distinct orbit labels. No equation repeats at two controls. The split
  counts `6,3,2,1` and captured ranks `5,3,2,1` therefore come from disjoint,
  rapidly thinning square-value events.
- The ten canonical rank-28 unexplained directions remain unexplained by the
  complete rational-bisection atlas and every tested trisection/quadrisection
  equation. Their quotient heights and their images under the second
  fibration show no low-complexity block or correlation.

The exact and numerical record is
[`../artifacts/generated-results/elliptic-curves/elkies_2026_exceptional_specialization_relations_v1.json`](../artifacts/generated-results/elliptic-curves/elkies_2026_exceptional_specialization_relations_v1.json).
Its SHA-256 is
`37566ee1c19a6afcd0d4d3f8ec243ba8388136588cc34fdc369c178f90c825b4`.

## What “divisor class” can mean here

The published parameter `t` is the base coordinate of one fixed elliptic K3
surface. It is not a parameter in a family of K3 surfaces. A rational point on
the elliptic fibre `E_t` is codimension two on the surface and has no intrinsic
class in `NS(X)`. It acquires a divisor class only after a multisection through
it has been chosen.

Consequently an isolated Noether--Lefschetz interpretation is not available:
the K3 surface and its geometric Picard rank remain fixed. Extra rational
points in `E_t(QQ)` do not create new K3 divisor classes. The relevant geometry
is instead splitting of existing multisections or arithmetic points confined
to individual fibres. This is consistent with the distinction in Theorem E of
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md).

For every displayed point `Q`, the calculation uses the canonical arithmetic
substitute

```text
delta(Q) = <Q,Q> - <Q,G> <G,G>^(-1) <G,Q>,
```

where `G=(P1,...,P17)` is the exact specialized generic subgroup and brackets
are Neron--Tate pairings on the number-field fibre. This is the squared norm
of `Q` in the real height quotient by `G`. The full Schur-complement Gram is
stored, not only its diagonal. These values are 90-digit PARI computations;
the integer embeddings of `G`, their primitivity, the public-point relations,
and the cover equations are exact.

## Per-direction results

“Removed” is the fraction of the raw canonical height removed by orthogonal
projection to the specialized generic subgroup. “Parent bits” and “parent
height” are measured after the exact q12/orbit5867 transport to the normalized
4A1/MW13 parent. “Unknown” means no absolute minimum is claimed.

| `t` | direction | public point | quotient defect | removed | parent bits | parent height | least known deformation |
|---|---:|---:|---:|---:|---:|---:|---|
| `-2/377` | Q1 | 1 | 10.677562 | 65.8% | 46 | 223.657 | degree 2, `orbit-1cb25` |
| `-2/377` | Q2 | 4 | 10.192867 | 69.0% | 84 | 316.907 | unknown |
| `-2/377` | Q3 | 6 | 8.344908 | 75.4% | 82 | 345.624 | degree 2, `orbit-051a1` |
| `-2/377` | Q4 | 10 | 11.758931 | 66.2% | 59 | 240.246 | unknown |
| `-2/377` | Q5 | 11 | 9.080505 | 74.0% | 91 | 355.988 | degree 2, `orbit-1d5bb` |
| `-2/377` | Q6 | 16 | 9.153663 | 74.5% | 53 | 193.668 | degree 2, `orbit-0d4ca` |
| `-2/377` | Q7 | 20 | 8.875768 | 77.1% | 104 | 437.741 | degree 4, `orbit-1cb25` and `orbit-0cff7` |
| `-2/377` | Q8 | 24 | 11.226782 | 72.6% | 72 | 327.179 | unknown |
| `-308/251` | Q1 | 1 | 15.206342 | 54.4% | 99 | 404.745 | unknown |
| `-308/251` | Q2 | 2 | 15.633060 | 53.7% | 73 | 304.709 | degree 4, `orbit-12c1b` and `orbit-1ea54` |
| `-308/251` | Q3 | 4 | 25.344305 | 26.2% | 78 | 315.273 | unknown |
| `-308/251` | Q4 | 5 | 14.463170 | 58.0% | 55 | 227.376 | unknown |
| `-308/251` | Q5 | 6 | 14.435902 | 58.2% | 66 | 255.642 | unknown |
| `-308/251` | Q6 | 8 | 12.678887 | 63.4% | 72 | 263.457 | unknown |
| `-308/251` | Q7 | 9 | 17.421453 | 49.8% | 91 | 367.449 | unknown |
| `-308/251` | Q8 | 10 | 13.258388 | 61.8% | 56 | 244.617 | unknown |
| `-308/251` | Q9 | 12 | 11.838899 | 65.9% | 68 | 274.597 | unknown |
| `2456/135` | Q1 | 1 | 13.622455 | 66.2% | 107 | 492.912 | unknown |
| `2456/135` | Q2 | 10 | 21.639900 | 49.6% | 86 | 399.037 | unknown |
| `2456/135` | Q3 | 11 | 23.782485 | 44.8% | 60 | 214.686 | unknown |
| `2456/135` | Q4 | 13 | 17.461561 | 60.1% | 95 | 434.254 | unknown |
| `2456/135` | Q5 | 18 | 18.887720 | 57.9% | 71 | 295.954 | unknown |
| `2456/135` | Q6 | 22 | 17.341111 | 62.4% | 133 | 631.967 | unknown |
| `2456/135` | Q7 | 24 | 21.905144 | 53.2% | 55 | 239.854 | unknown |
| `2456/135` | Q8 | 25 | 24.430377 | 47.8% | 77 | 313.402 | unknown |
| `2456/135` | Q9 | 26 | 20.739613 | 56.2% | 68 | 295.880 | unknown |
| `2456/135` | Q10 | 27 | 17.615988 | 64.9% | 82 | 365.217 | unknown |
| `-9529/5471` | Q1 | 1 | 16.530291 | 65.0% | 105 | 431.854 | unknown |
| `-9529/5471` | Q2 | 2 | 25.033733 | 47.1% | 79 | 351.557 | unknown |
| `-9529/5471` | Q3 | 3 | 27.069127 | 43.0% | 95 | 442.216 | unknown |
| `-9529/5471` | Q4 | 4 | 18.967883 | 60.1% | 113 | 472.844 | unknown |
| `-9529/5471` | Q5 | 7 | 20.092834 | 58.1% | 116 | 558.684 | unknown |
| `-9529/5471` | Q6 | 8 | 24.742883 | 48.5% | 85 | 387.604 | unknown |
| `-9529/5471` | Q7 | 9 | 20.620290 | 57.1% | 77 | 352.165 | unknown |
| `-9529/5471` | Q8 | 11 | 22.834605 | 52.5% | 101 | 457.818 | unknown |
| `-9529/5471` | Q9 | 15 | 21.128955 | 56.3% | 91 | 421.328 | unknown |
| `-9529/5471` | Q10 | 19 | 21.117864 | 56.6% | 66 | 319.355 | unknown |
| `-9529/5471` | Q11 | 22 | 18.869402 | 61.3% | 78 | 355.455 | unknown |

All 38 generic embeddings replay with selected exceptional complement index
one in the displayed public lattice. Thus the table is not using an
unsaturated coordinate projection inside the submitted subgroup.

## Exact algebraic traces and cover degrees

If `P_i` is one branch of the bisection with trace `tau_i`, then

```text
P_i + sigma_i(P_i) = tau_i.
```

For `R=sum c_i P_i` over the compositum of `k` independent quadratic fields,

```text
Trace(R) = 2^(k-1) sum c_i tau_i.
```

The six unit directions currently reached this way are:

| fibre | direction | branch combination | degree | trace in the published R17 basis |
|---|---:|---|---:|---|
| `-2/377` | Q1 | `orbit-1cb25` | 2 | `[-1,-1,0,0,0,-1,-1,0,0,0,0,0,0,0,0,1,0]` |
| `-2/377` | Q3 | `-orbit-051a1` | 2 | `[1,1,0,1,0,1,0,0,0,-1,-1,-1,0,0,1,1,0]` |
| `-2/377` | Q5 | `-orbit-1d5bb` | 2 | `[1,1,0,1,0,1,-1,0,1,-1,-1,0,0,-1,1,0,-1]` |
| `-2/377` | Q6 | `orbit-0d4ca` | 2 | `[0,0,0,-1,-1,1,1,1,1,0,0,0,-1,1,-1,-1,0]` |
| `-2/377` | Q7 | `orbit-1cb25 - orbit-0cff7` | 4 | `[-2,-2,2,0,-2,-2,-2,0,0,0,2,0,2,0,-2,2,0]` |
| `-308/251` | Q2 | `-orbit-12c1b + orbit-1ea54` | 4 | `[0,2,4,-4,4,-4,2,0,0,0,2,2,0,-4,-4,0,-4]` |

The quadratic classes are represented by distinct irreducible quadratics, so
the two degree-four rows genuinely require the displayed biquadratic
composita inside this mechanism. A different genus or cover construction may
still lower their absolute degree; that has not been excluded.

## Comparison of the specialization equations

Each rational bisection is defined by an irreducible quadratic

```text
u^2 = q_0+q_1 t+q_2 t^2.
```

The following primitive triples compare all equations that split on the four
controls. The full artifact retains the rational normalization, exact square
value, square root, generic trace, specialized point relation, and arithmetic
quotient height.

| fibre | cover | primitive `(q0,q1,q2)` | exceptional quotient vector |
|---|---|---|---|
| `-2/377` | `orbit-1cb25` | `(388152818881,1527533196332,1295099052676)` | `(1,0,0,0,0,0,0,0)` |
| `-2/377` | `orbit-0cff7` | `(144546477201,54191016960,7705261840)` | `(1,0,0,0,0,0,-1,0)` |
| `-2/377` | `orbit-1ea09` | `(1582589970256,1900664722216,514992145369)` | `(0,0,0,0,0,1,-1,0)` |
| `-2/377` | `orbit-051a1` | `(6396254569,2363218858,1578980689)` | `(0,0,-1,0,0,0,0,0)` |
| `-2/377` | `orbit-0d4ca` | `(11718940411984,7410964194920,1177213677145)` | `(0,0,0,0,0,1,0,0)` |
| `-2/377` | `orbit-1d5bb` | `(23743585636,12143358604,9658127089)` | `(0,0,0,0,-1,0,0,0)` |
| `-308/251` | `orbit-0da89` | `(17461276021225,17777564958608,3732461513872)` | `(1,1,1,0,1,2,-1,-2,0)` |
| `-308/251` | `orbit-12c1b` | `(-1591877453401,1064065327200,240228538560)` | `(1,0,0,0,1,2,-1,-1,-1)` |
| `-308/251` | `orbit-1ea54` | `(6320290432225,-2708966844430,471529050001)` | `(1,1,0,0,1,2,-1,-1,-1)` |
| `2456/135` | `orbit-195a4` | `(2094878738995344,1919852662499064,551005442820121)` | `(0,-1,0,0,0,0,0,1,1,0)` |
| `2456/135` | `orbit-00edf` | `(192763348439209,220152011867324,48346219485316)` | `(1,0,1,0,0,0,0,0,0,0)` |
| `-9529/5471` | `orbit-15a68` | `(45610752625,93849805300,39065678404)` | `(0,1,0,-1,1,0,0,-1,0,1,0)` |

There is no repeated cover label and no repeated quotient pattern across the
four controls. Primitive coefficient sizes range from 33 to 51 bits and the
specialized square roots from 47 to 75 projective bits, so the rank-28 hit is
not singled out by a larger equation. What changes is incidence: the same
complete atlas sees fewer independent square-value conditions as the public
exceptional rank rises.

## The ten rank-28 directions

For the canonical packet

```text
Q1,Q3,Q4,Q5,Q6,Q7,Q8,Q9,Q10,Q11,
```

the quotient-height diagonal ranges from `16.530291` to `27.069127`, with
median `20.869077`. The ten-dimensional Schur complement has eigenvalues

```text
1.625669, 2.714533, 6.437600, 11.789258, 14.108362,
19.538871, 22.916971, 32.071495, 35.802110, 64.969266.
```

Its largest absolute pair correlation is `0.4268`; in the submitted coordinate
basis it has no obvious repeated or nearly diagonal character block comparable
to the multiquadratic orthogonal blocks in Theorem F4. Under the second
fibration, the ten points land on ten
different parent fibres, with normalized base heights of 66--116 bits and
canonical heights `319.355--558.684`. Across these ten points the quotient
defect has Pearson correlation `-0.027` with parent height and `-0.185` with
parent-base bit height. The two parent complexity measures correlate `0.946`
with each other, but neither traces the R17 quotient defect.

Thus the current evidence rejects two proposed explanations and leaves one
open:

1. **Noether--Lefschetz:** not applicable to specialization in the base of a
   fixed K3 fibration.
2. **Known extension-defined low-degree divisors:** absent for the ten-class
   packet in the complete rational-bisection universe and the tested
   degree-three/four equations. The known bisection divisor is defined over
   `QQ`; only its two branches become sections over a quadratic function
   field.
3. **High-degree or higher-genus multisections:** still possible, but not
   positively indicated by the height or parent-fibration data.

The best current description is therefore **isolated arithmetic fibre points
relative to the tested divisor systems**, not a newly identified repeatable
ten-direction mechanism. This does not conflict with the general non-thin
rank-jump theorem for doubly elliptic K3 surfaces: that theorem supplies a
global abundance mechanism, not a common low-degree divisor for these ten
submitted points. See Pasten--Salgado,
[*Non-thin rank jumps for double elliptic K3 surfaces*](https://doi.org/10.1007/s00229-024-01554-2),
and Garbagnati--Salgado,
[*Rank jumps and multisections of elliptic fibrations on K3 surfaces*](https://arxiv.org/abs/2505.15159).

## Reproduction and boundary

```bash
python3 \
  elliptic-curves/cas/analyze_elkies_2026_exceptional_specialization_relations.py

python3 \
  elliptic-curves/cas/analyze_elkies_2026_exceptional_specialization_relations.py \
  --check
```

This report does not prove that the public point sets are full Mordell--Weil
bases, give exact ranks, exclude higher-genus bisections, exhaust degree three
or four, or produce divisor classes for isolated fibre points. The norm-26
trisection shell is complete only within its declared trace shell; the
norm-20 trisections and norm-34 quadrisections are deterministic samples. The
numerical height quotient is reproducible at 90 decimal digits but is not an
exact algebraic lattice certificate.
