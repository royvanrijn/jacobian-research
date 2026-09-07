# Reverse-engineering the rank-28 exceptional quotient (2026-09-02)

## Outcome

The eleven public directions at `t=-9529/5471` now have a reproducible
construction fingerprint relative to the specialized generic R17 subgroup
`M_t`.  The result is negative for the most obvious small operations but gives
one exact target object for a geometry-first continuation.

The cubic two-descent Kummer map

```text
Q |-> x(Q)-theta
```

exposes all eleven directions.  Each class has an exact norm square and an
explicit intersection of quadrics with verified rational witness
`[1:0:0:1]`.  Thus this is the first tested algebraic object that sees `11/11`,
where the complete rational-bisection atlas sees only the collective class

```text
Q2-Q4+Q5-Q8+Q10.
```

This is a representation mechanism, not yet a construction mechanism: the
Kummer class is computed from a point after that point is known.  It therefore
does not by itself turn an unknown rank-32 direction into a rational point.

Small division and isogeny operations do not split off a large common class.
For every `Q_i`:

- the `[2]` preimage polynomial is an irreducible quartic with Galois group
  `S4`;
- the `[3]` preimage polynomial is an irreducible nonic with PARI group
  `E(9):2S_4`, order `432` and transitive id `26`;
- within each level, all eleven polynomial discriminants have the same
  squareclass;
- the complete rational isogeny graph has one vertex and degree matrix
  `Mat(1)`.

The division fields therefore have the generic large first-level behavior for
all eleven directions.  Their exact factorization patterns at 32 common
unramified primes distinguish individual fields, but no common `[2]`, `[3]`,
or rational-isogeny operation explains seven to ten directions.

## Smallest relative representatives

Write the canonical-height Gram in the public rank-28 basis and project away
the real span of `M_t`.  For each public complement point the calculation
completely enumerates the first coset shell on integral Gram matrices rounded
at scales `10^4`, `10^5`, and `10^6`.  All three runs return the same result:
the public `Q_i` itself is already the shortest representative found in
`Q_i+M_t`; every generic correction is zero.

The first nonzero shell of the specialized generic subgroup has height
`47.878552386064374...`.  The public directions split as follows.  “First
shell” means no taller than this specialized generic first shell.

| direction | `x` height bits | raw height | quotient defect | shortest `Q_i+M_t` | first shell | parent base bits | parent height | exploratory cluster |
|---|---:|---:|---:|---:|:---:|---:|---:|---:|
| Q1 | 91 | 47.265917 | 16.530291 | 47.265917 | yes | 105 | 431.854 | 1 |
| Q2 | 91 | 47.341600 | 25.033733 | 47.341600 | yes | 79 | 351.557 | 1 |
| Q3 | 91 | 47.480937 | 27.069127 | 47.480937 | yes | 95 | 442.216 | 3 |
| Q4 | 91 | 47.561506 | 18.967883 | 47.561506 | yes | 113 | 472.844 | 2 |
| Q5 | 92 | 48.009635 | 20.092834 | 48.009635 | no | 116 | 558.684 | 2 |
| Q6 | 91 | 48.058291 | 24.742883 | 48.058291 | no | 85 | 387.604 | 3 |
| Q7 | 91 | 48.067204 | 20.620290 | 48.067204 | no | 77 | 352.165 | 3 |
| Q8 | 91 | 48.097760 | 22.834605 | 48.097760 | no | 101 | 457.818 | 3 |
| Q9 | 91 | 48.311435 | 21.128955 | 48.311435 | no | 91 | 421.328 | 3 |
| Q10 | 95 | 48.625264 | 21.117864 | 48.625264 | no | 66 | 319.355 | 5 |
| Q11 | 92 | 48.717476 | 18.869402 | 48.717476 | no | 78 | 355.455 | 4 |

For Q1--Q4 the smallest positive multiplier entering that shell is `n=1`.
For Q5--Q11 no positive multiplier enters it: `n=1` is already above the
threshold, while for every `n>=2` the continuous quotient lower bound
`n^2 ||Q_i+M_t||_quot^2` exceeds the first-shell height.  This conclusion is
relative to the stored 90-digit canonical-height form; it is numerically
stable, not an interval-certified height theorem.

A point on one elliptic fibre is codimension two on the K3, so it has no
intrinsic K3 divisor class or geometric intersection profile.  The certificate
records this as `NOT_INTRINSIC_FOR_A_FIBRE_POINT`.  It stores instead the
quotient-height pairings with all other `Q_j`, the exact bad-place Kummer code,
and the exact parent-fibration transport.  A genuine divisor class enters only
after choosing a multisection through the point.

## Exploratory clustering

An equal-weight consensus distance combines four views:

1. coordinate, height, Kummer, parent, and division-polynomial sizes;
2. the exact 53-bit bad-place Kummer code;
3. `[2]` and `[3]` Frobenius factor patterns at 32 common unramified primes;
4. angles in the quotient-height Gram.

Average-linkage clustering, with the number of clusters chosen by mean
silhouette among `k=2,3,4,5`, selects

```text
{Q1,Q2}, {Q4,Q5}, {Q3,Q6,Q7,Q8,Q9}, {Q11}, {Q10}.
```

The largest group has size five and the selected silhouette is only `0.2814`.
This is evidence against a visible seven-to-ten-direction small-operation
class in the chosen features, not evidence that no such geometric operation
exists.  The eleven points are a historically selected basis, feature weights
are choices, and singleton clusters must not be overinterpreted.

## Symbolic generic-R17 replay

On the short model `y^2=x^3+A*x+B`, the first two preimage objects are the
universal multiplication covers

```text
F_n(z;X) = numerator(x([n]R))-X*denominator(x([n]R)).
```

The exact quartic is

```text
z^4 - 4X*z^3 - 2A*z^2 - (4AX+8B)z + A^2-4BX.
```

The artifact also stores the exact degree-nine formula.  Both are irreducible
over `QQ(A,B,X)`.  After substituting the published R17 polynomials `A(t)` and
`B(t)`, they remain irreducible over `QQ(t,X)`; hashes pin the 3,132-byte and
32,399-byte expansions.  Specializing `t=-9529/5471` and `X=x(Q_i)` recovers
the twenty-two displayed defining polynomials.

This symbolic replay explains the negative result: `[2]` and `[3]` preimages
are universal covers in an independent point coordinate `X`.  They do not
produce a rational function `X(t)` and hence do not reproduce an exceptional
section on generic R17.

## Search consequence

The next geometry-first experiment should target the exact Kummer classes,
not another parameter score.  Concretely, search either

- higher-genus degree-two curves on the published R17 surface; or
- targeted higher-degree multisections,

subject to the constraint that their specialization at the rank-28 fibre maps
to one of the ten bisection-invisible classes `x(Q_i)-theta`.  The Kummer class
is then an exact terminal constraint for symbolic elimination.  A mechanism is
promoted only if it deforms over `QQ(t)`, specializes to several independent
target classes, and gives a reproducible rational curve or multisection.  The
present calculation does not authorize a Nagao sieve or a point-height search.

This viewpoint is compatible with the multisection/rank-jump framework of
Garbagnati--Salgado and with the use of explicit covering classes in descent;
see [Rank jumps and Multisections of elliptic fibrations on K3
surfaces](https://arxiv.org/abs/2505.15159) and Bruin--Dahmen,
[Visualizing elements of Sha[3] in genus 2
Jacobians](https://arxiv.org/abs/1001.5302).  Preimage-field factorization is
used here as a finite construction fingerprint in the spirit of the
preimage-tree literature; it is not an arboreal maximality theorem.

## Reproduction

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/analyze_elkies_2026_rank28_construction_fingerprints.py

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elliptic-curves/cas/analyze_elkies_2026_rank28_construction_fingerprints.py \
  --check
```

The generated certificate is
[`../artifacts/generated-results/elliptic-curves/elkies_2026_rank28_construction_fingerprints_v1.json`](../artifacts/generated-results/elliptic-curves/elkies_2026_rank28_construction_fingerprints_v1.json).
Its SHA-256 is
`a7261c221b174a5cca3fffdbb1dbbbdb040dea869a42617e0b393dad9c7edd64`.
