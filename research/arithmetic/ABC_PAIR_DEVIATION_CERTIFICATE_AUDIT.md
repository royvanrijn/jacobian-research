# Coverage gap and exact optimum of the abc modulus-locking program

The 20 November 2026 supplement to Zeng's modulus-locking preprint reports
\(5845/10009\) from an 84,147-leaf rational certificate. A point in one
ancestor region lies above that value and misses every child of the
pair-deviation branch. Completing the branch and replacing two other
leaves gives an 84,177-leaf certificate with maximum bound \(7/12\).
Together with Zeng's feasible point, this determines the optimum of the
FullOC exponent program. [20 November supplement](https://doi.org/10.5281/zenodo.24407593) [Original preprint and reproduction archive](https://doi.org/10.5281/zenodo.24354214)

## The missing branch

Put \(s_i=\alpha_{i,3}+\alpha_{i,4}\) for \(i\in\{a,b,c\}\) and
\(e=1/100\). The supplied verifier's pair-deviation case gives three
branches \(s_i-s_j\ge e\) for \((i,j)=(a,b),(a,c),(b,c)\), and one branch
with all \(|s_i-s_j|\le e\). The reverse inequalities
\(s_j-s_i\ge e\) are absent. The root ordering is on
\(\alpha_{i,4}+2\alpha_{i,5}\), a different linear form.

The tree reaches this node by the zero-based child indices
\[
(0,0,0,1,0,0,1,1,1,0,0,2,2,2,2,2,1,1,1,1,5).
\]
The following point satisfies all 115 root and ancestor rows. All entries
in the table are numerators over \(86149\), with last coordinate
\(D=50309/86149=7187/12307\).

| Letter | \(\alpha_{i,1}\) | \(\alpha_{i,2}\) | \(\alpha_{i,3}\) | \(\alpha_{i,4}\) | \(\alpha_{i,5}\) | \(\alpha_{i,6}\) | \(\alpha_i\) | \(F_i\) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| \(a\) | 2365 | 5523 | 11266 | 5500 | 3388 | 0 | 28042 | 11266 |
| \(b\) | 1573 | 5754 | 11266 | 5720 | 3278 | 0 | 27591 | 11474 |
| \(c\) | 5971 | 3388 | 11277 | 9829 | 51 | 0 | 30516 | 13217 |

Here \(s_a<s_b<s_c\), so the first three children fail. Since
\(s_c-s_a=20/397>1/100\), the fourth child fails too. The excess over
the announced bound is \(268/123180763\). The geometric inequality at
\((I_a,I_b,I_c)=(\{1\},\{1,2,4\},\{1,2,3,5\})\) cuts this witness, with
right-hand side \(7150/12307\). The point diagnoses a failure of coverage at
the selected node. [Geometric inequality, Proposition 2.4](https://doi.org/10.5281/zenodo.24354214)

## Exact upper bound

Append the three reverse inequalities \(s_b-s_a\ge e\), \(s_c-s_a\ge e\),
and \(s_c-s_b\ge e\) after the four original children. Every real triple
\((s_a,s_b,s_c)\) lies in one of the six directed branches or has all three
absolute differences at most \(e\). The three appended subtrees have
\(1,9,9\) leaves. Each of their leaf bounds is at most \(7/12\).

Among the 84,147 original leaves, exactly two exceed \(7/12\). The delta
replaces the leaves bounded by \(1689/2893\) and \(5845/10009\) with
subtrees of nine and four leaves, respectively. All other original leaves
retain their inherited inequalities and multipliers. Every further split
in the delta uses one of the original geometric or collision inequalities.
The 19 and 20 November checkers use identical root, geometry and collision
row formulas. Both supplied upper certificates use the 11 coprime pairs
of exponents at most six, a subset of the full 15-pair program, and every
third-letter modulus set used in the certificates omits index 1.

At each leaf the positive rational multipliers give coefficient \(1\)
on \(D\), nonnegative coefficients on all other variables, and a rational
right-hand side at most \(7/12\). The repaired tree has 5,330 geometric
nodes, 9,855 collision nodes, one pair-deviation node and 84,177 leaves.
Its maximum leaf bound is exactly \(7/12\). Zeng's Proposition 8.3 supplies
the matching feasible point for the full 15-pair family, including all
32 third-letter modulus sets. The result is the exact maximum of the
finite limiting program. [FullOC inequalities and lower obstruction, §§2–5 and Proposition 8.3](https://doi.org/10.5281/zenodo.24354214)

## Consequence for the exceptional-set count

Let \(N_\lambda(X)\) count ordered primitive positive triples
\(a+b=c\le X\) with \(\operatorname{rad}(abc)<c^\lambda\).
The compactness argument in Zeng's §7 applies with \(U=7/12\).
If for some \(\varepsilon>0\) the bound at \(\lambda=1\) failed, the
dyadic and anatomic decompositions would give boxes with exponents
\(D\ge D_0>7/12\) and errors tending to zero. The first six coordinates
and the three totals lie in a fixed compact set. Their limit satisfies
the original closed root, geometry and collision inequalities, and the
corrected tree covers the limit. This contradicts \(D\le7/12\).
Monotonicity in \(\lambda\) yields
\[
N_\lambda(X)\ll_\varepsilon X^{7/12+\varepsilon}
\qquad (0<\lambda\le1).
\]
The exponent-program feasible point establishes no corresponding lower
bound on the count of integer triples. [Counting reduction, §7](https://doi.org/10.5281/zenodo.24354214)

## Replay

The [certificate delta](abc_pair_deviation_7_12_delta.json) contains
the new subtrees, rational multipliers and two exact replacement paths.
From the repository root run the [standard-library verifier](verify_abc_pair_deviation_certificate.py).

```sh
python3 -B research/arithmetic/verify_abc_pair_deviation_certificate.py
```

The verifier checks the SHA-256 of the two upstream archive members,
replays the rational witness, applies the delta and checks all 84,177
leaves. Numerical optimization was used only to discover the replacement
subtrees and candidate multipliers.

| Original member | SHA-256 |
| --- | --- |
| [Compressed 84,147-leaf tree](sources/zeng-abc-modulus-locking-24407593/cert_fulloc_5840_improved_73_125.json.gz) | `a94bbf1c0a8bd1267f9be28cd862d980e2fe1eab7b74f2b6760d76d342ae8069` |
| [Row and tree verifier](sources/zeng-abc-modulus-locking-24407593/check_extended_pairdev.py) | `a1ac4cd83a296c55e2c9a245e845d66dc7935f8624febdd25a61a81104b31437` |
