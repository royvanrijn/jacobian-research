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

## Grouped-root certificate

The [grouped-root collision inequality](ABC_OMITTED_COORDINATE_GRAPH.md#collision-bounds-for-square-and-cube-roots-drawn-from-sparse-sets)
adds valid rows to the repaired program. These rows retain the six-coordinate
projection and use the same tail weight \(7\) in the lower bound for the
frozen modulus. Their diagonal term is affine. In the off-diagonal term
each positive part has a nonnegative coefficient. The off-diagonal term is
therefore the maximum of the affine expressions obtained by independently
selecting either zero or the argument of each positive part. At every such
split the diagonal branch and these off-diagonal branches cover all possible
profiles. The grouped residue constraint contributes at most one additional
positive part.

At budget one the retained [certificate delta](abc_grouped_fulloc_delta.json.gz)
replaces precisely **4,834** leaves whose exact original bounds exceed
\(1747/3000\). The remaining leaves keep their original multipliers. The
replacement subtrees use the original geometry bound and the new grouped-root
inequality. Each new leaf has a positive rational combination of the root and
branch rows with coefficient \(1\) on \(D\), nonnegative coefficients on the
other variables, and right-hand side at most \(1747/3000\).
The verified composite certificate has **106,845** leaves and gives the upper bound
\(1747/3000\). Floating-point optimization selected branches and candidate
multipliers during discovery. The verifier uses exact rational arithmetic.

From the repository root, run the [verifier](verify_abc_grouped_fulloc_certificate.py).

```sh
python3 -B research/arithmetic/verify_abc_grouped_fulloc_certificate.py
```

The verifier rebuilds and checks the previous \(7/12\) tree first. It computes
the complete set of paths that require replacement, reconstructs their
ancestor rows, and checks every branch and leaf of each replacement. It also
checks the \(37/66\) bound on the symmetric profile directly. The input tree,
upstream checker and the repair for the previous tree are the files already
identified above.

### Dependence on the radical budget

For a leaf with verified right-hand side \(b\) at budget one, write \(m\) for
its multiplier on the initial row \(L\le1\). At budget \(\lambda\) the
same multipliers give \(D\le b-m(1-\lambda)\). The maximum of the retained
and replacement leaf lines on \(0\le\lambda\le1\) is exactly
\[
 U_{\mathrm{cert}}(\lambda)=\max\left\{
 \frac{33678}{57833},
 \frac{65607}{413150}+\frac{87492}{206575}\lambda
 \right\}.
\]
The verifier checks that both lines occur and that every effective leaf
line is bounded by their maximum at zero, one and the single crossing
\(146664001/146664488\). Affinity gives the inequality on each intervening
interval. In particular the exact maximum at budget one is
\(240591/413150<1747/3000\).

## Exponent inequalities with both tall terms retained

For the height range \(Y/2\le c\le Y\), the larger addend is at least
\(Y/4\). After carrying the signs of the original equation through a
permutation, label the shorter addend \(a\) and the two terms comparable
with \(Y\) by \(b,c\), ordering the latter by
\(\alpha_{T,4}+2\alpha_{T,5}\). Thus of the two order inequalities in
the original FullOC root, retain only
\(\alpha_{c,4}+2\alpha_{c,5}\le\alpha_{b,4}+2\alpha_{b,5}\).
The height argument is independent of which of the two tall labels
is the original largest term.

A selection \(J\) on one letter may contain all higher coordinates
or none of them. Put \(\tau_T=L_T-\sum_{i\le6}\alpha_{T,i}\).
When \(J\) contains the entire higher part, its mass is
\(\tau_T+\sum_{i\in J,\ i\le6}\alpha_{T,i}\). Its upper weighted
height is \(1-\sum_{i\notin J,\ i\le6}i\alpha_{T,i}\). For either
tall letter the latter expression is also a lower bound for the
weighted height in the limiting program. For the shorter letter the
lower bound is \(7\tau_T+\sum_{i\in J,\ i\le6}i\alpha_{T,i}\).
A selection of only early coordinates has its ordinary weighted
height in both bounds. The bounds hold up to the vanishing
decomposition errors and allow the exceptional coefficients in the
frozen moduli.

These choices apply to the geometry estimate by releasing the product
of all selected coordinates. Its upper exponent uses the upper
weighted heights. In a collision estimate the modulus exponents
are lower weighted heights, while root and free-value exponents are
upper heights. A grouped root that includes the higher part can be
taken with exponent one, regardless of the individual degrees in
that part. The resulting root set is sparse at its mass scale.
Releasing the whole higher part of \(A\) leaves the off-diagonal
collision bound unchanged when the two root groups are kept fixed
and lowers the diagonal term. The direct congruence-energy estimate
uses the same lower bound for the modulus of its locked letter.
[Grouped collision and direct residue bounds](ABC_OMITTED_COORDINATE_GRAPH.md#direct-congruence-bound-from-two-residue-energies)


### A Fourier bound on aggregate tail mass

For one letter \(T\), let \(F_T\) be the full grouped Fourier
quantity, the maximum mass of degrees divisible by a common prime.
Write \(e_T=\sum_{i\le6}i\alpha_{T,i}\). For a prime \(q\) let \(m_q\)
be the smallest degree \(i\ge7\) with least prime factor \(q\).
Thus \(m_2=8\), \(m_3=9\), \(m_5=25\) and \(m_q=q\) for \(q\ge7\).
For any fixed \(t\ge7\), one has the affine inequality
\[
 e_T+t\tau_T
 \le 1+\sum_q(t-m_q)_+
   \left(F_T-\sum_{\substack{i\le6\\q\mid i}}\alpha_{T,i}\right).
\]
Subtract the complete weighted height from the left-hand side, discard
negative terms \((t-i)\alpha_{T,i}\) for \(i\ge7\), and assign each
remaining term to its least prime \(q\). For a given \(q\), its weight
is at most \(t-m_q\), and the early multiples of \(q\) already consume
part of its Fourier mass allowance. At \(t=7\) the formula gives the
usual tail height inequality. The additional inequalities can be
imposed using the current projected variables and apply uniformly
to unknown distributions of higher degrees.


### Exact certificates at fixed radical budgets

The additional rows give the following upper certificates for the
tagged program. Each tree starts from the three heights, the projected
pair and Fourier bounds, the budget inequality, and the single order
inequality for the two tall terms.

| Radical budget | Upper bound for \(D\) | Leaves |
| --- | --- | --- |
| \(7/10\) | \(2/5\) | 7,268 |
| \(4/5\) | \(9/20\) | 21,808 |
| \(9/10\) | \(13/25\) | 3,344 |
| \(1\) | \(11/20\) | 280,510 |

These trees start with the aggregate Fourier tail inequalities at
cutoffs \(8,9,10,12,15,18,20,24,32\). Their branches also use the
[interpolated energy bound and the reduced-fraction count](ABC_OMITTED_COORDINATE_GRAPH.md#uniform-energy-bounds-from-interpolation).
The archive retains the earlier \(43/100\) and \(49/100\) trees
for reweighting at other budgets.

At each split, the diagonal and off-diagonal terms, or the geometry
alternatives, give a disjunction of affine inequalities. For a sum
of positive parts with nonnegative coefficients, the verifier lists
every affine term of the maximum. Each leaf verifies a rational
linear combination of its ancestor and root rows, with coefficient
one on \(D\) and nonnegative coefficients on all other variables.
Floating-point optimization proposed the trees and multipliers.
[Certificate archive](abc_height_fulloc_certificate.json.gz) ·
[Exact verifier](verify_abc_height_fulloc_certificate.py)

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py
```

A leaf with bound \(b\) and multiplier \(m\) on the radical budget
row at budget \(\lambda_0\) gives the bound \(b+m(\lambda-\lambda_0)\)
at any fixed \(0<\lambda\le1\). Taking the largest of these bounds in
each tree and the smallest across the available trees gives an
upper exponent function. For example, \(\lambda=3/4\) gives \(1773/4012\).
The verifier option `--at-budget 3/4` evaluates this function exactly.

### Height refinements of the budget-one tree

At budget one the [height delta](abc_height_fulloc_delta.jsonl.gz)
refines the grouped-root tree at exactly **12,272** leaf bounds above
\(29/50\). It keeps the original root ordering. At a refined leaf,
the three possible positions of the shorter addend cover the boxes.
In each position, the other two letters have the tall modulus bounds
for frozen tails. Some replacements close directly with a tighter
rational combination of the existing rows and need no height case.
The composite tree has **368,294** effective leaves and exact upper
bound \(29/50\).

The [delta verifier](verify_abc_height_fulloc_delta.py) checks both
earlier layers, all replacement paths and all new branches. It checks
each rational dual using integer-scaled rows.

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_delta.py
```
