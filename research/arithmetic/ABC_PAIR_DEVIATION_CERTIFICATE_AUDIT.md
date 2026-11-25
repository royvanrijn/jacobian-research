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
| \(9/10\) | \(1/2\) | 39,968 |
| \(1\) | \(107/200\) | 1,822,745 |

These trees start with the aggregate Fourier tail inequalities at
cutoffs \(8,9,10,12,15,18,20,24,32\). Their branches also use the
[interpolated energy bound and the reduced-fraction count](ABC_OMITTED_COORDINATE_GRAPH.md#uniform-energy-bounds-from-interpolation).
The two newest budget-\(9/10\) trees start with the eleven unbranched
[residual-radical inequalities](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#residual-radical-bounds-for-hyperbolic-prime-power-signatures).
These rows retain divisibility information in the six exact valuation
classes. Schema 5 also adds the four [Euclidean signature inequalities](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#residual-radical-bounds-for-two-euclidean-signatures). Its available branch bounds include [separated residue collisions](ABC_OMITTED_COORDINATE_GRAPH.md#separating-repeated-residues-from-the-diagonal) and [the fixed-target cubic–quartic surface](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#a-fixed-target-surface-with-a-cubic-factor-and-a-quartic-power).
Schema 6 adds the [square–cube residual radical rows](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#averaging-the-squarecube-count-over-residual-radicals). Schema 7 also adds the [square–fifth rows](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#square-and-odd-prime-powers). The [repeated-power lattice inequalities](ABC_OMITTED_COORDINATE_GRAPH.md#collinear-congruence-lattices-for-repeated-power-residues) apply to the selected residue branches.
The archive retains the earlier \(43/100\), \(49/100\) and \(13/25\)
trees for reweighting at other budgets.

At each split, the diagonal and off-diagonal terms, or the geometry
alternatives, give a disjunction of affine inequalities. For a sum
of positive parts with nonnegative coefficients, the verifier lists
every affine term of the maximum. Each leaf verifies a rational
linear combination of its ancestor and root rows, with coefficient
one on \(D\) and nonnegative coefficients on all other variables.
Floating-point optimization proposed the trees and multipliers.
The verifier loads the earlier trees from the certificate archive and
concatenates the decompressed unit-budget fragments in numbered order.
[Certificate archive](abc_height_fulloc_certificate.json.gz) ·
[Unit-budget part 1](abc_height_fulloc_unit_certificate.part1.xz) ·
[Part 2](abc_height_fulloc_unit_certificate.part2.xz) ·
[Part 3](abc_height_fulloc_unit_certificate.part3.xz) ·
[Improved unit-budget part 1](abc_height_fulloc_improved_unit_certificate.part1.xz) ·
[Improved unit-budget part 2](abc_height_fulloc_improved_unit_certificate.part2.xz) ·
[Improved unit-budget part 3](abc_height_fulloc_improved_unit_certificate.part3.xz) ·
[Improved unit-budget part 4](abc_height_fulloc_improved_unit_certificate.part4.xz) ·
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

## Exact reweighting at two intermediate budgets

Only root row zero depends on the radical budget. If a leaf was checked
at budget \(\lambda_0\), with exact bound \(b\) and multiplier \(y_0\) for
that row, the same nonnegative dual gives
\[
D\le b+(\lambda-\lambda_0)y_0.
\]
The branch partition and every other coefficient stay valid. Taking the
maximum over all leaves gives the following exact bounds.

| Target \(\lambda\) | Source budget | Source bound | Leaves | Maximum reweighted bound |
| --- | --- | --- | --- | --- |
| \(3/5\) | \(7/10\) | \(2/5\) | 7,268 | \(1959/5015\) |
| \(3/4\) | \(4/5\) | \(9/20\) | 21,808 | \(1773/4012\) |

From the repository root replay the source rows and the reweighted
right-hand sides.

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 7/10 --at-budget 3/5
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 4/5 --at-budget 3/4
```

Monotonicity extends each bound \(N_\lambda(X)\ll_\varepsilon
X^{D+\varepsilon}\) to all positive budgets up to the indicated target.

## A sparse improvement at unit budget

The checked unit-budget tree with upper bound \(107/200\) has 1,822,745
leaves. Replace only the 4,586 leaves whose saved right-hand sides
exceed the new target. The compact [certificate delta](abc_unit_strict_delta.json.xz)
records their paths and replacement trees with exact rational multipliers.
The refined trees use fixed Hölder moments, the standard Hölder bound,
primitive ratio locks and short-addend ratio locks.

The [sparse verifier](verify_abc_unit_strict_delta.py) checks the SHA-256
of the decompressed baseline bytes, enumerates every baseline leaf and
checks that the replacements cover all leaves above the target. It
reconstructs the exact ancestor inequalities and checks every replacement
tree, including every affine branch. The prior exact replay of the pinned
baseline certifies the unaffected leaves. The resulting tree has 1,830,732
leaves. Their maximum bounds are
\[
 U_{\mathrm{replacement}}=\frac{7133338393}{13333343040},\qquad
 U_{\mathrm{untouched}}=\frac{5350003447}{10000006630}.
\]
The unit-budget bound is their maximum \(7133338393/13333343040\). This improves the
baseline bound by \(667/66666715200\).

From the repository root run the sparse verifier.

```sh
python3 -B research/arithmetic/verify_abc_unit_strict_delta.py
```

The optional flag `--full` replays every leaf in the reconstructed
tree with the standard verifier. The raw baseline remains available as
the four pinned archive fragments above.

## Sparse retargeting above budget \(9/10\)

The \(9/10\) source tree with 39,968 leaves has exact bound \(1/2\).
Increasing the budget changes only the first root row's right-hand side.
For a leaf with right-hand side \(b\) and multiplier \(y_0\) for that
row, its retargeted right-hand side is \(b+(\lambda-9/10)y_0\).
At each target budget replace every leaf that exceeds the new limit.
The preserved and replacement subtrees combine to give the following
fully replayed certificates.

| Budget | Exact exponent | Replaced leaves | Replayed leaves | Delta |
| --- | --- | --- | --- | --- |
| \(91/100\) | \(51/100\) | 6 | 39,969 | [abc_budget_910_delta.json.xz](abc_budget_910_delta.json.xz) |
| \(23/25\) | \(257500283/500000550\) | 502 | 41,663 | [abc_budget_920_delta.json.xz](abc_budget_920_delta.json.xz) |
| \(37/40\) | \(34666673/66666680\) | 316 | 40,328 | [abc_budget_925_delta.json.xz](abc_budget_925_delta.json.xz) |
| \(93/100\) | \(21/40\) | 270 | 40,132 | [abc_budget_930_delta.json.xz](abc_budget_930_delta.json.xz) |
| \(187/200\) | \(211000210979/400000400800\) | 500 | 40,318 | [abc_budget_935_delta.json.xz](abc_budget_935_delta.json.xz) |
| \(47/50\) | \(53/100\) | 839 | 41,152 | [abc_budget_940_delta.json.xz](abc_budget_940_delta.json.xz) |
| \(189/200\) | \(213/400\) | 1,275 | 51,090 | [abc_budget_945_delta.json.xz](abc_budget_945_delta.json.xz) |
| \(19/20\) | \(854400813351/1600001523200\) | 2,364 | 56,794 | [abc_budget_950_delta.json.xz](abc_budget_950_delta.json.xz) |
| \(191/200\) | \(1069/2000\) | 4,219 | 62,268 | [abc_budget_955_delta.json.xz](abc_budget_955_delta.json.xz) |
| \(24/25\) | \(106980070711/200000132200\) | 6,596 | 95,985 | [abc_budget_960_delta.json.xz](abc_budget_960_delta.json.xz) |
| \(193/200\) | \(52194191487/97561060000\) | 9,412 | 123,485 | [abc_budget_965_delta.json.xz](abc_budget_965_delta.json.xz) |
| \(97/100\) | \(267/500\) | 12,850 | 81,482 | [abc_budget_970_delta.json.xz](abc_budget_970_delta.json.xz) |
| \(39/40\) | \(213980180813/400000338000\) | 14,609 | 89,096 | [abc_budget_975_delta.json.xz](abc_budget_975_delta.json.xz) |
| \(391/400\) | \(85592071639/160000133920\) | 15,796 | 97,924 | [abc_budget_9775_delta.json.xz](abc_budget_9775_delta.json.xz) |
| \(49/50\) | \(1069/2000\) | 17,382 | 118,559 | [abc_budget_980_delta.json.xz](abc_budget_980_delta.json.xz) |
| \(1961/2000\) | \(3566601579/6666670000\) | 17,189 | 118,548 | [abc_budget_9805_delta.json.xz](abc_budget_9805_delta.json.xz) |
| \(981/1000\) | \(1069/2000\) | 17,834 | 128,858 | [abc_budget_981_delta.json.xz](abc_budget_981_delta.json.xz) |
| \(393/400\) | \(213885997/400001200\) | 18,323 | 149,081 | [abc_budget_9825_delta.json.xz](abc_budget_9825_delta.json.xz) |
| \(983/1000\) | \(534899976623/1000000003000\) | 18,391 | 149,187 | [abc_budget_983_delta.json.xz](abc_budget_983_delta.json.xz) |
| \(123/125\) | \(8916501783/16666670000\) | 18,736 | 165,396 | [abc_budget_984_delta.json.xz](abc_budget_984_delta.json.xz) |
| \(197/200\) | \(4279923132521/8000005856000\) | 19,130 | 211,584 | [abc_budget_985_delta.json.xz](abc_budget_985_delta.json.xz) |
| \(493/500\) | \(106998000419/200000000800\) | 19,526 | 301,028 | [abc_budget_986_delta.json.xz](abc_budget_986_delta.json.xz) |
| \(987/1000\) | \(534990173281/1000000324000\) | 19,888 | 407,079 | [abc_budget_987_delta.json.xz](abc_budget_987_delta.json.xz) |
| \(247/250\) | \(5349904295969/10000008030000\) | 20,282 | 482,561 | [abc_budget_988_delta.json.xz](abc_budget_988_delta.json.xz) |
| \(989/1000\) | \(21399616776161/40000031360000\) | 20,676 | 583,496 | [abc_budget_989_delta.json.xz](abc_budget_989_delta.json.xz) |
| \(99/100\) | \(26749499961979/50000000001500\) | 21,103 | 743,497 | [abc_budget_990_delta.json.xz](abc_budget_990_delta.json.xz) |

The [standard verifier](verify_abc_height_fulloc_certificate.py) pins
the source archive and selected tree hashes. It checks the root-row
change, applies the sparse replacements, and replays all branches and
rational leaf duals in the resulting tree. The replacement rows include
primitive ratio locks, fixed Hölder moments and, in some trees, the
[local Mordell packing bound](ABC_MORDELL_LOCAL_PACKING.md).

For example, replay the budget \(99/100\) tree from the repository root.

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 99/100
```

### Reweighting between and beyond the new budgets

The same \(47/50\) tree gives exact bounds \(263837/495000\) at
\(471/500\) and \(1411547/2640000\) at \(9429/10000\).
For a single certificate the reweighted bound is the maximum of the
affine leaf bounds \(b+y_0(\lambda-\lambda_0)\), hence is convex.
The straight lines between the certified points
\[
 (47/50,53/100),\quad
 (471/500,263837/495000),\quad
 (9429/10000,1411547/2640000)
\]
therefore give upper bounds throughout this interval. In particular the
reweighted count beats the unit-budget bound through \(9429/10000\).

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 47/50 --at-budget 471/500
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 47/50 --at-budget 9429/10000
```

### Reweighting the \(193/200\) certificate

Reweighting the verified \(193/200\) tree gives the exact
count exponent \(133749764986241/250000000750000\) at
\(\lambda=482503/500000\). At that budget the sparse strict unit tree
reweights to \(535000193634363/1000000462000000\), a strictly
larger exponent.

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 193/200 --at-budget 482503/500000
```

### Reweighting the \(97/100\) certificate

The \(97/100\) tree gives the exact upper exponents \(38861/72800\)
at \(\lambda=387/400\) and \(2606510113/4872000000\) at
\(\lambda=970691/1000000\). The sparse strict unit tree gives
\(1070000390799111/2000000924000000\) at the latter budget,
a strictly larger exponent.

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 97/100 --at-budget 387/400
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 97/100 --at-budget 970691/1000000
```

### Reweighting the \(39/40\) certificate

The \(39/40\) tree gives exact upper exponents
\(106950853/200000400\) at \(\lambda=389/400\) and
\(835936493193/1562500009375\) at \(\lambda=60941/62500\).
At the latter budget the sparse strict unit tree gives
\(66875024594361/125000057750000\), a strictly larger exponent.

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 39/40 --at-budget 389/400
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 39/40 --at-budget 60941/62500
```

### Reweighting the \(391/400\) certificate

The \(391/400\) tree gives exact upper exponents
\(7132667863/13333337200\) at \(\lambda=39/40\) and
\(53499982624411/100000000300000\) at \(\lambda=48877/50000\).
At the latter budget the sparse strict unit tree gives
\(53500019752617/100000046200000\), a strictly larger exponent.

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 391/400 --at-budget 39/40
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 391/400 --at-budget 48877/50000
```

### Reweighting the \(49/50\) certificate

The \(49/50\) tree gives the following exact upper exponents
at the indicated budgets.

| Radical budget | Count exponent |
| --- | --- |
| \(391/400\) | \(855199842919/1600000001600\) |
| \(979/1000\) | \(2137999607311/4000000004000\) |
| \(980011/1000000\) | \(534529057150099/1000000005000000\) |

At the last budget the sparse strict unit tree gives
\(1070000396586831/2000000924000000\).
The separate \(1961/2000\) retarget extends the strict
improvement farther. At its own budget the sparse strict unit
tree gives \(2140000793781/4000001848000\).

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 49/50 --at-budget 391/400
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 49/50 --at-budget 979/1000
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 49/50 --at-budget 980011/1000000
```

### Reweighting the \(393/400\) certificate

At \(\lambda=544721207/554858000\), the \(393/400\) tree gives the exact
upper exponent \(1069/2000\). At \(\lambda=196501/200000\),
its exact upper exponent is \(171793213/321280000\). The sparse
strict unit tree reweighted at that budget gives \(214000079627121/400000184800000\),
so the dedicated retarget has the smaller exponent. These are the
maxima of the rational reweighted leaf bounds.

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 393/400 --at-budget 544721207/554858000
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 393/400 --at-budget 196501/200000
python3 -B research/arithmetic/verify_abc_unit_strict_delta.py --at-budget 196501/200000
```

### Reweighting the \(983/1000\) certificate

At radical budget \(983003/1000000\), the \(983/1000\) retarget
gives exact upper exponent \(5802594141/10847000000\).
The sparse strict unit tree gives \(1070000398444863/2000000924000000\)
at the same budget. Both are maxima of the respective rational
reweighted leaf bounds.

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 983/1000 --at-budget 983003/1000000
python3 -B research/arithmetic/verify_abc_unit_strict_delta.py --at-budget 983003/1000000
```

### Reweighting the \(123/125\) certificate

At radical budget \(12300001/12500000\), the \(123/125\) retarget
gives exact upper exponent \(835921978968881/1562500009375000\).
The sparse strict unit tree gives \(13375004988300621/25000011550000000\)
at the same budget. Both are maxima of the respective rational
reweighted leaf bounds.

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 123/125 --at-budget 12300001/12500000
python3 -B research/arithmetic/verify_abc_unit_strict_delta.py --at-budget 12300001/12500000
```

### Reweighting the \(197/200\) certificate

At radical budget \(985001/1000000\), the \(197/200\) retarget
gives the exact upper exponent \(534991037592609/1000000003000000\).
The sparse strict unit tree gives \(1070000399685621/2000000924000000\)
at that budget. Both are maxima of the respective rational
reweighted leaf bounds.

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 197/200 --at-budget 985001/1000000
python3 -B research/arithmetic/verify_abc_unit_strict_delta.py --at-budget 985001/1000000
```

### Reweighting the \(493/500\) certificate

At radical budget \(986001/1000000\), the \(493/500\) retarget
gives the exact upper exponent \(267495625400507/500000002000000\).
The sparse strict unit tree gives \(1070000400306621/2000000924000000\)
at that budget. Both are maxima of the respective rational
reweighted leaf bounds.

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 493/500 --at-budget 986001/1000000
python3 -B research/arithmetic/verify_abc_unit_strict_delta.py --at-budget 986001/1000000
```

### Reweighting the \(987/1000\) certificate

At radical budget \(4935001/5000000\), the \(987/1000\) retarget
gives the exact upper exponent \(2431773951451/4545455000000\).
The sparse strict unit tree gives \(5350002004635621/10000004620000000\)
at that budget. Both are maxima of the respective rational
reweighted leaf bounds.

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 987/1000 --at-budget 4935001/5000000
python3 -B research/arithmetic/verify_abc_unit_strict_delta.py --at-budget 4935001/5000000
```

### Reweighting the \(247/250\) certificate

At radical budget \(9880003/10000000\), the \(247/250\) retarget
gives the exact upper exponent \(53499028958498837/100000000500000000\).
The sparse strict unit tree gives \(10700004015481863/20000009240000000\)
at that budget. Both are maxima of the respective rational
reweighted leaf bounds.

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 247/250 --at-budget 9880003/10000000
python3 -B research/arithmetic/verify_abc_unit_strict_delta.py --at-budget 9880003/10000000
```

### Reweighting the \(989/1000\) certificate

At radical budget \(4945001/5000000\), the \(989/1000\) retarget
gives the exact upper exponent \(441830515859/825857500000\).
The sparse strict unit tree gives \(5350002010845621/10000004620000000\)
at that budget. Both are maxima of the respective rational
reweighted leaf bounds.

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 989/1000 --at-budget 4945001/5000000
python3 -B research/arithmetic/verify_abc_unit_strict_delta.py --at-budget 4945001/5000000
```

### Reweighting the \(99/100\) certificate

At radical budget \(9900001/10000000\), the \(99/100\) retarget
gives the exact upper exponent \(213996134798720789/400000197200000000\).
The sparse strict unit tree gives \(10700004027900621/20000009240000000\)
at that budget. Both are maxima of the respective rational
reweighted leaf bounds.

```sh
python3 -B research/arithmetic/verify_abc_height_fulloc_certificate.py --source-budget 99/100 --at-budget 9900001/10000000
python3 -B research/arithmetic/verify_abc_unit_strict_delta.py --at-budget 9900001/10000000
```

### A twelve-leaf improvement near unit radical budget
At radical budget \(991/1000\), retarget the pinned strict unit tree and replace the
12 leaves whose exact reweighted bounds exceed the new ceiling. Each replacement
is a rational dual for the same leaf rows. The ordinary exact replay checked all
1,830,732 leaves and certified the upper exponent \(267494705491/500000001000\).
The 1,500-byte [overlay](abc_unit_budget_991_overlay.json.xz) has SHA-256
`48cba103fc4dad46df7f5e7c01de2e2fd9abbe33fa7411d718083182e560ea6b`.

Every verified leaf gives an affine bound in the budget, with its slope equal to
the multiplier of the radical budget row. The maximum of these lines is convex.
At \(57804903499/58300100000\) it equals \(53499/100000\),
and each line is nondecreasing, giving that ceiling at lower budgets.
The old tree has a leaf line \(267499947/500000231+(621/2000000924)\lambda\).
The new maximum is strictly below this line at \(991/1000\) and meets it at \(130801157300042477/130801188430064625\),
where the old maximum also equals the line. Convexity gives a strict improvement on
the whole intervening half-open interval.

| Radical budget | Overlay exponent | Same-budget unit exponent |
|---|---|---|
| \(9910001/10000000\) | \(891649018497667/1666666670000000\) | \(10700004034110621/20000009240000000\) |
| \(199/200\) | \(17833135833/33333333400\) | \(214000081179/400000184800\) |
| \(9999/10000\) | \(2674998941999/5000000010000\) | \(10700004089379/20000009240000\) |

```sh
python3 -B research/arithmetic/verify_abc_unit_strict_delta.py --overlay research/arithmetic/abc_unit_budget_991_overlay.json.xz --at-budget 9999/10000
python3 -B research/arithmetic/verify_abc_unit_strict_delta.py --at-budget 9999/10000
```

Add `--full` to replay all inherited rows at the source of the overlay.

### A strict improvement through unit radical budget

The [combined unit overlay](abc_unit_budget_1_interval_overlay.json.xz)
replaces 5,199 dual leaves in the existing strict unit tree and retains
all 1,830,732 leaf positions. Its SHA-256 is
`e6cd907dae663245adbdcd00fa0535c7844fba54eb8bcdb98fa664880360a441`.
At unit budget the full exact replay and sparse replay both give
\(42800023737/80000047360\), below the previous
\(7133338393/13333343040\).

The same overlay strictly improves the old unit tree on the closed interval \([991/1000,1]\).
The checker takes an ordered list of old knots. For each consecutive pair
it checks that a single old leaf line attains the old maximum at both ends.
Since the old maximum is convex and lies above that line, it agrees with
the line on that segment. The new maximum is also convex. Checking a
strict improvement at each knot proves it on the whole interval.

| Radical budget | Combined overlay exponent | Earlier unit tree |
|---|---|---|
| \(991/1000\) | \(267494705491/500000001000\) | \(1070000403411/2000000924000\) |
| \(130801157300042477/130801188430064625\) | \(13995732829209764109620851/26160250844612481064501275\) | \(93304840305054369/174401584573419500\) |
| \(171787695125876029/171787696585878054\) | \(45953234676129959392173581/85893900172823395935172308\) | \(245083774666788733/458100524229008144\) |
| \(24267451106599229/24267451126599274\) | \(25966186573679728296885923/48534930015162636829569456\) | \(9442244437864069/17649055364799472\) |
| \(648493562008930961/648493562183927025\) | \(173472120590860320060035111/324246966561122297103129150\) | \(138777619678434073/259397424873570810\) |
| \(453862750531578208/453862750531586879\) | \(19425336496091672059495110223/36309041537466815495954589440\) | \(173440404995357859/324187678951133485\) |
| \(1\) | \(42800023737/80000047360\) | \(7133338393/13333343040\) |

```sh
python3 -B research/arithmetic/verify_abc_unit_strict_delta.py \
  --overlay research/arithmetic/abc_unit_budget_1_interval_overlay.json.xz \
  --compare-unit-knots \
  991/1000 130801157300042477/130801188430064625 \
  171787695125876029/171787696585878054 24267451106599229/24267451126599274 \
  648493562008930961/648493562183927025 453862750531578208/453862750531586879 \
  1
```

Add `--full` to replay all inherited duals at unit budget as well.
