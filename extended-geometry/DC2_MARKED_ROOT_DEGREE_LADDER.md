# Marked-root degree ladder for the restricted `DC_2` search

> **Status and scope.** This note isolates exact structural patterns shared by
> the degree-six through degree-nine marked-root calculations.
> The support formulas are proved combinatorially. The rank formulas are exact
> computations in the four displayed degrees, not an all-degree theorem.
> Degrees six through eight are closed at order seven. Their section-zero
> components have residue degrees four, eight, and twelve; every complete
> inherited order-seven ideal is the unit ideal. Degree nine has passed the
> exact relative-rank calculation and multi-prime order-five discovery screen,
> but its characteristic-zero section-zero scheme is not yet reconstructed.

## 1. Uniform two-parameter marked-root slices

For `n>=7`, put `m=n-3` and

\[
H^{(n)}_{\sigma,\tau}(W)=W^2(W-1)P^{(m)}_{\sigma,\tau}(W),
\tag{1.1}
\]

where

\[
\begin{aligned}
P^{(m)}_{\sigma,\tau}(W)={}&W^m+\sigma W^{m-1}+\tau W^{m-2}+A_mW+B_m,\\
A_m={}&-\frac{2m+5}{2}-(m-1)\sigma-(m-2)\tau,\\
B_m={}&\frac{2m+1}{2}+(m-2)\sigma+(m-3)\tau.
\end{aligned}
\tag{1.2}
\]

All omitted intermediate coefficients are zero. Direct substitution gives

\[
P(1)=-1,\qquad P'(1)=-\frac52,
\tag{1.3}
\]

and hence

\[
H(0)=H'(0)=H(1)=0,qquad H'(1)=-1,qquad H''(1)=-9.
\tag{1.4}
\]

Thus every row lies in the same `kappa=-9` marked-root incidence class while
changing the classical degree. Exact homotopy interpolation gives quadratic
completing shears. The first two monic rows are

\[
q_7=\frac3{28028}(26104\sigma^2+21736\sigma\tau+134160\sigma
+4576\tau^2+56160\tau+75285)
\tag{1.5}
\]

and

\[
q_8=\frac3{28028}(71424\sigma^2+86112\sigma\tau+279300\sigma
+26104\tau^2+169352\tau+176269).
\tag{1.6}
\]

The next row is

\[
q_9=\frac3{476476}(2437664\sigma^2+3435360\sigma\tau
+8147828\sigma+1214208\tau^2+5770548\tau+5171465).
\tag{1.7}
\]

Each formula uses six exact samples and an unused seventh sample.

## 2. Closed support-count formula

The classical symbol bounds in degree `n` are

\[
(\deg_Z,\deg_B)(S)=(n,7n-6),\qquad
(\deg_Z,\deg_B)(T)=(n-1,7n-10).
\tag{2.1}
\]

At correction order `2r`, the bounds decrease by `(2r,4r)` and the weights
are `-2+6r` for `S_(2r)` and `-1+6r` for `T_(2r)`. Writing
`x-q-2z=nu`, the number of allowed monomials is therefore

\[
\boxed{
N^S_{n,r}=\sum_{z=0}^{n-2r}
\left(\left\lfloor\frac{7n-4-10r-5z}{2}\right\rfloor+1\right),}
\tag{2.2}
\]

and

\[
\boxed{
N^T_{n,r}=\sum_{z=0}^{n-1-2r}
\left(\left\lfloor\frac{7n-9-10r-5z}{2}\right\rfloor+1\right).}
\tag{2.3}
\]

These formulas are not extrapolations: they follow by writing
`x=q+nu+2z` and counting
`0<=q<=floor((deg_B-nu-5z)/2)`.

They also admit a compact closed form. Put `L=n-2r` and `K=n-1-2r`.
Because the numerator in each floor has the same parity as its upper summation
limit, the accumulated parity correction is `ceil(L/2)` or `ceil(K/2)`.
Consequently

\[
N^S_{n,r}=
\frac{(L+1)(9L+8r-4)-2\lceil L/2\rceil}{4},\qquad
N^T_{n,r}=
\frac{(K+1)(9K+8r)-2\lceil K/2\rceil}{4}.
\tag{2.4}
\]

Thus every correction-space size is an explicit quadratic
quasi-polynomial of period two in the classical degree. This is useful beyond
bookkeeping: if the conjectural kernel law in (3.1) holds, then the complete
rank of `D5` follows immediately from (2.4) as
`N^S_(n,2)+N^T_(n,2)-(2n-12)`.

For orders two and four they give

| degree | `#S2` | `#T2` | `#S4` | `#T4` |
|---:|---:|---:|---:|---:|
| 6 | 49 | 34 | 22 | 12 |
| 7 | 72 | 54 | 38 | 25 |
| 8 | 100 | 78 | 59 | 42 |
| 9 | 132 | 107 | 84 | 64 |

## 3. Exact rank ladder

At the rational point `(sigma,tau)=(1,0)`, exact characteristic-zero row
reduction gives:

| degree | `rank d3` | `dim ker d3` | `rank D5` | `dim ker D5` | `rank M5` | `rank [M5|O0]` | output |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 | 77 | 6 | 34 | 0 | 41 | 42 | 110 |
| 7 | 118 | 8 | 61 | 2 | 68 | 69 | 178 |
| 8 | 168 | 10 | 97 | 4 | 104 | 105 | 269 |
| 9 | 227 | 12 | 142 | 6 | 149 | 150 | 371 |

The four exact rows satisfy

\[
\dim\ker d_3=2n-6,qquad
\dim\ker D_5=2n-12,qquad
\operatorname{rank}M_5-\operatorname{rank}D_5=7.
\tag{3.1}
\]

Equation (3.1) is now a sharply formulated all-degree conjecture. Proving it
requires a symbolic pivot family or a representation-theoretic decomposition;
four successful rows alone are not a proof. The stable `+7` is especially
useful: although the number of quadratic lower-lift columns grows rapidly,
only seven new defect directions survive modulo current corrections in every
computed row.

In all four degrees, adjoining the constant raises the generic strong rank
by one. Thus the generic order-five obstruction persists as degree grows; the
search problem is consistently a section-zero problem on a proper locus, not
a generic repair problem.

## 4. Arithmetic and nonlinear pattern

| degree | section-zero component | nonlinear order-five scheme | order seven |
|---:|---|---|---|
| 6 | one irreducible quartic point | `A^4` over its residue field | unit ideal |
| 7 | one irreducible octic point | doubled `A^4` | unit ideal |
| 8 | one irreducible degree-twelve point | two-square thickening of `A^5` | unit ideal |
| 9 | not yet reconstructed; modular counts `2,0,1,0` at `p=23,29,31,37` | rank-six thickening of `A^6` at every recorded point | not authorized |

The first two residue degrees `4,8` suggested doubling. The exact degree-eight
component falsifies that extrapolation and replaces it by the first three
values of a linear `+4` pattern. Its stable characteristic-zero generators,
exact Buchberger completion, and unused-prime check prove residue degree 12.
Over that field, the nonlinear standard basis has three affine-linear
relations, two repeated-root quadratics with identically zero discriminant,
and one bridge relation. The reduced lift is therefore affine five-space with
coordinates `z0,z2,z4,z6,z8`; `z7,z9` carry nilpotent square directions.

The degree-nine modular fibres sharpen the next pattern. Each of the three
recorded survivors has radical
`A^6=Spec k[z0,z2,z4,z6,z8,z10]`; its nine-element standard basis makes the
full algebra finite free of rank six over that radical, and its univariate
`z11` cubic is a perfect cube. Consequently degrees seven, eight, and nine
have reduced lift `A^(n-3)` on even lower-lift coordinates and thickening
ranks `2,4,6`. The tempting law `2(n-6)` is exact in those three fibres, but
the degree-nine row is modular and the law is not promoted to a theorem.

At order seven, all four directions of `ker D5` are retained. The 220-node
cubic interpolation has an independent holdout, and the 152 projected
equations span six polynomial generators. Their exact ideal is `(1)`. See the
[degree-eight component note](DC2_DEGREE_EIGHT_MARKED_ROOT_SEARCH.md) for the
full reconstruction, valuation, and terminal descent ledgers.

Two arithmetic warnings are exact: `p=17` lowers the degree-eight
order-three rank from 168 to 167 at `(1,0)`. It is a bad-reduction prime and
is excluded from generic splitting inference. At degree nine, `p=19` lowers
the same rank from 227 to 226, while `p=17` is unavailable because it divides
the completing-shear denominator. The good primes `23,29,31,37` retain the
characteristic-zero signature. A direct cross-degree audit falsifies the
tempting rule “the bad prime is always `2n+1`”: degree six is good at 13 but
has a strong-order-five rank drop at 17, and degree seven is good at both 17
and 19. These are presentation-specific torsion primes, not a proved uniform
arithmetic law.

The degree-nine modular pivot chart supplies another exact finite-field
pattern. Its common denominator has degree

\[
72=2(9-3)^2,
\]

continuing the degree-seven/eight values `32,50`. On an independently checked
line, successive residual numerators have degrees
`74,76,78,80,82,84,86`, followed on the next weight row by
`73,75,77,79,81`. This is a degree-bound certificate for the pending
bivariate interpolation, not a characteristic-zero zero-scheme result.

## 5. Search consequence

Degrees six through eight now give three exact terminal rows. The next work
should not enlarge PBW support on any of their components: each order-seven
ideal is already a unit. Degree nine has reached the intended order-five gate.
Three genuinely new questions remain:

1. prove or disprove the all-degree rank laws (3.1), preferably by a symbolic
   pivot family or representation-theoretic decomposition;
2. reconstruct the degree-nine denominator-72 pivot chart over several large
   primes and certify its characteristic-zero section-zero scheme;
3. authorize nonlinear characteristic-zero and order-seven work only for the
   reconstructed degree-nine components.

The arithmetic pattern to test is successive `+4` residue-degree growth, not
the retired doubling guess. The linear pattern is still a conjecture from
three values; the stable rank-seven lower-variation image is likewise not yet
an all-degree theorem.

## 6. Literature boundary

The stable Jacobian--Dixmier comparison of
[Belov-Kanel--Kontsevich](https://arxiv.org/abs/math/0512171) proceeds by
reduction to finite characteristic; it does not supply a characteristic-zero
lift for one fixed noninjective polynomial Poisson symbol. Likewise, the
polynomial lifting programme of
[Kanel-Belov--Grigoriev--Elishev--Yu--Zhang](https://arxiv.org/abs/1707.06450)
starts with polynomial *symplectomorphisms* and constructs Weyl-algebra
automorphisms by tame approximation. The marked-root rows here are deliberately
noninjective Poisson maps, so they lie outside that automorphism-lifting input.
The restricted obstruction module is therefore a complementary diagnostic,
not a competing construction of the automorphism lift.

This scope distinction also explains why the exact-quantization controls in
the family census vanish identically while the marked-root rows do not: the
controls are already automorphisms, whereas the search targets the first locus
on which a noninjective classical symbol might nevertheless satisfy the Weyl
relations.

## 7. Reproduction

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_dc2_marked_root_degree_ladder.py \
  --output artifacts/generated-results/dc2_marked_root_degree_ladder.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/derive_degree_eight_marked_root_shear.py --jobs 7 \
  --output artifacts/generated-results/degree_eight_marked_root_shear.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_eight_relative_quantization_obstruction.py \
  --output \
  artifacts/generated-results/degree_eight_relative_quantization_obstruction.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/search_degree_seven_order_five_fitting_locus.py \
  --degree 8 --prime 31 --jobs 8 \
  --output artifacts/generated-results/degree_eight_order_five_scan_gf31.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/interpolate_degree_eight_order_five_chart.py \
  --prime 1009 --jobs 16 \
  --output artifacts/generated-results/degree_eight_order_five_chart_gf1009.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/reconstruct_degree_eight_order_five_rational_chart.py \
  --holdout-prime 1009 \
  --output artifacts/generated-results/degree_eight_order_five_rational_chart.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_eight_order_five_survivor.py \
  --output artifacts/generated-results/degree_eight_order_five_survivor.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/derive_degree_nine_marked_root_shear.py --jobs 7 \
  --output artifacts/generated-results/degree_nine_marked_root_shear.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_nine_relative_quantization_obstruction.py \
  --output artifacts/generated-results/degree_nine_relative_quantization_obstruction.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/search_degree_seven_order_five_fitting_locus.py \
  --degree 9 --prime 23 --jobs 8 \
  --output artifacts/generated-results/degree_nine_order_five_scan_gf23.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/probe_degree_nine_order_five_chart.py --prime 1009 --jobs 12 \
  --output artifacts/generated-results/degree_nine_order_five_chart_degree_probe.json
```

The degree-eight interpolation command is modular discovery only; its next two
commands provide the characteristic-zero reconstruction and terminal
order-seven certificate, both requiring Singular. All degree-nine commands
shown here stop before characteristic-zero component reconstruction.
