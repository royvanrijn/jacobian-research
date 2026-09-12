# Degree-nine marked-root obstruction search

> **Status and scope.** This is an exact characteristic-zero construction and
> relative-rank certificate, followed by modular order-five component
> discovery. It is not yet a characteristic-zero reconstruction of the
> Fitting locus. Consequently no degree-nine order-seven PBW computation,
> Ore-localized quantization, conductor gluing, Weyl endomorphism, or
> nonsurjectivity claim is made here.

## 1. Classical two-parameter family

Let

\[
H_{\sigma,\tau}(W)=W^2(W-1)P_{\sigma,\tau}(W),
\]

where

\[
P_{\sigma,\tau}(W)=W^6+\sigma W^5+\tau W^4
+\left(-\frac{17}{2}-5\sigma-4\tau\right)W
+\frac{13}{2}+4\sigma+3\tau.
\tag{1.1}
\]

Then `P(1)=-1` and `P'(1)=-5/2`, so the marked incidences remain
`H'(1)=-1` and `H''(1)=-9`. Six exact Hamiltonian-homotopy samples and an
unused seventh sample determine and verify the completing shear

\[
q_9(\sigma,\tau)=\frac3{476476}
\left(
2437664\sigma^2+3435360\sigma\tau+8147828\sigma
+1214208\tau^2+5770548\tau+5171465
\right).
\tag{1.2}
\]

This gives an explicit two-parameter family of noninjective polynomial
Poisson maps in the same normalized marked-root incidence class as the
degree-seven and degree-eight rows.

## 2. Relative restricted deformation complex

The inherited root-weight supports have dimensions

\[
\#S_2=132,\quad \#T_2=107,\qquad
\#S_4=84,\quad \#T_4=64.
\tag{2.1}
\]

The finite relative complex is organized as

\[
\text{order-three gauge/lift directions}\longrightarrow
\text{order-five corrections}\longrightarrow
\text{order-five defects}.
\tag{2.2}
\]

At `(sigma,tau)=(1,0)`, exact reduction over `Q` gives

| quantity | value |
|---|---:|
| order-three correction dimension | 239 |
| `rank d3` | 227 |
| `dim ker d3` | 12 |
| order-five correction dimension | 148 |
| `rank D5` | 142 |
| `dim ker D5` | 6 |
| quadratic lower-lift columns | 90 |
| `rank M5` | 149 |
| `rank [M5|O0]` | 150 |
| defect-space dimension | 371 |
| coherent obstruction-module dimension `dim coker M5` | 222 |

Thus the constant order-five class is generically nonzero. The same complete
signature is reproduced over `GF(23)`, `GF(29)`, `GF(31)`, and `GF(37)`.
At `p=19`, `rank d3=226` and `dim ker d3=13`; that reduction is excluded.
At `p=17`, the shear formula itself has a pole.

These values extend the exact degree ladder to a fourth row:

\[
\dim\ker d_3=2n-6,\qquad
\dim\ker D_5=2n-12,\qquad
\operatorname{rank}M_5-\operatorname{rank}D_5=7.
\tag{2.3}
\]

Equation (2.3) remains a four-case conjectural pattern, not an all-degree
theorem.

## 3. Modular Fitting and nonlinear screens

The order-five restricted Fitting condition is that adjoining the constant
defect does not raise the strong rank. Complete finite parameter-plane scans
give

| prime | strong-consistent base points |
|---:|---|
| 23 | `(10,8)`, `(14,14)` |
| 29 | none |
| 31 | `(9,22)` |
| 37 | none |

The nonlinear Kuranishi screen projects the genuine quadratic lower-lift
equations modulo every current correction. Singular finds a nonempty
six-dimensional fibre with a nine-element standard basis above all three
recorded base points. Its radical basis has six affine-linear generators,
so the reduced fibre is

\[
\mathbb A^6=\operatorname{Spec}k[z_0,z_2,z_4,z_6,z_8,z_{10}].
\tag{3.1}
\]

The leading monomial ideal proves that the full fibre is finite free of rank
six over this radical. In every fibre, the univariate cubic in `z11` is a
perfect cube. Thus the three points survive the genuine nonlinear order-five
test and carry a certified nonreduced rank-six thickening; they are not
artifacts of treating the quadratic lift monomials as independent columns.

Counts `2,0,1,0` are compatible with a finite arithmetic scheme having different
splitting types, but they do not determine its degree or prove that the three
fibres descend from one characteristic-zero component.

Across degrees seven through nine, the reduced nonlinear base now has the
uniform form

\[
\mathbb A^{n-3}=\operatorname{Spec}k[z_0,z_2,\ldots,z_{2n-8}],
\tag{3.2}
\]

and the known thickening ranks are `2,4,6`. The degree-seven/eight statements
are characteristic zero; the degree-nine statement is presently certified
only on the three modular fibres. The evident rank law `2(n-6)` is therefore
a new three-row conjecture, not yet a family-level theorem.

## 4. Pivot-chart degree law

A generic specialization selects 149 independent strong columns among 238
columns in a 371-dimensional defect support. Solving the first 149 pivot rows
and probing 166 points on the line `tau=0` over `GF(1009)` proves, with unused
samples in every rational interpolation, that the common chart denominator
has degree

\[
72=2(9-3)^2.
\tag{4.1}
\]

The first weight row of residual numerators has successive degrees

\[
74,76,78,80,82,84,86,
\tag{4.2}
\]

and the next row begins

\[
73,75,77,79,81.
\tag{4.3}
\]

Together with denominator degrees 32 and 50 at degrees seven and eight,
(4.1) establishes the three-row modular pattern `2(n-3)^2`. This is a sharp
interpolation bound and a computational optimization. It is not a proof of
the corresponding formula for arbitrary degree.

## 5. Authorization boundary

The next authorized calculation is the denominator-72 bivariate pivot-chart
interpolation over several large primes, followed by stable rational
reconstruction and an unused-prime check. Only a component reconstructed in
that way may receive characteristic-zero nonlinear lifting and order-seven
PBW computation.

The root-at-infinity valuation ledger itself is already fixed by the inherited
grading:

| coefficient order | `S` valuation | `T` valuation | support size |
|---:|---:|---:|---:|
| classical | -2 | -1 | classical symbol |
| `hbar^2` | 4 | 5 | `132+107` |
| `hbar^4` | 10 | 11 | `84+64` |
| `hbar^6` | 16 | 17 | `46+31` |

The last row records the complete possible inherited order-seven correction
support; it is not an order-seven consistency calculation.

In particular, the terminal ledger is currently:

| requested step | degree-nine status |
|---|---|
| reconstruct a component over `Q` | pending modular chart reconstruction |
| solve Ore-localized quantization | unauthorized before reconstruction |
| compute root-at-infinity valuations | complete inherited ledger through `hbar^6` above |
| conductor gluing | not reached |
| certify Weyl relations | not reached |
| certify nonsurjectivity | not reached |

This keeps the logical direction explicit: the computation is selecting
classical symbols on which the restricted obstruction class vanishes, not
attempting to repair the already closed degree-six/eight symbols.

## 6. Reproduction

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/derive_degree_nine_marked_root_shear.py --jobs 7 \
  --output artifacts/generated-results/degree_nine_marked_root_shear.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_nine_relative_quantization_obstruction.py \
  --output artifacts/generated-results/degree_nine_relative_quantization_obstruction.json

for prime in 23 29 31 37; do
  PYTHONPATH=scripts .venv/bin/python \
    scripts/search_degree_seven_order_five_fitting_locus.py \
    --degree 9 --prime "$prime" --jobs 8 \
    --output \
    "artifacts/generated-results/degree_nine_order_five_scan_gf${prime}.json"
done

PYTHONPATH=scripts .venv/bin/python \
  scripts/screen_degree_seven_order_five_survivors.py \
  artifacts/generated-results/degree_nine_order_five_scan_gf23.json \
  artifacts/generated-results/degree_nine_order_five_scan_gf31.json \
  --output artifacts/generated-results/degree_nine_order_five_nonlinear_screen.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/probe_degree_nine_order_five_chart.py --prime 1009 --jobs 12 \
  --output artifacts/generated-results/degree_nine_order_five_chart_degree_probe.json
```

The nonlinear screen requires Singular. Every scan and chart probe is modular
discovery; the exact relative-rank certificate is characteristic zero.

The prepared but not yet completed large interpolation batch is

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/interpolate_degree_nine_order_five_chart.py \
  --prime 1009 --jobs 12 \
  --output artifacts/generated-results/degree_nine_order_five_chart_gf1009.json
```

It evaluates a `150 x 84` grid, validates four residual numerators on held-out
grid points, saturates by the degree-72 denominator, and requires Singular.
