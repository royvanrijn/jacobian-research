# Degree-eight marked-root obstruction component

> **Status and scope.** This note gives an exact characteristic-zero
> calculation in the parity-preserving, root-weight-homogeneous restricted
> PBW complex. The degree-eight order-five section-zero scheme is one
> irreducible closed point of degree twelve. Its genuine reduced lift is an
> affine five-space over that residue field, with two additional nilpotent
> repeated-root directions. The complete inherited order-seven ideal is the
> unit ideal. Thus this component does not produce a Weyl-algebra
> endomorphism and does not settle `DC_2`.

## 1. Classical family and completing shear

Put

\[
H_{\sigma,\tau}(W)=W^2(W-1)P_{\sigma,\tau}(W),
\]

where

\[
P_{\sigma,\tau}(W)=W^5+\sigma W^4+\tau W^3
+\left(-\frac{15}{2}-4\sigma-3\tau\right)W
+\frac{11}{2}+3\sigma+2\tau.
\tag{1.1}
\]

Then `P(1)=-1` and `P'(1)=-5/2`, so the marked incidence has
`H'(1)=-1` and `H''(1)=-9`. Seven exact Hamiltonian-homotopy evaluations,
six for interpolation and one held out, give the completing shear

\[
q_8(\sigma,\tau)=\frac3{28028}
\left(
71424\sigma^2+86112\sigma\tau+279300\sigma
+26104\tau^2+169352\tau+176269
\right).
\tag{1.2}
\]

This produces an explicit noninjective two-parameter polynomial Poisson
family. The exact classical interpolation is independent of the PBW
obstruction calculation below.

## 2. Relative restricted complex

At orders two and four, the correction supports are

\[
\#S_2=100,\quad \#T_2=78,\qquad
\#S_4=59,\quad \#T_4=42.
\]

The relative complex is the finite complex

\[
\text{gauge}\longrightarrow
\text{corrections}\longrightarrow
\text{defects},
\]

with the order-three affine solution torsor supplying the lower-lift
coordinates and the cokernel of the current correction map supplying the
coherent obstruction module. At the rational point `(sigma,tau)=(1,0)`,
exact row reduction gives

| quantity | value |
|---|---:|
| `rank d3` | 168 |
| `dim ker d3` | 10 |
| `rank D5` | 97 |
| `dim ker D5` | 4 |
| `rank M5` | 104 |
| `rank [M5|O0]` | 105 |
| order-five output support | 269 |

The same signature holds at `19,23,29,31`. The reduction at `17` is bad:
`rank d3` falls to 167 there, so that prime is excluded from component
inference.

## 3. Order-five Fitting locus and reconstruction

Full-plane scans give the following section-zero points:

| prime | points `(sigma,tau)` |
|---:|---|
| 19 | none; bad model reduction for the eventual lex presentation |
| 23 | `(20,12)` |
| 29 | `(9,14)`, `(21,2)` |
| 31 | `(8,14)` |
| 37 | none |
| 41 | `(25,6)` |
| 43 | `(29,21)` |
| 47 | `(2,4)`, `(37,42)` |

On a fixed 104-column pivot chart, three saturated residual numerators
interpolate a five-polynomial standard basis of length twelve. The exact
descent uses three 31-bit images and three additional build primes. The first
three normalized standard-basis polynomials stabilize under removal of the
last build prime. Those three already generate the entire ideal over `Q`;
exact Buchberger completion recovers the two higher-height polynomials. The
completed basis reduces to every build image and to the unused `GF(1009)`
image exactly.

The lexicographic elimination polynomial is the irreducible polynomial

\[
\begin{aligned}
f(\sigma)={}&9699739875\sigma^{12}+41496410970\sigma^{11}
-846610243920\sigma^{10}-10549859894154\sigma^9\\
&-56081730021765\sigma^8-173027569716540\sigma^7
-339093184426920\sigma^6-467796452829460\sigma^5\\
&-560246361633175\sigma^4-553097854439102\sigma^3
-51771900582024\sigma^2+380750072251710\sigma\\
&-92072207284295.
\end{aligned}
\tag{3.1}
\]

The second lexicographic equation is linear in `tau`; its exact coefficients
are recorded in
`degree_eight_order_five_rational_chart.json`. Thus the Fitting locus is one
irreducible degree-twelve point over `Q`, not merely a compatible collection
of finite-field points. Exact rank 104 of the selected columns over this
residue field proves that the component lies inside the reconstruction chart.

## 4. Genuine nonlinear lift

Over `K=Q[sigma]/(f)`, the 172 projected quadratic Kuranishi equations span a
seven-dimensional polynomial space. Seven pivot equations therefore generate
the same ideal. Its reduced standard basis has six elements:

- three affine-linear relations eliminate `z1,z3,z5`;
- one quadratic in `z9` has zero discriminant;
- one quadratic in `z7` over `K[z8]` has identically zero discriminant;
- the bridge relation vanishes after substituting both repeated roots.

Consequently the reduced order-five lift is

\[
\mathbb A^5_K
=\operatorname{Spec}K[z_0,z_2,z_4,z_6,z_8],
\tag{4.1}
\]

while `z7,z9` are nilpotent repeated-root directions. The scheme is
nonreduced, but it introduces no further reduced residue-field extension.

## 5. Complete order-seven obstruction

The order-seven current correction supports and rank are

\[
\#S_6=28,\qquad \#T_6=16,qquad \operatorname{rank}D_7=44.
\]

All four free directions of `ker D5` are retained. Together with the five
coordinates in (4.1), the obstruction is cubic in nine variables. Exact
Newton interpolation uses all

\[
\binom{9+3}{3}=220
\]

degree-at-most-three nodes and an unused 221st point. The order-seven output
has dimension 196; projection modulo `D7` gives 152 equations, whose
polynomial span has six independent generators. Singular computes

\[
\operatorname{std}(I_7)=(1).
\tag{5.1}
\]

The unit result is first proved on the reduced lift. It also excludes the
nonreduced order-five thickening: if an ideal becomes the unit ideal modulo a
nilpotent ideal, then a relation `1+n` lifts, and `1+n` is a unit.

## 6. Valuations and terminal descent ledger

Every retained correction is homogeneous for the root-at-infinity valuation
`v(X)=1`, `v(Q)=-1`, `v(Z)=-2`:

| term | `S` side | `T` side |
|---|---:|---:|
| classical | -2 | -1 |
| `hbar^2` | 4 | 5 |
| `hbar^4` | 10 | 11 |
| `hbar^6` | 16 | 17 |

Hence every root-at-infinity valuation through the last solvable correction
order is explicit. The requested terminal steps resolve as follows:

| step | result |
|---|---|
| reconstruct component over `Q` | one irreducible degree-twelve point |
| solve Ore-localized quantization | exact through `hbar^6`; inconsistent at `hbar^7` |
| compute root-at-infinity valuations | complete table above |
| conductor gluing | not applicable: no order-seven local quantization |
| certify Weyl relations | impossible on this component at order seven |
| certify nonsurjectivity | not applicable: no Weyl endomorphism survives |

Thus the degree-eight family answers the revised question decisively: its
restricted obstruction class vanishes through order five on one arithmetic
component, but the coherent order-seven class does not vanish anywhere on
the genuine lift scheme.

## 7. Reproduction

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/derive_degree_eight_marked_root_shear.py --jobs 7 \
  --output artifacts/generated-results/degree_eight_marked_root_shear.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_eight_relative_quantization_obstruction.py \
  --output \
  artifacts/generated-results/degree_eight_relative_quantization_obstruction.json

for prime in 19 23 29 31 37 41 43 47; do
  PYTHONPATH=scripts .venv/bin/python \
    scripts/search_degree_seven_order_five_fitting_locus.py \
    --degree 8 --prime "$prime" --jobs 8 \
    --output \
    "artifacts/generated-results/degree_eight_order_five_scan_gf${prime}.json"
done

PYTHONPATH=scripts .venv/bin/python \
  scripts/screen_degree_seven_order_five_survivors.py \
  artifacts/generated-results/degree_eight_order_five_scan_gf*.json \
  --output \
  artifacts/generated-results/degree_eight_order_five_nonlinear_screen.json

for prime in 2147483647 2099999999 2049999979 1019 1021 1013 1009; do
  PYTHONPATH=scripts .venv/bin/python \
    scripts/interpolate_degree_eight_order_five_chart.py \
    --prime "$prime" --jobs 12 \
    --output \
    "artifacts/generated-results/degree_eight_order_five_chart_gf${prime}.json"
done

PYTHONPATH=scripts .venv/bin/python \
  scripts/reconstruct_degree_eight_order_five_rational_chart.py \
  --holdout-prime 1009 \
  --output \
  artifacts/generated-results/degree_eight_order_five_rational_chart.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_eight_order_five_survivor.py \
  --output \
  artifacts/generated-results/degree_eight_order_five_survivor.json
```

The reconstruction and survivor commands require Singular. The last command
is the expensive exact degree-twelve calculation.
