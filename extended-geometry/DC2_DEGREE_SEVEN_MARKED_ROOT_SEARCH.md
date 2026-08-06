# Degree-seven marked-root obstruction search

> **Status and scope.** This is the next genuinely noninjective classical
> branch after the closed degree-six `MR6` component. It constructs an
> explicit two-parameter degree-seven polynomial symplectic family and its
> complete inherited parity-preserving, root-weight-homogeneous order-five
> presentation. The generic strong obstruction is nonzero. Exhaustive scans
> over twelve finite fields find sixteen isolated section-zero points. Modular
> interpolation, a 151-bit CRT, and an independent holdout prime reconstruct
> their characteristic-zero closure as one irreducible degree-eight point.
> Over its residue field, the genuine order-five lift scheme is a doubled
> affine four-space. The complete order-seven equations, including both free
> order-five correction directions, generate the unit ideal. Thus this branch
> terminates at order seven; there is no quantized map on which conductor
> gluing, Weyl relations, or nonsurjectivity could be tested.

This continues the
[classical-symbol family search](DC2_CLASSICAL_SYMBOL_FAMILY_SEARCH.md)
without returning to the obstructed quintic or enlarging correction support
on the degree-six family.

## 1. Classical family

Fix `kappa=-9` and set

\[
\begin{aligned}
H_{\sigma,\tau}(W)=W^2(W-1)\bigg(&W^4+\sigma W^3+\tau W^2\\
&+\left(-\frac{13}{2}-3\sigma-2\tau\right)W
+\frac92+2\sigma+\tau\bigg).
\end{aligned}
\tag{1.1}
\]

Then

\[
H(0)=H'(0)=H(1)=0,
\qquad H'(1)=-1,
\qquad H''(1)=-9.
\tag{1.2}
\]

The leading coefficient is one, so every point has exact degree seven. This
is a two-dimensional affine slice of the three-dimensional fixed-`kappa`
marked-root family. The all-degree rank-two descent theorem applies. Exact
quadratic interpolation of its uniform four-residue functional gives the
completing shear

\[
\boxed{
s_2=\frac{3}{28028}\big(
26104\sigma^2+21736\sigma\tau+134160\sigma
+4576\tau^2+56160\tau+75285
\big).}
\tag{1.3}
\]

The interpolation uses six exact normalized Hamiltonian homotopies and an
independent seventh point. Its slow replay is
[`derive_degree_seven_marked_root_shear.py`](../scripts/derive_degree_seven_marked_root_shear.py).
The sparse constructor verifies the canonical relation `{S,T}=1` directly.
The classical symbol bounds are

\[
\begin{array}{c|cc}
&\deg_Z&\deg_B\\ \hline
S&7&43\\
T&6&39.
\end{array}
\tag{1.4}
\]

## 2. Relative restricted complex

The canonical relative complex remains

\[
C^0\longrightarrow C^1\longrightarrow C^2,
\tag{2.1}
\]

with the differentials and coherent obstruction module defined in the parent
family-search note. Intersecting it with the inherited root-boundary PBW
lattice gives the complete order-three and order-five correction summands

\[
\begin{array}{c|ccc|c}
&\deg_Z\le&\deg_B\le&\nu&\#\\ \hline
S_2&5&39&4&72\\
T_2&4&35&5&54\\
S_4&3&35&10&38\\
T_4&2&31&11&25.
\end{array}
\tag{2.2}
\]

The exact sparse presentation is
[`verify_degree_seven_relative_quantization_obstruction.py`](../scripts/verify_degree_seven_relative_quantization_obstruction.py).
At the rational point `(sigma,tau)=(1,0)` it gives

\[
\operatorname{rank}d_3=118,
\qquad \dim\ker d_3=8,
\tag{2.3}
\]

and

\[
\operatorname{rank}D_5=61,
\qquad \operatorname{rank}M_5=68,
\qquad \operatorname{rank}[M_5\mid O_0]=69.
\tag{2.4}
\]

The output support has dimension 178. A fixed 68-column pivot chart proves
exactly over `Q(sigma,tau)` that all remaining 39 strong columns are
dependent. Thus

\[
\boxed{\operatorname{rank}_{\mathbb Q(\sigma,\tau)}M_5=68.}
\tag{2.5}
\]

The exact rank-only replay is
[`reconstruct_degree_seven_order_five_zero_scheme.py`](../scripts/reconstruct_degree_seven_order_five_zero_scheme.py)
with `--rank-only`. Equation (2.4) at one rational point plus (2.5) proves
that the order-five section is nonzero on a nonempty characteristic-zero
open.

## 3. Multi-prime section-zero screen

The full-plane scanner specializes the complete presentation at every point
of `GF(p)^2`. Rank drop of `M5` and consistency after adjoining `O0` are
recorded separately. The results are:

| `p` | points scanned | exceptional signatures | section-zero points |
|---:|---:|---:|---:|
| 17 | 289 | 0 | 0 |
| 19 | 361 | 3 | 2 |
| 23 | 529 | 0 | 0 |
| 29 | 841 | 3 | 1 |
| 31 | 961 | 1 | 0 |
| 37 | 1369 | 2 | 2 |
| 41 | 1681 | 3 | 3 |
| 43 | 1849 | 3 | 1 |
| 47 | 2209 | 1 | 1 |
| 53 | 2809 | 4 | 3 |
| 59 | 3481 | 3 | 1 |
| 61 | 3721 | 3 | 2 |

There are sixteen recorded section-zero points. Every one has signature

\[
(118,8,61,68,68),
\tag{3.1}
\]

so the strong matrix does not drop rank there; the constant section becomes
dependent. The splitting counts

\[
0,2,0,1,0,2,3,1,1,3,1,2
\tag{3.2}
\]

are exactly the rational-point counts of the characteristic-zero scheme
reconstructed below. Low-degree lattice fits through the first primes fail
on later primes; the final reconstruction uses complete saturated
pivot-chart residuals, not a fit to these counts.

The scanner is
[`search_degree_seven_order_five_fitting_locus.py`](../scripts/search_degree_seven_order_five_fitting_locus.py).

## 4. Genuine nonlinear order-five gate

Strong consistency is only a necessary condition because the linear span
allows the quadratic lower-lift coefficients to vary independently. For each
of the sixteen section-zero points, project the genuine defect

\[
O_5(z)=O_0+\sum_{i=1}^8z_iO_i
+\sum_{1\le i\le j\le8}z_iz_jO_{ij}
\tag{4.1}
\]

modulo all 61 current corrections. This gives 117 quadrics in eight
variables. Singular reduces them, at every recorded point, to a four-element
Groebner basis of a dimension-four nonempty scheme. Hence all sixteen modular
points pass the nonlinear order-five gate.

The exact modular replay is
[`screen_degree_seven_order_five_survivors.py`](../scripts/screen_degree_seven_order_five_survivors.py).
This modular calculation is a discovery gate. The next section performs the
independent characteristic-zero reconstruction and exact nonlinear check.

## 5. Characteristic-zero section-zero scheme

On the fixed 68-column pivot chart, four residual rational functions suffice
to cut out the section-zero scheme. Their common denominator has total degree
32 and their numerator degree bounds are `34,36,33,33`. Nested univariate
interpolation on `76 x 44` grids reconstructs these functions independently
over each build prime. Saturating their numerator ideal by the common
denominator always gives a four-element standard basis of vector-space
dimension eight.

Coefficientwise CRT over the fifteen primes

\[
1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,
1069,1087,1091,1093,1097
\tag{5.1}
\]

has 151 bits. Balanced rational reconstruction is unchanged after dropping
the last build prime, and reduction agrees coefficientwise with the unused
prime `1103`. Exact standard-basis computation over `Q` gives dimension zero
and length eight. In lexicographic shape the component is

\[
\boxed{
\begin{aligned}
p(\sigma)={}&1687500\sigma^8+41047500\sigma^7
+462666750\sigma^6+3259667250\sigma^5\\
&+15781954748\sigma^4+53969799492\sigma^3
+126253770468\sigma^2\\
&+183369004011\sigma+119142437697=0,
\end{aligned}}
\tag{5.2}
\]

and

\[
\boxed{
\begin{aligned}
383180852409815403888\tau
&+123670522877062500\sigma^7
+2282786985300315000\sigma^6\\
&+19181170262239043250\sigma^5
+100324562410939815000\sigma^4\\
&+337020658576478999692\sigma^3
+599309543339487802104\sigma^2\\
&+840178458581066844732\sigma
+269809615785764223981=0.
\end{aligned}}
\tag{5.3}
\]

The polynomial `p` is irreducible over `Q`. Thus (5.2)--(5.3) is one closed
point of residue degree eight, not eight rational components. Its reductions
have exactly all sixteen full-plane section-zero points from Section 3,
including `(15,42)` over `GF(53)`, which lies outside the original pivot
minor after reduction. No modular survivor is discarded. Direct evaluation
of the defining `68 x 68` pivot minor over the octic residue field has rank
68, so the characteristic-zero component itself is genuinely inside the
chart; saturation did not create or lose it.

Let `K=Q[sigma]/(p)`, with `tau` given by (5.3). Exact recomputation gives

\[
(\operatorname{rank}d_3,\dim\ker d_3,
\operatorname{rank}D_5,\operatorname{rank}M_5,
\operatorname{rank}[M_5\mid O_0])=(118,8,61,68,68).
\tag{5.4}
\]

The 117 genuine quadratic Kuranishi equations over `K` have a four-element
standard basis: three affine-linear relations and one quadratic relation.
The quadratic discriminant is exactly zero. Consequently the order-five
lift scheme is a nonreduced double structure supported on

\[
\boxed{(\mathcal L_5)_{\mathrm{red}}\simeq\mathbb A^4_K,}
\tag{5.5}
\]

with free lower-lift coordinates `z0,z2,z4,z6`. This perfect-square check is
exact; it explains the modular dimension-four fibres without inventing a
spurious degree-sixteen extension.

## 6. Complete order-seven closure

The root-at-infinity valuation remains

\[
\nu(X,Q,Z)=(1,-1,-2),
\tag{6.1}
\]

and all correction weights through the terminal order are

| term | `S` | `T` | `S2` | `T2` | `S4` | `T4` | `S6` | `T6` |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `nu` | -2 | -1 | 4 | 5 | 10 | 11 | 16 | 17 |

The complete order-six correction summands have sizes

\[
\#S_6=14,\qquad \#T_6=6,
\qquad \operatorname{rank}D_7=20.
\tag{6.2}
\]

At order five, `D5` has a two-dimensional kernel. The verifier first solves
the current correction as a quadratic function on (5.5), using 15 Newton
nodes and the independent point `(3,1,2,4)`, and then retains both kernel
coordinates. The complete order-seven defect is therefore a cubic in six
variables: four lower-lift coordinates and two current-correction
coordinates. Its 124-dimensional output, modulo the 20 columns in (6.2),
gives 104 projected equations. Exact values at all 84 total-degree-at-most-
three Newton nodes reconstruct them over `K`; `(4,1,2,3,2,1)` is an unused
holdout.

Singular proves exactly

\[
\boxed{I_7=(1).}
\tag{6.3}
\]

Thus every reduced order-five lift is obstructed at `hbar^7`. The
Ore-localized quantization problem has been solved exactly through its first
terminal inconsistency. This also closes the doubled structure in (5.5): if
`N` is its square-zero nilradical, the image of the order-seven ideal in
`A/N` is `(1)` by (6.3). Hence the order-seven ideal contains an element
`1+n` with `n` nilpotent, and `1+n` is a unit. Therefore the order-seven
ideal is already `(1)` on the full nonreduced order-five scheme, not only on
its radical.

## 7. Terminal audit

| requested stage | result |
|---|---|
| reconstruct the component over `Q` | completed by (5.2)--(5.3) |
| solve Ore-localized quantization exactly | inconsistent at `hbar^7` by (6.3) |
| compute root-at-infinity valuations | every allowed summand through `S6,T6` is recorded in (6.1)--(6.2) |
| test conductor gluing | not applicable: there is no `hbar^7` local quantization to glue |
| certify Weyl relations and nonsurjectivity | not applicable: no quantized homomorphism survives |

The conditional downstream stages cannot be certified vacuously. Equation
(6.3) is precisely the certificate that there is no Weyl tuple in this
restricted filtration on which those tests could operate.

## 8. Reproduction

The completed certificates are produced by

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_seven_relative_quantization_obstruction.py \
  --output artifacts/generated-results/degree_seven_relative_quantization_obstruction.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/search_degree_seven_order_five_fitting_locus.py \
  --prime 19 --jobs 8 \
  --output artifacts/generated-results/degree_seven_order_five_scan_gf19.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/screen_degree_seven_order_five_survivors.py \
  artifacts/generated-results/degree_seven_order_five_scan_gf*.json \
  --output artifacts/generated-results/degree_seven_order_five_nonlinear_screen.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/reconstruct_degree_seven_order_five_zero_scheme.py \
  --rank-only \
  --output artifacts/generated-results/degree_seven_order_five_exact_strong_rank.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/interpolate_degree_seven_order_five_chart.py \
  --prime 1097 --jobs 8 \
  --output artifacts/generated-results/degree_seven_order_five_chart_gf1097.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/reconstruct_degree_seven_order_five_rational_chart.py \
  --holdout-prime 1103 \
  --output artifacts/generated-results/degree_seven_order_five_rational_chart.json

PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_seven_order_five_survivor.py \
  --output artifacts/generated-results/degree_seven_order_five_survivor.json
```

Repeat the interpolation command for the fifteen build primes in (5.1) and
for holdout prime `1103`. The reconstruction command requires all of those
artifacts. The final exact nonlinear and order-seven command requires
Singular and takes several minutes.

The slow independent classical shear replay is

```bash
.venv/bin/python scripts/derive_degree_seven_marked_root_shear.py \
  --output artifacts/generated-results/degree_seven_marked_root_shear.json
```
