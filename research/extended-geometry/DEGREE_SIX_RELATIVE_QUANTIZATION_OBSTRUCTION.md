# A parameter-uniform degree-six rank-two quantization obstruction

> **Status and scope.** On a nonempty Zariski open of the two-parameter
> degree-six weighted-symbol family with \(\kappa=-9\), the complete
> parity-preserving, root-boundary-weight-homogeneous quantization is
> obstructed at \(\hbar^5\). The calculation retains the complete
> six-dimensional \(\hbar^3\) lift torsor and all allowed \(\hbar^5\)
> corrections in the declared filtration. It is an exact
> characteristic-zero function-field theorem. It does not exclude
> corrections of other boundary weights, odd powers of \(\hbar\), wider
> filtrations, or different polarizations, and therefore does not prove
> anything about `DC_2`.

This is the first output of the
[boundary-selected classical-symbol census](DC2_CLASSICAL_SYMBOL_CENSUS.md).
It varies the classical degree-six symbol rather than enlarging correction
support over the former quintic.

## 1. The classical family

Fix \(\kappa=-9\), and let

\[
\begin{aligned}
H_{\sigma,\tau}(W)
=W^2(W-1)\bigg(
&\sigma W^3+\tau W^2
 \left(-\frac52-3\sigma-2\tau\right)W\\
&+\frac32+2\sigma+\tau
\bigg).
\end{aligned}
\tag{D6Q1}
\]

Then

\[
H(0)=H'(0)=H(1)=0,\qquad
H'(1)=-1,\qquad H''(1)=-9.
\tag{D6Q2}
\]

On \(\sigma\ne0\), this has exact degree six. The polynomial rank-two
completion uses the classical shear

\[
W=Z+s_2Q^2,
\qquad
s_2=
\frac{
6024\sigma^2+5016\sigma\tau+11088\sigma
+1056\tau^2+4752\tau-16929
}{2156}.
\tag{D6Q3}
\]

The resulting pair \(S_{\sigma,\tau},T_{\sigma,\tau}\) lies in
\(\mathbb Q[\sigma,\tau][X,Q,Z]\) and satisfies

\[
\{S,T\}=1
\tag{D6Q4}
\]

for the Ore-symbol bracket

\[
\{f,g\}
=f_Z\delta(g)-\delta(f)g_Z,\qquad
\delta=3X^2\partial_X+(2-6XQ)\partial_Q.
\tag{D6Q5}
\]

Its differential orders and Bernstein degrees are

\[
\begin{array}{c|cc}
&\deg_Z&\deg_B\\ \hline
S&6&36\\
T&5&32,
\end{array}
\qquad
\deg_B(X,Q,Z)=(1,1,3).
\tag{D6Q6}
\]

The exact sparse constructor is shared with the degree-five regression; it
builds \(H=W^2(W-1)P(W)\) coefficientwise and verifies (D6Q4) over any exact
coefficient field.

## 2. The complete restricted correction spaces

Give the root-at-infinity chart the valuation weight

\[
\nu(X,Q,Z)=(1,-1,-2).
\tag{D6Q7}
\]

Then \(\nu(S,T)=(-2,-1)\), while every Poisson contraction raises weight by
three. Use the parity-preserving ansatz

\[
S_\hbar=S+\hbar^2S_2+\hbar^4S_4,\qquad
T_\hbar=T+\hbar^2T_2+\hbar^4T_4.
\tag{D6Q8}
\]

The inherited order-lowering filtration and weight equations give the
complete relevant summands

\[
\begin{array}{c|ccc|c}
&\deg_Z\le&\deg_B\le&\nu&\#\text{ monomials}\\ \hline
S_2&4&32&4&49\\
T_2&3&28&5&34\\
S_4&2&28&10&22\\
T_4&1&24&11&12.
\end{array}
\tag{D6Q9}
\]

These are full weight summands, not sparse supports selected after solving.
Corrections of other weights are outside the theorem.

## 3. Uniform order-three lift

Let

\[
K=\mathbb Q(\sigma,\tau).
\]

The first equation is

\[
d_3(S_2,T_2)
=\{S_2,T\}+\{S,T_2\}
=-\frac1{24}\Pi^3(S,T),
\tag{D6Q10}
\]

where

\[
\Pi=\partial_Z\otimes\delta-\delta\otimes\partial_Z.
\]

On the \(49+34=83\) columns in (D6Q9), exact sparse row reduction over
\(K\) gives

\[
\operatorname{rank}d_3=77,
\qquad
\dim_K\ker d_3=6,
\tag{D6Q11}
\]

and (D6Q10) is soluble. Thus the complete restricted lower-lift scheme is
an affine six-space over the function field:

\[
\mathcal L_3\simeq\mathbb A^6_K.
\tag{D6Q12}
\]

Fix one solution and a kernel basis \(e_1,\ldots,e_6\), and write a general
lower lift as

\[
(S_2,T_2)=(S_2^0,T_2^0)+\sum_{i=1}^6z_ie_i.
\tag{D6Q13}
\]

## 4. The strong-cocycle and Fitting presentation

The order-five defect is

\[
\begin{aligned}
\mathcal O_5={}&
\{S_2,T_2\}
+\frac1{24}\Pi^3(S_2,T)
+\frac1{24}\Pi^3(S,T_2)\\
&+\frac1{1920}\Pi^5(S,T).
\end{aligned}
\tag{D6Q14}
\]

Expand (D6Q14) using (D6Q13):

\[
\mathcal O_5(z)
=O_0+\sum_i z_iO_i+\sum_{i\le j}z_iz_jO_{ij}.
\tag{D6Q15}
\]

There are six linear, six diagonal-quadratic, and fifteen cross-quadratic
coefficient columns, for a total of \(27\). Let \(D_5\) be the \(34\)-column
current-correction matrix from \(S_4,T_4\), and form

\[
M_5=
\left[
D_5\mid
(O_i)_i\mid
(O_{ij})_{i\le j}
\right].
\tag{D6Q16}
\]

All columns lie in the \(110\)-dimensional output support computed by the
exact expansion. Over \(K\),

\[
\operatorname{rank}D_5=34,\qquad
\operatorname{rank}M_5=41,
\qquad
\operatorname{rank}[M_5\mid O_0]=42.
\tag{D6Q17}
\]

The coherent consistency module and strong dual-cocycle module are

\[
E^{\mathrm{str}}_5=\operatorname{coker}M_5,
\qquad
\mathcal P_5=\ker M_5^\vee.
\tag{D6Q18}
\]

Generically \(E^{\mathrm{str}}_5\) and \(\mathcal P_5\) have rank \(69\).
In Fitting language, the \(41\)-minor locus of \(M_5\) is nonempty, every
\(42\)-minor of \(M_5\) vanishes over \(K\), and a \(42\)-minor becomes
nonzero after adjoining \(O_0\). Equivalently, there is

\[
\Lambda\in\mathcal P_5\otimes K
\tag{D6Q19}
\]

such that

\[
\Lambda(O_i)=\Lambda(O_{ij})=0,\qquad
\Lambda(D_5)=0,\qquad
\Lambda(O_0)\ne0.
\tag{D6Q20}
\]

Thus \(\Lambda(\mathcal O_5(z))\) is independent of every lower-lift
parameter and is nonzero.

At the boundary-clean rational point

\[
(\sigma,\tau)=(1,0),
\tag{D6Q21}
\]

the checker constructs an explicit \(30\)-term rational functional,
normalizes it by

\[
\Lambda(O_0)=1,
\tag{D6Q22}
\]

and verifies all \(61\) annihilation equations exactly over \(\mathbb Q\).
The complete coefficients are stored in the generated JSON certificate
rather than displayed here.

## 5. Parameter-uniform obstruction theorem

### Theorem

There is a nonempty Zariski open

\[
U\subset\mathbb A^2_{\sigma,\tau}
\tag{D6Q23}
\]

such that, for every geometric point of \(U\), no pair (D6Q8) in the
complete spaces (D6Q9) satisfies

\[
[S_\hbar,T_\hbar]=\hbar\pmod{\hbar^7}.
\tag{D6Q24}
\]

More precisely, every point of the six-dimensional order-three lift torsor
is obstructed at order five by a strong restricted dual cocycle.

### Proof

Equations (D6Q11) and (D6Q17) are exact ranks over
\(\mathbb Q(\sigma,\tau)\). Clear the finitely many pivot denominators and
retain the nonzero \(77\)-, \(41\)-, and augmented \(42\)-minors. Their
common nonvanishing defines a Zariski open \(U\). On \(U\), the order-three
solution scheme has the form (D6Q12), the strong module has constant rank,
and (D6Q20) specializes. Hence every evaluation of (D6Q15) has the same
nonzero \(\Lambda\)-value, while all current corrections have value zero.
The order-five equation is therefore inconsistent. The rational point
(D6Q21) has the same rank profile, proving that \(U\) is nonempty. QED

The point (D6Q21) is also the explicit Hessian-clean, boundary-clean witness
for the degree-six stable-moduli theorem. After shrinking \(U\) by that
marked open, the theorem applies to a two-dimensional family of pairwise
stably polynomially left-right inequivalent classical symbols.

This meets the census success criterion of a parameter-uniform obstruction
theorem. It does not imply that these symbols have no unrestricted
quantization: odd corrections, mixed boundary weights, enlarged
Bernstein/differential-order bounds, and different Weyl polarizations remain
outside the statement.

## Reproduction

Run:

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_six_relative_quantization_obstruction.py \
  --output \
  artifacts/generated-results/degree_six_relative_quantization_obstruction.json
```

The command performs the characteristic-zero function-field reductions,
constructs the rational strong cocycle, and repeats the ranks at
`31991`, `32003`, and `65521`. The good-prime calculations are regressions;
the theorem is proved by the exact function-field calculation.
