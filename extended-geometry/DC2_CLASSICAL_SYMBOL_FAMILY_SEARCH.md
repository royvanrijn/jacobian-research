# Search `(DC_2)` across classical-symbol families

> **Status and scope.** This note branches away from the obstructed
> degree-five symbol. It constructs four initial two-parameter classical
> polynomial symplectic families, attaches the relative
> `gauge -> corrections -> defects` complex, and screens the known
> order-five and order-seven restricted cocycles. Further marked-root slices
> in degrees seven and eight are then treated by the same pipeline. The
> characteristic-zero survivors have residue degrees four, eight, and twelve.
> All three have genuine order-five lifts, but their complete inherited
> order-seven equations generate the unit ideal. Consequently no component reaches conductor
> gluing, full Weyl relations, or nonsurjectivity. All quantum conclusions are
> internal to the declared parity-preserving, root-weight-homogeneous
> filtration and do not settle `DC_2`.

The calculation implements the change of question

\[
\text{“can one fixed symbol be repaired?”}
\quad\leadsto\quad
\text{“where does the restricted obstruction section vanish?”}
\]

without enlarging support over the old degree-five symbol.

## 1. Six two-parameter family rows

The census contains one genuine noninjective search branch, one transported
noninjective control, and two exact-quantization automorphism controls.

| key | construction | parameters | classical status | PBW screening status |
|---|---|---|---|---|
| `MR6` | normalized degree-six weighted marked-root incidence | `(sigma,tau)` | noninjective exact polynomial symplectic family | generically obstructed at order five; one quartic closed component reaches order seven |
| `MR7` | monic degree-seven weighted marked-root incidence | `(sigma,tau)` | noninjective exact polynomial symplectic family | one octic closed component has a doubled order-five lift space and is obstructed at order seven |
| `MR8` | monic degree-eight weighted marked-root incidence | `(sigma,tau)` | noninjective exact polynomial symplectic family | one degree-twelve closed component has a two-square thickening of `A^5` and is obstructed at order seven |
| `NF3` | normalized-factorization `c=-9` completion followed by two target symplectic shears | `(rho,eta)` | noninjective, but one transported target orbit | known order-five value remains `-49`; vanishing locus empty |
| `HS3` | cubic Hamiltonian suspension | `(lambda,mu)` | polynomial symplectic automorphisms | order-five and order-seven cocycles vanish identically |
| `RF3` | rank-one fibre shear over central `R`, with its connection coordinate | `(lambda,mu)` | polynomial symplectic automorphisms | order-five and order-seven cocycles vanish identically |

The `MR6`, `MR7`, and `MR8` rows are the nonautomorphism moduli branches. The other
rows are necessary controls: they distinguish a vanishing obstruction caused
by exact quantizability from a vanishing obstruction on a potentially
non-surjective classical symbol. The full higher-degree calculations are kept
in the [degree-seven marked-root note](DC2_DEGREE_SEVEN_MARKED_ROOT_SEARCH.md)
and [degree-eight marked-root note](DC2_DEGREE_EIGHT_MARKED_ROOT_SEARCH.md).

### 1.1 Marked-root family `MR6`

Fix `kappa=-9` and put

\[
\begin{aligned}
H_{\sigma,\tau}(W)=W^2(W-1)\bigg(&
\sigma W^3+\tau W^2
+\left(-\frac52-3\sigma-2\tau\right)W\\
&+\frac32+2\sigma+\tau\bigg).
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

The exact-degree chart is `sigma != 0`. The rank-two completion and its
quadratic shear are those of the
[degree-six relative obstruction calculation](DEGREE_SIX_RELATIVE_QUANTIZATION_OBSTRUCTION.md).
This is the marked-root incidence branch.

### 1.2 Normalized-factorization orbit `NF3`

Let `(R,T,D,S)` be the exact `c=-9` completion obtained from normalized
linear-quadratic factorization. Its canonical brackets are

\[
\{D,R\}=\{S,T\}=1
\]

and all mixed brackets vanish. Define

\[
\Phi_{\rho,\eta}
=\left(R,T,D+\rho R,S+\eta T\right).
\tag{1.3}
\]

The two added target shears are linear exact symplectic automorphisms. Hence
`Phi_(rho,eta)` is an explicit two-parameter family of noninjective
polynomial symplectic maps, but it is a single transported orbit rather than
new normalized-factorization moduli. Exact Weyl-linear target shears cannot
remove the transported `c=-9` restricted obstruction. Its order-five value
is still `-49` over `Q[rho,eta]`.

### 1.3 Hamiltonian suspension `HS3`

In canonical coordinates `(q_1,q_2,p_1,p_2)`, let

\[
K_{\lambda,\mu}
=\frac\lambda3q_1^3+\mu q_1q_2^2.
\tag{1.4}
\]

The time-one triangular map is

\[
\begin{aligned}
Q_1&=q_1,&
Q_2&=q_2,\\
P_1&=p_1+\lambda q_1^2+\mu q_2^2,&
P_2&=p_2+2\mu q_1q_2.
\end{aligned}
\tag{1.5}
\]

Because `q_1,q_2` commute in the Weyl algebra, the same formulas define an
exact Weyl automorphism. All odd Moyal terms of orders at least three vanish
coefficientwise. Thus both known cocycles vanish on the whole parameter
plane, but the maps are surjective controls.

### 1.4 Central rank-one fibre family `RF3`

Use canonical pairs `(v,P)` and `(R,U)`, with `R` central for the first
fibre. Set

\[
F_{\lambda,\mu}(v,R)
=\frac\lambda3Rv^3+\frac\mu2R^2v^2.
\tag{1.6}
\]

The relative fibre map and its Hamiltonian connection are

\[
\begin{aligned}
v'&=v,&R'&=R,\\
P'&=P+\lambda Rv^2+\mu R^2v,&
U'&=U+\frac\lambda3v^3+\mu Rv^2.
\end{aligned}
\tag{1.7}
\]

The `U` correction is forced by `partial_R F`; omitting it would destroy the
mixed canonical relations. Again `v,R` commute, so (1.7) is already an exact
Weyl automorphism. Its order-five and order-seven vanishing loci are the
entire parameter plane and are excluded from the nonsurjective search.

## 2. Relative restricted complexes

For every family, work in canonical fibre coordinates `(S,T)` over the
central parameter `R`. The common relative complex is

\[
C^0\xrightarrow{d_0}C^1\xrightarrow{d_1}C^2,
\tag{2.1}
\]

with

\[
d_0(h)=(-h_T,h_S,-h_R)
\tag{2.2}
\]

and

\[
d_1(s,t,a)
=\left(s_S+t_T,\ s_R-a_T,\ t_R+a_S\right).
\tag{2.3}
\]

The defects satisfy

\[
\partial_RF-\partial_SG-\partial_TH=0,
\tag{2.4}
\]

so `C^2` is the closed-defect module. In canonical polynomial degrees
`4 -> 3 -> 2`, the integral presentation has

\[
35\xrightarrow[\operatorname{rank}34]{d_0}
60\xrightarrow[\operatorname{rank}26]{d_1}26,
\qquad
(\dim H^0,\dim H^1,\dim H^2)=(1,0,0).
\tag{2.5}
\]

Consequently the coherent canonical obstruction module

\[
E^{\mathrm{can}}=\operatorname{coker}d_1
\tag{2.6}
\]

is zero and `Fitt_0(E^can)=(1)` after base change to each parameter ring.
For `HS3` and `RF3`, the global polynomial symplectic automorphism transports
this complex directly. For `NF3`, (2.5) is only the canonical local complex;
the root-boundary PBW lattice retains the transported order-five class.

The family-specific PBW complex is the pullback of (2.1) intersected with
the declared boundary-weight and filtration lattice. After Hamiltonian gauge
is fixed, its order-`m` differential is

\[
d_m(s_m,t_m)=\{s_m,T\}+\{S,t_m\}.
\tag{2.7}
\]

Let

\[
E_m=\operatorname{coker}d_m.
\tag{2.8}
\]

Lower-lift dependence is retained by adjoining every coefficient of the
Kuranishi section to `d_m`. If `M_m` is this strong presentation, then

\[
E_m^{\mathrm{str}}=\operatorname{coker}M_m,
\qquad
\mathcal P_m=\ker M_m^\vee.
\tag{2.9}
\]

The rank-drop Fitting locus of `E_m^str` and the zero scheme of the
obstruction section are related but not identical. The calculation below
keeps that distinction explicit.

## 3. Order-five Fitting screen on `MR6`

At order three, the complete correction module has `83` columns, rank `77`,
and a six-dimensional affine lift torsor. At order five:

\[
\operatorname{rank}D_5=34,
\qquad
\operatorname{rank}M_5=41,
\qquad
\operatorname{rank}[M_5\mid O_0]=42
\tag{3.1}
\]

over `Q(sigma,tau)`. Thus a nonempty open is obstructed.

Full scans over seven primes reconstruct the following primitive quartic as
the candidate rank-`40` Fitting divisor:

\[
\begin{aligned}
C(\sigma,\tau)={}&
2563\sigma^4+3954\sigma^3\tau+2319\sigma^2\tau^2
+608\sigma\tau^3+60\tau^4\\
&+7240\sigma^3+7200\sigma^2\tau+2280\sigma\tau^2
+200\tau^3\\
&+6970\sigma^2+3400\sigma\tau+250\tau^2+2250\sigma.
\end{aligned}
\tag{3.2}
\]

This divisor is a multi-prime reconstruction, not yet a complete
characteristic-zero equality with the ideal of all `41`-minors. It is not
used as a theorem below. What is verified exactly is the closed component
inside it.

The fibrewise consistency condition

\[
\operatorname{rank}[M_5\mid O_0]
=\operatorname{rank}M_5
\tag{3.3}
\]

has the splitting profile

\[
2,0,4,0,2,0,0,4,4,4
\tag{3.4}
\]

at

\[
p=17,19,23,29,31,37,41,139,167,197.
\]

The reconstructed characteristic-zero component is

\[
\boxed{
4\sigma^4+66\sigma^3+561\sigma^2+1260\sigma+900=0,
}
\tag{3.5}
\]

\[
\boxed{
4260\tau+52\sigma^3+498\sigma^2+12693\sigma+13050=0.
}
\tag{3.6}
\]

The quartic (3.5) is irreducible over `Q`. Substitution proves exactly that
(3.5)--(3.6) lies on (3.2). Let

\[
K=\mathbb Q[\alpha]/
(4\alpha^4+66\alpha^3+561\alpha^2+1260\alpha+900).
\tag{3.7}
\]

Over `K`, with `tau` given by (3.6), the exact ranks are

\[
\operatorname{rank}d_3=77,
\quad \dim\ker d_3=6,
\quad \operatorname{rank}D_5=34,
\quad \operatorname{rank}M_5
=\operatorname{rank}[M_5\mid O_0]=40.
\tag{3.8}
\]

This proves that the strong order-five cocycle disappears on the exact
quartic component. It does not yet prove that a nonlinear lower lift exists.

## 4. Genuine order-five Kuranishi scheme

Project the complete order-five defect modulo the `34` current-correction
columns. The result is `76` quadratic equations in the six order-three lift
parameters. Exact standard-basis computation over `K` proves that their
ideal is generated by two independent affine-linear forms. The verifier
checks both containments between the input ideal and these two generators.
Therefore

\[
\boxed{\mathcal L_5\simeq\mathbb A^4_K.}
\tag{4.1}
\]

This is a genuine lift scheme, not the linear-span relaxation used for
multi-prime discovery. The two exact generators have large quartic-field
coefficients and are retained in the verifier rather than copied into this
note.

## 5. Complete inherited order-seven calculation

The root-at-infinity valuation is

\[
\nu(X,Q,Z)=(1,-1,-2).
\tag{5.1}
\]

The induced weights are

\[
\begin{array}{c|rrrrrrrr}
&S&T&S_2&T_2&S_4&T_4&S_6&T_6\\ \hline
\nu&-2&-1&4&5&10&11&16&\text{absent}.
\end{array}
\tag{5.2}
\]

The complete order-six correction space in the inherited filtration is

\[
S_6\in\operatorname{span}
\{X^{16},X^{17}Q,X^{18}Q^2,X^{19}Q^3,X^{20}Q^4\},
\qquad T_6=0.
\tag{5.3}
\]

It has rank five. For every point of (4.1), the unique order-five current
correction is solved exactly and substituted into the order-seven defect

\[
\begin{aligned}
O_7={}&\{S_2,T_4\}+\{S_4,T_2\}\\
&+\frac1{24}\left(
\Pi^3(S_4,T)+\Pi^3(S,T_4)+\Pi^3(S_2,T_2)
\right)\\
&+\frac1{1920}\left(
\Pi^5(S_2,T)+\Pi^5(S,T_2)
\right)
+\frac1{322560}\Pi^7(S,T).
\end{aligned}
\tag{5.4}
\]

After quotienting the five columns (5.3), there are `60` projected equations
in the four affine parameters of (4.1). They have total degree at most three.
Exact values at the `35` Newton interpolation nodes reconstruct every
coefficient over `K`; an additional point `(4,1,2,3)` independently verifies
the reconstruction. Singular then proves

\[
\boxed{I_7=(1)\subset K[u_0,u_2,u_4,u_5].}
\tag{5.5}
\]

Hence every point of the complete affine four-space (4.1) is obstructed at
order seven.

## 6. Terminal audit

The requested post-survival stages are conditional. Their status here is:

| stage | result |
|---|---|
| reconstruct component over `Q` | completed by (3.5)--(3.6) |
| solve Ore-localized quantization exactly | completed through the terminal order-seven inconsistency (5.5) |
| compute root-at-infinity valuations | filtration weights (5.1)--(5.3) computed; no all-order correction exists to value further |
| test conductor gluing | not reached |
| certify Weyl relations | not reached |
| certify nonsurjectivity | not reached |

It would be incorrect to perform conductor or nonsurjectivity tests on this
component after (5.5): there is no Weyl tuple in the declared filtration to
glue or test.

The outcome of this census is therefore:

\[
\boxed{
\begin{array}{l}
\text{the two tame families have identically vanishing classes but are automorphisms;}\\
\text{the normalized-factorization orbit is uniformly obstructed at order five;}\\
\text{all three genuine marked-root survivors are obstructed at order seven.}
\end{array}}
\tag{6.1}
\]

There is no order-seven survivor in these six rows. The degree-seven
marked-root branch reconstructs an irreducible octic section-zero point whose
genuine order-five scheme is a doubled affine four-space; its 104 complete
order-seven equations in six lift coordinates again generate `(1)`. The
degree-eight branch reconstructs one irreducible degree-twelve point, whose
reduced lift is affine five-space and whose six independent order-seven cubic
generators generate `(1)`. The next search should therefore change the
noninjective classical construction, not enlarge the old quintic or any
now-closed marked-root component.

For comparison rather than support enlargement, the
[marked-root degree ladder](DC2_MARKED_ROOT_DEGREE_LADDER.md) now isolates the
closed PBW support formulas and the exact rank pattern through degree eight.
Its degree-eight row is now closed through the same characteristic-zero
reconstruction and terminal order-seven gates used here.

## 7. Reproduction

Construct the four families, verify their common relative complex, and check
the exact-quantization controls with

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_dc2_classical_symbol_families.py \
  --output artifacts/generated-results/dc2_classical_symbol_families.json
```

The complete modular scans used for discovery are produced by

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/search_degree_six_order_five_fitting_locus.py \
  --prime 23 --jobs 8 \
  --output \
  artifacts/generated-results/degree_six_order_five_fitting_GF23.json
```

Replace `23` by the other listed good primes. After reconstruction, the
candidate-divisor option restricts larger-prime checks to (3.2):

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/search_degree_six_order_five_fitting_locus.py \
  --prime 197 --jobs 8 --candidate-divisor \
  --output \
  artifacts/generated-results/degree_six_order_five_fitting_GF197.json
```

The characteristic-zero component, its exact order-five lift scheme, and the
order-seven unit ideal are checked by

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_six_order_five_survivor.py \
  --output \
  artifacts/generated-results/degree_six_order_five_survivor.json
```

The last command requires Singular. The generated JSON records the prime
splitting profile, exact ranks, interpolation data, valuation weights, and
SHA-256 hashes of both dynamically generated Singular programs.
