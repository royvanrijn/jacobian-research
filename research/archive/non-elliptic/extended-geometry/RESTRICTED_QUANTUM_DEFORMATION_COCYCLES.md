# Restricted quantum deformation cocycles

## 1. Outcome

The residue calculations at orders \(\hbar^5\) and \(\hbar^7\) admit a
precise gauge-invariant interpretation, provided the word *restricted* is
taken literally.

At a fixed correction order let

\[
 C^0\xrightarrow{d_0}C^1\xrightarrow{d_1}C^2
\tag{RDC1}
\]

be the filtered complex determined by the allowed Hamiltonians, correction
symbols, and defect symbols.  It is finite in every bounded ansatz; for the
all-pole Laurent space it is a filtered union, and the finite support of the
functional makes it factor through a finite quotient.

This is the Weyl instance of the
[unified gauge--correction--defect complex](UNIFIED_DEFORMATION_COMPLEX.md).
That note adds the classical-symbol base and the localization, saturation,
and conductor layers; the present note remains the canonical source for the
two exact strong-cocycle certificates.

The computed functionals are elements of

\[
 Z^2_{\mathrm{res},\vee}
 :=\ker\!\left(d_1^\vee:(C^2)^\vee\longrightarrow(C^1)^\vee\right).
\tag{RDC2}
\]

Thus they are dual \(2\)-cocycles, not merely row-reduction witnesses.  A
second, stronger condition makes their values independent of the chosen
lower-order lift.  It is naturally the kernel of a morphism of coherent
modules, and therefore includes Hamiltonian gauge invariance without
choosing a gauge slice.

For the existing degree-five calculations this gives exact
instantiations at orders five and seven.

1. At \((\kappa,\tau)=(0,1)\), the 16-term Laurent functional
   \(\Lambda_5\) is a strong restricted dual cocycle.  It annihilates every
   finite-\(X\)-Laurent order-five correction and is constant on an affine
   \(1075\)-dimensional superset of all finite-Laurent order-three lifts:
   \[
   \Lambda_5(O_5)
   =-\frac{47547660815739}{190658}\ne0.
   \tag{RDC3}
   \]
2. On the explicit \(K\)-rational fifth-order lift on the cubic exceptional
   component, where
   \[
   K=\mathbb Q[a]/(94a^3+335a^2+400a+160),
   \]
   coefficient extraction
   \[
   \Lambda_7=[X^{18}]
   \tag{RDC4}
   \]
   is a strong restricted dual cocycle on the complete
   \(20\)-dimensional fifth-order correction torsor.  Its value is
   \[
   \Lambda_7(O_7)=
   \frac{2189187}{83886080}
   \left(
   587583566a^2+1388701707a+831388850
   \right)\ne0.
   \tag{RDC5}
   \]
3. On the full reduced \(\mathbb A^{27}_K\) fifth-order lift component,
   eliminate the constant six-column pivot in the order-seven variation
   matrix.  The resulting Schur functional
   \(\Lambda^{\mathrm{rel}}_{18}\) is a global section of
   \(\mathcal P_7\), and its value is again the constant (RDC5).

The first class is global with respect to every finite Laurent lower lift at
the stated seed.  The raw coefficient extraction in the second item is
branchwise.  Its Schur correction in the third item globalizes it across
the complete reduced fifth-order lift component.

## 2. The restricted correction complex

Fix a classical canonical pair \(S,T\) with \(\{S,T\}=1\).  At one
correction order, let \(C^0\) be the declared filtered Hamiltonian space,
let \(C^1=U\oplus V\) be the two declared correction spaces, and let \(C^2\)
be the declared defect space.  Include the appropriate output projection
in the following formulas when \(C^2\) is a filtered quotient:

\[
\begin{aligned}
 d_0H&=(\{H,S\},\{H,T\}),\\
 d_1(u,v)&=\{u,T\}+\{S,v\}.
\end{aligned}
\tag{RDC6}
\]

Jacobi gives

\[
\begin{aligned}
d_1d_0(H)
&=\{\{H,S\},T\}+\{S,\{H,T\}\}\\
&=\{H,\{S,T\}\}=0.
\end{aligned}
\tag{RDC7}
\]

Consequently a functional \(\Lambda:C^2\to B\) is a dual cocycle exactly
when

\[
\Lambda\circ d_1=0.
\tag{RDC8}
\]

It then descends to

\[
\Lambda:H^2_{\mathrm{res}}
=\operatorname{coker}(d_1)\longrightarrow B.
\tag{RDC9}
\]

This is the first part of the desired invariantization.  It proves
independence of the correction at the current order.  By itself it does not
prove independence of lower-order choices, because the higher obstruction
\(O_n\) is a polynomial function of those choices.

## 3. Strong cocycles and lower-lift independence

Let \(B\) be a seed-parameter ring and let

\[
p:\mathcal L_{<n}\longrightarrow\operatorname{Spec}B
\]

be the affine scheme of lower-order lifts satisfying all equations before
order \(n\).  Put

\[
E_n=\operatorname{coker}(d_{1,n})
\]

and let

\[
\bar o_n\in
\Gamma(\mathcal L_{<n},p^*E_n)
\tag{RDC10}
\]

be the order-\(n\) defect modulo current-order corrections.  On an affine
chart, write \(R=p_*\mathcal O_{\mathcal L_{<n}}\).  Pairing with
\(\bar o_n\) defines

\[
\rho_n:E_n^\vee\longrightarrow R/B,\qquad
\Lambda\longmapsto
\langle\Lambda,\bar o_n\rangle\bmod B.
\tag{RDC11}
\]

Define the module of **strong restricted dual cocycles**

\[
\mathcal P_n:=\ker\rho_n.
\tag{RDC12}
\]

This definition has three immediate consequences.

* Since \(E_n^\vee=\ker d_{1,n}^\vee\), every section of \(\mathcal P_n\)
  is a dual cocycle of the restricted correction complex.
* Its evaluation on \(O_n\) is independent of every point of the
  lower-lift scheme, not only infinitesimal changes.
* Every allowed Hamiltonian gauge orbit lies inside the lower-lift scheme.
  Hence lower-lift independence implies gauge invariance without choosing
  a complement to \(\operatorname{im}d_0\).

In coordinates, expand

\[
O_n(z)=O_n^0+\sum_{\alpha\ne0}z^\alpha O_{n,\alpha}.
\tag{RDC13}
\]

Then (RDC12) is exactly the finite system

\[
\Lambda d_{1,n}=0,\qquad
\Lambda(O_{n,\alpha})=0\quad(\alpha\ne0).
\tag{RDC14}
\]

These are precisely the columns already used in the exact span
calculations: current corrections, followed by every nonconstant linear,
square, and cross-term coefficient of the lower-lift defect.  The former
rank-jump certificates were therefore computing fibers of
\(\mathcal P_n\), even though they were not named as cocycles.

There is a weaker gauge-only variant obtained by replacing \(B\) in
(RDC11) with the invariant subring \(R^{\mathrm{gauge}}\).  The present
certificates satisfy the stronger condition (RDC12), so no construction of
that invariant ring is needed at the certified fibers.

## 4. The all-pole \(\hbar^5\) cocycle

At \((\kappa,\tau)=(0,1)\), let \(C^1_5\) be the union of all allowed
finite-\(X\)-Laurent \(S_4,T_4\) corrections and let \(C^2_5\) contain the
finite defect support seen by the 16 monomials in the functional
\(\Lambda_5\).

The shift bounds prove that a negative-\(X\) correction monomial cannot
meet the support of \(\Lambda_5\).  Exact evaluation on every polynomial
column then proves

\[
\Lambda_5d_{1,5}=0
\quad\text{on the full finite-Laurent union.}
\tag{RDC15}
\]

For the lower lift, the finite cutoff argument produces a deliberately
larger affine scheme with 4065 ambient variables, rank 2990, and dimension
1075.  On its affine solution coordinates, the quadratic matrix of
\(\Lambda_5(O_5)\) is skew and every linear derivative is zero.  Since a
scalar quadratic polynomial uses the symmetric part of its matrix, this is
exactly

\[
\Lambda_5(O_{5,\alpha})=0
\quad\text{for every nonconstant coefficient.}
\tag{RDC16}
\]

Equations (RDC15)--(RDC16) put \(\Lambda_5\) in \(\mathcal P_5\).  Its
nonzero value (RDC3) is therefore invariant under:

* every allowed order-five finite-Laurent correction;
* every finite-Laurent solution of the order-three equation;
* every allowed Hamiltonian gauge change inside those spaces; and
* every choice introduced by the finite projected affine enlargement.

This locates the class in
\[
H^2(C^\bullet_X)^\vee,
\]
not in a connecting \(H^1\)-torsor: the obstruction remains nonzero after
localizing at \(X\).

The 15-term functional at \((\kappa,\tau)=(0,-3)\) satisfies the same
strong-cocycle conditions and gives value \(1\).  Its appearance only on
that fiber is the base-change failure of the coherent kernel
\(\mathcal P_5\), not a failure of gauge invariance.

## 5. The branchwise and global \(\hbar^7\) cocycles

On the cubic component, fix the sparse \(K\)-point of the genuine
fifth-order lift scheme displayed in
[Quantum residue obstruction](QUANTUM_RESIDUE_OBSTRUCTION.md#8-finite-laurent-bands).
The filtration permits no \(S_6,T_6\), so

\[
C^1_7=0,\qquad d_{1,7}=0.
\tag{RDC17}
\]

The remaining lower choice is the affine \(20\)-dimensional torsor of
solutions to the fifth-order correction equation.  The order-seven defect
is affine on this torsor.  Its 20 variation columns have rank six, and exact
calculation gives

\[
[X^{18}](O_{7,\alpha})=0
\quad\text{for all 20 variation columns.}
\tag{RDC18}
\]

Thus \([X^{18}]\in\mathcal P_7\) on this branch.  Evaluation gives (RDC5).
Nonvanishing is immediate because the factor

\[
587583566a^2+1388701707a+831388850
\]

has degree two and cannot vanish at a root of the irreducible cubic
defining \(K\).  As an independent arithmetic certificate, its norm is

\[
N_{K/\mathbb Q}
\left(587583566a^2+1388701707a+831388850\right)
=57870471643246100480\ne0.
\tag{RDC19}
\]

This proves that every one of the 20 fifth-order correction choices above
the fixed lower branch is obstructed at order seven by the same dual
cocycle.

The global extension uses the six variation columns numbered
\(9,11,13,15,17,19\) and the pure-\(X\) rows \(X^8,\ldots,X^{13}\).
Their pivot matrix \(P\) is constant and invertible over \(K\), and every
other variation column lies in their span over the full reduced
\(\mathbb A^{27}_K\).  If \(V_J\) denotes those pivot columns, then

\[
\Lambda^{\mathrm{rel}}_{18}
=e_{X^{18}}^\vee-(V_J)_{X^{18}}P^{-1}
 (e_{X^8}^\vee,\ldots,e_{X^{13}}^\vee)
\tag{RDC19a}
\]

annihilates all twenty columns identically.  Exact evaluation gives the
same constant (RDC5), independently of all 27 lower-lift parameters.
Moreover

\[
\begin{aligned}
&(32313555201a^2+79786133680a+49319661920)\\
&\quad\cdot(587583566a^2+1388701707a+831388850)\\
&-(201988446756823689a+256263622855091438)\\
&\quad\cdot(94a^3+335a^2+400a+160)
=1637349242961920.
\end{aligned}
\tag{RDC19b}
\]

Thus (RDC5) is a unit in \(K\).  The complete 401-generator consistency
ideal in ten effective parameters contains 27 nonzero constants, has exact
characteristic-zero Gröbner basis \((1)\), and admits a one-generator
degree-zero Nullstellensatz certificate.  Hence every point of the reduced
\(\mathbb A^{27}_K\) component is obstructed at order seven.

## 6. Exact scope

The cocycle statement is relative to the declared deformation problem:
Weyl ordering, the parity-preserving correction tower, the stated
filtration and differential-order bounds, and (for \(\Lambda_5\)) finite
Laurent localization at \(X\).  It does not exclude:

* odd-\(\hbar\) corrections;
* a different ordering prescription;
* higher differential order or relaxed filtration;
* a different polarization; or
* Hamiltonian reduction followed by quantization.

Accordingly this is a rigorous obstruction for the displayed symbol and
ansatz, not a proof of \(DC_2\).

## 7. First parameter-uniform presentation

Let \(\Sigma\) be the fixed 16-monomial support of the all-pole functional
at \((\kappa,\tau)=(0,1)\), and write
\(\mathcal P_{5,\Sigma}\) for strong cocycles supported on \(\Sigma\).
This is a support-restricted submodule of \(\mathcal P_5\), not yet the
complete dual obstruction sheaf.

Over \(\mathbb Q(a,\tau)\), all 41 nonzero conditions on this support have
rank 15.  The stable subsystem consists of ten order-five correction-image
rows and five lower-lift rows.  Clearing those rows primitively gives a
polynomial matrix

\[
M_\Sigma(a,\tau):B^{16}\longrightarrow B^{15},
\tag{RDC20}
\]

where the selected sparse lower-kernel basis is defined on

\[
B=\mathbb Q\!\left[
a,\tau,\frac1{a(a+1)\tau H(a,\tau)}
\right]
\]

and

\[
\begin{aligned}
H(a,\tau)={}&
4a^3\tau^2-24a^3\tau-72a^3
+8a^2\tau^2-54a^2\tau-216a^2\\
&+4a\tau^2-30a\tau-246a-105.
\end{aligned}
\tag{RDC21}
\]

The frozen presentation is
[`degree_five_qper_15x16_presentation.sing`](../artifacts/generated-results/degree_five_qper_15x16_presentation.sing).

Two maximal minors give explicit trivialization charts.  Deleting the first
and last functional coordinates respectively gives

\[
\begin{aligned}
\Delta_0&\doteq a^{12}(a+1)^{62}D_{34}(a,\tau),\\
\Delta_{15}&\doteq a^{10}(a+1)^{63}E_{35}(a,\tau),
\end{aligned}
\tag{RDC22}
\]

where \(D_{34}\) and \(E_{35}\) are irreducible, primitive, and have:

\[
\begin{array}{c|ccc}
&\text{terms}&\deg_{\rm total}&(\deg_a,\deg_\tau)\\ \hline
D_{34}&274&34&(22,12)\\
E_{35}&283&35&(23,12).
\end{array}
\tag{RDC23}
\]

The exact fraction-free gcd is

\[
\gcd(\Delta_0,\Delta_{15})
\doteq a^{10}(a+1)^{62}.
\tag{RDC24}
\]

Hence the rank-\(\le14\) Fitting locus of this presentation has no
codimension-one component on the selected localization
\(a(a+1)H\tau\ne0\).
Possible failure of the cocycle line is confined to a finite scheme inside
\(V(D_{34},E_{35})\).

The first finite-stratum audit is modular but stable across three good
primes.  At each of \(p=31991,32003,65521\):

\[
\begin{array}{c|cc}
\text{ideal}&\dim&\text{quotient length}\\ \hline
(D_{34},E_{35}):H^\infty&0&395\\
I_{15}(M_\Sigma):(a(a+1)H)^\infty&0&218.
\end{array}
\tag{RDC25}
\]

Thus all sixteen maximal minors reduce the two-chart complement from length
395 to a candidate Fitting stratum of length 218.  Eight exact points on
\(D_{34}\) over \(\mathbf F_{101}\) and \(\mathbf F_{103}\), all with
\(H\ne0\), retain rank 15 for the complete 41-condition fiber.  This
confirms directly that \(D_{34}\) is a normalization-chart boundary, not by
itself an extra-cocycle divisor.

The three saturated modular bases have not merely the same length but the
same leading ideal:

\[
\begin{split}
J_{\rm in}={}&
 (a^8\tau^{12},a^9\tau^{11},\ldots,a^{19}\tau,a^{20},\\
 &\tau^{21},a\tau^{20},\ldots,a^7\tau^{14}).
\end{split}
\tag{RDC26}
\]

The saturation exponent is 12 at all three primes.  The standard monomials
under the first row of (RDC26) contribute \(12+11+\cdots+1=78\);
those under the second contribute \(21+20+\cdots+14=140\).  This gives
\(78+140=218\) directly from the common staircase.  This remains modular
evidence: agreement of several good reductions does not prove that the
characteristic-zero leading ideal is (RDC26).

The agreement in (RDC25) is not yet a characteristic-zero length proof.
It does, however, permit a sharply bounded reconstruction.  The modular
lift in
[`compute_degree_five_qper_fitting.py`](../scripts/compute_degree_five_qper_fitting.py)
produces a 21-element rational candidate \(G_{\mathbf Q}\), stored in
[`degree_five_qper_fitting_basis_Q.sing`](../artifacts/generated-results/degree_five_qper_fitting_basis_Q.sing)
with SHA-256
`25788668021f563e17373b55703a08ef5693576077ebdbe53c4c3f2c659d98e6`.
The file is \(20{,}840{,}615\) bytes because its reduced coefficients have
thousands of digits.  The completed and timed checks are recorded separately
in
[`degree_five_qper_fitting_basis_certificate.txt`](../artifacts/generated-results/degree_five_qper_fitting_basis_certificate.txt).

Two properties of this rebuilt candidate now have independent exact
certificates.  First, the 20 \(S\)-polynomials between adjacent generators
in the two-variable staircase (RDC26) reduce to zero.  Buchberger's chain
criterion therefore proves that \(G_{\mathbf Q}\) is a Gröbner basis, without
checking all 210 pairs.  Second,

\[
1\in (G_{\mathbf Q},a(a+1)H),
\tag{RDC27}
\]

so \(G_{\mathbf Q}\) is already saturated by the selected boundary product.
Consequently its dimension-zero quotient has exact length 218.

This still does **not** identify \(G_{\mathbf Q}\) with
\(I_{15}(M_\Sigma):(a(a+1)H)^\infty\).  Exact input containment
\(I_{15}(M_\Sigma)\subseteq G_{\mathbf Q}\) and reverse membership
\((a(a+1)H)^{12}G_{\mathbf Q}\subseteq I_{15}(M_\Sigma)\) remain.  A direct
normal-form computation for the seventh input generator did not finish
within 30 minutes.  The correct next certificate is therefore a
fraction-free reconstruction of the quotient identities for these two
containments, followed by direct polynomial-identity checks; rerunning
coefficient-heavy rational normal forms is not the preferred route.

A first bounded input-containment reconstruction now makes that warning
quantitative.  Ordered modular division of all 16 primitive minors by
\(G_{\mathbf Q}\) has 336 nonzero quotient entries and 11,701 quotient
monomials.  Among 615 computed prime records, 613 have exactly the same
support and zero remainder and residual; the two support-unlucky primes are
70067 and 70099.  The dominant support has SHA-256
`dde5e0b72fad4ba9532e37e49361d48c6d001635d836bec1a5314ec674e2c548`.
The compact checkpoint is
[`degree_five_qper_input_quotients_modular.json`](../artifacts/generated-results/degree_five_qper_input_quotients_modular.json).

This is strong modular evidence for
\(I_{15}(M_\Sigma)\subseteq G_{\mathbf Q}\), but it is not an exact
containment certificate.  A CRT pool of 612 good primes has 18,116 bits,
with the distinct good prime 1000012337 held out.  Balanced rational
reconstruction produces 7,087 candidates, but only 30 of 11,701
coefficients agree with the held-out image.  The validated coefficients
already reach 7,660 numerator bits and 7,487 denominator bits.  Thus blindly
extending the canonical ordered-division quotients is not a reasonably
bounded route.  The next input-containment experiment should instead choose
a noncanonical lift modulo the syzygy module—equivalently, a fixed
fraction-free Macaulay/RREF pivot convention—designed to minimize coefficient
height.  Exact integer polynomial identities remain the terminal
certificate.

## 8. Reproduction and next step

The all-pole strong-cocycle calculation is:

```bash
.venv/bin/python \
  scripts/verify_degree_five_laurent_quantum_obstruction.py
```

The exact cubic branch and its one-term order-seven cocycle are:

```bash
.venv/bin/python \
  scripts/analyze_degree_five_cubic_fifth_order.py \
  --exact-cubic --seventh-order
```

The global order-seven section, generated unit program, and fast replay are:

```bash
PYTHONPATH=scripts .venv/bin/python \
  scripts/analyze_degree_five_cubic_fifth_order.py \
  --exact-cubic --seventh-component-elimination \
  --seventh-component-program-output \
    artifacts/generated-results/degree_five_cubic_h7_unit_certificate.sing
.venv/bin/python \
  scripts/verify_degree_five_cubic_h7_unit_certificate.py
```

The parameter presentation and modular Fitting audit are:

```bash
.venv/bin/python \
  scripts/compute_degree_five_fifth_order_function_field.py \
  --presentation-factor-output \
    artifacts/generated-results/degree_five_qper_pivot_D34.sing \
  --alternate-presentation-factor-output \
    artifacts/generated-results/degree_five_qper_pivot_E35.sing \
  --presentation-program-output \
    artifacts/generated-results/degree_five_qper_15x16_presentation.sing

PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_five_qper_presentation_strata.py

.venv/bin/python scripts/compute_degree_five_qper_fitting.py \
  --timeout 120

# Rebuild the rational candidate (about 11 minutes on the reference machine).
.venv/bin/python scripts/compute_degree_five_qper_fitting.py \
  --prime 0 --method modular-rebuild --timeout 900 \
  --basis-output \
    artifacts/generated-results/degree_five_qper_fitting_basis_Q.sing

# Exact properties of the rebuilt candidate.
.venv/bin/python scripts/verify_degree_five_qper_fitting_basis.py \
  --check shape
.venv/bin/python scripts/verify_degree_five_qper_fitting_basis.py \
  --check groebner --jobs 8 --timeout 600
.venv/bin/python scripts/verify_degree_five_qper_fitting_basis.py \
  --check boundary-unit --timeout 600
```

The next parameter-uniform steps are:

1. reconstruct fraction-free quotient identities proving both containments
   between the rebuilt \(G_{\mathbf Q}\) and the saturated maximal-minor
   ideal;
2. compare that finite scheme with the vertical fibers obtained from the
   other four supports;
3. lift the bounded-versus-Laurent equality of the 41-condition module over
   the parameter ring; and
4. start the same relative construction on a different classical-symbol
   family now that the cubic order-seven component is closed.
