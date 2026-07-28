# One deformation complex for gauges, corrections, and defects

> **Status and scope.**  This note defines a common computational interface
> and a research protocol.  It does not assert that left--right deformation,
> exact symplectic deformation, Weyl quantization, and constant-Hessian
> deformation are equivalent moduli problems.  The three-term complex below
> is exact formal packaging; every application must still specify its
> filtration, parameter base, gauge group, relation module, and boundary
> lattice.  The small checker verifies four exact field-level
> instantiations.  It does not replace the existing large certificates until
> their matrices have been rebuilt through this interface.

Several calculations in the repository repeatedly form a large matrix,
compute its kernel and left cokernel, then start again when the filtration,
prime, or classical seed changes.  Their common linear core is

\[
 C^0
 \xrightarrow{\,d_0\,}
 C^1
 \xrightarrow{\,d_1\,}
 C^2:
\qquad
 \text{gauges}\longrightarrow
 \text{corrections}\longrightarrow
 \text{relation defects}.
\tag{1}
\]

The useful replacement is not one universal matrix.  It is one functorial
complex whose modules and sparse differentials are assembled from the
chosen deformation problem.  Parameter components, higher Kuranishi maps,
localization, and global descent are then layers on the same presentation.

## 1. The relative correction complex

Let \(A\) be a characteristic-zero coefficient ring or a localization of
one.  Fix a classical object \(G\), a filtration order \(m\), and finite
free or perfect \(A\)-modules

\[
 C^0_{G,m},\qquad C^1_{G,m},\qquad C^2_{G,m}.
\]

They mean:

* \(C^0_{G,m}\): allowed infinitesimal gauge transformations;
* \(C^1_{G,m}\): allowed order-\(m\) corrections; and
* \(C^2_{G,m}\): coefficients of the order-\(m\) relation defects.

Differentiate the gauge action to obtain \(d_0\), and differentiate the
defining relations at \(G\) to obtain \(d_1\).  The gauge group must preserve
the declared relation problem.  Under precisely that hypothesis,

\[
 d_1d_0=0.                                             \tag{2}
\]

The cohomology has the operational interpretation

\[
\begin{aligned}
H^0&=\ker d_0
 &&\text{infinitesimal stabilizers},\\
H^1&=\ker d_1/\operatorname{im}d_0
 &&\text{homogeneous corrections modulo gauge},\\
H^2&=\operatorname{coker}d_1
 &&\text{obstruction classes}.
\end{aligned}                                         \tag{3}
\]

After lower orders have been chosen, the next equation is

\[
 d_{1,G,m}(\Delta_m)=-\mathcal O_m.                    \tag{4}
\]

It is solvable exactly when

\[
 [\mathcal O_m]=0\quad\text{in }H^2_{G,m}.             \tag{5}
\]

A row functional \(\Lambda\in(C^2)^\vee\) with
\(\Lambda d_1=0\) is a dual obstruction cocycle.  If
\(\Lambda(\mathcal O_m)\ne0\), it is an exact inconsistency certificate
without requiring a full row-reduced matrix in the final record.

This is the field-level calculation implemented in
[`deformation_complex.py`](../jcsearch/deformation_complex.py).  It checks
(2), computes the dimensions (3), returns a basis of dual obstruction
cocycles, evaluates (5), solves correctable defects with a normalized
particular correction, and compares exact ranks at good primes.

## 2. Lower lifts and Kuranishi maps

The linear complex alone does not describe the dependence of
\(\mathcal O_m\) on earlier choices.  Let

\[
 p:\mathcal L_{<m}\longrightarrow B
\]

be the scheme of lower-order lifts over a classical-symbol base \(B\).
The order-\(m\) obstruction is a section

\[
 \overline{\mathcal O}_m
 \in
 \Gamma\!\left(
 \mathcal L_{<m},p^*H^2_m
 \right).                                             \tag{6}
\]

Its zero scheme is the next Kuranishi locus.  Quadratic, cubic, and higher
Schur/Kuranishi equations are simply the homogeneous pieces of (6).
Changing a lower lift moves inside \(\mathcal L_{<m}\); it is not captured
by quotienting only the current correction space.

The strong dual-cocycle module from
[Restricted quantum deformation cocycles](RESTRICTED_QUANTUM_DEFORMATION_COCYCLES.md)
fits directly here.  A dual cocycle is **strong** when its evaluation on
(6) is pulled back from \(B\), so its value is independent of every
lower-lift choice.  Thus the existing rank-jump certificates have the form

\[
 \Lambda d_1=0,\qquad
 \Lambda(\mathcal O_{m,\alpha})=0
 \quad(\alpha\ne0),\qquad
 \Lambda(\mathcal O_{m,0})\ne0.                     \tag{7}
\]

Equation (7) should be the common certificate format for Weyl, Schur, and
coefficient-scheme obstruction calculations.

## 3. Vary the classical symbol, not only its correction support

A calculation at one point \(b\in B\) answers only whether that classical
symbol survives the chosen quantization problem.  It cannot settle a
universal conjecture whose putative counterexample may have another
principal symbol.

The correct family object is the complex \(C^\bullet_m\) over \(B\), not a
single specialized matrix.  Its coherent obstruction module is

\[
 E_m=\operatorname{coker}(d_{1,m}).                  \tag{8}
\]

The rank and Fitting strata of (8) identify where correction dimensions or
obstruction dimensions jump.  The section (6) then cuts out the locus of
classical symbols that continue to order \(m\).

At a zero of the obstruction section, the first absolute linearization can
be written locally as

\[
 C^0_{b,m}
 \longrightarrow
 C^1_{b,m}\oplus T_bB
 \xrightarrow{\ (d_{1,b,m},\,D_b\mathcal O_m)\ }
 C^2_{b,m}.                                          \tag{9}
\]

The second summand is the missing horizontal direction in a fixed-symbol
calculation.  Nonlinear components still require the full Kuranishi section
(6), not only its derivative (9).

For reconstruction across primes, the safe protocol is:

1. freeze an integral basis, support, filtration, and term order;
2. compute ranks, Fitting ideals, and component fingerprints at several
   good primes;
3. match only components whose dimensions, degrees, leading ideals, and
   incidence with named boundary divisors agree;
4. reconstruct candidate equations over \(\mathbb Q\); and
5. verify the reconstructed containments, saturations, and obstruction
   identities exactly in characteristic zero.

Agreement of ranks at several primes is a component-discovery screen, not a
characteristic-zero proof.  Small-characteristic Frobenius components must
be labelled vertical and must not be reconstructed as horizontal
characteristic-zero components.

## 4. Four instantiations

### 4.1 Commutative left--right and coefficient schemes

For a polynomial map \(F\), let \(C^1\) be a bounded space of map
perturbations \(G\).  For the determinant-one coefficient scheme, let
\(C^2\) be the declared coefficient space of
\(\det D(F+\epsilon G)-1\).  Then

\[
 d_1(G)
 =
 L_F(G)
 =
 \operatorname{tr}\!\left(
 \operatorname{adj}(DF)\,DG
 \right).                                             \tag{10}
\]

For admissible source and target vector fields \(V,U\),

\[
 d_0(U,V)=U\circ F+DF\cdot V.                         \tag{11}
\]

The admissibility condition

\[
 (\operatorname{div}U)\circ F+\operatorname{div}V=0   \tag{12}
\]

is exactly what makes (10)--(11) a complex.  If arbitrary left--right
fields are used, their determinant variation must be retained as an extra
relation coordinate rather than silently discarded.

At the explicit Keller maps, the unrestricted first-order source gauge
exhausts \(\ker L_F\).  Therefore the useful information begins with the
filtered complex and its quadratic second fundamental form, not the
ordinary tangent quotient.  This recovers the organization of
[Jelonek coefficient components](JELONEK_COEFFICIENT_COMPONENTS.md) and
[the torus-filtered LR module](TORUS_FILTERED_LR_MODULE.md) in one
presentation.

### 4.2 Exact symplectic maps

For a canonical pair \((S,T)\), take \(C^0\) to be allowed Hamiltonians and
\(C^1=U\oplus V\) to be pair corrections.  The first bracket defect is

\[
\begin{aligned}
d_0H&=(\{H,S\},\{H,T\}),\\
d_1(u,v)&=\{u,T\}+\{S,v\}.
\end{aligned}                                         \tag{13}
\]

Jacobi and \(\{S,T\}=1\) give \(d_1d_0=0\).  For several canonical pairs,
\(C^2\) contains all diagonal and mixed bracket defects.  In localized
Darboux coordinates, the \(R\)-connection equations are additional mixed
defect rows, not a separate kind of calculation.

### 4.3 Weyl commutators

For an exact Weyl tuple \(Y=(Y_i)\), let \(C^1\) contain the filtered PBW
corrections \(\Delta_i\), and let \(C^2\) contain one coordinate block for
each Weyl-relation defect.  For a relation
\([Y_i,Y_j]=\omega_{ij}\hbar\),

\[
 d_1(\Delta)_{ij}
 =
 [\Delta_i,Y_j]+[Y_i,\Delta_j].                      \tag{14}
\]

Inner filtered gauges give

\[
 d_0(H)=([H,Y_i])_i.                                  \tag{15}
\]

Jacobi and centrality of \(\omega_{ij}\hbar\) prove (2).  Formula (14) is
the operator-level version of the filtered equation

\[
 d_G(\Delta_m)=-\mathcal O_m.
\]

The Ore localization is not merely a faster basis choice.  It splits three
relations exactly and reduces (14) to a rank-one fiber Weyl complex plus
the Hamiltonian connection.  This should be used before any multi-PBW
support enumeration.

### 4.4 Hessian determinant constraints

For a potential \(\Phi\), let \(C^1\) be the selected potential corrections
\(h\), and let \(C^2\) contain the nonconstant coefficients of the Hessian
determinant defect.  The linear differential is

\[
 d_1(h)
 =
 \operatorname{tr}\!\left(
 \operatorname{adj}(\operatorname{Hess}\Phi)
 \operatorname{Hess}h
 \right).                                             \tag{16}
\]

The gauge module \(C^0\) must be the chosen normalized coordinate and
potential gauge that preserves the declared Hessian relation.  Schur
complements reduce (16) before expansion; the higher determinant layers
are the Kuranishi pieces of (6).  The same construction applies cellwise
to the Hessian--Ritt diagram, where its totalization is the larger
[Hessian--Ritt deformation complex](HESSIAN_RITT_DEFORMATION_COMPLEX.md).

## 5. Localization and boundary descent

The rank-two Weyl problem is already simplest after Ore localization.
Suppose \(A\) is the global coefficient ring, \(\mathfrak b\) is the
boundary ideal, and the relative complex becomes easy over
\(A[\mathfrak b^{-1}]\).  If the localized order-\(m\) equation is soluble,
then the global obstruction class

\[
 e_m=[\mathcal O_m]\in E_m=\operatorname{coker}d_{1,m}
\]

is supported on \(V(\mathfrak b)\):

\[
 e_m\in H^0_{\mathfrak b}(E_m).                      \tag{17}
\]

Consequently the support-saturation condition

\[
 H^0_{\mathfrak b}(E_m)=0,
\quad\text{equivalently}\quad
 \operatorname{im}d_{1,m}:_{C^2}\mathfrak b^\infty
 =\operatorname{im}d_{1,m},                          \tag{18}
\]

forces the obstruction to vanish globally.  This is precisely the
[support-saturation principle](../verified/SUPPORT_SATURATION_PRINCIPLE.md).
It proves existence of a global correction; it does not say that a
particular localized formula is already pole-free.

Apparent denominators must therefore be treated in the quotient torsor:

1. solve the localized correction equation;
2. compute boundary valuations modulo \(\ker d_1\) and
   \(\operatorname{im}d_0\), since a different homogeneous correction or
   gauge may cancel the displayed poles;
3. use (18) to rule out a defect class supported only on the boundary; and
4. glue solutions across normalization charts with the conductor square.

For a finite normalization \(A\subset\widetilde A\) with conductor
\(\mathfrak c=\operatorname{Ann}_A(\widetilde A/A)\), the basic gluing
diagram compares data over \(\widetilde A\), \(A/\mathfrak c\), and
\(\widetilde A/\mathfrak c\).  Module descent requires the corresponding
Tor or flatness hypotheses; these must be checked rather than inferred from
set-theoretic agreement.  The conductor records exactly the finite boundary
compatibility that localization forgets.

## 6. Incidence-first rank-two quantization

The primary \(DC_2\) branch should now run in the following order.

### Step A: quantize the natural incidence chart

Use the marked simple-root presentation

\[
 I^{\mathrm{simp}}\longrightarrow\mathbb A^3
\]

from [the marked-root model](../verified/MARKED_ROOT_MODEL.md), or the
equivalent normalized factorization chart, as the classical incidence
object.  For a family, let \(B\) parametrize the seed or classical symbol
and form the relative marked-root object over \(B\).  Do not specialize
immediately to the one degree-five rank-two symbol.  Keep the root, residual
factor, and seed parameters in the coefficient ring.

On the affine-root chart, the coordinates

\[
 u=y+x^{-1},\qquad v=x^{-1}
\]

solve the inverse/incidence problem.  On the second chart, the marked root
at infinity has the regular reconstruction in the coordinate \(s\).
The quantum calculation must record the exact bridge from these root
coordinates to the Ore variables; it should not rely on their similar
notation.  On the generic adapted chart this bridge is now exact.  If
\(W=Z+\psi(Q)\), \(Y=Q-XW/3\), and
\(y_{\rm seed}=-3Y/(2a)\), then

\[
\begin{aligned}
u&=y_{\rm seed}+X^{-1}
 =v-\frac{3Q}{2a}+\frac{Z+\psi(Q)}{2av},\\
v&=X^{-1},\qquad
Q=\frac{v(2-Rv)}3.                                  \tag{20}
\end{aligned}
\]

Moreover the argument at which the weighted seed polynomial is evaluated
is

\[
 (1+Xy_{\rm seed})\left(1-\frac32XQ\right)
 =\frac{Ru}{2}.                                     \tag{21}
\]

Thus on the root-at-infinity chart \(s=u^{-1}\), \(d=v/u\), a finite marked
argument \(w\) forces

\[
 R=2sw,\qquad
 1-\frac32XQ=dw.                                    \tag{22}
\]

These identities are checked by
[`verify_marked_root_ore_bridge.py`](../scripts/verify_marked_root_ore_bridge.py).
They give the first boundary valuation lattice; they do not yet prove that
the quantum corrections lie in it.

### Step B: solve fiber quantization and the connection

In the Ore chart, the exact Darboux variables satisfy

\[
 [P,v]=[U,R]=\hbar
\]

with all mixed commutators zero.  First solve the rank-one fiber equation

\[
 [S_\hbar,T_\hbar]=\hbar,
\]

over the incidence base.  Then solve

\[
\begin{aligned}
[A_\hbar,S_\hbar]&=-\hbar\partial_RS_\hbar,\\
[A_\hbar,T_\hbar]&=-\hbar\partial_RT_\hbar.
\end{aligned}                                        \tag{19}
\]

The formal-local contracting homotopy and the formal Hamiltonian primitive
show that both stages are locally unobstructed after completion.  The new
content is to construct them relatively over classical-symbol components
and in a boundary-controlled lattice.

The first relative complex is now explicit in canonical coordinates.  For
corrections \((s,t,a)\), its raw defect operator is

\[
 d_1(s,t,a)=
 \left(
 s_S+t_T,\,
 s_R-a_T,\,
 t_R+a_S
 \right).                                           \tag{23}
\]

Hamiltonian gauge is

\[
 d_0(h)=(-h_T,h_S,-h_R).                            \tag{24}
\]

The raw defects satisfy the Bianchi identity

\[
 \partial_R F-\partial_S G-\partial_T H=0.          \tag{25}
\]

Thus \(C^2\) is the closed-defect module, not three unrelated coefficient
blocks.  In polynomial degrees \(4\to3\to2\), the exact complex has
dimensions and ranks

\[
 35\xrightarrow[\operatorname{rank}34]{d_0}
 60\xrightarrow[\operatorname{rank}26]{d_1}
 26,
\qquad
 (\dim H^0,\dim H^1,\dim H^2)=(1,0,0).              \tag{26}
\]

The integer matrices base-change to every characteristic-zero seed algebra.
The construction and three-prime regression are
[`verify_relative_fiber_connection_complex.py`](../scripts/verify_relative_fiber_connection_complex.py).
Equation (26) proves local bounded exactness in the canonical chart; it does
not control the Ore filtration or boundary lattice.

### Step C: make polynomiality a boundary theorem

Compute valuations of the relative corrections on the root-at-infinity
chart.  Test cancellation after all homogeneous corrections and Hamiltonian
gauges, not term by term in one representative.  Present the remaining
class in \(E_m\), apply support saturation at every boundary divisor, and
use the normalization conductor to glue the two marked-root charts.

This is the quantum analogue of the classical construction: solve the easy
inverse/incidence problem first, then prove that the apparent denominators
cancel and that the result extends polynomially.

### Step D: branch across classical symbols

At each order retain the Kuranishi/Fitting locus in the incidence-symbol
base.  Reconstruct its horizontal components across several good primes and
verify the candidates over \(\mathbb Q\).  Sample different seed parameters,
degrees, and classical rank-two completions rather than only higher-support
corrections of one symbol.

Continuing the high-support \(\mathbf P^6\) of the current degree-five
symbol remains a valid bounded fiber calculation.  Eliminating it, or every
branch over that one symbol, would prove only nonquantizability of that
symbol in the declared filtration.  It would not prove \(DC_2\).

### Step E: separate global generation

Even a global filtered Weyl tuple is only an endomorphism candidate.
Nonsurjectivity still needs an Ore--Gröbner, module, center, or other exact
image-subalgebra certificate.  A non-Keller tuple of separate leading
symbols remains metadata, not a nonsurjectivity test.

## 7. Executed degree-five branch

The incidence-first workflow has now been run on the current degree-five
two-parameter symbol family.  This is a bounded exact calculation for that
family, not a result about all classical symbols.

1. The order-three correction system over \(\mathbb Q(a,\tau)\) has 42
   columns, rank 42, zero kernel, and zero residual.  Its denominators have
   only the recorded powers of \(a\); no additional parameter divisor is
   introduced.
2. At order five, four independent period supports were reconstructed
   exactly in characteristic zero after modular discovery.  Their primitive
   factors have total degrees \(41,32,38,39\).
3. Their exact common scheme has length 44, equal to the modular upper
   bound.  Hence the listed support is exhausted.  Away from the declared
   boundary strata, it contains the rational point
   \((a,\tau)=(-1/2,-3)\) and one residue-degree-three component
   \[
      94a^3+335a^2+400a+160=0,\qquad
      8\tau+658a^2+1593a+976=0.
   \]
4. The rational point has an exact all-pole order-five obstruction.  Its
   dual period annihilates every allowed finite Laurent correction, so
   localizing at the adapted coordinate does not remove it.
5. Over the cubic residue field there is an explicit order-five lift.  The
   induced order-seven extension has coefficient rank 6 and augmented rank
   7.  An exact one-term dual period therefore obstructs this branch at
   order seven.

Thus every nonboundary component found in this particular order-five
symbol scheme terminates before conductor descent.  There is no global Weyl
tuple from these branches, so the image-subalgebra/nonsurjectivity test in
Step E is not reached.  This closes the declared branches only; it is the
reason the next run must begin with different classical symbols.

The maximal-minor presentation also yields a 21-element exact Gröbner basis
candidate with quotient length 218 and a boundary-unit certificate.  Exact
identification of that candidate with the saturated original maximal-minor
ideal still requires both containment directions, so it is not used to
strengthen the component conclusion above.

## 8. What this replaces, and what it does not

The reusable engine should replace repeated implementations of:

* kernel and gauge-quotient construction;
* left-cokernel projection;
* dual obstruction-cocycle extraction;
* rank comparison at good primes; and
* assembly of sparse differentials from local relation formulas.

It does not replace:

* nonlinear component reconstruction;
* exact characteristic-zero verification after modular discovery;
* Kuranishi brackets and completed local algebras;
* support saturation or conductor descent; or
* an image-subalgebra certificate for a Weyl endomorphism.

Those are functorial layers attached to (1), not more rows in a single
matrix.

## Reproduction

The exact small regression is
[`verify_unified_deformation_complex.py`](../scripts/verify_unified_deformation_complex.py):

```bash
.venv/bin/python scripts/verify_unified_deformation_complex.py
```

It constructs four complexes directly from their defining operations:

* the derivative constraint on a determinant-one commutative
  left--right slice;
* Hamiltonian gauge and Poisson-bracket defects for a canonical pair;
* exact PBW multiplication with \([p,q]=1\) for a Weyl pair; and
* the linearized constant-Hessian determinant at a quadratic potential,
  with cubic corrections and quadratic defects retained as a restricted
  obstruction test.

The symplectic and exact PBW examples have the same dimensions, ranks, and
cohomology, while their differentials are independently assembled.  A
quadratic Hessian defect survives the restricted cubic-correction cokernel.
Three good-prime rank checks are included.  This is a regression of the
common linear interface, not a theorem about the nonlinear or global
comparison.

The executed degree-five branch is reproduced by the existing exact
function-field, period-component, Laurent-obstruction, and cubic-component
scripts.  In particular, after writing the four primitive period factors to
temporary files, use

```bash
.venv/bin/python scripts/analyze_degree_five_period_components.py \
  --exact-decomposition P41.txt Q32.txt R38.txt S39.txt
PYTHONPATH=scripts .venv/bin/python \
  scripts/verify_degree_five_laurent_quantum_obstruction.py \
  --rational-common-point
PYTHONPATH=scripts .venv/bin/python \
  scripts/analyze_degree_five_cubic_fifth_order.py \
  --exact-cubic --seventh-order
```

The factor files are generated inputs, not pinned certificates; the four
`compute_degree_five_fifth_order_function_field.py` support modes print the
exact output paths and reconstruction metadata.
