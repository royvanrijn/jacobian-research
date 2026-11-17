# Universal cubic quartic-kernel saturation frontier

## Status

The full 24-parameter relative-cotangent saturation question for the smooth
cubic symbol is proved.  The proof is a formal-gauge certificate, not a full
27-variable saturation:

\[
 H^0_{(x,y,z)}(\Omega_{B/R})=0.
\]

Consequently the universal annihilator equals the canonical different and
the actual-support Fittings are
\(\operatorname{Fitt}_6=(1)\) and
\(\operatorname{Fitt}_5=(0)\).

This note separates three levels which must not be conflated:

1. exact polynomial-family calculations on specified parameter subspaces;
2. exact calculations at isolated dense parameter points;
3. the universal formal-gauge calculation over
   \(\mathbb Q[u_1,\ldots,u_{24},x,y,z]\).

The earlier subspace calculations remain useful independent checks.  The
universal theorem is supplied by Section 5 below.

For the six singular squarefree symbols, the smooth formal-rigidity proof
does not apply, but the later route is also closed at its stated level.
`KDSQ6` proves cotangent saturation and a six-generated non-Cartier
different on every fiber of each complete quartic nongauge complement, and
`SSADPALL` propagates those statements through every compatible formal tail.
These results do not prove normality, algebraize an infinite formal gauge,
or construct a distinguished Keller open.  The double-line, triple-line,
and zero symbols remain outside the squarefree theorem and must first pass
generic-etaleness and Keller-compatibility gates.

## 1. Universal family

Let

\[
 A=\mathbb Q[x,y,z],\qquad
 M=\operatorname{coker}\left(
 A\xrightarrow{(z,-y,x)^{\mathsf T}}A^3
 \right).
\]

For a squarefree ternary-cubic representative \(h\), write \(\Phi_h\) for
the homogeneous generalized triple-cover tensor on \(M\).  The order-four
compatibility equations have the fixed primitive integral basis
\(\psi_1,\ldots,\psi_{24}\) constructed by
[`verify_cubic_symbol_double_saturation.py`](../scripts/verify_cubic_symbol_double_saturation.py).
Put

\[
 S=\mathbb Q[u_1,\ldots,u_{24}],\qquad
 R=S[x,y,z],
\]

and form the exact universal tensor

\[
 \Phi_{\mathrm{univ}}
 =\Phi_h+\sum_{i=1}^{24}u_i\psi_i.                 \tag{1.1}
\]

The checker reconstructs from (1.1) the associative multiplication on
\(R\oplus(M\otimes_A R)\), the 31-column presentation \(N\) of
\(\Omega_{B/R}\), the ramification-support module

\[
 T=B/\operatorname{Ann}_B(\Omega_{B/R}),
\]

and

\[
 E=\operatorname{Ext}^2_R(T,R).                    \tag{1.2}
\]

The universal theorem is:

\[
\begin{aligned}
&(N:(x,y,z)^\infty)/N=0,\\
&\sqrt{\operatorname{Fitt}_0^R(E)}=(x,y,z),\\
&E\text{ is }S\text{-flat of rank }6.
\end{aligned}                                      \tag{1.3}
\]

For the finite \(S\)-module obtained from the verified
\((x,y,z)^2E=0\) truncation, flat rank six is equivalent to

\[
 \operatorname{Fitt}_6^S(E)=(1),\qquad
 \operatorname{Fitt}_5^S(E)=(0).                   \tag{1.4}
\]

The first equality is proved in Section 5.  The canonical-different
calculation then supplies the last two assertions without a separate
parameter-discriminant elimination.

## 2. Exact full-support planes and lines

The strongest current mixed-support calculation is KDFP6.  With
\(\psi_+=\sum_i\psi_i\) and
\(\psi_-=\sum_i(-1)^i\psi_i\), it treats

\[
 \Phi_h+u\psi_++v\psi_-
\]

over \(\mathbb Q[u,v,x,y,z]\) for every squarefree symbol.  The cotangent
presentation is saturated and the pruned rank-three Ext presentation is
pulled back from the origin with multiplicity six.  On
\(D(u^2-v^2)\), all 24 fixed basis coordinates are nonzero.  Thus basis
sparsity is not a necessary condition for retaining the defect, although
one parameter plane is not an open subset of \(\mathbb A^{24}\).

Four fixed coefficient vectors with all 24 entries nonzero are recorded in
[`research_universal_cubic_quartic_kernel_saturation.py`](../scripts/research_universal_cubic_quartic_kernel_saturation.py).
For each vector \(a=(a_i)\) and each of the seven squarefree cubic-symbol
orbits, the checker works over \(\mathbb Q[t,x,y,z]\) with

\[
 \Phi_h+t\sum_i a_i\psi_i.                          \tag{2.1}
\]

All 28 polynomial families satisfy, exactly:

\[
\begin{array}{c|c}
\text{invariant}&\text{value}\\ \hline
\text{relative cotangent saturation quotient}&0\\
\operatorname{mult}(E)&6\\
t\text{-torsion of }E&0\\
\sqrt{\operatorname{Fitt}_0(E)}&(x,y,z)\\
\text{difference from the central Ext presentation}&0.
\end{array}                                         \tag{2.2}
\]

Consequently (1.4) holds after restriction to each recorded line, and the
universal exceptional locus has empty intersection with all 28 lines.
This is stronger than checking a dense endpoint: it excludes every scalar
on each line.  A finite set of lines is not Zariski dense enough to settle
the 24-parameter question.

## 3. Exact higher-dimensional subspace

For the smooth cubic symbol, the same checker treats the full coordinate
subspace on the first ten basis tensors:

\[
 \Phi_h+\sum_{i=1}^{10}p_i\psi_i
 \quad\text{over}\quad
 \mathbb Q[p_1,\ldots,p_{10},x,y,z].                \tag{3.1}
\]

The relative cotangent presentation is saturated, the radical support of
\(E\) is exactly the parameter ten-space at \(x=y=z=0\), its multiplicity
is six, and its pruned rank-three presentation is pulled back from the
parameter origin.  Hence the parameter discriminant restricts to the empty
scheme on this \(\mathbb A^{10}\).

This extends the previously recorded smooth coordinate-three-space result
for one specified nested subspace.  It does not cover all coordinate
ten-spaces or an arbitrary ten-dimensional linear subspace.

## 4. Universal cotangent input reduction

There is now an exact calculation over the full ring
\(\mathbb Q[u_1,\ldots,u_{24},x,y,z]\), but it is an input reduction rather
than the desired saturation theorem.  Comparing the universal 31-column
cotangent presentation \(N_{\mathrm{univ}}\) with its central value \(N_0\)
entry by entry gives only the parameter/collision bidegrees

\[
 (1,3),\qquad (1,5),\qquad (2,6).                 \tag{4.1}
\]

More precisely, 147 matrix entries change.  Their expanded terms comprise
1,019 terms of bidegree \((1,3)\), 933 of bidegree \((1,5)\), and 4,800 of
bidegree \((2,6)\).  In particular,

\[
 N_{\mathrm{univ}}\equiv N_0\pmod{(x,y,z)^3}.      \tag{4.2}
\]

The collision-degree bound follows structurally from the construction:
the quartic tensor changes the trace-free multiplication in degree three,
its cross term with the cubic multiplication enters the scalar part in
degree five, and its square enters in degree six.  Thus the same bound
applies to every cubic-symbol row, although the serialized matrix hash in
the artifact uses the smooth representative.

Before any standard-basis computation, six parameter-independent unit
pivots split off from the \(12\)-by-\(31\) presentation.  Exact elementary
row and column operations therefore replace it by a cokernel-equivalent
\(6\)-by-\(25\) presentation.  The six successive pivot values are
\(2,1,1,2,1,2\), hence are units over \(\mathbb Q\) on all of parameter
space.  The generated artifact records the pivot positions and a SHA-256
hash of the reduced universal matrix.

Equations (4.1)--(4.2) explain why the verified two-layer Ext module is
insensitive on many parameter slices: no quartic parameter occurs in the
two-jet of the cotangent presentation.  They do **not** prove saturation or
flatness.  Taking the annihilator of the cotangent module and then resolving
its support can convert higher-order relations into new low-order
syzygies.  Any exceptional locus must therefore enter through this
syzygy-lifting step, rather than through the raw two-jet.

### 4.1 The first filtered syzygy frontier

The central reduced \(6\)-by-\(25\) presentation can be resolved exactly.
For the smooth symbol, Singular produces a minimal cokernel resolution

\[
 0\longrightarrow R_0^7\longrightarrow R_0^{13}
 \longrightarrow R_0^6\longrightarrow\Omega_0\longrightarrow0,
\qquad R_0=\mathbb Q[x,y,z].
\tag{4.3}
\]

Thus \(\operatorname{pd}_{R_0}\Omega_0\leq2\), and the exact colon test
shows that \(x+y+z\) is \(\Omega_0\)-regular.  If either this resolution or
this regular element lifted over
\(\mathbb Q[u_1,\ldots,u_{24},x,y,z]\), the universal saturation would
follow by Auslander--Buchsbaum or directly by the regular-element
criterion.

There is, however, a nontrivial finite syzygy correction before that
argument can be made.  After unit pruning, the 150 matrix entries have the
following pairs

\[
(\text{central collision order},\text{perturbation collision order}):
\]

\[
\begin{array}{c|c}
\text{pair}&\text{entry count}\\ \hline
(1,\infty)&15\\
(2,3)&30\\
(3,4)&12\\
(4,5)&51\\
(\infty,3)&21\\
(\infty,4)&6\\
(\infty,5)&3\\
(\infty,\infty)&12.
\end{array}
\tag{4.4}
\]

Here \(\infty\) denotes a zero entry.  Hence every perturbation of a
nonzero central entry strictly raises collision order.  Nevertheless, let
\(S_0\) be the exact 24-column syzygy module of the central
\(6\)-by-\(25\) input and let \(D_{\mathrm{univ}}\) be the universal
matrix.  The exact reduction

\[
 D_{\mathrm{univ}}S_0
 \pmod{\operatorname{im}D_0}
\tag{4.5}
\]

has 12 nonzero columns.  Therefore the unchanged central syzygies do not
give a universal complex.  Entrywise low-jet agreement is not Rees
strictness.

Equation (4.5) is a frontier, not a counterexample to saturation.  Allowed
changes of relation generators and higher filtered corrections may kill
these classes.  Moreover, their span is a submodule of \(\Omega_0\), so
the regularity of \(x+y+z\) on \(\Omega_0\) implies regularity on this
span.  These remainders are therefore not collision-supported torsion
which support saturation could erase; they are a horizontal
presentation-gauge mismatch.  The next finite problem is to construct
corrected syzygies with the allowed generator changes.  The checker
[`verify_universal_cubic_filtered_syzygy_frontier.py`](../scripts/verify_universal_cubic_filtered_syzygy_frontier.py)
records the full twelve-column remainder and its hash.

## 5. Formal-gauge saturation theorem

This section gives the universal certificate.  Put
\(A=\mathbb Q[x,y,z]\), \(r=(z,-y,x)^{\mathsf T}\), and index the ten
components of a symmetric tensor by
\[
 000,001,002,011,012,022,111,112,122,222.
\]
Let \(K\subset A^{10}\) be the graded module of compatible tensor
corrections.  It is the kernel of the \(6\)-by-\(10\) matrix \(C\) encoding
\[
 zc_{0ij}-yc_{1ij}+xc_{2ij}=0
 \qquad(0\leq i\leq j\leq2).                         \tag{5.1}
\]

Write \(\phi_h=(c_{ijk})\) for the smooth Fermat cubic tensor.  A matrix
\(D\in\operatorname{Mat}_3(A)\) determines a simultaneous infinitesimal
change of the collision coordinates and the generators of the Koszul
module.  Define \(v(D)\) by
\[
 (v_z,-v_y,v_x)^{\mathsf T}=Dr.
\]
The induced tensor variation is
\[
 \delta_h(D)_{ijk}
 =
 v(D)\mathbin{\cdot}\nabla c_{ijk}
 +\sum_a\left(
 D_{ai}c_{ajk}+D_{aj}c_{aik}+D_{ak}c_{aij}
 \right)
 -\operatorname{tr}(D)c_{ijk}.                       \tag{5.2}
\]
Here the last term is essential: the generalized tensor takes values in
\(\det(M)\), and transport back to the original determinant line multiplies
by \(\det(I+D)^{-1}\).  Equivalently, the exact transformed tensor is
\[
 \det(I+D)^{-1}\,
 \sigma_D(\phi_h)\bigl((I+D)e_i,(I+D)e_j,(I+D)e_k\bigr),
                                                        \tag{5.2a}
\]
whose first-order term is (5.2).
For the nine matrix units, these columns form a homogeneous
\(10\)-by-\(9\) matrix \(G\).  Direct expansion gives
\[
 CG=0.                                               \tag{5.3}
\]
The checker does not merely assume (5.2): over
\(\mathbb Q[\epsilon]/(\epsilon^2)\), it expands the determinant-twisted
finite action (5.2a) for all nine matrix units and verifies that its
\(\epsilon\)-coefficient is exactly the corresponding column of \(G\).
It simultaneously checks the exact transformed Koszul relation.

There is one cubic direction not induced by a gauge.  Take \(\eta\) to be
the integral tensor attached to \(3XYZ\):
\[
\eta=
(0,-x^2y,-xy^2,-x^2z,0,y^2z,0,xz^2,yz^2,0)^{\mathsf T}.
                                                               \tag{5.4}
\]
Then \(C\eta=0\).  The exact module certificate is
\[
 \boxed{\quad
 \ker C=\operatorname{im}G+A\eta,\qquad
 (x,y,z)\eta\subseteq\operatorname{im}G.
 \quad}                                                   \tag{5.5}
\]
The second identity is witnessed without a standard basis:
\[
 GL=[x\eta\ \ y\eta\ \ z\eta],                         \tag{5.6}
\]
where, in the matrix-unit ordering
\((00),(01),(02),(10),(11),(12),(20),(21),(22)\),
\[
L=
\begin{pmatrix}
0&0&0\\
-y&0&0\\
0&-x&0\\
z&0&0\\
0&0&0\\
0&0&x\\
0&-z&0\\
0&0&-y\\
0&0&0
\end{pmatrix}.                                      \tag{5.7}
\]
For the first identity, the checker computes \(K=\operatorname{syz}(C)\)
over \(\mathbb Q[x,y,z]\) and reduces its ten exact generators by
\([G\ \eta]\), with zero remainder.  Independently,
\(\operatorname{coker}(\operatorname{modulo}(K,G))\) has vector-space
dimension one, \(\eta\) has nonzero remainder modulo \(G\), and
\((x,y,z)\) acts by zero.  Thus, with the grading in which the tensor
components have their polynomial degree,
\[
 \ker C/\operatorname{im}G\simeq\mathbb Q(-3).         \tag{5.8}
\]
In particular,
\[
 K_d=\operatorname{im}(G)_d\qquad(d\geq4).             \tag{5.9}
\]

For the first noncentral layer this surjectivity is also made explicit.
The 24 fixed quartic-kernel tensors \(\psi_1,\ldots,\psi_{24}\) have
coefficient rank 24 and span \(K_4\).  The degree-one gauge action has rank
24 and kernel dimension three.  Exact rational row reduction constructs a
linear-polynomial matrix
\[
 Q\in\operatorname{Mat}_{9\times24}(A_1)
 \quad\text{with}\quad
 GQ=[\psi_1\ \cdots\ \psi_{24}].                     \tag{5.9a}
\]
The generated certificate stores all entries of \(Q\), not only its rank
or a modular fingerprint.

Now complete \(R=S[x,y,z]\) along
\(\mathfrak m=(x,y,z)\).  Suppose inductively that a compatible tensor has
the form
\[
 \phi_h+\theta_d+\text{terms of collision degree \(>d\)}
 \qquad(d\geq4).
\]
By (5.9), choose a homogeneous matrix \(D_{d-3}\) with
\(\delta_h(D_{d-3})=-\theta_d\).  Put
\[
 (v_z,-v_y,v_x)^{\mathsf T}=D_{d-3}r.
\]
The coordinate change \(x\mapsto x+v_x\),
\(y\mapsto y+v_y\), \(z\mapsto z+v_z\), together with
\[
 e\longmapsto (I+D_{d-3})e,
\]
preserves the Koszul relation because the exact identity
\[
 (I+D_{d-3})r=r(x+v_x,y+v_y,z+v_z)                    \tag{5.10}
\]
holds.  Formula (5.2) says that this change kills \(\theta_d\), changes no
lower jet, and only creates terms of higher collision degree.  Iteration
converges in the \(\mathfrak m\)-adic topology.  It is valid over \(S\):
all right inverses used in (5.9) are rational linear maps, so the
successive coefficients remain polynomials in the universal parameters.
Functoriality of the generalized triple-cover construction therefore gives
a formal equivalence
\[
 \widehat B_{\mathrm{univ}}
 \simeq \sigma^*(B_0\widehat{\otimes}_{\mathbb Q}S)    \tag{5.11}
\]
for an \(S\)-linear formal coordinate automorphism \(\sigma\) preserving
\(\mathfrak m\).

The central smooth cotangent presentation is exactly saturated, as checked
independently by
[`verify_cubic_symbol_double_saturation.py`](../scripts/verify_cubic_symbol_double_saturation.py).
Flat scalar extension and (5.11) give
\[
 H^0_{\mathfrak m}(
 \Omega_{\widehat B_{\mathrm{univ}}/\widehat R})=0.
\]
Completion detects every \(\mathfrak m\)-power-torsion element: if
\(\mathfrak m^n q=0\), then
\((R/\mathfrak m^n)\otimes_R\widehat R=R/\mathfrak m^n\).
Hence
\[
 \boxed{H^0_{\mathfrak m}(\Omega_{B/R})=0.}            \tag{5.12}
\]

The exact matrices and (5.3), (5.5), and (5.6) are checked by
[`verify_universal_cubic_cotangent_saturation.py`](../scripts/verify_universal_cubic_cotangent_saturation.py).

### 5.1 Exact boundary of the formal-rigidity argument

The preceding proof is special to the smooth cubic symbol.  For a cubic
orbit representative \(h\), let \(G_h\) be the determinant-twisted gauge
differential obtained from (5.2), and put
\[
 Q_h=\ker C/\operatorname{im}G_h.
\tag{5.13}
\]
The same dual-number derivation of the finite action applies to every
orbit.  Exact graded module calculation gives the following complete
atlas:

\[
\begin{array}{c|c|c|c}
h&\operatorname{Hilb}_{Q_h}(t)&\operatorname{Ann}(Q_h)
&\dim_{\mathbb Q}(Q_h)_4\\ \hline
\text{smooth}&t^3&(x,y,z)&0\\
\text{nodal}&t^3/(1-t)^2&(x)&2\\
\text{cuspidal}&2t^3/(1-t)^2&(x^2)&4\\
\text{line + transverse conic}&2t^3/(1-t)^2&(yz)&4\\
\text{line + tangent conic}&3t^3/(1-t)^2&(y^3)&6\\
\text{triangle}&3t^3/(1-t)^2&(xyz)&6\\
\text{concurrent lines}&4t^3/(1-t)^2&(x^3)&8\\
\text{double line}&t^3(5-4t)/(1-t)^3&(0)&11\\
\text{triple line}&t^3(7-5t)/(1-t)^3&(0)&16\\
0&t^3(10-6t)/(1-t)^3&(0)&24.
\end{array}
\tag{5.14}
\]

In all ten rows, \(\ker C\) has ten homogeneous generators in collision
degree three.  The 24 fixed quartic tensors span the full compatible space
\((\ker C)_4\), so the last column records exact essential quartic moduli
for this gauge action.  Thus the smooth orbit is uniquely formally rigid
above degree three.  The six singular squarefree quotients are supported
on the displayed planes, unions of planes, and nonreduced plane
thickenings.  The last three quotients are faithful over \(A\), with
generic ranks one, two, and four, respectively.

For singular symbols, (5.14) is not a cotangent-saturation failure.  It
only proves that the deformation cannot be removed by the smooth proof's
coordinate/module gauge.  Indeed, the previously verified singular
parameter planes remain cotangent-saturated despite having nongauge tensor
directions.  At this stage their universal problem therefore had to retain
the deformation-dependent cotangent complex; formal reduction to the central
tensor is unavailable.  Sections 5.2 onward close the complete squarefree
formal-tail saturation and different-generator questions by a different
normal-form and strict-Rees argument.
The exact matrices and Hilbert numerators are checked by
[`verify_cubic_formal_gauge_cokernel_atlas.py`](../scripts/verify_cubic_formal_gauge_cokernel_atlas.py).

### 5.2 The nodal first formal slice

For the nodal symbol
\[
 h_{\mathrm{nod}}=Y^2Z-X^2(X+Z),
\]
the atlas can be sharpened from Hilbert data to an explicit cyclic
presentation.  Let \(\eta_{\mathrm{nod}}\) be the compatible cubic tensor
attached to \(Z^3\).  Exact module reduction gives
\[
 \ker C=\operatorname{im}G_{\mathrm{nod}}+A\eta_{\mathrm{nod}},
 \qquad
 \operatorname{Ann}_A\!\left(
 \ker C/\operatorname{im}G_{\mathrm{nod}}
 \right)=(x).
\tag{5.15}
\]
An explicit linear-polynomial column \(L_x\) satisfies
\[
 G_{\mathrm{nod}}L_x=x\eta_{\mathrm{nod}}.
\tag{5.16}
\]
Together with (5.14), this proves
\[
 \boxed{\quad
 \ker C/\operatorname{im}G_{\mathrm{nod}}
 \simeq A/(x)(-3)=\mathbb Q[y,z](-3).
 \quad}
\tag{5.17}
\]

In collision degree four, therefore,
\[
 (\ker C)_4
 =
 \operatorname{im}(G_{\mathrm{nod}})_4
 \oplus
 \mathbb Q\,y\eta_{\mathrm{nod}}
 \oplus
 \mathbb Q\,z\eta_{\mathrm{nod}},
\tag{5.18}
\]
with dimensions \(24=22+2\).  In the fixed primitive quartic basis
\(\psi_1,\ldots,\psi_{24}\), exact rational reduction gives
\[
 [\psi_1]=z\eta_{\mathrm{nod}},\qquad
 [\psi_2]=y\eta_{\mathrm{nod}},\qquad
 [\psi_i]=0\quad(3\le i\le24).
\tag{5.19}
\]
The certificate stores a \(9\)-by-\(24\) linear-polynomial matrix giving
all 22 gauge lifts in (5.19).

The previously checked full-support directions are also a transverse
slice.  For
\[
 \psi_+=\sum_i\psi_i,\qquad
 \psi_-=\sum_i(-1)^{i-1}\psi_i,
\]
one has
\[
 [\psi_+]=(y+z)\eta_{\mathrm{nod}},\qquad
 [\psi_-]=(-y+z)\eta_{\mathrm{nod}}.
\tag{5.20}
\]
The corresponding change-of-slice matrix is
\[
 \begin{pmatrix}1&-1\\1&1\end{pmatrix},
 \qquad\det=2.
\tag{5.21}
\]
Both the coordinate slice
\(\phi_{\mathrm{nod}}+u\psi_1+v\psi_2\) and the dense slice
\(\phi_{\mathrm{nod}}+u\psi_++v\psi_-\) have exact saturated cotangent
presentations and the central length-six Ext block.  The full-complement
certificate further gives a six-dimensional Nakayama quotient for the
Kähler different on every fiber, so the nodal quartic family is nowhere
Cartier-different at its collision.  This statement is quartic-order only.

<!-- status-consumer: KDSQ6 cd423f625f1f3cd2 -->

The non-Cartier calculation now continues through the next two complete
normal-form quotients.  In degrees five and six the compatible nodal spaces
decompose as

\[
 42=39+3,
 \qquad
 64=60+4,
\]

and the complementary spaces are

\[
 \langle y^2\eta,yz\eta,z^2\eta\rangle,
 \qquad
 \langle y^3\eta,y^2z\eta,yz^2\eta,z^3\eta\rangle.
\]

On the exact nine-parameter family formed with the two quartic slice
directions, strict Rees packets commute the annihilator with every fiber and
give

\[
 J/\mathfrak nJ\simeq
 \mathbb Q[p_0,\ldots,p_8]^{\oplus6}.
\]

Thus the Kähler different is non-Cartier throughout the complete nodal
order-six normal-form family.  This finite-jet calculation is retained as
an independent regression for the all-orders result below.

<!-- status-consumer: NSDP6 c5f68253995b7b6a -->

The cyclic presentation (5.17) also closes the formal tail.  Successive
homogeneous determinant-twisted gauge transformations put every compatible
tensor with fixed nodal leading symbol in the form

\[
 \phi_{\mathrm{nod}}+f(y,z)\eta,
 \qquad f\in(y,z)\mathbb Q[[y,z]].
\tag{5.21a}
\]

For the universal coefficient family `phi_nod+u*eta`, the multiplication
table is affine-linear in `u` and exact annihilator reduction gives

\[
 J/\mathfrak nJ\simeq\mathbb Q[u]^{\oplus6}.
\tag{5.21b}
\]

Assigning `u` collision weight one, the weighted Rees presentations of
`Omega` and `coker(B -> Omega^3)` are strict and have the central
presentations tensored with `Q[u]` as associated graded.  For every graph
`u -> f`, the initial form of `u-f` is monic in `u`, hence regular on both
packets.  The annihilator therefore commutes with the graph specialization,
and (5.21b) gives six minimal generators after every formal correction.
Thus the nodal Kähler different is non-Cartier to all formal orders.  This
does not prove normality or compatibility with a Keller open.

<!-- status-consumer: NADPALL 60218641ccdf6fac -->

The same construction now covers every singular squarefree symbol, not
only the cyclic nodal row.  Minimal degree-three tensor generators give
normal-coefficient counts `1,2,2,3,3,4` for nodal, cuspidal, transverse
line--conic, tangent line--conic, triangle, and concurrent lines.  The exact
small presentations retain annihilators
`(x),(x^2),(yz),(y^3),(xyz),(x^3)`.  On each universal coefficient family,
`J/nJ` is free of rank six and the weight-one Rees packets for `Omega` and
`coker(B -> Omega^3)` are strict with central initial modules.  Successive
monic graph equations therefore preserve the intrinsic annihilator and its
six generators after every compatible formal tail.  Thus all six singular
Kähler differents are non-Cartier to all formal orders.  The smooth quotient
has no positive-order normal coefficients, and its central different also
has six minimal generators.  Consequently every squarefree cubic formal
collision is non-Cartier.

<!-- status-consumer: SSADPALL 584a6e05374612ee -->

The next collision layer is also computable.  Fix the deterministic
rational lift in (5.19) obtained by setting all free row-reduction
variables to zero, apply its inverse first-order gauge, and project the
created degree-five term to
\[
 (Q_{\mathrm{nod}})_5
 =
 \left\langle
 y^2\eta_{\mathrm{nod}},\
 yz\eta_{\mathrm{nod}},\
 z^2\eta_{\mathrm{nod}}
 \right\rangle.
\tag{5.22}
\]
The resulting normal remainder
\[
 \kappa_5^{\mathrm{row}}(u_1,\ldots,u_{24})
 =
 (\kappa_{y^2},\kappa_{yz},\kappa_{z^2})
\tag{5.23}
\]
is an exact quadratic map.  Its three components contain respectively
14, 16, and 13 monomials, with 30 nonzero cross-parameter pairs.  All
coefficients are serialized in the certificate.  It vanishes identically
on the coordinate slice:
\[
 \kappa_5^{\mathrm{row}}(u_1,u_2,0,\ldots,0)=0.
\tag{5.24}
\]
On the dense plane \(u_i=a+(-1)^{i-1}b\), however, it is
\[
 \left(
 \frac{7(4a^2-13ab+b^2)}3,\
 -\frac{25a^2-193b^2}{12},\
 -\frac{55a^2-14ab+97b^2}{6}
 \right).
\tag{5.25}
\]
Thus a transverse saturated polynomial slice can acquire a nonzero higher
normal tail when compared with another slice.  This is not a cotangent
obstruction: the dense plane itself is saturated.

The dependence on the chosen quartic gauge lift can now be removed
exactly.  Let \(W\) be the kernel of the degree-four gauge-action matrix.
Then
\[
 \dim_{\mathbb Q}W=5.
\tag{5.26}
\]
The degree-five projections of every \(W\)-\(W\) pair and every pair of
\(W\) with one of the 22 removable directions vanish.  Only the pairings
with \(\psi _1,\psi _2\) remain.  In the ordered target
\((Q_{\mathrm{nod}})_5\oplus(Q_{\mathrm{nod}})_5\), their action is
\[
 W\longrightarrow\mathbb Q^6,\qquad
 \begin{pmatrix}
 1/3&0&0&0&1/2\\
 0&-7/3&0&4&0\\
 2&0&0&0&7/2\\
 0&-2&0&7/2&0\\
 7/3&0&0&0&4\\
 0&-1/3&0&1/2&0
 \end{pmatrix}.
\tag{5.27}
\]
This matrix has rank four.  Its cokernel is detected by the two rows
\[
 (-1,0,-1,0,1,0),\qquad (0,-1,0,1,0,1).
\tag{5.28}
\]
Consequently the slice--gauge part of the degree-five curvature has the
two lift-independent quotient coordinates
\[
 L_1=\frac34(u_3+2u_9+2u_{11}),\qquad
 L_2=\frac34(3u_6+2u_{10}+2u_{12}).
\tag{5.29}
\]
The pure removable-gauge curvature is already fixed by every lift change.
Its three components are
\[
\begin{aligned}
 Q_1&=-\frac14\big((u_3-3u_5)^2-9u_6^2\big),\\
 Q_2&=\frac32(u_3u_4+u_3u_6-3u_4u_5),\\
 Q_3&=\frac14(u_3^2-9u_4^2).
\end{aligned}
\tag{5.30}
\]
Thus (5.29)--(5.30), rather than the row-reduced representative (5.23),
are the intrinsic degree-five data.

There is also an exact scheme-theoretic description of the pure-curvature
zero locus in \(\mathbb A^4_{u_3,u_4,u_5,u_6}\).  Up to invertible scalar
factors its ideal is
\[
 I_5=\big(
 (u_3-3u_5)^2-9u_6^2,\
 u_3u_4+u_3u_6-3u_4u_5,\
 u_3^2-9u_4^2
 \big).
\tag{5.31}
\]
Its reduction is the union of the two rational planes
\[
\begin{aligned}
 P_+&=(u_3-3u_4,\ u_3-3u_5+3u_6),\\
 P_-&=(u_3+3u_4,\ u_3-3u_5-3u_6).
\end{aligned}
\tag{5.32}
\]
Moreover,
\[
 \sqrt{I_5}=I_5+(g),\qquad
 g=u_3u_5-3u_5^2+3u_4u_6+3u_6^2,
\tag{5.33}
\]
and \(g\notin I_5\), while
\((u_3,u_4,u_5,u_6)g\subset I_5\).  Hence the unreduced locus has exactly
one additional degree-two socle class supported at the origin.

This closes the quartic-lift ambiguity at degree five, but it is not the
universal nodal saturation theorem.  A gauge removing the other 22
quartic coordinates creates higher collision terms, and by (5.17) the
degree-\(d\) nongauge quotient has dimension \(d-2\).

The two reduced planes in (5.32) can now be continued one step further.
For \(\epsilon\in\{1,-1\}\), parameterize them by
\[
 (u_3,u_4,u_5,u_6)
 =
 (3\epsilon p,\ p,\ \epsilon(p+q),\ q).
\tag{5.34}
\]
Keep the row-reduced quartic lift used in (5.23), apply its exact finite
inverse gauge, and write the transformed tensor as
\[
 h_{\mathrm{nod}}+t^2R_5+t^3R_6+O(t^4).
\tag{5.35}
\]
The linear term vanishes identically.  On both branches, \(R_5\) lies in
the degree-five gauge image.  That gauge-action matrix has rank \(39\)
and a \(15\)-dimensional kernel.  An explicit quadratic matrix
\(E_\epsilon(p,q)\) in the certificate satisfies
\[
 G_{\mathrm{nod}}E_\epsilon=R_5.
\tag{5.36}
\]
Moreover every change of \(E_\epsilon\) by this \(15\)-dimensional kernel
has zero image in \((Q_{\mathrm{nod}})_6\).  Thus the next class is
independent of the degree-five correction.

In the ordered quotient basis
\[
 \langle
 y^3\eta,\ y^2z\eta,\ yz^2\eta,\ z^3\eta
 \rangle,
\]
exact expansion gives
\[
 [R_6-\delta_{\psi_\epsilon}(E_\epsilon)]
 =
 \frac{27}{8}
 \left(
 q^3,\ 3\epsilon pq^2,\ 3p^2q,\ \epsilon p^3
 \right).
\tag{5.37}
\]
Equivalently,
\[
 \boxed{\quad
 \kappa_{6,\epsilon}^{\mathrm{row}}(p,q)
 =
 \frac{27}{8}(qy+\epsilon pz)^3\eta.
 \quad}
\tag{5.38}
\]
Hence, for the declared quartic splitting, the common zero locus on each
reduced plane is only \(p=q=0\).

Equation (5.38) is independent from the quadratic correction
\(E_\epsilon\), but it is **not** independent from the earlier
five-dimensional quartic gauge lift.  On the plus branch, add
$pK_1$ to the stored quartic lift, where $K_1$ is the first recorded
basis vector of that kernel.  The quartic perturbation is unchanged, while
the exact degree-six quotient becomes

\[
 -\frac98
 \begin{pmatrix}
 q^2(4p-3q) & pq(8p-5q) & p^2(4p-q) & p^3
 \end{pmatrix}^{\mathsf T}.
\tag{5.39}
\]

This differs literally from (5.37).  Its zero scheme is nevertheless still
the origin: the last coordinate forces $p=0$, and the first then forces
$q=0$.  Thus the proposed literal lift-invariance statement is false and
is replaced by a concrete stabilizer-action problem.  The correct intrinsic
object is the orbit, or an invariant quotient, of the degree-six class under
the five-dimensional quartic stabilizer.

The calculation does not yet classify that full action, continue the
embedded socle class (5.33), or treat the full slice--gauge curvature locus.
Transversality at degree four alone does not identify the full universal
family with either two-parameter slice.  Exact computation 1.8f of the
frontend now proves cotangent saturation on every geometric fiber of the
full nodal quartic complement, so the former universal quartic cotangent
gap is closed; higher formal orders, normality, and Keller-open compatibility
remain open.  The exact certificate is
[`verify_nodal_cubic_formal_slice.py`](../scripts/verify_nodal_cubic_formal_slice.py).

## 6. Canonical-different complex and the Fitting reduction

Write the trace-free and scalar multiplication components as
\(\mu_{ij}\in M\) and \(s_{ij}\in R\).  Over the full 24-parameter ring,
form the seven-column canonical-different matrix

\[
 d_1=\left[
 (0,z,-y,x)^{\mathsf T},
 (s_{ij},2\mu_{ij})^{\mathsf T}_{0\leq i\leq j\leq2}
 \right].                                             \tag{6.1}
\]

The exact checker constructs a universal \(7\)-by-\(3\) matrix \(d_2\).
For \(r=(z,-y,x)\), associativity and the coefficient-module relation give

\[
 \sum_i r_i s_{ij}=0,\qquad
 2\sum_i r_i\mu_{ij}=q_jr.                            \tag{6.2}
\]

The top row of \(d_2\) is \((-q_0,-q_1,-q_2)\); its other six rows are the
fixed incidence coefficients expressing
\(\sum_i r_i(s_{ij},2\mu_{ij})\).  Thus \(d_1d_2=0\) identically.  Each
\(q_j\) has only parameter/collision bidegrees \((0,2)\) and \((1,3)\).
The fixed lower \(6\)-by-\(3\) block has maximal minors containing
\(x^3,y^3,z^3\), while \(d_1\) has a nonzero central \(4\)-by-\(4\) minor
for every squarefree symbol.  The Buchsbaum--Eisenbud grade criterion
therefore proves exactness of

\[
 0\longrightarrow R^3\mathop{\longrightarrow}^{d_2}R^7
 \mathop{\longrightarrow}^{d_1}R^4.                  \tag{6.3}
\]

Let \(T_\Delta=\operatorname{coker}(d_1)\).  Transposing the last
differential computes \(\operatorname{Ext}^2_R(T_\Delta,R)\).  The six
fixed linear rows alone generate a module \(L\subset R^3\) satisfying

\[
 \dim_{\mathbb Q}(R^3/L)=6,\qquad
 (x,y,z)^2(R^3/L)=0.                                  \tag{6.4}
\]

The varying top row lies in \((x,y,z)^2R^3\), so it is redundant.  Hence

\[
 \operatorname{Ext}^2_R(T_\Delta,R)
 \simeq (R^3/L)
 \simeq E_0\otimes_{\mathbb Q}
 \mathbb Q[u_1,\ldots,u_{24}].                        \tag{6.5}
\]

After truncation by \((x,y,z)^2\), the parameter module has 12 generators
and six independent constant relations.  It is free of rank six, proving
universally for the canonical-different support

\[
 \operatorname{Fitt}_6=(1),\qquad
 \operatorname{Fitt}_5=(0).                           \tag{6.6}
\]

The remaining issue is now an identification, not a Fitting-minor
calculation.  The ramification support in (1.2) uses
\[
 T=B/\operatorname{Ann}_B(\Omega_{B/R}),
\]
whereas (6.1) defines \(T_\Delta\) from the seven canonical different
generators.  On the full-support plane of Section 2, exact module reduction
proves that these seven generators span the complete annihilator for every
parameter and every squarefree symbol.  The resulting actual minimal
support resolution has tail

\[
 R^3\mathop{\longrightarrow}^{d_2}R^7
 \longrightarrow R^4\longrightarrow T\longrightarrow0.       \tag{6.7}
\]

After the exact minimal-resolution basis chosen by Singular, rows two
through seven of \(d_2\) are parameter-independent linear triples.  Let
\(L\subset R^3\) be the module they generate.  The checker verifies

\[
 \dim_{\mathbb Q}(R^3/L)=6,\qquad
 (x,y,z)^2(R^3/L)=0.                                  \tag{6.8}
\]

The remaining row lies in \((x,y,z)^2R^3\).  Its central part is quadratic,
its parameter-dependent part lies in \((x,y,z)^3R^3\), and that part is
linear in the two plane parameters.  Equation (6.2) therefore makes the
entire seventh row redundant:

\[
 \operatorname{coker}(d_2^{\mathsf T})
 \simeq R^3/L.                                        \tag{6.9}
\]

Thus (6.6) already closes the universal Fittings for \(T_\Delta\), and
(6.9) checks their identification with the requested Fittings on seven
full-support planes.  To close them globally for the actual \(T\), it is
enough to prove the universal annihilator--different equality

\[
 \operatorname{Ann}_B(\Omega_{B/R})
 =
 \left((0,z,-y,x),(s_{ij},2\mu_{ij})\right).          \tag{6.10}
\]

Equivalently, one must exclude additional annihilator generators supported
over a proper parameter locus.  No further determinant calculation is
needed once (6.10) is established.

There is now an exact conditional closure of this last step.  The universal
Deligne--Faddeev cubic algebra on a free trace-free rank-two module has

\[
 \operatorname{Fitt}_0^B(\Omega_{B/R})
 =
 \operatorname{Ann}_B(\Omega_{B/R}).                  \tag{6.11}
\]

The checker computes both ideals in the universal algebra over
\(\mathbb Q[a,b,c,d]\) and reduces them to the same three-generator ideal.
On \(D(x)\cup D(y)\cup D(z)\), the Koszul trace-free module is locally free,
so (6.11) identifies the canonical different in (6.1) with the actual
annihilator.

The exact complex (6.3) also shows that \(T_\Delta\) has projective
dimension at most two.  If a prime contains \((x,y,z)\), its height is at
least three; Auslander--Buchsbaum therefore gives depth at least one for
\((T_\Delta)_\mathfrak p\).  Hence

\[
 H^0_{(x,y,z)}(T_\Delta)=0.                           \tag{6.12}
\]

Assume now the requested universal relative cotangent saturation
\(H^0_{(x,y,z)}(\Omega_{B/R})=0\).  The canonical different annihilates
\(\Omega\) off the collision axis, so its action on \(\Omega\) is
\((x,y,z)\)-torsion and therefore vanishes globally.  Thus
\(\Delta\subseteq\operatorname{Ann}(\Omega)\).  Their quotient is supported
on the collision axis by (6.11) and injects into \(T_\Delta\); (6.12)
forces the quotient to vanish.  Consequently

\[
 H^0_{(x,y,z)}(\Omega)=0
 \quad\Longrightarrow\quad
 T=T_\Delta,\quad
 \operatorname{Fitt}_6(E)=(1),\quad
 \operatorname{Fitt}_5(E)=0.                         \tag{6.13}
\]

Thus there is no independent Ext-Fitting exceptional set: it is contained
in the cotangent-saturation failure locus.  Closing universal cotangent
saturation closes the requested Fittings at the same time.

The checker
[`verify_cubic_quartic_ext_tail_absorption.py`](../scripts/verify_cubic_quartic_ext_tail_absorption.py)
proves the actual-support statement on the seven full-support planes.
The checker
[`verify_universal_cubic_quartic_different_complex.py`](../scripts/verify_universal_cubic_quartic_different_complex.py)
proves (6.2)--(6.6) over all 24 parameters.  The checker
[`verify_universal_cubic_kahler_annihilator.py`](../scripts/verify_universal_cubic_kahler_annihilator.py)
proves (6.11).  The saturation theorem (5.12) makes (6.13)
unconditional, so
\[
 \operatorname{Ann}_B(\Omega_{B/R})=\Delta,\qquad
 \operatorname{Fitt}_6^S(E)=(1),\qquad
 \operatorname{Fitt}_5^S(E)=(0).                     \tag{6.14}
\]

### 6.1 The intrinsic localized six-charge

The repeated coefficient six in (6.4)--(6.9) is not a cubic-orbit
coincidence.  It is the second Chern number of the universal quotient on the
exceptional plane, equivalently the length of the transpose vertex defect of
the same tautological presentation.

Let \(V\) be the three-dimensional collision tangent space, put
\(E=\mathbb P(V)\), and write

\[
 0\longrightarrow\mathcal O_E(-1)
 \longrightarrow V\otimes\mathcal O_E
 \longrightarrow\mathcal Q\longrightarrow0           \tag{6.15}
\]

for the universal rank-two quotient.  Up to the orientation and harmless
diagonal conventions used for \(r=(z,-y,x)\), the fixed lower \(6\)-by-
\(3\) block of \(d_2\) is the matrix of tautological symmetric
multiplication

\[
 \sigma:\mathcal O_E(-1)\otimes V
 \longrightarrow \operatorname{Sym}^2(V)\otimes\mathcal O_E.
                                                               \tag{6.16}
\]

Indeed, the row indexed by \(i\leq j\) is the relation
\(r_i e_j+r_j e_i\), with the diagonal factor absorbed into the chosen
basis.  Fiberwise the image is the degree-two part of the ideal generated
by the tautological line.  Therefore (6.16) is injective and

\[
 \operatorname{coker}(\sigma)\simeq
 \operatorname{Sym}^2(\mathcal Q).                    \tag{6.17}
\]

Put \(H=c_1(\mathcal O_E(1))\).  From (6.15),

\[
 c_1(\mathcal Q)=H,\qquad c_2(\mathcal Q)=H^2.
\]

If \(a,b\) are the formal Chern roots of \(\mathcal Q\), the roots of its
symmetric square are \(2a,a+b,2b\).  Hence

\[
 \boxed{
 c_2\!\left(\operatorname{Sym}^2\mathcal Q\right)
 =2c_1(\mathcal Q)^2+4c_2(\mathcal Q)
 =6H^2.}                                               \tag{6.18}
\]

The affine transpose gives the same integer without intersection theory.
Let \(R_0=\operatorname{Sym}(V^\vee)\), let \(\mathfrak m=R_{0,+}\), and
let \(D\) be the fixed matrix in (6.16).  Then

\[
 E_0=\operatorname{coker}
 \left(D^{\mathsf T}:R_0(-1)^6\longrightarrow R_0^3\right)
                                                               \tag{6.19}
\]

has three generators in degree zero.  Its six linear relations kill the
symmetric part of \(V^\vee\otimes V^\vee\), so the degree-one layer is
\(\bigwedge^2V^\vee\), of dimension three.  Commutativity of \(R_0\) and
alternation of this first layer give
\(\mathfrak m^2E_0=0\).  Thus

\[
 \operatorname{gr}_{\mathfrak m}E_0
 \simeq V^\vee\oplus\bigwedge^2V^\vee(1),\qquad
 \boxed{\operatorname{length}_{R_0}(E_0)=3+3=6.}       \tag{6.20}
\]

Equations (6.18) and (6.20) are the projective and vertex forms of one
localized charge.  In \(K_0^{\{\mathfrak m\}}(R_0)\simeq\mathbb Z\), the
class is \(6[k]\); its localized top Chern character is the zero-cycle
\(6[0]\).  The relevant complex is Buchsbaum--Rim-type, but the integer in
(6.20) is a colength/localized-Chern charge.  It should not be silently
identified with the classical Buchsbaum--Rim multiplicity, which is defined
from the leading coefficient of a different Rees Hilbert polynomial.

Now return to the canonical-different complex.  Its only nonfixed row is in
\(\mathfrak m^2R_0^3\), so it maps to zero in \(E_0\) by (6.20).  Whenever
the Buchsbaum--Eisenbud grade conditions identify (6.3) and cotangent
saturation identifies the canonical different with the actual annihilator,

\[
 \boxed{
 \operatorname{Ext}_{R_0}^2(T,R_0)\simeq E_0,
 \qquad [\operatorname{Ext}_{R_0}^2(T,R_0)]_{\rm loc}=6[0].}
                                                               \tag{6.21}
\]

Consequently the length-six support-hull defect for a squarefree leading
symbol is forced by the primitive Koszul conormal and is insensitive to the
cubic orbit and to every compatible higher correction for which the current
cotangent-saturation theorem applies.  No orbit-by-orbit Ext calculation is
needed for the coefficient six.

There is a parallel first-Chern statement.  The discriminant of the binary
cubic on \(\mathcal Q\) is a section of
\((\det\mathcal Q)^{\otimes6}\simeq\mathcal O_E(6)\).  If the leading
ternary cubic is squarefree, this section is nonzero, so its divisor has
class \(6H\).  Thus the branch multiplicity six and the support-defect
length six are respectively

\[
 c_1\!\left((\det\mathcal Q)^{\otimes6}\right)=6H,
 \qquad
 c_2\!\left(\operatorname{Sym}^2\mathcal Q\right)=6H^2. \tag{6.22}
\]

This explains the persistent six intrinsically but does not make it vanish.
To classify boundary-minimal cubics one still needs a global
boundary-intersection or localized-Chern conservation theorem forcing the
local class (6.21) to be zero (equivalently \(T=T^{[2]}\)), followed by the
coefficient base-change rigidity theorem.

## 7. Superseded universal elimination bottleneck

The complete smooth universal input has 27 variables and approximately
0.5 MB of exact Singular source.  Two direct routes were tested:

- a polynomial ring in all 27 variables;
- the rational function field
  \(\mathbb Q(u_1,\ldots,u_{24})[x,y,z]\).

Both bottleneck before saturation or Ext, while standardizing the initial
31-relation cotangent module.  The polynomial-ring route did not reach its
first invariant after ten minutes.  The rational-function-field route
exceeded 4 GB resident memory before the same point.  A block order with
\((x,y,z)\) first reduced memory but did not remove the elimination
bottleneck.  On the smooth nested coordinate family, ten parameters
complete.  The corresponding first-eleven-parameter calculation did not
complete within a fresh 300-second bound; an earlier first-twelve-parameter
calculation did not complete within 900 seconds and reached approximately
1.27 GB resident memory.  These failed computational routes are superseded
by the formal-gauge proof in Section 5.

Unit pruning reduces the universal input to the \(6\)-by-\(25\) matrix of
Section 4 and approximately 0.2 MB of Singular source.  With
\((x,y,z)\) placed in the first elimination block, this reduced calculation
still did not return its first standard-basis diagnostic after six minutes
and reached approximately 1.75 GB resident memory.  This is another failed
direct route, not an exceptional parameter.

The canonical-different construction of Section 6 completes the
parameter-only matrix that this bottleneck originally suggested.  Its six
linear relation columns act on the nine-dimensional space

\[
 (x,y,z)\,S^3/(x,y,z)^2S^3.                         \tag{7.1}
\]

The Fittings of this matrix are closed by (6.6).  Before the formal-gauge
theorem, the next proposed certificate was the equality (6.10), formulated
as the finite quotient

\[
 \operatorname{Ann}_B(\Omega_{B/R})/\Delta,
 \qquad
 \Delta=((0,z,-y,x),(s_{ij},2\mu_{ij})),             \tag{7.2}
\]

It vanishes on every recorded line, plane, and the smooth coordinate
ten-space.  The theorem (5.12) and implication (6.13) now prove that it is
zero universally; no 24-parameter elimination remains.

## 8. Reproduction

Run

```bash
.venv/bin/python scripts/research_universal_cubic_quartic_kernel_saturation.py
.venv/bin/python scripts/verify_cubic_quartic_ext_tail_absorption.py
.venv/bin/python scripts/verify_universal_cubic_quartic_different_complex.py
.venv/bin/python scripts/verify_universal_cubic_kahler_annihilator.py
.venv/bin/python scripts/verify_universal_cubic_cotangent_saturation.py
.venv/bin/python scripts/verify_cubic_formal_gauge_cokernel_atlas.py
.venv/bin/python scripts/verify_nodal_cubic_formal_slice.py
.venv/bin/python scripts/verify_cubic_double_saturation_stratification.py
.venv/bin/python scripts/verify_nodal_sextic_different_persistence.py
.venv/bin/python scripts/verify_nodal_all_orders_different_persistence.py
.venv/bin/python scripts/verify_singular_squarefree_all_orders_different_persistence.py
```

The generated records are
[`universal_cubic_quartic_kernel_saturation_frontier.json`](../artifacts/generated-results/universal_cubic_quartic_kernel_saturation_frontier.json)
and
[`universal_cubic_cotangent_saturation.json`](../artifacts/generated-results/universal_cubic_cotangent_saturation.json),
together with the all-orbit boundary calculation
[`cubic_formal_gauge_cokernel_atlas.json`](../artifacts/generated-results/cubic_formal_gauge_cokernel_atlas.json)
and the first singular slice certificate
[`nodal_cubic_formal_slice.json`](../artifacts/generated-results/nodal_cubic_formal_slice.json).
The completed squarefree continuation is recorded in
[`cubic_double_saturation_stratification.json`](../artifacts/generated-results/cubic_double_saturation_stratification.json),
[`nodal_sextic_different_persistence.json`](../artifacts/generated-results/nodal_sextic_different_persistence.json),
[`nodal_all_orders_different_persistence.json`](../artifacts/generated-results/nodal_all_orders_different_persistence.json),
and
[`singular_squarefree_all_orders_different_persistence.json`](../artifacts/generated-results/singular_squarefree_all_orders_different_persistence.json).

The calculations require Singular 4.4.1.  The smooth checker uses Singular
only for the three-variable exact module identity (5.5) and the independent
central saturation check; it does not compute the 24-parameter
saturation.  The formal-gauge argument proves universal cotangent
saturation, and the canonical-different complex then identifies the actual
annihilator and closes the universal Fittings.

The canonical-different checker also verifies the tautological
symmetric-multiplication matrix, the `3+3` vertex layers, square
annihilation, and the Chern-root identity (6.18); this is the exact replay
of the localized charge in Section 6.1.

The atlas checker computes the exact three-variable graded cokernel for all
ten symbols and marks the
limit of that formal-rigidity argument.  The nodal checker resolves the
first singular row cyclically and splits the complete quartic space into
22 gauge directions and a certified saturated two-parameter slice.  The
later four checkers cover every singular-squarefree quartic nongauge
complement, the complete nodal degree-six normal form, the nodal formal
tail, and finally every squarefree compatible formal tail.  Their scopes
remain local/formal and do not certify normality or a Keller open.
