# Odd graded Long jets and nonlinear compatibility

The full homogeneous odd first-jet family at common spatial degree 23
has dimension six. Two directions admit exact uniformly bounded shears.
The complete nonlinear next equation forces the remaining \(H_6\)
parameter to vanish, while retaining actual second-jet points for every
choice of the three parameters \(Y,R^2T,H_4\). Their complete second
coefficients have two additional homogeneous freedoms, which are kept.
No evenization, higher compatibility or bounded full lift is proved.

## Complete first-jet problem

Use source order \((X,Q,P,Z)\), Poisson convention
\(\{P,X\}=\{Z,Q\}=1\), and source weights \((-1,1,4,2)\). Give \(h\)
weight 3. Let
\[
 f=(R,T,D,S)
\]
be the canonical polynomial Long images in
[long_images()](../dixmier-audit/verify_long_pbw.py), with image weights
\((-1,2,4,1)\). They satisfy the six full polynomial Poisson frame
relations. Work with Weyl symbols and the Moyal product.

For symbols
\[
 F_i=f_i+h g_i+O(h^2),
\]
the coefficient of \(h^2\) in each required commutator is
\[
 \{g_i,f_j\}+\{f_i,g_j\}=0,\qquad i<j.                 \tag{1}
\]
These are the homogeneous linearized Poisson equations. The nonlinear
equations at commutator order \(h^3\), including the cubic Moyal term,
are treated in the nonlinear compatibility section below.

Require every \(g_i\) to have ordinary source degree at most 23 and
weight \(w_i-3\). Writing its source monomial as \(X^aQ^bP^cZ^d\),
the complete admission condition is
\[
 a=b+4c+2d+3-w_i\ge0,\qquad a+b+c+d\le23.              \tag{2}
\]
Thus the component basis counts are
\[
              (80,114,140,102),\qquad 436\text{ in total}.       \tag{3}
\]
No positional, momentum-order or fixed-first-image restriction is
imposed.

## Six compact Hamiltonians

Put
\[
\begin{split}
 B&=3X^2P+2(1-3XQ)Z,\\
 A_0&=Q+3XQ^2-XB/3,\\
 t&=XQ,\qquad Y=RA_0-3t.
\end{split}                                                    \tag{4}
\]
These are commuting classical polynomial expressions. In the rational
Long chart \(a=A_0+X^{-1}\), the last expression is \(Y=Ra-2\),
since \(R/X=2-3t\). In particular \(Y\) is globally polynomial despite
this chart interpretation.

The complete solution space of (1)--(2) is
\[
 g_i=\{G,f_i\},\qquad
 G\in\operatorname{span}_K
       \{t,Y,Y^2,R^2T,H_4,H_6\},                               \tag{5}
\]
where
\[
 H_4=Y^4+4tY^2-\frac{16}{3}tY+10t^2,                           \tag{6}
\]
and
\[
 H_6=R^4D+
 \frac{9Y^5+45tY^3+45t^2Y-60tY^2+110tY-720t^2}{25}.            \tag{7}
\]
This holds over every characteristic-zero coefficient field \(K\).

The high source degrees of these Hamiltonians do not violate (2):
their variations have lower degrees through exact cancellations.

| Hamiltonian | Source degree | Degrees of its \(R,T,D,S\) variations |
|---|---:|---|
| \(t\) | 2 | \((-,15,23,11)\) |
| \(Y\) | 7 | \((-,8,16,4)\) |
| \(Y^2\) | 14 | \((-,15,23,11)\) |
| \(R^2T\) | 21 | \((-,-,18,6)\) |
| \(H_4\) | 28 | \((-,17,23,13)\) |
| \(H_6\) | 35 | \((12,17,23,20)\) |

Here a dash denotes a zero component. Only \(H_6\) changes the first
image, and
\[
                         \{H_6,R\}=R^4.                       \tag{8}
\]

For completeness, form the coefficient matrix of (1) in the entire
basis (2). It has 2494 nonzero-support rows and 68262 nonzero entries.
The [explicit witness](odd_graded_long_first_jets_witness.json)
selects a \(430\)-square minor whose determinant is
\[
                         288992\pmod {1000003}.                \tag{9}
\]
All denominators are invertible at this prime, so the rational rank is
at least 430. Exact substitution proves that the six variations in
(5) are kernel vectors with admitted support. A six-square coefficient
minor of these vectors, also specified in the witness, has determinant
\[
                              -10240/9.                        \tag{10}
\]
They are independent. Hence the rational rank is exactly 430 and the
kernel is precisely (5). This is a complete finite coefficient
classification, rather than a search among the displayed ansatz.

## Reconstruction by the full symplectic Jacobian

The Hamiltonian assertion can also be recovered directly from the
variations. The full polynomial symplectic Jacobian gives the source
one-form
\[
 dG=-g_D\,dR-g_S\,dT+g_R\,dD+g_T\,dS.                          \tag{11}
\]
It follows from the Poisson frame equations that this is the
Hamiltonian one-form of the corresponding source vector field.
The linearized frame equations make it closed. Its polynomial
primitive, normalized to vanish at the origin, is obtained by ordinary
Euler integration:
\[
 G(x)=\int_0^1\sum_{\ell=1}^4 x_\ell\,\omega_\ell(ux)\,du,
 \qquad \omega=\sum_\ell\omega_\ell(x)\,dx_\ell.                \tag{12}
\]
A coefficient of ordinary degree \(n\) in \(\omega_\ell\) is divided
by \(n+1\). The checker verifies (11)--(12) exactly for all six basis
vectors, recovering (5)--(7). Weighted Euler integration cannot be used
here: each \(G\) has weight zero.

## What can be removed without increasing the common cap

There is the exact classical identity
\[
                         RS=-6t-3Y^2-8Y.                      \tag{13}
\]
Consequently one may use
\(\{t,RS,Y,R^2T,H_4,H_6\}\) as a Hamiltonian basis.

The source Hamiltonian \(t=XQ\) gives the exact quantum shear
\[
 X\mapsto X,\quad Q\mapsto Q,\quad
 P\mapsto P-bhQ,\quad Z\mapsto Z-bhX.                          \tag{14}
\]
It is a linear substitution in source generators, so it preserves any
common ordinary spatial cap and the stated grading.

For any actual Weyl tuple, the target Hamiltonian \(RS\) similarly
gives
\[
 R\mapsto R,\quad S\mapsto S,\quad
 T\mapsto T+bhR,\quad D\mapsto D-bhS.                          \tag{15}
\]
This is linear in the four images, so it too preserves every common
cap and the grading. Equations (14)--(15) therefore justify removing
two first-jet parameters while remaining in the bounded problem.
A remaining first jet can be represented by the four Hamiltonians
\[
                              Y,\ R^2T,\ H_4,\ H_6.            \tag{16}
\]

The target Hamiltonian \(R^2T\) has the finite exact abstract flow
\[
 R\mapsto R,\quad T\mapsto T,\quad
 S\mapsto S-bhR^2,\quad D\mapsto D-2bhRT.                     \tag{17}
\]
Finiteness alone does not preserve the common source degree cap.
Suppose a first jet has nonzero \(H_6\)-coefficient \(c\).
Then \(R_1=cR_0^4\), whereas every admitted first variation \(T_1\)
has degree at most 17. The coefficient of \(h\) in the Weyl symbol
of \(R(h)T(h)\) is
\[
 R_1T_0+R_0T_1+\tfrac12\{R_0,T_0\}
                    =cR_0^4T_0+R_0T_1.                     \tag{18}
\]
The first term has degree 27, the second degree at most 20. Thus, if
\(bc\ne0\), (17) creates a degree-27 term in the \(h^2\) coefficient
of \(D\). No preexisting cap-23 \(D_2\) cancels it. This is a concrete
failure of uniform cap preservation, conditional only on that first
jet; it does not assert existence of a full corrected tuple.

Two natural classical Hamiltonian exponentials also show why a
first-order Hamiltonian description supplies no finite normalization:
\[
 \{Y,R\}=0,\qquad \{Y,X\}=-RX^3
 \quad\Longrightarrow\quad
 X(\tau)=\frac{X}{(1+2\tau RX^2)^{1/2}},                       \tag{19}
\]
and (8) gives
\[
                     R(\tau)=\frac{R}{(1-3\tau R^3)^{1/3}}.  \tag{20}
\]
Their formal expansions have unbounded spatial degrees in
characteristic zero. These statements concern the indicated natural
classical exponentials. They exclude no alternative nonlinear
parameter-dependent normalization with the same first derivative.

## Complete nonlinear compatibility at commutator order \(h^3\)

After only the two uniformly bounded gauges (14)--(15), write
\[
 g^{(1)}=\{aY+bR^2T+cH_4+dH_6,f\},
 \qquad F=f+h g^{(1)}+h^2g^{(2)}+O(h^3).
\]
Allow EVERY degree-23 polynomial \(g^{(2)}_i\) of weight \(w_i-6\).
Its complete component counts are \((54,80,102,71)\), totaling 307.
At commutator order \(h^3\), all six equations are
\[
 A g^{(2)}=b_0+\sum_{\mu\leq\nu}p_\mu p_\nu b_{\mu\nu},
 \quad(p_1,p_2,p_3,p_4)=(a,b,c,d),                         \tag{21}
\]
where
\[
 b_{0,ij}=-\Pi^3(f_i,f_j)/24,\qquad
 \sum_{\mu\leq\nu}p_\mu p_\nu b_{\mu\nu,ij}
                         =-\{g^{(1)}_i,g^{(1)}_j\}.
\]
The quadratic term is essential; this is not the even-jet problem.
No \(R^2T\), \(H_4\) or \(H_6\) parameter was discarded on the
basis of an infinitesimal gauge.

**Nonlinear compatibility theorem.** Over every characteristic-zero
coefficient field, (21) has an admitted solution if and only if
\(d=0\). The parameters \(a,b,c\) remain arbitrary, and actual
degree-23 polynomial second coefficients are retained for all
their values.

The row set has 2005 entries. The witness supplies a rational
functional \(\ell\), supported on 288 specified rows, satisfying
\[
 \ell A=0,\qquad \ell b_{dd}=1,\qquad
       \ell b_0=\ell b_{\mu\nu}=0
                   \quad((\mu,\nu)\ne(d,d)).              \tag{22}
\]
The replay checks (22) exactly on all 307 coefficient columns and
all eleven right-side coefficient vectors. Thus (21) forces
\(d^2=0\), hence \(d=0\) over a field. No modular-only inference
or evaluation at sampled parameter values is used.

For sufficiency, the witness retains seven exact polynomial tuples
\(P_0,P_{aa},P_{ab},P_{ac},P_{bb},P_{bc},P_{cc}\), with 346
nonzero rational coefficients in total, and verifies the corresponding
seven full matrix identities. Consequently
\[
 \Psi(a,b,c)=P_0+a^2P_{aa}+abP_{ab}+acP_{ac}
                    +b^2P_{bb}+bcP_{bc}+c^2P_{cc}         \tag{23}
\]
solves (21) for every \(a,b,c\) with \(d=0\). In fact \(P_{bb}=0\).
These are explicit finite polynomial points, not a further
formal-existence assertion.

The second matrix \(A\) has rational rank 305: the retained
305-square minor from the even first-correction certificate is
replayed, and its two exact independent kernels are
\[
 K_{R^3}=(0,0,-3R^2,0),\qquad
 K_{R^4S}=(0,R^4,-4R^3S,0).
\]
Thus the ENTIRE second-coefficient family is
\[
       g^{(2)}=\Psi(a,b,c)+uK_{R^3}+vK_{R^4S},
                           \qquad u,v\ \hbox{arbitrary}. \tag{24}
\]
The parameters \(u,v\) have not been normalized away in later jets.

Every tuple in (23), and both kernels in (24), has zero
\(R\)-component. Therefore the normalized admitted branch has
\(R_1=R_2=0\). Undoing the two allowed first gauges does not change
this: the source \(t\)-shear fixes positional \(R_0\) and
\(R_1=dR^4\), and the target \(RS\)-shear fixes \(R\) exactly.
Hence every degree-23 equivariant jet compatible through \(h^3\)
has \(R_1=R_2=0\), even before those two normalizations.
This does not say that a higher first-image correction vanishes.

The earlier degree-27 \(R^2T/H_6\) example remains a valid
first-order warning, but its nonzero \(H_6\) first jet cannot pass
(21). No uniform removal of the remaining odd directions is
thereby proved. There is also a concrete later warning for (24):
the \(c^2\) part of the retained \(S_2\) has degree 15, so
\(R^3S_2\) has a nonzero degree-24 part. Removing a nonzero
\(vK_{R^4S}\) by its natural order-\(h^2\) target shear therefore
generically adds an over-cap term to \(D_4\) when \(cv\ne0\).
The ordering contribution \(\Pi^2(R^3,S_0)/8\) has degree at most
16 and cannot cancel that part. No higher compatibility system
has been solved here.

Finally, the grading makes the unrestricted-parity cap-23 problem
finite: the coefficient of \(h^m\) has source weight \(w_i-3m\)
and is divisible by \(X^{\max(0,3m-w_i)}\).
Thus the maximal possible parameter indices for \((R,T,D,S)\)
are \((7,8,9,8)\). This is a finite polynomial-in-\(h\) coefficient
scheme in the stated grading, not a normalization of arbitrary
bounded lifts or a proof that these partial jets extend within
the same cap.

## Scope and replay

The result classifies all first odd homogeneous jets at common cap 23,
and the entire admitted second coefficient over its three surviving
parameters. The nonlinear \(h^3\) condition excludes precisely the
\(H_6\) first parameter over characteristic-zero fields, and forces
\(R_1=R_2=0\). It proves no \(h^4\) or \(h^5\) compatibility,
evenization theorem, full corrected finite lift or nongraded
nonexistence theorem. An even obstruction is not transferred to the
remaining odd family by normalization.

[The replay checker](check_odd_graded_long_first_jets.py) uses only
Python integers and Fraction, reusing
[the complete linear-column construction](check_graded_long_cap23.py)
and [the independent polynomial curvature engine](check_long_graph_curvature.py).
The canonical producer supplies input polynomials only. The
[retained receipt](odd_graded_long_first_jets_controls.json) pins the
two minors, all six exact kernel substitutions, ordinary Euler
reconstruction, shear identities and the degree-27 cap failure.
