# Finite matrix differential actions and primitive admission

15 November 2026. **QRA1.** Restriction of scalars along a finite free
map of smooth affine coordinate rings gives finite matrix differential
actions without an etaleness hypothesis. For a monic one-variable map
of any degree, the resulting Weyl generators have spatial degree at
most two. Their branch correction shows that a classical coordinate
factorization need not persist in the quantum action.

For the smooth elliptic isogeny of QIM1, an explicit Euler connection
extends the curve action to the whole varying-energy surface and all
six rank-two Weyl relations. The energy is noncentral in this enlarged
operator algebra. The physical-coordinate action still requires b inverse. The monic
coordinates x,w admit a different global polynomial full action. Both
have special left rank three, rather than the rank-one Long graph.

A separate all-degree admission criterion applies at b=0. If a special
coordinate endomorphism is injective but not surjective, it cannot
become a unit in the actual generic endomorphism algebra. An exact
primitive relation Y=XW then conflicts with the appropriate Weyl
relations. For the elliptic matrices det(X_0)=-4b, so retaining that
primitive does not repair the singular-energy extension, even after
dropping the explicit x inverse generator. Quantum changes which break
the exact primitive relation remain outside this obstruction.

## Restriction of differential operators along a finite free map

Let k have characteristic zero. Suppose A=k[u_1,...,u_n] and
B=k[z_1,...,z_n] is a finite free A-algebra through a polynomial map
u=F(z). Choose an A-basis beta_1,...,beta_m of B, and use columns.
Write M_g for multiplication by g in B. For a k-derivation delta of B
write delta(beta_j)=sum_i C_delta,ij beta_i. Chain rule gives
\[
 \rho(g)=M_g,\qquad
 \rho(h\delta)=\sum_{\alpha=1}^n
             M_{\delta(F_\alpha)}p_\alpha+h C_\delta,
 \qquad [p_\alpha,u_\gamma]=h\delta_{\alpha\gamma}.          \tag{1}
\]
Acting on coefficient columns, the two sides of (1) are the same
operator on B: differentiate both the basis and its A coefficients.
Products therefore define a faithful homomorphism from the source
Weyl algebra to Mat_m(A_n(k[h])). All source Weyl relations hold.
After replacing h by any nonzero quantum parameter a(t), the same
formulas hold over the finite PBW parameter algebra.

More generally, if A and B are smooth affine coordinate rings and
B is finite free over A, a differential operator of B of order d is
an A-differential operator of order at most d, since A is a subalgebra
of B. In a free basis this embeds the differential-operator algebra
of B in Mat_m(D(A)), preserving the order filtration. The same
statement for the Rees algebras follows by multiplying an order-d
operator by h^d. Faithfulness is faithfulness of its action on B.
For the differential-operator definitions and etale extension
comparison see [Stacks00UP](https://stacks.math.columbia.edu/tag/00UP).
The direct chain-rule construction here also allows ramified maps.

This is a full matrix action of the source Weyl algebra. It is not
an endomorphism of Mat_m(A_n), and the special left module has rank m.
The ordinary rank-one graph hypothesis of QMC1 is not being asserted.

## Ramified companion Weyl pairs in every degree

Let f(z) be monic of degree m>=1 over k, and set x=f(z).
Use the basis 1,z,...,z^(m-1) of k[z] over k[x]. Let Z be its
multiplication companion matrix, and let C be the constant matrix
with C_(j-1,j)=j for 1<=j<m. Formula (1) becomes
\[
              X=Z,\qquad
              T=f'(Z)p+h C,\qquad [T,Z]=hI_m.              \tag{2}
\]
Both sides represent multiplication by z and h*d/dz on k[z].
Thus (2) is an exact identity in every m and for every monic f,
rather than a finite degree census.

The entries of Z have x-degree at most one. To obtain f'(Z),
one reduces z^j f'(z), with degree at most 2m-2, modulo f(z)-x.
Each reduction introduces at most one x in every term: the x term
then has z-degree less than m, while the remaining terms have
constant coefficients and smaller z-degree. Consequently
\[
            \deg_{\rm spatial} Z\le1,\qquad
            \deg_{\rm spatial} T\le2,                     \tag{3}
\]
independently of m. Here h is a scalar, and x,p have spatial degree one.

For f(z)=z^m one also has Z^m=xI and
deg(Z^k)=ceil(k/m). In particular the powers grow unboundedly,
as required by [QCI1](CENTRALIZER_INDEX_AND_RATIONAL_POLARIZATION.md)'s
finite-lattice power criterion. This action is not a diagonal
representation of the rational power-cotangent pair excluded there.

For m>=2 the classical momentum factors through Z, but the actual
quantum momentum does not:
\[
 T_0=Z_0\bigl(mZ_0^{m-2}p\bigr),\qquad
 T=Z\bigl(mZ^{m-2}p\bigr)+h C,\qquad
 (Z^{-1}C)_{m-1,1}=1/x.                                  \tag{4}
\]
Indices start at zero. This is a positive control for a necessary
quantum correction which releases a classical primitive from the
actual polynomial matrix algebra. Classical factorization alone
is insufficient for the obstruction below. These examples have
arbitrary special rank m; no rank-two Dixmier counterexample follows
from their uniformly small generator degrees.

## The elliptic module admits an explicit energy connection

Use the family version of
[QIM1](ELLIPTIC_ISOGENY_MATRIX_MODULE.md):
\[
 \begin{split}
 A&=k[u,v,b,b^{-1}]/(v^2-u^3+27b),\\
 B&=k[x,x^{-1},y,b,b^{-1}]/(y^2-x^3-b)
      =A[w]/(w^3-3uw/4+v/4).
 \end{split}                                             \tag{5}
\]
Its basis is 1,w,w^2, with w=y/x. The two rings are smooth affine
surfaces; A is the open plane with b=(u^3-v^2)/27, and B is the
corresponding source open plane. The map A -> B is finite etale.

Alongside the invariant derivation delta, introduce the Euler
derivation E normalized by E(b)=1:
\[
 \begin{array}{lll}
 \delta(u)=v,&\delta(v)=3u^2/2,&\delta(b)=0,\\
 E(u)=u/(3b),&E(v)=v/(2b),&E(b)=1,\\
 \widetilde E(x)=x/(3b),&
 \widetilde E(y)=y/(2b),&
 \widetilde E(w)=w/(6b).
 \end{array}                                             \tag{6}
\]
The last row is the unique lift. The QIM1 lift of delta has
delta(x)=-y, delta(y)=-3x^2/2. On both surfaces
\[
 [E,\delta]=\delta/(6b),\qquad
 \det\begin{pmatrix}\delta(u)&E(u)\\
                     \delta(v)&E(v)\end{pmatrix}
   =-27/2.                                               \tag{7}
\]
Thus delta,E form a global basis for the tangent module over A;
the analogous determinant over B is -1/2.

Let q=h*delta and P=h*E. Retain QIM1's matrices
\[
 W=\begin{pmatrix}0&0&-v/4\\1&0&3u/4\\0&1&0\end{pmatrix},
 \quad C_\delta=\begin{pmatrix}0&u/2&v/2\\0&0&-u/2\\0&-1&0\end{pmatrix},
 \quad C_E=\operatorname{diag}(0,1/(6b),1/(3b)),
 \quad Z=qI+hC_\delta,\quad P_B=PI+hC_E.                   \tag{8}
\]
The energy connection is compatible in every operator degree:
\[
 \begin{split}
 E(W)+[C_E,W]&=W/(6b),\\
 E(C_\delta)-\delta(C_E)+[C_E,C_\delta]&=C_\delta/(6b),\\
 [P_B,Z]&=hZ/(6b),\qquad [P_B,bI]=hI.
 \end{split}                                             \tag{9}
\]
The second equation is the mixed curvature condition for (7).
Together with QIM1's delta identity it gives the full source
surface differential action. Equivalently, etaleness lifts both
target tangent generators, and the basis makes D_h(B) free of
rank three on both sides over D_h(A).

Set X=(4W^2-uI)/3 and Y=(2uW-vI)/3. All four source Weyl generators
can be written explicitly:
\[
 \begin{split}
 D_x&=-(Y/b)Z-3X^2P_B,\\
 D_y&=(2X/(3b))Z+2YP_B,\\
 [D_x,X]&=[D_y,Y]=hI,\quad
 [D_x,Y]=[D_y,X]=[D_x,D_y]=[X,Y]=0.
 \end{split}                                             \tag{10}
\]
These follow from
partial_x=-(y/b)*delta-3x^2*E and
partial_y=(2x/(3b))*delta+2y*E on B. They are an actual faithful
A_2 action by three-by-three differential matrices on the open
target plane. No infinite parameter series or momentum orders
are needed. In ordinary target Weyl coordinates,
q=v*p_u+(3u^2/2)*p_v and
P=u*p_u/(3b)+v*p_v/(2b).

The monic model across b=0 is
B_0=A_0[w]/(w^3-3uw/4+v/4)=k[x,w], with
u=4w^2-3x and v=8w^3-9xw. Formula (1) gives a different global
rank-two Weyl action with commuting positions X,W:
\[
 D_{x\mid w}=-3p_uI-9Wp_v,\qquad
 D_{w\mid x}=8Wp_u+(24W^2-9X)p_v+hC_0,\qquad
 (C_0)_{j-1,j}=j.
\]
All six relations for X,W,D_(x|w),D_(w|x) hold polynomially over k[u,v].
It has special left rank three. The physical function Y=XW is a
composite in this action: [D_(w|x),Y]=hX, rather than hI.

Its logarithmic energy operator is also polynomial:
\[
 L_B=bP_B=(bP)I+h\operatorname{diag}(0,1/6,1/3),\qquad
 bP=(u/3)p_u+(v/2)p_v,\qquad
 [L_B,b]=hb,\quad [L_B,Z]=hZ/6.                            \tag{11}
\]
This gives a global model with noncentral energy and a full Weyl
action in the monic coordinates x,w. It has an Euler relation with b,
rather than a Weyl relation, and does not give an independent physical
y derivative commuting with x. The derivatives (10) still require
their b denominators. Noncentrality of energy alone therefore does
not supply the missing physical canonical action.

This makes energy noncentral in the enlarged algebra, through
[P_B,b]=h. It quantizes the cotangent geometry of the varying
surface: x,y remain commuting positions. It does not change
their commutator to h as required for the physical isogeny pair
in [QEM1](ELLIPTIC_QUANTUM_ENERGY_MOTION.md), nor does it identify
the special rank-three module with the rank-one Long graph.
The poles in b are substantive.

## Inverse poles cannot hide an injective special coordinate

Let R=K[[t]], let M be an R-flat module, and let X be an R-linear
endomorphism. If X_0 on M/tM is injective and X becomes invertible
on M[t^-1], its inverse preserves M.

For m in M, put z=X^-1 m in M[t^-1]. If z is not in M, choose the
least integer j>=1 such that t^j z lies in M. Its class modulo t
is nonzero, while X(t^j z)=t^j m has zero class. This contradicts
injectivity of X_0. Hence the inverse is integral and
\[
 X_0\text{ injective but not surjective}
       \quad\Longrightarrow\quad
 X\text{ is not a unit in }\operatorname{End}(M[t^{-1}]).
                                                               \tag{12}
\]
This is an implication from an assumed actual inverse, not a
lifting theorem from a special inverse. For example 1-t*x is
special-invertible on R[x] but has no polynomial generic inverse.
The example X=tI shows why special injectivity is required to
remove inverse poles.

## Exact primitives obstruct global admission in any finite rank

Suppose the actual generic endomorphism algebra E embeds in
Mat_r(D) for a division ring D. It is directly finite: a right
inverse in E is also a left inverse. This applies to finite
torsion-free Weyl modules and to polynomial matrix actions.
Consider X,Y,W in E with Y=XW. Either of the following full
relation patterns forces X to be a unit:
\[
 \begin{array}{ll}
 [Y,X]=aI\ne0
   &\Longrightarrow\quad X[W,X]=aI,\\
 [T,X]=0,\quad [T,Y]=aI\ne0
   &\Longrightarrow\quad X[T,W]=aI.
 \end{array}                                             \tag{13}
\]
The scalar a is inverted in the generic field. This argument
uses actual multiplication and holds in every matrix rank and
every operator degree. In a finite full Weyl carrier, generic
torsion-freeness and projectivity are supplied by QMC1's stated
external finite-left bimodule theorem.

For the monic elliptic model over
A_0=k[u,v] with b=(u^3-v^2)/27, multiplication by x has
\[
                 \det X_0=-4b.                          \tag{14}
\]
After adjoining ordinary momentum symbols, this determinant is
still nonzero and is not a unit. Thus X_0 is injective but not
surjective on the free rank-three special module. Formula (11)
excludes a generic actual inverse.

Consequently no finite polynomial Weyl-matrix extension across
b=0 can retain the exact primitive Y=XW and either full pattern
(12) with this special module. This includes arbitrary quantum
corrections to X,Y,W and their partners, all momentum orders and
all ramified quantum parameters a(t). For the surface construction
(10) the second pattern is supplied by T=D_y. For an attempted
physical isogeny Weyl pair, the first pattern applies instead.
The gate is invariant under actual module isomorphisms which
transport every operator and preserve the identified special module.

Keeping x inverse would already fail (11). Dropping x inverse
but keeping the exact primitive W=y/x does not remove the failure.
The monic model's loss of that inverse is not enough to admit the
full physical Weyl action.

This leaves a precise possible escape: the primitive relation itself
must change, or W must leave the finite endomorphism algebra.
The genuine companion pairs (4) prove that such quantum corrections
can matter. Neither an exact-energy obstruction nor the classical
identity Y_0=X_0 W_0 can exclude all corrected matrix actions.
Special modules with a different determinant or rank, and genuinely
different quantum relations, need their own admission.

## Rational branch extraction does not change the generic rank

For the finite free polynomial-coordinate construction (1), let D_A
and D_B be the target and source differential-operator division rings.
The map of function fields has degree m and is separable. Lift target
derivations to the source function field. The finite differential-field
extension theorem used in QDI1 gives
\[
 \dim_{D_A}(D_B)=m
 \quad\text{on either side},\qquad
 0\ne V'\subset D_B,\ D_B V'\subset V'
       \quad\Longrightarrow\quad V'=D_B.                 \tag{15}
\]
The matrix action is the left regular action on the right D_A-basis
beta_1,...,beta_m. A right D_A-subspace invariant under every source
Weyl generator is invariant under inverses of their nonzero polynomial
operators: these act injectively on the whole division ring, hence
bijectively on a finite-dimensional invariant subspace. It is therefore
invariant under all D_B. A nonzero left ideal of D_B is the whole ring,
proving (15).

Dually the row representation has no proper nonzero invariant left
D_A-subspace. Its full-source commutant therefore has no nontrivial
idempotents. An extension of this same source action from A_n to a
matrix algebra Mat_k(A_n), k>1, would put nonzero orthogonal matrix
idempotents in that commutant, and is impossible. Thus no rational change of target matrix frame can extract
a rank-one full source action when m>1. For the monic elliptic model
this excludes rational compression from its generic rank three to a
rank-one carrier. It does not exclude changing the action or constructing
a separate Long carrier. Extension of the target coefficient field can
change this representation: power companions split after a suitable
coefficient cover, whose global finite admission must be checked anew.

## Verification and remaining goal

Run check_matrix_primitive_admission.py. It checks companion actions
for representative monic degrees, including ramification, and their
exact chain-rule matrix identities. The construction for every degree
is proved above. It verifies
the explicit elliptic Euler connection, mixed curvature, all six
physical surface Weyl relations modulo the curve equation, and all six
global monic-coordinate Weyl relations. It also checks
the norm, an interior zero-energy point, and the nonintegral primitive
correction in the power companions.

The restriction-of-scalars construction, uniform one-variable degree
bound, inverse-pole argument and exact-primitive obstruction are written
proofs. The examples do not establish an all-rank claim by census.
The physical-coordinate surface action is localized; the global monic
action changes the source coordinates. Both have special rank three,
without a matrix-algebra endomorphism or a full rank-one Long carrier.
No bounded full Long lift, all-representative nonlifting, DC2/JC2
solution, independent theorem replay, formal verification, external
review or novelty priority is claimed.
