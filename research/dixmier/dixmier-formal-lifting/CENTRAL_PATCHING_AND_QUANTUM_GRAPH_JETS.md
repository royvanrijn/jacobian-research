# Central patching and finite jet reduction of quantum graph carriers

15 November 2026. **QPG1.** For a full polynomial symplectic classical
map, the finite graph-carrier problem can be posed on the generic
Weyl algebra and the punctured formal Weyl algebra. A compatible
generic finite projective action and a formal rank-one graph glue
to an actual finite projective carrier. No spatial bound on the
overlap frame is required.

There are two further reductions. A generic rank-one free carrier
with a formal rank-one free completion is already free over the
parameter algebra. Also, full scalar formal lifts agreeing through
the first quantum order are conjugate by a regular formal unit.
Any remaining lower-jet difference can be matched by a finite word
of exact polynomial quantum shears.

Thus finiteness of a good algebraization is invariant among full
formal lifts of the same classical map, after finite exact target
changes and formal bimodule isomorphisms. This does not make the
displayed coefficients of every representative bounded. No generic
finite projective Long action is constructed here.

## Three rings and the overlap that must be used

Let K have characteristic zero, R=K[[t]], a(t)=t^e*c(t), e>=1,
c(0) nonzero, and let A be the finite PBW Weyl algebra over R:
[p_i,x_j]=a(t)*delta_ij. Set
\[
 A_L=A[1/t],\qquad B=\widehat A=\varprojlim A/t^N A,\qquad
 C=B[1/t].
                                                               \tag{1}
\]
A has a finite set of spatial monomials in each element.
B allows that set to grow with t-order. An element of C has a
uniform finite lower bound on t-order, with a polynomial spatial
coefficient at each retained order. It has no rational spatial
coefficient denominators.

Inside C, coefficientwise PBW comparison gives
\[
                A_L\cap B=A,\qquad A_L+B=C,\qquad
                A/t^N A=B/t^N B.
                                                               \tag{2}
\]
For the sum, retain the finitely many negative t-coefficients in
A_L and leave the rest in B. The intersection follows because an
element of A_L already has finite spatial support.

The t-valuation is multiplicative on nonzero elements of C:
the product of the leading coefficients is their product in the
domain S=K[x,p]. In particular
\[
 C^\times=t^{\mathbb Z}B^\times,\qquad
 B^\times=\{U:U_0\in K^\times\}.
                                                               \tag{3}
\]
If g*g^{-1}=1, its two leading polynomial coefficients multiply to
one, so both are nonzero constants. Conversely a constant leading
coefficient permits a geometric-series inverse in B.

The overlap C is smaller than the completed quantum division ring.
For example x_1^{-1} belongs to the latter but not to C.
A division-ring intertwiner alone is not patching data.

## A constructive patching lemma for free pieces

The matrix calculation is written for right modules and column
vectors. Applying it to the opposite ring proves the left-module
version used for graph carriers.

For g in GL_m(C), glue B^m to A_L^m by the relation b=g*d.
The resulting module is finite projective over A. Here is a
constructive proof that also explains why only a finite prefix of
an already admitted overlap matrix is needed.

With block elementary matrices U,V and the constant rotation J,
\[
 U(g)=\begin{pmatrix}I&g\\0&I\end{pmatrix},\quad
 V(h)=\begin{pmatrix}I&0\\h&I\end{pmatrix},\quad
 J=\begin{pmatrix}0&-I\\I&0\end{pmatrix},
\]
\[
       H=\operatorname{diag}(g,g^{-1})
               =U(g)V(-g^{-1})U(g)J.                           \tag{4}
\]
This identity holds over every associative ring.

Choose b>=0 with every entry of g and g^{-1} of valuation at
least -b. Retain their coefficients of t-order less than N to
obtain g_N,h_N in Mat_m(A_L). Take N>=5b+1 and put
\[
 H_-=U(g_N)V(-h_N)U(g_N)J,\qquad H_+=H H_-^{-1}.
                                                               \tag{5}
\]
All elementary factors of H_- and its inverse belong to A_L.
Telescoping the three nonconstant factors in (4) gives
v(H-H_-)>=N-2b, while v(H_-^{-1})>=-3b. Hence
\[
                  H_+-I\in t\operatorname{Mat}_{2m}(B),
           \qquad H_+\in GL_{2m}(B),\quad H=H_+H_-.
                                                               \tag{6}
\]
The inverse in B is the convergent geometric series. These are
valuation inequalities, not checks of a finite collection of b.

The module glued with H is A^{2m}: its coordinates are
b=H_+*v, d=H_-^{-1}*v. Projection onto the first m coordinates
on the gluing data becomes the same idempotent on both patches:
\[
 E=H_+^{-1} e H_+=H_- e H_-^{-1}\in\operatorname{Mat}_{2m}(A),
             \qquad e=\operatorname{diag}(I_m,0).
                                                               \tag{7}
\]
Equality uses that H commutes with e; membership in A uses (2).
Its image is exactly the module glued with g. This proves finite
projectivity and both base-change identifications.

Equation (7) constructs a finite spatial idempotent from admitted
data. It does not infer the existence of g or its exact inverse
from a finite jet.

## Projective generic pieces and their algebra actions

Let P be any finite projective left A_L-module and suppose an
actual C-linear isomorphism
\[
                   \alpha:C\otimes_{A_L}P\xrightarrow{\sim}C^r
                                                               \tag{8}
\]
is given. Glue it to B^r using alpha. To reduce to the preceding
free case, choose Q with P direct_sum Q=A_L^m and put
Q'=Q direct_sum A_L^r. Then P direct_sum Q' is free of rank m+r.
On C, the chosen alpha identifies Q'_C with C^m, using
Q_C direct_sum P_C=C^m. Thus this larger gluing has free pieces.
The compatible projector onto P and B^r descends by (7).
Its image M is finite projective and satisfies
\[
          M[1/t]=P,\qquad B\otimes_A M=B^r,\qquad
          M=\{(p,b):\alpha(p)=b\text{ in }C^r\}.
                                                               \tag{9}
\]
A compatible pair of module maps on the two patches preserves
this fiber product and therefore descends. The construction also
proves uniqueness: for a finite projective A-module, the
cartesian assertion follows from (2) on finite free modules and
then on direct summands.

Apply this to r=1. Give P a commuting full right A_L-action.
Give B its right action through a full scalar formal lift Phi
of a classical symplectic map psi. Require alpha to intertwine
those actions on C, not just their completed division actions.
The descended M then has all right Weyl relations, is R-flat,
and has identified special graph S_psi.

Conversely, every finite-left R-flat full graph carrier is
projective by QIR1, so supplies exactly these patching data.
Consequently
\[
 \begin{split}
 &\text{finite full graph carrier for a chosen formal graph}\\
 &\quad\Longleftrightarrow\
 \text{finite projective generic action with an admitted C comparison}.
 \end{split}                                                   \tag{10}
\]
The finite-corner certificate of QIR1 can be read from this
construction. Compatible right generators, extended by zero on
the complement, give the same finite matrices on both patches.
Their equality in C and (2) put them in the algebraic corner.

One may instead start with a C-free frame of P in which every
right generator image lies in B and reduces to psi. Those images
define their own full formal graph. No preselected formal
representative or explicit infinite intertwiner is needed.

This is the central noncommutative patching situation often called
a Karoubi square; compare
[Johnson, section 3, Proposition 1](https://doi.org/10.3390/math14050861).
The commutative counterpart is
[Beauville--Laszlo gluing](https://stacks.math.columbia.edu/tag/0BNI).
The proof above establishes the needed noncommutative case directly.

## Generic rank one freeness already supplies the integral frame

Suppose M is finite projective over A, its completion is B as a
left module, and M[1/t] is free of rank one. The overlap map
between those two free frames is multiplication by a unit of C.
By (3) it is t^j times a unit of B. Such a B unit changes no
formal lattice, so (2) and (9) give
\[
                          M\simeq t^j A\simeq A.
                                                               \tag{11}
\]
Thus a finite full rank-one graph carrier is integrally free if
and only if its generic left module is free. There is no extra
parameter-compatible frame obstruction after generic freeness.

In particular, for a generic free action, a frame that is
admissible over C makes the standard generic generator images
t-regular as well: conjugation by a C unit preserves valuation.
Those finite PBW elements therefore belong to A_L intersect B=A
and have the required classical symbols. This produces a scalar
bounded full lift.

For a generic projective module that is nonfree, patching still
produces a finite projective carrier and the matrix consequence
of QMC1. It does not remove that generic projective-class
obstruction. QPF1 has nonfree generic fiber and is consistent
with (11).

## Formal units preserve an intrinsic finite quantum jet

Let F_i be a scalar full formal tuple in B. For U in B^times,
divide its nonzero constant leading coefficient so that
U=1+tV. Every scalar algebra commutator is divisible by a(t).
Consequently
\[
                        U F_i U^{-1}-F_i\in t\,a(t)B.
                                                               \tag{12}
\]
The same holds for a C unit by (3). Thus the tuple modulo t*a(t)
is unchanged by an overlap frame change. This includes every
sub-quantum coefficient and the coefficient at the first
quantum order.

For example x -> x, p -> p+a(t)*x is an exact polynomial
quantum shear with the same classical identity map. It is not
conjugate to the identity tuple by a B or C unit, since its
difference is not in t*a(t)B. A formal canonical change and
an inner formal module-frame change are distinct notions.

For any two scalar full formal lifts F,G of the same polynomial
symplectic psi,
\[
 \exists U\in B^\times:\ F_i=U G_iU^{-1}\ \text{for all }i
       \quad\Longleftrightarrow\quad F_i=G_i\pmod{t\,a(t)}.
                                                               \tag{13}
\]
Only sufficiency needs proof. If the first remaining difference
is t^r*d_i with r>=e+1, the Weyl relations at order e+r give
the closed one-cochain condition for the Hamiltonian frame
{psi_i,-}. That frame is a polynomial basis of all derivations;
its cochain complex is the polynomial de Rham complex, as in
QFL1. Thus d_i={H,psi_i} for a polynomial H, up to a common
choice of sign. Conjugation by
1+t^(r-e)*H/c(0) removes that difference. Its inverse exists
in B and its effect starts at order r. Repeating gives a
convergent product U in 1+tB. No spatial bound on U is used.

This works even when psi is noninvertible: inverse Jacobian
matrices are polynomial, and no polynomial inverse map is
needed. Agreement modulo t*a(t) means agreement through t^e,
not merely equality of the classical symbols.

## The remaining finite jet is realized by exact polynomial shears

At the first differing order t^r of two full lifts, their
difference is Hamiltonian by the same H^1 argument. Every
polynomial Hamiltonian is a finite sum of powers of linear
forms. For a product of d linear forms the polarization identity is
\[
 \ell_1\cdots\ell_d
   =\frac1{d!}\sum_{I\subseteq\{1,\ldots,d\}}
       (-1)^{d-|I|}\left(\sum_{i\in I}\ell_i\right)^d .
                                                               \tag{14}
\]
Repeated linear forms give all monomials; constants have zero
Hamiltonian field.

For H=b*ell^d, the generator assignment
\[
               z_j\longmapsto
                  z_j+t^r b d\,\ell^{d-1}\{\ell,z_j\}
                                                               \tag{15}
\]
is an exact polynomial quantum automorphism over R. It fixes
ell, and its inverse replaces b by -b. The quadratic correction
commutators vanish because they are functions of the same ell;
the two linear corrections cancel. Equivalently the map is
exp((t^r/a)*ad(H)), whose displayed generator formula contains
no parameter denominator, including when r<e.

Match r=1,...,e by finitely many such shear words. Their
composition theta and its inverse have finite spatial support.
The adjusted lifts then agree modulo t*a(t), so (13) gives
\[
                          F=\operatorname{Ad}(U)\circ\theta(G),
              \qquad U\in1+tB .
                                                               \tag{16}
\]
There is no uniform degree bound on theta across all possible
finite jets, and its intermediate tuples need not fix a
selected image such as Long R.

Twisting an algebraic carrier by the actual automorphism theta
preserves finite generation, projectivity and freeness. Inner
conjugation by U is an isomorphism of its formal graph modules.
Therefore existence of a finite algebraization, or of a free
rank-one algebraization, is equivalent for any two full formal
lifts of the same psi. This concerns existence of an algebraic
model of the formal module, not boundedness of the displayed
series after an arbitrary formal gauge.

## Consequences and remaining admission tests

For the actual rank-one Long graph the remaining construction
can now be sought as a generic finite projective full action P
that becomes rank-one free over C, with t-regular generator
images and the right classical shadow. Central patching then
supplies an integral finite carrier; a separately guessed
integral lattice is unnecessary. If P is generically free,
(11) already upgrades the result to a scalar bounded lift.
If P is nonfree, the finite matrix consequence is the one in
QMC1.

C-freeness and regularity of the full action have not been
proved for a new Long candidate. Completed division-ring
equivalence alone does not establish them. QGM1 is coherent
and formally free but fails generic left finiteness; it
therefore fails the input to (10), despite having a different
good algebraization. QPF1 supplies projective left modules
but no admitted full right action. QIR1 prevents repairing
a fixed full rational action merely by changing its finite
generic lattice.

The finite-jet reduction (13)--(16) removes a requirement to
control infinitely many representative-dependent gauge
coefficients when comparing formal modules. It does not
supply a universal degree bound, turn fixed-image constraints
into free choices, or identify arithmetic Frobenius centers.
The QFL2 quantum shear control still separates the latter
from the classical t=0 map. QUR1 remains a different singular
inverse-normalization obstruction.

Run check_quantum_graph_patching.py. Exact controls check the
noncommutative elementary identity, a finite patched idempotent
and corner action from an unbounded overlap frame, both
directions of its splitting, the polarization identity and
mixed-position/momentum quantum shears. Patching, valuation,
generic-to-integral rank-one freeness, full-frame homotopy
and all-order convergence use the written proofs. No actual
Long carrier, scalar/matrix Dixmier counterexample, DC2/JC2
solution, independent theorem replay, formal verification,
external review or novelty priority is claimed.
