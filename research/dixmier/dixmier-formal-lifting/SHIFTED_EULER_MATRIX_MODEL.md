# Shifted Euler companion matrices and an exact rational Weyl tuple

The shifted five-by-five companion construction gives actual nonscalar
Weyl operators and a full rational rank-two Weyl tuple. It changes the
quantum identification with QCS1's fixed cubic transport: its target
quintic pair is diagonal with five different quantum entries. The
original pair occurs in the zeroth block. This is a matrix construction,
not a scalar fifth root or a bounded scalar Long lift.

The argument below uses actual quotient division rings and finite
matrices. Its irreducibility statements are scoped to the indicated
relative division ring, with a further exact-R extension proved below.
They are not nonexistence statements about arbitrary nonscalar carriers
or alternative quantum constructions.

## The elementary shifted companion identity

Let D be a characteristic-zero division ring, eta a nonzero central
scalar, and A in D^times, J in D satisfy [J,A]=-5eta A. Indices run
from zero to four. On the right D-vector space D^5 put

\[
 P=\begin{pmatrix}
 0&0&0&0&A\\
 1&0&0&0&0\\
 0&1&0&0&0\\
 0&0&1&0&0\\
 0&0&0&1&0
 \end{pmatrix},\quad
 N=\operatorname{diag}(0,1,2,3,4),\quad
 \widetilde J=JI-\eta N,\quad Q=\widetilde J P^{-1}.       \tag{1}
\]

Ordinary ordered matrix multiplication gives exactly

\[
 P^5=AI,\qquad [\widetilde J,P]=-\eta P,\qquad
 QP=\widetilde J,\qquad PQ=\widetilde J+\eta I,
                         \quad [P,Q]=\eta I.             \tag{2}
\]

The subdiagonal commutators use consecutive shifts by eta. At the
wrap entry, JA-A(J-4eta)=-eta A. The inverse has ones on its
superdiagonal and A^{-1} in its lower-left corner. Thus the nonconstant
entries of Q are J-i*eta on the superdiagonal, i=0,...,3, and
(J-4eta)A^{-1} in the lower-left corner. These are finite actual
division expressions. The same elementary construction works in
size n with [J,A]=-n*eta*A and shifts through (n-1)eta.

The shifts are essential. The unshifted entrywise amplification JI
would instead give

\[
                 [P,(JI)P^{-1}]=\operatorname{diag}(5\eta,0,0,0,0),
\]

not eta I. In QQL1's setting, an eigenoperator of ad(JI) with eigenvalue
-eta would have every entry zero. Formula (1) uses ad(tilde J), whose
matrix-index shifts change that equation. It does not contradict the
scalar spectrum theorem.

The scalar fifth-power linearization used in QQL1 also does not extend
to this matrix algebra. Indeed X=[N,P] is nonzero, but
\[
 \sum_{i=0}^4 P^i X P^{4-i}=[N,P^5]=[N,AI]=0.
\]
Thus the matrix companion root is compatible with the failure of the
scalar injectivity argument in this larger algebra.

## What survives and changes in the QCS1 transport

Use [QCS1, equations (9)--(12)](CLASSICAL_QUINTIC_COVER_SQUARE.md).
Let L=K_0((h)), k=L(r), c=4/r^2, eta=ch, and work first in the actual
relative source division ring

\[
 D=\operatorname{Frac}\bigl(k(v)[b;h\partial_v]\bigr).
\]

Its transported canonical pair U,Z satisfies [U,Z]=h and
J=5cUZ, so [J,U]=-5eta U. Fix a scalar beta and set

\[
 F_\beta(t)=t(t-\beta\eta),\qquad A=UF_\beta(J).           \tag{3}
\]

The scalar J is nonconstant, so these factors are nonzero and invertible
in D. In particular [J,A]=-5eta A, and (1)--(2) apply.
All source identities involving kappa,K,y,U,Z survive under their
entrywise embedding x to xI. The new target identification is different.
Define the target pair by the same formulas as QCS1:

\[
 \mathcal U_\beta=P^5
       [\widetilde J(\widetilde J-\beta\eta I)]^{-1},\qquad
 \mathcal Z_\beta=\frac1{5c}\mathcal U_\beta^{-1}\widetilde J.
                                                                  \tag{4}
\]

Then

\[
 [\mathcal U_\beta,\mathcal Z_\beta]=hI,\qquad
 \widetilde J=5c\mathcal U_\beta\mathcal Z_\beta,\qquad
 P^5=\mathcal U_\beta\widetilde J(\widetilde J-\beta\eta I).
                                                                  \tag{5}
\]

These are exact matrix target-quintic identities. The familiar ordering
beta=1 still gives mathcal U_1=P^3Q^{-2}. But the entries in (4) are

\[
 \begin{split}
 J_i&=J-i\eta,\\
 U_i&=UF_\beta(J)F_\beta(J_i)^{-1},\\
 Z_i&=\frac1{5c}U_i^{-1}J_i,\qquad
 \mathcal U_\beta=\operatorname{diag}(U_i),\quad
 \mathcal Z_\beta=\operatorname{diag}(Z_i).
 \end{split}                                             \tag{6}
\]

Each block is canonical. Block zero is exactly (U,Z); every i>0 has
U_i!=U. Indeed F_beta(J)-F_beta(J-i*eta) is the nonzero polynomial
2i*eta*J-i(i+beta)eta^2. The fixed scalar pair (UI,ZI) is therefore
not the target pair (4).

For QCS1's retained cubic coordinate y, y^3+Uy+U=0, whereas

\[
                     y^3+U_i y+U_i=(U_i-U)(y+1)\ne0
                                      \quad(i>0).        \tag{7}
\]

Thus the old scalar yI does not provide the cubic transport for all
five new target blocks. No replacement cubic branch operators are
supplied here. The exact relative degree-three division-field argument
in QCS1 cannot simply be repeated by counting dimensions in M_5(D),
which is not a division ring.

For h-regular QCS1 parameters, U_i and Z_i reduce to the same classical
u,zeta in every block, because all shifts i*eta vanish modulo h.
Also A_0=p_0^5 for QCS1's prescribed classical p_0. After adjoining
fifth roots of unity the classical companion has eigenvalues
zeta_5^i*p_0; one is the prescribed Long branch. This classical splitting
does not establish an actual scalar branch or a finite polynomial lattice.

## A full rational rank-two completion

Require the QCS1 source ordering parameters and beta to be independent
of r, with coefficients in K_0[[h]] when classical reductions are used.
Let H_s=vb and use QCS1's exact connection

\[
 \widetilde d=d+\frac{v^2}{3r^2},\qquad
 D_s=\widetilde d+\frac{H_s}{5r}
       =d+\frac{v^2a}{5r}+\frac{v^2}{5r^2}.
\]

It commutes with U,Z and has

\[
 [D_s,r]=h,\quad [D_s,J]=-2hJ/r,\quad
 [D_s,\eta]=-2h\eta/r,\quad [D_s,A]=-4hA/r.
\]

Inside M_5(Q_2), where Q_2 is the actual rank-two source Weyl division
ring over L, put

\[
                 \widehat D_s=D_s I-\frac{4h}{5r}N.       \tag{8}
\]

The subdiagonal entries and wrap entry give

\[
 [\widehat D_s,P]=-\frac{4h}{5r}P,\quad
 [\widehat D_s,\widetilde J]=-\frac{2h}{r}\widetilde J,\quad
 [\widehat D_s,Q]=-\frac{6h}{5r}Q.                       \tag{9}
\]

In the wrap entry the first identity is -4hA/r+16hA/(5r)=-4hA/(5r).
The derivative of eta in the middle identity is necessary.
Now define, with the indicated ordering,

\[
 \begin{split}
 R&=rI,\\
 S&=rP+\frac4{3r}I,\\
 T&=\frac r4Q-\frac S{6r}+\frac4{27r^2}I,\\
 F&=\frac{ST}{5r}-\frac{8T}{5r^2}-\frac{S^2}{20r^2},\\
 \mathbb D&=\widehat D_s+F.
 \end{split}                                             \tag{10}
\]

Exactly as the ordered QCS1 curvature calculation,

\[
 \begin{split}
 [\widehat D_s,S]/h&=S/(5r)-8I/(5r^2),\\
 [\widehat D_s,T]/h&=-T/(5r)+S/(10r^2),\\
 [F,S]/h&=-S/(5r)+8I/(5r^2),\\
 [F,T]/h&=T/(5r)-S/(10r^2).
 \end{split}
\]

Consequently the actual rational matrices satisfy the full relations

\[
 [\mathbb D,R]=[S,T]=hI,\qquad
 [R,S]=[R,T]=[\mathbb D,S]=[\mathbb D,T]=0.                \tag{11}
\]

This proves a finite rational matrix realization, with no scalar root
premise. It uses coefficients in L(r) and actual finite division
expressions; the larger field K_0(r)((h)) is not silently substituted.
For regular parameters it is h-regular in the rational Gauss valuation.
It does not assert polynomial admission in the original Weyl algebra.

## Irreducibility over the specified relative division ring

For the same relative D, [QQL1](QUINTIC_ORDERING_LAURENT_OBSTRUCTION.md)
proves that the scalar ad(J)-spectrum is exactly 5eta*Z.
If X commutes with P,Q, it commutes with tilde J=QP, giving
[J,X_ij]=(i-j)eta X_ij. Since |i-j|<5, QQL1 forces X diagonal.
The unit subdiagonal entries of P then force all diagonal entries
equal, while the wrap entry imposes commutation with A. Therefore

\[
 \operatorname{Cent}_{M_5(D)}(P,Q)
             =\{xI:[x,J]=[x,A]=0\}.                      \tag{12}
\]

It has no nontrivial idempotents. In fact the common action is
irreducible over D, not only without a direct-sum projector.
Let a nonzero common invariant right-D subspace contain a vector w
with the minimum possible number of nonzero coordinates. Choose
w_i!=0 and set a=w_i^{-1}(J-i*eta)w_i. The vector
z=tilde J*w-w*a is in that subspace and has coordinate i zero.
For j!=i with w_j!=0, the equality z_j=0 would give

\[
 [J,w_jw_i^{-1}]=(j-i)\eta\,w_jw_i^{-1},
\]

which QQL1 forbids. Thus if w had at least two nonzero coordinates,
z would be a nonzero vector of smaller support, a contradiction.
The subspace contains a coordinate line. The invertible companion P
cycles it through all five coordinates, so it is D^5. There is no
proper common invariant subspace or proper nonzero quotient of any
relative rank.

For (3), even P alone has no invariant right-D line or rank-one
quotient: either condition would produce a fifth root of A in D
by the companion recurrence. For example Pw=wp gives
A w_i=w_i p^5 at any nonzero component, hence
A=(w_i p w_i^{-1})^5, excluded by QQL1. For a rank-one quotient,
a row ell with ell*P=p*ell gives ell_i*A=p^5*ell_i at a nonzero
component, and hence A=(ell_i^{-1}p*ell_i)^5. This is a statement
about this actual companion realization, not all possible matrix carriers.

There is a further exact-R extension. Explicitly assume

\[
 Q_2=\operatorname{Frac}\bigl(D[e;h\partial_r]\bigr),\qquad
 e=\widetilde d,\quad [e,r]=h,\quad [e,v]=[e,b]=0.
\]

Then Cent_{Q_2}(r)=D. Indeed a nonzero element of e-order m has
leading commutator with r of order m-1 and coefficient m*h times
its leading coefficient. If m=0, subtract that coefficient in D;
a nonzero negative-order remainder again cannot commute with r.

A right-Q_2 subspace invariant under rI has a basis matrix B with
identity pivot rows. Invariance gives rB=BC, and reading those pivot
rows gives C=rI. Thus every entry of B commutes with r and lies in D.
If P,Q also preserve the subspace, their restriction matrices read at
the pivot rows also lie in D. Relative irreducibility above excludes
every proper such subspace. Hence no proper common rI,P,Q-invariant
right-Q_2 subspace exists. In particular no invariant scalar component
of the full tuple (10) can be extracted over this Q_2.
This uses the exact first image rI and this ambient division ring;
it is not a no-go theorem for unrelated models, different ambient
division rings or non-invariant extraction procedures.

## Finite operator admission in the two source charts

Let E_kappa=k(kappa)[K], with [kappa,K]=h, be the finite differential
algebra of the source quintic chart. Its coefficient field contains
U and U^{-1}; J has order one in K. More generally take A=UG(J),
0!=G in k(J). In the displayed companion basis the only nonconstant
entries are UG(J), the four J-i*eta, and
(J-4eta)G(J)^{-1}U^{-1}. They all lie in E_kappa exactly when

\[
                      G\in k[J],\qquad G\mid J-4\eta.    \tag{13}
\]

Thus G is a nonzero scalar or a scalar multiple of J-4eta.
Multiplication by U^{+/-1} preserves finite differential order, and
k(J) intersect E_kappa=k[J]: subtract a rational function's polynomial
part in J; a nonzero proper remainder has negative K-order.
This is only a displayed-basis criterion.

A rational gauge can improve it. For beta=1, G_1=J(J-eta), take

\[
 B=\operatorname{diag}(1,J,G_1,G_1,G_1),\qquad
                 P'=BPB^{-1},\quad Q'=BQB^{-1}.           \tag{14}
\]

The five forward edges of P', in cyclic order, are
(J,J-eta,1,1,U). The reverse edges of Q' are
(1,1,J-2eta,J-3eta,(J-4eta)U^{-1}). Thus both are in M_5(E_kappa),
and B^{-1}E_kappa^5 is a common stable free right-E_kappa lattice.
The fifth power after gauge is diag(B_i A B_i^{-1}); it is not being
asserted equal to AI.

This gauge also preserves finite differential-operator admission of
the full tuple in the localized (r,kappa) chart. Since
[D_s,B]B^{-1}=-(2h/r)diag(0,1,2,2,2), its transformed connection is

\[
 B\widehat D_s B^{-1}
   =D_s I+\frac hr\operatorname{diag}
                  (0,6/5,12/5,8/5,4/5).                 \tag{15}
\]

The transformed S',T',F' and mathbb D' are finite differential operators
in the rational position chart (r,kappa), with independent momenta
D_s,K. More precisely S',T' have K-order at most one and
mathbb D' has total differential order at most two in this localized
chart. These bounds describe the nonscalar matrix tuple only. They
still have rational position poles. This is not polynomial
admission in k[kappa,K] or in the original Long chart; the coordinate
change to kappa itself involves momentum denominators.

There is a scoped obstruction in the original finite-operator chart.
Suppose [b,v]=h, H=vb, kappa=c_v*v^5*F_s(H)^{-1}, and
 d=deg_infinity F_s>0, as for the quadratic QCS1 source. Ordinary
b-order gives ord(kappa)=-d, ord(U^{+/-1})=0, and
ord(J)=N=2d+1>1. No diagonal gauge diag(B_i(J)), B_i in k(J)^times,
can make both P and Q finite operators in M_5(k(v)[b]).
The two wrap edges after any such gauge have the form

\[
 p_4=UT(J),\qquad
 q_4=(J-4\eta)T(J)^{-1}U^{-1},\qquad 0\ne T\in k(J).
\]

If t=deg_infinity T, their orders are Nt and N(1-t). Both being finite
forces t=0 or t=1, so one edge is a nonzero order-zero finite operator,
hence a multiplication function f(v). But its exact commutator with J
has eigenvalue -5eta or +5eta. If f'!=0, that commutator has order
N-1>0; if f'=0 it vanishes. Neither case permits such a nonzero scalar
eigenvalue. This rules out every rational-J diagonal gauge in that
original finite-operator chart, and hence in its polynomial PBW subring.
It does not rule out arbitrary rational matrix gauges or other stable
lattices.

## A classical trace obstruction to global admission

The companion has trace zero, so its full second image has classical trace
\[
                    \operatorname{tr}(\overline S)=\frac{20}{3r}.
\]
The [valued matrix lattice theorem](VALUED_MATRIX_LATTICE_AND_SHIFTED_ADMISSION.md)
proves that this trace, and indeed the entire classical characteristic
polynomial, is unchanged by **any** rational matrix gauge whose output is
h-regular. The gauge and its inverse may have parameter poles or singular
classical reductions. Its finite elimination proof uses a central
h-uniformizer and ordinary matrix trace only after commutative reduction.

A globally admitted PBW matrix has polynomial classical entries and thus
a globally regular classical trace. In the original source,
\[
                           r=X(2-3XQ)
\]
has two disjoint reduced smooth zero components, and 20/(3r) has a simple
pole on each. Consequently no rational matrix gauge gives an h-regular
global polynomial second image for this fixed generic representation.
The same theorem excludes every finitely generated embedded full lattice
over the original integral PBW algebra stable even under S; projectivity
or a torsion-free classical reduction is not needed.

This closes the full-rank admission gate for this fixed representation,
including h-polar gauges and coherent lattice changes. It does not exclude
other generic representations, changed matrix images or selecting a
proper formal sector, which changes the generic representation.

## What the construction does and does not supply

Equations (1)--(11) supply an actual nonscalar shifted model and a full
rational Weyl tuple. The beta=1 gauge supplies finite differential
operators in one rational source chart. The changed diagonal transport
(6), original-chart polynomial admission, and scalar realization have
separate tests. Global full-rank admission of this fixed generic model is
excluded by the valued trace theorem, and the irreducibility statements
exclude invariant extraction over the indicated division rings. Other
generic representations and unrelated scalar realizations remain open.

No bounded scalar Long lift, original global PBW admission, DC1/DC2
resolution, formal verification, independent whole-proof replay,
external review or novelty claim is made.

[Exact matrix controls](check_shifted_euler_matrix.py), with
[their retained receipt](shifted_euler_matrix_controls.json), verify the
companion identities, changed canonical transport, shifted connection,
all rational Long commutators, the beta=1 finite-chart gauge and its
transformed connection, and a classical rank-one branch projector.
The irreducibility, universal diagonal-gauge admission and classical
trace obstructions use the written spectrum, differential-order and
regularity arguments above.
