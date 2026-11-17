# A quintic phase, finite matrix quantization and algebraic branch obstruction

14 November 2026. **QPH1.** A polynomial generating phase gives a
finite two-by-two full Weyl representation in a rational chart of
Long's map. Its classical limit contains two branches. The positive
branch is exactly Long, including its corrected fourth coordinate.

This construction does not yield a finite-left rank-one Long carrier.
Over every algebraic coefficient cover the finite matrix fiber has
no right-invariant rank-one left submodule or quotient. Keeping both
branches instead prevents global finite-left admission already
classically: the second branch has a pole at X=0.

The nonsplitting also covers arbitrary finite polynomial g_h deformations,
including higher momentum terms vanishing classically. A polynomial
connection proof establishes this without a degree cutoff.

Both branches do separate formally. The positive scalar branch is
another global polynomial-coefficient full Long lift. A Riccati
recurrence proves that its h^n coefficient in S has momentum order n
and total spatial degree 10n-1 for n>=2. Thus this explicit branch
has unbounded momentum as well as position degrees. No conclusion
about every other algebraic model or formal representative is drawn.

## Rational canonical coordinates matching the retained Long map

Write the original Weyl generators as X,Q,P,Z, with
[P,X]=[Z,Q]=h. Put
\[
 \begin{split}
 r&=X(2-3XQ),\quad v=X^{-1},\quad q=v^2/2,\\
 F&=3v-\tfrac53rv^2+\tfrac13r^2v^3,\\
 a&=-X^3P-\tfrac23X(1-3XQ)Z+F,\\
 G&=-\tfrac5{12}v^4+\tfrac2{15}rv^5,\qquad
 d=-Z/(3X^2)+G .
 \end{split}                                                    \tag{1}
\]
These are finite rational operators on X!=0. All coefficients are
on the left. The cotangent vector fields commute, and
G_v=v F_r, so the gradient shifts give exactly
\[
 [a,q]=[d,r]=h,\quad [a,v]=h/v,\quad
               [d,v]=[d,a]=[r,v]=[r,a]=0.                     \tag{2}
\]
No h correction is omitted in (1). First-order normal-operator
commutators equal their Poisson brackets.

Define g=4a-3ra^2. The retained classical Long fiber is
\[
 s_+=g+2v,\qquad
 t_+=(ra^3-a^2-va)/2,\qquad
 D_+=d+a^4/4-3ra^5/10+va^3/2.                                \tag{3}
\]
The canonical fourth-coordinate expression in (3) equals the
actual corrected Long D, not its uncorrected cubic D0.
Substitution in the retained formula verifies this identity;
it is also checked against QPL1's retained classical seed.

## The polynomial phase and its exact Fourier module

Use an auxiliary position s and the quintic polynomial
\[
 \mathcal F(r,s,a)
   =-\frac{as^2}{8}+\frac{a^2s}{2}-\frac{ra^3s}{4}
       -\frac{2a^3}{3}+\frac{3ra^4}{4}
       -\frac{9r^2a^5}{40}.
 \qquad
 \mathcal F_a=-\frac{(s-g)^2}{8}.                             \tag{4}
\]
For the phase mathcal F+a*q the critical equation is
q=(s-g)^2/8. On q=v^2/2 this has branches s=g+2v and s=g-2v.

There is no need to assume an oscillatory integral exists.
On finite polynomials in a,s with coefficients in K[r][[h]],
define left actions
\[
 a_{\rm left}=a,\qquad
 q_{\rm left}=-h\partial_a-\mathcal F_a,
 \quad
 s_{\rm right}=s,\qquad
 t_{\rm right}=h\partial_s+\mathcal F_s .                     \tag{5}
\]
Right operators compose in the opposite order as linear maps.
Thus (5) gives [a,q]=h on the left and [s,t]=h on the right;
the two actions commute by equality of mixed derivatives.

This is the partial Fourier transform of a polynomial exponential
module. The formal stationary-phase viewpoint is standard; compare
[Fang](https://msp.org/pjm/2011/254-1/pjm-v254-n1-p08-p.pdf) and
[Etingof--Kazhdan--Polishchuk, Section2.2](https://arxiv.org/html/math/0003009v2).
Those papers are context, not an algebraization or splitting theorem
being applied here. The formulas and subsequent proofs are direct.

The module is free of rank two on the left over the q,a Weyl algebra,
with basis 1,s. Indeed q raises s degree by two with nonzero leading
coefficient 1/8; its derivative term and every remaining term have
smaller s degree. Leading-degree induction gives generation by 1,s
and proves independence in all degrees.

## Four finite matrices satisfying every full relation

All matrix products below use the q,a Weyl product, [a,q]=h.
Abbreviate
\[
 c_0=8q-g^2,\quad t_0=a^2/2-ra^3/4,\quad
                       c=3a^4/4-9ra^5/20.
\]
The right fiber matrices and its parameter connection are
\[
 S=\begin{pmatrix}0&1\\c_0&2g\end{pmatrix},\quad
 T=\begin{pmatrix}
 t_0&-a/4\\h-ac_0/4&t_0-ag/2
 \end{pmatrix},\quad
 C=\begin{pmatrix}
 c&-a^3/4\\-a^3c_0/4&c-a^3g/2
 \end{pmatrix}.                                               \tag{6}
\]
Expressions such as ac_0 have the indicated noncommutative order.
The connection h*partial_r+mathcal F_r commutes with the left q,a
actions and the right s,t actions. In the basis 1,s its matrix
part is C. Equivalently,
\[
 [S,T]=hI,\qquad [C,S]=hS_r,\qquad[C,T]=hT_r.
\]
Adjoin the independent left canonical pair r,d and set
\[
                         R=rI,\qquad D=dI-C.                  \tag{7}
\]
Then [D,R]=[S,T]=hI and all four mixed commutators vanish exactly.

Equations (6)--(7) give a homomorphism A_2 -> Mat_2(A_2), with
polynomial entries in this chart. It is not an endomorphism of
Mat_2(A_2): no commuting matrix-unit system is supplied. On replacing
the chart by (1), its coefficients are finite Laurent operators in X.
Its left rank is two, rather than QMC1's required rank one.

## Why the algebraic fiber cannot separate its branches

Let L be any finite algebraic function-field extension of
K(r)((h))(v), with the derivation delta=(1/v)*partial_v extended
uniquely. Set O=L[a;h delta]. It is a left and right Euclidean
differential-operator ring. Work with the left module O^2 and its
right action (6).

Give its row basis weights 0,2. Right S has order at most two,
and right T has order at most three. Their leading graded matrices
are
\[
 \sigma_2(S)=
 \begin{pmatrix}0&1\\-9r^2\xi^4&-6r\xi^2\end{pmatrix},
 \quad
 \sigma_3(T)=
 \begin{pmatrix}-r\xi^3/4&-\xi/4\\
                 9r^2\xi^5/4&5r\xi^3/4\end{pmatrix}.
 \quad
 \det\sigma_2(S)=9r^2\xi^4,\quad
 \det\sigma_3(T)=r^2\xi^6/4.                                  \tag{8}
\]
The displayed degrees include the basis shifts. Both determinants
are nonzero over the generic field r!=0.

Suppose N is a nonzero right-invariant rank-one left submodule.
Euclidean ideal theory makes N=O*w for a nonzero row w.
Write wS=P*w and wT=Q*w with finite scalar operators P,Q in O.
Additivity of the weighted order gives orders at most two and three.
The invertible matrices (8) prevent the leading row symbol from
being killed, so these orders are exactly two and three.

The full right relation gives [P,Q]=h. This contradicts
[QSO1](ALGEBRAIC_CURVE_STRING_OBSTRUCTION.md), which excludes such
scalar operators over every algebraic coefficient curve and cover.
Therefore no such N exists. A rank-one quotient would have a
rank-one kernel and is equally excluded.

In fact this fiber bimodule is simple. A proper same-rank submodule
would have a torsion quotient over O, hence finite dimension over L.
Right S,T would act L-linearly on that quotient; taking traces in
[S,T]=hI forces its dimension to vanish. The only other possible
nonzero submodule rank was one, already excluded.

This is an all-order algebraic non-splitting statement, not a finite
search for a projector. This quadratic-symbol proof assumes r!=0;
the polynomial-connection proof below also covers r=0, where g=4a.

## The obstruction also excludes full rank-one extraction

In the full chart let B=K((h))(r)[v,v^-1][a;h delta] and
mathcal A=B[d;h partial_r]. Suppose a rank-one left submodule of
mathcal A^2 is stable under the right R,S,T (and hence possibly D).
Left r and right R=rI show it is stable under differentiation
with respect to d. In the d-left PBW order the finite Taylor
projection sum_k (-d)^k*(partial_d)^k/k! extracts the constant
coefficient and preserves the submodule. Applying it after each
d derivative extracts all coefficients.

Consequently the submodule is induced from its d-free part N0,
which is a rank-one left B-submodule stable under S,T.
Localizing the position coefficients gives the forbidden rank-one
O-submodule above. Thus full algebraic submodules or direct summands
cannot isolate the Long branch either. This assertion concerns this
particular full matrix model, not every algebraization of its formal
branch.

## Keeping both branches prevents global classical admission

In the chart X!=0 the second classical branch has
s_-=g-2v=s_+-4/X. The first branch s_+ is the actual polynomial Long S.
Thus the second branch has an unavoidable pole at X=0, and
\[
                    \operatorname{tr}S_0=2s_+-4/X             \tag{9}
\]
has a pole too.

There is no finite module over the normal ring K[X,Q,P,Z],
with this generic classical two-branch action and preserved by S.
The determinant trick would make both eigenvalues integral over
that ring. They belong to its fraction field, so normality would
make them polynomials, contradicting the pole of s_-.
This rules out a global finite-left flat carrier retaining both
classical branches of this model. It says nothing comparable about
a different carrier with only the Long classical graph.

## Formal splitting gives a different global scalar Long lift

Let a row branch generator be (-g+L,1). Its S-invariance is equivalent
to the Riccati identity
\[
             L^2+[g,L]=4v^2,\qquad
             L=2v+\sum_{n>=1}h^nL_n.                         \tag{10}
\]
At each new h order the unknown coefficient is multiplied by 4v,
a unit in the chart. Thus there is a unique polynomial-in-a
coefficientwise solution over K[r,v,v^-1][[h]].

The negative branch starts at -2v and is constructed identically.
Their two row generators have an invertible classical determinant
4v. They therefore form a formal basis diagonalizing S.
Its two diagonal reductions differ by 4v. The off-diagonal
equations from [S,T]=hI and [D,S]=0 have an invertible leading
coefficient, forcing T,D to be diagonal too. The positive scalar
full tuple is
\[
 \begin{split}
  R_h&=r,\qquad S_h=g+L,\\
  T_h&=(ra^3-a^2)/2-La/4,\\
  D_h&=d+a^4/4-3ra^5/10+La^3/4 .
 \end{split}                                                   \tag{11}
\]
It satisfies every full Weyl relation in all orders, and (3) is its
complete classical tuple.

A finite-order L would have order at most one: for order m>=2,
L^2 has order 2m while [g,L] has order at most m+1, so (10) is
impossible. But the first two recurrence steps give
\[
 L_1=(3ra-2)/v^2,\qquad
       [a^2]L_2=-45r^2/(4v^5)\ne0.                            \tag{12}
\]
Thus even allowing arbitrary completed position coefficients, this
normalized branch vector has unbounded momentum. Algebraic
non-splitting above remains stronger for arbitrary algebraic frames.

## All-order global admission and sharp growth

Write eta=ra and chi=rv. The recurrence in (10) gives
\[
 L_n=v^{1-3n}P_n(\eta,\chi),\qquad
 P_n\in\mathbb Q[\eta,\chi],\quad
 \deg P_n\le n,\quad\operatorname{ord}_a L_n\le n.              \tag{13}
\]
The polynomials are read with chi coefficients on the left of eta.
Here is an induction covering every order. In a normal product,
k contractions differentiate a k times and v by delta^k.
They contribute r^k*v^(-2k)=chi^k*v^(-3k).
For an h^n product the common factor is v^(2-3n).
Differentiating P in eta and multiplying chi^k does not increase
its degree. The first-derivative commutator with g raises degree
at most one from P_(n-1); its second derivative contributes one
chi from P_(n-2). Dividing the recurrence by 4v proves (13).

Both eta and chi are global polynomial operators in (1), and
v^(1-3n)=X^(3n-1). Hence every positive S correction is global.
Also a=X^-1+A with A a global first-order operator.
The normal coefficients of a^k have X pole order at most k:
the vector-field part of A has X component -X^3, so differentiating
a pole does not increase it.

The product P_n(eta,chi)*a^k has pole order at most n+k.
Its prefactor X^(3n-1) clears those poles for k=1,n>=1 and
for k=3,n>=2. This proves global admission of every positive T
term and all D terms from n=2 onward. The only remaining terms
are the finite base plus h*L1 terms in S,T,D-D0.
The retained Fraction seed checks their full normal-order
globality, including all their internal h corrections.
No later admission calculation is needed.

The coefficient of eta^n in P_n is c_n, where
\[
 c_1=3,\qquad
 c_n=-\frac14\left(
       \sum_{i=1}^{n-1}c_i c_{n-i}
           +6(3n-4)c_{n-1}\right)\quad(n>=2).                 \tag{14}
\]
The highest a-order terms in (10) give this recurrence.
Induction yields sign(c_n)=(-1)^(n-1), so c_n never vanishes.
In particular c2=-45/4 and c3=405/4.

In the original spatial Bernstein degree eta has exact degree seven,
chi degree two, and X degree one. Thus the unique eta^n term in
h^n L_n has degree 10n-1 and momentum order n. Every other monomial
of P_n has smaller degree. Earlier terms have smaller spatial
degree, apart from the finite base of degree eleven, and normal
ordering never increases either spatial degree or momentum order.
It follows that
\[
 \deg [h^n]S_h=10n-1,\qquad
       \operatorname{ord}_{P,Z}[h^n]S_h=n\quad(n>=2).           \tag{15}
\]
This is a written all-order growth proof, not extrapolation from jets.

## The phase obstruction survives every finite polynomial g deformation

The quadratic-order proof is not the limit of this construction.
Let k be any characteristic-zero constant field, h a nonzero element,
and let g(a) in k[a] have arbitrary degree m>=1. Write G'=g and H'=g^2,
and use
\[
       \mathcal F_g=-as^2/8+sG(a)/4-H(a)/8,\qquad
       (\mathcal F_g)_a=-(s-g(a))^2/8.                        \tag{16}
\]
The Fourier module (5) is again left free of rank two over k<q,a>.
It has no proper nonzero right-invariant left submodule after
q=v^2/2 and after any further finite algebraic coefficient cover.
This includes every finite-degree g_h with classical reduction
4a-3ra^2: higher a powers multiplied by arbitrary positive powers
of h do not evade this phase-model obstruction.

Here is an all-degree proof independent of QSO1. Over k[a,s] the
rank-one polynomial exponential connection is simple: undo its
polynomial gradient shift of the derivatives, and a nonzero
polynomial ideal stable under both derivatives contains 1.
Partial Fourier transformation is an algebra automorphism, so
the Fourier module is simple too.

Regard that module instead over k[q,s]. Its operator q raises
a degree by 2m with nonzero leading coefficient g_m^2/8; all
other terms, including -h*partial_a, have smaller a degree.
Thus 1,a,...,a^(2m-1) form a free k[q,s] basis. Multiplication
by a and h*partial_s+mathcal F_s give an integrable connection
with polynomial matrices on the whole affine (q,s) plane.
It has no finite singularities and is generically irreducible.

We use the following elementary characteristic-zero fact.
A generically irreducible nonsingular flat algebraic connection
on affine space remains irreducible over every finite algebraic
extension of its function field. Over C, an algebraic horizontal
subspace over such an extension can be transported using a global
analytic fundamental matrix: affine space is simply connected and
the flat connection has no finite singularities. Its resulting
Grassmannian-valued map is single-valued on the entire affine space.
Its Plucker-coordinate ratios are algebraic meromorphic functions
and have trivial monodromy. They are rational: the normalization
of their algebraic function field would otherwise give a connected
finite covering with more than one branch over the complement of
its branch locus, whereas a single-valued branch is fixed by all
monodromy. Connectedness of that covering forces a single branch.
The horizontal subspace therefore descends to the original rational
function field, contradicting irreducibility. General characteristic
zero follows by descent of the finite coefficients to a finitely
generated field and an embedding into C. This does not assume
regular singularity at infinity. For differential-module and
monodromy background see
[Singer, Sections3 and6](https://arxiv.org/html/0712.4124).
Section6, immediately after Proposition9, gives the same entire-function
argument for connectedness in one variable; the affine-plane subspace
argument used here is written out above.

Apply this fact to the polynomial connection above. The substitutions
q=v^2/2 and every subsequent finite coefficient cover are finite
algebraic extensions of k(q,s). A right-invariant left submodule
would be a subconnection there: right S supplies multiplication
by s, right T its connection derivative, and left a supplies the
v derivative. It is a finite module over the position ring because
the ambient connection is finite there. A proper submodule cannot
become the whole connection generically: its torsion quotient would
still carry a connection and differentiation of a minimal
annihilator in s would force that quotient to vanish. Hence
absolute irreducibility excludes it. The previous d-coefficient
argument also applies to full parameter families of this phase form.

Nonconstant g is essential: the constant-g control splits by the
explicit algebraic projector
I/2+S/(4v) after shifting away its constant value. The checker uses
g=0. In that case q no longer raises a degree, so the nonsingular
connection argument does not apply.

Non-splitting of the selected phase model is not an obstruction to
every quantization of its classical branch. For example for g=a,
over [a,v]=h/v the finite scalar pair
\[
       S_*=a+2v,\qquad
       T_*=-S_*^2/8+v^2/2,\qquad [S_*,T_*]=h                 \tag{17}
\]
has exactly the positive classical branch of this linear-g model.
The phase module itself still does not split algebraically.
Thus changing the algebraization can matter even when a finite
scalar alternative exists. This control prevents interpreting
the all-degree phase theorem as all-representative nonlifting.

## What this closes and what remains

The quintic phase gives an actual finite full matrix representation,
rather than a left module missing its right action. Its rank is two
and its second classical branch prevents global finite-left admission.
After the needed position cover, no algebraic rank-one extraction
exists in any finite coefficient cover. Formal extraction does exist
and yields the globally admitted but unbounded Riccati lift (11).

This closes the direct plan of extracting an algebraic rank-one
Long carrier from this phase form, including via a projective
summand and after any finite polynomial g_h deformation. The
nonconstant polynomial connection has an all-degree obstruction
independent of the special order-two/order-three classification. It does not exclude another algebraic model of the same
formal Long branch, just as QGM1 distinguishes a bad model from
intrinsic nonalgebraizability. The scalar and matrix targets in
[QMC1](FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md) remain open.

[Exact controls](check_phase_carrier.py) verify the phase derivatives,
all six matrix relations, the exact rational canonical chart,
the actual classical fourth image, the finite admission seed,
selected Riccati identities and leading growth coefficients, and
the split constant-g control, a higher-order phase deformation, and
a finite scalar alternative for the linear-g branch. Universal freeness,
absolute irreducibility of the nonsingular polynomial connection, formal splitting,
nonsplitting, global admission and growth use written proofs.
There is no independent theorem replay, formal verification,
external review, novelty priority, finite rank-one Long carrier,
bounded Long lift, matrix endomorphism counterexample or DC2 solution.
