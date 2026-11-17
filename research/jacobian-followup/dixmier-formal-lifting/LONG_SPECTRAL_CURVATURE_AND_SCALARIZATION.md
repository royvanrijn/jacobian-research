# Long joint centralizers, spectral curvature and matrix scalarization

15 November 2026. **QSC1.** Long's two commuting observables R,S
have joint polynomial Poisson centralizer K[R,S]. Consequently every
commuting formal quantum deformation of this pair has a commutative
joint centralizer, even if R deforms to a nonpositional operator.

This applies in an actual finite operator algebra, in finite
projective corners with the classical Long graph, and generically
in every finite-left full graph carrier satisfying QMC1. It does
not require a scalar free frame for that carrier.

For a fixed commuting pair, two canonical spectral translation
classes test existence of its individual conjugates. Once those
classes vanish, the remaining mixed relation is exactly one
algebraic curvature class in the commutative joint centralizer.
This is a finite-algebra criterion, not merely a formal criterion.

A concrete consequence is matrix scalarization: if both position
images are scalar matrices, every finite matrix completion traces
to a finite scalar completion, without increasing the degree bound.
Genuinely matrix-valued positions and nonfree projective frames
remain outside that conclusion.

## The classical joint centralizer

Use the actual global Long R,S and the retained rational canonical
chart from [QPH1](QUINTIC_PHASE_AND_BRANCH_OBSTRUCTION.md). Put
r=R, v=X^-1, with conjugate momentum a and parameter momentum d:
\[
 A_0[X^{-1}]=K[r,v,v^{-1},a,d],\quad
 \{d,r\}=1,\quad\{a,v\}=1/v,\quad
 s=S=g(r,a)+2v,\quad g=4a-3ra^2.                            \tag{1}
\]
The other chart brackets vanish. This is an equality of localized
polynomial algebras, not only an equality of their fraction fields.
The inverse expressions for X,Q,P,Z in QPH1 are Laurent in v.

Replacing v by (s-g)/2 gives the ring
K[r,s,a,d,(s-g)^-1]. Commutation with r first removes d.
On the remaining ring the Hamiltonian of s is
\[
                 \{s,-\}=-\frac4{s-g}\partial_a
                     \quad\hbox{with r,s held fixed}.     \tag{2}
\]
Its fraction-field kernel is K(r,s). The intersection
K[r,s,a,(s-g)^-1] intersect K(r,s) is exactly K[r,s]:
a reduced denominator in K[r,s] would have to divide a power
of s-4a+3ra^2. But no nonconstant element of K[r,s] divides
that polynomial, whose a-linear coefficient is the unit -4.
It follows that
\[
                       \operatorname{Cent}_{\{,\}}(R,S)
                                      =K[R,S].             \tag{3}
\]
No bound on a potential invariant, and no assumption on its
momentum order, enters this proof. This is a special feature
proved for this pair, not a theorem for all Lagrangian pairs.

## Every commuting quantum deformation has commutative centralizer

Let A_R be the finite PBW h-Weyl algebra over R0=K[[h]], and
Ahat its h-adic completion. Let F1,F2 in Ahat commute, with
classical reductions R,S. Their complete joint centralizer is
\[
 Bhat=\left\{\sum_{j\ge0}h^j P_j(F_1,F_2):
                                  P_j\in K[u_1,u_2]\right\}.
                                                               \tag{4}
\]
Indeed a centralizer element has classical reduction in (3).
Subtract its corresponding polynomial in the commuting F_i,
divide by h, and repeat. Algebraic independence of R,S makes
this representation unique. Conversely every such series
commutes with both F_i. In particular Bhat is commutative.

This argument permits arbitrary h-small corrections to both
images. It uses no positional normalization, scalar differential
order, finite Fourier comparison or prescribed arithmetic support.
The finite or generic algebraic centralizer can be strictly smaller
than (4), but it is commutative whenever it embeds faithfully there.
It is not automatically K[[h]][F1,F2].

## Finite projective carriers retain this commutativity

For a finite idempotent corner e*Mat_s(A_R)*e with e mod h=e11,
the completed left module is free of rank one. Its completed
corner identifies with Ahat. The finite corner embeds faithfully
in its completion, since its matrix coefficients are h-separated.
If two right images reduce to R,S, their finite joint centralizer
is therefore commutative by (4). No finite conjugating frame is
being inferred.

The same conclusion holds for the generic endomorphism algebra
of a finite-left full graph carrier from
[QMC1](FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md).
Here is the faithfulness point. QMC1 gives a rank-one torsion-free
generic left module P over the Weyl algebra A=A_R tensor K((h)).
An embedding P into A is obtained by a fractional-ideal realization.
After clearing h poles on finitely many generators of its lattice M,
there is an injection M into A_R. Thus M is h-separated.

Every generic endomorphism of P, multiplied by a sufficiently high
power of h, preserves M. It extends to its free completion. If that
extension is zero, separatedness makes it zero on M, hence on P.
Consequently the generic algebra E=End_A(P)^op embeds faithfully
in Ahat[h^-1]. Its centralizer of the two Long right images is
commutative by (4).

A partial finite-left module does not automatically meet the
rank-one, torsion-free or faithfulness hypotheses. Nor does
formal left freeness imply generic left finiteness; QGM1 remains
the countercontrol. The carrier consequence uses the stated QMC1
hypotheses and its external generic-projectivity dependency.

## Universal spectral translation classes in the actual algebra

Let E be a unital characteristic-zero algebra with central
nonzerodivisor h, and fix two commuting elements F1,F2.
The generic version has h invertible; integral applications
below use h-torsion-free algebras. Introduce central
spectral parameters u1,u2 and set
\[
 {\cal M}=E[u_1,u_2]\big/
       \big(E[u](F_1-u_1)+E[u](F_2-u_2)\big).               \tag{5}
\]
As a left E-module this is the regular module E: right evaluation
sends sum a_ij*u1^i*u2^j to sum a_ij*F1^i*F2^j.
The parameter ui acts by right multiplication by Fi.

The two commuting monic expressions Fi-ui have the free Koszul
resolution. Exactness follows successively by monic division in
the central variables, which gives injectivity at each step.
Applying Hom(-,M) gives the finite-algebra complex
\[
 E\xrightarrow{d_0}E^2\xrightarrow{d_1}E,\qquad
 d_0(T)=([F_1,T],[F_2,T]),\quad
 d_1(A,B)=[F_1,B]-[F_2,A].                                 \tag{6}
\]
Thus its cohomology computes the spectral self-Ext groups.
The closed vectors (h,0) and (0,h) define canonical classes alpha1,
alpha2 in H1. They vanish exactly when there are elements T1,T2
of the actual E satisfying
\[
                      [F_j,T_i]=h\,\delta_{ij}.            \tag{7}
\]
These are finite elements when E is the finite operator algebra
or finite corner; completing E changes the question.

Equivalently, right multiplication by Ti defines an E-linear
spectral connection with [nabla_i,ui]=h. This is a relative
connection interpretation of the explicit Koszul calculation.
For a finite projective left Weyl module P, use P[u] in the same
resolution. Its self-Ext complex has coefficients End_A(P)^op,
so the criterion is independent of a free left frame.

The connection viewpoint is standard context; compare the
[principal-parts sequence](https://stacks.math.columbia.edu/tag/09CN)
and the canonical first-order classes in
[Lauritzen--Thomsen, 2026](https://arxiv.org/html/2604.16636).
No identification with that paper's Hochschild class, or with QRF2,
is asserted. The calculation (5)--(7) is direct.

For one operator F, the same calculation gives
Ext1=M/(F-u)M=E/[F,E], and the class of h vanishes exactly
when a finite mate exists. This universal spectral family is
stronger than one selected fiber.

For example take E=K(x)[p;h*partial_x], F=p^2. There is no
finite mate by the constant-coefficient infinity lemma QVL1.
But
\[
                  [p^2,xp/2]=h p^2.                       \tag{8}
\]
The universal class of h is nonzero and is killed by u.
It vanishes after inverting u and in every nonzero spectral fiber.
Even at the zero fiber, the class is exact: for
T0=x^2/(2*h)-x^3*p/(3*h^2), one has
[p^2,T0]=h-(2*x^2/h)*p^2. Thus every individual spectral fiber
has zero translation class, while the universal class is nonzero.
One cannot infer its vanishing from all fibers or apply finite-module
Nakayama to this generally nonfinite spectral cohomology module.
A generic or fiberwise spectral connection is insufficient.

## The full completion obstruction is one curvature class

Assume alpha1=alpha2=0, choose T_i satisfying (7), and suppose
B=Cent_E(F1,F2) is commutative. This is now known for Long pairs
under the hypotheses above. Put
\[
 C=[T_1,T_2]\in B,\qquad
                 \delta_i(b)=[T_i,b]\quad(b\in B).         \tag{9}
\]
Jacobi gives preservation of B and
[delta1,delta2]=ad_C=0 on B.
Every other choice of individual mates is Ti+bi with bi in B.
Because B is commutative,
\[
 C\longmapsto C+\delta_1 b_2-\delta_2 b_1,\qquad
 [C]\in B/(\delta_1 B+\delta_2 B).                          \tag{10}
\]
The derivations on B are independent of the choice of Ti.
Hence the quotient and curvature class are canonical for the pair
once the two translation classes vanish.

All four Weyl relations are realizable in E exactly when
alpha1,alpha2 vanish and the class in (10) is zero. This is an
all-partner criterion: there is no degree, pole or order cutoff.
It does not supply the commuting pair F_i in the first place.

Actual algebra isomorphisms transport the complex, B, derivations
and curvature class. This includes finite Weyl coordinate changes,
with the images transported too, and changes of a projective frame.
The result uses an actual finite algebra, without selecting a
position-field or pseudodifferential completion.

If B=K((h))[F1,F2], curvature can always be removed finitely.
Here delta_i=-h*partial_(Fi), and choose
\[
         b_1=0,\qquad
         b_2=\frac1h\int C(F_1,F_2)\,dF_1.                 \tag{11}
\]
For polynomial degree at most N in the F_i this primitive has
degree at most N+1. If the F_i have spatial degree at most D,
the correction has spatial degree at most D*(N+1).
This is a sufficient bound given the polynomial representation
of C; it does not bound N from a formal representation alone.

For h-regular finite Long images and h-regular individual mates,
(4) makes (11) possible coefficientwise in the completion.
Commutators are divisible by h, so C/h has an integral formal
centralizer expansion. Generic mates with h poles use its Laurent
version; integrality is not asserted for those choices. If the
h-regular classical individual mates are already flat, the
correction can be h-small and preserve them.
The resulting series need not belong to the actual B.
Thus the algebraic curvature class vanishes after completion,
without being proved zero before completion.

## Matrix padding with scalar positions adds no completion freedom

Suppose Fi*I_m are the two position images and arbitrary finite
matrices Ui in Mat_m(E) satisfy
[ Fj*I_m,Ui ]=h*delta_ij*I_m and [U1,U2]=0.
Assume B=Cent_E(F1,F2) is commutative. Then the normalized traces
\[
             T_i=\frac1m\operatorname{tr}_{\rm matrix}(U_i)
                                                               \tag{12}
\]
are already a full scalar completion in E.

To prove the mixed relation, write Ui=Ti*I_m+Bi.
Every entry of Bi belongs to B, and tr(Bi)=0.
Taking normalized trace of [U1,U2] gives
\[
 0=[T_1,T_2]
   +\tfrac1m\operatorname{tr}([T_1I,B_2]+[B_1,T_2I])
   +\tfrac1m\operatorname{tr}[B_1,B_2]=[T_1,T_2].            \tag{13}
\]
The middle terms vanish by tracelessness; the last vanishes
because the matrix coefficients belong to the commutative B.
Ordinary matrix trace does not have this property for arbitrary
noncommutative coefficient algebras.

Consequently fixed Long scalar positions admit a matrix completion
if and only if they admit a scalar one. The normalized trace
does not increase spatial degree and preserves h regularity.
A bounded full formal matrix lift with scalar position images
and Long's classical graph would therefore yield a bounded scalar
lift, and QFL2 would give its nonsurjective Weyl realization.

For a projective carrier the scalarization takes place in E,
not automatically in a scalar Weyl algebra. Stafford's generic
two-copy frame may make both positions genuinely matrix-valued.
No free rank-one frame or scalar DC2 consequence is obtained
from the existence of an arbitrary projective carrier alone.
The QPH1 phase has a genuinely nonscalar S, so (12) does not
apply to it.

Identity Weyl stabilization also preserves completion of the
fixed positions. Add spectator positions y_j. Any mate of an
old Fi must commute with every y_j, hence has no spectator
momentum terms and belongs to E[y]. Evaluation y=0 is an
algebra homomorphism and gives an old flat completion.
Conversely an old completion extends by the identity pairs.
This does not normalize all deformed positions into a fixed pair.

## A boundary model where the algebraic curvature survives

Formal vanishing of (10) cannot replace its algebraic test.
Let qx=1-h*x, qy=1-h*y, and work over K[[h]] with
B_R=K[[h]][x,y,qx^-1,qy^-1].
Form the PBW algebra with [pi,xj]=h*delta_ij and
\[
                  [p_1,p_2]=\frac{h^2}{q_xq_y}.           \tag{14}
\]
It is associative: first adjoin p1 with derivation h*partial_x,
then p2 with derivation h*partial_y on B_R and
p2(p1)=-h^2/(qx*qy). This respects the first Ore relation.
The PBW basis proves h-flatness.

The positions x,y have individual mates Ti=-pi and joint
centralizer B_R, or B=B_R tensor K((h)) generically.
Any further finite mate differs by an element of B.
Their curvature class is nonzero, because
\[
 \operatorname{Res}_{x=1/h}\operatorname{Res}_{y=1/h}
         \frac{h^2\,dx\,dy}{q_xq_y}=1.                     \tag{15}
\]
Every rational exact two-form has zero such iterated residue,
so the class remains nonzero even in K((h))(x,y) and after
finite algebraic covers with normalized differential trace.

Nevertheless the h-adic completion has a full flat pair:
\[
 b_2=-\frac{\log(1-hx)}{1-hy},\quad b_1=0,\qquad
       \partial_x b_2=\frac{h}{q_xq_y}.                    \tag{16}
\]
The logarithm is its formal h series. Its coefficients are
polynomials, but their degrees grow. It kills the curvature
formally and preserves the identity classical reduction.
It is not rational, and the double residue forbids any
algebraic replacement satisfying (10).

This is an explicit boundary model in a twisted Ore algebra,
not a construction inside the finite Weyl algebra and not
a Long counterexample. Its pole divisors escape to infinity
as h tends to zero. It shows exactly why commutativity of a
completed polynomial centralizer cannot supply an algebraic
curvature primitive.

There is also a matrix boundary control when B is noncommutative.
Let E=A_2 tensor Mat_2(K), with position images x,y.
Choose noncommuting constant inner matrices A,B and let
U1=-p1*I_outer+diag(A,0),
U2=-p2*I_outer+diag(0,B).
These are a full flat pair, but their outer normalized traces
have curvature [A,B]/4, nonzero. Thus the centralizer
commutativity hypothesis in (12) is essential.

## Controls and next obligation

[Exact controls](check_spectral_curvature.py) verify the retained
Long chart and Hamiltonian kernels, its inverse localized ring
formulas, a globally finite nonpositional identity deformation,
the two Koszul maps and right evaluation,
the spectral-torsion identity, curvature change and finite polynomial
flattening, nontrivial scalar-position matrix completions and their
trace, the noncommutative-centralizer boundary, the PBW derivation
compatibility and escaping-boundary residue/formal primitive.
Universal centralizer exhaustion, all-deformation induction, faithful
carrier transport, exact spectral/curvature criteria, all matrix sizes
and identity stabilization use the written proofs.

The next actual carrier test is now precise: fix a finite algebra or
projective lattice and a commuting Long pair, compute the two spectral
translation classes there, then compute (10) if they vanish.
Proving a class nonzero excludes every conjugate choice in that
actual algebra. Proving it zero must exhibit finite representatives;
a completed primitive is insufficient.

No nonzero universal class has been computed for every Long
deformation, and no finite Long carrier is constructed here.
No bounded full lift, all-representative nonlifting, scalar or
matrix counterexample, independent theorem replay, formal
verification, external review, novelty priority or DC2 solution
is claimed.
