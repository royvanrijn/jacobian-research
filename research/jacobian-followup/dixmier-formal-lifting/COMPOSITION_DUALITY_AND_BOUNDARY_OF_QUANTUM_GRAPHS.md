# Composition duality and boundaries of finite quantum graph carriers

15 November 2026. **QTC1.** Forward tensor composition preserves finite
full quantum graph carriers, with an explicit block corner certificate
and multiplicative generic and completed indices. Generic index excess
cannot cancel under such compositions.

The opposite operation has a different finiteness condition. Every
finite graph carrier is an integral projective generator, but its dual
is finite on the opposite target side. For any noninvertible classical
symplectic map, composing the carrier with this dual over the source
algebra is not finite on the target left side. For Long the obstruction
is explicit: the off-diagonal classical self-correspondence is the
quadratic closure with its branch polynomial inverted.

Thus a hypothetical finite Long carrier would not automatically supply
the finite closure carrier excluded by QGC1. These results concern all
degrees and every quantum parameter ramification order.

## Full carriers and forward composition

Use R=K[[t]], L=K((t)) in characteristic zero and quantum parameter
a(t)=t^e u(t), e>=1, u a unit. Let A_R,B_R,C_R be separately labeled
copies of the rank-n finite-PBW Weyl algebra with this parameter.
Let M be an A_R-B_R full graph carrier with special graph S_psi and
N a B_R-C_R full graph carrier with special graph S_chi, as in
[QMC1](FINITE_GRAPH_CARRIERS_AND_MATRIX_REALIZATION.md).
All right Weyl relations, parameter flatness and left finiteness are
assumed. By [QIR1](INTEGRAL_GRAPH_PROJECTIVITY_AND_LATTICE_RIGIDITY.md),
both left modules are projective over their parameter Weyl algebras, and
\[
 P=M\otimes_{B_R}N,\qquad
 P/tP=S_\psi\otimes_{S_B}S_\chi
            \simeq S_{\chi\circ\psi}.                       \tag{1}
\]
Here the point-map convention is
(chi composed with psi)_i=chi_i(psi_1,...,psi_(2n)).
Since N is a summand of a finite free left B_R-module, P is a summand
of a finite power of M. It is finite projective on the left and R-flat.
The right C_R-action is defined on this actual module and satisfies all
Weyl relations. No completed tensor product is required.

For explicit certificates write
\[
 M=A_R^s E,\quad
 \rho:B_R\longrightarrow E\operatorname{Mat}_s(A_R)E,
 \qquad N=B_R^q F,\quad
 \sigma:C_R\longrightarrow F\operatorname{Mat}_q(B_R)F.      \tag{2}
\]
The map rho sends 1 to E. Apply rho to every entry of F and of the
right matrices sigma(c), and interpret the resulting blocks as
sq-by-sq matrices:
\[
 H=\rho^{(q)}(F),\qquad
 \tau(c)=\rho^{(q)}(\sigma(c)),\qquad
 P\simeq A_R^{sq}H.                                       \tag{3}
\]
The block matrix H lies in the corner with unit I_q tensor E, is
idempotent, and contains every tau(c). Their relations follow from
the two algebra homomorphisms. If E_0=e11 and F_0=e11, then H_0=e11
in the lexicographic block ordering and tau(z_i)_0 is
chi_i(psi)e11. Formula (3) is a finite normalized certificate.

If D_M>=1 bounds the spatial degrees of E and all rho(z_i), and D_N>=1
bounds those of F and all sigma(z_i), finite PBW substitution gives
\[
              \text{matrix size }sq,\qquad
              \deg H,\deg\tau(z_i)\le D_M D_N.             \tag{4}
\]
Constants in B_R map to scalar multiples of E, whose degree is included
in D_M. Matrix sums do not increase this bound. Scalar coefficient series
may be infinite in t, as in the original finite-PBW definition.

## Generic index excess cannot cancel

Choose generic rank-one fractional frames for M and N. Their right
division actions are embeddings rho_D,sigma_D:D -> D; the tensor
action is rho_D composed with sigma_D. Changing fractional frames
conjugates these embeddings and does not change their indices.
The tower law, separately on the two sides, gives
\[
 N_\ell(P)=N_\ell(M)N_\ell(N),\quad
 N_r(P)=N_r(M)N_r(N),\quad
 \delta(P)=\delta(M)\delta(N).                             \tag{5}
\]
The classical equality follows from the same tower law on function
fields. By [QDI1](QUANTUM_DIVISION_AND_COMPLETION.md), delta is also the
completed quantum degree. Finiteness and N_l,N_r>=delta are supplied by
QMC1/QDI1; no generic left/right index-symmetry theorem is needed.

For either side put epsilon(M)=N(M)-delta(M). Then
\[
 \epsilon(P)=\epsilon(M)\delta(N)+\delta(M)\epsilon(N)
                         +\epsilon(M)\epsilon(N)\ge0.     \tag{6}
\]
In particular N(P)/delta(P) is the product of the two ratios.
Equality N(P)=delta(P) holds exactly when equality holds for each
factor. Composing further full finite carriers cannot remove a strict
generic/completed index gap. This includes factors whose classical
limit is the identity; formal identity contact alone does not make
them invertible. The [QNI1 examples](IDENTITY_LIFTS_AND_UNBOUNDED_INVERSES.md)
provide such factors in rank at least three.

The actual rank-three QNI1 identity lifts have forward degree bound12,
generic indices3 and completed degree1. Their k-fold forward composition,
for every k>=1, therefore has both generic indices3^k, completed degree1
and forward degree at most12^k. It can retain identity contact of any
prescribed parameter order; for contact order at least two it is still
formally inner. Thus genuine polynomial full lifts have arbitrarily
large generic/completed index loss, even with the identity classical
limit. Identity padding gives every rank at least three. This supplies
no rank-two example or degree bound independent of k.

## Every full carrier is an integral projective generator

Set
\[
 M^\vee=\operatorname{Hom}_{A_R}(M,A_R),\qquad
 (b f)(m)=f(mb),\quad (f a)(m)=f(m)a,\qquad
 E_R=\operatorname{End}_{A_R}(M)^{op}.                      \tag{7}
\]
This is a B_R-A_R bimodule, finite projective on the target right
side. Integral projectivity makes dual formation commute with t=0;
there its right S-module is S, with left source action through psi.

Let I be the trace ideal generated by all f(m). The special free
module gives I+tA_R=A_R. The generic M_L is a nonzero projective over
the simple Weyl algebra, hence a generator; thus I[t^-1]=A_R[t^-1].
The quotient A_R/I is finite on the left, killed by one central power
t^b, and equal to t(A_R/I). Therefore
\[
                  I=A_R,\qquad M\text{ is a progenerator}.
                                                                  \tag{8}
\]
Here a progenerator means a finitely generated projective generator.
This conclusion uses bounded central torsion, not a Jacobson-radical
Nakayama argument. The usual Morita evaluation over E_R is an
isomorphism M tensor_(E_R) M^vee -> A_R.

The source B_R acts through E_R. Tensoring over that smaller algebra
instead gives
\[
 C=M\otimes_{B_R}M^\vee,\qquad
 \mathrm{ev}:C\twoheadrightarrow A_R,\quad
 m\otimes f\longmapsto f(m),\qquad K=\ker(\mathrm{ev}).      \tag{9}
\]
The map is balanced and is a bimodule map. Its image is exactly the
same trace ideal, so it is surjective. Finiteness of M^vee on the
target right side is not left finiteness over B_R.

## The classical self-correspondence is not finite

Let f:X=affine^m -> Y=affine^m be any characteristic-zero polynomial
etale map. The following equivalence holds:
\[
 X\times_Y X\longrightarrow X\text{ finite}
                    \quad\Longleftrightarrow\quad
 f\text{ is an automorphism}.                            \tag{10}
\]
The forward implication can be checked after algebraic closure, since
finiteness descends under field extension.

Write U=f(X), an open subset of Y. Its complement has codimension at
least two. Otherwise a nonconstant irreducible polynomial g vanishing
on a missing divisor would pull back to a nowhere-zero polynomial on
affine space, hence a nonzero constant. Dominance rules this out.

The map X -> U is surjective, flat and quasi-compact, hence fpqc.
If the projection in (10) were finite, fpqc descent of finiteness
would make X -> U finite etale; see
[Stacks02LA](https://stacks.math.columbia.edu/tag/02LA).
Every connected finite etale cover of U extends to a finite etale
cover of Y: normalize Y in its function field and use
[purity of branch locus](https://stacks.math.columbia.edu/tag/0BMB)
across the codimension-at-least-two complement. In characteristic zero,
affine space over an algebraically closed field has no nontrivial
connected finite etale covers; compare
[Milne, Section3](https://www.jmilne.org/math/CourseNotes/LECc.pdf#page=27).
Consequently X -> U is an isomorphism. Normality and the codimension
bound give Gamma(U,O)=Gamma(Y,O), so f^* on polynomial coordinate
rings is an isomorphism. This proves (10). The reverse implication
is immediate.

The diagonal of an etale separated map is open and closed. Thus the
affine self-correspondence decomposes into the diagonal and its
off-diagonal complement, and on coordinate rings
\[
 T=S\otimes_{S_B}S=S\times T',\qquad
 \ker(T\longrightarrow S)=T'.                             \tag{11}
\]
For a noninvertible f, T' is not finite over the first S by (10).
Its generic rank is delta-1. It is an etale, quasi-finite algebra,
which does not make it module-finite.

For a full quantum carrier, ordinary tensor base change and projective
dual base change identify C/tC with T. Since the target A_R of ev is
R-flat, reduction of its kernel sequence is still exact. Hence
\[
       C/tC=T,\qquad K/tK=T',\qquad
       f\text{ noninvertible}\Longrightarrow
       C,\ K\text{ are not finite as left }A_R\text{-modules}.
                                                                  \tag{12}
\]
No parameter-flatness assertion about C or K is required for this
implication. No assertion of generic left nonfiniteness follows just
from this special-fiber argument.

## Long has an explicit interior localization in the off-diagonal part

Use [QGC1](RAMIFICATION_AND_FINITE_CLOSURE_CARRIERS.md)'s actual source
chart and further invert r:
\[
 S_0=K[r,r^{-1},v,v^{-1},a,d],\quad
 s=2v+4a-3ra^2,\quad t_0=(ra^3-a^2-va)/2,\quad
 F(X)=rX^3-2X^2+sX+4t_0.                                 \tag{13}
\]
The output coordinates are r,s,t_0,d. Equality of outputs for a second
point with parameter b gives F(b)=0, and its v-coordinate is
v_b=F'(b)/2. Since r is a unit, the second point also has x_b invertible,
so the same source chart applies. Its affine admission requires v_b
to be a unit.

Factor F(X)=(X-a)G(X), where
G(X)=rX^2+(ra-2)X+(ra^2-2a+s) and G(a)=2v. The two factors
are comaximal over S_0. In the full root algebra the diagonal
idempotent is G(b)/(2v), and the off-diagonal component has G(b)=0.

Put A=3ra-2 and w=2rb+ra-2. On this component,
\[
 w^2=Q=A^2-8rv,\qquad
 b-a=(w-A)/(2r),\qquad
 v_b=w(w-A)/(4r),\qquad
 (w-A)(w+A)=-8rv.                                        \tag{14}
\]
The last expression is a unit. Thus b-a is already a unit and
v_b is invertible exactly when w is. The admitted off-diagonal ring is
\[
 T'_0=S_0[w,w^{-1}]/(w^2-Q)
     =\big(S_0[w]/(w^2-Q)\big)[Q^{-1}].                   \tag{15}
\]
It is finite etale of rank two over S_0[Q^-1]. It is not finite over
S_0: it contains 1/Q, which cannot be integral over S_0 because Q is
a nonunit. The divisor Q=0 lies inside this actual source chart.
For example r=a=1 and v=1/8 gives Q=0 while r,v are nonzero.
The same conclusion follows at its generic divisor from partial_v Q=-8r.

The finite quadratic normalization S_0[w]/(w^2-Q) retains the ramified
point w=0. The actual etale self-correspondence omits that point.
This explicitly distinguishes the finite ramified closure obstructed
by QGC1 from the nonfinite correspondence naturally obtained by
composition with a dual. Long's conclusion (12) thus also has this
local proof independent of the global covering argument.

## Arithmetic lifting does not remove this boundary

The primary-source review used two distinct lifting statements.
[Lauritzen--Thomsen, Theorem4.6](https://arxiv.org/html/2604.16636v1)
tests a first-order lift of an existing constant-degree Azumaya
endomorphism with formally smooth center by that center's Poisson bracket.
[Their Weyl specialization](https://arxiv.org/html/2601.23110v1)
also proves first-order Witt lifting over a perfect field of positive
characteristic when the endomorphism degree is less than the prime. Neither statement constructs a uniformly
bounded family across all orders and primes, or prescribes Long's
arithmetic center map for an arbitrary classical-limit lift.

The exact-center hypothesis in
[Kanel-Belov--Fedorov, Theorem4.2](https://arxiv.org/html/2610.33919v1)
therefore remains essential to its stated application. The bounded
quantum shear in [QFL2](BOUNDED_REALIZATION.md) already distinguishes
a classical limit from an arithmetic center map. First-order lifting
does not supply a new comparison theorem.

The constructive gate remains an actual new normalized finite corner.
Forward compositions of admitted corners are valid, with (4)--(6);
the arithmetic-compatible cotangent/shear words already excluded by
[QFL3](ARITHMETIC_COTANGENT_EXCLUSION.md) stay excluded. Dualizing and
forming a natural closure does not automatically stay in this finite
category. A changed quantum action or a different correspondence must
justify its own admission.

## Exact controls and scope

Run check_quantum_graph_composition.py. It verifies a block composition
of nonconstant corner certificates, all six full Weyl relations, its
special graph and the degree bound. It also verifies Long's cubic
factorization, diagonal idempotent, quadratic change of variable,
second-point chart admission and an actual interior branch point.
Generic-index arithmetic controls illustrate strict excess and its
failure to cancel; the tower law is the universal proof.

All projective-generator, composition, index and nonfiniteness statements
are written proofs in every rank, spatial degree and ramification order.
The general geometric equivalence (10) uses the stated external descent,
purity and characteristic-zero covering theorems. The explicit Long
localization does not require the global covering theorem. No actual
Long full carrier, bounded scalar lift, all-representative nonlifting,
DC2/JC2 solution, independent theorem replay, formal verification,
external review or novelty priority is claimed.
