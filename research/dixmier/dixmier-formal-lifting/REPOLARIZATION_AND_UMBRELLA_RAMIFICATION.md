# Repolarization and all ramification obstruction for the elliptic normalization module

15 November 2026. **QUR1.** The prescribed rank-three monic
inverse-isogeny module of QRA1 admits no regular formal physical Weyl
lift for any nonzero quantum parameter a(t) in t*k[[t]].
Arbitrary lower t-corrections and all spatial degrees are allowed.
The persistent-block hypothesis of QFR1 is unnecessary for this one
normalization module.

The proof changes the ambient polarization and obtains a rank-two
connection on the open Whitney umbrella. A regular-sequence argument
controls every commuting lower jet at once. At a moving common zero,
its commutator is too divisible by the quantum parameter to cancel
a nonzero curvature entry.

This does not exclude arbitrary higher-multiplicity special modules,
different extensions, or the ordinary rank-one forward Long graph.

## A rank two connection obstruction at a transverse common zero

Let k have characteristic zero, R=k[[t,q,z]], and let
a=a(t)=t^e*c(t), e>=1, c(0) nonzero. Suppose B_1,B_2 in Mat_2(R)
satisfy the flat connection equation
\[
 a(\partial_q B_2-\partial_z B_1)+[B_1,B_2]=0.
                                                               \tag{1}
\]
Write A_i=B_i-tr(B_i)*I/2. Scalar matrices commute, so the traceless
part of (1) is the same equation for A_1,A_2.

Assume the special matrices have the form
\[
 A_1^0=f_0 C_0,\qquad A_2^0=g_0 C_0,\qquad
 C_0=\begin{pmatrix}r_0&s_0\\1&-r_0\end{pmatrix},
 \quad f_0(0)=g_0(0)=0,
                                                               \tag{2}
\]
where (f_0,g_0) are formal coordinates at the origin. If
\[
                 (\partial_q g_0-\partial_z f_0)(0)\ne0,
                                                               \tag{3}
\]
there are no such B_1,B_2 for any e. A nonzero off-diagonal entry of
C_0(0) can first be made the lower-left unit by a constant frame
change, and its unit coefficient absorbed in f_0,g_0.

Here is the all-order proof. Set f=(A_1)_{21}, g=(A_2)_{21}.
The formal inverse function theorem makes (t,f,g) coordinates in R.
In particular, f,g are a regular sequence in R/(a), even though
that coefficient ring has nilpotents.

Equation (1) makes [A_1,A_2]=0 modulo a. For traceless two-by-two
matrices its lower-left and diagonal entries give
\[
 f(A_2)_{11}=g(A_1)_{11},\qquad
 f(A_2)_{12}=g(A_1)_{12}\pmod a.
                                                               \tag{4}
\]
The first Koszul syzygy of the regular sequence f,g says that
f*b=g*d implies d=f*h and b=g*h. Apply it to both equations.
There is a traceless C over R/(a), with C_{21}=1, for which
\[
                 A_1=fC,\qquad A_2=gC\pmod a.
                                                               \tag{5}
\]
This describes every sub-quantum commuting deformation. It does
not assume that a reducing projector lifts.

The equations f(t,q,z)=g(t,q,z)=0 have a unique formal solution
q=q(t), z=z(t), with q(0)=z(0)=0. Evaluate there. Any lift of C
in (5) gives A_1=aU, A_2=aV at this moving point, hence
\[
 [A_1,A_2](t,q(t),z(t))\in a^2\operatorname{Mat}_2(k[[t]]).
                                                               \tag{6}
\]
But the lower-left entry of the derivative term in (1) is exactly
a*(partial_q g-partial_z f). At this point its coefficient after
division by a is a unit by (3). Thus its t-valuation is e, whereas
the commutator in (6) has valuation at least 2e. This contradicts
(1). The proof uses the quotient by a itself, not a census of e.

The same argument works with extra commuting formal parameters:
set them to zero after taking the q,z derivatives. The coordinate
and regular-sequence assertions are relative to those parameters.

## The elliptic correspondence is an open Whitney umbrella

The special normalization is B=k[x,w], with
\[
 u=4w^2-3x,\qquad v=8w^3-9xw,\qquad y=xw.
                                                               \tag{7}
\]
The image ring is
k[u,v,x,y]=k[x,w^2,xw,w^3]. It omits exactly the linear monomial w
from B: its quotient in B is k*w. The conductor in B is (x,w^2).
This is a statement about the special correspondence ring, not
an identification with the whole correction cokernel of QFR1.

Suppose the physical quantum matrices of QFR1 existed over the
target [v,u]=27a, with [Y,X]=a*I. Let them act on columns of the
free right target module A^3. Their left matrix action commutes
with right multiplication by target u,v. Regard x=X,y=Y as
ambient operators. Right multiplication reverses the target
commutator, so the nonzero ambient relations are
\[
                    [v,u]=-27a,\qquad [y,x]=a.
                                                               \tag{8}
\]
All source-target mixed commutators vanish.

The following linear change is an exact quantum canonical change:
\[
 \begin{split}
 q&=(u+3x)/4,& z&=3x/5-2u/15,\\
 p&=-4v/45+8y/15,& r&=v/9+y.
 \end{split}
 \qquad [p,q]=[r,z]=a,
                                                               \tag{9}
\]
with all other brackets between these four generators zero.
There are no ordering corrections in a linear change.

On the normalization put z=x-8w^2/15. Then (9) becomes
\[
              q=w^2,\quad z=z,\quad
              p=\tfrac43 zw,\quad r=\tfrac89 w^3.
                                                               \tag{10}
\]
Rescaling q to q/2, z to 8z/3, p to 2p and r to 3r/8 gives
the standard parametrization
(q_1,q_2,p_1,p_2)=(w^2/2,s,sw,w^3/3).
This is precisely the first open Whitney umbrella in
[Bogaevski--Ishikawa, introduction, equation (1)](https://msp.org/pjm/2002/203-1/pjm-v203-n1-p06-p.pdf).
The algebra and canonical transformation above prove the
identification directly. No stability or quantization theorem
from that paper is imported.

## Completion changes the finite rank from three to two

The change of polarization mixes the source and target algebras.
It is a change for their combined module, not a claim that it is
a source-only or target-only automorphism.

Here are the completion details needed to use it. Complete the
putative free target module spatially at u=v=0 and in t. Its
underlying formal PBW module is
(k[[u,v]][[t]])^3 and is t-flat and t-complete. Its reduction is
the completion of B at x=w=0: this is the unique point over
u=v=0. The old ideal (u,v)B and the new ideal (q,z)B have the
same radical (x,w). Thus their spatial completions agree.

The commuting q,z operators in (9) are jointly topologically
nilpotent on the completed module. This follows first modulo t
from their action by w^2 and x-8w^2/15 on the finite special
module. It then follows in every t quotient: t is nilpotent
there, and each retained PBW star correction has finite
differential order. Consequently formal series in q,z act,
giving an action of k[[t,q,z]].

In this polarization the identified special module is
\[
 k[[z,w]]=k[[q,z]]\oplus k[[q,z]]w,\qquad w^2=q.
                                                               \tag{11}
\]
Lift 1,w to the completed module. The induced map
k[[t,q,z]]^2 -> M is an isomorphism modulo t. Recursively lifting
coefficients makes it surjective by t-completeness.
If its kernel contained a nonzero vector, divide by its least
t-power and reduce modulo t; t-flatness would contradict (11).
The map is therefore an isomorphism. This is a formal free
frame; no algebraic frame or spatial degree bound is asserted.

The other two canonical operators have to act in this frame as
\[
             p=a\partial_q+B_1,\qquad
             r=a\partial_z+B_2,\quad B_i\in\operatorname{Mat}_2(R).
                                                               \tag{12}
\]
Their commutators with q,z give the Leibniz rule, including on
formal series by continuity. Their mutual commutator is zero,
so their matrices satisfy exactly (1).

This argument also applies to a parameter-flat algebraic carrier
with the same identified special left module: its compatible
finite-parameter frames first give the putative formal matrix
action, as in QMC1 and QFR1. Changing that left model cannot
evade the local completion.

## The moving common zero gives a nonzero curvature entry

In the basis 1,w, let W denote multiplication by w. Formula
(10) makes the special connection matrices
\[
 W=\begin{pmatrix}0&q\\1&0\end{pmatrix},\qquad
 B_1^0=\tfrac43 zW,\qquad B_2^0=\tfrac89 qW.
                                                               \tag{13}
\]
Both are traceless. Here f_0=4z/3, g_0=8q/9 are transverse
coordinates and C_0=W has lower-left entry one. But
\[
 \partial_q B_2^0-\partial_z B_1^0
   =\begin{pmatrix}0&4q/9\\-4/9&0\end{pmatrix},
 \qquad \partial_qg_0-\partial_zf_0=-4/9.
                                                               \tag{14}
\]
The criterion (2)--(6) applies. Therefore the prescribed single
monic physical module has no regular formal lift for every
a=t^e*c(t), e>=1. The argument controls arbitrary lower
corrections, rather than requiring a fixed or persistent cusp
block. In particular there is no bounded lift of this module.

Identity spectator factors can also be included. For each such
factor, the combined source-target diagonal has a linear
polarization whose base coordinates are the averages of the
two position and the two momentum coordinates; its conjugate
momenta are the corresponding differences. Its special module
is rank one over that base. The completed total module is
still rank two over all the commuting base variables, so the
relative version of (2)--(6) applies.

Regular invertible quantum changes of frame or source and
target canonical coordinates transport this contradiction.
Polynomial symplectic automorphisms have formal quantum lifts
by QFL1, sufficient here because the exclusion is already
formal. A change allowing poles at the singular point changes
the coefficient category.

## Exact localized connections show where regularity matters

Let alpha=4/3, beta=8/9 and K=diag(1/6,-1/6).
On the two separate punctured charts there are exact solutions:
\[
 \begin{array}{c|cc}
  &B_1&B_2\\ \hline
 z\ne0&\alpha zW&\beta qW+aK/z\\
 q\ne0&\alpha zW+a\,\operatorname{diag}(0,1/(2q))&\beta qW .
 \end{array}                                                   \tag{15}
\]
Direct substitution in (1) gives zero for an indeterminate a.
Neither formula is regular at q=z=0. These are full localized
flat connections, not merely first-order corrections. No
rational or algebraic gluing of the two chosen frames is claimed.

The nonzero derivative hypothesis is substantive: if C is a
constant traceless matrix, B_1=zC and B_2=qC satisfy (1) for
every a. Their scalar factors are transverse, but their
derivative difference is zero.

The old persistent projector condition really can fail.
Replace v by v-6s in the commuting special matrices X_0,Y_0
and then replace Y_0 by Y_0-2sI. This is an exact commuting
polynomial family. At the old target cusp it has
X=4E_20/3+2s(E_01+E_12), Y=0. A projector reducing to E_11
cannot commute even modulo s^2: [P_1,E_20] has no E_01 or
E_12 term to cancel [E_11,2(E_01+E_12)]. Thus QUR1 does not
deduce automatic persistence; it avoids that false assertion.

## Scope and consequences for the main problem

QUR1 closes the higher-ramification escape for this single
prescribed monic inverse-isogeny normalization module, including
broken exact primitive relations and regular spatial formal
corrections. It does not quantify over different special
modules or arbitrary added matrix summands. The two-by-two
syzygy argument is not a proof for higher-multiplicity
normalization sheaves: their commutants are larger. QFR1's
finite-copy exclusions at quantum orders one and two remain
the established statements for those models.

The QRA1 surface/cotangent actions use a different Poisson
structure and remain valid. The ordinary scalar forward Long
graph is smooth and formally quantizable by QFL1; it is not
this singular inverse module. A finite normalized graph corner,
or a representative-independent obstruction to every bounded
Long lift, is still missing.

The geometric identification connects the computation to
singular Lagrangian modules. Smooth-Lagrangian quantization
results such as
[Baranovsky--Ginzburg--Kaledin--Pecharich](https://arxiv.org/html/1403.3493v2)
do not automatically quantize this normalization sheaf.
The obstruction here concerns the full simultaneous algebra
action, beyond quantizing its free underlying bundle.

Run check_umbrella_ramification.py. Its exact controls check
the linear canonical transformation, normalization and rank
change identities, traceless commutator equations, universal
factorization identity, both localized flat connections and
the explicit failed persistent projector. The regular-sequence,
moving-point and completion arguments are written proofs.
No degree or quantum-order census, full Long carrier,
scalar/matrix Dixmier counterexample, DC2/JC2 solution,
independent theorem replay, formal verification, external
review or novelty priority is claimed.
