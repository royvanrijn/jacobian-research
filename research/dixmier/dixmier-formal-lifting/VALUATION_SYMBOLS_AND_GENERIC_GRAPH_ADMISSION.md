# Valuation symbols and automatic formal frames for finite quantum actions

15 November 2026. **QVA1.** A finite generic rank-one full Weyl
action with Gauss-regular generator images automatically has the
rank-one formal frame needed by QPG1. Its polynomial symplectic
shadow is intrinsic to the rational action. The additional
C-freeness test is unnecessary under these hypotheses.

The argument uses the whole initial valuation ideal of a
fractional carrier. Divided left/right operators make its gcd
normalization stable under the residue Hamiltonians. A full
symplectic frame forces that normalized ideal to be the unit
ideal. Completeness then produces a formal generator.

Together with QPG1 this reduces the correct rank-one Long graph
construction to a finite generic rank-one action whose intrinsic
Gauss residues are Long. No such generic action is supplied here.

## The saturated valuation ideal

Use the rings A,A_L,B,C of QPG1, with a(t)=t^e*c(t), e>=1,
and S=K[x_1,...,p_n]. Let D=Frac(A_L), let Dhat be its
t-valuation completion, and denote the valuation by nu.
The residue field is k=Frac(S). The commutator calculus is
\[
 \nu([F,G])\ge\nu(F)+\nu(G)+e,\qquad
 \sigma([F,G]/a)=\{\sigma F,\sigma G\}
 \quad\text{when }\nu(F)=\nu(G)=0.
                                                               \tag{1}
\]
The bracket formula permits a zero right side. It follows from
PBW and the inverse quotient rule, then by continuity in Dhat.

Let I be a nonzero left ideal of C. Define
\[
 I_0=I\cap B,\qquad J=\sigma(I_0)\subset S,\qquad
                      I_0/tI_0=J.
                                                               \tag{2}
\]
The equality uses saturation: I_0 intersect tB=tI_0, since t
is a central unit of C. J is a nonzero polynomial ideal.
Any nonzero element of I can be shifted to valuation zero.
Let d be its polynomial gcd and put J'=J/d; thus J' has gcd one.

For every nonzero polynomial ideal H,
\[
                       (H:H)=S\quad\text{inside }k.
                                                               \tag{3}
\]
Divide H by its gcd. At each height-one prime some element of
the resulting ideal is a unit, so a rational multiplier preserving
it is regular there. The polynomial ring is normal, and the
intersection of those localizations is S.

If rho in Dhat right-multiplies I into itself, normalize it by
a central t-power to valuation zero. It then preserves I_0.
Its rational leading coefficient f satisfies fJ subset J,
so (3) makes f polynomial. Thus polynomiality of normalized
right multiplier residues follows from the ideal; it is not
an assumption about their rational formulas.

## Every residue Hamiltonian preserves the normalized ideal

Fix such a valuation-zero rho, with residue f in S. Choose a
polynomial PBW lift P of f. Initially nu(P-rho)>=1.
If nu(P-rho)>=j for j<e, the operator
\[
                  T_j=(L_P-R_\rho)/t^j:I_0\longrightarrow I_0
                                                               \tag{4}
\]
is well defined: both terms preserve I, their difference has
valuation at least j, and t is invertible in I.

Modulo t its commutator part vanishes because j<e.
It therefore acts on J by multiplication by
r_j=sigma((P-rho)/t^j). Equation (3) makes r_j polynomial.
Replace P by P-t^j*r_j. After finitely many steps,
\[
        \nu(P-\rho)\ge e,\qquad
        \nabla=(L_P-R_\rho)/a:I_0\longrightarrow I_0.
                                                               \tag{5}
\]
Here polynomial lifts of r_j are understood. This is a proof
for every e, not a sequence of numerical ramification tests.

Its residue on J is a scalar first-order operator:
\[
                  \overline\nabla=X_f+r,\qquad
                  X_f=\{f,-\},\quad r\in k.
                                                               \tag{6}
\]
After dividing J by d this becomes
X_f+r', with r'=r+X_f(d)/d, preserving J'.
At each height-one prime choose j in J' not in that prime.
Since X_f(j) is polynomial and (X_f+r')j is in J', r'
is regular there. Normality gives r' in S. Consequently
\[
                         X_f(J')\subset J'.                    \tag{7}
\]
No flatness of the scalar connections is needed.

This is the valuation analogue of the Bernstein-symbol argument
in QHI1. The two conclusions are different: (7) controls a
formal frame over B or C, whereas QHI1 controls a generic
finite-PBW frame using its spatial filtration.

The normalized J' is unchanged by a rational fractional-frame
change, up to a nonzero scalar. Right multiplication of a
fractional ideal by q in D^times scales its initial ideal by
the rational leading coefficient of q; gcd normalization
removes that factor whenever the two realizations lie in C.

## A full symplectic shadow forces a formal generator

Suppose I is preserved by a full rational Weyl tuple rho_i
whose images have nonnegative valuation and whose pairwise
commutators are a*J_ij. Equation (1) forces every canonical
generator to have valuation zero. The preceding multiplier
argument makes its residue psi_i polynomial, and (1) gives
\[
                         \{\psi_i,\psi_j\}=J_{ij}.               \tag{8}
\]
The Hamiltonians X_(psi_i) form a polynomial basis of all
derivations of S, by the symplectic Jacobian identity.
Equation (7) makes J' stable under every coordinate derivative.
A nonzero derivative-stable ideal in a characteristic-zero
polynomial ring is S: differentiate a nonzero polynomial of
least degree until a nonzero constant appears. Hence
\[
                             J'=S,\qquad J=(d).                 \tag{9}
\]

A principal J already gives an actual formal generator.
Lift d to q in I_0. For m in I_0, write its leading coefficient
as b_0*d and subtract a polynomial lift b_0*q. The remainder
is divisible by t inside I_0. Repeat after division by t.
The polynomial coefficients b_j define b=sum t^j*b_j in B,
and completeness gives m=b*q. Thus
\[
                            I_0=Bq,\qquad I=Cq.                 \tag{10}
\]
Conversely, a principal ideal Cq has principal initial ideal
generated by its normalized leading coefficient. Therefore
a nonprincipal initial valuation ideal is an exact obstruction
to C-freeness.

The new formal basis q makes every right image an element
F_i=q*rho_i*q^{-1} in B: right stability first gives F_i in C,
and valuation invariance gives nonnegative valuation.
Its residue remains psi_i, because the residue field k is
commutative. No spatial bound on q or F_i has been proved.

## The generic graph admission criterion

Let P be a finite projective left A_L-module of division rank
one, with a commuting full right A_L-action. Choose a rational
left frame to represent that action by rho_i in D. Its Gauss
valuations and rational leading coefficients do not depend
on this choice: conjugation in D preserves valuation and
acts trivially on the commutative residue field.

Clear a common right denominator to realize P as a left
ideal in A_L. Its extension P_C embeds as a nonzero left
ideal I of C. One can check injectivity in Dhat: P_C is a
direct summand of C^m, so its map into Dhat tensor_C P_C is
injective; the chosen rational frame identifies that extension
with a one-dimensional Dhat-space. The image of P_C itself
lies in C because the original fractional ideal lies in A_L.

If all rho_i are Gauss-regular, (8)--(10) give the admitted
C-frame automatically. Its formal full tuple has the intrinsic
polynomial shadow psi. QPG1 then patches a finite projective
integral full graph carrier. Conversely, a finite full graph
carrier has exactly this rank-one generic action and shadow
by QMC1 and QIR1. Therefore
\[
 \begin{split}
 &\text{a finite full graph carrier with special map }\psi\\
 &\quad\Longleftrightarrow\
 \text{a finite projective generic rank-one full action}\\
 &\hspace{39mm}\text{with Gauss-regular images and shadow }\psi .
 \end{split}                                                   \tag{11}
\]
The generic projectivity input for a finite-left Weyl bimodule
is the cited lemma already used in QMC1; no independent replay
of that external result is claimed.

If P is generically free, QPG1 gives an integral free carrier
and hence a scalar bounded lift. Otherwise its finite matrix
consequence remains the QMC1 statement. This proof does not
turn formal freeness into generic freeness.

For Long, the missing object is now the finite generic full
action with the required Gauss shadow. There is no separate
formal-frame hypothesis once those conditions hold.
A rational Weyl tuple with that shadow but no finite invariant
generic carrier still fails (11).

## Division intertwiners between regular full graphs have no spatial poles

Suppose two full scalar tuples F,G in B are intertwined by
U in Dhat^times. Their classical residues agree automatically.
Normalize U by a central t-power to valuation zero, and let
f=sigma(U) in k^times. The first possible change is at
quantum order, and
\[
 \sigma((F_i-G_i)/a)=\{f,\psi_i\}/f\in S.                        \tag{12}
\]
All logarithmic derivatives of f along the full Hamiltonian
frame are polynomial. Hence every coordinate logarithmic
derivative of f is polynomial. A nonconstant rational
function has a zero or pole on some irreducible divisor;
in characteristic zero its logarithmic differential has a
nonzero simple pole there. Therefore f is a nonzero constant.

It follows that F=G modulo t*a. QPG1 supplies V in B^times
intertwining F,G. The remaining V^{-1}U centralizes a full
formal tuple. Its normalized leading rational coefficient
is annihilated by every Hamiltonian field, hence is constant.
Subtract that scalar coefficient and repeat. Completeness
gives the centralizer K((t)), so
\[
                            U\in C^\times.                    \tag{13}
\]
This excludes hidden rational spatial poles in such a full
regular comparison. It does not hold for partial tuples.

Consequently a finite generic rank-one full action that has
a completed division comparison with a chosen regular full
formal graph has an actual C comparison as well: first use
(10), then apply (13) to the two regular formal frames.
QPG1 alone assumed C-admission; the additional full-action
and valuation-ideal argument here proves it. A bare division
intertwiner with no finite generic carrier is still insufficient.

## An interior projective defect fails in every degree

For one pair [p,x]=a, consider the generic left ideal
\[
                      I_{\rm int}=A_L x^2+A_L(xp+2a).
                                                               \tag{14}
\]
It is the charge-one projective ideal of QPF1 after rescaling
the momentum. Its relation (-p,x) has the generic splitting
column (-x/a,-p/a), since (px-xp)/a=1.

Its extension to C has initial valuation ideal
\[
                      J_{\rm int}=(x^2,xp),\qquad
                      J'_{\rm int}=(x,p).                     \tag{15}
\]
To see exhaustion, suppose the lowest coefficients of
F*x^2+G*(xp+2a) cancel. They satisfy F_0*x+G_0*p=0, so
(F_0,G_0)=h*(p,-x). Subtract a lift of the exact syzygy
p*x^2-x*(xp+2a)=0. This raises the minimum coefficient
valuation without changing the result. After finitely many
steps it reaches the valuation of any fixed nonzero result.
Its leading coefficient is therefore in (x^2,xp); the reverse
inclusion is immediate. This is not a sampled coefficient
calculation.

The normalized ideal is proper and is not stable under the
full derivative frame. Thus this generic projective module
has no Gauss-regular full rational right Weyl action, in any
spatial degree or parameter ramification. Spectator Weyl
factors preserve this conclusion. Even a multiplier with
valuation-zero residue x is impossible, since
{ x,p }=-1 does not belong to (x,p).

The moving-divisor QPF1 ideal is a different positive
left-module control: q=1-t^m*x has leading coefficient one,
so its initial valuation ideal is S and its completion is
free. This passes the present necessary test without
supplying its missing full right action.

## A complete rational logarithmic admission control

For rho_x=x and rho_p=p+a*r(x), r in K(x), every rational r
gives a full rational Weyl pair with identity Gauss shadow.
A finite generic rank-one carrier exists exactly when
\[
                 r(x)=H'(x)+f'(x)/f(x),\qquad
                 H\in K[x],\quad f\in K(x)^\times.              \tag{16}
\]
For necessity, (10) supplies a formal generator in a fractional
realization, whose leading coefficient g is rational in x,p.
The first quantum coefficient of the regular x image makes
partial_p(log g) polynomial. A rational logarithmic derivative
in p is proper at p infinity; being polynomial makes it zero.
Thus g lies in K(x). The first quantum coefficient of the p
image then makes r-g'/g polynomial. Integrate its polynomial
part in characteristic zero.

Conversely the free fractional module A_L*f is stable:
transporting its right images by f gives x and p+a*H'.
In particular r=lambda/x is admitted precisely for integer
lambda; higher poles such as x^{-2} fail. Over an algebraically
closed field the equivalent condition is that every finite
pole of r is simple with integer residue. This supplies exact
positive and negative controls showing that a full rational
pair and a polynomial symplectic residue do not by themselves
supply finite generic generation.

## Verification and remaining scope

The main argument is a t-valuation theorem. It does not
replace the generic Bernstein frame test of QHI1, identify
generic and completed division indices, require a Galois
residue cover, or assume a fixed reducing matrix block.

The residue-ideal method also suggests an exact certificate
for future candidates: compute the entire initial valuation
ideal and normalize its gcd. Valued Weyl Groebner methods
such as
[Hartanto--Ohara](https://doi.org/10.1016/j.jalgebra.2024.07.023)
provide relevant computational context, with homogenization
required by their algorithm. Their monomial term orders are
not automatically the whole coefficient ideal J used here;
an exhaustion argument or certified initial presentation is
still needed. No degree-box campaign is proposed.

Run check_valuation_graph_admission.py. It checks the interior
projective splitting and exact syzygy, the complete commutative
syzygy and nonprincipal symbol certificate, principal-divisor
connection cancellation, rational logarithmic positive and
negative controls, and an exact unbounded inner comparison.
The all-parameter leading-ideal, gcd, completion, generic
admission and intertwiner arguments are written proofs.
No new finite generic Long action, scalar/matrix Dixmier
counterexample, DC2/JC2 solution, independent theorem replay,
formal verification, external review or novelty priority
is claimed.
