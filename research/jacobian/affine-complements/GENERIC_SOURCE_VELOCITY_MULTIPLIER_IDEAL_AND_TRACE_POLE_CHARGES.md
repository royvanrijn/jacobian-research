# Generic source velocities have an exact target multiplier ideal

PF2RULED141. 16 November 2026. ALL actual complex plane Keller degrees,
in EVERY fixed polynomial source frame. This measures the full finite
normalization/trace-lattice landing left open by shared MS0/MT0. It does
not prove that landing or the full Wright/HC4/JC2 theorem.

## 1. Statement in the actual full finite order

Normalize Jac(P,Q)=1. Put B=C[P,Q] subset R=C[x,y], F=Frac B,
L=Frac R, d=[L:F], S the FULL finite normalization of B in L, and
Lambda its finite trace-dual. Retain ALL omitted primes D_i of
N=Spec S minus the WHOLE specified source X=A2. Write e_i for their
target normal ramification indices and C_i for their target images.

For a nonzero CONSTANT source direction v=(v1,v2), put

    delta_v=v1*partial_x+v2*partial_y,
    ell_v=v1*y-v2*x, a_v=delta_v P, b_v=delta_v Q.

There is a finite sufficient exceptional set of at most #D_i projective
directions such that, outside it,

    q_i=-ord_Di(ell_v)=-min(ord_Di x,ord_Di y)>0,
    min(ord_Di a_v,ord_Di b_v)=-(q_i+e_i) for EVERY i.       (1)

For each DISTINCT image curve C choose its reduced target equation H_C.
Define

    n_C=max_(i:C_i=C) ceil((q_i+e_i)/e_i)
       =1+max_(i:C_i=C) ceil(q_i/e_i),
    J=product_C H_C^n_C in B.                              (2)

Then the following THREE target multiplier ideals agree exactly:

    {f in B: f*a_v,f*b_v in S}
      ={f in B: (f*delta_v)(S) subset S}
      ={f in B: (f*delta_v)(Lambda) subset Lambda}
      =J*B.                                               (3)

There is no hidden finite-point condition in (3): normality accounts
for every finite point after the prime-divisor tests. This ideal is also
the ideal clearing ALL FOUR source Jacobian entries into S, and is
independent of the generic v within the fixed source frame. Polynomial
target changes preserve it as an ideal of B. Nonlinear source changes
require recomputing their own q_i; no frame-independent numerical
exponent is asserted.

Every n_C>=2. If C has ANY omitted unramified prime above it, then
n_C>=3. A single reduced nonproper-divisor factor therefore never
clears the generic first velocities; even its square cannot clear an
unramified omitted prime. These are actual full-order statements,
not merely pole bounds on a chosen Laurent open.

## 2. Generic source-coordinate poles and the Hamiltonian calculation

At every omitted divisor at least one of x,y has a pole. Otherwise its
DVR would define a morphism to the source A2 using x,y. Composing the
open source embedding into separated N agrees with its actual valuative
centre generically, hence everywhere, wrongly putting the omitted centre
inside X. This is the valuative argument already retained in KV0/PF90.

At a fixed D_i, a constant linear combination of x,y can cancel their
lowest-valued terms in at most ONE projective direction. If the leading
residue ratio is not constant, there is no such constant direction.
Exclude these finitely many directions. This proves the first part of
(1), with the displayed fixed q_i.

At a general smooth point of D_i and its image, write a normal parameter
u and a tangent parameter w. The actual target volume has the form

    Omega=dP wedge dQ=dx wedge dy
         =u^(e_i-1)*g(u,w)*du wedge dw,
    ell_v=u^(-q_i)*A(w)+higher u orders, g,A nonzero units. (4)

Residue degree need not be1; choose a general point away from its finite
ramification, or make an unramified residue extension. This does not
change the normal valuations e_i,q_i. The Hamiltonian identity is

    i_(delta_v) Omega=d ell_v.

Its du coefficient gives a tangential vector coefficient of EXACT order
-q_i-e_i; its dw coefficient gives normal order at least1-q_i-e_i.
In compatible local target coordinates t=u^e_i,z=w, the corresponding
target velocities have orders at least -q_i and exactly -q_i-e_i.
The change between (t,z) and (P,Q) has invertible regular Jacobian at
the generic target curve. It preserves the minimum valuation of a
two-entry vector, proving the second part of (1).

Thus the pole charge includes BOTH the source-coordinate pole q_i and
the full ramification index e_i. Keeping only the usual different
exponent e_i-1 would undercount these source velocities.

For e_i=1, retained PF46 Sections6--7 gives an intrinsic regular source
coordinate and a complementary pole q>=2, by exactness of the unit volume.
PF38's arbitrary-source-automorphism lemma makes every negative source
COORDINATE pole a positive multiple of that q. Hence the present q_i>=2
in EVERY source frame. This proves the stronger unramified bound n_C>=3.

## 3. Full-order and trace-lattice stabilization

For f in B, ord_Di(f)=e_i*ord_Ci(f). By (1), the pair f*a_v,f*b_v
is regular at D_i exactly when

    ord_Ci(f)>=ceil((q_i+e_i)/e_i).                         (5)

They are already polynomial on X. Normality of S makes these tests over
ALL D_i sufficient for membership in S. Since B is a UFD, (5) gives
the first ideal in (3), including repeated FINAL image supports.

At the same primes, multiplying the source vector coefficients in (4)
by such f makes its tangential coefficient regular and its normal
coefficient vanish to order at least1. Thus D=f*delta_v is a regular
derivation on N, tangent to EVERY omitted divisor. For every s in S,
D(s) is regular at all height-one primes and therefore belongs to S.
Conversely D(S) subset S puts D(P),D(Q) in S. This proves the second
equality in (3). No local nilpotence is inferred after multiplying by f.

For the last equality, reuse LO1's actual canonical-module identity
omega_N=Lambda*Omega. On the smooth locus the Lie derivative satisfies

    Lie_D(lambda*Omega)=(D lambda+delta_v(f)*lambda)*Omega. (6)

Because D is regular and tangent to the zero divisor of Omega,
its divergence relative to Omega is REGULAR at every height-one point:
in (4), D(u)/u is regular. Normality then gives delta_v(f) in S.
Lie derivative preserves regular canonical forms, and reflexive extension
across normal singular points preserves omega_N. Equation(6) therefore
gives D(Lambda) subset Lambda. Conversely this inclusion and the retained
multiplier identity (Lambda:Lambda)=S imply D(S) subset S by the
commutator [D,m_s]=m_(D s). This treats the divergence correction;
we do NOT incorrectly assume f*delta_v preserves Omega.

Two independent generic directions recover the four source partials by
constant linear combinations. Since each has the same ideal (2), that
ideal clears the entire source Jacobian matrix. A target automorphism
changes the two velocities by a matrix in GL2(B), proving the stated
target-coordinate invariance.

## 4. Exact recovery from finitely many actual field traces

For a rational target function z, let pole_C(z)=max(0,-ord_C z), with
pole_C(0)=0. The same exponent in (2) has the EXACT formula

    n_C=max_(1<=j<=d) {
       ceil(pole_C(Tr_(L/F)(a_v^j))/j),
       ceil(pole_C(Tr_(L/F)(b_v^j))/j)}.                   (7)

For curves outside the nonproper divisor the maximum is0. Indeed Newton
identities in characteristic zero say that a_v is integral over B_(C)
iff its first d power traces lie there. Apply this to f*a_v: its jth
power trace is f^j*Tr(a_v^j). Repeat for b_v. The least possible
ord_C(f) is therefore exactly (7), proving the formula. This is MS0's
finite trace test used PRIME BY PRIME, with the new charge (1).

Consequently the union of pole divisors of these 2d rational traces is
EXACTLY the reduced image union of ALL omitted primes, including wholly
unramified omissions. Every such curve has a trace pole of order >j
for some j<=d. If it has an unramified omitted prime, some order is >2j.
Trace-one cancellation, or a prefix independent of d, is not substituted
for this stated finite list. We do not identify these FIELD traces with
Gaussian/Haar ordinary moment functionals.

These pole divisors detect escaping sheets/nonproperness, NOT absence
of all source points over their target values. In particular they do
not solve PF140's general finite source-image point-gap question.

## 5. A bounded clearing power from the full primitive inverse

Choose ell_v additionally to be a SOURCE-FINITE primitive coordinate in
the sense of PF2OMIT3; this excludes only finitely many further linear
directions. Its primitive inverse equation has leading coefficient c_d
in B and reduced nonproper divisor H=product_C H_C. Retained PF2OMIT3
Section4's complete Newton ledger gives

    M_C=ord_C(c_d)=sum_(i:C_i=C) q_i*f_i,                 (8)

where f_i is the residue degree of D_i. This is POLE charge, not omitted
sheet count; all original primes are included in that identity. Since
M_C>=max_i q_i>=1, (2) implies

    n_C<=M_C+1<=2*M_C.

Consequently

    H*c_d belongs to J*B, c_d^2 belongs to J*B.            (9)

Thus the SQUARE of the actual primitive inverse's leading coefficient
clears ALL FOUR source Jacobian entries into the full S, and multiplying
each constant source direction by it stabilizes S and Lambda. The power2
does not grow with the field degree, boundary count or pole height.
This uses the complete source-finite inverse, not an arbitrary resolvent
or partial normalization. It does not make its raw order normal.

Equation(9) is a clearing result, not unscaled integrality. In particular
the next section shows why it supplies no locally finite source flow
under noninvertibility. No actual noninvertible Keller map is constructed
to claim attainment of a sharp bound on such maps.

## 6. Why target-denominator clearing does not supply a complete source flow

Under noninvertibility J is a nonunit, by (2). A nonzero f in J*B is
therefore a nonconstant target polynomial. In coordinates (s,ell_v)
adapted to delta_v=partial_s, its s-degree is at least2:

* Degree1 would restrict f to a nonconstant affine-linear function on
  some whole source line. PF52's retained one-line criterion forces
  the Keller inverse.
* Degree0 would give f=F_0(ell_v) nonconstant in B. Its Hamiltonian
  F_0'(ell_v)*delta_v is a nonzero source LND preserving B, since
  Hamiltonian brackets of target polynomials stay in B. Retained KX0
  target-field LND descent again forces inversion.

If deg_s f=a>=2, then (f*partial_s)^n(s) has degree
1+n*(a-1), with nonzero leading coefficient in characteristic zero.
These iterates are linearly independent. Thus NONE of these nonzero
target clearings produces a LOCALLY FINITE source derivation under
failure, even though it stabilizes the actual finite S and Lambda.
This is a single-direction rescaling obstruction, not a prohibition
on every possible combination of regular fields on N.

Allowing a rational target multiplier f in F does not enlarge polynomial
source rescalings: since delta_v is constant and nonzero, polynomial
coefficients of f*delta_v force f in R, and retained PF52 gives
F intersection R=B. Thus the same obstruction covers those rescalings.

Shared MT0 reconstructs the RAW polynomial source algebra in generic
directions, and MS0 identifies its finite B-module landing with
integrality of the first two velocities. The ideal here measures that
landing exactly. Clearing two generators by J does not make the entire
source derivative algebra a finite B-module. No inverse, contracting
source flow, complete action or regular full trace pairing follows.

## 7. Source audit and remaining goal

The full MS0, MT0 and LO1 proofs and the relevant PF38/46/PF52 inputs
are read and credited with their original classical-source limits.
PF2OMIT3's source-finiteness and FULL leading-coefficient pole-charge
ledger are reused only under their stated additional coordinate hypotheses.
No new theorem/classification import or predecessor mathematical replay
is used. Small local-power/trace/Newton/ideal controls support the written
ALL-DEGREE proof; they do not classify hypothetical Keller maps.

The original S, Lambda, ALL ramified and unramified primes, and all
final image identifications are retained. Automatic unscaled finite
landing, full trace-dual multiplication/marked-energy integrality,
cyclic/trivial-deck and non-equivariant cases, finite source-image point
gaps and the full Wright/HC4/JC2 theorem remain OPEN. Goal ACTIVE.
No novelty, formal verification, independent whole proof or external
review is asserted.
