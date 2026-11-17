# Scalar differential rigidity and one-component descent

Date: 2026-11-16. Claims HC4MYLE0 and HC4MYLE1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: an ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, equivalently regular inverse Gram entries;
the alternative full marked-energy integrality implication is also OPEN.
Neither is proved here. External review/outreach remains DEFERRED.

The preceding LD0/1 identifies why all-conjugate localization or finite
denominators cannot supply the missing whole-source regularity. Here we
strengthen a different, retained obstruction: complete source-flow descent.
For a source field of suitable divergence, descent of JUST ONE component
derivative already forces full target-field descent. No complete field or
such derivative landing is asserted to exist automatically.

## 1. LE0: no new scalar affine differential solutions

Let K=(P,Q):A2->A2 be an actual complex polynomial Keller map,
j=Jac(P,Q) in C*, R=C[x,y], B=C[P,Q], F=C(P,Q), L=C(x,y),
d=[L:F]. Let

    D=(1/j)*Jac(P,-), so D(P)=0,D(Q)=1.                 (1)

Its unique extension to L is the polynomial Nambu derivation. The retained
HX0 theorem gives C(P) relatively algebraically closed in L, because P
has no critical point. The retained KMFI1 theorem gives NO NONTRIVIAL
INTERMEDIATE FIELD GALOIS OVER F. Both hypotheses concern the ACTUAL
Keller extension; neither is assumed for arbitrary finite field extensions.

First its entire constants field is

    L^D=C(P).                                          (2)

Indeed for z in L with Dz=0, differentiate its monic minimal polynomial
over F. The resulting polynomial in z has smaller degree, so all its
coefficients vanish. The original coefficients therefore belong to
F^D=C(P). Thus z is algebraic over C(P); HX0 puts it in C(P).

We claim that for EVERY z in L and a,b in F,

    Dz=a*z+b  ==>  z in F.                             (3)

This is a rational-function statement, not merely a polynomial centralizer.
Trace and norm commute with D: extend D to a normal closure, where its
uniqueness makes it commute with all embeddings, then differentiate the
sum and product of conjugates. Put t=Tr_(L/F)(z)/d. Equation(3)'s hypothesis
gives Dt=a*t+b. Therefore w=z-t has

    Tr(w)=0, Dw=a*w.                                   (4)

If w!=0, its norm n=Norm_(L/F)(w) is nonzero and Dn=d*a*n. Hence

    D(w^d/n)=0, so w^d/n in C(P) by (2).                (5)

It follows that w^d is in F. Since all d-th roots of unity are in C,
F(w)/F is a finite GALOIS Kummer intermediate extension: its minimal
polynomial divides T^d-w^d and splits in F(w). KMFI1 forces F(w)=F.
But an element of F has trace d*w, so (4) then forces w=0, contradiction.
Thus w=0 and z=t is in F. This proves (3). The case d=1 is immediate.

When a=0, only the constants-field statement is needed: Dz in F implies
D(z-Tr(z)/d)=0, and relative closedness proves z in F directly. In
particular, by induction,

    D^n z in F, n>=1  ==>  z in F.                     (6)

For z in R, target-field membership is polynomial membership by the
retained PF52 identity F intersect R=B. Thus a polynomial outside B can
NEVER acquire membership in F under ANY positive iterate of D. It also
cannot solve a first-order affine equation over F for that derivation.
This includes nonzero logarithmic derivatives Dz/z in F as a special case.
The statement is valid after ANY polynomial target coordinate change and
for every nonzero constant Nambu direction, by choosing a complementary
linear target coordinate as P. There is NO assertion about arbitrary
two-variable differential equations or automatic finite differential orbits.

The elementary additive trace and multiplicative norm arguments are not
claimed new differential-algebra methods. Their exact use here combines
the already retained one-parameter regularity and Galois-intermediate
obstructions. Both are essential: in L=C(s,V), F=C(s^2,V), D=partial_V,
the new constant s has Ds=0. In L=C(U,z), F=C(U,V), V=z^e, e>=2,
the constants remain C(U) but

    D=(1/(e*z^(e-1)))*partial_z, Dz=z/(e*V).            (7)

This Kummer example violates KMFI1 and shows why relative closedness
alone does NOT prove the full affine equation assertion. Neither example
is a whole polynomial plane Keller map.

## 2. LE1: one derivative and divergence force full descent

Write an arbitrary polynomial source field as

    delta=A(x,y)*partial_x+C(x,y)*partial_y,
    div(delta)=A_x+C_y, H=delta(P), J=delta(Q).

In target Nambu coordinates let partial_P=(1/j)*Jac(-,Q) and
partial_Q=(1/j)*Jac(P,-). The exact divergence identity is

    div(delta)=partial_P(H)+partial_Q(J).               (8)

It follows either by expanding Jacobians or by taking the Lie derivative
of dP wedge dQ=j*dx wedge dy. No assumption of target-field descent
was used to obtain (8); H and J are functions in the full field L.

If div(delta) and H are in F, (8) puts partial_Q(J) in F. The ADDITIVE
part of LE0 then puts J in F. Conversely full target-field descent makes
both H and J and their target derivatives elements of F, so (8) puts
div(delta) in F. Consequently the exact criterion is

    delta(F) subset F
      <==> delta(P) in F AND div(delta) in F,            (9)

and symmetrically with Q in place of P. Here (9) uses only HX0's
constants-field statement and trace, not the Kummer part of LE0.
PF52 makes the descended H,J genuine target POLYNOMIALS whenever delta
is polynomial. No field denominators survive that conclusion.

Apply the retained KZ1 source-flow theorem. For ANY NONZERO HOLOMORPHICALLY
COMPLETE polynomial source field with div(delta) in F,

    delta(P) in F OR delta(Q) in F  ==> K inverse.      (10)

Thus under actual plane failure, BOTH component derivatives lie outside F
for every such field. The divergence hypothesis is retained explicitly.
In particular it always holds for a LOCALLY FINITE polynomial field:
its exponential gives polynomial automorphisms, with entire time-dependent
coefficients on a finite-dimensional invariant space containing x,y.
Their spatial Jacobians are scalar units. Differentiating those Jacobians
at time zero shows div(delta) is CONSTANT. Such fields are complete.
Therefore a SINGLE component derivative in the target field now suffices
for every nonzero locally finite source derivation, including every LND.
This strengthens KX0's full-field-stability premise; it does not assert a
new source symmetry automatically exists.

Some useful all-degree consequences in ANY polynomial source coordinates:

* If one Keller component, after a constant shift, is weighted homogeneous
  for ANY nonzero integer weights, of either sign or zero, K is inverse.
  Its source Euler field is locally finite and its derivative of that
  ONE component already belongs to B. No second homogeneous component
  or pre-existing invariant target algebra is assumed.
* For every nonzero polynomial h(x), h(x)*partial_y is locally nilpotent.
  Under failure h(x)*P_y and h(x)*Q_y BOTH lie outside F. No division by
  h(x) in F is used; exchange x,y and change polynomial coordinates freely.
* For every nonzero polynomial h, delta=h(x*y)*(x*partial_x-y*partial_y)
  has divergence ZERO and is holomorphically complete, with exact flow
  (x*exp(t*h(x*y)),y*exp(-t*h(x*y))). Thus (10) applies even when h is
  nonconstant and this field is NOT locally finite. Its iterates on x
  are h(x*y)^n*x, linearly independent when h is nonconstant. Replacing
  completeness by polynomiality of every time map is unnecessary.
* Every nonzero affine source field is locally finite, because it preserves
  bounded polynomial degree. Under failure the map from the SIX-dimensional
  affine-field Lie algebra to R/B, delta maps to delta(P) mod B, is INJECTIVE.
  Equivalently P_x,P_y,x*P_x,y*P_x,x*P_y,y*P_y are C-linearly independent
  modulo F, and the same holds for Q. This is NOT F-linear independence.

For a supplied field and component, (3) can also check the antecedent of
(10): a first-order affine equation over F for delta(P) under any Nambu
direction forces delta(P) in F. Such an equation is a certificate; it is
not supplied for every Keller pair.

## 3. The entire source group preserving one component pencil is finite abelian

Suppose K is noninvertible and set

    Gamma_P={sigma in Aut(A2):sigma(P)=a_sigma*P+b_sigma,
             a_sigma in C*,b_sigma in C}.              (11)

Then Gamma_P is FINITE ABELIAN. This includes possible discrete symmetries,
not just connected algebraic actions. Let c_sigma=Jac(sigma)/a_sigma.
Since source polynomial automorphisms have constant nonzero Jacobian,

    Jac(P,sigma(Q))=c_sigma*j.

Consequently sigma(Q)-c_sigma*Q is killed by Jac(P,-). The constants
calculation above puts it in C(P). PF52 and C(P) intersect B=C[P] give

    sigma(Q)=c_sigma*Q+h_sigma(P), h_sigma in C[P].      (12)

Thus sigma induces a TRIANGULAR target polynomial automorphism tau_sigma.
It preserves the ENTIRE nonproperness curve: equivariance by source and
target homeomorphisms preserves the proper locus. The map sigma to tau_sigma
is injective, because a source automorphism over the target is identity
by retained KMFI47. KY1 proves the full ambient target stabilizer of the
entire nonproperness curve is finite. Therefore Gamma_P is finite.

The two diagonal coefficients (a_sigma,c_sigma) define a character map
from this finite triangular group to (C*)^2. Its kernel is trivial:
if both are1, finite order first forces b_sigma=0 by iterating the
translation of P, then forces h_sigma=0 by iterating (12). Injectivity
over the target forces sigma=identity. Hence Gamma_P embeds in an
abelian group and is itself abelian. There is NO claim that every such
finite group is trivial or that any finite symmetry exists automatically.
No unrestricted finite-equivariance or odd-Keller inversion theorem is
imported. The exact same argument works for every polynomial target
coordinate chosen as the first component of a target automorphism.

## 4. Literature and limits of the fibre approach

The polynomial centralizer case Jac(P,r)=0 implies r in C[P] is the
classical Cheng--McKay--Wang theorem, not a new result here. Related
primary context was read in Vered Moskowicz, *The two-dimensional
Centralizer Conjecture*, arXiv:1802.04685v2, Section2.1 Theorem2.3 and
its complete proof, together with the introduction and Section3 limits:

    https://arxiv.org/html/1802.04685v2

The text credits Cheng, McKay and Wang and presents an independent proof
using Formanek's observation. We import NO new theorem from that paper:
our constants-field proof instead uses retained HX0, and the affine
first-order equation uses retained KMFI1 plus written trace/norm algebra.
No priority is claimed for the standard centralizer or differential methods.

The initial fibre investigation also checked Nguyen Van Chau,
*A note on the plane Jacobian conjecture*, arXiv:1005.3866v2, Theorem3
and Section2 proof. Its coordinate criterion requires ALL fibres to be
irreducible and of the SAME GENUS. Those hypotheses were not established
for an arbitrary Keller component. Generic parameter regularity from HX0
does not imply irreducibility of every special fibre, as its retained
submersion x+x^2*y already shows. The current v3 of Chau's
arXiv:0804.3172 only asserts inversion when BOTH components are rational
polynomials; the broader v2 abstract must not be substituted. No genus
constancy, single-component rationality theorem or new classification
is imported from these sources. The investigation changed direction to
the exact constants-field implication instead of assuming a missing fibre
invariance statement.

The previous goal turn is PROGRESS: LD0/1 changed the mathematical authority
and identified all-conjugate descent exactly. This turn supplies additional
all-degree scalar rigidity, one-component flow certificates and full finite
pencil-group constraints. The actual existence of a qualifying complete
field, one-component landing or affine differential certificate remains
UNPROVED. The last trace-dual/marked-energy implication and full HC4/JC2
remain OPEN. Published quartic retirement is unchanged; external review is
deferred. Independent replay, formal verification and external review are
false. Exact controls support the identities and hypothesis boundaries;
they do not formally verify the universal field and geometric proofs.
