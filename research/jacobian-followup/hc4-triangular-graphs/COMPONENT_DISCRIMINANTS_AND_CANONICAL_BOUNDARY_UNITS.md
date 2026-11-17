# Keller components: discriminants and canonical differentials at infinity

2026-11-16. HC4MYLQ0/1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, or full marked-energy integrality; OPEN.
External review/outreach is DEFERRED. The preceding LP0/1 turn made
PROGRESS: one-component geometric acyclicity was shown insufficient.
Here the ACTUAL second Keller component is used throughout.

## 1. LQ0: every fibre discriminant must be strictly above the minimum

Let (P,Q) be an actual complex plane Keller pair. In arbitrary polynomial
source coordinates choose a further generic LINEAR source projection
(X,Y) with deg_Y P=deg P=d. Its leading Y coefficient is a nonzero
constant. Put

    Delta_P(X,t)=disc_Y(P(X,Y)-t).

For d=1 the pair is already invertible: in coordinates (P,z), the
Keller equation makes Q=a*z+b(P), a!=0. Assume d>=2 below.

We explicitly import TWO published inputs from Gwozdziewicz--Ploski,
*On the singularities at infinity of plane algebraic curves*, Rocky
Mountain J. Math.32(2002),139--148:

* Theorem4: if deg_Y f=deg f=d and f admits a polynomial Jacobian mate,
  then f is a coordinate iff deg_X disc_Y(f)=d-1.
* Equation(4), printed pp145--146: for a polynomial with finitely many
  critical points and that same source projection,

      deg_X Delta_f(X,c)=mu(f)+lambda(f)-lambda_c(f)+d-1.  (1)

Here mu is the total finite Milnor number, lambda_c is the nonnegative
infinity Milnor jump at c, and lambda=sum_c lambda_c. The irregular
values are EXACTLY those with lambda_c>0. For the actual Keller P,
mu=0. The precise primary statement, its hypotheses, definitions and
complete relevant proof passages were read, including visual inspection
of printed pp139,140,142,144--146. Its proof credits Assi and Ephraim
for the one-irregular-value/coordinate step and Krasinski and
Cassou-Nogues for the discriminant formulas. These inputs are imported,
not independently reproduced or claimed new here.

Apply Theorem4 to P-c and its ACTUAL mate Q. A coordinate P makes the
whole pair invertible by the elementary degree-one argument in (P,z).
Consequently, in EVERY degree,

    (P,Q) invertible
      iff deg_X Delta_P(X,c)=d-1 for SOME c.          (2)

In the inverse case this equality holds for EVERY c. Under failure,
(1)--(2) give the strict simultaneous inequalities

    deg_X Delta_P(X,c)>=d for EVERY c,
    lambda-lambda_c>=1 for EVERY c,
    # {c:lambda_c>0}>=2.                             (3)

If there were no irregular value or only one, choose c carrying the
entire lambda and contradict (3). The last set is the component's
atypical-value set since it has no finite critical point. The primary's
Theorem1 explains the geometry in the singleton case: its fibre is
reducible with an A1 component. Retained GR0 also excludes that component
for a nonautomorphic Keller pair. The route through (1)--(2) suffices;
no additional classification import is required.

These statements hold for EVERY actual polynomial target-coordinate
component A(P,Q), with its polynomial-coordinate mate B(P,Q), and in
EVERY polynomial source frame followed by the stated monic linear
projection. Thus a failure cannot be hidden by changing a target
coordinate or by concentrating all component degenerations at one
finite value. No minimal discriminant or one-value condition is proved
automatic. This is a sourced classical inverse criterion and its strict
all-fibre consequence, not a new proof of JC2.

## 2. LQ1: generic units are exactly the reducible-fibre factors

Set k=C(P), A=k tensor_(C[P]) C[x,y], and let C=Spec A be the generic
curve. It is geometrically integral by retained HX0. Let Cbar be its
smooth projective completion and S=Cbar-C the FINITE set of punctures.
Write each reduced fibre as

    P-c=product_(j=1)^(r_c) f_(c,j),

with the factors irreducible in the WHOLE source polynomial ring.
Every unit of A is a quotient of source polynomials whose irreducible
factors divide a nonzero polynomial in P. Unique factorization therefore
gives the EXACT free group

    A*/k* = direct_sum_c (Z^(r_c)/Z*(1,...,1)),
    rank(A*/k*)=delta(P)=sum_c(r_c-1).                 (4)

Only finitely many fibres are reducible, as retained in LG1. Relations
in (4) are exactly the products P-c, which are base-field units. There
are no further relations by source unique factorization. In particular

    all fibres irreducible iff A*=k*.                (5)

This is the same integer delta that occurs in LG1's torsion and
saturated-horizontal calculation. We use the elementary unit formula
directly, without identifying raw M, its saturation and its localization.

### A pluricanonical section at infinity forces a reducible fibre

The ACTUAL mate makes dQ a nowhere-zero regular generator of Omega^1_C.
It is exact. On Cbar it is meromorphic and has a POLE at every pole of
the nonconstant function Q. Suppose, over the geometric generic field,
there exist m>=1 and a NONZERO regular pluricanonical section

    omega in H0(Cbar_(kbar),K_Cbar^m)

whose ZERO divisor is supported ENTIRELY on the punctures S_(kbar).
Its ratio

    u=omega/(dQ)^m

is then a unit on the affine generic curve. At a pole p of Q of order e,

    ord_p(dQ)=-e-1,
    ord_p(u)=ord_p(omega)+m*(e+1)>0.                  (6)

Hence u is nonconstant. There is a descent issue: u need not be defined
over k. Choose a FINITE coefficient extension k'/k defining it. The
constant extension A'=A tensor_k k' is a domain because C is geometrically
integral. The norm N_(A'/A)(u) is a unit in A. Every conjugate numerator
omega^sigma is still regular on Cbar, whereas dQ is fixed. Therefore
EVERY factor of that norm has strictly positive order at EVERY pole of
Q after extension to kbar. The norm cannot be a base-field scalar.
Equation(4) now gives

    such a geometric pluricanonical section ==> delta(P)>=1.   (7)

Norms are used only for a FINITE CONSTANT-FIELD extension of the generic
curve, not for the potentially polar full Keller trace pairing. No
integrality of that ordinary pairing is inferred.

This is an all-genus, all-degree boundary condition. If every fibre of
P is irreducible, NO positive canonical power can have a nonzero regular
section with all zeros among those punctures. It is NOT asserted that
such a section always exists. A canonical divisor represented by the
meromorphic dQ can have NEGATIVE coefficients at punctures; its existence
does not provide the EFFECTIVE divisor required in (7).

### Genus one and specified hyperelliptic classes

For generic genus1, a nonzero holomorphic differential has no zero at
all, so (7) applies with m=1. Thus any such actual Keller component has
a reducible fibre. This is the genus-one conclusion of Friedland's
Theorem2.11, credited and reproved here algebraically by (4)--(6).
The complete published argument on printed p231 was reread visually.
No claim that all special fibres are irreducible is introduced.

For geometric generic hyperelliptic genus g>=2, the same argument applies
whenever the punctures contain EITHER a Weierstrass point OR both points
of any fibre of its hyperelliptic degree-two map. Indeed

    K_Cbar ~ (g-1)*D

for a degree-two hyperelliptic fibre D. If D=2p at a Weierstrass point,
the corresponding holomorphic differential has divisor (2g-2)p.
For two distinct points p,q it has divisor (g-1)(p+q).
The local verification is explicit: place that fibre at infinity in
y^2=h(x). In odd degree2g+1, x=t^-2 and dx/y has order2g-2 at the one
point; in even degree2g+2 it has orderg-1 at each of the two points.
There are no affine zeros of dx/y. This supplies the section in (7).

Consequently the ENTIRE all-degree class with all fibres irreducible,
hyperelliptic generic curve and either stated puncture configuration is
excluded for an actual Keller pair. More generally a puncture p with
m*K_Cbar linearly equivalent to m*(2g-2)p for SOME m>=1 is excluded in
an all-irreducible component frame. A general hyperelliptic puncture
configuration or a general higher-genus curve is NOT excluded.
No unproved Kaliman change of component or preservation of these marked
punctures under such a change is assumed.

## 3. Exact controls and the still-open landing

For every d>=2 the familiar non-Keller component

    f_d(x,y)=y+x*y^d

is a whole-plane submersion, with zero fibre A1 union Gm. For t!=0 its
whole fibre has y!=0 and x=(t-y)/y^d, hence is Gm. It is trivial over
Gm_t, and its only atypical value is0. In monic source coordinates
X=x+y,Y=y it is

    f_d=Y+X*Y^d-Y^(d+1), total degree D=d+1.

The imported (1) gives discriminant degree D-1 at0 and a strictly larger
degree at each nonzero value, since the sole infinity jump is positive.
The new checker transcribes exact discriminants in selected new degrees,
where that larger degree is D, together with the whole-fibre
parametrization, the unit-relation lattices and the norm/canonical-order
identities. The general statements are written proofs, not extrapolations
from those samples. The family is NOT claimed new and has NO Keller
mate by (2); it supplies a check of the strict-inequality boundary.
No registered LP/LG/GR/known-map checker is replayed.

Full HC4/JC2 and the actual-source ordinary trace-dual multiplication/full
marked-energy gate remain OPEN. A proof must supply minimal discriminant
at one actual component value, or a qualifying canonical section in an
ACTUAL all-irreducible component frame, or another genuine source
landing. Neither condition is supplied automatically here. No formal
verification, independent whole-proof replay or external review is claimed.
