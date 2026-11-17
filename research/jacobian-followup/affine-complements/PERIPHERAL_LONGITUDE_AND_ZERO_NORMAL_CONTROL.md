# Peripheral longitude and a sharp zero-normal-degree control

15 November2026. Claim PF2RULED60. Programme status: archived.
The full research goal remains active.

## 1. A peripheral condition on the whole surviving quartic

Retain the FULL surviving quartic normalization of PF2RULED46, its
ramified line R, unramified line D, distinct image supports C,H,
immersive normalization C^nu=A1 and actual S4 monodromy. Let mu,lambda
be a meridian and preferred longitude of the target knot at infinity K.
For the actual monodromy rho supplied by the finite map, one has

    rho(mu) = a transposition,       rho(lambda)=1.          (1)

This holds for every source ruling, source pole order, target degree
and finite contact multiplicity in that quartic packet. The full
inverse image of K consequently has three components: the degree-one
ramification knot and two separate unramified degree-one knots. The
proof retains D and every other finite sheet.

Take the two generic unramified sheets above C as a quadratic etale
algebra over C(t), where t is the normalization parameter of C. Its
normalization over A1 is finite. At a double value with parameters
a_i,b_i and contact m_i, the selected branch has its own quadratic
fold and the other fold restricts to

    z^2=(t-a_i)^(m_i)*unit,

and similarly at b_i. Everywhere else this residual algebra is
unramified over the normalization. In particular intersections with
H add no ramification: pi is etale on D. Thus the squareclass is

    S(t)=product_(m_i odd) (t-a_i)(t-b_i).                  (2)

This is derived for the smooth two-prime packet directly from its
disjoint fold factors. The cuspidal, one-prime marking used in the
older PF2RULED33 is not substituted for this normalization.

The polynomial S has EVEN degree. Hence the completed residual cover
at infinity is split: its two sheets have trivial monodromy around
the large normalization circle. This also includes the globally split
case S=1. A longitude parallel to C at infinity therefore fixes each
of the two residual sheets. Changing that longitude by a meridian power
does not affect this assertion, since rho(mu) fixes those sheets.

The preferred longitude commutes with mu and bounds a Seifert surface
in the knot complement, so its image is an EVEN permutation. A
permutation centralizing (12) and fixing3 and4 is either1 or(12).
Evenness forces1. This proves (1) without assuming an unverified
zero twist for a normal framing. The cycles of the longitude on the
meridian orbits give the stated three components over K.

The quadratic residual rank is essential. In higher degrees, local
branch permutations of the d-2 residual sheets need not commute; an
even total contact count alone does not make their product trivial.
No general S_d longitude assertion is made here.

## 2. An exact additional linear lift condition

PF2RULED54 identifies S4 with V semidirect S3, V=F2^2. Fix a
nonconstant Fox3 representation phi of the knot group, and let z be
a relation cocycle. The meridian is a transposition exactly when

    (I+phi(mu))*z(mu)=0.

Since phi(lambda) commutes with a transposition and is even in S3,
phi(lambda)=I. Formula(1) adds the linear equation

    z(lambda)=0.                                          (3)

Coboundaries vanish on lambda, so this descends to a subspace of
the parabolic H1 of PF54. A surjective lift satisfying BOTH peripheral
conditions exists exactly when that subspace has a nonzero class:
the zero class gives an S3 complement, while every nonzero class
gives S4 by the same invariant-Klein-subgroup argument as PF54.

Consequently an infinity candidate is excluded if every nonconstant
Fox3 representation has zero constrained space. This applies to any
finite knot presentation with the actual preferred-longitude word,
with no cable-length or coefficient-depth cutoff. A passing infinity
representation must still extend to the affine complement with all
local commuting-fold relations and the marked plane source.

The equation is not a consequence of the meridian equation in general.
An exact group control uses the Coxeter presentation of S4 with
adjacent transpositions s1,s2,s3, meridian s1 and the commuting even
word ell=s1*s3. Its parabolic lift space has dimension1, but adding
z(ell)=0 reduces it to zero. This is a GROUP control; ell is not
claimed to be a preferred longitude of a knot with that presentation.

## 3. An unbounded Wronskian reduction for a two-stage binomial key

Let an immersive polynomial branch normalization have coordinates X(t),Y(t)
with degrees p*r,q*r, where p,q>=2 are coprime and r>=2. Suppose its
next and final infinity key is the actual binomial

    F(t)=X(t)^q-Y(t)^p,    deg F=k<p*q*r,   gcd(k,r)=1.     (4)

These are explicit hypotheses. A key with essential lower-weight
terms is not assumed transformable into this binomial. The pole
orders p*r,q*r,k are primitive, so the map is birational: the field
degree [C(t):C(X,Y)] divides all three and hence is1. The infinity
data are the first torus pair(p,q) and ordinary cable winding r,
parameter k. Neumann--Norbury Section4 and Section5.2 give

    eta=a_2-1=k+r*(p+q-p*q)-1.                            (5)

Define the weighted Wronskian

    W=q*X'*Y-p*X*Y'.

There is a universal identity

    Y*F'-p*Y'*F = X^(q-1)*W.                            (6)

Its left leading coefficient is proportional to k-p*q*r, which is
nonzero. Comparing degrees in (6) gives

    deg W=k+r*(p+q-p*q)-1=eta >=0.                       (7)

W cannot vanish identically, since that would make X^q/Y^p constant
and contradict the lower nonzero defect in (4). This is an all-degree
identity, not a Mason bound inferred from small examples.

If an actual surviving quartic branch has (4), PF2RULED59 also gives
eta<=0. Thus it MUST have eta=0 and W a nonzero constant. Negative
eta is excluded for this entire primitive binomial two-stage class,
with no bound on p,q,r. This does not exclude general two-stage keys
with lower-weight terms or every longer cabling sequence.

At eta0, X,Y,F are squarefree and X,Y are coprime. A common root
would force W=0 there. A repeated root of X or Y also forces W=0,
and a repeated root of F contradicts (6), since X,Y are nonzero
there. Moreover X',Y' have no common zero, so the normalization is
immersive. The pair X^q,Y^p attains

    deg F=(p*q*r+1)-(p*r+q*r),

the Davenport--Zannier bound with the corresponding uniform root
partitions (terminology of Pakovich--Zvonkin Definition1.2).
The rational function X^q/Y^p has critical values only0,1,infinity:
its derivative away from X*Y=0 is a nonzero constant times
X^(q-1)/Y^(p+1), and at infinity it approaches1 with order p*q*r-k.

There is also concrete geometry at eta0. The positive weighted Euler
field (p*x,q*y) is transverse to the affine immersion, since its
determinant with the tangent is W. If two branches are tangent at
distinct parameters a,b, their tangent velocities are EQUAL, not
merely proportional: applying this same determinant to both velocities
forces the proportionality factor to be1. These facts do not supply
a symmetry of the target FIELD induced from a source symmetry, or
an everywhere-smooth member of the canonical conic pencil.

## 4. A nodal eta0 curve passing the infinity conditions

Consider the explicit independently checked polynomials

    Y=t^6+4t^4+10t^2+6,
    X=t^9+6t^7+21t^5+35t^3+(63/2)*t.

They satisfy

    X^2-Y^3=-(27/4)*(4t^4+13t^2+32),
    2X'Y-3XY'=378.                                       (8)

Their primitive pole orders are6,9,4. Thus the image normalization is
finite birational A1, its infinity knot is the (winding3,parameter4)
cable of T(2,3), and a_2=4+3*(3+2-6)=1, eta=0. It is immersive
everywhere. This is a realized polynomial curve in the remaining
normal-degree range, rather than a numerical splice without a curve.

Here is an exact certificate of its complete affine singularity set.
Put A(t,s)=(X(t)-X(s))/(t-s), B(t,s)=(Y(t)-Y(s))/(t-s). Then

    Res_s(A,B)=(729/32)*Delta(t),
    Delta=(4t^2-6t+11)*(4t^2+6t+11)
           *(2t^8+12t^6+42t^4+70t^2+63).                 (9)

Delta is squarefree of degree12. The linear subresultant has leading
coefficient coprime to Delta, so every root has exactly one partner,
which in C[t]/(Delta) is

    s(t)=-t*(256t^8+1536t^6+5376t^4+8960t^2+8361)/297.     (10)

The checker verifies A(t,s(t))=B(t,s(t))=0 modulo Delta, s(s(t))=t,
and gcd(s(t)-t,Delta)=1. The tangent determinant
X'(t)Y'(s(t))-Y'(t)X'(s(t)) is a unit modulo Delta. Thus there are
exactly SIX ordinary nodes, with no cusp, tangency or triple value.
Monicity in s makes (9) a complete affine collision certificate,
with no discarded leading-coefficient locus.

The infinity determinant is9: the trefoil determinant3 is multiplied
by the odd winding3 at an even cable parameter4. Even the required
S4 infinity monodromy with trivial longitude EXISTS. To see this,
use the trefoil group <a,b | a^2=b^3> and set

    a=(12), b=(123), M=a*b^(-1)=(13), L=1.

Here M is its meridian and L its preferred longitude, since the
central fibre a^2=b^3=1 and L=M^(-6). Kill the regular fibre in the
(3,4) cable space. Its base is an annulus with one cone point of
order3. Send an outer boundary generator to beta=(13), the inner
meridian to tau=(34), and the cone generator to (beta*tau)^(-1).
The cone relation holds because beta*tau has order3. The outer
meridian is beta^3=M and the outer longitude is beta^(+/-4)=1,
so this matches the companion representation. The inner preferred
longitude is the fibre times tau^(-12), hence1. The image generated
by a,b,tau is S4. Ito's section/fibre bases on complete PDF5 fix
these slopes; orientation choices do not affect these order2/order4
evaluations. No cosmetic-surgery theorem is used.

This curve therefore passes the eta0, double-value, determinant-nine
and peripheral S4 infinity tests simultaneously. This DOES NOT give
an S4 representation of its affine complement, Picard3 or linking
certificates for its sign surface, a finite quartic normalization,
or the marked affine-plane source. None of those is inferred from
an infinity representation. In particular, a proposed universal
contradiction based only on those four passing conditions is false.

## Remaining geometric obligation

Use the whole canonical quartic and its actual normal cubic, shared
binary/cubic index, Cartier Picard/linking data and marked source.
The original goal is still to exclude or realize the full finite map.
The new longitude equation is a necessary global-to-peripheral gate.
The binomial-key reduction isolates an extremal eta0 class, and the
explicit nodal control prevents replacing that gate with an unfounded
strict-positivity or forced-triple assertion.

[The checker](verify_peripheral_zero.py) and
[source manifest](sources/peripheral-zero-control/manifest.json) retain
exact polynomial and permutation certificates and primary conventions.
Written unbounded proofs are distinct from the finite control. No full
quartic or higher-degree exclusion, formal verification, independent
full reproof, external review, novelty priority or goal completion
is claimed.
