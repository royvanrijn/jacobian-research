# Projective differential generators and inverse monodromy

Date: 2026-11-16. Claims HC4MYLH0 and HC4MYLH1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, equivalently regular inverse Gram entries;
the full marked-energy integrality alternative is also OPEN. Nothing below
supplies either regularity automatically. External review/outreach remains
DEFERRED. The previous goal turn was PROGRESS: LG0/1 gave the whole-source
uniform integration criterion and retained the raw relative quotient.

The earlier LE0 affine first-order criterion suggests trying a quadratic
equation instead. Quadratic equations do NOT give unconditional descent.
They give a different, exact restriction on the normal-closure ACTION.
The distinction between a full source field and a proper intermediate
field is essential when using the published quartic exclusion.

## 1. LH0: full Riccati systems are exactly constant conjugate configurations

Let F=C(t_1,...,t_n), with the full commuting target frame delta_i=partial_i.
Let z be algebraic over F, E=F(z), d=[E:F], and M its normal closure.
Target derivatives extend uniquely to M and commute there. Say that z
has a projective differential system if, for EVERY i,

    delta_i z=a_i*z^2+b_i*z+c_i, a_i,b_i,c_i in F.     (1)

For d>=3 this is equivalent to the following intrinsic condition:
all cross-ratios of four DISTINCT conjugates of z belong to C. When d=3,
the cross-ratio condition is empty and (1) is automatic by the power basis.

To prove the forward direction, write v_i(T)=a_i*T^2+b_i*T+c_i.
For two conjugates u,v, direct subtraction gives

    delta_i(u-v)/(u-v)=a_i*(u+v)+b_i.                 (2)

For four distinct conjugates, apply (2) to

    R(u,v,w,s)=(u-v)*(w-s)/((u-w)*(v-s)).

The logarithmic derivatives cancel, so every delta_i kills this ratio.
The COMMON constants of M are exactly C. Indeed for any common constant
h in M differentiate its monic minimal polynomial over F. Each derivative
is a relation of smaller degree, so all coefficients are common constants
of F, hence in C. Thus h is algebraic over C and belongs to C.
This uses the FULL frame, not algebraic closedness of constants for just
one target direction.

Conversely choose three distinct conjugates and interpolate the unique
quadratic v_i(T) over M agreeing with their delta_i derivatives. For any
fourth conjugate, differentiating its constant cross-ratio forces its
derivative to equal the same v_i, by (2). Thus v_i fits EVERY conjugate.
Every Galois conjugate of v_i fits them too, and uniqueness on three points
makes its coefficients Galois invariant, hence in F. This proves (1).

Normalize three conjugates to 0,1,infinity by a Mobius transformation over
M. The other normalized points are their constant cross-ratios, so all
lie in P1(C). Every sigma in G=Gal(M/F) permutes this finite configuration
by a UNIQUE CONSTANT Mobius map. A permutation preserves cross-ratios
because sigma fixes C; its Mobius map is determined by three points.
Composition agrees on these points, giving a faithful embedding

    G --> PGL2(C).                                   (3)

Faithfulness uses that ALL conjugates generate M. The stabilizer
H=Gal(M/E) is exactly the stabilizer of one point in this projective
action. In particular H is CYCLIC: move its fixed point to infinity;
a finite subgroup of affine transformations has a common fixed point by
averaging, and then embeds in C* by its linear coefficient.

The condition (1) is invariant under changing the full target derivative
basis over F and under ANY Mobius change z'=(a*z+b)/(c*z+d) over F with
ad-bc!=0. Differentiation and substitution give a quadratic in z' again;
the derivative of the rational gauge coefficients must be included.
Also F(z')=F(z). This is a field-generator condition, not a particular
choice of source formula.

There is an exact flatness check. For d>=3 commutation of delta_i,delta_j
gives the polynomial identity

    delta_i(v_j)-delta_j(v_i)+v_i*partial_T(v_j)
                                    -v_j*partial_T(v_i)=0.           (4)

The cubic terms cancel and its degree is at most TWO. Its evaluation at
z is zero; degree d>=3 makes the whole polynomial zero. The three
coefficients are

    delta_i(c_j)-delta_j(c_i)+c_i*b_j-c_j*b_i,
    delta_i(b_j)-delta_j(b_i)+2*(c_i*a_j-c_j*a_i),
    delta_i(a_j)-delta_j(a_i)+b_i*a_j-b_j*a_i.

Thus one obtains a genuine flat projective differential system, not
merely several unrelated quadratic identities in a finite algebra.

## 2. The finite projective classification, with its action data

For clarity we recall the elementary classification used in (3). A finite
subgroup of PGL2(C) is cyclic, dihedral, A4, S4 or A5. The point stabilizers
in a nontrivial orbit are crucial; the abstract group name alone is not
enough. Here is the unbounded classification argument rather than a finite
group search.

The quotient of P1 by a finite G is P1. A point stabilizer is cyclic, as
above. If its branch indices are e_1,...,e_s, Riemann-Hurwitz gives

    sum_j(1-1/e_j)=2-2/|G|.                          (5)

There are at most THREE branch values and not fewer than two for G!=1.
For two, since e_j<=|G|, (5) forces both equal |G| and G fixes two points,
hence is cyclic. For three, order e_1<=e_2<=e_3. The reciprocal sum is
greater than one, forcing exactly

    (2,2,m), (2,3,3), (2,3,4), (2,3,5).

Their orders from (5) are 2m,12,24,60. The inertia generators give the
orientation triangle presentations. For (2,2,m), two involutions whose
product has order m give the usual dihedral group of order 2m. The other
three spherical triangles tessellate the sphere: the orientation groups
have respectively 12,24,60 elements by their spherical area. The standard
permutation realizations are A4,S4,A5. For example generators
(12)(34),(123); (12),(234); and (12)(34),(135), respectively, have product
orders3,4,5 and generate groups of exactly these orders. Hence the
triangle quotients have the same order and are these groups. This proves
the needed classification. The spherical classification and Riccati
orbit degrees are classical; no novelty is claimed for them.

For a point outside branch orbits H=1 and its orbit has size |G|. In a
branch orbit H has exactly one of the orders in the corresponding
signature, and the orbit size is |G|/|H|. These are FULL point stabilizers,
not arbitrary cyclic subgroups of the abstract finite group.

## 3. LH1: the actual Keller constraint removes most projective actions

Now let F=C(K) subset L=C(x_1,...,x_n) be an ACTUAL complex polynomial
Keller extension. Suppose z in L minus F has the full system (1).
Retained KMFI1 applies to E=F(z) as well as to L: NO NONTRIVIAL
INTERMEDIATE FIELD GALOIS OVER F can be contained in E. Equivalently,
the conjugates of the point stabilizer H normally generate G.
Use this together with the faithful projective action in Section2.

* H=1 is forbidden: E/F would be nontrivial Galois.
* A cyclic group has no faithful nontrivial orbit with a nontrivial
  point stabilizer, so the cyclic case is excluded.
* In D_m of order2m, the rotation-fixed orbit has degree2 and a normal
  stabilizer. The remaining possible H is a reflection, with degree m.
  Writing D_m=<r,s:r^m=s^2=1,srs=r^-1>, the normal closure of s is
  <s,r^2>. It is all of D_m exactly when m is ODD. This is retained KMFI2,
  and the formula handles EVERY m without a degree bound.
* In A4 the signature is (2,3,3). An order2 stabilizer normally generates
  only V4 and is forbidden. An order3 stabilizer normally generates A4,
  giving degree4.
* In S4 the signature is (2,3,4). The order2 branch stabilizer is generated
  by a transposition and normally generates S4, giving degree12. The
  order3 stabilizer normally generates A4 and is forbidden, giving no
  degree8 case. The order4 stabilizer contains an odd4-cycle and normally
  generates S4, giving degree6. The order2 subgroup of an order4 point
  stabilizer is NOT a different full point stabilizer.
* In A5, which is simple, each branch stabilizer of order2,3,5 normally
  generates A5. Its degrees are30,20,12 respectively.

Thus the ONLY NECESSARY surviving actions for E/F are

| Normal-closure group | Point stabilizer | Degree [E:F] |
|---|---|---|
| D_m, m>=3 odd | reflection C2 | m |
| A4 | C3 | 4 |
| S4 | C4 | 6 |
| S4 | transposition C2 | 12 |
| A5 | C5 | 12 |
| A5 | C3 | 20 |
| A5 | C2 | 30 |

These are necessary conditions, NOT actual realization or sufficiency.
The table works for Keller maps in EVERY dimension. In particular it
does not turn the known cubic three-dimensional example into a
contradiction: its S3 natural degree3 action is D3 and survives.

### The extra plane consequence requires a primitive generator

For an actual PLANE Keller map, if z is PRIMITIVE, F(z)=L, published
Domrina/KR0 excludes the degree4 row. Hence a noninvertible plane map
admitting such a projective primitive generator must have generic degree

    an odd integer >=3, or 6,12,20,30,                (6)

with the ACTION specified in the table. In odd degree its group must be
the natural dihedral group, not merely any group of that odd degree.

Every source polynomial COORDINATE is primitive. Here is the standard
coordinate-line field argument. Write R=C[f,g] and restrict (P,Q) to the
generic f-line over an algebraic closure of C(f). A nonconstant polynomial
pair in g is finite onto its image. The normalization of that image is
A1, so factor the parametrization through a polynomial h(g). The Keller
differential makes P_g,Q_g nowhere simultaneously zero. Therefore h'
has no zero and h is linear. The generic line map is birational, giving
L=F(f). This is the classical coordinate-field phenomenon (also underlying
the retained KW0/PF2RULED53 line argument); it is not a new primitive-field
theorem claimed by priority.

Consequently for a full degree or full monodromy action outside (6), NO
polynomial source coordinate, in ANY polynomial source frame, can obey
quadratic rational equations along the FULL target derivative frame.
More generally no rational primitive generator has that property. In
particular this excludes natural S_d for d>=4 and natural A_d for d>=5
as projective-generator cases. It does NOT exclude those groups as Keller
monodromy without the differential hypothesis.

DO NOT delete the degree4 row for a PROPER intermediate E inside a larger
plane source field. Its finite normalization is not supplied with a
polynomial plane source; Domrina's full-source theorem does not apply
to that intermediate field. No quartic realization queue is reopened.

There is also a useful link to LE0. For nontrivial z in an actual plane
field, no nonzero constant Nambu direction can make the quadratic
coefficient in (1) vanish: that would be an affine first-order equation
and LE0 would force z in F. Thus a_1,a_2 must be C-LINEARLY INDEPENDENT
in F. This is NOT F-linear independence, and it does not assert the
existence of any projective element. An arbitrary rational linear
combination of target derivatives is not silently treated as a constant
Nambu direction.

## 4. Sharp controls and the full-frame boundary

The universal cubic f(T)=T^3+A*T+B gives a control of the surviving D3
row. Its normal source S=C[A,t], B=-t^3-A*t, is finite free of degree3
over C[A,B]. The polynomial is irreducible because the presentation
eliminates B, and its discriminant -4A^3-27B^2 is not a square. Thus its
group is S3 with reflection stabilizer, and it has no nontrivial Galois
intermediate. Implicit differentiation gives

    delta_A t=-t/(3t^2+A), delta_B t=-1/(3t^2+A).

Reducing either inverse modulo f gives a quadratic in t over C(A,B),
and their curvature (4) is zero. Nevertheless its Jacobian is
-(3t^2+A), NONCONSTANT. Its entire etale source inverts 3t^2+A and is
NOT whole A2. This is not a plane Keller or HC4 counterexample.
Indeed every separable cubic primitive element has quadratic derivative
representatives, so the degree3 table must not be advertised as an
inverse theorem.

Using only ONE derivative is insufficient for the general LH0 theorem.
There is an every-degree control: for every integer m>=4 take

    P=s^m-s, Q=t, L=C(s,t), F=C(P,Q), z=s.

Then delta_Q z=0, but delta_P z=1/(m*s^(m-1)-1) is NOT quadratic over F.
If it were a*z^2+b*z+c, multiply by the denominator and reduce using
z^m=z+P. The result would be

    m*c*z^(m-1)+(m-1)*a*z^2
       +(m*P*a+(m-1)*b)*z+(m*P*b-c)=1.              (7)

For m>=4 the distinct degrees force c=a=b=0, a contradiction. The
presentation is finite free of EXACT degree m. This is an unbounded
written proof, not extrapolation from a few m. The one-direction constants
are C(s), larger than C(P). The full common constants are still C.
These maps also have NONCONSTANT Jacobian, so they do not contradict
the actual-plane LE0 theorem. No single-direction constant-field assumption
is omitted from that retained theorem.

## 5. Literature, proof scope and what remains

The Riccati cross-ratio/superposition principle and the finite projective
groups are CLASSICAL. Primary context read: Jacques-Arthur Weil,
*Introduction to Differential Algebra and Differential Galois Theory*,
CIMPA-UNESCO-Vietnam lectures, Hanoi, December3-7,2001. The actual title
and complete PDF pages26-28 (printed24-26), including C.1.2/C.1.3,
were rendered and visually read. The binary tetrahedral/octahedral/
icosahedral groups there have orders24,48,120; their PROJECTIVE groups
here have orders12,24,60. We do not conflate them. The lecture's usual
one-derivation algebraically closed-constant setting is not imported as a
theorem about C(P); the FULL-frame common-constant proof is supplied here.
The finite classification used here is also derived in Section2 from
Riemann-Hurwitz and spherical triangles. No new major external theorem
or unchecked differential-Galois classification is imported.

The exact checker controls the cross-ratio identity, interpolation,
curvature, rational Mobius gauge, projective group point-stabilizer normal
closures, the cubic and the all-m coefficient obstruction. The universal
proof and its field hypotheses are written above; finite group calculations
are not an all-degree search or an existence certificate. All retained
proof/checker/source bytes remain unchanged; no registered mathematical
checker is rerun. Independent replay, formal verification, external review
and original goal completion remain false.

This closes a possible unsupported extension of affine differential
rigidity: full quadratic systems constrain monodromy but do not by
themselves force descent. No actual rational or polynomial projective
primitive element has been constructed universally. The remaining odd
dihedral and specified polyhedral actions are not realized or excluded
here. Actual uniform integration, unaveraged complete-field descent,
trace-dual pairing regularity, and full HC4/JC2 remain OPEN.
