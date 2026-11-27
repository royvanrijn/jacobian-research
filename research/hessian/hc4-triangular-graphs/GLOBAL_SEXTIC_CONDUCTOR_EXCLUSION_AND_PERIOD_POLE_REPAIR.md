# Global sextic conductors are excluded by repairing the period poles

Let K=(P,Q) be an actual complex polynomial plane Keller map in a
source frame where P is monic of positive degree in x. Let

    (C[P,Q,y]:C[x,y])=kappa*C[x,y]

be its exact coordinate-order conductor. Then

    deg_x kappa != 6.                                    (1)

This includes source fibres where its degree drops to four or two.
Together with the retained global degree-two and degree-four
exclusions, every noninvertible map in this frame has

    deg_x kappa >= 8, and this degree is even.             (2)

These are GLOBAL directional degrees. An individual special fibre
of a higher-degree conductor may still have degree two, four or six.
Higher conductors, automatic first-velocity or trace integrality,
HC4 and JC2 remain open.

The proof uses the complete primary genus-three classification.
It does not assume a generic family of target automorphisms.
The three ambient types are separated by finite filtered-order
conditions and by algebraic conditions on their physical action
vectors. Monicity and a nonzero low-defect action then remove every
possible pole of the rational canonical parameters.

## 1. Retained inputs and the pointwise genus-three classification

The [coordinate-order theorem](GENERIC_COORDINATE_ORDER_CONDUCTOR_AND_CLOSED_COLLISIONS.md)
gives, on EVERY source line y=c, the exact finite birational order

    A_c=C[P(x,c),Q(x,c)] subset C[x],

with immersed normalization C[x], specialized conductor kappa_c*C[x],
and defect delta(c)=deg_x(kappa_c)/2. Under noninvertibility every
specialized conductor is nonzero of positive even degree.

Assume global degree six. On the parameter open where its leading
coefficient a_6(y) is nonzero, the defect is three. Such an immersed
affine-line-normalized plane curve has one place at infinity.
The corrected Neumann--Rudolph knot-fibre lemma makes its defining
polynomial good at infinity. Smoothing the finite singularities
gives a once-punctured genus-three fibre: the local normalization
and smoothing count gives Euler characteristic 1-2*delta=-5.

The [primary genus-three audit](sources/neumann-genus-three-source-audit.json)
retains Neumann's COMPLETE three realizable ambient types:

    intrinsic quartic, intrinsic degree-six cable,
    intrinsic degree-seven hyperelliptic.

The classification is pointwise, up to polynomial ambient
automorphisms and affine changes of the defining value.
The [hyperelliptic sextic theorem](HYPERELLIPTIC_SEXTIC_CONDUCTORS_AND_THREE_PERIOD_RIGIDITY.md)
already excludes its generic hyperelliptic branch, including contacts
and degree drops.

Two pointwise action facts will also be used:

- Every immersed rational quartic-type fibre has a nonzero vector
  of matching-endpoint physical actions, by the
  [quartic chord theorem](QUARTIC_GENUS_THREE_CHORD_ACTIONS_AND_PERIOD_EQUATIONS.md).
- Every immersed rational cable-type fibre has a nonzero such vector,
  by the [weighted cable theorem](CABLE_GENUS_THREE_WEIGHTED_PERIODS_AND_NONVANISHING.md).
  This is a statement about the PHYSICAL cable action, not a
  Jacobian-preserving change to its intrinsic cubic/quartic functions.

For defect one or two the analogous nonzero action-vector statements
follow from the unrestricted nodal theorem and Sections1,2,5 of the
[constant degree-four proof](CONSTANT_DEGREE_FOUR_CONDUCTORS_AND_PERIOD_ODE.md).
They apply to an individual fibre: a node, two nodes or a tacnode.
No exclusion of an isolated low-degree fibre is imported.

## 2. Finite filtered descent replaces a family-type assumption

Put R=C[x,y], A=C[P,Q,y], k=C(y). Localize C[y] at a_6. The
normalized conductor is now monic of degree six, so R/(kappa) is
finite free of rank six over that localized PID.

The specialized inclusion A/(y-c) into R/(y-c) is injective by the
coordinate-order theorem. Thus R/A is torsion-free over C[y].
After this localization it is finite, hence free, of rank three.
The quotient maps from polynomials of source degree at most r,

    span{1,x,...,x^r} -> R/A,

are therefore finite matrices over the localized parameter ring.
Their generic kernels and their dimensions on a parameter open
are determined by ordinary linear algebra.

Every pointwise genus-three type has a nonconstant source polynomial
of degree at most three in its exact order. The hyperelliptic type
has one of degree two. The quartic and cable types have filtered
semigroup <3,4>: the cable square identity gives intrinsic generators
of degrees three and four, without asserting an ambient automorphism.

Consequently the generic least nonconstant source degree is two
or three; degree one would make the defect zero. If it is two,
the same finite matrix ranks give this least degree on a nonempty
parameter open. Neither quartic nor cable admits a degree-two
element, so that open is hyperelliptic, already excluded.

In the remaining case, on a further open the kernel in source degrees
at most three has dimension two and the kernel in degrees at most
four has dimension three: these are the filtered dimensions of both
remaining primary types. The matrix ranks therefore supply an
independent rational cubic and quartic. Normalize their
leading coefficients, depress the cubic by a rational source shift,
and remove the quartic cubic and constant terms:

    u=x-h(y),
    p=u^3+A*u, q=u^4+B*u^2+C*u,
    A,B,C,h in k.                                      (3)

The suborder k[p,q] has colength three: its k[p]-basis is
1,q,q^2 and the leading degrees 3i+4j, j=0,1,2, have distinct
residues modulo three. The missing degrees are 1,2,5.
Since k[p,q] is contained in A_k and both have colength three,

    A_k=k[p,q].                                         (4)

No coefficient in(3) is yet declared regular at a finite parameter.

The monic conductor of(3), calculated in the quartic chord theorem, is

    kappa_0=u^6+(A+B)u^4+2C*u^3+(A^2-AB+B^2)u^2
                    +C(A+B)u+A(A-B)^2+C^2.              (5)

Thus the actual conductor is a_6(y)*kappa_0(x-h), and
h=-[x^5]kappa/(6a_6).

## 3. Pair actions and algebraic selection of the period equations

Matched distinct parameters alpha,beta of(3) have sum sigma,
gap Delta=beta-alpha, and satisfy

    g(sigma)=sigma^3-(B-2A)*sigma-C=0,
    Delta^2=-K(sigma), K=3sigma^2+4A.                    (6)

Set

    U=2A-3B, V=2B-A, W=4A-3B,
    H=15sigma^2+6A-7B,
    R0=4A*U^2+27C^2=Res(g,K).                           (7)

The rational function R0 is nonzero by immersion.

For later use, persistent multiple roots of g can be disposed of
without assuming a family target automorphism. A double root
sigma!=0 gives parallel canonical tangent vectors whose first
components differ by 3*sigma*Delta!=0. Any two plane generating
pairs of the SAME order have an invertible change of tangent
coordinates at a singular image point. Indeed the inverse
polynomial substitutions compose to the identity modulo the
curve equation, whose gradient vanishes at that point.

For an actual persistent matched pair, parallel source tangent
vectors must be EQUAL. Differentiating the endpoint identities
along y and taking determinants shows that their proportionality
factor times the nonzero constant source Jacobian equals that
same Jacobian. The factor is one. The canonical tangent difference
therefore gives a contradiction.

The triple-root case is sigma=0, B=2A,C=0. It cannot be cable by
the cable normalization constraints below, and is the already
excluded persistent order-three-contact quartic branch. Thus g
may be assumed generically separable.

Let I_sigma be the actual matching-endpoint action integral P*dQ.
Its quotient by Delta is a symmetric divided difference of a
polynomial primitive in u, so it belongs to k[sigma]/(g).
Its three reduced coefficients form a vector over k.

At a quartic-type parameter that vector is proportional to

    R_Q(sigma)=K^2*H.                                   (8)

At a cable-type parameter the normalized coefficients necessarily
satisfy

    B=2A/3, C!=0,

and the PHYSICAL action vector is proportional to

    R_C(sigma)=C^2*(8A*sigma+9C).                        (9)

All nonzero fixed numerical factors are absorbed in the scale.
These precise statements are proved in the two pointwise gates.

Proportionality of two nonzero three-component vectors is the
vanishing of their two-by-two minors, an algebraic condition over
k. On a parameter open, the quartic proportionality condition
and the cable coefficient/proportionality conditions are CLOSED
conditions and cover all complex parameter points by the primary
classification. Restrict to C!=0 if C is not identically zero;
if C is identically zero the cable set is empty.

An irreducible parameter open cannot be the union of two proper
closed subsets. Hence one of these conditions holds identically.
Comparison of any nonzero reduced coefficient descends its common
scale Lambda to k, directly. This proves the following exhaustive
generic alternatives without a generic-family classification
assumption or descent of target automorphisms.

In the QUARTIC period alternative,

    I_sigma=Lambda*Delta*K^2*H,
    12V*Lambda'+42V'*Lambda=-3J*U/R0,                    (10)
    45C*Lambda'+105C'*Lambda=-9J*C/R0,
    (A-B)*(2W*Lambda'+7W'*Lambda)=0,

where J in C* is the actual source Jacobian.
If A!=B and W!=0, Lambda^2*W^7 is in C*.
The separate case W=0 remains allowed.

In the CABLE period alternative, B=2A/3,C!=0 and

    I_sigma=Lambda*Delta*C^2*(8A*sigma+9C),
    9C^3*Lambda'+30C^2*C'*Lambda=-J.                     (11)

For completeness, differentiating(9), using
sigma'=(3C'-4A'*sigma)/(9sigma^2+4A), and reducing modulo g
gives(11). It also gives

    4AC*Lambda'+10AC'*Lambda+5A'C*Lambda=0.

These are rational differential identities. No constant resultant
or absence of parameter poles has been assumed.

## 4. Constant fibre degree six is impossible

First suppose a_6 is a nonzero constant, and make kappa monic.
Then h is polynomial. In the centered conductor(5),

    S=A+B, T=A^2-AB+B^2, C

are polynomial, by its coefficients. Also AB=(S^2-T)/3 is
polynomial. The rational A,B are integral roots of

    z^2-S*z+(S^2-T)/3,

so they are polynomial by normality of C[y].

This polynomial descent also identifies EVERY specialized order.
Put Tcal=C[y,p,q]. Over C[y,p], it has basis 1,q,q^2.
The monic leading terms of p^i*q^j, j=0,1,2, have distinct
degrees 3i+4j. Leading-coefficient subtraction therefore proves

    Tcal_k intersection C[y,u]=Tcal.                    (12)

In detail, the highest coefficient of a polynomial belonging to
Tcal_k determines its unique basis monomial. It lies in C[y];
subtract that polynomial multiple and continue downwards.
The same finite leading-degree argument works in every fibre.

Actual P,Q are polynomial and belong to Tcal_k, hence to Tcal.
Each specialized Tcal_c has colength three with gap basis u,u^2,u^5.
The actual order A_c is contained in it and has the SAME colength.
Thus A_c=Tcal_c for every complex c. Actual immersion now implies

    R0 in C*.                                           (13)

Indeed a common zero of p_u,q_u would force P_u=Q_u=0.

Every such fibre order has no nonconstant element of degree two,
so no special fibre can be hyperelliptic. Each is quartic or cable,
and its actual physical action vector is nonzero by the pointwise
gates. The canonical vector in either generic alternative is also
nonzero at every finite parameter: for(8) this is the quartic gate;
for(9), equation B=2A/3 and(13) give R0=27C^2, so C is a unit,
and the cable gate applies.

The roots of monic g remain finite on each local parameter cover.
Equation(13) makes every pair gap a unit, including repeated
collision sums. The polynomial actual action integrals have finite
limits. A finite pole of the rational Lambda would contradict a
nonzero canonical limiting action. A finite zero of Lambda would
make ALL actual limiting actions zero, also impossible. Therefore

    Lambda in C*.                                       (14)

In the quartic alternative, the second equation of(10) now says
C'=k*C with a nonzero CONSTANT k. A polynomial solution is C=0.
Then R0=4A*U^2 is a unit, so A,U and hence B are constants.
The first equation of(10) has zero left side and nonzero right
side, a contradiction.

In the cable alternative C is constant nonzero, and(11) gives
Lambda'=-J/(9C^3)!=0, contradicting(14).
This proves the constant-fibre-degree case.

## 5. A finite unit action is available at every proposed degree drop

Suppose a_6 has a finite zero y=c. The specialized conductor has
degree TWO or FOUR, never zero or the zero polynomial.
Choose a special matching pair with NONZERO action, using the
pointwise low-defect inputs in Section1.

Lift a nonvertical branch of the conductor through one endpoint
after a finite Puiseux substitution. Such a branch exists because
a vertical conductor component would make the whole specialized
conductor zero. Its generic mate has the same P value. Monicity
of the ACTUAL P bounds that mate: it cannot escape while the
first endpoint stays finite.

The two endpoints cannot merge. Their identical P,Q values,
divided by their nonzero parameter difference, would give
P_x=Q_x=0 at a merged limit. The special node or tacnode has
exactly one distinct mate for the selected endpoint; thus the
generic mate tends to the selected OTHER endpoint. If the
special defect-two fibre has two nodes, choose the one with
nonzero action.

Consequently some generic matched pair has, over the local
Puiseux field,

    x_1,x_2 finite, x_2-x_1 a unit,
    I=integral_(x_1)^(x_2) P*dQ a unit.                  (15)

This is a NONZERO limiting action, not only a generically
nonzero one. It will control the period scale at the drop.

The ordinary triple-point alternative A=B can already be
removed. Here kappa=a_6*((x-h)^3+A*(x-h)+C)^2. All three
roots of the monic cubic are generic mates. Starting with
one finite conductor branch, monicity brings ALL three
mates to finite points. None can merge by the same divided
difference argument. The cubic therefore has regular
coefficients; a_6=0 would make the entire conductor fibre
zero. Thus no such zero exists, and the constant case
already excludes this alternative. Hence A!=B below.

## 6. Quartic period poles cannot survive a degree drop

Use the unit pair(15), set v=y-c, and put

    sigma=x_1+x_2-2h, K=3sigma^2+4A=-(x_2-x_1)^2.

Thus K is a unit. Suppose h has a pole of positive INTEGER order n.
Then sigma~-2h and A~-3h^2.

Write H in the useful exact form

    H=(35/2)*sigma^2-(5/6)*K+(7/3)*W.                    (16)

If W!=0, let k=-ord_v W and l=ord_v Lambda. The first integral
in(10) gives l=7k/2. Since I, the gap and K are units, l is
exactly the pole order of H.

If k>2n, equation(16) gives l=k, contradicting l=7k/2.
If k=2n, its pole order is at most 2n but l=7n, also impossible.
If k<2n, its pole order is EXACTLY 2n. Hence

    l=2n, k=4n/7.                                       (17)

If W=0 identically, the same l=2n follows directly from(16);
no first integral is asserted or needed in that case.

In every surviving case W is of smaller order than h^2, so

    A~-3h^2, B~-4h^2, C~-4h^3.                          (18)

The three roots of g, divided by h, have leading values
-2,1+i,1-i. The selected root has K a unit. At the other
two, K has exact pole order 2n with nonzero leading
coefficients. Thus

    ord_v R0=-4n.                                       (19)

Divide the second equation of(10) by C, which is nonzero in(18):

    45Lambda'+105(C'/C)Lambda=-9J/R0.

By(17),(18), its left side has nonzero leading coefficient
(90n-315n)*lambda_0 at order 2n-1. Its right side has order
4n. This is impossible. Therefore h has no pole at this drop.

Now sigma is regular, so K=3sigma^2+4A makes A regular.
If B had a pole of order k>0, W and H would both have that
pole. The unit action would give ord Lambda=k, while the
nonzero-W first integral gives 7k/2. This is impossible.
The case W=0 makes B=4A/3 regular directly. Thus B is regular,
and C=sigma^3-(B-2A)sigma is regular as well.

## 7. Cable period poles cannot survive a degree drop

Again use(15), with K=3sigma^2+4A a unit. In this alternative
B=2A/3 and the pair cubic gives the EXACT identity

    C=K*sigma/3.                                        (20)

If h has a pole of positive integer order n, sigma~-2h,
A~-3h^2, and C has pole order n. The physical cable action
factor can be rewritten exactly as

    C^2*(8A*sigma+9C)=K^2*sigma^3*(5K-6sigma^2)/9.

Its pole order is five n. The unit actual action therefore
forces ord Lambda=5n. In(11) the left side has nonzero leading
coefficient (45n-30n)*c_0^3*lambda_0 at order 2n-1.
The right side is the unit -J. Since n is a positive INTEGER,
2n-1>0, a contradiction.

Thus h is regular. Then sigma and K are regular,
A=(K-3sigma^2)/4 is regular, and(20) makes C regular.
The same holds for B=2A/3.

## 8. Global polynomial descent finishes the exclusion

Sections6,7 repair h,A,B,C at every zero of a_6. At any other
finite parameter, h and the normalized conductor coefficients
are already regular. Formula(5) gives C, A+B and A^2-AB+B^2
regular there; monic quadratic integrality makes A,B regular.
Thus all four rational functions h,A,B,C have no finite pole
and are POLYNOMIAL.

But now kappa=a_6*kappa_0(x-h) is a polynomial multiple of a
monic polynomial with polynomial coefficients. Any finite
zero of a_6 would make kappa_c identically zero, which the
exact source theorem forbids. Hence a_6 is a nonzero constant.
Section4 gives the final contradiction.

This exhausts all three pointwise genus-three types, their
persistent contact alternatives, every degree drop to defect
one or two, and every possible finite canonical parameter pole.
It proves(1). Combining it with the prior global degree-two/four
exclusions and generic even positivity proves(2).

## 9. Evidence and the remaining goal

The exact checker supports the monic conductor identities, degree
filtration and index, unit immersion resultant, physical cable
differentiation with variable coefficients, both pole-order
obstructions and the constant-case contradictions.
The unrestricted filtered descent, pair-lifting and valuation
arguments are the written proof, not finite parameter sampling.

The primary genus-three classification, corrected knot-fibre lemma
and low-defect pointwise action results remain explicit inputs.
No formal verification, independent whole-proof replay or external
human review is claimed.

This sextic exclusion gave even global coordinate-conductor degree at
least eight. The [current lower bound](GLOBAL_CONDUCTOR_TWENTY_EXCLUSION.md)
is twenty-two. Automatic trace/codifferent integrality, HC4 and JC2 remain open.
