# Kulikov's exact coframe and the plane-boundary weight filter

Claims HC4MYKT0/1. This is a sourced reconstruction of a KNOWN non-plane
counterexample to weakened Jacobian statements, with an explicit cohomology
comparison. It is not an HC4 counterexample or a proof of JC2. Full HC4/JC2
remains open.

## 1. Credit and the question this tests

Vik. S. Kulikov, *Generalized and local Jacobian problems*, Russian Acad.
Sci. Izv. Math.41:2 (1993),351--365,
[DOI](https://doi.org/10.1070/IM1993v041n02ABEH002265), constructs in
Section3 a simply connected surface with a nonbirational etale map to A2.
Printed pp357--358 identify its compactification with the blowup of P2
at one point and give the conic net used below. Adjamagbo,
[*The complete solution to Bass Generalized Jacobian Conjecture*](https://arxiv.org/abs/1210.5281),
pp3--5, records the same construction and its affineness, factoriality and
constant units. The construction and those properties are not new here.

The purpose is to test the proposed replacement of a polynomial source by
normality, factoriality, rationality, constant units, an exact differential
coframe, or quadratic trace-zero reconstruction. ALL these tests can hold
on an actual non-plane etale source. A global pair P,Q with
dP wedge dQ nowhere zero is therefore not by itself an obstruction.
The actual polynomial-source or plane-boundary hypothesis must be retained.

## 2. An explicit member of the published construction: HC4MYKT0

In homogeneous coordinates X,Y,Z put p=(1:1:1),

    Q_i = 3*X_i^2-X*Y-X*Z-Y*Z,   (X_1,X_2,X_3)=(X,Y,Z),
    R = sum_(i!=j) X_i^2*X_j - 6*X*Y*Z,
    C = Q_1+2*Q_2+4*Q_3.

The net defines a rational map [Q_1:Q_2:Q_3]. Its only base point is p:
the equations force X_i^2 equal, and the four possible projective sign
patterns leave only (1:1:1). Its linear terms span the cotangent space at
p, so one blowup resolves the net without an infinitely-near base point.
Let W=Bl_p(P2), with H the pullback line class and E the exceptional curve.
The resulting morphism f:W->P2 has

    f^*O(1)=2H-E,       (2H-E)^2=3.

The divisor 2H-E intersects E and every other irreducible curve positively:
for a curve of degree a and multiplicity b<=a at p, its intersection is
2a-b>0. Thus f contracts no curve. Properness makes it finite, and its
degree is THREE. This is a map of a rational surface, not a polynomial
plane Keller map.

The exact identity det(D(Q_1,Q_2,Q_3))=-36*R gives the ramification cubic.
Its unique singular point is p, an ordinary node. In local coordinates
X=1+a,Y=1+b,Z=1 the equation is

    a^2*b+a*b^2+2*a^2-2*a*b+2*b^2.

The nondegenerate quadratic part and the homogeneous gradient ideal
certify the node and uniqueness. A reducible cubic with an ordinary node
would have another intersection singularity; hence R is irreducible.
Its strict transform R' is a smooth P1 of class 3H-2E. The Jacobian
divisor on W is exactly R', with generic ramification index2. The conic
C is nonsingular: its symmetric matrix has determinant -127.

The strict transforms R',C' meet transversely at FOUR distinct points.
Here is an exact normalization certificate. Parametrize the nodal cubic
by b=t*a. Away from its node,

    a=-(2*t^2-2*t+2)/(t^2+t).

After removing the node factor, C|R has numerator

    -6*(2*t^4-4*t^3+7*t^2-2*t+1).

It is square-free, has degree4, and has no common root with the denominator.
The tangent direction of C at p is not a node direction, so no strict
intersection is hidden on E; the point t=infinity is not an intersection.
The [exact control](verify_kulikov_tate_control.py) checks these assertions.

Let L be the target line with equation U_1+2U_2+4U_3=0, and set

    Z_0=W-C'=f^-1(P2-L),
    S=W-(C' union R').

Z_0 is smooth, affine and finite of degree3 over A2=P2-L. Removing the
Cartier curve R' keeps S affine. The restricted map

    g:S -> A2,      P=Q_1/C, Q=Q_2/C

is ETALE of degree3. Ratios here mean the resolved sections on W; the
exceptional common factor is cancelled. Thus dP,dQ are an exact global
basis of Omega^1_S, and d(P*dQ)=dP wedge dQ is a nowhere-zero EXACT volume.
The surface is rational. Its coordinate ring A is factorial and A*=C*:
Pic(W)=Z H plus Z E, and the boundary classes (3,-2),(2,-1) form a
unimodular basis. Localization gives Pic(S)=Cl(S)=0. A unit would have
divisor supported on those two independent boundary classes, hence zero,
and is constant on the projective W. Kulikov additionally proves simple
connectedness for the transverse-line construction; that is a credited
published input, not re-proved by the finite algebra check.

The full finite normalization over the target is EXACTLY Z_0. Its one
omitted boundary curve is

    R' minus C' = P1 minus FOUR points.

It is not A1. It has e=2,f=1, leaving one actual unramified sheet over
the generic branch curve; monodromy is S3. In Cl(Z_0), E=2H, so

    [R']=-H,       [K_(Z_0)]=-H=[R'].

The canonical class is nonzero primitive. Thus the model also passes the
retained nontorsion-canonical and fixed-sheet necessary tests. It does not
evade KN by a torsion-canonical finite normalization.

### Even the quadratic trace-zero reconstruction holds

Let B=C[P,Q], L_0=C(S), and N=ker(Tr_(L_0/Frac(B)):A->Frac(B)).
The same localization/finite-IC argument as KQ applies to the SMOOTH pair
(Z_0,R'-C'): it does not need that this particular boundary curve be A1.
There is a semisimple full-support socle I=B plus I_0 and simple curve top
T=IC_branch. The S3 augmentation I_0 is simple of rank2. Trace on the top
is nonzero, so N=I_0. Hence N is SIMPLE and the source module has length3.

The quadratic image N*N over B is W-stable, has full generic rank3 and
therefore contains I, by the retained trace-hyperplane argument. In the
omitted ramified factor t^2=u, the square of the trace-zero pole t^-1 is
u^-1 and has nonzero boundary class. The simple top is therefore included,
and

    A = B[N] = image(Sym^2_B N -> A).

Consequently exact coframe plus factoriality plus constant units plus
simple trace-zero kernel plus complete quadratic reconstruction is STILL
not an unconditional affine-plane recognition theorem. KS0 correctly
assumed an ACTUAL plane-source open; this example does not contradict it.

## 3. The cohomology which the weakened tests miss

Write D=R' union C'. It is two P1 components meeting at four nodes.
Its dual graph has first Betti number3, hence H1(D,Q)=Q^3. Restriction
H2(W,Q)->H2(D,Q) has matrix

    [3 2]
    [2 1]

in the H,E and component bases, with determinant -1. The compact-support
exact sequence for S=W-D therefore gives H_c^2(S)=Q^3 and H_c^3(S)=0;
connectedness gives H_c^1(S)=0. Poincare duality yields

    H0(S,Q)=Q, H1(S,Q)=0, H2(S,Q)=Q^3,
    Hk(S,Q)=0 for k>2,       chi(S)=4.

Thus S is not A2. It cannot even CONTAIN an open A2. If such an open
existed, any boundary divisor would be principal by factoriality, with
its defining function a nonconstant unit on that plane, impossible.
With no boundary divisors, plane coordinates extend by normal Hartogs
and separatedness makes the open embedding an isomorphism, again impossible.

The weight is also explicit. Put q=uv in the compact-support Hodge--Deligne
polynomial. The newly screened standard MHS input is in the
[source context](sources/kulikov-tate-context.json): functorial mixed Hodge
structures, compact-support additivity and smooth Poincare duality.
Then

    E_c(W)=1+2q+q^2,
    E_c(D)=2(1+q)-4=2q-2,
    E_c(S)=q^2+3.

Only H_c^4=Q(-2) and H_c^2 contribute. The latter is therefore Q(0)^3,
and duality gives

    H2(S,Q)=Q(-2)^3: pure weight FOUR, type(2,2).       (1)

This is a boundary-graph contribution, not a Picard class: Pic(S)=0.
The standard sign/matrix computations here are independent of incidental
typographical examples in the expository MHS source.

## 4. An all-degree plane-completion filter: HC4MYKT1

For EVERY actual plane completion V used in KS0, or the maximal smooth
etale completion, the added boundary is c disjoint A1 curves over the
plane X=A2. The ordinary Gysin sequence gives

    H1(V,Q)=0, H2(V,Q)=Q^c, Hk(V,Q)=0 for k>2.

The Gysin generators are precisely the c boundary divisor classes. In
mixed Hodge structures,

    E_c(V)=q^2+c*q,
    H2(V,Q)=Q(-1)^c: pure weight TWO, type(1,1),
    c1:Pic(V) tensor Q -> H2(V,Q(1)) is an isomorphism. (2)

These follow from compact-support additivity and duality just as above;
the only nonzero ordinary cohomology groups are degrees0 and2, so there
is no cancellation hiding other Hodge types. This statement is independent
of the generic sheet degree and uses the actual disjoint-A1 boundary.

Equations(1) and(2) separate Kulikov's exact-coframe example from EVERY
such actual plane completion, even when the unweighted b2 ranks agree.
The boundary punctures, not a missing algebraic coframe, are the exact
failure in this model. Products with affine space retain this distinction.

This filter is AUTOMATIC for actual plane sources and is useful to reject
relaxed replacement algebras. It does not newly exclude the retained
one-boundary or higher-degree plane packets, which already satisfy their
own source-boundary hypotheses. Higher-exponent Danielewski completions
from UADANVOL1 also contain A2 charts and pass the weight-two test;
no etale coordinate pair on them is supplied. Another unmarked cohomology
budget is not a new HC4 obstruction. The remaining useful requirement is
an actual polynomial source or a constraint on its marked energy/field
multiplication, retaining the full source topology and boundary geometry.
