# An all-degree power-difference obstruction and the essential cubic index

15 November 2026. PF2RULED56 proves a general obstruction without a
coefficient bound: if p,q>=2 are coprime and an irreducible polynomial
G=A^p-B^q has zero curve with normalization A1, then A and B have an
affine common zero. Applied to the actual quartic's normal cubic
resolvent, this excludes **every global monogenic model** of that
algebra. All its primitive integral elements have a nonconstant index.
This is a restriction on the actual finite algebra, not a proof that
all the necessary Picard 3-torsion vanishes or that no quartic exists.

## 1. The power-difference theorem

Work over C, with A,B,G in C[X,Y], G=A^p-B^q irreducible, p,q>=2
coprime, and normalization of {G=0} equal to A1. Then

    V(A,B) is nonempty.                              (1)

Assume the contrary. On the normalization parameter t, the equality
A(t)^p=B(t)^q and coprimality imply A(t)=a*h(t)^q,
B(t)=b*h(t)^p, for nonzero constants a,b with a^p=b^q.
The polynomial h has no zero, so it is a nonzero constant. Consequently
there is a fixed smooth point P=(a0,b0), a0*b0!=0, of the torus curve

    D={f(u,v)=u^p-v^q=0}

such that H=(A,B) maps the entire curve {G=0} to P. In particular

    A=a0+G*R,   B=b0+G*S.                            (2)

The divisibility follows in the reduced curve ring C[X,Y]/(G), not
just pointwise at generic normalization parameters. We treat both
possible dimensions of the image of H; dominance is not assumed.

### 1a. A one-dimensional image is impossible

Suppose A,B are algebraically dependent. The normalization of the
affine image curve lifts the map from normal A2: its integral
coordinate functions lie in the integrally closed ring C[X,Y].
Restriction to a suitable affine line gives a dominant map from A1
to that normalized curve, so its function field is rational by
Luroth's theorem. The normalized curve is affine and its units inject
into C[X,Y]^*=C^*. A rational affine curve with r punctures has unit
rank r-1. Thus it has exactly one puncture, is A1, and

    A=P(h), B=Q(h) for some polynomial h in C[X,Y].

The polynomial P(s)^p-Q(s)^q is nonconstant. Since substituting h
gives irreducible G, its factorization over C can have just one
simple linear factor. Therefore

    P(s)^p-Q(s)^q=c*(s-s0), c!=0.                   (3)

The three terms in (3) are pairwise coprime. A shared zero of P,Q
would give a multiple zero of the linear right side. Mason-Stothers
gives, writing a=deg P, b=deg Q,

    max(p*a,q*b,1) <= a+b.

This is impossible: adding a<=(a+b)/p and b<=(a+b)/q contradicts
1/p+1/q<1 when a+b>0. If both degrees vanish, (3) is impossible too.
There is no assumption that H has two-dimensional image hidden here.

For completeness, the polynomial abc bound used here follows from
the elementary Wronskian proof in Snyder. If U+V=W are nonzero
pairwise coprime polynomials with some derivative nonzero, the nonzero
Wronskian U'V-UV' is divisible by
gcd(U,U')*gcd(V,V')*gcd(W,W'). Taking degrees bounds deg W by
deg rad(UVW)-1. Permuting the three terms gives the maximum bound.
This is a recalled classical proof, not a new arithmetic input.

### 1b. A dominant image is excluded by one logarithmic pluriform

Now assume H is dominant. The Taylor expansion of f at P and
f(H)=G imply, using (2),

    1 = f_u(P)*R + f_v(P)*S + G*(terms in R,S).

Hence (R,S) is the unit ideal and H^*I_P=(G). The universal property
of the blowup lifts H to a morphism

    Htilde:A2 -> Y=Bl_P(A2) minus the strict transform of D. (4)

Indeed the total transform of D is Dstrict+E, while both H^*D and
Htilde^*E are the divisor of G, so Htilde^*Dstrict=0 and its preimage
is empty. The exceptional curve E is retained in Y; it is not removed
with Dstrict. Smoothness of P and multiplicity one of G are essential.

Set N=p*q, k=p+q. Since p,q are coprime and >=2, N>k. We exhibit a
nonzero global logarithmic N-canonical form on a smooth completion
of Y with simple-normal-crossing boundary:

    theta = (du wedge dv)^N / (u^p-v^q)^k.            (5)

Here the superscript denotes a tensor power. A logarithmic
N-canonical form may have pole order at most N on each reduced
boundary component, and must be regular on the interior. The
following valuations check **all** boundary resolutions, not a
bounded list of fans.

At the affine origin resolve the binomial cusp by a regular toric
subdivision containing the primitive ray (q,p). For an exceptional
ray (a,b), a,b>0, the form numerator has order N*(a+b-1), and f has
order min(p*a,q*b). Its order plus N is

    N*(a+b)-k*min(p*a,q*b)
      = p*(q*b-p*a) if p*a<=q*b,
      = q*(p*a-q*b) if q*b<=p*a.                    (6)

Both expressions are nonnegative. Coprimality makes the strict
binomial transform smooth and transverse on the critical ray; all
other toric boundary crossings are normal crossings after a regular
subdivision. The strict D itself has pole order k<N. The coordinate
axes away from D are interior and (5) is regular there.

Interchange the two variables and exponents if needed so q>p. In
the infinity chart u=1/w, v=s/w, formula (5) becomes

    (ds wedge dw)^N * w^(q*k-3*N)/(w^(q-p)-s^q)^k.

The unique meeting point of the projective binomial and infinity
is resolved by a regular toric subdivision containing (q-p,q).
For a positive ray (a,b) in (s,w), its order plus N equals

    N*a+(q*k-2*N)*b-k*min((q-p)*b,q*a)
      = p*(q*a-(q-p)*b) if (q-p)*b<=q*a,
      = q*((q-p)*b-q*a) if q*a<=(q-p)*b.            (7)

Again both are nonnegative. On the generic line at infinity the
order is q*k-3*N>=-N, since q>=p. These subdivisions resolve the
entire boundary: the affine origin and the unique infinity point
are the only singular or nontransverse points before the blowup P.

At the additional blowup of the smooth point P, the Jacobian of
the blowup has order1 and f has order1 on E. Formula (5) has order
N-k>=0 there. Therefore it is regular along the retained interior
exceptional curve, as required in (4). This is the point for which
simply asserting a formula on A2 minus D would be insufficient.

A dominant morphism pulls logarithmic top forms back to logarithmic
top forms after resolving its compactification. The local reason is
that d(u)/u has at most a simple pole along a boundary valuation;
a wedge of two such forms still has at most a simple pole. Tensor
powers have at most N. This is the differential inclusion used in
Fujino Lemma2.12. It applies after a graph resolution whose added
source boundary lies over the line at infinity, since (4) is a
morphism on the whole affine source. No surjectivity of (4) is needed:
the extension between projective compactifications is dominant and
therefore surjective. Characteristic zero makes the pulled top form
nonzero.

One can see the contradiction directly on the original source.
The expansion (2) shows G divides J=Jac(A,B). Thus

    Htilde^*theta = G^(N-k)*(J/G)^N*(dX wedge dY)^N   (8)

is a nonzero polynomial N-canonical form. Such a form, with coefficient
of degree d>=0, has pole order 3*N+d at the generic line at infinity
of P2. A logarithmic N-canonical form can have pole order at most N.
Resolving finitely many source boundary points cannot change this
generic order. This contradiction proves (1).

This proof uses a general logarithmic obstruction, not a classification
of possible coefficients. Equations (6)-(7) prove the required
inequalities for every positive valuation ray and every coprime p,q.

## 2. The actual cubic resolvent is locally monogenic everywhere

Retain the full actual quartic packet and its S3 Galois quotient T
from [PF2RULED54](QUARTIC_COLOURING_AND_ETALE_TRIPLE_GATE.md). The
degree3 quotient Z=T/<reflection> has normal coordinate algebra C3
over A=C[X,Y]. This is a subfield of the S4 Galois closure, not a
cubic subfield of the degree4 source field. We use C3 here for an
algebra; the cyclic group of order3 is written Z/3 below.

The triple T->S is etale at every affine sign point. Over a completed
strict local target ring R along C, T consists of three copies of
S_R=R[z]/(z^2-G). A reflection fixes one copy and exchanges the other
two. Taking its invariants gives

    C3 completed locally = R x R[z]/(z^2-G).         (9)

The same description holds at a tangency value of C. The two local
S4 fold transpositions were disjoint, so their S3 images are the SAME
reflection. It is precisely this fact that prevents a local S3 image
or a triple-root fibre there. Off C, the algebra is etale of rank3.
Consequently C3 is finite locally free of rank3, all C fibres have
length partition2+1, and its trace discriminant divisor is exactly
the reduced curve C. These conclusions also follow directly from
the local product (9), whose trace discriminant is4G.

Each such local algebra is monogenic: near C, the element (1,z)
has polynomial (T-1)*(T^2-G), and the two factors are coprime where
1-G is a unit. Off C choose distinct residues in its three geometric
factors and use faithful flat descent/Nakayama for a local generator.
This local assertion does not provide one element generating globally.

## 3. No global generator; every primitive element has a divisor index

Suppose C3=A[t] for a globally integral generator. Translate by
Tr(t)/3 to write its characteristic polynomial as

    F(T)=T^3+c*T+d,   disc(F)=-4*c^3-27*d^2=u*G,

with u in C^*. Its fibre never has a triple root: by (9) its partition
along C is2+1 and off C it is reduced. Thus V(c,d) is empty. Rescale
c,d by nonzero complex constants to express G=A0^2-B0^3, with
V(A0,B0) empty. This contradicts Section1 with p=2,q=3 and the
actual normalization Ctilde=A1. Hence

    C3 is NOT globally monogenic.                   (10)

This rules out every monogenic normal cubic resolvent throughout the
surviving quartic packet, with no restriction on rulings, source pole
orders or target image degrees. Neither a free module basis nor the
primitive-element theorem for its fraction field changes (10).

There is a precise basis-independent index formulation. The trace
splits C3=A direct-sum E, E=ker(Tr), a rank2 vector bundle. The
determinant line det(C3) is trivial because Pic(A2)=0. Fix a generator
of it, and for t in C3 define the polynomial index I(t) by

    1 wedge t wedge t^2 = I(t)*(chosen determinant generator).

Locally, for a basis e1,e2 of E, I(u*e1+v*e2) is a binary cubic
in u,v. It is unaffected by adding a scalar in A to t. Its value
is a unit if and only if 1,t,t^2 is a global basis, equivalently
if t is a global algebra generator. For every t in C3 outside A,
the fraction-field degree is3, so I(t)!=0. Equation(10) therefore gives

    I(t) is nonconstant for EVERY t in C3 minus A,
    disc(F_t) = constant * G * I(t)^2.               (11)

The second identity is the trace Gram determinant under the basis
change. Its index-square factor is unavoidable. Its support may
coincide with C; no new distinct branch curve is asserted. Since I
is a nonunit polynomial, some target divisor witnesses the failure
of that order to be the normal cubic algebra. A chosen primitive
resolvent polynomial must retain its index/conductor data. Removing
them changes the algebra being tested.

## 4. Controls that distinguish the hypotheses

The smooth Oka branch G=(xu-2)^2-u^3 has V(xu-2,u) empty and a
monogenic cubic with discriminant -108G, as retained in PF54. Its
normalization is C*, with two ends, so it does not contradict (1).

Conversely a genuine A1-normalized immersive torus curve does exist:

    G0=x^2-(x+y^2)^3,
    x(t)=(1-t^2)^3, y(t)=t*(1-t^2).

Here x+y^2=(1-t^2)^2. Away from its one double value (0,0), the
inverse parameter is t=y*(x+y^2)/x. The parameter is integral, since
t^3-t+y=0. The derivatives of x(t),y(t) are coprime, so the
normalization is immersive. Only t=1,-1 are identified, and their
branches have contact3: writing h=1-t^2, one has x=h^3 and
y=t*h, so x/y^3 tends to1 and-1 on the two smooth branches. Thus
the affine singularity is A5.

Irreducibility here is over C, not just over Q: in coordinates
u=x+y^2 the parametrization has u=(1-t^2)^2, of degree4 over C(u).
The rational inverse shows C(u,y)=C(t), so the monic quartic
(u-y^2)^2-u^3 is the minimal polynomial of y over C(u).

The monic cubic T^3-3*(x+y^2)*T-2*x has discriminant -108G0, but
has a triple root over (0,0). Its local reflections are distinct
and generate S3. The corresponding cyclic triple over the A5 sign
surface is only quasi-etale and is not etale at that singular point.
Indeed the cubic quotient has a single point in that fibre and is
smooth there because its x-derivative is -2. An etale sign triple
would instead give the two local cubic factors (9), including the
singular double factor. Thus immersion and A1 normalization alone do not exclude a
monogenic S3 cubic; the everywhere-etale sign condition is essential.
This example is not a Keller map.

Finally p=q=2 is outside the coprime hyperbolic theorem:
(x+1)^2-(x-1)^2=4x is an irreducible A1 branch with no common zero
of the two coefficients. Reducible G also falls outside the theorem,
as x^2-1^3=(x-1)(x+1) shows. These are explicit failures of the
dropped hypotheses, not searched candidate families.

## 5. Evidence and remaining scope

[The checker](verify_power_difference.py) verifies the universal
piecewise valuation identities symbolically, the cubic trace/index
identity and the negative controls. Exact primary bytes and inspected
complete pages are pinned in [the manifest](sources/power-difference/manifest.json).
Snyder's complete proof and Fujino's complete differential-pullback
lemma were inspected. The torus log form, retained exceptional curve,
both image-dimension cases, and actual cubic-index consequence are
proved above. No general Iitaka subadditivity conjecture is imported.

Catanese-Perroni Section5.3 was consulted as primary background on
dihedral quotients and trace-zero modules, not as a monogenicity
theorem. Its PDF31 displayed expansion of (c+tau(c))^3 has a printed
coefficient2 where commutative expansion requires3; no numerical
identity from that display is used here. The trace splitting and
index determinant above are derived directly.

The new theorem does not show that the actual cubic algebra must be
monogenic. That missing implication would be needed to turn (10)
into a full quartic exclusion. Couple its unavoidable index to the
marked source and [PF53's trace-dual algebra](TRACE_DUAL_PENCIL_AND_EXCESS.md),
while retaining [PF55's Picard/linking obstruction](SIGN_SURFACE_PICARD_AND_SQUARE_DETERMINANT.md).
The power-difference theorem is valid in all polynomial degrees and
all coprime exponents. Applying it in another Keller field degree
requires an independently justified appropriate cover; none is
manufactured here. Independent full reproof, formal verification,
external review and novelty priority remain unclaimed. The full
research goal remains active.
