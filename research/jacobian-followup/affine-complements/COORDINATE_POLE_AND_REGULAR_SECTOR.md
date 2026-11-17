# Coordinate-invariant poles and the simple regular-coordinate sector

14 November 2026. PF2RULED38 replaces a coordinate-reduction guess by
an invariant and a complete marked normal form in one whole sector.
There is no bound on the field degree, polynomial degrees or q.

If an actual source coordinate y is regular on the normalization and
a complementary coordinate has pole q, every negative source-coordinate
pole is a multiple of q. The only regular source coordinates are affine
functions of y. Thus the q>=2 controls in PF2RULED37 really have no hidden
multiplicity-one coordinate ruling.

If y has a simple zero on the sole omitted line and the source volume
has divisor R, the entire smooth affine surface has the two-plane
normal form

    y=u,  x=u*t+h(u),
    h(u)=sum_(j=1)^q c_j*u^(-j), c_q!=0.                (1)

Its volume is exact exactly when c_1=0. The coefficients modulo the
two scalar changes below classify the marked surface in this sector.
Finiteness over a Keller target remains a separate condition.

An actual finite Keller map would additionally give a monic norm
family of degree d. Every nonzero parameter member is an immersed,
self-colliding A1-normalized curve; the zero member is the reduced
union of the boundary curve and one such immersed curve. This is a
constraint on the full finite algebra, not a constructed map.

The proofs use pinned primary automorphism and affine-line-fibration
theorems, as well as the one-line injectivity theorem. They are written proofs
with exact supporting controls, not formal verification or external review.

## 1. A valuation lemma covering every source automorphism

Let v be a valuation on C(x,y), trivial on C*, with

    v(x)=-q<0,  v(y)>=0.                                (2)

For every polynomial coordinate f of C[x,y],

    v(f)>=0  iff  f=a*y+b,  a!=0;
    { -v(f) : f a coordinate, v(f)<0 } = q*Z_(>0).       (3)

Here a polynomial coordinate means a member of an actual polynomial
automorphism, not a polynomial with merely a Jacobian mate.

Use the Jung--van der Kulk normal form, recorded in Neumann--Wightwick
Section1. Bruhat decomposition of an affine factor outside the triangular
subgroup writes it as triangular-affine, swap, triangular-affine.
Absorbing these factors gives an arbitrary automorphism in the form

    A_final o B_k o swap o ... o B_1 o swap o B_0,

where B_0 is arbitrary triangular, every B_j with j>=1 is triangular
of nonlinear degree n_j>=2, and A_final is affine. This is just the
alternating normal form; its length is arbitrary and finite.

Initially B_0 gives coordinates ax+p(y), by+c. Their pole orders are
q and nonpositive, respectively. The negative value cannot cancel
against p(y), whose valuation is nonnegative.

After the first swap and nonlinear triangular factor, the pole orders
are n_1*q and q, with the first strictly larger. At every further
step the nonlinear power of the larger-pole coordinate dominates
every other summand strictly. The next pair has pole orders n_j*M,M
if the preceding larger order was M. Thus both are positive multiples
of q and distinct. No coefficient or special residue choice permits
a cancellation at a strict inequality.

A final affine linear combination either keeps the larger pole,
keeps the smaller one, or, if no nonlinear step occurred and its
coefficient of ax+p(y) is zero, is affine in y. This proves necessity
in (3) and the claimed uniqueness of regular coordinates. Conversely
y+x^n is a coordinate of pole nq for each n>=1.

For a marked normalization (N,U,R), define

    mu(N,U,R)=min {-v_R(f): f a source coordinate, v_R(f)<0}. (4)

It is invariant under every isomorphism preserving the actual source
open and omitted valuation. Under (2), it equals q; changing the
complementary coordinate cannot lower it. This is a minimum over all
source automorphisms, not a bounded word search.

Primary import: [Neumann--Wightwick, Section1](https://arxiv.org/pdf/math/9911104v1).
Only its automorphism normal form is imported here, not its algorithms
or degree bounds.

## 2. Classify the sector v_R(y)=1 globally

Let N be a smooth affine surface with U=A2_(x,y), sole complement
R=A1, and

    div_N(dx wedge dy)=R,  v_R(y)=1.                    (5)

A shift has set the boundary value of y to zero. The normal-domain
valuation criterion makes y regular on N. Its divisor is

    div_N(y)=L+R,  L={y=0} in U.                        (6)

The retained line L is closed in N and disjoint from R. More generally
an embedded A1 in an affine open cannot acquire an endpoint inside
an affine surface: the normalization of its closure would be an
affine normal curve containing A1 and an additional point. The only
possible added point in its smooth completion is infinity; adding it
would make P1, contradicting affineness of that normalization.

Put W=N minus L. It is affine, since L is an effective Cartier divisor
on the smooth affine N. The morphism y:W->A1 is faithfully flat with
generic fiber A1; its fiber over zero is the reduced R and all other
fibers are source affine lines. It is smooth as well: y has a simple
zero along R, which is disjoint from L.

Miyanishi, TIFR60 ChapterI Theorem4.1.1, makes this an affine-line
bundle. Its base is A1, so both its linear bundle and its torsor class
are trivial. Therefore W=A2_(u,t) over u=y. All its hypotheses,
including affineness, reduced integral fibers and factorial base,
have been checked. [Primary lecture notes](https://mathweb.tifr.res.in/Documents/Publications/Lectures/tifr60.pdf).

On U intersect W={y!=0} the transition is

    x=a*u^k*t+h(u),  a!=0, k integer, h Laurent.

The volume divisor in (5) forces k=1. Scale t to set a=1. Removing
the nonnegative part of h by a source shear x->x-p(y) gives (1).
A negative term must remain: otherwise x and y would both be regular
on N, and A=C[x,y] intersect O_(N,eta_R) would equal C[x,y], making
the given nonempty complement impossible.

Thus q>=1 and c_q!=0. The entire U and W are plane charts, with
exactly their displayed overlap. No formal neighborhood has been
mistaken for an entire chart. For an actual one-index-two-prime
Keller normalization PF2OMIT7 supplies smoothness, and the volume
equation gives (5). The extra condition v_R(y)=1 is explicit;
a regular coordinate with higher vanishing order is not covered
by this classification.

## 3. An explicit affine ring and marked-isomorphism criterion

Conversely every principal part h in (1) gives a smooth affine
surface with the required marked plane, omitted line and volume
divisor. Put H(y)=y^q*h(y), so H(0)=c_q!=0, and set

    f=x*y^q-H(y),   a=x*f,   b=x^2*y^(q-1)*f,
    A_h=C[y,f,a,b] inside C[x,y].                       (7)

On W these are

    y=u, f=u^(q+1)*t,
    a=u*t*(H(u)+u^(q+1)*t),
    b=t*(H(u)+u^(q+1)*t)^2.                            (8)

There are three distinguished affine opens:

* at f, x=a/f recovers the source chart D_U(f);
* at f+H, t=b/(f+H)^2 and u=y recover D_W(f+H);
* at y, x=(f+H)/y^q and t=f/y^(q+1) recover the overlap.

Writing H(y)=c_q+y*H_1(y), the identity

    (f+H)-f-y*H_1=c_q

shows that these opens cover Spec(A_h). This proves affineness and
the full ring identification, uniformly in q and all c_j.
The special controls of PF2RULED37 are H=1.

By Section1 the regular source coordinates are exactly ay+b.
Two normalized models A_h and A_h' are isomorphic preserving
the plane open and omitted line if and only if they have the same
q and there are alpha,b in C* such that

    c_j'=b*alpha^j*c_j for every j.                     (9)

Indeed any marked isomorphism restricts to a source automorphism,
and regular-coordinate uniqueness forces y'=alpha*y (both boundary
values are zero). A complementary coordinate is x'=b*x+p(y).
For this to extend over R, its expression on W must satisfy

    u'=alpha*u,
    t'=(b/alpha)*t+p(u)/(alpha*u),
    h'(alpha*u)=b*h(u).

The negative coefficients give (9), and the constant term of p
must vanish. These conditions are sufficient as well: the displayed
polynomial affine transformation of W has a polynomial inverse.
Thus no additional marked identifications are being presumed.
This is not a classification after forgetting the source open.

In particular q>=2 models have mu=q and have no source coordinate
of pole one under any automorphism. A multiplicity-one source-coordinate
ruling would have precisely that pole. They cannot be identified
with the marked Wright model, whose source coordinate X has pole one.

## 4. Exact volume and a constructive primitive for every principal part

The two local primitives x dy and u*t du of dx wedge dy differ by
h(u) du. Its residue is c_1. The two-plane argument in PF2PERIOD1
therefore gives exactness if and only if

    c_1=0.                                             (10)

Here sufficiency can be constructed rather than inferred. For c_1=0 put

    K(u)=sum_(j=2)^q c_j/(1-j)*u^(1-j).

Choose a polynomial p(u) in u*C[u], of degree at most q-1, whose
negative Laurent part in p(u)*h(u) is K(u). Multiplication by h,
followed by negative-part projection on this finite space, is
triangular with diagonal c_q. It is invertible since c_q!=0;
solve starting with the u^(-(q-1)) coefficient. This is a uniform
triangular construction, not a coefficient search for a Keller map.

Then

    K_U=x*p(y),
    K_W=u*t*p(u)+(p(u)*h(u)-K(u))

are polynomial on their respective charts and K_U-K_W=K(u) on
the overlap. Hence the forms

    theta_U=x dy-dK_U,  theta_W=u*t du-dK_W              (11)

glue to a global regular primitive of dx wedge dy.

If a finite Keller map exists in this sector, its volume is exact,
so c_1=0, q>=2 and the marked invariant mu is at least two.
For every such q the models really exist as affine surfaces with
exact volume, including the H=1 controls. They still have no
alternative multiplicity-one coordinate ruling. This closes the
coordinate ambiguity explicitly left open in PF2RULED37; it does
not produce a finite Keller map or exclude one in this sector.

## 5. Necessary norm family for the actual finite algebra

Now suppose pi=(P,Q):N->A2 is an actual finite extension of a
polynomial Keller map, of field degree d, with the one-index-two
boundary and the simple regular coordinate above. Put B=C[P,Q],
K=Frac(B), and L=C(N).

For every beta in C the source line L_beta={y=beta} is closed in N.
Its finite image D_beta has normalization A1 and its normalization
map is immersive. To check degree one, factor its finite
parametrization through the normal curve D_beta^nu=A1. The
polynomial factor A1->A1 has no finite critical point, since
F is etale on L_beta. It is therefore affine linear.

Let H_beta be the reduced irreducible equation of D_beta and H_C
that of the omitted boundary image C. The latter has a finite
nonimmersion by PF2RULED37's all-degree Euler/Chau argument;
hence it is distinct from every immersed D_beta.

Because y is integral over B, all its field norms are polynomials.
Push forward the principal divisors along pi. Equations (6) and
div_N(y-beta)=L_beta for beta!=0 give, up to nonzero constants,

    Norm_(L/K)(y)=H_C*H_0,
    Norm_(L/K)(y-beta)=H_beta for beta!=0.               (12)

The coefficient of each image divisor is one: both finite curve
maps have residue degree one. In particular y is a primitive
element of L/K. Otherwise every displayed norm would be a proper
[L:K(y)]-th power by the tower formula, contradicting the reduced
irreducible divisor for beta!=0.

Its minimal polynomial chi(T)=Norm_(L/K)(T-y) is therefore monic
of degree d, with coefficients in B, and (up to constants)

    chi(0)=H_C*H_0,
    chi(beta)=H_beta for every beta!=0.                 (13)

Every nonzero parameter member is a reduced irreducible curve
with immersive A1 normalization; the zero member is a reduced
union of two distinct curves, one of them nonimmersed. This is
one monic polynomial family of plane curves. It is not a linear
pencil or a polynomial fibration of the target plane.

There is also a compulsory collision in every D_beta.
Gwozdziewicz, *Injectivity on one line* (1993), Theorem1.1 says that
a polynomial Keller map injective on one affine line is an automorphism.
Apply it in the actual source coordinates (x,y). It makes every
restriction to L_beta noninjective, since the actual map has a
nonempty omitted boundary. This imported theorem already includes
its affine-line embedding input; it is not inferred from a finite test.
[Primary two-page paper](https://arxiv.org/pdf/alg-geom/9305008v1).
Thus every member D_beta is self-colliding, including H_0.

One algebraic consequence is that S=B[y] is not normal: over a
self-collision on L_beta its finite birational normalization N->Spec S
identifies two distinct source points, since their y values agree.
Such collisions exist for every beta. The conductor support V(Ann_S(A/S)) has a
component dominating the y-line, since a finite-type support with
only vertical components could meet only finitely many fibers.
This is an additional global finite-algebra requirement. A primitive
integral source coordinate must not be confused with generation of
the normal ring by that coordinate.

## Scope and next step

The pole-spectrum theorem covers every valuation with one regular
source coordinate. The normal form, affine ring and marked-isomorphism
classification cover the full v_R(y)=1, canonical-divisor-R sector.
The norm-family and conductor constraints additionally require an
actual finite Keller map. All statements are uniform in field degree.

Regular coordinates of higher vanishing order, normalizations with
no regular source coordinate, extra omitted primes, higher-index
packets and the full finite-map problem remain open. In the classified
sector the remaining question is finiteness and the Keller algebra on
the genuine q>=2, c_1=0 marked models, not a further coordinate reduction.
Keep the original goal active. No formal verification, independent
full reproof, external review or novelty priority is claimed.

Primary inputs are pinned in
[sources/coordinate-pole-sector/manifest.json](sources/coordinate-pole-sector/manifest.json).
