# Adjoint-fibration exclusion of the two nonimmersed normal-r7 rows

14 November 2026. This is a written geometric proof with exact univariate
certificates. Its external Borisov inputs, including the source-model
finality assertion used in the proof of his Theorem 3.8, are explicit in
[SOURCE_AUDIT.md](SOURCE_AUDIT.md). No independent reproof of that paper,
formal verification, external review of this deduction or novelty priority
is claimed. Ordinary JC2 remains open.

## 1. Hypotheses and notation

Assume a complex plane Keller map has geometric degree six and its affine
finite normalization N has precisely one ramified prime E of generic
(e,f)=(2,1). Its image is one of the two nonimmersed curves in
PF2R7CUSP1, parametrized birationally by the retained monic polynomials
(p(t),q(t)) of degrees (6,10). All unramified omitted primes are retained.

The coefficient classification normalizes by affine parameter changes,
target translations/scalings and q->q+C*p+F. These are affine target
changes, so preserve the target line at infinity, the Keller property
and the geometric reduction. Once the retained coordinates are chosen,
the common target infinity point is [0:1:0] because q has degree10
and p degree6. No nonlinear target coordinate reduction is used.

Let W be the finite projective normalization over target P2, and
q:Y -> W its minimal smooth resolution. Write g:Y -> P2, H=g*L_infinity,
K=K_Y and R_tot=div(det(dg)). These are the actual models of JC2POS1--3.
The source open U=A2 is unchanged. Write C_i for the primes over generic
target infinity, and e_i,f_i for their normal and residue degrees.

We prove that these hypotheses are impossible for either coefficient
class. This does not classify immersed r7, normal r5 or nonnormal slices.

## 2. The adjoint pencil collapses the rank window to eight

JC2POS3 gives H^2=6, K.H=-4, gamma=2 and 8<=rho(Y)<=11.
Put D=K+H. Adjunction on a smooth general H-member and rationality give

    h0(D)=2, D.H=2, D^2=8-rho.

Remove the fixed divisor Z>=0 from |D|, leaving the complete movable
system |M|. Its nonzero class has M^2>=0 and H.M>0. Hodge index gives
M^2<= (H.M)^2/6 <=4/6, hence M^2=0. Two general movable members have
no common component and zero intersection, so |M| is basepoint-free.
After Stein factorization the base is P1 because Y is rational. Since
h0(M)=2, its general connected fibre F occurs with multiplicity one:
M=F, not a multiple of a fibre. Thus D=F+Z and H.F is either one or two.

H.F=1 is impossible. The restriction of g to a smooth general F would
be degree one onto a line, hence F=P1. Then (2K+3H).F=-4+3=-1,
contradicting JC2POS1. Consequently H.F=2 and H.Z=0. All components
of Z are g-contracted. On a relatively minimal smooth model adjunction
gives K.B>=0 for every contracted irreducible curve B: its square is
negative and the only negative-canonical case would be a smooth
rational (-1)-curve. Thus K.Z>=0. Since

    6-rho = K.D = K.F+K.Z >= -2,

we get rho<=8. Hence rho=8, K.F=-2 and K.Z=0. Also
F.Z=(K+H).F=0 and D^2=0 imply Z^2=0. Hodge index on H.Z=0 forces
Z=0. We have proved

    rho=8, K^2=2, |K+H| is a basepoint-free rational pencil,
    F=K+H, F^2=0, K.F=-2, H.F=2.                         (A)

The whole complete linear system and connected-fibre step are essential;
a formal pencil on a partial boundary would not establish (A).

## 3. Generic infinity and the cusp alternative

Borisov's Theorem 3.5(1) says f_i=1 for every type1 prime. Thus
JC2POS3's residue-degree sum two gives exactly two C_i. The former
one-prime (3,2) option is excluded. The remaining index pairs are
(1,5),(2,4),(3,3).

PF2R7BND1 gives t=3-epsilon additional unramified affine omitted primes,
where epsilon=0 or1 describes the cusp's local orbit size two or three.
At epsilon0 the normal degree-two local cover branched along an ordinary
cusp is analytically

    z^2=x^2-y^3, equivalently a*b=y^3.

This is a Du Val A2 singularity and its minimal resolution has two
(-2)-curves. The other finite normalization points are smooth: at every
two-branch value the separated ramified orbits give double covers of
one smooth target branch; all simple local factors are etale. At epsilon1
the cusp orbit is the smooth standard cubic cusp map
(s,z)->(s,z^3+s*z); the affine normalization is smooth everywhere.

Borisov's Theorem 3.4 requires W to be singular at its common infinity
point A, and all its other possible singular points are the affine cyclic
points just considered. A resolution needs at least one curve above A.
But at epsilon0 the whole-boundary count would give

    contracted curves = 8-(1+3+2)=2,

already exhausted by the affine A2 singularity. This is impossible.
Therefore epsilon1, t2, and exactly three contracted curves, all over A.
These counts are checks on the complete normalization, not on its
ramification divisor alone.

## 4. Transporting the type1 final-curve assertion

For each type1 curve, use the source model Z in Borisov's proof of
Theorem 3.8 on which that strict transform is a final smooth (-1)-curve.
The argument may be applied separately to each curve; a simultaneous
choice is not needed. The cited source asserts this choice; finality is NOT inferred from the
word 'minimal' for Y. The resolution Z -> W factors through Y, and
Z -> Y is a sequence of smooth blowdowns of g-contracted (-1)-curves.

For a type1 curve C on Z, H.C=f_i=1, K_Z.C=-1, so Q_Z.C=1 for
Q_Z=2K_Z+3H. At each such blowdown its intersection changes by

    Q_new.C_new = Q_old.C_old - 2*(C_old.E_exceptional).

The intersections subtracted are nonnegative integers. The final Q_Y is
nef by JC2POS1. Hence no positive intersection can have been subtracted:
otherwise 1 would become negative. No blowdown touches C. Thus on Y

    C_i is smooth, C_i^2=-1, H.C_i=1, (K+H).C_i=0.       (I)

The two infinity curves are vertical for the rational pencil (A). This
is the precise geometric bridge to the parabola calculation below.
Omitting the source-model assertion would leave additional horizontal
infinity cases, and the present conclusion could not be inferred.

## 5. General pencil members are parabolas with a fixed infinity point

Let B_j be the three contracted curves above A, let m_j>0 be their
coefficients in H=g*L_infinity, and put d_j=F.B_j>=0. The actual divisor
H consists of the two type1 primes and these contracted primes. By (I),

    2=H.F=sum_j m_j*d_j.                                 (B)

Thus at least one B_j is horizontal for the pencil. All its intersection
points with a general F map to the same target point a=g(A) at infinity.
The restricted map g_F has (degree of map)*(degree of image)=H.F=2.
There are two alternatives: a double cover of a line or an isomorphism
onto a smooth conic.

The line alternative is impossible. All image lines pass through a, so
they form a pencil of parallel affine lines. If eta:P1->P1 sends the
source pencil parameter to the image-line parameter, degree conservation
gives deg(eta)=3. The target pencil parameter is an affine linear
polynomial ell. On U,

    ell(F_source)=eta(pi).

Every pole of eta is omitted by pi(U). Choose one such omitted point as
infinity. Then pi|U is a nonconstant polynomial r on A2, hence surjective
onto A1. All poles of eta must consequently be at that one infinity
point, so eta is a cubic polynomial in r. Its derivative has a finite
root, where the differential of eta(r)=ell(F_source) vanishes. This
contradicts the Keller condition. Thus the general image is a smooth conic
and g_F is an isomorphism.

Since g_F is injective, only one point of a general F may map to a. In
(B) this forces exactly one horizontal B_j, with d_j=1 and m_j=2.
It is a section of pi and is smooth P1; adjunction gives B_j^2=-3.
The two other contracted curves are vertical, hence K.B_j=F.B_j=0
and are smooth (-2)-curves. Their incidence is not assumed to be SNC.

The general conic meets target infinity only at a, with multiplicity two.
The nonproper curve C also has its unique infinity point at a, by
Borisov's common-point theorem. In the retained coordinates this is
[0:1:0]. A smooth conic tangent to Z=0 at [0:1:0] has affine equation

    y=a(z)*x^2+b(z)*x+c(z),  a(z)!=0,                    (P)

where z is the rational pencil parameter and a,b,c are in C(z).
The section supplies a rational point of the generic conic, and x is
its affine coordinate. Hence C(Y)=C(z,x). The field degree of the
rational function in z on the right of (P), over C(x), is exactly six.

## 6. The ramified curve meets a general fibre once or twice

On the affine part of a general pencil fibre the only ramification
divisor is E, with coefficient one. In the coordinates (z,x), its
equation is

    a'(z)*x^2+b'(z)*x+c'(z)=0.

It is not identically zero (otherwise g is not dominant). All three
coefficients cannot have zero derivative. A vertical E would have
H.E<=H.F=2, contrary to H.E=deg(C)=10. Thus the map pi|E^nu has degree
one or two. Boundary coordinates and unramified omitted curves create
no additional critical divisor at a general finite x and generic z.

At a generic point t on E^nu, tangency in (P) gives

    a(z)*p(t)^2+b(z)*p(t)+c(z)=q(t),
    2*a(z)*p(t)+b(z)=s(t),  s(t)=q'(t)/p'(t).             (T)

All calculations use generic points where p'!=0; rational identities
extend to the whole function field.

After shrinking the pencil base, (z,x) are regular coordinates at every
finite-x point of a general fibre: a P1-bundle with its infinity section
removed has this affine coordinate on its generic fibre, and the generic
isomorphism spreads over an open base. Thus no coordinate Jacobian zero
has been counted as extra affine ramification. The two unramified omitted
curves remain regular finite-x points of this surface chart.

### 6.1 Degree one: the rational-function degree is too large

If pi|E^nu has degree one, t is a rational coordinate on the pencil
base. The derivative equation is linear in x: a second distinct root
would give another affine ramification point on the generic fibre, and
a repeated root would give coefficient two rather than E's coefficient
one. Thus a'(z)=0 and a is a
constant A (nonzero in the smooth-conic case). Formula (T) gives

    y=A*(x-p(t))^2+q(t)+s(t)*(x-p(t)).                   (L)

The two exact classes have gcd(p',q')=t-1. Put d=p'/(t-1), n=q'/(t-1).
Then deg(d)=4, deg(n)=8 and gcd(d,n)=1. The numerator of (L) is

    A*d*(x-p)^2+d*q+n*(x-p).

Over C(x), its gcd with d is one: modulo d the numerator is n*(x-p),
and neither n nor x-p vanishes identically at a root of d. If A!=0 the
numerator has degree16 in t. If A=0 it has degree14, with nonzero top
coefficient since 1-10/6=-2/3. Thus its rational-function degree is16
or14, contradicting six. The A0 boundary is included as a control.

### 6.2 Degree two: no possible involution

If pi|E^nu has degree two, its nontrivial deck automorphism is a Mobius
involution iota of P1_t. The paired points t and v=iota(t) must satisfy
the Hermite quadratic-contact identity

    (p(t)-p(v))*(s(t)+s(v))=2*(q(t)-q(v)).                (H)

Both p' polynomials are squarefree. Thus s has only simple finite poles.
If iota(infinity) were finite, s(iota(t)) would grow at most linearly
as t->infinity, while p(iota(t)) and q(iota(t)) stay bounded. In (H)
the coefficient of t^10 on the left minus right is 10/6-2=-1/3,
a contradiction. Hence iota fixes infinity, and has the form 2h-t.

Clear the slope denominators and form

    E_h(t)=(p(t)-p(2h-t))*(q'(t)*p'(2h-t)+q'(2h-t)*p'(t))
           -2*(q(t)-q(2h-t))*p'(t)*p'(2h-t).

For BOTH exact coefficient classes,

    coefficient(t^17,E_h)=3840*h^3.

Therefore an identity E_h=0 would force h=0. At h0 the leading
coefficient (at t^11) is -1280/243 in class I and -2048/6075 in class II.
It is nonzero in either case. The checker additionally replays exact
Bezout identities between these t17 and t11 coefficient polynomials,
using Fraction arithmetic independently of its polynomial gcd producer.
Thus no Mobius involution satisfies (H), closing the degree-two case.

## 7. Coverage and limits

The degree-two cusp orbit is excluded by the global rank/boundary count;
the degree-three cusp orbit is excluded by (A)--(P) and the exhaustive
degree-one/two ramification alternatives. Both retained coefficient
classes are therefore excluded under the explicit single-ramified-row
hypothesis and named external inputs.

No assumption that omitted primes equal ramified primes was introduced.
The excluded source has two unramified affine omitted primes in the
remaining cusp branch, until the final contradiction. The generic
infinity index pairs (1,5),(2,4),(3,3) are all covered by the same
verticality and field-degree argument; no graph search truncation occurs.

This does not prove all normal r7 is excluded: the immersed coefficient
family is separate. Nor does it exclude additional ramified affine rows,
normal r5, nonnormal terminal slices, the whole (75,125) candidate or JC2.
The geometric proof and cited source assertions are not certified by the
univariate checker. See SOURCE_AUDIT.md for the exact dependency boundary.
