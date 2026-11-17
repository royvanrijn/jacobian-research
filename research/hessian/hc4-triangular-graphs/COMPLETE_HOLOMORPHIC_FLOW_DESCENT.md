# Complete holomorphic flow descent and affine-line-normalized curves

Date: 2026-11-16. Claims: HC4MYKZ0 and HC4MYKZ1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: the actual polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, equivalently regular inverse Gram entries;
the alternative full marked-energy integrality implication is also OPEN.
Neither implication is proved here. External review/outreach remains DEFERRED.

The preceding [KY0/1 theorem](TARGET_BOUNDARY_SYMMETRY_AND_FINITE_STABILIZER.md)
excludes locally finite logarithmic target fields. We extend that exclusion
to EVERY polynomial target vector field with a holomorphic flow defined for
ALL complex times, even when its time maps are transcendental. We also
extend [KX0](ALGEBRAIC_FLOW_DESCENT_AND_INTEGRAL_COORDINATES.md) from locally
finite source fields to all holomorphically complete polynomial source
fields preserving the actual target field. No such field is shown to exist
for a hypothetical counterexample.

## 1. Exact classification input and the normalization observation

We adopt ONE NEW PROVED EXTERNAL INPUT: Brunella's classification of complete
polynomial vector fields on C^2, up to POLYNOMIAL coordinate changes,
Topology 43(2)(2004),433-445, DOI10.1016/S0040-9383(03)00050-8.
The original publisher page/abstract was checked, but its full PDF was not
accessible. The exact three normal forms below were read in the primary
research paper by Bustinduy, arXiv:0902.4726v2, Section1.2, printed/PDFp3,
with its Suzuki definitions on printed/PDFp2. That paper proves an extension
of the classification and explicitly credits these polynomial normal forms
to Brunella. The classification is ADOPTED AS REPRODUCED THERE; neither its
proof nor the inaccessible original full text is claimed independently
audited. The cached primary PDF and this precise access boundary are pinned.

All polynomials and constants below are complex. Write D for the field.
The forms, without any degree bound, are:

I. D=(c*x+d)*partial_x+(A(x)*y+B(x))*partial_y.

II. Put h=x^m*y^n, with m,n positive coprime integers. Then

    D(x)=n*lambda(h)*x,  D(y)=(a-m*lambda(h))*y,
    D(h)=n*a*h.                                           (1)

III. Put w=x^ell*y+p(x), h=x^m*w^n, with m,n,ell positive,
gcd(m,n)=1, deg(p)<ell, p(0)!=0. In the rational coordinates (x,w),

    D(x)=n*lambda(h)*x,  D(w)=(a-m*lambda(h))*w,
    D(h)=n*a*h,                                           (2)
    D(y)=[a*w-lambda(h)*((m+n*ell)*x^ell*y
                         +m*p+n*x*p')]/x^ell.

Here lambda is polynomial. The condition making III polynomial is

    lambda(h)*(m*p+n*x*p')-a*p in x^ell*C[x,y].             (3)

In particular a=m*lambda(0), by evaluation at x=0. The map (x,y)->(x,w)
is only a BIRATIONAL change, invertible on x!=0. It is never substituted
for an ambient polynomial automorphism at x=0.

Let C be an irreducible invariant algebraic curve whose normalization is
A^1. The complete holomorphic flow preserves C, and lifts to its
normalization, including all points above singularities. Indeed every time
map is an analytic isomorphism and hence lifts uniquely through the finite
normalization; the lifts vary holomorphically and satisfy the group law.
Every biholomorphism of the complex line is affine. Thus, in a normalization parameter t,
the restricted infinitesimal field is

    (alpha*t+beta)*partial_t.                             (4)

This includes alpha=beta=0 when D vanishes on C. The two coordinate
functions x(t),y(t) are POLYNOMIAL, since the normalization is A^1 and the
normalization morphism is algebraic. No smoothness, immersion, or injectivity
of that morphism is assumed.

For an invariant REDUCED finite union, each component is invariant: the
connected time group cannot permute finitely many components. The following
argument therefore applies to every component, in the SAME ambient normal
form. Independent coordinate changes for different components are not used.

## 2. KZ0: a complete field produces a locally finite symmetry of the curve

THEOREM. Let C be a nonempty reduced plane curve with A^1 normalization
for EVERY irreducible component. If a nonzero holomorphically complete
polynomial vector field is logarithmic for C, then C admits a NONZERO
LOCALLY FINITE polynomial logarithmic derivation. The latter field need
not equal the original one. The proof covers fields vanishing identically
on some or all components.

If the classified field is already locally finite, there is nothing to
prove. We treat all remaining cases.

### Form I

If c=d=0, vertical lines are invariant. Any nonvertical invariant C has
x(t) nonconstant, so D(x)=0 and (4) force the induced field to be zero.
Hence C is contained in A(x)*y+B(x)=0. If A=0 and B!=0, there is no such
curve and the original field is locally nilpotent. Otherwise divide A,B
by their gcd, obtaining coprime A0,B0. The unique possible nonvertical
curve has coordinate ring

    C[x,y]/(A0*y+B0)=C[x,1/A0].                            (5)

Bezout gives 1/A0 as a polynomial in x,y. Its normalization can be A^1
only if A0 is constant, because A^1 has no nonconstant units. Thus the
curve, if present, is a polynomial graph y=q(x). Scaling y-q(x), while
fixing x, preserves that graph AND every vertical line. If only vertical
lines occur, translation in y preserves their whole union. These supply
the required nonzero locally finite fields.

If c=0,d!=0, rescale so D(x)=1. Equation (4) applied to x(t) says
(alpha*t+beta)*x'(t)=1, whence alpha=0, beta!=0, and x is affine linear.
Every invariant component is a polynomial graph y=q(x). Two distinct
graphs would give a nonzero polynomial f satisfying f'=A*f. Comparing
degrees forces A=0 and f constant; in that case D=partial_x+B*partial_y
is locally nilpotent. Otherwise there is at most ONE graph and scaling
about it supplies a locally finite symmetry. No vertical curve is invariant.

If c!=0, translate x so D(x)=c*x. When A is constant, D is locally finite:
the span of y,1,x,...,x^deg(B) is stable. Assume A NONCONSTANT. An invariant
vertical component must be x=0. For a nonvertical component, (4) gives

    (alpha*t+beta)*x'(t)=c*x(t).

Polynomial roots and degrees imply alpha!=0 and, after translating t,
x(t)=C*t^k, c=k*alpha, C!=0, k>=1. Decompose

    y(t)=sum_(r=0)^(k-1) t^r*Q_r(t^k).

For r!=0 the y equation has no B contribution, so, with z=t^k,

    alpha*(r*Q_r+k*z*Q_r')=A(C*z)*Q_r.                    (6)

A nonzero Q_r is impossible: the right side has degree deg(Q_r)+deg(A)
and the left side at most deg(Q_r). Thus y belongs to C[t^k]. Birationality
of the normalization map then forces k=1: its coordinate field is C(t),
not the proper subfield C(t^k). The component is a graph y=q(x). Two
different graphs would give c*x*f'=A*f, again impossible for nonconstant
A by degrees. Therefore the ENTIRE curve is a subset of one graph and
the vertical line x=0. Scaling y-q(x) preserves both; if the graph is
absent, translation in y suffices. This completes I in all degeneracies.

### A common calculation for forms II and III

Use the coordinate pair (x,z), where z=y in II and z=w in III, so
h=x^m*z^n and (1)/(2) hold. If a!=0 and h(t) is nonconstant on an
A^1-normalized component, then (4) and D(h)=n*a*h imply, after translating
t, h(t)=C*t^k and the induced field is alpha*t*partial_t with k*alpha=n*a.
The factorization x(t)^m*z(t)^n=C*t^k in C[t] forces

    x(t)=c1*t^r, z(t)=c2*t^s, r,s>=0, m*r+n*s=k.

The x equation now reads r*alpha=n*lambda(C*t^k). This is impossible
when lambda is nonconstant, even if r=0. If h is a nonzero constant,
D(h)=n*a*h is already a contradiction. Consequently, when a!=0 and
lambda is nonconstant, every allowed component lies in h=0.

If a=0 and lambda!=0, h must be constant on every invariant component:
otherwise (4), D(h)=0, and nonconstant h force the induced field zero,
contradicting D(x)=n*lambda(h)*x since x is not identically zero when h is
nonconstant. A nonzero constant level h=C has no zeros of x or z on its
normalization. Both would be units in C[t], hence constant, making the
whole normalization map constant. Thus such a level CANNOT have A^1
normalization. This handles zero restrictions as well, including every
level where lambda(C)=0 and the ambient field vanishes identically.

### Form II

Constant lambda makes D a linear, hence locally finite, field. For
nonconstant lambda, the preceding calculation leaves only h=0: the two
axes. Every reduced subset of their union has a nonzero Euler symmetry.

### Form III

If lambda is constant lambda0, (3) first gives a=m*lambda0 and then
n*lambda0*x*p' divisible by x^ell. Since deg(x*p')<ell, either lambda0=0,
in which case D=0, or p is constant. In the latter case
D=n*lambda0*(x*partial_x-ell*y*partial_y), locally finite.

For nonconstant lambda, the common calculation permits only h=0. Its
components are x=0 and w=0. On w=0, x is a nonconstant invertible regular
function: p(0)!=0 and Bezout show its coordinate ring is C[x,x^-1]. Thus
w=0 has normalization Gm, not A^1, and is forbidden by the hypothesis.
Only the line x=0 remains, with a translation symmetry. There is no
forgotten component at the birational chart boundary.

The normal forms exhaust all complete polynomial fields by the explicit
classification input, and all changes used to transfer the final symmetry
back to the original plane are POLYNOMIAL. This proves KZ0 without a bound
on degrees, component counts, singularities, or normalization collisions.

## 3. KZ1: actual plane Keller consequences and complete source descent

Let K:X=A^2->Y=A^2 be complex Keller. Put R=C[X], B=C[K],
L=Frac(R), F=Frac(B), and let N be its entire reduced nonproperness curve.
Under noninvertibility N is nonempty and every component has A^1
normalization, by the retained polynomial-parametrization result PF2OMIT3
(also used in IY2). These are ACTUAL target components, not arbitrary
curves in an abstract normalization model.

If a nonzero holomorphically complete POLYNOMIAL target field preserves N,
KZ0 gives a nonzero locally finite logarithmic field for its WHOLE union.
KY1 excludes that. Therefore

    K noninvertible ==> N admits NO nonzero holomorphically
                        complete polynomial logarithmic field.    (7)

Now suppose delta is a nonzero holomorphically complete polynomial SOURCE
field and delta(F) subset F. The retained intersection identity F intersect
R=B (PF52) shows that delta induces a polynomial target field epsilon.
It is nonzero: a derivation killing F kills its finite separable extension L.
We prove that epsilon is complete; no polynomiality of its flow is assumed.

Let Phi:C times X->X be the complete holomorphic source flow. The actual
image V=K(X) is open and its complement has complex codimension at least2:
the retained actual-source sheet theorem HW0 supplies a source sheet over
EVERY target divisor. In the plane the complement is a finite set.

For y in V, choose x with K(x)=y and define

    Psi(t,y)=K(Phi(t,x)).                                 (8)

For two choices of x these are entire solutions to the SAME target ODE
with the same initial value. They agree near zero by uniqueness, then for
ALL t by the identity theorem. Thus (8) is well-defined. Local holomorphic
sections of the etale map K make it JOINTLY holomorphic on C times V.
Its values lie in V, and the source group law proves the target group law
and the inverse Psi(-t) on V.

Here the needed JOINT Hartogs extension has a direct Cauchy proof. Around
one missing point choose coordinates (z1,z2) with that point at zero and
no other missing point in the polydisc. For a coordinate function f of
Psi, set, for |z1|<r,

    ftilde(t,z1,z2)=(1/(2*pi*i))*integral_(|zeta|=r)
                         f(t,zeta,z2)/(zeta-z1) d(zeta).    (10)

The integrand is defined even at z2=0 since zeta!=0. The integral is
jointly holomorphic in t,z1,z2, locally on any time disc. For z2!=0 it
equals f by the one-variable Cauchy formula; continuity gives equality
also at z2=0,z1!=0. These local extensions agree uniquely and glue at all
finitely many missing points and all times. Thus each coordinate of Psi
extends to ALL of C times Y. The identity and group laws extend by density
on C times Y and C^2 times Y. Hence every Psi(t) is a global holomorphic
automorphism of Y and epsilon is complete. This is ANALYTIC extension,
not algebraic GAGA or local finiteness. Only the elementary one-variable
Cauchy formula is used; a separate general Hartogs theorem is unnecessary.

The equivariance K*Phi(t)=Psi(t)*K holds globally. Since the time maps
are homeomorphisms with continuous inverses, they are proper topologically.
They carry a sequence escaping to infinity with K-limit y to an escaping
sequence with K-limit Psi(t,y); inverses give equality. Thus Psi(t) preserves
the ENTIRE N, including components caused by unramified omissions. Formula
(7) gives the contradiction. We have proved in EVERY plane field degree

    nonzero complete polynomial source delta preserving F
                       ==> K has a polynomial inverse.           (9)

The stated classification and INVERSION CONSEQUENCE are PLANE-specific.
The known equivariant 3D Keller map in KW1 still has a
nonzero complete weighted Euler source field preserving B. It forbids
transferring (7)/(9) directly to a four-dimensional HC4 gradient. Use the
genuine plane base through retained Tang instead.

## 4. Hostile controls and what remains unsolved

Completeness is STRICTLY broader than local finiteness: x*y*partial_y has
entire flow (x,y)->(x,exp(t*x)*y), while delta^n(y)=x^n*y has unbounded
degree. It preserves every union of vertical lines and the graph y=0.
KZ0 constructs a DIFFERENT locally finite symmetry of that whole curve.

The A^1-normalization hypothesis cannot be replaced by rationality. FormIII
with m=n=ell=1,p=1,a=0,lambda(h)=h has

    h=x*(x*y+1),
    D=x*h*partial_x-(x*y+1)*(2*x*y+1)*partial_y.

It is complete by the classified form and D(h)=0. Every nonzero level
h=C has normalization Gm. Its x=0 boundary line has D=-partial_y.
The proof does not discard those Gm levels as if they were A^1 and does
not assume the birational (x,w) chart is an ambient automorphism.

A logarithmic polynomial field need NOT be complete. For the rational
nodal curve H=v^2-u^2*(u+1), use normalization
u=t^2-1, v=t*(t^2-1). The tangent Hamiltonian
D=H_v*partial_u-H_u*partial_v induces (t^2-1)*partial_t, which is NOT
complete: its solutions have finite complex poles. The two distinct
normalization points over the node cannot be moved by a complete affine
flow fixing their images. Merely writing down this Hamiltonian does not
supply the missing complete field on a possible N.

The proposed Lin-Zaidenberg next step was audited rather than promoted to
a new theorem: IY1 already imports its irreducible contractible-curve
classification, and the retained Chau math/0408048 primary paper Corollary3
already excludes a simply connected WHOLE N. A search hit claiming more
from Jelonek arXiv:2011.03472 is WITHDRAWN in v6 and is not imported here;
the repository already retains that version audit. No simple-connectedness,
one-branch, or symmetry-existence hypothesis is inferred from just one
point at infinity or a rational parametrization.

The full conjecture is still OPEN. No complete polynomial logarithmic
target field, or complete polynomial source field preserving F, is proved
to exist automatically. Noncomplete fields remain possible. The actual
trace-dual regularity and marked-energy integrality implications above
remain the last universal gaps. No whole surviving higher-degree packet
is excluded by this conditional obstruction alone.

The checker verifies exact algebraic transcriptions, degeneration controls,
source pins and scope boundaries. It is NOT a proof of Brunella's
classification, the Cauchy extension argument, global flow descent, or HC4. Independent
replay, formal verification and external review remain false.
