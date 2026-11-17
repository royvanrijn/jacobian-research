# Primitive quotients, the plinth obstruction and a global action certificate

15 November 2026. Claims **HC4MYJT0/1/2**. The K2 action obstruction
now concerns entire embedded invariant fields, not only displayed
derivations. Changing the speed of an action or the generators of
its invariant algebra cannot evade it. A different prospective
action must have a genuinely different rational quotient, even
after arbitrary K2-preserving source conjugacy.

There is also a coordinate-independent test: the action forced by
a hypothetical nonautomorphic unit-Hessian source member has a
proper nonprincipal plinth ideal, with no height-one prime over it.
The three-variable principal-plinth theorems must not be applied
to this four-variable action.

Finally the action approach has an exact global certificate. In
addition to an action, it requires a full polynomial coordinate
system with its inverse, a constant-Jacobian coefficient pair,
and equality of its kernel with the polynomial algebra generated
by the two base coordinates and the value parameter. These data
construct the actual unit-Hessian source function. Mere generic
freeness, a polynomial kernel or a matching first jet does not.

None of the needed nonstandard K2 actions or coordinate certificates
is constructed here. The normalized full-source questions, both
old endpoints and unrestricted HC4/JC2 remain open.

## 1. Statements

Let R=C[X1,...,Xn]. A nonzero polynomial derivation D is called
primitive here if the gcd of D(X1),...,D(Xn) is a unit. This is
equivalent to its image ideal having no common prime divisor.
It is preserved by polynomial source changes.

**JT0.** Let D be a primitive locally nilpotent derivation of R.
Every nonzero locally nilpotent derivation E with the same embedded
invariant field in Frac(R) has the form

    E=fD, f in ker D, f!=0.                            (1)

If E is primitive too, f is a nonzero constant. Equality of the
embedded invariant fields is equivalent to equality of the kernel
subalgebras in R. Abstract isomorphism of rational fields is not
this hypothesis. This is a reproduced elementary general fact,
not a novelty claim about locally nilpotent derivations.

For A=ker D define the plinth ideal

    pl(D)=A intersect D(R).

The action required by [JS0](GLOBAL_ADDITIVE_ACTIONS_AND_THE_K2_TRANSLATION_BARRIER.md)
for an N>=2 unit-Hessian source orbit is primitive. Its plinth
ideal is proper, nonprincipal and contained in no height-one
prime of A. No exact formula for that ideal is asserted.

**JT1.** For K2=s and its standard commuting actions D_u,D_v
from JS, take any nonzero module element

    D=a(x,s)D_u+b(x,s)D_v.

Divide a,b by their gcd in C[x,s]. The resulting module derivation
D0 is primitive in the entire four-variable source ring. Thus
every locally nilpotent derivation sharing D's embedded invariant
field is a polynomial invariant multiple of D0, and still has
fixed points on every s=kappa, kappa!=+-1.

The same holds for the y-standard module and for every embedded
field conjugate by an actual K2-preserving source automorphism.
No such invariant field can carry the action required by JS0.
This includes all invariant multiples, not just coefficients in
C[x,s], and all possible derivations with that same kernel.
Other rational quotient fields are not classified.

**JT2.** Let h in four variables have N>=2 finite Morse critical
points. Its full polynomial source orbit contains a unit-Hessian
member if and only if there are the following global data:

1. A full polynomial coordinate system (b1,b2,P,Q) of R, with
   an actual polynomial inverse.
2. Polynomials F,G in B=C[b1,b2] with J_(b1,b2)(F,G) in C*.
3. The derivation D=G partial_P-F partial_Q, with the equality
   of actual subalgebras ker D=B[h].

Under these hypotheses h=D0(b)+c(PF+QG), c in C*, in those
coordinates. Its Hessian determinant there is c^4 J(F,G)^2.
Its gradient is invertible exactly when the plane pair (F,G)
is invertible. In particular this is an exact global source
certificate, not a new sufficient test from an unqualified action.

Retain [HCTANG2](../tang-hc4-audit/STRUCTURAL_REVIEW.md),
[HG0](GLOBAL_VALUE_ORBITS.md), the JS uniform K2 value presentation,
complete invariant calculations and whole-fibre fixed-point theorem,
and the actual source-coordinate boundaries in
[JQ](TWO_NORMAL_MODELS_AND_RELATIVE_RANK_ONE_ALIGNMENT.md) and
[JR](ONE_HYPERBOLIC_COORDINATE_AND_INTRINSIC_CYLINDER_PLANE.md).
All new kernel, denominator, primitivity and plinth arguments
are proved below. No action classification is imported.

## 2. A primitive action is fixed by its embedded invariant field

Put A=ker D and K=Frac(A). Choose q with Dq=a!=0 in A:
take the last nonzero derivative before the kernel of any
noninvariant polynomial. In R_a put t=q/a, so Dt=1.
The finite Taylor operator

    pi(r)=sum_{j>=0} (-t)^j D^j(r)/j!

is D-invariant; successive Taylor subtraction expresses every
r in R_a as a polynomial in t with coefficients in A_a.
Its finiteness follows from local nilpotence. Thus
R_a=A_a[t], and after inverting all nonzero elements of A,

    R_(A minus 0)=K[t].

The constants in Frac(R)=K(t) are exactly K. Hence
A=R intersect K. Derivations with the same embedded invariant
field therefore have the same actual kernel A.

Such an E extends to a K-linear locally nilpotent derivation of
K[t]. It must be c partial_t with c in K*: if Et=p(t) had
positive degree, its iterated derivatives on t would never
vanish in characteristic zero. Since D=partial_t, E=cD.

Write c=u/v in Frac(R), with u,v coprime in the UFD R.
Polynomiality of E(Xi) makes v divide every D(Xi). Primitivity
makes v a unit, so c is in R intersect K=A. This proves (1).
If E is primitive, the common factor c must itself be a unit
of R, so c is a nonzero constant.

The primitive hypothesis cannot be dropped. For example
x partial_y and partial_y have the same invariant field, but
the latter removes the former's fixed hypersurface by dividing
by x. JT1 proves primitivity globally before using (1); it
does not silently assume fixed loci are invariant under every
rational change of an action's speed.

## 3. The necessary action has a nonprincipal plinth ideal

In the actual Tang coordinates the JS0 derivation is

    D=G(b) partial_P-F(b) partial_Q,
    A=C[b1,b2,h], J(F,G) in C*.

F and G are coprime in C[b1,b2], hence also in R. Otherwise
their common curve would make the Keller Jacobian vanish.
Thus D is primitive. Actual source conjugation preserves this.

The set pl(D) is an ideal of A: for a in A and r in R,
aD(r)=D(ar). It contains both F and G, since D(P)=G,D(Q)=-F.
It is proper: D vanishes at every point of each nonempty fixed
plane, so no D(r) can equal one. These fixed planes exist because
N>=2. Finally no height-one prime of the polynomial ring A can
contain both coprime F,G. Therefore no such prime contains pl(D).
In particular a principal proper plinth ideal is impossible.

This proves a global obstruction for every proposed action,
not merely the standard K2 modules. It does not identify the
plinth ideal with (F,G), assert purity of its support or classify
all actions having a nonprincipal plinth ideal. The already
retained native C2 example in JS has the same necessary action
data but nonunit native Hessian, so these tests are not sufficient.

## 4. Equidimensional K2 fibres control all hidden coefficient divisors

Use the exact JS notation s=K2,u=w,v=1-2xy-z^2 and

    2xy+(s+uv)^2+v-1=0, z=s+uv.

The polynomial map pi=(x,s):A4->A2 has every fibre of pure
dimension two. If x=c!=0, put r=2cy+z^2-1; the actual trace
coordinate is s=wr+z. The maps

    z=kappa-wr, y=(r-z^2+1)/(2c)

give a polynomial isomorphism of each pi fibre with A2_(r,w).
If x=0 the fibre is A1_y times the plane curve
w(z^2-1)+z=kappa. Every component of that nonzero plane
hypersurface has dimension one, including kappa=+-1.

Now let a,b in C[x,s] be coprime. Their common zero set in
A2 is finite or empty. Its inverse image under pi has dimension
at most two. Hence a,b cannot have a common prime divisor in
the four-variable source ring: a source prime hypersurface has
dimension three. This proves coprimality after pullback, rather
than assuming it from the abstract coefficient ring.

For D0=aD_u+bD_v the actual polynomial derivatives satisfy

    D0(u)=2xa, D0(v)=2xb.

Any common divisor of all D0(Xi) divides both displayed
derivatives, hence divides x. Since a,b are coprime they cannot
both be divisible by x. Let a0(s)=a(0,s),b0(s)=b(0,s); these
are not both zero. On x=0,s=kappa the y component, multiplied
by -(1-z^2), is exactly the JS polynomial

    2a0(kappa)z(1-z^2)^2
                       +b0(kappa)(z^2-2kappa z+1).       (2)

If D0(y) vanished identically on x=0, this polynomial would
vanish identically in z for generic kappa. Its degree-five
coefficient forces a0=0, and its remaining degree-two coefficient
forces b0=0, a contradiction. Thus x is not a common divisor.
D0 is primitive in the whole source ring.

For arbitrary nonzero a,b first remove their coefficient gcd;
this factor is invariant, so D and D0 have the same kernel and
the same embedded invariant field. JS already proves that D0
has an actual fixed affine line on every s=kappa,kappa!=+-1.
JT0 makes every other action sharing its invariant field an
invariant polynomial multiple fD0. This can only add fixed
points. If the prospective action is primitive it is even a
nonzero constant multiple. Interchanging x,y and transporting
by any K2-preserving source automorphism proves all of JT1.

This is a statement about the same embedded subfield of Frac(R),
including its value-preserving source-conjugacy saturation.
An abstractly isomorphic rational field or polynomial kernel
does not meet that hypothesis. A new candidate must change the
actual rational quotient, not just its generators or action speed.
No bound on a coefficient, replica or source degree is used.

## 5. Exactly what turns an action into a source certificate

Assume the JT2 global coordinate data are supplied. Write
H=PF+QG. The coprime-translation kernel proof in JS gives

    ker D=B[H].

The required equality with B[h] means both h is a polynomial
in H over B and H is a polynomial in h over B. Comparing
degrees in that one indeterminate, including the inverse
identity, gives

    h=D0(b)+c H, c a unit of B=C[b1,b2].

Every such unit is a nonzero complex constant. The full Hessian
in the actual new source coordinates has zero vertical block
and mixed block c D(F,G), so its determinant is c^4 J(F,G)^2.
The inverse of the supplied polynomial coordinate system is
the source map producing this actual unit-Hessian function.
No nonlinear Hessian invariance is asserted.

If (F,G) has a polynomial inverse, the last two gradient outputs
cF,cG recover b, after which the first two outputs solve P,Q
using the inverse Keller Jacobian. Conversely a polynomial
gradient inverse, specialized at first two outputs zero and
last outputs c times a free target pair, gives a polynomial
inverse to (F,G): every preimage b solves the first two equations
by a unique P,Q, so both inverse identities follow. This proves
the gradient-inverse equivalence.

For the reverse implication of JT2, start with an actual unit-
Hessian source member. N>=2 excludes an automorphic gradient,
so the retained Tang normal form supplies precisely these full
coordinates and the Keller coefficient pair. Transporting the
whole coordinate system and its inverse to R, JS gives the
exact equality ker D=B[h]. This proves the equivalence.

It is essential that the SAME base coordinate pair occurs in
the translation formula and in B[h]. A polynomial three-variable
kernel with h as one abstract coordinate does not by itself
supply that pair as coordinates of the entire source ring, nor
two global linear translation variables, nor the everywhere
constant coefficient Jacobian. Generic slices, local sections
and formal inverse jets cannot replace these global data.

## 6. The narrowed next question

For K2 any hypothetical unit-Hessian source orbit requires a
primitive action with proper nonprincipal plinth ideal and no
height-one prime over it. Its embedded invariant field must be
outside every standard-module field and every K2-preserving
source conjugate of those fields. An action that survives must
still satisfy the global coordinate-action certificate, or be
connected to an actual source map by another proved recovery.

These are necessary exclusions and an exact certificate format,
not a construction of that action, a classification of all
actions fixing x, or full K2 inverse closure. Full h0,h1 and
S_rho source questions, K2 one-coordinate unit-Hessian routes,
full K2/h1 equivalence, both old endpoints, actual global
embedding recovery and unrestricted HC4/JC2 remain open.

We inspected [Dasgupta--Gaifullin v2](https://arxiv.org/abs/2306.00510v2)
and the primary published discussion of rank and the principal
plinth theorem in THREE variables. Its full-text publisher
redirect was unavailable through the browser, so no theorem
from that endpoint is imported. We also inspected
[Shakhmatov v1](https://arxiv.org/abs/2605.00138v1) for the
geometric role of the plinth locus, as context only. No
three-variable kernel or principal-plinth theorem is transferred
to four variables. Every new assertion used above has its own
elementary proof. Exact finite controls support the written
all-degree argument; they are not independent replay, formal
verification or external review.
