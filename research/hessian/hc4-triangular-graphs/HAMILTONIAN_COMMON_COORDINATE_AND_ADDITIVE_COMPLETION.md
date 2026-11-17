# A common regular coordinate extends by its Hamiltonian field

HC4MYHCC0. 17 November 2026. A conditional bridge for actual complex
plane Keller completions, with no boundary-chart classification.

Let U=Spec(A) be a smooth integral AFFINE complex surface containing
the marked source X=Spec(R)=A2 as a dense open, so A is a subring of
R with the same fraction field. Assume a nowhere-zero regular two-form
Omega on U restricts to the source volume. Suppose ONE polynomial
source coordinate y has nonnegative valuation at EVERY boundary
prime divisor of U minus X.

Then y belongs to A. Its Hamiltonian field is a nonzero locally
nilpotent derivation on A and restricts to the source translation
partial_x in coordinates (x,y). Its kernel is exactly C[y], and
its additive action preserves the whole marked source X.

This conclusion does not need whole normalized root charts, boundary
components isomorphic to A1, disjointness of those components, simple
connectedness, or a global primitive of Omega. It does not supply the
missing common coordinate or establish unique LND kernels, all-action
source preservation, or cyclicity of the generic deck group.

## 1. Exact hypotheses and the weakened unresolved input

Fix source coordinates (x,y) with R=C[x,y]. Rescale the given volume
form by a nonzero constant so that

    Omega|X=dx wedge dy.

The hypotheses are precisely:

    U=Spec(A) is a smooth integral affine complex surface,
    X=Spec(R) is a dense open in U,
    A subset R, Frac(A)=Frac(R)=C(x,y),
    Omega is regular and nowhere zero on U,
    y is a polynomial coordinate of the MARKED source,
    v_D(y)>=0 for every prime divisor D contained in U minus X. (1)

The final clause is a DIVISORIAL common-coordinate condition.
No particular local equation, normalized cover, or boundary chart
is assumed at a closed boundary point.

For an actual plane Keller map, the pullback of the target volume is
a regular unit two-form on its etale locus. Thus, WHEN the affine
completion/open-source hypotheses in (1) are established for that
actual object, this bridge applies directly. It does not prove
affineness of the entire etale locus or any other completion assertion.

The missing divisorial common-coordinate existence is weaker than
the bundled HESS-COMMON-PENCIL-INPUT, which also asks for whole root
charts. This note proves only the implication from (1).

## 2. Normality makes the coordinate regular everywhere

A is a Noetherian normal domain. The standard height-one intersection
theorem gives

    A=intersection_(height p=1) A_p inside Frac(A).        (2)

This is the exact statement of [Stacks, Lemma10.157.6(2), Tag031T]
(https://stacks.math.columbia.edu/tag/031T), whose statement and proof
were inspected. It is used with the actual ring A.

At a height-one point whose generic point lies in X, y is already
regular because it is a polynomial on X. At every remaining height-one
point, (1) says precisely that y belongs to its DVR A_p. Equation (2)
therefore yields

    y in A.                                             (3)

This also handles every closed boundary point in codimension two:
normality leaves no extra isolated obstruction to regularity of a
rational function. Boundary components may meet, have singularities
as curves, or fail to be affine lines.

## 3. The unit two-form extends the source translation directly

On the smooth surface U, contraction with Omega is an isomorphism
of locally free sheaves

    T_U -> Omega_U^1, delta -> contraction_delta(Omega).

The inverse exists at EVERY point because Omega is nowhere zero.
Apply it to the regular one-form dy. This defines a global regular
vector field delta with

    contraction_delta(Omega)=dy.                         (4)

On X, contraction_(partial_x)(dx wedge dy)=dy, so uniqueness gives

    delta|X=partial_x.                                  (5)

Since U is affine, delta induces an algebraic derivation A->A.
Its extension to C(x,y) is exactly the source field derivation
partial_x. In particular no extension through a branched covering,
no formal vector field, and no codimension-two chart gluing is needed.

For any a in A, the source inclusion a in C[x,y] gives a finite
x degree. Thus

    delta^(deg_x(a)+1)(a)=0.                             (6)

Consequently delta is LOCALLY NILPOTENT on A. It is nonzero: if
it killed A, it would kill Frac(A)=C(x,y), contrary to delta(x)=1.

The kernel is exact:

    ker(delta|A)=A intersect ker(partial_x|R)
                =A intersect C[y]=C[y],                 (7)

where the last equality uses y in A.

## 4. Source preservation and every polynomial replica

The finite exponential on each a in A defines an algebraic Ga action
on U. For a polynomial source function a, its exponential equals the
ordinary substitution a(x+s,y). Hence the action restricted to
Ga times X is the polynomial source translation

    (s,x,y) -> (x+s,y).

Equality holds as morphisms to the affine surface U because their
pullbacks agree on every element of A. The action therefore preserves
the ENTIRE marked source X. Its complement is invariant as well.
Each irreducible boundary component is invariant, since the connected
group Ga cannot permute a finite set nontrivially.

For every polynomial q(y), the derivation

    q(y)*delta

preserves A and is locally nilpotent, because delta(q)=0 and
(q*delta)^n(a)=q^n*delta^n(a). Its source action is

    (x,y) -> (x+s*q(y),y).                               (8)

This proves the existence and explicit form of THESE additive actions.
It does not classify every LND on A, establish that every LND kills y,
or make a statement about rational deck transformations.

Affineness is needed when passing from the LND of global functions to
a complete action on U. For example the regular field partial_x on
A2 minus {(0,0)} has an LND on its global ring C[x,y], but its
translation orbit through (1,0) exits that open at parameter -1.
This example is not an object satisfying all of (1); it isolates the
role of the affine hypothesis.

## 5. Exact controls and barriers to stronger conclusions

A smooth affine control with a proper plane open is

    V: y*a=z*(z-1),
    z=x*y, a=x*(x*y-1) on the marked source A2_(x,y).

The source is the complement of the affine line

    D={y=0,z=1}, with free coordinate a.

Indeed on D(y) one recovers x=z/y, while on D(z-1) one recovers
x=a/(z-1). Their union is exactly V minus D and is the whole marked
source. The surface is smooth because the partial derivatives of
y*a-z*(z-1) have no common zero on V.

The two-form

    Omega=dz wedge dy/y
          =da wedge dy/(2z-1)

is regular and nowhere zero on the cover D(y), D(2z-1), and pulls
back to dx wedge dy. Its Hamiltonian for y is

    delta(y)=0, delta(z)=y, delta(a)=2z-1.

It is locally nilpotent and its exact action is

    y -> y,
    z -> z+s*y,
    a -> a+s*(2z-1)+s^2*y.                              (9)

These formulas preserve the defining equation, the source open,
and D, agreeing with x->x+s. No actual etale target pair or Keller
counterexample is asserted for this control.

It also shows the limit of the bridge. There is a second LND

    epsilon(a)=0, epsilon(z)=a, epsilon(y)=2z-1,

with a different kernel. Its action can move D into the source open.
Thus the present weaker hypotheses do not establish unique kernels
or preservation of the marked source by ALL additive actions.
Those stronger actual-Keller claims need their separate hypotheses.

Two other precise barriers are elementary. A regular function that
is not a source coordinate is insufficient: on A2 with its standard
volume, the Hamiltonian of x*y is x*partial_x-y*partial_y, which
is not locally nilpotent. A degenerate two-form is insufficient:
for Omega=u*du wedge dv, the Hamiltonian of v is u^-1*partial_u,
which has a pole at u=0.

## 6. Consequence for the common-pencil bridge and replay

For additive completion alone, whole-boundary root-chart hypotheses
may be replaced by the last, purely divisorial clause in (1).
Regularity at all closed boundary points then follows from normality
and the unit-volume Hamiltonian construction.

The existence of that common source coordinate is NOT proved here.
The missing PF46/50 valuation/classification content, any affineness
or full-completion prerequisites, uniqueness of LND kernels, free
deck actions and generic deck cyclicity remain separate matters.
No bundled HESS-MP theorem, HC4 or JC2 is declared proved by this note.

Replay the exact algebraic controls with

    .venv/bin/python research/hessian/hc4-triangular-graphs/verify_hamiltonian_common_coordinate.py

The checker validates contraction signs, coordinate conjugation,
nilpotence and exponential identities, the smooth affine control,
its proper plane chart and second action, and the explicit barriers.
These finite controls do not replace the normal-domain and vector-bundle
proof above. Internal derivation and replay are not formal verification
or external human review.
