# Source symmetry cannot descend to a nontrivial Keller target field

15 November 2026. PF2RULED52 closes a possible symmetry shortcut on the
PF2RULED51 affine models. It gives an all-degree obstruction already at
the target FUNCTION FIELD, explicit missing weight components, and a
trace/norm gate on those models. It does not exclude all Keller maps on
them. All other omitted primes are retained whenever a finite normalization
is used.

## 1. Recover the target algebra from the target field

Let F=(P,Q) be a complex plane Keller map, S=C[x,y], A=C[P,Q],
K=C(P,Q) and L=C(x,y). Then

    K intersect S=A.                                      (1)

Here is a proof that does not assume the desired finiteness of F.
Write a rational target function as g/h with coprime target polynomials.
If an irreducible H divides h, then H(P,Q) is nonconstant by algebraic
independence. It has a curve E of zeros. The Keller map is etale and
quasi-finite; thus E dominates the target curve H=0, and g(P,Q) does not
vanish identically on E. Consequently g(P,Q)/h(P,Q) has a pole on E.
No nonconstant denominator is possible for a source polynomial, proving
(1). The same proof works for a dominant quasi-finite polynomial map
of normal affine planes. A generically finite map with a contracted
curve does not suffice: (x,xy) has y in K intersect S but not in A.

We use Gwozdziewicz's Theorem1.1, *Injectivity on one line* (1993):
a Keller map injective on one affine line is an automorphism. The exact
two-page primary paper, including its proof and hypotheses, is retained
and visually checked in the [manifest](sources/source-symmetry-gate/manifest.json).
It was already used in PF2RULED38; no new line-embedding theorem is assumed.

In particular, if ANY element of A restricts to a nonconstant affine-linear
function on ONE source-coordinate line, F is an automorphism. The element
need not itself be a target coordinate. This criterion survives arbitrary
polynomial changes of source coordinates.

## 2. Additive and multiplicative source symmetries at field level

In any polynomial source coordinates, either of the following conditions
forces F to be an automorphism:

    partial_x(K) subset K;
    E_(a,b)(K) subset K, E_(a,b)=a*x*partial_x+b*y*partial_y,
      (a,b) in Z^2 minus {(0,0)}.                         (2)

For the additive condition, (1) gives partial_x(A) subset A. Choose
f in A of positive x degree m, which exists by algebraic independence.
Its (m-1)-st derivative is in A and has form a(y)*x+b(y), a!=0.
On a line y=beta with a(beta)!=0 this is affine-linear in x. Apply the
one-line theorem.

For the multiplicative condition, (1) first gives E_(a,b)(A) subset A.
For a source polynomial f with finite weight support W, each homogeneous
component is obtained by the exact spectral projector

    f_w = product_(v in W,v!=w) (E-v)/(w-v) applied to f.   (3)

Hence every component of every f in A belongs to A. Suppose a!=0.
The invertible differential of F at the source origin gives some f in
{P,Q} with a nonzero linear x term. Its weight-a component restricts
on y=0 to c*x, c!=0: the only surviving monomial of that weight on the
axis is x. This is an element of A, so F is injective on that line.
If a=0, use b!=0 and the y axis instead. This proof includes positive,
negative and zero weight patterns; it requires no linearization theorem
for an action on the target and no homogeneity of the chosen P,Q.

Thus a noninvertible Keller target field is preserved by neither a source
translation derivation nor a nonzero integral diagonal Euler derivation,
in ANY polynomial source coordinates. A regular translation or diagonal
C* action preserving K would preserve A by (1); differentiating its
polynomial action then gives (2), again impossible.

This is an obstruction to DESCENT. The source action itself may exist.
The geometric action on PF2RULED51's N_n is not ruled out. A finite map
alone does not make its target field invariant: the finite non-Keller map
(x,y)->(x+x^2,y) gives a concrete failure, detected by its deck involution
x->-1-x. Conversely a polynomial automorphism can have mixed-weight
coordinate expressions even though its target algebra is all S and hence
invariant. Do not identify invariant algebra with homogeneous generators.

## 3. Explicit field escape on every PF2RULED51 model

Keep the marked smooth affine surface

    B_n=C[h,x(h-C),yH(h),xyH(h)], h=x^2*y,
    H(h)=product_i(h-d_i), C,d_i nonzero and pairwise distinct.

Suppose pi:Spec B_n->A2 is an actual finite extension of a Keller map on
the original marked source. Its degree is arbitrary. The source action
has weights (x,y)=(1,-2), and E=x*partial_x-2y*partial_y preserves B_n.
An affine target change arranges F(0)=0 and DF(0)=identity. From the whole
weight modules of PF2RULED51 the indicated two components are

    u=P_1=x*F_1(h), (h-C) divides F_1(h), F_1(0)=1;
    v=Q_(-2)=y*G_1(h), H(h) divides G_1(h), G_1(0)=1.     (4)

Both u and v are regular on the ENTIRE N_n. Nevertheless

    u notin K, v notin K.                               (5)

If u were in K, (1) would put it in A, and its restriction u(x,0)=x
would make F injective on y=0. Similarly v(0,y)=y would force injectivity
on x=0. The map would then be an automorphism, incompatible with a finite
normalization having the displayed nonempty omitted boundary. This proves
(5); it is stronger than saying merely that u,v are absent from A.

At least one of E(P),E(Q) also lies outside K. If both belonged to K, the
derivation rule on rational expressions would give E(K) subset K,
contradicting (2). In particular the source Euler vector field cannot
be represented by a rational target vector field, even allowing target
poles. All these functions are integral over A when regular on finite N_n.

There is a separate support consequence. For EVERY target polynomial
coordinate Z=Z(P,Q), its restriction to every source-coordinate line has
degree at least2. A degree1 restriction directly gives injectivity. If
the degree is0, a complementary target coordinate W has everywhere
nonzero derivative along that line by the Keller differential; its
one-variable derivative is then a nonzero constant, again giving
injectivity. For weights (1,-2), the two axes therefore require

    max weight(Z)>=2, min weight(Z)<=-4.                 (6)

This holds after every target polynomial automorphism. It forbids all
one-sided and narrow-weight coordinate attempts, at arbitrary degrees;
it does not bound the weight width from above or force equivariance.

If the finite extension has NO proper intermediate fields, every u,v
in (5), and any Euler derivative just shown outside K, is a primitive
element of L/K. In particular this applies under PF2RULED46's retained
hypothesis: the sole ramified prime R has index2, residue degree1 and
immersive image normalization (zero excess), and other omitted primes
are unramified. The actual monodromy is then S_d. In the Galois closure
the subgroup for L is S_(d-1), maximal in S_d, so there is no intermediate
field. A short proof of maximality: if a subgroup contains S_(d-1) and
moves its fixed sheet, it contains two distinct point stabilizers; their
transpositions generate the whole S_d. Thus the norms of T-u and T-v
are their monic minimal polynomials of FULL actual degree d. Without the
no-intermediate-field hypothesis an outside element need not be primitive;
the fourth-power field control in the checker retains that failure.

The three primitive choices, when applicable, obey the exact coupling

    u^2*v=h*F_1(h)^2*G_1(h).                            (7)

PF2RULED51 already makes h primitive under its own verticality argument,
which does not require zero excess. Formula (7) supplies compatibility
between these choices, not a new claim that A[h,u,v] is normal.

## 4. Rational source-coordinate norms and a trace constraint

On N_n the full source-coordinate divisors are

    div(x)=L_x+sum_i D_i-R;
    div(y)=L_y-2*sum_i D_i+2R.                          (8)

Here L_x,L_y are the retained axes in h=0. They and the unramified D_i
map finitely, immersively and with residue degree1 to their target image
normalizations: a higher-degree finite etale map A1->A1 is impossible.
For this paragraph assume also residue degree1 on R (in particular the
zero-excess hypothesis above suffices). Let H_x,H_y,H_i,G_R define their
irreducible target images. Finite divisor pushforward gives, up to nonzero
constants, the identities of RATIONAL target functions

    Norm(x)=H_x*product_i H_i/G_R;
    Norm(y)=H_y*G_R^2/product_i H_i^2.                   (9)

Coincident global image supports are allowed and can cancel in (9).
These are not asserted to be polynomial norms of integral elements.
Multiplying gives Norm(h)=Norm(x)^2*Norm(y)=const*H_x^2*H_y, agreeing
with PF2RULED51's monic norm at zero. Equations (4),(7),(9) must hold in
ONE actual finite extension, with the same target image curves.

There is also an integrality statement despite x having a pole on R:

    Tr_(L/K)(x*B_n) subset A.                           (10)

At the generic point of the branch image, the R factor has tame index2
and x has pole order1. In a completed fold parameter w, its trace pairs
w with -w, cancelling the only negative Laurent term of x*f for every
regular f. Every other finite local factor has x*f regular. At other
target prime divisors all factors are regular. The trace therefore has
no codimension-one target pole, so lies in normal A. No assertion that x
itself is integral follows. The local trace of x^2 need not be regular;
this is an explicit negative control. Equivalently (10) places x in the
trace-dual module of the actual finite algebra.

## 5. What this closes and what remains

The source grading of the affine controls cannot descend, even rationally,
to a nontrivial Keller target field. The same holds for source translations.
Any search or argument that assumes such invariance has already excluded
the desired noninvertible map. The missing weight components (4), coupled
norms and trace-dual condition give explicit finite-algebra tests that do
not assume descent. They apply uniformly in n and in field degree.

The [checker](verify_source_symmetry_gate.py) writes the
[receipt](source_symmetry_gate_receipt.json). Exact controls cover spectral
projectors, source weight modules and axis restrictions, automorphic and
non-Keller counter-controls, primitive-field necessity, divisor/norm
identities and tame trace cancellation. The universal conclusions rest on
the written arguments and the imported one-line theorem, not on the finite
control set. The 2026 graded-Keller preprint is related context only; its
abstract is not imported as proof of field-level descent rigidity.

Common rulings, the remaining minimal complementary configurations and
the full quartic/higher-degree goal remain open. No finite Keller pair,
normal monogenic algebra, full exclusion, independent full reproof,
formal verification, external review or novelty priority is claimed.
