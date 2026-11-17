# A trace-dual pencil for the entire simple-pole ramified row

15 November 2026. PF2RULED53 strengthens PF2RULED52's trace inclusion to
an equality of modules and derives a symmetric matrix-pencil obstruction.
It keeps the actual finite map and every unramified omitted prime. The
result covers the whole complementary r2q1 row, not only the C*-invariant
B_n controls. It does not supply or exclude every finite Keller map.

## 1. Hypotheses and the coordinate primitive-element fact

Let pi:N->A2 be the original finite normalization of a complex plane
Keller map on U=A2_(x,y). Assume its sole ramified prime is R, of generic
index2. All other omitted primes D_1,...,D_n are unramified. Retain the
whole finite N; no primes are deleted. Suppose the complementary whole
coordinate charts have

    R: y=z^2, x=c/z+t, c!=0;
    D_i: x=lambda_i+w^s_i, with the retained full unramified chart.

Thus x has pole divisor exactly R, with a nonzero leading coefficient
at EVERY point of R, and is regular on all D_i. The D orders, poles,
boundary values and all coefficients are unrestricted. PF2OMIT7 and
the unramified boundary theorem make N smooth, all omitted curves A1,
and all residue degrees one. Put B=O(N), A=C[P,Q], K=C(P,Q), L=C(x,y),
d=[L:K]. A is the actual target polynomial algebra.

Every polynomial source coordinate is a primitive element of L/K. For
clarity here is the two-variable proof, rather than a generic-coordinate
assumption. On a line x=beta, F is a nonconstant polynomial map A1->A2,
finite over its image and immersive by the Keller differential. Its image
normalization is A1: the completed curve is rational, and finiteness with
the source's single point at infinity leaves exactly one target puncture.
The induced finite map A1->A1 is etale, hence has degree1, since a
higher-degree polynomial has a critical point. Apply this to the geometric
generic x fibre. It gives [L:K(x)]=1; likewise for y or any source coordinate.
This is the plane coordinate case of Formanek's field-of-fractions fact.
No symmetric-monodromy or zero-excess hypothesis is needed here.

The map x:N->P1 is globally regular: 1/x=z/(c+t*z) is a uniformizer
near every point of R, while x is regular everywhere else. Consequently
(pi,x) is finite onto its image W_x in A2 times P1 and is birational.
The normalization of W_x is N. Its finite-x chart has normalization
N minus R, which STILL includes every D_i. In particular x need not be
PF2OMIT3's source-finite admissible coordinate; its affine chart is not U.

## 2. The entire trace-dual module has two explicit B-generators

Identify the trace-dual with a fractional B-module in L:

    B^vee={q in L: Tr_(L/K)(q*B) subset A}.

Finite duality and the different identify it with O_N(R). One may check
this on the smooth finite map: the Jacobian different is exactly R, and
the inverse different is the trace-dual. Equivalently a regular canonical
form is q*pi*(dP wedge dQ), and the latter has divisor R. The primary
statements are [Stacks 0BW0](https://stacks.math.columbia.edu/tag/0BW0),
Lemma49.8.2, and [0BWB](https://stacks.math.columbia.edu/tag/0BWB),
Lemma49.12.6. They apply to this ACTUAL finite smooth normalization.

The sections 1,x generate this line bundle everywhere: 1 generates away
from R, and x, with its nonzero simple pole, generates at R. Since N is
affine, this proves the module equality

    B^vee=B+x*B.                                        (1)

This is stronger than Tr(x*B) subset A, and uses no source C* action.
It also gives the exact quotient

    B^vee/B = O_R(R),

as a B-module. This is an invertible module on R=A1, hence abstractly
free of rank1 there. Its pushforward to the target is the normalization
module of the branch curve C=pi(R), with a choice of trivialization.
The equality does NOT say x is integral, that B=A[x], or that N is the
monogenic inverse surface.

## 3. A symmetric pencil that must generate everywhere

B is finite locally free of rank d over A. On any target affine open
where an A-basis b_1,...,b_d is chosen, form

    G_0=(Tr(b_i*b_j)), G_1=(Tr(x*b_i*b_j)).              (2)

Both matrices have regular entries and are symmetric. This is a global
condition on the finite locally free module B, expressed in local bases;
no global freeness theorem is needed. Equation (1) is precisely

    [G_0 | G_1]: A^(2d) -> A^d is surjective.            (3)

In a trivialization, the maximal minors of this block matrix generate
the unit ideal. The condition is required at EVERY target point, including
singular image values. Merely having generic rank d does not suffice.

Let G_R define the reduced irreducible branch curve C. Since R is the
only ramified prime, of index2 and residue1, the normalized trace
discriminant is G_R. Thus, in a local basis,

    det G_0 = unit*G_R,
    det(T*G_0-G_1)=det(G_0)*Norm_(L/K)(T-x).             (4)

The second identity holds first over K, where G_0 is invertible, by the
characteristic polynomial of multiplication by x. Its left side has
regular coefficients. Because x is primitive, it is a degree-d equation
of x, not a proper power of an equation of smaller degree.

Write H_beta for the reduced image equation of the source line x=beta,
and H_i for pi(D_i). The source line is closed in affine N: adding its
only missing normalization point would complete A1 to P1 in an affine
surface. Hence it maps finitely with residue1 to its image, by the
argument in Section1. The full divisor is

    div(x-beta)=L_beta+sum_(lambda_i=beta) s_i*D_i-R.

Therefore the primitive polynomial F_x(T), normalized by leading
coefficient G_R, satisfies, up to nonzero constants,

    F_x(beta)=H_beta*product_(lambda_i=beta) H_i^s_i.     (5)

This allows globally coincident image supports and their multiplicities.
For beta outside the finite set of lambda_i the member is reduced
irreducible, with immersive A1 normalization. Its normalization still
must self-collide by PF2RULED38/52's one-line theorem. The family is not
asserted to be a linear pencil of target polynomials.

In particular

    det G_1 = unit*H_0*product_(lambda_i=0) H_i^s_i.      (6)

The trace Tr(x) is regular, so the T^(d-1) coefficient of F_x is divisible
by G_R. At the generic point of C exactly two roots go to infinity, and
the finite specialization has degree exactly d-2. The only pole is the
index-two R with x pole1; no other prime can be silently counted at infinity.

The projective coefficient ideal of F_x is the unit ideal, since W_x is
finite over the target. This and degree d-2 distinguish the full finite
model from an arbitrary polynomial with the same first and last coefficients.

## 4. The radical rank recovers all ramification excess

At a closed target value b, let

    V_b=ker G_0(b), rho_b=rank(G_1(b) restricted to V_b).

These quantities are basis-independent. The fibre algebra decomposes
into its local factors. At a point outside R the map is etale, so the
factor is C and the ordinary trace pairing has zero radical. At q in R
let ell_q be its length. Its local Artin algebra T_q is Gorenstein,
and the ordinary trace pairing has radical its maximal ideal m_q.

The functional f -> Tr(x*f), specialized to this factor, is a Frobenius
functional: x is a local generator of B^vee at q by its nonzero simple
pole and (1). In a local Gorenstein algebra with a Frobenius pairing,
m_q^perp is its one-dimensional socle. Since ell_q>=2, that socle lies
in m_q. The restriction of this pairing to m_q therefore has rank
dim(m_q)-1=ell_q-2. Summing the factors proves

    rho_b=sum_(q in R over b)(ell_q-2).                 (7)

In particular rho_b=0 precisely when every ramified point over b is a
fold. A self-collision of several immersive branches still gives rho_b=0;
it must not be confused with a nonimmersive point. A cusp of local length3
gives rank1, and length ell gives rank ell-2. This works at higher excess
without truncating the local model to an ordinary cusp.

The retained Orevkov Euler identity now becomes an all-degree matrix
constraint. With n unramified omitted primes,

    E_R=sum_b rho_b, d=n+3+sum_b rho_b.                 (8)

Only finitely many b contribute. Thus zero excess is equivalent to total
isotropy of ker G_0(b) for G_1(b) at EVERY target value. In the surviving
quartic packet n=1,d=4, this isotropy is compulsory everywhere. Equations
(3),(7),(8) are genuinely simultaneous finite-algebra tests; they do not
assume source symmetry or identify a target ruling.

## 5. Conductor compatibility with the same pencil

The raw binary discriminant/index identity and the reciprocal derivative
formula are already PF2OMIT4, not new theorems claimed here. Apply them
to W_x, retaining its finite-chart D_i. Write C_x for its divisorial
normalization conductor on N and J_x for its target index polynomial.
The graph is normal at generic R: in reciprocal coordinate 1/x, the
parameter has order1 and the projection has a generic fold. Hence
(C_x)_R=0. With f_x=F_x,T(P,Q,x), the reciprocal derivative formula gives

    div_N(f_x)=C_x-(d-3)*R,
    [C_x]=(d-3)[R], v_R(f_x)=3-d,                     (9)
    Disc_T(F_x)=constant*G_R*J_x^2,
    Norm(f_x)=constant*G_R^(3-d)*J_x^2.               (10)

Thus the same actual inverse polynomial, trace pencil and conductor must
have compatible pole order d-3. On the source f_x is a polynomial, with
its zero divisor the source part of the conductor. Under (8) its pole
order on R is n+E_R. The conductor can also contain D_i when inverse
coordinate values coincide generically over a shared image support.
It is not assumed to be supported only in U. All global support
coincidences remain allowed.

These formulas do not prove normality of A[x], finiteness of U over A[x],
or absence of conductor. In fact line-image collisions require conductor.
The earlier unregistered CENTRAL_INVERSE_CONDUCTOR draft in the omitted
divisor directory is not used as a proved input; only sealed PF2OMIT3/4
and the stated current claims are used.

## 6. Verification and the remaining finite-map gate

The [checker](verify_trace_dual_pencil.py) writes
[the receipt](trace_dual_pencil_receipt.json), with exact primary inputs
in [the manifest](sources/trace-dual-pencil/manifest.json).
It checks entire local trace pencils, the unit block-rank condition,
all-length Gorenstein radical-rank formulas, direct sums of several fold
and nonimmersive factors, normalized inverse poles and discriminants,
and the B_n weight-module description of B+xB. A non-Keller finite
control has the same two-generated trace-dual property and survives it;
this prevents treating the matrix condition as a Keller realization or
a universal exclusion. A generic-rank-only block matrix is another
explicit negative control.

The new equality and excess formula hold uniformly in field degree and
through the whole r2q1 complementary row, with arbitrary D orders and
coefficients. They are necessary constraints on the full finite algebra.
The other minimal R orders, common rulings and remaining quartic/higher
degree finite maps stay open. No growing Laurent search, actual Keller
pair, full exclusion, formalization, independent full reproof, external
review or novelty priority is claimed.
