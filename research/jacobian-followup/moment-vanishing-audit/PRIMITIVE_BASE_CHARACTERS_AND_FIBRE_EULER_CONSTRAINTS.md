# Primitive base characters and fibre Euler constraints

16 November 2026. Written all-degree proof, with no degree search or CAS
construction. Independent replay, formal verification and external review
are not claimed. Both FULL exceptional MZ statements remain **OPEN**;
external review is deferred. This supplies an input to the full
pure-to-safe/intrinsic-extension obligation, not its solution.

## 1. Every invariant phase has a canonical primitive affine-line base

Put X=SL2(C), A=C[X], and let Gamma be ANY finite subgroup acting by
right multiplication. For every nonconstant f in A^Gamma there are
regular q in A and nonconstant H in C[T] such that

    f=H(q),
    C(q) is relatively algebraically closed in C(X).       (1)

Thus q has geometrically integral generic fibre. The ring C[q] is the
integral closure of C[f] in the relative algebraic closure of C(f) in
C(X); in particular it is canonical. Its generator q is unique up to
an affine change. The following proof addresses the nonproperness of X:
Stein factorization is used only AFTER proper compactification.

First, X is smooth, rational and has only constant invertible regular
functions. For the last assertion use the projection to a nonzero bottom
row (z,w). Its fibres are affine lines, with local trivializations over
z!=0 and w!=0. A unit and its inverse are constant on each affine line,
and hence descend to a unit on A2 minus the origin. Regular functions
there are C[z,w], by intersecting the two localizations, whose only
units are constants.

Resolve a projective graph compactification of f:X->A1 to get a smooth
projective rational threefold Y containing X as an open set and a
morphism Y->P1 extending f. Apply proper Stein factorization:

    Y -> Cbar -> P1,

with geometrically connected fibres in the first map and a finite map
in the second. Here Cbar is a normal projective integral curve. Its
function field is the relative algebraic closure just described. These
are the Noetherian proper statements of
[Stacks, Theorem 37.53.4](https://stacks.math.columbia.edu/tag/03GX).
Graph compactification and resolution are retained characteristic-zero
inputs; no Stein theorem is being applied to the original affine map.

Because Y is rational, Cbar is dominated by a rational variety. Restrict
the resulting nonconstant rational map from A3 to a suitable affine line
to obtain a nonconstant rational map P1 to Cbar. It extends to a finite
morphism of smooth projective curves. Riemann-Hurwitz implies Cbar has
genus zero, so Cbar is P1. Only this elementary genus consequence of
[Riemann-Hurwitz](https://stacks.math.columbia.edu/tag/0C1B) is used.

Let B be Cbar with the finite map's poles removed. X maps dominantly
to B, since f is finite valued on X. At least one pole is removed. If
two or more were removed, B would have a nonconstant invertible regular
function (a ratio of linear coordinates for two removed points), whose
pullback would be a nonconstant unit on X. This is impossible. Exactly
one point is removed, so B=A1. A coordinate q on B pulls back to a
regular function on X, and the finite map B->A1 is H(q), a polynomial.
This proves (1) and identifies C[q] with the normalization ring above.
Equivalently the generic field extension over C(q) is regular; its
geometric generic fibre is integral. No claim about every special fibre
being connected or irreducible is made.

The map q:X->A1 is surjective: otherwise q-a would be an invertible
nonconstant regular function for an omitted value a, by the affine
Nullstellensatz. This also rules out using a rational function with a
finite pole as an unrecorded regular recomposition of q.

## 2. The only possible base action is a finite character

Since Gamma fixes f, it preserves the relative algebraic closure and
its normalization ring C[q]. It therefore acts on this ring by affine
transformations q -> alpha_gamma q+beta_gamma. A finite affine group
has a fixed point, obtained by averaging an orbit. After shifting q,
write q0=q-a; the action becomes

    gamma(q0)=chi(gamma) q0,
    image chi = mu_e,     K=ker chi.                      (2)

The image is cyclic, so [Gamma,Gamma] is contained in K. Since H(q)
is invariant, there is a polynomial G with

    f=G(q0^e),    r=q0^e in A^Gamma.                     (3)

In particular deg H is a multiple of e. The function r is primitive
on X/Gamma: over a generic value its pullback has e connected q0-fibres
which Gamma permutes transitively, and their quotient is connected.
Thus (3) gives the primitive affine-line factor on the ORIGINAL
quotient, not only on its covering space. No purity reflection from
G(r) to r is asserted.

For the two exceptional groups:

* Gamma=2I is perfect. Hence e=1: EVERY primitive root q can be chosen
  Gamma-invariant. A primitive phase on SL2/2I pulls back with connected
  generic fibre on SL2.
* Gamma=2O has commutator 2T and abelianization C2. Hence e=1 or2.
  Either q0 is invariant, or it transforms by the unique nontrivial
  character and q0^2 is the primitive phase on the quotient. A primitive
  quotient phase has at most TWO generic fibre components upstairs.

For completeness the needed commutator facts follow without a new
group computation. The central quotients are S4 and A5. The commutator
of 2O maps onto A4 and is contained in its order24 preimage 2T. If it
did not contain -I it would be an isomorphic A4 subgroup of SL2, which
is impossible: SL2 has a unique element of order2, namely -I, while
A4 has several involutions. Thus the commutator is 2T. Similarly the
commutator of 2I maps onto the perfect group A5; a putative index2
commutator would be an A5 subgroup avoiding -I, again impossible by
its involutions. Thus 2I is perfect. The facts S4'=A4 and A5'=A5
are elementary: 3-cycles generate their alternating groups, and in A5
each 3-cycle is a product of two double transpositions, so its
abelianization has both exponent dividing3 and generators killed by2.

## 3. Euler divisibility for the primitive fibres

For EVERY a in C, K acts freely on q0^{-1}(a): right multiplication
is free on all of SL2, and (2) makes K preserve this fibre. The affine
quotient is a finite etale covering of degree |K|. Consequently

    chi_c(q0^{-1}(a)) = |K| chi_c(q0^{-1}(a)/K).           (4)

This includes singular and reducible special fibres. The elementary
topological input is finite semialgebraic triangulation: triangulate
the quotient into finitely many locally closed contractible cells;
the covering is trivial over each cell and supplies |K| copies.
Additivity of compactly supported Euler characteristic proves (4).
No smoothness or properness of these affine fibres is required.

Thus EVERY primitive fibre Euler characteristic is divisible by24
for 2O, and by120 for 2I. In the trivial-character 2O branch the
stronger divisor is48. In general the divisor is |Gamma|/e, not an
unjustified |Gamma| for a root on which Gamma acts nontrivially.
The weaker statement that |[Gamma,Gamma]| divides the fibre Euler
characteristic is intrinsic and uniform in ALL degrees.

In particular, if a regular primitive q has ONE fibre of Euler
characteristic 1 or -1, no nonconstant H(q) can be invariant under
a noncyclic freely acting finite group Gamma: (4) forces K trivial,
and (2) then makes Gamma cyclic. This argument also applies to q
composed with ANY algebraic automorphism of SL2, without a tame
or coordinate-component hypothesis.

## 4. All automorphic recompositions of the small Long phase are excluded

Use Long's convention X=[[a,c],[b,d]], ad-bc=1, and put

    L=(1+c)(ad+b).

Set t=1+c. Then L=t+b t^2. For s!=0 the fibre L=s is exactly

    t in Cstar,    b=(s-t)/t^2,
    ad=((s+1)t-s)/t^2.                                   (5)

For s outside {0,-1}, the right side has exactly one simple zero,
at t=s/(s+1). Away from this point the fibre over t is Cstar and
has Euler characteristic0. At the zero it is the two coordinate axes
ad=0, with Euler characteristic 1+1-1=1. Hence

    chi_c(L^{-1}(s))=1,    s outside {0,-1}.              (6)

The function field is C(L)(t,a), with b,c,d rationally recovered by
(5), so C(L) is relatively algebraically closed in C(X). Thus L is
primitive. The preceding obstruction proves:

    H(L composed with phi) is NOT Gamma-invariant

for EVERY nonconstant polynomial H, EVERY algebraic automorphism phi
of SL2, and EVERY noncyclic finite Gamma acting freely on the right.
In particular this rules out the ENTIRE automorphic-recomposition
construction on both exceptional quotients, in all degrees.

The same proof for a matrix coordinate uses its nonzero fibres A2,
of Euler characteristic1. It extends the earlier coordinate obstruction
to all polynomial recompositions, rather than appealing to a degree
reduction theorem for automorphism tuples.

This is an exclusion of invariant functions in a specified architecture,
not a proof of purity or nonpurity for every invariant phase. It does
NOT apply to the regular Bezout maps Phi_b, whose Jacobians are
nonconstant and which are not automorphisms. The retained
[weighted transfer](REGULAR_BEZOUT_MAPS_AND_WEIGHTED_HAAR_TAILS.md) and
[opposite-weight Long closure](QUADRATIC_SUSPENSION_AND_ALL_LONG_RECOMPOSITIONS.md)
retain their exact, separate scopes. Mixed Bezout parameters are not
excluded by (6).

## 5. Safe controls and the remaining whole-quotient obligation

The Euler restriction allows the retained nonzero pure SAFE phases.
For a homogeneous binary invariant F_n with n distinct roots,
q=F_n(x,y) has generic fibres that are affine-line torsors over
F_n(x,y)=s in A2. The latter curve is an n-sheeted etale cover of
P1 minus the n roots. Thus its Euler characteristic, and that of the
threefold's q-fibre, is n(2-n). For the retained n=8 and n=12 forms
this is -48 and -120 respectively, satisfying (4). Their roots are
distinct: the degree8 form reduces to a quadratic in X^4/Y^4 with
distinct nonzero roots; the degree12 form has the two simple coordinate
roots and a quadratic in X^5/Y^5 with distinct nonzero roots.
Automorphism transport preserves these fibre Euler characteristics.

There is no contradiction with pure safety, and no Euler-sign argument
for intrinsic splitting has been established. Euler characteristic0
also gives no divisibility exclusion. These controls prevent promoting
(4) into a full MZ theorem.

Every remaining genuine failure must obey this primitive-base character
and commutator divisibility constraint AS WELL AS the retained capacity,
automorphism, image-interior and five/seven-weight constraints. The
all-degree pure-locus inclusion, equivalently the intrinsic unit
extension splitting, remains the target. The
[moment-ideal amplification route](PURE_MOMENT_CONORMAL_AMPLIFICATION.md)
is still unproved. No new decidability claim, degree campaign or
external review has been introduced.
