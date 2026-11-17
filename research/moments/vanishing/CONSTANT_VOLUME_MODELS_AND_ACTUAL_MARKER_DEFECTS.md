# Constant-volume models and actual marker defects

16 November 2026. Internal written proof; no independent replay, formal
verification or external review is claimed. Both FULL exceptional MZ
cases remain **OPEN**. Review is deferred and is not the next task.

A finite constant-volume correspondence transfers existence of a proper
model carrying the ORIGINAL top form in both directions. Applied to the
retained weighted-torus model, this gives such a model on the original
SL2 source, and hence on each original exceptional quotient. This is a
stronger geometric conclusion than existence of an auxiliary moment model.
It does not give extension of every original marker form.

There are actual exceptional phases for which the original-volume model
exists and the unit Cauchy submodule splits globally, while the FULL
original marked-function module has a nonsplit unit logarithm. The fixed
original marker is explicit. These phases are already known NOT to be
eventually pure. The examples rule out an unconditional replacement of
the full marked source by the proper-model regular-period image; they
neither give an MZ counterexample nor settle the purity-sensitive bridge.

See the [checker](verify_haar_source_models.py),
[receipt](haar_source_models_verification.json), and
[source audit](haar_source_models_source_audit.json).

## 1. A valuative criterion for the original form

Let X be a smooth irreducible quasi-projective complex variety of dimension
N, f a regular function, and Omega a nonzero regular algebraic top form.
Write V(X,f,Omega) for the following property: X embeds as an OPEN subset
of a smooth Y with Y -> A1 PROJECTIVE extending f, and the SAME rational
form Omega extends holomorphically on Y. The condition is equivalent to

    ord_v(Omega)>=0 for EVERY divisorial valuation v of C(X)
                         for which f belongs to O_v.             (1)

Necessity: projectivity gives a centre on Y for every such v. Extract v
on a common smooth resolution over Y. Pullback of a regular top form
along a morphism of smooth varieties is regular, so its order is
nonnegative. Sufficiency: take a projective graph compactification of f
and resolve its boundary, keeping X unchanged. Restrict over A1. Every
prime divisor has finite f and hence nonnegative Omega order by (1).
On this smooth variety, a rational section of the invertible canonical
sheaf is regular exactly when it has no divisorial pole. This proves
sufficiency. Characteristic-zero resolution and the proper valuative
criterion are the same standard inputs as in the retained graph proofs.

This is a model-independent criterion. A negative order cannot be repaired
by choosing a different birational compactification. No assertion about
all marker multiples Q Omega is included in (1).

## 2. Finite constant-volume correspondences

Let K0 subset K1 be a finite extension of complex function fields of equal
transcendence degree N. For nonzero top forms suppose

    Omega1=c*pullback(Omega0), c in C*,
    f1=pullback(f0).                                      (2)

A divisorial valuation v of K1 restricts to a positive integer multiple
e*w of a divisorial valuation of K0; conversely every divisorial w has
an extension v, again divisorial. For example, normalize a model
extracting w in K1. Finiteness preserves the codimension-one generic
points and the residue transcendence degree N-1. The same argument
starting from v gives the restriction assertion.

At these generic points the top-form formula is

    ord_v(Omega1)=e*ord_w(Omega0)+(e-1)
                =e*(ord_w(Omega0)+1)-1.                  (3)

Here is its local calculation. Choose lifts z1,...,z_(N-1) of a
separating transcendence basis of the residue field of w. They remain
separating after the finite residue extension. Choose uniformizers t,u
with t=u^e*a, a a unit. In the completed DVR the coefficient of du in
dt, modulo the dz_i terms, has order exactly e-1: its leading coefficient
is e times a nonzero residue, and e is invertible in characteristic0.
Wedging with dz1,...,dz_(N-1) removes the other terms. Multiplying the
form by t^k times a unit contributes e*k. This proves (3), including
arbitrary residue extensions in characteristic0. It is the familiar
local differential calculation in [Stacks, Section53.12, discussion
preceding Lemma53.12.1](https://stacks.math.columbia.edu/tag/0C1B),
with residue coordinates retained; no global curve-genus formula is used.

For integer k, e*(k+1)-1>=0 iff k>=0. Also f1 is finite at v iff f0
is finite at w. Combining (1)--(3), whenever smooth open sources with
regular original forms are specified, gives

    V(X1,f1,Omega1) iff V(X0,f0,Omega0).                  (4)

The correspondence may be rational away from its generic points; no
extension to a morphism between the chosen resolutions is required.
The criterion also applies to the smooth proper resolution of a singular
auxiliary model. It concerns the rational form on its function field.

In particular a LOGARITHMIC pole k=-1 stays EXACTLY -1 on EVERY finite
cover. A higher pole gets worse. Ramified base change cannot repair an
original-form pole by itself. Replacing the form by a Jacobian multiple
changes (2) and changes the unweighted Haar question. Statement (4)
alone does NOT transfer compact-cycle functionals or all marked periods;
that needs cycle information separately.

## 3. The weighted auxiliary model transfers to the ORIGINAL source

On SL2 let P=p+U+V with torus weights 0,a,-b, a,b>0, and U,V nonzero.
All three functions are invariant under the right exceptional group
Gamma=2O or2I. Set

    g=gcd(a,b), r=a/g, s=b/g, q=U^s V^r.

The quotient by the left torus is the smooth affine sphere Q2. In a
local fibre coordinate lambda the original volume is a fixed nonzero
constant times Omega2 wedge dlog(lambda). The retained
[weighted-cover theorem](WEIGHTED_TORUS_COVERS_AND_ALL_ZERO_WEIGHT_TRIPLES.md)
constructs a smooth projective model over A1_S of the auxiliary field

    W: A^s B^r=q(x), S=p(x)+A+B,
    Theta=Omega2 wedge dA/(r*A),                         (5)

carrying the SAME holomorphic Theta at every finite S. Its original
proof covers every degree and all collisions, not just smooth fibres.

There is a generically finite map over Q2

    psi: SL2 --> W, z |-> (pi(z), U(z), V(z)).           (6)

On a generic torus fibre, U=u0*lambda^a and V=v0*lambda^(-b).
Coprimality of r,s shows the fibre map onto (5) has degree g: its
kernel consists exactly of the g-th roots of unity. Thus (6) is dominant
and finite on function fields. Base differential terms vanish when
wedged with Omega2, giving

    psi*Theta=(a/r)*Omega2 wedge dlog(lambda)
              =g*Omega2 wedge dlog(lambda).             (7)

The ratio to the original SL2 Haar volume is a NONZERO CONSTANT. By
(4), the proper model of (5) yields V(SL2,P,Omega3). This argument does
not assume psi lifts as a morphism to the selected auxiliary resolution.

The finite etale quotient SL2 -> SL2/Gamma pulls the normalized quotient
volume back to a nonzero constant times Omega3. Applying (4) again gives
an ORIGINAL-volume proper model on each original exceptional quotient.
Every nonconstant polynomial H(P) has the same property: compose the
projective model map with the finite map H:A1 -> A1. This covers zero
and nonzero centres and every ramification multiplicity. It does NOT
infer eventual purity of P from purity of H(P).

The known all-H nonpurity of these opposite-charge phases is already a
retained result. The new conclusion here is the ORIGINAL-source geometry.
It does not cover arbitrary mixed phases: the retained
[global-bundle safe-pole controls](GLOBAL_TORUS_BUNDLE_AND_UNIT_INVISIBLE_INFINITY_POLES.md)
have permanent finite-value volume poles and therefore cannot admit
such proper models. Thus no universal mixed-phase model is asserted.

## 4. What the proper model says about the UNIT Cauchy submodule

For SL2=Q3 the dimension is odd. The retained
[proper-unit theorem](PROPER_VOLUME_UNIT_PERIODS_AND_COMPLEX_HOPF_CLOSURE.md)
therefore applies to the original model of Section3. BEFORE ANY PURITY,

    M1(P)=D_(C(s))<1,A1(s)>, A1(s)=h(1/(s-P)),

has a GLOBAL differential retraction to the unit, fixing literal1.
Here D_(C(s)) denotes closure under rational scalar multiplication and
differentiation. Finite quotient normalization leaves the Cauchy function
unchanged, so the assertion holds on the original exceptional quotients.
The same is true for EVERY nonconstant H(P).

The precise larger object in that proof is the compact-cycle orbit
period image of the regular de Rham cohomology of Y minus Y_s. It is a
global unit-retaining quotient of a pure polarizable source E_Y in odd
dimension; hence its unit splits globally. The period quotient kills
the ambient kernel on every continuation, and contains A1 via the
regular form Omega_Y/(s-f). This does not place every Q Omega/(s-f)
in that image. It does not identify M1 with the full marked module.

## 5. Actual exceptional proper-volume phases with nonsplit marked unit

Use the regular Bezout families on BOTH exceptional quotients. In the
notation of the retained
[uniform-marker theorem](COHOMOLOGICAL_PULLBACK_AND_A_UNIFORM_ORIGINAL_MARKER.md),
(n,m)=(8,12) or(12,20), L=n+m, and b is any invariant section of weight L.
Write Y_b=V+G0*b, of weight n, while F0 has weight -n. The phase is

    g_b=(1+Y_b)*(1+F0+Y_b*F0)
        =p+A+B,
    p=1+2*Y_b*F0,
    A=Y_b*(1+Y_b*F0), B=F0.                             (8)

These charged functions are nonzero: F0 is nonzero, Y_b cannot vanish
because the Bezout matrix has determinant1, and Y_b*F0 cannot be -1
since the only units of C[SL2] are constants and the weights are nonzero.
Thus Section3 applies with a=b_charge=n and gives proper ORIGINAL-volume
models for g_b and ALL nonconstant H(g_b). Section4 splits M1 globally.

Nevertheless the SAME original regular invariant marker

    C_b=W_b*(1+Y_b), h(W_b)=1,
    h(C_b*R(g_b))=integral_0^1 R(-t)dt                  (9)

works for every polynomial R by the retained uniform-marker identity.
In particular its actual Cauchy germ at infinity is

    A_C(s)=log((s+1)/s),
    A_C'(s)=1/(s+1)-1/s.                               (10)

A local unit retraction on the FULL original marked-function module
would send A_C to a single-valued meromorphic primitive of its derivative,
because it fixes1 and commutes with differentiation. This is impossible:
the derivative has residues -1 at0 and +1 at-1, whereas the derivative
of a meromorphic function has residue0. Hence the full marked module is
locally unit-NONSPLIT at BOTH0 and-1, before any purity assumption.
This agrees with the retained actual component-transfer theorem.

Thus actual exceptional geometry realizes

    globally unit-SPLIT submodule M1(g_b)
      inside a FULL marked module locally unit-NONSPLIT.            (11)

A submodule inclusion is not a global unit-retaining quotient. The
retained quotient-reflection lemma therefore creates no contradiction.
Nor should (6) be confused with the Bezout map Phi_b: for (8), psi has
generic fibre degree n, whereas Phi_b has ACTUAL integral H3 scalar2
and normalized nonconstant Jacobian 2*W_b. That Jacobian cannot be
dropped from a Haar integral. Section3 uses a different correspondence
with the constant ratio (7).

There is also a cohomological consequence, stronger than seeing a raw
pole of C_b Omega. For ANY original-volume proper model Y of g_b,
A_C cannot belong to its actual regular de Rham compact-orbit period
image described in Section4, even modulo exact forms and actual GLOBAL
period relations. Otherwise the rank2 differential module generated by
1,A_C would inherit that image's global unit retraction, contradicting
(10). The marker loss is real in the period image, not merely notation
for a pole that might be removed by an exact form.

All phases in this equivariant branch, and all their nonconstant
polynomial recompositions, have ALREADY been proved not eventually pure
by the retained all-zero-weight-triple/all-Long results. They are not MZ
counterexamples. The same marker in (9) gives
h(C_b*H(g_b)^j)=integral_0^1 H(-t)^j dt for every j; its non-eventual
vanishing and global nonsplitting for all H are retained results. No
claim that every such moment is nonzero, or that primitive purity
reflects from a composite, is made here.

## 6. The remaining FULL obligation

The universal global bundle preserves every actual original marked
function and relation, and removes finite-base poles of the UNIT volume.
Sections3--5 now show concretely that even GLOBAL extension of that unit
volume does not permit unconditional replacement of the full marked
module by its regular-period image. Any proposed full reduction must
control the missing marker directions using the actual purity identities
and global relations, rather than discard them because Omega extends.

This result does NOT localize the marker defect exclusively over finite
base or exclusively over base infinity. A raw marker residue on one
boundary divisor can be exact; it is not by itself an intrinsic unit
obstruction. The finite-base-specific source elimination question remains
OPEN. Under eventual purity, the full task remains to force the corrected
local component/torus-valuation vector to vanish for EVERY original
exceptional phase in EVERY degree, or produce an actual invariant pure
phase with one fixed original bad marker. Restricted nonpure examples
and unit-only splitting do not prove either alternative.

The Gaussian/formal Bernstein effectivity programme remains separately
owned. No new degree, support or parameter search was run. Both FULL
SL2/2O and SL2/2I MZ cases remain **OPEN**; the goal remains active.
