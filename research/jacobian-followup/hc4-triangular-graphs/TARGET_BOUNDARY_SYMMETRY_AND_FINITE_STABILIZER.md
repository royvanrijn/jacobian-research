# Target boundary symmetry forces inversion; the full stabilizer is finite

HC4MYKY0/1. Full HC4/JC2 remains OPEN. The result is all-degree: a
noninvertible complex plane Keller map can have only a FINITE polynomial
automorphism group preserving its ENTIRE reduced nonproperness curve.
The proof lifts target actions, rather than assuming a source symmetry
descends. It does not prove that every nonproperness curve has a symmetry.

## 1. A complete finite-cover action lifting theorem: KY0

Let U be a connected smooth complex algebraic variety and V->U a
CONNECTED finite etale cover of degree d. Suppose G=Ga or Gm acts
algebraically on U. Then:

    a Ga action lifts algebraically to V;
    a Gm action lifts after pullback by lambda->lambda^d. (1)

The lift is normalized to be the identity at the group identity, so it
is a genuine group action. No Galois or primitive-monodromy hypothesis
is required. More precisely, the needed torus power is the order m of
a particular CENTRAL monodromy permutation, and m DIVIDES d. For
center-free monodromy m=1. The bound d is uniform in every generic degree.
This theorem is valid in ALL dimensions over C.

Here is the complete lifting argument. View V(C)->U(C) as a finite
topological cover, with transitive monodromy rho on its d sheets. For a
point u and a path g(s) from the identity to g in G, lift the orbit path
a(g(s),u) starting at v above u. Define the lifted action by its endpoint.
For Ga independence of the path follows from simple connectedness of C.

For Gm an identity-based generator loop in G gives an orbit loop in U
and a sheet permutation sigma. For any loop gamma in U, the action
square (s,t)->a(g(s),gamma(t)) gives a homotopy commuting the orbit
loop with gamma. Thus sigma is in the CENTER of rho(pi1(U)). A central
permutation in a transitive action is semiregular: if sigma^k fixes one
sheet, commutation and transitivity make it fix every sheet. All cycles
therefore have length ord(sigma), so

    m=ord(sigma) divides d.                         (2)

Under lambda->lambda^d, every ambiguity of a path in Gm acts by sigma^d=1.
The endpoint construction is well defined and continuous on Gm times V.
It covers the dth-power target action. Group multiplication follows from
uniqueness of covering lifts on the connected product G times G times V;
the competing maps agree at the identity. This also proves the identity
and inverse laws, not just pointwise existence of lifts.

ALGEBRAICITY is essential and is not inferred from an analytic map being
regular. Form the ALGEBRAIC finite etale cover of W=G times V

    E=W times_U V -> W,

using the map (g,v)->a(g,pi(v)) for Ga or a(g^d,pi(v)) for Gm. The
topological lifted action is a SECTION of this finite cover. Full
faithfulness of generalized Riemann existence makes this section algebraic.
This is an explicit NEW PROVED EXTERNAL INPUT: SGA1 TheoremXII.5.1,
as stated in Poonen, Rational points on varieties, Theorem3.5.41,
printedpp77-78/PDFpp91-92. We adopt that theorem; we do not reprove it.
No nonproper GAGA principle or algebraization of an arbitrary infinite
analytic covering is used.

## 2. Extension over a finite normalization and every omitted prime

Now let Y=Spec(B) be a normal affine complex variety, Z=Spec(S) its
normalization in a finite separable field extension, and pi:Z->Y finite.
Let D be an invariant reduced closed set with U=Y minus D smooth,
pi etale over U, and Z_U connected. Let a Ga/Gm action of Y preserve D.
Apply Section1 to Z_U->U.

The lifted action extends ALGEBRAICALLY to ALL of Z. Indeed G times Z
has the normal coordinate ring S[t] or S[t,t^-1]. On its dense open
G times Z_U the action gives rational pullbacks of elements s in S.
Each s satisfies a monic polynomial over B. Equivariance makes its
pulled-back coefficients regular on G times Z. The rational pullback
is therefore integral over its normal coordinate ring, hence regular.
This gives the morphism G times Z->Z globally. Equivariance and the group
laws hold on the dense open and hence everywhere by separatedness.
No smoothness or finite flatness of Z is needed in this step.

Every prime component E of pi^-1(D) is invariant under the lift. The
image of the irreducible G times E is contained in the finite union of
components of pi^-1(D); its irreducible closure must lie in ONE component.
At the identity it contains E, so it lies in E. Inverse group elements
give equality. Therefore ANY open obtained by deleting selected full
prime components of pi^-1(D) is invariant as well.

For an ACTUAL plane Keller map K:X=A2->Y=A2, take its ENTIRE reduced
nonproperness curve D=N_K and its full finite normalization Z. Off D the
whole finite cover is the source and is etale. Every omitted prime lies
over D, and the retained pure-divisorial boundary theorem KMFI94 says
X is EXACTLY Z minus those omitted primes. Consequently

    every Ga symmetry of (Y,D) lifts to X;
    every Gm symmetry lifts to X after the dth-power isogeny. (3)

The source action is nontrivial whenever the target action is: the
target power map is surjective and K is dominant. Its polynomial source
derivation is nonzero and locally finite, and preserves the target field.
The all-degree plane obstruction KX0 therefore forces K to be an inverse.

Preserving merely the RAMIFIED BRANCH support is NOT the hypothesis.
Unramified omitted primes must also be preserved. Taking the ENTIRE D
makes them components of the invariant finite preimage and proves exactly
the source-open invariance that is otherwise missing.

## 3. Complete target-curve consequences: KY1

For any hypothetical noninvertible plane Keller map,

    Aut(A2,N_K) has NO nontrivial connected algebraic subgroup. (4)

Equivalently, its reduced equation H admits NO nonzero locally finite
polynomial logarithmic derivation delta, that is,

    delta locally finite, delta(H) in (H) ==> delta=0. (5)

To justify the equivalence without a new classification assumption, use
the Jordan/weight argument proved in KX0. A finite delta-orbit subspace
inside the stable ideal (H) is stable under both Jordan parts. A nonzero
nilpotent part integrates to a Ga action preserving (H). Otherwise the
finitely generated complex weight group supplies a nontrivial Gm subgroup
preserving (H). Section2 excludes either. Conversely, Lie elements of
connected algebraic actions are locally finite and logarithmic; a
nontrivial connected action in characteristic zero has a nonzero one.

Now import Blanc--Stampfli, Automorphisms of the plane preserving a curve,
Algebraic Geometry2(2)(2015),193-213, Theorem1 (printedp194,PDFp2;
proof conclusion printedp205,PDFp13), DOI10.14231/AG-2015-009. For EVERY
reduced plane curve, including reducible curves, its FULL polynomial
ambient stabilizer is an algebraic group unless a polynomial automorphism
turns it into a FENCE (a union of parallel coordinate lines). This is a
second NEW PROVED EXTERNAL INPUT, not independently reproved here.

A fence has a nontrivial Ga symmetry and is already excluded by (3).
Thus Aut(A2,N_K) is an algebraic group. Its identity component is trivial
by (4), so it has dimension zero and is FINITE. We obtain the full result

    K noninvertible ==> |Aut(A2,N_K)|<infinity.       (6)

This excludes infinite discrete polynomial stabilizers too; absence of
Ga/Gm alone would NOT justify that stronger conclusion without the
Blanc--Stampfli theorem. It does not force the finite group to be trivial
or give a uniform bound on its order.

In particular, in ANY polynomial target coordinates, the WHOLE D cannot
be a weighted-homogeneous divisor for ANY nonzero integer weights,
including mixed or zero weights. Its Euler field would contradict (5).
It cannot be a union of fibres of one target coordinate, either.
The condition is on the ENTIRE reduced D: individual components may have
symmetries broken by other components. No symmetry is asserted to exist
on every possible nonproperness curve. Thus (6) is an actual all-degree
necessary restriction, not a proof of full HC4 or a generic degree bound.

## 4. Sharp controls and failed substitutions

**The torus isogeny is essential in the general cover theorem.** For
V=Gm_t times A1_v -> U=Gm_u times A1_v, u=t^m, target scaling of u
does not lift with the original Gm parameter when m>1: a lifted action
on the unit t would have integer character a with m*a=1. After
lambda->lambda^m it lifts by t->lambda*t. Its monodromy is an m-cycle,
attaining the bound in (2). This is an etale cover of a punctured target,
not a whole-plane Keller map.

CONNECTEDNESS matters to the bound m divides d. On five disconnected
sheets a permutation with cycles2 and3 has order6, which does not divide5.
It cannot be used as a connected-cover counterexample to (2).

**Branch invariance alone does not preserve the source open.** Let Z=Y=A2
with the identity finite cover, but remove the unramified curve v-u^2=0
from the source. The finite branch set is empty and invariant under
v-translation. The deleted curve moves. The lift to Z therefore does
not restrict to the chosen source open. This source is not A2; the
example controls the exact failed inference about arbitrary opens.

**The Keller hypothesis is essential.** The retained quasi-finite maps
(x,y)->(u=x,v=x^a*y^2-y), a>=1, have invariant nonproper curve u=0
under v-translation. They are not Keller: J=2*x^a*y-1, and the true
finite branch curve 1+4*u^a*v=0 lies OUTSIDE the nonproper curve. The
cover over Y minus N is not etale, so Section2 is inapplicable. This
credits the prior PF2OMIT3 control rather than inventing a counterexample.

**A component symmetry is not a whole-curve symmetry.** The weighted
cusp v^2-u^3=0 is invariant under 2u*partial_u+3v*partial_v. Adjoining
v-u-1=0 breaks that action on the union. Conversely nonlinear target
coordinate changes preserve the obstruction: the curve
(v-u^2)^2-u^3=0 has logarithmic field
2u*partial_u+(3v+u^2)*partial_v and the same Gm action after conjugation.
No search through only linear target symmetries establishes (6).

**A tangent vector field is not an algebraic flow.** Every equation H
has its Hamiltonian logarithmic field, but it need not be locally finite.
For H=u^2*v^2+u+v, its highest Hamiltonian part acts on u^a*v^b by
2*(b-a)*u^(a+1)*v^(b+1). Starting at u, every highest iterate is
(-2)^n*u^(n+1)*v^n, with strictly increasing degree2n+1. Thus this
logarithmic field is not locally finite. It does not satisfy (5)'s
missing hypothesis merely by being tangent to a divisor.

The known whole3D Keller map is equivariant for Euler weights(1,-1,-2)
on its source and(-2,-1,1) on its target. Its full nonproperness set
is therefore invariant. This is consistent with the lifting theorem
in all dimensions, but shows that the PLANE conclusion (4)--(6) cannot
be transferred directly to arbitrary Keller maps or to an HC4 gradient.

## 5. Full target, dependencies and assurance

The new inputs are exactly generalized Riemann existence and the published
Blanc--Stampfli stabilizer theorem. Orbit-loop monodromy, its degree bound,
normalization extension and whole-source invariance are written directly.
The final plane inverse implication uses the RETAINED KX0 theorem, with
its existing Rentschler/Gutwirth/Gwozdziewicz assurance limits. No priority
claim is made for standard lifting facts or published classification.

FULL THEOREM TARGET remains full HC4/JC2 proof or an ACTUAL counterexample.
CURRENT LAST UNSOLVED IMPLICATION remains

    ACTUAL polynomial plane source
       ==> Tr(Lambda_L*Lambda_L) subset B.          OPEN

The alternative full marked canonical-energy integrality gate remains
open. No argument here supplies a target symmetry for every nonproperness
curve; curves with finite stabilizer remain possible. No full surviving
higher-degree packet is newly ruled out merely because it has no obvious
symmetry. Quartic realization remains retired by Domrina. External review
and outreach remain DEFERRED, and all assurance flags remain unchanged.
