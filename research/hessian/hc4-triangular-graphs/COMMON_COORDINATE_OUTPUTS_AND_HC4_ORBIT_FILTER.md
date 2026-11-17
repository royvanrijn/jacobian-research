# Common coordinate outputs, dimension descent, and the HC4 orbit filter

2026-11-16. Claims HC4MYLK0/1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, or full marked-energy integrality; OPEN.
External review/outreach is DEFERRED. The previous turn made PROGRESS
by proving the entire gamma source polynomially affine (LJ0/1).

The direct gamma-to-HC4 MAP conversion route is now CLOSED NEGATIVELY,
uniformly at every level and under arbitrary polynomial automorphisms
of source and target. Generic field degree, not the complexity of a
Darboux calculation, already obstructs it. Every global polynomial
2+2 mixed graph projection is excluded as well. A different construction
that changes the map/field degree is not excluded.

The general argument is all-degree and all-dimensional. It makes the
distinction between common coordinate outputs, parameter families, and
actual polynomial products explicit. The elementary coordinate-pair
restriction mechanism in BCR2 and the gradient fibre correspondence in
GS0 are PRIOR results; neither mechanism is claimed new here.

## 1. LK0: common coordinate outputs give genuine lower-dimensional Keller maps

Let K:A^n->A^n be a complex polynomial Keller map, of generic FIELD
degree d. Say that K has m common coordinate outputs if there are
polynomial source and target automorphisms putting it into the form

    K(z,t)=(G(z,t),t),
    z=(z_1,...,z_r), t=(t_1,...,t_m), r=n-m.          (1)

Equivalently, some m target functions are part of a polynomial coordinate
system, and their pullbacks are JOINTLY part of a polynomial source
coordinate system. A rational coordinate, primitive field element,
individually coordinate functions without joint completion, local slice,
or coordinate of an abstract isomorphic ring does not satisfy this.

Let kappa(K) be the largest such m, allowing m=0. This number is unchanged
by arbitrary polynomial source and target automorphisms. The notation is
only for an intrinsic maximum, not a proposed bounded coordinate search.

The block determinant in (1) gives

    det(D_z G)=j in C*.

Thus for EVERY lambda in C^m, G_lambda(z)=G(z,lambda) is an r-variable
polynomial Keller map. Its generic degree is d for all lambda in a
NONEMPTY ZARISKI OPEN set. This last qualifier must not be omitted.

Here is a complete specialization proof. Write B=C[v,t] and R=C[z,t]
via v=G(z,t). Since Frac(R)/Frac(B) is finite of degree d, invert one
nonzero h(v,t) so that R_h is finite free of rank d over B_h. Such h
exists directly: clear the coefficients of monic equations for the
finitely many source generators, then clear a field basis and its inverse
change-of-generators matrix. This uses the actual algebra R, not an
arbitrary rational model.

Outside the proper closed parameter set on which EVERY v-coefficient of
h(v,lambda) vanishes, this specialization is nonzero. The localized
base change is

    C[z]_[h(G_lambda(z),lambda)]
       finite free rank d over C[v]_[h(v,lambda)].

Its source is a domain and the localization is nonzero because the Keller
map G_lambda is dominant. Taking fraction fields proves degree d.
This is not a claim that every special parameter preserves degree.

Consequently, in ALL n and ALL d:

    n-r common coordinate outputs
       ==> d is realized by an actual complex r-variable Keller map. (2)

For r=1, the conclusion is stronger and global: partial_z G=j forces
G=j*z+b(t), with an explicit polynomial inverse. Therefore

    K nonautomorphic ==> kappa(K)<=n-2.              (3)

Use the already retained plane exclusions: PCC3/Orevkov excludes d=2,3;
KR0/Domrina excludes d=4. Applying (2) with r=2 yields

    d in {2,3,4}, n>=3 ==> kappa(K)<=n-3.            (4)

This holds with arbitrary nonlinear source/target coordinate changes,
without any bound on their degrees or word lengths. It is a lower-
dimensional spectrum obstruction, not an all-degree inversion theorem.

### Exact stabilization consequence for a three-dimensional seed

For ANY three-variable Keller map F of generic degree3 or4, and EVERY
s>=0,

    kappa(F times id_(A^s))=s.                     (5)

The literal identity factors give the lower bound s. Equation(4), with
n=3+s, gives the upper bound s. Thus arbitrary nonlinear recoordination
of any identity stabilization cannot reveal even one additional common
coordinate output. No stable coordinate-cancellation theorem is assumed.
The KNOWN degree3 map from KO2 is an actual instance; no existence of a
new degree4 seed is needed for the conditional degree4 statement.

For general maps only kappa(K times id_(A^s))>=kappa(K)+s follows directly.
We do NOT assert equality in all degrees or all dimensions.

## 2. LK1: the complete polynomial-map orbit criterion for HC4 gradients

Call two maps polynomially equivalent when one is obtained from the
other by an arbitrary polynomial automorphism in the source and an
arbitrary polynomial automorphism in the target, independently.

For a four-variable polynomial Keller map T, the following are equivalent:

    T is polynomially equivalent to a unit-Hessian gradient;
    T is polynomially equivalent to K times id_(A2)
       for some complex polynomial plane Keller map K.             (6)

The forward implication uses retained HCTANG2. If the gradient is an
automorphism, it is equivalent to id_(A4). Otherwise a CONSTANT LINEAR
source change puts its potential into

    f(z,w)=D(z)+w dot K(z),  z,w in A2, det(DK)=j in C*.

Its gradient is

    (xi,eta)=(grad D+(DK)^T*w, K(z)).

The source change (z,w)->(z,xi) is an ACTUAL polynomial automorphism:
its inverse is w=(DK)^(-T)*(xi-grad D), whose entries are polynomial
because det(DK) is a field unit. Permuting the outputs then gives

    (z,xi)->(K(z),xi).                              (7)

Conversely f=w dot K has constant Hessian determinant j^2, and the same
source change gives (7). This proves (6), not merely an abstract
equivalence between the two conjectures.

An arbitrary Keller map with two common coordinate outputs need only be
a PARAMETER FAMILY of plane maps, as in (1). It has not been proved to
be a polynomial product. Such a family is NOT silently substituted for
the stronger product condition in (6). Nevertheless if its d>1, (2)
supplies some actual noninvertible plane specialization and hence an
HC4 counterexample by the usual cotangent construction. Existence of
such a family is not supplied here.

Every nonautomorphic HC4 gradient has EXACTLY

    kappa(grad f)=2.                                (8)

Indeed (7) gives at least2, while (3) gives at most2. For automorphic
gradients kappa=4. In particular kappa=0,1,3 cannot occur for an HC4
gradient. The structural conclusion in this paragraph uses HCTANG2;
the dimension-descent theorem in Section1 does not.

### The entire marked field and normalization must descend as a product

Let B0=C[K_1,K_2], R0=C[z_1,z_2], F0=Frac(B0), L0=Frac(R0), and S0
the full normalization of B0 in L0. For (7), the embedded field pair is

    F0(xi_1,xi_2) subset L0(xi_1,xi_2).             (9)

This preserves generic degree and the WHOLE normal-closure permutation
action, including its full point stabilizer. If M0 is the normal closure
of L0/F0, its normal closure after adjoining the independent xi is exactly
M0(xi); a purely transcendental extension is linearly disjoint from the
finite field extensions. The subgroup actions and all intermediate
fields therefore extend without a new algebraic relation.

The full finite normalization and marked source are literally

    S=S0[xi_1,xi_2],  R=R0[xi_1,xi_2].              (10)

The first ring is finite over B0[xi], normal and has the required fraction
field, so it is the full normalization. The omitted divisors, branch
divisor and nonproperness locus are products with A2. After the polynomial
equivalences in (6), these are products in the corresponding target
coordinates. An unmarked rational field descent, an abstract product
normalization, or a standard coframe alone is not the marked product
assertion (10), nor a sufficient coordinate certificate.

This is an ALL-DEGREE necessary test for every proposed HC4 map conversion.
The plane exclusions apply to the FULL residual degree. A proper degree3
or4 intermediate field does NOT automatically have a polynomial plane
source, so those exclusions must not be applied to it.

## 3. Retire the direct gamma map-conversion and global mixed-chart queues

LJ0/1 proves that EVERY gamma level has an actual A4 source and a
noninjective polynomial etale projection Phi_c of generic degree3.
These facts and their assurance limits are reused; no old Jacobian,
collision, source-coordinate or reduction checker is rerun.

Generic FIELD degree is unchanged by polynomial automorphisms of source
or target. Equations(4),(6)--(8), together with PCC3, prove for EVERY c:

    kappa(Phi_c)<=1;
    NO polynomial source/target recoordination of Phi_c
       is a four-variable unit-Hessian gradient.                  (11)

This excludes ALL mixed four-coordinate source changes for that MAP
conversion, not only the two Darboux position pairs excluded in LJ1.
It does NOT decide whether the reduced symplectic form Omega_c itself
is polynomially standard. Even a successful standardization would not
turn this degree3 map into an HC4 gradient by separate source/target
automorphisms. The analogous direct degree4 conversion is excluded by
KR0. We adopt no fresh degree5 theorem.

There is a further useful obstruction before seeking a global generating
function. After ANY source and target polynomial recoordination, consider
the graph projection to TWO target coordinate functions and TWO source
coordinate functions. If this 2+2 projection were a polynomial
isomorphism of the graph onto A4, the two target pullbacks would form
a joint source coordinate pair. That would give kappa(Phi_c)>=2,
contrary to (11).

Thus EVERY such global polynomial 2+2 mixed graph chart is excluded,
including the usual globally invertible mixed projections used for
generating functions after a Darboux change. Local, formal, rational,
or finite non-isomorphic projections are NOT excluded by this argument.

### Limits of the retirement

The result concerns a FIXED MAP under separate source/target polynomial
automorphisms, and the specified globally invertible graph charts. It
does not exclude a new nonlinear restriction, a source-target mixing of
the graph outside these charts, multiplication by a variable matrix,
or a genuinely degree-changing construction. Nonlinear pullback of a
POTENTIAL changes its gradient by a variable Jacobian transpose; that
is not a separate target automorphism of the old map.

Compositions of k gamma maps have generic degree3^k, since field degree
multiplies. Only k=1 is excluded by the retained degree3 theorem; no
exclusion of all k>=2 or all proper low-degree intermediate fields is
claimed. No new composition/search queue is proposed in lieu of a
whole-source existence argument.

## 4. Provenance and remaining objective

The basic coordinate restriction in BCR2 and cotangent fibre-degree
correspondence in GS0 are retained and credited. The proof adds the
generic specialization argument, arbitrary-coordinate maximum and exact
stabilization bound, whole polynomial-map orbit criterion, and its
application to the gamma family/global mixed charts. No novelty priority
for standard block or transcendental-base-change facts is claimed.

PCC3 already imports Orevkov's theorem. The retained author-hosted PDF
was re-inspected at the COMPLETE first page, including the generic
multiplicity definition and Theorem1.1; it was rendered and visually
read. Its ten-page proof is NOT newly independently reproduced here.
KR0's published unrestricted degree4 exclusion and HCTANG2's existing
written reproduction/assurance boundaries are unchanged. There is
NO newly adopted external theorem and NO registered checker replay.

The finite checker tests the new block/inverse algebra, degree
specialization transcriptions and scope/provenance. Universal assertions
are the written arguments above, not inferred from bounded examples.
Independent replay, formal verification and external review remain false.

This closes an entire direct counterexample-conversion route, not HC4.
The main open implication remains actual polynomial PLANE source to
ordinary trace-dual pairing regularity, or full marked-energy integrality.
Further Darboux calculations for the fixed degree3 gamma projection are
not an HC4 completion gate. Full HC4/JC2 remains OPEN.
