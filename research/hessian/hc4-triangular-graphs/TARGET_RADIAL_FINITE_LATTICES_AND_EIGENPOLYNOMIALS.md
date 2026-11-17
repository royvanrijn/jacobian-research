# The lifted target radial field has no extra finite lattice or eigenpolynomial

16 November2026. HC4MYME0. ACTUAL complex polynomial PLANE Keller maps,
ALL DEGREES, every target centre, every polynomial target coordinate frame.
This is a written deduction from retained plane boundary and finite-cover
inputs, not a new external theorem import. The original HC4/JC2 and full
trace-pairing/fullmarkedenergy goal remain OPEN. External review DEFERRED.

## 1. Statement and distinction from source contraction

Let K=(P,Q):A2->A2 have nonzero constant Jacobian, B=C[P,Q] subset
R=C[x,y], F=Frac(B), L=Frac(R). For ANY b=(p0,q0) set

    rho_b=(P-p0) partial_P+(Q-q0) partial_Q,
    E_b=(P-p0) D_P+(Q-q0) D_Q.

The lifted derivatives are the polynomial source derivatives from the
Keller coframe. In particular E_b(R) subset R and E_b(B) subset B.
For EVERY finite B-submodule W subset R,

    E_b(W) subset W  ==>  W subset B.                    (1)

No freeness, projectivity, containment of1 or algebra structure on W is
assumed. Hence the largest finite B-module preserved by E_b is B itself.
Its locally finite part in the source is EXACTLY

    R_lf(E_b)={r:dim_C span(E_b^j r:j>=0)<infinity}=B.    (2)

All polynomial eigenfunctions of E_b therefore descend to B. Its
nonzero polynomial eigenvalues are nonnegative integers; an eigenfunction
with eigenvalue k is precisely a homogeneous degree-k target polynomial
in P-p0,Q-q0. In particular R^(E_b)=C, even though E_b has the rational
first integral (P-p0)/(Q-q0).

The [LN0/1 theorem](CONTRACTING_SOURCE_GRADINGS_AND_FINITE_LATTICES.md)
concerns a contracting SOURCE Euler field, which is locally finite on ALL
R and need not preserve B. Here the TARGET Euler field preserves B and
is lifted to R, where local finiteness is not assumed. The distinction is
essential. [LX1](NAMBU_SOBOLEV_LATTICE_AND_FINITE_DERIVATIVE_CORE.md) already
gives (1) for target additive fields by a one-variable Puiseux argument.
The new argument uses trace discriminants to treat the radial field.

## 2. Every omitted image is transverse to the target radial field

Use the full finite normalization S of B in L and its actual open source
Spec R. The omitted divisors and all retained source divisors together
account for all height-one points. The retained KU1 gives S subset R;
the source is an open subset of Spec S by Zariski Main.

The PREVIOUSLY IMPORTED Chau component Remark, recorded in JY0 and audited
in LX Section1, says that no irreducible nonproper-image component is
isomorphic to A1. This is the component Remark, not the theorem about
the whole exceptional set. No stronger reading is imported here.

An irreducible polynomial f defining a rho_b-invariant target curve must
satisfy f divides rho_b(f). Total degree in the shifted target coordinates
then gives rho_b(f)=c f for a scalar c. The homogeneous decomposition
makes f homogeneous. Over C an irreducible homogeneous polynomial in
two variables has degree one. Thus every invariant irreducible curve is
a LINE through b, hence an A1. No omitted-image component is invariant.
Therefore rho_b is generically transverse to EVERY omitted image,
including every branch-image component of S/B.

This works after ANY polynomial target automorphism: the embedded rings
B,S,R are unchanged, and an affine line in that coordinate frame is an
A1 component in the intrinsic sense excluded by the retained Remark.
No missing-value hypothesis is needed in (1).

## 3. Bounded derivative poles put a stable finite module inside S

Suppose W is finite over B and stable under E_b. At any omitted prime
E its elements have a common valuation lower bound: take the minimum
over finitely many B-generators, since B has nonnegative valuations there.

If w in W has negative valuation -q, choose a general point of E and
its target image. Both surfaces and the image are smooth there. After
a tame local change, target coordinates (t,z) and source coordinates
(u,v) give t=u^e,z=v. Generic transversality gives

    E_b=a(u^e,v)/(e*u^(e-1)) partial_u
                       +c(u^e,v) partial_v,   a(0,v)!=0.

The leading term of E_b^j(w) has exponent -q-j*e and coefficient

    leading(w)*a(0,v)^j*product_(i=0..j-1)(-q-i*e)/e^j. (3)

Every coefficient is nonzero in characteristic zero. Tangential
derivatives and derivatives of coefficient functions have larger
u-order, so cannot cancel this leading term. This includes e=1.
The valuations tend to minus infinity, contradicting the common bound.
Consequently every w in W is regular at EVERY omitted divisor. It is
already regular on the actual open source. Normality gives W subset S,
including at finite singular points. This is the valuation method of LX,
used here with the radial field and all boundary primes.

Thus A=B[W] is a finite B-algebra contained in S. It is a domain,
and E_b(A) subset A by the product rule. Put E=Frac(A); it is an embedded
intermediate field F subset E subset L. The next argument does NOT assume
that E or its normalization has a polynomial-plane model.

## 4. A stable finite order has logarithmically invariant discriminant

Fix an irreducible target divisor D, write O=B_(D), and choose a
uniformizer t. The derivation rho_b maps O to O. A_O is finite torsion-free
over this DVR, so is free, of some rank m, with basis a_1,...,a_m.
Stability gives a connection matrix C over O for E_b(a_i).

Let G=(Tr_(E/F)(a_i*a_j)). Its entries are in O because A is integral
over the normal base. Its determinant Delta is nonzero by separability.
Field trace commutes with the lifted derivation, so

    rho_b(G)=C^T G+G C,
    rho_b(Delta)/Delta=2 tr(C) in O.                    (4)

The second equality is Jacobi's determinant identity over F; it does
not require a unit discriminant. If v_D(Delta)=s>0 and rho_b(t) is a
unit, write Delta=t^s*h with h a unit of O. Then

    rho_b(Delta)/Delta=s*rho_b(t)/t+rho_b(h)/h

has a pole. The second summand is regular; the first has a nonzero
coefficient. This contradicts (4). Hence

    D transverse to rho_b ==> v_D(Delta)=0.             (5)

The normal finite cover in E has branch support contained in this
order discriminant. For example, trace duality over the DVR shows
that a unit trace determinant makes the order finite etale, so its
normalization cannot introduce ramification there. Equivalently use
the retained different/discriminant identities. This is a statement
about the WHOLE finite order, not only the diagonal of source base change.

Branch support of an embedded intermediate also lies in the original
branch support, by the valuation tower. Section2 makes every original
branch component transverse to rho_b; (5) removes ALL such components.
The intermediate therefore has no codimension-one branch. Retained
purity makes its finite normalization etale over A2, and connected
finite-etale triviality over complex A2 gives E=F. Since W subset S is
integral over B and now lies in F, normality of B gives W subset B.
This proves (1), in ALL degrees and without intermediate-plane recognition.

No normalization-stability theorem for derivations was needed: (4) is
proved on the stable order itself. In particular we do not extend a base
derivation regularly to an arbitrary finite normalization by fiat.

## 5. Locally finite elements, controls, and the remaining Keller mark

If r has a finite-dimensional C-linear E_b-orbit V, the B-module B*V
is finite and stable because E_b(B) subset B. Equation(1) puts r in B.
Conversely ordinary target Euler is locally finite on B. This proves (2)
and the eigenfunction assertions by its standard homogeneous decomposition.
The proof does not infer local finiteness on the source from two target
eigenfunctions or from their positive weights.

The WHOLE-plane hypothesis matters. Reuse the retained ALL-e Laurent
constant-Jacobian controls, without replaying their checkers:

    P=u^e, Q=v/u^(e-1), R_e=C[u,u^-1,v], e>=2.

For target centre zero the lifted radial field is

    E=(u/e) partial_u+((2e-1)/e)*v partial_v.

EVERY Laurent polynomial has a finite-dimensional orbit. Thus the locally
finite part is all R_e, and the finite normalization C[u,Q] is stable
and strictly larger than B. Its branch line P=0 is invariant, violating
Section2; the source has nonconstant units. This is NOT a whole-plane
Keller or HC4 counterexample. Fractional source eigenvalues occur here,
exactly where the actual-plane theorem disallows them.

For an actual omitted value b, [MD0](SMOOTH_PROJECTIVE_PENCILS_WITH_ARBITRARY_GENUS.md)
gives the zero-free E_b with divergence2, a smooth surjective projective
pencil, and rational first integral. The result here adds (2), but does
not exclude such a field. Its finite normalization and trace-dual are
not proved invariant. Indeed a transverse ramified divisor gives a
pole of order e-1 in E_b(u); the radial source field cannot extend
regularly across the full ramified boundary under failure. Completeness
and local finiteness are not automatic; KX/KZ already give stronger
inversion consequences IF those hypotheses are supplied.

The trace identity (4) is not an automatic argument for the full S or
Lambda: its regular matrix C exists ONLY for a supplied stable finite
order. It cannot be assumed to establish that very stability. The
original implication ACTUAL polynomial plane source ==> Tr(Lambda*Lambda)
subset B, or full marked-energy integrality, is still OPEN.

## 6. What the new binary-affine geometry does and does not add

This is a DIRECT APPLICATION of retained PF123 and PF59, not a new
nodal-deformation theorem. For a binary affine sheet-block quotient,
PF123 makes the full intermediate branch curve C A1-normalized with
either an everywhere immersive normalization or exactly one ordinary
cusp. It is noncontractible. Every double value has exactly two smooth
branches, arbitrary contact, and EMPTY original source fiber.

In the IMMERSIVE case, every local intermediate monodromy image is
abelian: the finite factors above different local branches have disjoint
supports and cyclic inertia. Therefore O(C) surjects onto the binary
affine group Q. This group is nonabelian: a faithful transitive abelian
group is regular, whereas its nonidentity transvections fix sheets.
PF59's general curve theorem now applies, with no intermediate A2
assumption. On the fixed efficient target infinity compactification,

    eta=deg N_(P1->X)=-K_X.Cbar-2 <=0.                  (6)

Positive eta is thus excluded in this profile for EVERY block degree
and block size, even when the original cover is imprimitive or has
visible original indices larger than2. This transfers an existing
all-curve theorem using newly established hypotheses; it is not a new
quartic campaign. It does not exclude eta<=0.

In the ONE-CUSP case the intermediate local image at the cusp is S3,
acting on b/4 triples and fixing b/4 sheets. Its commutator C3 is
nontrivial. Monodromy therefore does NOT factor through O(C), which
kills every local commutator. One cannot silently transfer (6) by
smoothing that cusp and keeping only the Orevkov invariant.

The primary Neumann--Norbury Definition1 and Theorem1, and Orevkov's
explicit negativity hypothesis, were checked in primary web text.
Orevkov's author PDF also explicitly withdraws an old general Jacobian
reduction; it is not an input here. The retained PF59 imports and sealed
source files are reused, not independently re-audited or replayed.
No unconditional abelianity theorem is inferred from polynomial
parametrization, immersivity, or preservation of a link alone.

The binary affine exclusion is still unproved. Natural alternating and
symmetric quotients and the full finite-algebra compatibility also remain
OPEN. New symbolic controls transcribe (3)--(4), field/coframe identities
and the Laurent distinction; universal geometry follows from the written
proof, not a bounded calculation. No independent replay, formal
verification, external review or full goal completion is claimed.
