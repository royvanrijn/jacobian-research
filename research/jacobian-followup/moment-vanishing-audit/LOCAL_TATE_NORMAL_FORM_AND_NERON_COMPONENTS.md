# Local Tate normal form and the Neron component criterion

16 November 2026. Internal written proofs, with named primary imports.
Independent replay, formal verification and external review are not claimed.
Both COMPLETE exceptional Mathieu--Zhao cases remain **OPEN**; review deferred.

For EVERY phase, the reduced intrinsic local obstruction is an ordinary
rational component class of its algebraic semiabelian section. It vanishes
exactly when that section specializes to a TORSION component of the Neron
lft model. Equivalently, local monodromy of the canonical boundary source
has no extra size-two Jordan block beyond its genus blocks. Under eventual
purity, either test is necessary and sufficient for ALL original markers.

The invariant-image qualification in the previous Tate-factor theorem can
also be removed for the REDUCED source: every invariant top Tate class
lifts to a local invariant of its genus--Tate quotient. This statement
does not concern the full source with compact weight-four terms.

See the [checker](verify_haar_neron.py),
[receipt](haar_neron_verification.json) and
[source audit](haar_neron_source_audit.json).

## 1. The genus--Tate quotient splits locally

Use the canonical admissible boundary source for any nonconstant regular
phase f on X=SL2/Gamma, ANY finite Gamma:

    0 -> L=Q(-2) -> E -> V -> 0,
    A=W5 E, P=GrW5 E, T=GrW6 E.                       (1)

The [canonical reduction](CANONICAL_BOUNDARY_REDUCTION_AND_LOWEST_WEIGHT_INVARIANTS.md)
has already removed the compact pure weight-four quotient: W4 E=L.
P is pure polarizable of weight5 with types(3,2),(2,3), and T is a finite
Tate local system of weight6. At EVERY puncture, including infinity, take
a finite local cover killing all semisimple monodromy, and write N for
the unipotent logarithm. The established strict limit MHS satisfies

    M3 E=0, M6 E=E, N(E) subset M4 E, N^2=0.          (2)

Its actual Hodge quotients and subobjects have strict limit filtrations.
For the pure genus part the monodromy filtration is centred at5, so

    M4 P=image N_P.                                   (3)

Indeed its Jordan blocks have size1 at grade5 or size2 at grades4,6.
No local semisimplicity of the genus part is imposed.

In 0 -> P -> V -> T -> 0 take any t in T and a lift v in M6 V=V.
Since N_T=0, Nv lies in P; by (2) and strictness it lies in M4 P.
Equation(3) supplies p in P with Np=Nv. Then v-p is invariant. Thus

    ker N_V -> T is SURJECTIVE.                       (4)

Choose such lifts linearly. This splits V=P direct_sum T as a LOCAL
SYSTEM on the cover. It need not split as a limit MHS, and says nothing
about GLOBAL splitting of V. A genuine global genus--Tate extension is
still allowed. Descending (4), average over the finite cyclic local
semisimple action on its original eigenvalue-one sector to obtain

    ker(T_V-I) -> T^1 is SURJECTIVE.                   (5)

Here T_V denotes original local monodromy, and T^1 is the original
eigenvalue-one part of the top Tate space, not the entire Tate space
when its finite monodromy is nontrivial. Therefore the former invariant
IMAGE J_boundary is exactly T^1. This strengthens the conservative
image-only statement in the
[full local theorem](FULL_INTRINSIC_LOCAL_TATE_CRITERION.md).
The image can still be proper for the UNREDUCED fibre quotient
H2(F)(-1), which retains a compact pure weight-four part. No result or
countercontrol about that full quotient is discarded.

## 2. Exactly one extra local Jordan block is the entire obstruction

The subextension 0 -> L -> A -> P -> 0 splits as a local system.
To see this, its connecting map ker N_P -> L(-1) is a morphism of limit
MHS. The source has weights at most5 by (3), whereas L(-1)=Q(-3) is
pure of weight6, so the map is zero. For a one-generator unipotent
representation this is precisely unit splitting: the extension row
annihilates ker N_P and is therefore a row multiple of N_P, removable
by a change of lift. Descent to original monodromy uses finite averaging.

Combine this splitting with (4). Lift the invariant T-summand of V to E.
Its logarithm lands in L. In a suitable rational local-system basis,

    E = P direct_sum E_T,
    N_E = N_P direct_sum [[0,beta],[0,0]],
    0 -> L -> E_T -> T -> 0, beta:T -> L.             (6)

The finite original semisimple action commutes with N. Splittings in
this construction may be averaged to respect it, so beta is zero on
all original non-one Tate eigenspaces. Its restriction to T^1 is the
full connecting row, not a functional on an unknown proper image.
Equation(6) is a decomposition of underlying local representations;
neither a Hodge splitting nor a global decomposition is asserted.

The unit lies in image N_E if and only if beta is nonzero. Hence

    rank N_E - rank N_P = rank beta in {0,1},
    xi_boundary,a=0 iff rank N_E=rank N_P.             (7)

In particular, a nonzero reduced local class means exactly ONE EXTRA
size-two Jordan block, however large the genus and top graph ranks.
For ANY positive integer e killing semisimple local monodromy, (2) gives

    N_E^(e)=log(T_E^e)=T_E^e-I,
    N_P^(e)=log(T_P^e)=T_P^e-I.                        (8)

These logarithms rescale with the cover; their ranks do not. Thus
ordinary rational matrix ranks after a finite cover decide (7).
Original non-one genus sectors contribute equally to the two ranks;
they need not be deleted to perform the test. No actual Haar orbit,
period relation kernel, or semisimplicity of V is required.

By the full intrinsic local theorem, under EVENTUAL PURITY,

    ALL original markers are safe
      iff rank(T_E^e-I)=rank(T_P^e-I) at a=0.           (9)

Without purity, (7) is still an exact source-local statement, while
(9) is not asserted as a safety criterion.

## 3. The local algebraic section detects exactly the same class

The [geometric 1-motive construction](ALGEBRAIC_BOUNDARY_ONE_MOTIVE_AND_CONSTANT_SECTION.md)
supplies, on an open algebraic parameter base and up to isogeny,

    M_f=[Z -> G_f], 1 |-> sigma_f,
    T_H(M_f)=(E(2))^dual.                             (10)

Here G_f is a semiabelian scheme, retaining its entire abelian genus
part, torus part, and extension between them. At a finite puncture use
z=s-a; at infinity use z=1/s. Put K=C((z)), R=C[[z]], and take the
Neron lft model mathcal G/R of the generic semiabelian variety G/K.
The standard Neron existence and component finiteness theorem gives

    Phi_G=pi_0(mathcal G_C), a finitely generated
    abelian group.                                   (11)

Precise primary references for this named input are Liu--Lorenzini,
J.reine angew.Math.532(2001), introduction and Section1.2, and Hertgen,
arXiv:1509.02797v1, introduction; both refer to BLR Chapter10 for
Neron lft existence. We import the standard finiteness theorem, not
the flawed intermediate lemma mentioned in Liu--Lorenzini Remark1.3.
Its cited correction and the later independent primary statement are
recorded in the source audit. No theorem about wild ramification is used.

The Neron property extends sigma to an R-point. Define its rational
component class by

    nu_a(sigma)=[sigma_C] in Phi_G tensor Q.           (12)

This is meaningful for the Neron LFT model; replacing it by a finite
type torus model can lose precisely the components that matter.
The claimed equivalence, BEFORE any purity hypothesis, is

    xi_boundary,a=0 iff nu_a(sigma_f)=0
      iff [m]sigma_C belongs to mathcal G_C^0
          for SOME nonzero integer m.                 (13)

It is torsion of the COMPONENT, not torsion of the section, and not
necessarily membership of sigma_C itself in the identity component.

Here is a Kummer proof covering arbitrary semiabelian degeneration.
Smoothness, completeness and the algebraically closed residue field
give an exact sequence of abstract abelian groups

    0 -> mathcal G^0(R) -> G(K) -> Phi_G -> 0.         (14)

For each integer n>0, mathcal G^0(R) is n-divisible. Its special fibre
is a connected smooth commutative group over C. Chevalley decomposition
reduces surjectivity of multiplication by n on this fibre to abelian
varieties, tori and unipotent groups, where it holds in characteristic0.
Given x in mathcal G^0(R), choose an n-th root of its reduction. The
differential of [n] at that root is multiplication by n on the Lie
algebra and is invertible over R. The henselian lifting property for
this etale map lifts the root to an R-point. No semistable reduction,
split torus, or properness of G is assumed.

Consequently projection induces, for every n,

    G(K)/nG(K) = Phi_G/nPhi_G.                         (15)

For injectivity, if a component is n times another, lift the latter to
G(K); the remaining point lies in mathcal G^0(R), where division by n
is possible. This proves (15), not just surjectivity.

For a fixed prime ell, the Kummer sequence injects the group in (15)
for n=ell^j into H1(K,G[ell^j]). The compatible injections give

    completed_ell Phi_G -> H1(K,T_ell G).             (16)

Indeed the finite groups H0(K,G[ell^j]) satisfy Mittag--Leffler, so
there is no inverse-limit H0 obstruction in (16). Since Phi_G is
finitely generated, tensoring with Q_ell shows that the rational
Kummer class of sigma vanishes EXACTLY when its component is torsion.
Finite component torsion is allowed and disappears rationally.

For completeness, this Kummer class is the same LOCAL obstruction as
in (10), rather than a different arithmetic invariant. The analytic
Betti realization of [Z -> G] is the inverse image of Z*sigma under
the fibrewise exponential Lie(G) -> G. A choice of logarithm of sigma
gives a lift of1; continuation around z=0 differs by an element of
H1(G,Z). Modulo ell^j this is exactly the torsor of ell^j-division
points of sigma, the Kummer cocycle in (16). Finite division covers
have the same local analytic and formal monodromy: over C((z)) they
are obtained by finite Puiseux extensions. Thus Betti and etale local
classes agree under the comparison of these finite covers.
Tensoring the Betti class with Q_ell is faithful on its finite
dimensional rational coinvariant space. Duality and the twist in
(10) identify its zero class with a local retraction onto L in E,
which is exactly xi_boundary,a=0. This proves (13).

Isogeny, integer rescaling used to construct (10), and finite local
base change preserve AND reflect zero in (13), either by rational
realization or restriction/corestriction. We do not claim a universal
coordinate scaling formula between arbitrary Neron component models.

## 4. Full safety, exact algebraic data, and every recomposition

Combining (9),(13) with the full local theorem gives, for EVERY
eventually pure phase on EVERY finite SL2 quotient,

    full ALL-original-marker safety
      iff sigma_f has a TORSION Neron component at0
      iff rank(T_E^e-I)=rank(T_P^e-I) at0
      iff xi_boundary,0=0.                            (17)

Under purity these are also equivalent to the previously proved
GLOBAL constant-TYPE section condition. The implication from one
local torsion component to that global condition USES purity; it is
not a general local-to-global theorem for semiabelian sections.

If a component group presentation Phi_G=Z^r/image C and a representative
v of the section component are supplied, then

    nu_a=0 iff v belongs to image(C tensor Q)
           iff rank_Q[C|v]=rank_Q C.                  (18)

This is finite exact integer/rational linear algebra, equivalently
Smith normal form. Formula(18) does not implement a general producer
of Neron models or section specializations. Finite torsion relations
must be retained before rationalizing.

For a split torus G=Gm^t this reduces to the FULL vector
(ord_a u_1,...,ord_a u_t)=0. A scalar norm can cancel opposite entries
and is insufficient. For a genuinely nonsplit semiabelian extension,
the actual component group and its relations are required: arbitrary
torus coordinates or independent norms need not encode (12).
If G is abelian, Phi_G is finite and (13) is automatically zero;
this recovers graph-free safety with arbitrary degenerating genus.
The constant section2 in Gm has zero component and is allowed even
though it is nontorsion and its Hodge extension is nonsplit.

Let Bad_Neron(f) be the finite set of finite primitive parameters with
nonzero rational component. Equation(13) identifies it exactly with
Bad_boundary(f). The retained all-branch theorem gives, for EVERY
nonconstant polynomial H of EVERY degree,

    IF H(f) is eventually pure, it is fully safe
      iff EVERY a with H(a)=0 has nu_a(sigma_f)=0.     (19)

Multiplicity and ramification are already accounted for by the
factor1 cyclic-coinvariant law; every primitive component is retained.
No primitive purity is inferred and no scalar norm replaces the tuple.
In particular the actual nonexceptional Long composite F(F+1) remains
a pure failure: the toric section s/(s+1) has nonzero valuations at
0 and -1, although the two recomposed branch values have scalar norm1.
No exceptional counterexample is produced by this control.

The whole exceptional problem is now equivalently: EVERY actual
eventually pure exceptional phase has TORSION section component at0,
or realize a nonzero component with an actual invariant pure phase
and an original fixed bad marker. Arbitrary top graph rank, degenerating
genus and globally nonsplit genus--Tate coupling are allowed. No such
universal torsion theorem, pure producer, coefficient-ring promotion,
universal machine proof or external-review gate is supplied here.
Both COMPLETE exceptional MZ cases remain OPEN; review stays deferred.
