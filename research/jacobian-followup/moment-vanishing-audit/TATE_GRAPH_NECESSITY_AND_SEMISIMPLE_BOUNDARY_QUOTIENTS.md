# Tate graph necessity and semisimple boundary quotients

16 November 2026. Internal written proofs. Independent replay, formal
verification and external review are not claimed. Both COMPLETE exceptional
Mathieu--Zhao cases remain **OPEN**; review stays deferred.

Every eventually pure failure, on SL2/Gamma for ANY finite Gamma, requires
NONZERO top Tate boundary graph cohomology. Arbitrary boundary genus,
including nonsemisimple local genus monodromy, cannot cause a failure by
itself. The statement retains all original markers and every polynomial
recomposition; it does not assert that all exceptional fibres have this
vanishing graph part.

More generally, if the reduced fibre quotient of the boundary source is
globally semisimple, safety under eventual purity is equivalent to ONE
intrinsic local unit test at zero. Actual relations cannot introduce a
second local unit obstruction in that case. This is a GLOBAL semisimplicity
hypothesis on a quotient, not a local eigenvalue semisimplicity assumption.

See the [checker](verify_haar_graph_necessity.py),
[receipt](haar_graph_necessity_verification.json) and
[source audit](haar_graph_necessity_source_audit.json).

## 1. A quotient lemma with a prescribed surviving unit

Let 0 -> L -> E -> V -> 0 be finite-dimensional local systems over a
connected base, with L a one-dimensional trivial system. Assume V is
GLOBAL semisimple over C. For ANY quotient q:E -> N retaining L, put
K=ker q. Its projection into V is injective; write Kbar for its image.
Because V is semisimple, choose a global projection r:V -> Kbar.
The composite

    E -> V -> Kbar -> K                                (1)

is identity on K and zero on L. Hence

    E=K direct_sum E',  L subset E',  E' ~= N.           (2)

No Hodge compatibility of K is required. In particular a retraction of
L on E for ANY subgroup of monodromy restricts to one on E', and a
retraction on N composes with q. Therefore

    L splits on E for that subgroup iff it splits on N. (3)

The subgroup can be global monodromy, one local cyclic monodromy, or a
finite-index subgroup. Importantly V may have NONSEMISIMPLE local
monodromy; global semisimplicity is what supplies (1).

## 2. A local unit functional suffices on an actual quotient

Let M be the full actual scalar Cauchy FUNCTION module, with its literal
rational unit u=1, BEFORE quotienting by C[s]. If a differential quotient
M -> N retains u and admits a LOCAL flat functional ell with ell(u)=1
at s=0, then regular singularity gives a differential morphism

    N -> K0=convergent meromorphic germs at0.            (4)

Compose M -> N -> K0 -> PP0=delta0. This fixes rational negative powers
of u, annihilates polynomials, and is Weyl-linear after polynomial
quotient. Thus it is well-defined on the ACTUAL Borel image, including
all actual relations. Under EVENTUAL PURITY that image contains delta0,
the Borel transform of the polynomial moment submodule C[t]. The map is
identity there, so it retracts the polynomial unit. The actual entire
moment image is indecomposable, forcing EVERY original marked moment
function to be polynomial.

This is exactly the local-functional proof in the
[intrinsic safety theorem](INTRINSIC_SOURCE_SPLITTING_AND_FULL_HAAR_SAFETY.md),
using only its unit covector, rather than the stronger sufficient
hypothesis that all local eigenvalue1 blocks are semisimple. N need not
embed in scalar functions or carry a mixed Hodge structure. A local
functional on a SOURCE need not descend to an arbitrary actual quotient;
Lemma1 is the additional argument that will ensure this below.

## 3. The reduced intrinsic local obstruction

For a regular phase f on X=SL2/Gamma retain the full regular base and
canonical admissible boundary source E. Write

    0 -> L=Q(-2) -> E -> V -> 0,
    V=(H2(F_s)/W2 H2(F_s))(-1),
    P=GrW5 E,  T=GrW6 E=GrW4 H2(F_s)(-1).             (5)

Thus V has only original weights5 and6. P is pure genus; T is the finite
monodromy Tate graph. Let N be the actual boundary quotient of E. The
original scalar module M surjects N, and u is retained, as proved in the
[canonical boundary theorem](CANONICAL_BOUNDARY_REDUCTION_AND_LOWEST_WEIGHT_INVARIANTS.md).

At any parameter a, including infinity, dualize (5). Pick a lift e of the
ambient generator in E^dual and let A_a act on V^dual. Its local class is

    xi_boundary,a=[(T_a-I)e] in V^dual/(A_a-I)V^dual.    (6)

This class is independent of the chosen lift. It vanishes exactly when
there is a local fixed ambient lift, equivalently a local unit retraction
on E. On a regular parameter it is zero. This is the intrinsic class of
the REDUCED source E, not an arbitrary physical Haar lift, and not an
unproved identification with the intrinsic class of the full source B.

If V is GLOBAL semisimple, Lemma1 applies to E -> N. Combining (3),(4)
and necessity from the global intrinsic safety theorem gives

    under EVENTUAL PURITY and V GLOBAL semisimple,
    ALL original markers safe iff xi_boundary,0=0.     (7)

Necessity: safety gives a global source-unit retraction and hence a local
one. Sufficiency: (3) transfers the source local retraction to N, and (4)
gives the actual Borel retraction. This does not infer finite GLOBAL
monodromy, finite LOCAL monodromy of V, or a Hodge-compatible actual kernel.

## 4. No top graph means full all-marker safety

Assume T=0, equivalently GrW4 H2(F_s)=0. Then V=P is PURE of weight5
and globally semisimple, by the audited pure-variation semisimplicity
theorem. There is NO assumption that P has trivial rank or semisimple
local monodromy, or that its global fixed part vanishes.

The intrinsic Tate-factor theorem from the
[unipotent extension proof](UNIPOTENT_HAAR_EXTENSIONS_AND_GENERAL_RECOMPOSITION.md)
applies to this geometric sequence. After unipotent base change its local
connecting morphism is

    delta:ker N_V -> L(-1)=Q(-3),                       (8)

and the actual strict limit MHS puts ker N_V in weights<=5. The target
has weight6. Hence delta=0, so xi_boundary,a=0 at EVERY puncture, INCLUDING
INFINITY. No assertion M=W is used: P can have nonzero local nilpotent
monodromy and a nontrivial relative monodromy filtration. Descent from
the finite local cover follows by averaging a fixed lift over its finite
cyclic quotient. Therefore (7) proves the general theorem

    GrW4 H2(F_s)=0 + eventual purity
       implies ALL-original-marker safety.             (9)

For a smooth SNC compactification of the generic affine surface F_s,
the top group in (9) is the first homology of the boundary dual graph,
with its Tate twist. Its vanishing means that graph is a FOREST; boundary
curves may have arbitrary genus and compact-middle cohomology may have
arbitrary rank. Disconnected fibres and all their components remain.

Consequently EVERY eventual-pure failure on EITHER exceptional quotient,
indeed on SL2/Gamma for ANY finite Gamma, has T!=0. Boundary weight3/genus
alone is now excluded as a failure mechanism. This is stronger than the
earlier pure-weight2 fibre theorem and closes every top-graph-free phase,
not merely phases with no boundary genus.

By the full source-splitting theorem, (9) also implies GLOBAL unit splitting
and global semisimplicity of E on this eventual-pure locus. In the
[algebraic 1-motive model](ALGEBRAIC_BOUNDARY_ONE_MOTIVE_AND_CONSTANT_SECTION.md),
G_f is then an abelian scheme with no torus and its section sigma_f is of
CONSTANT TYPE up to isogeny. This conclusion is forced by eventual purity;
an arbitrary abelian normal function without purity is not asserted constant.

## 5. EVERY H and an exact finite bad-value obstruction

Let H be ANY nonconstant polynomial. The established source law gives

    V_H(f)=pi_*V_f.                                    (10)

Finite-index induction and restriction preserve and reflect global
semisimplicity. In particular the hypothesis for (7) reflects EVERY H,
BEFORE purity. Top Tate rank becomes deg(H)*rank T_f; T=0 is likewise
preserved and reflected. Thus if T_f=0, EVERY eventually-pure H(f) is
fully safe, without assuming f or any shift f-a eventually pure.
If its actual primitive Cauchy rank is>1, NO eventually-pure H(f) exists:
otherwise safety would make composite actual rank1 and contradict the
retained actual rank injection. Rank1 has its unique safe centre c and
exactly H(c)=0 gives eventually-pure, fully safe composites.

For the more general globally semisimple V case define the FINITE set

    Bad_boundary(f)={a in C:xi_boundary,a!=0}.           (11)

Only full primitive bifurcations can occur in (11). For each root a of
H of multiplicity e, cyclic induction gives the canonical coinvariant
isomorphism from its e branch copies to primitive coinvariants. On (6),
the ambient diagonal lift changes in ONE wrap-around branch. The
isomorphism therefore sends xi_boundary,0 of H(f) to the tuple

    (xi_boundary,a)_(H(a)=0),                          (12)

with factor1, no cross-root cancellation, and ALL branches retained.
It is the linear-algebra proof already audited for full extensions,
applied to the reduced source whose induction/SUM law is geometric.

Equations(7),(10),(12) give the all-degree criterion

    if V_f GLOBAL semisimple and H(f) eventually pure,
    H(f) fully safe iff H^-1(0) misses Bad_boundary(f).  (13)

The source-global criterion gives a further purity exclusion: if the
primitive source unit is NONSPLIT and H^-1(0) misses (11), then NO such
eventually-pure H(f) exists. Otherwise (13) would give safety and hence
primitive source splitting, a contradiction. If Bad_boundary(f) is empty
but the source unit is globally nonsplit, NO eventually-pure composite
exists in ANY degree. Neither hitting a bad value nor global nonsplitting
is a producer of eventually pure phases. All implications retain that
quantifier; the primitive may be nonpure.

## 6. Exhaustive alternatives for every remaining pure failure

For ANY eventually-pure failing phase, T!=0 by (9). Moreover exactly the
following exhaustive alternatives are left:

* xi_boundary,0!=0: an intrinsic local SOURCE obstruction is present.
* xi_boundary,0=0: V is NOT globally semisimple, by (7).

In the second alternative the pure grades P and T are individually
globally semisimple, so their extension in V is genuinely nonsplit.
Both P and T are nonzero: a pure-only V would be semisimple. Thus a
failure invisible to the reduced intrinsic local unit test requires a
genuine GLOBAL genus--Tate coupling, in addition to the moving section
sigma_f. This keeps all nonsemisimple cases in scope. It does not prove
the coupling or the toric obstruction universally absent.

Without global semisimplicity the quotient lemma cannot be used. Even
admissible LOCAL limit data can have a split source unit but a quotient
with its unit in the monodromy image. A useful control has basis u,a,b,t,
original weights4,5,5,6, N b=a, N t=a, N u=N a=0. Its relative limit
weights are4 on u,a and6 on b,t. Take F3=span(b,t), F2=all; this is a
split Tate limit MHS and the original weight5 genus block is a standard
polarizable rank2 nilpotent orbit. The non-Hodge invariant kernel
span(a-u,b-t) retains u in the quotient but identifies N t with u.
This is LOCAL data only, not a supplied global algebraic phase or an
actual pure counterexample. It shows why a source fixed lift cannot
generally replace the actual quotient and why M=W would be wrong here.

The full target remains CONSTANT TYPE of the canonical algebraic section
for EVERY exceptional eventually pure phase, or an actual admissible
nonconstant-type pure failure with a fixed original bad marker. The new
theorem proves it whenever the top graph vanishes and confines all
remaining failures to the exhaustive alternatives above. It does NOT
assert that the full exceptional locus is top-graph-free: known actual
pure safe examples have unbounded nonzero Tate graph rank.

All universal arguments above are written proofs using the retained
geometric source/limit and actual Borel results. Finite linear-algebra
controls check their transcription, not universal machine verification.
No degree/support/prefix campaign, arbitrary-ring topological promotion
or external-review prerequisite is introduced. Both COMPLETE exceptional
MZ cases remain OPEN; review deferred.
