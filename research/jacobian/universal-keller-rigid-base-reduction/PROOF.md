# Rigid bases, the quartic automorphism obstruction, and Sylow degree budgets

15 November 2026. Characteristic zero unless explicitly stated otherwise.

These are written proofs with small exact controls, not formal verification or an independent review. The construction inputs and prior fixed-locus results are pinned in `sources.json`. No claim of publication priority is made. The original presentation-free affine Keller descent problem remains open.

## Conventions

A canonical natural symmetric datum of rank N consists of a dominant etale self-map, compatible source and target S_N actions, a fixed target point, and its **complete** reduced N-point fibre with the natural faithful S_N action. The generic degree M may exceed N. All geometric statements are over an algebraically closed field. Results over other characteristic-zero fields follow where explicitly stated by geometric descent. A geometric torus quotient is a smooth variety having a finite etale torus cover; this is the published notion of a Q-algebraic torus.

The previous quartic model (UAQUARTSELF1) is W=T/S_3, where

    T=Gm_x x Gm_z x A,
    A={(v1,v2,v3) in Gm^3 : v1 v2 v3=1},
    r(x,z;v)=(theta*x,z;v3,v1,v2),
    s(x,z;v)=(x^-1,-z;v2,v1,v3),  theta^3=1, theta!=1.

Over the original ground field the first Gm is the cyclotomic norm-one torus. Geometrically the displayed split coordinates apply. The free action has faithful linear part on the character lattice. Its fixed universal degree-four volume-preserving self-map Phi is induced by (x,z;v)->(x,z/2;v1^2,v2^2,v3^2). Its geometric monodromy S_4 is primitive, so UAFINATOM1 proves strong absolute stable atomicity. The results below concern the additional canonical action, not this established universality or atomicity.

## 1. UA1RIGIDCYL1: affine coordinates cannot move an A1-rigid base

Let Z be a smooth geometrically integral affine variety such that every geometric morphism A1->Z is constant. Then every morphism U x A^r -> Z, for a reduced variety U, factors uniquely through U. Indeed, its restriction to any affine line in a geometric affine fibre is constant, hence its restriction to the entire fibre is constant. It agrees on all geometric closed points with its restriction to the zero section composed with projection. Reducedness and the separated target give equality of morphisms. Descent gives the corresponding statement over the ground field.

In particular every endomorphism of Z x A^r has form

    G(z,v)=(f0(z),g(z,v)).                                      (1.1)

An automorphism induces an automorphism of Z: apply (1.1) to it and its inverse. Thus projection defines a homomorphism Aut(Z x A^r)->Aut(Z). If two A1-rigid bases have isomorphic affine cylinders, the induced base maps are inverse; fibre dimension then also identifies the numbers of affine coordinates. This is the standard relative A1-rigidity condition of Asok--Ostvaer, Definition4.1.6, proved here directly under the geometric hypothesis.

If G in (1.1) is dominant etale, its differential is block triangular. Both f0 and the vertical map

    H(z,v)=(z,g(z,v))

are etale, and G=(f0 x id) o H. If G is **finite**, UAFINATOM1 (finite composite between normal affine varieties) makes both factors finite. Each geometric vertical fibre of H is then a finite etale self-map of A^r. In characteristic zero A^r has no nontrivial connected finite etale cover; hence these fibres are automorphisms. The finite etale map H has degree one and is an automorphism. Consequently

    every finite etale G is (f0 x id) after a source automorphism. (1.2)

No assumption that the coefficients in g are constant was used. This closes arbitrary finite coefficient coupling on any such cylinder.

Every geometric torus quotient in characteristic zero satisfies the hypothesis. Pull back a finite etale torus cover along A1->Z. The pullback splits into copies of A1, giving a lift to the torus. All invertible regular functions on A1 are constant, so this lift, and hence the original map, is constant. Characteristic zero matters: A1 has wild finite etale covers in positive characteristic. The pseudo-planes in UAPSEUDOSYM1 have affine-line curves; this theorem does not settle their remaining projection-mixing cubic case.

## 2. UATORUSDOM1: dominant maps of torus quotients are automatically finite etale

Let X and Y be smooth geometric torus quotients of the same dimension in characteristic zero. Every dominant regular map f:X->Y is finite etale. Moreover its map on holonomy embeds the holonomy group of X into that of Y. In the self-map case the embedding is onto. Thus UAINFRASELF1 needs no separate finiteness assumption for dominant maps of the **base** itself; the statement about cylinders still distinguishes finite and nonfinite maps.

Here is a proof, including the characteristic-cover step used below. Work over C by spreading out and base change. Write

    1 -> Lambda_X -> Gamma_X -> F_X -> 1,
    1 -> Lambda_Y -> Gamma_Y -> F_Y -> 1,

for the universal affine deck groups, their translation lattices and linear holonomy groups. Each Lambda is Z^n and is the unique maximal abelian normal subgroup (UAINFRASELF1). The same description holds for a finite-index subgroup H of Gamma_Y, with maximal abelian normal subgroup H intersect Lambda_Y. One can see the uniqueness directly: for an element of an abelian normal subgroup with linear part L, commutators with lattice translations give (L-I)^2=0; finite-order semisimplicity forces L=I.

The image H=f_*(Gamma_X) has finite index in Gamma_Y. To justify this without assuming f finite, restrict to dense opens on which the dominant generically finite map is finite etale. The inclusion of each dense open in its smooth ambient variety induces a surjection on fundamental groups. The finite covering index then bounds the index of H.

The image f_*(Lambda_X) is abelian normal in H and therefore lies in H intersect Lambda_Y. The quotient H/f_*(Lambda_X) is finite, being a quotient of F_X. It follows that f_*(Lambda_X) has rank n and that f_* is injective on Lambda_X. The kernel of f_* consequently injects into finite F_X. The affine universal deck group Gamma_X is torsion free: a finite-order affine transformation of C^n has a fixed point, whereas a deck transformation acts freely. Thus the entire kernel is trivial. Since Gamma_X maps isomorphically onto H, uniqueness of the maximal abelian normal subgroup gives

    f_*(Lambda_X)=H intersect Lambda_Y,
    f_*^-1(Lambda_Y)=Lambda_X.                                  (2.1)

The finite etale torus cover of Y corresponding to Lambda_Y pulls back to a disjoint union of copies of the torus cover of X corresponding to Lambda_X. A connected component is the canonical torus cover of X. This identification is algebraic by the finite-cover correspondence. Its map to the torus cover of Y is dominant and regular. Such a torus map is a translation followed by a full-rank integral monomial map, hence a finite etale isogeny. The full pullback is finite etale, so f is finite etale by descent. The quotient of (2.1) embeds F_X into F_Y. In particular

    deg(f)=[F_Y : image(F_X)] [Lambda_Y : f_*(Lambda_X)].         (2.2)

For clarity, the canonical torus cover is classical: Jia--Shibata--Xie--Zhang, arXiv2104.05339v1, Lemma2.20. It is the translation-free Galois torus closure. We do not claim that cover as a new construction. Automorphisms lift to it because Lambda is characteristic. If X=T/F has faithful linear holonomy and T is this cover, then

    Aut(X)=N_Aut(T)(F)/F.                                      (2.3)

The finite-cover correspondence gives the algebraic lifts; their conjugation normalizes the deck group, and the ambiguity is exactly F. Conversely every normalizing automorphism descends. These assertions descend from C to every characteristic-zero field after geometric base change. Frobenius rules out the same assertion in arbitrary positive characteristic.

## 3. UAQUARTAUT1: the fixed universal quartic has no canonical S4 repair

For the specific W above, the complete geometric automorphism group is

    Aut(W) = C2 x (Gm semidirect C2) x (mu3 semidirect C2).       (3.1)

The involutions in the last two factors act by inversion. This group is metabelian (its second commutator subgroup is trivial), so it has no subgroup isomorphic to S_4, whose second commutator is V_4.

Proof. In logarithmic coordinates x,z,v1,v2 (v3=(v1*v2)^-1), the linear parts of r and s are

    R=diag(1,1,[[-1,-1],[1,0]]),
    S=diag(-1,1,[[0,1],[1,0]]).

The rational representation is the direct sum of the sign, trivial and two-dimensional standard representations of S_3. They are pairwise nonisomorphic and absolutely irreducible. Hence their rational linear centralizer consists of diag(a,b,c,c), and its integral unimodular subgroup has a,b,c in {+1,-1}. An automorphism of W lifts by (2.3). Every automorphism of S_3 is inner and S_3 has trivial centre; multiplying the lift by a unique deck transformation therefore gives a lift commuting with the actual affine S_3 action. The centralizer maps isomorphically to Aut(W).

Let the three linear signs be epsilon_x,delta,eta. Commutation with the translation x->theta*x forces epsilon_x=+1. Commutation with s forces the remaining x translation to be epsilon in mu2. In the A factor the translations fixed by both permutations are (zeta,zeta,zeta), with zeta in mu3. The z translation c is arbitrary. Every resulting map

    (x,z;v_i) -> (epsilon*x, c*z^delta; zeta*v_i^eta),
    epsilon in mu2, c in Gm, zeta in mu3, delta,eta in {+1,-1}, (3.2)

commutes with both r and s, and these are all the lifts. This proves (3.1), including its group law.

It follows that **no map having source W can carry a faithful natural S_N marked action for N>=4**, whatever its degree or target: faithfulness on the labels makes the homomorphism S_N->Aut(W) injective, while S_N is not metabelian.

More strongly, no finite etale endomorphism of W x A^r admits a canonical natural S_N datum for N>=4, for any r and arbitrary finite coefficient coupling. By Section1 it is (f0 x id) after a source automorphism. Its complete fibre projects bijectively onto the fibre of f0. Every source automorphism projects to an automorphism of W; the kernel of this projection fixes every point of the finite fibre, since each projection fibre contains only one marked point. Faithfulness of the natural label action would therefore give an injection S_N->Aut(W), impossible by (3.1).

This closes the canonical quartic gate on the constructed quotient, on every identity cylinder and on every finite coupled cylinder, including left-right equivalents. It does not exclude different torus quotients with larger automorphism groups. The geometric S_4 monodromy of Phi comes from the Galois closure T->T/S_4; it is not a faithful S_4 automorphism action on W. Thus universality and strong absolute stable atomicity remain intact.

## 4. UASYLOWDEF1: a classification-free degree budget for every transitive group

Let a finite group G act transitively on a set Omega of N>1 labels. If ell divides N, every Sylow ell-subgroup P of G has **no common fixed label**. Indeed, a label stabilizer has order |G|/N, whose ell-part is strictly smaller than |P|. No conjugate of that stabilizer can contain P. Individual elements of P may nevertheless have fixed labels; the assertion is about the whole subgroup.

Apply this observation to a canonical transitive polynomial etale datum on affine space, with complete marked rank N and generic degree M. UASMITHFIX1 gives smooth connected mod-ell acyclic whole P-fixed schemes and an etale restriction between them. The fixed target is present, but its complete fixed fibre has rank zero. If d_P is the generic degree of this restriction, the strict-henselian normalization budget of UAFIXDEFECT1 gives

    1 <= d_P <= M-N.                                          (4.1)

UAACYCLICRIG1 excludes d_P=1. In characteristic zero, if N is even choose ell=2; UAACYCLICGAL1 excludes d_P=2 as well. Hence, for **every finite transitive label group and arbitrary polynomial actions**,

    M >= N+2                       (characteristic zero),
    M >= N+3                       (characteristic zero, N even). (4.2)

Equivalently M>=N+2+1_(2 divides N). This extends UAPARITYDEF1 from even natural symmetric groups to all even-degree transitive groups. For odd N every 2-group action on the labels has a fixed point, so this particular empty-fibre argument supplies only the +2 bound. Attainment of either bound is not asserted.

The general transitive +2 argument no longer requires Fein--Kantor--Schacher's prime-power derangement theorem or the classification of finite simple groups: whole Sylow subgroups suffice. The earlier cyclic-derangement proof remains valid, but its classification-theoretic input is unnecessary for this conclusion.

In characteristic p>0, whenever N is not a power of p choose a prime ell|N different from p. Tame Smith theory, (4.1) and birational acyclic rigidity still give M>=N+2, for **every** transitive label group. Thus full-degree or one-surplus canonical transitive data can occur only in p-power rank. No positive-characteristic +3 assertion is made here; the quoted quadratic obstruction was proved in characteristic zero.

A useful exact control is the transitive six-point action of A_5 on cosets of a dihedral subgroup of order10. Its Sylow2 subgroup V_4 has no common fixed label, although each of its three nonidentity elements fixes two labels. The Sylow proof applies directly even though this subgroup contains no cyclic 2-power derangement. This is not a bounded-search proof of (4.2); it checks the distinction on which the general written proof relies.

## 5. UAVERTRED1: a rigid metabelian base returns the canonical problem to affine space

Let Z be smooth geometrically integral affine and geometrically A1-rigid in characteristic zero, with metabelian geometric Aut(Z). Let N>=4. Over an algebraically closed field, a canonical natural S_N datum on an arbitrary dominant etale endomorphism

    G:Z x A^r -> Z x A^r

restricts to such a datum on a single pair of affine fibres A^r->A^r. In particular, if Z has a rational point, existence of these split canonical data on Z x A^r is **equivalent**, at the same r, to existence on A^r itself. Over a nonclosed characteristic-zero field the same assertion holds for split natural S_N data and a rational point of Z. The reverse implication simply takes the product with id_Z and lets S_N act trivially on Z. This statement concerns canonical equivariant marked data; it does not by itself construct a family over the moduli stack B S_N.

Proof. By Section1 all automorphisms project to Z and G has form (1.1). Every homomorphism S_N->Aut(Z) kills the second commutator subgroup S_N''. This subgroup is V_4 for N=4 and A_N for N>=5, and is transitive on the natural labels. Therefore all N marked source points have the same projection z0 in Z. Since the full S_N permutes them, z0 is fixed by its projected source action. The target projection y0 is fixed by the projected target action. Source and target affine fibres over z0,y0 are invariant under the full group. The block differential of (1.1) makes the restriction

    g_z0:A^r -> A^r                                            (5.1)

etale. Its fibre over the marked target is exactly the complete original N-point fibre, with the same natural action. An etale polynomial self-map of affine space has nonzero constant Jacobian, which can be normalized to1 by a linear target change, transporting the target action. Thus (5.1) is the required ordinary polynomial Keller datum.

For descent, split natural data have all labelled points rational, so their common projection z0 is rational. Even when phrased as a geometric natural full S_N fibre with constant S_N action over the field, Galois permutes the labels through the centralizer of S_N in Sym(N), which is trivial for N>=3; this again makes the points rational. This argument is only about this split object, not a claimed canonical choice of a point on every input torsor.

The general mechanism is slightly broader: any normal transitive subgroup of a prescribed label group that is killed by the projected source action forces the same vertical restriction. Metabelianity is a sufficient condition for natural S_N, not a necessary condition for all groups.

### Dimension and degree consequences

Let c(N) be the retained arbitrary-action affine lower bound from UANONLINEARCOST1:

    c(N)=N+1  if N is a prime power;
         N    if N is a sum of two positive powers of one prime,
              and is not already a prime power;
         N-1  otherwise.

Every datum above requires **r>=c(N)**; it is the number of affine coordinates, not dim(Z)+r, which must pay the cost. In particular for the four-dimensional W, a quartic datum would require r>=5 and total dimension at least9. This is only a necessary bound, with no construction asserted at9.

There is also a multiplicative degree budget. In (1.1) put D=deg(f0), h=deg(H) and h0=deg(g_z0), so M=D*h. Then h0<=h. To see this, pass to a geometric generic point of the target affine fibre over z0. Each of the h0 points of (5.1) is an etale point of H over that point. An etale generically finite map of degree h between normal integral varieties has at most h geometric points in any fibre: Zariski Main embeds it in its finite normalization, and strict-henselian local factors at its etale points each consume one generic sheet. This is the same point-count argument as UAFIXDEFECT1; scheme-theoretic special-fibre length need not equal the generic rank. Applying Section4 to (5.1) gives

    M >= D * (N+2+1_(2 divides N)).                            (5.2)

For quartic data on W cylinders this gives M>=7D>=7. If G is finite, Section1 instead gives h0=1, an immediate contradiction. Thus even an arbitrary nonfinite coupled map can help only by containing the original higher-degree affine canonical datum in an invariant affine fibre. Merely enlarging the torus-quotient base or changing its coefficient functions does not remove that problem.

### A precise connection to JC(r)

If Z is a geometric torus quotient, Section2 makes every dominant base endomorphism f0 finite etale. If JC(r) holds in characteristic zero, every vertical map g_z in (1.1) is an automorphism. The etale map H is then geometrically bijective and hence an isomorphism (a surjective radicial etale map). Consequently every etale endomorphism of Z x A^r is finite. Conversely an etale nonautomorphism of A^r gives the nonfinite map id_Z x g, since a finite etale self-map of A^r would be an automorphism. Thus, over an algebraically closed characteristic-zero field,

    all etale endomorphisms of Z x A^r are finite <=> JC(r).    (5.3)

This is a conditional equivalence and reduction, not a proof or disproof of JC(r). Canonical natural symmetric data impose additional requirements. It does imply unconditional finiteness for r<=1 and rules out the canonical data in those cylinder dimensions.

## 6. What is closed and what remains

The full canonical S_4 action on the particular universally quartic W self-map is now excluded, including every finite coupled affine cylinder. Its universality and strong absolute stable atomicity survive. All dominant regular maps of equal-dimensional torus quotients are finite etale. Whole Sylow fixed schemes simplify the general degree argument and strengthen it for every even-degree transitive group. For any nonfinite attempt on a rigid metabelian base, the canonical problem reduces to an actual polynomial Keller datum on one affine fibre, with the same r and the budgets above.

Still open: canonical quartic self-receivers on different geometries; the original higher-degree nonfinite whole-affine-space boundary; effective simultaneous fixed-locus landing and gluing; atomicity for maps whose marked rank is smaller than their generic degree; the modular p-power cases; and the A1-rich pseudo-plane cubic mixing gate. Neither the exact automorphism calculation nor the Sylow simplification is a certificate that no similar result has appeared in the literature.
