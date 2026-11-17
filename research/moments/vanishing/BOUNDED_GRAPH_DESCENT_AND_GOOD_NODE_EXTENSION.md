# Bounded graph descent and good-node extension

16 November 2026. Written local geometric lemma for the ORIGINAL
Mathieu unit-extension calculation, in every degree. Both FULL
exceptional cases remain **OPEN**. External review remains deferred.

The usable bridge is precise: contraction of a good fibre subgraph
preserves the corrected component class EXACTLY when its good-cycle
holonomies are units. An invertible extension of the ACTUAL corrected
toric line bundle across a proper equisingular good boundary curve
proves that condition, even when its graph has cycles. A nonfree
node extension is an explicit obstruction. Original volume regularity
alone has not been shown to rule it out on every pure phase.

## 1. Input is the actual corrected original section

Fix an original phase f on X=SL2/2O or SL2/2I, with its original
volume and canonical geometric boundary 1-motive. At phase value 0
take a finite splitting/labeling cover s=u^e. Put K=C((u)), R=C[[u]].
Use the ACTUAL generic semiabelian decomposition

    0 -> T -> G -> A -> 0,
    t=[n]sigma-b in T(K),
    b in G^0(R), pi(b)=[n]pi(sigma).                  (1)

This is the identity-component correction of
[GHAARTORICCORR20261116](ACTUAL_TORUS_CORRECTIONS_AND_RAMIFIED_COMPONENT_VALUATIONS.md),
with n killing the abelian specialization component. Its invariant is

    v=val_T(t)/(e*n).                                (2)

No characters on all of G, good reduction of A, or splitting G=T times A
are assumed. T is the GENERIC toric subgroup, not a larger special-fibre
torus caused by abelian degeneration. The retained theorem identifies
v=0 with zero intrinsic reduced local unit obstruction; under eventual
purity this is equivalent to safety of ALL original markers.

The graph below must realize this ACTUAL torus, up to the fixed rational
comparison/isogeny, through the boundary of the ORIGINAL nearby fibre.
It is not the cut graph of an auxiliary {FG!=0} fibre. Clear isogeny
denominators once when necessary. All zero tests are rational and are
preserved and reflected by the finite cover and isogeny.

After labeling, let Gamma be this finite graph and identify

    Lambda=H_1(Gamma,Z),
    T_Gamma=Hom(Lambda,Gm).                           (3)

This is the usual toric subgroup of the graph Picard variety: component
normalizations have trivial line bundles, and node scalars are taken
modulo vertex changes of trivialization. In an actual further character
subquotient use its supplied lattice and the saturated good-character
sublattice; the lattice argument below applies unchanged. No comparison
with the intrinsic domain is inferred from a drawing alone.

Represent t by actual oriented node gluing scalars g_a in K*. Then

    t(z)=product_a g_a^(z_a),
    v(z)=(1/(e*n)) sum_a z_a ord_u(g_a)              (4)

for each cycle z. Multiplying vertex trivializations changes the edge
orders by a coboundary and leaves (4) unchanged. All components and
attachments are retained at this stage.

## 2. Exact bounded-descent lemma for graph contraction

Choose a good subgraph Gamma_good. Include every other vertex as an
isolated vertex and contract each connected good component to one
vertex. Retain all other edges, including loops and parallel edges.
Write p:Gamma->Gamma_c for this contraction. There is an INTEGRAL
exact sequence

    0 -> H_1(Gamma_good) -> H_1(Gamma)
                         -> H_1(Gamma_c) -> 0.        (5)

For a proof, lift a cycle in Gamma_c by joining its entry/exit vertices
along paths inside each connected good component. The result is a
cycle upstairs. Two choices differ by good cycles. A cycle killed
by contraction is a sum of cycles inside the good components. This
proves both surjectivity and the kernel assertion. The quotient is
free, so the good-cycle lattice is primitive.

Applying Hom(-,Gm) gives

    1 -> T_(Gamma_c) -> T_Gamma -> T_good -> 1.       (6)

In particular the contracted graph torus EMBEDS in the original torus.
Graph collapse alone is not a map sending an arbitrary actual section
to a contracted section. Its good holonomy has to be removed first.

The following conditions are equivalent:

    (a) ord_u t(z)=0 for every good cycle z;
    (b) the image of t in T_good(K) belongs to T_good(R);
    (c) there is a in T_Gamma(R) with
                         t/a in T_(Gamma_c)(K).      (7)

Condition (a) is exactly (b) in a character basis. For (b)->(c), the
primitive lattice inclusion in (5) makes T_Gamma(R)->T_good(R)
surjective: extend a basis of the good lattice to one of Lambda and
choose arbitrary unit values on the additional basis elements. Lift
the good image to a and divide. The converse is immediate from (6).

This also proves support localization at the cochain level. Write
w_a=ord_u(g_a). Vanishing on good cycles is equivalent to the existence
of integral vertex potentials p_v with w_a=p_head-p_tail on good
edges. Define p along a spanning tree in each good component; the
remaining edges satisfy the equation exactly because their cycle
sums vanish. A vertex gauge using powers of u then makes every good
edge order zero. All remaining nonzero orders are supported on the
uncontracted edges. This is a proved gauge operation on the ACTUAL
corrected gluing, not deletion of an untested boundary contribution.

Different bounded lifts a give the SAME class

    [t_c] in T_(Gamma_c)(K)/T_(Gamma_c)(R).            (8)

Indeed their quotient lies in the embedded contracted torus and is
bounded there: all its characters have valuation zero. Consequently
the normalized contracted vector is canonical and

    v=p^* v_c,
    v=0 iff v_c=0,                                   (9)

provided the equivalent conditions (7) hold. Here p^* is the injection
H^1(Gamma_c,Q)->H^1(Gamma,Q), and normalization is still 1/(e*n).
Changing b,n or a further parameter cover changes representatives
by bounded points and the established normalization only. Thus (8)
and (9) apply to the ORIGINAL intrinsic component class in (2).

This is LOCAL bounded descent. It does not produce a global quotient
of the algebraic variation on the punctured parameter line. The
global unit-retaining quotient theorem must not be applied to an
arbitrary local graph collapse. Local equivalence in (9) follows
instead from the actual corrected component criterion (2). A genuine
global section factorization through the embedded contracted torus
would supply the corresponding global 1-motive map, but is not assumed.

The explicit obstruction to this descent is

    o_good = v restricted to H_1(Gamma_good,Q).        (10)

It is an actual corrected cycle functional, not a residue vector of
the unit form. If o_good!=0, contracting that good subgraph erases
part of the obstruction and is invalid. If it is zero, (8) gives
the legitimate residual component calculation on Gamma_c.

No assertion of an algebraic surface blowdown is made. The contraction
here is the graph/Picard reduction needed for the component invariant.
An actual morphism contracting surfaces or positive-genus curves needs
its own geometric hypotheses.

## 3. A geometric criterion that PROVES good holonomies are units

Suppose the proposed good boundary curve has a proper flat model
mathcal C/R whose normalization components mathcal C_v are smooth
proper curves with geometrically connected fibres. Suppose its nodes
are disjoint pairs of sections of these normalizations: the generic
and special dual graphs are the same Gamma_good. Attachments to the
remaining boundary are kept as marked sections. These are explicit
equisingular-model hypotheses, not consequences of a regular volume.

Restrict the ACTUAL corrected toric point t to its good Picard torus.
Its corresponding generic line bundle L_K is trivial on every
normalization component. Then

    o_good=0
      iff L_K extends to an INVERTIBLE line bundle
          on the entire mathcal C.                   (11)

This is a necessary and sufficient extension criterion for this
equisingular model, rather than an assumption that Gamma_good is a tree.

To prove the reverse implication in (11), pull an invertible extension
back to each regular surface mathcal C_v. Its generic restriction is
trivial. A generic trivialization has a divisor supported on the
special fibre. That fibre is irreducible and is div(u), so this
line bundle is trivial on all of mathcal C_v. Choose such R-trivializations.
The node identifications of an invertible sheaf are now isomorphisms
of free rank-one R-modules. Every node scalar is a UNIT. Hence every
cycle product is a unit and o_good=0.

Conversely o_good=0 means the torus point is an R-point by (7)(b).
Represent it by unit node scalars after vertex gauge changes. Glue
trivial line bundles on the normalization components using these
units. At a node the module of pairs satisfying f(0)=g_a h(0) is
free rank one when g_a is a unit. Thus the gluing is an invertible
line bundle on mathcal C, with the prescribed generic fibre. This
proves the other implication.

Units need not equal 1 and need not be torsion: a constant gluing
parameter 2 is allowed. The statement tests the rational local
component, not splitting of a mixed Hodge structure.

Cycles, and genera of normalization components, are allowed in (11).
Its necessity uses the full supplied good Picard-torus realization.
If only a proper actual character domain is being tested, such an
invertible extension remains sufficient; no stronger necessity for
irrelevant extra graph directions is asserted.
The generic normalization bundles are trivial because the ACTUAL
abelian correction in (1) has already been made. One may not impose
this on sigma before subtracting b. If normalizations degenerate,
nodes collide, or extra vertical components appear, all those data
must be retained; (11) is not asserted for an unverified model.

## 4. Positive consequence even when the full fibre graph has cycles

If the extension condition (11) holds and Gamma_c is a forest, then
H_1(Gamma_c)=0 and (9) gives v=0. Under eventual purity every original
marker is safe by the intrinsic local theorem. This is an ALL-degree
geometric safety mechanism permitting arbitrarily many cycles INSIDE
the good subgraph. It does not require the full boundary graph to
be a tree.

More generally, after (11), the ONLY remaining local class is v_c.
Its pairing with every cycle of Gamma_c must still be computed. A
tree-shaped bad piece is not enough: attachments through one good
component can turn it into a loop after contraction.

The geometric gate supplied to the existing Mathieu calculation is
therefore: identify its ACTUAL corrected good line bundle, verify
the model hypotheses and invertibility at the good nodes, then use
(8)--(9) for the residual finite-fibre graph. This keeps the original
Haar unit and removes only a proved bounded component. No Jacobian
weight, volume replacement, marker loss or branch norm is involved.

## 5. Explicit obstruction: a nonfree node extension

Let mathcal C be a constant proper rational nodal curve over R with
one self-node; its normalization is P1_R with 0 and infinity identified.
The graph is one loop. The generic toric line bundle with gluing
parameter u^m, m>0, has good-cycle valuation m. It cannot extend
invertibly, by (11). This failure occurs even though the entire
curve and its normalizations have a constant smooth model away from
the persistent node.

The local defect is explicit. With

    A=R[[x,y]]/(xy),

the extension inside the two branch lattices is

    M_m={(f(x),g(y)):f(0)=u^m g(0)}.

It is R-flat (torsion free) and generically free rank one. It has
generators e_1=(x,0), e_2=(u^m,1) with

    x e_2=u^m e_1, y e_1=0.                          (12)

At the closed node both generators survive in M_m/(u,x,y)M_m. Its
fibre dimension is 2, so it is NOT invertible. The gluing-lattice
cokernel R/(u^m) has length m. Negative exponents exchange the two
branches. On a graph an edge exponent changes by a vertex coboundary;
the invariant defect is its cycle sum, not a chosen individual edge.

This is the Kummer 1-motive [Z->Gm],1->u^m, with nonzero rational
component and local logarithm m log u. It is a geometric obstruction
model, NOT an actual exceptional pure phase or an MZ counterexample.
It shows exactly what a proposed support-localization proof must
exclude in its ACTUAL node gluing. Extending merely as a torsion-free
sheaf is insufficient; invertibility is the crucial missing property.

## 6. Actual original-volume control: a valid contraction leaves a nonzero cycle

The retained nonexceptional Long phase makes the residual obstruction
geometric, not merely a free choice of torus section. On ORIGINAL SL2 put

    p=(1+y)(1+z+yz), v=1+y, s=p.

It is pure in every positive order and has the fixed original bad
marker Q=-y, with h(Q p^j)=(-1)^(j-1)/(j+1). Its actual boundary
1-motive has genus zero, torus rank one and Kummer section

    sigma(s)=s/(s+1),                                (13)

up to a nonzero constant and the retained rational isogeny normalization.
These are completed inputs from
[the boundary 1-motive proof](ALGEBRAIC_BOUNDARY_ONE_MOTIVE_AND_CONSTANT_SECTION.md),
Section 4. No new purity claim about an exceptional phase is made.

For s!=0,-1 the fibre has v!=0. Set X=xv,W=wv. The determinant equation
gives

    XW=(s+1)v-s,
    F_s ~= A2_(X,W) minus {XW=-s}.                   (14)

Starting from the ORIGINAL invariant volume dx wedge dy wedge dz/x,
one direct change of variables gives, up to its one fixed Haar scalar,

    Omega=(s+1) dX wedge dW wedge ds/(XW+s)^2.       (15)

This does not replace volume by a phase-dependent multiple. In the
generic compactification P1_X times P1_W, the boundary is the triangle
formed by H_s={XW=-s}, D_X={X=infinity}, D_W={W=infinity}. H_s has the
double volume pole. The other two components are volume-regular away
from their attachments to H_s; at their mutual node, with U=1/X,V=1/W,

    Omega=(s+1)dU wedge dV wedge ds/(1+sUV)^2,       (16)

which is holomorphic. The good subcurve D_X union D_W is a constant
proper two-P1 tree with persistent node and marked attachments. Its
torus is trivial, so the good gate is satisfied and contraction is valid.

The contracted GENERIC FIBRE graph has two vertices, the good vertex
and H_s, joined by TWO parallel attachment edges. Its cycle survives.
By (13), the actual corrected residual component at s=0 is nonzero:

    ord_0(sigma)=1.                                 (17)

Thus purity and a valid volume-good contraction do not alone make the
remaining coupling vanish on general SL2. Exceptional hypotheses must
enter the missing vanishing proof. They cannot be replaced by the
assertion that the bad piece was itself a single rational curve.

Equations (14)--(16) describe the nearby generic compactification and
the good neighbourhood, NOT a complete resolved original model over 0.
H_s degenerates at 0, and the original special fibre has an additional
v=0 component omitted by the generic v!=0 chart. Those components and
any required vertical resolution are NOT discarded: the component in
(17) is taken from the completed ACTUAL original source (13), which
retains them. The example does not infer monodromy from the partial
generic compactification alone.

There is a retained exceptional obstruction as well: the regular
Bezout maps Phi_b from BOTH exceptional quotients have nonzero actual
H^3 scalar 2. The completed
[cohomological pullback theorem](COHOMOLOGICAL_PULLBACK_AND_A_UNIFORM_ORIGINAL_MARKER.md)
therefore transfers the nonzero intrinsic component at 0 to EVERY
g_b=p composed Phi_b, with the recorded fixed original marker C_b.
Their unweighted purity remains unproved. Moreover Phi_b^*Omega equals
a nonconstant Jacobian times the original exceptional volume. Thus
neither the volume-good pieces in (15) nor unweighted purity are
transferred by dropping that factor. Only the already proved
cohomological obstruction and fixed-marker identity are being used.

## 7. An ambient tree cannot replace this fibre calculation

For example, in A3 take the SNC surfaces x=0,y=0 and phase
f=x+y+z^d, d>=2. The ambient incidence graph has one edge. In a
generic fibre their two connected curves meet in the d roots of
z^d=s, giving d parallel edges and first Betti number d-1. This
retained geometric example concerns the actual fibre intersection
graph, not an exceptional pure phase. Resolving a proper completion
retains these finite nodes, or subdivides their edges.

Even at the graph level, take two vertices joined by a good edge
and a bad edge. The good graph is a tree, but contracting it leaves
a bad loop. Unit good holonomy permits the contraction; it does
NOT make the remaining loop's coupling zero. Conversely good cycles
with unit holonomy are harmless for the component calculation even
though their graph is not a tree.

## 8. What is proved, and the remaining exceptional obligation

The lemma proves a legitimate local bridge from actual corrected
line-bundle extension to support localization, bounded descent and
good-component graph contraction. It also gives the precise obstruction
(10),(12) when that bridge fails. Its input is the original intrinsic
section and the actual uncut fibre comparison, not the nonzero
clutching-symbol cut graph of
[GHAARLOGCLUTCH20261116](ORIGINAL_LOG_SYMBOL_AND_TWO_ORBIT_BOUNDARY_SUPPORT.md).

The UNPROVED full step is that EVERY exceptional pure phase has the
invertible extension required on all proposed volume-good pieces,
with an appropriate verified good model, and that every residual
corrected cycle pairing v_c vanishes. Holomorphic original volume
on a piece, or the dominant-root pole theorem, has not been proved
to imply this invertibility. The explicit nonfree-node model is not
claimed to arise on an exceptional pure phase, so it does not refute
that desired implication either.

All support, capacity, automorphism and primitive constraints remain
inputs. The completed family-transfer theorem, all original markers
and relevant fibre components are retained. No ordinary-resolvent
algebraicity, full pure-to-safe inclusion, or actual pure counterexample
with a fixed original bad marker is obtained. Both FULL cases remain
OPEN; review remains deferred.
