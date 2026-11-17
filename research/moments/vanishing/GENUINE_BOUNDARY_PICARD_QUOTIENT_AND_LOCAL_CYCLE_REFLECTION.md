# Genuine boundary Picard quotient and local cycle reflection

16 November 2026. Written geometric deduction, using the published Hodge/Picard
comparison specified below. Both FULL exceptional Mathieu--Zhao statements
remain **OPEN**. Independent replay, formal verification and external review
are not claimed; external review remains deferred.

The genuine, UNCUT finite-fibre boundary supplies an actual Picard quotient
for EVERY original phase. Its kernel is proper, so it loses no rational
generic torus directions and reflects the local corrected component. This
removes an existence hypothesis from the preceding raw-degree calculation.
It does not compute the original Haar section from an arbitrary clutching
symbol, or prove that purity makes its degrees satisfy the required equations.

## 1. The quotient comes from the entire finite-fibre boundary

Let X=SL2/Gamma, with Gamma=2O or 2I, and let f be ANY nonconstant original
regular phase. Work on its full regular parameter locus, and shrink and take
a finite cover when needed to spread the following geometric data. For a
smooth fibre F choose a smooth PROJECTIVE surface compactification Y with
reduced SNC boundary C=Y minus F. Keep EVERY component of C. Disconnected
fibres are handled by disjoint unions throughout, with all branch groups
retained. No product norm replaces the branch tuple.

Write B=H^3(X minus F;Q), L=H^3(X;Q)=Q(-2). The retained canonical reduction
gives E=B/C_f, where W4 B=L direct_sum C_f, and

    0 -> L -> E -> (H^2(F)/W2 H^2(F))(-1) -> 0.       (1)

Set H=(E(2))^dual. The actual boundary-motive theorem supplies its algebraic
1-motive [Z -> G_f], up to isogeny, with its ORIGINAL ambient-unit section
sigma_f. Our new statement is

    G_f ~= G_C := Pic^0(C) / I, up to isogeny,
    I=im(Pic^0(Y) -> Pic^0(C)).                       (2)

The comparison in (2) is the full negative-weight mixed Hodge structure,
not merely an equality of genus and torus ranks. The section sigma_f is
additional extension data; (2) alone does NOT specify it.

To prove (2), use the MHS localization sequence

    H^2(Y) -> H^2(F) -> H^3_C(Y) -> H^3(Y).

Proper support duality on the smooth surface identifies

    H^3_C(Y)=H_1(C)(-2), H^3(Y)=H_1(Y)(-2).

The support map is induced by C -> Y. Deligne's compactification-image
description identifies im(H^2(Y)->H^2(F)) with W2 H^2(F). Hence, as FULL MHS,

    H^2(F)/W2 = K(-2),
    K=ker(H_1(C) -> H_1(Y)).                          (3)

Dualizing (1), twisting as prescribed, and using exactness gives

    W_-1 H = K^dual(1)
            =coker(H^1(Y) -> H^1(C))(1).             (4)

For the proper nodal curve C, the functorial Picard comparison is

    H_1(Pic^0(C);Q)=H^1(C;Q)(1),
    H_1(Pic^0(Y);Q)=H^1(Y;Q)(1).                     (5)

Here and below H_1 of a semiabelian variety means its Betti realization,
with its full MHS. For C one can see the entire comparison directly:
normalization with paired node branches computes both the simplicial
Hodge complex and line-bundle descent. The exponential sequence gives

    Pic^0(C)^an=H^1(C,O_C)/H^1(C,Z(1)).

The Hodge quotient H^1(C,C)/F^1 is H^1(C,O_C), as follows from the same
normalization complex. The lattice, Hodge filtration and node-torus weight
filtration are therefore precisely those of H^1(C)(1). The node identifications
retain the extension between the graph and normalization terms. Restriction
from Y uses the same exponential map, so this comparison is functorial.

For the algebraic structures and full comparison we use the already retained
published input [Barbieri-Viale--Rosenschon--Saito, Annals 158 (2003)](https://annals.math.princeton.edu/2003/158-2/p06),
Construction 3.1, Lemma 3.3, Corollary 5.2, Theorem 5.3 and comparison 5.4.
The normalization complex is a proper smooth simplicial presentation of C;
its identity Picard group is Pic^0(C). These statements include the extension,
not only the associated gradeds. Proposition 1.5 supplies point full
faithfulness up to isogeny. Selected complete statements/proofs have been
reread in the cached published PDF; the foundational Hodge theory is imported.

The image I is a connected proper abelian variety. Quotienting Pic^0(C) by
I exists as a semiabelian variety. Its H_1 is exactly (4), by (5) and exactness
of the rational realization. Point full faithfulness identifies it with G_f
up to isogeny. There is also a direct algebraic reading: use the relative
normalization complex for (Y,C) in BRS Construction 3.1 with r=2. Its
semiabelian cokernel is precisely Pic^0(C)/im Pic^0(Y), realizing the negative
part of H^2(Y,C)(1)=H_c^2(F)(1). Surface duality identifies that negative
part with (4). The Picard quotient thus comes from an algebraic relative
construction, with its full Hodge comparison. This is an identification between TWO already algebraically
constructed groups. In families use their natural horizontal comparison,
point full faithfulness and the Hom-rigidity/spreading argument of the
[boundary-motive note](ALGEBRAIC_BOUNDARY_ONE_MOTIVE_AND_CONSTANT_SECTION.md),
after shrinking and clearing denominators. No unrestricted relative essential
surjectivity theorem is assumed. Thus (2) is an algebraic generic-field
identification, not only a fibrewise equality of Betti dimensions.

## 2. The proper kernel preserves every rational graph direction

Let Gamma_C be the dual graph of the ORIGINAL generic finite-fibre boundary C.
Normalize C to obtain the actual semiabelian sequence

    0 -> T_Gamma -> J=Pic^0(C) -> product_v Jac(C_v^normal) -> 0,
    X^*(T_Gamma)=H_1(Gamma_C;Z).                      (6)

The graph includes all loops and parallel edges. Since I is proper,
I intersect T_Gamma is finite. The quotient map J -> G_C therefore induces
an ISOGENY on the generic tori. It is onto on the torus: its cokernel would
be both a quotient of the abelian J/T_Gamma and affine. Consequently

    X^*(T_(G_f)) tensor Q ~= H_1(Gamma_C;Q).           (7)

Choose an actual isogeny from G_C to G_f realizing the rational comparison,
and compose it with J -> G_C, obtaining an actual surjection

    psi:J -> G_f.                                    (8)

Its kernel has proper connected component I and finitely many components;
it is proper. It may contain finite torus torsion. No rational graph direction
is killed. This statement applies to C=Y minus F. An AUXILIARY boundary made
by deleting extra phase-dependent divisors can have a torus kernel, as the
already safe top-row Klein-form control demonstrates. Such cuts do not
satisfy this conclusion by merely having a visible graph.

Changing the smooth projective compactification changes J and I but not
their quotient up to its canonical realization (4). This follows directly
from (1)--(4), without identifying an ambient resolution tree with Gamma_C.

## 3. Every original section has an actual lift after finite cover

The fibre of psi over any positive integer multiple [a]sigma_f is a nonempty
finite-type scheme over the generic parameter field. A closed point gives
a lift after a FINITE field extension. Thus, on a finite cover and after
shrinking, there is an actual section ell of J such that

    psi(ell)=[a]sigma_f, a>0.                         (9)

The multiplier records clearing of rational lattices/isogenies. It is NOT a
generic degree of a map. The map of 1-motives with lattice map a and group
map psi preserves the original extension because (9) is imposed using the
ACTUAL sigma_f. A hypothetical line bundle with the same node orders, or
the original contact primitive on a cut open set, does not establish (9).

A J-valued point is a relative Picard class. After a further finite extension
it has a representing line bundle on C, by the fppf Picard definition and
a closed point in a finite-type local trivialization (equivalently one may
choose a rational point to rigidify the proper curve after finite extension).
Thus an actual generic line bundle L_K representing ell exists on the cover.
Its degree on EVERY generic normalization is zero. No generic normalization
triviality is asserted.

Near the original phase value 0 write s=u^e after completion and the finite
covers. The [actual torus-component theorem](ACTUAL_TORUS_CORRECTIONS_AND_RAMIFIED_COMPONENT_VALUATIONS.md)
identifies rational Neron component groups with the inertia-invariant generic
cocharacters. Equation (7) therefore makes (8) an isomorphism on these rational
component groups. With original-parameter normalization,

    psi_* nu_J(ell)=a nu_(G_f)(sigma_f),
    nu_J(ell)=0 iff nu_(G_f)(sigma_f)=0.              (10)

Different lifts in (9) have the same rational component: their difference
lies in the proper kernel, whose identity component is abelian and whose
component group is finite. Its image in a rational Neron component group
is zero. This can also be read off from the isomorphism in (10).

This is LOCAL rational-component reflection. It does not say that an
arbitrary lift ell has globally constant type whenever sigma_f does.
Nor does it reflect a local splitting through the discarded compact-middle
variation C_f. That quotient and the proper abelian kernel in (8) are
different operations with different local properties.

## 4. Uniform raw-degree input for the intrinsic Mathieu calculation

After further finite base change, take regular strict semistable models Y_v
of all generic normalizations, with disjoint node half-sections in their smooth
loci. Glue the persistent paired sections to a proper flat nodal curve over
R=C[[u]]. Keep the entire special normalization fibre sum_j E_(v,j).
The preceding [raw Picard lemma](RAW_PICARD_DEGREES_AND_A_CONSTRUCTED_IDENTITY_CORRECTION.md)
extends L_K invertibly after marked-point blowups, retaining all exceptional
degrees. Put

    delta_(v,j)=deg(L|E_(v,j)),
    M_v=(E_(v,j) dot E_(v,k)), M_v d_v=delta_v.       (11)

The solutions exist over Q and are unique modulo a constant on each connected
normalization special fibre. Intersection balancing and unit node gluing
construct the ACTUAL identity correction. Therefore, for every original torus
character chi, the intrinsic coupling is

    nu_f(chi)=1/(a e) sum_edges (psi^*chi)_b
                   (d_(tail,attachment)-d_(head,attachment)).  (12)

Unlike an arbitrary supplied Picard presentation, (7) says that the pullbacks
psi^*chi span ALL H_1(Gamma_C;Q). Consequently, under eventual purity,

    ALL original regular markers safe
      iff sum_edges z_b(d_tail-d_head)=0
          for EVERY z in H_1(Gamma_C;Q).              (13)

Equivalently sum_v delta_v dot q_v(z)=0 for every z, where M_v q_v(z) is the
attachment current of z. Equivalently SOME positive tensor power of this
genuine lift L_K has an invertible extension of degree zero on EVERY special
normalization component. The tensor power retains finite component torsion.
The equivalences use the retained full intrinsic local Tate theorem, not just
a sufficient generic-fibre semisimplicity test. Individual nonzero degrees
may cancel in the exact cycle pairing; requiring all given degrees zero is
only a sufficient test.

Equations (9),(11),(12) are a geometric interface for EVERY phase. They prove
existence of the correct presentation and lift, not an implemented producer
of the original Haar NS/Picard data and not a uniform purity-to-zero theorem.
The original volume is retained in sigma_f and in a; it is not replaced by
an arbitrary generator of graph cohomology.

## 5. Exactly what a good-component contraction must justify

There are three separate graphs: the ambient compactification/resolution
graph; the generic finite-fibre boundary graph Gamma_C in (6); and each
normalization's special-fibre intersection graph determining M_v in (11).
An ambient tree imposes no tree assertion on Gamma_C. An exceptional P1 in
a normalization model can carry essential degree and attachment information
even when that normalization's own intersection graph is a tree.

For a proposed good subgraph Gamma_good of Gamma_C, calculate (12) on
H_1(Gamma_good). The [bounded graph descent lemma](BOUNDED_GRAPH_DESCENT_AND_GOOD_NODE_EXTENSION.md)
then gives the exact conditional localization:

    all those corrected good-cycle periods vanish
      iff, modulo a bounded torus point, the corrected toric section
          comes from the torus of Gamma_C/Gamma_good.  (14)

All surviving attachment edges, loops and parallel edges remain after the
contraction. The contracted torus embeds; this is not a surface blowdown
or an automatic quotient of the full variation. A forest after this justified
contraction suffices for zero coupling. If cycles remain, their periods still
must be evaluated. Finite-base holomorphicity of the ambient volume by itself
has not been proved to imply the first condition in (14).

The actual nonexceptional Long phase supplies a sharp control: it is pure,
its good infinity components contract, yet the remaining attachments give
nu=1. On the degenerating hyperbola normalization the two-component
intersection graph is a tree, the raw degrees are (1,-1), and opposite
attachments lie on different components. Equation (12) retains the nonzero
period. Dropping the exceptional component or summing its degrees loses it.
This is NOT an exceptional invariant counterexample.

## 6. The outstanding exceptional implication

For the ENTIRE exceptional pure locus, the remaining obligation is to identify
the actual lift (9), or sufficient original NS/Picard data for it, and prove
the pairings (13) vanish. The dominant-left-root and fixed-polar-divisor
results remain inputs to locating these data. Neither has yet supplied
the necessary zero good-cycle periods in (14) for every pure phase.

Alternatively one must construct an actual original right-2O or right-2I
invariant pure phase with a NONZERO pairing and one fixed original regular
bad marker. An abstract Picard section, a nonzero cut-symbol residue, a
nonpure phase with nonzero coupling, or a finite list of degrees does not
meet that requirement. The family-transfer theorem, capacity, automorphism,
primitive-branch and fixed-marker constraints are preserved. Ordinary marked
resolvent algebraicity is not inferred from Gaussian exponential generating
functions. Both full exceptional cases and the all-pure vanishing remain OPEN.
