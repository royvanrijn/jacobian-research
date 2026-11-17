# Raw Picard degrees and a constructed identity correction

16 November 2026. Written all-degree geometric deduction with the
primary Picard inputs specified below. Both FULL exceptional
Mathieu--Zhao cases remain **OPEN**. External review remains deferred.

For a supplied ACTUAL Picard presentation, the abelian identity-component
correction no longer has to be supplied separately. Intersection balancing
constructs one. The corrected component can then be read directly from
the multidegrees of the UNCORRECTED line bundle, even when its generic
normalization classes are nontrivial. This is a geometric input-output
lemma, not a proof that purity forces its output to vanish.

## 1. Actual data and the scope of the presentation

Put R=C[[u]], K=C((u)), with original phase parameter s=u^e. Use a proper
nodal curve C_K with labeled generic graph Gamma. Its generic normalizations
have supplied proper regular strict semistable models Y_v/R, with geometrically
connected fibres and reduced special fibres sum_j E_(v,j). The marked half-node sections are
disjoint and meet the smooth locus of those fibres. Glue the sections to
form a proper flat nodal model C/R. Work componentwise if Gamma is
disconnected; no branch is dropped.

Write

    J=Pic^0(C_K), 0 -> T_Gamma -> J -> A -> 0,
    A=product_v Jac((Y_v)_K).                         (1)

Let ell in J(K) be represented by a generic line bundle L_K. It has degree
zero on every generic normalization, but those normalization bundles are
NOT required to be trivial. Suppose an INVERTIBLE extension L on C has
been supplied. Such an extension can be made after marked-point blowups
as in Section 6; its degrees must then retain all added components.

The direct application is J=G_f and ell=sigma_f in the actual geometric
boundary source. More generally Section 5 requires an ACTUAL algebraic
homomorphism from this J to G_f and a compatible lift of an integer multiple
of sigma_f. A same-rank torus, an auxiliary cut graph, or a Betti vector
alone does not supply these data. Original Haar normalization is retained
through the lattice multiplier in that actual map.

## 2. Why zero multidegree gives a Neron identity point

We first prove the boundedness fact used in the construction:

    If a line bundle B on C has degree zero on EVERY
    irreducible component of the special fibre, and
    B_K lies in J, then [B_K] lies in mathcal J^0(R).   (2)

Here mathcal J is the Neron LFT model of the GENERIC semiabelian J.
It is not a replacement by the special fibre's possibly larger torus.

The named primary imports are the Stacks Project:

* [Lemma 36.32.6, tag 0E0L](https://stacks.math.columbia.edu/tag/0E0L):
  universal identification of the structure-sheaf pushforward for proper
  flat families with geometrically reduced connected fibres;
* [Lemma 36.32.2, tag 0B9T](https://stacks.math.columbia.edu/tag/0B9T):
  local constancy and base-change compatibility of fibre Euler characteristics;
* [Proposition 99.11.8, tag 0D2C](https://stacks.math.columbia.edu/tag/0D2C),
  and [Section 108.9 with Lemma 108.9.5, tag 0DMD](https://stacks.math.columbia.edu/tag/0DMD):
  relative Picard representability and smoothness for the present curves.
  Its proof uses [Lemma 108.8.5, tag 0DM9](https://stacks.math.columbia.edu/tag/0DM9),
  whose obstruction argument uses vanishing of coherent H^2 on curves.

The complete relevant statements, setups and proofs were read and cached.
Only these Picard/cohomology facts are imported. In particular the separatedness
lemma requiring integral fibres is NOT used for our reducible fibres.

To apply them, the proper nodal fibres of C/R are geometrically reduced and
connected on each chosen connected piece. Thus the universal H^0 hypothesis
of Section 108.9 holds. Let P be its smooth relative Picard algebraic space.
The total degree on each Y_v is locally constant, by the Euler-characteristic
formula for degree and tag 0B9T. Let P' be the open-and-closed subgroup where
these degrees are all zero. Its generic fibre is exactly J.

The generic identity extends to a morphism P'->mathcal J: take smooth affine
etale charts of P', use the Neron mapping property on each chart, and glue
by uniqueness on their overlaps. Separatedness of the TARGET suffices;
P' itself need not be separated. The resulting morphism preserves the group
law, again by uniqueness from the generic fibre.

The special line bundle B_0 in (2) belongs to the connected identity subgroup
of Pic(C_0). Indeed degree-zero line bundles on each normalized special
component are parametrized by connected Jacobians, and their node gluing
is a connected torus. This normalization/gluing description proves
connectedness of the multidegree-zero subgroup directly. Its image in
mathcal J_0 is therefore in mathcal J_0^0. This proves (2).

This argument constructs a map from a smooth Picard space. It does not
assert that the entire Picard space IS the Neron model.

## 3. Construct the actual abelian correction from raw degrees

Let

    delta_(v,j)=degree(L|E_(v,j)),
    M_v=(E_(v,j) dot E_(v,k)).

Generic normalization degree zero and flatness imply sum_j delta_(v,j)=0.
The supplied reduced semistable model has the intersection property proved
in the [vertical-degree lemma](VERTICAL_MULTIDEGREES_AND_ATTACHMENT_PERIODS_OF_THE_HAAR_EXTENSION.md):
M_v is a negative graph Laplacian with kernel Q*1. Hence solve

    M_v d_v=delta_v over Q,
    D_v=sum_j d_(v,j) E_(v,j).                        (3)

Choose n>0 clearing all denominators of all D_v. On every normalization put

    B_v=(L|Y_v)^n tensor O(-nD_v).                    (4)

Its degree is zero on EVERY special component. Its generic class is the
n-th power of the generic normalization class of L. At each marked node
the two section fibres are free rank-one R-modules, so choose a unit
identification. Glue the B_v to an invertible B on C. Then B has zero
degree on every special component of C, and Section 2 gives

    b=[B_K] in mathcal J^0(R),
    pi(b)=[n]pi(ell).                                (5)

This is the ACTUAL identity-component lift required by the corrected
component theorem, constructed from the supplied raw Picard data. In
particular n automatically kills the relevant abelian specialization
component; no independent good-reduction assumption is introduced.

Different unit identifications, rational solutions modulo constants, or
integer multiples produce allowed identity lifts. The canonical rational
component does not depend on those choices.

## 4. The correction cancels out of the component formula

Set t=[n]ell-b in the ACTUAL generic torus of J. It is represented by
the generic restriction of

    L^n tensor B^-1,

an invertible line bundle on C whose normalization restrictions are O(nD_v).
Apply the preceding attachment formula to this GENERICALLY TRIVIAL
normalization bundle. The factor n cancels, giving, for every graph cycle z,

    nu_J(ell)(z)
      =(1/e) sum_a z_a
         [d_(tail,j(a,tail))-d_(head,j(a,head))].      (6)

Here nu_J is the original-parameter rational Neron component, not the
uncorrected torus coordinate of an artificially split J=T times A.
The d_v in (6) are obtained from RAW multidegrees of L by (3). No generic
normalization triviality of L and no separately supplied b are needed.

Equivalently, form the attachment current J_v(z) as before and solve
M_v q_v(z)=J_v(z). Then

    nu_J(ell)(z)=(1/e) sum_v delta_v dot q_v(z).        (7)

These formulas are independent of the extension L and its model, because
Sections 2--4 identify them with the intrinsic corrected component. An
alteration which moves degrees to exceptional curves must also retain
the changed attachment components. Omitting either changes the calculation.

There is an exact zero statement for the full J presentation:

    nu_J(ell)=0 iff SOME positive tensor power of L_K
    admits an invertible extension with degree zero
    on EVERY special normalization component.         (8)

The tensor power is essential: finite abelian component torsion is allowed.
For necessity, use (5). Zero rational component makes the corrected torus
point t bounded. It has a unit-glued extension with trivial normalization
bundles. Tensor that extension with B to obtain a zero-multidegree extension
of L_K^n. Sufficiency is (2), followed by division by that positive integer
in the rational component group. For any GIVEN extension, the exact test
is (6), not the demand that its individual degrees vanish.

## 5. Pass through the ACTUAL Picard quotient and retain its kernel

Suppose the intrinsic original boundary presentation supplies

    psi:J -> G_f,
    psi(ell)=[a]sigma_f, a>0.                        (9)

This is an actual semiabelian morphism and a compatible rank-one lattice
map. Its existence and its multiplier must be checked in the supplied
Picard/NS construction. Neither is inferred from an auxiliary graph.

Functoriality of Neron models sends b in (5) to an identity-component
point psi(b). The point

    [na]sigma_f-psi(b)=psi([n]ell-b)

lies in the ACTUAL generic torus of G_f. Thus for an actual character chi
of that torus, use its pullback cycle psi^*chi and obtain

    nu_f(sigma_f)(chi)
      =(1/(a*e)) sum_edges (psi^*chi)_a
           [d_(tail,attachment)-d_(head,attachment)]. (10)

All line-bundle degree inputs in (10) precede the abelian correction.
The kernel of psi is retained by testing ONLY pullbacks of actual
characters. One may not replace this by all cut-graph directions, nor
by a scalar norm cancelling opposite branches. The multiplier a is
the actual lattice multiplier, not an assumed generic degree of a map.

This proves a direct raw-Picard-degree input formula for the existing
intrinsic Mathieu calculation whenever (9) has been supplied. It does
not prove that an arbitrary candidate cut symbol has such a compatible
Picard lift, or implement a universal producer for (9).

## 6. Existence of a raw invertible extension does not prove vanishing

The extension assumption on L can be achieved after marked modifications.
First extend each generic normalization line bundle to the regular Y_v:
take a rational section on the generic curve and close its divisor;
regularity makes that divisor Cartier. At each marked node the desired
generic identification has some valuation m_a relative to the resulting
R-line lattices. Blow up its chosen tail marked point and twist that
normalization extension by m_a times the exceptional P1. The strict
marked section now meets that P1, so the remaining node identification
is a unit. Glue to obtain an invertible raw extension on the expanded C.

This is the earlier construction applied to arbitrary generic normalization
classes, not only toric ones. Every resulting exceptional degree enters
(3),(6),(10). The Long control has raw degrees +1,-1 and coupling 1;
invertibility by itself remains insufficient. An ambient resolution tree
does not control these actual attachment periods.

## 7. The remaining full exceptional proof

The usable new input is the ACTUAL uncorrected Picard line bundle and its
intersection degrees, together with the genuine map (9). The balancing
construction supplies the identity lift and computes the intrinsic local
component via (10). It is no longer necessary to guess an abelian
correction or assume it can be omitted.

To close either exceptional Mathieu--Zhao case, one still must show that
unweighted purity forces the raw degree vectors to pair to zero with
EVERY actual attachment current in (10), on the full pure locus. The
dominant-root and capacity/automorphism constraints remain inputs to that
geometric implication. It has NOT been proved by this lemma.

Alternatively an actual pure exceptional phase with a nonzero output
would require its fixed original bad marker to complete the counterexample.
No such phase is produced here. The original volume, all relevant
components and branches, the completed family-transfer theorem and the
ordinary/exponential generating-function distinction remain intact.
Both FULL cases remain OPEN. External review remains deferred.
