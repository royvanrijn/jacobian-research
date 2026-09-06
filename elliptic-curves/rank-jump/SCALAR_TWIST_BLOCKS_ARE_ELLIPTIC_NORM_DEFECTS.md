# Scalar-twist blocks are simultaneous elliptic norm defects

The three retained production CT blocks share a concrete arithmetic
object: the norm map from **E(Q(i)) to E(Q)**. They certify failures of
its local-global principle of dimensions at least **8, 6, 6**. The
corresponding independent Sha classes on the negative twists all become
soluble over the **same quadratic field Q(i)**.

There is also a geometric boundary. Each production cubic has Galois group
S3, so the standard identification is the **only rational 2-torsion gluing
to a scalar quadratic twist**. Its quotient, with the gluing polarization,
is a Weil restriction that is geometrically a product, not a smooth
genus-two Jacobian. The nontrivial cyclic relabelling in the small
self-gluing control is unavailable here.

These facts sharpen the existing
[production comparison](PRODUCTION_TWIST_INCIDENCE_AND_SOLUBILITY.md).
They do not compute new ranks, a full norm group, or new CT entries.

## A common norm sequence

Put

    E: y² = x³ + A x + B,
    E−: y² = x³ + A x − B,
    F = Q(i).

The isomorphism over F is

    ι: E− → E,     (x,y) ↦ (−x, i y).

Complex conjugation σ sends ι to −ι. Let R be the Weil restriction
Res_(F/Q)(E_F). Thus R(Q)=E(F), and over F it becomes E×E. The norm,
or trace in the elliptic group law, is N(T)=T+σT. There is an exact
sequence of abelian varieties

    0 → E− → R →[N] E → 0.

The first map is induced by ι; over F its image is the anti-diagonal.
On rational points and cohomology this gives

    E(F) →[N] E(Q) →[∂] H¹(Q,E−) → H¹(F,E_F).

The final map is restriction followed by ι, by Shapiro's lemma. In
particular, every class in the image of ∂ is killed by restriction to F.

Under the standard identification E[2]≅E−[2], ∂(P) is precisely the
image of the rational Kummer class δ_E(P) in H¹(Q,E−). To check this,
choose a half H of P over Qbar. The diagonal pair (H,H) in R(Qbar)
has norm P. Its coboundary is (gH−H,gH−H); these differences are
2-torsion, so the pair lies on the anti-diagonal and gives the labelled
twist Kummer cocycle.

Consequently, for any rational P on E,

    δ_E(P) has a rational representative on E−
        ⇔ ∂(P)=0
        ⇔ P belongs to N E(F).

All statements use the elliptic **group norm**, not the multiplicative
norm of a number field.

## The local-global defect is a restriction kernel in Sha

Let E(Q)_locN consist of the rational points that are norms locally at
every place. At a split place the local algebra is Q_v×Q_v and the norm
is addition, so it is surjective. At other places use the quadratic
completion. The same connecting-map argument works locally. Equivalently,
the Kummer image of local norms is the intersection of the original and
twist Kummer images; see
[Morgan–Paterson, Lemma 4.3(ii)](https://arxiv.org/pdf/2011.04374).

Restricting the global exact sequence to locally trivial torsors gives

    E(Q)_locN / N E(F)
        ≅ ker(Sha(E−/Q) → Sha(E−/F)).

This identifies the **whole** norm defect with a simultaneous solubility
obstruction in one quadratic extension. It does not require finiteness
of Sha. The defect is a finite elementary 2-group because

    2E(Q) ⊆ N E(F):    N(P)=2P for P∈E(Q).

In particular N E(F) has the same free rank as E(Q). A large norm index
is not an additional Mordell–Weil rank, and surjectivity of the norm by
itself is not a rank-jump event.

## The production blocks and the low-gain comparison

The retained strict space W∩U is represented by rational points on E
and is locally trivial at the complete bad-place set. It remains Selmer
on E−, so every displayed class is represented by a point in E(Q)_locN.
The existing scalar-cup certificate is the CT form on these twist classes.

Choose the retained nondegenerate subspace V of that form. If a nonzero
v∈V were a global elliptic norm, its twist torsor would be rational and
would pair trivially with every Selmer class. Nondegeneracy on V excludes
this. Thus V injects into the norm defect and into Sha(E−)[2].

| Original fibre | Generic / independent witness / observed quotient | Retained strict dimension | Norm-defect dimension at least | Norm-defect order at least | Retained norm-kernel dimension at most |
|---|---:|---:|---:|---:|---:|
| A1/MW16-05, 307/206 | 16 / 25 / +9 | 10 | 8 | 256 | 2 |
| R17, −2300/843 | 17 / 24 / +7 | 8 | 6 | 64 | 2 |
| R17, −1561/3133 | 17 / 17 / observed 0 | 6 | 6 | 64 | 0 |

The last column bounds the kernel only inside the retained strict space;
an unpaired retained vector may still be obstructed by an uncomputed
class. Full norm-defect dimensions and full ranks remain UNKNOWN.
The observed-zero label is still censored, not an exact rank assertion.

Every nonzero torsor in the selected subspaces has period **and index**
two: it is nontrivial of order two, and has a point over F, a degree-two
extension. Thus the same F splits every one of the 255, 63, or 63 nonzero
torsor classes in the respective displayed subspace. These are counts
inside certified subspaces, not the full Sha groups.

The low control's entire six-dimensional retained strict space comes from
its generic subgroup. It already has the same six-dimensional obstruction
as the R17 +7 control. A common quadratic splitting field and a substantial
norm defect therefore do not distinguish those observed jumps.

There is no conflict with the
[independent halving-field result](HALVING_FIELDS_AND_BLOCKS.md).
Making an underlying genus-one torsor soluble over F does not mean that
the associated 2-cover has a point above the elliptic origin, or that
halves of all its rational Kummer representatives lie in F. Splitting a
torsor and adjoining prescribed halves are different requirements.

## Why scalar gluing cannot use the small genus-two construction

The rational isogeny

    Φ: E×E− → R,       (P,Q) ↦ P+ι(Q)

becomes over F the matrix

    M = [[1,1],[1,−1]]

after identifying E− with E. Its kernel is the graph of the standard
2-torsion identification; it has order four. The conjugation matrices
satisfy swap·M=M·diag(1,−1), verifying descent. Also

    M²=2I,     MᵀM=2I.

Hence the product principal polarization on R over Qbar pulls back to
twice the product polarization on E×E−. It is exactly the principal
polarization furnished by the gluing. In particular this polarized
quotient is geometrically decomposable and is not the polarized Jacobian
of a smooth genus-two curve. This is the reducible case of the usual
2-gluing construction; see
[Howe–Leprévost–Poonen, Proposition 3](https://math.mit.edu/~poonen/papers/large.pdf).

Could another rational identification avoid this? For any scalar quadratic
twist, all identifications form a torsor under Aut_GQ(E[2]). The production
cubics have group S3=GL2(F2). Its centralizer in GL2(F2) consists only of
the identity, so there is only one rational identification. It is induced
by the geometric twist isomorphism and is always reducible in this sense.

The coefficient-only part of the experiment independently certifies the
three S3 groups: irreducible reductions occur at 29,17,17 respectively,
and each rational cubic discriminant is nonsquare. Exhausting the six
invertible binary matrices gives:

| Galois image on E[2] | Equivariant invertible identifications with a scalar twist |
|---|---:|
| S3, production controls | 1 |
| C3, small cyclic control | 3 |

For the small cubic T³−11T²−14T−1, the discriminant is 163². Its two
nonidentity centralizing rotations are the extra gluing choices; since
the elliptic j-invariant is neither 0 nor 1728, they do not come from
geometric elliptic automorphisms. This accounts for the smooth genus-two
objects with transported labels in the
[small Selmer-lifting control](JACOBIAN_SELMER_LIFTS_CAN_BE_SHA.md).

The exclusion concerns this **2-gluing to scalar twists with its induced
polarization**. It does not exclude gluing to a different 2-congruent curve,
other isogeny degrees, or other principal polarizations in an isogeny
class. The earlier fixed-cubic deformation Jacobians have different
elliptic targets and are unaffected.

## What would turn norms into a soluble rational block?

The norm interpretation also gives an explicit sufficient construction.
Suppose rational points P₁,…,P_k on E have independent Kummer classes,
and a construction supplies T_i∈E(F) with N(T_i)=P_i. Then

    Q_i = ι⁻¹(2T_i−P_i)

is rational on E−: the argument of ι⁻¹ is anti-invariant under σ.
The pair (P_i,Q_i) comes from the dual isogeny of Φ, so its two labelled
Kummer classes agree. Therefore the Q_i are independent modulo two.
For the S3 controls, E−(Q)[2]=0, and this would imply rank E−(Q)≥k.

The condition is substantive: taking T=P rational only gives N(T)=2P
and a trivial mod-two class. Nor does local norm solubility suffice;
the table certifies whole subspaces where it fails globally.

To explain a jump in an original family using this construction with the
roles reversed, the P_i and their common trace construction must be
available before its exceptional points are supplied. Independence must
then be checked modulo the original generic subgroup, not just in the
whole elliptic group. This turn supplies no such prospective construction.

## Ranked implications

1. **Supported solubility structure:** the retained production obstructions
   belong to one elliptic norm defect and are killed by one quadratic
   field. The blocks are not unrelated local failures.
2. **Excluded shortcut:** the cyclic small control's nonstandard scalar
   genus-two gluing cannot transfer to these S3 production curves.
3. **Weak rank explanation:** common splitting fields, norm-defect size,
   and product isogenies alone do not distinguish high gain from observed
   zero. Norm saturation changes a finite index, not free rank.
4. **Missing implication:** a point-independent construction must provide
   several global trace preimages, or another proof of simultaneous
   rational solubility, and then certify their relative independence.
   Vanishing on a retained CT subspace is insufficient.
5. **For Agent 1:** the centralizer is an equation-only **incidence gate**
   for a specific geometric construction. The norm defect is a
   **solubility diagnostic**; its present numerical bounds use retrospective
   points. No **point-search visibility** feature or candidate score is
   produced, and no search changes follow.

## Frozen evidence and replay

The [protocol](QUADRATIC_NORM_BLOCK_PROTOCOL.json) permits three cubics,
primes at most 503, six binary invertible matrices, and ten seconds per
command. No global point, norm-equation, class-group or descent campaign
ran. The exact CT values and independence of original witnesses are
inherited from their pinned certificates; this calculation rechecks their
symplectic masks and derives the norm consequences.

- [Inputs and source hashes](../../artifacts/generated-results/elliptic-curves/rank_jump_quadratic_norm_block_inputs_v1.json)
- [Norm-block certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_quadratic_norm_blocks_v1.json)
- [Independent verification](../../artifacts/generated-results/elliptic-curves/rank_jump_quadratic_norm_block_verification_v1.json)

    python3 elliptic-curves/rank-jump/quadratic_norm_blocks.py check
    sage -python elliptic-curves/rank-jump/verify_quadratic_norm_blocks.py check
    python3 -m unittest discover -s elliptic-curves/rank-jump -p test_quadratic_norm_blocks.py

The independent replay verifies the general twist identity and the Weil
restriction's affine equations symbolically, the descent and polarization
matrices, the centralizers, the cubic certificates, and all retained CT
block transformations. The cohomology exact sequence is the written
mathematical argument; no global norm group is claimed to have been
enumerated. Active search files and mathematical-status entries are
unchanged.
