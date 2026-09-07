# Fourth division of the generic points has no new Selmer classes

The [uniform congruence-filtration proof](THE_GENERIC_TWO_POWER_DIVISION_TOWER_IS_CLOSED.md)
now extends this exclusion to the entire generic two-power division tower,
without requiring the actual higher torsion images.

On all sixteen frozen fibres, the field obtained by adjoining all fourth
preimages of the marked generic points contains exactly the generic Kummer
classes and the derivative class. The latter is outside Selmer. Thus this
larger, noncentral construction supplies **zero new Selmer directions**.
Its translation group over Q(E[4]) is maximal on every row, high and low.
No exceptional point enters the result.

There is a useful refinement to the required class-creation mechanism:
new standard S3 layers in a field tower are not sufficient. They must
survive in the elementary abelian quotient detected by quadratic
characters. In this fourth-division tower the new standard layers lie
in doubles of the translation group and contribute nothing to H1.

## The field being tested

Let E/Q have full S3 image on V=E[2], and let rational points P_1,...,P_m
give independent classes spanning G in E(Q)/2E(Q). Define

    M = Q(E[2]),       F4 = Q(E[4]),
    L1 = Q(E[2], [2]^-1 P_i for every i),
    L2 = Q(E[4], [4]^-1 P_i for every i).

All preimages are included; both L1 and L2 are Galois. These are finite
extensions attached to the supplied generic points, not searches for
rational preimages. The usual affine Kummer representation is described
in [Lombardo–Tronto, Explicit Kummer Theory](https://arxiv.org/abs/1909.05376).
The particular capacity theorem below is proved here.

The retained [four-torsion calculation](FOUR_DIVISION_RAMIFICATION_CANNOT_SUPPLY_THE_JUMP.md)
proves, without assuming a full mod-four image,

    ker(H1(Q,V) -> H1(F4,V)) = <beta>,
    beta = -disc(f) f'(theta).

On every panel row beta is nonzero and outside Sel2(E). Since G is
rational, beta is independent of G. This local exclusion is the arithmetic
input needed for the whole generic fourth-division result.

## Full translation image follows from mod-two independence

Restriction to F4 remains injective on G: its kernel is G intersect
<beta>=0. Over F4 the V action is trivial. The image of the m restricted
Kummer characters is therefore a subgroup H of V^m, stable under the S3
action induced by Gal(F4/Q).

The standard V is simple with End_S3(V)=F2. A proper submodule of the
semisimple module V^m admits a nonzero S3-equivariant linear functional
V^m -> V vanishing on it. Such a functional would be a nonzero F2-linear
combination of the m restricted characters that vanishes identically,
contradicting their independence. Hence

    Gal(F4 L1/F4) = V^m.

Let T=Gal(L2/F4), viewed by translation of the chosen fourth preimages
as a subgroup of E[4]^m=(Z/4)^(2m). Restricting to their doubles gives
the reduction of T modulo two, whose image is the full V^m just proved.
For any subgroup T of A=(Z/4)^n with T+2A=A, the quotient A/T satisfies
A/T=2(A/T), and iterating gives A/T=4(A/T)=0. Thus T=A. In particular,

    Gal(L2/F4) = (Z/4)^(2m),       [L2:F4] = 16^m.

This is a proof of maximal translation image for the actual fields.
It does not infer the linear mod-four image from a full-group control.
Write that actual linear image as Gamma4=Gal(F4/Q). The joint affine
image projects onto Gamma4 and contains every translation, so it is the
entire semidirect product

    Gal(L2/Q) = (Z/4)^(2m) semidirect Gamma4.

The displayed splitting is a consequence of containing the full
translation subgroup in the affine group; no arithmetic splitting
assumption is added.

## Exact cohomological capacity

For A=(Z/4)^(2m), acting trivially on V, a cocycle on A semidirect Gamma4
restricts to an equivariant homomorphism A -> V. Conversely each such
homomorphism phi extends to the cocycle (a,g) -> phi(a). Subtracting it
leaves a cocycle inflated from Gamma4. Consequently

    H1(A semidirect Gamma4,V)
      = H1(Gamma4,V) direct sum Hom_Gamma4(A,V).

Every homomorphism A -> V kills 2A. Its equivariance therefore depends
only on the full S3 action on A/2A=V^m. The Hom space has dimension m,
and the coordinate cocycles are exactly the original generic Kummer
classes. The retained torsion theorem gives H1(Gamma4,V)=<beta>. Inflation
into H1(Q,V) now proves

    ker(H1(Q,V) -> H1(L2,V)) = G + <beta>,
    Sel2(E) intersect ker(res_L2) = G.

The strict part is exactly G intersect U. A class g+beta cannot become
locally admissible by adding g in Selmer, because the local point image
is a subgroup. No Cassels–Tate calculation is needed for this exclusion.

The subgroup Gal(L2/F4 L1)=2A is itself a sum of m standard V modules
under the S3 quotient. Nevertheless a cocycle restricted to A is a
homomorphism into the exponent-two group V, so it vanishes on 2A.
The new layers are in the Frattini subgroup of the abelian 2-group A.
They cannot give new quadratic Kummer characters. This is stronger
than the previous exclusion of *central* governing extensions, whose
kernels have trivial S3 action.

## Paired consequences and exact checks

The [sixteen-row certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_generic_fourth_division_v1.json)
binds the generic-only panel and both the prior torsion certificate and
its independent verification. It revalidates their direct source hashes
and checks each local exclusion before applying the theorem. It does not
repeat their already completed local arithmetic or generic independence
proofs. Outcome labels play no part in the calculation.

| Marked group | Rows | Relative degree over F4 | Split Kummer dimension | Split Selmer dimension | Additional Selmer dimension |
|---|---|---:|---:|---:|---:|
| MW17 | 00–10,13,14 | 2^68 | 18 | 17 | 0 |
| MW16 | 11,12,15 | 2^64 | 17 | 16 | 0 |

This includes the fresh 103b2 high/low pair, both 11952 matched pairs,
the A1/MW16 high/low pair, and historic +12/+14 controls. The gains remain
the frozen lower-bound labels, with low-gain results censored. There is
no translation-image collapse at this level distinguishing the highs.

The [bounded protocol](GENERIC_FOURTH_DIVISION_PROTOCOL.json) also runs
exact Cayley-graph cocycle calculations on full affine control groups:

| Generic points in abstract control | Group order | dim Z1 | dim B1 | dim H1 |
|---:|---:|---:|---:|---:|
| 0 | 96 | 3 | 2 | 1 |
| 1 | 1536 | 4 | 2 | 2 |
| 2 | 24576 | 5 | 2 | 3 |

Explicit generic coordinate cocycles, the retained derivative cocycle,
and coboundaries span every computed solution space. Separately, exhaustive
subgroup enumeration in (Z/4)^2 confirms that its only subgroup with full
mod-two image is the whole group; enumeration of equivariant maps to V
gives exactly zero and reduction modulo two. The general proof, rather
than these small controls, covers arbitrary m and proper linear images.
All checks finish within the 30-second cap.

## What is still needed

1. **Incidence:** an independently constructed S-split unramified standard
   block that survives the relevant elementary abelian quotient, outside
   the inherited class pool. Fourth division of the generic points does
   not provide it. Generic ideal relations and S-unit corrections remain
   possible sources; their missing relations have not been constructed here.
2. **Representation change:** actual higher torsion images remain UNKNOWN.
   This theorem is at level four. It is not a whole-tower exclusion, and
   full abstract mod-eight/mod-sixteen group calculations cannot replace
   those missing actual-image certificates.
3. **Solubility:** even a new surviving standard block still needs a
   rationality argument beyond local admissibility and vanishing CT.
   No new rational points or CT values were computed.

The structural chain has been narrowed, not completed: extra layers in
a generic division tower do not explain the exceptional dimensions here.
No search score or selector follows. Only rank-jump files are changed.

```sh
timeout 30 python3 elliptic-curves/rank-jump/generic_fourth_division.py check
```
