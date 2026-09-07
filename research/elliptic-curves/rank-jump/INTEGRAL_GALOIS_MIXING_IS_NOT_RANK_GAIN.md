# Integral Galois mixing measures common rational classes, not rank gain

The closed rank-one/rank-three rational/Sha control has **no integral
invariant/anti-invariant mixing** over Q(i). Its two-dimensional solubility
switch therefore needs no such mixing. Conversely, a prescribed point
(i,1+2i) gives an S3 example with one mixed block while both rational
elliptic curves have exact rank one.

The relevant invariant is now identified exactly: it is the dimension
of the intersection of the two **full rational Kummer images**, under
the standard scalar-twist labels. It measures a finite index and the
number of regular integral Galois summands, not additional free rank.

This develops the [elliptic norm interpretation](SCALAR_TWIST_BLOCKS_ARE_ELLIPTIC_NORM_DEFECTS.md)
without a production point search or another class-group calculation.

## The integral lattice and its common Kummer image

Let F/Q be quadratic with involution sigma, E− its associated twist,
and ι:E−_F→E_F the twist isomorphism. Assume E(F)[2]=0. This holds
for both small controls and the retained production cubics: a root of
an irreducible rational cubic cannot lie in a quadratic field.

Write M=E(F) tensor Z_2 and let M+ and M− be its invariant and
anti-invariant submodules. They are the completions of E(Q) and
ι(E−(Q)). Their ranks are r+ and r−. Put

    h = dim_F2 M/(M+ + M−).

The quotient is killed by two, since 2T=(T+sigma T)+(T−sigma T).
The trace map gives an isomorphism

    M/(M+ + M−) ≅ N E(F)/2E(Q).

For the kernel assertion, if N(T)=2P with P rational, then T−P is
anti-invariant. The converse is immediate. Odd torsion does not affect
these 2-primary statements.

Let G+ and G− denote the full rational 2-Kummer images of E and E−,
in the common H¹(Q,E[2]). The norm sequence from the preceding note gives

    N E(F)/2E(Q) ≅ G+ intersect G−.

In particular,

    h = dim(G+ intersect G−),
    [E(F) : E(Q)+ι(E−(Q))] = 2^h,
    [E(Q) : N E(F)] = 2^(r+−h).

The middle quotient is killed by two, so the displayed index is the
actual finite index, not just its 2-part. No rank or saturation assertion
about a finite displayed point list is substituted for the full groups.

## Exact decomposition into three types of integral summands

There is a Z_2[Gal(F/Q)]-module decomposition

    M ≅ Z_2^(r+−h)  ⊕  Z_2(sign)^(r−−h)  ⊕  Z_2[C2]^h.

Here a regular summand Z_2[C2] has a basis exchanged by sigma. Its
invariant and anti-invariant generators span a sublattice of index two.

For a direct proof, embed M between M+⊕M− and its half-lattice.
The image modulo M+⊕M− projects injectively to both M+/2M+ and
M−/2M−: a class projecting to zero on one side would already belong
to the other eigenspace. Hence it is the graph of an isomorphism between
two h-dimensional binary subspaces. Choose bases lifting those subspaces.
For each matched pair a_i,b_i the extra generator is (a_i+b_i)/2, whose
conjugate is (a_i−b_i)/2. They give one exchanged pair. Unmatched basis
vectors give trivial or sign summands. This proves the decomposition.

Equivalently,

    h = rank_F2(sigma−1 acting on E(F)/2E(F)).

This formula is for the rational Mordell–Weil image. Replacing it with
the action on Sel2(E/F) may add Sha contributions and is not an equality
without further evidence. Integral splitting and Selmer Galois action
are studied, in a different statistical setting, by
[Morgan–Paterson, §1.2](https://arxiv.org/pdf/2011.04374).
Their density theorem for full rational 2-torsion is not applied to our
S3 controls.

The companion twist interchanges the trivial and sign summands and
leaves the regular summands unchanged. In particular h is the same on
both sides of the scalar twist comparison. The rank difference is

    r+−r− = (number of trivial summands) − (number of sign summands),

with the mixed summands cancelling. A larger common rational Kummer
intersection alone cannot explain which side has greater rank.

## Closed rational/Sha control: a switch with h=0

Use the retained exact pair

    E0: y²=x³−11x²−14x−1,       rank 1,
    E1: y²=x³+11x²−14x+1,       rank 3.

In coordinates (beta0,beta1,u0,u1), with u0=−1−theta and
theta=u0*u1 modulo squares, the previously proved full spaces are

    G0 = <u0>,                 S0 = <beta0,beta1,u0>,
    G1 = <beta0,beta1,u0*u1>,    S1 = G1.

Thus G0 intersect G1=0, while S0 intersect S1 is the two-dimensional
strict block. This uses **standard scalar labels**. The order-three
transport used for the genus-two control is a different identification
and must not be inserted here.

There is a direct real-sign check of the zero rational intersection.
The roots of the common cubic lie in (−2,−1),(−1,0),(12,13).
The rational generator on E0 has Kummer sign pattern (+,−,−).
The three rational generators on E1 have patterns (−,−,+),
(−,−,+),(+,+,+), respectively. Their sign image does not contain the
first pattern. Since G0 is one-dimensional, its intersection with G1
is zero.

Consequently

    E0(F) tensor Z_2 ≅ Z_2 ⊕ Z_2(sign)^3,
    E0(F) = E0(Q)+ι(E1(Q)),

and the full group has rank four. For the reverse twist the trivial
and sign multiplicities swap. The norm indices and defects are:

| Norm target | Full norm index | Local norm codimension | Local-global norm-defect dimension |
|---|---:|---:|---:|
| E0(Q), rank 1 | 2 | 1 | 0 |
| E1(Q), rank 3 | 8 | 1 | 2 |

The two-dimensional defect in the second row maps to the already known
Sha(E0)[2]. Both strict classes are rational on E1 and Sha on E0, yet
h=0 throughout. This is a counterexample to requiring nonzero integral
mixing for a simultaneous rational/Sha switch. It is a controlled twist
rank difference, not a new specialization result for MW17 or MW16.

## Constructed control: h=1 with ranks only 1 and 1

Prescribe T=(i,1+2i), then solve the real and imaginary parts of
y²=x³+A x+B. They force A=5 and B=−3. This is one fixed
construction, not a coefficient sweep. Put

    E: y²=x³+5x−3,       E−: y²=x³+5x+3.

Both cubics are irreducible, with discriminant −743, so their 2-torsion
has S3 action. Exact group arithmetic gives

    T+sigma T = P=(4,−9),
    T−sigma T = ι(Q),     Q=(1,−3) on E−,
    2T = P+ι(Q).

The common Kummer class is also explicit. In K=Q(theta), where
theta³+5theta−3=0,

    (4−theta)(1+theta) = (theta²+2)².

At the good prime 3, both reductions have exactly the origin and
(0,0),(1,0),(2,0). Every double is the origin. Both P and Q reduce
to (1,0), so neither is globally divisible by two. Since there is no
rational 2-torsion, their nonzero Kummer classes also prove positive rank.

The two predeclared effort-zero PARI descents return

    [1,1,0,[(4,9)]],       [1,1,0,[(1,3)]].

The matching upper bounds close both ranks at one. PARI documents its
2-Selmer, rational 2-torsion and CT-based upper bounds as unconditional;
see [ellrank documentation](https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html#ellrank).
No analytic rank is used. The point lower bounds and trace identities
have a separate CAS-free replay; it does not independently redo the two
descent upper bounds.

The common nonzero class forces h≥1, and the exact ranks force h≤1.
Therefore

    E(F) tensor Z_2 ≅ Z_2[C2],
    [E(F):E(Q)+ι(E−(Q))]=2.

Both elliptic norm maps are surjective; their local-global defects vanish.
The point T supplies a representative of the nontrivial index-two coset.
This does not claim T and sigma T form an integral basis at every odd
prime. It is an exact 2-primary mixing statement with no large rank.

The scalar gluing is still the geometrically product Weil restriction.
The S3 obstruction to a smooth genus-two 2-gluing from the preceding note
does not prevent a mixed rational point on that abelian surface.

## Production implications without new twist searches

Let d be the codimension of globally rational points satisfying all local
norm conditions, and nu the local-global norm-defect dimension. The norm
quotients give

    r+−h = d+nu.

For the three production rows, the marked generic points already span
the entire relevant local quotient. Hence d is exactly the previously
certified local change 8,6,3, rather than an estimate from incomplete
local images. The CT norm blocks give nu≥8,6,6. Thus:

| Original control | d | nu at least | Trivial integral summands at least | Bound on h |
|---|---:|---:|---:|---|
| A1/MW16-05, 307/206, observed +9 | 8 | 8 | 16 | h≤rank(E/Q)−16 |
| R17, −2300/843, observed +7 | 6 | 6 | 12 | h≤rank(E/Q)−12 |
| R17, −1561/3133, observed 0 | 3 | 6 | 9 | h≤rank(E/Q)−9 |

If the known witness ranks 25,24,17 happened to be the full ranks, the
last bounds would be 9,12,8. That condition is not proved. The actual h
values, anti-invariant ranks, and full norm defects remain UNKNOWN.
These are abstract summand counts, not claims that designated generic
points form particular summands in a canonical basis.

## Ranked consequences and the next missing implication

1. **Exact solubility interpretation:** integral mixing equals the common
   full rational Kummer dimension. A trace preimage constructs a rational
   representative on the companion twist with matching class.
2. **Disproved necessity:** a simultaneous rational/Sha switch can occur
   with no mixed summands, as the closed two-class control demonstrates.
3. **Disproved sufficiency for high rank:** one mixed summand and a
   globally surjective norm coexist with ranks 1 and 1. Mixing measures
   an index and consumes existing invariant and anti-invariant dimensions.
4. **Missing production computation:** determine common rational images
   or construct several global trace preimages without exceptional-point
   input. A Selmer intersection or its Galois action alone leaves the
   rational/Sha distinction unresolved.
5. **For Agent 1:** this is a **solubility and saturation diagnostic**,
   not an **incidence predictor** or a **point-search visibility** score.
   No large Galois module should be counted as new rational rank before
   separating its rational image, its Sha quotient, and its generic
   subgroup. Current search choices remain unchanged.

## Evidence and replay

The [protocol](INTEGRAL_QUADRATIC_MIXING_PROTOCOL.json) fixes the two
small curves by the prescribed complex point, two effort-zero descents,
and a fifteen-second worker cap. The worker completes in under a second;
its raw checkpoint and log are retained under the ignored rank-jump
local-artifact directory.

- [Ranks, trace witnesses, complete small intersections, and production bounds](../../artifacts/generated-results/elliptic-curves/rank_jump_integral_quadratic_mixing_v1.json)
- [CAS-free witness and accounting verification](../../artifacts/generated-results/elliptic-curves/rank_jump_integral_quadratic_mixing_verification_v1.json)

    sage -python elliptic-curves/rank-jump/integral_quadratic_mixing.py check
    python3 elliptic-curves/rank-jump/verify_integral_quadratic_mixing.py check
    python3 -m unittest discover -s elliptic-curves/rank-jump -p test_integral_quadratic_mixing.py

The old full small-control spaces and production CT/local certificates are
explicit dependencies. Only the two prescribed low-rank curves received
new descent calls. No production point, parameter, class group, active
search file, or mathematical-status entry was changed.
