# The first singular genus-one layer is a finite tangent construction

**Written geometric theorem, with a symbolic-premise checker.** Fix an elliptic
K3 over a number field k with exactly24 geometric I1 fibres. There are only
finitely many k-defined geometrically integral bisections of arithmetic genus2
and normalization genus1, modulo translation by E(k(t)). Consequently these
bisections realize only finitely many literal quadratic extensions of k(t).
This is not an effective list or an emptiness theorem.

The result covers the complete first singular genus-one layer, not merely
the [frozen height-six chord nets](R17_ONE_NODE_CORRELATED_COVER_GATE_2026-09-12.md).
It supplies neither a matching pair nor a rational nontorsion point on any
cover. The [MW17 correlated-gain target](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
remains open. Higher arithmetic genus is outside this finiteness theorem.

## 1. Anti-trace height from the image's arithmetic genus

Let B be any geometrically integral bisection, C its normalization, and Q the
section obtained on the degree-two base change to C. Assume its branch values
are smooth parent fibres. Then chi(O of the pullback)=4 and all its fibres
are irreducible. Put T=Q+sigma(Q), V=Q-sigma(Q), z=B.O, and a=p_a(B).

In the Shioda decomposition, F^2=0, O^2=-2, O.F=1 and v_T is orthogonal to F,O
with v_T^2=-h(T). The generic degree-two divisor B has elliptic sum T, so

```
B=2O+(z+4)F+v_T,
2a-2=B^2=4z+8-h(T),
h(T)=4z+10-2a.
```

Pullback to the normalization preserves the intersection multiplicities with
O, including at singular points of B. Thus h(Q)=8+2z. Deck invariance gives

```
h(V)=4h(Q)-2h(T)=4(a+3).                         (1)
```

The two anti-invariant directions arising from images of arithmetic genera
a,b are independent if (a+3)/(b+3) is not a square in Q. Indeed a rational
linear dependence of nonzero directions forces their height ratio to be a
rational square. They are orthogonal to inherited sections, so this certifies
independence modulo the inherited subgroup. A square ratio is inconclusive.
The images must first be proved to have the same literal quadratic cover.

This includes the previous heights16 and20 for smooth genus-one and
arithmetic-genus-two images. It retains singular-image multiplicities and
places no coordinate denominator bound. Formula(1) does not apply unchanged
when the cover branches at singular parent fibres and reducible-fibre height
corrections appear.

## 2. The quotient attached to an arithmetic-genus-two image

Now assume a=2. The class calculation on the parent, which does not need the
smooth-branch assumption, gives h(T)=4z+6. This is2 modulo4. The height lattice
of a24-I1 K3 is even integral, so T cannot be twice a geometric section.
Its halving curve D_T is therefore connected of genus9 with full S4 monodromy,
by the [24-I1 halving argument](JOINT_HALVING_EXCEPTIONAL_BRANCH_LOCUS_2026-09-14.md).
No arithmetic saturation of a chosen finite subgroup may replace this
geometric primitivity check.

Consider eta(P)=T-P. It extends to an involution on the K3 and preserves B.
Its quotient Y is smooth and ruled over the original t-line, as in the
[ruled-quotient proof](Q80_HALVING_PENCIL_GEOMETRY_2026-09-15.md): at an I1 node
eta exchanges the local branches, whose quotient is smooth. The image of O
supplies a k-section, so Y is a split Hirzebruch surface F_n over k. Write C0
for its negative section and F for a ruling fibre.

The fixed curve maps isomorphically to the branch curve R=D_T. There are no
other fixed curves: none can be vertical, and generic fixed points are exactly
the connected four-half scheme. Hence R is smooth and irreducible, and

```
K_Y=-2C0-(n+2)F,       R=-2K_Y=4C0+(2n+4)F.
```

Irreducibility forces R.C0=4-2n>=0, so n<=2.

The image C of B is a section of Y -> P1_t. Its inverse image is precisely B,
so B^2=2C^2. Since B^2=2, a section C=C0+kF satisfies

```
1=C^2=-n+2k.
```

The only possibility with n=0,1,2 is n=1,k=1. Thus all the bisections in this
trace class lie over members of |C0+F| on F1.

## 3. Tangents and the finite arithmetic parameter set

Blow down C0 to p in P2. The curve R becomes a plane sextic of geometric
genus9. Members of |C0+F| which are sections correspond to lines not through
p. The lines through p give a reducible divisor C0+F_t on F1 and do not
produce an integral bisection of this class.

Along an eligible line, a contact of multiplicity m with R gives local equation
v^2=u^m times a unit. Its normalization delta is floor(m/2). Since B has
arithmetic genus2 and normalization genus1, its total geometric delta is1.
There is therefore exactly one contact of order2 or3 and no other repeated
contact. The point of contact is k-rational: the divisor of repeated contacts
is Galois-invariant and consists of this one geometric point. It lies on the
smooth curve R over k, not merely at a singular point of its plane model.

Conversely, a k-point of R away from C0 determines its unique tangent line.
When that line avoids p, has precisely the stated contact pattern, and gives
an integral inverse image, it produces one of the required bisections. These
open checks must be performed; a tritangent line can instead give a split
inherited divisor. Rational points over p must not be counted as ordinary
eligible tangent points. Their tangent lines through p are excluded.

Thus the eligible genus-one normalizations for fixed trace T inject into
R(k)=D_T(k). This set is finite by Faltings because g(D_T)=9. This use of
Faltings does not enumerate the points or give a height bound.

Finally E(k(t)) is finitely generated, and translation by S changes T to T+2S.
Choose the finitely many representatives of E(k(t))/2E(k(t)). Translating B
preserves its arithmetic genus, normalization and map to the original base.
For each representative the preceding argument gives finitely many members.
This proves finiteness modulo inherited translations, and finiteness of their
actual quadratic function fields. The latter includes constant squareclasses:
each member has its own fixed inverse-image equation, not an arbitrary
constant twist with the same branch divisor.

## 4. Consequences and limits for the active construction

The first singular layer can be approached through rational points on the
halving curves and their tangent lines. This replaces an unbounded search in
the coefficients and node positions of a fixed net by a precise arithmetic
curve problem. It does not authorize an enumeration of all trace parities.
The retained genus-nine curve controls the point field, tangency and residual
quartic; solving it would still require checking the actual cover and its
rational points.

To obtain two gains, one may try to match a resulting quartic with a cover
from a different image genus. Formula(1) makes independence automatic when
the genus ratio is nonsquare. Same-genus pairs still require their own
independence proof. No such match is supplied here, and the completed frozen
pair exclusions remain in force.

The [checker](scripts/verify_singular_bisection_geometry.py) verifies the
symbolic height identity, the three possible ruled-surface indices, the
sextic/genus intersection numbers and contact delta table. Its
[certificate](../artifacts/generated-results/elkies-k3-singular-bisection-geometry-v1/result.json)
is not an arithmetic point list or an independent proof of the geometry.
Replay:

```
.venv/bin/python research/elkies-k3/scripts/verify_singular_bisection_geometry.py
```

The quotient, tangent correspondence, finite generation and Faltings step are
written mathematics, without formal verification or external review. The
intersection methods are standard; no foundational novelty claim is made.
