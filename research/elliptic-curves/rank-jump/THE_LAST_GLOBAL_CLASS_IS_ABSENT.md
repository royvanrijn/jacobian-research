# The last global class is absent on the published R17 root curve

The remaining global incidence gate is now closed:

\[
 \boxed{J_C[2](\mathbf Q)=G,\qquad \dim G=17.}
\tag{1}
\]

Here C is the fixed published R17 cubic root curve, and G is generated
by the seventeen marked generic sections. The residual degree-four
étale scheme Ω defined in
[the preceding theorem](A_QUARTIC_GOVERNS_THE_LAST_GLOBAL_CLASS.md)
has **Galois group S4 and no rational point**. Its Frobenius at 131 is a
4-cycle. An explicit coefficient polynomial for Ω remains uncomputed;
it is unnecessary for this conclusion.

This eliminates the proposed extra common class in the retained native
construction. For **any subset of k of the 37 retained native supports**,
with their original scalar choices, the full multiquadratic pullback has
arithmetic generic rank **exactly 17+k**. Every singleton character has
rank one; every character involving two or more distinct supports has
rank zero. Arbitrary nonzero rational scalar choices give the upper
bound 17+k, without an assertion that the singleton bound is attained.

This is a fixed-family **incidence exclusion**, not a rank predictor in t.
It does not rule out additional rank at isolated rational points of any
of those parameter covers. It does not compute the new specialized
class groups or explain their growth. It does show that further attempts
to make the conjectural common class soluble on these native twists
cannot work: the class itself is absent.

## Retained evidence and new arithmetic

All new arithmetic uses only the published A(t), B(t), and seventeen
generic sections. No exceptional coordinates, specialization parameters,
gain labels, new points, or prospective selections enter it.

The original K3 surface at p=131 was already independently counted over
F131 and F131² in the
[Frobenius certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_native_twist_frobenius_v1.json)
and its
[independent replay](../../artifacts/generated-results/elliptic-curves/rank_jump_native_twist_frobenius_verification_v2.json).
The original-surface columns of the retained 8779-orbit ledger give
parabolic H1 traces 1884 and 319520. The native-twist columns play no role
here. Those fiber counts are not rerun.

The new certificate verifies the following exact data:

| Quantity | Result |
|---|---:|
| Discriminant degree, over Q and F131 | 24 |
| Repeated discriminant roots, over either field | 0 |
| Infinity fibre | smooth |
| Generic height Gram determinant, over either field | 948 |
| Smooth base values avoiding all marked y zeros | 119 |
| Degree-one cubic-root places at those values | 91 |
| Rank of their generic Kummer character matrix | 17 |
| Irreducible cubic specializations among those values | 52 |

For a pair of integral sections P_i=(x_i,y_i), P_j=(x_j,y_j), the
intersection count is

    (P_i.P_j) = deg gcd(x_i-x_j,y_i-y_j)
                + min(4-deg(x_i-x_j),6-deg(y_i-y_j)).

The second term uses the infinity coordinates u=1/t, X=u⁴x, Y=u⁶y.
It is the common vanishing order of the coordinate differences there.
Each section is disjoint from O. All singular fibres are I1, so no
reducible-fibre correction occurs. The height diagonal is 4, and the
off-diagonal entry is 2-(P_i.P_j). This is the standard intersection
height formula; see
[Schütt–Shioda, §11.8 and equation (22)](https://arxiv.org/pdf/0907.0298).
The full matrices over Q and F131 agree, not only their determinants.

At a retained place (t=a,theta=b), the character in column i is the
Legendre symbol of x_i(a)-b. These are evaluation characters of the
function-field Kummer classes x_i-theta. Their matrix has rank 17,
with seventeen independent rows explicitly retained. A function that is
a square has square residue at every such place. Hence the marked
sections are independent in E(F131(t))/2E(F131(t)). This is an exact
two-saturation witness once rank 17 is established below, not an
inference from a finite point-search miss.

## The full residual Frobenius polynomial from two moments

The nineteen divisor classes from the fibre, zero section and seventeen
generic sections form a nondegenerate subspace in H²: their determinant
has absolute value 948. Frobenius acts there by p. Its orthogonal
complement has dimension three and is a weight-two orthogonal
representation. Its eigenvalues have the form

    epsilon*p, alpha, p²/alpha,       epsilon in {+1,-1}.

Subtracting the seventeen known parabolic directions from the ledger
traces gives

    T1 = 1884 - 17*131 = -343,
    T2 = 319520 - 17*131² = 27783.

The necessary identity

    (T1-epsilon*p)² = T2+p²

holds only for epsilon=-1. Thus alpha+p²/alpha=-212, and the residual
characteristic polynomial in the eigenvalue variable is exactly

\[
 (X+131)(X^2+212X+17161).
\tag{2}
\]

Neither factor vanishes at X=131. Therefore the full H² p-eigenspace
has dimension nineteen and is spanned by the displayed divisor classes.
In particular Tate's divisor statement holds for this reduction, its
arithmetic Néron–Severi rank is nineteen, and E(F131(t)) has rank seventeen.
This argument does not assume an exact rank from a lower bound, or
identify geometric Picard rank with arithmetic Picard rank.

Let I be the index of the marked subgroup in E(F131(t)). There is no
torsion: on this nonconstant elliptic K3 with only irreducible fibres,
every nonzero section has height 4+2(P.O)>0. The two-saturation witness
proves I odd. The full height lattice is integral, since all fibres are
irreducible. Its determinant is 948/I², so I² divides 948=4*3*79. There
is no nontrivial odd square divisor. Consequently **I=1**, and the full
arithmetic Néron–Severi determinant has absolute value 948 as well.

## Artin–Tate removes finite-field Sha[2]

In the reciprocal variable T, equation (2) gives

    P2(T) = (1-131T)^19 (1+131T)(1+212T+17161T²).

Its Artin–Tate factor is

\[
 \lim_{T\to1/131}\frac{P_2(T)}{(1-131T)^{19}}
 =2\left(2+\frac{212}{131}\right)=\frac{948}{131}.
\tag{3}
\]

For a K3, the Artin–Tate formula is

    absolute leading factor = #Br(S) * abs(disc NS(S)) / p.

The Néron–Severi group is torsion-free and chi(O_S)-1=1. Tate's statement
has just been proved here from the explicit cycles and (2), so use of
Artin–Tate introduces no unproved conjecture. See
[Milne, the 1975 theorem and its stated consequences](https://www.jmilne.org/math/articles/1975a.html)
and the formula in
[Milne–Ramachandran, page 2](https://www.jmilne.org/math/articles/2013d.pdf).
Equations (3) and the determinant 948 give #Br(S)=1; only its trivial
two-primary part is needed. The elliptic surface has a section, so its
Brauer group identifies with the function-field Tate–Shafarevich group;
see
[Milne's comparison theorem](https://www.jmilne.org/math/articles/1982f.html).
The Kummer exact sequence therefore gives

\[
 \dim\operatorname{Sel}_2(E/\mathbf F_{131}(t))=17.
\tag{4}
\]

## Why this Selmer group is the root Jacobian's rational two-torsion

This step concerns finite residue fields and must not be replaced by the
weaker geometric-local definition of the earlier global pool. Put
k=F131, F=k(t), and K=k(C). The cubic descent identifies H¹(F,E[2]) with
the norm-square kernel in K*/K*². Let U be its subgroup represented by
functions with even valuation at every closed point of C.

For this 24-I1 model, the local Kummer condition at every place v of F
is precisely the even-valuation norm-square subgroup in the local cubic
algebra. Here is the local verification.

* At a good place, the cubic algebra is unramified. If it has s field
  factors, its unit squareclasses have dimension s and the norm to the
  base unit squareclasses is surjective. The norm kernel has dimension
  s-1. The local point image lies in this unit subgroup and has that same
  dimension: multiplication by two is invertible on the formal group,
  and the finite residue elliptic group has E/2 dimension dim E[2]=s-1.
* At an I1 place, the cubic algebra is a base-field factor times a
  ramified quadratic field. The even-valuation squareclasses have
  dimension two. Norm is the identity on the base unit squareclass and
  trivial on the ramified quadratic unit squareclass, so its kernel has
  dimension one. Every local point maps into it. To check this inclusion,
  an integral point cannot reduce to the node: over the maximal
  unramified extension the surface is locally uv=pi, so a section cannot
  have both u and v divisible by pi. At a simple two-torsion reduction,
  the equation makes the corresponding factor valuation even; the
  other factors are units. Near O the pole orders of x are even.
  The local point quotient has dimension one, because the component
  group is trivial and the residue identity component is a one-dimensional
  torus, whose finite group has even order in odd characteristic.
* Infinity is good after x is replaced by u⁴x. This coordinate change
  multiplies each Kummer representative by a square, so the first case
  applies without an extra local condition.

These exhaust the places. Thus Sel2(E/F)=U. Finally, write div(a)=2D
for a in U. The map a -> [D] identifies U with J_C[2](k). Its kernel
consists of constants, and the norm-square condition kills that kernel
because the cover has odd degree three. Conversely any rational
two-torsion line bundle has a rational divisor over a finite field
(Br(k)=0); choosing its square-trivializing function gives norm c*h(t)².
Multiplying the function by c makes the norm a square. This proves

\[
 \operatorname{Sel}_2(E/k(t))\simeq J_C[2](k).
\tag{5}
\]

The degree-24 squarefree discriminant and smooth infinity also certify
good reduction of the genus-ten root curve at 131: the finite
ramification is simple, and the weighted infinity cubic has distinct
roots. Prime-to-131 torsion therefore specializes injectively. Combining
(4)-(5) with the seventeen rational generic classes proves (1).

## The residual quartic and native capacity now have exact answers

The earlier theta-pairing theorem gives the Galois-module description

    J_C[2] = F2^16 direct-sum F2[Ω].

At 131 the fixed dimension is seventeen. The permutation on Ω therefore
has one orbit and is a 4-cycle. Together with the already certified
3-cycle at 211 and transposition at complex conjugation, this forces
Gal(Ω)=S4. In particular Ω has no rational point.

Write L for the original Q(t) global pool. The root-curve theorem gives
L=J_C[2](Q)=G. For each of the retained 37 branch polynomials q_i, the
existing exact branch computation gives

    G intersect Z_i = <w_i>,       w_i nonzero,

and the 37 lines are pairwise distinct. Every rational point class on
the d-twist, with d a scalar times a product of those q_i, belongs to
the intersection of its branch kernels. Thus a singleton has rank at
most one and every multiple-support twist has rank zero. There is no
rational two-torsion on these twists. The known native sections attain
each original singleton bound. Character decomposition proves

\[
 \boxed{\operatorname{rank}E\bigl(\mathbf Q(t)(\sqrt{q_i}:i\in I)\bigr)
          =17+|I|}
\tag{6}
\]

with the retained native scalars absorbed into q_i. The disjoint simple
branch supports guarantee degree 2^|I|. This uses the branch certificates
already linked in
[the common-class theorem](ONE_COMMON_CLASS_FOR_A_LARGE_NATIVE_BLOCK.md);
their generic section identities and local proofs remain unchanged.

## Consequence for the class-creation programme

The ranked structural lessons are now:

1. **Required specialized mechanism:** the successful fibres need new
   cubic-field unramified, S-split quadratic extensions, equivalently
   new standard S3 class-group blocks. The existing strict-capacity
   bounds on the frozen panel still force this. These classes have not
   been constructed point-independently in this calculation.
2. **Closed proposed mechanism:** hidden generic mixed-character blocks
   on the retained native covers do not exist. The former common-h
   solubility question is empty on this published model. The individual
   native directions remain valid; their number is not evidence for a
   larger hidden generic block.
3. **Separate second stage:** the six-direction rational/Sha switch
   remains a solubility control. It cannot supply the absent incidence
   dimensions, and the present finite-field Sha vanishing is not a CT
   computation for any rational specialization.
4. **Missing implication:** find an equation-only condition on t0 that
   constructs several independent strict classes outside the specialized
   global pool. Independence and S-splitting must be certified before
   asking whether their elliptic covers have rational points. There is
   still no such prospective selector here for Agent 1 to adopt.

This theorem is asserted for the fixed published R17 presentation and
the specified native supports. It does not silently sharpen every
other chart's pool bound or relabel the panel's censored low controls.
No active search, worker setting, candidate population or MATH_STATUS
entry changes.

## Reproduction and scope

The
[frozen protocol](RESIDUAL_QUARTIC_AT_131_PROTOCOL.json)
discloses the trace observation before the saturation computation. The
[new certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_residual_quartic_at_131_v1.json)
binds its script, generic input and retained counts by SHA-256. The
[independent verifier](../../artifacts/generated-results/elliptic-curves/rank_jump_residual_quartic_at_131_verification_v1.json)
uses SymPy polynomial arithmetic, exhaustive integer root enumeration,
Euler's criterion, and separate bit elimination. It recomputes all 272
height intersections, all 91 characters, the determinant and the trace
algebra. The original finite-field point recount remains supplied by
the retained independent 8779-orbit replay, not by a duplicate run.

```sh
timeout 60 sage -python elliptic-curves/rank-jump/residual_quartic_at_131.py check
timeout 60 sage -python elliptic-curves/rank-jump/verify_residual_quartic_at_131.py check
```
