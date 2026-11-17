# Primitive end classes and canonical contact constraints exclude triple torus packets

PF2RULED102, 16 November 2026. This note closes an actual higher-degree
boundary class, using the full marked peripheral data rather than an
abstract knot-cover substitution.

**Six-sheet consequence.** Let a complex plane Keller map have generic
field degree six. Suppose its actual finite normalization N is smooth,
and N minus the marked whole source plane consists of exactly one
embedded affine line D, of transverse ramification index three. Then
its branch-image knot at infinity is neither an unknot nor a torus knot.
All torus parameters and all finite singularities of the image are
covered; this is not a bounded search or a zero-excess assertion.

The general results used to obtain this exclusion hold in every field
degree and at every positive ramification index:

* The actual omitted end knot generates finite cyclic end homology
  Z/h; its order h has parity e-1, where e is the index.
* The boundary contact plane field has Gompf invariant
  theta=-1-(e-1)^2/h.
* If the oriented end is the boundary Y_h of the Euler-minus-h disk
  bundle over S2, h>=2, then h>=e+1. In particular an RP3 end forces
  e=1 and is impossible for an actual nonautomorphic one-prime Keller
  normalization.
* For a torus boundary image, the actual monodromy is centre-free,
  every peripheral longitudinal degree is one, and h divides the
  field degree d. Hence odd d and odd e are incompatible.

Smoothness, the sole embedded omitted line and the genuine finite
normalization are explicit hypotheses. No conclusion for arbitrary
singular or multiple-prime packets is asserted. The global
[Laurent-key/finite-algebra implication](LAURENT_FINITE_ALGEBRA_FINISH_LINE.md)
remains open. Published [Domrina2000 / HC4MYKR0](../../hessian/hc4-triangular-graphs/QUADRATIC_SOURCE_RECONSTRUCTION_AND_QUARTIC_RETIREMENT.md)
already settles degree four; it is not reopened here.

## 1. The actual source has a primitive end class

Write U=A2_(x,y), R=C[P,Q], Jac(P,Q)=c in C*, and let B be the integral
closure of R in C(x,y). Use the ACTUAL finite pi:N=Spec B→A2, with the
actual open immersion U⊂N and N\U=D≅A1. The finite field degree is d.
At the generic point of D, the source volume has divisor

    div_N(dx wedge dy)=(e-1)D,
    K_N=(e-1)[D].                                     (1)

There are no other components, since the source form is nonvanishing
on U. The finite normalization, its compactification over P2, resolution
only at infinity and the actual radial end are those proved in
[PF98](PROPER_FINITE_MARKED_END_AND_COMPONENT_OPEN_AUDIT.md).

Take a sufficiently large target ball. Its compact preimage X has
boundary W, and N retracts onto X. The proper radial function has no
critical points beyond this level: off D this follows from étaleness,
and on D from the proper polynomial parametrization of its image.
The proper line meets W in one circle K_D. No original finite-map
hypothesis is retained after deleting a boundary prime.

The relative Thom isomorphism for the smooth divisor D gives

    H_j(N,U;Z)=H_(j-2)(D;Z).

Since U and D are contractible, H_1(N)=H_3(N)=H_4(N)=0 and H_2(N)=Z.
The map H_2(N)→H_2(N,U)=Z is an isomorphism: its coordinate is
intersection with D. Thus c1(O(D)) is the primitive dual generator
of H^2(N)=Z, not a multiple of it. Equivalently the proper class of
D generates H_2(X,W), and its boundary is the ACTUAL K_D.

The intersection form on H_2(X) is negative. To see this without
assuming it for arbitrary Stein fillings, compactify pi as in PF98.
On the rational smooth projective surface S, the pullback H of a
target line at infinity is nef, H^2=d>0, supported outside N. The
image of the generator of H_2(N) is nonzero in H_2(S), because it
intersects the closure of D once. It is orthogonal to H. The Hodge
index theorem makes its square strictly negative. Write it as -h,
where h is a positive integer.

Poincare--Lefschetz duality and the pair sequence now give

    H_1(W)=coker(Z --[-h]--> Z)=Z/h,
    [K_D] generates H_1(W), chi(X)=2, sigma(X)=-1.      (2)

This proves the primitive marking, not merely the order of the group.
The characteristic property of c1(TN)=-(e-1)[D] gives

    h=e-1 modulo2.                                   (3)

Indeed c1 mod2 is w2; on a compact surface representing the generator,
the tangent Euler number is even and its normal Euler number is -h.
Thus the first Chern evaluation, of absolute value e-1, is -h modulo2.
No coordinate ruling or polynomial key hypothesis enters (2)-(3).

## 2. Use the genuine complex structure to obtain a contact invariant

The radial preimage W need not itself be strictly pseudoconvex. The
Levi form of |z|^(2e)+|w|^2 is degenerate at z=0 when e>1. We do NOT
declare that unperturbed boundary a contact manifold.

Choose an algebraic embedding of smooth affine N into affine space,
and let psi be its squared norm, a nonnegative strictly plurisubharmonic
function. On the compact radial core and its boundary collar, choose
epsilon>0 small. The function |pi|^2+epsilon psi is strictly
plurisubharmonic. Its sublevel at the chosen radius is a compact Stein
domain X_epsilon, smoothly isotopic to X by a C1-small collar
perturbation. Its complex tangency contact structure xi is positive
and holomorphically fillable, hence tight. The actual complex structure,
canonical class and all invariants in (1)-(2) are retained. This does
not require a classification of the finite branch singularities.

Gompf's Definition4.2 and Theorem4.5 define the homotopy invariant

    theta(xi)=c1(X_epsilon)^2-2chi(X_epsilon)-3sigma(X_epsilon).

It is defined here because H^2(W) is torsion. With intersection form
[-h] and Chern evaluation -(e-1), this becomes

    theta(xi)=-1-(e-1)^2/h.                           (4)

We import the contact theorem, not a new independent proof of it.
The complete relevant pages were rendered and read; see the
[primary-source manifest](sources/primitive-end-contact/manifest.json).
The formula is Gompf's theta; d3 would be theta/4.

## 3. A complete lens-end restriction, including RP3

Define Y_h to be the oriented boundary of the Euler -h disk bundle
over S2, h>=2, equivalently the usual negatively surgered L(h,1).
An oriented identification with Y_h is required; cyclic H1 alone
does not give such a lens-space identification.

Honda's Theorem2.1 says there are exactly h-1 positive oriented tight
contact structures on Y_h up to isotopy. Gompf's Corollary5.7 proof
constructs Stein contact structures on this SAME oriented circle bundle
with rotation numbers

    r=-h+2,-h+4,...,h-2,
    theta=-1-r^2/h,
    Gamma(xi,s)=(r-h)mu/2 in Z/h.                     (5)

Here s is the fixed spin structure specified by the characteristic
unknot in that surgery diagram. As r runs through the list, Gamma
runs through h-1 distinct classes. Definition4.7/Proposition4.8/
Corollary4.9 make Gamma an isotopy invariant for fixed s. The structures
in (5) therefore exhaust Honda's h-1 tight structures. This argument
also handles even h, when c1 alone would not distinguish them.

Consequently every tight contact structure on the oriented Y_h has
one of the theta values in (5). Comparing with (4) yields

    (e-1)^2=r^2, hence e-1<=h-2 and h>=e+1.             (6)

Honda's 2001 erratum corrects Proposition5.8 and Theorem2.2(1), concerning
nonrotative torus layers. It does not replace the selected lens-space
Theorem2.1. The erratum's full introductory page was checked; no
nonrotative-layer uniqueness is imported here.

At h=2, Y_2=RP3, the only theta value is -1. Formula (4) forces e=1.
For an actual nonautomorphic Keller normalization with just this one
prime, e>=2: otherwise purity makes the full finite map étale, hence
an isomorphism over complex A2. Its remaining plane open could not
omit a divisor, since a defining equation of that divisor would be a
nonconstant unit on A2. Thus

    an actual smooth one-prime Keller normalization has no RP3 end,
    at ANY field degree or ramification index.         (7)

The condition e>=2 is essential for the finite-Keller exclusion. The
actual affine quadric f V=W(W-1), with marked plane W=f x,
V=x(f x-1), has omitted line f=0,W=1 and index e=1. On its other chart
f=u,x=u^-1+t the volume is nonvanishing. This is the A1 affine quadric,
with RP3 end and trivial canonical class; it survives (4)-(6). No
finite Keller target is provided for that affine control.

## 4. Centre-free actual monodromy and the all-degree torus parity gate

Every field automorphism of C(x,y)/C(P,Q) extends to the actual finite
normalization. It preserves the sole ramified divisor D and hence U.
Its restriction is a polynomial deck automorphism of the Keller map.
A nonidentity deck automorphism has no fixed point on U: at a fixed
point the local inverse function theorem for F makes it the identity
on a neighborhood, hence everywhere. The finite deck group would
therefore act freely on A2. Its quotient is an affine variety of finite
CW type; Euler characteristic of that finite topological covering gives

    1=chi(A2)=|Deck|*chi(A2/Deck).

Thus Deck is trivial. If G is the faithful transitive Galois-closure
monodromy group and H its point stabilizer, Z(G) injects into
N_G(H)/H=Aut(C(x,y)/C(P,Q)). A central element in H fixes every sheet
and is the identity. Hence

    Z(G)=1.                                         (8)

This proof uses the actual marked source and sole ramified boundary;
it does not infer source preservation for arbitrary partial covers.

Suppose now C=pi(D) has torus knot T(p,q) at infinity, with positive
coprime p,q>1. The actual full infinity monodromy equals G by PF98's
retained affine-infinity epimorphism. The regular Seifert fiber h_0
of the torus-knot exterior is central, so maps to identity by (8).
Since h_0=mu^(pq)lambda, longitude monodromy is mu^(-pq), not necessarily
identity. It joins no meridian cycles. Every longitudinal degree is
one; in particular D→C has residue degree one at infinity and hence
generically. Each filled branch-preimage core is a regular fiber
upstairs, in the actual fiber framing. The ramified one is K_D.

The lifted Seifert-surface relation of PF101 is integral even with
nonzero shear. After filling it says sum_i e_i[K_i]=0. All [K_i]
are this regular fiber class, and sum_i e_i=d. The primitive marking
(2) therefore gives

    h divides d.                                    (9)

Together (3) and (9) exclude every such torus-at-infinity packet with
odd generic degree d and odd transverse index e. No smoothness or
immersion of the finite image curve is assumed. The source N and D
must still satisfy the stated smooth one-line hypotheses.
The unknot is excluded too: its cyclic exterior group could not have
the actual nontrivial transitive centre-free infinity image G.

## 5. Close every six-sheet triple-ramified torus possibility

Now d=6 and e=3. The preceding torus argument already forces residue
degree one, so a normal branch meridian is one 3-cycle with three
fixed sheets. In fact that residue degree also follows directly:
a generic point of C has at least one retained source preimage, since
a defining equation H_C(P,Q) is a nonconstant source polynomial and
the Keller map contracts no curve. Hence 3f<=5, so f=1.

The global group is A6. It is generated by conjugates of that single
3-cycle, since meridians normally generate the complement group over
simply connected A2. A transitive group generated by single 3-cycles
is alternating on its sheet set: join the supports as a hypergraph;
connected supports give one orbit. A new triple meeting the old set
in two points extends its alternating group by one point. If it meets
in one point, conjugating it by an old triple gives two triples meeting
in two points, and then extends by both new points. This inductive
argument begins with A3 and proves the assertion without a subgroup
census or a finite-degree inference. Thus G=A6, whose centre is trivial.

Lift the Seifert fibration of S3 with branch T(p,q) a regular fiber.
Trivial fiber monodromy means fiber degree one. A permutation cycle
of length l over the exceptional order p gives source order p/l;
similarly for q. The branch normal filling, in the fiber framing,
adds no exceptional order. Finite H1 forces the source base to have
genus zero. Let alpha,beta be the base permutations over the two
exceptional values and tau the branch 3-cycle. Riemann--Hurwitz gives

    c(alpha)+c(beta)+c(tau)=d+2,
    c(tau)=4, hence c(alpha)+c(beta)=4.               (10)

Both alpha,beta are even. Their cycle counts are positive even
integers, so both are two. Their only partitions of six are

    (5,1), (4,2), (3,3), with respective orders 5,4,3. (11)

The orders of alpha,beta divide coprime p,q, so their types differ.
Write p=o_alpha A, q=o_beta B, A,B positive. The source Seifert orders
including order-one entries, fiber order and full homology order are:

| Types | Source cone orders | Order of fiber | Order of H1(W) |
|---|---|---:|---:|
| (5,1),(4,2) | A,5A,B,2B | 3 | 3AB |
| (5,1),(3,3) | A,5A,B,B | 2 | 2AB |
| (4,2),(3,3) | A,2A,B,B | 1 | AB |

These are universal formulas in A,B, not enumerated parameter bounds.
For completeness, for arbitrary cycle lengths a_i,b_j, the lifted
rational Euler number has absolute value d/(pq), by the same
branched-base section-obstruction argument as PF34. The abelianized
Seifert determinant therefore gives

    |H1(W)|=(product source orders)*d/(pq).

Killing the regular fiber leaves the orbifold abelianization, of order
(product source orders)/(lcm source orders). Thus the fiber has order

    d*lcm(source orders)/(pq)
      =d/(gcd(a_i)*gcd(b_j)).                         (12)

Coprimality of p,q separates the two order lists and gives the table.
No source smooth-germ classification at the finite image singularities
is used in this end calculation.

But the ACTUAL K_D is that fiber and GENERATES H1(W) by (2). Each table
row therefore forces AB=1, i.e. A=B=1. All unbounded torus twists are
removed at once. Up to exchange, the only necessary base cases are

    T(5,4), H1=Z/3;
    T(5,3), H1=Z/2;
    T(4,3), H1=0.                                    (13)

For e=3, parity (3) requires h even. This excludes the first and third
cases. In the remaining T(5,3) case the source Seifert orders are
1,5,1,1. A closed orientable Seifert manifold over S2 with at most two
genuine exceptional fibers is a union of two solid tori, hence a lens
space (or S1xS2). Its finite H1 here has order two, so the manifold
is RP3. This is a Seifert-topology identification, NOT an inference
from cyclic H1 alone. It contradicts (7), completing the six-sheet
triple torus exclusion at every torus parameter.

The order and parity arguments also cover a mirror torus knot: use
the positive absolute values of its two parameters. The Euler sign
reverses, while all orders in (9)-(13) stay the same. The last case is
still RP3 with either orientation; quaternion conjugation on S3 descends
to an orientation-reversing diffeomorphism of RP3. It therefore admits
an orientation-preserving identification with Y_2 in either case, and
the positive-tight-contact comparison in (7) still applies. No oriented
Y_h identification for h>2 is inferred from this observation.

## 6. Controls and the exact remaining frontier

The [checker](verify_primitive_end_contact.py) verifies the source pins
and sealed predecessor receipts, exact characteristic/contact identities,
the complete partition argument and the universal two-parameter order
formulas. It also constructs actual connected A6 torus covers for all
three necessary base cases. Their full periphery and exact torsion
core orders, computed by PF101, agree with (13). The T(5,4) cover has
nonzero ramification shear; no longitude-one assumption is substituted.

Twisted actual covers such as T(5,6) have H1 order four but ramified
core order two. They show why homology alone is not the primitive
source marking. The affine index-one RP3 quadric is retained to show
why the finite-Keller and ramification hypotheses cannot be dropped.
Neither control is asserted to be an actual plane Keller counterexample.

The [receipt](primitive_end_contact_receipt.json) supports the written
unbounded arguments; finite controls are not their universal proof.
Honda's contact classification, Gompf's homotopy invariant and Stein
contact constructions are explicit external inputs, not independently
reproved here. No formal proof assistant or external review is claimed.

The next complete targets keep the original scope: exclude genuinely
iterated infinity knots for the remaining smooth triple-ramified packet,
handle singular or multiple-prime normalizations, or prove the global
actual Laurent-key/finite-algebra incompatibility. The all-degree
parity, primitive-mark and canonical-contact constraints are available
for those attempts. Do not replace the actual finite cover by a toy
cover, restart quartic realization or expand coefficient/parameter
campaigns. Closure comes before external-review expense; the broader
Wright/JC2 goal remains ACTIVE.
