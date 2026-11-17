# Quadric middle cycles close every constant-trace intrinsic free projection

Date: 2026-11-15. Claims: **HC4MYIW0/1/2**.
This note closes the noncritical rank-one case left open in
[IU0/1/2](FREE_PROJECTIONS_AND_CONSTANT_BASE_TRACE.md) and
[IV0](PRIMITIVE_HELPERS_AND_LEGENDRE_NORMALITY.md).
It does not solve unrestricted HC4 or JC2.

## Statement

Let K=(K1,K2): A2_z -> A2 be a complex plane Keller map. Suppose that

    j(z,w) = (w1,w2,X(z,w),Y(z,w),Z(z,w))

is an actual closed embedding into A5 and, for a nonzero complex number c,

    X^2 + w1 Y + w2 Z = c^2 + w1 K1(z) + w2 K2(z).          (1)

Then K is a polynomial automorphism, with no degree restriction on any helper.

Together with IU0/1, this proves the following complete constant-trace
free-projection result. Let Q be a nondegenerate quadratic polynomial on A5,
with centre o and critical value lambda. Let b1,b2 be independent affine ambient
linear functions with b(o)=0. If an actual closed embedding j:A4_(z,w)->A5
satisfies

    b composed with j = w,
    Q composed with j = mu + w dot K(z),

where mu is any constant and K is plane Keller, then K is an automorphism.
In particular the associated four-variable gradient is a polynomial
automorphism. The case mu != lambda is no longer an open higher-helper gate.

The new ingredient is an integral middle-dimensional cycle, rather than
helper rectification, a jet truncation, or normality transfer to a different
implicit image.

## External foundations and assurance

We use the following standard external facts explicitly, without claiming
that the symbolic checker proves them.

1. Principal divisors on integral subvarieties are zero in the ambient Chow
   group, by the definition of rational equivalence:
   [Stacks 02RV](https://stacks.math.columbia.edu/tag/02RV).
2. For complex varieties, the algebraic cycle map takes Chow classes to
   integral Borel-Moore homology classes. On a smooth oriented real manifold,
   these classes admit intersection with compact complementary cycles, by
   Poincare duality. See Burt Totaro,
   [PCMI notes 1: Chow groups](https://www.ias.edu/sites/default/files/pcmi1.pdf),
   Section 6, printed pages 9-10.
3. A proper homotopy of proper maps preserves the pushforward of a
   Borel-Moore fundamental class. The required special case also follows by
   pushing forward the locally finite fundamental chain of the source
   manifold times the compact interval: its boundary is the difference of
   the two endpoint chains.

The algebraic divisor calculations, all multiplicities, both proper
homotopies, transversality, and the nonzero integral pairing are supplied
below. No computation of the full Chow group of a quadric is imported.
The retained IU0/1 classification and HW0 Keller denominator argument are
the only earlier structural inputs. Source bytes and prior proofs are pinned
in [the source context](sources/quadric-middle-cycle-context.json).

The [checker](verify_quadric_middle_cycle.py) verifies exact algebraic
identities, local divisor charts, a genuine nonlinear closed-embedding
control, the linear rotation, transversality matrices and scope boundaries.
It is not a formal verification, an independent proof review, or a
replacement for the written all-degree proof.

## IW0: a nonzero common plane class on the affine quadric

Work first on the smooth complex affine fourfold

    V = {x^2 + w1 y1 + w2 y2 = 1}.

Let

    P = {w1=w2=0, x=1},       S = {y1=y2=0, x=1}.

Both are closed complex affine planes. Write [ ] for their integral
Borel-Moore classes in H_4^BM(V,Z).

**The plane class is of infinite order.** Inside V is the compact real sphere

    C = {x real, y=conjugate(w), x^2+|w1|^2+|w2|^2=1}.

This is S4, an oriented compact real four-manifold. Its intersection with P
is the single point x=1,w=y=0. At that point the tangent equation of V is dx=0.
The tangent of P consists of dw=0 and arbitrary complex dy; the tangent of C
consists of dy=conjugate(dw), with dw arbitrary complex. These real
four-dimensional subspaces are complementary in the real eight-dimensional
tangent of V. Their integral intersection number is therefore +1 or -1.
Consequently [P] has infinite order. No mod-two intersection is substituted.

**The two plane classes agree.** Put

    a_i=(w_i+y_i)/2,          b_i=(w_i-y_i)/(2i).

The quadratic equation becomes x^2+a1^2+a2^2+b1^2+b2^2=1.
Fix x,a1,a2 and rotate the pair (b1,b2) through the real angle theta,
0 <= theta <= pi. This is a continuous family of complex linear,
quadratic-preserving automorphisms. At theta=pi it exchanges w and y
simultaneously in both pairs; hence it takes P onto S.
The homotopy on P is proper: the interval is compact and the inverse linear
operators are uniformly bounded. Its endpoint maps are holomorphic, so
their induced orientations are the complex orientations. It follows that

    [S]=[P]=alpha,            alpha has infinite order.       (2)

This explicit rotation avoids any convention about quadric rulings.

**Every polynomial section with the same sign has this class.** Suppose

    sigma(w)=(w,x(w),y(w)),   x(0)=1,
    x(w)^2+w dot y(w)=1.

Define

    sigma_t(w)=(w,x(tw),t y(tw)),    0 <= t <= 1.             (3)

The quadric equation holds exactly, including t=0. At t=0 this is S, and at
t=1 it is the given section. The homotopy is proper because it retains w:
the inverse image of a compact target set has bounded w and lies in a
closed subset of that bounded set times the compact interval. Thus the
complex fundamental class of every such section is alpha.

The same conclusions hold for x^2+w dot y=c^2 with c != 0 and x(0)=c,
by scaling all five coordinates by c. Translations of y do not change them.
Both smoothness and the compact sphere argument require c != 0.

## IW1: two principal divisors count the generic sheets

Return to (1). On w=0, X(z,0)^2=c^2. The integral source ring implies that
X(z,0) is one constant sign. Replace c by -c if necessary, so X(z,0)=c.
The base change of the actual closed embedding along w=0 is a closed
embedding A2_z into the plane {w=0,X=c}. Its image is the whole plane:
a closed irreducible subset of the same dimension in that irreducible
plane cannot be proper. Thus the zero-section image is exactly P.

Let d=[C(z):C(K)]. Choose a complex target value beta=(beta1,beta2) in
an open set where its fibre consists of exactly d distinct points
z_1,...,z_d. Such an open set exists because a dominant generically finite
separable morphism is finite etale over some target open set.
Set f_i=K_i-beta_i. All fibres of f1 are smooth, since the Keller differential
has full rank everywhere. In particular f1 is square-free, and
gcd(f1,f2)=1 because their common zero set is finite.

Consider the smooth affine quadric

    V_beta = {X^2+w1(Y-beta1)+w2(Z-beta2)=c^2}.               (4)

Let E be the irreducible equation of the actual closed hypersurface j(A4).
It exists because the ambient polynomial ring is a UFD and the kernel of
the surjective coordinate-ring map has height one. On the source, the
equation (4) pulls back to

    F = w1 f1 + w2 f2.                                     (5)

This polynomial is irreducible: its coefficients f1,f2 are coprime, and over
C(z) it is a nonconstant linear polynomial in w. Hence the closed
threefold D=V_beta intersect V(E) is integral, with coordinate ring

    C[z,w]/(F).

In particular E does not vanish identically on V_beta. The closed
embedding gives an ambient polynomial L whose restriction to the source
is f1. This algebra lift, and not merely a rational lift, is essential.

On the integral threefold D take the two principal divisors of L and w2.
Their effective cycles, pushed into V_beta, are

    div_D(L)  = sum_{i=1}^d [S_i] + [B],
    div_D(w2) = [P] + [B],                                 (6)

where

    S_i = j({z_i} times A2_w),
    B   = j({f1=0,w2=0}).

Here [B] means the sum over its irreducible components. To see that these
are the complete lists of components, use the exact source ideals

    (F,f1) = (f1,w2 f2),
    (F,w2) = (w2,w1 f1).

Every component has dimension two. Their generic multiplicities are one:

- At S_i, w2 is a unit generically and (f1,f2) are local source coordinates.
  The equation F eliminates f2, and f1 cuts a reduced divisor.
- At B, f2 is a unit generically, while w1 is also a unit generically.
  The two principal divisors each cut the same reduced component.
  Distinct components of f1 remain reduced.
- At P, f1 is a unit generically, so F eliminates w1 and w2 is a local
  parameter of the divisor.

These are computations at generic points; intersections of the components
in smaller dimension do not add terms to their fundamental cycles.
There is no omitted component at infinity: D and both divisors are closed
in the affine V_beta, and (6) lists the complete supports there.

Both divisors in (6) are rationally zero after pushforward. Applying the
integral cycle map and subtracting gives

    sum_i [S_i] = [P] in H_4^BM(V_beta,Z).                  (7)

Translate (Y,Z) by beta. Each S_i is a polynomial section over the retained
w coordinates, and X(z_i,0)=c. By IW0 every [S_i] equals [P]=alpha,
with alpha of infinite order. Equation (7) gives

    (d-1) alpha=0,    hence d=1.                            (8)

This argument uses the actual entire sections; no bound on their degrees
or dependence on w appears.

Finally K is birational Keller. The retained HW0 coprime-pullback argument
removes denominators from its rational inverse. For completeness, if
s=a(K)/b(K) with coprime target polynomials a,b and nonconstant b, a source
prime divisor of b(K) would also divide a(K). A source divisor dominates
a target divisor under a quasi-finite Keller map, whereas the common zero
set of a,b has codimension at least two. This is impossible. Thus b is
constant, and both rational inverse coordinates are polynomial. Their
compositions are identities on dense opens and hence everywhere.

## IW2: all constant traces and constant affine base offsets

IU0 classifies the actual closed zero-section for an independent ambient
linear pair b vanishing at the centre of a nondegenerate Q. Its helper
quadric has rank one, two or three.

- Rank three has no A2 component: at the critical value it is a singular
  cone; away from that value its divisor class group is Z. Both alternatives
  contradict the whole A2 component supplied by the closed zero-section.
- Rank two permits only the critical-value reducible case. IU1 gives a
  Keller inverse from the resulting coordinate line.
- Rank one is Q=X^2+b1Y+b2Z after the retained normalization. At the
  critical value IU1 reads off the source coordinate pair directly.
  Away from it IW1 applies, closing the last case.

Thus every constant trace mu closes. This conclusion uses the actual
closed embedding at zero and globally, not just its differential.

The same result applies to a constant affine frame with source-dependent
offsets. Suppose in the actual cotangent expression

    Q composed with j = D(z)+w dot K(z),
    b composed with j = A w+v(z),       A in GL2(C),

and suppose

    D(z)-v(z)^T A^(-T) K(z) = mu

is a constant. Actual source coordinates (z,b) give the literal identity

    Q composed with j = mu+b dot K'(z),    K'=A^(-T)K.

The transformed plane pair is Keller. The constant-trace theorem makes K'
and hence K invertible. No Hessian transformation under a nonlinear
source change is asserted. In the original cotangent coordinates, inverse
recovery of the gradient is explicit: first recover z from K(z), and then
recover w from eta=grad D(z)+DK(z)^T w.

This closes the whole retained constant-trace free-projection family,
including arbitrary higher compensation terms from IV0. It does not
produce such a free projection or a constant trace for an arbitrary
actual HC4 embedding. Nonconstant traces, unaligned intrinsic projections,
actual endpoint representatives and unrestricted HC4/JC2 remain open.
The normality and source-divisor gates of IV1/2 remain valid tools on
those residual cases.

## Replay

    .venv/bin/python research/hessian/hc4-triangular-graphs/verify_quadric_middle_cycle.py

See [the replay receipt](receipt-quadric-middle-cycle.json) and
[scoped validation](validation-quadric-middle-cycle.json). The proof of the
integral cycle statement is the argument above together with the explicitly
listed foundational external inputs; exact symbolic controls alone do not
establish it.
