# Intermediate source algebras and replica fibre coverage

17 November 2026. Actual complex polynomial plane Keller maps, every
embedded intermediate field, and every explicitly given polynomial
source replica where stability is assumed.

Let
```
B=C[P,Q] subset R=C[x,y], Jac(P,Q) in C*,
F=C(P,Q), L=C(x,y), F subset E subset L.
```
Then the ACTUAL intersection algebra

    C=R intersect E

is finitely generated, normal and smooth, and W=Spec C is affine and
etale over Spec B. The specified source map X=A2->W is etale and its
open image has at most finitely many point gaps. No intermediate field
is assumed to be a polynomial plane.

If, in the GIVEN source coordinates,

    delta=p(y)*partial_x, p in C[y] nonzero,
    delta(E) subset E,

then there is a nonconstant polynomial h(y) such that

    ker(delta|C)=C[h(y)].

The actual source map X->W is then SURJECTIVE. Every reduced fibre of
h:W->A1 is a disjoint union of smooth affine lines. Each whole source
line y=c maps isomorphically onto one such component. Its transverse
multiplicity is exactly the multiplicity of c as a root of h(y)-h(c).

The proof uses the canonical etale model for E and explicit source
functions. It uses no common-pencil or boundary-chart theorem, no LA0
boundary classification, and no triangularization of arbitrary LNDs.
It repairs the stated intermediate-algebra and replica-fibre assertions
without invoking a classification of intermediate surfaces. It does not
prove generic first-velocity integrality or JC2.

## 1. Canonical etale model and the actual source image

Let S_E be the full normalization of B in E, let Lambda_E be its
trace-dual, and put A_E=B[Lambda_E]. KU0 and KU1, Sections1 and3 of
[canonical etale reconstruction](CANONICAL_TRACE_DUAL_ETALE_RECONSTRUCTION.md),
give the precise actual factorization

    X=A2 -> U_E=Spec A_E -> Y=Spec B,

where U_E is a smooth integral affine surface, both maps are etale,
and Frac(A_E)=E. The generic degree of the first map is [L:E].
The input is KU0/1; neither KU2 nor its unverified differential
source-layer input is used.

Write V for the image of X->U_E. Etale morphisms are open, so V is a
dense open subscheme and X->V is surjective etale. It is faithfully flat
and quasi-compact: all schemes here are of finite type and Noetherian.

The complement U_E minus V is closed. Write K_1,...,K_s for all its
one-dimensional irreducible components. Any remaining irreducible
components are finitely many isolated closed points. The K_i are
effective Cartier divisors on the regular surface U_E, even if they
are singular as reduced curves. Therefore

    W=U_E minus union_i K_i                              (1)

is affine. Here is the exact affineness argument. The sum of those
Cartier divisors is effective Cartier; on each affine open trivializing
its invertible ideal the complement is a principal affine open.
Its complement immersion is consequently affine, since affineness of
a morphism is local on its target. As U_E is affine, W is affine.

The surface W is smooth integral and etale over Y. It contains V, and

    W minus V is finite of codimension two.              (2)

Neither map X->V nor X->W has been assumed birational. Removing the
unreached curves in (1) does not replace either full finite
normalization S_E or S_L; W is a specified open model whose ring is
identified next.

## 2. Descent and Hartogs identify the actual intersection ring

We claim

    R intersect E=Gamma(V,O_V)=Gamma(W,O_W).              (3)

For the first equality, a regular function on V pulls back to a
polynomial on X and is a rational function in E, giving one inclusion.
Conversely take f in R intersect E. As a rational function on integral
V, f is regular on some dense open V_0. Its actual pullback to X is
regular everywhere by f in R.

The two regular pullbacks of f to X times_V X agree. Indeed this
fibre product is reduced, being etale over smooth X. Its map to V is
etale and open, so the inverse image of the dense open V_0 is dense:
every nonempty open subset has a nonempty open image in integral V.
The two pullbacks agree on that dense open, and therefore everywhere.
Faithfully flat descent of regular functions along X->V now gives a
unique regular function on V with pullback f. It is the specified
rational function f, as seen on V_0. The precise sheaf input is
[Stacks, Lemma35.13.7, Tag023Q](https://stacks.math.columbia.edu/tag/023Q),
applied to the represented functor of maps to A1; its displayed second
proof and equalizer argument were inspected.

For the second equality, W is normal and the complement (2) has
codimension two. Normal Hartogs extension gives identical rings of
regular functions; equivalently the same height-one local rings are
intersected inside E. Thus (3) holds.

Since W is affine, (3) identifies W with Spec(R intersect E) as an
ACTUAL embedded B-algebra model. This proves finite generation,
normality, smoothness and etaleness of C. The original source map has
image V and only the finite point gaps in (2). These conclusions hold
for EVERY E; no derivation-stability hypothesis was used.

## 3. A given replica makes the invariant curve an affine line

Now assume delta=p(y)partial_x preserves E. It preserves R and hence
C=R intersect E, and its restriction is locally nilpotent. It is nonzero
because the unit Jacobian prevents it from killing both P and Q.

Put D=ker(delta|C), a subring of C[y]. We give the required kernel
argument rather than invoke a surface classification. Choose g in C
and N>=1 with a=delta^N(g)!=0 and delta(a)=0, and set

    s=delta^(N-1)(g)/a in C[a^-1], delta(s)=1.

The finite Taylor-subtraction argument gives

    C[a^-1]=D[a^-1][s].                                  (4)

For an element of that localization, subtract its top delta derivative
times the corresponding power of s, and continue. Invariants in the
localization are exactly D[a^-1], since a is invariant. In particular
D has transcendence degree one and contains a nonconstant polynomial
d0(y).

The kernel D is normal in its fraction field. An element integral
over D is integral over the normal ring C, hence belongs to C; its
invariance then puts it in D. Also C[y] is finite over C[d0(y)], by
the monic polynomial for y furnished by d0. The ring D is a
C[d0]-submodule of that finite module. Noetherianity makes D finite
over C[d0], and therefore finitely generated. The same finite module
generators show that C[y] is finite over D.

For completeness, these facts force

    D=C[h(y)], h in C[y] nonconstant.                    (5)

Let T=Spec D, a smooth affine complex curve. The finite inclusion
D subset C[y] gives a finite surjective map A1_y->T. Take smooth
projective completions. The resulting finite map P1_y->barT has
degree n>=1. Riemann--Hurwitz gives

    -2=n*(2g(barT)-2)+sum_z(e_z-1),

so g(barT)=0 and barT=P1 over C. This uses the standard
characteristic-zero Riemann--Hurwitz formula for finite maps of smooth
complete curves, not an intermediate-surface theorem. The exact formula
and its differential proof were inspected in
[Stacks, Section53.12, Tag0C1B](https://stacks.math.columbia.edu/tag/0C1B),
including Lemma53.12.2. Both primary pages are cached with exact hashes
in the receipt.

Moreover the inverse image of T in P1 is precisely A1: C[y] is the
normalization of D in C(y), so the finite normalization over T is
already the given A1. Thus the nonempty finite set barT minus T has
all its preimages in the single point y=infinity. It consists of
one point. Hence T=P1 minus one point is A1. A generator of its
coordinate ring belongs to the embedded D subset C[y], giving (5).

No smoothness of a boundary of U_E, no root chart, and no Rentschler
normal form is assumed. Only the specified replica is used.

## 4. Whole source lines fill every possible point gap

Let u denote the generator h(y) of D. It defines W->A1_u. For any
a in C, its pullback fibre in the source is

    h(y)=a,

whose reduced support is the finite disjoint union of WHOLE source
lines L_c={y=c}, one for each root c of h(c)=a.

The fibre W_a is nonempty because h(y)-a has a source root. Its
support is pure of dimension one: on the smooth integral surface W,
u-a is a nonzero nonunit and defines an effective Cartier divisor.
Every point of this fibre lies on an irreducible curve component.

The generic point of each such component is in V, since W minus V
is finite. A source preimage lies on one of the lines L_c. The
restriction L_c->W is nonconstant: otherwise P,Q would both be
constant on L_c, contradicting the unit Jacobian.

This line image is CLOSED and the line maps finitely onto it.
Indeed choose one of P,Q with nonconstant restriction b(x) to L_c.
The image of C=O(W) in C[x] contains b(x). The equation b(T)-b(x),
divided by its nonzero constant leading coefficient, makes x integral
over C[b(x)]. Hence C[x] is finite over that polynomial ring and
over the intermediate image of C. This is precisely finiteness of
L_c->W onto its closed scheme-theoretic image.

That closed irreducible curve image contains the chosen generic
point, so it is the ENTIRE fibre component. Every component is
therefore covered by a whole source line. Since every point of
W_a lies on a component, W_a is covered. This applies to each a.
In particular no closed point of the finite set W minus V remains:

    X -> W is SURJECTIVE ETALE.                          (6)

This step uses the actual polynomial source, the affine model (3),
and the specified polynomial invariant h. It does not assume a
preexisting classification of the fibres of W.

## 5. Reduced fibres and their exact source multiplicities

Base change (6) to a fibre and then reduce. Reduction commutes with
etale base change, so

    disjoint-union_(h(c)=a) L_c -> (W_a)_red              (7)

is surjective etale. The source is a disjoint union of smooth A1s;
smoothness descends along a surjective etale cover. Thus (W_a)_red
is a smooth affine curve, whose irreducible components are disjoint.

The map from each L_c to its component is finite by Section4,
etale by (7), and surjective. If its degree is d, Euler
multiplicativity for this finite etale cover gives

    1=chi(A1)=d*chi(component).

The Euler characteristic on the right is an integer, so d=1.
A finite etale morphism of degree one is an isomorphism. Therefore
every component is A1, and each source line maps isomorphically to
its component, including at all its closed points.

Write the scheme fibre as

    W_a=sum_j e_j*C_j.

For a root c mapping to C_j, etaleness of X->W preserves the
transverse divisor order, so

    ord_c(h(y)-a)=e_j.                                   (8)

If n_j is the number of distinct source roots whose lines map onto
C_j, every n_j is positive and the polynomial degree identity gives

    sum_j n_j*e_j=deg h.                                 (9)

The fibre multiplicities and covering counts are thus constrained by
ONE actual polynomial h, not independent fibre passports. Neither
h nor the individual components are presumed target coordinates.

## 6. Scope and exact support

The first conclusion replaces the unverified smooth-intermediate-algebra
input for EVERY embedded E: it follows from KU0/1, source-image descent,
Cartier-complement affineness and Hartogs. It establishes finite point
gaps before adding stability.

The second conclusion fills those gaps for a GIVEN polynomial replica
preserving E and proves its entire reduced-fibre geometry. It does not
classify arbitrary source LNDs. No part requires the common-coordinate
completion input or an omitted-boundary root-chart theorem.

The construction does not prove E=F or E=L, make W a polynomial plane,
make X->W finite, or make source first velocities integral over B.
The full finite normalizations and every original prime remain in scope.
Trace perfection, generic first-velocity landing and HC4/JC2 remain open.

The checker supports the descent equalizer, normal Hartogs codimension
ledger, local slice and invariant degree identities, finite source-line
module reductions, reduction/ramification multiplicities, and finite
etale Euler arithmetic. It includes point-gap and nonstable-field
barriers. Its finite controls do not prove the unrestricted geometric
arguments. Riemann--Hurwitz, faithfully flat descent and the canonical
KU0/1 construction are foundational inputs used at the exact scopes
stated above. Internal replay is not formal verification or external review.

Replay:

    .venv/bin/python research/hessian/hc4-triangular-graphs/verify_intermediate_source_algebra.py
