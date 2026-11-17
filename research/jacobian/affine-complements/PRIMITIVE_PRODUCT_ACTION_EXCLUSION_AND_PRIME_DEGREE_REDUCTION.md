# Primitive product actions are excluded; the primitive and prime-degree reductions

PF2RULED120, 16 November 2026. No ACTUAL complex plane Keller map has
primitive monodromy preserving a nontrivial Cartesian product structure
with factors of size at least3. This all-degree exclusion is independent
of the classification of finite simple groups (CFSG). Combining it with
the published CFSG-dependent primitive-group classification reduces
primitive monodromy to alternating subset actions, binary affine actions,
and three binary classical actions. In PRIME generic degree the group
must be the natural alternating or symmetric group. The full finite-
normalization equality and the complete research goal remain OPEN.

## 1. Two strict support bounds for EVERY actual branch meridian

Let pi:N->A2 be the ORIGINAL finite normalization of a complex plane
Keller map, with the WHOLE marked open U=A2 and generic degree d>1.
Retain ALL omitted primes. For every DISTINCT irreducible FULL finite
branch-image component C, let sigma_C be a geometric meridian on the
d generic sheets and w_C its TOTAL moved support. Then

    w_C < 5d/6.                                       (1)

If sigma_C is an involution, there is the stronger bound

    w_C < 2d/3.                                       (2)

These are support bounds, not bounds on the largest cycle length.
As in PF118, w_C=sum_(D_i->C)e_i f_i, including every nontrivial
cycle and residue copy.

For (1), w_C<=d/2 is immediate. If 2w_C>d, PF112--113 make C
contractible, forbid unramified primes over C, and require another
RAMIFIED image. PF114's torus-cusp mass theorem applies to this C,
with precisely its one-Puiseux-pair hypothesis supplied by the
contractible-curve classification. In its notation

    6w_C<=5d-15-2F_C<5d.

The local excess charged to C also satisfies eta_C>=w_C/5.
Neither assertion is extended to an arbitrary singular image branch.

For (2), again only w_C>d/2 needs treatment. An involutory meridian
has only index2 ramified cycles. PF105 gives f_i=1 for these primes,
so there are exactly w_C/2 of them. The high-mass curve C is singular
and contractible, and its cusp costs at least one excess on EVERY
such prime (PF112, or the local PF105/PF119 calculation). PF113 also
requires another ramified image, costing at least2. Therefore

    d-1 >= w_C+w_C/2+2,
    w_C <= 2(d-3)/3 < 2d/3.

Other ramification indices elsewhere are unrestricted. We have not
assumed an entirely simple ramification divisor to prove (2).

## 2. A uniform wreath-product fixed-point dichotomy

Let Gamma have size b>=3 and s>=2, with Omega=Gamma^s and d=b^s.
Use the PRODUCT action of Sym(Gamma) wr Sym(s), whose top permutation
permutes the Cartesian coordinates. If x has nontrivial top permutation,
then exactly one of the following alternatives holds:

    (a) x is an involution and |Fix(x)|/d=1/b;
    (b) |Fix(x)|/d<=1/8.                               (3)

Here is the full calculation. Write x=(g_1,...,g_s;pi). Each cycle of
pi contributes the number of fixed points of the product of the g_i
around that cycle. Thus if pi has c cycles, at most b^c tuples are
fixed. If s-c>=2 the ratio is at most b^(-2)<=1/9.

The remaining case is one transposition and all other coordinates
fixed by pi. If its cycle product and every g_i on a fixed coordinate
are identities, x is conjugate by a base element to a pure coordinate
swap. This gives (a). Otherwise one cycle product is a nonidentity
permutation of b letters and fixes at most b-2 letters. The ratio is
at most (b-2)/b^2<=1/8, since

    b^2-8(b-2)=(b-4)^2>=0.

This proves (3) for all b,s without enumerating wreath products.
It distinguishes product actions of degree b^s from the ordinary
IMPRIMITIVE wreath action of degree bs.

## 3. The complete primitive product-action exclusion

Suppose the actual monodromy G is primitive and preserves the above
Cartesian product structure. Every branch meridian lies in the wreath
product. If its top permutation were nontrivial, alternative (3b)
would give w>=7d/8>5d/6, contradicting (1). Alternative (3a) would
give an involution with w=(1-1/b)d>=2d/3, contradicting (2).
Hence EVERY affine divisorial inertia has trivial top permutation.

Let K be the kernel of G->Sym(s), a NORMAL subgroup of G. All affine
inertia lies in K. The normalization of A2 in the corresponding
quotient of the Galois closure is therefore unramified in codimension
one. Purity makes it finite etale, and connected finite etale covers
of complex A2 are trivial. Thus the top quotient is trivial.
Equivalently, affine meridians normally generate the monodromy;
we have NOT assumed one chosen meridian per image generates G.

With trivial top action, the partition of Gamma^s according to its
first coordinate is G-invariant and nontrivial. This contradicts
primitivity. Thus ALL primitive Cartesian product-action plane Keller
monodromies with b>=3,s>=2 are excluded in every degree.

No intermediate field is assumed to have a polynomial-plane model,
and no generic deck transformation is assumed regular on the source.
This proof uses the actual whole-plane geometry only through (1)--(2)
and the finite-cover purity/normal-generation foundation. It uses NO
CFSG, genus census or extra knot parameters.

## 4. Published primitive classification: the remaining necessary containers

Now IMPORT the CFSG-dependent Theorem4 and Table2 of Burness--Guralnick,
*Fixed point ratios for finite primitive groups and applications*,
arXiv2112.03967v2. The inspected author-hosted PDF has 66 pages;
the COMPLETE printed pages4--5, especially Table2 and Theorem4 on p.5,
were rendered and visually checked. The classification is NOT proved
by our finite checker. The older Guralnick--Magaard classification
motivated this connection; we use the modern explicit table rather
than incomplete summaries or OCR-dependent subset parameter bounds.

PF118 gives mu(G)<=d/2<2d/3, so Theorem4 applies. Its product-type
case(vi) is excluded by Section3. Its equal-bipartition case(ii) has
mu/d=(1+1/(n-1))/2>1/2 and fails PF118; its M22:2 case(iii) has
mu=14,d=22 and also fails. Filtering the classical Table2 leaves
exactly the following necessary families, up to permutation isomorphism:

* Alternating subset actions: G=A_m or S_m acting on k-element
  subsets, with 1<=k<m/2. Retain the additional PF118 support test;
  this container list does NOT say every subset parameter passes.
  The exceptional outer A6/S6 descriptions give the same permutation
  containers up to the isomorphisms in the source theorem.
* Primitive binary affine actions V:H, V=(C2)^a, with irreducible
  H<=GL_a(2) containing a transvection. PF118--119 retain their
  one-image/two-cusp-alternative boundary restrictions.
* Sp_(2a)(2) on cosets of O_(2a)^-(2), a>=3; degree
  2^(a-1)(2^a-1).
* O_(2a)^+(2) on nonsingular points, a>=4; the same degree.
  Include its triality-equivalent degree120 descriptions for a=4.
* O_(2a)^-(2) on singular points, a>=3; degree
  (2^a+1)(2^(a-1)-1). The small U4(2).2 degree27 row is this same
  family under U4(2).2=O6^-(2). The a=3 plus-type group has alternating
  socle A8 and belongs in the alternating container instead.

The last three bullets are the single-factor binary classical
alternatives; their nontrivial product actions are ALREADY excluded.
No sporadic, diagonal, twisted-wreath, semilinear-projective or other
primitive action has been silently added by analogy. These are
NECESSARY containers from the published classification and our filter,
not a list of known Keller realizations or a proof that all imprimitive
actions are impossible.

For reproducibility, the discarded Table2 rows can be checked without
group calculations. Put A=2^(n/2-1) in the binary even-dimensional
rows. The signs of 2mu-d are:

    Sp_n(2)/O_n^eps(2):       eps*A;
    O_n^eps(2), singular:     eps*(A+eps);
    O_n^eps(2), nonsingular: -eps*A.

Thus only the minus symplectic, minus singular-orthogonal and plus
nonsingular-orthogonal rows survive. The L_n(2), L_n(3) and Sp_n(2)
projective-point rows fail PF118. U4(3) has d=112,mu=72 and fails.
For the ternary orthogonal rows, the exact positive expressions for
2mu-d are recorded in primitive_product_action.py. Their valid
almost-simple dimension ranges are retained; the formal n=4 plus-type
ternary equality has nonsimple socle and is not an extra exception.
The small isomorphisms and triality are retained explicitly.

## 5. All prime degrees, and a two-image reduction outside alternating socle

If d=p is prime, transitivity makes the action primitive. Cartesian
products and binary affine degrees>2 are composite. Each retained
binary classical degree is a product of two integers>1. Finally an
interior binomial coefficient C(m,k), 2<=k<=m-2, is composite:
m/gcd(m,k)>1 divides it, by k*C(m,k)=m*C(m-1,k-1), and C(m,k)>m,
so that divisor is proper. Therefore every hypothetical prime-degree
plane Keller map has

    G=A_p or S_p in its NATURAL p-letter action.       (4)

Degrees2 and3 are already impossible by PF118; degree4 is already
retired by published Domrina2000/HC4MYKR0. At p=5 PF118 also removes
A5, leaving S5 as the only necessary group container from this argument.
No existence, full prime-degree inversion, or new all-degree solution
is claimed by (4).

If the primitive group has NONALTERNATING socle, the surviving rows
give

    3d/7 <= mu(G) <= d/2.                             (5)

For the two even binary classical degrees, mu/d=(A-1)/(2A-1) with
A>=4, whose minimum is3/7. For the odd orthogonal degree it is
A/(2A+1)>=4/9. The binary affine minimum is1/2. Consequently there
are at most TWO distinct full branch-image components, by the exact
support ledger sum w_C+t0+Eta=d-1.

In fact EVERY branch meridian in this nonalternating primitive case
moves at most half the sheets. If one had w>d/2, PF113 would require
another ramified image of mass at least mu, while PF114 would charge
eta_C>=w/5. The complete ledger would then give

    d-1 >= (6/5)w+mu > (3/5+3/7)d=36d/35>d,

a contradiction. The stronger statement concerns ALL branch meridians,
whereas PF118 initially guaranteed only SOME meridian. If there are
two images, the exact remaining budget is

    t0+Eta <= d-1-2*ceil(3d/7).

No further ramification-cycle classification or global source recognition
is inferred here. In particular the binary classical families, the
alternating subset families and imprimitive actions remain genuine
research obligations.

## 6. What the abandoned blowup test establishes

Blowing up a target value missing from the open source does not preserve
the original affine-plane target hypotheses. Even the local fold
(x,y)->(u=x^e,v=y) produces on the v-chart x^e=v*s after target blowup,
with an A_(e-1) quotient point where the strict ramification curve meets
a new COMPACT exceptional curve. Those intersecting omitted curves
are not the disjoint affine boundary of the original N. Applying the
old quotient-order or affine-divisor theorem to that replacement would
be invalid. This exploratory route supplies no new contradiction and
is not a premise of Sections1--5.

The new primitive reductions leave actual embedded multiplication,
global trace-moment regularity and finite normalization B=C[x,y] OPEN.
Original N, all omitted primes, PF117's closed local trace gate,
PF118's earlier exclusions, PF119's exact ledger, shared priorities
and published quartic retirement are retained. No outreach, historical
priority, independent CFSG proof, formal verification or full goal
completion is claimed.
