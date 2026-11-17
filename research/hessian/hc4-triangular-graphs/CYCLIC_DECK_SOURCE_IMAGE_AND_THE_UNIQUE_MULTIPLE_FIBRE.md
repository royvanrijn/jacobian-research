# Cyclic deck source images have no point gaps and one multiple fibre

HC4MYMQ0. 16 November 2026. ALL actual complex polynomial plane Keller
degrees and every nontrivial generic rational deck SUBGROUP. Full
HC4/JC2 and the actual trace-dual/marked-energy implication remain OPEN.
No novelty, independent replay, formal verification or review is asserted.

This strengthens PF2RULED136's actual almost-surjective quotient: after
MP0's common-source-pencil reduction, its possible isolated missed points
do NOT occur. The quotient actually covered by the source is affine,
has one irreducible multiple fibre of order g, and reduced fibres elsewhere.
Its universal cover has a FREE source-translation action. These are actual
source/target statements, not just rational Darboux field coordinates.

## 1. Retained setup and the finite group on the common source pencil

Let K:X=A2->Y=A2 be an ACTUAL complex polynomial Keller map, with
R=C[x,y], B=C[K], L=Frac(R), F=Frac(B), full finite normalization N,
and ENTIRE etale completion U=Spec(A). Retain all omitted primes.
HC4MYLA0 gives U smooth affine, simply connected, with

    U minus X = D_1 disjoint-union ... disjoint-union D_b0,
    D_i=A1, chi(U)=1+b0,

and makes every generic deck subgroup act regularly and FREELY on U.
Ramified omitted primes remain outside U in the unchanged N.

HC4MYMP0, combining retained PF46/50, proves that all generic deck
groups are CYCLIC. If b0>0 there are source coordinates (x,y) with

    y in A, delta=partial_x preserving A,
    ML(A)=C[y], LND(A)=C[y]*partial_x.

The action delta preserves the WHOLE original X. Its zero order on
D_i is r_i-1, where y=lambda_i on D_i and
r_i=ord_(D_i)(y-lambda_i)>=1. The whole component charts prove this
at EVERY point of D_i. In particular delta is nonzero there exactly
when r_i=1. We use MP0's proof and its retained classical-source limits.

Fix ANY nontrivial deck subgroup G of order g>=2, put E=L^G, and let
q:U->T=U/G. Necessarily b0>0 by LA0's Euler divisibility. MP0 makes
G act FAITHFULLY on A1_y. Center its fixed point at0. A generator
sigma has

    sigma(y)=zeta*y, zeta of exact order g.                       (1)

On the whole generic A1 fibre, sigma(x)=a(y)*x+h(y), a!=0,
h in C(y). Preservation of the ACTUAL source volume dx wedge dy
forces a(y)*zeta=1. Hence

    sigma(x)=zeta^-1*x+h(y).                                    (2)

The rational h need not be polynomial: a deck map can move X.
Equations(1)--(2) are not a source automorphism or a target Ga action.

## 2. The distinguished fibre contains a free orbit of unit components

All fibres of y:U->A1 have the precise scheme divisor

    y^-1(lambda)=L_lambda+sum_(lambda_i=lambda) r_i*D_i,
    L_lambda={y=lambda} subset X is the ENTIRE source A1.        (3)

The source lines are closed in affine U, all omitted components are
vertical, and no extra components exist. The supports in (3) are
pairwise disjoint A1s. Multiplicities are retained.

An automorphism of U sends whole fibre components to whole components,
preserving multiplicity. A nonidentity element stabilizing a component
in y=0 would induce a finite-order automorphism of A1. Such an
automorphism has a fixed point; the identity restriction also fixes
points. Either contradicts the FREE deck action. Thus every component
orbit in y=0 has g components.

In particular the orbit of L_0 has g DISTINCT multiplicity-one A1s,
one of them L_0 and the other g-1 among the D_i. Therefore

    nontrivial order-g deck subgroup
      ==> at least g-1 added components with r_i=1 at y=0.        (4)

They all map to the SAME curve of the original plane target, since
G fixes B. In particular, if every added unramified component has
r_i>1, the generic deck group is TRIVIAL. This does not force the
original Keller map to be invertible.

## 3. The actual source image is affine and has NO isolated gaps

Put W=union_(gamma in G) gamma(X), an open subset of U. In general
intersections of translated boundary curves might leave isolated missed
points; PF136 deliberately did not discard that possibility. Here (3)
removes it: each gamma(D_i) is itself a WHOLE fibre component. Any two
such components are either equal or disjoint, even when they are source
lines rather than members of the original boundary. Consequently

    U minus W = intersection_gamma (U minus gamma(X))

is a finite union of ENTIRE original D_i, with NO extra point part.
It is a union of free G-orbits of curves. This argument uses MP0's
affine source-coordinate ruling, not an assumed permutation action on
the original boundary basis.

Let V=q(X)=q(W). Since q is finite etale, it is open and q^-1(V)=W.
The complement T minus V is exactly the quotient of those whole curve
orbits. Each orbit gives one embedded A1: its component stabilizer is
trivial by freeness and the finite-order A1 fixed-point argument.
Distinct such quotient curves are disjoint. They are Cartier on smooth
affine T. Their complement V is AFFINE. Similarly W is affine.

Thus PF136's largest almost-surjective open T' is EXACTLY V, and

    X=A2 -> V is SURJECTIVE ETALE of generic Galois degree g,
    W -> V is finite etale with group G,
    W=union_gamma gamma(X), pi1(W)=1, pi1(V)=C_g.                (5)

To check simple connectedness, W contains X and its complement over X
is still a union of disjoint smooth A1 divisors. Moving loops off them
makes pi1(X)->pi1(W) surjective. Hence the cover in (5) is universal.
It is the same actual source map and original target factorization
as in PF136, with no hypothetical chart substitution. Calling X->V
a Galois affine pseudo-covering is now justified by SURJECTIVITY.

We do NOT identify V with the ENTIRE etale locus V_E of the finite
intermediate normalization in E. Nor do we replace the full N or its
trace-dual by these selected opens. Those global distinctions remain.

## 4. Exact ruling and fibre counts on the source-image quotient

The invariant function u=y^g gives a SURJECTIVE A1-fibration

    rho:V->A1_u.

At u=0 the image of X is q(L_0), so the fibre has exactly ONE reduced
component C_0=A1. Its multiplicity is EXACTLY g: upstairs the divisor
of y^g has multiplicity g along every translate of L_0, and W->V is
etale. The map L_0->C_0 is an isomorphism because its stabilizer is
trivial. Thus

    rho^-1(0)=g*C_0.                                            (6)

For u=c!=0 the source fibre of y^g consists of g disjoint lines.
Each maps isomorphically to a fibre component of V. They may map to
the same component. Every component of V is reached, and all have
multiplicity one: their lifts in W are translates of those source
unit components. If s_c is the number of reduced components, then

    1<=s_c<=g;
    there are positive integers k_(c,1),...,k_(c,s_c),
    sum_j k_(c,j)=g,                                             (7)

where X->V has exactly k_(c,j) points over EVERY point of that
component. Over C_0 it has exactly ONE point, not g. Only finitely
many c have s_c>1. Thus V has ONE whole multiple fibre, irreducible
of order g, and all other fibres are REDUCED disjoint unions of A1s.

The scheme fibre y=0 in W consists precisely of the g translates of
L_0, all reduced. At every other value its components are again
translates of source unit components. Therefore W->A1_y has REDUCED
fibres everywhere. The source translation delta is nonzero on X and
on every added component of W, by its r_i=1 chart expression. It is
a FREE Ga action on W. This does not make W an affine plane: it can
have several components in a fibre and nonzero Picard group.

There is also a descended action on V. Using (1)--(2),
sigma*delta*sigma^-1=zeta*delta, so

    Delta=y^(g-1)*partial_x

is G-invariant and locally nilpotent on A. Its action preserves every
gamma(X), hence W, and descends to V. Its invariant ring on V is
EXACTLY C[u]: pulling an invariant function back to X gives a polynomial
in y, and G-invariance leaves C[y^g]. No target B-invariance is inferred.

Writing c'=chi(V), the actual fibre and universal-cover formulas give

    c'=1+sum_(c!=0)(s_c-1),
    chi(W)=g*c', #components(W minus X)=g*c'-1.                  (8)

PF136's deleted-orbit Euler budget is preserved. If all s_c=1 then V
is an affine pseudo-plane in the irreducible-fibre definition. Reducible
nonzero fibres are not excluded; surjectivity alone does not remove them.

## 5. The actual image algebra and a primitive quotient generator

There is no Hartogs ambiguity in this cyclic source image:

    O(V)=R intersect E.                                         (9)

One inclusion is pullback. Conversely a rational E-function regular
on X descends to regularity on V along the SURJECTIVE etale X->V.
This identifies KV1's maximal target-differential algebra in E with
the ACTUAL image, with no isolated points filled in to form its spectrum.

The general coordinate primitive-element lemma reproduced in PF53
Section1 gives L=F(y). Its proof is independent of PF53's later
special ramification row: the geometric generic source line maps
finitely and immersively onto an A1-normalized curve, so its polynomial
normalization factor has degree one. Consequently

    E=F(u), u=y^g, [E:F]=d/g.                                   (10)

Indeed [L:F(u)]<=g by T^g-u, while the g automorphisms fixing F,u
give the reverse inequality. This proves (10) for EVERY deck subgroup
without substituting a plane model for E. The quotient degree d/g is
greater than one by the retained Galois Keller obstruction.

The regular map V->Y is etale and the composite X->V->Y is the SAME
original K. The distinguished C_0 maps by the same polynomial immersion
as the source line L_0. Its normalization map is finite and birational;
the retained one-line injectivity theorem forces a self-collision under
nonautomorphy. Its g lifts in W provide a packet of g unramified
normalization points over each normalization branch of that target curve.
This is actual coupled source/target data, not a constructed target map.

Neither (9) nor (10) makes u integral over B. Poles at points of
V_E minus V, ramified omitted primes of N and the finite trace-dual
lattice remain relevant. Rational traces or a primitive field generator
do not establish regular multiplication in that lattice.

## 6. Source audit and the remaining gate

The new combination reads the complete retained MP0, LA0, KV1 and
PF53 canonical proofs and the concurrently registered PF136 in full.
PF136's curve-orbit quotient and Euler argument are credited; the new
steps use MP0 to eliminate isolated gaps, identify the multiple fibre,
force all covered components to have r_i=1, and obtain the free Ga
universal cover and the precise actual image algebra.

For a scope check, Xie's retained primary *Intersections of valuation
rings in k[x,y]*, complete PDF2--3 and32 were read; PDF3 was rendered
and visually checked. Theorem1.2 really uses the unnormalized dual
intersection matrix and applies to every finite set of divisorial
valuations centered at infinity. This corroborates PF50's input scope,
not an independent reproof of Xie's whole potential-theoretic argument.
Miyanishi's retained report PDF38--39 was also read in full: its precise
pseudo-cover definition and general ruling presentation are respected.
There are NO new external theorem imports, source downloads or caches.

The NEW preregistered checker supplies small exact coordinate, fibre,
orbit, multiplicity, Euler, source-seal and scope controls; written
arguments establish the unbounded conclusions. No predecessor
mathematical checker is executed or assurance flag upgraded.

The remaining nontrivial cyclic case therefore has an ACTUAL surjective
plane pseudo-cover onto the specified affine ruled surface V, followed
by an ACTUAL etale map V->A2. Establishing the latter's compatibility
with the full finite normalization/trace-dual remains the hard step.
The trivial-deck case also remains. Free Ga actions, surjectivity and
the unique multiple fibre do not prove inversion or target Ga descent.
Full HC4/JC2 stays OPEN, the published quartic PLANE exclusion remains
retained, and external review/outreach is DEFERRED.
