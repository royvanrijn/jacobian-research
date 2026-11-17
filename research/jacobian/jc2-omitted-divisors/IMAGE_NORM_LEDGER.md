# Explicit first-image equations and exceptional matching in both r7 charts

14 November 2026. PF2OMIT11 combines the complete coefficient charts
PF2R7FAM1/2 with the omitted-prime ledger PF2OMIT6/7. It computes the
known ramified image as a parameter family, not an actual second image
or a polynomial source realization. The norm matrices give exact branch
counts on both immersed charts without assuming a Laurent-family
equisingularity classification.

## 1. A six-by-six defining equation for the known image

Use p_lambda(t),q_lambda(t) from
[PF2R7FAM1](../jc2-r7-image-family/PROOF.md). In the polynomial chart L,
lambda=d with 3d^3+2 nonzero. In the Laurent chart K, lambda=x with
x*(3x^3+2)*(15x^3-4) nonzero. These are exactly the immersed domains by
PF2R7FAM2. The functions p,q have fixed monic degrees6,10 and give a
finite birational normalization A1_t of their image C.

Let U,V be normalized target coordinates. On the free rank-six algebra

    E_lambda=C[U,V,t]/(p_lambda(t)-U)

form the matrix M_lambda(U,V) of multiplication by q_lambda(t)-V in
the basis 1,t,...,t^5. The [exact matrix builder](image_norm.py) reduces
monic powers of t; the complete entries, p and q are recorded in
[image_norm_ledger.json](image_norm_ledger.json). In chart K their
coefficients are Laurent polynomials in x, with no denominator vanishing
on its domain. Then

    H_C(U,V;lambda)=det M_lambda(U,V)                 (N1)

is the reduced irreducible image equation for every allowed parameter.
Indeed the monic norm determinant equals Res_t(p-U,q-V). Its support is
the image, and its generic exponent is the normalization residue degree,
which is one by PF2R7FAM1. It is monic of degree6 in V, so there is no
parameter-dependent nonunit factor. The image has projective degree10:
the degree10 homogenizations of p,q,1 give a basepoint-free birational
map from P1. No expansion of the determinant is needed to specify its
equation exactly.

This replaces an unspecified first-image equation by two explicit
one-parameter determinantal equations, subject to the actual parameter
and the target coordinate changes of PF2R7FAM1. It does not identify
those changes for an unconstructed global map.

## 2. Exact branch and lost-length strata at every target value

At a closed target y=(U0,V0), the cokernel of M(y) is

    C[t]/(p(t)-U0,q(t)-V0).

Its dimension is deg gcd(p-U0,q-V0). Since the parametrization is
immersive, at each common root at least one polynomial has a simple
zero. The gcd is therefore squarefree, even when p alone has a repeated
root. Consequently

    c(y)=#(C^nu)_y=6-rank M(y).                       (N2)

The first ramified omitted prime R has (e,f)=(2,1), is omitted, and
has one length-two fold factor for each of these c(y) parameters. Thus

    L_R(y)=2c(y), with c(y) separate copies of C[w]/w^2. (N3)

Rank6 is off C, rank5 is its smooth locus, and rank<=4 is its singular
locus, set-theoretically on the immersed domain. More generally the
rank-(6-c) stratum is given by all (7-c)-minors zero and at least one
(6-c)-minor nonzero, with the usual endpoint conventions. This is an
explicit constructible description for every target in both charts.
It does not distinguish a node from a tangency with two branches, or
claim equality with a scheme-theoretic conductor or Jacobian ideal.

For a residual degree-six Keller candidate, c>=4 is impossible since
R alone would require length at least8. A rank<=2 target anywhere on a
proposed chart therefore excludes that parameter from these candidates.
PF2OMIT11 supplied this necessary gate. The subsequent
[PF2OMIT12 certificate](NO_FOUR_BRANCHES.md) proves that the locus is
empty on K; it removes no parameter. The other matching tasks remain.

The immersion premise cannot be omitted. For (p,q)=(t^2,t^3) at zero,
the same cokernel has dimension2 but there is only one normalization
point. Such a corank counts the nonreduced parametrization fibre, not
its number of branches. The nonimmersed parameters are explicitly
removed before (N2) or (N3) is used.

## 3. Complete matching at each singular value of C

Put b(y)=#(B^nu)_y for the unknown second ramified image and h(y) for
the optional unramified image H. The exact PF2OMIT6 formulas are

    A: retained=6-2c-2b-h,
    B: retained=6-2c-2b-1_(y=y*),
    C: retained=6-2c-3b.                             (N4)

Here y* is the unique ordinary-cusp value of B in row B. All retained
local factors have length one. If H=C then h=c; if H=B then h=b.
The [ledger](image_norm_ledger.json) enumerates every nonnegative
completion of (N4), keeping these global image groupings separate.
Its entries are necessary numerical completions, not curve realizations.

At a two-branch value c=2 (a node or tangency of C), the full list is:

| Row/grouping | Allowed additional branches | Retained length |
|---|---|---:|
| A, H distinct | b=0,h=0; b=0,h=1; b=0,h=2; b=1,h=0 | 2;1;0;0 |
| A, H=C | b=0,h=2 | 0 |
| A, H=B | b=0,h=0 | 2 |
| B | b=0 or1, neither at y* | 2 or0 |
| C | b=0 | 2 |

The original normalization factors at the first three A/distinct entries
are two folds, plus respectively zero, one or two unramified factors.
The last has three fold factors. No extra length is assigned because
the two image branches are tangent: the two omitted points remain
distinct in the original N.

At a three-branch value c=3 the first prime exhausts length6. Every
other omitted image must avoid it and the source fibre is empty, as
already proved in PF2R7COLL2. In particular the H=C grouping is impossible
for any first image having a triple value.

These inequalities also give global matching restrictions:

- In row C, B is disjoint from Sing(C). Every point of B intersect C
  has c=b=1 and exactly one retained point, hence is a singleton value.
  This imposes no transversality condition on the two target branches.
- In row B, its cusp value y* avoids Sing(C). If y* lies on C, then
  c=b=1 and it is a singleton value.
- In row A, if H shares either ramified image, B avoids Sing(C).
  When H=C, every two-branch value of C has empty source fibre, and
  C cannot have any triple value.

No surjectivity theorem is used to discard an empty fibre.

## 4. Apply the complete polynomial-family collision catalogue

PF2R7COLL1 supplies all singular values in L. They are exactly the rank4
and rank3 loci above, with their full counts:

| Immersed L stratum | Rank4 values | Rank3 values |
|---|---:|---:|
| Generic | 12 nodes | 0 |
| d=0 | 3 nodes | 3 ordinary triples |
| U(d)=0 | 9 nodes | 1 ordinary triple |
| M(d)=0 | 9 nodes | 1 ordinary triple |
| T(d)=0 | 10 nodes and1 tacnode | 0 |

The factors U,M,T and the exceptional target equations are imported
from [the certified collision classification](../jc2-r7-collisions/PROOF.md);
these parameter factors are distinct from the target coordinate U.
For example at d=0 the three triple values are

    (U,V)=(a^2+a,0), 81a^3+135a^2+63a+7=0.

Equivalently they are exactly V=0 with
6561U^3+2916U^2+378U+14=0.

There are three distinct such target values, as proved there. At each,
M has rank3, the three omitted factors have length2, and S and D must
avoid the value. The H=C grouping therefore excludes d=0 and all roots
of the degree9 U and degree3 M factors: thirteen immersed parameter
values in L. This excludes only that image-sharing grouping, not the
entire parameter, other groupings or all residual r7 possibilities.

The [checker](verify_image_norm.py) reconstructs both complete matrices
from the canonical p,q, verifies their exact multiplication relations,
checks resultant/corank controls and every matching row, and checks rank3
at all three d=0 triple values over their squarefree cubic algebra.
The universal norm, rank, immersion and conservation arguments remain
written proofs. The K collision taxonomy, actual second-image equations,
pole orders, global complement relations and source realization remain
unresolved. Explicit first-image families and constructible rank strata
do not supply those missing global data or prove singleton-zero properness.
