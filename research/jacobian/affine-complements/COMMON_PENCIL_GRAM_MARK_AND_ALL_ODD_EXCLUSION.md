# A common-pencil Gram obstruction closes the all-odd primitive-binomial class

15 November 2026. PF2RULED67 proves an all-degree SOURCE obstruction:
when a smooth affine plane extension has a common source-coordinate
pencil, a sole positive-volume omitted line and any number of
zero-volume omitted lines, the end knot of the positive-volume line
has NONZERO self-linking and hence nonzero homology class.

Applied to the actual quartic, this closes the common-pencil sector
left in PF66. Together PF65-PF67 exclude the ENTIRE primitive ACTUAL
binomial two-stage target class with coprime odd first exponents
p,q>=3, at all windings and all source orders/poles. This is an
exclusion as a branch of the surviving Keller packet. The polynomial
immersions of PF63 still exist. Even first pairs, general lower-weight
keys, longer sequences and the full quartic/higher-degree problem
remain open.

## 1. Explicit hypotheses of the general source theorem

Let N be smooth affine over C, with marked U=A2_(x,f) and

    N minus U = R disjoint-union D_1 ... disjoint-union D_n,

all of these curves isomorphic to A1. Assume:

1. f extends to a regular function on N and is constant, with finite
   value lambda_i, on each omitted line.
2. x extends to a morphism N->P1, has a pole on each omitted line,
   and maps it to infinity. Thus (x,f):N->P1_x times A1_f is a
   birational morphism extending the marked source coordinates.
3. The source volume omega=dx wedge df has divisor a*R on N,
   where a>=1. Its order on every D_i is zero.

The number n, ramification index a+1, source vanishing/pole orders,
Laurent supports and finite target degree have NO bounds here. This
is a theorem under the displayed source hypotheses, not a claim that
every higher-degree normalization has such a pencil or is smooth.

They hold for every common-pencilled actual sole-ramified-prime
packet with the whole charts retained in PF48/PF50. In those charts

    f-lambda=z^r,
    x=z^(e-r)*t+h(z),  h has a nonzero constant leading pole,

so 1/x is regular near the ENTIRE boundary line, including every
finite t. All unramified primes have e=1 and the sole ramified
prime has e=a+1>=2. The common pencil of the surviving quartic is
therefore included, with a=1,n=1. The source values need not coincide
for the general theorem.

Let Q be PF61's whole integral intersection lattice in the dual
omitted-line basis. Write G=-Q and let S be the Schur complement
of its D block, or G_RR if n=0. We prove

    S>1,
    lambda_end(K_R,K_R)=(G^(-1))_RR=1/S in (0,1).             (1)

The value in (1) is understood modulo Z. In particular it is
nonzero, so [K_R] cannot be zero. At a=1 it is also the canonical
boundary class. For a>1 the canonical class is a*[K_R]; its
nonvanishing is NOT asserted merely from (1).

## 2. A boundary blowup model that leaves N unchanged

Start with a smooth projective completion of N. Resolve the rational
map to P1_x times P1_f by point blowups OUTSIDE N, since (x,f) is
already a morphism on N. Stacks Lemma54.4.3, tag0C5H, explicitly
allows all centres outside the domain of definition. The resulting
proper birational morphism

    X -> P1_x times P1_f

is a sequence of point blowups by Stacks Lemma54.17.1, tag0C5R.
It is an isomorphism over the marked A2: there the original inverse
open embedding is defined, and a proper birational morphism to a
normal surface has connected fibres. A point in that open embedding
is isolated in its fibre, so no exceptional part can occur over it.
All blowup centres are therefore in the reduced source boundary.

The full boundary Gamma=X minus U is the reduced total transform
of the two infinity sections of P1 times P1. It is SNC, and its
components form a Z-basis of Pic(X). Put T=X minus N. The strict
closures of R,D_i are P1, and each meets T in exactly ONE point.
They cannot meet one another in N. Since Gamma is SNC, each of
these retained components is consequently a LEAF of its full graph.

A boundary blowup at a smooth point of one component creates a
leaf and increases its parent's degree by one. A blowup at a crossing
creates a degree2 vertex and leaves the parents' degrees unchanged.
Existing degrees never decrease. Hence each retained leaf was
created at a FREE centre (only one boundary component through it),
and no later free blowup was made on that retained component.
This is where retaining entire DISJOINT A1s is essential. An arbitrary
valuation cluster or an open curve with extra punctures need not
have this leaf property.

The pole order of omega on the original x-infinity section is -2.
For a boundary component E write b(E)=1+ord_E(omega). At blowups:

    free centre on E:       b(new)=b(E)+1;
    crossing E,F:           b(new)=b(E)+b(F).                 (2)

These are the ordinary Jacobian factors in the two surface blowup
charts. The retained R has b=a+1>=2 and its creation parent has
b=a>=1. Every retained D_i has b=1 and its creation parent has b=0.

## 3. Two ancestry facts, at every extraction depth

Call a retained valuation a descendant of another if its generic
curvette path includes the other's creation centre and continues
over the created component.

First, NO D_i is a descendant of R. A continuation on R must start
at its unique boundary point: a free centre elsewhere would destroy
the final leaf property of R. That boundary point is the intersection
with its creation parent, possibly separated by later blowups. At
the initial edge both endpoint b-values are positive. By (2), every
new component over that edge has positive b, and every new FREE
component over it has b>=2. A retained D_i would have to be such a
free component with b=1, impossible.

Second, the D_i are PAIRWISE incomparable. At D_i's unique boundary
point the initial b-values are (1,0). Under crossing blowups all new
components over that point have b>=1. A free component created on
any new component has b>=2. The only old b=0 endpoint has, over this
particular point, only intersections with the new exceptional chain;
a blowup there is a crossing blowup, not a free one. Thus no new
free b=1 leaf can be created over this point. A different D_j cannot
be a descendant. This is an induction under the two rules in (2),
not a bound on the number of centres.

If R is a descendant of some D_i, there are at least TWO new centres
on R's path after creation of D_i: the first centre at D_i's boundary
point is a crossing, while the final centre creating R is free.
No other D_j can visit these new centres, by pairwise incomparability.
Nor does D_i's own curvette visit them.

## 4. The entire source lattice is a negative Gram matrix

For each omitted line E, let m_E be its vector of positive integral
generic-curvette multiplicities at the successive point centres on
its path, padded with zeros at all other centres. Its final entry
at the creation centre of E is1. Use the orthogonal total-transform
exceptional basis of the blowup surface, with each exceptional square
-1. The dual class checkE has form

    checkE=q_E*F - sum_j m_(E,j)*E_j^*,                     (3)

where q_E=-ord_E(x)>0 and F is the fibre class of f on P1 times P1.
One verifies the base coefficient by pairing principal source
coordinate divisors: a generic x fibre pairs by q_E and a generic
f fibre by0. The exceptional coefficients are the usual curvette
proximity equations. Equivalently (3) follows by induction through
the two blowup charts, starting with the curvette and total-transform
bases, as in PF50's dual calculation. The base class has F^2=0 and
is orthogonal to every E_j^*. Therefore

    G_(E,E')=-Q_(E,E')=m_E dot m_E'.                         (4)

Thus the WHOLE lattice here, not merely a pairwise necessary matrix,
is a Gram matrix of nonnegative integral path vectors.

The final coordinate of each D_i is private among the D vectors:
its own value is1 and all other D_j values are0 by incomparability.
The D block is consequently positive definite. The final coordinate
of R is zero on all D_i, since none is a descendant of R.
In Euclidean vector notation the Schur complement is exactly

    S=dist(m_R, span_R{m_D_i})^2.                            (5)

We now prove the strict bound, including both ancestry possibilities.

If R descends from a D_i, the at least two new centres of Section3
have positive integral R multiplicities and zero values on EVERY
D vector. Their contribution to (5) is at least 1+1=2. Thus S>=2.

If R descends from no D_i, remove R's final coordinate, of value1,
and call the remaining vector v. It is nonzero: R is a free leaf
whose parent has b>=1, whereas the initial infinity section has
b=-1. At least one earlier exceptional centre is therefore needed.
The vector v is zero at every private D_i coordinate. If it were
in their span, those private coordinates would force every coefficient
to be zero, contradicting v!=0. R's final coordinate is orthogonal
to that whole span. Hence

    S=1+dist(v,span_R{m_D_i})^2>1.                           (6)

This also covers n=0, with the span zero. The strict inequality is
not inferred from matrix samples. Equations(2)-(6) prove it for
arbitrary extraction length, arbitrary Laurent characteristic data
and arbitrary numbers of unramified primes.

## 5. Consequences for the actual end

Equation(4) and the private R coordinate make G positive definite.
By PF61, H1(end)=coker Q, with outward linking -Q^(-1)=G^(-1)
modulo Z and the R knot represented by e_R. Block inversion gives
(G^(-1))_RR=1/S. Since S>1, its self-linking lies strictly between
zero and one. This proves (1), including nonvanishing of the mark.

For two omitted lines this is particularly concrete. With
Q=[[-n_R,c_src],[c_src,-n_D]] and D=det Q,

    D/n_D>1,   0<n_D/D<1.                                  (7)

The zero-mark requirement D|n_D of PF61 is impossible in the common
pencil. No inverse-coordinate coefficient bound is involved.

The disjoint-leaf hypothesis cannot be dropped. Take formal common
charts x=y*t+y^-2 and x=v+y^-2. The cluster Gram matrix is

    G=[[5,4],[4,4]],  S=1,

and its first inverse column is integral. But the second retained
component acquires degree2 in the full boundary: the first is created
by blowing up a finite point on it. The two full charts are not a
separated extension by disjoint entire affine lines. This is a useful
sharp countercontrol to any proof using just volume orders1,0 and
positive definiteness. In particular (7) is not a generic inequality
for arbitrary positive definite binary matrices.

## 6. The entire all-odd primitive-binomial quartic class is closed

Retain the FULL surviving quartic packet of PF46, including R,D,
zero excess and actual monodromy. Suppose its branch has an actual
primitive-binomial two-stage key, with coprime odd p,q>=3, winding
r>=2 and eta0. PF59/PF60/PF63 give c=pq-p-q and k=rc+1.

PF63's complete peripheral/tree/determinant gate leaves only C4 at
r=9 mod36, or C3 with4|r and9|k. PF65 excludes the latter when
3|pq, and proves [K_R]=0 in EVERY remaining representation.
PF50's intrinsic coordinate dichotomy is exhaustive:

* The complementary source is excluded throughout this target
  class by PF66's oriented linking congruence n_R=2 mod3.
* The common source satisfies Sections1-5, with a=1,n=1. Its
  R mark is nonzero, contradicting PF65.

There is no remaining source sector or representation in this class.
Thus NO actual quartic Keller normalization can have such an
all-odd primitive-binomial two-stage branch key. This conclusion
covers ALL windings and ALL source orders/poles; general lower-weight
keys and longer sequences are not silently made binomial/two-stage.

The Zannier polynomial eta0 immersions of PF63 are NOT invalidated.
They provide genuine branch-curve controls, sometimes passing every
coarse infinity test. They cannot extend to the full marked quartic
packet in this all-odd class. No finite self-intersection classification
or new polynomial coefficient construction is needed for this exclusion.

## Evidence and the next broader gate

[verify_common_pencil_mark.py](verify_common_pencil_mark.py) checks
actual P1 times P1 point blowups, free/satellite creation and final
leaf degrees, the entire negative-Gram identity, both ancestry cases,
the strict Schur/mark inequalities, arbitrary-count symbolic matrix
controls, and the nonseparated S=1 countercontrol. The proof above
is unbounded; the finite blowups do not certify all extraction trees.
[common_pencil_mark_receipt.json](common_pencil_mark_receipt.json)
pins the proof, checker, predecessor evidence and
[the primary manifest](sources/common-pencil-mark/manifest.json).
The Stacks resolution/factorization lemmas are imported from primary
HTML; no complete independent reproof of surface factorization is
claimed.

The useful next extension is to compute the ramification mark and
self-linking for OTHER actual infinity classes and compare them
with (1). Even first pairs can have noncyclic companion images and
different marks; PF65's all-odd formulas cannot be reused blindly.
Longer keys and higher degree still require their actual peripheral
data and source hypotheses. A single integral source coordinate
alone was not used as an automorphism theorem: the common-pencil
global morphism, disjoint affine lines, volume orders and end mark
are all essential to the argument given here.

No full quartic/higher-degree Keller exclusion, finite Keller map,
formal verification, independent full topological reproof, external
review, novelty priority or completion of the overall goal is claimed.
