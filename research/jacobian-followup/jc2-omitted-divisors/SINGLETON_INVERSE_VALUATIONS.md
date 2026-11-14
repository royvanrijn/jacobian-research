# Singleton pole packets force inverse coefficient and index orders

14 November 2026. PF2OMIT10 applies to the complete residual r7 singleton
packets of PF2OMIT7 and to the marked quartic cusp of PF2OMIT9. These
are necessary local constraints on the unknown global inverse, not a
construction or exclusion of any residual global map.

## 1. The reciprocal coordinate is a unit times a boundary power

Choose an admissible primitive source coordinate T=x+t*y by PF2OMIT3,
and let R_d=sum a_i T^i be its primitive inverse. The coefficient ideal
is one. The normalized projective root cover is the entire finite N,
and T:N->P1 is a morphism whose infinity divisor is exactly sum p_D D,
with positive integer pole orders p_D. This source-finiteness certificate
is essential: an arbitrary primitive rational function need not extend
to P1 at an exceptional boundary point.

At an omitted point q in the present packets, N is smooth and there is
only one original omitted boundary germ D=(s=0). Thus z=1/T is regular
there and

    z=unit*s^p_D.                                    (V1)

Indeed its Cartier zero divisor is p_D D and has no other local
component. In a regular local ring the remaining factor is a unit.
There is no untracked common numerator zero at q.

Take a generic target line through the singleton value, with parameter
epsilon. Avoid the finitely many tangent directions needed below. The
local models from PF2OMIT7 give the following valuations of the geometric
reciprocal roots over the Puiseux field of epsilon:

| Omitted block | Ambient length k | ord(s) on the normalized line | k reciprocal-root valuations |
|---|---:|---:|---|
| Unramified | 1 | 1 | p |
| Fold | 2 | 1/2 | p/2, p/2 |
| Cyclic index3 | 3 | 1/3 | p/3, p/3, p/3 |
| Ordinary cusp, original index2 | 3 | 2/3 | 2p/3, 2p/3, 2p/3 |

For the cusp use pi=(U,V)=(u,w^3+u*w). Its critical equation is
s=u+3w^2. On a generic line u=a*epsilon, V=b*epsilon with b!=0,
w has order1/3, and s has order2/3. The fold and cyclic models are
(u,w^k); a generic line gives w order1/k. The unramified image is
immersive, so its boundary equation vanishes to order one on a generic
line. Multiplication by the unit in (V1) changes none of these valuations.

The normalized pulled-back target line is smooth at all the relevant
points: the local differential has rank at least one and a generic
line equation has nonzero differential there. Its map to the target
line has index k at a length-k block. Therefore its normalized trace
discriminant has epsilon-order

    nu=sum_blocks (k-1).                             (V2)

This also equals the maximal-ideal order of the original normalized
target discriminant. Trace and its determinant commute with this base
change, and the pulled-back normal algebra is already normal near the
fibre. The retained singleton contributes zero to nu.

## 2. Coefficients and the raw discriminant

Let q_1<=...<=q_(d-1) be the reciprocal-root valuations from all omitted
blocks, repeated with their ambient lengths. Over the Puiseux field,

    R_d=g (T-t0) product_j(1-z_j T),

where g is a unit and t0 is the regular retained source root. Elementary
symmetric functions and positivity of all q_j give

    ord(a_i)>=ceil(sum_(j=1)^(i-1) q_j), 2<=i<=d,
    ord(a_d)=sum_j q_j.                              (V3)

The leading equality has no cancellation: it is a unit times the
product of all reciprocal roots. These are maximal-ideal orders at
the singleton target value, because a generic line detects the order
of each of the finitely many coefficient polynomials.

The homogeneous discriminant product formula has unit factors from
the retained root and from g, leaving

    ord Disc(R_d)=2 sum_(i<j) ord(z_i-z_j)
                 >=Delta(q):=2 sum_(i<j) min(q_i,q_j). (V4)

Equal leading terms can increase the order; no distinct-leading-term
assumption is needed for the lower bound. PF2OMIT4 separates normalized
discriminant from inverse-order index. If J is the index determinant,

    ord J >= ceil((Delta(q)-nu)/2).                   (V5)

For the packets here Delta is an integer. The ceiling also preserves
the discriminant parity when a larger pole order makes the displayed
bound have the wrong parity. No raw discriminant exponent is charged
as an extra omitted sheet.

## 3. Complete residual degree-six table

All p_D>=1. Replacing them by one can only lower the right sides of
(V3)-(V5). Thus every residual r7 singleton satisfies this table:

| Packet | Minimal reciprocal-root valuations | Lower orders of (a2,a3,a4,a5,a6) | nu | Raw discriminant order >= | Index order >= |
|---|---|---|---:|---:|---:|
| A0: five unramified | 1,1,1,1,1 | 1,2,3,4,5 | 0 | 20 | 10 |
| A1: one fold, three unramified | 1/2,1/2,1,1,1 | 1,1,2,3,4 | 1 | 13 | 6 |
| A2: two folds, one unramified | 1/2,1/2,1/2,1/2,1 | 1,1,2,2,3 | 2 | 10 | 4 |
| B: cusp and fold | 1/2,1/2,2/3,2/3,2/3 | 1,1,2,3,3 | 3 | 11 | 4 |
| C: cyclic3 and fold | 1/3,1/3,1/3,1/2,1/2 | 1,1,1,2,2 | 3 | 7 | 2 |

These coefficient entries are lower bounds, including for a6. The
exact leading order is obtained from (V3) with the actual p_D. Points
on the same original omitted prime must use the same pole order, even
when their image branches coincide. In B the fold can be another point
of the cusp's prime or belong to the other ramified prime. In A2 both
folds can be points of one prime or of the two distinct ramified primes.
These identifications do not alter the uniform table.

For example, with cusp pole p and fold pole r in B,

    ord(a6)=2p+r,
    Delta=4p+r+12 min(2p/3,r/2).

With cyclic3 pole p and fold pole r in C,

    ord(a6)=p+r,
    Delta=2p+r+12 min(p/3,r/2).

Both have nu=3 but their reciprocal valuations and index bounds differ.
An equality of the special Artin algebras does not identify these maps.

The normal discriminant orders also agree directly with the image
ledger. In A1 exactly one immersed ramified branch passes through the
value; A2 has two. In B a cusp contributes multiplicity2 and the fold
contributes one, whether or not they belong to one image curve. In C
the smooth cyclic image contributes its different weight2 and the fold
contributes one. Unramified omitted images add no normal different.

## 4. Marked quartic refinement and exact local controls

The PF2OMIT9 quartic singleton has just the cusp block. Its pole is
exactly one, so (V3)-(V5) give

    (ord a2,ord a3,ord a4)>=(1,2,2), ord a4=2,
    ord Disc(R4)>=4, ord J>=1.                       (V6)

This recovers its compulsory index zero and adds a coefficient order
test. It remains conditional on the actual marked model and finite map.

The [checker](verify_singleton_valuations.py) constructs monic reciprocal
root algebras for every packet over C[epsilon]. Their normal comparison
algebras are products of C[w] with epsilon=w^k, and the root coordinate
is a constant times w^(r*p), where r=1 except r=2 for a cusp. Further
controls use the unit(1+w) to separate conjugates at larger poles.
It computes the full inclusion matrix in the normal bases, its
determinant, both trace discriminants and all coefficient orders.
The one retained factor is separate and etale.

The p=1 controls attain every raw-discriminant and index bound in the
table. This proves sharpness for these local one-parameter root-algebra
controls, not for global plane Keller maps, source opens, r7 images or
all coefficient bounds. The general result is the written valuation
argument; a finite list of controls is not its proof. Actual image
equations, pole orders and exceptional target coordinates remain unknown.
