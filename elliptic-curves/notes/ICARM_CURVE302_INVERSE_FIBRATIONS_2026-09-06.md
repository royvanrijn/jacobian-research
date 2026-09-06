# Curve 302: inverse fibrations and a constructed nine-direction K3

An explicit K3 family through curve 302 now carries a **certified rank-nine
section sublattice**, with an injective specialization into the public
31-point group. This is a construction from selected known points, not
recovery of the discoverers' parent. No family explaining twelve or more
directions was found in the finite searches below. The full generic MW
group of the new surface, its exact overlap, and the original parent remain
`UNKNOWN`.

The authority is `EC-K3-CURVE302-INVERSE-NINE-DIRECTIONS` in
[`../../MATH_STATUS.json`](../../MATH_STATUS.json).

## What was searched

| finite search | candidates | result |
|---|---:|---|
| Complete minimum-norm-eight table on 103b2 | 63,925 | exact rational-parameter exclusion for every row |
| Complete minimum-norm-eight table on 08f72 | 63,917 | exact rational-parameter exclusion for every row |
| 31 cubic pencils, each with 207 pointed line covers | 6,417 | largest displayed section-span overlap is nine |

The first two runs extend the already completed 11952 atlas. They reuse the
398 residual-chord machinery without selecting a generic MW rank, include
every minimum-vector multiplicity from one through eight, and retain a
projective no-root prime for every row. The verifier reconstructs **all
127,842** comparison polynomials from the source equations and trace words.
It uses homogeneous degree 24, so a dropped leading coefficient retains a
possible point at infinity. These are table classes, not a deduplicated count
of fibrations across charts. No specialization means no defined overlap
score; these rows are not assigned score zero.

The point-directed search uses all 31 points symmetrically as cyclic
eight-point anchors. Each anchor pencil has nine rational basepoints; each
of its other 23 known points is tested on a line through each basepoint.
The discriminant of the residual quadratic gives a degree-two branch
polynomial in the pencil parameter. Covers are grouped by exact
`Q(u)` squareclass, including the constant twist. All 6,417 branch groups
are singletons within their respective pencils. Thus this tested mechanism
does not make several extra public directions generic on one common cover.
Cross-pencil isomorphism classification, conic/higher multisections, and
nontrivial combinations of the 31 points are outside this finite search.

Curve 398's existing exact control replay passes: its two presentations have
the same integral MW16 subgroup, intersection rank 16, and sum rank 16.

## The actual family

The complete rational coefficients, the nine cubic-model sections, the
Weierstrass coefficients, and the fibre isomorphism are in
[`../../artifacts/generated-results/elkies-k3-curve302-nine-direction-k3-v1.json`](../../artifacts/generated-results/elkies-k3-curve302-nine-direction-k3-v1.json).
Coefficients are stored in ascending powers; no interpolation or decimal
reconstruction is involved.

Let `P_i` denote the public points transported to the short model

```
E_short: Y^2 = X^3 + A*X + B,   A=-27*c4, B=-54*c6,
X=36*x+15,  Y=108*(2*y+x+1).
```

Put `F0=Y^2 Z-X^3-A X Z^2-B Z^3`. Eight independent linear conditions at
`P1,...,P8` give a pencil of cubics. The script chooses a deterministic
primitive integral second member `F1`. The nine basepoints are
`P1,...,P8,R`, where `R=-(P1+...+P8)`. All nine membership and transverse
intersection identities are checked over `Q`.

On the affine line

```
L(t) = P1 + t*(P9-P1)                 (affine coordinates, not group addition)
u(t) = -F0(L(t))/F1(L(t)) = N(t)/D(t),
```

cancel the common basepoint factor. Both `N,D` have degree at most two,
their ratio has degree two, and `u(1)=0`. The explicit pointed family is

```
D(t)*F0(X,Y,Z) + N(t)*F1(X,Y,Z) = 0,   zero section P1.
```

Its eight constant sections are `P2,...,P8,R`; its ninth is `L(t)`.
At `t=1`, translation by `-P1` on `E_short`, followed by the inverse short
model map, identifies its group with the public curve 302.

For the stored source Jacobian `y^2=x^3+a(u)x+b(u)`, write
`a(u)=sum a_i u^i`, `b(u)=sum b_i u^i`. The actual Weierstrass family is

```
y^2 = x^3 + A_K3(t)*x + B_K3(t),
A_K3(t) = sum_{i=0}^4 a_i*N(t)^i*D(t)^(4-i),
B_K3(t) = sum_{i=0}^6 b_i*N(t)^i*D(t)^(6-i).
```

Its coefficient degrees are at most `(8,12)`. The discriminant has degree
24 and is squarefree, with no common zero with `A_K3`; infinity is smooth.
Thus its minimal elliptic surface has 24 `I1` fibres and `chi=2`, hence is
K3. The source pencil similarly has 12 `I1` fibres. The rational
Weierstrass fibre at `t=1` is checked Q-isomorphic to 302. Sections are
explicit on the pointed cubic; generic conversion of their coordinates to
the displayed Jacobian Weierstrass gauge is not part of this certificate.

## Exact specialization score and lattice

In public-point coordinates, the nine specialized sections are

```
P2-P1, ..., P8-P1,
-2*P1-P2-...-P8,
P9-P1.
```

Their integer `9 x 31` matrix has rank nine and Smith factors
`1,1,1,1,1,1,1,1,9`. Every group-law identity is replayed. Independently,
32 exact quadratic-character rows over `F2`, using good primes at most 283,
have rank 31, certifying independence of the public inputs. Consequently
specialization is injective on the nine-section subgroup, and its rational
intersection with the known 31-dimensional span has dimension **nine**.
The displayed integral quotient is `Z^22 + Z/9`; no saturation is claimed.

The eight constant sections are pairwise disjoint and disjoint from zero.
The strict transform of `L` meets the zero exceptional divisor once and
none of the other eight; projection formula gives the same intersection
counts for its lifted section. With no reducible-fibre corrections,
Shioda's height formula gives

```
G = [ 2*I_8 + 2*J_8    3*1_8 ]
    [     3*1_8^T         6   ],      det(G)=4608.
```

This is a positive definite MW sublattice and determines an explicit
rank-eleven NS sublattice `U + (-G)`. It is not an isometric embedding into
the single-fibre numerical Néron–Tate lattice: specialization here certifies
a group injection and rank. The full generic arithmetic rank is between
nine and seventeen, and the overlap of the **whole** generic group is only
known to be at least nine. The upper bound is the arithmetic K3 bound over
`Q`; rank 18 or 20 would require leaving this K3-over-`Q` setting, for
example through further base change. This is a geometric constraint, not
a requested target-rank filter. See [Schütt–Shioda](https://arxiv.org/abs/0907.0298)
and [Schütt](https://arxiv.org/abs/0804.1558).

## Replay and limits

```sh
sage -python elkies-k3/scripts/construct_curve302_nine_direction_k3.sage --check
sage -python elkies-k3/scripts/certify_curve302_inverse_extension.sage
sage -python elkies-k3/scripts/search_curve302_cubic_pencil_overlap.sage --check
sage -python elliptic-curves/cas/verify_icarm_curve398_two_parent_collision.sage
```

The complete extension certificate is
[`elkies-k3-curve302-inverse-fibration-extension-full-v1.json`](../../artifacts/generated-results/elkies-k3-curve302-inverse-fibration-extension-full-v1.json).
Its compact prime vectors replay without local campaign checkpoints. The
cover-grouping certificate is
[`elkies-k3-curve302-cubic-pencil-overlap-v1.json`](../../artifacts/generated-results/elkies-k3-curve302-cubic-pencil-overlap-v1.json);
the script reconstructs its per-pencil local checkpoints. The first atlas
pilot, its original script snapshot, and the full raw run remain under
`artifacts/local/elkies-k3/curve302-inverse-*`.

The historical `check_icarm_curve302_rank31_pinned.py` byte replay failed
because its expected deterministic prime list differs from the list produced
by the current helper. That certificate and its expectations were not
changed. The new construction instead embeds a fresh exact character matrix
and independently checks its rank with Sage. This establishes the needed
31-point independence without claiming the historical byte replay passed.

The nine-direction construction is available for any suitably nondegenerate
point configuration of this form. It therefore supplies a reproducible
baseline, not evidence of 302's hidden provenance. An explained sublattice
of rank twelve or higher remains open after these bounded experiments.
