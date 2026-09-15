# The remaining polynomial genus-zero chart on Q80 is empty

On the literal direct11952 alternate-Q80 MW17 parent, **all five genus-zero
`k=0,j=0,...,4` common-quartic allocations are empty over Q_131**. There is
no coefficient-height or denominator bound. The proof includes the literal
cover scalar, higher contacts, infinity, repeated-quartic presentations and
conics without a rational point.

Together with the earlier `k=1` and `k=2` exclusions, this closes the full
polynomial genus-zero common-quartic chart over Q on this parent. The earlier
`k=2` aggregate retains its explicit norm-eight singular-pencil replay gap;
the present `k=0` component neither uses nor repairs that calculation.
**Five genus-one allocations remain open**, after the separate
[k2 orientation](Q80_GENUS_ONE_K2_NODAL_ORIENTATION_2026-09-14.md) and
[k1 cubic reciprocity](Q80_GENUS_ONE_K1_CUBIC_RECIPROCITY_2026-09-14.md) proofs.
Other parents, sections with
rational abscissas outside this polynomial chart, and covers ramified at bad
characteristic-zero fibres remain outside this result. The positive
[MW17 plus two-gain construction](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
remains unconstructed.

## 1. The exact local problem

Write the fixed equation as `E: y^2=f(x)=x^3+A(t)*x+B(t)` and put `p=131`.
The full genus-zero chart is

```
f(x_i)=d*q*s_i^2,       deg(x_i)<=4, deg(s_i)<=5, i=1,2,        (1)
```

where q is squarefree of degree two and d is its literal scalar squareclass.
The two branch fibres are smooth. For `k=0`, the two sections select
distinct nonzero 2-torsion points at each geometric branch value.
The [complete degree-two reciprocity gate](Q80_DEGREE_TWO_RECIPROCITY_2026-09-14.md)
allows q to be made monic and integral in this fixed coordinate, with

```
q mod p=(t-b)^2,
b in {2,9,12,27,32,35,38,46,48,63,74,75,103,105,107,115}.       (2)
```

Its proof includes split, unramified and ramified quadratic branch algebras,
arbitrary base denominators and infinity. All these sixteen fibres are
smooth and fully split. Distinct roots above a branch value have distinct
reductions, since the fibre discriminant is a unit.

The integral Q80 model is a smooth projective K3 over Z_p. Its reduced
elliptic fibration has squarefree degree24 discriminant, all geometric
singular fibres are I1, and all132 rational projective base fibres are
smooth. These are independently replayed properties of the literal model.

## 2. Denominators would force an elliptic component

For either solution in (1), the curve `B: x=x_i(t)` has class `2O+4F` and
arithmetic genus5. Its function field is `Q_p(t,sqrt(d*q))`: s_i is nonzero,
and the ordinate generates the quadratic extension. Thus B is geometrically
integral with geometrically rational normalization, even for an unpointed
conic.

Choose d with valuation epsilon in `{0,1}`. Suppose the Gauss coefficient
valuation of x_i is `-m<0`, and write `x_i=p^(-m)*X`, `s_i=p^(-n)*S`,
with X and S primitive integral polynomials. Gauss valuations in (1) give
`2n=3m+epsilon`. Consequently, for the unit `delta=d/p^epsilon`,

```
X^3+p^(2m)*A*X+p^(3m)*B = delta*q*S^2,
X_bar^3=delta_bar*(t-b)^2*S_bar^2.                            (3)
```

If r is the multiplicity of b in X_bar, (3) gives
`3r=2+2*ord_b(S_bar)`. Hence r is even and at least two.

The primitive section `p^m*x_coordinate-X(t)` of `O(2O+4F)` defines the
flat closure of B. Its special divisor, as an actual effective cycle, is

```
2O + sum_a ord_a(X_bar,h)*F_a,                               (4)
```

where `X_bar,h` is the homogeneous quartic, including its zeros at infinity.
In particular (4) contains the whole smooth elliptic fibre F_b.

This is impossible for the image limit of a geometrically rational curve.
After finite DVR base extension the normalization map extends to a stable
map of genus zero into the projective K3 model. This uses
[Abramovich--Oort, Theorem2.8 and section2.5](https://arxiv.org/pdf/math/9808074):
stable maps into a projective target are proper also in mixed characteristic.
The geometric special domain is a tree of rational components. Proper
pushforward of cycles commutes with specialization; its image cycle is
exactly (4), because the generic normalization map has degree one onto B.
No rational component can dominate a smooth elliptic curve, by Luroth's
theorem, including inseparable maps. Thus (4) cannot contain F_b.

It follows that **all coefficients of both x_i are p-integral**. This
argument compares the actual flat image cycles, not merely their numerical
classes; adding a vertical elliptic component cannot evade it.

## 3. The scalar and the complete reduced Mordell--Weil group

Integrality in (1) forces s_i integral as well: its Gauss valuation n satisfies
`epsilon+2n>=0`. If epsilon is one, reduction gives `f(x_bar_i)=0`, a
nonzero rational 2-torsion section. This contradicts the height formula on
a 24I1 elliptic K3: every nonzero section has
`height(T)=4+2*(T.O)>0`. Thus d is a unit.

The retained and independently replayed Frobenius certificate at131 gives

```
char_H2(T)=(T-p)^19*(T+p)*(T^2+212*T+p^2).                    (5)
```

See the retained [Frobenius proof](../elliptic-curves/rank-jump/THE_LAST_GLOBAL_CLASS_IS_ABSENT.md)
and its [receipt](../artifacts/generated-results/elliptic-curves/rank_jump_residual_quartic_at_131_verification_v1.json).
The new checker validates the retained bindings and factor arithmetic;
it does not rerun the old surface point count. Removing the zero/fibre
hyperbolic plane leaves a p-eigenspace of dimension17 and a minus-p
eigenspace of dimension one. The quadratic factor does not vanish at
minus p, since `2*p^2-212*p=6550`.

If d has nonsquare residue, the two reductions
`(x_bar_i,(t-b)*s_bar_i)` are height-four sections on the constant quadratic
twist `d_bar*Y^2=f(x)`. Twisting negates Frobenius on the parabolic
`H^1(P1,R^1*pi_*)` part; its p-eigenspace therefore has dimension at most
one. Shioda--Tate bounds the arithmetic twist MW rank by one. It has no
torsion, again by the 24I1 height formula. In a rank-at-most-one group two
height-four points are equal up to sign and have the same abscissa. But
`x_bar_1(b)` and `x_bar_2(b)` must be distinct roots. This excludes the
nonsquare unit case.

For a square unit d, absorb its Q_p square root into each s_i and set d=1.
The reductions

```
T_i=(x_bar_i,Y_bar_i),       Y_bar_i=(t-b)*s_bar_i              (6)
```

are polynomial sections, finite also in the weighted infinity chart, and
have height4. To enumerate all of them, the new check recomputes the
finite height Gram of the original17 Q80 basis sections and gets the
original integral matrix of determinant948. The rational-fibre character
matrix has rank17 over F_2, so this subgroup is 2-saturated. Formula (5)
bounds the full arithmetic rank by17. If its index is I, integrality of
the 24I1 height lattice gives `I^2 | 948`. Two-saturation makes I odd,
and948 has no odd square divisor. Thus I=1: this is the **full** reduced
arithmetic Mordell--Weil lattice, not just a selected subgroup.

## 4. Every simple ordinate zero has only the split lift

Fix one height-four point `(X,Y)` over F_p(t) with `Y(b)=0`. Put `u=t-b`
and `S=Y/u`. Linearize `f(X)=q*S^2` at `q=u^2`, with q monic. There are
13 coefficients: five in X, six in S and two in q. Their degree12 identity
has thirteen equations. Its tangent equation is

```
(3*X^2+A)*delta_X - 2*u^2*S*delta_S - S^2*delta_q = 0.        (7)
```

In homogeneous degrees8 and5, `a=3*X^2+A` and S are coprime. A common
zero would mean Y=a=0, impossible at a smooth fibre and impossible at an
I1 node for a section: the local equation `xy=t` prevents a section from
passing through the node. This includes infinity, where u does not vanish.
Thus (7) implies S divides delta_X. Degrees5 and4 force delta_X=0.
If Y has a simple zero at b, `S(b)!=0`, and (7) then implies
`u^2 | delta_q`. The leading coefficient of delta_q is zero because q is
monic, so delta_q=0 and delta_S=0. The coefficient Jacobian is invertible.
Homogeneous degrees here include lower-degree ordinates with a zero at
infinity; they are not silently removed from the chart.

Each reduced T is an integral word in the original basis by section3.
The corresponding characteristic-zero word has the same height4, hence
polynomial coordinates of degrees at most4 and6. Specialization at the
Gauss valuation is a group map on the good generic fibre and sends this
word to T. Since T is finite, the lifted polynomial coefficients are
p-integral. A simple root b of its ordinate lifts uniquely to a p-adic
root beta. Therefore there is an exact solution reducing to (7):

```
x=X_lift,    q=(t-beta)^2,    s=Y_lift/(t-beta).              (8)
```

Invertibility in (7) makes (8) the unique solution in its coefficient
residue ball. Indeed, divide the difference of any two distinct congruent
solutions by its least p-adic valuation and reduce; (7) would have a
nonzero kernel vector. But (8) has a double branch root. Consequently a
squarefree q solution in (1) requires **both** points (6) to have repeated
ordinate zeros at b.

## 5. Complete repeated-contact census

An independent exact rational LDL enumeration of the17-dimensional Gram
finds norm histogram `{0:1,4:2626}` up to norm4, visiting32104 nodes. It uses
exact rational bounds and integer square roots, not a floating radius.
The1313 vectors up to sign are all retained with their polynomial points.

The finite checker verifies every point equation, the degree bounds and
each word identity at seven distinct smooth fibres: if two height-four
sections were unequal, their difference would have height at most16 and
intersection with O at most6, so seven agreements prove equality. This
gives9191 independent finite group-law identities. Sign changes preserve
the zeros, so the1313 rows exhaust the full signed lattice for this test.

Across the sixteen sites in (2), there are443 ordinate-zero incidences.
All439 simple contacts have independently verified invertible13-by13
coefficient Jacobians. The only four repeated contacts are:

| Section index, starting at0 | b | Root X(b) | Ordinate order |
|---:|---:|---:|---:|
|274|2|58|2|
|738|2|58|2|
|496|27|36|2|
|1230|105|24|2|

At any one b, all repeated contacts select the same root. They cannot
provide the two distinct roots required by k=0. Together with sections2--4,
this proves the asserted local emptiness of all five allocations.

## 6. Evidence, replay and boundary

The [frozen finite input](../artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1/input.json)
binds the literal parent, basis, prior collision table, retained Frobenius
evidence and both finite implementations. The bounded
[Sage producer](scripts/certify_q80_genus_zero_collision_closure.sage) used
0.337 CPU seconds under a40-second/4GiB cap. Its
[result](../artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1/result.json)
binds the full lattice and polynomial-point records. The discovery probe
is preserved with its output, and no frozen script was retagged.

The canonical [complete independent checker](scripts/verify_q80_genus_zero_closure.py)
has a separate [frozen replay input](../artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1/complete-replay-input.json).
It supplies the independent norm-four enumeration and invokes the
[finite arithmetic checker](scripts/verify_q80_genus_zero_collision_closure.py).
The [complete receipt](../artifacts/generated-results/elkies-k3-q80-genus-zero-integral-reduction-v1/complete-independent-replay.json)
records PASS in2.219 seconds. The old norm-ten census is retained for
provenance, but its count is not needed as the sole completeness assertion.

```
python3 research/elkies-k3/scripts/verify_q80_genus_zero_closure.py
python3 -m unittest discover -s research/tests -p 'test_q80_genus_zero_collision_closure.py' -q
```

Nine targeted tests cover simple and repeated contacts, a degree drop at
infinity, invalid points and degree bounds, repeated-root label separation,
and exact enumeration with fractional centres and an invalid Gram.
The integrality, stable-map, Frobenius-twist and Hensel arguments are written
proofs, not formally verified. The independently replayed finite component
does not upgrade the old norm-eight singular-pencil assurance in the
earlier k=2 result. It does not construct an infinite rank-at-least19 source.
