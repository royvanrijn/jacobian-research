# Weil reciprocity closes six more common-quartic strata on Q80

On the actual direct11952 alternate-Q80 MW17 parent, two rational
anti-invariant sections on a quadratic cover branched only at smooth fibres
**cannot differ at exactly one geometric branch point**. Such a point would
be rational, and its 2-torsion residue violates an inherited-section
character at the good prime131.

All four genus-zero cases `k=1,j=0,...,3` and both genus-one cases
`k=3,j=0,1` are therefore empty. This component and the earlier same-branch
gates leave the following cases; the genus-zero row is now closed by the
separate [coefficient collision proof](Q80_GENUS_ZERO_COLLISION_CLOSURE_2026-09-14.md).
**Twelve genus-one cases remain open** on Q80:

| Base genus | Allocations outside this component | Branch degree outside agreement |
|---:|---|---:|
|0|`k=0,j=0,...,4`|2|
|1|`k=0,j=0,...,4`|4|
|1|`k=1,j=0,...,3`|3|
|1|`k=2,j=0,...,2`|2|

No positive [MW17 plus two-gain cover](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
is constructed. Other parents and disagreement over degree-two, three or
four branch divisors remain outside this exclusion. The new proof does
not use the older norm-eight singular-pencil replay or a Brauer bound.

## 1. A reciprocity condition on the cubic root curve

Let S be the smooth proper curve `theta^3+A(t)*theta+B(t)=0`. For an
inherited section `T=(X_T,Y_T)`, put `alpha_T=X_T-theta` in Q(S).
Every valuation of alpha_T is even. At a good fibre, a zero occurs at a
smooth 2-torsion point; the section identity writes alpha_T as Y_T squared
divided by a unit. At an I1 fibre a section cannot meet the node: the
local total-space equation `xy=u` precludes two vanishing section
coordinates. The other cubic root is simple, so the square/unit argument
still applies. Poles of X_T have even order. At infinity, multiplication
by `v^4` is a square and gives the same argument in the minimal chart.

If alpha has even divisor, beta has odd divisor exactly a reduced
effective divisor Z, and alpha is a unit at Z, then Weil reciprocity gives

```
Norm_Z(alpha) is a square in the ground field.         (1)
```

Indeed the product of the norms of the tame symbols of alpha,beta is1.
Outside Z both valuations are even, so those contributions are squares.
At Z the alpha valuation is zero and the beta valuation is odd, leaving
(1) modulo squares. This handles common support elsewhere. The norm form
over an arbitrary field is [Musicantov--Yom Din, Corollary2.5](https://arxiv.org/pdf/1410.5391).

Let b be a rational smooth base value, e a rational root of its cubic,
and Z the complementary degree-two divisor of the other roots. They need
not be individually rational. If beta has odd divisor exactly Z, then
for every inherited T with finite X_T(b) and nonzero Y_T(b),

```
Norm_Z(alpha_T)=Y_T(b)^2/(X_T(b)-e),
X_T(b)-e is a nonzero square.                          (2)
```

The argument works over Q and Q_131. A closed branch point of degree
greater than one is not a single rational point for this statement.

## 2. Every rational root residue violates (2)

Use the literal [direct11952 equation and sections](R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md).
Its degree24 discriminant is squarefree at131 and has **no zero on
P1(F_131)**, including infinity. The cubic root curve has good reduction.

| Cubic roots in F_131 | Base fibres | Root points |
|---:|---:|---:|
|0|54|0|
|1|62|62|
|3|16|48|
|Total|132|110|

At each of these110 root points `(b_bar,e_bar)`, the certificate supplies
an actual inherited section T with finite minimal-chart coordinates x,y,
nonzero y, and `x-e_bar` a nonzero nonsquare. All seventeen source sections
pass fresh exact characteristic-zero identities; nine are used in the
selected witnesses. At `(b_bar,e_bar)=(2,19)`, section1 gives `(86,77)`
and nonsquare67. At infinity, root124 and section1 give `(108,39)` and
nonsquare115 in weighted coordinates.

A hypothetical rational `(b,e)` reduces to this complete roster. A root
of a monic integral cubic is integral in the appropriate minimal base
chart. The selected T has unit coordinate denominators and nonzero
ordinate there. Thus `X_T(b)-e` is a p-adic unit with nonsquare residue,
contradicting (2). Properness of S and the weighted infinity calculation
include arbitrary denominators and b=infinity.

Hence **no rational function on S has odd divisor exactly the complementary
two roots above one rational smooth fibre**. This also holds over Q_131.
Neither a fully split cubic over Q nor an empty local root curve is assumed.

## 3. One disagreement produces the forbidden divisor

Let `L=Q(t,sqrt(D))` be a quadratic cover branched only at smooth fibres.
For a rational anti-invariant section `P=(x,r*sqrt(D))`, use its twist
Kummer representative

```
alpha_P=D*(x-theta),   Norm(alpha_P)=D^4*r^2.
```

Away from branch places it has even spectral valuations: the same
square/unit and even-pole arguments apply, and the unramified pullback at
an I1 fibre still prevents a section from meeting the node. Infinity is
treated in minimal weighted coordinates. At a branch place:

- If P specializes to nonzero 2-torsion e, the parity is0 at e and1 at
  each of the complementary roots.
- If P specializes to O, all three parities are0. The formal parameter
  `-x/y` has odd order in the ramification coordinate because P is
  anti-invariant. Thus x has odd pole order on the original base, and
  multiplication by D makes every spectral valuation even.

If P had a nonzero specialization at exactly one geometric branch point,
that point would be rational and alpha_P would have the forbidden odd
divisor of section2. Apply this to P1-P2 to prove the opening statement.
Poles are allowed; no formal halving of a Mordell--Weil word is used.

For the polynomial common-quartic sections, all branch specializations
are nonzero 2-torsion. Equality is exactly the common branch factor
`g=gcd(D,x1-x2)`; use the residual q in the full genus-zero chart.
The genus-zero k=1 and squarefree genus-one k=3 cases have exactly one
disagreement. This closes all six cases, including higher contacts,
repeated-quartic genus-zero presentations, infinity and unpointed conics.

One can also directly use `beta=(x1-theta)*(x2-theta)` in these polynomial
cases. Its odd divisor is the pair of selected roots at the disagreement
point, complementary to the third root.

## 4. Evidence and remaining arithmetic

The [input](../artifacts/generated-results/elkies-k3-q80-single-branch-reciprocity-v1/input.json)
binds the literal parent, all17 section records and both implementations.
The [producer](scripts/certify_q80_single_branch_reciprocity.sage) uses
Sage rational functions and square testing. The
[independent checker](scripts/verify_q80_single_branch_reciprocity.py)
uses integer/Fraction arithmetic: it clears denominators in all17 rational
section identities and rechecks the degree24 gcd, all132 projective fibres
and110 nonsquare witnesses. Its
[receipt](../artifacts/generated-results/elkies-k3-q80-single-branch-reciprocity-v1/independent-replay.json)
records0.154 seconds. Producer CPU was0.031 seconds under a40-second/4GiB
cap. No coefficient elimination, lattice reconstruction, section search
or new specialization-rank calculation was run.

```
python3 research/elkies-k3/scripts/verify_q80_single_branch_reciprocity.py
python3 -m unittest discover -s research/tests -p 'test_q80_single_branch_reciprocity.py' -q
```

Six tests accept a valid character and reject a square, incorrect root,
incorrect section point, zero ordinate and singular fibre. The written
reciprocity and specialization arguments are not formally verified.
The old norm-eight singular-pencil assurance gap remains in the earlier
all-k2 conclusion and is not upgraded by this independent gate.

The [small-prime local splitting witnesses](COMMON_QUARTIC_LOCAL_SPLITTING_GATE_2026-09-14.md)
remain valid. The present result uses the additional squareclass relations
forced by the actual section identities. At residual branch degree2,3
or4, reciprocity couples norms across several branch points; the
single-point contradiction does not apply. The genus-one systems and the
positive MW17 endpoint remain open.

The subsequent [degree-two calculation](Q80_DEGREE_TWO_RECIPROCITY_2026-09-14.md)
confines every remaining genus-zero branch pair to sixteen double-point
residues modulo131. It retains collisions and nodal quadratic reductions;
the subsequent [coefficient proof](Q80_GENUS_ZERO_COLLISION_CLOSURE_2026-09-14.md)
closes these polynomial genus-zero collision strata.
