# Direct bases for the square/cube transfer block

## Status

The all-`k` freeness statement remains a **conjecture**.  No counterexample
was found, but the calculations below do not prove full-base spanning or rule
out discriminant-supported torsion for arbitrary `k`.

What is proved directly, without the refuted conductor-ribbon norm, is:

1. a uniform elimination theorem reducing `U^2=V^3` to one explicit
   coefficient polynomial in `S,A`;
2. monic relative Groebner bases and free modules of ranks `2,4,8,16` for
   `k=1,2,3,4` over the full extracted `S`-base;
3. filtration-compatible bases indexed by subsets for those four cases; and
4. failure of the simplest monic-quadratic induction already at `k=1 -> 2`.

Exact searches find collided length `2^k` through `k=9`, the binomial
maximal-ideal filtration through `k=7`, and no failure on any coordinate
Artin base direction `S=Z^7+tau Z^(7-j)` with `tau^2=0`.  One order-three
test with `tau^3=0` also has the expected length.  These are bounded
regressions only.

Throughout, the ground field is `Q`.  Every displayed bounded presentation
therefore remains valid after extension to an arbitrary characteristic-zero
field.

## 1. Exact direct presentation

Put

\[
 S=Z^k+s_1Z^{k-1}+\cdots+s_k,
 \qquad
 A=a_{k-1}Z^{k-1}+\cdots+a_1Z+a_0,                         \tag{1}
\]

and divide by the monic polynomial `S`:

\[
 A^2=SQ+R,qquad \deg Q\le k-2,quad \deg R<k.              \tag{2}
\]

Define

\[
 \mathcal F_k(S,A)=9Q^2+8SAQ-48S^2R-64AR.                  \tag{3}
\]

### Theorem 1 (uniform elimination of `U`)

Let `V=S^2+A`.  In the formal neighborhood of `A=0`, there is at most one
monic degree-`3k` polynomial `U`, congruent to `S^3`, whose square is `V^3`.
It is

\[
 \boxed{U=S^3+\frac32SA+\frac38Q}.                         \tag{4}
\]

Moreover,

\[
 64(U^2-V^3)=\mathcal F_k(S,A).                            \tag{5}
\]

Consequently the completed monic square/cube factorization ring is the
completion at `(s_1,...,s_k,a_0,...,a_(k-1))` of

\[
 \mathcal T_k=
 \mathbb Q[s_1,\ldots,s_k,a_0,\ldots,a_{k-1}]
 /\bigl(\operatorname{coeff}_Z\mathcal F_k(S,A)\bigr).       \tag{6}
\]

If only completion along the transverse ideal `(a_0,...,a_(k-1))` is
wanted, use the corresponding adic completion of the same ring.  Formula
(6) is an equality of defining functors and rings, not a comparison of
length or tangent data.

### Proof

The high coefficients of `U^2=V^3` determine the coefficients of monic `U`
successively because their linear coefficient is two.  Equivalently, in the
Laurent series ring at `Z=infinity`, `U` must be the polynomial part of

\[
 S^3\left(1+{A\over S^2}\right)^{3/2}.
\]

For binomial degree `d>=3`, the Laurent degree is at most

\[
 (3-2d)k+d(k-1)=(3-d)k-d<0.
\]

Thus only `d=0,1,2` contribute to the polynomial part.  The polynomial part
of `A^2/S` is exactly `Q`, proving (4).  Substitute `A^2=SQ+R` into
`64(U^2-V^3)`; direct expansion cancels the `S^3Q` terms and gives (3).
This proves (5)--(6).  QED

This theorem is the conceptual gain of the direct approach: all `3k`
auxiliary `u_i` disappear before any Groebner calculation.

## 2. Uniform elimination order

Use graded reverse lexicographic order on the transverse variables

\[
 a_0>a_1>\cdots>a_{k-1}                                    \tag{7}
\]

as implemented by the ordered variable list `(a_0,...,a_(k-1))`, with
coefficient domain

\[
 \mathbb Q[s_1,\ldots,s_k].                                \tag{8}
\]

### Proposition 2 (bounded exact calculation)

For `k<=4`, the reduced relative Groebner basis of the coefficient ideal in
(6) is monic: every leading coefficient is `1`, not merely a nonzero element
of `Q(s_1,...,s_k)`.  Therefore polynomial division works over the full base.
Its standard monomials span and are linearly independent over (8).  The
quotient is free of rank `2^k`; base change and completion preserve this
basis, so there is no hidden torsion on the discriminant in these cases.

This last conclusion uses monicity over (8).  Generic leading terms over the
fraction field, or constant length in selected fibers, would not suffice.
The proposition is certified by exact rational polynomial division in
`scripts/search_transfer_basis.py`; it is deliberately quantified only for
the four displayed values of `k`.

## 3. Exact presentations and subset bases for `k<=4`

The coefficients `a_0,...,a_(k-1)` form a minimal transverse generating set:
at `S=Z^k,A=0`, all equations have order at least two, so their classes form
a basis of the `k`-dimensional transverse cotangent space.

The subset labels below index a basis; they do not assert Boolean
multiplication.  They are chosen so that the element labelled by `I` has
maximal-ideal filtration order `|I|` at the collided fiber.

### `k=1`

Write `S=Z+s_1` and `A=E`.  Then

\[
 \mathcal T_1=\mathbb Q[s_1,E]/(E^2),
\]

with

| subset | basis element |
|---|---|
| `empty` | `1` |
| `{1}` | `E` |

### `k=2`

Write

\[
 S=Z^2+s_1Z+s_2,qquad A=XZ+Y.
\]

The exact relative presentation is

\[
 \mathcal T_2=
 {\mathbb Q[s_1,s_2,X,Y]\over
 (X^3,\ 2XY-s_1X^2,\ Y^2-s_2X^2)}.                       \tag{9}
\]

The monic basis is

| subset | basis element |
|---|---|
| `empty` | `1` |
| `{1}` | `X` |
| `{2}` | `Y` |
| `{1,2}` | `X^2` |

At `s_1=s_2=0`, this is

\[
 \mathbb Q[X,Y]/(X^3,XY,Y^2).                              \tag{10}
\]

In particular, `X^2` is retained.  Any construction killing it cannot map
isomorphically to the factorization block.

### `k=3`

Write

\[
 S=Z^3+pZ^2+qZ+r,qquad A=XZ^2+YZ+T.
\]

The seven monic relative relations are

\[
\begin{aligned}
T^3={}&6T^2r^2+6X^2pr^3-12XYr^3,\\
T^2X={}&2T^2pr+2X^2p^2r^2-4XYpr^2,\\
TX^2={}&2T^2q+2X^2pqr-4XYqr,\\
X^2Y={}&4T^2p+4X^2p^2r-8XYpr,\\
X^3={}&6T^2+6X^2pr-12XYr,\\
2TY={}&-X^2pq+X^2r+2XYq,\\
Y^2={}&-2TX-X^2p^2+X^2q+2XYp.
\end{aligned}                                               \tag{11}
\]

A filtration-compatible subset basis is

| subset | basis element |
|---|---|
| `empty` | `1` |
| `{1}`, `{2}`, `{3}` | `X`, `Y`, `T` |
| `{1,2}`, `{1,3}`, `{2,3}` | `X^2`, `XY`, `TX` |
| `{1,2,3}` | `T^2` |

Although `T^2` is an ordinary quadratic monomial, (11) puts it in the third
maximal-ideal power at the collided fiber.  This is the first warning that
subset cardinality is a filtration index, not ordinary monomial degree.

### `k=4`

Write

\[
 S=Z^4+pZ^3+qZ^2+rZ+t,qquad
 A=AZ^3+BZ^2+CZ+D.
\]

An exact compact presentation is (6), with `Q,R` obtained by dividing this
`A^2` by `S`.  Its reduced relative Groebner basis has thirteen monic leading
monomials

\[
\begin{gathered}
A^2CD, D^3, CD^2, BD^2, B^2D, B^3, AD^2, ABD, AB^2,\\
A^2B, A^3, C^2, BC.                                    \tag{12}
\end{gathered}
\]

Together with (3), this specifies every relation exactly without abbreviating
the long base-dependent tails.  Running

```bash
.venv/bin/python scripts/search_transfer_basis.py --print-relations 4
```

prints all thirteen reduced polynomials over `Q[p,q,r,t]`.

The standard basis, arranged by collided filtration order, is

| subset size | subset labels and basis elements |
|---:|---|
| 0 | `empty -> 1` |
| 1 | `1->A`, `2->B`, `3->C`, `4->D` |
| 2 | `12->A^2`, `13->AB`, `14->B^2`, `23->AC`, `24->AD`, `34->BD` |
| 3 | `123->CD`, `124->D^2`, `134->A^2C`, `234->A^2D` |
| 4 | `1234->ACD` |

The ordinary standard-monomial list has eight quadratic and three cubic
monomials.  Relations such as `A^3=12CD` and `A^2B=2D^2` at the collided
fiber move `CD` and `D^2` into filtration order three, yielding the required
binomial distribution `(1,4,6,4,1)`.

## 4. The three proposed proof mechanisms

### A. Direct relative Groebner degeneration

This succeeds completely for `k<=4`.  It also succeeds in an additional
exact `k=5` run: the reduced basis has 25 monic relations and 32 standard
monomials.  That extra bounded case is evidence, not an induction step.

What is missing for all `k` is a formula for the leading ideal of (3) and a
uniform Buchberger reduction proving that every `S`-polynomial reduces with
unit leading coefficient.  The rapid growth

\[
 1,3,7,13,25,42,78,138,256
\]

in collided Groebner-basis sizes for `k=1,...,9` shows why a list-by-list
argument is not a proof strategy.

### B. Induction by a monic quadratic extension

The simplest version is impossible already from `k=1` to `k=2`.  In (10), a
general element of the maximal ideal is

\[
 e=\alpha X+\beta Y+\gamma X^2,qquad e^2=\alpha^2X^2.
\]

Thus every square-zero element lies in `span(Y,X^2)`, the socle, and
multiplication by it has rank at most one.  If (10) were free of rank two
over a copy of `Q[eta]/(eta^2)`, multiplication by `eta` would have rank two.
Therefore no such free quadratic-extension structure exists on the collided
fiber.  Any viable induction would need additional data and cannot be the
proposed monic quadratic tower.

### C. A filtered algebra with `2^k` graded classes

This is consistent with all computations.  At `S=Z^k`, exact maximal-ideal
linear algebra gives

\[
 \dim m^d/m^{d+1}={k\choose d}
\]

through `k=7`.  For `k<=4`, the tables above lift filtration-compatible
classes to a relative monic basis.

The missing theorem is a full-base filtered straightening law.  A Hilbert
function in the collided fiber gives neither spanning over
`Q[[s_1,...,s_k]]` nor independence away from that fiber, and it cannot rule
out discriminant-supported torsion.

## 5. Counterexample search and logical separation

The direct scripts distinguish the following statements.

| Property | Result |
|---|---|
| Generic rank | `2^k`, from direct one-root factorization after a separated-root etale base change |
| Relative spanning | proved by monic division for `k<=4`; additionally checked for `k=5` |
| Relative independence | proved by distinct standard remainders in the same bounded cases |
| Absence of discriminant torsion | follows from freeness only in those bounded relative cases |
| Collided-fiber length | exactly `2^k` for `k<=9` |
| Collided filtration | exactly `(1+t)^k` for `k<=7` |
| Coordinate dual-number base tests at `k=7` | length `256=2*2^7` in all seven directions |
| One order-three base test at `k=7` | length `384=3*2^7` |
| Uniform relative monicity | open for `k>=6` |

No rank jump or embedded Artin torsion was found.  Conversely, none of the
last four bounded tests proves flatness.  In particular, constant collided
length through `k=9` is not promoted to an all-`k` statement.

The search did find a failure of one proposed mechanism—the quadratic
induction—but not a counterexample to rank `2^k` itself.

## 6. Conjecture and next exact obligation

### Conjecture

For every `k>=1`, the ring (6), completed at the maximally collided point,
is a free `Q[[s_1,...,s_k]]`-module of rank `2^k`.  It admits a basis
`B_I`, indexed by subsets `I subset {1,...,k}`, such that the image of `B_I`
in the collided fiber has maximal-ideal order `|I|`.

No Boolean multiplication is asserted.

### Required next step

Derive directly from (3) either:

1. a uniform monic initial ideal with exactly `2^k` standard monomials and a
   complete all-`k` Buchberger reduction; or
2. a filtered rewriting system over `Q[[s_1,...,s_k]]` whose degree-`d`
   normal forms number `binom(k,d)` and whose overlap ambiguities close
   without dividing by a discriminant.

Until one of these is supplied, Failed All-k Transfer Block remains computational evidence and no
dependent Hensel-Product Affine Equalizer or Master Quotient Theorem claim is restored.  Factorization-Slice Obstruction's affine-complement theorems are
independent; only its explicitly suspended transfer connection touches Failed All-k Transfer Block.

## 7. Verification

Run the primary bounded audit:

```bash
.venv/bin/python scripts/search_transfer_basis.py
```

Extend the collided search to `k=9`:

```bash
.venv/bin/python scripts/search_transfer_basis.py \
  --collision-max 9 --filtration-max 7
```

Test every first-order coordinate direction at `k=7`:

```bash
.venv/bin/python scripts/search_transfer_basis.py \
  --collision-max 1 --filtration-max 1 --artin-k 7
```

Test the recorded order-three direction:

```bash
.venv/bin/python scripts/search_transfer_basis.py \
  --collision-max 1 --filtration-max 1 \
  --artin-k 7 --artin-order 3 --artin-direction 7
```

The calculations are exact over `QQ`.  Their bounds are part of their
statements.
