# Nine explicit chord fibrations on the new Mestre parent

The `u=11` determinant468 parent now has nine newly compiled elliptic
presentations, with exact birational maps, rational quartic sections and
complete generic-rank calculations. The initial `O+C0` pencil has arithmetic
generic rank **7**. Eight further pencils selected from a frozen finite
divisor dictionary have ranks **8–10**, below the original presentation's
exact rank11.

These are actual equations on the new parent, not hypothetical lattice
frames. Their fibre configurations differ from the original presentation.
We do not claim the nine are pairwise inequivalent. None satisfies this
experiment's higher-generic-rank construction gate, and no point search
follows. This does not prove inferior discovery yield, bound specialized
elliptic ranks, or exclude other fibrations on the parent. A useful rank16
fibration remains unconstructed.

## Why this pencil exists

The [full NS matrix](MESTRE_FULL_NS_GRAMS_2026-09-07.md) gives `O.C0=2`.
Consequently `D=O+C0` has square zero. It is effective and nef: it has
intersection zero with each of its two components and nonnegative
intersection with every other irreducible curve. Either finite I2
nonidentity component meets D once, so D is primitive and its genus-one
pencil has a rational section. Its degree relative to the old fibration
is `D.F=2`.

The pencil and section criterion is the standard K3 divisor construction;
see [Schütt–Shioda, Proposition12.10](https://arxiv.org/pdf/0907.0298).
Their Corollary6.13 supplies the Shioda–Tate rank formula used below.

More generally this construction applies to the retained rational sections
P with `O.P=2`. All selected sections meet both finite I2 nonidentity
components, and the compiled quartics explicitly recover the rational
sections supplied by those fibres.

Write the old model and pole section as

```
y^2 = x^3 + A(T)*x + B(T),
P = (a/q^2, b/q^3),   deg(q)=2,   gcd(a,q)=1.
```

Choose `c = b/a mod q^2`, with `deg(c)<4`, and put

```
m = (y+y(P))/(x-x(P)),
z = (q*m+c)/q^2,
W = (2*x+x(P)-m^2)/q.
```

The new generic fibre, with base z and fibre coordinate T, is

```
W^2 = ((q^2*z-c)^4 - 6*a*(q^2*z-c)^2 - 8*b*(q^2*z-c)
       - 3*a^2 - 4*A*q^4) / q^6.
```

Every coefficient is exactly divisible by `q^6`, and the resulting
right-hand side has degree four in T. The rational inverse is

```
m = q*z-c/q,
x = (q*W-x(P)+m^2)/2,
y = m*(x-x(P))-y(P).
```

Direct substitution verifies the old equation modulo the new quartic
equation, and the displayed forward coordinates recover z and W. These
are birational maps on the open sets where their denominators are nonzero.
At both old finite I2 base values the quartic value is an exact square in
Q[z], giving four displayed rational points on the new generic fibre.
Thus its Jacobian is an elliptic model of the same K3 surface over Q.

## Fixed extension after the first pencil

The initial `O+C0` compilation is a separate completed calculation. Its
low generic rank motivates a bounded extension, without changing any
previous point or score experiment.

Use the eleven proved basis sections
`C0,C1,C2,C3,C4,C5,C6,C7,C8,C9,C12`. The extension freezes all **781** words
with support at most three, nonzero coefficients ±1, and first nonzero
coefficient +1. Their exact heights and component classes give **151**
sections with `O.P=2`.

For each D, let r be the rank of the span of the explicitly known vertical
curves, including O and P. The fibre class lies in that span, so those
curves supply at least `r-1` rational reducible-fibre directions. Rational
Picard rank18 therefore bounds the new generic MW rank by `17-r`.

The complete finite accounting is:

| Upper bound from known vertical curves | Number of eligible words |
|---:|---:|
| 8 | 3 |
| 9 | 14 |
| 10 | 25 |
| 11 | 48 |
| 12 | 51 |
| 13 | 10 |

The protocol admits the first eight rows after descending bound, word and
component ordering, with no replacement after compilation failure. All
eight have preliminary bound13 and pass the quadratic-pole preparation.
Every admission is frozen before the first of these eight compilations.
All eight finish within their30-second,2-GiB individual caps and300-second
total ceiling, with one worker. They consume6.077452267 summed supervised
seconds. This dictionary is not a complete height ball, section group or
fibration enumeration. The other53 rows with preliminary bounds above11
remain uncompiled in this experiment.

## Complete fibre geometry and exact generic ranks

For each quartic the standard binary-quartic invariants give the Jacobian
`y^2=x^3-27*I*x-27*J`. Its coefficients have degrees8 and12, discriminant
degree22, and coprime c4 and discriminant. All fibres are multiplicative.
The complete finite factorization, together with I2 at infinity, accounts
for Euler number24.

| Presentation | Geometric fibre configuration | MW rank over Q(z) | MW rank over Qbar(z) |
|---|---|---:|---:|
| Initial O+C0 | 8 I1 + 6 I2 + I4 | 7 | 8 |
| Pencil0 | 14 I1 + 3 I2 + I4 | 10 | 11 |
| Pencil1 | 14 I1 + 3 I2 + I4 | 10 | 11 |
| Pencil2 | 10 I1 + 5 I2 + I4 | 8 | 9 |
| Pencil3 | 9 I1 + 6 I2 + I3 | 8 | 9 |
| Pencil4 | 10 I1 + 7 I2 | 9 | 10 |
| Pencil5 | 11 I1 + 5 I2 + I3 | 9 | 10 |
| Pencil6 | 13 I1 + 4 I2 + I3 | 10 | 11 |
| Pencil7 | 9 I1 + 6 I2 + I3 | 8 | 9 |

Every reducible fibre is at a rational base value. Each I3 or I4 is split,
verified by the node tangent squareclass. Each I2 is nonsplit, but its
single nonidentity component class is rational. Hence all geometric root
directions are rational. The previously proved rational/geometric Picard
ranks18/19 transport through the birational K3 identification, and
Shioda–Tate gives the two exact ranks in the table. No analytic rank or
bounded point-search miss is used. We have not constructed complete
Mordell–Weil point bases for these new models.

The extra reducible fibres were absent from the initial finite list of
known vertical curves. This explains why the preliminary rank13 bounds
were optimistic: they were valid upper bounds from incomplete component
data, not predicted ranks.

## Reproduction and preserved inputs

The [portable input](../../artifacts/generated-results/elliptic-curves/mestre_chord_pencils_v1.json)
contains all nine equations and maps, the complete151-row eligible roster,
selection protocol, exact fibre results, original NS marking and terminal
supervision records. The producer files are
[`compile_mestre_visible_two_neighbor.sage`](../cas/compile_mestre_visible_two_neighbor.sage),
[`prepare_mestre_chord_pencils.sage`](../cas/prepare_mestre_chord_pencils.sage)
and [`compile_mestre_chord_pencil.sage`](../cas/compile_mestre_chord_pencil.sage).

The [standalone verifier](../cas/verify_mestre_chord_pencils.sage) reconstructs
the complete781-word selection, checks the pole sections by exact group
law, verifies the rational inverse by polynomial reduction, checks all
four quartic sections per presentation, and verifies complete discriminant
factorizations by multiplication, squarefreeness and pairwise coprimality.
It then checks every split-component decision and exact generic rank.
It imports no compiler or repository modules.

The fresh-directory run passes in **0.998026215 seconds**, within120 seconds
and2GiB. The [replay record](../../artifacts/generated-results/elliptic-curves/mestre_chord_pencils_portable_replay_v1.json)
binds the copied inputs and source to the terminal transcript.

```sh
sage -python verify_mestre_chord_pencils.sage --input mestre_chord_pencils_v1.json
```

The constructor, fixed roster, eight compilations and geometry audit total
8.356831680 supervised seconds, excluding the standalone replay and prior
parent proofs. All stages complete. There is no live process or queued
point-search follow-up for this experiment.
