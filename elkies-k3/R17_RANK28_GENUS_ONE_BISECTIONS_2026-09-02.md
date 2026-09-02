# Barcode-targeted genus-one bisections on published R17 (2026-09-02)

<!-- status-consumer: EC-K3-R17-RANK28-GENUS1-BISECTION-PILOT 80fa6e59107cc9e6 -->

<!-- status-consumer: EC-K3-R17-RANK28-GENUS1-SIMULTANEOUS-SPLITTING-H10000 40fb0bc465e3e95c -->

<!-- status-consumer: EC-K3-R17-RANK28-GENUS1-MIXED-TRACE-SPLITTING-H10000 c7aa09836b842b60 -->

## Outcome

The first exact positive-control experiment succeeds for all eleven public
exceptional directions at

```text
t0=-9529/5471.
```

One equation-cheapest norm-eight R17 trace is enough:

```text
tau=-P2-P5.
```

For every `Q1,...,Q11`, the exact line-incidence equation selects a rational
member of the same genus-one bisection pencil.  All eleven normalized branch
polynomials are irreducible squarefree quartics over `QQ`, are coprime to the
degree-24 surface discriminant and to the trace denominator, and have the
displayed rational point above `t0`.  The lifted point specializes literally
to `Q_i`, so its cubic two-descent barcode is exactly

```text
x(P)-theta=x(Q_i)-theta,
```

not merely a matching local signature.

This is an exact equation-level construction and a strong validation of the
proposed template.  The `11/11` target incidence itself must not be
overinterpreted: a norm-eight class defines a genus-one pencil, so every
nondegenerate point of the K3 lies on one member.  The nontrivial checks are
that the selected members and their covers are rational over `QQ`, have the
required irreducible smooth quartic branch divisor, reproduce the Kummer
classes exactly, and give independent height-16 anti-invariant sections.

## Why norm eight is the correct first shell

For a degree-two divisor of arithmetic genus one on a rootless elliptic K3,
Proposition F5 of
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md)
gives the all-section threshold

```text
2*d^2-2*g+2 = 8.
```

The complete `R17/2R17` minimum distribution contains 63,925 classes whose
minimum is exactly eight, plus 43 deeper genus-one classes of minimum twelve.
The pilot completely enumerates and equation-ranks the 230,040 stored
minimum-norm-eight representatives and the resulting 63,925 translation
classes, but executes the equation construction only for the cheapest trace.
The selected published-basis vector is

```text
(0,-1,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0).
```

It costs one group addition and has two finite trace poles.

## Exact pencil formula

Write the height-eight trace as

```text
tau=(Nx/h^2,Ny/h^3),       deg(h)=2.
```

There is a unique polynomial `M0` of degree below four with

```text
M0*Nx+Ny == 0 mod h^2.
```

Every regular slope numerator in this genus-one pencil is

```text
M=M0+lambda*h^2.
```

For the line through `-tau`, put

```text
N=M^4-6*M^2*Nx-8*M*Ny-3*Nx^2-4*A*h^4.
```

The section identity and the congruence imply `h^6 | N`.  Since every term
of `N` has degree at most sixteen,

```text
q_lambda=N/h^6
```

has degree at most four.  The residual intersections satisfy

```text
x^2-(M^2-Nx)/h^2*x
  +((M*Nx+Ny)^2-B*h^6)/(h^4*Nx)=0,
```

whose discriminant is `h^2*q_lambda`.  Thus `s^2=q_lambda(t)` is the
normalized double cover.

For a target `Q=(x_Q,y_Q)` on the fibre at `t0`, exact line incidence gives

```text
lambda_Q =
  (h(t0)*(y_Q+y(tau(t0)))/(x_Q-x(tau(t0)))-M0(t0))/h(t0)^2.
```

The cover witness is then

```text
s_Q=(2*x_Q-(M^2-Nx)/h^2 evaluated at t0)/h(t0).
```

The replay checks `s_Q^2=q_lambda(t0)`, the residual quadratic, the line's
`y` equation, the short Weierstrass equation, and literal equality of the
specialized lifted point with `Q_i`.

## Certified common template

| invariant | value |
|---|---:|
| trace shell | R17 norm 8 |
| selected trace | `-P2-P5` |
| pole pattern | two finite poles, `deg(h)=2` |
| slope template | `M=M0+lambda*h^2` |
| targets attempted | 11 |
| exact successes | 11 |
| irreducible squarefree quartics | 11 |
| branch divisors avoiding singular fibres | 11 |
| exact Kummer barcode matches | 11 |
| anti-invariant height on each cover | 16 |

The degree-two base change has `chi=4`.  Each integral lift is disjoint from
zero and has self-intersection `-4`.  The two conjugates meet transversely at
the four branch points, so their difference has height

```text
2*(4-(-4))=16.
```

It is anti-invariant under the cover involution and is therefore independent
of the invariant R17 subgroup.  Each constructed cover consequently has
generic Mordell--Weil rank at least 18 over its genus-one base.

## First simultaneous-splitting gate

The eleven certified quartics were frozen and searched away from the fitted
fibre by two exact bounded procedures.  The compact projective scan enumerates
every primitive

```text
t=a/b,    |a| <= 10000,    1 <= b <= 10000.
```

It visits `121,589,943` rational parameters.  Sixteen good-prime
quadratic-residue masks reduce these to `7,889` possible simultaneous splits;
`15,861` exact integer-square tests leave exactly one: the original positive
control `t0=-9529/5471`, where all eleven quartics split.  Thus there is no new
simultaneous split in the complete displayed box.

Independently, each pointed quartic has a canonical non-torsion Jacobian point
coming from the opposite ordinate above `t0`.  Exact enumeration of multiples
`2P,...,30P` on all eleven Jacobians finds no cross-split.  These points reach
projective coordinate sizes between 386,035 and 682,732 bits, so this is not a
duplicate small-height search.

The negative result is bounded.  It does not prove that any pairwise fibre
product lacks further rational points, nor does it test other norm-eight trace
pencils.  Since no new exact split survived, specialization and quotient
independence are vacuous at this gate rather than silently assumed.

## Mixed-trace gate

The next experiment replaces the single trace by the seven
equation-cheapest distinct finite-pole norm-eight traces.  Each trace is fitted
through all eleven targets, producing `77` exactly verified irreducible
squarefree quartics.  A second complete scan of the same `121,589,943`
primitive parameters keeps a candidate only while square conditions from at
least two distinct trace pencils survive.

The mixed sieve leaves `5,179` candidates and `10,461` exact square tests again
leave only `t0`.  All `77` quartics split there, providing a strong positive
control, while no parameter away from `t0` splits quartics belonging to two
different selected traces.  This rules out a small mixed-trace collision in
the declared box; it does not cover the remaining 63,918 norm-eight trace
classes.

## Interpretation and next gate

This supplies the requested first mechanism template, but not a rank-32
specialization.  The eleven values of `lambda` were solved backward from the
already known `Q_i`; that step explains the positive control and does not
discover an unknown point.  The first non-tautological simultaneous-splitting
experiment and the first seven-trace expansion have now been run and are
negative in their declared regions.  A useful next expansion should change the
arithmetic source of rational points instead of merely enlarging the compact
box: certify additional independent generators on the pointed quartic
Jacobians, move to the norm-twelve reciprocal chart, or begin the proposed
barcode-targeted trisection fallback.  Any new simultaneous hit must still
pass exact specialization and quotient independence before promotion.

Alongside that arithmetic expansion, the equation experiment should be
extended in two controlled directions:

1. run a short prefix of distinct norm-eight trace classes and group them by
   finite/infinite pole pattern and coefficient complexity;
2. include the 43 norm-twelve deep genus-one classes using the reciprocal
   chart where necessary.

A rank-32 claim would still require at least fifteen independent directions
beyond the generic rank seventeen at one rational specialization, exact point
verification, and a specialization independence certificate.  Simultaneous
square conditions or a rational point on a fibre product are search inputs,
not rank proofs.

The geometric motivation is consistent with Garbagnati--Salgado,
[*Rank jumps and Multisections of elliptic fibrations on K3
surfaces*](https://arxiv.org/abs/2505.15159).  The present certificate is a
specific residual-chord calculation on the published R17 equation rather
than an application that strengthens their general theorem.

## Reproduction

Generate the deterministic certificate:

```bash
sage -python \
  elkies-k3/scripts/search_elkies_2026_rank28_genus_one_bisections.sage
```

Replay it byte for byte:

```bash
sage -python \
  elkies-k3/scripts/search_elkies_2026_rank28_genus_one_bisections.sage \
  --check
```

The generated result is
[`../artifacts/generated-results/elkies-k3-r17-rank28-genus-one-bisection-pilot-v1.json`](../artifacts/generated-results/elkies-k3-r17-rank28-genus-one-bisection-pilot-v1.json).

The certificate records every `lambda`, quartic, residual quadratic, lifted
section, rational cover witness, short-model target, Kummer generator and
input hash.  Its proof boundary is one equation-cheapest finite-pole trace;
it does not classify all genus-one pencils or promote a new rank record.

Run and replay the frozen-quartic simultaneous-splitting search with

```bash
sage -python \
  elkies-k3/scripts/search_elkies_2026_rank28_simultaneous_splitting.sage

sage -python \
  elkies-k3/scripts/search_elkies_2026_rank28_simultaneous_splitting.sage \
  --check
```

Its generated result is
[`../artifacts/generated-results/elkies-k3-r17-rank28-simultaneous-splitting-h10000-v1.json`](../artifacts/generated-results/elkies-k3-r17-rank28-simultaneous-splitting-h10000-v1.json).

Run and replay the seven-trace mixed search with

```bash
sage -python \
  elkies-k3/scripts/search_elkies_2026_rank28_mixed_trace_splitting.sage

sage -python \
  elkies-k3/scripts/search_elkies_2026_rank28_mixed_trace_splitting.sage \
  --check
```

Its generated result is
[`../artifacts/generated-results/elkies-k3-r17-rank28-mixed-trace-splitting-h10000-v1.json`](../artifacts/generated-results/elkies-k3-r17-rank28-mixed-trace-splitting-h10000-v1.json).
