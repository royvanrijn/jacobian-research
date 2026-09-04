# Comparative Kummer/class-group pressure theorem

<!-- status-consumer: EC-K3-R17-KUMMER-CLASSGROUP-PRESSURE-COMPARISON 74b1dae24470b531 -->

Date: 2026-09-04  
Status: exact six-fibre lower-bound theorem; no exact class groups or prospective predictor

## Result

The known-point Kummer calculation now covers the five published-R17 controls

```text
351 (+8), 356 (+12), 376 (+5), 377 (+6), 385 (+12)
```

and the native alternate-Q80 curve 12 (`+12`).  All six exceptional blocks add
zero bad-prime valuation rank modulo the specialized generic `MW17`.  After
quotienting by the generic contribution and allowing for every norm-positive
unit squareclass, the certified residual cubic class-image lower bounds are

```text
+5 -> 3,  +6 -> 5,  +8 -> 6,  +12 -> 10,11,11.
```

Thus the four observed jump strata are strictly separated by this lower bound.
The cross-frame curve-12 value is `11`, matching curve 356 and exceeding the
curve-385 value `10` only by the unavoidable signature/unit correction.

This localizes the class-group wall: the known exceptional points themselves
force a large residual 2-class image before a complete descent begins, and
they do so without creating new bad-prime valuation directions modulo the
generic subgroup.  It does **not** explain why those exceptional rational
points occur.  It is not yet a prospective rank predictor: the exceptional
points, and hence the known jump, are inputs to the invariant.

The complete machine-readable data set is
[`../artifacts/generated-results/elkies-k3-r17-kummer-classgroup-pressure-comparison-v1.json`](../artifacts/generated-results/elkies-k3-r17-kummer-classgroup-pressure-comparison-v1.json).

## The comparative invariant

Let (K=mathbf Q(\zeta)) be the irreducible completed-square 2-division
cubic, let (B\subset K^*/K^{*2}) be the binary span of all certified known
point classes, and let (G\subset B) be the specialized generic `MW17`
subspace.  Write

\[
 v:B\longrightarrow
 \bigoplus_{\mathfrak p\mid 2\Delta_E}\mathbf F_2
\]

for prime-ideal valuation parity.  Everywhere-even classes have the ideal
square-root map

\[
 c:\ker(v)\longrightarrow \operatorname{Cl}(K)[2].
\]

The precise known-point residual object is

\[
 \Pi(B,G)=
 \frac{c(\ker(v|_B))}{c(\ker(v|_G))}.
\]

This is the fail-closed version of “known Kummer half-ideal span modulo the
generic `MW17` contribution.”  The replay proves the general lower bound

\[
 \dim_{\mathbf F_2}\Pi(B,G)\ \geq\
 \dim(B/G)
 -\bigl(\operatorname{rank}v(B)-\operatorname{rank}v(G)\bigr)
 -(r_1+r_2-1).
\]

Indeed,

\[
 \dim\ker(v|_B)-\dim\ker(v|_G)
 =\dim(B/G)-\bigl(\operatorname{rank}v(B)-\operatorname{rank}v(G)\bigr),
\]

and the kernel of (c) is contained in the norm-positive unit squareclasses.
That space has dimension (r_1+r_2-1): the full unit squareclass space has
dimension (r_1+r_2), while `Norm(-1)=-1` in odd degree removes one direction.

No class-group computation is used in this theorem.  The values are lower
bounds, not exact dimensions of (Pi(B,G)).

## Exact six-fibre data

| curve | frame | known jump | signature ((r_1,r_2)) | `rank v(G)` | `rank v(B)` | residual valuation obstruction | unit ambiguity | proved `dim Pi(B,G)` lower bound | proved full `Cl(K)[2]` lower bound |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 351 | published R17 | 8 | `(3,0)` | 5 | 5 | 0 | 2 | **6** | 18 |
| 356 | published R17 | 12 | `(1,1)` | 7 | 7 | 0 | 1 | **11** | 21 |
| 376 | published R17 | 5 | `(3,0)` | 3 | 3 | 0 | 2 | **3** | 17 |
| 377 | published R17 | 6 | `(1,1)` | 7 | 7 | 0 | 1 | **5** | 15 |
| 385 | published R17 | 12 | `(3,0)` | 12 | 12 | 0 | 2 | **10** | 15 |
| 12 | alternate Q80 | 12 | `(1,1)` | 6 | 6 | 0 | 1 | **11** | 22 |

The total full-class-group lower bound does **not** track the jump: for
example curve 376 (`+5`) forces `Cl(K)[2]` rank at least 17, while curve 385
(`+12`) forces only 15.  The separation appears only after removing the
generic `MW17` contribution.  On this panel the residual valuation obstruction
is zero in every row, so the displayed residual bound simplifies to

\[
 \text{known jump}-(r_1+r_2-1).
\]

That identity is why the finite-panel tracking is exact and also why it cannot
yet be advertised as independent prediction.

## What the theorem localizes

The zero residual valuation-rank increment is the main structural observation,
not the subsequent ordering of the six lower bounds.  On this panel, once that
increment vanishes, the lower bound is forced to be

\[
 \operatorname{jump}-(r_1+r_2-1).
\]

Thus the calculation says where much of the already-known exceptional Kummer
information lives: after adjustment by generic point classes it is
everywhere-even, and—up to the unit ambiguity—maps to genuinely global cubic
2-class directions.  The large jumps are not being accounted for by an
expanding collection of bad-prime valuation-parity directions.

The converse is a different problem.  Since the 2-division cubic is
irreducible, `E(Q)[2]=0`.  For a known rank-`r` subgroup `G` whose Kummer image
has dimension `r`, the standard exact sequence gives

\[
 \dim_{\mathbf F_2}\bigl(\operatorname{Sel}_2(E)/\delta G\bigr)
 = (\operatorname{rank}E(\mathbf Q)-r)+\dim_{\mathbf F_2}\Sha(E)[2].
\]

Consequently even a complete large residual Selmer space would not determine
how many directions lift to rational points.  The theorem therefore explains
why a full-BNF front door encounters substantial global 2-class information;
it does not explain the production of Mordell--Weil directions rather than
Tate--Shafarevich classes.

Equivalently, the constructive target inside the residual Selmer quotient is
the globally soluble subspace

\[
 W_G=
 \frac{\delta(E(\mathbf Q)/2E(\mathbf Q))}{\delta(G/2G)},
 \qquad
 \dim_{\mathbf F_2}W_G=\operatorname{rank}E(\mathbf Q)-r.
\]

A Selmer class lies in this subspace exactly when its associated 2-covering
has a rational point.  The missing reverse implication is therefore to make
`W_G` large, not merely to make the ambient residual Selmer quotient large.

The continuation is now the canonical
[rational-solubility and residual-Selmer theorem package](RATIONAL_SOLUBILITY_AND_RESIDUAL_SELMER_THEOREMS.md).
It supplies the saturation correction, the precise Cassels--Tate radical
quotient `2 Sha[4]`, the distinction between midpoint charts and covering
maps, and exact soluble-cover certificates on the eleven requested fibres.

## Alternate-Q80 basis audit

Curve 12 requires a genuine basis change; its first seventeen published
points are not the alternate-Q80 generic basis.  The replay imports the exact
29-by-17 specialization matrix and adjoins the certified quotient basis

```text
P2, P11, P4, P3, P6, P8, P17, P10, P28, P24, P19, P15.
```

The resulting 29-by-29 matrix has determinant `-1`.  The 29 public
good-reduction signatures are transported through this matrix and retain full
binary rank 29.  The same half-ideal audit is then run on the ordered exact
basis

```text
G1,...,G17,Q1,...,Q12.
```

For all twelve `Q` directions the certificate records an explicit `G`
adjustment whose valuation parity is zero at every prime above every bad
rational prime.  It also stores the 29 exact integral point half-ideals

```text
A_P = (d^2*(4*x(P)-zeta), d^3*4*(2*y(P)+a1*x(P)+a3)),
d^2 = denominator(4*x(P)),
```

and verifies that each square correction is supported only above the declared
bad primes.

## Interpretation and next experiment

The theorem supplies a localization layer for the failed full-BNF front door.
A `+12` record does not merely coexist with a large auxiliary class group: its
certified exceptional Kummer classes already force ten or eleven independent
residual 2-class directions after the generic contribution is removed.  Native
alternate Q80 reproduces the phenomenon, so it is not an artifact of the
`074d9` chart.  The strict separation of the observed jump strata is largely a
formal consequence of the zero valuation-rank increment and the unit
correction; it is not a separate mechanism for the jumps.

A frozen prospective class/S-class feature could still be useful for candidate
scheduling, but another scalar class-group statistic cannot by itself supply
the missing implication.  The prospective small-field laboratory makes the
explicitly different localized choice

\[
 A_S=\operatorname{Cl}(K)/\langle S\rangle,
 \qquad
 A_S[2]/\langle c_S(G)\rangle.
\]

This is not `Pi(B,G)` above.  Reducing original class coordinates modulo two
does not compute it in the presence of higher 2-power torsion, and localizing
at `S` must not be identified with first taking the valuation kernel in the
full class group.  Any association found for this localized feature therefore
has its own interpretation.

The constructive experiment should instead:

1. compute a certified residual 2-Selmer basis modulo the specialized generic
   subgroup;
2. materialize basis classes and selected compatible combinations as explicit
   2-coverings, with complete local-solubility certificates;
3. use certified Cassels--Tate information as an obstruction and prioritization
   layer, without promoting its radical to rational points; and
4. search the compatible coverings for rational points, map every witness back
   to `E(Q)`, certify independence modulo the growing known subgroup, and
   iterate.

For an `MW17` rank-32 target, success means producing fifteen independent
rational directions, not merely exhibiting a residual Selmer space of
dimension at least fifteen.  The present six rows support neither that
construction nor a population-level correlation, Selmer upper bound, or
rank-search gate.

## Replay

```bash
sage -python \
  elkies-k3/scripts/certify_r17_kummer_classgroup_pressure_comparison.sage \
  --check
```

The replay first recomputes and compares the complete five-fibre v1 pressure
certificate, then performs the unimodular curve-12 basis change and the same
exact cubic/half-ideal calculation.  It uses Sage 10.9 and PARI through Sage.
