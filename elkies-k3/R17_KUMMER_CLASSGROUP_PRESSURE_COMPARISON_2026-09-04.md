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

This is a compelling explanation of the class-group wall: the known
exceptional points themselves force a large residual 2-class image before a
complete descent begins.  It is not yet a prospective rank predictor.  The
exceptional points, and hence the known jump, are inputs to the invariant.

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

## Interpretation and next test

The theorem supplies the missing explanatory layer for the failed full-BNF
front door.  A `+12` record does not merely coexist with a large auxiliary
class group: its certified exceptional Kummer classes already force ten or
eleven independent residual 2-class directions after the generic contribution
is removed.  Native alternate Q80 reproduces the phenomenon, so it is not an
artifact of the `074d9` chart.

What remains unproved is the predictive direction.  A valid next test must
construct a comparable class/S-class feature before revealing exceptional
points or the fibre's jump, then evaluate it on a frozen holdout panel.  The
present six rows cannot support a population-level correlation, a Selmer upper
bound, or a rank-search gate.

## Replay

```bash
sage -python \
  elkies-k3/scripts/certify_r17_kummer_classgroup_pressure_comparison.sage \
  --check
```

The replay first recomputes and compares the complete five-fibre v1 pressure
certificate, then performs the unimodular curve-12 basis change and the same
exact cubic/half-ideal calculation.  It uses Sage 10.9 and PARI through Sage.
