# Stratified adelic engineering in seed and target space

The fixed-seed adelic theorem prescribes a real chamber and finitely many
unramified factorization types by moving only the target coordinates `(s,t)`.
This note records the corresponding statement when the seed is also allowed
to move on a prescribed locally closed stratum. It applies to exact
omitted-value components, their complements, Hessian-symmetry strata, and
their intersections whenever the selected intersection has a rational chart.

The conclusion is deliberately local-to-global. It does not assert that an
arbitrary omitted-value component meets an arbitrary Hessian-symmetry stratum
or an arbitrary real chamber. Those are geometric nonemptiness questions.
Once compatible local points exist on one rational chart, the real and
finite-prime fiber conditions impose no further global coupling.

## 1. The universal seed--target space

Let

\[
 \mathcal A_n=\{H:H(0)=H'(0)=H(1)=0,\ H'(1)=-1,
                         \ H''(1)\ne-2\}.
\]

It is an open subset of affine `(n-3)`-space. Over
`mathcal A_n x A^2_(s,t)` put

\[
 E_{H,s,t}(W)=H(W)-sW+t.                              \tag{1}
\]

For a partition `lambda` of `n`, let `U_lambda` be the locally closed locus
where (1) is squarefree and has factor-degree partition `lambda` over the
residue field. Over the reals let `U_j(R)` be the open locus where (1) has
exactly `j` real roots. Squarefreeness makes the corresponding weighted
Keller fiber complete and regular, and its fiber algebra is

\[
 \mathbb Q[W]/(E_{H,s,t}).                            \tag{2}
\]

Only the conjugacy class of Frobenius is intrinsic without a labeling of the
roots. Thus prescribing a "Frobenius permutation" below means prescribing
its cycle type.

## 2. Stratum-wise weak approximation theorem

Let `X subset A_n` be a locally closed subvariety over `Q`. Suppose a dense
open `X^o` has a rational coordinate chart

\[
 \theta:V\buildrel\sim\over\longrightarrow X^o,
 \qquad V\subset\mathbb A^d_{\mathbb Q}\text{ open}.             \tag{3}
\]

Fix a parity-compatible integer `j`, a finite set `S` of rational primes,
and a partition `lambda_p` of `n` for every `p in S`. Assume:

1. there is a real point `(u_infty,s_infty,t_infty)` of
   `V x A^2` at which (1), with `H=theta(u_infty)`, has `j` simple real
   roots; and
2. after choosing one integral model of (3), every `p in S` has a residue
   point `(u_p,s_p,t_p)` at which the chart is defined and (1) is squarefree
   of factor-degree type `lambda_p`.

Then there are `u in V(Q)` and `(s,t) in Q^2` such that

\[
 H=\theta(u)\in X^o(\mathbb Q)                       \tag{4}
\]

has exactly `j` real sheets over the target `(t,s,1)` and has Frobenius cycle
type `lambda_p` at every `p in S`. At each selected prime the fiber algebra
(2) is etale over `Z_p`, so the displayed factor degrees are the degrees of
its unramified local factors.

### Constructive proof

Shrink around the real point until a rational box

\[
 B_\infty\subset V(\mathbb R)\times\mathbb R^2       \tag{5}
\]

lies in `U_j(R)` and avoids every denominator and deleted divisor of the
chart. Clear the finitely many chart denominators. Apply the denominator-one
CRT grid simultaneously to all `d+2` coordinates, using the residue vector
`(u_p,s_p,t_p)` at `p`. A sufficiently fine grid meets (5). Its common
denominator is `1` modulo every selected prime, so evaluation of the rational
chart commutes with reduction and recovers every prescribed residue point.

The real box preserves the root count. Reduction to a squarefree polynomial
preserves the modular factorization type and also proves that the lifted point
does not lie on any deleted divisor visible at that prime. If `S` is empty,
rational density in (5), followed by avoidance of finitely many proper closed
subsets, gives the same conclusion. This proves the theorem.

This is weak approximation on an explicit rational chart, not an appeal to
weak approximation for arbitrary varieties. The rational-chart hypothesis
is therefore load-bearing.

## 3. Geometric strata covered by the theorem

### Omitted-component membership

For `2a+3b=n`, the exact nonsurjective stratum of type `2^a3^b` has the
quotient-coordinate model

\[
 M=Q^2R^3,\qquad
 \Phi=M(1)-M(0)-M'(0)=0,                             \tag{6}
\]

with `deg Q=a`, `deg R=b`, followed by the open exactness and admissibility
conditions. Its seed is

\[
 H={-M+M'(0)W+M(0)\over M'(1)-M'(0)}.                \tag{7}
\]

On the exact stratum, unique factorization recovers `(Q,R)`, so (7) is a
rational chart wherever (6) is rationally parametrized. In the stable ranges
`a>=3` or `b>=3`, rationality is explicit: in endpoint coordinates the
equation is

\[
 X^2Y^3-x^2y^3-2xuy^3-3x^2y^2v=0,                  \tag{8}
\]

and one solves linearly for `u` or `v` on a dense open. The finitely many
low-rank cases may be handled by their individual charts. Therefore the
theorem constructs exact type-`2^a3^b` seeds with simultaneous real and
modular fiber data whenever the chosen chart has the required local points.

### Omitted-component nonmembership

The complement of a selected component is open in the affine seed space. It
is therefore already covered by the identity chart. The same applies to the
complement of the entire nonsurjective locus.

### Hessian symmetry

The exact Hessian-symmetry type is locally closed. For the coarse marked
strata whose symmetry fixes zero, it has the explicit form

\[
 H(W)=W^eR(W^q),\qquad R(1)=0,\qquad qR'(1)=-1,       \tag{9}
\]

with the larger-congruence subloci deleted to make the stabilizer exactly
`mu_q`. Equations (9) are affine-linear in the coefficients of `R`, hence
give rational charts. General centered symmetry strata are given by the same
Hessian-stabilizer equations with a moving center.

To impose omitted-component membership and a Hessian type simultaneously,
take `X` to be a rational locally closed chart in their intersection. The
theorem shows that the arithmetic fiber engineering is independent after
this geometric intersection is known to be nonempty. It does **not** turn an
empty or nonrational intersection into a nonempty one.

## 4. Explicit nonsurjective quintic with two Frobenius conditions

Consider

\[
 H(W)={1\over27}W^2(1-W)(W^2+8W+18).                \tag{10}
\]

It is a normalized admissible seed and

\[
 H(W)+W-1=-{1\over27}(W+3)^3(W-1)^2.                \tag{11}
\]

Thus its unique omitted target is `(s,t)=(-1,-1)` and its exact contact type
is `3+2`, the component `2^1 3^1`.

Its Hessian divisor is cut out, up to a nonzero scalar, by

\[
 (W+3)(5W^2+6W-3)=5W^3+21W^2+15W-9.                \tag{12}
\]

After translating its barycenter to zero, (12) becomes

\[
 5X^3-{72\over5}X-{64\over25}.                      \tag{13}
\]

Both lower coefficients are nonzero. A nontrivial affine stabilizer of a
reduced cubic divisor would force the centered cubic to have form `X^3+c`
(a `mu_3` symmetry) or `X(X^2+c)` (a `mu_2` symmetry). Hence this seed has
trivial Hessian stabilizer.

Prescribe the cycle types

\[
 \lambda_7=(5),\qquad \lambda_{11}=(2,2,1).          \tag{14}
\]

Local witnesses for the fixed seed are

\[
 (s,t)\equiv(0,4)\pmod7,\qquad
 (s,t)\equiv(1,1)\pmod{11}.                         \tag{15}
\]

A denominator-one CRT lift, followed by exact real-root certification, gives

\[
 \boxed{(s,t)=\left(-{161\over103},{48\over103}\right).}        \tag{16}
\]

In fact the same seed and the same two modular conditions realize every
allowed quintic real-root count:

| `j` | `s` | `t` |
|---:|---:|---:|
| `5` | `-161/103` | `48/103` |
| `3` | `-1435/963` | `-82/963` |
| `1` | `-1435/963` | `-223/107` |

Each row is squarefree and has exactly the displayed number of real roots.
For the five-real-root row, the primitive integral polynomial is, up to a
nonzero scalar,

\[
 P(W)=103W^5+721W^4+1030W^3-1854W^2-4347W-1296.     \tag{17}
\]

Modulo the selected primes, its squarefree factorizations have the forms

\[
 \begin{aligned}
 E(W)\bmod7&\sim W^5+3W^3+3W^2-3,
     &&\text{irreducible},\\
 E(W)\bmod11&\sim
 (W-1)(W^2-2W-1)(W^2-W-5).
 \end{aligned}                                      \tag{18}
\]

Therefore one nonsurjective seed of exact type `2^1 3^1` and trivial Hessian
symmetry has complete `5`-, `3`-, and `1`-real-sheet Keller fibers whose
Frobenius classes at `7` and `11` all have the prescribed cycle types. The
`5`-cycle at `7` also proves that every row defines a quintic field. Their
signatures are respectively `(5,0)`, `(3,1)`, and `(1,2)`; `11` has
unramified splitting type `(2,2,1)` in each field.

## Verification

Run

```bash
.venv/bin/python scripts/verify_stratified_adelic_engineering.py
```

The checker verifies normalization and admissibility, the exact omitted
factorization, trivial Hessian stabilizer, the three CRT lifts and exact real
root counts, both squarefree modular factorization types, and global
irreducibility.
