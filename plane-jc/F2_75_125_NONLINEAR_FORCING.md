# Exact nonlinear forcing on the F2 double-carrier component

## Result and claim boundary

The compiler
[`cas/compile_f2_75_125_nonlinear_forcing.py`](cas/compile_f2_75_125_nonlinear_forcing.py)
performs the nonlinear continuation left open by
[`F2_75_125_CARRIER_SPECIALIZATIONS.md`](F2_75_125_CARRIER_SPECIALIZATIONS.md).
It works on the two conjugate double-carrier points

\[
 R(w)=\frac{(w-\rho)^2}{25(1-\rho)^2},\qquad
 \rho^2-3\rho+1=0,                                      \tag{1}
\]

and presents the descent-eight ratio relatively by

\[
 27y^2-9y+1=0.                                           \tag{2}
\]

Thus the coefficient algebra is the quartic compositum
`QQ(rho)[y]/(27*y^2-9*y+1)`.  Equivalently, with
`theta=rho+y`, its primitive polynomial is

\[
729\theta^4-4860\theta^3+10341\theta^2
-7470\theta+1756.                                       \tag{3}
\]

The result is an exact arithmetic-circuit presentation of the component
ideal.  It is not a Gröbner calculation, a unit-ideal certificate, an F2
exclusion, or a Keller map.  The pinned machine record is
[`../artifacts/generated-results/jc2_f2_75_125_nonlinear_forcing.json`](../artifacts/generated-results/jc2_f2_75_125_nonlinear_forcing.json).

## 1. Source substitutions

The ten fixed-endpoint pivots are rebuilt from their exact `10 x 10` system
of determinant `75000`.  Their `1,489`-term solutions are stored as an
opaque, hash-consed arithmetic DAG rather than expanded into a multivariate
coefficient system.  The degree-seven follower lifts preserve `H(0)`; after
the complete band specialization the latter is identically zero.

Before the coupled rows, descents `1..11` are propagated by the exact tangent
solution

\[
 Q_{25-\delta}=-3C_0^2P_{15-\delta}.                    \tag{4}
\]

The resonant descents `5,10` retain the source centralizers `C0^4,C0^3`.
This is the already-certified rank-`134` upper elimination, now used as an
actual substitution in every later forcing vector.

## 2. The coupled 294+53 split

For each

\[
 \delta=12,13,\ldots,35,37,                             \tag{5}
\]

the compiler forms the complete bracket forcing after (4) and the ten
endpoint substitutions.  The top/new-`Q` image has a forced Laurent divisor.
Euclidean division by that divisor gives two canonical obstruction blocks:

1. the remainder, recording polynomial divisibility and local-jet
   compatibility; and
2. the quotient, projected into the pinned left cokernel of the normalized
   three-term operator.

The exact direct-sum ledger is

\[
 \boxed{294\text{ divisibility coordinates}
       +53\text{ quotient cokernel coordinates}=347}.   \tag{6}
\]

The numeral `53` alone is therefore not the full Laurent obstruction.  It is
the requested quotient block after forced-factor removal.  The other `294`
coordinates cannot be dropped: they say that the nonlinear forcing actually
admits that division.  The artifact pins the operator, cokernel-basis, and
equation digest separately at every descent.

## 3. Final functionals and incidence

For the target polynomial `G=J_4-1`, the compiler appends

\[
 G^{(j)}(\rho),\quad 0\le j\le4,                        \tag{7}
\]

and the two quotient residues of

\[
 \mathcal N(S)=5wE S'+(11E+22wE')S,qquad
 E=(w-1)(w-\rho).                                      \tag{8}
\]

These are the seven residual target coordinates.  The layer-zero block adds

\[
 H(\rho)-H(0),\qquad H^{(j)}(\rho),\quad1\le j\le5,     \tag{9}
\]

giving the promised `7+6=13` final functionals.

The descent-eight component itself is retained by five incidence equations:

\[
\begin{aligned}
&K_{P_7}(\rho)=K'_{P_7}(\rho)=K''_{P_7}(\rho)=0,\\
&K_{Q_1}(\rho)=0,\\
&K_{P_{-1}}(\rho)=25^3\rho(\rho-1)^6y
 \left(K'''_{P_7}(\rho)/6\right)^2.
\end{aligned}                                           \tag{10}
\]

Together with (2), the displayed circuit ideal has `347+13+5+1=366`
equations in `954` reachable source/circuit variables.  The exact DAG has
`710,439` reachable nodes and total degree at most seven on this carrier
specialization.  These are presentation counts, not dimension or
independence claims.

## 4. Squarefree first-defect routing

The rational carrier

\[
R_{\rm sf}(w)=\frac{w^2-3w+3}{25}                       \tag{11}
\]

has two simple roots and no movable double prime, so (10) is not its
component.  The compiler now routes it explicitly through every later
first-defect spacing `9..90`:

| spacing | primitive multiples strictly before target | count |
| --- | --- | ---: |
| `9..11` | `2,3` | 3 |
| `12..17` | `2` | 6 |
| `18..90` | none | 73 |

This is an exact `82`-row routing ledger, not an elimination of those rows.
Each spacing still needs its target-and-tail Fitting compiler; in particular,
failure of the double-carrier circuit ideal would not by itself exclude the
squarefree carrier.

## 5. Good-reduction tangent frontier

The follow-up driver
[`cas/probe_f2_75_125_nonlinear_modular.py`](cas/probe_f2_75_125_nonlinear_modular.py)
base-changes the full circuit to the split good reduction

\[
 \mathbf F_{31},\qquad \rho=14,\qquad y=3.             \tag{12}
\]

Unlike the component presentation above, this probe imposes the open
condition

\[
 a=K'''_{P_7}(\rho)/6\ne0                              \tag{13}
\]

by the explicit Rabinowitsch equation `a*localize_a_inverse-1=0`.  The seed
`K_P7=(w-rho)^3` has `a=1`.  Of the resulting `367` equations, exactly `41`
are initially nonzero: `35` coupled Laurent coordinates, five target/Hermite
coordinates, and one incidence coordinate.  The full sparse Jacobian has
rank `214`, and its inhomogeneous tangent system is consistent.

The minimal-support tangent correction exposes a much smaller natural chart:

\[
\begin{aligned}
P&:\ 11,7,3,-1,-5,-9,-13,-17,\\
Q&:\ 13,9,5,1,-3.
\end{aligned}                                           \tag{14}
\]

This spacing-four staircase has `169` coordinates including `y` and the
localizing inverse.  Its restricted Jacobian has rank `57`, so the consistent
affine tangent fiber has dimension `112`.  The pinned particular correction
uses only `45` coordinates.

A same-field Newton step is not a Hensel lift.  The driver therefore evaluates
the entire affine line through that particular correction.  All `366`
geometric equations restrict to polynomials of degree at most two on this
line; `64` restrictions are nonzero and their exact gcd in
`F_31[t]` is `1`.  Hence this particular line is disjoint from the nonlinear
zero locus even over the algebraic closure.  It is not merely a failed list
of `F_31` trials.

The corrected sparse solver back-substitutes pivots in descending column
order.  With that correction, the full Jacobian has a pinned
`153`-dimensional left cokernel and the particular line projects identically
to zero.  The complete directional compiler then parameterizes the affine
tangent fiber as `d0+u*k` and interpolates every one of its `112` coordinate
lines through the exact degree-seven bound.  All `112` projected restrictions
are zero.  Eight deterministic dense mixed lines are also identically zero
in all `153` coordinates.  Thus no first cokernel obstruction is detected on
the staircase chart.  Coordinate lines plus eight dense lines do not prove
that the full multivariate obstruction map vanishes, but they remove the
previously suspected spacing-eight quadratic block.

## 6. Fixed-Jacobian formal homotopy

The next compiler replaces same-field Newton iteration by the exact formal
homotopy

\[
 F(x(\lambda))=(1-\lambda)F(x_0).                       \tag{15}
\]

At order one this is exactly the inhomogeneous tangent equation.  At every
higher order the compiler evaluates the full arithmetic circuit in
`F_p[lambda]/(lambda^(N+1))`, projects the new forcing to the fixed Jacobian
cokernel, and solves against the same rank-`214` Jacobian only when that
projection vanishes.

Over `F_31`, with all `112` first-order free parameters set to zero, the path
lifts without obstruction through order `16`.  Every order-two-through-sixteen
forcing projects to zero in all `153` cokernel coordinates.  The localization
series is exactly

\[
 a(\lambda)=1.                                         \tag{16}
\]

The same construction over the independent good split reduction
`F_61`, `(rho,y)=(19,19)`, has the same Jacobian rank and cokernel dimension
and lifts through order `8`.  The smaller split prime `19` is rejected before
evaluation because a pinned rational circuit constant has denominator
divisible by `19`.

In the default pivot gauge, seven coefficient sequences exhibit the same
finite-prefix recurrence with denominator `(1-lambda)^2` at both good primes:
`P_-5_d0` and `P_7_d2..d7`.  This is a gauge effect.  Prescribing those seven
higher-order coefficients to be zero remains consistent and produces a
regular-gauge order-`16` jet with smaller correction supports.  Directly
substituting `lambda=1` into the truncated regular jet is not a solution; it
leaves `50` nonzero equations over `F_31`.  Thus the computation proves a long
unobstructed formal deformation, not convergence or a modular point.

This finite-field calculation proves neither emptiness nor nonemptiness of
the localized component.  It does change the computational frontier.  The
next target is continuation of the regular formal branch to `lambda=1`:
either reconstruct an algebraic/rational parameterization from longer jets,
or change the higher-order kernel gauge to obtain a branch regular at the
target.  Blind unit Newton iteration and a Gröbner basis in all `954`
reachable source variables remain unsupported by the data.

The pinned records are
[`../artifacts/generated-results/jc2_f2_75_125_modular_probe.json`](../artifacts/generated-results/jc2_f2_75_125_modular_probe.json),
[`../artifacts/generated-results/jc2_f2_75_125_tangent_obstruction.json`](../artifacts/generated-results/jc2_f2_75_125_tangent_obstruction.json),
[`../artifacts/generated-results/jc2_f2_75_125_formal_homotopy.json`](../artifacts/generated-results/jc2_f2_75_125_formal_homotopy.json),
[`../artifacts/generated-results/jc2_f2_75_125_formal_homotopy_regular_gauge.json`](../artifacts/generated-results/jc2_f2_75_125_formal_homotopy_regular_gauge.json),
and
[`../artifacts/generated-results/jc2_f2_75_125_formal_homotopy_mod61.json`](../artifacts/generated-results/jc2_f2_75_125_formal_homotopy_mod61.json).

## 7. Reproduction and next operation

Run

```bash
.venv/bin/python plane-jc/cas/compile_f2_75_125_nonlinear_forcing.py
.venv/bin/python plane-jc/cas/test_f2_75_125_nonlinear_forcing.py
.venv/bin/python plane-jc/cas/test_sparse_circuit_modp.py
.venv/bin/python plane-jc/cas/probe_f2_75_125_nonlinear_modular.py
.venv/bin/python plane-jc/cas/compile_f2_75_125_tangent_obstruction.py
.venv/bin/python plane-jc/cas/compile_f2_75_125_formal_homotopy.py
.venv/bin/python plane-jc/cas/compile_f2_75_125_formal_homotopy.py --regular-gauge --artifact artifacts/generated-results/jc2_f2_75_125_formal_homotopy_regular_gauge.json
.venv/bin/python plane-jc/cas/compile_f2_75_125_formal_homotopy.py --prime 61 --rho 19 --y 19 --maximum-order 8 --artifact artifacts/generated-results/jc2_f2_75_125_formal_homotopy_mod61.json
```

Intentional artifact regeneration uses `--refresh` on either artifact-producing
command.  The full directional and order-16 formal audits take several
minutes because they replay the complete arithmetic circuit many times.
