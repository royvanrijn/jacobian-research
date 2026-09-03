# Alternate Q80 Rootless Equation Handoff (2026-09-03)

<!-- status-consumer: EC-K3-H3-ROOTLESS-J2-MINIMAL-ACCESSIBILITY 631f50389e0a3283 -->

## Scope and status

This note pins the exact starting point for recovering the characteristic-zero
equation of one historically transported copy of the unique alternate
rootless determinant-948 `J2` frame.  It does not alter the completed two-class
rootless-frame classification and does not claim an equation-level fibration.

**Updated construction priority.**  The historical Q80 lift below is no longer
the preferred route to the frame class.  The exact classifier in
[`J2_GEOMETRIC_ACCESSIBILITY_2026-09-03.md`](J2_GEOMETRIC_ACCESSIBILITY_2026-09-03.md)
finds ten nef degree-two copies of the alternate class directly among the 43
norm-twelve genus-one bisections on the published R17 equation.  Each shares
the old zero and has zero-section degree one.  The cheapest stored witness is
`norm12-orbit-11952`.  Compile its pencil `|D|` before resuming the million-bit
third-q12 reconstruction recorded here.

The lattice handoff is complete and independently replayable.  Along this
historical route, the equation breakthrough (Level A) remains open because
there is no characteristic-zero equation for the immediate `A1/MW16` parent.
For the new direct route, the open gate is instead the resolved two-dimensional
Riemann--Roch pencil of one exact genus-one divisor.  The new Weierstrass model,
its 17 explicit sections, and all specialization comparisons remain open.

## Exact final divisor

Use row coordinates in

```text
NS_parent = U + (-M_A1/MW16).
```

The raw final `q=6`, `(a,b)=(2,3)` isotropic class is

```text
(2,3,0,-2,4,2,-1,2,1,-1,1,0,1,-1,1,0,0,0,0).
```

Reflection first in the old zero and then in the sole simple `A1` component
gives the physical nef fibre

```text
F_alt = (3,2,-1,-2,4,2,-1,2,1,-1,1,0,1,-1,1,0,0,0,0).
```

Applying the same Weyl reflections to the canonical child zero gives

```text
O_alt = (15,10,-5,-10,20,10,-5,9,5,-5,5,0,5,-5,5,0,0,0,0).
```

The old-horizontal section singled out by the divisor is

```text
S = (5,1,-1,-2,4,2,-1,2,1,-1,1,0,1,-1,1,0,0,0,0),
F_alt = O_old + S - F_old.
```

The checker proves exactly:

- `F_alt^2=0`, `F_alt` is primitive and nef;
- `O_alt^2=-2` and `F_alt.O_alt=1`;
- the minimum pairing with any old section is one;
- a negative bisection is impossible;
- the child frame is positive definite, rootless, rank 17, and determinant
  948;
- it has 2,626 oriented norm-four vectors, hence 1,313 pairs;
- it is not integrally isometric to the published R17 frame;
- the complete transport back to the pinned Q80 Neron--Severi basis is
  integral of determinant one.

The machine-readable source of these coordinates and matrices is
[`../artifacts/generated-results/elkies-k3-q80-alternate-final-divisor-handoff-v1.json`](../artifacts/generated-results/elkies-k3-q80-alternate-final-divisor-handoff-v1.json).

## What the Riemann--Roch result does and does not provide

The primitive nef isotropic K3 theorem gives

```text
h^0(X,O(F_alt)) = 2.
```

On an explicit equation for the parent, the divisor identity predicts the
formal generic-fibre basis

```text
1, (y+y(S))/(x-x(S)),
```

after the certified vertical `F_old` cancellation.  This is not yet a
resolved-surface valuation proof: neither the parent equation nor the
coordinates of `S` on it have been compiled.

The two-dimensional kernel previously constructed over `GF(73)` belongs to a
special CM24 shadow with child root type `A1+2A3+2A4/MW3`.  It must not be used
as the generic alternate rootless pencil.

## Correct upstream frontier

The equation pipeline must start earlier than the final q6 transition:

```text
exact third-q12 genus-one pencil
    -> exact third-q12 Jacobian and two-way maps
    -> fourth q12 equation
    -> alternate q4 A1/MW16 equation
    -> final q6 rootless equation.
```

At fixed `u=-2`, the third-q12 pencil is exact over the biquadratic field
`QQ(a,b)`, descends to `QQ(a*b)`, has a 63-term moving equation and a certified
smooth genus-one generic member.  The long Jacobian and maps are known over
many finite and p-adic precisions, but their attempted rational
reconstructions fail held-out primes.  No exact characteristic-zero child
Jacobian has passed the fail-closed gate.

An exact normalization audit now observes that for the reduced fraction

```text
q1*q2 = N/D,
```

`N` is a perfect square and `D` is not.  Therefore

```text
QQ(sqrt(q1*q2)) = QQ(sqrt(D)).
```

This permits the denominator-integral generator `delta^2=D` in place of
`omega^2=16*q1*q2`, without factoring the 33,894-bit denominator.  It is an
exact compiler normalization, not evidence that the resulting equation has
small coefficients.  The replay is
[`../artifacts/generated-results/elkies-k3-q80-third-q12-descent-field-normalization-v1.json`](../artifacts/generated-results/elkies-k3-q80-third-q12-descent-field-normalization-v1.json).

The proposed basis-height gate has also been run on all 63 coefficients of
the exact moving equation.  Before global model rescaling, the maximum
rational-coordinate heights are:

| basis | maximum bits | median bits |
|---|---:|---:|
| `theta^2` | 1,462,954 | 1,409,636 |
| `omega` | 1,456,899 | 1,403,570 |
| `delta` | 1,462,828 | 1,409,537 |

Relative to `omega`, `delta` is lower on 9 terms and higher on 54, increasing
the maximum by 5,929 bits.  After one rational projective normalization,
`delta` instead lowers the primitive maximum from 1,495,639 to 1,484,751 bits,
about 0.7 percent.  The cleared integer contents have only 7--12 bits, so
there is no large missed rational common factor in any of the three bases.
This is a bounded basis comparison, not an optimization over multiplication
by a general quadratic-field element, base `PGL2`, Weierstrass changes, or
integral ideals.  Its exact output is
[`../artifacts/generated-results/elkies-k3-q80-third-q12-pencil-basis-heights-v1.json`](../artifacts/generated-results/elkies-k3-q80-third-q12-pencil-basis-heights-v1.json).

<!-- status-consumer: EC-K3-Q80-THIRD-Q12-EXACT-LINEAR-CONDUCTOR 957479f39bedd57b -->

## Exact conductor progress

The first characteristic-zero factor of the third-q12 discriminant is now
recovered exactly.  Write the old base variables as `V,W`.  Rational
reconstruction from the `19^12288` factor lift gives a constant `r` with
6,059-bit numerator and 6,053-bit denominator.  A direct integral compilation
of the exact 63-term pencil in the shifted variable

```text
S = denominator(r)*W + numerator(r)
```

proves that the coefficients of `S^0,S^1,S^2` in the discriminant vanish
identically in `QQ(sqrt(q1*q2))[V]`, while the `S^3` coefficient is a nonzero
degree-eight polynomial.  Consequently

```text
(W+r)^3 | discriminant(V,W)
```

with exact multiplicity three over the full characteristic-zero base.  The
nonzero coefficient has nine terms, maximum coordinate height 6,138,563 bits,
and SHA-256
`9d936d08118c94e9d1404b6dca3b8f88500c7a3554a0a5407995c6a0f9915091`.
The exact replay is
[`../artifacts/generated-results/elkies-k3-q80-third-q12-exact-generic-linear-conductor-v1.json`](../artifacts/generated-results/elkies-k3-q80-third-q12-exact-generic-linear-conductor-v1.json).

At `V=0`, the exact discriminant has degree 15 and maximum coordinate height
5,783,744 bits.  Removing the proved cubic linear factor leaves degree 12 at
5,783,726 bits.  Both PARI's exact gcd and a Singular number-field gcd were
tried on this residual; neither completed, including a five-minute Singular
run.  This is the present exact arithmetic cost wall, not evidence against
the expected quartic square factor.

There is also a useful, deliberately unpromoted normalization for that
quartic.  In the `19^12288` lift, all four nonleading coefficients of the
monic quartic `Q` have one common monic linear denominator `H(V)`.  Its
constant has coordinate heights 10,904/10,897 bits in the rational coordinate
and 6/16,947 bits in the `delta` coordinate.  Direct reductions of the exact
117 MB pencil at the untouched inert primes 163, 191, and 199 each reproduce
the factor pattern `(1,3),(4,1),(4,2)` and the predicted `H` as all four
quartic denominators.  The audit is
[`../artifacts/generated-results/elkies-k3-q80-third-q12-quartic-denominator-candidate-v1.json`](../artifacts/generated-results/elkies-k3-q80-third-q12-quartic-denominator-candidate-v1.json).
This is three-prime evidence for the normalization only: exact
characteristic-zero `Q^2` divisibility, `Q`, the remaining factor, the
Jacobian, and its maps are still open.

<!-- status-consumer: EC-K3-Q80-THIRD-Q12-EXACT-SPECIALIZED-QUARTICS 725664f9e36ae8a7 -->

## Last exact-Q run checkpoint: 2026-09-03

The current `19^12288` quartic reconstruction space has now been exhausted
systematically in both the `omega` and `delta` bases.  The fail-closed probe

```text
probe_q80_third_q12_quartic_projective_reconstruction.sage
```

tests four reconstruction granularities against direct exact-pencil
factorizations at all of `p=163,191,199`:

| reconstruction block | result at both bases |
|---|---|
| each scalar, over its full asymmetric Euclidean-convergent chain | `0/16` coordinates accepted |
| each quadratic-coordinate pair | no LLL row accepted |
| each four-coordinate `W` coefficient | no LLL row accepted |
| all sixteen numerator coordinates jointly | no LLL row accepted |

The shortest rows occur at their random-lattice boundaries: about 34,800 bits
for a coordinate pair, 41,760 bits for one `W` coefficient, and 49,129 bits
for the joint block.  No row passes even the first held-out comparison.  Thus
there is no justification for promoting any current-precision rational or
projective reconstruction candidate.

The exact residual itself reduces to gcd degree four at each of the same three
untouched primes.  A custom quadratic-pair implementation of Brown's
subresultant PRS has an internal exact `Q^2 D` self-test that terminates at
degree four.  An earlier false degree-zero result was traced to Python `/` on
`gmpy2.mpz`, which produced approximate `mpfr` values; explicit `mpq`
inversion fixes that correctness bug.

The corrected characteristic-zero run had completed
the following exact nonzero subresultants after removal of `(W+r)^3` at
`V=0`:

| degree | maximum coordinate bits | preceding scaled-step time |
|---:|---:|---:|
| 10 | 16,894,494 | initial pseudo-remainder |
| 9 | 28,160,280 | 111.834 seconds |
| 8 | 39,426,096 | 328.980 seconds |
| 7 | 50,691,932 | 554.864 seconds |
| 6 | 61,957,773 | 725.906 seconds |
| 5 | 73,223,619 | not retained |
| 4 | 84,489,468 | 818.837 seconds |

The run was then observed at 100 percent CPU and about 1.90 GB RSS computing
the degree-six-to-five step.  At the next status check PID `1740880` no longer
existed, its detached terminal session was unavailable, and no result artifact
had been written.  There is no kernel OOM record for that PID, but its exit
status cannot be recovered.  Therefore nothing below degree six from that
particular run was retained.

The recovery worker now writes an atomic binary Brown-state checkpoint after
every completed exact subresultant, with a human-readable JSON companion and
input hashes for the operands, pencil, factor lift, base value, descent field,
and stripped residual.  A fresh replay has exercised this on the real
million-bit input and durably recovered the complete sequence through degree
four.  Its final checkpoint is
`artifacts/local/elkies-k3/q80-third-q12-exact-quartic-subresultant-checkpoint-v1.pickle`.
The resulting monic degree-four gcd `Q` has maximum coordinate height 320,859
bits.  Literal exact division gives a monic quartic `D` of maximum coordinate
height 1,735,258 bits and proves

```text
monic(Delta(V=0)/(W+r)^3) = Q(W)^2 D(W)
```

over the exact quadratic descent field.  The 18 MB certificate is
[`../artifacts/generated-results/elkies-k3-q80-third-q12-exact-discriminant-specialization-v1.json`](../artifacts/generated-results/elkies-k3-q80-third-q12-exact-discriminant-specialization-v1.json),
with SHA-256
`96025b0829943030150925d9911c753261a08c71e8ca48b4f6b878f34badded6`.
An independent `--check` replay loaded the completed checkpoint, repeated the
exact `Q^2` division, and returned `PASS_CHECK` with the same digest.  This is
an exact characteristic-zero theorem at `V=0`; it does not yet recover the
generic quartics over the full `V`-line or the Jacobian.

## Next gates

Retain the existing `omega` presentation for the current p-adic pipeline.
The small rational-projective benefit from `delta` does not compensate for
the worse raw coordinates or reduce the roughly 1.5-million-bit primitive
height to a tractable scale.  General two-coordinate `K^*` scaling, nine
projective coefficient charts, the `I6/I4` and `I2`-trace base gauges, and
joint evaluation lattices have already failed the independent `p=199` replay;
they should not simply be repeated in the `delta` basis.  Mere rational
content extraction is closed by the 7--12-bit content calculation.

The specialized degree-four gcd and exact square division are now certified.
The next justified compiler should use this exact `V=0` quartic together with
the candidate common denominator `H(V)` and the untouched-prime factorizations
to recover the generic `Q(V,W)` projectively.  It must finish by literal exact
division in the original two-variable characteristic-zero discriminant.  A
general expanded exact gcd and another unstructured high-precision `j`-map LLL
run have both reached their useful limits.  In parallel, the exact specialized
factors now provide a normalization anchor for the invariant-first computation

```text
exact factored pencil -> j -> (c4^3,Delta) -> minimal Jacobian -> maps.
```

An implementation should operate on factored or lazy rational coefficients,
not expand the 1.5-million-bit primitive coefficient vector, and must replay
the resulting invariant at the untouched primes `163`, `191`, and `199`.

After the exact third-q12 Jacobian is obtained, compile each remaining
neighbour in order.  At every stage retain the source equation, zero,
resolved fibre components, two-way maps, and NS marking.  The final rootless
equation must have root rank zero and geometric Mordell--Weil rank 17 before
section transport begins.

Hypotheses H1--H5 about multisection visibility, quotient minima, rank-28
exceptional directions, heuristic enrichment, and specialization tails are
not yet testable.  They require the equation and saturated 17-section basis;
no arithmetic-search conclusion is drawn here.

## Replay

```bash
/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/build_q80_alternate_final_divisor_handoff.sage --check

python3 \
  elkies-k3/scripts/audit_q80_third_q12_descent_field_normalization.py --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_q80_third_q12_pencil_basis_heights.py --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_q80_third_q12_exact_discriminant_specialization.sage \
  --base-value 0 --certify-generic-linear --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/probe_q80_third_q12_exact_discriminant_specialization.sage \
  --attempt-subresultant-prs \
  --output artifacts/generated-results/elkies-k3-q80-third-q12-exact-discriminant-specialization-v1.json \
  --check

/home/royvanrijn/.local/share/jacobian-sage-10.9/bin/python \
  elkies-k3/scripts/audit_q80_third_q12_quartic_denominator_candidate.sage \
  --check
```
