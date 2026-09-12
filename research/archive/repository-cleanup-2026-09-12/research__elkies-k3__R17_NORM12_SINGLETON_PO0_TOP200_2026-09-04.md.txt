# Alternate-Q80 direct-singleton `P.O=0` screen (2026-09-04)

<!-- status-consumer: EC-K3-R17-NORM12-11952-SINGLETON-PO0-TOP200 80ab545a98b4e2d7 -->

## Status

This is an exact bounded negative result. It does **not** prove generic rank
at least 19 or 20.

The complete 39,147-class smooth rational-bisection search gives one known
anti-invariant section on each direct quadratic twist. A finite-prime trace
census ranked those characters for a possible second section. For each of
the discovery top 200, two distinct usable finite fields were chosen and the
complete polynomial section box

\[
P\mathbin{.}O=0,\qquad \deg X\leq 6,\qquad \deg Y\leq 9
\]

was enumerated on the arithmetic-genus-three twist. The 400 selected shells
test 10,690,517,260 `X` polynomials and contain 51,410 representative-sign
solutions. After removing the reduction of the known section, 2,423 modular
points have full tangent rank eight. Exact digit-by-digit Hensel lifting
obstructs all 2,423: 655 first fail modulo `p^2`, 800 modulo `p^3`, 935 modulo
`p^4`, 31 modulo `p^5`, and 2 modulo `p^6`.

The consolidated two-prime certificate is
[`elkies-k3-r17-norm12-11952-singleton-twist-po0-two-prime-top200-audit-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-singleton-twist-po0-two-prime-top200-audit-v1.json).
The top-150 and top-20 certificates remain smaller positive-control and
regression artifacts.

## Exactness and positive control

For a twist of arithmetic genus `chi=3`, a section disjoint from the zero
section has polynomial coordinates in precisely the displayed degree box.
The finite-field enumerator fixes the smooth point met at infinity, recovers
the lower `Y` coefficients recursively, and exhausts the remaining `X`
coefficients. Each reported full-rank point is checked against every
coefficient equation over the integers. Failure to extend a solution from
`p^e` to `p^(e+1)` is therefore an exact local obstruction to a rational
section reducing to that point.

The first/second prime-pair multiplicities are 187 for `(17,23)`, 2 for
`(17,29)`, 1 for `(17,31)`, 4 for `(17,37)`, 1 for `(17,41)`, and 5 for
`(23,37)`. The campaigns hash-pin every export, known section, brute-force
shell, and Hensel audit, and checkpoint after every character. Eleven
characters have no isolated extra point at either prime; 107 have isolated
extras at both primes, all of which are locally obstructed.

As a positive control, the known section on `alternate-orbit-1c3d5` survives
511 successive lifting steps from `p=23` through `23^512`. Rational
reconstruction followed by substitution into the original twist equation is
literally exact. This shows that the lifting path recovers a genuine section
rather than merely rejecting spurious reductions.

## Heuristic holdouts

The first twenty characters were rescored on two disjoint 48-prime blocks,
`499--821` and `823--1151`. Their best weakest-block scores are respectively
`0.740349` and `0.784990`; none exhibits a stable rank-two signal. These
scores are target-selection heuristics and are not used in the exact local
exclusions.

## Replay

```bash
.venv/bin/python \
  elkies-k3/scripts/audit_r17_norm12_11952_singleton_po0_two_prime_top200.py \
  --check

.venv/bin/python \
  elkies-k3/scripts/audit_r17_norm12_11952_singleton_po0_top20.py \
  --check
```

## Boundary

This closes two distinct `P.O=0` finite-field shells for each of 200
heuristically selected characters. It does not exclude singular modular
branches: a characteristic-zero section could reduce into the known-section
or singular locus at both selected primes. It also does not exclude the
other 38,947 direct singleton characters, sections intersecting the zero
section, or characters not arising from one smooth rational bisection. In
particular, the requested second independent anti-invariant section remains
`UNKNOWN`.
