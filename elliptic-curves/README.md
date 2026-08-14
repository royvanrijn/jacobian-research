# Elliptic-curve rank/conductor programme

This is a research track separate from the Keller-map programme.  Its two
operational targets are:

1. exhibit an elliptic curve over \(\mathbb Q\) with **at least 21 independent
   rational points** and
   \(\log N<182.72\), where `log` is the natural logarithm; or
2. exhibit one with **at least 30 independent rational points**, improving the
   public rank-at-least-29 record.

The literal first cutoff is

```text
N <= 22609332114411420624526008180120083289443726642551045721326266745468787166869234.
```

The number `182.72` comes from the rounded value for Fermigier's rank-at-least-22
curve in Bober's table.  Its exact conductor has
`log(N) = 182.7249109506374287961...`, so that benchmark narrowly misses the
literal cutoff.  See [baselines and literature](notes/BASELINES_AND_LITERATURE.md).

## What is implemented

The first end-to-end calibration implements the proposed

```text
p-adic roots -> CRT -> exact 2D Gauss reduction -> small rational t
             -> integral curve -> minimal model/conductor/rank checks
```

for

\[
E_t:y^2=x^3-t^2x+t^2,
\qquad
\Delta(a,b)=-16a^4b^6(27b^2-4a^2),
\]

where \(t=a/b\).  Lifting roots at \(23^3,47^2,73^2\), considering all eight
root-sign combinations, and exhaustively searching the relevant rational-height
boxes gives

\[
t=-\frac{110627}{84367},\qquad
27b^2-4a^2=23^3\,47^2\,73^2.
\]

The resulting integral model is

\[
y^2=x^3-87109893594476435881x
 +620029989546545117143687312009.
\]

PARI/GP 2.15.4 verifies that the displayed model is already global minimal,
that the three shaped primes have conductor exponent one, that torsion is
trivial, and that `ellrank` returns bounds `[3,3]`.  The pinned manifest is
[here](../artifacts/generated-results/elliptic-curves/crt_lattice_calibration_v1.json).
This is an exact mechanism calibration, not a candidate for either target and
not an independent rank computation.

## First high-family seed

The programme now has an exact adapter for Fermigier's fixed-root Mestre
family.  Exact binary-quartic invariants connect its quartic construction to a
canonical model

```text
[1, a2(u), 1, a4(u), a6(u)],     u = s/2,
```

whose discriminant is one primitive even polynomial \(\Phi(u)\) of degree 20.
The apparent \(u^{12}\) and twelfth-power constant in the raw quartic
conversion are removable coordinate artifacts and are not used as shaping
factors.  Roots of \(\Phi\) at 89, 131, and 137 are simple and give split
multiplicative reduction.

Forcing all three primes to the second power, checking all eight root choices,
and exhaustively reconstructing rational parameters through the first occupied
height box gives

\[
u=\frac{673709}{29965},\qquad M=89^2 131^2 137^2.
\]

The homogeneous degree-20 discriminant factor has valuation exactly two at
each shaped prime.  PARI verifies split type \(I_2\), Tamagawa number 2, and
conductor exponent one at all three.  The
[pinned seed](../artifacts/generated-results/elliptic-curves/fermigier_crt_seed_v1.json)
does **not** claim a global conductor or a rank: its 469-bit remaining cofactor
is deliberately left unfactored.  Separately, the reconstructed thirteenth
quartic point and an exact finite-reduction certificate prove the twelve
generic section differences independent.  A second certificate proves all 22
published E22 points independent.  The
[pinned rank replay](../artifacts/generated-results/elliptic-curves/fermigier_rank_certificates_v1.json)
is a lower-bound certificate only: it supplies neither saturation nor an
unconditional rank upper bound.

There is also a reproducibility discrepancy in the source.  Fermigier prints
the symmetric shift `s=19754/39`, but that literal substitution gives a
different curve.  The displayed E22 curve is recovered exactly at
`s=39508/39`, or `u=19754/39`; the two literal specializations have different
\(j\)-invariants.  No erratum or derivation of this factor two was found, so the
coordinates remain explicitly separated rather than silently identified.

## Best low-conductor near miss

A bounded `ratpoints` search at

\[
u=\frac{28917}{20}
\]

found a point cloud from which exact finite-reduction arithmetic selects 20
independent points.  PARI gives the global minimal model

\[
y^2+xy+y=x^3+x^2
-4437412060110743641525245114305x
+3586842216822165612930264910099076801587288127
\]

with exact conductor

```text
2876153493562761211278364526603564191699143885403233935132057708367930
```

and `log(N) = 159.9348252255254533984...`, comfortably below the target
cutoff.  The [pinned near-miss certificate](../artifacts/generated-results/elliptic-curves/fermigier_rank20_near_miss_v1.json)
records the 58 searched abscissas, exact search bounds, 20 reduction-log rows,
and minimal-model data.  This is one independent point short of the target;
it is not a solution and no rank upper bound is claimed.

## Exact high-rank baselines

Two further replays keep the search calibrated against stronger public
families and records:

- Kihara's fully published family is reconstructed at (t=2).  Fourteen
  rational-function sections specialize to fourteen points with an exact
  full-rank certificate modulo 5.  This proves generic rank at least 14, but
  it supplies no rank-30 candidate and its simple specialization is not
  conductor competitive.
- The 29 displayed points on the Elkies--Klagsbrun record curve are checked
  exactly.  A full-rank matrix modulo 2 and an independent torsion witness
  prove rank at least 29.  This local replay does not claim exact rank 29 and
  does not produce a thirtieth point.

Their equations, sources, and certificate commands are recorded in
[baselines and literature](notes/BASELINES_AND_LITERATURE.md).

## Search drivers

`ecsearch/fermigier_score_sweep.cpp` is a deterministic staged
Mestre--Nagao-style ranking pass.  It precomputes local trace tables through
2000 and retains successively rarer rational parameters.  Its output is only
a heuristic ordering: every survivor still requires a bounded point search,
an exact independence certificate, global minimization, and an exact
conductor.  The current bounded sweep produced the rank-20 near miss above;
it did not produce a target curve.

## Important separation of local roles

If a prime \(p\) is forced into \(\Delta\), the specialized curve has bad
reduction at \(p\).  It therefore cannot simultaneously contribute an ordinary
good-reduction trace \(a_p\) to a Nagao score.  The search uses two roles:

- **shaping primes:** force prime powers in discriminant factors and score
  reduction type, local Euler factor, Tamagawa data, or root number;
- **scoring primes:** require good reduction and score finite-field traces.

The local constraints can still be assembled globally by CRT.  The final
specialization is always minimized before its conductor is accepted.

## Programme map

- [Mathematical pipeline](notes/CRT_LATTICE_PIPELINE.md) gives the exact
  congruence, lattice, primitive-pair, and conductor checks.
- [Baselines and literature](notes/BASELINES_AND_LITERATURE.md) records what is
  external, conditional, and current as of 2026-08-14, together with the exact
  Kihara rank-14 and public record rank-29 replays.
- [Fermigier reproduction](notes/FERMIGIER_REPRODUCTION.md) records the exact
  canonical-family bridge and unresolved factor-two source discrepancy.
- [`families/`](families/) contains the rank-at-least-two calibration family
  and the exact Fermigier--Mestre rank-at-least-twelve adapter.
- [`ecsearch/`](ecsearch/) contains dependency-free exact arithmetic.
- [`scripts/`](scripts/) contains generation and replay commands.
- [Reproduction commands](REPRODUCE.md) are the public entry points.

No target curve has been found in this repository yet.  The immediate search
priority is a twenty-first point on the rank-20 near miss or on a
conductor-qualified odd-root-number survivor, followed by broader rational
parameter sweeps and cofactor-aware constraint selection.
