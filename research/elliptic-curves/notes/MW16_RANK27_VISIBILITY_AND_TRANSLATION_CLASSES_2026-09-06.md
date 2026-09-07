# MW16 at -1867/270: exact visibility gain and distinct translation classes

The curve `new-20260906-90` retains certified rank at least **27**, with eleven
independent directions beyond this family's specialized generic sixteen.
This study explains the successful chart by exact arithmetic and proves a
useful consequence of its finite-reduction certificate. It does not establish
rank 28, exact rank, or a general theorem producing high-rank fibres.

The original equation, minimality, catalogue comparison and discovery history
remain in [the canonical rank-27 proof](NEW_MW16_RANK27_2026-09-06.md).
The later [specialized-parity experiment](SPECIALIZED_PARITY_CONTROL_2026-09-06.md)
and [2,3,5-saturation result](NEW27_SMALL_PRIME_SATURATION_2026-09-06.md)
also apply to this curve. None of those historical certificates is replaced.

## What made the discovered point visible

Let `H26` be the ordered 26-point subgroup immediately before the adaptive
follow-up, `G` its first sixteen generic generators, and `H27 = H26 + Z P`.
On the short search model the newly admitted point is

```
P = (-12447852056335152212800463/12,
     -20218538522161581965667523514701522728).
```

On the globally minimal integral model from the original certificate it is

```
P = (-1037321004694596017733372,
     -20218538522161063305165176216692656042).
```

The [visibility certificate](../../artifacts/generated-results/elliptic-curves/mw16_rank27_visibility_v1.json)
recomputes all centre combinations, exact rational chart transports, quartic
identities, coordinates and square witnesses. It also independently replays
the finite-reduction proof for all 27 points and their minimal-model transport.
Heights in the following table are primitive **search-coordinate heights**,
not canonical heights on the elliptic curve.

| Frozen roster | Best height for either P or -P | Chart, sign | Primitive coordinate | Signed presentations within height 100000 |
| --- | ---: | --- | --- | ---: |
| Original 43 generic-centre charts | 1379753695791601 | 22, P | [1379753695791601 : 1085900877990488] | 0 |
| First 301 adaptive charts | 79466 | 86, P | [-73571 : 79466] | 1 |

Thus these rosters differ by the exact factor
`1379753695791601 / 79466`, approximately **17.36 billion**, for this particular
signed point. Its adaptive height lies near the historical bound of 100000.
Merely extending the original charts to height one million would still not
expose either sign of this representative.

The winning centre has generic mask 585 and quotient word 672. The latter
uses quotient directions 6, 8 and 10 in one-based indexing; its combined
26-bit parity is 44040777. Its exact coefficient vector is

```
[-1,0,-2,1,0,0,1,0,0,-1,0,0,2,0,0,0,0,0,2,0,0,-1,2,-1,0,-1].
```

The original roster did **not** include generic mask 585, and the horizontal
maps also differ. This is a comparison of two completed rosters, not a
controlled causal ablation isolating quotient bits. The numerical height/CVP
policy chose the charts; it is not a certificate of lattice optimality.

## Translation-class lemma and certified corollary

For an elliptic curve `E/Q`, put `i_C(R) = C-R` and `tau_T(R) = R+T`.
Then

```
tau_T i_C tau_-T = i_(C+2T).
```

Consequently `i_C` and `i_D` are conjugate by translation by a rational point
if and only if `D-C` lies in `2E(Q)`. This follows by evaluating the displayed
identity; the converse uses a rational `T` with `2T = D-C`.

Suppose generators `B_1,...,B_r` of a subgroup `H` have independent columns
in a product of finite quotients `E(F_p)/2E(F_p)` at good primes. If
`sum a_i B_i = 2T` for `T in E(Q)`, reduction makes every `a_i` even. Hence
`H/2H -> E(Q)/2E(Q)` is injective. Distinct generator parities therefore
give distinct rational-translation classes of these involutions. This
argument needs neither a rank upper bound nor full saturation.

Apply this to the exactly certified `H27` and its first 26 and sixteen
generators. There are **1024 distinct lifts of each fixed generic parity**
inside `H26/2H26`, indexed by its ten quotient bits. The 43 original and 301
adaptive centres occupy **344 distinct classes in E(Q)/2E(Q)**. Every adaptive
centre has a nonzero quotient word, so none is rational-translation equivalent
to an original generic centre. The two rosters test only 344 of the
67,108,864 classes in `H26/2H26`; this is not a percentage of all search maps
or of all rational points. Different coordinates or rational translations
within one class can still change a finite box.

The [translation-class certificate](../../artifacts/generated-results/elliptic-curves/mw16_rank27_translation_classes_v1.json)
checks the injectivity witness and each parity/centre identity. This is an
elementary structural lemma with a new certified application to this example,
not a claim that the half-lattice method itself is new. For its established
high-rank context see [Elkies, Three lectures](https://arxiv.org/abs/0709.2908).

The invariant concerns pointed involutions under rational **translations**.
It does not classify abstract genus-one curves, prove an insoluble covering,
or supply a missing rational point. These pointed curves already have known
rational endpoints.

## Why the two quartic roots are one new direction

On the short model write `C=(x_C,y_C)` and

```
t_C(R) = (y_R + y_C)/(x_R - x_C).
```

The line through `-C` and `R` has third intersection `C-R`, by the chord law.
Thus `t_C(R) = t_C(C-R)`, with the usual tangent limits. Any invertible
horizontal coordinate change preserves this equality. For every chart and
both signs of `P`, the checker independently evaluates the rational point
`C-P` or `C+P` and verifies the identical primitive coordinate. Neither a
missed sign nor the other root of the same quadratic explains the initial
visibility gap.

Because `C` belongs to `H26`, the two partners represent opposite nonzero
classes modulo `H26 tensor Q`. They add precisely one direction together.
This is an exact quotient statement; it does not audit every translate or
combination representing that direction. A pointwise miss must not be
reported as absence of the entire quotient class.

## Bounded follow-up and research consequences

A separately frozen trial selects the first twelve charts in the already
fixed 49-chart numerical order on this curve's own 27-point subgroup. Each
receives height 1000000 and a 15-second cap, with one worker, at most 180
seconds of point invocations and a 300-second worker limit. It changes no
map or parameter and permits no refill. This tests a finite amount of
additional coordinate exposure, without extrapolating the known point's
height ratio to an unknown point.

The [trial certificate](../../artifacts/generated-results/elliptic-curves/mw16_rank27_visibility_followup_v2.json)
records completion versus censorship and the exact finite-reduction bounds
on the full seed-plus-output cloud. Version 1 stopped before any point
invocation because its guard compared exact rational values to their
serialized strings. Its source, protocol, failure log and supervision are
preserved. Version 2 corrects that type comparison, with identical
mathematical inputs and budgets.

All twelve invocations **timed out**; no box completed and no finite point
was returned. The retained cloud is exactly the original 27-point seed,
independently rechecked modulo 2, 3 and 5. This trial provides a cost/censorship
observation, not a completed-box null result or an improved rank bound.
No larger budget is launched automatically.

The research implications are specific:

- Rational parity classes provide a justified deduplication rule for centre
  involutions; a single assigned quotient word is not all lifts of a generic
  class. This says nothing about optimal representatives inside a class.
- Use exact post-discovery coordinates to distinguish a visibility failure
  from a failure to find additional rational directions. Keep such oracles
  out of prospective selection unless the experiment is explicitly retrospective.
- A future causal test should hold the generic mask and coordinate policy
  rule fixed while varying quotient words. This study identifies that missing
  control; it does not silently run or certify it.
- The known 2,3,5-saturation closes small-prime division as a repair of this
  particular subgroup. More independent directions and divisions at other
  primes remain open; bounded point searches provide no upper bound.

## Reproduction

The compact input retains exact geometry, centre coefficients, square hits,
point witnesses and the 27-point independence proof, with hashes identifying
the original transcripts. No new point enumeration is needed to replay the
visibility and translation results:

```sh
python3 elliptic-curves/cas/study_mw16_rank27_visibility.py check
python3 elliptic-curves/cas/certify_mw16_rank27_translation_classes.py --check
python3 elliptic-curves/cas/followup_mw16_rank27_visibility_v2.py check
python3 -m unittest discover -s elliptic-curves/tests -p test_mw16_rank27_visibility_study.py
```

The five targeted tests reject altered quartic coefficients, an off-curve
oracle and a wrong specialization, and distinguish a returned point from
its separately retained square witness.
