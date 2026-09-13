# Effective coordinate-height bounds on an actual V3 chart

Four models of the first pointed chart in the retained Curve302 M27 benchmark
now have independently replayed height bounds. These comprise the existing
factor-free and full-minimization outputs and two newly constructed minimal
models. A pair of minimal models has a better joint bound than either alone.
This is a concrete coordinate-coverage result. No point search ran, no new rank
direction was found, and no CPU advantage or cross-curve transfer is established.

The [portable packets](../../artifacts/generated-results/elliptic-curves/pointed_height_bounds_first_chart_v1/protocol.json)
contain the actual equations, complete rational parameter matrices, real bounds,
finite residue trees and independent receipts. The
[joint certificate](../../artifacts/generated-results/elliptic-curves/pointed_height_bounds_first_chart_v1/portfolio-verified.json)
records the local coupling and every subset of the three minimal models.

## Scope and inputs

Select the first centre in the first development case, `curve302-M27`, of the
[retained next-direction benchmark](NEXT_DIRECTION_RANK32_BENCHMARK_2026-09-12.md).
Reuse its two existing normalization maps. Model selection uses that equation,
anchor and certified discriminant-prime data; target points and acquired-direction
outcomes are not inputs to the bound or neighbour algorithms. This is retrospective
commissioning on a known control, not an independently held-out validation.

The earlier [anatomy note](CURVE302_V3_ANATOMY_AND_FINITE_ATLAS_2026-09-11.md)
already proves the pointed-chart identities and identifies the uncomputed local
constants. Its missing external measurement payload remains missing. This work
uses actual retained V3 normalization matrices and does not promote those external
measurements. The two blindly recovered strict covers motivate representation
choice but are not inputs to this calculation; their constructor is not reopened.

Fisher–Sills supplies the relevant framework: explicit local contributions to
covering-height bounds, and examples where searching different minimal models
helps. Sections 4–5 also explain why an elimination bound alone can be poor.
Our implementation below uses direct rational-polynomial bounds on the pointed
map; it does not claim to implement their full Tamagawa-distance algorithm.
[Fisher–Sills, 2012](https://arxiv.org/pdf/1103.4944).

## The checked height inequality

Fix the actual short Weierstrass equation and anchor `Q=(a,b)`. For a point `P`,
put `R=2P-Q`. Compose the exact identity `x(R)=N_Q(t)/F_Q(t)` from the anatomy
note with the retained rational parameter matrix. Clear denominators and remove
joint content to obtain integral homogeneous quartics `N(U,V),D(U,V)`.
The transformed discriminant quartic is denoted `q(U,V)`.

For primitive integers `(m,n)`, set `M=max(|m|,|n|)` and
`g=gcd(N(m,n),D(m,n))`. The identity is

\[
h(t)-\tfrac14h_x(R)
=-\tfrac14\log\!\left(\frac{\max(|N(m,n)|,|D(m,n)|)}{M^4}\right)
+\tfrac14\log g.
\]

Suppose the packet proves

\[
0<L\le\max(|N(u,v)|,|D(u,v)|)\le U
\quad(\max(|u|,|v|)=1,\ q(u,v)\ge0),\qquad g\mid C.
\]

Then valid constants are

\[
B_1=-\tfrac14\log U,\qquad B_2=\tfrac14\log(C/L).
\]

In multiplicative heights this gives the exact coverage implication

\[
H_x(2P-Q)\le H^4L/C\quad\Longrightarrow\quad H(t(P))\le H.
\]

This compares models of the **same marked chart** and the same elliptic equation.
It does not compare unrelated anchors by replacing their different `2P-Q` with a
common unknown quantity. It is a sufficient coordinate bound, not an assertion
that every visible point lies inside the guaranteed elliptic-height range.
It also does not guarantee that the point is independent of the starting subgroup.

Both projective parameter charts are included. At `D=0`, use projective
`[N:D]`, with `H_x(O)=1`. The two Bezout identities exclude simultaneous
vanishing, including at infinity. Execution must still enumerate the relevant
box and exceptional parameters completely; a timed-out backend is not complete
coverage merely because the coordinates satisfy a proved bound.

## Computing and replaying the constants

The producer uses exact real-algebraic arithmetic on `(t,1)` and `(1,t)` for
`-1<=t<=1`, restricted to `q>=0`. Extrema occur at endpoints, branch points,
crossings `N=+/-D`, or stationary points of the active polynomial. The resulting
extrema are rounded outward to rational `L,U`. The independent checker uses
Fraction arithmetic and simultaneous rational Sturm root isolation to check
every sign cell of the required inequalities. It is not a grid sample.

Two rational polynomial Bezout identities, in opposite parameter charts,
give an integer `C0` with `g|C0`. Specifically, after homogenizing their
degree-at-most-seven identities, a prime sees a unit in at least one of `m,n`.
Thus `g` divides the lcm of their two constants. This baseline needs no factoring.

At a previously certified prime dividing `C0`, a finite tree covers `Z_p` and
the projective infinity chart. It strips common coefficient valuations and
refines every simultaneous residue root of `N,D`. A whole ball is excluded only
when `q` has a fixed odd valuation or a fixed nonsquare unit. The test at2 uses
the unit modulo 8. The independent checker proves residue-list completeness by
polynomial gcd with `X^p-X`, rather than trusting a supplied root list.

Unresolved branches retain their original Bezout cap; unprocessed prime support
retains its unfactored integer contribution. Therefore an incomplete refinement
can weaken the bound but cannot silently strengthen it. These are valid upper
bounds on cancellation, not proofs that every retained local case is attainable.
All 20 contributing primes of the initial pair were covered, using the retained
unconditional prime certificate. No new large-integer factorization was needed.

## Actual models and results

The initial plain Bezout bounds gave `B2` approximately 74.03 and 71.59, which
give no useful guaranteed range at height 125000. Local refinement changes this:

| Representation | Approximate `B2` | Guaranteed upper `h_x(R)` at `H=125000` |
|---|---:|---:|
| Existing factor-free output | -13.83777 | 102.29535 |
| Existing full-minimization output | -15.38486 | 108.48371 |
| New minimal neighbour at 5, residue 1 | -14.81659 | 106.21062 |
| New minimal neighbour at 7, residue 0 | -14.80417 | 106.16096 |

Logs are natural. The exact inequalities use the rational `L,U,C` in the
packets, not rounded table entries. For the initial pair the cancellation
divisor drops by exactly `6^4`; the real contribution partially offsets that
improvement. The resulting sufficient parameter-height bound improves by
approximately 4.69779. This is not a measured point-search speedup.

The existing factor-free output is not minimal: its `c4` is `12^4` times that
of the minimal elliptic model. The full-minimization output has exactly the
minimal elliptic `c4,c6`. The new models are constructed in fixed order from
integral local changes

\[
x=r+px',\qquad y=s+py'.
\]

If the old equation is `y^2+Q(x)y=P(x)`, the new coefficients are

\[
Q'=(Q(r+px')+2s)/p,\qquad
P'=(P(r+px')-sQ(r+px')-s^2)/p^2.
\]

Require integrality and unchanged invariants, then apply `hyperellred`.
The first two admitted models occur at `(p,r)=(5,1),(7,0)`. Fraction-only
replay checks their equations, parameter maps, square transports and invariant
identities. All three minimal models have the discriminant of the verified
minimal elliptic equation. Their recorded height-box keys differ; a complete
classification of minimal models or global `GL2(Z)` orbits is not claimed.

## Why the pair helps

Let `a_p` indicate the residue ball `u=r v (mod p)` in the base minimal model.
The exact parameter matrices give, on primitive parameters,

\[
k_{\text{neighbour},p}=k_{\text{base},p}+2-4a_p.
\]

Outside the neighbour prime the local cancellations agree. The replayed base
caps are `(0,2)` outside/inside the selected ball at 5 and `(2,4)` at 7.
These produce four joint local cells. For each cell `s` and model `i`, let
`D_i(s)` be the rational bound satisfying
`H(t_i(P))^4 <= D_i(s) H_x(2P-Q)`.

For a set `T` of models the correct bound is

\[
D_T=\max_s\min_{i\in T}D_i(s).
\]

The model is chosen once for each global point. Summing independently minimized
contributions at different primes would not justify this conclusion.

The base minimal model together with the neighbour at 7 gives
`B2=(log D_T)/4` approximately **-15.77713**. Adding the neighbour at 5 does
not improve this worst-case bound. The selected pair reduces the sufficient
height per box by about **1.48033** relative to the best single model.

Two boxes must be charged. Under the deliberately limited rational-address
work proxy `number_of_models * H^2`, their work ratio is **0.91266** relative
to that single model. The exact subset comparison uses
`number_of_models^2 * D_T`, avoiding floating ranking. This roughly 9% address
advantage can disappear after model preparation, sieving and verification costs.
It is neither a CPU prediction nor evidence of an additional independent point.

## Cost and the next gate

The initial two bound computations took 0.325 CPU seconds in their arithmetic
kernels; their separate checkers took 0.452 seconds. Four fresh processes used
2.615 metered CPU seconds and 3.523 elapsed seconds. These component timings
exclude the historical generation of the retained maps. The neighbour factory's
successful kernel took 0.014 seconds; further bound/replay timings are retained
in `meter.json` and the individual receipts. Development included a constant-
polynomial root interface failure and two neighbour numeric-conversion failures.
Those diagnostics survive locally. Full development/commissioning CPU was not
completely metered, so no whole-pilot speed comparison is claimed.

The next gate is a frozen next-direction comparison using the existing control
banks and the same starting subgroups, with all new model construction, bounds,
failed attempts, search and verification charged. The coverage-work proxy must
earn its place by reducing CPU to a certified new direction. The comparison
must distinguish retained-bank timings from full cold setup;
unknown historical landscape costs remain unknown, rather than being charged
as zero. A later cross-curve test needs withheld directions sourced independently of the V3
clouds used to build these controls. Curve302's alternate subgroup is useful
validation on that curve, not a broad transfer test. That benchmark has not
been launched by this calculation; no parameter campaign is opened.

## Replay

From the repository root, using an output path that does not already exist:

```sh
sage -python research/elliptic-curves/cas/verify_pointed_height_bounds.py \
  --packet research/artifacts/generated-results/elliptic-curves/pointed_height_bounds_first_chart_v1/factor_free-bounds.json \
  --output /tmp/pointed-height-replayed.json
sage -python research/elliptic-curves/cas/verify_pointed_height_portfolio.py \
  --folder research/artifacts/generated-results/elliptic-curves/pointed_height_bounds_first_chart_v1 \
  --output /tmp/pointed-portfolio-replayed.json
sage -python research/elliptic-curves/tests/test_pointed_height_bounds.py -q
```

The portfolio checker replays the three minimal-model bound packets and both
neighbour transports before checking the joint bound. Primality is reused from
the bound retained certificate. Sage algebraic-number production and independent
Fraction/SymPy verification are distinct implementations, not formal verification.
