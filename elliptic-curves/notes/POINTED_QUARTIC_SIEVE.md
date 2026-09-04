# Exact pointed-quartic slope sieve

Implementation: [`half_lattice_pointed_sieve.py`](../cas/half_lattice_pointed_sieve.py)
and [`pointed_quartic_sieve.cpp`](../cas/pointed_quartic_sieve.cpp).
The current experiment and rank boundaries remain in the
[canonical MW16 ladder note](ICARM_MW16_BLIND_LADDER_AND_PROSPECTIVE_GATE_2026-09-04.md).

## Denominator transform

Start with an integral short curve `E: y²=x³+A*x+B` and
`Q=(a/d²,b/d³)`, with `d>0` and `gcd(a,d)=1`. The pointed chart is

```text
w² = f(t) = t⁴ - 6*x(Q)*t² - 8*y(Q)*t - 3*x(Q)² - 4*A.
```

Compute the centered residue `k = -b/a (mod d²)` by an extended-gcd
inverse, without factoring `d`. For `d=1`, take `k=0`. Substitute
`t=d*s+k/d`, `w=d*v`. This gives `v²=g(s)` with ascending coefficients

```text
[(k⁴-6*a*k²-8*b*k-3*a²-4*A*d⁴)/d⁶,
 4*(k³-3*a*k-2*b)/d⁴,
 6*(k²-a)/d²,
 4*k,
 d²].
```

Every coefficient is integral, and the implementation checks every division.
To see this without factoring, work in `Z[1/a]`, put `h=-b/a`, and use
`C=b²-a³=A*a*d⁴+B*d⁶`. Then

```text
h²-a = C/a²,
h³-3*a*h-2*b = -b*C/a³,
F(h) = 4*B*d⁶/a + C²/a⁴,
```

where `F(z)=z⁴-6*a*z²-8*b*z-3*a²-4*A*d⁴`. Taylor expansion at `h`,
with `k-h` divisible by `d²`, proves the respective divisibilities by
`d²`, `d⁴`, and `d⁶`. Since `a` is coprime to `d`, they hold in `Z`.
The resulting binary quartic has `I=-48*A` and `J=-1728*B` exactly.

For rational short input, a preliminary exact scaling
`(x,y) -> (u²*x,u³*y)` makes the coefficients integral. Only fixed small
prime valuations and exact integer-root tests are used. If a remaining
denominator is not a perfect power, the implementation uses it as a safe
clearing scale. Likewise, it removes a common twelfth-power invariant factor
only when verified. These choices do not assert minimality and can leave
large coefficients.

## Horizontal lattice

The numerator lattice of `t` has basis `(d²,0),(k,1)`. Exact Gauss reduction
uses the positive metric

```text
N² + (|a| + d²*(floor(sqrt(|A|))+1))*D².
```

The recorded unimodular matrix `U=(alpha,beta;gamma,delta)` changes
`s=(alpha*z+beta)/(gamma*z+delta)`. The implementation expands the binary
quartic exactly and checks its invariants. The independent replay checks all
five coefficients against direct substitution into the original rational
quartic, as well as the MW16 combination producing `Q`.

This metric balances slope numerators and denominators. It does not promise
optimal height compression or the same bounded box as PARI reduction. The
historical adaptive 54/55 recovery is not automatically inherited by this
backend. The separate initial-wave control ledger measures actual recovery.

## Modular search and certificates

The worker enumerates `-H<=n<=H`, `1<=d<=H`, with `gcd(n,d)=1`.
For each useful odd prime through 151 it precomputes the condition that
`F(n,d)` is a square modulo that prime, including zero. Homogeneous evaluation
handles `p|d` directly; there is no inversion of `d` modulo a sieve prime.
Uninformative all-pass primes are omitted. These are square-value filters,
not local-solubility or Selmer tests.

Eight filters operate on 64 candidates at once. Remaining primes filter the
survivors before any large-coefficient polynomial evaluation. GMP then tests
the surviving integer values for exact squares. Python independently checks
each square, transports both ordinate signs, and verifies the resulting
point on the original curve. The transformed point at slope infinity is also
checked. Points over original slope infinity are the already-known `O,Q`.

The worker emits completed-denominator and sieve-stage counts. A timeout is
never a completed box: the internal deadline stops between denominator rows,
and an external timeout discards the interrupted chunk. `search_box` accepts
denominator intervals for explicitly partitioned larger boxes. The campaign
checkpoints every chart atomically and replays incomplete charts on restart;
candidate checkpoints separately retain exact group classification. Complete
chart checkpoints are keyed by inputs, bounds, and source hashes.

Compile prerequisites are a C++17 compiler and GMP development headers. The
Python wrapper builds the worker into an ignored directory keyed by its source
hash. Sage/PARI are used by the existing lattice and group-classification
runner; the new quartic backend itself uses neither Sage nor PARI.
Construction of `Q` uses exact GMP rational group arithmetic with simultaneous
signed binary multiplication and a separate 60-second safety limit. This
avoids the enormous intermediate rational multiples in a serial Python sum.
Input sections and the resulting point are checked on the original curve;
regressions compare against an independent Python rational group law. The
chart wall timer includes group arithmetic and transformations, while the
two-second search budget applies to the modular worker. Compilation is cached.

## Reproduction

```bash
python3 -m unittest discover -s elliptic-curves/tests -p 'test_half_lattice_pointed_sieve.py'

sage -python elliptic-curves/cas/run_icarm_mw16_nagao_finalist_half_lattice.sage \
  --backend pointed-sieve --height-bound 10000 --timeout-seconds 2 \
  --output artifacts/local/elliptic-curves/mw16-pointed-sieve-h10000-all.json

sage -python elliptic-curves/cas/calibrate_icarm_mw16_pointed_sieve.sage

python3 elliptic-curves/cas/verify_icarm_mw16_pointed_sieve.py \
  --raw artifacts/local/elliptic-curves/mw16-pointed-sieve-h10000-all.json
python3 elliptic-curves/cas/verify_icarm_mw16_pointed_sieve.py --check --replay-charts
```

`--candidate-start` and `--maximum-candidates` bound independent campaign
slices. Re-running the same output resumes it; changed inputs or budgets
require a new output path. The deterministic gzip certificate retains the
full search evidence so checking does not depend on ignored local files.
No new parameter sweep, adaptive lift, unrestricted search, or Selmer
calculation is part of this replay.

<!-- status-consumer: EC-K3-ICARM-MW16-POINTED-SIEVE cb83c1afae1d0141 -->
