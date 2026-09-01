# Q80 third-q12 exact-lift checkpoint

## Status

The fixed `u=-2` construction has passed the horizontal lift gate.

- The exact `p=19` child remains fully certified: resolved genus-one pencil,
  Jacobian and maps in both directions, minimal `I6+I4+3I2+8I1` model, and
  transported `A5+A3+3A1` marking.
- The complete finite-field pipeline also passes at
  `p=61,67,83,89,103,131` with the same discrete profiles.
- Exact Hensel reconstruction proves that the characteristic-zero horizontal
  is defined over a biquadratic field, not one quadratic field.
- The connected resolved pencil nevertheless descends coefficientwise to the
  third quadratic subfield `QQ(a*b)`.
- The exact horizontal passes a direct test at the previously reserved good
  prime `p=71` for all four sign conjugates.
- All seven complete finite children are now uniquely aligned with that exact
  quadratic pencil at the generator and base-`PGL2` levels.

The complete connected correction, exact two-dimensional pencil, and generic
genus-one gate pass over the biquadratic field, with the pencil itself defined
over `QQ(a*b)`. The active gate is transport of the already computed finite
Jacobians and explicit maps into this common exact gauge, followed by a new
CRT/LLL lift. No characteristic-zero child equation or Mordell--Weil rank is
asserted at this checkpoint.

## Exact closure operands

The six characteristic-zero polynomial-closure equations are even in the
leading coordinate `l`. Replacing `l^2` by `q` leaves five rational unknowns.
At each of the two `p=19` operand branches, the resulting five-by-five
Jacobian is nonsingular. Newton--Hensel doubling and exact rational
reconstruction produce two literal solutions of the original six equations.

Write their reconstructed leading squares as `q1` and `q2`. Exact square
tests prove that `q1`, `q2`, and `q1/q2` are nonsquares over `QQ`. Hence the
two operands generate the degree-four field

```text
K = QQ(a,b),   a^2=q1,   b^2=q2.
```

The certificate includes exact substitution into all six closure equations
and literal reduction to the original `p=19` operand pair:

```text
artifacts/generated-results/
  q80-third-q12-um2-biquadratic-closure-operands-p19-hensel-qq.json
```

Replay it with:

```bash
sage elkies-k3/scripts/lift_q80_third_q12_closure_operands_p19_qq.sage \
  --biquadratic-operands
```

The largest reconstructed coefficient has 36,335 bits. These sizes explain
why a direct low-height rational search did not expose the lift.

## Exact horizontal and held-out prime

The two reconstructed operands can be composed exactly on the original
characteristic-zero Weierstrass model. In the basis `(1,a,b,a*b)`, the
result has

```text
x in span(1,a*b),   y in span(a,b).
```

Literal substitution proves the characteristic-zero Weierstrass identity.
This is a degree-four point in general, although its `x` coordinate descends
to the third quadratic subfield `QQ(a*b)`.

The same certificate reduces the exact object directly at the previously
reserved `p=71`. There both `q1` and `q2` are nonsquares while `q1*q2` is a
square. All four sign conjugates pass:

- `P.O=2`;
- height `8` by both height computations;
- the finite `I1*` identity component condition;
- the expected numerator and denominator degree profiles.

Modulo section sign, the four conjugates give two unsigned classes. The
certificate is:

```text
artifacts/generated-results/
  q80-third-q12-um2-biquadratic-horizontal-qq.json
```

Replay it with:

```bash
sage elkies-k3/scripts/certify_q80_third_q12_biquadratic_horizontal_qq.sage
```

This proves the exact horizontal and an independent held-out-prime
realization. The separate exact connected compiler below now proves its
resolved pencil, but not yet the child Jacobian.

## Why the earlier CRT route is retired

The exact square classes have the following local characters, listed as
`(q1,q2,q1*q2)`:

| prime | characters | inert operand |
| ---: | :---: | :---: |
| 19 | `(-,+,-)` | `q1` |
| 61 | `(+,-,-)` | `q2` |
| 67 | `(+,-,-)` | `q2` |
| 83 | `(-,+,-)` | `q1` |
| 89 | `(+,-,-)` | `q2` |
| 103 | `(-,+,-)` | `q1` |
| 131 | `(+,-,-)` | `q2` |

Thus the unique local quadratic target alternates between the two independent
global square classes. Local trace, norm, and coefficient discriminant are
generator-free at one prime, but they are not all reductions of a single
quadratic conjugation quotient across this mixed prime set. The old CRT
integers therefore cannot be used as rational reconstruction candidates.

The seven-prime artifact remains useful as a literal residue and branch
diagnostic. It is schema-versioned and status-labelled accordingly:

```text
artifacts/generated-results/
  q80-third-q12-um2-frobenius-crt-interface.json
```

```bash
python3 elkies-k3/scripts/compile_q80_third_q12_frobenius_crt_interface.py
```

It accumulates 1,947 ordered local slots modulo `7739891239523`, but its claim
boundary explicitly excludes rational reconstruction in its displayed local
generators and base gauges.

## Exact quadratic descent and repaired alignment

The exact 63-term moving equation reveals a smaller coefficient field than
the horizontal itself.  With `theta=a+b`, every coefficient has the form
`alpha*theta^2+beta` with `alpha,beta` rational.  Hence the connected pencil
descends exactly to

```text
F = QQ(theta^2) = QQ(a*b),
omega = 2*theta^2-2*(q1+q2) = 4*a*b,
omega^2 = 16*q1*q2.
```

This is the common global explanation of the alternating inert operand:
although the local extension is presented through `a` at some primes and
through `b` at the others, `q1*q2` is nonsquare at all seven primes.  The
local fields are therefore all reductions of `F`.

The new alignment certificate fixes more than the square class.  For each of
`p=19,61,67,83,89,103,131`, it finds the unique signed scale

```text
omega -> kappa_p*(2*r+ell_p)
```

and a three-by-three transformation on the quadratic `V` coefficients that
replays all 63 exact moving-equation coefficients.  The transformation
preserves the binary-quadratic discriminant conic up to a nonzero scalar, so
it is the symmetric-square action of a base `PGL2` change.  The certificate
is:

```text
artifacts/generated-results/
  q80-third-q12-um2-exact-quadratic-pencils-aligned.json
```

Replay it with:

```bash
python3 elkies-k3/scripts/align_q80_third_q12_exact_quadratic_pencil_primes.py --check
```

Its SHA-256 is
`d5237cb963f619a2cbee3d917c93fac1cdfd941a56fdb9341b43bb0974ce8331`.
This does not retroactively validate the old centered CRT integers: the seven
child models and maps must first be transported through the certified
generator scales and base gauges.

The legacy `p=19` resolved pencil uses a different local quadratic generator
and base gauge from the common-producer p=19 pencil.  A second run of the same
worker aligns that presentation explicitly; its artifact SHA-256 is
`e1ea9cc31518745ad569b0b8aebd23b74394b38ae33d22175674adf15581d429`.

The first transport stage is complete for all seven long Jacobians at
`p=19,61,67,83,89,103,131`.  Substituting the inverse certified base change
and rewriting every local coefficient in the basis `(1,omega)` gives one
common exact gauge.  The long degree profiles remain

```text
a1 2/2, a2 4/4, a3 4/4, a4 6/6, a6 8/8,
Delta 18/18, j 24/24.
```

At every prime, literal recomputation of `Delta` and `j` from the transported
Weierstrass coefficients matches the transported source invariants.  The
certificate is:

```text
artifacts/generated-results/
  q80-third-q12-long-jacobians-exact-quadratic-gauge.json
```

Its SHA-256 is
`38ca28d9fd643f2c2fd3b5aff706543e640f514a511ffc52023d29acb3675d07`.
Replay it with:

```bash
sage -python elkies-k3/scripts/transport_q80_third_q12_long_jacobians_exact_quadratic.sage --check
```

A fresh, valid CRT ledger now accumulates the resulting 292 rational-function
coefficient coordinates modulo `7739891239523`:

```text
artifacts/generated-results/
  q80-third-q12-long-jacobian-exact-quadratic-crt.json
```

Its SHA-256 is
`11343a3a9f9432987208e6b0b8238657564818e438a6c98ef9336d610682492f`.
All seven-prime reductions replay.  Apart from forced zero and monic slots, the
centered residues are already of generic modulus size, so the present modulus
does not support rational reconstruction; more aligned primes are required.

Both directions of the generic birational maps are now transported in the
same way at all seven primes.  Forward blocks are stored as rational
functions in exact `V` and old `W`; inverse coefficients retain the pinned
weighted monomials in the long Weierstrass `X,Y`.  The common transported
support and degree schema passes at every prime:

```text
artifacts/generated-results/
  q80-third-q12-birational-maps-exact-quadratic-gauge.json
```

Its SHA-256 is
`3e9f38ed42976fe3a554ba9ddac27e2d86c8505d944b44f03ac9844ca00447ef`.
The corresponding 3,484-slot CRT ledger, again modulo `7739891239523`, is:

```text
artifacts/generated-results/
  q80-third-q12-birational-maps-exact-quadratic-crt.json
```

Its SHA-256 is
`ca24c48f208fb39fa9817e39f4ab3eb7d1e9e44d6987608606ce6c57def1c114`.
Replay both with:

```bash
sage -python elkies-k3/scripts/transport_q80_third_q12_maps_exact_quadratic.sage --check
python3 elkies-k3/scripts/compile_q80_third_q12_maps_quadratic_crt.py --check
```

## Active exact-child path

The immutable `p=19` connected compiler uses only field arithmetic after it
loads the horizontal. The exact adapter now certifies the same sequence over
`K=QQ(a,b)`:

1. Smith saturation with degrees `(0,0,6)`;
2. the seven-dimensional shifted-Popov ambient;
3. the complete connected `D7` ideal;
4. the finite connected `D5` quotient;
5. the rank-five gate and two-dimensional kernel;
6. the moving equation of degrees `(2,9,3)`.

Run it with:

```bash
sage -python elkies-k3/scripts/compile_q80_third_q12_biquadratic_resolved_pencil_qq.py
```

The resulting 117 MB artifact has SHA-256
`ac67210166cd414945e1fa373e8f0d5829ff8231daf83c764c376a32ff4b641e` and
contains 63 exact nonzero moving-equation terms:

```text
artifacts/generated-results/
  q80-third-q12-um2-biquadratic-resolved-pencil-qq.json
```

Irreducibility of one good `p=19` reduction of the 63-term exact moving
equation proves characteristic-zero irreducibility. Together with the pinned
primitive isotropic divisor, complete pencil, old-fibre degree three,
Bertini, and K3 adjunction, this certifies a smooth generic genus-one member:

```bash
sage -python elkies-k3/scripts/verify_q80_third_q12_biquadratic_resolved_genus_qq.sage
```

The remaining ordered gates are:

1. compute the exact Jacobian and retain birational maps in both directions;
2. minimize the exact Jacobian;
3. factor its discriminant and certify `I6+I4+3I2+8I1`;
4. transport the old components and zero to certify `A5+A3+3A1`;
5. use finite fields only as independent replays, not to infer the
   characteristic-zero Mordell--Weil rank.

The failed first-marking genus-two-cover hypothesis remains closed: its
quadratic field splits at `p=19`, so it cannot supply the mandatory local
quadratic control there.

## Recommended Jacobian lift

The long Jacobian and two-way map transport stages, with valid new CRT
ledgers, are now done for all seven primes including the mandatory legacy
`p=19` control.  The next computational task is to add more complete primes
in the same exact gauge; the present modulus is tiny compared with the known
36,335-bit horizontal heights.
Minimal models should be derived after reconstructing the long model, so the
pinned simple-branch Laurent gauge is not obscured by prime-local minimal
scalings.

The alternative remains a coupled Hensel lift of the generic map identities,
seeded at `p=19`.  In either route, accept an exact reconstruction only after
literal substitution into the exact 63-term pencil over `QQ(omega)` and keep
`p=71` as an independent replay.  This avoids an expensive unspecialized
normalization over the original 36,335-bit biquadratic presentation.

## P-adic conductor path

The p-adic alternative has now passed its first two exact gates.  The exact
63-term pencil is compiled coefficientwise modulo `19^64` in the global basis
`(1,omega)`; every rational denominator is a 19-adic unit and the quadratic
extension is unramified:

```text
artifacts/generated-results/
  q80-third-q12-exact-pencil-p19-adic-precision64.json
```

SHA-256:
`fc8b86cd750a906804216d17fc348ab8c6bf7704b9c2036d78026dfb2e4dc945`.

The cubic discriminant has mod-19 factor shape `L^3 Q^2 D`, with degrees
`(1,4,4)`.  Treating the nine non-leading coefficients of `L,Q,D` as the
unknowns gives a rank-nine Hensel Jacobian.  A fixed-inverse digit lift proves
the factor identity through five p-adic digits, with valuation history
`1,2,3,4,5`.  The repeated cubic root modulo `LQ` lifts with the same history,
and produces the unique normalized candidate

```text
e=(z^2+A*z+B)/(L*Q),  deg_W(A,B)=(4,4),  B=0 mod L.
```

The certificate is:

```text
artifacts/generated-results/
  q80-third-q12-discriminant-factors-p19-adic-precision5.json
```

SHA-256:
`b98cb5dae2aa9dbfe05f02fd5e6bed7b23b17bd383b23156c17c1b383b5a705d`.

The proposed basis has now passed the generic modular integrality test.  A
dedicated exact worker represents every rational function as a numerator over
a power of one fixed global denominator `H(U)`.  In the multiplication matrix
of `z^2+A*z+B` on `(1,z,z^2)`, it verifies coefficientwise that the trace,
second symmetric coefficient, and determinant are divisible by `LQ`,
`(LQ)^2`, and `(LQ)^3`, respectively.  Thus `(1,z,e)` is an integral basis of
the generic cubic modulo `19^5` over the stated localization.  The certificate
is:

```text
artifacts/generated-results/
  q80-third-q12-integral-basis-mod19-power.json
```

SHA-256:
`18b2bc2dbb538a06c8f911f900f5966fa28a38ab5915fed620c2d38887e4d222`.
Replay with:

```bash
python3 elkies-k3/scripts/compile_q80_third_q12_exact_pencil_p19_adic.py --check
sage -python elkies-k3/scripts/lift_q80_third_q12_discriminant_factors_p19_adic.sage --check
python3 elkies-k3/scripts/verify_q80_third_q12_integral_basis_mod19_power.py --check
sage -python elkies-k3/scripts/compile_q80_third_q12_riemann_roch_p19_adic_sample.sage --check
```

This is a five-digit modular certificate, not characteristic-zero integrality.
The certified basis has also been used for the first p-adic Riemann--Roch
sample.  The exact-gauge value `U=16+7*omega` reduces to the pinned legacy
`T=1` sample.  Sage cannot directly construct the quadratic extension above
the unramified quadratic p-adic field, so the worker takes the norm quartic
over `QQ_19`, enforces regularity at both conjugate double branches, and then
checks that the normalized generators descend to `QQ_19(omega)`.  It obtains

```text
dim L(0P),...,dim L(3P) = 1,1,2,3,
ord_P(x),ord_P(y) = -2,-3
```

through all five digits.  The sample certificate is:

```text
artifacts/generated-results/
  q80-third-q12-riemann-roch-p19-adic-sample.json
```

SHA-256:
`9ec5c3558a461a9d1342a4c293dccf6949700f4ce1f05c241524c6dad5ed80fe`.

The sample now also contains the long Weierstrass relation and maps in both
directions.  The inverse weighted bounds are exactly the finite-field bounds
`4` for `W` and `10` for `z`; both identities are replayed coefficientwise
using the common conductor denominators `X/(LQ)^2` and `Y/(LQ)^3`.  Reduction
of the five long coefficients under `omega -> 2+13*r` literally reproduces
the pinned legacy `T=1` positive control.

The worker no longer hardcodes the legacy infinity roots.  It derives the
unique simple and double roots from the weighted infinity cubic modulo 19,
and has produced 95 accepted residue-distinct samples in the exact global
gauge.  Interpolation of the rational long coefficients gives degrees
`2/2,4/4,4/4,6/6,8/8` modulo `19^5`, using 93 samples and reserving two as
held-outs.  Every coefficient also reduces literally to the independently
transported generic p=19 long model.  The generic certificate is:

```text
artifacts/generated-results/
  q80-third-q12-long-jacobian-p19-adic-precision5.json
```

SHA-256:
`a2e65c630a9d79bf59f2d45ce9bac1ea69dc4d52bd8bdadaf45295f407cd9ba8`.

Replay with:

```bash
python3 elkies-k3/scripts/interpolate_q80_third_q12_long_jacobian_p19_adic.py --check
```

The generic maps now pass as well.  Direct scalar interpolation initially
exposed uncancelled conductor factors in the p-adic rational functions; a
literal function-field comparison showed that the pinned maps nevertheless
equal both positive controls.  The interpolator now canonicalizes every
sample by solving `n*D-d*N=0` at the transported `W`-degree bounds.  Of the 95
samples, 93 have the canonical generic degrees; 91 train the maps and two are
held out.  This covers the complete forward map, inverse `W` map of degree
`6/6`, and inverse `z` map of degree `40/40`.  Every scalar rational function
passes both held-outs and reduces literally to the independently transported
generic p=19 map.

Thus the complete long child and maps in both directions are now available
over `(Z/19^5)[omega](U)`.  The next gate is higher precision for rational
reconstruction, followed by exact minimization, fibres, and marking.  No
characteristic-zero Jacobian is claimed from five digits.

Intermediate rational compression after every Hensel digit removes the prior
expression-growth wall: the conductor factors and repeated root now lift with
valuation histories `1,...,16` in under twenty seconds.  Generic integrality
and the complete pinned child with maps both ways have independently replayed
through `19^16`.  The precision-16 artifacts are:

```text
q80-third-q12-discriminant-factors-p19-adic-precision16.json
q80-third-q12-integral-basis-mod19-power-precision16.json
q80-third-q12-riemann-roch-p19-adic-sample-precision16.json
```

Their SHA-256 hashes are, respectively,
`5d1d1999bc1d53c7fd20b3a5f5f904fb248a332b0871f712e69b580d54a8273f`,
`852fb6503c10696d38b41443137b9bcf392309f51cbe89064b48ecc4b61572ca`, and
`7a3ba5d0d4d8bb97e140b63fcd4c637183373e7ab06c82768a8d7ff075122098`.
The exact pencil was subsequently recompiled modulo `19^260`.  The complete
factor/root lift, generic integral-basis test, and pinned child with maps both
ways now pass through `19^256`.  The source, factor, basis, and child hashes
are, respectively:

```text
1ef2f30c59990caa8e71f3612ae563ac3de34a7b4128532cb26d7c6a824a0410
e13b1a0f7fd980bc121838ea97053724ba5f26ff9f8775c24e5df7d5bd533403
0183bfc64128e6ab3e686ec1b12cfdf26172041ceba1d9979b1338fc90077d41
c2986e59579b6a34c660bad3ca2be39866377dbf2e293089350d9156e0a12af8
```

Rational reconstruction of the pinned long coefficients is still
insufficient: most coordinates fail, and every returned fraction lies at the
543--544-bit boundary of the 1,088-bit modulus.  No such fraction is accepted
as exact.  The next target should be at least `19^1024`, preferably with
checkpointed factor lifting, before repeating generic interpolation.

Two prospective shortcuts were tested and rejected.  Reloading a completed
lower-precision JSON factor artifact into a larger p-adic field does not
preserve the normalized fraction-field representatives needed by the Hensel
system, while solving the full current 9-by-9 p-adic Jacobian is already too
slow at low precision.  Neither mode is retained in the certified worker; its
hash again agrees with the precision-5 and precision-256 artifacts.  A useful
checkpoint implementation must instead save native Sage state at the final
target precision.  The tempting first-marking genus-two shortcut is not an
alternative: the exact local-behaviour certificate already proves that field
splits at 19, while the target field is inert.  The exact fixed-`u` pencil is
already available over the different quadratic descent field `QQ(a*b)`; the
remaining issue is reconstructing its Jacobian and maps, not identifying its
coefficient field.
