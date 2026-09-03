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
checkpoint implementation must instead save native Sage state, or canonical
coefficient residues, at the final target precision.  The tempting
first-marking genus-two shortcut is not an
alternative: the exact local-behaviour certificate already proves that field
splits at 19, while the target field is inert.  The exact fixed-`u` pencil is
already available over the different quadratic descent field `QQ(a*b)`; the
remaining issue is reconstructing its Jacobian and maps, not identifying its
coefficient field.

## Precision 1024 gate

The exact 63-term pencil is now compiled modulo `19^1028`, and the conductor
factorization and repeated root pass through `19^1024`.  A restart test at
target precision 10 stopped once in the factor loop and once in the repeated-
root loop; the resumed output agrees coefficientwise with a fresh run.  For
the full lift, fixed-inverse digit corrections reached valuation 48.  Full-
precision 9-by-9 Newton solves at good `U` values, followed by interpolation
on the already certified rational-function supports and three held-outs,
then advanced

```text
48, 96, 192, 384, 768, 1024.
```

The repeated root similarly advances

```text
1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024.
```

At precision 16 the accelerated `L,Q,D,A,B` records are byte-for-byte equal
to the prior fixed-digit artifact.  The actual hybrid execution, final
checkpoint, canonical worker, and coefficientwise comparison of checkpoint
`L,Q,D` with the final artifact are separately pinned by

```text
q80-third-q12-discriminant-factors-p19-adic-precision1024-execution.json
```

with SHA-256
`47b02dc8483f068e82bcd5b32689a1c5d30fd8affc0231e4b0266f6d9b46485b`.
The source and final factor hashes are

```text
e0327dd535c4571a36d8b720e453090f7893087ae0c0c7aa1c65ed71d5964cbf
81757b810266bc6000959562e02dd80a276a1476d26536af1be83940e1ead961
```

Generic integral-basis divisibility also passes through `19^1024`; its hash
is
`c10bce1e647179f207c7cf291cadf9aa392f7ef27c9e1c69773d87b8f4ff06ea`.
The legacy-aligned sample independently recovers dimensions `1,1,2,3`, pole
orders `2,3`, the long Weierstrass equation, and maps in both directions at
the same precision.  Its hash is
`fff901bba6a6b8868edc9af724c12e407882e47c00c422f0d546354c211c0098`.

Coefficientwise rational reconstruction at this precision is still not a
certificate: successful fractions occur at the 2,175-bit square-root
boundary and some slots fail.  The next reconstruction therefore uses the
known rational-function supports, residue-distinct high-precision samples,
projective LLL per long coefficient, and independent reduction to all seven
aligned finite-prime models.  Literal characteristic-zero substitution into
the exact pencil remains mandatory before promoting the result from a
candidate equation to an exact Jacobian.

### Invariant-first reconstruction order

The immediate reconstruction gate is now the fixed-gauge degree-`24/24`
`j`-map, not the 3,484 coefficients of the two birational maps.  The cheaper
route needs only the existing 17 support-determining long-coefficient samples
and three held-outs: it interpolates the five small rational functions modulo
`19^1024`, then composes `c4`, `Delta`, and `j` symbolically.  Its numerator is
a scalar times the cube of a monic degree-8 polynomial, while the denominator
has squarefree degrees `8,3,1,1` at multiplicities `1,2,4,6`.  Reconstructing
those intrinsic blocks replaces an unstructured dimension-99 lattice by
dimension-17 and dimension-29 projective problems.

The resulting order is

```text
j -> (c4^3, Delta) -> minimal Jacobian -> birational maps.
```

The monic normalization of the `j` denominator pins the rational-function
representative but is not by itself a Weierstrass scaling.  The subsequent
`(c4^3,Delta)` step must therefore record its common scaling separately.
Neither modular reconstruction nor the exact numerator-cube and denominator-
factor identities alone prove literal characteristic-zero agreement with the
63-term pencil.

The sample and reconstruction commands are

```bash
sage -python elkies-k3/scripts/reconstruct_q80_third_q12_j_map_p19_adic.sage
```

The resulting candidate uses `19^1024` together with all six aligned
auxiliary primes `61,67,83,89,103,131`.  Its dimension-17 and dimension-29
primitive vectors have respectively 4,130 and 4,237 bits against random
boundaries 4,131 and 4,238.  All auxiliary primes were therefore consumed by
CRT; a new eighth aligned good prime is mandatory before promotion.  The
candidate artifact remains explicitly non-theorem status until that replay
or literal exact-pencil substitution succeeds.

### Blind eighth-prime rejection

The required clean replay was carried out at the next selected inert good
prime `p=163`.  The common producer was run before the characteristic-zero
candidate was reduced: the third-q12 horizontal, two-dimensional resolved
pencil, genus-one test, 72 finite Jacobian samples, generic interpolation,
exact quadratic-field alignment, and legacy base-gauge transport all pass.
The final independent inputs are pinned by

```text
resolved pencil: 11cb2c5493d814023e2697ded83db5fec15d34a02a11601ae0e053fad569ab67
interpolated Jacobian: 2d86b70cbad4537d1121f0308ac5faf5cbc8ea8a1635c504ee93fcdfcea73131
eight-prime alignment: 90f409015757a6ebca8214ca16cad28caa6879c1cdbc6b602db171776ba01751
eight-prime transport: f5acf0aeca875e4717bad8fc76b5b20eb00f299f9f8063b7c602c7c0df1e8c48
```

Only after those files existed was the candidate reduced modulo `163`.  It
fails literally in the pinned exact gauge: all 25 numerator coefficient
pairs and 24 of 25 denominator coefficient pairs differ.  The rejection
certificate is
`q80-third-q12-j-map-p163-heldout-replay.json`, SHA-256
`1da5975fc1f546b9e213fe2a44272e115c9442e78ba8030bf05336b866fe4e60`.
Thus `q80-third-q12-j-map-p19-adic-reconstructed-qq.json` is rejected; its
one-bit-below-random LLL vectors were spurious, and it must not be promoted by
characteristic-zero substitution.

This `j`-map belongs to the intermediate `A5+A3+3A1/MW6` fibration.  The
rank-29 and ICARM 398--400 recognition gate is defined for the final rootless
Q80 endpoint.  Comparing those curves with this intermediate map would test
the wrong fibration and is not a substitute for constructing the endpoint
`j`-map.

### Ninth/tenth-prime and intrinsic-basis audit

The adversarial sequence was continued rather than treating `p=163` as an
isolated accident.  Incorporating `p=163` produces the eight-prime candidate

```text
q80-third-q12-j-map-p19-adic-reconstructed-eight-prime-qq.json
sha256 e3625bac3e7c77b42e30eb432f5d6dd3094d6024f069d79fd240ed024e550e78
```

and the independently produced `p=191` model rejects it with the same literal
profile: all 25 numerator pairs and 24 of 25 denominator pairs disagree.  The
replay certificate has SHA-256
`8f3afbe80907698bc21c7255cac93021bc452c0911d1b17b840361c7145bb651`.
After incorporating `p=191`, the nine-prime candidate

```text
q80-third-q12-j-map-p19-adic-reconstructed-nine-prime-qq.json
sha256 dcbc03ee6e032fecbb0d4ede1af2e74939e754c7f4ee0771c5ebe7e81140d949
```

is independently rejected at `p=199`, again with all 25 numerator pairs and
24 of 25 denominator pairs unequal.  That replay has SHA-256
`198518bb83fa9f8667c70ab919cfdb19e25c0ec1f3eba46204c6cfdf4f6ace01`.
The ten-prime candidate including `p=199` remains only an unvalidated lattice
output.  Its degree-eight cube block is two bits below its random boundary and
its denominator-factor block is one bit below; there is no fresh-prime replay
for it.

The reconstruction worker now separates reconstruction primes from repeatable
`--holdout-prime` gates and supports four structured tests:

1. separate projective blocks for the rational and `omega` polynomial
   components;
2. coefficient-pair and scalar reconstruction;
3. exact interpolation after reconstructing the intrinsic degree-eight
   `c4` factor at eight residue-distinct `19^1024` sample values;
4. the base normalization `z=L6(U)/L4(U)` and its further scale fixing by the
   cubic `I2` trace coefficient.

Every alternative is rejected by the untouched `p=199` cube factor.  Thus the
failure is not just the common denominator of the original dimension-17
vector or the pinned base coordinate.  Multiple evaluation points at one
fixed modulus change the lattice basis and expose nonlinear structure, but
once they determine the same coefficient vector modulo `19^1024` they do not
multiply the available p-adic precision.

The long-model audit exposes one additional exact modular factorization:
through `19^1024`, and at all ten aligned finite primes,

```text
H(U) = (U+r)^2.
```

Reconstructing `r` rather than the expanded constant `r^2` is the correct
intrinsic long-model gate, but the nine-prime reconstruction is still rejected
at `p=199`.  This is a negative computation, not a characteristic-zero square
theorem.

### Complete 2048-digit audit and the next justified precision

The restructured `19^2048` run completed.  Its local certificate chain is
pinned by

```text
exact pencil mod 19^2051:                         35670a446fecba6e343f209ce4d6d9bdface5c07134d3c7af99831f7d2cdffc5
factor/root lift mod 19^2048:               33d1c6fc2fca171a13e218d75bd18bac94dbbf4002358d13a1d778c3f6dd28ea
generic integral basis:                     509229eb7aa1834541d7b8fd2013806800478c27876f75b36c53fa8ec59775e8
20-sample residue-distinct manifest:         33488fdc58178e0213c677f76db99bad4808259eac6e8a6a72e170a69bef3172
```

The first label means the exact-pencil artifact used for this run has 2051
p-adic digits; its filename is
`q80-third-q12-exact-pencil-p19-adic-precision2051.json`.  All 20 samples
independently pass the Riemann--Roch dimensions, long Weierstrass identity,
and two-way maps at 2048 digits.  The batch producer now defaults to the
support-determined `--limit 20`; the larger seed inventory is available only
through an explicit larger limit.  Its former human-readable `digits=5`
message was a display bug and now reports the runtime precision.

The full 2048-digit reconstruction audit still fails the untouched `p=199`
cube-factor replay.  The tested formulations include coefficient bundles,
separate rational and `omega` components, coefficient pairs, scalar
convergents, the `I6/I4` and `I2`-trace base normalizations, separately
reconstructed evaluations, one joint lattice across eight residue-distinct
evaluations, all nine projective coefficient charts, and a projective lattice
over the quadratic field with its correct two-coordinate `K^*` scale.  The
nine-prime component candidate has SHA-256
`71a3c70fb44e4bd5f0c5ae3359d117d416f876438e234a26902041a37cd99ea0`;
its independent `p=199` rejection certificate has SHA-256
`dbc37a03ebc9f1745088328c85f2a890040d934d2851e7511d4081a67012428a`.
The quadratic-projective diagnostic candidate has SHA-256
`41f49c133a0876a662f90163f205d00c138947c603fabfc15a0f573cdcf36c0d`.

The failure remains quantitatively random.  With combined reconstruction
modulus bit length 8754, the component lattice returns 7781-bit primitive
coordinates against a 7782-bit random boundary; the coupled quadratic-field
projective lattice returns 7780 bits against the same boundary.  In contrast,
the exact closure operands already record coefficient heights 33,926 and
36,335 bits.  A target of `19^8192` would place the relevant random boundary
near 30,900 bits and is therefore not justified.  The first precision with a
meaningful height margin is `19^12288`, whose boundary is about 46,400 bits.

The exact pencil has now been compiled at precision 12291, SHA-256
`c2df2b3011e1f75c84d359974bc929e51eb47a3aececcf2ffb510b72fa5f2ad4`.
The checkpointed pointwise factor/root lift to `19^12288` completed with
SHA-256
`ab776b418520968f14bc21bdc28ca3c53578e0b5e59c4a54f7a2eb271d3f878a`;
its final checkpoint has SHA-256
`f3c0c8ea90839f4cf26a74f16a6f64dee352f54952f830e2187c829b53f84f37`.
The generic integral-basis check at the five digits used by its theorem passes
with SHA-256
`bd25930a8e5aed900e2add74ae6b4a6b3006c1eea558a0271c233368af9fe502`.

The Riemann--Roch sample worker was then changed internally to replay the
Weierstrass relation in the existing compact conductor-power algebra instead
of repeatedly canonicalizing generic rational functions.  At precision 256
this reduces the same calculation from 212.6 seconds to 2.87 seconds; at
precision 2048 the new output is identical to the old certified sample apart
from the recorded worker hash.  This made the complete set of 20
residue-distinct `19^12288` samples practical.  Every sample passes the four
Riemann--Roch dimensions, long equation, and maps in both directions.  The
manifest has SHA-256
`a6a307360c358ff88be5c0f213ff02c319cfcd0d395d2411c9d128e3eab8e279`.

The extra precision does not validate the characteristic-zero reconstruction.
Both a separate-component candidate (SHA-256
`de1ad26a2c1fb749df82c0ebf205ce7a1a57ecda2694fc02eea820ee438470e7`)
and the coupled quadratic-projective candidate (SHA-256
`fb3ad92b66d3ca807f29483ac5f25bf91493ee9265b8378c02ccf523b76fa868`)
remain at the random-lattice boundary: the degree-eight factor uses 46,446-bit
primitive coordinates against a 46,448-bit boundary.  Independent reduction
at the untouched prime `p=199` rejects both, with all 25 numerator pairs and
24 of 25 denominator pairs unequal.  The replay certificates have SHA-256
`80432448d96e6168c09344bcd32df0483c75dc1a70490ec877c3e6e328385016`
and
`370b414a9d591e9099377a98c0fd25107953144fa180309beb8efc581088e2fa`.

This also shows why the previous height estimate was not a bound on the
normalized `j` coefficients.  In the exact coefficient field used by the
reconstructor, `omega^2=16q1q2` already has a 33,886-bit numerator and a
33,890-bit denominator.  The next reconstruction should therefore retain the
known field and denominator factors symbolically, or use a better integral
field/base normalization, rather than merely increase the same LLL modulus.
The order remains
`j -> (c4^3,Delta) -> minimal Jacobian -> birational maps`.

<!-- status-consumer: EC-K3-Q80-THIRD-Q12-EXACT-LINEAR-CONDUCTOR 163251202819137c -->

## 2026-09-03 exact conductor continuation

The exact 63-term pencil now yields one generic characteristic-zero conductor
factor without reconstructing the full `j`-map.  Rational reconstruction from
the `19^12288` lift gives `L=W+r`; an integral truncated expansion in
`S=denominator(r)W+numerator(r)` proves that the discriminant coefficients at
orders zero, one, and two vanish identically in the remaining base variable,
and that order three is nonzero.  Thus `L^3` divides with exact multiplicity
three.  The exact certificate is
[`../artifacts/generated-results/elkies-k3-q80-third-q12-exact-generic-linear-conductor-v1.json`](../artifacts/generated-results/elkies-k3-q80-third-q12-exact-generic-linear-conductor-v1.json).

The residual exact quartic gcd remains too expensive for the general PARI and
Singular number-field gcd paths tried here.  The pointwise lift nevertheless
reveals a common monic linear denominator `H(V)` in all four nonleading
coefficients of the expected monic quartic.  Direct exact-pencil reductions
at untouched inert primes 163, 191, and 199 reproduce both the discriminant
factor pattern `(1,3),(4,1),(4,2)` and this denominator.  That candidate audit
is
[`../artifacts/generated-results/elkies-k3-q80-third-q12-quartic-denominator-candidate-v1.json`](../artifacts/generated-results/elkies-k3-q80-third-q12-quartic-denominator-candidate-v1.json).
It is not a characteristic-zero proof of `Q`; the next implementation should
be a dedicated modular/subresultant degree-four recovery followed by exact
division, retaining the proved `L^3` and candidate `H` symbolically.
