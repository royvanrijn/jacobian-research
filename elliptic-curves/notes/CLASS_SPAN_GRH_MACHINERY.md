# Reusable class-group span certification under GRH

`MATH_STATUS.json`, entry `EC-CLASS-SPAN-GRH-MACHINERY-20260906`, is the
mathematical authority. This generalizes the implementation of the
[MW16 completion](SMALL_CONDUCTOR_CLASS_COMPLETION_PROOF_2026-09-06.md).
The frozen curve proof and its inputs remain unchanged.

The machinery verifies whether proposed ideal classes generate `Cl(K)/2`.
A positive interval test certifies a class-2-rank upper bound under GRH for
the relevant quadratic ordinary ideal-class characters. A valid input with
an inconclusive margin returns `UNKNOWN` and no upper bound. Invalid arithmetic
raises an error and produces no new certificate.

## General sufficient criterion

Let K have degree `n=r1+2*r2` and absolute discriminant D. Let H be the image
of a finite list of ideal classes in `Cl(K)/2`. Fix `T>1`, `L=log T`, and
`F(x)=max(1-x/L,0)`. Include every prime ideal of norm below T.

Exact principal relations certify that some of these prime classes belong to
H. Call these primes known; all others remain unknown. Define

```text
w(P,m) = log(NP) * (1 - m*log(NP)/L) / (NP)^(m/2),  (NP)^m < T,
S_min = sum_(P,m) w(P,m) * sign(P,m),
sign(P,m) = +1 if P is known or m is even; -1 otherwise,
C_upper = log(D) - n*(gamma+log(8*pi)) - r1*pi/2
          + (n*pi^2/2 + 4*r1*Catalan)/L.
```

**Criterion:** If an outward-rounded interval proves
`2*S_min-C_upper > 0`, then H is all of `Cl(K)/2`, under GRH for the
nontrivial quadratic ordinary ideal-class characters trivial on H.

Proof: if H were proper, a nonzero functional on `(Cl(K)/2)/H` would give
such a character chi. At known primes chi is +1; at unknown primes its odd
powers are at least -1 and its even powers are +1. Therefore `S_chi>=S_min`.
The unramified character explicit formula gives
`2*S_chi-C = -sum_rho Phi(rho)<=0` under GRH. The triangular Fourier transform
is nonnegative. Its archimedean constant satisfies `C<C_upper` by replacing
`1-F(x)` with `x/L` in the positive archimedean kernels. Positivity of the
displayed margin is a contradiction.

The source is the classical character explicit formula and triangular bound in
[Belabas–Diaz y Diaz–Friedman (2008), equations (3), (11), (13)](https://doi.org/10.1090/S0025-5718-07-02003-0).
The adaptation uses verified prime-class memberships and pessimistic signs for
the remaining primes. It is not a claimed new general explicit formula.

## Three interfaces

[`class_span_grh.py`](../cas/class_span_grh.py) provides:

- `RelationSpan(width)`: incremental F2 row insertion, followed by
  `analyze(anchor_rows)`. It preserves every prime coordinate, including those
  outside the analytic cutoff, and returns replayed anchor-coordinate witnesses.
- `quadratic_margin(discriminant, signature, prime_norms, known_columns,
  cutoff, bits)`: MPFI interval evaluation for arbitrary number-field signature.
  It also ranks unresolved prime ideals by their odd-power penalty.
- `verify_document(document)`: the arithmetic verifier that combines both
  components and returns `CERTIFIED_UNDER_GRH` or `UNKNOWN`.

The first two are components with mathematical input hypotheses. An arbitrary
parity matrix or a supplied list of known columns does not constitute a proof.
Existing audited relation collectors may call these components after replaying
their own witnesses. The MW16 regression uses that interface and retains the
complete frozen proof as a dependency.

For an incomplete principal-relation space R and anchor vectors A, the verifier
uses `dim((R+<A>)/R)` as an **upper bound** for `dim H`. If generation is
certified, this bounds the class-2-rank. Independence in the formal presentation
is never promoted to actual class independence. A matching lower-bound
certificate is a separate input to any exact-rank conclusion.

## Portable input and arithmetic audit

The JSON schema is `number-fields.class-span-input.v1`. A complete example is
[`imaginary_c2.input.json`](../../artifacts/generated-results/elliptic-curves/class-span-v1/imaginary_c2.input.json).
It contains:

- `field`: ascending coefficients of a monic integral irreducible polynomial,
  optional proved prime hints (an empty list is allowed), signed discriminant,
  signature `[r1,r2]`, and the maximal integral basis in polynomial coordinates;
- `cutoff` and `precision_bits`;
- `columns`: distinct prime ideals, specified by rational prime, ramification
  index, residue degree and HNF matrix in the stated integral basis;
- `anchors`: sparse ideal factorizations `[[column, exponent], ...]`;
- `relations`: nonzero principal generators in polynomial coordinates and their
  complete sparse prime-ideal factorizations. Rational element coefficients
  are exact strings; negative ideal exponents are supported.

The checker proves polynomial irreducibility and requires an empty
`nfcertify` result. It reconstructs every supplied prime ideal, checks all
rational-prime decompositions used, and verifies coverage of every prime ideal
with norm strictly below T. It checks each principal ideal by exact ideal
multiplication. Omitted outside factors cause rejection, even when an omitted
factor would have even exponent. Rational-prime relations are added only when
all prime factors of `(p)` are present. The
[PARI number-field documentation](https://pari.math.u-bordeaux.fr/dochtml/html-stable/General_number_fields.html)
specifies the maximal-order and ideal representations.

The production path uses no `bnfinit` or full class-group computation. Certified
small class groups are used only as independent calibration oracles in
`class_span_fixtures.py`. Large-field maximal-order initialization may still
require factorization work; no fixed runtime is promised.

## Commands and resource bounds

Run from the repository root with Sage on PATH:

```bash
python3 elliptic-curves/cas/verify_class_span_grh.py \
  artifacts/generated-results/elliptic-curves/class-span-v1/imaginary_c2.input.json \
  --output artifacts/local/class-span-example.json

python3 elliptic-curves/cas/verify_class_span_grh.py \
  artifacts/generated-results/elliptic-curves/class-span-v1/imaginary_c2.input.json \
  --output artifacts/local/class-span-example.json --check
```

The CLI supervises a single Sage process, defaults to 120 seconds and 1536 MiB
RSS, and retains a unique run directory with log and supervisor checkpoint.
`--wall-seconds`, `--rss-mib`, and `--run-dir` declare different limits or paths.
Existing certificates are preserved; `--check` recomputes and compares them.
Timeout or invalid arithmetic cannot produce a success certificate.

The fixed validation commands are:

```bash
sage -python elliptic-curves/cas/certify_class_span_machinery.sage launch-check
sage -python elliptic-curves/cas/certify_class_span_machinery.sage launch-check-mw16
```

The second command additionally needs the frozen MW16 evidence. Both stages
have fixed 600-second and 1536-MiB limits in the validation protocol.

## Validation and use in rank searches

The [validation certificate](../../artifacts/generated-results/elliptic-curves/class_span_machinery_v1.json)
records seven small-field cases of degrees 2, 3 and 4, with signatures
`(0,1)`, `(2,0)`, `(1,1)`, `(3,0)` and `(0,2)`, and certified
class-2-ranks 0, 1 and 2. Proper-span controls must remain `UNKNOWN`. A positive
control deliberately omits the relations at 97 in Q(i), leaving a nonzero
formal quotient while the analytic test certifies class-2-rank at most zero.
Tests also reject missing prime coverage, duplicate/forged primes, incorrect
bases, incomplete principal factorizations and inexact inputs. Random finite
row spaces are compared with Sage matrix ranks.

The [MW16 regression](../../artifacts/generated-results/elliptic-curves/class_span_mw16_regression_v1.json)
replays the existing arithmetic and checks that the general code reproduces
the 4,740 prime memberships below 50,000, upper bound 16 and positive margin.
The small standalone fixtures and general checker can be packaged separately
from that large historical dependency.

For another elliptic curve, feed audited relations from its own cubic field,
then combine a successful class bound with its separately certified local
Selmer correction, applicable Selmer parity, and independent rational points.
Matching bounds prove exact rank under the stated assumptions. If the bounds
do not match, the pipeline remains inconclusive. Unresolved-prime penalties
can guide the next bounded relation batch; they do not guarantee a successful
relation or a rational point. This is reusable rank-proof infrastructure,
not a universal fast or unconditional rank algorithm.
