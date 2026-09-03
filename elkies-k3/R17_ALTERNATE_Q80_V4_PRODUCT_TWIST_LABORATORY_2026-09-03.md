# Alternate-Q80 rational V4 bases and product-twist laboratory

<!-- status-consumer: EC-K3-R17-NORM12-11952-V4-RATIONAL-PAIR-SHORTLIST-64 e14368b602eebedb -->
<!-- status-consumer: EC-K3-R17-NORM12-11952-V4-BASE-RANK-SCREEN-64 f706a4396a0b13af -->

## Exact result

The cheapest-1,024 native alternate-Q80 bisection prefix contains many
arithmetically usable genus-one pair bases despite having no shared branch
fibres.

For a norm-ten bisection

```text
C_w = (2,2,w),
```

the direct-frame intersection is

```text
C_w.C_v = 8 - <w,v>.
```

Exact evaluation of all `523,776` unordered pairs finds `10,362` with
intersection number one.  Two rational curves with intersection number one
have a unique geometric intersection point.  It is fixed by Galois and gives
a rational point on the corresponding `V4` fibre product.

The production certificate selects the first 64 such pairs in the declared
priority/complexity order.  For every pair it also recovers the point from the
two explicit lifted-section formulas.  Substitution verifies

```text
s^2 = q_i(u),
t^2 = q_j(u),
(st)^2 = q_i(u)q_j(u),
```

and verifies that both lifts give the same point of the alternate elliptic
surface.  No bounded rational-point search is used in this proof.

All 64 bases are connected genus-one curves: their two irreducible quadratic
branch divisors are disjoint.  The two existing character sections therefore
give generic rank at least 19 over each base.  The shortlist certificate alone
does not determine the Jacobian ranks.

The exact artifact is
[`../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-pair-shortlist-64-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-pair-shortlist-64-v1.json),
SHA-256
`4739e0b24b00e276269b228d7c85a01743caa1dcbc62ef58ea59ca3152d116a0`.

## Exact base-rank screen

Every base Jacobian was isolated in its own Sage process and given an exact
PARI `ellrank` call with effort zero and a hard ten-second limit.  The bounded
screen gives:

| interval | bases |
|---|---:|
| `[1,1]` | 17 |
| `[0,2]` | 27 |
| `[0,4]` | 8 |
| `[1,3]` | 10 |
| timeout / `UNKNOWN` | 2 |

Thus seventeen shortlisted bases have certified Mordell--Weil rank exactly
one.  No ambiguous interval is treated as an exact rank.  The rank artifact is
[`../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-base-rank-screen-v1.json`](../artifacts/generated-results/elkies-k3-r17-norm12-11952-v4-base-rank-screen-v1.json),
SHA-256
`ed41526fca5724f529aaf73d3f695fe2a781eaf4ea6b6fb792a2af5676759f69`.

## Product-character bisection inversion

The unsliced product-twist `msolve` system is no longer the active route.
The structural inverse problem has now been completed on the entire
norm-eight/pole-order-zero bisection layer.

There are exactly `63,917` section-nonnegative degree-two isotropic
translation classes in the alternate lattice whose trace has norm eight.
For every class the regular chord pencil

```text
M=M0+lambda*h^2
```

was compared with all seventeen product quartics attached to the exact
rank-one pair bases.  A synthetic quartic constructed from the first known
norm-eight pencil is recovered at its declared `lambda=0` before the target
search.  The complete target pass then gives

```text
classes                              63,917
rank-one product targets                 17
modular survivors                         0
exact squareclass hits                     0
```

Prime `131` obstructs `63,915` classes; the two traces with bad reduction in
that calculation are both obstructed at `137`.  The comparison is projective
in the five quartic coefficients, includes `lambda=infinity`, and treats a
zero modular coefficient vector as a survivor.  Thus this is an exact finite
negative, not a bounded runtime observation.

The precise descent dictionary has an integral qualification.  A bisection
with lifts `P,sigma(P)` gives the twist section `P-sigma(P)`.  Conversely an
anti-invariant twist section `T` comes from a bisection exactly when
`T+tau` is divisible by two in the quadratic-cover Mordell--Weil group for
some generic section `tau`.  The cokernel is 2-primary and has not been
proved to vanish.  Therefore the computation excludes every minimal
product-character section in the integral bisection/coboundary image; it does
not yet exclude a non-coboundary height-eight twist section.

The full theorem, proof, certificates, and boundary are in
[`R17_ALTERNATE_Q80_PRODUCT_BISECTION_INVERSION_2026-09-03.md`](R17_ALTERNATE_Q80_PRODUCT_BISECTION_INVERSION_2026-09-03.md).
The older `msolve` export and timeout manifests are retained as historical
complexity evidence only.

## Replay

```bash
sage -python \
  elkies-k3/scripts/select_r17_norm12_11952_v4_pair_shortlist.sage

sage -python \
  elkies-k3/scripts/select_r17_norm12_11952_v4_pair_shortlist.sage --check

.venv/bin/python \
  elkies-k3/scripts/screen_r17_norm12_11952_v4_base_jacobian_ranks.py \
  --limit 64 --timeout 10 --jobs 2

.venv/bin/python \
  elkies-k3/scripts/screen_r17_norm12_11952_v4_base_jacobian_ranks.py \
  --check

sage -python \
  elkies-k3/scripts/rank_r17_norm12_11952_alternate_norm8_pencils.sage \
  --check

sage -python \
  elkies-k3/scripts/search_r17_norm12_11952_product_bisection_inversion.sage \
  --check
```

## Claim boundary

The 64 rational pair-base points, their genus-one geometry, and the seventeen
exact rank-one base Jacobians are established.  The complete norm-eight
bisection inversion has no product-character hit, so no generic
rank-at-least-20 `V4` base change is obtained.  A non-coboundary minimal
product-twist section remains `UNKNOWN` until the 2-primary integral descent
quotient is computed.
