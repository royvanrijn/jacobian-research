# A concrete star510 point for the isotropic-parent descent gate

The subsequent [level51 descent obstruction](STAR510_FRICKE_LIFT_EXCLUSION_2026-09-15.md)
proves that this non-CM point has no rational Fricke lift. The proposed
parent route below is closed; its input and earlier UNKNOWN states are
retained as the historical construction attempt.

**Exploratory arithmetic input; no new K3 parent is admitted.** The
[real-locus obstruction](REAL_LOCUS_C4_ANISOTROPY_OBSTRUCTION_2026-09-15.md)
requires isotropic T for the order-four mechanism on a fully rational
Picard19 K3. A newly checked modular point supplies a specific next descent
question for that route. The [correlated-gain objective](CORRELATED_QUADRATIC_GAINS_2026-09-12.md)
remains open.

Assaf–Hashimoto–Shnidman report an exceptional point on X0*(510) in
[their September2026 paper, Table1](https://arxiv.org/html/2609.00516v1#S2.T1).
Their [companion repository](https://github.com/sachihashimoto/AtkinLehnerQuotients)
at commit65c90892f05458396f7c0e649737e612b1f7a331 records the point(8:3:1).
Its non-CM identification is inherited from their work, not independently
reproved here.

## Correct coordinate model

The authors' general genus-three model file and point log use different
bases of differentials. Substituting(8,3,1) into the former quartic gives-198,
not0. This is a coordinate mismatch, not a refutation of their point.
From the cached three star forms, exact linear algebra instead recovers

```
2x³z - 2x²y² - 14x²z² + xy³ + 25xy²z - 17xyz² - xz³
-12y⁴ + 16y³z - 4y²z² = 0.
```

Here P=(8:3:1) lies on the curve and has gradient(72,-376,552).
On y=3z its equation factors as

```
2z (x-8z) (x²-8xz+36z²).
```

The [input](../artifacts/generated-results/elkies-k3-star510-fricke-input-v1/input.json)
retains the first865 coefficients of each author-cached form, upstream hashes,
commit, quartic and point. The [checker](scripts/verify_star510_fricke_input.py)
checks that the degree-four relation space from the first80 coefficients is
one-dimensional and that the displayed relation vanishes through coefficient864.
These are weight-eight forms on Gamma0(510), whose index1296 gives Sturm
bound864. Modularity and the meaning of the cached basis are inherited;
no new modular-form calculation is claimed. Unit Jacobian ideals in the affine
and infinity charts modulo7, plus the remaining projective point, certify
geometric smoothness of the quartic. A finite list of rational residue points
alone would not suffice for that assertion.

## The missing arithmetic lift

For squarefree510=2*3*5*17, the Atkin–Lehner group has order16. The map

```
X0(510)/<w510>  ->  X0*(510)
```

has degree8, with group(C2)^3. Over a non-CM point this is an etale torsor:
a fixed point of a nontrivial squarefree Atkin–Lehner involution would give a
non-scalar self-isogeny of nonsquare degree, hence CM. A rational point on
the star quotient does not assert that this torsor has a rational point.
The exact fibre above P must be computed and tested for a rational lift.
No such fibre algebra or lift is supplied by the retained input.

A rational Fricke lift would give a Galois-stable unordered pair of
non-CM elliptic curves joined by a cyclic510-isogeny. The
[Inose descent argument used at level311](DET622_INOSE_RATIONAL_MARKING_SOURCE_2026-09-15.md)
would then be a possible route to a fully rational marking of

```
NS = U + E8(-1)² + <-1020>,    T = U + <1020>.
```

This is an isotropic NS candidate, distinct from the already admitted
anisotropic determinant1020 NS. The equal determinants must not be used
to transfer a rootless-frame certificate.

Indeed the existing positive frame has cyclic discriminant generator value
10877/1020, whereas the desired frame genus has value1/1020. An isometry
would require a unit u with u²=10877 modulo2040. Reducing modulo5 gives
u²=2, impossible. The checker verifies this obstruction against the actual
retained frame certificate. A suitable rootless frame for the isotropic NS
therefore remains a separate UNKNOWN, even if the modular lift succeeds.

The next concrete step is the rational fibre of the degree-eight Fricke map
above this fixed point. Neither the full-quotient point alone nor the
wrong-genus determinant1020 witness authorizes an MW17 claim. Fixed-point
field, one new section and the infinite rational conic would still be needed
for the order-four two-gain construction.

## Replay and scope

The [result](../artifacts/generated-results/elkies-k3-star510-fricke-input-v1/result.json)
records the exact checks and the two UNKNOWN gates. Upstream data is retained
under its [MIT license](../artifacts/generated-results/elkies-k3-star510-fricke-input-v1/LICENSE.upstream.txt).
No author scripts, point searches, modular-curve rebuilds or frame campaigns
were run. The cached q-expansions were reused to repair the coordinate mismatch.
Replay takes under one second with20 CPU seconds and1GiB allowed:

```
.venv/bin/python research/elkies-k3/scripts/verify_star510_fricke_input.py
```
