# Complete two-prime cohomological triage of the seventeen product twists

**Arithmetic update (2026-09-05):**
[`19bad:083ad` has rank zero over `QQ(u)`](R17_PRODUCT_19BAD_083AD_ARITHMETIC_RANK_ZERO_2026-09-05.md),
proved from incompatible regulator squareclasses in its existing two
reductions. Four of the five Frobenius survivors remain arithmetic
candidates. The geometric classification table below is unchanged.

## Result

The full degree-28 cohomological quotient has now been computed for all
seventeen selected alternate-Q80 product twists.  Twelve have a good
reduction with no Tate factor and therefore have geometric Mordell--Weil rank
zero.  Five retain a degree-two Tate factor at both `p=131` and `p=137`; for
those five the unconditional conclusion is only

```text
0 <= rank E^(q_i*q_j)(QQbar(u)) <= 2.              (1)
```

They were the five cases from this cohomological gate justified for full
descent or class-sliced section construction; the arithmetic update above
removes `19bad:083ad`. A Tate factor at two primes does not prove
a characteristic-zero section.

The aggregate certificate is
[`elkies-k3-r17-all17-product-toric-frobenius-campaign-v1.json`](../artifacts/generated-results/elkies-k3-r17-all17-product-toric-frobenius-campaign-v1.json).

## Classification

The entries in the two prime columns are cyclotomic Tate degrees counted with
multiplicity.  A dash means the second prime was unnecessary.

| pair | `p=131` | `p=137` | unconditional geometric conclusion |
|---|---:|---:|---|
| `1463f:19bad` | 4 | 0 | rank 0 |
| `19bad:083ad` | 2 | 2 | rank at most 2 |
| `11ae6:0f82c` | 2 | 2 | rank at most 2 |
| `0f82c:025be` | 2 | 2 | rank at most 2 |
| `025be:13dbe` | 2 | 0 | rank 0 |
| `1ad20:1b24d` | 0 | — | rank 0 |
| `1a465:19b4e` | 2 | 0 | rank 0 |
| `19ead:146dc` | 0 | — | rank 0 |
| `19b4e:17b71` | 0 | — | rank 0 |
| `11ee2:0c36e` | 2 | 2 | rank at most 2 |
| `0c36e:02bf1` | 0 | — | rank 0 |
| `0fda0:1a6c8` | 4 | 0 | rank 0 |
| `1059f:1db8d` | 0 | — | rank 0 |
| `0fda0:1037d` | 0 | — | rank 0 |
| `0c10b:17a1a` | 2 | 2 | rank at most 2 |
| `1ede3:1c364` | 2 | 0 | rank 0 |
| `13dbe:1019b` | 0 | — | rank 0 |

The five persistent survivors are therefore

```text
alternate-orbit-19bad : alternate-orbit-083ad
alternate-orbit-11ae6 : alternate-orbit-0f82c
alternate-orbit-0f82c : alternate-orbit-025be
alternate-orbit-11ee2 : alternate-orbit-0c36e
alternate-orbit-0c10b : alternate-orbit-17a1a.
```

Four have normalized Tate factor `(Z-1)(Z+1)` at both primes.  The
`0f82c:025be` target has `(Z-1)^2` at `131` and `(Z-1)(Z+1)` at `137`.
These patterns remain upper-bound data only.

## Certificate chain

Every reduction uses the regular equation

```text
d(u)y^2=x^3+A(u)x+B(u)
```

and passes the exact degree, squarefreeness, discriminant-coprimality, Newton
support, and toric nondegeneracy gates.  From the raw controlled-reduction
output the independent verifier reconstructs

```text
P_boundary=P_Z/P_D,                           degree 8,
P_triv=(T-p)^2 P_D P_Z,                       degree 18,
P_H2=(T-p)^2 P_D^2 P_toric=P_triv P_E,        degree 46,
P_E=P_toric/P_boundary.                       degree 28.       (2)
```

It then checks integrality, the reciprocal functional equation, exact Weil
size by certified real-root isolation, and the independent fibrewise `n=1,2`
moments.  Finally it computes the gcd with every `Phi_m` satisfying
`phi(m)<=28`, including factor multiplicity.

For each of the twelve zero-Tate targets, (2) gives `rho(Fpbar)<=18`; the
explicit `U+4D4` lattice gives equality.  Specialization and Shioda--Tate then
prove product-twist rank zero without the Tate conjecture.  On the associated
torsion-free `48I1` quadratic base change this also gives

```text
A-=0,       Gamma_d=0,       Hhat^(-1)(<sigma>,A)=0.          (3)
```

Thus the earlier nonzero Tate-class loophole is gone for twelve targets, not
only for shortlist rank 55.  It remains genuinely open for the five
persistent survivors.

## Replay and next gate

```bash
JOBS=3 elkies-k3/scripts/run_r17_all17_product_toric_frobenius_campaign.sh
sage -python \
  elkies-k3/scripts/certify_r17_all17_product_toric_frobenius_campaign.sage \
  --check
```

The first command computes `p=131` for all targets, sends only cyclotomic
survivors to `p=137`, and then builds the aggregate.  The second is the cheap
fail-closed replay over the stored individual certificates.

The next expensive arithmetic work should be restricted to the five listed
survivors.  Shortlist rank is not a useful product-rank predictor: the
rank-55 target is zero, while survivors and zero targets are interleaved
throughout the stored order.

<!-- status-consumer: EC-K3-R17-PRODUCT-19BAD-083AD-ARITHMETIC-RANK-ZERO fe572bd5979b5d2c -->
