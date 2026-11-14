# Scoped exclusion of both nonimmersed normal-r7 classes

**Later route:** [PF2R7EX2](../jc2-r7-quartic-contact/README.md) now proves
the same scoped two-class exclusion without source finality. This
directory preserves the earlier source-relative proof and certificates.

14 November 2026. The written proof excludes both retained nonimmersed
normal-r7 coefficient classes **under the explicit degree-six,
single-ramified-affine-row hypothesis of PF2R7BND1**. Unramified omitted
primes are retained throughout. This does not close the immersed r7
family, other ramification patterns, r5, nonnormal slices or JC2.

- `PF2R7ADJ1`: the genus-two adjoint system forces **Picard rank8** and a
  basepoint-free rational pencil |K+H| with H-fibre degree2. The cusp's
  two-sheet local alternative is excluded by the exceptional-curve count.
- `PF2R7PAR1`: exact univariate certificates exclude every constant-leading
  parabola family and every Mobius two-contact involution for either
  retained parametrization. Fraction replay independently checks the
  Bezout identities.
- `PF2R7EX1`: the geometric bridge reduces the remaining cusp alternative
  to precisely those algebraic mechanisms, giving the scoped exclusion.

[Complete proof](PROOF.md) · [external source/model audit](SOURCE_AUDIT.md) ·
[author-side obligation review](REVIEW.md) · [exact certificates](certificates.json).

The external bridge uses Borisov's type1 residue-degree theorem and the
final-curve assertion in his proof of Theorem3.8. The latter is explicitly
cited, not independently reproved here. Its transport to our minimal
resolution is proved in the follow-up. Finite algebra does not certify
that geometric input. No formal verification, independent external review
or novelty priority is claimed.

```sh
timeout 60s .venv/bin/python research/jacobian-followup/r7-adjoint-exclusion/verify.py
```

[Checker](verify.py) · [actual receipt](receipt.json).
The small replay performs no new braid computation, graph enumeration or
large r7 elimination. Mathematical status is recorded in
[MATH_STATUS.json](../../MATH_STATUS.json).
