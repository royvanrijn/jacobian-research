# Minimal-resolution positivity and the r7 Picard-rank window

14 November 2026. **Subsequent result:** the
[adjoint-fibration proof](../r7-adjoint-exclusion/README.md)
sharpens the rank to eight and excludes both coefficient classes under
PF2R7BND1's hypothesis and the explicit Borisov source-model inputs. The bounds
below remain retained necessary conditions, not the current continuation gate.

Under the explicit degree-six, single-ramified-row hypothesis, the earlier
necessary conditions for either nonimmersed normal-r7 class were:

- general target-line pullback genus **two**;
- minimal-resolution Picard rank **8, 9, 10 or 11**;
- exactly four possible generic infinity decompositions;
- **two to seven** contracted boundary primes, after retaining the two or
  three unramified omitted affine primes already forced by PF2R7BND1.

[r7 proof and complete numerical table](R7_WINDOW.md) (`JC2POS3`). No r7
class or additional-image configuration is excluded outright.

The [general proof](PROOF.md) records `JC2POS1` and `JC2POS2`:

- reproduces nef/big full ramification and canonical adjoint on the actual
  minimal resolution;
- removes the original conditional smooth-finite/cubic-scroll donor
  dependency by a direct boundary/Euler/Riemann--Hurwitz argument;
- excludes general-line genus one and Picard rank one;
- gives exact smooth blowdown rules, retaining arithmetic genera after
  the boundary ceases to be SNC, and boundary intersection tests.

The campaign's original theorem and proof strategy retain their attribution.
Its source bytes and dependency notes are [pinned](sources/manifest.json).
This is written geometry supported by exact regressions, not formal
verification or independent external review; novelty priority is unassessed.

```sh
timeout 60s .venv/bin/python research/jacobian-followup/minimal-resolution-positivity/verify.py
```

[Checker](verify.py) · [implementation](positivity.py) · [actual receipt](receipt.json).
The receipt separates the conditional r7 deductions from applicability:
the live quartic/F2 packets still lack a complete global surface model;
the older 23-component transformed Laurent chart has Jacobian X^2 and is
ineligible as a Keller application. Its numerical reduction passes the
new positivity tests, so no historical exclusion is inferred.

[MATH_STATUS.json](../../MATH_STATUS.json) remains authoritative. No proof
of JC2, HC4, (75,125), normal r5, immersed r7 or nonnormal-slice exclusion
is claimed.
