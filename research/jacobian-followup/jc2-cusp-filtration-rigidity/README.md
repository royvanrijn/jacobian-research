# Cusp-filtration rigidity adopted into the planar programme

14 November 2026. A written audit of the public JC2 campaign's coprime
cusp-filtration theorem, followed by a finite recognition criterion and
a concrete necessary contact bound for the retained F2 row.

- **JC2CFR1:** reproduce the all-coprime same-filtration Euler/Poisson
  rigidity proof, preserving its source attribution and assurance limits.
- **JC2CFF1:** prove the sharp cusp-Jacobian degree inequality and reduce
  whole-ring scaling to a finite polynomial-degree certificate. The rigid
  classification also covers distinct source and target binomial cusps.
- **JC2CFF2:** in the retained F2 support envelope,
  delta_(2,5)(P^5+(125/729)*Q^3)>=456 for any actual Keller pair.
  An upper bound <=455 would exclude the candidate; that bound is unknown.

The [Wronskian follow-up](../jc2-cusp-contact/README.md) classifies the exact two-face F2 contact spectrum, whose minimum is507, and the adapted-quintic Euler obstruction. The earlier456/457/459 bounds remain valid but are weaker.

The [complete proof](PROOF.md) separates reproduction, new deductions and
the exact application gate. The [pinned source](sources/manifest.json)
retains the upstream MANUAL / PROVISIONAL status. The
[checker](verify.py) and [receipt](receipt.json) replay rational polynomial
identities, complete finite preserving-derivation kernels, graded symbols
with no eligible polynomial lift, and the two F2 support envelopes.

From the repository root, no third-party Python package is needed:

    python3 research/jacobian-followup/jc2-cusp-filtration-rigidity/verify.py

The checker includes a reusable exact sparse normal-form engine and
diagnose(P,Q,a,b,p,q): source polynomials are dictionaries from (x,y)
exponents to rational coefficients. Its output records component degrees,
cusp-difference degree, Jacobian degree and defect. Source and target cusp
exponents must be coprime and at least two. A zero Jacobian is explicitly
rejected as algebraic dependence.

The all-degree statements rely on written proofs, not finite regressions.
No formal verification, independent end-to-end implementation, external
review or novelty priority is claimed. The F2 contact upper bound, r7
source-contact estimate, arbitrary valuation reduction, JC2 and HC4 remain
open. No historical missing artifact was reconstructed. The [adoption audit](validation.json) records passing scoped checks and the pre-existing missing-artifact failures in repository-wide navigation validation.
[MATH_STATUS.json](../../MATH_STATUS.json) remains authoritative.
