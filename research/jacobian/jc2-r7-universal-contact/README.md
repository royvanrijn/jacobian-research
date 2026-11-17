# The entire gap-21 family: contact obstruction and three infinity escapes

**Subsequent refinement:** [cubic conic coefficients](../jc2-r7-cubic-cover/README.md)
exclude (2,4) and restrict the reciprocal case. The
[quartic-contact proof](../jc2-r7-quartic-contact/README.md) then excludes
both cusp classes without source finality. The
[universal extension](../jc2-r7-universal-quartic/README.md) closes the
entire quartic family, leaving only the special reciprocal cubic case.
The three-escape theorem below is preserved as a necessary-window dependency.

14 November 2026. [PF2R7PAR2](PROOF.md) extends the parabola algebra from
two cusp examples to every birational polynomial degree-(6,10), gap-21
pair. A literal gap certificate rules out fourth-order roots of p'; two
symbolic contact coefficients exclude every Mobius involution. The
constant-leading rational-function degree obstruction is also uniform.
The shifted-product pair (p-r,(p-r)*q) also has no contact involution.
Together with a pole-order bound, this removes two reciprocal-conic subcases.

PF2R7INF1 extends the rank-eight adjoint pencil to the whole family under
the explicit entire single ramified affine-prime (e,f)=(2,1), degree-six
Keller hypothesis. **Without source finality**, it excludes the balanced
infinity pair (3,3) and leaves exactly three horizontal-infinity budgets.
Their ramified projection degrees are now exactly 3,4,3; two of their
generic conic families and critical equations are explicit in the proof.
With the additional Borisov source-finality input already isolated in
the [source audit](../r7-adjoint-exclusion/SOURCE_AUDIT.md), those remaining
cases are excluded too, including the immersed family. That stronger
source-relative corollary is not an independent reproof of finality.

```sh
timeout 60s .venv/bin/python research/jacobian/jc2-r7-universal-contact/verify.py
```

[Proof and escape table](PROOF.md) · [exact constants](certificates.json) ·
[replay](receipt.json) · [validation](validation.json).

Additional ramification patterns, normal r5, nonnormal slices and JC2
remain unresolved. No formal verification, independent external review
or novelty priority is claimed.
