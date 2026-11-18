# A fixed hyperbolic three-cusped parent for uniform primeness

Date: 18 November 2026. The fixed three-component link below is explicitly defined and rigorously verified hyperbolic. The [threaded-curl diagram theorem](UNBOUNDED_PRIME_AUGMENTATION_GAP_FAMILIES.md) identifies the all-bridge parents with its fillings and proves prime augmentation gaps at every sufficiently large bridge number. This note supplies the indispensable fixed-parent geometry and disk-surgery marking.

Let `U` be the fourteen-letter prefix in [UNIFORM_PROTECTED_AUGMENTATION_GAPS.md](UNIFORM_PROTECTED_AUGMENTATION_GAPS.md), and set

\[
\begin{split}
R={}&\sigma_4\sigma_3\sigma_2^2\sigma_3\sigma_4\sigma_5^{-1}
          \sigma_4 U\sigma_4^{-1}\sigma_2\sigma_3\sigma_5,\\
D={}&\sigma_6\sigma_5\sigma_4^2\sigma_5\sigma_6^{-1},\qquad
N=S^3\setminus\widehat{RD}.
\end{split}
\]

This is a seven-strand 32-crossing presentation of a three-component link. Denote its components by the knot `K`, the old circle `C` and the new circle `D_0`.

## 1. Proved geometry and surgery of the fixed candidate

The word `R` is obtained from the bridge-five parent with its one full circle twist removed: start with

\[
Q'=\sigma_4U\sigma_4^{-1}\sigma_2\sigma_3,
\qquad P=\sigma_5\sigma_4\sigma_3\sigma_2^2\sigma_3\sigma_4\sigma_5^{-1},
\]

and cyclically cut `Q'P` just after the first suffix letter `sigma_5`. At that cut the old circle is in position five. The new circle word `D` is `A_(4,7) A_(5,7)`, so it encircles exactly the adjacent positions four and five: one knot arc and one old-circle arc. It bounds an embedded disk pierced once positively by each, and has no self crossings. In the retained PD component order `(D_0,K,C)`, the linking matrix is

\[
\begin{pmatrix}0&1&1\\1&0&3\\1&3&0\end{pmatrix}.
\]

The two selected-component roles are bound to the literal PD by every braid port and all seven closure identifications. The component edge counts are `6,48,10` in that order. Both circle projections have no self crossings.

The same geometric disk surgery argument as in the five-bridge proof gives, in its preferred meridian/disk-longitude marking,

\[
N\bigl(D_0:(1,-q)\bigr)
\cong S^3\setminus\widehat{R\sigma_4^{2q}},\qquad q\in\mathbb Z.
\]

For completeness, at the new-circle slab the two selected meridians `a,b` are conjugated by `c^{-1}`, and the circle output is `(abc)^{-1}c(abc)`. Closing it gives `[c,ab]=1`; filling `(1,-q)` gives `c=(ab)^q`, exactly the action of `q` positive two-strand full twists. The spanning-disk twist supplies the homeomorphism and fixes the old components outside the slab. The exact free-group control is retained, so the sign is not guessed from numerical cusp shapes.

The complete exterior of the literal three-component diagram has a retained 35-tetrahedron triangulation with three cusps and decorated signature

```text
JLLvLLvLwzPvQPPAAvwQPQQQccdiijnkoppstxzxrutyBwCBACBCEFIGFHIHIicgvseagggvftilbhbtabqfgivmvfnbgbggo_bcaCdbBBakbBbcB
```

The pinned Sage/SnapPy replay freshly verifies Krawczyk interval hyperbolicity at 100 bits on both the PD exterior and the reloaded triangulation. It also checks the face pairings, retained rational boxes, logarithmic gluing and cusp-completeness equations and positive orientations. No randomization was needed. Thus this `N` is a complete finite-volume hyperbolic three-cusped manifold under the stated software imports. Its approximate volume `31.0249194173386` is descriptive only.

## 2. The uniform identification

In the uniform all-bridge notation, remove one circle full twist and write

\[
Q'_b=\rho_bU\rho_b^{-1}\theta_d,
\quad \theta_d=\sigma_2\cdots\sigma_{b-2},\quad d=b-2.
\]

The resulting marked-circle parent is homeomorphic to `M_b`, with a changed circle peripheral marking. The [uniform diagram lemma](UNBOUNDED_PRIME_AUGMENTATION_GAP_FAMILIES.md) proves

\[
S^3\setminus\widehat{Q'_bP_b}
\cong
N\bigl(D_0:(1,-(b-5))\bigr),
\]

with the knot and old-circle components identified. An unmarked exterior homeomorphism suffices to establish parent hyperbolicity for all sufficiently large `b`. A fully marked simultaneous two-parameter surgery formula would additionally track the old-circle longitude correction; the prime-family theorem avoids that unproved marking assertion.

The main knot has the exact braid identity

\[
\sigma_3^{-1}(U\delta_b)\sigma_3
=(\sigma_3^{-1}U\sigma_2\sigma_3)
       \sigma_4\sigma_5\cdots\sigma_{b-1}.
\]

Thus it is an iterated positive Markov stabilization of one fixed four-strand knot. Every added strand pierces the old circle's disk. The port proof uses two finite oriented disk cells to turn these threaded curls into positive clasps between the circle and the knot. Its source-ladder formulas, stage induction and fixed four-port target contract prove the isotopy for every bridge number, preserving the entire fixed `U` block. This argument supplies the link identification; main-knot Markov equivalence and numerical manifold recognition alone would not suffice.

## 3. Retained evidence

The [input](protected_universal_parent_candidate_inputs.json), [checker](verify_protected_universal_parent_candidate.py), [finite receipt](protected_universal_parent_candidate_verification.json), [rigorous geometry receipt](protected_universal_parent_candidate_geometry_verification.json) and [source ledger](protected_universal_parent_candidate_source_ledger.json) retain the exact fixed link, component labels, disk-surgery convention and analytic certificate. Its coverage flag explicitly excludes certification of the uniform identification, which belongs to the separate port checker. Corrupted new-circle detours and component swaps are rejected.

```bash
.venv/bin/python research/meridional-rank/representations/verify_protected_universal_parent_candidate.py \
  --write-receipt research/meridional-rank/representations/protected_universal_parent_candidate_verification.json
~/.local/share/jacobian-sage-10.9/bin/python \
  research/meridional-rank/representations/verify_protected_universal_parent_candidate.py --sage \
  --write-receipt research/meridional-rank/representations/protected_universal_parent_candidate_geometry_verification.json
```

The fixed geometry checker certifies this `N` and its disk surgery, while the separate exact port proof establishes the all-bridge transfer. No formal verification, independent external replay, external review or novelty judgment is claimed.
