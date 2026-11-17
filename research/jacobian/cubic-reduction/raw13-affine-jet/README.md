# KCUB41 full affine-orbit controls

The [canonical proof](../RAW13_AFFINE_ORBIT_MINIMUM.md) proves raw13 sharp
for every affine source/output frame of the saved map,with equality only
at the original source centre and monomial linear frames. Combining KCUB40
extends the13 minimum to every affine frame of every degree<=7 diagonal-weight
collision. Global unweighted/nonlinear-coordinate optima remain open.

- [discovery.json](discovery.json): six exploratory full matrices and three
  preliminary two-variable rank-ideal unit answers,explicitly unused as proof.
- [input.json](input.json): six complete translated matrices,348 entries,
  and two selected-minor units using31 multiplier monomials.
- [checkpoint.json](checkpoint.json): completed production with source pins.
- [validation.json](validation.json): independent six-variable translation,
  all minors and unit multiplication,covariance,exact13/14 rank sums,
  eighteen nonzero-axis columns and the saved determinant/collision control.
- [integration-validation.json](integration-validation.json): current source,
  metadata,review and navigation-check boundary.

Reproduce from the parent directory with `python raw13_affine_jet.py` and
`python verify_raw13_affine_jet.py`. Neither makes solver calls. The two
short units replace preliminary Groebner answers. Universal invariance,
field coverage and UFD equality arguments are written,not formally verified
or externally reviewed. No new positive witness or global bound is claimed.
