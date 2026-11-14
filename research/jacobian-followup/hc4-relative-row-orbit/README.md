# The complete relative-linear row orbit and its missing-boundary gate

The unbounded relative-linear Schur route has an exact reformulation:
its possible plane coefficient pairs are precisely the ordered generators
of (x,y) with Jacobian determinant -2 at the origin. A successful Keller
pair therefore has a single reduced zero fibre. Proving invertibility
still requires control of the missing finite-normalization fibre at zero.

The [proof](PROOF.md) records two results:

- **HC4MYOR1:** starting from HC5T1's retained degree-ten normalizer, the
  full row orbit and its one-polynomial stabilizer are classified.
  Arbitrarily high matrix degrees can give exactly the
  same linear output pair.
- **HC4MYOR2:** the unit-Hessian/collision question in this relative-linear
  class is exactly a singleton-zero plane Keller problem. Properness at
  zero would close it. A normalized degree-99 sufficient gate strictly
  extends the old degree-89 matrix box, using only Moh's degree-100 theorem.

An explicit known three-variable Keller map demonstrates why a singleton
source fibre does not imply properness: two other sheets escape to infinity.
Its finite cubic normalization has fibre lengths 1+2 over the origin.
This is a control in dimension three, not a plane counterexample.

[Checker](verify.py) · [actual replay](receipt.json).

```sh
timeout 60s .venv/bin/python research/jacobian-followup/hc4-relative-row-orbit/verify.py
```

The unrestricted relative-linear class, the broader Schur problem and JC2
remain open. These results use no Tang theorem. No formal verification,
external review or novelty priority is claimed. Status is authoritative in
[MATH_STATUS.json](../../MATH_STATUS.json).
