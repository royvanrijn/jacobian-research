# F2 `(75,125)` common-edge boundary handoff

> **Status: exact contact classification and exact transfer obstruction.**
> The quadratic common-root family determines four contact packets on one
> selected toric divisor. It does not determine toroidal branch scales,
> dicritical components, or finite-normalization rows. Even the strongest
> naive promotion of contacts to ramification rows survives the available
> coarse boundary gates. This note therefore stops the F2 layer/boundary
> route; it does not exclude `(75,125)`.

The exact replay is
[`cas/audit_f2_75_125_boundary_handoff.py`](cas/audit_f2_75_125_boundary_handoff.py).
Its pinned artifact is
[`../artifacts/generated-results/jc2_f2_75_125_boundary_handoff.json`](../artifacts/generated-results/jc2_f2_75_125_boundary_handoff.json).

## 1. Input from the Laurent classification

The upper-band calculation in
[`F2_75_125_DERIVATION.md`](F2_75_125_DERIVATION.md) proves

\[
H(t)=(1+u+u^2+u^3+u^4)^2R(u^5),\qquad u=1+t,
\]

where

\[
R(v)=av^2+bv+\left(\frac1{25}-a-b\right),\qquad a\ne0.
\]

The common edge polynomial is therefore

\[
C_0(u)=t^7H(t)=(u-1)^5(u^5-1)^2R(u^5). \tag{1}
\]

It has degree 25. Since \(R(1)=1/25\), no root of `R` collides with the
fixed fifth-root-of-unity packet. Equation (1) always has:

- contact multiplicity seven at \(u=1\);
- contact multiplicity two at each of the other four fifth roots of unity.

The two roots of `R` provide the remaining contact degree ten.

## 2. Exhaustive root strata

Put

\[
c=\frac1{25}-a-b,\qquad \Delta=b^2-4ac.
\]

Because `R` is quadratic and \(R(1)\ne0\), the four rows below are disjoint
and exhaustive over an algebraically closed characteristic-zero field.

| `R` stratum | centers | contact partition of 25 |
| --- | ---: | --- |
| \(c\ne0,\ \Delta\ne0\) | 15 | \(7,2^4,1^{10}\) |
| \(c\ne0,\ \Delta=0\) | 10 | \(7,2^9\) |
| \(c=0,\ b\ne0\) | 11 | \(7,5,2^4,1^5\) |
| \(c=b=0,\ a=1/25\) | 6 | \(10,7,2^4\) |

For a nonzero root \(\rho\) of `R`, the equation \(u^5=\rho\) has five
distinct solutions; its multiplicity is the multiplicity of \(\rho\).
A root of order \(m\) at \(\rho=0\) instead gives contact order \(5m\) at
\(u=0\). This proves every row without numerical root approximation.

## 3. Why contacts are not ramification rows

The Newton/boundary dictionary identifies an edge-root multiplicity with
the contact multiplicity of the leading restriction on the selected toric
divisor. It does **not** identify it with:

- a toroidal branch scale;
- a source boundary prime in the finite normalization;
- a transverse ramification index \(e\);
- a residue degree \(f\);
- a puncture count;
- or a group of points in one target fiber.

This distinction is unavoidable. In smooth local coordinates `(s,z)`, all
germs

\[
g_q(s,z)=s^m+z^q,\qquad q=1,2,3,
\]

have the same edge restriction \(g_q(s,0)=s^m\), but their primitive
equality rays are

\[
\left(\frac m{\gcd(m,q)},\frac q{\gcd(m,q)}\right).
\]

Thus the first nonzero normal order—lower-band data—changes the entire
toroidal fan while leaving the contact packet unchanged. Feeding the
partition of 25 directly to the log-boundary compiler would invent its
input.

## 4. Strongest naive finite-normalization test

For diagnosis, the checker deliberately makes the unsupported assignments

\[
(e_i,f_i,s_i)=(m_i,1,1)
\]

for every contact multiplicity \(m_i\). Their boundary contribution is
always 25. Since a target nonproperness curve must retain a positive affine
sheet, the smallest compatible generic degree in this surrogate is

\[
d=25+1=26.
\]

All four resulting coarse signatures are residue-immersion-compatible. If
one additionally—and again without proof—places every center in a single
closed target fiber, its minimum packet length is 25, while the finite-flat
fiber length is 26. Hence every row receives the same exact verdict:

> fiber length permits this packet; further data are required.

The typed target-ledger audit remains `incomplete`, because no source
contact has been transferred to a normalization-boundary prime and neither
the boundary nor affine generic pullback ledger is exhaustive.

This surrogate is intentionally biased toward producing a contradiction.
Its survival proves that contact arithmetic alone cannot close F2 through
the current finite-normalization or conductor-packet inequalities.

## 5. Pivot verdict

The boundary route has not converted the F2 family into an exclusion:

- all four contact strata survive;
- the contact partitions do not determine local blowup scales;
- the coordinate degree 125 does not determine geometric degree \(d\);
- no target curve, target fiber, residue degree, or puncture profile is
  certified;
- consequently neither the log-boundary matrix nor the residual-different
  identity can be constructed honestly.

The missing datum is again the first nonzero normal term at each root
center. Recovering it requires the lower Laurent bands that the pivot was
intended to avoid. Continuing the remaining thirty layers merely to obtain
that datum would return to the rejected sequential strategy.

The appropriate conclusion is therefore:

> Stop the degree-specific F2 descent. A future bridge must be a
> degree-independent theorem converting a common-power edge packet into
> target-side finite-normalization data, or it must attack the canonical
> finite normalization without first reconstructing a Newton boundary.

## 6. Reproduction

```bash
.venv/bin/python plane-jc/cas/audit_f2_75_125_boundary_handoff.py
```

Intentional regeneration uses `--refresh`. The replay verifies all four
root strata, the local scale ambiguity, the conditional degree-26
signatures, the finite-flat packet budgets, and the typed refusal to infer
missing target data.
