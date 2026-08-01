# F2 `(75,125)` common-edge boundary handoff

> **Updated status: the raw contact census is cover-level, not the final
> branch list.**  The subsequent exact Kummer-orbit audit in
> [`F2_KUMMER_ORBIT_TRANSFER.md`](F2_KUMMER_ORBIT_TRANSFER.md) proves that
> every nonzero `u^5=rho` packet has one natural Newton polygon, excludes the
> two `R(0)=0` rows from the published F2 normal form, and shows that simple
> roots of `R` are not additional above-bisectrix F2 continuations.  The
> genuine remaining rows have one principal F2 chain when `disc(R)!=0` and
> two copies of the same chain when `disc(R)=0`.  The contact table below is
> retained as an exact factorization census and as a warning against treating
> one Kummer-cover chart as a finite-normalization ledger.

The original exact contact replay is
[`cas/audit_f2_75_125_boundary_handoff.py`](cas/audit_f2_75_125_boundary_handoff.py).
Its pinned artifact is
[`../artifacts/generated-results/jc2_f2_75_125_boundary_handoff.json`](../artifacts/generated-results/jc2_f2_75_125_boundary_handoff.json).
The orbit-transfer refinement is checked by
[`cas/verify_f2_kummer_orbit_transfer.py`](cas/verify_f2_kummer_orbit_transfer.py).

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

## 2. Exhaustive cover-level root strata

Put

\[
c=\frac1{25}-a-b,\qquad \Delta=b^2-4ac.
\]

Because `R` is quadratic and \(R(1)\ne0\), the four rows below are disjoint
and exhaustive as factorizations of the selected Kummer-cover restriction.
The last two rows are subsequently removed by the order-vertex gate
`A'_0=(1,0)`; they are shown here only to preserve the exact contact census.

| `R` stratum | cover centers | contact partition of 25 | F2 normal-form status |
| --- | ---: | --- | --- |
| \(c\ne0,\ \Delta\ne0\) | 15 | \(7,2^4,1^{10}\) | one admissible principal chain |
| \(c\ne0,\ \Delta=0\) | 10 | \(7,2^9\) | two copies of the same principal chain |
| \(c=0,\ b\ne0\) | 11 | \(7,5,2^4,1^5\) | excluded by `A'_0=(1,0)` |
| \(c=b=0,\ a=1/25\) | 6 | \(10,7,2^4\) | excluded by `A'_0=(1,0)` |

For a nonzero root \(\rho\) of `R`, the equation \(u^5=\rho\) has five
distinct solutions; its multiplicity is the multiplicity of \(\rho\).
A root of order \(m\) at \(\rho=0\) instead gives contact order \(5m\) at
\(u=0\). This proves every cover-level row without numerical root
approximation.

The Kummer transfer theorem adds the missing quotient information.  Every
nonzero fiber `u^5=rho` is one orbit, and all five centers have identical
natural Newton polygons.  At `rho=1` the complete terminal block transfers
exactly, including its bracket `X^4`.  The two simple `R` orbits in the first
row do not meet the Newton-step multiplicity bound `5/4<t_2`; a nonzero
double `R` root in the second row can be relabelled as the selected squared
factor and therefore gives the same unique principal chain.

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

This distinction remains unavoidable. In smooth local coordinates `(s,z)`,
all germs

\[
g_q(s,z)=s^m+z^q,\qquad q=1,2,3,
\]

have the same edge restriction \(g_q(s,0)=s^m\), but their primitive
equality rays are

\[
\left(\frac m{\gcd(m,q)},\frac q{\gcd(m,q)}\right).
\]

Thus a raw contact partition is insufficient.  What changes after the
Kummer audit is that the first normal data of the selected nonzero orbit are
no longer missing: every band has one Kummer character, so its exact order
transfers to all five natural charts.  The remaining target-side problem is
global gluing of one or two known principal chains and classification of the
simple spectator orbits, not reconstruction of fifteen unrelated local
fans.

## 4. Historical strongest-naive finite-normalization test

For diagnosis, the original checker deliberately made the unsupported
assignments

\[
(e_i,f_i,s_i)=(m_i,1,1)
\]

for every cover-level contact multiplicity \(m_i\). Their boundary
contribution was always 25. Since a target nonproperness curve must retain a
positive affine sheet, the smallest compatible generic degree in that
surrogate was

\[
d=25+1=26.
\]

All four coarse signatures were residue-immersion-compatible. If one
additionally—and again without proof—placed every cover center in a single
closed target fiber, its minimum packet length was 25, while the finite-flat
fiber length was 26. Hence every row received the same exact verdict:

> fiber length permits this packet; further data are required.

This remains a valid demonstration that contact arithmetic alone cannot
close F2.  It is **not** the current branch count: the zero-root rows are
normal-form-incompatible and each nonzero five-center packet is one Kummer
orbit.

## 5. Revised pivot verdict

The boundary route has not yet converted the F2 family into an exclusion,
but the obstruction is now narrower:

- the two `R(0)=0` rows are gone;
- the selected fifth-root packet is one orbit with a completely transferred
  terminal principal block;
- simple nonzero `R` roots are not additional admissible F2 continuations;
- a nonzero double root gives a second copy of the same known chain;
- the coordinate degree 125 still does not determine geometric degree `d`;
- no target curve, target fiber, residue degree, or puncture profile has yet
  been certified for the one-chain or two-chain global package.

The next honest tasks are therefore:

1. compile the global log boundary for one known F2 principal chain;
2. compile the overlap of two such chains on the double-root row;
3. identify the simple `R` orbits as affine spectators or boundary branches;
4. only then construct the finite-normalization and meridian ledgers.

This reopens the boundary route without reopening the thirty-layer sequential
descent.

## 6. Reproduction

```bash
.venv/bin/python plane-jc/cas/audit_f2_75_125_boundary_handoff.py
.venv/bin/python plane-jc/cas/verify_f2_kummer_orbit_transfer.py
```

The first checker verifies the cover-level factorization/contact census and
the deliberately unsupported coarse packet test.  The second verifies the
Kummer band transfer, the conjugate terminal bracket, the order-vertex
exclusion of `R(0)=0`, the Newton multiplicity gate, and the unique principal
endpoint filter.