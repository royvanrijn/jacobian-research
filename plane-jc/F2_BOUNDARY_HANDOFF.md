# F2 `(75,125)` boundary handoff

> **Updated status: one exact target-boundary row, global gluing open.**
> The raw contact census is cover-level rather than a list of independent
> boundary branches.  The Kummer-orbit audit reduces it to one known principal
> F2 chain, or two copies on the nonzero double-root row.  The terminal block
> then determines a genuine target extraction ray `(5,2)`, transverse index
> `1`, residue degree `6`, branch passport
> `(5,1)|(3,3)|(3,1,1,1)`, monodromy `A_6`, and an exact global meridian
> relation.  What remains is to glue this row to the full source completion and
> classify the simple spectator orbits; the degree pair is not excluded.

The three exact replays are:

- [`cas/audit_f2_75_125_boundary_handoff.py`](cas/audit_f2_75_125_boundary_handoff.py),
  for the original cover-level factorization/contact census;
- [`cas/verify_f2_kummer_orbit_transfer.py`](cas/verify_f2_kummer_orbit_transfer.py),
  for orbit transfer and normal-form filtering;
- [`cas/verify_f2_terminal_residue_cover.py`](cas/verify_f2_terminal_residue_cover.py),
  for the target row and meridian factorization.

The mathematical refinements are documented in
[`F2_KUMMER_ORBIT_TRANSFER.md`](F2_KUMMER_ORBIT_TRANSFER.md) and
[`F2_TERMINAL_RESIDUE_COVER.md`](F2_TERMINAL_RESIDUE_COVER.md).

## 1. Common edge and cover-level contact census

The upper-band calculation gives

\[
H(t)=(1+u+u^2+u^3+u^4)^2R(u^5),\qquad u=1+t,
\]

where

\[
R(v)=av^2+bv+\left(\frac1{25}-a-b\right),\qquad a\ne0.
\]

Hence

\[
C_0(u)=(u-1)^5(u^5-1)^2R(u^5). \tag{1}
\]

Put

\[
c=\frac1{25}-a-b,
\qquad
\Delta=b^2-4ac.
\]

As a factorization of the selected `X^5=x` cover restriction, the four
algebraic rows are:

| `R` stratum | cover centers | contact partition of 25 | F2 status |
| --- | ---: | --- | --- |
| \(c\ne0,\ \Delta\ne0\) | 15 | \(7,2^4,1^{10}\) | one principal chain |
| \(c\ne0,\ \Delta=0\) | 10 | \(7,2^9\) | two copies of the same chain |
| \(c=0,\ b\ne0\) | 11 | \(7,5,2^4,1^5\) | excluded |
| \(c=b=0,\ a=1/25\) | 6 | \(10,7,2^4\) | excluded |

The last two rows are incompatible with the required order vertex
`A'_0=(1,0)`: if `R(0)=0`, the approximate root has `y`-order at least five
and its cube has `y`-order at least fifteen.

## 2. Kummer-orbit transfer

Every Laurent coefficient on band `ell` comes from `k[X^5,y]` and has the
form

\[
f_\ell(t)=t^\ell u^{k_\ell}A_\ell(u^5),
\qquad
k_\ell\equiv-\ell\pmod5. \tag{2}
\]

At a nonzero conjugate center `mu^5=rho`, put `s=u-mu` and
`z_mu=X/s`.  Since `z=(s/t)z_mu`,

\[
t^\ell u^{k_\ell}A_\ell(u^5)z^\ell
=
s^\ell u^{k_\ell}A_\ell(u^5)z_\mu^\ell. \tag{3}
\]

Thus exact coefficient orders, Newton points, edges, and vertex
nonvanishing transfer to all five natural charts.  A nonzero fiber
`u^5=rho` is one Kummer orbit, not five unrelated scale problems.

The Newton-step inequality is

\[
\frac54<t_2\le4. \tag{4}
\]

Therefore a simple cofactor root is not an additional above-bisectrix F2
continuation.  The selected squared factor gives the unique `t_2=2`
principal row.  If `R` has a nonzero double root, either squared factor can be
selected and both selections give the same terminal chain.

## 3. The exact terminal target row

In Laurent coordinates `t=Xy,z=y^-1`, the terminal block is

\[
P=t^4z^3+t^{21}z^{15},
\]

\[
-Q=tz+3t^{18}z^{13}+\frac95t^{35}z^{25}. \tag{5}
\]

The support direction is `(17,12)`, with primitive normal

\[
\nu=(12,-17).
\]

It gives pole orders

\[
\nu(P)=-3,
\qquad
\nu(Q)=-5. \tag{6}
\]

At the `Q`-dominant target-infinity corner put

\[
a=(-Q)^{-1},
\qquad
b=P/(-Q).
\]

Their source orders are `(5,2)`, so the target extraction ray is `(5,2)`.
On the regular chart adjacent to `(3,1)`,

\[
\pi=b^3/a,
\qquad
\eta=a^2/b^5.
\]

The source orders are

\[
\nu(\pi)=1,
\qquad
\nu(\eta)=0. \tag{7}
\]

Hence the extracted source-to-target row has transverse index

\[
\boxed{e=1}. \tag{8}
\]

Writing `s=X^17y^5`, its residue map is

\[
\eta^{-1}
=\frac{P^5}{(-Q)^3}
=\frac{125s(s+1)^5}{(9s^2+15s+5)^3}. \tag{9}
\]

It has degree

\[
\boxed{f=6}. \tag{10}
\]

This is actual target-side toroidal data, not a contact-to-ramification
surrogate.

## 4. Braid and meridian data

The residue map in (9) has derivative

\[
h'(s)=\frac{625(s+1)^4}{(9s^2+15s+5)^4}.
\]

Its branch passport is

\[
(5,1),\qquad(3,3),\qquad(3,1,1,1), \tag{11}
\]

above branch values `0`, `infinity`, and `125/729`.  The total different is
`10=2*6-2`.  Exhaustive branch-cycle enumeration gives monodromy

\[
\boxed{A_6}, \tag{12}
\]

and the actual global meridian relation

\[
\sigma_0\sigma_\infty\sigma_{125/729}=1. \tag{13}
\]

This completes the local braid factorization for the terminal target row.

## 5. Why the old contact surrogate is retained

The original checker deliberately promoted each cover contact multiplicity
`m_i` to an unsupported row `(e_i,f_i,s_i)=(m_i,1,1)`.  Even that aggressive
promotion survived a degree-26 finite-flat packet budget.  The test remains a
valid warning that raw contact arithmetic alone cannot exclude F2.

It is not the current branch ledger.  The zero-root rows are impossible,
nonzero five-center packets are Kummer orbits, and the selected principal
orbit now has the certified row `(e,f)=(1,6)` rather than five rows derived
from contact multiplicities.

## 6. Remaining gap

The F2 route is reopened, but only at the global level.  The immediate tasks
are:

1. attach the source ray `(12,-17)` and target ray `(5,2)` to the original
   `A^2` and target completions;
2. determine how the one-chain row sits with the simple `R` spectator orbits;
3. on `Delta=0`, decide whether the two identical `A_6` residue covers land
   on the same target boundary component or on distinct components;
4. complete the source class-group/unit ledger and target canonical pullback;
5. identify which descendants, if any, map to affine nonproperness curves;
6. then run finite-normalization and global meridian filters.

The missing object is no longer the first normal order at fifteen centers and
no thirty-layer descent is required for the selected chain.  The unresolved
problem is the global source/target gluing of one or two explicit degree-six
residue packets.

## 7. Reproduction

```bash
.venv/bin/python plane-jc/cas/audit_f2_75_125_boundary_handoff.py
.venv/bin/python plane-jc/cas/verify_f2_kummer_orbit_transfer.py
.venv/bin/python plane-jc/cas/verify_f2_terminal_residue_cover.py
```
