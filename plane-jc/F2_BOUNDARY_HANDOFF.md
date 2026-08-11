# F2 `(75,125)` boundary handoff

> **Updated status: one exact target-boundary row, global gluing open.**
> The raw contact census is cover-level rather than a list of independent
> boundary branches.  The Kummer-orbit audit reduces it to one known principal
> F2 chain, or two copies on the nonzero double-root row.  The terminal block
> then determines a genuine target extraction ray `(5,2)`, transverse index
> `1`, residue degree `6`, branch passport
> `(5,1)|(3,3)|(3,1,1,1)`, geometric monodromy `A_6`, and an exact global
> meridian relation.  The `A_6` action is four-transitive and primitive, its
> target-fixed deck group is trivial, the transverse different is zero, and
> the residue map is uniform across the admissible cofactor strata.  Any
> global realization has geometric degree at least six; two distinct packets
> over the same target divisor force at least twelve.  The row is centered at
> target infinity, so no affine-sheet `+1` applies, and purity forces a
> separate affine ramification row.  The later carrier Wronskian calculation
> classifies the spectator ledger exactly.  The affine-purity frontier then
> proves that this row requires a new boundary component, raising the source
> floors to `28/49`; it does not determine the target curve or pullback
> factorization.  The target-curve atlas reduces the former to 24 singular
> normalization charts `(3k,5k)`, `1<=k<=24`, with a forced nonunit
> collision/critical ideal.  Its `k=1` chart is controlled by one exact
> quartic and generically has four affine nodes plus the fixed infinity cusp.
> Its puncture is transverse to `(5,2)`.  The apparent terminal
> `lambda=125/729,e=3` slot is only a formal compatibility: that neighborhood
> is already resolved, so the actual affine divisor must occur elsewhere.
> Lower-band realization and global gluing remain open, so the degree pair is
> not excluded.

The first global compiler is now implemented in
[`F2_75_125_GLOBAL_ATTACHMENT_COMPILER.md`](F2_75_125_GLOBAL_ATTACHMENT_COMPILER.md).
It tracks the translated terminal chart back to the nonmonomial original
valuation `nu(x),nu(y_old),nu(x*y_old^5-c)=(-25,5,12)`.  Six blowups extract
the carrier `(-5,1)` and six more at `v=xy_old^5=c` extract each principal
arm.  The logarithmic node-profile correction resolves each interior
attachment with two, not one, source blowups.  The subsequent carrier-profile
theorem adds the necessary common fans along each principal arm and the
carrier-local spectator/fivefold branches, strengthening the lower bounds
from `19/31` to `27/48` components.  All marked terminal, carrier-local,
aligned-arm, and spectator node log cokernels vanish.  Boundary
support also identifies the smooth endpoint cokernel as `R/(w^3)` with zero
normalization defect.  The compiler
completes these source class/unit/canonical skeletons,
the four-blowup target ledger, and the five marked neighbors per terminal
component.  The target valuation `(5,2)` is unique, so the two double-root
arms necessarily share it and force degree at least twelve.  The subsequent
[`carrier Wronskian classifier`](F2_75_125_CARRIER_WRONSKIAN_CLASSIFIER.md)
extracts `(5,36)`, forces `R=(v^2-3v+3)/25` in the squarefree case and
`rho^2-3rho+1=0` in the double case, and identifies the carrier residue maps
as a cyclic cubic and the terminal degree-six Belyi map.  Finally, the
[`upstream extraction profile`](F2_UPSTREAM_CARRIER_EXTRACTION_PROFILE.md)
finds the first forced nonzero matching class: the extraction-root cokernel is
`R/(W^3U^18)` and its branchwise quotient has length `54`.  The purity row,
global cancellation/Chern ledger, lower-band realization, and global
meridians were still missing at that stage.  The subsequent outgoing-tail
theorem closes the terminal continuation as unimodular log-etale, and the
affine-purity frontier forces one new component without determining its
target curve or pullback factorization.
<!-- status-consumer: PF2CW1 a7774b0fa736b64c -->
<!-- status-consumer: PF2GA1 57dea3062b1147fb -->
<!-- status-consumer: PF2LNP1 e4f0f231bf7494d5 -->
<!-- status-consumer: PF2CLP1 41625dd5d3f8f898 -->
<!-- status-consumer: PF2UCE1 7f15bc756cc73fff -->
<!-- status-consumer: LCBBC1 b3eb4679f781c55f -->

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
<!-- status-consumer: PF2GC1 6ba3fd9eb6a0bcdf -->

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
The natural degree-six `A_6` action is four-transitive and primitive, so this
residue cover has no nontrivial degree-two/degree-three factorization.  Its
centralizer in `S_6` is trivial, hence it has no target-fixed deck
transformation.  Moreover `e=1` contributes zero to the transverse different;
the residue-different coefficients are `(4,2,2,2)` and have total degree ten.
Formula (9) contains no `R` parameter, so this local packet is identical on
every admissible principal-chain stratum.

For the rational model, the discriminant of

\[
125s(s+1)^5-r(9s^2+15s+5)^3
\]

with respect to `s` is `5^17*r^4*(729*r-125)^2`.  Thus its arithmetic
monodromy over `Q(r)` is `S_6`, while its geometric monodromy is `A_6`; the
quadratic constant field of the Galois closure is `Q(sqrt(5))`.
After scaling by `729/125` the map is Belyi.  Its regular geometric `A_6`
closure has inertia signature `(5,3,3)` and genus `25`.

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

The F2 route is open only at the unresolved global support.  The immediate
tasks are:

1. retain the upstream extraction theorem's forced cyclic cokernel
   `R/(W^3U^18)`, and use the blowup-stable Cartier charge `27` rather than
   its model-dependent raw length `54`; the kernel-line theorem makes the
   actual cyclic contribution `deg(K_root)+27`, and contraction makes this
   `27-e_root<=27`; full tangential divisibility proves `e_root=0`, so the
   cyclic root term is exactly `27`;
2. solve or stratify the 24 purity-target collision charts, factor the chosen
   implicit equation, and locate the new component's proximity chain;
3. add any remaining global resolution centers to the source
   class-group/unit/canonical skeleton; and
4. run the localized-second-Chern, finite-normalization, and global meridian
   filters with the exact root term `27`, including every possible noncyclic
   or remaining-component cancellation.

<!-- status-consumer: LCCT1 2fc6ecea7a7c8b49 -->

<!-- status-consumer: LKGD1 8a357250b5005186 -->

<!-- status-consumer: LTKT1 32ac27318f16c20c -->

<!-- status-consumer: PF2OTT1 af25012e34020e11 -->

<!-- status-consumer: PF2APF1 192055eb737d3140 -->

<!-- status-consumer: PF2ATC1 9ab722c45c586b73 -->

<!-- status-consumer: PF2K1C1 358a6ba820e8b2f1 -->

<!-- status-consumer: PF2PPA1 b24c4d80c2f8230e -->

Thus the live global ledger has one squarefree case with one principal packet
and one double-root case with two identical packets over the unique target
component `(5,2)`.  The five local branch-cycle triples form one centralizer
orbit and do not add further monodromy cases.

Several consequences can already be entered before that gluing:

- the valuation equality over the extracted target divisor gives geometric
  degree `d>=6`;
- the two distinct double-root packets necessarily lie over that same divisor
  and give `d>=12`;
- the target pole orders `(3,5)` center this row at infinity, so the affine
  companion theorem does not strengthen these to `7` or `13`;
- the target toric nodes `h=0,infinity` have three distinct preimages in the
  source-divisor interior (`s=-1` and the two denominator roots), forcing
  three boundary-attachment points with different contributions `(4,2,2)`;
  the remaining contribution `2` lies at the source endpoint `s=infinity`
  over the smooth third branch value;
- `e=1` makes this row unramified in the normal direction, so purity requires
  a separate missing-boundary row with `e>1` over an affine nonproperness
  curve;
- the global geometric Galois group has a decomposition subgroup with
  quotient `A_6`; hence `A_6` is a section of the global group, the global
  group is nonsolvable, and its order is divisible by `360`; this section
  statement does not in general make `A_6` a composition factor; if `d=6`,
  the global group is `A_6` or `S_6`;
- a same-target isomorphism between the two identical local covers, if one
  exists, is unique because their deck group is trivial.

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
