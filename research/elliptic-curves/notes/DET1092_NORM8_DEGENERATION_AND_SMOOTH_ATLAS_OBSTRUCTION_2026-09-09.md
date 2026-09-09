# Norm-eight degenerations and the limit of the smooth bisection atlas

## Results and boundaries

**New deduction; independently verified application.** In the existing generic
norm-eight pencil, every singular member defined over `Q` is the union of two
original generic sections. The complete configuration is

\[
 5I_2+14I_1.
\]

All five `I2` parameters are rational; none of the fourteen irreducible nodal
`I1` parameters is rational. Consequently **degenerating any rational member
of this pencil cannot supply a new rational bisection seed cover**. This is
a pencil-wide obstruction, not a bounded list of missed members.

**New written corollary of completed results.** Even all40,917 smooth rational
bisection orbits together, including every generic translation, cannot exhaust
the parent's rational rank-jump parameters. Their combined incidence is thin,
whereas rank-at-least18 incidence is non-thin. Non-thin rank jumps remain after
removing this entire smooth degree-two atlas.

Neither result says whether302 belongs to that atlas. The first seed's
prospective member selection remains **UNKNOWN**. No point, parameter or
subgroup search, new CVP, integer factorization, class group, production
mutation or later cascade input was used.

## 1. A precise census boundary

The [complete degree-two census](CURVE302_LOW_DEGREE_MULTISECTIONS_2026-09-07.md)
classifies effective irreducible divisors with **arithmetic genus zero** and
original fibre degree two. On this smooth K3 these are smooth rational
`(-2)` curves. The subsequent
[index theorem](DET1092_RATIONAL_BISECTION_INDEX_AND_ODD_DIVISOR_CONSTRUCTION_2026-09-09.md)
proves all40,917 translation orbits are parametrizable over `Q`.

A singular curve of positive arithmetic genus can also have rational
normalization. Such a curve is not covered by that smooth rational census.
In particular, a nodal member of the already constructed genus-one pencil
would have arithmetic genus one and geometric genus zero. This calculation
tests that additional route on the fixed pencil, without selecting a new one.

## 2. Exact singular-member calculation

**Verified inputs.** Retain only the saved generic pencil, parent sections
and ten vertical section words from the completed full Picard-image frame.
The pencil has fibre class

\[
 D=2O+4F+\phi(w),\qquad w=e_{15}-e_{16},\quad w^2=8,
 \quad D^2=0,\quad D.F=2,
\]

using **one-based** original section indices. Its equation is
`C_z: W²=F_z(t)`, where the explicit quartic and maps are in the
[generic pencil construction](DET1092_GENUS_ONE_FIRST_SEED_COVER_2026-09-08.md).
The inherited section `S15` is a section for this second elliptic fibration.

Write `F_z(t)=a t⁴+b t³+c t²+d t+e`. The standard quartic invariants give

\[
 I=12ae-3bd+c^2,\qquad
 J=72ace+9bcd-27ad^2-27b^2e-2c^3,
\]
\[
 A=-27I,\quad B=-27J,\quad
 \Delta_z=-16(4A^3+27B^2).
\]

**Verified application.** Polynomial gcds, not factorization, give

\[
 \boxed{\Delta_z=c_0 L(z)^2R(z),\quad
 \deg L=4,\quad\deg R=14,\quad\deg\Delta_z=22.}       \tag{1}
\]

Here `L,R` are monic and squarefree, `gcd(L,R)=1`, and `gcd(A,Delta_z)=1`.
Also `deg A=8`, `deg B=12`. Thus the four finite repeated roots give `I2`
fibres, the fourteen simple roots give `I1` fibres, and infinity gives `I2`:
after the usual K3 rescaling there, `c4` is a unit and the discriminant has
order `24-22=2`. There are no other singular fibres or nonminimal locations.
The Euler sum is `5*2+14=24`.

**Established literature applied.** Unit `c4` and discriminant order `n`
give multiplicative type `I_n`; the Shioda--Tate formula then gives geometric
MW rank `19-2-5=12` for this alternate fibration. See
[Schuett--Shioda, Elliptic Surfaces](https://arxiv.org/abs/0907.0298),
the sections on singular fibres and the Neron--Severi group. This explains
the previously certified full carrier Picard-image rank12 geometrically.
It is not a rank claim for an individual specialized carrier.

### No rational irreducible singular member

**Verified application.** Every coefficient of the monic `R` is149-integral.
Its reduction has ascending coefficients

```text
[120,96,145,72,143,129,114,89,51,79,31,7,87,59,1].
```

The checker verifies all149 values are nonzero. Equivalently, over `F149`,
`gcd(R,z^149-z)=1`. The leading coefficient is1, so there is no projective
root at infinity either. A rational root of a monic149-integral polynomial
is149-integral and would reduce to a root. Therefore

\[
 R(z)\ne0\quad\text{for every }z\in\mathbf Q.         \tag{2}
\]

No irreducibility over `Q`, degree14 number field, or Galois group is inferred
from this root-free certificate. The prime was the first successful one in
the prefrozen pool of primes from3 through211; all earlier outcomes remain.

### Every rational degeneration is inherited

**Verified application.** The ten old vertical words pair into five divisors
`S_v+S_{w-v}=D`. Representatives `v` for the pairs are

| Pair | One component word `v`; the other is `w-v` |
| --- | --- |
|1|`e9-e16`|
|2|`e6+e12-e14`|
|3|`e11-e16`|
|4|`-e1+e3+e4+e6+e10-e11-e14`|
|infinity|`0`, giving `O+S_w`|

For the four finite pairs, the rational pencil coordinate computed on each
section is constant, the two constants agree, and the four distinct constants
give exactly the roots of `L`. The saved quartics there are squares in
`Q[t]`. The two signs map identically to the two listed generic sections;
the independent checker verifies their exact elliptic coordinates and words.
Infinity is identified by the divisor equality `D=O+S_w`.

Indeed the identities `v²=w.v` and `v²+(w-v)²=8` prove the divisor sum, and
the component intersection number is2. The `I2` classification confirms two
distinct transverse intersections. All rational singular fibres are now
accounted for, including infinity, without enumerating other parity cosets.

**New obstruction.** On any smooth original fibre, all points obtained from
these components lie in the specialized generic group `H`. This remains
true after any original generic translation. None can be the certified
nongeneric302 seed, or a new direction beyond full MW17 on any control fibre.
This conclusion needs no exceptional point coordinates.

For rational `z` outside these five values, the member is smooth of genus one
and has no nonconstant rational parametrization by `P1`, by Riemann--Hurwitz.
Its rational points may still give independent seeds; that smooth-member
route is not excluded. Curves obtained by moving in `z`, such as the old
alternate-section degree20 and58 constructions, are not contained in a fixed
member and are also not excluded.

**Number-field boundary.** At a root of `R`, defined over its residue number
field rather than `Q`, the irreducible nodal member has rational normalization:
the alternate section gives a smooth rational point, so its genus-zero
normalization is `P1` over that field. This observation supplies no `Q`-seed
and no descent of the new section to `Q`; no such field or maps were computed.

## 3. Why the complete smooth degree-two atlas cannot be universal

**New written deduction.** Let `S2` be the rational original parameters lying
under a rational point of *any smooth rational bisection*. The complete
quotient and index theorem give degree-two maps

\[
 T_i:\mathbf P^1_{\mathbf Q}\longrightarrow\mathbf P^1_t,
 \qquad 1\le i\le40917,
 \qquad S_2=\bigcup_i T_i(\mathbf P^1(\mathbf Q)).     \tag{3}
\]

Translations preserve the original parameter, so infinitely many translated
curves do not enlarge (3). Each image is thin and their finite union is thin.
In the fixed original rational coordinate its height count is `O(B)` for
`H(t)≤B`: for each degree-two morphism, `H(T_i(u))` is bounded below by a
positive constant times `H(u)²`, and `P1(Q)` has `O(B)` points of height at
most a constant times `sqrt(B)`. The finite sum only changes the constant.
No useful effective constant or full set of maps is computed here.

**Established literature, already verified application.** The parent has
nonconstant `j`, reduced fibres, a distinct rational elliptic fibration and
dense rational section curves. Thus its smooth rank-at-least18 parameter set
`R18` is non-thin by
[Pasten--Salgado, Theorem1.1](https://www.math.rug.nl/algebra/uploads/Main/PastenSalgado_2024.pdf),
as audited in the [two-fibration note](DET1092_TWO_FIBRATION_SEED_CONSTRUCTION_2026-09-08.md).
Consequently `R18 \ S2` is non-thin: otherwise its union with the thin set
`R18 ∩ S2` would make `R18` thin.

This is stronger than a limit on a short list of constructed conics: even
**all smooth rational degree-two multisections**, with unrestricted generic
translations, miss a non-thin supply of rank jumps. It does not locate302
inside or outside `S2`, assert positive density of `R18`, or exclude singular
rational curves of higher arithmetic genus or smooth positive-genus covers.
It also does not limit a pointed-quartic search to these curves: a parity
label does not force the chosen RR member to be its unique `(-2)` bisection.

## 4. Evidence, limits and next obstruction

The [immutable packet](../../artifacts/generated-results/elliptic-curves/det1092_norm8_singular_members_v1/)
contains the prefrozen protocol, full discriminant coefficients, no-root
certificate with all prime outcomes, all five parameter/component pairs,
and independent replay. Coefficient arrays are in ascending order.

The independent checker rebuilds quartic invariants, verifies (1) from the
product of the four displayed linear factors, checks squarefreeness and
minimality, independently evaluates the149 residues, and reconstructs every
component from the original generic elliptic group law. It imports no
producer and performs no factor discovery. Both computation and replay
complete below one second in Sage10.9, with25-second process limits:

```sh
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/audit_det1092_norm8_singular_members.sage
timeout 25s /home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/verify_det1092_norm8_singular_members.sage
```

The atlas corollary is a written implication of completed theorems, not a
new independent replay of their census or analytic proofs. Neither result
has formal verification or external review.

The original goal remains open. Generic data now exclude every rational
degeneration in this particular seed-carrying pencil, while the full smooth
rational-bisection universe is provably incomplete as a universal mechanism.
A prospective302 seed still requires an arithmetic selector for a useful
smooth member, another curve construction, or an exact obstruction covering
those remaining routes. No subsequent rank direction is investigated here.
