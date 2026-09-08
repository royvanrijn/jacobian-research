# A seed-specific genus-one covering and an obstruction to its continuation

## Exact result, with calibration exposed

**Verified application and new deduction.** The historical first302
elliptic Kummer class now has an explicit genus-one2-covering whose
**every rational point maps outside the generic MW17 rational span**.
This fixes the inherited-point ambiguity of the RR genus-two curves at
the level of this particular certified elliptic class.

The simplest norm-square continuation of that class to the parent is
also explicit. On all eight unchanged controls it fails a necessary local
condition: it is **not in the2-Selmer group**, so the corresponding cover
has no rational points. No integer factorization or class-group computation
is needed. The obstruction persists after multiplication by **any inherited
Kummer class**. Generically, the same lift ramifies at12 geometric places
where the elliptic family is smooth.

**Calibration limit.** This continuation literally retains the known seed
abscissa. It is not an oracle-free selector, and its302 solubility is not
a discovery or evidence that it predicts the first seed. The deliverable is
a seed-specific covering and an exact obstruction to this particular naive
continuation. Later points and V3 search artifacts are not inputs.

The [frozen protocol](../../artifacts/generated-results/elliptic-curves/det1092_seed_norm_lift_v1/protocol.json),
[expanded equations and maps](../../artifacts/generated-results/elliptic-curves/det1092_seed_norm_lift_v1/construction.json),
and [independent replay](../../artifacts/generated-results/elliptic-curves/det1092_seed_norm_lift_v1/replay.json)
retain all nine decisions and their exact witnesses.

**New verified diagnostic.** The [finite-place gate](#equation-only-finite-place-gate)
below removes the seed from construction of a *recognition test*: five local
characters separate the first seed from the generic rational span. This is
not yet a seed-producing construction. All eight controls have independently
replayed generic-only comparisons, including an exact blind spot.

**New verified obstruction.** The [polynomial-lift theorem](#polynomial-lift-obstruction-an-entire-family-not-one-copied-abscissa)
extends the ramification issue to every polynomial abscissa through degree24:
if its linear-norm lift is unramified at finite good parameter fibres,
it must already be generic. This is a construction-family obstruction,
not an exclusion of isolated seed specializations.

## What is new relative to the existing notes?

**Verified applications already established.** The
[marked transport](DET1092_MARKED_KUMMER_TRANSPORT_2026-09-08.md) identifies
the first elliptic class and obstructs the proposed ordinary RR-Jacobian
transfers. The [RR specificity controls](DET1092_RR_GENERIC_POINT_SPECIFICITY_2026-09-08.md)
show that a rational Jacobian class outside the Picard image can still come
from a generic elliptic point. The
[single-seed RR cover](DET1092_SINGLE_SEED_COVER_2026-09-08.md) supplies an
elliptic rank18 base change over a calibrated genus-two function field; its
whole-pencil gate identifies the genus-one chord curve with the original
elliptic fibre.

The present cover is different: it selects **one elliptic2-Kummer class**.
It maps with degree4 to a fixed elliptic fibre, not degree2 to the family
parameter line. Neither its curve nor its class is identified with the
genus-two RR Jacobian.

## The explicit class family

**Established literature.** For an odd cubic model `Y^2=f(X)`, write
`A=K[theta]/f(theta)`. The norm kernel in `A*/A*2` identifies with
`H^1(K,E[2])`; the rational Kummer map sends a nonbranch point to
`X-theta`. Norm-square representatives give explicit2-coverings, while
Selmer membership additionally requires local solubility everywhere. See
[Cremona–Fisher–O'Neil–Simon–Stoll, section3](https://www.mathe2.uni-bayreuth.de/stoll/papers/Explicit-Descent-III.pdf).

**Verified application.** Use the earlier exact coordinate change
`X=4x, Y=8y+4x+4`, so the determinant1092 parent is

\[
E_t:Y^2=f_t(X)=X^3+aX^2+b(t)X+c(t),
\]
\[
a=5,\qquad b=16a_4+8,\qquad c=64a_6+16.
\]

The calibrated constant is openly recorded:

\[
k=X_*=
\frac{3192654717569013009686535431767390554676}{1234321}.
\]

Set

\[
D(t)=f_t(k),\qquad
\boxed{\alpha_t=D(t)(k-\theta_t).}
\]

The determinant of multiplication by `alpha` proves

\[
\operatorname{Norm}(\alpha_t)=D(t)^4=(D(t)^2)^2.
\]

At302, `D(0)=Y_*^2`, so `[alpha_0]=[k-theta_0]` is exactly the known
seed class. Norm correction makes a cohomology class for every smooth
specialization with `D!=0`; it does **not** automatically make a Selmer class.

## Explicit genus-one model and map

**Verified construction.** For `gamma=z0+z1*theta+z2*theta^2`, expand

\[
(k-\theta)\gamma^2=R_0(z)+R_1(z)\theta+R_2(z)\theta^2.
\]

The covering is the intersection of two quadrics

\[
\mathcal T_t:\quad R_2(z)=0,\qquad w^2=-D(t)R_1(z).
\]

The first conic has the generic rational point `(1,0,0)`. A completely
explicit parametrization, in the affine coordinate `r`, is

\[
\begin{aligned}
z_0={}&(k+a)r^2+(2b-2a(k+a))r+(k+a)(a^2-b)-ab+c,\\
z_1={}&2r(r-k-a),\\
z_2={}&2(r-k-a).
\end{aligned}
\]

Substitution yields the binary-quartic model

\[
\boxed{w^2=Q_t(r):=-D(t)R_1(z(r)).}
\]

All five quartic coefficients are expanded in the certificate. Let
`N(r)=Norm(gamma(z(r)))`; the rational point map is

\[
\boxed{X=\frac{D R_0(z(r))}{w^2},\qquad
Y=\frac{D^2N(r)}{w^3}.}
\]

The exact norm identity verifies the elliptic equation. Cubic reduction
also gives

\[
X-\theta=\alpha_t\left(\frac{\gamma}{w}\right)^2.
\]

Thus every image has the specified Kummer class, not merely the inherited
trace class of an unmarked RR pair. The maps extend between the smooth
projective curves; the displayed affine formulas require `w!=0`.

**Verified invariants and new deduction.** With `c4,c6,Delta` in the
displayed cubic elliptic model,

\[
I(Q_t)=D^4c_4(E_t),\quad J(Q_t)=2D^6c_6(E_t),\quad
\operatorname{disc}(Q_t)=256D^{12}\Delta(E_t).
\]

Hence the cover is smooth of genus1 when `D*Delta!=0`, with Jacobian
Q-isomorphic to `E_t`. The rational function `X=-R0/R1` has degree4 in `r`,
with coprime numerator and denominator; since both curves have degree-two
maps to their respective coordinate lines, the covering map has degree4.
It is the2-covering of the displayed Kummer class, not a degree-four
base change of the original parameter line. No constant-field extension
is required to define the class, curve or maps.

## Why all its302 rational points are seeds

**Verified application.** At `r=k+a` the conic coordinates are `(D,0,0)`.
For `D(0)=Y_*^2`, the quartic point
`(r,w)=(k+a,D(0)*sqrt(D(0)))` maps to `(k,sqrt(D(0)))` on302.
A fresh finite-quotient calculation on the17 generic points plus this image
has binary rank18 and excludes rational2-torsion modulo31.

**New deduction.** A rationally soluble2-covering of class `delta(P_*)`
has image exactly

\[
P_*+2E_0(\mathbf Q).
\]

Let `M` be the generic subgroup and `M_sat` its rational-span saturation
inside `E_0(Q)`. The certified mod2 injection of its17 generators implies
that `M_sat/M` has odd index on the free part, so
`delta(M_sat)=delta(M)`. Equivalently, repeated division of any even
relation reduces it to an odd relation, using the mod2 injection and
absence of rational2-torsion. The retained separating character is zero
on `delta(M)` and nonzero on the seed class. Therefore

\[
\boxed{\pi(\mathcal T_0(\mathbf Q))
\cap(M\otimes\mathbf Q)=\varnothing.}
\]

This is the useful distinction from the RR curve with inherited points.
The seed class is represented by a rational point, belongs to Selmer,
and has zero Sha image. It is not a detected nonzero Sha class.
As an abstract pointed genus-one curve, `T_0` is Q-isomorphic to `E_0`;
the non-generic information is in its **covering map and Kummer class**.

## A factor-free local obstruction

**Established local fact, with an elementary proof.** At an odd prime
where the monic cubic is integral with unit discriminant, every local
rational Kummer image has even valuation in every factor of its unramified
cubic algebra. If `X` is nonintegral, the monic equation makes `v(X)`
even, and this is the valuation of every `X-theta_i`. If `X` is integral,
at most one factor has positive valuation: its reduction is the unique
simple root equal to `X mod p`, hence a linear factor. Its valuation is
even because the norm is `Y^2`. The identity and rational2-torsion cases
also have unramified images. This agrees with the general local Kummer
description in [Cremona et al., section2.2](https://www.mathe2.uni-bayreuth.de/stoll/papers/Explicit-Descent-III.pdf).

**New deduction.** Suppose in addition `k` is integral and `v_p(D)` is
odd. The factor of `k-theta` above its simple residual root has valuation
`v_p(D)`; the other factors have valuation0. Multiplication by `D` makes
the first valuation even but the remaining degree-two part odd. Consequently

\[
v_p(D)\text{ odd at such a good prime}
\Longrightarrow\alpha_t\notin\delta(E_t(\mathbf Q_p))
\Longrightarrow\mathcal T_t(\mathbf Q_p)=\varnothing.
\]

**Exact certificate without finding or factoring that prime.** Let `S`
be the product of2, the numerator and denominator of the cubic discriminant,
and the denominators of `k,D` and all cubic coefficients. Starting with
`N=abs(num(D))`, repeatedly divide by `gcd(N,S)` until the remainder `U`
is coprime to `S`. If `U` is not a square, some prime dividing `U` has
odd valuation, lies outside `S`, and therefore supplies the local obstruction.

The independent checker verifies each exact division, a Bezout identity
`aU+bS=1`, and integer inequalities `r^2<U<(r+1)^2`. Thus the existence
of an obstructing good prime is proved without identifying that prime,
factoring the large remainder, or assuming primality. This is a genuine
local obstruction, not a bounded failure to find a point.

## The frozen nine-address result

**Verified application.** On all eight unchanged controls the coprime
remainder is nonsquare. At302 it is square and the displayed rational lift
and independence certificate pass.

| Address | Bits of coprime remainder | Class conclusion |
|---|---:|---|
| scale-0131232 |1674|not Selmer; cover has no Q-points|
| scale-0257585 |1667|not Selmer; cover has no Q-points|
| scale-0487239 |1658|not Selmer; cover has no Q-points|
| scale-0177036 |1637|not Selmer; cover has no Q-points|
| scale-0043332 |1639|not Selmer; cover has no Q-points|
| scale-0590501 |1632|not Selmer; cover has no Q-points|
| scale-0290097 |1664|not Selmer; cover has no Q-points|
| scale-0748009 |1665|not Selmer; cover has no Q-points|
|302|292|rational non-generic seed class; Sha image zero|

Every member of an inherited elliptic Kummer subgroup is locally
unramified at these good primes. **New deduction:** multiplying the frozen
lift by any such inherited class or by a square cannot change its ramified
valuation pattern. In additive cohomology notation, for every control,

\[
\boxed{(\alpha_t+\delta(M_t))\cap\operatorname{Sel}_2(E_t)=\varnothing.}
\]

This conclusion covers the whole inherited coset without a parity census.
It does not exclude unrelated classes, another lift of the302 class, or
other exceptional points on these fibres. The earlier rational RR-Jacobian
classes remain rational with zero Sha image; they are different objects.

## Generic obstruction and next construction gate

**Verified application and new deduction.** The polynomial `D(t)` has
degree12, is squarefree, and is coprime to the parent's cubic discriminant.
At each of its12 geometric roots the same valuation argument ramifies
`alpha` while the elliptic fibre is smooth. The generic covering therefore
has no Q(t)-point; this norm-corrected lift is not an everywhere locally
soluble class over the geometric parameter field either. Multiplication by
generic rational Kummer classes cannot remove these good-fibre residues.

The obvious square-root construction `s^2=D(t)` has genus5, not genus0
or1. It does not supply a rational parametrization for additional addresses.
The earlier genus-two RR base cover through302 remains the better known
small-genus base change; these are not competing rank counts.

**Unresolved construction, not a conjectural result.** An oracle-free rule
must produce a different suitable class and establish its local conditions
before seeking a rational representative. Simply freezing the seed's
abscissa and multiplying its norm into the linear cubic element fails that
gate, and generic-class multiplication cannot fix it. We have not selected
such a replacement, computed full Selmer groups, proved a generic incidence
law, or obtained an independent control seed. The specific global-class
mechanism sought by the main goal remains open.

## Reproduction and limits

**Verified computation.** Constructor and independent checker each finished
in under one second, under separate25-second caps. All nine outcomes and
gcd exposures are retained. There was no point search, new address, full
integer factorization, class group, broad Selmer calculation, later-point
input, V3/V4 modification or background job.

```sh
sage -python research/elliptic-curves/cas/verify_det1092_seed_norm_lift.sage
```

The [checker](../cas/verify_det1092_seed_norm_lift.sage) does not import the
[constructor](../cas/construct_det1092_seed_norm_lift.sage). It uses direct
cubic reduction and resultants for the covering maps, independent integer
square bounds and Bezout witnesses for the local obstructions, and the
existing finite-group verifier for the302 independence certificate.

## Equation-only finite-place gate

### Exact first-seed characterization

**Verified application and new deduction.** Use the same cubic model
`E_0:Y^2=f_0(X)` and `K_0=Q[theta]/f_0(theta)` as above. For

\[
S_0=\{7,19,23,29,167\},
\]

each local point quotient `E_0(Q_p)/2E_0(Q_p)` has dimension one. Define
`epsilon_p(P)=0` when `X(P)-theta` is a square in `K_0 tensor Q_p`,
and `epsilon_p(P)=1` otherwise; set `epsilon_p(O)=0`. These are the unique
nonzero characters of those local point quotients, not arbitrary choices
of residue coordinates. Define

\[
\boxed{\chi(P)=\epsilon_7(P)+\epsilon_{19}(P)+\epsilon_{23}(P)
                   +\epsilon_{29}(P)+\epsilon_{167}(P)\pmod2.}
\]

Every one of the17 generic sections has `chi=0`. The historical first seed
has the exact character tuple

\[
(0,0,1,1,1),\qquad \boxed{\chi(P_{\rm unlock})=1.}
\]

The complete declared-place matrix has generic column rank17, and the
cubic is irreducible modulo31, excluding rational2-torsion. Thus
`H/2H -> E_0(Q)/2E_0(Q)` is injective. The rational saturation of `H`
has odd index over `H`, so every character killing `H` also kills its
rational saturation. Consequently

\[
\boxed{P\in E_0(\mathbf Q),\quad\chi(P)=1
       \ \Longrightarrow\ P\notin H\otimes\mathbf Q.}
\]

This is an exact sufficient recognition condition for a supplied rational
point. It does not assert that every point outside the rational span is
detected: for example, doubling an exceptional point kills its2-Kummer
class without making the point generic.

The17 characters at each of the five places, in the original generic
section order, are stored in the certificate. A nonzero local reference
is generic section0 at7,19,23,29 and generic section2 at167 (zero-based).
The independent square tests verify both the generic identities and the
first-seed tuple, without number-field local-character routines.

### What was frozen before the seed was evaluated?

**Verified application.** The constructor reads the parent equation and
17 generic sections, the equation-defined20 bad primes of302, and the
existing nine parameter addresses. It constructs the local images at
those same20 primes and infinity for *every* address. No control-specific
bad-prime expansion is made. Define intrinsically

\[
L_t=\prod_{v\in S}\delta_v(E_t(\mathbf Q_v)),\qquad
A_t=\operatorname{image}(H_t\longrightarrow L_t).
\]

All individual factors of `L_t` are proved to be generated by the generic
sections at all nine addresses. Local coordinates use independent generic
classes, and the constructor freezes the full annihilator `A_t^perp`.
The quotient `L_t/A_t` and its vanishing predicate are independent of
generic basis or local-coordinate changes. **These dimensions are finite
local-image dimensions, not Selmer dimensions.**

At302, `dim L_0=22`, `dim A_0=17`, giving five checks. Only after the
generic matrix is frozen is the historical first point evaluated; its
syndrome in the recorded check basis is `(1,0,1,1,1)`. The first seed
increases the joint image rank to18. The final independent replay confirms
that the earlier pre-evaluation matrix and the corrected integral-model
replay give exactly the same302 code.

**New deduction, explicitly retrospective compression.** Among the31
nonzero checks in the already frozen five-dimensional check space,16
separate the first seed. Their smallest place support is exactly five,
attained by `S_0` above. This proves that the seed matches some generic
combination at every subset of at most four places from the declared set,
but no one generic combination matches it at those five places together.
This minimum is limited to the declared places. Selecting the displayed
five-place presentation uses the seed for evaluation and is not represented
as a prospectively selected five-prime search policy.

### Unchanged controls and the exact limit of this test

**Verified application.** Each row below is equation-and-generic-section
only. A displayed generic image rank16 is a loss in this finite footprint,
not a claim that the certified MW subgroup has rank16.

| Address | `dim L_t` | `dim A_t` | Number of independent compatibility checks |
|---|---:|---:|---:|
| scale-0131232 |18|17|1|
| scale-0257585 |21|17|4|
| scale-0487239 |17|16|1|
| scale-0177036 |17|16|1|
| scale-0043332 |18|16|2|
| scale-0590501 |20|17|3|
| scale-0290097 |16|16|0|
| scale-0748009 |23|17|6|
|302|22|17|5|

**New deductions.** Single-place tests cannot separate any rational point
from the generic local image on this panel. On scale-0290097 even the
entire declared product is already generic: this particular footprint
cannot detect an exceptional point, whether or not one exists. The
scale-0748009 control has more checks than302. Thus the existence or count
of checks is not a discriminator for a rational seed. The prior bounded
null searches remain bounded nulls, not rank upper bounds.

For a general supplied parameter the equation-only code can be constructed
without an exceptional point. A nonzero syndrome of a *rational candidate*
certifies that candidate outside the rational span only after separately
certifying the generic subgroup's odd saturation index and excluding a
rational2-torsion contribution. On302 the rank17 local injection and the
modulo31 torsion exclusion supply this gate. A general rank17 row still
requires the corresponding torsion check; a rank16 row additionally needs
another saturation certificate. Zero
syndrome proves neither membership in the rational span nor nonexistence
of a seed.

### Independent finite-ring proof

**Established literature.** Odd-cubic2-descent identifies point classes
with `X-theta`; the local point quotient has dimension `m-1` at an odd
prime and `m` at2, where `m` is the number of local cubic factors.
The standard [explicit descent framework](https://www.mathe2.uni-bayreuth.de/stoll/papers/Explicit-Descent-III.pdf)
distinguishes rational Kummer classes from locally soluble covering classes
and their Sha images. The p-radical/multiplier test for p-maximal orders is
the classical Round2 criterion, illustrated in the
[algorithm documentation](https://docs.magma-maths.org/GlobalFields/NumberFieldsAndOrders/creation.html#example-round2).
That is a reference, not a computation using a licensed backend.

**Verified application and elementary deductions used by the checker.**
The constructor uses local PARI arithmetic only at the fixed primes. The
independent checker imports neither it nor its local-character helper:

1. It reconstructs every specialized equation and generic point, and checks
   the suggested cubic order by integral multiplication tables. Frobenius
   kernels give `rad(pO)`. A9-by3 linear system certifies
   `(rad(pO):rad(pO))=O`, hence p-maximality. No large discriminant is factored.
2. Each suggested square normalization is checked exactly in the cubic
   algebra. For an odd-prime unit `u`, put `F(a)=a^p` in `O/pO`. In
   `(O/pO)[z]/(z^2-u)`, Frobenius on the `z` summand is
   `M(u^((p-1)/2))*F`. Its fixed-space dimension equals the number of residue
   factors exactly when `u` is square in every factor; Hensel lifting then
   proves local squareness. Nilpotents do not alter these fixed-space counts.
3. At2 it tests all64 possible coefficient vectors modulo4 and their
   squares modulo8. For units, this is exact because `1+8O` consists of
   squares in the complete maximal order.
4. A normalized nonunit nonsquare is certified by showing that it generates
   a nontrivial radical ideal containing `pO`: its norm valuation matches
   the quotient dimension, and quotient Frobenius is bijective. Such an
   ideal cannot be generated by a square in a maximal order.
5. Exact Sturm chains verify real signs. Binary linear algebra then checks
   every local relation, every nonzero local-basis combination, the joint
   code, and all31 checks at302. This is not a subgroup-state census.

The [panel certificate](../../artifacts/generated-results/elliptic-curves/det1092_seed_local_code_v5/panel-replay.json)
and [first-seed replay](../../artifacts/generated-results/elliptic-curves/det1092_seed_local_code_v5/first-seed-replay.json)
contain the hashes and exact results. The canonical replay is:

```sh
for i in 0 1 2 3 4 5 6 7 8; do
  sage -python research/elliptic-curves/cas/verify_det1092_seed_local_code_v2.sage --case "$i"
done
sage -python research/elliptic-curves/cas/verify_det1092_seed_local_code_panel.sage
```

**Verified computational limits.** Every final independent case replay
finished in under one second, with a25-second per-call alarm. Two earlier
constructor calls hit that alarm because of unnecessarily large rational
scaling. A factor-free scaling-only correction completed both without
changing the addresses, places, generic code rule or time limit. Earlier
coercion, Sturm-method and denominator failures, sources and checkpoints
are retained in local-code versions1–4; version5 is the complete panel.
These version numbers do not refer to or modify the search pilots. No
point searches, later exceptional points, V3 artifacts, generic Selmer
dimensions, class groups or detached jobs were used.

### What is still missing for a parameter search?

**Established construction, limited to the calibrated class.** The genus-one
cover above remains the smallest-genus explicit curve currently isolating
this elliptic Kummer class. A rational curve cannot map nonconstantly to a
fixed elliptic curve, so genus one is minimal for that category of auxiliary
curves. This says nothing about minimal genus of covers of the *parameter*
line, and does not supersede the existing rational rank18 base change or
the calibrated RR genus-two construction.

**Exact next obstruction.** The new code is seed-free to construct but
does not construct a global norm-square cubic class with nonzero syndrome,
still less a rational point on its genus-one covering. Since each local
factor is full, nonzero local syndromes merely prescribe compatible-looking
local point classes; they are not automatically global Selmer classes or
rational classes. We need an equation-only source of a global class that
passes the appropriate local conditions and whose covering has a justified
rational point. The copied-abscissa continuation above fails precisely at
good-prime ramification on every control and at12 geometric good fibres
generically. Multiplying by inherited classes cannot fix it.

**Unresolved, not a claimed theorem or prediction.** No oracle-free parameter
selector or rational seed producer has been obtained. The RR-Jacobian class
`[P_unlock-P0]` is accessed here through its marked elliptic incidence point,
not through a nonexistent linear map from that Jacobian to the elliptic
curve. The seed's elliptic Kummer class is rational and has zero Sha image;
neither the compatibility quotient nor its nonzero abstract elements are
being labelled Selmer or Sha classes.

## Polynomial-lift obstruction: an entire family, not one copied abscissa

### The precise theorem

**New deduction, independently checked arithmetic.** Fix the original
parameter `t` and cubic model used above. For any `k(t) in Q[t]`, put

\[
D_k(t)=f_t(k(t)),\qquad
\alpha_k=D_k(t)(k(t)-\theta_t).
\]

If `deg k<=24` and the squareclass of `alpha_k` is unramified at every
finite **good parameter fibre**, then `D_k` is a square in `Q[t]` and
`alpha_k` is the Kummer class of a section in the full generic MW17 subgroup.
In particular, this family cannot specialize to the certified first seed
class while satisfying that generic unramifiedness condition.

**Corollary.** Odd degrees greater than four are impossible under the same
unramifiedness condition. Thus a non-generic polynomial linear-norm lift
of this kind would need even degree at least26. This is a lower bound in
the fixed polynomial-abscissa chart, not an intrinsic minimum over all
auxiliary varieties or coordinate systems.

**Scope distinction.** The theorem concerns ramification in the function
field `Q(t)` as `t` varies, not merely local solubility at rational primes
of one specialized elliptic curve. An isolated specialization of a
generically ramified family can still be a seed—the already known302
constant-abscissa example does exactly that. No elliptic fibre rank is
bounded above, and no other possible construction is excluded.

### 1. The branch polynomial has one rational orbit

**Verified application.** Write

\[
f_t(X)=X^3+aX^2+b(t)X+c(t),\qquad a=5,
\]

and let `Delta` be its monic discriminant polynomial. Its degree is24.
At the two already-used parent geometry primes149 and151, its complete
squarefree irreducible factor degrees are respectively

\[
(5,7,12),\qquad (2,2,6,14).
\]

The possible degrees of a rational factor must be subset sums of both
patterns. Their intersection is exactly `{0,24}`. Hence `Delta` is
irreducible over `Q`. The independent checker verifies every finite-field
factor by Frobenius powers and gcds and checks their product; it does not
ask a rational factorization routine to prove irreducibility.

In the field `L=Q[t]/Delta`, the cubic has one double root and one distinct
simple root. Their unique polynomial remainders of degree below24 are

\[
r\equiv\frac{9c-ab}{2(a^2-3b)}\pmod\Delta,\qquad
s\equiv-a-2r\pmod\Delta.
\]

Both remainders have degree **23**. Exact reduction verifies
`f(X)=(X-r)^2(X-s) mod Delta`, with `s-r` invertible in `L`.

### 2. Unramifiedness reduces the whole coefficient space to two cases

**New deduction.** At a finite good parameter place where `D_k` has
valuation `m>0`, the regular polynomial `k` meets one simple cubic root.
Over the unramified splitting algebra the valuations of `k-theta` are
`(m,0,0)`, so those of `alpha_k` are

\[
(2m,m,m).
\]

Thus an odd `m` ramifies the squareclass. Since `k` has no finite poles,
unramifiedness at every finite good fibre forces every odd factor of
`D_k` to divide `Delta`. Irreducibility gives exactly

\[
\boxed{D_k=d\,h(t)^2\quad\text{or}\quad
       D_k=d\,\Delta(t)h(t)^2,\qquad d\in\mathbf Q^*.}
\]

The cubic has no polynomial root: its specialization at0 is irreducible
modulo31. Therefore `D_k` is nonzero. No coefficient search or census is
used in this reduction.

For the second case, `k mod Delta` must equal `r` or `s` in the field `L`.
This is impossible for `deg k<=22`. At degree23, the degree of `D_k` is69,
whereas `d*Delta*h^2` has even degree, again impossible. At degree24 there
are just two one-parameter boundary forms:

\[
k=z+u\Delta,\qquad z\in\{r,s\},\quad u\in\mathbf Q^*.
\]

These are whole algebraic families, not a list of sampled coefficients.

### 3. Exact obstruction for both degree24 boundary forms

**New deduction and verified polynomial identities.** For either root
remainder `z`, define

\[
A=3z+a,\qquad B=3z^2+2az+b,\qquad C=f_t(z)/\Delta.
\]

Their degrees are23,46,45. Put `w=1/u`. The required square condition becomes

\[
\boxed{Q_z(t,w)=\Delta^2+wA\Delta+w^2B+w^3C
                    \quad\text{is a polynomial square in }t.}
\]

Indeed, `D_k/(u^3 Delta)=Q_z`. The polynomial `Q_z` is monic of degree48
in `t`, so any constant-times-square representation normalizes to an
actual monic polynomial square over `Q`. Its total degree in `(t,w)` is
at most48.

There is a unique monic degree24 polynomial `H_z(t,w)` whose square
matches the leading25 coefficients in `t`. Its coefficient of `t^i`
has `w`-degree at most `24-i`; this follows directly from the triangular
square-root recursion, which divides only by2. Consequently

\[
R_z=Q_z-H_z^2\quad\text{has }\deg_t R_z\le23,
\qquad \deg_w [t^{23}]R_z\le25.
\]

The residual is divisible over `Q` by `w^2` for the double-root template.
For the simple-root template it is divisible by `w^3`, because

\[
B-A^2/4\equiv0\pmod\Delta,
\]

so the first three terms of the formal square root are polynomials:

\[
\Delta+\frac{wA}{2}
  +w^2\frac{B-A^2/4}{2\Delta}.
\]

Let `e=2` or3 respectively, and write `R_z/w^e=sum_i P_i(w)t^i`.
At each of149 and151 the certificate provides an explicit identity

\[
\boxed{\sum_{i=0}^{23} U_i(w)P_i(w)=1\quad\text{in }\mathbf F_p[w].}
\]

It also verifies that `P_23` attains its exact characteristic-zero degree
bound:23 for the double-root template,22 for the simple-root template.
All rational input coefficients are p-integral, and the recursion divides
only by2. Hence `P_23` loses no degree modulo `p`.

**Why this proves a characteristic-zero obstruction.** A nonconstant
common factor over `Q` would divide the degree-preserving `P_23`. By
Gauss's lemma over `Z_(p)`, that factor has a nonconstant reduction which
divides all reduced `P_i`, contradicting the displayed Bezout identity.
Thus the rational polynomials have no common factor and no common
algebraic root. Neither boundary template has a nonzero square parameter
`w` in characteristic zero. The solution `w=0` is the infinite-`u` boundary,
not a finite degree24 abscissa. The `u=0` case already belongs to the
excluded degree23 case.

This argument uses four small finite-field identity certificates, not a
scan over values of `u` or `w` and not a heuristic modular miss.

### 4. The remaining constant twist is inherited

**Verified application of an already established theorem.** The remaining
case `D_k=d*h^2` gives the geometric section

\[
(k(t),\sqrt d\,h(t))\in E(\overline{\mathbf Q}(t)).
\]

The [completed parent proof](CURVE302_PARENT_SEARCH_AND_GEOMETRY_2026-09-07.md)
identifies the full torsion-free geometric MW group with the displayed
17-section rational basis. Every geometric section is therefore rational.
Since `h!=0`, this forces `d` to be a rational square, and the section is
in the generic MW17 subgroup. Its Kummer class is `alpha_k` because `D_k`
is a square. No generic rank, Selmer dimension or constant-field calculation
is redone here.

For any odd `deg k=n>4`, the leading term `k^3` gives `deg D_k=3n`, odd,
while both permitted square forms have even degree. This proves the
corollary excluding all such odd degrees, including25. Degree26 is merely
the first polynomial degree not excluded by this proof—not an existence
claim or a proposal to start a larger coefficient search.

### Consequences for302, controls and the next construction

**New deduction.** Any polynomial linear-norm lift of the non-generic first
seed class with degree at most25 must ramify at some finite good parameter
fibre. Changing the copied constant abscissa to a polynomial of those
degrees cannot remove all such ramification while retaining that seed.
This upgrades the earlier one-lift obstruction to a complete bounded
construction-family obstruction.

**Verified application.** All nine existing addresses, including the eight
unchanged null controls, are nonsingular specializations of this parent.
The same theorem applies uniformly: valid specializations of the
generically unramified families just classified are inherited. This is
not a distinction between their actual ranks and302's rank; exceptional
points and generically ramified families remain possible on every address.
No control exceptional point is read or sought.

**Unresolved constructive alternative.** A prospective producer must still
justify a rational seed. This calculation leaves open rational abscissas
with finite poles, other cubic-algebra representatives, larger-degree
families, and genuine base changes. It also leaves open using a ramified
family whose specially soluble fibres can be located without knowing the
seed. None of those alternatives is asserted to work. The established
genus-one seed covering and genus-two RR base change remain valid, with
their previously stated calibration limits.

### Reproduction, provenance and limits

**Verified computation.** The
[constructor](../cas/construct_det1092_polynomial_lift_gate.sage) and
[independent checker](../cas/verify_det1092_polynomial_lift_gate.sage)
each completed in under one second under separate25-second caps. The
[certificate](../../artifacts/generated-results/elliptic-curves/det1092_polynomial_lift_gate_v1/replay.json)
includes all four Bezout gates and the nine-address applicability check.
The checker uses finite-field Frobenius/gcd tests, expanded cubic identities
and supplied polynomial-square truncations; it does not import the
constructor or factor the discriminant over `Q`.

The [preflight record](../../artifacts/generated-results/elliptic-curves/det1092_polynomial_lift_gate_v1/preflight.json)
retains the earlier equation-only prime inspection and one rational
degree24 polynomial factorization. The frozen proof uses only149 and151;
the earlier inspection is not rewritten as a two-prime prospective test.
No exceptional coordinates, later points, V3 artifact, pilot mutation,
point search, new family parameter, class group, number-field computation,
Selmer-dimension calculation or detached job was used.

```sh
sage -python research/elliptic-curves/cas/verify_det1092_polynomial_lift_gate.sage
```
