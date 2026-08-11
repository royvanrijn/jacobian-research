# Carrier Wronskian classifier for F2 `(75,125)`

## Result and claim boundary

The exact checker
[`cas/verify_f2_75_125_carrier_wronskian.py`](cas/verify_f2_75_125_carrier_wronskian.py)
uses the generic source carrier, rather than a principal terminal arm, to
extract the first target equation not removable by target shears.  The result
is a finite classification of both live cofactor strata.

In the squarefree case one necessarily has

\[
\boxed{R(v)=\frac{v^2-3v+3}{25}}. \tag{1}
\]

The generic carrier maps with `(e,f)=(1,3)` to a target divisor of local ray
`(5,36)`.  After scaling its residue coordinate, the map is

\[
\boxed{g(v)=\frac{v(v^2-3v+3)}{(v-1)^3}
             =1+\frac1{(v-1)^3}}. \tag{2}
\]

Thus the two simple roots of `R` are unramified simple points over `g=0`.
They are no longer spectator points of unknown inertia.

In the double-root case, writing

\[
R(v)=\frac{(v-\rho)^2}{25(1-\rho)^2},
\]

one necessarily has

\[
\boxed{\rho^2-3\rho+1=0}. \tag{3}
\]

Put `alpha=3(rho+1)/5`.  The carrier residue map is, up to target scaling,

\[
g(v)=\frac{v(v-\alpha)^5}{(v-1)^3(v-\rho)^3}. \tag{4}
\]

Under `s=-v/alpha`, equation (3) gives the exact identity

\[
\boxed{
g(v)=\frac{729}{125}
\frac{125s(s+1)^5}{(9s^2+15s+5)^3}.
} \tag{5}
\]

Consequently (4) has degree six and passport

\[
(5,1)\mid(3,3)\mid(3,1,1,1),
\]

and is the same Belyi map already found on the principal terminal divisor.
The two double-root parameters are conjugate over `Q(sqrt(5))`.

This is an exact necessary condition for the F2 row.  It does not construct
the missing lower Laurent bands or a global Keller map, and it does not
exclude `(75,125)`.

The pinned output is
[`../artifacts/generated-results/jc2_f2_75_125_carrier_wronskian.json`](../artifacts/generated-results/jc2_f2_75_125_carrier_wronskian.json).

The argument is now proved in arbitrary primitive common-power degree in
[`COMMON_POWER_CARRIER_WRONSKIAN.md`](COMMON_POWER_CARRIER_WRONSKIAN.md).
That theorem explains the descent `36`, makes the fixed-carrier coefficient
step an explicit linear kernel, and derives both passports here directly from
the multiplicity partitions `(2,1,1,1)` and `(2,2,1)`.
<!-- status-consumer: PCW1 94b10929118f151d -->

## 1. The generic carrier expansion

Return to the original affine source coordinates

\[
q=y,\qquad v=xy^5,\qquad x=vq^{-5}.
\]

Then

\[
dx\wedge dy=-q^{-5}\,dq\wedge dv. \tag{6}
\]

The common top polynomial is

\[
C_R=x(v-1)^2R(v)=q^{-5}c(v),
\qquad c(v)=v(v-1)^2R(v). \tag{7}
\]

The certified common-power edge therefore reads

\[
P=q^{-15}c(v)^3+O(q^{-14}),
\qquad
-Q=q^{-25}\frac95c(v)^5+O(q^{-24}). \tag{8}
\]

Here `deg(c)=5`, `R(0)!=0`, and `R(1)=1/25`.

## 2. Unimodular target monomials

Set

\[
\pi=\frac{P^3}{(-Q)^2},
\qquad
h=\frac{P^5}{(-Q)^3}. \tag{9}
\]

The exponent matrix

\[
\begin{pmatrix}3&-2\\5&-3\end{pmatrix}
\]

has determinant one.  Hence

\[
d\pi\wedge dh
=\frac{\pi h}{P(-Q)}\,dP\wedge d(-Q). \tag{10}
\]

From (8),

\[
\pi=q^5U(v)+O(q^6),
\qquad
U(v)=\frac{25}{81c(v)},
\qquad
h=\frac{125}{729}+O(q). \tag{11}
\]

The right side of (10) has `q`-order `40` as a coefficient of
`dq wedge dv`.

## 3. Removing all earlier target shears

Suppose the first surviving coefficient of `h-125/729` occurs at descent
`d<36`:

\[
h-\frac{125}{729}=q^dH_d(v)+O(q^{d+1}).
\]

The coefficient below order `40` in (10) is zero, so

\[
5UH_d'-dU'H_d=0. \tag{12}
\]

Therefore `H_d^5/U^d` is constant.  Since `U` has order `-1` at `v=0`
and `H_d` is rational, equation (12) forces `5|d`.  In that case
`H_d=lambda*U^(d/5)` and is removed by the target shear

\[
h\longmapsto h-\lambda\pi^{d/5}.
\]

Iterating removes exactly descents

\[
5,10,15,20,25,30,35. \tag{13}
\]

After these shears, the first nonzero coefficient must occur at descent
`36`, because its wedge with the leading `q^5U` term has order `40`.

## 4. The forced Wronskian equation

Write the normalized target coordinate as

\[
w=q^{36}H(v)+O(q^{37}). \tag{14}
\]

The leading coefficient of (10) is

\[
\boxed{
5UH'-36U'H=\frac{5^6}{3^{12}}c(v)^{-9}.
} \tag{15}
\]

At a finite zero of `c` of multiplicity `m`, equation (15) forces

\[
\operatorname{ord}(H)=1-8m. \tag{16}
\]

The coefficient `5-4m` is nonzero for `m=1,2`, so there is no local
resonance.  At every other finite point equation (15) forbids a pole of
`H`.

At infinity there is a homogeneous resonance at order `36`; the particular
solution begins at order `39`.  This resonance is essential.  It implies
that the forced numerator of `H` has degree zero in the squarefree case and
degree one in the double-root case.

Let `d(v)` be the squarefree radical of `c(v)`.  The only possible form is

\[
H=\frac{d(v)N(v)}{c(v)^8}, \tag{17}
\]

where `deg(N)=0` or `1`.  Substitution reduces (15) to the low-degree
polynomial identity

\[
5(dN)'-4\frac{c'}c dN=\frac{5^4}{3^8}. \tag{18}
\]

This is the finite classifier.

## 5. Squarefree solution and spectator inertia

Write

\[
R(v)=av^2+bv+\left(\frac1{25}-a-b\right).
\]

For squarefree `R`, equation (17) has constant numerator `N=n`.  The three
coefficients of (18) have the unique solution

\[
a=\frac1{25},\qquad
b=-\frac3{25},\qquad
n=-\frac{5^6}{3^9}. \tag{19}
\]

This proves (1).  The target divisor selected by the generic carrier has
valuation orders

\[
\operatorname{ord}(\pi,w)=(5,36),
\]

so its transverse index is one.  Its residue coordinate is

\[
\zeta=\frac{w^5}{\pi^{36}}.
\]

Using (17)--(19) and ignoring a nonzero target scalar gives (2).  The zero
fiber consists of `v=0` and the two simple roots of `R`, all with local
degree one.  The pole at `v=1` and the point `v=infinity` both have local
degree three.

Thus the squarefree carrier row is `(e,f)=(1,3)`.  The old conditional
one-transposition and five-transposition spectator models do not describe
these two simple points: their actual carrier residue inertia is trivial.

## 6. Double-root solution and Belyi self-similarity

For the double-root cofactor, the numerator in (17) is `N=n0+n1*v`.
The two highest coefficients of (18) form a homogeneous `2 x 2` system.
Its determinant is a nonzero unit times

\[
\rho^2-3\rho+1.
\]

The constant coefficient makes `(n0,n1)` nonzero, proving (3).  Modulo (3),
one solution is

\[
n_0=\frac{625(3-\rho)}{3^8},
\qquad
n_1=\frac{625(4\rho-11)}{3^9}. \tag{20}
\]

The root of `N` is `alpha=3(rho+1)/5`.  Substitution in
`zeta=w^5/pi^36` gives (4).  Its zeros have multiplicities `(5,1)`, its
poles have multiplicities `(3,3)`, and its remaining index-three point is
`v=infinity`.  Equation (5) then identifies it with the terminal Belyi map.

This identity is geometric, not merely a matching passport.

## 7. Target extraction and revised remaining work

The ray `(5,36)` requires twelve point blowups above the smooth target point
`h=125/729`.  Its adjacent regular rays are

\[
(1,7),\qquad(4,29),
\]

and `zeta=w^5/pi^36` has orders `-1,0,1` on the left, carrier, and right
rays.  The checker records the complete thirteen-component boundary chain
with weights

```text
(0,-2,-2,-2,-2,-2,-2,-6,-1,-2,-2,-2,-2).
```

The global F2 problem is now finite at the carrier level:

1. squarefree `R` has the single value (1), with an exact cyclic cubic
   carrier cover and two unramified simple-root points;
2. double-root `R` has the two conjugate values (3), and its carrier cover is
   the same degree-six Belyi map as the terminal packet;
3. the exact
   [`carrier log-node profile`](F2_CARRIER_LOG_NODE_PROFILE.md) resolves all
   marked carrier-local points and the common fans along the principal arms.
   Their logarithmic exponent determinants are `1`, `3`, or `5`, so none
   supplies a normalization defect;
4. the subsequent
   [`upstream extraction profile`](F2_UPSTREAM_CARRIER_EXTRACTION_PROFILE.md)
   proves that the carrier-zero ladder is unimodular and identifies the
   extraction-root cokernel `R/(W^3U^18)`, whose branchwise matching quotient
   has length `54`;
5. the subsequent outgoing-tail theorem closes the terminal continuation as
   unimodular log-etale; the affine-purity frontier then forces a new
   affine-branch component and raises the source floors to `28/49`, while
   proving that coarse purity data do not determine its target curve or raise
   the degree floors.  The remaining geometric work is to recover that curve,
   factor its pullback, compile possible cancellation centers, and prove a
   global localized-Chern/descent identity.  The independent coefficient
   route still has to impose these rows on the lower Laurent system.

<!-- status-consumer: PF2CLP1 41625dd5d3f8f898 -->

<!-- status-consumer: PF2UCE1 7f15bc756cc73fff -->

<!-- status-consumer: PF2OTT1 af25012e34020e11 -->

<!-- status-consumer: PF2APF1 192055eb737d3140 -->

Until those steps are completed, the F2 row and `(75,125)` remain
unexcluded.

The first lower-system handoff is now implemented in
[`F2_75_125_CARRIER_SPECIALIZATIONS.md`](F2_75_125_CARRIER_SPECIALIZATIONS.md).
It specializes all exposed linear maps for the two double carriers over
`Q(sqrt(5))` and proves that the rational squarefree carrier is not a point
of that descent-eight component, so it must enter through the later-defect
branch compiler.
<!-- status-consumer: PF2CS1 666da98d2d24669e -->

## Reproduction

```bash
.venv/bin/python plane-jc/cas/verify_f2_75_125_carrier_wronskian.py
.venv/bin/python plane-jc/cas/test_f2_75_125_carrier_wronskian.py
```

Intentional artifact regeneration uses `--refresh` on the first command.
