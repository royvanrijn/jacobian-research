# A rational bisection independent in both fibrations

## Result and limits

**New deduction; verified application.** One of the three already constructed
generic conics, stored orbit47755, has degree two over **both** elliptic
parameter lines. Its moving point gives

\[
\boxed{\operatorname{rank} E_{T(u)}(\mathbf Q(u))\ge18,
\qquad \operatorname{rank} J_{Z(u)}(\mathbf Q(u))\ge13.}
\]

Here `E/Q(t)` is the original MW17 family, `J/Q(z)` the alternate MW12
family on the same K3, and `u` is a rational parameter on the conic.
These are two different elliptic curves/function-field inclusions, **not**
two extra directions on one original fibre. All maps are over Q; no
number-field extension or rational-point search is required.

**Verified application.** The rule “among the three existing conics take
the unique least alternate degree, then use its existing parameter `u=0`”
was frozen before candidate evaluation. It produces one rational fibre
with an independently certified rank-at-least18 subgroup. This is not302,
not an exact-rank assertion, and not a public novelty claim. Independence
on the alternate *specialization* at `u=0` was not tested.

**New exact obstruction, retrospectively evaluated.** All three source
curves have no rational point on the fixed first302 carrier `z=z_*`.
Consequently their entire generic alternate-translation/inversion orbits
miss both marked first-seed points. The two original translates of47755
that minimize its alternate degree are covered by this same obstruction.
Other original translates, other curves, alternating moves, and other
representatives of the302 residual direction remain outside this claim.

The prospective positive302 seed problem remains open.

## 1. Audit and input separation

**Verified prior applications, reused.** The
[conic-index construction](DET1092_RATIONAL_BISECTION_INDEX_AND_ODD_DIVISOR_CONSTRUCTION_2026-09-09.md)
already constructed these three smooth rational original bisections and
proved original generic rank18. This note does not present those conics
or that original function-field rank bound as new discoveries. The
[alternate geometry](DET1092_NORM8_DEGENERATION_AND_SMOOTH_ATLAS_OBSTRUCTION_2026-09-09.md)
already established `5I2+14I1`, geometric/rational MW rank12, and full
rational Picard marking. The
[original all-prime theorem](DET1092_ALL_PRIME_DIVISION_AND_MULTISECTION_THEOREM_2026-09-09.md)
has genus threshold9; that formula cannot simply be reused for this
different singular-fibre configuration.

**New scope.** We prove the corresponding alternate genus threshold4,
compose the three maps into the alternate parameter, certify one fixed
original specialization, and close their first-carrier incidence exactly.
These original **bisections** were not among the1,308 original **sections**
in the [previous obstruction](DET1092_INTRINSIC_LOW_DEGREE_SOURCE_OBSTRUCTION_2026-09-09.md).

The [generic protocol](../../artifacts/generated-results/elliptic-curves/det1092_bifibration_conic_sources_v1/generic-protocol.json)
reads only the parent, generic norm8 pencil, three old generic conics and
code. No exceptional point or carrier label enters its selection or seed
execution. This is a new equation-only replay of a retrospectively
developed route, not a claim that the historical derivation was blind.
The [separate incidence protocol](../../artifacts/generated-results/elliptic-curves/det1092_bifibration_conic_sources_v1/incidence-protocol.json)
then reads only the old first-carrier label and unchanged twelve-prime pool.
No later302 discovery or V3 search artifact is read.

## 2. Alternate low-genus saturation theorem

**New deduction from established height, modular and covering theory.**
Let `M=J(Qbar(z))`. For a smooth geometrically integral projective curve
`C` and any nonconstant base change `C -> P1_z` with `g(C)<=3`, the
inherited group is fully saturated and there is no torsion:

\[
 nP\in M\Longrightarrow P\in M\quad(n\ge1),
 \qquad J(\overline{\mathbf Q}(C))_{\rm tors}=0.
\tag{1}
\]

Thus any genuine multisection of degree greater than1 with normalization
of genus at most3 supplies an independent new direction. “Genuine” means
the normalization maps birationally onto its image on the elliptic surface;
a redundant parametrization of an old section does not qualify.

### Proof: isogenies and monodromy

**Established literature, applied here.** The height formula and `I2`
corrections are given in
[Schuett–Shioda, *Elliptic Surfaces*, section11.8](https://arxiv.org/pdf/0907.0298).
For every nonzero geometric section,

\[
 h(P)=4+2(P\mathbin{.}O)-\sum_v\operatorname{contr}_v(P)
      \ge4-5/2=3/2.
\]

Distinct effective sections have nonnegative intersection, each `I2`
correction is at most1/2, and `I1` contributes zero. Hence there is no
geometric torsion, in particular no rational2-torsion or degree2 isogeny.

For a prime `ell>=3`, an isogeny would induce a nonconstant map
`P1_z -> X0(ell)` factoring `j`. The cusp of width `ell` forces a pole
of `j` of order at least `ell`, whereas all poles here have order1 or2.
This is impossible. The cusp-width/moduli input is standard; see
[Milne, *Modular Functions and Modular Forms*](https://www.jmilne.org/math/CourseNotes/MF.pdf).
The argument uses surjectivity of a nonconstant map of proper curves,
not a bounded modular-polynomial test.

Therefore mod-ell monodromy is irreducible. An `I1` fibre supplies a
transvection. Conjugate it to move its fixed line; in the two fixed-line
basis these give upper and lower unipotents. Over the prime field their
powers generate all such unipotents, hence `SL2(F_ell)`. The Weil pairing
already places geometric monodromy inside that group.

### Proof: connected division covers and genus lower bounds

For `Z in M \ ell M`, the affine group of its division torsor projects
onto `SL2(F_ell)` with translation kernel either zero or all `F_ell^2`.
For odd `ell`, a zero kernel is the graph of a cocycle `c`. Central `-I`
gives `2c(g)=(I-g)c(-I)`, so the entire affine group fixes `c(-I)/2`.
That would be a generic division point, a contradiction. For `ell=2`,
a nontransitive four-point torsor has an orbit of size1 or2. Size1 gives
a generic half; a two-point orbit `{P,P'}` gives nonzero rational2-torsion
`P+P'-Z=P'-P`. Both are impossible.

Thus the primitive division cover is connected of degree `ell^2`.
The nonzero torsion cover is connected of degree `ell^2-1` by the
transitive linear action on nonzero vectors.

For `ell=2`, inertia at each of the14 `I1` fibres acts nontrivially on
the four division points, contributing at least1 to ramification.
Ignore the nonnegative contributions at the five `I2` fibres. By
[Riemann–Hurwitz](https://stacks.math.columbia.edu/tag/0C1B),

\[
g(T_{2,Z})\ge1-4+14/2=4.
\tag{2}
\]

For nonzero2-torsion, the three-sheet action has one transposition at
each `I1`; `I2` has trivial linear action. There is no ramification over
smooth fibres, so its genus is exactly5.

For odd `ell`, all19 singular fibres have nontrivial transvection linear
part (the coefficients1 and2 are nonzero modulo `ell`). An affine lift
`v -> Av+b`, with `A=I+N`, `N^2=0`, has order `ell`: both `ell*b` and
`ell(ell-1)Nb/2` vanish. It fixes at most `ell` vectors. Hence its
permutation defect is at least `(ell-1)^2`. This gives

\[
g(T_{\ell,Z})\ge1-\ell^2+\frac{19}{2}(\ell-1)^2
 =\frac{(\ell-1)(17\ell-21)}2\ge30\quad(\ell\ge3).
\tag{3}
\]

The last inequality follows from subtracting30 to get
`(ell-3)(17ell+13)/2`. The odd nonzero-torsion bound is one greater.
These are **lower bounds**, not claimed exact primitive-division genera.
In particular we do not assume affine inertia at `I2` always has a fixed point.

No curve of genus at most3 can map nonconstantly to any of these curves.
This excludes new prime-order torsion and hence all torsion. If `nP in M`,
take the least positive such `n`. For a prime `ell|n`, set
`P'=(n/ell)P`, `Z=ell P'`. If `Z in ell M`, absence of torsion forces
`P' in M`, contradicting minimality. Otherwise `P'` supplies a map from
`C` to the forbidden primitive division cover. This proves(1).

All arguments are geometric over `Qbar`; their obstruction is unchanged
by extending the constant number field. Rational rank12 on the original
alternate base, however, uses the already verified rational Picard marking.

## 3. Explicit maps and a rational source

**Verified application.** In one-based original MW coordinates let
`w=e15-e16` and `D=2O+4F+phi(w)` be the alternate fibre. An original
norm10 conic `C_a=2O+4F+phi(a)` satisfies

\[
 C_a\mathbin{.}F=2,\qquad C_a\mathbin{.}D=8-\langle a,w\rangle.
\tag{4}
\]

| Existing conic | Original word `a` | Degree over `t` | Degree over `z` |
|---|---|---:|---:|
|8044|`-e2+e10-e14`|2|6|
|47755|`e2+e12-e16`|2|2|
|103186|`e13-e16-e17`|2|5|

The packet stores exact rational-function coefficients for
[8044](../../artifacts/generated-results/elliptic-curves/det1092_bifibration_conic_sources_v1/map-8044.json),
[47755](../../artifacts/generated-results/elliptic-curves/det1092_bifibration_conic_sources_v1/map-47755.json),
and [103186](../../artifacts/generated-results/elliptic-curves/det1092_bifibration_conic_sources_v1/map-103186.json).
Each coefficient list is low-degree first; a rational function is its
`numerator` list divided by its `denominator` list. This specifies the
equations exactly without printing thousands of integer digits here.

To reproduce the construction from the old conic data, take its
parametrization `t=T(u)`, `W=W(u)`, `W^2=q(T)`, and maps

\[
x=x_0(T)+x_1(T)W,\qquad y=y_0(T)+y_1(T)W.
\]

These are **original Weierstrass** coordinates. Convert to the short model
by `X=x+b2(T)/12`, `Y=y+(a1(T)x+a3(T))/2`. From the fixed norm8 pencil
take `h,shift,cx=nx/h^2,cy=ny/h^3`, evaluated at `T`, and set

\[
 m=\frac{Y+c_y}{X-c_x},\qquad
 Z(u)=\frac{m+\mathrm{shift}/h}{h},\qquad
 V(u)=\frac{2X+c_x-m^2}{h}.
\tag{5}
\]

The checker verifies `V^2=F_Z(T)` for the alternate quartic and its
original elliptic equation exactly. The alternate origin is the old
degree-one section; the point `(T,V)` therefore represents a class on
`J_Z`. An explicit pointed-quartic/Weierstrass conversion is already in
[the two-fibration construction](DET1092_TWO_FIBRATION_SEED_CONSTRUCTION_2026-09-08.md).
No new Jacobian descent is being asserted.

The common fibre equations for `(T,W)` have gcd of degree1; with
`x1 != 0`, this proves the map is birational onto the actual conic.
Thus its normalization really is `P1_u`, with rational points furnished
by every `u in P1(Q)`. The three alternate degrees in the table are
verified after exact rational-function cancellation, independently of(4).
The original all-prime theorem and(1) prove generic independence in both
fibrations for **all three** sources.

**Explicit affine exclusions.** For any displayed affine evaluation,
exclude zeros of the stored denominators of `T,W,X,Y,Z`, of the
denominator of `V` computed in(5), and the
numerators of the original and alternate discriminants evaluated at
`T(u),Z(u)`. When using the factored formula(5), also exclude zeros and
poles of its intermediate denominators `h(T)` and `X-cx(T)`; cancelled
values may instead be evaluated from the reduced stored functions.
These specified rational functions give a finite, exact, conservative
chart-exclusion set. Proper maps may extend over some excluded chart
points. No claim is made that this chart set lists all specializations
where independence could fail. The fixed original `u=0` is checked directly.

## 4. Global minimum within the original translation orbit of47755

**New deduction, not a sampled CVP score.** For an original translation
by `q in M17`, the conic's alternate degree is

\[
D\mathbin{.}t_q(C_a)=\tfrac12\|2q+a-w\|^2-1.
\tag{6}
\]

This follows by translating its trace `a -> a+2q` and retaining
self-intersection `-2` and original degree2. For47755,
`||a-w||^2=10+8-2*6=6`. Every vector in this parity coset has norm
`2 mod4`; the rootless even original lattice has minimum4, so norm2
is impossible. The coset minimum is therefore exactly6, and(6) is at
least2 for **every** original translation.

There are only two minimum vectors, `a-w` and `w-a`. Indeed, if distinct
nonopposite norm6 vectors had the same parity, their half-sum and
half-difference would be two nonzero lattice vectors. Minimum4 would give
`||d+e||^2+||d-e||^2 >=32`, contradicting the parallelogram value24.
Thus the only minimizing translations are `q=0` and `q=w-a`.

The involution `P -> P_a-P` preserves the conic, since its two points on
each original fibre sum to `P_a`. Hence `t_(w-a)(C_a)` equals its image
under `P -> P_w-P`. This last involution is the deck involution of the
alternate quartic, hence inversion followed by a **generic alternate**
translation. It preserves `z`. Both minimum-degree conics consequently
have the same first-carrier incidence obstruction. This is a minimum
within this particular original translation orbit, not a global ranking
of all40,917 conics.

## 5. Exact first-carrier obstruction

**Verified retrospective application.** For each source write `Z=N/D`
in lowest terms. The polynomial `N-z_*D` has the full degree6,2,5,
respectively, so there is no preimage at infinity. The independently
replayed rational-root certificates exclude every finite rational root.
This is exact over Q, not a height-limited negative search.

The certificates use the unchanged old prime pool. Simple residue roots
are lifted to modulus `M>8H^2`, where the rational-root theorem bounds
numerator and denominator by `H`. Each supplied two-dimensional lattice
has determinant `M` and a verified Gauss-reduced basis. Either its shortest
vector is too long or the only possible bounded ratio is not a polynomial
root. The checker verifies the lifts, complete residue-root list, determinant,
reduction inequalities and final exclusions without importing the producer.

Generic alternate translations and inversion extend to automorphisms of
the smooth elliptic K3. They preserve `z` and rationality. A rational
point on any translated source at `z_*` would pull back to a rational
point on the original source there, contradicting the certificate.
This excludes both marked first-seed points without reading their coordinates.

It does **not** exclude any seed on the original302 fibre: that fibre can
meet a translated source at a different `z`, and its elliptic translates
can represent the same original nongeneric direction in other ways.

## 6. Fixed specialization and independent replay

**Verified application.** The selected source is47755, parameter `u=0`,
with no retries. Its original parameter and exact rational point are in
[the seed certificate](../../artifacts/generated-results/elliptic-curves/det1092_bifibration_conic_sources_v1/seed-result.json).
The generic-only footprint was checkpointed before candidate evaluation.
It uses the earliest complete good-prime prefix

`7,11,19,31,37,59,67,71,73,79,89,97,101,103,107,109`.

The independent checker rebuilds every `E(Fp)/2E(Fp)` using Sage point
groups, not the producer's manual finite-group implementation. The17
generic columns have rank17 over F2, and adding the constructed point
raises this to18. At7 the finite group has odd order13, excluding rational
2-torsion. Therefore these18 rational points are Z-independent: any
nonzero integral relation can be divided by the common powers of2 using
absence of2-torsion, giving a nonzero parity vector that contradicts
the injective mod2 image.
The stored independent quotient matrix makes this certificate inspectable.

Replay, with Sage10.9 and a25-second hard limit:

```bash
/home/royvanrijn/.local/bin/sage -python research/elliptic-curves/cas/verify_det1092_bifibration_conic_sources.sage
```

[Independent result](../../artifacts/generated-results/elliptic-curves/det1092_bifibration_conic_sources_v1/independent-replay.json):
`PASS_INDEPENDENT_BIFIBRATION_CONICS_AND_FIXED_M18_SEED`.
The all-prime proof is the written argument in section2; small-prime cycle
regressions check code, not the universal theorem. This is not formal verification.
An initial checker-only polynomial-coercion error is preserved with its
failed protocol and source; it produced no success certificate. The corrected
checker uses coefficient lists to distinguish the variable `v` from the
coefficient field `Q(u)`. No mathematical input or producer result was changed.

## Interpretation and next obstruction

**New deduction.** Generic independence is not the difficulty here: even
a rational curve of degree two in both geometries escapes both inherited
spans. Such a curve gives an equation-only source of seeds, but its rational
parameter image need not contain a prescribed fibre or carrier. Changing
the active Mordell–Weil subgroup cannot repair a proven rational-incidence
failure for a fixed source and a parameter-preserving orbit.

**Open/conjectural direction, not a result.** To reach the first302 seed by
this construction strategy requires a different rational source orbit, a
non-minimizing original translate, or a different representative/carrier.
No criterion selecting one of these from generic equations is proved here.
The long302 amplification chain remains an experimentally certified
visibility phenomenon; this note neither explains its length nor uses it
to infer seed incidence on the eight null fibres.

No point search, parameter sweep, complete atlas, class/unit-group job,
production mutation, or detached process was started. All computations ended.
