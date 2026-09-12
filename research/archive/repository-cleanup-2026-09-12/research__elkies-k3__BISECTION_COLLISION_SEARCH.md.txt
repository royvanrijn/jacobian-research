# Rootless MW17 bisection-collision search

This note defines the exact lattice gate for the proposed generic-rank-19
search. It does **not** claim a new elliptic surface, a quadratic base change,
or generic rank 19.

> **Current-chain warning (2026-08-23).** The canonical H3 entrance is
> `E7+E8/MW2 --q6--> E8+E6/MW3 --q8--> D13/MW4`; the q8 section has height
> `24`, collision degree `10`, and uses the full normalization
> `R*h*Dy == Ny*Dx mod Nx`. The later lattice/chamber chain reaches
> rootless/MW17 through eleven further nef degree-two neighbours, but those
> eleven equations have not been executed over characteristic zero. The
> degree-46 child-section and no-`Dx` q-normalizer calculations later in this
> note are historical diagnostics only. The authoritative repair ledger is
> [`H3_Q8_REAUDIT_2026-08-22.md`](H3_Q8_REAUDIT_2026-08-22.md).

## Lattice reduction

Let the exact rootless Mordell--Weil frame be `M`, stored in
[`data/lattice/rank17_gram.txt`](data/lattice/rank17_gram.txt), and write the
Neron--Severi lattice as `U + (-M)`. Use

```text
F = (1,0,0),     O = (-1,1,0).
```

A degree-two `(-2)` class is necessarily

```text
B_w = ((w.M.w-2)/4, 2, w),       w.M.w = 2 (mod 4).
```

Translation by the section indexed by `x in M` takes `w` to `w+2x`; hence the
section-translation orbit is a coset in `M/2M`. There are no reducible-fiber
component orbits in this rootless frame.

For the section `S_x=((x.M.x-2)/2,1,x)`, exact completion of the square gives

```text
B_w.S_x = (w-2x).M.(w-2x)/4 - 5/2.
```

Thus a degree-two class is nonnegative on every section precisely when its
coset has no representative of norm 2 or 6. On any K3 realization of this
rootless fibration, every survivor is in fact an irreducible smooth rational
bisection. Riemann--Roch and `B.F=2` make `B` effective. Rootlessness makes
every vertical effective component numerically a multiple of `F`. If
`B=C+kF` had one degree-two horizontal component, then `C^2=-2-4k`, so
adjunction forces `k=0`. If it had two degree-one components, they would be
sections `S1,S2`, and `B^2=-2` would give `S1.S2+2k=1`; then
`B.S1=-1-k<0`, contradicting the filter. Thus no reducible decomposition is
possible, and adjunction makes `B` smooth rational.

## Exact completed quotient enumeration

The checker
[`scripts/enumerate_rootless_bisection_orbits.sage`](scripts/enumerate_rootless_bisection_orbits.sage)
uses an explicitly verified unimodular short basis only to enumerate the
shell. It exports all representatives back in pinned `rank17_gram`
coordinates.

It establishes the following finite result.

| quantity | exact value |
| --- | ---: |
| cosets of `M/2M` | 131072 |
| cosets with norm `2 mod 4` | 65792 |
| removed by norm-6 representatives | 26672 |
| section-nonnegative bisection translation orbits | 39120 |
| minimum norm in every surviving orbit | 10 |
| norm-10 unoriented representatives across all surviving orbits | 806238 |

Each orbit consequently has a representative `B_w=(2,2,w)` with
`B_w^2=-2`, `B_w.F=2`, `B_w.O=0`, and minimum section intersection zero. The
complete orbit table is a generated result, rather than a height-cutoff: after
translation, this is the entire NS-lattice quotient for degree-two classes.

Fibrewise inversion sends `B_w` to `B_(-w)`, but this introduces no second
orbit: `-w-w=-2w` lies in `2M`. Thus the translation quotient already removes
the automatic inverse-bisection symmetry. An eventual equal-extension bucket
containing two *distinct* exported orbit masks is consequently not this
tautological inversion duplicate; its anti-invariant sections must still be
tested for independence by the height gate below.

Reproduce it with:

```bash
sage -python elkies-k3/scripts/enumerate_rootless_bisection_orbits.sage \
  --output artifacts/generated-results/elkies-k3-rootless-bisection-orbits.json \
  --orbits-output artifacts/generated-results/elkies-k3-rootless-bisection-orbits.tsv
```

The generated JSON summary has SHA-256
`43ec39e5a53ae562441c1ddc1dd9c647dceb6673239a42f7402ed57b3ac00b7a`; its
referenced complete orbit table has SHA-256
`fa43153cdb44bf8a82d8203616f0238617deba3049ed165bc3a7b0421000940b`.

### Disjoint-pair priority frontier

For two classes `B_w,B_v`, the intersection formula is
`B_w.B_v=(w-v).M.(w-v)/2-2`. Since translations independently change `w,v`
by `2M`, the minimum intersection of two translation orbits is zero precisely
when their orbit-mask XOR has a norm-four representative. The exact norm-four
shell has 1,311 masks, and all are active on the survivor set. It gives
8,895,801 unordered pairs of distinct bisection orbits that can be made
disjoint. This graph was the finite prioritization for equation matching; the
complete equation replay below now proves that none of its pairs collide.

```bash
sage -python elkies-k3/scripts/analyze_rootless_bisection_disjoint_frontier.sage \
  --output artifacts/generated-results/elkies-k3-rootless-bisection-disjoint-frontier.json
```

<!-- status-consumer: EC-K3-BISECT-DISJOINT-FRONTIER c7ad7497253ac0b3 -->

## Equation-level entrance prerequisite

### Exact unmarked Q80 first step

There is now one characteristic-zero equation-level bridge on the alternate
Q80 lattice route. The reconstructed rational Q80 coefficient curve lies
identically on the degree-eight collision divisor of the explicit first
`q=4` pencil. Thus

```text
U=(x-T)/T^2
```

defines over `QQ(u)` a first child with fibres `I5*+I5+8I1`, hence root
lattice `D9+A4`. The exact test is:

```bash
sage elkies-k3/scripts/verify_q80_unmarked_first_q4_collision_qq.sage
```

It uses a direct rational-function identity for the collision and an exact
`u=0` witness for the open nonvanishing conditions. This creates no global
marked sections and reaches neither the later `q12,q12,q4,q6` pencils nor a
rootless equation, so it supplies no bisection cover yet.

<!-- status-consumer: EC-K3-Q80-UNMARKED-FIRST-Q4-COLLISION d18185784da1e93d -->

The first source-side degree-two neighbor is the marked class
`D=O+(-P1)-F` on the explicit H3 `E7+E8/MW2` model. Its lattice shell name is
`q6`, but its old-fiber degree is two. It would be the natural first place to
derive an equation and then transport bisections toward the rootless frame.

The resolved-chart compiler supplies the local source data:

```bash
sage -python elkies-k3/scripts/verify_elliptic_neighbor_compiler.sage
sage -python elkies-k3/scripts/compile_h3_first_q6_preflight.sage
```

The first command tests the resolved-chart linear-algebra discipline. The
second is a local preflight whose corrected marked branch and valuation atlas
are consumed by the all-edge global q=6 RR cover. The complete E8 module is
`u*<1,Q>`, where `Q=u^2*(y-y(P1))/(x-x(P1))` is the integral chord. A Kodaira
label or a Smith-form saturation is not a replacement for these maps. The
four smooth `P1.O` collision blocks are also exact: their saturated frame is
`<1,(m-y(P1)/x(P1))/h>`. The simpler `h | b` rule applies only to its
base-regular sublattice.

The first exact q6 continuation builds the ten-term common ambient
`a=A/h^2`, `b=B/h`, imposes its rank-eight collision congruence, and produces
the genus-one pencil with `E8+E6` root data:

```bash
sage -python elkies-k3/scripts/assemble_h92_q6_global_rr.sage
sage -python elkies-k3/scripts/eliminate_h92_q6_global_pencil.sage
sage -python elkies-k3/scripts/certify_h92_q6_child_jacobian.sage
```

The actual-chart audit of the raw E7 chord proves `ord_Z(m/t)=-1` generically
on `E7_5`, and the repaired all-edge cover uses `Z*m/t=unit/W`.  It proves
that `m` has nonnegative exceptional orders on all seven actual E7 components
and is the marked local frame at `-P1`; therefore E7 adds no hidden row.  The
complete resolved matrix has rank eight on ten columns, so `h0(D)=2`.  Its
Jacobian has `E8+E6`, MW rank three, and the Weyl-transported sections have
the certified target Gram matrix.  Its height determinant is `316`; together
with root determinant `3` and trivial torsion/glue index, Shioda--Tate gives
absolute Néron--Severi discriminant `948`.
The final hop also pins the three source curves that meet the new fibre once,
the exact Néron--Severi old-E7 pairing rows of the third correction, and its
degree balance `4812-4811=1`; no child component is inferred from its Kodaira
symbol.  Those rows are not asserted to be a new resolved-chart trace.
For the two low-height sections, however, the affine E7 branch is transported
through its explicit blow-up chart to the two binary-quartic infinity points
and then to exact child-Jacobian coordinates; the E7_7 component is the
opposite signed point by the resolved E7 graph certificate.

There is no immediate bisection hidden among the three marked old sections:
the q6 parameter restricts with degrees `1,1,22` to the old zero, `-P1`, and
`P1`, respectively.  In particular, the degree-two neighbour does not itself
provide a quadratic cover from these sections.  Replay this narrow exclusion
with:

```bash
sage -python elkies-k3/scripts/audit_h92_q6_pencil_marked_section_degrees.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-marked-section-degrees.json
```

It excludes only those three sections, not arbitrary q6 multisections or any
rootless bisection class.
<!-- status-consumer: EC-K3-H3-Q6 177cd6e614c8b8e0 -->

The E8 module is independently replayed with:

```bash
sage -python elkies-k3/scripts/derive_h92_q6_e8_p1_branch_module.sage
sage -python elkies-k3/scripts/derive_h92_q6_smooth_po_module.sage
```

The explicit blow-ups also rule out a tempting shortcut: the marked chord
cannot be repaired by subtracting only two base Laurent terms. The terminal
infinity exceptional divisor has valuation `(u,X,Y)=(5,8,12)`, leaving a nonzero
`X/u^2` term of order `-2`:

```bash
sage -python elkies-k3/scripts/derive_h92_q6_e8_chord_obstruction.sage
```

This is an obstruction to base-only compensation, not the full E8 module and
not a neighbour equation. It writes
[`../artifacts/generated-results/elkies-k3-h92-q6-e8-chord-obstruction.json`](../artifacts/generated-results/elkies-k3-h92-q6-e8-chord-obstruction.json).

## Exact q=8 generic-fibre ambient

The selected lattice `q=8` representative is nef in the `q=6`-child chamber,
not initially in the original H3 `E7+E8` chamber.  Its source pullback has
122 deterministic fixed-component reflections.  Removing them gives the
source-nef class whose difference from `9*O+9*(-P1)` is the explicit vertical
divisor

```text
(-11,0,2,3,4,6,5,5,6,-4,-5,-7,-10,-8,-6,-4,-2,0,0)
```

in the pinned `[U,E7,E8,P1,P2]` basis. In particular, its old generic-fibre
restriction is exactly `9(O)+9(-P1)`. The compiler records that literal
repeated-section support, its simple-component vertical correction, and the
`-11F` twist before using any generic-fibre equivalence. Put

```text
m = (y-y(P1))/(x-x(P1)).
```

Eliminating `y` makes the function field quadratic over `QQ(t)(m)`, with
`x` satisfying a monic quadratic. Hence the exact 18-dimensional generic
Riemann--Roch ambient is

```text
1,m,...,m^9, x,x*m,...,x*m^7.
```

The pole orders show directly that these functions lie in
`L(9O+9(-P1))`; their quadratic-field independence and genus-one
Riemann--Roch prove that they are the complete generic-fibre space. Replay
this preparation with:

```bash
sage -python elkies-k3/scripts/derive_h92_q8_generic_rr_ambient.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-generic-rr-ambient.json
```

This does **not** yet construct the q=8 pencil. The remaining equation-level
gate is to derive the finite vertical and resolved E7/E8 quotient conditions
using this source-nef vertical divisor, and prove that their common kernel
has dimension two.  Applying q=6 local conditions to the pre-reflection class
or simply taking a ninth power is not a valid substitute.
<!-- status-consumer: EC-K3-H3-Q8-AMBIENT 2e14dd27b9a3dd79 -->

The actual q6 all-edge cover now also gives an exact generic-component layer
for the source q8 E7 condition.  For each endpoint term
`u^i/h^k*x^a*m^b`, the resolved residual order is evaluated on all seven E7
components from their actual H92 valuations and the Cartier twist.  On the
corrected 54-term seed this produces 139 negative-order groups; 41 singleton
groups give a rank-22 exact coordinate block.  The nonsingleton groups remain
for chart-function-field residue evaluation rather than being discarded as
component labels.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_all_component_generic_conditions.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-all-component-generic-conditions.json
```

<!-- status-consumer: EC-K3-H3-Q8-ALL-COMPONENT-GENERIC 6f60b40f3b99f693 -->

For the non-singleton generic groups, the compiler now pins actual H92
component frames before any residue calculation: every E7 component has a
selected blow-up chart, reduced equation, and normal `(Z,U,Y)` weight
reproducing its transported `(t,x,y)` orders.  The `Y=0` components retain
their resolved quadratic normal branch, so the next residue step can reduce
in the genuine component function field.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_generic_component_chart_frames.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-generic-component-chart-frames.json
```

<!-- status-consumer: EC-K3-H3-Q8-GENERIC-COMPONENT-CHART-FRAMES 48fc1dbd9486f9a2 -->

At the four smooth `O.(-P1)` collisions, this same vertical-difference
description identifies the actual q8 line-bundle lattice.  The E7/E8 terms
are supported at additive fibres and `-11F` may be represented away from a
chosen collision; locally the divisor is therefore just `9O+9(-P1)`.  In the
regular coordinates `q=(m-y(P1)/x(P1))/h` and `X=h^2*x`, the permitted lattice
is exactly `1,q,...,q^9,X,Xq,...,Xq^7`.  Thus every negative `h` principal
part in that frame must vanish.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_smooth_line_bundle_lattice.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-smooth-line-bundle-lattice.json
```

<!-- status-consumer: EC-K3-H3-Q8-SMOOTH-LINE-BUNDLE-LATTICE 340ac8bee0e38750 -->

The corrected marked-E7 normalization gives a 54-term endpoint seed.  Its
1080-by-54 smooth principal-part map has column rank 54 modulo `43`; the
resulting full-column-rank minor certifies injectivity over `QQ`.  Thus this
corrected seed has no nonzero fully smooth element.  The endpoint envelope
must therefore be enlarged; this is not a no-pencil or no-rank-19 result.

```bash
sage -python elkies-k3/scripts/probe_h92_q8_smooth_principal_parts_modp.sage \
  --prime 43 --output artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43.json
```

<!-- status-consumer: EC-K3-H3-Q8-ENDPOINT-SMOOTH-OBSTRUCTION d629df5d215d009c -->

There is a controlled enlargement: raise each endpoint denominator exponent
by a common `r`, retaining the corresponding full E8/E7 endpoint interval.
For the corrected `r=4` ambient, the smooth map has rank 342 on 342 columns
modulo `43`, so it has no modular survivor.  The earlier 335-column,
three-dimensional kernel used the invalid `m/t` marked-frame normalization
and is withdrawn.  This remains a bounded obstruction, not a
characteristic-zero q8 pencil statement.

```bash
sage -python elkies-k3/scripts/probe_h92_q8_smooth_principal_parts_modp.sage \
  --prime 43 --extra-h-power 4 --include-kernel \
  --output artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra4.json
```

<!-- status-consumer: EC-K3-H3-Q8-ENLARGED-ENDPOINT-SMOOTH-KERNEL 3a89050b29ce8eac -->

The marked E7 edge accepts every generator of this corrected `r=4` ambient:
the exact inequality `i<=4k+d(a,b)` remains true for all 342 labels, with
`d=6` for `m^b` and `d=8` for `x*m^b`.  Any new E7 condition must come from
the other five edges.

```bash
sage -python elkies-k3/scripts/certify_h92_q8_enlarged_endpoint_marked_e7.sage \
  --extra-h-power 4 --extra-e7-pole 0 \
  --output artifacts/generated-results/elkies-k3-h92-q8-extra4-marked-e7-cover.json
```

<!-- status-consumer: EC-K3-H3-Q8-EXTRA4-MARKED-E7-COVER 5deb19aa922fd23b -->

E8 also adds no row to this ambient: `h(0)=1`, while each generator retains
the same certified E8 `u`-floor.  Thus only five non-marked E7 edges remain.

```bash
sage -python elkies-k3/scripts/certify_h92_q8_enlarged_endpoint_e8_cover.sage \
  --extra-h-power 4 --output artifacts/generated-results/elkies-k3-h92-q8-extra4-e8-cover.json
```

<!-- status-consumer: EC-K3-H3-Q8-EXTRA4-E8-COVER 135e8e4da32d56b6 -->

### Exact degree-two marking on the q=6 child

The binary-quartic covariant used here is already the degree-two covering map
to the child Jacobian.  If `Pmap` and `Qmap` are the transported covariant
differences of the old `E7_7` and affine-`E7` components, the primitive marked
section is therefore

```text
S = Pmap + Qmap,
MW(S) = (-2,-2,0).
```

It has height `24`, meets the standard zero in degree `10`, and has `x` and
`y` denominators equal to the square and cube of the same reduced degree-ten
smooth collision divisor.  The formerly used section is exactly `2*S`; its
height `96` and degree-46 collision divisor explain the obsolete downstream
modules archived below.  The primitive divisor is `O+S`, with exact generic
Riemann--Roch basis `<1,m>`, `m=(y+y(S))/(x-x(S))`, and the checker records the
monic quadratic relation for `x` over `QQ(T)(m)`.  The complete corrected
pencil and its `D13/MW4` child are certified by the later characteristic-zero
checker; this subsection records only the marking identity.

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_marking.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-marking.json
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-MARKING 745bf011cb47e7f3 -->

### Physical E8+E6 component target

Changing from the canonical q6 zero to the explicit transported old zero is
also essential for the vertical divisor.  Relative to the finite root lattice
orthogonal to `<F,O+F>`, both `O` and `S` have degree zero on every finite
component.  The exact q8 vertical part is

```text
V = -F + (3,5,6,4,2,3)_E6 + (4,5,7,10,8,6,4,2)_E8.
```

In the pinned simple-root orders, its component degrees are
`(-1,-1,0,0,0,0)` on `E6` and `(−1,0,0,0,0,0,0,0)` on `E8`.  These labels are
lattice-only; matching them with actual resolved II*/IV* charts and compiling
the corresponding finite quotient modules remains an equation-level task.

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_physical_root_target.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-physical-root-target.json
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-PHYSICAL-ROOT-TARGET 064318c2afe537fd -->

The Weyl-nef q8 fibre differs from this dominant D13 representative by four
root reflections.  Its IV* ideal is the arm-invariant `(u^2,X,Y)`, and its
finite q-regular module is `<(1,lift(R/Nx)),(0,f_II^2*f_IV^2)>`; the infinity
lattice remains open.

<!-- status-consumer: EC-K3-H3-Q8-CHILD-NEF-LOCAL-MODULE e2887bd2bd4f6c27 -->

### A certified lattice bisection pencil on the explicit child

The abstract Weyl-nef representative is not component-nef for the transported
old zero.  Reflecting it in the actual `E6+E8` components (102 deterministic
reflections) produces a primitive isotropic class of old-fibre degree two.
Its simple-component degrees are `(1,0,0,0,1,0)` on `E6` and
`(1,0,0,0,0,0,0,0)` on `E8`; both affine degrees are zero.  The exact
short-coset enumeration has no negative section wall, and the standard
degree-two parity identity excludes a negative bisection.  Since any curve
negative on this effective class would be a fixed component and hence have
old-fibre degree at most two, this is a lattice proof that the class is nef.
It therefore defines a genus-one pencil whose generic member is an old-fibre
bisection.

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_physical_root_target.sage \
  --representative component-nef \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-component-nef-physical-root-target.json
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_section_walls.sage \
  --target artifacts/generated-results/elkies-k3-h92-q6-child-q8-component-nef-physical-root-target.json \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-component-nef-section-walls.json
```

This supplies one exact lattice bisection pencil, not its equation, branch
divisor, quadratic extension, collision, or a rank claim.  In particular,
the available standard-Weierstrass chord is for the divisor after translating
the transported old zero to the Weierstrass infinity section.  That exact NS
translation has not yet been transported for this component-nef class, so its
existing local chord modules cannot yet be used as this bisection equation.

<!-- status-consumer: EC-K3-H3-Q8-CHILD-COMPONENT-NEF-BISECTION-PENCIL 84142196f27e5e2d -->

The generic marking is now explicit.  If `P0` is the transported old zero as
a finite point on the standard Weierstrass child and `S` is the previous
marked point, put `Q=P0+S`. Translation by `-P0` takes the physical horizontal
divisor `P0+Q` to `O_standard+S`. Thus the desired chord is the pullback of
the standard chord: with `lambda=(y+y(P0))/(x-x(P0))`, set
`x'=lambda^2-x-x(P0)`, `y'=lambda*(x-x')-y`, and use
`(y'+y(S))/(x'-x(S))`.

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_component_nef_chord.sage
```

This is an exact generic-fibre function. Its translated II*/IV* chart modules
and infinity condition still need compilation before it becomes a global
pencil or produces a branch divisor.

<!-- status-consumer: EC-K3-H3-Q8-CHILD-COMPONENT-NEF-CHORD-TRANSPORT c3896078bcdd432d -->

At level zero this generic chord gives an exact quadratic function field. In
the translated coordinates it is the line `y'=-y(S)`; its fixed intersection
at `-S` is removed from the Weierstrass line intersection, leaving a monic
quadratic in `x'`. Pulling this equation back by `tau_-P0` gives a curve of
old-fibre degree two. Its discriminant is exactly
`(-3*Nx^2-4*A*Dx^2)/Dx^2`; `Dx` is a square, while the degree-192 numerator
is squarefree. Thus this uncorrected chord level has branch degree 192, so it
is *not* a rational bisection and must not enter the collision hash.

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_component_nef_bisection_branch.sage
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-COMPONENT-NEF-GENERIC-BISECTION-BRANCH 6d0dbb4b90b14710 -->

The translation centre itself presents no additive-fibre obstruction: at both
cusps the exact chord denominator is a unit, its translated cusp image has
the prescribed additive orders, and its tangent action is unipotent with
determinant one.  The exact finite point `P0` therefore specializes regularly
with nonzero Weierstrass gradient at both the II* and IV* fibres, so
translation by `-P0` extends on both Néron smooth loci. This is a
singular-germ and tangent prerequisite only; it is not a resolved blow-up-chart
pullback and does not identify transported local quotient modules.

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_component_nef_translation_additive_jets.sage
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-COMPONENT-NEF-TRANSLATION-ADDITIVE-JETS f76374817b493c9a -->

The II* vertical part now has an exact local ideal.  In the unit-normalized
germ `Y^2=X^3+u^4*a(u)*X+u^5*b(u)`, the physical E8 cycle maps to the
valuation cycle of `X`.  The required complete ideal is therefore
`(u^2,X,Y)`, with quotient basis `1,u` and colength two.  This does not yet
trivialize `<1,m>` near II*, so it is not by itself a q8 coefficient module.

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_iistar_vertical_ideal.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-iistar-vertical-ideal.json
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-IISTAR-VERTICAL-IDEAL ea029da7275e86da -->

### IV* vertical ideal pair

The IV* vertical condition is also explicit, up to the one remaining E6 arm
orientation.  In the unit-normalized germ

```text
Y^2 = X^3 + u^3*a(u)*X + u^4*b(u),   b(0)=c^2,
```

the physical E6 cycle has two chart realizations exchanged by the E6 diagram
involution.  They give the conjugate complete ideals

```text
(Y-c*u^2, u*X, X^2, u^3),
(Y+c*u^2, u*X, X^2, u^3),
```

each of colength four with quotient basis `1,u,u^2,X`.  Selecting between
them requires the still-missing attachment of the pinned physical E6 roots to
the resolved chart arms; neither choice may yet be imposed as the q8 global
condition.

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_ivstar_vertical_ideal.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-ivstar-vertical-ideal.json
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-IVSTAR-VERTICAL-IDEAL-PAIR 4300b94b65dabc64 -->

The arm is now oriented by the transported old `E7_7` section.  Its source NS
class meets physical E6 root five, and its exact IV* jet has
`Y/u^2=c`; it therefore meets the plus outer resolved branch.  This selects
the physical q8 IV* ideal

```text
(Y+c*u^2, u*X, X^2, u^3),
```

with chart cycle `(3,2,3,4,5,6)`, degrees `(0,0,-1,0,-1,0)`, and colength
four.  This is still only the vertical IV* condition, not a generic-chord or
global-pencil construction.

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_ivstar_orientation.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-ivstar-orientation.json
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-IVSTAR-ORIENTATION 4e931cdcb76201c3 -->

### Additive chord blocks

For the actual generic chord `m=(y+y(S))/(x-x(S))`, the selected additive
ideals can now be turned into coefficient conditions directly.  Reducing in
`(u^2,X,Y)` yields a rank-two block on the II* two-jets of `a+b*m`.  Reducing
in `(Y+c*u^2,u*X,X^2,u^3)` yields a rank-four IV* block on the three-jets and
the `X` coefficient.  The IV* residue includes the nonzero `u^2` correction
arising from `Y=-c*u^2`. Both blocks now use the compiler's resolved
marked-chord quotient adapter with their explicit quotient ideals and
normal-form bases. These are local blocks only: their common global ambient
and compatibility with the degree-46 smooth module remain unknown.

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_additive_chord_blocks.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-additive-chord-blocks.json
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-ADDITIVE-CHORD-BLOCKS bf6fc5b74f51fc0f -->

### Saturated global-ansatz reconnaissance

The local blocks identify a small, stable modular candidate window.  At good
primes `43`, `53`, and `59`, impose the smooth relation
`A*D+B*N=0 mod h^2`, restrict `B` to degree at most seven, and set the free
`h^2*C` correction to zero.  The eight-dimensional coefficient space has
six independent additive conditions, leaving a two-dimensional kernel.  The
adjacent windows of degrees six and eight leave dimensions one and three.
This motivates an exact global-envelope derivation, but cannot replace one:
the bounds are not yet consequences of the q8 divisor.

```bash
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_saturated_ansatz.sage \
  --prime 43 --max-b-degree 7 --max-c-degree -1 \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-saturated-ansatz-probe.json
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-SATURATED-ANSATZ-PROBE cdcfe4ca7989839d -->

### Branch screen of the saturated candidate pencil

The two-dimensional bounded slice can also be rejected at the level of its
actual finite pencil ratios.  For the deterministic kernel basis `(f0,f1)`,
each finite constant level `f0/f1=v` is substituted into the reduced
Weierstrass equation.  The discriminant in `x` of the resulting quadratic
has branch degrees `72:1, 74:2, 76:40` at `p=43`; `72:1, 74:4, 76:48` at
`p=53`; and `72:1, 74:6, 76:52` at `p=59`.  In particular, none has the
degree-four branch divisor of a genus-one double cover.  This rules out only
the finite levels of this deliberately bounded modular slice; it is not a
global q8-pencil nonexistence result.

```bash
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_saturated_pencil_modp.sage --prime 43
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_saturated_pencil_modp.sage --prime 53
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_saturated_pencil_modp.sage --prime 59
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-SATURATED-PENCIL-OBSTRUCTION 5a4beade3e2caff1 -->

### Complete finite coefficient module (modular)

Removing the arbitrary degree window, the smooth parametrization has two
polynomial parameters, `A=A0*B+h^2*C`.  Modulo the degree-five CRT modulus
supported at II* and IV*, the complete six-row condition matrix has rank six
on ten residues.  Its finite-compatible rank-two `F_p[T]` module has Smith
degrees `(1,5)` at each of the good primes `43`, `53`, and `59`.  This is an
intrinsic finite-base profile; the infinity trivialization needed for a q8
pencil remains open.

```bash
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_finite_module_modp.sage --prime 43
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_finite_module_modp.sage --prime 53
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_finite_module_modp.sage --prime 59
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-FINITE-MODULE-MODULAR 78fc7f298da9eaf6 -->

### Complete finite coefficient module (exact q frame)

The same finite gate is now established in characteristic zero without
computing the expensive global interpolation modulo `h^2`.  In the regular
smooth frame

```text
q=(m-p)/h,       coefficient = C+B*q,
```

`q` vanishes in the II* quotient and has the exact IV* residue
`c*u^2/(x(S)(0)h)+(-y(S)(0)/(x(S)(0)^2h))*X`. The six II*/IV* rows have a
nonzero exact minor in columns `B0,C0,C1,C2,C3,C4`. More strongly, their
kernel is exactly

```text
< (f_IV*,0), (0,f_II*^2*f_IV*^3) >,
```

with Smith degrees `(1,5)` and determinant degree six. This is the exact
counterpart of the modular finite profile, but it still carries no infinity
condition.

```bash
sage elkies-k3/scripts/derive_h92_q6_child_q8_finite_q_module_qq.sage
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-FINITE-Q-MODULE-QQ 83452e0f702d6d9a -->

The q-frame has a separate exact global pole profile.  Writing
`x(S)=Nx/Dx` and `y(S)=Ny/Dy`, clearing the chord shows that its generic
vertical pole divisor is exactly `Nx/gcd(Nx,Ny)`: it has degree `96`, is
coprime to the collision divisor, IV* factor, and discriminant, while `q`
has order `44` at infinity.  Thus a global pair `C+B*q` must either make `B`
vanish against this divisor (and `f_IV`) or cancel its base Laurent principal
parts through `C`; treating q as a globally base-regular coefficient is not
valid.

```bash
sage elkies-k3/scripts/derive_h92_q6_child_q8_q_pole_profile.sage
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-Q-FRAME-POLE-PROFILE 3685b882ed19702a -->

Modulo `43` and `59`, the unique base principal-part correction
`R=Ny*(h*Dy)^(-1) mod Nx` has degree `95`; it cancels the degree-96 pole
divisor and leaves normalized infinity order one. At `53` the leading
coefficient cancels and its degree drops to `94`, so that prime is unsuitable
for reconstructing the leading normalization coefficient.

```bash
sage elkies-k3/scripts/probe_h92_q6_child_q8_q_pole_normalization_modp.sage --prime 43
sage elkies-k3/scripts/probe_h92_q6_child_q8_q_pole_normalization_modp.sage --prime 59
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-Q-FRAME-NORMALIZATION-MODULAR 6bd8f54f7e8887a9 -->

The remaining exact normalization may be computed by deterministic CRT rather
than raw rational extended Euclid.  It uses 31-bit good primes and can save a
resumable local state of the 96 coefficient residues.  That state is only a
work checkpoint: it is not an artifact or certificate until rational
reconstruction, withheld-prime validation, and the exact `QQ[T]` congruence
all succeed.

```bash
sage elkies-k3/scripts/reconstruct_h92_q6_child_q8_q_pole_normalization_crt.sage \
  --prime-bits 31 --maximum-primes 4000 --minimum-primes 5000 \
  --reconstruct-every 5000 \
  --checkpoint artifacts/local/h92-q8-q-normalizer-4000.crt.json \
  --output artifacts/local/h92-q8-q-normalizer-4000.json
```

The exact continuation has now succeeded at 4,600 good 31-bit primes.  The
committed certificate records a 138,048-bit CRT modulus, five incorporated
and three withheld modular checks, and (decisively) the exact congruence
`R*h*Dy = Ny mod Nx` over `QQ[T]`.  It clears the generic q-frame principal
part; it does not yet supply the remaining additive/smooth gluing conditions
or a q8 pencil.  Reproduce it directly with:

```bash
sage -python elkies-k3/scripts/reconstruct_h92_q6_child_q8_q_pole_normalization_crt.sage \
  --prime-bits 31 --maximum-primes 4600 --minimum-primes 4001 \
  --reconstruct-every 100 --accepted-validation-primes 5
```

The resulting artifact is
[`../artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-pole-normalization-crt.json`](../artifacts/generated-results/elkies-k3-h92-q6-child-q8-q-pole-normalization-crt.json).

<!-- status-consumer: EC-K3-H3-Q8-CHILD-Q-FRAME-NORMALIZATION-CRT 6f6a4e097d4bddd5 -->

Subtracting this correction produces the globally base-regular frame
`q_regular=q-R/Nx`.  The checker expands it with the `Nx` factor cancelled,
checks the required coprimalities, and obtains infinity order one.  This
removes the generic-base gluing obstruction, but the finite additive and
smooth q8 modules must still be transformed and assembled separately.

```bash
sage -python elkies-k3/scripts/certify_h92_q6_child_q8_q_regular_frame.sage
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-Q-REGULAR-FRAME 1cead360048f47f6 -->

The complete finite II*/IV* condition transports exactly to this frame.  In
`(B,C)` coordinates for `C+B*q_regular`, it is generated by

```text
(f_IV*, lift(f_IV*R/Nx)),  (0, f_II*^2*f_IV*^3).
```

Here the lift is modulo `f_II*^2*f_IV*^3`; `Nx` is a unit there.  The module
still has finite codimension six and Smith degrees `(1,5)`.  This now aligns
the global base coordinate and finite additive constraints, but does not yet
impose the smooth/resolved q8 conditions needed for a pencil.

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_q_regular_finite_module_qq.sage
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-Q-REGULAR-FINITE-MODULE 794f629d93420ed4 -->

The degree-46 smooth collision module imposes no further row in the same
frame.  Locally at `h`, write `p=N/(h*D)`.  The exact transition from
`C+B*q_regular` to the old saturated coefficients is
`A=h^2*(C-B*R/Nx)-B*N/D`, so
`A*D+B*N=h^2*D*(C-B*R/Nx)`.  Conversely the old smooth congruence makes `C`
regular.  Hence `q_regular` is a full local frame at every smooth collision;
only the finite additive and infinity modules remain to be assembled.

```bash
sage -python elkies-k3/scripts/certify_h92_q6_child_q8_q_regular_smooth_frame.sage
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-Q-REGULAR-SMOOTH-FRAME b394a23706a914ce -->

The ratio of these two canonical normalized generators is also not the
pencil.  At every constant level over both `GF(43)` and `GF(59)`, its
quadratic discriminant over the old base has branch degree `488`, `490`, or
`492`, never the necessary degree `4`.  The bounded monomial family
`V=(T^d*g1+T^e*M)/M`, for `0 <= d,e <= 4` and every constant level, likewise
has no branch degree below `488` at either prime.  This rejects that bounded
finite-generator/infinity window, not general combinations or a nonconstant
construction.

```bash
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_q_regular_generator_modp.sage \
  --prime 43 --max-a-monomial-degree 4 --max-b-monomial-degree 4
sage -python elkies-k3/scripts/probe_h92_q6_child_q8_q_regular_generator_modp.sage \
  --prime 59 --max-a-monomial-degree 4 --max-b-monomial-degree 4
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-Q-REGULAR-GENERATOR-OBSTRUCTION 487f98d8a254b683 -->

The diagonal finite basis alone does not solve the infinity problem. Its
most direct ratio would be

```text
V=f_IV*(m-p)/(h*f_II^2*f_IV^3).
```

At good primes `43`, `53`, and `59`, eliminating the old Weierstrass
coordinate gives branch degree `484` for `V=1` (also `V=2` at `43`), rather
than the genus-one degree `4`. Thus this unnormalized diagonal candidate is
rejected; an actual pencil must use a nontrivial infinity normalization.

```bash
sage elkies-k3/scripts/probe_h92_q6_child_q8_diagonal_candidate_modp.sage --prime 43 --v 1
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-DIAGONAL-PENCIL-OBSTRUCTION ec76cc0097b76a4c -->

### Smooth zero-collision module

The degree-46 smooth collision divisor is now a complete local condition, not
merely a list of base points.  Write `P=-S` and
`m=(y-y(P))/(x-x(P))`.  If `h` is its collision divisor and
`p=y(P)/x(P)`, the exact saturated frame is

```text
<1, (m-p)/h>.
```

Consequently a base-regular `a+b*m` is regular precisely when `h | b`.  In
the saturated form `a=A/h^2`, `b=B/h`, the finite condition is

```text
A*(den(p)/h) + B*num(p) = 0 mod h^2.
```

The `h` and `h^2` quotients have dimensions 46 and 92.  This supplies every
smooth O.S collision condition for the selected q8 marking; the II*/IV*
modules and their common global kernel remain separate.

```bash
sage -python elkies-k3/scripts/derive_h92_q6_child_q8_smooth_collision_module.sage \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-q8-smooth-collision-module.json
```

<!-- status-consumer: EC-K3-H3-Q8-CHILD-SMOOTH-MODULE 559fade7fa0cb618 -->

### Superseded IV* branch calculations

The earlier depth-two IV* entrance and plus/minus component orientation were
for the incorrectly normalized q8 section.  The corrected section is smooth
at IV*, so those calculations are historical and must not be used for the q8
module or global pencil.

<!-- status-consumer: EC-K3-H3-Q8-CHILD-IVSTAR-ENTRANCE 86bf6cc2487fd0b4 -->
<!-- status-consumer: EC-K3-H3-Q8-CHILD-IVSTAR-COMPONENT ee2d661e3dc73ef0 -->

### Resolved q=8 E7 target

The q=8 source-nef class must still be converted from the source `E7`
numbering to the actual H92 blow-up charts.  The affine-component
normalization gives the exact map

```text
source E7_i -> resolved E7_j: (1,6,4,3,7,2,5).
```

In that resolved order the q=8 restriction degrees are
`(0,1,0,0,2,0,1)`.  If `c6` denotes the exceptional cycle of the certified
q=6 marked module `Z*J_-P1^dual`, the q=8 cycle is exactly

```text
c8 = 9*c6 + (2,5,6,4,6,3,5).
```

Thus the q=8 E7 condition is not a guessed higher-order analogue: it is the
ninth tensor power of the resolved q=6 non-Cartier module followed by this
integral exceptional twist.  Its marked-branch pole order is nine; reducing
only its E7 class group class must not replace it by a single branch factor.
Moreover, the twist has exceptional intersection degrees
`(0,1,0,0,-7,0,1)`.  It is therefore not anti-nef, so it cannot be replaced
by a single complete ideal in the singular E7 local ring.  The pending q=8
quotient compiler must retain resolved-chart trivializations and their
gluing.
Replay the target calculation with:

```bash
sage -python elkies-k3/scripts/derive_h92_q8_e7_local_target.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-local-target.json
```

The finite quotient map for this module, the E8 target, and the global
two-dimensional q=8 kernel remain separate gates.
<!-- status-consumer: EC-K3-H3-Q8-E7-TARGET a6dd94428dfa14e4 -->

The exact cleared ninth-power interface is also materialized. Set
`J_P1=(x-xP1,y-yP1)`; on every actual edge, multiply a q8 representative by
`(x-xP1)^9*g` and compare it with the ten generators of `t^9*J_P1^9`.
This retains the marked branch and pole order nine rather than treating the
E7 correction as one complete ideal. Replay the manifest with:

```bash
sage -python elkies-k3/scripts/derive_h92_q8_actual_e7_power_module.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-actual-e7-power-module.json
```

It is the precise input for the remaining resolved-chart quotient evaluator,
not that evaluator or a q8 pencil.

After the source envelope was deepened to 54 terms and `h^-15` smooth jets,
the declared `r=7` enlargement has 558 columns and a two-dimensional smooth
kernel modulo each of `43,53,59,89`. The generic point of the actual
`E7_4--E7_3` chart already excludes both directions: its two independent
unique negative leading labels have residual orders `-54` and `-51`. Replay
the mod-43 certificate with:

```bash
sage -python elkies-k3/scripts/probe_h92_q8_smooth_principal_parts_modp.sage \
  --prime 43 --extra-h-power 7 --include-kernel \
  --output artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra7.json
sage -python elkies-k3/scripts/probe_h92_q8_e3_generic_module_modp.sage \
  --kernel artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra7.json \
  --output artifacts/generated-results/elkies-k3-h92-q8-e3-generic-module-mod-43-extra7.json
```

This excludes only the bounded modular `r=7` smooth survivor space; it does
not give a q8 pencil or an equation-level bisection.
<!-- status-consumer: EC-K3-H3-Q8-EXTRA7-E3-GENERIC-OBSTRUCTION 41e4df0f0b84986e -->

The stronger all-component generic layer independently excludes the same
mod-43 smooth kernel before node calculations.  Its six singleton coordinate
rows have rank two on the two-dimensional smooth kernel, so no direction
survives.  This uses all actual E7 component valuations but still leaves the
nonsingleton-residue and overlap compiler work for larger ambients.

```bash
sage -python elkies-k3/scripts/assemble_h92_q8_endpoint_rr_ambient.sage \
  --extra-h-power 7 \
  --output artifacts/generated-results/elkies-k3-h92-q8-extra7-endpoint-rr-ambient.json
sage -python elkies-k3/scripts/derive_h92_q8_all_component_generic_conditions.sage \
  --ambient artifacts/generated-results/elkies-k3-h92-q8-extra7-endpoint-rr-ambient.json \
  --output artifacts/generated-results/elkies-k3-h92-q8-all-component-generic-conditions-extra7.json
sage -python elkies-k3/scripts/probe_h92_q8_all_component_singleton_modp.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-all-component-singleton-mod-43-extra7.json
```

<!-- status-consumer: EC-K3-H3-Q8-EXTRA7-ALL-COMPONENT-SINGLETON-OBSTRUCTION 2b07df0cd2f4e19f -->

At `r=10`, the current 774-column envelope has a 16-dimensional smooth
kernel modulo each of `43,53,59,89`; the generic unmarked-E7 screen retains
four directions. The actual `E7_4--E7_3` node then has four successive unique
Pareto-leading negative bidegrees and eliminates all four at every listed
prime. Reproduce the mod-43 three-stage bounded screen with:

```bash
sage -python elkies-k3/scripts/probe_h92_q8_smooth_principal_parts_modp.sage \
  --prime 43 --extra-h-power 10 --include-kernel \
  --output artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra10.json
sage -python elkies-k3/scripts/probe_h92_q8_unmarked_e7_generic_module_modp.sage \
  --kernel artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra10.json \
  --output artifacts/generated-results/elkies-k3-h92-q8-unmarked-e7-generic-module-mod-43-extra10.json
sage -python elkies-k3/scripts/probe_h92_q8_e7_4_3_node_module_modp.sage \
  --kernel artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra10.json \
  --generic artifacts/generated-results/elkies-k3-h92-q8-unmarked-e7-generic-module-mod-43-extra10.json \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-node-module-mod-43-extra10.json
```

This only obstructs the stated modular envelope; it neither completes
resolved-chart gluing nor produces a q8 pencil or bisection.
<!-- status-consumer: EC-K3-H3-Q8-EXTRA10-NODE-OBSTRUCTION 3152bae4a804d3d8 -->

The next three declared denominator enlargements are likewise excluded at the
same prime. The `(columns, rank, smooth-kernel, generic survivors)` tuples
for `r=13,16,19` are respectively `(990,952,38,26)`,
`(1206,1156,50,38)`, and `(1422,1372,50,38)`. In every case the actual
`E7_4--E7_3` node supplies exactly as many independent necessary leading rows
as remaining directions, so no mod-43 survivor remains in these envelopes.

```bash
for r in 13 16 19; do
  sage -python elkies-k3/scripts/probe_h92_q8_smooth_principal_parts_modp.sage \
    --prime 43 --extra-h-power "$r" --include-kernel \
    --output "artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra$r.json"
  sage -python elkies-k3/scripts/probe_h92_q8_unmarked_e7_generic_module_modp.sage \
    --kernel "artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra$r.json" \
    --output "artifacts/generated-results/elkies-k3-h92-q8-unmarked-e7-generic-module-mod-43-extra$r.json"
  sage -python elkies-k3/scripts/probe_h92_q8_e7_4_3_node_module_modp.sage \
    --kernel "artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra$r.json" \
    --generic "artifacts/generated-results/elkies-k3-h92-q8-unmarked-e7-generic-module-mod-43-extra$r.json" \
    --output "artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-node-module-mod-43-extra$r.json"
done
```

This is a one-prime bounded leading-term obstruction, not a complete q8
module or a rootless-equation result.

<!-- status-consumer: EC-K3-H3-Q8-EXTRA13-NODE-OBSTRUCTION 0cd4414b09443b94 -->

Changing the ambient shape by one unit of individual E7-pole slack does not
produce a low-window escape. At `r=4` the 360-column smooth block has kernel
zero. At `r=7`, its 9-dimensional smooth kernel is eliminated by unmarked
generic E7 rows. At `r=10`, the 25-dimensional smooth kernel leaves 13 after
those rows, and the actual `E7_4--E7_3` node eliminates all 13. These are
necessary unmarked/node conditions, so no unproved marked-E7 slack quotient
is being used to reject the windows.

```bash
for r in 4 7 10; do
  sage -python elkies-k3/scripts/probe_h92_q8_smooth_principal_parts_modp.sage \
    --prime 43 --extra-h-power "$r" --extra-e7-pole 1 --include-kernel \
    --output "artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra$r-e7slack1.json"
  sage -python elkies-k3/scripts/probe_h92_q8_unmarked_e7_generic_module_modp.sage \
    --kernel "artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra$r-e7slack1.json" \
    --output "artifacts/generated-results/elkies-k3-h92-q8-unmarked-e7-generic-module-mod-43-extra$r-e7slack1.json"
  sage -python elkies-k3/scripts/probe_h92_q8_e7_4_3_node_module_modp.sage \
    --kernel "artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra$r-e7slack1.json" \
    --generic "artifacts/generated-results/elkies-k3-h92-q8-unmarked-e7-generic-module-mod-43-extra$r-e7slack1.json" \
    --output "artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-node-module-mod-43-extra$r-e7slack1.json"
done
```

This is still a one-prime bounded screen, not a classification of slack
ambients or a q8-pencil construction.

<!-- status-consumer: EC-K3-H3-Q8-E7-SLACK1-OBSTRUCTION 419a01f2b59d42ea -->

The `E7_4--E7_3` node now also has its exact H92 local-module frame.  Its
actual equation is `Y^2-U*H(Z,U)=0` with `H(0)=1`; the P1 corrections are
strictly higher than `x` and `y` in the completed chart.  Hence the q6 module
is `t*R` and the q8 node condition is precisely
`Z^4*Y^6*f/t^9 in R`.  This validates the local frame used by the modular
screen, but does not yet expand its finite two-variable quotient.

```bash
sage -python elkies-k3/scripts/certify_h92_q8_e7_4_3_node_frame.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-node-frame.json
```

<!-- status-consumer: EC-K3-H3-Q8-E7-4-3-NODE-FRAME c3d6b194bb72dd59 -->

That principal-frame derivation now covers the other four unmarked edge
nodes as well. In edge order, their actual q8 Cartier factors are
`U^4Y^2`, `Z^4Y^6`, `U^5Y^6`, `Z^5Y^5`, and `Z^3Y^6`; every corresponding q6
module is exactly `t*R`. The `E7_2--E7_5` node remains outside this uniform
strict-order proof only because its genuine P1 leading cancellation requires
a separate exact calculation.

```bash
sage -python elkies-k3/scripts/certify_h92_q8_unmarked_e7_node_frames.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-unmarked-e7-node-frames.json
```

<!-- status-consumer: EC-K3-H3-Q8-UNMARKED-E7-NODE-FRAMES 59977f99242c92c6 -->

The sixth `E7_2--E7_5` node passes the corresponding cancellation-sensitive
audit: `x-x(P1)=Z^3U^2*unit` and `y-y(P1)=Z^3U^2Y*unit`, so its node module
is likewise `t*R` and its q8 factor is `Z^6Y^5`. This does not conflate the
node with the separate marked smooth point `-P1` on E7₅.

```bash
sage -python elkies-k3/scripts/certify_h92_q8_e7_2_5_node_frame.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-2-5-node-frame.json
```

<!-- status-consumer: EC-K3-H3-Q8-E7-2-5-NODE-FRAME eaff4f299d528774 -->

The six actual frames are now compiled into one two-parameter leading
principal-part template. It has 260 negative bidegree groups on the 54-column
seed, but only one independent exact initial coordinate row after Pareto
deduplication; this remains necessary node data only, since the other groups
need the finite local quotient before becoming matrix rows.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_all_e7_node_principal_bidegrees.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-all-e7-node-principal-bidegrees.json
```

<!-- status-consumer: EC-K3-H3-Q8-ALL-E7-NODE-PRINCIPAL-BIDEGREE-TEMPLATE f9faa97d398f8f02 -->

The one safe Pareto-leading node row is now a composable exact compiler block:
it has rank one on the 54-column seed. The residual node groups still require
their finite local quotient rather than a leading-order shortcut.

```bash
sage -python elkies-k3/scripts/compile_h92_q8_initial_node_conditions.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-initial-node-conditions.json
```

<!-- status-consumer: EC-K3-H3-Q8-INITIAL-NODE-CONDITION-BLOCK a2294ddca25d7344 -->

The E7₄–E7₃ quotient is now reduced to its actual principal frame: after
clearing only explicit node units, the 54-column seed requires common-numerator
membership in `t^17 R`. This is the structured finite-quotient input, rather
than the raw ten-generator singular-model power ideal.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_e7_4_3_principal_node_clearing.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-principal-node-clearing.json
```

<!-- status-consumer: EC-K3-H3-Q8-E7-4-3-PRINCIPAL-NODE-CLEARING dc4e3312a8eb2cef -->

The resulting `t^17` quotient is one-dimensional, not a finite node block:
in actual parameters it is `Z^51Y^68` times a unit.  The finite 3468-dimensional
rectangle `(Z^51,Y^68)` is a different corner jet and is explicitly rejected
as a substitute for product divisibility.

```bash
sage -python elkies-k3/scripts/certify_h92_q8_e7_4_3_node_divisibility_geometry.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-4-3-node-divisibility-geometry.json
```

<!-- status-consumer: EC-K3-H3-Q8-E7-4-3-NODE-DIVISIBILITY-GEOMETRY e287ab81214ee330 -->

For the recorded prime `43`, the stronger all-component generic screen uses
the audited `E7_5` order `ord(m)=0`, rather than inferring it from the
singular equation.  It successively removes every unique live negative-order
coefficient on all seven actual E7 components, giving 16 independent rows on
the 16-dimensional smooth kernel and no survivor.  Thus this one-prime
bounded envelope needs no node screen to be obstructed; groups with multiple
live residues still require the chart-function-field compiler for any larger
ambient.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_all_component_generic_conditions.sage \
  --ambient artifacts/generated-results/elkies-k3-h92-q8-extra10-endpoint-rr-ambient.json \
  --output artifacts/generated-results/elkies-k3-h92-q8-all-component-generic-conditions-extra10.json
sage -python elkies-k3/scripts/probe_h92_q8_all_component_generic_module_modp.sage \
  --kernel artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-mod-43-extra10.json \
  --conditions artifacts/generated-results/elkies-k3-h92-q8-all-component-generic-conditions-extra10.json \
  --output artifacts/generated-results/elkies-k3-h92-q8-all-component-generic-module-mod-43-extra10.json
```

<!-- status-consumer: EC-K3-H3-Q8-EXTRA10-ALL-COMPONENT-GENERIC-OBSTRUCTION 8fa617eca25dc71e -->

The first non-singleton residue stage now works in actual component function
fields, rather than merely comparing valuations. On `E7_4--E7_3` and
`E7_3--E7_7`, the resolved surface equations normalize `E7_4` and `E7_7` and
the strict marked-chord inequalities give the exact leading chord `m=y/x`.
For `r=10` this produces 124 and 131 residue rows. Modulo `43`, their
restriction alone has rank 16 on the 16-dimensional smooth kernel, leaving
no survivor before a node calculation. This is still only a two-component,
one-prime screen.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_e7_4_7_generic_residue_rows.sage \
  --conditions artifacts/generated-results/elkies-k3-h92-q8-all-component-generic-conditions-extra10.json \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-4-7-generic-residue-rows-extra10.json
sage -python elkies-k3/scripts/probe_h92_q8_e7_4_7_generic_residues_modp.sage \
  --residues artifacts/generated-results/elkies-k3-h92-q8-e7-4-7-generic-residue-rows-extra10.json \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-4-7-generic-residues-mod-43-extra10.json
```

<!-- status-consumer: EC-K3-H3-Q8-E7-4-7-GENERIC-RESIDUES f4ea27a1597d3df7 -->
<!-- status-consumer: EC-K3-H3-Q8-EXTRA10-E7-4-7-GENERIC-RESIDUE-OBSTRUCTION 146750c0a39e15f4 -->

The actual conic components now have the same treatment. On `E7_5` and
`E7_6`, set `Z=0`, clear the explicit t-and-chord denominator, and reduce in
the true coordinate ring `QQ(U,Y)/(F(0,U,Y))`. The seed calculation yields
228 E7₅ rows and 57 E7₆ rows; the E7₅ chord keeps the audited P1 cancellation
instead of substituting a guessed normal form.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_e7_5_6_generic_residue_rows.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e7-5-6-generic-residue-rows.json
```

<!-- status-consumer: EC-K3-H3-Q8-E7-5-6-GENERIC-RESIDUES de263de56ed5c2dc -->

The base 54-column seed now has an exact generic-residue cover of all seven
E7 components.  Actual Y-branch normalizations supply E7₁--E7₃, actual conic
coordinate rings supply E7₅--E7₆, and the normalized edge charts supply E7₄
and E7₇.  The row counts are `(42,189,391,42,228,57,34)`, hence 983 in total.
This is not yet a global q8 pencil: node, marked-branch, overlap, and gluing
conditions have not been folded into this certificate.

```bash
sage -python elkies-k3/scripts/derive_h92_q8_e7_1_3_generic_residue_rows.sage
sage -python elkies-k3/scripts/derive_h92_q8_e7_5_6_generic_residue_rows.sage
sage -python elkies-k3/scripts/derive_h92_q8_e7_4_7_generic_residue_rows.sage
sage -python elkies-k3/scripts/certify_h92_q8_all_generic_e7_residue_cover.sage
```

<!-- status-consumer: EC-K3-H3-Q8-ALL-GENERIC-E7-RESIDUE-COVER e3ea804e3f4dd675 -->

Stacking the 983 exact generic-residue rows with the 22 singleton
generic-component rows gives an actual 1005-by-54 q8 compiler block:

```bash
sage -python elkies-k3/scripts/compile_h92_q8_generic_component_conditions.sage
```

Its exact rational rank is 54, so the least endpoint ambient has zero kernel
already on the generic E7 component cover. This rejects only that bounded
ambient; it neither replaces node/overlap conditions nor rules out the
enlarged Riemann--Roch space needed for a q8 pencil.

<!-- status-consumer: EC-K3-H3-Q8-GENERIC-COMPONENT-CONDITION-BLOCK 33693f196eb13091 -->

The `r=7` endpoint enlargement is also excluded before node work. Its
actual smooth block has a two-dimensional kernel mod 43, but the full set of
2,487 actual generic E7 rows has restriction rank two on it. Thus the stacked
smooth-plus-generic matrix has full rank 558 mod 43, which proves the same
characteristic-zero rank. This rejects this one enlargement, not larger
Riemann--Roch ambients.

<!-- status-consumer: EC-K3-H3-Q8-EXTRA7-SMOOTH-GENERIC-REJECTION e0e67861f23d4f24 -->

The next `r=8` enlargement is excluded by the same transported-chart block:
its 630-column smooth matrix has mod-43 rank 624, and all actual generic E7
residue rows have rank six on the resulting six-dimensional smooth kernel.
The stacked matrix is therefore full rank modulo 43 and hence over QQ. This
remains an ambient rejection before nodes, marked branches, overlaps, or E8;
it does not construct a q8 pencil or rule out larger Riemann--Roch spaces.

<!-- status-consumer: EC-K3-H3-Q8-EXTRA8-SMOOTH-GENERIC-REJECTION 0d66d55dd153a589 -->

At `r=9` the same check gives smooth rank 690 and generic restriction rank 12
on the twelve-dimensional kernel, so the 702-column stacked block is full rank
modulo 43 and over QQ. This rejects just that transported-chart ambient before
the node, marked-branch, overlap, and E8 layers.

<!-- status-consumer: EC-K3-H3-Q8-EXTRA9-SMOOTH-GENERIC-REJECTION 8035f0b465038995 -->

The next nested enlargement has the same outcome.  At `r=11`, the 846-column
smooth block has rank 822 modulo `43`; the complete resolved E7 residue cover
has rank 24 on its 24-dimensional kernel, so the stacked block is full rank
over both `F_43` and `QQ`.  This is another endpoint-envelope rejection only,
before node, marked-branch, overlap, and E8 conditions.

<!-- status-consumer: EC-K3-H3-Q8-EXTRA11-SMOOTH-GENERIC-REJECTION c755bb2f34d68f41 -->

The four smooth P1.O collision fibres now contribute an exact q8 condition
block for the least endpoint seed.  In the transported frame
`q=(m-y(P1)/x(P1))/h`, `X=h^2*x`, all negative `h`-principal parts give a
1080-row condition template on 54 coefficients. Its finite ambient image is
compiled through the shared condition interface, and its mod-43 reduction has
full column rank, which proves characteristic-zero rank 54; hence this seed
has no smooth-compatible direction.  This is only a local obstruction for
the least ambient, not a completed q8 cover or a nonexistence theorem for
enlargements.

```bash
sage -python elkies-k3/scripts/compile_h92_q8_smooth_principal_parts_exact.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-smooth-principal-parts-exact.json
```

<!-- status-consumer: EC-K3-H3-Q8-SMOOTH-PRINCIPAL-PART-CONDITION-BLOCK 6888f78205d642be -->

### Source q=8 E8 target

The q=8 E8 component degrees in the pinned source order are
`(1,0,0,0,0,0,0,0)`.  Because the E8 Cartan matrix is unimodular, they give a
unique integral exceptional cycle:

```text
(-4,-5,-7,-10,-8,-6,-4,-2).
```

Equivalently, the q=8 E8 module is the ninth tensor power of the certified
q=6 module `u*<1,Q>` followed by that exceptional twist.  This target is
replayed by:

```bash
sage -python elkies-k3/scripts/derive_h92_q8_e8_local_target.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e8-local-target.json
```

The source cycle is now assigned to the eight ordinary E8 blow-up components:
in chart order `(B1,B2,B3,B4,N3,N40,N4B,N4inf)` it is
`(-2,-4,-6,-10,-4,-7,-5,-8)`, with degree vector
`(0,1,0,0,0,0,0,0)`.  It still does not supply the finite local quotient
condition.
<!-- status-consumer: EC-K3-H3-Q8-E8-TARGET 7f4bac8ed72db930 -->

### Complete q=8 E8 module

The chart-cycle has an unexpectedly small complete ideal.  The exact
exceptional valuations of `u`, `X`, and `Y` prove that its integral twist is

```text
(u^2, X, Y),
```

with quotient basis `1,u` and colength two.  As the q=6 unit chord `Q` is a
unit at the E8 singularity, the full q=8 E8 module is
`u^9*(u^2,X,Y)`.  Reproduce it with:

```bash
sage -python elkies-k3/scripts/derive_h92_q8_e8_complete_module.sage \
  --output artifacts/generated-results/elkies-k3-h92-q8-e8-complete-module.json
```

The bounded global ambient and the E7 finite quotient remain to be assembled.
<!-- status-consumer: EC-K3-H3-Q8-E8-MODULE 74327bc7489c8ca6 -->

### Bounded q=6-child polynomial reconnaissance

The selected `q=8` class has old-fibre degree two on the exact `E8+E6/MW3`
q=6 child, even though its source-side ambient has degree 18.  As a narrow
reconnaissance test, the script below first rejects primes that do not retain
the child's `II*` and `IV*` valuations.  At the good prime `43` it then
exhausts the 79,507 degree-`<=4` polynomial `x` coordinates on the explicitly
declared singular `IV*` branch `f_IV^2 | x`, requiring `f_IV^2 | y` and
`deg(y)<=6`:

```bash
sage -python elkies-k3/scripts/search_h92_q6_child_polynomial_sections_modp.sage \
  --prime 43 --max-x-degree 4 --require-iv-star-singular \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-polynomial-sections-mod-43-iv-singular.json
```

The mod-43 output has three `x` coordinates and their six `y` sign choices,
all with recorded IV* orders `(2,2)` and smooth II* specialization.  This is
an exhaustive finite-field sample only for that stated ansatz and local
branch.  It neither establishes a characteristic-zero lift nor identifies a
q=8 pencil; in particular it must not be read as a no-collision or no-rank-19
result.

The accompanying Hensel diagnostic currently finds coefficient-Jacobian rank
11 rather than 12 for each of the six mod-43 residues (and likewise for the
six corresponding mod-53 residues), so they are not unique p-adic lifts:

```bash
sage -python elkies-k3/scripts/lift_h92_q6_child_polynomial_sections.sage \
  --input artifacts/generated-results/elkies-k3-h92-q6-child-polynomial-sections-mod-43-iv-singular.json \
  --precision 64 \
  --output artifacts/generated-results/elkies-k3-h92-q6-child-polynomial-sections-hensel-43.json
```

This is a guard against incorrectly promoting the modular residues to the
three characteristic-zero MW generators.  A deformation-aware local analysis
is required before using them in a q8 construction.

The older exploratory
[`scripts/build_h92_q6_chord.sage`](scripts/build_h92_q6_chord.sage) must not
be used as an alternative certificate. On the current pinned inputs its final
squarefree-cover assertion fails: the raw odd factor has degree 21 rather than
the asserted degree 3 or 4. Thus it produces neither a valid first neighbor
equation nor a branch divisor for this collision search.

## Direct rootless bisection compiler and first equation batch

The compact published rootless equation and its seventeen exact section/chord
records remove the former equation-level entrance obstruction.  If a
height-ten trace section is written

```text
tau=(Nx/h^2,Ny/h^3),       deg(h)=3,
```

then the residual line through `-tau` is recovered without a nonlinear solve.
There is a unique polynomial `M` of degree below six such that

```text
M*Nx+Ny == 0 mod h^2.
```

For the slope `m=M/h`, exact substitution gives the residual quadratic

```text
x^2-s(t)*x+p(t)=0
```

with

```text
s=(M^2-Nx)/h^2,
p=((M*Nx+Ny)^2-B*h^6)/(h^4*Nx).
```

The elliptic equation forces

```text
s^2-4p = h^2*q(t),
q = (M^4-6*M^2*Nx-8*M*Ny-3*Nx^2-4*A*h^4)/h^6,
deg(q)=2.
```

The proof, including the `h^6` divisibility, is Proposition F1 of
[`RANK_MUTATION_AND_LIFT_THEOREMS.md`](RANK_MUTATION_AND_LIFT_THEOREMS.md).
Thus the squareclass is read directly from `q`, while the two points over
`u^2=q(t)` have exact coordinates

```text
x=s/2 + (h/2)u,
y=y0(t) + (M/2)u.
```

The exact regressions in
[`scripts/verify_elkies_2026_rank18_first_cover.sage`](scripts/verify_elkies_2026_rank18_first_cover.sage)
and
[`scripts/verify_elkies_2026_rank19_paired_cover.sage`](scripts/verify_elkies_2026_rank19_paired_cover.sage)
recover the two published chord numerators from this single congruence before
checking their displayed cover sections.

### Exact equation-complexity priority

The complete survivor set is ranked in the published Mordell--Weil basis by a
reproducible tuple rather than a floating heuristic:

```text
(group-addition upper bound,
 support size,
 recursive section/chord dependency count,
 exact serialized coordinate-bit cost,
 maximum coefficient,
 coefficient L1 norm,
 published-basis vector).
```

The scalar-multiplication component uses the standard binary double-and-add
upper bound.  Both orientations and every norm-ten representative of an orbit
are tested, so the retained representative is the cheapest for this declared
score, not merely the lexicographic short-basis representative.  This is an
exact arithmetic-cost proxy; it is not a theorem that coefficient height or
wall-clock time is monotone in the ordering.

```bash
sage -python elkies-k3/scripts/rank_elkies_2026_bisection_orbits.sage \
  --pool-size 39120 \
  --output artifacts/generated-results/elkies-2026-bisection-equation-priority-full.json \
  --table-output artifacts/generated-results/elkies-2026-bisection-equation-priority-full.tsv \
  --pairs-output artifacts/generated-results/elkies-2026-bisection-equation-priority-disjoint-pairs-full.tsv
```

The complete ranking has 77 orbits of group-addition cost two, followed by
547 of cost three.  The first 1,000 equation-cheapest orbits contain 11,823
exact disjoint-priority pairs.  The two previously published trace vectors
occur much later in this score, so they are validation anchors rather than
appropriate starting points for the batch search.  The full priority table
contains all 39,120 classes and its full graph table contains all 8,895,801
disjoint pairs.

### Complete 39,120-orbit equation batch and injectivity certificate

The batch compiler forms each trace by exact group law, solves the linear
double-pole congruence, verifies the quadratic relation and both Weierstrass
coefficient identities over `u^2=q(t)`, and emits the pinned orbit mask and
vector together with the exact relation accepted by the squareclass checker.

```bash
sage -python elkies-k3/scripts/construct_elkies_2026_bisections.sage \
  --priority-table artifacts/generated-results/elkies-2026-bisection-equation-priority-full.tsv \
  --limit 39120 \
  --output artifacts/generated-results/elkies-2026-equation-bisections-full.json \
  --orbits-output artifacts/generated-results/elkies-2026-equation-bisections-orbits-full.tsv

.venv/bin/python elkies-k3/scripts/hash_bisection_extensions.py \
  --compact \
  --input artifacts/generated-results/elkies-2026-equation-bisections-full.json \
  --output artifacts/generated-results/elkies-2026-equation-bisection-collisions-full-compact.json
```

All 39,120 records pass the exact quadratic relation, lifted-section,
two-branch, and orbit-attachment gates.  Exactly one trace, orbit `0x0c54f`,
has a finite denominator of degree below three because a pole lies at infinity;
the compiler performs the same construction after `t -> 1/t` and transports
the certified relation back.  The other 39,119 records use the published
affine chart.  Every branch quadratic is coprime to the degree-24 surface
discriminant, so every cover branches at two smooth fibres.

Exact normalization gives 39,120 distinct elements of
`QQ(t)^*/QQ(t)^{*2}` and no collision.  Since the input coverage gate proves
that there is exactly one equation record for every surviving `R17/2R17`
translation orbit, this is a complete finite injectivity certificate for the
rootless bisection-to-quadratic-extension map on this surface.  In particular,
none of the 8,895,801 disjoint pairs has a common quadratic cover, and this
entire norm-ten rootless-bisection mechanism cannot produce the desired
rank-two anti-invariant collision or a new generic rank-19 family.  It does
produce 39,120 explicit smooth quadratic covers, each with its exact split
bisection section, for generic-rank-at-least-18 base changes.

The pinned replay hashes are:

```text
full priority summary       c1ce639e21f773148b800c3b905cb87d118d31881c159380bdb5c60f4e58d480
full priority table         1296c8a81c8df49757a4308f7abb087035507e16c5f41cbd2a257343fa3eb166
full disjoint-pair table    15c00a0374c683def5c88a77f130ad651db0f48467535c08cfb54bb5ecc5a3e2
full equation batch         78e037dc4170955b8f79ddce4d1d3e0c0d3e9bb8f9614644c59ccc7d605226c4
full orbit attachment table 6edfb82b8dad06b5e3a0c26a0045b999efaf4d34698387224d1d44800c18a85b
compact collision report    8e81214901bb75760f662fc87d6d30c1d1941c41ff1de8954ad1113449ee3d19
extension manifest          e2419266fc527090a23a6ab4d7bee8f3ca37a0f6364cb482b3192f883d88ca73
```

<!-- status-consumer: EC-K3-BISECT-EQUATION-BATCH a0570a5a4ea8e02b -->

## Collision and height gate

The quotient enumeration supplies a finite input set, not a squareclass hash.
Mapping a lattice class to a quadratic extension requires the explicit
rootless characteristic-zero fibration and an equation for the corresponding
bisection. Those inputs now exist for every one of the 39,120 priority classes.
The checker:

1. derive and squarefree-normalize the branch divisor for every realized
   orbit representative;
2. hash its quadratic extension by squareclass (or unordered branch pair);
3. for each hash collision, construct the anti-invariant sections on the
   common quadratic twist; and
4. compute their Shioda height matrix and test the rank increment.

The deterministic post-processing component is already available as
[`scripts/hash_bisection_extensions.py`](scripts/hash_bisection_extensions.py).
It accepts exact numerator and denominator coefficient lists for a cover
`z^2=f(t)`, first cancels common factors, retains the rational constant
squareclass, factors the odd monic polynomial support, and groups exactly
equal classes in `QQ(t)^*/QQ(t)^{*2}`. It verifies that the resulting cover
has exactly two geometric branch points (including infinity when appropriate),
as required for a connected smooth rational bisection; split and higher-genus
quadratic covers are rejected. It also accepts a bisection already eliminated
to an exact quadratic relation `a(t)z^2+b(t)z+c(t)=0` and derives its branch
squareclass from `b^2-4ac`.  When the resolved record also supplies `h` and
`q`, it verifies `b^2-4ac=h^2*q` exactly before normalizing the smaller
quadratic.  Compact mode retains the complete label/orbit/extension-digest
manifest and all collision provenance, pins the input SHA-256, and omits only
bulky singleton provenance. When a collision is supplied with coordinates in
a declared twist-height lattice, it also returns the exact height matrix and
its rational rank. It does not manufacture those coordinates from a collision.
Whenever lattice masks are supplied, each mask must be globally unique across
the full input: section translates induce the same quadratic cover and the
same anti-invariant direction even if inconsistent branch data would hash
them differently.

For a complete rootless run, `required_lattice_orbits.table` is set to the TSV
above (with its SHA-256 pinned) and every record has its `lattice_orbit_mask`.
Also give each record the 17 integral `pinned_rank17_w` coordinates of the
realized bisection class.  The checker verifies its norm is at least 10 and
is `2 mod 4`, then verifies that it is congruent modulo `2M` to the
representative of the claimed orbit.  It rejects omitted, duplicate,
extraneous, or mis-attached orbit records. This binds the complete
equation-level table to all 39,120 translation orbits while still allowing
any section-translate as the equation's representative.

For the alternate q80 rootless frame, declare the alternate TSV together
with its frame artifact and use `alternate_rank17_w` as the record key.  The
same exact coverage gate then validates the vector norm and its class modulo
twice that alternate lattice; the synthetic regression exercises both frame
modes.

Its independent synthetic regression is:

```bash
.venv/bin/python elkies-k3/scripts/hash_bisection_extensions.py --self-test \
  --require-collision-heights --require-rank-at-least 19
```

The gate validates positive definiteness of every declared rational twist
height lattice; a merely symmetric or indefinite matrix cannot be used to
manufacture a rank contribution.

There is also a geometric route once a common cover has been certified to
branch only at smooth fibres and its rootless pullback is known.  For lifted
sections `P_i` and the deck involution `tau`, the checker verifies
`P_i^2=-4` and `P_i.tau(P_i)=2`, then computes

```text
<P_i-tau(P_i), P_j-tau(P_j)> = 2*(P_i.tau(P_j)-P_i.P_j).
```

Thus the equation-level intersection calculation can supply the exact height
matrix directly; if it is supplied alongside twist coordinates, both matrices
must agree.

<!-- status-consumer: EC-K3-BISECT-EXTENSION-PROTOCOL 90dc72ea57ae22dc -->

For a production collision run, add `--require-collision-heights`.  It rejects
any collision bucket lacking either declared twist-height coordinates or the
certified geometric lift-intersection data, rather than silently returning a
collision with missing height data.  A computed bucket records its exact height-matrix rank;
that rank is the certified anti-invariant rank contribution of the supplied
sections.  Add `invariant_mw_rank: 17` to a rootless production input to have
the tool report the base-changed rank lower bound as `17` plus that rank. The
deck involution splits the rational Mordell--Weil space into invariant and
anti-invariant eigenspaces, so a collision height matrix of rank two then
certifies the lower bound `19`. The optional
`--require-rank-at-least 19` rejects a run unless at least one collision meets
that declared lower bound.

Whenever lattice masks are supplied, each `lattice_orbit_mask` must be unique
across the entire run. Mixing masked and unmasked records inside one
squareclass bucket, or listing two section translates from the same `M/2M`
orbit even with inconsistent branch data, is rejected before hash grouping.
This rules out a duplicate geometric bisection masquerading as a rank-two
collision.

All 39,120 pinned classes now have exact equation-level records and pairwise
distinct squareclasses.  Thus the pinned map is injective on its complete
survivor set and supplies no common quadratic cover with a rank-two
anti-invariant height matrix.  The 39,147 alternate classes are still only
exact lattice candidates because that alternate finite-field endpoint has not
been lifted to a characteristic-zero rootless surface.

Injectivity does not end the paired-cover programme.  Taking two distinct
quadratic extensions gives a biquadratic base whose two new sections have
different Galois characters.  The complete follow-up classification in
[`BISECTION_PAIR_COVER_GEOMETRY_2026-08-31.md`](BISECTION_PAIR_COVER_GEOMETRY_2026-08-31.md)
proves that all 39,120 conics are rational, every distinct pair has a genus-one
base and exact height matrix `diag(24,24)`, and the first new low-complexity
pair has base Jacobian rank at least 3.  Neither that pair nor the published
rank-19 pair lies in the norm-four disjointness graph, so that graph is a
priority heuristic rather than a valid hard filter for character independence.
<!-- status-consumer: EC-K3-BISECT-BIQUADRATIC-R19 707bffd8b85f8f3e -->
<!-- status-consumer: EC-K3-BISECT-ORBIT 81da2fd80c3623b6 -->

The alternate q80 q6 endpoint supplies a second, nonisometric rootless
rank-17 lattice that is closer to the available finite-field equation route.
The streamed exact quotient has 39,147 section-nonnegative bisection orbits,
with 805,466 unoriented norm-ten representatives.  It is independently
cross-checked against PARI's signed short-vector count through norm ten;
the finite-field endpoint has not been lifted to a characteristic-zero
rootless surface, so these are not yet branch-cover records.

```bash
sage -python elkies-k3/scripts/enumerate_rootless_bisection_orbits.sage \
  --frame-artifact artifacts/generated-results/q80-alternate-fifth-q6-rootless-transport.json
```

<!-- status-consumer: EC-K3-ALT-BISECT-ORBIT eca5fc0bfee5038d -->
