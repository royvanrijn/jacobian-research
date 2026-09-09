# Generic rational curves: six constructions and an infinite-orbit obstruction

The subsequent [intrinsic source-degree theorem](DET1092_INTRINSIC_LOW_DEGREE_SOURCE_OBSTRUCTION_2026-09-09.md)
excludes all original section sources of alternate degree at most2 and
their entire indicated orbits. The present higher-degree displayed-source
exclusions and six construction maps remain separate valid evidence.

## Main result and scope

**New verified obstruction.** Fix the smooth norm-eight carrier `C_(z*)`
through the historical first302 seed. Start with `O` and both signs of
the17 displayed original generic section curves. No composition of
translations by arbitrary generic sections of the **alternate fibration**,
and its elliptic inversion, can move a rational point from this source set
to the first seed or its carrier conjugate.

This excludes an infinite family of rational curves, not just finitely
many small multiples. It does **not** exclude all original MW section
curves, automorphisms alternating between fibrations, or other points
representing the same original exceptional quotient direction.

**New construction over Q.** A separately frozen, entirely generic rule
constructs six rational covers of degrees13,14,15,25,26,27, each with a
certified independent extra section after pullback, hence generic rank
at least18. All six have no rational preimage of302. They remain useful
rank18 families; they do not construct the requested302 seed.

**Verified closure of an earlier diagnostic.** Both previously constructed
degree20 and58 covers miss all nine fixed302/null-control parameters.
The five earlier UNKNOWN entries are now excluded by exact rational-root
certificates. Neither the covers, the panel, nor the old prime pool was
changed. The initial bounded records remain intact.

## Why the signed-source obstruction holds

**Established setup, reused.** The norm-eight elliptic pencil has parameter
`z` and origin `S14`; write `J/Q(z)` for its Jacobian with this origin.
On the first-seed member, the full Picard image

\[
 H_{z_*}=\rho_{z_*}(\operatorname{NS}(X))
 \subset J_{z_*}(\mathbf Q)
\]

has rank12, and the first marked class is outside its rational span.
Both statements, including the conjugate's exclusion, were already
proved in the [full-Picard comparison](DET1092_GENUS_ONE_PICARD_SPECIFICITY_2026-09-08.md).
Every generic section of `J` specializes inside this Picard image.

**Verified source incidence.** For an original section `S_i`, use its short
coordinates `X_i,Y_i` and the generic pointed-quartic data `cx,cy,h,shift`:

\[
 m_i^\pm=\frac{\pm Y_i+cy}{X_i-cx},\qquad
 z_i^\pm(t)=\frac{h m_i^\pm+\mathrm{shift}}{h^2}.
\]

The exact degree equals the lattice intersection
`D.(+/-S_i)=G_ii -/+ (Gw)_i`, with `w=S14-S15`. We check this on all34
curves, then solve `z_i^pm(t)=z*` over the rational projective line.

| Source curves | Exact outcome on `C_(z*)` |
|---|---|
| `+S14`, `-S15` | Degree-one alternate sections, hence inherited classes |
| Other32 signed displayed sections | No rational intersection |
| `O` | Disjoint from the smooth member because `D.O=0` |

All32 absences include the point at infinity. The independent checker
verifies60 Hensel/Gauss certificates; it does not invoke factorization,
root finding, or lattice reduction. Exact functions and certificates are
in the [signed-source packet](../../artifacts/generated-results/elliptic-curves/det1092_signed_source_orbits_v1/).
The only seed-derived input is the already fixed carrier label `z*`;
later points and search outputs are not read.

**New deduction; proof.** A fibrewise translation is a rational
automorphism preserving `z`. Its inverse takes a rational point of a
translated source curve on `C_(z*)` to a rational point of the source
curve on the same carrier. For32 sources no such point exists. For the
remaining two, its class belongs to `H_(z*)`; arbitrary generic
translations and inversion keep it in that rational span. The first
seed and its conjugate are outside the span. This proves the exclusion
for the entire indicated group, without enumerating any translation orbit.

**Control interpretation, not a predictor.** In the nine old generic-point
controls, the carrier label was defined as `z_0(tau)`. Those marked
points lie on source curve `S0` by construction. This explains how an
originally inherited point can nevertheless give the observed alternate
carrier rank increase12 to13. On the first-seed carrier, even that signed
source set is inaccessible. The two label selections differ, so this is
not a candidate-free discriminator on the null fibres.

## Six generic rational-curve constructions

**Frozen rule.** Among the displayed original section curves, take all
with smallest alternate-fibration degree greater than one. The degrees
are computed from generic lattice data, before any specialization:

```text
4,2,4,4,3,3,3,2,3,2,3,3,3,8,1,7,3.
```

The selected sections are `S1,S7,S9`, all degree two. Translate each
exactly once by both signs of the already inherited alternate section
`B`, the conjugate of the origin. No orbit or sign is chosen from outcomes.
Unlike the earlier `2B,3B` route, these sources are multisections of the
alternate fibration, not its sections.

**Explicit construction.** Let `u` parametrize an original selected section.
Its known generic functions give `z=z_i(u)` and `W=W_i(u)`. Substitute this
`z` into the pointed alternate quartic and its inherited origin:

\[
 F_z(t_0+v)=s^2+q_1v+q_2v^2+q_3v^3+q_4v^4.
\]

The monic Jacobian is

\[
 y^2=x^3+q_2x^2+(q_1q_3-4s^2q_4)x
       +s^2q_3^2+q_1^2q_4-4s^2q_2q_4.
\]

Set `v=u-t0` and

\[
 x=\frac{2s(W_i+s)+q_1v}{v^2},\qquad
 y=\frac{(x^2-4s^2q_4)v-q_1x-2s^2q_3}{2s}.
\]

The inherited point has
`bx=q1^2/(4s^2)-q2`, `by=-(q1*bx+2s^2*q3)/(2s)`.
For `epsilon=+/-1`, the ordinary chord formula gives

\[
 \ell=\frac{\epsilon by-y}{bx-x},\quad
 x'=\ell^2-q_2-x-bx,\quad y'=-y+\ell(x-x').
\]

Finally put

\[
 v'=\frac{2sy'+q_1x'+2s^2q_3}{x'^2-4s^2q_4},\quad
 T(u)=t_0+v',\quad
 W'(u)=\frac{x'v'^2}{2s}-s-\frac{q_1v'}{2s}.
\]

Then `W'(u)^2=F_(z_i(u))(T(u))`. Map back to the original elliptic family
by the already verified inverse chart

\[
 m=h(T)z_i(u)-\mathrm{shift}(T)/h(T),\quad
 X=(h(T)W'-cx(T)+m^2)/2,\quad
 Y=m(X-cx(T))-cy(T),
\]

and undo the original short-model coordinate change. These are explicit
rational maps over Q; all coefficients of `T,W',z_i` are saved in
[map-00 through map-05](../../artifacts/generated-results/elliptic-curves/det1092_translated_sections_v1/).

**Verified independent maps and rank.**

| Source | Sign | Degree of `T` | Certifying prime | Smooth ramification support lower bound |
|---|---:|---:|---:|---:|
| S1 | - |27|127|52|
| S1 | + |13|127|23|
| S7 | - |26|131|50|
| S7 | + |14|127|26|
| S9 | - |25|149|48|
| S9 | + |15|127|28|

The independent checker rebuilds the source from the original generic
coordinates, verifies the manual chord and inverse chord, and checks the
quartic and parameter maps exactly. Translation extends to an
automorphism of the minimal elliptic K3; its image of the original section
is therefore a smooth rational curve, and the inverse chord verifies
birationality of the displayed parametrization.

**Established literature, verified application.** Ramification above at
least one smooth original fibre proves the new section independent of
the pulled-back MW17; see
[Garbagnati--Salgado, Lemma2.9](https://www.cambridge.org/core/journals/forum-of-mathematics-sigma/article/rank-jumps-and-multisections-of-elliptic-fibrations-on-k3-surfaces/602807FB00055D106E3CEA5418DE08F7).
For `T=N/D`, the certificates preserve the degrees of `N,D,N'D-ND'`
and the original discriminant. Modulo the indicated prime, the critical
radical has a factor outside `D*D^24*Delta(N/D)`. If every characteristic-zero
critical point were bad, the critical polynomial would divide a power
of this bad support, contradicting that reduction. This proves the
required smooth ramification.

No no-bad-ramification height formula is assumed for these six covers;
the lower bounds above suffice for independence. They give six separate
rank-at-least18 families over `Q(u)`, not a combined rank claim.

**Rational source and exclusions.** The base is `P1_Q`, so rational `u`
are available directly. Exclude zeros of the displayed chart denominators,
`h(T(u))`, and the original discriminant, when exporting affine packets.
These are explicit finite algebraic sets. Specialization independence
also has a finite exceptional set which is not made effective here;
individual exported points would still require rank certification.
No new specialization is exported. All six have no rational projective
preimage of `T=0`; the first eligible old-prime test at53 certifies this.

## Exact closure of the old degree20/58 control cases

**Verified application.** The initial12-prime screen excluded13 of18
cover/control pairs and left five UNKNOWN. Exact Q-factorization returned
no linear factor in each remaining case. A separate elementary proof
now excludes rational projective roots in all five, without factorization.
Thus both old covers miss302 and all eight unchanged null controls.
The [initial evidence](../../artifacts/generated-results/elliptic-curves/det1092_pencil_multiples_v2/controls.json)
and [exact follow-up](../../artifacts/generated-results/elliptic-curves/det1092_pencil_exact_incidence_v1/)
remain separate.

**Elementary certificate, verified application.** For primitive integral
`f` with nonzero constant coefficient, every rational root `a/b` in lowest
terms satisfies `|a|<=|f(0)|`, `|b|<=|lc(f)|`. Put `H` equal to the larger
bound. Choose the first prime from the same pool with nonzero leading
coefficient and only simple roots. Lift every root to a modulus `M=p^k`
with `M>8H^2`. Hensel uniqueness captures every possible rational root.

For each lift `r`, the congruence lattice consists of pairs satisfying
`a-rb=0 mod M`. The certificate supplies an exact Gauss-reduced basis
`v,w`: determinant of absolute value `M`, `||v||<=||w||`, and
`2|v.w|<=||v||^2`. These inequalities prove that `v` is a shortest vector.
If `||v||^2>2H^2`, no rational-root pair is possible. Otherwise every
candidate pair is collinear with `v`, since two such independent vectors
would have determinant at least `M` but absolute determinant at most
`2H^2`. Only one rational ratio remains to test exactly.

For302, the chosen prime is131. Its two roots modulo131 lift to moduli
of17,345 decimal digits; both shortest-vector bounds exclude rational
roots. Several earlier primes in the old screen had a degree drop or
multiple roots. The old small-prime survival therefore was not evidence
of rational incidence. Five positive/negative synthetic regressions also
verify the method's handling of genuine rational roots.

## Reproducibility and what remains open

**Subsequent new deduction.** The first operation involving both fibrations
now has an [exact degree barrier](DET1092_TWO_FIBRATION_DEGREE_GAP_2026-09-09.md):
from the cheapest degree13 curve, a nonzero arbitrary original translation
followed by the fixed alternate move has minimum degree272 or442. This
narrows that switching route; it is not a general automorphism-orbit exclusion.

**Verified checkpoints.** Sage10.9/PARI2.17.3 were used. All jobs completed
within their25-second process caps. No detached job, new prime pool,
parameter scan, class-group/unit-group calculation, point search or
production change occurred. The norm-equation principality fork remains
untouched and unresolved.

- [Infinite-orbit independent replay](../cas/verify_det1092_signed_source_orbits.sage)
- [Six-cover independent replay](../cas/verify_det1092_translated_sections.sage), one `--case 0..5`
- [Old-cover rational-root replay](../cas/verify_det1092_pencil_rational_preimages.sage), one `--case 0..4`

**Open constructive boundary.** The original target is not achieved.
Increasing alternate-fibration multiples of the fixed signed sources
cannot reach the known first seed. A new route must change the source
curve outside that set, change the alternate-fibration parameter through
another operation, or construct the nongeneric arithmetic class directly.
The theorem does not assert that any such change will succeed, and it
does not exclude other seeds on302. Prospective seed construction and a
candidate-free distinction from the null fibres remain open.
