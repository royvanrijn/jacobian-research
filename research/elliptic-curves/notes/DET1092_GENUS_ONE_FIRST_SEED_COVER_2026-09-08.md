# A positive-rank genus-one carrier through the first302 unlock

## Result and calibration boundary

**New constructive deduction, independently verified.** A generic-selected
norm-eight pencil on the determinant1092 parent contains an explicit
degree-two **genus-one** cover through the historical first302 unlock.
Its base curve has an explicit nontorsion rational point and rational maps
to and from a Weierstrass model. The pulled-back original elliptic family
has a certified subgroup of rank at least18.

Thus the first seed has a verified positive-rank genus-one carrier, not
only the previously constructed genus-two carrier. This changes the
available rational-point source: multiples on an explicitly given elliptic
curve supply infinitely many rational base points.

**Selection is still retrospective.** The pencil is chosen using only the
generic lattice and sections. The particular member through the first seed
is then fitted using that permitted first witness. Removing its coordinates
from subsequent execution does not make the member selection prospective.
The generic fixed member `z=0` misses302 and all eight controls; the fitted
member hits302 and misses all eight controls.

**Verified equation-only execution.** A separate worker reads only sanitized
equations, maps, generic sections and the unchanged address roster. It
reconstructs a point over302 and independently proves rank18 using complete
finite group quotients. No saved exceptional coordinate is an execution
oracle. This recovers the known first direction, not a new302 point or an
arithmetic predictor.

## Audit: what existed, and what was missing?

**Previously established applications.** The parent has full geometric
MW17, determinant1092 and24 irreducible nodal fibres. Its complete old
degree-two quotient and the norm-ten orbit8044 rational conic already exist.
The historical first norm-ten RR pencil has no rational lower-genus member
through its first witness. None of those statements excludes a genus-one
member with a different centre.

**Implementation audit.** A committed generic norm-eight pencil constructor
was present, but its result artifact and theorem entry were absent. Its
bounded replay completed the generic algebra, then attempted to construct
a singular finite reduction before testing its discriminant. The original
failure is retained in
[version1](../../artifacts/generated-results/elliptic-curves/det1092_norm8_seed_cover_v1/failure.json).
Version2 moves that good-reduction check before the finite-curve constructor.
The pencil, fixed member, prime pool and limits are unchanged.

**New verified application.** The frozen generic selector minimizes
`(l1, linfinity, orbit)` among the already enumerated norm-eight rows. It
chooses orbit20124 and the trace

\[
Z=P_{14}-P_{15},\qquad \langle Z,Z\rangle=8,
\]

with zero-based indices in the inherited basis. The single parity replay
uses655 nodes within its200000-node limit. The independent proof checks
the stored selector, the exact norm, the equations and smooth members;
it does not repeat a subgroup census.

## Explicit pencil and its elliptic lift

**Verified construction.** In the short coordinates of the original family,
write

\[
X=x+b_2/12,\quad Y=y+(a_1x+a_3)/2,\quad
Y^2=X^3+aX+b,
\]

where `a=-c4/48`, `b=-c6/864`. Write the trace coordinates as `(cx,cy)` and
put

\[
h^2=\operatorname{den}(x_Z),\qquad N_x=h^2c_x,\quad N_y=h^3c_y.
\]

Here `h` has degree two. Let `shift` be the unique degree-below-four
remainder satisfying `shift*Nx=Ny mod h^2`, and put `M_z=h^2*z-shift`.
The exact branch quartic is

\[
\boxed{F_z(t)=
\frac{M_z^4-6N_xM_z^2-8N_yM_z-3N_x^2-4ah^4}{h^6}.}
\]

The numerator is divisible by `h^6`; the quotient is polynomial of degree
four in `t`. The auxiliary curve and degree-two base map are

\[
C_z:W^2=F_z(t),\qquad C_z\longrightarrow\mathbf P^1_t.
\]

With `m=M_z/h`, the extra point on the pulled-back original family is

\[
X=\frac{hW-c_x+m^2}{2},\qquad Y=m(X-c_x)-c_y,
\]
\[
x=X-b_2/12,\qquad y=Y-(a_1x+a_3)/2.
\]

The [generic certificate](../../artifacts/generated-results/elliptic-curves/det1092_norm8_seed_cover_v2/generic.json)
contains all polynomial coefficients. The checker independently verifies
the quartic identity and both coefficients of the elliptic equation modulo
`W^2-F_z(t)`, including the original rational twist factor.

**Inherited rational points.** The old sections `P14` and `-P15` each meet
a pencil member once. Their exact restrictions are `(t0(z),s(z))` and
`(t0(z),-s(z))`, with `t0(z)` fractional linear. They are available from
generic data alone. The old zero section has intersection zero with these
members; this is a different inherited-point source from `D.O=1` in the
genus-two RR net.

## Two fixed members, no member search

**Generic choice.** Member00 is precisely `z=0`, frozen by the old script.

**Retrospective choice.** If `(X*,Y*)` is the historical first point in the
short coordinates at `t=0`, the unique member through it is

\[
z_*=
\frac{h(0)(Y_*+c_y(0))/(X_*-c_x(0))+\operatorname{shift}(0)}{h(0)^2}.
\]

The exact rational value, openly retaining its calibration, is

```text
1018595891632634549154748055919128315329729591633806393135864425102630466026311226984286139170415482796192040753045215403074194221499999977424460241888261147/126812491621022165092163677308365336127830599682204216711429022109284058607666372800652948683672856076106132814795031434807177478522000000
```

The immutable [calibration record](../../artifacts/generated-results/elliptic-curves/det1092_norm8_seed_cover_v2/calibration.json)
identifies the sole exceptional input. No later point or V3 artifact is used.
Both fixed quartics are squarefree and coprime to the original elliptic
discriminant, so both are smooth genus-one covers unramified over all24
nodal fibres.

## An explicit infinite rational-point source

**New explicit application.** For either member, use its inherited point
`(t0,s)`, and expand

\[
F_z(t_0+u)=s^2+q_1u+q_2u^2+q_3u^3+q_4u^4.
\]

Its pointed Weierstrass model is

\[
J_0:Y^2=X^3+q_2X^2+(q_1q_3-4s^2q_4)X
       +s^2q_3^2+q_1^2q_4-4s^2q_2q_4.
\]

The origin maps to `(t0,s)`. For other points the inverse map is

\[
\boxed{u=\frac{2sY+q_1X+2s^2q_3}{X^2-4s^2q_4},\quad
t=t_0+u,\quad W=\frac{Xu^2}{2s}-s-\frac{q_1u}{2s}.}
\]

The forward formulas on `u!=0` are

\[
X=\frac{2s(W+s)+q_1u}{u^2},\qquad
Y=\frac{(X^2-4s^2q_4)u-q_1X-2s^2q_3}{2s}.
\]

The other inherited point `(t0,-s)` extends to

\[
\overline P=(\overline X,\overline Y),\qquad
\overline X=\frac{q_1^2}{4s^2}-q_2,\quad
\overline Y=-\frac{q_1\overline X+2s^2q_3}{2s}.
\]

**Verified arithmetic.** This point is nontorsion on both members. The
checker verifies the birational identities and uses elementary complete
finite group enumeration on the rationally isomorphic short model.
For member00, the group orders at47 and53 are61 and68, with gcd1.
For the calibrated member, the orders at47,53,61 are50,60,63, again
with gcd1. Good odd reduction therefore forces rational torsion to be
trivial. The displayed finite point is not the origin.

All nonzero multiples `n*Pbar` consequently give infinitely many rational
points on `C_z` through the explicit inverse map. On both fixed curves
`q4` is positive and nonsquare, so the denominator `X^2-4s^2q4` cannot
vanish at a rational finite point. These are elliptic-group multiples,
not a rational parametrization by `P1`. Their coefficients and exact
generator are exported in the
[independent replay](../../artifacts/generated-results/elliptic-curves/det1092_norm8_seed_cover_v2/replay.json).

## Rank18 after pullback

**Established literature applied to exact geometry.** We use the
intersection height pairing and base-change scaling from
[Schuett--Shioda, sections5 and11](https://arxiv.org/abs/0907.0298).

**New deduction.** The primitive RR line has coefficient degree bounds
`(8,4,2)`. Its full divisor is `3O+8F`, with no vertical component, including
at infinity. Removing the known section `-Z` gives

\[
D=2O+4F+\phi(w),\qquad D^2=0,\qquad D.O=0.
\]

The squarefree quartic normalization has genus one, equal to `p_a(D)`.
The pulled-back surface has48 `I1` fibres and `chi=4`; the extra section
has zero intersection with the zero section and height8. Its trace is
`Z`, and deck symmetry gives the full height Gram

\[
\mathcal G=\begin{pmatrix}2G&Gw\\w^tG&8\end{pmatrix},\qquad
\operatorname{Schur}(\mathcal G)=4,
\]
\[
\boxed{\det\mathcal G=4\cdot2^{17}\cdot1092=572522496>0.}
\]

Hence both pulled-back families have rank at least18 over `Q(C_z)`.
The anti-invariant section has height16. There is one new quotient
direction; its conjugate is not a second independent one.

**Established specialization theorem, verified application.** The original
`j(t)` is nonconstant. The specialization theorem for nonconstant elliptic
families over a curve therefore makes the certified subgroup injective
outside finitely many rational base points; see
[Silverman's theorem, pages3--4 of his presentation](https://www.math.brown.edu/johsilve/Presentations/JA09Talk.pdf).
Combined with the explicit infinite point source and degree-two map, this
gives infinitely many distinct rational original parameters with rank at
least18. The exceptional set or an effective bound was **not computed**.
Unlike the separate conic progression theorem, this is not a uniform rank
certificate for every integer multiple of the base generator.

## Exact panel replay and remaining obstruction

**Verified application.** On the unchanged eight controls and302, member00
is nonsplit at every address. Member01 is split precisely at302. Exact
integer square-root inequalities retain all18 cover/address outcomes.
Every nonsplit result concerns that cover above that address, not absence
of exceptional elliptic points or rational points elsewhere on the cover.

The [equation-only member01 payload](../../artifacts/generated-results/elliptic-curves/det1092_norm8_seed_cover_v2/equation-only-cover-01.json)
contains no exceptional point coordinates or control outcomes. The separate
[equation-only replay](../../artifacts/generated-results/elliptic-curves/det1092_norm8_seed_cover_v2/equation-only-replay.json)
extracts the square at zero, maps it back, and certifies generic rank17 and
rank18 with its image using complete finite `E(Fp)/2E(Fp)` quotients.
A root-free cubic modulo31 excludes rational2-torsion. Only afterward does
the independent geometric checker compare with the historical first point,
up to the inherited centre-complement.

**Exact chart exclusions.** The affine maps exclude zeros of `h`, the
elliptic discriminant, the branch polynomial and every denominator in the
stored maps, parent equation and generic section coordinates. The worker
derives this finite list directly from those equations and checks all nine
addresses. Proper smooth maps may extend beyond this affine patch; no
undefined specialization is assigned a rank.

**Subsequent verified application.** The
[full Picard-image comparison](DET1092_GENUS_ONE_PICARD_SPECIFICITY_2026-09-08.md)
now proves inherited rank12 and rank13 with the first class. Thus the
inherited rational-point source cannot generate either point above302.
However, nine generic-old-point controls also give carrier rank12 to13
while remaining generic on the original elliptic fibres; this condition
fails specificity.

**Unresolved.** Choosing `z_star` from generic data without the first point
remains open. No full Selmer dimension, strict-class discriminator or
nonzero Sha event is inferred here.
The genus-zero possibility through this particular seed remains unknown.
The earlier genus-two minimum in its fixed norm-ten pencil is unchanged.

## Reproduction and resource record

**Verified computation.** Each stage finishes within two seconds under a
25-second cap. The original failed attempt is retained, and no selection
or budget was retuned. No point search, original-parameter sweep, later
point, V3 input, class group, paid backend, pilot mutation or detached job
was used. The two positive-rank auxiliary bases are explicit constructions;
only the nine original addresses receive the comparison/rank replay.

```sh
sage -python research/elliptic-curves/cas/verify_det1092_norm8_first_unlock_cover.sage
sage -python research/elliptic-curves/cas/replay_det1092_genus_one_seed_equations.sage
```
