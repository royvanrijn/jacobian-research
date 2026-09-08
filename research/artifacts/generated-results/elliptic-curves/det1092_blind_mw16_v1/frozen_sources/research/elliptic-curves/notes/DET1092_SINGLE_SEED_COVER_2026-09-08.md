# A degree-two cover certifying one determinant-1092 seed

## Strongest result and its boundary

**Verified application and new deduction.** The existing genus-two RR
member through the historical first302 unlock is now an explicit base-change
construction of an independent eighteenth **elliptic** section. Over its
function field the section has height10, its anti-invariant difference has
height20, and the inherited-plus-new height matrix has Schur complement5
and determinant715653120. All equations and maps are over **Q**.

**Verified application.** A frozen, equation-only execution at the nine
existing parameter addresses constructs and independently certifies an
eighteenth point at302, and proves this cover **nonsplit at every one of the
eight controls**. The first historical point's coordinates are read only in
a subsequent comparison, after execution, which confirms exact agreement.
No V3 search artifact or later exceptional point is an input.

**Important limit.** The cover itself was selected using the historical
unlock. Its302 success is therefore calibrated, **not a held-out discovery**.
The second, already source-fixed RR member `B` also gives rank at least18
over its function field, but is nonsplit at all nine addresses. The control
comparison proves exact facts about these two covers; it does not validate
a general rank-incidence predictor or prove the controls have rank exactly17.

The [frozen protocol](../../artifacts/generated-results/elliptic-curves/det1092_single_seed_covers_v1/protocol.json),
[expanded successful cover and maps](../../artifacts/generated-results/elliptic-curves/det1092_single_seed_covers_v1/cover-01-input.json),
and [independent panel replay](../../artifacts/generated-results/elliptic-curves/det1092_single_seed_covers_v1/panel-replay.json)
are authoritative numerical evidence. No point search or parameter expansion
was run. The seeded amplifier remains outside this work.

## Audit: which existing results are being reused?

**Verified applications already established, not re-proved as new results.**

- The parent has full geometric and arithmetic MW rank17, with height
  matrix `G` of determinant1092, and24 nodal fibres.
- The [first RR net](DET1092_FIRST_UNLOCK_RR_NET_2026-09-08.md) gives exact
  generic polynomials `A,B` and the whole historical chord chart. The
  successful coordinate was `-1714/2373`; the first transition was17 to18.
- The [Jacobian certificate](DET1092_FIRST_UNLOCK_JACOBIAN_CLASS_2026-09-08.md)
  proves `xi=[P_unlock-P0]` independent of the full rank17 Picard restriction
  image. Its true Kummer class is rational, hence has zero Sha image.
- The [specificity controls](DET1092_RR_GENERIC_POINT_SPECIFICITY_2026-09-08.md)
  prove exactly the same Jacobian independence for generic elliptic points
  on all eight null fibres. Thus this Jacobian property alone does not
  certify an amplifier seed.
- The [marked transport gate](DET1092_MARKED_KUMMER_TRANSPORT_2026-09-08.md)
  identifies the actual elliptic cubic Kummer class once a rational branch
  is resolved. It excludes the proposed sextic-to-cubic norm and ordinary
  Jacobian/2-torsion-module transfer routes, not nonlinear incidence.
- The [whole witness-pencil gate](DET1092_FIRST_WITNESS_PENCIL_GENUS_GATE_2026-09-08.md)
  excludes every rational genus-zero or genus-one member in this particular
  RR pencil, not in every divisor system.
- The [orbit8044 conic](DET1092_RANK18_BASE_CHANGE_AND_INITIAL_UNLOCK_2026-09-08.md)
  already gives a rationally parametrized rank18 base change, with infinitely
  many rational base points. It is nonsplit above302 and does not explain
  this seed. The present genus-two result does not supersede that construction.

The new work is the explicit elliptic height certificate for the two existing
RR sextics, a witness-free execution of their point maps, and exact splitting
obstructions at the unchanged eight control addresses.

## Explicit one-dimensional auxiliary cover

**Verified application.** Use the historical generic centre `Z=P_w`, where

```text
w=(1,-1,-2,1,0,-1,1,-1,1,-1,1,0,1,0,0,-1,0),  w^t G w=10.
```

On the literal parent `y^2+a1*x*y+a3*y=x^3+a2*x^2+a4*x+a6`, put

\[
c_x=x_Z+b_2/12,\qquad c_y=y_Z+(a_1x_Z+a_3)/2,
\qquad a=-c_4/48,\qquad h^2=\operatorname{den}(x_Z).
\]

For either frozen value of `u`, set `f_i=B_i+u A_i`, with `v=0`, and

\[
m=-f_1/f_2+a_1/2,\quad
F_u(t)=\frac{f_2^4}{h^6}
 \left(m^4-6c_xm^2-8c_ym-3c_x^2-4a\right).
\]

The checker proves that `F_u` is an exact rational polynomial of degree6,
squarefree and coprime to the elliptic discriminant. No constant squareclass
is discarded. The auxiliary curve and parameter map are

\[
\boxed{C_u:s^2=F_u(t),\qquad C_u\longrightarrow\mathbf P^1_t,\ (t,s)\mapsto t.}
\]

Its elliptic point is given by

\[
X=\frac{m^2-c_x+h^3s/f_2^2}{2},\qquad
Y=m(X-c_x)-c_y,
\]
\[
\boxed{x=X-b_2/12,\qquad y=Y-(a_1x+a_3)/2.}
\]

Every coefficient of `F_u`, and the four rational functions in the linear
expressions `x=x0+x1*s, y=y0+y1*s`, is expanded in the equation-only payload.
The checker independently substitutes them into the original Weierstrass
equation modulo `s^2-F_u`, and verifies the residual quadratic elimination.
The deck involution satisfies

\[
P+\sigma(P)=Z,\qquad \sigma(t,s)=(t,-s).
\]

**Frozen choices.** Cover00 is the already fixed `B`, namely `u=0`.
Cover01 uses the previously recorded value

```text
u0=-24866539961702974106913661426507481835384195926726682303377531616158216282536772141507729874986070411728576375513733832475829768144420211347951753983533120422763938994630741972201949175154930592306424929883241467002148851014200/139537.
```

No new member is selected. Cover01 is exactly the earlier curve carrying
the certified Jacobian class `xi`, including its rational twist factor.

**Excluded parameters.** The displayed affine certifier excludes zeros of
the elliptic discriminant, `F_u,f2,h`, and every denominator in the point
maps, parent equation and17 generic sections. The input retains the four
polynomial factors; the checker constructs the denominator list directly
from the equations. None of the nine tested addresses is excluded.
These are chart/certifier exclusions, not a claim that every omitted smooth
point fails to extend on the proper curves.

## Independence over the cover: a geometric certificate

**Established literature.** Use the elliptic-surface intersection height
formula, degree scaling under base change, and positivity modulo torsion
from [Schuett--Shioda, sections5 and11](https://arxiv.org/abs/0907.0298).
The six simple branch points give genus2 by
[Riemann--Hurwitz](https://stacks.math.columbia.edu/tag/0C1B).

**New deduction applied to verified equations.** The line is primitive, has
coefficient degrees `(10,6,4)`, and has no vertical fibre component, including
at infinity. Removing its trace section `-Z` gives

\[
D=3O+10F-(-Z)=2O+5F+\phi(w),\quad D^2=2,\quad D.O=1.
\]

Squarefree sextic normalization has genus2, equal to `p_a(D)`, so the
irreducible residual curve is smooth. The base change is unramified over
all24 nodal fibres. Its elliptic surface has48 `I1` fibres and `chi=4`.
Projection formula gives `P.O=1`, hence

\[
\langle P,P\rangle=2\chi+2P.O=10.
\]

The inherited trace has height20 after pullback. Deck symmetry and
`P+sigma(P)=Z` give the complete height Gram

\[
\mathcal G=
\begin{pmatrix}2G&Gw\\w^tG&10\end{pmatrix},\qquad
10-(Gw)^t(2G)^{-1}(Gw)=5.
\]

Thus `det(Gcal)=5*2^17*1092=715653120>0`, certifying rank at least18 over
`Q(C_u)` for **both** existing covers. In particular

\[
R=P-\sigma(P)=2P-Z,\qquad \langle R,R\rangle=20,
\]

is a nonzero anti-invariant direction, orthogonal to the inherited group.
This is an elliptic function-field result, not a genus-two Jacobian rank.
It requires no number-field extension of constants.

**New deduction: minimality with a stated scope.** A degree-one cover of
the parameter line cannot add a section, since the full original generic
rank is17. Degree2 is therefore the smallest possible base degree for an
extra direction. Genus2 is the minimum in the *fixed next RR pencil through
this witness*, by the existing modulo191 genus gate. Global minimum genus
among other centres, translates or divisor systems through this witness is
**UNKNOWN**; the independent orbit8044 conic shows genus0 is possible for
other seeds on this parent.

## The actual seed criterion, distinct from curve solubility

**Verified application / exact sufficient criterion.** For a supplied
rational address `t` in the displayed open set:

1. Test exactly whether `F_u(t)` is a rational square. A nonsquare excludes
   this degree-two cover above that address; no elliptic rank conclusion follows.
2. If it is square, extract its rational square root and construct `P` by
   the displayed map. This is direct equation evaluation, not a point search.
3. Certify independence of the17 generic sections together with `P`. The
   replay uses exact finite groups `E(F_p)/2E(F_p)` at the18 inherited fixed
   primes through197. Binary column rank18, together with a root-free
   2-division cubic modulo31, proves independence by infinite descent.
   Failure of this finite certificate remains **UNKNOWN**, not dependence.

In symbols, with `M_t` the elliptic generic specialization subgroup,

\[
F_u(t)\in\mathbf Q^{\times2}\quad\text{and a verified rank18 certificate}
\quad\Longrightarrow\quad P_t\notin M_t\otimes\mathbf Q.
\]

**New deduction / necessary distinction.** Splitting alone is insufficient.
Every one of these curves already has an inherited rational pair over the
unique root of `f2/h`. Its two elliptic images, after extending the affine
map, are `O` and `Z`. Generic positive height therefore cannot replace the
specialized independence check. Also `R=2P-Z` has inherited elliptic
2-Kummer image `delta(Z)`; use `P`, not `R`, in the mod2 independence test.

For the RR Jacobian, the two marked branch classes satisfy
`xi(P)+xi(sigma(P))=[K_C-2P0]`, an inherited class. For the elliptic fibre,
`P+sigma(P)=Z` is inherited. Thus either branch gives the same non-generic
direction up to sign in the respective rational quotient. Neither fact
identifies the two ambient groups or induces a Jacobian-to-elliptic homomorphism.

## Exact unchanged control comparison

**Verified application.** All nine original addresses are copied unchanged
from the previous generic-input control protocol. Small-prime obstructions
below certify failure of the square condition at a **specified address**;
they do not say the whole genus-two curve is locally insoluble. Each curve
has an inherited rational point elsewhere.

In the table `(p,v,r)` records `v=v_p(F_u(t))` and the residue `r` of its
`p`-adic unit part. Odd `v`, or a nonsquare `r` when `v` is even, is exact.

| Address | Source-fixed `B` | Calibrated first-unlock member |
|---|---|---|
| scale-0131232 | `(3,11,1)` | `(3,11,1)` |
| scale-0257585 | `(17,-6,7)` | `(17,-6,7)` |
| scale-0487239 | `(3,13,1)` | `(3,13,1)` |
| scale-0177036 | `(3,8,2)` | `(3,8,2)` |
| scale-0043332 | `(37,0,24)` | `(29,22,8)` |
| scale-0590501 | `(3,-3,2)` | `(3,-3,2)` |
| scale-0290097 | `(3,-3,2)` | `(3,-3,2)` |
| scale-0748009 | `(3,3,1)` | `(3,3,1)` |
| 302, `t=0` | `(5,64,3)` | rational split; independently certified17 to18 |

The cover01 execution reads no exceptional coordinates and reconstructs
exactly the historical replay's elliptic point. This identifies the seed;
it does **not** make the cover's historically fitted coefficient `u0`
independent of that seed. The eight controls were never used to retune `u`.

## Whole-pencil calibration gate

**New deduction, independently verified application.** The eight fixed-cover
obstructions do **not** extend to the whole already selected pencil
`B+(u0+v*t)A`. At every control address its rational points include an
explicit point obtained entirely from generic section0. At302 the affine
equation instead becomes the same nonzero constant square for every finite
`v`, because this pencil was fitted through that marked fibre point.

This distinction has an exact invariant formula. Define the nonzero
**generic-input constant**

\[
\kappa=\frac{A_2B_1-A_1B_2}{h^3}\in\mathbf Q^\times.
\]

The new checker proves that this quotient is constant; it is independent of
the fitted `u0`. For this paragraph write `F(t,v)=F_(u0,v)(t)`, where the
same branch polynomial is degree6 in the parameter-line coordinate `t`
and degree4 in the pencil coordinate `v`. With
`m_num=-f1+a1*f2/2`, the determinant of the slope substitution is

\[
\det\begin{pmatrix}[v]m_{\rm num}&[1]m_{\rm num}\\
                         [v]f_2&[1]f_2\end{pmatrix}
=\kappa t h^3.
\]

**Verified polynomial identities.** In the conventional binary-quartic
normalization, exact coefficient and universal covariance checks give

\[
\boxed{I(F(t,v))=(\kappa t)^4c_4(E_t),\qquad
J(F(t,v))=2(\kappa t)^6c_6(E_t),}
\]
\[
\boxed{\operatorname{disc}_v F(t,v)
=256(\kappa t)^{12}\Delta(E_t).}
\]

The independent checker first reconstructs the quartic by **literal
Weierstrass elimination**, not the constructor's short-coordinate formula.
It then proves the binary-quartic covariance law over a symbolic polynomial
ring and the discriminant identity by a separate universal resultant.
No interpolation or numerical height calculation enters these identities.

**New deduction.** For `t!=0` in the indicated smooth chart, the slope
transformation is invertible. The curve `s^2=F(t,v)`, with `v` now its
coordinate, is the pointed genus-one chord curve, birational over Q to
the **original elliptic fibre** `E_t`. Its Jacobian short equation
`y^2=x^3-27I*x-27J` is scaled from the parent's short equation by
`6*kappa*t`. This is not a second genus-two Jacobian or a new arithmetic
cover class to compare with `E_t`.

At `t=0`, the raw quartic invariants vanish to orders4 and6 and its
discriminant to order12, while `Delta(E_0)!=0`. Dividing out the displayed
Weierstrass scaling recovers the smooth original elliptic model. The
constant-square affine equation is therefore a **pencil-chart collapse**,
not evidence of bad elliptic reduction or a special Selmer dimension.
The projective pencil member at infinity is the previously classified
`C_min+F_0`; its vertical component must not be counted as a new point source.

For an arbitrary marked address `tau0`, the same pencil construction
`B+(u0+v*(t-tau0))A` replaces `t` by `t-tau0` in the determinant and scaling.
Thus the factor `(t-tau0)^12` is forced by the choice of a pencil through
that fibre, not by a rank-jump theorem at the address.

**Generic-only counterwitnesses on the controls.** Let the already frozen
section0 define `r0(t)=-B(S0(t))/A(S0(t))`. Then

\[
v_{\rm gen}(t)=\frac{r_0(t)-u_0}{t},\qquad
s_{\rm gen}(t)=
\frac{(2x_{S_0}^{\rm short}+c_x)f_2(t,v_{\rm gen})^2
                -m_{\rm num}(t,v_{\rm gen})^2}{h^3}
\]

are explicit rational functions. Their complete expressions are in the
[pencil certificate](../../artifacts/generated-results/elliptic-curves/det1092_seed_pencil_gate_v1/construction.json).
The map `t -> v_gen(t)` has degree7 and a simple pole at302. Exact elimination
proves `s_gen^2=F(t,v_gen)` over Q(t), and independent evaluation at all eight
controls maps the pair back to **exactly `S0(t)`**. Its class modulo elliptic
MW17 is zero. The whole-pencil rational-point condition is consequently
insufficient even on the same eight addresses that the fixed member excludes.

These are diagnostic section evaluations, not retuned production members,
point searches, or new seed discoveries. The two earlier frozen covers and
their negative exposure remain unchanged. No Jacobian rank/Selmer-dimension
panel or later-direction analysis is added.

The [independent gate replay](../../artifacts/generated-results/elliptic-curves/det1092_seed_pencil_gate_v1/replay.json)
passes under a25-second cap. It took under one second; the separately capped
constructor also completed. Reproduce this supplement with

```sh
sage -python research/elliptic-curves/cas/verify_det1092_seed_pencil_gate.sage
```

## Next obstruction to a useful parameter search

**Established consequence and new deduction.** No nonconstant rational or
genus-one curve can map to either smooth genus-two cover, by
Riemann--Hurwitz. Consequently the known inherited points and302 witness
do not give a rational parametrization or an elliptic group-law generator
for additional addresses on this fixed cover.

The missing step is now precise: supply additional rational points of this
cover whose mapped elliptic points pass independence, or give a generic-input
rule for a different cover with a justified source of such points. Varying
the whole pencil without such a rule restores the original elliptic fibre;
the exact invariant identities and generic control witnesses above now
verify this boundary. A single rational-function formula defined for generic
`t` would be an original generic section and hence lie in rank17; a successful
construction must restrict addresses or use a genuine base change.

**UNKNOWN.** There is no new source of rational cover points, oracle-free
member selector reproducing302, smaller-genus cover through this witness
outside the classified pencil, or newly constructed control seed. The
literal RR-Jacobian criterion outside its Picard image already holds for
generic-point controls and is not sufficient for the amplifier. No claim
about subsequent exceptional directions or their amplification is made.

## Replay and limits

**Verified computation.** Each command took under one second locally,
under an explicit25-second cap. All18 cover/address decisions, unsuccessful
small-prime checks, point-map identities and finite rank rows are retained.
No long-running process, parameter sweep, global Selmer/class-group job,
point search or V3/V4 mutation was started.

```sh
sage -python research/elliptic-curves/cas/verify_det1092_single_seed_panel.sage
```

The [per-cover checker](../cas/verify_det1092_single_seed_covers.sage) is
independent of the [input constructor](../cas/prepare_det1092_single_seed_covers.sage).
The [panel checker](../cas/verify_det1092_single_seed_panel.sage) re-runs both,
checks the source bindings, and only then compares with the historical witness.
