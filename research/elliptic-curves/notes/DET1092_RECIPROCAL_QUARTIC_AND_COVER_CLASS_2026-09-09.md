# Which Jacobian class can actually detect a seed?

## Result and literature boundary

**Verified application of established binary-quartic descent.** For any
smooth elliptic curve `E:Y^2=X^3+aX+b` in characteristic zero and inherited
centre `Z=(c,d)`, the pointed slope quartic

\[
\mathscr D_Z:\quad w^2=q(m)=m^4-6cm^2-8dm-3c^2-4a
\]

is isomorphic to `E` with a rational origin. Its conventional degree-four
covering map, with the signs used here, is

\[
\boxed{\Phi(P)=Z-2P.}
\tag{1}
\]

Consequently its class as a2-covering is **`delta(Z)`**, in the inherited
Kummer subgroup, and its Sha image is zero. This is not the Kummer class
of a newly found point `P`. See
[Cremona--Fisher, section6 and Theorem13](https://johncremona.github.io/papers/quartequiv.pdf)
for the general covering-class correspondence. No literature novelty is
claimed for that correspondence.

**New programme obstruction, verified generically.** Transposing the
existing norm-eight pencil to fix the original parameter `t` and vary its
member label `z` gives precisely such a quartic. Its Jacobian is the
original `E_t`, not a smaller auxiliary arithmetic problem. Its covering
class is inherited on every smooth fibre in the chart, including302 and
the null controls. Testing whether this class survives descent, has zero
Sha image, or lies outside the generic image cannot distinguish a seed:
the first two hold automatically and the last never holds.

This does **not** rule out useful chart reduction, marked-point descent,
other Selmer classes, or the existing constructive fixed-member base changes.

## Three objects which must not be conflated

**Verified distinction.** The same marked surface point lies on two
different curves in the pencil, with different origins and Jacobians:

| Object | Fix / vary | Relevant class | What is established |
|---|---|---|---|
| Carrier `C_z` | Fix `z`, vary `t` | `[A-B]` in `Jac(C_z)` | First302 gives inherited12 to13; the nine old-point controls do too. |
| Reciprocal fibre `D_t` | Fix `t`, vary `z` | `[A-infinity+]` in `Jac(D_t)=E_t` | Exactly the original elliptic point; this is new iff that point is outside the original generic rational span. |
| Degree-four covering of `D_t` | Retain the covering map as part of the object | `delta(Z_t)` in `H^1(Q,E_t[2])` | Always inherited, with zero Sha image, independently of which rational point is marked. |

The earlier
[Picard-image specificity obstruction](DET1092_GENUS_ONE_PICARD_SPECIFICITY_2026-09-08.md)
and
[nonisogeny/torsion-transfer obstruction](DET1092_GENUS_ONE_TRANSPORT_OBSTRUCTION_2026-09-08.md)
concern the **first row**. They are unchanged. In particular, the new
identity in the second row supplies no map from that fixed carrier Jacobian
to the marked original elliptic fibre.

## Exact universal identities and the sign of the covering map

**Verified algebra.** The curve isomorphism on the affine patch is

\[
X=(w-c+m^2)/2,\qquad Y=m(X-c)-d,
\]
\[
m=(Y+d)/(X-c),\qquad w=2X+c-m^2.
\tag{2}
\]

The two rational points at infinity, `w/m^2=+1,-1`, extend to `O,Z`.
For the first, `X` has a double pole; for the second,
`w=-m^2+3c+4d/m+O(m^{-2})`, giving `(X,Y)->(c,d)`.
The quartic involution `w->-w` is therefore `P->Z-P`. In particular,

\[
\boxed{w=0\quad\Longleftrightarrow\quad 2P=Z.}
\tag{3}
\]

**New consequence for seed incidence.** Any rational branch point is in
the rational span of the inherited centre, even when it enlarges an
integral subgroup by index two. It cannot be a rank-increasing seed.
Rational points away from the branch locus can be either inherited or new.

**Verified invariant identities.** In the conventional binary-quartic
normalization,

\[
I=-48a,\qquad J=-1728b,\qquad \operatorname{disc}(q)=256\Delta(E).
\tag{4}
\]

Thus `Jac(D_Z):v^2=u^3-27Iu-27J` is carried to `E` by
`X=u/36`, `Y=v/216`, over the original ground field.

For an independently checkable covering map put

\[
G=cm^4+4dm^3+(6c^2+4a)m^2+4cdm+c^3+4b,
\quad H=(q'G-qG')/2.
\]

The verifier proves the polynomial identity

\[
H^2=4G^3+4aGq^2+4bq^3.
\]

Hence the covariant map is

\[
\Phi(m,w)=\left(G/q,\ H/(2qw)\right).
\]

**Verified identification.** Formula(2) pulls `dX/(2Y)` to `dm/w`;
the covariant map pulls it to `-2dm/w`. Its value at positive infinity is
`Z`, since `G` has leading coefficient `c` and `H` has degree six with
leading coefficient `2d`. Translate its image by `-Z`. A pointed morphism
between elliptic curves is a homomorphism, and in characteristic zero
the differential determines that homomorphism. It is therefore `[-2]`,
proving(1). Equivalently the degree-four map is a translate of doubling,
not the degree-one isomorphism(2).

The rational point above `O` makes the covering class a rational Kummer
class. For every rational marked point,

\[
\Phi(P)-Z=2(-P),\qquad\delta(\Phi(P))=\delta(Z).
\]

This conclusion holds regardless of the rank of the original fibre.

## Generic determinant1092 application

**Verified generic application, no exceptional input.** Reuse the
norm-eight centre `Z=P14-P15` and its generic polynomials `h,nx,ny,shift`.
With `c=nx/h^2`, `d=ny/h^3`, write the existing pencil as

\[
C_z:\quad W^2=F_z(t),\qquad
m=h(t)z-\mathrm{shift}(t)/h(t),\quad w=h(t)W.
\]

Fixing `t` gives

\[
F_z(t)=h(t)^{-2}q\bigl(h(t)z-\mathrm{shift}(t)/h(t)\bigr).
\]

The affine substitution and ordinate scaling cancel in the binary-quartic
invariants. The checker verifies over **Q(t)**, not just at a few fibres,

\[
I_z(F)=-48a(t),\qquad J_z(F)=-1728b(t),\qquad
[z^4]F=h(t)^2.
\]

Thus the two reciprocal infinities are always rational. Choosing `z` so
that the member splits above a prescribed `t` is exactly finding a rational
point on this alternative model of `E_t`. Relative to positive infinity,
the Abel--Jacobi class is exactly `P`; subtracting another inherited
basepoint changes it only by an inherited point.

**Patch and scope.** Require `h(t)!=0`, `Delta(E_t)!=0` and the generic
coefficient denominators to be defined. Affine formulas also omit their
displayed denominators; the smooth projective curve isomorphism extends
over removable chart exceptions. Degenerate fibres are not assigned a
Jacobian or rank by these identities.

## Independent old-control replay

**Verified retrospective application.** After freezing the universal and
generic proof, a separate script reconstructs covariants via the quartic
Hessian, then compares them with Sage elliptic group operations. It does
not import the universal checker. Its32 points are:

- the two old split branches on each of302 and four old conic fibres;
- one inherited section on each of those five fibres;
- the nine completed blinded generic reconstructions, assessed over full17;
- one inherited section on each of the eight unchanged V3/V4-null fibres.

Eight displayed branches are already certified independent over full17;
24 points are inherited, including both branches of the dependent conic.
The nine masked recoveries are still new over their M16 cores, but old
over full17. The replay checks `Phi(P)=Z-2P`, the exact half `-P` of
`Phi(P)-Z`, and `P+conjugate(P)=Z` in **all32 cases**.

The common result is inherited covering class and zero Sha image, despite
the different marked-point rank outcomes. These are32 diagnostic records,
not32 different points claimed independent. No later rank21 or302 point
was used; the rank21 control uses only its previously certified first seed.
The original independence decisions remain bound to the separate
[halving-or-cycle certificates](DET1092_SPLIT_SEED_HALVING_DICHOTOMY_2026-09-08.md).

## What remains constructive, and what is now excluded

**Verified application.** The conic progression and the new
[two-fibration rational covers](DET1092_TWO_FIBRATION_SEED_CONSTRUCTION_2026-09-08.md)
still construct independent sections after nonconstant base change. Those
are global curves crossing many original fibres, not merely new models
of a prescribed original fibre.

**New programme conclusion.** Using the known-centre quartic's covering
class as a putatively nongeneric Selmer/Sha event is excluded for this
entire class of charts. Using the reciprocal quartic's full Jacobian
descent is instead descent on the original elliptic curve. The distinction
between old and new lives in the **marked rational point's class**, and
its divisibility relative to the generic subgroup, not in that covering
label. This extends the earlier trace/twist observation to the actual
pointed-quartic covariant map, including its normalization and origin.

**Unresolved.** A generic-only selection of a cover that produces a302
seed is still open. The above identification is an obstruction to one
proposed auxiliary discriminator, not a solution of rational-point
construction on every fibre or an explanation of subsequent amplification.

## Replay and resources

**Verified resources.** Both commands complete in under one second under
25-second caps, using Sage10.9. No point search, parameter sweep, class
group, Selmer campaign, production change, or background job is started.

```text
timeout 25s sage -python research/elliptic-curves/cas/verify_det1092_reciprocal_quartic.sage
timeout 25s sage -python research/elliptic-curves/cas/replay_det1092_reciprocal_quartic_controls.sage
```

The [generic proof](../../artifacts/generated-results/elliptic-curves/det1092_reciprocal_quartic_v1/generic.json),
[input protocol](../../artifacts/generated-results/elliptic-curves/det1092_reciprocal_quartic_v1/protocol.json)
and [control replay](../../artifacts/generated-results/elliptic-curves/det1092_reciprocal_quartic_v1/control-replay.json)
are immutable and hash-bound.
