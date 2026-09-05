# Shared quadric surfaces, distinct ruling base changes

All twenty inherited anchor Kummer classes have the same rational
quadric total space, up to an explicit Q-isomorphism:

\[
 XY=Z^2-BW^2.
\]

Their rulings are defined over `F=Q(sqrt(B))`. This is the same field
that supplies the fixed-cubic pencil's single geometric generic
Mordell–Weil direction. But the forty ruling lines through the anchor
points give forty **distinct** quadratic base changes. The common
quadric does not make their level curves soluble together.

The [protocol](QUADRIC_RULING_PROTOCOL.json) and
[certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_quadric_rulings_v1.json)
retain 23 determinant checks, twenty explicit isomorphisms, forty
section equations and their branch values.

    python3 elliptic-curves/rank-jump/quadric_rulings.py check
    sage -python elliptic-curves/rank-jump/quadric_rulings.py verify

## The total space of a fixed Kummer class

Write `z=a+bθ+cθ²` and let `Q_j(z)` be the coefficient of `θ^j`
in `βz²`. The projective two-cover of `E_u` is

\[
 Q_1(z)+w^2=0,\qquad Q_2(z)+uw^2=0.
\]

On `w≠0` its rational image has `x=Q_0(z)/w²`. The first equation
defines a fixed quadric surface; `u=-Q_2/w²` is a rational function
on that surface, and its level curves carry the remaining solubility
problem.

The coefficient pairing `(r,s)↦[θ](rs)` has determinant `B` in the
basis `1,θ,θ²`. Multiplication by `β` has determinant `N(β)`.
Consequently

\[
 \det(Q_1)=B\,N(\beta).
\]

For a norm-square Kummer class the four-variable quadric has
discriminant squareclass `B`. The certificate checks this for all
twenty anchor basis classes and the three previously fixed CT chain
classes. It requires no full class group or new descent.

For `β=p-θ` with `P=(p,q)∈E_0(Q)`, put `H=p²+A` and
`t=a-pb-Ac`. Then

\[
 Q_1=-t^2+Hb^2+2Bbc-pBc^2,
\qquad
 H(Hb^2+2Bbc-pBc^2)=(Hb+Bc)^2-Bq^2c^2.
\]

The invertible linear change

\[
 X=H(t+w),\quad Y=t-w,\quad Z=Hb+Bc,\quad W=qc
\]

gives the displayed common quadric. The retained points all have
`Hq≠0`. The anchor `(z,w)=(1,1)` maps to `[1:0:0:0]`.
Projection from this rational point makes the surface rational over Q.
Its two lines through that point require `sqrt(B)`; the fixed
nonsquare `B` prevents either ruling from being defined over Q.

## The ruling lines and their branch values

For `ε=±1` set

\[
 b_\epsilon=\frac{-B+\epsilon q\sqrt B}{p^2+A},\qquad
 v_\epsilon=(pb_\epsilon+A)+b_\epsilon\theta+\theta^2.
\]

Then `[θ](βv_ε)=0` and `[θ](βv_ε²)=0`. The line
`z=1+s v_ε,w=1` stays on the quadric. It gives

\[
 u(s)=-2[\theta^2](\beta v_\epsilon)s
          -[\theta^2](\beta v_\epsilon^2)s^2,
\]
\[
 x(s)=[1](\beta(1+sv_\epsilon)^2),\quad
 y(s)=q\,N(1+sv_\epsilon).
\]

Here `[1]` denotes the constant coefficient. The script verifies the
curve identity coefficientwise for all forty lines over `F[s]`.

The quadratic map has branches at infinity and

\[
 \boxed{
 u_\epsilon(P)=
 \frac{p^2}{Ap+2B+2\epsilon q\sqrt B}
 =\frac1{x(P-\epsilon T)},\qquad T=(0,\sqrt B)\in E_0(F).
 }
\]

Thus the branch values are controlled by translations on one common
auxiliary elliptic curve, rather than forty unrelated formulas. Both
expressions agree exactly in the certificate; all forty finite values
are distinct and avoid `D=0`.

There is also a structural reason for distinctness. Equality of two
reciprocal x-coordinates forces
`P_i-εT=±(P_j-ηT)`. Unless `P_i=±P_j`, this would make `2T`
rational. Conjugation sends `T` to `-T`, so that would force `4T=0`.
But `y(2T)=-(A³+8B²)/(8B sqrt(B))≠0` for the pinned anchor.
The independent anchor points have no equal or opposite pairs.

## What the constructions prove

Each quadratic base change over F has **exact generic rank two over
F(s)**. The original generic section gives one direction. The displayed
section is not invariant under the quadratic deck transformation:
its other value above `u=0` differs from the anchor point, as checked
exactly. Its anti-invariant difference is nonzero and vanishes at the
finite branch point, which is a good fibre. A nonzero torsion section
cannot specialize to zero at a good characteristic-zero fibre.
The anti-invariant direction is therefore non-torsion and independent.

For the upper bound, the original `3I_2+I_0^*` configuration becomes
six `I_2` fibres. Ramification at infinity removes `I_0^*`, and the
other branch is smooth. The new minimal elliptic surface over `P¹`
has Euler number 12 and is geometrically rational. Shioda–Tate gives
`10-2-6=2`. Thus the two explicit directions attain the upper bound.
This is a function-field result over F; it gives no new Q-curve rank.

Choosing one ruling for each of the twenty classes gives twenty private
finite branch points and one shared point at infinity. The connected
compositum therefore has genus `1+2^18(21-4)=4,456,449`.
Shared branching reduces the genus but does not yield a low-genus
twenty-direction base.

These are **incidence constructions over F**. Shared quadric rationality
is only a prerequisite for **solubility** of its fibres. In particular,
the earlier CT-obstructed Q-fibres lie on these same fixed total spaces;
their obstruction is not removed by exhibiting a rational surface.

The stronger Q-rational construction found while analyzing these
quadrics is the [linear-twist soluble block](LINEAR_TWIST_SOLUBLE_BLOCKS.md).
It gives a common affine Kummer correction and a proved rank bound,
and identifies a genus-two component-splitting event at the anchor.

## Superseded first test

The proposed product-character test on the 38 target-fitted historic
quartics did not need a new computation. The pinned
[visibility filtration](../../artifacts/generated-results/elkies-k3-r17-multisection-visibility-filtration-v1.json)
already records factor degrees `[4]` for every quartic, and its source
checker requires irreducibility. None can be a product of the atlas's
rational quadratic branch polynomials. This excludes those
presentations only; an exceptional point need not have a unique
cover through it.
