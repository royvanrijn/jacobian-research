# The whole block lifts through Jacobian Selmer and still remains Sha

The two strict classes in the small self-gluing control **do admit ordinary
Jacobian 2-Selmer lifts**. None of those lifts can be rational. Its
Jacobian has exact rank two, 2-Selmer dimension six, and four-dimensional
`Sha[2]`. Thus the missing implication in this control lies after the
Jacobian Selmer condition, not merely after the unrestricted norm equation.

| Curve / Jacobian | Exact rank | `dim Sel2` | `dim Sha[2]` | Strict block admitted by Jacobian Selmer projection | Strict block admitted by rational projection |
|---|---:|---:|---:|---:|---:|
| `C+: Z^2=X^6-14X^4+11X^2+1`, `J+` | 4 | 4 | 0 | 0 | 0 |
| `C-: Z^2=-X^6-14X^4-11X^2+1`, `J-` | 2 | 6 | 4 | 2 | 0 |

These are the same two fixed genus-two controls from the
[nonscalar norm calculation](NONSCALAR_CUP_BLOCK_AND_SELF_GLUING.md).
No new curve, descent, norm equation, class group or point search was run.
The [protocol](SMALL_JACOBIAN_SELMER_PROTOCOL.json) declares the finite
local checks and their confirmatory scope.

## Keep the second-quotient labels

Put `f(T)=T^3-11T^2-14T-1` and
\[
C_s:\ Z^2=f(sX^2-1),\qquad s=\pm1.
\]
The two degree-two maps are
\[
\pi_0(X,Z)=(sX^2-1,Z),\qquad
\pi_s(X,Z)=(X^{-2},sZ/X^3).
\]
Their elliptic targets are `E0: y^2=f(x)` and
`E_s: y^2=x^3+11s*x^2-14x+s`. Thus `E_-=E0` and
`E_+=E0^(-1)`. The sextic discriminant is
`-64*s*163^4`, so these are smooth genus-two curves with no bad prime
outside 2 and 163 in the displayed model.

Let `tau(theta)=theta^2-12theta-2=-1/(theta+1)` be the cubic automorphism.
At the Weierstrass pair above a root `theta`, the second quotient's
2-torsion abscissa is
\[
\frac{s}{\theta+1}=-s\tau(\theta).
\]
Therefore the second quotient's standard cubic squareclasses pull back
by **`tau`**. Replacing this by the identity would give the wrong common
Selmer space for the positive sign.

Use global norm-kernel coordinates `(beta0,beta1,u0,u1)`, with
`U=<beta0,beta1>`, `u0=-1-theta`, `u1=tau(u0)`, and `u0*u1=theta`.
The [exact elliptic calculation](NORM_LIFTS_CAN_BE_ENTIRELY_SHA.md) gives
\[
S_0=U+\langle u_0\rangle,\qquad
S_+^{\rm standard}=U+\langle\theta\rangle.
\]
Since `tau(theta)=u0^(-1)`, the correctly pulled-back second spaces are
\[
S_+^{\rm pulled}=U+\langle u_0\rangle=S_0,\qquad
S_-^{\rm pulled}=U+\langle u_1\rangle.
\]
Their common spaces `S_C` consequently have dimensions three and two.

## Complete local conditions and the summed space

At 2 the cubic algebra is the unramified cubic extension. All three
nonzero classes `u0,u1,u0*u1` are nonsquare there. Exhausting the 512
residues modulo 8 proves this; an independent replay already excludes
square roots among the 64 residues modulo 4. The power basis is integral
and étale at 2, so a local unit square root would reduce to one of these
residues. Thus `u0,u1` are independent local squareclasses.

The elliptic local point image has dimension one at 2. Its nonzero
class is certified by a global point: `u0` for `E0`, and `theta` for
the standard twist. Pullback by `tau` therefore makes the two local
lines equal for `C+`, and distinct for `C-`. The same distinction holds
at infinity, using their real signs. At 163 the irreducible local cubic
has no rational 2-torsion, so the odd-adic Kummer dimension is zero.
At every other prime both local conditions are the unramified image.

Let `C_v` and `D_v` be the intersection and sum of the two local images.
Let `t_v` be the connecting rank from local rational 2-torsion in the
Jacobian kernel sequence. The complete differing-place table is

| Sign | Place | `dim C_v` | `dim D_v` | `t_v` | Ordinary Jacobian Kummer dimension |
|---|---|---:|---:|---:|---:|
| `+` | 2 | 1 | 1 | 0 | 2 |
| `+` | 163 | 0 | 0 | 0 | 0 |
| `+` | infinity | 1 | 1 | 1 | 1 |
| `-` | 2 | 0 | 2 | 0 | 2 |
| `-` | 163 | 0 | 0 | 0 | 0 |
| `-` | infinity | 0 | 2 | 2 | 0 |

The middle dimension is `dim D_v + dim C_v - t_v`; dropping the real
connecting term would overcount. The real dimensions are independently
checked using the even-subset model of `J[2]`: conjugation fixes a
three-dimensional space for `J+` and a two-dimensional space for `J-`.
Subtracting the genus gives real component dimensions one and zero.
Both Jacobians have no rational 2-torsion, also checked in that model.

For the positive sign, all local lines coincide, so `S_D=S_C=S0` has
dimension three. For the negative sign, the summed conditions admit the
entire four-dimensional norm-kernel space `V_S=U+<u0,u1>`. Hence
\[
(\dim S_C,\dim S_D)=(3,3)\text{ for }J_+,
\qquad(2,4)\text{ for }J_-.
\]
The finite checks certify the local dimensions, not particular cocycle
representatives of every Jacobian local Kummer image.

## Ordinary Selmer lifts and the two sources of Sha

Write `Phi:E0 x E_s -> J_s` and `Psi:J_s -> E0 x E_s` for the dual
degree-four isogenies, with composites multiplication by two. The
[ordinary Kummer diagram](JACOBIAN_LOCAL_CONDITIONS_AND_CT.md) gives
\[
0\longrightarrow S_D\longrightarrow\mathrm{Sel}_2(J_s)
\longrightarrow\operatorname{rad}(\Delta|_{S_C})\longrightarrow0,
\]
where `Delta` is the sum of the two correctly labelled elliptic CT forms.
The Selmer pairing kernel and additivity are
[Morgan–Smith, Theorem 1.3 and Proposition 4.4](https://arxiv.org/pdf/2103.08530).

For `J+`, the twist's full Selmer group is rational and its CT form is
zero. In basis `(beta0,beta1,u0)`, the difference is `H + 0^1`, so its
radical is `<u0>`. This gives `dim Sel2(J+)=3+1=4`.
For `J-`, both forms on `U` are the same nondegenerate alternating
plane after transport: every invertible map of a two-dimensional
F2-space preserves it. Their sum is zero. Thus all of `U` survives,
and `dim Sel2(J-)=4+2=6`.

There is a useful general refinement. Let `G0,G1` denote the **full
rational elliptic Kummer images** in the common cohomology space.
Assume the common 2-torsion module has no rational invariants, as here.
The rational-point Kummer diagram gives
\[
0\longrightarrow G_0+G_1\longrightarrow J(\mathbb Q)/2J(\mathbb Q)
\longrightarrow G_0\cap G_1\longrightarrow0.
\]
Indeed, a rational Jacobian point projects to a pair with equal Kummer
class. Conversely, equality makes the `Psi`-Kummer obstruction zero,
so the pair lifts rationally. The kernel consists of `Phi` images,
whose 2-Kummer classes are the sums of the two elliptic classes.
The absence of rational kernel invariants makes the left map injective.

Quotienting this row from the Selmer row proves the exact F2-sequence
\[
\boxed{
0\longrightarrow \frac{S_D}{G_0+G_1}
\longrightarrow\Sha(J)[2]
\longrightarrow\frac{\operatorname{rad}(\Delta|_{S_C})}{G_0\cap G_1}
\longrightarrow0.}
\]
The left quotient measures excess in the summed local conditions; the
right measures common Selmer lifts without common rational representatives.
This is an exact-sequence consequence, not a claim of a canonical direct
sum or a new theorem that a zero CT radical produces points.

Here the two quotient dimensions are `(0,0)` for `J+` and **`(2,2)`
for `J-`**. In the latter case the sequence has the form
`0 -> U -> Sha(J-)[2] -> U -> 0`. Its two-dimensional quotient comes
from Selmer lifts of the very same strict block. No rational Jacobian
point projects to a nonzero class of that block.

The Jacobian ranks equal the sums of the elliptic ranks under the
isogenies: four and two. Since rational Jacobian 2-torsion is zero,
the Selmer dimensions independently give Sha[2] dimensions zero and four.
No assertion about higher 2-primary Sha is needed.

## Rational lifting, halving, and killing elliptic Sha

For `J+`, `G0=<u0>` and the pulled-back second rational image fills `S0`.
Their intersection is one-dimensional. It has an explicit rational
Jacobian representative
\[
j=[\infty_+-(0,1)],\qquad
\Psi(j)=(-P,R),\quad P=(-1,1),\ R=(0,1).
\]
The point `infinity+` is the rational branch with `Z/X^3 -> +1`.
The maps above verify its projections, and `u0*tau(theta)=1` verifies
the matching Kummer class. The relation
`2j=Phi(-P,R)` gives index two for the product's rational image in
`J+(Q)`. This is a halving contribution, not an additional free rank.

For `J-`, the two rational images are `<u0>` and `<u1>`, with zero
intersection. Therefore `Phi(E0(Q) x E0(Q))=J-(Q)`. Isomorphic elliptic
quotients do not make this a split labelled extension: the gluing acts
on 2-torsion by `tau`, rather than the identity.

There is also an instructive contrast of maps. The original elliptic
Sha[2] is killed in `J+`, whose Sha[2] is zero, but its strict classes
do not lift through the right-hand Jacobian Selmer projection. In `J-`
they do lift through that projection, and every lift remains non-rational.
Killing a torsor in an auxiliary abelian variety and lifting a common
Selmer class are different assertions.

## Ranked implications for rank-jump analysis

1. **Solubility:** a common arithmetic object carries the whole obstruction
   block; here the same two-dimensional strict space accounts for both
   terms of a four-dimensional Jacobian Sha[2] space. This is structural
   block behaviour, not independent point-search failure.
2. **Disproved sufficiency:** even ordinary Jacobian Selmer lifting,
   after all local conditions, need not yield a rational direction.
   The earlier explicit norm witnesses can be corrected to some Selmer
   lifts in the negative case, by the pairing-kernel theorem. This does
   not assert that the stored norm witnesses already give the required
   local cocycle representatives without correction.
3. **Incidence:** the exact sequence distinguishes the summed Selmer
   excess and the common CT-radical excess. It applies to production
   comparisons only after using the correct labels and full common
   Selmer space; a radical computed on a retained subspace is insufficient.
4. **Weak rank explanation:** self-gluing, a smaller Galois closure,
   and a larger Selmer lift space do not force higher rational rank.
   The positive sign's extra rational lift is an index-two saturation
   event with no additional free rank.
5. **Missing production implication:** an equation-defined incidence
   block still needs common rational representatives, followed by
   independence modulo the original generic subgroup. The groups `G0,G1`
   in the formula are not point-blind selectors. Agent 1 can use the
   exact sequence as an obstruction/accounting framework; it supplies
   no validated candidate score or chart policy.

The [calculation](../../artifacts/generated-results/elliptic-curves/rank_jump_small_jacobian_selmer_v1.json)
and [independent replay](../../artifacts/generated-results/elliptic-curves/rank_jump_small_jacobian_selmer_verification_v1.json)
retain labels, local dimensions, torsion checks and both Sha quotients.
No active-search output or mathematical-status entry was changed.

The sealed producer's direct `--check` compares integer dictionary keys
with JSON string keys and fails on that representation difference. The
replay below recomputes the same output and compares canonical JSON,
preserving the original source and certificate. The independent verifier
checks the mathematical inputs separately.

```sh
python3 elliptic-curves/rank-jump/replay_small_jacobian_selmer.py
sage -python elliptic-curves/rank-jump/verify_small_jacobian_selmer.py --check
python3 -m unittest discover -s elliptic-curves/rank-jump -p test_small_jacobian_selmer.py
```
