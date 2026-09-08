# Two determinant-1092 covers: rank19 and the genus-one base obstruction

## Verified application

The two already constructed rational bisections, orbits **8044 and 127449**,
each give one explicitly verified extra section on their own quadratic base
change. On their common degree-four cover the two directions are independent,
giving a subgroup of **rank at least19**. The common base has **genus one**,
so simultaneous splitting cannot be parametrized by a rational parameter line.
The existence of a rational point on that common base remains `UNKNOWN`.

This follows the independently checked
[rank18 construction and initial-unlock obstruction](DET1092_RANK18_BASE_CHANGE_AND_INITIAL_UNLOCK_2026-09-08.md).
The present calculation uses precisely those two covers; it selects no new
orbits, runs no rational point search, and does not change the eight-fibre V3
pilot or interpret its running checkpoints as final results.

The complete [certificate](../../artifacts/generated-results/elliptic-curves/det1092_two_cover_rank19_genus_gate_v1.json)
and [independent replay](../../artifacts/generated-results/elliptic-curves/det1092_two_cover_rank19_genus_gate_replay_v1.json)
contain the exact equations, second conic parametrization, lifted section,
trace identity, and positive-definite 19-section height matrix.

## Explicit common base

Let `q1(t)` be the quadratic for orbit8044 in the rank18 certificate. The
second bisection discriminant is exactly `h2(t)^2*q2(t)`, with

```text
q2(t) = alpha2*t^2 + beta2*t + gamma2
alpha2 = 442672757925512203739677418804643828333652510289675767653949339776753457
beta2  = 97444477872846123883554680615682804397661851854905069920795047645996400
gamma2 = 10437655460487956834254111412184389987866303895084604973565973964840000.
```

Both quadratics are irreducible and squarefree over Q, and their gcd is one.
Both are coprime to the parent discriminant. Their geometric branch pairs
are disjoint. In particular `q1`, `q2` and `q1*q2` are nonsquares even after
extending the constant field to the algebraic closure.

The smooth common base `B` has function field

\[
 \mathbf Q(B)=\mathbf Q(t,s_1,s_2),\qquad
 s_1^2=q_1(t),\quad s_2^2=q_2(t).
\]

It is geometrically connected and degree four over `P1_t`, with deck group
`(Z/2)^2`. If `t=T(u)=N(u)/D(u)` is the already displayed parametrization of
the first conic, an equivalent explicit genus-one model is

\[
 v^2=D(u)^2q_2(T(u)),\qquad t=T(u),\quad s_2=v/D(u).
\]

The equations and maps are defined over Q. Each individual conic has a
verified rational point and parametrization; this does not imply that the
common genus-one curve has a rational point. Its rational points would give
the simultaneous square conditions `q1(t),q2(t) in Q^2`, outside the recorded
chart exclusions. No such point is asserted here.

## Independence is proved using sections, not only characters

**Verified application of established height theory.** On each quadratic
cover, the residual-quadratic and RR-line formulas construct a point `Qi`.
The checker verifies the elliptic equation, the exact conjugate sum
`Qi+sigma_i(Qi)=Wi`, and polynomial coordinates of degrees at most8,12 on
the `chi=4` model. Thus `Qi.O=0`, its height is8, and

\[
 V_i=2Q_i-W_i,\qquad \langle V_i,V_i\rangle=12.
\]

On the common cover these heights double to24. The two nonzero sections
transform under different deck characters, `(-1,+1)` and `(+1,-1)`.
Applying the two character projections to a dependence proves their
independence from each other and from the invariant rank17 subgroup.
Their heights establish that neither character contributes only torsion.

For the fixed trace words, the exact parent pairing is `w1^t G w2=0`.
The 19-section Gram has blocks

\[
 H=\begin{pmatrix}
 4G & 2Gw_1 & 2Gw_2\\
 2w_1^tG &16&0\\
 2w_2^tG &0&16
 \end{pmatrix}.
\]

Its Schur complement is `diag(6,6)`, so it is positive definite and has
determinant `36*4^17*1092`. This proves rank at least19 over **Q(B)**. It
does not give a rationally parametrized rank19 family or certify rank19
at a rational specialization.

## Exact obstruction to a rational base

**Established literature.** Use the characteristic-zero
[Riemann--Hurwitz formula](https://stacks.math.columbia.edu/tag/0C1B).

**Verified application.** The four geometric branch points each have two
points above them with ramification index2. Hence

\[
 2g(B)-2=4(-2)+4\cdot2=0,
 \qquad g(B)=1.
\]

Any nonconstant rational parametrization by `P1` would extend to a morphism
of smooth projective curves. Riemann--Hurwitz would give `-2=0+R` with
`R>=0`, a contradiction. This obstruction persists after any constant-field
extension; it is not removed by finding a point over a number field.

This does not obstruct rational points or even infinitely many rational
points on `B`. Those require further arithmetic of this specific genus-one
curve. The retained primitive simultaneous-square checks at 11 fixed prime
powers through `2^10`, `3^6`, `5^4` and `31^2` all have survivors. They
exclude an obstruction at those finite moduli only; they prove neither
everywhere local solubility nor global solubility.

## New deduction: an exact construction gate

For `k` independent quadratic characters whose `2k` simple geometric branch
points are pairwise distinct, the common degree-`2^k` cover satisfies

\[
 2g-2=-2\,2^k+2k\,2^{k-1},\qquad
 g=1+2^{k-1}(k-2).
\]

Thus the base genera for one through four such covers are `0,1,5,17`.
Independently proved extra sections can coexist with rapidly increasing
base genus. This supplies an exact gate before spending specialization
budget: check branch overlap and the genus of the simultaneous splitting
base, and require a justified rational-point source there.

To retain a rational parameter line, the alternatives include specially
shared branch loci or multiple independently proved sections on one
quadratic character. This is a necessary geometric direction, not an
assertion that either option has been constructed on this parent. A branch
collision or repeated squareclass still needs actual section identities and
an independence proof.

The gate helps interpret completed negative exposure: a fixed generic-start
chart search does not force these square conditions. Conversely a construction
that forces several conditions must still solve the arithmetic of their
common base. No inference about the full ranks of the sampled pilot fibres
follows from this calculation.

## Bounded replay

The additional constructor has a 40-second execution cap in the recorded
run and completed within it. The independent replay uses no conic solver
and has a 30-second cap. No detached job is left by this experiment.

```sh
sage -python research/elliptic-curves/cas/verify_det1092_two_cover_gate.sage
```

The constructor `construct_det1092_two_cover_gate.sage` refuses to overwrite
its certificate. The previous rank18 theorem and both original bisection
certificates remain intact. The equations and computations are labelled
**verified application**, the general genus calculation **new deduction**,
and the height/Riemann--Hurwitz tools **established literature**. Rational
solubility of the common base remains unresolved rather than conjecturally
promoted.
