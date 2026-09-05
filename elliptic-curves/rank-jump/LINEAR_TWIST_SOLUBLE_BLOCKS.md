# An explicit soluble block and a splitting event

There is an exact simultaneous-solubility mechanism in the fixed-cubic
pencil. For independent anchor points `P_i=(a_i,b_i)`, the conditions

    1-a_i*u = z_i^2 in Q*,   i=1,...,k

produce explicit rational points on `E_u` whose Kummer classes have
rank at least `k-1`. This is a proved rank lower bound, not a visibility
score. All conditions refer to one genus-two correspondence. At `u=0`
that correspondence splits into two copies of the anchor elliptic curve,
and all the original directions survive together.

This explains the engineered rank-at-least-20 anchor of this
**generic-rank-zero pencil**. It does not yet explain the large jumps
in the published R17 or A1/MW16 families. The anchor points are public
exceptional data used retrospectively; they must not enter Agent 1's
prospective selector.

## The point formula

Use the family and assumptions of the
[rank-jump reassessment](../notes/RANK_JUMP_REASSESSMENT_2026-09-05.md):

\[
 f(a)=a^3+Aa+B,\quad D(u)=1+Au^2+Bu^3,
\]
\[
 E_u:\ y^2=x^3+2Au x^2+(A+3Bu+A^2u^2)x+B+ABu^2-B^2u^3.
\]

Assume the anchor cubic is irreducible, `B≠0`, and its discriminant is
nonzero. Let `P=(a,b)∈E_0(Q)`, with `b≠0`. Set `g=1-au`. The
quadratic twist by `g` has the section

\[
 X=a+(Aa+B)u^2,\qquad Y=bD(u).
\]

Indeed the universal polynomial identity is

\[
 X^3+2Au\,gX^2+(A+3Bu+A^2u^2)g^2X
 +(B+ABu^2-B^2u^3)g^3=f(a)D(u)^2.
\]

Thus whenever `g=z²≠0` and `D(u)≠0`,

\[
 \boxed{\quad
 P_u=\left(\frac{a+(Aa+B)u^2}{z^2},
                 \frac{bD(u)}{z^3}\right)\in E_u(\mathbf Q).
 \quad}
\]

At `u=0,z=1` this is exactly `P`. No point search or genus-one
solubility assumption occurs after the square test.

## The generic incidence locus for linear twists

The formula gives more than a sufficient construction. For any rational
`a` with `f(a)≠0`,

\[
 \operatorname{rank}E_u^{\,1-au}/\overline{\mathbf Q}(u)=1,
\qquad
 \operatorname{rank}E_u^{\,1-au}/\mathbf Q(u)
 =\begin{cases}1&f(a)\in\mathbf Q^{\times2},\\0&\text{otherwise}.\end{cases}
\]

For `a=0` this is the reassessment's original theorem. For `a≠0`,
twisting removes the `I_0^*` fibre at infinity and puts an `I_0^*` at
`u=1/a`. The three `I_2` fibres remain; `f(a)≠0` keeps `1/a` away
from them. The minimal surface is still rational, so Shioda–Tate gives
geometric rank `10-2-4-3=1`.

For completeness, take `v=u/(1-au)` and the integral coordinates
`x_v=(1+av)^2X`, `y_v=(1+av)^3Y` on the twisted model. Its section is

\[
 x_v=a+2a^2v+f(a)v^2,\quad
 y_v=\sqrt{f(a)}\,[1+3av+(3a^2+A)v^2+f(a)v^3].
\]

It is disjoint from the zero section, meets the nonidentity component
at each of the three `I_2` fibres, and meets the identity component
at the `I_0^*` at infinity. Its height is `2-3/2=1/2`. It therefore
spans the geometric rank-one space, on which Galois acts through
`sqrt(f(a))`. This proves the arithmetic rank formula.
The surface and height tools are the same as those used in the
reassessment; see [Schütt–Shioda, Sections 8 and 11](https://arxiv.org/abs/0907.0298).

Consequently the anchor curve `b²=f(a)` is the exact arithmetic
incidence locus for these rank-one linear twists. This statement concerns
function-field ranks, not upper bounds on specialized curves.

## A common Kummer shift gives the block bound

Let `K=Q(θ)` with `f(θ)=0`, `α=θ+uθ²`, and

\[
 \beta_i=a_i-\theta,\quad
 \gamma=1-u\theta,\quad
 \kappa=1+u\theta+u^2(A+\theta^2),\quad
 \eta=D(u)\gamma.
\]

Direct multiplication gives `γκ=D(u)` and, for every transported point,

\[
 x(P_{i,u})-\alpha
   =\frac{\beta_i\kappa}{z_i^2}
   =\frac{\beta_i\eta}{(z_i\gamma)^2}.
\]

Therefore its Kummer class is `[β_i]+[η]`. The shift is the same for
every point; `N(η)=D(u)^4` is a square. This transporter does **not**
preserve each inherited class separately. That distinction matters when
comparing it with the previously obstructed inherited CT classes.

Suppose the `k` classes `[β_i]` are independent. Every even-cardinality
sum of transported classes equals the corresponding sum of anchor
classes. Thus the `k-1` differences

    P_(i,u)-P_(1,u),   i=2,...,k

have independent Kummer classes `[β_iβ_1]`. Since `D(u)≠0` preserves
the irreducible cubic algebra, `E_u(Q)[2]=0`. Hence

\[
 \boxed{\operatorname{rank}E_u(\mathbf Q)\ge k-1.}
\]

More precisely, the transported classes have rank `k` unless `[η]`
belongs to their anchor span with odd coordinate sum; in that case
their rank is `k-1`. At `u=0` the shift is zero and all `k` classes
remain independent. This proof does not rely on generic independence
surviving specialization, the failure measured in the
[previous cover comparison](BRANCH_BLOCKS_AND_SPECIALIZATION.md).

## One genus-two curve controls simultaneous solubility

For fixed `u≠0` with `D(u)≠0`, consider

\[
 C_u:\qquad b^2=f(a),\qquad z^2=1-au.
\]

Over the `a`-line the two double covers share the branch point at
infinity. Their union has five branch points, so the connected
normalization has genus two. The two elliptic quotients are `E_0` and

\[
 H_u:\quad r^2=(1-au)f(a),\qquad r=bz.
\]

The second quotient is birational over Q to `E_u`:

\[
 x=\frac{a+(Aa+B)u^2}{1-au},\qquad
 y=\frac{rD(u)}{(1-au)^2}.
\]

Away from their poles the inverse is

\[
 a=\frac{x-Bu^2}{1+Au^2+xu},\qquad
 r=\frac{yD(u)}{(1+Au^2+xu)^2}.
\]

The identities extend to the smooth projective curves. Thus the usual
bielliptic construction gives a degree-four isogeny

\[
 \operatorname{Jac}(C_u)\longrightarrow E_0\times E_u.
\]

One can see the splitting on differentials: the two degree-two quotient
maps give complementary eigenspaces for the involutions. Their norm and
pullback maps compose to multiplication by two on each elliptic factor.
This is the classical `(2,2)` decomposition, not a new general theorem;
compare [Wetherell, Chapter 2, Section 4](https://swc-math.github.io/aws/1999/99WetherellThesis.pdf).

At `u=0` the equation `z²=1` splits. Affinely there are two copies of
`E_0`, each mapping isomorphically to `E_u=E_0`. In the finite double
cover compactification over `E_0` they meet at the point at infinity.
For a local parameter `t` there, writing `w=tz` gives
`w²=t²-u(t²a)`, where `t²a` is a unit with value one.
The central fibre has the node `(w-t)(w+t)=0`.

This is an actual component-splitting event, visible from the equations
before any exceptional point is supplied. It releases an entire copy of
the anchor Mordell–Weil group at once. With the pinned anchor it gives
twenty independent rational directions against generic rank zero.
The construction deliberately starts from that anchor; it is not evidence
that unrelated high-gain fibres share this degeneration.

## Frozen retrospective experiment

The [protocol](LINEAR_TWIST_BLOCK_PROTOCOL.json) uses exactly twenty
already pinned anchor points and `u=-3,-2,-1,0,1,2,3`.
The [certificate](../../artifacts/generated-results/elliptic-curves/rank_jump_linear_twist_blocks_v1.json)
retains all 140 square tests, the twenty polynomial identities, all
transported points and their exact Kummer square corrections.

| u | -3 | -2 | -1 | 0 | 1 | 2 | 3 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Soluble transports from the fixed basis | 0 | 0 | 0 | 20 | 0 | 0 | 0 |

The six misses exclude only this fixed basis transporter. Combinations
of anchor points, other points on `C_u`, and rational points of `E_u`
that do not lift to `C_u` remain outside the test. No upper rank follows.

    python3 elliptic-curves/rank-jump/linear_twist_blocks.py check
    sage -python elliptic-curves/rank-jump/linear_twist_blocks.py verify

The second command independently checks the universal point formula and
inverse quartic identity symbolically. The rank and degeneration statements
are proved above; they are not inferred from the 140 tests.

## A representative-independent obstruction to large transported blocks

Let `W` be the twenty-dimensional inherited anchor Kummer space and
`W_u=W∩Sel_2(E_u)` the previously computed complete local intersection.
Let `R_u` be the radical of the certified CT pairing on `W_u`.

If `k` independent classes in `W` have rational point representatives
on the anchor that satisfy the transport square conditions, their
`k-1` even differences are rational Kummer classes on `E_u`. They
therefore belong to `W_u` and pair trivially with every element of
`W_u` under CT. Thus

\[
 k-1\leq\dim R_u.
\]

This does not require the affine shift `η` itself to belong to W or
to the Selmer group. It bounds a whole independent anchor-class block,
not just the twenty chosen coordinate representatives.

| u | dim W_u | CT rank on W_u | dim R_u | Maximum independent anchor classes that can satisfy transport together |
|---|---:|---:|---:|---:|
| -3 | 17 | 16 | 1 | 2 |
| -2 | 13 | 12 | 1 | 2 |
| -1 | 18 | 16 | 2 | 3 |
| 0 | 20 | 0 | 20 | 20, limited by dim W |
| 1 | 13 | 12 | 1 | 2 |
| 2 | 13 | 12 | 1 | 2 |
| 3 | 15 | 14 | 1 | 2 |

These caps survive changes of representatives and replacement by any
other independent basis inside W. They do not bound the total number
of transported points, the full curve rank, or the rank of an
unsaturated subgroup whose Kummer images are dependent. The radicals
are necessary receptacles for rational classes, not proofs of solubility.

## What this changes, and what remains missing

1. **Incidence, proved:** the linear-twist rank-one locus is exactly the
   anchor elliptic curve `f(a)=square`.
2. **Solubility, proved:** simultaneous squares `1-a_i u` give rational
   points, and the common Kummer shift certifies a block of dimension
   at least `k-1`. At the reducible member `u=0` the full anchor block
   survives. At the six retained nonzero members, CT forbids a large
   independent inherited block through this transporter, regardless
   of the chosen point representatives.
3. **Weak explanation:** a shared cubic field, shared quadric or shared
   ruling field by itself. The companion
   [quadric analysis](SHARED_QUADRICS_AND_RULINGS.md) shows the common
   geometry but also forty distinct base changes.
4. **Missing:** a comparable point-independent correspondence or
   degeneration in the R17/MW17 and A1/MW16 high/low pairs. Their
   quotient dimensions must be measured relative to their own generic
   subgroups, not this generic-rank-zero pencil.

Agent 1 can eventually use the **form of the certificate**: a fixed
auxiliary construction, explicit rational-solubility conditions, and a
uniform Kummer correction giving a quotient-independence bound. The
public anchor-point square tests here are not an authorized selector.

The next useful bounded experiment is to determine whether the affine
coset `η_u+W` meets `Sel_2(E_u)` at each of the six retained nonzero
parameters. This is a new incidence question: the transported classes
are shifted, whereas the old descent computed only the inherited
linear space W_u. Solve the affine local conditions without a point
search. If a nonempty coset extends W_u, the next CT calculation needs
only one additional representative against the existing basis.

This would locate the precise missing class and its solubility
obstruction. It would not authorize treating local compatibility,
or a radical, as a rational point.
