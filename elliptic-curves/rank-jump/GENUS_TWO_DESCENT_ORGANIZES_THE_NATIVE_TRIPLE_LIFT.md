# The native triple lifts through a genus-two isogeny condition

The successful `01333,0b2d0,19e45` triple has a smaller **unramified
descent base** than its genus-five simultaneous-solubility carrier:

\[
\boxed{H:\ y^2=h(t)f(t)g(t),\qquad
C_3\longrightarrow H,\quad(t,a,u,v)\longmapsto(t,auv).}
\]

This map is an étale degree-four cover. Its exact lifting condition is
membership in the image of a specific degree-four isogeny of abelian
surfaces. Consequently the unresolved global task can be posed on a
**genus-two curve with a specified descent class**, rather than by searching
the rank-three elliptic pair carrier for further squares.

The computations also give

\[
\operatorname{rank}J(C_3)(\mathbf Q)
=8+\operatorname{rank}J(H)(\mathbf Q)\geq8.
\]

Thus the ordinary rank-below-genus Chabauty criterion cannot apply directly
to C3. The rank of J(H), and the complete set of its curve's rational
points in the required descent class, remain **UNKNOWN**.

This explains the arithmetic object controlling a retained **two-direction
subblock** of the +7 fibre `08234-003`. It does not explain all seven
quotient directions. The generic subgroup has rank 17; the retained
independent subgroup has rank 24. The triple's specialized quotient rank
is exactly two, as proved by the
[existing carrier and intersection certificate](MINIMAL_CARRIER_AND_RATIONAL_SPLITTING_OF_A_TWO_DIRECTION_BLOCK.md).

## Equations and the precise condition on t

Use the fixed primitive native forms

\[
\begin{aligned}
h(t)&=3255283501715844+1254304425186516t+125950947365881t^2,\\
f(t)&=16496921457+11037654810t+1654807609t^2,\\
g(t)&=4897771825+2856794060t+580965316t^2.
\end{aligned}
\]

Their roots are simple and pairwise disjoint, and all three quadratics
are irreducible over Q. Hence no finite rational t is a branch value.
The full marked carrier is

\[
C_3:\ a^2=h(t),\quad u^2=f(t),\quad v^2=g(t).
\]

For finite rational t the necessary and sufficient condition is

\[
\boxed{\exists y\in\mathbf Q:\quad y^2=h(t)f(t)g(t)
\quad\text{and}\quad([h(t)],[f(t)])=(1,1).} \tag{1}
\]

The classes belong to Q*/Q*². When (1) holds, choose rational a and u
and set v=y/(au). This constructs all three native points through the
already certified native point maps. Independence after specialization
is a separate check: it holds with quotient rank two at the retained
parameter, and is not inferred for arbitrary solutions of (1).

The extra content of the genus-two formulation is that the last two
squareclasses are the restriction of an **isogeny descent homomorphism
on J(H)**. They are not unrelated tests fitted to the exceptional points.

For an explicit normalization take the retained point

\[
P_0=\left(-\frac{288}{65},
\frac{147350257714085825322}{21125}\right)\in H(\mathbf Q).
\]

Write J=J(H). Let K be the rational isotropic subgroup of J[2] defined
by the three quadratic factors, let A=J/K, and let
\(\phi:J\to A\), \(\psi:A\to J\) be the quotient and dual
isogenies, normalized so that \(\psi\phi=[2]\). Then (1) is equivalently

\[
\boxed{\exists y\in\mathbf Q:\quad y^2=h(t)f(t)g(t),\qquad
[(t,y)-P_0]\in\psi\bigl(A(\mathbf Q)\bigr).} \tag{2}
\]

The factorization and the cover exist before any exceptional point is
supplied. P0 only identifies their particular soluble twist with a
pointed pullback. This use of the known point is explicitly retrospective;
it is not an admissible prospective selector.

## Why this is a degree-four isogeny pullback

Let Wh, Wf, Wg denote the respective degree-two Weierstrass divisors on
H, and put I=∞+ + ∞−, a rational divisor even if its points are not
individually rational. On H,

\[
\operatorname{div}(h)=2(W_h-I),\quad
\operatorname{div}(f)=2(W_f-I),\quad
\operatorname{div}(g)=2(W_g-I).
\]

The three divisor classes Dh, Df, Dg are the nonzero elements of a
two-dimensional rational subgroup K of J[2]. Their sum is zero because
\(W_h+W_f+W_g-3I=\operatorname{div}(y)\). In the even-subset
description of hyperelliptic 2-torsion, they are the three disjoint pairs
in the six branch points. They are distinct nonzero classes modulo
complement and have trivial pairwise Weil pairing. Thus K is maximal
isotropic and A carries the quotient principal polarization.

The character group of ker ψ is K. Pulling ψ back along
\(i:H\to J,\ P\mapsto[P-P_0]\) therefore gives the two quadratic
extensions defined by

\[
\sqrt{h(t)/h(t_0)},\qquad\sqrt{f(t)/f(t_0)}.
\]

Both normalizing constants are rational squares. Since y supplies
\(\sqrt g=y/(\sqrt h\sqrt f)\), this is exactly Q(C3). The fibre
over P0 is split, so no constant torsor twist is left undetermined.
Equivalently,

\[
C_3\simeq H\times_{J,i,\psi}A.
\]

The connecting homomorphism for ψ, restricted to i(H), evaluates the
functions h and f on P−P0. Its kernel is ψ(A(Q)), proving (2).
The general isogeny descent construction for a genus-two sextic split
into three quadratics is described by
[Arnth-Jensen and Flynn, §2](https://people.maths.ox.ac.uk/flynn/genus2/af/artlong.pdf).
The broader unramified-cover framework and its local/global distinction
are developed by
[Bruin and Stoll](https://mathe2.uni-bayreuth.de/stoll/papers/twocoverdescent.pdf).

This is the nondegenerate Richelot situation: the coefficient determinant
of h,f,g, using columns 1,t,t², is

```text
2756819365126924444728271453799958.
```

The three derivative brackets are quadratic, squarefree and mutually
coprime. Their coefficients are retained in the certificate. No genus-two
Selmer computation or explicit Richelot point correspondence has been
claimed merely from these algebraic checks.

## The two notions of minimality

The two marked directions B=`0b2d0` and D=`19e45` have the minimal
function-field carrier

\[
C:\ u^2=f(t),\quad v^2=g(t),\qquad g(C)=1.
\]

The three marked native point maps require all three independent
quadratic characters, so their carrier has degree eight and genus five.
Their specialized relation does not remove a generic character. These
minimality assertions concern these specified native maps, not every
possible construction of the observed quotient.

H is minimal in a different, useful sense: it is the unique quotient by
a maximal freely acting subgroup of the native deck group. In
G=(Z/2)³, the three inertia elements are the individual sign changes.
The even-parity subgroup

\[
G_0=\{000,011,101,110\}
\]

avoids all three and acts freely, including at infinity. It is the unique
order-four subgroup with this property. Its invariant root is auv,
giving H. Riemann–Hurwitz reads \(8=4(2g(H)-2)\), so g(H)=2.
No nonconstant unramified map from C3 to a genus-zero or genus-one curve
can exist, again by Riemann–Hurwitz. This makes H a minimal-genus
unramified descent base, not a substitute for the full carrier's rational
lifting condition.

## Why the elliptic base is a harder global formulation

Forgetting a instead gives C3→C, a ramified double cover. Over each of
the two roots of h, the pair carrier has four geometric points; h has
simple zeros there. Its poles at the four infinity points have even order
two. Thus the relative ramification divisor has eight points.

The exact algebra

\[
\mathbf Q[t,u,v]/(h(t),u^2-f(t),v^2-g(t))
\]

is one degree-eight field. To certify this without a guessed Galois
group, the quadratic field Q[t]/h has nonsquare norms for each of
f, g and fg. Thus their three nontrivial squareclasses remain nontrivial
there. A primitive element t+u+v is also checked by independent resultant
elimination against its multiplication matrix. The eight ramification
points do not split into smaller rational components.

C has exact Jacobian rank three and hence infinitely many rational
points. C3 has genus five and a known point, so its rational points,
and their images in C(Q), are finite and nonempty by Faltings' theorem.
See [Stoll's survey](https://www.mathe2.uni-bayreuth.de/stoll/papers/JA2009-paper.pdf)
for the rational-point finiteness and ordinary Chabauty framework.
It follows that the liftable subset of C(Q) cannot be a union of cosets
of any finite-index subgroup: every nonempty such coset is infinite.
In particular, a criterion solely on C(Q)/nC(Q) cannot characterize it.

Even a fixed finite collection of local lifting tests on C cannot suffice.
Near the known unramified lift, h remains a nonzero local square. For
any finite set of places, compactness of the product of the local elliptic
groups gives infinitely many multiples of a nontorsion rational point
returning to a sufficiently small neighborhood of zero. Translating by
the known base point gives infinitely many rational points passing every
one of those local tests, but only finitely many lift globally. This
statement concerns a fixed finite set of places while t varies. For a
specified rational t, squareness at **all** places is sufficient for a
rational number to be a square.

There is no contradiction with (2): the finite-index condition there
holds in the genus-two Jacobian and must be intersected with the embedded
curve H. That curve-point requirement is the essential global constraint.

## Jacobian factors and exact bounded descents

The four nonzero differential character spaces of C3 give

\[
J(C_3)\sim_{\mathbf Q}
E_{hf}\times E_{hg}\times E_{fg}\times J(H),
\]

where Eq is the Jacobian of y²=q(t). The pair-product characters each
have one differential, the triple-product character has two, and the
three singleton quotients have genus zero. Pullback on differentials
therefore proves that the resulting homomorphism of equal-dimensional
abelian varieties is an isogeny. The Prym of C3→C has dimension four
and is isogenous to Ehf × Ehg × J(H).

Initial effort-zero descents on the product models did not determine
their ranks. Using the retained lift to parameterize one conic gave the
two-isogenous pair models and sharpened all bounds:

| Character | Initial product rank interval | Exact rank via pair model | Pair-model Sha[2] dimension | Product-model Sha[2] dimension |
|---|---:|---:|---:|---:|
| hf | [1,3] | 3 | 0 | 14 |
| hg | [0,4] | 2 | 2 | 16 |
| fg | [1,3] | 3 | 0 | 10 |

Each degree-two isogeny is checked exactly. The computed pair models,
isogeny kernels, raw PARI outputs and returned points are retained.
The rank bounds and dimensions use
[PARI's documented ellrank semantics](https://pari.math.u-bordeaux.fr/dochtml/html-stable/Elliptic_curves.html#ellrank).
Replaying the same PARI implementation is a reproducibility check, not
an independent descent algorithm. The product hg model also has
dim(2Sha[4])=2; its initial unresolved interval was consistent with this.

These large auxiliary Sha groups coexist with the known rational lift
on C3. They must not be labelled an obstruction to that already soluble
carrier. They demonstrate why isogenous-model Selmer dimensions alone
would misrepresent this problem.

The elliptic rank sum is eight. Ordinary Chabauty's rank<genus hypothesis
fails on C3 regardless of the remaining genus-two rank. This does not
rule out Chabauty on H, covering methods or quadratic Chabauty.

Independent finite-field checks include all branch and infinity points:

| Prime | #H(Fp) | H points that lift | #C3(Fp) |
|---|---:|---:|---:|
| 131 | 142 | 29 | 116 |
| 137 | 124 | 30 | 120 |

Each soluble fibre has four points. These checks calibrate the exact
lift gate; they are not frequencies for rational t or rank predictions.

## Mechanism, missing implication and next falsifiable computation

The strongest current mechanism is **rational splitting of a specific
arithmetic cover/intersection**. For the retained generic relation, the
previous exact degree-twelve intersection has factorization 1+11, and
its unique rational parameter satisfies **65t+288=0**. That certifies one
soluble two-direction subblock. It does not enumerate H(Q) or C3(Q).

The global carrier problem is now precisely

\[
i(H(\mathbf Q))\cap\psi(A(\mathbf Q)).
\]

Knowing that J(H) or its isogeny Selmer group has many classes is
insufficient: those classes must meet the embedded curve and the
specified isogeny image. The missing implication is a complete rational
point calculation on H in that class, followed by original-curve
independence checks for each parameter obtained.

A small falsifiable next experiment is a **single fixed genus-two
Jacobian descent** on this H and, if needed, its Richelot quotient A.
Predeclare one worker, at most 120 seconds per model, checkpoints before
each descent, and no new original-parameter point search. The test is
whether an unconditional rank upper bound at most one makes ordinary
Chabauty available on H. An exact lower bound at least two refutes that
route; a wider interval or timeout remains UNKNOWN. Actual Chabauty
would then require a separately certified subgroup and its needed
saturation. This experiment is designed here, not executed.

For eventual use by Agent 1:

1. **Solubility:** preserve the specific cover class and rational splitting
   data. A rational point on H in the trivial two-squareclass fibre supplies
   simultaneous native points.
2. **Incidence:** certify how many of those native maps remain independent
   modulo the original generic subgroup after specialization. A third
   soluble cover can add zero quotient rank.
3. **Weak explanation:** auxiliary rank, raw Selmer size, small trace norms
   and a fixed finite local panel do not supply the missing global curve
   points. The matched degree-twelve controls already separate identical
   lattice capacity from rational splitting.
4. **Visibility:** no new visibility feature or search policy follows from
   this calculation. No active search file or candidate population changed.

## Reproduction

Frozen inputs and outputs are under
`artifacts/generated-results/elliptic-curves/rank_jump_native_genus_five_lift*`,
`rank_jump_native_pair_factor_descent*`, and
`rank_jump_native_genus_two_lift_gate_v1.json`.
The two descent protocols bound five new auxiliary elliptic descents to
one worker and sixty seconds each. All completed. The previous fg pair
descent was reused. The algebra-only genus-two gate enumerates the 256
subsets of the eight-element native deck group and checks the two fixed
finite fields; it performs no descent or parameter search.

```sh
sage -python elliptic-curves/rank-jump/verify_native_genus_five_lift.py check
sage -python elliptic-curves/rank-jump/native_genus_two_lift_gate.py check
```
