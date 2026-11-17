# Parametric inverse Galois theory for quintic Keller fibers

## Status

This is a research audit and terminology correction.  Its main Keller
theorem is the absolute quintic corollary proved in
[`verified/UNIVERSAL_RELATIVE_KELLER_MAP.md`](verified/UNIVERSAL_RELATIVE_KELLER_MAP.md):
one explicit determinant-one map of `A^5_Q` realizes every quintic number
field as a complete fiber.  Thus the unrestricted fixed-map minimum is
already one, not two.

For the particular split-seed map of `A^3`, universality remains open.  The
existing exact work proves generic families and infinitely many fields, not
field-surjectivity.  For one-parameter target curves, published local
specialization theorems rule out every finite collection for `A_5` and
`S_5`.  No corresponding conclusion is asserted here for `C_5`, `D_5`, or
`F_{20}`.  The separate
[`KELLER_BECKMANN_BLACK_SPECIALIZATION.md`](KELLER_BECKMANN_BLACK_SPECIALIZATION.md)
studies the stronger requirement that the fixed map itself retain the
proper generic group `G`.

## 1. Three questions that must be separated

The fixed split-seed map has inverse polynomial

\[
 f_{\Pi,B,C}(T)
 =T^5-5T^3-2\Pi BT^2+4\Pi^3T-2\Pi^5C.                 \tag{1.1}
\]

There are three different notions:

1. A standard `G`-parametric extension is a regular Galois extension of
   `Q(t)` with group `G` whose specializations include every `G`-extension
   of `Q`.
2. A rational target curve
   `\gamma:P^1\dashrightarrow A^3_{\Pi,B,C}` gives such a one-parameter
   question only after taking the splitting field of
   `f_{\gamma(t)}` and checking that its generic group is `G`.
3. The full target of the fixed map is three-dimensional.  Its generic
   splitting group is `S_5`, so it is not a one-parameter extension and,
   for `G<S_5`, is not a `G`-extension at all.

The quintic fiber itself is a degree-five field.  Except for `C_5`, that
field is not Galois; `G` refers to its Galois closure in the natural
degree-five action.

These distinctions are standard in the definition of a parametric
extension; see
[Legrand, *Parametric Galois Extensions*](https://arxiv.org/abs/1310.6682)
and
[Legrand, *On parametric extensions over number fields*](https://arxiv.org/abs/1602.06706).

## 2. The unrestricted fixed-map minimum is one

The universal relative quadratic-gauge construction has two seed
parameters in degree five.  Promoting them to unchanged source and target
coordinates gives the single polynomial map

\[
 {\cal U}_5:\mathbb A^5_{\mathbb Q}\longrightarrow\mathbb A^5_{\mathbb Q}
\]

displayed explicitly in Section 3.1 of the universal note.  It has
determinant one and geometric degree five.  At target
`(u,v,\pi,b,c)`, its inverse polynomial is

\[
 v\pi^5S^5+u\pi^4S^4+\pi S^3+bS^2+S-\frac c2.         \tag{2.1}
\]

For every monic squarefree quintic `P(T)`, choose a rational origin `a`
away from the finitely many zeros of `P'P^{[3]}` and normalize
`P(a+S)/P'(a)`.  Its five non-linear coefficients give
`(u,v,\pi,b,c)` by the explicit formulas in the universal note, and the
complete fiber is `Spec(Q[T]/(P))`.

Therefore:

\[
\boxed{\text{minimum number of fixed Keller maps realizing every
quintic field}=1}                                      \tag{2.2}
\]

provided the ambient affine dimension is allowed to be five.  A minimum of
zero is impossible, so this is exact.

The interesting restricted question is now:

> Does one fixed degree-five Keller map of `A^3_Q` realize every quintic
> field, and if not, how many such threefold maps are required?

For the split-seed map, the only remaining intrinsic descent condition is
the rational-point problem on the Kummer threefold

\[
 {\cal V}_A:\quad
 \operatorname{Tr}(\eta)=0,\quad
 \operatorname{Tr}(\eta^2)=10,\quad
 \operatorname{Tr}(\eta^4)=50-16\Pi^3
\]

attached to the quintic algebra `A`; see Section 10 of
[`FIXED_QUINTIC_MODULI_DOMINANCE.md`](FIXED_QUINTIC_MODULI_DOMINANCE.md).
No local obstruction, universal rational point, or finite-cover theorem for
these threefolds is currently proved in the repository.

There is also a family-level strengthening.  Any monic separable quintic
`P_{\mathbf r}(T)` over a rational parameter variety can be normalized by
one rational origin and inserted coefficientwise into the same fixed
`A^5` map.  The pulled-back inverse cover is exactly
`Spec(P_{\mathbf r})`, so its splitting field is preserved.  In particular,
every two-parameter generic polynomial for a transitive quintic group gives
a rational target surface in this single fixed map.

For example,

\[
 T^5+rT^3+sT+s
\]

gives the `S_5` target surface

\[
 (u,v,\pi,b,c)
 =
 \left(0,\frac{s^4}{r^5},\frac rs,0,-2\right),
\]

and Brumer's generic `D_5` polynomial gives another small surface.
Section 3.2 of the universal note now also displays the Lecacheux `F_{20}`,
Buhler `A_5`, and Hashimoto--Tsunogai `C_5` surfaces explicitly and verifies
all five coefficient identities.  Thus the one fixed map contains
multi-parameter fiber-parametric families for all five transitive quintic
groups, although its ambient generic group is `S_5`.

## 3. What the fixed threefold map currently gives

The exact state is:

| group | family inside the fixed target | proved arithmetic output | one-parameter status |
|---|---|---|---|
| `S_5` | rational curve and dominant rational surface | infinitely many pairwise nonisomorphic fields; local refinements | no finite set of curves can be parametric |
| `F_{20}` | rational surface | infinitely many pairwise nonisomorphic fields; local refinements | open |
| `D_5` | rational curve | infinitely many pairwise nonisomorphic fields; local refinements | open |
| `A_5` | explicit oriented-discriminant and smaller descent incidences | certified fibers in both allowed real signatures | no finite set of curves can be parametric |
| `C_5` | explicit automorphism incidence and a cubic cover of the Lehmer line | one certified field | open |

The formulas and proofs are in
[`FIXED_QUINTIC_GALOIS_STRATIFICATION.md`](FIXED_QUINTIC_GALOIS_STRATIFICATION.md).
In particular, the existence of a rational curve with generic group `G`
and infinitely many `G`-specializations is much weaker than the assertion
that every `G`-extension is a specialization.

## 4. One rational space versus one rational curve

Hashimoto and Tsunogai constructed generic polynomials over `Q` with two
parameters for all five transitive subgroups of `S_5`; see
[their degree-five theorem](https://doi.org/10.3792/pjaa.79.142).
Consequently one rational space of dimension two suffices for each of
`C_5,D_5,F_{20},A_5,S_5`.

This is compatible with the classification of one-parameter generic
extensions.  Over a characteristic-zero field `k`, a one-parameter generic
extension exists only for the cyclic and odd-dihedral cases satisfying the
stated root-of-unity conditions.  Over `Q`, those conditions fail for
`C_5` and `D_5`, and `F_{20},A_5,S_5` are not on the list; see
[Dèbes--König--Legrand--Neftin, Theorem 2.5](https://arxiv.org/abs/2102.07465).
Failure of genericity does not by itself prove failure of
`Q`-parametricity.

Krashen and Neftin define `R`-equivalence on `H^1(Q,G)` by chains of
`G`-extensions over `Q(t)`.  Their odd-cyclic-kernel theorem gives

\[
 H^1(\mathbb Q,C_5)/R
 \cong H^1(\mathbb Q,1)/R,\quad
 H^1(\mathbb Q,D_5)/R
 \cong H^1(\mathbb Q,C_2)/R,
\]

and

\[
 H^1(\mathbb Q,F_{20})/R
 \cong H^1(\mathbb Q,C_4)/R.
\]

The right sides are trivial over `Q`, so `C_5,D_5,F_{20}` are
`R`-trivial.  The two-parameter generic polynomials similarly make
`A_5` and `S_5` `R`-trivial.  Thus the `R`-equivalence count is one for all
five groups.  It counts rational parameter spaces, not one-parameter
curves.  See
[Krashen--Neftin, *On rational connectedness and parametrization of finite Galois extensions*](https://arxiv.org/abs/2502.15674).

## 5. The decisive local obstruction for `A_5` and `S_5`

Let `E_1/Q(t),...,E_r/Q(t)` be any finite collection of regular
`G`-extensions, with `G=A_5` or `S_5`.  There are infinitely many primes
`p` which split in the residue fields of every branch point of every
`E_i`.  Outside a finite exceptional set, every completion at such a prime
of every specialization of every `E_i` has cyclic decomposition group.

Both `A_5` and `S_5` contain `V_4=C_2\times C_2`.  The required weak
Grunwald statements produce global `G`-extensions whose completion at one
of these primes has group `V_4`.  Such an extension cannot specialize from
any `E_i`.  Hence

\[
\boxed{\text{no finite collection of rational one-parameter families is
parametric for }A_5\text{ or }S_5.}                    \tag{5.1}
\]

For `A_5`, the published proof uses Mestre's specialization theorem.  For
`S_5`, the generic polynomial supplies the needed Grunwald property.
The general cyclic-versus-`V_4` criterion is Theorem 7.2 and its examples
in
[König--Legrand--Neftin, *On the local behaviour of specializations of function field extensions*](https://arxiv.org/abs/1709.03094).

This exactly implements the proposed “find a local extension that cannot
occur” attack, but only for a finite collection of target curves.  It
cannot obstruct the full five-dimensional universal map, and it does not
by itself obstruct the full three-dimensional coefficient pencil.

A naive finite-field test on the projective descent model also cannot prove
the latter obstruction: the `(2,4)` compactification always contains the
singular rational boundary point
`[\eta:s:W]=[0:0:1]`.  More strongly, the blowup calculation in Section
10.5.1 of
[`FIXED_QUINTIC_MODULI_DOMINANCE.md`](FIXED_QUINTIC_MODULI_DOMINANCE.md)
shows that this boundary produces primitive realization points over
`\mathbb Q_p` at every sufficiently large good prime.  It uses an isotropic
direction `e` for the trace quadric with
`\operatorname{Tr}(e^4)\ne0`, followed by the exact scaling

\[
 \eta=\lambda e+\lambda^7mf,\qquad s=\lambda^4\sigma.
\]

The rescaled two-equation Jacobian is a unit, so the point lifts into
`sW\ne0`.  Thus any local counterexample must be sought among the finite
small, ramified, or otherwise exceptional primes for the chosen field;
simply finding no affine point modulo a large good prime is insufficient.

The groups `C_5`, `D_5`, and `F_{20}` have no `V_4` subgroup, so this
particular local argument is unavailable.  Their `R`-triviality points
toward rational spaces of dimension at least two, while one-parameter
`Q`-parametricity remains a separate problem.

## 6. A Keller version of `R`-equivalence

Fix a Keller map `F:A^m_Q->A^m_Q` and its squarefree degree-five target
open `U`.  For quintic fields occurring as complete fibers, define direct
Keller linkage by

\[
 K_0\sim_F^{\,1}K_1
\]

if a rational map `\gamma:P^1\dashrightarrow U` has endpoint fibers
`K_0,K_1` and its generic inverse polynomial is irreducible.  Let
`\sim_F` be the equivalence relation generated by direct linkages.  If the
generic splitting group along the curve is `G`, taking splitting fields
maps a Keller linkage to classical `R`-equivalence in `H^1(Q,G)`.

The universal map `{\cal U}_5` is compatible with arbitrary monogenic
one-parameter quintic families: over `Q(t)`, choose an origin away from the
two jet divisors and apply the same coefficient formulas as in (2.1).
Thus every such family lifts to a rational target curve of this one fixed
map.  Subject to choosing a model and primitive generator regular at the
endpoints, `{\cal U}_5` reproduces the classical rational linkage rather
than introducing a new arithmetic obstruction.

For the split-seed map of `A^3`, the relation is genuinely finer: a
classical rational family lifts only when its primitive generators solve
the square/cube descent conditions, equivalently the Kummer condition on
`\mathcal V_A`.

## 7. Corrected headline targets

The unrestricted headline is settled:

> **One explicit Keller map of `A^5_Q` realizes every quintic field.**

Two sharper frontiers remain:

1. **Threefold universality.**  Decide whether the fixed split-seed map of
   `A^3` realizes every quintic field by proving or disproving rational
   solubility of every `\mathcal V_A`.
2. **Curve-parametricity.**  No finite collection can work for `A_5` or
   `S_5`.  For `C_5,D_5,F_{20}`, determine whether a rational target curve
   can be `Q`-parametric, or prove non-parametricity by a different local
   invariant.

The most concrete next computations are therefore:

- search for a `Q_p`-point obstruction on `\mathcal V_A` at the finite set
  of small and ramified primes of explicit quintic fields not already
  realized by the split-seed map;
- compute the branch-cycle and residue-field data of the explicit `D_5`
  curve and selected curves on the `F_{20}` surface;
- use the three free trace-quadric parameters in the `C_5` incidence before
  imposing the cube, rather than only the current Lehmer section;
- search the smaller `A_5` descent surface for rational curves, while
  recognizing that no finite list of resulting curves can be globally
  parametric.
