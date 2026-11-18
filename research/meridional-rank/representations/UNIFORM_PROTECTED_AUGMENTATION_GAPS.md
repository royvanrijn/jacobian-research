# Uniform protected augmentation gaps and prime six-bridge families

Date: 18 November 2026.

The root/covector obstruction and alternating colorings admit one formula in every bridge number. This note gives the uniform two-cusped parent construction and a direct bridge-six certificate. The [threaded-curl transfer theorem](UNBOUNDED_PRIME_AUGMENTATION_GAP_FAMILIES.md) now identifies all parents with fillings of one fixed hyperbolic manifold, proving infinite prime augmentation-gap families at every sufficiently large bridge number.

Write

\[
U=\sigma_1\sigma_2\sigma_3^{-2}\sigma_2\sigma_3^{-2}
       \sigma_2\sigma_1^2\sigma_2^3\sigma_3^{-1}.
\]

For `b>=5`, put `d=b-2` and

\[
\begin{split}
\rho_b&=\sigma_{b-1}\sigma_{b-2}\cdots\sigma_4,\\
\delta_b&=\rho_b^{-1}(\sigma_2\sigma_3\cdots\sigma_{b-2})\rho_b\\
&=\sigma_2(\sigma_3\sigma_4\sigma_3^{-1})
                   \sigma_5\cdots\sigma_{b-1},\\
r_{b,k}&=(d+1)(1+d^2k),\qquad
K_{b,k}=\widehat{U\delta_b^{r_{b,k}}}\quad(k\ge0).
\end{split}
\]

The second expression for `delta_b` follows from the adjacent braid relation and commuting nonadjacent generators; equivalently `rho_b` identifies the protected generators with the adjacent generators on positions two through `b-1`.

**Theorem 1 (uniform augmentation gap).** For every `b>=5` and `k>=0`, `K_(b,k)` is a knot and

\[
\operatorname{ar}_{\mathbb C}(K_{b,k})\le b-1<
\operatorname{mr}(K_{b,k})=\operatorname{br}(K_{b,k})=b.
\]

Its group maps onto `A_(2b)` with every meridian a double transposition. This statement makes no unrestricted primeness assertion.

**Theorem 2 (prime six-bridge tail).** In particular, put

\[
\delta_6=\sigma_2\sigma_3\sigma_4\sigma_3^{-1}\sigma_5,
\qquad K_{6,k}=\widehat{U\delta_6^{5+80k}}.
\]

All sufficiently large `k` give hyperbolic prime knots. Infinitely many different knot types occur among these prime members, and every member has

\[
\operatorname{ar}_{\mathbb C}(K_{6,k})\le5<
\operatorname{mr}(K_{6,k})=\operatorname{br}(K_{6,k})=6.
\]

No effective exceptional threshold or pairwise distinctness for every parameter is asserted. The corresponding prime-tail result at bridge five is proved in [PRIME_FIVE_BRIDGE_AUGMENTATION_GAP_FAMILY.md](PRIME_FIVE_BRIDGE_AUGMENTATION_GAP_FAMILY.md).

## 1. Explicit alternating seeds in every degree

On `2b` symbols, use

\[
\begin{split}
x_0&=(01)(23),&x_1&=(14)(23),\\
x_2&=(05)(16),&x_3&=(12)(67),\\
x_4&=(08)(79),&&\\
x_j&=(2j-2\ \ 2j)(2j-1\ \ 2j+1)&& (5\le j\le b-1).
\end{split}
\]

The first four generate `A_8`. Each later generator has two distinct old symbols and two new ones. The two-new-symbol extension lemma in [ALTERNATING_TARGET_PRIME_FAMILIES.tex](ALTERNATING_TARGET_PRIME_FAMILIES.tex), Theorem `thm:targets`, therefore proves that the full tuple generates `A_(2b)`, for every `b`.

Use the positive Hurwitz rule `T(a,b)=(b,b^{-1}ab)`. At the end of `U`, the first four labels are

\[
x_0,\quad a=(05)(14),\quad h=(16)(23),\quad C=x_3;
\]

the new labels are unchanged. Frame the middle tuple as

\[
(a,h,z_1,\ldots,z_{d-2}),\qquad
z_i=Cx_{i+3}C^{-1}.
\]

The protected-tail identity in the audit identifies `delta_b` with the ordinary Hurwitz Coxeter move on this `d`-tuple. Its ordered product is explicitly

\[
p=(0\ 8\ 10\ \cdots\ 2b-2\ 5)
  (1\ 6\ 9\ 11\ \cdots\ 2b-1\ 4)(2\ 3)(7).
\]

These disjoint cycles have lengths `d,d+1,2,1`; since `d(d+1)` is even, the order is exactly `d(d+1)`.

Put `X=x_1`, `Y=x_2`. The cycle formula gives

\[
pXp^{-1}=h,\qquad z_i=p^iYp^{-i},\qquad
p^dYp^{-d}=a.
\]

For example, `p` sends the endpoints of `X=(14)(23)` to those of `(16)(23)`. The first new framed label is `(08)(69)`, and each successive iterate introduces the next two symbols. The last equality follows because `p^d` fixes the first `d`-cycle and acts by one backward step on the `(d+1)`-cycle, sending `(16)` to `(14)`.

For any `d`-tuple `(a_1,...,a_d)` with product `p`, the Coxeter move `R=T_1...T_(d-1)` satisfies

\[
R(a_1,\ldots,a_d)=(a_2,\ldots,a_d,p^{-1}a_1p).
\]

It follows by iteration that `R^d` conjugates every entry by `p^{-1}`, and

\[
R^{d+1}(a_1,\ldots,a_d)=
(p^{-1}a_2p,\ldots,p^{-1}a_dp,p^{-2}a_1p^2).
\]

Applying this to the displayed prefix tuple and using the three identities yields exactly

\[
(X,Y,z_1,\ldots,z_{d-2}),
\]

which is the initial framed tuple. Hence `U delta_b^(d+1)` returns the full alternating seed. Moreover `R^(d^2(d+1))` fixes the prefix tuple because `p` has order `d(d+1)`. Thus return holds for every exponent

\[
(d+1)+d^2(d+1)k=r_{b,k}.
\]

This proves the full-alternating quotient uniformly, without enumerating any growing alternating group.

## 2. Knot permutation, bridge number and augmentation

The permutation of `U` is `(1 2)(3 4)` in one-based strand labels. The permutation of `delta_b` is a `d`-cycle on the protected positions `{2,3,5,...,b}`, fixing positions one and four. Since `r_(b,k)` is congruent to one modulo `d`, its permutation does not change with `k`. Multiplying by the two transpositions of `U` inserts the two outer fixed positions into that cycle, at positions two and three respectively. The resulting permutation is a `b`-cycle, proving that every closure is a knot.

At most `b-1` double transpositions have at most `2b-2` transposition edges on `2b` symbols. Their graph is disconnected, and their generated group cannot be transitive. The full alternating quotient therefore forces `mr>=b`; the `b`-strand braid gives `mr<=br<=b`.

The expression for `delta_b` belongs to the protected subgroup generated by

\[
\sigma_2,\quad \sigma_3\sigma_4\sigma_3^{-1},
\quad\sigma_5,\ldots,\sigma_{b-1}.
\]

The all-field paired-span Theorem 3 of [UNBOUNDED_PRIME_AUGMENTATION_GAP_AUDIT.md](UNBOUNDED_PRIME_AUGMENTATION_GAP_AUDIT.md) applies to every power. It gives `ar_C<=b-1`. This proves Theorem 1. Its augmentation inference uses the full root/covector obstruction and Cornwell's correspondence, not just finite permutation class rank.

## 3. One marked parent for each bridge number

Let

\[
t_b=\sigma_2\cdots\sigma_{b-2},\qquad
Q_b=\rho_b U\rho_b^{-1}t_b^{d+1}.
\]

The closure of `Q_b t_b^(d q)` is the closure of `U delta_b^(d+1+d q)` by braid conjugation. Add one last strand and define

\[
P_b=\sigma_b\sigma_{b-1}\cdots\sigma_3\sigma_2^2
                \sigma_3\cdots\sigma_{b-1}\sigma_b^{-1}.
\]

Let `M_b` be the exterior of the two-component link `closure(Q_b P_b)`. The added last component is a circle bounding a disk pierced by the `d` consecutive positions two through `b-1`.

Here is a uniform geometric description, including the slope marking. Factor

\[
P_b=A_{2,b+1}A_{3,b+1}\cdots A_{b-1,b+1},
\qquad
A_{i,b+1}=\sigma_b\cdots\sigma_{i+1}\sigma_i^2
                       \sigma_{i+1}^{-1}\cdots\sigma_b^{-1}.
\]

These are the last-puncture point-pushing loops around the indicated adjacent punctures with a common stem. Closing the last strand and pulling its cap down the stem gives the oval around precisely those `d` strands. Its detour across the outside `b`-th strand has opposite crossing signs with the circle over in both crossings, and cancels. The remaining circle has no self crossing and is over all `d` selected strands on the outward trip and under them on return. It bounds the transverse disk, with linking number `d` and its preferred zero-framing longitude supplied by that disk.

For the group sign control, write the slab input as `(x,a_1,...,a_d,e,c)` and put `w=a_1...a_d`. Exactly as in the three-strand disk proof in the five-bridge note, the output is

\[
\left(x,c^{-1}a_1c,\ldots,c^{-1}a_dc,
 c_{\rm out}c^{-1}ec\,c_{\rm out}^{-1},c_{\rm out}\right),
\qquad c_{\rm out}=(wc)^{-1}c(wc).
\]

This also follows by inducting through the displayed point-pushing factors. The circle closure gives `[c,w]=1` and restores `e`. Circle filling `mu-q lambda` imposes `c=w^q`, hence conjugates the selected meridians by `w^{-q}`. A positive `d`-strand full twist `t_b^d` has that action for `q=1`. The disk-twist homeomorphism of the unknot exterior identifies the filling geometrically with that full twist and preserves the old knot meridian. Therefore, for every integer `q`,

\[
M_b(1,-q)\cong
S^3\setminus\widehat{Q_b t_b^{d q}}.
\]

In particular

\[
M_b(1,-d(d+1)k)\cong S^3\setminus K_{b,k}.
\]

**Conditional consequence.** If the specified parent `M_b` is hyperbolic, then all sufficiently large `k` give hyperbolic prime `K_(b,k)`, and infinitely many knot types occur. The proof is the same fixed-parent filling and strict-volume-squeeze argument as in Section 5 of the five-bridge note: the slope lengths diverge; Thurston's Theorem 5.8.2 gives hyperbolicity; Futer–Kalfagianni–Purcell Theorem 1.1 and Thurston Theorem 6.5.6 force the filled volumes to approach the parent volume from below. Mostow–Prasad rigidity then excludes a finite collection of exterior types. Exact primary source statements and versions are pinned in the new ledger.

The conditional statement is completed uniformly by the [threaded-curl transfer theorem](UNBOUNDED_PRIME_AUGMENTATION_GAP_FAMILIES.md): the parents are hyperbolic at every sufficiently large bridge number. The argument uses a fixed three-cusped parent and a uniform diagram isotopy.

## 4. The six-bridge parent is hyperbolic

For `b=6`, `d=4`, the base braid and circle suffix are

\[
Q_6=\sigma_5\sigma_4 U\sigma_4^{-1}\sigma_5^{-1}
                     (\sigma_2\sigma_3\sigma_4)^5,
\qquad
P_6=\sigma_6\sigma_5\sigma_4\sigma_3\sigma_2^2
                        \sigma_3\sigma_4\sigma_5\sigma_6^{-1}.
\]

Their seven-strand closure has a connected spherical 43-crossing PD projection. The added circle component has ten edges and no self crossings; the other component has 76 edges. The oriented linking number is four. Its complete two-cusped exterior has a retained 34-tetrahedron triangulation with decorated signature

```text
ILLvvPvLzLQAzLMPzQAAAMQkbeilmotrssstnqvuuBAzwADDzEDCFHGHGHGhhsaaalaaltsvbafflhgoorsgancqldwswr_badbCBBbrzRy
```

The literal braid is bound to the PD by all incoming/outgoing ports and all seven closure identifications. A fresh exterior constructed from this PD matches the retained triangulation. With the pinned Sage 10.9/SnapPy 3.3.2/Spherogram 2.4.1 environment, both that exterior and the independently reloaded triangulation pass rigorous Krawczyk interval hyperbolicity verification at 100 bits, without randomization. The retained face pairings, exact rational interval bounds, logarithmic gluing and completeness equations and positive orientations are checked. This proves `M_6` hyperbolic under the stated software imports.

The filling slopes are `(1,-20k)`. The conditional consequence therefore proves Theorem 2. The approximate parent volume `29.9270860877851` is descriptive and is not used as a numerical inequality.

## 5. Evidence and limits

The new [checker](verify_protected_six_bridge_uniform_family.py) checks the uniform seed, product, return and knot-permutation formulas for `4<=b<=24` as bounded controls of the written all-`b` proof. Its exact six-bridge controls bind the PD, circle-component identification, point-pushing factorization and universal surgery action. It rejects an altered detour and a swapped circle component. These computations do not replace the written proof of the unrestricted formulas.

The [input](protected_six_bridge_uniform_family_inputs.json) retains the full literal two-component diagram, component marking, complete triangulation, face pairings, logarithmic equations and rational interval boxes. The [finite receipt](protected_six_bridge_uniform_family_verification.json), [rigorous geometry receipt](protected_six_bridge_uniform_family_geometry_verification.json) and [source ledger](protected_six_bridge_uniform_family_source_ledger.json) state their exact scopes and hashes.

```bash
.venv/bin/python research/meridional-rank/representations/verify_protected_six_bridge_uniform_family.py \
  --write-receipt research/meridional-rank/representations/protected_six_bridge_uniform_family_verification.json
~/.local/share/jacobian-sage-10.9/bin/python \
  research/meridional-rank/representations/verify_protected_six_bridge_uniform_family.py --sage \
  --write-receipt research/meridional-rank/representations/protected_six_bridge_uniform_family_geometry_verification.json
```

The rigorous computation is an internal interval replay with pinned source and kernel editions. No formal verification, independent replay, external review or novelty judgment has been performed. The [unbounded-bridge theorem](UNBOUNDED_PRIME_AUGMENTATION_GAP_FAMILIES.md) proves the prime-family result for every sufficiently large bridge number; it does not compute the finite exceptional bridge set. The direct bridge-five and bridge-six certificates are independent of that non-effective threshold.
