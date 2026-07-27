# Adversarial audit of universal Keller-fiber multiplicity

## Verdict

The theorem in
[UNIVERSAL_KELLER_FIBER_MULTIPLICITY.md](UNIVERSAL_KELLER_FIBER_MULTIPLICITY.md)
survives the checks below:

\[
\boxed{
\begin{aligned}
&|\mathcal R_K(A)|=\infty
&&\text{for number fields }K\text{ and ranks }N\ge4,\\
&|\mathcal R_K(A)|=\infty
&&\text{for characteristic-zero }K\text{ and ranks }N\ge5.
\end{aligned}}
\]

The audit found and corrected one field-extension issue in the quartic clean
locus: the quadratic

\[
 4\alpha^2+4\alpha+3
\]

has no rational root but can have roots over a number field containing
`\sqrt{-2}`.  Those are only two additional excluded parameter values, so
the infinitude proof is unchanged.

The theorem is not independent of the repository's stable-normalization
machinery.  Its stable-inequivalence conclusion imports:

1. exact stable orbit classification on the quadratic-gauge coefficient
   torus;
2. intrinsic weighted boundary selection and selected-root Torelli.

Everything else in the proof is elementary polynomial algebra, standard
quadratic-form arithmetic, or exact full-fiber reconstruction.

## 1. Attack on the generator choice

### Possible failure

The proof requires a primitive trace-zero generator `eta` satisfying
`Tr(eta^2)!=0`.  It would fail if the primitive open were disjoint from the
nonisotropic open.

### Resolution

For a finite etale algebra `A/K`, primitive elements form a nonempty Zariski
open.  Translating a primitive element by `Tr(theta)/N` preserves the
generated algebra and puts it in

\[
 A_0=\ker(\operatorname{Tr}_{A/K}).
\]

The trace pairing restricts nondegenerately to `A_0`, because

\[
 A=K\cdot1\perp A_0,\qquad\operatorname{Tr}(1)=N\ne0.
\]

Thus `Tr(eta^2)=0` is a proper quadric in `A_0`.  Both complements are
nonempty opens in one affine space.  Since every characteristic-zero field
is infinite, their intersection has a `K`-point.

Newton's identity then gives

\[
 c_{N-2}=-\frac12\operatorname{Tr}(\eta^2)\ne0.
\]

No Hilbertian or local-global hypothesis enters this step.

## 2. Attack on the clean coefficient torus

### Possible failure

Translation might force one derivative coefficient to vanish identically,
leaving the locus on which the stable-moduli theorem applies.

### Resolution

For

\[
 G_s(S)=P(s+S)-P(s)=\sum_{j=1}^Ng_j(s)S^j,
\qquad
 g_j(s)=\frac{P^{(j)}(s)}{j!},
\]

every `g_j` is a nonzero polynomial: its leading term comes from the monic
term `T^N`, and characteristic zero prevents its coefficient from vanishing.
Hence

\[
 g_1g_3\cdots g_N
\]

is a nonzero polynomial.  Only finitely many translations are excluded.
The omitted `g_2` is exactly the coefficient removed by a polynomial target
shear.

This verifies all hypotheses of the clean coefficient-torus theorem, not
only `g_1g_3g_N!=0`.

## 3. Attack on the stable invariants

### Quintic

The normalized coefficient weights are

\[
 w_3=(-2,-1),\quad w_4=(-3,-4),\quad w_5=(-4,-5).
\]

Their primitive relation is

\[
 -w_3-6w_4+5w_5=0,
\]

so

\[
 I=\frac{a_5^5}{a_3a_4^6}
\]

is invariant under both independent source--target scalings.  For a centered
quintic,

\[
 g_3=10s^2+c_3,\qquad g_4=5s,\qquad g_5=1.
\]

Since `c_3!=0`, `g_3` is a unit at zero.  The coefficient of `s^2` in
`g_1` is `3c_3`, so `ord_0(g_1)<=2`.  Therefore

\[
 I(s)=\frac{g_1(s)^2}{g_3(s)g_4(s)^6}
\]

has order at most `4-6=-2` at zero.  Cancellation cannot make it constant.

### Degrees at least six

For `j>=4`,

\[
 w_j=(1-j,-j).
\]

The relation

\[
 w_{N-2}+w_N-2w_{N-1}=0
\]

is exact.  Normalization by `g_1` cancels from the corresponding monomial,
giving

\[
 J_N=\frac{a_{N-2}a_N}{a_{N-1}^2}
    =\frac{g_{N-2}g_N}{g_{N-1}^2}.
\]

The top three derivative jets are

\[
 g_N=1,\qquad g_{N-1}=Ns,\qquad
 g_{N-2}=\binom N2s^2+c_{N-2},
\]

and therefore

\[
 J_N(s)=\frac{N-1}{2N}+\frac{c_{N-2}}{N^2s^2}.
\]

The second summand is nonzero.  Thus `J_N` is nonconstant even after every
lower coefficient is specialized adversarially.

The argument needs only a separating invariant: it does not claim that this
single monomial generates the whole stable quotient when `N>5`.

## 4. Attack on infinitude after deleting exceptions

A nonconstant rational function on `\mathbb A^1_K` over an infinite field
has infinite image.  Otherwise a finite product of its differences from the
finitely many image values would be the zero rational function.

Deleting the finite zero set of `g_1g_3...g_N`, together with poles, cannot
make that image finite.  Thus infinitely many clean translations have
pairwise distinct invariant values.

The stable-moduli theorem is geometric over an algebraic closure.  If two
maps were stably equivalent over `K`, they would be equivalent after scalar
extension.  Distinct values of `I` or `J_N` therefore obstruct equivalence
over `K` as well.

## 5. Attack on the common full fiber

At

\[
 y_s=\left(1,0,-\frac{2P(s)}{g_1(s)}\right),
\]

the inverse equation is identically

\[
 G_s(S)-\frac{g_1(s)}2
       \left(-\frac{2P(s)}{g_1(s)}\right)
=P(s+S).
\]

Because `eta` is primitive, its characteristic polynomial has degree `N`;
because `A` is etale, it is squarefree.  Exact quotient-ring reconstruction
therefore identifies the entire scheme fiber with

\[
 \operatorname{Spec}K[S]/(P(s+S))
 \simeq\operatorname{Spec}A.
\]

The inverse polynomial has degree `N`, equal to the geometric degree of the
map.  Thus the fiber is full, not only a collection of selected inverse
points.  The fixed target scaling that changes determinant `-2` to `1`
leaves the displayed second target coordinate zero and changes none of these
identities.

## 6. Attack on the quartic local-global step

The quartic trace-chord form

\[
 \mathcal Q_A=\operatorname{Tr}_{A/K}(\eta^2)-2e^2-4u^2
\]

is nondegenerate of dimension five.  At a real place, the three possible
signatures are

\[
 (3,2),\qquad(2,3),\qquad(1,4),
\]

so it is isotropic.  Complex places are automatic.  At a nonarchimedean
completion of a number field, the `u`-invariant is four, so every
five-dimensional form is isotropic.  Hasse--Minkowski supplies a global
`K`-point.

The resulting smooth projective threefold quadric is `K`-rational, hence has
Zariski-dense `K`-points.  None of the excluded loci contains the quadric:

- `e=0`, `d=e-2u=0`, and fixed `alpha` values are hyperplane sections;
- nonprimitive elements lie in a finite union of proper finite-etale
  subalgebras, hence proper linear subspaces in `A_0`;
- a nondegenerate rank-five quadric cannot be contained in any such proper
  linear subspace.

The normalized seed parameter

\[
 \alpha=u/e-\frac12
\]

is nonconstant, because otherwise the quadric would lie in one hyperplane.
It therefore assumes infinitely many values on the good open.

This arithmetic step is genuinely number-field-specific.  The explicit
anisotropic example in
[LOW_RANK_MULTIPLICITY_BOUNDARIES.md](LOW_RANK_MULTIPLICITY_BOUNDARIES.md)
shows that the same trace-chord argument does not work over every
characteristic-zero field.

## 7. Attack on quartic stable separation

The normalized weighted seed is

\[
 H_\alpha(W)=W^2(W-1)(\alpha W-\alpha-1).
\]

On the declared clean open, the full Hessian Fitting divisor and the marked
boundary points determine the seed up to rerooting at a nonzero primitive
root.  The intrinsic affine root-one sheet kills every nontrivial rerooting.
Thus stable equivalence recovers the exact normalized seed, not merely its
unmarked Hessian divisor.

Consequently two different values of `alpha` cannot become equivalent after
an affine root change, a polynomial left--right transformation, or identity
stabilization.  This is precisely where the affine mark is load-bearing.
Without it, the quartic Hessian reflection would identify rerooted seeds and
the argument would fail.

## 8. Dependency and evidence ledger

| step | status | source |
|---|---|---|
| primitive centered generator | elementary proof | Sections 1--2 above |
| translated inverse identity | exact algebra | universal multiplicity checker |
| full scheme fiber | theorem | finite-etale reconstruction |
| quintic invariant | exact algebra | quintic checker |
| all-`N>=6` invariant | exact algebra | higher-degree checker |
| quadratic stable separation | theorem | quadratic-gauge stable moduli |
| quartic isotropy over number fields | external arithmetic theorem | Hasse--Minkowski and local `u=4` |
| weighted stable separation | theorem | decorated normalization and selected-root Torelli |
| degree `4,5,6` concrete regressions | exact algebra | multiplicity witness cards |

No bounded search is used to prove a universal quantifier.  The witness
cards are regressions and exposition; the universal proof remains the
symbolic and arithmetic argument above.

## 9. Reproduction

Run

```bash
.venv/bin/python scripts/verify_quadratic_gauge_stable_moduli.py
.venv/bin/python scripts/verify_universal_quartic_fiber_multiplicity.py
.venv/bin/python scripts/verify_universal_quintic_fiber_multiplicity.py
.venv/bin/python scripts/verify_universal_higher_degree_fiber_multiplicity.py
.venv/bin/python scripts/verify_universal_multiplicity_witness_cards.py
.venv/bin/python scripts/verify_low_rank_multiplicity_boundaries.py
```
