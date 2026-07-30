# Two marked fibers recover the power-shift gauge

The whole-plane stable-multiplicity theorem shows that one complete inverse
plane at `P=1` does not determine a Keller map: every power-shift exponent
gives the same marked cover there while the stable Fitting area changes.
The next reconstruction step is sharp.  A second marked fiber at a
non-torsion value of `P` recovers the exponent.  At a root of unity it does
not, and the residual congruence ambiguity is exact.

This is a theorem about the normalized quadratic-gauge power-shift family,
not an unrestricted Torelli theorem for arbitrary Keller maps.

## 1. Marked-fiber data

Let `K` have characteristic zero, let `N>=4`, and let

\[
 G(S)=\sum_{j=1}^N g_jS^j,\qquad
 g_1g_3g_4\cdots g_N\ne0.
 \tag{1.1}
\]

For `m>=0`, the power-shift inverse equation is

\[
\begin{aligned}
 E_m(P,B,C;S)
 ={}&g_1S+P(g_2S^2+g_3S^3)
       +\sum_{j=4}^Ng_jP^{j+m}S^j\\
    &-\frac{g_1}{2}(BS^2+C).
\end{aligned}
\tag{1.2}
\]

Fix a target `y=(p,B,C)` with `p!=0` for which (1.2) is squarefree.
The corresponding **root-marked fiber** is

\[
 \left(A_{m,y},s_{m,y}\right)
 =
 \left(K[S]/(E_m(p,B,C;S)),\ S\bmod E_m\right).
 \tag{1.3}
\]

The mark here retains the distinguished inverse-root generator, not merely
the abstract finite etale algebra or one primitive idempotent after
splitting.  Consequently its monic annihilator

\[
 M_{m,y}(S)=\frac{E_m(p,B,C;S)}{g_Np^{N+m}}
 \tag{1.4}
\]

is intrinsic to the marked pair (1.3).  Squarefreeness ensures that this is
the full annihilator even when the fiber algebra is disconnected.

Write

\[
 \lambda_{m,p}=[S]M_{m,y}.
 \tag{1.5}
\]

Neither `B` nor `C` changes the linear or leading coefficient of (1.2), so

\[
 \boxed{
 \lambda_{m,p}=\frac{g_1}{g_Np^{N+m}}.
 }
 \tag{1.6}
\]

This one coefficient is the entire reconstruction calculation.

## 2. Two-fiber reconstruction theorem

At `P=1`, the monic marked polynomial is

\[
 M_{m,(1,B,C)}
 =\frac{1}{g_N}
 \left(G(S)-\frac{g_1}{2}(BS^2+C)\right),
 \tag{2.1}
\]

independent of `m`.  It recovers the normalized seed `G/g_N`: its linear
coefficient is `g_1/g_N`, its quadratic coefficient becomes `g_2/g_N`
after adding `B(g_1/g_N)/2`, and all coefficients of degrees at least three
are already `g_j/g_N`.  Overall scaling of `G` does not change the
quadratic-gauge map, which depends on the ratios `g_j/g_1`.

Now add any squarefree marked fiber at `P=c`, where `c in K^times`.
The two linear coefficients satisfy

\[
 \boxed{
 \frac{\lambda_{m,1}}{\lambda_{m,c}}=c^{N+m}.
 }
 \tag{2.2}
\]

Thus the pair first recovers the normalized seed and then recovers the
character value

\[
 c^m=c^{-N}\frac{\lambda_{m,1}}{\lambda_{m,c}}.
 \tag{2.3}
\]

### Theorem 2.1

Let `c in K^times` have infinite multiplicative order.  On the normalized
degree-`N` power-shift locus, one root-marked squarefree fiber over `P=1`
and one over `P=c` recover the seed up to its irrelevant common scalar and
recover the exponent `m` exactly.  Hence they recover the stable
power-shift class and its normalized Fitting Newton area

\[
 \boxed{2N-3+(N-2)m.}
 \tag{2.4}
\]

The targets on the two planes may have different `B,C` values; formula
(1.6) is independent of them.

The choice

\[
 \boxed{c=2}
 \tag{2.5}
\]

works over every characteristic-zero field.  Indeed its prime field
contains a copy of `Q`, in which no positive power of `2` is one.

### Corollary 2.2

Adding the whole inverse plane `P=c` is more data than necessary: any one
squarefree root-marked fiber on it already supplies (2.2).  Accordingly the
two framed inverse planes `P=1,c` recover `m` exactly when `c` is
non-torsion.  For the restricted polynomial maps themselves, the same
criterion follows from

\[
 F_m=(P,B_{\mathrm{low}}+P^mH_B,
         C_{\mathrm{low}}+P^mH_C),
\tag{2.6}
\]

because the nonzero higher decoration is multiplied on `P=c` by `c^m`.

For reconstruction in the full normalized family, two is minimal in the
following sense.  The `P=1` marked fiber recovers the seed but is identical
for every `m`.  A single fiber at `P=c` with the seed unknown confounds the
factor `c^m` with the seed: replacing `m` by `m'` and replacing every
`g_j`, `j>=4`, by `c^{m-m'}g_j` leaves the specialized equation at `P=c`
literally unchanged.  The pair at `1,c` removes exactly that ambiguity.

## 3. Sharp torsion counterexamples

Suppose `c` has finite multiplicative order `d`.  For a fixed seed,

\[
 E_m(c,B,C;S)=E_{m'}(c,B,C;S)
 \quad\Longleftrightarrow\quad
 m\equiv m'\pmod d.
 \tag{3.1}
\]

The forward implication follows from the nonzero top coefficient; the
reverse implication is immediate from (1.2).  The fibers at `P=1` also
agree.  Therefore the two marked fibers retain infinitely many stable
classes:

\[
 m,\ m+d,\ m+2d,\ldots,
 \tag{3.2}
\]

whose Fitting areas are strictly increasing.

More generally, after the seed has been recovered at `P=1`, sample finitely
many supplementary nonzero planes `P=c_1,...,c_r`.  If every `c_i` is
torsion of order `d_i`, the complete marked data are unchanged after

\[
 m\longmapsto m+\operatorname{lcm}(d_1,\ldots,d_r).
 \tag{3.3}
\]

If at least one `c_i` has infinite order, that supplementary sample
separates the exponent.  Thus the finite-plane criterion for a recovered
seed is exact:

\[
\boxed{
\text{the exponent is recovered iff }
\mathbb Z\longrightarrow(K^\times)^r,\quad
n\longmapsto(c_1^n,\ldots,c_r^n)
\text{ is injective},
}
\tag{3.4}
\]

equivalently, iff some `c_i` is non-torsion.

The value `c=0` is excluded: the displayed inverse equation loses degree
there and does not give a rank-`N` finite etale fiber.

## 4. One transverse inverse-cover line

There is an equally sharp one-family version.  Let

\[
 \ell(u)=(p(u),B(u),C(u)),\qquad p(u)=a+bu,\quad b\ne0,
 \tag{4.1}
\]

be an affine target line transverse to the planes `P=constant`.  Delete
`p=0` and the repeated-root discriminant.  On the resulting open line the
root-marked inverse cover has monic annihilator whose linear coefficient is

\[
 \lambda_m(u)
 =\frac{g_1}{g_Np(u)^{N+m}}.
 \tag{4.2}
\]

At the unique point `u_0=-a/b` where the completed line meets `P=0`,

\[
 \boxed{\operatorname{ord}_{u_0}\lambda_m=-(N+m).}
 \tag{4.3}
\]

### Theorem 4.1

One root-marked transverse inverse-cover family, together with its known
base map to the `P`-line, recovers the gauge exponent as

\[
 m=-\operatorname{ord}_{u_0}\lambda_m-N.
 \tag{4.4}
\]

In particular two distinct power shifts cannot have the same root-framed
cover on a transverse line.  Geometrically, the transverse family sees the
degree-drop valuation at `P=0` that the common plane `P=1` misses.

This does not say that the abstract etale cover of the punctured line,
after forgetting the inverse-root generator, the `P`-coordinate, and the
missing-point valuation, determines `m`.  Those forgetful variants remain
separate reconstruction questions.

## 5. Reconstruction hierarchy

For this family the proposed hierarchy now has an exact first answer.

| retained data | power-shift exponent |
|---|---|
| abstract finite algebra on `P=1` | invisible |
| root-marked fiber on `P=1` | invisible, although the normalized seed is recovered |
| full root-marked inverse plane `P=1` | invisible |
| marked fibers on `P=1` and `P=c`, `c` non-torsion | recovered |
| any finite set of torsion planes | invisible modulo a nonzero period |
| one root-marked transverse line with its `P=0` valuation | recovered |

The result is deliberately relative to the canonical quadratic-gauge
coordinates.  Promoting it to a map-intrinsic stable Torelli theorem would
require recognizing the root mark and the base character `P` from the
finite-normalization package.  The stable-moduli theorem recognizes `P` up
to scalar only after retaining the ordered boundary and Fitting data; that
is the next, stronger reconstruction layer.  Allowing arbitrary polynomial
gauge multipliers changes the finite-sampling answer completely: the
[finite marked-plane nonreconstruction theorem](FINITE_MARKED_PLANE_NONRECONSTRUCTION.md)
constructs infinitely many stable classes agreeing on any prescribed
finite set of marked inverse planes.

## 6. Exact regression

Run

```bash
.venv/bin/python scripts/verify_two_marked_fiber_gauge_reconstruction.py
```

The checker verifies the monic linear-coefficient formula in degrees four
through eleven, reconstructs every normalized seed coefficient at `P=1`,
recovers `m` and the stable area from `P=2`, checks the exact periodic
counterexamples on finite collections of cyclotomic planes, and computes
the transverse-line pole order `N+m`.  It also checks representative
squarefree fibers on both reconstruction planes.
