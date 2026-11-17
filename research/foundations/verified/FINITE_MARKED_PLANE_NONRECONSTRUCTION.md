# No finite collection of marked fibers reconstructs a stable class

The two-marked-fiber theorem recovers the exponent inside the monomial
power-shift family.  That positive result is locus-relative.  Once the same
quadratic gauge is allowed its natural polynomial freedom, no finite
collection of marked fibers—and not even finitely many complete marked
inverse planes—determines the stable class.

The obstruction is interpolation in the first target coordinate.  An
arbitrary polynomial multiplier can equal one on every sampled plane while
acquiring arbitrarily many new degree-drop boundary components elsewhere.

## 1. Polynomially interpolated quadratic gauges

Let `K` have characteristic zero, let `N>=4`, and fix

\[
 G(S)=\sum_{j=1}^Ng_jS^j,\qquad
 g_1g_3g_4\cdots g_N\ne0.
\tag{1.1}
\]

For a nonzero polynomial `R(P) in K[P]`, put

\[
\boxed{
 G_{P,R}(S)
 =g_1S+P(g_2S^2+g_3S^3)
  +R(P)\sum_{j=4}^Ng_jP^jS^j.
}
\tag{1.2}
\]

As usual, set

\[
 t=1+xy,\qquad
 q=t^2z+\frac{g_1}{g_3}y^2(1+3t),\qquad P=tq.
\tag{1.3}
\]

The denominator-free map `F_R=(P,B_R,C_R)` is

\[
\begin{aligned}
B_R={}&y+3\frac{g_3}{g_1}xq+2\frac{g_2}{g_1}tq\\
 &+R(P)\sum_{j=4}^N
      j\frac{g_j}{g_1}t^2x^{j-2}q^j,\\
C_R={}&x(5-3t)-\frac{g_3}{g_1}x^3z\\
 &-R(P)\sum_{j=4}^N
      (j-2)\frac{g_j}{g_1}x^jq^j.
\end{aligned}
\tag{1.4}
\]

These are polynomials for every polynomial `R`.  On `t!=0`, use

\[
 S=\frac xt,\qquad Q=y+xq,\qquad D=1-SQ+PS^2=\frac1t.
\tag{1.5}
\]

The inverse equation is

\[
\boxed{
 E_R(P,B,C;S)
 =G_{P,R}(S)-\frac{g_1}{2}(BS^2+C)=0.
}
\tag{1.6}
\]

The general marked-line gauge already allows arbitrary polynomial
coefficient functions in `P`.  At fixed `P`, its slope/intercept
calculation gives

\[
 \partial_SE_R=g_1D,\qquad
 \det\frac{\partial(B_R,C_R)}{\partial(S,Q)}=-2D.
\tag{1.7}
\]

Differentiating `R(P)` in the three-dimensional map only adds multiples of
`dP` to `dB_R,dC_R`, which disappear from
`dP\wedge dB_R\wedge dC_R`.  The reciprocal chart contributes `D^{-1}`.
Hence

\[
\boxed{\det DF_R=-2.}
\tag{1.8}
\]

The inverse equation is primitive and linear in `C`, has `S`-degree `N`
over `K(P,B,C)`, and its generic simple roots reconstruct.  Thus every
nonzero `R` gives a geometric-degree-`N` Keller map.

## 2. Interpolation on arbitrary finite samples

Fix any finite set of plane values

\[
 \mathcal C=\{c_1,\ldots,c_s\}\subset K.
\tag{2.1}
\]

Adjoin `1` to this set if it is not already present, and discard
repetitions.  Put

\[
 H(P)=\prod_{c\in\mathcal C}(P-c).
\tag{2.2}
\]

For every `k>=1`, consider

\[
 R_{k,\tau}(P)=1+\tau P^kH(P),\qquad \tau\in K^\times.
\tag{2.3}
\]

It satisfies

\[
 R_{k,\tau}(c)=1\quad(c\in\mathcal C),\qquad
 R_{k,\tau}(0)=1,
\tag{2.4}
\]

and has degree `k+s`.

There is a choice of `tau` for which `R_{k,tau}` is squarefree.  Indeed let
`f=P^kH`.  If `1+tau f` has a repeated root `alpha`, then

\[
 f'(\alpha)=0,\qquad
 \tau=-\frac1{f(\alpha)}.
\tag{2.5}
\]

Only finitely many values of `tau` arise from the finite critical locus of
the nonconstant characteristic-zero polynomial `f`.  Since `K` is
infinite, they can be avoided.

Choose one such `tau_k` for each `k` and write `R_k=R_{k,tau_k}`.  Formula
(1.4) gives

\[
\boxed{
 F_{R_k}|_{P=c}=F_1|_{P=c}
 \quad\text{for every }c\in\mathcal C.
}
\tag{2.6}
\]

This is coordinatewise equality of the restricted polynomial maps on the
source divisors, not merely an abstract isomorphism of inverse algebras.
Consequently every squarefree target fiber on any sampled plane is
literally common, with the same inverse-root coordinate and reconstruction.

## 3. Boundary created by a multiplier root

Let `rho` be a root of a chosen squarefree `R_k`.  It is nonzero by (2.4).
Use `u=P-rho` at the generic point of the target plane `P=rho`.

In (1.6), the coefficients through `S^3` are units and every coefficient
of `S^j`, `j>=4`, has valuation one.  The lower Newton polygon is

\[
\boxed{
 (0,0)\longrightarrow(3,0)\longrightarrow(N,1).
}
\tag{3.1}
\]

The horizontal block gives three finite simple roots.  Since `rho!=0`,
they reconstruct to affine source points.  The second block has horizontal
length `N-3`, height one, and therefore supplies one boundary prime with

\[
\boxed{(e,f)=(N-3,1).}
\tag{3.2}
\]

Its target image is the plane `P=rho`.  Distinct roots of `R_k` give
distinct target boundary components.

Because `R_k(0)=1`, the `P=0` Newton ledger is unchanged from the minimal
quadratic gauge:

\[
 (0,0)\longrightarrow(2,0)\longrightarrow(3,1)
 \longrightarrow(N,N).
\tag{3.3}
\]

It retains the canonical target image `P=0`.  Away from `P=0`, `R_k=0`,
and the repeated-root locus, every inverse root is simple and reconstructs
regularly.  The remaining boundary image is the irreducible ramified
discriminant, whose generic label is `(2,1)`.  This exhausts the
degree-`N` normalization.

Hence the complete geometric list of target images of boundary divisors
has cardinality

\[
\boxed{
 1_{\rm discriminant}+1_{P=0}+\deg R_k
 =k+s+2.
}
\tag{3.4}
\]

Equivalently, after deleting all boundary intersections, the normalized
ramified stratum has coordinate ring

\[
 \overline K[P^{\pm1},R_k(P)^{-1},r^{\pm1}],
\tag{3.5}
\]

whose unit rank modulo constants is `deg R_k+2`.

The canonical finite-normalization boundary is functorial under polynomial
left--right equivalence.  Identity stabilization takes every boundary
component times affine space and introduces no new units.  Therefore both
the component count (3.4) and the unit rank from (3.5) are stable
invariants.  Since they grow strictly with `k`, the maps `F_{R_k}` are
pairwise stably inequivalent.

## 4. Finite marked-fiber nonreconstruction theorem

### Theorem 4.1

Let `K` have characteristic zero and `N>=4`.  Given any finite set of
target planes `P=c_1,...,c_s`, there are infinitely many pairwise stably
inequivalent geometric-degree-`N` Keller maps that agree coordinatewise on
every corresponding source divisor.  They therefore share the entire
squarefree root-marked inverse cover over every sampled plane.

In particular, given any finite collection of squarefree marked target
fibers of one such map, infinitely many stable classes realize exactly the
same marked fibers at the same targets.

Thus there is no universal finite integer

\[
 r(N)<\infty
\tag{4.1}
\]

such that `r(N)` marked fibers recover the stable class of an arbitrary
degree-`N` Keller map, even inside the polynomial quadratic-gauge family.
The same negative statement holds if each sampled fiber is enlarged to its
complete two-dimensional inverse plane.

This does not contradict the two-marked-fiber theorem.  That theorem fixes
the rigid monomial locus `R(P)=P^m`; evaluating at a non-torsion `c`
determines its one exponent.  Polynomial interpolation leaves the sampled
values fixed while changing `R` away from them.

## 5. What a transverse line sees

The negative result is genuinely finite-sampling.  A root-marked inverse
cover on a transverse affine line retains a polynomial function of `P`, not
just finitely many values.

Normalize by `R(1)=1`, so that a marked fiber at `P=1` recovers `G/g_N`.
For a transverse line with nonconstant affine coordinate `p(u)`, the monic
annihilator has linear coefficient

\[
 \lambda_R(u)
 =\frac{g_1}{g_Np(u)^NR(p(u))}.
\tag{5.1}
\]

Therefore

\[
\boxed{
 R(p(u))
 =\frac{g_1}{g_Np(u)^N\lambda_R(u)}.
}
\tag{5.2}
\]

Since `p(u)` is an affine coordinate, one root-framed transverse family
recovers the complete polynomial multiplier `R` in this enlarged gauge
family.  This remains a coordinate-relative statement: it retains the
inverse-root generator and the base character `P`.

## 6. Answer to the reconstruction hierarchy

For degree at least four, the current sharp picture is:

| retained data | arbitrary polynomial quadratic gauges |
|---|---|
| any finite collection of unmarked fibers | insufficient |
| any finite collection of root-marked fibers | insufficient |
| any finite collection of complete parallel inverse planes | insufficient |
| `P=1` plus one non-torsion marked fiber, restricted to `R=P^m` | recovers `m` |
| one root-framed transverse inverse-cover line with known `P` | recovers `R` |
| full finite normalization plus boundary | detects at least `deg R` and all multiplier-root boundary planes |

Finite marked sampling is therefore impossible, but the final decorated
layer is complete.  The
[polynomial-gauge decorated Torelli theorem](POLYNOMIAL_GAUGE_DECORATED_TORELLI.md)
shows that the intrinsic boundary recovers the `P`-character and the
unmarked ramified Fitting divisor recovers the seed and multiplier up to
exactly the ordinary source--target scalings.  More generally, the full
finite-normalization morphism plus reconstruction boundary determines the
map by restriction to its distinguished affine open.

## 7. Exact regression

Run

```bash
.venv/bin/python scripts/verify_finite_marked_plane_nonreconstruction.py
```

The checker verifies the arbitrary-multiplier marked-line Jacobian, a
direct denominator-free quartic determinant, inverse, and reconstruction
identity, exact interpolation on three sample planes, squarefree
interpolants of growing degree, the multiplier-root Newton polygon and
ramification index, the unchanged `P=0` ledger, and the strictly increasing
boundary count `deg(R)+2`.
