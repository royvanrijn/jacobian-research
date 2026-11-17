# Three-real-variable weighted Gaussian search (archived)

> **Archive notice.**  This note predates the five-term explicit
> three-variable witness, which settles GMC negatively in every real Gaussian
> dimension \(n\geq3\) after adjoining unused coordinates.  Its finite
> weighted-family exclusions are retained for provenance only; the surviving
> ansatzes below are not an active research frontier.  See the
> [current GMC(2) program](../../extended-geometry/GMC2_RESEARCH_PROGRAM.md).

## 1. Outcome and scope

Let `h in C[z]`, `h(0)=1`, put

\[
 D=h-zh',
 \qquad g=u h(g),
\]

and retain the weighted bridge's observable `Q=Z`.  The target in this note is
the strong, seed-preserving version of a three-real-variable bridge:

\[
 \mathbb E(e^{uP})=1,
 \qquad
 \mathbb E(Ze^{uP})=g(u).                           \tag{1.1}
\]

Here `Z=(X+iY)/sqrt(2)`, `W=(X-iY)/sqrt(2)`, and `T` is a third independent
standard real Gaussian.  Merely assigning the same known three-variable
witness to every seed is not called a representation of the weighted family;
(1.1) preserves its injective mixed-moment fingerprint.

The present result is a **real obstruction, but not a complete
three-variable nonexistence theorem**:

1. every correction separated from `W` is classified by one exact scalar
   Gaussian identity;
2. every such correction of degree at most two in `T` is impossible for every
   nonlinear `h`, even with linear and constant terms allowed;
3. the Lagrange--Good determinant architecture itself cannot acquire a second
   fixed-point coordinate from a single unpaired real Gaussian without extra
   self-contractions;
4. exact Groebner calculations exclude several small `W`-dependent supports;
5. higher `T`-degree searches survive nontrivial finite truncations, so the
   unrestricted polynomial question remains open.

Thus the old failure was not peculiar to the displayed numerator of a few
weighted seeds: it extends to the entire separated quadratic class and every
nonlinear polynomial `h`.  It still does not justify an obstruction to an
arbitrary polynomial in `(Z,W,T)`.

## 2. Clean-room audit of the Gaussian--Lagrange lemma

The standalone
[formal Gaussian--Lagrange lemma](../../extended-geometry/FORMAL_GAUSSIAN_LAGRANGE_LEMMA.md) was
treated as read-only.  A clean-room pass reconstructed its proof in the
following order.

- The inverse of `F_u(z)=z-u Phi(z)` exists in the `(u,z)`-adic complete local
  ring because `F_u` is the identity modulo `u`; evaluating its inverse at
  the moving section `y=0` is therefore legal even when `Phi(0) != 0`.
- The coefficientwise Laurent module used in the note is large enough for
  `(z_i-u Phi_i)^{-b_i-1}`: at fixed `u`-degree it has a finite lower Laurent
  bound.  The homotopy `z-su Phi(z)` preserves this module coefficientwise.
- Differentiating the pullback along that homotopy gives a relative exact top
  form.  Its formal residue is zero, so the moving-section residue change is
  valid; no illicit substitution into an ordinary `z`-adic Laurent ring is
  being made.
- Expanding the changed denominators gives the stated multivariate binomial
  coefficient, and choosing `B=A/det(I-u JPhi)` cancels the determinant with
  the displayed, not reciprocal, orientation.
- Finally `E(Z^beta W^alpha)=delta_(alpha,beta) alpha!` converts the
  coefficient identity into the Wick identity.  Each `u` coefficient uses
  finitely many multi-indices.

No defect was found.  The existing nonlinear two-variable exact regression
was rerun successfully through `u^6`, as was the independent affine
three-real-variable Wick regression.  This is an internal clean-room check,
not the missing external specialist review, and the lemma file was not
modified.

## 3. Exact reduction for every separated correction

Consider the broad separated ansatz

\[
 P=W h(Z)+V(Z,T),                                   \tag{3.1}
\]

where `V` is any polynomial.  Define

\[
 F(u,z)=\mathbb E_T\!\left(e^{uV(z,T)}\right)
\]

coefficientwise.  Apply the one-variable Gaussian--Lagrange identity to each
`u` coefficient of the `u`-dependent prefactor `F` and sum.  This gives

\[
 \mathbb E(e^{uP})
 =\frac{F(u,g(u))}{1-u h'(g(u))}.                   \tag{3.2}
\]

Since `u=g/h(g)` and `1-u h'(g)=D(g)/h(g)`, (3.2) equals one if and only if

\[
 \boxed{
 \mathbb E_T\!\left[
  \exp\!\left(\frac{z}{h(z)}V(z,T)\right)
 \right]=\frac{D(z)}{h(z)}.}                       \tag{3.3}
\]

This is necessary and sufficient for all pure moments to vanish.  Moreover,
once (3.3) holds, insertion of `Z` in the same calculation multiplies the
right side of (3.2) by `g`; hence the second identity in (1.1) follows
automatically.  Equation (3.3), not the special half-pair formula, is the
classification of all `W`-separated one-complex-plus-one-real corrections.

It also gives a precise higher-chaos search problem.  Writing

\[
 V(z,T)=\sum_{j=0}^d v_j(z)\operatorname{He}_j(T)
\]

turns (3.3) into a prescribed Gaussian exponential integral.  For `d>=3`
this is a finite-rank integration-by-parts connection, but it is no longer an
algebraic square-root identity.  Clearing its connection denominators can
create spurious components, so chart saturation is essential in any
Groebner proof.

## 4. Complete obstruction through quadratic real chaos

Take the most general separated polynomial of `T`-degree at most two,

\[
 V(z,T)=c(z)+a(z)T+b(z)T^2.                         \tag{4.1}
\]

At `u=z/h`, elementary Gaussian integration turns (3.3) into

\[
 \frac{
  \exp\!\left(uc+\dfrac{u^2a^2}{2(1-2ub)}\right)}
 {\sqrt{1-2ub}}
 =\frac{D}{h}.                                      \tag{4.2}
\]

The exponential of a nonconstant rational function cannot be algebraic over
`C(z)`: a pole, including the pole at infinity, would become an essential
singularity.  The exponential factor in (4.2) must therefore be one (its
value at `z=0` is one).  Consequently

\[
 \boxed{
 b=-\frac{h h'(2h-zh')}{2D^2},
 \qquad
 c=-\frac{z a^2D^2}{2h^3}.}                         \tag{4.3}
\]

The first formula is independent of the optional linear term `a`.  Thus a
linear Gaussian correction cannot repair the half-pair denominator.

There is in fact no nonlinear exceptional polynomial `h`.  Suppose
`deg(h)=d>=2`.  Then `deg(D)=d` and `D(0)=1`.  If every root of `D` were a
root of `h`, it would be a multiple root of `h`, because `D(alpha)=h(alpha)=0`
implies `h'(alpha)=0`.  At a root of `h` of multiplicity `m`, the multiplicity
in `D=h-zh'` is exactly `m-1` (all roots are nonzero because `h(0)=1`).  Hence
the total multiplicity of all roots of `D` would be at most

\[
 \sum_{h(\alpha)=0}(m_\alpha-1)<d,
\]

contradicting `deg(D)=d`.  Therefore `D` has a root `beta` with
`h(beta) != 0`.  At this root,

\[
 h'(\beta)=h(\beta)/\beta\ne0,
 \qquad 2h(\beta)-\beta h'(\beta)=h(\beta)\ne0.
\]

So the numerator of `b` in (4.3) is nonzero at `beta`; even `D`, let alone
`D^2`, does not divide it.

### Theorem 4.1

For every nonlinear polynomial `h` with `h(0)=1`, there are no polynomials
`a,b,c` for which (4.1) produces (1.1).  For affine `h=1+az`, `D=1` and the
known half-pair construction is recovered.  This completely classifies the
separated `T`-quadratic case.

## 5. Why a nontriangular fixed point does not come for free

The determinant formula used by the four-real-variable bridge is based on a
Gaussian polarization

\[
 \mathbb E(Z_iZ_j)=\mathbb E(W_iW_j)=0,
 \qquad \mathbb E(Z_iW_j)=\delta_{ij}.              \tag{5.1}
\]

Each fixed-point coordinate has a distinct Wick-dual source coordinate.  The
span of the variables in (5.1) is hyperbolic of even dimension `2r`; a
Lagrange--Good system with no self-contracted leftover directions is therefore
even-dimensional.  Three real Gaussians complexify as one such hyperbolic pair
plus an anisotropic line `T` with `E(T^2)=1`.

If one writes a putative second source term `T Psi(Z,T)`, contractions of `T`
occur both with `Psi` and with other source copies of `T`.  These are the
self-contractions absent from the proof of the Gaussian--Lagrange lemma, and
they generate additional terms rather than a second row of
`det(I-u JPhi)`.  Making `Phi` nontriangular does not remove them.

Therefore an exact two-component Good fixed-point system cannot simply be
installed on one complex plus one real Gaussian.  A successful nontriangular
three-real construction would need a genuinely different identity that
controls the anisotropic self-contractions.  This is an obstruction to the
current determinant-cancelling architecture, not to arbitrary polynomials.

## 6. Corrections involving both `Z` and `W`

Once the correction contains `W`, (3.2) no longer applies: powers of the
correction change the Wick-source degree as well as the target polynomial.
There is no divisor theorem presently covering this class.  Exact finite
support searches are nevertheless decisive when their moment ideal is the
unit ideal.

The script
[`search_three_real_weighted_ansatz.py`](../tooling/search_three_real_weighted_ansatz.py)
uses the first nonlinear canonical seed

\[
 h=1+z^2-z^3
\]

and forms the exact rational moment ideal for

\[
 P=Wh(Z)+\sum c_{kij}T^kW^iZ^j.                    \tag{6.1}
\]

Including both the pure moments and the target mixed moments, its rational
Groebner basis is the unit ideal in the following cases:

| real powers | `W` degree | `Z` degree | moments used |
|---|---:|---:|---:|
| `T^2` | 1 | 0 | `m<=5` |
| `T^2` | 1 | 1 | `m<=5` |
| `T^2` | 1 | 2 | `m<=4` |

These are exact exclusions of the displayed finite supports: an all-order
witness would in particular solve those finite equations.  They are not an
exclusion of arbitrary `W` degree or arbitrary support.

As a positive control, selecting the affine seed, `W` degree zero, `Z` degree
two, and `T^2` recovers the unique coefficients
`(-1,-3/2,-1/2)` of Long's correction through the same moment equations.

For comparison, the `W`-free support with `T,T^2,T^3`, coefficient degree at
most one in `Z`, has a nonunit ideal through `m=6`.  Thus higher real chaos has
enough freedom to survive meaningful truncations.  Survival is not an
all-order construction; it explains why extrapolating the quadratic divisor
argument would be unsafe.

The separate cubic-Hermite connection search
[`explore_three_gaussian_cubic_chaos.py`](../tooling/explore_three_gaussian_cubic_chaos.py)
is useful on the more normalized ansatz without an independent zeroth
Hermite term.  Over `Q`, all genuine-cubic charts with coefficient degree at
most one give a unit compatibility ideal for `h=1+z^2-z^3`.  Degree two is
currently only modular/low-coefficient reconnaissance; its full
characteristic-zero elimination was not completed.

## 7. Reproduction and open frontier

Run:

```bash
.venv/bin/python scripts/verify_formal_gaussian_lagrange.py
python3 scripts/verify_long_gaussian_moments.py
.venv/bin/python archive/tooling/search_three_real_weighted_ansatz.py \
  --z-degree 1 --w-degree 1 --moment-order 5
.venv/bin/python archive/tooling/search_three_real_weighted_ansatz.py \
  --seed affine --z-degree 2 --w-degree 0 --moment-order 5 --print-basis
.venv/bin/python archive/tooling/search_three_real_weighted_ansatz.py \
  --z-degree 2 --w-degree 1 --moment-order 4
.venv/bin/python archive/tooling/search_three_real_weighted_ansatz.py \
  --z-degree 1 --w-degree 0 --t-degrees 1,2,3 --moment-order 6
.venv/bin/python archive/tooling/explore_three_gaussian_cubic_chaos.py \
  --degree 1 --chart 0 --prime 0
.venv/bin/python archive/tooling/explore_three_gaussian_cubic_chaos.py \
  --degree 1 --chart 1 --prime 0
```

The most useful next target is not another unstructured moment expansion.  It
is an all-orders analysis of (3.3) for `T`-degree three including an
independent `He_0` coefficient, followed by a saturated characteristic-zero
elimination.  In parallel, `W`-dependent corrections need a normal form under
Gaussian linear changes before support-by-support searches can become a
classification.

At the time of this search, its internal answer was:

\[
 \boxed{\text{no separated quadratic bridge for any nonlinear seed; the
 unrestricted weighted-family ansatz was not closed by these calculations.}}
\]

That historical search gap is now moot for GMC witness existence in dimension
three: the five-term explicit witness settles it without using this weighted
family.
