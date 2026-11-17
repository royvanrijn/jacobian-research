# Complete integer parameters on the fixed Keller target line

This note counts every signed integer parameter giving a Hasse failure on
the target line of the [fixed-map manuscript](../../unpublished/fixed-map-hasse-failures/sections/03-local-global.tex).
The earlier sufficient positive integer family has a smaller leading
constant. The rational-line count and its ambient interpretation remain
separate results.

## 1. Integer height and the exact multiplicative variable

Retain the fixed map `Phi` and targets

\[
 y_a=(-1,32a/9,(8a+1)/3),\qquad a\in\mathbb Z.
\]

For `a != 0`, the [rational local classification](../../unpublished/fixed-map-hasse-failures/sections/03-local-global.tex)
says that the fibre is everywhere locally soluble and has no rational
point precisely when `a` is not a cube and

\[
 3\mid v_p(a)\quad(p\equiv2\pmod3),\qquad
 3\mid v_3(a),\qquad
 a/3^{v_3(a)}\equiv\pm1\pmod9.                    \tag{1.1}
\]

The target at `a=0` has no real point, so it contributes nothing.
For every nonzero integer `a`, the primitive-coordinate calculation is

\[
 [1:y_a]=[9:-9:32a:24a+3],\qquad
 H(y_a)=\frac{32|a|}{\gcd(a,3)}.                    \tag{1.2}
\]

Indeed the content is `gcd(a,3)`, and `32|a|` is the largest absolute
coordinate. In particular height and condition (1.1) are unchanged by
`a -> -a`; the integer count has exact sign symmetry.

Define the positive height index

\[
 h=H(y_a)/32=|a|/\gcd(a,3).
\]

For parameters satisfying (1.1), its exponent at three is

\[
 v_3(h)\in\{0,2,5,8,\ldots\}.                       \tag{1.3}
\]

Conversely, from any positive integer `h` with (1.3), recover the positive
parameter by `a=h` if `3` does not divide `h`, and by `a=3h` otherwise.
The inert-prime conditions pass unchanged between `a` and `h`. If
`h_0=h/3^{v_3(h)}`, the last condition in (1.1) becomes
`h_0 == +/-1 (mod 9)`. This gives a bijection before excluding cubes.

## 2. The complete integer count theorem

Put

\[
 I_\Phi(B)=\#\{a\in\mathbb Z:H(y_a)\le B,
                 \ \Phi^{-1}(y_a)\text{ fails the Hasse principle}\}.
\]

Let `chi_{-3}` be the quadratic Dirichlet character modulo three, and set

\[
\begin{aligned}
 G_3(s)&=\left[L(s,\chi_{-3})(1-3^{-s})
                   \prod_{p\equiv2(3)}(1-p^{-2s})\right]^{1/2},\\
 J(s)&=\prod_{p\equiv2(3)}(1-p^{-3s})^{-1},\\
 A_3(s)&=1+\frac{3^{-2s}}{1-3^{-3s}},\\
 Q(s)&=G_3(s)J(s)A_3(s).                            \tag{2.1}
\end{aligned}
\]

The square root is the branch positive on the real line immediately to the
right of one. These functions are holomorphic near one: the inert Euler
products converge absolutely there, and `L(1,chi_{-3})` is nonzero.
All logarithms below are natural.

**Theorem 2.1 (complete signed integer count).** As `B -> infinity`,

\[
\boxed{
 I_\Phi(B)=C_{\rm int}\frac{B}{\sqrt{\log B}}
       \left(1+\frac{\kappa_{\rm int}}{\log B}\right)
       +O\!\left(\frac{B}{(\log B)^{5/2}}\right),}
                                                               \tag{2.2}
\]

where

\[
\boxed{
 C_{\rm int}=\frac{29G_3(1)}{1248\sqrt\pi}
                 \prod_{p\equiv2(3)}(1-p^{-3})^{-1},\qquad
 \kappa_{\rm int}=\frac{\log32}{2}+\frac12-\frac\gamma4
                         -\frac{Q'(1)}{2Q(1)}.}       \tag{2.3}
\]

Here `gamma` is Euler's constant. The positive integer count is exactly
`I_Phi(B)/2`, including at finite height. In particular the leading-term
version has error `O(B/(log B)^{3/2})`.

### Proof of Theorem 2.1: character projector and analytic hypotheses

Let `f(h)` be the multiplicative indicator imposing (1.3) and the inert
valuation conditions, without the residue condition. Its Dirichlet series
is

\[
 F(s)=A_3(s)\prod_{p\equiv1(3)}(1-p^{-s})^{-1}
                 \prod_{p\equiv2(3)}(1-p^{-3s})^{-1}
     =\zeta(s)^{1/2}Q(s).                            \tag{2.4}
\]

Choose a character `chi` of order three on the units modulo nine. Extend
it to a completely multiplicative function `chi_*` on positive integers
by defining `chi_*(3)=1`, and using `chi_*(p)=chi(p)` at all other primes.
Thus `chi_*(h)=chi(h_0)`. This extension is essential: the ordinary
Dirichlet character would vanish on the valid height indices divisible by
three. Its kernel on the units is `{1,-1}`, so

\[
 \mathbf1_{h_0\equiv\pm1(9)}
     =\frac13\sum_{j=0}^2\chi_*(h)^j.                \tag{2.5}
\]

Apply de la Breteche--Tenenbaum,
[Theorem 1.2, formulas (1.10)--(1.12), author PDF p. 3](https://tenenb.perso.math.cnrs.fr/PPP/On-SD.pdf),
to `f(h)chi_*(h)^j`. Every prime-power value has absolute value at most one.
The theorem's class `F(r,sigma)` applies with `r=1/2`, `sigma=3/4`:
the absolute prime values are the indicator of the primes `1 (mod 3)`,
and both required higher-power series converge. Fixed-modulus prime
distribution modulo nine gives the prime-mean error with every fixed
logarithmic power. The prime means are `rho=1/2` for `j=0` and `rho=0`
for `j=1,2`, because either nontrivial character sums to zero on
`{1,4,7}`. The altered factor at three affects none of these prime means.

This is the same inspected analytic theorem and fixed-modulus prime input
as in the manuscript's [source audit](../../unpublished/fixed-map-hasse-failures/LITERATURE_AUDIT.md).
For `rho=0`, all asymptotic coefficients vanish by formula (1.10).
Taking its prime-error exponent `A=3` gives errors smaller than
`X/(log X)^{5/2}`, after absorbing the theorem's fixed iterated-log factor.
For the principal function, the coefficient of order two can itself be
absorbed into `O(X/(log X)^{5/2})`.

To compute the first correction with its sign, write `z=s-1`. The
coefficient-generating expression in (1.10) is

\[
 \frac{z^{1/2}F(1+z)}{1+z}
 =Q(1)+\left[Q'(1)+(\gamma/2-1)Q(1)\right]z+O(z^2).
                                                               \tag{2.6}
\]

The constant term is divided by `Gamma(1/2)`, and the linear coefficient
by `Gamma(-1/2)=-2Gamma(1/2)`. Consequently

\[
 \sum_{h\le X} f(h)
 =\frac{Q(1)X}{\sqrt\pi\sqrt{\log X}}
  \left[1+\frac{1/2-\gamma/4-Q'(1)/(2Q(1))}{\log X}\right]
  +O\!\left(\frac{X}{(\log X)^{5/2}}\right).          \tag{2.7}
\]

Equations (2.5) and (2.7) supply a factor `1/3`, and the two signs supply
a factor two. Substitute `X=B/32`. The expansion of
`(log(B/32))^{-1/2}` adds `(log32)/2` to the correction. Finally

\[
 A_3(1)=1+\frac{1/9}{1-1/27}=\frac{29}{26},\qquad
 \frac{2}{3\cdot32\sqrt\pi}Q(1)
       =\frac{29G_3(1)J(1)}{1248\sqrt\pi}.
\]

The calculation so far includes global cubes. If `a=b^3`, (1.2) and
`H(y_a)<=B` imply `|b|^3<=3B/32`; hence only `O(B^{1/3})` such parameters
occur. Removing them changes neither displayed coefficient nor the stated
error, and the local classification then identifies exactly the Hasse
failures. This proves the theorem.

## 3. Comparison with the sufficient positive family

The [earlier count](../../unpublished/fixed-map-hasse-failures/sections/04-counting.tex)
uses positive parameters supported entirely on split primes, with residue
one modulo nine. Its leading target-height constant is
`C_suff=G3(1)/(96sqrt(pi))`. Therefore

\[
 \frac{C_{\rm int}}{C_{\rm suff}}
     =\frac{29}{13}\prod_{p\equiv2(3)}(1-p^{-3})^{-1}.
                                                               \tag{3.1}
\]

For positive parameters alone the ratio is half of (3.1). The extra
factors account for inert prime cubes, powers of three together with their
target content, and, in (3.1), the second sign.

## 4. Verification scope

Run `python3 research/arithmetic/verify_fixed_map_integer_count.py` from
the repository root. The [checker](verify_fixed_map_integer_count.py) uses
only Python's standard library. It compares the height-index bijection
with the original integer local criterion, recomputes primitive heights,
checks the finite character projector and local Euler-factor algebra, and
compares direct and indexed bounded-height counts. Its
[receipt](fixed_map_integer_count_verification.json) records the exact
bounds, checker and proof hashes, and outcomes. These are finite internal
controls, not verification of an asymptotic by finite extrapolation.

The analytic conclusion uses the written proof and the imported
Selberg--Delange theorem. Formal verification, independent whole-proof
replay, external mathematical review and publication are not asserted.
