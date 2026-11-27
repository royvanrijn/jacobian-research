# An unrestricted three-cell comparison through correlation 317/500

Every ordinary balanced Gaussian three-cell partition in dimension
at least two has noise stability at most that of the simplex for
$0<\rho\le317/500$. On the added interval $[3/5,317/500]$, equality
requires a simplex up to rotation, relabeling and null sets. The
proof handles all planar circular profiles and all angular shapes.
Existence and dimension reduction for maximizers apply as in the
[comparison through $3/5$](THREE_CELL_QUARTIC_EXTENSION_THROUGH_THREE_FIFTHS.md).

The [directed checker](check_three_cell_excess_purity_extension.py)
covers the closed interval from $3/5$. The comparison uses internal
verification only.

## Combining unconditional variance and conditional purity

Use $d=1-\rho^2$, $\kappa=\rho^2/d$ and $d\nu(t)=e^{-t}dt$. Let
$\mu_1$ have density $\sqrt t\,e^{-\kappa t}/Z_1$ relative to $\nu$,
with $Z_1=\Gamma(3/2)d^{3/2}$. Write $z=p-\mathbf1/3$, let
$P_1=\mathbb E_{\mu_1}(1-|p|^2)$ be the purity factor and let
$V_1=\operatorname{Var}_{\mu_1}(p)$. Set $\omega_1=\rho d/(2\pi)$.

For the sharp aligned-arc scalar form $q$ of (GH4), its total
Taylor slack has both lower bounds

\[
               \frac{12}{25}\omega_1 V_1,
               \qquad \omega_1P_1V_1.
\]

The first is (SV1) and the second follows from (QE1)–(QE2).
Taking the larger of these bounds gives the comparison

\[
 S_\rho(A)-S_\rho(T)
 \le q^\#(z)-\omega_1(P_1-12/25)_+ V_1,
 \qquad q^\#(z)=q(z)-\frac{12}{25}\omega_1V_1.             \tag{EP1}
\]

Thus the conditional correction charges only the excess over the
unconditional variance credit. [Unconditional variance bound](THREE_CELL_SHARP_VARIANCE_ABSORPTION.md)
[Quartic purity bound](THREE_CELL_QUARTIC_EXTENSION_THROUGH_THREE_FIFTHS.md)

In the positive scalar series for $q^\#$, the only change is to
replace $M_1$ by $M_1+24/(25\pi^2)$ and keep $P_1$ in the angular
series equal to the original $P_{\sharp,1}$. Indeed, in (SV2) the
added row and kernel coefficients are $2d\beta_1/Z_1$ and
$\beta_1/Z_1^2$, where $\beta_1=(12/25)\omega_1$. Substitution in
(GH3) and (GH5) gives the stated common increment to the degree-one
moment. In the remainder of this note the purity factor is denoted
by $\mathcal P=\mathbb E_{\mu_1}(1-|p|^2)$ to distinguish it from
the angular series moments.

## Fixed slab operator and trial data

On each rational slab $[l,u]$ use the anchored comparison
(QE3)–(QE4) for this amended series. In particular the row at $u$
is mixed by $(l/u)^n$ in degree $n$, the finite kernel is at $u$,
and the balance leakage is bounded against the anchored mean-zero
projection. The unamended full row $J$ bounds the tail by (SV5).
The additional row term of (SV2) is positive, so the same tail
bound applies to the amended form. The checker uses $m=54$ row
terms, $N=168$ kernel degrees and the enclosed linear-interpolation
moments of the variance checker on $[0,14]$ with 4096 panels.

For clarity, in the notation of (QE4) the trial coefficients in the
centered features $g_i$ are $\alpha_i=\sqrt{k_i}\,r_i$. A floating-
point diagonalization of the midpoint Gram matrix chooses the vector
$r$, and the checker rounds its coordinates to multiples of
$2^{-46}$. For these exact dyadic coefficients it encloses the trial
norm, Rayleigh quotient, complement Hilbert–Schmidt norm, residual,
and mean-zero leakage. The resulting quadratic bound is

\[
                    q^\#(z)\le eA^2+2RAG-\delta G^2.         \tag{EP2}
\]

The physical reference and its pole scale as in (QE4)–(QE5).
On the anchor core $0\le\sqrt{t_0}\le1/4$ the checker bounds
$\sqrt{t_0}\phi(t_0)$ below by $c-C\sqrt{t_0}$ with $C=13/10$.
The positive first row term gives the lower weight $j x$ used in
the origin clip over $0\le x\le1/4$.

For the physical $\mu_1$ density relative to the mixed weighted
measure one can use

\[
 B=\frac{8\pi}{l(1-u^2)(11/2+24/25)}.
\]

This is the density formula in the quartic checker with the
amended degree-one coefficient. Remainder orthogonality and
(SR15)–(SR19) apply in the mixed weighted space. The checker
bounds the reference variance below by $V$ from the core
$0\le\sqrt{t_0}\le4$ using the two-sample variance formula and
the scale factor $(1-u^2)/(1-l^2)$. It bounds the raw second
moment above by $M_\phi$ on the anchor space. The latter is also
an upper bound on the physical raw second moment because the
scale factor is at most one. [Variance and moment inequalities](THREE_CELL_SHARP_SCALAR_REMAINDERS.md#rank-one-bounds-for-the-quartic-penalty)

For the raw second moment the checker integrates on $[0,12]$.
Writing $F(x)=e^{-\kappa_u x^2}\sum_i\alpha_i x^i$ and $c_0$ for
its centering constant, its integrand is proportional to
$e^{-(1-\kappa_u)x^2}(F(x)-c_0)^2/Q(x)^2$ with the positive row
polynomial $Q$. For $x\ge L=12$, use $Q(x)\ge Q(L)$ and

\[
 |F(x)|\le e^{-\kappa_u x^2}(x/L)^N\sum_i|\alpha_i|L^i,
 \qquad (F(x)-c_0)^2\le2F(x)^2+2c_0^2.
\]

The two Gaussian tails then have the elementary bound used in
(GH8). This verifies an upper raw moment over the whole half-line.

## Directed amplitude test

If the scalar upper bound in (EP2) can be nonnegative and $A>0$,
put $r=G/A$. For an upper bound $r_+$ on $r$, use the pole
cap $A\le A_{r_+}$ from (QE5). If $M_\phi<B/(1+r_+^2)$, use
the upper raw coefficient $(\sqrt{M_\phi}+r_+\sqrt{B-M_\phi})^2$
from (SR19). Otherwise use $B(1+r_+^2)$. Calling this coefficient
$D_+$, the purity and variance bounds are

\[
 \mathcal P\ge 2/3-A_{r_+}^2D_+,
 \qquad
 V_1\ge A^2(\sqrt V-r_+\sqrt{B-V})_+^2 .                  \tag{EP3}
\]

On each of 128 equal intervals $[r_-,r_+]$ from zero to the
positive scalar root, the checker computes the maximum of
$e+2Rr-\delta r^2$. It either verifies a negative scalar upper
bound or compares it to the extra term in (EP1), using (EP3)
and $u(1-u^2)/(2\pi)\le\omega_1$. This proves strictness for
every nonuniform mass profile. If $A=0$, (EP2) is strictly
negative unless the mass profile is uniform. The equal-circle
rigidity argument then gives the simplex as the only tie.

The checker uses 120-bit Arb enclosures and 136 adjacent rational
correlation slabs of width $1/4000$. The replay command is

    OPENBLAS_NUM_THREADS=1 .venv/bin/python research/plurality-stability/check_three_cell_excess_purity_extension.py

## Discrete plurality

For every fixed $0<\rho\le317/500$ and $\varepsilon>0$ there
exist $D$ and $\tau>0$ such that for every $m\ge1$, any balanced
$f:[3]^m\to\Delta_3$ with $\operatorname{Inf}_i^{\le D}(f_j)\le\tau$
for all $i,j$ satisfies

\[
 S_\rho(f)\le\lim_{k\to\infty}S_\rho(\operatorname{PLUR}_{k,3})+\varepsilon.
\]

As in the preceding transfer, the noise independently retains each
coordinate with probability $\rho$ and otherwise resamples it
uniformly. The Gaussian realization preserves balance and the
discrete plurality limit is the simplex stability.
[Isaksson–Mossel, Theorem 7.1 and Lemma 7.3](https://arxiv.org/pdf/0903.3362)
