# Leading constant for Hasse failures on the second fixed-map target line

This note sharpens the two-sided order in [Theorem FMA6](FIXED_MAP_AMBIENT_LOCAL_CLASSIFICATION.md) to an asymptotic with an explicit positive constant. All targets belong to the same fixed degree-five determinant-one Keller map `Phi` over `Q`. The analytic input is exactly the coprime-pair Selberg–Delange estimate already proved in [the rational-counting section](../../unpublished/fixed-map-hasse-failures/sections/04-rational-counting.tex), equation `eq:rational-polygon-count`; no additional analytic theorem is asserted. The proof and new height calculation are internal written mathematics, with exact controls in [the checker](verify_fixed_map_second_line_count.py). Formal verification and external review are not claimed.

Set

\[
 \gamma=\frac3{32},\qquad k=\gamma^3=\frac{27}{32768},\qquad
 z_r=\left(2,\frac{56}3-\frac{1024}9r,-\frac76-\frac{64}3r\right).
                                                        \tag{1}
\]

Let `H` be projective multiplicative height, so a primitive integer representative has height the maximum absolute coordinate. For `r in Q`, write `H(z_r)=H([1:z_r])`. Define

\[
 N_2(B)=\#\{r\in\mathbb Q:H(z_r)\le B,
                         \Phi^{-1}(z_r)\text{ is Hasse-failing}\}.
\]

Here Hasse failure means local points at every completion and no rational point. FMA4 gives the exact classification: such targets are precisely `r=a+k/a`, where

\[
 \mathcal H=\{a\in\mathbb Q^\times\setminus\mathbb Q^{\times3}:
 3\mid v_p(a)\ (p\equiv2\bmod3),\ a\in\mathbb Q_3^{\times3}\}.
                                                        \tag{2}
\]

Every target in this image has exactly two parameters, `a` and `k/a`, both in `H`. The parameter map is not injective and is not surjective onto the rational line.

## 1. Theorem SLC1: the second-line leading constant

Put `D=32768`, `K=27`, and define the positive homogeneous height function

\[
 Q(x,y)=\max\left\{
 9216|xy|,
 16|5376xy-Dx^2-Ky^2|,
 3|1792xy+Dx^2+Ky^2|
 \right\}.
                                                        \tag{3}
\]

Let

\[
 \Omega_2=\{(x,y)\in\mathbb R^2:y>0,\ Q(x,y)\le1\},\qquad
 \mathcal A_2=\operatorname{area}(\Omega_2).
\]

Then

\[
 \boxed{N_2(B)\sim C_2\frac{B}{\log B}},\qquad
 \boxed{C_2=\frac{256\mathcal A_2}{13\sqrt3\,\zeta(2)}
       \prod_{p\equiv2\ (3)}\frac{1+p^{-3}}{1-p^{-3}}}.
                                                        \tag{4}
\]

The product converges absolutely and `C_2>0`. An explicit one-dimensional integral for the archimedean factor is

\[
 \mathcal A_2=\frac12\int_{-\infty}^{\infty}\frac{dt}{Q(t,1)}.
                                                        \tag{5}
\]

For orientation, numerical evaluation gives `A_2` approximately `0.000099194120294697` and `C_2` approximately `0.0008977265476`; the proof uses the exact integral, not numerical quadrature.

## 2. Lemma SLC2: exact primitive height and its local content

Write `a=m/n`, with `m!=0,n>0,gcd(m,n)=1`, and put `A=Dm^2+Kn^2`. A projective integer representative for `z_(a+k/a)` is

\[
 [4608mn:9216mn:16(5376mn-A):-3(1792mn+A)].                \tag{6}
\]

Its positive common divisor is exactly

\[
 g(m,n)=2^{\min(2v_2(n),15)}3^{\min(2v_3(m),3)}.
                                                        \tag{7}
\]

Consequently

\[
 H(z_{a+k/a})=\frac{Q(m,n)}{g(m,n)}.                     \tag{8}
\]

**Proof.** Substituting `r=A/(Dmn)` in (1), clearing denominators by `18Dmn` and dividing all four coordinates by `128` gives (6). A common prime divisor other than two or three would divide `mn` through the first coordinate. If it divides `m`, the last coordinate reduces to `-81n^2`; if it divides `n`, the third reduces to `-524288m^2`. Coprimality excludes both cases.

For the valuation at two, there is no common factor if `n` is odd, because the last coordinate is then odd. If `e=v_2(n)>0`, `m` is odd. The three summands in `1792mn+A` have valuations `e+8,15,2e`. Their minimum is unique and equals `min(2e,15)`: the possible ties between `e+8` and another valuation are above the third valuation, and `15=2e` is impossible. The third coordinate has valuation `4+min(2e,15)` and the first has valuation `9+e`, so the last coordinate realizes the common minimum in (7).

For three, there is no common factor if `m` is a unit at three, because the third coordinate is then nonzero modulo three. If `f=v_3(m)>0`, `n` is a unit. The three summand valuations inside the third coordinate are `f+1,2f,3`; those inside the last coordinate, before its outside factor three, are `f,2f,3`. When `f=1` the last coordinate has valuation two and the third has at least two. When `f=2` the last has valuation three and the third has at least three. When `f>=3` the third has valuation exactly three and the last at least four. The first coordinate has valuation `2+f`. These cases prove the exponent `min(2f,3)`, including possible cancellations in the other coordinate. ∎

The second coordinate in (6) is twice the first; its contribution `9216|mn|` is part of the projective height and is included in (3).

## 3. The coprime-pair analytic input

Let `g_0` be the multiplicative indicator on positive integers defined by

\[
 g_0(u)=1\quad\Longleftrightarrow\quad
 (u,6)=1\text{ and }3\mid v_p(u)\text{ for every }p\equiv2\pmod3.
\]

The existing rational-counting proof supplies constants `K_0>0` and `E(2)>0` such that, for every fixed bounded Jordan measurable region `R` in the positive quadrant,

\[
 \#\{(u,v)\in T R\cap\mathbb Z^2:
 g_0(u)g_0(v)=1,\ (u,v)=1,\ u/v\equiv\pm1\pmod9\}
 \sim\frac{K_0^2E(2)}3\operatorname{area}(R)
                  \frac{T^2}{\log T}                    \tag{9}
\]

when its area is positive, with a corresponding little-oh bound when it is zero. Its exact constant is

\[
 K_0^2E(2)=\frac7{36\sqrt3\,\zeta(2)}
       \prod_{p\equiv2\ (3)}\frac{1+p^{-3}}{1-p^{-3}}.   \tag{10}
\]

It also supplies the uniform upper bound

\[
 \sum_{u\le X}g_0(u)\ll X/\sqrt{\log(2X)}.               \tag{11}
\]

The source proves (9) first for rectangles, with the cubic-character projector modulo nine and the absolutely convergent coprimality convolution, then by finite grid approximation. Thus its use for the curved but Jordan measurable regions below requires no new uniform asymptotic theorem.

## 4. Proof of Theorem SLC1

Every reduced parameter satisfying the local conditions in (2), including global cubes temporarily, has a unique expression

\[
 m=\varepsilon 8^i27^j u,\qquad n=8^k27^l v,
 \quad \varepsilon\in\{1,-1\},\quad i,j,k,l\ge0,
 \quad ik=jl=0,                                         \tag{12}
\]

where `g_0(u)g_0(v)=1`, `(u,v)=1`, and `u/v=+-1 mod9`. The powers of eight and the sign are cubes modulo nine; the powers of three have been removed before applying the unit condition.

Set `r_0=8^i27^j`, `s_0=8^k27^l`. Lemma SLC2 makes the content on this stratum constant:

\[
 g_{i,j,k,l}=2^{\min(6k,15)}27^{\mathbf1_{j>0}}.           \tag{13}
\]

For each fixed stratum and each sign, (8) becomes

\[
 (\varepsilon r_0u,s_0v)\in\sqrt{g_{i,j,k,l}B}\,\Omega_2.
\]

The corresponding positive-quadrant regions in `(u,v)` have combined area `g_(i,j,k,l) A_2/(r_0s_0)`. Using (9) with `T=sqrt(B)` therefore gives the fixed-stratum parameter count

\[
 \frac{2K_0^2E(2)}3\,
 \frac{g_{i,j,k,l}\mathcal A_2}{r_0s_0}
                  \frac{B}{\log B}(1+o(1)).             \tag{14}
\]

The factor two is `log B/log sqrt(B)`, before quotienting the two-to-one parametrization.

We justify summing all strata. If `F_2,F_3` are the last two entries of (6), then

\[
 F_2+16F_3=-64(Dm^2+Kn^2),\qquad
 Q(m,n)\ge\frac{64}{17}(Dm^2+Kn^2).                     \tag{15}
\]

Also `g<=DK` by (7). Hence a height-bounded pair satisfies `|m|,n<=C sqrt(B)` with an absolute constant, so `u<=C sqrt(B)/r_0` and `v<=C sqrt(B)/s_0`. For `r_0,s_0<=B^(1/8)`, (11), after dropping coprimality and the unit restriction, bounds its count by

\[
 \ll\frac{B}{r_0s_0\log B}.
\]

The sum of `1/(r_0s_0)` over the permitted indices converges. For the remaining strata, the trivial bound is `O(B/(r_0s_0))`, and

\[
 \sum_{\max(r_0,s_0)>B^{1/8}}\frac1{r_0s_0}
 \le B^{-1/16}\sum_{i,j,k,l\ge0}(8^{i+k}27^{j+l})^{-1/2}
 \ll B^{-1/16}.                                        \tag{16}
\]

Thus those strata contribute `o(B/log B)`. Dominated convergence applies to (14).

The valuation weight factors exactly:

\[
 \begin{aligned}
 W&=\sum_{ik=jl=0}\frac{g_{i,j,k,l}}{8^{i+k}27^{j+l}}\\
 &=\left(1+\sum_{i\ge1}8^{-i}
          +\sum_{k\ge1}2^{\min(6k,15)}8^{-k}\right)
   \left(1+27\sum_{j\ge1}27^{-j}+\sum_{l\ge1}27^{-l}\right)\\
 &=\left(1+\frac17+8+64+\frac{512}7\right)
       \left(1+\frac{27}{26}+\frac1{26}\right)
   =\frac{1024}7\frac{27}{13}=\frac{27648}{91}.           \tag{17}
 \end{aligned}
\]

Global cubes contribute `O(B^(1/3))`: (15) bounds both `|m|` and `n` by `C sqrt(B)`, and a reduced cube has numerator and denominator cubes. Removing them is negligible. FMA4 then identifies the remaining parameters with Hasse failures, exactly twice per target. Thus the parameter-to-target factor one-half cancels the factor two in (14). Combining (10) and (17) yields

\[
 C_2=\frac{K_0^2E(2)}3\mathcal A_2W
     =\frac{256\mathcal A_2}{13\sqrt3\,\zeta(2)}
       \prod_{p\equiv2\ (3)}\frac{1+p^{-3}}{1-p^{-3}},
\]

as claimed. The region is bounded by (15), has positive area, and its boundary is a finite union of quadratic arcs and coordinate lines. It is therefore Jordan measurable. ∎

## 5. Lemma SLC3: an explicit archimedean integral

Set

\[
 R_0=768\sqrt6,\qquad L=\frac{80640}{19},\qquad
 U=\frac{91392}{13},
\]

and define

\[
 M(R)=\begin{cases}
 16(5376-R),&R_0\le R\le L,\\
 3(R+1792),&L\le R\le U,\\
 16(R-5376),&U\le R.
 \end{cases}
\]

Then an alternative exact area formula is

\[
 \mathcal A_2=\int_{R_0}^{\infty}
 \left(\frac1{16(R+5376)}+\frac1{M(R)}\right)
                 \frac{dR}{\sqrt{R^2-R_0^2}}.           \tag{18}
\]

The first summand is the negative-parameter area `A_(2,-)` and the second is the positive-parameter area `A_(2,+)`. In particular both signs contribute positive constants, obtained from (4) by replacing `A_2` with the corresponding area.

**Proof.** Substituting `x=ty` in the area integral gives (5). For `t<0`, the second expression in (3) dominates the other two and equals `16(Dt^2+K+5376|t|)`. For `t>0`, set `R=Dt+K/t`; the arithmetic-geometric mean inequality gives `R>=2sqrt(DK)=R_0`. The third expression in (3) dominates the first because `R_0>1280`. Comparing the remaining two expressions gives breakpoints `16(5376-R)=3(R+1792)` and `16(R-5376)=3(R+1792)`, which are precisely `L,U`. Each `R>R_0` has two positive inverse branches, interchanged by `t -> k/t`, and on either branch

\[
 \frac{|dt|}{|t|}=\frac{dR}{\sqrt{R^2-R_0^2}}.
\]

The two branches cancel the factor one-half in (5), proving (18). ∎

The integral is elementary. For a closed expression, put

\[
 q=\sqrt{\frac{5376-R_0}{5376+R_0}},\quad
 u_L=\sqrt{\frac{L-R_0}{L+R_0}},\quad
 u_U=\sqrt{\frac{U-R_0}{U+R_0}},\quad
 \eta=\sqrt{\frac{R_0-1792}{R_0+1792}}.
\]

Substituting `u=sqrt((R-R_0)/(R+R_0))` in each interval of (18) gives

\[
 \mathcal A_2=
 \frac1{12288\sqrt{43}}
  \log\left(\frac{u_U+q}{u_U-q}\frac{q+u_L}{q-u_L}\right)
 +\frac{\arctan(u_U\eta)-\arctan(u_L\eta)}{384\sqrt5}.
                                                        \tag{19}
\]

Here `0<u_L<q<u_U<1`, so the logarithm is real. Formula (18) remains the definition needed for (4); the checker compares (18), (19), and a polar-coordinate evaluation of the original area as numerical controls, separately from the exact symbolic breakpoint checks.

## 6. Corollary SLC4: local solubility has the same leading constant

Let `L_2(B)` count all everywhere locally soluble fibres on (1) in height at most `B`, and let `G_2(B)` count fibres having a rational point. Then

\[
 G_2(B)=O(B^{2/3}),\qquad
 L_2(B)=N_2(B)+G_2(B)\sim C_2B/\log B.                  \tag{20}
\]

**Proof.** The fixed quadratic factor has no rational root. Thus the cubic factor `U^3-3gamma U-r` has a rational root exactly when `r=t^3-(9/32)t` for `t in Q`. Write `t=m/n` in lowest terms, `n>0`. The integer pair representing `r` is

\[
 [32m^3-9mn^2:32n^3].                                  \tag{21}
\]

Its content divides `32`: an odd prime dividing the denominator divides `n` and cannot divide the numerator; the two-adic valuation is at most five because the possible valuations `5` and `2v_2(n)` cannot tie when `n` is even, and if `n` is odd the denominator has valuation exactly five. With `M=max(|m|,n)`, the maximum absolute entry in (21) is at least `23M^3` and at most `41M^3`. For `n>=|m|` use the denominator for the lower bound; for `|m|>n`, use `|m|(32m^2-9n^2)>=23|m|^3`. Hence `H(r)` is comparable to `H(t)^3`. The fixed linear embedding (1) makes height at most `B` imply `H(t)<<B^(1/3)`, giving `O(B^(2/3))` rational parameters and therefore at most that many targets. This estimate also includes `t=0`.

Every rational-point fibre is locally soluble; by definition the remaining locally soluble fibres are precisely Hasse failures. Since `B^(2/3)=o(B/log B)`, (4) gives (20). ∎

## 7. Theorem SLC5: global-solubility leading constant and conditional frequency

Define the homogeneous cubic height function

\[
 Q_3(x,y)=\max\left\{36y^3,
 |336y^3-2048x^3+576xy^2|,
 |21y^3+384x^3-108xy^2|\right\}\quad(y\ge0),
\]

and put

\[
 \Gamma_3=\{(x,y):y>0,Q_3(x,y)\le1\},\qquad
 \mathcal A_3=\operatorname{area}(\Gamma_3)
    =\frac12\int_{\mathbb R}Q_3(t,1)^{-2/3}\,dt.
                                                        \tag{22}
\]

Then the rational-point count on the entire second line has the asymptotic

\[
 \boxed{G_2(B)=C_{G,2}B^{2/3}+O(B^{1/3}\log B)},\qquad
 \boxed{C_{G,2}=\frac{\mathcal A_3}{3\zeta(2)}
         (1+2^{5/3})(3+3^{2/3})}.                       \tag{23}
\]

Consequently the conditional frequency of Hasse failures among everywhere locally soluble targets is

\[
 \frac{N_2(B)}{L_2(B)}
 =1-\frac{C_{G,2}}{C_2}\frac{\log B}{B^{1/3}}(1+o(1)).  \tag{24}
\]

**Proof.** The rational-point parameterization from SLC4 is `r=f(t)=t^3-(9/32)t`. It is injective on `Q`. If `f(t)=f(s)` for distinct rational `t,s`, division by `t-s` gives

\[
 t^2+ts+s^2=9/32.
\]

The left side has even two-adic valuation. To see this directly, extract `2^e`, where `e=min(v_2(t),v_2(s))`; the normalized pair has two-adic integral coordinates, at least one odd, and its quadratic form is odd modulo two. The valuation is therefore `2e`. This contradicts the valuation `-5` of the right side. The case of a zero coordinate is covered by taking the minimum of the two nonzero/finite valuations. Thus every rational-point target has exactly one rational parameter `t`, and no division by three is appropriate.

For `t=m/n` in lowest terms, `n>0`, a projective integer representative of its target is

\[
 [18n^3:36n^3:
   336n^3-2048m^3+576mn^2:
   -21n^3-384m^3+108mn^2].                              \tag{25}
\]

Its common divisor is

\[
 d(m,n)=2^{e(n)}3^{\mathbf1_{3\mid m}},\qquad
 e(n)=\begin{cases}
 0,&n\text{ odd},\\
 3,&v_2(n)=1,\\
 7,&v_2(n)\ge2.
 \end{cases}                                          \tag{26}
\]

Here the divisibility convention includes `3|0`. For `m=0` coprimality forces `n=1`; the representative has content three and its target has height `112`, contributing only one target.

We verify (26) for `m!=0`. A common prime other than two or three would divide `n` through the first coordinate and would then fail to divide the third coordinate, which reduces to `-2048m^3`. At three, if `m` is a unit the third coordinate is a unit; if `3|m`, `n` is a unit and both last coordinates have valuation exactly one. At two, if `n` is odd the last coordinate is odd. If `v_2(n)=1`, the last coordinate has valuation three, less than the valuations of the others. If `v_2(n)>=3`, its `-384m^3` term has unique least valuation seven. At the remaining case `n=4s`, `m,s` odd, divide the last coordinate by `64`; the result is

\[
 -21s^3-6m^3+27ms^2,
\]

which is even, while the first coordinate has valuation exactly seven and the third has valuation at least ten. Thus the common valuation is seven here as well. This proves (26), including the extra factor two from that cancellation.

The target height is exactly `Q_3(m,n)/d(m,n)`. The region `Gamma_3` is bounded: the term `36y^3` controls `y`, and when `|x|>=y` the absolute value in the third coordinate is at least `1136|x|^3`. Thus `Q_3(x,y)>=36max(|x|,y)^3`. Its fixed algebraic boundary is a finite union of coordinate-monotone arcs; its integer lattice count in any fixed union of residue classes is area times the lattice density plus `O(T+1)` after dilation by `T`. This follows by counting unit grid squares intersected by the finitely many monotone arcs; all other squares are entirely inside or outside.

There are three local cases at two and two at three, all determined modulo twelve. Each such local stratum has constant `d`, so its height region is `(dB)^(1/3) Gamma_3`. The local strata impose coprimality at two and three. Use Möbius inversion only for common divisors coprime to six. Multiplication by such a divisor preserves each local stratum and permutes its residue classes. For each divisor `h`, the lattice error is `O(B^(1/3)/h+1)`; the bounded region allows truncation at `h<=CB^(1/3)`. Summation gives `O(B^(1/3)log B)`. Extending the main sum to infinity changes it by `O(B^(1/3))`, and its density factor is

\[
 \prod_{p\ge5}(1-p^{-2})=\frac3{2\zeta(2)}.
\]

For clarity, the two-adic local densities and weights `d_2^(2/3)` are

\[
 \begin{array}{c|c|c}
 \text{local stratum}&\text{density}&\text{weight}\\\hline
 n\text{ odd}&1/2&1\\
 v_2(n)=1,\ m\text{ odd}&1/8&4\\
 v_2(n)\ge2,\ m\text{ odd}&1/8&2^{14/3}.
 \end{array}
\]

Their weighted sum is `1+2^(5/3)`. At three the cases `3 does not divide m` and `3|m,3 does not divide n` have densities `2/3,2/9` and weights `1,3^(2/3)`, respectively. Their weighted sum is `(2/9)(3+3^(2/3))`. Multiplying these two independent local sums, the outside-six primitive density and the area gives (23).

Finally `L_2=N_2+G_2`, and SLC1/SLC4 give `L_2(B)~C_2B/log B`. Division of (23) by this asymptotic proves (24). ∎
## 8. Verification and remaining scope

Run from the repository root:

```sh
.venv/bin/python research/arithmetic/verify_fixed_map_second_line_count.py --output research/arithmetic/fixed_map_second_line_count_checks.json
```

The checker verifies the rational coordinate identity, breakpoint and constant arithmetic, and bounded arithmetic controls of the content formula; it numerically compares independent archimedean area integrals. The written valuation proof covers all coprime pairs. Neither bounded arithmetic controls nor numerical quadrature prove the unrestricted asymptotic. The count imports exactly the previously scoped coprime-pair analytic estimate and the exact second-line classification.

This counts targets on one specified line, with the two rational parameters identified. It does not count all rational targets of `Phi` or distinct fibre-algebra isomorphism classes on this line. Those ambient and isomorphism-class problems retain separate hypotheses and height questions. The new constant does not strengthen the old first-line constant: the target lines and height growth laws are different.
