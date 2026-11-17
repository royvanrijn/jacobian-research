# A Hasse-principle failure for a complete Keller fiber

This note gives one explicit rational target of one explicit polynomial
Keller map of affine three-space which lies in the image over every
completion of `Q`, but not over `Q`.  It is an arithmetic specialization of
the [tangent-map and weighted-suspension theorem](TANGENT_MAP_CORE.md).

## 1. The theorem

There is a polynomial map

\[
 F:\mathbb A^3_{\mathbb Q}\longrightarrow\mathbb A^3_{\mathbb Q}
\]

with constant Jacobian determinant `-38` and a target

\[
 \boxed{q=(12138,-308652,1)}                              \tag{1}
\]

such that

\[
 F^{-1}(q)(\mathbb Q)=\varnothing,
\]

but

\[
 F^{-1}(q)(\mathbb R)\ne\varnothing,
 \qquad
 F^{-1}(q)(\mathbb Q_p)\ne\varnothing
 \quad\text{for every prime }p.                         \tag{2}
\]

The complete fiber is the finite etale scheme

\[
 \boxed{F^{-1}(q)\simeq\operatorname{Spec}\mathbb Q[W]/(P(W))}, \tag{3}
\]

where

\[
\begin{aligned}
 P(W)={}&(W^2-2)(W^2-17)(W^2-34)
          (157W^2-267W+399)\\
 ={}&157W^8-267W^7-7922W^6+14151W^5+85613W^4\\
 &\quad-181560W^3+89828W^2+308652W-461244.
                                                               \tag{4}
\end{aligned}
\]

In particular, (2) is a Hasse-principle failure for a complete regular
zero-dimensional Keller fiber.

## 2. The intersective polynomial

Put

\[
 Q(W)=(W^2-2)(W^2-17)(W^2-34),
 \qquad
 R(W)=157W^2-267W+399.                                  \tag{5}
\]

### No rational root

None of `2`, `17`, and `34` is a square in `Q`.  Moreover,

\[
 \operatorname{disc}(R)=(-267)^2-4\cdot157\cdot399
 =-179283<0.                                             \tag{6}
\]

Thus neither `Q` nor `R` has a rational root, and hence `P=QR` has no
rational root.  All six roots of `Q` are real, while `R` is positive
definite, so `P` has exactly six real roots.

### A root over every `Q_p`

At `p=2`, the unit `17` is a square in `Q_2` because

\[
 17\equiv1\pmod 8.                                      \tag{7}
\]

At `p=17`, the integer `2` is a square because

\[
 6^2\equiv2\pmod {17},                                  \tag{8}
\]

and the root lifts by Hensel's lemma.

Now let `p` be an odd prime different from `17`.  If either `2` or `17` is
a quadratic residue modulo `p`, the corresponding quadratic factor in (5)
has a simple root.  If both are nonresidues, multiplicativity of the Legendre
symbol gives

\[
 \left(\frac{34}{p}\right)
 =\left(\frac2p\right)\left(\frac{17}p\right)=1.         \tag{9}
\]

Thus one of `2`, `17`, and `34` has a simple square root modulo `p`, which
again lifts to `Q_p`.  This proves the nonarchimedean part of (2).  These
roots are all `p`-adic integers, so the same argument plus the Chinese
remainder theorem shows that `P` has a root modulo every positive integer.

The factors in (4) are separable and pairwise coprime.  The polynomial `P`
is therefore squarefree.

## 3. Forcing the admissible normalization

The quadratic `R` was selected to impose the tangent identity

\[
 P(1)-P(0)=P'(0).                                       \tag{10}
\]

Indeed,

\[
 P(0)=-461244,
 \qquad P'(0)=308652,
 \qquad P(1)=-152592.                                   \tag{11}
\]

Remove the affine part:

\[
\begin{aligned}
 H(W)&=P(W)-P(0)-P'(0)W\\
 &=157W^8-267W^7-7922W^6+14151W^5+85613W^4
   -181560W^3+89828W^2.                                 \tag{12}
\end{aligned}
\]

Then

\[
 H(0)=H'(0)=H(1)=0,
 \qquad H'(1)=38.                                       \tag{13}
\]

Take

\[
 c=-38.
\]

The remaining weighted-chart condition is also nondegenerate:

\[
 H''(1)=160590,
 \qquad
 \kappa=\frac{H''(1)}c=-\frac{80295}{19}\ne-2.          \tag{14}
\]

Consequently `H` is an admissible weighted primitive.  Its source parameter
is

\[
 a_0=-\frac{1+\kappa}{2+\kappa}
 =-\frac{80276}{80257}.                                 \tag{15}
\]

The choice of `R` is elementary.  Since

\[
 Q(0)=-1156,\qquad Q(1)=-528,
\]

a quadratic `rW^2+lW+m` makes `QR` satisfy (10) precisely when

\[
 -528r+628(l+m)=0.                                      \tag{16}
\]

The solution `(r,l,m)=(157,-267,399)` has negative discriminant and also
makes the resulting Jacobian constant the small integer `-38`.

## 4. The explicit polynomial map

Write `h=H'` and set

\[
 v=xy,qquad S=x^2z,qquad u=1+v,
\]

\[
 \gamma=1-\frac{80276}{80257}v+S,
 \qquad W=u\gamma.                                      \tag{17}
\]

Define

\[
\begin{aligned}
 F_1&=\frac{1}{x^2}
 \left(u+\frac{Wh(W)-H(W)}{c\gamma^2}\right),\\[2mm]
 F_2&=\frac{1}{x}
 \left(c+\frac{h(W)}\gamma\right),\\
 F_3&=x\gamma,
 \qquad F=(F_1,F_2,F_3).                               \tag{18}
\end{aligned}
\]

Formula (18) is a compact polynomial formula, not merely a rational map.
Since `H` begins in degree two and `h` begins in degree one, substitution of
`W=u gamma` first cancels every displayed power of `gamma`.  The numerator
of `F_2`, regarded as a polynomial in `(v,S)`, lies in `(v,S)` and hence is
divisible by `x`.  The numerator of `F_1` lies in `(v^2,S)` and hence is
divisible by `x^2`.  The vanishing of its coefficient of `v` is exactly
equation (15).

Thus all three coordinates belong to `Q[x,y,z]`.  Their total degrees and
numbers of monomials are respectively

\[
 (32,82),\qquad(31,76),\qquad(4,3).                     \tag{19}
\]

The weighted-suspension determinant identity gives

\[
 \boxed{\det DF=c=-38}.                                 \tag{20}
\]

The exact checker also expands (18) and verifies (20) directly.

## 5. Identification of the complete fiber

For a target `(A,B,C)`, the inverse pencil of (18) is

\[
 E_{A,B,C}(W)=H(W)-BCW+cAC^2.                           \tag{21}
\]

At the target (1),

\[
\begin{aligned}
 E_q(W)
 &=H(W)+308652W-38\cdot12138\\
 &=H(W)+308652W-461244\\
 &=P(W).                                                \tag{22}
\end{aligned}
\]

Because `P` is squarefree, every root is in the regular reconstruction
chart.  Explicitly, for a root `r` in any characteristic-zero field, put

\[
 \gamma_r=-\frac{P'(r)}c=\frac{P'(r)}{38},
 \qquad x_r=\frac1{\gamma_r},
 \qquad y_r=r-\gamma_r,                                 \tag{23}
\]

and

\[
 z_r=\frac{\gamma_r-1-a_0(r/\gamma_r-1)}{x_r^2}.        \tag{24}
\]

Then `F(x_r,y_r,z_r)=q`.  Conversely, a point of `F^{-1}(q)` has

\[
 F_3=x\gamma=1.                                         \tag{25}
\]

Hence neither `x` nor `gamma` vanishes, and the suspension identities force
`W` to satisfy (22).  Equations (23)--(24) then recover the source point
uniquely.  This proves the complete scheme identification (3).

Equation (25) is also the boundary-completeness argument: no `x=0`,
`gamma=0`, or `C=0` boundary point can contribute another preimage.

Finally, a rational source point would give a rational value of the source
polynomial `W`, hence a rational root of `P`, which does not exist.  Every
local root constructed in Section 2 has nonzero derivative and therefore
reconstructs a point over the same completion.  This completes the proof of
(2).

## 6. Why the special fiber is reducible

Reducibility is necessary here, rather than a weakness of the example.  If a
finite etale `Q`-scheme were connected of degree greater than one, its
transitive Galois action would contain a derangement.  Chebotarev would then
produce unramified primes at which the scheme has no local point.  Thus an
everywhere locally soluble finite etale scheme without a rational point must
be disconnected.

This does not change the generic monodromy: the universal pencil theorem
applied to (12) gives generic monodromy `S_8` for `H(W)-sW+t`.  Only the
distinguished Hasse fiber (22) is reducible.

## Verification

Run

```bash
.venv/bin/python scripts/verify_hasse_keller_fiber.py
```

The checker audits the residue-class covering, exceptional Hensel lifts,
absence of rational roots, squarefreeness, normalization, admissibility,
polynomial coordinate expansion, constant Jacobian, suspension square,
inverse pencil, and quotient-ring reconstruction identities.
