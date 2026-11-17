# Complete all-degree minimal-excess closure in the scalar `HC4` packet

## Status

This note upgrades the tail computations of `HC4RSD44--50` to a single
all-degree induction.

> **Theorem HC4RSD51 — complete all-degree `h=0` closure.**
> In the synchronized scalar reverse-Schur packet, every first transverse
> face of **minimal transverse excess**
> \[
> h=0
> \]
> has a fixed ruling.  Hence every such packet reduces to `HC2` or to the
> exact `JC2` cotangent endpoint.  There is no bound on the leading degree,
> no bound on the number of repeated roots, and no root-partition census.

Together with `HC4RSD25--26`, this means that any still-live scalar packet
must simultaneously have positive transverse excess `h>=1` and a genuinely
repeated leading-root architecture.

## 1. Setup

Write

\[
 f=A^2B,
 \qquad o=\deg B,
 \qquad k=\deg A=\frac{d-o}{2},
\]

where `B` is the squarefree product of roots having odd multiplicity.  At
`h=0`,

\[
 e=k+o,
\]

and admissibility `e<=d-2` gives

\[
 k\ge2.                                             \tag{1.1}
\]

Write the complete same-weight potential as

\[
 c=\sum_{i=0}^{J}\frac{z^i}{i!}R_i(x,y),           \tag{1.2}
\]

where

\[
 R_0=A^2B,
 \qquad R_1=aAB,
 \qquad R_2=\kappa a^2B,
\]

(up to the harmless normalization of `\kappa`) and, for every `i`,

\[
 \deg R_i=d-ik=o-(i-2)k.                           \tag{1.3}
\]

The case `o=0` is already HC4RSD45, so suppose `o>0`.

## 2. Two derivative-count identities

Let

\[
 \mathcal U(c)=\nabla c^{\mathsf T}
 \operatorname{adj}(\operatorname{Hess}c)\nabla c.
\]

Every monomial of `\mathcal U(c)` is quartic in copies of `c` and has two
properties:

1. the total number of `z` derivatives across the four copies is exactly two;
2. the total number of `y` derivatives across the four copies is exactly two.

For example, writing

\[
 p=c_x,\ q=c_y,\ r=c_z,
\]

and

\[
 A_0=c_{xx},\ B_0=c_{xy},\ C_0=c_{xz},
\ D_0=c_{yy},\ E_0=c_{yz},\ F_0=c_{zz},
\]

one has

\[
\begin{aligned}
\mathcal U(c)={}&A_0D_0r^2-2A_0E_0qr+A_0F_0q^2-B_0^2r^2
 +2B_0C_0qr+2B_0E_0pr-2B_0F_0pq\\
&-C_0^2q^2-2C_0D_0pr+2C_0E_0pq+D_0F_0p^2-E_0^2p^2.
\end{aligned}                                      \tag{2.1}
\]

Thus the coefficient of `z^N` only involves quadruples
`(R_{i_1},...,R_{i_4})` satisfying

\[
 i_1+i_2+i_3+i_4=N+2.                              \tag{2.2}
\]

## 3. Positive-degree highest tail

Suppose the highest nonzero tail `R_J` has degree `r>0`.  HC4RSD48 gives,
after a constant linear change,

\[
 R_J=x^r.                                           \tag{3.1}
\]

> **Tail-degree lemma.** For every `i<J`,
> \[
> \deg_y R_i\le J-i.                               \tag{3.2}
> \]

### Proof

Proceed downward by induction on `i`.  Put `\ell=J-i` and inspect the
coefficient

\[
 [z^{3J+i-2}]\mathcal U(c).                        \tag{3.3}
\]

By (2.2), the term containing the new coefficient `R_i` is necessarily the
linear term with three copies of `R_J`.  HC4RSD48 computes it exactly, up to
factorial scale, as

\[
 -Jr(r+J)x^{3r-2}(R_i)_{yy}.                       \tag{3.4}
\]

Every other contribution uses only indices larger than `i`.  If
`\delta_a=J-i_a`, then their deficits sum to `\ell`.  By induction their
`y`-degrees sum to at most `\ell`; by (2.1) exactly two `y` derivatives are
applied, so every nonlinear contribution has `y`-degree at most
`\ell-2`.

If `deg_y R_i>\ell`, the top `y`-term of `(R_i)_{yy}` has degree greater than
`\ell-2` and cannot cancel.  Since the coefficient in (3.4) is nonzero,
(3.2) follows.  QED.

## 4. Scalar highest tail

Suppose instead that the highest tail is a nonzero scalar

\[
 R_J=t.
\]

Let `I<J` be the first lower nonzero coefficient of positive degree.  In the
coefficient with two copies of `t` and two copies of `R_I`, the leading term
is

\[
 J^2t^2\det\operatorname{Hess}R_I.                 \tag{4.1}
\]

Therefore `R_I` has zero binary Hessian and, after a constant change,

\[
 R_I=x^r.
\]

Here (1.3) and `deg R_J=0` give

\[
 r=(J-I)k\ge2.                                     \tag{4.2}
\]

For any lower `R_i`, the coefficient with two scalar top copies, one `R_I`
copy and one `R_i` copy contains the nonzero linear term

\[
 J^2t^2r(r-1)x^{r-2}(R_i)_{yy}.                    \tag{4.3}
\]

The same deficit/y-derivative count now gives

\[
 \deg_yR_i\le I-i.                                 \tag{4.4}
\]

Thus scalar top tails do not create an exceptional branch; they merely reset
the induction at the first positive-degree descendant.

## 5. Contradiction when `kappa != 0`

If `\kappa\ne0`, then `R_2` is a nonzero scalar multiple of the squarefree
binary form `B`.

In the positive-degree-top case, (3.2) gives

\[
 \deg_yB\le J-2.                                   \tag{5.1}
\]

Because `deg R_J>0`, equation (1.3) gives

\[
 (J-2)k\le o-1.
\]

Together with `k>=2`,

\[
 J-2\le\left\lfloor\frac{o-1}{k}\right\rfloor<o-1
\]

whenever a higher tail exists.  But after any projective normalization a
squarefree binary form of degree `o` has

\[
 \deg_yB\ge o-1,                                   \tag{5.2}
\]

because at most one of its distinct linear factors can be `x`.  Equations
(5.1)--(5.2) contradict each other.

For a scalar top, (4.4) is stronger:

\[
 \deg_yB\le I-2\le J-3=\frac ok-1<o-1.
\]

Thus no higher tail exists.  The resulting quadratic-in-`z` face is exactly
the tail-free packet closed by HC4RSD46.

## 6. Contradiction when `kappa = 0`

Now `R_2=0`, but

\[
 R_1=aAB\ne0.
\]

The Pell sieve HC4RSD47 says that the odd-root count is

\[
 o=1,8,49,288,\ldots .                             \tag{6.1}
\]

If `o=1`, then `k>=2` means (1.3) admits no higher tail, so HC4RSD46 already
closes the packet.  Hence suppose `o>=8`.

For a positive-degree top, (3.2) gives

\[
 \deg_y(AB)\le J-1
 \le1+\left\lfloor\frac{o-1}{k}\right\rfloor
 \le1+\frac{o-1}{2}<o-1.                           \tag{6.2}
\]

But `B|AB` and `B` is squarefree of degree `o`, so

\[
 \deg_y(AB)\ge\deg_yB\ge o-1,
\]

contradiction.

For a scalar top, (4.4) gives

\[
 \deg_y(AB)\le I-1\le J-2=\frac ok\le\frac o2<o-1,
\]

again impossible.

This proves HC4RSD51.

## 7. Consequence for the scalar frontier

The scalar branch should no longer be organized by total degree.  The exact
remaining hierarchy is now:

1. `h=0`: **closed in every degree** by HC4RSD51;
2. `h>=1` with squarefree top: already closed by HC4RSD25;
3. `h>=1` with exactly one double root: already closed by HC4RSD26;
4. the live scalar frontier therefore requires `h>=1` **and** at least two
   double roots or one root of multiplicity at least three, subject to the
   resonance square equation of HC4RSD44.

This is a qualitative reduction: increasing the total degree alone no longer
creates new minimal-excess cases.

## 8. Verification

Run

```bash
.venv/bin/python scripts/verify_hc4_minimal_excess_complete.py
```

The checker verifies the derivative-count structure, the universal positive
highest-tail mixed coefficient, and both scalar-top starter coefficients.
The induction itself is symbolic degree bookkeeping and uses no fixed degree.
