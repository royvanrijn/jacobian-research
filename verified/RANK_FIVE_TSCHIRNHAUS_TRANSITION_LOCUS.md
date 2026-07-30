# Rank-five ambient and marked Tschirnhaus transition loci

This note begins the marked-versus-ambient continuation of the generic
Tschirnhaus non-descent theorem.  It makes the rank-five transition loci
exact and reduces every still-possible marked nonprojective transport to a
self-equivalence orbit problem for one fixed quintic Keller map.

Work over an algebraically closed characteristic-zero field.  Fix the split
presentation

\[
 r=(1,2,3,4,5)
\]

with root polynomial

\[
 R(T)=\prod_{i=1}^5(T-i)
 =T^5-15T^4+85T^3-225T^2+274T-120.                    \tag{1.1}
\]

All statements are on the clean locus where the root polynomial is
squarefree and the linear, cubic, and quartic coefficients are nonzero.

## 1. The ambient stable-equivalence hypersurface

Let

\[
 Q(T)
 =T^5+A_4T^4+A_3T^3+A_2T^2+A_1T+A_0.                 \tag{1.2}
\]

For a monic quintic the compiler-slice stable invariant is

\[
 I(Q)=\frac{A_1^2}{A_3A_4^6}.                          \tag{1.3}
\]

For (1.1),

\[
 I(R)=\frac{75076}{968203125}.                         \tag{1.4}
\]

The exact quadratic-gauge stable-moduli theorem therefore gives

\[
\boxed{
 F_Q\sim_{\mathrm{stable}}F_R
 \quad\Longleftrightarrow\quad
 968203125A_1^2-75076A_3A_4^6=0.
}                                                       \tag{1.5}
\]

Equivalently,

\[
 A_3=\frac{968203125A_1^2}{75076A_4^6}.                \tag{1.6}
\]

Thus the clean ambient stable-equivalence locus is a rational
four-dimensional hypersurface with coordinates
`(A_0,A_1,A_2,A_4)`.

In ordered-root coordinates `u=(u_1,\ldots,u_5)`, write `e_j=e_j(u)`.
Then

\[
 Q(T)=T^5-e_1T^4+e_2T^3-e_3T^2+e_4T-e_5,
\]

and (1.5) becomes

\[
\boxed{
 968203125e_4^2-75076e_2e_1^6=0.
}                                                       \tag{1.7}
\]

## 2. The labelled projective locus

The unique fractional-linear transformation matching

\[
 1\mapsto u_1,\qquad2\mapsto u_2,\qquad3\mapsto u_3
\]

also matches the fourth and fifth roots exactly when

\[
\begin{aligned}
 D_4={}&u_1u_2-4u_1u_3+3u_1u_4
       +3u_2u_3-4u_2u_4+u_3u_4=0,\\
 D_5={}&2(u_1u_2-3u_1u_3+2u_1u_5
       +2u_2u_3-3u_2u_5+u_3u_5)=0.                    \tag{2.1}
\end{aligned}
\]

These are the two framed residuals from the all-rank projective theorem.
At `u=r`, their Jacobian has rank two.  The Jacobian of (1.7), `D_4`, and
`D_5` has rank three there.  Hence:

\[
\begin{array}{c|c|c}
\text{locus in the five-dimensional ordered-root chart}
&\text{local codimension at }r&\text{local dimension}\\ \hline
\text{ambient stable equivalence}&1&4\\
\text{projective Tschirnhaus transport}&2&3\\
\text{their intersection}&3&2.
\end{array}                                             \tag{2.2}
\]

In particular, neither the ambient stable-equivalence locus nor the
projective locus contains the other.  There is no quotient
`\mathcal E_5/PGL_2` to take without first imposing an actual group action
preserving `\mathcal E_5`.

## 3. The explicit ambient equivalence

Normalize (1.2) by its linear coefficient.  Its compiler target and seeds
are

\[
\begin{aligned}
 \pi_Q&=\frac{A_3}{A_1},&
 b_Q&=\frac{A_2}{A_1},&
 c_Q&=-\frac{2A_0}{A_1},\\
 v_4&=\frac{A_4A_1^3}{A_3^4},&
 v_5&=\frac{A_1^4}{A_3^5}.
                                                               \tag{3.1}
\end{aligned}
\]

Let `(u_4,u_5)` be the seeds compiled from `R`.  On (1.5), put

\[
\boxed{
 \alpha=-\frac{1275A_1}{274A_3A_4}.
}                                                       \tag{3.2}
\]

Direct substitution of (1.6) gives

\[
 v_4=\alpha^5u_4,\qquad v_5=\alpha^6u_5.               \tag{3.3}
\]

Thus the stable equivalence is already the ordinary coefficient-torus
equivalence.  In source and normalized target coordinates it is

\[
 \sigma_\alpha(x,y,z)
 =(\alpha x,\alpha^{-1}y,\alpha^{-2}z),                \tag{3.4}
\]

\[
 \tau_\alpha(\pi,b,c)
 =(\alpha^{-2}\pi,\alpha^{-1}b,\alpha c),              \tag{3.5}
\]

and

\[
 F_Q\circ\sigma_\alpha=\tau_\alpha\circ F_R.           \tag{3.6}
\]

## 4. When the canonical equivalence carries the marked fibre

The selected target for (1.1) is

\[
 y_R=
 \left(
 \frac{85}{274},
 -\frac{225}{274},
 \frac{120}{137}
 \right).                                               \tag{4.1}
\]

The canonical equivalence carries the selected fibre of `R` to that of `Q`
exactly when

\[
 y_Q=\tau_\alpha(y_R).                                 \tag{4.2}
\]

In elementary-symmetric root coordinates, (3.2) and (4.2) are

\[
\begin{aligned}
 274\alpha e_1e_2&=1275e_4,\\
 274\alpha^2e_2&=85e_4,\\
 274\alpha e_3&=225e_4,\\
 274e_5&=120\alpha e_4.                                \tag{4.3}
\end{aligned}
\]

Together with the ambient equation (1.7), these imply

\[
\boxed{
 e_1=15\alpha,\quad
 e_2=85\alpha^2,\quad
 e_3=225\alpha^3,\quad
 e_4=274\alpha^4,\quad
 e_5=120\alpha^5.
}                                                       \tag{4.4}
\]

Indeed the first two equations give `e_1=15*alpha`.  The second equation
gives `e_4=(274/85)alpha^2e_2`; after these substitutions, (1.7) is a
nonzero clean-locus multiple of

\[
 e_2(e_2-85\alpha^2).
\]

The last two equations then recover `e_3` and `e_5`.

Consequently

\[
 Q(T)=\prod_{i=1}^5(T-\alpha i).                       \tag{4.5}
\]

On the labelled frame, `u_i=\alpha i`.  Hence:

> **Canonical marked-transport theorem.**
>
> On the clean rank-five ambient stable-equivalence hypersurface, the
> explicit coefficient-torus equivalence carries the selected complete
> fibre and its ordered roots exactly on the one-dimensional scaling locus
> `u_i=alpha*i`.

The five equations consisting of (1.7) and (4.3), in the six variables
`(u_1,\ldots,u_5,alpha)`, have Jacobian rank five at `(r,1)`.  The scaling
locus is therefore locally the complete marked locus for the canonical
transport.

## 5. Exact reduction of arbitrary marked transport

The preceding theorem classifies the marked behaviour of the canonical
ambient equivalence, not every possible stable self-equivalence of the
quintic map.

For a point `Q` of (1.5), pull its selected target back to the fixed map:

\[
\boxed{
 z_Q=\tau_\alpha^{-1}(y_Q)
 =\left(
 \alpha^2\pi_Q,\,
 \alpha b_Q,\,
 \alpha^{-1}c_Q
 \right).
}                                                       \tag{5.1}
\]

Let `Stab_st(F_R)` be the stable polynomial left--right self-equivalence
group of `F_R`, and let `Stab_st^t(F_R)` denote its induced target action.
Composing any proposed equivalence `F_R -> F_Q` with the inverse of
(3.6) gives:

\[
\boxed{
 \begin{array}{c}
 \text{a stable marked equivalence from }(F_R,y_R)
 \text{ to }(F_Q,y_Q)\\
 \Longleftrightarrow\\
 z_Q\in Stab_{\mathrm{st}}^t(F_R)\cdot y_R.
 \end{array}
}                                                       \tag{5.2}
\]

Equation (5.2) is the exact remaining gate.  On the canonical scaling locus,
`z_Q=y_R`.  Away from it, every marked nonprojective lift requires a
nontrivial stable target self-equivalence of the one fixed map `F_R`.

The quadratic-gauge stable-moduli theorem shows that the displayed
coefficient-torus stabilizer is trivial in rank five: the seed weights are
five and six.  It does not classify the kernel of all polynomial stable
self-equivalences acting trivially on the intrinsic normalized boundary.
That kernel, and its target orbit of `y_R`, are the next problem.

## 6. All-rank form of the reduction

For every `N>=5`, equality of the complete fingerprint `Phi_N` is equivalent
to membership in one residual scaling orbit on the compiler seed torus.
Faithfulness of the weights makes the scaling parameter unique.  The same
target calculation shows:

\[
 H_u(S)=\alpha H_r(S/\alpha)                            \tag{6.1}
\]

exactly when the canonical ambient equivalence carries the selected target.
Thus the canonical marked locus is again the root-scaling locus.

For arbitrary marked transport, pull the second selected target back through
the unique canonical seed equivalence and test its orbit under the stable
self-equivalence group of the first map.  Rank five is the smallest case
with a positive-dimensional ambient stable quotient and is therefore the
correct stabilizer testbed.

The marked and unmarked dimension counts, and the resulting span over
`BS_N`, are organized in the
[clean quadratic-gauge decorated receiver](QUADRATIC_GAUGE_DECORATED_RECEIVER.md).

## Exact regression

Run

```bash
.venv/bin/python scripts/verify_rank_five_tschirnhaus_transition_locus.py
```

The checker verifies the ambient hypersurface, the two projective residuals,
all three Jacobian ranks in (2.2), the explicit seed scaling (3.3), the
canonical target equations, and the rank-five marked scaling locus.
