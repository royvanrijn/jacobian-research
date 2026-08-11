# F2 `k=1` fixed-coordinate Keller-pullback theorem

> **Status.**  The affine target changes used to obtain the four-parameter
> `k=1` normal form can be undone exactly with five scalar parameters.  This
> gives a denominator-free quintic in the fixed F2 target coordinates, an
> explicit carrier residue and eight-jet test, and four fixed-coordinate
> node fibers.  For every plane Keller map its pullback is reduced, and its
> complete affine singular and conductor schemes are the étale base changes
> of the target schemes.  Thus no unknown affine local type or local
> conductor length remains in the generic `k=1` packet; only four finite
> fiber counts remain to be evaluated from the source map.  The boundary
> factors and boundary point corrections still require the unresolved
> lower-Laurent map.

The coordinate identities and finite jet calculation are replayed by
[`verify_f2_affine_k1_keller_pullback.py`](../scripts/verify_f2_affine_k1_keller_pullback.py).

## 1. Restoring the fixed target coordinates

Let

\[
 p_0(t)=t^3+at,
 \qquad
 q_0(t)=t^5+bt^4+ct^2+dt                         \tag{1.1}
\]

be the normal form of the generic `k=1` target.  Undoing the target
translations, independent scalings, and triangular shear used in its
construction gives every fixed-coordinate parametrization in the form

\[
 \boxed{
 \begin{aligned}
 p(t)&=P_0+A p_0(t),\\
 q(t)&=Q_0+B q_0(t)+\Gamma p_0(t),
 \end{aligned}
 \qquad AB\ne0.}                                  \tag{1.2}
\]

The affine change of normalization parameter has already been absorbed in
`a,b,c,d`.  Conversely, (1.2) is taken back to (1.1) by

\[
 U=\frac{P-P_0}{A},\qquad
 V=\frac{A(Q-Q_0)-\Gamma(P-P_0)}{AB}.             \tag{1.3}
\]

Let `F(U,V)` be the twelve-support implicit quintic from
[`F2_AFFINE_TARGET_K1_IMPLICIT_CONDUCTOR.md`](F2_AFFINE_TARGET_K1_IMPLICIT_CONDUCTOR.md).
Then

\[
 \boxed{
 G(P,Q)=A^5B^3 F\!\left(
 \frac{P-P_0}{A},
 \frac{A(Q-Q_0)-\Gamma(P-P_0)}{AB}
 \right)}                                         \tag{1.4}
\]

is a polynomial, not merely a rational expression, and
`G(p(t),q(t))=0`.  Its degree is five and its top homogeneous part is
`B^3P^5`.  Formula (1.4) closes the inverse target-normalization step: a
source compiler can work directly in the fixed F2 coordinates without
recomputing an implicit equation.

## 2. Exact carrier residue and jet interface

At the target puncture put `u=1/t`.  Equation (1.2) gives

\[
\begin{aligned}
p&=Au^{-3}\left(1+au^2+\frac{P_0}{A}u^3\right),\\
-q&=-Bu^{-5}\left(
1+bu+\frac{\Gamma}{B}u^2+cu^3
+\left(d+\frac{\Gamma a}{B}\right)u^4
+\frac{Q_0}{B}u^5\right).
\end{aligned}                                      \tag{2.1}
\]

Hence the fixed carrier residue is

\[
 \boxed{\lambda=\frac{A^5}{(-B)^3}.}              \tag{2.2}
\]

Writing

\[
 \pi=\frac{P^3}{(-Q)^2},\qquad
 h=\frac{P^5}{(-Q)^3},\qquad
 w=h-\lambda-\sum_{j=1}^7c_j\pi^j,               \tag{2.3}
\]

all coefficients of `w` through `u^8` are explicit rational polynomials in

\[
 A^{\pm1},B^{\pm1},P_0,Q_0,\Gamma,a,b,c,d,c_1,\ldots,c_7. \tag{2.4}
\]

In particular,

\[
 [u]w=-3\lambda b-c_1\frac{A^3}{B^2}.             \tag{2.5}
\]

On `lambda=125/729`, the contact number is the first nonzero coefficient in
this eight-term list, truncated at eight.  Thus restoring the fixed target
chart introduces no infinite or unspecified jet problem.

## 3. Universal Keller-pullback theorem

Let `k` be algebraically closed of characteristic zero and let

\[
 \Phi=(P_s,Q_s):\mathbb A^2_{x,y}\longrightarrow\mathbb A^2_{P,Q}
\]

have nonzero constant Jacobian.  Put

\[
 H(x,y)=G(P_s(x,y),Q_s(x,y)),\qquad D=V(H).        \tag{3.1}
\]

### Theorem 3.1

For every reduced target curve `C=V(G)`:

1. `Phi` is étale and `D=C times_(A2) A2` is reduced; equivalently, `H` is
   squarefree in `k[x,y]`;
2. the chain rule and the Keller determinant give the equality of ideals
   \[
   \boxed{(H_x,H_y)=
   (G_P(P_s,Q_s),G_Q(P_s,Q_s));}                  \tag{3.2}
   \]
3. scheme-theoretically,
   \[
   \boxed{\operatorname{Sing}(D)=
   \operatorname{Sing}(C)\times_C D;}            \tag{3.3}
   \]
4. normalization, its quotient, and the conductor commute with this étale
   base change.

#### Proof

The constant nonzero Jacobian makes `Phi` étale by the Jacobian criterion.
Étale base change preserves reduced schemes, proving (1).  The chain rule is

\[
 \binom{H_x}{H_y}=
 \begin{pmatrix}P_{s,x}&Q_{s,x}\\P_{s,y}&Q_{s,y}\end{pmatrix}
 \binom{G_P}{G_Q}\!\bigg|_{(P_s,Q_s)}.            \tag{3.4}
\]

The matrix in (3.4) is invertible over `k[x,y]` because its determinant is a
unit, giving (3.2).  Reducing (3.2) modulo `H` proves (3.3).  Normalization
commutes with smooth, hence étale, base change.  Flatness of the étale map
then base-changes the normalization exact sequence; locally the annihilator
of its finite quotient is the base-changed conductor.  This proves (4).
\(\square\)

## 4. Consequence for the generic four-node quintic

Let `r` run through the roots of the collision quartic

\[
 R(r)=r^4+br^3+ar^2+(2ab-c)r-(a^2+d).            \tag{4.1}
\]

The four target nodes in normalized coordinates are

\[
 X(r)=-r(r^2+a),\qquad
 Y(r)=(r^2+a)(r^3+2ar+ab-c).                     \tag{4.2}
\]

Their fixed-coordinate values are therefore

\[
 \boxed{
 (P_r,Q_r)=
 (P_0+AX(r),\ Q_0+BY(r)+\Gamma X(r)).}           \tag{4.3}
\]

For a generic four-node target, every affine source singularity of `D` is an
ordinary node lying in one of the four explicit fibers

\[
 P_s=P_r,\qquad Q_s=Q_r.                         \tag{4.4}
\]

Each geometric solution of (4.4) has normalization-quotient length one and
conductor-divisor degree two.  There are no other affine singular or
conductor local types.  The total normalization-defect contribution is the
number of geometric solutions of the four fibers in (4.4), counted once
each.  These fibers may be counted directly from the eventual source
polynomials; factorization of `H` is not needed for that finite calculation.

There is also an unconditional finite-flat bound.  Let `B` be the
normalization of `k[P,Q]` in `k(x,y)`.  Its finite cover has rank equal to the
geometric degree `d`, and the affine source is its étale open subset.  For a
node value `y_r`, write `N_r` for the number of affine solutions of (4.4)
and `L_r^partial` for the scheme length in the missing boundary.  Flatness
and affine étaleness give

\[
 \boxed{d=N_r+L_r^\partial.}                     \tag{4.5}
\]

Every point of the nonproperness curve has a nonempty boundary fiber, so
`L_r^partial>=1`.  Consequently

\[
 \boxed{0\le N_r\le d-1,\qquad
 \ell_{\rm aff}^{\rm node}=\sum_{r=1}^4N_r\le4(d-1).} \tag{4.6}
\]

At the squarefree/double degree floors this gives affine nodal
normalization-defect bounds `20/44` and conductor-divisor bounds `40/88`.
Formula (4.5) also identifies the missing datum exactly: the four affine
counts are complementary to four boundary-fiber lengths.

This is an important separation.  The still-unknown point-supported terms
in the logarithmic boundary calculation are genuinely boundary terms, not
uncontrolled affine node terms.

The subsequent
[`all-stratum conductor-conservation theorem`](F2_AFFINE_TARGET_K1_CONDUCTOR_CONSERVATION.md)
extends this conclusion across the collision discriminant.  It proves that
every `k=1` degeneration has affine delta four and the same degree-eight
conductor divisor, while replacing the four node counts by finitely many
delta-weighted singular-fiber counts.

<!-- status-consumer: PF2K1CC1 f152c82ef2d54c32 -->

The
[`affine strict-log-étale resolution theorem`](AFFINE_KELLER_STRICT_LOG_ETALE_RESOLUTION.md)
identifies the precise role of those counts.  After any embedded resolution,
the Keller pullback is a strict étale map of curve-log pairs, so its relative
logarithmic cokernel is zero for every affine singularity type.  The
normalization/conductor lengths above constrain finite fibers and boundary
escape; they are not finite point terms of the relative log module.

<!-- status-consumer: PAER1 60eb24b2232d159e -->

## 5. Remaining source input

The generic `k=1` pullback compiler now has a fixed input contract:

1. supply the complete fixed-coordinate Laurent polynomials `P_s,Q_s`;
2. form the single denominator-free expression (1.4);
3. compute its valuations and factors at the unresolved source boundary;
4. use (4.4) for all affine node/conductor lengths; and
5. determine the boundary component's `(e,f,E^2)` and boundary-local
   nonunit-`Fitt_1` corrections.

The currently certified F2 materials do not supply item 1: they give corner
envelopes, selected Laurent bands, and conditional/common-power systems, but
not the complete normalized polynomial pair.  Therefore this theorem closes
the inverse-normalization and affine-correction gaps, but it does not locate
the purity divisor, factor the boundary pullback, exclude `(75,125)`, or
prove `JC(2)`.

## Sources

- [Stacks Project, Lemma 35.18.1](https://stacks.math.columbia.edu/tag/034E),
  for preservation of reducedness under smooth, hence étale, maps.
- [Stacks Project, Lemma 37.19.2](https://stacks.math.columbia.edu/tag/03GV),
  for normalization commuting with smooth base change.

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_affine_k1_keller_pullback.py
```
