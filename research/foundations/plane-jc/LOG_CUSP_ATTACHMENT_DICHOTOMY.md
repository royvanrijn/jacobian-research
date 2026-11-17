# Cusp boundary attachments: smooth folds versus SNC nodes

> **Status.**  Exact local theorem and exact countermodel to a universal
> `2f` cusp charge.  At a boundary point above a unibranch target value of
> multiplicity `m_C` and residue index `q_p`, let `epsilon_p=0` when the
> point is smooth on the reduced source boundary and `epsilon_p=1` when it
> is an SNC boundary node.  If the logarithmic `Fitt_1` defect is isolated,
> its colength is at least
> `q_p*m_C-1+epsilon_p`.  Over a complete residue-degree-`f` fiber with `h`
> geometric points, `c` of them boundary nodes, the resulting lower ledger
> is `m_C*f-h+c`.  For an ordinary cusp this lies between `f` and `2f`.
> The upper endpoint is the node-saturated packet of
> [`LOG_UNIBRANCH_ATTACHMENT_FITTING.md`](LOG_UNIBRANCH_ATTACHMENT_FITTING.md).
> The lower endpoint is attained by an exact Keller-compatible fold whose
> cusp pullback factors as a doubled boundary branch plus a tangent affine
> companion and whose point quotient has length one.  This theorem computes
> local Fitting packets; inserting them into the global logarithmic `ch_2`
> budget still requires the common perfect-complex filtration.

The identities and fiber arithmetic are replayed by
[`verify_log_cusp_attachment_dichotomy.py`](../scripts/verify_log_cusp_attachment_dichotomy.py).

## 1. The general logarithmic order gate

Let `f:X->Y` be a morphism of smooth characteristic-zero surfaces.  At a
point `p` of a reduced SNC source boundary, choose a boundary component

\[
 E=(r=0)
\]

that maps nonconstantly to a reduced unibranch plane curve `C`.  Write the
completed source ring as `R=k[[r,t]]`.  The point is of one of two types:

\[
\begin{array}{c|c|c}
\text{type}&D_X&\epsilon_p\\ \hline
\text{smooth boundary point}&(r=0)&0,\\
\text{boundary node}&(rt=0)&1.
\end{array}                                                    \tag{1.1}
\]

Let `tau` be a normalization parameter of `C`, let

\[
 m_C=\min(\operatorname{ord}_\tau x,
          \operatorname{ord}_\tau y),                         \tag{1.2}
\]

and suppose the residue map on `E` has local index `q_p`, so

\[
 \tau=t^{q_p}\cdot\text{unit}.                                \tag{1.3}
\]

Use `(dlog r,t^(1-epsilon_p)dlog t)` as the source logarithmic basis.  Thus
the second basis element is `dt` for `epsilon_p=0` and `dlog t` for
`epsilon_p=1`.  In target coordinates `(x,y)`, the logarithmic differential
matrix is

\[
 \Theta_p=
 \begin{pmatrix}
  r x_r&r y_r\\
  t^{\epsilon_p}x_t&t^{\epsilon_p}y_t
 \end{pmatrix}.                                                \tag{1.4}
\]

Modulo `r`, the two entries in the second row are the derivatives of the
composed branch normalization.  Characteristic zero and (1.2)--(1.3) give

\[
 \min\bigl(\operatorname{ord}_t x_t(0,t),
           \operatorname{ord}_t y_t(0,t)\bigr)
 =q_p m_C-1.                                                   \tag{1.5}
\]

Every first-row entry is divisible by `r`.  Consequently

\[
 \boxed{
 \operatorname{Fitt}_1(\operatorname{coker}\Theta_p)
 \subseteq (r,t^{N_p}),\qquad
 N_p=q_p m_C-1+\epsilon_p.}                                   \tag{1.6}
\]

If the `Fitt_1` locus is isolated at `p`, its ideal is primary to `(r,t)`,
and (1.6) gives the exact lower bound

\[
 \boxed{
 \ell_p R/\operatorname{Fitt}_1(\operatorname{coker}\Theta_p)
 \ge q_p m_C-1+\epsilon_p.}                                  \tag{1.7}
\]

Equality holds precisely when the Fitting ideal itself is
`(r,t^(N_p))`.  If `Fitt_1` is not primary, the failure is not a smaller
point term: it is a positive-dimensional noncyclic packet that must be kept
separately.

The later cyclic-submodule positivity theorem supplies the missing K-theory
sign for the isolated case: the scalar determinant module injects into the
actual cokernel with an effective finite quotient, and that quotient length
dominates the colength in (1.7).

<!-- status-consumer: LCSP1 8658eebeb1d65671 -->

This bound requires neither transverse unimodularity nor the minimal
Jacobian hypothesis of the two-boundary theorem.  It is only a statement
about the first nonzero tangential derivative and the actual reduced source
log boundary.

## 2. A universal monomial fold

The lower bound is sharp.  For coprime `1<m<n` and residue index `q`, put

\[
 \boxed{
 x=t^{mq}+r,\qquad
 y=t^{nq}+\frac nm r t^{(n-m)q}.}                           \tag{2.1}
\]

On `E=(r=0)` this is the monomial branch `(t^(mq),t^(nq))`.  The linear
term in the pullback of its equation cancels:

\[
 y^m-x^n
 =-\frac{n(n-m)}{2m}
   r^2t^{mq(n-2)}+O(r^3).                                    \tag{2.2}
\]

Thus the boundary factor has transverse multiplicity two.  The ordinary
Jacobian factors exactly as

\[
 \boxed{
 J_f=\frac{n(n-m)q}{m}
     r t^{(n-m)q-1}.}                                        \tag{2.3}
\]

There are now two cases.

1. If `(n-m)q=1`, necessarily `q=1` and `n=m+1`.  The ramification divisor
   is only `E`, so `p` may be a smooth source-boundary point and
   `epsilon=0`.
2. If `(n-m)q>1`, the residual divisor `t=0` occurs in (2.3).  When the map
   is étale off the compactification boundary, that divisor must be another
   boundary component.  Thus `p` is an SNC boundary node and `epsilon=1`.

In either case the column operation

\[
 C_2\longmapsto C_2-\frac nm t^{(n-m)q}C_1                 \tag{2.4}
\]

reduces (1.4) to

\[
 \begin{pmatrix}
 r&0\\
 mq\,t^{mq-1+\epsilon}&
 \dfrac{n(n-m)q}{m}
 r t^{(n-m)q-1+\epsilon}
 \end{pmatrix}.                                               \tag{2.5}
\]

Therefore

\[
 \boxed{
 \operatorname{Fitt}_1=(r,t^{mq-1+\epsilon}),\qquad
 \ell R/\operatorname{Fitt}_1=mq-1+\epsilon.}                \tag{2.6}
\]

The fold family attains the lower bound (1.7) in both boundary-incidence
types.  The order shift by one is not an arbitrary convention: it records
whether `dt` or `dlog t` belongs to the source log basis.

## 3. The ordinary-cusp fold has point charge one

The case `(m,n,q)=(2,3,1)` is the decisive local model:

\[
 x=t^2+r,\qquad
 y=t^3+\frac32rt.                                             \tag{3.1}
\]

It has exact factorizations

\[
 \boxed{
 y^2-x^3=-r^2\left(r+\frac34t^2\right),\qquad
 J_f=\frac32r.}                                               \tag{3.2}
\]

Thus:

- `E=(r=0)` occurs twice in the cusp pullback, so its transverse index is
  `e=2`;
- the companion `r+3t^2/4=0` is an affine pullback branch, not a second
  boundary component;
- the companion meets `E` with contact two; and
- the surface map is étale away from `E` in this completed neighborhood.

With the actual source boundary `D_X=E`, the logarithmic matrix is

\[
 \Theta_{\mathrm{fold}}=
 \begin{pmatrix}
 r&\frac32rt\\
 2t&3t^2+\frac32r
 \end{pmatrix}
 \sim
 \begin{pmatrix}r&0\\t&r\end{pmatrix}.                       \tag{3.3}
\]

Hence

\[
 \det\Theta_{\mathrm{fold}}\sim r^2,
 \qquad \operatorname{Fitt}_1=(r,t).                         \tag{3.4}
\]

If `M_fold=coker(Theta_fold)`, the first target generator gives the exact
sequence

\[
 \boxed{
 0\longrightarrow R/(r^2)
 \longrightarrow M_{\mathrm{fold}}
 \longrightarrow R/(r,t)
 \longrightarrow0.}                                         \tag{3.5}
\]

The scalar determinant packet `R/(r^2)` therefore misses one positive point
class.  This is an exact local logarithmic correction of length one.  It is
also a concrete reason why the transverse divisorial index `e=2` must not be
multiplied blindly into the node formula `q_p*m_C`.

## 4. Complete-fiber ledger

Let a boundary normalization cover the target normalization with residue
degree `f`.  Over one unibranch value, write the local residue indices as
`q_1,...,q_h`, so

\[
 \sum_{p=1}^h q_p=f.                                          \tag{4.1}
\]

Let `c` of these `h` points be source-boundary nodes.  Summing (1.7) gives

\[
 \boxed{
 B_C\ge\sum_p(q_p m_C-1+\epsilon_p)
       =m_C f-h+c.}                                           \tag{4.2}
\]

Since `h<=f`, every isolated complete fiber obeys

\[
 \boxed{B_C\ge(m_C-1)f.}                                    \tag{4.3}
\]

For an ordinary cusp,

\[
 \boxed{B_{A_2}\ge2f-h+c,\qquad f\le2f-h+c\le2f.}          \tag{4.4}
\]

The two exact endpoints have different geometry:

- `h=f,c=0`: all points are unramified smooth-boundary folds, and the exact
  fold model gives total point correction `f`;
- `c=h`: every point is an SNC boundary attachment, and the exact minimal
  packet of `LUAF1` gives total `2f`.

Intermediate residue partitions and boundary incidences give intermediate
ledgers.  Therefore cusp multiplicity and residue degree alone do not select
one number.  The source compiler must additionally record whether each cusp
preimage is a smooth boundary point or a boundary node.

## 5. Consequence for the F2 `k=1` cusp face

On the generic nonimmersion face the target packet is one ordinary cusp plus
three nodes.  The affine strict-log-etale theorem still removes all affine
preimages from the relative logarithmic ledger.  At the compactification
boundary, however, there are now two exact minimal possibilities.

If the common `ch_2` filtration identifies the local Fitting quotients with
its finite point pieces, the cusp subledger is

\[
 B_{\mathrm{cusp}}=2f-h+c                                  \tag{5.1}
\]

for a fiber assembled from the exact minimal packets above.  Thus the
previous subtraction `2f` is valid only at the node-saturated endpoint.  At
the all-unramified smooth-fold endpoint it must be replaced by `f`.

For the minimal arithmetic signature `(e,f,n)=(2,1,1)`, the doubled virtual
point numerators before the cusp packet are

\[
 12-4b-s_X,\qquad49-4b-s_X.                                \tag{5.2}
\]

At the smallest compatible values `s_X=0/1`, the unidentified residuals are

\[
\begin{array}{c|c|c}
&\text{smooth fold, }B=1&\text{boundary node, }B=2\\ \hline
\text{squarefree}&5-2b&4-2b\\
\text{double}&23-2b&22-2b.
\end{array}                                                   \tag{5.3}
\]

This does not select the actual F2 packet.  It does close a logical gap: a
target cusp does not by itself force the stronger `2f` subtraction.  The
next source calculation is now sharply finite—locate the `h` points, compute
their residue indices, and mark the `c` points that are actual boundary
nodes.  Higher local expansion is needed only where equality in (1.7) fails.

## Reproduction

```bash
.venv/bin/python scripts/verify_log_cusp_attachment_dichotomy.py
```
