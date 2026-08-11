# Unibranch boundary attachment and the logarithmic `Fitt_1` point class

> **Status.**  Exact local theorem under an explicit minimal-transverse SNC
> attachment hypothesis.  If a source boundary component maps to a target
> plane branch of multiplicity `m_C` with residue ramification index `q_p`, and
> the forced residual ramification branch is a second reduced SNC boundary
> component, then the completed logarithmic cotangent cokernel is
> `R/(r) direct-sum R/(t^(q_p*m_C))`.  Its scalar determinant module is
> `R/(r*t^(q_p*m_C))`, and the normalization/Smith discrepancy is the point
> module `R/(r,t^(q_p*m_C))` of length `q_p*m_C`.  For an ordinary cusp this is
> `2q_p`, not the corank one of the carrier-jet parameter map.  The theorem does
> not prove that a particular unresolved F2 attachment satisfies the stated
> transversality or SNC hypotheses.

The normal forms, Fitting ideals, and point lengths are replayed by
[`verify_log_unibranch_attachment_fitting.py`](../scripts/verify_log_unibranch_attachment_fitting.py).

## 1. The local geometric setup

Let `f:X->Y` be a morphism of smooth surfaces over a characteristic-zero
field.  Suppose `p` lies on two reduced SNC source-boundary components

\[
 E=(r=0),\qquad F=(t=0)                         \tag{1.1}
\]

in the complete regular local ring

\[
 R=k[[r,t]],                                    \tag{1.2}
\]

and suppose `f(p)` lies in the affine interior of `Y`, so there is no target
logarithmic divisor near `f(p)`.

Assume that `E` maps nonconstantly to a reduced unibranch plane curve `C`.
Let `tau` be a normalization parameter of `C`, let

\[
 m_C=\min(\operatorname{ord}_\tau x,
          \operatorname{ord}_\tau y)            \tag{1.3}
\]

be the branch multiplicity in regular target parameters `(x,y)`, and let the
residue map on `E` have ramification index `q_p`:

\[
 \tau=t^{q_p}\cdot\text{unit}.                  \tag{1.4}
\]

Then the tangential derivative pair has common order

\[
 \rho=q_p m_C-1.                                \tag{1.5}
\]

More invariantly, (1.5) says that in `R/(r)=k[[t]]`,

\[
 (\partial_t x,\partial_t y)=(t^\rho).          \tag{1.6}
\]

The theorem needs two minimality hypotheses at `p`:

1. the transverse derivative pair `(partial_r x,partial_r y)` is unimodular;
2. the ordinary Jacobian has exactly the forced tangential order,

\[
 \det\frac{\partial(x,y)}{\partial(r,t)}
 =t^\rho\cdot\text{unit}.                       \tag{1.7}
\]

Condition (1.7) makes `F` the reduced residual ramification branch.  If the
order is larger or the residual divisor is nonreduced/non-SNC, the same
calculation gives a lower bound but not the diagonal normal form below.

## 2. Why a singular branch forces an attachment

If `m_C>1`, the normalization map of `C` is nonimmersive at the marked point,
so `rho>=1`.  Along `E`, (1.7) forces the ordinary Jacobian to vanish at `p`.
Its zero set is divisorial.  When `f` is étale off the source boundary, that
residual divisor must also lie in the boundary.  Thus a nonimmersive image
point cannot be represented by an isolated rank drop on a lone boundary
component: a second boundary branch (or a non-SNC/nonreduced replacement)
must pass through it.

This is the geometric distinction between two statements that otherwise
look numerically similar:

- the raw carrier-jet map ramifies in the **parameter space** of target
  curves; and
- the logarithmic cotangent matrix degenerates on the **source boundary**.

Only the second contributes to the localized Chern complex.

## 3. Exact logarithmic matrix normal form

Use `(dlog r,dlog t)` as the basis of
`Omega_X^1(log(E+F))`.  The map from ordinary target differentials is
represented by

\[
 \Theta=
 \begin{pmatrix}
 r x_r&r y_r\\
 t x_t&t y_t
 \end{pmatrix}.                                 \tag{3.1}
\]

The transverse unimodularity permits an invertible target-column operation
making the first row `(r,0)`.  After that operation write

\[
 \Theta\sim
 \begin{pmatrix}r&0\\h&g\end{pmatrix}.          \tag{3.2}
\]

Equations (1.6) and (1.7) imply

\[
 g=t^{\rho+1}\cdot\text{unit},\qquad
 h\in(r,t^{\rho+1}).                            \tag{3.3}
\]

Row and column operations first remove the `r`-part and then the
`t^(rho+1)`-part of `h`.  Hence

\[
 \boxed{
 \Theta\sim\operatorname{diag}(r,t^{\rho+1})
 =\operatorname{diag}(r,t^{q_p m_C}).}          \tag{3.4}
\]

Consequently

\[
\boxed{\begin{aligned}
T_f^{\log}&=\operatorname{coker}\Theta
 \cong R/(r)\oplus R/(t^{q_p m_C}),\\
\operatorname{Fitt}_0(T_f^{\log})&=(r t^{q_p m_C}),\\
\operatorname{Fitt}_1(T_f^{\log})&=(r,t^{q_p m_C}).
\end{aligned}}                                  \tag{3.5}
\]

The first line contains the gluing information that the scalar determinant
forgets.

## 4. The exact point correction

For arbitrary positive `N`, the standard intersection/sum sequence for the
ideals `(r)` and `(t^N)` is

\[
0\longrightarrow R/(rt^N)
\longrightarrow R/(r)\oplus R/(t^N)
\longrightarrow R/(r,t^N)
\longrightarrow0.                              \tag{4.1}
\]

The last module has basis

\[
 1,t,\ldots,t^{N-1}                             \tag{4.2}
\]

over `k`, and therefore length `N`.  Taking `N=q_p*m_C` gives

\[
 \boxed{
 [T_f^{\log}]=[R/(r t^{q_p m_C})]
              +q_p m_C[k(p)].}                 \tag{4.3}
\]

Thus after the generic divisorial determinant/Smith packet has been booked,
the point-supported correction is positive and equals

\[
 \boxed{\ell_p\bigl(R/\operatorname{Fitt}_1(T_f^{\log})\bigr)
 =q_p m_C.}                                     \tag{4.4}
\]

This is exactly the generalized version of

\[
0\to R/(uv)\to R/(u)\oplus R/(v)\to R/(u,v)\to0,
\]

with the second branch thickened to order `q_p*m_C`.

## 5. Model and the ordinary-cusp specialization

The monomial-transverse model

\[
 x=t^m,qquad y=t^n+r,qquad 1\le m<n           \tag{5.1}
\]

has boundary `rt=0` and logarithmic matrix

\[
 \Theta_{m,n}=
 \begin{pmatrix}0&r\\m t^m&n t^n\end{pmatrix}. \tag{5.2}
\]

The column operation

\[
 C_2\longmapsto C_2-\frac nm t^{n-m}C_1         \tag{5.3}
\]

followed by a column swap gives `diag(r,m*t^m)`.  Since `m` is a unit,
(3.4)--(4.4) hold with point length `m`.

For the ordinary cusp `(m,n)=(2,3)` and residue ramification `q_p`, the
composed normalization has minimum order `2q_p`.  Therefore a minimal
transverse SNC boundary attachment has

\[
 \boxed{
 T_f^{\log}\cong R/(r)\oplus R/(t^{2q_p}),
 \qquad \ell_p R/\operatorname{Fitt}_1=2q_p.}  \tag{5.4}
\]

For the generic F2 `k=1` nonimmersion stratum, the raw seven-jet map has
corank one, but (5.4) gives point length `2q_p`.  The two numbers measure
different maps and must not be identified.  When `q_p=1`, the point length two
happens to equal the cusp conductor exponent; that equality is special to
the ordinary cusp and is not the general formula for a plane branch.

## 6. Consequence for the global `ch_2` budget

This theorem turns one formerly qualitative gap into a finite local test.
For every boundary point above a nonimmersive target-curve value, the source
compiler should record:

1. the residue ramification index `q_p` on the boundary normalization;
2. the target branch multiplicity `m_C`;
3. whether the transverse derivative pair is unimodular;
4. whether the residual Jacobian branch is reduced SNC and has the minimal
   order (1.7).

If all four hold, the point correction is exactly `q_p*m_C`; no higher Laurent
descent is needed.  If they fail, the failure itself identifies the required
non-toroidal or higher-order packet.

For `(75,125)`, the generic `k=1` cusp stratum would therefore contribute
`2q_p` at each minimal boundary attachment above its cusp value.  If the
boundary normalization covers the target normalization with residue degree
`f` and every point over the cusp is minimal transverse SNC, then the degree
formula for a finite map of smooth curves gives

\[
 \boxed{\sum_{p\mapsto\mathrm{cusp}}2q_p=2f.}   \tag{6.1}
\]

More generally, a complete minimal fiber over a unibranch value contributes
`m_C*f`.  This weight is branchwise and produces the following distinct F2
patterns:

- an ordinary cusp has `m_C=2` and contributes `2f`;
- the monomial `(3,5)` cusp has `m_C=3` and contributes `3f`;
- an ordinary node has two immersive smooth branches and forces no residual
  attachment by this theorem; and
- the three branches of an ordinary triple point are likewise individually
  immersive, so their total point multiplicity three must not be substituted
  for `m_C`.

Thus neither delta invariant nor conductor degree determines the logarithmic
point charge.  They can remain constant while the branch-derivative
attachment ledger changes.

The two-boundary hypothesis is essential.  The complementary
[`cusp attachment dichotomy`](LOG_CUSP_ATTACHMENT_DICHOTOMY.md) proves the
general isolated-Fitting lower exponent

\[
 q_p m_C-1+\epsilon_p,
\]

where `epsilon_p=0` at a smooth source-boundary point and `epsilon_p=1` at
an SNC boundary node.  In particular the exact ordinary-cusp fold

\[
 x=t^2+r,\qquad y=t^3+\frac32rt
\]

has transverse index two, only one source-boundary component, and point
correction one.  Hence the total `2f` in (6.1) is the node-saturated endpoint,
not a consequence of the target cusp and residue degree alone.

<!-- status-consumer: LCAD1 7b9c15d3dfae0337 -->

This remains
conditional because the complete fixed-coordinate source pair has not yet
located those points, determined their residue indices, or proved the SNC
transversality.  Affine source preimages of the cusp contribute nothing:
`K_f` is already exact on `A^2`.  Only boundary attachments are counted.

## Reproduction

```bash
.venv/bin/python scripts/verify_log_unibranch_attachment_fitting.py
```

The checker verifies the ordinary-versus-logarithmic order shift, diagonal
reduction, Fitting ideals, ideal intersection/sum sequence, and point length
for model multiplicities and residue indices through a bounded exact range.
