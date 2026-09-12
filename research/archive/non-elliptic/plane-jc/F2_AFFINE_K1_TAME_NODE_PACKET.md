# F2 `k=1` tame-node logarithmic packet theorem

> **Status.**  A target node does not by itself force a logarithmic
> `Fitt_1` or localized-`ch_2` correction.  Every fine-and-saturated tame
> Kummer toroidal pullback of the two target branches is log-étale in characteristic
> zero.  This includes both a split smooth source point
> `(x,y)=(r^e,t)` and the collided cyclic model `z^e=xy`; after the minimal
> toric resolution their logarithmic cotangent cokernels vanish.  Therefore
> the degree-eight target conductor and its affine normalization defects
> cannot be inserted as positive boundary-Chern lengths.  Any nonzero
> boundary point correction must come from non-toroidal branch gluing,
> residue ramification, or a genuinely nonunit `Fitt_1` matrix.  More
> generally, an arbitrary SNC monomial-with-unit pullback has zero cokernel
> whenever its exponent matrix has full rank.  In rank one its cokernel is
> cyclic, and two explicit first-jet equations are necessary for a singular
> determinant support and hence for a positive normalization mismatch.

The lattice matrices, the full `A_(e-1)` resolution fan, and the vanishing
claim are replayed by
[`verify_f2_affine_k1_tame_node_packet.py`](../scripts/verify_f2_affine_k1_tame_node_packet.py).

## 1. The split smooth packet

Let the completed target surface at a node have parameters `(x,y)`, so the
two completed branches of the target curve are the coordinate axes.  At a
source point where one ramified boundary branch and one companion branch
remain separated, the saturated tame model is

\[
 x=r^e,\qquad y=t,qquad e\ge2.                  \tag{1.1}
\]

Using the full reduced pullback divisor `rt=0`, the logarithmic differential
matrix is

\[
 \begin{pmatrix}e&0\\0&1\end{pmatrix}.           \tag{1.2}
\]

In characteristic zero, `e` is a unit.  Thus (1.2) is invertible,

\[
 \operatorname{coker}\theta^{\log}=0,
 \qquad \operatorname{Fitt}_1=\mathcal O,
 \qquad ch_2^{\rm loc}=0.                       \tag{1.3}
\]

The ordinary curve normalization quotient at the node still has delta
length one.  Equation (1.3) says that this ordinary conductor length is not
an additional logarithmic differential defect.

## 2. The collided cyclic packet

The opposite extreme is the normal cyclic cover

\[
 X_e=\operatorname{Spec}
 k[[x,y,z]]/(z^e-xy)\longrightarrow
 \operatorname{Spec}k[[x,y]].                   \tag{2.1}
\]

This is the toric `A_(e-1)` surface.  In additive character notation, let

\[
 M_Y=\mathbb Zx\mathbin\oplus\mathbb Zy,
 \qquad M_X=M_Y+\mathbb Zz,
 \qquad ez=x+y.                                    \tag{2.2}
\]

Thus `M_X/M_Y` is cyclic of order `e`.  In the source-lattice basis put

\[
 m_1=x,\qquad m_2=z=\frac{x+y}{e}.
\]

Then

\[
 x=m_1,\qquad y=e m_2-m_1,                      \tag{2.3}
\]

so the logarithmic character matrix is

\[
 A_e=\begin{pmatrix}1&-1\\0&e\end{pmatrix},
 \qquad \det A_e=e.                             \tag{2.4}
\]

Its Smith form over the integers is `diag(1,e)`.  Thus the integral group
cokernel is exactly `Z/e`, while `A_e` becomes invertible over the
characteristic-zero ground field and its entries generate the unit ideal.

In the dual lattice, the cone is spanned by `(0,1)` and `(e,1)`.  Its minimal
regular subdivision has rays

\[
 v_j=(j,1),\qquad 0\le j\le e.                  \tag{2.5}
\]

Adjacent determinants are one in absolute value.  Every interior ray obeys

\[
 \operatorname{ord}_{v_j}(x,y)=(j,e-j),          \tag{2.6}
\]

so it is contracted to the target node, and

\[
 v_{j-1}+v_{j+1}=2v_j                            \tag{2.7}
\]

shows that the `e-1` exceptional curves form the usual `(-2)` chain.
Toric subdivision is itself log-étale.  Hence the resolved map still has

\[
 \boxed{
 \mathcal T_f^{\log}=0,
 \quad \operatorname{Fitt}_1=\mathcal O,
 \quad ch_2^{\rm loc}=0.}                       \tag{2.8}
\]

The exceptional chain changes the ordinary boundary graph but contributes
no logarithmic cokernel term.

## 3. General fs tame Kummer toroidal packet

Both examples are instances of a rank-two Kummer injection between fine
and saturated characteristic monoids.  Here **toroidal packet** includes the
essential chart hypothesis: after completion, the source is étale over the
toric base change defined by the monoid map.  On group completions there is
an inclusion

\[
 M_Y\hookrightarrow M_X                            \tag{3.1}
\]

with finite cokernel of order prime to the characteristic.  Logarithmic
differentials tensor this inclusion with the ground field.  The finite
cokernel is killed, so

\[
 M_Y\otimes k\mathrel{\cong}M_X\otimes k.        \tag{3.2}
\]

Kato's chart criterion now applies: the comparison morphism is étale by the
toroidal-packet hypothesis, and the group cokernel has invertible order.
The packet is therefore log-étale.  Every regular toric subdivision is a log
blowup and remains log-étale.  This proves the vanishing in (2.8) for every
fs tame Kummer toroidal node packet, independently of the Kummer index and
the chosen regular fan.

The toroidal chart hypothesis cannot be omitted.  A completed source map
with the same boundary exponents may contain transverse or branch-gluing
terms not étale over the toric base change; those terms are precisely where
a nonzero logarithmic cokernel can still occur.

## 4. The general SNC rank and first-jet gate

There is a useful calculation that does not assume the toroidal chart
hypothesis.  Let `(u,v)` be regular parameters on a completed smooth source
with reduced SNC divisor `uv=0`, and suppose

\[
 x=u^a v^b\alpha(u,v),\qquad
 y=u^c v^d\beta(u,v),                            \tag{4.1}
\]

where `alpha,beta` are units.  In the target basis `(dlog x,dlog y)` and
source basis `(dlog u,dlog v)`, the logarithmic matrix is

\[
 \Theta=
 \begin{pmatrix}
 a+u\partial_u\log\alpha&c+u\partial_u\log\beta\\
 b+v\partial_v\log\alpha&d+v\partial_v\log\beta
 \end{pmatrix}.                                  \tag{4.2}
\]

Its residue matrix and determinant are

\[
 A=\begin{pmatrix}a&c\\b&d\end{pmatrix},
 \qquad \det\Theta(0,0)=ad-bc.                  \tag{4.3}
\]

This gives an immediate rank gate.

- If `ad-bc` is nonzero, then `det Theta` is a unit and the logarithmic
  cokernel vanishes.  No higher unit jet can change this conclusion; in
  particular, full-rank exponent data are enough even when the map has not
  been put in an exact toric chart.
- If `A` has rank one, some exponent is a nonzero characteristic-zero
  scalar.  Hence `Fitt_1(coker Theta)=R`; provided `det Theta` is not
  identically zero, the cokernel is cyclic, `R/(det Theta)`.

Write

\[
 A_u=(\partial_u\log\alpha)(0),\quad
 A_v=(\partial_v\log\alpha)(0),\quad
 B_u=(\partial_u\log\beta)(0),\quad
 B_v=(\partial_v\log\beta)(0).
\]

Under the rank-one equation `ad-bc=0`, direct expansion gives

\[
 \det\Theta\equiv
 u(dA_u-bB_u)+v(aB_v-cA_v)\pmod{(u,v)^2}.        \tag{4.4}
\]

If either displayed coefficient is nonzero, the determinant support is
smooth at the point.  The cyclic module is then already a rank-one module
on its normalization, so its conductor-gluing quotient and its
codimension-two normalization/Smith defect vanish.  A positive point defect
can occur only behind the simultaneous gate

\[
 \boxed{dA_u-bB_u=0,\qquad aB_v-cA_v=0.}         \tag{4.5}
\]

These equations are necessary, not sufficient: once they hold, the
quadratic and higher unit jets of `det Theta` and the full matrix still have
to be calculated.  A nonunit `Fitt_1` is even more restrictive.  It cannot
occur in the monomial-with-unit form (4.1) at a point over `x=y=0`, because
at least one exponent is nonzero.  Such a locus certifies that the chosen
log divisor or completed monomial chart is incomplete.

Thus the next source-side input is finite: at each unresolved attachment,
extract `(a,b,c,d)` and the four logarithmic unit derivatives in (4.4).
Full rank ends the packet immediately; rank one with a nonzero first jet
ends its point-defect calculation; only rank one satisfying both equations
in (4.5) requires a higher-order local expansion.

## 5. Consequence for the F2 `k=1` budget

The all-stratum conductor theorem proves that the target affine
normalization quotient has total length four and its conductor divisor has
degree eight.  Those are invariants of the pulled-back curve.  They are not,
by themselves, lengths of the logarithmic cotangent cokernel.

This theorem closes the automatic-node-correction route:

1. a split node packet contributes zero;
2. collision of the two Kummer branches and the entire resulting toric
   resolution chain still contribute zero; and
3. no positive point term may be charged merely because the target has a
   node or a nonzero conductor divisor; the same applies to a resolved cusp
   only after its source packet satisfies the toroidal chart hypothesis.

Therefore a nonzero boundary-local term in the `k=1` Chern sieve must be
certified by additional data: a non-toroidal completed pullback matrix,
ramification of the residue map along a conductor point, collision of
distinct source packets, or a nonunit `Fitt_1` ideal.  The complete source
Laurent pair is still needed to decide whether any of these occurs.  This
theorem does not determine `(e,f,E^2)`, exclude `(75,125)`, or prove
`JC(2)`.

The complementary
[`unibranch attachment theorem`](LOG_UNIBRANCH_ATTACHMENT_FITTING.md)
computes one such nonunit `Fitt_1` packet.  When a boundary component maps to
a plane branch of multiplicity `m_C` with local residue index `q_p` and meets the
forced residual ramification branch minimally and transversely, its log
matrix is `diag(r,t^(q_p*m_C))`.  The point correction is `q_p*m_C`; over a
complete residue-degree-`f` ordinary-cusp fiber the total is `2f`.  This is
conditional on the attachment chart and
does not charge affine cusp preimages automatically.

<!-- status-consumer: LUAF1 b0279670ffbd3fa5 -->

The general input is Kato's chart criterion for logarithmic étaleness
([Theorem 3.5](https://www.math.brown.edu/dabramov/LOGGEOM/Kato-log.pdf));
the preservation under subdivision uses the fact that log blowups are
log-étale, as in F. Kato's treatment of
[log modifications](https://arxiv.org/abs/math/9907124).  The calculations
above identify the relevant charts and do not replace those general results.

## Reproduction

```bash
.venv/bin/python scripts/verify_f2_affine_k1_tame_node_packet.py
```
