# Posterior orientation and the bisector of two Gaussian triple points

Fix \(0<r<1\). Let \(A\) be an ordinary balanced planar
Gaussian noise-stability maximizer. Put
\(U_i=T_r1_{A_i}\), write \(V_i=U_i-\beta_i\) for the
adjusted scores with \(\sum_i\beta_i=0\), and use their strict
winner sets \(O_i=\{V_i>V_j\text{ for }j\ne i\}\).
The strict winners agree with the partition up to null sets.
A maximal triple point is a point where all three adjusted
scores are equal. Call it transverse if the differential of
\[
                     F=(V_1-V_3,V_2-V_3)
\]
is invertible there. The sign of its determinant specifies
one of the two orientations in a fixed planar coordinate system.
[Fixed-mass assignment and winner equations](MAXIMIZER_PERTURBATIVE_RIGIDITY.md#assignment-and-the-exact-information-supplied-by-global-maximality)

Three restrictions follow. A bounded strict-winner component
has zero signed count of regular preimages of each score value
in its winner region. Two triple points force every label
to cross their contracted bisector. If they have equal radii,
their centered bisector reflection itself must realize all
three unordered mismatch types.

## Degree zero inside each bounded winner component

Let \(C\) be a bounded connected component of \(O_i\).
Let \(W_i\) be its open winner wedge in the plane of two-score
differences. For example \(W_1=\{(u,v):u>\max(v,0)\}\).
Then \(F(C)\subset W_i\) and \(F(\partial C)\subset\partial W_i\).
Indeed continuity gives the corresponding nonstrict winner
inequalities at a point of \(\partial C\). If all those
inequalities were strict there, a ball in \(O_i\) would meet
\(C\), making the point interior to \(C\).

At every value \(w\in W_i\) regular for \(F|_C\),
\[
              \sum_{\substack{x\in C\\F(x)=w}}
                       \operatorname{sgn}\det DF(x)=0.     \tag{DB1}
\]
Here regular means that the differential is invertible at
every preimage in \(C\). The preimages are finite by
compactness of \(\overline C\), the boundary exclusion, and
local invertibility. The assertion uses only the strict-winner
representation and smoothness of the score map.

For a direct degree proof without regularity assumptions on
\(\partial C\), take any \(h\in C_c^\infty(W_i)\) and extend
it by zero to the plane. Define the one-form
\[
 \alpha(w)=\left(\int_0^1 t h(tw)\,dt\right)
                   (w_1\,dw_2-w_2\,dw_1).
\]
All rays through the angular directions of the support stay
strictly inside \(W_i\), and \(h\) is supported away from the
origin. Thus \(\alpha\) vanishes in a neighborhood of
\(F(\partial C)\), including the origin if it occurs.
Differentiation under the integral gives
\(d\alpha=h(w)\,dw_1\wedge dw_2\), by integrating the derivative
of \(t^2h(tw)\) in \(t\). The pullback \(F^*\alpha\) on \(C\)
extends by zero to a smooth compactly supported one-form,
since \(F^*\alpha=0\) near the compact boundary. Its exterior
derivative integrates to zero. Therefore
\[
                        \int_C h(F(x))\det DF(x)\,dx=0.
\]
For \(h\) supported sufficiently close to a regular value \(w\),
every preimage in \(C\) lies in one of the inverse-function
charts of the finite fiber over \(w\). This follows by
compactness and the boundary exclusion. Changing variables
on these charts gives (DB1).

If every triple point in \(\overline C\) is transverse, it
follows that \(C\) has equally many incident triple points
of the two orientations. In fact the triple points in this
compact set are finite, since \(F^{-1}(0)\) is closed and each
of these points has a locally invertible score map. Near an
incident triple \(p\) the inverse chart identifies \(O_i\) with
the local connected wedge \(W_i\). This sector belongs to
\(C\) since \(p\in\overline C\). Thus each incident triple
supplies one preimage in \(C\) of every sufficiently small
\(w\in W_i\). There are no additional preimages for small
\(w\) by compactness and the list of zeros of \(F\) in
\(\overline C\). Apply (DB1) and use the determinant sign in
each inverse chart. In particular a bounded component with
one incident transverse triple cannot have that triple as
its only incident triple under these hypotheses.

## Two triple points and the contracted bisector

Let \(p,q\) be distinct maximal triple points. Set
\[
 m=(p+q)/2,\qquad e=(p-q)/|p-q|,\qquad
 \sigma^2=1-r^2,\qquad d=r|p-q|/2,\qquad \tau=d/\sigma^2,
\]
and let \(\mu\) be the Gaussian with mean \(rm\) and
covariance \(\sigma^2 I\). For \(y\) in the source plane
write \(z=(y-rm)\cdot e\). Then each label satisfies
\[
                      \int_{A_i}\sinh(\tau z)\,d\mu(y)=0. \tag{DB2}
\]
Indeed \(\sum_i V_i=1\), so \(U_i(p)=U_i(q)=1/3+\beta_i\).
The conditional source measures at \(p,q\) have densities
\(e^{-d^2/(2\sigma^2)}e^{\tau z}\) and
\(e^{-d^2/(2\sigma^2)}e^{-\tau z}\), respectively, relative
to \(\mu\). Subtract their probabilities of \(A_i\).
Since each label has positive area and the bisecting line
is null, (DB2) forces positive area of each label on both
sides of \(y\cdot e=rm\cdot e\).

This straddling has a uniform positive margin when \(r\) lies
in a compact positive band and \(p,q\) lie in a fixed disk.
In fact let \(\theta_*>0\) be the score lower bound (PP1)
on the band at winning points. The conditional density at
\(p\) relative to \(\mu\) has squared \(L^2(\mu)\) norm
\(e^{d^2/\sigma^2}\). Hence Cauchy--Schwarz gives
\[
              \mu(A_i)\ge\theta_*^2e^{-d^2/\sigma^2}
                                      \ge c_0>0
\]
uniformly. The centered scalar \(z\) has variance \(\sigma^2\)
bounded away from zero. Choose \(a>0\) uniformly small
so that \(\mu(|z|\le a)<c_0/2\), by the uniform density
bound. Thus \(\int_{A_i}|z|\,d\mu\ge a c_0/2\).
Equation (DB2) and \(|\sinh(\tau z)|/\tau\ge|z|\) give
an integral of at least \(a c_0/4\) for
\(|\sinh(\tau z)|/\tau\) on each half of \(A_i\).
The integrands have a uniform \(L^2(\mu)\) bound, since
\[
              |\sinh(\tau z)|/\tau\le |z|e^{\tau|z|}
\]
and \(\tau\) is bounded on the stated parameter sets.
Another Cauchy--Schwarz application lower-bounds the
\(\mu\)-mass of each half by a uniform constant. Removing
a sufficiently small strip about \(z=0\) retains a uniform
positive fraction of this mass by the same density bound.
In particular there exist \(\kappa,c>0\) independent even
of the separation of \(p,q\) such that
\[
        \mu(A_i\cap\{z\ge\kappa\})\ge c,\qquad
        \mu(A_i\cap\{z\le-\kappa\})\ge c
                     \quad(i=1,2,3).                     \tag{DB3}
\]
[Uniform winner score lower bound, PP1](THREE_CELL_PRESSURE_CONTRACTION_PROPAGATION.md#a-positive-own-label-score-on-its-winning-set)

For example, if a bounded face of a compact theta belongs
to one label and is that label's entire cell, and the face
lies in the strip between the chord projections of \(p,q\),
then (DB3) enforces
\[
                   |(1-r)m\cdot e|
                               \le |p-q|/2-\kappa .
\]
This consequence requires the projection-strip containment.
Neither (DB2) nor (DB3) needs that containment, a theta
topology, or transversality of the triples.

## The targeted reflection at equal vertex radii

Suppose \(|p|=|q|\). The reflection \(R\) through the line
orthogonal to \(p-q\) and through the origin exchanges
\(p,q\). Let \(H\) be its open half-plane containing \(p\)
and put \(\delta_i(y)=1_{A_i}(y)-1_{A_i}(Ry)\) there.
If at most two unordered mismatch types occur on positive
measure in \(H\), the strict two-type orientation comparison
(DG12) makes every \(\delta_i\) of one weak sign almost
everywhere on \(H\). To see the sign implication, orient
the one or two active types away from a common label.
The maximizing signs are then all equal, so each label's
coordinate has a fixed sign on its support.
[Exact two-type reflection comparison, DG12](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md#exact-comparisons-within-reflection-orbits)

On the other hand the equality of all unadjusted scores at
the two triples gives
\[
 0=U_i(p)-U_i(Rp)
   =\int_H [M_r(p,y)-M_r(p,Ry)]\delta_i(y)\,d\gamma_2(y).
\]
The bracket is strictly positive for \(y\in H\), since the
normal coordinates of \(p,y\) have the same positive sign.
Each \(\delta_i\) would therefore vanish almost everywhere.
The three cells would all be fixed by \(R\) up to null sets,
so their Gaussian centroids would be in the reflecting
line. This contradicts full active-plane centroid rank.
Hence every one of the three mismatch types has positive
measure for this prescribed reflection.
[Full active-plane centroid rank](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-active-dimension-is-impossible)

The quantifier here differs from (DG13). That result supplies
some reflection with three types for every nonsimplex actual
maximizer. Here the two triple locations prescribe the
reflection, and their shared conditional score vector forces
the three types there. The argument does not conclude that
equal-radius triples occur or that their presence is impossible.
[Existential reflection-triangle witness, DG13](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md#reflection-order-and-a-necessary-three-type-witness)
