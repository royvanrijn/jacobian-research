# A quadratic dual bound for central radial mass

The undamped central-profile comparison uses
\(\ell(p)=\min(1/4,p_1,p_2,p_3)\) and
\(e_0=\int \ell(p)\,d\mu_1\).
The mass lower bound below keeps the squared residual budget and the
exact least coordinate of the clipped reference. It applies at each
amplitude and label-direction bin of the variance bootstrap.
[Undamped central variance](THREE_CELL_UNDAMPED_CENTRAL_VARIANCE.md)

## The pointwise minimum

Let \(\Delta_3\) be the probability simplex. For \(q\in\Delta_3\)
write \(m=\min_i q_i\). Fix \(w\ge0\) and \(\lambda>0\), and set
\(\delta=\min(m,w/(3\lambda))\). Then
\[
 \begin{split}
 H_\lambda(m,w)
 &:=\inf_{y\in\Delta_3}
       \bigl[w\ell(y)+\lambda|y-q|^2\bigr]\\
 &=\min\left\{\frac w4,\,
             w(m-\delta)+\frac{3\lambda}{2}\delta^2\right\}.
                                                               \tag{CM1}
 \end{split}
\]
Indeed the minimum of the four affine branches in \(\ell\) can be
taken after minimizing each branch. The constant branch gives \(w/4\).
In branch \(i\) it suffices to take \(y_i=q_i-d\) with
\(0\le d\le q_i\), since any increase above \(q_i\) costs at least
\(wq_i\). At a fixed decrease \(d\) the other coordinates both increase
by \(d/2\) in the least-distance choice. This choice is feasible and
has squared distance \(3d^2/2\). The minimizing decrease is
\(\min(q_i,w/(3\lambda))\). The resulting branch value is increasing
in \(q_i\): it is quadratic up to \(w/(3\lambda)\) and affine with
slope \(w\) thereafter. Thus the least reference coordinate attains
the smallest affine-branch value.

The resulting \(H_\lambda\) is nonnegative and increasing in \(m\)
and \(w\). The monotonicity in \(w\) also follows from the pointwise
infimum, since \(\ell(y)\ge0\). These facts allow lower panel enclosures
of both arguments when integrating.

## Integrating against a residual budget

Write \(dJ=J'_m\,d\nu\) for the weighted measure of a spectral slab,
and put \(W=d\mu_1/dJ\). In the usual trial decomposition take
\(\bar p=p_*+A v\phi\), \(q=\Pi_{\Delta_3}(\bar p)\), and \(k=p-q\).
If \(K^2\) bounds \(\int |k|^2dJ\), then (CM1) gives
\[
       e_0\ge\int H_\lambda(\min_i q_i,W)\,dJ-\lambda K^2 .
                                                               \tag{CM2}
\]
This bound holds for every \(\lambda>0\), and taking a maximum over
a finite set of positive choices preserves it. The nonnegative
integrand can be omitted outside any set of radial panels. The
budget follows from the projection inequality
\(\|k\|_J^2\le G_0^2-\|q-\bar p\|_J^2\).

Optimizing (CM2) over \(\lambda\) includes the Cauchy–Schwarz
bound on every set \(E\) with \(\int_E W^2dJ<\infty\).
Indeed \(\ell\) is \(r_*=\sqrt{2/3}\)-Lipschitz on the simplex,
so completing a square in \(|y-q|\) gives
\[
 H_\lambda(\min_i q_i,W)
       \ge W\ell(q)-\frac{r_*^2 W^2}{4\lambda}.
\]
Integrating this inequality on \(E\) and optimizing the two
terms involving \(\lambda\) gives the lower bound
\(\int_E\ell(q)d\mu_1-r_*K(\int_E W^2dJ)^{1/2}\).

For a signed trial with \(v\) sorted as in (SA1), the minimum
coordinate of \(q\) is
\[
     \min_i q_i=(1/3-A\,a_{\operatorname{sgn}\phi}|\phi|)_+,
     \qquad a_+=-v_3,\qquad a_-=v_1 .
                                                               \tag{CM3}
\]
When all trial coordinates are nonnegative the projection is the
identity. Otherwise at least one coordinate of the projection is
zero, since it subtracts a nonnegative threshold and clips at zero.
Thus on a direction bin, upper bounds on \(a_+\), \(a_-\), \(A\)
and \(|\phi|\) supply a directed lower bound on this minimum.
[Exact signed clipping](THREE_CELL_SIGNED_AXIS_EXTENSION.md)

The same pointwise argument permits a linear multiplier for forced
alignment. If \(\langle k,v\phi\rangle_J\ge L\), choose \(\eta\ge0\).
Replacing \(q\) in the quadratic distance by
\(q+\eta v\phi/(2\lambda)\) and minimizing the four affine branches
by Euclidean projection gives
\[
\begin{split}
 e_0\ge{}&
 \int \inf_{y\in\Delta_3}
   \left[W\ell(y)+\lambda
            \left|y-q-\frac{\eta}{2\lambda}v\phi\right|^2\right]dJ\\
 &-\lambda K^2+\eta L-\frac{\eta^2}{4\lambda}.
\end{split}                                                    \tag{CM4}
\]
Here \(\|v\phi\|_J=1\). The shifted reference need not belong to
\(\Delta_3\), so the simpler expression (CM1) is used only in the
unshifted bound.

These inequalities supply lower bounds on \(e_0\). An unrestricted
correlation interval still requires a complete spectral and
amplitude cover.

## Upper tilted moment from joint clipping and trial orthogonality

Let \(C=\Delta_3-p_*\), \(y=Av\phi\) and \(z=p-p_*=y+g\) with
\(\langle g,y\rangle_J=0\) and \(\|g\|_J^2\le A^2r^2\). For
\(\lambda>0\) and \(\zeta\le\lambda\), put
\[
 B_{\lambda,\zeta}(W,y)=
 \sup_{x\in C}
  \bigl[(W-\lambda)|x|^2+2(\lambda-\zeta)\langle x,y\rangle
                          +(2\zeta-\lambda)|y|^2\bigr].
\]
Subtract the quadratic distance and the axial pairing from the
moment, then restore their integrated bounds. This gives
\[
 \int |z|^2d\mu_1
   \le \lambda A^2r^2+\int B_{\lambda,\zeta}(W,y)dJ .           \tag{CM5}
\]
For \(W<\lambda\), completing the square evaluates the supremum as
\[
 \left(W+\frac{(W-\zeta)^2}{\lambda-W}\right)|y|^2
   -(\lambda-W)\operatorname{dist}^2
           \left(\frac{\lambda-\zeta}{\lambda-W}y,C\right).
                                                               \tag{CM6}
\]
The signed ray of the shifted point agrees with that of \(y\)
when the point is nonzero, so the distance is given by (SA1).
For \(W\ge\lambda\) the quadratic is convex, so the maximum is
attained at a vertex \(e_i-p_*\). Its value is
\[
   \frac23(W-\lambda)+2(\lambda-\zeta)\max_i y_i
                             +(2\zeta-\lambda)|y|^2 .
                                                               \tag{CM7}
\]
If an interval enclosure for \(W\) contains \(\lambda\), a
division-free upper bound follows by replacing \(W-\lambda\)
in (CM7) by \((W-\lambda)_+\). This uses
\(|x|^2\le2/3\) and \(\langle x,y\rangle\le\max_i y_i\).
The two multipliers in (CM5) allow the tilt and the signed
clipping geometry to enter the same bound.
