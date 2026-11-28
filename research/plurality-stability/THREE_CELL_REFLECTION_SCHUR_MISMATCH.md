# An odd-kernel reflection test at equal posterior scores

Let \(A\) be an ordinary balanced Gaussian three-cell maximizer
at \(0<r<1\) in the plane. Fix a reflection \(R\) through the
origin and one of its open half-planes \(H\). Write
\(F=(1_{A_i})_{i=1}^3\),
\(\delta(x)=F(x)-F(Rx)\) on \(H\), and
\(D=\{x\in H:\delta(x)\ne0\}\). With the Gaussian-relative
Mehler kernel \(M_r\), define the odd kernel
\[
 W(x,y)=M_r(x,y)-M_r(x,Ry),\qquad x,y\in H .
\]
For \(p\in H\), let \(w(y)=W(p,y)\) and put
\[
 \begin{split}
 L&=\int_D w(y)\,d\gamma_2(y),\\
 J&=\int_D w(y)\delta(y)\,d\gamma_2(y),\\
 B&=\int_D\sqrt{W(p,p)W(y,y)-w(y)^2}\,d\gamma_2(y).
 \end{split}
\]
Then
\[
                    \|J\|^2+2B^2\ \ge\ \frac89 L^2.       \tag{SR1}
\]
For any balanced ordinary partition violating (SR1), a specified
reflection-pair relabeling strictly increases stability. The
comparison uses all three mismatch types.

If \(U_i(p)=U_i(Rp)\) for every \(i\), where
\(U_i=T_r1_{A_i}\), then \(J=0\) and (SR1) becomes
\[
                                 B\ge \frac23 L.          \tag{SR2}
\]
This applies to any two distinct equal-radius maximal triple
points, using the centered reflection that exchanges them.
The earlier targeted reflection test requires all three
mismatch types at that reflection. Inequality (SR2) also
constrains their spatial distribution.
[Targeted reflection at equal radii](THREE_CELL_POSTERIOR_DEGREE_AND_TRIPLE_BISECTOR.md#the-targeted-reflection-at-equal-vertex-radii)

## The bound from the odd kernel

For \(x,y\in H\), \(W(x,y)>0\), since the normal
coordinates of \(x,y\) are positive. The Mehler kernel is
positive semidefinite on finite point sets. Indeed it is the
kernel \(\exp(c\,x\cdot y)\), \(c=r/(1-r^2)>0\), multiplied
by a positive factor in each variable. The exponential has
a positive semidefinite power series of tensor inner products.
Pairing the antisymmetric point vectors for \(x,Rx\) and
\(y,Ry\) gives \(2W(x,y)\), so \(W\) is positive semidefinite.
In particular \(W(p,p)>0\) and the Schur complement
\[
                      W_\perp(x,y)=W(x,y)-\frac{w(x)w(y)}{W(p,p)}
\]
is positive semidefinite. This follows by adding the anchor
\(p\) to any finite point set and minimizing the resulting
nonnegative quadratic form over the coefficient at \(p\).
Hence
\[
                  |W_\perp(x,y)|
                     \le\sqrt{W_\perp(x,x)W_\perp(y,y)} .
\]
The terms are integrable because
\(\sqrt{W(y,y)}\le\sqrt{M_r(y,y)}\), which is a constant
times \(\exp[r|y|^2/(2(1+r))]\) and lies in \(L^1(\gamma_2)\).
The positive semidefinite kernels have nonnegative quadratic
forms on bounded measurable functions here. For example, compact
truncation, finite partitions, and continuity give this first on
bounded sets, and the diagonal bound justifies passage to the
limit.

Integrating over the two halves, the only part of \(S_r\) affected
by reflection-pair swaps is one half of
\[
 E(\delta)=\iint_{H^2} W(x,y)\delta(x)\cdot\delta(y)\,
                              d\gamma_2(x)d\gamma_2(y).
\]
Indeed the complementary term depends only on \(F(x)+F(Rx)\)
and the even kernel \(M_r(x,y)+M_r(x,Ry)\). Each nonzero
\(\delta\) has Euclidean norm \(\sqrt2\). The rank-one part
of \(W\) contributes \(\|J\|^2/W(p,p)\), and the diagonal
bound on the residual gives
\[
                         E(\delta)\le
                              \frac{\|J\|^2+2B^2}{W(p,p)} .
                                                               \tag{SR3}
\]

## A feasible finite relabeling

Fix cyclic label vectors
\(d_1=e_1-e_2,\ d_2=e_2-e_3,\ d_3=e_3-e_1\), and let
\(t_i\) be the integral of \(w\) over the mismatch type
\(\{\delta=\pm d_i\}\). For signs \(\epsilon_i\in\{-1,1\}\),
orient every reflection pair of type \(i\) to give
\(\widetilde\delta=\epsilon_i d_i\) on its part of \(H\).
This is an ordinary measurable relabeling preserving every
Gaussian mass, since each orbit retains its two labels and
\(R\) preserves Gaussian measure.

The function
\[
                f(t)=\max_{\epsilon}
                           \left\|\sum_{i=1}^3\epsilon_i t_i d_i\right\|
\]
is convex and symmetric in \(t\). Convexity follows from
the triangle inequality. For symmetry, the three \(d_i\)
are an equilateral triple in the centered label plane,
and permuting them is an orthogonal map there. The average
of all permutations of \(t\) is \((L/3,L/3,L/3)\), so
\[
                     f(t)\ge f(L/3,L/3,L/3)
                                 =\frac{2\sqrt2}{3}L .
\]
For the equality computation, the common-sign choices give
zero because \(\sum_i d_i=0\), and any other choice leaves
twice one of the three vectors up to sign. Choose maximizing
signs. Positive semidefiniteness of the Schur residual gives
\[
                  E(\widetilde\delta)\ge
                          \frac89\,\frac{L^2}{W(p,p)} .
\]
Together with (SR3) this proves (SR1) at an actual maximizer,
and gives the gain lower bound
\[
              S_r(\widetilde A)-S_r(A)\ge
                 \frac{(8/9)L^2-\|J\|^2-2B^2}{2W(p,p)}
\]
for any balanced ordinary partition. This relabeling is the
same feasible reflection-pair operation used in the exact
two-type comparison.
[Reflection-pair comparison, DG12](THREE_CELL_DIRECT_GLOBAL_INEQUALITY.md#exact-comparisons-within-reflection-orbits)

Finally, odd reflection gives \(J=(U_i(p)-U_i(Rp))_{i=1}^3\).
For a geometric form of (SR2) put
\[
                   h_p(y)=\frac{W(p,y)}{\sqrt{W(p,p)W(y,y)}} .
\]
Write \(u_x>0\) for the normal coordinate and \(v_x\) for the
tangential coordinate of \(x\in H\). The formula for Mehler gives
\[
 h_p(y)^2=e^{-c(v_p-v_y)^2}
       \frac{\sinh^2(cu_pu_y)}
            {\sinh(cu_p^2)\sinh(cu_y^2)} .                  \tag{SR4}
\]
The quotient extends continuously to \(u_p=0\) or \(u_y=0\)
by \(\sinh(t)/t\to1\). If \(D\) has positive measure and
\(h_p(y)>3/\sqrt{13}\) almost everywhere on \(D\), then
\(\sqrt{1-h_p(y)^2}/h_p(y)<2/3\) there. This contradicts
(SR2) when the anchored posterior vectors agree. The
conclusion excludes mismatch support confined to this
neighborhood in the normalized odd kernel. It leaves open
whether an unrestricted maximizing pair of junctions has
mismatch support this concentrated.


## A conditional slice witness without equal radii

For any distinct planar observations \(p,q\) with \(U_i(p)=U_i(q)\)
for all labels, put \(e=(p-q)/|p-q|\). For each pair \(i\ne j\),
a set of positive transverse measure in \(e^\perp\) consists
of offsets \(z\) for which the line \(z+\mathbb R e\)
meets both \(A_i\) and \(A_j\) in sets of positive line measure.
The assertion also applies to two distinct maximal triple points
of arbitrary radii.

To prove it, in \(\mathbb R^3\) consider the cylinder
\(\widehat A_i=A_i\times\mathbb R\). It is a global maximizer:
dimension reduction places the optimum for three cells in at
most two directions, and the cylinder has the planar value.
Its scores are \(U_i(x)\), independent of the added coordinate.
[Gaussian dimension reduction and existence](GAUSSIAN_NEAR_MAXIMIZER_COMPACTNESS.md)

Choose \(0<\alpha<1\), put \(b=\sqrt{1-\alpha^2}\), and let
\(n=(\alpha e,b)\) be the normal of a centered reflection
\(\widehat R\). For \(m=(p+q)/2\) and \(a=|p-q|/2\), set
\[
 \widehat m=(m,-\alpha (m\cdot e)/b),\qquad
 \widehat p=\widehat m+(a/\alpha)n,\qquad
 \widehat q=\widehat m-(a/\alpha)n .
\]
Then \(\widehat m\cdot n=0\), so \(\widehat R\widehat p=\widehat q\),
and the projections of these observations are \(p,q\). Their
cylindrical posterior vectors agree. If at most two mismatch
types occurred for \(\widehat R\), the two-type comparison
(DG12) would force a fixed sign for each mismatch coordinate
on the half-space containing \(\widehat p\). The strictly
positive odd kernel there would then force every mismatch
coordinate to vanish, since all the score differences between
\(\widehat p\) and \(\widehat q\) vanish. Reflection invariance
would put the cylinder's centroids in \(n^\perp\). They already
lie in the base plane, so their planar projections would all
be orthogonal to \(e\), contradicting their planar rank two.
[Full centroid rank](THREE_CELL_GLOBAL_OBSTRUCTIONS.md#one-active-dimension-is-impossible)

Thus each unordered pair \(i,j\) occurs on a set of reflection
orbits of positive Gaussian measure. In coordinates \(x=z+te\)
and added height \(h\), the reflection keeps \(z\) fixed and
changes the planar coordinate \(t\) to
\[
                         t'=(1-2\alpha^2)t-2\alpha b h.
\]
The map \((z,t,h)\mapsto(z,t,t')\) is invertible and carries
Gaussian measure to a measure equivalent to Lebesgue measure.
The half-space containing \(\widehat p\) becomes \(t>t'\).
For a given pair, Fubini yields a positive-measure set of \(z\)
for which the pairs \((t,t')\) with these two labels in either
order have positive two-dimensional measure. This requires
both line sections at that \(z\) to have positive measure,
as asserted.
