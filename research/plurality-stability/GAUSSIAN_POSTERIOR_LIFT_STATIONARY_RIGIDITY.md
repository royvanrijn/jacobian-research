# Only simplex posterior lifts can be stationary

Fix \(k\ge3\), a finite base dimension \(d\), and parameters
\(0<s,\rho<1\). For an ordinary balanced Gaussian partition
\(A\) in \(\mathbb R^d\), put \(p(x)=T_s1_A(x)\).
Let \(H=\mathbf1^\perp\subset\mathbb R^k\), and let
\(b(p)\in H\) be the centered price vector for which the
winner of \(z_i-b_i(p)\), for standard Gaussian \(z\) on
\(H\), has masses \(p\). The posterior simplicial lift
\(F=\mathcal L_s A\) assigns \((x,z)\) this winning label at
\(p=p(x)\).

If \(F\) is stationary for Gaussian noise stability at
correlation \(\rho\), then \(A\) is the centered regular
simplex, up to orthogonal transformations, relabeling and
null sets. The stationarity hypothesis here is equality of
the Gaussian score differences to fixed label-price
differences on all regular shared facets. No stability
or maximality assumption on \(F\) is needed. For \(d<k-1\)
this stationarity is impossible. Conversely a centered simplex
lifts to a centered simplex, which is stationary.
[The posterior simplicial lift and its centroid identities](GAUSSIAN_WEIGHTED_CENTROID_INEQUALITY.md#a-sharp-posterior-entropy-and-centroid-tradeoff)

## Global Lipschitz bounds for the model prices

Write \(B(x)=b(p(x))\), and let \(L_p\) be the facet-area
Laplacian of the translated regular simplex with masses
\(p\). Since every cell has positive mass, \(p(x)\) has
positive entries everywhere. The model price map is smooth
in the interior mass simplex and its mass Jacobian is
\(-L_p/\sqrt2\). This follows from the displacement of the
pair hyperplane \(z_i-z_j=b_i-b_j\) and strict positivity
of every pair facet. [Price derivative and Gaussian
multi-bubble profile](GAUSSIAN_SCORE_SUPPORT_FLOW.md#a-graph-formula-for-the-functional-perimeter-inequality).

For a unit observation direction \(v\), put
\(q_i=-\partial_v B_i/\sqrt2\). The chain rule gives
\(\partial_v p=L_p q\). Sort the entries of \(q\).
For an adjacent gap \(g\ge0\) in this list, let \(S\)
be the labels above it. If \(A_{ij}(p)\) is the model
facet area for \(i,j\), then
\[
 g\sum_{i\in S,\ j\notin S}A_{ij}(p)
 \ \le\ \sum_{i\in S}\partial_v p_i
 \ \le\ \frac{s}{\sqrt{1-s^2}}\,
             \phi\!\left(\Phi^{-1}\!\Big(\sum_{i\in S}p_i\Big)\right).
 \tag{PL1}
\]
The first inequality follows by summing the Laplacian
equations on \(S\). For the second, differentiating the
posterior probability of the union of the original labels
in \(S\) gives \(s/\sqrt{1-s^2}\) times a conditional
one-dimensional first moment. For a Gaussian set of
fixed measure, this directional moment is bounded above
by the moment of the upper halfspace of that measure:
subtract the threshold from the linear score and compare
the two indicators pointwise. The halfspace moment is
the displayed scalar profile.

The model union of labels in \(S\) has this same mass
and Gaussian perimeter \(\sum_{i\in S,j\notin S}A_{ij}(p)\).
Gaussian isoperimetry bounds the scalar profile by this
perimeter. Thus every adjacent gap is bounded by
\(s/\sqrt{1-s^2}\), and every price contrast obeys
\[
 \operatorname{Lip}(B_i-B_j)\le
                   \frac{\sqrt2(k-1)s}{\sqrt{1-s^2}} .
 \tag{PL2}
\]
The cut weights are positive for every nontrivial cut,
so division in (PL1) is valid. The estimate holds in
every observation direction and at every finite point.
[Gaussian isoperimetry as the two-cell
case](https://arxiv.org/html/1805.10961v3).

## Stationarity gives a scalar conditional median equation

Fix labels \(i\ne j\) and put \(B_{ij}(x)=B_i(x)-B_j(x)\).
The \(ij\) interface in the \(z\) fiber over any \(x\)
contains a nonempty relatively open subset of the hyperplane
\(z_i-z_j=B_{ij}(x)\). Its graph in \((x,z)\) is regular
there since the \(z\) normal is nonzero, and the two labels
strictly exceed every other label. If stationarity gives
facet score equality only almost everywhere on the regular
interface, continuity gives it on this full open portion.
For fixed \(x\), the score \(U_{ij}=T_\rho(1_{F_i}-1_{F_j})\)
is real analytic in \(z\), so the equality to the pair price
holds on the whole supporting hyperplane in that fiber.

Let \(P_H\) be projection onto \(H\) and
\(\tau=P_H(e_i+e_j)\). This direction is nonzero because
\(k\ge3\). It favors labels \(i,j\) equally over all other
labels. On the supporting hyperplane take observations
\(z_0+R\tau\) and \(z_0-R\tau\), with \(R\to\infty\)
and \(z_{0,i}-z_{0,j}=B_{ij}(x)\). Put
\(\sigma=\sqrt{1-\rho^2}\). The conditional base observation
is \(Y=\rho x+\sigma G\) with standard Gaussian \(G\).
For each fixed conditional draw of the base and label
noises, all model prices at \(Y\) are finite. In the
reverse direction neither \(i\) nor \(j\) eventually wins,
so dominated convergence makes the pair price zero. In the
forward direction the winning label is eventually one of
\(i,j\); the conditional difference in their raw label
coordinates is Gaussian with mean \(\rho B_{ij}(x)\) and
variance \(2\sigma^2\). Its winning threshold at the source
is \(B_{ij}(Y)\). Thus, writing \(b=B_{ij}\), one obtains
\[
 \mathbb E\,\psi\bigl(b(\rho x+\sigma G)-\rho b(x)\bigr)=0,
 \qquad \psi(t)=2\Phi\!\left(\frac{t}{\sqrt2\,\sigma}\right)-1 .
 \tag{PL3}
\]
The sign of the forward score is the negative of this
expectation. Only its vanishing is used in (PL3).

## A Lipschitz Liouville argument

Let \(b:\mathbb R^d\to\mathbb R\) be smooth and globally
Lipschitz and satisfy (PL3).

Differentiate (PL3) in a fixed direction \(v\), using the
same standard Gaussian \(G\) at all \(x\). The bounded
derivatives of \(b\) and bounded \(\psi'\) justify the
differentiation by domination. Since \(\rho>0\), the
bounded function \(f=\partial_v b\) satisfies
\(f(x)=\int f(y)K(x,dy)\) for the Markov kernel
\[
 K(x,dy)=
 \frac{\psi'(b(y)-\rho b(x))\,p_\sigma(y-\rho x)\,dy}
      {\int\psi'(b(z)-\rho b(x))\,p_\sigma(z-\rho x)\,dz},
                                                               \tag{PL4}
\]
where \(p_\sigma\) is the density of covariance \(\sigma^2 I_d\).
The kernel in (PL4) is reversible for the finite measure
\[
  e^{-b(x)^2/4}m(x)\,d\gamma_d(x),\qquad
  m(x)=\int\psi'(b(y)-\rho b(x))p_\sigma(y-\rho x)\,dy .
\]
Indeed the Gaussian transition is reversible for \(\gamma_d\),
and the additional factor in the reversible pair density is
\[
 e^{-b(x)^2/4}\psi'(b(y)-\rho b(x))
   =C_\rho\exp\left(-\frac{b(x)^2+b(y)^2-2\rho b(x)b(y)}
                            {4(1-\rho^2)}\right),
\]
which is symmetric. The finite invariant measure has a strictly
positive density and the transition kernel is everywhere positive.
Since \(f\) is bounded and \(Kf=f\), its reversible Dirichlet
energy is zero. Thus \(f(x)=f(y)\) for almost every pair
under the positive joint density, so \(f\) is constant.
Apply this to every direction \(v\) to conclude \(b\) is affine.

## Identification of the original partition

Apply the Liouville argument to every pair in (PL3).
All contrasts of the centered price map \(B\) are affine,
so \(B\) is affine. The lifted winner is an affine fan
with affinely independent sites, because the label
coordinate blocks of the sites are the projected vectors
\(e_i\) on \(H\). Its one-correlation stationary score
equations satisfy the hypotheses of affine fan rigidity.
Hence the lift is a centered regular simplex.
[Rigidity of stationary full-rank affine fans](AFFINE_SIMPLEX_FAN_RIGIDITY.md)

The affine site equalities now have constant pair
distances, and all their supporting hyperplanes contain
the origin. Thus \(B(0)=0\), and the row differences of the
linear price map all have the same length, since the
label coordinate blocks already have equal pair distances.
The map is centered, so the price matrix has row Gram
equal to a scalar multiple of \(P_H\). A zero scalar
would make \(T_s1_A\) constant, contradicting injectivity
of \(T_s\). Write \(B(x)=-\lambda O x\) where
\(\lambda>0\) and \(O:\mathbb R^d\to H\) is a coisometry.
The conditional label formula gives
\(p=T_r1_{T_O}\) for the centered simplex with active map
\(O\) and \(r=\lambda/\sqrt{1+\lambda^2}\).
Since \(p=T_s1_A\), semigroup injectivity implies \(r=s\):
if they differ, factor the operator with the smaller parameter
through the one with the larger and cancel the latter. An
ordinary nontrivial indicator would then equal a strictly
smoothed indicator. Cancelling at \(r=s\) gives \(A=T_O\).

The argument from (PL3) through affine fan rigidity also
applies to a winner of \(z_i-B_i(x)\) given any smooth
globally Lipschitz price map, without a posterior source.
The posterior structure supplies the Lipschitz estimate
(PL2) and the final identification of \(A\).
