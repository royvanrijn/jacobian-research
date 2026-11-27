# Polynomial covariance control for nearly equal Gaussian cells

For \(k\) labels, put \(n=k-1\), \(P=I-\mathbf1\mathbf1^T/k\),
\(Q_*=P/n\), \(u=\mathbf1/k\), and \(F_k=(\mathbb E M_k)^2/n\),
where \(M_k\) is the maximum of \(k\) independent standard Gaussians.
For \(p\) a probability vector and \(Y\sim N(0,Q)\) define the
prescribed-mass score value
\[
 C_p(Q)=\inf_{\lambda\in\mathbb R^k}
   \left(p\cdot\lambda+\mathbb E\max_i(Y_i-\lambda_i)\right).
\]
We prove a polynomial quantitative covariance bound for this value
using the Gaussian multi-bubble perimeter theorem and Zhang's
[score covariance comparison](sources/zhang-equal-gaussian-cells-ebfd905.md).
The new estimate comes from a lower bound on every price Hessian
weight along the initial covariance path.

**Covariance bound.** Fix \(B_0\ge1\). There exist \(c,C,\epsilon>0\)
depending only on \(B_0\) and an integer \(k_0(B_0)\) such that for
\(k\ge k_0(B_0)\), centered positive semidefinite \(Q\) of trace one
with \(Q_{ii}\le B_0/k\), and \(\max_i|kp_i-1|\le\epsilon\),
\[
 C_p(Q)^2\le F_k + C(\log k)^2\|p-u\|_2^2
   -c\left(\|\operatorname{diag}(Q-Q_*)\|_2^2+
                         \frac{\|Q-Q_*\|_F^2}{k}\right).
 \tag{VC1}
\]
The mass perturbation is quadratic. The relevant coercivity has
different scales on diagonal and off-diagonal entries.

## Price weights at the start of the covariance path

Write \(z=\Phi^{-1}(1-1/k)\), \(H=Q-Q_*\), and \(Q_t=Q_*+tH\).
Add an independent common Gaussian coordinate to the scores and
scale by \(\sqrt n\) to use a Gaussian vector \(X\) with covariance
\[
 R_t=nQ_t+\mathbf1\mathbf1^T/k
      =(1-t)I+t(nQ+\mathbf1\mathbf1^T/k).
 \tag{VC2}
\]
An additive common score leaves the maximizing label unchanged.
The covariance in parentheses is positive semidefinite with
diagonal bounded by \(B_0+1\). Thus all diagonal entries of
\(R_t\) are \(1+O_{B_0}(t)\) and off-diagonal entries \(O_{B_0}(t)\).
The same bounds hold for conditional variances and covariances on
remaining coordinates after fixing one or two coordinates. This
follows from the Schur complement because the conditioning covariance
matrix has dimension at most two and is \(I+O_{B_0}(t)\).
Conditional means are bounded by \(C_{B_0}t\) times the sum of the
absolute values of the fixed coordinates.

Let \(P_i(\lambda)\) be the probability that \(X_i-\lambda_i\) wins.
For the objective \(p\cdot\lambda+\mathbb E\max_i(X_i-\lambda_i)\),
the gradient is \(p_i-P_i(\lambda)\). Its Hessian is a Laplacian
with pair weights
\[
 \widehat w_{ij}(t,\lambda)
 = \int_{\mathbb R}
 f_{ij}(x+\lambda_i,x+\lambda_j)\,
 \Pr(X_l\le x+\lambda_l\ \forall l\notin\{i,j\}\mid
             X_i=x+\lambda_i,\ X_j=x+\lambda_j)\,dx ,
 \tag{VC3}
\]
where \(f_{ij}\) is the bivariate marginal density. For \(i\ne j\)
this follows by conditioning on \(X_i\) in \(P_i\) and differentiating
the \(j\)-th upper limit. The diagonal Hessian follows from the
invariance of winners under a common price shift.
The relevant integral derivatives are dominated in a neighborhood
of any fixed price vector by an integrable Gaussian along the
two-coordinate slice.

For all \(0\le t\le 1/z^2\) and all \(|\lambda_i|\le 1/z\), we have
\[
 \widehat w_{ij}(t,\lambda)\ge c_0 z/k^2
 \tag{VC4}
\]
for some \(c_0>0\) depending only on \(B_0\), once \(k\) is large
enough. To prove this, retain \(x\in[z+L/z,z+(L+1)/z]\) in (VC3).
Here \(L\) is a sufficiently large constant depending only on \(B_0\).
The conditional mean and variance bounds imply, uniformly in \(l\),
that the standardized threshold \(x+\lambda_l\) for the remaining
conditional variable is at least \(z+2/z\). Indeed the possible
loss relative to \(x\) is at most \(1/z+C_{B_0}t(x+1/z)\).
The union bound and \((1-\Phi)(z+2/z)\le e^{-2}/k\) give a
conditional probability at least \(1-e^{-2}\) in (VC3).
The bivariate density is at least \(c_{B_0}\phi(z)^2\) on the same
interval. Its covariance is \(I+O_{B_0}(t)\), and its inverse
has the same bound, so its exponential factor relative to
\(\phi(z)^2\) differs by at most a bounded factor since
\(t z^2\le1\) and \(|x+\lambda_i-z|=O_{B_0}(1/z)\).
The interval has length \(1/z\), and \(\phi(z)\ge z/k\) by the
Mills bound. This proves (VC4) independently of the smaller
parameter chosen next.

## The balancing prices stay in the box

For any prescribed \(\eta>0\), there are \(a=a(B_0,\eta)>0\) and
\(k_1(B_0,\eta)\) such that
\[
 \max_i|kP_i(0)-1|\le\eta
     \qquad (0\le t\le a/z^2,\ k\ge k_1).
 \tag{VC5}
\]
To prove this, condition on \(X_i=x\) and standardize the remaining
coordinates. Their correlations are \(O_{B_0}(t)\). The zero-price
thresholds are \(b_l=x+O_{B_0}(t|x|)\).

For the upper estimate on the contribution of \(x\ge0\), increase
all thresholds to \((1+C_{B_0}t)x\) and compare their conditional
correlation matrix to the equicorrelation matrix with off-diagonal
entry \(\theta=C_{B_0}t\). Take \(C_{B_0}\) to dominate all the
preceding constants. Entrywise increase of off-diagonal entries
increases a Gaussian lower-orthant probability. One way to check
the comparison is to interpolate the two positive definite matrices,
differentiate the Gaussian density in covariance, and integrate by
parts on the orthant. The derivative with respect to one symmetric
off-diagonal entry is the bivariate threshold density times the
remaining conditional probability. All terms in the interpolation
are therefore nonnegative for entrywise increase. Both matrices
here are positive definite for \(t\le a/z^2\) and large \(k\).

The equicorrelated maximum has law
\(\sqrt{1-\theta}M_{k-1}+\sqrt\theta\,G\) for a standard normal
\(G\) independent of \(M_{k-1}\). Integrating the upper comparison
over all \(X_i\), which enlarges its integral over \(X_i\ge0\),
gives at most
\[
 \mathbb E\big[1-\Phi(c_t M_{k-1})\big],\qquad
 c_t=
 \frac{\sqrt{1-\theta}}
      {((1+C_{B_0}t)^2\operatorname{Var}(X_i)+\theta)^{1/2}}
 \ge 1-C'_{B_0}t .
\]
On \(0\le M_{k-1}\le2z\), the normal hazard rate bound
\(\phi(y)/(1-\Phi(y))\le C(1+y)\) for \(y\ge0\) shows that
\(1-\Phi(c_t M_{k-1})\le e^{C_{B_0}a}(1-\Phi(M_{k-1}))\).
The contribution on \(M_{k-1}<0\) is at most \(2^{-(k-1)}\).
For \(M_{k-1}>2z\) it is bounded by \(1-\Phi(2z(1-C_{B_0}'t))\),
which is \(o(1/k)\) since \(z^2\sim2\log k\).
The exact identity \(\mathbb E[1-\Phi(M_{k-1})]=1/k\) now
bounds this part of \(P_i(0)\) by \(e^{C_{B_0}a}/k+o(1/k)\).

For \(x<0\) the conditional standardized thresholds are nonpositive
for large \(k\). Replacing them by zero and using the same comparison
gives a contribution at most
\[
 \Pr(M_{k-1}\le z/2)+
 \Pr(\sqrt\theta\,G\le -\sqrt{1-\theta}z/2)=o(1/k),
\]
where the second term is zero when \(\theta=0\). For the first term,
write \(Q_{\rm N}=1-\Phi\) for the normal tail. Translation in the
normal tail integral gives \(kQ_{\rm N}(z/2)\ge e^{3z^2/8}\), so
\(\Pr(M_{k-1}\le z/2)\le e^{-(k-1)Q_{\rm N}(z/2)}\).
The second term is at most \(e^{-c_{B_0}z^4/a}\) by the Gaussian
tail bound. Taking \(a\) small proves the upper half of (VC5).

For the lower estimate choose a fixed compact window
\([u_0,u_1]\subset(0,\infty)\) with
\(\int_{u_0}^{u_1}e^{-v}dv\) as close to 1 as desired.
Restrict the \(X_i\) integration to \(v=kQ_{\rm N}(x)\in[u_0,u_1]\).
The Mills bounds give \(x=z+O_{u_0,u_1}(1/z)\) uniformly there.
Thus
\[
 \frac{\phi_{\operatorname{Var}(X_i)}(x)}{\phi(x)}
      =1+O_{B_0,u_0,u_1}(a)
\]
for \(a\le1\), where the numerator denotes the centered Gaussian
density of the displayed variance. Similarly
\(Q_{\rm N}(b_l)=\big(1+O_{B_0,u_0,u_1}(a)\big)v/k\)
by the hazard rate bound. The independent lower-orthant probability
at the thresholds \(b_l\) is \(e^{-v}+O_{B_0,u_0,u_1}(a)+o(1)\)
uniformly over the window.

The actual orthant probability differs from the independent one by
at most \(C_{B_0,u_0,u_1}a\). Indeed, use the same Gaussian
covariance derivative along the segment to independence. Every
bivariate threshold density is at most
\(C_{B_0,u_0,u_1}\phi(z)^2=O_{B_0,u_0,u_1}(z^2/k^2)\)
using \(b_l=z+O_{B_0,u_0,u_1}(1/z)\) and conditional correlations
\(O_{B_0}(t)\). The sum of the absolute increments of correlations
is \(O_{B_0}(k^2t)\), and \(z^2t\le a\).
The substitution \(dv=-k\phi(x)dx\) in the restricted integral
now gives the lower half of (VC5), by first choosing the window,
then \(a\), then the lower bound on \(k\).

Choose \(\eta<c_0/4\) in (VC5) and restrict \(\epsilon<c_0/4\)
in the hypotheses of (VC1). Minimize the scaled score objective
on the compact box \(|\lambda_i|\le1/z,\ \sum_i\lambda_i=0\).
Its gradient at zero has each entry bounded in modulus by
\((\epsilon+\eta)/k\). If a minimizer had \(\lambda_i=1/z\),
then integrating the \(i\)-th gradient along the segment from zero
using (VC4) gives
\[
 (p_i-P_i(\lambda))-(p_i-P_i(0))
   =\int_0^1\sum_{j\ne i}\widehat w_{ij}(t,s\lambda)
                                   (\lambda_i-\lambda_j)\,ds
   \ge c_0/k .
\]
In particular the gradient there is positive. At any lower boundary
\(\lambda_i=-1/z\) the analogous calculation gives a negative
gradient. Since the gradient sums to zero, moving an upper-boundary
coordinate downward and another coordinate with negative gradient
upward strictly reduces the objective. The negative-gradient coordinate
cannot be at the upper boundary. The case of a lower boundary is
symmetric. The minimizer is therefore interior in the box and all
coordinates of its gradient vanish. Convexity makes it a global
minimizer. Hence the prices balancing to \(p\) satisfy (VC4)
through the entire early path \(0\le t\le a/z^2\).

## The covariance comparison at prescribed masses

Use the unscaled scores \(Q_t\) and denote their price Laplacian
by \(L_t\). The scaling in (VC2) gives
\(w_{ij}(t)=\sqrt n\,\widehat w_{ij}(t,\lambda_t)\).
Indeed the objective is multiplied by \(\sqrt n\) when prices
are multiplied by \(\sqrt n\), and the common Gaussian component
makes zero contribution to its expectation. Let
\(C(t)=C_p(Q_t)\), \(C_*=C_p(Q_*)\), and let \(S_*\) denote
the Gaussian perimeter of the shifted regular fan with masses \(p\).
Here are the derivative and flux identities for the prescribed
mass objective. For \(Q\) positive on \(E\), factor \(Q=BB^T\)
with \(B:\mathbb R^n\to E\) invertible. The rows \(b_i\) are
affinely independent scores. The minimizing prices exist:
in the gauge \(\min_i\lambda_i=0\), testing the largest
score against any zero-priced label gives an objective lower
bound \(p\cdot\lambda\ge(\min_i p_i)\max_i\lambda_i\).
Ties are Gaussian-null and differentiation in the prices
gives the winning masses \(p\). All pair facets have positive
area for every finite price vector, since the score
differences against one label give invertible coordinates.
Thus the price Hessian is strictly positive on the mean-zero
gauge. The prices are unique and locally smooth there by the
implicit function theorem. The needed smoothness of the
winner probabilities follows by expressing the winning regions
in the nonsingular score-difference coordinates. Their cones
are fixed in these coordinates, and local parameter derivatives
of the Gaussian densities admit integrable polynomial-Gaussian
majorants.

Write \(L_Q\) for the price Laplacian, with weights
\(w_{ij}=A_{ij}/\ell_{ij}\), where \(A_{ij}\) is the ordinary
Gaussian interface area of the winning cells and
\(\ell_{ij}=|b_i-b_j|\).
Gaussian integration by parts gives their moment matrix
\(M=L_Q B\). Indeed each inward unit normal of cell \(i\)
is \((b_i-b_j)/\ell_{ij}\) at its \(ij\) interface.
A ball cutoff removes the unbounded-boundary term by Gaussian
decay. The balancing prices cancel in the objective, giving
\(C_p(Q)=\langle M,B\rangle_F=\operatorname{tr}(L_Q Q)\).
For any centered covariance variation \(D=BRB^T\) with \(R\)
symmetric, realize \(Q+sD\) by \(B(I+sR)^{1/2}\).
The derivative at \(s=0\) of the score objective is the
pairing of the winning moment matrix with \(BR/2\).
The price derivative cancels because the masses equal \(p\).
Consequently
\[
 dC_p(Q)[D]=\frac12\operatorname{tr}(L_Q D).
 \tag{VC5a}
\]
This reconstructs the identities in
[Zhang's covariance method](sources/zhang-equal-gaussian-cells-ebfd905.md)
directly for the masses used here.

Since \(tH=Q_t-Q_*\), (VC5a) yields
\(2tC'(t)=C(t)-\operatorname{tr}(L_t)/n\). Thus for \(0<t<1\),
\[
 t(C^2-C_*^2)'=C^2-C_*^2-D(t),\qquad
 D(t)=C(t)\operatorname{tr}(L_t)/n-C_*^2 .
\]
At \(Q_*\) all edge lengths are \(\sqrt{2/n}\), so by flux
\(C_*^2=2 S_*^2/n\). For the winning cells at \(Q_t\), let
\(S(t)=\sum_e w_e(t)\ell_e(t)\) be their Gaussian perimeter.
The shifted regular fan minimizes perimeter among partitions
with these masses by
[Milman–Neeman, Theorem 1.1](https://annals.math.princeton.edu/2022/195-1/p02).
Using \(C(t)=\sum_e w_e(t)\ell_e(t)^2\),
\(\operatorname{tr}L_t=2\sum_e w_e(t)\), and \(S(t)\ge S_*\),
the weighted square expansion gives
\[
 D(t)=\frac2n\left[
 S(t)^2-S_*^2
 +\sum_{e<f}w_e(t)w_f(t)(\ell_e(t)-\ell_f(t))^2
 \right]\ge
 \frac2n\sum_{e<f}w_e(t)w_f(t)(\ell_e(t)-\ell_f(t))^2 .
\]
Unlike the equal-mass argument, a linear term remains at the initial
covariance. Integrating with this term retained gives
\[
 C_p(Q)^2
   = C_*^2+C_*\operatorname{tr}(L_0H)
                       -\int_0^1 \frac{D(t)}{t^2}\,dt .
 \tag{VC6}
\]
Smoothness at \(Q_*\) identifies the \(t\downarrow0\) limit of
\((C(t)^2-C_*^2)/t\) with \(C_*\operatorname{tr}(L_0 H)\).
Continuity at \(Q\) holds even if \(Q\) is singular by the
Gaussian square-root coupling, uniformly over prices.
These are the only endpoint limits in (VC6).

Put \(d_{ij}=H_{ii}+H_{jj}-2H_{ij}\), so that
\(\ell_{ij}(t)^2=2/n+td_{ij}\).
Since \(\max_i (Q_t)_{ii}=O_{B_0}(1/k)\), each edge length is at
most \(C_{B_0}/\sqrt k\). Consequently
\[
 (\ell_e(t)-\ell_f(t))^2\ge c_{B_0} k t^2(d_e-d_f)^2 .
\]
The zero trace and row sums of \(H\) give
\(\sum_e d_e=0\) and
\[
 \sum_{i<j} d_{ij}^2
     = k\|\operatorname{diag}H\|_2^2+2\|H\|_F^2 .
\]
For the second identity, expand the squares.
The diagonal-square term is \((k-2)\sum_i H_{ii}^2\) since
\(\sum_i H_{ii}=0\). The mixed term is \(4\sum_i H_{ii}^2\) by
the row sums. The remaining term is \(2(\|H\|_F^2-\sum_i H_{ii}^2)\).
With \(q=k(k-1)/2\), this yields
\(\sum_{e<f}(d_e-d_f)^2=q\sum_e d_e^2\).
On \(0<t\le a/z^2\), (VC4) now gives
\[
 \int_0^1 \frac{D(t)}{t^2}\,dt\ge
 c_{B_0}\left(\|\operatorname{diag}H\|_2^2+
                                      \|H\|_F^2/k\right).
 \tag{VC7}
\]

To bound the initial slope in (VC6), put
\(\delta=p-u\) and let \(\lambda_0\) be the centered, scaled
balancing price vector at \(Q_*\). At \(t=0\), the scaled scores
are independent standard Gaussians, so the objective gradient at
zero is \(\delta\). On the mean-zero gauge (VC4) gives a Hessian
lower bound \(c_0 z/k\) throughout the segment from 0 to
\(\lambda_0\). Since the gradient at \(\lambda_0\) vanishes,
\[
 \|\lambda_0\|_2\le C k\|\delta\|_2/z .
\]
For independent scores the scaled price weights are
\[
 \widehat w_{ij}(0,\lambda)=
  \int \phi(x+\lambda_i)\phi(x+\lambda_j)
                          \prod_{l\ne i,j}\Phi(x+\lambda_l)\,dx .
\]
Uniformly on \(|\lambda_i|\le1/z\), their derivatives in
\(\lambda_i,\lambda_j\) have modulus at most \(Cz^2/k^2\),
and their derivatives in other \(\lambda_l\) at most \(Cz^2/k^3\).
For these bounds put \(y=x+1/z\). Each density factor is bounded
by \(\phi(y)e^{2|y|/z}\), each CDF factor by \(\Phi(y)\).
On \(|y|\le2z\), the exponential factor is bounded, and a density
derivative adds a factor at most \(2z+2/z\).
The required integrals are then controlled by \(I_2(k)\) and
\(I_3(k)\) of [the facet asymptotic (VH2)](NEGATIVE_VOTING_HESSIAN.md#the-asymptotic-gap).
For \(y>2z\), the Gaussian powers \(\phi(y)^j e^{2j|y|/z}\)
and possible additional factor \(1+|y|\) have negligible
integrals relative to \(z^{j-1}/k^j\), even without the
CDF powers. For \(y<-2z\), the integrals are exponentially
small in \(k\) from \(\Phi(y)^{k-j}\le2^{-(k-j)}\).

Integrating those derivative bounds along the price segment and
summing squares over the pairs yields
\[
 \|\widehat w(0,\lambda_0)-\widehat w(0,0)\|_{\ell^2(e)}
       \le C\frac{z^2}{k^2}\sqrt k\,\|\lambda_0\|_2 .
\]
To see the scaling, the pairwise bound is
\(Cz^2/k^2\) times
\(|\lambda_i|+|\lambda_j|+\|\lambda\|_1/k\).
The edge norm of that factor is at most \(C\sqrt k\|\lambda\|_2\).

The zero-price Laplacian at regular scores is scalar on \(E\),
and hence its pairing with \(H\) vanishes. The bound
\(C_*\le\sqrt{F_k}=O(z/\sqrt k)\) follows by testing zero prices
in \(C_p(Q_*)\). Cauchy–Schwarz on the edges gives, with the
\(\sqrt n\) scaling of \(L_0\),
\[
 \begin{split}
 |C_*\operatorname{tr}(L_0H)|
 &\le C\frac{z}{\sqrt k}\sqrt n\,
      \|\widehat w(0,\lambda_0)-\widehat w(0,0)\|_{\ell^2(e)}
      \left(\sum_e d_e^2\right)^{1/2} \\
 &\le C z^2\|\delta\|_2
     \left(\|\operatorname{diag}H\|_2^2+
                                  \|H\|_F^2/k\right)^{1/2}.
 \end{split}
\]
Apply Young's inequality to absorb half of (VC7) in (VC6),
and use \(z^4=O((\log k)^2)\) to prove (VC1).

## Consequences for Gaussian assignments

Let \(f\) be a measurable ordinary or fractional Gaussian
assignment with masses \(p\). Let \(A\) have rows
\(a_i=\int xf_i(x)d\gamma_d(x)\), and put \(W=\|A\|_F^2\).
For \(W>0\), \(Q=AA^T/W\) is a centered trace-one covariance
matrix. Testing the weighted pointwise maximum against \(f\)
gives \(W\le C_p(AA^T)=\sqrt W C_p(Q)\). Consequently (VC1)
bounds \(W\) itself.

Suppose \(f\) is a possible improver over the simplex at \(-r\)
within \(0<r\le1/(k-1)\). The
[mass estimate (M1)](NEGATIVE_VOTING_MINIMIZERS.md#rebalancing-possible-improvers)
gives \(\|p-u\|_2=O(r\log k/\sqrt k)\), and hence
\(\max_i|kp_i-1|\le C r\sqrt{k}\log k=o(1)\) uniformly.
The necessary bounds in the
[negative voting comparison](NEGATIVE_VOTING_INTERVAL.md#v4-first-moments-of-possible-improving-assignments)
give \(W\ge F_k+\|p-u\|_2^2/r-O(r/k)\).
In particular \(W\gtrsim\log k/k\).
The halfspace bound \(|a_i|\le\phi(\Phi^{-1}(p_i))\) and
\(p_i=O(1/k)\) give \(Q_{ii}\le B_0/k\) for an absolute
\(B_0\). Fix \(B_0\) at this value before choosing \(\epsilon\)
in (VC1). Combining (VC1) with the necessary bound absorbs
the \(C(\log k)^2\|p-u\|_2^2\) term into
\(\|p-u\|_2^2/r\) uniformly throughout the voting interval.
The result is
\[
 \frac{\|p-u\|_2^2}{r}
   +\|\operatorname{diag}(Q-Q_*)\|_2^2
   +\frac{\|Q-Q_*\|_F^2}{k}\ \le\ C\frac r k .
 \tag{VC8}
\]
In particular the centroid effective rank
\(\operatorname{tr}(Q)^2/\|Q\|_F^2\) is at least \(c k\),
since \(\|Q\|_F^2=1/(k-1)+\|Q-Q_*\|_F^2=O(1/k)\).
This applies to balanced and free-mass possible improvers.

Combining (VC1) with the signed second-degree estimate gives a
[sharper label localization](NEGATIVE_VOTING_INTERVAL.md#dimension-uniform-localization-on-the-voting-scale).
This estimate is used in the [exact large-alphabet same-partition theorem](NEGATIVE_VOTING_INTERVAL.md#exact-gaussian-simplex-minimization).

## Numerical constants for the covariance estimate

For \(B_0=2000\), the constants in (VC1) may be taken as
\[
 k_0=\lceil e^{1000}\rceil,\qquad
 \epsilon=e^{-50},\qquad c=e^{-1201},\qquad C=e^{1500}.
 \tag{VC9}
\]
The inequalities below only use \(\log k\ge1000\). Put \(a=e^{-1000}\)
and restrict the early path to \(t\le a/z^2\).
Mills' inequalities give
\[
 \log k\le z^2\le2\log k,\qquad z/k\le\phi(z)\le 2z/k .
\]
For the lower bound on \(z^2\), at \(y=\sqrt{\log k}\) use
\(Q_{\rm N}(y)\ge \phi(y)y/(1+y^2)\ge 1/k\).
The remaining bounds follow from the Gaussian Chernoff and Mills
bounds.

The entries of \(R_t-I\) in (VC2) have absolute value at most
\(2001t\). Conditioning on one or two coordinates uses a covariance
block with inverse norm at most two, since the block differs from
identity in operator norm by at most \(4002t\). Its conditional
mean in any remaining coordinate is at most \(10^5t\) times the
largest absolute value of the fixed coordinates, and its
conditional variances differ from one by at most \(10^5t\).
After conditioning on one coordinate the remaining standardized
correlations have absolute value at most \(10^5t\).
These bounds follow directly by bounding the conditioning row
in Euclidean norm by \(\sqrt2\cdot2001t\). For a conditional
threshold of the form \(x\) after fixing \(X_i=x\), its
standardized value differs from \(x\) in absolute value by at
most \(10^5t|x|\), using the one-coordinate Schur complement.

For \(|\lambda_i|\le1/z\) in (VC3) restrict the tie height to
\([z+10/z,z+11/z]\). The preceding conditional estimates
give other standardized thresholds at least \(z+7/z\).
The two fixed coordinates have squared radial length at most
\(2z^2+49\), and the corresponding inverse covariance differs
from identity by at most \(10^4t\). The density is therefore
at least \(e^{-28}\phi(z)^2\). The conditional union bound is
at most \(e^{-7}\), using
\(Q_{\rm N}(z+7/z)\le e^{-7}/k\).
Integrating on the height interval proves (VC4) with \(c_0=e^{-40}\)
throughout this early path.

For zero prices, one can take \(\eta=e^{-60}\) in (VC5).
For its upper estimate use equicorrelation
\(\theta=10^5t\) and conditional threshold \((1+10^5t)x\)
on \(x\ge0\). The resulting \(c_t\) is at least \(1-10^7t\).
On \(0\le M_{k-1}\le2z\), the normal hazard rate bound by
\(2(1+y)\) shows
\(Q_{\rm N}(c_t M_{k-1})\le e^{10^9 a} Q_{\rm N}(M_{k-1})\).
The contribution for \(M_{k-1}<0\) is at most \(2^{1-k}\)
and for \(M_{k-1}>2z\) is at most \(e^{-1.9 z^2}\).
The contribution from \(x<0\) is bounded by
\[
 \exp\bigl(-\tfrac12e^{3z^2/8}\bigr)+
 \exp\bigl(-z^2/(16\theta)\bigr),
\]
with the second term zero at \(\theta=0\).
Indeed \(kQ_{\rm N}(z/2)\ge e^{3z^2/8}\) and the independent
common normal would have to fall below
\(-\sqrt{1-\theta}z/(2\sqrt\theta)\) otherwise.
Each of these four exterior contributions is at most
\(e^{-65}/k\) for \(\log k\ge1000\), while
\(e^{10^9a}-1\le e^{-65}\).
Using \(\mathbb E Q_{\rm N}(M_{k-1})=1/k\) gives
\(kP_i(0)\le1+e^{-60}\).

For the lower estimate restrict to \(v=kQ_{\rm N}(x)\) in
\([e^{-80},100]\). Gaussian tail translation locates \(x\) between
\(z-20/z\) and \(z+100/z\).
The conditional standardized thresholds all lie in
\([z-21/z,z+101/z]\). In this box the bivariate densities
that enter the orthant comparison are at most
\(e^{50}\phi(z)^2\), since their covariance differs
from identity by at most \(10^5t\) in operator norm.
The sum of their products with absolute correlation
increments is at most \(2\cdot10^5 e^{50}a\).
For the independent part, the change from the threshold \(x\)
to each of the standardized conditional thresholds changes
its normal tail by a relative factor between
\(e^{-10^7a}\) and \(e^{10^7a}\), again by the hazard bound.
The product changes in absolute value by at most \(10^{10}a\)
because \(v\le100\). Also
\(\left|(1-v/k)^{k-1}-e^{-v}\right|\le 3\cdot10^4/k\)
by the Taylor bound for \(\log(1-v/k)\).
The density ratio \(\phi_{\operatorname{Var}(X_i)}(x)/\phi(x)\)
is between \(e^{-10^7a}\) and \(e^{10^7a}\) throughout the window.
Inserting these bounds in the \(dv=-k\phi(x)dx\) integral
and using
\(\int_{e^{-80}}^{100}e^{-v}dv\ge1-e^{-80}-e^{-100}\)
proves \(kP_i(0)\ge1-e^{-60}\).

The proof of the balancing box now applies because
\(\epsilon+\eta<e^{-40}\). The unscaled facet weights on the
early path are at least \(\sqrt n\,e^{-40}z/k^2\).
All squared length differences satisfy
\[
 (\ell_e(t)-\ell_f(t))^2\ge
 \frac{k}{32000}t^2(d_e-d_f)^2
\]
by the bound on the diagonal of \(Q_t\). Using
\(q=k(k-1)/2\ge k^2/4\) in the edge variance computation and
integrating \(t\) over \([0,a/z^2]\) gives
\[
 \int_0^1D(t)t^{-2}\,dt
 \ge e^{-1200}
       \left(\|\operatorname{diag}H\|_2^2+\|H\|_F^2/k\right).
\]

For the slope estimate at \(Q_*\), use the weight derivatives
at scaled prices in the balancing box. Their absolute bounds
are \(e^{100}z^2/k^2\) for the two tied-coordinate derivatives
and \(e^{100}z^2/k^3\) for each other derivative.
For a check of the constants, the Mills substitution in (VH2)
gives \(I_j(k)\le e^{20}z^{j-1}/k^j\) for \(j=2,3\).
On \(y\ge0\) one can use
\(\phi(y)\le2(1+y)Q_{\rm N}(y)\) and
\((1+y)/z\le 4(1+\sqrt{\log_+(1/v)})\) after
substituting \(v=kQ_{\rm N}(y)\).
Then \((1-v/k)^{k-j}\le e^{-v/2}\) and
\[
 \int_0^\infty v^{j-1}
       (1+\sqrt{\log_+(1/v)})^{j-1}e^{-v/2}\,dv\le 20 .
\]
The negative half contributes at most \(2^{-(k-j)}\)
to the original \(I_j\). For the shifted price integrals use
\(y=x+1/z\). On \(|y|\le2z\) each density factor is at most
\(e^4\phi(y)\) and a tied-coordinate derivative adds at most
\(3z\). These contributions fit within \(e^{100}\) times
the required powers of \(z\) and \(k\). For \(y>2z\),
\(6y/z\le y^2/16\), so the product of density factors is
bounded by \(e^{-(j/2-1/16)y^2}\).
Integrating with a possible factor \(1+y\) gives at most
\(100e^{-(2j-1/2)z^2}\) by retaining \(e^{-y^2/16}\)
for integration. This is negligible compared with
\(z^{j-1}/k^j\). For \(y<-2z\) use
\(\Phi(y)^{k-j}\le2^{-(k-j)}\) and integrability of the
shifted Gaussian density products.

The balanced price norm therefore obeys
\(\|\lambda_0\|_2\le e^{40}k\|\delta\|_2/z\) and
\(\|\widehat w(0,\lambda_0)-\widehat w(0,0)\|_{\ell^2(e)}
 \le e^{142}z\|\delta\|_2/\sqrt k\).
Also \(C_*\le 2z/\sqrt k\) using
\(\mathbb E M_k\le z+1/z\), from integrating the union bound on
the maximum above \(z\).
The initial slope in (VC6) thus has absolute value at most
\[
 e^{144}z^2\|\delta\|_2
          \left(\|\operatorname{diag}H\|_2^2+\|H\|_F^2/k\right)^{1/2}.
\]
Absorb this by Young's inequality with half the early-path
integral. Since \(z^4\le4(\log k)^2\), the result gives
the constants in (VC9).
