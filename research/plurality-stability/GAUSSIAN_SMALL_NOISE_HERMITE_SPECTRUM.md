# Hermite spectrum of moving partitions near correlation one

For a sequence that attains the sharp Gaussian perimeter limit at one
noise scale, the limit of its rescaled Hermite spectrum is determined by
the normals of the limiting interfaces. This applies to actual
maximizers and to sequences with disagreement within \(o(\sqrt{1-\rho})\)
of optimal in a fixed admissible dimension. The statement resolves the
joint distribution of degrees in the individual coordinates and keeps
the matrix of label correlations.

## Interface normals determine the rescaled degree law

Fix \(k,d\) and the coordinate axes in \(\mathbb R^d\). Suppose \(t_j>0\)
tends to zero and ordinary partitions \(A^{(j)}\) converge in Gaussian
label distance to a finite-perimeter partition \(A\). Write
\(N_\rho=1-S_\rho\), \(P_\gamma(A)=\sum_{i<l}\nu_{il}(\Sigma_{il})\), and take
\[
 \lim_j\frac{N_{e^{-t_j}}(A^{(j)})}{\sqrt{t_j}}
                       =\frac{2}{\sqrt\pi}P_\gamma(A).       \tag{HD1}
\]
Here \(\nu_{il}\) is the Gaussian surface measure of the reduced interface
\(\Sigma_{il}\). Let \(n_{il}\) be its unit normal in either orientation and
\(q(n)=(n_1^2,\ldots,n_d^2)\). For a multi-index \(\alpha\), use the
orthonormal product Hermite basis and write \(v_{j,\alpha}\in\mathbb R^k\)
for the vector of coefficients of the \(k\) cell indicators. Put
\[
 \eta_j=\frac{1}{\sqrt{t_j}}\sum_{\alpha\ne0}
     (1-e^{-t_j|\alpha|_1})\,v_{j,\alpha}v_{j,\alpha}^{\mathsf T}
                               \delta_{t_j\alpha}.
\]
These finite positive semidefinite matrix measures converge weakly on
\([0,\infty)^d\) to the measure characterized by the following formula
for bounded continuous functions \(g\):
\[
 \int g\,d\eta
  =\frac{1}{2\pi}\sum_{i<l} (e_i-e_l)(e_i-e_l)^{\mathsf T}
        \int_{\Sigma_{il}}\int_0^\infty
             g(uq(n_{il}))(1-e^{-u})u^{-3/2}\,du\,d\nu_{il}.   \tag{HD2}
\]
In particular, writing \(W_m(A^{(j)})=\sum_{|\alpha|_1=m}|v_{j,\alpha}|^2\),
for every \(b>0\) one has
\[
 \lim_j\frac{1}{\sqrt{t_j}}\sum_{t_jm>b} W_m(A^{(j)})
                  =\frac{2P_\gamma(A)}{\pi\sqrt b}.          \tag{HD3}
\]
For a simplicial limit, (HD2) is supported on the rays given by the
finitely many facet normals in the chosen coordinates. For fixed
partitions, lower bounds of order \(m^{-1/2}\) for unweighted
Hermite tails of indicators appear in
[Kindler and O'Donnell, Theorem 2.9](https://www.cs.cmu.edu/~odonnell/papers/gaussian-noise-sensitivity.pdf).
The statement here concerns moving partitions at the scale of (HD1)
and determines the limit including its constants, label matrix and
coordinate directions.

### Midpoint measures with the label pair kept

Set \(\rho_j=e^{-t_j}\) and \(e_j=\sqrt{2(1-\rho_j)}\). Use the midpoint
\(w\) and the difference variable \(z\) with endpoints \(w\pm e_jz/2\).
Let \(\mu_{j,il}\) be the joint law of \((w,z)\) for covariance
\(\rho_j I\), restricted to the event that the endpoint labels are
\(i,l\) in either order, and divided by \(e_j\). Then
\(\mu_{j,i}:=\sum_{l\ne i}\mu_{j,\min(i,l),\max(i,l)}\) is the
difference-quotient measure of the \(i\)-th indicator, as in (P2).
The measure lower bound (P3) in fixed dimension gives
\[
 \mu_{j,i}\ \rightharpoonup\
 \phi(w)|z\cdot n_i(w)|\,d|D1_{A_i}|(w)\,d\gamma_d(z)=:\mu_i .
\]
Indeed any vague cluster limit dominates \(\mu_i\). The sum of their
total masses tends to \(2\sqrt{2/\pi}P_\gamma(A)\) by (HD1), since
\(e_j/\sqrt{t_j}\to\sqrt2\). This is the sum of the limiting masses.
The domination is therefore equality for every \(i\), with tightness.
The proof of (P3) uses local \(L^1\) convergence at a fixed
\(z\) to obtain distributional convergence of signed quotients, then
lower semicontinuity and Fatou to pass to nonnegative joint tests.
[Midpoint measure lower bound (P2)–(P3)](THREE_CELL_NEAR_ONE.md#compactness-and-the-sharp-perimeter-lower-bound)

The unordered pair measures also converge narrowly. Any joint cluster
limit of the finitely many measures is dominated in its \(il\)-entry
by both \(\mu_i\) and \(\mu_l\). The reduced interfaces give a disjoint
decomposition up to surface-null sets, so the \(il\)-entry is supported
on \(\Sigma_{il}\). The sum over the neighbors of \(i\) equals \(\mu_i\).
Restricting the sum to each interface then gives
\[
             \mu_{j,il}\ \rightharpoonup\
           |z\cdot n_{il}(w)|\,d\nu_{il}(w)\,d\gamma_d(z).
                                                               \tag{HD4}
\]

### Changing the covariance at the same scale

For a fixed symmetric positive definite \(H\), let \(R_{j,H}=e^{-t_jH}\)
be the cross covariance, with both marginals standard Gaussian.
With the same coordinates \((w,z)\), the two variables are independent,
and their covariances are
\[
             C_{j,H}=\frac{I+R_{j,H}}2,\qquad
             B_{j,H}=\frac{I-R_{j,H}}{1-\rho_j}.
\]
Let \(D_{j,H}\) be the ratio of the joint densities in these coordinates
for \(H\) and \(I\). On compact sets \(D_{j,H}\) converges uniformly to
the ratio \(D_H(z)\) of the \(N(0,H)\) and \(N(0,I)\) densities:
\(B_{j,H}\to H\) and \(C_{j,H}\to I\).

The reweighted integrals against (HD4) converge even when this ratio
grows as \(|z|\to\infty\). To prove it, write
\(M_j=\sqrt{K\log(1/e_j)}\) for any fixed \(K>4\), taking \(j\) large so
that \(e_j<1\). On \(|w|>M_j\), the total reweighted mass is at most
\(e_j^{-1}\Pr(|W_{j,H}|>M_j)\), where \(W_{j,H}\) has covariance \(C_{j,H}\).
Since \(C_{j,H}\preceq I\), this is at most \(2^{d/2}e_j^{K/4-1}\).
For \(|w|\le M_j\), choose \(c>\max(1,\lambda_{\max}(H))\).
The difference-variable contribution to \(D_{j,H}/D_{j,cI}\)
is at most a constant times \(e^{-\kappa|z|^2}\) for a
\(\kappa>0\) and all sufficiently large \(j\), because
\(B_{j,H}^{-1}-B_{j,cI}^{-1}\to H^{-1}-c^{-1}I\succ0\).
The midpoint covariance exponents change by \(O(t_j M_j^2)=o(1)\)
and the determinant ratios stay bounded. Thus
\[
         D_{j,H}(w,z)\le C e^{-\kappa|z|^2}D_{j,cI}(w,z)
                           \quad (|w|\le M_j).
\]
The sum of the \(D_{j,cI}\)-weighted masses equals
\(N_{e^{-ct_j}}(A^{(j)})/e_j\), bounded by
\(c N_{\rho_j}(A^{(j)})/e_j\) since
\(1-e^{-ct_j m}\le c(1-e^{-t_jm})\) for \(m\ge1\).
This controls the \(|z|\) tails uniformly. On bounded \(|z|\) and
\(|w|\le M_j\) the ratio \(D_{j,H}\) itself is bounded uniformly in \(j\).
Tightness in \(w\) now follows from (HD4). Compact approximation proves
the convergence, and
\[
 \lim_j\int D_{j,H}\,d\mu_{j,il}
  =\sqrt{\frac2\pi}\int_{\Sigma_{il}}
                    \sqrt{n_{il}^{\mathsf T}Hn_{il}}\,d\nu_{il}. \tag{HD5}
\]
The last step integrates the absolute normal displacement under
\(N(0,H)\).

For the diagonal matrix \(H=\operatorname{diag}(h)\), \(h_a>0\), the
multivariate Mehler multiplier on \(\alpha\) is \(e^{-t_j h\cdot\alpha}\).
The matrix with entries
\(\langle 1_{A^{(j)}_i},(I-T_{R_{j,H}})1_{A^{(j)}_l}\rangle\)
also equals half the expectation of the outer product of the difference
of the endpoints' label vectors. Using (HD5) and \(e_j/\sqrt{t_j}\to\sqrt2\),
\[
 \lim_j\int \frac{1-e^{-h\cdot x}}{1-e^{-|x|_1}}\,d\eta_j(x)
 =\frac{1}{\sqrt\pi}\sum_{i<l}(e_i-e_l)(e_i-e_l)^{\mathsf T}
          \int_{\Sigma_{il}}\sqrt{h\cdot q(n_{il})}\,d\nu_{il}. \tag{HD6}
\]
Each \(\eta_j\) puts zero mass at \(0\). For identification, subtract
the version of (HD6) for \(h\) from that for \(h+\mathbf1\).
The integrand difference is \(e^{-h\cdot x}\). The right sides are
the corresponding integrals for (HD2), since
\(\int_0^\infty (1-e^{-a u})u^{-3/2}\,du=2\sqrt{\pi a}\),
by integration by parts. The total masses also converge by (HD6)
at \(h=\mathbf1\).

For completeness, tightness follows by taking the trace Laplace
transforms at \(h=c\mathbf1\). For any \(R,c>0\) the mass on
\(\{|x|_1>R\}\) is at most
\((1-e^{-cR})^{-1}\int(1-e^{-c|x|_1})\,d\operatorname{tr}\eta_j\).
Its upper limit can be made arbitrarily small by first choosing \(c\)
small and then \(R\) large, using total-mass and transform convergence.
The candidate (HD2) is finite, so dominated convergence justifies the
limit as \(c\downarrow0\). The trace controls the total variation
of each matrix entry. Every subsequential weak limit has the same
Laplace transforms for strictly positive \(h\), and then also for
\(h\ge0\) by bounded convergence. These determine the finite measure:
under \(y_a=e^{-x_a}\) they give all monomial integrals on
\([0,1]^d\), and polynomials are uniformly dense there.
This proves (HD2). Divide its radial projection by \(1-e^{-u}\)
on \(u>b\), where the reciprocal is bounded and the boundary \(u=b\)
has zero limiting measure. Taking traces yields (HD3).

## Consequences for maximizing sequences

Let \(d\ge k-1\), let the mass vectors converge to an interior
vector \(p\), and suppose the disagreement at \(\rho_j=e^{-t_j}\)
is within \(o(\sqrt{t_j})\) of the optimum for those masses.
The midpoint compactness and recovery by the translated simplex,
combined with Gaussian multi-bubble perimeter minimality and equality,
give convergence up to orthogonal alignment to the simplicial
partition at \(p\) and give (HD1). These are the compactness and
energy steps in the small-noise price proof. For actual global
maximizers in varying admissible dimensions, first reduce to
dimension \(k-1\) by the positive-correlation dimension theorem.
The statement about multi-indices uses the aligned active coordinates.
The total-degree statement is invariant under rotations and adding
unused coordinates.
[Moving masses and the Gaussian multi-bubble limit](GAUSSIAN_ALL_CELL_SMALL_NOISE_PRESSURES.md#perimeter-convergence-at-moving-masses)
[Heilman–Tarter, Theorem 1.9](https://arxiv.org/pdf/2011.05583v2)

The noise profile and every fixed positive-order derivative at
proportional scales converge. In particular
\[
 \frac{1-S_{e^{-ct_j}}(A^{(j)})}{\sqrt{t_j}}
          \longrightarrow \frac{2P_\gamma(A)}{\sqrt\pi}\sqrt c,
 \qquad
 \sqrt{t_j}\,S'_{e^{-ct_j}}(A^{(j)})
          \longrightarrow \frac{P_\gamma(A)}{\sqrt{\pi c}}.     \tag{HD7}
\]
For the first limit take \(H=cI\) and sum (HD5). For the derivatives,
differentiate the scalar expression in (HD6) with \(h=c\mathbf1\).
For order \(m\ge1\) the differentiated function of \(u=|x|_1\) is
\((-1)^{m-1}u^m e^{-cu}/(1-e^{-u})\). It extends continuously
to zero and is bounded. The bound and continuity in \(c\) hold in
the uniform norm on compact positive \(c\)-intervals, by the same
bound for the next derivative. Weak convergence in (HD2) gives
locally uniform convergence of each positive-order derivative
to the derivative of the right side, and the chain rule gives
the displayed \(S'\) limit. Mixed derivatives in positive diagonal
\(h\) follow with the product of coordinate powers in place of \(u^m\).
The finite-measure convergence gives leading orders at
\(t_j|\alpha|_1\) bounded away from zero. Exact finite-correlation
optimality requires a deficit estimate beyond these common limits.
