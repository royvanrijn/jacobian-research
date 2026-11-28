# Gaussian perimeter decreases under the posterior simplex lift

Let \(A\) be an ordinary Gaussian partition with \(k\ge2\)
positive masses in any finite dimension. Fix \(0<s<1\) and
write \(p(x)=T_s1_A(x)\) for its posterior mass vector.
The posterior simplex lift \(\mathcal L_sA\) uses a fresh
standard Gaussian on the label hyperplane and assigns the
label of the model simplicial cluster of masses \(p(x)\)
in that fiber. It has the same unconditional label masses
as \(A\). Write \(P_\gamma\) for Gaussian cluster perimeter,
counting each interface once, and let \(I(p)\) be the model
perimeter. Then
\[
 P_\gamma(\mathcal L_s A)
 \ \le\ \frac{\mathbb E I(p(X))}{\sqrt{1-s^2}}
 \ \le\ P_\gamma(A).                                        \tag{PL1}
\]
The lift always has finite perimeter. If \(A\) has finite
perimeter, equality of its perimeter with the lift occurs
exactly for regular simplicial clusters, allowing cylindrical
extensions, relabeling and null sets.

In particular for every fixed nonsimplicial ordinary \(A\)
and fixed \(s\in(0,1)\), there is a \(\rho_0(A,s)<1\) such that
\[
                   S_\rho(\mathcal L_s A)>S_\rho(A)
             \qquad(\rho_0(A,s)<\rho<1).                    \tag{PL2}
\]
This is a fixed-input interval, with no uniform assertion as
the partition or the lift parameter varies.
[Posterior simplex lift and its first moments](GAUSSIAN_WEIGHTED_CENTROID_INEQUALITY.md#a-sharp-posterior-entropy-and-centroid-tradeoff)

## The graph and the conditional perimeter

Set \(\sigma=\sqrt{1-s^2}\) and write \(q(p)=\nabla I(p)\)
on the label hyperplane. The lifted pair interfaces are
the graphs of \(q_i(p(x))-q_j(p(x))\) in the normalized
pair-score coordinate of the fresh Gaussian. If \(A_{ij}(p)\)
is the model facet area, their Gaussian area is
\[
 P_\gamma(\mathcal L_s A)
    =\mathbb E\sum_{i<j} A_{ij}(p(X))
        \sqrt{1+|\nabla_x(q_i(p(X))-q_j(p(X)))|^2}.
\]
The graph formula is valid for all positive-mass ordinary
partitions, since \(p\) is smooth and interior at every
observation point. The area formula on bounded sets and
exhaustion allow infinite area before a bound is established.
[Gaussian graph perimeter formula, (PF1)–(PF2)](GAUSSIAN_SCORE_SUPPORT_FLOW.md#a-graph-formula-for-the-functional-perimeter-inequality)

At each \(x\), the sum of the squared graph slopes weighted
by \(A_{ij}(p)\) is
\[
    \operatorname{tr}\big[(Dp)^{\mathsf T}(-D^2I(p))Dp\big]
        \le \frac{s^2}{1-s^2}I(p).
\]
The identity uses \(D^2I(p)=-L_p^{-1}\), and the upper
bound is the sharp weighted centroid inequality applied
to the conditional partition. Cauchy–Schwarz for the
weighted square roots now bounds the graph area integrand
by \(I(p)/\sigma\). Since \(I\) is bounded, this proves the
first inequality and finiteness of the lifted perimeter.
For a nonsimplicial \(A\) the pointwise weighted centroid
comparison is strict, because conditional partitions are
similarities of \(A\). The first inequality is then strict.
[Weighted posterior gradient bound, (WC4)](GAUSSIAN_WEIGHTED_CENTROID_INEQUALITY.md#posterior-gradients-in-the-perimeter-metric)

If \(A\) has finite Gaussian perimeter, the standardized
conditional partition is
\((A-sx)/\sigma\), of masses \(p(x)\). The minimum perimeter
\(I(p(x))\) is at most its perimeter. If the original
dimension is less than \(k-1\), embed the conditional
partition cylindrically before invoking the perimeter
theorem. Under this similarity the conditional perimeter
is
\[
 P_\gamma((A-sx)/\sigma)
   =\sigma\int p_{\sigma^2}(y-sx)\,d\mathcal H^{d-1}(y),
\]
where \(p_{\sigma^2}\) is the centered Gaussian Lebesgue
density of variance \(\sigma^2 I\) and the integral is
over the original reduced interfaces, counted once.
This is the surface Jacobian \(\sigma^{-(d-1)}\) together
with Gaussian density rescaling. The transformation applies
to finite-perimeter representatives. Finite Gaussian
perimeter gives locally finite Euclidean perimeter, and
the nonnegative conditional integral can be averaged
without a prior finiteness bound. Since \(sX+\sigma Z\) is
standard Gaussian for independent standard \(X,Z\),
Tonelli's theorem gives
\[
           \mathbb E P_\gamma((A-sX)/\sigma)=\sigma P_\gamma(A).
\]
This proves the second inequality of (PL1). If the original
perimeter is infinite that inequality is automatic.
[Gaussian multi-bubble perimeter and its equality cases](https://annals.math.princeton.edu/2022/195-1/p02)

For a regular simplicial cluster in active coordinates,
write its label as the winner of \((Ox)_i-b_i\), where
\(O\) is a coisometry onto the label hyperplane. Its
conditional winner is the model translated simplex of
masses \(p(x)\), with price \((b-sOx)/\sigma\).
Thus the lift assigns the winner of \(sOx+\sigma z-b\).
This is the original simplicial rule composed with a
coisometry on the product space. Perimeter and noise
stability are both preserved in this case.

## Heat content and strict improvement near one

For an ordinary Gaussian partition \(C\) with finitely many cells, write
\(t=-\log\rho\). The Ornstein–Uhlenbeck heat-content formula
gives the extended-real limit
\[
      \lim_{\rho\uparrow1}\frac{1-S_\rho(C)}{\sqrt{-\log\rho}}
                   =\frac{2}{\sqrt\pi}\,P_\gamma(C).
                                                               \tag{PL3}
\]
Indeed \(1-S_\rho(C)=\sum_i\int_{C_i^c}T_\rho1_{C_i}\,d\gamma\).
Equation (3.1) of Ambrosio–Miranda–Pallara gives the
corresponding limit for each finite-perimeter cell with
constant \(1/\sqrt\pi\), and their finite-liminf criterion
gives divergence for an infinite-perimeter cell.
Summing over finitely many cells counts each interface
twice in the finite case.
[Ambrosio–Miranda–Pallara, §3, equation (3.1) and its converse](https://cvgmt.sns.it/media/doc/paper/2601/finepropertiesfin.pdf)

If \(A\) has finite perimeter and is nonsimplicial,
apply (PL3) to the two fixed partitions and use the strict
gap in (PL1). If \(A\) has infinite perimeter, its rescaled
noise deficit tends to infinity and the lift's tends to a
finite value. Both cases give (PL2).
