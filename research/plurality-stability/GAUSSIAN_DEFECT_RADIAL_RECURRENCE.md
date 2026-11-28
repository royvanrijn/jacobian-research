# Radial recurrence of defects with exact reference score planes

Let \(k\ge3\), let \(T\) be the centered regular simplex in
its active space \(\mathbb R^{k-1}\), and let \(A\) be an ordinary
Gaussian partition obeying the winner rule for its own
adjusted scores \(V_i=T_\rho 1_{A_i}-\beta_i\).
Assume every adjusted contrast vanishes on its full
reference plane \(H_{ij}\). A discrepancy from \(T\)
has a definite amount of Euclidean mass near its
inward Gaussian source point. In particular a
noncompact discrepancy, under a small pressure bound,
has a thick sequence of defects along a single geometric
ray. A sufficiently clean thick annulus already forces
an exact simplex exterior, even with arbitrary changes
inside the annulus.

## Each defect requires a defect at the inward source scale

Fix a compact band \(J=[a,b]\subset(0,1)\). Write
\(D=\{x:A(x)\ne T(x)\}\), understood up to null sets,
\(\sigma=\sqrt{1-\rho^2}\) and
\(\sigma_b=\sqrt{1-b^2}\). The saturated simplex
margin gives \(\vartheta>0\), depending only on \(k,J\),
such that
\[
 (T_\rho1_{T_i}-T_\rho1_{T_j})(x)
                  \ \ge\ \vartheta\min(t,1),
 \qquad x\in T_i,\quad t=\operatorname{dist}(x,H_{ij}).
\]
[The saturated simplex margin (FG1)](GAUSSIAN_FINITE_SIMPLEX_GERM_PROPAGATION.md#a-comparison-in-a-growing-cylinder)

Assume for the moment
\(\max_{ij}|\beta_i-\beta_j|\le\vartheta/2\).
There are \(C,\delta>0\), uniform over \(\rho\in J\),
for which
\[
       |D\cap B_C(\rho x)|\ \ge\ \delta
                       \qquad\text{for a.e. }x\in D.       \tag{DR1}
\]
Here \(|\cdot|\) means Euclidean volume. To prove this,
suppose \(T(x)=i\ne j=A(x)\) at a point where the winner
rule holds and the reference winner is strict. Define
\[
 f_{ij}=1_{A_i}-1_{A_j}-1_{T_i}+1_{T_j},\qquad
 d_{ij}=T_\rho f_{ij}-(\beta_i-\beta_j).
\]
Then \(|f_{ij}|\le2\,1_D\), \(d_{ij}=0\) on \(H_{ij}\),
and \(d_{ij}(x)\le-\vartheta\min(t,1)\).
If \(t\ge1\), this forces
\(|T_\rho f_{ij}(x)|\ge\vartheta/2\). If \(0<t<1\),
the mean value theorem along the normal segment from
\(H_{ij}\) to \(x\) gives
\(|\partial_\nu T_\rho f_{ij}(x')|\ge\vartheta\)
at a point \(x'\) with \(|x-x'|\le1\), where
\(\nu\) is the unit normal to the plane.

For either test, use \(y=\rho x'+\sigma Z\) with \(Z\)
standard Gaussian and \(x'=x\) in the first case.
The conditional score and its derivative are integrals
of \(f_{ij}(y)\) and
\((\rho/\sigma)(Z\cdot\nu)f_{ij}(y)\), respectively.
The contribution of \(|Z|>M\) is uniformly small as
\(M\to\infty\) in both tests, using \(b/\sigma_b\) for
the prefactor. In the remaining integral the conditional
density is bounded by \((2\pi\sigma_b^2)^{-(k-1)/2}\)
and \(|Z\cdot\nu|\le M\). Thus a uniform positive lower
bound on either test forces a positive amount of Euclidean
volume of \(D\) in \(B_M(\rho x')\). Enlarging this ball
by one gives (DR1).

## A single geometric ray supports recurrent defects

If \(D\) is not essentially bounded under the pressure
hypothesis, there are a unit vector \(\omega\), a radius
\(r_*\in[1,\rho^{-1}]\) and \(C_*>0\), with \(C_*\) depending
only on \(k,J\), such that
\[
       |D\cap B_{C_*}(r_*\rho^{-q}\omega)|\ge\delta
                         \qquad(q=0,1,2,\ldots).            \tag{DR2}
\]
Indeed start at any typical defect \(x_0\) and use (DR1)
to select typical defects \(x_{\ell+1}\in B_C(\rho x_\ell)\).
Then \(|x_\ell-\rho^\ell x_0|\le C/(1-\rho)\).
In addition, (DR1) gives volume at least \(\delta\)
in the ball of radius \(C\) at \(\rho x_{\ell-1}\).
Choose a sequence of starting defects \(x_0^{(n)}\) with
radii tending to infinity, and write
\[
     |x_0^{(n)}|=r_n\rho^{-m_n},\qquad 1\le r_n<\rho^{-1}.
\]
Take a subsequence with \(r_n\to r_*\) and
\(x_0^{(n)}/|x_0^{(n)}|\to\omega\). For each fixed
\(q\ge0\), use step \(\ell=m_n-q\) for large \(n\).
Its source ball lies within distance
\(2C/(1-b)+1\) of \(r_*\rho^{-q}\omega\) for all
sufficiently large \(n\). The same enlarged radius is
valid for every \(q\). This proves (DR2).
In particular no noncompact discrepancy can have finite
Euclidean volume, and in every sufficiently large
logarithmic radial window the volume grows at least
linearly in the number of these geometric scales.

## A clean thick annulus forces exterior agreement

Assume the exact adjusted planes, but do not impose
the small pressure bound in advance. There are constants
\(M_0,C_0,\eta>0\), uniform on \(J\), with the following
property. For some \(R\ge1\) and \(S\) with
\[
     S\ \ge\ (R+C_0)/a+C_0 ,
\]
suppose
\[
  |D\cap B_{M_0}(y)|\le\eta
            \qquad\text{whenever } R\le |y|\le S.           \tag{DR3}
\]
Then \(A=T\) outside a compact set, and all prices are
equal. This includes exact agreement on a sufficiently
thick annulus after a fixed enlargement of its endpoints.

To prove this, choose \(M_0\) at least as large as \(C\)
in (DR1) and large enough that the conditional score
tail outside \(B_{M_0}(\rho x)\) is at most
\(\vartheta/4\). At a point on \(H_{ij}\) with
\(R\le|\rho x|\le S\), the equality
\(T_\rho f_{ij}(x)=\beta_i-\beta_j\) and (DR3) bound the
remaining score contribution by \(\vartheta/4\) when
\(\eta\) is small enough. This verifies the pressure
bound for all pairs. Choose also \(\eta<\delta\).
If a typical defect lies at radius larger than \(S/\rho\),
iterate \(x_{\ell+1}\in B_C(\rho x_\ell)\) until the
first step with \(\rho |x_\ell|\le S\). There is such a
step, since the recurrence forces the limiting upper
radius to be at most \(C/(1-\rho)\) and the constants
can be chosen so \(S/\rho\) exceeds it. At the first
step \(\rho |x_\ell|\ge\rho S-\rho C>R\), by the preceding
radius and the assumed bound on \(S\). The lower bound
(DR1) now contradicts (DR3). Thus the discrepancy has
bounded support. In this case \(T_\rho f_{ij}(x)\) tends
to zero along each unbounded reference plane, so the
exact adjusted equality gives \(\beta_i-\beta_j=0\).

In particular the same conclusion holds if the
discrepancy volume in balls of fixed radius tends
uniformly to zero at infinity. It also holds when
\(D\) has finite Euclidean volume. The proof uses
both exact reference planes and the stationarity
winner rule. It leaves open whether an arbitrary
distant maximizing competitor has those planes or
whether defects recurring as in (DR2) can occur.
