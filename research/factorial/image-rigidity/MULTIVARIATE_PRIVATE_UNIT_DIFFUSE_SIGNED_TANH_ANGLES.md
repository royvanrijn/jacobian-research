# All degree private unit angular separation for diffuse signed tails

Consider an ordinary two-monomial profile with one private unit
row in each orientation, any shared mirrored core, and any collection
of additional shared rows. Index the additional rows by \(i\) and
write \(\ell_i\ge1\) for the smaller exponent, \(d_i\ne0\) for the
signed gap, and \(M_i=\ell_i+|d_i|/2\). Set
\[
 \begin{gathered}
 \delta=\max_i |d_i|/M_i,\qquad S=\left|\sum_i d_i^3/M_i^2\right|,
 \qquad D=\sum_i d_i^4/M_i^3,\\
 \kappa=\sum_{\rm core\ rows}\frac{d_i^2}{\ell_i+|d_i|+1/2}
                   +\sum_{\rm additional\ rows}\frac{d_i^2}{M_i+1/2},
 \qquad b=1-(1+\kappa/4)^{-1}.
 \end{gathered}
\]
Let \(J\ge4\) be an integer at least the
[absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).

**Theorem.** Suppose \(\kappa>0\) and \(\delta\le1\). The profile
has disjoint consecutive root arguments at every order if
\[
 \begin{gathered}
 R\ge\max(2,4000 b\sqrt J),\qquad H\ge\max(2,30000000 b^2J),\\
                 R(S/4+5D/24)+HD/6\le b^5/2^{16}              \tag{1}
 \end{gathered}
\]
for some real \(R,H\). This proves physical SFC(2) for the
stated class over every characteristic-zero field. Positive
rescalings may be made independently in the two degrees.
At fixed positive reference curvature, the sufficient scales
are \(S=O(J^{-1/2})\) and \(D=O(J^{-1})\), with constants
specified in (1). There is no constraint on the total signed gap
or on the number of diffuse rows.

These quantities use the exact joint signed tanh blocks of
the additional rows. At each slot put
\(v=(n/T'-(n-1)/T)^{-1}=|d_i|/(2M_i)\), and retain the slot's
orientation \(\epsilon=\operatorname{sgn} d_i\).
Then \(\delta=2\max v\), \(S=4|\sum\epsilon v^2|\),
\(D=8\sum v^3\), where the sums count slots. The split into
core and additional rows can be read from the recovered
ordinary Gauss grids. The quantity \(\kappa\) is a rational
lower bound on the reference curvature.
[Reconstruction of ordinary rows](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)

## A linked Gaussian reference for the exact midpoint clocks

The private pair cancels the binomial coefficient. Write \(L_m(v)\)
for the endpoint-normalized log symbol per degree of the additional
rows, with \(v=k/m\), so its contribution to the normalized
log coefficient is \(mL_m(k/m)\). At \(v=1/2+w\) the exact
curvature is
\[
 K_m(v)=L_m''(v)
       =\sum_i d_i^2\int_0^\infty
                   e^{-(M_i+d_i w)u}h(u/m)\,du,
                 \qquad h(t)=\frac{t}{e^t-1}.                \tag{2}
\]
The representation follows by changing variables in the
trigamma series or its positive integral formula. Set
\(\tau_m=K_m(1/2)\). The Gaussian reference multiplier of
degree \(m\) is \(\exp[-\tau_m k(m-k)/(2m)]\), together with the
actual mirrored core and private pair.

The midpoint clock \(\tau_m\) is nondecreasing in \(m\) by
(2), since \(h\) is decreasing. Also \(\tau_m/m\) is
nonincreasing by the original formula
\(\tau_m/m=\sum_i d_i^2\psi'(m M_i+1)\).
Thus with \(N=n+1\), the reference Gaussian variances satisfy
\[
              (n/N)^2\,\tau_n/n\ \le\ \tau_N/N\ \le\ \tau_n/n .
\]
The first inequality already follows from
\(\tau_N\ge\tau_n\). Centered Gaussian densities with these
variances have exactly the folded likelihood-ratio order used
in the mirrored comparison. The folded kernels of the
mirrored core preserve that comparison under convolution.
Accordingly the proofs of the matched phase inequality and the
quantitative reference gap apply with these order-dependent
variances. Both use only the folded comparison for adjacent
orders and the resolvent estimates in each individual order.
The reference gap is greater than \(1/(150000000N^2)\).
[Folded order with a Gaussian factor](MULTIVARIATE_PRIVATE_UNIT_MIRRORED_ROOT_MONOTONICITY.md#adding-the-gaussian-symbol)
[Quantitative mirrored angular reserve](MULTIVARIATE_PRIVATE_UNIT_MIRRORED_UNIFORM_ANGULAR_RESERVE.md)

The log reference coefficient is attenuated by curvature
at least \(\kappa/m\) in the index \(k\). The core rows
satisfy the trigamma bound in the
[near mirror derivative estimate](MULTIVARIATE_PRIVATE_UNIT_NEAR_MIRROR_QUADRATIC_ANGLES.md#the-mirrored-reference-slope).
For the Gaussian, the curvature is \(\tau_m/m\) in \(k\), and
\(\tau_m\ge\sum_i d_i^2/(M_i+1/2)\) for \(m\ge2\) by
\(\psi'(y)\ge1/y\). Hence the reference coefficients satisfy
(2) of that estimate with the present \(b\) and
\(h_k=\min(k,m-k)\). The same slope proof gives
\(\lvert P_m'(\zeta)\rvert\ge mb/8\) at every reference root.
On the disjoint circles of radii
\(\epsilon_m=b^2/(128m)\) about these roots,
\(\lvert P_m\rvert\ge b^3/2048\).

## Even and odd curvature of the residual

The Taylor coefficient of order \(p\) of \(K_m(1/2+w)\)
at zero has absolute value at most
\(\sum_i |d_i|^{p+2}/M_i^{p+1}\).
In (2), use \(0<h\le1\), differentiate, and integrate
\(u^p e^{-M_i u}\). This also proves absolute convergence
for \(|w|\le1/2\) under \(\delta\le1\).
The positive and negative reflections are the even and odd
parts about \(v=1/2\). Summing the even powers from order two
gives
\[
                         \|K_m^+-\tau_m\|_\infty\le D/3 .
\]
For the first odd power, use \(0\le1-h(t)\le t/2\) for
\(t\ge0\); the upper bound is equivalent to
\((t/2)\coth(t/2)\ge1\), which follows from \(\tanh x\le x\).
Thus
\[
               |K_m'(1/2)|\le S+\sum_i |d_i|^3/(mM_i^3)
                             \le S+D/2 .
\]
The last inequality uses \(m\ge2\) and the integer gap
\(|d_i|\ge1\). The odd powers from order three are at most
\(\sum_i |d_i|^5/[8M_i^4(1-\delta^2/4)]\le D/6\)
on the interval. Combining with the first power gives
\[
                        \|K_m^-\|_\infty\le S/2+5D/12 .
\]

Subtract the Gaussian log symbol with curvature \(\tau_m\)
from \(L_m\). The resulting difference has zero endpoints,
as do its even and odd parts. The Dirichlet chord bound
for a function of second derivative bounded by \(C\) is
\(C v(1-v)/2\). Multiplying by \(m\) and using
\(m v(1-v)\le h_k\), the even and odd parts \(e_{m,k},o_{m,k}\)
of the log coefficient error satisfy
\[
         |e_{m,k}|\le h_k D/6,\qquad
         |o_{m,k}|\le h_k(S/4+5D/24).                         \tag{3}
\]
The signed cubic mass appears only in the odd error.
[Matched curvature and higher moments](MULTIVARIATE_MATCHED_TANH_ANGULAR_HIERARCHY.md)

## Angular displacement from the two errors

Multiply reference coefficient \(k\) by
\(\exp[t o_{m,k}+u e_{m,k}]\). At \((t,u)=(1,1)\) this is the
ordinary normalized polynomial. For \(|t|<R\) and \(|u|<H\),
(1) and (3) give a log coefficient error at most
\(Mh_k\), with \(M\le b^5/2^{16}\).
The curvature-attenuated coefficient sum bounds the
evaluation difference on the reference circles by
\(16M/b^2\le b^3/4096\), using the same geometric series
as in the [near mirror root chart](MULTIVARIATE_PRIVATE_UNIT_NEAR_MIRROR_QUADRATIC_ANGLES.md#an-analytic-deformation-respecting-reversal).
Rouché gives one holomorphic simple root in each disk
with \(\big|\log(z_\zeta(t,u)/\zeta)\big|\le2\epsilon_m\).

For \(u=0\), reversing the polynomial sends \(t\) to \(-t\).
Thus the angular displacements for real \(t\) are the
imaginary parts of the even pieces of the analytic root logs.
By Schwarz the displacement at \(t=1,u=0\) is at most
\(2\epsilon_m/R^2\). At \(t=1\), subtract the root log at \(u=0\)
from its value at general \(u\). This is analytic and bounded
by \(4\epsilon_m\) and vanishes at \(u=0\).
Schwarz gives at most \(4\epsilon_m/H\) further displacement
at \(u=1\).

For \(2\le n<N<J\), the sum of angular errors in both degrees
is less than
\[
                   \frac{b^2}{32nR^2}+\frac{b^2}{16nH}.
\]
Each summand is below half the reference gap since
\(N^2/n\le 3J/2\) and (1) gives
\(225000000<16\cdot16000000\) and
\(225000000<8\cdot30000000\).
This proves the angle separation below the absolute cutoff.
Terminal two follows by strict positivity of the second moment
on the real axis and the real root of the first moment.
Terminals \(N\ge J\) use the absolute theorem.

For an unbalanced example with intermediate curvature, let \(L\)
be even and take \(L/2+1\) increasing rows \((L,L+1)\) and
\(L/2\) decreasing rows \((2L+1,2L)\), together with the private
pair. The total signed shared gap is one, and the shared absolute
curvature tends to \(3/4\). For this family \(S=O(1/L)\) and
\(D=O(1/L^2)\), so (1) holds for all sufficiently large even
\(L\) at each admissible cutoff \(J\).
