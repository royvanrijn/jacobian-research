# Matched second-order angular interlacing at every degree

Let \(n\ge2\) and \(N=n+1\). Match the Gauss factors of an ordinary
fully shared profile as in the
[matched tanh moment lemma](MULTIVARIATE_MATCHED_TANH_MOMENTS.md).
Let \(X=\sum_i x_{n,i}>0\) be the order-\(n\) reciprocal sum.
Here a single index \(i\) ranges over the matched factors, including
all factors within a row. For \(j=n,N\) put
\[
 \xi_{j,i}=x_{j,i}/X,\quad R_j=\sum_i\xi_{j,i},\quad
 V_{j,m}=\sum_i\xi_{j,i}^m,\quad D_j=\sum_i\sigma_i\xi_{j,i}^2,
\]
where \(\sigma_i=1\) for an increasing row and \(-1\) for a reflected
row. All moments here use the common scale \(X\). In the deformation below,
the tanh coordinate is
\[
 T_{j,i}(\lambda)=
 \frac{(j-1)\lambda^2\xi_{j,i}}{2+(j-1)\lambda^2\xi_{j,i}}.
\]
In particular \(\sum_i\sigma_iT_{j,i}^2\) has leading term
\((j-1)^2\lambda^4 D_j/4\). The unsigned tanh moments at powers two and
three likewise start with the corresponding \(V_{j,m}\).

For a fixed list \(\xi_{j,i}\), deform the reciprocals to
\(\lambda^2 \xi_{j,i}\).
If \(iy\) is an imaginary Hermite anchor of degree \(j\) and variance
\(R=R_j\), the signed angular displacement from the negative real
axis, divided by \(\lambda\), has the fixed-degree expansion
\[
 \Theta_{j,y}(\lambda)
      =y+\lambda^2 y B_j(y)+O(\lambda^4),
 \tag{1}
\]
where the positive upper-half-plane convention gives \(y>0\), and
\[
 \begin{aligned}
 B_j(y)={}&-\frac{jR}{24}-\frac{(j-1)V_{j,2}}{4R}
       +\frac{V_{j,3}}{12R^2}\left(2j-3-\frac{y^2}{R}\right)\\
 &\hspace{14mm}
       -\frac{D_j^2}{36R^3}\left(\frac72j-5-\frac{2y^2}{R}\right).
 \end{aligned}
 \tag{2}
\]

**Theorem.** If \(nX\le1/100\), the positive heights of the two
second-order models
\[
                       y\ \longmapsto\ y(1+X B_j(y))
 \tag{3}
\]
strictly interlace, with the order-\(N\) height first in decreasing
order. Conjugation gives the negative heights and the odd degree
has the zero height. This is a theorem for the model (3). A uniform,
linked error estimate for the exact angular heights on the scale of
the adjacent Hermite gaps is open in this argument.

## Exact angular coefficient

Write \(k-j/2=v\). In an increasing factor the log of the normalized
coefficient with index \(k\) is
\[
 \sum_{h=0}^{k-1}\log(1+h\lambda^2\xi)
       -\frac{k}{j}\sum_{h=0}^{j-1}\log(1+h\lambda^2\xi).
\]
For a reflected factor replace \(v\) by \(-v\). Let \(F_k(\lambda^2)\)
be the product of these normalized coefficients. Use the exponential
variable \(t=-e^{-\lambda z}\), so that the imaginary part of \(z\)
is exactly the divided angle near the negative axis. The transformed
polynomial, with a nonzero exponential factor removed, is
\[
 C_j(z,\lambda)
    =\lambda^{-j}\sum_{k=0}^j(-1)^k\binom jk F_k(\lambda^2)
                                             e^{-\lambda vz}.
 \tag{4}
\]
The logarithmic coefficient of \(\lambda^{2r}\) in \(F_k\) has
degree at most \(r+1\) in \(v\). The centered alternating moments come
from
\[
       \sum_{k=0}^j (-1)^k\binom jk e^{v s}
                 =(-2\sinh(s/2))^j
                 =(-1)^j s^j(1+js^2/24+O(s^4)).
 \tag{5}
\]
All negative powers of \(\lambda\) in (4) cancel. At degree \(j\)
in the alternating sum, the leading moment is \((-1)^j j!\).
At degrees \(j+1\) and \(j+2\) the moments are \(0\) and
\((-1)^j(j+2)!\,j/24\), respectively.

Set \(w=\lambda v\), and abbreviate
\(\beta=-(j-1)V_{j,2}/4\), \(\gamma=V_{j,3}/12\), \(d=D_j/6\).
In the bookkeeping by powers of \(\lambda\) and \(w\), the logarithm
of \(F_k\) is, through degree two in \(\lambda\),
\[
   Rw^2/2-\lambda d w^3
          +\lambda^2\bigl(-Rj^2/8+\beta w^2+\gamma w^4\bigr).
\]
These terms follow from the first three power sums in \(h\).
For the second log term, the centered cubic changes sign under
reflection and has coefficient \(-D_j/6\). Its centered quadratic
coefficient is \(-(j-1)V_{j,2}/4\) at both orientations.
For \(r\ge4\), the degree bound \(r+1\) gives at least
\(\lambda^{2r-(r+1)}\), which cannot affect the displayed orders.
Writing \(H_j(z)=j![s^j]\exp(zs+Rs^2/2)\), equations (4) and (5) give
\[
 \begin{aligned}
 C_j(z,\lambda)={}&H_j(z)+\lambda d H_j'''(z)\\
 &+\lambda^2\left(
 -\frac{Rj^2}{8}H_j(z)+\beta H_j''(z)+\gamma H_j^{(4)}(z)
 +\frac{d^2}{2}H_j^{(6)}(z)+\frac{j}{24}H_{j+2}(z)\right)
 +O(\lambda^3).
 \end{aligned}
 \tag{6}
\]
The \(H_{j+2}\) term uses the same variance \(R\).

At a zero \(z_0=iy\) of \(H_j\), put
\(K_m=H_j^{(m)}(z_0)/H_j'(z_0)\).
The zeros are simple by the tridiagonal Hermite recurrence.
Differentiating \(R H_j''+zH_j'-jH_j=0\) gives
\(RK_{m+2}+z_0K_{m+1}+(m-j)K_m=0\).
Also \(H_{j+2}(z_0)/H_j'(z_0)=Rz_0\).
The first displacement of the zero is \(-dK_3\), which is real.
The part of its second displacement proportional to \(d^2\) is
\[
 -d^2(K_6/2-K_3K_4+K_3^2K_2/2)
     =-d^2\left((7j/2-5)\frac{z_0}{R^3}
                          +2\frac{z_0^3}{R^4}\right).
\]
The other terms are \(-\beta K_2-\gamma K_4-jRz_0/24\).
Taking imaginary parts gives (2).
The real coefficients and
\(C_j(z,-\lambda)=(-1)^jC_j(-z,\lambda)\) imply that the
analytic root branch obeys \(z(-\lambda)=-\overline{z(\lambda)}\)
for real \(\lambda\) near zero. Its imaginary part is even, giving
the remainder in (1) at fixed degree.

## Comparison across the two orders

The exact reciprocal link satisfies
\[
        \frac nN\le\frac{x_{N,i}}{x_{n,i}}
              =\frac nN+\frac{h_i}{N(Na_i+h_i)}
              \le\frac nN+\frac XN .
\]
In particular \(R_n=1\), \(n/N\le R_N\le(n+X)/N\), and \(R_N\le1\).
Write \(a_m=V_{j,m}/R_j^m\) and \(b=D_j/R_j^2\) within an order.
Then
\[
 B_j(y)=R_j(j c+c_0)+y^2c_1,\quad
 c=-1/24-a_2/4+a_3/6-7b^2/72,\quad
 c_0=a_2/4-a_3/4+5b^2/36,\quad
 c_1=-a_3/12+b^2/18.
\]
The ratios of the normalized factor weights \(\xi_i/R_j\) at the
two orders differ from one by at most \(X/n\).
Using \(X\le1\), the total changes of \(a_2,a_3,b\) are at most
\(8X/n,12X/n,8X/n\), respectively. Indeed the difference of the
\(m\)-th powers of two weights is at most
\(m(1+X/n)^{m-1}(X/n)\) times the lower-order weight to power \(m\),
and these lower-order powers sum to at most one.
It follows that the changes of \(c,c_0,c_1\) are at most
\((50/9)X/n,(65/9)X/n,(17/9)X/n\).
Using \(|c|\le5/9\), \(|c_0|\le23/36\),
\(0\le NR_N-n\le X\), and \(0\le1-R_N\le1/N\) gives
\[
             |B_N(y)-B_n(y)|\le25X+\frac1n
             \qquad (y^2\le4N).
 \tag{7}
\]
In this last bound \(4N\le6n\).
On the same interval, (2) gives
\[
 |B_j(y)+yB_j'(y)|\le 5j/9+(5/12)y^2\le(20/9)N.
 \tag{8}
\]

Let \(Y_{j,k}\) be the positive Hermite zeros in decreasing order.
Their Jacobi matrices have consecutive off-diagonal entries
\(\sqrt{R_j\,\ell}\) for \(1\le\ell\le j-1\).
Each height is at most \(2\sqrt N\) by the row sums.
The leading \(n\)-by-\(n\) block of the order-\(N\) matrix is
\(\sqrt{R_N}\) times the order-\(n\) matrix, and Cauchy interlacing
gives \(Y_{N,k+1}\le\sqrt{R_N}Y_{n,k}\) whenever the two heights exist.
The trailing block has consecutive ratios
\(\sqrt{R_N(\ell+1)/\ell}\) to the order-\(n\) off-diagonals.
They decrease and are at least \(c_*=(1-1/n^2)^{-1/2}\).

For completeness, write \(r_\ell\) for these ratios divided by \(c_*\).
Choose a diagonal \(D_*\) by \(D_{*,n}=1\) and
\(D_{*,\ell}D_{*,\ell+1}=r_\ell\).
The logarithms of its entries are alternating sums of decreasing
nonnegative numbers, hence each diagonal entry is at least one.
The trailing block is \(c_*D_*J_nD_*\).
Under the bipartite decomposition of the path matrix \(J_n\), the
positive eigenvalues are the singular values of its off-diagonal
rectangular block. Multiplication of that block on either side by
a diagonal matrix with entries at least one increases each singular
value in the weak sense. This follows from the positive-semidefinite
comparison of the appropriate Gram matrices. Cauchy interlacing
now gives \(Y_{N,k}\ge c_*Y_{n,k}\).

Both adjacent gaps thus admit relative lower bounds:
\[
 Y_{N,k}-Y_{n,k}\ge\frac{Y_{n,k}}{2n^2},\qquad
 Y_{n,k}-Y_{N,k+1}\ge\frac{Y_{n,k}}{4n^2}.
\]
For the second bound, use \(1-\sqrt{R_N}\ge(1-X)/(2N)\ge1/(4n^2)\)
when \(nX\le1/100\).
By (8) the maps (3) are increasing with slope at least \(1/2\).
Their difference at the same height is at most
\(X(25X+1/n)y\), which is at most \(y/(80n^2)\) under \(nX\le1/100\).
It is smaller than half either displayed gap, so the same signs
hold after the correction (3).

The [exact algebra checker](verify_multivariate_matched_second_order_angles.py)
reconstructs the first three coefficients of (4) directly from
truncated Gauss products for several rational factor lists and degrees.
It also checks the angular formula by exact polynomial reduction
modulo each corresponding Hermite polynomial. Replay with

~~~sh
python3 -O research/factorial/image-rigidity/verify_multivariate_matched_second_order_angles.py
~~~
