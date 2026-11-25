# Consecutive factorial moments at uniformly large curvature

Let \(a,b\) be distinct nonnegative integer exponent columns in any
finite dimension. For a positive terminal \(N\), set
\[
 S=\sum_{i:\,\max(a_i,b_i)>0}
       \frac{(b_i-a_i)^2}{\max(a_i,b_i)} .
\]
For every \(N\ge256\) with
\[
                         S\ge 2\log N+5,                   \tag{1}
\]
the ordinary factorial moment polynomials of orders \(N-1,N\)
are coprime. The criterion is uniform in the dimension, the
exponents, and the characteristic-zero coefficient field. In
particular it also applies when \(\sum_i a_i=\sum_i b_i\).

## An endpoint coefficient estimate

For \(m=N-1,N\), define
\[
 C_{m,k}=\prod_i((m-k)a_i+kb_i)!,\qquad
 G_m=(C_{m,m}/C_{m,0})^{1/m},\qquad
 T_m(z)=\sum_{k=0}^m c_{m,k}z^k,\quad
 c_{m,k}=\binom mk\frac{C_{m,k}}{C_{m,0}G_m^k}.
\]
The physical scalar ratio \(t\) satisfies \(z=G_m t\), with
\(G_m>0\), and both endpoint coefficients of \(T_m\) equal one.

The [discrete curvature estimate](MULTIVARIATE_FIXED_TERMINAL_CURVATURE.md#large-curvature-and-the-endpoint-radii)
gives
\[
 0<c_{m,k}\le \binom mk
       \exp\left(-\frac{S k(m-k)}{2m}\right).               \tag{2}
\]
For completeness, in a single row with \(0\le a_i<b_i\) and
\(d=b_i-a_i\), the second difference of the logarithm of its
factorial is
\[
 \sum_{j=1}^d \log\left(1+\frac d{ma_i+(k-1)d+j}\right)
   \ge \frac{d^2}{mb_i}\qquad (1\le k<m).
\]
The bound uses \(\log(1+x)\ge x/(1+x)\) and
\(ma_i+kd+j\le mb_i\). Reversing the row keeps the bound
with \(b_i\) replaced by \(\max(a_i,b_i)\). After summation,
subtracting a quadratic of second difference \(S/m\)
leaves a convex sequence. Comparing it with its endpoint
chord yields (2).

Under (1), both orders satisfy the sharper geometric estimate
\[
                  c_{m,k}\le 8^{-\min(k,m-k)}
                  \qquad (0\le k\le m).                   \tag{3}
\]
Put \(A=\log N+5/2\). The function \(2(\log x+5/2)/(x-1)\)
decreases for \(x\ge256\). Since \(\log 2<7/10\), this implies
\(2A/m<1/12\). By symmetry of the bound (2), it suffices to
treat \(k\le m/2\). For \(k=1\) the bound is at most
\(\exp(-5/2+A/m)<1/8\). For \(k\ge2\), the elementary
estimate \(\binom mk\le (em/k)^k\) bounds the logarithm of
(2), divided by \(k\), by
\[
 g(k/m),\qquad g(w)=1-\log w-A(1-w).
\]
This function is convex on \([2/m,1/2]\). Its endpoint values
satisfy
\[
 g(2/m)\le -\frac{17}{12}-\log2<-3\log2,\qquad
 g(1/2)\le-\frac14-3\log2.
\]
Thus it is below \(-\log8\) on the entire interval.

## Angular slope and the common-root contradiction

For either order write \(h=\lfloor m/2\rfloor\) and
\[
 L_m(w)=\sum_{k=0}^{h}c_{m,k}w^k,\qquad
 U_m(w)=\sum_{j=0}^{m-h-1}c_{m,m-j}w^j,\qquad
 T_m(z)=L_m(z)+z^m U_m(1/z).
\]
Set \(R=e^{1/m}\). For either half \(H=L_m,U_m\), (3) and
\(e^{1/m}<8/7\) imply, throughout \(|w|\le R\),
\[
 |H(w)-1|<\sum_{j\ge1}7^{-j}=\frac16,\qquad
 |wH'(w)|<\sum_{j\ge1}j7^{-j}=\frac7{36},\qquad
 \left|\frac{wH'(w)}{H(w)}\right|<\frac7{30}.              \tag{4}
\]
In particular each half has positive real part on the disk,
so its argument is continuous there and equals zero on the
positive real radius. On \(|z|=R\), \(z^m U_m(1/z)\)
dominates \(L_m(z)\) in modulus by (4) and \(R^m=e\).
The dominating polynomial has all its zeros in \(|z|<R\),
since \(U_m\) is zero-free for \(|w|\le R\). Rouché's
theorem puts all zeros of \(T_m\) inside that disk. Applying
the same domination to the reciprocal polynomial, with
\(L_m,U_m\) exchanged, puts all root moduli above \(R^{-1}\).

Let \(z=r e^{i\theta}\) be a root with \(0<\theta\le\pi\).
The positive coefficients exclude \(\theta=0\). By (4),
the derivatives of \(\arg L_m(r e^{i\theta})\) and
\(\arg U_m(r^{-1}e^{-i\theta})\) with respect to \(\theta\)
both have absolute value less than \(7/30\). Integrating from
zero and writing
\[
 \delta_m=\arg L_m(r e^{i\theta})-\arg U_m(r^{-1}e^{-i\theta})
\]
gives \(|\delta_m|<7\theta/15\). The exact root equation
gives \(m\theta=(2j_m+1)\pi+\delta_m\) for some integer \(j_m\).

If \(t\) were a common root of the physical polynomials,
conjugation would give \(0<\arg t\le\pi\). The two normalized
roots have the same angle \(\theta=\arg t\), since \(G_m>0\).
Subtracting their phase equations would give
\[
 2\pi(j_N-j_{N-1})=\theta-\delta_N+\delta_{N-1}
       \ \in\ (\theta/15,29\theta/15)\ \subset\ (0,2\pi),
\]
which is impossible. Coprimality over \(\mathbb C\) gives
rational Bézout coefficients for the two ordinary moment
polynomials, whose original coefficients are integers. The
same Bézout identity works over every characteristic-zero field.
