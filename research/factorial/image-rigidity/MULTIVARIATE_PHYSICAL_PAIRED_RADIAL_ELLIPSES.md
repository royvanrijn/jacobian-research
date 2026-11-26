# Physical radial separation by reflected coefficient pairs

For a factorial binomial, the consecutive moment polynomials share the
same physical ratio of monomial coefficients. A circle certificate exploits
this constraint directly. It applies to two arbitrary orders, every row
dimension, and any private exponents. A sharper endpoint-scale comparison
is available when all unequal rows are shared.

## A circle certificate using rational arithmetic

Let \(Q(z)=\sum_{k=0}^m q_kz^k\) have nonnegative real
coefficients with \(q_0,q_m>0\). At radius \(R>0\), write \(t_k=q_kR^k\).
Assume \(t_0\ne t_m\), and define
\[
 B_Q(R)=
 \sum_{1\le k<m/2}
 \max\left\{\frac{t_k+t_{m-k}}{t_0+t_m},
             \frac{|t_{m-k}-t_k|}{|t_m-t_0|}\right\}
 +\begin{cases}
   t_{m/2}/(t_0+t_m),&m\ \text{even},\\
   0,&m\ \text{odd}.
  \end{cases}
\]
If \(B_Q(R)<1\), the circle \(|z|=R\) contains no root. It encloses
all \(m\) roots when \(t_m>t_0\), and none when \(t_m<t_0\), counting
multiplicities.

Indeed, multiply \(Q(Re^{i\theta})\) by \(e^{-im\theta/2}\) and
divide the real and imaginary parts by \(t_0+t_m\) and \(t_m-t_0\),
respectively. The two endpoint terms then contribute the unit vector
\((\cos(m\theta/2),\sin(m\theta/2))\).
A reflected pair contributes
\[
 \left(
  \frac{t_k+t_{m-k}}{t_0+t_m}\cos((m/2-k)\theta),
  \frac{t_{m-k}-t_k}{t_m-t_0}\sin((m/2-k)\theta)
 \right).
\]
Its Euclidean norm is at most the corresponding maximum in \(B_Q\).
The middle term in even degree is entirely real. Scaling the interior
coefficients by a common factor in \([0,1]\) therefore keeps the circle
free of zeros. The argument principle preserves the root count from
\(q_0+q_mz^m\) along this homotopy, giving the claim.

For the ordinary factorial binomial with exponent columns \(a,b\), use
\[
 q_{m,k}=\binom mk\prod_i((m-k)a_i+kb_i)! .
\]
In particular, a certificate at orders \(n,N=n+1\) consists of a single
positive rational radius \(R\) such that
\[
 B_{Q_n}(R)<1,\qquad B_{Q_N}(R)<1,\qquad
 (q_{n,n}R^n-q_{n,0})(q_{N,N}R^N-q_{N,0})<0.
\]
Every comparison is rational. The roots lie on opposite sides of the
same circle in the variable \(z\), so the resultant is nonzero. For
example, when only the constant and leading coefficients of \(Q_m\)
are initially used to choose test radii, the search interval is the
open interval between their two endpoint radii. The certificate applies
over any characteristic-zero field because the moment polynomials have
rational coefficients.

## Endpoint-normalized root bounds

Set \(G_m=(q_{m,m}/q_{m,0})^{1/m}\) and
\(P_m(w)=Q_m(w/G_m)/q_{m,0}\). Write \(p_k\) for its coefficients.
For a reflected pair put
\[
 g_k=\sqrt{p_kp_{m-k}},\qquad
 s_k=\tfrac12\log(p_{m-k}/p_k),\qquad \nu_k=m/2-k .
\]
When \(x\ne0\), the pair contributes to \(B_{P_m}(e^x)\) by
\[
 g_k\max\left\{
 \frac{\cosh(s_k+\nu_k x)}{\cosh(mx/2)},
 \frac{|\sinh(s_k+\nu_k x)|}{|\sinh(mx/2)|}
 \right\}.
\]
The middle term contributes \(p_{m/2}/(2\cosh(mx/2))\) in even
degree. Thus if \(x_m^-<0<x_m^+\) both pass the certificate,
all normalized log root moduli lie in \((x_m^-,x_m^+)\).
A physical common root requires, with \(C=\log(G_N/G_n)\),
\[
                        x_N^- -x_n^+ < C < x_N^+-x_n^- .
\]
In the frequently useful case
\[
 \sum_{1\le k<m/2}g_k
       +\mathbf 1_{2\mid m}\,p_{m/2}/2 < 1,
 \qquad
 b_m=\max_{1\le k<m/2}|s_k|/k ,
\]
all normalized log root moduli belong to \([-b_m,b_m]\), using
\(b_m=0\) for an empty sum of pairs. This follows by applying the
certificate at \(x=\pm b\) for every \(b>b_m\), since
\(|s_k+\nu_k x|\le mb/2\). When \(b_m>0\) the same argument applies
at \(b=b_m\), giving strict bounds. The symmetric case \(b_m=0\)
recovers the strict real reciprocal case of
[Lakatos and Losonczi, Theorem 1(ii)](https://publi.math.unideb.hu/paper/985/download/10_5486_PMD_2004_3250.pdf).
The criterion can also certify exclusion at a nonzero radius when the
unweighted sum of interior terms is too large, because each reflected
pair is divided separately by the two endpoint axes.

## A signed estimate for the physical shift

Assume all unequal rows are shared. For each row let
\(\ell=\min(a_i,b_i)\), \(d=|b_i-a_i|\), and
\(\epsilon=\operatorname{sgn}(b_i-a_i)\). Delete equal rows. Let
\(\Delta=\sum_i\epsilon_i d_i\). Put
\[
 L_\pm=\sum_{\epsilon_i=\pm1}\log(1+d_i/\ell_i),\qquad
 H_\pm=\sum_{\epsilon_i=\pm1}
                    \left(\frac1{\ell_i}-\frac1{\ell_i+d_i}\right).
\]
Then
\[
 \begin{split}
 C&=\Delta\log(N/n)-\frac{L_+-L_-}{2nN}+E,\\
 -\frac{1/n^2-1/N^2}{12}\,H_-
     &\le E\le
       \frac{1/n^2-1/N^2}{12}\,H_+ .
 \end{split}
\]
This estimate retains cancellation in the signed logarithmic term.

To prove the bound, for \(u>0\) use the trigamma series to write
\[
 \psi'(1+u)=\int_0^\infty e^{-ut}\frac{t}{e^t-1}\,dt .
\]
For \(t>0\) one has
\[
           1-\frac t2\le\frac{t}{e^t-1}\le1-\frac t2+\frac{t^2}{12}.
\]
Indeed \(t/(e^t-1)=x\coth x-x\) for \(x=t/2\). The lower estimate
uses \(x\cosh x-\sinh x\ge0\) by its power series. For the upper
estimate, the coefficient of \(x^{2k+1}\) in
\((1+x^2/3)\sinh x-x\cosh x\) is
\(4k(k-1)/(3(2k+1)!)\), which is nonnegative for \(k\ge0\).
Integration gives
\[
          \frac1{2u}-\frac1{6u^2}
                    \le 1-u\psi'(1+u)\le\frac1{2u}.
\]
Also \(1-u\psi'(1+u)\ge0\), by the decreasing-function upper integral
for the trigamma series. An increasing row contributes
\(\int_\ell^{\ell+d}\psi(1+jx)\,dx\) to \(\log G_j\).
Its derivative is \(d/j\) minus
\[
 \frac1j\int_\ell^{\ell+d}(1-jx\psi'(1+jx))\,dx .
\]
Thus its integrated correction to \(d\log(N/n)\) is
\(-\log(1+d/\ell)/(2nN)\) plus a nonnegative term at most
\((1/n^2-1/N^2)(1/\ell-1/(\ell+d))/12\).
A decreasing row reverses both signs. This proves the interval.

There is also a uniform sharpening of the
[physical grid at fixed reflected mean](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md#the-physical-grid-at-fixed-reflected-mean).
For \(n\ge2\), define \(L_n^*=\log(2n/(n+1))\). Let \(U_i\) denote
the outer loss of a row at order \(n\), and \(U_\pm\) their sums by
orientation. Then
\[
 \Delta\log(N/n)-\frac{\log 2}{nN L_n^*}\,U_+
 \ \le C\le\
 \Delta\log(N/n)+\frac{\log 2}{nN L_n^*}\,U_- .
\]
For a row with \(y=d/\ell\), the definition of \(U_i\) gives
\[
 U_i\ge \frac d2 \log\frac{n(1+y)}{n+y}.
\]
The ratio \(\log(1+y)/\log(n(1+y)/(n+y))\) increases with \(y>0\).
In fact, at \(v=\log(1+y)\) its denominator is
\(v-\log((e^v+n-1)/n)\). The second term is convex and vanishes at
\(v=0\), so its quotient by \(v\) increases. Since \(y\le d\), this
gives
\[
 \frac{\log(1+y)}{2U_i}
 \le\frac{\log(1+d)}{d\log(n(1+d)/(n+d))}
 \le\frac{\log 2}{L_n^*}.
\]
For the last inequality both \(\log(1+d)/d\) and the reciprocal of
\(\log(n(1+d)/(n+d))\) decrease for \(d\ge1\). The row correction
has its stated sign and magnitude at most
\(\log(1+d/\ell)/(2nN)\). Summing over orientations proves the bound.
The uniform factor \(\log 2\) is sharp for this form of estimate over
all orders. For one row \((1,2)\), \(U_i=L_n^*/2\) and the signed
expansion gives \(nN(\log(N/n)-C)\to(\log 2)/2\) as \(n\to\infty\).


## Exact example and scope

For the one-row profile \((a,b)=(1,4)\), take \(n=2\) and
\(R=1/300\). The coefficient lists are
\[
             (2,240,40320),\qquad (6,2160,1088640,479001600).
\]
The endpoint differences have opposite signs. The certificate values
are \(50/153\) and \(2010/2473\), respectively. At order three the
ratio of the sum of the individual interior magnitudes to the minimum
endpoint magnitude is \(2010/1223>1\), so the pair estimate is needed
to certify this particular circle by endpoint dominance.

A profile whose exponent columns are exchanged by a row permutation
has both endpoint scales equal to \(1\) at every order. Some adjacent
moments of such profiles share root moduli. Already for rows
\((1,2),(2,1)\), the first two moment polynomials are proportional to
\(1+z\) and \(2+3z+2z^2\), respectively. All their roots have modulus
\(1\), and the quadratic is nonzero at \(-1\). This gives an explicit
case requiring information about arguments for a physical-root proof.

