# Higher label separation by linked endpoint noise and integer drift

Consider an ordinary two-monomial profile with one private unit
row per orientation and arbitrary shared rows. For \(N=n+1\),
any common physical nonreal root would have the same positive
angle-box label \(k\) in the two degrees after conjugation if
needed.

**Theorem.** Common physical roots are excluded in the following
ranges, for all such profiles:
\[
                   n\ge56,\quad k\ge2,
       \qquad\text{or}\qquad n\ge24,\quad k\ge3 .
\]
The proof uses the exact adjacent Beta link and has no eventual
cutoff. In these ranges a same-angle contact of the normalized
roots would be too close in radius to absorb any nonzero
integer physical drift. The balanced Euclidean theorem covers
zero drift.

For high shared curvature the
[every-degree curvature theorem](MULTIVARIATE_BOTH_PRIVATE_UNIT_CURVATURE_15_OVER_8.md)
already covers the complete angle range. We therefore take
\(0<S<15/8\) where \(S=\sum_i d_i^2/(\ell_i+|d_i|)\)
on shared rows, \(d_i=b_i-a_i\), and \(\ell_i=\min(a_i,b_i)\).
Put \(V=\sum_i d_i^2/\ell_i\) and \(\Delta=\sum_i d_i\).
Exact ordinary row packing gives \(V<61/12\) and
\(|d_i|/\ell_i\le2\).
[Absolute imbalance bound and the packing](MULTIVARIATE_PRIVATE_UNIT_IMBALANCE_TWENTY_FIVE.md#shared-curvature-and-the-endpoint-laws)
[Angle boxes and matching labels](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md)
[Balanced Euclidean theorem](MULTIVARIATE_BOTH_PRIVATE_EUCLIDEAN_DESCENT.md)

## Exact signed log noise

The private unit factorials cancel the binomial coefficient.
Write \(G_m>0\) for the endpoint-normalizing argument scale.
The normalized coefficients are the moments of a positive
endpoint variable \(T_m\) with \(\mathbb E T_m^m=1\).
In the given orientation it is the product of \(X_i^{d_i}\)
divided by \(G_m\), with independent unit-scale Gamma variables
\(X_i\) of shape \(m a_i+1\).
The reverse law uses the opposite exponents and scale \(1/G_m\).

The beta-gamma product identity gives an independent link in
distribution \(\log T_n\stackrel{\rm law}=\log T_N+Z\), where
\[
     Z=\log(G_N/G_n)+\sum_i d_i\log Y_i,\qquad
                   Y_i\sim {\rm Beta}(na_i+1,a_i)
\]
independently. Indeed multiplying an independent Gamma variable
of shape \(Na_i+1\) by \(Y_i\) gives shape \(na_i+1\);
this follows by changing the variables of two independent
Gamma inputs to their sum and fraction. Both endpoint
orientations have the same estimates below.
The link is the law version of the
[positive Mellin coefficient quotient](MULTIVARIATE_BOTH_PRIVATE_LINKED_COEFFICIENT_ORDER.md#positive-mellin-representation-for-the-quotient).

Put \(J(x)=\psi(Nx+1)-\psi(nx+1)\). The exact endpoint
normalization is
\(\log(G_N/G_n)=\sum_i d_i\int_0^1J(a_i+d_i u)\,du\).
The beta mean gives
\[
 \mathbb E Z=\sum_i d_i\int_0^1[J(a_i+d_i u)-J(a_i)]\,du.
\]
The trigamma Laplace formula gives
\[
 0\le J'(x)=\int_0^\infty e^{-xy}\big[h(y/N)-h(y/n)\big]\,dy
                  \le1/(2nNx^2),\qquad h(v)=v/(e^v-1).
\]
Here \(h\) is decreasing and \(-h'\le1/2\), by its convexity
and the derivative at zero. Thus
\[
                        0\le\mathbb E Z\le V/(4nN).           \tag{1}
\]
[Exact trigamma transport kernel](MULTIVARIATE_BOTH_PRIVATE_LINKED_COEFFICIENT_ORDER.md#comparing-a-fixed-index-at-adjacent-degrees)

Let \(M(t)=\mathbb E e^{tZ}\) and \(\Lambda=\log M\).
They are finite for \(|t|<n/2\), since \(na_i+d_i t>0\)
using \(|d_i|/a_i\le2\).
Twice differentiating the beta integral and integrating the
trigamma difference in real order gives
\[
 \begin{aligned}
 \Lambda''(t)
   &=\sum_i d_i^2[\psi'(na_i+d_i t+1)-\psi'(Na_i+d_i t+1)]\\
   &\le \frac{V}{nN(1-2|t|/n)^2}.
 \end{aligned}                                                \tag{2}
\]
For the last bound use
\(|\psi''(1+x)|\le x^{-2}\) for \(x>0\) and
\(ma_i+d_i t\ge ma_i(1-2|t|/n)\) for \(n\le m\le N\).
In particular
\[
 |\Lambda'(t)|\le V/(4nN)
              +\frac{V|t|}{nN(1-2|t|/n)} .
                                                               \tag{3}
\]

## Resolvent variation at any horizontal position

Let \(0<\theta\le\pi/2\) and \(F(s)=(1-e^{s+i\theta})^{-1}\).
On the real \(s\)-line,
\[
 |F'|=\operatorname{Im}F/\sin\theta,\qquad
 |F''|\le\operatorname{Im}F/[\sin\theta\sin(\theta/2)],\qquad
 |(\log\operatorname{Im}F)'|\le1/\sin\theta.
\]
The first two follow by differentiating and using
\((1+e^s)/|1-e^{s+i\theta}|\le1/\sin(\theta/2)\).
For the third, use \(\operatorname{Im}F=
\sin\theta/[2(\cosh s-\cos\theta)]\) and
\((\cosh s-\cos\theta)^2-\sin^2\theta\sinh^2s
=(\cos\theta\cosh s-1)^2\).
In particular \(|\partial_x\log\mathbb E F(x+\log T_m)|
\le1/\sin\theta\), with no restriction on the real shift \(x\).

Put \(b=n\sin\theta\), \(\alpha=1/\sin\theta\), and take \(b>2\).
For a common real shift \(x\), write \(E_n=\mathbb E F(x+\log T_n)\)
and \(E_N=\mathbb E F(x+\log T_N)\).
Taylor expand \(F(x+\log T_N+Z)\) in \(Z\) independently of \(T_N\).
The first variation costs
\(|\mathbb E Z|\operatorname{Im}E_N/\sin\theta\).
Along the segment from zero to \(Z\), the second derivative
is bounded by \(\operatorname{Im}F(x+\log T_N)e^{\alpha|Z|}
/[\sin\theta\sin(\theta/2)]\).
Thus \(|E_n-E_N|\le\varepsilon |E_N|\), where
\[
 \varepsilon=\frac{\mathbb E Z}{\sin\theta}
         +\frac{M''(\alpha)+M''(-\alpha)}{\sin^2\theta}.
                                                               \tag{4}
\]
We used \(2\sin(\theta/2)\ge\sin\theta\) and
\(e^{\alpha|Z|}\le e^{\alpha Z}+e^{-\alpha Z}\).
The same estimate holds for the reverse endpoint packet with
conjugate angle. In (4), \(\alpha<n/2\) by \(b>2\).

## Subtracting the two modulus equations

Write \(K_m(x,\theta)\) for the log modulus of the forward
endpoint packet minus that of the reverse packet, evaluated
at the arguments \(e^{x+i\theta},e^{-x-i\theta}\), respectively.
The finite geometric sum gives \(m x_m=K_m(x_m,\theta)\)
at a normalized root of log radius \(x_m\). If \(\varepsilon<1\)
bounds (4) for both orientations, then
\[
 |\partial_xK_m|\le2/\sin\theta,\qquad
 |K_N-K_n|\le 2[-\log(1-\varepsilon)]
\]
at fixed \(x\). Hence
\[
 N|x_N-x_n|
  \le \frac{b}{b-2}
        \left(\frac{11}{n}+2[-\log(1-\varepsilon)]\right)
                                                               \tag{5}
\]
at two normalized roots sharing an angle, using
\(|x_n|\le11/n\) and \(N-2/\sin\theta\ge(1-2/b)N\).
[Uniform normalized radial wall](MULTIVARIATE_PRIVATE_UNIT_IMBALANCE_TWENTY_FIVE.md)

At a common physical root, \(x_N-x_n=\log(G_N/G_n)\).
The exact endpoint drift bound differs from
\(\Delta\log(N/n)\) by at most \(V/(nN)\).
Since \(\log(N/n)>1/N\), the right side of (5) plus \(V/n\)
is a strict upper bound for \(|\Delta|\) whenever \(\Delta\ne0\).
[Endpoint drift estimate](MULTIVARIATE_ALL_DEGREE_SIGNED_IMBALANCE.md#an-explicit-linked-radial-condition)

For \(n\ge56,\ k\ge2\), the angle boxes give
\(b\ge n\sin(3\pi/n)>28/3\).
The inequality follows from
\(\sin u\ge u-u^3/6\), \(157/50<\pi<22/7\).
On \(|t|\le\alpha\), (2) gives
\(\Lambda''\le(13/8)V/(nN)\), (3) gives
\(|\Lambda'|\le V/(7N)\), and (1)–(2) give
\(\Lambda\le1/16\) by Taylor's bound.
Thus \(M\le16/15\) and
\[
 M''(\alpha)+M''(-\alpha)
       \le(32/15)(13/8+6/49)V/(nN)<(15/4)V/(nN).
\]
Equation (4) gives \(\varepsilon\le9/40\), using \(V<61/12\).
Use \(-\log(1-\varepsilon)\le
\varepsilon+\varepsilon^2/[2(1-\varepsilon)]\).
The resulting integer bound is
\[
      \frac{14}{11}\left(\frac{11}{56}+\frac9{20}
           +\frac{(9/40)^2}{1-9/40}\right)
                            +\frac{61}{12\cdot56}<1 .
\]

For \(n\ge24,\ k\ge3\), use \(b\ge n\sin(5\pi/n)>14\)
by the same sine estimate. On \(|t|\le\alpha\),
\(\Lambda''\le(11/8)V/(nN)\),
\(|\Lambda'|\le V/(10N)\), and \(\Lambda\le1/25\).
Then \(M\le25/24\) and the sum of the second derivatives
of \(M\) is at most \(3V/(nN)\).
Here (4) gives \(\varepsilon\le1/12\), and the integer bound is
\[
                  \frac76\left(\frac{11}{24}+\frac2{11}\right)
                             +\frac{61}{12\cdot24}<1.
\]
Thus any common physical root in either range must have
\(\Delta=0\), where the balanced theorem excludes it.
