# Transport of the full even and odd tanh kernels

Let \(n\ge2\), \(N=n+1\), \(R=n^2/(n^2-1)\) and
\(\chi=1-1/n-1/N\). Match the Gauss factors of an ordinary
two-monomial exponent profile. Use the lower and upper tanh
coordinates \(T_{i,h},T'_{i,h}\) and outer losses \(U,W\) from the
[excess-loss theorem](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md).
The signs \(\epsilon_i\) are constant within each row. Equal rows
can be omitted.

On a row with smaller exponent \(\ell_i\ge0\) and gap \(d_i>0\), put
\(L_i=\ell_i+d_i\) and \(A=(n-1)/2\). The path joining the lower
coordinate to the contracted upper coordinate is
\[
 F_{i,h}(\lambda)
 =\frac{A}{n\ell_i/d_i+A+1/2+\lambda(h/d_i-1/2)},
 \qquad \frac nN\le\lambda\le1 .
\]
It starts at \(F_{i,h}(n/N)=T'_{i,h}/R\) and ends at
\(F_{i,h}(1)=T_{i,h}\). Choose \(b_i>0\) bounding all these paths
and write \(B=\max_i b_i\), \(E=\max_i L_i b_i/d_i\).
One may choose \(b_i=\max_{h,\lambda}F_{i,h}(\lambda)<1\), including
private rows. It suffices to check \(h=1\) at the two path endpoints, since the
denominator increases with \(h\) and is affine in \(\lambda\).
These exact path maxima also satisfy \(E\le1\), using \(A\le n\).
Define
\[
    C_+=\frac{4N}{n-1}E^2,\qquad
    C_-=\frac{2(2n+1)}{\chi(n-1)^2}E^3 .
\]

**Theorem.** For any complex-valued \(\phi\in C^2([0,B])\),
\[
 \left|\sum_{i,h}\epsilon_i
             \bigl(\phi(T'_{i,h}/R)-\phi(T_{i,h})\bigr)\right|
 \ \le\ (W-RU)
       \max\left\{C_+\|\phi'\|_\infty,\,
          C_-\|2\phi'+t\phi''\|_\infty\right\}.                 \tag{1}
\]
The unsigned statement follows by using \(\epsilon_i=1\).
The estimate applies to a full logarithmic kernel and retains
cancellation between its powers.

Indeed, fix one row and put \(T_h=F_{i,h}(\lambda)\),
\(v_h=(2h-d_i-1)/(2d_i)\), and \(m=1/(2d_i)\). Its difference is
\[
 \frac1A\int_{n/N}^1\sum_h (m+v_h)T_h^2\phi'(T_h)\,d\lambda.
\]
The constant part is bounded in modulus by
\(b_i^2\|\phi'\|_\infty/(n^2-1)\).
As a function of the affine shift in the denominator, the
derivative of \(T^2\phi'(T)\) in modulus is at most
\(b_i^3\|2\phi'+t\phi''\|_\infty/A\).
Pairing \(v_h\) and \(-v_h\) and using
\(\sum_h v_h^2=(d_i^2-1)/(12d_i)\) bounds the remaining part by
\[
   \frac{(2n+1)(d_i^2-1)b_i^3}{6d_i(n^2-1)^2}
               \|2\phi'+t\phi''\|_\infty .
\]
These two bounds equal the row versions of the entries of
the maximum in (1) times, respectively,
\[
    G_{+,i}=\frac{d_i^2}{4N^2L_i^2},
    \qquad
    G_{-,i}=\frac{\chi d_i^2(d_i^2-1)}{12N^2L_i^3}.
\]
The [integer-grid outer-loss gain](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md#two-contributions-to-the-row-gain)
gives \(\sum_i (G_{+,i}+G_{-,i})\le W-RU\), completing the proof.
For a noninteger power \(p\ge1\) the argument can be run on the
compact positive range of the finitely many paths and gives back
the power bound of the excess-loss theorem.

## Two full kernels for the factorial coefficients

For \(-1\le z\le1\) set
\[
 \begin{aligned}
 E_z(t)&=\tfrac12\log(1-z^2t^2)+\tfrac12 z^2t^2,\\
 O_z(t)&=\operatorname{artanh}(zt)-z\operatorname{artanh}t .
 \end{aligned}
\]
The even kernel has the leading quadratic term removed. Its
higher even powers give the signed skew beyond \(M_2\).
The odd kernel has the linear contribution removed and gives the
unsigned reflected mean defects. Both are applied to the
centered index grid of the [exact Cayley expansion](MULTIVARIATE_CAYLEY_OPERATOR_RAY_TEST.md#ordinary-factorial-moments).

If \(RB<1\), apply (1) to \(t\mapsto E_z(Rt),O_z(Rt)\).
Writing \(b=RB\), the following are upper bounds for the two norms
in (1), with a common factor \(R\):
\[
 \begin{array}{c|cc}
 & \|\phi'\|_\infty/R & \|2\phi'+t\phi''\|_\infty/R\\ \hline
 E_z(Rt) &
 \dfrac{|z|^4 b^3}{1-z^2 b^2} &
 \dfrac{|z|^4 b^3(5-3z^2 b^2)}{(1-z^2 b^2)^2}\\[2mm]
 O_z(Rt) &
 \dfrac{|z|(1-z^2)b^2}{(1-b^2)^2} &
 \dfrac{4|z|(1-z^2)b^2}{(1-b^2)^3}
 \end{array}
\]
Thus the signed even and unsigned odd kernels are transported from
\(RT_{i,h}\) to \(T'_{i,h}\) with the indicated loss and vanishing
orders. The condition \(RB<1\) is needed only for this dilation of
the kernels. Formula (1) itself also works when \(RB\ge1\).

For the even row use
\(E_z'(t)=-z^4t^3/(1-z^2t^2)\) and
\[
  2E_z'(t)+tE_z''(t)
    =-z^4t^3\frac{5-3z^2t^2}{(1-z^2t^2)^2}.
\]
The last quotient increases with \(t\), since it equals
\(3/(1-z^2t^2)+2/(1-z^2t^2)^2\). For the odd row use
\[
  O_z'(t)=\frac{z(z^2-1)t^2}{(1-z^2t^2)(1-t^2)} .
\]
Each denominator is at least \(1-b^2\). The derivative of the
reciprocal product, multiplied by \(t\), is at most
\(4t^2/(1-b^2)^3\). In \(2O_z'+tO_z''\) this gives the second
odd bound by combining
\(4/(1-b^2)^2+4b^2/(1-b^2)^3=4/(1-b^2)^3\).
The extra factor \(R\) in both rows comes from dilation.

This compares the logarithmic kernels at a common grid location.
The two adjacent orders have different centered grids. A
gap-scale angular comparison still needs to track their
alternating coefficient sums after this transport.
