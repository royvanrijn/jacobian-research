# A root-force reserve for the shared mirrored gap-two pair

The exact ordinary pair \((1,3),(3,1)\) has
[circular interlacing in every degree](MULTIVARIATE_MIRRORED_CRITICAL_GAP_TWO_INTERLACING.md).
At adjacent orders \(n,N=n+1\), the smallest circular distance between
their root arguments is greater than \(10^{-10}N^{-2}\). This also
gives a neighborhood of ordinary mirrored profiles having disjoint
consecutive root arguments in every degree.

## The force sign for a real shape parameter

For \(B\ge n\) real, set
\[
 F_{n,B}(z)=\sum_{k=0}^n\binom nk
       \Gamma(B+2k+1)\Gamma(B+2n-2k+1)z^k,
 \qquad f_{n,B}(\theta)=e^{-in\theta/2}F_{n,B}(e^{i\theta}).
\]
All zeros are simple and circular by the
[factorial circle-root theorem](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md#mirrored-profiles-at-every-fixed-order),
using the two Gamma factors with common shape \(B+1\) and mirrored
private exponents two. At a zero in \(0<\theta<\pi\), let
\(r=f_{n,B}''/f_{n,B}'\) and \(t=\tan(\theta/2)\).
Define \(S=\sum_{j\ne i}\csc^2((\theta_i-\theta_j)/2)\) as in the
[exact ordinary proof](MULTIVARIATE_MIRRORED_CRITICAL_GAP_TWO_INTERLACING.md#a-mesh-bound-rules-out-a-wrong-sign).
The coefficient recurrence used there, with \(B\) in place of the
baseline \(n\), gives
\[
 S=r^2-\frac{b}{3}tr-\frac{K+1}{3},\quad
 h=B+n+1,\quad b=4h+2-2n,\quad K=h(h+1)-n(2h+1).               \tag{1}
\]
Indeed the centered differential equation has real part
\(nh(h+1)/2+b\partial^2\) and imaginary part \(4\partial^3-K\partial\).
Since \(K=(B+1)^2+(B+1)-n^2>0\), equation (1) precludes \(r=0\).
The positive-half angular roots vary continuously with \(B\).
They stay simple, \(1\) is never a root, and a fixed root at
\(\theta=\pi\) occurs precisely when \(n\) is odd, by reciprocity
and simplicity. Thus the sign of \(r\) is constant along each
positive-half root branch.

At fixed \(n\) and \(B\to\infty\), put \(x=\sqrt B\cot(\theta/2)\).
The beta integral for the two Gamma factors gives the centered
polynomial up to a positive factor as
\(\sin(\theta/2)^n q_{n,B}(x)\), with \(q_{n,B}\) proportional to
\[
 \int_{-1}^{1}(1-u^2)^B
       \big(x(1+u^2)+2i\sqrt B\,u\big)^n\,du .
\]
Multiplying this integral by \(\sqrt{B/\pi}\) makes its coefficients
converge to those of
\(q_n(x)=\pi^{-1/2}\int_{\mathbb R} e^{-w^2}(x+2iw)^n\,dw\).
In fact \(w=\sqrt B\,u\) gives
\((1-w^2/B)^B\le e^{-w^2}\) on the integration interval and
\(1+w^2/B\le2\), providing an integrable coefficientwise bound.
The generating function \(e^{xt-t^2}\) for \(q_n/n!\) makes \(q_n\)
the Hermite polynomial with variance two, with simple real zeros,
and gives \(2q_n''-xq_n'+nq_n=0\).
At a positive zero of \(q_n\), the derivative ratio is \(x/2>0\).
Because the roots and the derivative ratios converge at these
simple roots, the angular identity
\[
 r=(n-1)\cot(\theta/2)
       -\frac{\sqrt B}{2}\big(1+\cot^2(\theta/2)\big)
          \frac{q_{n,B}''(x)}{q_{n,B}'(x)}
\]
shows that \(r<0\) for large \(B\), hence for all \(B\ge n\).

This sign also gives an exact comparison within degree \(n\).
The identity \(f_{n,B+1}=(h^2+4\partial^2)f_{n,B}\) follows from the
two Gamma increments. The sign of \(f_{n,B+1}\) at each positive-half
root of \(f_{n,B}\) is therefore opposite that of \(f_{n,B}'\).
In increasing order, each positive-half root moves strictly toward
\(\pi\) under this unit change of shape, remaining below the next
positive-half root of the lower shape whenever there is one.
Indeed the signs alternate at the lower roots. The sign just below
\(\pi\) of either polynomial is \((-1)^{\lfloor n/2\rfloor}\),
since both are positive at \(\theta=0\) and each has exactly
\(\lfloor n/2\rfloor\) zeros in \((0,\pi)\). This gives the last
interval and exhausts the positive-half roots.

## An angular bound at the ordinary shape

Return to \(B=n\), and write \(f=f_{n,n}\),
\(g=f_{N,N}\). As in the
[exact ordinary proof](MULTIVARIATE_MIRRORED_CRITICAL_GAP_TWO_INTERLACING.md#the-central-root-and-the-cut),
the angular mesh estimate gives
\[
 S\le \pi^2 L/12<5N^2/3,\qquad
 L=2n^2+\tfrac32 n+\tfrac12<2N^2 .
\]
At a root in \(0<\theta<\pi\), (1) now reads \(S=r^2-2Ntr-N\).
Since \(r<0\) and \(t>0\), this gives
\(|r|\le\sqrt{S+N}<2N\). The same absolute bound holds at the
reflected roots, and \(r=0\) at any central root by parity.

The same identity yields
\((4n+3)+2tr\ge4n+2-S/N>2n\) on the positive half.
The angular mesh through the cut at \(1\) is twice the smallest
positive argument and is at least \(2\pi/\sqrt L\).
Consequently \(\sin(\theta/2)\ge1/\sqrt L>1/(2N)\) at every
positive-half root. Equation (3) of the
[exact ordinary proof](MULTIVARIATE_MIRRORED_CRITICAL_GAP_TWO_INTERLACING.md#exact-signs-in-terms-of-the-lower-roots)
then gives \(|g/f'|>12nN\ge6N^2\).
For a central root, its separate formula in that proof gives
\[
 \frac{g(\pi)}{f'(\pi)}
  =\frac{4(2h+3)}{6N+2}
       \big(h((4n+8)h-2n)+2(3n+2)\big)>6N^2 .
\]
Here \(h=2n+1\ge N\), \((4n+8)h-2n\ge4Nh\), and
\(6N+2\le7N\). Reflection covers the other roots.

To control \(g\) near any lower root \(\theta_i\), put
\(y_j=\cot((\theta_i-\theta_j)/2)\) for \(j\ne i\). The centered
sine factorization gives
\[
 \frac{f(\theta_i+u)}{f'(\theta_i)}
       =2\sin(u/2)\prod_{j\ne i}
               \big(\cos(u/2)+y_j\sin(u/2)\big),
 \qquad \sum_{j\ne i} y_j=r,\qquad
 \sum_{j\ne i} y_j^2\le S .
\]
Let \(\delta=1/(10N)\) and allow \(u\) to be complex with
\(|u|\le\delta\). Write a factor in the product as \(1+w_j\).
The bounds \(|\cos(u/2)-1|\le |u|^2/4\) and
\(|\sin(u/2)|\le|u|\) give
\[
 \max_j|w_j|<1/4,\qquad
 \left|\sum_{j\ne i}w_j\right|\le N\delta^2/4+2N\delta<1/4,
 \qquad
 \sum_{j\ne i}|w_j|^2\le N\delta^4/8+2S\delta^2<1/20 .
\]
Using \(\log|1+w|\le\Re w+|w|^2\) on \(|w|\le1/4\)
bounds the product by \(e^{3/10}<2\). Thus
\(|f(\theta_i+u)|\le4\delta |f'(\theta_i)|\).
Cauchy's estimate on radius \(\delta/2\) gives, for real
\(|\theta-\theta_i|\le\delta/2\) and \(0\le j\le5\),
\[
  |f^{(j)}(\theta)|\le 8j!(20N)^{j-1}|f'(\theta_i)| .
\]

In equation (2) of the exact ordinary proof, write
\(g=\sum_{j=0}^4 a_j f^{(j)}\). Expanding \(q=h^2 f+4f''\) gives
\[
 (a_0,a_1,a_2,a_3,a_4)=
 \big(2c(h+1)(h+2)h^2,\ 4s(2h+3)h^2,\
             8c(3h+2),\ 16s(2h+3),\ -32c\big).
\]
On the real line the absolute values are bounded in sequence by
\(48N^4,88N^3,56N,88N,32\). The corresponding derivative
bounds are half as large. The five product-rule terms, with the
Cauchy bounds above, give \(|g'(\theta)|<10^{10}N^4 |f'(\theta_i)|\)
for \(|\theta-\theta_i|\le\delta/2\).
Since \(|g(\theta_i)|>6N^2 |f'(\theta_i)|\), no upper root
lies at circular distance at most \(10^{-10}N^{-2}<\delta/2\).
This proves the stated reserve.

## Ordinary mirrored tails

Take the core \((1,3),(3,1)\) and add ordinary shared mirrored
pairs \((\ell_i+d_i,\ell_i),(\ell_i,\ell_i+d_i)\) with
\(\ell_i,d_i\ge1\). Let \(J\ge3\) be any integer at least the
[computable absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
If
\[
                 2\sum_i d_i^2/\ell_i<\frac{10^{-20}}{\pi^2J^5},
\]
the full profile has disjoint consecutive root arguments in all
degrees. For \(N<J\), the
[same-degree angular transport of a mirrored tail](MULTIVARIATE_MIRRORED_UNIT_CORE_ANGULAR_TRANSPORT.md#transfer-to-the-ordinary-tail)
moves both orders by a total of at most
\(\pi\sqrt{2J\sum_i d_i^2/\ell_i}<10^{-10}J^{-2}\).
The exact core reserve therefore implies that these pairs still
strictly interlace in the angular ordering cut at \(1\). Both
polynomials stay real reciprocal with simple circular roots along
the symmetric interpolation used in the transport, so the cut
cannot be crossed by a root. The absolute theorem covers \(N\ge J\).
