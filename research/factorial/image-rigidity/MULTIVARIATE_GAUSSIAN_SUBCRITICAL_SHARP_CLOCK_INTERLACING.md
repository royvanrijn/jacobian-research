# All-degree Gaussian angular interlacing through critical coupling

Let
\[
 Q_j(z;t)=\sum_{k=0}^j\binom jk e^{-t k(1-k/j)}z^k.
\]
**Theorem.** For every \(0<t\le2\) and every \(N=n+1\ge2\), the
unit-circle roots of \(Q_n(\cdot;t)\) and \(Q_N(\cdot;t)\) strictly
interlace. In raw heat times, this gives interlacing for every lower
time \(t/n\) and upper time in \([t/N,t/n]\).
The same \(t\) in the two polynomials is the sharp physical Gaussian
comparison. The theorem is about the Gaussian reference.

We use the [exact saddle phase](MULTIVARIATE_GAUSSIAN_EXACT_SADDLE_PHASE.md).
At every positive root angle, \(0<\theta<\pi\), the angle lies
above the zero-free edge. In its notation the right contour
\(z=x+i y(x)\), \(x>0\), starts at \(i(\pi-\theta)\), ends at
infinity, and is a strictly decreasing graph through the saddle
\(z_*=x_*+i y_*\).
The complex phase
\[
 f(z)=-z^2/(4t)+\log\cosh((z+i\theta)/2)
\]
has constant imaginary part \(\phi\in(0,\pi/2)\) on the graph.
The real part \(U(x)=\Re f(x+iy(x))\) strictly increases to
\(U_*=U(x_*)\) and then strictly decreases. This follows from
\(U'=I_y(1+y'^2)\), where \(I=\Im f\). Write
\(\beta(x)=\arctan(-y'(x))\) and \(\beta_*=\beta(x_*)\).
The exact degree-phase correction is
\[
  \chi_j=\arg(a_j-i b_j),\qquad
  a_j=\int e^{jU}\,dx,\qquad b_j=\int e^{jU}(-y')\,dx,
\]
and the phase at the \(i\)-th positive root is
\(j\phi+\chi_j=(i-\tfrac12)\pi\).
Since the previous saddle phase theorem treats \(\phi\ge\pi/4\),
throughout the proof of the new estimate we may take \(\phi\le\pi/4\).

## On the right of the saddle

Put \(C=\cosh x\), \(S=\sinh x\), and \(c=\cos(\theta+y)\).
On the right branch one has \(I_{xy}<0\). In fact
\[
 2I_{xy}=-1/t+(1+C c)/(C+c)^2 .
\]
For \(C\ge2\), the last rational function is increasing in \(c\)
on \([-1,1]\) and is at most \(1/(C+1)<1/t\).
For \(1<C<2\), it has its maximum at
\(c_*=C-2/C\), with value \(C^2/(4S^2)\).
When the minimizing line \(y=u(x)-\theta\) is in the strip, the
right branch has \(c\ge c_0=tS/x-C\).
If \(c_0\ge c_*\), compare with \(c=c_0\), where differentiating
\(I_y(x,u(x)-\theta)=0\) gives \(I_{xy}=-I_{yy}u'<0\).
If \(c_0<c_*\), then \(t<2x S/C\), so even the unconstrained
maximum satisfies \(t C^2/(4S^2)< xC/(2S)<1\).
For \(x\) beyond the interval where the minimizing line is
interior, \(t\le x\tanh(x/2)<2x S/C\), so the last bound still
applies if \(C<2\).

The same bound holds in the region between the right graph and
the horizontal line through \(y_*\), because \(u\) increases.
At the saddle \(I_{xy}<0\) gives \(0<\beta_*<\pi/4\).
Moreover \(y'\) never crosses down through \(-1\) on the right.
Indeed the twice-differentiated level equation gives
\(y''=2I_{xy}/I_y>0\) at such a point, since \(I_y<0\).
Thus the right abscissa and ordinate differences from the saddle
satisfy \(|y(x)-y_*|\le x-x_*\).

For a log-modulus deficit \(v=U_*-U(x)>0\), let \(R(v)=x-x_*\)
on this branch. Along the chord from \(z_*\) to \(z\),
\(f''\) has negative real and imaginary parts. The real sign
uses the bound on \(I_{xy}\), and the imaginary sign follows
from \(\Im f''=-I_{yy}<0\).
Taylor's integral for \(f(z)-f(z_*)\) uses the weight \(1-s\)
on this chord, while that for \(f'(z)\) uses weight one.
The former integral has modulus at most the latter, by
componentwise comparison in this quadrant. Hence
\(v/|f'(z)|\le |z-z_*|\).
Since \(|U'|=|f'(z)|\sqrt{1+y'^2}\) and \(R'(v)=1/|U'|\),
we obtain
\[
                         v R'(v)\le\sqrt2\,R(v).                  \tag{1}
\]

## Convexity on the left of the saddle

Put \(A=\pi-\theta\), \(w=z-iA=x-is\), and \(s=A-y(x)>0\).
Write
\[
 \gamma=-\arg\sinh(w/2),\qquad
 \kappa=\arg\coth(w/2),\qquad
 \zeta=\arg(-f''(z)).
\]
Here \(0<\kappa<\gamma<\pi/2\).
In fact \(\coth(w/2)\) has positive real and imaginary parts,
and its argument is \(\gamma+\arg\cosh(w/2)\), with the last
argument negative.
Since
\[
 -4 f''=\coth^2(w/2)+(2/t-1),\qquad
                  f'''=\tfrac14\operatorname{csch}^2(w/2)\coth(w/2),
\]
one has \(\zeta\le2\kappa\), and the unwrapped argument of \(f'''\)
is \(2\gamma+\kappa\).

We first show \(\zeta>2\gamma-2\pi/3\) in the left region and at
the saddle. This is immediate if \(\gamma\le\pi/3\), since
\(\zeta>0\).
Otherwise put \(q=S/x\), \(B=C-1\).
The condition \(I_y\ge0\) on this branch gives
\(C-\cos s\le tq\le2q\). Also \(q-1\le B/3\) by the power series.
Therefore \(C+1-2q\ge B/3\). The formula
\(\tan^2\gamma=\frac{C+1}{C-1}\frac{1-\cos s}{1+\cos s}\)
gives
\[
 (q-1)^2\tan^2\gamma
   \le \frac{(C+1)(2q-B)}3
   \le \frac{(q+1)^2}3 < \frac{4q^2}3 .
\]
Writing \(z_0=\tan\gamma\), the resulting bound
\((q-1)z_0<2q/\sqrt3\) implies
\(q\sin(2\gamma-2\pi/3)<\sqrt3/2\), by the double-angle identity.
Now \(-f''\) is the sum of positive real \(1/(2t)\) and
\(\tfrac14\operatorname{csch}^2(w/2)\) of argument \(2\gamma\).
The former's modulus divided by the latter's is
\((C-\cos s)/t\le q\).
Thus the imaginary part of \(-f''\) after rotation through
\(-(2\gamma-2\pi/3)\) is positive, proving the claim.

On the left, \(f'\) has argument \(\beta\), while \(z'\) has
argument \(-\beta\). Differentiating yields
\[
 \beta'=-P\sin(\zeta-2\beta),
 \qquad P=|f''z'/f'|>0 .
\]
At a stationary value of \(\beta\) on the left, \(\zeta=2\beta\).
The same equality holds at the saddle, since the tangent is
in the direction of steepest descent there.
At any such point the angle bounds just proved give
\[
 0<2\gamma+\kappa-3\beta<\pi .
\]
For the lower bound use \(\beta\le\kappa<\gamma\).
For the upper bound use \(\beta>\gamma-\pi/3\) and
\(\kappa<\gamma\).
Consequently \(\zeta'=\Im(f'''z'/f'')<0\) at those points.
An interior zero of \(\beta'\) on the left therefore has
\(\beta''>0\). At the saddle, the Taylor expansion of
\(f'(z(x))\) gives
\[
 2\beta'_*=\Im(z''_*/z'_*+f'''(z_*)z'_*/f''(z_*))
          =-\beta'_*+\zeta'_*.
\]
Hence \(\beta'_*<0\).
There can be no interior zero of \(\beta'\) on the left, since
every such zero would cross from negative to positive while
\(\beta'<0\) near the saddle. Thus \(\beta\) strictly decreases
and \(y\) is strictly convex on the left.

## The left width elasticity

For \(x<x_*\), put \(\ell=x_*-x\) and
\(b=\arctan((y-y_*)/\ell)\).
By the left convexity, \(\beta_*\le b\le\beta(x)\).
Let
\[
 E_1=\int_0^1 f''(z_*+s(z-z_*))\,ds,\qquad
 E_2=\int_0^1(1-s)f''(z_*+s(z-z_*))\,ds .
\]
The unwrapped arguments of \(-E_1,-E_2\) in \((0,\pi)\) are
\(\beta+b,2b\), respectively. This follows from
\(f'(z)=(z-z_*)E_1\), \(f(z)-f(z_*)=(z-z_*)^2 E_2=-v\),
and \(\arg(z-z_*)=\pi-b\).
Since \(\Im(-f'')>0\) on the chord,
\(\Im(-E_2)\le\Im(-E_1)\).
For the left abscissa gap \(L(v)=\ell\), this yields
\[
 \frac{v L'(v)}{L(v)}
   \le \frac{\cos\beta}{\cos b}\,
                       \frac{\sin(\beta+b)}{\sin(2b)} .          \tag{2}
\]
Here \(|U'|=|f'|\sec\beta\), as on the right.
For \(\beta\ge b\), the right side is at most \(1\) if
\(b\ge\pi/6\) and at most \(1/[4\sin b(1-\sin b)]\) if
\(b\le\pi/6\).
Indeed \(2\cos\beta\sin(\beta+b)=\sin(2\beta+b)+\sin b\).
The maximum is at \(\beta=b\) in the first case and at
\(2\beta+b=\pi/2\) in the second.

It remains to bound the saddle tangent below when \(\phi\le\pi/4\).
At the saddle, \(I_{xy}=-I_{yy}u'(x_*)\), so
\(2\beta_*=\arctan(1/u'(x_*))\).
In the next formulas \(x=x_*\), \(u=u(x_*)\),
\(q=\sinh x/x\), and
\[
 z_1=\tanh(x/2)\tan(u/2),\qquad \gamma_1=\arctan z_1 .
\]
Direct substitution of \(t=x(\cosh x+\cos u)/\sinh x\) gives
\[
 u'=\frac{(q+1)z_1+(q-1)/z_1}{2},\qquad
 \phi=\gamma_1-\frac{\sin(2\gamma_1)}{2q}.
\]
Also \(x_*\le x_c<5/2\) because \(e^{5/2}>9\), and
\(\sinh(5/2)<25/4\) gives \(q<5/2\).
Both exponential inequalities follow by summing the series for
\(e^{5/2}\) through order seven and bounding its positive tail
above by the first omitted term times \(18/13\).
Since \(1-\cos u\ge x^2/6\) for \(t\le2\), we have
\(u\ge x/\sqrt3\). Thus
\[
 (q-1)/z_1
   \le \frac{(\cosh x-1)/3}{\tanh(x/2)\,x/(2\sqrt3)}
   =\frac{2\sqrt3}{3}q<\frac76q.
\]
The phase assumption gives
\(\gamma_1\le\pi/4+1/(2q)\).
Use \(\tan \xi\le \xi/(1-\xi^2/2)\) for \(0\le\xi\le1/2\)
and the tangent addition formula. The following upper bounds for
\(z_1\) and \(u'\) use the two endpoints of each \(q\)-interval.
\[
 \begin{array}{c|c|c}
 q\text{ interval}&z_1&u'\\
 [1,5/4]&11/3&233/48\\
 [5/4,3/2]&33/13&421/104\\
 [3/2,2]&23/11&142/33\\
 [2,5/2]&39/23&2443/552
 \end{array}
\]
All entries in the last column are below \(49/10\).
Thus \(2\beta_*>\arctan(10/49)>2\arcsin(1/10)\), the last inequality
following from \(\tan(2\arcsin(1/10))=\sqrt{99}/49\).
Equation (2) consequently gives
\[
                              v L'(v)\le\frac{25}{9}L(v).         \tag{3}
\]

## The adjacent phase increment

Let \(X(v)=L(v)+R(v)\) be the total horizontal width at deficit
\(v\). Equations (1) and (3) give \(vX'(v)\le(25/9)X(v)\).
After removal of \(e^{mU_*}\), the real part of the right
contour amplitude at every real degree \(m>0\) is
\(\widetilde a_m=\int_0^\infty e^{-mv}X'(v)\,dv\).
Integration by parts gives
\(\widetilde a_m=m\int_0^\infty e^{-mv}X(v)\,dv\).
The horizontal width is \(O(\sqrt v+1)\) at infinity by Gaussian
decay along the contour, and tends to zero at zero.
The mean deficit in this integral with density proportional to
\(e^{-mv}X'(v)\) is therefore at most \(25/(9m)\).
The imaginary amplitude in absolute value has the analogous
positive Laplace representation in the total vertical width at
deficit \(v\). Its mean deficit is nonnegative.
Differentiating the argument
\(\chi_m=\arg(\widetilde a_m-i\widetilde b_m)\) hence gives
\(\chi_m'\ge-25/(18m)\), using
\(\widetilde a_m\widetilde b_m/
(\widetilde a_m^2+\widetilde b_m^2)\le1/2\).

At every degree-\(n\) root, \(\phi>\pi/(2n)\) by the exact
phase label. For \(\phi\le\pi/4\),
\[
  \psi_{n+1}-\psi_n
       \ge\phi-\frac{25}{18}\log(1+1/n)>0 ,
\]
because \(25/9<\pi\).
The increment is also less than \(\pi\), since
\(\phi<\pi/2\) and both corrections are in \((-\pi/2,0)\).
This gives the strict adjacent interlacing at those roots.
The larger phases were already covered by the exact saddle phase
theorem. The extension to larger physical upper times follows from
outward motion and the
[inner-gap comparison](MULTIVARIATE_GAUSSIAN_REFERENCE_ANGULAR_INTERLACING.md).
