# Strict Gaussian interlacing for all consecutive degrees and physical times

Let
\[
 Q_j(z;t)=\sum_{k=0}^j\binom jk e^{-t k(1-k/j)}z^k.
\]
For every \(t>0\) and \(n\ge1\), \(Q_n(\cdot;t)\) and \(Q_{n+1}(\cdot;t)\)
strictly interlace in symmetric angular order. The conclusion holds
for lower raw heat time \(t/n\) and upper raw time anywhere in
\([t/(n+1),t/n]\).

At every positive degree-\(n\) root, the exact consecutive phase
increment satisfies
\[
 \psi_{n+1}-\psi_n>\frac{\pi/2-\log4}{n}>0 .
 \tag{1}
\]
The main step is a concavity bound on the log modulus of the left
contour in the coordinate \(x^2\). It controls a ratio of Laplace
amplitudes at orders \(m\) and \(m/2\). For \(2<t<4\), the proof
uses an exact sign check over three continuous contour parameters.
The other ranges use analytic inequalities.

## Saddle phase and left curvature

Take the decreasing graph \(z(x)=x+iy(x)\) from the
[exact saddle phase construction](MULTIVARIATE_GAUSSIAN_EXACT_SADDLE_PHASE.md),
for \(\theta\) between the zero-free edge and \(\pi\).
All positive root angles are in this interval. With \(I=\Im f\)
and \(U=\Re f\), put
\[
 f(z)=-z^2/(4t)+\log\cosh((z+i\theta)/2),\qquad
 I(x,y(x))=\phi,\qquad z_*=x_*+iy_* .
\]
The graph runs from \(i(\pi-\theta)\) to infinity and passes through
the nondegenerate saddle \(z_*\). On the left of \(z_*\) one has
\(I_y>0\) and \(U'=I_y(1+y'^2)\). For \(t>2\) let \(x_0>0\) solve
\(x_0\coth(x_0/2)=t\). The saddle obeys \(x_*>x_0\), and the phase
satisfies \(\phi_t(0)=0\) and \(\phi_t'(\theta)=x_*/(2t)\). In
particular \(\phi>x_0\theta/(2t)\).

Write \(w=z-i(\pi-\theta)\), \(\beta=\arctan(-y')\), and
\[
 \gamma=-\arg\sinh(w/2),\quad
 \eta=-\arg\cosh(w/2),\quad
 \kappa=\gamma-\eta=\arg\coth(w/2),\quad
 \zeta=\arg(-f'').
\]
Here \(0<\kappa<\gamma<\pi/2\) and \(0<\eta<\pi/2\).
The identities
\[
 -4f''=\coth^2(w/2)+(2/t-1),\qquad
 f'''=\tfrac14\operatorname{csch}^2(w/2)\coth(w/2)
\]
give \(2\kappa<\zeta<\pi\) when \(t>2\), and the unwrapped
argument of \(f'''\) is \(3\gamma-\eta\).

For \(t\le2\) the [left convexity argument](MULTIVARIATE_GAUSSIAN_SUBCRITICAL_SHARP_CLOCK_INTERLACING.md#convexity-on-the-left-of-the-saddle)
gives \(\beta'<0\) on the left. That argument uses \(I_y\ge0\)
and \(t\le2\) throughout the strip and has no restriction on the
saddle phase. For \(t>2\), we proceed as follows.

On the left, \(f'(z)z'(x)=U'(x)>0\), so
\[
 \beta'=-|f''z'/f'|\sin(\zeta-2\beta).
\]
At a zero of \(\beta'\), or at the saddle, \(\beta=\zeta/2\).
At any such point \(3\gamma-\eta-3\beta<2\eta<\pi\). Its positivity
is equivalent to
\[
  2(a+b)\cos(2\eta/3)>d,\qquad
  a=\sinh^2(x/2),\quad b=\sin^2((\theta+y)/2),\quad d=t/2-1 .
 \tag{2}
\]
Indeed, putting \(H=\cosh^2(w/2)\), the identity
\(-f''=\operatorname{csch}^2(w/2)(H+d)/(2t)\) gives
\(\zeta=2\gamma-\alpha\) with \(\alpha=-\arg(H+d)\in(0,2\eta)\).
Thus the angle in question is \(3\alpha/2-\eta\). After rotation
through \(2\eta/3\), \(H+d\) has negative imaginary part precisely
when \(2|H|\cos(2\eta/3)>d\), and \(|H|=a+b\).

Failure of (2) implies \(a+b\le d\). The power series for
\((1+s^2/3)\sinh s-s\cosh s\) has nonnegative coefficients,
so \(x_0\coth(x_0/2)\le 2+x_0^2/6\). Hence
\(x^2\le4a\le4d<12d\le x_0^2\). In particular (2) holds at the saddle.

## Excluding a left curvature zero with nonpositive rotation

Suppose (2) fails at a left zero of \(\beta'\). Write
\[
 u=\theta+y,\quad S=\sinh x,\quad D=\cosh x+\cos u,\quad
 F=1+a-b,\quad T=\tan(u/2),\quad p=-y'.
\]
The level equation gives
\[
 p=\frac{t\sin u-yD}{tS-xD}.
 \tag{3}
\]
Set \(A=a/d\), \(B=b/d\), \(s=A+B\), and \(e=1-s\).
Here \(0<s<1\) and \(\sin^2\eta=A(1-dB)/s\).
If \(s\ge1/2\), failure of (2) is equivalent to
\(\cos(2\eta/3)\le1/(2s)\). The cubic \(4q^3-3q\) increases
on \([1/2,1]\). Applying the triple-angle identity gives
\[
 s^2\big(3+2s-4A(1-dB)\big)\le1 .
\]
For \(s\le1/2\) the same inequality follows from
\(s^2(3+2s)\le1\). Since \(1-s^2(3-2s)=(1-s)^2(1+2s)\),
both cases give
\[
 4B(1+dA)s^2\le(1-s)^2(1+2s).
 \tag{4}
\]
In particular \(B\le1/2\). This is immediate for \(s\le1/2\),
and otherwise
\((1-s)^2(1+2s)/(4s^2)\le1/2\), since its numerator decreases
and its denominator increases for \(1/2\le s\le1\).
Write \(B=e^2 V\). Then \(V\le2\) by \(B\le s\) for
\(s\le1/2\) and by (4) for \(s\ge1/2\).

At the curvature zero, differentiating the level equation twice gives
\[
 p^2-2\frac{I_{xy}}{I_{yy}}p-1=0 .
\]
Put \(p=(\sin u/S)M/e\), \(c=1+2eV\), and
\(J=Fc-2d(1+d)AeV\). Using \(D=2F\) in this equation gives
\[
 V(1-dB)M^2-\frac{J}{1+d}M-A(1+dA)=0.
 \tag{5}
\]
Indeed \(I_{xy}/I_{yy}=(1+\cosh x\cos u-D^2/t)/(S\sin u)\),
and the numerator equals \(2deJ/(1+d)\). Equation (5)
has a single positive root because \(A,V>0\) and \(dB=b<1\).

There are two upper bounds for \(M\). The first follows from
\(y\ge0\) in (3), using
\(1+d-(x/S)F\ge1+d-F=d(1-A+B)=dec\).
The second uses the phase of \(\cosh((z+i\theta)/2)\),
\[
 \eta_0=\arctan\big(\tanh(x/2)T\big)\le xT/2,\qquad
 \phi=\eta_0-xy/(2t)>x_0(u-y)/(2t).
\]
These give \(y\ge tT-(tT-u)/(1-x/x_0)\). In addition,
\[
 \frac{1}{1-x/x_0}\le\frac{1}{1-\sqrt{A/3}}
      \le\frac{24}{19-10A},\qquad
 \frac{u}{2T}=\int_0^1\frac{d\lambda}{1+T^2\lambda^2}
      \ge\frac{1}{1+T^2/3}=\frac{1-b}{1-2b/3}.
\]
The rational bound follows from
\((A+1/2)^2\ge2A\) and \(5/12>1/\sqrt6\). The integral bound
is Jensen's inequality applied to \(\lambda^2\).
Using these in (3) gives
\[
 M\le M_{\mathrm{height}}=\frac{1+d}{dc},\qquad
 M\le M_{\mathrm{phase}}=\frac{g}{h(1-dB)c},\quad
 h=(19-10A)(3-2dB),\quad
 g=24F(3-2dB+B)-h(1+d)A .
 \tag{6}
\]
In the second substitution \(tT-u\ge0\) since \(T\ge u/2\) and \(t>2\).
Both denominators in (6) are positive.

We show that at least one of the following polynomials is strictly
positive on the possible parameters:
\[
 \begin{split}
 d c\big[J+A(1+dA)dc\big]-V(1-dB)(1+d)^2,\qquad\\
 Jg h c +(1+d)\big[A(1+dA)(hc)^2(1-dB)-Vg^2\big].
 \end{split}
 \tag{7}
\]
For \(d\ge1\) the first polynomial suffices. To see this,
put \(r=\sqrt b\) and \(v=a+b\) in (4), obtaining
\[
        4b(1+a)v^2\le(d-v)^2(d+2v).
\]
If \(d\le a+r\), then \(a\ge1-r\). The right side increases
for \(d\ge v\), so substitution of \(a+r\) and division by \(b\)
give
\[
 4(1+a)(a+r^2)^2\le(1-r)^2(3a+r+2r^2)
                 \le a^2(3a+r+2r^2).
\]
But \(4(1+a)>3a+r+2r^2\) for \(0<r<1\). Thus \(d>a+r\).
Write \(q=(d-a+b)/(1+d)\). The first polynomial in (7),
multiplied by \(de^2/(1+d)^2\), is
\[
   ((1+a)q-b)^2-b(1-q)^2 .
\]
The smaller factor in this difference of squares is
\((1+a)q-b-r(1-q)=F(d-a-r)/(1+d)>0\), since \(q<1\).
This proves the sign.

For \(0<d<1\) the exact polynomial certificate uses the box
\(0\le d\le1,\ 0\le s\le1,\ 0\le V\le2\), with
\(e=1-s,\ B=e^2V,\ A=s-B\). Restrict to \(A\ge0,\ B\le1/2\)
and \(4V(1+dA)s^2\le1+2s\). The last inequality is (4) divided
by \(e^2\) and retains its limiting domain at \(e=0\). On this
closed domain (7) holds.
The [integer interval checker](verify_multivariate_gaussian_supercritical_phase_reserve.py)
bisects the box using outward dyadic arithmetic. It excludes
1,691 leaves by the domain tests, accepts 6,084 leaves by the
first polynomial and 5,162 by the second, and reaches depth 24.
All comparisons use integer endpoints. Every bisection retains
its midpoint in both children, and each rejection test is strict.

For the actual parameters \(d,A,V>0\), the two expressions in (7)
are the negative of the left side of (5) at the bounds in (6),
multiplied by positive denominators. If a bound is negative then
(6) already contradicts \(M>0\). Otherwise the positive certificate
puts that bound below the positive root of (5), also a contradiction.
This proves (2) at all left zeros of \(\beta'\).

We conclude \(\beta'<0\) on the left. At each hypothetical zero,
(2) gives \(\zeta'=\Im(f'''z'/f'')<0\) and therefore
\(\beta''>0\). At the saddle, expanding \(f(z(x))\) to third
order along the constant phase gives \(3\beta'_*=\zeta'_*<0\).
Every interior zero of \(\beta'\) would cross upward, so there can
be no such zero with \(\beta'<0\) near the saddle.

## Dyadic width and the phase reserve

Write \(v=U_*-U(x)>0\). On the left let
\(L(v)=x_*-x\), and on the right let \(R(v)=x-x_*\).
For \(r=x^2\) on the left, for any \(t>0\),
\[
 U_r=\frac{1+y'^2}{4}
 \left(\frac{\sinh x}{x(\cosh x+\cos(\theta+y))}-\frac1t\right)>0.
\]
The bracket decreases strictly. Indeed \(D'\ge\sinh x\) along
the decreasing graph and
\[
 \frac{d}{dx}\log(\sinh x/x)
 =\coth x-1/x<\frac{\sinh x}{\cosh x+1}\le D'/D .
\]
As \(\beta'<0\) also makes \(1+y'^2\) decrease, \(U_r\) decreases.
It follows that \(L_2(v)=x_*^2-x^2\) is concave in \(v\), since
\(L_2'(v)=1/U_r\) decreases as \(v\) increases. Thus
\(L_2(2v)\le2L_2(v)\), and \(L=L_2/(x_*+x)\) implies
\[
 L(2v)\le4L(v).
\]

On the right the bound \(vR'(v)\le\sqrt2R(v)\) from the
[right-branch argument](MULTIVARIATE_GAUSSIAN_SUBCRITICAL_SHARP_CLOCK_INTERLACING.md#on-the-right-of-the-saddle)
continues to hold. For \(t\le2\) this is the original right-branch argument. For
\(t>2\), in its case \(\cosh x\ge2\), \(x\ge x_*>x_0\) gives
\(\cosh x+1>\cosh x_0+1=t\sinh x_0/x_0>t\).
For \(1<\cosh x<2\) the minimum-height curve and both rational
bounds there are unchanged. These inequalities hold in the region
between the right branch and the horizontal through the saddle,
as required by the chord estimate. Consequently
\(R(2v)\le2^{\sqrt2}R(v)\), and the total width \(X=L+R\) obeys
\(X(2v)\le4X(v)\).

For real \(m>0\), remove the common factor \(e^{mU_*}\) from the
complex contour amplitude and write the result \(a_m-i b_m\).
Both \(a_m,b_m\) are positive. The former is
\[
 a_m=\int_0^\infty e^{-mv}\,dX(v)
     =m\int_0^\infty e^{-mv}X(v)\,dv.
\]
The second equality uses decay at infinity and \(X(0)=0\).
The width bound gives \(a_{m/2}\le4 a_m\).
Log-convexity of a positive Laplace transform then yields
\[
 E_{a_m}[v]=-(\log a_m)'
      \le\frac{2}{m}\log(a_{m/2}/a_m)
      \le\frac{2\log4}{m}.
\]
The vertical width gives the analogous positive transform for \(b_m\),
so \(E_{b_m}[v]\ge0\). Thus for \(\chi_m=\arg(a_m-i b_m)\) we have
\[
 \chi'_m=\frac{a_m b_m}{a_m^2+b_m^2}
        \big(E_{b_m}[v]-E_{a_m}[v]\big)\ge-\log4/m .
\]
At a positive degree-\(n\) root the exact label gives
\(\phi>\pi/(2n)\). Integrating in \(m\) from \(n\) to \(n+1\) proves
(1). The increment is less than \(\pi\), since \(\phi<\pi/2\)
and both corrections lie in \((-\pi/2,0)\). The root labels now
give strict adjacent interlacing, with the endpoint \(-1\) handled
by parity. Outward root motion and the
[equal-raw-time inner-gap comparison](MULTIVARIATE_GAUSSIAN_REFERENCE_ANGULAR_INTERLACING.md)
extend the result to the whole physical interval.

The Gaussian polynomials satisfy \(Q_j(z;t)=(-1)^jH_j(-z;2t/j)\)
in the unitary Hermite convention. Their location and the limiting saddle geometry
are treated by [Kabluchko](https://arxiv.org/html/2203.05533v3).
The monotone heat motion is given by
[Hou, Jiang and Newman](https://arxiv.org/html/2208.00917).
The interlacing comparison here concerns finite consecutive
degrees at their different physical clocks. Applying it to the
ordinary factorial coefficients requires a separate angular
perturbation bound.
