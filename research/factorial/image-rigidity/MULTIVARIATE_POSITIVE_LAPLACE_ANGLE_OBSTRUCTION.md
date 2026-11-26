# Inward angular motion under positive Laplace curvature

Positive symmetric Laplace curvature can move an outer Gaussian
reference root inward and subsequently send roots off the unit
circle. Both effects occur at degree four. Thus the positivity
of the matched curvature increment alone gives no outward angular
comparison for interlacing. The
[finite small-movement theorem](MULTIVARIATE_MATCHED_LAPLACE_NONLINEAR_TRANSFER.md)
uses the positive increment with its quantitative mesh and size
conditions.

## An exact Gaussian deformation

Set \(q=99/100\), and choose \(h>0\) so that
\(\cosh(h/4)=2\). Take the normalized log-coefficient direction
\[
                     D(s)=\cosh(h(s-1/2))-\cosh(h/2).
\]
Its curvature is the Laplace transform of two positive atoms at
\(\pm h\). The endpoint values of \(D\) are zero, and
\(D(1/4)=D(3/4)=-5,\ D(1/2)=-6\).
The reciprocal quartic obtained by exponentiating the direction
is
\[
       P_t(z)=1+4q^3e^{-5t}z+6q^4 e^{-6t}z^2
                          +4q^3 e^{-5t}z^3+z^4 .
\]
At \(t=0\) this is the quadratic Gaussian reference with
parameter \(c=-4\log q>0\) in the form
\(\binom4k e^{-c k(4-k)/4}\). Its roots are simple and on the
circle. In the coordinate \(x=z+z^{-1}\), the quotient by
\(z^2\) equals
\[
                   H_t(x)=x^2+4q^3e^{-5t} x+6q^4e^{-6t}-2 .
\]
Its larger root at \(t=0\) is
\[
                   x_+=-2q^3+(1-q^2)\sqrt{4q^2+2}<-9q/5 .
\]
For the strict inequality, the square-root summand is less
than \(1/20\), whereas \(2q^3-9q/5\) is greater than \(1/20\).
It follows that
\(\partial_tH_t(x_+)|_{t=0}=-20q^3x_+-36q^4>0\), and
\(H_0'(x_+)>0\). Thus \(x_+'(0)<0\). Since the positive
displacement \(\theta\) from the negative ray has
\(x=-2\cos\theta\), this is inward angular motion. Both roots
of \(H_t\) stay simple and in \((-2,2)\) for small positive \(t\).

For an exact finite deformation, take \(e^{-t}=p=999/1000\).
The discriminant of \(H_t\), divided by eight, is
\[
                       2q^6p^{10}-3q^4p^6+1.
\]
Its numerator over \(10^{42}\) is
\[
                       -325267580761146894693907050513937721198<0.
\]
The two \(x\) roots are then nonreal, so the quartic has no
unit-circle root at this parameter.

## Directions from an actual matched symbolic increment

A nonzero positive symmetric Laplace curvature measure that gives
a finite endpoint-normalized direction \(D\) has
\[
                       1<\frac{|D(1/2)|}{|D(1/4)|}<\frac43
\]
provided the measure assigns positive mass to nonzero Laplace
coordinates.
For a pair of atoms at \(\pm h\ne0\), writing
\(r=\cosh(h/4)>1\) makes the ratio \(2(r+1)/(2r+1)\).
The ratio for a measure is its weighted average, with weights
given by the positive quantity \(|D(1/4)|\) for each summand.
A zero-coordinate curvature atom contributes the limiting
ratio \(4/3\).

In particular the exact matched increment of any nonconstant
fully shared mirrored profile has the strict inequality, by
the [rowwise positivity theorem](MULTIVARIATE_MATCHED_LAPLACE_NONLINEAR_TRANSFER.md#positive-curvature-increment-for-ordinary-shared-rows).
Apply its direction to a degree-four Gaussian reference with
\(q<1\) sufficiently close to one. As \(q\uparrow1\) the larger
\(x\) root tends to \(-2\), while the inward derivative condition
is \(x_+<-(3/2)q\,|D(1/2)|/|D(1/4)|\). Hence that matched
direction can also move a Gaussian reference root inward.
Comparing the actual consecutive factorial polynomials requires
the additional relation between the root configuration and
the symbol that produces the direction.
