# Derivative signs at same-height zeros of the half-line cubic

For \(T\in\mathbb R\), let
\[
                 F(T,H)=\int_0^\infty e^{-u^3/3-Tu^2/2+Hu}\,du .
\]
Suppose \(H_i=x_i+i y\), \(i=1,2\), are two zeros for the
same \(T\), with \(x_1<x_2\) and \(y>0\). Then their left
member is simple and
\[
                         \Im F_H(T,H_1)>0.                 \tag{1}
\]
If \(x_1\ge0\), the right member is also simple, with
\(\Im F_H(T,H_2)<0\). In this range there are at most two
zeros at a given positive imaginary height. When two
occur, their local zero branches satisfy
\[
          \frac{d}{dT}\bigl(\Im H_2(T)-\Im H_1(T)\bigr)>0   \tag{2}
\]
at the crossing. For \(T\ge0\) both signs and (2) apply
to every horizontal pair. Conjugation gives the
lower-half-plane version.

## Resolvent representation on the interval between the zeros

Put \(a=T/2\), \(I=(x_1,x_2)\), \(f(x)=e^{a x}\), and
\(U(x)=e^{a x}F(T,x+i y)\). The cubic integral gives
\(F_{HH}+T F_H-H F=1\). Therefore, with \(L\) the real
Dirichlet operator on \(I\),
\[
          (L+i y)U=-f,\qquad
          L=-\frac{d^2}{dx^2}+x+a^2,\qquad
          U(x_1)=U(x_2)=0.
\]
At every upper zero \(x_i+a^2>0\). Indeed the real
density in the defining transform is decreasing for
\(x_i\le-a^2\) if \(T<0\), and already for \(x_i\le0\)
if \(T\ge0\). Pairing its positive and negative sine
half periods gives a strictly positive imaginary part in
those ranges. Thus the potential on \(I\) is positive
and the Dirichlet spectrum strictly positive. Hence
\[
                U(x)=-\int_0^\infty e^{-i y s}u(s,x)\,ds,
                \qquad u(s,\cdot)=e^{-sL} f .             \tag{3}
\]
First suppose \(x_1\ge0\). The positive solution \(u\) of the
Dirichlet heat problem is strictly decreasing in its time
parameter in the interior. Indeed the stationary function \(f\) is a
supersolution with \(Lf=x f\ge0\), and it is positive at
the endpoints. The maximum principle therefore gives
\(u(h,x)<f(x)\) in \(I\) for \(h>0\). For \(s,h>0\),
\[
          u(s,x)-u(s+h,x)=e^{-sL}(f-u(h,\cdot))(x)
\]
is again the heat evolution of a strictly positive interior
function. At either endpoint its inward normal derivative is
strictly positive, by the parabolic Hopf lemma. Thus both
inward fluxes
\[
               A_1(s)=u_x(s,x_1),\qquad
               A_2(s)=-u_x(s,x_2)
\]
are positive and strictly decreasing for \(s>0\).

For the remaining case \(x_1<0\) the zero localization
gives \(a<0\) and \(x_1+a^2>0\). Write \(L_0=x_2-x_1\).
For each \(s>0\), make the heat rescaling
\[
 w_s(t,r)=e^{s(x_1+a^2)t-a x_1}\,
             u(st,x_1+\sqrt{s}\,r).
\]
On \(0<r<L_0/\sqrt{s}\) this solves
\[
  \partial_t w_s=\partial_r^2 w_s-s^{3/2}r w_s,\qquad
  w_s(0,r)=e^{a\sqrt{s}r},
\]
with Dirichlet endpoints. Increasing \(s\) shrinks the
domain, increases the potential, and decreases the
initial data because \(a<0\). The maximum principle
therefore orders the solutions at \(t=1\) on the smaller
domain, hence orders their inward derivatives at \(r=0\).
Since
\[
 A_1(s)=e^{a x_1-s(x_1+a^2)}s^{-1/2}\partial_r w_s(1,0),
\]
the positive left flux is strictly decreasing in this
case as well.

The fluxes are integrable in \(s\). Near zero, comparison
with the half-line heat evolution of \(\|f\|_\infty\)
bounds each by \(\|f\|_\infty/\sqrt{\pi s}\): the two
comparison functions are the corresponding error functions
of the distance to the endpoint divided by \(2\sqrt s\).
For large \(s\), the positive solution at \(s=1\) is bounded
by a finite multiple of the first Dirichlet sine function,
since its endpoint derivatives are finite. Evolving this
upper bound with the zero-potential Dirichlet heat equation
gives an exponential flux bound, by the same comparison.
These bounds also dominate the difference quotients of the
heat solutions at the endpoints, by the error-function
and sine comparisons. This justifies differentiating
(3) there.

For a positive strictly decreasing integrable function
\(A\), pairing successive half periods gives
\[
 \int_0^\infty A(s)\sin(y s)\,ds
   =\sum_{k\ge0}\int_{2k\pi/y}^{(2k+1)\pi/y}
      \bigl(A(s)-A(s+\pi/y)\bigr)\sin(y s)\,ds>0 .
\]
Consequently (3) gives \(\Im U'(x_1)>0\) always,
and \(\Im U'(x_2)<0\) when \(x_1\ge0\).
At a zero, \(U'(x_i)=e^{a x_i}F_H(T,H_i)\), giving
the claimed field derivative signs. Three same-height
zeros with nonnegative real parts would give opposite
signs at the middle zero from the adjacent intervals,
so at most two can occur in this range.

At a simple zero, \(F_T=-F_{HH}/2\) and the field
equation give \(H_T=(1/F_H-T)/2\). Thus
\(\Im H_{1,T}<0\) always, and \(\Im H_{2,T}>0\) when
\(x_1\ge0\), proving (2).

## Two identities for individual zeros

At a simple branch let \(P(T)=1/F_H(T,H(T))\).
Differentiating the two equations for \(F\) gives the
closed system
\[
                  H_T=(P-T)/2,\qquad
                  P_T=\frac{P}{2}(H+TP-P^2).              \tag{4}
\]
Indeed \(F_{HH}=1-TF_H\) at a zero and
\(F_{HHH}=(H+T^2)F_H-T\) there. In particular
\(\partial_T\Im P=(\Re P)(\Im H)/2\) wherever \(P\) is
real. A speed zero \(2H-T H_T=0\) for nonzero \(T\)
requires \(P=T+4H/T\). Every upper zero satisfies
\(\Re H>-T^2/4\) for \(T<0\), again by the decreasing
density sine transform. Thus a speed zero for \(T<0\)
would require both the real and imaginary parts of \(P\)
to be negative. For \(T>0\) it requires both to be positive.

There is also a boundary-slope formula from the Airy
equation. Write \(Z=H+a^2\) and \(\omega=e^{2\pi i/3}\).
For any upper zero, simple or multiple,
\[
 F_H(T,H)
  =-\omega \int_0^\infty e^{a\omega s}
              \frac{\operatorname{Ai}(\omega^2 Z+s)}
                   {\operatorname{Ai}(\omega^2 Z)}\,ds .  \tag{5}
\]
The denominator is nonzero since the Airy zeros are
negative real and \(Z\) is in the upper half-plane.
To prove the formula, let
\(V(s)=e^{a(H+\omega s)}F(T,H+\omega s)\) and
\(A(s)=\operatorname{Ai}(\omega^2 Z+s)\).
Then \(V(0)=0\), and
\[
             (A V'-A'V)'=\omega^2 A(s)e^{a(H+\omega s)} .
\]
The Wronskian tends to zero at infinity by Airy decay
and the endpoint expansion of \(F\) along the ray of
argument \(2\pi/3\). Evaluation at zero gives (5).
The needed Airy zero location and sectorial decay are
the [Airy zero formulas](https://dlmf.nist.gov/9.9) and
[Airy asymptotics](https://dlmf.nist.gov/9.7). The endpoint
expansion is obtained by the same rotated-ray estimate
as in the
[half-line cubic comparison](MULTIVARIATE_HALF_LINE_CUBIC_ZERO_SEPARATION.md#an-exact-airy-decomposition).
