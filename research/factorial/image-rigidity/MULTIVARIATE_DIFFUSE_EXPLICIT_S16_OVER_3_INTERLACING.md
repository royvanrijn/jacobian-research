# An explicit Gaussian complex-disk comparison down to \(S=16/3\)

For \(B_j(z;S)=\sum_k\binom jk e^{-S k(j-k)/(2j)}z^k\),
the roots of \(B_n\) and \(B_{n+1}\) strictly interlace on the
unit circle whenever \(S\ge16/3\) and \(n\ge10^7\).
Consequently all physical upper time ratios give Gaussian
angular separation, by the
[inner comparison and outward motion](MULTIVARIATE_GAUSSIAN_REFERENCE_ANGULAR_INTERLACING.md).

The proof gives a nonvanishing polynomial half and analytic logarithm
on \(|z|\le1001/1000\) in the new range \(16/3\le S\le6\).
For \(S\ge6\) the
[earlier half comparison](MULTIVARIATE_DIFFUSE_EXPLICIT_S6_INTERLACING.md)
has a larger disk and lower degree threshold.
The new endpoint estimate uses both the modulus and phase of
the saddle, including near the negative \(z\) direction.

## Saddle and normalized contour

Put \(a=S/2\in[8/3,3]\), \(R_0=1001/1000\), and use
\(G(x,z)=e^{-x^2/(4a)}(1+ze^{x-a})\) and its half \(I_j\)
from the \(S\ge6\) proof. The saddle equation in \(q\) is
\[
   q=z e^{-a}\exp(2a q/(1+q)),\qquad x_*=\frac{2a q}{1+q}.
\]
At the Rouché radius \(q_0=13/100\) use the maximum
\(\Re(q/(1+q))=13/113\). Put \(\kappa=87/113\) and
\(Q=129/1000\). Since \(e^{-a\kappa}\le e^{-232/113}<Q\)
and \(R_0Q<q_0\), there is one simple analytic saddle on the disk.
It satisfies \(|q/z|\le Q\), with its continuous value at zero,
\(\Re x_*<7/10\), and \(|x_*|<9/10\).
Monotonicity of \(a e^{-\kappa a}\) gives
\(2a|q|/(1-|q|)^2 \le (16/3)R_0 Q/(1-q_0)^2<91/100\).
Writing \(F_*=\log G(x_*,z)\), \(A_*=-F_{xx}(x_*,z)\), and
\(C_*=(2a A_*)^{-1/2}\), we have
\[
 |F_*|<11/50,\qquad |zF_*'(z)|<3/20,\qquad
 \Re A_*\ge3/200,\qquad |C_*|>7/10,\qquad |\log C_*|<6 .
\]
For example \(|F_*|\le3(q_0/(1-q_0))^2+q_0/(1-q_0)\),
and the last two bounds follow from \(|2a A_*-1|\le91/100\).

Put \(b=a/2\). On the horizontal through \(x_*\) with
\(\Re x\le b\), the modulus of \(p=ze^{x-a}\) is at most
\(R_0e^{-4/3}<53/200\).
The maximum of \(\Re[p/(1+p)^2]\) at fixed modulus \(r\)
is \(r/(1+r)^2\); the same coefficient check as in the
\(S\ge6\) proof works up to \(r=53/200\), since
\(1-2r-6r^2-2r^3>0\).
Thus \(2a\Re[p/(1+p)^2]\le
6(53/200)/(1+53/200)^2<199/200\).
The horizontal phase has real second derivative at most
\(-1/(400a)\le-1/1200\).

For \(b\le\Re x\le a\), the normalized modulus is
nonincreasing in \(\phi=\arg z\) at fixed \(|z|\) and
on the horizontal through \(\Im x_*\).
This is the horizontal angular comparison from the
\(S\ge6\) proof: \(u=\phi+\Im x_*\) increases on \([0,\pi]\),
and \(|q|\le |p|\le1/|q|\), so the two values of
\(\Im[p/(1+p)]\) have the required ordering.
At real positive \(z\), the spatial maximum occurs at the
endpoints since \(1/x+1/(2a-x)\le1\) on \([a/2,a]\).
At \(x=b\), the real saddle concavity gives log loss at least
\((19/30)^2/2400>1/6400\).
At \(x=a\), the normalized bound is at most
\[
         \frac{e^{-a/4}(1+R_0)}{1+R_0 e^{-a}}
           \ \le\ \frac{e^{-2/3}(1+R_0)}{1+e^{-3}}
           \ <\ \frac{6499}{6500}.
\]
The same upper bound follows at \(b\) from the log loss.

The vertical connector at \(x=a+iv\) decreases in modulus
for \(0\le v\le y=\Im x_*\), exactly as in the \(S\ge6\) proof:
the required inequality is
\(|q/z|\le |1+q|^2/(1+R_0)^2\), since
\(\sin(\phi+v)\ge(v/y)\sin u\).
For clarity, here is the sharper normalized comparison of its
real endpoint at \(v=0\).
Write \(r=|q|\), \(x=\cos u\), and \(D_q=1+2r x+r^2\).
Then \(y=c\sin u\), with \(c=2ar/D_q\le91/100\).
The normalized endpoint modulus is nonincreasing in \(\phi\)
provided
\[
 \frac{D_q}{1+2|z|\cos\phi+|z|^2}
                   \frac{\sin\phi}{\sin u}\ \ge\ |q/z| .
\]
For \(x\ge0\), use
\(\sin(u-y)\ge\sin u(1-cx-c^2(1-x^2)/2)\).
This bounds the numerator \(D_q\sin\phi/\sin u\) below by
\(1-Bx-C(1-x^2)\), with \(B=431/1000\), \(C=238/1000\).
Indeed \(2(a-1)r\le B,\ 2a^2r^2\le C\) by decreasing
\((a-1)e^{-\kappa a}\) and \(a e^{-\kappa a}\) on this range.
Thus the quotient in the required inequality is at least
\((1-C-B^2/(4C))/(1+R_0)^2>Q\).
For \(x\le0\), use \(\sin\phi/\sin u\ge\cos(91/100)\)
and \(\cos\phi\le\sin(91/100)\).
The resulting lower bound
\[
   \frac{(1-q_0)^2\cos(91/100)}
          {1+R_0^2+2R_0\sin(91/100)}>Q
\]
follows from the alternating Taylor bounds through degrees
six and five, respectively. This proves the connector
comparison, including by continuity at the endpoints of
the angular interval. Therefore the normalized far
horizontal and connector moduli are at most \(6499/6500\)
over a combined length less than four.


## Explicit adjacent phase estimate

On \(x=x_*+t\) with \(|t|\le1/20\), one has
\(|p|<(13/100)e^{1/20}<7/50\).
The third and fourth derivatives of \(F=\log G\) are
bounded in modulus by \(3/10\) and \(1/2\), respectively,
by the same rational expressions as in the \(S\ge6\) proof.
Put \(F(x_*+t)-F_*=-A_*t^2/2+E(t)\),
\(c_3=F_{xxx}(x_*)/6\). Then
\[
 |E(t)|\le |t|^3/20,\quad |E(t)-c_3t^3|\le |t|^4/48,\quad
    (\Re A_*)t^2/2-|E(t)|\ge t^2/200 .
\]
Subtract the cubic term in the symmetric local Gaussian integral.
The normalized integrated residual is bounded by the
Gaussian fourth and sixth moments with coefficients
\(j/48\) and \(j^2/800\), using the damping \(e^{-jt^2/200}\).
It is less than \(87500/j\), since the full Gaussian mass
factor is \(\sqrt{200/(4a)}<9/2\) and the two moment
coefficients after division by \(1/j\) are \(625\) and \(18750\).

At \(|t|=1/20\) the real phase loss is at least \(1/80000\),
and the outward derivative has magnitude at least \(1/3000\)
with decreasing phase. For the derivative use
\((3/200)/20-(3/20)/20^2>1/3000\).
Both bounds then persist by left concavity to the end
of the left piece. The matching Gaussian also satisfies
these two loss bounds on its two tails. Their normalized
integral bounds, divided by \(|C_*|\), contribute at most
\(3100e^{-j/80000}/\sqrt j\) in total.
With the far-contour bound and the local residual above,
the analytic half is
\[
 I_j=e^{jF_*}C_*(1+\eta_j),\qquad
 |\eta_j|\le \frac{125000}{j}
       +\frac{3100}{\sqrt j}e^{-j/80000}
       +6\sqrt j(6499/6500)^j<1/50
                           \qquad(j\ge10^7).               \tag{1}
\]

The adjacent estimate uses the real-order continuation of the
normalized local integral after substituting \(y=\sqrt j\,t\),
as in the \(S\ge6\) proof. The derivative of the exponent is
\(E(t)-tE'(t)/2\), of modulus at most \(|t|^3/8\).
The scaled Gaussian majorant is \(e^{-y^2/200}\).
The differentiated local integral is therefore bounded by
\(1300/j^{3/2}\) plus the moving-endpoint bound
\(e^{-j/80000}/(50\sqrt j)\).
Integrate between \(n\) and \(n+1\), then add at both orders
the nonlocal left horizontal and bounded far-contour costs.
The Gaussian completion cancels in this adjacent comparison.
The result is
\[
 |\eta_{n+1}-\eta_n|
 \le \frac{1300}{n^{3/2}}
       +\frac{3100}{\sqrt n}e^{-n/80000}
       +12\sqrt n(6499/6500)^n<10^{-6}
                            \qquad(n\ge10^7).               \tag{2}
\]
The inequalities in (1) and (2) follow by evaluating the decreasing
upper terms at \(10^7\), using \(e>2\) and
\((6499/6500)^j\le e^{-j/6500}\).

Put \(H_j=jF_*+\log C_*+\log(1+\eta_j)\).
The logarithm of the amplitude cancels in the adjacent comparison.
By (1)–(2), the difference of the two logarithmic errors
is less than \(1/500000\) on \(|z|\le R_0\).
Cauchy with radius \(1/1000\), and \(|zF_*'|<3/20\),
give on \(|z|=1\)
\[
      |H_{n+1}'-H_n'|<1/5,\qquad
      |H_j'|\le (3/20)j+1000(6+1/49)<j/2.
\]
The strict phase comparison for the symmetric halves is then
the same as in the \(S\ge6\) proof.

