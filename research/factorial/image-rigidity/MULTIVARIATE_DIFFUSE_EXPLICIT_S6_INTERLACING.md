# A half-saddle phase bound for Gaussian coupling \(S\ge6\)

For
\[
 B_j(z;S)=\sum_{k=0}^j\binom jk e^{-S k(j-k)/(2j)}z^k
\]
the two orders \(j,j+1\) have simple and strictly interlacing
unit-circle roots whenever \(S\ge6\) and \(j\ge10000\).
It suffices to prove \(6\le S\le8\), by the
[previous half-saddle proof](MULTIVARIATE_DIFFUSE_EXPLICIT_S8_INTERLACING.md).
The same Gaussian reference order also holds throughout the physical
upper time range by outward motion and the
[inner comparison](MULTIVARIATE_GAUSSIAN_REFERENCE_ANGULAR_INTERLACING.md).

In the new range, we prove the analytic half estimate on
\(|z|\le21/20\). The far-contour bounds compare its modulus to the
saddle value at each complex \(z\). This avoids bounding the
unnormalized contour and its saddle separately.

## Analytic saddle and left concavity

Put \(a=S/2\in[3,4]\), \(R_0=21/20\),
\(G(x,z)=e^{-x^2/(4a)}(1+ze^{x-a})\), and
\[
 I_j(z)=\sqrt{\frac{j}{4\pi a}}\int_{-\infty}^a G(x,z)^j\,dx .
\]
Then \(B_j=I_j+z^jI_j(1/z)\), and \(I_j\) has positive coefficients.
As in the \(S\ge8\) proof, use the saddle variable
\[
 q=ze^{-a}\exp(2a q/(1+q)),\qquad x_*=\frac{2aq}{1+q}.
\]
On \(|q|=9/100\) the maximum of \(\Re(q/(1+q))\) is \(9/109\).
Indeed its formula at fixed modulus \(r\) is
\((r\cos u+r^2)/(1+2r\cos u+r^2)\), increasing in \(\cos u\)
for \(r<1\). Put \(\kappa=91/109\).
The right-hand side of the saddle equation is bounded on that
circle by \(R_0 e^{-\kappa a}<R_0/12<9/100\), using
\(e^{273/109}>12\). Rouché gives one simple analytic solution.
Its symmetry under conjugation follows from uniqueness.
In particular
\[
 |q|<9/100,\quad |q/z|\le e^{-\kappa a}<1/12,\quad
 |x_*|<4/5,\quad \Re x_*\le2/3 .
\]
The bound on \(q/z\) has its continuous meaning at \(z=0\).
Monotonicity of \(a e^{-\kappa a}\) gives
\(2a|q|\le 2\cdot3 R_0/12=21/40\).

Write \(F=\log G=-x^2/(4a)+\log(1+ze^{x-a})\) where
the series converges, \(F_*=F(x_*,z)\), and
\(A_*=-F_{xx}(x_*,z)\). Then
\[
 \|F_*\|<3/20,\qquad |zF_*'(z)|=|q/(1+q)|<1/10,\qquad
 |2aA_*-1|\le16/25,\quad \Re A_*\ge9/200.                  \tag{1}
\]
For the first bound use \(a(9/91)^2+9/91<3/20\).
The last bounds use \(2a|q|/(1-|q|)^2<16/25\).
Thus \(C_*=(2aA_*)^{-1/2}\) has an analytic logarithm
with \(|C_*|>3/4\) and \(|\log C_*|<1\).

On the horizontal through the saddle and to the left of
\(b=a/2-1/5\), put \(p=ze^{x-a}\).
We have \(|p|\le p_0=R_0e^{-a/2-1/5}<1/5\).
For \(|p|\le1/5\), the maximum of
\(\Re[p/(1+p)^2]\) at a given modulus \(r\) is
\(r/(1+r)^2\), as in the \(S\ge8\) proof.
The same derivative test applies on \([0,1/5]\), since
\(1-2r(1+r^2)-8r^2>0\).
Moreover \(2a p_0/(1+p_0)^2\) decreases with \(a\ge3\):
its logarithmic derivative is at most \(1/3-1/2+1/6=0\).
As \(p_0<1/5\) at \(a=3\), this quotient is at most \(5/6\).
Consequently \(\Re F_{xx}\le-1/(12a)\) on the left horizontal.

## Normalized far-contour bounds

Fix \(\rho=|z|>0\) and \(0<\phi=\arg z<\pi\) and put
\(u=\phi+y\) with \(y=\Im x_*(z)\).
Differentiating the saddle relation gives
\[
        u'(\phi)=\Re\left(1-\frac{2a q}{(1+q)^2}\right)^{-1}>0 .
\]
Thus \(u\) increases from zero to \(\pi\) as \(\phi\) does,
by conjugation at both real endpoints. Writing \(r=|q|\), put
\[
  D_q=1+2r\cos u+r^2,\qquad
  y=2a\frac{r\sin u}{D_q}=c\sin u .
\]
Then \(0<y\), and \(c\le21/40\) for \(u\le\pi/2\),
while \(c\le16/25\) for \(u\ge\pi/2\).

For a fixed real \(x\in[b,a]\), put \(p_x=z e^{x+iy-a}\).
The normalized horizontal modulus satisfies
\[
 \partial_\phi\log\frac{|G(x+iy,z)|}{|G(x_*,z)|}
       =(1+y')\left(\Im\frac{q}{1+q}
                        -\Im\frac{p_x}{1+p_x}\right)\le0 .
\]
To verify the sign, \(p_x,q\) have common argument \(u\),
and \(r\le|p_x|\le1/r\) since \(x\ge b>\Re x_*\)
and \(|p_x|\le R_0\). At fixed \(u\in(0,\pi)\) the function
\(s\sin u/(1+2s\cos u+s^2)\) increases for \(s<1\)
and is unchanged under \(s\mapsto 1/s\).
The displayed derivative uses the saddle identity
\(y/(2a)=\Im[q/(1+q)]\).

It follows that the maximum normalized modulus on this piece
at fixed \(\rho\) occurs for positive \(z=\rho\).
For that real \(z\), the spatial maximum on \([b,a]\)
is attained at an endpoint: the derivative sign is nonpositive
precisely when \(\rho e^{-a}\le e^{-x}x/(2a-x)\),
and the right side decreases on the interval because
\(1/x+1/(2a-x)\le600/611<1\).
At \(x=b\), left concavity and \(b-\Re x_*>3/5\) give
\[
 \log |G(b,z)/G(x_*,z)|<-\frac{(3/5)^2}{24a}<-1/300
 \qquad (z=\rho).
\]
At \(x=a\), use \(G(x_*,\rho)\ge G(0,\rho)\ge1\) and
\(G(a,\rho)\le(41/20)e^{-3/4}<399/400\).
This bounds the whole far horizontal piece by \(399/400\)
relative to the saddle.

For the connector at \(x=a+iv\), \(0\le v\le y\), the derivative
of log modulus is \(v/(2a)-\Im[p/(1+p)]\) with
\(p=\rho e^{i(\phi+v)}\). By concavity of sine on \([0,\pi]\),
\(\sin(\phi+v)\ge(v/y)\sin u\). Also
\[
       r/\rho\le1/12
                < (91/100)^2/(1+R_0)^2
                \le D_q/(1+\rho)^2 .
\]
These inequalities and \(y/(2a)=r\sin u/D_q\) give
\(\Im[p/(1+p)]\ge v/(2a)\).
Thus the largest connector modulus is at \(v=0\).

At \(v=0\) its normalized modulus decreases with \(\phi\) also.
Its angular derivative is
\(-\Im[z/(1+z)]+\Im[q/(1+q)]\).
Since \(\sin\phi\ge\sin u-y=(1-c)\sin u\), it suffices
to show \((1-c)D_q/(1+2\rho\cos\phi+\rho^2)\ge r/\rho\).
For \(u\le\pi/2\) the left side is at least
\((19/40)/(1+R_0)^2>1/12\).
For \(u\ge\pi/2\) we have \(\phi\ge\pi/2-16/25\),
so \(\cos\phi\le16/25\). This gives lower bound
\[
       \frac{(9/25)(91/100)^2}
                  {1+R_0^2+2R_0(16/25)}>1/12 .
\]
The inequalities at zero modulus follow by continuity, and the
lower complex semicircle follows by conjugation.
The connector therefore satisfies the same normalized bound
\(399/400\), using its real endpoint at positive \(z\).
Together the far horizontal part and connector have length below \(4\).


## Relative adjacent saddle estimate

Deform the original contour to the horizontal through \(x_*\)
followed by the connector at \(a\). On the symmetric local
interval write \(x=x_*+t\) with \(|t|\le1/4\).
Here \(|ze^{x-a}|\le(9/100)e^{1/4}<3/25\).
The differentiated logarithm gives \(|F_{xxx}|<1/5\)
and \(|F_{xxxx}|<8/25\): the respective modulus bounds
at a given \(|p|\) are
\(|p|(1+|p|)/(1-|p|)^3\) and
\(|p|(1+4|p|+|p|^2)/(1-|p|)^4\).
Put \(F(x_*+t,z)-F_*=-A_*t^2/2+E(t)\) and
\(c_3=F_{xxx}(x_*,z)/6\).
Then \(|E(t)|\le |t|^3/30\),
\(|E(t)-c_3 t^3|\le |t|^4/75\), and
\[
 \frac{\Re A_*}{2}t^2-|E(t)|\ge t^2/75
\]
since \(9/400-1/120>1/75\).

The integral of the cubic term against the matching Gaussian on
the symmetric local interval is zero. The remaining local
difference from the Gaussian is bounded in modulus by the
integral of
\((j|t|^4/75+j^2 |t|^6/1800)e^{-jt^2/75}\).
After normalization by \(\sqrt{j/(4\pi a)}\), the Gaussian
fourth and sixth moments give less than \(1300/j\).
At either local endpoint the real phase loss is at least
\(1/1200\), and the derivative in the outward direction is
at most \(-1/200\): use
\(-\Re A_*/4+(1/5)(1/4)^2/2\le-1/200\).
Left concavity extends the linear loss outward as far as
the left horizontal portion extends.
The same loss bound works for the two tails of the Gaussian.
Each pair of tails after normalization costs at most
\((400/\sqrt{4\pi a})e^{-j/1200}/\sqrt j\).
Using the far-contour bound and \(|C_*|>3/4\) now gives
\[
 I_j=e^{jF_*}C_*(1+\eta_j),\qquad
 |\eta_j|\le \frac{1750}{j}
       +\frac{178}{\sqrt j}e^{-j/1200}
       +6\sqrt j(399/400)^j<1/5\quad(j\ge10000).            \tag{2}
\]

For the adjacent difference use the local integral normalized
by \(e^{jF_*} C_*\) and change variables to \(y=\sqrt j\,t\).
Continue this local expression in real \(j\), with endpoints
\(\pm\sqrt j/4\). The derivative of its exponent in \(j\) is
\(E(t)-tE'(t)/2\), bounded by \(|t|^3/12\).
Using the scaled Gaussian majorant \(e^{-y^2/75}\)
gives the derivative bound \(105/j^{3/2}\) for the
integrand contribution. The summed moving-endpoint
contribution is bounded by \(e^{-j/1200}/(18\sqrt j)\).
Integrate these bounds between \(n\) and \(n+1\), then add
the two orders' actual nonlocal horizontal and far pieces.
There is no Gaussian completion error in this adjacent difference.
The result is
\[
 |\eta_{n+1}-\eta_n|
 \le \frac{105}{n^{3/2}}
       +\frac{180}{\sqrt n}e^{-n/1200}
       +12\sqrt n(399/400)^n<1/1000
                          \quad(n\ge10000).                \tag{3}
\]
All three upper terms decrease on this range.
For (2) and (3) at \(n=10000\), one can use
\(e^{-n/1200}<(7/19)^8\) and
\((399/400)^n\le e^{-25}<2^{-25}\).

## Angular interlacing

Let \(H_j=jF_*+\log C_*+u_j\), where \(u_j=\log(1+\eta_j)\)
is the analytic logarithm defined by (2).
All logarithms are real on the positive radius and have
conjugation symmetry. Their norms satisfy
\(\|u_j\|<1/4\) and \(\|u_{n+1}-u_n\|<1/800\)
on \(|z|\le21/20\) by (2) and (3).
The amplitude is common at fixed \(S\).
On \(|z|=1\), (1) and Cauchy's estimate with radius \(1/20\)
give
\[
 |H_{n+1}'-H_n'| < 1/10+1/40,\qquad
 |H_j'|< j/10+25<j/2 .
\]
As in the
[alternating phase comparison](MULTIVARIATE_DIFFUSE_EXPLICIT_S8_INTERLACING.md#phase-comparison),
the phases \(T_j(\theta)=j\theta-2\Im H_j(e^{i\theta})\)
are strictly increasing from \(0\) to \(2j\pi\),
and \(0<T_{j+1}(\theta)-T_j(\theta)<2\pi\) at every interior
\(\theta\). The half decomposition gives one root at every
odd multiple of \(\pi\) in \(T_j\), proving simplicity and
strict interlacing.

