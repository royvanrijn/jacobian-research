# Sharpness of the private-unit angle boxes for fixed exchangeable bases

Fix positive degrees \(n,N=n+1\) and an angle
\(0<\theta<\pi\) that belongs to a positive-angle box
\[
       I_{m,k}=\left((2k-1)\pi/m,\ 2k\pi/(m+1)\right)
\]
at each of \(m=n,N\), with a label available at both degrees.
There exists a fixed positive exchangeable base \((U,V)\)
with just three atoms such that, after adding independent
unit-rate exponential variables \(X,Y\) as a private unit
pair, both
\[
                     \mathbb E[(UX+zVY)^m],\qquad m=n,N
\]
vanish at \(z=e^{i\theta}\). Thus the overlapping boxes in
the [private-unit angle theorem](MULTIVARIATE_PRIVATE_UNIT_ASYMMETRIC_ANGLE_BOXES.md)
are sharp for consecutive coprimality over positive exchangeable
bases. A common angle is possible already at degrees three
and four. The base used here has finite support and is outside
the ordinary exponent-row class.

## Construction at any overlapping angle

For \(0<\varepsilon<1\) put
\[
     G_m(z)=\sum_{k=0}^m z^k,\qquad
     H_{m,\varepsilon}(z)
        =\frac12\sum_{k=0}^m(\varepsilon^k+\varepsilon^{m-k})z^k.
\]
At \(z=e^{i\theta}\) define the centered real values
\(g_m=e^{-im\theta/2}G_m(e^{i\theta})\) and
\(h_m=e^{-im\theta/2}H_{m,\varepsilon}(e^{i\theta})\).
The ordinary geometric sum gives
\[
       g_m=\frac{\sin((m+1)\theta/2)}{\sin(\theta/2)},\qquad
          h_m\longrightarrow\cos(m\theta/2)
                                      \quad(\varepsilon\to0).
\]
On \(I_{m,k}\) the first sign is \((-1)^{k-1}\) and
the limiting second sign is \((-1)^k\). They are strict.
Choose one sufficiently small positive \(\varepsilon\) for
both degrees. Then \(R_m=-g_m/h_m>0\) at both orders.

Put \(B=R_N/R_n>0\) and choose \(p\in(0,1)\) by
\[
                        \frac{1-p}{p}=R_n B^{-n}.
\]
Give \((U,V)\) probability \(p\) at \((1,1)\) and
probability \((1-p)/2\) at each of
\((B,\varepsilon B),(\varepsilon B,B)\).
The same base is used in both moments.
Indeed integration over \(X,Y\) cancels the binomial
coefficient and yields
\[
   \frac1{m!}\mathbb E[(UX+zVY)^m]
                    =pG_m(z)+(1-p)B^m H_{m,\varepsilon}(z).
\]
By construction \((1-p)B^m/p=R_m\) at both orders.
The claimed two zeros follow. This argument explains why using
one fixed law at both degrees still permits independent
adjustment of the relative masses of the two base orbits.

## An exact rational example

Take \(n=3\), \(\varepsilon=1/10\), and \(z=(5+12i)/13\).
Its angle lies in \(I_{3,1}\cap I_{4,1}\) since its cosine
is strictly between \(\cos(2\pi/5)=(\sqrt5-1)/4\) and \(1/2\).
Writing \(x=\cos(\theta/2)\), the centered values are
\[
            g_3=\frac{20}{13}x,\qquad
            h_3=-\frac{1573}{13000}x,\qquad
            g_4=\frac{61}{169},\qquad
            h_4=-\frac{1107569}{1690000}.
\]
Use \(R_3=20000/1573\) and \(R_4=610000/1107569\) in the
construction. The resulting atom coordinates and probabilities
are all positive rationals. No limiting zero or numerical
solution is involved.
