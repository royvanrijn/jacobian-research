# A quadratic pencil for the actual linked beta root geometry

The diagonal matrix in the external
[Vishnyakova proof, Sections 1–3](https://github.com/1337isnot1337/solved-conjectures/blob/758066055ef49c2c66cb1087918ae62ad75e7ac1/results/vishnyakova-root-radius/PROOF.md)
transfers to the actual binary linked moment polynomial.
The transfer gives an exact root modulus in terms of a Rayleigh
quotient, a sector bound that is uniform over moment orders at
fixed positive exponents, and explicit necessary conditions on
a consecutive collision. These are constraints on
`GAP-FACTORIAL-BINARY-LINKED-BETA-COPRIMALITY`.

## 1. The physical polynomial in the quadratic pencil

For \(\alpha,\gamma>0\), put \(\beta=\alpha+\gamma\) and define
\[
\begin{aligned}
 G_m(s;\alpha,\gamma)
    &=\sum_{k=0}^m\binom mk(\alpha)_{m-k}(\gamma)_k s^k,\\
 J_m(t;\alpha,\gamma)
    &=\sum_{k=0}^m(-1)^k\binom mk
                          \frac{(\alpha)_k}{(\beta)_k}t^{m-k}.
\end{aligned}
\]
Here \((a)_k\) is the rising factorial. In the actual linked problem,
\(\alpha=mr+1,\gamma=mu+1\) with fixed \(r,u\in\mathbb N_0\).
The polynomial \(G_m\) is the ordinary factorial moment of
\(X^rY^u(X+sY)\) divided by \((mr)!(mu)!\).
The [linked beta identity](canonical-density/moment_pair_obstructions.md)
or its binomial expansion gives
\[
 J_m(t;\alpha,\gamma)
   =\frac{(t-1)^m}{(\beta)_m}
                 G_m\left(\frac{t}{t-1};\alpha,\gamma\right).
\]
The endpoints are nonzero. All coefficients of \(G_m\)
are positive, so no root has \(s\in[0,\infty)\).

Vishnyakova’s original weighted normalization appears in
[arXiv:2608.03723v1, Open Problem 7.2](https://arxiv.org/pdf/2608.03723v1).
We use the diagonal recurrence of the pinned external proof,
where its z-variable is the argument of a falling factorial.
For temporary variables \(z,t\) define
\[
 p_j(z;t)=(\beta)_j\sum_{k=0}^j\binom jk
                        \frac{z^{\underline k}}{(\beta)_k}t^{j-k}.
\]
Summing first over \(j-k\) gives its exponential generating
series \((1+(1-t)v)^z(1-tv)^{-z-\beta}\).
Multiplying its derivative by
\((1+(1-t)v)(1-tv)\) and comparing coefficients gives
\[
 p_0=1,\quad p_1=z+\beta t,\qquad
 p_{j+1}=[z-j+(2j+\beta)t]p_j
                 +t(1-t)j(j+\beta-1)p_{j-1}.
\]
Thus with \(D=\operatorname{diag}(0,\ldots,m-1)\), and
\(K\) real symmetric with zero diagonal and
\(K_{j-1,j}=K_{j,j-1}=\sqrt{j(j+\beta-1)}\)
for \(1\le j<m\), the leading-block determinant recursion gives
\[
 p_m(z;t)=\det[zI-D+t(2D+\beta I)-i\sqrt{t(1-t)}K].
\]
The square-root choice does not affect the determinant, since
conjugation by the alternating diagonal sign matrix sends
\(K\) to \(-K\). At \(z=-\alpha\),
\(p_m(-\alpha;t)=(\beta)_m J_m(t;\alpha,\gamma)\).
Take \(t=s/(s-1)\) and use
\(\sqrt{t(1-t)}=i\sqrt{s}/(s-1)\).
Multiplication by \((s-1)^m\) now gives the exact quadratic pencil
\[
 G_m(s;\alpha,\gamma)
     =\det\!\bigl[\alpha I+D+\zeta K+s(\gamma I+D)\bigr],
                    \qquad \zeta^2=s.                              \tag{1}
\]
The identity extends to \(s=1\) as a polynomial identity. The disk
theorem for zeros in Vishnyakova’s falling-factorial variable is
distinct from (1). Here we use the underlying matrix before
applying its numerical-range calculation in the other variable.

## 2. Root radius and a pairwise obstruction certificate

Let \(s\) be a root and choose a unit null vector \(v\) of the
matrix in (1). Set \(d=v^*Dv\), \(k=v^*Kv\); both are real.
For \(m>1\) the first and last coordinates of \(v\) are nonzero:
a zero end coordinate would force its neighbor to vanish by
the end row, then every coordinate by induction, since each
offdiagonal entry \(\zeta K_{j-1,j}\) is nonzero.
Hence \(0<d<m-1\), with \(d=0\) at \(m=1\).
The scalar equation is
\[
               (\gamma+d)\zeta^2+k\zeta+\alpha+d=0.
\]
A real \(\zeta\) would give \(s\ge0\), already excluded.
The imaginary part of the equation gives
\(\Re\zeta=-k/(2(\gamma+d))\). Substitution of this into the
real part gives the exact radius
\[
                   |s|=\frac{\alpha+d}{\gamma+d}.                 \tag{2}
\]
For the physical parameters this becomes
\[
 |s|=\frac{r+\delta_m}{u+\delta_m},
                    \qquad \delta_m=\frac{1+d}{m}.
\]
At \(m\ge2\), \(1/m<\delta_m<1\).
For \(r<u\), every root therefore satisfies
\[
          \frac{mr+1}{mu+1}<|s|<\frac{r+1}{u+1}.
\]
For \(r>u\), reverse the endpoint inequalities; for \(r=u\),
all roots have modulus one. For \(r<u\), the inverse chart
\(t=s/(s-1)\) sends the upper radius
\(R=(r+1)/(u+1)<1\) to the disk
\(\lvert t\rvert<R\lvert t-1\rvert\), with rightmost coordinate
\(R/(1+R)=(r+1)/(r+u+2)<1/2\).
Reflection gives the analogous strict leftmost bound when
\(r>u\).

For an unequal physical profile, the function
\(\delta\mapsto(r+\delta)/(u+\delta)\) is strictly monotone.
If two consecutive orders share a root, their two actual null
vectors must therefore satisfy
\[
 \frac{1+d_m}{m}=\frac{1+d_{m+1}}{m+1}=\delta,\qquad
 \frac{k_m}{m}=\frac{k_{m+1}}{m+1}
                    =-2(u+\delta)\Re\sqrt{s},                     \tag{3}
\]
using one square root at both orders. The radii intervals are
nested as the order grows, so the radius bounds allow an
overlap. The remaining comparison must constrain actual
consecutive null vectors beyond their separate numerical
ranges. Real common roots are already excluded by strict
positivity of the even-order \(J_m\) on the real axis.

## 3. Order-uniform angle and compactness

For \(\alpha,\gamma\ge1\) and \(m\ge2\), the same null vector
gives a stronger angular bound. The identity
\[
 (j+\alpha)(j-1+\gamma)-j(j+\beta-1)=\alpha(\gamma-1)
\]
implies, for \(1\le j<m\),
\[
 j(j+\beta-1)\le (1-c_1)(j+\alpha)(j-1+\gamma),
 \qquad c_1=\frac{\alpha(\gamma-1)}
                  {(\alpha+m-1)(\gamma+m-2)}.
\]
Apply Cauchy–Schwarz separately to the two weighted coordinate
lists in \(v^*Kv\). The first omits coordinate zero and the
second omits coordinate \(m-1\), which both have positive
mass. Thus
\[
 |k|<2\sqrt{(1-c_1)(\alpha+d)(\gamma+d)}.
\]
Using \((j+\gamma)(j-1+\alpha)\) instead gives the same bound
with \(c_2=\gamma(\alpha-1)/((\gamma+m-1)(\alpha+m-2))\).
With \(c=\max(c_1,c_2)\), equation (2) and the formula for
\(\Re\zeta\) yield
\[
 \sin^2\!\left(\frac{\arg s}{2}\right)
 =1-\frac{k^2}{4(\alpha+d)(\gamma+d)}>c .                           \tag{4}
\]
For every fixed positive \(r,u\), each of \(c_1,c_2\) exceeds
\(c_\infty=ru/((r+1)(u+1))>0\).
For example, the difference \(c_1-c_\infty\) is
\[
         \frac{u(r+u+1)}
          {(r+1)(u+1)(m(u+1)-1)}>0 .
\]
Thus the forbidden sector in the native coefficient variable
has angle bounded away from zero uniformly in \(m\).
If one exponent is zero and the other positive, (4) still has
\(c>0\) at each \(m\), though this uniform lower bound is zero.

The [Jacobi comparison for the physical pair](LINKED_BETA_UNEQUAL_PROFILE_PIVOTS.md)
uses \(\xi=(1+t)/(1-t)=1-2s\). For positive \(r,u\), set
\(\eta=\min\{r/u,(r+1)/(u+1)\}>0\) and
\(R=\max\{r/u,(r+1)/(u+1)\}\).
All transformed roots at orders \(m\ge2\) lie in the image
of the fixed compact annular sector
\(\eta\le |s|\le R\) with
\(\sin^2(\arg s/2)\ge c_\infty\).
In particular they remain uniformly away from the weight
cut \((-\infty,1]\) in that Jacobi coordinate, and
\[
             |\xi-1|\ge2\eta,\qquad
             |\xi+1|\ge4\sqrt{\eta c_\infty}.
\]
The second estimate follows from
\(|1-s|^2=(1-|s|)^2+4|s|\sin^2(\arg s/2)\).
This bounds the root search in a fixed region away from
the auxiliary real cut. The moving-arc and Airy endpoint
comparisons for adjacent powers remain separate obligations.

## Verification

For each of the eight profiles
\((0,0),(0,2),(3,0),(1,1),(1,2),(2,1),(3,11),(17,20)\)
and orders \(1,\ldots,9\), the exact determinant recurrence
\[
 D_{j+1}(s)=[\alpha+j+s(\gamma+j)]D_j(s)
                      -s j(j+\beta-1)D_{j-1}(s)
\]
from (1), with fixed physical shapes for the selected order,
matched every coefficient of the direct binomial expansion.
These 72 checks used exact integer polynomial arithmetic.
As a separate numerical diagnostic, 360 roots computed at
60 decimal digits were used to reconstruct a tridiagonal
null vector recursively starting with its first coordinate
equal to one. The last-row residual and the exact radius
identity agreed at the stated test tolerances below. Every
tested root of order at least two satisfied the strict
angular bound. These finite checks supplement the determinant
identity and eigenvector proof.

```python
from math import comb, prod
from sympy import Poly, Symbol
import mpmath as mp

s = Symbol("s")
mp.mp.dps = 75
profiles = [(0,0),(0,2),(3,0),(1,1),(1,2),(2,1),(3,11),(17,20)]

def rising(a, n):
    return prod(a+j for j in range(n))

for r,u in profiles:
    for m in range(1, 10):
        a,b = m*r+1,m*u+1
        beta = a+b
        original = Poly(sum(comb(m,j)*rising(a,m-j)*rising(b,j)*s**j
                            for j in range(m+1)),s)
        lower,upper = Poly(1,s),Poly(a+b*s,s)
        for j in range(1,m):
            lower,upper = (
                upper, Poly(a+j+(b+j)*s,s)*upper
                       - j*(j+beta-1)*Poly(s,s)*lower
            )
        assert upper == original
        for root in original.nroots(n=60,maxsteps=500):
            re,im = root.as_real_imag()
            tau = mp.mpc(str(re.evalf(60)),str(im.evalf(60)))
            zeta = mp.sqrt(tau)
            vector = [mp.mpc(1)]
            for j in range(m-1):
                row = (a+j+(b+j)*tau)*vector[j]
                if j:
                    row += zeta*mp.sqrt(j*(j+beta-1))*vector[j-1]
                vector.append(-row/(zeta*mp.sqrt((j+1)*(j+beta))))
            last = (a+m-1+(b+m-1)*tau)*vector[-1]
            if m>1:
                last += zeta*mp.sqrt((m-1)*(m+beta-2))*vector[-2]
            assert abs(last) < mp.mpf("1e-46")*max(map(abs,vector))
            norms = [abs(w)**2 for w in vector]
            d = sum(j*w for j,w in enumerate(norms))/sum(norms)
            assert mp.almosteq(abs(tau),(a+d)/(b+d),rel_eps=mp.mpf("1e-50"))
            if m>1:
                c = max(mp.mpf(a)*(b-1)/((a+m-1)*(b+m-2)),
                        mp.mpf(b)*(a-1)/((b+m-1)*(a+m-2)))
                assert 0 < d < m-1
                assert mp.sin(mp.arg(tau)/2)**2 > c
```
