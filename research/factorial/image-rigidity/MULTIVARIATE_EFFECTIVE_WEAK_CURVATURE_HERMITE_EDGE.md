# A scaled Hermite spacing bound at the outer edge

Let \(N=n+1\ge3\) and \(\rho\in[n/N,1]\).
The distance between the roots of \(\mathrm{He}_n\) and the
points \(\sqrt\rho\,y\) with \(\mathrm{He}_N(y)=0\) is at least
\[
                           \frac1{3200N^{7/6}}.              \tag{1}
\]
Together with the
[even angular displacement estimate](MULTIVARIATE_EFFECTIVE_WEAK_CURVATURE_ANGULAR_PARITY.md),
this gives disjoint root arguments and simplicity for the two
Gamma moment polynomials whenever
\(\sigma\le(54\cdot10^6 N^{11/6})^{-1}\).
For ordinary factorial moments, it suffices that
\(S\le10^{-16}N^{-8/3}\).

## The largest lower-order root

Let \(J_n\) be the size-\(n\) real symmetric Hermite Jacobi matrix
with zero diagonal and off-diagonal entries \(\sqrt{k}\),
\(1\le k<n\). Its characteristic polynomial is
\(\mathrm{He}_n\) by the Hermite recurrence.
For a real unit vector \(u\), add \(u_0=u_{n+1}=0\) and
write \(a_k=\sqrt k\) for \(0\le k\le n\).
Then
\[
 2\sqrt N-u^\mathsf T J_n u
    =\sum_{k=0}^n a_k(u_{k+1}-u_k)^2
     +\sum_{i=1}^n(2\sqrt N-a_{i-1}-a_i)u_i^2.
\]
All terms on the right are nonnegative.
Set \(L=\lfloor N^{1/3}\rfloor\), so \(N^{1/3}/2\le L\le N/2\).
For \(i\le n-L\), the potential in the second sum is at least
\(2(\sqrt N-\sqrt i)\ge L/\sqrt N\).
If at least half the squared norm is at these indices, the
energy is at least \(L/(2\sqrt N)\ge N^{-1/6}/4\).

Otherwise at least half the squared norm is on the last \(L\)
indices. Telescoping from \(u_{n+1}=0\) and applying
Cauchy–Schwarz gives
\[
       \sum_{i=n-L+1}^n u_i^2
          \le L^2\sum_{k=n-L+1}^n (u_{k+1}-u_k)^2 .
\]
Here \(a_k\ge\sqrt{N-L}\ge\sqrt{N/2}\).
The energy is therefore at least
\(\sqrt{N/2}/(2L^2)\ge N^{-1/6}/(2\sqrt2)\).
The variational formula proves that every lower-order root
\(x\) satisfies
\[
                      x\le2\sqrt N-\tfrac14N^{-1/6}.        \tag{2}
\]

## The adjacent scaled gap

Let \(N\ge4096\) and \(x>0\) be a lower-order root.
It is a stationary point of \(\mathrm{He}_N\), at which
\(\mathrm{He}_N(x)\ne0\) by interlacing. Let \(y>x\) be
the nearest larger-order root to its right.
With \(t=\sqrt N(\xi-x)\), normalize
\[
          g(t)=\mathrm{He}_N(x+t/\sqrt N)/\mathrm{He}_N(x).
\]
Then \(g(0)=1,\ g'(0)=0\), and
\(g''-c(t)g'+g=0\) for \(c(t)=x/\sqrt N+t/N\).
For \(0\le t\le1\), (2) gives \(c(t)\le2\).
Before its first positive zero one has \(g'(t)<0\):
it is negative near zero as \(g''(0)=-1\), and at any
critical point with \(g>0\) the derivative of \(g'\)
is \(-g<0\).
Thus \(p=-g'/g\ge0\) there and
\(p'=1+c(t)p+p^2\). The constant coefficient solution for
\(c=2\) is \(g_*(t)=(1-t)e^t\), giving \(p_*(t)=t/(1-t)\).
While both solutions are defined for \(0\le t<1\),
\[
          (p-p_*)'\le (2+p+p_*)(p-p_*)
\]
and \(p(0)=p_*(0)=0\). Integration gives \(p\le p_*\).
The first positive zero of \(g\) therefore cannot be less
than one, since at a smaller first zero \(p\) would diverge
while \(p_*\) is bounded. Thus \(y-x\ge1/\sqrt N\).

Put \(s=\sqrt{n/N}\). Scaling upwards from \(s\) only
increases the outer gap. By (2),
\[
 \sqrt\rho\,y-x
    \ge s\,y-x
    \ge\frac{s+s^2-x/\sqrt N}{(1+s)\sqrt N}
    \ge\frac{\tfrac14N^{-2/3}-2/N}{2\sqrt N}
    \ge \frac1{16N^{7/6}},
\]
using \(N^{1/3}\ge16\).
For the inner, zero, and negative neighbors, the
[earlier Hermite comparison](MULTIVARIATE_EFFECTIVE_WEAK_CURVATURE_IMPROVEMENT.md#the-hermite-gap)
gives a bound at least \(1/(40\sqrt N)\).
Reflection covers the negative lower-order roots.
For \(N<4096\), that earlier argument gives
\(1/(200N^{3/2})\ge1/(3200N^{7/6})\).
This proves (1) for all \(N\).

## Applying the even displacement bound

The perturbation proof in the even angular displacement note
uses individual disks of radius \(r=1/(80\sqrt N)\).
It bounds the imaginary displacement from each Hermite root by
\(r/(L_v^2-1)\), where \(L_v=(6\cdot10^6\sigma N^{3/2})^{-1}\)
is the complex perturbation disk radius.
Under \(\sigma\le(54\cdot10^6N^{11/6})^{-1}\),
one has \(L_v\ge9N^{1/3}\), and hence \(L_v^2-1>80N^{2/3}\).
Thus
\[
                \frac{r}{L_v^2-1}<\frac1{6400N^{7/6}}.
\]
The sum of the two possible angular displacements is smaller
than the Hermite gap (1). The same midpoint argument therefore
gives disjoint physical root arguments and simplicity.
For ordinary factorial moments under \(S\le10^{-16}N^{-8/3}\),
\(\sigma^2\le2S/n\le3\cdot10^{-16}N^{-11/3}\), which implies
the stated bound on \(\sigma\).
