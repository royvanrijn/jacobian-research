# A direct coefficient-ratio criterion for the first sector

For distinct nonnegative integer exponent columns \(a,b\), write
\[
P_m(z)=\sum_{k=0}^m \binom{m}{k} F_{m,k}z^k,\qquad
F_{m,k}=\prod_i((m-k)a_i+kb_i)!.
\]
Call the profile fully shared if \(a_i,b_i>0\) in every row after omitting zero rows. For a fully shared profile set \(d_i=b_i-a_i\) and \(V=\sum_i d_i^2/\min(a_i,b_i)\).

**Theorem.** Orient a fully shared profile so \(\Delta=\sum_i d_i\ge0\), let \(N=n+1\), and suppose
\[
 n\ge2,\qquad
 V\le\frac{4}{(1+1/n+1/(2N))^2}.
 \tag{1}
\]
If \(P_n(z)=P_N(\kappa z)=0\) and \(\kappa\ge1\), then \(|\arg z|>2\pi/N\).
The same is true under the alternative sufficient condition
\[
 \left((N/n)^2+\frac1{2n}\right)V-2(N/n)\sqrt V\le\Delta.
 \tag{1a}
\]
Thus in the first sector a same-ray upper root must have strictly smaller radius in this orientation. In particular, fully shared profiles with \(V\le16/5\) have no common physical root in the closed first sector in any degree. The first-sector width used here is \(2\pi/N\). The theorem gives no exclusion of the other sectors.

For the last assertion, (1) applies when \(n\ge13\) because its right side increases with \(n\) and \(4(364/405)^2>16/5\). For \(n\le12\), the physical exclusion follows from the [ordinary moment theorem through terminal thirteen](MULTIVARIATE_TERMINAL_THIRTEEN_ALL_PROFILES.md). The degree-dependent claim (1) does not use those finite results.

## One sine crossing for the direct coefficient ratio

Let \(A(z)=\sum_{k=0}^n A_k z^k\) and \(B(z)=\sum_{k=0}^N B_k z^k\) with all coefficients positive. If \(B_k/A_k\) is nondecreasing for \(0\le k\le n\), then \(A,B\) have no common root with \(0<|\arg z|<2\pi/N\). If the quotient is strictly increasing, the exclusion includes \(|\arg z|=2\pi/N\).

Indeed, by conjugation take a root \(z=r e^{i\theta}\) of \(A\) with \(0<\theta\le2\pi/N\). The positive term at \(k=1\) in the imaginary-part sum forces a negative term among \(k\le n\), so \(\pi<n\theta<2\pi\). As the nonzero sines change sign once, from positive to negative, take \(c\) between the last \(B_k/A_k\) with positive sine and the first with negative sine. Then
\[
\Im B(z)=\sum_{k=0}^n\left(\frac{B_k}{A_k}-c\right)
           A_k r^k\sin(k\theta)+B_N r^N\sin(N\theta).
\tag{2}
\]
The sum is nonpositive, and it is strictly negative in the strictly increasing case, by taking \(c\) interior to the interval. The last term is negative for \(\theta<2\pi/N\) and zero at equality. No positive real root is possible.

## The exact factorial ratio

For \(A=P_n\) and \(B(z)=P_N(\kappa z)\), the coefficient ratio in (2) is
\[
 \frac{B_k}{A_k}=\frac{N\kappa^k}{N-k}\,g(k),\qquad
 g(x)=\prod_i\prod_{j=1}^{a_i}(na_i+d_i x+j).
 \tag{3}
\]
The factor \(1/(N-k)\) is retained here. The [polar derivative criterion](MULTIVARIATE_POLAR_DERIVATIVE_FIRST_SECTOR.md) instead compares \(g(k)\) directly. For example, at \(a=(2,3)\), \(b=(3,2)\), the signed total is zero and \(g(n)<g(n-1)\) in both orientations for \(n\ge2\). Indeed, \(g(n)/g(n-1)=(3n+2)(2n+1)/[(3n)(2n+4)]<1\). The direct criterion below applies to this example.

Put \(L=\log(N/n)\) and \(W_+=\sum_{d_i>0}d_i/a_i\). Each row's contribution to \((\log g)'(x)\) is
\[
 d_i R_i,\qquad
 R_i=\frac1{a_i}\sum_{j=1}^{a_i}
        \left(n+\frac{d_i x}{a_i}+\frac{j}{a_i}\right)^{-1}.
\]
Let \(I_i\) denote the integral of \((n+d_i x/a_i+t)^{-1}\) over \(0\le t\le1\). Right endpoint sums for this decreasing convex integrand give
\[
 I_i-\frac{1}{2a_i}\left[
       \left(n+\frac{d_i x}{a_i}\right)^{-1}
       -\left(n+\frac{d_i x}{a_i}+1\right)^{-1}\right]
 \le R_i\le I_i .
\]
Convexity makes the trapezoidal sum an upper bound on the integral, and its difference from the right sum is the displayed endpoint correction. For \(d_i>0\) the correction multiplied by \(d_i\) is at most \(d_i/(2a_i nN)\). For \(d_i<0\), the upper bound on the right sum gives \(d_iR_i\ge d_i I_i\).

For either sign, direct subtraction gives
\[
 d_i(I_i-L)=-\frac{d_i^2 x}{a_i}
       \int_0^1\frac{dt}{(n+d_i x/a_i+t)(n+t)}
       \ge-\frac{d_i^2 x}{n^2\min(a_i,b_i)}.
\]
Here \(0\le x\le n\) gives \(n+d_i x/a_i\ge n\min(a_i,b_i)/a_i\). Summing the bounds proves
\[
 \left(\log\frac{g(x)}{N-x}\right)'
 \ge \Delta L+\frac1{N-x}-\frac{Vx}{n^2}-\frac{W_+}{2nN}.
 \tag{4}
\]
Because every positive \(d_i\) is an integer, \(W_+\le V\). With \(y=x/n\), the last three terms of (4) are nonnegative as soon as
\[
 V\left(y+\frac{1}{2N}\right)\left(1+\frac1n-y\right)\le 1.
\]
The product of the two positive factors is at most one quarter of the square of their sum. Condition (1) therefore suffices. Even at equality in (1), the derivative bound is strictly positive except possibly at the unique point where the two factors are equal. This makes \(g(x)/(N-x)\) strictly increasing on \([0,n]\). Since \(\kappa\ge1\), (3) is strictly increasing, and (2) proves the claim.

To retain the gain from the signed degree difference, set \(t=x/N\),
\(U=V(N/n)^2>0\), and \(C=V/(2n)-\Delta\). Multiply (4) by \(N-x\)
and use \(N\log(N/n)>1\) and \(W_+\le V\). The resulting lower
bound is \(1-(Ut+C)(1-t)\). If \(U+C\le0\), this is at least
one for \(0\le t\le n/N\). If \(U+C>0\), the concave quadratic
\((Ut+C)(1-t)\) has maximum \((U+C)^2/(4U)\) on the real
line. Condition (1a) gives \(U+C\le2\sqrt U\). The derivative is
then positive except possibly at one point, again giving the
strictly increasing ratio in (3).
