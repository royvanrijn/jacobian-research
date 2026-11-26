# Physical coprimality for linked Gamma moments and a single shared unit gap

For a positive real shape \(u\), let \(X_u\) denote a unit-scale Gamma random variable, and write

\[
Q_m(w;u)=\mathbb E(w+X_u)^m
        =\sum_{k=0}^m\binom mk(u)_k w^{m-k}.
\]

**Theorem.** For every integer \(n\ge1\), every real \(u>0\), and every integer \(h\ge0\), the polynomials \(Q_n(w;u)\) and \(Q_{n+1}(w;u+h)\) have no common complex zero. The same holds if both Gamma variables have the same arbitrary positive scale.

In particular, all ordinary two-monomial profiles with exactly one unequal exponent row and that row of gap one satisfy SFC(2) in every degree. If the unequal row is \((a,a+1)\) with \(a\ge0\), reversing the factorial moment polynomial and stripping its common factorial factors gives \(Q_m(w;ma+1)\). The theorem applies with \(u=na+1\) and \(h=a\). Swapping the two columns reverses the moment polynomials. All equal exponent rows contribute only nonzero scalar factors. The result concerns common physical zeros. It does not compare roots having the same angle at different radii.

## The numerical range of the upper moment

Put \(N=n+1\), \(b=u+h-1\), \(M=b+N\), and \(Y(w)=Q_N(w;u+h)\). A common zero must be nonreal because an even-order Gamma moment is strictly positive at every real \(w\). The constant coefficients are positive, so every candidate is nonzero.

The coefficient formula gives the monic generalized Laguerre representation

\[
Y(w)=(-1)^N N!L_N^{-M-1}(w),\qquad
wY''-(M+w)Y'+NY=0.
\tag{1}
\]

Fix the Laguerre parameter \(-M-1\) at every matrix size. The monic recurrence has diagonal \(2k-M\) and successive off-diagonal products \(-k(M+1-k)\). Thus \(Y\) is the characteristic polynomial of the complex symmetric tridiagonal matrix with diagonal \(2k-M\) for \(0\le k<N\) and entries \(i\sqrt{k(M+1-k)}\) linking \(k-1\) and \(k\) for \(1\le k<N\). The square roots in this matrix are the positive real ones. The recurrence and differential equation follow from the finite Laguerre coefficient formula. The classical formulas are recorded in [DLMF §18.5](https://dlmf.nist.gov/18.5.E12), [DLMF Table 18.8.1](https://dlmf.nist.gov/18.8.T1), and [DLMF Table 18.9.1](https://dlmf.nist.gov/18.9.T1).

At a zero \(w=-s\) of \(Y\), choose a unit eigenvector \(v\) and put \(\mu=\sum_{k=0}^{N-1} k|v_k|^2\). Taking the Hermitian inner product with \(v\) gives

\[
x:=\operatorname{Re}s=M-2\mu>M-2N+2,\qquad
 |\operatorname{Im}s|
 \le 2\sqrt{\mu(M-\mu)}.
\tag{2}
\]

The strict first bound holds because equality would force the eigenvector to be supported only at \(N-1\), whose off-diagonal entry is nonzero. For the second bound, apply Cauchy–Schwarz to \(\sqrt k|v_k|\) and \(\sqrt{M+1-k}|v_{k-1}|\), using

\[
\sum_{k=1}^{N-1}(M+1-k)|v_{k-1}|^2
\le M-\mu.
\]

Hence \(|s|^2\le (M-2\mu)^2+4\mu(M-\mu)=M^2\).

## The linked lower moment at an upper zero

The exponential generating function \(e^{w t}(1-t)^{-u}\) shows that

\[
Q_n(w;u)=\frac{1}{N}(1-\partial_w)^hY'(w)
        =\frac{(-1)^h}{N}e^w g^{(h)}(w),
\qquad g(w)=e^{-w}Y'(w).
\tag{3}
\]

At the nonzero root \(w=-s\), the differential equation (1) gives \(Y'(w)\ne0\) by uniqueness for a second-order ordinary differential equation. It also gives \(g'(w)/g(w)=M/w\). Thus \(h=0,1\) already give nonzero lower values.

For \(h\ge2\), differentiate (1) and substitute \(Y'=e^w g\) to obtain

\[
w g''+(w-M+1)g'-bg=0.
\]

At \(w=-s\), set \(H_j(s)=w^j g^{(j)}(w)/g(w)\). The differentiated equations give

\[
H_0=1,\qquad H_1=M,\qquad
H_{j+2}=(s+M-1-j)H_{j+1}-s(b-j)H_j.
\tag{4}
\]

Consequently \(H_h(s)=\det(D+sT-qC)\), where \(q^2=s\), \(D=\operatorname{diag}(M-j)_{j=0}^{h-1}\), \(T=\operatorname{diag}(0,1,\ldots,1)\), and \(C\) is the real symmetric tridiagonal matrix with zero diagonal and entries \(\sqrt{b+1-j}\) linking \(j-1\) and \(j\) for \(1\le j<h\). These square roots are positive because \(b+1-j=u+h-j>0\). Either square root \(q\) gives the same determinant by the recurrence (4).

Suppose \(H_h(s)=0\) and let \(v\) now denote a unit kernel vector. Set \(t=\sum_{j=1}^{h-1}|v_j|^2\), \(\ell=\sum_{j=1}^{h-1}j|v_j|^2\), and \(c=v^*Cv\in\mathbb R\). The equation of its quadratic form is

\[
M-\ell+st-qc=0.
\]

Here \(t>0\), since \(t=0\) would give \(M=0\). Dividing by \(q\) and taking imaginary parts gives \(M-\ell=|s|t\) because \(s\) is nonreal and hence \(\operatorname{Im}q\ne0\). Taking real parts then gives \(c=2t\operatorname{Re}q\). Cauchy–Schwarz on the off-diagonal terms shows

\[
 x+|s|=2(\operatorname{Re}q)^2=\frac{c^2}{2t^2}
 \le \frac{2((b+1)t-\ell)}{t^2}.
\]

The upper-zero bound \(|s|\le M\) implies \(\ell\ge M(1-t)\). Therefore

\[
x\le \frac{(2b+2-M)t-\ell(2-t)}{t^2}\le 2b+2-M=M-2N+2.
\]

For the middle inequality, substitute \(\ell\ge M(1-t)\) and subtract \((2b+2-M)t^2\) from the numerator. The result is at most \((1-t)((2b+2)t-2M)\le0\), since \(b+1=u+h>0\), \(b+1\le M\) and \(0<t\le1\). This contradicts the strict inequality in (2). Thus (3) is nonzero at every possible common zero. Common positive scale only rescales \(w\) and multiplies each moment by a nonzero scalar.
