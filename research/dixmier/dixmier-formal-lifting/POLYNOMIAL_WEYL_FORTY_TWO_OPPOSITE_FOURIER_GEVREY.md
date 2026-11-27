# Divergence in the opposite Fourier equations at degree forty-two

**Theorem PolynomialWeylFortyTwoOppositeFourierGevrey.**
Take the centered hypothetical degree-forty-two pair \(P,Q\) with
\([Q,P]=1\). Represent it in the opposite Fourier coordinate by
\(X=-D_z,\ Y=z\), where \(D_z=d/dz\). For every
\(\lambda\in\mathbb C\), the two opposite Fourier spectral
equations have formal Taylor solution dimensions \(2\) and \(7\)
at \(z=0\). Every nonzero Taylor solution has exact Gevrey order
\(1/2\), and neither equation has a nonzero meromorphic solution
germ at zero. Their polynomial cokernels have dimensions
\(18\) and \(63\), respectively, at every spectral value.

Here order \(1/2\) means that the coefficient \(a_n\) is bounded
in absolute value by \(C A^n(n!)^{1/2}\) for some \(C,A>0\),
and admits no such bound with an exponent smaller than \(1/2\).

## The two ends in the Fourier coordinate

Both operators can be treated together. Let \(T\) denote \(P\)
or \(Q\) in the opposite representation, and put
\[
   (c,d,M)=(2,18,12)\quad\text{or}\quad(7,63,42)
\]
respectively. Here \(M\) is the differential order in \(z\).
Write \(p_{ij}\) for its coefficients before the representation
in the ordered monomials \(X^iY^j\), and put \(C_z=zD_z\).
With the falling factorial notation \((C_z+j)_i\),
\[
          T=\sum_{i,j}(-1)^i p_{ij} z^{j-i}(C_z+j)_i .
                                                        \tag{1}
\]
By the [crossing-face theorem](POLYNOMIAL_WEYL_BERNSTEIN_FORTY_TWO_REAL_BOUND.md),
\(3i-j\le 3c\) on the support, and the coefficient \(p_{c0}\)
is nonzero. If \(i-j\ge c\), then \(i\ge c\) and
\(3i-j=2i+(i-j)\ge3c\). Equality requires \(i=c,\ j=0\).
Thus \((c,0)\) supplies the unique maximal lowering of the
Taylor exponent. The triangular recurrence with diagonal
\((-1)^c p_{c0}(n)_c\) for \(n\ge c\) leaves the first
\(c\) Taylor coefficients free. The rows for \(n<c\)
vanish because \(T-\lambda\) acts on the formal Taylor ring.

The [regular finite-pole theorem](POLYNOMIAL_WEYL_FORTY_TWO_REGULAR_FINITE_POLE.md)
gives \(j-i\le d\). The support corner \((M,3(M-c))\)
satisfies \(j-i=d\). It is the only support point with
\(i=M\), since the crossing bound gives \(j\ge3(M-c)\)
there and the support rectangle gives \(j\le3(M-c)\).
Consequently the highest differential coefficient of \(T\) is
a nonzero scalar times \(z^{3(M-c)}\), so zero is its only
possible finite singularity. The equation at infinity is
regular singular. Indeed \(z^{-d}(T-\lambda)\) in (1) has
coefficients holomorphic in \(1/z\), and its \(C_z^M\)
coefficient has a nonzero value at \(1/z=0\).

## Polynomial eigenvectors and meromorphic germs

On polynomials, \(T-\lambda\) raises degree by at most \(d\)
and by exactly \(d\) above some finite integer \(N_0\).
Indeed the coefficient for the shift \(d\) in (1) is a
nonzero polynomial in \(n\) of degree \(M\), with nonzero
leading term supplied by the corner. In particular every
polynomial killed by a power of \(T-\lambda\) has degree
at most \(N_0\), since a polynomial of higher degree increases
in degree by \(d\) under every iteration.

The Weyl mate acts on polynomials and has a nonzero scalar
commutator with \(T\). If \(v\ne0\) were a polynomial
eigenvector, repeated applications of the mate \(S\) would
produce an infinite linearly independent family
\(v,Sv,S^2v,\ldots\) in the generalized eigenspace. In fact if
\([T,S]=a\ne0\) then \((T-\lambda)S^k v=ka S^{k-1}v\).
This contradicts the degree bound, so \(T-\lambda\) is
injective on polynomials. For \(N\) sufficiently large its
image intersects the polynomials of degree at most \(N+d\)
in the image of polynomials of degree at most \(N\).
Every polynomial of higher degree reduces to degree at most
\(N+d\) by subtracting images of monomials. Thus the cokernel
dimension is \((N+d+1)-(N+1)=d\).

A nonzero meromorphic germ at zero would have initial
exponent \(n\) with \((n)_c=0\) by (1), hence would be
holomorphic there. Continuation gives an entire solution
because zero is the only possible finite singularity and
the germ is single-valued around it. Regular singularity at
infinity gives a polynomial growth bound for this entire
solution. For example the Euler companion system has a
bounded coefficient matrix for the radial Euler derivative
on \(|z|\) sufficiently large. Radial integration gives a
power bound uniformly in angle using a fixed compact circle
of initial data. Cauchy estimates then make the solution
a polynomial, which has just been excluded.

## Exact Gevrey order

For a formal Taylor solution \(\sum_{n\ge0}a_n z^n\)
compare the coefficient at \(z^{n-c}\), and normalize
\(b_n=a_n/\sqrt{n!}\). A term \((i,j)\) uses the index
\(n-\delta\), where \(\delta=c+j-i\). Include the spectral
term by replacing the coefficient of \(X^0Y^0\) with
\(p_{00}-\lambda\). Every term other than \((c,0)\)
has \(1\le\delta\le D=c+d=2(M-c)\).
For sufficiently large \(n\), its normalized factor before
multiplication by the fixed operator coefficient is
\[
       \frac{(n-\delta+j)_i}{(n)_c}
                    \sqrt{\frac{(n-\delta)!}{n!}}
        = n^{i-c-\delta/2}\bigl(1+O(n^{-1})\bigr).
                                                        \tag{2}
\]
The crossing inequality gives \(\delta\ge2(i-c)\).
Thus all recurrence coefficients in
\(b_n+\sum_{\delta=1}^D h_\delta(n)b_{n-\delta}=0\)
have finite limits. For lag \(D\) the only term with
equality in (2) is the corner with \(i=M\). Its coefficient
is nonzero, so the limiting \(D\)-step companion matrix
is invertible.

Both the companion matrices and their inverses therefore
have bounded norms for sufficiently large \(n\).
The vectors \((b_n,\ldots,b_{n-D+1})\) are nonzero
throughout this tail for a nonzero solution, since a zero
state would give an eventually zero series and a polynomial
eigenvector. Their norms are bounded above by an
exponential and below by a positive constant times a
decaying exponential, by iterating the forward and inverse
norm bounds. The upper bound gives Gevrey order \(1/2\).
Any estimate \(|a_n|\le C A^n(n!)^\theta\) with
\(\theta<1/2\) would instead make the entire state
decay faster than every exponential. This contradicts the
lower bound and proves exactness.

More explicitly, write the original crossing symbol as
\(\alpha (X L(XY^3))^c\) for the relevant nonzero scalar
\(\alpha\). The limiting lag polynomial
\(1+\sum_{\delta=1}^D h_\delta(\infty)\xi^\delta\)
is \(L(-\xi^2)^c/L(0)^c\). Indeed equality in (2)
occurs on the crossing line and gives
\((i,j)=(c+b,3b)\), with lag \(2b\) and sign \((-1)^b\)
relative to the leading Taylor term. The two nonzero
crossing roots therefore give the four nonzero roots of
this lag polynomial with their multiplicities.
