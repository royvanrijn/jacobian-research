# All-order linked beta coprimality on the one-private boundary

Let \(a\) be a nonnegative integer and let \(\gamma>0\).
For every \(j\ge1\) define the truncated negative-binomial
polynomial
\[
                 P_j(z)=\sum_{k=0}^j
                             \frac{(aj+\gamma)_k}{k!}\,z^k.
\]
Then \(P_n,P_{n+1}\) are coprime at every order \(n\ge1\).
At \(\gamma=1\) this proves SFC(2) for
\[
                         U X^a(CY+E X),
\]
where \(U\) is an arbitrary monomial in variables independent
of \(X,Y\), and \(C,E\) are coefficients in any
characteristic-zero field. This closes the boundary \(r=0\)
or \(u=0\) of the
[linked binary beta problem](LINKED_BETA_VISHNYAKOVA_PENCIL.md#1-the-physical-polynomial-in-the-quadratic-pencil).
For general \(\gamma>0\) the shared variable \(X\) has Gamma
shape \(\gamma\) and the private variable \(Y\) has shape one.

Indeed, for independent unit-rate Gamma variables of these shapes,
\[
 \mathbb E\bigl(X^a(Y+zX)\bigr)^j
      =j!(\gamma)_{aj}\,P_j(z).
\]
The factorial of \(Y\) cancels the corresponding binomial
denominator. Multiplication by \(U\) only introduces a positive
factor at each order. The even-order polynomial has no real root
since it is the expectation of a nonzero even power of a real
linear combination times a strictly positive factor.

## A difference of two beta survival polynomials

If \(a=0\), the adjacent polynomials are consecutive truncations
of one series, so \(P_{n+1}-P_n\) is a nonzero multiple of
\(z^{n+1}\). Neither polynomial vanishes at zero.

Assume \(a\ge1\). For \(A>0\) and a positive integer \(b\) put
\[
 S_{A,b}(z)=(1-z)^A
        \sum_{k=0}^{b-1}\frac{(A)_k}{k!}z^k.
\]
Use the branch on \(\operatorname{Re}z<1\) that is positive
for \(1-z>0\). On the real interval \(0<z<1\) this is the
beta survival function
\[
 S_{A,b}(z)={1\over \mathrm B(A,b)}
                    \int_z^1 t^{b-1}(1-t)^{A-1}\,dt,
\]
as seen by differentiating the finite sum and comparing the
endpoint at \(z=1\). Here \(\mathrm B\) is the complete beta
integral. The two contiguous identities are
\[
 \begin{aligned}
 S_{A+1,b}-S_{A,b}
       &=-\frac{z^b(1-z)^A}{A\,\mathrm B(A,b)},\\
 S_{A,b+1}-S_{A,b}
       &= \frac{z^b(1-z)^A}{b\,\mathrm B(A,b)}.
 \end{aligned}
\]
The first follows by differentiation and the beta recursion
\(\mathrm B(A+1,b)=A\,\mathrm B(A,b)/(A+b)\), and the
second by the additional term in the finite sum.
They extend to the indicated complex half-plane by
analytic continuation.

Set \(A=an+\gamma\), \(b=n+1\), and shift \(A\) up \(a\)
times before increasing \(b\) by one. With \(v=1-z\), the
difference becomes
\[
 S_{A+a,b+1}(z)-S_{A,b}(z)
  = z^b v^A
       \left[\frac{v^a}{b\,\mathrm B(A+a,b)}
           -\sum_{r=0}^{a-1}\frac{v^r}{(A+r)\mathrm B(A+r,b)}\right].
                                                               \tag{1}
\]
If \(P_n,P_{n+1}\) share a nonreal root in the domain,
the square bracket vanishes there. Dividing by \(v^a\)
makes this the equality of a positive real number and a
linear combination of \(v^{-1},\ldots,v^{-a}\) with
strictly positive real coefficients.

## The root-radius bound forbids that equality

The consecutive coefficient ratios of \(P_n\) are
\(k/(an+\gamma+k-1)\) for \(1\le k\le n\). These increase
because \(an+\gamma>1\), and the Eneström--Kakeya upper
radius is
\[
             R_n=\frac{n}{(a+1)n+\gamma-1}<\frac1a .
\]
For completeness, the positive coefficient ratios give
this radius by scaling \(z=R_n w\), after which the
coefficients are nondecreasing. Multiplication by \(1-w\)
gives one negative top term whose magnitude strictly
dominates the other terms at \(|w|>1\).

In particular \(R_n<1\), so all roots are in the branch
domain of (1). If \(|z|\le R_n\), the disk for \(v=1-z\)
has \(|\arg v|\le\arcsin R_n\). Convexity of arcsine on
\([0,1]\) gives
\[
                   a|\arg v|\le a\arcsin R_n<\pi/2 .
\]
When \(z\) is nonreal, every \(v^{-k}\), \(1\le k\le a\),
has imaginary part of the same nonzero sign. Thus the
positive combination in (1) cannot be real. Real common
roots were excluded by the even moment, which proves
coprimality for all consecutive orders. In the ordinary
case the polynomials have rational coefficients, and the
resulting rational Bézout identity gives the field claim.
