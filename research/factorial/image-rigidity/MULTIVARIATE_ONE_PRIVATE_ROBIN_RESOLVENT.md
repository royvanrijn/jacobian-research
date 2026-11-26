# A positive resolvent test for an ordinary shared unit gap and one private unit

For \(N\ge2\), \(\alpha>0\), write
\[
 B_{N,\alpha}(x)=
     \sum_{k=0}^N\binom Nk(\alpha)_k k!\,x^{N-k}
     =\mathbb E(x+G_\alpha Y)^N ,
\]
where \(G_\alpha\) is unit-rate Gamma of shape \(\alpha\) and \(Y\)
is an independent unit exponential. Suppose \(\alpha>N-1\).
Then \(B_{N,\alpha}\) is squarefree. Moreover,
\(B_{N-1,\alpha}\) and \(B_{N,\alpha+r}\) are coprime for
\(r=1,2,3\). They are also coprime for \(r=4\) on the ordinary
path \(\alpha=4(N-1)+1\). Thus the original binomial
\(X^a(x+XY)\), with \(a=1,2,3,4\), satisfies SFC(2) in every
degree. The same test below is defined for every positive integer
shape advance \(r\).

The comparison uses a positive selfadjoint operator to control the
value of \(B_{N,\alpha}\) at a nonreal zero of its derivative. It
retains the exact shared Gamma shape and the private exponential.

## The Robin form

Set \(B=B_{N,\alpha}\) and \(\rho=N+\alpha>2N-1\).
The coefficient recurrence of \(z^N B(1/z)\) is
\(c_{k+1}=(N-k)(\alpha+k)c_k\). Equivalently
\[
                 xB''+(1-\rho)B'+B=x^N.                       \tag{1}
\]
On the Hilbert space \(\mathcal H=L^2([1,\infty),t^{-\rho}dt)\)
use the form, linear in its first argument,
\[
  \mathfrak a(v,w)=\int_1^\infty t^{1-\rho}v'(t)\overline{w'(t)}\,dt
                                      +N v(1)\overline{w(1)}.
\]
Its domain consists of locally \(H^1\) functions in \(\mathcal H\)
with finite displayed derivative integral. This is a dense, closed,
nonnegative form. Completeness in the form norm follows by local
\(H^1\) completeness and convergence in the two weighted spaces;
near \(1\) the local \(H^1\) convergence also gives trace convergence.
Smooth functions of compact support on \([1,\infty)\) are a form
core. Indeed cutoff at \(R,2R\) with derivative bounded by \(C/R\)
costs at most \(O(R^{-1})\int_R^{2R}t^{-\rho}|v|^2dt\) in the
derivative term, and compactly supported functions can then be
smoothed in the local Sobolev norm.

The form therefore defines a nonnegative selfadjoint operator
\(A\), with formal expression
\[
             A_{\rm diff}v=-t^\rho(t^{1-\rho}v')'.
\]
In particular a polynomial of degree at most \(N-1\) satisfying
\(v'(1)=N v(1)\) belongs to \(\operatorname{Dom}(A)\) and \(Av\)
is given by this expression. This follows by integration by parts
against the form core and the representation theorem, since the
formal derivative is again in \(\mathcal H\). The sign of the
Robin condition follows from the lower integration boundary.
The kernel of \(A\) is zero, because zero form energy forces a
constant with value zero at \(1\).
[Kato, chapter VI, §2, first representation theorem](https://link.springer.com/book/10.1007/978-3-642-66282-9)

Put \(h(t)=1/N\), so \(\mathfrak a(v,h)=v(1)\) throughout the
form domain. For each \(1\le k\le N-1\), \(t^k\) belongs to that
domain and integration by parts gives
\[
 \mathfrak a(t^k,w)
   =k(\rho-k)\langle t^{k-1},w\rangle
                         +(N-k)\overline{w(1)}.
\]
Hence \(t^k-(N-k)h\) is in the operator domain, with image
\(k(\rho-k)t^{k-1}\). In particular \(t^{k-1}\) is in the
range of \(A\). By injectivity and functional calculus this gives
\[
 t^k=f_k(A)h,\qquad
 f_0(\lambda)=N,\qquad
 f_k(\lambda)=N-k+\frac{k(\rho-k)}{\lambda}f_{k-1}(\lambda)
                 \quad(\lambda>0).
\]
Each equality is in \(\mathcal H\). The possibly unbounded inverse
is applied only to a vector already shown to lie in its domain.
Write \(f=f_{N-1}\). It has positive coefficients as a polynomial
in \(\lambda^{-1}\), so \(f(\lambda)/\lambda\) is strictly decreasing
on \((0,\infty)\). If \(\mu\) is the spectral measure of \(h\), then
\[
       \int\lambda\,d\mu(\lambda)=\mathfrak a(h,h)=1/N,\qquad
       \int f(\lambda)\,d\mu(\lambda)
                           =\langle t^{N-1},h\rangle=1/(N\alpha).
                                                               \tag{2}
\]
The measure has no atom at zero because \(A\) is injective,
and it is not concentrated at one point. Indeed such a
concentration would put \(h\) in the operator domain. Testing
the form of \(h\) against functions supported away from \(1\)
would give \(Ah=0\), whereas \(\mathfrak a(h,h)>0\).

## The value at a critical point

Let \(\xi\) be any nonreal root of \(B'\). Since \(B\) is monic,
\[
 v(t)=\frac{B(\xi t)}{\xi^N}-t^N
\]
has degree at most \(N-1\) and is in the form domain. Put
\(y=B(\xi)/\xi^N\). Then \(v(1)=y-1,\ v'(1)=-N\). By (1),
\[
        (\xi-A_{\rm diff})v=N\alpha t^{N-1}=F_0.
\]
The boundary defect \(N v(1)-v'(1)\) equals \(\delta=Ny\).
Thus \(v-\delta h\) is in the operator domain, with image
\(A_{\rm diff}v\), by the same integration by parts. As
\(A\) is selfadjoint, \(R_\xi=(\xi-A)^{-1}\) exists, and solving
the resulting equation gives
\[
                   v=R_\xi F_0-\delta A R_\xi h.
\]
Both terms on the right are in the form domain. In the second
term \(A R_\xi\) is a bounded spectral multiplier, so it preserves
the form domain of \(h\). Evaluation at \(1\) using
\(\mathfrak a(v,h)\) now gives
\[
 y-1=N\alpha\int\frac{\lambda f(\lambda)}{\xi-\lambda}\,d\mu(\lambda)
                -Ny\int\frac{\lambda^2}{\xi-\lambda}\,d\mu(\lambda).
\]
Both integrals converge by the form norm and \(f(A)h\in\mathcal H\).
Using (2) and \(\xi\ne0\), this simplifies to
\[
          y=\alpha\,
             \frac{\int f(\lambda)(\xi-\lambda)^{-1}\,d\mu(\lambda)}
                  {\int \lambda(\xi-\lambda)^{-1}\,d\mu(\lambda)} .
                                                               \tag{3}
\]
The denominator has nonzero imaginary part. Multiply the
numerator by the conjugate of the denominator. Its imaginary
part is \(\Im\xi\) times
\[
 \frac12\iint\!
 \frac{(\nu-\lambda)[f(\lambda)\nu-f(\nu)\lambda]}
      {|\xi-\lambda|^2|\xi-\nu|^2}
       \,d\mu(\lambda)\,d\mu(\nu).
\]
This integral is strictly positive. The integrand is
nonnegative by the strict decrease of \(f(\lambda)/\lambda\),
and it is positive for distinct positive \(\lambda,\nu\).
The support property of \(\mu\) ensures positive product mass
off the diagonal. The finite first moments in (2) justify
the product integrals. Thus \(y\) and \(\xi\) have imaginary
parts of the same strict sign.

In particular \(\xi\) cannot be a root of \(B\). A real repeated
root is impossible as well, since the even polynomial \(B\)
or the even polynomial \(B'/N\) is a strictly positive even
moment on the real line.

## A finite sign criterion for a shape advance

Write \(E=x\partial_x\). Up to a nonzero scalar, the upper
shape-\(\alpha+r\) polynomial is \((\rho-E)_r B\), whereas the
lower polynomial is \(B'/N\). The coefficient identities are
the Gamma shape advances in the
[linear-monomial Euler bridge](MULTIVARIATE_LINEAR_MONOMIAL_EULER_BRIDGE.md).
At a lower root \(\xi\) the upper zero equation has the form
\[
                 P_r(\xi)y+\xi V_r(\xi)=0 ,
\]
with \(P_0=1,\ P_1=\rho,\ V_0=V_1=0\), and
\[
 \begin{split}
 P_{j+2}(x)&=(\rho+j+1)P_{j+1}(x)-xP_j(x),\\
 V_{j+2}(x)&=(\rho+j+1)V_{j+1}(x)-xV_j(x)+(\alpha)_j .
 \end{split}                                                    \tag{4}
\]
Indeed \((1-u)^{-\rho}B((1-u)x)\) has \(j\)-th \(u\)-derivative
at zero equal to \((\rho-E)_j B(x)\). Equation (1) transforms
into
\[
 (1-u)\partial_u^2 F-(\rho+1)\partial_u F+xF
                                  =x^{N+1}(1-u)^{-\alpha},
 \qquad F(u)=(1-u)^{-\rho}B((1-u)x).
\]
Taking derivatives at zero gives (4) at a critical point.

The determinant of the size-\(r\) tridiagonal matrix with diagonal
\(\rho,\rho+1,\ldots,\rho+r-1\) and off-diagonal entries
\(-\sqrt{x}\) is \(P_r(x)\). Conjugating by the inverse square
root of the positive diagonal makes the determinant proportional
to that of \(I-\sqrt{x}J\), where \(J\) is real symmetric,
tridiagonal with zero diagonal and positive off-diagonal entries.
Its eigenvalues are simple by the three-term recurrence, and
occur in opposite pairs by alternating sign conjugation.
Consequently the zeros of \(P_r\) are positive and simple.
Also \(\deg(x V_r)\le\lfloor r/2\rfloor=\deg P_r\) by (4).
If every residue of \(-x V_r/P_r\) at the zeros of \(P_r\) is
nonnegative, this rational function has imaginary part of the
opposite weak sign to \(\Im x\). Its polynomial part is real
constant. Equation (3) then excludes a nonreal common root.
Even-moment positivity excludes a real one.

For \(r=1\) the ratio is zero. For \(r=2,3\), (4) gives
\(-xV_r/P_r=\kappa x/(x-c)\), with \(\kappa,c>0\).
In detail, \(c=\rho(\rho+1)\) and \(\kappa=1\) for \(r=2\),
and \(c=\rho(\rho+2)/2\) and
\(\kappa=1-N/[2(\rho+1)]\) for \(r=3\). These pass the
residue test.

For \(r=4\), put \(C=(\alpha+2)(\alpha+\rho+2)+(\rho+1)(\rho+2)\).
Then \(V_4=C-x\) and
\[
 P_4=x^2-3(\rho+1)(\rho+2)x+\rho(\rho+1)(\rho+2)(\rho+3).
\]
The residues of \(x(x-C)/P_4(x)\) are positive if \(C\) lies
strictly between the two roots of \(P_4\). On the ordinary
path \(\alpha=4(N-1)+1\), direct substitution gives
\(P_4(C)=-N^2(229N^2-115N+9)<0\), so the criterion applies.

For the phase \(X^a(x+XY)\), division of each degree-\(m\)
moment by \((ma)!\) gives \(B_{m,ma+1}(x)\). Thus the relevant
arguments above are \(\alpha=a(N-1)+1\) and \(r=a\).
Pure monomials at either coefficient endpoint have nonzero
moments. Rational Bézout identities give the original results
over every characteristic-zero field.
