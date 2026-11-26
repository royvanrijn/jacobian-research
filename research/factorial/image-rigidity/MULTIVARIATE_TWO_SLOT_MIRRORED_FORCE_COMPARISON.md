# Two-slot mirrored force and the upper polar comparison

Consider either two shared mirrored unit-gap pairs or one
shared mirrored gap-two pair. Each fixed-order polynomial is
a circular coefficient convolution of exactly two Gegenbauer
arrays. At each adjacent order, the actual upper polynomial
strictly interlaces the lower-order polynomial with its
two Gegenbauer shapes advanced to the upper values. Each
native symmetric shape advance moves positive-half lower
roots toward \(\pi\), yielding one side of all the desired
ordinary adjacent angular inequalities.

## Negative force for a two-slot convolution

For real \(\lambda_1,\lambda_2>0\), write
\[
 G(z)=\sum_{k=0}^n\binom nk
       \prod_{j=1}^2
         \frac{(\lambda_j)_k(\lambda_j)_{n-k}}{(\lambda_j)_n}z^k .
\]
The [Gegenbauer circle theorem](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md#mirrored-profiles-at-every-fixed-order)
gives simple unit-circle roots. For \(n\ge2\) and
\(\lambda_1\lambda_2\ge n^2/4\), we claim that
\[
          f''(\theta)/f'(\theta)<0 \quad\text{at each root in }(0,\pi),
       \qquad f(\theta)=e^{-in\theta/2}G(e^{i\theta}).
                                                               \tag{1}
\]
Let \(K=n/2\), \(A_j=K+\lambda_j\),
\(B=K+\lambda_1+\lambda_2\), and \(C=K^2-\lambda_1\lambda_2\).
The exact consecutive-coefficient ratio gives, with
\(E=z\partial_z\),
\[
 \left[E\prod_{j=1}^2(n+\lambda_j-E)
       -z(n-E)\prod_{j=1}^2(\lambda_j+E)\right]G=0 .
\]
In centered variables \(E=K-i\partial_\theta\).
The first term without \(z\) therefore acts as
\(K A_1 A_2+B\partial_\theta^2+
 i(C\partial_\theta+\partial_\theta^3)\).
The second term has the conjugate centered operator.
Subtracting the two terms and using \(z=e^{i\theta}\)
gives the real equation
\[
   \cos(\theta/2)(C f'+f''')
           =\sin(\theta/2)(K A_1 A_2 f+B f'').
                                                               \tag{2}
\]
At a root in \((0,\pi)\), put
\(r=f''/f'\), \(s=f'''/f'\), and
\(S=\sum_{\theta_j\ne\theta}\csc^2((\theta-\theta_j)/2)\),
where the sum is over the other circular roots.
The sine product gives
\[
                      s=\frac34(r^2-S)-\frac14 .
\]
Indeed the factor corresponding to the given root has
normalized expansion \(u-u^3/24+O(u^5)\), and the log
of the product of the remaining normalized factors is
\(r u/2-Su^2/8+O(u^3)\).
If \(r=0\), then \(s<0\), whereas (2) gives \(s=-C\ge0\).
Thus \(r\) never vanishes on a positive-half root branch
in the region \(\lambda_1\lambda_2\ge K^2\).

To determine its sign deform both shapes to \(\lambda_j+v\)
with \(v\to\infty\) and \(n\) fixed. Roots stay simple,
and the positive cut and the fixed central parity prevent
a positive-half root from leaving \((0,\pi)\).
Represent the Gamma moments by the independent products
\(U=X_1X_2\), \(V=Y_1Y_2\), with \(X_j,Y_j\) independent
Gamma variables of shape \(\lambda_j+v\), each divided
by its mean. The common mean divisions only normalize
the polynomial by a scalar. The centered polynomial
is, up to a positive scalar,
\[
 [\sin(\theta/2)]^n v^{-n/2}
             \mathbb E[(x(U+V)+i\sqrt v(V-U))^n],
                 \qquad x=\sqrt v\cot(\theta/2).
\]
The two individually centered and scaled product variables
\(\sqrt v(U-1),\sqrt v(V-1)\) converge in moments to
independent normals of variance \(2\).
For an individual normalized Gamma factor the
centered moment generating function is
\(e^{-t\sqrt v}(1-t\sqrt v/(\lambda_j+v))^{-(\lambda_j+v)}\),
which converges locally analytically to \(e^{t^2/2}\).
Expanding the two-factor products gives the claim,
including convergence of the fixed-degree mixed moments.
Thus the polynomial in \(x\) converges coefficientwise
to \(2^n\) times the monic Hermite polynomial with
generating function \(e^{xt-t^2/2}\).
Its roots are simple, and at each positive root its
second-to-first derivative ratio equals \(x>0\).
The roots and the derivative ratios of the polynomials
in \(x\) converge there.

For \(f=(\sin(\theta/2))^n p(x)\), the angular ratio at a root is
\[
      r=(n-1)\cot(\theta/2)
         -\frac{\sqrt v}{2}\bigl(1+\cot^2(\theta/2)\bigr)
                                      \frac{p''(x)}{p'(x)} .
\]
It is negative along every positive-half root for
sufficiently large \(v\), and hence throughout the
deformation. This proves (1). By reflection, the sign
of \(r\) reverses on \((\pi,2\pi)\).

## Application to the native Euler factors

For a shared mirrored pair of gap \(d\), the normalized
factor at degree \(m\) has Gegenbauer shapes
\((m\ell+j)/d\), \(1\le j\le d\), with \(\ell\ge1\).
Hence two shared unit-gap pairs have two shapes exceeding
\(m/2\), as does a single shared gap-two pair.
The same is true during the symmetric Gamma-shape stages
in the [exact polar factorization](MULTIVARIATE_MIRRORED_POLAR_EULER_FACTORS.md#factorization-and-the-attached-shapes).
Therefore (1) applies both to the ordinary lower polynomial
and to every intervening symmetric stage.

Let \(G_n\) be the advanced degree-\(n\) polynomial
with both Gegenbauer shapes at their upper values.
The upper polar half is
\[
                  H_N=(1-\beta_1 E)(1-\beta_2 E)G_n,
                   \qquad \beta_j=(n+\lambda_{N,j})^{-1}.
\]
At a root of \(G_n\), the
[two-residual sign formula](MULTIVARIATE_MIRRORED_POLAR_EULER_FACTORS.md#two-decreasing-factors-can-fail-the-sign)
gives the centered upper value divided by twice the
centered lower derivative as
\[
          (a_1\beta_2+a_2\beta_1)\sin(\theta/2)
                        -\beta_1\beta_2 r\cos(\theta/2),
                 \qquad a_j=1-\beta_j n/2 .
\]
This is positive on both halves by (1) and reflection.
At the central root the cosine term vanishes and the sine
term is positive. The lower-root sign test proves strict
interlacing of \(G_n,F_N\).

For the successive symmetric stages themselves,
one native update acts on the centered polynomial as
\(a+b\partial_\theta^2\) with \(a,b>0\).
At each old positive-half root the new centered value
has the sign opposite to the old derivative by (1).
Thus each positive-half root of the new stage lies
strictly to the right of the corresponding old root
and before the next old root, when there is one.
For the last interval use the sign
\((-1)^{\lfloor n/2\rfloor}\) of every centered
polynomial just below \(\pi\), since all stages are
positive at the cut and have the same parity and root count.
This establishes the claimed motion without assuming a
continuous phase comparison between the native stages.

Write \(l_i,a_i,u_i\) for the positive-half
root arguments of the ordinary lower polynomial, the
advanced lower polynomial, and the actual upper polynomial.
The estimates give
\[
                      l_i<a_i,\qquad u_i<a_i<u_{i+1}
\]
when \(u_{i+1}\) exists, with the central endpoint used
for the last interval when the upper degree is odd.
In particular \(l_i<u_{i+1}\).
For the ordinary consecutive pair the remaining
inequalities are \(u_i<l_i\). The sign of the
advanced-core comparison alone does not bound this
remaining symmetric displacement.
