# The initial Strong Factorial window on coalescing ternary rays

Fix distinct ordinary exponent columns \(b_1,\ldots,b_s\in\mathbb
Z_{\ge0}^r\) with \(s\le3\), and a strictly positive integer vector
\(a=(a_\rho)\). For a positive integer \(T\) give term \(j\) the
exponents \(Ta+b_j\).

**Theorem.** For all sufficiently large \(T\), the first \(s\)
positive factorial moments have no simultaneous nonzero projective
coefficient zero. This holds over every characteristic-zero field.
The same holds for every subset of the columns with its own support
size. The translation preserves the exponent differences and hence
can yield fully shared homogeneous triples with arbitrarily small
relative differences.

For any support size, if the total degrees of the \(b_j\) are
pairwise distinct, the same conclusion holds for any prescribed
\(s\) distinct positive orders, with a threshold depending on those
orders. For \(s=3\) with repeated degrees, the theorem concerns
the initial window at orders \(1,2,3\).

## Exact tilted Gamma normalization

Write \(d_j=\sum_\rho b_{\rho j}\) and set \(t=1/T\). The ordinary
Gamma integral for \(\Lambda\) tilts coordinate \(\rho\) at order
\(m\) to a unit-rate Gamma law of shape \(mTa_\rho+1\).
For independent variables \(X_\rho\) of these shapes, put
\(Y_\rho=X_\rho/(mTa_\rho)\). Removing the common nonzero factorial
\(\prod_\rho (mTa_\rho)!\) and changing coefficient coordinates by
\(C_j=c_jT^{d_j}\prod_\rho a_\rho^{b_{\rho j}}\), the normalized
moment is exactly
\[
      M_m(C,t)=\mathbb E\!\left[
        \left(\sum_j C_j m^{d_j}Y^{b_j}\right)^m\right].
 \tag{1}
\]
The coefficient change is independent of \(m\). We use the Gamma
moment and cumulant formulas, which follow by differentiating the
identity \(\mathbb E e^{zX_\rho}=(1-z)^{-mTa_\rho-1}\).

Put \(\delta_\rho=Y_\rho-1\) and \(\epsilon=t/m\). Its first
cumulant is \(\epsilon/a_\rho\), and for \(k\ge2\) its \(k\)-th
cumulant is
\[
      (k-1)!\left(1+\epsilon/a_\rho\right)
                 (\epsilon/a_\rho)^{k-1}.
 \tag{2}
\]
All mixed cumulants of independent coordinates vanish. The
moment-cumulant partition formula shows that each joint moment
of degree \(h\) in \(\delta\) is \(O(t^{\lceil h/2\rceil})\) for
fixed \(m\). It also shows convergence of all joint moments of
\(\delta/\sqrt t\) to those of independent centered real Gaussians
with variances \(1/(m a_\rho)\). In the proof all moment orders
being considered are fixed, and all coefficients are bounded
after normalization, so the errors in finite polynomial
expansions are uniform.

If the asserted threshold fails, choose \(T\to\infty\) and
projective common zeros with \(\max_j|C_j|=1\), and take a convergent
subsequence. Coefficient zeros over an arbitrary
characteristic-zero field imply zeros over the complex numbers:
the coefficient system has rational coefficients and a nonempty
affine chart after imposing one coefficient equal to one.
The Nullstellensatz gives a point in that chart over
\(\overline{\mathbb Q}\) if it has a point over any extension.

## Distinct total degrees

For fixed \(m\), (1) tends to
\((\sum_j C_j m^{d_j})^m\). If all degrees are distinct, any \(s\)
distinct positive orders force the limiting \(C\) to vanish:
a nonzero polynomial using only \(s\) distinct integer degrees
has at most \(s-1\) distinct positive real roots by Descartes'
rule of signs. For complex coefficients apply the real rule to
a nonzero real or imaginary part. This contradicts projective
normalization. This argument works for arbitrary \(s\) in the
distinct-degree case.

## Three equal total degrees

Remove the factor \(m^{m d_1}\) in (1), and write
\[
 S=\sum_j C_j,\qquad
 L_\rho=\sum_j C_j b_{\rho j},\qquad
 H_{\rho\sigma}=\sum_j C_j
       (b_{\rho j}b_{\sigma j}-\mathbf1_{\rho=\sigma}b_{\rho j}),
 \qquad
 V=\sum_\rho L_\rho^2/a_\rho .
\]
Here \(S,L,H\) are the constant, gradient and Hessian of
\(\sum_j C_jY^{b_j}\) at \(Y=\mathbf1\). From \(M_1=0\) and the
moment bound after (2) we get \(S=O(t)\). With this estimate,
\[
             M_2/2^{2d_1}=S^2+\tfrac t2 V+O(t^2)
\]
gives \(V=O(t)\). Thus at a projective limit \(S=V=0\).
For two columns, these limiting conditions are already
impossible: with their coefficients summing to zero, \(L\)
is a nonzero complex scalar times a nonzero real vector,
and has nonzero bilinear variance \(V\). The single-term case
is immediate.

### Independent offsets and the cubic obstruction

For noncollinear triples, \(L\ne0\) at a nonzero limiting \(C\)
with \(S=0\). Put \(D_a=\operatorname{diag}(1/a_\rho)\). The
only order-\(t^2\) terms of \(M_3/3^{3d_1}\) after using
\(S,V=O(t)\) are the cube of the linear term and three times
the square of the linear term times the quadratic term.
Indeed the remaining terms either contain \(S\) or have
degree at least five in \(\delta\), and the term
\(3S\,\mathbb E(L\cdot\delta)^2\) is \(O(t^3)\) since
\(\mathbb E(L\cdot\delta)^2=tV/3+O(t^2)\).

At \(m=3\), the cube of the linear form contributes
\[
   \epsilon^2\left(2\sum_\rho L_\rho^3/a_\rho^2
                                  +3V\sum_\rho L_\rho/a_\rho\right)
       +O(\epsilon^3),
\]
by (2). The term from the quadratic form contributes
\[
  \tfrac32\epsilon^2
     \left(V\,\operatorname{tr}(H D_a)+
          2(D_a L)^{\mathsf T}H(D_a L)\right)
           +O(\epsilon^3)
\]
by the three pairings in the fourth joint moment. Consequently
vanishing of the third moment requires the limiting identity
\[
            3\sum_j C_j\big(b_j\cdot D_a L\big)^2
                         -\sum_\rho L_\rho^3/a_\rho^2=0.
 \tag{3}
\]
The subtraction includes the diagonal correction in \(H\).

To exclude (3), put \(y_{\rho j}=b_{\rho j}/\sqrt{a_\rho}\) and
\(\widetilde L_\rho=L_\rho/\sqrt{a_\rho}\). The \(y_j\) form a
real nondegenerate triangle of circumradius \(R\). Since
\(\widetilde L\) is nonzero and isotropic and lies in the complex
span of its edge plane, write \(\widetilde L=\ell(e+if)\)
with orthonormal real vectors \(e,f\) in that plane and
\(\ell>0\). Write \(z_j=(y_j-y_1)\cdot(e+if)\). We have
\(\sum_j C_j z_j=0\) and \(\sum_j C_j\overline z_j=2\ell\).
Shifting the origin to \(y_1\) does not alter the sum of squares
in (3), since \(\sum_j C_j=V=0\). Solving for \(C_2,C_3\)
with \(z_1=0,z_2=u,z_3=v\) gives
\[
       \left|\sum_j C_j z_j^2\right|
       =\frac{2\ell|u v(u-v)|}{|v\overline u-u\overline v|}
       =2\ell R ,
\]
where the denominator is four times the triangle area. The
first term in (3) therefore has modulus \(6\ell^3 R\).

Let \(w_\rho=e_\rho+i f_\rho\). If \(w_\rho\ne0\), at least
one edge has nonzero projection to that coordinate. The
integrality of \(b\) gives that edge a coordinate difference
of magnitude at least \(1/\sqrt{a_\rho}\). All projective edge
lengths in the coordinate \(\rho\) are bounded above by
\(2R |w_\rho|\), since the unrestricted edge length is at most
\(2R\) and the coordinate projection on the plane has norm
\(|w_\rho|\). Hence
\[
  \left|\sum_\rho L_\rho^3/a_\rho^2\right|
  \le \ell^3\sum_\rho |w_\rho|^3/\sqrt{a_\rho}
  \le 2R\ell^3\sum_\rho |w_\rho|^4
  \le 4R\ell^3 .
\]
The last inequality uses \(|w_\rho|^2\le1\) and
\(\sum_\rho|w_\rho|^2=2\). Thus (3) is impossible in every
positive integral direction \(a\).

### Collinear offsets and a quadratic Gaussian limit

Write \(b_j=b_0+\lambda_j u\) for three distinct real numbers
\(\lambda_j\) and a nonzero real direction \(u\). Set
\(J_k=\sum_j C_j\lambda_j^k\), so \(J_0=S\) and
\(L=b_0 S+u J_1\). Since \(S,V=O(t)\) and
\(\sigma^2=\sum_\rho u_\rho^2/a_\rho>0\), the expansion
\(V=\sigma^2 J_1^2+O(t)\) gives \(J_1=O(\sqrt t)\).
At the coefficient limit \(J_2\ne0\), by the Vandermonde
matrix on the three distinct \(\lambda_j\).

Pass to limits of \(S/t\) and \(J_1/\sqrt t\). The Hessian
tends to \(J_2 uu^{\mathsf T}\), since all its other parts
contain \(S,J_1\), or \(L\). The rescaled polynomial
\(t^{-1}\sum_j C_jY^{b_j}\) now tends in moments to a common
quadratic \(A+B Z+C Z^2\) with \(C\ne0\), evaluated with
\(Z\) centered real Gaussian of variance \(1/m\) at order
\(m\). This follows by applying (2) to the Taylor expansion:
terms of degree \(h\ge3\) are \(t^{h/2-1}\) times bounded
polynomials in \(\delta/\sqrt t\). After dividing out
\(C^m\), the limiting polynomial is \(Z^2+\beta Z+\gamma\).
Its first moment forces \(\gamma=-1\), and at order two the
moment is \(3/4+\beta^2/2\), giving \(\beta^2=-3/2\).
At order three the Gaussian has variance \(1/3\), and the
moment is \(-4/9\), since the \(\beta^2\) contribution is
\(3\beta^2(\mathbb E Z^4-\mathbb E Z^2)=0\).
Thus the three moments cannot vanish.

## Two equal degrees and one distinct degree

Let \(d_1=d_2=d\), \(d_3=d+k\) with a nonzero integer \(k\).
Put \(S=C_1+C_2,\ q=C_3\), and remove \(m^{md}\) from each
normalized moment. At order one,
\(S+q=O(t)\). Order two yields
\((S+2^k q)^2=O(t)\) by (2), so \(q=O(\sqrt t)\) and
\(S=O(\sqrt t)\). The coefficient limit has
\(C_1=-C_2\ne0\), and we can take \(A=\lim q/\sqrt t\).

After dividing the polynomial inside the moment by \(\sqrt t\),
its limit in moments at order \(m\) has mean
\(A(m^k-1)\) and a centered Gaussian linear term with
bilinear variance \(W/m\), where
\[
                    W=C_1^2\sum_\rho
                        (b_{\rho1}-b_{\rho2})^2/a_\rho\ne0 .
\]
Indeed the constant comes from \(S+q=O(t)\), the Taylor
linear term tends to \(C_1(b_1-b_2)\cdot(\delta/\sqrt t)\),
and the higher terms tend to zero in moments. Vanishing at
order two requires \(A^2(2^k-1)^2=-W/2\), so \(A\ne0\).
The third moment requires \(A^2(3^k-1)^2=-W\), since
\(3^k-1\ne0\). It would follow that
\(((3^k-1)/(2^k-1))^2=2\), contrary to rationality of the
ratio. This completes the initial ternary theorem.

## Comparison

The support-size formulation is [Edo and van den Essen,
*The Strong Factorial Conjecture*, Conjecture 2.8](https://arxiv.org/abs/1304.3956).
The [saddle-envelope theorem](SFC_DILATED_SUPPORT_SADDLE_ENVELOPE.md)
treats large dilation with separated normalized offsets and
distinct degree slopes. Here the normalized columns coalesce
along a common ray, and the first nonvanishing order in the
Gamma expansion is crucial at an isotropic coefficient limit.
The [existing-row peeling theorem](SFC_EXISTING_ROW_PEELING_CERTIFICATE.md)
varies separate row baselines and yields certificates on finite
collections of order sets. The present initial theorem works
on every prescribed strictly positive common-height direction.
