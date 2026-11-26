# Universal linear Cayley compression misses an actual factorial pair

Fix \(N\ge6\), \(n=N-1\), and use the Cayley transform
\[
     {\cal C}_j P(x)=(x+i)^jP\bigl((x-i)/(x+i)\bigr).
\]
A real palindromic polynomial of degree \(j\) transforms
to a real polynomial of parity \(j\). For positive
coefficients the degree is preserved, since the leading
coefficient is \(P(1)>0\). The unit circle away from the
positive cut maps to the real line, preserving the
interlacing order up to reversal.

**Theorem.** Suppose \(T\) is a real linear map from the
degree-\(N\) parity space to real polynomials of degree at
most \(n\). Assume that for every positive palindromic
degree-\(N\) polynomial \(P\) with simple unit-circle
roots and strictly log-convex binomial-normalized
coefficients, \(T({\cal C}_N P)\) has degree \(n\) and
strictly interlaces \({\cal C}_N P\).
For each integer \(b\ge1\), let
\[
          U=X^{b+1}Y^b,\qquad V=X^bY^{b+1},
\]
with ordinary factorial moments, and write \(P_j\)
for the moment polynomial of \(U+zV\) at order \(j\).
Then \(T({\cal C}_N P_N)\) cannot be a nonzero scalar
multiple of \({\cal C}_n P_n\).
The conclusion permits \(T\) to be chosen separately
for each degree and profile.

This addresses exact input-output interlacing by a
universal linear compression. General zero-set
preservers are classified by
[Borcea and Brändén](https://doi.org/10.4007/annals.2009.170.465).
The repeated-root tests below use the stronger
requirement that the output interlace its own input.
The actual moment pair in the theorem itself strictly
interlaces by the
[single mirrored unit-gap theorem](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md#consecutive-interlacing-for-a-single-mirrored-unit-gap).

## Compression on the parity spaces

Write \(y=x^2\), \(N=2p+\epsilon\) with
\(\epsilon\in\{0,1\}\), and take the repeated inputs
\(x^\epsilon(y-u)^p\) for small real \(u>0\).
Their inverse Cayley images, up to positive scalars,
are
\[
                (1+z)^\epsilon\big((1+z)^2-tz\big)^p,
                         \qquad t=\frac{4u}{1+u}.
\]
Their binomial-normalized coefficients at small \(t\)
are
\[
                  1-\frac{p\,k(N-k)}{N(N-1)}\,t+O(t^2).
\]
The second differences of their logarithms are
\(2pt/[N(N-1)]+O(t^2)>0\) at every internal index
for sufficiently small \(t>0\).
Splitting the repeated positive \(u\)'s into nearby
distinct positive numbers keeps all coefficients
positive and preserves these strict logarithmic
inequalities. It gives simple unit-circle roots,
including the odd structural root. The repeated
inputs are therefore in the closure of the
assumed universal interlacing class.

For a sequence tending to a repeated input,
the interlacing output roots remain in a fixed
bounded interval. If the leading output coefficient
tends to zero, all its coefficients tend to zero
by the elementary symmetric root formulas.
Otherwise the limiting output has the repeated
factors dictated by weak interlacing.
In either case it vanishes to order at least
\(p-1\) at both \(\sqrt u\) and \(-\sqrt u\).
Its projection to the target parity retains
those factors. Use this projection, denoted
\(T_{\rm par}\), throughout the classification.

If \(N=2p\), the projected output on the repeated
input is \(c(u)x(y-u)^{p-1}\).
The coefficient of \(x^{2p-1}\) shows \(c(u)\)
is polynomial in \(u\) of degree at most \(p\).
The coefficient of \(x\) is \(c(u)(-u)^{p-1}\)
and also has degree at most \(p\), so \(c(u)=A+Bu\).
The powers \((y-u)^p\), as \(u\) ranges over
any open interval, span the degree-\(\le p\)
space. Hence for every \(f\) in that space,
\[
       T_{\rm par}(f(y))
          =x\left[\frac{A+By}{p}f'(y)-B f(y)\right].
                                                               \tag{1}
\]
One can check (1) directly on the repeated inputs.

If \(N=2p+1\), the even projection on the
repeated inputs is \((y-u)^{p-1}(a(u)y+b(u))\).
The two highest coefficients show that \(a,b\)
are polynomials in \(u\). The constant coefficient
then gives \(\deg b\le1\), and the coefficient
of \(y\) gives \(\deg a\le1\).
As \(u\) decreases to zero, weak interlacing
forces the limiting output on \(x^{2p+1}\)
to be a scalar multiple of \(x^{2p}\),
possibly zero. Thus \(b(0)=0\), so
\(a(u)=A+Bu\), \(b(u)=Cu\).
Again the repeated inputs span the source space.
Substitution yields
\[
      T_{\rm par}(x f(y))
        =\frac{(A+C)y+By^2}{p}f'(y)-(By+C)f(y).
                                                               \tag{2}
\]
The comparison with an actual target of the required
parity may therefore use (1) or (2), even when
\(T\) has an additional component of the other parity.

## Exact obstruction on the factorial coefficients

For the stated unit-gap pair, the Cayley moments
at order \(j\) are scalar multiples of the monic
polynomial \(F_{j,\lambda_j}\), where
\(\lambda_j=jb+1\) and
\[
 F_{j,\lambda}(x)=\sum_{h=0}^{\lfloor j/2\rfloor}
          c_{j,h}(\lambda)x^{j-2h},\qquad
 c_{j,h}(\lambda)=
       \frac{(-1)^h j!}{h!(j-2h)!\,4^h(\lambda+1/2)_h}.        \tag{3}
\]
To verify (3), the binomial-normalized factors at
fixed \(\lambda\) are \((\lambda)_k(\lambda)_{j-k}\).
The Cayley generating function with division by
\(j!\) is
\[
  \left(1-2xv+(x^2+1)v^2\right)^{-\lambda}
    =\sum_{h\ge0}\frac{(-1)^h(\lambda)_h}{h!}
                      v^{2h}(1-xv)^{-2\lambda-2h}.
\]
Divide the degree-\(j\) coefficient by its leading
coefficient and use
\((2\lambda)_{2h}=4^h(\lambda)_h(\lambda+1/2)_h\).

Put \(\lambda=\lambda_N\), \(a=\lambda+1/2\),
and \(c=a-b>0\), so \(\lambda_n=\lambda-b\).
The lower factorial coefficient ratio relative
to the fixed-upper-parameter polynomial is
\[
       \frac{c_{n,h}(\lambda-b)}{c_{n,h}(\lambda)}
                         =\frac{(a)_h}{(c)_h}.               \tag{4}
\]
For even \(N=2p\), substituting (3) in (1)
shows that the coefficient at \(x^{n-2h}\),
after division by \(c_{n,h}(\lambda)\), is
\[
                        A+\frac{B(n-2h)}{2(a+h)}.
\]
In detail, the \(A\) term is
\(A(1-h/p)c_{N,h}(\lambda)=A c_{n,h}(\lambda)\).
The \(B\) term is
\(-B(h+1)c_{N,h+1}(\lambda)/p\), giving the
displayed ratio by (3). An exact match with
(4), up to a scalar, would force
\[
                        R_h=\frac{(a+1)_h}{(c)_h}
\]
to be affine in \(h\). But \(0,1,2\) are
available indices and
\[
                       \Delta^2R_0=\frac{b(b+1)}{c(c+1)}>0.
\]

For odd \(N=2p+1\), substitution in (2)
instead gives the divided coefficients
\[
  \frac{N}{N-2h}\left(A-\frac{(A+C)h}{p}\right)
                    +\frac{B N(N-2h-1)}{4p(a+h)}.
\]
Here the \(B\) term is again obtained from
\(-B(h+1)c_{N,h+1}(\lambda)/p\), and the
other terms use
\(c_{n,h}(\lambda)=(N-2h)c_{N,h}(\lambda)/N\).
Matching (4) would force \((N-2h)R_h\)
to be quadratic in \(h\). The indices
\(0,1,2,3\) are available. Set \(d=b+1\);
the rising ratios give
\[
 \Delta^2 R_1=\frac{d(d-1)(a+1)}{c(c+1)(c+2)},\qquad
 \Delta^3 R_0=\frac{d(d-1)(d-2)}{c(c+1)(c+2)}.
\]
The discrete product formula yields
\[
 \Delta^3\big((N-2h)R_h\big)_{h=0}
      =\frac{d(d-1)[N(d-2)-6(a+1)]}{c(c+1)(c+2)}<0,
\]
because \(a=Nb+3/2\), so the bracket is
\(-5Nb-N-15\). This excludes exact matching
in both parities, for every \(b\ge1\).
