# Sharp curvature triangle at stationary Gaussian contacts

Fix consecutive orders \(n,N=n+1\) and a finite ordinary factorial
core. The [signed jet tests](MULTIVARIATE_POSITIVE_DIFFUSE_EXCEPTION_SIGNED_JETS.md)
use an angular expansion in the strength of a sufficiently diffuse
tail and the next log-coefficient terms. At a positive stationary
contact, the former feasibility test was
\(0\le u\le W_*,\ 0\le v\le W_*u\), for the normalized quartic
and squared cubic responses \(u\) and \(v\). The ordinary row
integrality and the limiting signed bias \(\lambda_*\) give the
smaller sharp triangle
\[
  0\le u\le M_*,
  \qquad \max(0,2M_*u-M_*^2)\le v\le M_*u,
  \qquad M_*=W_*-|\lambda_*|.                                      \tag{1}
\]
In particular, the triangle is a single point at saturated bias.

## Exact finite-profile constraints

For each unequal shared ordinary row set
\[
  m_i=(a_i+b_i)/2,\qquad d_i=b_i-a_i,\qquad
  x_i=d_i/m_i,\qquad w_i=d_i^2/m_i .
\]
Use \(W=\sum_i w_i\), \(\lambda=\sum_i x_i\),
\(\eta=\max_i |x_i|\), and \(Z=\sum_i x_i^2>0\).
The signed jet tests use
\[
 K=\sum_i w_i x_i,\qquad B_4=\sum_i w_i x_i^2,\qquad
 u=B_4/Z,\qquad k=K/\sqrt Z,\qquad v=k^2 .
\]
Put \(y_i=x_i/\sqrt Z\). The integer inequality \(|d_i|\ge1\)
allows the nonnegative trimmed weights
\[
   \omega_i=w_i-|x_i|,\qquad
   A=\sum_i\omega_i=W-\sum_i|x_i|\le M=W-|\lambda|.
\]
Define \(u_0=\sum_i\omega_i y_i^2\) and
\(k_0=\sum_i\omega_i y_i\). Then
\[
 \begin{aligned}
   0\le u-u_0&=\sum_i |x_i| y_i^2\le\eta,\\
   |k-k_0|&=\left|\sum_i|x_i|x_i/\sqrt Z\right|\le\sqrt Z
                   \le\sqrt{W\eta}.
 \end{aligned}                                                     \tag{2}
\]
The last inequality follows from
\(\sum_i x_i^2\le\eta\sum_i|x_i|\le\eta W\).

The pair \((u_0,k_0^2)\) satisfies (1) with \(M_*\) replaced by
either \(A\) or \(M\). For \(M>0\), put the extra weight \(M-A\)
at an amplitude of zero. Write \(\kappa=k_0/M\) for the mean
and \(U=u_0/M\) for the second moment under the resulting
probability weights. Every amplitude is in \([-s,t]\)
where \(s,t\ge0\) and \(s^2+t^2\le\sum_i y_i^2=1\).
The inequality for the two extreme amplitudes implies
\(y_i^2\le(t-s)y_i+st\). Hence
\[
   2U-\kappa^2
      \le s^2+t^2-(\kappa-t+s)^2\le1 .
\]
The other bound is \(\kappa^2\le U\) by Cauchy–Schwarz.
The same argument works at mass \(A>0\). If either chosen
mass is zero the two trimmed moments vanish. Thus, for
bounded \(W\), the untrimmed \((u,v)\) lies within a uniformly
bounded \(O_W(\sqrt\eta)\) distance of the triangle at \(M\).
For instance \(|k|,|k_0|\le W\), so
\(|k^2-k_0^2|\le2W\sqrt{W\eta}\).

## Sharpness with arbitrary fixed integer signed gap

Fix \(W_*>0\), \(|\lambda_*|\le W_*\), and any integer \(h\).
Every point of (1) is the limit of \((B_4/Z,K^2/Z)\) along
finite ordinary shared tails for which
\[
       \eta\longrightarrow0,\qquad
       W\longrightarrow W_*,\qquad
       \lambda\longrightarrow\lambda_*,\qquad
       \sum_i d_i=h .
\]
One can prescribe either sign for the limiting \(K/\sqrt Z\).

For \(M_*>0\) write \(U=u/M_*\) and choose
\(\kappa=\pm\sqrt v/M_*\). Take the roots \(y_1,y_2\) of
\[
                    Y^2-\kappa Y+\kappa^2-U=0.
\]
They are real and the interval between them contains
\(\kappa\), because the left side at \(\kappa\) is
\(\kappa^2-U\le0\). Choose probabilities \(p_1,p_2\) on the
two roots giving mean \(\kappa\). Their second moment is \(U\)
by the same quadratic equation, and
\(y_1^2+y_2^2=2U-\kappa^2\le1\). If the two roots are equal,
use the single point at zero. Omit any zero probability.
The fixed target row weights are \(w_i^*=M_*p_i>0\).

For a positive parameter \(\varepsilon\to0\) and nonzero
\(y_i\), choose an integer \(d_i\) of sign \(\operatorname{sgn}(y_i)\)
whose magnitude is nearest to \(w_i^*/(\varepsilon |y_i|)\).
Choose \(m_i\) nearest to \(d_i^2/w_i^*\) in
\(\mathbb Z+d_i/2\). Then \(m_i\pm d_i/2\) are positive
integers for small \(\varepsilon\),
\(w_i\to w_i^*\), and \(x_i/\varepsilon\to y_i\).
If \(y_i=0\), choose \(d_i=\lceil\varepsilon^{-2}\rceil\)
and the same midpoint rounding to obtain \(x_i=o(\varepsilon)\).
For the remaining variance slot
\(c^2=1-\sum_i y_i^2>0\), add a single increasing unit row
with midpoint nearest \(1/(\varepsilon c)\) in
\(\mathbb Z+1/2\). It contributes
\(Z/\varepsilon^2\to c^2\) and has vanishing weight.
In using this slot, the sum of squared amplitudes includes
only the retained positive-probability nodes.

If \(M_*=0\) omit all prescribed positive-weight nodes and
use the full remaining variance slot \(c=1\).
In every case, add
\(\lfloor|\lambda_*|(H+1/2)\rfloor\) unit-gap rows with
midpoint \(H+1/2\) and the sign of \(\lambda_*\), where
\(H=\lceil\varepsilon^{-4}\rceil\). Omit these rows when
\(\lambda_*=0\). Their combined curvature tends to
\(|\lambda_*|\), their bias tends to \(\lambda_*\), and
their contribution to \(Z\) is \(O(\varepsilon^4)\).
If the signed gap is now \(h-D\), add \(|D|\) unit-gap rows
of the sign of \(D\) at an integer smaller exponent exceeding
\((1+|D|)\varepsilon^{-4}\). Their total curvature and
absolute bias are at most \(\varepsilon^4\).
The full profile has \(\sqrt Z/\varepsilon\to1\) and the
required limiting strengths. Its limiting \(k\) and \(u\)
are \(\sum_i w_i^*y_i=M_*\kappa\) and
\(\sum_i w_i^*y_i^2=M_*U\). This proves sharpness in the
moment coordinates. The construction prescribes no exact
common root.

## Uniform stationary-contact exclusion

Use the notation of the
[stationary signed jet test](MULTIVARIATE_POSITIVE_DIFFUSE_EXCEPTION_SIGNED_JETS.md#stationary-contacts-with-a-common-phase-sign).
Suppose the two contacting roots are simple at a positive
physical reference contact \((W_*,\lambda_*)\). The phase
difference is \(g(W)=a(W-W_*)^{2r}+O((W-W_*)^{2r+1})\)
with \(a\ne0\), \(r\ge1\). Assume its first signed cubic
response vanishes through order \(r-1\), and let \(b\)
denote its coefficient of order \(r\). Denote the limiting
clock, quartic and second cubic responses by \(B,C,D\)
respectively, as in that test. A sufficient condition to
exclude the contact for all sufficiently diffuse ordinary
tails within a neighborhood is
\[
 \operatorname{sgn}(a)
       \left[B+C u+\left(D-\frac{b^2}{4a}\right)v\right]>0         \tag{3}
\]
at the three points
\[
                 (0,0),\qquad (M_*/2,0),\qquad (M_*,M_*^2).
\]
The middle point comes from two equally weighted opposite
amplitudes. The last point comes from one amplitude of modulus
one. This criterion uses no bound on the tail row count, gaps
or heights. All hypotheses concern the two fixed reference
polynomials and this one contact.

Indeed the imaginary part of the signed jet expansion first
forces \(Y=(W-W_*)^r/\sqrt Z\) to stay bounded along a putative
contact sequence. The bounds \(K^2\le W^2 Z\) and
\(B_4\le W Z\), together with the vanishing to order \(r\)
of the cubic response, give
\(|W-W_*|^{2r}=O(\sqrt Z\,|W-W_*|^r+Z)\).
After division by \(Z\), the contact equation is then
\[
    0=aY^2+b kY+B+C u+D k^2+o(1).
\]
Completing the square gives the sum of a term with sign
\(\operatorname{sgn}(a)\) or zero and the affine expression
in (3). The exact bounds (2) and
\(M\to M_*\) put the limiting pair in the triangle (1).
Strict positivity on the three vertices therefore excludes
the sequence. This argument is uniform and effective using
the signed jet Taylor bounds, (2) and the strict vertex
margins. In particular at saturated bias it suffices to test
\(\operatorname{sgn}(a)B>0\).

The signs in (3) can be decided exactly. All the response
coefficients are algebraic at the exceptional clock, and
\(M_*=-2nN\log(q_*v_*^{-\epsilon})\) with
\(q_*=e^{-W_*/(2nN)}\), \(v_*=e^{-\lambda_* /(2nN)}\) both
positive algebraic and \(\epsilon=\operatorname{sgn}(\lambda_*)\).
If \(M_*=0\) use the single vertex directly. Otherwise
\(M_*\) is transcendental by the same Hermite–Lindemann
argument as in the signed jet test, and each vertex value is
a polynomial in \(M_*\) with algebraic coefficients. An
identically zero polynomial fails strict positivity. Every
other sign is decided by convergent interval evaluation.
No such sign assumption is made about transverse interior
contacts, where the stationary test does not apply.
