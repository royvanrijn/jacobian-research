# All-degree coprimality for a shared unit gap and one private unit

For every integer \(a\ge1\), the ordinary moments of
\[
                          X^a(x+XY)
\]
are consecutively coprime in every degree, over every
characteristic-zero field. Equivalently, the exponent rows are
\((a,a+1),(0,1)\), with arbitrary equal rows allowed. The argument
uses a single estimate for every \(a\) and every terminal order.

## The exact comparison at a lower root

Use \(B_{j,\alpha}\) from the
[Robin resolvent comparison](MULTIVARIATE_ONE_PRIVATE_ROBIN_RESOLVENT.md).
For terminal \(N\ge2\) set \(m=N-1\),
\(\alpha=am+1\) and \(\rho=N+\alpha\). After scalar normalization
the two physical moments in the same variable \(x\) are
\(B_{m,\alpha}\) and \(B_{N,\alpha+a}\). At any nonreal zero
\(\xi\) of \(B_{m,\alpha}=B'_{N,\alpha}/N\), the Robin comparison
shows that
\[
                   y=\frac{B_{N,\alpha}(\xi)}{\xi^N}
\]
has imaginary part of the same strict sign as \(\Im \xi\).
The exact upper zero equation from that comparison is
\[
                       P_a(\xi)y+\xi V_a(\xi)=0,               \tag{1}
\]
where
\[
 \begin{split}
 P_0=1,\quad P_1=\rho,\qquad
 &P_{j+2}=(\rho+j+1)P_{j+1}-xP_j,\\
 V_0=V_1=0,\qquad
 &V_{j+2}=(\rho+j+1)V_{j+1}-xV_j+(\alpha)_j .
 \end{split}
\]
When \(a=1\), equation (1) would give \(y=0\). We now prove
that for every \(a\ge2\), the term \(\xi V_a(\xi)/P_a(\xi)\)
also has the same strict imaginary sign as \(\xi\).

The coefficients of \(B_{m,\alpha}\) in descending order are
\(c_k=(m)_{\underline k}(\alpha)_k\), where the underline denotes
a falling factorial. Their consecutive ratios satisfy
\[
             \frac{c_{k+1}}{c_k}=(m-k)(\alpha+k)\le m\alpha.
\]
This bounds the modulus of any root by \(R=m\alpha\).
For completeness, scaling a positive-coefficient polynomial
\(f(x)\) by its largest ascending coefficient ratio makes the
ascending coefficients of \(f(Rz)\) nondecreasing. In
\((1-z)f(Rz)\) the last coefficient's absolute value equals
the sum of all preceding coefficients, which are nonnegative.
For \(|z|>1\) the last term strictly dominates the sum of the others.
In our case \(R\) is an upper bound on that largest ratio, so the
same proof applies.

## Trailing determinant phases

Write \(\Delta_i(x)\), \(0\le i<a\), for the trailing tridiagonal
determinant with diagonal
\(\rho+i,\rho+i+1,\ldots,\rho+a-1\) and off-diagonal elements
\(-\sqrt{x}\), and put \(\Delta_a=1\). These determinants are
polynomials in \(x\), and \(\Delta_0=P_a\). The recurrence for
\(V_a\) gives
\[
           V_a=\sum_{k=0}^{a-2}(\alpha)_k \Delta_{k+2}.
                                                               \tag{2}
\]
Indeed the forcing term \((\alpha)_k\) enters at index \(k+2\)
and propagates to index \(a\) by the homogeneous determinant
recurrence from the successive diagonal entries. The first
two values of that propagation are \(1\) and \(\rho+k+2\).

Set \(b_i=\Delta_{i+1}/\Delta_i\). For each trailing determinant,
congruence by the inverse square root of the diagonal gives
a positive scalar times \(\det(I-\sqrt{x}J)\), with \(J\) a
real symmetric tridiagonal matrix of zero diagonal and off-diagonal
entries \(1/\sqrt{(\rho+j)(\rho+j+1)}\). Its operator norm is
at most \(2/\rho\) by the row sums. Its eigenvalues are simple,
are symmetric about zero, and strictly interlace those of its
next trailing block. Simplicity and strictness follow from the
nonzero off-diagonals and the determinant recurrence. Passing
from its positive eigenvalues to their inverse squares shows
that each pole of \(b_i\) is positive and at least \(\rho^2/4\),
and is followed by a zero before the next pole. The last pole
may be paired with infinity when there is one fewer zero.

For \(|x|<\rho^2/4\), choose the logarithm of
\(b_i(x)/b_i(0)\) vanishing at zero. Pairing each pole \(p\)
with its following zero \(z\) gives the exact identity
\[
 \log \frac{b_i(x)}{b_i(0)}
                       =\sum_{(p,z)}\int_p^z \frac{x}{t(t-x)}\,dt .
                                                               \tag{3}
\]
For \(i=a-1\) the ratio is constant and the integral is empty.
For the other ratios, differentiating at zero, or using the
exact recurrence \(b_i(x)^{-1}=\rho+i-xb_{i+1}(x)\), gives
\[
  \sum_{(p,z)}\int_p^z \frac{dt}{t^2}
                          =\frac1{(\rho+i)(\rho+i+1)} .
                                                               \tag{4}
\]
If \(\Im x>0\), all the unwrapped phases in (3) are nonnegative.
These formulas use only the determinants on the actual, common
degree-linked pencil.

## The half-plane sign on the root disk

Since \(\rho=(a+1)m+2\) and \(\alpha=am+1\), put
\[
 K=R\left(\frac1{\rho}-\frac1{\rho+a-1}\right)
                    <\frac{a(a-1)}{(a+1)^2}<1,\qquad
 u=\frac{4R}{\rho^2}<\frac{4a}{(a+1)^2}<1 .
\]
The bounds follow from
\[
       a\rho^2-m(am+1)(a+1)^2
                  =m(a+1)(3a-1)+4a>0 .
\]
For \(x=|x|e^{i\theta}\) with \(0<\theta<\pi\), \(|x|\le R\),
consider the phase of a summand
\(x\Delta_{k+2}(x)/P_a(x)=x\prod_{i=0}^{k+1}b_i(x)\) in (2).
Its phase lift is \(\theta\) plus the sum of the logarithmic
imaginary parts in (3). If \(\Re x\le0\), then
\(|t-x|^2\ge t^2\), so (3)–(4) bound the increment by
\(K\sin\theta<\pi-\theta\).

For \(\Re x\ge0\), write \(w=|x|/t\) on the intervals in (3).
Then \(0\le w\le u\). The Cauchy–Schwarz inequality
\[
          2w\cos\theta+(1-w^2)\sin\theta\le1+w^2
\]
gives
\[
       \frac{\sin\theta}{1-2w\cos\theta+w^2}
                                             \le\frac1{1-w^2}.
\]
Using this in the imaginary part of (3) bounds the increment
by \(K/(1-u^2)\). For \(a\ge3\) this is less than one, since
\[
    (3a+1)(a+1)^2-16a^2=(a-1)(3a^2-6a-1)>0 .
\]
For \(a=2\) it is less than \(18/17<\pi/2\). Thus the phase
lift of every summand is strictly between zero and \(\pi\) in
the upper half of the root disk. All the weights \((\alpha)_k\)
in (2) are positive. Consequently \(xV_a(x)/P_a(x)\) has
strictly positive imaginary part there. Conjugation gives the
corresponding negative sign in the lower half.

At a nonreal lower root, equation (1) would therefore add two
numbers with the same strict imaginary sign and give zero.
No real common root exists by strict positivity of the even
moment on the real axis. At the two coefficient endpoints the
phase is a single monomial. The resulting rational Bézout
identities transfer coprimality to every characteristic-zero
field.
