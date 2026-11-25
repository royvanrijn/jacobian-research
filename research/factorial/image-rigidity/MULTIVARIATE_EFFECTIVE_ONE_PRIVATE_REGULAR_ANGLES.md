# Effective angular separation on the regular one-private bands

Orient the private rows of an ordinary two-monomial profile as
\((0,e_i)\). On unequal shared rows let \(d_i=b_i-a_i\). Put
\[
 P=\sum_{\mathrm{private}}e_i,\qquad p=\#\{\text{private rows}\},
 \qquad \sigma=\sum_{\mathrm{shared}}d_i^2/a_i,\qquad
 S=P+\sum_{\mathrm{shared}}\frac{d_i^2}{\max(a_i,b_i)}.
\]

**Theorem.** Given rational caps \(K>0\) and \(\kappa>0\), there is
a computable cutoff above which every such profile with \(S\le K\)
and either \(P\ge2\) or \(P=1,\ \sigma\ge1+\kappa\) has disjoint
consecutive root arguments. The exclusion holds under independent
positive rescalings of the polynomials.

We use the compact families and certified one-term and two-term
expansions of the
[effective regular one-private comparison](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_REGULAR_COMPLEMENT.md).
Those expansions are local in the field and independent of any
relation between its real parts at the two orders. Equal arguments give
fields \(\xi_n,\xi_N\), \(N=n+1\), of the same height in the corresponding
centered logarithms. By conjugation we may take the height in \([0,\pi]\).
The one-term witnesses already exclude either root individually.

## The regular modulus comparison

On a common two-term chart, the normalized exact polynomial is
\[
 L(\xi)(1+e_j(\xi))
       +j^\gamma A(\xi)e^{jD(\xi)}(1+f_j(\xi)),
       \qquad \gamma=p/2,\qquad |e_j|+|f_j|\le K_0/j.        \tag{1}
\]
Here \(L,A,D\) are independent of the order for a fixed profile.
In a convex subchart with \(\Re D'\ge c>0\), choose one analytic
logarithm \(C=\log(L/A)\). Let
\[
                  E_j=\log(1+e_j)-\log(1+f_j).
\]
For large \(j\), \(E_j\) has computable uniform bound \(2K_0/j\).
At an exact root,
\[
 j\Re D(\xi_j)=-\gamma\log j+\Re C(\xi_j)+\Re E_j(\xi_j).
                                                               \tag{2}
\]
Let \(\delta=\xi_N-\xi_n\), which is real.
By (2), \(|\Re D(\xi_n)|\le(\gamma\log n+\|C\|+2K_0)/n\).
Subtracting (2) gives, with computable \(K_1\),
\[
       (Nc-\|C'\|)|\delta|
                       \le K_1(1+\log n)/n.
\]
To obtain this bound, write the left phase difference as
\(N\Re(D(\xi_N)-D(\xi_n))+\Re D(\xi_n)\).
The first difference has modulus at least \(Nc|\delta|\).
On the right use \(\log(N/n)\le1/n\), the derivative bound for
\(C\), and the two error bounds. Thus the computable estimate
\[
                         |\delta|\le K_2(1+\log n)/n^2        \tag{3}
\]
holds uniformly above a terminal threshold.

At a limiting equal-height pair, either one field has a one-term
witness or both fields are at the same balance point, by the uniqueness
proved in the
[regular comparison](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_REGULAR_COMPLEMENT.md#linked-comparison-and-termination).
At a balance point \(\Re D'= \Re q>0\), supplying the convex
neighborhood of (1) and the positive real slope in (2).

## Interior heights and the negative real boundary

For all the limiting profiles of the theorem, the upper saddle
\(q(\xi)=D'(\xi)\) has \(\Re q>0\) in \(0<\Im\xi<\pi\), by the
[axis exclusion and continuation](MULTIVARIATE_UNIFORM_ONE_PRIVATE_REGULAR_COMPLEMENT.md#the-positive-real-boundary-and-the-adjacent-comparison).
Its imaginary part is positive there by the one-sheet construction.
On \(\Im\xi=\pi\), the saddle continues to real \(q>1\) in the
right gap. The upper logarithm formula for its phase gives
\(\Im D=\pi q-\pi(q-1)=\pi\) at this boundary. All other terms
in the phase are real at such \(q\). In particular, for each fixed
real \(x\) the identity
\(\partial_y\Im D(x+iy)=\Re q(x+iy)>0\) gives \(\Im D(x+iy)<\pi\)
at every interior height.

Let \(x_0\) be the positive-real coexistence field in the
centered logarithm coordinate. The simple right saddle there
has \(D(x_0)=0\) and continues from the upper sheet. For \(0<y<\pi\),
\(\partial_y\Re D(x_0+iy)=-\Im q(x_0+iy)<0\),
so \(\Re D(x_0+iy)<0\).
Strict increase of \(\Re D\) with the real coordinate shows
that any balance at this height has real coordinate \(x>x_0\).
For every \(x>x_0\) the upper sheet continues at \(y=0\) to the
simple real right maximizer. This follows from the real fold
geometry in the regular proof and the implicit function theorem,
since \(X'(q)<0\) at that maximizer, for \(X+\xi=\Phi_q\).
At these real fields \(D(x)\) is real. Integrating
\(\partial_y\Im D=\Re q>0\) now shows that every interior
balance has \(0<\Im D<\pi\).

Thus every interior or negative-boundary balance has
\(e^D\ne1\). Such a strict margin can be checked by
interval subdivision on a smaller chart. Dividing the two
exact zero equations of (1) gives
\[
 \exp\bigl(D(\xi_n)+N(D(\xi_N)-D(\xi_n))\bigr)
    =(n/N)^\gamma\exp\bigl(C(\xi_N)-C(\xi_n)
                                  +E_N(\xi_N)-E_n(\xi_n)\bigr).
\]
By (3) and the derivative bounds, this implies
\(|e^{D(\xi_n)}-1|\le K_3(1+\log n)/n\) with computable \(K_3\).
For sufficiently large \(n\) this contradicts the chart's strict
positive lower bound on \(|e^D-1|\).

## Positive-real coexistence

At the remaining balance point \(x_0\), the
[positive-boundary expansion](MULTIVARIATE_EFFECTIVE_ONE_PRIVATE_REGULAR_COMPLEMENT.md#positive-boundary-witnesses)
gives (1) on a symmetric field disk. We make the analytic errors
real on the real diameter. For \(L(1+e_j)\), take the exact
short initial lattice sum. Its complex uniform approximation
by \(L\) is the boundary estimate of the cited expansion.
Take the remaining polynomial as the other packet, including
the exponentially small complement outside the retained saddle
interval. The saddle lattice expansion and strict exterior
loss are uniform on a sufficiently small disk, giving the stated
relative \(O(j^{-1})\) bound for \(f_j\). Both packets are
holomorphic and real on the real diameter because they are
finite sums of the original normalized terms with real coefficients
in \(e^\xi\). The leading stationary value is real on that diameter
and the leading Gaussian amplitude \(A\) there is positive.
The endpoint leading amplitude \(L\) is positive as well.
Thus the branches of \(C\) and \(E_j\) can be chosen real on
the real diameter. The strict margins used for these decompositions
hold on a larger disk, so Cauchy estimates on a smaller symmetric
rectangle give computable constants \(C_2^\ast,D_2^\ast,E_1^\ast\)
bounding \(|C''|,|D''|,j|E'_j|\), respectively. All bounds persist
on a small profile cell where the real branches continue.

Write \(\xi_j=x_j+iy\) in that rectangle. Real roots on the positive
ray are excluded by the positive coefficients. For \(y>0\), the
reality properties and integration along vertical segments give
\[
 \begin{aligned}
 c y&\le\Im D(x_j+iy)\le\|D'\|y,\\
 |\Im D'(x_j+iy)|&\le D_2^\ast y,\qquad
 |\Im C'(x_j+iy)|\le C_2^\ast y,\qquad
 |\Im E_j(x_j+iy)|\le E_1^\ast y/j .
 \end{aligned}
\]
The derivative bounds hold at every point of the horizontal joining
segment. Taking imaginary parts of the logarithmic zero equations
and subtracting shows that
\[
 \begin{split}
 2\pi(m_N-m_n)
   = \Im D(\xi_n) + N\Im\bigl(D(\xi_N)-D(\xi_n)\bigr)\\
       -\Im\bigl(C(\xi_N)-C(\xi_n)+E_N(\xi_N)-E_n(\xi_n)\bigr)
 \end{split}                                                    \tag{4}
\]
for integers \(m_n,m_N\).
The last two terms combined have modulus at most
\[
          y\bigl((ND_2^\ast+C_2^\ast)|\delta|
                       +2E_1^\ast/n\bigr).
\]
By (3) that coefficient is less than \(c/2\) above a computable
order. Shrink the rectangle height so that
\(y(\|D'\|+c/2)<2\pi\). Then the right side of (4) lies
strictly between zero and \(2\pi\), excluding the pair.

## Computable cutoff

The same compact profile nets, field rectangle, and candidate
contours used for the effective regular comparison enumerate
the witnesses above. At a positive-real balance use the real
interval and saddle witnesses in a symmetric disk, the positive
leading amplitudes, and the derivative bounds on the inner
rectangle. The smaller chart can be tested by computable strict
inequalities. At other balances use the nonzero margin for
\(e^D-1\) alongside the existing contour witnesses. Pair
uniqueness and the one-term exclusions leave no other cases in
these regular bands. All witnesses apply on profile cells and
field pairs by strict margins, and compactness gives a finite
subcover. The maximum of their explicit local order thresholds
is a computable uniform cutoff.
