# Joint positive moment shapes do not exclude consecutive common roots

For every adjacent pair \(n,N=n+1\) with \(n\ge4\), there are
successive moment polynomials
\[
                  P_j(z)=\mathbb E(U+zV)^j,\qquad j=n,N ,
\]
of the same finite, strictly positive joint law of \(U,V\),
with the following simultaneous properties.

- Their coefficients and the probabilities of the law are
  rational. They share an irreducible quadratic with nonreal
  roots arbitrarily close to the negative real axis.
- On the full fractional grid, both normalized log moment
  symbols are strictly convex, and the upper-minus-lower
  symbol is strictly convex. Their curvatures can all be
  made arbitrarily small.
- In either orientation the direct coefficient ratio is
  strictly increasing and strictly log convex. The
  normalized curvature \((N-x)^2 f''(x)\) of its logarithm
  is strictly increasing on \(0\le x<N\), as in the
  [factorial direct-ratio comparison](MULTIVARIATE_DIRECT_RATIO_CURVATURE_AND_SECTOR_COLLAR.md).

For every fixed \(R\ge1\), and all sufficiently large adjacent
orders, one can also impose positivity on all even derivatives
of orders \(2,4,\ldots,2R\) of both normalized symbols and
of their difference throughout the fractional grid.
This rules out an implication from those finite derivative
signs, the direct ratio shape, and the common positive moment
law alone. The laws here do not come from ordinary factorial
exponent rows.

## Symmetric atoms with a prescribed moment zero

For \(j\ge4\) let \(p_j(t)\) be the real part of \((t-i)^j\)
when \(j\) is even and its imaginary part when \(j\) is odd.
A symmetric measure on \(t\in\mathbb R\) has
\(\mathbb E(t-i)^j=0\) precisely when \(\mathbb E p_j(t)=0\),
since the opposite component cancels under \(t\mapsto -t\).
The real function of the positive magnitudes to integrate is
\(p_j(t)\) on \(t>0\), or equivalently its value on one
representative of each symmetric pair.
More explicitly the average of the two powers at \(t,-t\)
is \(p_j(t)\) if \(j\) is even and \(i p_j(t)\) if \(j\) is odd.

Put \(\vartheta=\arctan(1/t)\) for \(t>0\). Apart from the
positive modulus factor, the relevant polynomial is
\(\cos(j\vartheta)\) in even degree and
\(-\sin(j\vartheta)\) in odd degree. Thus its positive
zeros are simple, and occur at magnitudes
\[
  t=\cot\frac{(2k+1)\pi}{2j}
       \quad(j\ \mathrm{even}),\qquad
  t=\cot\frac{k\pi}{j}
       \quad(j\ \mathrm{odd}),
\]
retaining angles strictly between zero and \(\pi/2\).
The smallest positive zero is
\(s_j=\tan(\pi/(2j))<1/2\). The largest zero \(T_j\) is
\(\cot(\pi/(2j))>1\) for even \(j\) and
\(\cot(\pi/j)>1\) for odd \(j\ge5\).

Set \(v_n=1\) and \(v_N=1-1/(2N)\).
A symmetric law supported on magnitudes \(s_j,T_j\)
can realize variance \(v_j\) by assigning probability
\[
                       p_j^*=\frac{v_j-s_j^2}{T_j^2-s_j^2}
\]
to the larger magnitude. Both weights are positive.
It has the desired complex moment zero because both
magnitudes are zeros of the real polynomial above.

The law may be made rational while keeping that moment
equation and the variance exact. At each of \(s_j,T_j\)
choose two nearby rational magnitudes bracketing its simple
zero. Their \(p_j\) values have opposite signs, so a
positive rational mixture on each pair makes its mean
\(p_j\) exactly zero. The resulting second moments
straddle \(v_j\) when the brackets are sufficiently small.
Mix these two laws with positive rational weights to get
variance \(v_j\). Their four support columns for the
functions \(1,t^2,p_j(t)\) have rank three. Indeed the
two line segments in the \((t^2,p_j(t))\)-plane cross
the horizontal axis at two distinct points and each has
a nonzero vertical component.

At each \(j\), add small positive rational mass on all
remaining magnitudes from the union of the two supports.
Correct the original weights at three independent columns
to preserve the mean of \(1,t^2,p_j(t)\). The corrections
tend to zero with the added mass, so all weights stay
positive. Give the positive and negative representatives
of each magnitude equal weight. This yields symmetric
rational laws \(\nu_n,\nu_N\) with the same finite support
and satisfies
\(\mathbb E_{\nu_j}(t-i)^j=0\), \(\mathbb E_{\nu_j}t^2=v_j\).
The construction can approximate the initial two-magnitude
laws as closely as desired in every fixed moment.

## The same positive joint law and its exact coefficients

Choose a small positive rational \(\varepsilon\) so that
\(X(t)=1+\varepsilon t>0\) on the common support.
Write \(p_{j,t}>0\) for its mass at \(t\) in law \(\nu_j\).
On the same support put
\[
   U(t)=p_{N,t}/p_{n,t},\qquad V(t)=X(t)U(t),\qquad
   \mu(t)=\frac{p_{n,t}U(t)^{-n}}{\sum_s p_{n,s}U(s)^{-n}} .
\]
This is a single positive rational joint law of \(U,V\).
With the common positive scalar
\(C=(\sum_s p_{n,s}U(s)^{-n})^{-1}\) its two successive
moment polynomials are exactly
\[
           P_j(z)=C\,\mathbb E_{\nu_j}(1+zX(t))^j,
                  \qquad j=n,N .
\]
Put \(z_*=-1/(1+i\varepsilon)\).
Then \(1+z_* X(t)=z_*\varepsilon(t-i)\), so
\(P_n(z_*)=P_N(z_*)=0\). Rationality and conjugation
give the common quadratic
\[
                        (1+\varepsilon^2)z^2+2z+1,
\]
which is irreducible over \(\mathbb Q\).
Every coefficient is strictly positive and the even
polynomial is strictly positive for real \(z\), because
its support has more than one value of \(X\).

## Convexity, the paired increment, and the direct ratio

Define the unnormalized log coefficient symbols
\[
 F_j(x)=\log\mathbb E_{\nu_j}X(t)^x,\qquad
 L_j(s)=j^{-1}\bigl[F_j(js)-sF_j(j)\bigr] .
\]
The variables \(x,s\) are real. Symmetry gives, uniformly
with any fixed number of \(x\) derivatives on
\([-1,N]\),
\[
                  F_j(x)
                  =\frac{\varepsilon^2 v_j}{2}x(x-1)
                                                +O(\varepsilon^4).
\]
Indeed the sum before the logarithm has constant term
one, vanishing odd powers of \(\varepsilon\), and the
displayed second-order term from the binomial expansion.
The finite support gives the claimed uniform Taylor
remainders on this compact \(x\) interval.

Thus \(F_j''(x)=\varepsilon^2 v_j+O(\varepsilon^4)\)
and \(F_j'''(x)=O(\varepsilon^4)\).
In particular \(L_j''(s)>0\) and
\[
                 L_N''(s)-L_n''(s)
                       =\tfrac12\varepsilon^2+O(\varepsilon^4)>0
\]
uniformly for \(0\le s\le1\).
For the direct coefficient ratio,
\[
            f(x)=\log\frac{N}{N-x}+F_N(x)-F_n(x),
                       \qquad 0\le x<N .
\]
The difference of the two \(F_j''\) is
\(-\varepsilon^2/(2N)+O(\varepsilon^4)\).
Therefore both \(f'\) and \(f''\) are strictly positive
for sufficiently small \(\varepsilon\), as their respective
binomial terms are at least \(1/N\) and \(1/N^2\).
The derivative of the normalized curvature is
\[
 \begin{aligned}
  \frac{d}{dx}\bigl[(N-x)^2 f''(x)\bigr]
     &=(N-x)
       \left[-2(F_N''-F_n'')+(N-x)(F_N'''-F_n''')\right]\\
     &=(N-x)\bigl[\varepsilon^2/N+O(\varepsilon^4)\bigr]>0 .
 \end{aligned}
\]
Reversal replaces the difference of \(F_j(x)\) by
\(F_N(N-x)-F_n(n-x)\) and at most a constant.
The same uniform estimates still hold, including the
signs and the normalized curvature derivative.
All smallness statements apply to positive rational
\(\varepsilon\) in some nonempty interval.

## Positivity through an arbitrary finite even order

Fix \(R\) and let \(n,N\) be sufficiently large. For this
part use a different large zero \(T_j\) if needed.
When \(n\) is even use \(T_n=\cot(3\pi/(2n))\) and
\(T_N=\cot(\pi/N)\). When \(n\) is odd use
\(T_n=\cot(\pi/n)\) and \(T_N=\cot(\pi/(2N))\).
Keep the same small zero and target variances as above.
These are still zeros of the appropriate real polynomial
and straddle the target variances. The ratios \(T_N/T_n\)
tend to \(3/2\) and \(2\) in the two cases.

For the initial two-magnitude law, its \(2\ell\)-th
moment is \(v_j T_j^{2\ell-2}(1+O_\ell(j^{-2}))\)
for \(\ell\ge2\). Indeed \(T_j\) is comparable to \(j\),
\(s_j=O(j^{-1})\), and the larger-magnitude weight is
\((v_j-s_j^2)/(T_j^2-s_j^2)\).
Every term in the \(2\ell\)-th cumulant other than
the \(2\ell\)-th moment is a product of at least two
positive-order even moments. This follows from expanding
the logarithm of the moment generating series and using
symmetry. Each such term is \(O_\ell(T_j^{2\ell-4})\)
because the \(2r\)-th moment is \(O_r(T_j^{2r-2})\).
Writing \(\kappa_{j,2\ell}\) for the even cumulants gives
\[
 \kappa_{j,2\ell}=v_j T_j^{2\ell-2}
                         (1+O_\ell(j^{-2})),\qquad \ell\ge2 .
\]
For every \(1\le\ell\le R\), both cumulants are positive
and
\(N^{2\ell-1}\kappa_{N,2\ell}>n^{2\ell-1}\kappa_{n,2\ell}\)
eventually. For \(\ell=1\) this is exactly the
chosen variance relation. For \(\ell\ge2\) it follows
from the strict limiting ratio of the large zeros.
The rational bracket and common-support corrections can
be chosen to preserve all these strict inequalities.

Finally \(F_j(x)\) is the log moment generating function
at \(x\) of \(\log(1+\varepsilon t)\).
Its \(2\ell\)-th derivative has uniform expansion
\[
         F_j^{(2\ell)}(x)
               =\varepsilon^{2\ell}\kappa_{j,2\ell}
                                      +O(\varepsilon^{2\ell+2}).
\]
For the leading term, note that \(\varepsilon^{-1}
\log(1+\varepsilon t)\) tends to \(t\) on the finite
support, so differentiation of its log moment generating
function at \(\varepsilon x\) tends to the indicated
cumulant. Analyticity in \(\varepsilon\) and symmetry
remove the next odd term. Since
\(L_j^{(2\ell)}(s)=j^{2\ell-1}F_j^{(2\ell)}(js)\),
the positive cumulant margins imply all the claimed
derivative signs for small positive rational
\(\varepsilon\), uniformly over \(0\le s\le1\).
The same exact common quadratic persists throughout.

## An absolutely continuous variant

The same obstruction also holds for an absolutely continuous
positive joint law with rational polynomial coefficients.
Replace each common positive magnitude \(t\) by a uniform
law on a sufficiently short rational interval around \(t\),
and replace its negative magnitude by the reflected law.
Choose the intervals disjoint. The averaged feature columns
for \(1,t^2,p_j(t)\) are rational and converge to the
original columns. Apply the rank-three correction to the
interval weights to preserve those three expectations.
All weights remain positive, and every strict even cumulant
condition persists by continuity. The density ratio
\(d\nu_N/d\nu_n\) is the ratio of the two weights inside
each interval, so it remains positive and constant there.
The same construction of \(U,V,\mu\) works with densities.
Finally multiply both variables by an independent positive
uniform random variable \(W\) on \([1,2]\). This multiplies
each moment polynomial by a positive rational scalar,
preserving all roots and shape conclusions.
Inside each interval the map from \((t,W)\) to
\((WU(t),WV(t))\) has Jacobian modulus
\(\varepsilon W U(t)^2>0\). The resulting law is
absolutely continuous, as claimed.

## A simple common pair with an exact affine Euler factor

For every terminal \(N\ge4\) there is a second obstruction with
an exact polar derivative relation. Write \(n=N-1\) and
\(\mathcal D=z\,d/dz\). There are positive rational moment
polynomials \(A,B\) in degrees \(n,N\), from one finite
positive joint law, satisfying
\[
    B-\frac1N\mathcal D B=(1+\alpha\mathcal D)A,\qquad
                    \frac1N B'=R(\mathcal D)A                    \tag{7}
\]
for a positive rational \(\alpha\) and a real polynomial
\(R\). The function \(R(k)/(1+\alpha k)\) strictly increases
on \(k=0,\ldots,n\). Nevertheless there is a common nonreal
conjugate pair, simple in both parents. Their matched symbol
curvature is strictly positive. The direct coefficient ratios
have the strictly increasing normalized curvature of the
analytic interpolation
\(\log[N/(N-x)]+\log(1+\alpha x)\).
Equation (7) matches one exact affine factor and the monotone
quotient in the
[ordinary Euler relations](MULTIVARIATE_LINKED_GAMMA_ANGLE_COMPACTIFICATION.md#euler-operator-relations-for-integer-exponent-rows).
It does not require the complementary multiplier \(R\) to
factor into the prescribed affine row factors.

To construct the pair, write \(N=2r+\delta\) with
\(\delta\in\{0,1\}\). The symmetric density on \(\mathbb R\)
proportional to \((1+t^2)^{-(r+\delta+1)}\) has odd moments
zero, and
\[
                  \mathbb E T^{2k}
                          =\frac{\binom r k}{\binom N{2k}}
                          \qquad(0\le k\le r).
\]
In fact substituting \(u=t^2\) in the integral gives the
ratio of beta integrals
\(B(k+1/2,r+\delta+1/2-k)/B(1/2,r+\delta+1/2)\).
Simplification using the half-integer gamma recurrence
gives the displayed binomial ratio. All moments up to \(N\)
are finite.

This moment vector has a rational quadrature with positive
weights and full affine rank. For clarity, the expectation
of \((t,\ldots,t^N)\) lies in the closed convex hull of this
curve by truncation and approximation. It lies in its
interior because any nonzero polynomial of degree at most
\(N\) that is nonnegative on \(\mathbb R\) has positive
expectation under the strictly positive density.
The rational points on the curve have the same closed
convex hull. Approximate the vertices of a small simplex
surrounding the expectation by points in the convex hull
of rational curve points. This places the moment vector
in the interior of a finite rational-node hull. A small
positive initial weight on every chosen node leaves a
renormalized target still inside that hull. Thus the
original target has a representation with all weights
positive, and the weights can be chosen rational by
Gaussian elimination and rational approximation in the
solution space. Reflecting and averaging makes the
quadrature symmetric. Its affine rank is preserved.

For sufficiently small rational \(\varepsilon>0\), all
transformed nodes \(X=1+\varepsilon T\) are positive.
They still give full moment rank, since this is an
invertible affine change of variable. Hence all sufficiently
near finite moment data of orders through \(N\) can be
represented by small weight corrections on these same
positive nodes, preserving strictly positive weights.
The base upper polynomial given by this quadrature is
\[
       \mathbb E(1+zX)^N=(1+z)^\delta Q(z)^r,\qquad
                      Q(z)=(1+z)^2+\varepsilon^2 z^2 .
\]
This follows by binomial expansion and the exact \(T\)
moments. Put \(z_*=-1/(1+i\varepsilon)\), a simple root of
\(Q\). If \(r>2\), add a sufficiently small nonzero rational
multiple of \(zQ(z)^2\) to the base polynomial. Its order
at \(z_*\) becomes exactly two, with constant coefficient
one and all stripped moments still in the interior
just proved. If \(r=2\), no perturbation is needed.
Denote the resulting double-root polynomial by \(B_0\).

For a variable \(\alpha\) near zero, correct the
coefficients at \(z,z^2,z^3,z^4\) of \(B_0\) to obtain
\(B_\alpha\) with
\[
 B_\alpha(z_*)=A_\alpha(z_*)=0,\qquad
   A_\alpha=(1+\alpha\mathcal D)^{-1}
                          (B_\alpha-\mathcal D B_\alpha/N).
\]
The inverse multiplies coefficient \(k\) by \(1/(1+\alpha k)\).
At \(\alpha=0\) the linear system on the four real
corrections is equivalent to imposing value and first
derivative at \(z_*,\bar z_*\). A polynomial in the four
chosen monomials is \(z\) times a cubic, so it cannot vanish
twice at both distinct nonzero points unless it is zero.
The system is therefore nonsingular at zero and remains
so for small \(\alpha\). The corrected polynomials depend
analytically on \(\alpha\), with rational coefficients
when \(\alpha\) is rational. Their constant terms stay one.

A dot denotes the derivative at \(\alpha=0\).
Differentiating \(B_\alpha(z_*)=0\) and the displayed
equation for \(A_\alpha\) gives
\[
   0=\dot A(z_*)=-\frac{z_*}{N}\dot B'(z_*)
                           +\frac{z_*^2}{N}B_0''(z_*).
\]
Thus \(\dot B'(z_*)=z_*B_0''(z_*)\ne0\), so \(B_\alpha\)
has a simple root at \(z_*\) for all small nonzero \(\alpha\).
Also \(A_0'(z_*)=-z_* B_0''(z_*)/N\ne0\), so the lower
root stays simple. Take \(\alpha>0\) rational and small.

For the stripped upper coefficients
\(b_k=[z^k]B_\alpha/\binom Nk\), the stripped lower
coefficients are \(b_k/(1+\alpha k)\).
Both arrays can be realized by positive rational
weights on the fixed positive nodes, by the full
moment rank and closeness to the original common data.
The previous construction of a joint law from the
two weight lists then represents \(A_\alpha,B_\alpha\)
up to a common scalar. The matched symbol increment
remains strictly positive. Indeed for the unperturbed
common node law \(F(x)=\log\mathbb E X^x\),
\(F''(x)=\varepsilon^2\mathbb E T^2+O(\varepsilon^4)\)
uniformly on \([0,N]\). Thus
\(N F''(Ns)-n F''(ns)>0\) throughout the grid
for small \(\varepsilon\). The subsequent two
coefficient perturbations preserve this strict margin
through their small weight corrections.

The complementary multiplier in (7) is determined on
the sample indices by
\[
                    R(k)=(1+\alpha k)\frac{b_{k+1}}{b_k}.
\]
Strict log convexity of the positive nonconstant moment
law for \(b_k\) makes \(b_{k+1}/b_k\) strictly increasing.
Finally the exact direct coefficient ratio at \(k\le n\)
is \(N(1+\alpha k)/(N-k)\). Its stated interpolation has
normalized curvature
\(1-\alpha^2(N-x)^2/(1+\alpha x)^2\), which strictly increases
on \([0,N)\). Thus the common simple roots survive an exact
affine polar factor together with the positive joint law.
