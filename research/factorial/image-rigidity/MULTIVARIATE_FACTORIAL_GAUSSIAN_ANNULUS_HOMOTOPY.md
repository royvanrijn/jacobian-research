# A radial bound preserved by the factorial-to-Gaussian path

Consider two distinct ordinary exponent columns and a consecutive
pair \(n,N=n+1\) with \(n\ge2\). Endpoint-normalize their factorial
moments. Write \(U_0,W_0\) for the lower and upper reflected outer
losses and \(R=n^2/(n^2-1)\). The result below concerns positive
Gamma shapes along the path. At its starting point the columns
are ordinary.

**Theorem.** There is a continuous family of positive normalized
coefficient pairs from the factorial pair to the Gaussian pair
with losses \(U_0,RU_0\). It keeps the lower outer loss \(U_0\)
fixed. The upper loss \(W\) satisfies \(W\ge RU_0\). Every
initial lower reciprocal cap remains valid for the surviving
Gauss factors, and the usual matched upper ratio bound for that
cap remains valid for the total losses. At each stage and both
orders, every root of the normalized polynomial satisfies
\[
                     u_j^2\le |z|\le u_j^{-2},
 \qquad u_j=\exp(-U_j)
\]
for its own outer loss \(U_j\). If all surviving tanh parameters
at order \(j\) are at most \(b\), the stronger annulus is
\([u_j^{\,1+b},u_j^{-1-b}]\).

## Factor bounds and endpoint normalization

For an increasing Gauss factor at order \(j\) let its shape be
\(A>0\). The moments before normalization have multiplier
\((A)_k\). Their roots have moduli between \(1/(A+j-1)\) and
\(1/A\) by the [single-factor Gamma annulus](MULTIVARIATE_PAIRED_GAMMA_ROOT_ANNULUS.md#matched-factor-bounds-for-the-ordinary-profiles).
Endpoint normalization multiplies root moduli by
\(\mu=(A)_j^{1/j}\in[A,A+j-1]\). Thus this factor contributes
normalized radii between \(e^{-2s}\) and \(e^{2s}\), where
\(s=\tfrac12\log((A+j-1)/A)\). Reversal gives the same bound for
a decreasing factor.

Concavity of the logarithm also yields
\[
 s\le \log(\mu/A)\le s+\log\cosh s.
\]
The left inequality averages the chord joining the first and
last Gauss terms. The right inequality takes the logarithm of
the arithmetic mean of the \(j\) terms \(A,A+1,\ldots,A+j-1\).
Now \(\log\cosh s\le s\tanh s\le sb\). The same normalized
factor bounds are therefore within
\([e^{-(1+b)s}, e^{(1+b)s}]\), in either orientation.
Binomial-basis composition multiplies the lower and upper radii
by the [disc composition lemma](GAP_TWO_GRACE_GLOBAL_LOCALIZATION.md#1-a-self-contained-closed-disc-convolution-lemma).
Applying it also after reversal gives the product annulus.
Gaussian factors have unit-circle roots and so add no radial
width. [Gaussian reference proof](MULTIVARIATE_GAUSSIAN_REFERENCE_ANGULAR_INTERLACING.md#root-motion-on-the-upper-semicircle)

## A path with fixed lower outer loss

On each unequal row write \(\ell\) for the smaller exponent and
\(d\) for its positive integer gap. For \(v\ge0\) replace \(\ell\)
by \(\ell+v\). The associated shape at order \(j\) is
\((j(\ell+v)+h)/d\) for \(1\le h\le d\). Keep the same orientations
and match the same integer shift \(h\). The arguments for the
[outer ratio](MULTIVARIATE_MATCHED_OUTER_RATIO_MONOTONICITY.md)
and the [excess-loss transport](MULTIVARIATE_TANH_EXCESS_LOSS_TRANSPORT.md)
use \(\ell\ge0\) and integer \(d,h\), so they hold on this path.
Every row reciprocal decreases as \(v\) increases.

Let \(U(v),W(v)\) be the two outer losses from the surviving
Gauss factors. Compose at lower order with the normalized Gaussian
of outer loss \(U_0-U(v)\) and at upper order with the normalized
Gaussian of outer loss \(R(U_0-U(v))\). A zero-loss Gaussian here
means \((1+z)^j\). In binomial basis, the reflected outer losses
add and the Gaussian contributes only the Gaussian logarithmic
direction. Thus the total losses are
\[
       U_0,\qquad W(v)+R(U_0-U(v)).
\]
Their excess is \(W(v)-R U(v)\), exactly the surviving factorial
excess. The total ratio is a weighted mean of \(R\) and \(W(v)/U(v)\).
This proves the claimed lower and cap-based upper ratio bounds.
The surviving signed and unsigned higher moments obey the
matched bounds with \(U(v),W(v)\). Upper estimates that increase
with the outer losses remain valid with the total losses. For
the Jensen lower estimates the total losses enlarge the
denominators. The coupled lower ratio bound follows with the same
signed moment since \(U(v)\le U_0\).

The radius estimate from the preceding section uses the sum of
the surviving Gauss losses. Enlarging this to the total loss
retains the bound at either order. As \(v\to\infty\) the normalized
Gauss factors tend to the binomial multipliers and their losses
tend to zero. The Gaussian multipliers tend to those for losses
\(U_0,RU_0\), completing the path.

A separating-ray proof can use the tighter annulus without
changing the angular counts argument. After endpoint
normalization a column reversal reciprocates the roots and
preserves the angular questions under conjugation. For the Cayley
inverse image of a separator with \(\Re z<0\), it therefore
suffices to check the half corresponding to \(|z|\le1\) and to
discard its portion with \(|z|<u_j^2\). The starting integer-row upper ratio bounds also remain valid
along this path by the following monotonicity.

## The row ratio decreases under inflation

For fixed integer \(d\) and \(0<y\le1/2\), the integral \(J\)
in the matched outer-ratio proof can be written
\[
 J(y,\ell)=\int_0^\infty e^{-\ell v}A_y(v)\,dv,\qquad
 A_y(v)=
 \frac{2e^{-(d+y)v/2}\sinh((1-y)dv/2)\sinh(ydv/2)}
      {(1-y)dv\,\sinh(yv/2)}.
\]
This follows by integrating \(e^{-v(1-y)ds}\) over \(0\le s\le1\)
and summing \(e^{-vyh}\) for \(1\le h\le d\). Let
\(H(z)=\coth z-z\operatorname{csch}^2 z\). It is positive
and increasing on \(z>0\), since \(H(0+)=0\) and
\(H'(z)=2\operatorname{csch}^2 z\,(z\coth z-1)>0\).
Direct differentiation gives
\[
 \partial_v\partial_y\log A_y(v)
 =-\frac12-\frac d2
   \left[H((1-y)dv/2)-H(ydv/2)\right]-\frac12 H(yv/2)<0 .
\]
Hence \(A_{1/N}(v)/A_{1/n}(v)\) is strictly increasing in \(v\).
The quotient \(J(1/N,\ell)/J(1/n,\ell)\) is the expectation of
this increasing function under the probability density proportional
to \(e^{-\ell v}A_{1/n}(v)\). Its \(\ell\)-derivative is minus
the covariance of this function with \(v\), which is strictly
negative by the two-variable covariance identity. The Laplace
densities and their first moments decay exponentially also at
\(\ell=0\), justifying differentiation. The row ratio \(W/U\)
is \(R J(1/N,\ell)/J(1/n,\ell)\), so it decreases under inflation.
In particular the upper ratio bound obtained from the original
integer rows on a fixed lower-mean band persists along the entire
path. Any refined reciprocal cap valid at the start likewise
persists because the reciprocals decrease.
