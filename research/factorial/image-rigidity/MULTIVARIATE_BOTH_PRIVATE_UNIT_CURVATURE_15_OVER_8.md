# Private unit pairs at shared curvature fifteen eighths

An ordinary two-monomial profile with exactly one private
unit variable on each side has disjoint consecutive root arguments
in every degree if its shared absolute curvature satisfies
\[
                        S=\sum_{\rm shared}
                     \frac{(b_i-a_i)^2}{\max(a_i,b_i)}\ge15/8.
\]
The shared rows may have arbitrary orientations.
We strengthen the [private unit threshold \(77/40\)](MULTIVARIATE_BOTH_PRIVATE_UNIT_CURVATURE_77_OVER_40.md)
using some exact shared-row coefficient attenuation before applying
its paired estimates. The same initial terminal results apply.

## A stronger coefficient cap from one row

Let a shared row have smaller exponent \(\ell\), gap \(d\ge1\),
and larger exponent \(b=\ell+d\). At order \(m\) its
endpoint-normalized log factor at \(k=mw\) is the negative
Green integral of its continuous second derivative. For
\(0<w\le1/2\), write
\[
                  G(w,x)=\min(w,x)(1-\max(w,x)).
\]
Under the probability density \(2G(w,x)/[w(1-w)]\) on \([0,1]\),
the mean of \(x\) is \((1+w)/3\), by integrating the two
linear pieces. Thus the mean affine row exponent is at most
\(b-d/3\) in either orientation. The trigamma is convex and
decreasing, and \(\psi_1(y)\ge1/y\) by integral comparison.
Jensen's inequality in the Green integral therefore bounds the
row's normalized coefficient by
\[
              \exp\left(
                 -\frac{d^2}{b-d/3+1/m}\frac{k(1-k/m)}{2}
                  \right).
\]
For \(m\ge13\), the replacement curvature
\[
                             \frac{d^2}{b-d/3+1/13}             \tag{1}
\]
is thus valid for the two coefficient halves. The derivative
estimate for the actual linked coefficients still uses
\(V=\sum d_i^2/\ell_i\); (1) changes only the upper caps.

## The remaining band \(15/8\le S<77/40\)

If a row has \(d>\ell\), it is \((\ell,d)=(1,2)\) or
\((2,3)\), as in the \(77/40\) proof. In these cases (1)
increases its curvature by \(46/141\) or \(108/265\),
respectively. Thus the existing \(s=77/40\), \(v=16/3\),
\(C=21/20\) cap and derivative calculation applies, with
\(s\) the coefficient-envelope curvature, \(v\) the
upper bound on \(V\), and \(e^{C/m}\) the radius.
The packing bound for \(V\) is unchanged, since \(S<2\).

Now assume every ratio \(d/\ell\le1\). If some ratio is at least
\(1/2\), (1) has gain at least \(1/40\). For \(d=1\), the
possible smaller exponents are one and two, and the gains
are \(5/68\) and \(10/321\), respectively.
For \(d\ge2\), the absolute row curvature is at least \(2/3\).
The denominator in (1) divided by \(b\) is at most
\(1-1/9+1/39=107/117\), using \(d/b\ge1/3\) and \(b\ge3\).
This gives gain at least \(20/321>1/40\).
Consequently we can use the cap curvature \(s=19/10\),
with \(V\le2S<4\). If every ratio is below \(1/2\),
we have \(V\le3S/2<3\) and can use \(s=15/8\) directly.

These two cases use the following scalar rows in the
\(77/40\) proof's replay, with all angle and order ranges
unchanged.
\[
\begin{array}{c|c|c}
 s&v&C\\ \hline
 15/8&3&11/10\\
 19/10&4&43/40
\end{array}
\]
In particular the same Beta quotient orders the first-box
phases before their midpoint correction. The proof uses the
new cap symmetrically for both orientations and the unchanged
linked derivative and quotient estimates.

## Exact replay

Execute the scalar replay in
[MULTIVARIATE_BOTH_PRIVATE_UNIT_CURVATURE_77_OVER_40.md](MULTIVARIATE_BOTH_PRIVATE_UNIT_CURVATURE_77_OVER_40.md#exact-scalar-replay)
with its final band list replaced by
\[
 (s,v,C)\in
 \left\{(15/8,3,11/10),\ (19/10,4,43/40)\right\}.
\]
Both rows have each of the squared Rouché ratio, the bounded-order
angle factor, and the uniform large-order angle factor strictly
below one.
