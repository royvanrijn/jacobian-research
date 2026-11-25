# Matched shared-coefficient rays at terminal nine

Let all unequal rows of distinct ordinary exponent columns be shared.
Write \(\ell_i=\min(a_i,b_i)\), \(d_i=|a_i-b_i|\), and let \(u\)
be the normalized order-eight outer reflected mean. With
\[
       \frac9{100}\le u\le\frac{44}{100},
       \qquad \max_i \frac{d_i}{8\ell_i+1}\le\frac18,
\]
the eighth and ninth factorial moment polynomials have disjoint root
arguments under independent positive rescalings. In particular they
are coprime over every characteristic-zero field.

## Matched skew parameter

For a matched Gauss factor let \(x=d/(8\ell+h)\),
\(x'=d/(9\ell+h)\), \(a=h/d\), with \(1\le h\le d\).
At order \(j\) write \(s_j=\log(1+(j-1)x_j)/2\) and
\(T_j=\tanh s_j=(j-1)x_j/(2+(j-1)x_j)\), where
\(x_8=x,\ x_9=x'\). The reciprocal caps are \(1/8,8/71\),
because \(x'=8x/(9-ax)\). Thus the matched tanh values satisfy
\[
 \frac{T_9}{T_8}
 =\frac{64}{63}\left(1+\frac{(2a-1)x}{18+(64-2a)x}\right).
\]
The correction in parentheses differs from one by at most
\(\eta=1/206\), using \(0\le a\le1\) and \(x\le1/8\).
Put \(R=64/63\), \(U_j=\sum s_j\), and denote by \(T_{*,j},s_{*,j}\)
the corresponding values at the reciprocal cap.

Fix a choice of sign for the orientation of each shared row so
that the signed skew at the reflected coefficient with lower
index is given by the signed Taylor vectors in the
[private linear proof](MULTIVARIATE_TERMINAL_NINE_PRIVATE_LINEAR.md#grouping-the-shared-logarithmic-corrections).
These orientation signs \(\epsilon_i\) agree at both orders for
matched factors. Define their first signed moments
\[
 P_j=\frac12\sum_i\epsilon_i T_{j,i}^2,
 \qquad
 P_{*,j}=\frac{T_{*,j}^2}{2s_{*,j}} .
\]
For every \(q\ge2\) the function \(T^q/\operatorname{artanh}T\) is
increasing on \(0<T\le T_{*,j}\).
Indeed the derivative sign is that of
\(q(1-T^2)\operatorname{artanh}T-T\), which is nonnegative since
\(\operatorname{artanh}T\ge T\) and \(T_{*,j}^2<1/2\).
Consequently \(\frac12\sum T_{j,i}^2\le P_{*,j}U_j\).
The squared transfer bounds give
\[
 \left|P_9-R^2(1+\eta^2)P_8\right|
       \le 2R^2\eta P_{*,8}U_8 .
\]
Exchanging the exponent columns changes the signs of both \(P_j\).
It replaces each root \(z\) by \(1/z\) up to positive rescaling,
which preserves the set of upper root arguments by conjugation.
We can therefore use the shared scalar
\(\lambda=P_8/(P_{*,8}U_8)\in[0,1]\).

For a higher skew vector indexed by \(l=2,3\), its scalar is
\(\sum_i\epsilon_i T_{j,i}^{2l}/(2l)\). Put
\(c_{j,l}=T_{*,j}^{2l-2}/(2l)\). The same scalar equals \(c_{j,l}P_j\)
plus an error with absolute value at most \(c_{j,l}P_{*,j}U_j\).
In fact \(T_{j,i}^{2l-2}\) differs from \(T_{*,j}^{2l-2}/2\)
by at most \(T_{*,j}^{2l-2}/2\). Add the two leading terms
\(c_{j,l}P_j\) to the first skew vector, using one scalar \(P_j\)
for the combined direction and two bounded symmetric residuals.

The first two reflected mean-defect vectors have scalars
\(D_{j,l}=\sum_i T_{j,i}^{2l+1}/(2l+1)\), for \(l=1,2\).
Their upper bounds are \(U_jT_{*,j}^{2l+1}/((2l+1)s_{*,j})\).
For the lower bound Jensen applied to the convex function
\(y\mapsto y^{2l}\) with weights \(T_{j,i}/\sum_i T_{j,i}\)
gives
\[
 D_{j,l}\ge\frac{(2|P_j|)^{2l}}{(2l+1)U_j^{2l-1}} .
\]
Here \(\sum T_{j,i}\le U_j\) and \(\sum T_{j,i}^2\ge2|P_j|\).
The grouped checker uses the implied lower bound whenever the
interval for \(P_j\) has positive lower endpoint.
The bounds on the omitted skew and mean tails are those of the
[geometric series estimate](MULTIVARIATE_TERMINAL_NINE_PRIVATE_LINEAR.md#grouping-the-shared-logarithmic-corrections).

## Exact projection and continuation

The [checker](verify_multivariate_terminal_nine_shared_matched_light.py)
uses the exact mean-ratio bounds for
\(\rho=U_9/U_8\) from the
[outer-mean lemma](MULTIVARIATE_TERMINAL_NINE_PRIVATE_ROWS.md#shared-outer-mean-boxes).
On a closed rectangle in \(u,\lambda\), insert the reference
Gaussian coefficients at the center of \(u\) and the nominal
combined skew at the centers of \(u,\lambda\). The shifted
order-nine nominal skew is \(R^2(1+\eta^2)\) times the order-eight one.
The exponential of each nominal shift is enclosed by its Taylor
series through degree twelve, using the geometric tail
\(A^{13}/(13!(1-A))\) for absolute shift at most \(A<1\).
Include the difference from this nominal skew as a grouped scalar,
together with the common Gaussian shift, the two bounded skew
residuals and two mean defects. Use the quadratic exponential
model with explicit cubic and omitted logarithmic errors from the
[private linear proof](MULTIVARIATE_TERMINAL_NINE_PRIVATE_LINEAR.md#grouping-the-shared-logarithmic-corrections).
The same model applies to envelopes larger than the actual profiles,
with its coefficient-vector contributions multiplied by independent
intervals for each scalar and each quadratic scalar product.

The seven cosines are fixed piecewise affine functions of \(u\).
Subdivide \(u,\lambda\) and the radial interval until a strict
Bernstein projection sign is verified on each cell. Cosine
intervals use the polynomial sine-ratio recurrence and exact
Bernstein checks from the
[terminal-eight paired bridge](MULTIVARIATE_TERMINAL_EIGHT_PAIRED_BRIDGE.md#paired-projections).
Root radii remain between \(1/12\) and \(12\), by leading and
constant coefficient domination using the independent coefficient
caps at the upper endpoint of the proved range.

At fixed \(u,\rho\) first scale every actual skew and mean defect
by a common factor \(t\) decreasing from one to zero. The transfer,
skew, and tail upper bounds remain valid, with \(\lambda\) replaced
by \(t\lambda\). Each lower mean-defect bound remains valid because
\(D_{j,l}\) is scaled by \(t\) and its lower bound at \(tP_j\)
is scaled by \(t^{2l}\). At zero skew and defects, move \(\rho\)
to \(64/63\) in its permitted interval. The resulting Gaussian
references then move along the range to \(u=9/100\).
The independent coefficient caps hold throughout, since the first
homotopy interpolates the coefficient logarithms between the actual
shape and the Gaussian shape. The annulus and ray exclusions
therefore give the root counts of the Gaussian anchor.

At \(u=9/100\), exact rational Rouché disks isolate the
roots of the two Gaussian reference arrays. They give one root in
each upper sector, with orders \(9,8,9,8,9,8,9,8\) in increasing
argument, and one additional negative real root of order nine.
Join the first and last sectors to their conjugates. Only the
single ninth-order root in the last joined sector is real, because
the actual eighth-order factorial moment is strictly positive on
the real axis by its Gamma integral. This excludes common root
arguments. Rational Bézout transfers to all characteristic-zero
fields.

Replay from the repository root with

~~~sh
python3 research/factorial/image-rigidity/verify_multivariate_terminal_nine_shared_matched_light.py
~~~

The replay verifies 20,351 radial cells and both Gaussian anchor arrays.
