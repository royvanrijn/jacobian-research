# Higher lifting obstructions for homogeneous bivariate heights

Fix \(s\ge2\) distinct offsets \((\lambda_j,D-\lambda_j)\) of
one total degree and a set \(I\) of \(s\) distinct positive integer
orders. Add independent common heights \(L_x,L_y\). Put
\[
 \bar\lambda=s^{-1}\sum_j\lambda_j,\quad
 \mu_j=\lambda_j-\bar\lambda,\quad e_3=e_3(\mu_1,\ldots,\mu_s),
 \qquad H_x=L_x+\bar\lambda,\quad H_y=L_y+D-\bar\lambda .
\]
Here \(e_3=0\) for \(s=2\). Write \(\ell=s-1\),
\(\sigma^2=H_x^{-1}+H_y^{-1}\) and \(r=H_x^{-1}/\sigma^2\).
The Gaussian system consists of
\(G_m(F)=\mathbb E F(Z_m)^m\) for \(m\in I\), where
\(Z_m\sim\mathcal N(0,1/m)\) and \(\deg F\le\ell\).
At a common projective zero of full degree \(\ell\) where the
projective Jacobian has rank \(s-1\), choose a left null vector \(u\)
and set
\[
 \Theta=\sum_{m\in I}u_m \mathbb E[F(Z_m)^{m-1}F'(Z_m)],\qquad
 \Phi=\sum_{m\in I}u_m \mathbb E[Z_m^2 F(Z_m)^m].
\]
Use the monic chart for \(F\).

**Theorem.** If \(\Theta\ne0\), any family of exact common zeros
with both heights tending to infinity and centered Gaussian
coefficient limit \(F\) obeys
\[
 \frac{r-1/2}{\sigma}\longrightarrow\frac{9\Phi}{64\Theta}.
 \tag{1}
\]
In particular a nonreal \(\Phi/\Theta\) excludes such a family.
If \(\Phi=0\), the sharper necessary condition is
\[
 L_x-L_y\longrightarrow
 D_*:=D-2\bar\lambda+\frac{12e_3}{s(s+1)(s+2)}.
 \tag{2}
\]
At integer heights this excludes the family if \(D_*\) is nonintegral,
and otherwise confines it eventually to \(L_x-L_y=D_*\).

These are fixed-support, fixed-window statements with complex
coefficients. Section 4 gives a finite decision under additional
hypotheses, with a nonvanishing conclusion over every
characteristic-zero field at ordinary integer heights.

## 1. Centered phase and change of law

In the factorial integral at order \(m\), absorb the common
factor \(X^{mH_x}Y^{mH_y}\) into the density on the positive
orthant. The coordinates are then independent gamma variables
of shapes \(mH_x+1,mH_y+1\) and unit scale. Divide by the
nonzero gamma factors and use the random variable
\[
 Z=\sigma^{-1}\bigl(\log(X/(mH_x))-\log(Y/(mH_y))\bigr).
\]
The residual monomial of index \(j\) equals
\((H_x/H_y)^{\mu_j}e^{\sigma\mu_j Z}\). Absorb the first
factor into its coefficient. This change is independent of \(m\)
and valid also for real heights in the positive corner.

Choose \(w_{kj}\) by the invertible Vandermonde system
\(\sum_j w_{kj}\mu_j^t/t!=\delta_{kt}\) for \(0\le t\le\ell\).
In the exact basis \(E_k(z;\sigma)=
\sigma^{-k}\sum_j w_{kj}e^{\sigma\mu_j z}\), let
\(F(z)=\sum_{k=0}^{\ell}a_k z^k\) denote the limiting phase for
coefficients \(a\). At fixed \(a\) the phase has expansion
\[
 F+\sigma^2 F_2+\sigma^3 F_3+O(\sigma^4),\qquad
 \deg F_i\le\ell+i,\qquad
 [z^{\ell+3}]F_3=\frac{a_\ell e_3}{s(s+1)(s+2)}.
 \tag{3}
\]
For \(k<\ell\) the linear correction is zero by construction.
For \(k=\ell\), the coefficient of \(v^{s-1}\) in the
reduction of \(v^s\) modulo \(\prod_j(v-\mu_j)\) is
\(\sum\mu_j=0\). The same reduction of \(v^{s+2}\) has
coefficient \(e_3\), since \(v^s\equiv -e_2v^{s-2}
+e_3v^{s-3}+\cdots\) for \(s\ge3\). For \(s=2\),
\(v^4\) reduces to a constant. The factorial ratio from the
Vandermonde normalization gives the coefficient in (3).

Let \(D_z\) denote differentiation. To second order, the
change from the Gaussian expectation on polynomial tests
is substitution of
\[
 1+\sigma T_1+\sigma^2(T_2+T_1^2/2),\qquad
 T_1=(2r-1)\left(\frac{D_z}{2m}-\frac{D_z^3}{6m^2}\right),\quad
 T_2=-\frac{r^2+(1-r)^2}{4m^2}D_z^2
        +\frac{r^3+(1-r)^3}{12m^3}D_z^4 .
 \tag{4}
\]
The third-order law term on each polynomial test has
coefficients polynomial in \(r\) and vanishes at \(r=1/2\).
For details, write \(C=mH_i\). The log-gamma cumulants
before scaling are \(\psi(C+1)-\log C\) for order one
and \(\psi_{k-1}(C+1)\) for order \(k\ge2\).
The [Stirling expansion with remainder](https://dlmf.nist.gov/5.11)
and the [polygamma asymptotics](https://dlmf.nist.gov/5.15)
give \(1/(2C)+O(C^{-2})\), \(C^{-1}-1/(2C^2)+O(C^{-3})\),
\(-C^{-2}+O(C^{-3})\), and \(2C^{-3}+O(C^{-4})\),
respectively, at orders one through four. Expanding orders
one and three one term further and including order five
gives the terms of degree three in \(\sigma\), with
polynomial coefficients in \(r\) after \(C^{-1}=\sigma^2r_i/m\)
where \(r_i=r,1-r\). Subtract the second coordinate for odd
orders and add it for even orders. All odd cumulants vanish
at \(r=1/2\), while even cumulants have only even powers
of \(\sigma\). Orders \(k\ge6\) contribute first at
\(\sigma^{k-2}\) or later. The finite moment-cumulant
formula now gives (4) and the asserted third-order property.

For fixed orders the remainder is uniform for
\(0<r<1\) and bounded \(a\). For the polynomial expectations
this follows from \(C^{-1}\le\sigma^2/m\) and the cumulant
expansion. To use (3) inside the integral, the gamma ratio
formula for \(\mathbb E e^{tZ}\) gives uniform bounds at
every fixed real \(t\). Indeed the log ratio for each
coordinate has a linear term \(O(\sigma)\) after subtracting
\(t\log C/\sigma\), and the Taylor remainder of degree one
is \(O(\sigma^{-2}\sup_{|v-C|\le |t|/\sigma+1}v^{-1})=O(1)\)
by the trigamma estimate. Thus polynomial powers times
\(e^{c|Z|}\) have bounded expectations for each fixed power
and \(c\), using the bounds at \(t=\pm(c+1)\). Taylor
expansion of the exponential basis through degree \(k+3\)
then gives the stated remainder in every fixed moment norm.

## 2. Projection at a simple zero

Let \(F\) now be the limiting zero. Define the linear functional
\(\mathcal J(P)=\sum_m u_m m\mathbb E[F^{m-1}P]\),
with the law at each summand understood. The null-vector
condition and Euler's identity give \(\mathcal J(P)=0\)
for \(\deg P\le\ell\). Gaussian integration by parts gives
\[
\begin{aligned}
 \mathcal J(zF)&=\mathcal J(F')=0,\\
 \mathcal J(z^2F)&=\sum_m u_m G_m(F)+\mathcal J(zF')=0,\\
 \mathcal J(z^3F)&=2\Theta+\mathcal J(z^2F')=2\Theta,\\
 \mathcal J(z^4F)&=3\Phi+\mathcal J(z^3F')=3\Phi .
\end{aligned}
\]
The first two lines extend the kernel to all degrees at most
\(\ell+2\), since \(F\) has degree \(\ell\). This explains
both cancellations in the last two lines.

For \(f=F^m\), integration by parts rewrites \(\mathbb E T_1 f\)
as \((2r-1)(\mathbb E[Z_m f]-(m/6)\mathbb E[Z_m^3 f])\).
Its projection by \(u\) is \((2/3)(2r-1)\Theta\).
At \(r=1/2\) the second-order law term is
\(-\mathbb E f''/(8m^2)+\mathbb E f^{(4)}/(48m^3)\).
Using \(G_m(F)=0\), this is
\(-\mathbb E[Z_m^2f]/4+(m/48)\mathbb E[Z_m^4f]\).
Its projection is \(-3\Phi/16\).
The phase adds \(\mathcal J(F_2)=0\).
At \(r=1/2\) the law has no third-order term, and the
first-order law cannot act on \(F_2\) there. By (3) the
third-order projected phase is
\[
 \mathcal J(F_3)=\frac{2e_3}{s(s+1)(s+2)}\Theta,
\]
since replacing the top term of \(F_3\) by the corresponding
multiple of \(z^3F\) changes it by a polynomial of degree
at most \(\ell+2\).

## 3. Necessary limits for an escaping family

Use the monic chart and let \(\delta\) be the coefficient deviation
from the limiting zero. The rank assumption gives a local lower
bound \(c\|\delta\|\) for the norm of the selected \(s-1\) Gaussian
equations, by Taylor expansion with invertible linear term.
At exact zeros (3)–(4) therefore give
\(\|\delta\|=O(\sigma|2r-1|+\sigma^2)\).
Project the full equations by \(u\); the Gaussian part contributes
\(O(\|\delta\|^2)\). The first-order projected coefficient
is \((2/3)(2r-1)\Theta\) up to an
\(O(|2r-1|\|\delta\|)\) error. Hence \(r-1/2=O(\sigma)\)
when \(\Theta\ne0\), and then \(\|\delta\|=O(\sigma^2)\).
Put \(h=r-1/2\). The second-order coefficient at the
limiting zero is an even polynomial in \(h\) by (3)–(4).
The third-order coefficient at \(h=0\) was computed in
Section 2, and its change for \(h=O(\sigma)\) is \(O(\sigma)\).
The projected equation becomes
\[
 0=\frac43\sigma h\Theta-\frac3{16}\sigma^2\Phi
       +\sigma^3\frac{2e_3}{s(s+1)(s+2)}\Theta
       +O(\sigma^4).
 \tag{5}
\]
For the \(\sigma^4\) remainder, the error from the first-order
coefficient is \(O(\sigma|h|\|\delta\|)\), and the second-order
errors are \(O(\sigma^2(h^2+\|\delta\|))\).
This proves (1). If \(\Phi=0\), (5) gives
\(h/\sigma^2\to-3e_3/(2s(s+1)(s+2))\). Now use
\[
 H_x-H_y=\frac{1-2r}{\sigma^2r(1-r)}
\]
and the definitions of the centered heights to get (2).

For a nonzero real ratio in (1), put \(c=9\Phi/(64\Theta)\),
\(S=L_x+L_y+D=H_x+H_y\), and \(\Delta=L_x-L_y\).
Equation (5) still gives the next term in \(h\), so that
\[
 \Delta=D_*-\frac{8c}{\sigma}+O(\sigma).
\]
Since \(r(1-r)=1/4-h^2\) and
\(\sigma^{-2}=S r(1-r)\), we have
\(\sigma^{-1}=\sqrt{S}/2+O(1/\sqrt S)\). In particular,
\[
  (\Delta-D_*)^2-\frac{81}{256}(\Phi/\Theta)^2 S=O(1).
 \tag{6}
\]
When the squared ratio is rational, the left side lies in
a fixed rational lattice at ordinary integer heights.

## 4. Arithmetic consequences and finite decision

For a trinomial let \(a=\lambda_2-\lambda_1\) and
\(b=\lambda_3-\lambda_2\) in increasing order.
Then \(3\mu_j\) are \(-2a-b,a-b,a+2b\), and
\(D_*=D-2\bar\lambda+e_3/5\). If \(a\not\equiv b\bmod3\),
each of the three integers \(3\mu_j\) is a unit modulo 3.
The term \(e_3/5\) has denominator divisible by \(27\),
whereas the other term has denominator dividing 3.
Thus \(D_*\) is nonintegral. If \(a\equiv b\bmod3\), the
centered exponents are integers and \(D_*\) is integral
exactly when one of them is divisible by 5. These tests
apply on the branch with \(\Phi=0\).

**Finite decision.** Suppose all projective Gaussian zeros have
full degree, rank \(s-1\), and \(\Theta\ne0\). Suppose further
that at each zero the ratio \(\Phi/\Theta\) is nonreal or has
rational square. Then there is a finite decision for the existence of
a complex common factorial zero, at any nonnegative integer heights,
for the prescribed support and orders.

The Gaussian zeros are isolated by the rank condition and
form a finite set by projective compactness. Compute them
over the algebraic numbers, together with \(\Phi/\Theta\).
The ratio is independent of the scaling of the null vector
and the phase. Each nonreal ratio is excluded by (1).
Choose a rational radius for a strip about \(\Delta=D_*\)
less than one and, for noninteger \(D_*\), less than its
distance to \(\mathbb Z\). A root with zero ratio leads to
this strip by (2). For each nonzero real ratio, put
\(\alpha=81(\Phi/\Theta)^2/256\), a positive rational
number under the hypothesis. Restrict to a strip
\(\lvert(\Delta-D_*)^2-\alpha S\rvert\le B\)
with a separate bound for each such \(\alpha\).
Compactness of projective coefficient space and the uniform
normalization imply that some finite \(B\) and a cutoff
cover all real-height common zeros with both heights above
the cutoff. The proof of (5) uses bounds uniform in a fixed
neighborhood of each root. Thus the \(O(1)\) in (6) is
uniform near that root for the real-height common zeros.

For an effective search for the bounds and cutoff, use real
quantifier elimination in the complement of these strips,
trying increasing positive integers for the cutoff and \(B\).
The stripped factorial forms represent the same equations
for real heights: after removing
\(\Gamma(mL_x+1)\Gamma(mL_y+1)\), the remaining gamma ratios
are polynomial rising factorials. Equivalently one can use
the [exact polynomial normalization](SFC_BIVARIATE_TRANSLATION_GAUSSIAN_JETS.md#1-exact-graded-normalization).
An empty union of strips simply requires a corner exclusion.

The uncovered integer heights lie on finitely many bounded-height
rows and columns. The zero-ratio strip contributes at most the
diagonal \(\Delta=D_*\). Each positive-ratio strip contributes
finitely many rational parabolas
\((\Delta-D_*)^2-\alpha S=b\), since the possible values \(b\)
lie in a lattice of known denominator and satisfy \(|b|\le B\).
Parameterize each parabola by the integer \(\Delta\), using
\(S=((\Delta-D_*)^2-b)/\alpha\). The two heights are rational
polynomials in \(\Delta\), so their integrality is a finite
congruence test after clearing denominators. Their
nonnegativity is a univariate inequality test.
Parameterize the rows, columns and diagonal by their free
integer height. On each resulting univariate family,
restrict the homogeneous resultant, and test the integer
zeros or, if it vanishes identically, the full admissible
parameter set. For exact support one can instead eliminate
on the coefficient torus; its univariate parameter image
is constructible and hence finite or cofinite, with a
computable list of exceptional parameters. At ordinary
integer heights all forms are rational, so a complex
solution exists exactly when an algebraic one does.
Consequently any nonvanishing certificate holds over
every characteristic-zero field.

The general odd-start quadratic Gaussian problem remains open.
A deficient Jacobian, vanishing \(\Theta\), or a nonzero real
ratio with irrational square can prevent this decision from applying.
