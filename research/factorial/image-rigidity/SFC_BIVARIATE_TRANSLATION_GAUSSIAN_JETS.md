# Gaussian limits and lifting tests for bivariate common translations

Fix any bivariate monomial offsets \(b_j=(\lambda_j,d_j-\lambda_j)\),
with \(s\) distinct columns. Add the same independent heights \(L_x,L_y\)
to all columns. Fix a set \(I\) of \(s\) distinct positive moment orders.
For a total-degree class \(g\), let \(d_g\) be its degree and \(r_g\)
its number of columns.

**Theorem.** There are rational homogeneous forms
\(\mathcal F_m(\sigma,r,a)\) in \(s\) coefficient variables \(a\),
with coefficients polynomial in \(\sigma,r\), that represent the
original common zero problem for
\[
 \sigma^2=L_x^{-1}+L_y^{-1},\qquad
 r=\frac{L_x^{-1}}{L_x^{-1}+L_y^{-1}}
\]
when \(L_x,L_y>0\). At \(\sigma=0\) they are independent of \(r\)
and are the Gaussian moment forms
\[
 G_m(a)=\mathbb E P_m(Z_m)^m,\qquad
 P_m(z)=\sum_g\sum_{k=0}^{r_g-1}a_{g,k}m^{d_g}z^k,\qquad
 Z_m\sim\mathcal N(0,1/m).
 \tag{1}
\]
Thus if the Gaussian forms have no common projective zero, all
independent height failures of the prescribed window lie on
finitely many computable horizontal or vertical fibers.
The exact resultant test on those fibers gives a finite decision.

At a Gaussian zero where the projective Jacobian has rank \(s-1\),
the first perturbation gives an additional lifting condition.
For a homogeneous support of size \(s\), at a simple Gaussian
zero of degree \(s-1\) with nonzero transverse scalar as specified
below, any escaping family converging to it must have
\(L_x/L_y\to1\). These are statements at fixed offsets and fixed
orders, over every characteristic-zero field at the ordinary
integer heights. The Gaussian gate is an explicit finite
polynomial system for each prescribed window.

## 1. Exact graded normalization

Use the
[Gamma normalization of the independent heights](SFC_INDEPENDENT_COMMON_HEIGHT_DECISION.md#1-exact-normalization-and-the-case-of-distinct-degrees),
with \(Y_i\) of mean \(1+1/(mL_i)\) at order \(m\).
For each degree class choose a basis \(U_{g,k}\), \(0\le k<r_g\),
in its monomial span, satisfying
\[
 U_{g,k}(Y)=
   Y_2^{d_g-k}(Y_1-Y_2)^k+
   \sum_{\ell=k+1}^{d_g}
           c_{g,k,\ell}Y_2^{d_g-\ell}(Y_1-Y_2)^\ell .
 \tag{2}
\]
Indeed, express \(Y_1^\lambda Y_2^{d_g-\lambda}\) in the powers
of \(Y_1-Y_2\). The first \(r_g\) coefficients are
\(\binom{\lambda}{\ell}\), \(0\le\ell<r_g\). Their
evaluation matrix at the distinct \(\lambda\)'s is invertible
by the Vandermonde determinant. Therefore that initial block
can be made the identity. For \(\ell<r_g\) every off-diagonal
entry in (2) is zero.

Change the normalized phase to
\(\sum_{g,k}a_{g,k}m^{d_g}U_{g,k}(Y)/\sigma^k\).
This is an invertible coefficient change independent of \(m\)
for every \(\sigma>0\). Write
\[
 Z=\frac{Y_1-Y_2}{\sigma},\qquad W=\frac{Y_2-1}{\sigma}.
\]
The divided bases are polynomials in \(\sigma,Z,W\),
since (2) contains only powers \(\ell\ge k\).
For \(i=1,2\), let \(r_1=r,r_2=1-r\). The first cumulant
of \((Y_i-1)/\sigma\) is \(\sigma r_i/m\), and its cumulant
of order \(\ell\ge2\) is
\[
         (\ell-1)!m^{1-\ell}\sigma^{\ell-2}r_i^{\ell-1}
                       (1+\sigma^2r_i/m).
\]
Mixed cumulants of the coordinates vanish by independence.
All polynomial moments in \(Z,W\) are therefore polynomial
in \(\sigma,r\) by the finite moment-cumulant formula.
This gives the exact forms in the theorem, even at the
boundary values \(r=0,1\) by polynomial evaluation. At
\(\sigma=0\), the joint Gaussian has
\[
 \operatorname{Var}(Z)=1/m,\qquad
             \operatorname{Cov}(W,Z)=-(1-r)/m .
\]
Consequently the constant term in \(\sigma\) is (1).

If the Gaussian gate has no projective zero, compactness of
the unit sphere in the coefficient variables and of
\(0\le r\le1\) gives a positive interval of \(\sigma\) on
which the exact forms remain nonvanishing. One can search
for a rational interval using real quantifier elimination
on real and imaginary parts. In two variables the complement
in the ordinary height plane is covered by finitely many
fibers as in the
[independent-height decision](SFC_INDEPENDENT_COMMON_HEIGHT_DECISION.md#4-effective-search-in-two-independent-heights).
It is enough to test over \(\mathbb C\), since a common
direction for the original rational moment equations over
any characteristic-zero field gives one over
\(\overline{\mathbb Q}\) by the Nullstellensatz in a
coefficient chart.

## 2. The first lifting obstruction

Define
\[
 V_m(z)=\sum_{g,k}(d_g-k)a_{g,k}m^{d_g}z^k,\qquad
 U_m(z)=\sum_{g,k}c_{g,k,k+1}a_{g,k}m^{d_g}z^{k+1},
\]
taking a missing coefficient as zero. The first derivative
of the exact form at \(\sigma=0\) is
\[
 \mathcal B_m(r,a)
 =m\,\mathbb E\!\left[
   P_m(Z_m)^{m-1}
   \left(U_m(Z_m)-(1-r)Z_m V_m(Z_m)
            +\frac{2r-1}{3}Z_m^3 P_m(Z_m)\right)\right].
 \tag{3}
\]
In particular it is affine in the relative height \(r\).

Indeed the divided phase before taking the moment is
\(P_m(Z)+\sigma(W V_m(Z)+U_m(Z))+O(\sigma^2)\), using (2).
For the first correction one substitutes the limiting
Gaussian and uses \(\mathbb E(W\mid Z)=-(1-r)Z\).
For the change of law in the leading term
\(f(Z)=P_m(Z)^m\), the first and third cumulants of \(Z\)
are \(\sigma(2r-1)/m\) and
\(2\sigma(2r-1)/m^2+O(\sigma^3)\), respectively.
The variance has no linear term, and higher cumulants
are \(O(\sigma^2)\). The moment-cumulant formula, or the
formal exponential of the corresponding differential
operators, gives the first change
\[
      (2r-1)\left(\frac{\mathbb E f'(Z_m)}m+
                   \frac{\mathbb E f'''(Z_m)}{3m^2}\right)
      = \frac{m(2r-1)}3\,\mathbb E[Z_m^3 f(Z_m)] .
\]
The equality follows by Gaussian integration by parts.
This proves (3).

Suppose \([a_*]\) is a common Gaussian zero and that the
Jacobian of \((G_m)_{m\in I}\) in a projective chart has
rank \(s-1\). Write \(u\) for a nonzero left null vector of
that Jacobian. If exact common zeros approach \([a_*]\) as
all heights grow, with \(r\to r_*\), then
\[
                           u\,\mathcal B(r_*,a_*)=0.          \tag{4}
\]
To prove it, choose \(s-1\) equations with invertible
Jacobian in the chart. By the implicit function theorem
they fix the nearby coefficient point as an analytic
function of \((\sigma,r)\), equal to \(a_*\) at
\(\sigma=0\). Its deviation is \(O(\sigma)\) uniformly
near \(r_*\). Multiply the first-order expansion of the
full system by \(u\), divide by \(\sigma\), and pass to
the limit. This proves (4). In particular, if all Gaussian
zeros have the stated rank and none satisfies (4) for a
real \(r_*\in[0,1]\), then the same finite independent-height
decision applies. The rank and last nonvanishing condition
are polynomial tests over real algebraic numbers.

## 3. A homogeneous support and the central ratio

Suppose there is just one degree class and the limiting
Gaussian polynomial \(F\) has degree \(s-1>0\).
Remove the nonzero common order factors so that
\(P_m=F\) at every order. In this normalization write
\[
 \Theta=\sum_{m\in I}u_m\,\mathbb E[F(Z_m)^{m-1}F'(Z_m)] .
\]
At a simple common zero as in Section 2, condition (4) becomes
\[
                              (2r_*-1)\Theta=0.              \tag{5}
\]
For the proof let \(\ell=s-1\) and set
\[
 \mathcal J(P)=\sum_{m\in I}u_m m\,\mathbb E[F(Z_m)^{m-1}P(Z_m)] .
\]
The projective null-vector condition and Euler's identity
at a homogeneous zero imply that \(\mathcal J\) kills all
polynomials of degree at most \(\ell\), including the
coefficient direction omitted by the chart.
Gaussian integration by parts gives
\(\mathcal J(zF)=\mathcal J(F')=0\).
The leading coefficient of \(F\) is nonzero, so
\(\mathcal J\) also kills every polynomial of degree at most
\(\ell+1\). This eliminates \(U_m-(1-r)z V_m\) from the
projected perturbation once the factors \(m^{d_g}\) have
been removed. Using integration by parts on \(z^2 F^m\) gives
\[
       \mathcal J(z^3F)
          =2\sum_{m\in I}u_m \mathbb E[Z_m F(Z_m)^m]
                         +\mathcal J(z^2F')
          =2\Theta .
\]
In the last equality the \(\mathcal J\) term has degree at
most \(\ell+1\), and \(m\,\mathbb E[Z_m F^m]=
m\,\mathbb E[F^{m-1}F']\).
Substitute in (3) to obtain (5).
Thus \(\Theta\ne0\) forces \(r_*=1/2\).
In fact near this root the analytic implicit solution of
the \(s-1\) selected equations gives, after projection, the
equation \(\sigma(2(2r-1)\Theta/3+O(\sigma))=0\) with
uniform error for \(r\) near \(1/2\). Hence an escaping
family has \(r-1/2=O(\sigma)\), or
\[
                   |L_x-L_y|=O\big(\sqrt{L_x+L_y}\big).
\]
The radius of this strip depends on the root, offsets and orders.

For homogeneous trinomials, a consecutive-window Gaussian
zero necessarily has degree two by
[scaled Hermite isolation](SFC_SCALED_HERMITE_ROOT_ISOLATION.md)
for the lower-degree possibilities. The first lifting test
therefore applies to hypothetical odd-window roots of the
[quadratic Gaussian gate](SFC_ODD_GAUSSIAN_LOCAL_WINDOW_TESTS.md)
that satisfy the rank and scalar tests. In the equal-ratio
direction (5) leaves a higher-order problem.
