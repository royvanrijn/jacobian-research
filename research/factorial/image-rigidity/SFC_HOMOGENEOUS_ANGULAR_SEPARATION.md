# Consecutive factorial windows for angularly separated homogeneous supports

Let \(D\) be a common total degree in \(k\) variables. For two exponent
columns \(e,f\in\mathbb R_{\ge0}^k\) of total \(D\), use the distance
\[
 \mathfrak d_D(e,f)
 =2\sqrt{D+k}\arccos\!\left(
      \frac{\sum_{\ell=1}^k\sqrt{(e_\ell+1)(f_\ell+1)}}{D+k}
                         \right).
 \tag{1}
\]
This is spherical distance on the square roots of the shifted and
normalized exponent vectors. For \(k=2\) and \(e=(u,D-u)\), put
\[
 \xi_D(u)=2\sqrt{D+2}\arcsin\sqrt{\frac{u+1}{D+2}} .
\]
Then (1) is \(|\xi_D(u)-\xi_D(v)|\).

**Trinomial sign theorem.** For homogeneous bivariate trinomials with
real nonnegative common translation heights, the sign \(R<0\) of
the [initial signed remainder](SFC_HOMOGENEOUS_TRINOMIAL_SHORT_SPANS.md)
holds whenever both adjacent distances in (1) are at least \(12\).
In particular the sign holds if the successive exponent gaps
\(a,d-a\) satisfy \(\min(a,d-a)\ge6\sqrt{D+2}\), where \(D=L+H+d\).

**Separated window theorem.** For each fixed number of variables \(k\ge2\),
support size \(r\ge2\), and starting order \(n\ge1\), there is a computable
\(T=T(k,r,n)\) with the following property. A homogeneous \(r\)-term
polynomial with integer exponent columns separated pairwise by at
least \(T\) in (1) has a nonzero ordinary factorial moment among
\(n,\ldots,n+r-1\), provided all exponent entries are at least \(T\).
In two variables the same conclusion holds for a threshold \(T(2,r,n)\)
with no lower bound on exponent entries. The assertion about
nonvanishing is over any characteristic-zero coefficient field.
The proof gives a procedure to compute the threshold from \(k,r,n\).

The angular reduction is the one used in
[Liu and Sun, Theorem 2.6](https://doi.org/10.1017/S0004972719000546),
where the homogeneous bivariate theorem assumes vanishing at all
positive orders. The
[uniform trinomial corner](SFC_UNIFORM_HOMOGENEOUS_TRINOMIAL_INITIAL_CORNER.md)
allows coalescing angular nodes for three initial moments and proves
nonvanishing under a different large-exponent condition. The results
here use pairwise separation to obtain the specific remainder sign
and to handle arbitrary fixed window length and start.

## 1. The mixed-moment estimate

On the simplex \(t_\ell\ge0,\ \sum t_\ell=1\), with measure given
by the first \(k-1\) coordinates, put
\[
 Z_m(e)=\int t^{me}\,dt
       =\frac{\prod_\ell\Gamma(me_\ell+1)}{\Gamma(mD+k)},\qquad
 f_e=\frac{t^e}{Z_1(e)},\qquad J_m(e)=\frac{Z_m(e)}{Z_1(e)^m}.
\]
These formulas also hold for real nonnegative exponents. The
change \(x_\ell=\rho t_\ell\) in the factorial integral removes a
positive row factor \(\Gamma(mD+k)\). Rescaling the column
coefficients by their first-order simplex integrals leaves the
functions \(f_e\) and preserves any common zero of the moment forms.

For a tuple \(e_1,\ldots,e_m\) of total-\(D\) columns with diameter
\(W\) in (1), the mixed moment satisfies
\[
 \int\prod_{i=1}^m f_{e_i}
 \ \le\ e^{-W^2/8}\prod_{i=1}^m J_m(e_i)^{1/m}.
 \tag{2}
\]
Indeed, if \(\bar e=m^{-1}\sum_i e_i\), the logarithm of the
left side divided by the product on the right is
\(\log Z_m(\bar e)-m^{-1}\sum_i\log Z_m(e_i)\).
The [trigamma series](https://dlmf.nist.gov/5.15.E1) gives \(\psi_1(u)\ge1/u\) for \(u>0\)
by comparison with its tail integral. Thus the Hessian of
\(\log Z_m\), computed before restriction to the fixed-total
hyperplane with its denominator held constant, dominates \(m\)
times the Hessian of
\[
 H(e)=\sum_\ell(e_\ell+1)\log(e_\ell+1).
\]
Consequently the negative of the log ratio is at least
\(\sum_i H(e_i)-mH(\bar e)\). This last quantity is the sum of
the divergences from the tangent plane to \(H\) at \(\bar e\).
For \(p=(e+1)/(D+k),\ q=(f+1)/(D+k)\), each such divergence is
\((D+k)\sum_\ell p_\ell\log(p_\ell/q_\ell)\).
Concavity of log bounds the relative entropy below by
\(-2\log\sum_\ell\sqrt{p_\ell q_\ell}\). If \(\theta\) is the
corresponding angle, \(-2\log\cos\theta\ge\theta^2\), by
integrating \(\tan\theta\ge\theta\).
The divergence is therefore at least
\(\mathfrak d_D(e,f)^2/4\). For two columns realizing the
diameter, the sum of their squared distances to \(\bar e\) is
at least \(W^2/2\) by the triangle inequality. This proves (2).

## 2. The signed trinomial remainder

For \(k=2\), the normalized diagonal moments satisfy a useful
uniform comparison. Write
\[
 w(u)=\sqrt{\frac{(u+1)(D-u+1)}{D+2}},\qquad
 C_m=\frac{m^{m-1}}{(m-1)!}.
\]
At any fixed \(D,m\) there is a positive number \(c_{m,D}\)
independent of \(u\) such that
\[
 c_{m,D}w(u)^{1-m}\le J_m(u)\le C_m^2 c_{m,D}w(u)^{1-m}.
 \tag{3}
\]
To prove this, the [gamma multiplication formula](https://dlmf.nist.gov/5.5#E6)
factors \(\Gamma(mu+1)/\Gamma(u+1)^m\) into a constant
times \(m^{mu}\prod_{j=1}^{m-1}\Gamma(u+j/m)/\Gamma(u+1)\).
The exponential factor cancels its \(u\)-dependence when combined
with the other exponent. For \(0<h<1\), log convexity of Gamma
on the intervals with endpoints \(u+h,u+h+1\) and \(u+1,u+2\)
can be used in the form
\[
 (u+1)^{h-1}\le\frac{\Gamma(u+h)}{\Gamma(u+1)}
          \le h^{-1}(u+1)^{h-1}.
\]
Specifically \(\Gamma(v+1-h)\le\Gamma(v)v^{1-h}\) at
\(v=u+h\) gives the first inequality, and
\(\Gamma(u+h+1)\le\Gamma(u+1)(u+1)^h\) gives the second by
recurrence. The products of these bounds give (3).

As a function of \(\xi_D\), \(w=\sqrt{D+2}\sin(\xi_D/\sqrt{D+2})/2\)
has derivative of magnitude at most \(1/2\) and \(w\ge1/2\).
Thus (3) bounds the ratio of the largest to the smallest pure
moment among nodes of diameter \(W\) by \(C_m^2(1+W)^{m-1}\).
Together with (2), for a mixed tuple at order \(m=2,3\) this gives
\[
 \int\prod_i f_{e_i}\ \le\
     21(1+W)^2 e^{-W^2/8}\min_i J_m(e_i).
 \tag{4}
\]
When every adjacent distance is at least \(12\), every mixed
tuple has \(W\ge12\). The factor on the right decreases there
and is at most \(3549e^{-18}<1/1024\).

Here is the stability estimate for the signed remainder. Label
the three unit-integral densities \(f_0,f_1,f_2\) in exponent order.
Write \(x,y,z\) for their pure quadratic moments, \(p,q,r\) for
their pure cubic moments, and \(\varepsilon=1/1024\).
The pure moments satisfy
\[
 y\le\max(x,z),\qquad
 x^2\le p\le2x^2,\quad y^2\le q\le2y^2,\quad z^2\le r\le2z^2.
 \tag{5}
\]
For the first inequality, duplication expresses \(J_2(u)\)
as a positive constant times
\(\Gamma(u+1/2)\Gamma(D-u+1/2)/(\Gamma(u+1)\Gamma(D-u+1))\).
Each gamma ratio as a function of its exponent is a beta moment
and is log convex by Hölder. The product is log convex in \(u\),
giving the endpoint bound. For the other inequalities, the lower
bound follows by Cauchy-Schwarz since \(\int f_i=1\).
Each \(f_i\) is log concave and attains its maximum at some mode
\(t_i\). The inequality \(f_i((t+t_i)/2)^2\ge f_i(t)f_i(t_i)\)
gives \(\int f_i^2\ge f_i(t_i)/2\) by a linear substitution.
This proves the upper bound on \(\int f_i^3\).

Use the moment coefficients \(A,B,C,D_3,F,G\) of
\(X(f_0-f_2)+Y(f_1-f_2)\), where \(D_3\) denotes the pure
\(X^3\) coefficient. Set \(A_0=x+z,\ C_0=y+z,\ B_0=2z\).
With all mixed moments set to zero, the negative remainder is
\[
 S_0=C_0^2p+2A_0zq+
        [xz+(3x+z)y-y^2]r .
\]
Define
\[
 E_p=C_0^2p,\qquad E_q=A_0zq,\qquad E_r=z(x+y)r,\qquad
 E=E_p+E_q+E_r .
\]
Then \(S_0\ge E/8\). If \(y\le z/2\), the cubic bracket is
at least \(xz+yz/2\). Otherwise the bracket is at least \(xz\)
by (5), and \(q\ge r/8,\ y\le A_0\) give
\(A_0zq\ge yz r/8\).

Bound the errors from the true mixed terms using (4).
In order two,
\[
 \begin{split}
 &0\le A_0-A\le2\varepsilon\min(x,z),\qquad
  0\le C_0-C\le2\varepsilon\min(y,z),\\
 &2(1-2\varepsilon)z\le B\le2(z+\varepsilon y),\qquad
 |B-B_0|\le2\varepsilon(\min(x,y)+\min(x,z)+\min(y,z)).
 \end{split}
\]
In order three write \(D_3=p-r+d_1,\ G=q-r+g_1,\ F=-3r+f_1\).
Then
\[
 |d_1|\le6\varepsilon p,\qquad
 |g_1|\le6\varepsilon\min(q,r),\qquad
 |f_1|\le3\varepsilon\bigl(3\min(p,q)+\min(p,r)+3\min(q,r)\bigr).
\]
The last bound follows by expanding
\((f_0-f_2)(f_1-f_2)^2\) and using the minimum participating
pure moment in each term.

For clarity, the losses from the pure and mixed contributions to
\[
 S=C^2p+ABq+(3AC-C^2-AB)r+C^2d_1+ABg_1-ACf_1
\]
are bounded by \(\varepsilon\) times the following quantities:
\[
\begin{array}{c|c}
 {\rm contribution} & {\rm loss}/\varepsilon\ {\rm at\ most}\\ \hline
 C^2p & 3E_p\\
 ABq & 6E_q\\
 (3AC-C^2-AB)r & 32E_r+4E_q\\
 C^2d_1 & 6E_p\\
 ABg_1 & 36E_q\\
 ACf_1 & 30E_p+27E_q
\end{array}
\]
For the third row, expand the three product differences. Use
\(C_0\min(x,z),\, A_0\min(y,z)\le z(x+y)\),
\(C_0\min(y,z)\le2yz\), \(A_0\min(x,z)\le2xz\),
and \(rA_0\min(x,y)\le E_r+2E_q\).
The last estimate follows from \(\min(x,y)\le y\) and the separate
cases \(y\le z\) and \(y>z\). In the second case \(yr\le2zq\).
The same cases give \(y\min(q,r)\le2zq\), which proves the fifth
row. For the last row, the three bounds are
\[
 A_0C_0\min(p,q)\le2E_p,\qquad
 A_0C_0\min(p,r)\le4E_p,\qquad
 A_0C_0\min(q,r)\le3E_q.
\]
For the first compare \(A_0/C_0\) with \(\max(1,x/y)\) and
use \(q/p\le2(y/x)^2\) when \(x>y\). For the second compare
\(A_0/C_0\) with \(1+x/z\) and use \(r/p\le2(z/x)^2\)
when \(x\ge z\). The third follows from \(C_0=y+z\) and
\(y\min(q,r)\le2zq\). The other table rows follow directly
from the quadratic perturbations and \(A\le A_0,\ C\le C_0\).

Thus \(S\ge E(1/8-73\varepsilon)>0\).
The angular chart here differs from the factorial chart by positive
column and row rescalings, which preserve the sign of the remainder.
Finally \(\xi_D'(u)=\sqrt{(D+2)/((u+1)(D-u+1))}\ge2/\sqrt{D+2}\),
proving the stated sufficient ordinary gap bound.

## 3. Effective separation for arbitrary fixed windows

Write \(M=n+r-1\), put
\[
 h(e)=\frac{\prod_\ell\sqrt{e_\ell+1}}{(D+1)^{k-1/2}},\qquad
 \tau_m(e)=h(e)^{m-1}J_m(e),\qquad
 \mu_m=(2\pi)^{(k-1)(1-m)/2}m^{-(k-1)/2}.
\]
For \(1\le m\le M\), the quantities \(\tau_m(e)\) have common
computable rational bounds \(0<\tau_-\le\tau_m(e)\le\tau_+\),
uniformly over all exponents. Moreover they tend uniformly to
\(\mu_m\) as the least exponent entry tends to infinity,
with an effective error bound. Here are direct estimates.
Set
\[
 U_m(u)=(u+1)^{(m-1)/2}
             \prod_{j=1}^{m-1}\frac{\Gamma(u+j/m)}{\Gamma(u+1)}.
\]
The gamma bounds in Section 2 give \(1\le U_m(u)\le C_m\).
Keeping the factor \((u+1)/(u+j/m)\) in that upper bound
gives \(U_m(u)\le(1+1/H)^{m-1}\) for \(u\ge H\ge1\).
Multiplication and recurrence for Gamma give
\[
 \frac{\tau_m(e)}{\mu_m}
 =\frac{\prod_\ell U_m(e_\ell)}{U_m(D)}
   \prod_{j=1}^{k-1}
    \left(\frac{D+j}{D+1}\right)^{m-1}
    \frac{m(D+j)}{mD+j}.
 \tag{6}
\]
For each \(j\), the last pair of factors is between \(1\)
and \(m j^{m-1}\). This gives rational enclosures for
\(\tau_m\) after enclosing \(\mu_m\). If all entries are
at least \(H\), the last pair is at most \((1+k/H)^m\).
In that case (6) puts \(\tau_m/\mu_m\) between
\((1+1/H)^{-(m-1)}\) and
\((1+k/H)^{m(k-1)+k(m-1)}\).

We use two finite-dimensional margins. Fix \(C\ge1\), a block size
\(s\le r\), and \(s\) consecutive orders in \([n,M]\).
For complex \(z_i\) with \(\max_i|z_i|=1\) and real
\(\lambda_i\in[C^{-1},C]\), the forms
\(\mu_m\sum_i\lambda_i z_i^m\) cannot vanish simultaneously.
Group the equal nonzero \(z_i\). Their combined weights are
positive, and the resulting weighted Vandermonde system excludes a
common zero. If \(s\le2\), the forms
\(\sum_i\lambda_i \tau_{mi}z_i^m\) with arbitrary
\(\tau_{mi}\in[\tau_-,\tau_+]\) also cannot vanish simultaneously.
For two nonzero \(z_i\), dividing the two consecutive
equations forces \(z_1/z_2\) to be positive real, incompatible
with either equation. A single nonzero \(z_i\) is immediate.

By compactness, in either case some row has absolute value at least
\(2\delta(C)>0\), uniformly over the finite list of possible blocks.
A positive rational \(\delta(C)\) is computable by real quantifier
elimination on the coefficient and weight domains. In the first
case divide out \(\mu_m\) in each row and then use a rational
lower bound for the \(\mu_m\). The normalization
\(\max_i|z_i|=1\), the weight intervals, and the bounds on the
absolute values of the forms are real polynomial conditions.
Choose \(H(C)\) so large that
\[
 |\tau_m(e)-\mu_m|<\delta(C)/(rC)
 \quad\text{whenever all entries of }e\text{ are at least }H(C).
 \tag{7}
\]
Equation (6) gives an explicit procedure for choosing \(H(C)\).

Choose \(r+1\) increasing levels \(K_0,\ldots,K_r\) beginning with
\(K_0=2\). At level \(j\) put \(C_j=K_j^r\) and
\(\delta_j=\delta(C_j)\). Pick \(R_j>1\) with
\(r C_j\tau_+/R_j<\delta_j/4\) and then take
\(K_{j+1}>\max(K_j,R_j^{2M+1})\).
Among the \(r\) intervals \([K_j,K_{j+1})\), at least one contains
none of the \(r-1\) adjacent width ratios when the widths are
ordered decreasingly. Use this level \(j\). Split the ordered
columns whenever that ratio is at least \(K_{j+1}\).
Inside a group all width ratios are bounded by \(C_j\);
between groups the ratio of any earlier to any later width
exceeds \(R_j^{2M+1}\).

For \(k>2\) choose the coordinate cutoff at least
\(\max_j H(C_j)\). For \(k=2\) a group containing a column whose
least exponent is less than \(H(C_j)\) has at most two columns,
provided the angular separation is large enough. Indeed, the ratio
of widths \(h(e)\) in this dimension agrees with the ratio of
the widths \(w(u)\) in Section 2, and
\((\min(u,D-u)+1)/2\le w(u)^2\le\min(u,D-u)+1\).
Thus every column of that group has least exponent plus one
at most \(2C_j^2(H(C_j)+1)\).
A column with such a bound lies at angular distance at most
\(\pi\sqrt{2C_j^2(H(C_j)+1)}\) from its nearer axis, by
\(\arcsin v\le\pi v/2\). Choose \(T\) larger than twice this
quantity for every level. Then at most one column of the group
lies near each axis. This gives the small group case even without
a coordinate cutoff.

Suppose a common projective zero in the window exists. In the
unit-integral angular chart let its coefficients be \(c_i\), and
put \(v_i=c_i/h_i\). Give each group \(g\) a representative
width \(h_g\) and put \(N_g=\max_{i\in g}|v_i|\).
Its row weight is \(W_g(m)=h_g N_g^m\).
If a later group \(b\) is within a factor \(R_j\) of the maximum
at row \(m\), then for any earlier nonzero group \(a\),
\[
 (N_b/N_a)^m\ge(h_a/h_b)/R_j>R_j^{2M}.
\]
Consequently \(W_b(m+1)/W_a(m+1)>R_j\). The
[ordered assignment lemma](SFC_SEPARATED_DEGREE_CLUSTER_WINDOWS.md#1-the-ordered-assignment-of-moment-rows)
now supplies one group uniquely active throughout its assigned
consecutive window of length equal to its number of columns.
At these rows its weight is the maximum and every other group
has weight below this maximum divided by \(R_j\).

Divide the assigned rows by the chosen group weight. Its pure
terms have \(\lambda_i=h_i/h_g\) and \(z_i=v_i/N_g\).
If all entries in the group are at least \(H(C_j)\),
(7) gives a row with pure contribution of modulus at least
\(\delta_j\). Otherwise it is the bivariate group of size at
most two, and its margin is at least \(2\delta_j\).
The other pure terms have total modulus below
\(r C_j\tau_+/R_j<\delta_j/4\).
By (2) every ordered mixed tuple contributes after normalization
at most \(e^{-T^2/8}C_j\tau_+\). Indeed, before the
exponential factor its absolute value is bounded by the
geometric mean of the pure bounds
\(\tau_+ h_i|v_i|^m\le \tau_+ C_j W_g(m)\).
Choose \(T\) large enough that
\(r^M e^{-T^2/8} C_j\tau_+<\delta_j/4\) at every level.
In the row with margin \(\delta_j\), the full moment then
cannot vanish.

All choices of \(T\) above admit rational upper bounds. The
factorial moment forms at integer exponents have rational
coefficients. As usual, absence of a common complex zero on
each projective coefficient chart transfers by the Nullstellensatz
to every characteristic-zero field.

The sign theorem controls actual real heights in a separated regime.
The unrestricted sign and the coefficientwise normalized gap
increments remain open.
