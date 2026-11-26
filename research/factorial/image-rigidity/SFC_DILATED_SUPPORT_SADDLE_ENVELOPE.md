# All-order Strong Factorial after dilation with distinct weighted degrees

Let \(b_1,\ldots,b_s\) be fixed nonnegative integer exponent columns in
\(r\) variables. Write \(d_i=\sum_\nu b_{\nu i}\), and assume that the
\(d_i\) are pairwise distinct. For a positive integer \(D\) put
\[
                         F_D(X)=\sum_i c_i X^{D b_i}.
\]
Let \(\Lambda(X^\alpha)=\prod_\nu \alpha_\nu!\).

**Theorem.** There is a computable threshold in the exponent columns
alone such that for every \(D\) above it and every nonzero coefficient
vector, the number of positive integers \(m\) at which
\(\Lambda(F_D^m)=0\) is at most the number of its nonzero coefficients
minus one. The field may be any characteristic-zero field. The orders
need not be consecutive.

Every fixed ordinary support admits such a dilation in the same
variables, with the same number of monomials and the same zero pattern
in its exponent matrix. For example, choose \(B\) larger than every
entry and multiply row \(\nu\) by \(B^{\nu-1}\) before uniformly
multiplying all rows by \(D\). Its column degrees are distinct by
base-\(B\) uniqueness. One can choose \(B,D\) to be powers of the same
prescribed prime. The [prime-power row comparison](SFC_PADIC_ROW_INFLATION_HIGH_PRECISION.md)
then also applies at that prime under its hypotheses. The theorem
does not assert removal of the dilation.

The proof separates the monomial peaks of the ordinary Gamma
integral. A peak contains mixed terms whose total contribution can
be large at high orders. Its complex saddle supplies a nonzero local
amplitude and an effective exponential rate. The rate functions have
strictly ordered logarithmic slopes. Each zero moment then requires
a distinct switch in the maximizing peak. The proof includes peaks
on coordinate faces and allows the constant monomial.

## Exact integral and separation of the peaks

Work over \(\mathbb C\), discard zero coefficients, and assume at
least two remain. The case of one monomial is immediate. All
constants below can be chosen from the original columns, so they
apply to any of the finitely many active subsets. For \(t>0\) define
\[
 \begin{split}
 u_i(y)&=e^{-\sum_\nu y_\nu}\prod_\nu y_\nu^{b_{\nu i}},&
 P_i(y)&=u_i(y)/u_i(b_i),\\
 q_i(t)&=c_i(Dt)^{D d_i}u_i(b_i)^D,&
 Q(t)&=\max_i |q_i(t)|.
 \end{split}
 \tag{1}
\]
Here \(0^0=1\) in the evaluation of \(u_i(b_i)\).
The functions \(P_i\) are entire in \(y\), nonnegative on the
positive orthant, and each has its unique maximum \(1\) there at
\(b_i\), including its boundary. This follows by maximizing
\(y^a e^{-y}\) for each nonnegative integer \(a\). For every
positive integer \(m\), set \(\lambda=mD\). Scaling \(X_\nu=\lambda
y_\nu\) in the product-exponential integral gives the exact formula
\[
 \Lambda(F_D^m)=\lambda^r
       \int_{\mathbb R_{\ge0}^r}\left(\sum_i q_i(m)P_i(y)^D\right)^m\,dy .
 \tag{2}
\]

Take a small positive rational \(\gamma\) and small disjoint real
boxes \(V_i\) about the peaks, restricted to the orthant. In
a coordinate where \(b_{\nu i}=0\) the interval is \([0,\epsilon_i]\),
and in the other coordinates it is
\([b_{\nu i}-\epsilon_i,b_{\nu i}+\epsilon_i]\). Choose
\(\epsilon_i\) so that \(P_i\) is nonvanishing on a complex
polydisc centered at \(b_i\) containing \(V_i\), and, on that
polydisc,
\[
                  |P_j/P_i|\le e^{-4\gamma}\qquad(j\ne i).
 \tag{3}
\]
Use a smaller concentric polydisc still containing \(V_i\) for
derivative bounds. These choices exist since \(P_j(b_i)<1\)
for \(i\ne j\). Off the union of the boxes there is an absolute
\(\delta>0\) with \(P_j\le e^{-\delta}\) for every \(j\).
Indeed, lying outside a box around \(b_j\) forces one
factor of \(P_j\) to be bounded below its maximum by a fixed
factor. In each coordinate that bound can be evaluated at the
endpoints of the corresponding interval by monotonicity.

At order \(t\) call \(i\) eligible if
\(|q_i(t)|\ge e^{-\gamma D}Q(t)\). Its set of eligible times is
an interval in \(\log t\), since each \(\log|q_i(t)/q_j(t)|\) is
affine in \(\log t\). At an eligible peak write the local phase as
\[
 \Phi_i(y,t)=\log P_i(y)+\frac1D\log
 \left(1+\sum_{j\ne i}\frac{q_j(t)}{q_i(t)}
                         \big(P_j(y)/P_i(y)\big)^D\right).
 \tag{4}
\]
Use the analytic logarithm of \(P_i\) with value zero at \(b_i\),
and the power-series logarithm for the perturbation. By (3)
the perturbation is \(O(D^{-1} e^{-3\gamma D})\), uniformly on the
polydisc and for all eligible parameters when \(D\) is large.
Cauchy estimates give the same bound for any fixed number of
\(y\)-derivatives on the smaller polydisc. For the partial derivative
in \(\log t\) the bound is \(O(e^{-3\gamma D})\), since each ratio
has logarithmic derivative \(D(d_j-d_i)\) times itself.
The same constructions extend to a neighborhood of each
finite eligibility endpoint, using \(|q_j/q_i|\le e^{2\gamma D}\)
there.

For noneligible \(i\), the integral on \(V_i\) has absolute value at
most \(\operatorname{vol}(V_i)Q(m)^m(s e^{-\gamma D})^m\).
This uses \(|q_i|<e^{-\gamma D}Q\) and (3) to bound all other terms
there. Off the union, the integrand is bounded by
\[
      (sQ(m))^m e^{-\lambda\delta/2}
                           \max_j P_j(y)^{\lambda/2}.
\]
Each \(P_j\) is integrable on the orthant, so for \(\lambda\ge2\)
the integral of the last maximum is bounded by
\(\sum_j\int P_j\). Increasing the threshold in \(D\), the
combined omitted contribution to the integral in (2) is
\(O(Q(m)^m e^{-c\lambda})\) for a fixed \(c>0\), uniformly in all \(m\).

## Uniform local saddle and endpoint estimate

Write \(A_i=\{\nu:b_{\nu i}>0\}\) and \(E_i=\{\nu:b_{\nu i}=0\}\).
Keep the coordinates in \(E_i\) at zero, and let
\((\theta_i(t),0)\) denote the solution of
\(\nabla_{A_i}\Phi_i=0\) near \(b_i\). For the pure phase
\(\log P_i\) the active Hessian at the center is
\(-\operatorname{diag}(1/b_{\nu i})_{\nu\in A_i}\).
For example, adding \(\operatorname{diag}(b_{\nu i})\,
\nabla_{A_i}\Phi_i\) to the active coordinate map gives a map
whose derivative vanishes at the pure center. On a fixed sufficiently
small active polydisc the derivative norm is at most \(1/2\) for
large \(D\); its displacement at the center is
\(O(D^{-1}e^{-3\gamma D})\). It is a contraction into that polydisc.
Thus the saddle exists, has this same displacement bound, and
depends differentiably on \(t\) on its interval by the implicit
function theorem. If \(A_i\) is empty, all these assertions mean
evaluation at the endpoint \(y=0\).

Let \(H_i\) be the active Hessian and \(v_{i,\nu}\) the \(E_i\)-partial
derivatives of the phase at the saddle. Put
\[
 \rho_i=|A_i|/2+|E_i|,\qquad
 K_i(t)=(2\pi)^{|A_i|/2}\det(-H_i)^{-1/2}
                       \prod_{\nu\in E_i}(-v_{i,\nu})^{-1}.
 \tag{5}
\]
Use the square-root branch continuing the pure positive
determinant, and the empty determinant and product conventions.
These amplitudes are bounded above and away from zero uniformly.
Indeed \(H_i\) is uniformly close to the pure diagonal Hessian
and \(v_{i,\nu}\) is close to \(-1\).

The following estimate at integer \(t=m\) holds uniformly at each
eligible peak with \(\lambda=mD\):
\[
 \int_{V_i}\left(\sum_j q_j(m)P_j(y)^D\right)^m dy
 =q_i(m)^m e^{\lambda\Phi_i(\theta_i(m),0;m)}
     \lambda^{-\rho_i}\big(K_i(m)+O(\lambda^{-1/2})\big).
 \tag{6}
\]
Here are details of the uniform estimate. In the integral of
\(\exp(\lambda\Phi_i)\), translate the active intervals by
\(\theta_i-b_{A_i,i}\), retaining the endpoint intervals.
Cauchy's theorem applied one coordinate at a time gives side
integrals whose absolute values are \(O(e^{-c_1\lambda})\).
In a side integral at least one active coordinate is within a
small complex displacement of a box edge. The pure phase has a
strict negative real bound there, uniform in the other variables
on their untranslated real intervals. The small translations
and the bound on the perturbation preserve a fixed negative
bound; side contour lengths and remaining areas are uniformly
bounded.

On the translated box write the active displacement as
\(x\in\mathbb R^{|A_i|}\) and the inactive variables as
\(y\in\mathbb R_{\ge0}^{|E_i|}\). After subtracting the saddle
value the phase has real part at most
\(-c_2(|x|^2+\sum y_\nu)\), with \(c_2>0\) uniform. To see
this, along the active real directions in the translated box at
\(y=0\) the real Hessian is uniformly negative definite, and
the active gradient at its center is zero. With active
coordinates fixed, the real part of each inactive partial
derivative is at most \(-1/2\) throughout the box. Both estimates
hold by the pure phase formulas, the small translation and the
uniform Cauchy bounds. The quadratic-linear approximation is
\(\frac12 x^{\mathsf T}H_i x+\sum v_{i,\nu}y_\nu\). It has
the same real bounds with a possibly smaller \(c_2\).
The Taylor error is
\[
                 O\left(|x|^3+|x|\sum y_\nu+(\sum y_\nu)^2\right).
\]
For the difference of the two exponentials use the integral of
the derivative along the straight segment between the two
phases. Its absolute value is at most \(\lambda\) times this
Taylor error times \(e^{-c_2\lambda(|x|^2+\sum y_\nu)}\).
Scaling \(x\) by \(\lambda^{-1/2}\) and \(y\) by
\(\lambda^{-1}\) makes the integrated error
\(O(\lambda^{-\rho_i-1/2})\). Extending the quadratic-linear
integral to the full active real space and inactive orthant has
exponentially small error. The resulting Gaussian integral
continues holomorphically from real definite matrices, and each
inactive integral is that of an exponential with strictly
negative real rate. They give (5). The contour-side errors
remain exponentially small relative to the main scaling, since
\(\Phi_i(\theta_i,0;t)=O(D^{-1}e^{-3\gamma D})\).
This proves (6), including the case that every coordinate
is inactive.

## Counting switches at all orders

For \(t\ge1/2\) and eligible \(i\), define the real function
\[
 G_i(t)=\log |q_i(t)|+
         D\operatorname{Re}\Phi_i(\theta_i(t),0;t)
                    -\frac{\rho_i}{t}\log(tD).
 \tag{7}
\]
The magnitudes in (6) are \(e^{mG_i(m)}\) times amplitudes
bounded above and away from zero when \(D\) is large. The
following estimates are uniform for eligible \(i\) and \(t\ge1/2\):
\[
 \begin{split}
     G_i(t)-\log |q_i(t)|&=O(\log(2D)),\\
 \frac{dG_i(t)}{d\log t}&=D d_i+O\big(D e^{-3\gamma D}+\log(2D)\big).
 \end{split}
 \tag{8}
\]
The derivative of the critical value uses only the partial
\(\log t\)-derivative of \(\Phi_i\), because the active gradient
vanishes at the saddle and the boundary coordinates stay zero.
The last term of (7) and its derivative are bounded as in (8)
because \((1+|\log t|)/t\) is bounded on \([1/2,\infty)\).

The maximum of the eligible functions \(G_i\) is continuous.
When a function enters or leaves eligibility, its raw
\(\log|q_i|\) is below the raw maximum by \(\gamma D\), which
dominates the adjustment \(O(\log(2D))\). It cannot take part in
a switch there. In the intersection of their eligible intervals
every pairwise difference has derivative in \(\log t\) with
the sign of \(d_i-d_j\) and magnitude at least \(D/2\), by
(8) and the distinct integral degrees. Thus each pair crosses
at most once. At every switch the active slope index increases
in the ordering by \(d_i\). There are at most \(s'-1\) such
switches, where \(s'\) is the number of active coefficients.

Suppose an integer order \(m\) gives a zero of (2). Then at least
two eligible functions exist, and a runner-up differs from the
maximum by at most \(C_0/m\) for a constant \(C_0\) independent of
\(D,m\) and the coefficients. Otherwise the largest eligible local
integral strictly dominates the sum of the others and the
omitted integral. For the eligible integrals this follows
from (6) with uniform upper and lower amplitude bounds. For
the omitted integral, the eligible maximum is at least the
adjusted value of a raw maximizing \(q_i\), and its exponential
\(e^{m\max G_i(m)}\) is at least
\[
              Q(m)^m e^{-C m e^{-3\gamma D}}
                                   (mD)^{-\max_i\rho_i}.
\]
The omitted bound above is \(O(Q(m)^m e^{-c mD})\), hence
uniformly small relative to this lower bound for \(m\ge1\)
as \(D\) grows.

Take a fixed \(h\in(0,1/3]\) so small that
\(2h\max_{i,j}|d_i-d_j|\le\gamma/4\). At a zero order the
maximal and runner-up raw logarithms differ from \(\log Q(m)\)
by at most \(O(\log(2D))+C_0\). For \(D\) sufficiently large
this is less than \(\gamma D/4\). For \(t\in[m-h,m+h]\)
their deficits can increase by at most
\[
      D\max_{i,j}|d_i-d_j|\,|\log(t/m)|\le \gamma D/4 .
\]
They thus stay strictly eligible. Their difference at \(m\)
has magnitude at most \(C_0/m\), while its monotone change
on either side, in the required direction, has magnitude at
least \(D h/(4m)\) before the endpoint of the interval.
For large \(D\) the runner-up overtakes the maximizing
function in one of the two directions. The maximum has a
switch in \((m-h,m+h)\). This also follows directly when the
two are tied at \(m\), since their derivatives differ and both
are eligible near \(m\). The intervals for distinct integers
are disjoint. Thus there are at most \(s'-1\) zero orders,
as claimed.

All threshold choices can be made effectively. The boxes and
the positive margins in (3) are tested by elementary
exponential and monomial bounds on small complex polydiscs.
The pure Hessian eigenvalues, inactive slopes and all the
pure derivative bounds on the prescribed discs have explicit
bounds in the exponents. Cauchy estimates, the contraction
test, and integrals of Gaussian and exponential polynomials
bound the uniform constants in (6) and (8). The remaining
conditions require only fixed inequalities between \(D\),
\(\log(2D)\), and decaying exponentials. No numerical
threshold has been evaluated. The field extension follows by
embedding the field generated by the finitely many
coefficients into \(\mathbb C\).

The ordinary exponential integral is the gamma identity on
monomials. The local argument above is the active-saddle and
linear-endpoint version of the methods compared in
[NIST DLMF, §§2.4(i), 2.4(iii), 2.4(iv)](https://dlmf.nist.gov/2.4).
Its uniformity in the coefficients comes from (3) and the
eligibility restriction. [Rocks, Theorem 3.10](https://people.math.sc.edu/rocks/BradyRocksSFC2015.pdf)
gives all-order results for sums of private prime-power
monomials by a different arithmetic argument. The hypotheses
here allow arbitrary exponent sharing, including fully shared
supports. They use a sufficiently large dilation of the given
columns and distinct weighted total degrees.
