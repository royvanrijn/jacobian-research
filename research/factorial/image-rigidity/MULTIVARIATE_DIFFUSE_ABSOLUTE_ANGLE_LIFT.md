# A degree-independent diffuse tolerance for ordinary factorial angular separation

For distinct nonnegative integer exponent columns, delete equal rows
and put \(S=\sum_i(b_i-a_i)^2/\max(a_i,b_i)\) and
\(\eta=\max_i|b_i-a_i|/\max(a_i,b_i)\).
The ordinary factorial moment polynomials at any two consecutive
orders have disjoint root arguments under independent positive
rescalings whenever
\[
                         S\ge13,\qquad S\eta\le10^{-9}.      \tag{1}
\]
The two moments are coprime over every characteristic-zero field.
The tolerance is independent of the order.

## General analytic lift

Fix the Gaussian parameter \(S\ge13\), \(n\ge2\) and \(N=n+1\).
Suppose the two positive polynomial arrays are obtained from
\(B_j(z;S)=\sum_k\binom jk e^{-S k(j-k)/(2j)}z^k\) by multiplying
each coefficient at \(z^k\) by \(e^{E_j(k/j)}\).
Assume \(E_j\) is real-symmetric and analytic on \(|v|\le5\),
\(E_j(0)=E_j(1)=0\), and
\[
 \|E_j''\|_{|v|\le5}\le10jU,\qquad
 \|E_N''-E_n''\|_{|v|\le5}\le13U,\qquad U\le10^{-9}.
                                                               \tag{2}
\]
Then the two arrays have disjoint root arguments. The proof compares
the two analytic halves of each moment separately, as in the
[polynomial diffuse lift](MULTIVARIATE_DIFFUSE_POLYNOMIAL_ANGLE_LIFT.md).
Here the constant and linear parts of the multiplier's logarithm are
applied exactly around an implicit fractional index.

For \(n\le21\), that lift applies directly to (2).
The proof there only needs \(\|E_j''\|\le10jU\) and the two
endpoint values, which give \(\|E_j\|\le150jU\), and \(NU\le10^{-6}\)
in its half-log comparison. Use \(n\ge22\) below.
For either choice \(\Phi_j(v)=E_j(v)\) or \(E_j(1-v)\), (2) implies
\[
 \Phi_j(0)=0,\quad |\Phi_j''(v)|\le10jU\ (|v|\le4),\quad
 |\Phi_j'(v)|\le50jU\ (|v|\le1/50),\quad
 |\Phi_N'(v)-\Phi_n'(v)|\le72U\ (|v|\le4).
\]
The derivative bounds follow by integrating twice with the two
endpoint values. For example
\(\|E_N'-E_n'\|_{|v|\le5}\le(1/2+5)13U<72U\).

## Centering the half-multiplier

Use \(I_j=e^{jF_*}(1+\omega_j)\) from the
[Gaussian half-integral proof](MULTIVARIATE_DIFFUSE_EXPLICIT_S13_INTERLACING.md).
It has nonnegative coefficients, with \(|F_*|<1/250\)
and \(|\omega_j|<21/100\) on \(|z|\le2\).
Here \(F_*\) is common to the two degrees at fixed \(S\).
Write \(\sigma(z)=zF_*'(z)\) and \(D=z\,d/dz\).
Cauchy's bounds give \(|\sigma(z)|\le3/250\) and
\(|z\sigma'(z)|<1/10\) on \(|z|\le3/2\).

For \(|z|\le7/5\) define \(w_j=w_j(z)\) by
\[
              w_j=\sigma(z e^{h_j}),\qquad h_j=\Phi_j'(w_j)/j.
\]
On \(|w|\le1/50\) the right side has modulus at most \(3/250\),
since the shifted argument is less than \(3/2\) in modulus.
Its derivative in \(w\) has modulus at most \(U\), by (2).
The fixed point is unique, analytic in \(z\), and real-symmetric.
Put \(z_j=z e^{h_j}\) and let \({\cal R}_j(v)\) be
\(\Phi_j(v)-\Phi_j(w_j)-\Phi_j'(w_j)(v-w_j)\).
For each \(z\), hold \(w_j,h_j\) fixed when applying
\[
 \widetilde I_j(z)
 =e^{\Phi_j(w_j)-w_j\Phi_j'(w_j)}
       \left[e^{{\cal R}_j(D/j)}I_j\right](z_j).             \tag{3}
\]
This identity follows termwise, since the combined action
multiplies each coefficient at \(z^k\) by \(e^{\Phi_j(k/j)}\).
The two choices of \(\Phi_j\) give the two actual halves.

For \(|Z|\le3/2\), hold \(w=\sigma(Z)\) constant and write
\(D_{j,m}(Z)=[(D/j-w)^m I_j](Z)/I_j(Z)\).
For \(1\le m\le j\), apply Cauchy in \(h\) to
\(e^{-jwh}I_j(Ze^h)/I_j(Z)\) on \(|h|\le\sqrt{m/j}/4\).
The first derivative of \(F_*(Ze^h)-F_*(Z)-wh\) is zero
at \(h=0\). Its second derivative in \(h\) has modulus
less than \(13\): on \(|Ze^h|<39/20\), Cauchy gives
\(|F_*'|\le2/25\) and \(|F_*''|\le16/5\).
The quotient in the Cauchy estimate is bounded by
\(2e^{13m/32}\). With \(m!\le m^m\) this gives
\[
                   |D_{j,m}(Z)|\le2(7\sqrt{m/j})^m .
\]
For \(m>j/64\), coefficient positivity instead gives
\(|D_{j,m}(Z)|\le (1+|w|)^m I_j(|Z|)/|I_j(Z)|\le2(9/4)^m\).
Here \(I_j(|Z|)/|I_j(Z)|\le2e^{j/125}\), \(j<64m\),
and \((51/50)e^{64/125}<9/4\).

Set \(B=5jU\) and let \(c_m\) be the coefficients in the series
of \(e^{{\cal R}_j(v)}-1\) about \(w_j\).
It starts at \(m=2\), and \(|{\cal R}_j(v)|\le B|v-w_j|^2\)
on \(|v-w_j|\le3\). To bound the relative remainder in (3):

- For \(2\le m\le j/64\) with \(m\le9B\), Cauchy at radius
  \(\sqrt{m/B}\) gives \(|c_m|\le e^m(B/m)^{m/2}\).
  Hence \(|c_m D_{j,m}(z_j)|\le2(47\sqrt U)^m\)
  and this sum is less than \(5000U\).
- For \(2\le m\le j/64\) with \(m>9B\), Cauchy at radius
  three and \(e^{9B}\le e^m<3^m\) give \(|c_m|\le9B\).
  The sum is at most
  \(18B(49/j)\sum_{m\ge2}m(7/8)^{m-2}=317520U\).
- For \(m>j/64\), the same circle gives
  \(|c_m D_{j,m}(z_j)|\le18B e^{9B}(3/4)^m\).
  Since \(9B\le2880Um\), this is at most \(18B(4/5)^m\).
  The sum is at most \(90B(4/5)^{j/64}\le144000U\),
  using \(j e^{-j/320}\le320\).

The relative error is less than \(470000U\) and its logarithm
has modulus at most \(10^6 U\). Moving \(\log(1+\omega_j)\)
from \(z_j\) back to \(z\) costs less than \(38U\), by
Cauchy on \(|z|\le3/2\) and \(|h_j|\le50U\).
Thus the actual half has analytic logarithm on \(|z|\le7/5\)
\[
 H_j+W_j+e_j,\qquad
 W_j=j(F_*(z_j)-F_*(z))+\Phi_j(w_j)-w_j\Phi_j'(w_j),\qquad
 |e_j|\le\varepsilon:=1100000U,                              \tag{4}
\]
with \(H_j=\log I_j\). All these functions have conjugation
symmetry and are real on the positive radius.


The leading correction satisfies
\[
                  |W_j|\le3jU,\qquad |W_N-W_n|\le7U.        \tag{5}
\]
The first inequality uses \(|w_j|\le1/50\),
\(|\Phi_j'(w_j)|\le50jU\), and
\(|F_*(z_j)-F_*(z)|\le(3/250)|h_j|\le3U/5\).
For the order comparison, at a common \(w\) the log dilations
\(\Phi_j'(w)/j\) differ by at most \(122U/N\).
The derivative of \(\sigma\) with respect to logarithmic argument
is less than \(1/10\). Combining this with the contraction
gives \(|w_N-w_n|\le26U/N\) and \(|h_N-h_n|\le123U/N\).
The first term of \(W_N-W_n\) is bounded by
\((3/250)(50+123)U<3U\).
The other two terms are bounded respectively by
\((72/50)U+1300U^2<2U\) and
\((72/50)U+(1300+26/5)U^2<2U\), using the symbol derivative bounds.
This proves (5) for each orientation of the half.

## Coupled radial and angular comparison

Use the \(+\) and \(-\) half logarithms (4) given by
\(\Phi_j=E_j\) and \(E_j(1-\cdot)\), respectively.
The bounds \(\|H_j\|_{|z|\le2}\le j/18\) and (4)–(5)
show that the actual roots are in \(|\log|z||<1/8\),
by the outer and reciprocal Rouché argument of the
polynomial diffuse lift. Indeed
\(2(1/18+3U+\varepsilon/j)<1/8\).

Let \(J_j(q,\theta)\) be the actual plus half logarithm
at \(e^{q+i\theta}\) minus the minus half logarithm at
\(e^{-q-i\theta}\). Let \(J_j^0\) use \(H_j\) at both
points. For a root with \(z=e^{q_j+i\theta}\) and
\(0<\theta\le\pi\), the equation is
\[
 j(q_j+i\theta)=J_j(q_j,\theta)+(2k_j+1)i\pi .
\]
The evaluation points have modulus below \(6/5\).
The baseline derivative bounds from the polynomial lift are
\(|\partial_qJ_j^0|\le j/6\) and
\(|\partial_q\partial_\theta J_j^0|\le2j/3\).
For either half correction bounded by \(M\) on \(|z|\le7/5\),
Cauchy at radius difference \(1/5\) bounds \(z\) times its
first derivative by \(6M\), and \(z\) times its first derivative
plus \(z^2\) times its second derivative by \(78M\).
The two corrections together contribute at most \(12M\)
and \(156M\) to those respective derivatives of \(J_j\).
With \(M=3jU+\varepsilon\), this gives
\(|\partial_q J_j|\le j/5\) and
\(|\partial_q\partial_\theta J_j|\le j\).

Since \(\Re J_j^0(0,\theta)=0\), the real modulus equation
also implies \(|q_j|\le(12/5)(3U+\varepsilon/j)\).
Suppose the two roots have a common argument \(\theta\).
At a common \(q\), the real part of the baseline difference
vanishes at \(q=0\) and its \(q\)-derivative is at most \(3/2\)
because \(\|H_N-H_n\|_{|z|\le2}<12/25\).
Using (5), the real parts at the actual radii give
\[
 N|q_N-q_n|
 \le \tfrac54\bigl(\tfrac52|q_n|+14U+4\varepsilon\bigr)
 \le 40U+6\varepsilon .
\]
The first inequality uses \(|\partial_q J_N|\le N/5\)
to absorb the radial difference on the right-hand side.

Compare the imaginary order difference to the baseline at \(q=0\).
By symmetry the imaginary parts vanish at \(\theta=0\) for
real \(q\). The mixed derivative at order \(N\) bounds the
error of changing \(q_N\) to \(q_n\) by
\((40U+6\varepsilon)\theta\).
The baseline order difference's mixed derivative is at most
\(12\|H_N-H_n\|_{|z|\le2}<6\), so changing \(q_n\) to zero
there costs at most \((44U+\varepsilon)\theta\).
At a common radius the two half differences of the leading
corrections (5) contribute at most \(84U\theta\), and the
two differences of \(e_j\) contribute at most \(24\varepsilon\theta\).
Thus the total imaginary error is at most
\((168U+31\varepsilon)\theta<7\theta/200\).

At \(q=0\) the baseline imaginary order difference is
\(2\Im(H_N-H_n)(e^{i\theta})\).
Its modulus is at most \(24\theta/25\) because the derivative
norm on \(|z|\le1\) is less than \(12/25\) and both functions
are real at \(z=1\). Subtracting the two root equations now places
an integer multiple of \(2\pi\) strictly between
\(\theta/200\) and \(399\theta/200\), with \(0<\theta\le\pi\).
This is impossible, proving the general analytic lift.

## Verification for ordinary factorial coefficients

Under (1) every unequal row is shared. Take the Gamma interpolant
\(E_j\) of the polynomial diffuse lift for real \(j\in[n,N]\);
it has the required exact sampled values and endpoint zeros.
That proof gives \(\|E_j''\|_{|v|\le5}\le10jS\eta\).
It remains to verify the adjacent bound in (2).

Write \(x_i(v)=a_i+(b_i-a_i)v\), \(L_i=\max(a_i,b_i)\), and
\(q=jx_i(v)\). The same sum–integral comparison as for the
trigamma gives
\[
 |\psi''(q+1)+1/q^2|\le2/(\Re q)^3 .
\]
Indeed the difference between the sum and integral for
\(-2(q+t)^{-3}\) is bounded by the integral of
\(6/(\Re q+t)^4\). Differentiating \(E_j''(v)\) in \(j\)
gives row terms proportional to
\(2j\psi'(q+1)+j^2 x_i(v)\psi''(q+1)\), whose integral
approximation is \(1/x_i(v)\).
Use \(|x_i(v)-L_i|\le6|b_i-a_i|\),
\(\Re x_i(v)\ge(1-6\eta)L_i\), and \(1/L_i\le\eta\)
for the unequal integer rows.
On \(|v|\le5\) the sum, after subtracting \(S\), is bounded by
\[
 U\left(\frac6{1-6\eta}+\frac2{j(1-6\eta)^2}
                  +\frac{2(1+6\eta)}{j(1-6\eta)^3}\right)<13U .
\]
Integrating between \(n\) and \(N\) supplies (2).
Endpoint normalization preserves arguments, so the general
lift proves the result for \(n\ge2\). At \(n=1\), the first
moment has only a negative real root and the second moment
is positive on the real line by the integral of a nonzero
real binomial squared. Rational Bézout gives the stated
characteristic-zero conclusion.

