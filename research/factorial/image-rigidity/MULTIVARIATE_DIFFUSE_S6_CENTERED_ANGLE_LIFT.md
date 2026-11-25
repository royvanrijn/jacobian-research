# Exact factorial angular separation down to Gaussian coupling six

For distinct ordinary exponent columns, omit equal rows, let
\(d_i=b_i-a_i\), \(L_i=\max(a_i,b_i)\),
\(S=\sum_i d_i^2/L_i\), and \(\eta=\max_i|d_i|/L_i\).

**Theorem.** The ordinary factorial moment polynomials of orders
\(n,N=n+1\), \(n\ge10000\), have disjoint root arguments whenever
\[
                   6\le S\le8,\qquad S\eta\le4\cdot10^{-15}.
                                                               \tag{1}
\]
The result also holds under the signed matched-tanh conditions
in (5) below. Positive radial rescalings are independent at the
two orders. Consequently the original integer polynomials are
coprime over every characteristic-zero field.

The [Gaussian saddle on the smaller disk](MULTIVARIATE_DIFFUSE_EXPLICIT_S6_INTERLACING.md)
allows the [centered coefficient lift](MULTIVARIATE_DIFFUSE_S8_CENTERED_ANGLE_LIFT.md)
to work in this range. To obtain radial confinement on the
smaller disk, compare the two half moduli by integrating their
radial derivative from the unit circle.

## Centered transport inside the smaller disk

Write the symmetric Gaussian reference as
\(I_j(z)+z^jI_j(1/z)\) with positive half coefficients.
The saddle estimate on \(|z|\le21/20\) gives
\[
 H_j=\log I_j=jF+C+u_j,\quad
 \|F\|\le3/20,\quad \|DF\|\le1/10,\quad
 \|C\|\le1,\quad \|u_j\|\le1/4,\quad \|u_N-u_n\|\le1/800,
                                                               \tag{2}
\]
where \(D=z\,d/dz\). Here \(\|F\|\le3/20\) follows from
\(F=-a(q/(1+q))^2+\log(1+q)\), \(a\le4\),
\(|q|\le9/100\), and \(|q/(1+q)|\le1/10\).
The leading functions \(F,C\) are common to the orders.

Suppose the coefficient at index \(k\) is multiplied by
\(e^{E_j(k/j)}\), where \(E_j\) is holomorphic on \(|v|\le2048\),
has conjugation symmetry and vanishes at \(v=0,1\). We shall use
\[
        \|E_j''\|\le j\omega,\quad
        \|E_N''-E_n''\|\le\omega,\qquad \omega\le10^{-11}.       \tag{3}
\]
Write \(\Phi_j=E_j(\cdot)\) or \(E_j(1-\cdot)\) for the two
halves. On \(|z|\le26/25\) define \(w_j,h_j,Z_j\) as in the
[curvature-eight lift](MULTIVARIATE_DIFFUSE_S8_CENTERED_ANGLE_LIFT.md#centered-multiplier-on-the-smaller-saddle-disk).
Here the fixed-point disk is \(|w|\le1/8\). The endpoint zeros
and second-derivative bound give \(|h_j|\le\omega\) and
\(|Z_j|\le1041/1000\). Cauchy's bound \(|D^2F|\le22\) on
\(|z|\le209/200\) gives contraction constant at most \(22\omega\).
Use the exact centered multiplier identity (4) of the
[curvature-eight lift](MULTIVARIATE_DIFFUSE_S8_CENTERED_ANGLE_LIFT.md).
Its centered exponent \({\cal R}_j\) has bound \(B|v-w_j|^2\)
on \(|v-w_j|\le1024\), where \(B=j\omega/2\).

For the normalized centered derivatives,
the Cauchy radius in the log dilation is \(\sqrt{m/j}/512\)
when \(m\le j\). The evaluation points remain within \(209/200\).
The first-order saddle term cancels and the quadratic bound
is \(11|h|^2\). The quotient of amplitude terms costs at most
\(e^{5/2}<16\), so
\(|D_{j,m}|\le16(513\sqrt{m/j})^m\).
For \(m>j/4\) use the positive coefficients instead.
The half quotient is at most \(16e^{3j/10}\), and
\((11/10)e^{6/5}<4\) gives \(|D_{j,m}|\le16(4)^m\).

Let \(c_m\) be the coefficients of \(e^{{\cal R}_j}-1\)
in powers of \(v-w_j\). As in the curvature-eight proof,
sum \(|c_m D_{j,m}|\) for \(m\ge2\) using the following bounds.

- If \(m\le j/4\) and \(m\le1024^2j\omega\), Cauchy's circle
  at radius \(\sqrt{m/(j\omega)}\) gives
  \(16(855\sqrt{\omega})^m\), with sum less than
  \(12000000\omega\).
- If \(m\le j/4\) and \(m>1024^2j\omega\), use radius 1024
  with \(|e^u-1|\le |u|e^{|u|}\). The first two powers
  of \(513\sqrt{m/j}/1024\) yield a factor \(m/j\);
  the remaining powers, combined with \(e^{m/2}\), have
  ratio at most \(21/50\). The sum is at most
  \[
       24\cdot513^2\omega
             \sum_{m\ge2}m(21/50)^{m-2}<31000000\omega .
  \]
- If \(m>j/4\), the same radius gives
  \(8\cdot1024^2j\omega(1/200)^m\), since
  \(e^{(1024^2/2)j\omega}\le (5/4)^m\).
  Using \(j<4m\), the sum is less than \(2000\omega\).

The case \(\omega=0\) has no correction. Otherwise the
relative error is less than \(44000000\omega\).
Taking its logarithm costs less than \(45000000\omega\),
and moving \(C+u_j\) back from \(Z_j\) to \(z\) costs less
than \(145\omega\) by Cauchy's first derivative bound.
Hence on \(|z|\le26/25\) each actual half has an analytic log
\[
 \log\widetilde I_j=H_j+W_j+e_j,\qquad
 \|W_j\|\le j\omega,\quad\|W_N-W_n\|\le2\omega,\quad
 \|e_j\|\le\varepsilon:=50000000\omega .                       \tag{4}
\]
The leading correction is
\(W_j=j(F(Z_j)-F(z))+\Phi_j(w_j)-w_j\Phi_j'(w_j)\).
At fixed \(w\) the log shifts differ by at most \(2\omega/N\).
The contraction gives \(|w_N-w_n|\le45\omega/N\) and then
\(|h_N-h_n|\le3\omega/N\). Together with \(\|DF\|\le1/10\)
and \((\Phi_j(w)-w\Phi_j'(w))'=-w\Phi_j''(w)\) these prove
the bounds on \(W_j\) in (4).

## Radial confinement and paired phases

Use the half-log difference \(J_j(q,\theta)\) and its Gaussian
version \(J_j^0\) as in the curvature-eight proof.
For \(|q|\le1/100\) the evaluation points have modulus
below \(51/50\). On this disk, (2) and Cauchy give
\[
 |D H_j|\le11j/100,\qquad |D^2 H_j|\le4j.
\]
Indeed \(|D^2F|\le34/10\) and \(\|C+u_j\|\le5/4\).
The bounds on the corrections in (4) cost \(51M\) and \(5253M\)
for these derivatives, respectively, where \(M=j\omega+\varepsilon\).
Thus
\[
          |\partial_q J_j|\le j/4,\qquad
          |\partial_q\partial_\theta J_j|\le9j .
\]
On \(q=1/100\), conjugation at \(q=0\) and the bound on
\(D H_j\) give \(\Re J_j^0(q,\theta)\le(22/100)jq\).
Adding \(2(j\omega+\varepsilon)\) still gives strictly less
than \(jq\). Rouché for the two halves on this circle, and
for the reciprocal polynomial, confines every actual zero to
\(|q|<1/100\).

At any actual root, the real part of the exact phase equation
gives \(|q_j|\le3(\omega+\varepsilon/j)\).
At a common \(q\), the baseline order difference is formed
from \(H_N-H_n=F+u_N-u_n\). Its \(q\)-derivative has modulus
less than one and its mixed \(q,\theta\) derivative has modulus
less than fourteen. The former bound follows from
\(2(1/10+34/800)<1\), and the latter from
\(2(34/10+(34+2\cdot34^2)/800)<14\).
At \(q=0\) its real part vanishes. Consequently, at roots
with a common argument \(\theta\), the real difference
of their two equations gives
\[
 N|q_N-q_n|
 \le \tfrac43(2|q_n|+4\omega+4\varepsilon)
 \le14\omega+6\varepsilon .
\]
The imaginary difference from the baseline at \(q=0\) is
bounded by
\[
 \bigl(9(14\omega+6\varepsilon)+42(\omega+\varepsilon/n)
                    +204\omega+204\varepsilon\bigr)\theta .
\]
These terms respectively change the upper radius, reset the
common baseline radius, compare the leading corrections, and
compare the individual errors. Conjugation gives zero imaginary
parts at \(\theta=0\), so each derivative bound integrates with
a factor \(\theta\). The total is less than \(7\theta/50\).
On \(|z|=1\), \(|D(F+u_N-u_n)|\le1/10+20/800=1/8\).
The baseline imaginary difference is therefore bounded by
\(\theta/4\). With \(0<\theta\le\pi\), subtraction of the
two exact root equations would put a multiple of \(2\pi\)
strictly between zero and \(2\pi\). This proves separation
under (3).

## Ordinary and signed applications

Use the exact Gamma interpolation in the curvature-eight
proof. On \(|v|\le2048\), \(|x_i-L_i|\le2049|d_i|\),
\(\Re x_i\ge(99/100)L_i\) under (1), and \(1/L_i\le\eta\).
The same trigamma integral gives
\[
 \|A_j''-S\|
  \le\left(\frac{2049}{99/100}+\frac{1}{2j(99/100)^2}\right)S\eta,
 \qquad
 N\|A_N''-A_n''\|\le \frac{S\eta}{2n(99/100)^2}.
\]
These prove (3) with \(\omega=2072 S\eta\).

For the signed criterion, let \(r_i=nT_i/(n-1)\), \(A=\sum_i r_i\),
\(P=\sum_i\epsilon_i r_i^2\), and \(V_p=\sum_i r_i^p\).
A sufficient condition at \(n\ge10000\) is
\[
      6\le2A\le8,\qquad
      8194|P|+45000000 V_3+4V_2/n\le10^{-11}.                 \tag{5}
\]
The proof of the signed curvature-eight criterion applies
with \(|g|=|2v-1|\le4097\). The cubic moment bound gives
every lower midpoint parameter below \(1/(8\cdot4097)\).
The upper matched parameters are less than twice as large,
so \(\Re(1+\epsilon_i r g)\ge3/4\) along the parameter paths.
The signed curvature bound uses \(2|g|\le8194\) and
\(2\cdot4097^2/(3/4)\) for the reciprocal expansion.
The same-\(r\) grid drift and the adjacent parameter drift
then give (3), since
\(2\cdot4097^2/(3/4)+2(1/4+5/12)/(3(3/4)^3)<45000000\).
