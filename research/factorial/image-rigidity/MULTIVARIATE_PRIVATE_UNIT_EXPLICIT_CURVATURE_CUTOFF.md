# An explicit positive curvature angular cutoff for private unit pairs

Consider a two-monomial ordinary factorial profile with exactly
one private unit row in each orientation and arbitrary shared rows.
Discard equal rows. For each remaining shared row let \(\ell_i\)
be its smaller exponent and \(d_i>0\) its absolute gap. Put
\[
 V=\sum_i d_i^2/\ell_i,\qquad
 \kappa=\sum_i d_i^2/(\ell_i+d_i+1/2).
\]
The case \(\kappa=0\) gives consecutive geometric sums, whose
nontrivial roots have arguments with coprime successive
denominators and are disjoint. Assume \(\kappa>0\).
Define
\[
 \begin{gathered}
 q=(1+\kappa/8)^{-1},\quad s=1-q,\quad A=16/s^2,\\
 U_0=(2V+4)/s^2,\qquad U_1=(4V+4)/s^3,\qquad
 L_0=4U_0/s,\qquad
 L_1=4U_1/s+16U_0/s^4,\qquad L_2=24/s^6 .
 \end{gathered}
\]

**Theorem.** At consecutive orders \(n,N=n+1\), the root
arguments are disjoint whenever
\[
 n\ge\max\left\{4,\ \frac{8A}{\kappa},\
         \frac{16V}{\kappa s^3},\
         \frac{16}{s^3},\ 8L_1,\
                   \sqrt{16 L_2(A+2L_0)}\right\}.             \tag{1}
\]
The bound allows independent positive rescalings of the two
moment polynomials. It is uniform on every region
\(\kappa\ge\kappa_0>0,\ V\le V_0\), by substituting those bounds
into its right side. It can also replace the absolute cutoff
in the near-mirror and signed diffuse private-unit angular
criteria for a specific profile. The proof uses two elementary
geometric sums and the exact linked coefficient order.
[Near-mirror criterion](MULTIVARIATE_PRIVATE_UNIT_NEAR_MIRROR_QUADRATIC_ANGLES.md)
[Signed diffuse criterion](MULTIVARIATE_PRIVATE_UNIT_DIFFUSE_SIGNED_TANH_ANGLES.md)

Both \(\kappa\) and \(V\) are exact functions of the
reconstructed ordinary rows from any one joint matched tanh pair.
[Matched Gauss grid](MULTIVARIATE_MATCHED_TANH_INTEGER_MULTIPLICITY.md#recovering-the-gauss-grid)

## Nonvanishing endpoint halves and root radii

Cancel the private unit pair against the binomial coefficient
and normalize the constant and leading coefficients to one
by positive scalings. Write \(c_{m,k}\) for the resulting
ordinary moment coefficient. For integer \(m\ge2\), continuous
interpolation in \(k\) by the log-Gamma factors gives
\[
       \kappa/m\le \partial_k^2\log c_{m,k}\le V/m.
\]
Indeed \(\psi'(x+1)=\sum_{j=1}^\infty(x+j)^{-2}\) lies between
\(1/(x+1)\) and \(1/x\) by integral comparison, and every
row argument \(x\) lies in \([m\ell_i,m(\ell_i+d_i)]\).
The endpoints of \(\log c_{m,k}\) are zero.
The Green function for the chord gives
\[
 -V k(m-k)/(2m)\le\log c_{m,k}\le-\kappa k(m-k)/(2m).
                                                               \tag{2}
\]
Also \(\partial_k\log c_{m,k}|_{k=0}\le-\kappa/2\), by
integrating the lower second derivative against \(1-k/m\).
These same facts hold in the reversed orientation.

Split the polynomial exactly as \(F_m(z)=E_m(z)+z^m D_m(1/z)\).
The first half uses the coefficients at \(0\le k\le\lfloor m/2\rfloor\).
The other uses the reversed coefficients at
\(0\le k\le\lfloor(m-1)/2\rfloor\).
Take \(m\) to be \(n\) or \(N\). By (2), each half has terms whose
absolute values on \(|z|\le e^{A/n}\) are at most \(q^k\).
Indeed the log coefficient is at most \(-\kappa k/4\) on
either half, and \(A/n\le\kappa/8\).

To bound each half below, first truncate it at
\(K_m=\lfloor\kappa m/(4V)\rfloor\), which lies in both halves
for \(m\ge4\).
For \(k\le K_m\), the log coefficient slope is at most
\(-\kappa/2+Vk/m\le-\kappa/4\).
On any radius \(r\le e^{A/n}\), the consecutive terms
\(a_k=c_{m,k}r^k\) of the truncated half have ratios at most \(q\).
For the truncated polynomial in \(w=z/r\), \(|w|=1\),
multiplication by \(1-qw\) gives constant term one and
nonpositive remaining coefficients. Its value at \(w=1\)
is \(s\sum_{k=0}^{K_m}a_k\ge s\). The triangle inequality
therefore bounds the absolute value of the truncated half
below by \(s/2\) on that radius. At \(r=0\) the bound follows
from the constant coefficient. The truncated remainder has
absolute value at most \(q^{K_m+1}/s\le s/4\), since
\[
 K_m+1\ge4/s^3,\qquad q^h=(1-s)^h\le(1+hs)^{-1}
\]
for positive integer \(h\).
This proves \(|E_m|,|D_m|\ge s/4\) throughout the common
disk. The upper bound for both is \(1/s\).

On \(|z|=e^{-A/m}\), the reflected term \(z^mD_m(1/z)\)
has modulus at most \(e^{-A}/s<s/4\).
Rouché excludes roots inside that circle, since \(E_m\)
has no zeros on the common disk. Apply the same argument
to the reversed polynomial with \(D_m\) as the constant
half to exclude roots outside \(e^{A/m}\).
Thus every normalized root has \(|\log|z||\le A/m\).

## Linked variation of the half phases

Let \(H_m\) denote either half in its nonnegative-power
orientation, and \(D_z=z\,d/dz\).
On the common disk, the geometric sums give
\[
 |H_m|\le 1/s,\qquad |D_z H_m|\le1/s^2,\qquad
                                      |D_z^2 H_m|\le2/s^3 .
\]
At an index \(k\le n/2\) in both half arrays, the
[linked coefficient order](MULTIVARIATE_BOTH_PRIVATE_LINKED_COEFFICIENT_ORDER.md)
and (2) give
\[
   0\le c_{n,k}-c_{N,k}\le c_{n,k}\,2Vk/n
\]
in either orientation. Indeed the exponent in the lower
bound of that comparison exceeds one by at most \(4/n\),
and \(-\log c_{n,k}\le Vk/2\).
Each half can gain at most one coefficient with \(k\ge n/2\).
For that index \(kq^k\le4/(n s^2)\), since
\(\sum_{j\ge0}(j+1)q^j=1/s^2\) and its first \(k+1\)
terms sum to at least \((k+1)(k+2)q^k/2\).
Consequently
\[
  |H_N-H_n|\le U_0/n,\qquad |D_z(H_N-H_n)|\le U_1/n .
\]
Both half polynomials admit logarithms on the disk, real
on its positive radius. The lower bound \(s/4\) and the
preceding estimates give
\[
 \begin{aligned}
 |\log|H_N|-\log|H_n||&\le L_0/n,\\
 |D_z\log H_N-D_z\log H_n|&\le L_1/n,\\
 |D_z\log H_m|&\le4/s^3,\qquad |D_z^2\log H_m|\le L_2 .
 \end{aligned}
\]
The first line uses the Lipschitz bound of the real
logarithm above \(s/4\). The second follows by subtracting
the two ratios \((D_z H_m)/H_m\). In the last line,
the second derivative costs at most \(8/s^4+16/s^6\).

Suppose there are roots with a common upper-half argument
\(\theta\in(0,\pi]\), and set \(x_m=\log|z_m|\) at those roots.
Write
\[
 K_m(x,\theta)+i\Phi_m(x,\theta)
        =\log E_m(e^{x+i\theta})-\log D_m(e^{-x-i\theta}).
\]
The modulus equation is \(m x_m=K_m(x_m,\theta)\).
Subtract the two equations at the same \(\theta\),
using \(|\partial_x K_m|\le8/s^3\),
\(|K_N-K_n|\le2L_0/n\) at fixed \(x\), and \(|x_n|\le A/n\).
Because \(N-8/s^3\ge N/2\), the radius difference satisfies
\[
                      |x_N-x_n|\le 2(A+2L_0)/(nN).
\]
At fixed \(x_n,x_N\), both half log differences are real
at \(\theta=0\). The derivative in \(\theta\) of their
imaginary difference is bounded by
\(2L_1/n+2L_2|x_N-x_n|\le1/2\), by (1).
Integrating from zero gives
\[
                   |\Phi_N(x_N,\theta)-\Phi_n(x_n,\theta)|
                               \le\theta/2.
\]
The root phase equations would therefore put
\(\theta-[\Phi_N(x_N,\theta)-\Phi_n(x_n,\theta)]\)
in \(2\pi\mathbb Z\). This number lies in
\([\theta/2,3\theta/2]\subset(0,2\pi)\), a contradiction.
The positive real ray is root-free by positivity and the
lower half-plane follows by conjugation.
