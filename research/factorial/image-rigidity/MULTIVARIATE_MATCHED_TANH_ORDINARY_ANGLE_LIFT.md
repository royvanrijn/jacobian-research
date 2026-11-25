# A matched tanh condition for ordinary angular separation

Use the [matched tanh factors](MULTIVARIATE_MATCHED_TANH_MOMENTS.md)
of an ordinary two-monomial profile at \(n\ge2,\ N=n+1\).
The index ranges over every Gauss factor, and \(\epsilon_i=1\)
or \(-1\) according to its row orientation. Put
\[
 \lambda=\frac{n}{n-1},\qquad
 r_i=\lambda T_i,\qquad r_i'=\lambda T_i'/R,\qquad
 A=\sum_i r_i,\qquad P=\sum_i\epsilon_i r_i^2,\qquad
 V_p=\sum_i r_i^p,\quad R=\frac{n^2}{n^2-1}.
\]

**Theorem.** The two exact ordinary moment polynomials have disjoint
root arguments under independent positive rescalings if
\[
       2A\ge13,\qquad
       22|P|+324V_3+\frac{4V_2}{n}\le10^{-8}.               \tag{1}
\]
They are therefore coprime over every characteristic-zero field.
This is a sufficient condition in matched moments. In particular,
the signed second moment can cancel between the two orientations.

The proof supplies the complex analytic log-coefficient bounds
required by the [degree-independent phase lift](MULTIVARIATE_DIFFUSE_ABSOLUTE_ANGLE_LIFT.md).
The common Gaussian reference parameter is \(2A\).

## Exact curvature of the assembled symbol

At degree \(j\), a single increasing row factor has midpoint
parameter \(r=jT/(j-1)\). Put \(c=(1-r)/(2r)\).
The contribution to the normalized log coefficient divided by \(j\)
is the centered chord of
\(j^{-1}\log\Gamma(j(c+s)+1/2)\) at \(s=k/j\).
Indeed the factor at midpoint \((l+1/2)/j\) in the Gauss product is
\[
 1-r+2r(l+1/2)/j=(2r/j)\big(jc+l+1/2\big).
\]
The \(2r/j\) terms cancel in the chord. A decreasing factor
uses \(s\mapsto1-s\).
The midpoint parameters at the two orders are exactly \(r_i,r_i'\)
because \(\lambda/R=N/n\).

For either orientation write \(g=2s-1\),
\(x=(1+\epsilon r g)/(2r)\), and
\[
  k_j(s;r,\epsilon)=j\psi'(jx+1/2)
      =\int_0^\infty e^{-xu}\,D(u/j)\,du,\qquad
  D(v)=\frac{v}{2\sinh(v/2)} .
\]
This is the exact second \(s\)-derivative of the normalized chord.
The integral follows by expanding the geometric series for
\(1/(1-e^{-u/j})\) in the trigamma integral.
Let \(L_j(s)\) be the sum of the analytic normalized chords.

For real \(v\ge0\), \(0<D(v)\le1\) and
\(0\le-dD/d(v^2)\le1/24\).
Indeed put \(y=v/2\). The upper inequality is equivalent to
\(3(y\cosh y-\sinh y)\le y\sinh^2 y\).
The coefficient of \(y^{2m+1}\) on the left, multiplied
by \((2m)!\), is \(6m/(2m+1)\), and on the right it is
\(2^{2m-1}\). They agree at \(m=1\), and the latter is larger
for \(m\ge2\). The lower bound follows from
\(y\cosh y-\sinh y\ge0\).
Consequently throughout \(\Re x>0\),
\[
 \left|k_j-\frac1x\right|
      \le \frac1{12j^2(\Re x)^3},\qquad
 |k_N-k_n|\le
       \frac{1/n^2-1/N^2}{12(\Re x)^3}                     \tag{2}
\]
for a fixed parameter \(r\). Apply the slope inequality to the
integrands and use the integral of \(u^2e^{-(\Re x)u}\).


## Size and adjacent change of the log-coefficient correction

The pointwise matched identity gives
\[
 r_i'=\frac{r_i}{1-\xi_i},\qquad |\xi_i|\le\frac{r_i}{nN}.
\]
In (1), \(324V_3\le10^{-8}\) forces every \(r_i<1/1000\).
It follows that \(r_i'<2r_i\), so both parameter lists are
bounded by \(1/44\).
On the disk \(|s|\le5\), \(|g|\le11\), and hence
\(\Re x\ge d/(2r)\) with \(d=3/4\).
Use the expansion
\[
 \frac1x=2r-2\epsilon r^2 g+
                 \frac{2r^3g^2}{1+\epsilon r g}.
\]
Summing and using (2) gives
\[
       \|L_n''-2A\|_{|s|\le5}
              \le22|P|+\left(\frac{242}{d}
                          +\frac{2}{3n^2d^3}\right)V_3.     \tag{3}
\]
The signed \(r^2\) term is summed before its absolute value is taken.

For the changing parameters, differentiate the integral for \(k_N\)
in real \(r>0\), obtaining
\[
       |\partial_r k_N|
          \le\frac1{2r^2}\int_0^\infty u e^{-(\Re x)u}\,du
          \le\frac2{d^2}.
\]
The matched identity gives
\(|r_i'-r_i|\le r_i^2/(nN-1/44)\).
Together with the fixed-parameter difference in (2), this gives
\[
 \|L_N''-L_n''\|_{|s|\le5}
 \le \frac{2V_2}{d^2(nN-1/44)}
       +\frac{2 V_3}{3d^3}(1/n^2-1/N^2).                  \tag{4}
\]
A profile-specific alternative for the first term in (4) follows
from the [full-kernel transport theorem](MULTIVARIATE_MATCHED_TANH_KERNEL_TRANSPORT.md).
Apply it to \(t\mapsto k_N(s;\lambda t,\epsilon)\) row by row.
The second derivative combination cancels the first derivative
term of the kernel exactly:
\[
 (2\partial_r+r\partial_r^2)k_N
       =\frac1{4r^3}\int_0^\infty u^2e^{-xu}D(u/N)\,du .
\]
Thus the first term in (4) can also be bounded by
\((W-RU)\max\{2C_+\lambda/d^2,4C_-\lambda/d^3\}\), with
\(C_+,C_-\) as in that theorem. Only the compact positive
range of the paths is needed for this application.

Set \(E_j(s)=j[L_j(s)+A s(1-s)]\).
It has zero endpoint values, conjugation symmetry, and gives the
exact logarithmic multiplier relative to \(B_j(z;2A)\)
on the fractional grid \(s=k/j\).
For \(n\ge2\), direct comparison gives
\[
 \frac{2N}{d^2(nN-1/44)}\le\frac4n,\qquad
 N(1/n^2-1/N^2)\le5/12,\qquad
 \frac{242}{d}+\frac{2}{3n^2d^3}+\frac{2}{3d^3}\frac5{12}
 \le\frac{968}{3}+\frac{32}{81}+\frac{160}{243}<324 .
\]
From (3)–(4) we obtain
\[
 \max\left\{\|E_n''\|/n,\ \|E_N''\|/N,\
                 \|E_N''-E_n''\|\right\}
           \le22|P|+324V_3+4V_2/n .
\]
All norms are on \(|s|\le5\).
Condition (1) gives the two analytic norm bounds of the
degree-independent lift with tolerance \(10^{-9}\)
and reference parameter \(2A\ge13\).
That lift compares the exact exponentiated coefficients and
the phases after allowing independent positive rescalings.

For example, take the two reflected rows
\((\ell,\ell+d_0),(\ell+d_0,\ell)\) with
\(\ell=600000000000\) and \(d_0=2000000\).
Then \(P=0\) and each \(r_i=d_0/(2\ell+d_0+(2h-d_0)/n)\)
lies between \(d_0/(2\ell+2d_0)\) and \(d_0/(2\ell)\).
Consequently
\[
 A\ge d_0^2/(\ell+d_0)>13/2,\qquad
 V_2\le1/90000,\qquad V_3\le1/54000000000.
\]
Condition (1) holds for every \(n\ge15000\).
Here the unsigned relative gap condition \(S\eta\le10^{-9}\)
of the original diffuse lift fails.

