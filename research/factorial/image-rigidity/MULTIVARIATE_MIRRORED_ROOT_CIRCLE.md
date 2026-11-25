# Simple unit-circle roots for mirrored binomial moments

For mirrored exponent profiles with equal Gamma shapes within each
paired row, every fixed-order moment polynomial has simple unit-circle
roots. Private exponents are unrestricted. A pair of such polynomials
can share a root only on the circle, and the successive orders remain a
separate coprimality question. In particular, applying Schur–Szegő
composition separately at each degree does not in general preserve
coprimality of two successive polynomials, even when the input pairs
strictly interlace and have positive log-convex binomial-normalized
coefficients.

## Mirrored profiles at every fixed order

**Theorem.** Fix integer profiles \(b_i\ge0\) and \(d_i\ge1\),
\(1\le i\le r\), \(r\ge1\), and positive real Gamma shapes
\(\gamma_i\). Give the two independent variables \(X_i,Y_i\) common
shape \(\gamma_i\) and rate one. For
\[
 U=\prod_{i=1}^r X_i^{b_i+d_i}Y_i^{b_i},\qquad
 V=\prod_{i=1}^r X_i^{b_i}Y_i^{b_i+d_i},
\]
the polynomial \(Q_m(t)=\mathbb E[(U+tV)^m]\) has \(m\) distinct
zeros on \(|t|=1\) for every \(m\ge1\). Any additional independent
Gamma coordinates with equal exponents in both monomials contribute
only a nonzero scalar to \(Q_m\).

The coefficient of \(t^k\), divided by \(\binom mk\), is
\[
 \prod_i (\gamma_i)_{m b_i+d_i k}
          (\gamma_i)_{m b_i+d_i(m-k)}.
\]
Grouping the factors of a rising factorial by their residues modulo
\(d_i\) gives
\[
 (A)_{d_i k}
   =d_i^{d_i k}\prod_{j=0}^{d_i-1}
              \bigl((A+j)/d_i\bigr)_k .
\]
Here \(A=\gamma_i+m b_i>0\). After removing factors independent of
\(k\), the normalized coefficient array of \(Q_m\) is therefore the
entrywise product of the arrays
\((\lambda)_k(\lambda)_{m-k}\), one for each
\(\lambda=(\gamma_i+m b_i+j)/d_i>0\).

For any \(\lambda>0\) use the generating series
\[
 (1-2xv+v^2)^{-\lambda}=\sum_{m\ge0} C_m^\lambda(x)v^m .
\]
The corresponding binomial polynomial
\[
 P_{m,\lambda}(t)
   =\sum_{k=0}^m\binom mk(\lambda)_k(\lambda)_{m-k}t^k
\]
satisfies
\(P_{m,\lambda}(e^{2i\theta})
   =m! e^{im\theta} C_m^\lambda(\cos\theta)\).
This follows by factoring the generating function into the two
factors with arguments \(e^{i\theta}v\) and \(e^{-i\theta}v\).
For completeness, the standard orthogonality argument here works
for all positive \(\lambda\), including \(\lambda<1\). Apply
\((1-x^2)\partial_x^2-(2\lambda+1)x\partial_x+
(v\partial_v)^2+2\lambda v\partial_v\) to the generating function
to get zero. The coefficient \(y=C_m^\lambda\) then satisfies
\[
 \bigl((1-x^2)^{\lambda+1/2} y'\bigr)'
 +m(m+2\lambda)(1-x^2)^{\lambda-1/2}y=0.
\]
The eigenvalues \(m(m+2\lambda)\) are distinct and the boundary
terms vanish at \(\pm1\). Integration by parts gives orthogonality
for the positive integrable weight \((1-x^2)^{\lambda-1/2}\).
Since the polynomials have exact degrees \(m\), orthogonality to
every lower degree forces \(m\) sign changes in \((-1,1)\). Indeed,
multiplying by the linear factors at all sign-changing zeros would
otherwise give a nonzero integrand of one sign with a multiplier
of degree less than \(m\). Hence every \(P_{m,\lambda}\) has simple
unit-circle roots.

Use the degree-\(m\) Schur–Szegő composition
\[
 \Big(\sum_k\binom mk a_k t^k\Big)\star_m
 \Big(\sum_k\binom mk b_k t^k\Big)
 =\sum_k\binom mk a_k b_k t^k .
\]
Here is its disk argument with the hypotheses made explicit.
Suppose the first factor \(F\) has no zero in \(|t|<1\) and the
second is monic with unit-circle zeros \(\zeta_1,\ldots,\zeta_m\).
Then their composition is
\[
 (-1)^m\Big(\prod_{i=1}^m\zeta_i\Big)
 \widetilde F(-t/\zeta_1,\ldots,-t/\zeta_m),\qquad
 \widetilde F(x_1,\ldots,x_m)=\sum_{k=0}^m a_k e_k(x_1,\ldots,x_m),
\]
where \(e_k\) denotes the elementary symmetric polynomial.
This is the symmetric multi-affine polarization of \(F\). For
\(|t|<1\) its displayed arguments all belong to the open unit disk.
The Grace–Walsh–Szegő coincidence theorem for this convex circular
domain makes its value equal to \(F(\eta)\) at some \(|\eta|<1\).
Thus the composition has no interior root. For real palindromic
inputs of exact degree \(m\) the composition remains palindromic
of degree \(m\), so reciprocity excludes exterior roots as well.
The version of the coincidence theorem used here is Theorem 11,
pp. 478–479 in [Borcea and Brändén, Pólya–Schur master theorems for circular domains and their boundaries](https://annals.math.princeton.edu/wp-content/uploads/annals-v170-n1-p14-p.pdf).

Simplicity follows by openness. In the real palindromic space of
exact degree \(m=2l\), write \(F(t)=t^l H(t+t^{-1})\), where \(H\)
is real of degree \(l\). At \(m=2l+1\), first remove \(t+1\) and
use the same expression for its palindromic quotient. The
circle-rooted set corresponds exactly to \(H\) having all roots
in \([-2,2]\), and its relative interior corresponds exactly to
distinct roots in \((-2,2)\). At an endpoint one can move a root
outside the interval, and at a repeated root one can perturb two
equal real factors to a nonreal conjugate pair. In the other
direction simple interior roots persist as real interior roots
under small real coefficient changes. Distinct interior roots of
\(H\) are equivalent to simple unit-circle roots of \(F\), including
the simple structural factor \(t+1\) at odd degree. For any positive
palindromic circle-rooted \(P_{m,\lambda}\), composition by it is
an invertible linear map on the real palindromic space. The disk
argument shows it maps the circle-rooted set into itself. An
invertible linear map is open, so it also maps the interior into
itself. Composing all factors \(P_{m,\lambda}\) now proves the
location and simplicity claimed for \(Q_m\).

For example \(U=X^3Y,\ V=XY^3\) has no private exponents. The
theorem still gives simple circular roots at all orders. The
[private-exponent theorem](SFC_OVERLAPPING_SWAP_PAIRS.md)
proves the stronger statement of strict interlacing for its
specified parameter range.

## Obstruction to varying the composition degree

For \(j=1,2\), form monic pairs
\[
 f_j(t)=(1+t)(t^2+b_jt+1),\qquad
 g_j(t)=(t^2+a_j t+1)(t^2+c_j t+1)
\]
with
\[
 (a_1,b_1,c_1)=(-1/2,11/6,15/8),\qquad
 (a_2,b_2,c_2)=(-1/2,7/4,433/234).
\]
For each \(j\), we have \(-2<a_j<b_j<c_j<2\). The identity
\(t^2+a_jt+1=0\) on \(t=e^{i\theta}\) gives
\(a_j=-2\cos\theta\), which increases strictly on \(0<\theta<\pi\).
Together with the cubic root \(-1\) and conjugate symmetry this
proves strict interlacing in each input pair. All coefficients
are positive. In binomial normalization the quartics have
endpoint and interior arrays
\[
 (1,11/32,17/96,11/32,1),\qquad
 (1,79/234,503/2808,79/234,1).
\]
Each array is log-convex because its middle coefficient \(v\)
satisfies \(u^2\le v\le u\) for the adjacent coefficient \(u\).
The cubics likewise have normalized arrays \((1,u,u,1)\)
with \(0<u<1\).

Exact coefficient multiplication gives
\[
\begin{aligned}
 f_1\star_3 f_2 &= (1+t)(t^2+\tfrac{115}{72}t+1),\\
 g_1\star_4 g_2 &=
     (t^2+\tfrac{115}{72}t+1)(t^2-\tfrac{707}{624}t+1).
\end{aligned}
\]
For the cubic the interior coefficient is
\((17/6)(11/4)/3=187/72\). The quartic's coefficients at
\(t\) and \(t^2\) are \((11/8)(158/117)/4=869/1872\) and
\((17/16)(503/468)/6=8551/44928\), which agree with the
displayed factors. Both outputs remain circle-rooted, but they
share the nonreal pair from the common quadratic. This example
concerns arbitrary input polynomial pairs. The Gamma coefficient
relations across orders would require their own argument.
