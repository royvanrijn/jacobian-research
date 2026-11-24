# The convergent core of arbitrary angular Gaussian monomial fibres

Let \(Z,W\) be a normalized circular Wick pair, so
\(\mathbb E[Z^uW^v]=\delta_{u,v}u!\). Let \(Y_1,\ldots,Y_s\) be independent
standard ordinary Gaussians, independent of the Wick pair. Assume \(s\ge1\)
and \(k_i\ge1\). Set
\[
 B(Y)=\prod_{i=1}^sY_i^{k_i},\qquad P_0=Wc(Z)B(Y),\qquad c\in\mathbb C[Z].
\]
For a phase \(P\), write \(\Phi_P[A]=\mathbb E[Ae^{tP}]\) for the actual
formal Wick period, including for polynomial-time amplitudes. Put
\[
 M_P=\mathbb C[t]\Phi_P[\mathbb C[Z,W,Y]],\qquad
 V_P=\mathbb C(t)M_P,\qquad
 {\cal A}=\mathbb C\{t\}[t^{-1}].
\]
Thus \({\cal A}\) is the field of convergent meromorphic germs at zero.
The generic actual space is finite-dimensional and derivative stable by
[Gaussian holonomicity](WHOLE_GAUSSIAN_SOURCE_HOM_AND_EFFECTIVE_CUTOFFS.md).

For every polynomial \(c\), the exact convergent spaces are
\[
 \begin{aligned}
 M_{P_0}\cap\mathbb C\{t\}&=\mathbb C[t],&
 V_{P_0}\cap{\cal A}&=\mathbb C(t),\\
 M_{P_0+a_0}\cap\mathbb C\{t\}&=e^{a_0t}\mathbb C[t],&
 V_{P_0+a_0}\cap{\cal A}&=e^{a_0t}\mathbb C(t).
 \end{aligned}                                                    \tag{1}
\]
The largest globally regular differential submodule is the constant line
for \(a_0=0\) and zero for \(a_0\ne0\).

This gives an obstruction in
\(\texttt{GAP-GAUSSIAN-REGULAR-TRANSCENDENTAL}\) for all angular polynomials
in the monomial fibre, including non-binomial \(c\), any number of ordinary
factors, and arbitrary positive exponents. The
[explicit binomial-fibre result](GAUSSIAN_NONLINEAR_MONOMIAL_FIBRE_FULL_IMAGE.md)
additionally computes the ranks and full lattices in its special family.
Here the proof uses paritywise moment normalization and algebraicity of
the normalized angular coefficients.

## 1. The normalized angular germ

For a monomial amplitude \(A=Z^aW^b\prod_iY_i^{d_i}\), the original
Wick evaluation gives
\[
 [t^n]\Phi_{P_0}[A]
 =\frac{(n+b)!}{n!}\,[z^{\,n+b-a}]c(z)^n\,
                       \prod_i\mathbb E[Y_i^{k_i n+d_i}].          \tag{2}
\]
Coefficients of negative powers of a polynomial are zero. For
\(n=2q+e\) with \(e\in\{0,1\}\), use the common positive weights
\[
 \epsilon_{i,e}\equiv k_i e\pmod 2,\quad \epsilon_{i,e}\in\{0,1\},
 \qquad
 w_{q,e}=\prod_i(2k_i q+k_i e+\epsilon_{i,e}-1)!!.                 \tag{3}
\]
We use \((-1)!!=1\). An inadmissible parity gives zero in (2). On an
admissible parity, \(d_i-\epsilon_{i,e}\) is an even nonnegative
integer. Dividing the moments in (2) by (3) gives a polynomial in \(q\).
The Wick factor \((2q+e+b)!/(2q+e)!\) is also polynomial in \(q\).

For any integer \(h\), define
\[
 A_h(X)=\sum_{n\ge0}[z^{\,n-h}]c(z)^n X^n.
\]
This series is algebraic over \(\mathbb C(X)\). If \(c(0)\ne0\), on a
fixed small circle \(|z|=r\) and for sufficiently small \(X\), termwise
integration of the geometric series gives
\[
 A_h(X)=\frac1{2\pi i}\int_{|z|=r}\frac{z^h}{z-Xc(z)}\,dz .        \tag{4}
\]
For \(X\ne0\) near zero the denominator has exactly one zero inside the
circle by Rouché's theorem, the simple algebraic germ \(z_*(X)\) with
\(z_*=Xc(z_*)\) and \(z_*(0)=0\). Its residue is
\(z_*^h/(1-Xc'(z_*))\). If \(h<0\), the pole at zero contributes a
rational function in \(X\), since it uses only finitely many Taylor
coefficients of \((z-Xc(z))^{-1}\) at \(z=0\). Consequently the entire
germ (4), where the individual Laurent terms are added, is algebraic.
The same selected-branch identity is the classical Lagrange inversion
formula with \(\psi(z)=z^h\).
[Gessel, equation (2.2.2) and §4.1](https://people.brandeis.edu/~gessel/homepage/papers/Lagrange.pdf)

If \(c=z u\) with \(u(0)\ne0\), the coefficient in \(A_h\) equals
\([z^{-h}]u^n\). It vanishes for \(h>0\), and for \(h\le0\) a finite
binomial expansion gives \(u(0)^n\) times a polynomial in \(n\).
Thus \(A_h\) is rational. If \(c=0\) or \(c\) is divisible by \(z^2\),
the same coefficient vanishes for all sufficiently large \(n\), so
\(A_h\) is polynomial. If \(c\) is constant, it is also polynomial,
since the index \(n-h\) must be zero.

For each parity,
\[
 \sum_{q\ge0}[z^{\,2q+e-h}]c(z)^{2q+e}x^q
\]
is algebraic by applying the even or odd projector to \(A_h(X)\) and
substituting \(X=\sqrt{x}\), with a factor \(X^{-e}\). Applying a
polynomial in \(x\partial_x\) preserves algebraicity by differentiation
of a separable algebraic relation. In conjunction with (2)–(3), this
proves that the normalized parity germ of every monomial period is
algebraic.

## 2. Euler denominator clearing for the full actual core

Let \(H\in M_{P_0}\) be convergent. It is a finite sum of terms
\(t^j\Phi_{P_0}[A]\) with \(j\ge0\) and monomial \(A\), including scalar
coefficients. Fix an output parity \(e\). In such a term the source
coefficient at output index \(2q+e\) has index \(2(q+v)+e'\) for fixed
\(v\in\mathbb Z\), \(e'\in\{0,1\}\). For all sufficiently large \(q\),
the normalized output coefficient is the coefficient of an algebraic
parity germ shifted by \(v\), multiplied by
\[
 \frac{w_{q+v,e'}}{w_{q,e}}.                                      \tag{5}
\]
For each \(i\), the two double-factorial arguments in (5) differ by a
fixed even integer. Their ratio, for large \(q\), is a finite product
or reciprocal of affine factors in \(q\). Thus (5) is rational in
\(q\). A shift of an algebraic germ changes it by a monomial factor
and a finite Taylor subtraction, so the underlying shifted series
is still algebraic.

Choose a nonzero polynomial \(Q_e(q)\) clearing all the finitely many
rational denominators for output parity \(e\). Define
\[
 U_e(x)=\sum_{q\ge0}\frac{[t^{2q+e}]H}{w_{q,e}}x^q .
\]
Then \(Q_e(x\partial_x)U_e\) is a finite sum of polynomial Euler
operators on algebraic germs up to a finite polynomial correction.
The correction accounts for the indices where a shifted source
coefficient is absent or (5) was used only as a tail identity. Hence
\(Q_e(x\partial_x)U_e\) is algebraic.

Convergence of \(H\) bounds its Taylor coefficients exponentially.
Since \(s\ge1\), \(k_i\ge1\) and \(w_{q,e}^{1/q}\to\infty\), \(U_e\)
and \(Q_e(x\partial_x)U_e\) are entire. An entire algebraic function
is polynomial. Indeed, outside a large disk the leading polynomial
coefficient of an algebraic relation has no zero. The root bound
for its values gives polynomial growth, and Cauchy's estimates
then annihilate all sufficiently high derivatives at zero. Since
\(Q_e(q)\ne0\) for all sufficiently large \(q\), \(U_e\) itself is
polynomial. Both parities now give \(H\in\mathbb C[t]\).

The time-polynomial angular Jacobian marker gives
\[
 \Phi_{P_0}[1-tc'(Z)B(Y)]=1.                                     \tag{6}
\]
At \(n>0\) the cancellation follows from
\([z^n]c^n=[z^{n-1}]c'c^{n-1}\), by differentiation. The ordinary
Gaussian moment accompanying both coefficients is \(\mathbb E[B^n]\).
This proves the reverse inclusion for the polynomial-time core.
Clearing a nonzero polynomial time denominator proves the generic
claim, including meromorphic germs with a finite pole at zero.

A constant phase term \(a_0\) multiplies all periods by \(e^{a_0t}\),
a convergent unit at zero, establishing the two twisted claims in (1).
Any globally regular actual submodule consists of convergent Laurent
germs by the Frobenius theorem, applied to scalar cyclic operators.
Thus it lies in the displayed rank-one line. The untwisted line is
regular globally. When \(a_0\ne0\), the connection one-form on the line
is \(a_0\,dt\), which has a double pole at infinity. A rational basis
change adds a logarithmic derivative, which has at most a simple pole.
The twisted line is irregular, and has no nonzero proper differential
submodule. [Fuchs–Frobenius theory](https://dlmf.nist.gov/2.7.i)

## 3. Exact diagonal cases

Consider
\[
 P=P_0+a_0+\sum_i a_iY_i^2,\quad D_i=1-2a_i t,\quad
 G(t)=e^{a_0t}\prod_i D_i^{-1/2},\quad
 T=t\prod_i D_i^{-k_i/2}.
\]
All powers use their branches equal to \(1\) at zero. For every
monomial \(A=Z^aW^b\prod_iY_i^{d_i}\), formal coefficient expansion gives
\[
 \Phi_P[A](t)=G(t)\left(\prod_iD_i^{-d_i/2}\right)
                       \Phi_{P_0}[A](T).                         \tag{7}
\]
Indeed, for even \(r=2v\) the ordinary Gaussian coefficients satisfy
\[
 \mathbb E[Y^{2v} e^{taY^2}]
   =(2v-1)!!(1-2at)^{-v-1/2},
\]
by the binomial series, and the analogous odd moments are zero.
Insert this formula with \(r=k_i n+d_i\) into (2) to obtain (7).
The time map is an analytic coordinate at zero, since \(T'(0)=1\).

For arbitrary \(c\) the polynomial-time amplitude in (6), now
multiplied by \(\prod_i Y_i^{2r_i}\) with \(r_i\ge0\), gives the exact
period
\[
 G(t)\prod_i (2r_i-1)!!D_i^{-r_i}.                                \tag{8}
\]
The circular cancellation used in (6) holds coefficientwise before the
ordinary Gaussian expectation in this instance. It leaves precisely
the diagonal Gaussian moments. Therefore
\(G(t)\mathbb C[t,D_1^{-1},\ldots,D_s^{-1}]\subset M_P\) for every \(c\).

If \(c\) is constant or divisible by \(Z^2\), all its unperturbed
monomial periods are polynomial by (2). Every surviving power \(T^n\)
in (7) pairs with \(\prod_iD_i^{-d_i/2}\) to give
\(t^n\prod_iD_i^{-(k_i n+d_i)/2}\), with integer nonnegative
exponents in the denominators. Together with (8), this gives the
full equalities
\[
 M_P=G(t)\mathbb C[t,D_1^{-1},\ldots,D_s^{-1}],\qquad
 V_P=G(t)\mathbb C(t).                                           \tag{9}
\]

For arbitrary \(c\), the exact convergent-core result holds also whenever
the total active exponent \(K_{\rm act}=\sum_{i:a_i\ne0}k_i\) is \(1\) or
\(2\). The diagonal coefficients can be distinct. Put
\[
 x=T^2=\frac{t^2}{\prod_iD_i^{k_i}},\qquad
 R=\mathbb C[t,D_1^{-1},\ldots,D_s^{-1}].
\]
Write the original period as
\(\Phi_{P_0}[A](T)=h_{A,0}(x)+T h_{A,1}(x)\) and let \(N\) be the
\(\mathbb C[x]\)-span of all these parity germs. Whenever \(h_{A,e}\ne0\),
each \(k_i e+d_i\) is even, so the relative prefactor of that parity in
(7) is the rational element
\[
 t^e\prod_iD_i^{-(k_i e+d_i)/2}\in R.
\]
Consequently \(M_P/G\subset R\,N\) inside \(\mathbb C((t))\).

Taking the even part of (6) gives \(1\in N\). The proof of §2, applied to
these parity series in \(x\), gives \(N\cap\mathbb C\{x\}=\mathbb C[x]\).
For this application use \(w_{q,0}\) as a single baseline. The factor
\(w_{q,e}/w_{q,0}\) is rational in \(q\), since each double-factorial
argument changes by a fixed even integer, and polynomial \(x\)-time
shifts contribute only further ratios of the same type.
It follows by clearing denominators that
\[
 (\mathbb C(x)N)\cap\mathbb C\{x\}[x^{-1}]=\mathbb C(x).
\]

The rational map \(x(t)\) has degree \(2\) globally and vanishes to
order \(2\) at \(t=0\). Write \(\prod_iD_i^{k_i}=1-B_1t+B_2t^2\).
The nontrivial deck transformation is
\(t'=-t/(1-B_1t)\). It is analytic at zero and fixes \(x\). Over
\(\mathbb C(x)\) use the basis \(1,t\) of \(\mathbb C(t)\).
Thus every element of \(V_P/G\) is of the form \(u(x)+tv(x)\) with
\(u,v\in\mathbb C(x)N\). If this is a convergent Laurent germ in \(t\),
composing with \(t'\), subtracting, and dividing by \(t-t'\) shows
that both \(u(x(t))\) and \(v(x(t))\) are convergent Laurent germs.
In the analytic local coordinate \(T\), their Laurent powers are even.
So \(u,v\) are convergent Laurent germs in \(x=T^2\), and the preceding
identity gives \(u,v\in\mathbb C(x)\).

For the polynomial-time lattice use the same degree equality.
The multiplication map
\[
 \mathbb C(t)\otimes_{\mathbb C(x)}\mathbb C((x))
       \longrightarrow \mathbb C((t))
\]
is injective because \(1,t\) remain independent over \(\mathbb C((x))\):
all nonzero elements of the latter field have even \(t\)-valuation.
Also \(R\) is torsion-free, hence flat over the principal ideal domain
\(\mathbb C[x]\). Localizing \(R\) along \(\mathbb C[x]\setminus\{0\}\)
gives \(\mathbb C(t)\), since it already contains the algebraic
generator \(t\) over \(\mathbb C(x)\). Tensor the intersection
\(N\cap\mathbb C(x)=\mathbb C[x]\) inside \(\mathbb C((x))\) with \(R\).
Exactness and the injective multiplication map then identify
\((R\,N)\cap\mathbb C(t)=R\) inside \(\mathbb C((t))\).
Together with (8), this proves
\[
 M_P\cap\mathbb C\{t\}=G(t)R,\qquad
 V_P\cap{\cal A}=G(t)\mathbb C(t)
       \qquad (1\le K_{\rm act}\le2).                             \tag{10}
\]
In (9) and (10), the global regular core is the line
\(G(t)\mathbb C(t)\) if \(a_0=0\), and zero otherwise. Its logarithmic
derivative is \(a_0+\sum_i a_i/D_i\). This has only logarithmic
poles as a connection form, apart from the double pole of
\(a_0\,dt\) at infinity when \(a_0\ne0\).

For general diagonal coefficients and angular \(c\), (7) still
places the untwisted, rescaled generic image in the actual
algebraic-coefficient span of the unperturbed periods after the
inverse time change. Formula (5) concerns finite time shifts.
An algebraic time-dependent multiplier can combine infinitely
many such shifts, so (5) does not establish an upper bound for
that larger convergent core. For \(K_{\rm act}>2\) the degree of the rational map exceeds its local
degree at zero, so the quadratic descent used for (10) also ceases to
apply. The general diagonal full-image exclusion remains open here.
