# A computable cutoff on bounded moving private cubic windows

Consider ordinary two-monomial profiles with exactly one private
linear row \((0,1)\), all other unequal rows shared. Keep the notation
of the [moving cubic model](MULTIVARIATE_MOVING_PRIVATE_CUBIC_MODEL.md)
and write
\[
 \sigma=\sum_i d_i^2/a_i,\quad
 c_2=\tfrac12(1+\sum_i d_i^3/a_i^2),\quad
 \ell=c_2^{-1/3},\quad \epsilon_j=j^{-1/3},\quad
 T_j=\ell^2(1-\sigma)/\epsilon_j,\quad H_j=\ell\xi_j/\epsilon_j^2 .
\]
Here \(\xi_j=\log(j^\Delta t)+\sum_i d_i\log a_i\) on the local
positive-endpoint branch, so
\(\xi_{n+1}-\xi_n=\Delta\log(1+1/n)\) at common physical \(t\).

**Theorem.** Given rational \(B_T,B_H\ge1\), there is an algorithm
producing an integer \(N_0(B_T,B_H)\) with the following property.
Suppose \(3/4\le\sigma\le5/4\) and both orders \(j=N-1,N\) have
\[
                         |T_j|\le B_T,\qquad |H_j|\le B_H.
\]
If \(N\ge N_0(B_T,B_H)\), their moment polynomials cannot both
vanish at the physical \(t\). The test uses a real trapezoid
estimate for the exact lattice sum and compact interval tests
for two model residuals. It does not require model zero labels or
a numerical cutoff for the escaping temperature and field regimes.

## Uniform real lattice expansion

Put \(L=1/1024\), \(c_*=7/64\), and \(r_i=d_i/a_i\).
The bounds \(\sum_i|r_i|\le\sigma\) and \(|r_i|\le\sigma\)
follow from integral nonzero gaps. For a negative row, if
\(|d_i|=1\) then \(a_i\ge2\), and otherwise
\(|r_i|\le(5/4)/2\). Thus \(r_i\ge-5/8\), which gives
\[
 c_*\le c_2\le41/32,\qquad 1/2<\ell<3.
\]
The first lower bound follows from \(2c_2\ge1-(5/8)\sigma\).
The normalized exact summand \(f_j(q)\) at \(q=k/j\) is the
Gamma quotient in equation (1) of the moving cubic model.
Write
\[
 \Phi(q,\xi)=\xi q+
    \int_0^q\left(\log(1-z)+\sum_i d_i\log(1+r_i z)\right)dz .
\]
For \(|q|\le1/4\) one has
\[
 \Phi(q,\xi)=\xi q-(1-\sigma)q^2/2-c_2q^3/3+q^4 E(q),
 \qquad |E(q)|\le8 .
\]
Indeed \(\sum_i|d_i||r_i|^k\le(5/4)^k\) for \(k\ge1\), so the
degree-\(k\) coefficients in the logarithmic derivative are bounded
by \(2(5/4)^k/k\). The sum of their degree-at-least-three
part is bounded by \(8|q|^3\) for \(|q|\le1/4\). This also gives
\(|\Phi_q-\xi+(1-\sigma)q+c_2q^2|\le8q^3\) on the real interval.

For sufficiently large \(j\), the field bounds imply
\[
   |\sigma-1|\le c_*L/64,\qquad |\xi_j|\le c_*L^2/512,
\]
and a threshold for these inequalities follows directly from
\(|\sigma-1|\le4B_T\epsilon_j\), \(|\xi_j|\le2B_H\epsilon_j^2\).
At \(q_0=L/2\), insertion of these bounds and \(8L<c_*/12\)
gives \(\Re\Phi(q_0,\xi_j)\le-\gamma\),
\(\Re\Phi_q(q_0,\xi_j)<0\), where \(\gamma=c_*q_0^3/5\).
For \(0\le q\le1\), set
\(R(q)=(1-q)\sum_i d_i^2/(a_i+d_i q)\).
Because \(3/8\le 1+r_i q\le9/4\) and \(b_i/a_i\ge3/8\),
\[
 R'(q)=-\sum_i d_i^2 b_i/(a_i+d_i q)^2
                 \le-(2/27)\sigma\le-1/18 .
\]
For \(q\ge q_0\) this gives \(R(q)<1\). Since
\(\Phi_{qq}=(R(q)-1)/(1-q)\) for \(q<1\), the phase bound
\(\Re\Phi\le-\gamma\) persists throughout \([q_0,1]\).

For real \(q\in[0,L]\) write the exact interpolation as
\[
 f_j(q)=A_0(q)e^{R_j(q)}e^{j\Phi(q,\xi_j)},\qquad
 A_0(q)=
  \frac{\exp(\tfrac12\sum_i\log(1+r_i q))}{\sqrt{1-q}} .
\]
Binet's formula for the remainder in
\(\log\Gamma(x+1)=(x+\tfrac12)\log x-x+
 \tfrac12\log(2\pi)+{\cal E}(x)\)
gives \(|{\cal E}'(x)|\le x^{-2}\) and
\(|{\cal E}''(x)|\le2x^{-3}\). Differentiate its real integral
\(2\int_0^\infty\arctan(t/x)(e^{2\pi t}-1)^{-1}\,dt\)
and use \(2\int_0^\infty t(e^{2\pi t}-1)^{-1}\,dt<1\).
For example this last bound follows by splitting at one and using
\(e^{2\pi t}-1\ge6t\) below one and \(e^{2\pi t}-1\ge e^{6t}/2\)
above one. [DLMF, Binet's formula](https://dlmf.nist.gov/5.9.E10)

In the Stirling quotient,
\[
 R_j(q)={\cal E}(j)-{\cal E}(j(1-q))
       +\sum_i[{\cal E}(j(a_i+d_i q))-{\cal E}(ja_i)] .
\]
Its value at zero is zero. The two derivative bounds and
\(\sum_i |d_i|/a_i^2,\sum_i d_i^2/a_i^3\le\sigma\) give
\(|R_j'(q)|\le12/j,\ |R_j''(q)|\le100/j\) on \([0,L]\).
They give a bound \(10/j\) for the shared first derivative
on \([0,1]\).
At integer grid points through \(q=1\) the falling-factor part
can instead be compared to the integral of \(\log(1-q)\).
Its left Riemann-sum error is between zero and \(\log j+1\),
by telescoping upper and lower rectangle sums through \(q=1-1/j\)
and integrating the last interval separately.
Since \(|\frac12\sum_i\log(1+r_i q)|\le(4/3)\sigma\) for
\(0\le q\le1\), the summands for \(q\ge q_0\) are bounded
in absolute value by \(C j e^{-\gamma j}\) with an explicit
absolute constant, for example \(C=e^{13}\).
This proves an effective tail bound for the exact lattice sum.

Fix \(T=T_j,H=H_j,\epsilon=\epsilon_j\), and put
\(q=\ell\epsilon y,\ h_j=\epsilon^2/\ell\).
Normalize the moment divided by its constant coefficient by
a further factor \(1/(j\ell\epsilon)\). The result is
\(h_j\sum_{k=0}^j f_j(k/j)\). Truncate at \(m=\lfloor jL\rfloor\).
For \(j\ge2/L\), \(m/j\ge L/2\), so the tail bound just proved
applies. Its contribution is at most a computable multiple
of \(\epsilon^4\), for example use
\(e^{-\gamma j}\le3!/(\gamma j)^3\) against its
\(O(j^{4/3})\) prefactor. For
\(g(y)=f_j(\ell\epsilon y)\) the trapezoid estimate gives
\[
 \left| h_j\sum_{k=0}^m g(kh_j)-\int_0^{mh_j}g(y)\,dy
              -\tfrac12 h_j(g(0)+g(mh_j))\right|
            \le \tfrac18 h_j^2\int_0^{mh_j}|g''(y)|\,dy .
\]
This follows by integrating the quadratic Peano kernel on
each mesh interval and bounding its absolute value by \(h_j^2/8\).
Here \(g(0)=1\), and the right endpoint term is exponentially
small by the same lattice tail bound.

After rescaling, the exponent equals
\[
   \Psi_0(y) + \epsilon\ell^4 y^4 E(\ell\epsilon y),\qquad
                    \Psi_0(y)=-y^3/3-Ty^2/2+Hy .
\]
For \(0\le e\le\epsilon\) and \(0\le y\le mh_j\),
\[
 \Re\left(\Psi_0(y)+ e\ell^4y^4 E(\ell e y)\right)
 \le -y^3/4+B_T y^2/2+B_H y.                              \tag{1}
\]
Indeed \(|\ell e y E(\ell e y)|\le 8L<c_*/12\) and
\(\ell^3 c_2=1\).
On \(|q|\le1/4\), \(|A_0(q)|\le2\) by the log series
and \(\sum_i|r_i|\le5/4\).
For real \(q\le L\), Cauchy on radius \(1/8\) gives
\[
   |A_0^{(k)}(q)|\le 2 k!8^k,\qquad
   |E^{(k)}(q)|\le 8 k!8^k .
\]
Together with the explicit \(R_j\) bounds these give an
upper bound for \(|g''|\) by a computable fixed polynomial
in \(y\) times the exponential in (1).
For instance differentiate \(g=A_0(q)e^{R_j(q)+P(y)}\)
twice, replace \(\ell,\epsilon,|T|,|H|\) by
\(3,1,B_T,B_H\) in absolute bounds, and apply the displayed
derivative constants. This bounds the integral in the
trapezoid error uniformly.

Moreover replacing \(e^{R_j(q)}\) by one under the integral
costs at most a computable multiple of \(\epsilon^4\):
\(|R_j(q)|\le12q/j\le36\epsilon^4 y\) by \(R_j(0)=0\).
In the remaining integrand, expand
\[
                 A_0(\ell e y)
                \exp\bigl(e\ell^4y^4 E(\ell e y)\bigr)
\]
to degree three at \(e=0\). The fourth-order Taylor remainder
after multiplication by \(e^{\Psi_0(y)}\) is bounded by
\(\epsilon^4\) times a polynomial in \(y\) times (1), with
computable coefficients. To obtain this polynomial, apply
the product and chain rules four times using the above
Cauchy bounds and \(\ell\le3\). Keep the exponential intact
and bound it combined with \(e^{\Psi_0}\) by (1).
The four Taylor coefficients are profile-dependent polynomials
\(P_p(y)\), fixed for both orders, with \(P_0(y)=1\). Extending each coefficient
integral to infinity costs another computable multiple of
\(\epsilon^4\): on its added tail, \(y\ge mh_j\ge L/(6\epsilon)\),
so \(1\le(6\epsilon y/L)^4\).

We obtain, writing \({\cal Q}_j\) for the nonvanishingly
normalized exact moment,
\[
 {\cal Q}_j
  =F(T_j,H_j)+\sum_{p=1}^3\epsilon_j^p C_p(T_j,H_j)
                 +{\cal R}_j,\qquad
 |{\cal R}_j|\le K_0\epsilon_j^4.                           \tag{2}
\]
Here the \(C_p\) are integrals of the three Taylor polynomials
against the cubic density, with \(1/(2\ell)\) added to \(C_2\)
for the left trapezoid endpoint. The same polynomial moments
bound all needed fixed derivatives in \(T,H\), by multiplying
by powers of \(y\) and \(-y^2/2\).
Every constant is computable using only \(B_T,B_H,L,c_*\).
For example all moment bounds of (1) can be computed by cutting
at an integer \(Y\) with \(B_T/(2Y)+B_H/Y^2\le1/8\), bounding
on \([0,Y]\) by a finite enclosure, and using
\(e^{-y^3/8}\) beyond \(Y\).
This real-lattice expansion uses only the exact linked coefficients
and avoids comparing independent \(O(\epsilon)\) errors.

## Reduced model residuals and termination

Write \(\epsilon=n^{-1/3}\), \(\lambda=\log(1+1/n)\),
\(T=T_n\), and \(H=H_n\). The exact relation is
\[
  T_N=e^{\lambda/3}T,\qquad H_N=e^{2\lambda/3}H+v,\qquad
  \epsilon_N=e^{-\lambda/3}\epsilon,\qquad
  v=\ell N^{2/3}\Delta\lambda\in\mathbb R .
\]
If \(\Delta\ne0\), then
\(|v|\ge\frac12N^{-1/3}\ge\epsilon/3\), using
\(\lambda\ge1/N\) and \(\ell>1/2\).
For large \(n\) the two bounded fields also give
\(|v|\le2B_H+1\), since
\(|(e^{2\lambda/3}-1)H|\le1\) for a computable threshold.

Define the continuous divided difference
\[
       B(T,H,v)=\int_0^1 F_H(T,H+s v)\,ds .
\]
For \(v\ne0\), this is \((F(T,H+v)-F(T,H))/v\).
Put \(D=B_T+2B_H\), and take a rational bound \(U\ge1\) for
\(K_0\) and the values and first two derivatives of \(F,C_1,C_2,C_3\)
throughout \(|T|\le B_T+1,\ |H|\le4B_H+4\).
All these bounds follow from the moments used in (2).
Since \(0<\lambda\le1/n\), one has
\[
 |T_N-T|+|H_N-H-v|\le D\lambda,\qquad
 |\epsilon_N^p-\epsilon^p|\le\lambda\epsilon^p
                     \quad (1\le p\le3).
\]
At an exact common root, (2) gives \(|F(T,H)|\le4U\epsilon\).
If \(\Delta\ne0\), divide the exact consecutive difference by
\(v\). The model difference differs from \(B\) by at most
\(UD\lambda/|v|\le3UD\epsilon^2\).
Each correction difference before division is at most
\(U\epsilon^p(|v|+(D+1)\lambda)\).
The two remainders after division contribute at most
\(6U\epsilon^3\). Summation gives
\[
                         |B(T,H,v)|\le(12D+18)U\epsilon.
\]

For \(\Delta=0\) use the other model residual
\(A(T,H)=T F_T(T,H)+2H F_H(T,H)\).
Divide the exact consecutive difference by \(\lambda/3\).
Differentiating the field path
\(F(e^{s/3}T,e^{2s/3}H)\) twice for \(0\le s\le\lambda\)
bounds its second derivative by \(U(D^2+D)\).
Thus the divided model term differs from \(A(T,H)\) by at most
\(2U(D^2+D)\lambda\).
Each divided correction costs at most \(3U(D+1)\epsilon^p\)
by the same bound with \(v=0\). Since \(\lambda\ge1/N\),
the remainders cost at most \(12U\epsilon\).
In this case
\[
                          |A(T,H)|\le(2D^2+11D+21)U\epsilon .
\]

On the required compact rectangles, both residual sums
\[
             |F(T,H)|+|B(T,H,v)|,\qquad |F(T,H)|+|A(T,H)|
                                                               \tag{3}
\]
have positive infima. For the first, at a zero of \(F\) a
vanishing divided difference at \(v\ne0\) would give two
same-height zeros. At \(v=0\) it would give a multiple zero.
Both are excluded by the
[Riccati derivative comparison](MULTIVARIATE_HALF_LINE_CUBIC_RICCATI_POSITIVITY.md).
For the second, at any model zero \(A=F_H(2H-T H_T)\).
This is nonzero by the same Riccati result for \(T\le0\),
the [bounded temperature certificate](MULTIVARIATE_HALF_LINE_CUBIC_EXTENDED_TEMPERATURES.md)
for \(0<T\le2\), and the
[exterior and compact speed comparison](MULTIVARIATE_HALF_LINE_CUBIC_POSITIVE_SPEED_CERTIFICATE.md)
for \(T\ge2\). Real fields have positive \(F\).

The lower bounds in (3) can be computed without knowing the
zeros. Subdivide the rational parameter rectangles and evaluate
interval enclosures at cell centers by real quadrature on
\([0,Y]\) and the elementary cubic tails. For \(B\), use the
exact rational center value of \(v\) in the divided difference,
or the field derivative when \(v=0\). Bound the variation
through each cell using the integral for \(B\) and the uniform
moment bounds on derivatives. The same procedure works for \(A\).
Refine the mesh and the quadrature tolerances to zero until the
resulting lower bounds \(\mu_1,\mu_0\) for (3) are positive.
Continuity, compactness, and the proved positivity give
termination. Finally enforce
\((12D+22)U\epsilon<\mu_1\) and
\((2D^2+11D+25)U\epsilon<\mu_0\), together with the elementary
thresholds used above. These explicit tests produce
\(N_0(B_T,B_H)\) and contradict either value of \(\Delta\).
