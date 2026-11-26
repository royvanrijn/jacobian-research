# Isolated two-clock exclusion at a positive exceptional Gaussian strength

Fix an ordinary factorial core \(C\) and adjacent orders \(n,N=n+1\),
\(n\ge2\). Write its original factorial moments as \(Q_j^C\).
For independently variable positive Gaussian strengths define
\[
 F_j(\tau_j,z)=\sum_{k=0}^j [z^k]Q_j^C(z)\,
                   \exp[-\tau_j k(j-k)/(2j)]\,z^k,
                   \qquad j=n,N .
\]
At a positive common strength \(\tau_*\), suppose
\(F_n(\tau_*,z)\) and \(F_N(\tau_*,z)\) have a common root in
the same physical variable. Set \(q_*=\exp[-\tau_* /(2nN)]\).

Use the two complex clock coordinates
\((X,Y)=(\exp[-\tau_n/(2n)],\exp[-\tau_N/(2N)])\).
The two polynomials have coefficients polynomial in their
respective clock coordinates. Say the contact is isolated if
\((q_*^N,q_*^n)\) is an isolated complex point of the locus where
those polynomials have a common factor of degree at least two.
This is a finite algebraic test on the Sylvester matrix.
The condition permits multiple roots and a singular contact.

Take a fixed finite list of reflected row pairs with integer
gaps \(d_i>0\) and positive integer multiplicities \(c_i\).
The pairs have midpoints \(m_i=b_i/u+a_i\), with fixed positive
rational \(b_i\), rational \(a_i\), and variable \(u>0\).
Repeat the entire list \(r\) times. For ordinary tails the endpoints
\(m_i\pm d_i/2\) and \(r\) must be nonnegative integers.
Put \(A=2\sum_i c_i d_i^2/b_i\) and \(\tau=Aru\).

**Theorem.** If the two-clock contact is isolated in the above
sense, there are positive computable \(\delta,\epsilon\) such that
every indicated ordinary reflected tail with
\(0<u<\epsilon\) and \(|\tau-\tau_*|<\delta\) has coprime original
consecutive moments at \(n,N\).
This applies in particular to repeated copies of one reflected
fixed-gap pair at any sufficiently large common midpoint.

A finite list of such exclusions is uniform over the finitely
many positive exceptional strengths of a fixed core whenever
their two-clock tests pass. The proof gives an exact transfer
through positive exceptions without a bound on contact order.

## The two independent clock directions

Let \(R_{j,k}(d,m)\) be the endpoint-normalized coefficient
multiplier of one reflected pair, with the binomial factor stripped.
Its endpoints are both one, and direct factorial cancellation gives
\[
 R_{j,k}(d,m)=
 \prod_{h=1}^{dk}
 \frac{j(m-d/2)+h}{j(m+d/2)-dk+h}.
\]
Write \(F_j^0(k)=k(k-j)/(2j)\). Expansion of this finite product's
logarithm gives, with \(x=1/m\),
\[
 \log R_{j,k}(d,m)=
 2d^2 x F_j^0(k)-\frac{d^2}{j}x^2 F_j^0(k)+O(x^3).
\]
Indeed if \(\alpha_h=-d/2+h/j\) and
\(\beta_h=d/2-dk/j+h/j\), then
\(\alpha_h-\beta_h=-d(1-k/j)\) and
\(\sum_{h=1}^{dk}(\alpha_h+\beta_h)=dk/j\).
These identities give the first and second log coefficients.
Higher coefficients are rational by the log series.

Define
\[
 f_{j,k}(u)=\frac1{Au}\sum_i c_i\log R_{j,k}(d_i,m_i).
\]
This is analytic at zero with rational Taylor coefficients.
It satisfies \(f_{j,k}(0)=F_j^0(k)\) and
\[
 f_{j,k}'(0)=-\left(\beta+\frac{\gamma}{2j}\right)F_j^0(k),
 \quad
 \beta=\frac2A\sum_i\frac{c_i d_i^2 a_i}{b_i^2},\qquad
 \gamma=\frac2A\sum_i\frac{c_i d_i^2}{b_i^2}>0 .
\]
The actual repeated tail multiplies the \(k\)-th coefficient by
\(\exp(\tau f_{j,k}(u))\). Mirror symmetry gives the same
physical argument in both orders, up to nonzero scalar factors.

Consider instead the analytic arrays with coefficients
\[
 [z^k]Q_j^C \,e^{\tau_*F_j^0(k)}
  \exp\left(v f_{j,k}(u)+t\frac{f_{j,k}(u)-F_j^0(k)}{u}\right).
                                                               \tag{1}
\]
The quotient is continued at \(u=0\). The actual tails use
\(v=\tau-\tau_*\) and \(t=\tau_*u\).
On the face \(u=0\), the arrays are Gaussian with strengths
\[
                 \tau_j=\tau_*+v-\left(\beta+\frac{\gamma}{2j}\right)t .
                                                               \tag{2}
\]
The linear map \((v,t)\mapsto(\tau_n-\tau_*,\tau_N-\tau_*)\)
is invertible, since \(\gamma>0\) and \(n\ne N\).

## Exact elimination of an accumulating contact

At the positive common strength, the even reference is strictly
positive on the real axis, by the
[lognormal averaging representation](MULTIVARIATE_DIFFUSE_PHYSICAL_EXCEPTION_POINTS.md#the-finite-set-of-physical-reference-contacts).
Hence a physical contact contributes a conjugate pair.
Also the ordinary same-strength Gaussian resultant is nonzero
as a polynomial in \(q=\exp[-\tau/(2nN)]\). At \(q=0\),
one has two endpoint binomials, which have no common root
since one degree is even and the other odd and all four
endpoint coefficients are positive. In fact common roots would
give incompatible signs for the two evaluations of \(z^{nN}\).
The polynomial has rational coefficients, so \(q_*\) is algebraic.
It follows from Hermite–Lindemann that \(\tau_*>0\) is
transcendental. [Hermite–Lindemann input](MULTIVARIATE_ARBITRARY_CORE_DIFFUSE_EXTENSION.md#an-angular-eliminant-for-the-exact-core)

For the arrays (1), let \(D(u,t,v)\) be the Sylvester
determinant and let \(M_i(u,t,v)\) be all its minors of size
\(n+N-1\). Their common zero locus is the locus where the
polynomials have a common factor of degree at least two.
This is the Sylvester kernel dimension formula, with leading
coefficients nonzero near the point in question.
The analytic functions \(D,M_i\) have algebraic Taylor
coefficients, since \(e^{\tau_* F_j^0(k)}\) is an integer
power of \(q_*\).

The function \(D(0,0,v)\) vanishes at zero to some finite order
\(b\), since it is the same-strength resultant just considered.
Choose a small circular contour in \(v\) containing only this
zero. For small \(u,t\), the contour encloses \(b\) zeros
\(v_\ell(u,t)\) of \(D\), with multiplicities. For formal weights
\(s_i\), form
\[
 {\cal E}(u,t;s)=
       \prod_{\ell=1}^b \left(\sum_i s_i M_i(u,t,v_\ell(u,t))\right).
\]
Its coefficients in \(s_i\) are analytic in \((u,t)\), even
at multiple zeros. To see this, obtain the power sums of the
product factors by integrating
\((\sum_i s_i M_i)^k\,\partial_v D/D\) around the contour
and dividing by \(2\pi i\). Newton's identities recover the
product. Expanding the contour integrals in \(u,t\) gives
residues of Laurent series in \(v\) with algebraic coefficients.
Thus every Taylor coefficient of \({\cal E}\) is algebraic.
For fixed \(u,t\) the polynomial in \(s\) vanishes identically
exactly when one contour zero makes every minor vanish:
a product of linear forms is zero as a polynomial precisely
when at least one factor is identically zero.

If there were an accumulating sequence of actual contacts
with \(u\to0\), \(\tau\to\tau_*\), the ordinary even moment
would exclude real roots. Thus each contact would make
\({\cal E}(u,\tau_*u;s)\) identically zero in \(s\).
The identity theorem makes this true along the entire germ
\(t=\tau_* u\). Each coefficient of \(u^a\) in that restriction
is a polynomial in \(\tau_*\) with algebraic coefficients,
so transcendence makes every such polynomial identically zero.
It follows that \({\cal E}(u,t;s)\) vanishes on the full germ.

In particular on \(u=0\), each small \(t\) gives a contour zero
where all the minors vanish. The zeros all approach zero as
\(t\to0\) because \(D(0,0,v)\) only vanishes at zero there.
Equation (2) now gives distinct arbitrarily close complex
two-clock contacts, contrary to isolation. This proves exclusion.

For effectivity, algebraic isolation is a real quantifier
elimination test on the real and imaginary parts of both complex
clocks and the Sylvester minors. A rational contour radius and
a complex parameter polydisc can be certified by Rouché,
using the finite leading order of \(D(0,0,v)\) and analytic
bounds from (1). The contour formula bounds \({\cal E}\) on
a smaller polydisc and computes each Taylor coefficient by
finitely many algebraic operations. Search for the first nonzero
coefficient of some weight monomial on \(t=\tau_*u\). The
isolation proof guarantees that the search terminates.
The test for coefficient zero is the algebraic polynomial test
before evaluation at \(\tau_*\). Interval enclosures at that
explicit logarithm and a Cauchy estimate for the subsequent
terms give a rational punctured disk where the coefficient stays
nonzero. This yields the claimed ordinary neighborhood.
