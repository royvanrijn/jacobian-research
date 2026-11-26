# All-order angular separation after a diffuse extension of an arbitrary factorial core

Fix any finite ordinary two-monomial exponent profile \(C\). Call its
rows the core, and allow the empty core. Append any finite set of unequal
ordinary rows \(E\), with arbitrary orientations. For these added rows
put
\[
 \tau=\sum_{i\in E}\frac{d_i^2}{M_i},\qquad
 \eta=\max_{i\in E}\frac{|d_i|}{M_i},\qquad
 d_i=b_i-a_i,\quad M_i=\max(a_i,b_i).
\]

**Theorem.** For each core there is a finite, effectively described set
\(\mathcal T_C\subset(0,\infty)\). Every element is a nonzero logarithm
of a positive algebraic number times a nonzero rational number, and
is transcendental. On any compact \(K\subset(0,\infty)\setminus\mathcal T_C\)
there is a single positive tolerance \(\eta_C(K)\) such that, if
\(\tau\in K\) and \(\eta\le\eta_C(K)\), the extended ordinary profile
has disjoint consecutive root arguments in every degree. The
tolerance is computable, for example, when \(K\) is given as a finite
union of compact rational intervals disjoint from \(\mathcal T_C\).
In particular, each positive rational \(\tau\) has such a tolerance.

There is also a computable \(a_C>0\), \(\tau_C>0\), and a single
absolute integer \(L>0\), such that
\[
             0<\tau\le\tau_C,\qquad \eta\le a_C\tau^L
\]
suffice at all degrees. The constants use the unevaluated absolute
ordinary angular cutoff. The conclusions allow arbitrary core
curvature, private rows in the core, and off-circle roots. Both orders retain
the exact factorial coefficients of the core. In characteristic zero the
extended profiles satisfy SFC(2) by the rational Bézout identity.

## An angular eliminant for the exact core

Fix \(n\ge2\) and \(N=n+1\). Write \(B_j\) for the unnormalized
order-\(j\) factorial polynomial of the core, with coefficients
\(b_{j,k}\). Its leading coefficient is positive. Form the monic polynomials
\[
 H_j(q,z)=\frac1{b_{j,j}}\sum_{k=0}^j
       b_{j,k}\,q^{(nN/j)k(j-k)}z^k,\qquad j=n,N.
\]
For \(q=e^{-\tau/(2nN)}\) these put an exact common Gaussian
strength \(\tau\) on top of the core. This also makes sense for
any pair of algebraic real polynomials of the indicated degrees
with positive endpoints and nonnegative coefficients, provided
the even polynomial is nonnegative on the real axis. In the
factorial case this last property follows by writing the even
moment as \(\mathbb E[(U+zV)^j]\) for positive monomials \(U,V\).
The empty core gives \(B_j=(1+z)^j\).

At \(\tau>0\) the real even polynomial in the deformed pair is
strictly positive. Indeed, if \(X\) is a real centered Gaussian
with variance \(\tau/j\), then
\[
 b_{j,j} H_j(e^{-\tau/(2nN)},z)
                =\mathbb E\, B_j(e^{X-\tau/2}z).
\]
For real \(z\ne0\), the argument has a density on its ray, so it
misses the finite zero set of the nonzero even polynomial with
probability one. At zero its constant coefficient is positive.

For any roots \(z_1,\ldots,z_j\) of \(H_j\), including
multiplicities, put
\[
 A_j(q,u)=\prod_{r\ne s}(u-z_r/z_s),\qquad
 D_j(q,u)=(u-1)^{\boldsymbol{1}_{j\ {\rm odd}}}A_j(q,u).
\]
These are monic polynomials in \(u\) over \(\mathbb Q[q]\) for
an ordinary core. More explicitly, with \(c_j=H_j(q,0)>0\)
constant in \(q\), the product formula for the resultant gives
\[
 A_j(q,u)=\frac{(-1)^j\,\operatorname{Res}_z(H_j(q,z),H_j(q,uz))}
                   {c_j^j(u-1)^j}.
\]
The division by the monic factor \((u-1)^j\) is exact, as the
diagonal terms of the root product provide the factor.

Every root \(z\) gives \(z/\bar z\) as a root of \(D_j\). For
a nonreal \(z\) the conjugate root occurs at a different index
and yields the required quotient. For a real \(z\) only the odd
degree is possible, where the added factor supplies \(1\).
Consequently two roots of the same argument force
\[
                    R_n(q)=\operatorname{Res}_u(D_n,D_N)=0.
\]
The eliminant may vanish for other reasons as well.

The specialization \(q=0\) makes \(H_j=z^j+c_j\) and hence
\(A_j(0,u)=((u^j-1)/(u-1))^j\). The nontrivial \(n\)-th and
\(N\)-th roots of unity are disjoint, and exactly one extra
factor \(u-1\) was inserted. In fact
\[
                   R_n(0)=j_{\mathrm{even}}^{j_{\mathrm{even}}},
 \qquad j_{\mathrm{even}}\in\{n,N\}\text{ the even integer}.
\]
To check the value,
at a nontrivial \(n\)-th root \(\xi\) one has
\((\xi^N-1)/(\xi-1)=1\). If \(n\) is even the remaining
product is \(\prod_{\xi\ne1}(\xi-1)^n=n^n\). If \(N\)
is even the added lower root \(1\) gives \(N^N\).

Thus each \(R_n\) is a nonzero rational polynomial. Any
positive \(\tau\) with \(e^{-\tau/2}\) transcendental gives
disjoint arguments for the Gaussian-augmented core in every
degree, without an eventual threshold. Indeed
\(q^{nN}=e^{-\tau/2}\), so every such \(q\) is transcendental.
The [Hermite–Lindemann theorem](https://isa-afp.org/entries/Hermite_Lindemann.html)
applies in particular to every nonzero algebraic \(\tau\).
The same conclusion holds for algebraic core coefficients in
the wider positive-even setting above, with \(R_n\) over the
algebraic numbers.

The pair of orders one and two already separates by the strict
positivity of the even reference on the real axis.

The algebraic-strength conclusion also allows independent, arbitrary
positive algebraic strengths \(C_n,C_N\). Use independent variables
\(q_n,q_N\) with exponents \(k(j-k)\) in the respective \(H_j\)'s,
then evaluate \(q_j=e^{-C_j/(2j)}\).
The same two-variable ratio resultant has the same nonzero
constant \(R_n(0,0)\). Its value is an algebraic linear combination
of exponentials with algebraic exponents
\(-aC_n/(2n)-bC_N/(2N)\), \(a,b\ge0\).
After grouping equal exponents the constant cannot cancel because
it is the unique term at exponent zero. The
[Lindemann–Weierstrass linear-independence theorem](https://isa-afp.org/entries/Hermite_Lindemann.html)
gives a nonzero value. This separates arguments for any ratio
between the two positive algebraic strengths.

This includes the cores at every algebraic point of a physical
correlated-Gaussian row interpolation. The
[correlated factorial kernel](SFC_PHYSICAL_MIRRORED_ROOT_CIRCLE.md)
is a polynomial in squared correlation at the fixed integer
exponents, and even moments are nonnegative on the real axis.
The same algebraic-strength Gaussian addition therefore gives
all-order separation at those parameter values, including
the ordinary endpoint cores.

## Uniform transfer to the added ordinary rows

Choose a computable integer \(J\ge4\) from the
[absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
For the extended ordinary profile all terminals \(N\ge J\)
are covered by that theorem, independently of the core and
the new rows. At terminal two the even ordinary moment is
strictly positive on the real axis because its two exponent
columns are distinct. Only \(3\le N<J\) are needed below.

Define \(\mathcal T_C\) as the union, over those finitely many
pairs, of \(-2nN\log q\) for real roots \(q\in(0,1)\) of \(R_n\).
It is effective by algebraic root isolation and has the stated
arithmetic restriction by Hermite–Lindemann. Here and below a
compact rational interval in \(\tau\) can be mapped into a
slightly enlarged rational interval in \(q\), refining the
exponential bounds until it remains disjoint from the isolated
roots of \(R_n\).

For a quantitative gap, normalize both endpoints to one by
a positive scalar and a positive rescaling of \(z\). The normalized
coefficients of the references are bounded by \(\binom jk\),
by log convexity of the core coefficients after stripping binomial
coefficients and by the attenuation. Their roots have moduli
between \(2^{-j}\) and \(2^j\) by the monic root estimate and
reversal. The ratios in \(D_j\) accordingly have modulus at most
\(2^{2j}\). Write \(m_j=j(j-1)+\boldsymbol{1}_{j\ {\rm odd}}\).
For any cross-degree reference roots \(z,w\), the product formula
for \(R_n\) gives
\[
 \left|\frac{z}{|z|}-\frac{w}{|w|}\right|
 \ge \frac{\left|z/\bar z-w/\bar w\right|}{2}
 \ge \frac{|R_n(q)|}{2(2^{2N+1})^{m_nm_N-1}}.                 \tag{1}
\]
The first inequality follows by factoring the difference of
squares of the unit projections. A positive lower bound \(g\)
in (1), uniform over the finite terminals and over \(K\),
can be computed by rational polynomial bounds on the compact
intervals in \(q\).

Let \(P_j\) be the endpoint-normalized ordinary extension
and \(G_j\) the endpoint-normalized reference. Suppose
\(\eta\le1/2\). Applying the
[joint factorial curvature bound](MULTIVARIATE_JOINT_CURVATURE_PROFILE_CLOSURE.md)
row by row to the added rows gives midpoint-weighted second
differences, at \(u+v+2=j\), between
\(\tau/(j(1+\eta/2))\) and \(\tau/(j(1-\eta))\).
For the lower bound use \(1/M_i\le\eta\) and \(j\ge2\).
They differ from the diffuse second difference \(\tau/j\)
by at most \(2\tau\eta/j\). Comparing with the quadratics
of constant second difference \(2\tau\eta/j\) and its
negative gives
\[
 \max_k\left|\log\frac{[z^k]P_j}{[z^k]G_j}\right|
                                  \le j\tau\eta/4 .
\]
The endpoint normalizations factor rowwise, so the exact
core cancels from this comparison. The binomial-stripped
logarithms of both coefficients are nonpositive. Since
the exponential has derivative at most one on this range,
\[
                  |P_j(z)-G_j(z)|
                          \le (j\tau\eta/4)(1+|z|)^j.
\]
At a root of \(P_j\) one has \(|z|\le2^j\). Factoring the
monic \(G_j\) there exhibits a reference root \(\zeta\) with
\[
 |z-\zeta|\le(1+2^j)(j\tau\eta/4)^{1/j}.
\]
Taking unit projections multiplies this bound by at most
\(2/|\zeta|\le2^{j+1}\), by the triangle inequality.
Thus for some computable \(B_J\) the angular chord displacement
of any such root is at most \(B_J(J\tau\eta/4)^{1/(J-1)}\)
provided \(J\tau\eta/4\le1\). Taking it strictly below \(g/3\)
proves the uniform tolerance on \(K\). This argument does not
need simple or circular reference roots.

## An entire one-sided neighborhood of zero added curvature

Write \(D=J^6\). At the pair \(n,N\), the maximum exponent of
\(q\) in \(H_j\) is at most \(nNj/4\). The Sylvester determinant
giving \(A_j\) has size \(2j\), so its \(q\)-degree, and that of
\(D_j\), is at most \(nNj^2/2\). Exact division by a monic
polynomial in \(u\) does not increase this degree. The resultant
in \(u\) is homogeneous of degrees \(m_N,m_n\) in the two
coefficient lists. Therefore
\[
 \deg_q R_n\le \frac{nN}{2}(n^2m_N+N^2m_n)
                     \le n^3N^3<D .
\]
In particular \(\mathcal T_C\) has fewer than \(J^7\) elements
regardless of the core. Factor \(R_n(q)=(q-1)^{h_n}S_n(q)\) with \(S_n(1)\ne0\).
For sufficiently small positive \(\tau\), exact rational
bounds on \(S_n\) near \(1\), together with
\(1-e^{-\tau/(2nN)}\ge\tau/(4nN)\), give
\(|R_n(e^{-\tau/(2nN)})|\ge b_n\tau^{h_n}\).
These are finite computable choices, and \(h_n<D\).
Taking \(0<\tau\le\tau_C\le1\) and the minimum over terminals
in (1) gives a bound \(b_C\tau^D\) on the reference chord gaps.
The root-displacement estimate then proves the theorem with
\(L=D(J-1)-1\), after choosing \(a_C\) sufficiently small.

For example, at any prescribed rational \(\tau>0\) one may
append equally many rows \((A-1,A)\) and \((A,A-1)\),
using a total of \(\tau A\) rows, with \(A\) a large integer
multiple of twice the denominator of \(\tau\). Their exact
added curvature is \(\tau\), their relative gap is \(1/A\),
and they leave the signed total degree difference unchanged.
The theorem actually covers every sufficiently diffuse
addition at that curvature, independently of the choice of
gaps and orientations.

For a finite family of fixed cores, such as the finitely many
partitions of a given row list into merger blocks, take the union
of the exceptional sets and the minimum of the tolerances.
The same ordinary diffuse tail then gives angular separation
in all degrees for every prescribed core in the family.
