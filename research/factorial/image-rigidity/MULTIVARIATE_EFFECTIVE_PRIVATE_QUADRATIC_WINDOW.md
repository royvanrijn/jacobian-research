# Effective exclusion on bounded private quadratic fields

Use ordinary profiles with exactly one private linear row \((0,1)\)
and shared gaps \(d_i=b_i-a_i\), and put
\[
       \sigma=\sum_i d_i^2/a_i,\qquad \mu=1-\sigma .
\]
Let \(0<\mu\le1\). Use the physical fields \(\xi_j\) of the
[moving private model](MULTIVARIATE_MOVING_PRIVATE_CUBIC_MODEL.md),
with the same local logarithm at the two orders. Write
\[
                H_j=\sqrt{j/\mu}\,\xi_j .
\]

**Theorem.** Given a rational \(B\ge1\), an algorithm produces
a constant \(K(B)\) with the following property. If \(N=n+1\),
\(n\mu^3\ge K(B)\), and both \(|H_n|,|H_N|\le B\), then the
two moment polynomials cannot share the physical root \(t\).
There is no restriction on the integer degree difference.
This makes the
[bounded quadratic escaping window](MULTIVARIATE_PRIVATE_ENDPOINT_UNIFORM_THERMAL_TAILS.md#the-bounded-quadratic-window)
effective. It applies also on fixed strictly subcritical margins.

## A real-order interpolation with a fourth-order lattice error

As in the
[effective cubic expansion](MULTIVARIATE_EFFECTIVE_BOUNDED_PRIVATE_CUBIC.md),
put \(L=1/1024\), \(q_0=L/2\), and \(r_i=d_i/a_i\).
Here \(|r_i|\le1\) and \(\sum_i|r_i|\le1\), and for negative
rows \(|r_i|\le1/2\). The latter follows from \(a_i\ge2\)
if \(d_i=-1\), and otherwise from
\(|r_i|\le\sigma/|d_i|\). The phase and leading amplitude
satisfy
\[
 \Phi(q,\xi)=\xi q-\mu q^2/2-c_2 q^3/3+q^4 E(q),\qquad
 1/4\le c_2\le1,\qquad |E(q)|\le8,\quad |A_0(q)|\le2
\]
on \(|q|\le1/4\), with both analytic functions real for real \(q\).
The explicit derivative bounds on \(0\le q\le L\) are
\(|E^{(k)}|\le8k!8^k\), \(|A_0^{(k)}|\le2k!8^k\).

For real \(j\ge1\), set
\[
 \alpha_j=(j\mu^3)^{-1/2},\qquad
 r_j=(j\mu)^{-1/2},\qquad h_j=(\mu/j)^{1/2}.
\]
They obey \(r_j=\mu\alpha_j\) and \(h_j=\mu^2\alpha_j\).
At real \(y\) with \(q=r_j y\in[0,L]\), let
\[
 g_j(y,H)=A_0(q)e^{R_j(q)}
   \exp\!\left(-y^2/2+Hy-c_2\alpha_j y^3/3+
                          \alpha_j r_j y^4 E(r_j y)\right)
\]
be the exact Gamma interpolation with physical
\(\xi=H\sqrt{\mu/j}\). The function \(R_j\) is the sum
of the four types of Stirling remainder terms in the
effective cubic expansion and has \(R_j(0)=0\).
For all \(k\ge1\) needed below, Binet's real integral
gives \(|{\cal E}^{(k)}(x)|\le k!x^{-k-1}\) for the
single Gamma remainder. In fact, write
\(\arctan(t/x)=[\log(x+it)-\log(x-it)]/(2i)\).
After \(k\) derivatives in \(x\), take the difference of
the two reciprocal powers along the vertical segment
between \(x\pm it\). Its bound is \(k! t x^{-k-1}\).
Use \(2\int_0^\infty t/(e^{2\pi t}-1)\,dt<1\).
Summing the Gamma remainders then gives computable
\(|\partial_q^k R_j|\le C_k/j\) on \(0\le q\le L\),
since \(\sum_i |d_i|^k/a_i^{k+1}\le\sum_i|r_i|\le1\).
In particular \(|R_j(q)|\le C_1 q/j\).
[Binet's formula](https://dlmf.nist.gov/5.9.E10)

At integer \(j\), a nonzero scalar times the exact
moment is \(h_j\sum_{k=0}^j g_j(kh_j,H)\), where outside
the interpolation interval \(g_j\) denotes the original
summand. For uniformly bounded \(|H|\) and sufficiently
large \(j\), \(|\xi|=|H|\sqrt{\mu/j}\) is uniformly
small. As in the effective cubic expansion,
\(\Re\Phi\le-\gamma\) on \(q_0\le q\le1\) for a fixed
absolute \(\gamma>0\). Indeed at \(q_0\) the cubic
bound gives strict loss and negative derivative,
while the quadratic term is nonpositive and
\(\Phi_{qq}=(R(q)-1)/(1-q)\le0\) thereafter.
Here \(R(0)=\sigma\le1\) and \(R(q)\) is nonincreasing.
One can use \(|\xi|\le(7/64)L^2/512\) and the \(\gamma\)
from the effective cubic expansion.
The same falling-factor Riemann estimate bounds the
exact summands at \(q\ge q_0\) by \(C j e^{-\gamma j}\).

In the local integral the perturbation of
\(-y^2/2+Hy\) is real and nonpositive:
\(|r_j y E(r_j y)|\le8L<c_2/6\).
For bounded \(H\) and \(\alpha_j\le1\) the function \(g_j\)
and its first four \(y\)-derivatives are bounded by
explicit polynomials in \(y\) times \(e^{-y^2/2+|H|y}\).
This follows by differentiating the displayed exponent
and amplitude with \(r_j,\alpha_j\le1\) and the stated
analytic derivative bounds. In particular the fourth
derivative has a computably bounded absolute integral.

On a mesh interval \([a,a+h]\) the corrected trapezoid
rule is
\[
 \frac h2(g(a)+g(a+h))-\frac{h^2}{12}(g'(a+h)-g'(a)).
\]
Its difference from the integral annihilates polynomials
of degree at most three. Taylor's integral remainder
at \(a\), in the integral and the two right endpoint
terms, gives an error of modulus at most
\(h^4\int_a^{a+h}|g^{(4)}|/3\).
Apply this on the mesh \(h_j\) through
\(m=\lfloor jL\rfloor\). The correction at the right
endpoint and the short extension to \(L/r_j\) are
exponentially small by the physical loss at \(q\ge q_0\).
This holds also for the right endpoint derivative,
whose prefactor there is at most a fixed polynomial in \(j\).
The outer exact summands are exponentially small as above.
Thus the real-order interpolation
\[
  J_j(H)=\int_0^{L/r_j} g_j(y,H)\,dy
                     +h_j/2-h_j^2 g_j'(0,H)/12
\]
differs from the normalized exact moment at integer
order by at most \(C(B)j^{-2}\).
Here \(h_j^4\le j^{-2}\), and the exponential errors
are bounded by explicit powers using
\(e^{-\gamma j}\le k!/(\gamma j)^k\) for integer \(k\).

## Order derivative and model residuals

With \(G(H)=\int_0^\infty e^{-y^2/2+Hy}\,dy\), the
following bounds are uniform on each bounded \(H\)-set:
\[
        |J_j-G|+|\partial_H(J_j-G)|+|j\partial_j J_j|
                               \le C(B)\alpha_j.            \tag{1}
\]
The comparison in value uses \(A_0(r_j y)-1=O(r_j y)\),
\(R_j(r_j y)=O(r_j y/j)\), and the nonpositive real
perturbation of the Gaussian exponent with modulus at
most \(C\alpha_j y^3\).
The field derivative inserts a factor \(y\).
For the missing Gaussian tail, \(y\ge L/r_j\) implies
\(1\le r_j y/L\), giving \(O(r_j)\) with the same moment
bounds. The two mesh corrections are \(O(h_j)\), also
after a field derivative.

For the real-order derivative at fixed \(H\), one has
\(j r_j'=-r_j/2\), \(j\alpha_j'=-\alpha_j/2\), and
\(j h_j'=-h_j/2\). The quadratic Gaussian exponent is
independent of \(j\). The derivative of its displayed
perturbation is \(O(\alpha_j)\) times a fixed polynomial
in \(y\). At fixed \(q\) we have
\(j\partial_j R_j(q)=O(q/j)\) by the mean-value theorem
on each row difference: the derivative in \(x\) of
\(x{\cal E}'(jx)\) is bounded by \(3/(j^2x^2)\).
At \(q=r_jy\) the extra \(q\)-derivative contributes
\(O(r_j y/j)\). Differentiating the leading amplitude
gives \(O(r_j y)\). The moving endpoint of the integral
contributes \(O(r_j)\), since \(y=L/r_j\) there and
\(y^2 |g_j|\) is uniformly bounded by the Gaussian majorant.
Finally
\[
      g_j'(0,H)=H+r_j[A_0'(0)+\partial_q R_j(0)].
\]
Its order derivative multiplied by \(j\) is \(O(r_j)\),
by the same Gamma derivative bounds. This proves (1).
The constants in all these statements are computable:
integrals and suprema of the finitely many Gaussian
polynomial majorants can be enclosed by a finite interval
and the tail \(e^{-y^2/4}\) beyond \(y\ge4|H|\).

At the consecutive orders,
\[
 \lambda=\log(N/n),\qquad
 H_N=e^{\lambda/2}H_n+v,\qquad
 v=\sqrt{N/\mu}\,\Delta\lambda\in\mathbb R.
\]
An exact common root makes \(J_j(H_j)=O(j^{-2})\)
at both orders. In particular \(|G(H_n)|\le C(B)\alpha_n\)
by (1). If \(\Delta\ne0\), \(|v|\ge(N\mu)^{-1/2}\)
and \(|v|\le2B+1\) for large \(n\).
Define
\[
           B_G(H,v)=\int_0^1G'(H+s v)\,ds .
\]
Divide the consecutive difference by \(v\).
At order \(n\), comparing the field difference at
\(H_n,H_n+v\) to the model costs \(O(\alpha_n)\) using
(1) for the field derivative. Changing the order at
fixed field costs \(O(\lambda\alpha_n/|v|)\) by (1).
Changing the final field by \((e^{\lambda/2}-1)H_n\)
costs \(O(\lambda/|v|)\), and \(\lambda/|v|\le\sqrt{\mu/N}\).
The exact-lattice errors divided by \(|v|\) are
\(O(n^{-3/2})\). Since \(n^{-1/2}\le\alpha_n\),
all these constants are uniform and give
\[
                         |B_G(H_n,v)|\le C(B)\alpha_n.
\]

If \(\Delta=0\), divide the consecutive difference by
\(\lambda/2\). Along the linked path
\(j=n e^s,\ H=H_ne^{s/2}\), \(0\le s\le\lambda\), the
divided change is the average of
\(2j\partial_j J_j+H\partial_H J_j\).
It differs from \(H_nG'(H_n)\) by \(O(\alpha_n+\lambda)\)
by (1) and the bounded Gaussian moments.
The two exact-lattice errors now give \(O(1/n)\),
using \(\lambda\ge1/N\). Hence
\[
                         |H_nG'(H_n)|\le C(B)\alpha_n .
\]

To finish, the compact residual sums
\(|G(H)|+|H G'(H)|\) and \(|G(H)|+|B_G(H,v)|\)
have positive infima for the stated field and drift caps.
Indeed \(G'=1+HG\) gives simplicity, and \(G\) is
positive on real fields. For the divided difference,
distinct zeros of \(G\) have distinct imaginary heights.
The [half-Gaussian proof](MULTIVARIATE_PRIVATE_ENDPOINT_UNIFORM_THERMAL_TAILS.md#the-bounded-quadratic-window)
uses the Gaussian Fourier representation to show
\(\Re G(-H)>0\) for \(\Re H>0\), so all zeros have
\(\Re H>0\). On each such horizontal line
\[
 \partial_x\log |G(-H)e^{-H^2/2}|
                       =-\Re(1/G(-H))<0 ,
 \qquad x=\Re H.
\]
At a zero of \(G\) the modulus is \(\sqrt{2\pi}\)
by the whole Gaussian integral, which proves uniqueness
at each height. At \(v=0\) the divided difference is
\(G'\), already nonzero at a zero.

Subdivision with interval quadrature for \(G\) and
its derivatives yields positive rational lower bounds
for these residuals, as in the
[bounded cubic certificate](MULTIVARIATE_EFFECTIVE_BOUNDED_PRIVATE_CUBIC.md#reduced-model-residuals-and-termination).
The computed constants in the two root inequalities and
\(\alpha_n\le K(B)^{-1/2}\) then determine \(K(B)\).
Since \(\mu\le1\), choosing \(K(B)\) large also enforces
every preliminary lower bound on \(n\).
