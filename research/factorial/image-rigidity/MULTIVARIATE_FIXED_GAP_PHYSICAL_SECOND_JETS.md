# Physical second-order jets for fixed-gap height tails

Fix adjacent positive orders \(n,N=n+1\) and a fixed ordinary
binomial core \(C\). Its moment polynomial is
\[
 P_j^C(z)=\sum_{k=0}^j\binom jk
      \prod_{(a,b)\in C}((j-k)a+kb)!\,z^k .
\]
Rows in the core are counted with multiplicity. Append rows
\((A_i,A_i+d)\) with fixed common gap \(d\ge1\), positive
heights \(A_i\), and positive multiplicities \(r_i\). Write
\(R_0=\sum_i r_i\), \(D=dR_0\), and
\(B_j(z)=P_j^C(j^Dz)\). These are the two limiting
polynomials in a branch of the
[height reduction](MULTIVARIATE_THREE_HEIGHT_EFFECTIVE_DECISION.md).
The result below gives an effective second-order test at their
simple common roots. It uses both coordinates of the logarithmic
root displacement. The
[endpoint-normalized weak-tail test](MULTIVARIATE_WEAK_TAIL_SECOND_ORDER_CONTACT_GATE.md)
uses the angular coordinate at second order for variable signed
gaps.

## Physical response at a simple contact

Set \(\varepsilon_i=A_i^{-1}\), \(H_1=\sum_i r_i\varepsilon_i\)
and \(H_2=\sum_i r_i\varepsilon_i^2\). Dividing the tail
constant terms out of the full polynomials and scaling \(z\)
by the common factor \(\prod_i A_i^{-d r_i}\) gives, exactly,
\[
 F_j(\varepsilon,z)=\sum_{k=0}^j[z^k]B_j\,z^k
      \prod_i\prod_{h=1}^{dk}(1+h\varepsilon_i/j)^{r_i}.
                                                               \tag{1}
\]
Define \(\theta=z\partial_z\) and two polynomials in the index:
\[
 U_j(k)=\frac{dk(dk+1)}{2j},\qquad
 V_j(k)=\frac{dk(dk+1)(2dk+1)}{12j^2}.
\]
If \(\alpha\) is a simple common root of \(B_n,B_N\), then
\(\alpha\ne0\) by the constant terms. For the analytic root of
\(F_j\) near \(\alpha\), let \(y_j=\log(z_j/\alpha)\) with
\(y_j=0\) at \(\varepsilon=0\). Its expansion is
\[
 y_j=H_1a_j+H_1^2b_j+H_2c_j+O(H_1^3),                         \tag{2}
\]
where all evaluations below are at \(\alpha\):
\[
 a_j=-\frac{U_j(\theta)B_j}{\theta B_j},\qquad
 b_j=-\frac{(a_j\theta+U_j(\theta))^2B_j}{2\theta B_j},\qquad
 c_j=\frac{V_j(\theta)B_j}{\theta B_j}.                         \tag{3}
\]
Indeed, the logarithm of each coefficient multiplier in (1) is
\(H_1U_j(k)-H_2V_j(k)+O(H_1^3)\). In the multivariate Taylor
series at \(\varepsilon=0\), the omitted terms have total degree
at least three. The root is analytic by the simple-root implicit
function theorem. Substituting \(z=\alpha e^{y_j}\) applies
the commuting operator \(e^{y_j\theta}\) to (1) evaluated at
\(\alpha\). Its linear term in \(y_j\) is
\(y_j\theta B_j(\alpha)\) with nonzero coefficient. Comparing
the terms of degrees one and two gives (2)–(3). On the
nonnegative orthant, all coordinates satisfy
\(0\le\varepsilon_i\le H_1\), so the analytic remainder is
\(O(H_1^3)\) uniformly as the heights diverge.

Write \(\Delta a=a_n-a_N\), and similarly for \(b,c\). Set
\(\tau=H_2/H_1^2\) and \(r_{\min}=\min_i r_i\). Cauchy–Schwarz
gives \(1/R_0\le\tau\), and
\[
 \sum_i r_i\varepsilon_i^2
   = \sum_i (r_i\varepsilon_i)^2/r_i \le H_1^2/r_{\min}.
\]
Thus \(\tau\in I=[1/R_0,1/r_{\min}]\). The necessary condition for a
sequence of collisions tending to \(\alpha\) is that
\[
             \Delta a=0,\qquad
             \Delta b+\tau_*\Delta c=0
             \quad\text{for some }\tau_*\in I.              \tag{4}
\]
This follows from (2) by division first by \(H_1\) and then
by \(H_1^2\), passing to a subsequence with convergent
\(\tau\) for the second division. At a contact in an imaginary
quadratic field, if \(\Delta c\ne0\) and the ratio
\(-\Delta b/\Delta c\) is real, that ratio is rational since
a nonreal quadratic field has real subfield \(\mathbb Q\).

If every common root of \(B_n,B_N\) is simple in each
polynomial and fails (4), the test gives an effective bound on
the minimum height in any collision on this branch. To see
this, the polynomial \(\operatorname{Res}_z(F_n,F_N)\) has
a root product in a neighborhood of \(\varepsilon=0\).
The root pairs that limit to distinct points give a factor
tending to a nonzero constant, even if some nonmatching roots
are repeated. At each matching point \(\alpha\), the root
difference is \(\alpha(y_n-y_N)(1+O(H_1))\). Hence the first
nonzero homogeneous layer of the resultant has total degree
\[
 e=\#\{\alpha:\Delta a\ne0\}+2\#\{\alpha:\Delta a=0\}.
\]
On the nonnegative orthant it is \(H_1^e\) times a fixed
nonzero constant and the product of
\(\Delta b+\tau\Delta c\) over the \(\Delta a=0\) contacts.
The leading coefficient in this description follows, for
example, by substituting \(\varepsilon=t s\) for generic \(s\)
with all entries positive and comparing powers of \(t\).
All lower layers vanish on the open positive orthant and
therefore vanish as polynomials. The product is bounded away
from zero on \(I\), and its minimum absolute value can be
bounded below by positive rational arithmetic using algebraic
root isolation. The remainder is bounded by
\(H_1^{e+1}\) times the sum of the absolute values of the
higher-layer coefficients for \(H_1\le1\). Computing the
resultant and these bounds gives a positive rational cutoff
for \(H_1\), and \(H_1\le R_0/\min_i A_i\) converts it into
a height bound. For the case with no common roots the
resultant has a nonzero constant term and the same argument
has \(e=0\).

## Two retained unit-gap rows

Take \(d=1\) and retained rows \((A,A+1),(B,B+1)\) with
\(A,B\ge0\), each occurring once. The two rows may coincide.
Up to a constant factor the comparator is
\[
 B_j(z)=\sum_{k=0}^j\binom jk(jA+1)_k(jB+1)_k(j^Dz)^k .
\]
At a simple stationary contact \(\alpha\), meaning
\(\Delta a=0\), introduce
\[
  R=\frac{\theta^2 B_j(\alpha)/\theta B_j(\alpha)+1}{j},
  \qquad Y=A+B+R,\qquad X=(A+R-1)(B+R-1),
  \qquad d_i=n^{-D-i}-N^{-D-i}.
\]
The value \(R=-2a_j\) is independent of the comparator.
Then the second-order differences are
\[
 \Delta b=\frac{YX}{8}-\frac{(Y-1)d_1+2d_2}{8\alpha},
 \qquad
 \Delta c=\frac{2-R}{12nN}-\frac{d_2}{6\alpha}.                \tag{5}
\]

For the derivation, the coefficient ratio gives
\[
 \theta B_j=zj^D(j-\theta)(jA+1+\theta)(jB+1+\theta)B_j .
\]
Put \(s_l=\theta^l B_j(\alpha)/\theta B_j(\alpha)\) and
\(\kappa=j^{-D}/\alpha\). The right-hand operator polynomial
has coefficients \(-1,C,E,F\) in descending order, with
\[
 C=j(1-A-B)-2,\quad E=j^2(A+B-AB)+j(2-A-B)-1,\quad
 F=j(jA+1)(jB+1).
\]
At the root its value is \(\kappa\theta B_j(\alpha)\).
Applying \(\theta\) to the operator equation with
\(j^{-D}z^{-1}\theta B_j\) on one side gives
\[
 s_2=jR-1,\qquad s_3=Cs_2+E-\kappa,\qquad
 s_4=Cs_3+(E-\kappa)s_2+F+\kappa.
\]
Formula (3) now reads
\[
 b_j=-\frac{s_4-2s_2s_3+s_2^3}{8j^2},\qquad
 c_j=\frac{2s_3+3s_2+1}{12j^2}.
\]
Substitution gives the expansions
\[
\begin{split}
 s_4-2s_2s_3+s_2^3 &=j^3YX+j^2X+j\kappa(Y-1)+2\kappa,\\
 2s_3+3s_2+1 &=2j^2[(1-A-B)R+A+B-AB]+j(2-R)-2\kappa .
\end{split}
\]
Taking differences proves (5).

If \(n\ge5\), \(\alpha\) is nonreal and \(R\) is real, this
contact fails (4) for every \(\tau\in I\). Indeed, the
purported second tie gives
\[
 \alpha\left(\frac{YX}{8}+\frac{\tau(2-R)}{12nN}\right)
   =\frac{(Y-1)d_1}{8}+\frac{(3+2\tau)d_2}{12}.
\]
Both brackets outside \(\alpha\) are real, so each vanishes.
The right-hand side gives
\(Y=1-2(3+2\tau)d_2/(3d_1)\). The ratio \(d_2/d_1\) is at
most \((D+2)/((D+1)n)\) and is strictly less than that bound:
write numerator and denominator as the integrals of
\((D+2)t^{-D-3}\) and \((D+1)t^{-D-2}\) over \([n,N]\).
Consequently \(0<Y<1\) for \(n\ge5\), since \(D\ge1\)
and \(\tau\le1\). Both factors in \(X\) are strictly negative,
and \(R\le Y<1\). The left bracket is strictly positive, a
contradiction. Reality of the stationary velocity is an
additional hypothesis. The argument gives no proof that it
holds for all contacts of these two retained rows. Each
specific comparator pair can instead use (5) directly in
the algebraic test (4).

## Binomial Newton face tests

Fix a prime \(p\ge5\) with \(p\nmid nN\) and an integer
\(a\ge1\), and embed the root field of \(\alpha\) into a
finite extension of \(\mathbb Q_p\), with \(v_p(p)=1\).
For \(j=n,N\), replace \(B_j(z)\) by \(B_j(\alpha x)\)
and divide by any nonzero scalar having the minimum valuation
of its coefficients. The two normalized polynomials vanish
at \(x=1\). Suppose each has precisely two terms modulo
elements of valuation at least \(a\), supported at
\(\ell_j<u_j\), both with unit coefficients. Suppose also
that \(p\nmid w_j=u_j-\ell_j\). This gives a unit Euler
derivative at \(x=1\), because the two surviving coefficients
are opposite modulo elements of valuation at least \(a\).
The jet formulas (3) reduce in the same quotient and agree
with those of two exact binomials whose two coefficients
are opposite. In fact, division by the Euler derivative
in (3) is division by a unit, and all integer denominators
in (3) are units by the assumptions on \(p\).

For one such exact binomial the displacement is the
negative log ratio of the two coefficient multipliers in
(1), divided by \(w_j\). It follows that
\[
 a_j\equiv-\frac{d[d(\ell_j+u_j)+1]}{2j},\qquad b_j\equiv0
 \pmod{v_p\ge a}.                                             \tag{6}
\]
If \(p\nmid d\) and \(p^a>2dnN\), stationarity thus forces
the exact integer equality
\[
        N[d(\ell_n+u_n)+1]=n[d(\ell_N+u_N)+1].                 \tag{7}
\]
Indeed, each side is positive and at most \(2dnN\), and
their difference is divisible by \(p^a\). By coprimality of
\(n,N\), there is a common integer \(m\) such that
\(d(\ell_j+u_j)+1=jm\). If \(d\ge2\), the two
congruences \(nm\equiv Nm\equiv1\pmod d\) imply both
\(m\equiv0\pmod d\) and \(nm\equiv1\pmod d\), which
is impossible. Thus one qualifying prime suffices to
exclude stationarity for a fixed common gap \(d\ge2\).

For \(d=1\), (7) permits \(m=1,2\). The first case has
\(\ell_j+u_j+1=j\) at both orders. The second has
\(\ell_j=j-1,\ u_j=j\) at both orders. For \(d=1\)
the exact binomial formula also gives
\[
 c_j \equiv \frac{3(\ell_j+u_j+1)^2+w_j^2-1}{24j^2}
                  \pmod{v_p\ge a},                           \tag{8}
\]
since \(c_j=(V_j(u_j)-V_j(\ell_j))/w_j\) for
an exact binomial. In the first case \(\Delta c\)
reduces to the rational number
\[
 T=\frac1{24}\left(\frac{w_n^2-1}{n^2}
                        -\frac{w_N^2-1}{N^2}\right)\ne0 .
\]
If it vanished, \(n^2\) would divide \(w_n^2-1\) and
\(N^2\) would divide \(w_N^2-1\). Both widths are at most
their respective orders minus one, so both widths would
be one. This contradicts \(w_j\equiv j-1\pmod2\) at
adjacent orders. If \(T\) is a \(p\)-adic unit, an exact
second tie with rational \(\tau_*\) in (4) yields
\(v_p(\tau_*)\ge a\), since \(v_p(\Delta b)\ge a\) and
\(\Delta c\) is a unit. Consequently, over any finite
collection of distinct primes meeting the centered case
and unit condition at the same contact, the numerator
of a rational \(\tau_*\) in lowest terms is divisible by
the product of their \(p^a\). Its denominator is at
least as large because \(0<\tau_*\le1\).
These constraints apply to the coefficient in the
limiting jet test. A weak-tail sequence can approach such
a rational value through different finite height ratios.

In the final width-one case, the two surviving positions
have the exact multiplier ratio
\(\prod_i(1+\varepsilon_i)^{r_i}\) at either order.
Thus the binomial faces alone yield identical displacement
series at the two orders. The deeper coefficients govern
any further exclusion in this case. This section gives
conditional prime tests, and existence of a prime meeting
their hypotheses for arbitrary comparators remains open.
