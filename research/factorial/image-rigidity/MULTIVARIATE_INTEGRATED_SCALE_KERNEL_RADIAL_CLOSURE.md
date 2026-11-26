# Physical radial separation by an integrated factorial scale kernel

An ordinary diffuse tail on a mirrored core can satisfy SFC(2)
uniformly in every degree even when the signed endpoint polynomial
and its rational transform both change sign. The exact physical
radial drift admits a positive kernel with a uniform linear
approximation. The odd part of the factorial coefficients admits
a separate bound that retains cancellations between unequal gaps.

## A sufficient all-degree test

Fix any ordinary mirrored core, including private reflected pairs.
The added rows are fully shared, with integer lower heights \(a_i\ge1\),
upper heights \(b_i=a_i+d_i\), and gaps \(d_i\ge1\). Set
\(\epsilon_i=1\) for increasing rows and \(-1\) for decreasing rows.
Assume \(\sum_i\epsilon_i d_i=0\). Define
\[
 \begin{split}
 T(X)&=\sum_i\epsilon_i(X^{a_i}-X^{b_i}),&
 L&=\sum_i\epsilon_i\log(b_i/a_i),\\
 R(x)&=\sum_i\epsilon_i\bigl((a_i+x)^{-2}-(b_i+x)^{-2}\bigr),&
 V&=\int_0^\infty |R(x)|\,dx .
 \end{split}
\]
For \(m_i=(a_i+b_i)/2\) also define
\[
 \begin{split}
 W(u,y)&=\sum_i\epsilon_i d_i^2 e^{-m_i u}\sinh(d_i y u),\\
 M&=\int_0^\infty\sup_{|y|\le1/2}|W(u,y)|\,du .
 \end{split}
\]
Write \(c=1/3\). Either of the following rules can supply a
positive lower-bound parameter \(D\):
\[
 \begin{array}{ll}
 D=\left|\displaystyle\sum_i\epsilon_i\log
                       \frac{b_i+c}{a_i+c}\right|,
   &\text{if \(R\) has one sign on \([0,\infty)\) and is nonzero},\\[2mm]
 D=|L|-c V,
   &\text{if \(|L|>c V\)} .
 \end{array}
\]
Use the absolute constants \(J,B,\eta_1\) from the
[exact skew test for mirrored cores](MULTIVARIATE_ORDERED_OPPOSITE_PAIR_RADIAL_CLOSURE.md#an-exact-skew-test-for-any-diffuse-tail-of-a-mirrored-core).

**Theorem.** If either positive \(D\) is available and
\[
        \max_i d_i/m_i\le\eta_1,\qquad M\le D/(2B J^3),            \tag{1}
\]
the ordinary extension satisfies SFC(2) in every degree over
every characteristic-zero field. The conclusion uses the
physical coefficient ratio. The test is uniform in the core
and in the number of tail rows.

The first sign condition is a rational-function sign test
with a positive denominator. For \(V\), split at its numerator
roots on \([0,\infty)\) and integrate on the sign intervals
using the rational antiderivative. The logarithms in \(D\) and
\(L\) are of positive rational numbers. To bound \(M\) more
simply, let \(P_d(X)=\sum_{d_i=d}\epsilon_i X^{a_i}\).
The inequality \(|e^{-du/2}\sinh(dyu)|\le du/2\)
for \(|y|\le1/2\) gives
\[
       M\le\frac12\sum_d d^3\int_0^\infty u|P_d(e^{-u})|\,du .  \tag{2}
\]
The grouped integral has the elementary finite evaluation by
root isolation in the [signed endpoint test](MULTIVARIATE_SIGNED_ENDPOINT_LAPLACE_RADIAL_CLOSURE.md).
For the ungrouped \(M\), at fixed \(u\) substitute
\(r=e^{-u/2}\), \(s=e^{uy}\). The integrand is the
maximum of the absolute value of a Laurent polynomial in \(r,s\)
over \(r\le s\le1/r\). Real polynomial inequalities therefore
give enclosures on compact subintervals of \(0<r<1\).
The bound
\(\sum_i d_i^2 e^{-a_i u}\) and continuity at zero give
effective control of the two endpoints of the integral.
This gives computable enclosures for all strict instances
of (1).

## The kernel and the two lower bounds

Fix \(N=n+1\), \(n\ge2\). Write \(q(u)=u/(e^u-1)\) and let
\(C=\log(G_N/G_n)\), where \(G_j\) is the ratio of the leading
and constant coefficients at order \(j\), raised to \(1/j\).
The [exact scale identity](MULTIVARIATE_SIGNED_ENDPOINT_LAPLACE_RADIAL_CLOSURE.md#the-physical-radial-scale)
gives
\[
 C=-\int_0^\infty [q(u/N)-q(u/n)]\,T(e^{-u})\,\frac{du}{u^2}.
\]
Introduce the piecewise linear function
\[
 F(x)=\int_0^x\left(\frac{\lfloor Nt\rfloor}{N}
                         -\frac{\lfloor nt\rfloor}{n}\right)\,dt,
                  \qquad x\ge0 .
\]
Integrating each floor as a sum of step functions gives
\(\int_0^\infty e^{-ux}F(x)\,dx=[q(u/N)-q(u/n)]/u^3\).
The identity \(R(x)=\int_0^\infty u e^{-xu}T(e^{-u})\,du\)
now gives
\[
                              C=-\int_0^\infty F(x)R(x)\,dx.    \tag{3}
\]
The exchange is absolute: \(F(x)=O(1+x)\),
\(T(e^{-u})=O(u^2)\) at zero by the signed gap balance,
and \(T(e^{-u})\) decays exponentially at infinity.

On the first period, summing at the knots \(k/N\) gives
\[
                       F(k/N)=\frac{k(k-1)}{2n^2N^2},
                              \qquad 0\le k\le N.
\]
Indeed, the integrals of the first and second floor terms
there are \(k(k-1)/(2N^2)\) and
\(k(k-1)(n-1)/(2n^2N)\), respectively. At \(k=N\) the last
step of the second floor contributes zero. Between successive
\(N\)-knots the derivative of \(F\) is nonincreasing, since
only the subtracted floor can jump there. Concavity on these
intervals gives \(F\ge0\). The two floors' difference is
periodic, and \(F(x+1)=F(x)+1/(2nN)\), so \(F\ge0\)
everywhere.

The residual \(F(x)-x/(2nN)\) is periodic and piecewise linear.
At the knots of either floor in the first period it equals
\[
 \begin{array}{ll}
 -k(N-k)/(2n^2N^2),&x=k/N,\quad 0\le k\le N,\\
 \ell(n-\ell)/(2n^2N^2),&x=\ell/n,\quad 0\le\ell\le n.
 \end{array}
\]
The first formula follows from the values of \(F(k/N)\).
For the second, at \(x=\ell/n\) the two floor integrals
are \(\ell^2/(nN)-\ell(\ell+1)/(2N^2)\) and
\(\ell(\ell-1)/(2n^2)\); the final step at \(x=1\) contributes
zero as before. Thus the absolute residual is bounded by
\(\lfloor N^2/4\rfloor/(2n^2N^2)\). In view of \(F\ge0\), this gives
\[
 \left|F(x)-\frac{x}{2nN}\right|\le \frac{c}{2nN},
                 \qquad
 F(x)\ge\frac{(x-c)_+}{2nN}.                                    \tag{4}
\]
Indeed \(\lfloor N^2/4\rfloor/(nN)\le1/3\) for \(n\ge2\):
when \(n=2\) this is an equality, and when \(n\ge3\) use
\(N/(4n)\le1/3\). The constant \(c\) in (4) is attained at
\(n=2,\ x=1/3\).
If \(R\) has one sign, use the second bound in (4) in (3).
For any shift \(c_0\ge0\),
\[
    \int_{c_0}^\infty (x-c_0)R(x)\,dx
                =\sum_i\epsilon_i\log\frac{b_i+c_0}{a_i+c_0}.
\]
For each summand this follows by integrating by parts and
then integrating \((a_i+x)^{-1}-(b_i+x)^{-1}\).
The boundary term vanishes by the paired difference at
infinity. With \(c_0=c\) this yields \(|C|\ge D/(2nN)\)
for the first rule. Using the first bound in (4) and taking
\(c_0=0\) in the integral identity gives
\[
 |C|\ge |L|/(2nN)-cV/(2nN)=(|L|-c V)/(2nN).
\]
This is the second rule and allows \(R\) to change sign.

Under the first rule, any one-sign nonzero \(T\) is already
sufficient for the sign of \(R\), by its Laplace integral.
The converse fails. For example the gap-one block with
signed multiplicities \(-1,3,-2\) at heights \(H,H+1,H+2\)
has \(T=X^H(1-X)^2(2X-1)\). Write \(f(s)=s^{-2}\) and
\(\Delta f(s)=f(s+1)-f(s)\). For \(s=H+x\), its \(R(x)\)
is \(2\Delta^2 f(s+1)-\Delta^2 f(s)\). In the double integral
of \(f''\) over the unit square this is positive when
\((1+1/H)^4<2\), for instance \(H\ge6\).
The one-sign centered-kernel criterion therefore extends
the one-sign endpoint condition to actual integer tails.

## The odd coefficient and the physical exclusion

At fixed order \(j\), let \(J_j^o(v)\) denote the odd part of
the normalized log coefficient multiplier at \(v=k/j\),
continued to \(0\le v\le1\). From the trigamma integral,
as in the [exact coefficient identity](MULTIVARIATE_SIGNED_ENDPOINT_LAPLACE_RADIAL_CLOSURE.md#the-coefficient-asymmetry),
\[
        (J_j^o)''(v)=-j\int_0^\infty q(u/j) W(u,v-1/2)\,du.
\]
This retains the sum over different gaps before taking an
absolute value. Since \(q\le1\) and \(J_j^o\) has zero
endpoint values, the Dirichlet Green kernel, of maximum
integral \(1/8\), gives \(\delta_j\le jM/8\) for the skew
norm in the exact skew test.

For any possible collision with \(3\le N<J\), the latter
test gives \(|C|\le2B(\delta_n+\delta_N)\). Under (1),
\[
       2B(\delta_n+\delta_N)\le BJ M/2
                  \le D/(4J^2)<D/(2nN)\le |C|.
\]
The contradiction gives separation below \(J\). The cited
test covers the other orders by the absolute eventual
angular cutoff and terminal-two positivity. As the supplied
polynomials are rational, Bézout over \(\mathbb Q\) gives
the characteristic-zero assertion.

## A family in which both signs change

For an integer \(H\ge7\), take a unit-gap block with signed
lower-height polynomial
\[
             P_1(X)=X^H(1-X)[(H+3)X-H].
\]
This is realized by the signed multiplicities
\(-H,2H+3,-(H+3)\) at \(H,H+1,H+2\). They sum to zero.
The endpoint polynomial \(T=(1-X)P_1\) changes sign.
Its rational transform is
\[
 R(x)=(H+3)\Delta^2 f(H+x+1)-H\Delta^2 f(H+x),
 \qquad f(s)=s^{-2}.
\]
At \(x=0\) the integrand in the double-difference representation
of \(R\) is negative. For \(0\le t\le2\),
\((1+1/(H+t))^4\ge1+4/(H+2)>1+3/H\).
As \(x\to\infty\) the integrand is eventually positive uniformly
in \(t\), since the prefactors differ by three. Therefore \(R\)
changes sign on the physical transform domain.

From \(|(H+3)e^{-u}-H|\le3+Hu\) and \(1-e^{-u}\le u\), (2) gives
\(M\le6/H^3\). Also
\[
 V\le\int_0^\infty |T(e^{-u})|\,du\le12/H^3,
\]
where the first bound follows from integrating
\(u e^{-xu}|T(e^{-u})|\) in \(x\). Finally,
\[
 H^2 L=H^2\int_0^\infty T(e^{-u})\,\frac{du}{u}
              \longrightarrow \int_0^\infty e^{-w}w(3-w)\,dw=1
\]
by \(w=Hu\) and dominated convergence, with majorant
\(e^{-w}w(3+w)\). Thus \(D=|L|-c V\) is positive and of
order \(H^{-2}\), while \(M\) is of order at most \(H^{-3}\).
The theorem holds for this family at all sufficiently large
\(H\), uniformly in the core and any integer multiplicity of
the whole block. A cutoff can be obtained by bounding the
exponential remainder in the last integral, or by the
explicit bound below. If \(z=w/H\), then
\(0\le w-H(1-e^{-z})\le w^2/(2H)\). Hence
\[
 \begin{split}
 \left|H^2(1-e^{-w/H})^2-w^2\right|&\le w^3/H,\\
 \left|(H+3)e^{-w/H}-H-(3-w)\right|&\le(3w+w^2/2)/H .
 \end{split}
\]
Inserting these two bounds gives \(|H^2L-1|\le21/H\).
Every quantity in (1) can therefore be bounded uniformly
for the whole integer ray.

## Common translations with cancellations across gaps

Start with any fixed finite template of rows as above, now
allowing offsets \(a_i=0\). Add a shared integer height \(H\)
to both offsets in every row, with the original orientations
and gaps. The resulting rows are fully shared for \(H\ge1\).
Suppose the template is gap-balanced and \(T\ne0\), and let
\(p\ge2\) be the order of vanishing of \(T(e^{-u})\) at zero.
Choose the orientation so
\[
             T(e^{-u})=t_p u^p+O(u^{p+1}),\qquad t_p>0.
\]
Assume the Taylor coefficients of template \(W(u,y)\) through
order \(p-1\) vanish identically in \(y\). These are finite
rational identities in the integer offsets and gaps. They
can mix different gap groups.

Then every sufficiently large \(H\) satisfies (1), uniformly
in the mirrored core and the multiplicity of the entire
template block. Indeed, Taylor's formula on \(u\ge0\) gives
\[
 |T(e^{-u})-t_p u^p|\le C_T u^{p+1},\qquad
 \sup_{|y|\le1/2}|W(u,y)|\le C_W u^p ,
\]
with \(C_T=\sum_i(a_i^{p+1}+b_i^{p+1})/(p+1)!\) and
\(C_W=\sum_i d_i^2 b_i^p/p!\). The derivative bounds use
only that all relevant exponents in the finite exponential
sums are nonnegative. The translated transforms have the
additional factor \(e^{-Hu}\). Integrating the bound for
\(T\) with \(u e^{-(H+x)u}\) yields
\[
 R_H(x)\ge \frac{t_p(p+1)!}{(H+x)^{p+2}}
              -\frac{C_T(p+2)!}{(H+x)^{p+3}}.
\]
For \(H\ge2C_T(p+2)/t_p\), the first term supplies at least
half its own coefficient as a lower bound. Thus the first
rule applies, and integration against \((x-c)_+\) gives
\[
          D_H\ge \frac{t_p(p-1)!}{2(H+c)^p},\qquad
          M_H\le \frac{C_W p!}{H^{p+1}} .
\]
Both the relative-gap and the quantitative condition in
(1) hold for all sufficiently large integer \(H\).

For \(p=2\) the only nontrivial odd cancellation is
\(\sum_i\epsilon_i d_i^3=0\). An explicit template uses five
increasing gap-one rows of lower offset \(1\), one increasing
gap-three row of lower offset \(1\), one decreasing gap-two
row of lower offset \(0\), and three decreasing gap-two rows
of lower offset \(2\). Its signed gap and signed gap cube
are zero:
\[
                5+3-8=0,\qquad 5+27-32=0.
\]
Its signed endpoint square difference
\(\sum_i\epsilon_i(b_i^2-a_i^2)=-10\), so \(p=2\) and \(t_p=5\).
The endpoint polynomial is negative near \(X=0\) and
positive near \(X=1\) before the common translation. Its
three gap groups are individually unbalanced. The right side
of (2) has order \(H^{-2}\) for this template, the same
order as \(D_H\), whereas the exact joint bound gives
\(M_H=O(H^{-3})\). The cancellation is in the original
factorial coefficient symbol across the unequal gaps.
