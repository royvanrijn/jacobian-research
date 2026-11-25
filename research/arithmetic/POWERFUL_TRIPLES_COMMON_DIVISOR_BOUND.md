# Powerful triples of large order with arbitrary common divisor

An integer is \(s\)-full if each of its nonzero prime valuations
is at least \(s\).

## Theorem. The sharp exponent with arbitrary common divisor

For all sufficiently large fixed \(s\), the
number of ordered positive \(s\)-full triples \(a+b=c\le X\),
with arbitrary common divisor, is \(O_{s,\varepsilon}(X^{1/s+\varepsilon})\).
Its exponent is sharp, since
\[
(a,b,c)=(2^s t^s,2^s t^s,2^{s+1}t^s)
\]
gives \(\gg_s X^{1/s}\) distinct triples for positive integers \(t\).
Any triple of orders \(p,q,r\) with minimum \(s\) belongs to
this count.

## Lemma. Keeping the deficient primes after division

Write \(g=\gcd(a,b,c)\) and \(u=a/g,v=b/g,w=c/g\), so \(u+v=w\)
is primitive. If a prime \(\ell\) divides \(g\), its valuation in
\(g\) is the minimum of three positive valuations at least \(s\).
Thus \(g\) is \(s\)-full. Let \(D\) be the product of primes whose
positive valuation on one of \(u,v,w\) is less than \(s\). Such a
prime must divide \(g\), as otherwise its valuation stayed
unchanged. We also have
\[
 \operatorname{rad}(uvw)\le D(uvw)^{1/s}<D w^{3/s}.
 \tag{1}
\]

For squarefree \(D\), the number of \(s\)-full \(g\le G\) divisible
by \(D\) is
\[
 \ll_{s,\alpha}G^{1/s}D^{-1+\alpha}\qquad(\alpha>0).
 \tag{2}
\]
Write \(g=\prod_{e=s}^{2s-1}z_e^e\) by the standard powerful
factorization. For a prime with valuation \(qs+j\), \(0\le j<s\),
use exponent \(q\) in \(z_s\) when \(j=0\), and otherwise use
exponent \(q-1\) in \(z_s\) together with one factor in \(z_{s+j}\).
In any representation each prime of \(D\) divides
some root \(z_e\). Assign it to one such root. Write \(D_e\) for
the product of the primes assigned to \(e\), and put \(z_e=D_e y_e\).
For a fixed assignment the number of tuples is at most
\[
 \left(\frac{G}{\prod_e D_e^e}\right)^{1/s}
      \prod_{e=s+1}^{2s-1}\left(\sum_{n\ge1} n^{-e/s}\right)
 \ll_s G^{1/s}/D .
\]
The series converge, and \(s^{\omega(D)}\ll_{s,\alpha}D^\alpha\)
bounds the number of assignments. This proves (2).
[Canonical powerful factorization](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#fourier-counting-in-a-box)

## Excluding the remaining limiting profile

Suppose the stated upper bound fails for an unbounded sequence of
orders \(s\). For each fixed such \(s\), some \(\varepsilon_s>0\)
gives failures at unbounded \(X\). Dyadically split into bins
\(g\asymp G, w\asymp Y, D\asymp H\). Let \(M_s(H,Y)\) be the count
of primitive quotients in the last two bins. Equation (2) bounds
the number of triples in the bin by
\(\ll_{s,\alpha}G^{1/s}H^{-1}X^\alpha M_s(H,Y)\), since \(D\le g\le X\).
Choosing \(\alpha\) small relative to \(\varepsilon_s\) and absorbing
the \(O((1+\log X)^3)\) bins gives bins with
\[
 M_s(H,Y)\gg_s Y^{1/s}H X^{\varepsilon_s/2},
\]
after reducing \(\varepsilon_s\) if necessary. Here \(GY\ll X\).
In particular \(Y\) tends to infinity on a subsequence. At fixed
\(s\) further select a subsequence with \(h_s=\lim \log_Y H\), which
exists since \(1\le H\ll Y^3\). Along this subsequence the quotient
count has limiting exponent at least \(1/s+h_s\).

For every fixed \(\rho>0\), (1) puts these quotients in the \(abc\)
count with height \(2Y\) and budget \(3/s+h_s+\rho\) for sufficiently
large \(Y\).
Using the \(2\lambda/3\) bound and then taking \(\rho\) down to zero
at fixed \(s\) gives \(h_s\le3/s\). Next, for each fixed \(\eta>0\)
and all sufficiently large \(s\), the limiting one-third count
applies at budgets just above \(3/s+h_s\). It gives
\[
 1/s+h_s \le (1/3+\eta)(3/s+h_s).
\]
For \(0<\eta<2/3\) this implies
\(s h_s\le3\eta/(2/3-\eta)\), and hence \(s h_s\to0\).
[The elementary and limiting small-budget bounds](ABC_EXCEPTIONAL_SET_AND_POWERFUL_TRIPLES.md#exceptional-set-estimates-and-dyadic-coordinates)

Now form anatomic boxes for the primitive quotients, with truncation
loss \(o(1/s)\). The compression index of the anatomic factorization
can be chosen at least \(s\). Thus the roots at indices \(e<s\)
group exactly the primes in \(D\), by their actual valuations in
the quotients. Their total logarithmic root mass in the box is
at most \(h_s+o(1/s)\). The total root mass \(L_s\) of a populated
box is at most \(3/s+h_s+o(1/s)\) by (1) and the anatomic radical
comparison. After paying the dyadic and truncation losses at
fixed \(s\) and then sending \(s\to\infty\), a box attaining the
quotient count has exponent at least \(1/s+h_s-o(1/s)\).
The limiting upper bound \(D_{\mathrm{box}}/L_s\le1/3+o(1)\) forces
\(D_{\mathrm{box}}/L_s\to1/3\) and \(sL_s\to3\).

Consequently the total normalized mass below degree \(s\) tends to
zero. At every other coordinate the scaled degree \(eL_s\) is at
least \(sL_s\), so all three weak-limit measures are supported
on \([3,\infty)\). Balance fixes each mass at \(1/3\), and the
numerical height bounds each first moment by one. Thus all three
are point masses at three, contradicting the endpoint lower bound
\(1/2\) on their summed nonatomic mass.
[Anatomic factorization, Lemma 2.2](sources/bernert-browning-lichtman-teravainen-2410.12234v2.pdf) · [The endpoint diffusion theorem](ABC_WIDE_COFACTOR_CONGRUENCES.md#necessary-conditions-at-normalized-exponent-13)

## Asymptotic expansion for comparable orders

Let the fixed orders be \(s_A,s_B,s_C\), with
\(m=\min s_i\), \(M=\max s_i\). There is an absolute \(\theta>0\) such that for
\(M/m\le1+\theta\) and large \(m\), the count with arbitrary
common divisor has an asymptotic expansion
\[
 A_{\boldsymbol s}(X)=
 \sum_{j=0}^J C_{\boldsymbol s,j}X^{1/(M+j)}
 +O_{\boldsymbol s,J,\varepsilon}\bigl(X^{1/(M+J+1)+\varepsilon}\bigr),
 \qquad (-1)^j C_{\boldsymbol s,j}>0 .
 \tag{PA}
\]
Here \(A_{\boldsymbol s}(X)\) counts ordered positive triples
with the three prescribed fullness orders and \(c\le X\).
Every fixed \(J\) is allowed for sufficiently large orders.
In fact, with the constant \(\kappa\) below, it suffices that
\(M\ge16(J+1)^2\) and \(M\ge100(J+1)/\kappa\), in addition
to an absolute lower cutoff on \(m\). In particular the number
of terms allowed is \(\gg\sqrt M\). The implied constant need
not be uniform in the orders or \(J\).

This makes the maximum order \(M\) govern the leading exponent
even when the orders are unequal. All three equal orders give
the count in the opening theorem. The single-integer
asymptotic underlying the calculation is classical. The proof
uses a coprimality-uniform form and a summable estimate for
primitive additive quotients. [Cohen–Davis, §§6–7](https://acta.bibl.u-szeged.hu/14284/1/math_031_fasc_003_004_363-371.pdf)

### Uniform bound for primitive quotients

For a primitive quotient \(u+v=w\), let \(P=\operatorname{rad}(uvw)\)
and let \(D\) be the squarefree product of primes whose valuation
on their quotient letter is below its required order. For
\(w\asymp Y,D\asymp H\), there are absolute
\(\kappa,\theta>0\) such that the number of these quotients,
at sufficiently large comparable fixed orders, is
\[
 \ll_{\boldsymbol s,\epsilon}
 Y^{(1-\kappa)/M+\epsilon}H^{m/M}.
 \tag{PQ}
\]
Indeed otherwise choose counterexample orders with
\(m\to\infty,\ M/m\to1,\ \kappa\to0\). At each fixed
triple of orders, failure for some fixed positive exponent
allowance gives arbitrarily large \(Y\). Pass to a subsequence
on which \(h=\lim\log_Y H\) exists, since \(1\le H\ll Y^3\). The quotient exponent is at least
\((1-\kappa)/M+(m/M)h\). Write \(\lambda_0=\sum_i1/s_i\).
The radical bound is \(\operatorname{rad}(uvw)\le D w^{\lambda_0}\).
The \(2\lambda/3\) bound, at budgets just above \(h+\lambda_0\),
first gives \(h=O(1/M)\). Then the limiting
\((1/3+\eta)\lambda\) bound gives
\[
 (m/M-\tfrac13-\eta)Mh
 \le(\tfrac13+\eta)M\lambda_0-(1-\kappa).
\]
Let the orders tend to infinity and then \(\eta\downarrow0\),
to obtain \(Mh\to0\).

Use anatomic boxes with truncation loss \(o(1/M)\) and
compression index above all the required orders. Every
coordinate on a letter below its order consists only
of its deficient primes, so their total logarithmic
root mass is at most \(h+o(1/M)\). The total root mass
is at most \(\lambda_0+h+o(1/M)\), while the chosen box
exponent is at least \(1/M-o(1/M)\). The limiting
one-third upper coefficient forces \(ML\to3\) and
normalized exponent \(d\to1/3\). The absence of
normalized low root mass then forces all three weak
measures to be supported in \([3,\infty)\).
They each have mass \(1/3\) and first moment at most
one, hence all are point masses at three. This
contradicts endpoint diffusion. All height limits and
the choice of \(h\) precede passage to varying orders,
so the contradiction gives fixed positive \(\kappa,\theta\).

### Counting the coprime part of the common divisor

Write the common divisor as \(g=f b\) with all prime factors
of \(f\) dividing \(P\) and \(\gcd(b,P)=1\). The latter
factor \(b\) is \(M\)-full. At each prime \(\ell\mid P\), let
\(K_\ell=\max_i(s_i-v_\ell(u_i))\), with
\((u_A,u_B,u_C)=(u,v,w)\). A positive exponent of \(f\)
must be at least \(K_\ell\), and exponent zero is allowed
exactly when \(\ell\nmid D\). Primitivity gives \(K_\ell\ge m\).
Thus, for any
fixed \(\sigma>0\), its allowed series satisfies
\[
 \sum_f f^{-\sigma}
 \le D^{-m\sigma}\prod_{p\mid P}(1-p^{-\sigma})^{-1}
 \ll_{\sigma,\epsilon}D^{-m\sigma}w^\epsilon.
 \tag{SF}
\]
This follows by comparing each local factor outside a
finite set with \(p^{\epsilon/3}\), and using \(P\le w^3\).

Let \(B_M(Y;P)\) count \(M\)-full integers at most \(Y\)
coprime to \(P\), and put \(a_i=M+i\) for \(0\le i\le J\).
Write \(\zeta_P(t)=\zeta(t)\prod_{p\mid P}(1-p^{-t})\).
The Dirichlet series of its indicator is
\[
 \prod_{p\nmid P}\left(1+\sum_{\ell\ge M}p^{-\ell s}\right)
   =G_M(s;P)\prod_{i=0}^J\zeta_P(a_i s),
\]
where \(G_M\) has absolutely summable Dirichlet coefficients
on each real \(\sigma>1/(M+J+1)\), uniformly in \(P\).
Indeed the local product with \(\prod_i(1-p^{-a_i s})\)
cancels all positive exponents below \(M+J+1\), for
\(M>J\). The local factors and hence \(G_M\) are positive
at each of the real points \(1/a_i\).

For completeness, if \(I\subset\{0,\ldots,J\}\) is nonempty,
put \(d=|I|\), \(A=\sum_{i\in I}a_i\). The count of
\(\prod_{i\in I}x_i^{a_i}\le Y\) with all \(x_i\) coprime
to \(P\) is
\[
 \rho_P\sum_{i\in I}Y^{1/a_i}
       \prod_{j\in I\setminus\{i\}}\zeta_P(a_j/a_i)
 +O_{M,J}\left(2^{d\omega(P)}
                   Y^{e_I}\right),
 \tag{DH}
\]
where \(\rho_P=\prod_{\ell\mid P}(1-1/\ell)\) and the
error exponent in (DH) can be taken to be
\(e_I=(d-1)(1+J/M)/A\).
For \(d=1\) the elementary coprime count gives error
\(O(2^{\omega(P)})\). For \(d\ge2\), put \(U=Y^{1/A}\).
At least one coordinate is at most \(U\), so apply
inclusion-exclusion over the nonempty sets \(K\) of these
coordinates. For proper \(K\), apply induction to the
complement \(B=I\setminus K\). The resulting error sums
converge with the required exponent since \(a_j e_B<1\)
for each \(j\in K\). Indeed, for \(b=|B|\),
\[
 a_j e_B\le(1-1/b)(1+J/M)^2<1
\]
under \(M\ge16(J+1)^2\). Summing these errors over
the small coordinates gives \(U^{|K|+e_B\sum_{i\in B}a_i}\),
at most \(U^{(d-1)(1+J/M)}\), with the stated
\(2^{d\omega(P)}\) factor.

For each main term use, at real \(t>0,\ t\ne1\),
\[
 \sum_{\substack{x\le U\\(x,P)=1}}x^{-t}
   =\zeta_P(t)+\rho_P\frac{U^{1-t}}{1-t}
       +O_t(2^{\omega(P)}U^{-t}).
\]
The unrestricted weighted sum follows by partial summation
with the zeta constant continued to \(t<1\).
Inclusion-exclusion over the divisors of \(P\) gives the
display uniformly. Extending the divisor sums from
\(\ell\le U\) to all divisors costs
\(O_t(2^{\omega(P)}U^{-t})\) in both leading terms.

For a fixed pole index \(i\), put
\(Z_j=\zeta_P(a_j/a_i)\) and
\(R_j=\rho_P U^{1-a_j/a_i}/(1-a_j/a_i)\) for \(j\ne i\).
The alternating sum over nonempty \(K\subseteq I\setminus\{i\}\)
of the products of \(Z_j+R_j\) on \(K\) and \(Z_j\) elsewhere
is \(\prod_{j\ne i}Z_j+(-1)^d\prod_{j\ne i}R_j\).
The first term gives the coefficient in (DH).
Upon multiplication by \(Y^{1/a_i}\), the second has
coefficient
\((-1)^d\rho_P^d a_i^{d-1}/\prod_{j\ne i}(a_i-a_j)\)
on \(U^d\). Lagrange interpolation gives
\(\sum_i a_i^{d-1}/\prod_{j\ne i}(a_i-a_j)=1\).
Thus these terms cancel the all-small box, whose leading
term is \((-1)^{d+1}\rho_P^d U^d\).
In any product using a remainder in the \(j\)-th slot,
the factor \(U^{-a_j/a_i}\) cancels that slot from the
initial \(U^{A/a_i}\). Each of the other \(d-1\) slots
then contributes an exponent at most \(1+J/M\).
This gives the same \(U^{(d-1)(1+J/M)}\) bound for
the remainder in the main-term calculation. This proves (DH).
Finally \(e_I<1/(M+J+1)\) for the full \(I\) under
\(M\ge16(J+1)^2\), using
\((1-1/d)(1+(J+1)/M)^2<1\).

Convolve (DH) with the absolutely summable coefficients
of \(G_M\), including main term tails beyond \(Y\).
For \(1/(M+J+1)<\sigma<1/(M+J)\), the result is
\[
 B_M(Y;P)=\sum_{j=0}^J b_j(P)Y^{1/a_j}
                 +O_{M,J,\sigma}(2^{(J+1)\omega(P)}Y^\sigma),
\]
where
\[
 b_j(P)=\rho_P G_M(1/a_j;P)
                \prod_{i\ne j}\zeta_P(a_i/a_j).
\]
Exactly \(j\) of the zeta arguments lie in \((0,1)\);
their values are negative and the other factors are
positive, so \((-1)^j b_j(P)>0\).
The coefficients are uniformly bounded in \(P\),
and the estimate extends to \(0<Y<1\), where the
actual count is zero and \(Y^{1/a_j}\le Y^\sigma\).

### Summing the quotient and common factors

The sought count is \(\sum_{u+v=w\ \mathrm{primitive}}
\sum_{f\ \mathrm{allowed}} B_M(X/(wf);P)\).
To justify summing its error and main terms, choose
\(\sigma>1/(M+J+1)\) sufficiently close to that
endpoint, and put \(h=\log_Y H\). By (SF), the quotient
summability condition for the error is bounded by
the dyadic sums of \(Y^{-\sigma+\epsilon}H^{-m\sigma}\)
times the quotient count. For \(h\le20/M\), (PQ) supplies
exponent at most
\[
 -\kappa/M+(1/M-\sigma)(1+mh)+O(\epsilon)
 \le-\kappa/M+\frac{21(J+1)}{M^2}+O(\epsilon)<0 .
\]
For \(h>20/M\), use the elementary two-smallest-radicals
count \(O_\epsilon(Y^{2\lambda_0/3+\epsilon}H^{2/3})\).
This is uniform in \(H\): the product of the two smallest
radicals is at most \((D w^{\lambda_0})^{2/3}\), and fixing
their complete values determines the third. In that branch
the dyadic exponent is at most
\[
 (2/3)\lambda_0-\sigma-(m\sigma-2/3)h+O(\epsilon)<0
\]
with a fixed margin of order \(1/M\), since
\(\lambda_0\le3.3/M\) and \(m\sigma-2/3\ge1/6\).
The latter follows from \(M/m\le1.1\),
\(M\ge16(J+1)^2\), and \(\sigma>1/(M+J+1)\).
At each fixed order take \(\epsilon\) sufficiently small.
The logarithmically many \(H\)-bins do not affect
summability. The same bound implies absolute convergence
of all main coefficients, since \(1/a_j>\sigma\).
Explicitly,
\[
 C_{\boldsymbol s,j}
 =\sum_{u+v=w\ {\rm primitive}} b_j(\operatorname{rad}(uvw))
                       \,w^{-1/a_j}\sum_{f\ {\rm allowed}} f^{-1/a_j}.
 \tag{PC}
\]
The terms have a common strict sign since every quotient
allows positive choices of \(f\), and the quotient \(1+1=2\)
already gives a nonzero term. Choosing \(\sigma\) as close
to \(1/(M+J+1)\) as needed proves (PA).

In particular, after dividing by \(X^{1/M}\), the count
coming from any fixed primitive quotient converges to its
term in \(C_{\boldsymbol s,0}\). The limiting distribution
on primitive quotients has total mass one after normalization
by that constant. Its height tail is a power saving by the
same dyadic summability argument.
