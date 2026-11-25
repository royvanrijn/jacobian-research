# Entropy decrement for power-rough twin pairs

Fix \(0<\eta<1/2\) and write \(y=X^\eta\). Let \(P_y\) be the product of odd primes at most \(y\), and set
\[
 f_y(n)=\mu(n)\mathbf1_{n\ {\rm odd}}\mathbf1_{(n,P_y)=1}.
\]
Thus \(f_y\) is supported on odd squarefree integers with no prime factor at most \(y\). For a real weight \(|G(a,b)|\le 1\), put
\[
 C(G)=\sum_{X<n\le2X}\frac{f_y(n)f_y(n+2)G(n,n+2)}n,
 \qquad
 U(G)=\sum_{X<n\le2X}\frac{|f_y(n)f_y(n+2)|G(n,n+2)}n.
\]
Write \(G_3(a,b)=\mathbf1_{\omega(a)+\omega(b)\le3}\) and
\(C_2(X)=\sum_{X<r\le2X,\ r,r+2\ {\rm prime}}1/r\).
On the support of \(G_3\) in these sums, the only configurations are two primes or a prime and a squarefree semiprime. Their parity products are \(+1\) and \(-1\), respectively. Hence, for large \(X\),
\[
 C(G_3)+U(G_3)=2C_2(X).                                           \tag{1}
\]
A bound \(C(G_3)\ge -U(G_3)+c/(\log X)^2\) with a fixed \(c>0\) would give \(\gg X/(\log X)^2\) twin primes in \((X,2X]\). The factor-count restriction matters. Two squarefree semiprimes also have parity product \(+1\).

## Divisibility data and opposite prime counts

The classical Selberg example uses \(1+\lambda(n)\) to weight integers
of even prime-factor parity. Here the weights also match the local
density \(1/\varphi(d)\) for odd divisors of \(p+2\). A correction at
\(3\) and \(5\) makes the global signed sums agree as well.
[Ford, p. 1664, for the classical example and the twin-prime local density](https://www.ford126.web.illinois.edu/wwwpapers/asymsieve.pdf)

Put \(B_X=X/(\log X)^2\), \(L_X=\log\log X\), and write
\(\mathfrak S_2=2\prod_{q>2}(1-(q-1)^{-2})\), with \(q\) prime.
For odd \(n\), set
\[
 f(n)=\prod_{q\mid n}\frac{q-1}{q-2},
 \qquad C=\frac12\prod_{q>2}\left(1+\frac1{q(q-2)}\right)
         =\mathfrak S_2^{-1}.
\]
Define \(\chi(n)=1\) on \((n,15)=1\), \(\chi(n)=-9/10\) on
\((n,15)>1,\ 9\nmid n,\ 25\nmid n\), and \(\chi(n)=0\) otherwise.
For \(n\in I_X=(X+2,2X+2]\) odd, use
\[
 w_X^\pm(n)=\frac{f(n)}{C\log(n-2)}
                      \bigl(1\pm\chi(n)\lambda(n)\bigr),
\]
where \(\lambda(n)=(-1)^{\Omega(n)}\). Extend both weights by zero
outside the stated support.

**Theorem.** For all sufficiently large \(X\), the weights \(w_X^\pm\)
lie in \([0,1]\). For every fixed \(K>0\), there is \(A_K\) such that,
uniformly for odd \(d\le X\exp(-A_K L_X^2)\) and
\(X+2\le U\le2X+2\),
\[
 \begin{split}
 \sum_{n\le U,\ d\mid n} w_X^\pm(n)
     &=\frac{\operatorname{Li}(U-2)-\operatorname{Li}(X)}{\varphi(d)}
                        +O_K\!\left(\frac{f(d)X}{d(\log X)^K}\right),\\
 \sum_{n\le U,\ d\mid n}\mu^2(n)w_X^\pm(n)
     &=(\operatorname{Li}(U-2)-\operatorname{Li}(X))\,g_{\rm sf}(d)
                        +O_K\!\left(\frac{f(d)X}{d(\log X)^K}\right),
 \end{split}                                                     \tag{P1}
\]
where \(g_{\rm sf}(d)=0\) if \(d\) is not squarefree and otherwise
\[
 g_{\rm sf}(d)=
 \left(\prod_{q>2}\left(1-\frac1{q(q-1)}\right)\right)
        \prod_{q\mid d}\frac{q-1}{q(q-1)-1}.
\]
Both \(\sum_{n\le U}w_X^\pm(n)\lambda(n)\) and
\(\sum_{n\le U}w_X^\pm(n)\mu(n)\) are \(O_K(X/(\log X)^K)\)
uniformly over these prefixes. The plus model gives zero weight to
every prime in \(I_X\). The minus model has prime weight
\((2\mathfrak S_2+o(1))B_X\).

For the unfiltered sums, write \(n=dm\) and expand
\(f(dm)=f(d)\sum_{e\mid m,\,(e,d)=1}h(e)\), where \(h\) is supported
on odd squarefree integers and \(h(q)=1/(q-2)\). Put \(Y=X/d\).
Counting odd multiples of \(e\) gives the unsigned main term
\[
 f(d)\frac{U-X-2}{2d}
     \prod_{\substack{q>2\\q\nmid d}}(1+h(q)/q)
 = C\frac{U-X-2}{\varphi(d)} .
\]
Here \(f(q)/(1+h(q)/q)=q/(q-1)\). The floor and product-tail error is
\(O(f(d)Y^{3/4})\) since \(\sum_e h(e)/e^{3/4}<\infty\).
With an inserted \(\lambda(n)\), complete multiplicativity replaces
the inner count by an odd Liouville sum. For \(e\le\sqrt Y\) apply
the classical estimate \(O(t\exp(-c\sqrt{\log t}))\) for the
summatory Liouville function, also valid on odd integers by separating
even terms. For \(e>\sqrt Y\) bound absolutely using the convergence
at \(3/4\). Together these give
\(O(f(d)Y\exp(-c'\sqrt{\log Y}))\) uniformly in the prefix. The
corresponding bounds for \(\mu\) and \(\lambda\) follow from the
classical zero-free region.
[Montgomery and Vaughan, (6.17) and Exercise 6.2.1.11](sources/montgomery-vaughan-prime-number-theorem-ch6.pdf)

For the squarefree sums, \(d\) must be squarefree and
\((m,2d)=1\). For the two sequences
\(\mu^2(m)f(m)\mathbf1_{(m,2d)=1}\) and
\(\mu(m)f(m)\mathbf1_{(m,2d)=1}\), factor their Dirichlet series as
\(\zeta(s)J_d(s)\) and \(H_d(s)/\zeta(s)\), respectively. At
\(q\nmid2d\) the factors of \(J_d,H_d\) are
\[
 (1+f(q)q^{-s})(1-q^{-s}),\qquad
 \frac{1-f(q)q^{-s}}{1-q^{-s}} .
\]
At \(q\mid2d\) they are \(1-q^{-s}\) and \((1-q^{-s})^{-1}\).
At \(\sigma=1-1/L_X\) the absolute coefficient norms of the two
correction series are \(L_X^{O(1)}\), uniformly in \(d\le X\).
The factors away from \(2d\) have bounded absolute product for
\(\sigma\ge3/4\). For the remaining factors,
\(\sum_{q\mid 2d}q^{-\sigma}=O(\log L_X)\) by splitting at
\((\log X)^2\) and using Mertens below the split.

Convolution with \(1\) gives the unsigned main term
\(f(d)(U-X-2)J_d(1)/d\) and an error
\(O(f(d)Y^\sigma L_X^{O(1)})\). Convolution with \(\mu\) and splitting
at \(\sqrt Y\) gives for the signed sum
\[
 O\!\left(f(d)Y L_X^{O(1)}
     \left(e^{-c'\sqrt{\log Y}}+Y^{-1/(2L_X)}\right)\right).
\]
This uses the Möbius estimate up to the split and the weighted
absolute norm beyond it. The Euler identities
\[
 \frac{J_1(1)}C=\prod_{q>2}\left(1-\frac1{q(q-1)}\right),\qquad
 \frac{f(d)J_d(1)}{dJ_1(1)}
      =\prod_{q\mid d}\frac{q-1}{q(q-1)-1}
\]
give the asserted squarefree density.
These are also the shifted-prime local factors. Conditional on
\(d\mid p+2\) for squarefree odd \(d\), excluding \(q^2\mid p+2\)
retains a factor \(1-1/q\) for \(q\mid d\) and
\(1-1/\varphi(q^2)\) for \(q\nmid d,\ q>2\). Multiplication by
\(1/\varphi(d)\) gives \(g_{\rm sf}(d)\).

Both signed errors give arbitrary powers of \((\log X)^{-1}\)
when \(Y\ge\exp(A_K L_X^2)\) with \(A_K\) large. Inserting \(\chi\)
is a fixed linear combination of divisibility tests at moduli
\([d,e]\) for \(e\mid225\), and increasing \(A_K\) accommodates these
moduli. Partial summation with \(1/(C\log(n-2))\) then gives (P1).
Mertens gives \(f(n)\ll L_X\) on \(I_X\), so \(w_X^\pm(n)\le1\)
for large \(X\).

For the global signed sums the term without the perturbation
already has the signed error just proved. In the perturbation,
\(\lambda^2=1\) and \(\lambda\mu=\mu^2\). The required unsigned
averages of \(\chi\), with or without the squarefree mask, are zero
in the local main term. Indeed, the local mass at \(3,5\) for
both valuations zero is \(3/8\); the mass for both at most one
and at least one positive is
\((1-1/6)(1-1/20)-3/8=5/12\). Thus \(3/8-(9/10)(5/12)=0\).
The valuations at other primes leave this identity unchanged.
Finally \(\chi(p)=1,\ \lambda(p)=-1\) at primes \(p>5\), so the
claimed prime weights follow from the prime number theorem.

The averaged error in (P1) is \(O_K(X(\log X)^{-K})\) even after
multiplying the error at \(d\) by any fixed divisor-function power,
on choosing a larger \(A_K\). The unsquarefree divisor data therefore
give the same small-divisor main term \(\mathfrak S_2X+o(X/L_X)\)
in [equation (B1)](TWIN_CORRELATION_DOUBLING_AND_DILATION.md#the-large-divisor-sum-in-the-twin-prime-estimate)
for both models. With \(\log(n-2)\Lambda(n)\) as the full test,
the two models give \(o(X/L_X)\) and \(2\mathfrak S_2X+o(X/L_X)\),
respectively. Proper prime powers are negligible. Their
large-divisor remainders are therefore
\(-\mathfrak S_2X+o(X/L_X)\) and \(+\mathfrak S_2X+o(X/L_X)\).
The models constrain what can be deduced from divisor totals and
global parity averages. The additional support condition that
\(n-2\) is prime is absent from these test sequences.

## A cutoff variation on power-rough shifted primes

For fixed \(1/4<\eta<1/3\), use the actual prime set
\[
 {\cal R}_\eta(X)=\{X<p\le2X:\ p\ {\rm prime},\
                  \mu^2(p+2)=1,\ P^-(p+2)>X^\eta\}.
\]
Fix any \(A>0\) and take \(D_-=\sqrt X/(\log X)^{A+1}\),
\(D_+=\sqrt X/(\log X)^A\). Define
\[
 J_D(n)=\sum_{\substack{d\mid n\\d\le D}}\mu(d)\log(D/d),\qquad
 V_X(n)=\frac{J_{D_+}(n)-J_{D_-}(n)}{\log(D_+/D_-)}.
\]

**Theorem.** If \(N_2(I)\) counts the twins with lower endpoint
in \(I=(X,2X]\), then
\[
 \sum_{p\in{\cal R}_\eta(X)}
       \bigl(1-\mu(p+2)+V_X(p+2)\bigr)
   =3N_2(I)+O_A(B_X L_X/\log X).                               \tag{P2}
\]
Thus a positive lower bound of order \(B_X/L_X\) for the signed
sum in (P2) along an unbounded sequence suffices. Both cutoffs
in the divisor statistic are within the Bombieri–Vinogradov
range when \(A\) is taken sufficiently large.

Indeed, \(V_X(n)\) averages \(\sum_{d\mid n,d\le D}\mu(d)\) over
\(D_-\le D\le D_+\) with measure \(dD/(D\log(D_+/D_-))\).
A squarefree \(X^\eta\)-rough \(n\in(X+2,2X+2]\) has at most three
prime factors. At a prime, \(V_X(n)=1\) and \(\mu(n)=-1\).
For a product of three primes, each prime is
\(O(X^{1-2\eta})<D_-\) and each pairwise product is
\(>X^{2\eta}>D_+\) for large \(X\). Hence \(V_X(n)=-2\)
and \(\mu(n)=-1\), giving zero total. A semiprime with smaller
factor \(q\le D_-\) has \(V_X(n)=0\) and \(\mu(n)=1\).
The remaining semiprimes have \(0\le V_X(n)\le1\).

For each smaller prime factor
\(D_-<q\le\sqrt{2X+2}\), the dimension-two upper sieve on the
forms \(r\) and \(qr-2\) gives \(O(X/(q(\log X)^2))\) possible
other prime factors \(r\). The Euler factor at the prime \(q\)
is bounded uniformly, and a small fixed-power sieve cutoff
gives a negligible CRT error on intervals of length at least
a constant times \(\sqrt X\). The reciprocal-prime sum over
this range of \(q\) is \(O_A(L_X/\log X)\) by Mertens.
This gives the leakage bound in (P2).

The analytic estimate in (P2) is for the divisor
variation and parity on the same power-rough support.
The signed shifted-prime estimate in Lichtman's Theorem 1.2
averages over all primes in its interval, and the proof's strong
cancellation on integers with occupied subpower prime bands
does not cover this support.

## The large-divisor estimate after removing small factors

In [the large-divisor reduction](TWIN_CORRELATION_DOUBLING_AND_DILATION.md#the-large-divisor-sum-in-the-twin-prime-estimate)
one can remove the short \(m\)-ranges without changing the required
precision. Put \(D=X^{1/3}\) and \(z=X^{c/L_X}\) for a sufficiently
small fixed \(c>0\), and let \(P(z)\) be the product of odd primes at
most \(z\).

**Theorem.** With \(I=(X,2X]\) and \(N_2(I)\) its twin-prime count,
\[
 (\log X)^2 N_2(I)=\mathfrak S_2 X+
 \sum_{\substack{d>D,\ m\ge2\\X<md-2\le2X\\md-2\ {\rm prime}\\
                 (md,P(z))=1}}
       \mu(d)\log m\,\log(md-2)
       +o(X/L_X).                                               \tag{P3}
\]
In the remainder, both \(m\) and \(d\) are \(z\)-rough, with
\(z<m\le(2X+2)/D\) and \(D<d\le(2X+2)/z\).

To prove this, insert \(\mathbf1_{(p+2,P(z))=1}\) before expanding
\(\Lambda(p+2)\). Every prime \(p+2\) survives, and proper prime powers
and the conversion to \((\log X)^2 N_2(I)\) cost \(o(X/L_X)\).
For \(d\le D\) only \((d,P(z))=1\) is needed. Apply upper and
lower Bonferroni truncations to the roughness condition at consecutive
even and odd orders of size \(K L_X\), with \(K\) a sufficiently large
fixed constant. Their moduli \(de\), with \(e\mid P(z)\), are at most
\(X^{5/12}\) when \(c\) is small enough. The prime progression main
term at each \(de\) has density \(1/(\varphi(d)\varphi(e))\).
The product
\[
 V(z)=\prod_{3\le q\le z,\ q\ {\rm prime}}
                         \left(1-\frac1{q-1}\right)
\]
differs from either truncated main factor by
\(O((\log X)^{-M})\) for any prescribed \(M\) when \(K\) is large.
Indeed the degree-\(j\) product tail is at most
\((\sum_{3\le q\le z}1/(q-1))^j/j!\), and the reciprocal sum over
primes is \(L_X+O(\log L_X)\).

The progression errors, summed over \(d,e\), have multiplicity at
most \(\tau(de)\). The divisor-weighted Bombieri–Vinogradov bound
here follows from Cauchy–Schwarz, unweighted Bombieri–Vinogradov
with a sufficiently large logarithmic saving, and the trivial error
\(O(X/\varphi(de))\) for prime counts at these moduli.
Partial summation introduces at most \(O((\log X)^2)\) for the
logarithmic weights. Both this error and the summed Bonferroni
factor error are \(o(X/L_X)\).

It remains to evaluate the two divisor sums in the small-divisor
main term with \((d,P(z))=1\) and a factor \(V(z)\).
Extend \(a(d)=\mathbf1_{d\ {\rm odd}}\mu(d)/\varphi(d)\) by zero
on even \(d\). The restricted sequence is \(a*g_z\), where \(g_z\)
is multiplicative, with \(g_z(q^j)=(q-1)^{-j}\) for odd primes
\(q\le z\) and zero on other prime powers. Thus
\(V(z)\sum_e g_z(e)=1\). Rankin's bound with exponent \(1/\log z\)
gives
\[
 V(z)\sum_{e>\sqrt D}g_z(e)\ll
    \exp\!\left(-\frac{\log D}{2\log z}\right)
    =(\log X)^{-1/(6c)} .
\]
For the bounded product factor in this bound, take the ratio of
the tilted and untilted Euler products. Its logarithm is
\(O(\sum_{q\le z}(q^{1/\log z}-1)/q)=O(1)\) by Mertens.

For \(e\le\sqrt D\) the cutoff \(D/e\) in the \(a\)-sum is at least
\(\sqrt D\), so the estimates in (B1) give limits \(0\) and
\(-\mathfrak S_2\) for \(\sum a\) and \(\sum a\log\), with
arbitrarily large logarithmic savings. For \(e>\sqrt D\) use the
absolute bounds for these sums and the Rankin tail, taking \(c\)
small enough to pay the extra logarithmic factors.
The two normalized restricted sums thus have the same limits to
the required accuracy. The small-divisor main term is
\(\mathfrak S_2X+o(X/L_X)\), proving (P3).

The remaining signed estimate is a bilinear sum with a prime
condition on \(md-2\) and with both factors exceeding the indicated
cutoffs. In particular, sieve information at the prime bands
below \(z\) is already absorbed into the evaluated term.

## Entropy comparison at the rough-pair scale

Let \(H\) be an integer tending to infinity with \(H\le\sqrt{\log X}\), and let
\[
 \mathcal P=\{p\ {\rm prime}:H/2<p\le H\},\qquad
 R=H\sum_{p\in\mathcal P}\frac1p\asymp\frac H{\log H},\qquad
 b=\frac{R}{(\log X)^2}.
\]
Choose \(N\) in \([\lceil XH/3\rceil,\lfloor3XH\rfloor]\) with probability \(1/(ZN)\), where \(Z\asymp1\) normalizes the mass. Let \(Y_p=N\bmod p\) and \(V_p\) denote the uniform law on residues modulo \(p\). For \(m=N+j\), \(1\le j\le H\), let \(m^{(p)}=m/p^{v_p(m)}\) and define
\[
 \begin{split}
 f_{y,p}(m)&=\mu(m)\mathbf1_{m\ {\rm odd}}\mathbf1_{(m,P_y/p)=1},\\
 A_{j,p}&=\mathbf1_{pX<m\le2pX} f_{y,p}(m)f_{y,p}(m+2p)
                    G\big(m^{(p)},(m+2p)^{(p)}\big).
 \end{split}
\]
The array allows \(p\) at either endpoint and has no explicit divisibility indicator for \(p\). Write \(A_p=(A_{j,p})_{j=1}^H\) and set, with logarithms in nats,
\[
 d=\sum_{p\in\mathcal P}
   \left(I(A_p;Y_p)+D_{\rm KL}(\mathcal L(Y_p)\Vert V_p)\right).
\]

**Theorem.** Uniformly for these weights,
\[
 C(G)=\frac{Z}{R}\,\mathbb E\sum_{p\in\mathcal P}\sum_{j=1}^H\frac{A_{j,p}}p
   +O_\eta\left(\frac{\sqrt{db}+d}{R}+\frac1X\right).             \tag{2}
\]
In particular \(d=o(b)\) suffices for an error \(o((\log X)^{-2})\). Each marginal relative entropy in \(d\) is \(O((p/(XH))^2)\).

For the exact change of variables, note that \(f_{y,p}(pn)=-f_y(n)\). If \(p\mid N+j\) and both factors in \(A_{j,p}\) are nonzero, stripping \(p\) gives \(n,n+2\). With \(N=pn-j\), the interval for \(N\) contains all the terms and \(1/(pn-j)=(1+O(X^{-1}))/(pn)\). Consequently
\[
 \mathbb E\sum_{p,j}A_{j,p}\mathbf1_{Y_p=-j\bmod p}
       =(R/Z)C(G)+O(R/X).                                       \tag{3}
\]
An identical identity with absolute values applies to \(U(G)\) for \(G\ge0\).

For each fixed column, take \(v\) uniform modulo \(p\) and put
\(S_p(v)=\sum_j A_{j,p}(\mathbf1_{v=-j\bmod p}-1/p)\).
Its mean is zero, \(|S_p(v)|\le4\), and
\(\mathbb E_v S_p(v)^2\le(2/p)\sum_j |A_{j,p}|^2\). At most two indices \(j\) occupy any residue because \(p>H/2\). Taylor's inequality therefore bounds
\(\log \mathbb E_v e^{sS_p(v)}\le C s^2\sum_j|A_{j,p}|^2/p\)
for \(|s|\le1/4\), with an absolute \(C\).

The unsigned sieve bound gives
\[
 \mathbb E\sum_{p,j}|A_{j,p}|^2/p\ll_\eta b.                     \tag{4}
\]
Indeed, \(m,m+2p\) avoid the two distinct classes \(0,-2p\bmod q\) at each odd \(q\le y\) other than \(p\), and \(N\) has mass \(O(1/(XH))\) at each point. The dimension-two Selberg bound on an interval of length \(O(XH)\) is \(O_\eta(XH/(\log X)^2)\), uniformly in \(p\). One may run the [Selberg argument for fixed forbidden classes](SQUAREFREE_TWIN_CORRELATION_NEAR_ZERO.md) up to a sufficiently small fixed power of \(X\) below \(y\), omitting \(p\). Omitting one prime changes the product by a bounded factor.

Apply the relative-entropy inequality to the conditional law of \(Y_p\) given \(A_p\) against \(V_p\) in both signs. Averaging the conditional logarithmic moment bound and summing over \(p\) gives \(s|\mathbb E\sum_p S_p(Y_p)|\le d+C_\eta s^2b\) for \(0<s\le1/4\). Choosing \(s\) of order \(\min(1,\sqrt{d/b})\), with the case \(d=0\) obtained by taking \(s\) to zero, proves (2) using (3). Summation by parts with \(1/N\) shows that each residue modulo \(p\) has mass \(1/p+O(1/(XH))\). A Taylor bound for relative entropy gives the marginal estimate. The argument uses separate conditional laws for each prime and requires no independence between the true \(Y_p\).

For \(G(a,b)=\mathbf1_{a,b\ {\rm prime}}\), one has \(C(G)=C_2(X)\) and an explicit expression for the reference average. The nonzero entries have either \(p\nmid m\), in which case \(m,m+2p\) are prime, or \(p\mid m\), in which case the stripped pair consists of twin primes. Both cases give \(A_{j,p}=1\). Averaging over \(N\), and replacing \(1/(m-j)\) by \(1/m\) as above, gives
\[
 \begin{split}
 (1-\rho_H)C_2(X)
  &=\frac{\displaystyle\sum_{p\in\mathcal P}p^{-1}
          \sum_{\substack{pX<r\le2pX\\ r,r+2p\ {\rm prime}}}r^{-1}}
         {\displaystyle\sum_{p\in\mathcal P}p^{-1}}
    +O_\eta\left(\frac{\sqrt{db}+d}{R}+\frac1X\right),\\
 \rho_H&=\frac{\sum_{p\in\mathcal P}p^{-2}}{\sum_{p\in\mathcal P}p^{-1}}
           \le\frac2H .
 \end{split}                                                     \tag{5}
\]
Thus one analytic target for the independent-residue average is a lower bound of order \((\log X)^{-2}\) for prime pairs with shifts \(2p\) averaged over the small primes \(p\) in \(\mathcal P\).

## Information revealed by divisibility

Write \(Y=(Y_p)_{p\in\mathcal P}\) and \(V=\prod_p V_p\). Define the array recording actual divisibility by the permitted prime by
\[
 B_{j,p}=\mathbf1_{A_{j,p}\ne0}\mathbf1_{Y_p=-j\bmod p},\qquad
 B_p=(B_{j,p})_{j=1}^H .
\]

**Theorem.** For \(G=1\) and fixed sufficiently small \(\eta>0\), the marked columns satisfy
\[
 \sum_{p\in\mathcal P} I(B_p;Y_p)\gg_\eta b\log H.               \tag{6}
\]
For every fixed \(\eta\) in the original range and any weight for which
\(\mathbb E\sum_{p,j}|A_{j,p}|^2/p\gg b\), the full permissive array satisfies
\[
 I(A;Y)\gg b\,\frac H{\log H}.                                 \tag{7}
\]
The hypothesis of (7) holds for \(G=1\) and sufficiently small fixed \(\eta\). It also holds for the prime-pair projection in (5) if that reference average has the proposed lower bound.

For (6), the dimension-two lower-bound sieve gives \(U(1)\gg_\eta(\log X)^{-2}\). The [lower bound for rough pairs](SQUAREFREE_TWIN_CORRELATION_NEAR_ZERO.md) applies on a fixed-ratio interval with this cutoff. Squares of primes \(q>y\) remove \(O(X/y+\sqrt X)\) pairs. By (3) for absolute values and the upper sieve, \(M=\mathbb E\sum_{p,j}B_{j,p}\asymp_\eta b\). Conditional on a fixed marked column, \(\sum_j B_{j,p}\mathbf1_{v=-j\bmod p}\) for uniform \(v\) is at most two and is nonzero with probability at most \(\sum_j B_{j,p}/p\). Under its true residue this sum is \(\sum_j B_{j,p}\). The relative-entropy inequality for the column therefore implies, after summation,
\[
 \sum_p\big(I(B_p;Y_p)+D_{\rm KL}(\mathcal L(Y_p)\Vert V_p)\big)
     \ \ge\ sM-(e^{2s}-1)(2/H)M .
\]
Take \(s=\tfrac14\log H\) and use the negligible marginal terms.

For (7), let \(K=\sum_{p,j}\mathbf1_{A_{j,p}\ne0}\) and \(L=|\mathcal P|\). Its mean is \(\gg H L/(\log X)^2\) by the hypothesis and \(p\asymp H\). For \(G=1\) and small fixed \(\eta\), the uniform lower bound at the shifts \(2p\) follows by applying the same dimension-two lower sieve on each interval \(pX<m\le2pX\) of length \(pX\). Each odd prime other than \(p\) removes two classes. Both conditions \(p^2\mid m\) and \(p^2\mid m+2p\) imply \(p\mid m\). The upper-bound sieve for the other odd primes on this progression gives \(O_\eta(XH/(p(\log X)^2))\). This is \(o(XH/(\log X)^2)\). Squares at \(q>y\) cost \(O(XH/y+\sqrt{XH})\). For the prime-pair projection the assertion follows directly from the average in (5) before the entropy comparison.

For two distinct entries of \(A\), their endpoint offsets are \(\{j,j+2p\}\) and \(\{i,i+2q\}\). There are \(O(HL^2)\) pairs sharing an endpoint and \(O(H^2L^2)\) pairs with four distinct endpoints. All offsets lie in \([1,3H]\), so any sieve prime greater than \(3H\) gives distinct roots for distinct endpoints. For \(k=3,4\) offsets, dropping primes up to \(3H\) and using the dimension-\(k\) upper sieve bounds the joint probability by \(O_\eta((\log H/\log X)^k)\). For example the same Selberg proof as for (4) works with \(k/q\) and the product is bounded using Mertens' estimate. Hence
\[
 \mathbb E K(K-1)\ll_\eta HL^2\left(\frac{\log H}{\log X}\right)^3+
                    H^2L^2\left(\frac{\log H}{\log X}\right)^4
             =o\left(\frac{HL}{(\log X)^2}\right).
\]
Here \(L\asymp H/\log H\) and \(H\le\sqrt{\log X}\). The first two inclusion–exclusion terms give \(\Pr(K>0)\gg HL/(\log X)^2\).

Given any array with \(K>0\), choose its first nonzero \((j,p)\) in a fixed ordering. Its true conditional law for \(Y\) excludes two residues \(-j,-j-2p\) at every \(q\in\mathcal P\setminus\{p\}\). The product-uniform probability of those simultaneous exclusions is \(\prod_{q\ne p}(1-2/q)\le \exp(-c/\log H)\). A probability law supported on an event of reference mass \(u\) has relative entropy at least \(\log(1/u)\). Average over \(A\) to obtain a lower bound \(\gg HL/((\log X)^2\log H)\) for \(I(A;Y)+D_{\rm KL}(\mathcal L(Y)\Vert V)\). By the Chinese remainder theorem, the marginal of \(Y\) corresponds to a residue modulo \(Q=\prod_{p\in\mathcal P}p=\exp(O(H))\). Its relative entropy against \(V\) is \(O((Q/(XH))^2)\) by the same residue mass estimate. Since \(b\asymp L/(\log X)^2\), this proves (7).

## Relation to Chowla and weighted dilation

For \(p\le y\), the original rough weight vanishes at \(pn\). Allowing \(p\) restores an exact multiplicative change of variables. Recording its occurrence defeats even the column-wise entropy budget by (6). Leaving \(p\) optional still makes the full-array budget fail under (7), because nonzero entries reveal restrictions at the other primes. Neither lower bound controls the column-wise information \(d\) of the permissive array.

Tao's logarithmically averaged two-point Chowla proof uses entropy decrement with genuinely multiplicative inputs and a dense error scale. Formula (2) gives the smaller information budget needed here. Lichtman's November 2026 preprint reports a binary Möbius bound \(O_\epsilon(X/(\log X)^{1-\epsilon})\) for a fixed nonzero shift. Its proof obtains much stronger savings on integers with factors in several subpower prime bands, and treats the contribution with both endpoints outside that part by an unsigned sieve. Power-rough pairs miss all those bands for every fixed \(\eta>0\). [Tao, The logarithmically averaged Chowla and Elliott conjectures for two-point correlations](https://arxiv.org/abs/1509.05422v4), [Lichtman, Improved bounds for the binary Chowla conjecture, Theorem 1.1 and §2](https://arxiv.org/pdf/2611.28392v1).

The [weighted-dilation audit](TWIN_CORRELATION_DOUBLING_AND_DILATION.md) identifies the actual marked prime factors required by the weighted-dilation theorem. If the marking primes are below \(y\), power-rough pairs are outside those marked sums. The variance estimate (4) uses only an unsigned sieve. The signed independent-residue average in (2) still needs an analytic bound.

**Open question.** Obtain the one-sided parity margin in (1) at scale \((\log X)^{-2}\). The column-wise route requires a bound \(d=o(b)\) for the permissive array and a lower bound of order \((\log X)^{-2}\) for the reference average in (5), or corresponding bounds directly for the signed and unsigned degree-three projection in (1).
