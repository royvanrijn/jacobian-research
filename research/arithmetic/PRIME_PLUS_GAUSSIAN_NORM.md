# A prime plus a Gaussian norm value

Let \(b\) be the indicator of positive sums of two integer squares
and let \(K\) be the Landau–Ramanujan constant. The pointwise target
is the count \(R(N)=\sum_{p<N}b(N-p)\). Here we evaluate the arithmetic
model for every target, prove the asymptotic outside a set of at most \(O_A(X(\log X)^{-A})\)
exceptions per dyadic block, and estimate the target-grid major arcs of the
Jiang–Lau–Lü mean square under ZF. The all-target minor-arc problem
is stated in Section 3.

Jiang–Lau–Lü, Theorem 1.1, evaluates \(\sum_{p\le x}b(p-1)\)
with constant \(K\prod_{q\equiv3\ (4)}(1-1/[q(q-1)])\).
Their variable in Theorem 4.1 has the fixed shift one.
Linnik's classical existence and representation-weighted asymptotic
for a prime and two squares are a separate benchmark for \(R(N)\).
The distinction is explicit in the opening of Assing–Blomer–Li,
who count with \(r_2\). [Jiang–Lau–Lü v2](https://arxiv.org/abs/2611.12467v2)
· [Assing–Blomer–Li, §1](https://arxiv.org/abs/2005.13915)

## 1. The uniform model and an almost-all result

Write \(X<N\le2X\), \(L=\log X\), and use the model \(D_X\) of
(CK.7), with \(y=\exp(L^{2/5})\) and sufficiently deep fixed order.
Its outer parameters stay fixed throughout the dyadic block.
Define
\[
 C(N)=K
   \prod_{\substack{q\equiv3\ (4)\\q\nmid N}}
                 \left(1-\frac1{q(q-1)}\right)
   \prod_{\substack{q\equiv3\ (4)\\q\mid N}}(1+1/q).
 \tag{PN.1}
\]
The product has a positive absolute lower bound. It is comparable
to \(J(N)=\prod_{q\equiv3\ (4),\,q\mid N}(1-1/q)^{-1}\)
with absolute constants.

**Theorem ReflectedPrimeNormAlmostAll.** Uniformly for \(X<N\le2X\),
\[
 \sum_{p<N}D_X(N-p)
  =C(N)\frac{N}{(\log N)^{3/2}}\bigl(1+O(L^{-1/2})\bigr),
 \qquad R(N)\ll C(N)\frac{N}{(\log N)^{3/2}} .
 \tag{PN.2}
\]
For every fixed \(A>0\), all but \(O_A(XL^{-A})\) integers in that
block satisfy
\[
 R(N)=C(N)\frac{N}{(\log N)^{3/2}}
                                  \bigl(1+O_A(L^{-1/2})\bigr).
 \tag{PN.3}
\]
This uses classical prime distribution, sieve estimates, and the
one-variable model estimates of
[Connors–Keating, Section 2](CONNORS_KEATING_INDICATOR_CORRELATIONS.md#2-the-indicator-density-and-both-principal-branches).
Those estimates use the classical moderate-height zero-free region
and Siegel’s theorem at polylogarithmic conductors.

**The local factors.** For the prime variable \(p\), put normalized
Haar measure on \(\mathbb Z_q^\times\). At an inert prime \(q\mid N\),
the norm variable \(N-p\) is a unit. At an inert \(q\nmid N\), the
odd valuations of \(N-p\) have probability
\(q/(q^2-1)\) on the prime units, since their ambient mass is
\(1/(q+1)\) and lies entirely over prime units. The marginal
even-valuation mass is \(q/(q+1)\). Dividing the allowed mass by
that marginal gives \(1+1/q\) or \(1-1/[q(q-1)]\), respectively.
At two the norm variable is uniform on one parity coset. The odd
part is \(1\bmod4\) on half that coset, for either parity. Thus its
normalized factor is one. Split odd primes have no restriction.
All local masses are positive. The comparison to \(J(N)\) follows
from convergence of \(\prod_{q\equiv3\ (4)}(1-q^{-2})\) and of the
unramified factors in (PN.1). In the fixed predecessor problem,
\(N-p\) is replaced by \(p-1\) and the normalized odd factors are
all \(1-1/[q(q-1)]\), agreeing with Jiang–Lau–Lü.

**Evaluating the model.** Put
\[
 U_y(n)=\nu_2(n)\prod_{\substack{q\le y\\q\equiv3\ (4)}}
                      \mathbf1_{2\mid v_q(n)},\qquad
 \lambda_y=\frac{2K}{\sqrt L}
                      \prod_{\substack{q\le y\\q\equiv3\ (4)}}(1+1/q).
\]
On \(n\ge X/L^B\), outside small-prime parts exceeding
\(\exp(\sqrt L)\), (CK.14) gives
\(D_X(n)=\lambda_y U_y(n)(1+O(L^{-1/2}))\).
Indeed \(\log n=L+O_B(\log L)\) and the logarithm of the rough
cofactor differs by \(O(\sqrt L)\). The omitted integers up to
\(2X\) cost \(O(XL^{-B}+XL^{O(1)}\exp(-cL^{1/10}))\), for any
fixed \(B\), by the count of the small endpoints and the Rankin
bound on smooth parts.

We evaluate the sum of \(U_y(N-p)\) using the fundamental lemma
and Bombieri–Vinogradov. [Bombieri, On the large sieve](https://doi.org/10.1112/S0025579300005313)

Discard \(p\le y\). Retain inert square
prefixes \(d^2\) with \(d\le L^B\) and dyadic prefixes \(2^e\le L^B\).
The integer error is \(O(X/L^B+\sqrt X)\) by summing
\(X/d^2+1\) over discarded \(d\) and counting multiples of the
first discarded dyadic power. The condition
\(N-p\equiv 2^e\bmod 2^{e+2}\) and \(d^2\mid N-p\) gives one
prime residue modulo \(Q=2^{e+2}d^2\). Discard the class if it
is nonreduced, since \(p>y\). On the cell, exclude every inert
prime at most \(y\) from \((N-p)/d^2\). For a squarefree sieve
divisor \(k\), the excluded-root conditions give one residue
modulo \(Qk\). Relative to the mass \(1/\varphi(Q)\) of the
base cell, its prime local densities are multiplicative, with
factor \(1/q\) for \(q\mid d\), \(1/(q-1)\) for \(q\nmid dN\),
and zero for an inert \(q\mid N\). These follow either from
the totient ratio or from the absence of any reduced prime
class when \(p\equiv N\equiv0\bmod q\).

The dimension is at most \(1/2\) by the inert-prime product
with local densities \(1/(q-1)\). Use the upper and lower
fundamental-lemma weights at level \(D_0=X^{1/100}\), with
\(\log D_0/\log y\asymp L^{3/5}\).
The relative sieve error is smaller than every fixed logarithmic
power. For each base cell \(Q\), the moduli \(Qk\) are distinct
and at most \(4L^{3B}D_0\). The prime disagreement summed over
these classes is \(O_A(XL^{-A})\) by Bombieri–Vinogradov,
after choosing its accuracy to pay for \(O(L^B(1+\log L))\)
base cells. The \(y\) small primes reintroduced by using the
usual progression counts cost at most \(y\) per term, also
negligible at this level.

The main masses are disjoint Haar probabilities in the product
of the prime unit groups. Their omitted Haar square tail is
bounded by \(\sum_{d>L^B}1/\varphi(d^2)\ll L^{-3B/4}\), using
\(\varphi(d^2)^{-1}\ll d^{-7/4}\). Their dyadic tail is
\(O(L^{-B})\). Consequently the count has main term
\(\operatorname{Li}(N)\) times the product of the local masses
through \(y\), with error smaller than any fixed inverse log
power times \(X\), by choosing \(B\) sufficiently large.
Multiplying by \(\lambda_y\) gives (PN.2) for the model.
The inert primes \(q>y\) dividing \(N\) have reciprocal sum
at most \(\log(2X)/(y\log y)\), so restoring their factors and
the absolutely convergent \(O(q^{-2})\) tails is uniform in \(N\).

**Uniform upper bound.** For \(p>w\), where \(w=X^{1/1600}\),
sieve out the residue \(p\equiv0\bmod q\) for every prime
\(q\le w\). At inert \(q\nmid N\) additionally exclude
\(v_q(N-p)=1\), using \(q-1\) classes modulo \(q^2\).
These tests are disjoint. Their combined allowed fraction is
\((1-1/q)^2\) at these inert primes, and \(1-1/q\) at the
other primes. The product is \(\ll J(N)/L^{3/2}\).
The Selberg denominator argument in the
[quadratic-indicator envelope](QUADRATIC_CHARACTER_INDICATOR_PAIRS.md#4-verifying-the-shared-comparison-hypotheses)
with sieve index \(z=X^{1/100}\) applies here too, since
\(w=z^{1/16}\) and each forbidden fraction is at most \(2/q\).
A least-common-multiple atom \(k\) has at most
\(2^{\omega(k)}k\) CRT roots. The \(O(z^{4+o(1)})\)
error and the discarded \(p\le w\) fit (PN.2).

**Passing to almost all targets.** Put \(\xi_X=b-D_X\). Use a
smooth cutoff in \(n\) that is one on \([2X/L^B,2X]\) and supported
on \([X/L^B,3X]\), with bounded logarithmic derivatives. Let
\(\widehat\xi(\alpha)\) denote the Fourier sum with this cutoff
and let \(P(\alpha)=\sum_{m\le2X}\Lambda(m)e(m\alpha)\).
Their convolution differs from the full weighted discrepancy
by \(O(XL^{1-B})\), uniformly in the target.

For any requested \(H>0\) choose a fixed large \(B_0\) and
take major arcs \(q\le L^{B_0}\), \(|\alpha-a/q|\le L^{B_0}/X\).
On the complement, the Vinogradov–Vaughan bound after Dirichlet
approximation with ceiling \(X/L^{B_0}\) gives
\(\lvert P(\alpha)\rvert\ll_H XL^{-H}\). On the major arcs,
the recordwise centering (CK.7) and Section 2 of the
Connors–Keating note give \(\lvert\widehat\xi(\alpha)\rvert
\ll_D XL^{-D}\) for every prescribed fixed \(D\), at suitable
fixed model depth. In applying it, dyadically partition the
cutoff and write \(n=rv\). The rational phase is periodic
in \(v\) and the remaining phase has polynomial-logarithmic
smoothness on each dyad. Drop \(r>\exp(\sqrt L)\) by Rankin
and sum \(1/r\) over the retained records, costing a fixed
power of \(L\). The major-arc measure is \(O(L^{3B_0}/X)\),
so its contribution to any convolution coefficient is \(O(XL^{-2})\)
on choosing \(D\) after \(B_0\). On the minor arcs, Parseval
bounds the sum of squares of all convolution coefficients by
\[
 \sup_{\alpha\ {\rm minor}}|P(\alpha)|^2
                \int_0^1|\widehat\xi(\alpha)|^2\,d\alpha
 \ll_H X^3L^{-2H}.
\]
The last integral is \(O(X)\) by boundedness. The prime
bound used here is Vaughan's rational-approximation estimate.
[Vaughan, On the distribution of \(\alpha p\) modulo 1](https://doi.org/10.1112/S0025579300009025)
 Hence only
\(O_H(XL^{2-2H})\) targets have a weighted minor-arc error
exceeding \(X/L\). On the remaining targets, remove prime
powers at cost \(O(\sqrt X L^2)\), omit primes smaller than
\(X/L^B\) at cost \(O(XL^{1-B})\), and use
\(\log p=L+O_B(\log L)\) on the remaining positive prime
sums for \(b\) and \(D_X\). The weighted error \(O(X/L)\)
therefore gives (PN.3) from the evaluated model, since
\(C(N)\gg1\). This proves the theorem.

## 2. Major arcs for the reflected mean square

Use \(T\asymp X\) for the target, \(M,R\ge X^\delta\) with
\(MR\asymp X\), and \(Z=\exp(L^{0.24})\).
Use the actual progression banks of Jiang–Lau–Lü, Section 2,
with all primes dividing \(T\) deleted. The reciprocal mass
lost to the deletion is \(O(L\exp(-L^{0.1}))\). Let
\(\mathcal A\) consist of products of one prime from each of the \(k_0\) banks
and \(\mathcal V\) be the product of the fixed number of bank
normalizers. Write \(\eta(a/Y)\) for the smooth auxiliary-product
cutoff on \([Y,4Y]\), with \(Y\le\exp(O_{k_0}(L^{0.2}))\) for a nonempty cutoff.
Let \(\alpha_m=\mu(m)m^{iv}\) on the \(Z\)-rough integers in
a subinterval of \([M,2M]\) coprime to \(T\), and zero elsewhere,
with \(|v|\le L^C\). Let \(|\beta_n|\le L^C\) be supported
on the \(Z\)-rough integers in \([R,2R]\) coprime to \(T\). In an extracted Type II sum this last
restriction costs an arbitrarily strong log saving. A prime
\(\ell\mid T\) that divides a rough variable exceeds \(Z\);
counting its multiples costs in total
\(O(XL/Z+L(M+R))\). The original auxiliary weight is bounded
by a fixed power of \(L\), since each integer of size \(O(X)\)
has at most \(O(L)\) prime divisors, and the bank normalizers
are bounded below.

For a fixed smooth cutoff \(\Phi(mn/T)\) supported in
\(0<mn/T<1\), set
\[
 Q_Y(T)=\sum_{k\ge1}
 \left|\frac1{\mathcal V}
  \sum_{\substack{a\in\mathcal A,\,m,n\\T-mn=ak}}
            \eta(a/Y)\alpha_m\beta_n\Phi(mn/T)\right|^2 .
\]
Pairing \(T-mn=ak\) and \(T-rs=ck\) gives
\(cmn-ars=T(c-a)\). For coprime \(a,c\) the converse holds on the
support of the cutoffs, by divisibility and positivity.

Let \(\mathfrak M\) consist of the arcs \(q\le L^B\),
\(|\theta-u/q|\le L^B/Y\) for a fixed \(B>0\). For a smooth
even compactly supported \(\psi\) equal to one on \([-4,4]\)
put
\[
 \begin{split}
 S_j(\theta,\xi)&=
 \sum_{c\in\mathcal A}\sum_{m,n}
   \eta(c/Y)\alpha_m\beta_n\Phi(mn/T)
      e\!\left(\frac{j+\theta+\xi/Y}{T}cmn-\theta c\right),\\
 Q_Y^{\mathfrak M}(T)&=
 \frac1{T\mathcal V^2}\sum_{j\bmod T}
 \int_{\mathfrak M}\int_{\mathbb R}
          \widehat\psi(\xi)|S_j(\theta,\xi)|^2\,d\xi\,d\theta .
 \end{split}                                                    \tag{PN.4}
\]
Here \(\frac1T\sum_j e(j(cm n-ar s)/T)\) imposes divisibility
by \(T\). The remaining phases detect
\((cmn-ars)/T=c-a\) and cut off its size on scale \(Y\).
This is the unrestricted major approximation, prior to the
label restrictions used for the minor-arc transfer. It is
the target-grid expression (LC.22) of the
[earlier determinant audit](PRIMITIVE_NORM_LEADING_CONSTANT.md#7-the-variable-target-issue-in-the-prime-plus-norm-method)
with the prime applications' coefficients.

**Proposition ReflectedPrimeNormMajorArcs.** Under the retained
fixed Dirichlet zero-free input ZF, for the coefficients and
banks just specified and any fixed \(D>0\),
\[
                       |Q_Y^{\mathfrak M}(T)|\ll_D XY L^{-D}.
                                                               \tag{PN.5}
\]
The bound is uniform in \(T\asymp X\), including targets with
arbitrary prime divisors.

First, in the average over \(j\) in (PN.4), additive orthogonality
groups the products \(cmn\) modulo \(T\). Each residue has
\(O(Y)\) integers in the product range. At each integer there
are at most \(\tau_3(cmn)\) contributing tuples with each
coefficient bounded by a fixed power of \(L\).
Cauchy–Schwarz on each residue and the divisor second moment
therefore bound the average by \(O(XY^2 L^{C_1})\),
independently of \(\theta,\xi\). This pays for truncation of
\(\xi\) at a sufficiently large fixed power of \(L\), since
\(\widehat\psi\) is Schwartz and
\(|\mathfrak M|\ll L^{3B}/Y\).

Fix \(\theta=u/q+\gamma\) in one arc. Set \(Q=qT\).
The additive numerators \(qj+u\) are distinct modulo \(Q\).
Replacing their average by the average over all numerators
costs at most \(q\). The remaining product phase is
\[
 e((Y\gamma+\xi)\,cmn/(TY)).
\]
On a fixed compact set in \(cmn/(TY)\), separate that phase
by smooth Mellin inversion, and likewise separate
\(\Phi(mn/T)\). Their Mellin \(L^1\) norms and fixed derivative
seminorms cost powers of \(L\) on the truncated range of
\(\xi\). Restrict the Mellin variables to fixed powers of \(L\)
by Schwartz decay. The same additive residue bound just
used, now for separated coefficients of the same absolute
size, justifies these tails in the mean square.

All three variables \(c,m,n\) are units modulo \(Q\), since
the least prime factors of \(c\) exceed \(L^B\) and those of
\(m,n\) exceed \(Z\). For each Mellin term,
additive orthogonality over all \(Q\) numerators and then
multiplicative orthogonality on \((\mathbb Z/Q\mathbb Z)^\times\)
identify the resulting average with
\[
 \frac1{\varphi(Q)}\sum_{\chi\bmod Q}
       |C_\chi A_\chi B_\chi|^2 .
                                                               \tag{PN.6}
\]
Here \(A_\chi\) and \(B_\chi\) are the \(\alpha_m,\beta_n\)
character polynomials with Mellin twists of polynomial-logarithmic
height, and \(C_\chi\) is the corresponding \(c\) polynomial
with \(\eta(c/Y)e(-\theta c)\).
In particular, there is no sum of separate absolute errors
over the \(T\) additive classes.

For characters whose primitive conductor does not divide
\(\operatorname{lcm}(q,4)\), one has
\(|C_\chi|\ll_E Y L^{-E}\) for every fixed \(E>0\).
To see this, expand \(e(-uc/q)\) on the prime units using
characters modulo \(q\), at a polynomial-logarithmic cost.
Write \(c=p c'\) with \(p\) from the first bank. With \(c'\)
fixed, the bank restriction \(p\equiv1\bmod4\) uses only
characters modulo four. No resulting character on \(p\)
is principal in this case. The sharp prime estimate under ZF
in the [pointwise proof](../../results/pointwise-binary-splitting-semigroups/prime-estimates.tex)
permits modulus \(\operatorname{lcm}(Q,4)\), and all \(p\)
are larger than \(\exp(cL^{0.1})\). Partial summation with
the cutoff and \(e(-\gamma p c')\) costs fixed powers of
\(L\) on \(p c'\asymp Y\). The deleted target divisors cost
\(O(L)\) per prime sum. The resulting error is smaller than
every log power times \(Y/c'\); summing \(1/c'\) over the
other banks costs \(O_{k_0}(1)\).

For these characters, bound the mean of \(|A_\chi B_\chi|^2\)
by \(O(XL^{C_2})\). The product variable \(mn\) ranges
over \(O(X)\) integers, with at most \(O(1)\) in a class
modulo \(Q\), and the coefficient at each integer is bounded
by \(L^C\tau(mn)\). This proves the claimed bound by
character orthogonality and the ordinary divisor second moment.
There are only \(O(q)\) exceptional characters, all induced
from the fixed polylogarithmic modulus \(\operatorname{lcm}(q,4)\).
For them use Jiang–Lau–Lü, Lemma 7.2, and partial summation
on the subinterval defining \(A_\chi\). This gives
\(|A_\chi|\ll_E ML^{-E}\). The frequencies here are
polylogarithmic, so the same Perron shift with ZF suffices:
on \(\Re s=1-c/\log Z\), the fixed-strip logarithmic derivative
bound and the Euler bound at \(1+c/\log Z\) bound the
reciprocal \(L\)-function by a power of \(L\), with the
principal pole regularized near one. The finitely removed
Euler factors cost powers of \(L\), and the shift saves
\(\exp(-c' L^{0.76})\). The restriction \((m,T)=1\)
only subtracts rough multiples of target primes \(>Z\),
of total size \(O(ML/Z+L)\). Use
\(|B_\chi|\ll RL^C,\ |C_\chi|\ll Y\) and
\(\varphi(Q)\gg Q/L\) to bound the exceptional contribution
to (PN.6) by \(XY^2\) times any fixed inverse log power.
This also bounds the nonexceptional contribution by the
estimate just proved.

Absorb the \(q\) and Mellin costs in these arbitrary log
savings, then integrate over \(\theta\). The factor \(1/Y\)
in the measure bound removes one \(Y\) from \(XY^2\).
This proves (PN.5).

### Removing the auxiliary weight uniformly

Take \(J\) disjoint families of \(k_0\) banks from Jiang–Lau–Lü,
Section 2, trimmed by \(q\nmid T\), with \(1\le J\ll_{k_0}\log L\)
as in their grid. In \(W_j\), normalize the number of bank
prime divisors of the norm by \(\sum_q1/(q-1)\) in that bank.
Then for \(\overline W=J^{-1}\sum_jW_j\), uniformly in \(T\asymp X\),
\[
 \left|R(T)-\sum_{p<T}b(T-p)\overline W(T-p)\right|
            \ll_{k_0} \frac{C(T)T}{L^{3/2}\sqrt J}.
 \tag{PN.7}
\]
Under the independent bank law conditioned on prime units,
the events \(q\mid T-p\) have probability \(1/(q-1)\).
Their band second moments are bounded, so
\(\mathbb E(\overline W-1)^2\ll_{k_0} J^{-1}\).
Expand \(Z_W(n)=(\overline W(n)-1)^2=\sum_Q c_Q1_{Q\mid n}\).
As in Jiang–Lau–Lü, (8.4)–(8.5), \(Q\) is squarefree
on at most \(2k_0\) bank primes,
\(\sum_Q|c_Q|\le\exp(O_{k_0}(L^{0.2}))\),
\(\sum_Q|c_Q|/\varphi(Q)\ll_{k_0}1\), and
\(\sum_{Q:q\mid Q}|c_Q|/\varphi(Q)\ll_{k_0}1/q\) for each bank
prime \(q\). For the last two bounds, use the bounded band
second moments on the constant, one-weight and two-weight terms
before collecting coefficients. Prescribing \(q\) in one band
costs \(O(1/q)\), including repeated selection of that prime.

Apply the upper sieve of (PN.2) to \(Z_W(T-p)\), with \(p\)
temporarily an integer in \((w,T)\). For a sieve atom \(r\),
put \(g(r)=\prod_{q\mid r}\rho_q\), its forbidden fraction.
Since the bank primes are split and avoid \(T\), \(Q\mid T-p\)
is incompatible with the atom if \((Q,r)>1\). Otherwise CRT
gives main density \(g(r)/Q\) and error \(O(r^2g(r))\).
Replacing \(\mathbf1_{(Q,r)=1}/Q\) by \(1/\varphi(Q)\) in
the signed sum costs at most
\[
 \ll_{k_0} \sum_{q\ {\rm bank}}q^{-2}
              +\sum_{\substack{q\mid r\\q\ {\rm bank}}}q^{-1}
 \ll_{k_0} L\exp(-cL^{0.1}),
\]
using the two reciprocal-totient bounds and \(r\le z^2\).
The cost in the quadratic sieve is negligible, since
\(\sum g([d_1,d_2])\le\prod_{q\le w}(1+3\rho_q)\ll L^6\).
The CRT errors give \(z^{4+o(1)}\exp(O_{k_0}(L^{0.2}))\).
The primes \(p\le w\) cost \(O(wL^{2k_0})\).
The main term is the same Selberg quadratic form as in (PN.2),
multiplied by \(\mathbb E Z_W\). This bounds the weighted
variance by \(O_{k_0}(C(T)T/(J L^{3/2}))\). Cauchy–Schwarz and
(PN.2) give (PN.7).

## 3. Primitivity and the open minor-arc transfer

On the unit-supported comparison, write \(t=cmn-ars\).
In the determinant kernel the only terms have \(T\mid t\)
and \(|t/T|\ll Y\). When the auxiliary products \(a,c\)
are coprime, \(t=0\) would force \(a\mid mn\), contrary to
roughness. The possible endpoint gcds \((m,s)\), \((r,n)\)
are coprime to \(T\), divide \(t/T\), and are \(Z\)-rough.
They therefore equal one. The conclusion holds for these
kernel terms as well as the original paired equations.

Overlapping auxiliary products have negligible cost. For fixed
\(a,c\asymp Y\), the sum over paired \(k\) costs at most
\(XL^{O(1)}/Y\). For example, on the reduced progression
\(T-ak\), use \(\tau(n)^2\le\tau_4(n)\). In every ordered
four-factorization, the product of the three factors other than
a largest one is \(O(X^{3/4})\). Counting the fourth factor
in one class modulo \(a\) bounds the second moment by
\((X/a)L^{O(1)}+X^{3/4}L^{O(1)}\). Cauchy–Schwarz then
gives the two-progression bound. For the major kernel,
\(T\mid c(mn)-a(rs)\) leaves \(O(1)\) choices of \(rs\)
in a fixed residue class for each \(mn\), and conversely,
since \((ac,T)=1\). Its absolute cost is
\(XL^{O(1)}|\mathfrak M|\ll XL^{O(1)}/Y\) per \((a,c)\).
There are \(O(Y^2\sum_{q\ {\rm bank}}q^{-2})\) pairs sharing
a bank prime. Repeated shared labels are rare by the reciprocal
bank law, with repetition probability \(O_{k_0}(J_0^2/q_{\min})\)
when \(J_0\) labels are inserted per bank. The same per-pair
bound pays this cost. Both errors are smaller than every log
power times \(XY\).

With the shared auxiliary product \(D_{\rm sh}\) coprime to
\(T\), use the physical points
\(P=(D_{\rm sh}cm,s)\), \(Q'=(D_{\rm sh}ar,n)\).
When the auxiliary labels are distinct, the preceding gcd
calculation makes both primitive in the ambient lattice.
Their coordinates are units modulo \(T\). They belong to
the same lattice
\[
 \Lambda_\lambda=\{(u,v)\in\mathbb Z^2:v\equiv \lambda u\bmod T\},
       \qquad \lambda=P_2P_1^{-1}\pmod T .
\]
The coordinate map \((u,v)\mapsto(u,(v-\lambda u)/T)\)
preserves primitivity here, because \((u,T)=1\), and divides
the determinant by \(T\). It also leaves the first coordinate
divisibility by bank primes intact.

The remaining mean-square target is
\[
 |Q_Y(T)-Q_Y^{\mathfrak M}(T)|\ll_D XY L^{-D}.
                                                               \tag{PN.8}
\]
The transported physical region has area of order
\(D_{\rm sh}Y=X^{o(1)}\) in each lattice, since its original
area is of order \(D_{\rm sh}YMR\) and the index is \(T\).
Jiang–Lau–Lü, Lemmas 6.3–6.4 and Proposition 6.5, use
polynomial-size physical sides to average the congruence
conditions along paths. Rescaling the determinant leaves
that averaging step to be proved across the projective
lattices. The proof of (PN.5) estimates the unprojected
major expression and does not supply (PN.8).
A full adaptation also needs to complete the reflected
Type I reduction. The fixed-shift argument in their Section 3
uses the reduced class \(-1\bmod d\). The class \(T\bmod d\)
can be nonreduced when \((d,T)>1\).
The model evaluation in Section 1 determines the candidate
uniform main term independently of that reduction.

## 4. The coefficient detector and the progression benchmark

For \(f_z(n)=z^{\Omega(n)}\), the bank identity is exact:
\(\xi_z(D_{\rm sh}n)=z^{\Omega(D_{\rm sh})}\xi_z(n)\) for
any candidate model with the same eigenvalue as \(f_z\).
For \(0<|z|\le1\) a possible normalization on the exact bank
count is to divide the endpoint by \(z^{Mb}\), where there
are \(b\) active banks and \(M\) exact marks per bank. This
is bounded for \(f_z\), and for a model bounded in absolute
value by a constant times \(|z|^{\Omega_{\mathcal B}(n)}\).
After removing the shared labels the remaining factor is
\(z^{-\sum_g t_g}\), with \(\Omega(D_{\rm sh})=Mb-\sum_g t_g\).
For fixed probe complexity \(\sum_g t_g=O(\log\log X)\),
independent of the eventual mark length. This gives only a
fixed power of \(\log X\) on compact annuli bounded away
from zero. The normalized model bounds, rough means and
sparse envelope for these weights would still need to be
proved before invoking a comparison theorem.

For the unweighted finite generating sum over \(a+b=N\)
with \(a,b\ge2\), the coefficient of \(zw\) counts ordered
prime pairs. A contribution that imposes at least two
bank-prime divisors on one endpoint has zero \(zw\) coefficient.
If an analytic remainder is \(E_N(z,w)\), the coefficient
error is bounded by
\[
 |[zw]E_N|\le
 (rs)^{-1}\sup_{|z|=r,\ |w|=s}|E_N(z,w)|.
\]
Recovering a Goldbach-scale asymptotic by this inequality
requires its right side to be \(o(N/(\log N)^2)\) on a
contour where the main term's coefficient is also controlled.

For twins, Lü–Wang–Wu, Theorem 1, gives
\(\sum_{p\le X}\mu(p+2)=o(\pi(X))\), equivalent by partial
summation and removal of proper prime powers to \(M_2(X)=o(X)\).
Use this as the global input. In Murty–Vatwani, Theorem 1.1(a)
and the remark immediately after the theorem, the restricted
residue \(n\equiv-2\bmod q\) suffices. Thus with
\(0<\varepsilon<1/2\) a sufficient separate estimate is, for
every \(A>0\),
\[
 \sum_{\substack{q\le X^{1/2+\varepsilon}\\(q,2)=1}}
 \max_{y\le X}
 \left|\sum_{\substack{n\le y\\n\equiv-2\bmod q}}
                 \Lambda(n)\mu(n+2)-\frac{M_2(y)}{\varphi(q)}\right|
       \ll_A X(\log X)^{-A}.
\]
In their theorem choose \(\theta=1/2-\varepsilon/2\). The
ordinary prime distribution at \(X^\theta(\log X)^C\) is
then classical Bombieri–Vinogradov, and the hypothesis
above covers their shifted Möbius range \(X^{1-\theta}\).
This would give the full weighted twin-prime asymptotic.
[Lü–Wang–Wu, Theorem 1](https://arxiv.org/abs/2610.20242)
· [Murty–Vatwani, Theorem 1.1(a)](https://mrammurty-ai.github.io/TwinPrimes-Parity.pdf)
