# An order-uniform obstruction to exact-power repair

Multiplying a primitive \(abc\) triple by a suitable common integer
can make its three letters exact powers of prescribed coprime
orders. On triples whose valuation masses nearly agree at the
individual orders, the height of the least such repair has a
lower bound uniform in the three chosen orders. The unrestricted
count by this repaired height has a matching obstruction at the
scale needed for a strict one-third saving.

## The least exact-power multiplier

Let \(U_1=a,\ U_2=b,\ U_3=c\), where \(a+b=c\) and
\((a,b)=1\). Fix integers \(q_i\ge2\) and let
\(Q=\operatorname{lcm}(q_1,q_2,q_3)\).
A common multiplier \(M\ge1\) repairs the triple if each
\(M U_i\) is a \(q_i\)-th power. For a prime of exponent \(e\)
on letter \(i\), its exponent in \(M\) must satisfy
\[
        r\equiv-e\pmod{q_i},\qquad
        r\equiv0\pmod{q_j}\quad(j\ne i).
 \tag{EP1}
\]
Primitivity ensures that a prime occurs on at most one letter.
The congruences at that prime are compatible exactly when
\(\gcd(q_i,\operatorname{lcm}_{j\ne i}q_j)\mid e\).
Equivalently, a repair exists exactly when each letter \(U_i\)
is an exact \(g_i\)-th power, where
\(g_i=\gcd(q_i,\operatorname{lcm}_{j\ne i}q_j)\).
Here \(g_i=1\) imposes no condition.
Whenever they are compatible throughout the triple, write
\(r_i(e)\in[0,Q)\) for the least residue in (EP1).
At primes outside \(abc\) the least repair uses exponent zero.
This defines a unique least multiplier \(M_{\boldsymbol q}\).
Every other repair is \(M_{\boldsymbol q}t^Q\) for an integer
\(t\ge1\).
If \(x_i^{q_i}=M_{\boldsymbol q}U_i\), then
\(\operatorname{rad}(\gcd(x_1,x_2,x_3))=\operatorname{rad}(M_{\boldsymbol q})\),
since \(\gcd(M_{\boldsymbol q}U_1,M_{\boldsymbol q}U_2,M_{\boldsymbol q}U_3)
=M_{\boldsymbol q}\).
Thus the roots are primitive in the usual generalized
Fermat sense exactly when \(M_{\boldsymbol q}=1\).
Set
\[
                    H_{\boldsymbol q}=c M_{\boldsymbol q}.
\]
The repair exists for every primitive triple when the orders
are pairwise coprime. Conversely, if a repair exists and
the positive valuation orders on each of two letters have
gcd one, the three \(q_i\) must be pairwise coprime. This
applies when each of the two letters has at least two
different prime valuation orders. A common divisor of any
pair of orders divides every valuation on either affected
letter by (EP1).
The effective-order reduction and the construction for
pairwise coprime orders are classical.
[Darmon–Granville, §1.1, p. 519](https://www.math.mcgill.ca/darmon/pub/Articles/Research/12.Granville/pub12.pdf)
· [Bennett–Mihăilescu–Siksek, §4, p. 22](https://personal.math.ubc.ca/~bennett/BeMuSi-Springer-2016.pdf)

## The height lower bound at matching valuation masses

For each integer \(e\ge1\), let
\[
 u_{i,e}=\frac{\log\!\left(\prod_{p:\,v_p(U_i)=e}p\right)}{\log c},
 \qquad m_e=\min_{1\le i\le3}u_{i,e},\qquad
 L=\frac{\log\operatorname{rad}(abc)}{\log c}.
\]
The normalization uses the actual height \(c\).
Suppose the repair exists. For every \(e\) with \(m_e>0\),
all three residues in (EP1) are defined and
\(\sum_i r_i(e)\equiv-e\pmod Q\). Nonnegativity gives
\(\sum_i r_i(e)\ge Q\lceil e/Q\rceil-e\). It follows that
\[
\begin{aligned}
 \frac{\log_c H_{\boldsymbol q}}{Q}
  &= \frac1Q+\sum_{i,e}\frac{r_i(e)}Q u_{i,e}\\
  &\ge \sum_e m_e\left\lceil\frac eQ\right\rceil
          +\frac{1-\sum_e e m_e}{Q}\\
  &\ge \sum_e m_e .
\end{aligned}
 \tag{EP2}
\]
Here terms with zero mass can be omitted, and
\(\sum_e e m_e\le\sum_e e u_{i,e}=\log_c U_i\le1\)
for every \(i\). No hypothesis on the sizes of the \(q_i\)
or the exponents \(e\) is needed.

Write \(u_i=(u_{i,e})_{e\ge1}\).
At each \(e\), the mean of the three masses exceeds their
minimum by at most one third of the sum of their pairwise
absolute differences. Thus
\[
 \frac{\log_c H_{\boldsymbol q}}{Q}
       \ge \frac L3-\frac13\sum_{i<j}\|u_i-u_j\|_1.
 \tag{EP3}
\]
In particular, if \(L\ge(1-\eta)\lambda\) and every
\(\|u_i-u_j\|_1\le\eta\lambda\), then the weighted repair
height exponent is at least \((1/3-4\eta/3)\lambda\).
The same local formula gives
\[
 \log_c\operatorname{rad}(M_{\boldsymbol q})
     =L-\sum_i\sum_{e:\,q_i\mid e}u_{i,e},
\]
because \(r_i(e)=0\) exactly when \(q_i\mid e\).
For example, if all valuation orders are prime and every
individual \(u_{i,e}\) is at most \(\eta\lambda\), this
repair radical has exponent at least \(L-3\eta\lambda\)
uniformly in the orders. Each \(q_i\ge2\) divides at most
one prime valuation order.

The height comparison concerns the masses at the original integer
valuation positions. Weak convergence to a common measure
of scaled orders alone does not imply the required
\(\ell^1\) proximity. For example, identical weak limits can
come from disjoint dense sets of positions.

For a formal contrast, choose distinct primes \(q_i>5\) of size
comparable to \(L^{-1/2}\), with \(L\to0\), so
\(Q=q_1q_2q_3\) has size comparable to \(L^{-3/2}\).
On letter \(i\), put mass \(L/3\) on the degrees \(e\) with
\[
          2\le Le\le4,\qquad
          e\equiv -Q/q_i\pmod{q_i},\qquad e\equiv i\pmod5.
\]
The three sets of positions are disjoint. On each letter the
mesh of the scaled positions is \(5Lq_i\to0\), so equal
weights give the same uniform weak limit after normalization.
The weighted mean of \(Le\) is \(3+o(1)\). An \(o(1)\)
redistribution of the normalized weights toward an endpoint
can make this mean exactly \(3\), corresponding to numerical
height one on every letter. Every occupied position on
letter \(i\) has \(r_i(e)/Q=1/q_i\) by (EP1). The resulting
formal weighted repair cost is
\[
                  \frac1Q+\frac L3\sum_i\frac1{q_i}
                        =O(L^{3/2})=o(L).
\]
This is a valuation-profile example without an assertion of
realization by primitive additive triples.

## A lower bound for the unrestricted repaired-height count

For pairwise coprime \(q_i\), define
\[
 B_{\boldsymbol q}(Y)
   =\#\{(a,b,c):a+b=c,\ (a,b)=1,\ a,b>0,\
                  H_{\boldsymbol q}(a,b,c)\le Y\}.
\]
This counts primitive quotients once, using only their least
repairs. Put \(s_i=r_i(1)\). Each \(s_i\) is positive and
divisible by \(Q/q_i\), so \(s_i\le Q-Q/q_i\).
Also \(\sum_i s_i\equiv-1\pmod Q\), whence
\[
                   \sum_i s_i=dQ-1,\qquad d\in\{1,2\}.
\]
For a primitive triple with all three letters squarefree,
\(M_{\boldsymbol q}=a^{s_1}b^{s_2}c^{s_3}\) and therefore
\(H_{\boldsymbol q}\le c^{dQ}\).
The prime exponent labels in \(H_{\boldsymbol q}\) are
\(s_1,s_2,s_3+1\) according to the original letter. They
are all positive and pairwise distinct. Modulo \(q_1\),
\(s_2\) differs from both \(s_1\) and \(s_3+1\). Modulo
\(q_2\), \(s_1\) differs from \(s_3+1\).
Thus \(H_{\boldsymbol q}\) itself recovers the ordered
squarefree triple by its prime factorization.

There are \(\gg Z^2\) such primitive squarefree triples with
\(c\le Z\). In the box \(1\le a,b\le T\), every local condition
modulo \(p^2\) excludes a proportion \(O(p^{-2})\).
The conditions require \(p^2\nmid a,b,a+b\) and
\(p\nmid(a,b)\). The class \((1,1)\bmod p^2\) is
admissible for every prime. Consequently the product of the
local densities is positive. For a fixed cutoff \(P\),
the Chinese remainder theorem gives the product density
over \(p\le P\) in the limit \(T\to\infty\).
The number of pairs excluded by some \(p>P\) is
\[
            O\left(T^2\sum_{p>P}p^{-2}+T^{3/2}\right).
\]
For each of \(p\mid(a,b)\), \(p^2\mid a\), and \(p^2\mid b\),
use \(O(T^2/p^2)\). For \(p^2\mid a+b\), there are
\(O(T^2/p^2+T)\) pairs per prime, and \(p\le\sqrt{2T}\).
Taking \(P\) large and fixed gives the lower density.
Use \(T=\lfloor Z/2\rfloor\).

Taking \(Z=Y^{1/(dQ)}\) now gives
\[
                  B_{\boldsymbol q}(Y)\gg Y^{2/(dQ)}
                                             \ge Y^{1/Q}.
 \tag{EP4}
\]
The lower constant and the threshold for \(Y^{1/(dQ)}\)
can be absolute. The lower bound counts distinct repaired
heights as well as distinct primitive quotients. Each
quotient uses its least repair.

## The exact counting exponent on the squarefree benchmark

Put \(\alpha=s_1,\ \beta=s_2,\ \gamma=s_3+1\). For the
squarefree restriction of \(B_{\boldsymbol q}(Y)\), the exact
exponent is
\[
   E_{\boldsymbol q}^{\rm sf}
      =\max\left\{\frac2{\alpha+\beta+\gamma},
                   \frac1{\alpha+\gamma},\frac1{\beta+\gamma}\right\}.
                                                               \tag{EP6}
\]
More precisely the count lies between
\(\gg_{\boldsymbol q}Y^{E_{\boldsymbol q}^{\rm sf}}\) and
\(O_{\boldsymbol q}(Y^{E_{\boldsymbol q}^{\rm sf}}(1+\log Y)^2)\).
The same holds for the count of distinct heights.

On \(a\le b\), the exact height
\(a^\alpha b^\beta(a+b)^\gamma\) is comparable, with a
fixed constant, to \(a^\alpha b^{\beta+\gamma}\). In a
dyadic box with lower endpoints \(A\le B\) the number of
tuples is at most \(O(AB)\). The conditions \(0\le\log A\le
\log B\) and \(\alpha\log A+(\beta+\gamma)\log B\le\log Y\)
give the maximum of \(\log(AB)/\log Y\) at the balanced or
one-sided vertex of this triangle. These exponents are
\(2/(\alpha+\beta+\gamma)\) and \(1/(\beta+\gamma)\).
Interchanging \(a,b\) gives the other one-sided exponent.
There are \(O_{\boldsymbol q}((1+\log Y)^2)\) boxes.

The lower exponent \(2/(\alpha+\beta+\gamma)\) comes from
the two-variable sieve above. For the other two lower
bounds fix one summand to one. Consecutive positive
integers are simultaneously squarefree with positive
density, by the finite CRT sieve and the tail bound
\(O(T\sum_{p>P}p^{-2}+\sqrt T)\) among integers at most \(T\).
These pairs are automatically primitive. For example
\(b=1\) gives height at most \(2^\gamma a^{\alpha+\gamma}\).
The injectivity of the height labels proved above applies
throughout. If \(d=2\) then \(\alpha+\beta+\gamma=2Q\),
and all three individual weights are less than \(Q\).
Thus in this case (EP6) gives exactly \(1/Q\).

## How low the weighted repair height can be made

The minimum scale in (EP3) has a complementary upper
construction. Given any fixed \(\lambda,\epsilon>0\), choose
large coprime \(q_1,q_2\). Choose \(q_3\ge2\) in the simultaneous
classes
\[
 q_3\equiv-q_2^{-1}\pmod{q_1},\qquad
 q_3\equiv-q_1^{-1}\pmod{q_2}.
\]
Then all three orders are pairwise coprime and
\(Q/q_i\equiv-1\pmod{q_i}\) for \(i=1,2\). Hence
\[
 \frac{r_i(e)}Q=\frac{[e]_{q_i}}{q_i}\le \frac e{q_i}
        \quad(i=1,2),\qquad \frac{r_3(e)}Q<1,
\]
where \([e]_q\) is the least nonnegative residue.
If \(\rho_i=\log_c\operatorname{rad}(U_i)\), then
\[
             \frac{\log_c H_{\boldsymbol q}}Q
                   \le \rho_3+\frac1Q+\frac1{q_1}+\frac1{q_2}.
 \tag{EP5}
\]
This holds with no truncation of the valuation orders.
Choose the orders so that the last three terms sum to at
most \(\epsilon\lambda\), and assign slot \(3\) to a least
radical letter. The choice of the three orders is fixed
at the budget, with only the three possible positions of
the designated letter. At total budget \(L\le\lambda\),
the weighted repair height is then at most
\((1/3+\epsilon)\lambda\). For these orders
\(\sum_i s_i=Q-1\), so their unrestricted lower bound in
(EP4) is in fact \(Y^{2/Q}\).

## Consequence for counting the balanced sector

Consider dyadic height families \(X<c\le2X\) with
\(L=(1-o(1))\lambda\), pairwise \(\ell^1\) differences
\(o(\lambda)\) at the individual valuation orders, and
the prime-order condition above. The losses can be read
in the iterated limit, first taking \(X\to\infty\) at
fixed \(\lambda\), then sending the normalized tolerances
to zero with \(\lambda\). Alternatively (EP3) gives the
corresponding finite-tolerance statement.
For any compatible exact-power orders, even if chosen
from the profile, (EP3) forces any repaired-height cap covering
the family to satisfy
\[
             \frac{\log Y}{Q\log X}
                       \ge (1/3-o(1))\lambda .
\]
By (EP4) the corresponding unrestricted count
\(B_{\boldsymbol q}(Y)\) is itself at least
\(X^{(1/3-o(1))\lambda}\). Thus bounding the family solely
by the unrestricted count at its repaired height cannot
prove a strict one-third saving in this sector. This holds
even when every common multiplier is chosen minimally and
quotients are counted exactly once.

The statements apply to actual triples satisfying their
hypotheses. They do not assert the existence of a sequence
of actual triples saturating the small radical budget with
these valuation profiles. A count retaining the original
three radical caps, the original height, or the distribution
of the repair factors is outside the obstruction.
