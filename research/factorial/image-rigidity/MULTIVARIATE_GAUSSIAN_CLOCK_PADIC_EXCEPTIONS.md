# Prime-adic restrictions on exceptional Gaussian clock parameters

Fix an ordinary binary factorial core and consecutive degrees
\(n\ge2\), \(N=n+1\). Write \(b_{j,k}>0\) for the original
moment coefficients, \(j=n,N\), and put
\[
 H_j(q,z)=\sum_{k=0}^j b_{j,k}q^{K_j k(j-k)}z^k,\qquad
                             K_j=nN/j .
\]
For \(q=e^{-\tau/(2nN)}\) this is the common positive Gaussian
strength \(\tau\) on the exact core, as in the
[core angular eliminant](MULTIVARIATE_ARBITRARY_CORE_DIFFUSE_EXTENSION.md).
For every algebraic \(q\in(0,1)\), the bounds below restrict
each valuation of a possible exceptional parameter.

For a prime \(p\), let \(a_{j,k}=v_p(b_{j,k})\), and define
\[
 \begin{aligned}
 \lambda_{j,k}&=a_{j,0}+\frac{k}{j}(a_{j,j}-a_{j,0}),\\
 P_p&=\max\left(0,\max_{\substack{j=n,N\\1\le k<j}}
             \frac{\lambda_{j,k}-a_{j,k}}{K_j k(j-k)}\right),\\
 d_{j,k}&=a_{j,k}-a_{j,k-1}\qquad(1\le k\le j),\\
 V_p&=\frac12\sum_{j=n,N}
                     \left(\max_k d_{j,k}-\min_k d_{j,k}\right).
 \end{aligned}
\]
**Theorem.** If the two Gaussian-deformed polynomials at
an algebraic \(q\in(0,1)\) have equal root arguments, then for
every embedding at every \(p\) one has
\[
                             -V_p\le v_p(q)\le P_p.             \tag{1}
\]
The valuation is normalized by \(v_p(p)=1\).
At a good prime, where every \(b_{j,k}\) is a unit, an exceptional
clock therefore has \(v_p(q)=0\). This excludes both numerator
and denominator primes. In particular all rational exceptional
\(q\)'s belong to the following explicitly finite list, which
requires no angular resultant to generate:
\[
   q=\prod_{p\le P_*}p^{e_p}\in(0,1),\qquad
       -\lfloor V_p\rfloor\le e_p\le\lfloor P_p\rfloor ,
 \quad P_*=\max\bigl(N,N\max_i\max(a_i,b_i)\bigr).
\]
Here \(a_i,b_i\) in the maximum are the exponent columns
of the core, and the maximum of an empty list is zero.
Primes above \(P_*\) are good because all factorial arguments
and both binomial degrees are smaller than such a prime.
This is a necessary constraint on exceptions, not a sufficiency
test for collision.

The core can instead be any pair with positive rational
coefficients whose even-degree polynomial is nonnegative on the
real axis, using its own good primes. For the ordinary core the
even reference at positive strength is strictly positive on the
real axis. Indeed Gaussian attenuation is the average of
\(B_j(e^{X-\tau/2}z)\) with \(X\) centered real Gaussian of variance
\(\tau/j\). For real \(z\ne0\) the argument misses the finite zero
set of the nonzero even polynomial almost surely, and the constant
coefficient handles zero. Thus equal arguments can occur only
at nonreal roots \(z,w\). They imply the algebraic identity
\(z/\bar z=w/\bar w\) between two off-diagonal root quotients.
It suffices to exclude such a quotient equality at any one place.

## Positive clock valuation

Suppose \(v_p(q)>P_p\). For each degree make a local variable
change \(z=\xi_j x\) with \(\xi_j^j=b_{j,0}/b_{j,j}\), then divide
by \(b_{j,0}\). This changes no root quotients.
The endpoints become one. Every internal coefficient has
strictly positive valuation by the definition of \(P_p\),
so the reduction is \(1+x^j\). All roots are units. For at
least one of \(j=n,N\), \(p\nmid j\), so that reduction is
squarefree in that degree. Its off-diagonal quotients reduce
to nontrivial \(j\)-th roots of unity. Quotients in the other
degree reduce to roots of unity of that other order,
possibly to one if its reduction is inseparable.
No equality is possible since \(\gcd(n,N)=1\).
No assumption on \(p\) dividing the other degree is needed.

## Negative clock valuation

Suppose \(v_p(q)=-s\) with \(s>V_p\).
The valuation at coefficient index \(k\) is
\(a_{j,k}-sK_jk(j-k)\). Its second difference is
\(d_{j,k+1}-d_{j,k}+2sK_j>0\) by the definition of \(V_p\).
Thus all lower Newton polygon edges have width one.
In index order \(1\le k\le j\), the root valuations are
\[
                 sK_j(j+1-2k)-d_{j,k}.
\]
This is the Newton polygon root rule with coefficient indices
read in increasing order. The rule follows by factoring into
linear factors in a valued algebraic closure. At each weight the
minimum weighted coefficient valuation of a product is the sum
of the minima for its factors, since their initial polynomials
multiply over the residue field.
Valuations of individual roots are all distinct within an
order. Valuations of their quotients need not be distinct.
If a quotient at indices \(k,l\) in order \(n\) equals one
at indices \(u,v\) in order \(N\), valuation equality requires
\[
 2s[N(l-k)-n(v-u)]
            =(d_{n,k}-d_{n,l})-(d_{N,u}-d_{N,v}).                 \tag{2}
\]
The bracket on the left is a nonzero integer: a zero would imply
\(n\mid l-k\), impossible for two different lower indices.
The right side has modulus at most \(2V_p\), which contradicts
\(s>V_p\). This proves (1).

Even for \(s\le V_p\), if the actual ordinates are strictly convex,
(2) enumerates a finite list of candidate negative valuations
from the core first differences. At any such valuation, a further
necessary congruence for an equal quotient is
\[
   \frac{A_{n,k-1} A_{n,l}}{A_{n,k} A_{n,l-1}}
       \bigg/
   \frac{A_{N,u-1} A_{N,v}}{A_{N,u} A_{N,v-1}}
                       \equiv1\pmod{\mathfrak m},\qquad
                 A_{j,k}=b_{j,k}q^{K_j k(j-k)} .
\]
Here \(\mathfrak m\) is the maximal ideal in a valued extension
containing \(q\). To check it, at the width-one edge from \(k-1\)
to \(k\), divide the root by \(-A_{j,k-1}/A_{j,k}\).
After rescaling the polynomial by the lowest power in its
valuation group, exactly the two edge terms survive. The
rescaled root therefore reduces to one. Taking ratios gives
the congruence even when several quotients have the same
valuation. A shared residue need not lift to an equal quotient.

The same negative-valuation test allows independent clock
parameters \(q_n,q_N\) at a good prime, using the coefficient
multipliers \(q_j^{k(j-k)}\) in degree \(j\). If their valuations
are both negative, equality requires the ratio of their
absolute valuations to lie among \(b/a\) with
\(1\le a\le n-1\), \(1\le b\le N-1\).
If just one valuation is negative and the other is nonnegative,
the off-diagonal quotients of the negative-clock polynomial
have nonzero valuations and the other-clock quotients are units.
This excludes common arguments and gives coprimality. For common strength as above,
the unequal factors \(K_n=N,K_N=n\) enforce precisely the
arithmetic separation in (2).

These comparisons concern the exact Gaussian-deformed pair.
The zero-strength ordinary specialization \(q=1\) does not
satisfy either strict valuation hypothesis. Nor do these
inequalities eliminate all possible positive algebraic clock
parameters at once.
