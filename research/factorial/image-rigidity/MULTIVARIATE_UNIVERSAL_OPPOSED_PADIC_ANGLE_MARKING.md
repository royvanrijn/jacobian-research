# Universal all-degree angular separation after adding two opposed rows

Fix any ordinary two-monomial exponent profile \(C\), with no restriction
on the dimension, signs, or private rows. Let \(J\ge4\) be any integer
at least the [effective absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md).
Fix integer gaps \(d_+,d_-\ge J-1\), and write
\(D=\max(d_+,d_-)\), \(d_0=\min(d_+,d_-)\).
Adjoin the rows
\[
             (A_+,A_++d_+),\quad (A_-+d_-,A_-),\qquad
                           A_+,A_-\ge0.
\]

**Theorem.** For any fixed prime \(p\) there is a computable
\(E=E(C,d_+,d_-,p,J)\) such that every pair of nonnegative integer
bases \(A_\pm\equiv-1\pmod{p^E}\) gives disjoint consecutive root
arguments in every degree. This holds after independent positive
rescalings at the two orders, so in particular the extended
profile satisfies SFC(2) over every characteristic-zero field.
If \(M\) is the largest core exponent and \(p>\max(D^2,DM)\),
one can take \(E=1\).

Along the common baseline \(A_+=A_-=A\), every sufficiently
large nonnegative integer \(A\) works. The threshold is
computable, and every exceptional integer base has \(A+1\)
supported only at primes up to \(\max(D^2,DM)\).

The extension can add any prescribed integer signed gap
\(\Delta=d_+-d_-\) at arbitrarily small curvature. At equal gaps
the rows form a mirrored pair, so the signed total gap and every
endpoint normalization of the argument are preserved exactly.
The thresholds use the existing absolute angular theorem.
The progression proof does not require searches through
the exact zeros of individual terminals.

## A prime-adic coefficient profile valid across degrees

Write \(b_{m,k}>0\) for the integer coefficient of \(z^k\) in the
core degree-\(m\) moment. Remove the two row factorial scalars
\((mA_+)!\) and \((mA_-)!\) to obtain the \(k\)-th coefficient
\[
 F_{m,k}(A_+,A_-)=b_{m,k}
       (mA_++1)_{d_+k}(mA_-+1)_{d_-(m-k)},                   \tag{1}
\]
where the rising product of length zero is one.
Fix a prime \(p\), write \(e_\pm=v_p(A_\pm+1)\), and let
\(L=\max_{1\le h\le D^2}v_p(h)\).
For \(2\le m<J\) and \(1\le h\le mD\), the identity
\(mA_\pm+h=m(A_\pm+1)+(h-m)\) gives, when \(e_\pm>L\),
\[
 v_p(mA_\pm+h)=
 \begin{cases}
  e_\pm+v_p(m),&h=m,\\
  v_p(h-m),&h\ne m .
 \end{cases}                                                 \tag{2}
\]
For the second case the nonzero integer \(h-m\) has absolute value
at most \(D^2\), and the first summand has larger valuation.
Since both gaps are at least \(m\), the factor with \(h=m\)
appears only in the decreasing row at \(k=0\), only in the
increasing row at \(k=m\), and in both at every internal index.

Let \(c_{m,k}\) be \(v_p(b_{m,k})\) plus the sum of the valuations
\(v_p(h-m)\) over the two products in (1), omitting \(h=m\) from
both sums. These numbers are fixed independently of the bases.
With \(t=k/m\) and \(a_{m,k}=v_p(F_{m,k}(A_+,A_-))\),
equation (2) gives for \(0<k<m\)
\[
 a_{m,k}-(1-t)a_{m,0}-t a_{m,m}
 =(1-t)e_++t e_-+v_p(m)+c_{m,k}-(1-t)c_{m,0}-t c_{m,m}.       \tag{3}
\]
Choose a positive integer \(E>L\) strictly larger than
\[
       \max_{\substack{2\le m<J\\0<k<m}}
       \left(\frac{m-k}{m}c_{m,0}+\frac{k}{m}c_{m,m}
                    -c_{m,k}-v_p(m)\right).
\]
Then (3) is strictly positive whenever both \(e_\pm\ge E\).
At \(p>\max(D^2,DM)\) the core factorials and binomial coefficient
are units, and so are all the factors in (2) with \(h\ne m\).
Hence \(L=0\), all \(c_{m,k}=0\), and \(E=1\) suffices.

## Binomial reductions separate root quotients

For each such order \(m\), independently choose a \(p\)-adic
algebraic change of variable \(z=\xi_m x\), with
\(\xi_m^m=F_{m,0}/F_{m,m}\), and divide by \(F_{m,0}\),
suppressing the bases in the coefficients.
Equation (3) says the resulting polynomial has integral coefficients
and reduces to \(1+x^m\). Its roots are units because it is monic
integral with unit constant coefficient. Therefore every quotient
of two roots is a unit and reduces to an \(m\)-th root of unity.
When \(p\nmid m\), the reduction is squarefree. In that case
quotients of distinct roots reduce to nontrivial \(m\)-th roots
of unity. The variable changes cancel from the quotients.

In consecutive degrees \(n,N=n+1\), at least one degree is coprime
to \(p\). A quotient of distinct roots from that degree cannot
equal a quotient of distinct roots from the other degree.
In fact its reduction would be a nontrivial root of unity of order
dividing both \(n\) and \(N\). This argument also applies when the
other reduction is inseparable.

An equality of physical arguments between nonreal roots \(z,w\)
would give \(z/\bar z=w/\bar w\). These are quotients of distinct
roots because the polynomials have real coefficients. The identity
can be tested at a \(p\)-adic place of the same splitting field.
A real argument cannot be common either, since the even factorial
moment is strictly positive on the real line. It integrates the
even power of a nonzero real binomial against the product
exponential measure. This establishes disjoint arguments for
\(2\le n<N<J\). The pair starting at one is also excluded by the
even moment's strict positivity, and the absolute angular theorem
covers \(N\ge J\). The progression assertion follows.

At the degrees covered by the prime-adic proof, at most one distinct
root can be shared under any nonzero complex rescaling between the
two polynomials. Indeed, two distinct shared roots would give a
common quotient of distinct roots. This statement also holds for a
prime dividing one of the degrees.

## A finite bound on exceptional bases

For the algebraic dependence along the common baseline
\(A_+=A_-=A\), further divide the coefficient array (1) by
\((mA+1)_{d_0\lceil m/2\rceil}\). One of the two rising products
always contains it. The resulting polynomial \(H_m(A,z)\) has
integer polynomial coefficients in \(A\) of degree at most
\(q_m=Dm-d_0\lceil m/2\rceil\). At equal gaps \(d\), this becomes
\(d\lfloor m/2\rfloor\). Neither endpoint coefficient vanishes
for \(A\ge0\).

Following the [angular root-quotient eliminant](MULTIVARIATE_ARBITRARY_CORE_DIFFUSE_EXTENSION.md#an-angular-eliminant-for-the-exact-core),
define the polynomial
\[
 D_m(A,u)=(u-1)^{\epsilon_m}
       \frac{(-1)^m\operatorname{Res}_z(H_m(A,z),H_m(A,uz))}
            {(u-1)^m},\qquad
 \epsilon_m=\begin{cases}1,&m\text{ odd},\\0,&m\text{ even}.\end{cases}
\]
The division by the monic polynomial \((u-1)^m\) is exact
over \(\mathbb Z[A,u]\). In the root product for the resultant
it removes the \(m\) diagonal pairs. With nonzero endpoints,
the remaining roots in \(u\) are the ordered off-diagonal
quotients, counted with multiplicity, and \(1\) with the added
odd-degree multiplicity. In particular the degree in \(u\) is
\(r_m=m(m-1)+\epsilon_m\). The Sylvester determinant for the
resultant in \(z\) gives
\(\deg_A D_m\le2m q_m\).
The exact division by a fixed monic polynomial in \(u\) does
not increase this bound.

For \(2\le n<N=n+1<J\), set \(R_n(A)=\operatorname{Res}_u(D_n,D_N)\).
Each \(R_n\) is nonzero. To prove this, choose a prime
\(p>\max(D^2,DM)\) and evaluate at \(A=p-1\).
Both residual binomials are squarefree, so the complex-root
polynomials are simple as well. Their off-diagonal quotient
sets are disjoint by the prime-adic comparison, and the
extra root \(1\) in the odd degree cannot occur among the
even degree's off-diagonal quotients. Since the endpoints
are nonzero, neither \(D_m\) loses degree at this evaluation.
This proves \(R_n(p-1)\ne0\).

Every angular exception for \(A\ge0\) at \(n,N\) comes from
a nonreal common ray and makes \(R_n(A)\) vanish at least to
second order as a polynomial in \(A\). To see the multiplicity,
the shared quotient \(u=z/\bar z\) has modulus one and
\(u\ne1\). If \(u\ne-1\), both \(u\) and \(u^{-1}\) are
distinct common roots of the quotient polynomials. If \(u=-1\),
the two ordered pairs \((z,\bar z)\) and \((\bar z,z)\)
give multiplicity at least two in each quotient polynomial.
Thus the gcd in \(u\) has degree at least two, giving
Sylvester matrix corank at least two. The determinant of a
polynomial matrix with corank at least two at a parameter
value vanishes there to order at least two. This follows by
constant invertible row operations making two entire rows zero
at that value, then factoring the parameter difference from
each of those rows.

Applying the Sylvester degree bound in \(u\) shows that the number
of exceptional nonnegative integer bases over all orders is at most
\[
 \sum_{n=2}^{J-2}
     \big(nq_n r_{n+1}+(n+1)q_{n+1}r_n\big).           \tag{4}
\]
One can bound all possible exceptional bases by taking a root bound
for the finite product of the nonzero \(R_n\).
The same proof of cofiniteness works on every ray
\(A_\pm=h_\pm T-1\) with fixed positive integers \(h_\pm\) and
integer \(T\ge1\). Use the undivided array (1) in the resultant,
and take \(T=p\) at any sufficiently large prime. The degree of
each coefficient in \(T\) is at most \(Dm\), so (4) holds there
with \(q_m=Dm\).

## Effect on the ordinary core

In the common-base family, the added curvature is
\[
                  \frac{d_+^2}{A+d_+}+\frac{d_-^2}{A+d_-}.
\]
The signed total gap increases by exactly \(\Delta=d_+-d_-\).
Let \(G_m\) be the positive scale of the argument that sets both
endpoint coefficients to one after a scalar normalization.
The scale for the extended core, divided by that of the core, is
\[
       \left(\frac{(mA+1)_{md_+}}{(mA+1)_{md_-}}\right)^{1/m}.
\]
For \(A>0\) its logarithm differs from \(\Delta\log(mA)\)
by at most \(|\Delta|\log(1+D/A)\). Cancel the common factors
in the quotient and use \(h/m\le D\) for each factor that remains.
The adjacent relative center is therefore multiplied in the
limit by \((N/n)^\Delta\). At \(\Delta=0\) even the individual
centers are exactly unchanged for every ordinary \(A\).

The two row log-factorials as a function of \(k\) have
total second derivative at most \((d_+^2+d_-^2)/(mA)\) for
\(A>0\), from \(\psi'(x+1)\le1/x\).
Subtract the endpoint chord and compare with the quadratic
of this constant curvature. The added log multiplier of the
normalized coefficient at \(k\) lies between
\(- (d_+^2+d_-^2)k(m-k)/(2mA)\) and zero.
The multipliers thus converge to one at each fixed order,
uniformly on a finite set of orders, while the added curvature
tends to zero.

In particular, taking equal gaps shows that ordinary profiles
with all-degree angular separation are dense in the ordinary
coefficient closure at fixed signed total gap, including the
relative endpoint centers. With prescribed nonzero signed gap,
the two rows can make arbitrarily weak extensions even at
any of the integer twists in the
[exceptional set of the original core](MULTIVARIATE_ARBITRARY_CORE_TWISTED_TAIL_STABILITY.md#the-finite-set-and-the-stability-constant).
A limit of separated extensions can still have a common
argument. The theorem concerns the exact extensions of every
core and does not assert coprimality of the original core.
