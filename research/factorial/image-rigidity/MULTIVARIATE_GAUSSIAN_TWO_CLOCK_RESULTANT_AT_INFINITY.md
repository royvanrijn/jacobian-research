# Independent Gaussian clocks and the angular resultant at infinity

Let \(n\ge2\), \(N=n+1\). In each degree \(j=n,N\) start with
a polynomial of degree \(j\) with positive real algebraic
coefficients \(b_{j,k}\), and form
\[
    F_j(q,z)=\sum_{k=0}^j b_{j,k}q^{k(j-k)}z^k .
\]
The input can be any exact ordinary factorial core, and the two
inputs need not be the same core. The Gaussian strengths giving
these clocks are \(-2j\log q\) for \(0<q\le1\). This note
determines the highest-weight part of the angular resultant in
the two clocks. It also separates nonreal root arguments for
mixed arithmetic types of clock, including a zero strength in
one order and a positive algebraic strength in the other.

For roots \(z_{j,r}\) of \(F_j(q,\cdot)\), with multiplicities,
define the polynomial in the unordered quotient coordinate \(v\):
\[
 B_j(q,v)=\prod_{r>s}
  \bigl(v-z_{j,r}/z_{j,s}-z_{j,s}/z_{j,r}\bigr),\qquad
 \mathcal R(X,Y)=\operatorname{Res}_v(B_n(X,v),B_N(Y,v)).
\]
These are polynomials over the algebraic numbers. Indeed, after
multiplication by \((\prod_r z_{j,r})^{j-1}\), the coefficients
of \(B_j\) are symmetric polynomials in its roots. The multiplier
and the leading coefficient of \(F_j\) are nonzero constants.
A nonreal same-angle pair forces \(\mathcal R=0\), since for
a root \(z\) its conjugate root supplies the coordinate
\(z/\bar z+\bar z/z\). The use of unordered pairs is equivalent
to passing from the ordered quotient polynomial to \(v=u+u^{-1}\).
If one uses the ordered polynomial in both degrees, its
resultant is \(\mathcal R^2\). This follows from the ordered
factorization \(u^{j(j-1)/2}B_j(q,u+u^{-1})\) and the product
over reciprocal pairs.
[Ordered quotient construction](MULTIVARIATE_ARBITRARY_CORE_DIFFUSE_EXTENSION.md#an-angular-eliminant-for-the-exact-core)

## Exact leading terms and the two-clock degrees

At \(q\to+\infty\), the roots can be indexed as
\[
 z_{j,k}(q)=-\frac{b_{j,k-1}}{b_{j,k}}\,
                q^{2k-j-1}\bigl(1+O(q^{-2})\bigr),
                         \qquad 1\le k\le j.                 \tag{1}
\]
To see this, substitute \(z=q^{2k-j-1}w\) and divide
by \(q^{k(k-1)}\). The summand with index \(\ell\) then has
exponent \(-(\ell-k)(\ell-k+1)\) in \(q\). Only \(\ell=k-1,k\)
have exponent zero, giving \(w^{k-1}(b_{j,k-1}+b_{j,k}w)\).
Its nonzero root is simple, so a disk around that root gives
one nearby root with error \(O(q^{-2})\). Their scales are
different for all \(k\), accounting for every root.

Put \(m_j=j(j-1)/2\) and
\[
 D_j=\frac{j(j^2-1)}3,\qquad
 C_j=(-1)^{m_j}\frac{\bigl(\prod_{k=1}^{j-1}b_{j,k}\bigr)^2}
                              {(b_{j,0}b_{j,j})^{j-1}} .
\]
By (1), an unordered quotient at \(r>s\) is asymptotic to
\(c^j_{rs}q^{2(r-s)}\), where
\(c^j_{rs}=b_{j,r-1}b_{j,s}/(b_{j,r}b_{j,s-1})\).
Consequently \(B_j(q,v)\) has exact degree \(D_j\) in \(q\),
and its leading coefficient in \(q\) is \(C_j\), independent
of \(v\). Indeed the exponent is
\(2\sum_{r>s}(r-s)=D_j\) and the corresponding coefficient
is the product of \(-c^j_{rs}\). This holds for every fixed
complex value \(v\). Comparing at infinitely many \(v\)
gives the claim in the polynomial ring.

It follows that \(\mathcal R(X,Y)\) has \(Y\)-degree
\(m_nD_N\) with leading coefficient \(C_N^{m_n}\), constant
and nonzero in \(X\). Evaluate the resultant at the roots
of the monic \(B_n(X,v)\), and use the leading term just
computed for \(B_N\) at each such root. In the opposite
direction the \(X\)-degree is \(m_ND_n\), with leading
coefficient \((-1)^{m_nm_N}C_n^{m_N}\), constant in \(Y\).
These identities include clocks where an input polynomial
has repeated roots. Each fixed value of either clock thus
leaves only finitely many exceptional values of the other.

More generally, at positive weights \(a,b\) the leading
weighted form of \(\mathcal R\) is exactly the product, over
\(r>s\) and \(u>v\), of the leading weighted forms of
\[
       c^n_{rs} X^{2(r-s)}-c^N_{uv}Y^{2(u-v)}.                \tag{2}
\]
Here “leading” means largest weight, retaining both terms
when their weights agree. In fact substitute
\(X=X_0 t^{-a}, Y=Y_0 t^{-b}\) with \(X_0,Y_0>0\) and
\(t\downarrow0\), and apply (1) to the resultant product.
Take \(X_0,Y_0\) away from the finite list of binomial
equalities in (2), so every leading factor is nonzero.
The resulting asymptotic identifies the weighted degree
and leading coefficient for all such \(X_0,Y_0\), an open
set. Equality of polynomials gives (2), including its
values on the excluded binomial curves. The exposed
Newton faces at positive weights are thus determined
by the degree pair, while their coefficients record the
actual adjacent ratios of the core.

For the common strength \(X=q^N,Y=q^n\), the exponents
in each factor of (2) are distinct. In fact \(Nd=ne\)
would force \(n\mid d\), where \(d=r-s\) is in
\(\{1,\ldots,n-1\}\). Therefore the exact common-clock
resultant degree is
\[
 \begin{aligned}
 \deg_q\mathcal R(q^N,q^n)
   &= 2\sum_{d=1}^{n-1}\sum_{e=1}^n
                   (n-d)(N-e)\max(Nd,ne)\\
   &=\frac{n(n-1)(n+1)(14n^3+22n^2-n+2)}{60}.
 \end{aligned}                                               \tag{3}
\]
The evaluation of the finite sum separates at \(e=d\),
since \(Nd>ne\) exactly for \(e\le d\).
This gives an explicit upper bound for the number of
candidate positive common strengths at every finite core
and terminal. Actual same-angle pairs are a subset of
the resultant witnesses.

## Mixed clock noncollision

If either \(X\) or \(Y\) is algebraic and the other is
transcendental, their \(\mathcal R(X,Y)\) is nonzero. This
follows directly from the nonzero finite-degree fiber
over the algebraic numbers proved above.

Suppose instead that \(X,Y>0\), \(X\) is transcendental,
and \(Y=\gamma X^{r}\), where \(\gamma>0\) is algebraic
and \(r\) is rational. The same conclusion holds. For
\(r=0\) it is the previous case. Otherwise write \(r=u/v\)
with \(v>0\) and put \(t=X^{1/v}\), which is transcendental.
The polynomial \(\mathcal R(t^v,\gamma t^u)\) is a Laurent
polynomial over the algebraic numbers, and is nonzero:
for \(u>0\), both clocks tend to zero as \(t\downarrow0\).
The binomial quotient sets at zero are disjoint because
\(n,N\) are coprime. In the unordered coordinate, equality
of the values \(\xi+\xi^{-1}\) only identifies \(\xi\)
or \(\xi^{-1}\), so the same disjointness holds. For
\(u<0\), the lower quotient coordinates stay bounded
while all upper ones diverge by (1). Thus no root
quotients overlap for small enough positive \(t\) in
either case. Nonzero Laurent polynomials with algebraic
coefficients cannot vanish at a nonzero transcendental
point.

For \(0<X,Y\le1\) these conditions exclude every nonreal
common argument. They exclude the real directions as well
if the even input is an ordinary factorial moment on
distinct monomials. The even input is positive on the
real axis by the even-power expectation, and attenuation
preserves this positivity by lognormal averaging. More
generally, a nonnegative even input on the real axis is
strictly positive after any strictly positive attenuation.
These conclusions do not assume simple roots in the
actual clock range.

This comparison concerns fixed algebraic inputs, including
arbitrary finite ordinary cores. No claim is made that a
countably retained compactified core has algebraic
coefficients.

