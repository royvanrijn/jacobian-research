# Adjacent collisions for smooth positive moment laws

For every integer \(N\ge5\), every \(c>0\) rational and
every \(\epsilon_0>0\), there is a positive random variable \(X\)
with a smooth compactly supported strongly log-concave density,
supported within \(\epsilon_0\) of \(c\), such that the polynomials
\[
                         B_j(z)=\mathbb E(1+zX)^j
\]
have rational coefficients through \(j=N\), and \(B_{N-1},B_N\)
have exactly one common irreducible quadratic over \(\mathbb Q\).
Its roots lie in \(\mathbb Q(i)\) and are simple at order \(N-1\)
and double at order \(N\). In particular the relative variance
can be made arbitrarily small. The same positive law supplies
the moments at both orders and at all the intervening orders.

The weighted Hankel forms of this moment sequence are all
strictly positive: for every nonzero real polynomial \(A\) and
every integer \(r\ge0\), \(\mathbb E[X^r A(X)^2]>0\).
Consequently strict moment positivity, strong log-concavity,
concentration, and rational coefficients together do not give
coprimality for general positive integral analogues of the
factorial moment. The construction uses the adjacent derivative
identity from the identical positive law. This is a different
relaxation from the
[common Laplace-symbol collision](MULTIVARIATE_MATCHED_TANH_ANGULAR_HIERARCHY.md#an-exact-collision-for-an-unrestricted-common-laplace-symbol).

## Losing imaginary-rootedness along a log-concave path

For a symmetric real random variable \(Y\) put
\(A_j(w)=\mathbb E(w+Y)^j\). The exact identity is
\(A_N'=NA_{N-1}\). Write \(m=\lfloor N/2\rfloor\), and put
\(e=N-2m\). By symmetry
\[
               A_N(w)=w^e Q(w^2),\qquad
 Q(y)=\sum_{j=0}^m\binom N{2j}\mathbb E[Y^{2j}]\,y^{m-j}.
                                                               \tag{1}
\]
For a positive symmetric density on an interval all coefficients
of \(Q\) are positive and its leading coefficient is one.

For \(Y\) uniform on \([-L,L]\), one has
\[
     A_N(w)=\frac{(w+L)^{N+1}-(w-L)^{N+1}}{2L(N+1)} .
\]
Its \(N\) roots come from
\((w+L)/(w-L)\) equal to the nontrivial \((N+1)\)-st roots
of unity. They are distinct and purely imaginary. The root at
zero occurs only when \(N\) is odd, and is then simple.
Thus \(Q\) has only simple strictly negative roots.

For the symmetric Laplace density proportional to \(e^{-|y|}\)
on the real line, \(\mathbb E[Y^{2j}]=(2j)!\), and the coefficient
of \(y^k\) in \(Q\) is \(N!/(2k+e)!\). If all roots of \(Q\)
were real, their negativity and Newton's inequality at index
\(m-1\) would require
\[
     \frac{N(N-1)}{(N-2)(N-3)}\ \ge\ \frac{2m}{m-1}.           \tag{2}
\]
Indeed if the roots were \(-b_i\) with \(b_i>0\), the two
coefficients below the monic term would be \(\sum_i b_i\) and
\(((\sum_i b_i)^2-\sum_i b_i^2)/2\). Cauchy–Schwarz gives
(2). The strict reverse inequality holds for every
\(N=2m\) with \(m\ge3\) and every \(N=2m+1\) with \(m\ge2\).
The reductions are \(2m-1<2(2m-3)\) and
\(2m+1<2(2m-1)\), respectively. Therefore this \(Q\)
has a nonreal root.

To join the two examples inside the required smooth class,
choose \(L\) so large that the Laplace law truncated to
\([-L,L]\) and renormalized still violates (2).
This is possible by convergence of its finitely many moments.
On \((-L,L)\) use the normalized densities
\[
 f_{a,\epsilon}(y)\ \propto\
      \exp\!\left(-a\sqrt{y^2+\epsilon^2}
                   -\epsilon y^2
                   -\frac{\epsilon}{L^2-y^2}\right),
        \qquad 0\le a\le1 .
\]
Extend them by zero outside that interval. They are smooth
and symmetric. The negative logarithm has second derivative
at least \(2\epsilon\) throughout the support interior,
since the other two displayed summands of the potential are
convex. For \(\epsilon\downarrow0\) the endpoint \(a=0\)
tends in moments to the uniform law, and the endpoint \(a=1\)
to the truncated Laplace law. Fix a sufficiently small positive
\(\epsilon\). The first endpoint still has simple negative
roots of \(Q\), by continuity of simple real roots, while the
second still violates (2).

Follow the connected component of the first endpoint in the
parameter set where all roots of \(Q\) are simple and negative.
At its last parameter the roots are all real by continuity of
the root multiset. They are strictly negative because the
coefficients of \(Q\), including its constant term, are positive.
At least one root is repeated, since simple negative roots
persist on a neighborhood. This supplies a symmetric smooth
strongly log-concave density with a double or higher root
of \(Q\) at some \(-a_0<0\).

## Rational moments and one common pair

One can perturb the last density while staying within this
class and prescribe all its first \(m\) even moments in a
neighborhood. Multiply the density by
\(\exp(\sum_{j=1}^m t_j y^{2j})\) and renormalize.
At \(t=0\) the Jacobian of those moments is the covariance
matrix of \(Y^2,\ldots,Y^{2m}\). It is positive definite:
a nontrivial linear combination cannot be constant on an
interval with positive density. The inverse function theorem
therefore gives the asserted neighborhood. Small tilts preserve
the strong convexity bound on the negative log density,
since the second derivatives of the tilting polynomials are
bounded on the support.

By (1), these moment coordinates freely parametrize the
nearby monic polynomials \(Q\). Write the original polynomial
as \((y+a_0)^2R(y)\). Approximate it by
\[
                         \widetilde Q(y)=(y+t^2)^2\widetilde R(y),
          \qquad t\in\mathbb Q_{>0},
\]
where \(t^2\) is close to \(a_0\) and \(\widetilde R\) is monic
rational, squarefree, and nonvanishing at \(0,-t^2\).
Such choices are dense near the original factorization.
The moment neighborhood realizes this exact rational polynomial.
The resulting \(A_N\) has exactly two repeated roots,
\(\pm it\), both of multiplicity two. In the odd case
its factor \(w\) is simple because \(\widetilde Q(0)\ne0\).
Thus the adjacent Appell gcd is exactly \(w^2+t^2\).

Now set \(X=c+dY\) for a sufficiently small positive rational
\(d\). The support lies in \((0,\infty)\) and within
\(\epsilon_0\) of \(c\). Its density remains smooth and strongly
log-concave. All moments of \(Y\) through \(N\) are rational,
so the same holds for \(X\). The affine transform carries
the repeated Appell pair to \(-c\pm idt\), and reciprocating
the Appell variable gives \(B_j\). This recovers the single
common quadratic
\[
                        1+2cz+(c^2+d^2t^2)z^2.
\]
It is irreducible over \(\mathbb Q\). Its multiplicities are
preserved through the affine and reciprocal transformations.
The assertion about the Hankel forms follows from positivity
of \(X\) and of its density on an interval.

## The first possible degree in the symmetric construction

For a symmetric nondegenerate log-concave density the fourth
moment satisfies \(\mathbb E Y^4\le6(\mathbb E Y^2)^2\).
For completeness, compare the density of \(|Y|\) on the
positive half-line with an exponential density of the same
second moment. The difference is positive on an interval
and nonpositive outside, since the logarithm of the density
ratio is concave. If the densities differ, matching their mass
and second moment forces two crossings. With just one crossing
at \(a\), multiplication of the difference by \(y^2-a^2\)
would have one strict sign in the integral, contradicting
those equalities. With the two crossing points \(a,b\),
multiplying by \((y^2-a^2)(y^2-b^2)\) gives a nonpositive
integral. The equalities of mass and second moment leave
the desired fourth-moment comparison. The exponential
fourth moment is six times the square of its second moment.

For \(N=4\), repetition of \(Q(y)=y^2+6\mathbb E[Y^2]y+
\mathbb E[Y^4]\) would require
\(\mathbb E[Y^4]=9(\mathbb E[Y^2])^2\). This is impossible
under the comparison. At \(N\le3\) the symmetric
Appell polynomial is simple directly from (1).
Thus \(N=5\) is the first possible collision for the
translated symmetric log-concave construction.
