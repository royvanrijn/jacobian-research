# Cofinite all-degree angular separation on every ordinary row fiber

Fix any ordinary two-monomial profile, choose one unequal row, and
hold every other row fixed. On the chosen row fix the positive gap
\(d\) and its orientation, and let the smaller exponent \(L\)
vary over the nonnegative integers. Write \(J\ge4\) for any
[absolute ordinary angular cutoff](MULTIVARIATE_EFFECTIVE_ABSOLUTE_EVENTUAL_ANGLES.md),
so terminals \(N\ge J\) have disjoint consecutive root arguments.
Put \(M=J-2\).

**Theorem.** At most
\[
          B_d=d\,\frac{M(M+1)(M-1)(M+2)(2M+1)}{10}
\]
values of \(L\) fail all-degree angular separation. At most
\(d M(M+1)(M+2)/6\) values fail SFC(2). Both bounds are uniform
in all the other exponent rows, including private rows. The value
of \(J\) is computable and unevaluated.

More precisely, for each pair \(2\le n<N\), the equality of a
nonreal root argument forces the vanishing of a nonzero integer
polynomial in \(L\) of degree at most \(d n(n-1)N(N-1)\).
For opposite parity the nonreal arguments give the entire
angular comparison. The same polynomial remains nonzero modulo
any prime \(p>\max(N^2 d,N a_*)\), where \(a_*\) is the largest
entry among the fixed rows, or zero if none remain.

In particular an extension of any fixed core by a single shared
unit-gap row is angularly regular at every order for all sufficiently
large smaller exponents. Its added curvature tends to zero. A fixed
row of opposite gap and orientation can be included in the core
when an extension preserving the original total-degree difference
is desired.

## Degenerating the upper moment at a negative baseline

Swap the two columns if needed so that the varying row is
\((L,L+d)\). At nonnegative integer \(L\), strip the nonzero
row scalar \((jL)!\). Use the resulting polynomial family
\[
 F_j(L,z)=\sum_{k=0}^j b_{j,k}(jL+1)_{dk}z^k,\qquad
 (x)_e=x(x+1)\cdots(x+e-1).
\]
Here \(b_{j,k}>0\) are the fixed core moment coefficients, including the
binomial factor. Set \(L_0=-1/N\), and use a formal parameter \(s\)
by \(L=L_0+s^N\). At degree \(N\), the factor \(NL+1=Ns^N\)
occurs once in each coefficient with \(k\ge1\). Thus, putting
\(z=u/s\),
\[
 F_N(L_0+s^N,u/s)
 = b_{N,0}+D_Nu^N+sD_{N-1}u^{N-1}+O(s^2),
 \qquad D_k=N b_{N,k}(dk-1)! .
\]
The constant polynomial in \(s\) has \(N\) simple nonzero roots
\(\lambda_r\). Solving the power series equation successively at
each root gives
\[
        u_r(s)=\lambda_r+c\,s+O(s^2),\qquad
        c=-\frac{D_{N-1}}{ND_N}\ne0 .
\]
For distinct \(r,t\) the ratio \(u_r(s)/u_t(s)\) has a nonzero
coefficient of \(s\), namely \(c(\lambda_t-\lambda_r)/\lambda_t^2\).

The lower polynomial \(F_n(L_0,z)\) has nonzero endpoints because
\(1-n/N>0\). Every lower root \(\xi\) of \(F_n(L_0+s^N,z)\)
has valuation zero for the valuation \(v_s(s)=1\), since its
endpoints are units and its coefficients are integral. Work
in a valued algebraic closure of the Laurent series field
over the algebraic closure of the coefficient field. Reduction
gives a nonzero root \(\xi_0\) of \(F_n(L_0,z)\). If its
multiplicity is \(m\le n\), then
\[
 v_s(\xi-\xi_0)\ \ge\ N/m .
\]
Indeed \(F_n(L_0,\xi)=O(s^N)\) and its other constant root
factors have valuation zero at \(\xi\). Every ratio of two
lower roots therefore differs from its constant value by
an element of valuation strictly greater than one. It
cannot agree with any off-diagonal upper ratio, including
when the constant ratios agree. At positive characteristic
the identical argument holds provided the displayed endpoints
and \(ND_N D_{N-1}\) are nonzero and \(N\) is invertible.
The stated prime bound ensures these conditions. It bounds
all the factorial arguments in the fixed coefficients and
all positive numerators \(Nh-n\) for \(1\le h\le nd\).

## An integer polynomial for the ratio obstruction

For degree \(j\) let \(f_k=[z^k]F_j(L,z)\) and let \(z_r\)
be its roots counted with multiplicity. Put \(m_j=j(j-1)/2\).
Use the reciprocal expression for unordered root pairs
\[
 A_j(L,v)=f_j^{\,2(j-1)}
       \prod_{r<t}\big(v z_r z_t-z_r^2-z_t^2\big).
\]
Its leading coefficient in \(v\) equals
\((f_0 f_j)^{j-1}\). In each coefficient the symmetric root
polynomial is homogeneous of degree \(j(j-1)\) and has degree
at most \(2(j-1)\) in each root. When expressed in elementary
symmetric functions, it thus has weighted degree \(j(j-1)\)
and at most \(2(j-1)\) elementary factors per term.
To see the second bound, successively subtract the product
of elementary functions with a given leading lexicographic
root monomial. If its nonincreasing exponent vector is
\(\lambda\), this product uses \(\lambda_1\) elementary factors,
each obtained from one column of the partition \(\lambda\).
The subtracted product has degree at most \(\lambda_1\)
in each root, so the per-root cap persists throughout.

Substitute \(e_\ell=(-1)^\ell f_{j-\ell}/f_j\), where \(e_\ell\)
is the elementary symmetric function of degree \(\ell\).
Because \(\deg_L f_k=dk\), the prefactor clears all
denominators and gives
\[
        A_j\in\mathbb Z[L,v],\qquad
        \deg_v A_j=m_j,\qquad
        \deg_L A_j\le d j(j-1).
\]
For the last bound each term has the same bound
\(2dj(j-1)-d\sum\ell=dj(j-1)\), with \(\ell\) counted by
its multiplicity among the elementary factors.

Two inverse pairs have the same value of \(v=q+1/q\)
exactly when their ratios agree or are reciprocal. Thus the
negative-baseline calculation shows that
\[
                \mathcal R_{n,N}(L)=
                      \operatorname{Res}_v(A_n,A_N)
\]
is nonzero, including at the stated good primes. Its degree
in \(L\) is at most \(d n(n-1)N(N-1)\) by resultant
homogeneity. At a nonnegative ordinary baseline, roots \(z,w\)
of the two polynomials on the same nonreal ray yield the same
ordered ratios \(z/\bar z=w/\bar w\). So the resultant
vanishes. Opposite parity precludes a common real ray,
since the even moment is strictly positive on the real line
as the expectation of an even power of a nonzero sum of
distinct real monomials.

For adjacent pairs \(N=n+1<J\) the degree bound is
\(d(n^4-n^2)\). The degree-one pair has disjoint arguments
by the same even-moment positivity. Sum the degree bounds
and use the absolute cutoff for later pairs to obtain \(B_d\).
There is also a modular version. For any prime
\(p>\max(J^2d,J a_*,B_d)\), all but at most \(B_d\)
residue classes of \(L\) modulo \(p\) give all-degree
ordinary angular separation whenever every ratio resultant
is nonzero on that class. Each of the finitely many
polynomials remains nonzero over \(\mathbb F_p\), so their
combined zero set has the stated size.

## The smaller bound for actual common roots

The ordinary resultant \(\operatorname{Res}_z(F_n,F_N)\)
is nonzero already by specialization at \(L_0\). The upper
polynomial becomes a nonzero constant, and the lower one
retains its full degree. The homogeneous resultant of degrees
\(n,N\) at this specialization is the constant to the
\(n\)-th power times the lower leading coefficient to the
\(N\)-th power, up to sign.

This also gives a direct certificate for a given ordinary baseline.
If there is a prime \(p\) such that
\[
 p\mid NL+1,\qquad
 p\nmid N b_{N,0} b_{n,n}\prod_{h=1}^{nd}(Nh-n),
\]
then the homogeneous resultant is a unit modulo \(p\), since
the same constant-versus-full-degree specialization holds in
that residue field. No degrees or bases other than those of
the tested pair enter this certificate.

The degree in \(L\) is at most \(dnN\). Indeed each
resultant monomial has coefficient-index weight \(nN\)
under simultaneous scaling \(z\mapsto vz\) of both inputs.
For consecutive degrees the resultant is nonnegative on
a real open interval about each ordinary baseline, since
the even polynomial there has positive leading coefficient,
no real roots, and a conjugate pairing of all roots.
The factorial representation on that interval uses real
exponents \(L,L+d\) for one independent exponential
variable, with each moment integral finite when \(NL>-1\).
The two monomials remain distinct because the selected row
has gap \(d>0\). The resultant product over the roots of
the even polynomial is therefore a product of squared
absolute values times a positive leading factor.

Every zero on the ordinary fiber has multiplicity at least
two. Summing half the degree bounds for the remaining
consecutive pairs, with the degree-one pair already excluded,
gives at most \(d M(M+1)(M+2)/6\) exceptions.
The uniform top-degree proof for translating all rows and the
leading coefficient for partial translations are separate
[exceptional-locus estimates](MULTIVARIATE_ALL_ORDER_EXCEPTIONAL_LOCUS.md#partial-translations-and-the-smaller-support).
The fiber argument here holds with the other rows fixed
arbitrarily, even if the leading partial-translation
coefficient vanishes.

## Varying an entire reflected row pair

One may instead vary two rows \((L,L+d),(L+d,L)\) together
while fixing an arbitrary core. At most
\[
                  d\,\frac{M(M+1)(M+2)}{3}
\]
ordinary values of \(L\) fail SFC(2) at some degree. If the core
is mirrored, the same bound holds for all-degree angular
separation. This version preserves exact mirroring as well as
the difference of the total degrees.

At degree \(j\) the appended coefficient multiplier, after
stripping the two factorial scalars, is
\((jL+1)_{dk}(jL+1)_{d(j-k)}\).
In any adjacent pair let \(E\) be its even member, assume
\(E\ge4\), and put
\[
                        h=dE/2+1,\qquad L_0=-h/E .
\]
At order \(E\), the factor \(EL+h\) kills every term except the
one with \(k=E/2\). That term survives, so the specialized
homogeneous polynomial is a nonzero multiple of
\(x^{E/2}y^{E/2}\), whose projective roots are the coordinate
axes. The other order \(j\) is relatively prime to \(E\).
Also \(h\) is not a multiple of \(E\), since its residue is
\(1\) or \(E/2+1\). Therefore \(jL_0\) is nonintegral,
so both pure coefficients at order \(j\) are still nonzero
and the two forms have no common projective point.

This proves that the ordinary consecutive resultant is nonzero
on every reflected-pair fiber for \(E\ge4\). The pairs with
\(E=2\) are already regular for all ordinary baselines by the
[endpoint nondivisibility argument](MULTIVARIATE_TWO_MONOMIAL_NONDIVISIBILITY.md#an-inequality-using-the-endpoint-coefficients).
Each coefficient on the fiber has \(L\)-degree at most \(dj\).
Homogeneity gives resultant degree at most \(2dnN\), and
even-moment positivity gives real zero multiplicity at least two
as above. The degree sum and the absolute cutoff prove the
stated bound.

For a mirrored core, the
[fixed-degree mirrored circle theorem](MULTIVARIATE_MIRRORED_ROOT_CIRCLE.md)
places every ordinary root on the unit circle. Equality of
arguments then entails equality of roots, which proves the
angular conclusion with the same smaller bound.
