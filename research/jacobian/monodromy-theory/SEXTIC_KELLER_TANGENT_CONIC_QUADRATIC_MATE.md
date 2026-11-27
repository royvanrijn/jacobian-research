# Low free-variable degrees in the quartic tangent-conic scalar

**Theorem SexticTangentConicScalarLowMateExclusion.**
Let \((P,Q,R)\) be a characteristic-zero polynomial Keller
triple of ordinary degree at most six with
\(R=g(x)+x^2(xz+y^2)\). On substituting \(z=v-y^2/x\),
both mates have \(y\)-degree at least three.

Work over \(\mathbb C\) by embedding the coefficient field,
and subtract \(g(0)\) so \(g(0)=0\). The submersion condition
gives \(g'(0)\ne0\). Write \(r=g(x)+x^3v\).
The substitution is an isomorphism on \(x\ne0\), with determinant
one. In particular coefficients after substitution belong to
\(\mathbb C[x,x^{-1},v]=\mathbb C[r,x,x^{-1}]\).
The values at \(y=0\) in the new frame come from polynomials
in the original frame, since \(v=z\) there. The same holds for
the first \(y\)-derivative evaluated at \(y=0\).
On \(x=0\), the original mates give a plane Keller pair of
degree at most six with nonzero determinant. Indeed, \(dR\)
there is \(g'(0)\,dx\). This pair is a plane automorphism by
[Moh's bound](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md).

## Support and the affine cases

Group the original monomials \(x^{N+b}y^{d-2b}z^b\) by
\((N,d)\), leaving \(x^Ny^d S(t)\) with \(t=xz/y^2\).
The allowed exponent interval in \(S\) is
\(\max(0,-N)\le b\le\lfloor d/2\rfloor\), and \(N+d\le6\).
The substitution replaces \(t\) by \(xv/y^2-1\).
If the new \(y\)-degree is \(m\), a nonzero group must be
divisible by \((1+t)^k\) for \(k=\max(0,\lceil(d-m)/2\rceil)\).
Different groups give distinct powers of \(x\) for each fixed
pair of exponents of \(y,v\), so this vanishing occurs in each
group. Since \((1+t)\) is coprime to \(t\), the exponent interval
implies \(N\ge-\lfloor m/2\rfloor\). A leading \(y^m\)-coefficient
term of \(v\)-degree \(k\), coming from \(d=m+2k\), has
\(x\)-exponent between \(k-\lfloor m/2\rfloor\) and \(6-m-k\).

Therefore a mate of \(y\)-degree at most one restricts on
\(x=0\) to a polynomial in \(y\) alone, necessarily linear
with nonzero slope because it is a coordinate. If its
\(y\)-degree is zero, its original polynomial is even in
\(y\) on \(x\ne0\), hence everywhere, a contradiction.

For an affine mate \(Q=\alpha y+B\), replace \(y\) by \(Q\)
over the generic \(r\)-fibre after localizing at \(\alpha\).
The [affine-mate comparison](KELLER_REPEATED_PARALLEL_SCALAR_AFFINE_MATE.md#localization-by-the-affine-mate)
gives \(P=\Psi(Q,r)+V\) and \(\partial_x V=\kappa/(x^3\alpha)\)
at fixed \(r\), with \(\Psi\in\mathbb C(r)[Q]\).
Since \(V=P|_{y=0}-\Psi(B,r)\) is Laurent in \(x\) over
\(\mathbb C(r)\), \(\alpha=e(r)x^j\) for \(e\in\mathbb C[r]\).
The exponent bound gives \(4\deg e+j\le5\).
Also \(\alpha(0,v)\) is the nonzero slope of \(Q(0,y,z)\),
so \(j=-\operatorname{ord}_0 e\). These bounds imply
\(\deg e\le1\) and \(j=0\) or \(-1\).
Integration makes \(V\) an invariant plus a nonzero multiple
over \(\mathbb C(r)\) of \(x^{-2-j}\).
The generic degree of the original triple is thus \(2+j\).
Degree one contradicts the reducible fibre \(R=0\) by the
birational Keller theorem. Degree two is excluded by the
[classical Galois case](../../../unpublished/common-arithmetic-fibers/LITERATURE_AUDIT.md#classical-rank-two-exclusion).

## The quadratic case

Suppose now \(Q=qy^2+By+C\).
The [quadratic-mate argument](KELLER_REPEATED_PARALLEL_SCALAR_QUADRATIC_MATE.md#regularity-on-the-generic-scalar-fibre)
through its generic-fibre regularity step takes place over
\(\mathbb C(r)[x,x^{-1}]\). It applies to the Laurent
coefficients here, with \(q=e(r)x^{2j}\) and \(e\in\mathbb C[r]\).
The integration and unit comparison in the
[monomial-coefficient sheet formula](KELLER_PARABOLIC_SHEAR_QUADRATIC_MATES.md#generic-sheet-degree)
likewise give \(K=\tau+u x^{-h}\), \(sh=2+j\), \(cu^s=\lambda\)
and sheet degree \(|h|(2s-1)\), with \(K=C-B^2/(4q)\).
For \(t=\deg e\), the support bounds at the highest \(v\)-degree
and at \(v=0\) give
\[
                 4t+2j\le 4,\qquad \operatorname{ord}_0 e+2j\ge-1 .
\]
Thus \(t\le1\). If \(t=1\), we have \(j=0\), or \(j=-1\)
and \(e\) proportional to \(r\). If \(t=0\), then \(j=0,1,2\).
The case \(j=-1\) has sheet degree one, already excluded.

For the remaining cases \(q\) is regular along \(x=0\).
By the group bound, no term involving \(z^2\) survives in
\(Q(0,y,z)\). If its coefficient of \(z\) were nonconstant in
\(y\), a root of that coefficient would give a reducible
fibre of this plane coordinate. A nonzero constant
coefficient \(a\) of \(z\) at \(x=0\) would give a term
\(-a x^{-1}\) in \(q\) at \(v=0\): this is the group with
\(N=-1,d=2\), in which \(S\) is a multiple of \(t\).
Therefore \(Q(0,y,z)\) is independent of \(z\), and is
linear with a nonzero slope. In particular \(B\) is a unit
at the generic point of \(x=0\).

Write \(H=r/x\). Off \(x=0\), the divisor \(H=0\) has
parameter \(x\) and uniformizer \(r\). If \(e(0)\ne0\),
\(q\) is a unit there. Regularity of \(K\) gives
\(\operatorname{ord}_0\tau,\operatorname{ord}_0u\ge0\).
Since the centered coordinate \(Z=x^j(y+B/(2q))\) is regular
with unit slope there, the leading odd coefficient of \(P\),
proportional to \(c e^{s-1}\) in \(Z\), is regular. Thus
\(cu^s=\lambda\) gives \(\operatorname{ord}_0u=0\).
At \(x=0\), evaluated with \(v\) fixed, \(K\) then has
pole order \(h>0\). The formula \(C-B^2/(4q)\) with
\(B\) a unit forces \(h=2j>0\). Since \(s h=2+j\)
and \(0\le j\le2\), only \(j=2,s=1\) is possible.
Here \(e\) is constant. The coefficients of the even
part \(\Psi(Q,r)\) are polynomial in \(r\), since
the centered calculation uses \(\mathbb C[r,x,x^{-1}]\),
takes \(e\) constant, and its invariant subring is
\(\mathbb C[r]\). The odd part is \(c Z\).
At \(y=0\) it equals \(c B/(2e x^2)\), which has a
pole at \(x=0\) because \(c(0)\ne0\).
This contradicts polynomiality of \(P|_{y=0}\).

It remains that \(q=rx^0\) up to rescaling and \(j=0\).
In this case \(s=1\) would have degree two, so \(s=2,h=1\).
The odd part is \(2c((Q-\tau)Z-rZ^3/3)\).
For an odd mate degree it forces \(m=3\) and leading
\(y\)-coefficient proportional to \(cr\). This belongs
to \(\mathbb C[r]\) and \(4\deg(cr)\le6-3\), so it is
constant. The resulting odd valuation of \(c\) contradicts
\(cu^2=\lambda\).

For an even mate degree \(m=2b\ge4\), the leading
coefficient is \(\psi_b r^b\) and is polynomial in \(r\).
The support bound gives \(4\deg(\psi_b r^b)\le6-2b\),
so \(\psi_b=\gamma r^{-b}\) with \(\gamma\ne0\).
Now use the
[parabolic two-component pole comparison](KELLER_PARABOLIC_SHEAR_QUADRATIC_MATES.md#the-remaining-exponent-in-degree-seven)
with the same \(s=2,h=1\), but \(Z=y+B/(2q)\).
Its argument along \(H=0\) is unchanged because \(x\) is
a unit there. In particular \(K\) and \(\tau\) have a
simple pole, \(b_0=\operatorname{ord}_0 u\ge0\), and the
negative Laurent powers in \(\Psi+O_*\) cancel. The pole
order \(N=2b_0+2\) of \(O_*\) must satisfy \(N\ge2b\)
by the coefficient \(r^{-b}Q^b\), so \(b_0\ge1\).
At \(x=0\) make the comparison only after setting \(y=0\),
where both ordinary mates are polynomial in \(x,v\).
Then \(U=u/x\) is regular, \(r\) has order one, and (3) of
the parabolic comparison again gives the nonzero order-two
pole for matching signs and an order-\(N\) pole for
opposite signs. This excludes the final quadratic case.
