# Quadratic mates and the parabolic quartic scalar

**Proposition KellerMonomialCoefficientQuadraticSheetDegree.**
In a characteristic-zero polynomial Keller triple with
\(r=g(x)+x^a v\), \(a\ge2\), suppose a mate is quadratic in
the free variable \(y\), with leading coefficient
\(e(r)x^{2j}\), where \(j\in\mathbb Z\) and \(e\ne0\).
Then \(a+j-1\ne0\). For some positive divisor \(s\) of
\(|a+j-1|\), the generic sheet degree is
\[
                 \frac{|a+j-1|}{s}(2s-1).
\]
In particular for \(a=2\), the exponent \(j=0\) is impossible,
and \(j=1\) forces \(s=2\) and sheet degree three.

**Theorem KellerParabolicShearQuadraticMateExclusion.**
Let \((P,Q,R)\) be a characteristic-zero polynomial Keller
triple of ordinary degree at most seven, with quartic component
\(R=g(x)+x^2(z+y^2)\). In the polynomial source coordinates
\(x,y,v=z+y^2\), both mates have \(y\)-degree at least three.

## Generic sheet degree

Work over \(\mathbb C\) by embedding the finitely generated
coefficient field. With \(k=\mathbb C(r)\), use the
[quadratic-mate formula](KELLER_REPEATED_PARALLEL_SCALAR_QUADRATIC_MATE.md).
In its notation \(w=x^j\), \(Q=eZ^2+K\),
\(J_T=F\), and \(\partial_x J(K)=\kappa x^{-a-j}\) at fixed \(r\).
Here \(K\in k[x,x^{-1}]\). The derivative of a Laurent
polynomial has zero \(x^{-1}\)-coefficient, so \(a+j-1\ne0\).
Integration gives \(J(K)-U=\lambda x^{1-a-j}\) with
\(U\in k\) and \(\lambda\ne0\).

For each root \(\tau\) of \(J(T)-U\) over an algebraic closure
of \(k\), the factor \(K-\tau\) is a unit of the Laurent
polynomial ring. Two such units with nonzero constant difference
would both be constant, contrary to the derivative equation.
Thus this polynomial has one distinct root, and, with
\(s=\deg J\), we have
\[
 J(T)-U=c(r)(T-\tau(r))^s,\qquad
 K=\tau(r)+u(r)x^{(1-a-j)/s},\qquad c u^s=\lambda .
\tag{1}
\]
The root and the displayed coefficients belong to \(k\) by
coefficient comparison in characteristic zero. The degree of
\(K\) as a rational map in \(x\) is \(|a+j-1|/s\).
The quadratic-mate sheet-degree formula now gives the claim.

The scalar \(r\) has a reducible fibre at \(g(0)\), as
\(g'(0)\ne0\) by the submersion condition and
\(r-g(0)=x((g(x)-g(0))/x+x^{a-1}v)\).
Consequently the triple is not an automorphism. Generic degree
one is excluded by the birational Keller theorem, and generic
degree two by the classical Galois case. The latter applies
because a quadratic extension in characteristic zero is Galois.
These are the same inputs used in the
[rank-two Keller exclusion](../../../unpublished/common-arithmetic-fibers/LITERATURE_AUDIT.md#classical-rank-two-exclusion).
For \(a=2,j=0\) the only sheet degree is one. For \(a=2,j=1\)
the possibilities are two for \(s=1\) and three for \(s=2\).

## Support after the parabolic shear

In the parabolic case subtract \(g(0)\) and normalize \(g'(0)=1\)
by scaling the scalar and the transverse variables. Absorb the
remaining univariate terms into \(v=z+y^2+h(x)\), with
\(\deg h\le2\), so \(r=x+x^2v\). The \(y\)-degrees are
the same as in \(v=z+y^2\). The
[affine-mate theorem](KELLER_REPEATED_PARALLEL_SCALAR_AFFINE_MATE.md)
excludes degrees zero and one.

Let a polynomial of original degree at most \(D\) have free
\(y\)-degree \(m\) after the shear. Its leading \(y\)-coefficient
has weight at most \(D-\lceil m/2\rceil\) for weights \(1,2\)
on \(x,v\). To check this first take \(h=0\). Fix an \(x\)-exponent
\(n\), and group terms \(x^n y^{d-2b} z^b\) by \(d\).
After factoring \(x^ny^d\), the group is a polynomial \(S(t)\)
in \(t=z/y^2\) with exponents between
\(\max(0,n+d-D)\) and \(\lfloor d/2\rfloor\). Under the shear,
\(t=v/y^2-1\). To contribute to the leading \(y^m\)-coefficient
at \(v^k\) one needs \(d=m+2k\) and \((t+1)^k\mid S(t)\).
For fixed \(n\) the different \(d\)'s give distinct monomials
after expansion and impose this condition separately.
Since \(t\) and \(t+1\) are coprime, the length of the allowed
exponent interval is at least \(k\). This implies
\(n+2k\le D-\lceil m/2\rceil\). The univariate shift
\(v\mapsto v-h(x)\) preserves the weight bound.

If \(Q\) is quadratic in \(y\), the quadratic-mate formula
gives \(q=e(r)x^{2j}\). The rational factor \(e\) is polynomial,
by testing denominators on the component of each \(r\)-fibre
where \(x\ne0\). The
[full exceptional-fibre theorem](KELLER_QUADRATIC_MATE_EXCEPTIONAL_FIBRE.md)
gives \(e(0)=0\). Put \(t=\deg e\). The coefficient of \(v^t\)
in \(q\) is proportional to \(x^{2t+2j}\), and polynomiality
at \(x=0\) gives \(2j+\operatorname{ord}_0 e\ge0\).
The support bound therefore gives
\[
                     4t+2j\le D-1\le6,\qquad 2j\ge-t.
\]
As \(t\ge1\), these leave \(t=1,\ j=0,1\) or \(t=2,\ j=-1\).
The last is excluded by the derivative in (1).
The case \(j=0\) was excluded by the sheet degree.
Thus \(q=rx^2\) after rescaling, and \(D=7\) is necessary.
For this case (1) forces \(s=2\). We will use
\(K=\tau+u/x\) and \(c u^2=\lambda\).

## The remaining exponent in degree seven

The quadratic-mate formula is \(P=\Psi(Q,r)+O\), where
\[
             O=2c\bigl((Q-\tau)Z-rZ^3/3\bigr),\qquad
             Z=x\left(y+\frac{B}{2q}\right).
\tag{2}
\]
The odd part has degree three in \(y\).
If \(P\) has odd \(y\)-degree, it is three and its leading
coefficient is proportional to \(c(r)rx^3\).
Polynomiality forces \(cr\in\mathbb C[r]\), testing on the
nonvertical part of the scalar fibres. The support bound gives
\(4\deg(cr)+3\le5\), so \(c\) has order \(-1\) at \(r=0\).
This is incompatible with \(cu^2=\lambda\).

For even \(y\)-degree \(m=2b\ge4\), let \(\psi_b(r)\) be the
coefficient of \(Q^b\) in \(\Psi\). The leading coefficient
\((\psi_b r^b)x^{2b}\) forces \(\psi_b r^b\in\mathbb C[r]\).
If its \(r\)-degree is \(d\), the bound is
\(4d+2b\le7-b\). Thus \(b=2\), \(d=0\) and
\(\psi_2=\gamma r^{-2}\), with \(\gamma\ne0\).

Write \(H=1+xv\), so \(r=xH\).
The centred formula \(K=C-B^2/(4q)\) gives at most a simple
pole along \(H=0\). If \(K\) were regular there, comparison
with \(K=\tau+u/x\) over \(\mathbb C(x)\) at \(r=0\)
would make both \(\tau\) and \(u\) regular. Then \(B\)
vanishes along \(H=0\), and \(Z\) is a unit-slope affine
coordinate over its local ring. The even and odd parts
must separately be regular. If \(u\) vanishes at zero,
\(c r\) has negative order by \(cu^2=\lambda\), contrary
to the leading odd coefficient. Otherwise \(Q|_{H=0}=
\tau(0)+u(0)/x\) is nonconstant, preventing cancellation of
any maximal coefficient pole of \(\Psi\) after substitution.
In particular the pole of \(\psi_2\) is impossible.
Thus \(K\) has a simple pole, and \(B|_{H=0}\ne0\), so
\(Q|_{H=0}\) is nonconstant.

Put \(b_0=\operatorname{ord}_0 u\). The simple pole implies
\(b_0\ge-1\) and \(\operatorname{ord}_0\tau\ge-1\).
If \(b_0=-1\), \(c\) has order two and \(Z\) has at most a
simple pole along \(H=0\), so both terms of (2) are regular
there. They cannot cancel the poles of \(\Psi\), because
\(Q|_{H=0}\) is nonconstant. Hence \(b_0\ge0\) and
\(\tau\) has a simple pole with nonzero leading coefficient.

For a formal \(T\) define
\[
 Z_*^2=(T-\tau)/r,\qquad
 O_*=\tfrac43\,c\,r^{-2}\,[r(T-\tau)]^{3/2},
\]
choosing the same square root in both formulas. The power
series \(rZ_*\) has nonzero constant term. The pole order of
\(O_*\) is \(N=2b_0+2\), and its coefficient at order
\(-N+k\) has \(T\)-degree at most \(k\).
Put \(U=u/x\). At fixed \(T=Q\), \(Z^2=(T-\tau-U)/r\),
so differentiating (2) with respect to \(U\) gives
\(-cU/(r Z)\). Where \(U/(T-\tau)\) has positive order, the
expansion with compatible signs is
\[
       O-O_*=-\frac{\lambda x^{-2}}{2rZ_*}
                         \left(1+O\!\left(\frac{U}{T-\tau}\right)\right).
\tag{3}
\]
This holds along \(H=0\), where the right side is regular.
Use the sign of \(Z_*\) selected by \(Z\) there. Regularity of
\(P\) at \(H=0\) forces all negative Laurent coefficients of
\(\Psi(T,r)+O_*(T,r)\) to vanish, since the first nonzero one
would remain nonzero on substituting the nonconstant \(Q|_{H=0}\).
The \(T^2r^{-2}\) coefficient therefore gives \(N-2\ge2\),
so \(b_0\ge1\).

Along \(x=0\), \(r=xH\) has order one and \(H\) is a unit.
Thus \(U=u/x\) is regular and (3) applies. If the sign there
were opposite, the pole of \(2O_*\) would have order \(N\ge4\)
and the correction order two. For the matching sign,
\(\Psi+O_*\) is regular, and (3) supplies a nonzero pole of
order two since \(rZ_*\) is a unit. Both cases contradict
polynomiality of \(P\).

## The full sextic exclusion

**Corollary SexticParabolicQuarticScalarExclusion.**
There is no characteristic-zero polynomial Keller triple of
ordinary degree at most six with quartic component
\(R=g(x)+x^2(z+y^2)\).

Use the normalized coordinates \(r=x+x^2v\) above, and write
\(m\ge n\ge3\) for the free degrees and \(p,q\) for their leading
coefficients. The bracket gives \(p^n/q^m\in\mathbb C(r)^\times\).
The [leading-coefficient residue](SEXTIC_KELLER_TRIPLE_LINE_CUBIC_MATE.md#leading-coefficients-and-equal-degrees)
excludes a constant leading coefficient. It also shows that an
invariant leading coefficient must vanish at \(r=0\).
The same residue applies here because it uses only the polynomial
ring in the sheared variables and the equation for its bracket.

The support weights of \(p,q\) are at most
\(6-\lceil m/2\rceil,6-\lceil n/2\rceil\), respectively.
In particular both are at most four. A full fibre factor \(r-a\)
costs weight four, so any coefficient containing one is
proportional to it. Stripping full fibre factors from \(p,q\)
leaves their power ratio constant: a zero or pole at \(r=a\)
would force another full \(r-a\) factor in the corresponding
coefficient by squarefreeness of the fibre.
For equal degrees, if \(p/q\) were nonconstant, one of them
would contain such a factor and both would be invariant.
Then the residue forces both to be proportional to \(r\).
Thus the ratio is constant, and a constant target subtraction
decreases one free degree without changing the ordinary bound.
After reordering, the degrees are unequal and still at least
three, since degrees at most two were excluded above.

Put \(d=\gcd(m,n)\). The stripped leading coefficients are,
up to scalars, \(\ell^{m/d},\ell^{n/d}\) for a common polynomial
\(\ell\). If \(\ell\) is constant, \(p,q\) are invariant and
necessarily multiples of \(r\), so \(m=4,n=3\).
Dividing by the invariant leading coefficients and centering
the cubic gives the monic Laurent setting of the
[coprime pole calculation](SEXTIC_KELLER_TRIPLE_LINE_COPRIME_POLE_LEMMA.md),
with bracket a nonzero invariant instead of \(\kappa/x\).
After dividing the bracket by \(x\) its right side has pole order
at most one. Any pole of the centered cubic would give leading
bracket order at least \((4+3-1)/3=2\), so the same leading
power argument excludes it at both Laurent ends. The
coefficient-drop argument then makes both mates invariant,
a contradiction.

For nonconstant \(\ell\) neither coefficient can have a full
fibre factor by the weight bounds. We have
\[
             (m/d)\operatorname{wt}(\ell)
                           \le 6-\lceil m/2\rceil .
\]
Since \(m/d\ge2\), \(m\le8\). Checking the unequal pairs
\(3\le n<m\le8\) leaves \((4,3),(6,3),(6,4),(8,4)\),
and in each \(\ell\) is linear in \(x\).
For \((6,3),(8,4)\), \(q\) is proportional to \(\ell\)
and has a simple zero in \(\mathbb C(r)[x]\).
The identity \(p^n/q^m\in\mathbb C(r)\) then requires the
degree of every replacement mate to be divisible by \(n\).
Successively subtract invariant multiples of powers of \(Q\)
over \(\mathbb C(r,x)\). Each replacement keeps the bracket
constant, and the process eventually makes a mate independent
of \(y\). Its bracket with \(Q\) then has degree \(n-1\)
or is zero, a contradiction.

For \((4,3)\), the residue for the cubic leading coefficient
requires a primitive of \(dx/(x^2\ell)\) over \(\mathbb C(r)\).
The root of \(\ell\) must therefore be zero. In the centered
coordinate \(Z=xy+[y^2]Q/(3x^2)\) the monic Laurent bracket
is \(\kappa/x\). The coprime pole calculation gives a single
simple pole \(\eta/x\) in its constant term, while the other
centered coefficients of \(Q\) are polynomials in \(r\), by
Laurent comparison in \(\mathbb C[r,x,x^{-1}]\).
At \(y=0\), regular \(Z\) leaves that pole and a pole of \(Z\)
instead makes \(Z^3\) dominate. Both contradict polynomiality.

The pair \((6,4)\) has leading coefficients proportional to
\(\ell^3,\ell^2\) with \(\ell\) linear in \(x\). It is excluded
by the
[centered sextic-quartic calculation](SEXTIC_KELLER_TRIPLE_LINE_SEXTIC_QUARTIC_MATES.md),
which uses the polynomial mates in \(x,v,y\) and their leading coefficients.
