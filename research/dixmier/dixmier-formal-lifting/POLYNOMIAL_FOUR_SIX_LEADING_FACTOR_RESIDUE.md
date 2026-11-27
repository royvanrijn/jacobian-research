# A residue constraint on the common leading factor

**Theorem MonicScalarCommutatorResidue.** Let \(K\) be a
characteristic-zero differential field with derivation \(\delta\).
If a monic operator \(U\in K[E;\delta]\) of positive order has a mate
\(V\in K[E;\delta]\) with \([U,V]=\kappa\ne0\) in the constant field,
then \(1\in\delta K\).

**Corollary PolynomialFourSixLeadingFactorResidue.** Let \(k\) have
characteristic zero and suppose \([P,Q]=\kappa\in k^\times\) for
\(P,Q\in k[x][D]\) of orders four and six, where \([D,x]=1\).
The leading coefficients are \(\alpha h^2,\beta h^3\) for
\(h\in k[x]\) and \(\alpha,\beta\in k^\times\). Over
\(\overline{k}\), the differential \(dx/\sqrt h\) is exact in
\(\overline{k}(x,\sqrt h)\). The polynomial \(h\) has no exact double
root and, if \(\deg h\ge 2\), has a root of multiplicity at least
three.

## The residue identity

Use formal pseudodifferential operators with coefficients on the
left and product rule
\[
 E^i a=\sum_{j\ge0}\binom{i}{j}\delta^j(a)E^{i-j},
 \qquad i\in\mathbb Z.
\]
Each series has powers of \(E\) bounded above. Write
\(\operatorname{res}_E\) for the coefficient of \(E^{-1}\). For a pair
of terms \(aE^i,bE^\ell\), a residue in their commutator can occur
only for \(j=i+\ell+1\ge0\), and is then
\[
 \binom{i}{j}a\delta^j(b)-\binom{\ell}{j}b\delta^j(a).
\]
Modulo the additive subspace \(\delta K\), the Leibniz rule gives
\(a\delta^j(b)=(-1)^j b\delta^j(a)\).
Since \(\ell=j-i-1\), the binomial identity
\(\binom{\ell}{j}=(-1)^j\binom{i}{j}\) shows that this residue lies in
\(\delta K\). In any two series only finitely many pairs of terms
can satisfy \(i+\ell\ge-1\). Thus the residue of a commutator is in
\(\delta K\), even when both series have infinite negative parts.

For \(n=\operatorname{ord}_E U\), let \(L=E+O(E^0)\) be the monic
formal \(n\)-th root of \(U\). It is constructed by comparing
coefficients of \(L^n\) in descending order. Each unknown root
coefficient enters the relevant equation multiplied by \(n\).
The inverse \(L^{-1}\) has leading term \(E^{-1}\), by the same
coefficient recursion for \(LL^{-1}=1\). The identity
\[
 [U,V]=\sum_{i=0}^{n-1}L^i[L,V]L^{n-1-i}=\kappa
\]
gives \([L,V]=(\kappa/n)L^{1-n}\). Indeed this expression solves
the linear equation for \([L,V]\). On a nonzero difference of
two solutions, with leading term \(aE^j\), the left side has leading
term \(naE^{j+n-1}\), a contradiction. Hence
\[
 [L^{n-1},V]=\frac{(n-1)\kappa}{n}L^{-1}.
\]
For \(n>1\) the residue is a nonzero constant. Division by that
constant proves \(1\in\delta K\). For \(n=1\), \(U=E+u\) is
itself an Ore generator with derivation \(\delta\). Express \(V\)
in powers of \(U\). Then the commutator acts coefficientwise by
\(\delta\), and its order-zero coefficient is \(\kappa\).

For algebraic function fields, the retained
[algebraic-time calculation](ALL_COPRIME_ALGEBRAIC_STRING_OBSTRUCTION.md#2-the-commutator-itself-supplies-an-algebraic-time)
also gives this exactness by taking the residue of \([P,QW^{-1}]\).
The finite-root conclusions below use the polynomial leading factor.

## The leading factor at its finite roots

The order-nine coefficient of the original commutator is
\(4p_4q_6'-6q_6p_4'\). Its vanishing makes \(q_6^2/p_4^3\)
constant. Factorization in \(k[x]\) gives the claimed leading
coefficients. Extend constants to \(\overline{k}\) and write \(k\)
for the extended field in the rest of the proof. With
\(\lambda^2=h\) and \(E=\lambda D\), work over \(K=k(x,\lambda)\).
There \(E\) has derivation \(\delta=\lambda\,d/dx\).
After a constant scaling \(P\) is monic in \(E\), with the
commutator still a nonzero constant. The theorem gives a
primitive \(f\in K\) with \(f'=1/\lambda\).

If \(h=H^2\) is square in \(k[x]\), choose \(\lambda=H\) and
\(f\in k(x)\). At an exact double root of \(h\), the function
\(1/H\) has a simple pole with nonzero residue, which no rational
derivative can have. More generally, if \(h\) is nonconstant,
let \(r\) be the number of distinct roots of \(H\) and \(N=\deg H\).
All these roots have multiplicity at least two in \(H\). A root
of multiplicity \(v\) gives a pole of order \(v-1\) for \(f\),
and \(f\) has no other finite poles. At infinity \(f\) is bounded
and \(f-f(\infty)\) has vanishing order \(N-1\), since \(f'=1/H\)
has degree \(-N\) there. The denominator of \(f-f(\infty)\)
has degree \(N-r\), so this vanishing order is at most \(N-r\).
It follows that \(r=1\).

If \(h\) is nonsquare, use the basis \(1,1/\lambda\) over
\(k(x)\) to write \(f=A+B/\lambda\). The two components of
\(f'=1/\lambda\) give \(A'=0\) and
\[
                  B'-\frac{h'}{2h}B=1.                    \tag{1}
\]
At any finite point \(b\), let \(m\) be its multiplicity in
\(h\), taking \(m=0\) away from the roots. If \(B\) has a
pole with leading exponent \(j<0\) in \(x-b\), the leading
coefficient on the left is that of \(B\) multiplied by
\(j-m/2\), at exponent \(j-1\), which is impossible.
Thus \(B\in k[x]\). Every root \(b\) of \(h\) must be a
root of \(B\) to remove the simple pole in (1). A zero of
order two or more in \(B\) gives zero constant term in (1).
A simple zero instead gives
\[
                         B'(b)(1-m/2)=1.                  \tag{2}
\]
This excludes exact double roots of \(h\). At any remaining
root of \(B\) which is not a root of \(h\), equation (1) gives
\(B'(b)=1\). In particular \(B\) is squarefree.

When \(h\) has at least two distinct roots, let \(d=\deg h\)
and \(M=\deg B\ge2\). The coefficient of \(x^{M-1}\) in (1)
is \((M-d/2)\) times the leading coefficient of \(B\).
It must be zero, so \(M=d/2\). In particular no squarefree
\(h\) of degree at least two can occur, since the number of
its roots would exceed \(\deg B\). The residue sum of
\(1/B\) yields the same degree relation. By (2) the residues
at the roots of \(h\) are \(1-m/2\), and those at all other
roots of \(B\) are one. Their sum is \(M-d/2\).

These cases prove the root-multiplicity claim. Combining them
with the [higher-root local gate](POLYNOMIAL_FOUR_SIX_HIGHER_ROOT_LOCAL_GATE.md)
also shows that any candidate of degree at least two has even
degree. In the square branch there is only one distinct root.
In the nonsquare branch there are several distinct roots and a
positive even number of simple roots. Indeed a nonsquare with
a single root would have odd multiplicity, and the local gate
gives even multiplicity for every root of multiplicity at least
three. With several roots the degree is even and the simple
roots are exactly those of odd multiplicity.
