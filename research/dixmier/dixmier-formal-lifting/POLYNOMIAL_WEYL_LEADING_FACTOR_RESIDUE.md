# Leading-factor restrictions for polynomial Weyl pairs in arbitrary orders

**Theorem PolynomialWeylLeadingFactorResidue.** Let \(k\) have
characteristic zero, and let \(P,Q\in k[x][D;\partial_x]\) have
positive differential orders \(n,m\) with
\([P,Q]=\kappa\in k^\times\). Put \(g=\gcd(n,m)\). Their leading
coefficients are \(\alpha h^{n/g},\beta h^{m/g}\), with
\(\alpha,\beta\in k^\times\) and \(h\in k[x]\). Over an algebraic
closure of \(k\), for any choice of \(\lambda^g=h\),
\(dx/\lambda\) is exact in \(k(x,\lambda)\), where \(k\) now
denotes the extended constant field. Equivalently there is a
polynomial \(B\) in that extended \(k[x]\) satisfying
\[
                        B'-\frac{h'}{g h}B=1.                       \tag{1}
\]
Every root of \(h\) is a simple zero of \(B\), and no root of \(h\)
has multiplicity \(g\). If \(h\) has at least two distinct roots,
\[
    \deg B=\frac{\deg h}{g},\qquad
    g\mid\deg h,\qquad \deg h\ge2g,
\]
and at least one root of \(h\) has multiplicity greater than \(g\).

Write \(p,q\) for the leading coefficients of \(P,Q\). The
order-\((n+m-1)\) commutator term gives \(n p q'-m q p'=0\).
Therefore \(q^{n/g}/p^{m/g}\) is constant, and factorization in
\(k[x]\) gives the displayed leading factor. Extend constants
to an algebraic closure, still denoted by \(k\). In
\(K=k(x,\lambda)\), use \(E=\lambda D\) as the derivative
generator. Its derivation is \(\lambda\,d/dx\). The leading
coefficient of \(P\) in \(E\) is the nonzero constant \(\alpha\).
After scaling, the
[monic residue theorem](POLYNOMIAL_FOUR_SIX_LEADING_FACTOR_RESIDUE.md#the-residue-identity)
produces \(f\in K\) with \(f'=1/\lambda\).

For \(d=[K:k(x)]\) set \(B=d^{-1}\operatorname{Tr}_{K/k(x)}(\lambda f)\).
Differentiation commutes with the trace, since every embedding of
the finite separable extension into an algebraic closure commutes
with the unique extension of \(d/dx\). Since
\((\lambda f)'=1+(h'/gh)\lambda f\), this gives (1) in \(k(x)\).
Conversely \((B/\lambda)'=1/\lambda\) for any solution of (1).

At a finite point \(x_0\), let \(\mu\ge0\) be the multiplicity of
\(h\). If \(B\) had a pole there, its leading exponent \(j<0\)
in \(x-x_0\) would give a nonzero coefficient proportional to
\(j-\mu/g\) at exponent \(j-1\) on the left of (1). This contradicts
(1). Hence \(B\) is a polynomial. At any root of \(h\), \(B\) has
to vanish to remove the simple pole. Its zero is simple, since a
higher zero would give zero constant term on the left of (1),
and
\[
                        B'(x_0)(1-\mu/g)=1
\]
there. In particular \(\mu\ne g\).
Every extra root of \(B\) is simple as well, with slope one by (1).

If \(h\) has at least two distinct roots, \(M=\deg B\ge2\).
In (1) the coefficient of \(x^{M-1}\) is proportional to
\(M-\deg h/g\) and must vanish. This proves the degree identity
and lower bound. If all the root multiplicities of \(h\) were less
than \(g\), then \(\deg h/g\) would be strictly less than the
number of distinct roots, although \(B\) vanishes at every root.
The multiplicity restriction follows.

For instance in the remaining polynomial order-\((6,9)\) case,
exact triple roots of \(h\) are excluded. A common leading factor
with at least two distinct roots has degree divisible by three
and at least six, and must have a root of multiplicity at least
four. These are necessary conditions on the leading coefficients.
