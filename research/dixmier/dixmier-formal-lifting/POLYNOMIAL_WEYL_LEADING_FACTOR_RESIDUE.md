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

## Perfect-power leading factors

**Corollary PolynomialWeylPerfectPowerLeadingFactor.** Under the
hypotheses of PolynomialWeylLeadingFactorResidue, if \(h\) is a
\(g\)-th power over the algebraic closure, it is constant or has
one geometric root of multiplicity \(gj\) for some \(j\ge2\).

Work over the algebraic closure and write \(h=f^g\). Equation (1)
gives \((B/f)'=1/f\). The
[rational-primitive obstruction](POLYNOMIAL_WEYL_CONDITIONAL_THROUGH_FIFTY_ONE_FROM_F7.md#an-integer-multiplier-obstruction)
excludes two or more distinct roots of \(f\), and a simple root
is excluded by the pole of order one in \(1/f\). Thus a nonconstant
\(f\) has exactly one root, of multiplicity \(j\ge2\).

## A rational map from the residue equation

**Corollary PolynomialWeylLeadingFactorRamification.** Suppose \(h\)
has \(r\ge2\) geometric roots of multiplicities \(m_i\) and admits
the polynomial \(B\) of (1). Put \(M=\deg h/g\), and let \(\Delta\)
be the gcd of the integers \(|g-m_i|\), including \(g\) as an
additional entry if \(M>r\). Over the algebraic closure there is
a rational map \(S\) with \(S^\Delta=B^g/h\) and degree
\[
                  N=\frac{\sum_{m_i>g}(m_i-g)}{\Delta}.
\]
Its finite zeros have orders \((g-m_i)/\Delta\) at \(m_i<g\) and
\(g/\Delta\) at the \(M-r\) extra roots of \(B\). Its poles have
orders \((m_i-g)/\Delta\) at \(m_i>g\). Its only possible branch
values are \(0,\infty,S(\infty)\), with \(S(\infty)\ne0,\infty\).
The fibre over \(S(\infty)\) contains infinity with ramification
index \(M-1\), and every other point of that fibre is unramified.
In particular \(N\ge M-1\).

The residue theorem gives \(\deg B=M\), and every zero of \(B\) is
simple. It follows that the finite divisor of \(B^g/h\) has the orders
used to define \(\Delta\). At infinity it has order zero since
\(\deg B^g=\deg h\). Factorization over the algebraic closure
gives a rational \(\Delta\)-th root \(S\), up to a root of unity.
Logarithmic differentiation of (1) gives
\[
                         \frac{S'}{S}=\frac{g}{\Delta B}.
\]
At finite points off the zeros and poles of \(S\), the derivative
is nonzero because the zeros of \(B\) are exactly the points in
its divisor. In the coordinate \(z=1/x\) near infinity,
\(dS/dz\) has order \(M-2\), since \(S\) has finite nonzero value
there and \(B\) has degree \(M\). In characteristic zero this
makes \(S-S(\infty)\) vanish to order \(M-1\). The degree follows
by summing the pole orders, and the index at infinity cannot
exceed that degree.

For instance, if \(1\le\ell<g\), a pattern with \(t\) roots of
multiplicity \(g-\ell\) and \(t\) of multiplicity \(g+\ell\) gives
\(M=r=2t,\ \Delta=\ell,\ N=t\). Thus \(t\ge2\) is impossible
for any choice of the roots. At \(t=1\), the residue equation
does have a solution. For \(h\) proportional to
\((x-a)^{g-\ell}(x-b)^{g+\ell}\) with \(a\ne b\), it is
\[
                      B=\frac{g(x-a)(x-b)}{\ell(a-b)} .
\]
This is a statement about the leading-factor equation.

## Real roots and parity

**Corollary PolynomialWeylLeadingFactorRealParity.** Suppose
\(h\in\mathbb R[x]\) has at least two geometric roots and a
polynomial solution of (1). If every real root of \(h\) has
multiplicity less than \(g\), then \(h\) has at most one real
root and \(M=\deg h/g\) is odd. For odd \(g\), there is exactly
one real root and its multiplicity is odd. For even \(g\), any
sole real root has even multiplicity.

Take \(B\) real by averaging a polynomial solution with its complex
conjugate. The degree and simple-zero conclusions of the residue
theorem apply to this solution. Off the real roots define
\(V(x)=B(x)/|h(x)|^{1/g}\). Each interval has
\[
                           V'(x)=|h(x)|^{-1/g}>0 .
\]
At a real root of multiplicity \(m<g\), the simple zero of \(B\)
makes \(V\) tend to zero from both sides because \(1-m/g>0\).
Thus \(V\) extends continuously and strictly increasingly across
the real line. It can vanish at at most one real root of \(h\).
If \(A,b\) are the leading coefficients of \(h,B\), the limits at
\(+\infty\) and \(-\infty\) are \(b/|A|^{1/g}\) and
\((-1)^M b/|A|^{1/g}\), respectively. Since the function is
strictly increasing, the limits differ and \(M\) is odd.
Nonreal roots contribute even total multiplicity, which gives the
last two conclusions from \(\deg h=gM\).
