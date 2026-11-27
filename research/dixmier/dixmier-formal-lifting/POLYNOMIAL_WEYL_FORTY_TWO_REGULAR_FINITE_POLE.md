# Regularity of the finite spectral pole in the degree-forty-two Weyl residual

**Theorem PolynomialWeylFortyTwoRegularFinitePole.**
Let \(P,Q\) be a hypothetical degree-forty-two proper pair over
\(\mathbb C\), in the centered Fourier presentation of the
[finite-pole theorem](POLYNOMIAL_WEYL_FORTY_TWO_FINITE_POLE.md).
Every monomial of \(P\) satisfies \(j-i\le18\) for its
\(X^iY^j\) exponents, and every monomial of \(Q\) satisfies
\(j-i\le63\). The finite Euler parameters are \(k=9,\ \ell=15\).
Both spectral equations are regular singular at their only possible
finite pole. For every spectral value they have at least \(18\) and
\(63\) entire solutions, respectively.

## Exclusion of intermediate weight faces

Work with the weight \((\rho,\sigma)=(-a,b)\) for coprime integers
\(0<a<b\). Put \(s=a/b\). The support of \(P\) is in
\([0,12]\times[0,30]\) and includes \((12,30)\). Any maximizer
\((i,j)\) for this weight satisfies
\[
             j-si\ge30-12s,\qquad
             j-i\ge30-12s-(1-s)i\ge18.
\]
The weight of \(P\) is strictly positive. By
[Guccione, Guccione and Valqui, Theorem 4.1 and part (4)](sources/guccione-guccione-valqui-1111.6100v3.pdf),
its weighted leading symbol \(R\) admits a homogeneous
polynomial bracket mate \(H\) of weight \(b-a\), with
\(\{R,H\}=cR\) for a nonzero scalar \(c\). Here the bracket
uses the usual commutative symbols. Changing the bracket convention
changes only the sign of \(c\).

Write \(T=X^bY^a\) and \(R=X^uY^v L(T)\) with \(L(0)\ne0\).
Since \(H\) is polynomial and the integer exponent pairs of weight
\(b-a\) are \((1,1)+n(b,a)\), it has the form \(H=XY G(T)\).
Indeed \(b\ge2\), so negative \(n\) would give a negative
\(X\)-exponent. Computing the bracket and cancelling \(X^uY^v\)
gives, with \(\{X,Y\}=1\),
\[
 (u-v)LG+(b-a)TL'G+(au-bv)TLG'=cL .
\]
Let \(M=\deg L\) and \(N=\deg G\).
If \(N>0\), the coefficient of \(T^{M+N}\) on the left is
the product of the two leading coefficients times
\[
                  (u+bM)-(v+aM)+(au-bv)N .
\]
The first difference is at most \(-18\) by the maximizer
inequality. Also \(au-bv<0\) by positivity of the weight of \(R\).
The coefficient therefore cannot vanish, contradicting the
degree \(M\) of the right side. Thus \(G\) is a nonzero constant.
Comparing the constant and all positive-degree coefficients of
\(L\) now forces \(L\) to be constant because \(b-a>0\).
Hence \(P\) has a unique maximizing monomial at every rational
\(s\in(0,1)\).

## The limiting Euler weight

At \(s=0\) the centered leading differential coefficient gives a
unique maximizer, the corner \((12,30)\). A first change of
maximizer as \(s\) increases would tie two integer support points
at a rational \(s\in(0,1)\). The corner therefore keeps
maximizing all the way to \(s=1\), where the weights attain
\(18\). The finite-pole theorem shows that the weights of \(P,Q\)
at \(s=1\) are \(2k,7k\), respectively, and that \(k=9\) forces
\(\ell=15\). Its regular singularity and entire-solution statements
now apply throughout the degree-forty-two residual. In particular
the root-pair formula applies with six excess pairs and the
finite-to-infinity coefficient map uses only the regular Euler
grades.
