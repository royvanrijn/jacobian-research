# Two coefficient slices of quartic-leading winding-one cubics

Let \(K\) have characteristic zero. On \(t^2+2xy=1\), put
\(U=1+x\), \(v=Ut\), and
\[
 f=\frac{U}{x}P(x,v),\qquad P(x,v)=v^3+d(x)v^2+b(x)v+c(x).
\]
The sphere functional \(\Lambda\) takes the constant term in \(x\)
of normalized integration in \(t\) over \([-1,1]\). Regularity of
\(f\) is equivalent to \(P(0,\pm1)=0\), so for \(d_0=d(0)\) the
polynomial at \(x=0\) is \(P_0(v)=(v^2-1)(v+d_0)\).
This is the normalized \(E=U^4\) branch of the
[quadruple-root reduction](WINDING_ONE_CUBIC_LOGARITHMIC_BOUNDARY.md#a-quadruple-root).

**Proposition 1.** Suppose \(\deg d\le3\), \(b=-1\), and either
\(c=-d_0\) or \(c=-d(x)\). Then all positive moments of \(f\) vanish if
and only if \(d\) is constant. In particular no phase of least sphere
degree seven in either slice is eventually pure. The first slice covers
every regular phase with constant \(b,c\).

**Proposition 2.** Suppose \(P(x,v)=P_0(v)+x^3Q(v)\) for a polynomial
\(Q\) of degree at most two. Then all positive moments of \(f\) vanish
if and only if \(Q=0\). In particular no phase of least sphere degree
seven in this slice is eventually pure.

## The fixed interval formula

Changing the interval variable to \(v\) gives, for \(m\ge1\),
\[
 \Lambda(f^m)=[x^m]U^{m-1}\frac12\int_{-U}^{U}P(x,v)^m\,dv
             =[x^m]U^{m-1}\frac12\int_{-1}^{1}P(x,v)^m\,dv.       \tag{1}
\]
For the second equality substitute \(v=1+xz\) on the first endpoint
strip and \(v=-1-xz\) on the other. The polynomial on each strip is
divisible by \(x\) because \(P(0,\pm1)=0\). The Jacobian supplies one
additional factor of \(x\), so the strip contributions have order at
least \(m+1\). This is a polynomial identity over \(K\).

The moment computations below use vanishing of the odd monomials and
the beta integral
\[
 \int_{-1}^1 v^{2j}(1-v^2)^a\,dv
 =\frac{2^{a+1}a!}{(2j+1)(2j+3)\cdots(2j+2a+1)},\qquad j,a\ge0.  \tag{2}
\]
These rational identities follow by integrating the binomial expansion,
or by substituting \(v^2\) and applying the beta recurrence.

## Constant or matched constant terms

For Proposition 1 write \(P=P_0+h(x)Q_\rho(v)\), where
\(h=h_1x+h_2x^2+h_3x^3\) and \(Q_\rho(v)=v^2-\rho\). Here \(\rho=0\)
gives \(c=-d_0\) and \(\rho=1\) gives \(c=-d(x)\). The normalization
and (2) give the following values, with every integral over \([-1,1]\).

| Integral | \(\rho=0\) | \(\rho=1\) |
|---|---|---|
| \(\tfrac12\int Q_\rho\,dv\) | \(1/3\) | \(-2/3\) |
| \(\int P_0Q_\rho\,dv\) | \(-4d_0/15\) | \(16d_0/15\) |
| \(\int P_0^2Q_\rho\,dv\) | \(\tfrac{16}{105}(d_0^2+1/3)\) | \(-\tfrac{32}{35}(d_0^2+1/9)\) |
| \(\int P_0^3Q_\rho\,dv\) | \(-\tfrac{32}{315}d_0(d_0^2+9/11)\) | \(\tfrac{256}{315}d_0(d_0^2+3/11)\) |

The first moment is \(h_1\) times the first row, so purity forces
\(h_1=0\). The second moment is then \(h_2\) times the second row.
If \(d_0\ne0\), purity forces \(h_2=0\).

If \(d_0=0\), the polynomial \(P_0\) is odd. In the fourth moment the
linear perturbation integrates to zero because \(P_0^3Q_\rho\) is odd.
The quadratic perturbation contributes
\[
 3h_2^2\int_{-1}^{1}P_0(v)^2Q_\rho(v)^2\,dv,
\]
and higher perturbations have \(x\)-order greater than four. At \(d_0=0\)
the integral is a positive rational number, as its integrand is the
square of a nonzero real polynomial. Thus \(h_2=0\) in this case too.

Now \(h=h_3x^3\) and the third and fourth moments are \(3h_3/2\) and
\(6h_3\) times the entries in the last two rows, respectively. In each
column the last two polynomials in \(d_0\) have no common root. Thus
purity requires \(h_3=0\). Conversely, \(P=P_0\) has zero positive
moments in (1) because \(U^{m-1}\) has degree \(m-1\).

## Perturbations at order three in the base variable

In Proposition 2 the third through fifth moments have only the linear
contributions of \(x^3Q\). Their vanishing gives
\(\int_{-1}^1P_0^jQ\,dv=0\) for \(j=2,3,4\), since each coefficient
\([x^{m-3}]U^{m-1}\) is a nonzero binomial coefficient. Put \(u=d_0^2\).
In the coefficient order \((v^2,v,1)\), the three rows of this linear
system, after scaling individually by nonzero rational numbers, are
\[
 \begin{pmatrix}
 3u+1 &6d_0 &3(7u+1)\\
 d_0(11u+9)&3(11u+1)&33d_0(3u+1)\\
 13u^2+18u+1&4d_0(13u+3)&143u^2+78u+3
 \end{pmatrix}.                                                  \tag{3}
\]
The first coordinates of the cross products of rows \((1,2)\) and
\((1,3)\) are \(-9(11u^2-4u+1)\) and
\(-6d_0(39u^2-10u+3)\). If \(d_0=0\) the first is nonzero.
Otherwise a common zero would have \(46u-6=0\) by subtracting the
suitably scaled quadratics, and \(11u^2-4u+1\) is nonzero at
\(u=3/23\). Hence the matrix has rank at least two. Since
\(\int_{-1}^1P_0^j P_0'\,dv=0\), its kernel is exactly the line spanned
by \(P_0'\). Write \(Q=\kappa P_0'\).

If \(\kappa\ne0\), the linear terms in the sixth and seventh moments
again integrate to zero. The quadratic terms in (1) have nonzero
binomial factors and are proportional to
\[
 \begin{aligned}
 \int_{-1}^{1}P_0^4(P_0')^2\,dv&=\frac{1024}{765765}J_6(u),\\
 \int_{-1}^{1}P_0^5(P_0')^2\,dv&=-\frac{2048}{765765}d_0J_7(u),\\
 J_6(u)&=221u^3+425u^2+119u+3,\\
 J_7(u)&=85u^3+187u^2+75u+5.
 \end{aligned}
\]
These evaluations and (3) follow from (2). Vanishing forces
\(d_0\ne0\) since \(J_6(0)=3\), and then requires a common root of
\(J_6,J_7\). The exact identity
\[
 -(4505u^2+8976u+2343)J_6(u)
 +17(689u^2+1182u+161)J_7(u)=6656
\]
excludes that possibility. Thus \(Q=0\), and its converse follows
from (1).

## Eventual purity in least degree seven

Both slices have \(E=U^4\), so \(E(0)=1\), actual height three and
negative winding, and \(x=-1\) is the unique leading root of
multiplicity four. Their numerator coefficients satisfy
\(U^3\mid D\), \(U^2\mid B\), and \(U\mid C\). The
[degree-seven boundary corollary](WINDING_ONE_CUBIC_LOGARITHMIC_BOUNDARY.md#consequences-in-degree-seven)
therefore makes eventual purity imply purity at every order in
exact least degree seven. In either proposition the pure members
have \(P=P_0\) and least sphere degree six by the
[winding-one degree formula](MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md).
Hence the propositions exclude eventual purity in exact degree seven
on the two displayed coefficient slices.
