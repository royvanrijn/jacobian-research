# Polynomial Ore intertwiners and subprincipal descent

A nonzero polynomial differential intertwiner between two actual
polynomial Weyl pairs is a bijection on \(k[x]\). If one image has a
centered monomial rectangle with position degree greater than its
differential order, every nonconstant such intertwiner strictly changes
its subprincipal corner by a negative integer. In particular it cannot
connect two pairs with zero subprincipal corner in that rectangle.

This is an actual-global admission theorem for scalar polynomial Ore
intertwiners. It does not classify rational Weyl-division intertwiners,
matrix intertwiners, or nonpolynomial formal gauges, and does not force
the subprincipal corner of an unknown pair to be zero.

Let \(k\) have characteristic zero and
\(A_1(k)=k[x][D]\), \([D,x]=1\). Assume
\[
 P,Q,P',Q'\in A_1(k),\qquad
 [Q,P]=[Q',P']=\kappa\in k^\times,
\]
and a nonzero \(S\in A_1(k)\) satisfies
\[
                          PS=SP',\qquad QS=SQ'.          \tag{1}
\]
No surjectivity of either Weyl algebra endomorphism is assumed.

## Every nonzero polynomial Ore operator is Fredholm on polynomials

Write \(S=\sum a_{ij}x^iD^j\) and define
\[
 m=\max\{i-j:a_{ij}\ne0\},\qquad
 \sigma_S(n)=\sum_{i-j=m}a_{ij}(n)_j,                    \tag{2}
\]
where \((n)_j=n(n-1)\cdots(n-j+1)\). Distinct \(j\) give
linearly independent falling-factorial polynomials, so \(\sigma_S\)
is nonzero. For all sufficiently large integers \(n\),
\[
                     \deg S(x^n)=n+m.                   \tag{3}
\]
Consequently \(S:k[x]\to k[x]\) has finite kernel and cokernel,
and its index is
\[
                       \dim\ker S-\dim\operatorname{coker}S=-m.
                                                               \tag{4}
\]
Indeed choose \(N\) beyond all nonnegative integer roots of
\(\sigma_S\), large enough that \(N+m\geq0\). The kernel is
contained in \(F_N=k[x]_{\leq N}\). Repeated leading-term subtraction
using (3) represents every cokernel class in \(F_{N+m}\), and
\(\operatorname{im}S\cap F_{N+m}=S(F_N)\): an input of degree
greater than \(N\) has output degree greater than \(N+m\).
The finite-dimensional calculation gives
\((N+1)-(N+m+1)=-m\). It also proves both finiteness assertions.

Under (1), the kernel is invariant under \(P',Q'\), and the image
is invariant under \(P,Q\). Thus both finite-dimensional kernel and
cokernel carry the respective scalar commutator \(\kappa\).
Their traces force both dimensions to be zero. Therefore
\[
 S\text{ is a bijection on }k[x],\qquad m=0.              \tag{5}
\]
In particular \(S\) preserves every degree filtration space \(F_N\).
Its restriction to each \(F_N\) is injective and hence bijective,
so the diagonal polynomial in (2) satisfies
\[
                         \sigma_S(n)\ne0\quad(n\geq0).  \tag{6}
\]
These are necessary conditions for an actual pair intertwiner, not
sufficient conditions for polynomial conjugation.

## Exact outer faces and leading-coefficient restriction

Assume further that \(P\) is supported in
\(0\leq i\leq a,0\leq j\leq b\), with both outer faces precisely
\[
                \operatorname{lc}_D(P)=x^a,\qquad
                \operatorname{lc}_x(P)=D^b,
                \qquad a>b>0.                          \tag{7}
\]
Differential order and position degree are additive on products of
nonzero scalar Weyl operators: their associated leading-coefficient
rings are domains. Hence (1) gives
\(\operatorname{ord}_D P'=b\) and \(\deg_xP'=a\).
At the top differential order, cancellation of the nonzero leading
coefficient of \(S\) gives \(\operatorname{lc}_D(P')=x^a\).
At the top position degree, the leading coefficient ring is \(k[D]\);
its commutativity and lack of zero divisors give
\(\operatorname{lc}_x(P')=D^b\). Thus both centered outer faces
of \(P'\) are exactly the same as those of \(P\).

Put \(r=\operatorname{ord}_D S\), and write its leading coefficient
as \(s_r(x)\ne0\). Denote the \(D^{b-1}\) coefficients of
\(P,P'\) by \(p_{b-1},p'_{b-1}\). Comparing the coefficient of
\(D^{b+r-1}\) in (1) gives the **exact** identity
\[
 p'_{b-1}=p_{b-1}
                +b x^a\frac{s_r'}{s_r}-ar x^{a-1}.      \tag{8}
\]
Both sides are polynomial. At a nonzero root \(\alpha\) of
\(s_r\), of multiplicity \(e>0\), the middle term would have
nonzero residue \(b\alpha^a e\). Therefore all roots of \(s_r\)
are at zero, after algebraic constant extension, and
\[
                  s_r=\lambda x^c,
                  \quad \lambda\in k^\times,\ c\geq0.   \tag{9}
\]
Moreover (5) says every monomial of \(S\) has position degree at
most its differential degree, so \(c\leq r\). Equations (8)--(9)
give
\[
                   p'_{b-1}=p_{b-1}+(bc-ar)x^{a-1}.      \tag{10}
\]
No common-centralizer hypothesis enters this argument.

## Consequence for the Weyl subprincipal corner

Let \(A,A'\) be the respective normal PBW coefficients of
\(x^{a-1}D^{b-1}\), and let
\(A_W=A-ab/2\), \(A'_W=A'-ab/2\). Then
\[
                         A'_W-A_W=bc-ar.                 \tag{11}
\]
If \(S\) is nonconstant, (5) forces \(r>0\): an order-zero
polynomial with \(m=0\) is constant. Since \(c\leq r\) and
\(a>b\),
\[
                   ar-bc\ \geq(a-b)r>0                 \tag{12}
\]
is a positive integer. Thus the shift in (11) is a strictly negative
integer, understood in the embedded copy of \(\mathbb Z\subset k\).

Consequently, if \(A'_W=A_W\), including when both are zero,
the only possible \(S\) is a nonzero scalar. In that case (1)
also gives \(P=P'\) and \(Q=Q'\). This rigidity is restricted
to the stated polynomial intertwiner and centered rectangle. It is
not bounded formal-gauge rigidity or a general division-ring theorem.

## A near-identity auxiliary-parameter corollary

Keep the basic relation \([D,x]=1\) fixed, and introduce an independent
central parameter \(\varepsilon\). Suppose the four pair members and
\(S_\varepsilon\) have finite PBW support, with coefficients in
\(k[[\varepsilon]]\), and satisfy (1) with a common
\(0\ne\kappa(\varepsilon)\in k[[\varepsilon]]\). Suppose
\[
 S_\varepsilon\equiv1\pmod\varepsilon,
\]
and the fixed centered outer faces of \(P_\varepsilon\) are
\(\alpha(\varepsilon)x^aD^b\), with \(a>b>0\) and
\(\alpha(\varepsilon)\in k[[\varepsilon]]^\times\). Then
\(S_\varepsilon\) is scalar; hence both pairs are identical.
The faces of \(P'_\varepsilon\) automatically agree by the preceding
leading-coefficient argument.

Indeed work over the generic field \(k((\varepsilon))\) and divide
\(P_\varepsilon,P'_\varepsilon\) by the shared unit \(\alpha\).
If \(A,A'\) are their *unnormalized* mixed PBW coefficients, use
\[
 \mathcal A_W=A/\alpha-ab/2,\qquad
 \mathcal A'_W=A'/\alpha-ab/2.
\]
For a nonconstant \(S_\varepsilon\), (11)--(12) give
\(\mathcal A'_W-\mathcal A_W=bc-ar\), a nonzero integer independent
of \(\varepsilon\). On the other hand reducing (1) modulo
\(\varepsilon\) gives \(P_0=P'_0\). Unit \(\alpha\) implies
\(\mathcal A'_W-\mathcal A_W\in\varepsilon k[[\varepsilon]]\),
a contradiction. Equivalently
\(A'-A=\alpha(bc-ar)\) has valuation zero, contradicting
\(A'-A\equiv0\). The highest derivative coefficient of
\(S_\varepsilon\) may itself start at positive parameter order;
its scalar parameter factor cancels from \(s_r'/s_r\).

This is an auxiliary parameter statement with fixed Weyl parameter.
It does not directly apply when \(\varepsilon\) is the physical
\(\hbar\) in \([D,x]=\hbar\). In that convention (8)'s shift
contains an extra factor \(\hbar\), compatible with reduction
modulo \(\hbar\); equivalently normalizing to \(D/\hbar\)
makes a positive-order corner coefficient nonunit. It also does not
apply to corners first appearing at positive \(\varepsilon\)-order,
intertwiners with unbounded PBW support, or general rational Ore
intertwiners.

For the necessity of the unit-corner argument, a finite **one-operator**
control has
\[
 S=1+\varepsilon xD,\quad
 P=\varepsilon x^3D+(1+2\varepsilon)x^2,\quad
 P'=\varepsilon x^3D+x^2.
\]
It satisfies \(PS=SP'\), \(S\equiv1\), and \(P\equiv P'\),
but its generic highest corner coefficient is the nonunit
\(\alpha=\varepsilon\). This is not a full-pair construction or
a counterexample to a stronger pair statement; it demonstrates the
failed valuation step when the unit hypothesis is removed.

## Finite controls and why each hypothesis matters

[The exact checker](check_polynomial_ore_intertwiner.py) verifies finite
Ore and polynomial-module controls. They supplement the universal
written proof; their finite ranks are not extrapolated to all operators.

* \(S=D,x,xD-2,D^2+xD+1\) give the predicted degree shifts and
  finite-index controls. The first two have nonzero index, and
  \(xD-2\) has a nonzero finite kernel and cokernel.
* \(S=xD+1\) is a bijection on polynomials, with diagonal \(n+1\).
  It nevertheless cannot satisfy \(xS=SP'\) for any polynomial
  \(P'\): order and leading coefficients would force \(P'=x\),
  whereas \(xS-Sx=-x\). Bijectivity alone is insufficient.
* A one-operator Darboux identity uses \(S=xD-2\),
  \(P=x^4D^2\),
  \(P'=x^4D^2-2x^3D+2x^2\). It verifies (1)'s first equation
  and the shift \(-2\), but the finite kernel forbids completing
  both sides to pairs satisfying all of (1).
* \(S=x^2D+x\) commutes with \(P=S^2\), whose subprincipal
  corner is zero. Here \(m=1\); its nonzero cokernel precludes
  an actual full-pair intertwiner. The second intertwining equation
  is therefore indispensable.

[The receipt](polynomial_ore_intertwiner_controls.json) pins this note
and checker. All actual pair claims are precisely the written theorem;
the one-operator examples are explicitly not full Weyl-pair
constructions. No old source or metadata is changed, and no formal
verification, external review or DC1 resolution is claimed.
