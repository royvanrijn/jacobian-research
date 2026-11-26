# The Davenport–Zannier bound is strict in the marked quadratic model

Let \(M,\Sigma\in\mathbb C[u]\), with \(M\ne0\),
\(m=\operatorname{ord}_0M\ge2\), and \(\Sigma(0)\ne0\). Consider
the marked surface

\[
V=\{w^2=u\Sigma(u)^2+M(u)z\}
\]

and regular functions \(P,Q\in\mathcal O(V)\). Use the same symbols
for their expressions in \(K[w]\), where \(K=\mathbb C(u)\) and
\(z=(w^2-u\Sigma(u)^2)/M(u)\). Suppose

\[
P_uQ_w-P_wQ_u=\frac1M.
\tag{1}
\]

Write \(r=\deg_wP=dr_0\), \(s=\deg_wQ=ds_0\),
where \(d=\gcd(r,s)\) and \(r_0,s_0\ge2\). Put \(n=rs/d\).
There is a unique \(\lambda\in\mathbb C^*\) cancelling the leading
terms of \(E=Q^{r_0}-\lambda P^{s_0}\). Then

\[
\deg_w E\ge n-r-s+2.
\tag{2}
\]

In particular, \((r,s)=(8,12)\) gives \(\deg_w(Q^2-\lambda P^3)\ge6\).
These are degrees in \(w\) over \(K\). If \(m=2\) and both restrictions
of \(P,Q\) to the central line \(u=w=0\) are nonzero polynomials in
\(z\), then the stronger bound

\[
\deg_w E\ge n-r-s+3
\tag{2a}
\]

holds. In particular \((8,12)\) gives degree at least seven under
these extra hypotheses. The reduced frame in the
[common leading divisor reduction](NORMAL_FORM_COVERAGE_AND_COMMON_LEADING_DIVISORS.md)
has \(m=2\), and both central restrictions have degree at least two.
Coverage of plane Keller maps by such models remains open.

## 1. The two derivative identities

Let \(\ell_P,\ell_Q\) be the leading coefficients. The coefficient of
\(w^{r+s-1}\) in (1) gives \(s_0\ell_P'/\ell_P=r_0\ell_Q'/\ell_Q\),
where prime denotes \(d/du\) on \(K\). Consequently
\(\ell_Q^{r_0}/\ell_P^{s_0}=\lambda\in\mathbb C^*\).
If integers \(i,j\) satisfy \(ir_0+js_0=1\), then
\(H=\ell_P^i\ell_Q^j\) gives
\(\ell_P=\alpha H^{r_0}\), \(\ell_Q=\beta H^{s_0}\) with
\(\alpha,\beta\in\mathbb C^*\).

Set

\[
W=r_0PQ_w-s_0P_wQ,\qquad D=r_0PQ_u-s_0P_uQ.
\]

For either \(i=u\) or \(i=w\) we have

\[
P E_i-s_0P_i E=Q^{r_0-1}(r_0P Q_i-s_0P_iQ).
\tag{3}
\]

If \(E=0\), both \(D\) and \(W\) vanish by (3), which contradicts
\(W P_u-D P_w=r_0P/M\) from (1). Thus \(e=\deg_wE<n\) is defined.
In (3) for \(w\), the coefficient of \(w^{r+e-1}\) on the left is
\((e-n)\ell_P\ell_E\ne0\). For \(u\), the degree on the left is
at most \(r+e\). Hence, with \(\Delta=e-(n-r-s+1)\),

\[
\deg_w W=\Delta,\qquad \deg_w D\le\Delta+1.
\tag{4}
\]

All roots of \(P\) and \(Q\) over an algebraic closure of \(K\)
are simple, and their root sets are disjoint. Indeed, along an
algebraic root \(\theta\) of \(P\), differentiation of
\(P(u,\theta(u))=0\) gives \(P_u=-\theta'P_w\) there.
A repeated root would make \(P_u=P_w=0\). At a shared root the
gradient of \(Q\) would also satisfy \(Q_u=-\theta'Q_w\).
Both alternatives contradict (1). A root of \(E\) has \(P,Q\ne0\),
so \(dE=r_0Q^{r_0-1}dQ-\lambda s_0P^{s_0-1}dP\) is nonzero there
by (1). The same differentiation argument makes every root of
\(E\) simple.

## 2. Equality fixes the root configuration

Suppose \(\Delta=0\). Equation (4) makes \(W\in K^*\) and
\(v=D/W=Aw+B\) with \(A,B\in K\). Then

\[
r_0P(Q_u-vQ_w)=s_0Q(P_u-vP_w).
\]

Work in a finite splitting field of \(PQ\), with the unique extension
of \(d/du\). At each root \(\theta\) of either polynomial the last identity
and implicit differentiation give \(\theta'=-v(\theta)\). Choose two
distinct roots \(\theta_1,\theta_2\) of \(P\). Each ratio
\((\theta-\theta_1)/(\theta_2-\theta_1)\) has zero derivative since
each difference satisfies \(y'=-Ay\). The constants of any finite
extension of \(K\) are \(\mathbb C\). Indeed, if an algebraic
element has zero derivative, differentiation of its monic minimal
polynomial over \(K\) gives a relation of smaller degree, so all
coefficients of the minimal polynomial are constants.

Put \(\gamma=1/(\theta_2-\theta_1)\), \(\delta=-\gamma\theta_1\),
and \(T=\gamma w+\delta\). The leading coefficient identity in
Section 1 now gives

\[
P=h^{r_0}p(T),\qquad Q=h^{s_0}q(T),\qquad h=H\gamma^{-d},
\tag{5}
\]

where \(p,q\in\mathbb C[T]\) have degrees \(r,s\), leading
coefficients \(\alpha,\beta\) and disjoint root sets. Substitution in
\(W\) shows that \(r_0pq'-s_0p'q=c\in\mathbb C^*\), with primes
on \(p,q\) denoting differentiation in \(T\). In the determinant
(1), the terms involving \(T_u\) cancel, giving

\[
c\gamma h^{r_0+s_0-1}h'=\frac1M.
\tag{6}
\]

The classical interpretation is the Davenport–Zannier equality
case. The powers \(Q^{r_0}\) and \(\lambda P^{s_0}\) have degree
\(n\), disjoint roots with multiplicities \(r_0\) and \(s_0\)
respectively, and difference of degree \(n+1-(r+s)\).
Pakovich and Zvonkin's *Minimum degree of the difference of two
polynomials over Q, and weighted plane trees*, Corollary 2.3,
Propositions 2.4 and 2.11–2.12, and Theorem 2.24 give the bound,
the three branch value criterion and the weighted tree description.
There are finitely many trees at fixed passport. The nonsimple
point at \(w=\infty\) above the leading cancellation value marks
infinity, reducing equivalence of source coordinates to affine
changes. The differential argument above obtains (5) over the
finite splitting field directly. [Pakovich–Zvonkin, 2014](https://arxiv.org/abs/1306.4141)

## 3. The two marked arcs contradict equality

Extend the splitting field to include \(\sqrt u\) and choose a place
over \(u=0\). Use its valuation \(\nu\) normalized by \(\nu(u)=1\).
For a local parameter \(t\) we can write \(u=t^N\) after absorbing
the unit into \(t\) in the completed local field. The residue field
is \(\mathbb C\). Differentiating the Laurent series of an element
therefore decreases its order by one in the \(\nu\) normalization
if the order is nonzero. For a unit with nonzero derivative, the
derivative has order strictly greater than \(-1\).

The points \(w_\pm=\pm\sqrt u\,\Sigma(u)\), \(z=0\), give two
regular arcs through the marked locus. Both evaluations of each of
\(P,Q\) have nonnegative order. Their affine coordinates
\(T_\pm=\gamma w_\pm+\delta\) satisfy

\[
\nu(T_+-T_-)=\nu(\gamma)+\tfrac12.
\tag{7}
\]

Put \(k=\nu(h)\). We use all three possible signs of \(k\).

- If \(k>0\), (6) gives
  \(\nu(\gamma)+\tfrac12=\tfrac32-m-(r_0+s_0)k<0\).
  By (7), at least one arc has \(\nu(T_\pm)\) at most this value.
  There the leading term of \(p(T_\pm)\) has uniquely least order.
  As \(r=dr_0\), we obtain
  \[
  \nu(P(u,w_\pm))=r_0(k+d\nu(T_\pm))
    \le r_0\bigl(d(\tfrac32-m)+(1-d(r_0+s_0))k\bigr)<0.
  \]
- If \(k=0\), (6) gives \(h'\ne0\) and
  \(\nu(\gamma)=-m-\nu(h')<1-m\).
  Again (7) gives a pole of \(T_\pm\) on at least one arc, and
  \(\nu(P(u,w_\pm))=r\nu(T_\pm)<0\) there.
- If \(k<0\), consider either arc. If \(T_\pm\) has a pole, so
  does \(h^{r_0}p(T_\pm)\). If \(T_\pm\) has nonnegative order,
  at least one of \(p,q\) is nonzero at its residue. The corresponding
  value in (5) then has negative order.

This proves (2) using the actual two marked arcs and the full
Jacobian determinant.

## 4. Excess one at a double zero of \(M\)

Assume both central restrictions are nonzero. Put \(a=r_0\),
\(b=s_0\) and suppose \(\Delta=1\). Write

\[
W=A(u)(w-\theta(u)),\qquad A\in K^*,\quad \theta\in K,
\qquad L=r+s-2.
\]

Let \(f=Q^a/(\lambda P^b)\). At \(w=\infty\), its Laurent expansion
begins with \(1+f_0 w^{-L}\), where \(f_0\ne0\) because
\(\deg_w E=n-L\). The identity \(f_w=Q^{a-1}W/(\lambda P^{b+1})\)
gives \(f_0\in\mathbb C^* A H^{-(a+b)}\). The Jacobian gives

\[
P_u-\frac{f_u}{f_w}P_w=\frac{aP}{MW}.
\tag{8}
\]

Since \(f_u/f_w=-(f_0'/Lf_0)w+O(1)\) at infinity, comparison in
degree \(r\) in (8) gives \(aH'/H+(r/L)f_0'/f_0=0\). It follows
that \(A'/A=(2/d)H'/H\).

To determine the order of \(A\), write the regular representatives
of \(P,Q\) as \(P_0(u,z)+wP_1(u,z)\) and \(Q_0(u,z)+wQ_1(u,z)\),
and write \(p_0=P_0(0,z)\), \(q_0=Q_0(0,z)\). Holding \(u\)
fixed on \(M\ne0\) gives \(M\partial_w=M\partial_w|_z+2w\partial_z\)
on those representatives. In the completion at the generic point
of the central line, \(w\) is the local parameter and
\(u=w^2/\Sigma(0)^2+O(w^4)\). The intrinsic Jacobian on the
central line is

\[
P_1(0,z)q_0'(z)-p_0'(z)Q_1(0,z)=\frac1{\Sigma(0)^2}.
\tag{9}
\]

Indeed, differentiation of the surface equation gives
\(du\wedge dw/M=dw\wedge dz/((u\Sigma^2)'+M'z)\) near the
line. The factor in parentheses is a unit there and \(du\)
restricts to zero at \(w=0\) in local coordinates \(w,z\). The
normal first jet of \(MW\) is therefore

\[
MW=2w(a p_0q_0'-b p_0'q_0)+O(w^2).
\tag{10}
\]

The coefficient of \(w\) is nonzero. If the bracket vanished,
\(q_0^a/p_0^b\) would be constant, and the coprime integers \(a,b\)
would give common degree powers of a polynomial with exponents
\(a,b\). Either both \(p_0,q_0\) would be constant, or their
derivatives would both vanish at a root of that polynomial. Each
possibility contradicts (9). Now the two summands of
\(MW=MAw-MA\theta\) have odd and even orders in \(w\)
respectively, with distinct orders whenever both are nonzero.
Equation (10) forces

\[
\operatorname{ord}_0 A=-m,\qquad \operatorname{ord}_0\theta\ge1
\tag{11}
\]

with infinite order allowed when \(\theta=0\). In particular the
bracket in (10) is a nonzero constant. The leading coefficient
relation gives \(\operatorname{ord}_0 H=-dm/2\). An odd ruling
degree would make the coefficient of \(z^{\lfloor r/2\rfloor}\)
in \(P_1\), or its analogue for \(Q\), have order \(-m/2\),
which is impossible for a regular representative. Both ruling
degrees are even, and the leading \(z\) coefficients of \(P_0\)
and \(Q_0\) are units at zero. Thus the central degrees are
\((r/2,s/2)\).

Now adjoin \(\sqrt A\) and use the affine coordinate
\(x=\sqrt A(w-\theta)\). In this coordinate the leading coefficients
of \(P,Q\) are constant by \(A'/A=(2/d)H'/H\). The transformed
polynomial \(aP Q_x-b P_x Q\) is \(x\), and
\(P_uQ_x-P_xQ_u=j=1/(M\sqrt A)\), with \(u\)-derivatives now
holding \(x\) fixed. Write the top two terms of \(P\) as
\(\ell(x^r+p x^{r-1})\) and of \(Q\) as
\(k(x^s+q x^{s-1})\), with \(\ell,k\) nonzero constants.
The order \(L>1\) at infinity of \(f-1\) gives \(a q=b p\).
Consequently, for a nonzero constant \(C\),

\[
f_x=Cx^{-L-1}\left(1-\frac{r+s}{r}\frac{p}{x}+O(x^{-2})\right).
\]

Integrating in \(x\) with constant term \(f(\infty)=1\) and
differentiating in \(u\) gives
\(f_u/f_x=(r+s)p'/(r(L+1))+O(x^{-1})\). In the transformed
equation (8), the coefficient of \(x^{r-1}\) is now

\[
\left(1-\frac{r+s}{L+1}\right)p'=a j,
\qquad p'=-\frac{a(r+s-1)}{M\sqrt A}.
\tag{12}
\]

If \(m=2\), (11) makes \(du/(M\sqrt A)\) have a nonzero residue
at every place over \(u=0\): one can choose a local branch of
\(\sqrt A\) in \(\mathbb C((u))\) of order \(-1\). Equation (12)
would make this differential a constant multiple of the exact
differential \(dp\), whose Laurent expansion has zero residue.
This excludes \(\Delta=1\) and proves (2a).

## 5. Branch data and monodromy for excess one

Define \(f=Q^{r_0}/(\lambda P^{s_0})\). Over an algebraic closure of
\(K\) it has degree \(n\). The fibres of \(0,\infty\) have uniform
multiplicities \(r_0,s_0\) respectively, and the fibre of \(1\)
has \(e\) simple finite points and the point \(w=\infty\) of
index \(n-e\). Away from these fibres the finite ramification is
given by the zeros of \(W\), counted with multiplicity. Indeed
\(f_w=fW/(PQ)\). Such zeros are disjoint from the finite points in
the three fibres by the simplicity and disjointness in Section 1.
Thus a finite critical point contributes its ramification excess
to \(\Delta\). Distinct critical points can have the same value
when \(\Delta>1\).

If \(\Delta=1\), write \(\theta\in K\) for the unique root of the
linear polynomial \(W\). Its branch value \(f(u,\theta(u))\) has
nonzero derivative, since

\[
  \frac{d}{du}f(u,\theta(u))=\frac{fD}{PQ}(u,\theta(u)),\qquad
  -D P_w(u,\theta(u))=\frac{r_0P(u,\theta(u))}{M}.
\]

For \((r,s)=(8,12)\), the degree is \(24\). The partitions over
\(0,\infty,1\) and the moving value are respectively
\(2^{12},3^8,(18,1^6),(2,1^{22})\). Roberts develops the four-value
Hurwitz framework and uses ramification partitions to test for
intermediate covers in Sections 2 and 4 of a different family.
[Roberts, 2018](https://www.numdam.org/articles/10.5802/pmb.22/)

For the displayed passport, (1) forces geometric monodromy
\(S_{24}\) even before the double-zero condition is used. An
imprimitive block system would give a decomposition
\(f=g\circ h\) with inner degree \(k\) and outer degree \(24/k\).
This follows from the correspondence between block stabilizers
and intermediate fields, and Lüroth's theorem. Consider the
branch permutation over \(1\). Its 18-cycle cannot remain within
a single block, since a proper block has size at most 12. Its
moved blocks consist entirely of points in the 18-cycle. Hence
\(k\) divides both 18 and 24, giving \(k\in\{2,3,6\}\).
The index of \(g\) at \(h(\infty)\) is \(18/k\), so \(h\) has
index \(k\) at infinity and can be chosen polynomial. The
transposition over the fourth value fixes all blocks, so it
uses one of the \(k-1\) finite ramification indices counted
with their excess for \(h\).

Every index-one point of \(g\) over zero requires a full fibre
of index-two points for \(h\), using \(k/2\) finite excess.
Every index-one point of \(g\) over infinity similarly requires
\(2k/3\) finite excess. These are finite fibres, since
\(g(h(\infty))=1\). All other outer indices over zero are two,
and all other outer indices over infinity are three.

- For \(k=2\) the fourth value uses the entire finite excess.
  Thus \(g\) is a quotient of a square by a cube, up to a
  constant. The coprimality of \(P,Q\) makes each a scalar
  multiple of a polynomial in \(h\).
- For \(k=3\), the outer degree is eight. At least two outer
  index-one points over infinity are needed to supplement the
  index-three points, costing four units of finite excess.
  The available excess is two.
- For \(k=6\), the outer degree is four. One index-one point
  \(c\) over infinity is required, using four units of finite
  excess. The fourth value uses the fifth, leaving no index-one
  point of \(g\) over zero. The other point \(d_1\) over infinity
  has index three. Therefore \(g\) is a constant multiple of
  \(q_1(y)^2/((y-c)(y-d_1)^3)\) with \(q_1\) quadratic, and
  \(h-c\) is a constant multiple of the cube of a quadratic
  \(h_1\). It follows from coprimality that \(P,Q\) are scalar
  multiples of \(h_1(h-d_1)\) and \(q_1(h)\), respectively.
  Both are polynomials in \(h_1\).

In each surviving composition case the Jacobian would be divisible
by the nonconstant derivative \(T_w\) of a quadratic polynomial \(T\).
If \(P=U(u,T), Q=V(u,T)\), the terms involving \(T_u\) cancel
and the determinant is \(T_w(U_u V_T-U_T V_u)\). All
coefficients can be taken in a finite algebraic extension of
\(K\), where \(d/du\) extends. This contradicts (1). Thus the
monodromy is primitive. The conjugates of a transposition form
the edges of a graph whose connected components are blocks.
Primitivity makes the graph connected, and its transpositions
generate \(S_{24}\).

For \(m>2\) with nonzero central restrictions, (12) still imposes
exactness and excludes any case in which the differential has a
nonzero residue. The constant bracket in (10) also makes
\(q_0^a\) and \(\lambda p_0^b\) a central Davenport–Zannier pair.
They have the same leading coefficient by the leading relation
in the generic ruling, and their difference has extremal degree
by (3) in the variable \(z\). At \((r,s)=(8,12)\), their ratio
defines a degree-12 cover with partitions \(2^6,3^4,(9,1^3)\).
The higher-order cases that pass these tests and the larger
excesses remain open here.
