# Cubic and quadratic mates of the squarefree quartic scalar

Let \(v=xy-1\) and \(R=xvz+b(x,y)\), with the
[squarefree quartic normal form](QUARTIC_AFFINE_SOURCE_COMPONENTS_THROUGH_DEGREE_SEVEN.md)
for \(b\). In particular \(b(0,y)=\lambda_0y+\mu_0\)
with \(\lambda_0\ne0\). By the
[cubic leading reduction](SEXTIC_KELLER_SQUAREFREE_REDUCED_CUBIC_CASES.md),
the two remaining patterns for cubic and quadratic source-direction
mates \(P,Q\) are \(p=q=xv\), or \(p=x^3\) and
\(q=x^2(\gamma v+\delta)\). Here \(p\) and \(q\) are the leading
coefficients in \(z\), after rescaling.

**Theorem SexticKellerSquarefreeQuarticCommonCubicQuadraticLeadingExclusion.**
In characteristic zero and ordinary degree at most six, the leading
pair \(p=q=x(xy-1)\), with mate \(z\)-degrees \(3,2\), is impossible
for the squarefree quartic normal form.

**Theorem SexticKellerSquarefreeQuarticMixedCubicQuadraticLeadingExclusion.**
In the same setting, the leading pair \(p=x^3\),
\(q=x^2(\gamma(xy-1)+\delta)\) is impossible if
\(\gamma\ne0\) and \(\gamma\ne2\delta\).
Thus the remaining quadratic leading coefficients in this pattern,
up to nonzero scalars, are \(x^2\) and \(x^2(2xy-1)\).

## Expansion in the target parameter

Work over \(\mathbb C\) by embedding the finitely generated
coefficient field. In the first pattern put \(t=xv\), and in the
second put \(t=v\). In both put \(s=1/x\) and use \(r=R\) as a
coordinate on \(xv\ne0\).
In \(t,s,r\) the target Jacobian is
\(\{P,Q\}_{t,s}=\kappa/t\), with \(\kappa\in\mathbb C^\times\).
Indeed, \(y=ts^2+s\) in the first case and \(y=(t+1)s\)
in the second. The respective Jacobians of \((x,y)\) in \((t,s)\)
are \(1\) and \(1/s\), and \(z_r\) is \(1/t\) and \(s/t\),
respectively.

Write \(\epsilon=r^{-1}\) and set \(U=Q/r^2,\ V=P/r^3\). At
fixed \(r\),
\[
 U=\ell+\epsilon u+\epsilon^2 w,\quad
 V=m+\epsilon V_1+\epsilon^2 V_2+\epsilon^3 V_3,\qquad
 \{V,U\}_{t,s}=\epsilon^5\kappa/t .
\]
Here \((\ell,m)=(1/t,1/t^2)\) in the common-leading case and
\((\ell,m)=((\gamma t+\delta)/t^2,1/t^3)\) in the other case.
Put \(K=\mathbb C(t)\), and let \(\partial=(d\ell/dt)^{-1}d/dt\)
on \(K\). We use primes for this derivation. For \(A\in K\)
the notation \(A[U]\) denotes its Taylor series
\(\sum_{i\ge0} A^{(i)}(U-\ell)^i/i!\). The formal chain rule
gives \(\{A[U],U\}=0\). For instance, \(dA=A' d\ell\), so
differentiating term by term yields \(d(A[U])=A'[U]\,dU\).

Successively comparing the first four coefficients of the bracket
gives \(A,B,C,D\in K\) such that
\[
 V=m[U]+\epsilon A[U]+\epsilon^2 B[U]
            +\epsilon^3 C[U]+\epsilon^4 D[U]\pmod{\epsilon^5}.
 \tag{1}
\]
At each step, the coefficient of the residual has zero bracket
with \(\ell\), so is independent of \(s\). This proves the assertion
without any choice of root of \(\ell\) or convergence assumption.
The coefficient of \(\epsilon^4\) on the right is zero. Write
\(G(u,w)\) for that coefficient and \(H(u,w)\) for the coefficient
of \(\epsilon^5\). Taking the bracket of the residual
\(-\epsilon^5 H+O(\epsilon^6)\) with \(U\) gives
\[
 G(u,w)=0,\qquad H(u,w)=\frac{\kappa}{t(d\ell/dt)}s+E(t),
 \quad E\in K .
 \tag{2}
\]
In particular \(H(u,w)\) is a rational function of degree one in \(s\).

## Both leading coefficients equal to \(xv\)

Write \(Q=tz^2+q_1z+q_0\) and
\(P=tz^3+p_2z^2+p_1z+p_0\). Then \(tu=q_1-2b\).
The coefficient of \(\epsilon\) in (1) gives
\(t^2 A=p_2-2q_1+b\). This belongs to \(\mathbb C[t]\), since
it is in \(K\) and, after substituting \(x=1/s,\ y=ts^2+s\),
in \(\mathbb C[t,s,s^{-1}]\). Its ordinary degree in \(x,y\)
is at most five. As \(\deg t=3\), it is affine in \(t\).
Thus \(a_0=A''\) is a scalar and \(A'''=0\).

The two coefficients in (2) are
\[
 \begin{aligned}
 G &=w^2+(a_0u+B')w+B''u^2/2+C'u+D,\\
 H &=a_0w^2/2+B''uw+B'''u^3/6+C'w+C''u^2/2+D'u .
 \end{aligned}
\]
Complete the square by putting \(Z=w+(a_0u+B')/2\). Then
\[
 Z^2=T(u)=L u^2+M u+N,\qquad
 L=a_0^2/4-B''/2,\quad M=a_0B'/2-C',\quad N=B'^2/4-D .
\]
Modulo this relation, \(H=-(2Lu+M)Z+J(u)\). The coefficients of
\(u^3,u^2\) in \(J\) are \(B'''/6\) and
\(a_0^3/4-3a_0B''/4+C''/2\).

If \(u\) is constant in \(s\), the equation for \(Z^2\) gives
the same for \(Z\), contradicting (2). At any pole of \(u\)
over \(\overline K\), the pole order of \(Z\) is at most that
of \(u\). If \(B'''\ne0\), the cubic term of \(J\) gives \(H\)
a pole of order at least three, contradicting (2). Hence
\(B''\) and \(L\) are scalars.

If \(L\ne0\) and \(T\) is nonsquare over \(\overline K\), the
quadratic extension \(\overline K(u,Z)\) has two places over
\(u=\infty\). At each, \(u\) has a simple pole, and \(Z/u\)
has opposite nonzero limits \(\pm\sqrt L\). The leading
coefficient of \(H/u^2\) is different on the two places, so
\(H\) has at least one double pole. This contradicts (2),
because (2) gives \(\overline K(u,Z)=\overline K(s)\).
If \(T\) is square, \(Z=\sigma u+\theta(t)\), with \(\sigma^2=L\).
Here \(\theta\in K\) since \(u,Z\in K(s)\) and \(\sigma\) is
a scalar. The same form with \(\sigma=0\) holds if \(L=M=0\),
as in that case \(Z\) is independent of \(s\). In either
case \(H\) is a polynomial in \(u\). Degrees of rational
functions multiply under composition, so (2) forces \(u\)
to be affine in \(s\), and \(w=k u+j(t)\) for a scalar \(k\).

This last possibility violates polynomiality in the original
coordinates. Write \(tu=A_1(t)s+A_2(t)\). Since \(tu=q_1-2b\)
is polynomial, \(A_1,A_2\in\mathbb C[t]\). Regularity along
\(x=0\) forces \(A_1(0)=0\), so \(tu\) restricts there to the
constant \(-A_1'(0)+A_2(0)\), where this derivative is in \(t\).
But
\[
 tQ=R^2+(R+k)(tu)+t j(t).
\]
The last term is polynomial in the source, hence polynomial
in \(t\) by the same Laurent-ring intersection. At \(x=0\) this
identity gives a monic quadratic relation over \(\mathbb C\)
for the nonconstant \(\lambda_0y+\mu_0\), an impossibility.

Finally let \(L=0,\ M\ne0\). Now \(u=(Z^2-N)/M\), and (2) and
the formula \(H=-MZ+J(u)\) force \(J\) to be constant in \(u\).
Indeed any positive power of \(u\) in \(J\) contributes an
uncancelled even degree in \(Z\), larger than one.
Since \(B''=a_0^2/2\), vanishing of the \(u^2\)-coefficient
gives \(C''=a_0^3/4\), so \(M'=0\). Hence \(M\) is a
nonzero scalar, and \(Z\) is a nonzero scalar multiple of
\(ts=v\) plus an element of \(K\), by (2).
Writing \(B'=a_0^2/(2t)+b_0\) with \(b_0\) constant gives
\[
 Q=(R^2-a_0^2/4)/t+(R-a_0/2)u+Z-b_0/2 .
\]
On \(R=a_0/2\) this restricts to \(Q=\alpha v+j(t)\),
\(\alpha\in\mathbb C^\times\).

On any fibre \(R=r_0\), the open part \(tv\ne0\) has coordinates
\(t,v\), with \(s=v/t\). There \(\{P,Q\}_{t,v}=\kappa/t^2\).
If \(Q=\alpha v+j(t)\), holding \(Q\) fixed and integrating in
\(t\) gives \(P=-\kappa/(\alpha t)+g(Q)\) with
\(g\in\mathbb C(Q)\). The function \(j\) is regular on
\(t\ne0\), since all pairs \(t,v\ne0\) occur on this open
part. The fibre also contains the line \(x=0,\ b(0,y)=r_0\).
At its generic point \(t\) has a simple zero because
\(R_y=\lambda_0\ne0\) there and \(v=-1\). Regularity of \(Q\)
there gives regularity of \(j\) at zero. The value
\(c=j(0)-\alpha\) occurs in the open part: choose \(t\ne0\)
with \(j(t)\ne c\), and take \(v=(c-j(t))/\alpha\).
Thus \(g\) is regular at \(c\), since \(P\) and \(1/t\)
are regular on that open part and \(Q_v=\alpha\ne0\).
On the line \(x=0\), \(g(Q)\) then cannot cancel the simple
pole of \(1/t\). This completes the common-leading exclusion.

## A quadratic leading coefficient with nonzero \(\gamma\)

Let \(f=\gamma t+\delta,\ d=\gamma t+2\delta\), \(\gamma\ne0\).
Here \(G\) is quadratic in \(w\), with nonzero leading
coefficient \(m''/2\). Its terms of weight four
for \({\rm wt}(u,w)=(1,2)\) are
\[
 m''w^2/2+m'''u^2w/2+m^{(4)}u^4/24.
\]
The weight-five part of \(H\) is
\(m'''uw^2/2+m^{(4)}u^3w/6+m^{(5)}u^5/120\).
For \(w=u^2 h\) set \(I=\delta t^2-d^2h\) and
\(J_0=2\delta^2+2\delta\gamma t-\gamma^2 t^2\).
These parts, after division by \(u^4,u^5\) respectively, are
\[
 \frac{3t}{d^7} I(fI-\gamma^2t^4),\qquad
 \frac{3t^3}{d^9} I(-J_0 I+2\delta\gamma^2t^4).
 \tag{3}
\]
The first has two distinct roots. Only \(I=0\) is common
to the two expressions, since
\(2\delta f-J_0=\gamma^2t^2\).

By (2), \(u\) is nonconstant in \(s\), since otherwise the
quadratic equation would also make \(w\) constant. If \(G\)
were irreducible over \(\overline K(u)\), there would be two
places over \(u=\infty\), with \(w/u^2\) tending to the two
distinct roots of the first polynomial in (3). This follows
by solving for \(w/u^2\) over \(\overline K[[1/u]]\) at
each simple root. At one of the places \(H\) has pole order
five by (3), contradicting the fact that
\(\overline K(u,w)=\overline K(s)\) by (2).

Consequently \(w\) is polynomial of degree at most two in
\(u\). Indeed a rational root of the monic quadratic over
\(\overline K[u]\) is polynomial, and its degree bound follows
from the weights of \(G\). Its coefficients descend to \(K\)
since \(u,w\in K(s)\) and \(u\) is transcendental over \(K\).
Again by (2), \(u=\alpha(t)s+\beta(t)\) with \(\alpha\ne0\).
Writing \(w=\rho u^2+\nu u+\omega\), the coefficient of \(u^5\)
in \(H\) must vanish, so (3) gives
\[
                         \rho=\delta t^2/d^2 .
\]
All coefficients of \(U,V\), as Laurent polynomials in \(s\),
belong to \(\mathbb C[t,t^{-1}]\). In particular \(\alpha,\beta\)
are regular at \(t=-1\).

Suppose now \(d(-1)\ne0\). If \(\alpha(-1)\ne0\), the
coefficients \(\nu,\omega\) are also regular at \(-1\), by
changing from \(s\) to \(u\) in \(w\). At the generic point of
the original divisor \(x=0\), \(t=-1+xy\) and \(s=1/x\), so
\(u\) has a simple pole while \(R\) is regular. If \(\delta\ne0\)
then \(Q=\ell R^2+Ru+w\) has an uncancelled double pole
from \(\rho u^2\). If \(\delta=0\) it has an uncancelled
simple pole from \((R+\nu)u\), since \(R|_{x=0}\) is
nonconstant. Both contradict polynomiality. Therefore
\(\alpha(-1)=0\), which makes \(u\) regular along \(x=0\).
The expression \(w=Q-\ell R^2-Ru\) is regular there as well.

Now \(m'\) and \(m''\) are regular in \(t\) at \(-1\)
because \(d(-1)\ne0\). In (1), \(V_1=m'u+A\) shows
that \(A\) is regular in \(t\) there, using the coefficients
in the Laurent variable \(s\). Then \(A'\) is also regular.
Similarly \(V_2=m'w+m''u^2/2+A'u+B\) shows that \(B\)
is regular there. Thus \(V_1,V_2\) are both regular along
\(x=0\) in the original source.

Write \(P=x^3z^3+p_2z^2+p_1z+p_0\). Substituting
\(z=(r-b)/(xt)\) into \(P/r^3\) gives
\[
 V_1=p_2/(x^2t^2)-3b/t^3,\qquad
 V_2=p_1/(xt)-2b\,p_2/(x^2t^2)+3b^2/t^3 .
\]
Since \(t=-1\) is a unit along \(x=0\), regularity of \(V_1\)
first gives \(x^2\mid p_2\). Regularity of \(V_2\) then gives
\(x\mid p_1\). Moreover \(u=q_1/(xt)-2fb/t^2\) for \(Q\)'s
linear \(z\)-coefficient \(q_1\), so \(x\mid q_1\).
Consequently \(P_z,Q_z,R_z\) all vanish on \(x=0\),
contradicting the nonzero Jacobian determinant.
