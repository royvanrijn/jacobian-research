# Degree-five obstruction for squared quadratic-height profiles

Let \(K\) be a field of characteristic zero and let

\[
S=K[x,y,t]/(t^2+2xy-1),\qquad
\Lambda(f)=\operatorname{CT}_x\frac12\int_{-1}^1
         f\left(x,\frac{1-t^2}{2x},t\right)dt.
\]

The interval integral is the rational polynomial functional. A phase is pure if every positive power has zero trace. It is safe if for every fixed marker \(a\in S\), the moments \(\Lambda(af^m)\) vanish for all sufficiently large \(m\). All degrees use the polynomial filtration of \(S\).

**Theorem.** Let \(f=R(x)g^2\), where \(R\in K[x]\) and the literal Laurent expression for \(g\in S\) has \(t\)-degree at most two. If \(f\) is pure and has least degree five, it is safe. If \(g\) has a pole at \(x=0\) and \(R(0)\ne0\), the degree-five phase \(f\) cannot be pure.

In particular, this excludes the regular quintic sector of literal height four and pole order two whose quartic numerator has two double roots over \(\overline{K(x)}\). In that sector write the numerator as \(R(x)G(x,t)^2\), with \(R\in K[x]\) and \(G\in K[x,t]\) by taking the primitive polynomial part over \(K[x]\). Regularity forces the constant \(x\)-coefficient of the numerator to be a nonzero multiple of \((1-t^2)^2\), so \(R(0)\ne0\) and \(g=G/x\) is regular. Height four with a generically nonsquare numerator is a separate case.

## The diagonal and the trace

The normal monomials are \(y^jt^r\) for \(j\ge0\) and \(x^it^r\) for \(i\ge1\). Their leading parts are still independent in the associated graded domain \(K[x,y,t]/(t^2+2xy)\). A coefficient of \(x^{-j}\), \(j\ge1\), in a regular Laurent expression is divisible by \((1-t^2)^j\). Thus a \(g\) of literal height at most two has at most a simple pole at \(x=0\), and with a pole it is

\[
g=x^{-1}(A-Bt^2+Ct),\qquad A(0)=B(0)\ne0,\quad C(0)=0,
\qquad A,B,C\in K[x].
\]

Scale the phase and root to arrange \(A(0)=B(0)=1\), and also \(R(0)=1\) when it is nonzero. Put

\[
D=A+\frac{C^2}{4B},\qquad L=\sqrt{D/B}\in K[[x]],\quad L(0)=1.
\]

Write \(a=\operatorname{ord}_{x=0}R\), \(R=x^a R_0\), \(R_0(0)\ne0\). Completing the square with \(z=(t-C/(2B))/L\) gives

\[
\Lambda(f^m)=c_{2m}[x^{(2-a)m}] (R_0D^2)^m L,
\qquad
c_{2m}=\frac12\int_{-1}^1(1-z^2)^{2m}\,dz\ne0.                 \tag{1}
\]

Indeed, the two transformed endpoints are congruent to \(1,-1\) modulo \(x\). Each endpoint displacement changes the primitive of \((1-z^2)^{2m}\) by a multiple of \(x^{2m+1}\). This gives (1), with negative requested exponents interpreted as zero. Formal changes of variables in these integrals mean substitution in a polynomial primitive.

Suppose \(a=0\). Set

\[
u=\frac{x}{\sqrt R\,D},\qquad \sigma=u^2=\frac{x^2}{RD^2},
\qquad \rho=\frac{2\sigma}{x\sigma'}L.
\]

The square root of \(R\) has constant term one. The map \(u\) has formal inverse \(v\) with \(v(0)=0\). Formal residue substitution in (1) gives

\[
[x^{2m}](RD^2)^m L=[u^{2m}]\rho(v(u)).                          \tag{2}
\]

In fact \(u'(x)/u(x)=\sigma'(x)/(2\sigma(x))\), so the substituted multiplier for the coefficient residue is \(L(v(u))\,u v'(u)/v(u)\). The substitution follows on monomials from vanishing residues of derivatives and \(\operatorname{Res}(u'/u)dx=1\).

Purity and \(\rho(0)=1\) show that the two inverse branches \(x=v(u),v(-u)\) of \(\sigma=u^2\), with the chosen evaluations of \(L\), satisfy

\[
\rho(v(u))+\rho(v(-u))=2.                                      \tag{3}
\]

Both terms of (3) are conjugates over \(K(\sigma)\) in a normal closure of \(K(x,L)/K(\sigma)\). Taking the normalized field trace gives normalized trace one for \(\rho\). If \(L\notin K(x)\), the quadratic conjugation over \(K(x)\) sends \(L\) to \(-L\) and therefore sends \(\rho\) to \(-\rho\). Transitivity of trace gives trace zero, a contradiction. It follows that \(BD\) is a polynomial square \(J^2\), with \(J(0)=1\), because \(BD=AB+C^2/4\in K[x]\) and \(BL=\sqrt{BD}\) is rational. In particular \(L=J/B\).

Let \(d_\sigma=[K(x):K(\sigma)]>0\). Equation (3) now gives the trace of a one-form on the projective \(x\)-line:

\[
\omega=\frac{J}{xB}\,dx,qquad
\operatorname{Tr}_{K(x)/K(\sigma)}(\omega)
             =\frac{d_\sigma}{2}\frac{d\sigma}{\sigma}.          \tag{4}
\]

Here the trace of \(h\,d\sigma\) means the field trace of \(h\) times \(d\sigma\). Indeed \(\rho=2\sigma\,\omega/d\sigma\). After extending the constants to an algebraic closure, the residue of the trace at \(\sigma=0\) is the sum of the residues of \(\omega\) at the points above zero. To see this locally, choose a formal coordinate \(s\) at any such point with \(\sigma=s^e\), where \(e\) is the ramification index. The \(e\) inverse branches substitute \(s=\zeta\sigma^{1/e}\). Each contributes \(1/e\) of the residue from the term \(ds/s\). Terms of other powers of \(s\) contribute no residue. Thus if all the poles of \(\omega\) map to zero, the residue theorem on the projective line makes the residue zero, contradicting (4).

## The two degree-five cases with nonzero \(R(0)\)

Degree is additive because the associated graded ring is a domain. With a pole, \(\deg g\ge1\), so \(\deg f=5\) gives the two possibilities \((\deg g,\deg R)=(1,3),(2,1)\). The normal monomials give the bounds

\[
\deg A\le\deg g+1,\qquad \deg B\le\deg g-1,
\qquad \deg C\le\deg g.
\]

* For \(\deg g=1\), we have \(B=1\) and \(\deg J\le1\). Now \(\sigma=x^2/(RJ^4)\). Since \(\deg R=3\), infinity maps to zero. The only possible poles of \(\omega=(J/x)dx\) are zero and infinity, and both map to zero. This contradicts (4).
* For \(\deg g=2\) we have \(\deg B\le1\), \(\deg J\le2\), and \(\deg R=1\). Now \(\sigma=x^2B^2/(RJ^4)\). If \(\deg J\ge1\), infinity maps to zero because the numerator has degree at most four and the denominator has degree at least five. A finite pole of \(\omega\) other than zero can occur only at a root \(\beta\) of the linear \(B\) with \(J(\beta)\ne0\). Then \(\operatorname{ord}_\beta\sigma=2-\operatorname{ord}_\beta R\ge1\), so \(\beta\) also maps to zero. This contradicts (4). If \(J=1\) and \(B\) is linear, \(\omega\) has poles only at zero and \(\beta\), again both over zero. The only remaining case is \(J=B=1\). Then \(D=1\) and \(A=1-C^2/4\). Since \(\deg A\le3\), we have \(\deg C\le1\) and \(g\) has degree one. This contradicts \(\deg g=2\).

## The singular charts and safety

If \(g\) has no pole, \(f\in K[x,t]\). For such a phase, every power moment equals the interval moment of \(f(0,t)^m\). The interval theorem of [Müger–Tuset, Theorem 2.1](https://arxiv.org/abs/1910.12148v1) makes purity possible only when \(f(0,t)=0\). The theorem applies after embedding the finitely generated coefficient field in \(\mathbb C\) and affinely transforming the interval. The phase then has strictly positive \(x\)-exponents and is safe, since a fixed marker has bounded negative \(x\)-exponents.

Suppose \(g\) has a pole and \(a\ge2\). If \(a=2\), (1) has a nonzero first moment. If \(a\ge3\), \(f\) has only strictly positive Laurent \(x\)-exponents and is safe.

It remains to treat \(a=1\). Set \(M=R_0D^2\). Equation (1) shows that purity is equivalent to \([x^m]M^mL=0\) for every \(m\ge1\). These conditions determine \(L\) uniquely among series with constant term one, since the next coefficient has multiplier \(M(0)^m\ne0\). The identity \((M^m)'=mM^{m-1}M'\) shows that the unique solution is \(L=1-xM'/M\).

* If \(\deg g=1\), then \(B=1\), \(D\) is polynomial and \(L^2=D\) forces the rational \(L\) to be polynomial. A root of the polynomial \(M\) would give a pole of \(1-xM'/M\), since \(M(0)\ne0\). Thus \(M=R_0D^2\) is constant, contradicting \(\deg R_0=2\).
* If \(\deg g=2\), then \(R_0\) is constant. Squaring the formula for \(L\) gives \(B=D^3/(D-2xD')^2\). At a finite zero or pole of \(D\) of order \(r\ne0\), \(D-2xD'\) has order \(r-1\) and \(B\) has order \(r+2\). All these points are nonzero, and every other finite point has order zero in \(B\) because \(B\) is polynomial. If \(n\) is the degree of \(D\) at infinity, the leading term of \(D-2xD'\) is \(1-2n\) times that of \(D\), so \(\deg B=n\). The finite valuations instead sum to \(n+2s\), where \(s\) is the number of distinct finite zeros and poles of \(D\). It follows that \(D\) is constant and then \(B=1\). As before, \(A=1-C^2/4\) gives \(\deg C\le1\), contradicting \(\deg g=2\).

This proves the theorem over every characteristic-zero field. The normalization and transfer of safety to arbitrary fixed original markers are those of the [quadratic-height classification](QUADRATIC_HEIGHT_HOPF_CLASSIFICATION.md#5-exact-least-degree-parity-and-the-image-normalization). Its endpoint calculation for a single quadratic phase is the closest earlier case. Here (2) supplies only the even diagonal, so (3) and the trace-residue argument provide the additional obstruction. No finite moment computation is needed for the theorem.
