# Hall's conjecture through cubic fields and Padé approximation

Write \(k=y^2-x^3=dm^2\ne0\), with \(x,y,m>0\) and \(d\) signed
squarefree. The cubic \(P(T)=T^3-3xT-2y\) has discriminant \(-108k\).
Throughout, primitive means \(\gcd(x,y)=1\). If \(k\) is squarefree,
primitivity is automatic since \(\gcd(x,y)^2\mid k\).

Three bounds isolate parts of the height problem. In the primitive
irreducible squarefree case the index of the displayed generator divides
\(6\). A sufficiently good point in a quadratic descent branch gives an
effective height ceiling for further points in the same branch and angular
sector. A reduced binary cubic in the same root field gives an effective
bound for the first point with separate dependence on \(d\) and \(m\).
The desired uniform power bound on that first height remains open.

## Reducible cubics and the cyclic squareclass

An integer root \(r\) gives \(2y=r^3-3xr\), so
\[
4k=(r^2-x)^2(r^2-4x).
\]
For \(a=r^2-x\), the cases \(a\le-1\) and \(a\ge3x+1\) give
\(|a^2(a-3x)|\ge3x-1\) directly. On \(1\le a\le3x-1\),
\(a^2(3x-a)\) has its minimum at an endpoint. Thus
\[
 |k|\ge(3x-1)/4.
\]
This requires no coprimality. Equality occurs at
\(x=4s^2-1,\ y=8s^3-3s,\ r=-2s\), with \(s\ge1\).
These examples are primitive since \(y=s(2x-1)\).

For \(d\ne1\), put \(E=\mathbb Q(\sqrt d)\) and
\(\eta=y+m\sqrt d\). The cubic is reducible exactly when \(\eta\)
is a cube in \(E\). If \(\eta=\gamma^3\), then \(N\gamma=x\) and
\(\gamma+\bar\gamma\) is a root. Conversely, if \(r\) is a rational
root, take \(z,w\) with \(z+w=r,\ zw=x\). The unordered pair of their
cubes is \(\{\eta,\bar\eta\}\). Since these cubes generate \(E\), the
quadratic polynomial with roots \(z,w\) already splits over \(E\).
This improves the earlier \(x/32\) bound wherever the descent proves
that \(\eta\) is a cube.

In fact the entire primitive cyclic squareclass \(d=-3\) is reducible.
Here \(x\) is odd, since \(x\) even and \(y\) odd make
\(y^2+3m^2\) nonzero modulo \(8\). The ideals of \(\eta\) and
\(\bar\eta\) are coprime in the Eisenstein integers. A common odd
prime would divide \(y\) and \(x\), and \(2\) is excluded by the odd
norm. The class number is one, so \(\eta=\epsilon\gamma^3\).
Modulo \(2\) the residue field is \(\mathbb F_4\).
There \(\eta=y+m=1\) and \(\gamma^3=1\). The Eisenstein units
congruent to \(1\) modulo \(2\) are \(\pm1\), both cubes.
An irreducible cubic with square discriminant would require \(d=-3\).
Thus a primitive irreducible \(P\) always has non-Galois root field.

## Field discriminant and the generator index

Suppose \(P\) is irreducible. Let \(L=\mathbb Q(\theta)\) for a root
\(\theta\), let \(D_L\) be its signed field discriminant, and write
\(I=[{\cal O}_L:\mathbb Z[\theta]]\). Put \(K=\mathbb Q(\sqrt{-3d})\)
and write \(D_K\) for its fundamental discriminant. Use \(D_K=1\)
if \(-3d\) is a square. The discriminant identity gives
\[
 D_L=D_K f^2,\qquad
 I f=m\sqrt{108|d|/|D_K|}.
\]
The square root on the right lies in \(\{3,6,9,18\}\).
For a primitive point, the local calculation below gives
\[
 f\mid18,\qquad m\mid I,\qquad I\mid18m.
\]
For squarefree \(k\) it gives \(I\mid6\). In particular
\(|D_L|\asymp |d|\) uniformly in the primitive irreducible case.

At \(p\ge5\) with \(p\mid k\), primitivity gives \(p\nmid xy\).
Modulo \(p\), \(P=(T-r)^2(T+2r)\) with \(r=-y/x\).
The simple factor lifts and the remaining algebra is quadratic.
Thus \(v_p(D_L)=v_p(D_K)\), which equals \(1\) if \(p\mid d\)
and \(0\) otherwise. When \(p\nmid k\), the polynomial discriminant
is already prime to \(p\).

At \(2\), odd \(x\) gives \(P\equiv T(T+1)^2\) and the same
quadratic-algebra argument gives \(v_2(D_L)=v_2(D_K)\).
Even \(x\) gives odd \(y\) and an Eisenstein cubic with tame
discriminant exponent \(2\). Then \(k\equiv1\bmod8\) and \(D_K\)
is odd.

At \(3\), take \(s\equiv-y\bmod3\), \(w=s^2-x\), \(C=P(s)\), and
\(j=v_3(C)\ge1\). The translated polynomial is
\(T^3+3sT^2+3wT+C\).

- If \(3\nmid k\) then \(w\) is a unit. For \(j=1\) this is
  Eisenstein and its derivative has uniformizer valuation \(3\).
  For \(j\ge2\), \(P(s+3U)/9\) is linear with unit derivative
  modulo \(3\) and gives a local root.
- If \(3\mid k\), then \(s\) is a unit and \(3\mid w\). For \(j=1\),
  the Eisenstein derivative has uniformizer valuation \(4\).
  For \(j=2\), the translated root \(\rho\) has valuation \(2\)
  in the totally ramified cubic, normalized by \(v(3)=3\).
  The uniformizer \(\pi=\rho^2/3\) has Eisenstein polynomial
  \[
  T^3-(3s^2-2w)T^2+(w^2-2sC/3)T-C^2/27 .
  \]
  Its derivative has valuation \(3\).
  For \(j\ge3\), \(P(s+3U)/27\) is monic with quadratic
  coefficient \(s\not\equiv0\bmod3\). Its reduction is separable
  or has a simple linear factor.

Whenever there is a local linear factor, the remaining quadratic
discriminant is the resolvent discriminant, so \(v_3(D_L)=v_3(D_K)\).
When \(3\mid k\) one can choose \(s=\pm1\), write \(x=1+3a\), and
see that \(j=1, j=2, j\ge3\) correspond to \(v_3(k)=1,2,\ge3\).
Indeed \(C=-2(s+y)-9as\). For \(j\ge2\), writing \(y=-s+9b\)
gives \(k/9\equiv-s(as+2b)\bmod3\).
These cases and the values of \(108|d|/|D_K|\) give the index
divisibilities above. In particular when \(v_3(k)\le1\) the index
valuation at \(3\) is at most one.

An additional effective exclusion applies when \(x\) is odd and
\(P\) has a \(\mathbb Q_3\)-root. In that case \(f=1\).
For irreducible \(P\), its Galois closure would be an everywhere
unramified cubic extension of \(K\). If the classes over split
primes \(p\ge5\) dividing \(ym\) span \(\mathrm{Cl}(K)/3\mathrm{Cl}(K)\),
this is impossible. For \(p\mid y\) the reduction has a linear
factor and is separable, and for \(p\mid m\) the simple factor in
the double-root reduction lifts. In both cases splitting in \(K\)
then forces complete local splitting. The resulting order-three
class character would kill the spanning classes. This test includes
\(\mathrm{Cl}(K)[3]=0\), and applies to the positive squareclasses
as well. The condition at \(3\) holds, for example, when \(27\mid k\)
and \(x\) is odd.

For a cubic index form \(J\) of the maximal order and integral
coordinates \(u,v\) of \(\theta\) modulo integers, one has
\(|J(u,v)|=I\) and \(\operatorname{Hess}(J)(u,v)=9x\), with the
Hessian normalized to have leading coefficient \(b^2-3ac\) for
\(aU^3+bU^2V+cUV^2+eV^3\). In the basis \(1,\theta,\theta^2\),
the signed index form relative to that order is
\(\prod_i(U-V\theta_i)=P(U/V)V^3\), and its Hessian at
\((1,0)\) is \(9x\). Changing to a maximal-order basis multiplies
the form by \(I\) after substitution of a change of determinant
\(\pm1/I\). The two factors cancel in the Hessian transformation.
Thus squarefree Hall here asks for
\(\operatorname{Hess}(J)(u,v)\ll_\varepsilon |D_L|^{2+\varepsilon}\)
on actual solutions of \(|J|=1,2,3,6\) with the Hall local data.
This is a height estimate on the specified generator.
[Bennett–Ghadermarzi, Section 2](https://arxiv.org/abs/1311.7077) ·
[Belabas, Theorem 1.1 and Lemma 1.6](https://www.math.u-bordeaux.fr/~kbelabas/research/mathcomp97.pdf) ·
[Gaál–Pohst–Pohst, Lemmas 2–3](https://agag-fieker.math.rptu.de/~fieker/pub/gaal13.pdf)

## Reduced descent forms and an effective first height

Write \((\eta)=\mathfrak b\mathfrak a^3\) by reducing ideal valuations
modulo \(3\). For a fixed \(\mathfrak b\) and ideal class of
\(\mathfrak a\), choose an integral representative \(\mathfrak a_0\).
Fix representatives of units modulo cubes. This gives branches
\(\eta=\beta\gamma^3\) with \(\gamma\in\mathfrak a_0^{-1}\) and
\((\beta)=\mathfrak b\mathfrak a_0^3\). In an integral basis
\(\gamma=\gamma(u,v)\) of that fractional ideal put
\[
 F(u,v)=\frac{\beta\gamma^3-\bar\beta\bar\gamma^3}{\sqrt d},
 \qquad Q(u,v)=(N\beta)^{1/3}N\gamma .
\]
The cube root in \(Q\) is signed. Then \(F,Q\) have integral
coefficients, \(F(u,v)=2m\) and \(Q(u,v)=x\) at the Hall point.
Indeed \(F\) is twice the \(\sqrt d\) coefficient of an algebraic
integer, and \(Q\) is \(\pm(N\mathfrak b)^{1/3}N\mathfrak a_0\,N\gamma\),
an integer multiple of an integral ideal norm form. The norm
\(N\mathfrak b\) is a cube since \(N\eta=x^3\).

Let \(e=D_E/d\in\{1,4\}\) and \(r=(N\mathfrak b)^{1/3}\). Scaling the
imaginary quadratic ideal lattice by \(\beta^{1/3}\) gives covolume
\(r\sqrt{|D_E|}/2\). In the real case scaling each embedding by its
signed cube root gives covolume \(r\sqrt{|D_E|}\). Applying the
binary cubic discriminant and Hessian transformation rules to
\(2(3X^2Y-Y^3)/\sqrt{|d|}\) in the imaginary case and
\((X^3-Y^3)/\sqrt d\) in the real case gives
\[
 \operatorname{disc}(F)=-27r^6e^3d,\qquad
 \operatorname{Hess}(F)=9r^2 e\, Q.
\]

The root fields of \(F\) and \(P\) coincide. Multiplying the
quadratic coordinates by the actual \(\gamma\) changes the root
equation to \(\eta h^3=\bar\eta\bar h^3\). Write
\(h=1+t\sqrt d\), choose \(z^3=\eta,\ w^3=\bar\eta,\ zw=x\), and
put \(U=z+w\). Then
\[
 t=\frac{w-z}{U\sqrt d}
   =-\frac{2m}{U(U^2-x)},\qquad
 U^2=\frac{4x}{1-dt^2}.
\]
These formulas recover each quantity rationally from the other.
The denominators are nonzero at these roots because \(y,m>0\).
Consequently Bugeaud–Győry's regulator bound for solutions of
\(F(u,v)=2m\) uses the regulator \(R_L\) of the field in the
generator formulation.

For a primitive point \(r\le2\). Reducing \(F\) by a unimodular
substitution gives a coefficient height polynomial in \(D=|d|\).
The transformed Hessian has polynomially bounded coefficients too,
so \(x\ll D^{O(1)}\max(1,|u|,|v|)^2\) in those coordinates.
Theorem 3 of Bugeaud–Győry now gives
\[
 \log(2x)\ll \log(2D)
  +(1+R_L)\log(2+R_L)\bigl(R_L+\log(2Dm)\bigr).
\]
The class number formula and a cubic zeta residue bound give
\(R_L\ll\sqrt D\log^2(2D)\). With \(L_D=\log(2D)\) this yields
\[
 \log(2x)\ll D L_D^5+\sqrt D\,L_D^3\log(2m).
\]
This effective bound applies to the first point and retains the
square part separately. For bounded \(m\), von Känel's general
Mordell height bound has the stronger logarithmic dependence.
[von Känel, Corollary 7.4](https://academic.oup.com/tlms/article/1/1/85/2841113)

For a nonprimitive point put \(G=\operatorname{rad}(\gcd(x,y))\).
The common ideal factors lie above \(2G\), so \(r\le2G\).
At \(p\ge5\), a cubic field discriminant exponent is at most \(2\);
an extra exponent \(2\) over the quadratic resolvent can only
occur at \(p\mid G\), by the simple-factor argument. Thus
\(|D_L|\ll D G^2\). The same proof bounds \(\log(2x)\) by
\[
 D G^2\log^5(2DG)+\sqrt D\,G\log^3(2DG)\log(2m)
\]
up to an effective absolute factor. Both regulator bounds concern
irreducible \(P\) and \(d\ne1\).
[Cremona on binary cubic reduction](https://doi.org/10.1112/S1461157000000073) ·
[Bugeaud–Győry, Theorem 3](https://matwbn.icm.edu.pl/ksiazki/aa/aa74/aa7437.pdf#page=3) ·
[Louboutin and the class number formula](https://www.numdam.org/article/JTNB_2005__17_2_559_0.pdf)


A stronger effective estimate holds in a monic Thue sector.
Suppose \(d<0\), \(d\equiv2,3\bmod4\), \(m=1\), and the reduced
branch has a representation of \(\pm2\) at primitive coordinates
of size at most \(D^A\) for a fixed \(A\). The form \(F\) is
even because \({\cal O}_E=\mathbb Z[\sqrt d]\). Dividing by \(2\)
and moving the known representation to \((1,0)\) by a unimodular
change gives an integral monic form of height \(D^{O_A(1)}\),
up to sign. Thus \(u-\alpha v\) is a unit in the totally real
cubic field for any other representation, where \(\alpha\) is
a root of the monic form.

For its conjugates \(\alpha_i\), put \(\beta_i=u-\alpha_i v\).
When \(|v|\) exceeds a fixed power of the form height \(H\),
the nearest root satisfies \(|\beta_i|\ll H^4/|v|^2\) by the
unit norm and the root spacing \(\gg H^{-2}\). For \(j,k\ne i\),
the nonzero quantity
\[
 \frac{\alpha_j-\alpha_i}{\alpha_k-\alpha_i}
 \frac{\beta_k}{\beta_j}-1
 =
 \frac{\alpha_j-\alpha_k}{\alpha_k-\alpha_i}
 \frac{\beta_i}{\beta_j}
\]
is \(O(H^7/|v|^3)\). The ratio on the left before subtracting
\(1\) is positive in the range where the error is small.
A reduced basis of the rank-two unit log lattice has product
of heights \(O(R_L)\), and fixed-degree Northcott gives an
absolute lower bound for its shortest vector. The unit exponents
are therefore \(O(\log(2|v|)+\log(2H))\). Matveev applied in
the totally real splitting field gives, for \(Y=\log(2|v|)\),
\[
 Y\ll (1+R_L)\log(2H)\log(2+Y+\log(2H)).
\]
The range of small \(|v|\) follows directly from the unit norm
since \(|u|\ll(1+H)|v|+2\). Inverting and using the regulator
bound gives \(\log(2x)\ll_A\sqrt D\,\log^4(2D)\), hence
\(D\gg_A(\log x)^2/(\log\log x)^8\) for large \(x\). The short
representation is a hypothesis on the branch.
[Matveev](https://www.mathnet.ru/links/5646f8921f46364ead641a660a874d86/im314_eng.pdf)

## Padé height ceilings inside an occupied branch

Fix \(0<\eta_0<1/2\). Suppose a primitive descent branch for \(d\ne1\)
has a point \((x_0,y_0)\) with \(x_0\ge \exp(C/\eta_0)\) and
\(0<|k_0|\le x_0^{1/2-\eta_0}\). Then any point in the same branch
and angular sector, with \(x\ge x_0\) and \(0<|k|\le\sqrt x\),
satisfies \(x\le x_0^{C'/\eta_0}\), for absolute effective \(C,C'\).
The real case has just the positive diagonal sector.

For the proof write \(a=y_0+m_0\sqrt d,\ b=\bar a\), and
\(\eta=a\xi^3\). The ideal of \(\xi\) is
\(\mathfrak a/\mathfrak a_{\rm seed}\), so \(x_0\xi\) is integral
because \(N\mathfrak a_{\rm seed}\mid x_0\).
Apply Voutier's Theorem 2.1 with degree \(3\).
For \(d<0\) use the imaginary quadratic base, \(\beta_1=1,\beta_2=-1\),
\(\gamma_1=-a,\gamma_2=-b\), evaluation point \(0\), and \(g=1\).
The target is \((\alpha-1)/(\alpha+1)\) with
\(\alpha=(a/b)^{1/3}\) on the principal branch.
For \(d>1\) use the rational base, \(\beta_1=\sqrt d,\beta_2=-\sqrt d\),
\(\gamma_1=-b,\gamma_2=-a\), evaluation point \(0\), and
\(g=d\sqrt d\). In that case use \(h_j=1\) for even \(j\) and
\(h_j=\sqrt d\) for odd \(j\), as \(h_j/g^j\in\mathbb Q\).
The target is \(\sqrt d(\alpha-1)/(\alpha+1)\) with
\(\alpha=(b/a)^{1/3}\). In the rational base the theorem takes
\(\tau=1\).

Put \(S=x_0^{3/2}+y_0\). The squared small and large square-root
differences of the normalized pair \(a,b\) are \(2|k_0|/S\)
and \(2S\). This follows from their real part or their product,
according to the sign. Voutier's parameters therefore have the form
\[
 E=\frac{S}{2J|k_0|},\qquad Q_{\rm P}=2JS,\qquad
 \kappa=\frac{\log Q_{\rm P}}{\log E},
\]
where \(s\) is Voutier's content parameter and
\(J=D_3/N_{s,3}\) is bounded above and below absolutely.
The theorem's prefactor is \(x_0^{O(1)}\) in this range.
For \(d<0\) approximate \(\alpha\) by \(\bar\xi/\xi\) in the
given angular sector. For \(d>1\) use \(\xi/\bar\xi\).
The corresponding Möbius transforms have integral numerator and
denominator in the required base and denominator
\(O(\sqrt{x_0x})\). In the real case multiplication of the
difference of the two integral conjugates by \(\sqrt d\) makes
the numerator an integer. Their errors are
\(O(\sqrt{|k|}/x^{3/2})\), with an additional factor \(\sqrt d\)
in the real case. Equality would give \(k=0\).
Absorbing \(|d|\le |k_0|\) into the seed-dependent factor, the
irrationality estimate gives
\[
 |k|\gg x_0^{-C_1}x^{2-\kappa}.
\]
The seed hypothesis gives \(\kappa\le3/2-c\eta_0\), which proves
the ceiling. Even \(|k_0|\le c_0\sqrt{x_0}\) for a sufficiently
small absolute \(c_0>0\) gives
\(3/2-\kappa\gg1/\log(2x_0)\) and hence
\(\log x\ll(\log(2x_0))^2\).

Combining the ceiling with the area determinant for two distinct
points in a branch bounds the number over all heights by
\(O((1+\log(1/\eta_0))h_3(d))\) per nonsplit squareclass in
\(0<|k|\le x^{1/2-\eta_0}\). In fact the determinant gives
\(x_{i+1}\gg |d|x_i^{3/2+\eta_0}\) in one angular sector, so the
large heights grow at least by the power \(4/3\). This bounds
the number up to \(\exp(C/\eta_0)\) and the number between the
first point above that threshold and its ceiling. Summing the
quadratic class-group first moment gives a bound linear in the
squareclass cutoff, uniformly over all \(x\). The ceiling starts
from an occupied branch. It leaves the least occupied height as
the pointwise problem. [Voutier, Theorem 2.1](https://arxiv.org/pdf/0805.3176#page=5) ·
[Bhargava–Shankar–Tsimerman](https://arxiv.org/abs/1005.0672)

## The split squareclass

For \(k=m^2\) let \(g=\gcd(y+m,y-m)\). Write
\((y+m)/g=aU^3\) and \((y-m)/g=bV^3\), with coprime cube-free
\(a,b\). Since \(g^2ab U^3V^3=x^3\), there is an integer
\(\delta>0\) with
\[
 g=ab\delta^3,\qquad x=ab\delta^2UV,\qquad
 m=\tfrac12ab\delta^3|aU^3-bV^3|.
\]
If the cube-free part of \(g\) is a power of \(2\), then
\(\{a,b\}\) is one of \(\{1,1\},\{1,2\},\{1,4\}\).
This includes every primitive point and all points with
\(2\)-primary \(\gcd(x,y)\). The bound for the irrationality
measure of \(\sqrt[3]2\) transfers to \(\sqrt[3]4\) by inversion.
For comparable \(U,V\) it gives
\(|aU^3-bV^3|\gg (UV)^{0.5675/2}\) after factoring the
difference of cubes. Widely separated \(U,V\), and \(a=b=1\),
follow by integer spacing. Therefore \(m^2\gg x^{0.5675}\)
with an effective absolute constant in this sector. For primitive
odd \(x\), the factors are coprime cubes and \(m^2\gg x^2\).
[Voutier, Corollary 2.2](https://www.numdam.org/item/JTNB_2007__19_1_263_0.pdf)


Voutier's Theorem 2.1 also treats some moving cube-free parts.
Suppose \(a>b\) after interchanging the factors and
\(a-b\le a^{1/6}/100\). Its Padé parameters satisfy
\(Q_{\rm P}\le 4e^{0.911}a\) and
\(E\ge e^{-0.911\cdot3^{3/2}}a/(a-b)^2>Q_{\rm P}^{2/3}\).
The irrationality exponent is at most \(5/2\) and the prefactor
in the denominator is \(O(a)\). Applying the bound to \(V/U\)
and using \(b\asymp a\) gives
\(|aU^3-bV^3|\gg (UV)^{1/4}\), with widely separated variables
again elementary. Thus \(m^2\gg (ab)^{3/2}\delta^5\sqrt{x}\)
uniformly in this nonprimitive sector.
[Voutier, Theorem 2.1](https://www.numdam.org/item/JTNB_2007__19_1_263_0.pdf#page=5)

## A resultant and divisibility bound near an integer square

In the Hall range, translation by an integer nearest \(-\sqrt x\)
gives a monic cubic of coefficient height \(\asymp\sqrt x\) with
two nearby roots and the third at distance \(\asymp\sqrt x\).
Thus \(|k|\asymp x^2\operatorname{sep}(P)^2\). This is the relevant
normalization in the Bugeaud–Mignotte equivalence.
[Bugeaud–Mignotte, Theorem 4](https://irma.math.unistra.fr/~bugeaud/travaux/PolSurvRev1.pdf)

For a direct divisibility argument put \(n\) nearest \(\sqrt x\),
\(t=x-n^2\), and \(h=2y-2n^3-3nt\). Also put \(s=4nh-3t^2\)
and \(g_2=\gcd(n,t^2)\). The normalized derivative resultant is
\(\operatorname{Res}(P,P'/3)=4k\). After translation by \(-n\),
the quadratic truncation is \(-3nT^2-3tT-h\), of discriminant
\(-3s\), and the resultant identity gives
\[
 8k=2n^2s+t^3+3ts+2h^2,\qquad g_2\mid s .
\]
Assume \(|t|\le n\) and \(|t|^3\le c g_2n^2\), with a
sufficiently small absolute \(c\). The second Taylor approximation
to \(x^{3/2}\) is
\[
 Y_2=n^3+\frac32 nt+\frac{3t^2}{8n}.
\]
It lies in \((g_2/8n)\mathbb Z\), and the error is
\(O(|t|^3/n^3)\). If \(y\ne Y_2\), this gives
\(|k|\gg g_2 n^2\). If \(y=Y_2\), then \(s=0\) and
\(k=t^3/8+9t^4/(64n^2)\). For \(k\ne0\), \(t\ne0\) and
\(4n\mid3t^2\), so \(|k|\gg n^{3/2}\). Thus the stated
sector has the uniform pointwise bound \(|k|\gg x^{3/4}\).
There is no primitivity condition.

The exponent \(3/4\) is attained in this sector. Take
\(n=b^2,\ t=4b,\ h=12\) for positive \(b\) coprime to \(6\).
Then \(x=b^4+4b,\ y=b^6+6b^3+6\), with \(k=8b^3+36\) and
\(\gcd(x,y)=1\). Primes dividing \(b\) do not divide \(y\),
and a common prime dividing \(b^3+4\) would divide \(2\).
The second Taylor expansion and its equality branch already
appear in Zagier. The useful extra parameter here is the
denominator reduction by \(\gcd(n,t^2)\).
[Zagier, Section 2](https://people.mpim-bonn.mpg.de/zagier/files/doi/10.2307/2007900/fulltext.pdf#page=8) ·
[Repository resultant valuation mechanism](../factorial/image-rigidity/SFC_CUBIC_UNIT_DEFORMATION_TRANSFER.md)

A coefficient-independent bound on the local resultant valuation
also fails for these cubics. Fix \(j\ge1\), take \(n\equiv1\bmod35\),
and put \(x=n^2,\ y=n^3+5^j z\). Choose \(z\) coprime to \(5n\)
and with \(y\equiv4\bmod7\) by the Chinese remainder theorem.
Then \(\gcd(x,y)=1\) and \(v_5(k)=j\). The cubic is irreducible
modulo \(7\), since \(T^3-3T\) never takes the value \(1\) there.
Thus even primitive irreducible Hall cubics can have arbitrarily
deep resultant divisibility at a fixed prime.

The remaining issue for the resultant route is to control the
cancellation when \(|t|^3\) exceeds the spacing supplied by
\(\gcd(n,t^2)n^2\). Local valuation or exact nonvanishing of a
resultant alone does not give an Archimedean lower bound of order
\(x^{1/2-\varepsilon}\).

## Squarefree gaps excluded by a unit in the real resolvent

Let \(T>1\) be squarefree, \(3\nmid T\), and \(k=-3T\). If an
integer point exists, it is primitive. Its cubic \(P\) is irreducible.
Indeed in the reducible identity \(4k=j^2(j-3x)\), the condition
\(3\mid k\) gives \(3\mid j\), and then \(27\mid k\). This contradicts
\(v_3(k)=1\). By the local calculations the cubic is totally
ramified at \(3\) with exponent \(4\) in its field discriminant.
Its quadratic resolvent is \(K=\mathbb Q(\sqrt T)\), unramified at
\(3\), and the conductor of the normal closure over \(K\) is
\(9{\cal O}_K\) if \(x\) is odd and \(18{\cal O}_K\) if \(x\) is
even. The even case requires \(T\equiv5\bmod8\) by reduction of the
Mordell equation modulo \(8\).

A necessary condition for the cubic is that its anti-invariant
order-three ray character be nontrivial on the trace-zero line
of \(1+3{\cal O}_K\) modulo \(9\). Indeed this group modulo
\(1+9{\cal O}_K\) is the additive space \({\cal O}_K/3\), and
the residue unit group has order prime to three. Inversion by
quadratic conjugation forces the character to vanish on the
invariant line. Exact conductor \(9\) at \(3\) forces the character
to be nonzero on the other line. The norm of any global unit is
\(\pm1\), so its image modulo local cubes lies in that other line.
A global unit noncube modulo \(9\) therefore rules out the
conductor when \(x\) is odd. The same exclusion works when \(x\)
is even if that unit is a cube modulo \(2\), since the dyadic
character then vanishes on it. This uses the ray class exact
sequence and the conductor formula for non-Galois cubic fields.
[Gaál–Pohst–Pohst, Lemmas 2–3](https://agag-fieker.math.rptu.de/~fieker/pub/gaal13.pdf) ·
[Cohen–Thorne](https://arxiv.org/abs/1301.3563) ·
[Breuer–Punch](https://doi.org/10.1080/10586458.2026.2633713)

This proves insolubility over \(\mathbb Z\) in two explicit families.

- If \(T=n^2+1\) is squarefree and \(n\not\equiv0,\pm2\bmod9\), use
  the unit \(n+\sqrt T\) of norm \(-1\). If \(x\) is even, then
  \(T\equiv5\bmod8\) and \(n\) is even. In the dyadic residue
  field \(\sqrt T\equiv1\) since \((1+\sqrt T)/2\) is integral.
  Thus the unit is \(1\) modulo \(2\).
- If \(T=n^2-1\) is squarefree and \(3\parallel n\), use the unit
  \(n+\sqrt T\) of norm \(1\). Squarefreeness forces \(n\) even
  and \(T\equiv3\) or \(7\bmod8\), so \(x\) is odd.

For the mod \(9\) condition, if \(n+\sqrt T\equiv(a+b\sqrt T)^3\)
then \(N(a+b\sqrt T)\equiv N(n+\sqrt T)\pmod3\). The constant
coefficient of the cube is \(4a^3-3aN(a+b\sqrt T)\). Its residues
are \(0,\pm2\) for norm \(-1\) and \(0,\pm1\) for norm \(1\).
The specified \(n\) fail these necessary conditions. At \(3\) the
coordinates \(a,b\) can be taken modulo \(9\) even when the ring
of integers has a basis with denominator \(2\).
Each family contains infinitely many \(T\). For \(n^2+1\) this
follows from the quadratic squarefree-values theorem in
\(n\equiv1\bmod18\). For \(n^2-1\), sieve the coprime odd linear
factors \(n-1,n+1\) in \(n\equiv6\bmod18\). Each \(p\ge5\) excludes
two residues modulo \(p^2\) and the tail for \(p>Z\) up to
\(n\le X\) costs \(O(X/Z+\sqrt X)\). [Ricci](https://doi.org/10.1007/BF03017586)

These families are everywhere locally soluble. At \(2\) choose
\(y=1\) for even \(k\) or \(y=0\) for odd \(k\), and take an odd
\(2\)-adic cube root of \(y^2-k\). At \(3\), \(x=1\) gives a
unit square root of \(1+k\). For \(p\ge5\) dividing \(k\),
\((1,1)\) is a nonsingular point on the reduction. At the other
primes the Hasse bound supplies an affine point on the smooth
reduction and Hensel lifts it. Real points exist for large \(x\).

### A density consequence of conductor \(9\)

Among negative squarefree residuals divisible by \(3\), the
lower relative density of insoluble residuals is at least \(5/9\).
For \(T\not\equiv5\bmod8\) in \(k=-3T\), the cubic field would
have exact conductor \(9\) over the real quadratic resolvent.
The mean number of such fields per real quadratic resolvent
unramified at \(3\) is \(1/3\). To obtain it, subtract the average
size \(4/3\) of the ordinary \(3\)-torsion from the mean \(2\)
for the odd part of the ray \(3\)-torsion modulo \(9\), then
divide by \(2\) for the two characters per field. The ordinary
average holds with \(3\) specified split or inert separately.

The mean \(1/3\) persists on fixing the quadratic \(2\)-adic
algebra. In the conductor-\(9\) field count the cubic algebras
at \(2\) are non-totally-ramified, with total mass \(3/2\).
The split quadratic resolvent has mass \(1/6+1/3=1/2\) from
cubic types \(111,3\), and the unramified quadratic field has
mass \(1/2\) from type \(21\). A ramified quadratic algebra
\(E_2\) has mass \(1/(2\operatorname{Disc}_2(E_2))\) from
\(\mathbb Q_2\times E_2\). After normalization these are precisely
the quadratic proportions. In each of these local types ordering
by \(T\) or by the fundamental discriminant differs by a fixed
scale. Thus the first moment rules out at least \(2/3\) of the
residuals in these types. The class \(T\equiv5\bmod8\) has
relative proportion \(1/6\) among squarefree \(T\) prime to \(3\).
This gives \((2/3)(5/6)=5/9\).
[Varma, Propositions 2.10, 4.4 and 4.5](https://arxiv.org/pdf/1609.02292) ·
[Bhargava–Taniguchi–Thorne, Theorem 1.3 and Section 8](https://arxiv.org/pdf/2107.12819)

For \(T\equiv5\bmod8\) either conductor \(9\) or conductor
\(18\) is possible. The conditional first moments there are
\(1/3\) and \(2/3\), respectively. For the second, Varma's
conductor-\(18\) average over resolvents unramified at \(2,3\)
gives one field in total on average for conductor dividing
\(18\). Subtract the averages \(1/6,1/6,1/3\) for exact
conductors \(1,2,9\). The \(1/6\) for exact conductor \(2\)
averages one half of \(1/3\) on the inert dyadic resolvents and
zero on the split ones by Proposition 4.5. Exact conductor
\(18\) also occurs only on the inert half, giving \(2/3\) there.
Thus the leading cubic count for the two candidate conductors
is the resolvent count in this residue class. Each cubic count
has a strictly negative secondary term of order \(X^{5/6}\) by
Bhargava–Taniguchi–Thorne. The sign comes from \(\zeta(1/3)<0\),
and their formula for the local secondary factors is positive
for these specifications. The quadratic squarefree progression
count has error \(O(\sqrt X)\). Subtracting the two field tallies
therefore rules out \(\gg X^{5/6}\) residuals of this remaining
class up to \(X\).


## Padé exclusion at the small root of a varying descent cubic

Fix a sufficiently small absolute \(c>0\). Let \(t\) be large
and let \(u\ne0\) be an integer with
\[
 (t,u)=1,\qquad 3\nmid u,\qquad |u|\le c t^{3/8},\qquad
 D=4t^3-u^2>0\text{ squarefree}.
\]
The branch considered here is the one generated by
\(\beta=(u+\sqrt{-D})/2\) in \(E=\mathbb Q(\sqrt{-D})\). The
squarefree equation \(y^2-x^3=-D\) in this branch has no point
corresponding to
\[
 s^3-3t sv^2-uv^3=2,\qquad |s|\le\sqrt t\,|v| .
\]
This is the region about the small real root of the descent
form. The other two real roots are outside this exclusion.

In fact \(N\beta=t^3\). The ideals of \(\beta\) and \(\bar\beta\)
are coprime since their sum element \(u\) is prime to \(t\).
Thus \((\beta)=\mathfrak a_0^3\) for
\(\mathfrak a_0=(t,\beta)\) of norm \(t\). The conjugate ideal
has integral basis \(t,\bar\beta\). For a point in this branch
write \(y+\sqrt{-D}=\beta\xi^3\) and
\(\xi=s+v\bar\beta/t\), since \(t\xi\in\bar{\mathfrak a}_0\).
Norm and imaginary part give the cubic equation above and
\[
 x=t s^2+usv+t^2v^2 .
\]
Conversely these formulas come from the norm and imaginary part
of the same expression. In the small root region \(v\ne0\), and
the cubic gives \(2t|s|v^2\le|u||v|^3+2\).
Therefore \(\xi/v\) has real part \(O((|u|+1)/t)\) and
imaginary part of order \(\sqrt t\).

Set \(a=-\beta,\ b=\bar\beta,\ q=t\xi,\ p=-\bar q\). Then
\(p/q\) is close to \(1\) and
\((p/q)^3=(a/b)(y-\sqrt{-D})/(y+\sqrt{-D})\). For the principal
cube root \(\alpha=(a/b)^{1/3}\) this gives
\[
 |\alpha-p/q|\ll t^{3/2}x^{-3/2},\qquad |q|=\sqrt{tx}.
\]
Let \(X_r(z)={}_2F_1(-r,-r-1/3;2/3;z)\), and let \(e_r\) clear
its denominators. Use its homogeneous Padé pair
\(P_r=e_rb^rX_r(a/b)\), \(Q_r=e_ra^rX_r(b/a)\).
The Padé remainder and denominator estimates give
\[
 |Q_r|\ll C^r t^{3r/2},\qquad
 |Q_r\alpha-P_r|\ll C^r |u|^{2r+1}t^{-3(r+1)/2}.
\]
These are the equal-modulus case of Voutier's estimates with
\(|a|=|b|=t^{3/2}\) and \(|a-b|=|u|\). All constants may be
chosen absolute.
[Voutier, Lemmas 7.1–7.4](https://arxiv.org/pdf/0805.3176#page=22)

For odd \(r\), conjugation interchanges \(P_r,Q_r\) with a minus
sign. Consequently \(\Delta_r=P_r q-Q_r p\) is a purely
imaginary algebraic integer. Its absolute value, if nonzero,
is at least \(\sqrt D/2\). It is nonzero for every odd \(r\).
Indeed the hypergeometric coefficient of \(z^j\) is
\[
 (-1)^j\binom rj\prod_{i=0}^{j-1}\frac{-3r-1+3i}{2+3i}.
\]
Thus \(X_r(z)\equiv(1-z)^r\bmod3\) and \(e_r\) is a \(3\)-adic
unit. In particular \(Q_r\equiv e_r(a-b)^r\bmod3{\cal O}_E\)
is a rational integer unit modulo \(3\), so \(Q_r^3\equiv
\pm1\bmod9{\cal O}_E\). If \(\Delta_r=0\), conjugation forces
\(q/Q_r=\lambda\in\mathbb Q\), and the \(\sqrt{-D}\)
coefficient \(m_r\) of \(\beta Q_r^3\) satisfies
\(1=(\lambda/t)^3m_r\). But \(2m_r\equiv\pm1\bmod9\).
To see the latter, subtract the conjugate of \(\beta Q_r^3\)
and use that differences of integral conjugates are integral
multiples of \(\sqrt{-D}\). Hence \(m_r\equiv\pm4\bmod9\),
incompatible with a rational cube of \(3\)-adic valuation zero.

On the other hand,
\[
 |\Delta_r|\ll C^r
 \left(\frac{t^{3r/2+2}}x+
       |u|^{2r+1}\sqrt x\,t^{-3r/2-1}\right).
\]
The right-hand side is \(<\sqrt D/2\) throughout intervals with
lower endpoint \(C_0 C^r t^{3r/2+1/2}\) and upper endpoint
\(C_0^{-2}C^{-2r}|u|^{-4r-2}t^{3r+5}\), after increasing
absolute constants. Successive odd intervals overlap.
The ratio of the upper endpoint for \(r\) to the lower endpoint
for \(r+2\) is at least a constant times
\(C^{-3r}c^{-4r-2}t^{3/4}\), which is \(>1\) for small enough
\(c\) and large \(t\). The odd intervals start at \(O(t^2)\)
and their upper endpoints tend to infinity. Thus \(x\ll t^2\).
The norm expression now gives \(|v|=O(1)\), and
\(2t|s|v^2\le |u||v|^3+2\) forces \(s=0\) for large \(t\).
Then \(-uv^3=2\) is impossible since \(u\) is odd.

For example \(u=1\) and \(t\) a multiple of \(6\) gives
infinitely many squarefree \(D\) by the cubic squarefree-values
theorem. The ideal class of \(\mathfrak a_0\) is nontrivial
for large \(t\). If it were principal, its generator would
have absolute value \(\sqrt t\) but any nonreal quadratic
integer has imaginary part of size at least \(\sqrt D/2\).
A real generator would make \((\beta)\) a rational ideal,
contradicting coprimality with its conjugate. These branches
survive the immediate mod \(2\) and mod \(9\) conditions on
the cubic equation. The result thus excludes the small root
in genuinely nonprincipal branches.
[Hooley](https://doi.org/10.1515/crll.1968.229.147) ·
[Bombieri–van der Poorten–Vaaler, Section 9](https://www.numdam.org/item/ASNSP_1996_4_23_2_211_0/)

## Pointwise bounds for sparse normalized index forms

Let \(k\) be squarefree and suppose the maximal-order index form
in some unimodular coordinates is
\[
 J(U,V)=AU^3+BU^2V+CUV^2+EV^3,\quad A\ne0 .
\]
Put \(a=3(3AC-B^2)>0\),
\(b_0=2B^3-9ABC+27A^2E\), and \(b=|b_0|>0\).
If \(a\ge360b^4\) and \(27A^2 I\le a+b+1\), then the
uniform bound \(|k|\gg x^{3/4}\) holds in this sector. The
sufficient condition \(a\ge162A^2\) covers all four levels
\(I=1,2,3,6\).

Indeed \(W=3AU+BV\) gives
\(G(W,V)=W^3+aWV^2+b_0V^3=27A^2J(U,V)\). Also
\[
 \operatorname{disc}(G)=729A^2D_L,\qquad
 H_G(W,V)=729A^2x.
\]
Write \(h=\gcd(W,V)\). Since \(\gcd(U,V)=1\), one has \(h\mid3A\).
Apply Wakabayashi's classification to \((W/h,V/h)\).
The coordinate and mixed directions have \(x\ll a^2\) by
the Hessian and \(h\mid3A\), while \(|k|\gg a^3/A^2\gg a^2\).
The exceptional direction is \(W/V=-b_0/a\) and the value
\(|G(W,V)|=|W|^3\) forces
\[
 |A|=I s^3,\qquad |W|=3I s^2,\qquad
 |U|=\left|\frac{27A^2E-B^3}{s b_0}\right|.
\]
If \(e=v_p(s)>0\) and \(v_p(B)=\ell<2e\), the numerator
has valuation \(3\ell\). Then \(v_p(b_0)\le3\ell-e\).
The only term in \(b_0\) that could have such a low valuation
is \(-9ABC\), whose valuation is at least \(3e+\ell\).
This is impossible. Thus \(v_p(B)\ge2e\), giving
\(v_p(a)\ge3e+v_p(3)\) and \(v_p(b_0)\ge5e\).

For \(p\ge5\) the identity
\(4a^3+27b_0^2=729A^2|D_L|\) now forces \(v_p(k)\ge3e\).
At \(2\) the numerator has valuation at least \(9e+1\),
which forces \(v_2(k)\ge3e-1\). Both contradict squarefreeness.
At \(3\), the bound on \(B\) gives \(v_3(a)\ge3e+2\) and
\(v_3(b_0)\ge5e+1\). The numerator has valuation at least
\(9e+6\), giving \(v_3(k)\ge3e-3\) and thus \(e\le1\).
Hence \(s\mid3\). In particular \(A\) is absolutely bounded
on the exceptional direction. There \(x\ll a^4\) by the
Hessian, while \(|k|\gg a^3\) by the discriminant.
[Wakabayashi, Theorem 3](https://doi.org/10.1016/S0022-314X%2802%2900010-0)

For the directly depressed monic form \(U^3+aUV^2-bV^3\) with
\(a\ge360b^4\), the squarefree points can be classified exactly.
The axis \((1,0)\) has negative Hessian and \((0,1)\) would
give \(x=a^2/9,\ k=b^2(a/3)^3+b^4/4\). Integrality forces
\(3\mid a\) and \(2\mid b\), contradicting squarefreeness.
The \((\pm1,1)\) values exceed \(6\). The exceptional direction
has index \((b/\gcd(a,b))^3\), so \(b\mid a\) and \(I=1\).
Since \(108k=b^2(4b(a/b)^3+27)\), squarefreeness forces
\(b=2\) or \(6\). The two surviving families are
\[
\begin{aligned}
 (x,y,k)&=(36t^4+4t,\ 216t^6+36t^3+1,\ 8t^3+1),\\
 (x,y,k)&=(4t^4+4t,\ 8t^6+12t^3+3,\ 8t^3+9),
\end{aligned}
\]
at their squarefree values and in the stated parameter range.
They come from \(a=6t\). The corresponding defining cubics
are Eisenstein at \(2\). The \(b=6\) cubic is Eisenstein at \(3\).
For \(b=2\), translating by \(-1\) gives an Eisenstein polynomial
at \(3\) since squarefreeness excludes \(t\equiv1\bmod3\).
At primes \(p\ge5\) the discriminant exponent is at most one.
Thus the displayed base orders are maximal. The squarefree
cubic-values theorem supplies infinitely many values in each
family, so \(3/4\) is attained in this sector. [Hooley](https://doi.org/10.1515/crll.1968.229.147)


## Finite cross-check

A PARI check on \(1\le x\le100\), \(1\le y\le130\) compared the
local formulas with the field discriminant for all \(7{,}750\)
primitive pairs with nonzero residual and irreducible \(P\).
All checks passed. Every predicted conductor in
\(\{1,2,3,6,9,18\}\) and every predicted squarefree-case index
in \(\{1,2,3,6\}\) occurred. This checks the local case split on
a finite range. The proofs above establish its unrestricted scope.

The continued squarefree checks used separate exact inputs. An exhaustive
mod \(9\) cube table gives \(n\equiv0,\pm2\bmod9\) in the
\(n^2+1\) unit family and \(n\equiv0\bmod9\) among \(3\mid n\)
in the \(n^2-1\) family. The hypergeometric coefficient congruence
was checked through degree \(60\). In the two Wakabayashi families,
PARI found the asserted maximal discriminant at all \(340\) and
\(376\) squarefree values, respectively, for \(1\le t\le600\).
An integer scan of the \(136\) explicit forbidden squarefree
residuals with \(2\le n\le200\) found no point for \(1\le x\le5000\).
These are finite cross-checks of the proofs.

