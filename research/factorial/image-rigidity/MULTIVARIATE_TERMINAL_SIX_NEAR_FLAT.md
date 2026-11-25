# Explicit near-flat separation at terminal six for shared support

Suppose all unequal exponent rows of an ordinary two-monomial profile
are shared, so both exponents are positive. Let \(u\) be the
order-five normalized outer geometric mean, after dividing out the
binomial coefficients and normalizing both endpoint coefficients to
one. If \(u\ge24/25\), the order-five and order-six polynomials
have disjoint root arguments after independent positive rescalings.
The profile may have arbitrary dimension.

## Scaling and the order ratio

Use the quintic and sextic notation from the
[terminal-six coefficient links](MULTIVARIATE_TERMINAL_SIX_COEFFICIENT_LINKS.md).
The normalized multiplier arrays are
\[
 (1,u e^\alpha,v e^\chi,v e^{-\chi},u e^{-\alpha},1),\qquad
 (1,w e^\beta,y e^\eta,z,y e^{-\eta},w e^{-\beta},1).
\]
Put \(p=1-u=\lambda^2\), \(q=1-w=r p\), so \(0<\lambda\le1/5\).
The strict inequality \(p>0\) holds for distinct profiles. We use
\[
 519/500\le r\le21/20 .
 \tag{1}
\]

For the lower bound, the [outer-mean comparison](MULTIVARIATE_LINKED_GAMMA_ANGLE_COMPACTIFICATION.md)
gives \(w\le u(1+(1-u^2)/24)^{-1/2}\).
For \(t=(1-u^2)/24\le1/300\),
\(1-(1+t)^{-1/2}\ge t/[2(1+t)^{3/2}]\).
Therefore \(r\ge1+u(1+u)/[48(1+t)^{3/2}]\ge519/500\),
using \(u\ge24/25\).

For the upper bound, one Gauss factor of an ordinary shared row
with smaller exponent \(a\) and gap \(d\) has reciprocal
\(x_5=d/(5a+h)\), \(1\le h\le d\). Every such \(x_5\) is at most
\((u^{-2}-1)/4\le49/2304\). As \(x_5\ge d/(5a+d)\),
we get \(d/a\le245/2255<1/8\).
Hence the factorwise ratio \(5x_6/(4x_5)\) is no larger than
\(\frac54(5+1/8)/(6+1/8)=205/196<21/20\).
The product comparison gives \(w\ge u^{205/196}\).
The tangent inequality \(u^\gamma\ge1-\gamma(1-u)\) for
\(\gamma\ge1\) yields the bound on \(r\).

## Stable Cayley coefficients

Make the substitution
\[
 t=\frac{\lambda X-1}{\lambda X+1}.
\]
For order \(j\), multiply its normalized polynomial by
\((1+\lambda X)^j/(2^j\lambda^j)\). For the quintic set
\[
 a=u\cosh\alpha,\quad b=v\cosh\chi,\quad
 c=u\sinh\alpha,\quad d=v\sinh\chi .
\]
Before division by \(2^j\lambda^j\), and in the variable
\(Z=\lambda X\), its coefficients from degree five to zero are
\[
 2(5a+10b+1),\quad10(3c+2d),\quad20(a-2b+1),\quad
 -20(c+2d),\quad-10(3a-2b-1),\quad-10(c-2d).
 \tag{2}
\]
For the sextic set \(a=w\cosh\beta,\ b=y\cosh\eta,\ e=z\),
\(c=w\sinh\beta,\ d=y\sinh\eta\). The coefficients from
degree six to zero are
\[
\begin{gathered}
 2(6a+15b+10e+1),\quad 12(4c+5d),\quad
 30(2a-b-2e+1),\quad -120d,\\
 -30(2a+b-2e-1),\quad-12(4c-5d),\quad
 -2(6a-15b+10e-1).
\end{gathered}                                                       \tag{3}
\]
Both expressions follow directly from the binomial expansion of
\((Z-1)^k(Z+1)^{j-k}\) in the Cayley transform.

Below, decimals are exact rationals. The parameter bounds
\(p\le.04,\ q\le.05\) imply
\[
 x_5^*:=(u^{-2}-1)/4\le .55p,\qquad
 x_6^*:=(w^{-2}-1)/5\le .45q .
\]
The sum of all Gauss reciprocals at each order is at most
its corresponding starred quantity, by expanding the product
that gives the outer geometric mean. The
[quintic skew estimates](MULTIVARIATE_TERMINAL_FIVE_ANGULAR_SEPARATION.md)
give \(|\alpha|\le .31p^2\),
\(|\alpha-2\chi|\le .055p^4\), and hence \(|\chi|\le .16p^2\).
For the sextic, the local outer skew is bounded by
\(5x_6^2/3\): write its numerator as
\((g_1+g_4-g_5)+(g_2+g_3-g_5)\) and bound their logarithms by
\(4x_6^2\) and \(6x_6^2\). The
[sextic link](MULTIVARIATE_TERMINAL_SIX_COEFFICIENT_LINKS.md#within-order-bounds)
then gives
\[
 |\beta|\le .34q^2,\quad |4\beta-5\eta|\le .25q^4,\quad
 |\eta|\le .274q^2 .
\]
For all these small skews we use
\(|\sinh t|\le1.001|t|\), \(|\sinh t-t|\le|t|^3\), and
\(\cosh t-1\le .501t^2\). These follow by bounding the tails
of the defining series geometrically on \(|t|\le.01\).

Here are all the quantities substituted into (2) and (3)
to bound their scaled even coefficients. Write
\[
\begin{aligned}
 a_5&=1-p+A_5p^4,&
 b_5&=1-\tfrac32p+R_5p^2+B_5p^4,\\
 a_6&=1-q+A_6q^4,&
 b_6&=1-\tfrac85q+\tfrac{12}{25}q^2+(R_y+H)q^3+B_6q^4,\\
 e_6&=1-\tfrac95q+\tfrac{18}{25}q^2+(R_z+\tfrac32 H)q^3+Dq^4 .
\end{aligned}
\]
Their enclosure is
\[
\begin{array}{c|c@{\qquad}c|c}
 A_5&[0,.05]& B_5&[0,.013]\\
 R_5&[3/8,3/8+.26p]&H&[0,.235]\\
 A_6&[0,.059]&B_6&[0,.038]\\
 R_y&[.064,.064+.03q]&R_z&[.048,.048+.02q]\\
 D&[-.078,.078]&&
\end{array}                                                        \tag{4}
\]
The \(A,B\) bounds follow from the hyperbolic estimates above.
The geometric-mean defects \(V=\log(v/u^{3/2})\) and
\(Y=\log(y/w^{8/5})\) come from sums of
\(G_a(s)=as-\operatorname{artanh}(a\tanh s)\) with \(a=1/2,3/5\).
As \(G_a'(s)\le a s^2\), we get
\(V\le [-\log(1-p)]^3/6\) and \(Y\le [-\log(1-q)]^3/5\).
Since \(-\log(1-h)\le h/(1-h)\) and \(e^h\le(1-h)^{-1}\)
for \(0\le h<1\), this implies
\[
 0\le v-u^{3/2}\le .19p^3,\qquad
 0\le \delta_y:=y-w^{8/5}\le .235q^3 .
\]
The Taylor coefficients of \(u^{3/2},w^{8/5},w^{9/5}\)
give the stated \(R_5,R_y,R_z\). For \(R_5\) the pure cubic
remainder is bounded by \(.07p^3\); for \(R_y,R_z\) the next
remainders are bounded by \(.03q^4,.02q^4\), respectively,
by the corresponding derivatives on \(u\ge.96,w\ge.95\).

For the last entry in (4), let
\[
 J=-\log z+\tfrac32\log y-\tfrac35\log w,\qquad
 z=w^{-3/5}y^{3/2}e^{-J}.
\]
The link gives \(0\le J\le(6/5)(x_6^*)^5\le .023q^5\).
Expanding \(w^{-3/5}y^{3/2}\) at \(y=w^{8/5}\)
yields the linear term \(\frac32 w^{1/5}\delta_y\),
with second-order error at most \(\frac38 w^{-7/5}\delta_y^2\).
For \(w\ge.95\), use \(|w^{1/5}-1|\le.21q\),
\(w^{-7/5}\le1.1\), \(w^{-3/5}\le1.04\).
Together with \(1-e^{-J}\le J\) and \(y\le1\), these give
\[
 \left|z-w^{9/5}-\tfrac32\delta_y\right|\le .078q^4.
\]
This proves the bound on \(D\). The coefficient calculations
using (4) preserve cancellations in the first-order terms
and in the central sextic coefficient.

For the scaled odd-part perturbations, the checker uses the
following absolute bounds on the coefficients at the indicated
degrees. They are the coefficients with parity opposite to
the order of the polynomial.
\[
\begin{array}{c|ccc}
 j=5,\ \text{degree}&4&2&0\\
 \text{bound}&.4\lambda^3&.4\lambda&.052\lambda\\ \hline
 j=6,\ \text{degree}&5&3&1\\
 \text{bound}&.515r^2\lambda^3&.516r^2\lambda&.158r^3\lambda .
\end{array}                                                        \tag{6}
\]
The entries except the constant of the quintic and the linear
term of the sextic follow by substituting the separate skew
bounds into (2) and (3). For the quintic constant, write
\(c-2d=u(\sinh\alpha-2\sinh\chi)+2(u-v)\sinh\chi\).
We have \(0\le u-v\le p/2\) from the shape bounds and the
tangent to \(u^{3/2}\). The linked skew and the hyperbolic cubic
error then bound this difference by
\(.055p^4+.04p^6+.161p^3\).
For the sextic linear term write
\(4c-5d=w(4\sinh\beta-5\sinh\eta)+5(w-y)\sinh\eta\).
Here \(0\le w-y\le3q/5\), and the separate skew and cubic
error bounds give \(.25q^4+.27q^6+.825q^3\).
Division by the indicated scaling in (2) and (3) gives (6).

## Exact ray exclusions and the flat anchor

The [checker](verify_multivariate_terminal_six_near_flat.py)
substitutes (4) and (6) into (2) and (3) by rational interval
arithmetic. Leading-term domination on \(|X|=2\) shows that
all Cayley roots have \(|X|<2\). In this disk the original
argument is near the negative real direction, since
\[
 t=\frac{\lambda^2|X|^2-1+2i\lambda\,\Im X}{|1+\lambda X|^2}.
\]
For a positive ray parameter \(b\), consider the original
argument \(\theta=\pi-2\arctan(\lambda b)\).
Writing \(X=x+iY\), the ray equation is exactly
\[
 (Y-b)(1+\lambda^2 bY)=-b\lambda^2 x^2,\qquad 0\le Y\le b .
 \tag{7}
\]
The checker starts from this interval for \(Y\), then
encloses it by two applications of (7). It is enough to
consider \(0\le x\le2\): applying \(X\mapsto-\overline X\)
to a root also changes the sign of the opposite-parity
coefficients, and their enclosing intervals in (6) are
symmetric about zero.

The four ray parameters, in increasing original angle, are
\(527/500,21/25,13/25,33/100\). In every rational box for
\(\lambda,r,x\), the checker proves a strictly signed
projection of the form \(\Im(P_j(X)\overline X^{\,h})\)
for some integer \(h\in[0,j]\), using the scaled polynomial
\(P_j\). It uses
\(\Im(X^k\overline X^{\,h})
 =|X|^{2h}\Im(X^{k-h})\) for \(k\ge h\), and the conjugate
identity for \(k<h\). The first replay used respectively
\(203,43,273,43\) boxes at order five and
\(20515,4527,12933,1743\) at order six.

At \(\lambda=0\), the scaled limits given by (2)--(4) are
\[
 X^5+\tfrac54 X^3+\tfrac{15}{64}X,\qquad
 X^6+\tfrac32 rX^4+\tfrac9{20}r^2X^2+\tfrac3{200}r^3 .
\]
The positive squared imaginary coordinates of the quintic
are \((5\pm\sqrt{10})/8\). Those of the sextic are \(r/4\)
times the three real roots of
\[
 A^3-6A^2+\tfrac{36}{5}A-\tfrac{24}{25}.
\]
The latter lie in \((.152,.153),(1.427,1.428),(4.420,4.421)\)
by the signs at these six rational endpoints.
For every \(r\) in (1) the imaginary coordinates alternate in
decreasing order around the four ray parameters. Thus the
five upper sectors contain sextic, quintic, sextic, quintic
and sextic roots, respectively, at the anchor.

For a given actual profile, multiply every shared Gauss
reciprocal at both orders by the same parameter decreasing
from one to zero. The cross-order factor ratios remain
fixed, and the derivation of (1) still applies: its lower
bound uses the generic linked-factor comparison and its
upper bound uses only the unchanged factor ratios.
Therefore the scaled arrays stay in the verified boxes.
For positive parameter the transformations are continuous,
and at zero the scaled coefficients tend to the displayed
limit with \(r\) given by the ratio of the first derivatives
of \(1-w\) and \(1-u\) along the path. The angular rays are
represented in the scaled plane by (7), which becomes
\(Y=b\) at the anchor. The root counts are consequently
preserved. In the last sector about the negative real axis
the quintic's one root is real by conjugation, and the
actual sextic is strictly positive on the real axis.
The arguments are disjoint, as claimed.
