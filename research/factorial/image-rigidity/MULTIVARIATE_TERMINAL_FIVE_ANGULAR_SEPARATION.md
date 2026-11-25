# Angular separation at terminal five for arbitrary exponent rows

For distinct nonnegative integer vectors \(a,b\), let
\(Q_n(t)=\Lambda((X^a+tX^b)^n)\), where
\(\Lambda(X^c)=\prod_j c_j!\). Then \(Q_4\) and \(Q_5\) are coprime.
The result applies over every characteristic-zero field. The angular
argument also applies when the factorial functional is replaced by
independent Gamma variables with arbitrary fixed positive shapes.

The proof uses three rational rays in the upper half-plane. In ascending
angular order, the first sector contains one quintic root, the second one
quartic root, the third one quintic root, and the last one quartic root.
The quintic's remaining root is negative real. Thus the quartic and quintic
arguments are disjoint after positive rescaling of either polynomial.
The middle and near-flat ranges use exact rational inputs and outward
interval arithmetic in the [replay](verify_multivariate_terminal_five.py).
Run the replay with Python and python-flint. It verifies a contiguous
cover of 49 slabs on \([1/500,13/20]\), 30 slabs on \([13/20,97/100]\),
and the six near-flat ray exclusions on \([.9676,1)\).

## The normalized coefficient ranges

Divide each moment polynomial by its constant term and positively rescale
its argument so that its leading coefficient is one. The resulting
polynomials are
\[
\begin{aligned}
 P_4(z)&=1+4u e^\alpha z+6v z^2+4u e^{-\alpha}z^3+z^4,\\
 P_5(z)&=1+5w e^\beta z+10y e^\chi z^2+10y e^{-\chi}z^3
                     +5w e^{-\beta}z^4+z^5.
\end{aligned}
\]
Each distinct exponent row contributes its own normalized coefficient
factors. In a row with increasing exponents \(a_0,a_0+d\) and shape
\(\gamma>0\), Gauss multiplication reduces these contributions to factors
\[
 F_{n,k}=\frac{(T)_k}{(T)_n^{k/n}},
 \qquad T=\frac{\gamma+n a_0+j}{d},\quad 0\le j<d.
\]
A decreasing row reverses the array. Equal rows cancel. For \(n=4\) and
\(5\), the parameters in one factor satisfy \(T\le T'\le5T/4\).
Write \(x=1/T\) and \(g_j=\log(1+jx)\), using the appropriate \(T\) at
each order. The squared local outer geometric means are \(1/(1+3x)\)
and \(1/(1+4x)\). Dividing the local inner mean by the outer one gives
\(\sqrt{(1+x)/(1+2x)}\) at order four and
\(\sqrt{(1+x)/(1+3x)}\) at order five. Multiplying the local inequalities
gives
\[
 u^{4/3}\le v\le u\sqrt{\frac{1+2u^2}{2+u^2}},\qquad
 w^{3/2}\le y\le w\sqrt{\frac{1+3w^2}{3+w^2}},\qquad
 u^{4/3}\le w\le \sqrt{\frac{15u^2}{16-u^2}}<u<1.                 \tag{1}
\]
For the first two bounds use
\(F_c(s)=((1+c)s+1-c)/((1-c)s+1+c)\) with \(c=1/3,1/2\).
For \(0<c<1\) and \(0<s,t\le1\), one has
\(s^c\le F_c(s)\) and \(F_c(s)F_c(t)\le F_c(st)\) by the logarithmic
derivative and the Cayley coordinate \((1-s)/(1+s)\), respectively.
For the cross-order lower bound use
\(1+4/T'\le(1+3/T)^{4/3}\). For the upper bound,
\[
 \frac{1}{1+4/T'}\le\frac{15s}{16-s},\qquad s=\frac{1}{1+3/T}.
\]
For \(f(s)=s/(\nu+(1-\nu)s)\), where \(\nu>1\), the product of the
denominators at \(s,t\) minus the denominator at \(st\) equals
\(\nu(\nu-1)(1-s)(1-t)\). Hence \(f(s)f(t)\le f(st)\).

Set \(\kappa=\beta-2\chi\). With \(x=(u^{-2}-1)/3\) for the quartic
and \(x=(w^{-2}-1)/4\) for the quintic, the skew bounds used in the
middle range are
\[
\begin{aligned}
 e^{|\alpha|} &\le \left(\frac{(1+x)(1+2x)}{1+3x}\right)^{1/4},\\
 e^{|\beta|} &\le
  \left(\frac{((1+x)(1+2x)(1+3x))^2}{(1+4x)^3}\right)^{1/10},\\
 e^{|\chi|} &\le
  \left(\frac{(1+2x)^4}{(1+x)(1+3x)(1+4x)}\right)^{1/10},\\
 e^{|\kappa|} &\le
  \left(\frac{((1+3x)(1+x))^4}{(1+2x)^6(1+4x)}\right)^{1/10}.
\end{aligned}                                                       \tag{2}
\]
For a single increasing factor, the signed logs are
\(\alpha=(g_3-g_1-g_2)/4\),
\(\beta=(3g_4-2g_1-2g_2-2g_3)/10\),
\(\chi=(g_1+g_3+g_4-4g_2)/10\), and
\(\kappa=(g_4-4g_3+6g_2-4g_1)/10\). For a decreasing
factor all four signs reverse. To sum the bounds across factors,
put \(s=\log(1+(n-1)x)\), \(q=e^s\), and
\(g_j(s)=\log((n-1-j+jq)/(n-1))\). Each log on the right of (2) is
convex, vanishes at zero, and has derivative zero there. For the
first two logs convexity follows by summing the terms
\(g_j''=j(n-1-j)q/(n-1-j+jq)^2\) with positive coefficients.
For the last two, the second derivatives of ten times the logs are
\[
 \frac{2q(3q^4+72q^3+170q^2+72q+3)}
 {(q+1)^2(q+3)^2(3q+1)^2},\qquad
 \frac{6q(q-1)^2(11q^2+26q+11)}
 {(q+1)^2(q+3)^2(3q+1)^2}.
\]
A convex function vanishing at zero is superadditive on the nonnegative
axis. This proves (2), since the factorwise \(s\)'s sum to
\(-2\log u\) or \(-2\log w\).

The near-flat range uses the following simpler estimates, with the
same conventions for \(x\):
\[
 |\alpha|\le x^2/2,\qquad |\beta|\le x^2,\qquad
 |\kappa|\le3x^4/5.                                                   \tag{3}
\]
For the quartic use the integral of \(-g''\) over the rectangle
\([0,1]\times[0,2]\) to bound \(g_1+g_2-g_3\) by \(2x^2\).
For the quintic use
\(g_j-(j/4)g_4\le j(4-j)x^2/2\), which follows by comparison
with a quadratic with second derivative \(-x^2\). Finally,
the fourth difference of \(g\) is an integral of \(g^{(4)}\)
over the unit four-dimensional cube and
\(|g^{(4)}|\le6x^4\). These give the three single-factor estimates.
For their sums use
\(\sum x_i\le(\prod_i(1+(n-1)x_i)-1)/(n-1)\) and
\(\sum x_i^k\le(\sum x_i)^k\).

## The small outer mean

With \(T_*=1+(n-1)x\), one has
\(F_{n,k}\le T_*^{-\min(k,n-k)/n}\).
To see this, take logarithms, set \(G_j=\log(1+jx)\), and split
\(\sum_{j=1}^{n-1}G_j\) before \(k\). The terms with indices below
\(k\) are no larger than the later terms. When \(k\le n/2\), compare
\((n-k)\sum_{j<k}G_j\) with \(k\sum_{j=k}^{n-2}G_j\), using
\((n-k)(k-1)\le k(n-1-k)\). When \(k>n/2\), add
\((2k-n)G_{n-1}\) to the right; now the total weights agree.
Reversal preserves the bound.

If \(u\le1/500\), apply this bound to both coefficient arrays.
The sum of normalized interior coefficients is at most
\(4(\sqrt u+u^{3/2})+6u<.193\) at order four and
\(5(w^{2/5}+w^{8/5})+10(w^{4/5}+w^{6/5})<.497\) at order five.
These bounds use (1) and the pair products \(u^2,w^2,y^2\), since
if \(A,B\le M\) and \(AB\le p\) with \(p\le M^2\), then
\(A+B\le M+p/M\). For the rational comparisons one can use
\(\sqrt u<.045\) and \(w^{2/5}<.084\).

If that sum is \(\epsilon<1\), any root has angular distance at
most \(\arcsin(\epsilon)/n\) from an odd multiple of \(\pi/n\).
For a root of modulus at most one, take the argument of the
disk \(|1+z^n|\le\epsilon\); for larger modulus apply the same
argument to the reciprocal polynomial. These angular discs
can be followed from \(1+z^n\) by multiplying the interior
coefficients by a scalar in \([0,1]\), so each contains one
root. The closest odd multiples at orders four and five
differ by \(\pi/20\). The combined errors are less than
\(.154<\pi/20\), by
\(\arcsin x\le x+x^3/(6(1-x^2))\). To obtain this upper bound,
integrate \((1-t^2)^{-1/2}\le1+t^2/(2(1-x^2))\) for \(0\le t\le x\).

## Exact ray checks on compact intervals

On a ray of angle \(\theta\in(0,\pi)\), put \(Y=\cos\theta\) and
\(R=|z|\). The bounds \(|z|\le16,32\), respectively, follow from
the root equation because the endpoint-normalized coefficients are at
most the binomial coefficients. Applying the same argument to the
reciprocals gives \(R\ge1/16,1/32\). These bounds apply throughout the
coefficient range (1)–(2) if one also retains that elementary
coefficient cap.

For the quartic, write \(A=e^\alpha\). Reversal preserves its set
of upper-half-plane arguments and permits \(A\ge1\). The imaginary
part of \(P_4(z)/z^2\), divided by \(\sin\theta\), gives
\[
 R/A-A/R=D,\qquad D=-Y(R^2-R^{-2})/(2u).
\]
It determines \(A=R(\sqrt{4+D^2}-D)/2\). The real part then gives
\[
 (R^2+R^{-2})(2Y^2-1)+4uY\sqrt{4+D^2}+6v=0.                \tag{4}
\]
The replay excludes the pair of necessary conditions by enclosing
both sides on dyadically subdivided boxes in \(u,R\).

For the quintic, put \(B=e^\beta\ge1\) after reversal and \(C=e^\chi\).
Set \(E=8Y^2-2\) and \(L=4Y\). The real and imaginary parts of
\(P_5(z)/z^{5/2}\), divided by \(2\cos(\theta/2)\) and
\(2\sin(\theta/2)\), give, on adding and subtracting and then
multiplying by \(R^{5/2}\),
\[
\begin{aligned}
 P&=ER^5+20Y(w/B)R^4-10wBR-L=-20yR^3/C,\\
 M&=-LR^5-10(w/B)R^4+20YwBR+E=-20yR^2C.
\end{aligned}                                                       \tag{5}
\]
In particular \(R M-C^2P=0\) and \(PM=400y^2R^5\). Use (1)–(2)
to enclose \(y\) and \(C^2\), along with \(C^2/B=e^{-\kappa}\).
This excludes a root when any necessary equality in (5) fails.
It is enough to subdivide the three variables \(w,R,B\).
An equivalent enclosure of \(RM-C^2P\) uses
\((RM-BP)+B(1-C^2/B)P\), with the first term expanded at \(R=1\).

At each rational slab of \(u\), use an exact anchor with
\(u\) at the midpoint, \(w=u^{16/15}\), \(v=u^{4/3}\),
\(y=w^{3/2}\), and zero skews. These values satisfy (1).
After substituting \(X=z+z^{-1}\) and removing \(z+1\) from
the quintic, the two quadratics are
\[
 X^2+4uX+6v-2,\qquad X^2+(5w-1)X+10y-5w-1.
\]
Their signs at the three values \(X=2Y\) give the four
sector counts. The endpoint signs at \(X=-2,2\) locate the roots
in \((-2,2)\). The certificate checks the signs at \(-2\)
and at \(2Y\); the signs at \(2\) are immediate from the formulas.

The coefficient ranges are connected to this anchor. At fixed
\(u,w,v,y\), multiply the skews by a common scalar in \([0,1]\).
This preserves (2), the cap by the binomial coefficients, and
the corresponding symmetric conditions. Then move \(v,y\) to
their lower bounds, move \(w\) continuously to \(u^{16/15}\)
keeping \(y=w^{3/2}\), and finally vary \(u\) to the slab
midpoint. For the cap during skew scaling, each geometric
mean is at most one and the absolute value of each skew
decreases. This path allows the sector root counts to be
continued because no root can meet a separating ray, or
meet zero. Count with multiplicity in the full complex
sectors bounded by the rays and their reflections. In particular,
the last sector meets the negative axis and contains the last
quartic conjugate pair and one quintic root. The latter is
real by conjugation. The actual quartic is strictly positive
on the real axis by its Gamma integral, so the remaining
quartic pair has distinct nonreal arguments.

## Stabilized equations near the negative axis

For \(13/20\le u\le97/100\), use \(C=\cosh\log R\),
\(S=\sinh\log R\), and \(Y=\cos\theta\). Reversing a polynomial
replaces \(R\) by \(R^{-1}\) and reverses all skews, so it is enough
to check \(S\ge0\) with both skew signs allowed. The radius bounds
give \(S<8,16\), respectively. From the imaginary part of the
quartic and the real part after eliminating \(\sinh\alpha\) one gets
\[
\begin{aligned}
 S(CY+u\cosh\alpha)&=uC\sinh\alpha,\\
 4Y^2-2-4S^2+8uY\cosh\alpha/C+6v&=0.                        \tag{6}
\end{aligned}
\]
For the quintic put
\(p=5w\cosh\beta-1\), \(q=10y\cosh\chi-5w\cosh\beta-1\),
\(d=5w\sinh\beta\), and \(\Delta=10y\sinh\chi-d\). Dividing its
reciprocal expression by \(2\cosh(\xi/2)\) for
\(\xi=\log R+\mathrm i\theta\) gives
\(4\cosh^2\xi+2p\cosh\xi+q-2d\sinh\xi-\Delta\tanh(\xi/2)\).
The imaginary part and the real part after eliminating \(\Delta\) are
\[
\begin{aligned}
 2(C+Y)\{S(4CY+p)-dC\}&=\Delta,\\
 4Y^2+2pY+q
 +S^2\{8Y^2-4-8CY-2p+2pY/(C+1)\}
 +2dS(C-Y)&=0.                                              \tag{7}
\end{aligned}
\]
Both use \(C^2-S^2=1\). On an upper-half-plane ray \(C+Y>0\).
The linked skew improves the interval for \(\Delta\). In fact,
\(0\le w-y\le(1-w)/2\) by (1), and
\[
 |\Delta|\le 5(1-w)\sinh|\chi|
 +5w\bigl\{|\kappa|\cosh(\max(2|\chi|,|\beta|))
                  +2\sinh|\chi|(\cosh|\chi|-1)\bigr\}.       \tag{8}
\]
This follows by writing
\(2\sinh\chi-\sinh\beta
=(2\sinh\chi-\sinh2\chi)+(\sinh2\chi-\sinh\beta)\).

For numerically stable interval evaluations, write
\(C-1=S^2/(C+1)\) and
\[
 v=1-\tfrac43(1-u)+J_4(1-u)^2,\qquad
 y=1-\tfrac32(1-w)+J_5(1-w)^2.                             \tag{9}
\]
Here \(J_4\ge2/9\) and \(J_5\ge3/8\) by the lower bounds in (1)
and the integral form of the Taylor remainder. For an upper bound
use \(v\le F_2(u)\), \(y\le F_3(w)\), where
\(F_j(t)=t\sqrt{(1+jt^2)/(j+t^2)}\).
The quotients
\((F_j(t)-1+F_j'(1)(1-t))/(1-t)^2\) are decreasing for
\(1/2\le t<1\). Indeed, the third derivative of \(F_j\) has the
sign of \(h^4-8h^3-10h^2-2h+1\) when \(j=2\) and
\(3h^4-36h^3-30h^2-4h+3\) when \(j=3\), where \(h=t^2\).
These polynomials decrease on \([1/4,1]\) and are negative at
\(1/4\). Thus the integral Taylor remainder is largest at the
lower interval endpoint. In the two-dimensional ray checks, expand
(6)–(7) using \(C-1,1+Y,1-u,1-w\) and (9). Both necessary
conditions are enclosed on boxes in either \(u,S\) or \(w,S\).

## The near-flat interval

It remains to cover \(u\ge1-.18^2=.9676\). Put
\(\lambda=\sqrt{1-u}\), \(r=(1-w)/\lambda^2\), and set
\(z=-1+\lambda\zeta\). The bounds (1) give
\(21/20\le r\le4/3\). For the lower bound use the derivative at
\(u=.9676\) of \(\sqrt{15}\,u/\sqrt{16-u^2}\), which increases in
\(u\) and exceeds \(21/20\). Use the same \(J_4,J_5\) as in (9),
and set
\[
 A=\alpha/\lambda^4,\quad B=\beta/\lambda^4,\quad
 K=\kappa/\lambda^8,\quad D=(B-K\lambda^4)/2=\chi/\lambda^4.
\]
Let \(S_1(t)=\sinh(t)/t\), \(E_2(t)=(\cosh t-1)/t^2\),
\(S_3(t)=(\sinh t-t)/t^3\), with their continuous values at zero.
The polynomials after dividing by \(\lambda^4,\lambda^5\) are
\[
\begin{aligned}
 \widehat P_4={}&\zeta^4-4z\zeta^2+6J_4z^2
 +4uA S_1(\alpha)\lambda z(1-z)\zeta
 +4uA^2E_2(\alpha)\lambda^4(z+z^3),\\
 \widehat P_5={}&\zeta^5-5rz\zeta^3+10r^2J_5z^2\zeta
 +5wB S_1(\beta)\lambda z(1-z)\zeta^2
 +\lambda L z^2(1-z)\\
 &+\lambda^4\{5wB^2E_2(\beta)z(1-z+z^2)
                 +10yD^2E_2(\chi)z^2\}\zeta,\\
 L={}&10(-r/2+r^2\lambda^2 J_5)D S_1(\chi)
 +5w\{-K\lambda^2+\lambda^6(2D^3S_3(\chi)-B^3S_3(\beta))\}.
\end{aligned}                                                       \tag{10}
\]
For the quintic odd terms, the cancellation is
\(z-z^4+z^2-z^3=z(1-z)(1+z)^2\).
The leftover coefficient \(10y\sinh\chi-5w\sinh\beta\) divided by
\(\lambda^6\) is exactly \(L\). This follows from
\(y-w=-r\lambda^2/2+r^2J_5\lambda^4\)
and \(2\chi-\beta=-K\lambda^8\).

The replay uses (3) and (9) on each box. Globally,
\(0\le A\le.3\), \(0\le B\le2/3\) after reversing as needed,
\(|K|\le.21\), \(|D|\le.334\), \(J_4\le.3\), \(J_5\le.5\).
These follow already for \(\lambda\le1/4\) by using
\(u\ge15/16,w\ge11/12\) in (3), the Taylor bounds, and
\(r\le4/3\). The arguments of the hyperbolic functions have
absolute value below \(.01\). Their power series give
\(1\le S_1\le1.001\), \(1/2\le E_2\le.501\),
and \(1/6\le S_3\le.167\). Each upper estimate follows by
bounding the series after its constant term with a geometric
series of ratio at most \(.01^2\).

The scaled root radius is below \(20\). If \(T=|\zeta|\ge20\),
then \(|z|\le.30T\), \(|1-z|\le.35T\), and
\(|1-z+z^2|\le.1075T^2\), using \(\lambda\le1/4\).
The displayed scalar bounds imply \(|L|<2.5\).
Substitution in (10) bounds the sum of the absolute values of
the nonleading terms by \(.07T^4\) at order four and
\(.11T^5\) at order five. The replay verifies these last
rational comparisons as well as the scalar guards.

On the three rays \(z=\rho(-1+\mathrm i\lambda c)\), \(\rho>0\),
use \(c=2,3/2,4/5\). In the scaled variable write
\(\zeta=X+\mathrm i c(1-\lambda X)\). It is enough to check
\(-20\le X\le20\), and the replay certifies nonzero values
of (10) throughout the corresponding parameter boxes.
At \(\lambda=0\), the polynomials in (10) become
\[
 \zeta^4+4\zeta^2+4/3,\qquad
 \zeta(\zeta^4+5r\zeta^2+(15/4)r^2).
\]
The positive imaginary root magnitudes squared are
\(2\pm\sqrt{8/3}\) and \(r(5\pm\sqrt{10})/2\), respectively.
The three horizontal cuts at \(c=2,3/2,4/5\) separate them
in the desired order throughout \(21/20\le r\le4/3\).
For any fixed positive \(\lambda\), move to the symmetric
anchor as in the middle range, then decrease \(\lambda\)
to zero along \(w=u^{16/15}\). The corresponding roots
converge to the simple limiting roots, and the moving rays
never meet one. This gives the required sector counts in the
near-flat range.

The three ranges cover \(0<u<1\). The positive rescalings used
in normalization preserve root arguments. The original
ordinary factorial polynomials lie in \(\mathbb Q[t]\),
so their complex coprimality gives a rational Bézout identity.
A zero monomial coefficient leaves a single monomial whose
moments are nonzero.
