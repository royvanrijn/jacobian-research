# Exclusion of the conductor-twenty-two sequence \((12,8,6,11)\)

**Theorem ConductorTwentyTwoTwelveEightSixElevenExclusion.** Let
\((P,Q)\) be a complex polynomial plane Keller map in a source frame
with \(P\) monic of positive longitudinal degree. Its generic exact
directional order cannot have reduced characteristic sequence
\((12,8,6,11)\).

## Parameter identities

Over \(K=\mathbb C(y)\), use monic generators \(W,V\) of degrees eight
and twelve. The first approximate root gives
\[
                       V^2=W^3+kW+R,\qquad \deg R=6 .
\tag{1}
\]
Indeed that root is monic cubic in \(W\), and the lower powers of \(W\)
have distinct residues of their degrees modulo twelve. Complete the
square in \(V\), translate \(W\) to depress the cube and absorb the
constant into the degree-six root, as in the
[three-stage root identity](CONDUCTOR_FOURTEEN_THREE_STAGE_EMBEDDING.md#the-two-approximate-root-identities)
using [Assi–García-Sánchez, Propositions 2 and 4](sources/assi-garcia-sanchez-1407.0490v1.pdf).

After the source translation over \(K\) specified below, take
\(\rho^{20}=\lambda\) with \(\lambda\in K^\times\) the target
determinant before scaling the source by \(\rho\) and the
generators by \(\rho^{-8},\rho^{-12}\). This preserves (1) and sets
the target determinant to one.

Write \(W=z^8+wz^7+\cdots\) and \(R=\sigma z^6+\theta z^5+\cdots\)
with \(\sigma\ne0\). In the
[parameter identity](CONDUCTOR_TWELVE_CABLE_EXCLUSIONS.md#1-common-source-and-target-identities)
\[
                 \operatorname{Jac}_{z,y}(W,V)=J\rho+\partial_z H(W,V)
\]
use \(V=W^{3/2}+kW^{-1/2}/2+RW^{-3/2}/2+O(z^{-20})\).
The polynomial Jacobian is the sum of the polynomial parts of
\(k'(\sqrt W)_z\) and \(\tfrac12 W^{-3/2}(W_z R_y-W_y R_z)\).
Its degree is at most three. Since \(\langle6,8,11\rangle\)
has no positive degree at most four, the evaluated \(H\) is constant.
The degree-three and degree-one equations give \(k'=\sigma'=0\).
After \(\sigma'=0\), the degree-twelve coefficient of the numerator
of the second term is \(8\theta'-6\sigma w'\), giving
\[
                  (4\theta-3\sigma w)'=J\rho .
\tag{2}
\]
The nonzero constant \(\sigma\) has weight eighteen, so
\(\rho^{18}\) and \(\rho^{20}\) lie in \(K\). Thus \(\rho^2\in K\).

## Exact square-root coefficients

Write \(W=L^2+M\) with \(L\) monic quartic and \(\deg M\le3\).
The square root \(S=(W^3+kW)^{1/2}\) differs from \(V\) by a
term of degree exactly \(-6\). Expanding in \(M/L^2\) gives
\[
 V=L^3+\tfrac32 LM+\tfrac38\operatorname{quot}(M^2,L),\qquad
 S-V=\frac{6\operatorname{rem}(M^2,L)+8k}{16L}
                -\frac{M^3}{16L^3}+O(z^{-8}).
\]
This is the same expansion as in the
[degree-five approximate-root calculation](CONDUCTOR_TWENTY_TWELVE_EIGHT_FIVE_EXCLUSION.md#the-approximate-root-equation).
Thus \((6\operatorname{rem}(M^2,L)+8k)L^2-M^3\) has degree six.
If \(\deg M\le2\), the first term vanishes. Then \(M\) has
degree two, and \(M^2\) is a constant multiple of \(L\) plus a
constant. Both \(W,V\) would be polynomial in \(M\), contrary
to the full function field.

Thus \(M=r m\) with \(r\ne0\) and \(m\) monic cubic. The form
\(6\operatorname{rem}(M^2,L)+8k\) must be linear with leading
coefficient \(r^3\) to cancel the degree-nine term of \(M^3\).
Translate its root to zero before scaling. Its scaled value is
\(r^3z\). Put \(m=z^3+az^2+bz+c\) and \(\Delta=a^2-4b\).
The coefficients in degrees eight and seven of \(m^3-zL^2\)
determine the top coefficients of \(L\):
\[
 L=z^4+\tfrac32 a z^3+(\tfrac32 b+\tfrac38a^2)z^2+
        (\tfrac32 c+\tfrac34 ab-\tfrac1{16}a^3+j)z+\ell.
\]
The degree-three and degree-two coefficients of
\(\operatorname{rem}(m^2,L)\) vanish. The
[exact checker](verify_conductor_twenty_two_twelve_eight_six_eleven_exclusion.py)
gives \(c=-a\Delta/8+2j\) and
\(\ell=-(5a^4-16a^2b-128aj-16b^2)/64\).
The remainder has linear coefficient \(\Delta(a\Delta-64j)/128\)
and constant coefficient \((\Delta^3+2(a\Delta-32j)^2)/512\).
The degree-six coefficient of \(m^3-zL^2\) is \(-2j\).
Comparing the linear form and the leading coefficient of
\(V^2-S^2=-(S-V)(S+V)\) gives
\[
 r=\frac{3\Delta(a\Delta-64j)}{64},\quad
 k=-\frac{3r^2}{2048}\big(\Delta^3+2(a\Delta-32j)^2\big),\quad
 \sigma=-r^3j/4 .
\tag{3}
\]

## Rational parameter curve

Let \(X=r^2,\ Y=X^2/\Delta\), both nonzero by (3). These
belong to \(K\) because \(r\) has weight five, \(\Delta\) has
weight two and \(\rho^2\in K\). Put
\(A=64/3,\ B=32jr^3,\ C=-2048k/3\), so \(A,B\) are nonzero
constants and \(C\) is constant. Equations (3) give
\[
 r(a\Delta-32j)=(AY+B)/X,\qquad
 X^9=Y^3\big(CX^2-2(AY+B)^2\big).
\tag{4}
\]
If \(Y\) is constant then \(X\) is constant by the monic
equation in (4). This fixes \(r,\Delta,j,a\) by (3), and
hence \(W,V\) are constant in the parameter. This contradicts
(2). It remains to exclude a nonconstant rational map into
the curve (4).

At infinity over the \(Y\)-line, normalize the valuation by
\(v(Y)=-1\). In the polynomial
\(F=X^9-CY^3X^2+2Y^3(AY+B)^2\), the constant term has
valuation \(-5\). Every root has valuation \(-5/9\), because
at that value the middle term has strictly greater valuation,
and at any other value a unique term would have the least
valuation. Hence the degree-nine polynomial is irreducible
over \(\mathbb C(Y)\), and its smooth projective normalization
has ramification index nine at infinity.

If \(C=0\), the same valuation calculation at \(Y=-B/A\)
gives index nine and at \(Y=0\) gives index at least three.
These contributions give genus at least one by Riemann–Hurwitz.
If \(C\ne0\), away from \(Y=0\) the nonzero critical roots over
the \(Y\)-line satisfy
\[
 9X^7=2CY^3,\qquad X^2=\frac{18(AY+B)^2}{7C}.
\]
Choose \(c_0^2=18/(7C)\) in \(\mathbb C\). The two signs
\(X=\pm c_0(AY+B)\) give equations
\((AY+B)^7=\pm dY^3\) for a nonzero constant \(d\).
Their roots are disjoint and neither zero nor \(-B/A\).
A multiple root must satisfy \(4AY=3B\), as seen by
differentiating \((AY+B)^7/Y^3\). The numerator of that derivative
has a simple zero there, so the only possible multiple root is
double and lies in at most one of the two equations.
Thus at least twelve distinct simple roots remain.

At each simple root \(F_{XX}=14CY^3\ne0\).
Also \(F_Y\ne0\), by differentiating \(F(\pm c_0(AY+B),Y)\)
at the root, where \(F_X=0\). Along this substitution,
\(F=X^2(X^7-2CY^3/9)\), whose derivative is nonzero by
simplicity. These points are smooth and have ramification
index two over the \(Y\)-line. Including infinity therefore
gives \(2g-2\ge-18+8+12>0\).
In both cases the irreducible curve has positive genus, so a
nonconstant rational map from \(\mathbb P^1_y\) is impossible
by Riemann–Hurwitz. This contradicts (4).
