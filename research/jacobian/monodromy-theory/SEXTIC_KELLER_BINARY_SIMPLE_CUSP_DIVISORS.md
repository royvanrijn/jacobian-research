# Ordinary cusps of binary leading tuples

**Theorem SexticKellerBinaryOrdinaryCuspExclusion.**
Let \(F\) be a characteristic-zero three-variable Keller map
of ordinary component degrees \((6,5,4)\) or \((6,5,5)\).
Suppose its highest forms \(H_i\) belong to \(k[x,y]\)
and have generic derivative rank two as a map of the binary
plane. Over the algebraic closure consider each line
dividing the gcd of the \(2\times2\) binary derivative minors.
At the generic point of such a line require that a component
\(H_a\) be nonzero there and that the following condition
hold. On the generic fixed level \(H_a=S\), expand the other
two \(H_i\) in a local parameter transverse to the line.
After subtracting the constant values, their linear terms
vanish and their quadratic and cubic coefficient vectors
are linearly independent. Then no such Keller map exists.

This condition describes an ordinary cusp of the ordered
pair on a level of \(H_a\). It is independent of the choice
of local parameter with nonzero linear part. For a line
normalized to \(x=0\), one can test it by setting
\(t=x/y\), \(h_i(t)=H_i(t,1)\), and using the two series
\(h_i(t)h_a(t)^{-d_i/d_a}\) for \(i\ne a\), where \(d_i=\deg H_i\).
The factors \(S^{d_i/d_a}\) on the fixed level do not
affect the independence test. In particular any remaining
binary tuple in these two profiles needs a common derivative
line where this condition fails.

## Reference descent and the local tangent parameter

Work over an algebraic closure. Let the reference be
\(G_i(\widetilde x,\widetilde y)\) with highest forms \(H_i\),
where \(\widetilde x,\widetilde y\) have leading forms \(x,y\).
Let \(r\ge1\) be the first degree loss in the error and
\(E_i\) its parts of degrees \(d_i-r\).
Rescale all components by their degrees with \(\epsilon\).
The determinant has order \(N=\sum d_i-3\), which is
twelve or thirteen. For any first error between positive
degrees, \(r\le5\), so \(2r<N\).

Put \(n=H_x\times H_y\). The linear determinant equation
at loss \(r\) gives \(\partial_z(n\cdot E)=0\).
This is the calculation in the
[saturated leading-pencil descent](SEXTIC_KELLER_MIXED_DEGREE_LEADING_REDUCTIONS.md#descending-through-the-lower-terms),
now allowing nonprimitive minors. Absorb \(E(x,y,0)\)
into the binary reference. The remaining \(E\) is tangent
to the highest binary map and vanishes at \(z=0\).
Over the fraction field write \(E=H_x U+H_y V\).
The scalars \(U,V\) are polynomial away from potential
poles on lines of the common derivative divisor, since
a nonzero minor solves for them. We show they are
regular also at the generic point of every such line.

Normalize one of the lines to \(x=0\), and choose \(H_a\)
as in the hypothesis. Work over an algebraic extension
\(K\) of \(k(S)\) containing the chosen value of \(y\)
at \(x=0,H_a=S\). The derivative in \(y\) of \(H_a\)
is nonzero there. At fixed rescaled value
\(F_{a,\epsilon}=S\), use \(t=\widetilde x_\epsilon/
\widetilde y_\epsilon\) and \(z\) as formal coordinates,
solving for the radial variable. The binary reference
on its own first level \(G_{a,\epsilon}=S\) then has
other components \(b(t,\epsilon),c(t,\epsilon)\).
Choose the component \(b\) with a nonzero quadratic
coefficient in the leading level. Center its critical
point and take a formal square root so that \(b=b_0(\epsilon)+s^2\),
where \(s\) is a regular parameter at the leading line.
After replacing \(c\) by \(c-c(0,\epsilon)-\beta(\epsilon)s^2\),
where \(\beta\) is its quadratic coefficient before this
subtraction, its leading term has order three in \(s\)
with nonzero coefficient \(\gamma\). These choices use
the ordinary cusp hypothesis. The formal critical point
exists by the nonzero leading second derivative of \(b\);
the square root divides out the resulting square of its
centered coordinate.

Make the same output subtractions in the actual chart.
Write the resulting components as \(s^2+U_*\) and
\(c(s,\epsilon)+V_*\), with errors starting at order \(r\).
At first order along the level, fixing \(F_{a,\epsilon}=S\)
subtracts the radial component of the tangent error \(E\).
Thus the residual tangent parameter \(\tau\) in the
\(s\)-direction satisfies
\[
              (U_*)_r=2s\tau,\qquad (V_*)_r=c_0'(s)\tau.
 \tag{1}
\]
Both coefficients in (1) are regular in \(s\), and their
values at \(z=0\) vanish, since the original \(E\) vanishes
there. The parameter \(\tau\) at this stage is only
known to be rational along the generic chart.

## The pole comparison eliminates the local denominator

On the generic Laurent chart, replace \(s\) by
\(\sqrt{s^2+U_*}=s+O(\epsilon)\).
As an \(\epsilon\)-formal substitution this is invertible
over the \(s\)-Laurent field. In the resulting coordinates
the first two actual outputs are fixed by \(S\) and the
new \(s\), and the determinant implies that the derivative
of the third component in \(z\) vanishes modulo
\(\epsilon^N\). The Jacobian of the coordinate change
is an \(\epsilon\)-unit at generic \(s\), so the order
bound is preserved. The normalized coefficients are
therefore Laurent series in \(s\) over \(K\).
Consequently there is a scalar \(R(s,\epsilon)\) with
\[
       c+V_*=\sum_{m\ge0}\frac{U_*^m}{m!}\delta^m R
                 \pmod{\epsilon^N},
          \qquad \delta=(2s)^{-1}\partial_s .
 \tag{2}
\]
The notation for the root and Taylor series is the
formal one used in the
[tangent-conic comparison](SEXTIC_KELLER_TANGENT_CONIC_CUBIC_OUTER.md#comparing-the-first-two-possible-poles).
By (1) the components \(R_i\) for \(i\le r\) are exactly
the reference coefficients \(c_i\).

Write \(\alpha_i=[s]c_i\) and \(u=(U_*)_r|_{s=0}\).
If there is a first nonzero \(\alpha_p\) for \(p<r\), the
only possible pole at order \(r+p\) besides the scalar
\(R_{r+p}\) has residue coefficient \(\alpha_p u/2\).
Thus \(u\) is constant in \(z\), hence zero by its
value at \(z=0\).
If all \(\alpha_i\) for \(i<r\) vanish, at order \(2r\)
the only possible pole besides the scalar \(R_{2r}\) has
coefficient
\[
                            \alpha_r u/2+3\gamma u^2/8 .
\]
The linear term comes from \(\delta c_r\), the quadratic
from the leading cubic term in \(\delta^2 c_0/2\).
All other terms are regular because lower reference
linear coefficients vanish and the quadratic coefficient
was removed. Since \(\gamma\ne0\), regularity makes \(u\)
algebraic over \(K\). It is a rational function of the
free \(z\), so is constant and vanishes at \(z=0\).
All the orders used are below \(N\).

In either case \(s\mid (U_*)_r\), giving regularity of
\(\tau\) by (1). The leading first level \(H_a\) together
with \(s\) gives regular coordinates at the generic
point of the line after the indicated scalar extension,
so the radial coefficient from \(E_a\) and the transverse
coefficient \(\tau\) together give regular \(U,V\).
The chart has ramification index one in the line
parameter, so this regularity excludes poles in the
original height-one localization as well.

Every potential denominator has now been excluded.
Thus \(U,V\) are polynomials of degree \(1-r\) by
homogeneity, and can be added to the two source
parameters of the reference. This absorbs the first
error. After repeating through the positive degrees and
adding the constant errors, the map factors through
two polynomials, a contradiction.

## Quadratic critical points with higher odd terms

**Lemma SexticKellerBinaryQuadraticCriticalPoleBound.**
Use the binary source and reference setup from the ordinary-cusp
proof, with rescaled determinant of order \(N\). On a fixed
generic level of a nonzero leading component, suppose the
chosen second component has a nonzero quadratic term at a
common derivative line. Write \(s\) for the centered parameter
in which this component is a scalar plus \(s^2\). Let \(m\)
be the first odd exponent in the leading series of the third
component in \(s\). At first error loss \(r\), the tangent
parameters are regular at that line if \(r(m+1)/2<N\).
They are also regular for \(m=5,r=4,N=12\) when the
two components used for the chart have degrees at most five.

**Theorem SexticKellerBinaryShortOddContactExclusion.**
Let a characteristic-zero three-variable Keller map have ordinary
profile \((6,5,4)\) or \((6,5,5)\), with binary highest tuple
of generic derivative rank two. At each common derivative line,
suppose the two lower-degree components can be used as the first
and second components of the preceding chart, with the first
nonzero generically on that line and the second having a nonzero
quadratic term on the fixed level. If the first odd exponent of
the sextic component in the centered leading parameter is at most
five at every such line, there is no Keller completion.

To prove the lemma use the reference descent and coordinates of
the ordinary-cusp proof. Write the actual second component as a
scalar plus \(s^2+U\), where \(U\) has order at least \(r\) in
\(\epsilon\). In the rectified generic chart
\(s'=(s^2+U)^{1/2}\) the third component is
\(R(s',\epsilon)\) modulo \(\epsilon^N\), where
\(R_i\in K((s'))\) and \(K\) is the algebraic coefficient
field of the fixed level. Indeed at fixed first component
and \(s'\) the \(z\)-derivative vanishes at orders below \(N\)
by the determinant. The leading series \(R_0\) has its first
odd term at exponent \(m\). Put \(u=U_r|_{s=0}\).
As in the ordinary-cusp proof, \(u\) is in \(K(z)\) and
vanishes at \(z=0\) since the binary reference has absorbed
\(E(x,y,0)\). If nonzero, it is transcendental over \(K\).

Give \((\epsilon,s)\) weights \((2,r)\). If \(u\ne0\),
the minimum-weight part of \(s^2+U\) is \(s^2+\epsilon^r u\).
Among terms \(\epsilon^i s^\ell\) of \(R\) with \(i<N\)
and \(\ell\notin 2\mathbb Z_{\ge0}\), let \(W\) be the
least weight. It exists since the term at \(i=0,\ell=m\)
occurs. Each Laurent coefficient has a finite lower bound on
its exponents, so \(W\) is finite and \(W\le rm\).
Remove terms of weight at most \(W\) with nonnegative even
\(\ell\) by subtracting their evaluations in \(s^2+U\).
These evaluations are regular. At weight \(W\) the term
\(\epsilon^i s^\ell\) has the binomial expansion with
terms proportional to
\[
        \binom{\ell/2}{k}\epsilon^{i+rk}s^{\ell-2k}u^k .
\]
If \(\ell<0\), use \(k=0\); otherwise \(\ell\) is odd and
use \(k=(\ell+1)/2\). Both give a pole, at an order less
than \(N\): in the second case it is
\((W+r)/2\le r(m+1)/2\). Terms at the same \(\epsilon\)-order
and pole order contributed by different pairs \((i,\ell)\)
have distinct powers of \(u\), since \(i+rk\) fixes \(i\)
once \(k\) is known and \(2i+r\ell=W\) fixes \(\ell\).
The resulting pole coefficient is a nonzero polynomial in
the transcendental \(u\), contradicting regularity.
Thus \(u=0\). This makes \(U_r/(2s)\), the first
transverse tangent coefficient, regular. The radial
coefficient is then also regular because the derivative
of the fixed leading component in the radial direction is
a unit on the chart.

For the endpoint \(m=5,r=4,N=12\), compute also at
order \(N\). In the rectified chart the second actual
coordinate is a scalar on the chosen first level plus
\((s')^2\). By the chain rule the derivative of \(R_N\)
in \(z\) is
\[
       \partial_z R_N=\frac{\kappa}{2s'\sigma(s')},
\]
up to an overall sign, where \(\sigma\) is the leading
Jacobian of the chart from the binary source variables to
the first component and \(s'\). It is a scalar unit at
\(s'=0\), independent of \(z\), by transversality and the
nonzero radial derivative. Thus \(R_N\) is a scalar
Laurent series plus a term linear in \(z\) with at most
a simple pole. This \(z\)-dependent term has weight at
least \(2N-r=20\). Include the scalar part of \(R_N\)
in the minimum-weight comparison. If \(W<20\), the same
argument gives a nonzero pole coefficient through order
\(N\), with no interference from the \(z\)-dependent term.
If \(W=20\), the coefficient of
\(\epsilon^{12}s^{-1}\) has the form
\[
 a_{0,5}\binom{5/2}{3}u^3+a_{4,3}\binom{3/2}{2}u^2
          +a_{8,1}\binom{1/2}{1}u+a_{12,-1}+c z .
\]
Here the \(a_{i,\ell}\) are scalar coefficients of the
normalized series, with \(a_{0,5}\ne0\). Since both
components of the fixed-level chart have degree at most
five, their error parts at loss four have degree at most
one. After the subtraction at \(z=0\) they are scalar
multiples of \(z\). The first-order change from fixing
the scalar level has a binary coefficient, so \(u\) is a
scalar multiple of \(z\) as well. If \(u\ne0\), the
coefficient displayed above has a nonzero cubic term
in \(z\), again a contradiction. This gives regularity
also at the endpoint.

For the theorem, \(m\) is either three or five because
the common derivative line kills the linear term on the
fixed level. At losses \(r\le4\), the bound applies
directly or by the endpoint calculation. At losses
\(r\ge5\), the first error in either lower-degree component
has degree at most zero, so it vanishes after subtracting
the binary value. Hence \(U_r=0\) and there is no transverse
tangent pole. Every common derivative line therefore
admits regular tangent parameters at every loss.
The polynomial tangent absorption and descent in the
ordinary-cusp proof now give a factorization through two
source polynomials and a zero determinant.

For example, take the binary highest forms
\[
 A=y^6+xy^5,\qquad B=x^5+x^3y^2+xy^4,\qquad C=x^4.
\]
The two minors against \(C\) leave only \(x=0,y=0\) as possible
common derivative lines, since \(x^2+2y^2\) and \(6y+5x\)
are coprime. The minor of \(A,B\) is nonzero generically on
\(x=0\). On \(y=0\), fixing a generic value of \(C\) fixes
\(x\ne0\). The centered \(B\)-parameter is odd in \(y\) with
nonzero linear coefficient, and \(A\) has first odd exponent
five in that parameter. Thus the theorem excludes this tuple
with a common lower leading factor.
