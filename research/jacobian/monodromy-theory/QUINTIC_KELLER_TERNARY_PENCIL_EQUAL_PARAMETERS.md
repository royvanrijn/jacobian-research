# Weighted restrictions in the equal-parameter ternary pencil

**Proposition QuinticKellerTernaryPencilEqualParameters.** Use the coordinates
of the [first-equations theorem](QUINTIC_KELLER_TERNARY_PENCIL_FIRST_EQUATIONS.md).
In a Keller map with \(a=0\), use the
[unequal-parameter exclusion](QUINTIC_KELLER_TERNARY_PENCIL_UNEQUAL_PARAMETERS.md)
to put \(v_1=u_5=v_2=v_5=0\) and \(u_0=u_4=k\).
With weights \((1,2,0)\), put \(s=q/y\) and \(h=v_4-v_0\). The weight-four
part of \(F_2\) has the form \(y^2 B(s)\), where \(B(s)=k s^2+f s+B_0\).
Then either \(B=0\), or \(B=k(s+h)^2\) with \(k\ne0\).
In the latter case \(F_1\) has weight two and its leading part is
\[
            y(s+h)\,L\!\left(\frac{z(s+h)^2+C_0(s)}{B(s)}\right),
\]
where \(L\) is a nonconstant univariate polynomial of degree at most two
and \(C_0(s)\) is defined below. If \(L\) is linear then \(s+h\mid C_0(s)\).
If \(L\) is quadratic then \((s+h)^2\mid C_0(s)\).
The [lower determinant equations](QUINTIC_KELLER_TERNARY_PENCIL_TRANSVERSE_PARAMETER.md) exclude both alternatives.
The argument extends from \(\mathbb C\) to characteristic zero by the
same coefficient-field embedding as in the preceding exclusion.

## Leading weighted forms and rational invariants

The [exact checker](verify_quintic_ternary_pencil_equal_parameters.py)
uses the rank-twenty cubic equation to verify that \(F_2,F_3\) have
weight at most four and that their weight-four parts are
\[
        y^2 B(s),\qquad y^2 C(s,z),\qquad
        C(s,z)=z(s+h)^2+C_0(s),\quad C_0(s)=2v_0s^2+g s+C_* .
 \tag{1}
\]
Here \(f,g\) are the coefficients of \(x^2y\) in the second and third
cubic components, and \(B_0,C_*\) are the coefficients of \(y^2\) in
the second and third quadratic components. The checker also verifies
that \(F_1\) has weight at most four with weight-four part a constant
multiple of \(y^2\).

Assume \(B\ne0\). In the chart \((x,y,z)=(tX,t^2,z)\), let \(s=X^2+z\),
\(r=C(s,z)/B(s)\) and \(w=t^4B(s)\). The two leading functions are
\(w,rw\). If the leading weighted form of \(F_1\) is \(t^j A(X,z)\),
its determinant with those leading functions vanishes. Indeed the
coefficient has weight \(j+5>0\), while the Keller determinant is
constant. Thus the same wedge-product argument as in the
[ordinary conic pencil](QUINTIC_KELLER_CONIC_PENCIL_FIRST_INTEGRALS.md)
gives \(A^4/B^j\) algebraic over \(\mathbb C(r)\).
This function is even in \(X\) and therefore belongs to
\(\mathbb C(s,z)=\mathbb C(r)(s)\). It follows that \(A^4/B^j\in\mathbb C(r)\).

No such leading form can have odd weight \(j\). In that case \(A\) is odd
in \(X\), so \(A^4/B^j\) vanishes at \(X=0\), or \(z=s\), for generic \(s\).
On this curve the ratio \(r\) equals
\[
                     \frac{s(s+h)^2+C_0(s)}{B(s)},
\]
a nonconstant function because the numerator has degree three and the
denominator degree at most two. A nonzero rational function of \(r\)
cannot vanish on this curve. The argument works for the leading weighted
form of any constant linear combination of \(F_1,F_2\) in place of
\(F_1\), provided it is nonzero.

## The nonconstant case

If \(B\) is nonconstant, \(F_1\) cannot have weight four. Its leading
form would be a nonzero scalar times \(y^2\), and the determinant
condition would give \(B'(s)(s+h)^2=0\) by differentiating (1).
Weight zero would make \(F_1\) a univariate polynomial in \(z\)
of degree five, contrary to its being critical-point-free.
Thus \(F_1\) has weight two.

For its weight-two form \(t^2 A(s,z)\), the polynomial \(A(s,z)\)
is affine in \(s\). The same wedge-product argument applied to
the ratio \(A^2/B\) instead gives \(A^2/B\in\mathbb C(r)\).
Every root of \(B\) in the rational-function variable \(s\) must then
have even multiplicity, by taking valuations in \(\mathbb C(r)(s)\).
Hence \(B=k(s+s_0)^2\) with \(k\ne0\). It follows that
\(A/(s+s_0)\in\mathbb C(r)\), using again relative algebraic
closedness. Write \(A=(s+s_0)L(r)\). The rational function \(L\)
is polynomial since \(r\) is affine and nonconstant in \(z\) over
\(\mathbb C(s)\) and \(A\) is polynomial in \(z\).
Its degree is at most two because the original degree-five leading
form is \(z^5\) and a weight-two monomial in \(F_1\) has ordinary
degree at most four.

If \(L\) were a nonzero constant, the leading weighted part of \(F_1\)
would be a nonzero multiple of \(x^2+y(z+s_0)\). Its gradient has
a nondegenerate zero at \(x=y=0,\ z=-s_0\). The gradient of
\(t^{-2}F_1(tx,t^2y,z)\) then has a formal zero over
\(\mathbb C[[t^{-1}]]\), found recursively using the invertible
leading Hessian. This contradicts the Keller equation after
substitution, since the Jacobian of the scaling is invertible
over \(\mathbb C((t^{-1}))\).
Thus \(j_0=\deg L\) is one or two. The coefficient of \(z^{j_0}\)
in \(A=(s+s_0)L(r)\) is a nonzero scalar times
\((s+h)^{2j_0}/(s+s_0)^{2j_0-1}\). Since \(A\) is polynomial in \(s,z\),
we must have \(s_0=h\).

Put \(M=s+h\). If \(L\) is linear, the constant coefficient in \(z\)
of \(A=M L(r)\), with \(r=(z+C_0/M^2)/k\), forces \(M\mid C_0\).
If \(L\) is quadratic, the coefficient of \(z\) first forces
\(M\mid C_0\). If that divisibility were exact to order one,
the constant coefficient would have a simple pole from \(C_0^2/M^3\)
that cannot cancel against the lower terms. Thus \(M^2\mid C_0\).

## Excluding a nonzero constant

Suppose \(B\) is a nonzero constant. If \(F_1\) has weight four,
subtract the multiple \(\gamma F_2\) that cancels its weight-four
part. Call the resulting first component \(P\), taking \(\gamma=0\)
if no subtraction is needed. The preceding argument excludes odd
leading weight for \(P\). Weight zero is impossible. If \(\gamma=0\)
it would leave \(F_1\) univariate of degree five, and if \(\gamma\ne0\)
the term \(-\gamma z^3 x^2\) from the ordinary degree-five leading
part would have weight two. Thus \(P\) has weight two.

Write its leading weighted form as \(t^2 A(s,z)\), again affine
in \(s\). The ratio \(A^2/B\) belongs to \(\mathbb C(r)\). Since
\(B\) is constant, \(A=L(r)\) for a univariate polynomial \(L\).
If \(L\) had positive degree \(i\), the highest \(z\) coefficient
of \(A\) would have degree \(2i\) in \(s\) from (1), contrary to
affineness. Hence \(A\) is a nonzero constant. The monomial
\(-\gamma z^3x^2\) once again shows that \(\gamma=0\). We obtain
\[
                          F_1=c y+x\,S(z)+T(z),
\]
where \(c\ne0,\ \deg S\le3\) and \(T\) has leading term \(z^5\).

On each fibre \(F_1=u\), substitute
\(y=(u-xS(z)-T(z))/c\). The resulting planar pair \(F_2,F_3\)
in \(x,z\) has constant nonzero Jacobian by the original determinant
and the constant \(y\) coefficient of \(F_1\). Its degrees are
exactly ten and thirteen. For \(F_2\) the only weight-four
term is \(B y^2\), giving degree ten after substitution. A monomial
of smaller weight contains at most one \(y\), and the ordinary
degree bound then gives substituted degree at most nine.
For \(F_3\) the term \(z^3 y^2\) gives degree thirteen.
No other term attains that degree, since every term has weight
at most four and the other \(y^2\)-terms have smaller \(z\)-degree.

The [retained planar degree bound](../../BREAKTHROUGH_TARGET_ASSESSMENT.md#four-variable-hessian-conjecture-and-planar-degree-bounds)
makes this pair a planar polynomial automorphism. The component
degrees of a planar polynomial automorphism divide one another,
by the [division lemma for Jung's theorem](https://arxiv.org/pdf/math/0408077).
Ten and thirteen violate it. This excludes every nonzero constant \(B\).
