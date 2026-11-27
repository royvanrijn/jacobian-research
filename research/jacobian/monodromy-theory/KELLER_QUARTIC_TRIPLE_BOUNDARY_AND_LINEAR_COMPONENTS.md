# Triple points and repeated support in quartic Keller boundaries

Let \(F=(P,Q,R)\) be a complex three-variable polynomial Keller map.
Write \(R_i\) for the homogeneous degree-\(i\) part of \(R\). For
a quartic \(R\), the plane leading curve is \(C=\{R_4=0\}\subset
\mathbb P^2\), the section at infinity of every compactified fibre.

**Theorem KellerQuarticOrdinaryTripleBoundaryExclusion.** The curve
\(C\) cannot be reduced with all singularities ordinary nodes,
cusps, tacnodes, or ordinary triple points.

**Theorem KellerQuarticTransverseRepeatedSupportExclusion.** Suppose
that \(R_3\) is nonzero at every singular point of \(C_{\rm red}\),
and that every zero of \(R_3\) on a repeated component of \(C\)
is a simple zero of its restriction to that component. Then \(R\)
cannot be a quartic Keller component. These hypotheses permit
nonreduced leading curves.

Here the first theorem extends the [mild quartic boundary exclusion](KELLER_QUARTIC_COMPONENT_BOUNDARY_OBSTRUCTION.md)
at every ordinary triple point, including a vanishing projective
quadratic jet.

## A residue restriction for components linear in a source variable

**Proposition KellerLinearSourceCoefficientResidue.** Suppose
\(R=a(x,y)z+b(x,y)\) with \(a\ne0\). On the normalization of every
multiplicity-one component of \(a=0\), the residue of
\(dx\wedge dy/a\) has a rational primitive polynomial in the
restriction of \(b\). In particular, the residue is regular at
every point of that normalization lying over the affine plane.

Put \(\kappa=\det JF\). On the open \(a\ne0\) of the fibre \(R=t\)
use \(x,y\) as coordinates, with \(z=(t-b)/a\). Its volume form is
\(\alpha=dx\wedge dy/a\), with the polynomially supplied primitive
\(\gamma_t=\kappa^{-1}P|_{R=t}\,d(Q|_{R=t})\), by the same
identity as [Wang–Yan, equation (3.28)](https://arxiv.org/pdf/2610.36495v1).
The coefficients of \(\gamma_t\) on this open are polynomial in
\(t\), with rational functions of \(x,y\) as coefficients. All
these primitives have derivative \(\alpha\).

On the smooth locus of a multiplicity-one component, away from
the other components, \(da\) is nonzero. The differentials \(da\)
and \(db\) are independent there. Otherwise one can choose \(z\)
so that \(z\,da+db=0\) at a point of \(a=0\), giving a critical
point of \(R\). Thus \(b\) is nonconstant on this component, and
\(s=a,\ u=b\) are local plane coordinates along its general point.

Expand \(\gamma_t=A(s,u,t)\,ds+B(s,u,t)\,du\) along \(s=0\),
and write \(h_t\) for the coefficient of \(s^{-1}\) in \(A\).
These coefficients are rational functions on the normalized
component. For instance, at its generic point the completion is
\(L[[s]]\), with \(L\) the component function field and coefficient
field lifting the separable extension \(L/\mathbb C(b)\) while
fixing \(u=b\). Thus the coefficients of \(h_t\) are in \(L[t]\).
Equivalently one can take the expansions in the analytic
coordinates \(a,b\) and follow each branch over the \(b\)-line.
In \(d\gamma_t=(\partial_s B-\partial_u A)\,ds\wedge du\),
the coefficient of \(s^{-1}\) in \(\partial_s B\) is zero.
It follows that \(dh_t=-\operatorname{Res}(\alpha)\) on the
component. In particular all the positive-degree-in-\(t\)
coefficients of \(h_t\) are constant on the component. Hence
\(h_t=h_0+c(t)\) with \(c\in\mathbb C[t]\).

At a general point of the component with \(b=t_*\), restrict
\(\gamma_{t_*}\) to the local normal slice \(b=t_*\). In the fibre
this is the regular curve \(z=0,\ b=t_*\). The original primitive
is regular there, so its normal residue \(h_{t_*}\) at the chosen
point vanishes. Evaluating this at all general points gives
\(h_0=-c(b)\) in the function field of the component. Its derivative
is regular on the whole affine normalization because \(b\) is
regular there. This proves the proposition.

For reduced \(a\) of degree at most three this excludes all affine
singularities of \(a=0\). When distinct smooth branches meet,
the residue on one smooth branch has a pole at their meeting:
express \(a\) as a local equation of that branch times the other
factors and a unit. For the only singular analytic branch of a
reduced plane cubic, an ordinary cusp, the same claim follows
by writing its equation as \(v^2-u^3=0\). Up to a unit the residue
is \(du/(2v)\), with a double pole after \(u=\tau^2,\ v=\tau^3\).
The cubic singularity cases also appear in
[Wang–Yan, immediately after Lemma 3.13](https://arxiv.org/pdf/2610.36495v1).

## Local analysis at an ordinary triple point

Choose projective source coordinates so that the triple point is
\([0:0:1:0]\). Write \(G_t(u,v,w)\) for the local equation with
\(w=0\) the plane at infinity. Its restriction to \(w=0\) starts
with a binary cubic having three distinct factors. If \(G_t\)
has a linear term, the fibre is smooth there. If its quadratic
jet has rank at least two, the double-point splitting argument
in the [quartic boundary proof](KELLER_QUARTIC_COMPONENT_BOUNDARY_OBSTRUCTION.md#ordinary-double-points-of-the-plane-leading-curve)
gives a holomorphic extension of the residue to a resolution.

A nonzero quadratic jet of rank one must be a multiple of \(w^2\),
since the quadratic jet vanishes on \(w=0\). Split off \(w\) by the
analytic Morse lemma with parameters \(u,v\). The shift solving
the \(w\)-critical equation has order at least two, so the binary
cubic of the remaining equation is still the original binary
cubic of the leading curve. Blow up this double point. The
exceptional line has possible singularities only at the three
roots of that cubic. They are simple, giving surface ordinary
double points by the \(U^2+s(V-V_0)\) quadratic jet. The
double-point blowup and these ordinary-double-point resolutions
extend the residue holomorphically, as in the preceding
quartic argument.

If the linear and quadratic jets both vanish, then \(R\) has
degree at most one in the affine variable \(z\). Indeed, the
homogenized term of \(x^i y^j z^k\) has local degree \(4-k\)
at this point in the chart with projective \(z=1\). Thus
\(R=a(x,y)z+b(x,y)\) with \(a\) a cubic whose binary degree-three
part has three distinct factors. The coefficient \(a\) is
squarefree, and the projective closure of \(a=0\) is smooth at
infinity because its intersection with the line at infinity
has three simple points. The residue proposition excludes
every affine singularity of this closure. It would therefore
be a smooth projective cubic. The affine residue of
\(dx\wedge dy/a\) extends to a nonzero holomorphic one-form on
that cubic by adjunction. Since it is the derivative of a
rational function, the rational function can have no pole,
and is constant on the compact curve. This is impossible.

When all singularities of a reduced leading quartic are of the
listed types, general \(t\) works simultaneously at the finitely
many boundary points. At the triple points the argument either
already gives a contradiction or extends the fibre form
holomorphically. At the other points the mild-boundary proof
gives either a local nonzero period or a holomorphic extension.
Thus the same global obstruction proves the first theorem.

## Transversality on repeated support

At a smooth point of the reduced support use local plane
coordinates \(u,v\), with the support given by \(u=0\).
If the component has multiplicity at least two, \(G_t\) has a
singularity there only when \(R_3\) vanishes. At such a point its
quadratic term has the form
\[
 c u^2+w(d u+e v)+f w^2 .
\]
Here \(e\ne0\) by the simple-zero hypothesis on \(R_3\). The
restriction to the \(v,w\) plane is a nondegenerate quadratic
form, so the germ has Hessian rank at least two. Points on the
singular locus of the reduced support are smooth on the projective
fibre by the nonzero \(R_3\) hypothesis. Simple smooth points
of the leading curve are smooth points of the projective fibre
by its tangential derivative.

All singularities of the compactification are therefore isolated.
The affine part is smooth, and the first derivatives at infinity
can all vanish only at the finitely many transverse intersections
just analyzed. Normality follows from the hypersurface depth
condition and regularity in codimension one. The rank-two or
rank-three double-point argument extends the residue form
holomorphically, contradicting the [global exactness obstruction](KELLER_QUARTIC_COMPONENT_BOUNDARY_OBSTRUCTION.md#residues-on-a-resolution).
This proves the second theorem.

Both exclusions apply to any quartic constant target-linear
combination of degree-five Keller components. In characteristic
zero the coefficient field embeds into \(\mathbb C\) after
restriction to a finitely generated subfield, and the geometric
hypotheses are preserved under algebraically closed base extension.

## Invertibility with a reduced quadratic coefficient

**Theorem QuinticKellerLinearSourceQuadraticCompletion.** Let \(F\)
have ordinary degree at most five and let a quartic constant target
component be \(R=a(x,y)z+b(x,y)\), with \(a\ne0\) squarefree of
degree at most two. Then \(F\) is invertible. In particular, this
applies when \(R_4=0\) has an ordinary quadruple point, the projective
linear and quadratic jets of the fibre equation vanish there,
and the resulting coefficient \(a\) in the direction of that
point is nonzero and squarefree.

It suffices to give a coordinate completion of \(R\) with inverse
degree at most eight. The
[bounded coordinate-completion theorem](LOW_DEGREE_COORDINATE_KELLER_MAPS_IN_THREE_VARIABLES.md)
then gives invertibility, since the resulting plane-map degree
is at most forty.

If \(a\) is constant, the graph completion \((x,y,R)\) has inverse
degree at most four. If \(a\) is linear or is a product of two lines,
the lines must be disjoint in the affine plane. In fact the
residue proposition excludes an affine intersection. After an
affine change write \(a=a(x)\), with \(d=\deg a\) equal to one
or two. At each root \(x=\xi\), the restriction \(b(\xi,y)\)
has nowhere vanishing derivative because \(da,db\) are independent
along this line. Thus that restriction is affine linear with
nonzero slope. Interpolation gives
\[
 b=f(x)y+g(x)+a(x)c(x,y),\qquad
 \deg f,\deg g<d,\qquad \deg c\le4-d,\qquad \gcd(f,a)=1 .
\]
Set \(Z=z+c\). Choose \(r,s\in\mathbb C[x]\) with
\(rf+sa=1\), \(\deg r<d\), and \(\deg s<\deg f\) when \(f\)
is nonconstant. If \(f\) is constant take \(r=1/f\), \(s=0\).
Then \((x,R,M)\) with \(M=-sy+rZ\) is an automorphism, with inverse
\[
 y=r(R-g)+aM,\qquad Z=s(R-g)+fM,\qquad z=Z-c(x,y).
\]
The degrees of \(y,Z\) in the inverse are at most \(d+1\) and
\(d\), respectively. The degree of \(z\) is at most six, using
\((4-d)(d+1)\le6\).

For irreducible quadratic \(a\), its projective conic is smooth.
If it meets the line at infinity at two points, an affine change
and scaling give \(a=xy-1\). The residue along the conic is, up
to a nonzero constant, \(dx/x\). It cannot have a rational
primitive, by its nonzero residue at \(x=0\) on the normalization.
This contradicts the proposition.

The remaining conic is tangent to the line at infinity. An affine
change and scaling give \(a=y-\phi(x)\) with \(\deg\phi=2\).
The restriction of \(b\) to this graph is again affine linear
with nonzero slope. Write \(b=\lambda x+\mu+ac(x,y)\), with
\(\lambda\ne0\), \(\deg c\le2\), and take \(Z=z+c\).
The completion \((R,a,Z)\) has inverse
\[
 x=(R-\mu-aZ)/\lambda,\qquad y=a+\phi(x),\qquad z=Z-c(x,y).
\]
The inverse degrees are at most two, four, and eight.

For the quadruple point application, use the projective chart
of the triple point proof. The vanishing jets imply \(R=az+b\),
and the leading plane quartic has no local cubic term at that
point. Therefore \(\deg a\le2\). Over a general characteristic-zero
field the same proof applies after complex embedding of a finitely
generated coefficient field. The formal inverse then descends as
in the bounded coordinate-completion theorem.
