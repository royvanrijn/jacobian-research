# Generic exclusion of quartic-leading winding-one septics

Use the normalized quadruple-root phase
\[
 f=\frac{1+x}{x}(v^3+d(x)v^2+b(x)v+c(x)),\qquad v=(1+x)t,
\]
on \(t^2+2xy=1\), with regularity \(b(0)=-1\), \(c(0)=-d(0)\)
and bounds \(\deg d\le3\), \(\deg b\le5\), \(\deg c\le7\).
Set \(Y=v+d/3\), and write the depressed cubic as \(Y^3+p(x)Y+q(x)\).
Use \(d_3,b_5,c_7\) for the three coefficients of maximal allowed
degree in \(d,b,c\).

**Theorem HighPoleLefschetzExclusion.** The
[Lefschetz exclusion](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_LEFSCHETZ_EXCLUSION.md)
continues to hold if its degree conditions on \(p,q\) are
replaced by either of the following conditions.

1. \(\deg d\le2\) and \(\deg p=5\).
2. \(d_3\ne0\) and \(b_5^2-4d_3c_7\ne0\).

The other hypotheses of that theorem retain their stated form.
In particular they require \(d(0)^2\ne1\), the two nonvanishing
conditions at \(x=-1\), and distinct nonzero interior Morse
critical values. Under these hypotheses, infinitely many positive
pure moments are nonzero.

## The finite model at the pole of the profile

Work over \(\mathbb C\) by embedding the finitely generated
coefficient field. Use \(u=x/(1+x)\) and the value parameter
\(\zeta\) of the Lefschetz exclusion. Its proof at \(u=0,\infty\)
is unchanged. At \(u=1\) put \(\tau=u-1\).

If \(\deg p=5\) and \(\deg d\le2\), then \(\deg q\le7\).
Write \(A(\tau)=\tau^5p\) and \(B(\tau)=\tau^7q\) in this local
coordinate, with \(A(0)\ne0\). After \(\tau=\eta^2\) and
\(Y=\eta^{-5}T\), the equation is
\[
 T^3+A(\eta^2)T+\eta B(\eta^2)=\zeta(1+\eta^2)\eta^{15}.
\]
At \(\eta=0\) the roots are simple. The involution
\((\eta,T)\mapsto(-\eta,-T)\) exchanges the two nonzero-root
sheets and fixes the zero-root sheet. The first two sheets give
one smooth quotient chart, with projection \(\tau=\eta^2\).
On the zero-root sheet \(T\) is odd as a function of \(\eta\),
so the quotient is itself a smooth chart with coordinate \(\tau\)
and \(Y\) meromorphic in \(\tau\). This second chart allows finite
\(Y\) when the degree of \(q\) is at most five. The local model
has two disjoint sections over the finite value line.

If \(d_3\ne0\), write \(D(\tau)=\tau^3d\),
\(B(\tau)=\tau^5b\), \(C(\tau)=\tau^7c\). In the undepressed
coordinate \(W=\tau^3v\), the equation is
\[
 W^3+D(\tau)W^2+\tau B(\tau)W+\tau^2 C(\tau)
     =\zeta(1+\tau)\tau^9.
\]
The root \(-D(0)\) is simple at \(\tau=0\). For the two other
branches put \(W=\tau T\). The equation becomes
\[
 \tau T^3+D(\tau)T^2+B(\tau)T+C(\tau)
     =\zeta(1+\tau)\tau^7.
\]
The quadratic at \(\tau=0\) has discriminant
\(b_5^2-4d_3c_7\), since \(x=-(1+\tau)/\tau\).
It has two simple roots. Thus there are three smooth sections
over finite \(\zeta\), with projection unramified to \(\tau\)
on each.

These local analytic sections, with their indicated ramification,
describe the normalization of the graph near \(u=1\). For
\(\tau\ne0\) they cover all three cubic roots. Each chart is
finite over a neighborhood in \((\tau,\zeta)\), and its
meromorphic root defines an analytic map to the projective
\(Y\)-line. The map to the graph is proper and quasi-finite
locally, hence finite. Its restriction off \(\tau=0\) is
birational onto the corresponding graph germ. This identifies
the union of the charts with the normalization there. The poles
of \(Y\,du/u\) over \(u=1\) lie on these disjoint sections,
which are smooth and avoid all finite nodes.

## Connectivity of the inner part of a large fiber

Let \(D_1\) be a small disk centered at \(u=1\). The proof
of the Lefschetz exclusion needs the preimage of \(D_1\) in a
large regular fiber to be connected. The preimage of its
complement is connected, as before, by the unique totally
ramified point over \(u=\infty\).

In the degree-five \(p\) case scale
\(\tau=\zeta^{-2/15}s\), \(Y=\zeta^{1/3}Z\) along large positive
real \(\zeta\). The limiting curve on \(s\ne0\) is
\[
 Z^3+A(0)s^{-5}Z=1.
\]
It is irreducible. Over \(\mathbb C(Z)\) the equation for \(s\)
is \(s^5=A(0)Z/(1-Z^3)\), with a simple zero at \(Z=0\), so
Eisenstein's criterion at that zero applies. This gives
irreducibility of the curve and thus a transitive monodromy
action on the three roots over the punctured \(s\)-line.

In the degree-three \(d\) case scale \(\tau=\zeta^{-1/9}s\)
in the equation for \(W\). The limiting curve is
\[
 W^2(W+D(0))=s^9.
\]
Again it is irreducible, by the simple zero at \(W=-D(0)\) in
the equation for \(s\). Its degree-three projection to \(s\)
has transitive monodromy.

In each case choose finitely many based loops away from zero
and the limit branch values whose permutations generate a
transitive subgroup. The coefficients of the scaled equation
converge uniformly on their compact images. Since the cubic
discriminant is nonzero there, the permutations persist under
this perturbation. For all sufficiently large positive \(\zeta\)
these loops are inside \(D_1\) in the original coordinate.
Every connected component of the preimage of \(D_1\) maps
onto it by properness and openness of a finite map of curves.
The transitive permutations therefore put all components in
one. The fixed circle itself has three distinct lifted
circles by the large-value inverse expansion. Since its two
sides have connected preimages, the selected lift is nonzero
in compact homology, exactly as in the original Lefschetz proof.

Lemma FiniteNodesDetectInversePeriod in that proof now applies.
The normalized total surface is rational, the finite nodes are
unchanged, and the pole sections stay away from the twisting
annuli. The lemma rules out the zero inverse period and proves
the theorem.

## Existence of a dense open exclusion

**Corollary GenericQuarticLeadingSepticExclusion.** A nonempty
Zariski-open dense subset of the normalized regular quartic-leading
septic coefficient space has infinitely many nonzero positive
pure moments. This also holds in each of the three coefficient
strata consisting of \(d_3\ne0,\ b_5^2-4d_3c_7\ne0\), or
\(d_3=0,\ b_5\ne0\), or \(d_3=b_5=0,\ c_7\ne0\).

Indeed least degree seven in this normalized family means
\((d_3,b_5,c_7)\ne(0,0,0)\), by the
[winding-one degree formula](MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md).
In the third stratum \(\deg p\le4,\ \deg q=7\) and the
original Lefschetz result applies. In the second stratum
\(\deg p=5\), giving the first case of the theorem. The first
stratum is its second case and is dense in the entire septic
coefficient space.

The two nonvanishing conditions at \(x=-1\) and \(d(0)^2\ne1\)
define a nonempty open in each stratum. For example take
\(d(0)=0\) and choose \(b_1,b_2\) to prescribe
\(p(-1)=0\), \(p'(-1)\ne0\). Then choose \(c_1\) so that
\(q(-1)=-1\). These changes leave the three top coefficients
fixed. Since \(x=1/(z-1)\) has nonzero derivative at \(z=0\),
they ensure \(\Delta_-(0)\ne0\) and \((p_-/\Delta_-)'(0)\ne0\).

For the interior critical points use the coordinate \((x,v)\)
with \(x\ne0,-1\). On this domain each prescribed second jet
of \(f\) is obtainable by varying the coefficients of
\(x,x^2,x^3\) in \(c\), the coefficients of \(x,x^2\) in \(b\),
and the coefficient of \(x\) in \(d\). Their variations are
\(U\) times the six polynomials \(1,x,x^2,v,xv,v^2\),
which span the second jets since \(U\ne0\). Thus the incidence
of a degenerate critical point imposes codimension three at
each fixed point, from two first derivatives and a singular
Hessian. The incidence of a zero critical value has the same
codimension. Allowing the point to vary adds dimension two,
so the closures of the projected bad sets have codimension at
least one in the coefficient stratum.

For two interior points with distinct \(x\), varying the
coefficients of \(x,\ldots,x^4\) in \(c\) and of \(x,x^2\) in
\(b\) gives arbitrary first jets at both points. First prescribe
the two \(v\)-derivatives using the \(b\)-variation and then
the two remaining values and \(x\)-derivatives using Hermite
interpolation in the \(c\)-variation. The equal-critical-value
incidence has codimension five at fixed points, so its projection
also has closure of codimension at least one. Two distinct
points with the same \(x\) cannot share a critical value, since
that would give two distinct double roots of a cubic in \(v\).
All these variations preserve the top coefficients and the
regularity constraints. The coefficient strata are irreducible.
Intersecting the complements of the three projected closures
with the boundary nonvanishing conditions gives the claimed
nonempty dense opens.
