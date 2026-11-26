# A Lefschetz obstruction for quartic-leading septic phases

Let \(K\) have characteristic zero. Consider a regular winding-one phase
on \(2xy+t^2=1\) of least sphere degree seven with normalized quadruple
leading root and the [boundary divisibilities](WINDING_ONE_CUBIC_LOGARITHMIC_BOUNDARY.md#a-quadruple-root).
Thus \(U=1+x\), \(v=Ut\), and
\[
 f=\frac{U}{x}(v^3+d(x)v^2+b(x)v+c(x))
   =\frac{U}{x}(Y^3+p(x)Y+q(x)),\qquad Y=v+d(x)/3.
\]
Regularity gives \(b(0)=-1\), \(c(0)=-d(0)\). The sphere functional
\(\Lambda\) takes the angular constant term and integrates in \(t\)
over \([-1,1]\) with measure \(dt/2\).

**Theorem DegreeFourSepticLefschetzExclusion.** Suppose the following
conditions hold over an algebraic closure of \(K\).

1. \(\deg p\le4\), \(\deg q=7\), and \(d(0)^2\ne1\).
2. Put \(u=x/(1+x)\), \(z=1/u\), \(p_-(z)=p(1/(z-1))\),
   \(q_-(z)=q(1/(z-1))\), and \(\Delta_-=-4p_-^3-27q_-^2\).
   Then \(\Delta_-(0)\ne0\) and \((p_-/\Delta_-)'(0)\ne0\).
3. On the chart \(u\ne0,1,\infty\), every critical point of
   \[
    F(u,Y)=\frac{Y^3+p(u/(1-u))Y+q(u/(1-u))}{u}
   \]
   is Morse, and their critical values are nonzero and pairwise distinct.

Then \(\Lambda(f^m)\) is nonzero for infinitely many positive \(m\).
The hypotheses admit varying depressed invariant, as shown in §5.

The nearest hyperelliptic result concerns a vanishing cycle of a Morse
point for \(Y^2+g(u)\). The family here is trigonal and the cycle in the
inverse-period condition is a lift of a circle about a fixed pole of
the profile. [Christopher and Mardešić, Theorems 1.5–1.6](https://arxiv.org/pdf/0803.2110).

## 1. The inverse-branch cycle and the finite model

Embed a finitely generated coefficient subfield into \(\mathbb C\).
Use the fibers
\[
 C_\zeta:\quad Y^3+p(u/(1-u))Y+q(u/(1-u))=\zeta u
\]
and the form \(\alpha=Ydu/u\). For large \(\zeta=w^3\), let \(\gamma\) be
the lift of a small positively oriented circle about \(u=1\) on the
branch with leading term \(u^{1/3}w\), choosing \(u^{1/3}=1\) at \(u=1\).
The inverse expansion converges on the fixed circle for large \(w\).
Here \(t\,dx/x=(Y-d/3)du/u\), the coefficient of \(w\) in \(\alpha\) is
regular at \(u=1\), and the depressed inverse has no constant term in
\(w\). Eventual purity would therefore imply the zero period germ
\(\int_\gamma\alpha\). [Prime-height inverse-period condition](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_INVERSE_PERIOD.md).

Normalize the graph closure in
\(\mathbb P^1_u\times\mathbb P^1_Y\times\mathbb A^1_\zeta\), giving a
proper surface \(S^\circ\) over the finite value line. Its boundary
charts are as follows.

- At \(u=0\) the cubic in \(v\) is \((v^2-1)(v+d(0))\), with three
  distinct roots. The derivative of the curve equation in \(Y\) is
  nonzero, giving three smooth sections. For finite \(u\ne1\) the monic
  curve equation excludes \(Y=\infty\).
- At \(u=1\) put \(u-1=\eta^3\) and \(Y=\eta^{-7}Z\). The equation is
  \[
   Z^3+\eta^{14}p(u/(1-u))Z+\eta^{21}(q(u/(1-u))-\zeta u)=0.
  \]
  The second coefficient is divisible by \(\eta^2\). The last term
  has a nonzero constant coefficient in \(\eta\) independent of finite
  \(\zeta\). Thus the three local sheets for \(Z\) are simple. For a
  primitive cube root \(\omega\), the action
  \((\eta,Z)\mapsto(\omega\eta,\omega^7 Z)\) permutes the sheets freely
  and transitively at \(\eta=0\). The quotient is the smooth normalized
  chart. Indeed \(\eta,Z\) are integral over the graph coordinates near
  the missing point, with \(Z\) a unit and \(1/Y=\eta^7/Z\). The
  quotient is then finite and birational to the graph closure. The
  pole of \(q\) excludes finite \(Y\) over \(u=1\). The result is a
  smooth section and a totally ramified point in each finite fiber
  for projection to the \(u\)-line.
- At \(u=\infty\), put \(P_-(z,Y)=Y^3+p_-(z)Y+q_-(z)\). The finite
  \(Y\) chart is the graph \(\zeta=zP_-(z,Y)\). Its only critical
  points at \(z=0\) are the three nodes \((0,y_i)\), with \(y_i\) the
  simple roots of \(P_-(0,Y)\). Their local coordinates are \(z,P_-\).
  At \(Y=\infty\) put \(h=1/Y\). The graph equation
  \(\zeta h^3=z(1+p_-(z)h^2+q_-(z)h^3)\) has unit derivative in \(z\)
  at \(h=z=0\), giving smooth coordinates \((h,\zeta)\). For
  \(\zeta\ne0\) the unique point over \(u=\infty\) is totally ramified
  of degree three. For \(\zeta=0\) the component \(z=0\) is a projective
  line meeting the rest of the fiber in the three nodes.

The total surface is smooth and rational. Its finite singular fibers
have the three nodes over zero and the interior nodes of hypothesis 3.
The possible poles of \(\alpha\) on smooth fibers lie on marked
sections over \(u=0,1,\infty\). These sections extend across the finite
critical values away from the nodes. In particular over \(u=\infty\)
the marked section is \(h=z=0\), away from the three finite roots.

For large regular \(\zeta\), projection from \(C_\zeta\) to the
\(u\)-sphere is finite of degree three. Both sides of the small
circle about \(u=1\) have connected preimages, since every component
maps onto its side and each side contains a unique totally ramified
point. The inner point is over \(u=1\) and the outer one is over
\(u=\infty\). The lifted circle is a union of three disjoint circles
by the uniform large-\(w\) expansion. Cutting along just \(\gamma\)
leaves a connected surface because the other two circles still join
both sides. Thus \([\gamma]\ne0\) in compact \(H_1(C_\zeta;\mathbb Q)\).

## 2. The vanishing cycles span compact homology

**Lemma FiniteNodesDetectInversePeriod.** Consider the same trigonal
pencil and inverse-branch period, with the conditions on the cubic
at \(u=0,\infty\) and on the interior critical points in
Theorem DegreeFourSepticLefschetzExclusion. Suppose that the
normalized finite model over \(u=1\) is smooth, without critical
points, and the poles there lie on disjoint sections over the
finite value line. If the preimages of both sides of the inverse
circle are connected for a large regular value, the inverse-branch
period is nonzero as a germ.

The charts at \(u=0,\infty\), the three nodes at value zero, and
the description of the three lifted circles depend only on the
stated shared conditions. Under the lemma's hypotheses the same
gluing shows that the selected circle has nonzero compact homology.
The following spanning and period arguments prove the lemma and
complete the announced dominant constant-term case.

Compactify \(S^\circ\) to a smooth projective rational surface \(S\)
with simple normal-crossing boundary \(D\), by compactifying the value
line and resolving over infinity. The restriction of
\(H^1(S^\circ;\mathbb Q)\) to a compact smooth finite fiber is zero.
Indeed the support sequence gives
\[
 0=H^1(S;\mathbb Q)\longrightarrow H^1(S^\circ;\mathbb Q)
   \longrightarrow \bigoplus_i\mathbb Q[D_i]
   \longrightarrow H^2(S;\mathbb Q).
\]
The degree-two group with support in \(D\) is the sum of the
groups for its irreducible components, by the Thom isomorphism and
because their intersections have real codimension four. The last
map sends each generator to its divisor class. On the rational
surface \(S\), the divisor class map from
\(\operatorname{Pic}(S)\otimes\mathbb Q\) is injective, by
\(H^1(\mathcal O_S)=0\) and the exponential sequence. After clearing
denominators and torsion, a relation in the kernel is the divisor
of a rational function with zeros and poles only in \(D\). Its
logarithmic differential divided by \(2\pi i\) realizes the
corresponding meridian values. These classes span
\(H^1(S^\circ;\mathbb Q)\), and every such function is invertible
regular on \(S^\circ\) and constant on the compact fiber. This proves
the zero restriction, so the fiber also has zero image in
\(H_1(S^\circ;\mathbb Q)\).

The kernel of inclusion on \(H_1\) is the span of the finite vanishing
cycles transported to a common smooth reference fiber. To use the
Lefschetz handle description, take a disk containing the finite
critical values with the reference value on its boundary. Outside
the disk the proper smooth fibration adds a collar. Over the disk,
the description starts with a regular fiber times a disk and
attaches one two-handle per node along its vanishing cycle. For the
three nodes over zero, use disjoint Morse balls. Small smooth shifts
of the values supported in the balls make those values distinct,
without new critical points on the cutoff annuli. The attaching
annuli agree in a nearby regular fiber outside the modified value
cluster. Therefore the same kernel formula applies to the three
cycles over zero. The zero image proved above now makes the kernel
all of the compact \(H_1\). The handle description and positive Dehn
twist in the distinct-value case are given in
[Apostolakis, Piergallini and Zuddas, §6 and Proposition 9](https://arxiv.org/pdf/1104.4536).

## 3. Testing periods at the finite nodes

At an interior node a vanishing cycle \(\delta\) has a nonzero period
germ \(\int_\delta\alpha\). Its derivative in \(\zeta\) is the
Gelfand–Leray integral of \(d\alpha=dY\wedge du/u\). This follows by
differentiating along a lift of the parameter vector field, since
the exact term in Cartan's formula integrates to zero. The two-form
is nowhere zero at the node. In local Morse coordinates with
\(\zeta-\zeta_*=AB\), write it as a unit times \(dA\wedge dB\).
On \(B=(\zeta-\zeta_*)/A\), the relative form is, up to sign, that
unit times \(dA/A\). Its integral on the annular core tends to a
nonzero multiple of the value of the unit at the node.

Transport the putative zero period of \(\gamma\) to a near regular
value and around an interior critical value. The local monodromy is
a positive Dehn twist along \(\delta\). Its variation on homology is
\(\langle\gamma,\delta\rangle\delta\), up to the common convention
for intersection signs. This follows by counting the oriented arcs
crossing the twisting annulus. It applies in punctured homology as
well, since the marked sections miss that annulus. The zero germ
stays zero under continuation. The node period is nonzero, so
\(\gamma\) has zero compact intersection with every such \(\delta\).

At \(\zeta=0\) fix a small circle \(|z|=\epsilon\) and let
\(Y_i(z,w)\) solve \(P_-(z,Y_i)=w\) analytically near \((0,0)\) with
initial value \(y_i\). For small \(\zeta\ne0\) the lifts with
\(w=\zeta/z\) represent the three vanishing cycles \(\delta_i\) of
the nodes \(\zeta=zP_-\). Orient them by the common positive
direction in \(z\). Their sum is zero in compact homology, since
they bound the part of the smoothed fiber near the rational
component. The monodromy is the product of the three disjoint
positive twists. Its period variation is, up to sign,
\(\sum_i a_i I_i(\zeta)\), where
\(a_i=\langle\gamma,\delta_i\rangle\), \(\sum_i a_i=0\), and
\(I_i(\zeta)=\int_{\delta_i}\alpha\).

## 4. Independence at the three nodes over zero

The integrals \(I_i\) extend analytically to zero along the fixed
circles. Substituting \(\alpha=-Y_i(z,\zeta/z)dz/z\) gives
\[
 \frac{I_i(0)}{-2\pi i}=y_i,\qquad
 \frac{I_i'(0)}{-2\pi i}=
 B_i:=\left.\frac{d}{dz}\frac{1}{P_{-,Y}(z,Y_i(z,0))}\right|_{z=0}.
\]
Both vectors have coordinate sum zero, by summing the inverse
branches of the depressed cubic and differentiating. The quadratic
coefficient of the polynomial of degree at most two interpolating
\(B(y_i)=B_i\) is, by Lagrange interpolation,
\[
 \begin{aligned}
 \sum_i B_i/P_{-,Y}(0,y_i)
 &=\frac12\left.\frac{d}{dz}\sum_i P_{-,Y}(z,Y_i(z,0))^{-2}\right|_{z=0}\\
 &=-3\left.\frac{d}{dz}\left(\frac{p_-(z)}{\Delta_-(z)}\right)\right|_{z=0}.
 \end{aligned}
\]
Multiplying each inverse squared derivative by the discriminant
leaves the squared difference of the other two roots. The sum of
the three squared differences is \(-6p_-(z)\), giving the last
equality. By hypothesis 2 this coefficient is nonzero, so
\((y_i)\), \((B_i)\) are independent in the trace-zero plane. Together
with the constant vector they form a basis of \(\mathbb C^3\).
Hence \(\sum_i a_i I_i=0\) and \(\sum_i a_i=0\) force every \(a_i=0\).

Thus the compact class of \(\gamma\) pairs to zero with every finite
vanishing cycle. These span compact \(H_1\) by §2, so nondegeneracy
of the intersection form forces \([\gamma]=0\), contrary to §1.
This rules out eventually vanishing pure moments.

## 5. An exact varying-invariant example

Take \(d=0\), \(b=-(1+x)\), \(c=x^7\). It is regular and has sphere
degree seven by the [degree formula](MIXED_CUBIC_ZERO_CENTER_AND_CONDITIONAL_OFFCENTER.md).
Here \(p=b\), \(q=c\), with nonconstant ratio \(q^2/p^3\). At \(z=0\),
\(p_-=0\), \(q_-=-1\), \(p_-'=1\), so the two nonvanishing
conditions in hypothesis 2 hold.

For the interior critical points, the critical equation in \(Y\) gives
\(x=3Y^2-1\). Put \(D=3Y^2-1\) and restrict \(F\) to this curve. Its
value and derivative are
\[
  R(Y)=3Y^2D^6-\frac{6Y^5}{D},\qquad
  R'(Y)=\frac{6Y}{D^2}H(Y),\qquad
  H(Y)=D^7(21Y^2-1)-Y^3(9Y^2-5).
\]
The interior chart excludes \(YD=0\). The polynomial \(H\) has degree
16 and has no root with \(YD=0\), by substitution. At such a critical
point \(F_{YY}=6Y/u\ne0\). Since \(dx/dY=6Y\ne0\) on the critical
curve, a simple root of \(H\) is exactly a nondegenerate critical
point of \(F\). Indeed eliminating the \(Y\)-critical equation by
the implicit function theorem gives the Schur complement of
\(F_{YY}\) in the Hessian. It remains to check that the roots are
simple and the corresponding critical values are distinct and
nonzero.

The following finite arithmetic certifies these assertions. Reduce
modulo 19, where \(H\) retains degree 16. For \(H,H'\) the degrees of
the successive nonzero Euclidean divisors, starting with \(H'\), are
\[
 15,14,13,12,8,7,6,5,4,3,2,1,0.
\]
In \(\mathbb F_{19}[Y]/(H)\), the class of \(R\) in increasing
degree order has coefficients
\[
 (18,0,4,14,3,0,18,0,18,0,9,0,18,0,16).
\]
To reproduce this row, multiply it by \(D\) and reduce modulo \(H\),
giving the reduction of \(3Y^2D^7-6Y^5\). Multiplication by \(D\) is
invertible as \(H\) and \(D\) have no common root modulo 19. The traces
of multiplication by the first 16 positive powers of this class are,
in order,
\[
 (14,18,7,12,0,10,18,10,5,6,5,10,11,2,8,4).
\]
Each trace is computed by summing, for \(0\le j<16\), the coefficient
of \(Y^j\) in the remainder of \(Y^jR^k\) modulo \(H\). Newton's
identities give the monic characteristic polynomial \(L\) in
increasing degree order as
\[
 (3,18,16,14,0,13,15,17,1,5,12,6,13,2,13,5,1).
\]
For \(L,L'\) the Euclidean divisor degrees are
\(15,14,13,12,11,10,9,8,7,6,5,4,3,2,1,0\). Since the leading
coefficient of \(H\) and the resultant of \(H,D\) are units modulo
19, the multiplication matrix is defined over \(\mathbb Z_{(19)}\)
and reduces to this calculation. Squarefreeness of \(H\),
squarefreeness of \(L\), and \(L(0)\ne0\) consequently hold in
characteristic zero. Over the algebraic closure each root of \(H\)
contributes its critical value as one eigenvalue of multiplication by
\(R\). Thus hypothesis 3 holds exactly for this example.
