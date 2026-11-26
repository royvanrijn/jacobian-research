# A Lefschetz obstruction at a second-order repeated center

Use the regular quartic-leading winding-one septic setup of
[DegreeFourSepticRepeatedCenterExclusion](WINDING_ONE_CUBIC_REPEATED_CENTER_LEFSCHETZ_EXCLUSION.md).
Thus \(u=x/(1+x)\), the depressed pencil is
\(P(u,Y)=\zeta u\) with \(P(u,Y)=Y^3+p(u/(1-u))Y+q(u/(1-u))\),
and the form used for the inverse period is \(\alpha=Ydu/u\).
Assume \(d(0)^2=1\). Let \(y_0=-2d(0)/3\) be the double root of
\(P(0,Y)\), let \(h(0)=y_0\) on the analytic branch
\(P_Y(u,h(u))=0\), and let \(\phi(u)=P(u,h(u))/u\) and
\(\zeta_c=\phi(0)\). The initial zero of the numerator gives the
extension of \(\phi\) to \(u=0\).

**Theorem DegreeFourSepticSecondOrderCenterExclusion.** Suppose
\(d(0)^2=1\), use one of the
[smooth pole models of the first-order obstruction](WINDING_ONE_CUBIC_REPEATED_CENTER_LEFSCHETZ_EXCLUSION.md),
and let hypotheses 2 and 3
of [DegreeFourSepticLefschetzExclusion](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_LEFSCHETZ_EXCLUSION.md)
hold. Suppose \(\phi'(0)=0\), \(\phi''(0)\ne0\), and \(\zeta_c\)
is nonzero and distinct from all interior critical values. Then
the sphere moments \(\Lambda(f^m)\) are nonzero for infinitely many
positive orders. The conditions admit a varying-invariant example.

## 1. The scaled elliptic piece

Work over \(\mathbb C\) after embedding a finitely generated field
of coefficients. As in the first-order center, the local exact
equation in an analytic coordinate \(w\) is
\(w^2=u(\zeta-\phi(u))\). Since \(\phi(u)-\zeta_c\) vanishes to
order exactly two, take a square root \(r\) of it with a simple
zero at zero. Write \(s=\zeta-\zeta_c\) and rescale \(w\) by the
unit \(\sqrt{u/r}\) to obtain the local analytic equation
\[
                         w_1^2=r(s-r^2).
\]
On the local double-cover piece near \(y_0\), the inverse analytic
coordinate change has expansion
\(Y=h(u)+g(u)w_1+O(w_1^2)\) with \(g(0)\ne0\).

For small nonzero \(s\), let \(L_s\) be the part of this double
cover over the closed disk \(|r|\le\epsilon\). It has three simple
branch points at \(r=0,\pm\sqrt{s}\). The boundary is connected,
since a loop about the disk acts by transposition on the two
sheets, and the Euler characteristic is \(2-3=-1\). Thus \(L_s\)
is a genus-one surface with one boundary circle. It has the single
marked point \(r=w_1=0\). The other sheet near the simple root of
\(P(0,Y)\) lies outside \(L_s\). The complement of the interior of
\(L_s\) in the compact smooth fiber is connected and forms a smooth
family with boundary and marked sections across \(s=0\), since the
only singularity at that value is the center. It trivializes over
the parameter disk. The two compact pieces meet in one separating
circle, so compact first homology splits as the sum of their first
homology groups. On punctured homology, Mayer–Vietoris shows that
the change of any class under monodromy around \(s=0\) is represented
by a class from \(L_s\) minus its marked point. On the compact
homology of \(L_s\) this change is \((T-1)\) applied to the local
component of the class, where \(T\) is the local monodromy.

On writing \(s=\tau^4\), \(r=\tau^2 Z\) and \(w_1=\tau^3 W\), one
gets the exact scaled curve
\[
                         E:\quad W^2=Z(1-Z^2).
\]
For \(\tau\) small, any fixed compact cycle in this curve avoiding
\(Z=0\) and infinity lies in the scaled local piece. Capping the
one boundary of \(L_s\) at infinity identifies its compact-piece
homology with that of the smooth projective elliptic curve \(E\).
A circuit of \(s\) induces \((Z,W)\mapsto(-Z,\pm iW)\) on this
homology, by the change \(\tau\mapsto i\tau\). Its square is the
elliptic involution, which acts by minus one on first homology.
Indeed it changes the sign of the holomorphic differential and
hence of both parts of \(H^1(E;\mathbb C)\). In particular \(T-1\)
is invertible and \(T\) acts irreducibly over \(\mathbb Q\) on
this two-dimensional space.

## 2. Detection of rational local cycles by the period

On any such fixed cycle the expansion is
\[
 \alpha=h(u)\frac{du}{u}+\tau^3g(0)\frac{W\,dZ}{Z}+O(\tau^5).
\]
This follows from \(du/u=(1+O(\tau^2))\,dZ/Z\) and the local
inverse-coordinate expansion. The first term has constant periods
on the scaled cycles, since it differs from \(y_0\,du/u\) by an
analytic exact form. In particular the marked-point meridian
contributes no \(\tau^3\) term. The selected differential on \(E\)
satisfies
\[
                  \frac{W\,dZ}{Z}=\frac23\,dW+\frac23\frac{dZ}{W}.
\]
The second term is a nonzero holomorphic differential, and its
period functional on rational homology has zero kernel. To see
this, the differential is an eigenvector under the order-four
monodromy and its kernel is therefore a rational \(T\)-stable
subspace. It is proper because a holomorphic differential with
all periods zero on a compact curve would have a global
holomorphic primitive. Irreducibility then forces that kernel
to be zero.

If the inverse period of \(\gamma\) were identically zero, its
monodromy variation around \(s=0\) would be zero. Represent this
variation in local punctured homology as in §1 and take its
\(\tau^3\) coefficient. The preceding calculation forces its
compact local class to vanish. As \(T-1\) is invertible, the
local component of \([\gamma]\) is zero. Thus \([\gamma]\) has
zero compact intersection with the entire local homology.

## 3. The compact-fiber proof after local resolution

The surface \(w_1^2=r(s-r^2)\) is an ordinary double point at the
origin. Blow it up. On the \(r\)-chart \(w_1=rW_1,s=rS\), the
equation yields \(S=W_1^2+r\), so the surface is smooth and the
map near the only critical point in this chart is
\[
 s=r(W_1^2+r)=(r+W_1^2/2)^2-W_1^4/4.
\]
The exceptional conic outside this chart is covered by the
\(s\)-chart \(r=sR,w_1=sW_1\) with equation \(W_1^2=R-sR^3\).
The marked section is \(R=W_1=0\), where the map is smooth.
Thus the rational finite model after this blowup is smooth, with
one isolated \(A_3\) critical point over \(\zeta_c\).

For the compact-homology spanning argument, perturb the local map
\(A^2-W_1^4/4\), \(A=r+W_1^2/2\), to
\(A^2-W_1^4/4+\varepsilon W_1\), using a cutoff in a small
neighborhood of the critical point. In the inner chart the three
critical points have \(A=0\), \(W_1^3=\varepsilon\), and distinct
critical values \(3\varepsilon W_1/4\). They are Morse. For sufficiently
small \(\varepsilon\), submersion persists on the cutoff annulus and
no other critical point is affected. This perturbation changes
neither the total topological space nor a reference fiber outside
the small value disk. Its three new vanishing cycles are supported
in the local piece, because the perturbation is inside that piece
and the complement remains trivial over the disk. The same
two-handle description applies to the perturbed map, so those
cycles together with the other finite vanishing cycles span compact
first homology by the
[rational-surface inclusion argument](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_LEFSCHETZ_EXCLUSION.md#2-the-vanishing-cycles-span-compact-homology).
[Apostolakis, Piergallini and Zuddas, §6 and Proposition 9](https://arxiv.org/pdf/1104.4536).

The models over \(u=1,\infty\), the ordinary critical values and the
simple center root are unchanged. The original branch circle is
therefore still nonseparating. If the pure moments vanished
eventually, the [inverse-period condition](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_INVERSE_PERIOD.md)
would make its period vanish. The original arguments at the
ordinary critical values and the three nodes over zero give
zero compact intersection with their cycles. Section 2 gives
zero compact intersection with all cycles supported in the
center piece, including the three cycles used in the spanning
argument. This contradicts the nonzero class of the branch circle
and the nondegenerate compact intersection form.

## 4. An exact example on the second-order locus

Take \(d=1\), \(p=-4(1+x)/3\), and
\(q=x^7-16/27-10x^2/9\). In the original cubic this gives
\(b=-1-4x/3\), \(c=x^7-1-4x/9-10x^2/9\). In particular the phase
is regular of least sphere degree seven and the ratio \(q^2/p^3\)
varies. At \(z=1/u=0\), one has \(p_-=0\), \(q_-=-73/27\), and
\(p_-'=4/3\), giving both required nonvanishing conditions at that
boundary. On the center branch \(h=-2(1+x)^{1/2}/3\) one finds
\[
 P(u,h)=\frac89x-\frac89x^2-\frac1{27}x^3+O(x^4),\qquad
 \phi(u)=\frac89-\frac{25}{27}u^2+O(u^3).
\]

As in the first-order example, write \(V=3Y/2\), \(D=V^2-1\) on
the critical branch in \(Y\), so that \(x=D\), \(u=D/V^2\). The
critical value function here is
\[
 R(V)=V^2D^6-\frac{16}{27}\frac{V^2(V^2-V+1)}{V-1}
                       -\frac{10}{9}V^2D.
\]
Its derivative is \(2V(V+1)H(V)/(27(V-1)^2)\), where
\[
 H(V)=\frac{27(V-1)^2 D^5(7V^2-1)-8(3V^3-6V^2+4V-2)
                   -30(V-1)^2(2V^2-1)}{V+1}.
\]
This polynomial has degree 13 and is nonzero at \(V=0,\pm1\), as
also verified by the reduction below. Its interior simple roots
give Morse points because \(P_{YY}\ne0\) and \(du/dV\ne0\) there.

Over \(\mathbb F_{29}\), \(H\) retains degree 13 and has values
\(H(0)=15\), \(H(1)=4\), and \(H(-1)=23\). For \(H,H'\) the
nonzero Euclidean divisor degrees from the derivative onwards are
\((12,11,10,9,6,5,4,3,2,1,0)\). In the quotient by \(H\), the
class of \(R\) has coefficients in increasing degree order
\[
 (27,6,16,1,3,25,23,2,26,0,27,23,28).
\]
This row is checked by multiplication by \(27(V-1)\) and reduction
modulo \(H\). For multiplication by its first 13 positive powers,
sum the coefficients of \(V^j\) in the remainders of \(V^jR^k\)
over \(0\le j<13\). The traces are
\((0,26,19,28,22,20,20,1,7,28,20,0,19)\). Newton's identities give
the monic characteristic polynomial \(L\) in increasing degree order
\[
 (24,28,6,26,20,5,5,9,18,5,13,16,0,1).
\]
For \(L,L'\) the nonzero divisor degrees are
\((12,11,10,9,8,7,5,4,3,2,1,0)\), and \(L(0)=24\), \(L(8/9)=8\).
The preserved degree and invertibility of \(27(V-1)\) in the
localized quotient allow all these nonvanishing tests to lift to
characteristic zero. After splitting \(H\), the eigenvalues of its
quotient multiplication by \(R\) are exactly the interior critical
values. Thus these values are distinct and avoid zero and
\(\zeta_c=8/9\), giving the stated example.
