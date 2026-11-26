# Tangent low-branch collision at the high pole

Keep the notation and the hypotheses of the
[low-branch collision exclusion](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_LOW_BRANCH_COLLISION.md)
at \(u=0,\infty\) and at the interior critical points. In particular,
the normalized winding-one cubic is of exact least sphere degree
seven, with quartic leading coefficient and \(d_3\ne0\). Suppose
\(k\ge7\). The centered low quadratic at \(u=1\) has the form
\(S^2=\tau^7K(\tau,\zeta)\), with \(\tau=u-1\), and \(\zeta_*\) is the
unique zero of \(K(0,\zeta)\).

**Theorem TangentLowBranchCollisionExclusion.** Write
\[
 K(\tau,\zeta_*)=b\tau^2+c\tau^3+O(\tau^4),\qquad
 K_\zeta(\tau,\zeta_*)=a_0+a_1\tau+O(\tau^2),
 \qquad a_0\ne0.
\]
If
\[
 b\ne0,\qquad \kappa=\frac{9a_1}{4a_0}-\frac{7c}{4b}-1\ne0,
\]
then infinitely many positive pure sphere moments are nonzero. The
assertion includes coincidence of \(\zeta_*\) with zero or with one
interior critical value.

## The local smoothing

Work over \(\mathbb C\) by embedding a finitely generated coefficient
field as in the low-branch exclusion. The normalized chart
\(L^2=\tau K\) resolves in the chart \(L=\tau h\) with
\(\tau h^2=K\). Put \(w=\zeta-\zeta_*\) and write
\(K=K_0(\tau)+w a(\tau,w)\), where \(K_0=K(\tau,\zeta_*)\) and
\(a(\tau,0)=K_\zeta(\tau,\zeta_*)\). The analytic coordinates
\[
 \xi=\tau\sqrt{\frac{K_0(\tau)}{b\tau^2}\frac{a_0}{a(\tau,w)}},
 \qquad H=h\sqrt{\frac{\tau}{a(\tau,w)\xi}},
 \qquad A=-\frac{b}{a_0}
\]
give the exact equation \(w=\xi(H^2+A\xi)\). Both square roots
are units, and the coordinate change is invertible near the critical
point. Completing the square in \(\xi\) gives
\(A(\xi+H^2/(2A))^2-H^4/(4A)\), the isolated \(A_3\) model.
The rational exceptional curve is \(\xi=0\) on the special fiber.
It meets the strict fiber only at this point. The marked section
meets the exceptional curve away from it, where the projection is
smooth, in the other resolution chart described in the low-branch
exclusion.

With \(w=s^4\), \(\xi=s^2 Z\), and \(H=sW\), the compact core of the
local smoothing is the elliptic curve
\[
             Z(W^2+AZ)=1,\qquad R=ZW,\qquad R^2=Z-AZ^3 .
\]
The local Milnor fiber is this curve with small disks removed at
\(Z=0,\infty\). Indeed, near the critical point both solutions for
\(\xi\) at a fixed small \(H\) remain small, and the two points
removed are the simple poles of \(W\), beyond \(|W|=\epsilon/|s|\).
The attachment charts use \(t=s/H\), \(\xi=H^2 J(t)\), where
\(J(1+A J)=t^4\). The branches \(J(0)=0,-1/A\) are simple. They
attach to the exceptional and strict components respectively, and
the equation \(s=Ht\) gives one positive twist along each attachment
circle after the base change. The explicit elliptic scaling
trivializes the core on the \(s\)-disk, and away from the two
attachments the local fibration has no further variation on that
disk. The base monodromy induces \((Z,W)\mapsto(-Z,-iW)\), or its
inverse, on elliptic homology. Its square sends \(R\) to \(-R\) at
fixed \(Z\) and hence acts by \(-1\).

The attachment to the exceptional curve bounds the smoothed piece
of that curve in the compact fiber. The piece is a disk, since the
exceptional curve has one attachment and no other singular point.
The other attachment circle is opposite to it in the homology of
the two-boundary local fiber. Both circles are consequently zero
in compact homology. In the punctured fiber, the first circle is
peripheral around the marked section, and the second has the
opposite peripheral class. Each pairs to zero with every closed
cycle. A Dehn twist changes a homology class by its intersection
with the twist circle times that circle. Thus the two twists act
trivially on the punctured homology as well. This proves that the
fourth power of monodromy supported in this local neighborhood
acts trivially on the whole punctured first homology.

## Two period coefficients

On the resolved chart the form \(\alpha=Y\,d\tau/(1+\tau)\) is a
common term \(F(\tau,w)d\tau\) plus
\(\tau^2 h\,d\tau/(1+\tau)\), as in the low-branch exclusion. Here
\(F\) is analytic in \(w\) and meromorphic in \(\tau\) with a fixed
finite polar order. Let
\(C_1=(c/b-a_1/a_0)/2\). At fixed \(w\), the coordinate formulas give
\[
 \begin{aligned}
 \xi&=\tau(1+C_1\tau+O(\tau^2,w)),\\
 h&=\sqrt{a_0}H\bigl(1+(a_1/a_0+C_1)\xi/2+O(\xi^2,w)\bigr),\\
 \frac{\tau^2h\,d\tau}{1+\tau}
   &=\sqrt{a_0}\,\xi^2H\bigl(1+\kappa\xi+O(\xi^2,w)\bigr)d\xi.
 \end{aligned}
\]
In the last line the coefficient of \(\xi\) is
\((a_1/a_0+C_1)/2-4C_1-1=\kappa\).
On fixed elliptic cycles avoiding the two attachments the common
term has an expansion with even powers of \(s\), since \(\tau\) depends
only on \(\xi=s^2Z\) and \(w=s^4\). The odd terms at degrees seven
and nine are therefore exactly
\[
          \sqrt{a_0},Z^2 W\,dZ,
          \qquad \sqrt{a_0}\,\kappa Z^3 W\,dZ .
\]
These expansions are uniform on the chosen cycles. With
\(\eta=dZ/R\), differentiation gives
\[
 d(Z^jR)=\bigl((j+\tfrac12)Z^j-(j+\tfrac32)A Z^{j+2}\bigr)\eta.
\]
Taking \(j=0,2\) and \(j=1,3\), respectively, gives
\[
 [Z^2 W\,dZ]=\frac{2}{21A}[\eta],\qquad
 [Z^3 W\,dZ]=\frac{2}{15A}[Z\eta].
\]
The first basis differential is holomorphic and the second has
one double pole at infinity with zero residue. A meromorphic
primitive with a single simple pole is impossible on an elliptic
curve, so the classes are independent in compact de Rham
cohomology. The two odd period coefficients detect the entire
elliptic quotient of local first homology when \(\kappa\ne0\).

On either boundary circle, use a fixed small nonzero \(|H|\). The
solutions for \(\xi\), and hence for \(\tau\) and \(h\), are analytic
in \(w\) along that circle. The integral of the common term is
\(2\pi i\) times its analytic residue at \(\tau=0\) times the fixed
winding number in \(\tau\). The other term is analytic in \(w\)
directly on the circle. Thus boundary periods contribute only
multiples of four to the powers in \(s\). If there are ordinary
nodes at the same finite value, their periods are analytic in \(w\)
by the fixed-circle representations of the
[Lefschetz calculation](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_LEFSCHETZ_EXCLUSION.md).

## Detecting the local variation

If the moments vanish eventually, let \(\gamma\) be the transported
inverse-branch circle, with identically zero period by the
[inverse-period condition](WINDING_ONE_CUBIC_QUADRUPLE_ROOT_INVERSE_PERIOD.md).
Let \(T_*\) be the monodromy supported in the \(A_3\) neighborhood,
and \(T_{\rm old}\) the product of any ordinary twists over
\(\zeta_*\), supported away from it. Both supports can be chosen
disjoint and the actual monodromy over \(\zeta_*\) is their product.
Its variation on \(\gamma\) is the sum of the separate local
variations. A representative of \((T_*-1)\gamma\) in the local
homology is the sum of a lift from elliptic homology and a boundary
class. The period of the whole variation is zero. The two odd
coefficients calculated above force the local elliptic projection
to vanish, even in the presence of coinciding ordinary nodes.

The remaining local variation is then a multiple of the boundary
class. The map \(T_*\) fixes that class, hence
\((T_*^4-1)\gamma=4(T_*-1)\gamma\) in punctured homology. The
attachment twist calculation makes the left side zero, so the
entire \(A_3\) variation vanishes there. On the elliptic quotient
\(T_*^{-1}-1\) is invertible. For each local cycle \(\delta\),
symplecticity gives
\[
 0=\langle(T_*-1)\gamma,\delta\rangle
   =\langle\gamma,(T_*^{-1}-1)\delta\rangle .
\]
The local boundary pairs to zero with \(\gamma\), so \(\gamma\)
is orthogonal to all the local vanishing cycles. The remaining
variation at \(\zeta_*\) is that of the ordinary nodes. Its zero
period then gives orthogonality to those nodes as in the usual
single-node and three-node tests. All other finite values use the
same tests unchanged.

The resolved total space and the connected inner and outer parts
of the large fiber are the ones in the low-branch exclusion.
For spanning, the \(A_3\) chart can be perturbed to three local
Morse points by adding \(\varepsilon H\) to its completed-square
form, as in the
[double-root Lefschetz proof](WINDING_ONE_CUBIC_REPEATED_MINUS_ONE_LEFSCHETZ_EXCLUSION.md#compact-spanning-and-the-parameter-conditions).
Their cycles lie in the original local vanishing homology. The
finite handle argument then gives compact spanning by that
homology and the other finite vanishing cycles. Their joint
orthogonality to \(\gamma\) contradicts its nonseparating compact
class and rules out eventual purity.

## Higher even orders of contact

**Theorem EvenOrderLowBranchCollisionExclusion.** Keep the global
hypotheses of TangentLowBranchCollisionExclusion, and suppose
\(K_0(\tau)=K(\tau,\zeta_*)=b\tau^m+O(\tau^{m+1})\), with \(b\ne0\)
and even \(m\ge4\). Write \(K=K_0+w a(\tau,w)\) as above and define
\[
 \xi=\tau\left(\frac{K_0(\tau)}{b\tau^m}\frac{a_0}{a(\tau,w)}\right)^{1/m},
 \qquad H=h\sqrt{\frac{\tau}{a(\tau,w)\xi}}.
\]
Use the unit with constant term one in the definition of \(\xi\),
and use the same square root as in the order-two chart. In the
inverse coordinate \(\tau=\tau(\xi,w)\) put
\[
 J(\xi,w)=\left(\frac{\tau}{\xi}\right)^2
        \sqrt{\frac{a(\tau,w)\xi}{a_0\tau}}
        \frac{\partial\tau}{\partial\xi}\frac1{1+\tau},
 \qquad J(\xi,0)=\sum_{j\ge0}J_j\xi^j .
\]
Assume that \(J_1,\ldots,J_{m-1}\) are nonzero. Then infinitely
many positive pure sphere moments are nonzero, including when
\(\zeta_*\) coincides with zero or one interior critical value.

The same coordinate substitution gives \(w=\xi H^2+A\xi^m\) with
\(A=-b/a_0\). For \(w=s^{2m}\), use \(\xi=s^2Z\) and
\(H=s^{m-1}W\). The compact completion of the scaled curve is
\[
          Z(W^2+AZ^{m-1})=1,\qquad R=ZW,\qquad R^2=Z-AZ^{m+1}.
\]
The polynomial on the right is squarefree of odd degree \(m+1\),
so its hyperelliptic curve has genus \(m/2\) and a single point at
infinity. Near \(Z=0\), \(W\) has a simple pole. At infinity,
\(Z\) and \(W\) have poles of orders two and \(m-1\), respectively.
There are no other poles of \(Z\) or \(W\). The original local fiber
in a small rectangle in \(\xi,H\) is the scaled curve with two small
disks removed at those points. One may take the \(H\)-radius larger
than \(2\sqrt{|A|}\) times the \(\xi\)-radius to the power
\((m-1)/2\). For small \(s\), the \(H\)-cutoff is then the active
bound near \(Z=0\) and the \(\xi\)-cutoff the active bound at infinity.

At \(Z=0\) use \(t=1/W\) and write \(Z=t^2 j(t)\), where
\(j+A t^{2m}j^m=1\) with \(j(0)=1\). The gluing to the outer
coordinate \(H\) has equation \(Ht=s^{m-1}\). At infinity use
\(Z=t^{-2}\) and \(W=t^{-(m-1)}V(t)\), where
\(V(t)^2=t^{2m}-A\) and \(V(0)\ne0\). Setting \(r t=s\) gives
\(\xi=r^2\), \(H=r^{m-1}V(t)\). At \(s=0\), \(r\) parametrizes the
single strict branch \(H^2=-A\xi^{m-1}\). The two gluing equations
give \(m-1\) positive twists toward the exceptional curve and one
toward the strict branch. The rational exceptional curve has one
attachment and carries the marked section. Both circles are again
peripheral, by the boundary of the two-ended local smoothing. The
base-changed twists act trivially on punctured homology by the
intersection formula. Base monodromy fixes each boundary class.
On compact core homology its \(m\)-th power is the hyperelliptic
involution, as seen by replacing \(s\) by a primitive \(2m\)-th
root of unity times \(s\) in the scaling. Hence it has no fixed
nonzero class in that quotient.

The part of the period form that depends only on \(\tau,w\) yields
even powers in \(s\) on fixed compact core cycles. On the two outer
circles its period is an analytic residue times a fixed winding,
as before. The remaining part is exactly
\(\sqrt{a_0}\,\xi^2H J(\xi,w)d\xi\). Its first \(m\) odd coefficients,
at the degrees \(m+5+2j\) for \(0\le j<m\), are
\(\sqrt{a_0}J_j Z^{j+2}W\,dZ\). Contributions from the \(w\)-dependence
of \(J\) start at degree \(m+5+2m\). The boundary periods are analytic
in \(w\). On the exceptional side this follows by solving for
\(\xi\) on a fixed small \(H\)-circle, where the derivative in
\(\xi\) at \(w=0,\xi=0\) is \(H^2\). On the strict side use a fixed
small \(r\)-circle, where
\(H=r^{m-1}\sqrt{-A+w/r^{2m}}\). As before, the form of its common
term permits a residue even if the projected \(\tau\)-circle
shrinks. The ordinary nodes sharing this value also contribute
only analytic functions of \(w\).

Put \(\eta=dZ/R\). The derivative identity
\[
 d(Z^nR)=\bigl((n+\tfrac12)Z^n-(n+\tfrac{m+1}{2})A Z^{n+m}\bigr)\eta
\]
shows that
\([Z^n W\,dZ]=\frac{m}{2n+m+1}[Z^n\eta]\), and that
\([Z^{n+m}\eta]\) is a nonzero multiple of \([Z^n\eta]\) for
\(n\ge0\). Thus the first \(m\) odd coefficients give nonzero
multiples of the classes \([Z^r\eta]\) for \(0\le r<m\), up to
permutation, since \(J_0=1\). All these forms are regular at finite
points and have zero residue at infinity by the local coordinate
\(Z=t^{-2}\). To see independence, an exact combination would
have a meromorphic primitive with no finite poles. The forms
change sign under \(R\mapsto-R\), so the anti-invariant part is
itself a primitive. Its ratio with \(R\) is a rational function of
\(Z\), with
no finite pole, since the anti-invariant primitive vanishes at
every finite branch point. Hence that ratio is a polynomial \(B\).
If \(B\ne0\), the displayed identity
makes the highest degree in the derivative at least \(m\). The
combination has degree below \(m\), so \(B=0\) and the combination
is zero. These \(m\) independent classes form a basis of compact
de Rham cohomology for the genus \(m/2\) curve.

The first \(m\) odd coefficients apply the same period test to the
genus \(m/2\) compact quotient. The local variation of the inverse-branch
circle has zero compact projection, hence is peripheral. Since that
variation is fixed by local monodromy and the \(2m\)-th power acts
trivially in punctured homology, the variation is zero there.
Invertibility of local monodromy minus the identity on the compact
quotient gives orthogonality to all local cycles, and the remaining
ordinary variations are detected as before.

The local cycles here also generate the vanishing directions in
the compact handle argument. Add \(\varepsilon\xi\) to the local
function. Its two critical points with \(\xi=0\) have
\(H^2=-\varepsilon\) and value zero. Its other \(m-1\) critical
points have \(H=0\), \(mA\xi^{m-1}=-\varepsilon\) and the distinct
nonzero values \((1-1/m)\varepsilon\xi\). All \(m+1\) points are
Morse by the Hessian determinant. A small further term \(\delta H\)
separates the first two values because their derivatives in
\(\delta\) are the two distinct \(H\)-coordinates. Cutting off
both perturbations away from the isolated critical point leaves
the differential nonzero on the cutoff annulus. The new local
vanishing cycles lie in the original smoothing region, so the
same compact spanning contradiction applies.
