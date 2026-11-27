# Degree obstructions for separate coordinate changes of the Alpöge–Gallagher map

**Theorem AlpogeGallagherSeparateCoordinateDegreeObstructions.**
Let \(F=(A,B,C)\) be the Alpöge–Gallagher degree-seven Keller map
over a characteristic-zero field. Any polynomial automorphism
\(\Theta\) of the target satisfies \(\deg(\Theta\circ F)\ge7\).
Any dominant polynomial map \(\Psi\) of the source satisfies
\(\deg(F\circ\Psi)\ge7\). The statements allow an additional
affine automorphism on the opposite side.

Use the formulas in [Gao, Theorem 3.3 and its displayed expansions](sources/gao-2608.00222v2.pdf).
Set \(u=1+xy\), \(g=2-3xy-x^2z\), \(w=gu\), and write
\[
 C=xg,\qquad
 BC=-3w^2+4w+2g,\qquad
 AC^2=-w^3+w^2+gw.
\]
The identity \(x^2 A=u^2(u+1-gu)\) follows.

## Polynomial changes of the target

Give \(x,y,z\) weights \(1,-1,-2\) respectively. Then \(A,B,C\)
have weights \(-2,-1,1\), so substitution preserves the weight
decomposition for this target grading. Each weight part of a
polynomial of ordinary degree at most six again has ordinary
degree at most six.

Use \(X=BC,\ Y=AC^2\). Any target polynomial \(V\) of weight
\(-2\) satisfies \(C^2V=J(X,Y)\) for some polynomial
\(J\in(Y,X^2)\). Indeed, in a monomial \(A^iB^jC^k\) of that
weight, \(k+2=2i+j\). The resulting monomial after multiplication
by \(C^2\) is \(Y^iX^j\), with \(i\ge1\) or \(j\ge2\).
A nonzero \(V\) gives a nonzero \(J\). Here and below, the equations
in \(A,B,C\) are polynomial identities before evaluation at \(F\).

Let \(\Delta(X,Y)\) be the discriminant of
\(T^3-2T^2+XT-2Y\). It is
\[
 \Delta=4X^2-4X^3-64Y-108Y^2+72XY.
\]
After pullback, \(T=w\) is a root, the derivative at that root is
\(2g\), and the quotient quadratic is
\(T^2+(w-2)T+(w^2-2w+X)\). Its discriminant gives
\[
 \Delta(B C,A C^2)=4g^2(9w^2-12w+4-8g).
\]

In the source's ordinary degree, the leading forms of \(g,w\)
are \(-x^2z,-x^3yz\), of degrees three and five. The exponent
vectors of these two monomials are independent over \(\mathbb Q\).
Consequently distinct monomials of a polynomial in \(g,w\) have
distinct leading monomials, and its source degree is its weighted
degree in \(g,w\) with weights three and five. In those weights
\(X,Y,\Delta\) after pullback have leading terms proportional to
\(w^2,w^3,g^2w^2\), respectively.

Since \(\Delta\) has degree two in \(Y\) with constant nonzero
leading coefficient, \(\mathbb C[X,Y]\), or the same ring over
the coefficient field, is free over the polynomial ring
in \(X,\Delta\), with basis \(1,Y\). This follows by division
by \(\Delta(X,Y)-D\) in \(Y\), treating \(D\) as an independent
variable. The pullbacks of \(X^iY^j\Delta^k\), \(j=0,1\), have
distinct maximal monomials \(w^{2i+3j+2k}g^{2k}\): the \(g\)-exponent
fixes \(k\), and then parity of the \(w\)-exponent fixes \(j\).
Thus the ordinary degree of the pullback of \(J\) is the maximum
of \(10i+15j+16k\) over its nonzero coefficients in this basis.

If \(V\circ F\) had degree at most six, the degree of \(J\circ F\)
would be at most \(6+2\deg C=14\), where pullback of \(J\) uses
\(X=BC,Y=AC^2\). The bound leaves only \(1,X\) in the expansion
of \(J\), which is impossible for nonzero \(J\in(Y,X^2)\).
In particular every target polynomial with pullback of degree at most
six has zero weight-\(-2\) part and hence zero coefficient of \(A\)
in its linear part. A polynomial automorphism has invertible
linear part at the target origin, which rules out three such outputs.

## Polynomial changes of the source

Let \((x_1,y_1,z_1)\) be a dominant polynomial substitution,
with strictly positive degrees \(d_x,d_y,d_z\).
Put \(u_1=1+x_1y_1\) and \(g_1=2-x_1(3y_1+x_1z_1)\).
Dominance rules out \(3y_1+x_1z_1\) being constant, so
\(d_g=\deg g_1\ge d_x+1\).
The term \(g_1u_1\) strictly dominates \(u_1+1\) by ordinary
degree. The identity for \(x^2 A\) therefore gives
\[
                     \deg(A\circ\Psi)=d_x+3d_y+d_g .
\]
If \(d_y\ge2\), this degree is at least nine. If \(d_y=1\),
then \(d_x+d_z>1\), so \(x_1z_1\) strictly dominates
\(3y_1\) and \(d_g=2d_x+d_z\). The displayed degree is then
\(3d_x+3+d_z\ge7\). An affine automorphism on the other side
preserves the ordinary degree of the vector map.

The proof addresses separate polynomial changes and leaves simultaneous
nonlinear changes of both source and target open.
