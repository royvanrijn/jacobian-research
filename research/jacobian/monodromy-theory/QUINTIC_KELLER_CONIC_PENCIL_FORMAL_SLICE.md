# A formal slice obstruction for the quintic square pencil

**Theorem QuinticKellerConicPencilFormalSlice.** Let \(\ell\) be nonzero linear, and let \(q\) be a homogeneous quadratic form in three variables over \(\mathbb C\). If \(q\) has rank two on the plane \(\ell=0\), there is no degree-five Keller map with homogeneous leading part
\[
                         (\ell^5,\ \ell^3 q,\ \ell q^2).
\]
The same conclusion holds over characteristic-zero fields after extending the finitely generated coefficient field to \(\mathbb C\).

This treats the rank-two restriction in the conic-pencil case with \(\ell^2\) in the pencil, beyond the scope of the
[first-integral exclusion](QUINTIC_KELLER_CONIC_PENCIL_FIRST_INTEGRALS.md).
Wang and Yan use the same rank distinction to normalize the square pencil in the quartic case, [printed pages 12–13, equation (3.12)](https://arxiv.org/pdf/2610.36495v1).

Take \(\ell=z\). A binary linear change, followed by adding multiples of \(z\) to \(x,y\), makes \(q=xy+\gamma z^2\). Replacing \(q\) by \(q-\gamma z^2\) corresponds to a constant linear target change on the three Veronese components. Thus assume \(H_1=z^5,\ H_2=z^3xy,\ H_3=zx^2y^2\).

Let \(F=(F_1,F_2,F_3)\) have these homogeneous leading parts and determinant \(\kappa\ne0\). There is a formal substitution
\[
                 z=Z(x,y,t)=t+c+O(t^{-1})
\]
over \(\mathbb C[x,y]((t^{-1}))\) with \(F_1(x,y,Z)=t^5\) and \(c\) constant. Indeed set \(Z=t U\) and solve \(t^{-5}F_1(x,y,tU)=1\) for \(U=1+O(t^{-1})\). The derivative of the leading equation \(U^5=1\) at \(U=1\) is \(5\), so the equation determines each successive coefficient by division by \(5\). The coefficient of \(t^{-1}\) in \(U\) is \(-1/5\) times the coefficient of \(z^4\) in \(F_1\), which is constant because \(H_1=z^5\).

The coefficient of \(z^4\) in \(F_2\) is also constant, say \(b\), and its coefficient of \(z^3\) is \(xy\) plus an affine binary polynomial. Consequently
\[
              t^{-3}\bigl(F_2(x,y,Z)-b t^4\bigr)
                =xy+a x+d y+e+O(t^{-1})
\]
for constants \(a,d,e\). At \(t^{-1}=0\) its gradient in \(x,y\) vanishes at \((-d,-a)\), with invertible Hessian. The two gradient equations therefore have a solution \(x(t),y(t)\) in \(\mathbb C[[t^{-1}]]\) with those constant terms. This follows recursively because at each order the unknown pair of coefficients is multiplied by the same invertible Hessian.

After the substitution, the determinant in \(x,y,t\) is \(\kappa Z_t\). Since the first component is \(t^5\), the binary Jacobian of \(F_2(x,y,Z), F_3(x,y,Z)\) at fixed \(t\) is \(\kappa Z_t/(5t^4)\). At \(x(t),y(t)\) the left side vanishes by the two gradient equations. The right side is nonzero there, with leading term \(\kappa/(5t^4)\) because \(Z_t=1+O(t^{-2})\). This contradiction proves the theorem.
