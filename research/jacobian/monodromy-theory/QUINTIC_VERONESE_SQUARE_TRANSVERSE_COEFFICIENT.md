# The quadratic coefficient of the second component in the quintic Veronese square case

**Theorem QuinticVeroneseSquareBinaryCoefficient.** Under the hypotheses of
[QuinticVeroneseSquareAffineFirstComponent](QUINTIC_VERONESE_SQUARE_RESIDUAL.md),
the coefficient \(b=[z^2]F_2\) depends only on \(x\). Thus in the weighted square form \(N(x^2,xy+\beta)\), the coefficient \(b\) is a quadratic in \(x\) with leading term \(\beta x^2\). The [full square exclusion](QUINTIC_VERONESE_SQUARE_EXCLUSION.md) uses this restriction.

## Roots of the first component's linear coefficient

Scale \(y\) by \(\beta\) and act on the target by the inverse of \(\operatorname{diag}(1,\beta,\beta^2)\). This gives \(\beta=1\) and preserves the nonzero constant determinant condition. If \(b\) had a nonzero \(y\)-coefficient \(\gamma\), take \(s=b=x^2+a x+\gamma y+b_0\) and \(v=x\) as source coordinates. By the preceding theorem \(p=[z]F_1\) is a polynomial in \(s\). Its binary degree is four and its leading term is \(x^4\), so
\[
                         p=s^2+\eta s+\theta .
\]
Use the depressed form of that theorem, with \(b=s,\ \delta=1\) and \(A=p'q+d(s)\). Here \(d\) belongs to \(\mathbb C[s]\), as \(A,p',q\) are polynomials.

At a simple root \(s_0\) of \(p\), the polynomial \(q(s_0,v)\) has degree one in \(v\). In fact, if \(q_v(s_0,v_0)=0\), the first component at \(s=s_0,v=v_0\) has zero \(v\) and \(z\) derivatives. Its \(s\) derivative is linear in \(z\) with nonzero coefficient \(p'(s_0)\). Choosing \(z\) to make it vanish contradicts the Keller condition. Hence \(q_v(s_0,v)\) has no zero and is a nonzero constant.

At a double root \(s_0\) of \(p\), the first row of the derivative matrix on that plane is \((2q+d',0,0)\). The constant nonzero determinant makes \(2q(s_0,v)+d'(s_0)\) a unit in \(\mathbb C[v,z]\), so \(q(s_0,v)\) is constant.

## The cubic coefficient on the planes \(b=s\)

Work first in the original binary variables \(x,y\). The ordinary degree-four term \(H_4\) is tangent to \(N(x,y)\) by the
[degree-eleven determinant identity](QUINTIC_VERONESE_LINEAR_Z_WEIGHTED_REDUCTION.md#from-ordinary-to-weighted-degree).
Write \([z]H_4=DN(x,y)(P_2,Q_2)\) with \(P_2,Q_2\) homogeneous quadratics. The cubic term of \(p\) is \(2xP_2\). From \(p=s^2+\eta s+\theta\) it follows that \(P_2=x(ax+\gamma y)\). Write \(Q_2=Lx^2+Mxy+N_2y^2\), and let \(D,E\) be the \(y^2,xy\) coefficients in the binary quadratic term of the original \(z\)-coefficient of \(F_2\). This original coefficient is
\[
 q_{\mathrm{old}}=x^3y+xQ_2+yP_2+D y^2+E xy+\text{terms in }x^2,\ x,\ y,\ 1.
\]
The original \(z^2\)-coefficient \(c\) of \(F_3\) is \(2xy\) plus an affine polynomial. Depressing the third component gives \(q=q_{\mathrm{old}}-2s c/3\).

Substitute \(x=v,\ y=(s-v^2-av-b_0)/\gamma\). The degree-five \(v\)-coefficient of \(q\) is \(N_2/\gamma^2\), which vanishes because \(q(s_0,v)\) has degree at most one at every root of \(p\). With \(N_2=0\), combine
\[
 x^3y+xQ_2+yP_2=xy(s-b_0)+Lx^3+M x^2 y .
\]
The degree-four \(v\)-coefficient of \(q\) is then \(-M/\gamma+D/\gamma^2\), also zero. With \(D=\gamma M\), the terms \(M x^2y+D y^2\) become \(M y(s-a x-b_0)\). Consequently
\[
 q=xy(-s/3-b_0-a M+E)+Lx^3+\text{terms of \(v\)-degree at most two}.
\]
In particular \([v^3]q\) is affine in \(s\) with slope \(1/(3\gamma)\). It must vanish at every root of \(p\), so \(p\) has a double root \(s_0\).

## The double root

On the plane \(s=s_0\), the first row is \((2q+d',0,0)\), with \(2q+d'\) a nonzero constant. The planar determinant in \(v,z\) of the last two components is therefore a nonzero constant. Write \(b=s_0\), \(q=q(s_0)\), and view \(B,\rho,C\) as polynomials in \(v\) on this plane. The coefficients at \(z^2,z^1,z^0\) of that determinant are
\[
 3B_v-2b\rho_v,\qquad -2b C_v-q\rho_v,\qquad B_v\rho-q C_v.
\]
If \(b\ne0\), the two vanishing coefficients make the constant coefficient equal to
\[
 \rho_v(2b\rho/3+q^2/(2b)).
\]
This cannot be a nonzero constant. Such a product would force both factors to be units, whence \(\rho\) would be constant and \(\rho_v=0\).

Thus \(s_0=b=0\) and \(p=s^2\). The same planar equations first give \(B_v(0,v)=0\), then \(q(0)\ne0\) and \(\rho_v(0,v)=0\). The integrated \(z^3\) equation of the preceding theorem becomes
\[
                   2s B+q^2+d'q-s^2\rho=e_0(s)
\]
with \(e_0\in\mathbb C[s]\), since the left side is polynomial. Differentiate in \(s\) and set \(s=0\). Because \(B(0,v)\) and \(q(0,v)\) are constant in \(v\), this makes \((2q(0)+d'(0))q_s(0,v)\) constant in \(v\). Its first factor is nonzero, while \([v^3]q_s(0,v)=1/(3\gamma)\). This contradiction proves the theorem.

The [symbolic checker](verify_quintic_veronese_square_transverse_coefficient.py) replays the coefficients on the planes \(b=s\) and the planar determinant calculation. The root argument and the use of the determinant equations are in the written proof.
