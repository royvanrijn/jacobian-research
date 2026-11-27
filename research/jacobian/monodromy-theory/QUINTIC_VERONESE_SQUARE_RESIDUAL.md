# The affine first component in the quintic Veronese square case

**Theorem QuinticVeroneseSquareAffineFirstComponent.** Suppose that a complex polynomial map \(F\) in three variables has degree five, nonzero constant Jacobian determinant and leading part
\[
H_5=(zx^2+g_0(x,y))(x^2,xy,y^2),
\]
where \(g_0\) is a homogeneous binary cubic. Suppose that the coefficient of \(t^6\) in \(F(tx,ty,t^2)\) is \(N(x^2,xy+\beta)\), where \(N(u,v)=(u^2,uv,v^2)\) and \(\beta\ne0\). Then \(x\mid g_0\) and the first component \(F_1\) is affine in \(z\).

Write \(\{\cdot,\cdot\}\) for the binary Jacobian in \(x,y\). With \(p=[z]F_1\) and \(b=[z^2]F_2\), one has \(\{p,b\}=0\), with leading binary parts \(p_4=x^4\) and \(b_2=\beta x^2\). In polynomial source coordinates \(s,v,z\), the map has the form
\[
 F_1=p(s)z+A(s,v),\quad
 F_2=b(s)z^2+q(s,v)z+B(s,v),\quad
 F_3=\delta z^3+\rho(s,v)z+C(s,v),\quad \delta=\beta^2 .
\]
In \(\mathbb C(s)[v]\), put \(f=p'/b'\), with primes denoting differentiation in \(s\). Then there exist \(d,e\in\mathbb C(s)\) with
\[
 A=fq+d,\qquad
 3\delta\bigl(p' B+\tfrac12 f' q^2+d' q\bigr)+(p b'-2b p')\rho=e.
\]
The hypotheses describe the non-scalar branch of the
[linear-\(z\) weighted reduction](QUINTIC_VERONESE_LINEAR_Z_WEIGHTED_REDUCTION.md).
They apply after extending the finitely generated coefficient field of a characteristic-zero example to \(\mathbb C\).

## Two weighted determinant coefficients

Set \(u=x^2,\ w=xy+\beta,\ \nu=(w^2,-2uw,u^2)\). Write
\[
 F(tx,ty,t^2)=t^6N(u,w)+t^5 M_1+t^4 M_2+\cdots .
\]
The determinant in \(x,y,t\) is \(2\kappa t^3\), where \(\kappa=\det JF\). Invert \(x\) and work in \(R=\mathbb C[x,x^{-1},w]\). The variables \(u,w\) have Jacobian \(2x^2\) in \(x,y\), so the determinant in \(u,w,t\) has order \(t^3\). The derivations in \(u,w\) extend to \(R\). In particular \(E=u\partial_u+w\partial_w=\tfrac12x\partial_x+w\partial_w\), where the derivative in \(x\) holds \(w\) fixed.

If \(t^{6-r}M\) is the first correction to \(t^6N(u,w)\), the determinant coefficient at \(t^{17-r}\) is
\[
 (12-2r)\nu\cdot M-6\nu\cdot E M
   =(24-2r)S-6E S,\qquad S=\nu\cdot M.
\]
This is the
[mixed determinant identity in the weighted reduction](QUINTIC_VERONESE_LINEAR_Z_WEIGHTED_REDUCTION.md#excluding-a-simple-common-zero-of-the-quadratic-arguments),
using \(E\nu=2\nu\). The eigenvalues of \(E\) on \(R\) are half-integers because \(x^i w^j\) has eigenvalue \(i/2+j\). At \(r=1,2\) the equation \((24-2r)S-6ES=0\) therefore forces \(S=0\).

For \(r=1\) we obtain \(\nu\cdot M_1=0\). This equality holds in \(\mathbb C[x,y]\). The polynomials \(u,w\) are coprime since \(\beta\ne0\). Reduction of the equation modulo \(u\) and \(w\) shows that the first entry of \(M_1\) is divisible by \(u\) and its third entry by \(w\). Substitution then gives
\[
 M_1=DN(u,w)(U,V),\qquad U,V\in\mathbb C[x,y].
\]
The ordinary degree-eleven determinant equation in the weighted reduction expresses \(H_4\) as a polynomial combination of the columns of \(DN(x,y)\). Write
\[
 M_1=g_0N(x,y)+DN(x,y)(P_2,Q_2)+L_1
\]
where \(P_2,Q_2\) are homogeneous of degree two and \(L_1=[z^2]H_3\) is homogeneous of degree one. This is the decomposition by ordinary degree of the coefficient at \(t^5\). In its first coordinate the equality with \(DN(u,w)(U,V)\) gives
\[
 2x^2 U=x^2g_0+2xP_2+(L_1)_1.
\]
Divisibility by \(x^2\) and the distinct binary degrees give \((L_1)_1=0,\ x\mid P_2\) and \(U=g_0/2+P_2/x\). At \(x=0\) the second coordinate gives \(\beta U(0,y)=(L_1)_2(0,y)\). The latter is homogeneous linear, so the binary cubic \(g_0\) vanishes on \(x=0\).

Use the \(t^{-1}\)-adic translation in \(x,y\) by \(-t^{-1}W\), where \(W=D(u,w)^{-1}(U,V)\) has entries in \(R\). It removes \(M_1\) and multiplies the determinant by a unit \(1+O(t^{-1})\). The new second correction is
\[
 M_2'= M_2-dM_1[W]+\tfrac12d^2(N(u,w))[W,W].
\]
By the chain rule and \(D(u,w)W=(U,V)\), the last two terms together are congruent to \(-N(U,V)\) modulo the image of \(DN(u,w)\). Indeed \(dM_1[W]\) contributes \(2N(U,V)\) modulo this image, and the second derivative term contributes \(N(U,V)\). The order-\(t^{15}\) coefficient with \(r=2\) now gives \(\nu\cdot(M_2-N(U,V))=0\), again a polynomial equality in \(x,y\). At \(x=0\) it gives \((M_2)_1(0,y)=U(0,y)^2\), a multiple of \(y^2\).

The first coordinate of \(M_2\) is the sum of the first coordinates of \([z^0]H_4,\ [z]H_3,\ [z^2]H_2\). At \(x=0\) the first summand vanishes by tangency to \(N(x,y)\), the second is a multiple of \(y^2\) and the last is constant. Thus \(([z^2]H_2)_1=0\). The coefficient \(([z^2]H_4)_1\) vanishes by the weighted leading form and \(([z^2]H_3)_1=(L_1)_1=0\) was proved above. The first coordinate has zero cubic coefficient in \(z\) by the weighted leading form and \([z^3]H_4=0\) by the preceding reduction. This proves the affine assertion.

## Equations in polynomial source coordinates

Initially \(F_1=pz+P,\ F_2=bz^2+qz+Q,\ F_3=\delta z^3+cz^2+rz+T\), with all coefficients binary polynomials except for constant \(\delta\). The \(z^5\) coefficient of the determinant is \(3\delta\{p,b\}\). Here \(b=\beta x^2+\ell(x,y)+b_0\), with \(\ell\) linear. If its \(y\)-coefficient is nonzero, use \((s,v)=(b,x)\) as polynomial coordinates. If that coefficient is zero, use \((s,v)=(x,y)\). Each change has nonzero constant determinant, and in each choice \(\{p,b\}=0\) now says \(p\) depends only on \(s\). The substitution \(z\mapsto z-c/(3\delta)\) depresses the third component and preserves the leading coefficients \(p,b,\delta\). This gives the asserted form. Both \(p'\) and \(b'\) are nonzero.

The \(z^4\) and \(z^3\) coefficients of the determinant of that form are
\[
 \begin{aligned}
 &3\delta(p' q_v-b'A_v),\\
 &3\delta(p' B_v+A_s q_v-A_v q_s)+(p b'-2b p')\rho_v.
 \end{aligned}
\]
Both vanish. The first gives \(A=fq+d\). Substitution into the second gives the \(v\)-derivative of the asserted expression, so that expression belongs to \(\mathbb C(s)\).

The [symbolic checker](verify_quintic_veronese_square_residual.py) replays the Taylor identities for the second correction and the determinant coefficients. The divisibility and weighted arguments are in the written proof. The [full square exclusion](QUINTIC_VERONESE_SQUARE_EXCLUSION.md) uses these equations.
