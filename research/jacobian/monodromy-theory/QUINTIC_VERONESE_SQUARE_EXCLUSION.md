# Excluding the weighted square case of the quintic Veronese class

**Theorem QuinticVeroneseWeightedSquareExclusion.** There is no complex polynomial map \(F\) of degree five with nonzero constant Jacobian determinant, ordinary leading part
\[
 H_5=(zx^2+g_0(x,y))N(x,y),\qquad N(x,y)=(x^2,xy,y^2),
\]
and coefficient \(N(x^2,xy+\beta)\) at \(t^6\) in \(F(tx,ty,t^2)\), where \(g_0\) is homogeneous cubic and \(\beta\ne0\). The same conclusion holds in characteristic zero by embedding the finitely generated coefficient field in \(\mathbb C\).

Together with the
[weighted classification](QUINTIC_VERONESE_LINEAR_Z_WEIGHTED_REDUCTION.md)
and the [scalar exclusion](QUINTIC_VERONESE_SCALAR_WEIGHTED_EXCLUSION.md),
this excludes the entire case of \(z\)-degree one for the homogeneous cubic factor of the quintic Veronese leading part.

## The coefficients on a plane \(x=a\)

By the
[affine first component theorem](QUINTIC_VERONESE_SQUARE_RESIDUAL.md)
and the [quadratic coefficient restriction](QUINTIC_VERONESE_SQUARE_TRANSVERSE_COEFFICIENT.md),
\(F_1\) is affine in \(z\), \(p=[z]F_1\) depends only on \(x\), and \(b=[z^2]F_2\) is quadratic in \(x\). Scale \(y\) by \(\beta\) and act on the target by the inverse of \(\operatorname{diag}(1,\beta,\beta^2)\) to make \(\beta=1\). Translate \(x\) to make \(b=x^2+b_0\). This translation preserves the ordinary and weighted leading parts. Here and below the nonzero constant Jacobian determinant may be changed by a nonzero scalar in the coordinate choices.

The ordinary degree-four term is tangent to \(N(x,y)\). Write
\[
 [z]H_4=DN(x,y)(P_2,Q_2),\qquad L_1=[z^2]H_3 .
\]
The proof of the affine first component theorem gives
\(g_0N(x,y)+DN(x,y)(P_2,Q_2)+L_1=DN(x^2,xy+1)(U,V)\), where
\(U=g_0/2+P_2/x\) and \(U,V\) are polynomials. The linear part in the second coordinate of the right side is \(P_2/x\), since the second coordinate is \(x^2V+(xy+1)U\). Thus \(P_2/x=(L_1)_2\), the linear part of \(b\). It follows that \(P_2=0\). In particular \(p\) is monic quartic with zero cubic coefficient.

Depress the third component by replacing \(z\) by \(z-c/3\), where \(c\) is its original quadratic coefficient in \(z\). Use the notation
\[
 F_1=pz+A,\qquad F_2=bz^2+qz+B,\qquad F_3=z^3+\rho z+C .
\]
The coefficients other than \(p,b\) are binary polynomials and \(c=2xy+\text{affine terms}\). The two equations from the affine first component theorem give, in \(\mathbb C(x)[y]\),
\[
 A=fq+d,\qquad
 3\bigl(p'B+\tfrac12 f'q^2+d'q\bigr)+\theta\rho=e,
 \qquad f=p'/b',\quad \theta=p b'-2b p',
 \tag{1}
\]
with \(d,e\in\mathbb C(x)\) and primes denoting \(x\) derivatives. The polynomial \(q\) is obtained by subtracting \(2bc/3\) from the original linear coefficient \(x^3y+xQ_2+\text{terms of binary degree at most two}\). Thus \(q\) has degree at most two in \(y\), the coefficient of \(y^2\) is affine in \(x\), and the coefficient of \(y\) has degree three with leading term \(-x^3/3\).

At a simple nonzero root \(a\) of \(p\), the coefficient \(f(a)\) is nonzero. If \(q_y(a,y)=0\) at some \(y\), the \(y,z\) derivatives of \(F_1\) vanish there, and its \(x\) derivative can be made zero by choosing \(z\), since \(p'(a)\ne0\). Hence \(q(a,y)\) is affine and nonconstant.

At a repeated nonzero root \(a\) of \(p\), \(f\) vanishes at \(a\) and the first row of the derivative matrix on \(x=a\) is \((f'q+d',0,0)\). This also holds at \(a=0\) if \(p=x^4\), since then \(f=2x^2\). All terms of (1) are regular at these points. The determinant condition makes \(f'q+d'\) a nonzero constant on the plane and the planar Jacobian of \(F_2,F_3\) a nonzero constant. The \(z^3\) coefficient of that planar Jacobian is \(3q_y\), so \(q(a,y)\) is constant. The next coefficients are
\[
 3B_y-2b\rho_y,\qquad -2b C_y-q\rho_y,\qquad B_y\rho-qC_y .
\]
If \(b(a)\ne0\), the last coefficient becomes
\(\rho_y(2b\rho/3+q^2/(2b))\), which cannot be a nonzero constant. In fact its two factors would be units, making \(\rho\) constant and \(\rho_y=0\). Therefore \(b(a)=0\), and the same equations show that \(q(a)\ne0\) and \(B_y(a,y)=\rho_y(a,y)=0\).

At any of these repeated roots, \(p(a)=p'(a)=b(a)=0\) gives \(\theta(a)=\theta'(a)=0\). Differentiate (1) in \(x\) and evaluate at \(a\). All terms except \(3(f'q+d')q_x\) are constant in \(y\) or vanish there, because \(q(a,y), B(a,y)\) are constant. Since \(f'q+d'\) is a nonzero constant on the plane, \(q_x(a,y)\) is constant in \(y\).

These facts force the coefficient of \(y^2\) in \(q\) to vanish identically. A repeated nonzero root, or \(p=x^4\) at zero, supplies a double zero of this affine coefficient. In the remaining case there are at least two simple nonzero roots of \(p\), because the roots sum to zero. Each gives a zero of the coefficient. Thus write
\[
 q=q_0(x)+k(x)y,\qquad k=-x^3/3+O(x^2).
\]
For any of the repeated roots just treated, taking coefficients of \(y\) in (1) gives
\[
 3p'[y]B+3(f'q_0+d')k+\theta[y]\rho=0 .
 \tag{2}
\]
At \(x=a\), the first and last coefficient polynomials \([y]B,[y]\rho\) vanish. Also \(f'q_0+d'\) is nonzero there, and \(\theta\) vanishes to at least the order of \(p'\), since \(\theta=p b'-2b p'\). It follows that \(k\) vanishes at \(a\) to at least the multiplicity of \(a\) in \(p\). In particular \(p=x^4\) is impossible, and \(p\) cannot have two distinct repeated nonzero roots.

## The determinant on a first-component fibre

The vanishing of the \(y^2\) coefficient of \(q\) implies that \(Q_2\) has no \(y^2\) term, since that coefficient was \(x[y^2]Q_2\) plus a constant. The original \(z\) coefficient of \(F_3\) was
\(x^2y^2+2yQ_2+\text{terms of binary degree at most two}\). Depressing subtracts \(c^2/3\) from it. Thus \(\rho\) has degree at most two in \(y\), and \(m=[y^2]\rho=-x^2/3+O(x)\). Equation (1) then gives \(\deg_y B\le2\) and
\[
 B_2=[y^2]B=-\frac{f'k^2}{2p'}-\frac{\theta m}{3p'} .
 \tag{3}
\]
The original constant coefficient of \(F_3\) as a polynomial in \(z\) has \(y\)-degree at most four since \(x\mid g_0\), so \(g_0y^2\) has \(y\)-degree at most four. The other contributions to the depressed constant coefficient are \(2c^3/27\) and \(-c/3\) times the original linear coefficient in \(z\), each of \(y\)-degree at most three. Thus \(\deg_y C\le4\).

Over \(\mathbb C(x)\), use \(u=F_1\) and \(y\) as coordinates on the remaining two variables. Write \(h=fk,\ d_0=fq_0+d\), and set
\[
 z=(u-d_0)/p-\alpha y,\qquad \alpha=h/p .
\]
Let \(G,H\) be \(F_2,F_3\) after this substitution. If \(\kappa\) is the nonzero constant Jacobian before substitution, then the binary Jacobian \(\{G,H\}_{x,y}\) at fixed \(u\) equals \(\kappa/p\). Write
\[
 G=L y^2+(M_1 u+M_0)y+J(x,u),\qquad
 L=b\alpha^2-k\alpha+B_2,\qquad M_1=(k-2b\alpha)/p .
\]
Here \([u^2]J=b/p^2\). The coefficient of \(u\) in \([y^2]H\) is
\[
 T_2=(3\alpha^2+m)/p .
\]
The coefficient of \(y^4\) in \(H\) is a polynomial \(C_4(x)=[y^4]C\), and that of \(y^3\) is independent of \(u\). At \(x=\infty\), the leading terms \(b\sim x^2,\ p\sim x^4,\ f\sim2x^2,\ k\sim-x^3/3,\ m\sim-x^2/3\) give
\[
 \alpha\sim-2x/3,\qquad M_1\sim x^{-1},\qquad T_2\sim x^{-2},
 \qquad L=O(x^3).
 \tag{4}
\]
For the last bound, (3) has leading term \((-1/18-1/6)x^4=-2x^4/9\), which cancels the \(2x^4/9\) from \(b\alpha^2-k\alpha\). Orders at infinity here refer to the Laurent series in \(x^{-1}\).

If \(C_4\ne0\), the \(u y^4\) coefficient in \(\{G,H\}_{x,y}\) implies \(4M_1' C_4-M_1 C_4'=0\). Thus \(C_4/M_1^4\) is a nonzero constant. This contradicts the polynomiality of \(C_4\) and \(M_1\sim x^{-1}\). Hence \(C_4=0\), and \(H\) has joint degree at most three in \(u,y\).

## When \(G\) is quadratic in \(y\)

Suppose \(L\ne0\). Take the homogeneous parts of \(G,H\) in \(u,y\) of degrees two and three, then set \(u=1\), to obtain \(Q,T\in\mathbb C(x)[y]\). The degree-four term in the binary Jacobian gives \(\{Q,T\}_{x,y}=0\). In particular
\[
 Q=L y^2+M_1 y+b/p^2,\qquad T(x,0)=1/p^3 .
\]
The second evaluation follows from the depressed monic cubic in \(z\), since its cubic part at \(y=0\) is \((u/p)^3\).

Use \(v=Q\), \(w=2Ly+M_1\), and \(D=M_1^2-4Lb/p^2\). Then \(w^2=4Lv+D\). This is a quadratic field extension of \(\mathbb C(x,v)\), because the right side is linear in \(v\) with nonzero leading coefficient. Reducing \(T\) in this extension gives \(T=R(v)+S(v)w\), with \(R,S\in\mathbb C(x)[v]\) of degrees at most one. The bracket condition says \(T_x=0\) at fixed \(v\). The two basis elements \(1,w\) therefore give \(R\in\mathbb C[v]\) and \(S^2(4Lv+D)\in\mathbb C[v]\).

If \(S=0\), then \(T=c_1 Q+c_0\) with constants \(c_1,c_0\). Evaluating at \(y=0\) gives \(1=c_1 bp+c_0p^3\), contrary to the existence of a root of \(p\). So \(S\ne0\). The unique root of odd multiplicity of \(S^2(4Lv+D)\) over an algebraic closure of \(\mathbb C(x)\) is \(-D/(4L)\). This is a complex constant \(\sigma\), since the polynomial belongs to \(\mathbb C[v]\). Its leading coefficient is \(4L\) times a rational square, so \(L\) itself is a square in \(\mathbb C(x)\). Consequently
\[
 b-\sigma p^2=\frac{p^2 M_1^2}{4L}
\]
is a polynomial square, since a polynomial that is a rational square is a polynomial square. If \(\sigma\ne0\) and \(R_0^2=b-\sigma p^2\), the product \((R_0-\sqrt{-\sigma}p)(R_0+\sqrt{-\sigma}p)=b\) forces both factors to have degree at most two, while their difference has degree four. If \(\sigma=0\), then \(L=p^2 M_1^2/(4b)\sim x^4/4\), contrary to (4). This excludes \(L\ne0\).

## When \(G\) is linear in \(y\)

Let \(L=0\), so \(G=M y+J\) with \(M=M_1u+M_0\) and \(M_1\ne0\). In the coordinates \(u,v=G,x\), solving \(y=(v-J)/M\) gives \(H\) rational in \(x,u,v\). Its derivative in \(x\) at fixed \(u,v\) is \(-\kappa/(pM)\). A rational derivative has zero residue at every finite pole. If the rational function \(-M_0/M_1\) in \(x\) were nonconstant, a generic complex \(u\) would give a simple zero of \(M\) at a finite point away from the zeros of \(p\) and from the zeros and poles of \(M_1\). Indeed the nonconstant rational function has only finitely many critical values and values at these exceptional points. Such a zero would give a nonzero residue, which is impossible. Thus \(M=M_1(u-u_0)\) for a complex constant \(u_0\).

At \(u=u_0\), \(G=J\) is independent of \(y\). The determinant equation becomes \(J_x H_y=\kappa/p\), which makes \(H\) affine in \(y\) there. As its \(y^3\) coefficient was independent of \(u\), it is identically zero. Its \(y^2\) coefficient is then \(T_2(u-u_0)\). For variable \(u\), the coefficient of \(v^2\) in \(H(x,u,(v-J)/M)\) is \(T_2/(M_1^2(u-u_0))\). Its \(x\) derivative is zero because the derivative of \(H\) is independent of \(v\). The asymptotics (4) imply \(T_2=M_1^2\).

Put \(\lambda=b f/p\). Substitution in \(T_2=M_1^2\) gives
\[
 m=k^2(1-2\lambda)^2/p-3f^2k^2/p^2 .
\]
Use this in (3) and in \(0=L=b\alpha^2-k\alpha+B_2\), and divide by \(k^2\). The resulting identity is
\[
 (1-2\lambda)^3+3\lambda^3 p/b^2+\frac{3f'}{2b'}=0 .
 \tag{5}
\]
In detail, \(\theta/p'=(p/f)(1-2\lambda)\), so after substitution \(L/k^2=-b f^2/p^2-f'/(2p')-(1-2\lambda)^3/(3f)\), which yields (5).

Every root of \(p\) must now be a root of \(b\). At a nonzero root \(a\) with \(b(a)\ne0\), \(\lambda=b p'/(p b')\) has a simple pole. The first term of (5) has a pole of order three, while the second has lower pole order because \(p(a)=0\) and the third is regular. If zero is a root of \(p\) of multiplicity \(j\) and \(b(0)\ne0\), then \(\lambda\) has a pole of order two, so the first term has a pole of order six. The second has lower pole order since \(j>0\), and \(f'/b'\) has pole order at most three because \(f=p'/(2x)\) and \(b'=2x\). This is again impossible.

If \(b=x^2\), the inclusion of roots gives \(p=x^4\). If \(b=x^2+b_0\) with \(b_0\ne0\), its roots are nonzero opposites. The zero cubic coefficient of \(p\) makes the multiplicities of those roots in \(p\) equal, giving \(p=b^2\). Both cases contradict the repeated-root restrictions from (2). This excludes \(L=0\) and proves the theorem.

**Corollary QuinticVeroneseLinearZExclusion.** In characteristic zero
a three-variable degree-five polynomial Keller map cannot have
leading part \(g(x,y,z)N(x,y)\) with a homogeneous cubic \(g\) of
\(z\)-degree one.

The [weighted classification](QUINTIC_VERONESE_LINEAR_Z_WEIGHTED_REDUCTION.md)
leaves only the scalar form and the normalized square form. The
[scalar exclusion](QUINTIC_VERONESE_SCALAR_WEIGHTED_EXCLUSION.md)
and the theorem above cover them respectively.

The [symbolic checker](verify_quintic_veronese_square_exclusion.py)
replays the determinant coefficient, the fibre-coordinate identities
and equation (5). The root, degree and pole arguments are in the proof.
