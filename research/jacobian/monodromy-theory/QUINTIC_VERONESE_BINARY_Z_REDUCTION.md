# The quadratic \(z\)-coefficient for a binary quintic Veronese factor

**Theorem QuinticVeroneseBinaryZReduction.** Let \(F\) be a three-variable Keller map over \(\mathbb C\) of degree five whose ordinary leading part is \(H_5=g(x,y)N(x,y)\), where \(g\ne0\) is homogeneous cubic and \(N=(x^2,xy,y^2)\). Then
\[
                       F=C(x,y)+z B(x,y)+z^2 A(x,y),
\]
where \(\deg B\le3\) and \(A\) is nonzero affine. Up to affine changes of the two binary source variables and a linear target change, \(A\) is one of
\[
                          (0,0,h(x,y)),\qquad (x,1,0),
\]
with \(h\) a nonzero affine polynomial. In the second case the third entry of \(B\) depends only on \(x\). The [rank exclusions](QUINTIC_VERONESE_BINARY_RANK_EXCLUSIONS.md) leave only nonconstant rank one. The result extends to characteristic zero by embedding the finitely generated coefficient field in \(\mathbb C\).

## The high determinant coefficients

Write \(\nu=(y^2,-2xy,x^2)\) and let \(\Delta_d\) be the degree-\(d\) homogeneous part of \(\det JF\). The
[adjugate coefficient identity](QUINTIC_VERONESE_LEADING_RESIDUAL.md#the-maximal-z-degree-case)
gives \(5g^2\partial_z(\nu\cdot H_j)\) for the contribution of type \((5,5,j)\), since \(g_z=0\). Thus \(\Delta_{11}=0\) makes every positive-\(z\) coefficient of \(H_4\) tangent to \(N\). Such a column is \(W=DN(x,y)(U,V)\) by the
[polynomial syzygy](QUINTIC_VERONESE_LINEAR_Z_WEIGHTED_REDUCTION.md#from-ordinary-to-weighted-degree).
In particular \([z^4]H_4=0\).

Put \(D=xV-yU\) for a tangent column \(W\). Directly from the formulas for \(N,\nu\),
\[
 \begin{aligned}
 \nu\cdot W_x&=-2yD,& \nu\cdot W_y&=2xD,\\
 (gN)_x\times W&=(2gV+g_xD)\nu,&
 (gN)_y\times W&=(-2gU+g_yD)\nu .
 \end{aligned}
\]
Consequently
\[
 \det((gN)_x,W_y,W)+\det(W_x,(gN)_y,W)=-10g D^2 ,
 \tag{1}
\]
using \(xg_x+yg_y=3g\).

In \(\Delta_{10}\) the other contribution besides \((5,5,3)\) is \((5,4,4)\). Its \(z^5\) coefficient is three times (1) for \(W=[z^3]H_4\). Here \(U,V\) are constant, so \(D=0\) forces \(W=0\). The coefficient of \(z^3\) is now twice (1) for \(W=[z^2]H_4\). Here \(U,V\) are linear and \(D=0\) gives \(U,V\) proportional to \(x,y\) by a common constant. Hence \([z^2]H_4=cN\). There is no \((5,5,3)\) contribution at either of these powers because it has \(z\)-degree at most two.

At \(z^2\), the \((5,4,4)\) terms use one \(cN\) and one copy of \(W=[z]H_4\). Terms with an \(x\) or \(y\) derivative on \(cN\) are all tangent to \(N\). In the other two terms the cross products \((gN)_x\times cN=cgy\nu\) and \((gN)_y\times cN=-cgx\nu\), paired with \(\nu\cdot W_y=2xD,\ \nu\cdot W_x=-2yD\), give cancelling contributions. Therefore \(\Delta_{10}\) at \(z^2\) gives \(15g^2\nu\cdot[z^3]H_3=0\). The column \([z^3]H_3\) is constant, so it vanishes. This yields the quadratic form in \(z\), with degree-two part \(A_2=cN\).

If \(c\ne0\), the \(z^5\) coefficient of the full determinant gives \(\det(A_x,A_y,A)=0\). Homogenize \(A\) to degree two and apply Wang–Yan [Lemma 2.4, printed pages 6–7](https://arxiv.org/pdf/2610.36495v1). The components are independent at infinity, and the independent homogeneous degree-two factorization is a basis of binary quadratics in two linear forms. Comparing at infinity with \(cN\) gives \(A=cN(x+\alpha,y+\beta)\). Translate \(x,y\) to make \(A=cN\). The \(z^4\) coefficient is then
\[
                    c^2(6-2E)(\nu\cdot B),\qquad E=x\partial_x+y\partial_y.
\]
Indeed it comes from two uses of \(cN\) and one use of \(B\), with coefficients given by \(N_x\times N_y=2\nu,\ N_y\times N=-x\nu,\ N\times N_x=-y\nu\). The degree-two part of this coefficient forces \(\nu\cdot B(0,0)=0\), hence \(B(0,0)=0\). This makes \(F_z(0,0,z)=0\), a contradiction. Thus \(c=0\).

One further equation in the original variables is available. With \(B_3=DN(P_2,Q_2)\) and \(A_1\) the linear homogeneous part of \(A\), the \(z\) coefficient in \(\Delta_{10}\) is the sum of \(10g^2\nu\cdot A_1\) and (1) for \(W=B_3\). It gives
\[
                          g\,\nu\cdot A_1=(xQ_2-yP_2)^2 .
 \tag{2}
\]

## The affine-in-\(z\) alternatives

Suppose \(A=0\). Since \(F_z=B\), the three entries of \(B\) have no common affine zero. Their span as polynomials has dimension one, two or three. In dimension one, \(B\) is a constant vector, since a nonconstant binary polynomial over \(\mathbb C\) has a zero. Put \(B=(0,0,1)\) by a target change. The planar determinant of the first two entries of \(C\) would be constant and nonzero, but their independent degree-five homogeneous leading parts have nonzero binary Jacobian. Indeed the Jacobian of two homogeneous forms of the same positive degree vanishes only if their ratio is constant, by Euler's identity. This excludes dimension one.

In dimension two, put \(B=(p,q,0)\), so \(p,q\) are coprime and \(\deg p,\deg q\le3\). The third entry \(C_3\) of \(C\) has degree five because the entries of \(H_5\) are independent. It has no critical point by the Keller condition. The coefficient of \(z\) in the determinant is \(p\{q,C_3\}-q\{p,C_3\}=p^2\{q/p,C_3\}\), where the brackets denote the binary Jacobian. In characteristic zero its vanishing means \(q/p\) is algebraic over \(\mathbb C(C_3)\).

A complex polynomial without critical points is non-composite. In fact if \(C_3=S(T)\) for a polynomial \(T\) and a univariate polynomial \(S\) of degree at least two, \(T\) attains a root of \(S'\) and the gradient of \(C_3\) vanishes there. Arzhantsev and Petravchuk [Theorem 1, (i) and (iii)](sources/arzhantsev-petravchuk-math-0608157v2.pdf) show that a non-composite polynomial is closed, meaning that the polynomial algebra it generates is integrally closed in the ambient polynomial ring. Petravchuk and Iena [Remark 2](sources/petravchuk-iena-math-0701588v2.pdf) identify this with relative algebraic closedness of the generated rational-function field in the ambient rational-function field. Thus \(q/p\in\mathbb C(C_3)\). Write it as \(J(C_3)/K(C_3)\) for coprime univariate polynomials \(J,K\). A Bézout identity shows that the substituted numerator and denominator are still coprime. If the ratio is nonconstant, one of \(p,q\) must therefore have degree at least \(\deg C_3=5\), which is impossible. A constant ratio contradicts dimension two.

In dimension three, the \(z^2\) coefficient gives \(\det(B_x,B_y,B)=0\). Homogenize to the actual degree \(d\le3\) and apply Wang–Yan Lemma 2.4 again. Its common homogeneous factor specializes to a constant at the affine chart because \(B\) has no common zero. A positive-degree homogeneous factor with constant specialization would be a power of the homogenizing variable, which is impossible at the actual degree. The degree constraint \(em=d\) with \(e\ge2\) then gives \(m=1\). The two affine-linear arguments have dependent linear parts, since otherwise they would have a common affine zero. Hence the entries of \(B\) depend on a single affine-linear variable. If \(d=3\), the tangent column \(B_3\) is a constant vector times the cube of a binary linear form, and \(\nu\cdot B_3=0\) would force that vector to be zero. Thus \(d=2\), and a linear target change gives \(B=(1,x,x^2)\) in an affine binary frame.

The \(z\) coefficient of the determinant in this frame is the \(y\) derivative of \(C_3-2xC_2+x^2C_1\), so this scalar equals \(T(x)\). With \(v=C_1+z\) and \(W=C_2-xC_1\), the map becomes
\[
                   (v,\ xv+W,\ x^2v+2xW+T(x)).
\]
The Jacobian in \(x,y,v\) is \(-(2W+T')W_y\), which cannot be a nonzero constant. It would make both factors units, so \(W_y=0\). Thus \(A=0\) is impossible in every case.

## The two remaining affine forms

The \(z^5\) coefficient \(2\det(A_x,A_y,A)\) excludes dimension three for the span of the affine entries of \(A\). If the span has dimension one, a target change gives \(A=(0,0,h)\) as claimed.

For a two-dimensional span, write \(A=(a,b,0)\) and let \(r\) be the third entry of \(B\). The \(z^4\) coefficient is
\[
                         r\{a,b\}+2a\{b,r\}-2b\{a,r\}.
\]
If the linear parts of \(a,b\) are independent, choose \(A=(x,y,0)\). The coefficient becomes \((1-2E)r\), forcing \(r=0\). Then \(F_3=C_3(x,y)\) and the \(z^3\) coefficient is \(2E C_3\), contrary to \(\deg C_3=5\). If the linear parts are dependent, the span has a constant direction and one may choose \(A=(x,1,0)\). The \(z^4\) coefficient is then \(-2r_y\), making \(r\) a polynomial in \(x\).
