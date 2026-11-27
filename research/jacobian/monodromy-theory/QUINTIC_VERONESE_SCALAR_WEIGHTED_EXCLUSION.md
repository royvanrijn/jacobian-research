# Exclusion of the scalar weighted linear-in-z Veronese form

**Theorem QuinticVeroneseScalarWeightedExclusion.** The scalar form \((q_2+c)N\) in [QuinticVeroneseLinearZWeightedReduction](QUINTIC_VERONESE_LINEAR_Z_WEIGHTED_REDUCTION.md) is impossible. Thus a degree-five polynomial Keller map with \(H_5=(zq_2+r_3)N(x,y)\), \(q_2\ne 0\), must have the square-factor weighted form \(zN(x^2,xy+\beta z)\) in the stated coordinates, with \(\beta\ne0\).

Work over \(\mathbb C\), as permitted by the weighted reduction. In the scalar form \(F=z^2 A+zB+C\), with \(A_2=cN\) and \(B_4=q_2N\), where subscripts denote binary homogeneous parts and \(N=(x^2,xy,y^2)\).

Suppose first that \(c\ne0\). The coefficient of \(z^5\) in the determinant is \(2\det(A_x,A_y,A)\). Homogenize \(A\) to degree two with homogenizing variable \(s\). Euler's identity makes its homogeneous Jacobian singular. Its three components are independent by restriction to \(s=0\). [Wang–Yan, Lemma 2.4](https://arxiv.org/pdf/2610.36495v1) factors it as a common homogeneous factor times three independent binary forms in a pair of coprime ternary forms. If the degrees of the factor, the binary forms, and the pair are \(r,e,m\), respectively, then \(r+em=2\), with \(e\ge2\), \(m\ge1\). Hence \(r=0,e=2,m=1\). The pair restricts to independent linear forms at \(s=0\). After choosing the basis of the pair, comparison at \(s=0\) shows \(A=cN(x+a,y+b)\). Translate and rescale to \(A=N(x,y)\). With \(\nu=(y^2,-2xy,x^2)\) and \(E=x\partial_x+y\partial_y\), the \(z^4\) determinant coefficient is
\[
  (A_x\times A_y)\cdot B+2(A_y\times A)\cdot B_x+
       2(A\times A_x)\cdot B_y=(6-2E)(\nu\cdot B).
\]
Indeed \(A_x\times A_y=2\nu\) and \(2A=xA_x+yA_y\). The vanishing forces the degree-two part of \(\nu\cdot B\) to be zero, which means \(B(0,0)=0\). Thus the \(z\)-column vanishes along \(x=y=0\).

If \(c=0\), the weighted reduction gives \(\deg_{(x,y)} A\le1\). If \(A\) has rank one as a tuple of polynomials, write it in target coordinates as \(A=(0,0,h)\) and \(B=(b_1,b_2,b_3)\). The degree-four binary leading forms of \(b_1,b_2\) are independent and their Jacobian bracket \(\{b_{1,4},b_{2,4}\}\) is nonzero. The top binary part of the \(z^3\) determinant coefficient is
\[
  2h_d\{b_{1,4},b_{2,4}\}+b_{1,4}\{b_{2,4},h_d\}
     -b_{2,4}\{b_{1,4},h_d\}
   = (2-d/4)h_d\{b_{1,4},b_{2,4}\},
\]
using Euler's identity for the top homogeneous part \(h_d\), of degree \(d=0\) or 1. This part cannot vanish. Rank three is impossible by the \(z^5\) coefficient, since independent affine components have invertible Jacobian after degree-one homogenization.

For rank two, choose \(A=(a,b,0)\) and put \(r=B_3\), where the subscript here denotes the third coordinate. The \(z^4\) coefficient is \(r\{a,b\}+2a\{b,r\}-2b\{a,r\}\). If the linear parts of \(a,b\) are independent, its degree-four part is \((1-8)r_4\{a,b\}\), a contradiction. Otherwise affine source and linear target changes give \(A=(x,1,0)\). Now \(r_y=0\) from the same coefficient. Write
\[
 (P,Q,R)=(xz^2+pz+C_1,\ z^2+qz+C_2,\ rz+C_3).
\]
The \(z^3\) coefficient gives
\[
 (C_3)_y=r' p_y+(r/2-xr')q_y .
\]
The left side has binary degree at most four, and the first summand on the right has degree at most six. Since \(r_4\) is a nonzero multiple of \(x^4\), the degree-seven equation forces \((q_4)_y=0\), contrary to the independence of \(p_4,q_4,r_4\).

It remains to treat \(A=0\). Then \(\det(B_x,B_y,B)=0\) from the \(z^2\) coefficient, and \(B\) has no common zero by the \(z\)-column. Homogenize \(B\) to degree four and apply the same factorization. Its homogeneous common factor restricts at \(s=1\) to a constant because \(B\) has no common zero in the affine plane. The factor is then a scalar, since a positive power of \(s\) would make the degree-four part zero. Consequently \((m,e)=(1,4)\) or \((2,2)\). Two affine linear polynomials with independent gradients have a common zero. If their gradients are dependent, all three degree-four parts of \(B\) are proportional. This excludes \(m=1\).

When \(m=2\), the gcd of the three degree-four parts is the square of the gcd of the two quadratic arguments at \(s=0\). Hence \(q_2\) is a square and the two arguments at \(s=0\) are a common linear factor times a basis of binary linear forms. Choose this common factor as \(x\), and compare degree-four parts after a change of basis in the pair. As in the weighted reduction this gives \(B=N(p,q)\), with \(p=x^2+ax+by+c_0\) and \(q=xy+dx+ey+f\). If \(b\ne0\), substitution of \(p=0\) yields a cubic with nonzero leading coefficient in the equation \(q=0\), giving a common zero. If \(b=0\), every root of \(p\) must be \(-e\), since at any other root the equation \(q=0\) is linear in \(y\) with nonzero slope. Translating \(x,y\) gives \(p=x^2,\ q=xy+c_1\) with \(c_1\ne0\).

On \(x\ne0\) set \(\lambda=q/p\), so \(y=\lambda x-c_1/x\) and \(B=p^2(1,\lambda,\lambda^2)\). The \(z\) coefficient in \((x,\lambda,z)\) of the Jacobian determinant of \(F\) still vanishes, since the change of variables is independent of \(z\). Because \(B_x\) at fixed \(\lambda\) is proportional to \(B\), this coefficient is a nonzero scalar multiple of
\[
 p^4\,\partial_x\big(C_3-2\lambda C_2+\lambda^2 C_1\big)\big|_\lambda .
\]
Thus the expression in parentheses lies in \(\mathbb C[\lambda]\), by the Laurent formula for \(y\). Multiplying by \(p^2=x^4\) gives a polynomial in \(x,y\). Since \(q\) is a unit along \(x=0\) and \(\lambda\) has order \(-2\) there, the polynomial in \(\lambda\) has degree at most two. Subtract a constant target vector to make it zero. Then \(q^2 C_1-2pq C_2+p^2 C_3=0\). Coprimality gives \(C=(2pU,qU+pV,2qV)\) with polynomial \(U,V\).

Set \(v=2pU+zp^2\) and \(W=pV-qU\). Now \(F=(v,\lambda v+W,\lambda^2v+2\lambda W)\), whose determinant in \((v,W,\lambda)\) is \(2W\) up to sign. By the chain rule its original determinant is, up to sign, \(2W p^2\{W,\lambda\}\). The second factor \(p^2\{W,\lambda\}=p\{W,q\}-q\{W,p\}\) is a polynomial. If the product is a nonzero scalar, \(W\) is constant, making it zero. This finishes the affine case.
