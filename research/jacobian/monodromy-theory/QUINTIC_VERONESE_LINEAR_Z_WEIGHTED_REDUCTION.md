# Weighted reduction of the linear-in-\(z\) quintic Veronese class

**Theorem QuinticVeroneseLinearZWeightedReduction.** Suppose that a complex polynomial map \(F\) in three variables has nonzero constant Jacobian determinant, ordinary degree five, and leading part \(H_5=g(x,y,z)(x^2,xy,y^2)\) with \(g\) a homogeneous cubic of \(z\)-degree one. Put \(q_2=[z]g\ne0\) and \(N=(x^2,xy,y^2)\). Then \([z^3]H_4=0\), and the coefficient of \(t^6\) in \(F(t x,t y,t^2)\) lies in one of the following two classes.

- It is \((q_2+c)N\) for a constant \(c\). In this class \(F\) has \(z\)-degree at most two and \([z^2]H_4=cN\).
- After a binary linear change and its corresponding linear action on the three target coordinates, it is \(N(x^2,xy+\beta)\), with \(q_2=x^2\) and \(\beta\ne0\). In this class
  \[
   [z^2]H_4=(0,\beta x^2,2\beta xy),\qquad
   [z^3]H_3=(0,0,\beta^2).
  \]
  In particular the case of nonzero \([z^3]H_3\) is possible only for square \(q_2\).

Any example over a characteristic-zero field gives such a complex example by embedding the finitely generated subfield of coefficients. The [scalar exclusion](QUINTIC_VERONESE_SCALAR_WEIGHTED_EXCLUSION.md) and [square exclusion](QUINTIC_VERONESE_SQUARE_EXCLUSION.md) rule out both classes.

## From ordinary to weighted degree

Write \(\nu=(y^2,-2xy,x^2)\) and use the
[adjugate coefficient formula](QUINTIC_VERONESE_QUADRATIC_Z_EXCLUSION.md#bounding-the-whole-maps-z-degree)
\(g(5g S_z-jg_z S)\), where \(S=\nu\cdot H_j\), for the contribution of type \((5,5,j)\). If \(S\ne0\) has \(z\)-degree \(e\), the coefficient at \(z^{e+1}\) is \(q_2^2(5e-j)[z^e]S\). The ordinary degree-eleven determinant part shows that \(\nu\cdot H_4=0\). The same polynomial syzygy gives \([z^3]H_4=DN(a,b)\) for constants \(a,b\), and \(H_4\) has no \(z^4\) term.

At ordinary degree ten, type \((5,5,3)\) has \(z\)-degree at most four. In type \((5,4,4)\), \(z^6\) comes only from \(zA\) and two uses of \(z^3B\), where \(A=q_2N\), \(B=DN(a,b)\). Its coefficient is
\[
 3\det(A_x,B_y,B)+3\det(B_x,A_y,B)+\det(B_x,B_y,A)
   =-22q_2(bx-ay)^2 .
\]
Here \(B\) has binary degree one, \(A\) has binary degree four and \(B_x\times B_y=2(b^2,-2ab,a^2)\). Since the coefficient vanishes, \(a=b=0\).

In \(F(t x,t y,t^2)\) the highest \(t\)-degree is now six, with coefficient
\[
 P=q_2N+[z^2]H_4+[z^3]H_3 .
\]
The determinant in \((x,y,t)\) is \(2\kappa t^3\), where \(\kappa=\det JF\ne0\). The coefficient at \(t^{17}\) gives \(6\det(P_x,P_y,P)=0\). Homogenize \(P\) to degree four in \((x,y,s)\), obtaining \(T\). Euler's identity gives \(\det JT(x,y,1)=4\det(P_x,P_y,P)=0\). The restriction \(T(x,y,0)=q_2N\) has linearly independent components.

Apply Wang–Yan [Lemma 2.4, printed pages 6–7](https://arxiv.org/pdf/2610.36495v1). Write \(T_i=h\Phi_i(p,q)\) with \(p,q\) coprime forms of common degree \(m\), the \(\Phi_i\) coprime binary forms of common degree \(e\ge2\), and \(h\) homogeneous of degree \(r\) with \(r+em=4\). The restrictions \(p_0,q_0\) to \(s=0\) are independent, since otherwise \(T(x,y,0)\) would have proportional components.

If \(m=1\), the specialized forms \(\Phi_i(p_0,q_0)\) are coprime and their common factor in \(T(x,y,0)\) is exactly \(h(x,y,0)=q_2\) up to a scalar. Thus \(r=e=2\). Changing basis in \(p,q\) and comparing at \(s=0\) gives \(P=h(x,y,1)N(x+\alpha,y+\beta)\). But \(P\) is even under \((x,y)\mapsto(-x,-y)\). The projective ratio of its first two components is the ratio of the two affine arguments of \(N\). Proportionality of these arguments with those after sign change forces \(x\beta-y\alpha=0\) identically, so \(\alpha=\beta=0\). Therefore \(h(x,y,1)\) is even and \(P=(q_2+c)N\), giving the first class. The absence of a \(z^3\) term in both \(H_3\) and \(H_4\) makes the full map quadratic or affine in \(z\).

If \(m=2\), then \(r=0,e=2\). The common factor of \(\Phi_i(p_0,q_0)\) is exactly \(\gcd(p_0,q_0)^2\) up to a scalar. Indeed after dividing \(p_0,q_0\) by their gcd, no projective binary point can be a simultaneous zero of the divided pair or of all three substituted \(\Phi_i\). Thus \(q_2=\ell^2\) for a binary linear form \(\ell\), and \(p_0,q_0\) are \(\ell\) times a basis of binary linears. Normalize that basis and compare at \(s=0\), giving \(P=N(p(x,y,1),q(x,y,1))\) with top parts \(\ell x,\ell y\). Evenness of \(P\) forces both arguments to be even. For example, the first coordinate is the square of the first argument, so that argument after sign change is equal to itself or to its negative, and its nonzero even quadratic part rules out the negative. After a binary linear change putting \(\ell=x\), followed by the induced inverse action on the target Veronese, the arguments are \(x^2+\alpha,xy+\beta\). If both constants vanish, this is included in the first class. The next argument rules out \(\alpha\ne0\).

## Excluding a simple common zero of the quadratic arguments

Suppose \(\alpha\ne0\). The two arguments \(x^2+\alpha,xy+\beta\) have a common zero with \(x\ne0\), and their Jacobian determinant there is \(2x^2\ne0\). Take these arguments as formal local coordinates \(u,v\) at that point using the invertible formal Jacobian. The transformed family \(F(t x,t y,t^2)\) now has leading term \(t^6N(u,v)\) in \(\mathbb C[[u,v]]((t^{-1}))\). Its Jacobian determinant in \((u,v,t)\) is \(2\kappa t^3\) times a unit whose constant in \(u,v\) is nonzero.

For a coefficient \(t^{6-r}M(u,v)\), the linear contribution with two leading terms to the determinant at order \(t^{17-r}\) is
\[
 (12-2r)\nu\cdot M-6\nu\cdot(u\partial_u+v\partial_v)M,
 \qquad \nu=(v^2,-2uv,u^2).
\]
This follows from the same three cross products for \(N\). On a Taylor-homogeneous piece \(M_j\) of degree \(j\), the coefficient is \((12-2r-6j)\nu\cdot M_j\). For \(r=1,2,3,4,5\), the only resonance \(12-2r-6j=0\) occurs at \(r=3,j=1\). In this range two uses of any retained term with index three cannot contribute, since their order corresponds to index at least six. Inductively remove each nonresonant piece by the completed syzygy \(M_j=DN(U_j,V_j)\) and the near-identity precomposition
\[
 (u,v,t)\longmapsto(u-t^{-r}U,\ v-t^{-r}V,\ t).
\]
Here \(U,V\) are the sums of the homogeneous syzygy coefficients of the nonresonant part, convergent in the \((u,v)\)-adic topology. Substitution is also defined \(t^{-1}\)-adically. It changes the source Jacobian by a unit of the form \(1+O(t^{-1})\), preserving the leading determinant order at most \(t^3\). These steps give
\[
             t^6N(u,v)+t^3L(u,v)+C(u,v)+O(t^{-1}),
\]
with \(L\) homogeneous linear and \(C\) formal.

Now substitute \(u=t^{-3}U,\ v=t^{-3}V\). The resulting family has coefficients in \(\mathbb C[U,V][[t^{-1}]]\): each Taylor degree \(j\) in a coefficient at \(t^{-i}\) contributes only to order \(t^{-i-3j}\), so only finitely many terms affect an individual coefficient. Its leading coefficient is \(A=N(U,V)+L(U,V)+C(0,0)\). Its determinant has leading term \(c_0t^{-3}\) with \(c_0\ne0\), because the additional source determinant is \(t^{-6}\) and the former local unit now has nonzero constant. If \(D_1\) is the coefficient at \(t^{-1}\), the determinant coefficient at \(t^{-2}\) gives \(\det(A_U,A_V,D_1)=0\).

The differential of every affine quadratic map \(A=N+U b+V d+c\) has a finite, nonempty critical set. To see this, test \(\lambda A_U+\mu A_V=0\). If \(d_1=0\), the column \(A_V\) vanishes at \(U=-d_2,\ V=-d_3/2\). If \(b_3=0\), \(A_U\) vanishes at \(U=-b_1/2,\ V=-b_2\). For \(\lambda=1,\ \mu=s\ne0\), the first and last coordinates give
\[
 U=-\frac{b_1+s d_1}{2},\qquad V=-\frac{b_3/s+d_3}{2}.
\]
The middle coordinate is equivalent to
\[
 -d_1s^3+(2d_2-b_1)s^2+(2b_2-d_3)s-b_3=0.
\]
If \(b_3,d_1\) are both nonzero, this has a nonzero complex root. If the cubic vanishes identically, then \(b_3=d_1=0\) and every \(s\ne0\) gives the same point. Otherwise there are only finitely many roots. This proves the claim.

In particular the minors of \((A_U,A_V)\) have no common irreducible divisor. The equation on \(D_1\) first gives \(D_1=A_U E+A_V G\) over the fraction field of \(\mathbb C[U,V]\). At each irreducible divisor some minor is a unit in the associated discrete valuation ring, so Cramer's rule places both \(E,G\) in that ring. Since \(\mathbb C[U,V]\) is a UFD, the coefficients are polynomials. Precompose at order \(t^{-1}\) by \((U-t^{-1}E,V-t^{-1}G,t)\) to remove \(D_1\). The source Jacobian is \(1+O(t^{-1})\). Now the determinant coefficient at \(t^{-3}\) is \(-2\det(A_U,A_V,D_2)\) for a polynomial \(D_2\), which vanishes at a critical point of \(A\). This contradicts \(c_0\ne0\), finishing the exclusion of \(\alpha\ne0\).

The [symbolic checker](verify_quintic_veronese_linear_z.py) replays the mixed determinant identities, the finite resonance calculation and the critical-point cubic with SymPy 1.14.0. The weighted factorization, formal substitutions and denominator argument are in the written proof.

~~~sh
.venv/bin/python research/jacobian/monodromy-theory/verify_quintic_veronese_linear_z.py
~~~
