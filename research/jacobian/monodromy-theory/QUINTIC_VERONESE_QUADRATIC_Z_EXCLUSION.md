# Excluding the quadratic-in-\(z\) quintic Veronese class

**Theorem QuinticVeroneseQuadraticZExclusion.** There is no three-variable polynomial map over a characteristic-zero field with nonzero constant Jacobian determinant, ordinary degree five and leading part \(H_5=g(x,y,z)(x^2,xy,y^2)\) with \(g\) a homogeneous cubic of \(z\)-degree two.

It suffices to work over \(\mathbb C\) since the finitely generated subfield containing the coefficients of any example embeds in \(\mathbb C\). Write \(N=(x^2,xy,y^2)\), \(\nu=(y^2,-2xy,x^2)\), and \(\ell=[z^2]g\ne0\), a binary linear form. Write \(DN(P,Q)=P N_x+Q N_y\).

## Bounding the whole map's \(z\)-degree

For a homogeneous column \(H_j\) of degree \(j\), the contribution to the determinant using twice the leading part and once \(H_j\) is
\[
 \operatorname{tr}(\operatorname{adj}(JH_5)JH_j)
      = g(5g S_z-jg_z S),\qquad S=\nu\cdot H_j .
\]
Indeed the adjugate rows from \(N_x\times N_y=2\nu\), \(N_y\times N=-x\nu\), \(N\times N_x=-y\nu\) give \(g(5g e_3-g_z(x,y,z)^{\mathsf T})\nu^{\mathsf T}\). Euler's identity on \(H_j\) yields the formula. A nonzero \(S\) of \(z\)-degree \(e\) has top contribution \((5e-2j)\ell^2[z^e]S\) in degree \(e+3\) if \(5e-2j\ne0\).

The ordinary degree-eleven determinant part is entirely of type \((5,5,4)\). Since \(5e-8\ne0\), it gives \(\nu\cdot H_4=0\). The syzygy \(\nu\cdot H=0\) implies \(H=DN(P,Q)\): reduction modulo \(x\) and \(y\) gives \(H_1=2xP,\ H_3=2yQ\), and then \(H_2=yP+xQ\). Hence \(H_4\) has \(z\)-degree at most three and \(B=[z^3]H_4=DN(a,b)\) with constants \(a,b\).

At ordinary degree ten, type \((5,5,3)\) has \(z\)-degree at most six by the adjugate formula. Type \((5,4,4)\) has degree at most seven, with the coefficient at \(z^7\) furnished by \(z^2A\), \(A=\ell N\), and two uses of \(z^3B\). This coefficient is
\[
 3\det(A_x,B_y,B)+3\det(B_x,A_y,B)+2\det(B_x,B_y,A)
   =-7 A\cdot(B_x\times B_y)
   =-14\ell(bx-ay)^2 .
\]
For the first equality, \(B=xB_x+yB_y\) and \(x A_x+y A_y=3A\). The last equality uses \(B_x\times B_y=2(b^2,-2ab,a^2)\). Vanishing forces \(a=b=0\), so \(H_4\) has \(z\)-degree at most two.

Now the \((5,4,4)\) contribution has \(z\)-degree at most five. The \((5,5,3)\) formula in degree six therefore forces \([z^3](\nu\cdot H_3)=0\). Since \([z^3]H_3\) is a constant vector and the entries of \(\nu\) are linearly independent, \([z^3]H_3=0\). All homogeneous parts now have \(z\)-degree at most two.

## Centering the quadratic coefficient

Write the full map as \(F=z^2 A_2(x,y)+zB_1(x,y)+C_0(x,y)\). The cubic leading part of \(A_2\) is \(\ell N\). The coefficient at \(z^5\) in \(\det JF\) is \(2\det(A_{2,x},A_{2,y},A_2)\), hence it vanishes. Let \(T(x,y,t)\) be the homogenization of \(A_2\) to degree three. At \(t=1\), Euler's identity gives \(\det JT=3\det(A_{2,x},A_{2,y},A_2)=0\), so \(\det JT=0\) identically by homogeneity. The components are linearly independent from their restrictions \(T(x,y,0)=\ell N\).

Apply Wang–Yan [Lemma 2.4, printed pages 6–7](https://arxiv.org/pdf/2610.36495v1). Its common-degree statement gives \(T_i=d\,\Phi_i(p,q)\), with \(p,q\) coprime homogeneous forms of a common degree \(m\ge1\), \(\Phi_i\) coprime binary forms of a common degree \(e\ge2\), and \(d\) homogeneous of degree \(r\) with \(r+em=3\). Thus \(m=1\). The restrictions of \(p,q\) to \(t=0\) are independent because the \(T_i(x,y,0)\) are independent and would all be proportional if \(p,q\) restricted to proportional linears. Therefore the restricted \(\Phi_i(p,q)\) are coprime, and \(d(x,y,0)\) is a nonzero scalar multiple of the common factor \(\ell\). Hence \(r=1\) and \(e=2\). Change basis in the binary arguments so that \(p=x+\alpha t,\ q=y+\beta t\), and absorb the leading scalar of \(d\) into the binary forms. Comparing \(T(x,y,0)=\ell N(x,y)\) then gives \(T=(\ell(x,y)+\gamma t)N(x+\alpha t,y+\beta t)\). Translating \(x,y\) after \(t=1\) makes
\[
 A_2=hN,\qquad h=\ell+c
\]
for a constant \(c\), without changing the Keller condition.

## The vanishing derivative

Let \(S=\nu\cdot B_1\) and \(E=x\partial_x+y\partial_y\). The coefficient of \(z^4\) in the determinant is
\[
 \det(A_{2,x},A_{2,y},B_1)+2\det(B_{1,x},A_{2,y},A_2)
       +2\det(A_{2,x},B_{1,y},A_2)
   =h\big((7\ell+6c)S-2hES\big).
\]
Indeed the three contributions are \(h(2h+Eh)S\), \(-2h^2 x \nu\cdot B_{1,x}\) and \(-2h^2 y \nu\cdot B_{1,y}\). Use \(E\nu=2\nu\) to combine them. If \(S\ne0\), its top homogeneous part \(S_s\) of degree \(s\) contributes \((7-2s)\ell S_s\) to the expression in parentheses in degree \(s+1\). No integer \(s\) makes this vanish. Thus \(S=0\). The polynomial syzygy gives \(B_1=DN(P,Q)\), which vanishes at \(x=y=0\). Also \(A_2(0,0)=0\), so \(F_z(0,0,z)=0\). This contradicts the constant nonzero Jacobian determinant.

The [symbolic checker](verify_quintic_veronese_quadratic_z.py) replays the adjugate and two mixed determinant identities with SymPy 1.14.0 over rational coefficients. The degree arguments and the use of homogeneous factorization are in the written proof.

~~~sh
.venv/bin/python research/jacobian/monodromy-theory/verify_quintic_veronese_quadratic_z.py
~~~
