# Linear arguments in independent quintic Keller leading parts

**Theorem QuinticKellerLinearArgumentLeadingExclusion.** Let \(F\)
be a degree-five three-variable Keller map in characteristic zero.
Its leading part cannot be \(g\Phi(x,y)\) in linear source
coordinates, where \(\Phi\) consists of linearly independent
coprime homogeneous binary forms of degree \(e\ge2\) and \(g\)
is homogeneous of degree \(5-e\).

Work over \(\mathbb C\) by embedding the finitely generated
coefficient field. The [Veronese exclusion](QUINTIC_VERONESE_BINARY_RANK_EXCLUSIONS.md)
covers \(e=2\). The [independent binary exclusion](QUINTIC_KELLER_INDEPENDENT_BINARY_LEADING_EXCLUSION.md)
covers \(g\) independent of \(z\), including \(e=5\).
It remains to consider \(\deg_z g=1\) for \(e=3,4\)
and \(\deg_z g=2\) for \(e=3\).

## Mixed highest coefficients

Write \(n=\Phi_x\times\Phi_y\), and set \(S=n\cdot H_j\)
for the degree-\(j\) homogeneous part of \(F\). The trace
term in the determinant using two copies of the leading
Jacobian is
\[
                  \frac{g}{e}(5g S_z-jg_z S).
                  \tag{1}
\]
Indeed Euler's identity for \(\Phi\) gives the rows of the
adjugate of \(J(g\Phi)\) in the form
\[
       \frac{g}{e}(-xg_z,-yg_z,5g-zg_z)^{\mathsf T} n^{\mathsf T}.
\]
For \(j=4\), the integer \(5\deg_z S-4\deg_z g\) cannot
be zero since \(\deg_z g\) is one or two. Thus the
degree-eleven determinant coefficient forces \(n\cdot H_4=0\).
In particular \([z^4]H_4=0\).

Let \(s=\deg_z g\), \(P=[z^s]H_5\), and \(W=[z^3]H_4\).
The \(z^{s+5}\) term of the degree-ten determinant comes only
from the triple \((5,4,4)\); the \((5,5,3)\) triple has
\(z\)-degree at most \(2s+2<s+5\). Its coefficient is
\[
 3\det(P_x,W_y,W)+3\det(W_x,P_y,W)+s\det(W_x,W_y,P)
                =(4s-15)P\cdot(W_x\times W_y)
\]
by the mixed Euler identity for linear \(W\) and \(\deg P=5-s\).
A rank-one \(W\) cannot be tangent to \(n\), since \(n\) has
no constant annihilator. A rank-two \(W\) gives a nonzero
constant \(W_x\times W_y\) and the displayed equation would
contradict independence of the entries of \(P\). Thus \(W=0\).

## The quadratic \(z\)-factor

If \(s=2\), by a linear source shear in \(z\) and scaling
take \(g=z^2+g_2(x,y)\). The normalization preserves the
previous argument. The \(z^6\) term of the degree-ten
determinant now comes only from \((5,5,3)\) and by (1)
is a nonzero multiple of \(n\cdot[z^3]H_3\). Hence
\([z^3]H_3=0\), since this column is constant.
Write \(F=z^2 A+zB+C\), where \(A\) has primitive cubic
leading part \(\Phi\), and \(\deg B\le3,\deg C\le5\).

The \(z^5\) determinant coefficient is
\(2\det(A_x,A_y,A)=0\). Homogenize \(A\) to degree three
in \((x,y,w)\). Euler's identity converts this equation
to a singular Jacobian for the independent homogeneous
tuple. [Wang and Yan, Lemma 2.4](https://arxiv.org/pdf/2610.36495v1)
factors it as a common factor times independent binary
forms in coprime homogeneous forms. By primitiveness
at \(w=0\), the common factor is constant. The remaining
degrees must be three in linear arguments, independent
at \(w=0\). Comparing there with \(\Phi\) and translating
the binary source variables gives \(A=\Phi(x,y)\).

The \(z^4\) coefficient is then
\[
                    n\cdot B-\frac23 n\cdot E B,\qquad
                    E=x\partial_x+y\partial_y .
\]
This follows by substituting \(A=(xA_x+yA_y)/3\) in the
two determinant terms using \(B_x,B_y\), with
\(n=A_x\times A_y\). Its binary degree-four part is
\(n\cdot B(0,0)\). Thus \(B(0,0)=0\) by the no-constant-
annihilator property, and \(F_z\) vanishes on \(x=y=0\).
This contradicts the Keller condition.

## The linear \(z\)-factor

Let \(s=1\). Set \(z=t^2\) after scaling \(x,y\) by \(t\).
The coefficient of the highest power \(t^6\) is
\[
       T=P+[z^2]H_4+[z^3]H_3,\qquad P=[z]H_5=g_1\Phi .
\]
It has independent components and is even in \(x,y\). All
higher powers are absent by the preceding mixed coefficients.
The Jacobian in \((x,y,t)\) is \(2\kappa t^3\).
Its \(t^{17}\) coefficient gives \(\det(T_x,T_y,T)=0\).
Homogenize \(T\) to degree four and apply Lemma 2.4 again.
Its binary leading part is \(P\), with common-factor degree
\(4-e\), zero or one.

If the arguments in the homogeneous factorization have
degree two, the outer forms are a basis of binary quadrics
and there is no common factor before restriction to infinity.
The binary leading arguments have independent restrictions
since the components of \(P\) are independent. Their common
factor contributes its square to the common factor of
\(P\). Hence this case requires \(e=4\) and
\(P=N(p_0,q_0)\), with \(p_0,q_0\) coprime binary quadrics.
After a linear shear in \(z\), the original quintic leading
part is \(zN(p_0,q_0)\). This is excluded by the
[external binary pencil theorem](QUINTIC_KELLER_EXTERNAL_BINARY_CONIC_PENCIL_EXCLUSION.md).

Otherwise the arguments have degree one and independent
restrictions at infinity. Choosing a basis of these arguments
and comparing the highest binary terms gives
\[
                         T=h(x,y)\Phi(x+a,y+b)
\]
with \(\deg h=4-e\). Evenness makes the projective directions
of \(\Phi(x+a,y+b)\) and \(\Phi(x-a,y-b)\) coincide.
If \((a,b)\ne(0,0)\), change binary coordinates so that
\(a\ne0,b=0\). The two slopes \(y/(x+a)\) and \(y/(x-a)\)
are algebraically independent: they recover \(x,y\)
rationally. A nonconstant projective map in one slope cannot
give the same value at these independent slopes. This
forces \(a=b=0\). Then \(h\) has the same parity as
\(\Phi\), so it is homogeneous because \(\deg h=0\) or one.
It follows that \(T=P\).

The last equality kills \([z^2]H_4\) and \([z^3]H_3\),
which have different binary degrees. Now \(F=z^2 A+zB+C\)
with \(A\) affine and \(B_4=P\) linearly independent.
We include the affine-coefficient argument since it uses
only this independence.

An affine tuple \(A\) with three independent entries is
excluded by the \(z^5\) coefficient \(\det(A_x,A_y,A)=0\).
If \(A=(0,0,h)\) has rank one, let \(p,q\) be the first
two entries of \(B\). The top binary part of the \(z^3\)
coefficient is
\((2-d/4)h_d\{p_4,q_4\}\), with \(d=\deg h\le1\),
using binary Euler's identity. It is nonzero. If the
rank is two, the independent-linear-part case normalizes
\(A\) to \((x,y,0)\). The \(z^4\) coefficient is then
\((1-2E)r\), where \(r\) is the third entry of \(B\),
forcing \(r=0\), contrary to its nonzero quartic part.
With dependent linear parts, take \(A=(x,1,0)\).
Now \(z^4\) forces \(r_y=0\), and \(z^3\) gives
\[
                   (C_3)_y=r' p_y+(r/2-xr')q_y .
\]
The left side has degree at most four and the first
summand at most six. The degree-seven term therefore
forces \((q_4)_y=0\) since \(r_4\ne0\) is proportional
to \(x^4\). This contradicts independence of the
quartic parts of \(p,q,r\).

Thus \(A=0\). The tuple \(B\) has no common zero,
and \(\det(B_x,B_y,B)=0\) by the \(z^2\) coefficient.
Homogenization and Lemma 2.4 apply once more. The
common homogeneous factor is constant since its
affine specialization has no zero and \(B\) has
actual degree four. With linear arguments the
degrees of the outer forms must be four. The leading
arguments at infinity are independent by independence
of \(B_4\), so their affine forms have a common zero,
a contradiction. With quadratic arguments the outer
forms are a basis of binary quadrics. The common
factor of \(B_4\) is again the square of the common
factor of the two leading arguments, so its degree
cannot be one. In the degree-zero case the two
quadratic arguments have no common point at infinity.
By Bézout their homogenizations must meet in the
affine plane, again giving a common zero of \(B\).
This finishes the proof.
