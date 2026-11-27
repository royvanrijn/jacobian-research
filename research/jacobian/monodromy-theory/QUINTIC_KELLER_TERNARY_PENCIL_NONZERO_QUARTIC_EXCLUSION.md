# A quartic coefficient obstruction in the ternary quintic pencil

**Theorem QuinticKellerTernaryPencilNonzeroQuarticExclusion.**
Let \(F=\sum_{j=0}^5 H_j\) be a polynomial map over a
characteristic-zero field with \(H_j\) homogeneous of degree \(j\) and
\[
 H_5=zN(p,q),\qquad p=z^2,\quad q=x^2+yz,\quad N(p,q)=(p^2,pq,q^2).
\]
If the degree-eleven through degree-eight terms of \(\det JF\)
vanish and the coefficient of \(x^4\) in the first entry of \(H_4\)
is nonzero, then \(\det JF=0\). In particular, any Keller completion
in this leading form has that coefficient zero.

## Normalization and a matching model

Use the notation \(U,V,a,b,c,d,e\) of the
[first two determinant equations](QUINTIC_KELLER_TERNARY_PENCIL_FIRST_EQUATIONS.md).
The coefficient in the statement is \(a\). Over the algebraic closure
scale the source by \((x,y,z)\mapsto(\lambda x,\lambda^2y,z)\)
and the target by \(\operatorname{diag}(1,\lambda^{-2},\lambda^{-4})\).
This preserves \(H_5\) and sends \(a\) to \(\lambda^4a\), so take \(a=1\).

The quartic changes from translating the three source coordinates are
the directional derivatives of \(H_5\). Writing \(DN\) for the
derivative of \(N\) in its two arguments, the derivatives in \(x,y,z\)
are \(DN(p,q)(0,2zx)\), \(DN(p,q)(0,z^2)\), and
\(DN(p,q)(5p/2,q/2+yz)\). A translation in \(z\) first removes
the \(z^2\) coefficient of \(U\); translations in \(x,y\) then
remove the \(xz,z^2\) coefficients of \(V\).
The degree-ten conditions now give \(U=u q\), \(V=v x^2+w yz\).
Set \(Q=q+(w-v)y\). There is a model in \(k[z,Q]^3\) with
the prescribed homogeneous parts in degrees five and four:
\[
 K=
 \begin{pmatrix}
 z^5+Q^2+(2u+b)pQ+cp^2\\
 z^3Q+uQ^2+vpQ+dp^2\\
 zQ^2+2vQ^2+ep^2
 \end{pmatrix}.
\]
The degree-four changes from using \(Q\) in the quintic terms are
\(0,(w-v)z^3y,2(w-v)zqy\), which give the required terms in \(V\).

## Comparing determinant equations

Put \(n=(q^2,-2pq,p^2)\) and \(D=-z\partial_x+2x\partial_y\).
For a homogeneous vector \(R\) of degree \(j\), the contribution
with two columns from \(JH_5\) to determinant degree \(j+7\) is
\(L(R)=5z^3D(n\cdot R)\). As in the first-equations calculation,
\(\ker D=k[z,q]\) and the vectors with zero contraction by \(n\) are
\(DN(p,q)(U',V')\) for homogeneous \(U',V'\) of degree \(j-2\).

Since \(K\) has zero Jacobian determinant, the difference in degree
three lies in \(\ker L\). The contractions of \(\langle z^3,zq\rangle^3\)
span \(z\langle p^3,p^2q,pq^2,q^3\rangle\), the possible invariants
\(n\cdot R\) in degree seven. Tangent vectors with arguments proportional
to \(z\) already lie in \(\langle z^3,zq\rangle^3\). Add arbitrary
multiples of \(z^3,zQ\) to each entry of \(K\) to make the cubic
difference a linear combination of the four tangent directions
\[
 DN(p,q)(x,0),\quad DN(p,q)(y,0),\quad
 DN(p,q)(0,x),\quad DN(p,q)(0,y).
\]
The difference of the degree-nine equations then consists of \(L\)
on the quadratic difference and the mixed determinants with one
column from each of \(JH_5,JH_4\) and the cubic difference.
The [exact checker](verify_quintic_ternary_nonzero_quartic.py)
computes four rational linear functionals on degree-nine monomials
annihilating the quadratic image of \(L\).
On the displayed tangent directions the mixed determinants give
the lower triangular matrix
\[
 \begin{pmatrix}
 -2&0&0&0\\
 0&4&0&0\\
 9b+30u&0&30&0\\
 0&-6b-20u&0&-20
 \end{pmatrix}.
\]
Its determinant is \(4800\). The cubic parts therefore agree.

Now the degree-nine comparison puts the quadratic difference in
\(\ker L\). In degree two this kernel is \(\langle p,q\rangle^3\):
its contractions span all binary cubics in \(p,q\), and the
zero-contraction tangent vectors have constant arguments.
Adding multiples of \(p,Q\) to the model matches degree two.
The degree-eight comparison likewise puts the linear difference in
\(\ker L\). In degree one this kernel is \(\langle z\rangle^3\):
its contractions span \(z\langle p^2,pq,q^2\rangle\), and there
are no nonzero tangent vectors. Adding multiples of \(z\) and
constants matches the remaining terms. Hence the normalized map is
in \(k[z,Q]^3\), and its determinant vanishes identically.
