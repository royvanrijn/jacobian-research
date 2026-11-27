# The tangent double conic in a mixed rank-two quintic pencil

**Theorem QuinticKellerRankTwoTangentDoubleConic.** No Keller map in three
variables over a field of characteristic zero has homogeneous leading parts
\[
             P_5=zq^2,\qquad Q_5=z^3q,\qquad R_4=q^2,
             \qquad q=x^2+yz .
\]
The same conclusion holds when \(q\) is a geometrically irreducible ternary
quadratic whose projective conic is tangent to \(z=0\). Indeed over
\(\mathbb C\) its binary restriction is a nonzero rank-one form. After making
that restriction \(x^2\), the coefficient of \(yz\) is nonzero by
nondegeneracy, and redefining \(y\) yields the displayed normal form.
It suffices to prove the theorem over \(\mathbb C\), by embedding the finitely
generated coefficient field.

This treats the tangent case with both pencil parameters zero in the
[generic conic-argument pencil](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_GENERIC_PENCIL.md).
Wang and Yan also use the double-conic normal form in the proof of their
degree-four classification. The Keller equations below use degrees
\((5,5,4)\). [Wang and Yan, §3.4.3](https://arxiv.org/pdf/2610.36495v1).

## Determinant equations

Use source coordinates \((z,q,x)\) on \(z\ne0\), where
\(y=(q-x^2)/z\) and the determinant of the coordinate change is \(-z\).
At \(z=1\) put \(r=q=x^2+y\), so \(y=r-x^2\). Let \(H_i\) be the
dehomogenization of \((P_{5-i},Q_{5-i},R_{4-i})\) in \((x,r)\).
In particular \(H_0=(r^2,r,r^2)\). If
\(D_i=\operatorname{diag}(5-i,5-i,4-i)\), the determinant component of
ordinary degree \(11-n\) at \(z=1\), up to the common sign, is
\[
         J_n=\sum_{i+j+k=n}
         \det\big(D_iH_i,\partial_rH_j,\partial_xH_k\big).
 \tag{1}
\]
Indeed the Euler column \(z\partial_z+2q\partial_q+x\partial_x\) in the
coordinates \((z,q,x)\) can replace the first column at \(z=1\), since
its additional two terms are multiples of the other columns.
The Keller equations require \(J_n=0\) for \(n<11\) and \(J_{11}\ne0\).
Terms with constant homogeneous components make no contribution.

The coefficient of the unknown \(H_n=(p,t,s)\) in \(J_n\) is
\[
          \partial_x\big(r^2(6p-2rt-5s)\big).
 \tag{2}
\]
The two highest rows used in this cofactor calculation are represented
by \(D_0H_0=(5r^2,5r,4r^2)\) and \(\partial_rH_0=(2r,1,2r)\).
Equation (2) also follows before dehomogenization from the cofactor column
\((6z^2q^2,-2q^3,-5z^3q^2)\) for the \(\partial_x|_{z,q}\) derivative.

## Normalizing the degree-four correction

Write the ring uniquely in the basis \(x^i y^j\), \(i=0,1\), \(j\ge0\),
over \(\mathbb C[q,z]\), by monic division using \(x^2=q-yz\).
This basis also works with \(z=1\) over \(\mathbb C[r]\).
At \(n=1\), equation (2) gives
\[
             6z^2P_4-2qQ_4-5z^3R_3\in\mathbb C[q,z].
\]
For each positive-degree basis monomial, let \(p,t,s\) denote its
coefficients in \(P_4,Q_4,R_3\). Then \(q\) divides \(6p-5zs\)
in \(\mathbb C[q,z]\). The coefficient triple therefore has the form
\[
       (q,3z^2,0)\frac{p-5zs/6}{q}+(z,0,6/5)\frac{5s}{6}.
\]
Consequently, up to terms in \(\mathbb C[q,z]\), the correction
\((P_4,Q_4,R_3)\) is \((q,3z^2,0)U+(z,0,6/5)V\) for homogeneous
\(U,V\) of degrees two and three.

At \(n=2\), reduce the determinant at \(z=1\) modulo \(r\). The term
(2) vanishes. With the actual first column \(\partial_z|_{q,x}\),
the two differentiated leading columns on \(r=0\) are
\(0\) and \((0,1,0)\), respectively. Thus the quadratic term reduces
to the determinant in the first and third entries of the
\(\partial_z|_{q,x}\) and \(\partial_x|_{z,q}\) columns of the
correction. Let \(v(x)=V(x,-x^2,1)\). If \(A\) and \(C\) are the
constant terms of the radial corrections to \(P_4,R_3\) on this
chart, the two entries before taking the derivatives are \(v+A\)
and \(6v/5+C\). The determinant is, up to sign,
\[
              \big(6v/5+24A/5-3C\big)v'.
\]
Here the \(x\partial_x\) Euler contributions to the first column
cancel in the determinant. This product vanishes, so \(v\) is
constant. Homogeneity gives \(V=qL\) modulo \(\mathbb C[q,z]\)
for a linear \(L\). The exact identity
\[
    (z,0,6/5)qL=\tfrac35(2zq,z^3,2q)L-\tfrac15(q,3z^2,0)zL
\]
places the first term in the direction of the \(q\)-derivative
of the leading triple. A source translation induces any
prescribed linear correction to \(q\) and can remove this term.
The correction from changing \(z\) in the leading forms depends
only on \(q,z\).

Absorb the part of \(U\) in \(\mathbb C[q,z]\) into the radial
terms. Subtract multiples of \(R\) from \(P,Q\) to eliminate
the \(q^2\) coefficients of \(P_4,Q_4\). At \(z=1\) the result is
\[
    H_1=(rU+a+br,\ 3U+d+er,\ g+hr),\qquad
    U=A_1xy+B_1y^2+\alpha x+\beta y .
 \tag{3}
\]
The exact degree-nine equation, upon integrating in \(x\) at
fixed \(r\), is
\[
 r^2(6p_3-2r t_3-5s_2+6U^2)
  +r(15g-24a)U+r^2(-12b+8d+5h+4er)U\in\mathbb C[r],
 \tag{4}
\]
where \(H_2=(p_3,t_3,s_2)\). The primitive of the quadratic
contribution from \(H_1\) in (4) is verified by differentiating it
and using (1). If \(15g-24a\ne0\), divide (4) by \(r\) and reduce
modulo \(r\). The positive-degree basis coefficients of \(U\)
are scalars at \(z=1\), and the same basis on \(r=0\) is
linearly independent. Thus \(U=0\) in that case. If \(U\ne0\)
then \(g=8a/5\), and the \(y^4\) basis coefficient in (4) divided
by \(r^2\) is \(6B_1^2\). Therefore \(B_1=0\).

## Exact lower coefficient witnesses

The coefficient witnesses are replayed in
[verify_quintic_rank_two_tangent_double_conic.py](verify_quintic_rank_two_tangent_double_conic.py).
Here is the rational elimination used by the replay. In step \(n\),
choose for each component of \(H_n\) the monomials
\(x^i(r-x^2)^j\), in increasing \(i\) and then \(j\),
with \(i+j\) at most its ordinary homogeneous degree before
dehomogenization. Integrate the known part of (1) in \(x\).
Discard the \(x\)-independent part of this primitive, and apply the
same projection to the three cofactor products in (2). Set the
resulting coefficients of \(x^i r^j\) to zero in descending
lexicographic order.

The last step produces a linear system for the unknown coefficients
of \(H_n\) with a rational matrix. Use the first independent columns
and their first independent rows to solve for a particular
solution. The remaining rows and monomials give necessary
compatibility equations. Add the ordered rational nullspace basis
with free parameters \(f_i\) for \(H_2\) and \(\ell_i\) for \(H_3\).
A compatibility polynomial displayed below is a coefficient
equation at such a step, after imposing the substitutions stated
before it. These are exact identities over \(\mathbb Q\).

### A nonzero correction \(U\)

Take \(g=8a/5\) and \(B_1=0\). The degree-nine system for \(H_2\)
has no compatibility conditions. The degree-eight system for \(H_3\)
has the witness \(10A_1^3\), so \(A_1=0\). Its witnesses include
\(12a\beta^2\), and when \(\beta=0\) they include \(12a\alpha^2\).
Thus \(a=g=0\) if \(U\ne0\).

If \(\beta\ne0\), the shear
\((x,y,z)\mapsto(x+\lambda z,y-2\lambda x-\lambda^2z,z)\)
preserves \(q,z\). It removes \(\alpha\) by taking
\(\lambda=\alpha/(2\beta)\), with any radial terms absorbed as in (3).
Rescaling \(y,z\) inversely and rescaling targets then makes
\(\beta=1\). If \(\beta=0\) and \(\alpha\ne0\), rescaling \(x,y\)
by \(\alpha,\alpha^2\) scales \(q\) by \(\alpha^2\).
After rescaling targets the correction has \(\alpha=1\).
Hence the remaining normalized choices at \(z=1\) are \(U=y\)
and \(U=x\).

For each choice, the degree-eight compatibility equations solve
linearly for \(f_6,f_8,f_9,f_{10},f_{11}\). The replay verifies a
square subsystem with nonzero constant determinant. Substituting
this forced solution into the degree-seven system for \(H_4\) gives
the compatibility coefficient \(35/3\) in both cases.

### The radial correction \(U=0\)

Use (3) with \(U=0\), and retain \(a,b,d,e,g,h\) without imposing
any relations. The degree-seven compatibility equations
successively contain
\[
   \frac5{12}f_8^2,\qquad \frac5{12}f_{10}^2,\qquad
   \frac23(f_3-f_5)^2,\qquad \frac5{12}(f_{11}-f_7)^2.
\]
At each step set the squared expression to zero before
checking the next. With these substitutions, the degree-eight
and degree-seven compatibility equations respectively contain
\[
     \frac12 f_9(8a-5g),\qquad
     \frac1{12}\big(5f_9^2-6(8a-5g)\ell_5\big).
\]
Multiplying the second equation by \(f_9\) and using the first
gives \(f_9=0\). Thus
\[
 H_2=(f_0+f_1r+f_4rx/3,\ f_2+f_4x+f_5r,\ f_6+f_7r).
\]
Recompute the parameterized \(H_3,H_4\) with this \(H_2\).
The degree-five determinant \(J_6\) has coefficients
\[
     [x^3]J_6=\frac53\ell_5^2,\qquad
     [x^2r]J_6=\frac{10}{9}f_4(3e\ell_5-f_4^2).
\]
The first gives \(\ell_5=0\) and the second gives \(f_4=0\).
For those two substitutions, formula (1) gives \(J_{11}=0\)
identically. This contradicts the Keller condition and finishes
the exclusion.
