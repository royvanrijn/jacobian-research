# The tangent pure-fifth double-conic pencil

**Theorem QuinticKellerRankTwoTangentPureFifthDoubleConic.** Over
characteristic zero there is no three-variable Keller map with
homogeneous leading parts
\[
        P_5=z^5,\qquad Q_5=z^3q,\qquad R_4=q^2,
        \qquad q=x^2+yz .
\]
Equivalently, the projective conic of the irreducible quadratic \(q\)
can be any conic tangent to the repeated pencil line \(z=0\). Over
\(\mathbb C\), its binary restriction can be made \(x^2\), and the
remaining terms are \(z\) times a linear form with a nonzero \(y\)
coefficient. Make that form the new \(y\). It suffices to prove
the result over \(\mathbb C\) by embedding the finitely generated
coefficient field.

This excludes the tangent case of the
[pure-fifth double-conic residual](QUINTIC_KELLER_RANK_TWO_DOUBLE_CONIC_PENCIL_RESIDUAL.md#removing-the-nonzero-coefficient).
The normalization of the conic and use of the polynomial invariant ring
parallel the tangent case in Wang and Yan's quartic proof. The cofactor
and lower equations here use degrees \((5,5,4)\).
[Wang and Yan, §§3.4.1 and 3.4.3](https://arxiv.org/pdf/2610.36495v1).

## Upper corrections

On \(z=1\) write \(r=q=x^2+y\), so the dehomogenized polynomials
are written in \((x,r)\) with \(y=r-x^2\). For each drop \(n\)
let \(H_n\) be the dehomogenization of
\((P_{5-n},Q_{5-n},R_{4-n})\), and put
\(D_n=\operatorname{diag}(5-n,5-n,4-n)\).
As in the [mixed tangent calculation](QUINTIC_KELLER_RANK_TWO_TANGENT_DOUBLE_CONIC.md#determinant-equations),
the degree-\((11-n)\) Jacobian at \(z=1\) is, up to the overall sign,
\[
   J_n=\sum_{i+j+k=n}
          \det(D_iH_i,\partial_rH_j,\partial_xH_k).
 \tag{1}
\]
Here the determinant of \((z,q,x)\) relative to the original
coordinates is \(-z\). The Euler column at \(z=1\) gives
the columns \(D_iH_i\). Thus \(J_n=0\) for \(n<11\) and \(J_{11}\)
must be nonzero. All constant homogeneous components have zero
derivatives.

Since \(H_0=(1,r,r^2)\), the cross product of \(D_0H_0\) and
\(\partial_rH_0\) is \((6r^2,-10r,5)\). For \(H_n=(p,t,s)\),
the part linear in this new correction in (1) is consequently
\[
                 \partial_x(6r^2p-10rt+5s).
 \tag{2}
\]
In particular the degree-ten equation is
\(6r^2p-10rt+5s\in\mathbb C[r]\) for \(H_1=(p,t,s)\).

Use the free basis \(x^i y^j\), \(i=0,1\), \(j\ge0\), over
\(\mathbb C[q,z]\) from monic division by \(x^2=q-yz\).
The same basis works over \(\mathbb C[r]\) at \(z=1\).
The degree-ten equation modulo \(r\) forces every positive-degree
basis coefficient of \(s\) to be divisible by \(r\).
Homogenization gives \(R_3=qL\) modulo \(\mathbb C[q,z]\) for
a linear \(L\). A source translation gives any desired linear
correction to \(q\) and removes \(qL\) from \(R_3\) through
\(R_4=q^2\). With \(s\in\mathbb C[r]\) the equation says
\(t=3rp/5\) modulo \(\mathbb C[r]\). For the positive-degree
basis coefficients, multiplication by \(q\) adds two to the
degree in the homogeneous coordinates. Because \(Q_4\) has
degree four, the correction to \(P_4\) has the form of a
homogeneous quadratic multiplied by \(z^2\), modulo
\(\mathbb C[q,z]\).

Subtract multiples of \(R\) from \(P,Q\) to remove their radial
\(q^2\) terms in degree four. We obtain
\[
        H_1=(U+a+br,\ 3rU/5+d+er,\ g+hr),
        \qquad U=Axy+By^2+\alpha x+\beta y .
 \tag{3}
\]
Here \(a,b,d,e,g,h\) are free scalars. The terms of \(U\)
depending only on \(r\) have been absorbed into them.

## Exact coefficient elimination

The witnesses below are verified in
[verify_quintic_rank_two_tangent_pure_fifth.py](verify_quintic_rank_two_tangent_pure_fifth.py).
The replay shares the determinant and monomial routines of the
[mixed tangent checker](verify_quintic_rank_two_tangent_double_conic.py).
For each unknown \(H_n\), use the basis \(x^i(r-x^2)^j\)
with \(i+j\) bounded by the respective homogeneous degree,
in increasing \(i\) and then \(j\). Integrate the known terms of
(1) in \(x\) at fixed \(r\), remove the part of the primitive
depending only on \(r\), and set its sum with the similarly
projected expression (2) to zero. Order the rows by descending
lexicographic order on the \(x,r\) exponents. The resulting
coefficient equations are linear in \(H_n\) with a rational
matrix. Use the first independent columns and their first
independent rows to parametrize all possible \(H_n\), by adding
the ordered nullspace basis to the particular solution. The
remaining rows and monomials give necessary compatibility
conditions. The free parameters in \(H_2\) are denoted \(f_i\).

For (3), the degree-nine system for \(H_2\) contains
the compatibility coefficient \(24B^2/5\). With \(B=0\),
it contains \(12A^2/5\), so \(U=\alpha x+\beta y\).
If \(\beta\ne0\), the shear
\((x,y,z)\mapsto(x+\lambda z,y-2\lambda x-\lambda^2 z,z)\)
preserves \(q,z\) and eliminates the \(\alpha x\) term modulo
radial terms. A uniform source dilation and componentwise
target rescaling restoring the three leading forms then normalize the
coefficient of \(y\) to one. The same dilation normalizes the
coefficient of \(x\) if \(\beta=0\) and \(\alpha\ne0\).
It remains to consider \(U=y,\ x,\ 0\).

For \(U=y\), the degree-eight system for \(H_3\) contains the
nonzero compatibility coefficient \(-2/25\).

For \(U=x\), the degree-seven system for \(H_4\) gives the
following three compatibility equations after removing nonzero
rational constants from their coefficients:
\[
      b f_4=0,\qquad
      625f_4^2+75f_4+3=0,\qquad
      375f_1f_4+1250f_4^2+3=0.
\]
The first two imply \(b=0\) since \(f_4\ne0\). The degree-eight
compatibility equation
\[
                   20be-50bf_5-15f_1+4=0
\]
then gives \(f_1=4/15\). Subtract twice the first quadratic from
the second after this substitution to obtain \(-50f_4-3=0\).
But the first quadratic has value \(3/4\) at \(f_4=-3/50\).

For \(U=0\), the degree-seven system gives the coefficient
\(-10 f_4^2/3\). With \(f_4=0\), recompute \(H_3,H_4\) using
the same linear procedure. Equation (1) then gives the exact
coefficient
\[
                         [x^2r^2]J_6=10 f_5^3/9 ,
\]
where \(J_6\) is the degree-five equation. With \(f_5=0\) as well,
all five pieces \(H_0,\ldots,H_4\) lie in \(\mathbb C[r]^3\).
Their \(\partial_x\) columns vanish, so (1) gives \(J_{11}=0\).
This contradicts the Keller condition.
