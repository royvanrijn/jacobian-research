# The quartic pencil of a double conic in the rank-two quintic case

**Theorem QuinticKellerRankTwoDoubleConicQuarticPencil.** Let \(q\) be a
geometrically irreducible ternary conic and \(z\) a source linear
coordinate over a characteristic-zero field. There is no three-variable
Keller map \(F=(P,Q,R)\) with homogeneous leading parts
\[
                     P_5=z^5,\qquad Q_5=zq^2,\qquad R_4=q^2 .
\]
Thus for the double conic the quartic-argument case \(m=4\) in the
[rank-two quartic pencil](QUINTIC_KELLER_RANK_TWO_QUARTIC_PENCIL.md)
is impossible. The proof works over \(\mathbb C\), and applies in
characteristic zero by embedding the finitely generated coefficient
field.

## The first mixed equations

Put \(S=Q-zR\), of degree at most four, and \(C=q^2\). Let
\(\{\cdot,\cdot\}\) denote the binary Jacobian bracket in \(x,y\).
Writing \(\kappa=\det JF\), the exact determinant identity is
\[
       \kappa=P_z\{S,R\}-S_z\{P,R\}+R_z\{P,S\}-R\{P,R\}.
 \tag{1}
\]
It follows by expanding \(dP\wedge d(S+zR)\wedge dR\).
Use \(r=q(x,y,1)\) for the dehomogenized conic, which has degree two
and is irreducible since \(q\) is irreducible and \(z\nmid q\).
The polynomial constants of the binary Jacobian derivation of \(r\)
are exactly \(\mathbb C[r]\), by the closed-polynomial and
differential criteria of Wang and Yan,
[Proposition 2.2 and Lemma 2.3, printed pages 5–6](https://arxiv.org/pdf/2610.36495v1).
In particular an ordinary homogeneous form whose dehomogenization
belongs to \(\mathbb C[r]\) belongs to \(\mathbb C[q,z]\).
We call such forms forms in \(q,z\).

In degree ten, (1) gives
\(\{5z^4 S_4-C P_4,C\}=0\). Dehomogenize at \(z=1\), denoting
the values of \(P_i,S_i,R_i\) by \(p_i,s_i,t_i\).
It follows that \(5s_4-r^2p_4\in\mathbb C[r]\).
Reduction modulo \(r^2\) and the degree bound four show
\(s_4\in\mathbb C[r]\), since the remainder of a polynomial in \(r\)
modulo \(r^2\) has degree at most two in \(x,y\).
Then \(p_4\) also belongs to \(\mathbb C[r]\).
Subtract suitable constant multiples of \(R\) from \(P,Q\)
to remove the \(q^2\) terms of \(P_4,S_4\).
They now have the form \(a z^2 q+c z^4,\ b z^2q+d z^4\).

If \(a=b=0\), the degree-nine equation gives
\(5s_3-r^2p_3\in\mathbb C[r]\), and the same degree argument
makes \(P_3,S_3\) forms in \(q,z\).
If they too depend only on \(z\), the degree-eight equation
makes \(P_2,S_2\) forms in \(q,z\).
If these again depend only on \(z\), the degree-seven
equation makes \(P_1,S_1\) univariate, and \(P\) is univariate
of degree five with a critical point. Thus take the first
index \(k\) among \(4,3,2\) with a nonzero coefficient of \(q\)
in at least one of \(P_k,S_k\), and write
\[
          P_k=a z^{k-2}q+\text{univariate term},\qquad
          S_k=b z^{k-2}q+\text{univariate term}.
\]
All higher terms of \(P,S\) are univariate.
The degree-\((k+5)\) equation of (1) gives
\[
    2r(5s_{k-1}-r^2 p_{k-1})-(5b-a r^2)t_3\in\mathbb C[r].
 \tag{2}
\]
Indeed its only possible nonzero contributions use \(P_{k-1},S_{k-1}\)
against \(C\) and \(P_k,S_k\) against \(R_3\), multiplied by
the highest \(P_z\) or \(R\) terms in (1). Both highest
quartic brackets with a form in \(q,z\) vanish.
Terms with \(S_z\) or \(R_z\) have lower degree.

If \(b\ne0\), reduce (2) modulo \(r\) to obtain \(t_3=rL+c_0\)
with \(L\) affine in \(x,y\). Subtract the part with \(c_0\),
divide by \(r\), and reduce modulo \(r^2\).
The degree bound \(\deg s_{k-1}\le3\) forces
\(s_{k-1}-(b/2)L\in\mathbb C[r]\), then
\(p_{k-1}-(a/2)L\in\mathbb C[r]\).
Lift \(L\) homogeneously to degree one. Replacing \(q\) by
\(q+L/2\) absorbs precisely the nonradial terms just found:
\(R_4+R_3\) is then a square plus a univariate cubic modulo
degree at most two. The nondegeneracy of the irreducible
ternary conic makes the linear correction removable by
translation, since the linear parts of its translations span
all linear forms. After that translation \(q+L/2\) becomes
\(q\) plus a constant. To preserve \(S=Q-zR\) as the translated
polynomial \(S\), subtract the translation multiple of \(R\)
from \(Q\). Thus \(R_3,P_{k-1},S_{k-1}\) are now forms in
\(q,z\), with the same highest \(a,b\) coefficients.

## The case \(b=0\)

Here \(a\ne0\). Dividing (2) by \(r\) gives
\[
                      10s_{k-1}+a r t_3-2r^2p_{k-1}\in\mathbb C[r].
\]
Reduce its polynomial in \(r\) modulo \(r^2\).
The degree of \(10s_{k-1}+a r t_3\) is at most five, so
the remaining quotient by \(r^2\) is affine in \(x,y\), say \(a L\).
Consequently \(t_3=rL+T\), where \(\deg(rT)\le k-1\).
Also \(p_{k-1}-(a/2)L\) is a polynomial in \(r\), and
\(s_{k-1}=-(a/10)rT+\text{a polynomial in }r\).
If \(k\le3\) then \(T\) is constant. Complete \(q+L/2\) and
translate as in the generic case to obtain the same forms
in \(q,z\).

For \(k=4\), that translation instead gives
\[
       R_3=z^2 U+f z^3,\qquad P_3=u zq+v z^3,\qquad
       S_3=-\frac{a}{10}qU+w zq+e z^3 ,
\]
with \(U\) binary linear. The degree-eight equation in (1),
at \(z=1\), gives
\[
     2r(5s_2-r^2 p_2)+a r^2t_2-\Theta(r)U+\frac{3a}{4}U^2
                    \ \in\ \mathbb C[r],
 \tag{3}
\]
where
\[
                         \Theta(r)=5w-4ad-af+
                                  (4ac/5-u)r^2+2a^2r^3/5 .
\]
For verification the terms involving \(U\) outside the linear
equation for \(p_2,s_2,t_2\) are multiples of \(\{r,U\}\).
Their coefficients come from \(5z^4\{S_3,R_3\}\),
\(P_{4,z}\{S_3,C\}\), \(-S_{4,z}\{P_4,R_3\}\),
\(-C\{P_3,R_3\}\), \(-R_3\{P_4,R_3\}\), and
\(C_z\{P_4,S_3\}\). Dehomogenized, the terms with \(q_z\)
cancel and the sum of the coefficients is
\(\Theta(r)-(3a/2)U\). This is the bracket of the last
two terms in (3) with \(r\).

Reducing (3) modulo \(r\) makes
\(-\Theta(0)U+(3a/4)U^2\) constant on the irreducible conic
\(r=0\). Since \(a\ne0\), factor the resulting univariate
quadratic over \(\mathbb C\); in the domain \(\mathbb C[x,y]/(r)\)
it forces \(U\) to be constant. A binary linear polynomial
can restrict to a constant on this conic only if it is
identically zero. Thus \(U=0\), completing the reduction
to \(R_3,P_{k-1},S_{k-1}\) forms in \(q,z\) in every case.

## Finishing the descent

Take \(j=2,1\) in order and write \(\ell=k+j-4\).
The next relevant degree \(k+j+2\) of (1) has the same
binary equation as (2), now with \(p_\ell,s_\ell,t_j\):
\[
        2r(5s_\ell-r^2p_\ell)-(5b-a r^2)t_j\in\mathbb C[r].
 \tag{4}
\]
Terms \(p_\ell,s_\ell\) with \(\ell\le0\) may be omitted from
this equation since their brackets vanish. The degree assertion
follows as in (2). In particular \(P_i,S_i\) above \(\ell\)
and \(R_i\) above \(j\) are already forms in \(q,z\), and
\(P_i,S_i\) above \(k\) depend only on \(z\).
The term with \(R_z\) in (1) has degree at most
\(k+\ell+1< k+j+2\); the term with \(S_z\) also has lower degree.

If \(\ell\le0\), equation (4) immediately gives \(t_j\in\mathbb C[r]\)
since \(5b-a r^2\ne0\) as a polynomial. If \(\ell\ge1\) and
\(b\ne0\), first reduce modulo \(r\). With \(\deg t_j\le2\) this
already gives \(t_j\in\mathbb C[r]\), and reducing the remaining
relation modulo \(r^2\) gives \(s_\ell,p_\ell\in\mathbb C[r]\).
Finally if \(\ell\ge1\) and \(b=0\), divide by \(r\) as in (3)
with \(U=0\). Since \(\deg t_j,\deg s_\ell\le2\), the quotient
of \(10s_\ell+a r t_j\) by \(r^2\) after subtracting a remainder
in \(\mathbb C[r]\) is constant. Then \(t_j\) is a linear
polynomial in \(r\) since its remaining nonconstant part,
multiplied by \(r\), would have degree greater than two.
The same equation puts both \(s_\ell,p_\ell\) in \(\mathbb C[r]\).

Thus all terms of \(P,S,R\) are forms in \(q,z\). The same
is true of \(Q=S+zR\), so the full map factors through the two
polynomials \(q,z\). Its determinant is zero, a contradiction.
