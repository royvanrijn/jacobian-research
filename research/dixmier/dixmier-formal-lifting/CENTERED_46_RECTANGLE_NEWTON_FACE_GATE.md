# Necessary Newton faces for the centered order-(4,6) rectangle

This is an exact necessary-condition calculation and a method control,
not an exclusion of the order-(4,6) case or a claim that it lies beyond
published degree bounds. It proves no new DC1 counterexample or global
normalization theorem.

Let \(k\) be a characteristic-zero field, \([D,x]=1\), and let
\(P,Q\in A_1(k)\) satisfy \([Q,P]=\kappa\in k^\times\).
Write their Weyl symbols as
\[
 p=x^8\xi^4+\sum_{i<8,j<4}p_{ij}x^i\xi^j,\qquad
 q=x^{12}\xi^6+\sum_{i<12,j<6}q_{ij}x^i\xi^j.          \tag{1}
\]
All exponents are nonnegative. Thus both outer edges are single
centered monomials. The extra subprincipal gate is
\(p_{7,3}=0\). The
[mixed-corner commutator identity](SCALAR_GAUGE_SUBPRINCIPAL_AND_POLYNOMIAL_ADMISSION.md)
then gives \(q_{11,5}=0\). These are Weyl-symbol coefficients;
in normal PBW form the corresponding coefficients are 16 and 36.

## Literature context and the role of this rectangle

Guccione--Guccione--Valqui, [arXiv:1111.6100v3](https://arxiv.org/abs/1111.6100v3),
Corollary 7.4 on printed page 42, gives a lower bound \(B>15\),
where \(B\) is the minimum gcd of the **total degrees** over
counterexamples. The [cached source](sources/guccione-guccione-valqui-1111.6100v3.pdf)
and [selected inspection record](sources/guccione-guccione-valqui-1111.6100v3-source-audit.json)
fix that scope. The present corners have total degrees \(12,18\)
and gcd \(6\). They are a local method/control case, not a proposed
new counterexample beyond that published bound.

For general coprime \(n,m\), the gcd of the proportional corners' total
degrees is \(d+l\). In the one-root ratio-two family \(d=2l\),
this is \(3l\); exceeding 15 requires \(l\ge6\). Thus the first
larger tests are substantially higher differential-order cases.
The selected source statement and its final table were inspected;
its full earlier proof dependencies are not imported by this note.
The local necessary-condition proof below is independent.

## Polynomial coefficient vanishings

Independently of this extra subprincipal gate, (1) forces
\[
 \begin{split}
 p_{7,0}=p_{7,1}&=0,\\
 q_{11,0}=q_{11,1}=q_{11,2}=q_{11,3}
                    =q_{10,0}=q_{10,1}&=0.              \tag{2}
 \end{split}
\]
Use weights \(\operatorname{wt}(x,\xi)=(1,-r)\), successively
\(r=1/4,1/3,1/2\). An Ore contraction lowers weight by \(1-r>0\).
Thus the top nonzero part of the Weyl commutator is the Poisson
bracket of the top faces; higher Moyal terms have strictly lower
weight. All top bracket weights used below are positive, whereas
\(\kappa\) has weight zero. Consequently those Poisson brackets
must vanish. This argument uses the full finite rectangle support,
not a restriction of the unknown coefficients to selected monomials.

At \(r=1/4\), the potentially dominant terms of \(q\) are first
\(q_{11,0}x^{11}\), then \(q_{11,1}x^{11}\xi\). Their brackets
with the unavoidable corner of \(p\) have respective nonzero
coefficients \(-44q_{11,0}\) and \(-36q_{11,1}\). Both vanish.
The remaining top faces are
\[
 p_*=x^8\xi^4+A x^7,\qquad q_*=x^{12}\xi^6+B x^{11}\xi^2.
\]
Their bracket is
\[
 14(3A-2B)x^{18}\xi^5+14ABx^{17}\xi.
\]
It forces \(A=B=0\) over a field.

At \(r=1/3\), the remaining top faces are
\[
 p_*=x^8\xi^4+A x^7\xi,\qquad
 q_*=x^{12}\xi^6+B x^{11}\xi^3+C x^{10}.
\]
Their bracket has coefficients
\(10(3A-2B),10(AB-4C),-10AC\). They imply
\(A=B=C=0\). This proves all but the last vanishing in (2).
At \(r=1/2\), a remaining \(q_{10,1}x^{10}\xi\) would be
uniquely dominant. Its bracket has the nonzero coefficient
\(-32q_{10,1}\), proving the last vanishing.

## The first surviving face is a square and a cube

At \(r=1/2\), after (2), the full top-face spaces are
\[
 p_*=x^6F(u),\quad q_*=x^9G(u),\quad u=x\xi^2,
\]
where
\(F=u^2+A u+B\), \(G=u^3+C u^2+D u+E\).
The bracket is
\(6x^{15}\xi(2FG'-3GF')\). Its four coefficient equations are
\[
 3A-2C=0,\quad AC+6B-4D=0,\quad
 4BC-AD-6E=0,\quad 2BD-3AE=0.                           \tag{3}
\]
The first three determine \(C,D,E\); after substitution the last
is \(3(B-A^2/4)^2=0\). Over fields this yields a parameter
\(t=A/2\) such that
\[
          p_*=x^6(u+t)^2,\qquad q_*=x^9(u+t)^3.          \tag{4}
\]
No nonreduced parameter-scheme assertion is made.

## The next face under the zero-subprincipal gate

Set \(s=x^3(u+t)\), so (4) is \(p_*=s^2,q_*=s^3\).
Because \(p_{7,3}=q_{11,5}=0\), the next possible faces have
weights \(11/2,17/2\) and exactly the forms
\[
 p_1=Kx^6\xi,\qquad q_1=B x^{10}\xi^3+C x^9\xi.
\]
The next commutator equation has weight 14. Third Moyal contractions
of the top faces have weight \(27/2<14\), so it is precisely
\[
                   \{s^3,p_1\}+\{q_1,s^2\}=0.
\]
Equivalently \(\{s,3sp_1-2q_1\}=0\). Put
\(M=3K-2B\), \(N=3Kt-2C\). The latter bracket is
\[
 8Mx^{13}\xi^4+(14N-9Mt)x^{12}\xi^2-3Ntx^{11}.
\]
Thus \(M=N=0\), or
\[
                              q_1=\tfrac32s p_1.        \tag{5}
\]
The parameter \(K\) remains free at this step. These equations
give no contradiction and do not force a full operator square or
cube, an adjoint symmetry, or a finite mate.

[The exact checker](check_centered_46_newton_face.py) enumerates every
allowed rectangle monomial at each face and verifies all displayed
polynomial identities with independent Fraction arithmetic.
[The receipt](centered_46_newton_face_controls.json) pins the note
and checker. The result is scoped to (1), with the stated extra
subprincipal gate only for (5). No source theorem is imported as a
premise of this necessary-condition calculation.
