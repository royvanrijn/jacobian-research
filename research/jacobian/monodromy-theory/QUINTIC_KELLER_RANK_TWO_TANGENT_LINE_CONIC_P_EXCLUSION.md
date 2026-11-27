# Exclusion of the pure-fifth span in the tangent double-line case

**Theorem QuinticRankTwoTangentLineConicPureFifthExclusion.**
There is no complex three-variable Keller map of degree five
with \(R=z+z^2(q+c)\), where \(q=x^2+zy\) and \(c\) is constant,
whose quintic leading span is \(\langle z^5,z^3q\rangle\).
Consequently the tangent irreducible conic branch of the
degree-two primitive pencil is excluded.

Set \(p=z^2\). Normalize the quintic parts to
\(P_5=z^5/5,\ Q_5=z^3q\). The
[pure-fifth restriction and its exceptional-plane corollary](QUINTIC_KELLER_RANK_TWO_TANGENT_LINE_CONIC_P_SPAN.md)
give \(g=0\) in the stated normal form. We write its full
quartic and cubic terms as
\[
\begin{aligned}
 P_4&=p z L+a p^2+b p q, & Q_4&=q z L+e p^2+f p q,\\
 P_3&=pV+2zL^2+4a p L+i z^3+kzq, &
 Q_3&=4e p L+qV+j z^3+n zq ,
\end{aligned}
\]
where \(L=A x+B y\) is nonzero and \(V=C x+D y\).
Any pure \(z\) term in the linear form \(V\) is absorbed into
the displayed invariant cubic terms. The quadratic and linear
terms of \(P,Q\) are arbitrary at this stage. Write
\(u_1=[xy]P_2\).

Let \(\Delta=\det J(P,Q,R)\).
The [exact coefficient calculation](verify_quintic_tangent_line_conic_p_exclusion.py)
gives
\[
                         [x^3 y^2 z^3]\Delta=-12 B^3 .
\]
All nonconstant coefficients must vanish. Thus \(B=0\), so
\(A\ne0\). The exceptional-plane corollary gives \(D=0\)
if \(V\ne0\), and this also holds if \(V=0\).
With these two substitutions the same exact calculation gives
\[
 \begin{aligned}
 [x^4 z^4]\Delta &= 6A^3-Ak+2u_1,\\
 [y^2 z^6]\Delta &= -Ak-u_1,\\
 [x^5 z^2]\Delta &= -4A u_1 .
 \end{aligned}
\]
The first two equations force \(k=2A^2\) and \(u_1=-2A^3\).
The third coefficient is then \(8A^4\ne0\), a contradiction.
The calculation is an identity with every lower coefficient
independent, so no specialization of the surviving parameters
evades these equations.

The
[two-span reduction](QUINTIC_KELLER_RANK_TWO_TANGENT_LINE_CONIC_SPANS.md)
and the
[divisible-span exclusion](QUINTIC_KELLER_RANK_TWO_TANGENT_LINE_CONIC_Q_SPAN.md)
had left only the span treated here. Together the results
exclude the tangent irreducible conic branch.
All conclusions hold in characteristic zero by embedding the
finitely generated coefficient field.
