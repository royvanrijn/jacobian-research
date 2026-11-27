# Closure of a primitive cubic pencil containing the quartic member

**Proposition QuinticRankTwoClosedCubicPencil.** Let \(F=(P,Q,R)\)
be a three-variable complex Keller map of degree five whose leading
forms in some linear coordinates satisfy
\[
       P_5=z^5,\qquad Q_5=z^2 T,\qquad R_4=zT,
\]
where \(T\) is a homogeneous cubic coprime to \(z\) and
\(r=T(x,y,1)\) is a closed polynomial. Then \(F\) is invertible.
In particular, the \(m=3\) branch of the
[rank-two repeated-line reduction](QUINTIC_KELLER_RANK_TWO_REPEATED_LINE_PENCIL.md)
can survive only with \(R_4=z^4\).

The hypothesis on \(r\) means that \(\mathbb C(r)\) is relatively
algebraically closed in \(\mathbb C(x,y)\). Hence the polynomial
constants of the binary derivation \(\{\cdot,r\}\) form
\(\mathbb C[r]\): in characteristic zero a rational function with
zero bracket is algebraic over \(\mathbb C(r)\), and a reduced
fraction in \(r\) that is polynomial must have constant denominator
by a univariate Bézout identity after substitution. Also
\(\deg r=3\) because \(z\nmid T\).

In the
[rank-two cubic calculations](QUINTIC_KELLER_RANK_TWO_QUARTIC_PENCIL.md#further-restrictions-in-the-cubic-pencil-case),
take \(S=Q-zR\).
The determinant identity (3) and its leading equation (4) hold
under the displayed leading-form assumptions. The coefficient
calculations (5), the degree-nine equation, and the weighted
determinant step use only \(\deg r=3\), \(T(x,y,0)\ne0\), and the
description of the polynomial constants. They give a binary linear
form \(L\), constants \(\alpha,\beta,\gamma,\delta,\epsilon\),
and a binary linear \(L_1\) such that after removing \(\alpha R\)
from \(P\),
\[
  S_4=TL+\delta z^4,\qquad P_4=5z^3L+\beta z^4,\qquad
  P_3=\gamma T+10zL^2+z^2L_1+\epsilon z^3.
\]
Furthermore \(L=0\) or \(T(x,y,0)\) is proportional to \(L^3\).

For \(L=0\), use the
[descent for the vanishing binary linear term](QUINTIC_KELLER_RANK_TWO_QUARTIC_PENCIL.md#vanishing-binary-linear-term).
If \(P_3\) is univariate in \(z\), the bounded coordinate completion
there gives invertibility. Otherwise its final divisibility argument
gives \(T_y=0\) after a binary linear source change. The divisibility
uses only that \(T+\xi z^3\), for a constant \(\xi\), is coprime to
\(z\) and has degree three. Now \(r=T(x,1)\) would be a univariate
polynomial of degree three. This contradicts closedness because the
remaining binary coordinate \(x\) is algebraic of degree three over
\(\mathbb C(r)\).

For \(L\ne0,\gamma=0\), the
[zero-\(\gamma\) proof](QUINTIC_KELLER_RANK_TWO_CUBIC_PENCIL_ZERO_GAMMA.md)
depends only on the resulting normal form of \(P\) and the Keller
condition. It gives bounded coordinate completion or fibrewise
planar inversion. For \(L\ne0,\gamma\ne0\), take \(x=L\) and use
the
[critical-point and residue proof](QUINTIC_KELLER_RANK_TWO_CUBIC_PENCIL_CLOSURE.md).
It begins with its formula (1), which follows here from the same
\(P_4,P_3\) and \(T(x,y,0)\) above, and uses only that formula and
the Keller condition to prove invertibility. Thus all the
calculations after the normal form apply to closed reducible cubic
members as well.

In the \(m=3,\ R_4=zJ\) branch of the repeated-line theorem, take
\(T=J\) and choose the target quintic span basis \(z^5,z^2J\).
The pencil has basis \(z^3,J\), so \(J(x,y,1)\) generates the closed
pencil field, and the proposition applies. The same result holds
over characteristic zero by the coefficient-field embedding and
inverse descent.
