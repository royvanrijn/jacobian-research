# Exclusion of the quintic span divisible by the tangent conic

**Proposition QuinticRankTwoTangentLineConicDivisibleSpanExclusion.**
Let \(F=(P,Q,R)\) be a complex three-variable Keller map of ordinary
degree five with \(R=z+z^2(q+c)\), where \(q=x^2+zy\) and \(c\)
is constant. Its quintic leading span cannot be
\(\langle z^3q,zq^2\rangle\).

Use \(p=z^2\), and the derivations \(E,E_*,\delta\) and
notation of the
[two-span proof](QUINTIC_KELLER_RANK_TWO_TANGENT_LINE_CONIC_SPANS.md).
We may normalize \(P_5=z^3q,\ Q_5=-zq^2/3\) so that
\(E(P_5)=zqp,\ E(Q_5)=zq^2\). The common factor in that proof
is \(H=q\). For the first noninvariant layer \(j=3,4\), we
have \(P_j=pW+I_j,\ Q_j=qW+J_j\), where \(W\) is homogeneous
of degree \(j-2\) and noninvariant, and \(I_j,J_j\in\mathbb C[z,q]\).
The class of \(W\) modulo \(q\) remains noninvariant by
localization, since a homogeneous linear or quadratic polynomial
invariant on \(q=0\) is by homogeneity a multiple of \(z^{j-2}\)
there and differs from that multiple by an invariant multiple of
\(q\).

## First noninvariance in degree four

Set \(u=x/z\) on \(q=0,\ z\ne0\), and denote by a bar reduction
on that surface. Equation (4) of the two-span proof implies that
\(J_4\) is divisible by \(q\). Indeed, modulo \(q\) its
only potentially nonconstant term in \(u\) is
\(-pE(J_4)W\), and \(E(J_4)\) modulo \(q\) is four times the
\(z^4\) term of \(J_4\). Thus \(Q_4\) is divisible by \(q\).
After dividing that equation by \(q\), reduce once more on
\(q=0\). The result and (3) of the two-span proof give
\[
                 \overline P_4=z^4 A(u),\qquad
                 \overline Q_3=z^3(-2A(u)^2+d),
\]
where \(A\) is nonconstant and \(d\) is a scalar.
In detail, the \(W^2\) term after division is \(2p\overline W^2\),
and the linear term is \(4a p^2\overline W\) with \(a z^4\)
the pure \(z\) term of \(I_4\). Since \(J_4\) is a combination
of \(pq,q^2\), \(E(J_4)/q\) vanishes modulo \(q\).

In degree seven in the bracket \(E_*P\,\delta Q-\delta P\,E_*Q\),
the only surviving term modulo \(q\) comes from \(P_4,Q_3\).
The sources involving \(E(P_5),E(Q_5)\) or \(Q_4\) vanish there.
The contribution from \(E_*-E\) that lowers degree two involves
a quintic and a quartic. It vanishes there since \(\delta Q_4\)
and \((Q_5)_y\) are both divisible by \(q\). The contribution
that lowers degree three uses only quintic invariants.
On this surface \(E=z\partial_z-2u\partial_u\).
The remaining term in degree seven is
\(z^7(-10A^2-3d)A'\), which cannot vanish for nonconstant \(A\).

## First noninvariance in degree three

The quartic terms \(P_4,Q_4\) are invariant. The degree-seven
bracket equation is
\[
 zq(p\delta Q_2-q\delta P_2)+
                  (qE(P_4)-pE(Q_4))\delta W=0.
\tag{1}
\]
Reducing modulo \(q\) shows \(q\mid Q_4\). Write
\[
 P_4=a p^2+b pq+d_4 q^2,\qquad Q_4=e pq+f q^2 .
\]
After canceling \(q\), (1) becomes
\[
 z(p\delta Q_2-q\delta P_2)
                       +4(a p^2+f p q-d_4 q^2)\delta W=0.
\tag{2}
\]
As in the two-span proof, the \(u\)-degrees of \(P_2,Q_2\)
are at most \(\deg_u W\). Otherwise their leading
coefficients would have ratio \(p/q\), but the \(u\)-leading
coefficient of a homogeneous quadratic with positive \(u\)-degree
is a scalar multiple of \(z^2\). In particular \(P_2,Q_2\)
have \(u\)-degree at most two, so each is a linear combination
of \(p,q,zx,zy\). The coefficient of \(q^2\) in (2), in
coordinates \(u,z,q\), forces \(d_4=0\).
The remaining terms give
\[
                 P_2=4fzW+I_2,\qquad Q_2=-4azW+J_2,
                       \qquad I_2,J_2\in\mathbb C[z,q].
\tag{3}
\]

In degree six in the bracket, the term from \(P_3,Q_3\)
integrates under \(\delta\) to \(2pqW^2\) plus an invariant
multiple of \(W\). The other terms after integration are
invariant multiples of \(P_1,Q_1,P_2,Q_2,W\).
The contribution involving \(E_*-E\) uses quintics against
cubics and is also linear in \(W\), since the quintics are
invariant. If \(W\) had a nonzero \(y\)-coefficient, \(W^2\)
would have \(u\)-degree four, while every other term would
have \(u\)-degree at most two. Thus, after absorption of a
pure \(z\) term into \(I_3,J_3\), \(W=A x\) for a
nonzero scalar \(A\).
The only quadratic term in \(u\) from the cubic pair is now
\(2pq A^2 x^2\). Replacing \(x^2\) by \(q-zy\) and comparing
the \(y\)-coefficient with the quintic and linear contribution
\(zq(pQ_1-qP_1)\) gives \((P_1)_y=0,\ (Q_1)_y=2A^2\).

Write \(J_3=j_0 z^3+n zq\). In degree six modulo \(q\),
the terms from the quartic and quadratic layers and from
the cubic pair give
\[
                   (-16a^2-3j_0)z^5\,\delta W=0.
\]
All other terms vanish modulo \(q\). In particular, the
degree-two correction involving \(E_*-E\) vanishes because
\(\delta Q_3\) and \((Q_5)_y\) are divisible by \(q\).
This gives \(j_0=-16a^2/3\).

A constant rescaling preserves the scalar shape and normalizes
\(A=1\). Namely replace \(z,x,y\) by
\(\lambda z',\lambda^{-1/2}x',\lambda^{-2}y'\) with
\(\lambda^{1/2}=A\), and rescale \(R,P,Q\) by
\(\lambda^{-1},\lambda^{-2},\lambda\), respectively.
The full forms, up to constant terms within \(P,Q\), are now
\[
\begin{aligned}
 P={}&z^3q+a p^2+b pq
      +p x+i z^3+k zq+4fzx+g p+h q+d z+r x,\\
 Q={}&-zq^2/3+e pq+f q^2
       +q x-\tfrac{16}{3}a^2 z^3+n zq
       -4azx+\ell p+m q+2y+v x+w z .
\end{aligned}
\]
In their exact Jacobian determinant with \(R\), the coefficient
of \(x^2z^3\) minus the coefficient of \(yz^4\) is \(-10\).
The [exact symbolic checker](verify_quintic_tangent_line_conic_q_span.py)
verifies this identity with all the displayed constants arbitrary.
It contradicts the constant determinant condition and eliminates
\(j=3\). The result extends to characteristic zero by embedding
a finitely generated field of coefficients.
