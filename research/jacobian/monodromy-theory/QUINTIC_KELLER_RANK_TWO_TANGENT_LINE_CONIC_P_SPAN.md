# The pure-fifth span in the tangent double-line case

**Proposition QuinticRankTwoTangentLineConicPureFifthRestriction.**
Let \(F=(P,Q,R)\) be a complex three-variable Keller map of
degree five with \(R=z+z^2(q+c)\), \(q=x^2+zy\), \(p=z^2\),
and \(P_5=z^5/5,\ Q_5=z^3q\). Then for some binary linear
form \(L(x,y)\ne0\), a homogeneous linear form \(V(x,y,z)\),
scalars \(a,b,e,f,g\), and homogeneous cubics
\(I_3,J_3\in\mathbb C[z,q]\), the lower terms satisfy
\[
\begin{aligned}
 P_4&=p z L+a p^2+b p q,&
 Q_4&=q z L+e p^2+f p q+g q^2,\\
 P_3&=p V+2zL^2+(4a p+4g q)L+I_3,&
 Q_3&=4e pL+q V+J_3 .
\end{aligned}
\]
In particular the first noninvariant homogeneous layer for
\(\delta=z\partial_x-2x\partial_y\) is quartic.

Use \(E,E_*,\delta\) from the
[two-span proof](QUINTIC_KELLER_RANK_TWO_TANGENT_LINE_CONIC_SPANS.md).
The common factor of the corresponding binary quadratics is
\(H=p\), since \(E(P_5)=z p^2,\ E(Q_5)=z p q\).
That proof gives \(j=3\) or \(4\) for the first noninvariant
layer, and
\[
                 P_j=pW+I_j,\qquad Q_j=qW+J_j,
\tag{1}
\]
where \(W\) is noninvariant of degree \(j-2\) and
\(I_j,J_j\) are invariant. Write the invariant quartic
terms in either case as
\[
                    a p^2+b p q+d_4 q^2,\qquad
                    e p^2+f p q+g q^2 .
\]

## Excluding the cubic first layer

If \(j=3\), the degree-seven equation in the bracket
\(E_*P\,\delta Q-\delta P\,E_*Q\) is
\[
 z p(p\delta Q_2-q\delta P_2)
           +4(a p^2 q-d_4 q^3-e p^3+g p q^2)\delta W=0 .
\tag{2}
\]
Here \(\delta W\) is a nonzero linear combination of \(z,x\),
so it is divisible by at most one power of \(z\).
All terms apart from the \(d_4\)-term are divisible by
\(z^2\). Thus \(d_4=0\).
Cancel \(p\) from (2). After inverting \(z\), the \(u=x/z\)
degree of \(P_2,Q_2\) cannot exceed that of \(W\);
otherwise the highest coefficients would have ratio \(p/q\),
while every leading coefficient of a homogeneous quadratic
of positive \(u\)-degree is a scalar multiple of \(z^2\),
as proved in the two-span argument. In particular
\(P_2,Q_2\) are combinations of \(p,q,zx,zy\).
The \(q^2\) coefficient after cancellation in (2) forces
\(g=0\), since \(\delta P_2,\delta Q_2,\delta W\) are independent
of \(q\) in the localized coordinates \(z,q,u\).
The other coefficients give
\[
           P_2=4azW+I_2,\qquad Q_2=4ezW+J_2
\]
for invariant \(I_2,J_2\).

In degree six, integrate the bracket equation under \(\delta\).
The contribution from the two cubic terms is \(2pqW^2\)
plus an invariant multiple of \(W\). The other terms are
invariant multiples of \(P_1,Q_1,P_2,Q_2,W\). The corrections
from \(E_*-E\) involving quintics against cubics are still
linear in \(W\). If \(W\) had a nonzero \(y\)-coefficient,
its square would have \(u\)-degree four, while the other
terms have degree at most two. Absorb a pure \(z\) term
into the invariants in (1) and write \(W=A x\), \(A\ne0\).
The only quadratic term in \(u\) from the cubic pair is
\(2pq A^2 x^2\). Substitute \(x^2=q-zy\) and compare with
the contribution \(zp(pQ_1-qP_1)\). This yields
\[
                        (P_1)_y=-2A^2,\qquad (Q_1)_y=0 .
\]

On \(z=0\) the scalar \(R\) has \(R_x=R_y=0,\ R_z=1\).
The displayed forms give \(P_y=-2A^2,\ Q_y=0\) there, since
the quartic invariants have \(d_4=g=0\) and the quadratic
terms restrict to multiples of \(q|_{z=0}=x^2\).
The restrictions of the cubics there are \(0\) and \(A x^3\),
respectively. Hence the determinant on that plane is
\(2A^2 Q_x\), whose \(x^2\) coefficient is \(6A^3\ne0\).
This rules out \(j=3\).

## The necessary quartic and cubic terms

With \(j=4\), equation (4) of the two-span proof reads
\[
 z p(pQ_3-qP_3)+2pqW^2+
          4(a p^2 q-d_4 q^3-e p^3+g p q^2)W
                                                \in\mathbb C[z,q].
\tag{3}
\]
We can choose \(W\) with only \(zx,zy,xy,y^2\) terms,
absorbing the invariants \(p,q\) into \(I_4,J_4\).
If its restriction to \(z=0\) is nonzero, it depends on \(y\).
The restriction of (3) to that plane therefore gives
\(d_4=0\). Otherwise \(W=zL\) for a nonzero binary linear
\(L\). The invariant on the right of (3) vanishes at \(z=0\).
It is homogeneous of degree eight, so all powers of \(z\)
in its expression in \(z,q\) are even, and it is divisible
by \(z^2\). Dividing (3) by \(z\) and restricting to the
plane again gives \(d_4=0\).

In either case we can now divide (3) by \(p\) in the invariant
ring. Its restriction to \(z=0\) becomes
\[
                2x^2 W(x,y,0)^2+4g x^4 W(x,y,0),
\]
which can be independent of \(y\) only if \(W(x,y,0)=0\)
with our chosen form of \(W\). Thus \(W=zL\) in all cases.
Divide the equation once more by \(z\) to obtain
\[
 pQ_3-qP_3+2zq L^2+(4a p q-4e p^2+4g q^2)L
                                               \in\mathbb C[z,q].
\tag{4}
\]
Modulo \(q\), after inverting \(z\), the class of
\(Q_3-4e pL\) is invariant and homogeneous of degree three.
It is therefore a multiple of \(z^3\), and we may write
\(Q_3=4e pL+qV+J_3\) for a homogeneous linear \(V\)
and invariant \(J_3\). Substituting this into (4) then
gives \(P_3=pV+2zL^2+(4a p+4g q)L+I_3\). This proves
the stated restriction. The result transfers to characteristic
zero by embedding the finitely generated coefficient field.

## Restriction on the exceptional plane

**Corollary QuinticRankTwoTangentPureFifthExceptionalPlane.**
In the normal form of the proposition, \(g=0\).
If \(V(x,y,0)\ne0\), both \(V(x,y,0)\) and \(L(x,y)\)
are multiples of \(x\), and \(P(x,y,0)\) is affine linear
with linear part proportional to \(x\).
If \(V(x,y,0)=0\) and \(P(x,y,0)\) is quadratic, then \(L\)
is proportional to \(x\) and the quadratic leading term
of \(P(x,y,0)\) is proportional to \(x^2\).

Let \(P_0^*=P(x,y,0)\) and \(Q_0^*=Q(x,y,0)\).
Since \(R_x=R_y=0,\ R_z=1\) there, the pair
\(P_0^*,Q_0^*\) has constant nonzero planar Jacobian.
It is a planar automorphism by the
[retained degree-100 planar theorem](../../BREAKTHROUGH_TARGET_ASSESSMENT.md#four-variable-hessian-conjecture-and-planar-degree-bounds),
as both degrees are at most five.
If \(g\ne0\), the leading forms on that plane are
\(4g x^2 L\) for \(P_0^*\) and \(g x^4\) for \(Q_0^*\).
Degrees three and four contradict the
[plane automorphism division lemma](https://arxiv.org/pdf/math/0408077v1).
Thus \(g=0\).

In this frame \(P_0^*\) has degree at most two, and the
cubic part of \(Q_0^*\) is \(x^2 V(x,y,0)\).
If \(V(x,y,0)\ne0\), the division lemma forces \(P_0^*\)
to be linear up to its constant term. The highest planar
Jacobian equation makes \(x^2 V(x,y,0)\) a cube of that
linear part, so \(V(x,y,0)\) and the linear part of \(P_0^*\)
are proportional to \(x\).

In the full determinant, \(R_x,R_y\) are divisible by
\(z^2\) and \(R_z=1+2z(x^2+c)+O(z^2)\). Hence the
coefficient of \(z\) gives the binary bracket identity
\[
       \{P_z|_{z=0}, Q_0^*\}
       +\{P_0^*, Q_z|_{z=0}\}+2(x^2+c)\kappa=0 .
\]
The cubic part of \(Q_z|_{z=0}\) is \(x^2L\), and the
quadratic part of \(P_z|_{z=0}\) is \(2L^2+\xi x^2\)
for a scalar \(\xi\).
When \(V(x,y,0)\ne0\), the only degree-three term of this
identity comes from the bracket of the latter quadratic
with a nonzero multiple of \(x^3\) in \(Q_0^*\).
Its vanishing gives \(L_y=0\).
When \(V(x,y,0)=0\) and \(P_0^*\) is quadratic,
the only degree-three term is the bracket of the leading
quadratic of \(P_0^*\) with \(x^2L\).
Two nonzero homogeneous forms of degrees two and three with
zero bracket are powers of a common linear form, by Euler's
identity. Thus \(x^2L\) is proportional to \(x^3\), which
gives the remaining assertions.
