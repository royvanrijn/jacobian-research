# Actual alternating monodromy on affine three-space

**KMFI12--13.** Over an algebraically closed characteristic-zero field,
the natural actions of A_(4j+1), j>=1, and A_12 occur as generic monodromy
of polynomial Keller self-maps of ordinary affine three-space. These
are attributed constructions, not new realizations claimed here.

Source: K. Kistner and T. Shaska,
[Orbits and fields of definition for graded Keller maps,
arXiv:2608.02863v1](https://arxiv.org/html/2608.02863v1), Sections 4--6.
The versioned [source manifest](sources/alternating-realizations-manifest.json)
pins the retained HTML. The algebra and monodromy arguments below are
rederived, with exact controls in
[the checker](verify_alternating_realizations.py) and
[its receipt](alternating_realizations_receipt.json).
No formal verification or independent whole-proof review is claimed.

## 1. A determinant identity for a global polynomial lift

Fix p>=1, set u=x^p y, v=x^p z and L=1+u. Let A,B be polynomials in
L,v which vanish at (1,0). Then

\[
 F=(A(L,v)/x^p,\ B(L,v)/x^p,\ xL)
 \tag{1}
\]

is polynomial on all of A^3. Indeed, after L=1+u, every monomial
u^i v^j of A or B has i+j>=1 and becomes divisible by x^p.
There are no deleted source divisors in this argument.

Put Q=L^p A, P=L^p B, C=xL. The changes from (x,y,z) to (L,v,x)
and from (F_1,F_2,C) to (P,Q,C) give

\[
 \det JF=-\frac{\operatorname{Jac}_{L,v}(P,Q)}{L^{2p-1}}.
 \tag{2}
\]

This equality holds first in the rational function field; once the
coordinates (1) are polynomial it is an identity on all affine space.
Moreover

\[
 k(x,y,z)=k(L,v,C),\qquad
 k(F_1,F_2,C)=k(P,Q,C).
 \tag{3}
\]

The rational inverse is x=C/L, y=(L-1)/x^p, z=v/x^p.
Here C is independent over k(L,v), and is independent over k(P,Q)
when P,Q are algebraically independent. Consequently both generic degree
and the normal-closure permutation group of the plane field pair
survive unchanged in the actual polynomial threefold map.

## 2. The beta--Hermite family, including an actual A5 map

For an integer p>=2 define

\[
 M_p(t)=\int_1^t(\theta^2-1)^{p-1}\,d\theta.
\]

It has degree 2p-1 and a zero of exact order p at 1.
Its values at 0 and -1 are nonzero: the real integrand has constant
nonzero sign on (0,1) and (-1,1). Thus there is a root
zeta outside {-1,0,1}; it is simple. Put delta=zeta^2-1, which is nonzero.
Choose a polynomial g(L) with

\[
 g(L)^2\equiv1+\delta L\pmod{L^p},\qquad g(1)=\zeta.
 \tag{4}
\]

The truncated square-root series through degree p-1, with a suitable
multiple of L^p added, gives such a polynomial. Define

\[
 \rho=g(L)+L^p v,\qquad
 Q=\frac{\rho^2-1-\delta L}{2},\qquad
 P=\frac{2p}{\delta^p}
     \int_1^\rho(\theta^2-1-2Q)^{p-1}\,d\theta.
 \tag{5}
\]

The integral is a polynomial expression, treating Q as constant before
substituting its value. Equation (4) gives L^p|Q.
To see L^p|P, expand the integrand in Q: the term of degree j in Q
contains an integral vanishing to order p-j at rho=1.
Since rho-1 is divisible by L, and Q by L^p, its L-order is at least
p-j+pj>=p. Write Q=L^p A and P=L^p B.
At (L,v)=(1,0), rho=zeta, Q=0 and P is proportional to M_p(zeta)=0.
Thus A and B vanish there, and Section 1 proves global polynomiality.

Direct differentiation gives

\[
 \operatorname{Jac}_{L,v}(\rho,Q)=\frac{\delta}{2}L^p,\quad
 \left.\frac{\partial P}{\partial\rho}\right|_Q
       =\frac{2p}{\delta}L^{p-1},\quad
 \operatorname{Jac}_{L,v}(P,Q)=pL^{2p-1}.
 \tag{6}
\]

The resulting polynomial map (1) therefore has **det JF=-p**.

The reconstruction

\[
 L=(\rho^2-1-2Q)/\delta,\qquad
 v=(\rho-g(L))/L^p
 \tag{7}
\]

proves k(L,v)=k(rho,Q). At fixed Q, P has degree 2p-1 in rho,
so (3) gives actual generic degree 2p-1.

Over the algebraic closure of k(Q), the derivative of this polynomial
is a nonzero constant times
(rho^2-1-2Q)^(p-1). It has exactly two distinct critical points.
Their values are distinct: their difference is proportional to
h^(2p-1) times the nonzero integral of (t^2-1)^(p-1) from -1 to 1,
where h^2=1+2Q. The finite inertia consists of two p-cycles, and
infinity has one (2p-1)-cycle.

These two finite inertia generators generate a transitive group.
Each moves p letters; transitivity forces their supports to cover all
2p-1 letters and hence meet in exactly one letter.
Two p-cycles meeting in one letter generate A_(2p-1) if p is odd and
S_(2p-1) if p is even. Here is an elementary group proof.
Label them a=(0,1,...,p-1) and b=(0,p,...,2p-2).
The commutator aba^-1 b^-1 is the three-cycle (0,1,p).
Its conjugates by powers of b have connected overlapping supports
covering supp(b) union {1}, and then their conjugates by powers of a
cover every letter. The three-cycle lemma in Section 3 gives A_(2p-1).
Parity of a decides between the alternating and symmetric groups.
For p=2 the two transpositions give S3 directly.

This first computes the group after algebraic closure of k(Q).
The group over k(P,Q) can in principle be larger. If p is even it is
already the full symmetric group. If p is odd, the derivative in (6),
as a polynomial in rho at fixed Q, is a constant times a square.
The polynomial discriminant of P(rho,Q)-T is therefore a square in
k(Q,T), up to a nonzero constant in k; this constant is a square since
k is algebraically closed. Thus the group over k(P,Q) is contained in
A_(2p-1) as well. Equation (3) completes the proof for the actual map.

**KMFI12.** For every p>=2 the family realizes S_(2p-1) for even p
and A_(2p-1) for odd p. In particular it realizes every natural
A_(4j+1), j>=1, absolutely on A^3. The exact checker instantiates p=3:

\[
 \zeta=\frac{-9+\sqrt{-15}}6,\quad \delta=\zeta^2-1,\quad
 g=1+\tfrac{\delta}{2}L-\tfrac{\delta^2}{8}L^2+
       (\zeta-1-\tfrac{\delta}{2}+\tfrac{\delta^2}{8})L^3.
\]

Equations (1),(5) are an explicit A5 Keller map over Q(sqrt(-15)),
with determinant -3 and degree five. The
[model](alternating_realizations_model.json) retains its exact A,B
polynomials and the checker verifies lift divisibility, quotient
Jacobian and the rational inverses. The proof for all p is the
written argument above; finite controls are not its substitute.

## 3. Three-cycle generation lemma

If a group is generated by three-cycles, its orbits are exactly the
connected components of the hypergraph of their three-letter supports.
On each component it is the full alternating group.

For completeness, start with one triple S. Suppose A(S) is available
and a further triple meets S. If it adds one letter, conjugation by
A(S), with inversion when needed in |S|=3, supplies the three-cycles
on the enlarged set. If it meets S in one letter and adds two letters
d,e, conjugate (a,d,e) to (b,d,e), with a!=b in S.
Their product, with one inverted, is a three-cycle on a,b and one of
d,e. Add that letter by the previous case, then add the other.
Induction along connected triples proves the claim.
This elementary lemma avoids needing a classification of primitive groups.

## 4. The degree-twelve alternating shear

Let

\[
 c=\frac{9-3i}{8},\quad d=\frac{1664+1152i}{125},\quad
 L=1+xy,\quad w=xz-1,\quad
 a=1+cLw^3,\quad b=w+dLa^3,\quad A=ab.
\]

The exact polynomial B, in these factored variables, is

\[
\begin{split}
B={}&a^2+2cLab^3-\tfrac{c^2}{3}L^2b^6
 +6cdL^2a^4b^2+6c^2dL^3a^3b^5-2cd^2L^3a^7b\\
 &+15c^2d^2L^4a^6b^4+\tfrac25cd^3L^4a^{10}
 -\tfrac{20}{3}c^2d^3L^5a^9b^3+3c^2d^4L^6a^{12}b^2\\
 &-\tfrac67c^2d^5L^7a^{15}b+\tfrac19c^2d^6L^8a^{18}.
\end{split}
\tag{8}
\]

Here d(1-c)^3=1 and 32c^2-72c+45=0.
At the lift point (L,w)=(1,-1), a=1-c, b=0,
A=0 and B=(32c^2-72c+45)/45=0. Hence

\[
 F_{12}(x,y,z)=(A/x,\ B/x,\ xL)
 \tag{9}
\]

is polynomial on ordinary affine three-space: after substituting L,w,
A and B vanish at x=0 identically in y,z, so both are divisible by x.

Set X=La^2, Y=Lab=LA and

\[
\begin{split}
\mathcal P(X,Y)={}&X+2cX^{-1}Y^3-\tfrac{c^2}{3}X^{-3}Y^6
 +6cdXY^2+6c^2dX^{-1}Y^5-2cd^2X^3Y\\
 &+15c^2d^2XY^4+\tfrac25cd^3X^5
 -\tfrac{20}{3}c^2d^3X^3Y^3+3c^2d^4X^5Y^2\\
 &-\tfrac67c^2d^5X^7Y+\tfrac19c^2d^6X^9,\qquad
 \Phi=X^2-c(Y-dX^2)^3 .
\end{split}
\tag{10}
\]

Substitution gives P=LB=mathcal P(X,Y), Q=LA=Y and the exact identities

\[
 \operatorname{Jac}_{L,w}(X,Y)=X,\quad
 \Phi(X,Y)=L^2a^3,\quad
 \mathcal P_X=\Phi^2/X^4.
 \tag{11}
\]

Thus Jac(P,Q)=Phi^2/X^3=L. Equation (2), with p=1,
proves **det JF_12=-1 globally**. This is a full polynomial Jacobian
proof through the chain rule, without expanding degree-eighty expressions.

There is also an explicit inverse field change:

\[
 L=\Phi^2/X^3,\qquad w=X(Y-dX^2)/\Phi.
 \tag{12}
\]

Consequently k(L,w)=k(X,Y), not a degree-two intermediate cover.
The double-cover shear used in the source paper has introduced no
extra sheets into the actual degree computation.

At generic Y!=0, mathcal P has poles of orders 3 at X=0 and 9 at
infinity, so its degree is twelve. Its finite critical points are the
six simple zeros of Phi, since

\[
 \operatorname{disc}_X\Phi=
 64Y^3c^4d^9(27Y^2cd+4)^2 .
\]

All have local degree three. Their critical values are generically
distinct, as follows from an exact nonzero specialization.
The scaling X=sX', Y=ds^2Y', mathcal P=s mathcal P',
where cd^3s^4=1, reduces any nonzero c,d to c=d=1.
For c=d=Y=1 the resultant of X^3(mathcal P-T) and Phi is a
degree-six polynomial in T, with gcd one with its derivative.
The exact rational polynomial is retained in the model. Therefore the
generic resultant over k(Y) also has six distinct roots. This uses the
specialization to prove a polynomial nonvanishing statement, not to
infer a generic group from one numerical fibre.

The geometric curve group is transitive and generated by its six
finite inertia generators, each a single three-cycle; the infinity
generator is their inverse product and has type (9,3).
Section 3 gives A12, and Riemann--Hurwitz is exact:
6*2+8+2=22.

Again parameter-field enlargement must be checked. Let
D=X^3(mathcal P-T), whose leading and constant coefficients are
c^2d^6/9 and -c^2Y^6/3. Equation (11) gives

\[
 XD'=\Phi^2+3D,\qquad
 \operatorname{disc}_X D
 =-\frac{27}{c^4d^6Y^6}\operatorname{Res}_X(D,\Phi)^2 .
 \tag{13}
\]

The right side is a square over k(Y,T) because k is algebraically
closed. Thus the full plane-field group lies in A12 and contains
the geometric curve A12. Equation (3) now proves:

**KMFI13.** Equation (9) is an actual polynomial Keller map on A^3
over Q(i), with generic degree twelve and geometric normal-closure
group A12. As a direct noninjectivity control, exact directional
derivatives at x=0 give

\[
 F_{12}(1,-1,1)
 =F_{12}(0,10-3i,\tfrac72-\tfrac32i)=(0,1,0).
\]

Equation (13) has square class -3 over the coefficient field Q(i);
over that field the arithmetic group is S12, with constant extension
Q(i,sqrt(-3)). This distinction follows from (13) and geometric A12.
It does not assert that the actual chosen threefold map is defined
over Q: the Laurent quotient with normalized c=d=1 is defined over Q,
where its arithmetic constant extension is Q(sqrt(-3)).

## 5. Consequences and scope correction

[KMFI14--16](ALL_ALTERNATING.md) extend the occurrence result to
every natural A_n, n>=5, and give exact independent-curve
presentations for the odd family. The even power-two examples are
proved to have no such presentation.

A5 is no longer merely an allowed action awaiting an affine completion
in the occurrence table. Both A5 and A12 are actual simple-group
realizations in a proper subgroup of the symmetric group of the same
degree. The claim that all such realizations remain missing is too broad.
The family gives infinitely many more A_(4j+1) realizations.

This does not realize the natural degree-twelve PSL2(11) action:
A12 has order 12!/2, whereas PSL2(11) has order 660.
It also does not realize either exceptional degree-eleven action.
Those specified cases and the full spectrum remain open.

The parameter-dependent plane covers above cannot simply be called
independent stabilizations of a fixed one-variable cover. The
[global curve-field obstructions](CURVE_DESCENT_LCM.md) were proved for
that specified stabilization class. Their extension to arbitrary
higher-dimensional field pairs has not been established here.
No sufficiency theorem for the action-level necessary conditions follows.
