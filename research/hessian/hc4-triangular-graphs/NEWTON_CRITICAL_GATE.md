# Newton faces and finite ramification close whole nonlinear graph families

Let `y=R(x,p)` in the original Meng--Yang potential, with `L*N != 0`.
The three retained duals may undergo any polynomial automorphism over K[x],
with polynomial inverse. All statements are in characteristic zero and have
no bound on either polynomial degree or automorphism word length.
Roots and factors are taken over the algebraic closure. The arguments may
be proved over C: any proposed characteristic-zero example and its actual
dual inverse use finitely many coefficients whose field embeds in C.

**HC4MYGN2. Every separated graph `y=a(x)P(p)+b(x)`, with a nonzero and
P nonconstant, forbids a nonzero constant actual four-variable Hessian after
every such automorphism.** The proof covers all degrees of P at once. It
uses the elementary critical-Hessian invariant, formal leading terms and
a polynomial derivative lemma; it needs no Tang or plane classification.
If a=0 or P is constant the graph is base-only, and the distinct unresolved
[unit-ideal endpoint](RELATIVE_COORDINATE_ENDPOINT.md) remains essential.

Three reusable gates below apply to general R as well: an empty critical
scheme with noncoordinate generic fibres, a simple root on the initial
Newton face, and a smooth finite ramification point. They narrow the
remaining genuinely nonlinear coefficient dependence, rather than opening
a sequence of increasing degree searches.

## 1. An empty critical scheme is also an obstruction (HC4MYGN0)

For any polynomial f(s,z1,z2,z3), suppose its entire relative critical
scheme is empty and its generic three-dual value fibre has compactly
supported Euler characteristic different from one. Then f cannot have
unit Hessian. Both hypotheses transport through actual relative dual
automorphisms.

Here is the complete argument. If its gradient were an automorphism, the
inverse image of the target line `(tau,0,0,0)` would be A1, not empty.
Otherwise the retained [Tang alternative](../tang-hc4-audit/STRUCTURAL_REVIEW.md)
supplies a jointly affine constant two-plane. Its projection to the single
free base has rank zero or one, as in [HC4MYGR1](CRITICAL_COMPONENTS.md).
In rank zero the form is `D(s,U)+P F(s,U)+Q G(s,U)` with (F,G) Keller.
The zero fibre of (F,G) is finite, so for generic s its two polynomials in
U are coprime. A Bezout completion of the linear dual form makes every
generic value fibre A2, of Euler one. In rank one the form is
`D(U,Z)+s A(U,Z)+V C(U,Z)` with (A,C) Keller. Its relative critical
scheme is isomorphic to the entire plane curve C=0: the other two
equations solve uniquely and polynomially for s,V using the constant
Jacobian. This curve is nonempty, since C is a nonconstant polynomial
over C. Both ranks contradict a hypothesis. No theorem about a single
rational fibre is substituted for this rank argument.

For the present graphs put

\[
                         F(x,p)=3+2xR(x,p).
\]

Assume deg_p R>0. For generic fixed x, the rational dual coordinates
`(p,T=A,q)` identify the actual three-dual potential with
`L T^2+c(p)T+e(p)+d(p)q`, where d has r>=1 distinct roots. A generic
value fibre has Euler `(1-r)+2r=1+r>=2`: off the roots it is
`(A1 minus r points) times A1`, while each root supplies two affine
lines. This is the previously proved [Euler stratification](CRITICAL_COMPONENTS.md),
valid for arbitrary root multiplicities.

If R0=R(0,p) is nonzero, has no simple root, and
`R0(-M/(2L)) != 0`, the boundary has no relative critical point. Indeed
its equations are `N R0=0` and `H'(p)+N R0'(p)q=0`; every root of R0
is multiple and avoids the zero of H'. By [HC4MYGK1](RELATIVE_CRITICAL_LOCUS.md)
the rest of the critical scheme is exactly `F=0` on `x*R_p != 0`.
Consequently if **every p-dependent irreducible factor of F over C(x)
has multiplicity at least two**, the whole critical scheme is empty
and unit Hessian is impossible. Vertical factors do not supply points
on this open set. This closes all such entirely repeated curves,
not merely a bounded collection of perfect powers.

For example `R=3p^2+(3/2)xp^4`, with M nonzero, satisfies
`F=3(1+xp^2)^2` and has no critical points anywhere. In contrast
`R=3p+(3/2)xp^2` has a simple boundary root and an actual critical
line; it is outside this empty-scheme theorem. Its earlier gradient
invertibility conclusion is retained without claiming a nonexistent
second component or universal unit-Hessian nonexistence.

## 2. A simple initial Newton root fixes the slope (HC4MYGN1)

Write `R=sum_j a_j(x)p^j` and define, over positive j with a_j nonzero,

\[
 \lambda=\min_j\frac{1+\operatorname{ord}_0 a_j}{j}>0,\qquad
 \mathcal P(z)=3+2\sum_{1+\operatorname{ord}_0a_j=j\lambda}
                         a_{j,\operatorname{ord}_0a_j}z^j.        \tag{1}
\]

The constant coefficient a0 of R does not enter: `x a0(x)` has positive
order. Thus P has constant term 3 and all its roots are nonzero.
For any simple root zeta of P choose positive integers e,m with
lambda=m/e, set `x=t^e`, `p=t^-m z(t)`, and lift z(0)=zeta in
`F(t^e,t^-m z)=0`. All powers of t are nonnegative. At each successive
order the unknown coefficient of z(t) has nonzero linear coefficient
P'(zeta), so recursion constructs the formal root. This elementary
implicit lifting suffices; no classification of all Puiseux branches
is needed. It is an actual algebraic branch, since F is polynomial.

On this branch

\[
 R_p=\tfrac12\mathcal P'(\zeta)x^{\lambda-1}(1+o(1)),\qquad
 p=\zeta x^{-\lambda}(1+o(1)).                                  \tag{2}
\]

R_p is nonzero, so this gives a genuine open relative critical branch.
The [exact native critical-Hessian formula](RELATIVE_CRITICAL_LOCUS.md) is

\[
 K_f=-32LN^2 x^6 R_p^2 D'',\quad
 D=-\frac{Np}{4x^2}-\frac{(M-13N/(2x^2))^2}{4L},                  \tag{3}
\]

where prime is differentiation along the branch with respect to x.
Any allowed unit-Hessian pullback requires the same nonzero constant
K_f on every branch, by HC4MYGK0.

For lambda<2 the universal leading term of D is
`-169*N^2/(16*L)*x^-4`. Hence

\[
 K_f=1690N^4\mathcal P'(\zeta)^2
                         x^{2\lambda-2}(1+o(1)).                 \tag{4}
\]

For lambda>2 the p term dominates D and (3) has positive order lambda.
For lambda=2, D has pole order at most four, even if its two leading
terms cancel, and (3) has order at least two. Thus a simple root on
the initial face **forces lambda=1**. At that slope it also forces

\[
                 K_f=1690N^4\mathcal P'(\zeta)^2                 \tag{5}
\]

for every simple face root.

If P is squarefree, (5) forces deg P<=2. Indeed, for a squarefree
polynomial P with P'(zeta)^2=c nonzero at every root, write
`(P')^2-c=P Q`. Differentiation and evaluation at all roots give
`Q=2P''`, since their difference has degree less than deg P. A second
differentiation gives `P P'''=0`, so P'''=0. This is a proof in all
degrees, not a check of several small face polynomials.

In particular, if the highest p coefficient a_n has a_n(0) nonzero
and deg_p R=n>=2, (1) is `3+2a_n(0)z^n` with lambda=1/n. It has
simple roots and unit Hessian is impossible. This includes every
monic-in-p graph in every p degree at least two. More generally, if
R0 is nonconstant of degree d>=2, lambda<=1/d<1; the mere existence
of **one** simple initial face root already forbids unit Hessian.

## 3. Close every separated polynomial dependence (HC4MYGN2)

Let `R=a(x)P(p)+b(x)`, with `n=deg P>=1`, `k=ord_0 a` and leading
coefficient P_n. Its initial face is exactly

\[
 \lambda=(k+1)/n,\qquad \mathcal P(z)=3+2a_k P_n z^n.
\]

All roots are simple. Section 2 therefore requires k=n-1 and n<=2.
The n=1 case is already excluded in every coefficient degree by the
elementary [HC4MYGK2](RELATIVE_CRITICAL_LOCUS.md), without Tang.

For n=2, k=1, there are two formal critical roots p_+,p_- with
`p_+-p_-=c/x+o(x^-1)`, c nonzero. For a quadratic polynomial in p
the squares `R_p(x,p_+)^2=R_p(x,p_-)^2` are exactly the same rational
discriminant, not just equal leading coefficients. Equation (3)
with the same constant on both branches therefore requires

\[
                  \left((p_+-p_-)/x^2\right)''=0.
\]

But its leading term is `12*c*x^-5`, nonzero. This contradiction
closes n=2 and proves the whole separated family. Zeros of a away
from zero, higher coefficients of a,b, repeated roots of P, and
arbitrary degrees do not affect this argument. The only roots required
simple are those of the initial binomial face, which are always simple.

## 4. Smooth finite ramification gives a pole (HC4MYGN3)

There is another general gate, even when the initial face has only
multiple roots. Suppose F=0 has a smooth point `(x0,p0)` with x0
nonzero, F_p=0, and its local irreducible component is nonvertical.
Then F_x is nonzero and we may use `t=p-p0` as local parameter,
with `x=x0+alpha*t^e+...`, alpha nonzero, e>=2 finite. A component
on which x is identically constant is explicitly excluded.

On this germ, R_p has order e-1, since
`2x R_p=F_p=-F_x dx/dt`. In (3) the critical value is
`D=D0-N*t/(4*x0^2)+O(t^2)`. The leading term of D'' has order
`1-2e`. Hence K_f has a **simple pole**, with leading coefficient

\[
              2LN^3 x_0^2 F_x(x_0,p_0)^2(1-e)\,t^{-1},           \tag{6}
\]

which is nonzero. A constant critical Hessian is impossible. The point
itself has R_p=0 and is not a relative critical point; the obstruction
is the pole along the genuine open critical branch approaching it.
No nonexistent finite critical point is inferred from this limit.

For instance `R=3p^2+(3/2)xp^4+(3/2)xp^2` has a repeated initial
face, while F is generically squarefree. At x=-4,p=i/2 it has
F=F_p=0, F_x=3 and F_pp=-96. Thus e=2 and this otherwise surviving
leading-face pattern is excluded by the finite gate.

## 5. What is still required

For a remaining nonautomorphic unit-Hessian p-only y graph, the
previous boundary constraints still hold: R0 is nonzero, has no
simple root and avoids `p=-M/(2L)`. Now F must also have a
p-dependent factor of multiplicity one. Its critical branches must
avoid the finite smooth ramification obstruction. If R0 is nonconstant,
every root of its initial face must be multiple. If R0 is constant,
a face with any simple root must have lambda=1; if the whole face is
squarefree, it is necessarily `3+c*z^2`. Mixed multiple/simple faces
and later splitting of repeated leading roots are not declared closed.
Singular finite branch points and branches escaping to p=infinity
at finite x require their own analysis; this note proves no exhaustion
of those cases. Separated coefficient dependence is fully closed.

[verify_newton_critical.py](verify_newton_critical.py) and its
[receipt](receipt-newton-critical.json) check exact coefficients, quadratic
discriminant identities, repeated-curve controls and finite-ramification
examples. The branch lifting, all-degree derivative lemma and degree
exhaustion are written proofs. The empty-scheme theorem imports retained
Tang; the Newton, finite-pole and separated-family arguments do not.
No new rational-fibre theorem, formal verification, independent external
review or novelty priority is claimed. The unit-ideal plane-surjectivity
and separate singleton-zero properness endpoints, q/r-dependent base
graphs, free-base/dual mixing, and unrestricted HC4/JC2 remain open.
