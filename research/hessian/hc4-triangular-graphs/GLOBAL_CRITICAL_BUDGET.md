# A global Euler budget for every remaining p-only base graph

Let `y=R(x,p)` in the Meng--Yang potential, `L*N != 0`, with any actual
polynomial automorphism of the three retained duals over K[x]. All degrees
and word lengths are unrestricted. Roots and curves may be taken over C;
a characteristic-zero candidate and its inverse descend through their
finitely generated coefficient field as in the preceding proofs.

This note combines the entire critical curve, its value functions and
generic dual fibres into one global test. It does not enumerate further
singularity types or assume that local critical-Hessian tests suffice.

Write `F=3+2xR`, and let r be the number of distinct roots of F over the
algebraic closure of C(x). Every p-dependent irreducible factor is counted
once, including factors whose multiplicity in F is greater than one.
Let S be the product of factors of multiplicity exactly one, as in
[HC4MYGM1](CRITICAL_NORMALIZATION.md). The other factors are called
inactive only with respect to this relative critical scheme.

We consider the remaining nonautomorphic boundary: R0=R(0,p) is nonzero,
has no simple root and `R0(-M/(2L)) != 0`. If deg_p R>0 and an allowed
pullback had unit Hessian, S=0 would be the entire nonempty smooth critical
curve C, with `(S,F_p)=(1)` and x etale to C*. Put

\[
 D=-\frac{Np}{4x^2}-\frac{(M-13N/(2x^2))^2}{4L},\quad
 \tau=\frac{dD}{dx},\qquad E=D-x\tau.                            \tag{1}
\]

These are regular functions on C. Let h be the generic degree of the
function `E+s*tau` on C, summed over its components, where s is a generic
constant independent of the curve coordinate x. Equivalently h counts
the points of `E+s*tau=lambda` for generic s,lambda.

## 1. Compare both generic value fibres (HC4MYGW0)

The necessary [entire plane Keller-fibre realization](CRITICAL_NORMALIZATION.md)
from HC4MYGM3 has the rank-one Tang form
`E(U,Z)+s A(U,Z)+V B(U,Z)`, with (A,B) Keller and an isomorphism
of the whole C with B=0 carrying tau to A. This uses retained Tang;
the rank-zero form was excluded by the native generic dual Euler value.

Fix generic s. On the plane complement of B=0, the value equation solves
uniquely for V. This stratum has Euler `1-chi(C)`. On B=0 the equation
is `E+s A=lambda`; its h points each supply a full affine V line.
Consequently the generic three-dual value fibre has Euler

\[
                            1-\chi(C)+h.
\]

The same fibre before the relative dual automorphism has Euler `1+r`:
the model `L T^2+c(p)T+e(p)+d(p)q` has an open stratum of Euler1-r
and two affine lines over every distinct root of d, irrespective of
its multiplicity. This is the retained [HC4MYGR2 calculation](CRITICAL_COMPONENTS.md).
Actual dual automorphisms preserve the entire value fibres. Thus

\[
                              h=r+\chi(C).                      \tag{2}
\]

The full curve is essential. Counting only one critical branch or omitting
inactive roots from r would give the wrong identity. Generic s is a scalar
fibre parameter, not the function x on C.

## 2. Eliminate the genus and sum the boundary (HC4MYGW1)

Normalize a projective closure of every component of C. Let the boundary
points over finite nonzero x where p escapes to infinity have local
orders `x=x0+t^e+...`, `p=t^-m*(unit)`, with e,m positive. Write
`epsilon=sum e` over all these escape points, not just base values.

The etale projection x on the closed affine C gives

\[
                             \chi(C)=-\epsilon.                 \tag{3}
\]

Here is an elementary Euler proof. Its generic degree is `a=deg_p S`.
Remove the finite set of escape base values from C*. Over that complement
the projection is a finite unramified a-sheeted cover. If there are b
removed values the Euler contribution is -ab. Over a removed value the
finite points are unramified and number `a-sum e` at that value; adding
these fibres gives -epsilon. There are no further missing finite points,
by HC4MYGM1. This argument needs neither the genus nor a classification
of the projective boundary. It does not assert that x is finite everywhere.

Combining (2) and (3) gives the global degree budget

\[
                               h+\epsilon=r.                    \tag{4}
\]

There is an explicit pole version. At any boundary point put
`b_q=max(0,-ord_q E,-ord_q tau)`. For generic s no leading cancellation
in E+s*tau occurs, so `h=sum b_q`. At a finite-base escape, D has
pole m, tau has pole m+e, and E has the same pole as tau because
x0 is nonzero. Thus `b_q=m+e`. All remaining boundary points lie over
x=0 or x=infinity. Equation (4) becomes

\[
 r=\sum_{\text{finite-base escapes}}(m_q+2e_q)
       +\sum_{x(q)=0\text{ or }\infty}
                      \max(0,-\operatorname{ord}_q E,-\operatorname{ord}_q\tau).
                                                                    \tag{5}
\]

Every summand is nonnegative. The earlier even-pole/conductor law still
holds at every escape; in particular each escape costs at least four
units of r. Repeated factors contribute to r only once and cannot buy
unbounded capacity merely by increasing their multiplicity.

At x=0 all branches have p tending to infinity, since F(0,p)=3.
Write x=t^e, p=beta*t^-m+... there. Unless `m=2e` and
`beta=-169N/(4L)`, the critical value D has pole at least4e.
Then tau has pole at least5e and the last sum in (5) costs at least5e.
Let `a_good` be the sum of e over these noncancelling active branches.
In particular

\[
                   r\ge5a_{\rm good}
                          +\sum_{\text{finite escapes}}(m+2e).  \tag{6}
\]

If all p-dependent F factors are simple and none of the x=0 branches
has the exceptional leading cancellation, r=a_good and (6) is impossible.
This closes that whole configuration in arbitrary degree. Exceptional
cancellation branches and inactive factors are retained in the general
budget (5), rather than silently discarded.

## 3. An algebraic differential descent gate (HC4MYGW2)

There is a complementary global test requiring no Tang. Let an irreducible
active component have function field L=C(x,p). If `R_p^2` on that
component belongs to C(x), a constant critical Hessian would force
`D'' in C(x)` by HC4MYGK1, hence `(p/x^2)'' in C(x)`.

An algebraic function u over C(x) whose second derivative is rational
is itself rational. To prove this, extend the derivation uniquely to a
normal closure. For every conjugate u_i, `(u_i-u)''=0`; its difference
is alpha_i*x+beta_i, with constants in C. The constant field is C:
a meromorphic function on a characteristic-zero normalized curve with
zero derivative is constant. Averaging the conjugates gives a rational
trace and then u rational. Applying this to u=p/x^2 shows that L=C(x).
Thus an active nonrational cover is impossible whenever its derivative
square descends to the base. This is a field condition on the actual
full R_p, not a claim about a factor derivative alone.

One complete consequence is that a squarefree quadratic F in p admits
no unit-Hessian pullback, with arbitrary polynomial coefficients. Its
two roots have the same rational R_p squared, so an irreducible cover
is excluded. If it splits into rational roots p_+,p_-, (1)'s common
constant Hessian forces `((p_+-p_-)/x^2)''=0`. Their difference is
therefore `x^2(alpha*x+beta)`. Both roots must have a pole at x=0;
they have the same principal part and pole order b>=1. Since the
constant coefficient of F is a unit there, its leading p coefficient
has order2b. Hence R_p has order at least2b+1 on both branches.
But D'' has pole order at most max(b+4,6), so
`ord K_f >= 4b+8-max(b+4,6)>0`, a contradiction.

For every graph with deg_p R=2, either this squarefree case applies or
F is entirely repeated in p and has no open critical curve. In the latter
case a simple boundary root gives gradient invertibility by HC4MYGR2;
a zero trace, a critical multiple root, or the remaining empty-scheme
case instead forbids unit Hessian by HC4MYGU0/HC4MYGR2/HC4MYGN0.
Thus **all quadratic p dependence has invertible gradient under actual
unit Hessian**, in all coefficient degrees and after every retained-dual
automorphism. The entirely repeated simple-boundary case is not
incorrectly relabeled universal unit-Hessian nonexistence. Base-only
graphs retain their separate unresolved unit-ideal endpoint.

## 4. An all-degree family passes the local gates but fails the global budget (HC4MYGW3)

Set L=N=1,M=0. For every odd k>=3 put l=(k-1)/2 and

\[
 G=3+\frac34x^{2k-2}p+\frac{507}{16}x^{2k-4},\quad
 Z=-\frac14x^{2k-3}p-\frac{169}{16}x^{2k-5},\quad
 Q=Z^l+(1+x)G/3,\qquad F=GQ^2,\qquad R=(F-3)/(2x).             \tag{7}
\]

These are actual polynomial graphs: G(0)=3, Z(0)=0, Q(0)=1.
Their boundary R0 is the nonzero constant -459/16 for k=3 and3
for every odd k>=5. On G=0, xZ=1 and Q=x^-l. Consequently G is
the only active factor and is the whole critical curve C*, with

\[
 p=-4x^{-(2k-2)}-\frac{169}{4}x^{-2},\quad
 R_p=\frac38x^{k-2},\quad D=x^{-2k},\quad
 K_f=-9k(2k+1)\ne0.                                            \tag{8}
\]

It passes the entire finite Bezout test: `(4/3)Z^(k-1)` represents
the inverse of F_p modulo G. There are no finite-base escapes.
All initial face roots are repeated: its face is
`3*(1+(-z/4)^l)^2` at slope2k-3. Thus the family satisfies the
constant critical-Hessian, finite normalization and local boundary gates
in every degree. Its critical target-line degree is2k+1>=7, outside
the previously used degree-one and degree-three exclusions.

For generic x, Q has l distinct roots in p and is coprime to G, so
`r=l+1`. To check squarefreeness without a bounded search write v=G/3;
then Q is `x^-l(1-v)^l+(1+x)v`. Its derivative obeys
`l*Q+(1-v)*Q_v=(1+x)(1+(l-1)v)`. A multiple root for l>1
requires v=-1/(l-1) and `x^l(1+x)=l^l/(l-1)^(l-1)`, only
finitely many x. The l=1 case is a generic linear polynomial.

But on C, `tau=-2k*x^(-2k-1)` and `E=(2k+1)*x^(-2k)`.
The generic function E+s*tau has degree h=2k+1 and epsilon=0.
The budget would require `2k+1=(k+1)/2`, impossible. Therefore
**every member of (7) forbids unit Hessian after every retained-dual
automorphism**, despite its nonzero constant native critical determinant.
The global comparison closes the family; passing the local tests was
never an existence claim. The native boundary Hessian already vanishes,
but that observation alone would not exclude nonlinear dual pullbacks.

## Verification and remaining obligation

[verify_global_critical_budget.py](verify_global_critical_budget.py) and its
[receipt](receipt-global-critical-budget.json) check the value-stratum
identities, pole counts, differential descent controls and exact family
identities. The all-degree Euler comparison, covering calculation and
conjugate argument are written proofs. Tang is an explicit retained
input for the global budget and the repeated-boundary quadratic
invertibility conclusion; the differential descent itself is elementary.
No new one-special-fibre rationality theorem, formal verification,
independent external review or novelty priority is claimed.

A remaining p-only graph must satisfy the full budget (5), as well as
the earlier finite Bezout and local conductor conditions, and must avoid
derivative-square descent on each nonrational active component. Those
conditions are still necessary, not sufficient for the entire abstract
plane Keller-fibre realization. Repeated-root splitting, exceptional
leading cancellation and complicated boundary configurations satisfying
this budget have not been exhausted. The two distinct plane endpoints,
q/r-dependent omitted bases, actual free-base/dual mixing and general
HC4/JC2 remain open.
