# Every centered quadratic-composite graph returns to the old endpoint

15 November 2026. Claims **HC4MYJE0/1**. This closes the nonzero-Jordan
frontier left by JC0/1 and JD0. The outer polynomial, the graph degree
and every polynomial source parameterization have unrestricted degree.
The common-center hypothesis and the stated nondegenerate ambient
frame remain essential. The singleton-zero endpoint and unrestricted
HC4/JC2 remain open.

## 1. Statements and exact retained inputs

Let q_A(x)=x^T A x be nondegenerate in four complex variables,
q_B(x)=x^T B x, and H in C[t] have degree d>=2. Put T=A^(-1)B.

**JE0.** If T is not nilpotent, then every polynomial source
automorphism Phi for which g=(q_A+H(q_B)) composed with Phi has
nonzero constant Hessian determinant gives a polynomial gradient
inverse for g. For diagonalizable nonzero T, JC0 gives the stronger
nonexistence of unit Hessian. The new argument handles all remaining
nonnilpotent Jordan configurations.

**JE1.** In the graph family

    Q(u,X)=q_A(u)+alpha X^2+beta X+gamma,
    det A!=0, alpha!=0,          X=P0(q_B(u)),            (1)

all parameterizations with nonnilpotent T close to gradient
invertibility under unit Hessian. Nilpotent T has exactly the old
nondegenerate quadratic source orbit by JD0. Constant P0 has that
same orbit. Hence the universal inverse question for this entire
commonly centered graph-image family, including its certified
quadratic-similitude saturation, is exactly the old singleton-zero
plane Keller endpoint. It creates no further endpoint. This is not
a proof of that endpoint.

The retained inputs are:

- [HG0](GLOBAL_VALUE_ORBITS.md#1-global-invariants-and-the-whole-value-profile-hc4myhg0):
  if a four-variable unit-Hessian source reparameterization has
  nonautomorphic gradient, the generic whole value fibre has
  E=L^3-nL, L=uv, where n is its number of full critical points.
  At every critical point a source automorphism transforms the
  Hessian by congruence with one constant determinant factor.
- The additive/product Hodge--Deligne invariant and smooth affine
  curve formula, explicitly sourced in
  [HC0](HODGE_TRACE_CLOSURE.md#1-the-additional-invariant-and-its-exact-input).
- The retained characteristic-zero
  [Riemann--Hurwitz formula, Stacks 0C1B](https://stacks.math.columbia.edu/tag/0C1B),
  already used in IY0. The curve genera below are derived from it.
- The semisimple exclusion [JC0](SEMISIMPLE_QUADRATIC_PENCILS_AND_ALL_SOURCE_ORBITS.md),
  the nilpotent normal form and quadratic endpoint [JD0](NILPOTENT_QUADRATIC_PENCILS_AND_ENDPOINT_ORBITS.md),
  and the literal ambient transfers HP0/HM0.

No new unproved external theorem is imported. No multiplicativity
for an arbitrary finite cover is used. All affine-line factors below
come from solving an equation on explicit disjoint strata.

Write U=u+v, epsilon=d mod 2,

    rho=2-epsilon,      a=(d+epsilon-2)/2.

Choose a generic value kappa avoiding H(0) and all the finitely many
critical values of the one-variable polynomials which occur.

## 2. Two curve calculations, including every missing point

For lambda!=0 let F_lambda(r)=r+lambda H(r). The affine double cover

    D_lambda: v^2=F_lambda(r)-lambda kappa

is smooth, with d simple finite branch points. At infinity it has
epsilon ramification points and rho punctures: odd degree has one
ramified point; even degree has two unramified points. The local
unit factors have square roots in C[[1/r]], so these are all the
infinite places. Riemann--Hurwitz gives genus a. Thus

    E(D_lambda)=M:=L-aU+1-rho.                           (2)

Removing v=0 removes exactly d affine points.

For the plane curve

    C_{b,c}: z^2+w^2+H(b z^2+c w^2)=kappa,              (3)

the following exhaustive facts will suffice:

- If b=c=0, E(C_{0,0})=L-1.
- If b=c!=0, the equation in s=z^2+w^2 has d distinct
  nonzero roots. Thus E(C_{b,b})=d(L-1), by disjoint conics.
- If b!=c and at least one is nonzero, C_{b,c} is smooth
  connected and its completion has positive genus.

Here is a proof of the last assertion. If c=0,b!=0, write
w^2=kappa-z^2-H(b z^2). This is a squarefree polynomial of
degree 2d in z for generic kappa. Its double cover has genus
d-1>=1 and two punctures, so E=L-(d-1)U-1.

If bc!=0 and b!=c, put r=b z^2+c w^2. Then

    z^2=(r+c H(r)-c kappa)/(b-c)=A_1(r),
    w^2=(b kappa-b H(r)-r)/(b-c)=A_2(r).                 (4)

Both degree-d polynomials have simple roots; they have no common
root for generic kappa, since a common root forces r=0 and
kappa=H(0). Their nonempty disjoint simple root divisors make their
square classes independent, giving a connected degree-four cover
of the r-line. At the 2d finite branch points the total
ramification contribution is 4d. At infinity it is 2epsilon,
with 4-2epsilon punctures. Indeed the two odd-valuation classes
give the same local quadratic extension when d is odd; for even
d both split locally. Riemann--Hurwitz gives

    G=2d-3+epsilon>=1,
    E(C_{b,c})=L-GU+1-(4-2epsilon).                      (5)

The affine double-square equations in (4) are the whole original
curve: the linear relation eliminates r. They are smooth because
the two squarefree branch divisors are disjoint. Thus no finite
normalization points are inadvertently removed.

Finally consider the two-variable Jordan-block potential

    J_lambda(v,s)=2vs+H(2lambda vs+v^2).

For lambda=0 the shear in JD0 gives a nondegenerate quadratic,
so its generic fibre has E=L-1. If lambda!=0, its v=0 stratum
is absent for generic kappa. On v!=0 set r=2lambda vs+v^2.
The fibre is exactly D_lambda with v=0 removed, with inverse
s=(r-v^2)/(2lambda v). Therefore

    E(J_lambda=kappa)=M-d       (lambda!=0).              (6)

In particular the omitted branch points count d, not zero.

## 3. The shared off-slice calculation

For a length-two block of nonzero eigenvalue lambda, write

    q_A=2xy+q_C(z,w),
    q_B=2lambda xy+x^2+q_D(z,w).
    f(z,w)=q_C(z,w)+H(q_D(z,w)).

Assume the binary form lambda q_C-q_D is nondegenerate. On x!=0,
put r=q_B and eliminate y. The fibre equation is

    x^2=F_lambda(r)+lambda q_C(z,w)-q_D(z,w)-lambda kappa.

A constant binary change identifies the last quadratic with pq.
For the auxiliary threefold x^2-pq=F_lambda(r)-lambda kappa,
stratifying by p!=0 and p=0 gives E=(L-1)L^2+LM.
Its x=0 part is a Danielewski surface pq=-F_lambda(r)+lambda kappa.
Stratifying by p!=0 and the d simple roots gives E=L^2+(d-1)L.
Subtracting it gives

    E(off x=0)=L^3-L^2-aLU+(2-rho-d)L.                   (7)

On x=0 the original fibre is {f=kappa} times the free y-line.
Consequently

    E(h=kappa)=L^3-L^2-aLU+(2-rho-d)L
                         +L E(f=kappa).                 (8)

Every multiplication in this calculation comes from an actual
product stratum; the double cover itself is accounted for by
(2), not by treating it as two copies of its base.

## 4. Complete nonsemisimple block audit

Distinct generalized eigenspaces of the A-self-adjoint T are
A-orthogonal. For example T-lambda is nilpotent on its lambda
primary space and invertible on another primary space; moving a
sufficiently large power across the pairing proves orthogonality.
Their restrictions are nondegenerate. On a primary space write
T=lambda I+N and apply the elementary nilpotent chain normalization
proved in JD0. Each block has q_B=lambda q_A+q_N.

The four possible nonsemisimple partitions of four are therefore
2+1+1, 2+2, 3+1, and 4, with arbitrary eigenvalues on their blocks.
The following covers zero eigenvalues as well as nonzero ones.

Assume, towards contradiction, a unit-Hessian source
reparameterization with nonautomorphic gradient. Put b=H'(0).
The origin is critical, so

    det Hess h(0)=16 det(A) det(I+bT)!=0.                 (9)

In particular 1+lambda b!=0 for every eigenvalue.

We will use one critical-count observation. If every eigenspace for
a nonzero eigenvalue is totally A-isotropic, the origin is the only
critical point. Indeed grad h=2A(I+H'(q_B)T)x. A nonzero critical x
must be in a nonzero-eigenvalue eigenspace, hence q_B(x)=lambda
q_A(x)=0. Then (9) makes I+bT invertible, a contradiction.
Thus the origin is one reduced point and n=1 in HG0. This applies
where explicitly invoked below.

### Partition 2+1+1

Use

    q_A=2xy+z^2+w^2,
    q_B=2lambda xy+x^2+mu z^2+nu w^2.

If lambda=0, nonnilpotence means (mu,nu)!=(0,0).
On x!=0 the y-equation gives E=L^3-L^2; on x=0 use
C_{mu,nu} times A1. If mu=nu!=0, the whole E has coefficient
d-1>0 at L^2. If mu!=nu, it has a negative, nonzero LU term
from the positive genus in (3). Both contradict HG0.

Suppose lambda!=0. If mu=lambda or nu=lambda, choose a root r of
H'(r)=-1/lambda. By (9), r!=0. The corresponding eigenspace
contains the isotropic eigenline of the length-two block and a
nonisotropic simple direction. Its q_B=r locus has a free
isotropic coordinate, giving a positive-dimensional critical
variety. Unit Hessian is already impossible.

Otherwise lambda q_C-q_D is nondegenerate and (8) applies with
f defining C_{mu,nu}. If mu=nu!=0 the L^2 coefficient is d-1.
If mu!=nu with at least one nonzero, the LU coefficient is
-(a+G), with G>=1. These contradict HG0. If mu=nu=0, (8) is

    E=L^3-aLU-(d+rho-1)L.                               (10)

All nonzero eigenspaces are isotropic, so n=1. For d>=3,
a>=1; for d=2, a=0,rho=2 and (10) is L^3-3L.
Neither is L^3-L.

### Partition 2+2

Use

    q_A=2xy+2zw,
    q_B=2lambda xy+x^2+2mu zw+z^2.

If lambda!=mu, choose the nonzero one as lambda.
The binary form 2(lambda-mu)zw-z^2 is nondegenerate, and (8)
applies with f=J_mu. If mu=0, (6) gives exactly (10).
If mu!=0, (6) gives

    E=L^3-2aLU+(3-2rho-2d)L.                            (11)

Again every nonzero eigenspace is isotropic, so n=1. For d>=3
the LU term is nonzero; for d=2, (11) is L^3-5L.
Both alternatives contradict HG0.

If lambda=mu!=0, put r=q_B. Project the fibre to (x,z,r):

    x^2+z^2=F_lambda(r)-lambda kappa.

This Danielewski surface has E=L^2+(d-1)L. Its d points with
x=z=0 do not lift: in the original space x=z=0 forces
q_A=q_B=0 and kappa=H(0), excluded. Over the complement the
equation 2xy+2zw=kappa-H(r) has one affine-line fibre.
Partition into x!=0 and x=0,z!=0 to get actual products.
Thus the whole fibre has

    E=L^3+(d-1)L^2-dL,                                 (12)

contradicting HG0. The case lambda=mu=0 is nilpotent and
belongs to JD0, rather than this nonnilpotent theorem.

### Partition 3+1

Use

    q_A=2xy+z^2+w^2,
    q_B=lambda(2xy+z^2)+2xz+mu w^2.

If lambda=0, solve y on x!=0. If lambda!=0, put r=q_B and solve

    z=(F_lambda(r)+(lambda-mu)w^2-lambda kappa)/(2x),
    y=(r-lambda z^2-2xz-mu w^2)/(2lambda x).

Both give E=L^3-L^2 on that open. On x=0 the fibre is
C_{lambda,mu} times the y-line. Nonnilpotence means
(lambda,mu)!=(0,0). Equal nonzero eigenvalues give the forbidden
positive L^2 coefficient d-1; distinct ones give a nonzero
negative LU coefficient by (3). This covers the case of a
nilpotent length-three block paired with a nonzero simple one.

### Partition 4

Here lambda!=0 and

    q_A=2xw+2yz,
    q_B=lambda(2xw+2yz)+2xz+y^2.

On x!=0 put r=q_B and solve

    z=(F_lambda(r)-y^2-lambda kappa)/(2x),
    w=(r-2lambda yz-2xz-y^2)/(2lambda x).

The open has E=L^3-L^2. At x=0 the fibre is
{J_lambda(y,z)=kappa} times the w-line. Formula (6) yields
exactly (10). Its single reduced critical point gives n=1,
and the same genus/coefficient contradiction applies.

These are all the nonsemisimple partitions. Together with JC0
they prove JE0 for every nonnilpotent T and every d>=2.

## 5. Image-family conclusion and what stays open

For nonconstant P0 in (1), put
H=alpha P0^2+beta P0+gamma; its degree is 2 deg P0>=2.
Every closed parameterization of that graph is the graph map
composed with a polynomial source automorphism. JE0 therefore
closes its nonnilpotent cases. The nilpotent cases have exactly
the nondegenerate quadratic orbit by JD0. A constant graph is
already that orbit.

One fixed nilpotent example, q_A=2xy+z^2+w^2 and q_B=x^2
with any chosen nonconstant P0, realizes the whole quadratic
source orbit. By HG2 the universal inverse statement for all
these images is consequently equivalent, in both directions,
to the existing singleton-zero plane Keller endpoint.

The same equivalence holds after certified polynomial quadratic
similitudes, and after the literal HM0 transfer to the original
Meng--Yang potential: the actual source functions differ only
by a nonzero scalar and an additive constant. No normality
transfer, arbitrary ambient rectification or nonlinear Hessian
invariance is assumed.

Independent centers, more general graph geometry, nonconstant
higher helper traces and intrinsic ambient alignment remain
separate all-embedding gates. The present image family is fully
reduced to the old endpoint; that endpoint and unrestricted
HC4/JC2 are not resolved. The finite checker controls the exact
strata, inverse formulas, ramification arithmetic and profile
contradictions. The written identities and complete block audit
prove the unrestricted-degree statements.
