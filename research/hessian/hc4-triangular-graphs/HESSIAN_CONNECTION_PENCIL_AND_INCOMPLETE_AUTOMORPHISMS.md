# Hessian connection pencils and incomplete polynomial automorphisms

HC4MYMI0. 16 November 2026. The family below has an explicit polynomial
gradient inverse and arbitrary polynomial degree. Nevertheless EVERY
nonendpoint member of its natural connection pencil is curved and
geodesically incomplete. In particular its Hessian metric is incomplete.
These are countercontrols to a proposed completeness argument, not
counterexamples to HC4. Full HC4/JC2 remains OPEN.

## 1. Three connections that must be distinguished

Let h be a complex polynomial in n variables with constant nonzero Hessian
determinant. Set H=Hess(h), Phi=grad(h), and write d for the standard
source affine connection. Define

    C_i=H^(-1)*partial_i H,
    nabla^(lambda)=d+lambda*C,              lambda in C.       (1)

All coefficients are polynomial, since H^(-1) is polynomial. Every member
is torsion-free: C^k_(ij)=sum_l(H^-1)_(kl)*h_(lij) is symmetric in i,j.
Every member preserves the standard source volume, because

    tr(C_i)=partial_i log(det H)=0.

The endpoints d and nabla^(1) are flat. Indeed nabla^(1) is the pullback
by Phi of the flat target connection, and the Maurer--Cartan identity gives

    R^(lambda)=lambda*(lambda-1)*C wedge C.                   (2)

The Levi-Civita connection of the complex bilinear Hessian metric H is
nabla^(1/2), by the usual coordinate formula and symmetry of h_(ijk).
It is NOT the flat target-pullback connection nabla^(1).

For a primary terminology comparison, Suzuki,
[*An explicit formula for Hessian potentials on warped product manifolds*](https://link.springer.com/article/10.1007/s41884-025-00173-4),
Section2, Proposition2.3, records the standard dual-connection identity
D'=2*nabla^LC-D. That section was read. The elementary matrix derivation
above proves the identities needed here directly; no warped-product,
positivity, compactness or completeness theorem is imported.

## 2. The all-degree family and its explicit gradient inverse

Choose ANY a,b in C*, ANY D in C[x,y], and put

    u=x+a*y^2,       v=y+b*u^2,
    K(x,y)=(u,v),    h(x,y,r,s)=D(x,y)+r*u+s*v.                (3)

This is a composition of two elementary shears. Its Jacobian A=DK has
determinant1 and its inverse is

    y=v-b*u^2,       x=u-a*(v-b*u^2)^2.                       (4)

Write z=(x,y), w=(r,s). The gradient outputs are

    eta=K(z),       xi=grad(D)(z)+A(z)^t*w.

Recover z by(4) and then w=A(z)^(-t)*(xi-grad(D)(z)). These are POLYNOMIAL
formulas. Thus grad(h) is an automorphism for EVERY D, with no degree
bound. This is the elementary cotangent construction already used in KO3,
not a new construction of a Keller counterexample.

The Hessian and its inverse have block forms

    H=[ C    A^t ],   C=D''+r*u''+s*v'',
      [ A     0  ]

    H^-1=[ 0       A^-1            ].
         [ A^-t   -A^-t*C*A^-1     ]

Consequently det H=1. For real data its bilinear metric has signature(2,2)
everywhere: the lower block is a maximal null plane and the cross block
is invertible. Positive-definite completeness or rigidity theorems cannot
be applied to it. Arbitrarily large deg D gives arbitrary-degree examples;
the inverse and the subsequent base equations are unaffected by that degree.

## 3. The projected geodesic equation is independent of D

Let z(t),w(t) be a geodesic of nabla^(lambda). The upper block of H^-1
and the identity h_(w_i z_j z_k)=(K_i)_(z_j z_k) give

    z''+lambda*A^-1*K''[z',z']=0.                            (5)

There are NO w or D terms in(5). All remaining geodesic equations form
a linear inhomogeneous second-order system for w,w', with coefficients
depending on z,z'. To see this directly, the lower equations are

    w''=-lambda*A^-t*( D'''[z',z']
                +sum_i w_i*K_i'''[z',z']
                +2*sum_i w_i'*K_i''*z'
                -C*A^-1*K''[z',z'] ).                      (6)

Here C is affine-linear in w. Thus every base solution of(5) on a simply
connected parameter domain lifts with arbitrary initial w,w'. Elementary
existence for a linear holomorphic ODE applies on every compact subinterval
of a nonsingular base path. This is not a completeness claim for the base.

Use the GLOBAL polynomial base coordinates (u,y), rather than (u,v).
Differentiating u=x+a*y^2 and using(5) yields exactly

    u''=2*a*(1-lambda)*(y')^2,
    y''=-2*b*lambda*(u')^2.                                (7)

Both equations, including their coefficients, are independent of D.
They describe a connection on the entire (u,y)-plane, and the projection
of the full connection to this plane is affine. Equation(6) proves the
lift statement even when the zero section w=0 is not preserved.

## 4. Every nonendpoint parameter is curved and incomplete

Fix ANY lambda different from0 and1. Set

    alpha=2*a*(1-lambda),       beta=-2*b*lambda.

Choose c in C* with c^3=beta/alpha and put kappa=alpha*c^2, which is nonzero.
On a disc with the analytic branch L(t)=-log(1-t), L(0)=0, put

    u(t)=L(t)/kappa,       y(t)=c*L(t)/kappa,
    x(t)=u(t)-a*y(t)^2.                                    (8)

Then u'=1/[kappa*(1-t)], y'=c/[kappa*(1-t)], and direct substitution,
using beta=alpha*c^3, proves BOTH equations(7). This is a written
ALL-PARAMETER solution; no numerical integration or parameter search is used.

As real t tends to1 from below, |u(t)| tends to infinity. Since u is a
polynomial of the source coordinates, this base geodesic cannot extend
holomorphically across t=1. Equation(6) supplies a full geodesic on a
neighborhood of every compact portion of [0,1), consistently by uniqueness.
That lift also cannot extend across1, since its polynomial base coordinate
u diverges. Hence the full connection is geodesically INCOMPLETE for
EVERY a,b!=0, EVERY D and EVERY lambda not in{0,1}.

For real a,b,lambda with lambda not in{0,1}, choose the real cube root c;
real initial dual data give the same finite-time REAL incompleteness.
For the Levi-Civita case a=b=1,lambda=1/2, choose c=-1,kappa=1. If D=0,
the zero section is preserved and a particularly short full geodesic is

    (x,y,r,s)=(L-L^2,-L,0,0),       L=-log(1-t).             (9)

It is a geodesic of the everywhere nondegenerate unit-Hessian polynomial
metric, despite the explicit polynomial gradient inverse in Section2.

Curvature is just as uniform. In the (u,y) coordinates the only nonzero
base Christoffel symbols are

    Gamma^u_(yy)=-alpha,       Gamma^y_(uu)=-beta.

Their constant connection matrices have commutator

    R_(uy)=4*a*b*lambda*(1-lambda)*diag(1,-1).               (10)

It is nonzero for EVERY nonendpoint parameter. Since the full connection
projects to this base connection, it too is curved, for every D. In
particular constant Hessian determinant does not imply a flat Levi-Civita
metric. Polynomial coefficients, preserved volume, contractible whole
source and a polynomial gradient inverse do not prevent its incompleteness.

The endpoints are complete in this family: nabla^(0) is the standard
connection, while nabla^(1) is flat in the GLOBAL polynomial coordinates
Phi=grad(h). Their midpoint is therefore incomplete even though both
endpoints are complete. Completeness cannot be interpolated in this pencil.

## 5. The flat dual endpoint retains the original global issue

For any polynomial unit-Hessian h, geodesic completeness of nabla^(1)
is equivalent to invertibility of grad(h). This is the standard developing-
map completeness criterion, not a new proposed independent HC4 criterion.
Here is the direct argument to avoid replacing it by a metric theorem.

The polynomial fields V_i=H^-1*e_i are parallel for nabla^(1), commute,
and satisfy V_i(Phi_j)=delta_ij. If all dual geodesics are entire, every
integral curve of each V_i is entire, since it is a geodesic with that
initial tangent. Their commuting complete flows give a holomorphic C^n
action. Its orbits are open because the fields span every tangent space;
connectedness makes the action transitive. Equivariance under Phi sends
this action to target translation, so its stabilizer is zero. Phi is
therefore a biholomorphism. A bijective etale polynomial morphism is an
algebraic isomorphism (equivalently use the usual birational Keller/Zariski
main argument). Conversely a polynomial inverse transports all target
straight lines to entire dual geodesics.

This argument applies to the FULL complex source in any dimension. It
does not follow from completeness of the standard source connection,
from the Levi-Civita connection, or from a positive-definite comparison.
The retained KZ complete-flow theorem is a different PLANE-specific
boundary theorem; its classification is not extended here to dimension4.

## 6. Verification and the remaining objective

The checker supplies exact symbolic transcription controls for the Hessian
blocks, gradient inverse, projected full geodesic equations, arbitrary-D
independence, the explicit real Levi-Civita geodesic, and the all-parameter
velocity/curvature identities. It does not prove completeness by a finite
simulation. The all-degree and all-parameter conclusions follow from the
written block and ODE arguments above. No registered checker is replayed.

No source realization is created for a hypothetical finite Keller cover.
No automatic completeness of the dual endpoint is proved. These examples
have degree-one field extension, so they are not new nonproperness packets
or counterexamples to the retained trace/energy criteria. This result
retires automatic Hessian-metric flatness/completeness and interpolation
as universal shortcuts; it does not reject every geometric approach.

HC4 and the plane Jacobian conjecture remain open.
CURRENT LAST UNSOLVED IMPLICATION remains

    ACTUAL polynomial plane source
       ==> Tr(Lambda_L*Lambda_L) subset B.                 OPEN

Full marked-energy integrality is also OPEN. Published generic-degree-four
plane exclusion remains retained. No independent replay, formal verification
or external review is claimed; review and outreach remain DEFERRED.
