# Semisimple quadratic pencils forbid unit Hessians in every source orbit

15 November 2026. Claims **HC4MYJC0/1**. The graph degree, the degree
of the outer polynomial and the degree of every source reparameterization
are unrestricted. This gives an image obstruction without requiring
the intrinsic cotangent plane to be represented by ambient projections.
It does not classify nonsemisimple pencils or resolve HC4/JC2.

## 1. The general potential theorem (JC0)

Let q_A(x)=x^T A x be a nondegenerate complex quadratic form in n>=1
variables, q_B(x)=x^T B x a nonzero quadratic form, and H in C[t]
with deg H>=2. Assume

    T=A^(-1) B is diagonalizable over C.                 (1)

Then no polynomial automorphism Phi of A^n makes

    h composed with Phi,       h=q_A+H(q_B),             (2)

have a nonzero constant Hessian determinant.

Both quadratic forms are homogeneous about the same source origin.
Translations of the complete configuration and arbitrary constant
linear changes are allowed. Independent linear terms in q_A or q_B
are not included without a common centering certificate.

Only the retained critical-scheme transformation of
[HG0](GLOBAL_VALUE_ORBITS.md#1-global-invariants-and-the-whole-value-profile-hc4myhg0)
is needed. If det D Phi=delta in C*, the chain rule at each critical
point gives

    Hess(h composed with Phi)=D Phi^T Hess(h) D Phi.      (3)

Thus unit Hessian in (2) would force every critical point of h to be
nondegenerate and every native critical Hessian determinant to have
the same nonzero value. This is a statement on the critical scheme;
no nonlinear Hessian invariance elsewhere is assumed.

### Simultaneous normalization and repeated eigenvalues

The endomorphism T is self-adjoint for A. Its distinct eigenspaces
are A-orthogonal, since
lambda A(v,w)=A(Tv,w)=A(v,Tw)=mu A(v,w).
Diagonalizability and nondegeneracy of A make the restriction to
each eigenspace nondegenerate. Diagonalize that restriction and
rescale its basis over C. A constant source change therefore gives

    q_A=sum_j x_j^2,       q_B=sum_j lambda_j x_j^2.      (4)

This constant change multiplies all critical Hessian determinants by
the same nonzero factor and does not affect the required equality.

Write p=H', m=deg p>=1, and b=p(0). The origin is always critical.
Its Hessian is 2 diag(1+b lambda_j). If some 1+b lambda_j=0,
the origin is degenerate and the theorem follows.

Otherwise take a nonzero eigenvalue lambda. Every root r of
p(r)=-1/lambda is nonzero. On its eigenspace, the quadric q_B=r
consists of critical points of h, with all other coordinates zero.
If that eigenspace has dimension at least two, this is a nonempty
positive-dimensional critical variety, which is incompatible with
unit Hessian after any automorphism.

We may therefore assume every nonzero eigenvalue is simple. Denote
them lambda_1,...,lambda_k, k>=1; zero eigenvalues may repeat.
Set

    P(s)=product_{i=1}^k (1+lambda_i s),
    C_i=product_{j!=i}(1-lambda_j/lambda_i),
    s_i=-1/lambda_i.                                    (5)

All C_i and P(b) are nonzero. For each root r of p(r)=-1/lambda_i
there are two critical points on the i-th axis,
x_i^2=r/lambda_i. A multiple root would give a degenerate critical
point, so unit Hessian would make all m roots simple.

### Equal critical Hessians force an impossible scalar identity

Direct differentiation gives

    Hess h=2 diag(1+lambda_j p(q_B))
           +4 p'(q_B) (lambda_j x_j)_j (lambda_j x_j)_j^T.

At the origin its determinant is 2^n P(b); at the critical points
associated with i,r it is

    2^n * 2 lambda_i r p'(r) C_i.                        (6)

Equality with the origin therefore says r p'(r)=c_i for all roots
of p(r)+1/lambda_i, with
c_i=P(b)/(2 lambda_i C_i). The degree-m polynomial
t p'(t)-c_i vanishes at those m distinct roots. Comparing leading
coefficients yields

    t p'(t)-c_i = m (p(t)+1/lambda_i).                   (7)

All intermediate coefficients of p vanish. In particular

    p(t)=a t^m+b,  a!=0,
    c_i=-m(b+1/lambda_i).                               (8)

This conclusion already follows from one nonzero eigendirection.
The remaining eigendirections and (6) now require

    P(b)=2m(s_i-b)P'(s_i)          for every i.           (9)

The degree-k polynomial

    R(s)=2m(s-b)P'(s)-P(b)

has all k roots of P, and leading coefficient 2mk times that of P.
Consequently R=2mk P. At s=b this gives

    -P(b)=2mk P(b),

contradicting P(b)!=0 and characteristic zero. This proves JC0
in every n and every degree of H, without Tang, a Keller conjecture,
properness, an inverse-degree bound or a finite search.

For m=1 this is already a polynomial-in-a-quadratic obstruction.
For arbitrary m the same scalar identity suffices; increasing H's
degree creates no new case.

## 2. An image-only HC4 graph obstruction (JC1)

In an affine frame of A5 let

    Q(u,X)=q_A(u)+alpha X^2+beta X+gamma,
    det A!=0, alpha!=0,              u in A4.            (10)

This is a nondegenerate ambient quadratic. For any nonconstant
P0 in C[t], take the actual graph image

    S={X=P0(q_B(u))},                                  (11)

where B!=0 and A^(-1)B is diagonalizable. Every polynomial closed
embedding j:A4->A5 with image S has nonconstant or zero actual
Hessian determinant for Q composed with j; it never has unit Hessian.

Indeed the graph map g(u)=(u,P0(q_B(u))) and its inverse projection
are polynomial. The given j equals g composed with a polynomial
automorphism Phi of A4, since both identify A4 with that same
closed image. The native potential is

    Q composed with g=q_A+H(q_B),
    H(t)=alpha P0(t)^2+beta P0(t)+gamma,

whose degree is 2 deg P0>=2. JC0 applies. Thus no preferred source
parameters or affine pullbacks need survive. The graph and source
parameterization degrees are both unrestricted.

The conclusion persists on the entire polynomial quadratic-similitude
saturation of these images. If an ambient automorphism U obeys

    Q composed with U=lambda_*+c(Q-lambda_*), c!=0,

then every embedding into U(S) is U composed with g composed with
some Phi. Its actual potential is a nonzero scalar multiple of
the former one plus a constant; its four-variable Hessian determinant
is scaled by c^4. The retained
[HP0](QUADRATIC_ISOMETRY_SATURATION.md#1-saturation-under-every-polynomial-quadratic-similitude-hc4myhp0)
therefore transfers the nonexistence result. The literal ambient
normalization in [HM0](AMBIENT_QUADRIC_GRAPHS.md#1-a-literal-ambient-automorphism-hc4myhm0)
also transfers it to genuine embeddings for the original Meng--Yang
potential, with the actual source potential unchanged.

The earlier native-graph theorem already gives gradient invertibility
under unit Hessian in native graph coordinates. The new conclusion
here excludes unit Hessian throughout this whole source orbit and
quadratic-similitude image saturation. Nonlinear source changes have
not been treated as Hessian invariant.

## 3. Why the spectral hypothesis is essential

In four variables put

    q_A=2xy+z^2+w^2,       q_B=x^2.

Then A^(-1)B is a nonzero square-zero endomorphism, so it is not
diagonalizable. For every H, including arbitrary nonlinear H,

    h=2xy+z^2+w^2+H(x^2)

has Hessian determinant -16 and an explicit gradient inverse.
The polynomial shear

    y -> y+(H(x^2)-H(0))/(2x)

is legitimate because the numerator is divisible by x^2; it identifies
h with q_A+H(0). For a gradient target (a,b,c,d), recover

    x=b/2, y=(a-2x H'(x^2))/2, z=c/2, w=d/2.

This is an actual unit-Hessian positive control outside (1).
It forbids dropping diagonalizability or treating every nonzero
quadratic inner form as a nonexistence obstruction. Constant P0
also lies outside JC1: its image is an affine hyperplane and retains
the separately classified endpoint possibilities.

Nonsemisimple pencils, independently shifted quadratics, arbitrary
graph images and full genuinely non-descended embeddings remain
outside the theorem. Nothing says those classes contain a
nonautomorphic example.

## 4. Literature and route boundaries

The selected primary introduction of
[Ni, arXiv:2608.14217v1](https://arxiv.org/html/2608.14217v1),
Theorem 1.2 and its coordinate-change convention, concerns native
potentials of degree four. No theorem from that paper is imported:
arbitrary nonlinear source automorphisms here have unbounded degree,
and the argument above proves a different source-orbit nonexistence
statement. The polynomial graph result likewise permits arbitrary
P0 degree. No priority claim beyond this audit is made.

The rank-three helper identity suggested a matrix-factorization route
for D+K1K2-UV. The inspected abstract of
[Leuschke, arXiv:2609.37839v3](https://arxiv.org/abs/2609.37839v3)
states a stable-category equivalence over a complete hypersurface
ring. No passage from such an equivalence to polynomial coordinates,
actual closed-image recovery or ambient symmetry rectification has
been established here; no such implication is imported.

The only retained proof inputs for JC0/1 are the critical chain rule,
the literal quadratic normal form, and polynomial quadratic-similitude
transfer. Simultaneous normalization and the all-degree scalar
contradiction are proved above. The checker supplies finite exact
controls, not a proof by checking degrees, eigenvalues or source maps.
