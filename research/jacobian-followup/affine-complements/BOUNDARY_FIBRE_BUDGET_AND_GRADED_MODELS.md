# A boundary-fibre budget and an unbounded affine test family

15 November 2026. PF2RULED51 uses the FULL valuation matrix to strengthen
PF2RULED50 beyond its pairwise tests. It also constructs actual affine
minimal complementary models with arbitrarily many unramified pieces,
and gives a primitive-element and norm-family gate for any finite Keller
map on those models. No such finite map is constructed.

## 1. One boundary value, with one explicit exception

Keep an actual plane Keller finite normalization, every omitted prime,
an index-two prime R and all unramified primes D_i. Suppose R has a
coordinate pencil distinct from their shared pencil. PF2RULED50 makes
them complementary; use source coordinates (x,y), y|R=0, v_R(y)=r,
v_R(x)=-q. Write

    x|D_i=lambda_i, v_(D_i)(x-lambda_i)=s_i,
    v_(D_i)(a source coordinate complementary to x)=-p_i.

The boundary constants lambda_i are finite. Let n_R and n_i be the
primitive diagonal sizes of PF2RULED50. Their bounds and cross terms are

    0<n_i<=p_i^2, checkR.checkD_i=q*p_i.

If lambda_i!=lambda_j, then

    checkD_i.checkD_j=0.                              (1)

Indeed take P1_y times P1_x. The two extractions lie over distinct finite
x values on the section y=infinity. After extracting D_i, its dual has
coefficient zero on that section: pairing the principal divisor x-beta
with it, for a generic beta, gives -v_(D_i)(x-beta)=0. Pushing the dual
of D_j back multiplies that coefficient by p_j, exactly as in the cross
pairing proof of PF2RULED50. This proves (1) at arbitrary orders and with
all Laurent coefficients present.

Group the D_i by lambda. Write D_lambda for their primitive dual block,
P_lambda=diag(p_i), and

    A_lambda=-P_lambda^-1*D_lambda*P_lambda^-1,
    E_lambda=1^t*A_lambda^-1*1.

The full matrix must be negative definite by Xie's theorem. Its principal
blocks give positive definite A_lambda with diagonal n_i/p_i^2<=1.
The Schur complement at R is the exact strict condition

    n_R/q^2 > sum_lambda E_lambda.                    (2)

Cauchy--Schwarz in the A_lambda inner product gives

    E_lambda >= max_(i in lambda) p_i^2/n_i >=1.       (3)

Thus if k is the number of distinct boundary values,

    k<n_R/q^2.                                       (4)

PF2RULED50's minimal R possibilities now show:

* r=1,q=2 gives n_R/q^2=5/4, so all lambda_i coincide;
* r>=3,q=r-1 gives (r^2-r+1)/(r-1)^2<2, again one value;
* only r=2,q=1,n_R=3 can have two distinct values; three are impossible.

These are all-degree restrictions on a genuine finite normalization.
They need no simple-zero assumption on D_i and no target positivity.
They bound distinct coordinate values, NOT the number of omitted primes.
For R sharing the unramified pencil, no such bound is asserted.

## 2. Exact grouped matrices in the simple pole-two exception

In the r=2,q=1 row, suppose every D_i has s_i=1,p_i=2. At a fixed
value lambda the whole chart is y=v_i+d_i*(x-lambda)^-2. Distinct primes
have distinct nonzero d_i: the same coefficient would give the same
valuation. Two such primes share the two initial monomial blowups and
have cross pairing -2; the individual diagonals are -4. Hence a group
of size m has

    -D_lambda=2*(I_m+J_m),
    E_lambda=2m/(m+1).                               (5)

With group sizes m_1,...,m_k, the sharp necessary condition is

    sum_j 2*m_j/(m_j+1)<3.                            (6)

For two groups either one is a singleton, or one has size2 and the other
has size at most4. This follows directly from
1/(m_1+1)+1/(m_2+1)>1/2; it is not a numerical cutoff experiment.
For one group every number of components passes (6). The next section
shows that actual smooth affine exact-volume geometry exists for all
of those component counts.

There is also a component-count consequence outside the exception.
Suppose every unramified D_i has simple zero and R has order r>=3.
PF2RULED50 forces q=r-1 and p_i=2; Section1 puts all D_i over one
coordinate value. If n is their number, the full matrix requires

    2n/(n+1)<(r^2-r+1)/(r-1)^2,
    n<(r^2-r+1)/(r^2-3r+1).                          (6a)

Thus n<=6 for r=3, n<=2 for r=4 and n<=1 for every r>=5. These are
conditional component-count bounds. The unrestricted number of unramified
primes remains unbounded in the r2q1 affine control family below.

## 3. An explicit affine family with every number of unramified pieces

Choose n>=1 and distinct nonzero complex numbers C,d_1,...,d_n, with
C!=d_i. Put H(h)=product_i(h-d_i), T(h)=h*H(h). Inside the marked source
C[x,y], define

    h=x^2*y, a=x*(h-C), b=y*H(h), z=x*y*H(h),
    B_n=C[h,a,b,z].                                   (7)

This is the full polynomial valuation-intersection ring for a component
of orders v_R(x)=-1,v_R(y)=2 and n components of orders
v_(D_i)(x)=1,v_(D_i)(y)=-2. On the first component h-C has order1;
on D_i, h-d_i has order2. The complete weight decomposition, with
weights deg(x)=1,deg(y)=-2, is

    (B_n)_k=C[h]*a^k, k>=0;
    (B_n)_(-2l)=C[h]*b^l, l>=1;
    (B_n)_(-2l-1)=C[h]*z*b^l, l>=0.                   (8)

To verify completeness, a source weight-k polynomial is x^k*f(h).
For k<0, polynomiality requires h^ceil(-k/2) dividing f. Regularity
on R requires (h-C)^max(k,0) dividing f; regularity on every D_i requires
(h-d_i)^max(ceil(-k/2),0) dividing it. These conditions give exactly (8).
Weight components cannot cancel in a boundary valuation: their leading
coefficients have distinct powers of the independent boundary parameter.

There is the explicit determinantal presentation

    a*b=(h-C)*z, a*z=(h-C)*T(h), z^2=T(h)*b.          (9)

It is an integral smooth affine surface N_n=Spec B_n. Here is a global
proof, not just formal boundary charts. The two principal opens h-C!=0
and T!=0 cover it, because T(C)!=0. On them (9) is respectively

    x^2*b=T(h), x=a/(h-C);
    a*z=(h-C)*T(h), b=z^2/T(h).                       (10)

The first is smooth since T has only simple roots; the second is smooth
since at a=z=0 one has h=C and its h derivative is -T(C)!=0. Both
opens are integral and share a nonempty open, proving integrality of
(9) and its identification with the subring (7).

The marked map A2_(x,y)->N_n is an open embedding. Where x!=0 and
h!=C its inverse is y=h/x^2. Near h=0 it is y=b/H(h). On T!=0,z!=0
it is x=T(h)/z, y=z^2/(h*H(h)^2). These formulas cover the source.
The exact complement consists of the disjoint affine lines

    R: h=C,z=b=0, a arbitrary;
    D_i: h=d_i,a=z=0, b arbitrary.                    (11)

At R take an etale square root rho of h, rho|R=c with c^2=C. Then
w=z/(rho*H(h)) satisfies y=w^2, x=rho/w. The parameter
t=(rho-c)/w is regular, with a|R=2C*t. Thus the whole local embedding
is x=c/w+t,y=w^2, giving index-two volume order1 and minimal r2q1.
At D_i, x is a uniformizer and y=d_i*x^-2+v with v regular and a
nonzero linear coordinate along D_i. Thus each D_i has simple x zero,
pole2 and volume order0. The component opens are affine: deleting any
of these effective Cartier divisors from the smooth affine N_n is an
affine morphism, locally the complement of one equation.

The explicit rational one-form on the common function field

    theta=-2y*dx-x*dy, dtheta=dx wedge dy              (12)

is regular globally. In the displayed R and D_i coordinates it becomes
-2w^2*dt-2w*t*dw and -2v*dx-x*dv. Its divisor of volume is exactly R.
The grading is a genuine C* action on N_n with weights (h,a,b,z)=(0,1,-2,-1).
It is a SOURCE action; no action or ruling on a hypothetical target is inferred.

The full primitive matrix is

    M_RR=-3, M_Ri=2, M_ii=-4, M_ij=-2 (i!=j).

It is negative definite for every n, with

    det(-M)=2^n*(n+3),
    Schur complement=3-2n/(n+1)=(n+3)/(n+1)>0.        (13)

The Euler characteristic is n+2 by (11). If a finite Keller map on N_n
has zero excess on its sole ramified prime, the retained Euler budget
would require field degree n+3; with positive excess it would be larger.
This is only a conditional compatibility count. These surfaces supply
NO finite target algebra, target branch parametrization or Keller pair.

There is a genuine geometric exclusion at the omitted resonance C=d_i:
the source curve x^2*y=C is C*. Its two ends extend to finite points
of R and D_i, giving a nonconstant P1 map into affine N, impossible.
The distinct-value hypothesis in this construction cannot be removed.

## 4. A primitive monomial gate that DOES use the full finite map

More generally let a hypothetical plane Keller finite normalization have
a globally regular primitive source monomial h=x^a*y^b, a,b>0,gcd(a,b)=1,
with EVERY omitted prime vertical over a finite h value. Then h is a
primitive element of the actual extension C(x,y)/C(P,Q).

Proof. A generic fibre h=lambda is wholly the source C*, is closed in N
and maps finitely and immersively into the target. The induced finite map
to the image normalization is etale. Its completion is a map P1 to a
genus-zero curve; Riemann--Hurwitz and the two deleted source points make
the affine normalization C*. Thus the map is, after coordinates over
C(lambda), a power map on C*. The two target punctures are rational over
C(lambda), since they are the images of the rational source punctures.

If its degree were k>1, its deck subgroup mu_k would act on the generic
source fibre by

    (x,y)->(zeta^b*x,zeta^(-a)*y).

Both P and Q would be invariant on that generic fibre, and therefore as
source polynomials. The origin is fixed. The invertible Keller differential
there forces zeta^a=zeta^b=1, hence zeta=1 by coprimality, a contradiction.
So the curve residue degree is one. Since h is regular on finite N, its
norm polynomial chi(U)=Norm(U-h) is monic with target-polynomial coefficients.
For generic lambda, chi(lambda) is reduced and irreducible. A proper
intermediate field C(P,Q,h) would make it a proper power by the norm tower
formula. Therefore C(P,Q,h)=C(x,y), as claimed.

For (7), h=x^2*y satisfies every verticality hypothesis. If a finite Keller
pair P,Q exists in B_n, chi has degree the FULL Keller degree and its
specializations must satisfy, up to nonzero constants,

    chi(0)=H_x^2*H_y;
    chi(d_i)=J_i*H_i^2;
    chi(C)=G_R*J_C^f, f in {1,2}.                    (14)

Here H_x,H_y are the images of the retained source axes; H_i the image
of D_i; J_i the image of the other component of h=d_i; and J_C the image
of the other component of h=C. Coincident GLOBAL supports are not silently
excluded; the local branches at the crossings over 0,d_i are transverse.

Indeed div(h)=2L_x+L_y, div(h-d_i)=S_i+2D_i and div(h-C)=S_C+R.
The curves S_i,S_C are A1, as follows from (10). All curves in the first
two rows map birationally and immersively: they are finite A1 maps in the
etale locus, and a higher-degree polynomial residue map has a critical point.
For S_C, the only possible critical point is its crossing with R. If its
residue degree f exceeds one, the polynomial residue map has that sole
critical point and local degree f there. At a nonimmersive R point the
kernel is tangent to R (PF2OMIT7), so transverse S_C is immersive and f=1.
At a fold, a transverse S_C has a target coordinate of order at most2.
Hence f<=2, proving the last row. Generic members are immersed C* curves;
no assertion that this is a linear pencil or a normal monogenic algebra is made.

## 5. Verification and the failed compact-curve shortcut

The [checker](verify_boundary_fibre_models.py) writes
[the receipt](boundary_fibre_models_receipt.json), with inputs pinned in
[the manifest](sources/boundary-fibre-models/manifest.json).
It checks grouped primitive matrices by actual point blowups, Schur
complements, full polynomial presentations and both localized smooth charts,
source inverses, boundary valuations, the explicit volume primitive and
weight decompositions. These controls support the universal written proofs.

A tempting polynomial for minimal r>=3 and simple pole2 D is
G=x^r*y^(r-1)-d*(xy)^(r-2)-C. Its R and D chart roots do not suffice:
putting u=xy gives x=(C+d*u^(r-2))/u^(r-1). Every root of
C+d*u^(r-2)=0 gives another end x=0,y=infinity, with finite nonzero u.
Those ends do NOT reach D, where u diverges. This invalidates the proposed
complete-curve exclusion, including when its two selected chart tests pass.
The extra ends are retained as negative controls, not hidden or discarded.

The new fibre budget is an exclusion; the affine family is actual geometry;
(14) is a necessary finite-algebra gate. Neither affine controls nor passing
matrices construct a finite Keller map. Common rulings, the remaining minimal
complementary cases and the full quartic/higher-degree problem remain open.
Independent full reproof, formal verification, external review and novelty
priority are not claimed.
