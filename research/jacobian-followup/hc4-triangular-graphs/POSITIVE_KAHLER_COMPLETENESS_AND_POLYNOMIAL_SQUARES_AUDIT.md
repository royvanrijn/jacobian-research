# Positive Kähler completeness and polynomial squares

16 November 2026. ROUTE AUDIT ONLY. Full HC4/JC2 remains OPEN.
No new ledger theorem, external theorem adoption, primary cache or
mathematical checker invocation. The arguments below screen a proposed
global positive-metric route in every degree, including arbitrarily many
polynomial helper functions. Classical pullback and Gram-matrix arguments
are used without a novelty claim.

## 1. The positive metric is different from the retained Hessian metric

For an ACTUAL polynomial plane Keller map K=(P,Q), normalized by det DK=1,
put

    u=|P|^2+|Q|^2,       g=K^*g_Euclidean.

In complex Hessian conventions, g_(i bar j)=sum_a (partial_i K_a)
overline(partial_j K_a); equivalently its usual Hermitian quadratic form
is the Gram form of DK. Thus g is positive definite, its coefficients
are polynomial in z and bar z, and det g=1. It is FLAT because K is a
local holomorphic isometry to Euclidean space. In particular its Kähler
Ricci tensor vanishes automatically.

This is NOT the complex bilinear Hessian metric of a cotangent potential
h=D+w.K studied in HC4MYMI0/HC4MYMJ0. Its positivity does not repair that
metric's incomplete Levi-Civita geodesics. Its automatic Kähler Ricci zero
does not supply MJ0's DIFFERENT Levi-Civita tensor hypothesis.

For this pullback metric,

    g complete  <==>  K has a polynomial inverse.              (1)

Indeed a complete local isometry from a connected complete Riemannian
manifold to connected Euclidean space is a covering onto that space.
One can see surjectivity by lifting straight segments: a finite-length
lift is Cauchy at a finite endpoint, so completeness and the local inverse
continue it. The usual covering argument, and simple connectedness of
Euclidean space, then give one sheet. A bijective étale polynomial map
is an algebraic isomorphism. Conversely an inverse makes g globally
isometric to the Euclidean target. These are standard completeness and
étale-isomorphism facts; no new completeness theorem is asserted.

The missing polynomial input cannot be dropped. The ENTIRE HOLOMORPHIC,
NONPOLYNOMIAL map

    F(z,w)=(exp(z),exp(-z)*w),       det DF=1

gives a positive flat metric of determinant one on all C2. Along
gamma(t)=(-t,0), t>=0, its length is integral_0^infinity exp(-t) dt=1
(up to a harmless fixed metric-convention factor). The path escapes,
so the metric is incomplete. This map has image C* times C and is NOT
a polynomial Keller/HC4 counterexample.

## 2. Growth and convexity hypotheses are substantial

For this particular potential, ANY coercive lower bound u(z)->infinity
as |z|->infinity already makes K proper. A proper étale polynomial map
to C2 is a finite covering of simply connected C2, hence an isomorphism.
Thus deriving such a lower bound from the Keller condition would already
solve the original inversion problem, before a Monge–Ampère theorem is
used. Adding a constant to u makes no difference.

Conversely known inverses need not obey near-quadratic Euclidean growth.
For EVERY m>=2 the shear

    K_m(z,w)=(z+w^m,w),        u_m=|z+w^m|^2+|w|^2

has an explicit polynomial inverse and a complete positive pullback
metric, but u_m(0,r)=r^(2m)+r^2. On the escaping curve
(z,w)=(-r^m,r), u_m=r^2 while |(z,w)|^2=r^(2m)+r^2. Consequently both
the near-quadratic upper bound and the near-quadratic coercive lower
bound, even with fixed logarithmic losses, fail in these coordinates.
This does not deny properness: r^2 still tends to infinity.

The same potentials are NOT real convex, even outside any compact set.
On the real slice, at w=1,

    partial_w^2 u_m
       =2*m^2+2*m*(m-1)*(z+1)+2.

It is negative for real z sufficiently negative, with z arbitrarily far
out. Strict plurisubharmonicity therefore supplies neither real convexity
nor real convexity outside a compact set. Linear coordinate changes do
not repair real nonconvexity. Choosing a polynomial target automorphism
that turns K into an affine map would already exhibit an inverse for K.

There is a further elementary obstruction to using a pluriharmonic
gauge to manufacture the needed upper growth. Let f_a be finitely many
holomorphic polynomials, let h be ANY entire holomorphic function, and
let u=sum_a |f_a|^2+Re(h). If

    u(z) <= C*(1+|z|^2)*(1+log(e+|z|^2))^beta

globally for a fixed finite beta>=0, then EVERY f_a has degree at most
one. To prove this, average over the sphere |z|=r. Re(h) has its constant
mean Re(h(0)). Homogeneous pieces of different degrees of each f_a are
orthogonal under scalar rotations of the sphere. The mean of the sum
of squares is therefore sum_k a_k*r^(2k), with a_k>=0 and a_k>0 for
each degree actually occurring. The displayed upper bound prohibits
any k>=2. This argument requires no determinant equation or PDE estimate.
For polynomial squares the cited upper growth already excludes nonlinear
Keller coordinates; a gauge cannot cancel their positive spherical mean.

## 3. Arbitrarily many polynomial squares retain a square Keller map

Here n and N>=n are ANY finite integers. Let f=(f_1,...,f_N) be polynomial
on C^n and assume the positive Gram form of Df has determinant c>0
constant. Cauchy–Binet gives the polynomial identity on the real locus

    sum_(I subset {1,...,N}, |I|=n) |det Df_I(z)|^2=c.       (2)

Every maximal minor is bounded, hence constant. Choose I0 with nonzero
minor and set K=f_I0. This is a square polynomial Keller map. For every
remaining row, Cramer's rule writes

    df_a=sum_j c_aj*dK_j,          c_aj in C,               (3)

because every coefficient is a ratio of constant maximal minors with a
fixed nonzero denominator. Subtracting the right side and differentiating
shows f_a-sum_j c_aj*K_j is constant. Hence

    f=C*K+b,           rank C=n,          g=K^*H,
    H=C^*C a fixed positive definite Hermitian form.         (4)

In particular this metric is already flat, and its completeness is
equivalent to inversion of K by (1), in ANY dimension. In dimension two
this is exactly the square Keller endpoint. Arbitrarily many helpers,
arbitrary polynomial degrees and positive constant weights (absorbed into
f) do not add a global inversion constraint. An added Re(h) does not
change the metric or (2). Conversely any square Keller map gives such
a constant-determinant positive metric, so this class retains the entire
original problem. This is a classical Gram/Cauchy–Binet reduction, NOT
a new proof of inversion or a claim that all polynomial positive
Monge–Ampère potentials have such a representation.

For example, trying to add epsilon*|z|^2 to the plane Keller potential
gives g_epsilon=g+epsilon*I and

    det(g_epsilon)=1+epsilon*tr(g)+epsilon^2.

For epsilon>0 this is nonconstant unless K is affine: constant tr(g)
is a constant sum of squares of the polynomial entries of DK and makes
all of them constant. The added source term ensures metric completeness,
but loses the constant-density equation precisely in the nonlinear case.
If the original K is retained among a larger finite polynomial tuple,
(3) also says every added helper must be affine in that SAME K for the
Gram determinant to stay constant.

## 4. Primary-source screen and its limits

The following sources supply CONTEXT, not adopted HC4 inputs.

* Hu–Sheng, [arXiv:2610.28375v1](https://arxiv.org/html/2610.28375v1),
  Section1.3 Theorem1.1, full Section7 proof, and AppendixB PropositionB.1
  and its full proof were read. Their weak quadratic hypothesis is
  c*rho/psi^alpha <= u <= C*rho*psi^(1-alpha), where rho=e+|z|^2,
  psi=1+log(rho), 0<=alpha<=1. The proof invokes their prior estimates;
  the interior regularity/Calabi foundations were NOT independently
  reproduced. AppendixB gives a smooth positive nonquadratic solution
  with quadratic upper growth, demonstrating why a lower condition
  matters. No such lower condition is established here.
* Chen–Hu–Sheng, [arXiv:2609.18265v1](https://arxiv.org/html/2609.18265v1),
  introduction/main statements and full Section4 proof were read. Their
  Liouville result assumes REAL convexity, also allowing convexity outside
  a compact set. The proof uses their global Hessian bound and
  Evans–Krylov regularity; those foundations were NOT independently
  reproduced. The shears in Section2 fail the required convexity.
* Pan–Zhang, [arXiv:2504.17625v1](https://arxiv.org/html/2504.17625v1),
  introduction and Section4 Theorems4.1/4.2 and their complete displayed
  proofs were read. Their entire nonquadratic flat examples and restricted
  ansatz do not establish completeness for polynomial Keller pullbacks.
  Theorem4.1 refers back to earlier flatness computations; Osserman's
  theorem used in Theorem4.2 was NOT independently reproduced. Neither
  theorem is imported here. Broader Calabi questions are not adjudicated.

No full-paper or foundational reproof, independent replay, formal
verification, or external review is claimed. No primary files cached.
The shared PF2RULED130 ledger scope was screened for changed global
geometry; its checker was NOT run and its proof was not newly imported
into this audit.

## 5. Remaining gate

FULL THEOREM TARGET: prove or break full HC4/JC2.
CURRENT LAST UNSOLVED IMPLICATION remains

    ACTUAL polynomial plane source
       ==> Tr(Lambda_L*Lambda_L) subset B.                 OPEN

Full marked-energy integrality remains OPEN. Positive Kähler flatness
is automatic, completeness retains inversion, coercive growth already
implies properness, and finite polynomial-square enlargement retains a
square Keller map. Revisit this route only with a new global polynomial
input proving completeness/properness, or a different metric construction
with its global hypotheses established. No unchecked SOS factorization,
automatic convexity/growth, or Ricci-tensor substitution is admissible.
Published generic-degree-four PLANE exclusion is retained. Review and
outreach remain DEFERRED.
