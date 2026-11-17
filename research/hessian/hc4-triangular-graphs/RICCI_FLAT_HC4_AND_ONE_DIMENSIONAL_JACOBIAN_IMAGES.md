# Ricci-flat HC4 and one-dimensional Jacobian images

HC4MYMJ0. 16 November 2026. ALL complex polynomial degrees.

**Theorem.** If a four-variable polynomial has nonzero constant Hessian
determinant and its complex bilinear Hessian metric has zero LEVI-CIVITA
RICCI TENSOR, its gradient is a polynomial automorphism. In particular this
holds if the metric is flat. The universal four-variable conclusion USES
the retained Tang theorem HCTANG2. The plane lemma below is elementary
and does not use that classification.

This settles a strict subclass, not HC4. Constant Hessian determinant does
NOT make this Ricci tensor zero. Scalar curvature alone is insufficient.
This is neither a Hermitian metric nor the affine/Kahler Ricci expression
-partial_i partial_j log(det H), which is automatically zero here.
No automatic full trace-pairing or marked-energy implication is proved.

## 1. A global plane lemma, with an explicit inverse

Let K=(P,Q):A2 -> A2 be complex polynomial with det A=j in C*, where
A=DK. Regard the FOUR entries of A as a polynomial map A:A2 -> A4.
Suppose its differential has rank at most ONE over C(x,y). Then

    K(z)=c+L*(z+v*f(ell(z))),                              (1)

where L is a constant invertible matrix, ell is a nonzero constant linear
form, v is a constant vector with ell(v)=0, and f is a polynomial. The
constant-matrix case is included by f=0. Thus K is affine-equivalent to
one shear, with explicit inverse

    z=L^-1*(K-c)-v*f(ell(L^-1*(K-c))).                     (2)

Here is the GLOBAL proof; no generic local parameter is promoted to a
polynomial source coordinate without justification.

If A is constant the result is immediate. Otherwise let Z be the reduced
closure of its image. The characteristic-zero rank/dimension theorem makes
Z an irreducible affine curve. Restrict A to an affine source line on which
it is nonconstant. This gives a nonconstant polynomial parametrization
gamma:A1 -> Z. It is finite onto its image: one nonconstant coordinate
gamma_i(t) makes t integral over C[gamma_i], and the same finite generators
work over the intermediate image coordinate ring. Hence the image is closed
and is ALL of Z.

The normalization Z^nu is A1. Indeed its projective completion is dominated
by P1, so has genus zero; every omitted point has all its preimages at the
single source infinity. There is exactly one omitted point. This is the
elementary curve argument already retained in MH0, before its separate
etale degree-one step. No etale or injectivity property is needed for gamma.

Write O(Z^nu)=C[t]. Its generator t is in Frac(O(Z)) and integral over
O(Z). Pulling it into C(x,y), it is integral over C[x,y], so NORMALITY of
C[x,y] makes it a polynomial t(x,y). Consequently

    A=B(t(x,y)),          B(T) in Mat_2(C[T]).              (3)

No claim that t itself is a source coordinate is made. This step retains
the full normalization, rather than a rational parametrization with poles.

Integrability of each row of A gives

    b_i1'(t)*t_y=b_i2'(t)*t_x.                             (4)

Choose a row for which (b_i1',b_i2') is not identically zero. For generic
c, that vector at c is nonzero and the fibre t=c is reduced and smooth.
There are only finitely many exceptional critical values: t is constant
on each irreducible component of its critical locus. Every fibre is
nonempty, since a nonconstant restriction to a suitable source line is a
surjective one-variable polynomial.

On every irreducible component of a generic fibre, (4) makes the linear
function b_i1'(c)*x+b_i2'(c)*y have zero differential. In characteristic
zero it is constant there. Thus that closed one-dimensional component
is an ENTIRE affine line. Generic fibres can have several components;
neither connectedness nor geometric irreducibility is assumed.

Lines in different fibres are disjoint. Two nonparallel complex affine
lines intersect, so the normals (b_i1'(c),b_i2'(c)) have the SAME projective
direction for all generic c. This is a constant direction, not a moving
normal form on different fibres. Returning to (4) gives a nonzero constant
linear form ell with t=F(ell) for a polynomial F. Therefore every entry
of A is a polynomial in the SINGLE LINEAR SOURCE FORM ell.

Now write A=C(ell). Row integrability says C_i'(T) is a multiple of ell.
Integrating in T gives

    A=L+G'(ell)*ell,       G(0)=G'(0)=0,
    K(z)=c+L*z+G(ell(z)),  L=C(0), det L=j.                (5)

Here ell is a row and G' a column. The determinant lemma and det A=j give

    ell*L^-1*G'(T)=0,
    ell*L^-1*G(T)=0.                                      (6)

Thus L^-1 G(T)=v*f(T) in the one-dimensional kernel of ell. This proves
(1) and (2), with no degree bound and no line-injectivity theorem invoked.

## 2. The actual Levi-Civita Ricci tensor

Let H=Hess(h) be nondegenerate with constant determinant, in any dimension,
and put C_i=H^-1 partial_i H. Use curvature convention

    R_ij=partial_i Gamma_j-partial_j Gamma_i
                       +[Gamma_i,Gamma_j],
    Ric_(bj)=sum_a (R_aj)^a_b.

For a Hessian metric Gamma_i=C_i/2. Direct contraction gives

    Ric_(ij)=1/4*tr(C_i*C_j).                              (7)

For clarity, no positivity argument is hidden here. The coordinate Ricci
formula, tr Gamma_i=0 and symmetry of h_(ijk) reduce it to

    1/2*sum_(k,l) H^kl*h_(lkij)-1/4*tr(C_i*C_j).

The derivative term in partial_k H^kl vanishes: by symmetry it is the
contraction of H^-1 with a derivative of log(det H). Differentiating
tr(H^-1 partial_i H)=0 once more identifies the remaining fourth-derivative
contraction with tr(C_i C_j), proving (7). Equivalently it follows by
contracting the standard Hessian curvature identity.

Primary comparison: Amari--Armstrong,
[*Curvature of Hessian manifolds*, arXiv:1312.1103v1](https://arxiv.org/html/1312.1103v1),
Section3, Proposition1 and its complete proof were read. Its algebraic
curvature identity supplies context; (7) is derived here directly for the
complex bilinear metric. No positive-definite, compactness, characteristic-
class or completeness theorem is imported, and no full-paper independent
reproduction or novelty priority is asserted.

## 3. The cotangent block and the three-dimensional trace form

Let h(z,w)=D(z)+w dot K(z), with z=(x,y) and det DK in C*. Set A=DK and
L_i=A^-1 partial_i A, for i=x,y. The Hessian blocks and inverse from MI0
apply to ANY plane Keller pair, not just its two-shear controls. Multiplying
these blocks shows

    tr(C_i*C_j)=2*tr(L_i*L_j),             i,j in{x,y},
    tr(C_i*C_j)=0 if i or j is a w index.                  (8)

All D and w terms are in strictly lower off-diagonal blocks and do not
contribute. Thus the ONLY possibly nonzero Ricci block is

    Ric_(ij)=1/2*tr(L_i*L_j),             i,j in{x,y}.       (9)

In particular the scalar curvature is ALWAYS ZERO for EVERY cotangent
plane Keller potential: the upper-left block of H^-1 is zero. This scalar
identity does not imply that the Ricci block (9) vanishes.

Since det A is constant, L_x,L_y are traceless 2-by-2 matrices over the
field F=C(x,y). The trace form (U,V)->tr(UV) is nondegenerate on sl_2(F),
which has dimension THREE. If Ric=0 then

    tr(L_x^2)=tr(L_x L_y)=tr(L_y^2)=0.                     (10)

The span W of L_x,L_y is totally isotropic. Nondegeneracy gives
W subset W^perp and dim W+dim W^perp=3, so dim W<=1. This argument is
valid over F itself; algebraic closure and positive definiteness are not
assumed. Hence A_x and A_y are linearly dependent and the polynomial map
defined by the entries of A has differential rank at most one.

Section1 now makes K an explicit affine shear automorphism. The cotangent
gradient inverse follows by recovering z from K(z) and then

    w=A(z)^(-t)*(xi-grad D(z)).                            (11)

This proves the Ricci-flat cotangent result in EVERY degree, for arbitrary
D and every plane Keller pair meeting (10).

## 4. Universal HC4 consequence, using Tang

Apply HCTANG2 to an arbitrary complex four-variable unit-Hessian h. If its
gradient is already an automorphism there is nothing to prove. Otherwise
Tang supplies its cotangent form after a CONSTANT invertible linear source
change. Hessian metrics under this change are actual pullback metrics;
vanishing of their Levi-Civita Ricci tensors is invariant. Section3 then
makes the plane pair, and hence the full gradient, invertible, contradicting
the nonautomorphic alternative. This proves the theorem stated above.

The use of Tang is essential to this universal claim. Without it Sections1
and3 concern the stated plane/cotangent classes. No arbitrary-dimensional
HC theorem or classification-free four-variable conclusion is asserted.

## 5. Sharp distinctions and retained hostile controls

MI0 already supplies polynomial gradient automorphisms for which (9) is
NONZERO. For its a=b=1,D=0 example, at the source origin

    L_x=2*E_21, L_y=2*E_12,
    Ric on the base = [[0,2],[2,0]].                       (12)

The scalar curvature is nevertheless zero. Therefore neither constant
Hessian determinant nor scalar flatness supplies the hypothesis here.
Ricci flatness is a sufficient condition, not a necessary condition for
gradient inversion or a renamed equivalence with HC4.

Ricci flatness is also weaker than full metric flatness within the actual
polynomial class. Take

    K=(x,y+x^2),       D=y^3,
    h=y^3+r*x+s*(y+x^2).                                  (13)

It has determinant1 and an explicit gradient inverse. Its L_x is a multiple
of E_21 and L_y=0, so its ENTIRE Ricci tensor is zero. But at the origin
the curvature operator R_xy has entry (row r,column y)=-3 with the above
convention. Thus the metric is curved. No flatness or completeness was
silently inferred from the Ricci condition.

The rank-one image lemma also requires the Keller determinant. For example
K=(x^2,y) has a Jacobian matrix varying in one direction and is not an
automorphism. The normality and global polynomial-image steps in Section1
must not be replaced by a rational local factorization: the original
source and the entire image curve are essential.

## 6. Verification and remaining scope

The checker supplies exact tensor, matrix, shear inverse, full coordinate
Ricci and hostile-control transcriptions. The normalization, fibre-line,
and trace-form dimension arguments are WRITTEN ALL-DEGREE PROOFS, not
extrapolations from those controls. Prior registered checkers are not run.

This result excludes every hypothetical HC4 counterexample with zero
Levi-Civita Ricci tensor. It does not prove that tensor zero automatically;
(12) explicitly forbids that unrestricted inference. A remaining
nonautomorphic Tang cotangent potential would be scalar-flat with a
nonzero base Ricci block and a genuinely two-dimensional Jacobian-matrix
image. No such polynomial Keller source is constructed here.

FULL THEOREM TARGET: prove or break full HC4/JC2.
CURRENT LAST UNSOLVED IMPLICATION remains

    ACTUAL polynomial plane source
       ==> Tr(Lambda_L*Lambda_L) subset B.                 OPEN

Full marked-energy integrality is also OPEN. No independent replay, formal
verification or external review is claimed. Review and outreach remain
DEFERRED. Published generic-degree-four PLANE exclusion is retained; no
intermediate is assumed to be a plane and no missing artifact is rebuilt.
