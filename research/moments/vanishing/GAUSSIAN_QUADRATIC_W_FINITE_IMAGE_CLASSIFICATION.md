# Complete finite images in the quadratic-W, constant-a architecture

Let K have characteristic zero. Use the normalized Gaussian with
E[ZW]=E[Y^2]=1 and all other basic pairings zero. Consider

    P=ell(Z)W^2+W(h(Z)+c(Z)Y)+aY^2+u(Z)Y+k(Z),      (C1)

where a is constant and ell,h,c,u,k are arbitrary polynomials over K.
Write Phi_P[Q]=E[Q exp(tP)], V_P=Phi_P[K[Z,W,Y]] and M_P=K[t]V_P.
Safety retains EVERY original polynomial marker; purity means Phi_P[1]=1.

## 1. The exhaustive finite-image table

The following cases are disjoint and exhaustive. The stated condition is
necessary and sufficient for the FULL actual image M_P to be finite over
K[t]. No coefficient degree or root multiplicity is otherwise restricted.

| Case | Necessary and sufficient condition | Primitive exponent psi |
|---|---|---|
| c nonconstant, any ell | a=0, ell in Z^3K[Z], h,c in Z^2K[Z] | k(0)t+u(0)^2t^2/2 |
| c constant, deg ell>=2 | c=a=0, ell in Z^3K[Z], h in Z^2K[Z] | k(0)t+u(0)^2t^2/2 |
| c constant, deg ell=1 | a=0 and h,u,k constant | k(0)t+u(0)^2t^2/2 |
| ell,c constant, (ell,c)!=(0,0) | P quadratic and B=Gamma A nilpotent, in (C4) | Formula (C5) |
| ell=c=0 | a=0 and either h is nonzero constant or h in Z^2K[Z] | Formula (C6) |

In EVERY finite case,

    Phi_P[1]=exp(psi), M_P=exp(psi)K[t].             (C2)

Thus every finite original vacuum in this entire architecture is
primitive. For nonconstant P one has literally V_P=M_P; for the constant
phase P=k0, the necessary exception is V_P=K exp(k0t), while
M_P=exp(k0t)K[t]. The zero phase is included in this exception.

Safety is precisely the appropriate table condition together with

    psi=0.                                         (C3)

Consequently EVERY safe member of (C1) is pure. This is a scoped
safety-to-purity result. It neither proves that implication for arbitrary
Gaussian phases nor classifies every nonfinite image or its generic rank.
Purity alone is not substituted for the whole-image premise in rows1-2.

The proofs and their exact assurance boundaries are:

* Row1: the [unrestricted mixed theorem](GAUSSIAN_QUADRATIC_W_MIXED_FINITE_IMAGES.md).
* Row2: the [arbitrary-c variable-square theorem](GAUSSIAN_VARIABLE_W_SQUARE_FINITE_IMAGES.md).
  A nonzero constant c cannot belong to Z^2K[Z], so only c=0 survives.
* Row3: the [affine-square theorem](GAUSSIAN_AFFINE_W_SQUARE_FINITE_IMAGES.md).
* Row4: the [constant nonzero W-square theorem](GAUSSIAN_CONSTANT_W_QUADRATIC_FINITE_IMAGES.md)
  when ell!=0, and the [constant nonzero mixed-source theorem](GAUSSIAN_W_AFFINE_FINITE_IMAGE_OBSTRUCTIONS.md#6-a-companion-constant-c-primitive-vacuum-result)
  when ell=0,c!=0.
* Row5: the [zero-mixed-source criterion](GAUSSIAN_W_AFFINE_FINITE_IMAGE_OBSTRUCTIONS.md#8-complete-finite-image-criterion-when-c0)
  and the direct scalar argument in Section3 below.

This table combines complete internal written proofs and separately
audited imports. It adds no claim of external review or formal verification.

## 2. Exact quadratic and zero-mixed-source exponents

In row4 write X=(Z,W,Y)^T and

    P=X^T A X+l^T X+k0,
    Gamma=[[0,1,0],[1,0,0],[0,0,1]], B=Gamma A.       (C4)

Thus h=h0+h1Z, u=u0+u1Z, k=k0+k1Z+k2Z^2, and

    A=[[k2,h1/2,u1/2],[h1/2,ell,c/2],[u1/2,c/2,a]],
    l=(k1,h0,u0)^T.

If B is nilpotent, B^3=0 and the primitive exponent is

    psi=k0t+(t^2/2)l^T(I-2tB)^(-1)Gamma l.          (C5)

The inverse is a polynomial matrix of degree at most two, so deg psi<=4.
The safety condition is k0=0 and l^T B^j Gamma l=0 for j=0,1,2.

In row5 the exponents are

    h=h0 nonzero constant:
      psi=t k(h0t)+t^2u(h0t)^2/2;
    h in Z^2K[Z]:
      psi=t k(0)+t^2u(0)^2/2.                       (C6)

For nonzero constant h0, safety is equivalently the exact polynomial identity

    k(Z)=-Z u(Z)^2/(2h0).                           (C7)

The original degree filtration is retained. In the SAFE shear subbranch
(C7), let
p=deg u if u!=0. The known Gaussian-preserving triangular shear gives
cutoff C=2p+1, attained by W^D for each D when u!=0; when u=0 the cutoff
is C=1, attained by Z^D. The proof and exact monomial costs are in
[the original shear formula](GAUSSIAN_AFFINE_ISOTROPIC_CHART.md#safety-and-cutoffs-in-the-first-two-branches).

## 3. A scalar criterion for every W-affine phase with c=0

The row5 finite-image condition follows already from the ACTUAL scalar
identity Phi_P[1]=p(t)exp(psi(t)), with p polynomial, p(0)=1 and
psi in tK[t]. This strengthens the formerly unstated c=0 scalar case;
constant a is essential.

First work over C. Let zeta be the selected formal root zeta=t h(zeta),
and put

    J=1-t h'(zeta), d=1-2at,
    R=t k(zeta)+t^2u(zeta)^2/(2d).

Original circular Lagrange contraction followed by ordinary Gaussian
completion gives exactly

    Phi_P[1]=exp(R)/(J sqrt(d)).                    (C8)

The quantities zeta,J,R,d are algebraic over C(t), with R(0)=0,J(0)=1.
Squaring the assumed scalar identity gives

    exp(2(R-psi))=p^2 J^2 d.                        (C9)

An algebraic function with algebraic exponential is constant. On a
smooth projective curve carrying both functions, a pole of the first
would give a pole of order at least two in its differential, while the
logarithmic differential of the second has at most simple poles. Thus
there are no poles and the first function is constant. The t=0 value
in (C9) consequently yields

    R=psi, p^2J^2d=1.                              (C10)

If h(0)=0, the selected root is identically zero. Then J=1-t h'(0), so
the polynomial product (C10) forces a=h'(0)=0,p=1. This is exactly the
h in Z^2K[Z] case.

Suppose h(0)!=0 and deg h=D>=2. The root curve z-t h(z)=0 is irreducible
because gcd(z,h)=1 and the polynomial is primitive and linear in t.
Its rational normalization is t=z/h(z), and the selected local branch
gives an injective embedding into C((t)). Hence (C10) is an identity
on that entire curve. There

    J=(h-z h')/h.

Set N=h-z h'. It has degree D. At a root alpha of h of multiplicity m,
alpha!=0 and N has EXACT multiplicity m-1: the term -z h' has that
order while h has order m. Therefore N/gcd(h,h') has positive degree
equal to the number of distinct h-roots and has no root in common with h.
Choose a root beta of this quotient. Then h(beta)!=0, t_beta=beta/h(beta)
is finite, and J(beta)=0. Both p(t_beta) and d(t_beta) are finite,
contradicting (C10). This eliminates every D>=2 without a moment cutoff.

If h=h0+h1Z with h0!=0, then zeta=h0t/(1-h1t) and J=1-h1t.
The same polynomial product forces a=h1=0,p=1. Thus only a nonzero
constant h survives, with exponent (C6).

Conversely either row5 condition gives the displayed primitive scalar
and the whole original finite image. For h in Z^2 use positive angular
balance after translating Y. For h=h0!=0 and a=0, the exact auxiliary
source is

    E exp(tP+sZ+vW+qY)
      =exp(s(h0t+v)+q^2/2+t k(h0t+v)
                 +t q u(h0t+v)+t^2u(h0t+v)^2/2).   (C11)

Every fixed s,v,q coefficient is exp(psi) times a polynomial in t.
Thus all ORIGINAL marked periods are controlled. The nonconstant-phase
original-image theorem gives V=M, with the constant exception already
stated after (C2). Finite-coefficient embedding into C descends necessity
to any characteristic-zero field; these source formulas prove sufficiency
algebraically.

For spatially variable a(Z), d(t,zeta) can have poles when the selected
root escapes to infinity at finite time and cancel J^2. No such variable-a
scalar theorem is asserted here.

## 4. A sharper cutoff for safe quadratic phases in three variables

This argument uses the exact quadratic Gaussian formula and does not
assume general safety implies purity. It applies to any nonsingular
symmetric three-dimensional Gaussian covariance Gamma: no special
entry of the displayed Witt matrix is used. If P is quadratic and safe, its
unit f=p is polynomial. The quadratic source gives

    f=det(I-2tB)^(-1/2) exp(R),
    R=k0t+(t^2/2)l^T(I-2tB)^(-1)Gamma l.

Rational-exponential rigidity applied to exp(2R)=det(I-2tB)p^2 implies
R=0, det(I-2tB)=1 and p=1. Hence B^3=0. Put v=Gamma l and use the
nondegenerate bilinear form with matrix Gamma^(-1). The operator B is
self-adjoint. Since R=0,

    <v,B^jv>=0, j=0,1,2.

The cyclic span U=span(v,Bv,B^2v) is totally isotropic, because
<B^iv,B^jv>=<v,B^(i+j)v>, and B^3=0 handles all larger sums.
Thus U is contained in its perpendicular and 2dim U<=3. It has dimension
at most one, and nilpotence then forces Bv=0.

The tilted mean is therefore t v, of degree one, while its covariance
(I-2tB)^(-1)Gamma has degree at most two. A Wick term in a marker of
degree D with r pairs and D-2r single means has time degree at most
2r+(D-2r)=D. Consequently

    E[Q P^m]=0 whenever m>deg Q.                    (C12)

Every safe row1-3 member and every safe row5 zero-root member has the
same C=1 by angular balance. Thus the nonlinear constant-h shear (C7)
is the only table branch needing a larger uniform original cutoff.
The retained mixed-cutoff theorem gives r-jet budget (Cd+1)r for
degree-at-most-d perturbations at each of these safe phases, and the
corresponding certified nilpotence tail budget (Cd+1)(nu-1).

## 5. Exact safe recompositions in the entire architecture

The [general safe-recomposition theorem](GAUSSIAN_SAFE_RECOMPOSITION_CLASSIFICATION.md)
now becomes an explicit test in every row. Let S_P={f in K[T]:f(P) is safe}.
If the relevant table condition fails, S_P=(0). If it holds, compute
the primitive exponent psi from the table. Then

    S_P=(T-lambda) if psi(t)=lambda t;
    S_P=(0) otherwise.                             (C13)

Indeed every period in a finite case is exp(psi) times a polynomial,
so a safe shift exists exactly when psi is linear, and its shift is
lambda. The general recomposition theorem handles every nonconstant
polynomial f, even when f(P) has higher W degree than (C1). Constant
phases have the same formula directly: f(P)=0 is the condition f(k0)=0.
No single common bad marker is asserted for all rejected recompositions.

In the constant-h row5 branch, a safe shift exists precisely when
k(Z)+Z u(Z)^2/(2h0) is a constant lambda. In rows1-3 and the zero-root
row5 branch, it requires u(0)=0 and lambda=k(0). For row4 it requires
l^T B^j Gamma l=0, j=0,1,2, and lambda=k0.

## Evidence and remaining scope

The table and its all-original-marker primitivity conclusion use the
complete written branch proofs. The new scalar row5 argument is an
elementary ramification/Jacobian obstruction, independent of the marked
curve-volume estimates. The new quadratic cutoff follows from the
three-dimensional isotropic-span bound.

The branch checkers retain their precise bounded algebraic scopes;
their passing does not certify unrestricted comparison theory. Those
imports retain their current primary-source audits and pending
independent-review status. No external review, formal verification,
publication or priority claim is made here. General Gaussian
safety-to-purity, arbitrary nonfinite image ranks and higher W/fibre
architectures remain open.
