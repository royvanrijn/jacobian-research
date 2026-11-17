# Gaussian moment-preserving Keller endomorphisms and the inverse boundary

15 November 2026. Internal written deductions, with exact auxiliary certificates. Independent replay, external review and formal verification of these new claims remain pending. No novelty claim is made. Existing foundational/cotangent maps and moment-map identities are credited explicitly below.

The [previous automorphism theorem](POLYNOMIAL_GAUSSIAN_SYMMETRIES_AND_INVARIANT_TIMES.md) assumes a polynomial inverse. Its sufficient half needs less. This note proves that extension, gives a universal quadratic-preserving lift of any normalized Keller map, and applies the Gaussian criterion to the already known weighted six-dimensional cotangent collision. The extension preserves full polynomial moment data while leaving a real distinction between automorphisms and endomorphisms.

## 1. Sufficient Gaussian criterion without an inverse

Let A be any commutative Q-algebra, B an invertible symmetric covariance matrix, S=X^T B^(-1)X/2, and L_B the Wick functional. Let Phi be a polynomial map whose Jacobian determinant j is a UNIT in A[X]. No polynomial inverse is assumed. Put

    ell=log(j/j(0)),    C=S composed Phi-S-ell.

As in the previous theorem, ell is a finite polynomial on the nilpotent nonconstant coefficients of the unit j. If C is constant in the coordinates, then

    L_B(F composed Phi)=L_B(F) for every polynomial F.         (1)

Proof. Let D_i be the derivations given by the columns of (Jac Phi)^(-1), which is polynomial because j is a unit. They satisfy D_i(phi_k)=delta_ik whether or not Phi is globally invertible. Piola's cofactor identity and Gaussian integration by parts give

    L_B(D_i G)=L_B((D_i S+D_i ell)G).

Constancy of C makes D_i S+D_i ell=(B^(-1)Phi)_i. Taking G=F composed Phi shows that F->L_B(F composed Phi) satisfies the normalized Gaussian IBP relations. Their uniqueness proves (1). No surjectivity of the pullback algebra, global inverse or analytic change of variables is used.

Thus a constant-unit-Jacobian map preserving S is sufficient in every degree and dimension over every Q-algebra. Over fields this includes quadratic-preserving Keller maps. The proof is only sufficient for endomorphisms: the necessity proof in the previous theorem used a polynomial inverse to test against EVERY source marker. We do not assert a complete logarithmic classification without that hypothesis.

The pullback homomorphism Phi*:A[Y]->A[X] is injective under (1). If Phi*F=0, then L_B(FG)=L_B(Phi*F Phi*G)=0 for every G; the separating Hermite pairing from the previous proof gives F=0. This is injectivity of an algebra homomorphism, not injectivity of the point map. A collision below makes the distinction concrete.

## 2. A general all-degree quadratic-preserving Keller suspension

Let F:A^n->A^n be polynomial with F(0)=0, JF(0)=I and det JF=1. Any constant-unit-Jacobian Keller map can be put in this form by translation and invertible linear output normalization. Set

    F_t(x)=t^(-1)F(tx),   A_t=J_x F_t,
    Y=A_t^(-T)y,
    H_F(x,y,t,z)=(F_t(x),Y,t,z+(partial_t F_t)^T Y).           (2)

All expressions are polynomials. The constant term of F is zero, so F_t has no negative powers of t; det A_t=1 and A_t^(-1)=adj(A_t). No inverse of F is needed.

The map H_F in dimension2n+2 preserves

    S=x dot y+t z,

has Jacobian determinant1, and hence preserves every polynomial Gaussian moment for the split covariance with pairings L(x_i y_j)=delta_ij and L(tz)=1.

To prove quadratic preservation, each degree-d part of F contributes t^(d-1)F_d(x). Weighted Euler gives

    A_t x=F_t+t partial_t F_t.

Pair this with Y, using A_t^T Y=y, and obtain

    F_t dot Y+t(z+(partial_t F_t)^T Y)=x dot y+t z.

For the determinant order coordinates as (x,y,t,z). After the row/column for the fixed t is removed, the final z coordinate has coefficient1, while the x,y block is triangular with diagonal blocks A_t,A_t^(-T). Their determinants multiply to1. The t-derivatives introduce no extra determinant factor.

The construction is the signed cotangent lift of the scaling family (x,t)->(F_t(x),t). It preserves the one-form y dot dx-z dt exactly. Cotangent lifting is standard prior work, retained in [SW1](../../archive/non-elliptic/extended-geometry/SYMPLECTIC_WEYL_LIFT.md); the weighted Euler argument explains this particular Gaussian quadratic. No new cotangent-lift mechanism is claimed.

H_F is a polynomial automorphism if and only if F is. If F has a polynomial inverse, the inverse family t^(-1)F^(-1)(tX) is polynomial, and its cotangent lift inverts (2). Conversely H_F preserves t, and target y=0,z=0 forces source y=0,z=0, because A_t^T is invertible. Restricting a hypothetical polynomial inverse to this zero-momentum section at t=1 yields a two-sided polynomial inverse of F.

Every collision F(a)=F(b) lifts to

    H_F(a,0,1,0)=H_F(b,0,1,0).

Over characteristic-zero fields the generic degree of H_F equals that of F: on t!=0 the position map F_t is linearly conjugate to F, and the momentum/last-coordinate transformations are invertible linear/triangular maps over the position function field. No statement about a special fibre degree being constant is made. In fact at t=0 the base F_t is identity, while H_F restricts to a triangular automorphism of the fibre; this does not imply generic inversion.

## 3. Weighted cotangent maps and the existing six-dimensional collision

There is a smaller route when F commutes with a nonzero-weight Euler field. Let W=diag(w_1,...,w_n) with every w_i a unit in A and suppose

    JF(x) W x=W F(x).                                       (3)

For any Keller F with constant-unit Jacobian, its usual cotangent map

    K_F(x,p)=(F(x),(JF(x))^(-T)p)                            (4)

has determinant1 and preserves S_W=x^T W p. Indeed (3) implies

    F^T W (JF)^(-T)p=x^T W p.

Here symmetry/commutativity of the diagonal W is used. The Hessian matrix of S_W is invertible, with inverse covariance

    B_W=[[0,W^(-1)],[W^(-1),0]].

Section1 therefore gives preservation of the ENTIRE Wick functional, not just the quadratic moment. The collision, generic degree and noninvertibility in (4) are precisely the standard cotangent facts already retained in SW1.

For the foundational F1 map write u=1+xy and normalize its outputs to

    f1=x-3x^2y/2-x^3z/2,
    f2=y+3xu^2z+3xy^2(4+3xy),
    f3=u^3z+y^2u(4+3xy).                                   (5)

This is (F1_3/2,F1_2,F1_1), so Jf(0)=I and det Jf=1. It obeys (3) for weights (1,-1,-2), directly monomial by monomial. The three rational source points

    a0=(0,0,-1/4), a+=(1,-3/2,13/2), a-=(-1,3/2,13/2)

have common image (0,0,-1/4). Thus K_f on six variables has the three distinct zero-momentum preimages (a0,0),(a+,0),(a-,0), preserves S_W, has determinant1, and preserves every polynomial Gaussian moment for B_W. It has no polynomial inverse.

The map F1, its collision and its cotangent lift are existing results, not new counterexamples. Source credit: Alpoge's foundational map is retained as F1; the cotangent lift is credited to W.G.P.Mayner in the contemporaneous source and SW1; the exact weighted moment-map identity is credited to Daniel Kirtchakov's TheoremC(i), at the [pinned primary source](https://github.com/05oz/certify/blob/dd607119e6a198df39bf0b9960f46d240a5c3404/paper/preprint-dixmier-poisson.tex). That source uses different source/target weight orderings. The normalization (5) makes the same quadratic appear on both sides. We rederive the weight and inverse-Jacobian identities. The Gaussian functional implication has the separate proof in Section1. No claim of priority, independent theorem replay or new formal verification is made.

To express this as ordinary variance-one Gaussian moments, use the fixed invertible complex linear change for each pair

    x_i=U_i+iV_i,   p_i=(U_i-iV_i)/(2w_i).

Then S_W=sum_i(U_i^2+V_i^2)/2 and B_W is the corresponding covariance. Conjugating K_f by this change gives a complex polynomial map of six variables with constant Jacobian1 and an explicit three-point complex collision, preserving every ordinary Gaussian polynomial moment. This is a holomorphic polynomial-moment assertion for a complex-coefficient map; it is not a real Gaussian probability-space bijection or an assertion about moments involving independent conjugate output functions.

## 4. What transfers without an inverse

For ANY moment-preserving polynomial map Phi, no Jacobian assumption is needed for these identities:

    L((P composed Phi)^m)=L(P^m),
    L((Q composed Phi)(P composed Phi)^m)=L(QP^m).            (6)

Consequently:

* all pure sequences and their literal ordinary coefficient ideals are preserved in a shared coefficient family;
* every fixed failure witness Q pulls back to Q composed Phi, retaining every marked moment;
* if P composed Phi is eventually safe against all SOURCE markers, then P is eventually safe against all TARGET markers;
* a linear cutoff C for P composed Phi gives coefficient C deg Phi for P.

The last point uses deg(Q composed Phi)<=deg Phi*deg Q. Conversely, no proof here transfers target safety to source safety: a general source marker need not be in the pullback algebra, and a polynomial inverse is absent. The preceding automorphism theorem retains its exact two-way result.

In particular the field first94 pure test for a ternary cubic transports to its pullback by a moment-preserving map (after harmless target padding if needed). This yields pure-null structured unbounded-degree families. Their full mixed safety is NOT thereby proved. Likewise the harmonic-core theorem cannot be transported as a full safe class merely from (1).

The remaining structural question is registered as OP-GMC-ENDOMORPHISM-SAFETY: for a Gaussian moment-preserving Keller endomorphism over C and a pure-null safe target phase P, must P composed Phi remain safe against every polynomial source marker? Automorphisms are closed positively by the previous theorem. Equivariant strict angular-separation cases remain positive by their own weight argument. This note gives no counterexample to or general proof of the implication. It also does not assert that the sufficient logarithmic criterion is necessary for every moment-preserving endomorphism.

The normalized homogeneous construction (2) is equivariant for weights x:+1,y:-1,t:-1,z:+1. For example, its individual pulled-back isotropic coordinate phases remain strictly one-sided and safe. A noninvertible state-preserving map is therefore not by itself a scalar GMC counterexample. The known lower-dimensional scalar GMC failures and all protected cutoffs remain unchanged.

## 5. Checkers and assurance boundary

The [checker](verify_gaussian_keller_endomorphisms.py) implements a sufficient logarithmic residual certificate without requesting an inverse, an explicit weighted cotangent constructor, and the scaling-family constructor (2). It verifies the foundational normalized determinant, Euler identities, inverse-Jacobian identities, moment-map equality and exact collision, plus selected finite moments. Determinant1 of the large lift follows from the written block identity; no expensive unneeded full eight-variable expansion is used. The six-dimensional lift is specified by its small base map and inverse-Jacobian certificate in the [receipt](gaussian_keller_endomorphism_verification.json).

The [source audit](gaussian_keller_endomorphism_source_audit.json) binds this note, checker, receipt, pinned primary source, F1/SW1 dependencies and retained evidence. The inherited formal assurances on F1 are preserved; they are not transferred to these new Gaussian deductions. General inverse recognition, a complete classification of moment-preserving endomorphisms, arbitrary safe-phase pullback and general moment solvers remain open.
