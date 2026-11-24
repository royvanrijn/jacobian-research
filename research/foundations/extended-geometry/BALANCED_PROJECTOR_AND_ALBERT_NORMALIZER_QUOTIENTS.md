# Balanced projector and Albert full-normalizer quotients

This proof gives explicit quotient phases on the balanced unitary and quaternionic Grassmann models and on compact and reduced Albert full-normalizer quotients. Its claims concern these actual rational models. The known endpoint and CDF mechanisms retain their attribution; no general arbitrary-phase descent theorem is used.

## 1. Inputs, exact model scope, and trace restriction

Use the actual rational models in [the classical-model proof](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md) (SHA256 `77e79160d3cc1478cabd9eea55fdc4c0374ac937312157a59fb90c99d047b9de`), Sections 2.2, 2.3, 3, 4, 5, 7, and 13. For positive integer k let n=2k. In AIII the group is SU_n(Q(i),identity Hermitian form), acting on rank-k Hermitian projectors P; in CII it is the quaternionic unitary group on n Hamilton-quaternion coordinates, acting on rank-k quaternionic Hermitian projectors. The stabilizers H are respectively S(U_k times U_k) and Sp(k) times Sp(k). The Q-varieties X are smooth, geometrically integral and affine. Their coordinate D=P11 has the exact beta polynomial functional with parameters (alpha,alpha), where alpha=k in AIII and alpha=2k in CII. These are positive integers, including k=1.

For the Albert assertions use [the Albert-model proof](ALBERT_EIV_AND_REDUCED_NORM_ONE_PRODUCT_COUNTEREXAMPLES.md) (SHA256 `297d5c3f2ab8439039644a83f423005c8bdc18fb8b71b64dc42884d0242ed1af`), Sections 1--8. This distinguishes the compact unitary Q-model X_c=G_c/H_c from an arbitrary reduced Albert norm-one variety Y_J over a characteristic-zero field. Its geometric model, full F4 stabilizer and exact beta recurrence are inputs, not conclusions inferred from a finite moment prefix.

If a finite subgroup of equivariant automorphisms acts freely on one of these affine X, its affine quotient is smooth and geometrically integral, and the map is finite etale. Equivalently the quotients below are homogeneous G/N with reductive finite extension N of H. Coordinate pullback is injective. The normalized invariant trace on the quotient equals the restriction of the upstairs trace, by invariance, normalization and uniqueness. Thus a phase that is ACTUALLY invariant has the same every-order identities downstairs. This proves no descent for arbitrary upstairs phases.

**Theorem 1 (balanced full-normalizer products).** For each k>=1, the stated balanced AIII and CII full-normalizer quotients have, over every characteristic-zero field, the pure seed p and fixed unit marker t of Section 4, with exact marked value c_m at every positive order, sharp monomial-marker equivalences, absolute unit-nonpower and primitive Laurent support. Over every formally real field both individual factors are Mathieu and their product is not.

## 2. Balanced projectors: actual full normalizer

Let P0=diag(I_k,0_k) and

    J=[[0,I_k],[-I_k,0]].

This is Q-rational in both groups, has J squared=-I in H, and sends P0 to I-P0. Its determinant is 1 in the AIII model. The right action of its coset on G/H is the equivariant involution

    sigma(P)=I-P.

There are no geometric fixed points: P=I-P would give P=I/2, which is not idempotent in characteristic zero. Consequently X/sigma is a genuine free finite-etale quotient; it is not a replacement by an unrelated projector scheme.

In AIII the connected center of H acts on the standard representation through its two distinct scalar weight spaces, the two k-dimensional blocks. A normalizer preserves that center and hence permutes these two spaces. The block-preserving part is H; the block-interchanging coset is JH. In CII the standard representation over an algebraic closure is the direct sum of the two nonisomorphic irreducible H-modules from the two symplectic factors. A normalizer again permutes these two isotypic spaces, and the block-preserving part is exactly H. This remains valid at k=1, since the two factors act on different modules. Therefore in both actual models

    N_G(H)=H union JH,   N_G(H)/H=constant Z/2,
    G/N_G(H)=X/{P maps to I-P}.

The geometric calculation and the Q-rational representative J establish the statement over Q, not merely over the real points. Centers of the acting SU or quaternionic unitary groups already act trivially on projectors, so passing to a central quotient of the acting group changes neither this actual model nor these functions.

## 3. Symmetric CDF and the invariant radius

For a positive integer alpha define the rational polynomial

    v(u)=B(alpha,alpha)^(-1) integral_0^u x^(alpha-1)(1-x)^(alpha-1) dx.

The integrand is symmetric under x maps to 1-x, and its integral over [0,1] is 1 after normalization. Polynomial differentiation and evaluation at 0 give the exact polynomial identity

    v(1-u)=1-v(u).

The beta functional and polynomial chain rule give lambda(v(D)^j)=1/(j+1) for every j>=0. Put

    S=(2v(D)-1)^2.

Then S is a Q-regular sigma-invariant function, so belongs to the ACTUAL quotient coordinate algebra. For every polynomial g,

    lambda(g(S))=integral_0^1 g((2u-1)^2) du
                =integral_0^1 g(z^2) dz.

The second equality follows by splitting at 1/2 and replacing z by its negative on one half. Over a general characteristic-zero field these equalities mean rational polynomial integration, followed by coefficient extension; no order or measure is assigned to that field. Thus every invariant S-moment is exactly 1/(2j+1) at every j.

An actual rational projector with D=1/2 exists. Take a rank-one projection onto e1+e_(k+1) in the plane spanned by these two vectors, using its matrix with all four entries 1/2, and add the projections onto e2,...,ek. This construction is Hermitian over Q in both models and has rank k. At its image in the quotient v(D)=1/2 and S=0.

## 4. A seed that is already pure and absolute nonpower

On either full-normalizer quotient times split Gm put

    T=1+t,
    p=T/t [1-T^2 S]^2,   H_marker=t.

This is an explicit rational regular pair. For k>=1 define the strictly positive rational

    c_k=integral_0^1(1-z^2)^(2k)dz
       =2^(4k)((2k)!)^2/(4k+1)!.

To show the every-order statement, let J_k(T)=integral_0^T(1-w^2)^(2k)dw. For every integer ell>=-k, coefficient extraction yields

    Lambda(t^ell p^k)=[t^(k-ell)] T^(k-1) J_k(T).

The derivative J_k'(T) vanishes to order 2k at T=1, so J_k(1+t)-J_k(1) is divisible by t^(2k+1). In the stated domain the required degree is at most 2k. Hence the exact identity is

    Lambda(t^ell p^k)=c_k binom(k-1,k-ell), ell>=-k.

In particular for EVERY m>=1, Lambda(p^m)=0 and Lambda(t p^m)=c_m, which is nonzero over EVERY characteristic-zero field. The restriction on ell is retained; for a fixed negative ell it is valid eventually.

Absolute unit-nonpower is especially simple here. Along the geometrically prime divisor t+1=0, the bracket 1-T^2S equals 1, t is a unit, and p has valuation exactly 1. Thus p cannot be a unit times an e-th power for any e>=2 after any algebraically closed extension. This argument is on the quotient itself, and needs no transport of an arbitrary upstairs divisor. At the actual S=0 point, p specializes to 1+t^-1, so its coefficients at exponents 0 and -1 are both nonzero. Its Laurent support is not contained in a translate of dZ for d>=2.

For completeness, every positive power p^a is pure. For that phase a fixed integer monomial t^ell is unsafe iff ell>=1, and it survives EVERY positive order iff 1<=ell<=a. These follow from the same binomial formula for k=am; for ell<1 take m large, and for ell>=1 take m large enough that k>=ell. If ell>a the first marked trace is zero, establishing necessity for the every-order statement. For p itself t is the unique every-order monomial marker. No finite-prefix inference enters these equivalences.

Over every formally real field the compact quotient factor is Mathieu: pullback preserves nonzero functions and square trace; use strict compact Haar positivity over R and the existing rational Gram/LDL argument. The scalar split torus is Mathieu over every characteristic-zero field. The explicit p,H pair therefore proves genuine two-factor Mathieu instability for all formally real fields.

**Theorem 2 (Albert full-normalizer and projective cubic products).** The compact rational EIV full-normalizer quotient and, separately, the projective cubic open U_J of every reduced Albert algebra over any characteristic-zero field have the explicit CDF phase E and unit marker t^-1 of Sections 5--6. Their marked values are (-1)^(m-1)/(m+1) at every positive order. Each phase is absolute unit-nonpower with primitive Laurent support. Only the compact rational quotient has the stated formally-real two-factor Mathieu instability conclusion.

## 5. EIV: full normalizer and unchanged invariant witness

Write J for the positive rational Albert algebra, e its unit, C=g e, and C# its quadratic adjugate. In the compact unitary model C#=tau C, where tau conjugates scalar i only. The existing regular coordinate and CDF are

    D=C11(C#)11=C11 tau(C11),   v=1-(1-D)^8.

We first prove the full normalizer assertion. The fixed vectors of the full stabilizer H=Aut(J) in the geometric Albert module form precisely the scalar line through e. This does not require an uninspected irreducibility claim: use the existing Peirce commutators Delta=2[L_(ei-ej),L_Hij]. For an invariant vector their diagonal coordinate equation Delta w=2d z forces every off-diagonal Peirce coordinate z=0; their Peirce coordinate equation Delta z=w_j-w_i then forces all diagonal coordinates equal. The Peirce coefficients d are nonzero. The argument holds after algebraic closure and descends.

If a normalizes H, it preserves this fixed line, so a e=zeta e. Cubic-norm preservation gives zeta^3=1. The scalar zeta lies in the geometric scalar center; after removing it a fixes e, hence belongs to H. Conversely H and that center normalize H, and their intersection is trivial because H fixes e. Thus N_G(H)=H Z(G), with geometric quotient N/H=mu3. For the compact Q-group Z(G) is the corresponding finite unitary form; we do NOT identify it with the split Q-group scheme mu3. This group-scheme assertion is checked geometrically and descends. Its action on X is free, since its intersection with every conjugate of H is trivial.

Under a central scalar, C maps to zeta C and C# maps to zeta^2 C#. Hence D is invariant. Equivalently in the compact real description zeta tau(zeta)=1, so C11 tau(C11) is unchanged. Thus the already canonical phase

    E=(1-t^-1)(1-v+vt),   H_marker=t^-1

is an actual function on G/N_G(H) times Gm. Its normalized quotient trace is the restriction of the old trace, so for every m>=1 it satisfies

    Lambda(E^m)=0,   Lambda(t^-1 E^m)=(-1)^(m-1)/(m+1).

All quantities are rational and give the same statement after EVERY characteristic-zero scalar extension. The prime divisor t=1 has valuation 1 in E, since the other factor restricts to 1, proving absolute unit-nonpower directly on the quotient. At the image of C=e, v=1, so the coefficients at exponents 1 and 0 are nonzero, proving primitive Laurent support. The compact quotient factor is Mathieu over all formally real fields by the same square-trace injection and rational Gram proof. No arbitrary witness or global Mathieu-descent theorem is invoked.

## 6. Reduced Albert projective cubic complements

There is a further separately scoped consequence over ANY characteristic-zero field K. Let J=Herm_3(O,gamma) be any REDUCED Albert algebra with its given frame, and let

    U_J=P(J) minus {N=0}.

The map Y_J={N(C)=1} to U_J sending C to its projective class is a mu3 torsor. Indeed over a projective chart where a homogeneous coordinate is nonzero, a representative has N invertible, and scaling to N=1 amounts to adjoining a root of a cubic equation with invertible constant and derivative; all choices differ by mu3. It is therefore finite etale of degree 3 and surjective, also when the required cubic root is not in K. In particular U_J is smooth geometrically integral affine: it is the standard projective open D_+(N), and smoothness follows by finite-etale descent from Y_J. It is Str_0(J)/(Aut(J) times mu3), its full normalizer quotient by the preceding fixed-line argument.

The actual degree-zero regular coordinate on this projective open is

    Dbar=C11(C#)11/N(C).

It pulls back to D on Y_J and is regular without making an arbitrary choice of projective representative. Put vbar=1-(1-Dbar)^8 and Ebar=(1-t^-1)(1-vbar+vbar t). The canonical reduced-Albert beta recurrence, restricted along injective pullback, proves the same every-order zero and marked trace above for Ebar on U_J times Gm. The same direct divisor t=1 and point [e] give absolute unit-nonpower and primitive Laurent support. There is NO compact-factor Mathieu conclusion for arbitrary reduced norm-one forms or these generally noncompact projective opens. Division Albert algebras without a frame are not included by this coordinate construction.

## 7. Literature and remaining inventory boundary

The normalizer circuit and high-order endpoint method are classical inputs already attributed in the canonical SCCP proof to Long's spherical draft, exact Sep 22 2026 revision e9edac262b4b64eb68feb0e6e80496cefe4a1dee, Section 5 and complete Proposition 5.1 proof. Its geometric normalizer example shows why quotient invariance must be checked; it does not supply an automatic arbitrary-phase descent theorem. The beta laws and full Albert model geometry are reused from the exact canonical documents named in Section 1, with their retained primary scopes.

As a narrow bibliographic comparison, R. Skip Garibaldi, The Rost invariant has trivial kernel for quasi-split groups of low rank, dated 30 October 2000, Example 3.4 printed page 12, explicitly identifies the projective stabilizer in the SPLIT Albert E6 representation as F4 times mu3. The actual whole example and its cited Section 2.4 setup were inspected. Its Jacobson transitivity and stabilizer references remain imports; no whole-paper or Rost-invariant proof audit is claimed. Retained [primary](sources/garibaldi-rost-kernel-2000.pdf), SHA256 `2eaacc14cb1f53bf7b765ec6a95173d8b3ca752e659c5d6cb6da0b7bf23cb75c`. It is comparison only, not a new necessary premise for the elementary fixed-line calculation here.

The two canonical parent scopes establish connected-stabilizer projector models and compact EIV, not explicitly these balanced full-normalizer quotients. The theorems above establish those named actual-coordinate descent cases only. A broader inventory, including other classical disconnected fixed groups, global central forms, exceptional cases and low-rank coincidences, requires separate exact matching. No absence of an equivalent literature result is inferred from targeted searches.
