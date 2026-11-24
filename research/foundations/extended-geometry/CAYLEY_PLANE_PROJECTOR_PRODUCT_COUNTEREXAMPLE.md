# Cayley-plane projector products: a direct Albert-algebra proof

The diagonal coordinate of a primitive Albert idempotent has a beta law obtained directly from sixteen Jordan derivations. Its polynomial cumulative distribution gives an explicit pure phase and a fixed unit marker on the product with a split torus.

## 1. Exact statement and actual coordinates

**Theorem 1.**

Let O_Q be the standard Cayley-Dickson octonion algebra over Q with its usual basis1,e1,...,e7, conjugation, and norm sum of eight squares. Put

    J_Q=Herm_3(O_Q), X_Q={P in J_Q : P o P=P, tr(P)=1},
    P o Q=(PQ+QP)/2.

X_Q means the affine primitive-idempotent variety. It is smooth and geometrically integral, geometrically homogeneous F4/Spin9, of dimension16. Its real form X_R is the compact Cayley plane. The projective complex E6/P1 Cayley plane is a different variety; we use its affine trace-normalized idempotent chart. Section2 specifies the classical primary geometry input and the equation-scheme tangent check.

Let K be any characteristic-zero field and X_K=X_Q tensor K. Its normalized Aut(J_K)-invariant trace lambda_K exists and is unique. Define the actual scalar regular coordinate

    S(P)=P11.

For EVERY integer j>=0,

    lambda_K(S^j)=(4)_j/(12)_j.                       (1.1)

Equivalently, for every polynomial g over K,

    lambda_K(g(S))=1320 integral_0^1 u^3(1-u)^7 g(u)du. (1.2)

These integrals are finite polynomial integration over Q, extended K-linearly. Over R, S really has the Beta(4,8) distribution under compact Haar measure.

Define the integer polynomial

    v(S)=330S^4-1848S^5+4620S^6-6600S^7
          +5775S^8-3080S^9+924S^10-120S^11.          (1.3)

On the ACTUAL product X_K x G_m, with t the split radial coordinate, put

    F=(1-t^-1)(1-v+v t), H=t^-1.                    (1.4)

With Lambda=lambda_K tensor CT_t, for EVERY m>=1,

    Lambda(F^m)=0,
    Lambda(H F^m)=(-1)^(m-1)/(m+1) !=0.             (1.5)

The SAME fixed regular unit marker H works at every order. F is not a unit times any proper power even after algebraic closure. Thus the product trace kernel is not Mathieu for every such base-changed Q-model over any characteristic-zero K.

For K=R, and more generally every formally real K, each factor X_K and G_m separately DOES have the Mathieu property. Thus these give genuine product instability, including the exceptional compact rank-one symmetric space over R. This does not assert the same theorem for an arbitrary Albert form, an arbitrary exceptional homogeneous variety, or a general real descent from complex failure.

## 2. Classical affine geometry and the trace

The needed classical facts are the following.

* The algebra J_R is the Euclidean Albert Jordan algebra, with positive trace form <A,B>=tr(A o B).
* Its connected algebraic automorphism group is the compact real form of F4. It acts transitively on primitive idempotents; the stabilizer of diag(1,0,0) is Spin9.
* Over an algebraic closure of a characteristic-zero field, the trace-one idempotents are the smooth affine homogeneous variety F4/Spin9.

The classical inputs are Brada, Elements de la geometrie des octaves de Cayley (1986), Sections9-11 and the trace/real-diagonalization portions of12/14 specified in Section8; and Petrov-Semenov, Rost motives, affine varieties, and classifying spaces, exact arXiv1506.07788v1, Section4 printed7 BEFORE Theorem4.2. The affine F4/Spin9 identification is taken from the latter source independently of Brada's stabilizer topology argument. Petrov-Semenov explicitly identifies the16-dimensional affine variety of primitive idempotents as a twisted form of F4/Spin9. Their following motive theorem assumes a first Tits construction; that theorem and that restriction are NOT imported here.

For clarity about the equation scheme, at E=diag(1,0,0), differentiating P o P=P gives

    2E o delta P=delta P.

The Peirce1 diagonal line and Peirce0 lower Hermitian2x2 block have zero tangent component. The two octonionic off-diagonal blocks involving row1 are free, giving16 dimensions; their trace is zero automatically. This equals the orbit dimension52-36. Geometric transitivity gives the same calculation at every geometric point. Thus the equations define the reduced smooth orbit, without hidden nilpotent components.

All multiplication constants are rational. The standard Q-algebra J_Q is an Albert algebra and these classical geometric statements apply after algebraic closure. In particular X_Q is smooth geometrically integral. Connected reductivity and local finiteness/complete reducibility give the unique normalized Reynolds functional on its coordinate ring, and this functional commutes with every characteristic-zero scalar extension. The canonical normalized trace over R agrees with compact Haar averaging by invariance and uniqueness. No assertion that every real point of an arbitrary form lies in one compact orbit is used.

## 3. Why the required Jordan operators are genuine Lie derivations

For a Jordan algebra let L_A(B)=A o B. Every commutator [L_A,L_B] is a derivation. The short algebraic proof is included to avoid an associative-matrix assumption.

The full polarization of the Jordan identity gives

    ((a o b)o d)o c + ((b o c)o d)o a + ((c o a)o d)o b
      =(a o b)o(d o c)+(b o c)o(d o a)+(c o a)o(d o b).

Subtract this identity with triple(b,c,z),middle a from the identity with triple(a,c,z),middle b, or take its negative. The paired product terms cancel by commutativity; the remaining equality is exactly

    [L_a,L_b](c o z)=([L_a,L_b]c)o z+c o([L_a,L_b]z).

Thus these operators lie in Lie Aut(J). This holds over every characteristic-zero field. The trace form is associative, so L_A,L_B are self-adjoint and their commutator is skew-adjoint. Since the commutator kills the unit, it also kills trace. Its linear action therefore preserves the idempotent variety and induces a derivation on its coordinate ring. The invariant Reynolds functional annihilates every such derivation.

This general Jordan argument uses no false associative identity for octonion matrices. Brada also gives a direct proof for the relevant zero-trace skew-Hermitian commutator operators at printed30-33; it agrees with the specific operators used next.

## 4. Sixteen coordinate derivations and the exact recurrence

Write E_i for the three real diagonal matrix units. For i=2,3 and an orthonormal octonion basis unit e, let H_i(e) be Hermitian with e at(1,i), conjugate(e) at(i,1), and other entries0. Set

    A_i=E_1-E_i,
    D_(i,e)=2[L_(A_i),L_(H_i(e))],
    x_(i,e)=Re(P1i conjugate(e)).

Direct calculation of ONLY the indicated matrix entries gives

    D_(i,e) S=2 x_(i,e),
    D_(i,e) x_(i,e)=Pii-S.                          (4.1)

Here D acts on a coordinate function f by the derivative of f(exp(tD)P) at0. Choosing the opposite conventional sign for the group action changes both terms and leaves invariance unchanged.

For completeness, L_A scales a(1,i) entry by(1-1)/2=0. The(1,1) entry of L_H P is x, while the same entry of L_H L_A P is0. Therefore2[L_A,L_H]S=2x. The(1,i) entry of L_H P is (Pii+S)e/2; after applying L_A it is0. The corresponding entry of L_H L_A P is(S-Pii)e/2. This proves the second formula in(4.1), because e has norm1. No product of three arbitrary octonions occurs.

The(1,1) idempotent equation gives

    sum_(i=2,3;e) x_(i,e)^2=S-S^2.                  (4.2)

The trace-one equation gives

    sum_(i=2,3;e)(Pii-S)=8(1-3S).                   (4.3)

Let mu_j=lambda_K(S^j), mu_0=1. For j>=1, apply invariant lambda_K to the sum of D_(i,e)(S^(j-1)x_(i,e)). Equations(4.1)-(4.3) yield

    0=2(j-1)(mu_(j-1)-mu_j)+8(mu_(j-1)-3mu_j)
     =(2j+6)mu_(j-1)-(2j+22)mu_j.

At j=1 the first differentiated term is zero before substitution, so there is no inversion of S. Consequently

    (j+11)mu_j=(j+3)mu_(j-1),
    mu_j=(4)_j/(12)_j.

This proves(1.1) in all orders over K itself. The number8 is the dimension of the octonion Peirce space; the two choices i=2,3 supply exactly16 derivations.

## 5. Beta law and the actual polynomial CDF

Elementary beta integration gives

    1320 integral_0^1 u^(j+3)(1-u)^7du
      =(4)_j/(12)_j,

since B(4,8)=3!7!/11!=1/1320. Hence(1.2) follows by polynomial linearity. Over R, equation(4.2) shows0<=S<=1. Its compact-Haar pushforward and the beta probability measure have the same integrals of every polynomial; polynomial density in continuous functions on[0,1] then gives equality of the probability measures. This identifies the claimed regular diagonal coordinate, without a polar-distance or octonionic homogeneous-vector parametrization.

The explicit polynomial(1.3) satisfies

    v'(u)=1320 u^3(1-u)^7, v(0)=0, v(1)=1.

For every j>=0, polynomial antiderivatives therefore give

    lambda_K(v(S)^j)
      =integral_0^1 v(u)^j v'(u)du=1/(j+1).          (5.1)

All exchanges are finite polynomial identities; these statements do not need any analytic integration over K.

## 6. The interval-Laurent pair, all orders, and absolute nonpower

Apply lambda_K coefficientwise to F^m from(1.4), using(5.1). The formal uniform variable v gives

    lambda_K(F^m)
      =((t-1)/t)^m integral_0^1[1+(t-1)v]^m dv
      =(t-1)^(m-1)(t^(m+1)-1)/[(m+1)t^m].          (6.1)

The first summand has Laurent exponents1,...,m; the second has exponents-m,...,-1. Its constant term is zero. After multiplication by H=t^-1, the first contributes(-1)^(m-1)/(m+1), and the second still has strictly negative exponents. This proves BOTH identities(1.5) at every positive order, with one H fixed.

For absolute nonpower, complexify or pass to any algebraic closure of K. The divisor t-1 is geometrically prime because X is geometrically integral, and its local ring at the generic point is a discrete valuation ring because X x G_m is smooth. In

    F=(t-1)/t * (1-v+v t),

the second factor equals1 modulo t-1. Thus ord_(t-1)(F)=1. Units have valuation0 and e-th powers have valuation divisible by e. No unit times a proper power can equal F. The marker H itself is a regular unit.

No finite cover, marker averaging, phase descent or complex-failure-to-real-failure inference is used. The phase is written over Q directly in the single regular diagonal coordinate S and the actual split radial coordinate t.

## 7. Why each factor is Mathieu over R and every formally real field

Over R a nonzero real regular f has lambda_R(f^2)>0. Indeed its pullback to the algebraic group is nonzero; real compact group points are Zariski dense, so it is nonzero somewhere on the compact orbit, and continuity supplies a positive-Haar-measure neighborhood where f^2>0. Thus the only pure real phase is0, which is safe.

The formally real extension also has a direct finite Gram proof. Fix a nonzero f in O(X_Q) tensor K, and write it as a finite sum sum_i c_i b_i with c_i in K and b_i Q-linearly independent. Faithful scalar extension makes the b_i R-linearly independent. The rational matrix

    G_ij=lambda_Q(b_i b_j)

is positive definite over R by the preceding Haar-square argument. Rational LDL^t elimination diagonalizes it with strictly positive rational pivots d_i. Each positive rational d=a/b is a sum of ab copies of(1/b)^2; thus d is a sum of rational squares without adjoining roots. Over a formally real field, a sum of squares can vanish only when each summand is zero. It follows that the same rational positive-definite Gram form is anisotropic over K. Hence lambda_K(f^2)!=0, and again only0 is pure. This is a fixed finite-span argument, not generic global Mathieu ascent.

The split radial torus has the Mathieu property over every characteristic-zero field by the Laurent constant-term theorem: restrict a given finite coefficient field to C, apply the Duistermaat--van der Kallen Newton theorem, then use strict support separation to overwhelm each fixed marker. The existing sharp-quadratic proof already gives this all-field argument. Thus for every formally real K the two separate factors are Mathieu while their product fails by(1.5).

For non-formally-real K the product bad pair remains valid by the direct algebraic trace recurrence, but this proof makes no claim that the Cayley factor is Mathieu.

## 8. Primary inputs and the attribution boundary

**Brada.** Claude Brada, Elements de la geometrie des octaves de Cayley, Publications du Departement de Mathematiques de Lyon,1986,fascicule2D,1-90. The trace form and Albert product are in Sections9-10, printed28-32; the derivation characterization is in Section11-1/11-2, printed32-34; preservation of the trace/Euclidean form is in Section12-1, printed46; the real diagonalization theorem and its proof are in Section12-2, printed47-50.

[Original primary](https://www.numdam.org/article/PDML_1986___2D_A1_0.pdf), [retained copy](sources/brada-cayley-geometry-1986.pdf).

A specific source error is excluded: printed58 asserts that the stabilizer is globally homeomorphic to S^8 x Spin8. That global product/trivialization assertion is erroneous and is not a premise here. The sphere bundle Spin8 -> Spin9 -> S^8 is locally trivial, not the asserted global product. The trace, derivation and real diagonalization facts used above are independent of that assertion. We use the affine homogeneous/stabilizer identification from Petrov--Semenov separately.

**Petrov--Semenov.** Victor Petrov and Nikita Semenov, Rost motives, affine varieties, and classifying spaces, exact arXiv1506.07788v1, Section4 printed7, geometric prelude to Theorem4.2. For an Albert algebra in characteristic different from2,3, the paragraph identifies the affine16-dimensional primitive-idempotent variety as a twisted F4/Spin9, distinct from the projective E6/P1 variety. The following first-Tits/Rost-motive theorem is not used.

[Exact primary version](https://arxiv.org/pdf/1506.07788v1), [retained copy](sources/petrov-semenov-1506.07788v1.pdf).

**Normalized trace and reusable functional lemmas.** Characteristic-zero rational complete reducibility and the same-field Reynolds projection are as in Milne, Algebraic Groups, Theorem22.42/Corollary22.43 and Section12l. Uniqueness and base change are explained in Section2 above. The general compact-Haar positivity, rational Gram, polynomial-beta CDF and interval--torus lemmas also appear in [Standard compact classical product counterexamples](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md), Sections3,5 and13. Their hypotheses are verified here by the actual Albert geometry and recurrence; this does not import a theorem restricted to the seven classical families.

The scalar torus-factor assertion uses J.J. Duistermaat and W. van der Kallen, Constant terms in powers of a Laurent polynomial, Indag.Math.9(1998),221-231, [DOI](https://doi.org/10.1016/S0019-3577(98)80020-7). Only scalar Laurent coefficients are involved. No function-valued version of that theorem or generic global Mathieu ascent is used.

**Known interval construction.** The interval--torus seed in Section6 is Christopher D. Long's known example, Counterexamples to the xz-Conjecture and the Mathieu Conjecture for SU(2), arXiv2607.19012v1, Theorem2.1 and Proposition3.1. The uniform coordinate is supplied here by the displayed regular polynomial CDF, and the complete finite polynomial calculation is reproduced above. No priority is claimed for the beta law, the seed or the classical structural inputs.

The proof concerns the standard Q-defined Albert model and its scalar extensions. It does not classify arbitrary Albert forms, arbitrary exceptional homogeneous varieties or arbitrary finite quotients. The direct Jordan--Lie recurrence supplies the coordinate law; no radial-metric identification or complex-failure-to-real-descent inference is required.
