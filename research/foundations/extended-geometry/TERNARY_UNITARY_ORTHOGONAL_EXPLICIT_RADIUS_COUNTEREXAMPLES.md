# Explicit radius counterexamples for ternary unitary and orthogonal quotients

The classical HCIZ integral gives a complete two-piece radius law. It yields a rational all-order phase, sharp pure-power and scalar-monomial marker profiles, and direct geometric certificates on the specified rational homogeneous models and their full normalizer quotients. This is an explicit branch, with no classification of all rational forms of the complex AI pair.

## 1. Specified rational data and theorem

Let d>0 be squarefree, E=Q(sqrt(-d)), and q any positive definite rational symmetric 3-by-3 matrix. Define the actual Q-groups

    G=SU(E/Q,q),       H=SO(q).

Here H is embedded by real matrices. Choose P in GL_3(Q) such that q0=P^t q P is positive diagonal. Put D=diag(-1,0,1), D_q=P D P^-1, and

    x=sqrt(-d)*D_q in Lie(G)(Q).

The matrix D_q is q-selfadjoint with trace zero, so x satisfies the unitary Lie equation and trace-zero condition. Its entries are rational coordinates relative to the fixed Q-basis 1,sqrt(-d) of E. This is an actual rational Lie point, not an assumed rational coweight.

For g in G set

    Z(g)=Ad(g^-1)x,       K(g)=(Z(g)+bar(Z(g)))/2,
    S(gH)=-(2/d)*tr(K(g)^2).                              (1)

The bar is the Q-linear involution of E. Section 2 proves that S is an actual regular Q-function on Y=G/H. It is also invariant under the full normalizer N_G(H), and thus defines the same function on Y_N=G/N_G(H).

**Theorem 1.** On either Y or Y_N, take the native normalized rational Reynolds trace lambda and compose it with Laurent constant term to obtain Lambda. Put

    T=1+t,    a(T)=(1-T^2*S)*(1-T^2*S/4),
    p=T*t^-1*a(T)^2,       F_a=p^a*(1+p),                 (2)

where a is a positive integer. Then the following assertions hold after every characteristic-zero scalar extension of the specified Q-model.

* The seed power p^a is pure if and only if a>=3.
* For every a>=3, both p^a and F_a are pure. For either phase, the integer monomial marker t^ell is unsafe if and only if ell>=2.
* For a>=3, t^ell marks EVERY positive order of p^a nontrivially if and only if 2<=ell<=a. It marks EVERY positive order of F_a nontrivially if and only if 2<=ell<=a+1.
* Each F_a with a>=3 is absolutely not a unit times a proper power and has primitive actual Laurent support.

In particular F=F_3 and the SAME Laurent UNIT H_marker=t^2 satisfy, for every m>=1,

    Lambda(F^m)=0,
    Lambda(t^2 F^m)=sum_(h=0)^m binom(m,h)
       *3*(81/16)^(3m+h)*((6m+2h)!)^2/(12m+4h+1)! >0.     (3)

The displayed marked values are positive rational numbers before scalar extension and remain nonzero over every characteristic-zero field. Over every formally real extension, each compact quotient factor and the split Laurent factor separately has a Mathieu kernel, while their product does not. The individual-factor scope and the all-field pair have different quantifiers.

Purity means Lambda(Q^m)=0 for every m>=1. Unsafe means nonzero marked moments for infinitely many orders. The theorem classifies scalar Laurent monomials for these displayed phases, not arbitrary function-valued markers, optimal degrees or every compact rational AI form.

## 2. Actual quotient, projection and rational trace

The actual compact real group is a form of SU(3), and H is a compact SO(3). Connected reductivity, the affine direction of Matsushima and the represented homogeneous quotient give Y smooth, geometrically integral and affine. Its identity coset is a rational point. The native normalized rational trace is the invariant projection to constants, commutes with characteristic-zero scalar extension, and over R equals compact Haar integration after pullback to G. Only these generic quotient/trace facts are used, as recorded in [RHT], [HQP] and [SCCP]. No earlier counterexample theorem is a premise.

Write a unitary Lie element as A+sqrt(-d)*B with rational real matrix coordinates. Its equations give A in so(q) and B q-selfadjoint with trace zero. These summands are orthogonal for the negative ambient Killing form, by tr(AB)=0. Thus the real part K in (1) is precisely the rational Killing-orthogonal projection to h=so(q). Projection is equivariant under FULL H. Therefore tr(K^2), and hence S, is right-H-invariant and descends regularly by the quotient torsor. Individual entries of K are not claimed to descend.

The rational basis change g->P^-1 g P reduces q to q0. Let C be the positive real diagonal square root of q0. Conjugation by C identifies the real group pair with the standard (SU3,SO3). It commutes with D, sends x to sqrt(d) times iD, and sends K to sqrt(d) times the projected standard matrix. The factor 1/d in (1) therefore makes S exactly the standard normalized radius square. This is an actual real isomorphism applied AFTER defining the rational function on the specified Q-model. It does not infer a Q-form from a real label.

For the standard normalization let M=Ad(SU3)(iD), of real dimension 6. If K is a real antisymmetric matrix representing a physical vector of length r0, then tr(K^2)=-2r0^2. Thus S=v^2 with v=2r0. Native inverse-Haar measure has the same law as Haar, so the projection law on M gives the polynomial functional of S on Y. No pointwise identification Y=M is asserted: Y has dimension 5, while M has dimension 6.

For comparison, zero is a GLOBAL regular value of the h-moment map on M. A zero-fiber point is iC0 with C0 real symmetric and distinct eigenvalues -1,0,1. Orthogonal diagonalization puts it in the SO3 orbit of iD; a determinant-minus-one diagonalizer can be adjusted by a diagonal reflection commuting with D. Its so3-centralizer is zero, since a real antisymmetric matrix commuting with a distinct diagonal matrix is zero. This also verifies the one-point input of the general nonabelian theorem. The explicit construction below does not import that theorem's density or Chevalley-lifting argument.

## 3. Classical scalar law and the exact radius functional

Take a real rotation generator J in so3 with eigenvalues i,-i,0 and put s=-tr(zJ) for z in M. The Hermitian matrices -iz and -iJ have eigenvalues -1,0,1, so s has endpoints -2 and 2. Normalized SU3 and U3 conjugation integrals agree, because scalar matrices act trivially by conjugation.

Use the classical HCIZ formula [McS, equation (2)] with A=D and B=aD, initially a>0. Their Vandermondes are 2 and 2a^3, and the normalization product is 1!*2!=2. The exponential determinant is

    exp(2a)-2exp(a)+2exp(-a)-exp(-2a).

Therefore the entire moment-generating function, with its removable value 1 at a=0, is

    E exp(a*s)=[exp(2a)-2exp(a)+2exp(-a)-exp(-2a)]/(2a^3).    (4)

Analytic continuation gives every complex a. Direct integration shows that (4) is the Laplace transform of the even probability density

    rho_A(s)=(2-s^2)/4       for |s|<=1,
             (2-|s|)^2/4   for 1<=|s|<=2,
             0             otherwise.                       (5)

Its integrals sum to 1. Equality of compactly supported measures follows from Fourier uniqueness, or from all polynomial moments and polynomial density in continuous functions.

The projected so3 law is rotation invariant. Conditional on a nonzero radius v, the scalar component s is uniform on [-v,v], by the ordinary round-sphere latitude formula. If nu is the law of v>=0, then for s>0

    rho_A(s)=integral_[s,infinity) (1/(2v)) dnu(v).           (6)

An atom at v=0 would produce a scalar atom at 0, absent in (5). An atom at a positive radius would give a jump in rho_A, also absent. Differentiating (6) as measures gives dnu(v)=-2v rho_A'(v)dv. Hence the EXACT normalized radius density is

    rho(v)=v^2            for 0<v<1,
           v*(2-v)       for 1<v<2,
           0             otherwise.                         (7)

This uses the actual invariant S and sphere mixtures; it does not average a phase while assuming its purity unchanged. The law on Y for every specified d,q is the same by Section 2. Elementary integration gives, for EVERY j>=0,

    lambda(S^j)=integral_0^1 v^(2j+2)dv
                 +integral_1^2 v^(2j+1)*(2-v)dv
               =(4^(j+1)-1)/[(j+1)*(2j+3)].                 (8)

These values are also native rational trace values, so the complete polynomial functional commutes with every characteristic-zero extension. As a normalization check, the mean is 3/2; no bounded list is used to infer the law.

## 4. Complete cancellation and exact exceptional moments

Put f(v)=(1-v^2)*(1-v^2/4). Its nonzero boundary roots are 1 and 2, and f(0)=1. For j=1,2 and each interval [alpha,beta] equal to [0,1] or [1,2], define

    J_jk(T)=T^(j+1)*integral_alpha^beta v^j*f(Tv)^(2k)dv.

This is a polynomial in T. Direct differentiation gives

    J_jk'(T)=T^j*[v^(j+1)*f(Tv)^(2k)]_alpha^beta.

At a nonzero endpoint f(v)=0, so f(Tv)^(2k) is divisible by (T-1)^(2k). At the radial endpoint 0 the factor v^(j+1) kills the flux. Thus J_jk(1+t)-J_jk(1) is divisible by t^(2k+1). Applying (7) coefficientwise gives the full exact mixed identity

    Lambda(t^ell p^k)
      = binom(k-3,k-ell)*integral_0^1 v^2*f(v)^(2k)dv
        +integral_1^2 [2v*binom(k-2,k-ell)
                       -v^2*binom(k-3,k-ell)]f(v)^(2k)dv,   (9)

for EVERY k>=1 and integer ell>=-k. Here binom(n,j)=[t^j](1+t)^n for integer n, and is zero when j<0. Indeed p^k contributes t^-k times T^k and the discarded Taylor error begins at degree 2k+1; the requested degree k-ell is at most 2k. Negative powers of T are formal units at t=0. This is a complete polynomial argument, not an asymptotic or finite moment-prefix inference.

For k>=3 the upper binomial indices in (9) are nonnegative and strictly below k. Consequently

    Lambda(p^k)=0 for EVERY k>=3.                            (10)

The only lower exceptional pure values are

    Lambda(p)=77/480,       Lambda(p^2)=-1/15.               (11)

An exact finite-expansion formula specifying this arithmetic is: if c_kj is the coefficient of u^j in (1-5u/4+u^2/4)^(2k), then

    Lambda(p^k)=sum_(j=0)^(4k) c_kj
                  *(4^(j+1)-1)/[(j+1)*(2j+3)]*binom(k+2j,k).

Taking k=1,2 gives (11). Therefore, for a positive integer a,

    p^a is pure IFF a>=3.                                  (12)

Sufficiency uses the whole tail (10); necessity is the first pure moment in (11), not a sampled inference. For a>=3, expansion of F_a^m into seed powers p^(am+h), 0<=h<=m, proves its purity as well.

## 5. Sharp unsafe and every-order scalar-monomial profiles

At ell=2 and k>=3, (9) becomes a strictly positive beta integral:

    B_k:=Lambda(t^2*p^k)
       =2*integral_1^2 v*f(v)^(2k)dv
       =integral_1^4 [(1-z)*(1-z/4)]^(2k)dz
       =3*(81/16)^k*((2k)!)^2/(4k+1)! >0.                  (13)

Use z=v^2 followed by z=1+3u to obtain the last equality. Expanding F_a^m gives sum_(h=0)^m binom(m,h)B_(am+h), proving (3) at a=3 and a positive rational marked value at EVERY order for every a>=3.

For every fixed integer ell<=1, (9) vanishes for all sufficiently large k: its lower binomial index exceeds both nonnegative upper indices once the stated domain applies. Thus such markers are eventually safe against every p^a and F_a with a>=3.

For ell>=3 and k>=max(3,ell), put B0=binom(k-3,ell-3)>0. The Origin contribution in (9) is B0 times a strictly positive integral. Its Away integrand is

    v*B0*[2*(k-2)/(ell-2)-v]*f(v)^(2k).

It is nonnegative on 1<v<2 because k>=ell. Hence the complete value is positive. Along with (13), EVERY ell>=2 is eventually positive. For the displayed pure seed and modified families,

    t^ell is unsafe IFF ell>=2.                             (14)

The every-order statements are sharper. Fix a>=3. If 2<=ell<=a, then every k=am is at least max(3,ell), so every t^ell-marked moment of p^a is positive. If ell>a, its first marked moment is zero because p^a has no Laurent exponent below -a. If ell<=1, eventual vanishing already excludes every-order survival. Thus

    t^ell marks EVERY order of p^a IFF 2<=ell<=a.            (15)

For F_a the same argument gives positivity at every order when 2<=ell<=a. At ell=a+1, the first order has zero contribution from p^a and a positive contribution from p^(a+1). For m>=2 every exponent am+h is at least a+1, so all its contributions are positive. If ell>a+1, the first marked moment is zero because F_a has no exponent below -(a+1). Therefore

    t^ell marks EVERY order of F_a IFF 2<=ell<=a+1.          (16)

In particular the least unsafe exponent is 2. The first local homogeneous density degree on the Origin interval is j0=2, whereas the Away density has a nonzero lower-degree term 2v. This explicitly shows why a local leading degree alone does not imply least unsafe=r+j0: the general leading-degree marker would be 3 here, while the actual least exponent is 2. No classification of arbitrary function-valued markers is asserted.

## 6. Full normalizer and actual adjoint quotient

We prove the FULL equality, scheme-wise,

    N_G(H)=H*Z(G).                                          (17)

Over Qbar, identify G with SL3 and H with SO3 in its defining representation, after changing q to the identity. Write L_ij=E_ij-E_ji. A matrix commuting with all L_ij also commutes with their squares. The squares for (1,2) and (1,3) have three distinct joint coordinate eigenspaces, so that matrix is diagonal. Commutation with the L_ij then makes all diagonal entries equal. Thus the full matrix commutant of h is scalar, and the center of H is trivial: a scalar in H has square 1 and determinant 1 in odd dimension, hence is 1.

Here is a local three-dimensional proof that every automorphism of h=so3 is inner, sufficient for the normalizer argument. An automorphism preserves its nondegenerate Killing form B_h and the nonzero alternating three-form omega(u,v,w)=B_h([u,v],w). It therefore belongs to SO(h,B_h). The adjoint map H->SO(h,B_h) is faithful by the scalar-commutant calculation. Its closed image has dimension 3, equal to the dimension of the connected target SO(h,B_h), so the map is onto. Thus every Lie automorphism is Ad(h0) for some h0 in H.

For n normalizing H, conjugation on h is such a Lie automorphism. Choose h0 with the same action. Then h0^-1*n commutes with h, hence is scalar. Determinant 1 puts it in the center of G. This proves (17) over Qbar; conversely H and the center plainly normalize H. The two Q-defined smooth subgroups are equal after faithful scalar extension, so (17) descends to Q. No rational representative of a central component and no automorphism classification by name is required.

The intersection H with Z(G) is trivial, by the same scalar square/determinant argument. Hence N_G(H)/H is finite etale (geometrically of order 3). Its right action on Y is free. The actual quotient Y_N=G/N_G(H) is affine, smooth and geometrically integral, and Y->Y_N is finite etale. It is equivalently the corresponding homogeneous quotient of the adjoint unitary group G/Z(G).

For a central z, Ad((gz)^-1)x=Ad(g^-1)x. Thus S and every phase coefficient in (2) are ACTUALLY invariant under the full normalizer and descend. The native normalized traces satisfy exact restriction under Y->Y_N, so every identity (8)--(16) remains true downstairs. This is descent of a specifically invariant phase, not a theorem that arbitrary upstairs witnesses descend.

## 7. Direct geometric certificates and separate factors

On both Y and Y_N, the actual rational identity image has S=0. There p=1+t^-1. At t=-1/2 the modifier 1+p has derivative -4 and p=-1 is a unit. Smoothness supplies a geometric prime divisor with valuation one in F_a for every a>=1. A unit times an e-th power has valuations divisible by e at every geometric prime, so F_a is absolutely not a unit times a proper power. The same actual point and derivative persist after every algebraically closed extension. No arbitrary nonpower descent is used.

At that actual point

    F_a=(1+t^-1)^a*(2+t^-1).

Its coefficients at adjacent exponents 0,-1 are 2,2a+1, nonzero in every characteristic-zero field. This proves actual Laurent-support primitivity, in particular coefficients 2 and 7 for the main F_3.

Every nonzero real regular function on either quotient pulls back injectively to the connected compact group, giving strictly positive square trace on its dense real orbit. Finite rational Gram/LDL and sums of rational squares make this square trace anisotropic after every formally real extension; only zero is pure on that factor. The scalar Laurent factor is separately Mathieu over every characteristic-zero field by scalar Duistermaat--van der Kallen and the finitely generated coefficient-field argument. These are the generic facts in [SCCP]. Thus the all-order pair proves genuine two-factor instability over formally real fields, while the pair itself remains valid over all characteristic-zero extensions.

## 8. Exact primary scope, dependencies and attribution

[McS] Colin McSwiggen, A new proof of Harish-Chandra's integral formula, Communications in Mathematical Physics 365 (2019),239--253, DOI 10.1007/s00220-018-3259-9. Exact arXiv:1712.03995v4, 27 April 2020; [primary](https://arxiv.org/abs/1712.03995v4), [retained PDF](sources/mcswiggen-harish-chandra-1712.03995v4.pdf), SHA256 dc50b3c21828d4f3a0099d95b8961ba5002ceca996776eea9c95449d14519a5d. Imported scope: equation (2) and U(N)/SU(N) discussion printed 2; Section 2 notation and Theorem 1 printed 3--4; ENTIRE Section 3 heat-flow proof printed 5--10, including Lemmas 2,5--7 and normalization 3.4. The source expressly imports radial differential-operator existence (Theorem 3) and the radial Laplacian identity (Proposition 4) from Helgason/Harish-Chandra. Their earlier proofs, ordinary heat-kernel uniqueness and steepest descent remain classical analytic inputs. No whole-paper or historical-proof certification is claimed. Its later symplectic localization Sections 4.2--4.3 and original DH Section 4 are not premises here.

[RHT] [Normalized reductive homogeneous trace](REDUCTIVE_HOMOGENEOUS_AUTOMORPHISM_INVARIANT_TRACE.md) and [HQP] [homogeneous quotient traces](HOMOGENEOUS_QUOTIENT_TRACE_AND_COUNTEREXAMPLE_PULLBACK.md) supply only represented affine quotient/native normalized trace/basechange scopes. [UCQ] [compact rational coadjoint products](UNIVERSAL_COMPACT_RATIONAL_COADJOINT_PRODUCT_COUNTEREXAMPLES.md), Sections 2--3, supplies generic compact orbit/Haar facts. [SCCP] [standard compact classical products](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md), Sections 3 and 13, supplies rational Gram/SOS and scalar Laurent-factor facts. No prior model-specific moment law or counterexample is used.

The cellwise endpoint-normalizer method is classical, with Long attribution recorded in [the radial invariant-test criterion](RADIAL_WALL_INVARIANT_TEST_NORMALIZER_CRITERION.md) and its predecessors. Section 4 writes the complete two-interval Euler calculation locally. The HCIZ scalar law is classical; the actual radius mixture, native rational normalization, complete marked beta identity and sharp profiles are derived here. This explicit branch needs no CDKW derivative principle or Chevalley restriction theorem.

Known complex all-order pairs already imply complex Laurent-unit/nonpower/support conclusions by P=f(1+th), as recorded in [the nonabelian zero-submersion theorem](COMPACT_RATIONAL_NONABELIAN_ZERO_SUBMERSION_COUNTEREXAMPLES.md), Section 10, with Liu's antecedent. The present scope is the explicitly specified E,q rational models and full normalizer quotients, their native rational values and characteristic-zero scalar extensions. The parameterization is not asserted to exhaust every abstract compact rational AI form. No new complex-existence, classification or worldwide-priority claim is made.
