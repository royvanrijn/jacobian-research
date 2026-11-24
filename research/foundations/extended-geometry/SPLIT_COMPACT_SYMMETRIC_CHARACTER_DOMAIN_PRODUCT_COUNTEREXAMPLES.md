# Character-domain product counterexamples for compact split-type symmetric spaces

The construction below applies to specified compact rational Chevalley models and their central quotients. It uses the [uniform bounded-domain normalizer criterion](UNIFORM_BOUNDED_DOMAIN_NORMALIZER_CRITERION.md). Polar integration, the character Jacobian and the one-dimensional endpoint mechanism have classical antecedents, recorded with their exact scopes in Section 9.

## 1. Models and statement

Fix a split simply connected simple Chevalley group G_sp overQ, of rank r>=1, with its standard torus, simple roots and root maps. Let omega be the Chevalley involution

    omega(H_i)=-H_i, omega(E_alpha)=-E_-alpha.

Put L=Q(i). Descend G_sp,L using the semilinear involution omega composed with the nontrivial automorphism of L/Q. Denote this compact Q-form by G. Its real compact Lie algebra has rational basis

    iH_i, E_alpha-E_-alpha, i(E_alpha+E_-alpha).

Because omega commutes with this descent datum, it defines a Q-involution of G, again denoted omega. Put H=G^omega, the FULL fixed subgroup, and X=G/H. Define the actual Cartan map

    c(gH)=g omega(g)^-1.

Section 2 proves that X is smooth geometrically integral affine and that c is Q-regular. Let Phi=(Phi_1,...,Phi_r) be the real rational fundamental-character coordinates constructed there. There is an explicit polynomial d in Q[v_1,...,v_r] whose restriction to a maximal torus is

    d(Phi(a))=product_(alpha root)(1-alpha(a)).

The compact character domain Omega_bar=Phi(G(R)^0) has a bounded nonempty interior Omega. Its boundary lies in d=0, and the normalized trace of P(Phi(c)) is uniform Lebesgue integration over Omega.

For ANY q in Q^r intersect Omega with d(q)!=0, put

    y=Phi(c)-q,
    f_q(y)=d(q+y)/d(q), T=1+t,
    p=T/t f_q(Ty)^2,
    F=p^r(1+p), H_marker=t^r.                              (1.1)

Such rational q exist. They are actual real image values of Phi(c), not merely points in an abstract ambient coordinate space.

**Theorem 1.** For the specified X and EVERY characteristic-zero field K, the same Q-regular functions in (1.1) satisfy, for EVERY integer m>=1,

    Lambda_K(F^m)=0,
    Lambda_K(H_marker F^m)
      =sum_(h=0)^m binom(m,h) mu_(rm+h) !=0,
    mu_k=lambda_Q(f_q(y)^(2k)) in Q_(>0).                 (1.2)

Here Lambda=lambda_X tensor Laurent constant term. The marker is a unit. F is not a unit times any proper power over Kbar, and its nonzero Laurent coefficient functions have exponent differences with gcd1. Over every formally real K, each individual factor X_K and splitGm has Mathieu trace kernel, while their product does not.

The seed has sharp converses: for every positive integer a, p^a is pure IFF a>=r. For a>=r, a fixed monomial t^ell is unsafe for p^a IFF ell>=r, and it survives every positive order IFF r<=ell<=a. For F, the unsafe monomial markers are exactly ell>=r, and every-order survival occurs exactly for r<=ell<=r+1.

For r=6,7,8 the corresponding compact split-type exceptional quotients are EI, EV and EVIII, respectively. Section 8 specifies their actual finite-center isotropy groups. Section 7 gives a distinguished rational q and proves direct invariance of the CONSTRUCTED functions under every central isogeny quotient, including any additional fixed-subgroup components downstairs.

## 2. Rational quotient, trace and character coordinates

The descent above is explicit finite Galois descent of the split affine group. OverL it is split, and its compact real Lie algebra is exactly the displayed Chevalley Q-span. Let U=G(R)^0, a connected compact real group. The fixed subgroup H is Q-defined and smooth in characteristic0. Its identity component is reductive: the unipotent radical is defined overR by uniqueness and descent, but a nontrivial real unipotent group has nonidentity real unipotent points, while a compact matrix group has none. Thus H is reductive, allowing finite components. Matsushima over Qbar and affine descent give X affine. Homogeneous-space smoothness and surjectivity from the geometrically integral G give smooth geometric integrality of X.

The formula for c is unchanged by g->gh for h fixed by omega, so it descends to the FULL quotient, not just G/H^0. Normalized Reynolds projection exists, is unique, and commutes with all field extensions. After R-extension its restriction to the dense compact orbit is Haar probability, by the generic trace lemma in [Standard compact classical products](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md), Section 3. No assertion that all real points lie in this orbit is needed.

OverL, let chi_i be the irreducible fundamental characters. Steinberg's character theorem says that they freely generate the algebra of regular conjugation-invariant functions of the simply connected group and separate its semisimple conjugacy classes. Under the compact descent, Galois conjugation sends chi_i to its dual chi_(i*). For each exchanged pair use

    (chi_i+chi_(i*))/2, (chi_i-chi_(i*))/(2i),

and for each self-dual character use chi_i itself. These r functions form Phi. They descend to Q-regular class functions and are real on U. The change from the chi_i is invertible and linear overL. Hence the algebra of Q-regular class functions is Q[Phi_1,...,Phi_r]. This proves rationality of the coordinates without assuming that every fundamental representation is itself realized overQ.

Let d be the polynomial in these coordinates obtained as follows. In det(zI-Ad_g), take the coefficient of u^r after setting z=1+u. On a torus this is product_(alpha root)(1-alpha(a)), so the resulting regular class function determines a unique d in Q[Phi]. On the compact torus it is nonnegative and equals the squared absolute Weyl denominator. It is nonzero precisely at regular semisimple torus points.

## 3. The exact uniform character-domain law

Choose the compact torus A whose Lie algebra is spanned by the iH_i. It lies in the -1 eigenspace of omega, and omega(a)=a^-1 on A. This is a maximal torus of U and a maximal abelian subgroup in the symmetric directions. Consequently the restricted root system is the FULL root system and every root multiplicity is1.

Compact polar integration applies to U/(U^omega)^0. It also computes the Haar pushforward on X for the functions considered here, since they are invariant under the full fixed subgroup. For a=exp(iY), with alpha(Y) REAL root angles, its density is

    product_(alpha>0)|sin alpha(Y)|.

There is no suppressed complex-sine convention in this expression. Locally it follows directly: on the real root plane with basis E_alpha-E_-alpha and i(E_alpha+E_-alpha), Ad_(a^-1) rotates by alpha(Y); projection from the fixed direction onto the symmetric direction has determinant factor sin alpha(Y). The torus directions contribute1. Global surjectivity and the constant polar multiplicity are the compact Cartan/polar integration formula. The later even-multiplicity hypotheses of the cited Paley--Wiener paper are not used.

For h fixed by omega,

    c(haH)=h a^2 h^-1,

so every Phi_i(c(haH)) equals Phi_i(a^2). If Delta denotes the ordinary Weyl denominator, then

    |Delta(a^2)|=2^(number of positive roots)
                         product_(alpha>0)|sin alpha(Y)|.

Squaring pushes Haar measure on A to Haar measure on A. Thus the joint law of Phi(c) is proportional to the pushforward of |Delta(a)| da under Phi:A->R^r.

The fundamental-character Jacobian is exactly a nonzero constant times Delta. One short proof also explains its all-type scope. The Jacobian is Weyl anti-invariant. Its highest weight is omega_1+...+omega_r=rho, with coefficient1, and all other terms have lower weight. Expanding into alternating Weyl sums leaves only the rho term: a strictly dominant integral weight strictly below rho cannot occur. The Weyl denominator identity identifies this alternating sum with Delta. Replacing dual pairs by their real and imaginary parts changes the Jacobian by only a nonzero constant. This is also Steinberg Lemma8.2 and its complete proof, with the Haar-volume observation8.3.

Characters separate compact conjugacy classes; a regular class meets A in exactly one Weyl orbit of size|W|. The ordinary change-of-variables/area formula therefore cancels |Delta| against the character Jacobian and gives CONSTANT Lebesgue density on the character image. Normalize its mass to1. We obtain, for every polynomial P,

    lambda_R(P(Phi(c)))=vol(Omega)^-1 integral_Omega P(v)dv.  (3.1)

No averaging of a proposed pure phase has occurred; this computes the law of actual regular coordinates.

## 4. Domain boundary and an actual rational interior fiber

The compact image Omega_bar is closed and bounded. At a regular torus point the real character Jacobian is nonzero, so its image contains an open neighborhood. Thus Omega=int(Omega_bar) is nonempty. If a point of the image is on its boundary, no torus preimage can be regular; consequently d vanishes on the whole boundary. Regular torus points are dense, so Omega_bar is the closure of Omega. The boundary is contained in the zero set of a nonzero polynomial and has Lebesgue measure0. The image is semialgebraic, but no smooth-boundary or finite-perimeter theorem is needed below.

There is a rational point q in the nonempty open regular image, by density of Q^r. In fact Phi(c(U/(U^omega))) equals Phi(U): by polar decomposition every c-value is conjugate to a square in A, and every a in A is a square of another element of A. Thus q has an ACTUAL compact real fiber in X. The fiber over this rational q is a nonempty finite-type Q-scheme, so it has a Qbar-point. Translating gives an actual geometric common zero y=0. The polynomial f_q has f_q(0)=1, vanishes on the whole translated domain boundary, and has rational coefficients. These facts close the model hypotheses of the uniform-domain normalizer.

## 5. Complete multivariable endpoint argument and sharp reverse statements

The applicable general result is [the uniform bounded-domain criterion](UNIFORM_BOUNDED_DOMAIN_NORMALIZER_CRITERION.md), Sections 1--5. For clarity the analytic/algebraic calculation is included here. A bounded open domain and a polynomial f vanishing on its whole topological boundary suffice. If a C1 polynomial vector field V vanishes on that boundary, then integral_Omega divV=0: slice by each coordinate line, write the bounded open slice as countably many intervals, use the fundamental theorem at their boundary endpoints, and sum by absolute integrability and Fubini.

Work in the translated domain, still denoted Omega, with0 interior and f=f_q. Write

    J_k(T)=vol(Omega)^-1 integral_(T Omega) f(w)^(2k)dw
          =T^r vol(Omega)^-1 integral_Omega f(Tu)^(2k)du.

Initially T>0; the second expression defines J_k as a polynomial for formal T. Differentiation gives

    J_k'(T)=T^(r-1)/vol(Omega)
              integral_Omega div_u[u f(Tu)^(2k)]du.

In the coefficient of (T-1)^j for j<2k, the vector polynomial u f(Tu)^(2k) is divisible by f(u), so vanishes on the boundary. Its divergence integral is0 by the preceding zero-flux argument. Multiplication by T^(r-1) cannot lower the vanishing order. Hence

    J_k(1+t)-J_k(1) belongs to t^(2k+1) R[t].

Let mu_k=J_k(1)>0. Positivity follows on a neighborhood of0, where f0=1. Since all actual trace coefficients are rational, mu_k=lambda_Q(f(y)^(2k)) is rational. Coefficientwise substitution gives, for every k>=1 and integer ell>=-k,

    Lambda(t^ell p^k)=mu_k binom(k-r,k-ell).                (5.1)

Indeed Lambda(t^ell p^k) extracts coefficient t^(k-ell) in T^(k-r)J_k(T). The error starts in degree2k+1 even when k-r<0, because T is a unit in R[[t]]. The domain ell>=-k is essential. A negative lower binomial index means0.

For k>=r the value at ell0 is0 and the value at ellr is mu_k. Expanding F^m yields (1.2). All nonzero marked values are positive rational numbers, so remain nonzero in every characteristic-zero field.

For a<r the first trace of p^a is mu_a(-1)^a binom(r-1,a), nonzero; for a>=r all its pure moments vanish. For k>=r, the ordinary binomial in(5.1) is nonzero exactly when r<=ell<=k. This proves the unsafe and every-order marker classifications for p^a. For F, every sufficiently large summand k=rm+h has the same positive binomial sign whenever ell>=r and k>=ell; ell<r is eventually safe. At m1 only k=r,r+1 occur, so every-order survival requires ell<=r+1, and for r<=ell<=r+1 the h=m term survives at every m. Fixed negative ell eventually lie in the exact domain(5.1), which suffices for their eventual safety.

## 6. Actual nonpower, support and factor comparison

At a geometric point y=0,

    p=(1+t)/t, 1+p=2+t^-1.

At t=-1/2, p=-1 and partial_t(1+p)=-4. Smooth geometric integrality of X gives an actual prime divisor with valuation1 on1+p and valuation0 onp. Therefore F is not a unit times any proper power over an algebraic closure. Also

    F|_(y=0)=(1+t^-1)^r(2+t^-1),

with adjacent coefficients2 and2r+1 at exponents0 and-1. Actual Laurent support differences have gcd1. The same geometric argument works over every characteristic-zero extension, or may be replayed after embedding Qbar in Kbar.

Normalized trace/base change transfers all rational identities to every such K. Dense compact Haar gives strict real square positivity for every nonzero regular function on X. The finite rational Gram/LDL lemma in [Standard compact classical products](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md), Section 13.1, then gives square anisotropy over EVERY formally real K, with no embedding of the whole field in R. Only0 is pure on each individual X_K, so its trace kernel is Mathieu. The scalar Laurent constant-term kernel is Mathieu by the ordinary all-field Laurent theorem in the same source, Section 13.2. The pair(1.2) proves failure for their product. No complex-square positivity or unrestricted Mathieu ascent/descent is asserted.

## 7. A specified rational parameter and central-isogeny quotients

A free rational interior parameter is enough for Theorem1. There is also a finite root-theoretic choice that yields explicit invariance under every center quotient. In the split group choose the usual simple-root Weyl representatives

    n_i=x_alpha_i(1)x_-alpha_i(-1)x_alpha_i(1), n=n_1...n_r.

In the root SL2 this is the matrix with rows(0,1),(-1,0). It is fixed by the Chevalley involution and has rational split coefficients; therefore each n_i, and n, is a Q-point of the COMPACT descended group. It is in U, as a product of points of the connected compact rank-one subgroups.

Steinberg Theorem1.4/8.1 says that every element of the product of root cosets X_alpha_1 n_1 ... X_alpha_r n_r is regular in the simply connected group. In particular n is regular. As a compact group element it is semisimple. Thus

    q_*=Phi(n) in Q^r

is a regular interior point of the compact character domain. This gives a specified q in(1.1), using only finite root-group operations, fundamental characters and the adjoint discriminant polynomial. It does NOT assert that the entire Steinberg cross-section is defined over the compact Q-form. That would conflict with the Q-Borel obstruction in Section9 of Steinberg; only this single point is descended here.

Let Z be the geometric center. The image w of n in the Weyl group is a Coxeter element. The map a->a w(a)^-1 on the split maximal torus is surjective: on its character lattice 1-w is injective, as follows successively from the independence of the simple roots. This is Steinberg Lemma7.6 and its proof. For every z in Z choose a with a w(a)^-1=z. Then

    a n a^-1=z n.

The center acts linearly on fundamental-character coordinates, since each irreducible representation has a scalar central character. Denote the corresponding linear action in the real rational coordinates by R_z. The preceding conjugacy gives R_z q_*=q_*. Also Ad(zg)=Ad(g), so

    d(R_z v)=d(v).

This is a polynomial identity, because the character map is dominant. Consequently

    f_(q_*)(T R_z y)=f_(q_*)(T y)

as a polynomial identity in T,y. Thus p and F are invariant under central multiplication of the Cartan-map argument, although the individual y_i need not be invariant.

Now let pi:G->G'=G/C be ANY Q-defined central isogeny quotient, and let H'=(G')^omega be the FULL fixed subgroup. The involution descends because it inverts the center. If a in Gbar satisfies pi(a) in H', then omega(a)=a c_a for some c_a in C. Hence

    c(g a)=c_a^-1 c(g).

The just-proved center invariance makes every coefficient of p and F invariant under the full preimage of H'. They therefore descend to ACTUAL regular functions on X'=G'/H'. This also handles possible additional components of H' beyond the image of H. It is not a general descent theorem for arbitrary phases.

The induced dominant map X->X' has injective pullback on coordinate algebras, and normalized traces agree on its image by uniqueness. All pure/marked identities and sharp marker converses descend for this particular pair and seed. If the downstairs F were a unit times a proper power, its pullback would have that property upstairs; this contradicts Section6. Nonzero coefficient functions have the same Laurent support under the injective pullback, so support primitivity also descends. Compact square positivity and rational Gram give the individual formally-real factor comparison downstairs as well. All these statements commute with every characteristic-zero base extension of the specified Q-models.

## 8. The three exceptional split-type applications

The fixed Lie algebra of the Chevalley involution has one basis vector E_alpha-E_-alpha for every positive root, and no Cartan direction. Thus its dimensions in E6,E7,E8 are36,63,120. The compact symmetric tangent dimensions are respectively42,70,128, and the symmetric ranks are the full ranks6,7,8. This is the compact dual construction of the SPLIT real groups: the split Cartan decomposition has fixed directions E_alpha-E_-alpha and symmetric directions H_i,E_alpha+E_-alpha; multiplying the latter by i gives exactly the compact Q-basis in Section1. In particular the restricted multiplicities1 used above are proved by this construction, not inferred from a list or extrapolated from another real form.

The conventional real compact quotients, with their actual finite-center isotropy notation, are

    EI:    E6/(Sp(4)/{+1,-1}), rank6, dim42;
    EV:    E7/(SU(8)/{+I,-I}), rank7, dim70;
    EVIII: E8/Ss(16),          rank8, dim128.

Here the numerator is the compact simply connected group. Yokota's explicit full fixed-subgroup results identify these compact duals; the source details below distinguish the precise global groups from merely their Lie algebras. The Q-model in this theorem is the explicit Chevalley descent of Section1. No identification with every Q-form of those real homogeneous spaces is asserted.

For each of these three models choose any rational regular interior q, or the specified q_* of Section7, use the corresponding r=6,7,8 in(1.1), and use H_marker=t^r. The seed power cutoffs are EXACTLY6,7,8. For the nonpower phase F the every-order unsafe monomial markers are exactly t^r and t^(r+1), whereas every fixed exponent ell>=r is eventually unsafe and every ell<r is eventually safe. The cutoff is independent of the degree of the discriminant polynomial. Each pair is a real rational regular pair before arbitrary characteristic-zero base change. The three compact complex symmetric products already fall under known complex existence theorems; this proof does not claim new complex existence.

Together with other separately proved families these cases may support a later comparison across standard compact symmetric models. No such global classification is assumed or proved in this note.

### Semisimple product corollary

The same conclusions hold for a nonempty product of the specified simply connected simple Q-groups, with the involution acting factorwise, total rank r, and every Q-defined central quotient, including quotients by diagonal central subgroups. Indeed the actual character domain is the product of the component domains, its normalized pushforward is uniform product Lebesgue measure, and the adjoint discriminant is the product of the component discriminants. The bounded-domain theorem therefore applies with total rank r and gives precisely the same sharp cutoffs. Take n to be the product of the component compact rational Coxeter representatives. Then q_* is the tuple of the component parameters and 1-w is invertible on the direct sum of their cocharacter spaces. Every central element fixes this q_*, so the full fixed-subgroup descent in Section 7 applies without a componentwise restriction on the finite central subgroup. Smoothness, the actual common zero, absolute unit-nonpower, Laurent support and the formally real factor argument are unchanged. This is a factorwise corollary, not an assertion for arbitrary Q-forms or arbitrary symmetric involutions.

## 9. Sources and attribution

Thomas Branson, Gestur Olafsson and Angela Pasquale, *The Paley-Wiener Theorem and the Local Huygens' Principle for Compact Symmetric Spaces*, exact arXiv:math/0411383v1, is retained as [the polar-integration primary](sources/branson-olafsson-pasquale-math0411383v1.pdf), SHA256 `fc7352545c3aba0e3fcdeb02a9cf5efda4da75fa91c637ae450ef01b54808f7f`. The imported scope is whole Section 1, printed 3--5, especially equation 1.1 and Lemma 1.1 with its polar-Jacobian justification. This section has no even-multiplicity restriction. The later Paley--Wiener and shift-operator assertions require even multiplicities and are excluded. Section 3 above supplies the unambiguous real root-angle convention.

Robert Steinberg, *Regular elements of semi-simple algebraic groups*, Publ. Math. IHES 25 (1965), 49--80, is retained as [the character primary](sources/steinberg-regular-elements-1965.pdf), SHA256 `bda20c4a4b3cf4a7e9e5b17a4eacd0147d0c9370063ab3a3778114a53c8349b9`. The character-algebra and semisimple-separation scopes are Theorem 6.1, Lemmas 6.2--6.4 and their proof, and Corollary 6.6 with proof, printed 62--63. The Jacobian scope is Lemma 8.2 and its entire proof, with the volume-form observation 8.3, printed 70--71. The optional specified q_* and central-invariance branch use Lemma 7.6 with proof, printed 67, and Theorem 1.4/8.1 with its Section 8 argument, printed 70--71. Earlier structural lemmas cited there remain classical imports. Only the single compact rational point n is descended; no rational compact-form Steinberg cross-section is asserted, respecting the Section 9.1 obstruction.

Ichiro Yokota, *Exceptional Lie groups*, exact arXiv:0902.0431v1, is reused as [the global-model comparison](sources/yokota-exceptional-lie-groups-0902.0431v1.pdf), SHA256 `a8226e37e80cc203113f2fc0bb17ccbe0773cfa12918252f8e3bc06f4436ecc9`. The relevant read scopes are entire Section 3.12, printed 95--100, including Lemma 3.12.1 and the complete Theorem 3.12.2 surjectivity/kernel proof; entire Section 4.12, printed 142--146, including Proposition 4.12.3 and Theorem 4.12.5 with both surjectivity arguments; and entire Section 5.8, printed 176--184, including the so16 bracket calculation and global Ss16 argument. The latter retains its cited representation/irreducibility inputs. Split-real polar identifications 3.15.1, 4.15.1 and 5.13.1, printed 105, 156 and 199, give the compact-dual nomenclature. This is a comparison of conventional real global models, not a whole-book classification audit or a classification of all Q-forms. The all-type multiplicity-one calculation and actual coordinate law are proved independently in Sections 2--4 above.

The normalized Reynolds projection, same-field base change, compact square positivity, rational Gram/LDL and scalar Laurent Mathieu theorem are used only in the generic scopes of [Standard compact classical products](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md), Sections 3 and 13. Matsushima affineness and the supporting characteristic-zero group facts keep the precise existing sources recorded in [the coadjoint proof](PIECEWISE_DENSITY_AND_COMPACT_COADJOINT_ALL_ORDER_COUNTEREXAMPLES.md), Section 12. No full-automorphism or Duistermaat--Heckman theorem is a premise of the present construction.

The closest identified normalizer antecedent is Christopher D. Long's spherical draft dated 22 September 2026, exact repository revision `e9edac262b4b64eb68feb0e6e80496cefe4a1dee`, retained [PDF](sources/long-spherical-e9edac2.pdf), SHA256 `b2bf25fbc2d6211333d4ff3e1f6de0b6562990e869c6078a2e219ef746f1026b`. Whole Section 5 and the complete Proposition 5.1 proof, printed 5--6, give the known one-dimensional interval normalizer and high-order endpoint cancellation, including its marker tower. Theorem 1.1 already gives complex pure/unsafe existence on these nontrivial spherical products. The multivariable boundary cancellation and actual rational character-domain realization are proved above and in the separate bounded-domain criterion. No worldwide novelty follows from this comparison or from a search not locating an equivalent formulation. The classical normalizer, polar and character mechanisms retain their attribution.
