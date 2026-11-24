# Compact rational homogeneous spaces with a nonabelian submersion over zero

The theorem gives a sufficient class of actual rational homogeneous spaces. A single submersion point over zero suffices for an all-order pair; it does not make zero a regular value of the entire moment map. Individual Cartan coordinates are not asserted to be functions on the homogeneous source.

## 1. The actual-model theorem

Let G/Q be connected semisimple with compact G(R), and let H be a closed Q-subgroup with h=Lie(H^0) of positive dimension. Take B=-Killing_g and suppose there is an actual nonzero rational vector

    x in h(Q)^perp,              h intersect C_g(x)=0.          (1)

Let M=Ad(G(R)^0)x have dimension 2D, put r=rank H^0, and let R_H be the number of positive roots of H^0. Then G/H is an actual smooth geometrically integral affine Q-variety. Its normalized rational trace lambda is compatible with every characteristic-zero scalar extension. For Lambda=lambda composed with Laurent constant term, there are Q-defined functions on (G/H) times split G_m of the form

    a(T)=theta(f(Tu)),       T=1+t,
    p=T*t^-1*a(T)^2,        F=p^L*(1+p),        H_marker=t^e,    (2)

where f is a rational invariant polynomial and L is sufficiently large, such that for every m>=1

    Lambda(F^m)=0,           Lambda(H_marker F^m)>0.            (3)

All marked values in (3) are positive rational numbers. The exponent satisfies

    r<=e<=r+2R_H=dim H^0.                                    (4)

The same phase and fixed scalar Laurent unit have nonzero marked moments at EVERY positive order after EVERY characteristic-zero extension of Q. The phase is not a unit times a proper power over the geometric coefficient field and has primitive actual Laurent support.

Over every formally real extension each separate factor, G/H and split G_m, has a Mathieu trace kernel, while the product fails. The all-field pair and the formally-real factor statements have different quantifiers.

If zero is separately assumed to be a GLOBAL regular value of the moment map Phi below, the marker can be chosen explicitly as

    H_marker=t^(dim H^0).                                    (5)

Hypothesis (1) alone gives one submersion point over zero, not the global hypothesis. Neither version asserts a least unsafe exponent, a smallest L or failure for every proper subgroup.

## 2. Geometry, native trace and the actual zero point

The positive definite rational form B makes ad(x) skew after extension to R, hence semisimple over C. The closed-coadjoint-orbit argument recorded in [UCQ, Sections 2--3] supplies the actual smooth geometrically integral affine orbit X=Ad(G)x, its rational full connected reductive centralizer, and its native normalized trace. That argument applies to arbitrary compact rational forms and nonperiodic rational x; no coweight or circle assumption is added here.

The identity component H^0 is reductive. A positive-dimensional real unipotent radical would have unbounded real points inside compact G(R), impossible in characteristic zero. Matsushima's affine direction over an algebraic closure, followed by affine descent, gives G/H^0 affine. The finite component group H/H^0 acts freely on this quotient; its quotient is the represented G/H and is affine. Quotient torsors in characteristic zero make G/H smooth; its geometric irreducibility follows from surjectivity of the connected geometric G. The identity coset eH is a rational point. These assertions also hold after every scalar extension. Only the affine direction of Matsushima is needed.

Complete reducibility gives the normalized rational Reynolds trace on O(G/H), unique among G-invariant functionals with lambda(1)=1. It commutes with field extension and equals the restriction of the normalized trace on O(G). Over R the latter is compact Haar integration. Density of G(R)^0 and invariance identify it with integration over that connected compact group; on X it pushes to normalized KKS/Liouville probability on M. These are the generic normalized-trace and quotient facts of [RHT], [HQP] and [UCQ], not any earlier counterexample phase.

Projection pr_h for B is Q-linear and equivariant under FULL H, since H normalizes its identity-component Lie algebra and preserves B. On M define

    Phi(z)=pr_h(z).                                          (6)

It is the Hamiltonian moment map for K_H=H(R)^0, using B to identify h with h*. Equation (1) gives Phi(x)=0. At x its differential has annihilator h intersect C_g(x), because tangent vectors are [a,x] and invariance of B identifies the annihilation condition with [x,v]=0. Therefore dPhi_x is surjective.

## 3. Absolute continuity and the local lower bound

The maximal minors of dPhi, in algebraic orbit coordinates, are real analytic. At least one does not vanish identically because of the preceding submersion point. The critical locus is a proper real-algebraic subset of the connected analytic orbit, hence has Liouville measure zero. On its complement, countably many submersion charts and ordinary Fubini show that Phi_* of Liouville probability is absolutely continuous on h_R.

A small product chart above x has the form (v,w) with Phi(v,w)=v. Restrict w to a smaller compact box and v to a ball about zero. The positive smooth source density has a positive lower bound there. Integrating over that box proves that the COMPLETE pushed probability has density at least c>0 almost everywhere on that ball. Other parts of the zero fiber may be critical and can add singular density behavior. This proves a lower bound, not global smoothness at zero.

The measure is K_H-invariant. Its density can be replaced by its compact-group average, retaining an invariant almost-everywhere representative and the same local lower bound on a smaller invariant ball.

## 4. Abelian polynomial law and finite algebraic walls

Choose a maximal Q-torus T_H of H^0 and a maximal Q-torus T_G of G containing it. For example choose a maximal torus in the connected reductive centralizer of T_H; it is maximal in G. The compact torus T_H(R)^0 has dimension r. Its action on M has no infinitesimal kernel because the H-action is locally free at x. Dividing a possible finite kernel makes the action effective without changing the real moment coordinates except for lattice normalization.

The Abelian moment map is pr_tH. The original Duistermaat--Heckman Section 3 theorem [DH] therefore gives an absolutely continuous probability law, polynomial of degree at most D-r on regular-value components. Only the variation/reduced-volume theorem is used; no statement from original Section 4 is imported.

Here is a finite algebraic hyperplane arrangement sufficient to refine those components, even when x is not regular in G. Expand a point z in the weight decomposition of g_C for T_H. Its stabilizer Lie algebra in t_H is an intersection of the kernels of the finitely many nonzero restricted roots that occur in z. Thus only finitely many stabilizer Lie spaces can occur. They are defined by torus characters and have real algebraic bases in any rational coordinate system.

For each nonzero possible stabilizer space a choose a nonzero real algebraic vector h_a in it. A point fixed by the corresponding connected subtorus A can be put into t_G by conjugation in C_G(A). Indeed C_G(z) is connected reductive by the full-centralizer argument in [UCQ]. A lies in it; extend A to a maximal torus T_z of C_G(z). The semisimple element z is central in this Lie algebra and belongs to Lie(T_z). A semisimple-element centralizer has full G-rank, so T_z is also a maximal torus of G. It lies in C_G(A), and compact maximal-torus conjugacy there carries T_z to T_G. This conjugation fixes h_a and preserves the h_a moment component.

The orbit intersection X_C intersect t_G,C is the finite Weyl orbit of a representative x_t. It is a finite scheme defined over Q, so its geometric points are algebraic. Real points lie in the compact real orbit as needed. Every critical value is consequently contained in a hyperplane

    B(v,h_a)=B(w x_t,h_a).                                   (7)

There are finitely many, with real algebraic coefficients. Their constants MAY vanish. Every connected component of the complement is an open convex polyhedron. On it the Abelian density is a polynomial of degree at most D-r. If the chamber is unbounded, compact support forces this polynomial to vanish on a nonempty open tail and hence identically. Bounded chambers suffice for every nonzero density piece. This argument uses no regular-G-orbit hypothesis and no genericity that would remove through-zero walls.

## 5. Exact nonabelian derivative scope and absence of wall mass

The nonregular elements of the compact reductive Lie algebra h_R form the zero set of a nonzero invariant discriminant polynomial; when H^0 is a torus the nonregular set is empty. They have Lebesgue measure zero. Section 3 therefore proves that the actual invariant pushed measure nu=Phi_*probability satisfies

    nu(h_R minus Ad(K_H)t_H,regular)=0.                       (8)

This is the strong support hypothesis in [CDKW, Theorem 3.2], not a consequence asserted for arbitrary Hamiltonian measures.

Let tau send a regular element to the chosen open positive Weyl chamber. CDKW definitions (1.1)--(1.3), Theorem 3.2 and Corollary 3.3 give, with consistent measure normalizations,

    tau_*nu = p_H(v) * product_(alpha>0) partial_(-alpha)
                         (pr_tH)_*nu                        (9)

on that OPEN chamber. The factor p_H is the orbit-volume polynomial, linear in each positive root and of degree R_H. It is NOT the squared discriminant from ordinary Lie-algebra Lebesgue integration.

On every chamber refining the Abelian hyperplanes (7) and the H-root hyperplanes, (9) gives a polynomial density. Differentiation lowers degree by R_H and multiplication by p_H raises it by R_H, so the degree is still at most D-r.

The radial probability tau_*nu is itself absolutely continuous on t_H: nu is absolutely continuous on h, the radial quotient is a smooth submersion on regular h, and the singular locus has nu-mass zero. The finitely many critical and root hyperplanes thus have zero radial mass. Off-wall polynomial densities determine the complete measure almost everywhere. No singular wall terms need be discarded by a formal computation. In particular, the warning in CDKW Remark 3.6 about positive critical-wall mass is respected by a separate proof of absolute continuity in this application.

Extend the radial probability by zero outside the chosen positive chamber. It now has a finite polyhedral polynomial decomposition in R^r with nonnegative density pieces, algebraic supporting hyperplanes and degree bound b=D-r. All nonzero pieces are bounded by compact support. There is no claim that this density is globally continuous.

## 6. The squared Lie-algebra Jacobian and the marker range

For completeness the precise ordinary Weyl Jacobian needed here has a local proof. Let K be any compact connected reductive group, T a maximal torus, and t_reg its regular Lie algebra set. The map

    psi:(K/T) times t_reg -> k_reg,
         ([g],v) -> Ad(g)v                                  (10)

is a |W_0|-sheeted local diffeomorphism. Every regular element has centralizer T, so its preimages correspond to the Weyl group. At (eT,v) the differential is

    (a mod t,w) -> [a,v]+w.

On the real root plane for alpha, ad(v)=alpha_R(v)J, where J is orthogonal and J^2=-1. Its determinant in absolute value is alpha_R(v)^2; the t summand has determinant 1. Thus the Jacobian is

    Delta_H(v)^2=product_(alpha>0)alpha_R(v)^2.                (11)

Ordinary change of variables on local charts of the finite covering, followed by compact averaging, proves that an invariant absolutely continuous h-measure with density a has actual radial probability density

    C*Delta_H(v)^2*a(v)                                      (12)

on the positive chamber, for a fixed C>0 recording normalization. Nonregular sets have measure zero, so they do not affect the formula. This proves the required Weyl integration assertion directly and distinguishes it from the LINEAR p_H in (9).

The local lower bound in Section 3 implies

    rho(v)>=c'*Delta_H(v)^2                                 (13)

almost everywhere near zero in the positive chamber. Hence some nonzero full-dimensional cell has zero in its closure. Write its and the other such cells' polynomials into homogeneous pieces, and let j0 be the least degree occurring among them. Along an open cone the polynomial lower bound (13) forces j0<=2R_H. Nonnegativity forces each first local homogeneous term to be nonnegative as in [RIT]. Therefore

    e=r+j0 lies in [r,r+2R_H]=[rank H^0,dim H^0].              (14)

If zero is a GLOBAL regular value, compactness makes the moment map a proper submersion above a sufficiently small ball. Integration in submersion charts gives a smooth positive invariant h-density there. Formula (12) has first homogeneous term C*a(0)*Delta_H^2 with nonzero coefficient, so j0=2R_H. A single submersion point over zero does not justify this last step.

## 7. Rational Chevalley restriction and full components

The exact invariant-polynomial input is [Etingof, Theorem 10.1 and Remark 10.2(2)]: restriction from a complex reductive Lie algebra to a Cartan is a graded algebra isomorphism onto Weyl invariants. The complete published lecture proof first uses dense semisimple elements for injectivity, then triangular inversion of highest-weight characters into Weyl orbit sums. Powers of dominant integral weights span each symmetric power by their Zariski density, giving surjectivity. The preceding highest-weight and semisimple-density facts remain classical inputs; no Chevalley--Shephard--Todd theorem is used.

Apply this theorem to h_C and t_H,C. In each degree the invariant spaces are kernels of rational linear equations. Connected characteristic-zero group invariance agrees with Lie invariance; a vector with zero differential action has stabilizer containing the connected group. Kernels commute with field extension. The rational restriction map becomes an isomorphism after C-extension, so it already is a graded Q-isomorphism

    Q[h]^(H^0) -> Q[t_H]^(W_0).                              (15)

This descent of the restriction map does not require rational representations affording every highest weight.

The normalizer of T_H in FULL H has identity component T_H. Its quotient by T_H is a finite etale Q-group, possibly with kernel in its action on t_H. Let W_full be its finite geometric matrix image. It is Galois-stable and preserves the positive rational form B|t_H.

Every geometric component element h of H carries T_H to another maximal torus of H^0. Multiplication by a suitable element of H^0 adjusts it to normalize T_H. For an H^0-invariant polynomial this adjustment does not change the invariance question. Chevalley injectivity now proves the graded full-component isomorphism

    Q[h]^H -> Q[t_H]^(W_full).                              (16)

No rational representative of a geometric component is assumed. Finite Galois-stable matrix image, rather than an arbitrary abstract finite algebraic matrix set, is the rational symmetry used below.

## 8. Actual invariant-test functions, origin evaluation and all-order conclusion

For P in the right-hand side of (16), let P_tilde be its unique full-H-invariant polynomial lift. Define

    theta(P)(gH)=P_tilde(pr_h(Ad(g^-1)x)).                    (17)

This is an actual unital Q-algebra homomorphism to O(G/H). Projection is full-H-equivariant and P_tilde is full-H-invariant, so (17) descends regularly from G by its quotient torsor. Individual Cartan coordinates are not asserted to descend.

Inverse Haar measure on G(R)^0 has the same law as Haar. The native trace of (17) is consequently the integral of P_tilde against the actual h-moment law nu. Chevalley restriction then identifies it with the integral of P against the radial density of Section 5. Thus EVERY hypothesis of the invariant-test conditional theorem [RIT] holds, with the actual algebra O(G/H), W_full, b=D-r, the native scalar-compatible trace and the density just proved.

The finite W_full/Galois orbit product supplies the rational polynomial killing every NONRADIAL supporting hyperplane; through-zero faces have zero Euler flux. Its invariant positive quadratic damping gives the f of [RIT]. Choose L>=max(D,K) as there and apply its complete coefficient and positivity proof. Equations (2)--(4) follow. Their full marked values are positive rational numbers by native trace, not by a claim that every density coefficient is rational.

At the ACTUAL rational identity coset,

    theta(P)(eH)=P(0),                                      (18)

because pr_h(x)=0. The geometric conclusions therefore hold DIRECTLY downstairs by [RIT, Section 7]: p=1+t^-1 at that point, the modifier has simple root t=-1/2 with derivative -4 and p=-1, giving an actual valuation-one prime divisor; coefficients 2 and 2L+1 at exponents 0,-1 give primitive Laurent support. Neither an arbitrary nonpower pullback principle nor an imaginary tuple of Cartan coordinates is used.

## 9. Formally-real factors and field quantifiers

Every nonzero real regular function on G/H pulls back to a nonzero function on G. The compact real identity component is Zariski dense, so this function is nonzero on a positive-Haar-measure open set. Its square has strictly positive trace. Every finite rational trace Gram matrix is thus positive definite over R. Rational LDL elimination has positive rational pivots; each pivot a/b is the sum of ab copies of (1/b)^2 for positive integers a,b. Over ANY formally real extension, such a finite sum-of-squares form is anisotropic. Each element uses only a finite rational span, hence lambda_K(a^2)=0 forces a=0. Only zero is pure on the compact factor, so its trace kernel is Mathieu.

The split scalar Laurent factor is independently Mathieu over every characteristic-zero field by the scalar Duistermaat--van der Kallen theorem and its coefficient-field argument [SCCP, Section 13.2]. It is not an algebra-valued constant-term theorem. Over formally real fields both individual factors are therefore Mathieu, but their product fails by (3).

All-field statements concern scalar extensions of THIS actual Q-model, its native trace and its actual rational phase. No arbitrary twisted form is introduced. No global Mathieu ascent theorem is inferred. The phase has the stated divisor and support properties after every extension by its direct actual certificate, independently of the positivity of that extension field.

## 10. Classical complex antecedent and elementary unitization

A useful comparison is elementary. Suppose a commutative characteristic-zero algebra has a linear functional ell and an ALREADY GIVEN pair f,h with ell(f^m)=0 and ell(h f^m)!=0 for every m>=1. On its Laurent extension put

    P=f*(1+t*h),                 U=t^-1.

The binomial expansion immediately gives

    Lambda(P^m)=0,               Lambda(U P^m)=m*ell(h f^m).   (19)

If the algebra is the coordinate ring of a smooth geometrically integral affine variety, the geometric open D(fh) is nonempty. At a point there, t=-1/h is nonzero and 1+t*h has nonzero t-derivative h while f is a unit. Thus P already has an actual order-one prime divisor and cannot be a unit times a proper power. Its nonzero coefficients f and fh at adjacent exponents 0,1 give primitive Laurent support.

Liu's deposited preprint version 2.2 [Liu], Theorem 3.1 and Corollary 3.2, gives complex finite-type compact-group functions, under its stated independent adjacent-weight hypotheses,

    f=(1+Bbar)*(B-A*Abar),       h=Bbar,
    integral f^m=0,              integral h f^m=1/dim V_(m lambda).

Its coherent-state proof uses the Cartan tensor component, the adjacent root-string vector and Schur orthogonality; the factorial radial circuit is proved there. Equation (19) therefore already gives the complex Laurent-unit and geometric-strength conclusions on the corresponding affine complex group product. These conclusions alone are not a new complex-existence claim.

That input does not supply a native Q-defined pure f on an arbitrary specified compact Q-form. Indeed such a nonzero pure f on the compact factor would contradict strict positivity of its real square trace. The present theorem instead constructs the rational phase on the product by an actual invariant-test moment law, and obtains a witness over every characteristic-zero extension of that Q-model. Genuine instability of two separately Mathieu factors over formally real fields is also a separate conclusion. No priority assertion follows from this comparison.

## 11. A finite rational input gate and split-type involutions

The geometric input (1) has a finite exact algebraic test. Put d=dim h>=1 and p_orth=h^perp. Invariance of B gives [h,p_orth] contained in p_orth. Choose rational bases a_1,...,a_d of h and v_1,...,v_n of p_orth. For x(c)=sum_i c_i v_i form the n-by-d rational matrix

    M(c):h->p_orth,             a -> [a,x(c)].                 (20)

Its entries are homogeneous linear polynomials in c, and the required zero Lie stabilizer is EXACTLY rank M(c)=d. Thus the following are equivalent:

1. Some x in p_orth,C has h_C intersect C_gC(x)=0.
2. Some maximal d-minor of M(c) is a nonzero rational polynomial.
3. Some point c in the finite rational grid {0,1,...,d}^n has rank M(c)=d.
4. Some rational x in p_orth has the actual submersion-zero property (1).

Only 2=>3 needs more than the rank criterion. A d-minor has total degree d, hence degree at most d in each variable. A polynomial with these separate degree bounds vanishing on the entire grid is zero, by univariate interpolation in one variable followed by induction on the others. Therefore a nonzero minor is nonzero at a grid point. If n<d, no maximal minor exists and rank d is impossible; the test fails correctly. Since d>0, a successful x is automatically nonzero.

This finite test certifies an INPUT rank condition, not infinitely many moments. It may require up to (d+1)^n rational rank computations and makes no efficiency or analytic-amplifier bound claim. The nonempty full-rank locus is invariant under complex isomorphism of the pair (g,h), since the Killing orthogonal complement and isotropy representation transport with that isomorphism. Consequently every ACTUAL compact rational form of a given complex pair with an infinitesimally generically free isotropy representation passes this rational test. A label or an unverified branching table alone does not identify that complex pair.

There is an intrinsic useful sufficient case without a type classification. Suppose the actual compact Q-group has a Q-involution theta_0, let h=g^(theta_0), and suppose its minus eigenspace p_minus,C contains a Cartan subalgebra c of g_C. B-invariance of theta_0 gives p_minus=h^perp. Choose a complex regular element z in c. Its Lie centralizer is c, contained in p_minus,C, so h_C intersect C_gC(z)=0. The finite-grid argument gives an ACTUAL rational x satisfying (1). Hence the theorem applies to every FULL closed Q-subgroup H with this Lie algebra h, under the stated positive-dimension assumption. In particular the conclusion is not restricted to a standard Chevalley Q-model of that complex split-type pair.

This corollary imports no list of symmetric types and gives no global regular-zero conclusion. Failure of the finite test is not a converse obstruction to arbitrary Mathieu counterexamples. It tests precisely the sufficient-class input used in Sections 2--8.

## 12. Exact sources and limits of the imports

[CDKW] M. Christandl, B. Doran, S. Kousidis and M. Walter, Eigenvalue Distributions of Reduced Density Matrices, exact arXiv:1204.0741v2, retained [primary PDF](sources/christandl-doran-kousidis-walter-1204.0741v2.pdf), SHA256 88a81ab03ab006491e2cfdad7dc5966125730b1f126ca3e61c9b45422e5a5683. Actual read scope: Definitions (1.1)--(1.3), Assumption 1.1 and its discussion on printed 5--6; Section 3 printed 15--17, including Theorem 3.2 statement and COMPLETE proof, Corollary 3.3 COMPLETE proof and Remark 3.6. The imported theorem requires the measure to give zero mass to nonregular coadjoint orbits, as proved here in Section 5. Its proof uses the earlier maximal-orbit Fourier formula (3.1)--(3.3) and an imported Harish-Chandra orbital formula; that historical proof is not independently audited here. No complete-paper certification is claimed. The previously retained Wolf use of this same physical PDF was context; the present derivative use has its own mathematical hypotheses.

[Etingof] P. Etingof, MIT 18.757, Fall 2023, Lecture 10, [official source](https://ocw.mit.edu/courses/18-757-representations-of-lie-groups-fall-2023/mit18_757_f23_lec10.pdf), retained [PDF](sources/etingof-chevalley-restriction-lecture10-2023.pdf), SHA256 629d82eeb9649aee60c6ad23558292c2465cc5e52790d8dec7163cfcb065c417. Actual scope is the ENTIRE Theorem 10.1 statement and proof on printed 54--55/PDF 1--2 plus Remark 10.2(2). Highest-weight existence/triangularity and semisimple-density facts from earlier representation theory remain classical inputs. Neither the rest of the lecture course nor Chevalley--Shephard--Todd is imported. The rational/full-component extension and the squared Weyl Jacobian are proved here.

[DH] J. J. Duistermaat and G. J. Heckman, On the variation in the cohomology of the symplectic form of the reduced phase space, Invent. Math. 69 (1982), 259--268, [retained PDF](sources/duistermaat-heckman-1982-original.pdf), SHA256 e314c2033972ab39d244b8bc98ddfb37d8ade96c8c36bfa95e05d82fc1033be5. The retained exact scope is the variation theorem and complete Section 3 reduced-volume/polynomial-density argument, printed 260--266. Finite isotropy/orbifold reductions are included; no global-freeness assumption is made. No original Section 4 statement is used.

[RHT] [Normalized reductive homogeneous trace](REDUCTIVE_HOMOGENEOUS_AUTOMORPHISM_INVARIANT_TRACE.md), Sections 2 and 14, provides only affine-quotient/native-Reynolds geometry used here. Its Matsushima input is Arzhantsev, Invariant ideals and Matsushima's criterion, arXiv:math/0506430v1, Section 1 and Theorem 2.1/Proposition 2.2; the affine direction descends from an algebraic closure. Quotient torsors and complete reducibility use the retained Milne Algebraic Groups (2017) scopes 5.25--5.28, 12l/12.53 and 22.42--22.43. No all-automorphism theorem is needed.

[UCQ] [Universal compact rational coadjoint products](UNIVERSAL_COMPACT_RATIONAL_COADJOINT_PRODUCT_COUNTEREXAMPLES.md), Sections 2--3, supplies the generic arbitrary-compact-Q semisimple closed-orbit/full-centralizer and Haar/Liouville trace lemmas. Its generic rational torus/centralizer inputs are the retained Milne 17.37--17.40, 17.59, 17.76, 17.82--17.83 and 17.90, with the stated proof scopes. Neither its scalar-height density theorem nor its earlier phase is a premise here.

[HQP] [Homogeneous quotient traces](HOMOGENEOUS_QUOTIENT_TRACE_AND_COUNTEREXAMPLE_PULLBACK.md), Sections 2--3 and 7, records generic native trace restriction, finite-component reductivity and compact rational Gram/SOS factor arguments. [SCCP] [Standard compact classical products](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md), Sections 3 and 13, records the same generic trace/Gram lemmas and the scalar Laurent theorem, from the retained Duistermaat--van der Kallen PDF. No seven-family model result is required.

[RIT] [Radial invariant-test criterion](RADIAL_WALL_INVARIANT_TEST_NORMALIZER_CRITERION.md) is the separate complete analytic theorem actually applied in Section 8.

[Liu] Z. Liu, Coherent-State Circuit Cancellation and the Compact-Group Mathieu Conjecture, deposited preprint version 2.2, [primary record](https://zenodo.org/records/21689760), DOI 10.5281/zenodo.21689760. The comparison uses Theorem 3.1 and Corollary 3.2 within their stated highest-weight hypotheses. Retained comparison PDF SHA256 132a9bea591006d5cf2fd3e55524c8bff76fcf41e29cc58f20eb421da809bcee. This is a comparison, not a premise of (1)--(18), not a native-Q-form theorem and not a certification of every auxiliary classification in that preprint.

There is no claim for all proper H. The proof supplies no lower-dimensional radial-law criterion, no phase descent by arbitrary averaging and no conclusion from a finite prefix. Finite stabilizers have rank zero and are outside this theorem. A subgroup may fail (1), even when some other counterexample construction exists.
