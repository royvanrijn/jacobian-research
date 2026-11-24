# Compact rational torus normalizer product counterexamples

## 1. The model, joint moment law and counterexamples

Let G be a connected semisimple algebraic group over Q with compact real points, and let S be ANY positive-dimensional Q-torus in G. Put s=Lie(S), r=dimS, B=-Killing_g, and K=G(R)^0. Then there is a maximal Q-torus T containing S and a regular y in t(Q), t=Lie(T), with the following properties.

Write M=Ad(K)y, dimM=2D, and identify s_R with its dual by B. The actual moment map

    nu:M -> s_R,       nu(z)=pr_s(z)

has 0 as a NONEMPTY regular value. Its normalized Liouville pushforward has a nonnegative probability density rho, polynomial of degree at most D-r on every chamber of a finite arrangement of REAL ALGEBRAIC affine hyperplanes. None of the hyperplanes contains0. The central chamber is bounded and rho(0)>0. Every unbounded chamber has identically zero density. Thus the measure is a finite sum of polynomial densities on bounded full-dimensional open convex polyhedra, with every supporting hyperplane avoiding0. No global continuity across the walls is needed for this assertion.

On the ACTUAL affine Q-variety Y=G/C_G(S), the Q-regular tuple

    mu(g C_G(S))=pr_s(Ad(g^-1)y)

has exactly this joint polynomial functional under its native normalized rational trace. In particular mu=0 has an ACTUAL real point and a geometric Qbar-point. All polynomial moments are rational and commute with characteristic-zero scalar extension.

Finally N=N_G(S) has N^0=C_G(S), and W=N/C_G(S) is finite etale overQ, with a faithful Q-linear action on s. The actual quotient Z=G/N is finite-etale covered by Y. There exists an ACTUAL rational polynomial f on s, with f(0)=1, invariant under the FULL geometric W action, vanishing on every supporting hyperplane of the above bounded cells. This statement includes nonconstant finite-Q-group schemes and does not presume rational normalizer representatives.

**Theorem 1 (torus-normalizer product).** For every G and S above, there are a W-invariant rational polynomial f with f(0)=1 and an integer L>=D such that the ACTUAL functions

    T=1+t, p=T*t^-1*f(Tmu)^2, F=p^L(1+p), H_marker=t^r

are regular on Z times splitGm and satisfy Lambda(F^m)=0 and Lambda(H_marker F^m)>0 inQ for EVERY m>=1. The same fixed UNIT marker and the same Q-defined phase work after EVERY characteristic-zero extension. The phase is absolutely not a unit times a proper power and has primitive Laurent support. The integer r=dimS is the LEAST unsafe integer Laurent-monomial marker exponent for this phase. Over every formally real field, the quotient factor and split torus factor separately have Mathieu kernels, while their product does not.

The expression f(Tmu) denotes its ACTUAL descended invariant coefficients. Individual mu coordinates need not descend from Y to Z. The all-order analytic input is [the polyhedral density criterion](POLYHEDRAL_PIECEWISE_DENSITY_NORMALIZER_CRITERION.md); all its geometric hypotheses are proved here. There is no sharp minimal-L or all-larger-marker classification, phase averaging, or arbitrary-form inference.

## 2. Rational torus, projection and the actual native trace

C=C_G(S) is connected reductive by the torus-centralizer theorem. Choose a maximal Q-torus T of C. Such a torus exists overQ and remains maximal after extension, by Milne17.82--17.83. The central torus S of C is contained in T: otherwise ST would be a strictly larger torus in C. Any torus of G containing T also centralizes S and hence lies in C, so T is maximal in G.

For an explicit affine-model proof, choose x0 in s(Q) outside the kernels of every NONZERO restricted root. Then Lie(C_G(x0))=Lie(C_G(S)) by the root decomposition. There is an inclusion C_G(S) contained in C_G(x0), and both FULL centralizers are connected: the latter follows from the already proved compact semisimple-element/exponential-torus argument. Equal Lie algebras therefore give C_G(x0)=C_G(S). The existing arbitrary compact-Q coadjoint closed-orbit proof identifies G/C_G(S) with the ACTUAL closed affine orbit Ad(G)x0. In particular affineness here needs no new Matsushima premise; x0 is a rational Lie element, not necessarily a coweight.

The form B is Q-valued, invariant and positive definite on g_R. Fix a Q-basis e_i of s with Gram matrix B_s. Then

    pr_s(v)=sum_(i,j)e_i (B_s^-1)_ij B(e_j,v)

is an ACTUAL Q-linear map. For n in N_G(S), invariance of B and preservation of s show that pr_s commutes with Ad(n). Consequently

    mu(gn)=Ad(n)^-1 mu(g),

and mu(gc)=mu(g) for every c in C, scheme-wise. Thus all coordinates of mu descend to Q-regular functions on the actual affine quotient Y=G/C. No rational basis of integral real circle generators is asserted.

The normalized Reynolds trace on O(Y) is the restriction of the normalized G-trace under pullback to O(G). OverR it is integration on the dense compact orbit K/(K intersect C(R)). Inversion preserves Haar probability on K. Hence for EVERY polynomial P,

    lambda_Y(P(mu))
       =integral_K P(pr_s Ad(g^-1)y)dg
       =integral_K P(pr_s Ad(g)y)dg
       =integral_M P(nu(z))dnu_M(z).                    (2.1)

Here dnu_M is normalized K-invariant probability. It equals the positive normalized KKS Liouville measure. The actual coordinate functions and trace are rational, so their polynomial moments belong toQ; this is not an assertion that the density coefficients separately are rational.

The orbit M and its affine geometric Q-orbit exist by semisimplicity of y. Compactness makes every rational Lie element semisimple, and regularity gives centralizer T. The actual closed-orbit and normalized-trace arguments are proved in [the arbitrary compact-Q coadjoint theorem](UNIVERSAL_COMPACT_RATIONAL_COADJOINT_PRODUCT_COUNTEREXAMPLES.md), Sections2--3; only those model/trace scopes are used here, not its all-order phase construction. No compact real label alone is used to specify an unrelated Q-form.

## 3. A finite list of all possible isotropy Lie spaces

Work first overR and then complexify. Let R_G be the roots of (g_C,t_C). Their restrictions to s_C are the nonzero S-weights of g_C; the remaining summand is the zero-weight Lie algebra of C. If z=sum_lambda z_lambda in that finite weight decomposition, then

    s_z={h in s_R:[h,z]=0}
       =intersection_(z_lambda !=0) ker(lambda|s_R).     (3.1)

Thus the possible s_z lie in a FINITE collection A of intersections of restricted-root kernels. Include s itself and omit the zero space when selecting walls. These real spaces are algebraic: relative to Q-bases, the roots and their restrictions have algebraic coefficients. Compactness makes their values on real Lie vectors purely imaginary, so multiplication by i gives real algebraic equations. Each such space is the Lie algebra of a connected closed subtorus of S_R: it is the identity component of an intersection of character kernels. This does not require the subtorus itself to be Q-defined.

For each nonzero a in A choose ONE nonzero vector h_a in a with real algebraic coordinates. This is possible by linear algebra over the real algebraic numbers. Let W_G be the full ordinary Weyl group of the compact pair (K,T(R)^0). It is finite, acts orthogonally on t_R and has algebraic matrix entries relative to t(Q).

Choose y in t(Q) which avoids BOTH the root hyperplanes and every hyperplane

    B(wy,h_a)=0,       a in A nonzero, w in W_G.         (3.2)

These are finitely many proper real linear hyperplanes. They are proper because B is positive definite on t_R, h_a!=0, and w is invertible. Rational points are dense in t_R, so a rational y avoiding their union exists. This is a choice in a Q-vector space, not an assertion that y is periodic.

## 4. All critical values lie on finitely many walls avoiding zero

For the Hamiltonian S-action on M, the annihilator of the differential image is

    (im dnu_z)^perp = {h in s_R:[h,z]=0}=s_z.           (4.1)

Indeed B(h,[a,z])=B([z,h],a). Therefore z is critical precisely when its isotropy Lie space a=s_z is nonzero.

Let A_a be the associated connected compact subtorus and L=C_K(A_a)^0. If z is fixed by A_a, then the regular centralizer torus of z contains A_a. This torus and T(R)^0 are maximal tori of L: the former centralizes A_a, and the latter contains it. Conjugacy of maximal compact tori in L therefore puts z into t_R. The usual Weyl conjugacy inside K gives

    z in Ad(L)(w y) for some w in W_G.                 (4.2)

Only this finite covering of the fixed set by such orbits is required; no claim about exact connected-component labels is needed.

L centralizes h_a, so on the orbit in (4.2),

    B(nu(z),h_a)=B(z,h_a)=B(wy,h_a).

Every critical value is consequently contained in one of the affine hyperplanes

    H_(a,w)={v in s_R:B(v,h_a)=B(wy,h_a)}.              (4.3)

Their equations have real algebraic coefficients. By (3.2) EVERY constant on the right is nonzero, so NONE of these entire supporting hyperplanes passes through0. This stronger genericity is useful: merely requiring pr_a(wy)!=0 would not ensure that a preselected scalar covering hyperplane avoids0.

The collection (4.3) is a finite affine hyperplane arrangement H. Its complement has finitely many full-dimensional open convex polyhedral chambers. All its points are regular values of nu whenever they lie in its image; empty fibers count as regular. Artificial cuts of the actual regular-value regions are harmless.

## 5. Zero existence without convexity or localization

Here is a complete elementary lemma, valid for ANY adjoint orbit of a compact connected semisimple group and any abelian subalgebra s.

Minimize F(z)=||pr_s z||^2 on the compact orbit, and at a minimum put eta=pr_s z. Differentiating along Ad(exp(ta))z gives

    F'(0)=2B(eta,[a,z])=2B([z,eta],a).

Thus [z,eta]=0. Suppose eta!=0. These commuting compact Lie elements lie in a common compact Cartan algebra t'. On a real root plane for a positive root alpha choose the orthogonal complex structure J so that

    [H,a]=alpha_R(H)J a,  J^2=-1, H in t'.

Whenever alpha_R(eta)!=0, this root plane is B-orthogonal to ker ad(eta), which contains s because eta belongs to the abelian s. Hence pr_s[a,z]=0 on this plane. Invariance then gives the SECOND derivative

    F''(0)=2||pr_s[a,z]||^2+2B(eta,[a,[a,z]])
          =-2 alpha_R(eta)alpha_R(z)||a||^2 >=0.

Therefore alpha_R(eta)alpha_R(z)<=0 for every root with alpha_R(eta)!=0; the product is zero for the other roots. The exact trace definition of the Killing form on the real root planes gives

    B(eta,z)=2 sum_(alpha positive)alpha_R(eta)alpha_R(z)<=0.

But the orthogonal-projection identity gives B(eta,z)=B(eta,eta)>0. This contradiction proves eta=0 and hence 0 belongs to nu(M).

For our selected y, Section4 shows that0 is not a critical value. It is therefore a NONEMPTY regular value, and a submersion chart at any point above0 proves that0 lies in the image interior. The same chart, integrating the strictly positive smooth Liouville density over a small fiber box, gives strictly positive density at0. Proper submersion over a neighborhood makes that density smooth there.

This proof uses no convexity theorem and no stationary-phase formula. As a CLASSICAL comparison, Kostant1973 Theorem8.2 proves the stronger equality pr_t Ad(K)y=conv(W_G y), after identifying its p and a with i g_R and i t_R by complexification. Its Weyl-average0 gives another zero-existence proof. That theorem is not a logical premise here.

## 6. Exact DH polynomial law and bounded cells

The torus S_R acts Hamiltonianly on M, with moment map nu after the B identification. Its infinitesimal action is effective: if h in s acts trivially, then at y it acts trivially on every root plane. Regularity of y makes every root plane tangent, so alpha(h)=0 for every root, forcing h=0 by semisimplicity. The group kernel is therefore finite. Divide it out to obtain an effective compact torus of the SAME dimension r. No individual rational Lie vector is alleged to generate a closed circle.

M is compact, so nu is proper. The actual Section3 Duistermaat--Heckman theorem applies: the normalized pushforward is globally absolutely continuous and its density is polynomial of degree at most D-r on each connected component of regular values. Regular reductions may have finite stabilizers and be V-manifolds; globally free action is not required. Coordinate changes on s_R and probability normalization do not change polynomial degree. Here D is the number of positive roots for the regular orbit; D>=rankG>=r, so the bound is nonnegative.

Every chamber of H is a region of regular values, so rho is one polynomial on that chamber. It is zero outside the bounded compact image. If a chamber is unbounded, it contains an exterior open set where the measure is zero; a polynomial zero on that open set is identically zero throughout that chamber. Thus every unbounded chamber has zero density. The chamber containing0 has rho(0)>0 and is consequently bounded.

Let Omega be the union of ALL bounded chambers, including any on which the density happens to be zero. There are finitely many, and their union is bounded and open. The pushforward measure is exactly the sum over these bounded cells of rho_cell(v)dv; the hyperplanes themselves have measure zero by absolute continuity. Every cell is a bounded full-dimensional open convex polyhedron, every boundary lies in H, and every full supporting hyperplane of this chosen arrangement avoids0. On the central cell the constant coefficient is rho(0)>0. The coefficients on the other cells may be arbitrary real numbers; their actual polynomial densities are nonnegative on those cells.

There is no need to identify the moment polytope's facets, to require density continuity across walls, or to assume a global polynomial density. This is the exact cellwise input for the separate analytic criterion.

## 7. Normalizer components, rational group action and wall polynomial

The normalizer N_G(S) is Q-defined and smooth in characteristic zero. Its Lie algebra consists of a with [a,s] contained in s. The S-weight decomposition shows that every nonzero-weight component of such an a is zero; conversely the zero-weight subalgebra is Lie(C_G(S)). Thus

    LieN=LieC,   N^0=C,   W=N/C finite etale overQ.      (7.1)

The action on s is faithful: an automorphism of a torus with identity differential is the identity in characteristic zero, because its integral character-lattice matrix is determined by that differential. Its kernel in N is exactly C. Therefore W is an actual finite Q-group scheme acting Q-linearly on s. Its finite set of geometric matrices is Galois-STABLE. It need not be constant overQ, and no rational normalizer representative is used.

The free right W-action on Y=G/C has an actual affine quotient Z=G/N. The map Y->Z is finite etale and surjective; smoothness and geometric integrality descend. The coordinate equivariance from Section2 holds for the FULL group scheme action, so any W-invariant polynomial in mu descends coefficientwise.

Normalize each hyperplane equation (4.3) as

    ell_(a,w)(v)=1-B(v,h_a)/B(wy,h_a),

and let p0 be their product. Choose a finite Galois extension E/Q containing all coefficients and all geometric W-matrix entries. Form

    q(v)=product_(sigma in Gal(E/Q))sigma(p0)(v),
    f(v)=product_(u in W(Qbar)) q(u v).                 (7.2)

Then q has rational coefficients, f(0)=1, and f is W-invariant. The Galois stability of the finite W-set gives sigma(f)=f, so f also has rational coefficients. Every original hyperplane is a factor of f, so f vanishes on the boundary of EVERY bounded cell. Extra conjugate or group-translated zeros cause no problem. This construction would not be justified for an arbitrary finite algebraic matrix group whose Galois conjugates generate an infinite group; the actual Q-group-scheme action in (7.1) is the necessary finiteness gate.

The rational positive quadratic v->B(v,v) is invariant under the same full W action. It is available for the separate analytic amplification step, without any need for a rational/real representative of a component or for a density symmetry assumption.

## 8. Apply the analytic theorem on the actual quotient

Take the actual algebra A=O(Y), the rational coordinates mu and the normalized native trace from Section2. Section6 gives a finite bounded polyhedral-cell law of degree at most b=D-r, with0 inside a cell of positive constant density. Section7 gives its Galois-stable finite W-action and the positive invariant rational quadratic B(v,v), and the required full supporting hyperplanes all avoid0. These discharge EVERY functional and symmetry hypothesis of [the polyhedral analytic criterion](POLYHEDRAL_PIECEWISE_DENSITY_NORMALIZER_CRITERION.md).

That theorem therefore selects a rational W-invariant boundary polynomial, its rational quadratic amplifier and an integer L>=r+b=D large enough that, on the ACTUAL Y product,

    p=T*t^-1*f(Tmu)^2, F=p^L(1+p), H_marker=t^r

have zero pure moments and strictly positive rational marked moments at EVERY order. It also proves eventual vanishing for each fixed integer marker exponent ell<r. This is an unrestricted tail-and-binomial proof, not a sampled moment inference. No density symmetry under W is required.

Since f(TAd(n)^-1 mu)=f(Tmu) for every geometric n inN, the ACTUAL phase coefficients and scalar Laurent unit descend to Z=G/N. Its normalized trace restricts to the Y-trace, so for EVERY m>=1 and integer ell,

    Lambda_Z(t^ell F^m)=Lambda_Y(t^ell pi^*(F)^m).

Thus the SAME unit t^r survives every order downstairs, while every integer ell<r eventually vanishes. This proves the least-unsafe-exponent assertion. No individual coordinate or arbitrary upstairs phase is presumed to descend. Rational positivity before scalar extension gives the identical nonzero rational values over EVERY characteristic-zero field.

### Direct downstairs prime divisor and support

Section5 gives g inK with pr_s Ad(g^-1)y=0. The corresponding ACTUAL Q-defined mu0 fiber on Y is real-nonempty, hence has a Qbar-point y0. Its image z0 onZ is an actual geometric point. Even though mu itself may not descend, the descended polynomial family has f(Tmu)(z0)=f(0)=1 for allT. Consequently

    p(z0,t)=T/t,   (1+p)(z0,t)=2+t^-1.

At t=-1/2 the latter has derivative-4 and p=-1 is a unit. Smooth geometric integrality of Z times Gm gives an ACTUAL prime divisor through that point on which ord(1+p)=1 and ord(p)=0. Therefore F has valuation1 and is absolutely not a unit times any proper power. This is a direct downstairs certificate, valid over every algebraically closed extension.

At the same actual point,

    F(z0,t)=(1+t^-1)^L(2+t^-1),

whose coefficients at exponents0,-1 are2,2L+1. They are nonzero in every characteristic-zero field, proving actual primitive Laurent support. These geometric conclusions come from an ACTUAL zero image, not from an abstract probability law.

### Individual factors

Every nonzero real regular function onZ has strictly positive square trace: pull it back injectively toG and integrate its nonzero square on the dense compact group orbit. Rational Gram/LDL with positive rational pivots written as sums of rational squares makes the square trace anisotropic after EVERY formally real extension, with no embedding of that whole field intoR. Only the zero phase can be pure, so the compact quotient factor is Mathieu. The scalar split-torus factor is Mathieu over every characteristic-zero field by the scalar Laurent theorem. The displayed pair proves their product fails over every formally real field; the pair itself holds over ALL characteristic-zero fields.

## 9. Homogeneous sources and full centralizer normalizers

**Corollary 2.** Let H be ANY closed Q-subgroup of G normalizing a positive-dimensional Q-torus S. The ACTUAL affine homogeneous source G/H has the same conclusions after pullback of the displayed phase: one scalar Laurent unit survives EVERY order over EVERY characteristic-zero extension, the phase is absolute unit-nonpower with primitive support, and its least unsafe integer monomial exponent is r=dimS. Over every formally real field its compact factor and split torus are Mathieu, while their product is not.

Indeed H is contained in N_G(S). Compactness excludes a positive-dimensional unipotent radical in H^0, so H^0 is reductive and G/H is affine. The actual map G/H->G/N_G(S) is smooth and surjective. [The homogeneous-quotient trace theorem](HOMOGENEOUS_QUOTIENT_TRACE_AND_COUNTEREXAMPLE_PULLBACK.md) pulls back this GIVEN pair with its scalar unit and exact marked values. Injectivity preserves Laurent support, and the ACTUAL order-one divisor proved above supplies the strong nonpower pullback certificate. Compact Haar positivity and rational Gram prove the factor conclusion separately; no general global-Mathieu-ascent theorem is used.

If H^0 has positive-dimensional center, S=Z(H^0)^0 is a nontrivial Q-torus normalized by the FULL H. Thus Corollary2 includes this case even when the full stabilizer permutes all rational lines of that center.

**Corollary 3 (all full centralizer normalizers).** For ANY0!=x in g(Q), let C=C_G(x). Then S=Z(C)^0 is a positive-dimensional Q-torus,

    C_G(S)=C,       N_G(C)=N_G(S).

Therefore Theorem1 applies to the FULL centralizer-normalizer quotient G/N_G(C), regardless of the dimension of its central Lie algebra. It also applies by Corollary2 to every closed Q-subgroup normalizingC, including every full subgroup preserving the rational lineQx.

To prove this, x is semisimple by compactness and the FULL C is connected reductive by the established semisimple-Lie-centralizer argument. Its central Lie algebra contains the nonzero x. In characteristic zero this central Lie algebra is Lie(Z(C)^0), so S is a positive-dimensional Q-torus and x belongs toLieS. Centrality gives C contained inC_G(S). Conversely centralizingS fixesx by differentiation, so C_G(S) is contained inC_G(x)=C. Every normalizer ofC preserves its characteristic connected centerS, and every normalizer ofS preservesC_G(S)=C; hence the full normalizers are equal, scheme-wise. If H preservesQx, invariance of the nonzero Killing norm gives Ad(h)x=+/-x and thus H normalizesC.

These are NEW existence pairs selected from the joint S-moment law, with marker t^(dimS). They do not identify the phase or its marked values with the older scalar-height construction whose marker was t^-D, nor assert that an arbitrary old witness descends to a larger normalizer. A stabilizer with no normalized positive-dimensional rational torus is outside this corollary.

## 10. Precise sources and boundaries

The actual torus-centralizer and maximal-Q-torus scopes are Milne17.37--17.40,17.59 and17.82--17.83, already retained with canonical PDF hash503f1a7d511632af05541c52fb2bf6fc36f5a75b0323f2f9110f03dcce11fd67. The compact coadjoint affine/trace/KKS ingredients are the exact already established arbitrary compact-Q orbit arguments. Basic compact root-plane decomposition and maximal-torus conjugacy are used explicitly in Sections3--5; no root-system classification is imported.

For the density law the necessary primary is Duistermaat--Heckman1982, original retained PDF hash e314c2033972ab39d244b8bc98ddfb37d8ade96c8c36bfa95e05d82fc1033be5: full Section3 printed264--266, especially Lemma3.1, Proposition3.2 and Corollary3.3, with the variation theorem and proof printed260--264. It is used for the ACTUAL compact torus S, dimensionr, not for an unproved periodic rational vector. No original/corrected Section4 stationary-phase formula, torus localization or global C^q premise is used.

The classical convexity comparison is B.Kostant, On convexity, the Weyl group and the Iwasawa decomposition, Ann.Sci.ENS6(1973),413--455, DOI10.24033/asens.1254. [Primary](https://www.numdam.org/item/10.24033/asens.1254.pdf). Actually inspected: the Section3.1 convex-hull definition and the ENTIRE Theorem8.2 statement and proof printed452--454, both text and rendered pages. Earlier normal-form/reflection-closure results used in that published proof remain imports; no whole43-page audit is claimed. Its compact application uses the complexification with p=i k and a=i t, so it is not obtained by applying a noncompact p-projection theorem to a compact group with p=0. The elementary zero lemma here makes this reference comparison only.

The wall arrangement, stronger genericity, nonempty regular zero, positive central density, bounded-cell reduction, full Q-group-scheme action and rational invariant wall polynomial are proved here. The separate complete analytic theorem is then applied in Section8. Affine coadjoint model/trace, rational Gram and the scalar Laurent/homogeneous-quotient facts are used only in their established generic scopes. The result concerns these ACTUAL G/Q,S/H data and their scalar extensions. It is not a classification of all homogeneous varieties, an arbitrary-phase descent statement or a claim of priority for moment-map and normalizer mechanisms.
