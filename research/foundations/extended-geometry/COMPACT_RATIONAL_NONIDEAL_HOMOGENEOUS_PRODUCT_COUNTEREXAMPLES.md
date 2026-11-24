# Compact rational homogeneous products with nonideal stabilizer Lie algebra

A principal Weyl-face cross-section gives an actual rational normalizer circuit whenever the stabilizer Lie algebra is not an ideal of the acting semisimple Lie algebra. The construction uses a torus moment law on the cross-section and its linear-root Liouville weight. It does not require a regular-H-image or a full-dimensional radial density.

## 1. The actual rational theorem

Let G/Q be connected semisimple with compact G(R), and let H be a closed Q-subgroup with reductive identity component. Suppose h=Lie(H^0) is NOT an ideal of g=LieG. Throughout put B=-Killing_g and p_orth=h^perp. Then there is an actual rational phase F on (G/H) times split G_m and ONE fixed scalar Laurent unit t^e such that, for every m>=1,

    Lambda(F^m)=0,          Lambda(t^e F^m)>0.                (1)

Here lambda is the native normalized rational trace and Lambda is lambda composed with Laurent constant term. Complete marked values are positive rational numbers. Thus the SAME phase and fixed unit work after every characteristic-zero scalar extension of this Q-model. The phase has an actual geometric prime divisor of valuation one and primitive actual Laurent support. Over every formally real extension the two separate factors are Mathieu, while their product fails.

The construction starts by choosing x in h(Q)^perp, x nonzero, whose projection mu:M->h_R on M=Ad(G(R)^0)x is not identically zero. Section10 proves that existence of THIS INPUT is equivalent to h not being an ideal, and supplies a finite rational basis-bracket certificate. If a particular rational x with these properties is already given, the same construction applies to that orbit. Choose the KKS sign so that i_(a_M)omega=-d B(mu,a); its normalized positive Liouville probability is the same invariant Haar probability for either sign. Put dimM=2D. The fixed marker satisfies

    1<=e<=D.                                                (2)

No least unsafe exponent or smallest normalizer power is claimed.

The analytic input is [the positive-stratified invariant-test theorem](POSITIVE_POLYHEDRAL_STRATIFIED_INVARIANT_TEST_COUNTEREXAMPLES.md). Only its linear-stratum Euler/leading-term/damping construction is needed. Its separate normal-delta decomposition and CDKW derivative-principle corollary are NOT premises of the present application.

Actual geometry is separate from the moment construction. Compactness makes x semisimple; the generic compact-Q coadjoint lemmas give its closed smooth geometrically integral affine orbit and full connected reductive centralizer. Matsushima's affine direction for H^0, finite components and affine descent give the actual smooth geometrically integral affine quotient G/H. Complete reducibility gives its unique Q-valued normalized invariant trace, compatible with every scalar extension. Over R it equals Haar on the Zariski-dense compact identity component and restricts from O(G); on M its inverse-Haar pushforward is normalized KKS probability. These generic facts are recorded in [UCQ](UNIVERSAL_COMPACT_RATIONAL_COADJOINT_PRODUCT_COUNTEREXAMPLES.md), Sections2--3, [RHT](REDUCTIVE_HOMOGENEOUS_AUTOMORPHISM_INVARIANT_TRACE.md), Sections2/14, and [HQP](HOMOGENEOUS_QUOTIENT_TRACE_AND_COUNTEREXAMPLE_PULLBACK.md), Sections2--3. No earlier counterexample phase or all-automorphism trace theorem is a premise.

## 2. Principal face and the local symplectic cross-section

Take an ACTUAL maximal Q-torus T_H of H^0 from the outset, and let T=T_H(R)^0 be the corresponding maximal compact torus of K=H(R)^0. Choose a positive Weyl chamber C in its real Lie algebra dual, identified by B=-Killing_g. The principal-cross-section theorem [LMTW, Theorem3.1] supplies a unique relatively open Weyl face sigma such that the radial image is contained in closure(sigma), together with a connected principal cross-section

    Y=mu^-1(sigma),       L=K_sigma,       S=Z(L)^0.

The image mu(Y) is dense in the full radial image. The induced space K times_L Y identifies with an open dense subset of M [LMTW, Remark3.9(b)]. Its FULL-MEASURE property is proved separately in Section3; density alone is not substituted for that property.

The local symplectic assertion and its volume data have a direct check. The natural slice U_sigma is open in l^* for l=LieL. Principal-image containment gives mu^-1(U_sigma)=Y. At eta in sigma, ad(eta) is invertible on l^perp. Fundamental directions from l^perp map under dmu onto that complement, so mu is transverse to the slice. Hence Y is an actual real-analytic submanifold. Its tangent space is the preimage of l^* under dmu and is symplectically orthogonal to the l^perp fundamental directions. On that complement the KKS form is nondegenerate because ad(eta) is invertible. Thus the tangent space splits symplectically and Y is symplectic.

The restricted moment map is L-equivariant and takes values in sigma, contained in z(l)^*. All Hamiltonian components for [l,l] are identically zero on Y, so their vector fields vanish. The connected derived group [L,L] acts trivially. The remaining action is that of the compact centre torus S, up to finite isogeny. The source theorem supplies connectedness of Y and existence/uniqueness of the principal face; the displayed local argument writes the required transversality/KKS assertion directly.

## 3. Why the principal open part has full measure

Let c=dim l. Every element of compact h is semisimple. The invariant polynomial

    P_c(eta)=[w^c]det(wI-ad_h(eta))

is nonzero on sigma, where the centralizer has dimension c, and zero on each strict boundary face of closure(sigma), where centralizer dimension is larger. Since Y is nonempty, P_c(mu) is a nonzero analytic function on connected M. Its proper analytic zero set contains M minus K.Y. Thus this complement has Liouville measure zero.

This argument permits a positive centre of h. If sigma is the central face and the image never leaves it, K.Y=M directly; the leading characteristic coefficient is1. No conclusion from density alone, no regular-H-orbit support hypothesis and no full-dimensional radial absolute continuity are used.

## 4. The effective torus and the actual linear moment subspace

The S-action on connected real-analytic Y has a global connected kernel with Lie algebra k. Each component B(mu,xi) for xi in k is constant on Y, so its image lies in an affine flat

    A=c0+k^perp inside s^*,        s=Lie(S).

This is its actual affine hull. To see this without an unstated generic-freeness assertion, restrict the linear adjoint representation of S on g. Stabilizer Lie algebras on M, hence on Y, are intersections of kernels of subsets of finitely many weights. There are only finitely many such spaces. For every stabilizer space a strictly larger than k, choose xi_a in a minus k. Its vector field on connected analytic Y is not identically zero, so its zero set is proper analytic and null. Outside the finite union of these sets, the stabilizer is exactly k. The derivative of the torus moment map consequently has image k^perp there. Its image therefore contains a relatively open subset of A.

The principal radial image is dense in the full radial image by the principal-face theorem and continuity. Since mu(x)=0, zero belongs to closure(mu(Y)). Hence c0=0 in the quotient by k^perp, and A is the LINEAR flat

    A=k^perp.                                                (3)

The connected kernel is a subtorus, so its Lie algebra is rational relative to the torus character lattice. All groups and root coordinates in this real construction can be taken over the real algebraic numbers; thus A is an algebraic linear flat. It need not be a Q-defined subspace of the chosen rational Cartan, and need not be preserved by the full geometric W.

Let n=dim A. Nonzero moment projection implies n>=1: if n=0 then the connected cross-section image is the single point0, so its dense radial closure would make all of mu zero. The effective compact torus S/Ker(S acting onY)^0 has dimension n; remove any remaining finite global kernel without changing its Lie moment coordinates. Its moment map is mu_Y into the convex relatively open set sigma intersect A. It is proper INTO THAT OPEN SET: compact target subsets lie inside sigma, and their preimages are closed compact subsets of the original compact M. No global properness into all of s^* is asserted.

## 5. Finite algebraic critical hyperplanes on the noncompact cross-section

Choose a maximal torus T_G of G containing S. Its Cartan coordinates and roots are algebraic. The finite S-weight decomposition of g gives finitely many possible stabilizer Lie spaces a on Y. A critical point of the effective torus moment map has a strictly larger than k. Choose, once for each such a, a nonzero real algebraic xi_a in a minus k. Such a choice is possible because the subtori and k have algebraic Lie algebras.

For y fixed by the corresponding subtorus A_a, place y in t_G by conjugation inside C_G(A_a): extend A_a to a maximal torus in the connected reductive centralizer of the semisimple y, then use compact torus conjugacy. The conjugation fixes xi_a. The point in t_G belongs to the finite Weyl orbit of x. Its coordinates are real algebraic, since x is rational semisimple. Thus every critical value lies in one of finitely many hyperplanes

    B(eta,xi_a)=B(w x_t,xi_a).                               (4)

The restriction of this equation to A is a GENUINE hyperplane, because xi_a is not in k. Add the homogeneous root hyperplanes bounding sigma and finitely many algebraic equations defining A. No uncounted family of critical walls is introduced by noncompactness of Y.

## 6. Local proper-torus DH and global finiteness

Put q=(dim K-dim L)/2, so dim Y=2(D-q). Here is the complete local proper-torus degree argument, including finite stabilizers. Let an effective compact n-torus act Hamiltonianly on a symplectic Z of dimension2d, with moment map psi and convention i_(xi_Z)omega=-d<psi,xi>. Take a small ball U of regular values over which psi is proper. An empty tube gives zero density; otherwise shrink around a nonempty compact fibre. A proper submersion has locally constant image, so every nearby fibre occurs.

The annihilator of dpsi is the stabilizer Lie algebra. Thus the tube is locally free, even though finite stabilizers may remain. An invariant metric constructs an invariant Lie(T)-valued connection alpha with alpha(xi_Z)=xi by projection onto vertical tangent directions and inversion of the injective infinitesimal action. Average horizontal lifts to be T-equivariant and integrate over compact target paths. Properness ensures existence of these flows. The resulting Ehresmann trivialization is T-equivariant, and the quotient tube is an orbifold product Q times U, where Q=psi^-1(eta0)/T is compact.

Put theta_conn=<psi,alpha>. Invariance gives i_xi dtheta_conn=-d<psi,xi>; hence beta=omega-dtheta_conn is closed, invariant and horizontal. The torus curvature F_conn=dalpha is closed and basic as well. On a moment level eta,

    omega_eta=beta_eta+<eta,F_eta>,
    [omega_eta]=c0+sum_i eta_i e_i.                          (5)

In the first line F_eta denotes the restriction of F_conn to that quotient fibre, and the pairing is between eta and this Lie(T)-valued curvature. The classes in the second line are constant under the product identification. Thus orbifold Stokes and the finite expansion of omega_eta^(d-n)/(d-n)! give a polynomial reduced volume of degree at most d-n. Orbifold Stokes follows in finite quotient charts with invariant partitions of unity and reciprocal finite-group weights.

The upstairs Liouville density is a fixed positive constant times this reduced volume. Choose n vertical action vectors, n transverse vectors with dpsi values the coordinate basis, and a reduced tangent complement. The vertical/transverse symplectic block has identity pairing up to orientation. Its top exterior power factors into reduced Liouville volume, moment coordinate differentials and the n connection forms. Fibre integration gives the fixed torus Haar volume. Finite quotient chart weights are exactly the reciprocal stabilizer factors in orbifold integration; removing a finite global kernel changes only the constant.

Apply this argument to Z=Y, for which d=D-q and properness is into sigma intersect A as proved in Section4. Its UNNORMALIZED local Liouville pushforward therefore has degree at most

    D-q-n.                                                  (6)

No finiteness of the entire Liouville measure of Y, compactness of Y or globally proper moment map into all of s^* is assumed. This is the connection-curvature variation mechanism of DH Section3, written over compact regular-value tubes. The finite stabilizer/orbifold version is explicitly included; no free-action theorem is substituted.

The complement of the critical locus is full measure on analytic Y by Section4. Thus no additional critical-value singular measure occurs for the Y torus law. On every connected chamber in A intersect sigma cut out by(4), the local polynomials agree on overlaps and hence form ONE polynomial of degree(6). Each chamber is convex, so connected. There are finitely many chambers because the hyperplane arrangement is finite. This supplies a finite global polynomial description without assuming Y compact or assuming a globally bounded torus moment map alone proves finiteness.

## 7. The precise Liouville Jacobian on the induced open part

At y in Y, the tangent space of the open part K times_L Y splits into T_yY and the fundamental orbit directions from l^perp in h. These are symplectically orthogonal: for a in l^perp and v in T_yY,

    omega(a_M,v)=-B(dmu(v),a)=0,

because dmu(v) lies in s^* inside l^*. The orbit block is the KKS form of the coadjoint orbit through eta=mu(y):

    omega(a_M,b_M)=-B(eta,[a,b]) under the stated moment convention.

With the corresponding symplectic orientation on K/L, the positive volume of this block is a fixed positive constant times the product of the positive root linear forms that DO NOT vanish identically on sigma. There are q such roots. Therefore for every invariant test P the actual Liouville integral on M is, up to a single positive normalization constant,

    integral_Y P(mu(y)) V_sigma(mu(y)) dvol_Y,               (7)
    V_sigma(eta)=C product_(alpha>0,alpha|sigma not0)alpha(eta).

Here C>0 is independent of eta, V_sigma is polynomial of degree q and is strictly positive on sigma. This is a LINEAR root product, not a squared Euclidean orbit Jacobian. No integration over individual non-invariant Cartan coordinates downstairs is inferred.

This is the standard symplectic-induction volume mechanism. The argument above writes the needed orthogonal-block/Pfaffian statement locally, so no globally proper induction theorem is applied outside its properness domain.

## 8. The required finite positive polyhedral radial law

Multiply each polynomial torus density from Section6 by the degree-q polynomial V_sigma. The actual radial probability of mu on the whole closed Weyl chamber equals this weighted law, by full measure of K.Y and(7). Thus it has nonnegative polynomial densities on finitely many relatively open polyhedral cells of the SAME linear flat A, with

    n+degree(rho_cell)<=D.                                  (8)

Only the positive radial probability is normalized; no separate finiteness assumption on unweighted Y measure is needed. Compactness of M makes the radial support compact. A nonzero polynomial density on an unbounded chamber would vanish on a nonempty open tail outside that compact support and hence vanish identically. Delete it. All retained cells are bounded. Their walls and flat A are real algebraic. Boundary/root-wall mass not represented by the principal law is zero by Section3, rather than discarded via the CDKW regular-support hypothesis.

Zero belongs to the radial support since mu(x)=0 and source probability has full support. There is no atom at zero: mu is not identically zero, so some analytic coordinate has a proper null zero set containing mu^-1(0). Finiteness of the cell decomposition now gives a positive-dimensional nonzero cell whose closure contains0. These are precisely the positive-stratified construction's hypotheses. No open-regular-H-chamber requirement remains.

## 9. Actual rational invariant tests, phases and fields

Use the rational Lie coordinates of the actual maximal Q-torus T_H chosen in Section2 for the ambient radial space. The finite full-component Weyl matrix image W_full of N_H(T_H)/T_H is Galois-stable and preserves B. It need not preserve A or an individual cell. The conditional positive-stratified theorem requires only W-invariant POLYNOMIAL TESTS, not those extra symmetries.

The connected Chevalley restriction theorem with its reductive extension, degreewise Q descent and full-component adjustment gives

    Q[h]^H = Q[t_H]^(W_full).

For each invariant P let Ptilde be its unique H-invariant lift. The actual homomorphism

    theta(P)(gH)=Ptilde(pr_h Ad(g^-1)x)

descends regularly on the actual represented affine quotient. Native inverse Haar has exactly the weighted radial law(7). At the actual rational identity image theta(P)=P(0).

Apply the positive-stratified theorem using(8). Its ambient rational W/Galois boundary polynomial f has f(0)=1 even if A is not Q-defined or W-stable. For a sufficiently large integer L0,

    p=(1+t)t^-1 theta(f((1+t)u))^2,    F=p^L0(1+p)

has all unmarked moments zero and one fixed unit t^e surviving every order, with1<=e<=D. Rationality is native; individual cell coefficients need not be rational. At identity p=1+t^-1, and the modifier has simple zero -1/2, derivative -4, with p=-1. This gives direct geometric order1 and adjacent support coefficients2,2L0+1, including full-component quotient cases.

## 10. An exact nonideal-Lie criterion for existence of input

If h is an ideal of g, then p_orth is its complementary B-orthogonal ideal. Every x in p_orth has its entire G-orbit in p_orth, so mu is identically zero. Hence THIS projection method has no input in that case.

Conversely, if h is not an ideal, then [h,p_orth] is nonzero. Indeed reductivity and B-invariance give [h,p_orth] subset p_orth; if it were zero, then [g,h] subset h and h would be an ideal. Choose actual rational x in p_orth and a in h with [a,x] nonzero. Put b=[x,a]. At x the tangent vector [b,x] has projected pairing B(a,[b,x])=B([x,a],b)=B(b,b)>0 with a. Thus the moment projection is nonconstant, providing the input used above. The input can be found by a finite rational basis bracket calculation.

Consequently existence of an x satisfying the projection method's input is equivalent to h NOT being an ideal. This is an equivalence about THIS INPUT, not a converse to Mathieu product failure.

For an absolutely simple compact G and any positive-dimensional proper H^0, h is a nonzero proper subalgebra and cannot be an ideal. Thus all such actual rational quotients, including full disconnected H with that identity Lie algebra, satisfy the all-order conclusion by this route. Positive-dimensional centres of H cause no additional principal-face obstruction, and their reductive invariant restriction and trace have the same full-component and field quantifiers proved in Sections9 and11.

The ideal case, finite H, H=G and any attempted classification of all compact homogeneous quotients are outside this input theorem. Some have other existing group/coadjoint or sphere constructions. No failure of this input proves failure of every construction.

## 11. Fields, geometric certificates and separately Mathieu factors

Every full marked value belongs to Q because the actual phase and normalized trace are rational. The analytic proof makes it positive over R; injectivity of field extension preserves its nonzero rational value over EVERY characteristic-zero extension. This is propagation of one given rational phase on one specified Q-model, not inference of arbitrary twisted forms from a real label.

At identity the actual invariant-test evaluation is P(0), so the simple modifier certificate and adjacent Laurent coefficients in Section9 are on G/H itself. Over Qbar, smoothness makes the derivative -4 a transverse local parameter and supplies an actual valuation-one prime divisor. Its valuation excludes every expression unit*Q0^a with a>=2. The same geometric point and parameter persist after algebraically closed extension. Coefficients2 and2L0+1 at adjacent exponents prove literal Laurent-support primitivity.

For the compact factor, a nonzero real regular function pulls back nontrivially to the Zariski-dense compact real identity component of G, so its square has strictly positive Haar trace. Finite rational Gram matrices are positive definite. Rational LDL with positive rational pivots, each a finite sum of rational squares, makes the Gram form anisotropic over every formally real field. Every individual element uses a finite rational span, so only zero is pure and its trace kernel is Mathieu. The split scalar Laurent factor is independently Mathieu by the scalar DvK theorem and its coefficient-field argument. These generic facts are recorded in [SCCP](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md), Sections3 and13, and [HQP](HOMOGENEOUS_QUOTIENT_TRACE_AND_COUNTEREXAMPLE_PULLBACK.md), Section7. No global Mathieu ascent or algebra-valued Laurent theorem is used.

Thus over every formally real extension the two separate factors are Mathieu and their product fails. The pair itself works over all characteristic-zero extensions, where no compactness or positivity assertion about that field's rational points is required.

## 12. Exact sources, antecedents and method boundary

[LMTW] E. Lerman, E. Meinrenken, S. Tolman and C. Woodward, Non-abelian convexity by symplectic cuts, Topology37(1998),245-259, DOI10.1016/S0040-9383(97)00030-X. The exact retained [author PDF](sources/lerman-meinrenken-tolman-woodward-nonabelian-cuts-author.pdf), from https://sites.math.rutgers.edu/~ctw/lmtw.pdf, has SHA256 333a4c9b0c4fee0473736f137b6a277a77e32fa4412c1f842b6bff1a1a6ee849. The imported principal-face/connected-cross-section statement is author Theorem3.1, with complete proof and supporting Theorem3.8, Lemmas3.10--3.11 and Remark3.9(b), printed4--8. Its numbering differs from earlier versions. The principal-orbit/slice facts used there remain classical inputs; Proposition3.2 cites the orbifold slice adaptation rather than proving all prior slice theory. No whole-paper certification is claimed. The full-measure characteristic-coefficient argument, effective-kernel affine-hull proof, critical-wall finiteness, proper local DH with finite stabilizers and root-volume Pfaffian are written here.

The ordinary DH mechanism uses the retained original variation/Section3 source, printed260--266, SHA256 e314c2033972ab39d244b8bc98ddfb37d8ade96c8c36bfa95e05d82fc1033be5, with the precise local/proper finite-stabilizer extension proved in Section6. No original DH Section4 formula is used. LMTW's symplectic cuts and proper-open-target torus discussion give a separate classical route to compact target localization, but are not substituted for the displayed local proof.

Rational invariant lifting uses the retained Etingof Lecture10 Theorem10.1 ENTIRE proof and Remark10.2(2), printed54--55, SHA256 629d82eeb9649aee60c6ad23558292c2465cc5e52790d8dec7163cfcb065c417. The graded restriction map descends degreewise from C because the finite-dimensional invariance kernels and restriction maps are Q-defined. Full components adjust to the normalizer by maximal-torus conjugacy; no rational representatives or rational highest-weight models are needed. The exact same-field normalized trace/Matsushima/quotient and compact Gram/DvK scopes are the generic ones cited in Sections1 and11. No seven-family model theorem or earlier constructed bad phase is a premise.

The analytic [positive-stratified invariant-test theorem](POSITIVE_POLYHEDRAL_STRATIFIED_INVARIANT_TEST_COUNTEREXAMPLES.md) proves the complete coefficient, damping and origin-certificate steps. Only its linear-stratum construction is used. Its positive-distribution lemma and optional CDKW derivative-principle corollary are independent results; neither CDKW nor regular-H-orbit support is a premise of the present principal-face route.

Long's interval-Laurent normalizer and related compact-product circuits are antecedents. Complex all-order pairs also admit an elementary unitization: if ell(f^m)=0 and ell(h f^m)!=0 for every m, then P=f(1+t h) has Lambda(P^m)=0 and Lambda(t^-1 P^m)=m ell(h f^m). On a smooth geometric domain its modifier has a transverse zero on D(fh), and coefficients f,fh give primitive support. Liu's deposited preprint v2.2, Theorem3.1/Corollary3.2, supplies such complex compact-group pairs within its stated coherent-state hypotheses, https://zenodo.org/records/21689760. Thus complex Laurent-unit/geometric strength alone is not a new complex-existence assertion. Native Q-form construction, the same rational values after all extensions, and genuine formally-real two-factor instability are separate statements proved here. No priority claim is made.

The exact equivalence in Section10 is about existence of INPUT for THIS moment-projection construction, not about all possible counterexamples. Finite H and ideal h can fail that input; no general conclusion about their products is inferred. The proof gives no classification of all compact homogeneous spaces and no estimate from a finite moment prefix.
