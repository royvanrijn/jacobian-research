# Positive polyhedral strata and invariant-test counterexamples

This conditional theorem gives one fixed Laurent unit marker at every positive order from a finite positive polyhedral stratification and an actual invariant-test algebra. A separate elementary lemma decomposes positive polyhedral derivative distributions into polynomial measures on faces. Group applications must supply their own moment law and actual algebraic homomorphism.

## 1. Stratified measure and actual algebraic input

Let r>=1. Let Omega_s be finitely many pairwise disjoint bounded relatively open convex polyhedral cells in R^r, with affine hull L_s defined over the real algebraic numbers. Put n_s=dim L_s. The reference measure on each cell is relative Euclidean Lebesgue measure; in dimension zero it is counting measure at the point. Let rho_s be a nonzero nonnegative real polynomial on L_s of degree at most b_s, and assume

    nu=sum_s rho_s(u)du_(L_s)|Omega_s

is a probability measure. Overlapping finite descriptions can first be refined into disjoint polyhedral cells and their density polynomials added. Delete zero-density cells.

Assume nu has NO atom at zero and at least one positive-dimensional cell has zero in its closure. These are automatic from zero in the support together with no atom at zero for a finite such decomposition: some nonzero cell closure must contain zero, and a zero-dimensional cell closure can do so only when that cell is the forbidden atom.

Let W be a finite Galois-stable geometric linear matrix group over Qbar preserving a positive definite rational quadratic form q. Density and cells need not be W-invariant. Put R=Q[u_1,...,u_r]^W. Suppose an ACTUAL unital Q-algebra homomorphism theta:R->A and normalized Q-linear lambda:A->Q satisfy

    lambda(theta(P))=integral P(u)dnu(u) for EVERY P in R.      (1)

Assume compatible scalar extension of the actual algebra and functional. Individual coordinates u_i need not descend to A.

Call a cell Origin if zero lies in its closure. Its affine hull necessarily contains zero and is a linear subspace. On every linear L_s decompose

    rho_s=sum_j rho_sj

into polynomials homogeneous about zero within L_s. Define

    d_star=min_Origin(n_s+j_s),                              (2)

where j_s is the first nonzero homogeneous degree of rho_s. Since no Origin cell is zero-dimensional, d_star>=1.

There exist a rational W-invariant f with f(0)=1 and an integer L>=1 such that, with T=1+t,

    a(T)=theta(f(Tu)),       p=T*t^-1*a(T)^2,
    F=p^L*(1+p),             H=t^d_star,                     (3)

the native Laurent functional Lambda=lambda composed with constant term satisfies

    Lambda(F^m)=0,           Lambda(H F^m)>0 for EVERY m>=1.   (4)

The complete marked values are positive rational numbers, so the SAME actual Q-phase and fixed Laurent unit have the corresponding zero and nonzero values after every compatible characteristic-zero extension.

If A=O(X) for a smooth geometrically integral affine Q-variety and there is an ACTUAL geometric point xi with theta(P)(xi)=P(0) for every invariant P, then F has an actual geometric prime divisor of order one and primitive Laurent support. Hence it is not a unit times any proper power after every characteristic-zero extension. The proof of these additional conclusions is in Section 7.

## 2. A rational invariant polynomial killing the required flats and faces

If L_s does not contain zero, choose an algebraic affine linear form ell_s with ell_s(0)=1 and ell_s identically zero on ALL of L_s. This is possible because at least one equation of that closed affine flat has nonzero constant term. It also handles a zero-dimensional cell at any nonzero point.

If L_s is linear and a relative facet hyperplane does not contain zero, choose its relative affine equation normalized to constant term1 and extend its linear part to the ambient vector space. Algebraic linear algebra gives an ambient algebraic affine form with value1 at zero and zero on that relative facet. If the relative facet hyperplane DOES contain zero, no factor is needed.

Take all these finitely many normalized forms, close the collection under Galois, multiply it, and then multiply its transforms under the finite Galois-stable W. This gives a rational W-invariant polynomial f_base with f_base(0)=1. It kills every affine L_s not through zero and every nonradial facet of every linear L_s. Empty products are1. Finiteness uses the Galois stability of the geometric W matrix set; no rational representatives of its elements are needed.

Any rational W-invariant multiplier with value1 at zero preserves these properties. A positive rational quadratic damping factor will later suppress the remaining Away linear cells.

## 3. Affine flats not through zero are invisible in the needed coefficient range

For any f divisible by f_base and u on a flat L_s not through zero, f(u)=0. Thus the polynomial f(Tu) has a factor T-1 pointwise on the whole flat. Integration of the finite polynomial gives

    integral_Omega_s rho_s(u)f((1+t)u)^(2k)du_(L_s)
          =O(t^(2k)).                                      (5)

Here O means exact polynomial divisibility. The contribution to Lambda(t^ell p^k) is the coefficient at t^(k-ell) of T^k times (5). It is ZERO whenever

    k>=1 and ell>-k.                                      (6)

In particular it is zero for every unmarked moment and every fixed positive marker, at every k>=1. At the excluded boundary ell=-k the coefficient at order2k can survive; the domain (6) is intentionally strict. This erasure uses no Euler calculation, no dimension bound, no sign cancellation and no global affine-shape assumption.

## 4. Relative Euler identity on the linear strata

For a linear L_s of positive dimension n_s and a homogeneous rho_sj, define

    J_skj(T)=T^(n_s+j)*integral_Omega_s rho_sj(u)f(Tu)^(2k)du_(L_s),
    nu_skj=integral_Omega_s rho_sj(u)f(u)^(2k)du_(L_s).        (7)

Differentiation and the ordinary divergence theorem within L_s give

    J_skj'(T)=T^(n_s+j-1)*integral_relativeBoundary
                  (u dot n)rho_sj(u)f(Tu)^(2k)dS.             (8)

On a radial relative facet, u dot n is identically zero for all T. On a nonradial facet f(u)=0, so the integrand is divisible by (T-1)^(2k). Hence

    J_skj(1+t)=nu_skj+O(t^(2k+1)).                           (9)

Every homogeneous component f_d is W-invariant, because the action is linear. Thus f(Tu) belongs to R[T] and theta extends coefficientwise; every coefficient of f(Tu)^(2k) is an invariant test in (1). Taking the constant term proves the COMPLETE mixed identity

    Lambda(t^ell p^k)
       =sum_(linear L_s,j) nu_skj * binom(k-n_s-j,k-ell)      (10)

for k>=1 and integer ell>-k. The affine-flat contributions vanished in Section 3. For the linear terms, the requested coefficient order is at most2k, so (9) has no remainder contribution. As usual binom(n,d)=[t^d](1+t)^n for integer n, and is zero for d<0. If no affine flats occur, the linear calculation also permits ell=-k.

Put

    B_star=max_(linear L_s)(n_s+b_s).                       (11)

The maximum is nonempty because there is an Origin stratum. Every such n_s>=1; a zero-dimensional linear stratum would be the excluded atom. For k>=B_star, every upper argument k-n_s-j is a nonnegative integer strictly below k. Therefore

    Lambda(p^k)=0 for every k>=B_star.                      (12)

For every fixed ell<=0, the same formula is eventually zero: k becomes large enough for (6), and k-ell exceeds every upper argument. More refined lower-marker statements are given in Section 8.

## 5. The first effective homogeneous degree is positive

For an Origin cell, convexity and zero in its closure imply epsilon*u belongs to the same cell for 0<epsilon<=1. Let j_s be its first nonzero homogeneous degree. Dividing rho_s(epsilon*u)>=0 by epsilon^j_s and letting epsilon tend to zero proves rho_s,j_s(u)>=0 throughout the cell. It is strictly positive on a nonempty relatively open subset.

At effective degree d_star, the only possible Origin contribution is rho_s,d_star-n_s. If n_s>d_star the index is negative and no such term occurs. If n_s+j_s>d_star it is zero. Otherwise it is exactly the first homogeneous piece and is nonnegative. At least one stratum has a nonzero positive contribution. All Origin pieces of effective degree n_s+j<d_star vanish identically.

Only Away LINEAR cells can now interfere. Their compact closures are separated from zero. Choose rational c>0 so small that 0<1-c*q<=1 on all cell closures. If these Away cells are nonempty, let B=max|f_base| on the closures and gamma=max_AwayLinear(1-c*q)<1. Choose an integer M with B*gamma^M<1/4, and set

    f=f_base*(1-c*q)^M.                                    (13)

If there are no Away linear cells use M=0. Affine flats not through zero remain killed identically; their density signs never enter this choice.

After fixing M choose a sufficiently small relative ball INSIDE one positive Origin cell, sufficiently near zero that |f|>1/2 and its first effective homogeneous density is strictly positive on the ball closure. Such balls exist along a relative cone of positive directions. If A_plus is its positive density integral, the total Origin contribution at effective degree d_star is at least A_plus*4^-k.

For k>=B_star and ell=d_star, formula (10) gives

    Lambda(t^d_star p^k)
      =sum_(Origin,n_s+j=d_star) nu_skj
       +sum_(AwayLinear,n_s+j<=d_star)
             nu_skj*binom(k-n_s-j,d_star-n_s-j).              (14)

Every omitted Origin piece of lower effective degree is zero. Higher effective degrees have vanishing binomial coefficient. On Away linear cells, |f|<1/4; their polynomial pieces are integrable and there are finitely many. The complete second sum has absolute value at most C*(1+k)^d_star*16^-k for a finite C independent of k. Thus

    Lambda(t^d_star p^k)
       >=4^-k*[A_plus-C*(1+k)^d_star*4^-k]>0                 (15)

for every k>=K for some finite integer K. No bounded prefix is used. With arbitrary real density coefficients this is an existence argument, not an effective amplifier-selection algorithm.

## 6. One fixed unit at every order and scalar fields

Choose L>=max(B_star,K,1). Expand

    F^m=sum_(h=0)^m binom(m,h)p^(Lm+h).

All exponents lie in the pure and positive marked tails (12) and (15), proving (4) for the SAME unit t^d_star at EVERY order. Each complete marked value belongs to Q because (1), theta and the phase are native rational. It is positive as a real number and remains nonzero in every compatible characteristic-zero extension. No rationality of the individual density coefficients is required.

The proof uses only coefficients of f(Tu)^(2k), all of which lie in the invariant test ring. It never asserts an actual image for an individual non-invariant coordinate or for a homogeneous density piece.

## 7. Direct actual origin certificate and the atom obstruction

Under the geometric evaluation hypothesis, a(T)(xi)=f(0)=1 for every T. Therefore p(xi,t)=1+t^-1. At t=-1/2 the modifier1+p vanishes, its t-derivative is-4 and p=-1 is a unit. On the actual smooth geometric product this is a transverse local parameter. Some prime divisor through the point has valuation one in F, excluding a unit times any proper power. The same point and derivative persist after algebraically closed extension.

At that point

    F=(1+t^-1)^L*(2+t^-1)

has nonzero coefficients2 and2L+1 at adjacent exponents0,-1. The actual Laurent support is primitive. These certificates must be on the actual algebraic source; formal radial coordinates are insufficient.

The excluded atom at zero is a genuine obstruction to THIS circuit. If it had weight w>0, it would contribute w to Lambda(p^k) at every k, because f(0)=1. All other linear strata are eventually pure by (12), and affine flats are erased by (5), so the complete tail would be w, not zero. Likewise it would contribute w*2^m to Lambda(F^m) once all non-atomic exponents are in the tail. This is a limitation of this construction, not an obstruction theorem for arbitrary Mathieu witnesses.

## 8. Scope of lower-marker statements

All fixed integer ell<=0 are eventually zero by Section 4. If EVERY homogeneous piece on EVERY linear stratum with effective degree n_s+j<d_star vanishes, including Away linear strata, then formula (10) proves eventual zero for every ell<d_star. Since ell=d_star is nonzero at every order of F, it is then the least unsafe integer monomial exponent.

Without this GLOBAL effective-degree vanishing assumption, Away linear pieces of smaller effective degree may give smaller unsafe exponents. The local minimum d_star alone does not rule them out. No classification of all exponents above d_star or smallest L is claimed.

## 9. Positive polyhedral-distribution decomposition: statement

The following elementary lemma supplies a separate structural gate. Fix a finite arrangement of real algebraic affine hyperplanes in R^r and its finite polyhedral cell complex. A polynomial face-measure is P(u)du_F on a relatively open face F, where du_F is relative Lebesgue measure, including counting measure when dimF=0. Faces may initially be unbounded.

Let D be a distribution obtained by a finite sum of constant-direction derivatives of polynomial face-measures, followed by multiplication by polynomials. Suppose its restriction to an open polyhedral region U is a nonnegative Radon measure. Refine the arrangement by the boundary hyperplanes of U. Then on U the measure has a finite decomposition into nonnegative polynomial densities on relatively open polyhedral faces. NO normal derivatives of delta measures remain.

If the actual measure has compact support, every retained nonzero face is bounded. If an actual measure on the closure of U gives no mass to its boundary, this decomposition extends by zero across that boundary. Positivity is used only where the actual Radon measure is known; a restriction of a distribution to an open chamber is not multiplied by a discontinuous chamber indicator.

## 10. Finite normal form by ordinary integration by parts

For each affine face F, fix tangent and normal coordinates, taking polynomial coefficients to be constant in the normal directions. Every such D has a finite normal form

    D=sum_(F,alpha) partial_normal^alpha(P_F,alpha du_F).     (16)

To see closure under a constant directional derivative, split the direction into its tangent and normal parts at the affine hull of F. The normal part adds a normal derivative. The tangent part obeys ordinary integration by parts on the polyhedron:

    D_v(P du_F)=(D_v P)du_F
                 -sum_relativeFacets (v dot n)P du_facet.   (17)

This is an equality of distributions tested on compactly supported smooth functions, so no infinity boundary term occurs for unbounded faces. Existing normal derivatives commute with the new derivative. Normals to F remain normal to each facet, so the expression stays of the form (16).

Multiplication by a polynomial also preserves (16), by the normal Leibniz formula. In normal coordinates z, multiplication of partial_z^alpha delta(z) by Q gives a finite sum with derivatives partial_z^(alpha-beta)delta and coefficients (-1)^|beta|binom(alpha,beta)*(partial_z^beta Q)|F. Tangential polynomial coefficients multiply as usual. Refinement into the arrangement's relatively open faces is finite. Algebraic flats and algebraic derivative directions give algebraic affine hulls throughout; arbitrary real density coefficients are permitted.

## 11. Positivity removes all normal derivatives, dimension by dimension

Collect equal faces and equal normal multi-indices in (16). Start with the largest remaining face dimension. At a point in a relative interior F, take a small neighborhood in U avoiding every other remaining face of dimension no larger than dimF. After higher-dimensional measure parts have been removed, the distribution there is supported on F and has the form

    sum_|alpha|<=m P_alpha(y) partial_z^alpha delta(z)dy,      (18)

with y tangent and z normal coordinates.

Suppose m>=1 is the highest nonzero normal order. Test with

    phi_epsilon(y,z)=epsilon^m psi(y)chi(z/epsilon),

where psi,chi are compactly supported smooth functions in a fixed small coordinate box. Their supports stay in one fixed compact set and their sup norms are O(epsilon^m). A Radon measure applied to them therefore tends to zero. In (18), terms of order below m tend to zero, while the order-m terms equal the fixed linear combination of the derivatives of chi at zero against P_alpha psi. Arbitrary prescribed finite jets of chi isolate each alpha. Thus every highest-order P_alpha vanishes as a distribution in y, and hence as a polynomial on F. Descending normal order removes all alpha!=0.

The remaining polynomial P_0 is nonnegative on F, since the residual distribution is a positive measure. Its face-measure is exactly that measure's restriction to the relative interior F. These interiors are disjoint. Subtracting their nonnegative measures leaves the original positive measure restricted to lower-dimensional strata, hence remains positive. This is the ordinary decomposition by Borel strata; it is not subtraction of an arbitrary positive distribution whose remainder might change sign.

Proceed downward through the finitely many face dimensions. At dimension zero the same test removes delta derivatives and leaves only nonnegative atoms. This proves Section 9. Polynomial coefficients that vanish on a face relative interior vanish identically on that affine hull, so their entire face distributions, including their closure supports, vanish; no hidden boundary derivative is lost by the induction.

If a retained face were unbounded while the actual measure had compact support, its polynomial density would vanish on a nonempty relatively open tail. It would therefore vanish identically. Deleting such faces leaves finitely many bounded nonzero strata.

## 12. Effective-degree bookkeeping for derivative-principle inputs

There is a useful precise degree bound. In a normal-form term partial_normal^alpha(P du_F), attach the upper weight

    weight=dimF+degree(P)-|alpha|.                           (19)

A constant derivative lowers this upper weight by at least one: a tangent derivative lowers polynomial degree; a boundary term lowers face dimension; a normal derivative raises |alpha|. Multiplication by a polynomial of degree q raises the upper weight by at most q. In the normal Leibniz expansion, taking beta normal derivatives of the multiplier lowers both its degree and the remaining delta-derivative order by |beta|, preserving this bound.

Consequently, start with polynomial-density strata with dimF+degree(P)<=D. Apply R constant derivatives and multiply by a polynomial of degree at most R. After positivity removes normal derivatives, EVERY nonzero final stratum satisfies

    n_s+b_s<=D.                                             (20)

This bound is about finite polynomial distributions and their positive decomposition, not a claim of full-dimensional radial absolute continuity. In that setting the analytic theorem permits L>=max(D,K,1), and its selected marker exponent satisfies1<=d_star<=D.

## 13. Exact derivative-principle application gates, not a universal subgroup theorem

Suppose an actual compact Hamiltonian H^0-action has an Abelian DH law given by a finite algebraic polyhedral polynomial decomposition. If the torus has a connected ineffective kernel, use its effective quotient: the law is a polynomial density on the corresponding affine moment subspace, and is a face-measure in the ambient Cartan. If an actual zero moment value exists, that ineffective affine subspace passes through zero. The ordinary Section3 DH theorem applies to the effective torus; no original DH Section4 input is needed.

Suppose separately that the invariant h-valued pushed measure gives zero mass to singular H^0-coadjoint orbits. Then CDKW Theorem3.2 applies on the OPEN positive Weyl chamber: the actual radial probability there is p_H times the positive-root directional derivative of the Abelian measure. Its orbit-volume factor has degreeR_H and the derivative order isR_H. The actual radial measure is nonnegative Radon. Sections9-12 therefore give a finite positive polynomial stratification WITH possible lower-dimensional wall measures. Positivity removes normal delta DERIVATIVES, not positive wall measures themselves.

The singular-H support hypothesis gives zero mass to the root-chamber boundary. Thus the actual radial measure is the zero extension of the positive-chamber decomposition; no product of an arbitrary distribution with a chamber indicator is used. Finite algebraic wall arguments and algebraic root directions ensure the resulting flats are algebraic. Compact support makes nonzero cells bounded. If the Abelian initial effective-dimension/degree bound is D, then (20) holds.

An actual zero fiber and no atom at zero then supply the Origin hypothesis of Section1. For a compact connected real-analytic coadjoint orbit with full-support smooth probability, an actual zero point makes zero belong to moment-pushforward support. If the projection is not identically zero, one nonzero analytic coordinate makes its zero fiber measure zero, so there is no atom at zero. If its image contains an H^0-regular point, the nonzero discriminant composed with the analytic projection shows singular-H preimage has measure zero. These are separate statements; neither says that the h-valued map is a submersion or that its radial measure is absolutely continuous in the full Cartan dimension.

Finally an application must supply the ACTUAL rational full-component invariant-test homomorphism theta, native trace compatibility, finite Galois-stable symmetry, invariant positive rational q and actual origin evaluation. The previously proved full-H Chevalley restriction and inverse-Haar projection construction are a possible route. They are not replaced by imaginary Cartan coordinates.

The exact classical source for the distributional derivative input is CDKW, arXiv1204.0741v2, retained SHA88a81ab03ab006491e2cfdad7dc5966125730b1f126ca3e61c9b45422e5a5683, Definitions1.1-1.3/Assumption1.1, Theorem3.2 and complete proof, Corollary3.3 and Remark3.6, within their already reviewed scopes. Its regular-orbit-support hypothesis remains mandatory. The finite positive-distribution decomposition above is proved locally rather than attributed to an uninspected structure theorem.

For example, a lower-dimensional radial wall measure inside the OPEN regular chamber is now compatible with this analytic construction. A projection confined to the singular-H discriminant still fails the stated CDKW gate. No claim is made that all actual subgroups or every abstract rational form passes these gates. No finite-prefix calculation supplies any unrestricted moment or group-application assertion.
