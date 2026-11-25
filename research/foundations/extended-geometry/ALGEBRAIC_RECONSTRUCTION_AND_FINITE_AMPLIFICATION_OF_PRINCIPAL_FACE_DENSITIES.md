# Algebraic reconstruction and finite amplification of principal-face densities

Let G/Q be a connected semisimple algebraic group with compact real points, and let H be a closed subgroup with reductive identity component. Fix x in Lie(G)(Q), and let nu be the radial probability of the actual map g -> projection_h(Ad(g^-1)x) under normalized compact Haar measure, using the positive form -Killing_G. Fix a maximal Q-torus T_H of H0, its full finite geometric Weyl image W, the connected Weyl subgroup W0, and a closed positive chamber C. Suppose nu on C has a finite algebraic polyhedral description with polynomial density degree bounds and realizes native rational moments for every test in Q[u]^W.

Then all density coefficients in an algebraic cell description are real algebraic. Given exact algebraic cells, finite actions and exact native moment access, finitely many such moments reconstruct the coefficients by an invertible algebraic matrix. Under the additional hypotheses of the positive-stratified invariant-test normalizer criterion, exact finite searches also certify amplifier parameters for all later orders. The proof supplies genuine full-component symmetry; Section4 shows why invariant moments alone do not suffice for an arbitrary nonsymmetric measure.

## 1. Every geometric component has a real representative

Let G/R be an affine connected algebraic group with compact real points and a faithful real representation. Average a positive definite form over G(R). Its invariance extends to the algebraic group by Zariski density of the smooth real identity neighborhood. After a real basis change G is a closed subgroup of O_N. Let J be any R-defined algebraic subgroup of G, possibly disconnected.

For g in J(C), g*=bar(g)^transpose belongs to J(C), because J is real-defined and orthogonal. Hence a=g*g is positive Hermitian and belongs to J(C). We claim every real power a^t belongs to J(C)^0. Diagonalize a with positive eigenvalues lambda_i. For any polynomial equation Q vanishing on J, Q(a^t) is a finite sum of real-positive exponentials beta^t with complex coefficients, where beta are products of the lambda_i. It vanishes for every nonnegative integer t because a^n belongs to J. Group equal beta and use the finite Vandermonde matrix at integer t=0,...,s-1: every grouped coefficient is zero. Thus Q(a^t)=0 for ALL real t. The continuous path a^t from identity shows a^t belongs to the identity component. The matrix u=g*a^(-1/2) lies in J(C), in the same component as g, and is unitary. It is also complex orthogonal, so u*=u^transpose and therefore bar(u)=u. Thus u is in J(R).

Consequently EVERY geometric component of J contains a real point. Apply this to J=N_H(T_H), inside an actual compact G/R. Its identity component is T_H, so every geometric Weyl-component matrix has a real representative. Since these matrices were already algebraic over Q, their entries are REALALG. No rational representative is asserted. The finite Galois-stable full matrix set W remains defined over Q.

## 2. Genuine outer symmetry of the actual radial probability

Let W0 be the connected Weyl group of H0, choose its closed positive chamber C, and put

    Omega={w in W : w C=C}.

Every w permutes W0 chambers; unique connected-Weyl correction gives W=W0 semidirect Omega. On the closed chamber, two points are in the same W-orbit exactly when they are in the same Omega-orbit: a connected Weyl orbit meets the closed chamber in one point, including its boundary.

The actual h-valued moment probability is invariant under every real H-component. This follows from normalized group-trace bi-invariance and inverse Haar: translating the group variable by h in H(R) conjugates the projected vector by h and leaves all polynomial moments unchanged. Compact support makes those moments determine the measure. Hence its radial probability on C is Omega-invariant. Section1 supplies real representatives for ALL geometric components used here; geometric invariance is not being substituted for a real action.

In the principal-face setting the principal face and the actual linear support flat are Omega-stable by uniqueness of the face and affine hull. Refine the finite algebraic arrangement by all Omega-translates. It is still finite and algebraic, and Omega permutes its bounded cells. Use the invariant positive algebraic B-metric for relative face volumes. Equivalently, ordinary Euclidean reference volumes have algebraic change-of-volume constants. Degree bounds are preserved.

## 3. Finite-dimensional invariant moment reconstruction

Fix a compact Omega-stable union of finitely many algebraic polyhedral cells in C, with a finite degree bound on each cell. Before forming V, enlarge the degree bound on each Omega-orbit of cells to the maximum of the given bounds on that orbit; this finite enlargement retains every measure admitted by the original bounds and makes the face-measure space Omega-stable. Let V be the real vector space of signed measures having polynomial densities of those degrees on these cells, using relative face volume. Include lower-dimensional faces if present. The cells are disjoint; choose intrinsic algebraic coordinates on each affine hull. Polynomial face-measure basis elements are linearly independent and defined over the real algebraic numbers. Let V_alg be their REALALG-linear span, so V=V_alg tensor_REALALG R. The finite Omega-action preserves V_alg and has algebraic coefficients. Choose a REALALG basis mu_1,...,mu_s of (V_alg)^Omega. Extension of scalars makes the same list a real basis of V^Omega.

Suppose an actual Omega-invariant probability nu lies in V and every native test moment for

    R=Q[u]^W

is rational. These tests separate Omega-orbits in C. Indeed real W-invariant polynomials separate W-orbits, for example by products of B-squared distances to a finite orbit. Moreover the degreewise W-Reynolds operator preserves Q-polynomials: averaging coefficients over W is Galois-fixed. Its image over Q spans the invariant-polynomial space after real scalar extension. Therefore some Q-invariant test separates any two distinct Omega-orbits in C.

The real algebra generated by these rational tests contains constants and separates points on the compact quotient by Omega. Stone-Weierstrass makes it uniformly dense in continuous quotient functions. If an Omega-invariant signed measure annihilates every test, it annihilates every invariant continuous function, then every continuous function by finite-group averaging, and hence is zero. Thus invariant test moments determine V^Omega injectively.

Finite-dimensional linear algebra now selects s rational invariant polynomials P_1,...,P_s for which

    M_ji=integral P_j dmu_i

is invertible. Every entry is REALALG. To see this, triangulate each bounded algebraic polytope into algebraic simplices and integrate its polynomial exactly; simplex Jacobians and monomial integrals are algebraic. Relative metric volume contributes only an algebraic square root of a Gram determinant.

Writing nu=sum_i c_i mu_i gives M c=b, where b_j is the native rational moment. Hence c=M^-1 b is REALALG. Consequently ALL density coefficients in the refined cell description are algebraic. Coefficients in the original cells are algebraic as well, by polynomial equality on any nonempty relatively open subcell.

This proves algebraicity, not rationality. No positivity is needed in the reconstruction after the actual symmetric probability has been supplied.

## 4. Why component symmetry cannot be omitted

Take W={+1,-1} acting on R, C=R, and a transcendental real tau with0<tau<1. The probability density

    rho(u)=(1+tau*u)/2 on(-1,1)

has algebraic polyhedral support and polynomial degree1. Every invariant rational test is even, and its moment is rational when its coefficients are rational. Nevertheless the density coefficient tau/2 is transcendental. It is not Omega-invariant. Thus the generic invariant-test normalizer theorem, which does not assume density symmetry, cannot be upgraded to algebraic coefficients merely from native rational invariant moments. The genuine full-component symmetry proved in Sections1--2 is essential to the present actual-model application.

## 5. Finite exact recovery, with explicit data access

Given the exact algebraic cells, degree bounds and finite W/Omega matrices, compute an algebraic basis of V^Omega. Enumerate rational invariant polynomials by applying the finite Reynolds operator to rational monomials. Compute their algebraic moment rows on the basis measures and retain independent rows until rank s. Section3 proves termination; no uniform degree or polynomial-time bound is asserted. Evaluate the finitely many selected native rational moments and solve the algebraic system.

This is an exact finite procedure PROVIDED those native moments are available exactly. For explicitly presented matrix groups and the actual invariant lift, the usual finite rational linear-algebra route is available: Chevalley lifting in each fixed degree solves a rational restriction matrix; compose with the rational adjoint/projection map; compute normalized group trace on a finite stable rational coordinate-module by Lie-invariance equations and lambda(1)=1. Semisimplicity and transitivity ensure its coinvariant space is the constants. Coordinate-ideal reduction is a separate standard computational-algebra input. No implementation or complexity bound is claimed here.

## 6. Exact amplifier parameters once coefficients are algebraic

For this section additionally assume the hypotheses of the [positive-stratified invariant-test criterion](POSITIVE_POLYHEDRAL_STRATIFIED_INVARIANT_TEST_COUNTEREXAMPLES.md), Sections1--8: finite positive bounded strata, no atom at0, a positive-dimensional Origin piece, and an actual unital invariant-test homomorphism with its compatible native functional. Use its rational invariant face polynomial f_base, effective degree d_star>=1, and B_star=max(n_s+b_s) over the relevant linear strata. Actual-origin evaluation and smooth geometric integrality remain extra hypotheses for the geometric conclusions, independently of density reconstruction.

All relevant polyhedral bounds can now be chosen algebraically. The finite W/Galois orbit wall polynomial is rational. Choose rational c so1-c*q is strictly positive on the compact cell closures. Polynomial maxima on compact algebraic polyhedra are algebraic and can be bounded/compared exactly by real-algebraic decision procedures. Thus a bound B for |f_base| and gamma<1 for Away-cell damping can be chosen algebraically. Search integer M until B*gamma^M<1/4; termination follows from gamma<1. If there are no Away linear cells, take M=0 and omit B and gamma; the negative Away bound below is then C=0.

After M is fixed, choose a SMALL ALGEBRAIC SIMPLEX inside one positive Origin cone where |f|>1/2 and its first effective homogeneous density is positive. Such a simplex exists by relative openness. Its integral A_plus of that homogeneous polynomial is positive algebraic, by simplex integration. Using a simplex avoids importing a possibly transcendental ball-volume constant.

Bound Away homogeneous polynomial pieces by algebraic suprema times algebraic polytope volumes; exact integrals of absolute values are unnecessary. The finite binomial coefficients admit a bound

    C*(1+k)^d_star*16^(-k)

with algebraic C. For k>=2*d_star, the sequence (1+k)^d_star*4^(-k) is strictly decreasing: its successive ratio is at most (1+1/(2*d_star))^d_star/4<1/2, by the binomial-series bound below2. Search integer K>=max(B_star,2*d_star,1) until

    C*(1+K)^d_star*4^(-K)<A_plus.

This exact inequality then certifies ALL later orders. Take L>=max(B_star,K,1). For p=(1+t)t^-1 theta(f((1+t)u))^2, the criterion then gives F=p^L(1+p) and the same fixed unit marker t^d_star at every order. These are finite exact searches with termination proofs; they are not finite-prefix moment tests. Standard exact algebraic-number operations and semialgebraic feasibility/optimization are explicitly the decision-procedure inputs.

## 7. Scope

The algebraicity conclusion applies to actual principal-face radial laws once the genuine real full-component symmetry is established, using the already proved finite algebraic cell and degree bounds plus native rational invariant moments. It is not a theorem for arbitrary nonsymmetric invariant-test measures. The exact recovery/parameter procedure additionally assumes an exact presentation of cells, finite actions and native moment evaluation. There is no complexity, new worldwide-priority or arbitrary-model classification claim.
