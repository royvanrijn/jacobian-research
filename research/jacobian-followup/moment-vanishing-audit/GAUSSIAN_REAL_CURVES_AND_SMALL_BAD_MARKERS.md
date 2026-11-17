# Gaussian real curves and small fixed bad monomials

Date: 15 November 2026. Internal written deductions; independent review, external review and formal verification remain pending. The accompanying exact controls do not prove the universal assertions. Mathematical status is in the live claim registry.

This audits Section 2 of the [forwarded draft](../cross-programme-connections/PROOF_DRAFT.tex). Its normalization argument extends the algebraic part of [the compact normalization note](COMPLEX_HEIGHT_IMAGES_AND_CAUCHY_RIGIDITY.md), with a separate proof below because that theorem's full hypotheses include a compact core. The [finite-module marker lemma](FINITE_CURVE_MARKERS_AND_POSITIVE_LOOP_MODEL.md) then applies to the Gaussian restricted pairing. No novelty or exhaustive literature coverage is asserted.

## 1. Polynomial normalization of a real curve image

Let A=C[X_1,...,X_n], n>=1, and let sigma conjugate coefficients while fixing every X_i. For nonconstant P, the following are equivalent:

1. P(R^n) has empty interior in C.
2. All two-by-two minors of the real differential of (Re P,Im P) vanish everywhere.
3. Every polynomial minor P_i(sigma P)_j-P_j(sigma P)_i is zero.
4. P and sigma P are algebraically dependent.
5. P=p(h), where h is a nonconstant polynomial with real coefficients and p is a nonconstant complex univariate polynomial.

For (1) to (2), a point with real differential rank two gives an open image near that point by the submersion theorem. The minors in (2) are real polynomial functions; real Zariski density and the invertible linear change between (Re P,Im P) and (P,sigma P) give (3). In characteristic zero the rank of the differentials in C(X_1,...,X_n) is their transcendence degree, giving (4). One can see this differential criterion directly by extending a transcendence basis: separability identifies the differentials of the basis as linearly independent in the rational function field.

For (4) to (5), set A0=C[P,sigma P], a one-dimensional affine domain, and let B be its normalization in Frac(A0). Normalization is finite, so B is again an affine normal curve ring ([Stacks, Lemma 33.27.1](https://stacks.math.columbia.edu/tag/0BXQ)). Since every b in B is integral over A0 and lies in Frac(A), normality of A embeds B in A. Therefore B^* is contained in A^*=C^*.

The inclusion B in A gives a morphism A^n to Spec B. Choose an affine line on which P is nonconstant; its composite A^1 to Spec B is nonconstant and hence dominant. It extends to a nonconstant morphism P^1 to the smooth projective completion Cbar of Spec B. Properness of the target extends the rational map at the missing points of the nonsingular source. In characteristic zero this map is separable. Riemann--Hurwitz gives -2=e(2g(Cbar)-2)+R with e>0 and R>=0, forcing g(Cbar)=0 ([Stacks, Section 53.12](https://stacks.math.columbia.edu/tag/0C1B)). Over C a smooth projective genus-zero curve is P^1: a degree-one divisor gives an isomorphism by [Stacks, Lemma 53.10.2](https://stacks.math.columbia.edu/tag/0C6L).

Thus Spec B is P^1 with a finite nonempty set removed. It cannot be all of P^1 because it is affine and nonconstant. If at least two points were removed, a rational function with a zero at one and a pole at another would be a nonconstant unit of B. Therefore exactly one point is removed and B=C[h0].

The ring B is sigma-stable, since A0 and its fraction field are. Write sigma(h0)=alpha*h0+beta, where alpha!=0. The involution gives bar(alpha)*alpha=1 and bar(alpha)*beta+bar(beta)=0. Choose c!=0 with bar(c)*alpha=c: use c=1+alpha if alpha!=-1, and c=i otherwise. Then bar(c)*beta is purely imaginary. Put d=bar(c)*beta/2, so d-bar(d)=bar(c)*beta. The coordinate h=c*h0+d satisfies sigma(h)=h. Hence h is a real polynomial, B=C[h], and P=p(h).

Finally, (5) implies (4), and the polynomial parametrized real path p(R) has no interior: its real and imaginary coordinates are algebraically dependent, so its image lies in the zero set of a nonzero real polynomial in two variables. This also proves (1). The image is semialgebraic, so empty interior here is equivalent to real dimension at most one. Constants have rank zero and are treated separately. No compactness is used.

## 2. A fixed bad monomial of degree at most deg P minus one

Let L be the ordinary positive real Gaussian expectation, with nonsingular real covariance (standard independent variance-one coordinates suffice). Suppose P!=0, L(P^m)=0 for every m>=1, and P(R^n) has empty interior. Then there is a FIXED nonconstant monomial X^a, 1<=|a|<=deg P-1, such that L(X^a P^m) is nonzero for infinitely many m. In particular P fails eventual mixed safety.

Proof. A pure-null constant must be zero, so Section 1 supplies P=p(h), with h real nonconstant. On B=C[h] the pairing L(ab) is nondegenerate. Indeed for b=q(h)!=0, choose a=bar(q)(h); then L(ab)=L(|q(h)|^2)>0. Nonzero polynomial functions cannot vanish on all R^n, and the Gaussian has positive density on every open set. Every polynomial integral is finite.

Put d=deg p and H=deg h. The ring B is free over C[P] with basis 1,h,...,h^(d-1). Division by the relation p(h)-P gives generation. Independence follows because the highest h-degrees of distinct summands h^j r_j(p(h)) are distinct modulo d.

If every basis generator had an eventually zero tail, choose M beyond all their cutoffs. For any b=sum_j h^j r_j(P), every L(P^M h^j P^k), k>=0, vanishes, so L(P^M b)=0. Nondegeneracy gives P^M=0, impossible in B because P!=0. This is exactly GFINITEMODULEMARKER20261115 with its restricted-pairing hypothesis verified. The marker 1 is excluded by purity; therefore d>=2 and some fixed h^j, 1<=j<d, has infinitely many nonzero marked moments.

Highest homogeneous parts give deg P=dH and deg h^j=jH<=deg P-H<=deg P-1. Expand this ONE h^j as a finite monomial sum. If every monomial occurring in it had an eventually zero tail, their sum would also. Thus one fixed monomial of that degree bound has infinitely many nonzero moments. Its degree is positive because the constant marker is pure-zero. This is not a claim of nonvanishing at every index, nor a construction of a nonzero pure-null thin-image phase.

The same argument applies to every positive measure with full support on R^n and finite polynomial moments. Positivity alone does not rule out all pure thin-image phases for arbitrary other supports; the retained positive polynomial-loop model is a boundary, not a Gaussian realization.

## 3. Safe classes and finite witness selection

Every nonzero pure-null Gaussian phase that is eventually safe has image with interior. Equivalently P and sigma P are algebraically independent and at least one minor in Section 1 is nonzero. Hence n>=2. In particular this applies to every nonzero pure-null ternary cubic by G3CUBICALL20261115, and to all admitted polynomial Gaussian shear/orthogonal orbits of those cubics, even in unbounded degree, by GMCCUBICORBIT20261115. See the [shear transfer note](GAUSSIAN_SHEAR_NORMAL_FORM_AND_TRANSFERS.md). The corollary uses safety preserved by those maps; it does not assert invariance of real-image geometry under a complex shear.

The converse fails: the retained quartic GMC counterexample has image with interior. Rank two therefore filters geometric architectures; it does not certify safety. No all-degree Gaussian classification, improved cubic prefix, or general vanishing theorem follows.

For coefficients in an exactly presented number field stable under conjugation, one can test all minors exactly. Given all-order purity and a nonzero thin-image input of degree D, enumerate the finitely many monomials of degrees 1 through D-1 and decide polynomiality of each L(X^a exp(tP)) by the retained fixed-amplitude holonomic polynomiality theorem MVHOLPOLY20261115. A nonpolynomial series is a bad fixed monomial; Section 2 guarantees that at least one exists. This is a theoretical terminating witness selection algorithm under those hypotheses. It does not require computing the normalization h. The general holonomic endpoint algorithm and this witness selection procedure are NOT implemented here; finite moments cannot certify an infinite bad tail.

The accompanying checker implements ONLY the exact minor/rank filter over Q(i), coefficient conjugation fixing the coordinate symbols, and the conditional marker-degree bound. It tests constants, unary phases, nonlinear real-seed compositions, safe shear examples and an interior-image failure control. It never treats a thin nonpure input as a GMC counterexample. It does not compute h, test all-order purity, decide safety or select an actual infinite-tail witness.

## 4. Source and assurance boundary

The [source audit](gaussian_real_curve_source_audit.json) binds this proof, its [checker](verify_gaussian_real_curves.py), the [small-control receipt](gaussian_real_curve_verification.json), retained inputs and primary Stacks pages. Section 2 of the forwarded unregistered draft is adopted only in the precise Gaussian scope above; its other sections remain separate owners' work. The small monomial strengthening is proved by the finite-support argument, not by enumerated moments. All universal/topological conclusions remain written internal proofs pending independent review.
