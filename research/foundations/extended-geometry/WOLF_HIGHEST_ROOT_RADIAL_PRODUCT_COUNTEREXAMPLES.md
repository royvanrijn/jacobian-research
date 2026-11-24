# Highest-root Wolf quotients: radial product counterexamples over all characteristic-zero fields

The polynomial endpoint construction is the sheared circuit in [Piecewise densities and compact coadjoint products](PIECEWISE_DENSITY_AND_COMPACT_COADJOINT_ALL_ORDER_COUNTEREXAMPLES.md), Sections 1--4. The rank-one projection and derivative mechanism is classical. The construction below uses an actual rational regular coordinate on the highest-root quotient. Moment identities and nonpower are proved on that quotient.

## 1. Actual rational models and the all-order pair

**Theorem 1 (highest-root quotients and their split-torus products).** Let G_Q be the standard compact rational Chevalley form of a simple adjoint algebraic group of type other than A1. Use the compact Q-Lie algebra construction in [the coadjoint proof](PIECEWISE_DENSITY_AND_COMPACT_COADJOINT_ALL_ORDER_COUNTEREXAMPLES.md), Section 8: its real form has negative Killing form and its Q(i)-complexification is split. Fix a highest LONG root theta. Let U_Q be the compact rank-one subgroup associated to the root pair +/-theta, with Lie algebra u_Q. It is the image of the usual root SL2 after compact Q-descent; finite central kernels do not affect the adjoint sphere used below.

Let T0 be its one-dimensional torus, x=i H_theta a nonzero rational generator, and set

C=C_G(T0)=C_G(x),
K0=C_G(U)^0,
H=K0 U,
Z_Q=G_Q/C, X_Q=G_Q/H.

Section2 proves the full group equalities and reductivity needed for these actual affine quotients. Their compact real orbits are the usual highest-root coadjoint/twistor orbit and the compact Wolf space.

Write < , > for any positive rational invariant form on the compact Lie algebra, for example minus the Killing form, and let pi_u be its rational orthogonal projection onto u_Q. The function

D(gH)=<pi_u Ad_(g^-1)x, pi_u Ad_(g^-1)x>/<x,x>             (1.1)

is an actual Q-regular function on X_Q. On the compact real orbit0<=D<=1, and r=sqrt(D) denotes only its nonnegative REAL random radius. No global square root is asserted in the coordinate ring.

Put

L=(dim X+2)/2=dim Z/2, b=L-1=dim X/2,
T=1+t,
p=T t^-2 [(1-T^2 D)(1-4T^2 D)^2]^2,
F=p^L(1+t^L p), H_marker=t^(-(L-1)).                    (1.2)

These are Q-regular functions on the actual affine product X_Q times SPLIT Gm, and H_marker is a unit. There is a positive rational constant A_X such that, for every m>=1,

Lambda(F^m)=0,
Lambda(H_marker F^m)
 = A_X 18^(2Lm) (2Lm)!/(2Lm+L-1)! !=0.                  (1.3)

The trace Lambda is the normalized invariant trace on X times Laurent constant term in t. The same identities hold after EVERY characteristic-zero base change. Positivity in(1.3) refers to the rational number at the real model; it implies nonvanishing over every such field.

Moreover F is not a unit times a proper power even after algebraic closure. Its nonzero Laurent coefficient functions have exponent differences with gcd1. For every formally real field K, the individual X_K and split Gm factors are Mathieu, while their product is not. This is a theorem about these standard rational models and their scalar extensions, not arbitrary twisted Wolf forms.

## 2. Full subgroup geometry, descent, and affineness

The compact rational root subgroup U exists explicitly: its Lie algebra is spanned overQ by i H_theta, E_theta-E_-theta and i(E_theta+E_-theta). These are the usual compact rational Chevalley generators. The corresponding split root group over Q(i) is preserved by compact conjugation, so it descends toQ. The centralizer groups and identity component K0 are therefore also defined overQ.

For any root alpha orthogonal to theta, alpha+/-theta cannot be a root: theta is a longest root and the squared norm of such a sum is ||alpha||^2+||theta||^2>||theta||^2. Hence its root groups commute with U. Root decomposition gives

Lie C=Lie K0 direct-sum Lie T0.                         (2.1)

The groups K0 and T0 commute. Their product is a closed connected algebraic subgroup of C, and C is connected because it centralizes a torus in a connected group. Equality of dimensions in(2.1) gives the FULL group equality

C=K0 T0.                                               (2.2)

It also proves K0 reductive: K0 is normal in C, its unipotent radical is characteristic in K0 and therefore normal in C, and C is reductive by the torus-centralizer theorem. Thus that radical is trivial. H is the finite-central image of the product of reductive groups K0 times U, since K0 intersect U lies in the finite center of U. Therefore H is connected reductive and closed.

The equality C_G(x)=C_G(T0) is also a full group statement. OverC, fixing x is equivalent to commuting with every matrix exp(sx), s real. These matrices are Zariski dense in the complex torus T0,C. Conversely centralizing T0 fixes x by differentiation. All algebraic subgroups are smooth in characteristic0, so complex-point equality yields scheme equality and descends toQ.

In fact H is the full normalizer of U: every automorphism of the rank-one semisimple group U over an algebraic closure is inner, so N_G(U)=C_G(U)U. Any element of C_G(U) belongs to C_G(T0)=K0 T0, and T0 lies in U. Thus C_G(U)U=K0 U. No unnoticed finite component of the stabilizer changes the quotient.

The highest-root rank-one group is proper when G is not typeA1. Its symmetric quotient has positive dimension. The affine quotients G/C and G/H are smooth geometrically integral by homogeneous-space geometry and Matsushima affineness for reductive stabilizers in characteristic0. The real compact group orbit is Zariski dense and the group pullback is dominant. Same-field normalized Reynolds traces exist, are unique and commute with all field extensions, by the generic normalized-trace lemma in [Standard compact classical products](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md), Section 3.

Formula(1.1) is Q-regular before descent: inversion and the adjoint representation are regular on G, and pi_u is given by a constant rational Gram matrix inverse. Every h in H normalizes u and preserves the form, so pi_u is H-equivariant and the squared norm is unchanged under g->gh. Thus D descends to the actual affine quotient X. Only rational nonzero constants are divided out.

## 3. Deriving the radial law by right Haar invariance

Let K=G(R)^0 and let its probability Haar measure push forward to the actual compact orbits of X and Z. On Z, define the scalar regular height

phi(gC)=<x,Ad_g x>/<x,x>.

This is a normalized component of the coadjoint moment map. The compact coadjoint orbit has dimension2L. Its normalized symplectic Liouville measure is K-invariant and therefore equals the probability pushforward of Haar. The circle T0 acts Hamiltonianly and nontrivially. The Duistermaat--Heckman theorem gives an absolutely continuous law g(v)dv for phi, with polynomial pieces of degree at most L-1.

For any continuous test function a, right invariance of Haar gives the EXACT identity

integral_K a(phi(g))dg
 =integral_K integral_U a(phi(gu))du dg.                 (3.1)

For fixed g, Ad_u x is uniform on the ordinary2-sphere of radius||x|| in the3-dimensional compact Lie algebra u. Its scalar product with the fixed projection pi_u Ad_(g^-1)x is uniform on[-r(g),r(g)], where r(g)=sqrt(D(gH)). Hence, if nu denotes the probability law of r on X,

g(v)=integral_{r>=|v|} (2r)^-1 dnu(r), v!=0.             (3.2)

The same sphere mixture gives the exact TWO-WAY moment bridge, for every integer j>=0,

lambda_X(D^j)=(2j+1)lambda_Z(phi^(2j)),
lambda_Z(phi^(2j+1))=0.

Indeed a uniform variable on[-r,r] has even moment r^(2j)/(2j+1) and odd moment0. Thus the D-moment sequence and the scalar coadjoint moment sequence determine each other. Over the real compact models this also determines both probability laws: compactly supported measures on[-1,1] are determined by all polynomial moments, and r->r^2 is a homeomorphism of[0,1]. Both regular coordinates and normalized traces are rational, so the displayed identities extend to every characteristic-zero field; no probability measure over that field is asserted. As a normalization check, the average of the rank-one tensor Ad_g x tensor Ad_g x is ||x||^2/dim(g) times the identity: it is invariant and self-adjoint, and its eigenspaces would otherwise split the simple adjoint representation. Taking the rank3 orthogonal projection gives lambda_X(D)=3/dim(g). This check is not a premise of the endpoint argument.

No assertion that averaging preserves purity is made. Equation(3.1) only computes the scalar height law. It bypasses any possible real-component issue in a quotient-fiber identification. The usual compact fiber H/C is a2-sphere by(2.2), but that fact is not required to justify the averaging calculation.

The radius is nonconstant. It takes value1 at the identity coset. If it were identically1, projection equality would put the entire adjoint orbit of x inside u. The linear span of a nonzero adjoint orbit is an invariant ideal in the simple Lie algebra g, hence is all of g. This would force g=u and G to have typeA1, excluded. Therefore D is a nonconstant real analytic regular function on the connected compact orbit.

Every level set of a nonconstant real analytic function on a connected real analytic manifold has volume zero. A short proof is useful here: at each zero of a nonzero analytic function some finite-order derivative is nonzero. Take one of least order; a derivative of order one lower vanishes there and has nonzero gradient. The zero set is therefore covered by a countable union of regular hypersurfaces of derivatives, each of measure zero. The impossible alternative of infinite vanishing order would make the analytic function locally, and then globally, zero.

Apply this to D-c for each fixed real c. Thus nu has no atom at any radius, including0 and1.

For0<a<b<1, formula(3.2) gives

g(a)-g(b)=integral_[a,b) (2r)^-1 dnu(r).

Equivalently, on(0,1) as measures,

dnu(r)=-2r dg(r).                                      (3.3)

Since nu is non-atomic, g is continuous there. The DH density is piecewise polynomial with only finitely many knots; continuity removes point-mass derivatives at those knots. Therefore nu has the ordinary nonnegative piecewise-polynomial density

rho(r)=-2r g'(r),                                      (3.4)

of degree at most L-1. There is no residual mass at0 or1, so this density has total mass1.

## 4. The knot set and right-endpoint order

The height phi is critical exactly when Ad_g x commutes with x. Commuting conjugates reduce inside the compact centralizer to a fixed Cartan, so its critical values are the normalized Weyl pairings

<theta,wtheta>/<theta,theta>.

Both theta and wtheta are LONG roots of the same length. Crystallographic root integrality and the Cauchy--Schwarz bound imply

critical values are contained in {-1,-1/2,0,1/2,1}.      (4.1)

Some of these values may be absent. For instance0 need not occur. Extra roots in the polynomial circuit cause no difficulty. It follows that rho has knots contained in{0,1/2,1}; an interval can have identically zero density.

The height phi has its UNIQUE maximum1 at the orbit point x, because

1-phi(xi)=||x-xi||^2/(2||x||^2).

The smooth positive local volume argument in [the coadjoint proof](PIECEWISE_DENSITY_AND_COMPACT_COADJOINT_ALL_ORDER_COUNTEREXAMPLES.md), Section 6, gives

Prob(phi>1-epsilon) ~ C_tail epsilon^L, C_tail>0.

Indeed in smooth local coordinates the squared distance has a positive-definite quadratic leading term; rescaling coordinates by sqrt(epsilon) and applying dominated convergence gives the stated exponent. Compactness and the unique maximum put the entire small superlevel set in that chart.

The rightmost DH density is a nonzero polynomial of degree at most L-1. Its tail asymptotic therefore forces vanishing order L-1 at1. Because that equals its degree bound, the polynomial is EXACTLY

g(v)=C_Z(1-v)^(L-1) near1, C_Z>0.                       (4.2)

By(3.4),

rho(r)=2 C_Z(L-1)r(1-r)^(L-2) near1.                   (4.3)

Consequently its vanishing order is q=L-2, and

A_X=(-1)^(L-2)rho^(L-2)(1)
   =2 C_Z (L-1)! >0.                                  (4.4)

If A_Z=(-1)^(L-1)g^(L-1)(1), then A_X=2A_Z. This endpoint calculation uses the scalar coadjoint law on Z, not a new metric normal form or the invalid uncorrected stationary-phase formula.

## 5. Even-polynomial transfer and the all-order marker

Take the universal EVEN polynomial

f(u)=(1-u^2)(1-4u^2)^2.

It has f(0)=1, a simple zero at the selected right endpoint1, and double zeros at1/2 (and the irrelevant negative partner). Its derivative at1 is -18. It vanishes doubly at every other possible nonzero radial knot. The knot0 is fixed by multiplication byT and needs no root.

Apply the sheared endpoint coefficient theorem in [the piecewise-density proof](PIECEWISE_DENSITY_AND_COMPACT_COADJOINT_ALL_ORDER_COUNTEREXAMPLES.md), Sections 1--4, to the abstract probability functional on R[r] with density rho. Its degree bound b=L-1 gives exponent L=b+1, and Section4 gives q=L-2. That theorem proves purity and

CT_t integral rho(r) t^(-(L-1)) F(r,t)^m dr
 =A_X 18^(2Lm)(2Lm)!/(2Lm+L-1)! >0.

This moment calculation really occurs in the actual coordinate algebra. Every coefficient of p and F is an EVEN polynomial in r, because f is even. The algebra homomorphism R[r^2]->R[X], r^2->D, identifies its probability functional with the actual normalized trace on those coefficients. Therefore(1.3) holds for the actual regular phase(1.2) and its scalar Laurent marker. No square root r is adjoined to the coordinate algebra, and no arbitrary phase or marker is transported.

This step transfers the moment identities only. Absolute nonpower is established separately on the actual quotient in Section6; it is not inferred from nonpower in an abstract interval algebra.

All coefficients of D,p,F and the marker are rational. Normalized trace commutes with scalar extension. Since(1.3) at m=1 is a positive nonzero rational and the remaining factors are rational and nonzero, A_X itself is a positive rational. The all-order identities and nonvanishing therefore remain valid over every characteristic-zero field.

## 6. Actual geometric zero and absolute nonpower

We need a geometric zero of D, not necessarily a zero on the compact real orbit. Such a zero always exists.

Over an algebraic closure, every regular unit on a connected semisimple algebraic group is constant. Indeed restrict the unit to any right coset of any root subgroup Ga. This gives an invertible polynomial in one variable, hence a constant. The unit is invariant under all right root unipotents; they generate the semisimple group, so the unit is constant. Pullback by the dominant quotient map G->X then proves that every unit of O(X) is constant as well.

D is nonconstant, so it is not a unit. Its principal ideal is proper, and the weak Nullstellensatz gives a point x0 of X over the algebraic closure with D(x0)=0. This remains true if the real radius has a positive lower bound.

At this ACTUAL point,

p(x0,t)=(1+t)/t^2,
1+t^L p(x0,t)=1+t^(L-2)+t^(L-1).

For L>=3 its numerator is monic and squarefree: a common root with its derivative would be the nonintegral rational number -(L-2)/(L-1), impossible for a monic integer polynomial. It has no root0 or-1. Choose any root alpha. At(x0,alpha), p is a unit and the t derivative of the second factor is nonzero.

Smoothness of X times Gm therefore gives a prime divisor with valuation1 on1+t^L p and valuation0 onp. The valuation of F is1. It cannot be a unit times a proper power. The same geometric argument works after every characteristic-zero base change.

The phase also has primitive support in the specified Laurent coordinate t. At the same geometric point D=0,

F(x0,t)=t^(-2L)(1+t)^L+t^(-L-2)(1+t)^(L+1).

For L>=4 the coefficients at exponents -2L and -2L+1 are 1 and L. For L=3 they are 1 and 4, because the second summand also contributes at -5. In every case these adjacent coefficients are nonzero. Therefore the nonzero coefficient functions of the actual F have exponent differences with gcd1. In particular no Laurent-unit multiple of F is the pullback of a Laurent polynomial under t->t^r for an integer r>1. This concerns the chosen split-torus coordinate only; it is not a statement about arbitrary morphisms or factorization of the quotient.


The nontrivial Wolf quotients have L>=3. The formal cases L=1,2 of the general endpoint theorem are not needed for these quotients.

## 7. Genuine formal-real factor instability and exceptional cases

Compact real group points are Zariski dense, and dominant pullback identifies the invariant trace with Haar averaging. Every nonzero real regular a on X has positive square trace. The generic finite rational Gram/LDL argument in [Standard compact classical products](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md), Section 13.1, therefore gives lambda_K(a^2)!=0 for every nonzero a and every formally real extension K ofQ, with no whole-field embedding intoR. Thus only0 is pure and the individual factor X_K is Mathieu. The scalar split torus factor is Mathieu over every characteristic-zero field by the ordinary Laurent constant-term theorem, in the all-field form proved in [Standard compact classical products](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md), Section 13.2. Formula(1.3) gives genuine product failure over every formally real K.

The standard compact exceptional Wolf quotients, with customary finite-center notation, illustrate the construction:

G2/SO4 has dimX8, L5, marker t^-4;
F4/(Sp3 Sp1) has dimX28, L15, marker t^-14;
E6/(SU6 Sp1) has dimX40, L21, marker t^-20;
E7/(Spin12 Sp1) has dimX64, L33, marker t^-32;
E8/(E7 Sp1) has dimX112, L57, marker t^-56.

For each, D is precisely the rational projection-norm coordinate(1.1) and the SAME fixed universal polynomial in(1.2) applies with the indicated exponent. These are the specified standardQroot models; no arbitrary exceptional real form or global classification of all homogeneous spaces is asserted.

## 8. Primary antecedents and exact imported scopes

The highest-root geometry is classical. Joseph A. Wolf, *Complex homogeneous contact manifolds and quaternionic symmetric spaces*, J. Math. Mech. 14 (1965), 1033--1047, is retained as [the original author PDF](sources/wolf-contact-quaternionic-1965.pdf), SHA256 `89dd24bf382f53b7fc5900b82d67c821f3c19f80f5a1b7caa20b7910943fdc3d`. The imported scopes are Theorem 4.2 and its complete proof, printed 1042--1043; Theorem 5.4 and its complete proof, printed 1043--1044; and Theorem 6.1 with its complete compact-fibration proof, printed 1044--1045. These identify the compact highest-root quotient and twistor sphere fibration. Section 2 supplies the full algebraic subgroup equalities, reductivity and rational descent directly. The noncompact continuation and any whole-paper classification audit are not imported.

The derivative principle has classical antecedents. Christandl--Doran--Kousidis--Walter, *Eigenvalue Distributions of Reduced Density Matrices*, exact arXiv:1204.0741v2, is retained as [the context PDF](sources/christandl-doran-kousidis-walter-1204.0741v2.pdf), SHA256 `88a81ab03ab006491e2cfdad7dc5966125730b1f126ca3e61c9b45422e5a5683`. The inspected scope is Example 3.1, Theorem 3.2 and its entire proof, Corollary 3.3 and proof, and Remark 3.6, printed 16--17. Their theorem requires an invariant Radon measure concentrated on maximal-dimensional coadjoint orbits. For the present rank-one U, the only singular coadjoint orbit is0 and Section 3 proves the needed zero-mass hypothesis. The theorem must not be applied with ambient G to the singular highest-root orbit Z. The direct Haar mixture and measure differentiation in Section 3 prove the required normalization independently, so this reference is method attribution, not a mathematical premise. No wall-free assertion is imported.

The conventional exceptional names in Section 7 are compared with Corollary 1.3 of Goertsches--Roschig--Stecker, *Revisiting the classification of homogeneous 3-Sasakian and quaternionic Kähler manifolds*, DOI `10.1007/s40879-023-00601-8`, retained as [the names context](sources/goertsches-roschig-stecker-quaternionic-2023.pdf), SHA256 `d26f205a8e37a44541219c10e8e26d5f586000f4f7d283cb6b57382887088b43`. Only that statement, printed 3, is used for nomenclature. No full classification proof is needed or claimed.

The necessary polynomial-density input is [Duistermaat--Heckman's original 1982 paper](sources/duistermaat-heckman-1982-original.pdf), Section 3 Lemma 3.1, Proposition 3.2 and Corollary 3.3, with the complete Theorem 1.1 proof. The exact retained source and hypotheses are recorded in [the piecewise/coadjoint proof](PIECEWISE_DENSITY_AND_COMPACT_COADJOINT_ALL_ORDER_COUNTEREXAMPLES.md), Section 12. Only this density statement is imported; the uncorrected nonisolated stationary-phase claim in Section 4 is excluded. The degree bound is applied to the actual compact Hamiltonian circle, including its finite-isotropy orbifold regular levels.

The same source section records Kerr--Robles, Section 6.3, printed 38--39, for the compact rational Chevalley construction; Milne's torus-centralizer theorem for connected reductivity; and Arzhantsev's Matsushima criterion for affine reductive quotients. The exact physical sources are reused. Normalized trace/base change, rational Gram/LDL and the scalar Laurent Mathieu theorem are used only in their generic scopes in [Standard compact classical products](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md), Sections 3 and 13. No full-automorphism theorem, finite-cover descent of arbitrary phases, or family-specific classical moment calculation is needed.

The endpoint mechanism, rank-one derivative method and Wolf geometry retain these classical antecedents. The result constructs this specific pair overQ and its scalar extensions. It does not assert worldwide priority, a classification of arbitrary twisted forms, or a product theorem for every compact homogeneous variety.
