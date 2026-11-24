# Symmetric normalizers, invariant phases and the EVII quotient

## 1. Scope

Let K be a field of characteristic zero, G a connected semisimple algebraic K-group, theta an involution, and H0=(G^theta)^0. The group-theoretic statements below apply over K and every extension. Descent of functions is conditional on the actual phase and marker being invariant under the indicated finite stabilizer extension. No averaging of a phase is used. In particular an arbitrary function on G/H0 need not descend to the quotient by the full normalizer.

Write N=N_G(H0), Z=Z(G), and let theta_ad be the induced involution on G_ad. We use normalized invariant traces only on actual affine homogeneous quotients for which they exist. For affine reductive homogeneous spaces these are the established normalized Reynolds traces. Section7 uses the specific compact rational EVII coadjoint model; arbitrary twisted forms are not asserted.

## 2. The full normalizer

**Theorem 1.** There are equalities of closed K-subgroups

    N = inverse_image_G((G_ad)^theta_ad)
      = {a in G : a theta(a)^(-1) belongs to Z}.

Here the fixed subgroup on the adjoint group is the FULL fixed subgroup. Moreover N^0=H0, so N/H0 is finite etale.

**Proof.** It suffices first to work over an algebraic closure. Write g=h direct_sum p for the +1 and -1 eigenspaces of dtheta. The Killing form B is nondegenerate, theta preserves B, and B(h,p)=0. Consequently its restrictions to h and p are nondegenerate. If a normalizes H0, Ad(a) preserves h and its Killing orthogonal complement p. Hence Ad(a) commutes with dtheta. Since

    Ad(theta(a)) = dtheta Ad(a) dtheta^(-1),

this is equivalent to a theta(a)^(-1) in ker(Ad)=Z.

Conversely this centrality condition implies that Ad(a) preserves h. The connected groups a H0 a^(-1) and H0 have the same Lie algebra, hence are equal in characteristic zero: their intersection is smooth, has the same dimension as each, and its identity component is an open connected subgroup of each. Thus a normalizes H0. Algebraic groups and closed subgroup schemes in characteristic zero are smooth, so equality on geometric points gives the closed-group equality; it descends to K.

Differentiating the centrality condition gives (1-dtheta)X in Lie Z=0. Therefore Lie N=h. Since H0 is connected and contained in N, N^0=H0. The quotient N/H0 is finite smooth, hence finite etale. This also proves finiteness without assuming that the full fixed group G^theta is connected. QED.

**Corollary 2.** The map

    z:N -> Z,  z(a)=a theta(a)^(-1)

is a K-homomorphism, its kernel is the full group G^theta, and theta(z(a))=z(a)^(-1).

Indeed z(ab)=a z(b) theta(a)^(-1)=z(b)z(a), because z(b) is central. The assertions about its kernel and theta follow directly.

## 3. The precise Cartan-coordinate gate

The Cartan map c:G->G, c(g)=g theta(g)^(-1), is right H0-invariant and satisfies

    c(g a) = z(a)c(g)       (a in N).

Thus a regular function written as F0(g,t)=Ftilde(c(g),t) is right N-invariant if Ftilde is invariant under multiplication by the actual subgroup z(N) of the center. Invariance under the full geometric center is a sufficient condition, and a geometric identity between K-defined functions descends to K.

This condition applies to the constructed phase itself, including every coefficient in its Laurent variable t. It does not follow merely from invariance of some auxiliary coordinate, nor merely from the existence of a non-invariant pure phase upstairs. Also c can identify points of G/H0 when G^theta is disconnected; arbitrary H0-invariant functions need not be functions of c. Direct right N-invariance of an actual phase is an alternative sufficient gate.

## 4. Finite stabilizer extensions and exact traces

Assume X0=G/H0 is affine. Let H be a K-defined closed subgroup with H0 contained in H contained in N. Then H0 is normal in H, H/H0 is finite etale, and its right action on X0 is free. The quotient Y=G/H is affine, and

    pi:X0 -> Y

is an H/H0-torsor, in particular finite etale and surjective. One way to see the scheme statements is to pass to an algebraic closure, take the invariant ring of the free finite group action on the affine scheme X0, and use etale-local triviality and faithfully flat descent. This construction descends over K, including when H/H0 is not a constant group scheme.

Suppose an actual Laurent phase F0 and a marker M0 on X0 times Gm are right H/H0-invariant. They descend uniquely to regular functions F and M on Y times Gm. A scalar Laurent monomial marker t^ell always descends and remains a unit. If lambda_0 and lambda_Y are normalized G-invariant traces, uniqueness gives

    lambda_0(pi^*b)=lambda_Y(b).

Indeed the left side is normalized and G-invariant as a functional of b. Taking the constant term in t gives, for EVERY m>=1,

    Lambda_Y(F^m)=Lambda_0(F0^m),
    Lambda_Y(M F^m)=Lambda_0(M0 F0^m).

Consequently an invariant pure phase and the actual same invariant fixed marker retain every marked value, every zero or nonzero order, and every explicit monomial-marker threshold. This is coefficient descent, not averaging a phase. The assertion concerns these given functions; it makes no global Mathieu-ascent claim.

If the upstairs phase is absolutely not a unit times a proper power, the downstairs phase has the same property: a downstairs unit-times-power factorization would pull back to one upstairs. More precisely an upstairs geometric prime divisor of order one maps to a downstairs prime divisor with the same order because pi is etale (ramification index one). Thus an actual simple-divisor proof descends. Injectivity of pi^* also preserves exactly the nonzero Laurent coefficients and hence Laurent support and radial-primitivity statements.

For compact-real applications, factor positivity on X0 passes to Y by trace restriction and injectivity. For a formally real extension, the same statement follows from the already established rational Gram/SOS argument whenever the model and its trace have that argument. Neither step asserts positivity for arbitrary twisted forms or for noncompact norm-one models.

## 5. Compatible central isogenies

Let pi_G:G->G'=G/C be a central isogeny. To discuss a descended involution theta' one must require theta(C)=C, or explicitly provide theta' with pi_G theta=theta' pi_G. Under this compatibility the image of H0 is (G'^theta')^0, since their Lie algebras agree and both groups are connected. The inverse image of any stabilizer between (G'^theta')^0 and its normalizer lies between H0 and N. The preceding actual-invariance argument therefore handles its quotient.

Without theta-stability, an arbitrary involution need not descend: it can permute central factors of a semisimple product. For the compact Chevalley involution the action on the center is inversion, so every central subgroup is stable. Even in that case a witness still needs its actual central/normalizer invariance proved.

## 6. Exact family gates and the EIII sign obstruction

The preceding statements identify a gate, not a classification of symmetric global forms. A Cartan-coordinate function must actually be invariant under the central action; an arbitrary phase on a connected-stabilizer cover is not automatically such a function. Balanced projector complements can require a new symmetric-CDF phase, and a selected-endpoint coadjoint phase can require a new sign-invariant phase. The concrete balanced AIII/CII and Albert quotient constructions are separately proved in [the balanced-projector and Albert theorem](BALANCED_PROJECTOR_AND_ALBERT_NORMALIZER_QUOTIENTS.md). This proof does not reuse their model-specific results as premises.

The EIII and EVII standard coadjoint models illustrate the sign issue precisely. Their full connected stabilizers H=C_G(x) have one-dimensional central Lie algebra generated by x. A normalizer preserves that line. Since the invariant form has (x,x)!=0, it can send x only to x or -x; the plus case is exactly H. Sign minus occurs exactly when -x is conjugate to x, equivalently when the corresponding coweight -h is in the Weyl orbit of h.

In EIII take the specified minuscule E6 coweight h. The Weyl group acts trivially on the coweight lattice modulo the coroot lattice: each simple reflection changes a coweight by an integral multiple of a coroot. The group P^vee/Q^vee has order detCartan=3, and the chosen h represents a NONZERO class. For example its fundamental-coweight simple-coroot coordinates have the nonintegral diagonal entry 4/3 in the inverse Cartan matrix, so h is not in Q^vee. If -h were Weyl-conjugate to h, its class would equal its negative, impossible for this nonzero order-three class. Thus N_G(H)=H in this specified EIII model. This is a root-coset argument, not a conclusion from only a list of numerical height values.

For EVII the sign minus DOES occur. In the standard E7 realization as the roots in E8 perpendicular to e7+e8, the seven roots

    e1-e2, e1+e2, e3-e4, e3+e4, e5-e6, e5+e6, e7-e8

are mutually orthogonal and span the seven-dimensional root space. The product of their reflections is -1. Thus -h belongs to the Weyl orbit of the specified E7 minuscule coweight. Consequently

    N_G(H)/H = Z/2,

acting on the actual adjoint orbit by x' maps to -x'. This quotient group scheme is constant overQ: its faithful action on the rational central line identifies it with the full subgroup {+1,-1}. This argument does not assert that its nontrivial coset has a Q-rational representative. As a finite subgroup scheme it acts freely on G/H by right translation, so its actual affine quotient is smooth, geometrically integral and finite-etale covered by G/H.

## 7. The actual standard EVII coordinate and its sign-invariant law

Use the [compact rational coadjoint theorem](PIECEWISE_DENSITY_AND_COMPACT_COADJOINT_ALL_ORDER_COUNTEREXAMPLES.md), Sections8--11, for the standard compact Chevalley E7 group, the minuscule coweight of node5 in its stated numbering, x=ih, and H=C_G(x). The model X=G/H is an actual smooth geometrically integral affine Q-variety of complex dimension54. Its rational normalized trace is scalar-compatible and agrees overR with compact Haar. The transported positive coweight pairing is (h,h')_cw=(ih,ih')_compact; no negative Killing form is applied to a split h without this transport.

The established self-height and its symmetric version are

    w=(1+v)/2,    v=(x,Ad_g x)/(x,x).

Both w and v are Q-regular functions on the ACTUAL orbit. The Weyl orbit contains -x, so the nontrivial normalizer action is v->-v and w->1-w. Put

    S=v^2.

This is a regular function on the actual quotient Y=G/N_G(H). The asymmetric endpoint phase on X is not claimed to descend unchanged.

The cited theorem gives w in[0,1] a finite piecewise polynomial density of degree at most26, with critical knots0,1/3,2/3,1. Its unique endpoint maximum has the norm-distance quadratic defect and the positive-volume tail exponent27. Hence its right density has EXACT vanishing order26. The density g(v) of v is symmetric; the law of r=|v| has density rho(r)=2g(r) on[0,1], with knots0,1/3,1, degree at most26, and exact right vanishing order26. Write

    A=(-1)^26 rho_right^(26)(1)>0.

If A_w denotes the corresponding signed derivative of the density of w, then A=A_w/2^26. This follows from g(v)=rho_w((v+1)/2)/2. Density coefficients need not be asserted rational.

The radius r need not be a regular function. Only its EVEN polynomial algebra is used: the actual trace of any polynomial in S equals the radius-law integral of the same polynomial in r^2. No square root is adjoined to the quotient coordinate algebra.

## 8. The invariant phase and every marked order

Let T=1+t and define ACTUAL regular functions on Y times splitGm by

    p = T t^(-2) [(1-T^2 S)(1-9 T^2 S)^2]^2,
    F = p^27(1+t^27 p),       M=t^(-27).

Here M is a unit. The endpoint polynomial f(r)=(1-r^2)(1-9r^2)^2 has f(0)=1, a simple zero at1, a double zero at1/3, and f'(1)=-128. Its negative roots introduce no additional knots inside[0,1].

The precise endpoint theorem in Sections1--4 of the cited coadjoint proof says: for a finite piecewise polynomial density of degree b, selected right endpoint1 of order q<=b, f simple there and double at all other nonzero knots, put L=b+1 and p=(1+t)t^-2 f((1+t)r)^2. Then F=p^L(1+t^L p) is pure, and the SAME fixed marker t^(-(q+1)) has value

    A_q f'(1)^(2Lm) (2Lm)!/(2Lm+q+1)!

at EVERY m>=1, where A_q=(-1)^q rho_right^(q)(1)>0. The negative endpoint part, other-knot remainders and all modified-phase binomial terms are included in that all-order theorem. It is not a finite-prefix assertion or a purity-root converse.

Apply it with b=26,L=27,q=26. The phase uses only r^2=S, so its functional identity is exactly the trace identity in the actual coordinate algebra:

    Lambda_Y(F^m)=0,
    Lambda_Y(M F^m)=A * 128^(54m) * (54m)!/(54m+27)! > 0

for EVERY m>=1. Since F,M and the normalized Q-trace are rational, the marked value at m=1 is rational; dividing by its positive rational factorial factor proves A itself is a positive rational. Thus these are nonzero rational values after EVERY characteristic-zero extension, with no ordering or embedding of the whole extension intoR. Coefficient descent also shows that these are the same values under the finite-etale pullback to X.

## 9. Direct absolute nonpower and primitive support on the quotient

The connected compact orbit contains points with v=1 and v=-1, hence a point with v=0 by continuity. Its image has S=0 on Y. The Q-fiber S=0 is therefore a nonempty finite-type scheme and has a geometric Qbar-point y0. At this ACTUAL point

    p(y0,t)=T/t^2,
    1+t^27 p(y0,t)=Q(t)=1+t^25(1+t).

This polynomial is squarefree. Indeed Q'(t)=t^24(25+26t). A common root cannot be0; at its only other possible value,

    Q(-25/26)=1-25^25/26^26>0.

Also Q(0)=Q(-1)=1. Thus choose any root alpha inQbar: alpha is nonzero, alpha+1 is nonzero, and Q'(alpha)!=0. At (y0,alpha), p is a unit, while 1+t^27 p has nonzero t-derivative. Smoothness of the ACTUAL quotient product gives a geometric prime divisor through that point on which its order is1 and p has order0. Therefore F is not a unit times a proper power over any algebraically closed extension. This proves the assertion directly downstairs; it is not an inference that arbitrary nonpower phases survive a morphism.

The same ACTUAL fiber gives

    F(y0,t)=t^-54 T^27+t^-29 T^28.

Its coefficients at exponents -54 and -53 are1 and27. Hence the actual Laurent support has adjacent exponents and its support differences have gcd1. These simple-zero and support arguments replay after every characteristic-zero extension through an embedding of Qbar in its algebraic closure.

## 10. Factor property, central acting groups and scope

For every formally real extension K/Q, the compact standard source X_K has anisotropic square trace by the established rational Gram/LDL argument. Injective pullback Y_K->X_K preserves square traces, so the quotient factor Y_K also has only the zero pure phase and its normalized trace kernel is Mathieu. The scalar split torus factor is Mathieu over every characteristic-zero field by the scalar Laurent theorem. Their product fails for every formally realK by the explicit pair F,M; the pair itself works over EVERY characteristic-zero field. This is genuine two-factor instability with a specified actual phase, not a global Mathieu-ascent statement.

The center of the acting E7 group acts trivially on its adjoint orbit. Passing through a compatible central isogeny changes neither this actual orbit quotient nor its invariant phase. No arbitrary Q-form of EVII is inferred from a compact real label or from complex failure.

The abstract normalizer theorem concerns any stated semisimple involutive model with its actual invariance hypotheses. The explicit second result concerns this standard rational EVII quotient. The nonsymmetric flag G/T is not covered by the involutive normalizer argument merely because it has a Weyl-group normalizer. Broader family inventories still require the actual coordinate-invariance and compatible-model gates; they are not automatic consequences of either result.

## 11. Mathematical inputs and attribution

The abstract group identity is proved by Killing orthogonality and the central cocycle. The coefficient, trace and divisor descent arguments are written in Sections3--5; the generic normalized-trace restriction agrees with [the homogeneous-quotient theorem](HOMOGENEOUS_QUOTIENT_TRACE_AND_COUNTEREXAMPLE_PULLBACK.md), Sections1--3. A full fixed group must not be replaced by its identity component without the explicit finite quotient argument.

The actual EVII rational model, full connected torus centralizer, affine orbit, exact DH-law hypotheses and local endpoint order are reused ONLY in their stated scopes from [the coadjoint theorem](PIECEWISE_DENSITY_AND_COMPACT_COADJOINT_ALL_ORDER_COUNTEREXAMPLES.md), Sections6,8--11. That proof retains the precise Duistermaat--Heckman, Kerr--Robles, Milne and Matsushima primary scopes. No new primary classification is imported here. Its complete distinguished-endpoint theorem supplies the all-order circuit calculation; the sign-invariant regular S-coordinate and direct downstairs certificates are proved above. The compact factor's formal-real argument and scalar-torus theorem use their generic scopes in [the classical-model theorem](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md), Sections3,13. The method retains the endpoint attribution recorded in those proofs. There is no priority claim for the normalizer or endpoint mechanism and no inference from a finite collection of moments.
