# Homogeneous quotient pullback of traces and fixed witnesses

Let K be a characteristic-zero field, G a connected reductive K-group, and L contained in H closed K-subgroups whose identity components are reductive. Assume that the represented quotients

    X=G/L,   Y=G/H

are affine. These are scheme quotients; no transitivity of G(K) on quotient K-points is assumed. Write q:X->Y for the natural map. Both varieties are smooth and geometrically integral. Let lambda_X and lambda_Y be their unique normalized G-invariant traces. A phase f is pure if lambda(f^m)=0 for every m>=1, and safe if for every fixed regular marker a, lambda(a f^m)=0 for all sufficiently large m. Unsafe means not safe. The kernel of lambda is Mathieu if every pure phase is safe.

Write A=O(X), B=O(Y), and identify B with q*B inside A. For the products with the SPLIT multiplicative group, write

    A_t=A[t,t^-1], B_t=B[t,t^-1],
    Lambda_X=lambda_X tensor CT_t,
    Lambda_Y=lambda_Y tensor CT_t.

All the statements below hold after every field extension of K.

## 1. Statement of the bridge

**Theorem.** The map q is smooth and surjective, and q*:B->A is injective. Its normalized traces satisfy

    lambda_X(q*b)=lambda_Y(b).                         (1)

There is a canonical K-linear B-linear retraction E:A->B, commuting with extension of the ground field, such that

    lambda_X(a q*b)=lambda_Y(E(a)b)                    (2)

for every a in A and b in B. Its coefficientwise extension E_t:A_t->B_t gives

    Lambda_X(a q*f^m)=Lambda_Y(E_t(a) f^m)             (3)

for every fixed a in A_t, every f in B_t and every integer m>=0. In particular, for EACH GIVEN downstairs phase f, the following are equivalent in both directions: purity of f and q*f; all-marker safety of f and q*f; and unsafety of f and q*f. A single upstairs marker pushes down by E_t with exactly the same marked value at every order. The pushed-down marker need not be a unit.

Conversely, any fixed downstairs marker h pulls back to q*h with exactly the same marked values. If h is a unit, q*h is a unit. Hence a pure phase with one fixed unit marker nonzero at EVERY positive order pulls back to such a pair. Pullback also preserves the exact nonzero Laurent support in the displayed variable t, so it preserves any proved Laurent-support primitivity or exact Laurent-monomial marker cutoff.

For absolute nonpower a geometric divisor hypothesis is stated separately. Suppose a regular phase F on Y times G_m has a geometric prime divisor D with ord_D(F)=1. Then q*F has a geometric prime divisor with valuation1, and therefore cannot be a unit times a proper power over the algebraic closure. More generally this preservation of a simple divisor holds under ANY smooth surjective morphism of smooth geometrically integral varieties. No claim about nonpower under arbitrary pullback is made.

Finally,

    ker(Lambda_X) Mathieu  ==>  ker(Lambda_Y) Mathieu. (4)

The same implication holds without the Laurent factor. Its converse is not asserted. The given-phase equivalences in (3) concern phases already in B_t, not arbitrary phases of A_t.

## 2. Quotient geometry and trace restriction

For a closed subgroup J of a smooth algebraic group, the represented quotient G/J exists and G->G/J is a faithfully flat J-torsor [M,5.25-5.28]. In characteristic zero these groups are smooth. Consequently the quotients are smooth. Over an algebraic closure the image of the connected, hence irreducible, group is the whole quotient, so X and Y are geometrically integral.

After the faithfully flat base change G->G/H, the map q becomes

    G times (H/L) -> G.

Explicitly, (g,hL) maps to (g,ghL) in G times_Y X; this is well-defined by the quotient sheaf property and is an isomorphism, as can be checked after lifting an X-point to G. The quotient H/L is nonempty and smooth, since H is smooth. Thus this projection is smooth and surjective, and smoothness and surjectivity descend to q. In particular q is faithfully flat. A regular function vanishing after q* vanishes on the target, so q* is injective. Quotients and this diagram commute with ground-field extension. This proves the geometry assertions as scheme statements, including when H or L is disconnected.

There is also a direct tangent description over an algebraic closure: the derivative at the identity coset is the surjection

    g/l -> g/h.

Translation gives the same rank at every point. It is consistent with, and does not replace, the smooth quotient argument.

The composite lambda_X q* is normalized and G-invariant. Uniqueness of normalized trace on the affine homogeneous target gives (1). Equivalently, identify

    O(X)=O(G)^(right L), O(Y)=O(G)^(right H).

These identifications follow by faithfully flat torsor descent. The normalized left-G trace lambda_G on O(G) restricts to lambda_X and lambda_Y, because each restriction is normalized and left-G-invariant. Only the generic same-field Reynolds trace and its uniqueness are used; no all-automorphism trace theorem is needed.

## 3. Ambient right-H Reynolds averaging

The averaging is defined in O(G), where right H acts. It is NOT necessary that right H preserve O(G)^(right L), and generally it does not.

First, H is linearly reductive in characteristic zero, even if it is disconnected. Here is the required finite-component argument. A finite-dimensional rational H-module and H-submodule W admit an H^0-equivariant projection p onto W by complete reducibility of the connected reductive H^0. Over an algebraic closure average the maps h p h^-1 over representatives of the finitely many components H/H^0. The average is independent of the representatives, because p commutes with H^0; it is Galois invariant, because p is defined over K; it restricts to identity on W and is H-equivariant. It therefore descends to the desired H-equivariant projection over K. This proves complete reducibility for H; it is also the proof of [M,Corollary22.43].

Every rational H-module is the union of its finite-dimensional submodules. Its invariant summand and sum of nontrivial simple summands give the canonical Reynolds projection R_H. This projection is natural for H-module homomorphisms: such homomorphisms preserve the invariant part and send the nontrivial simple part into the nontrivial part or zero. Apply it to the right-H rational module C=O(G):

    R_H:C -> C^(right H)=B.

For b in B, multiplication by b is right-H-equivariant. Naturality gives R_H(ba)=b R_H(a), so R_H is B-linear. Left G commutes with right H; hence R_H is left-G-equivariant. The normalized left-G-invariant functional lambda_G R_H must equal lambda_G. The uniqueness assertion applies to the whole O(G), not just to B.

Now restrict the DOMAIN of R_H to A=C^(right L), and denote this restriction by E. Its output lies in B, which is a subalgebra of A because L is contained in H. It is a B-linear retraction and

    lambda_Y(E(a)b)=lambda_G(R_H(a)b)
                  =lambda_G(R_H(ab))
                  =lambda_G(ab)=lambda_X(a q*b).

This proves (2). The projection and the identity commute with every field extension. Indeed invariants of a rational module commute with field extension, as follows from their description as the kernel of the coaction-minus-trivial-coaction map. On each finite-dimensional module the canonical invariant projection is uniquely characterized by equivariance and its restriction to invariants; that characterization survives field extension, using complete reducibility over the enlarged field. Passing to the locally finite union proves the assertion for O(G) and E.

Thus (2) is a relative Reynolds projection formula. It does not require finite fibers or an algebra trace of a finite morphism. In finite-etale situations another construction may supply the same marked identity, but no finite-etale theorem is an input here.

## 4. Exact moment identities and both fixed-phase directions

Extend E coefficientwise to E_t. It is B_t-linear. For a=sum a_i t^i and b=sum b_j t^j, constant term uses only finitely many products a_i b_j, so (2) gives

    Lambda_X(ab)=Lambda_Y(E_t(a)b).

Taking b=f^m proves (3).

With a=1, every pure moment is identical upstairs and downstairs. This proves BOTH purity directions. It also proves that a target fixed marker and its pullback have the same surviving orders, including nonzero values at every positive order.

If f is safe downstairs, take ANY upstairs marker a. The ONE downstairs marker E_t(a) is fixed independently of m. Its marked values eventually vanish by downstairs safety, and (3) makes the upstairs values vanish at exactly those same orders. Thus q*f is safe upstairs. Conversely, if q*f is safe upstairs, test it on pullbacks of each downstairs marker. Equation(1), coefficientwise, gives downstairs safety. This proves BOTH safety directions.

If f is unsafe downstairs, a fixed downstairs marker with infinitely many surviving orders pulls back with the same orders. If q*f is unsafe upstairs, choose an upstairs marker a with infinitely many surviving orders. The fixed marker E_t(a) has those same marked values downstairs. This proves BOTH unsafety directions directly. The identical argument transports a marker surviving at EVERY positive order; no pigeonhole choice of a new marker by order is needed.

Reynolds averaging can destroy invertibility: it is a linear projection, not an algebra homomorphism. For example, for the right action of H=G_m on C=K[z,z^-1], its Reynolds operator is constant term, and sends the unit z to0. Accordingly unit markers are asserted to pull BACK, but arbitrary upstairs unit markers are only pushed down as regular markers. No unit pushdown claim is used.

## 5. Laurent support and smooth simple divisors

If F=sum_i f_i t^i belongs to B_t, injectivity of q* gives

    {i:f_i!=0}={i:q*f_i!=0}.

Consequently the Laurent support is literally unchanged. Any statement that its differences have gcd1 is preserved. Every Laurent-monomial moment Lambda(t^ell F^m) is also unchanged, so an exact unsafe-monomial threshold already proved for the target phase is inherited with BOTH directions of that threshold.

For the divisor assertion, work over an algebraic closure and let D be a prime divisor of Y times G_m with ord_D(F)=1. Shrink to a nonempty open neighborhood U of the generic point of D on which D is a smooth Cartier divisor, F=u a, a is a local defining equation for D, and u is a unit. Such an open exists because the target is smooth and the ground field has characteristic zero. Remove the other zeros of F if needed.

The inverse image of D in (q times id)^-1(U) is nonempty by surjectivity. It is an effective Cartier divisor because q times id is flat. It is smooth over D because the morphism is smooth, hence reduced. At a generic point of any of its irreducible components, its local equation a has valuation1 in the source DVR. The pulled-back u is still a unit. The closure of such a component is a geometric prime divisor D' of X times G_m with

    ord_D'(q*F)=1.

Any unit has valuation0 there, whereas the valuation of an e-th power is divisible by e. Thus q*F cannot be a unit times a proper e-th power for any e>=2. Splitting of the divisor into several components causes no problem; every component on the chosen open has multiplicity1. The same proof works after any further extension of K.

The divisor hypothesis is essential to the stated mechanism. The argument does not infer it merely from an abstract nonpower assertion, nor does it prove nonpower for arbitrary dominant or ramified pullbacks. For instance the smoothness hypothesis fails at0 for A1->A1, x->x^2, and the simple phase y pulls back to x^2.

Even a smooth surjective map can destroy an abstract unit-nonpower when there is no simple-divisor certificate. Here is an example OUTSIDE the homogeneous hypotheses. Over an algebraically closed characteristic-zero field take

    Y_0: y^2=x^3-x, X_0: v^2=u^4-1,
    q_0(u,v)=(u^2,uv).

These are smooth affine integral curves. The map is finite, since u and v satisfy monic equations u^2=x and v^2=x^2-1. Over x!=0 it is the etale quadratic algebra adjoining u with u^2=x and derivative2u invertible; over x^2-1!=0 it is the etale quadratic algebra adjoining v with v^2=x^2-1 and derivative2v invertible, with u=y/v. These two opens cover Y_0. Each algebra has rank2, so q_0 is etale and surjective. In the smooth projective closure of Y_0, the unique point at infinity has ord(x)=-2 and ord(y)=-3. Every regular function is uniquely a(x)+y b(x); its possible nonconstant pole orders are even2deg(a) or odd3+2deg(b), which cannot cancel. Thus no nonconstant regular function has pole order1, and all units are constant, since a nonconstant function and its regular inverse would both have a pole at infinity. If x=c g^r with c a unit and r>=2, its pole order forces r=2 and g to have pole order1, impossible. Hence x is an absolute unit-nonpower on Y_0, while q_0*x=u^2 is a proper power. Its zero divisor downstairs has multiplicity2, so it does not satisfy the simple-divisor hypothesis above.

## 6. Global Mathieu descent and its precise direction

If ker(Lambda_X) is Mathieu and f is pure downstairs, then q*f is pure upstairs. Every pulled-back downstairs marker is eventually annihilated upstairs. Its marked values equal the downstairs ones, so f is safe. Thus ker(Lambda_Y) is Mathieu. The argument with no Laurent variable is the same.

This is the standard preimage property of Mathieu subspaces, applied to the algebra embedding q*:B_t->A_t and the identity

    (q*)^-1 ker(Lambda_X)=ker(Lambda_Y).

For completeness, for ANY algebra homomorphism phi:B->A and Mathieu subspace M of A, if b^m lies in phi^-1(M) for all m>=1, then phi(b)^m lies in M. For each fixed c in B, the Mathieu property gives phi(c)phi(b)^m in M eventually, which means c b^m lies in phi^-1(M) eventually. This is the entire proof. It is classical: see Zhao [Z,Proposition4.9(b) and Remark4.10], restated for all sidedness choices in [Z2,Proposition2.5]. The quotient-by-an-IDEAL iff in [Z2,Proposition2.7] is a different theorem; a geometric quotient pullback here is an embedding of coordinate algebras, not that surjective algebra quotient.

There is no converse deduced from fixed-phase safety. Arbitrary phases of A_t may have coefficients outside B. In particular, the point has a Mathieu trace kernel, while the known compact-factor-times-split-torus counterexamples supply homogeneous affine sources mapping smoothly to a point with non-Mathieu trace kernel. The theorem does not descend arbitrary phases, nor assert global Mathieu invariance under homogeneous quotient maps.

## 7. Compact rational applications and genuine factor instability

Now assume G,L,H and the quotients are defined over Q, and G(R) is compact. Let K be any formally real characteristic-zero field. For any nonzero real regular function a on X_R, its pullback to G_R is nonzero. Smooth real group points are Zariski dense, and the compact Haar integral of the nonnegative nonzero square is strictly positive. This Haar functional is the normalized invariant trace on X_R by uniqueness. Thus lambda_X,R(a^2)>0.

For any finite Q-independent set a_1,...,a_s in O(X_Q), its rational trace Gram matrix is positive definite overR. Rational LDL elimination diagonalizes it with positive rational pivots. A positive rational d=a/b, with positive integers a,b, is the sum of ab copies of (1/b)^2. Therefore after extension to a formally real field, lambda_X,K(a^2) is a finite sum of squares after an invertible rational linear change. Vanishing would force every summand zero, and hence a=0. This uses the definition of formal reality, not an embedding of K intoR or a compactness assertion about its K-points. Every individual nonzero element lies in a finite such rational span. Consequently the only pure phase on X_K is0, and its trace kernel is Mathieu. The scalar split-torus factor is independently Mathieu over any characteristic-zero field by scalar Laurent constant-term theory.

Suppose one of the stated target models Y_Q has a displayed rational pure pair F,H on Y_Q times G_m, with H a unit and every marked order nonzero, together with a geometric simple-divisor certificate and primitive Laurent support. Then the SAME formulas in the pulled-back target coordinates give all those conclusions on X_K times G_m over EVERY characteristic-zero field K. For every formally real K, the individual factors X_K and G_m have Mathieu kernels but the product does not. This is genuine factor instability for the specified rational homogeneous source, not a general statement about arbitrary forms.

This conditional application includes L=1, so the group G itself times the split torus inherits every such quotient witness. It also includes every closed reductive L contained in the target stabilizer, whether the resulting quotient is spherical or nonspherical. No spherical classification is needed.

## 8. Applications to the explicit quotient constructions

The bridge applies to each EXISTING construction only with its own stated field, group, and model hypotheses:

* The compact coadjoint-orbit theorem overR pulls back through every actual affine G/L->G/H for that same real acting group and nested reductive stabilizers. It remains a real theorem unless a rational model and rational formula have separately been supplied.
* The standard compact Chevalley rational coweight-orbit constructions, including their EIII, EVII and regular G2 cases, supply rational target formulas. Their same all-order marked values, simple-divisor certificate, and primitive radial support pull back over every characteristic-zero extension ofQ.
* The highest-root Wolf quotient constructions use their actual regular squared-radius coordinate. Their direct smooth-point valuation1 certificate, not an abstract square-root radius, is what is transported here.
* The compact Cayley-plane projector and compact unitary EIV quotient constructions use their stated positive Albert rational models. Their t=1 divisor has valuation1, so the strong nonpower certificate pulls back directly. For the Cayley-plane phase, at the actual projector P=e1 its diagonal coordinate is1 and its polynomial CDF has value v=1. Thus F specializes to t-1, so its coefficient functions at exponents0 and1 are nonzero. This supplies its radial-support assertion locally, independently of any existing support-node label.
* The reduced Albert norm-one construction, when its specified K-defined frame exists, also pulls back by the algebraic theorem under its same E6 acting group. Its all-field pair transports; no compact-factor Mathieu assertion is added for noncompact norm-one forms.

For PWC and Wolf phases the existing direct nonpower proofs exhibit a smooth point with a simple zero of the modifying factor and nonzero remaining factor. This supplies exactly a prime divisor of valuation1, even when that divisor is not globally irreducible. For CDF phases the t=1 certificate is immediate. The present theorem constructs no initial bad pair on a target where none has been proved.

No finite quotient phase averaging, general morphism functoriality, generic preservation of nonpower, or converse global Mathieu implication is inferred. The finite-component average in Section3 constructs a Reynolds projection on rational modules; it does not replace a phase by an average or claim purity of such an average.

## 9. Primary inputs and antecedent scope

[M] J. S. Milne, Algebraic Groups (2017), [exact retained author PDF](sources/milne-algebraic-groups-2017.pdf), SHA256 503f1a7d511632af05541c52fb2bf6fc36f5a75b0323f2f9110f03dcce11fd67. Quotient existence/torsor/fpqc statements:5.25-5.28, printed105-106. Complete reducibility and the finite-component proof:22.42-Corollary22.43, printed477-478. Locally finite semisimple modules and canonical Reynolds projection:12.53 and12l, printed248-249. The proof above gives the finite-component and relative projection arguments explicitly. The universal same-field normalized-trace and rational Gram/SOS lemmas are also recorded in [the standard compact-product proof](STANDARD_COMPACT_CLASSICAL_PRODUCT_COUNTEREXAMPLES.md), Sections3 and13.1. Its Section13.2 records the scalar split-torus Mathieu theorem over every characteristic-zero field, from Duistermaat and van der Kallen's scalar Laurent theorem ([retained primary PDF](sources/dvk-powers-1998.pdf)); only the individual torus-factor assertion is used here. No new whole-book theorem or coordinate-algebra-valued Laurent theorem is imported.

[Z] Wenhua Zhao, Generalizations of the Image Conjecture and the Mathieu Conjecture, J. Pure Appl. Algebra214(2010),1200-1216, arXiv:0902.0212v3. The exact primary preprint, https://arxiv.org/pdf/0902.0212v3, Proposition4.9(b), full proof and Remark4.10, printed29-30, proves the classical inverse-image property. Its entire short proof is reproduced in Section6. No other conjecture or integral-positivity assertion of that paper is imported.

[Z2] Wenhua Zhao, Mathieu subspaces of associative algebras, J. Algebra350(2012),245-272, DOI10.1016/j.jalgebra.2011.09.036, arXiv:1005.4260v2, https://arxiv.org/pdf/1005.4260v2. Proposition2.5 and Corollary2.6, printed11, restate preimage/intersection; Proposition2.7 and its full proof on the same page concern a surjective algebra quotient by an ideal. These narrow comparisons distinguish the valid implication from an inapplicable quotient iff. They are antecedent credit; the proof above is independent of an external Mathieu-subspace lemma.

The trace, projection, divisor, and witness assertions here are written proofs. No bounded moment experiment, finite group surrogate, or executable checker is a proof premise. No novelty or priority claim is asserted for Reynolds averaging or the classical preimage property.
