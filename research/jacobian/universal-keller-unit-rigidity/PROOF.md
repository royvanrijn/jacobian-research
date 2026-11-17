# Primitive units close the nonconstant-unit intermediate branch

15 November 2026. Two written general claims with selected exact controls. No independent whole-proof replay, formal verification, external review or novelty certificate. OP-KDESC remains open. This strengthens the [intermediate-surface constraints](../universal-keller-intermediate-surfaces/PROOF.md): a PROPER intermediate has no nonconstant units, even geometrically.

Work in characteristic zero. Put

    t=1+xy, Q=t^2*z+y^2*(1+3t), sigma=t*Q.

Consider ANY whole polynomial Keller map

    F=(u,sigma,B,C)

with any number of retained identity variables u. B,C need NOT belong to the native module or satisfy the native source jets. Every native map is covered. For ANY whole polynomial self-map factorization F=E composed G, set

    p=(u,sigma), q=(E_u,E_sigma), K0=k(u,a),
    U_eta=Spec R, R=K0[S,Y,r^(-1)], r=1-SY=1/t,
    V_eta=Spec A, A=O(V_eta) subset R,
    M=Frac A.

UAINTERSURF1 supplies the actual generic etale diagram U_eta->V_eta->A2 for native maps. Its smoothness, geometric integrality, cofinite-image and localization arguments use only the retained source tuple, so they apply to the broader F above. The carried A1 and the plane target-factor criterion use the native witness separately.

## 1. Root descent is a published input, not a new general theorem

A whole Keller polynomial subring is ROOT CLOSED:

    H(P)^m in k[G] implies H(P) in k[G], m>=1.          (1.1)

This property belongs to the established square-free/root-closed literature. [de Bondt--Yan, Irreducibility properties of Keller maps, version2 (2014)](https://arxiv.org/html/1304.0634v2), Corollary2.2, states square-free preservation. [Jedrzejewicz--Zielinski, Analogs of Jacobian conditions for subrings, version1 (2016)](https://arxiv.org/html/1601.01508v1), Theorems2.4,3.4,3.6, relates the Jacobian condition, square-factorial closure and root closure. We give the short polynomial-ring argument needed here and claim no novelty for (1.1).

Suppose P(G)=H^m. Factor P into distinct irreducible target polynomials, with exponents e_i. Their pullbacks are square-free and pairwise coprime: the product of any two distinct factors is square-free and remains so. Unique factorization in the SOURCE polynomial ring forces every e_i divisible by m. Hence P=lambda*P0^m. Unique factorization again makes H a scalar in k times P0(G). Thus H belongs to k[G].

Together with UAQFSAT1 this gives the useful rational-field form

    H(P)^m in k(G) implies H(P) in k[G]                 (1.2)

for a SOURCE POLYNOMIAL H. First apply polynomial rational saturation to the power, then (1.1). This is not factorial closure: the target polynomial a pulls back under the known cubic to sigma=t*Q while t itself is not recovered by that cubic. A product can descend without its separate factors.

One must apply (1.1) on the WHOLE polynomial affine space. An arbitrary etale map on an open surface can have a nonprimitive unit pullback, as the all-degree examples below show.

## 2. UAPRIMUNIT1: a primitive unit makes a cofinite etale surface map invertible

Let K be ANY characteristic-zero field, put

    U=Spec K[S,Y,(1-SY)^(-1)], r=1-SY.

Let g:U->V be a dominant etale map to a smooth geometrically integral AFFINE surface, whose image misses at most finitely many points. If some unit v in O(V)^* satisfies

    v composed g=r,                                   (2.1)

then g is an isomorphism. Multiplying v by a scalar or taking its inverse gives the equivalent primitive-unit versions. No prescribed degree or assumed affine-plane structure of V is used.

### The generic unit fibre is a split multiplicative curve

Work first after an algebraic closure of K; the final isomorphism descends. Over H=K(r), the generic source fibre of r is Gm with coordinate S and Y=(1-r)/S. The generic fibre W of v on V is a smooth geometrically integral affine curve. The restricted g is a SURJECTIVE etale map

    Gm_H -> W.

Indeed any point omitted by g has unit value algebraic over K, so no such point remains on the generic r fibre. Geometric integrality follows by the dominant injection into H[S,S^(-1)] and flat scalar extension.

Complete the curves to obtain a nonconstant finite map P1_H->Wbar. The classical Riemann-Hurwitz formula gives genus zero for Wbar, and the image of S=1 gives an H-rational point. Thus Wbar is P1_H by the standard genus-zero criterion. Both facts are recorded in [Stacks Riemann-Hurwitz](https://stacks.math.columbia.edu/tag/0C1B) and [genus zero](https://stacks.math.columbia.edu/tag/0C6L).

The boundary Wbar-W is nonempty because W is affine. Every boundary point has all preimages among0 and infinity, since the map from Gm lands in W. Therefore there are at most TWO boundary points, both H-rational (they are images of0 or infinity).

If W were A1, its coordinate would pull back to a Laurent polynomial f(S) with derivative a Laurent unit. In characteristic zero this forces

    f(S)=c*S^e+b, c!=0, e!=0:

the derivative has only one monomial, and an S^(-1) derivative term cannot be integrated to a Laurent polynomial. This map misses b, contradicting surjectivity. Hence W is Gm_H, with a coordinate pulling back as c(r)*S^e. Its degree is d=abs(e).

### A cyclic deck action cannot be etale at the retained node

The total function-field extension is H(S)/H(S^d), up to the coefficient c(r) and inversion. Over the algebraically closed K, its cyclic deck action is

    (S,Y)->(zeta*S,zeta^(-1)*Y), zeta^d=1.              (2.2)

These are REGULAR automorphisms of the ENTIRE U, fixing r. Since they fix the field of V, g composed with each is g as an actual morphism. At o=(0,0), every such automorphism fixes o and acts on its tangent space by diag(zeta,zeta^(-1)). Differentiating g composed h=g there gives

    Dg_o * diag(zeta,zeta^(-1))=Dg_o.

Etaleness makes Dg_o invertible, so zeta=1. Therefore d=1.

Now g is birational and etale, hence an open immersion by the usual birational quasi-finite form of Zariski's Main Theorem. Its complement is finite. Normal affine Hartogs extension gives O(V)=O(U), hence g is an isomorphism. The height-one intersection fact is [Stacks tag031T](https://stacks.math.columbia.edu/tag/031T). All reasoning is geometric, so faithful-flat descent returns the isomorphism over the original K.

The node o must remain in the actual source. Removing it would remove the fixed-point contradiction. The image condition and PRIMITIVITY in (2.1) are also genuine hypotheses.

## 3. UAUNITFACT1: nonconstant intermediate units occur exactly at an invertible first stage

Return to the whole factorization F=E composed G. Then the following are equivalent:

1. G is a polynomial automorphism.
2. V_eta has a nonconstant unit over K0.
3. V_eta has a nonconstant GEOMETRIC unit.
4. t belongs to the embedded field k(G).
5. E_sigma is reducible on the WHOLE intermediate affine space.
6. The zero fibre of E_sigma has exactly two geometrically irreducible components.

When these hold, there are whole polynomial A1,A2 on intermediate space with

    A1 composed G=t, A2 composed G=Q,
    E_sigma=A1*A2.                                    (3.1)

The factors are irreducible and comaximal, and their generic unit pullback generates the FULL source lattice Z. No proper finite-index unit sublattice occurs for a whole factor. Consequently every PROPER factorization, with neither whole stage invertible, has

    O(V_eta_bar)^*=K0_bar^*.                            (3.2)

This holds after arbitrary field extension and arbitrary identity stabilization/intermediate mixing. No native-jet, nonnegative-r or degree hypothesis occurs.

### Units force primitive whole-source recovery

The source unit group is R^*=K0^* times r^Z. If v in O(V_eta)^* is nonconstant, its pullback is c(u,a)*r^d, d!=0. Invert v if necessary. The field k(G) contains all retained values (u,sigma), hence c. Thus a POSITIVE power of t belongs to k(G). Equation (1.2) makes t itself a whole polynomial in G. Since Q=sigma/t is a SOURCE polynomial rational in G, UAQFSAT1 makes Q a polynomial in G as well. Dominance gives (3.1).

In particular the intermediate generic unit A1^(-1) pulls back to r EXACTLY. UAINTERSURF1 gives g=G_eta etale with cofinite image. UAPRIMUNIT1 makes G_eta invertible; the whole G is therefore birational and Keller, hence a polynomial automorphism by the criterion used in UAALLATOM1. This proves2 implies1, including every stabilization.

Conversely a whole automorphism G makes its generic fibre map an isomorphism U_eta->V_eta, so V_eta has the source nonconstant unit. Condition4 directly gives the same argument without first taking a power.

A geometric unit descends far enough for this proof. Its coefficients and inverse lie over a finite extension K1/K0. Its field norm is a unit of O(V_eta); if the geometric pullback is c*r^d, that norm pulls back to Norm(c)*r^(d*[K1:K0]), still nonconstant. Thus3 implies2. The source unit lattice is split and fixed by this Galois action; cancellation of the exponents is impossible.

### The complete reducible-fibre pattern on intermediate affine space

The source polynomials t and Q are geometrically irreducible, and sigma-c is geometrically irreducible for EVERY nonzero constant c. This is the primitive linear-in-z argument retained in UASIGFIELD1.

If E_sigma factors, dominance pulls its nonconstant factors to nonconstant polynomials with product t*Q. Hence there are exactly two factors, pulling back to scalar multiples of t and Q; rescale to obtain (3.1). Reducibility therefore gives condition4 and G invertible. Conversely an invertible G carries the two source zero components isomorphically to two components of E_sigma=0. Condition6 implies the same conclusion by applying the factor argument after an algebraic closure of k and descending the unique inverse.

The two factors are comaximal: E_sigma is a Keller component, so its gradient is unimodular. If E_sigma=A1*A2, every gradient entry lies in (A1,A2), making that ideal the unit ideal. Multiplicities are one. Their geometric irreducibility follows by pulling back through the automorphism G.

For c!=0, a factorization of E_sigma-c would pull back to one of sigma-c, impossible. Consequently, if G is NOT invertible, EVERY constant-level polynomial

    E_sigma-c

is geometrically irreducible, including c=0. This assertion concerns constant levels on the WHOLE intermediate affine space over algebraic extensions of k; it does not assert geometric irreducibility of every further special fibre after fixing all E_u.

Thus a proper factor must land in the all-irreducible-fibre branch of the component problem. Published irreducibility reductions, including the discussion in de Bondt--Yan, do not make that branch invertible. No general JC(2), higher-dimensional inverse theorem or full atomicity is inferred.

At any cut of a finite polynomial factor chain the same theorem applies to its total source prefix G. A noninvertible whole Keller stage has field degree greater than one, and these degrees multiply. Thus the generic intermediate unit lattice is the FULL Z while the prefix is invertible, drops to ZERO at the first noninvertible stage, and never reappears. There is no proper finite-index intermediate unit lattice and no later recovery of a lost unit. This is an all-length chain rule, not a bound on the number of later unit-free factors.

## 4. Another necessary ring constraint, with the field boundary kept explicit

The intermediate generic ring is actually a UFD over K0:

    O(V_eta)=k[Z_1,...,Z_n] localized at all
              nonzero polynomials in (E_u,E_sigma).    (4.1)

Localization preserves unique factorization. This is an ARITHMETIC statement over K0; scalar extension to K0_bar does not automatically preserve factoriality. Geometric factoriality is not proved.

For a proper factor of a NATIVE map, or another map satisfying the UARELPLANE1 curve witness, the retained surface constraints therefore sharpen to:

- A NON-PLANE smooth geometrically integral affine complete intersection, with free tangent bundle and the whole carried A1 component in the native case.
- CONSTANT geometric units, not merely unit rank at most1.
- The exact embedded field, both saturation identities and at most finitely many missing points of both etale arrows.
- Arithmetic factoriality over K0, with geometric divisor-class behaviour still unresolved.
- Both whole field degrees greater than one, and a whole intermediate signal component with no reducible constant fibre.

The unit, arithmetic-factorial and whole-component conclusions apply to the broader F as well. The NON-PLANE conclusion uses the embedded-curve witness and is not asserted for arbitrary B,C without it. These are necessary conditions. They do not construct such a proper intermediate, exclude every unit-free non-plane surface, or make an abstract surface model into a whole polynomial factor.

## 5. All-degree sharp controls and the remaining task

The UAINTERSURF1 family

    J_d(S,Y)=(S,Y*(1+r+...+r^(d-1))), d>=1,
    r composed J_d=r^d, det DJ_d=d*r^(d-1),

is etale and surjective on U in every degree. At S=0 it maps Y to dY; on S!=0 solve r^d=r_target. For d>1 its unit map has index d, so it violates the primitive-unit hypothesis and does not contradict UAPRIMUNIT1. Whole root descent excludes that index in actual whole Keller factorizations. Its canonical whole lift still has the earlier source pole and wrong volume.

There is a second exact model for the node obstruction. The invariants of (2.2) are

    K[U,V,r,r^(-1)]/(U*V-(1-r)^d),
    U=S^d, V=Y^d.                                     (5.1)

Indeed invariant monomials satisfy i-j divisible by d, and SY=1-r. For d>1 this quotient is singular at U=V=0,r=1, exactly where the cyclic group fixes the source node. The quotient is etale away from that point; treating only that open part would miss the obstruction. This is an explicit calculation, not an imported quotient-classification theorem.

The all-degree nonconstant-unit intermediate branch is therefore CLOSED for every one-block retained-signal Keller F. For native maps the remaining proper generic surface is unit-free and non-plane. The explicit nonlinear native polynomial integration, complete canonical Spec(A), whole-source landing, effective gluing, absolute stable atomicity and improved costs remain open. Multi-block canonical receivers require their own two-dimensional source-fibre hypotheses; this theorem does not collapse their free dimensions.

Square-free preservation/root closure, curve Riemann-Hurwitz, genus-zero classification and normal Hartogs extension are classical published inputs. The primitive-unit and whole-factor deductions are written proofs with selected exact controls; independent proof review and publication priority remain unverified.
