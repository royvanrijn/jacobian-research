# Embedded curves and every intermediate generic surface

15 November 2026. Three written general claims with selected exact controls. No independent whole-proof replay, formal verification, external review or novelty certificate. OP-KDESC remains open. This strengthens the [plane-completion theorem](../universal-keller-plane-completion/PROOF.md): the former boundary-pole alternative for signal-retaining factors is now closed without a Laurent-degree restriction.

Work over a characteristic-zero field k. Identity variables u may have any finite length. Put

    t=1+xy, q=t^2*z+y^2*(1+3t), sigma=t*q, X=x*q,
    K=(sigma,B0,C0),
    B0=-y/2-3*x*q/2, C0=x*(5-3t)-x^3*z.

A native map T=(u,sigma,B,C) adds finitely many polynomial source multiples of

    v_j=(0,-(j/2)*sigma^2*X^(j-2),-(j-2)*X^j), j>=4

to the last two outputs. UANATIVE2GEN1 gives two unique effective controls. The results below require only the displayed retained coordinates, the Keller condition, and the embedded-curve witness; arbitrary native degrees and positive or negative Laurent-r dependence are allowed.

## 1. UARELPLANE1: an embedded curve cancels a plane Keller postfactor

The published inputs are HC4MYGL0, Gwozdziewicz's [Injectivity on one line](https://arxiv.org/html/alg-geom/9305008v1), Theorem1.1, and Abhyankar--Moh rectification of a polynomially embedded affine line, explicitly invoked in that source. They were already imported in the HC(4) research. We reproduce their elementary consequence, not a proof of these classical theorems.

**Plane cancellation lemma.** Let E:A2->A2 be Keller. If a polynomial curve gamma:A1->A2 makes E composed gamma a CLOSED EMBEDDING, then E is an automorphism.

Indeed, with s the curve parameter,

    k[s]=k[(E composed gamma)_1,(E composed gamma)_2]
          subset k[gamma_1,gamma_2] subset k[s].

Thus gamma is itself a closed embedding. Over an algebraic closure, rectify gamma(A1) by a polynomial plane automorphism. The resulting plane Keller map is injective on a straight line, hence invertible by HC4MYGL0. Descent gives the unique inverse over k.

In particular gamma need not be a straight line, and a map into gamma may have poles AWAY from that curve. The proof uses only its polynomial restriction to the curve. This removes the extension assumption in the earlier boundary-pole discussion.

### The general relative theorem

Let F:A^n->A^n be Keller, n>=2. Split its TARGET coordinates as

    F=(p_1,...,p_(n-2), f,g), p=(p_1,...,p_(n-2)),
    K0=k(p_1,...,p_(n-2)).

For n=2 the tuple p is empty and K0=k. Suppose the generic source fibre X_eta of p contains a morphism i:A1_K0->X_eta such that (f,g) composed i is a closed embedding in A2_K0.

Take ANY whole polynomial self-map factorization

    F=E composed G.

Both factors are Keller by dominance and the constant chain-rule product. On the INTERMEDIATE affine n-space define

    q=(E_1,...,E_(n-2)), V_eta=q^(-1)(generic point).

Then

    E is an automorphism
       iff the GEOMETRIC generic fibre V_eta is A2.     (1.1)

For the forward direction, q is projection after the automorphism E. For the reverse, base-change E to the generic point of the target p-space. Its restriction E_eta:V_eta->A2_K0 is etale, being a base change. After an algebraic closure and an identification V_eta with A2, it is a polynomial plane Keller map: its Jacobian is a unit in a polynomial plane ring, hence constant. The curve G composed i maps into V_eta and its E_eta image is the given closed embedding. Plane cancellation makes E_eta invertible after this scalar extension. Isomorphism descends faithfully flatly, so E is birational, and the classical birational-etale affine self-map criterion used in UAALLATOM1 makes E an automorphism.

No coordinate system on V_eta over K0, global fibre trivialization, or affine-plane-form classification is assumed. Geometric generic A2 is enough. The result applies after field extension and identity stabilization whenever the retained tuple leaves exactly TWO generic coordinates and the stated curve witness exists. It is not an unrestricted affine-space inverse theorem in higher residual dimension.

For example the original three-dimensional K is noninvertible despite mapping the curve (0,s,0) to the embedded parabola (4s^2,-s/2,0). The plane dimension in the cancellation lemma is essential.

## 2. UASIGFACT1: all signal-retaining native factorizations are trivial on the target side

For a native T, on x=0 one has

    (u,sigma,B,C)=(u,z+4*y^2,-y/2,0).

Fixing generic retained values (u,a), the ACTUAL source curve

    i(s)=(u,0,s,a-4*s^2)                               (2.1)

lies in sigma=a, and its last outputs are (-s/2,0). It is a closed embedded target line, with s=-2B. This curve is in the original generic source, away from the deleted hyperbola; it does not depend on a plane completion.

Consequently, in every whole polynomial factorization

    T=E composed G,
    G=(u,sigma,Gb,Gc), E=(u,a,Eb,Ec),                    (2.2)

the target factor E is an automorphism. Its generic retained fibre is literally A2, so UARELPLANE1 applies. There is NO need to assume nonnegative r powers of the final OR intermediate outputs. This closes the former UAPLANEATOM1 boundary-pole alternative in the signal-retaining class, in ALL degrees.

The same conclusion holds if the retained target tuple (u,a) is merely a PART OF A POLYNOMIAL COORDINATE SYSTEM on the intermediate affine space, after a field extension if necessary. Change intermediate coordinates to complete that tuple, then apply (1.1). More generally, (1.1) closes the case where that tuple has geometric generic fibre A2 even without a global coordinate completion.

This is absolute and identity-stable relative to the specified retained tuple; it does not assert full absolute stable atomicity of every native T. The extra identity coordinates are included in that tuple, while the intermediate factors in the next section may mix them arbitrarily. A nontrivial factor must express the retained tuple by functions with a genuinely non-plane geometric generic fibre.

The argument applies more broadly to any whole polynomial Keller map with first outputs (u,sigma) and this embedded-curve witness. Native module membership is a sufficient way to obtain it, not a hidden necessity. Mutually changing blocks with more than two remaining generic coordinates require their own analysis.

## 3. UAINTERSURF1: every stabilized factor passes through a constrained non-plane surface

Allow ANY number of identity variables u and ANY polynomial self-map factorization T=E composed G of a native T. The factors need NOT retain u or sigma, and may mix all coordinates. Let n=length(u)+3 and set

    p=(u,sigma), q=(E_u,E_sigma), K0=k(u,a),
    U_eta=Spec R, R=K0[S,Y,r^(-1)], r=1-SY.             (3.1)

UAM1 and UARPLANE1 identify this generic source leaf; the coordinate change is S=x/t,Y=y,r=1/t. Set

    V_eta=Spec B=q^(-1)(generic retained values),
    A=K0[b,c], L=Frac B.

There is an actual diagram of etale dominant maps

    U_eta --G_eta--> V_eta --E_eta--> A2_K0,             (3.2)

with composition the generic T pair. Both maps are base changes of the whole Keller factors. The following are necessary conditions on EVERY intermediate model, including after arbitrary stabilization mixing.

### Geometric integrality, complete intersection and tangent bundle

The map q is smooth of relative dimension2: it is target projection composed with the etale E. Thus V_eta is a smooth affine surface, cut out by the length(u)+1 equations q_i=generic_value_i in intermediate affine n-space. These equations form a complete intersection; smoothness supplies independent differentials.

Dominance gives B injective into R. After an algebraic closure of K0, flat scalar extension preserves this injection into

    R_bar=K0_bar[S,Y,(1-SY)^(-1)],

which is a domain. Hence V_eta is GEOMETRICALLY INTEGRAL, and G_eta stays dominant. It is unirational because it is dominated by the rational open plane U_eta. No rational-surface classification is imported.

Since E_eta is etale, its differential identifies the tangent bundle of V_eta with the pullback of the free tangent bundle of A2. In particular

    T_(V_eta) is free of rank2, omega_(V_eta) is trivial. (3.3)

These necessary conditions do not identify V_eta with affine plane.

### The carried curve is a whole component of one target-line preimage

The curve gamma=G_eta composed i has E_eta image (-s/2,0). Therefore B->K0[s] is surjective: its image contains s=-2*E_b. This proves gamma is a CLOSED EMBEDDING even though G_eta need not extend across r=0.

In E_eta^(-1)({c=0}), gamma gives a section over the whole affine b-line. The preimage is etale and separated over that line; its section is open and closed. Thus it contains a whole component isomorphic to A1 on which E_eta is an isomorphism. Other components have not been excluded. The standard section fact is [Stacks tag024T](https://stacks.math.columbia.edu/tag/024T).

### Units and saturated polynomial rings

The UFD localization R has

    R^*=K0^* times r^Z.

Dominant pullback therefore injects B^*/K0^* into Z. Its unit rank is at most1; the same assertion holds geometrically. No assumption that all intermediate units are constant is justified.

The full polynomial saturation UAQFSAT1 localizes to BOTH arrows of (3.2):

    B=R intersection L,
    A=B intersection K0(b,c).                           (3.4)

For the first identity, write an element of R as H(P)/h(p), with h a nonzero polynomial in retained values. If it is a rational function f(G) on the intermediate space, h(q)*f pulls back to a polynomial H. Saturation for the WHOLE G makes h(q)*f a polynomial in intermediate coordinates, so f is in B. The converse is immediate. Apply the same argument to E, localizing target denominators only in the retained tuple, for the second identity.

In fact coprime TARGET polynomials in A remain without a common source divisor in B. This follows either by localizing the full no-contracted-divisor argument, or by etaleness and divisor-image dimension. A nonconstant target polynomial cannot pull back to a unit in B: if it did, its reciprocal would be in the intersection in (3.4), contradicting polynomiality on A2. Hence E_eta meets EVERY target prime divisor generically, including after algebraic scalar extension. It is not merely an arbitrary open etale map to the plane.

Here is a precise geometric-base-change argument. The image of each WHOLE Keller factor is open and its missing closed set has codimension at least2, by the divisor argument of UAQFSAT1. Such a set has dimension at most n-2. Taking its generic fibre over the retained base of dimension n-2 leaves dimension at most0. Hence BOTH arrows of (3.2) have open image with at most finitely many missing target points. Image formation commutes with field extension by the nonemptiness criterion for fibres, and a zero-dimensional missing set stays zero-dimensional. Thus both arrows meet every geometric target divisor as well. Normality and the same divisor-valuation argument give the geometric versions of (3.4). No surjectivity at those possible missing points or finiteness is inferred.

Thus, for any actual factor, its intermediate coordinate ring is determined by its embedded function field L as R intersection L. An arbitrary abstract surface or intermediate field is not automatically such a whole polynomial factor.

### The exact plane obstruction for unrestricted factors

By UARELPLANE1,

    E is an automorphism iff V_eta_bar is A2.           (3.5)

A proper factorization with neither stage invertible consequently requires a NON-PLANE geometric generic intermediate surface satisfying all preceding constraints. Its field degrees must also make both whole factors nonbirational; a degree-one whole Keller factor is an automorphism. This is an all-factor necessary obstruction after arbitrary identity stabilization, not a proof that no such surface exists or that every model satisfying it algebraizes. The endpoint V_eta=U_eta with G identity satisfies the listed conditions; a PROPER intermediate must make both field degrees greater than one.

In particular intermediate boundary poles alone can never evade the obstruction when the intermediate retained fibre is a plane. The former pole route is closed; the open geometry is the non-plane generic fibre of q, not merely a choice of coordinates on the old completed plane.

### Sharp generic controls: the source surface itself has many etale endomorphisms

The abstract surface U_eta cannot simply be declared an atom. For EVERY integer d>=1 put

    H_d(r)=1+r+...+r^(d-1),
    J_d(S,Y)=(S,Y*H_d(r)).
    r composed J_d=r^d,
    det DJ_d=d*r^(d-1).                                (3.6)

This is an etale endomorphism of U_eta of generic degree d. On S!=0 use coordinates (S,r): it is (S,r)->(S,r^d). At S=0 one has r=1 and Y->dY, so d>1 is nonfinite: a generic degree-d fibre loses sheets there. It still maps the distinguished affine line isomorphically. This is an all-degree positive control, not a whole affine-space Keller construction.

Its canonical lift back to the SAME original sigma leaf has

    x'=x*t^(d-1),
    y'=y*(1+t^(-1)+...+t^(-(d-1))).

For d>1 the second expression has an uncancelled t=0 pole. The canonical source volume r*dS wedge dY pulls back with ratio

    d*r^(2d-2),                                       (3.7)

so this lift is not a constant-Jacobian whole map either. UASIGFIELD1 independently excludes nontrivial whole Keller lifts normalizing that original signal field. These examples do not rule out other non-plane intermediate models with their required volumes and whole-source realization.

Likewise the open immersion U_eta->A2 is etale and contains the embedded line but FAILS (3.4): 1/(1-bc) is regular on the source and rational but nonpolynomial on the target. It cannot serve as the generic target stage of a whole affine-space Keller factor. These boundaries explain why surface shape, a generic etale formula and an embedded curve are not sufficient construction data.

## 4. Current gap and verification boundary

The all-degree relative factor problem is now closed for the entire native class, including negative Laurent-r outputs. Every arbitrary stabilized factorization has the surface diagram (3.2); a proper intermediate must be a non-plane smooth geometrically integral affine complete intersection with free tangent bundle, the carried whole A1 component, the unit bound, at most finitely many missing target points, and BOTH saturation identities.

This does not construct the nonlinear native integration or the atomic canonical receiver. The fixed canonical Spec(A), whole-source polynomial landing at every signal, effective gluing, absolute stable atomicity and improved costs remain open. The canonical trace receiver has several changing blocks; this single-block witness must not be applied to it without the two-dimensional fibre and curve hypotheses.

The plane line theorem, embedded-line rectification and standard etale base-change/section facts are published inputs. The compiler deductions and surface constraints are written general arguments; selected exact controls do not constitute independent whole-proof replay, formal verification, external review or a publication-priority certificate.
