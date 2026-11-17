# Native nonlinear feedback: exact profile orbits and cubic descent certificates

15 November 2026. Three written general claims with selected exact controls. No independent whole-proof replay, formal verification, external review or novelty certificate. OP-KDESC remains open. This extends the [two-control normal form](../universal-keller-feedback-normal-form/PROOF.md) to an all-degree coordinate-orbit classification and an exact old-factor test.

Except for the general quasi-finite lemma, work over a characteristic-zero field k. Parameters u are retained; source and target coordinate operations involving them must fix u. Put

    t=1+x*y, q=t^2*z+y^2*(1+3*t), sigma=t*q, X=x*q,
    K=(sigma,B0,C0), B0=-y/2-3*x*q/2, C0=x*(5-3*t)-x^3*z.

A NATIVE map is K plus finitely many polynomial source coefficients times

    v_j=(0,-(j/2)*sigma^2*X^(j-2),-(j-2)*X^j), j>=4.

UANATIVE2GEN1 gives two unique effective controls at indices4 and5, without restricting their degrees. Every such correction has B-component in (x^2) and C-component in (x^4). Use UAM1 rational proof coordinates

    S=x/t, r=1/t, X=sigma*S,
    x=S/r, y=(1-r)/S,
    z=sigma*r^3-r*(r+3)*(1-r)^2/S^2.                       (0.1)

All statements concern actual polynomial maps on the WHOLE source. Rational and formal charts are used only to prove them.

## 1. UAQFSAT1: quasi-finite polynomial maps have polynomial rational descent

Let F:A^n_k->A^n_k be a dominant quasi-finite polynomial map, over ANY field. Then

    k[P] intersection k(F_1,...,F_n)=k[F_1,...,F_n].      (1.1)

In other words, a rational target function with polynomial source pullback is already polynomial on the target. Quasi-finiteness suffices; no flatness, etaleness, separability, finiteness or surjectivity assumption is needed. The earlier UAJOINTKER1 saturation argument is the etale special case.

To prove this, let h=A(Y)/B(Y) with coprime target polynomials A,B. If a nonconstant irreducible q divides B, dominance makes q(F) a nonconstant nonunit source polynomial. Choose an irreducible source divisor D above it. Its dimension is n-1. The restriction D->V(q) is quasi-finite, so the closure of its image also has dimension n-1 and is the whole irreducible V(q). Thus D dominates V(q). If h(F) were polynomial, D would divide A(F); density would then force q to divide A, contradiction. Therefore B is constant.

This is an elementary divisor/dimension consequence of the classical quasi-finite fibre-dimension facts, not a claim of a new general geometric theorem. The hypotheses matter. For F(x,y)=(x,xy), the rational target function Y2/Y1 pulls back to y, but F has a positive-dimensional fibre at the origin and is not quasi-finite. Constant affine units matter as well: the open immersion Gm->A1 allows1/Y to become regular. For every dominant quasi-finite affine-space self-map, each target prime divisor is therefore met generically, even when some higher-codimension fibres are missing.

## 2. UANATIVEORBIT1: exact all-degree normal form for profile-coordinate feedback

Let T=(sigma,Bout,Cout) be a polynomial native Keller map. Write K0=k(u,sigma). The following are EQUIVALENT:

1. Cout is affine in r over K0(S).
2. Some polynomial plane automorphism over K0 sends (Bout,Cout) to a pair affine in r.
3. There exist a UNIQUE retained-signal profile map H=(sigma,Bh,Ch), with higher profiles polynomial in u,sigma, and a UNIQUE h in k[u,a,c], such that

       T=E_f composed H,
       E_f(a,b,c)=(a,b+f(u,a,c),c),
       f(u,a,c)=a^5*c^3*h(u,a,c).                        (2.1)

4. After some identity stabilization, T is polynomial left-right equivalent to a retained-signal profile map, with the SOURCE coordinate automorphism normalizing the field k(u,sigma) and all coordinate operations fixing u.

Condition4 is a restricted stable equivalence, not unrestricted stable equivalence. In condition2 the plane-coordinate coefficients may be rational in u,sigma; the conclusion (2.1) nevertheless algebraizes them globally. No degree bound and no tameness assumption over k[u,a] is imposed.

Consequently all native feedback with an affine third output is a coordinate change of an existing absolutely stable atomic profile. Arbitrarily nonlinear r-dependence in the other output can be a coordinate artifact. A genuinely new map for this signal-respecting profile-orbit question must have a third output nonlinear in r, and cannot be made affine in r by a generic polynomial plane coordinate change. Neither property alone proves full-map atomicity.

### From a generic affine pair to a profile

Assume condition2. Rescale one target coordinate over K0 to make the automorphism Jacobian one. Let its affine pair be

    p(S)*r+b(S), q(S)*r+c(S).

It is regular on the generic sigma hyperbola-complement leaf. Its coefficients are Laurent polynomials in S. The Laurent Wronskian/source-regularity argument in UAAFFINER1, BEFORE imposing native jets, gives

    p'*q-p*q'=-1/S,
    p,q in span{S,S^(-1)}.

Indeed Riemann-Hurwitz first gives powers S^n,S^(-n), and source regularity at r=1-S*y forces n=1. A determinant-one linear target change over K0 puts p=1/(2S),q=S. Integrating b'*S-c'/(2S)=0 and translating target constants gives a normalized pair

    Bh=(r-1)/(2S)+sum_(n>=1) beta_n*S^n,
    Ch=(r+1)*S+sum_(n>=1) (2n/(n+2))*beta_n*S^(n+2).      (2.2)

At x=0 both this pair and T restrict to (-y/2,0). The residual plane automorphism over K0 therefore fixes the line c=0 POINTWISE. Its c image is an irreducible coordinate divisible by c, hence gamma*c. Its b image is epsilon*b+f(c). The fixed line and determinant one give epsilon=gamma=1 and f(0)=0. Thus Cout=Ch and Bout=Bh+f(Ch), with f in K0[c].

The native C jet forces beta_1=-3*sigma/2. All other terms of Ch have the retained-profile form with indices>=4. Their rational sigma profiles are in fact polynomial. To see this using Ch ALONE, exclude every finite nonzero coefficient pole at sigma=alpha: on that irreducible source leaf X is nonconstant, so distinct powers X^j cannot cancel the leading residue. At sigma=0 use the t=0 component, where X=1/x is again nonconstant; every j-2 is nonzero in characteristic zero, so negative sigma powers cannot cancel either. Thus all profiles are polynomial over k(u). Comparing coefficients of x^j on y=0,z=sigma makes their u coefficients polynomial as well. The full H is therefore an actual retained-profile polynomial Keller map by UAPROFILE1, with the UAALLATOM1 degree and atomicity properties.

Now f(sigma,Ch)=Bout-Bh is polynomial. UAQFSAT1 for H with identities makes f a polynomial in u,a,c. The remaining native landing test determines exactly which f are allowed.

### The exact target-coordinate landing ideal

For any retained-profile H and ANY target polynomial automorphism E with E composed H native, dominance first gives E_a=a. On x=0, H covers the entire target plane c=0 and native corrections leave its values unchanged. Hence E fixes that plane pointwise. Coordinate irreducibility gives E_c=gamma*c. Since (gamma-1)*Ch would have source x-order1 while a native C correction has order at least4, gamma=1. It follows that E_b=b+f(u,a,c).

Using the two-control basis, an output correction (f(sigma,Ch),0) has effective controls

    delta A=-3*f(sigma,Ch)/(sigma^2*X^2),
    delta B= 2*f(sigma,Ch)/(sigma^2*X^3).                 (2.3)

They are polynomial exactly when f belongs to the principal target ideal

    (a^5*c^3) in k[u,a,c].                               (2.4)

For necessity, source x valuation of f(sigma,Ch) equals target c valuation: Ch=x*L with L modulo x equal2, and sigma on x=0 is freely variable. Source q valuation equals target a valuation, since t is a unit on q=0 and Ch there is nonconstant (on y=z=0 it is2x). The denominator sigma^2*X^3=t^2*x^3*q^5 therefore requires c^3 and a^5. Conversely if f=a^5*c^3*h and Ch=x*L, (2.3) becomes

    delta A=-3*t^3*x*q*L^3*h(u,sigma,Ch),
    delta B= 2*t^3*L^3*h(u,sigma,Ch),                     (2.5)

both whole-source polynomials. Uniqueness of H follows from Ch and the independent powers X^j in its profile expansion; f and h are then unique by dominance. This proves condition2 implies3, and3 gives condition1.

### One affine coordinate is already enough

Suppose a nonzero K0-linear combination U of Bout,Cout is affine in r, U=p*r+b. The native source jets show p!=0 and b'!=0: if the B coefficient is nonzero, p,b have opposite nonzero S^(-1) poles; if only C occurs, p=S times that coefficient plus higher terms, and b has the same nonzero linear term. Choose a complementary linear target coordinate V. If V had a negative lowest r power m, the Jacobian would have a nonzero unique term m*b'*V_m*r^(m-1), impossible. Hence V is polynomial in r.

If its highest degree N>=2 has coefficient v_N, the coefficient of r^N in the Jacobian gives

    N*p'*v_N-p*v_N'=0,
    v_N=constant_in_K0 * p^N.

Subtract that constant times U^N from V. This is a polynomial target shear over K0 and lowers N. Repeating finitely gives an affine pair. In particular condition1 implies2. This argument imposes no bound on the other output's degree and allows its original Laurent expression; the negative powers are ruled out by the exact Jacobian equation.

### Identity stabilization does not hide a different target orbit

Under condition4, UASIGFIELD1 makes the source map on P a weighted scaling; its remaining identity coordinates are a vertical polynomial automorphism. For a profile H0 this scaling satisfies

    H0(eta*x,eta^(-1)*y,eta^(-2)*z)=L_eta composed H1(P),
    L_eta=diag(eta^(-2),eta^(-1),eta),
    a_j^H1(sigma)=eta^(-j-1)*a_j^H0(eta^(-2)*sigma).

The compiler outputs of the target coordinate map cannot depend on the added identity coordinates: over k(u,P), the vertical source automorphism makes those coordinates freely variable, while T is independent of them. Dominance makes this independence a polynomial identity on the target. Thus T=E composed H1 for a polynomial target map E. The full target automorphism has block-triangular Jacobian, so E is Keller; its inverse supplies a polynomial right inverse for E by setting its target identity coordinates to zero. Generic degrees multiply to one, so E is birational and hence an automorphism by the affine birational-etale lemma used in UAALLATOM1. No JC(m) assertion is needed.

The preceding target-landing classification now proves (2.1). Conversely (2.1) gives condition4 without stabilization. The same argument works relatively with retained u; polynomial source scalings and their inverses force eta to be a unit of k[u], hence a field constant.

As a corollary, two retained-profile maps H_f,H_g are equivalent in this signal-respecting stable category exactly when, for some eta in k*,

    g_j(u,sigma)=eta^(-j-1)*f_j(u,eta^(-2)*sigma)

for every j. The weighted source formula gives sufficiency. For necessity, the unique shear in (2.1) must vanish because both maps have affine B outputs while every nonzero permitted shear has c-degree at least3. This is a RESTRICTED stable profile-moduli classification, not a classification under source coordinates moving the signal field.

### A finite degree-reduction gate, not an unbounded coefficient search

If both outputs are POLYNOMIAL in r with positive degrees m,n and m+n-1>1, their leading coefficients p_m,q_n satisfy

    n*p_m'*q_n-m*p_m*q_n'=0.

Whenever m divides n, q_n/p_m^(n/m) is in K0. A polynomial target shear lowers n; similarly when n divides m. These operations strictly lower the degree sum and therefore terminate. A target coordinate independent of r cannot occur: its restriction to the native image line c=0 would be constant, so coordinate irreducibility would make it a constant plus a scalar times Cout; the native C jet forbids r independence. Reaching degree1 invokes the preceding argument and proves (2.1).

A polynomial-r candidate surviving this elementary reduction must stop with two degrees at least2, neither dividing the other. This is a remaining gate, not proof that such a candidate exists or is atomic. In particular the WHOLE degree-at-most-two polynomial-r sector reduces: if both degrees are2, the top coefficients are K0-proportional, yielding an affine linear combination. Formula (2.1) then permits no nonzero shear of degree at most2, because every nonzero f has c-degree at least3. Such a map is IDENTICALLY a retained-profile map. Negative Laurent-r sectors remain separate.

The quadratic conclusion uses native MODULE membership, not only the source jets. For example the whole polynomial target shear b->b+a^5*c^2 gives a quadratic-r Keller map satisfying the source jets, but its effective quintic control has an x pole and is not polynomial. On y=0,z=1 that control is2*(2-x^2)^2/x. Thus the broader jet hypothesis used in the earlier affine-r theorem cannot replace native module membership in the present nonlinear classification.

### Genuine nonlinear-r examples, and their exact marked contact

For every m>=3, retained profile H and ell>=1, take

    f=lambda*a^5*c^m*I(a)^ell,
    I(a)=a*(3*a-1)*(4*a+1).

Equations (2.5) give explicit polynomial controls with r-degree m, exact determinant one and generic degree equal to H. They retain its entire restricted map over V(I(a)^ell), including the complete filter fibres when H is the cubic. These are positive whole-source nonlinear-r solutions, but literal target-coordinate changes: they produce no new stable class and preserve every composition factor of a full receiver.

More precisely, inside this target orbit, equality with H over the entire thickening V(I(a)^ell) holds exactly when I(a)^ell divides f. The exact ideal is

    (a^5*c^3) intersection (I(a)^ell)
      =(a^max(5,ell)*(3*a-1)^ell*(4*a+1)^ell*c^3).         (2.6)

The equivalence follows by valuations on nonzero signal fibres and on both components t=0,q=0 at zero; Ch is nonconstant on each relevant component. Thus arbitrarily high marked contact and nonlinear controls alone do not certify a new map class.

## 3. UACUBICDESC1: exact all-degree old-factor certificates in one cubic algebra

Fix the original cubic K, not an arbitrary unknown profile. Let F0=k(u,a,b,c). Its full source function field is the degree-three field

    L=F0[S]/(E), E=a*S^3+b*S^2+S-c/2,
    r=E'=1+2*b*S+3*a*S^2,
    x=S/r, y=(1-r)/S,
    z=a*r^3-r*(r+3)*(1-r)^2/S^2.                         (3.1)

UAM1 supplies the reconstruction and generic degree. E is irreducible and separable. S and r are invertible in this generic field; use Euclid modulo E, with S^(-1)=2*(1+b*S+a*S^2)/c. Every source polynomial therefore has a UNIQUE reduction

    p0(u,a,b,c)+p1(u,a,b,c)*S+p2(u,a,b,c)*S^2.            (3.2)

It is a polynomial in the old outputs exactly when p1=p2=0. The forward direction is clear; the converse is UAQFSAT1 for K with identities. Thus the test requires no target-degree ansatz or coefficient cutoff. Polynomial ring division cancels the apparent target denominators after descent is certified. The cubic algebra has fixed dimension three, although coefficient sizes may grow with the input.

For any polynomial candidate T, applying (3.2) to every new output detects EXACTLY whether an old-output polynomial factor T=E_target composed K exists. In the retained-signal native class only the last two outputs need checking. If T is Keller, dominance and the chain rule make E_target Keller. Generic degrees satisfy

    deg T=3*deg E_target.

If deg E_target>1 this is an actual nontrivial polynomial factorization. If deg E_target=1, the target map is an automorphism and the native landing theorem gives exactly

    T=(sigma,B0+sigma^5*C0^3*h(u,sigma,C0),C0).           (3.3)

Consequently among output-descended native Keller maps, self-atomicity is equivalent to degree3 and to the cubic-coordinate form (3.3). Degree1 etale affine self-maps are automorphisms; no JC(2) inverse assertion is used. Failure of this old-factor test does NOT prove atomicity: another intermediate map, another retained profile or a source change moving the signal field may remain.

The scalar formal tangent R_4*v_4 is not output-descended, as witnessed by its unequal values at the two equal-K points in UAETALEGAUGE1. For any hypothetical finite polynomial parameter integration with that tangent, reduction (3.2) has coefficients polynomial in the parameter over F0. Its nonzero first derivative therefore proves generic failure of descent through the old K. Such an integration would remove that old factor and escape the classified signal-respecting profile orbit. Every map in that orbit has exactly the fixed r coefficient S in its C output, so a parameter tangent there has no r dependence in that output; the displayed tangent does. It would still need a separate proof that NO other factors remain. No such polynomial integration is constructed here.

## 4. Literature and the remaining construction

The quasi-finite/divisor lemma extends the retained etale saturation argument without importing a general algebraization theorem. The generic coordinate reduction uses only explicit elementary shears, the UAAFFINER1 Laurent Wronskian argument and elementary irreducibility of a plane coordinate. We do not assume an arbitrary Keller plane map invertible, and do not identify the hyperbola source leaf with affine plane.

The tame/symplectic approximation paper [Kanel-Belov et al., 2017](https://arxiv.org/html/1707.06450v1) separates power-series approximation from the unresolved polynomiality of a limit. Its introduction and main approximation statements were reviewed as related context, not as a theorem dependency solving our nonlinear determinant equation. The exact nonlinear coordinate examples above illustrate why polynomial success must also be checked against coordinate orbits and actual intermediate fields.

The remaining native candidate must have Cout genuinely nonlinear in r, survive the polynomial-r degree-reduction gate or use an unclassified Laurent sector, and fail all actual-factor tests relevant to its proposed atomicity. It must also keep the full canonical marked scheme, whole-source landing and effective gluing. Changed first outputs and genuinely changing multi-block settings remain other departures. Canonical absolute atomicity, improved costs and independent proof/novelty review remain open. No ever-deeper coefficient, source-degree or formal-order sweep is scheduled.
