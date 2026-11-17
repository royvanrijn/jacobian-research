# Plane charts, carried curves and the Danielewski gluing barrier

15 November 2026. Two written general claims with selected exact controls. Independent whole-proof replay, formal verification, external review and publication priority remain unverified. OP-KDESC and JC(2) remain open.

This extends the [log-volume/Picard obstruction](../universal-keller-log-picard/PROOF.md). Passing exact volume is insufficient even at the generic-surface level: a higher-power Danielewski target would give compatible JC(2) counterexamples on its classical affine-plane atlas. The carried native affine line must cross the boundary of EVERY plane chart.

Work in characteristic zero, geometrically unless otherwise stated. In an actual retained-signal whole factorization F=E composed G, UAINTERSURF1 gives

    U --g--> V --e--> A2, U=Spec K0[S,Y,(1-SY)^(-1)].

For native F, or any F with the UARELPLANE1 witness, there is a CLOSED curve Gamma=A1 in V on which e is a closed embedding, namely the carried whole component over the target line. Its entire image, not a punctured piece, is the relevant datum.

## 1. UAPLANECHART1: the carried curve cannot fit inside a plane chart

Let V be a smooth geometrically integral surface, e:V->A2 etale and separated, and j:A2->V an open immersion. Put H=e composed j.

Then H is a polynomial plane Keller map. Indeed it is etale, and its Jacobian is an everywhere invertible polynomial, hence a nonzero constant. If H is an automorphism, define

    s=j composed H^(-1):A2->V.

This is a section of e. A section of a separated unramified morphism is open and closed (Stacks tag024T). Since V is connected, its image is all of V. Thus

    H an automorphism implies e an isomorphism and j surjective. (1.1)

This is an elementary section argument, not a new general inverse theorem.

Now suppose Gamma:A1->V is a closed embedding, e composed Gamma is a closed embedding, and its ENTIRE image is contained in j(A2). Its inverse chart parametrization gamma=j^(-1) composed Gamma is polynomial. UARELPLANE1 applies to H composed gamma: published affine-line rectification and Gwozdziewicz line injectivity make H an automorphism. Equation (1.1) makes e an isomorphism.

For an actual whole native factor, a generic isomorphism e makes the whole E birational, hence invertible by the birational-etale affine self-map criterion in UAALLATOM1. Therefore a proper intermediate necessarily satisfies

    Gamma intersects V \ j(A2) for EVERY open plane chart j.    (1.2)

It is not enough to check a selected chart. Conversely, merely crossing all selected chart boundaries is not proved sufficient for a proper factor.

### The precise JC(2) boundary

If e is not an isomorphism, H is a NONAUTOMORPHIC plane Keller map on EVERY open plane chart. Thus any such e with even one plane chart would disprove JC(2) over the geometric retained field. Over characteristic zero this also gives a complex counterexample: the polynomial coefficients descend to a finitely generated field, embed in C, and nonbirational generic degree is preserved. A birational plane Keller map is already invertible.

In particular, assuming JC(2), any smooth integral non-plane surface containing an open A2 admits no etale map to A2. This is CONDITIONAL. We neither assume JC(2) in (1.1)-(1.2), nor infer it from the surface constraints.

The open-plane reduction is established literature: Dubouloz's Basel2009 talk explicitly uses this setting and explains why elementary differential obstructions do not decide x^2*z=y^2-1. Our compiler application uses the additional carried WHOLE affine line to obtain (1.2).

### Formal plane neighborhoods do not supply a global plane chart

Rectify the embedded target e(Gamma) to the line c=0 and identify Gamma with that line. For every n>=0 let

    L_n=Spec k[b,c]/(c^(n+1)).

The initial section L_0->V lifts UNIQUELY over e to s_n:L_n->V, by classical formal etaleness (Stacks tag025H). The lifts are compatible by uniqueness. Each s_n identifies L_n with the open-and-closed component supported on Gamma in V times_A2 L_n. Equivalently, the formal completion of V along Gamma is the standard formal plane along c=0. Since Gamma is a whole open-and-closed component of e^(-1)(c=0), the ideal (c) agrees near Gamma with its Cartier ideal, so this really is completion along that curve, not along the other components.

For a noninvertible native target, (1.2) nevertheless forbids any GLOBAL open A2 chart containing Gamma. Thus EVERY finite-order neighborhood, and the full formal neighborhood, can be standard while the required global chart is impossible. This is a classical formal-etale application, extending the UAETALEGAUGE1 formal-versus-global boundary to the actual intermediate carried curve. No independent new formal-etale theorem is claimed.

A control is e:(x,y)->(x,y+y^2) on {1+2y!=0}, with Gamma={y=0}. The quadratic branch y=(-1+sqrt(1+4c))/2 supplies every infinitesimal lift and its formal limit, but no polynomial inverse. The source has a nonconstant unit and misses the discriminant line in its target, so it is NOT an actual whole-factor intermediate. It only displays the local/global distinction.

## 2. UADANCHART1: the all-exponent atlas and its exact Keller compatibility problem

Let P(z) be square-free of degree d>=2 and m>=1. Over an algebraic closure normalize P to monic by rescaling v, and write its distinct roots alpha_1,...,alpha_d. Let

    W_m: u^m*v=P(z),
    D_i={u=0,z=alpha_i},
    U_i=W_m \ union_(j!=i) D_j.

The following CLASSICAL atlas is reproduced explicitly, not claimed as a new construction. For each i,

    j_i:A2_(u,w_i)->U_i,
    z=alpha_i+u^m*w_i,
    v=w_i*product_(j!=i)(alpha_i-alpha_j+u^m*w_i).       (2.1)

For u!=0 the inverse is w_i=(z-alpha_i)/u^m. Near D_i it is equivalently v/product_(j!=i)(z-alpha_j), whose denominator is nonzero there. These two descriptions cover U_i and agree; hence j_i is an isomorphism. The U_i cover W_m.

On the common open u!=0,

    w_j=w_i+(alpha_i-alpha_j)/u^m.                       (2.2)

The canonical volume Omega_m=du wedge dz/u^m pulls back on EVERY chart to du wedge dw_i. These are genuine A2 charts; the transition has a boundary pole. The torsor/multiple-origin atlas is explicitly present in Dubouloz2004, introduction and Section2, not new here.

### An exact reformulation of the target-map problem

An etale e=(B,C):W_m->A2 is equivalent to polynomial pairs

    H_i=(B_i(u,w_i),C_i(u,w_i))

such that every H_i has a nonzero constant Jacobian, and the identities

    H_i(u,w)=H_j(u,w+(alpha_i-alpha_j)/u^m)              (2.3)

hold in k[u,u^(-1),w] for all i,j. Necessity is restriction; sufficiency is Zariski gluing and locality of etaleness. The chart transitions have Jacobian one, so the constants agree.

UADANVOL1 gives constant units and a nonzero geometric Picard group for W_m, so W_m is not A2. UAPLANECHART1 therefore makes EVERY H_i a JC(2) counterexample if such a compatible family exists. One cannot glue plane AUTOMORPHISMS into the desired non-plane etale target. A regular primitive for Omega_m only solves a sum-of-differentials equation; it gives none of the pairs or compatibility (2.3).

For m=1 UADANVOL1 already excludes e by the volume period. For m>=2 the atlas identifies a genuine JC(2)-level problem, not an etale construction. Neither a cofinite g:U->W_m nor a different whole smooth spread is constructed.

### A terminating all-degree obstruction along the A1 ruling

Here is an elementary plane Keller reduction, used without a coefficient bound. If H=(f,g) in k[u,w]^2 has constant nonzero Jacobian, put p=deg_w f, q=deg_w g. Both cannot be zero.

If p=0, f=f(u) and f'(u)*g_w is a nonzero constant. Both factors are constants, so H is triangular and invertible. The case q=0 is symmetric.

If p,q>0 and q=l*p, let A(u),B(u) be the respective leading w-coefficients. The coefficient of w^(p+q-1) in Jac(H) is zero:

    q*A'(u)*B(u)-p*A(u)*B'(u)=0.

Thus the rational function B/A^l has zero derivative; characteristic zero makes it a scalar lambda. The target shear g->g-lambda*f^l strictly decreases q and preserves the Jacobian.

Swap the outputs when necessary and repeat whenever the larger positive w-degree is divisible by the smaller. The sum p+q strictly decreases, so this process terminates. If either degree is at most1, it reaches degree0 and proves H invertible.

For (2.3), every chart has the SAME w-degrees and leading coefficients, since (2.2) is a translation in w. Consequently each reducing target shear works simultaneously on the ENTIRE atlas. A hypothetical e therefore must terminate at

    p,q>=2, p does not divide q, q does not divide p.    (2.4)

Any target-coordinate presentation with a ruling degree at most1 is excluded. This is a necessary obstruction, not a classification of the nondividing sector or a proof of JC(2). The argument uses standard leading-coefficient calculus; no novelty is asserted for that elementary plane-map lemma.

### The carried curve must meet at least two root components

Suppose an actual native intermediate is geometrically W_m. Write the carried curve as

    gamma(s)=(U(s),V(s),Z(s)),
    U(s)^m*V(s)=P(Z(s)).

By UAPLANECHART1 its image cannot fit in any U_i. If the curve meets at most one D_i, it IS contained in the corresponding U_i (or in every U_i if it meets none), a contradiction. Thus it must meet at least TWO distinct root components.

In particular U is nonconstant and has at least two distinct zeros with different root labels Z(s)=alpha_i. A curve in u=0 lies in one component; a nonzero constant U meets none. Both are excluded.

Let N=deg U and r be its number of distinct zeros. At a zero s_0 of U of multiplicity e_0 with Z(s_0)=alpha_i, square-freeness gives

    ord_(s_0)(Z-alpha_i)=m*e_0+ord_(s_0)V>=m*e_0.

The labels include two distinct roots, so Z is nonconstant. Summing the resulting derivative orders gives the UNBOUNDED contact constraint

    N>=2,
    deg Z-1>=m*N-r,
    deg Z>=m*N-r+1>=(m-1)*N+1.                         (2.5)

For m>=2 this includes deg Z>=2m-1. These are intrinsic contact conditions for the chosen ruling and the ENTIRE carried curve. They do not assert a finite bound on possible curves or embedding degrees.

### Sharp contact alone does not produce the embedded curve

For EVERY m>=2, take P=z(z-1) and define

    c_m=(2m-1)!/((m-1)!)^2,
    Z_m(s)=c_m*integral_0^s t^(m-1)*(1-t)^(m-1) dt,
    U(s)=s(s-1),
    V_m(s)=Z_m(s)*(Z_m(s)-1)/U(s)^m.                   (2.6)

Z_m has order m at0, Z_m-1 has order m at1, so V_m is polynomial. These maps meet both root components and attain equality deg Z_m=2m-1 in (2.5).

Nevertheless they are NEVER geometrically injective. The symmetry Z_m(1-s)=1-Z_m(s) implies

    Z_m(s)-1/2=(s-1/2)*H_m(U(s))

for a polynomial H_m of degree m-1. Here H_m(0)=1 and

    H_m(-1/4)=Z_m'(1/2)=c_m*4^(1-m)!=0.

Thus H_m has a root lambda different from0 and -1/4. The two distinct solutions of s(s-1)=lambda have the same image

    (u,v,z)=(lambda,-1/(4*lambda^m),1/2).              (2.7)

The all-degree contact construction therefore fails the CLOSED-EMBEDDING requirement; it cannot be the carried native curve. For m=2, Z=3s^2-2s^3, V=4U-3 and lambda=1/2 give an exact quadratic collision. No degree sweep underlies the general argument.

This does not classify all embedded affine lines in W_m. Other curves meeting multiple components remain unclassified here.

## 3. Scope, literature and next gate

The prior UALOGVOL1/UAGEOPIC1 restrictions remain: a proper native intermediate is unit-free, negative-log and geometrically nonfactorial, with nonzero geometric Picard group and no retained-field invariant class. The present chart obstruction adds a condition on the WHOLE carried curve and exposes the JC(2) burden of the higher-power family.

The actual next gate is a compatible nonautomorphic Keller atlas (2.3), a carried CLOSED embedded line crossing at least two root components and satisfying (2.5), both cofinite etale arrows, exact saturated fields, and a DIFFERENT globally smooth retained-base spread. UADANVOL1 already excludes the raw spread a=z^d-u^m*v by its critical zero fibre. Passing any one test does not prove a proper whole factor.

The recent Borisov--Gabber--Vasiu2026 paper discusses sphere-like/generalized Danielewski surfaces and etale maps, largely in positive characteristic. We reviewed its Section2.2 discussion and retained source; its positive-characteristic constructions are NOT transplanted to this characteristic-zero problem. Dubouloz2009 already documents the precise difficulty of the exponent-two target; this is a priority boundary, not evidence of a new solution.

Finite nonlinear native integration, canonical complete Spec(A), whole-source landing, effective gluing, improved costs and absolute stable atomicity remain open. The multi-block canonical receiver requires its own two-dimensional residual fibre and curve witness before these single-block deductions apply. The plane-chart reduction, classical atlas and elementary degree lemma overlap established mathematics; compiler applications and the new curve controls require independent review and a fuller novelty audit.
