# Scalar self-feedback: whole-source rigidity and a formal integration barrier

15 November 2026. Three written general claims with selected exact controls. No independent whole-proof replay, formal verification, external review or novelty certificate. OP-KDESC remains open. This addresses the additional source-channel gate after the [root-feedback theorem](../universal-keller-root-feedback/PROOF.md), for one affine coefficient direction, without restricting the control to a root-polynomial subring.

Work in characteristic zero. Identity parameters u are allowed. Retain the actual entire source triple P=(x,y,z) and

    t=1+x*y, q=t^2*z+y^2*(1+3*t), sigma=t*q, X=x*q.

For finite polynomial lists a_j^0(u,sigma),c_j(u,sigma), j>=4, with some c_j nonzero, take an ARBITRARY polynomial g(u,x,y,z), and set

    a_j(P,u)=a_j^0(u,sigma)+c_j(u,sigma)*g(P,u).             (0.1)

Thus the whole coefficient vector moves along one specified affine line. Several entries can vary, but all use the same scalar g AFFINELY. The source map is the retained polynomial compiler

    Phi=(sigma,
      -y/2-3*x*q/2-sum_j (j/2)*a_j*t^2*x^(j-2)*q^j,
      x*(5-3*t)-x^3*z-sum_j (j-2)*a_j*(x*q)^j),             (0.2)

with u kept unchanged. This is globally polynomial for every g. The theorem decides exactly when it is Keller.

## 1. The exact nonlinear response equation

Use the UAM1/UAROOTDAG1 proof coordinates

    S=x/t, r=1/t,
    x=S/r, y=(1-r)/S,
    z=sigma*r^3-r*(r+3)*(1-r)^2/S^2.                        (1.1)

Put K=k(u,sigma). The entire polynomial source ring embeds in the Laurent ring K[S,S^-1,r,r^-1]. No source locus is removed in the construction; this embedding is only used to prove identities and exclude poles.

In K[S] define

    h0(S)=S+sigma*S^3+sum_j a_j^0*sigma^j*S^j,
    V(S)=sum_j c_j*sigma^j*S^j,
    A=S*h0''-h0', B=S*V''-V', L=S*V'/V.                    (1.2)

V is nonzero and has order e>=4 at S=0. At fixed sigma,u, the outputs are

    Bout=(r-h0'-g*V')/(2*S),
    Cout=r*S+2*h0-S*h0'+g*(2*V-S*V').                      (1.3)

When g is an independent frozen parameter the two-dimensional source Jacobian is -r/S, as for the cubic. After substituting the actual source function g, the rank-one determinant lemma gives the EXACT identity

    det DPhi = 1 + V/(r*S) * E(g),
    E(g)=S*g_S+((1-L)*r+A+B*g)*g_r.                        (1.4)

The g term in the response is essential. A first integral of the response at g=0 need not solve this nonlinear equation.

Any constant nonzero determinant in (0.2) is one. Indeed at x=0 every parameter-response vector has value zero: its B component contains x^(j-2), its C component x^j. The additional rank-one dg term vanishes there, while the frozen-parameter determinant remains one by UAPROFILE1. Thus (0.2) is Keller exactly when E(g)=0.

## 2. UASELFLAUR1: complete Laurent solutions and their only exceptions

The following field statement isolates the mechanism. Let K be any characteristic-zero field, take A in K[S] with A(0)=-1, and a nonzero V in K[S] with order e>=3 and degree d. Use B,L as in (1.2). Every solution

    g in K[S,S^-1,r,r^-1], E(g)=0                         (2.1)

is either constant in K, or V is a MONOMIAL v*S^j, j>=3, and g has the explicit form

    a0=2*(j-1)/(j*(j-2)*v), A(S)=sum_(n>=0) A_n*S^n,
    g=a0*S^(1-j)*r
       -a0*sum_(n>=0) A_n/(n+j-1)*S^(n+1-j)
       +b0*S^(2-2*j),                                   (2.2)

where b0 is arbitrary in K. Conversely all (2.2) solve (2.1). These are rational chart solutions, not globally polynomial source controls.

### Laurent degree in r

B is nonzero since V has order at least three. If the lowest r exponent m of g were negative, the unique lowest term in E(g) would come from B*g*g_r, with exponent2m-1 and coefficient m*B*g_m^2. This is nonzero in characteristic zero, contradiction. If the highest r exponent n were at least two, B*g*g_r would similarly have unique highest exponent2n-1, coefficient n*B*g_n^2, again impossible. Thus

    g=a(S)*r+b(S), a,b in K[S,S^-1].                      (2.3)

If a=0, the equation reduces to S*b'=0 and g belongs to K. Otherwise the r coefficient and constant coefficient are

    S*a'+(1-L)*a+B*a^2=0,
    S*b'+a*(A+B*b)=0.                                   (2.4)

### One reciprocal and a degree/divisibility argument

Set w=1/a. The first equation in (2.4) becomes

    derivative_S(w*V/S)=V*B/S^2.

Let W be the polynomial primitive of V*B/S^2 with W(0)=0. Then

    a=V/(S*(C+W)), C in K.                              (2.5)

The order of W is2e-2 and its degree is2d-2. Both leading coefficients are nonzero: the contribution of v*S^j is j*(j-2)*v^2/(2j-2), with j>=3.

If C!=0, C+W has nonzero constant coefficient and degree2d-2>d, so it cannot divide V even in K[S,S^-1]. Hence (2.5) is not Laurent polynomial. If C=0, write V=S^e*P and W=S^(2e-2)*Q with P(0),Q(0) nonzero. Laurent polynomiality requires Q to divide P. But deg Q=2(d-e), deg P=d-e, impossible unless d=e. Therefore V is monomial. Substituting it into (2.5) gives a=a0*S^(1-j); integrating the second equation in (2.4) gives exactly (2.2). There is no logarithm: all denominators n+j-1 are nonzero. This proves the full field classification, with no degree bound or bounded coefficient search.

## 3. UASELFCTRL1: an entire polynomial control can only use the retained signal

For the polynomial control (0.1), the following are equivalent:

1. Phi is Keller.
2. det DPhi=1.
3. g belongs to k[u,sigma].                              (3.1)

In fact a stronger whole-map statement holds. Allow ANY rational source control g in k(u,x,y,z), and ask that the RATIONAL expression (0.2) extend polynomially to the whole original affine source. Such an extension is Keller exactly when

    g in k(u,sigma), and every effective
    a_j^0+c_j*g belongs to k[u,sigma].                    (3.2)

A rational control need not itself be polynomial when a common factor in the direction cancels it; for instance c_4=sigma,g=1/sigma gives effective coefficient1. The conclusion concerns the actual whole polynomial map and its effective profiles. This also excludes repairing a genuinely source-dependent rational scalar control merely by cancellation against the native response factors.

### Rational-control landing first forces determinant one

Let vB=-V'/(2*S), vC=2*V-S*V'. Whole polynomial landing makes g*vB and g*vC Laurent polynomials in S,r. Since vB is a nonzero function of S alone, g lies in K(S)[r,r^-1]. If the determinant is a nonzero constant c, (1.4) says

    E(g)=(c-1)*r*S/V.

The same lowest/highest r-degree argument as before gives g=a(S)*r+b(S), now with a,b rational in S. Its equations are

    S*a'+(1-L)*a+B*a^2=(c-1)*S/V,
    S*b'+a*(A+B*b)=0.                                   (3.3)

If a=0, the r coefficient gives c=1. Suppose c!=1 and write V=v*S^e+... at zero. Balancing orders in the first equation forces a=alpha*S^(1-e)+..., alpha!=0. On the ACTUAL plane y=0, r=1,S=x,sigma=z. Since vB has order e-2, regularity of g*vB there forces b=-alpha*S^(1-e)+... . The leading terms in (3.3) are

    -2*(e-1)*alpha+e*(e-2)*v*alpha^2=(c-1)/v,
    (e-2)*alpha*(1-e*v*alpha)=0.

Hence alpha=1/(e*v), and the first equality gives c=0, contradiction. Thus every entire Keller extension in this rational-control class has determinant one as well.

### Whole-map landing eliminates all nonconstant-r branches

With c=1, the reciprocal argument gives a=V/(S*(C+W)). Landing implies both a*V and a*V' are Laurent polynomials in S, by combining the r coefficients of g*vB and g*vC.

If C!=0, the denominator C+W is coprime to S. It would have to divide P^2, where V=S^e*P and P(0)!=0. But deg(C+W)=2d-2>2(d-e)=deg P^2. This is impossible. If C=0, write W=S^(2e-2)*Q. Landing of a*V implies Q divides P^2. Their degrees are equal, so Q is a nonzero constant multiple of P^2. Landing of a*V' then forces P to divide S*P', hence P is constant in characteristic zero. Thus V is monomial, and the remaining g is exactly the exceptional branch (2.2).

Its pole cannot be hidden by the native response. On y=0, the b0 term gives a B-output pole S^(-j), strictly worse than all others, so b0=0. The remaining B-output has residue

    -(j*v/2)*(j*a0/(j-1))=-j/(j-2) != 0                 (3.4)

at S=0. The polynomial background cannot cancel this pole. Therefore a=0, E(g)=S*b'=0, and g belongs to K=k(u,sigma).

UAPROFILE1 over k(u) now says that every effective profile is polynomial in sigma. To descend its parameter coefficients as well, restrict the C output to y=0,z=sigma and compare the independent powers x^j: sigma^j*a_j is in k[u,sigma]. Together with a_j in k(u)[sigma], this gives a_j in k[u,sigma]. The converse in (3.2) is the retained-profile determinant theorem.

For an initially polynomial g, the polynomial section x=y=0,z=sigma further gives g in k[u,sigma], proving (3.1). The effective direction must be nonzero; a specialization on which every c_j vanishes makes the control irrelevant and is not a classification of that degenerate control.

Therefore arbitrary polynomial dependence on t, q or ANY other own-source expression cannot supply a new scalar self-feedback construction along one affine coefficient line. Every allowed map is already in the retained profile family of UAALLATOM1: effective degree at least three, arithmetic/geometric symmetric monodromy and absolute stable self-atomicity in the existing scope. No fixed degree across specializations or common stable class for different profiles is asserted. This is a one-block statement with any other coordinates retained as identity parameters. It does not classify simultaneous mutually changing blocks, arbitrary curved coefficient families or several affine coefficient directions.

In particular multiplying an arbitrary self-source control by a marked factor I(sigma)^ell does not solve the Jacobian condition. It preserves the marked restricted morphism when the first output is retained, but an entire polynomial Keller map still requires the resulting effective scalar control to lie in k[u,sigma]. The marked-ideal identity of UAMARKEDDEF1 and Keller self-consistency are separate requirements.

## 4. UASELFFORMAL1: all-order formal Keller controls with no polynomial integration

The global rigidity is sharp. In the cubic background, choose one controlled coefficient j>=4. Let

    R_j=X^(j-1)*(q+sigma/(j-1))-3*X^(j+1)/(j+1),
    D_j=sigma^2*X^(2j-2), c_j=j*(j-2)/(2*(j-1)).          (4.1)

These are POLYNOMIAL source functions. In (S,r) coordinates,

    R_j=sigma^j*S^(j-1)*(r+1/(j-1)-3*sigma*S^2/(j+1)),
    D_j=sigma^(2j)*S^(2j-2).

R_j is a first integral of the response with g=0, but the actual self-control g=lambda*R_j has

    det DPhi - 1
      =lambda^2*j*(j-2)*R_j*sigma^3*X^(3j-3)*t,           (4.2)

which is not zero. This is an exact false-lead certificate, not just a failed numerical search.

Instead define the rational control

    g_lambda=lambda*R_j/(1+c_j*lambda*D_j).               (4.3)

It solves E(g_lambda)=0 exactly. One direct check uses

    E(R_j)=g*B*partial_r R_j,
    E(R_j-c_j*g*D_j)=-c_j*D_j*E(g),

where E denotes the derivation with the specified current g. The relation g=lambda*(R_j-c_j*g*D_j) then implies E(g)=0 whenever 1+c_j*lambda*D_j is invertible. Equivalently substitute (4.3) into (1.4). This rational Keller family has a genuine source pole; it is not an entire polynomial map.

### Every finite order lifts, and the full formal series exists

Expand

    g_lambda=sum_(n>=0) (-c_j)^n*lambda^(n+1)*R_j*D_j^n.  (4.4)

This is a formal lambda-adic series with POLYNOMIAL source coefficients of unbounded degree. Its truncation through lambda^L defines an actual polynomial Keller self-map over k[lambda]/(lambda^(L+1)), for every L>=1. The full series defines a morphism of lambda-adic formal affine spaces, with exact determinant one. It generally does not belong to k[[lambda]][x,y,z], which requires a finite uniform source-degree bound; it belongs to k[x,y,z][[lambda]]. This distinction prevents treating the formal solution as a finite polynomial formula.

There is NO polynomial family g in k[lambda,x,y,z] solving the same one-coefficient-line Keller problem with first derivative R_j at lambda=0: UASELFCTRL1, treating lambda as an identity parameter, forces g in k[lambda,sigma], while R_j is not a function of sigma. Thus these infinitesimal solutions are unobstructed to every finite order and formally, but do not integrate to an entire finite polynomial family inside this specified architecture. This is not a claim that broader multi-direction deformations cannot integrate them.

### The marked restrictions can be fixed at the same time

Replace lambda everywhere by lambda*I(sigma)^ell, ell>=1. Both the rational identity and formal construction remain valid, since the feedback derivation fixes sigma. Each finite-order polynomial map over k[lambda]/(lambda^(L+1)) agrees with the cubic on the ENTIRE target thickening V(I(a)^ell), not merely on listed points. Take I(T)=T*(3*T-1)*(4*T+1): all three filter planes and the complete empty/singleton/triple fibres remain unchanged. The same first-output/divisor argument as UAMARKEDDEF1 excludes hidden preimages. Hence exact finite-order Keller identities together with exact marked restrictions still do not certify global polynomial integration.

These auxiliary deformations may be inserted in UACANFIB1 over nilpotent parameter bases, or formally, keeping the trace selector. The idempotent relation beta^2=beta on the first-stage fibre and sigma_aux=beta/3 make I(sigma_aux)=0, so its complete marked fibre remains Spec(A). The old auxiliary-first/compiler-second factors remain. No new canonical atom, constant generic degree over arbitrary formal parameter fields, or improved cost is claimed.

## 5. Literature and the remaining general gap

[De Bondt, Proposition1.3](https://arxiv.org/html/1501.05168v1) gives the classical self-invariance criterion for quasi-translations. Its hypotheses concern maps x+H and the actual displacement H; it does not allow first integrals of a frozen response to be transplanted into an arbitrary Keller coefficient line. The explicit correction term B*g in (1.4) is why (4.2) fails.

The scalar first-order equation has a Burgers/Riccati form. [Lu2020](https://www.numdam.org/articles/10.5802/crmath.136/) treats analytic entire and meromorphic generalized inviscid Burgers equations; only its abstract/scope is used as related literature, not as a theorem dependency. Our argument is elementary Laurent algebra over every characteristic-zero field, with a direct source-plane pole test. Novelty and independent proof review remain open.

A remaining constructive candidate must escape both the earlier descended-time and root-polynomial factor classes. This theorem additionally shows that merely replacing a retained-signal scalar by an arbitrary own-source polynomial along one affine coefficient line supplies no new entire map. Curved coefficient families, multiple affine directions and genuine multi-block interaction remain outside the scalar theorem. The explicit rational denominator in (4.3) is now a concrete integration obstruction to address: cancel it through a proved broader whole-source construction, while retaining the complete canonical scheme and proving absolute atomicity. No degree sweep inside the closed scalar line is scheduled.
