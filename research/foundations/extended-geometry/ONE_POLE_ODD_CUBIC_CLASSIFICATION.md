# One-pole odd-cubic-height purity classification

Date: 19 November 2026. The classification concerns the stated odd-cubic-height class; it makes no global claim about arbitrary prime-degree Image phases.

## 1. Scope and exact classification

Let K be any characteristic-zero field. Work in the sphere coordinate ring
`K[x,y,t]/(t^2+2xy-1)` with its normalized Gaussian/sphere angular functional Lambda. Consider the regular unit-center odd-cubic-height phase

    p=t/x*(C(x)+t^2 D(x)),  C,D in K[x], C(0)=-D(0)!=0.

Normalize by dividing the phase by C(0), so C(0)=1. Assume C has at most one geometric zero. If it has a zero then C=(1+lambda*x)^c with lambda in K*, c>=1. This representation descends to K: if C=sum_j c_j*x^j has degree c, its unique root is exactly -c_(c-1)/(c*c_c). Set U=1+lambda*x.

**Theorem OnePoleOddCubicPurity.** All positive pure moments vanish if and only if one of the following holds.

1. C=1,D=-1. Here p=2ty is one-sided and eventually safe for every fixed polynomial marker; its least degree is 2.
2. C=U^(N+1), N>=1, and the inversion branch holds:

       P(U)=(U^(2N)-1)/(U-1),
       D= -U^(N+3)*P(U)^2/(4N^2).

   Its exact least sphere degree (equivalently least radial-primitive balanced lift degree) is 5N+3.
3. C=U^(N+1), N>=1, and for some divisor d|N, q=N/d, the fixed-zero branch holds:

       P(U)=((2U^d-1)^q-1)/(U-1),
       D= -U^(N+3)*P(U)^2/(4N^2).

   Its exact least degree is 3N+3.

Both quotients P are literal polynomials. There are no pure phases with C=U. The nonconstant branches are unsafe: each has a radial-primitive non-proper-power original two-pair lift and a single fixed original marker nonzero at infinitely many positive even exponents. For the inversion branch marker Z already works at EVERY positive even exponent, by the retained hyperbolic family. For the fixed-zero branch one may always take `(R+lambda Z)^N`; no assertion that marker Z works is made.

Consequently odd prime degrees in this entire unit-center, one-pole, cubic-height class are exactly the primes congruent to 3 modulo 10 greater than 3. They arise from case 2 with N=(ell-3)/5 even. This constrains this class only. The higher odd-height family already realizes prime degrees 41 and 59 outside the cubic-height class.

## 2. Retained all-order endpoint and rational-descent inputs

The endpoint calculation and CharacterSpanDescent in [ODD_CUBIC_HEIGHT_LOW_DEGREE_CLASSIFICATION.md](ODD_CUBIC_HEIGHT_LOW_DEGREE_CLASSIFICATION.md) imply the following facts for every all-pure regular unit-center phase in this class. Choose

    X=sqrt(-D/C), X(0)=1; A=C/X; u=x/A=xX/C.

Then u and X are rational over K, u(0)=0,u'(0)=1, and

    D=-C^3*u^2/x^2.                                      (O1)

Write v(u) for the formal inverse branch x with v(0)=0. Let V'(x)=1/C(x), V(0)=0, a formal primitive, not necessarily rational. Purity is exactly

    V(v(u))-V(v(-u))=2u.                                 (O2)

Indeed Phi(u)=v'(u)/C(v(u)) has no positive even powers iff (O2) holds. Odd original powers vanish by t parity. Conversely rational u with u'(0)=1, (O1) polynomial, and (O2) supplies a pure phase.

For any fixed polynomial G(x) and even m>=2 the marked angular moment is

    Lambda(G p^m)=kappa_m [u^m]G(v(u))*Phi(u),
    kappa_m=int_0^1 t^m(1-t^2)^m dt !=0.                 (O3)

These are all-order identities, not consequences of a finite prefix.

## 3. Rational primitive and induced field involution

Assume C=U^(N+1), N>=1. Then

    V=(1-U^(-N))/(N*lambda).

Work temporarily over the algebraic closure k of K, and put E=k(u,V), as a subfield of k(U). Embed E in k((u)) using the specified inverse branch. Formal replacement u by -u preserves E, since (O2) sends the generators

    u -> -u, V -> V-2u.

This defines a genuine automorphism sigma of E: any algebraic relation is preserved by the injective formal embedding, and the replacement has square the identity. It is nontrivial since u is nonzero and the characteristic is zero.

The extension k(U)/k(U^N)=k(U)/k(V) is cyclic of degree N, with automorphisms U->zeta U for zeta^N=1. By the Galois correspondence its intermediate field E is k(W), W=U^d, for some d|N. This can also be checked directly: the subgroup of order d has invariant field k(U^d), because that rational map has degree d. Let q=N/d.

Every automorphism of k(W) is a fractional-linear transformation. An elementary proof uses [k(W):k(f(W))]=deg(f): if f has a rational inverse, the two degrees multiply to 1. The degree formula follows because for coprime polynomials A,B, A(T)-zB(T) is irreducible over k(z), by its linearity in z and coprimality. Thus write sigma(W) as a Mobius function.

On the formal branch W(0)=1, sigma(W)(u)=W(-u), so sigma(1)=1. Moreover

    u=(sigma(W)^(-q)-W^(-q))/(2N*lambda).                 (O4)

The pole restriction in (O1) says that u can have no finite U-poles except U=0. If a=sigma^(-1)(0) is a nonzero finite W-value, the first term in (O4) has a pole at every U with U^d=a, whereas the second is regular. Multiplication by U^(3N+3) in (O1) cannot remove that pole. This contradicts polynomial D. Hence a is 0 or infinity.

A nontrivial Mobius involution fixing 1 with sigma^(-1)(0)=0 is exactly

    sigma(W)=W/(2W-1).

With sigma^(-1)(0)=infinity it is exactly

    sigma(W)=1/W.

For completeness: an involution matrix can be taken with trace zero. In the first case its numerator is aW, denominator cW-a, and fixing 1 gives c=2a. In the second case its numerator is constant, and trace zero makes the denominator proportional to W; fixing 1 gives the displayed map.

The inversion branch in (O4) gives u=(U^N-U^(-N))/(2N*lambda), independently of d. In fact E=k(U^N), so its intermediate-field parameter is d=N. The fixed-zero branch gives

    u=((2U^d-1)^q-1)/(2N*lambda*U^N).                   (O5)

Both formulas have u=x+O(x^2). Substituting in (O1) gives exactly the two polynomial D formulas in the theorem. Conversely each displayed involution fixes the formal point W=1 and sends u to -u and V to V-2u, which proves (O2) and sufficiency. All displayed formulas have coefficients in K, so the algebraic-closure classification descends without an additional parameter extension.

## 4. Constant and logarithmic boundary cases

For C=1, equation (O1) makes u polynomial: it has no finite poles. Here V=x and E=K(u,V)=K(x). As above purity induces an involution of K(x) sending x to x-2u(x). This is polynomial, so as a fractional-linear map it is affine. It fixes 0 and has derivative -1 there. Hence it is -x, giving u=x and D=-1. The resulting p=2ty has strictly negative x/y rotation weight, so every fixed marker is eventually annihilated. The standard bilinear transfer gives the corresponding statement for arbitrary original markers.

For C=U, V=log(U)/lambda formally. Let x_+=v(u),x_-=v(-u). The rationality of u implies that these are algebraic over K(u). Inside their finite common field L, differentiating (O2) gives

    dlog((1+lambda*x_+)/(1+lambda*x_-))=2lambda du.

In characteristic zero the derivation extends uniquely through finite algebraic extensions and commutes with embeddings. Take the field norm to K(u). By its product-over-embeddings formula, g=Norm((1+lambda*x_+)/(1+lambda*x_-)) in K(u)* satisfies

    g'/g=2lambda [L:K(u)] !=0.

A nonzero rational function g has g'/g=O(1/u) at infinity, so it cannot have a nonzero constant logarithmic derivative. This contradiction excludes C linear. It uses no analytic branch or completeness premise.

## 5. Exact degrees and original lifts

The general least-degree formula is max(deg C,deg D+2). Indeed the unique normal form with t-degree at most 1 is

    p=t*(C+D)/x-2tyD.

Its two summands have different y-degrees, and replacing t^2 by 1-2xy never increases total degree. Thus this normal form realizes the least possible degree, exactly the displayed maximum. In the inversion branch deg P=2N-1, hence deg D=5N+1. In the fixed-zero branch deg P=N-1 with nonzero leading coefficient 2^q, hence deg D=3N+1. Thus the exact degrees in the theorem follow. No leading coefficient can cancel in characteristic zero.

For either nonconstant branch put e=deg D, b=e+2, c=deg C. With C_h,D_h the degree-c and degree-e homogenizations in R,Z, the exact original polynomial lift is

    F=T/Z*(R^(b-c) C_h(R,Z)+T^2 D_h(R,Z)).             (O6)

Here R=xi1*z1+xi2*z2, Z=xi1*z2, W=2xi2*z1, T=xi1*z1-xi2*z2, and T^2+2ZW=R^2. The quotient in (O6) is a literal polynomial after this relation: the constant-in-Z part is R^(b-2)(R^2-T^2)=2ZW R^(b-2), and all remaining terms are divisible by Z. It dehomogenizes exactly to p, without the optional clearing scalars used for the hyperbolic family.

At R=0 its value is the nonzero top D term, T^3 D_e Z^(e-1), so R does not divide F. Modulo T, F/T agrees on the nonempty open set RZ!=0 with R^(b-c)C_h/Z, which is not identically zero: R,Z are independent in the prime hypersurface T=0. Thus v_T(F)=1. Since T=xi1*z1-xi2*z2 is irreducible, this excludes every scalar proper power over field extensions. This simple T-order argument applies here because the starting height power is 1.

## 6. One fixed unsafe marker for every fixed-zero branch

Put z=U^(-d)=1+s. The fixed-zero coordinate and Jacobian factor become

    u(s)=((1-s)^q-(1+s)^q)/(2N*lambda),
    Phi=2(1+s)^(q-1)/((1-s)^(q-1)+(1+s)^(q-1)).

The polynomial u(s) is odd with nonzero linear coefficient -q/(Nlambda), hence has an odd formal inverse s(u). For G=U^N=(1+s)^(-q), the even part of G(v(u))*Phi(u), expressed in s, is

    E(s)=((1+s)^(-1)+(1-s)^(-1)) /
          ((1-s)^(q-1)+(1+s)^(q-1)).                   (O7)

It has simple poles at s=1 and s=-1: the denominator there is 2 for q=1 and 2^(q-1) for q>1, nonzero. If only finitely many even coefficients in (O3) were nonzero, E(s(u)) would be a polynomial in u^2. Composing with the polynomial u(s) would make E(s) a polynomial in s, contradicting those poles. Hence infinitely many positive even exponents have nonzero marked moments. The fixed original marker is G_h=(R+lambda Z)^N. Homogeneous factorial/angular normalization multiplies its angular moments by the nonzero factorial (bm+N+1)!, so they survive in the original Image contraction as well. Odd marked moments vanish by t parity.

For N=2,d=1, the formula gives C=U^3,D=-U^7,u=x/U, recovering the retained degree-9 odd-composite witness. In this case v*Phi=u, so marker x has zero moments at EVERY positive exponent. This is why the theorem uses the uniform marker U^N, rather than extrapolating the inversion marker.

## 7. Rational-primitive corollary and a squarefree obstruction

**Corollary RationalPrimitiveOddCubicPurity.** The preceding theorem is also a complete classification of regular unit-center pure odd-cubic phases for which the differential dx/C has a rational primitive. In particular the only possible odd prime degrees in this rational-primitive class are the primes congruent to 3 modulo 10 greater than 3.

To prove the scope equivalence, suppose V'=1/C is rational. If C is nonconstant of degree M, its zeros must have multiplicities m_i>=2, because a derivative of a rational function has no simple poles. At each such zero V has pole order m_i-1, and no other finite poles. The rational function V has a finite limit V_infinity and V-V_infinity has zero order M-1 at infinity, by differentiating its Laurent expansion. If h is the number of distinct roots of C, the degree of V is therefore sum_i(m_i-1)=M-h. The order of a fiber at one point cannot exceed the map degree, so M-1<=M-h, giving h<=1. (Equivalently, count numerator zeros after subtracting V_infinity.) A single root of multiplicity at least 2 integrates rationally, and constant C does too. This proves the equivalence over the algebraic closure and hence over K. No assertion that every pure odd-cubic phase has a rational primitive is made.

**Lemma SquarefreeCenterObstruction.** No regular unit-center all-pure odd-cubic phase has a nonconstant squarefree C, even when C has arbitrarily many roots.

Proof. Over the algebraic closure write

    1/C(x)=sum_i a_i/(x-r_i),  a_i=1/C'(r_i)!=0.

Rational-u descent still applies. In a finite common field L over k(u) containing x_+=v(u) and x_-=v(-u), differentiation of (O2) gives

    sum_i a_i dlog((x_+-r_i)/(x_--r_i))=2du.

Take the field trace, using the norm-product identity proved in the logarithmic case. This gives

    sum_i a_i g_i'(u)/g_i(u)=2[L:k(u)],
    g_i=Norm((x_+-r_i)/(x_--r_i)) in k(u)*.

Every logarithmic derivative on the left is O(1/u) at infinity, whereas the right is a nonzero constant. This contradiction proves the lemma. It extends the linear-center exclusion but does not classify mixed simple and repeated roots.

## 8. Sources and finite-control scope

The necessary rational-u descent and endpoint/residue identities are exactly the retained results in [ODD_CUBIC_HEIGHT_LOW_DEGREE_CLASSIFICATION.md](ODD_CUBIC_HEIGHT_LOW_DEGREE_CLASSIFICATION.md). The inversion branch and its exact every-even marker formula are already retained in [HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md](HYPERBOLIC_ODD_HEIGHT_IMAGE_WITNESSES.md). The present classification adds the exhaustive one-pole argument, boundary cases, and fixed-zero all-order marker-pole argument.

J.S. Milne, Fields and Galois Theory v5.10 (2022), Theorem 3.17 (printed p39) and Remark 5.47 (printed p84), [author PDF](https://www.jmilne.org/math/CourseNotes/FT.pdf), retained as [milne_fields_galois_v510.pdf](sources/milne_fields_galois_v510.pdf), were inspected for the precise finite Galois correspondence and separable norm product formula. The Mobius and logarithmic-derivative steps above are proved directly. No novelty or priority assertion is made.

The separate standard-library checker constructs the actual C,D coefficients, checks exact degrees and involution equations, compares raw interval moments with residue coefficients, and computes original four-variable factorial controls. It deliberately labels all finite checks as controls; the written field and pole arguments supply the universal quantifiers.


The [source ledger](one_pole_odd_cubic_source_ledger.json) pins the exact retained endpoint/descent and inversion-family inputs and the primary Galois/norm source. The [checker](verify_one_pole_odd_cubic.py) and [receipt](one_pole_odd_cubic_verification.json) bind the finite controls to these proof and source bytes.

```bash
python3 research/foundations/extended-geometry/verify_one_pole_odd_cubic.py \
  --write-receipt research/foundations/extended-geometry/one_pole_odd_cubic_verification.json
```
