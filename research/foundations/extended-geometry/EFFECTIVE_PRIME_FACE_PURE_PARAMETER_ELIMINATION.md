# Effective parameter elimination for pure spherical phases

Put A_K=K[x,y,t]/(t^2+2xy-1), for a characteristic-zero field K, and let

    Lambda(f)=CT_x (1/2) integral[-1,1]
                  f(x,(1-t^2)/(2x),t) dt.

The integral is the rational polynomial rule. A phase f is PURE if Lambda(f^m)=0 for every integer m>=1. Its moments are denoted M_m. This note combines two existing results: the computable all-field sphere prefix and prime-face parameter rigidity. It computes necessary parameter polynomials and, in the one-parameter pure cases below, the exact finite GEOMETRIC parameter set. No cutoff, exceptional polynomial or candidate list is actually computed here.

## 1. Exact input and the two reused results

An input coefficient family is given by a finite presentation

    C=P/J0,    P=Q[T_1,...,T_r,c_1,...,c_n],

and an explicit polynomial class f in C[x,y,t]/(t^2+2xy-1), of degree at most d in x,y,t. The T_i are the retained parameters. Any finite list of required coefficient inversions is encoded by extra variables: h!=0 becomes wh-1=0. Likewise a regular parameter expressed as a rational function is encoded by its graph equation and the required inverse. Thus all computations below take place in a polynomial ring over Q; a boundary chart is not silently discarded.

Choose polynomial representatives mu_m in P of Lambda(f^m). They are exactly computable, for example by expanding powers and using

    Lambda(x^i y^j t^k)=0 unless i=j and k=2a,
    Lambda(x^i y^i t^(2a))
      =2^i*i!*(2a)!*(i+a)!/[a!*(2i+2a+1)!].              (1)

Changing representatives changes mu_m only modulo J0 and has no effect on any ideal used below.

We use the following previously proved facts with their precise scopes.

(U) GLOBAL_SPHERE_PURITY_RECOGNITION.md, Theorem ComputableUniformSpherePrefix and Section4: a terminating exact algorithm computes B(d), such that over EVERY characteristic-zero field, every sphere phase of degree at most d is pure if and only if its first B(d) moments vanish. This is a computable number, not a presently supplied numerical bound. The same source proves the equivalent reduced moment-ideal statement.

(P) PRIME_FACE_PARAMETER_RIGIDITY.md, Theorem1(b,c) and its specified applications: for a promised prime-face family, the parameter image of the PURE locus, and separately that of each FIXED tail locus, is contained in some proper rational algebraic subset. For r=1 these images are finite. The promise is the actual residue hypothesis there, not merely a statement that a generating series is D-finite.

We also use effective polynomial-ideal elimination over Q, proved in the existing narrow node EPA-GROEBNER-ELIMINATION-AND-NORMALFORM-TERMINATION, and the classical Nullstellensatz in the existing source EXT-MOMENT-LOCAL-HEIGHT-AND-NULLSTELLENSATZ. No ordinary-prime conjecture is used.

## 2. Necessary polynomials from a terminating prefix search

For an integer m0>=1 write

    I_infty,m0=J0+(mu_m : m>=m0) subset P,
    I_b,m0=J0+(mu_m : m0<=m<=b),    b>=m0.                (2)

Suppose the image of the fixed-tail locus is contained in a proper rational algebraic subset: there is a NONZERO g in Q[T_1,...,T_r] which vanishes at the retained parameter values of EVERY characteristic-zero field point killing I_infty,m0. This is exactly the relevant conclusion(P).

**Theorem 1.** The following algorithm terminates on this promised input. For b=m0,m0+1,..., compute the elimination ideal

    J_b=I_b,m0 intersect Q[T_1,...,T_r].                  (3)

Stop at the first nonzero J_b and output any nonzero H in it. The output is a necessary polynomial for every phase in the specified fixed-tail locus, over every characteristic-zero field. In particular m0=1 gives a necessary polynomial for all pure phases. There is no uniform assertion over an unbounded varying tail threshold.

Proof. Consider any prime ideal p containing I_infty,m0. The phase over Frac(P/p) is an actual characteristic-zero field point of this locus, so g belongs to p. Intersecting all such primes gives

    g in radical(I_infty,m0).

If I_infty,m0 is the unit ideal this conclusion is automatic, with g=1. By the definition of a radical, some positive power g^e lies in I_infty,m0. Membership in an ideal uses a FINITE sum. Consequently g^e already lies in I_b,m0 for some finite b. It then belongs to the elimination ideal(3), which is therefore nonzero. Exact Groebner elimination is a terminating computation at every fixed b, proving termination of the outer search.

If H is returned, H belongs to I_b,m0 and hence vanishes at every field point satisfying all the indicated tail equations. This proves necessity without a finite-prefix purity inference. QED.

This argument handles nonreduced ideals: it uses the power g^e, not an unsupported assertion g in I_infty,m0. It does not need(U), a D-finite recurrence, or a numerical bound for b. A finite ideal-membership certificate can be retained in the form

    H=sum_i A_i*F_i + sum_(m=m0)^b Q_m*mu_m,              (4)

where the F_i generate J0. Groebner computations can track these combinations. Alternatively, once H is found, enumerate coefficient-degree bounds for A_i,Q_m and solve the corresponding rational linear systems; membership guarantees that one of these finite searches succeeds. The certificate is checked by literal polynomial multiplication. Rational denominators may be cleared.

The roots or zeros of this first H may strictly overinclude actual pure parameters. It is a necessary condition, not a recognition theorem, and its first stopping index need not be a valid purity cutoff.

## 3. One complete pure prefix with arbitrary coefficient constraints

Compute B=B(d) by(U), and set

    I_B=J0+(mu_1,...,mu_B),
    I_infty=J0+(mu_m : m>=1).                            (5)

**Lemma 2.** These ideals have the same radical, even with the extra coefficient relations and inverse variables in J0. The field points of V(I_B) are exactly the pure phases in this coefficient family, over every characteristic-zero field.

Proof. The field-point assertion is precisely(U) applied after any coefficient evaluation; the degree in x,y,t stays at most d. For radical equality, take any prime ideal p containing I_B. The phase over Frac(P/p) has zero first B moments and is therefore pure. Hence every mu_m belongs to p, so I_infty is contained in radical(I_B). The opposite radical inclusion follows from I_B subset I_infty. QED.

This proof does not infer that ordinary moment ideals are equal. It also avoids assuming that a radical identity survives an arbitrary nonflat coefficient specialization: the all-FIELD statement(U) is applied directly to every prime quotient. Coefficient inversion is exact because wh=1 forces h to remain a unit at every such point.

## 4. Elimination computes the geometric pure-image closure

Compute

    J=I_B intersect Q[T_1,...,T_r].                      (6)

Use a block monomial order placing all eliminated c variables before all T variables. A Groebner basis then yields generators for J. To see the elimination step directly, a polynomial in J has leading monomial involving only T; a basis leading monomial dividing it also involves only T. For a block elimination order, the corresponding entire basis polynomial lies in Q[T]. Division of the original polynomial consequently uses only those retained basis polynomials. The reverse inclusion is immediate. The existing terminating Groebner completion theorem supplies the finite basis.

**Lemma 3.** Over any algebraically closed characteristic-zero field Omega, V_Omega(J) is the Zariski closure of the retained-parameter image of the pure coefficient locus. Thus(6) computes its exact closed algebraic envelope; it does not generally assert that the image itself is closed.

Proof. Over Omega, Lemma2 identifies pure phases with V(I_B). By the Nullstellensatz, the ideal of polynomials vanishing on its parameter image is

    radical(I_B*Omega[c,T]) intersect Omega[T].

Radicals commute with contraction. Also

    (I_B*Omega[c,T]) intersect Omega[T]=J*Omega[T].        (7)

Indeed J is the kernel of Q[T]->P/I_B; tensoring this linear map with the field Omega is exact and gives exactly(7). The vanishing ideal is therefore radical(J*Omega[T]), whose zero set is V_Omega(J). QED.

The ideal J itself need not be radical. Its computed generators nevertheless define the correct closed set. If a reduced defining ideal is desired, one may compute its radical by the existing exact characteristic-zero ideal algorithms. A single nonzero member of J is already a valid necessary parameter polynomial.

## 5. Exact finite sets in one parameter

Suppose r=1 and(P) proves that the pure parameter image is finite. Then J in(6) is NONZERO. Compute its monic generator h in Q[T]. If J=(1), output h_sf=1. Otherwise compute

    h_sf=h/gcd(h,h'),                                    (8)

with monic normalization. This is a finite rational polynomial computation.

**Theorem 4.** The roots of h_sf in Qbar are EXACTLY the retained parameter values realized by pure phases over Qbar in the specified coefficient family. Every such value has a realization with all phase coefficients in a number field. The polynomial h_sf=1 describes the empty set. For arbitrary characteristic-zero K, every pure K-phase satisfies h_sf(T)=0; this statement does not promise a K-rational realization of every root lying in K.

Proof. Lemma3 identifies V(h) with the closure of the geometric image. The image is finite by hypothesis and therefore already closed, so it equals V(h). Passing to h_sf changes neither the roots nor the zero set.

A Qbar point of the finite-type coefficient scheme has finitely many algebraic coordinates, hence lies in some number field. Conversely, if a pure phase originally has transcendental lower coefficients and an algebraic parameter value alpha, its coefficient fiber over Q(alpha) is nonzero. After any embedding Q(alpha)->Qbar, faithful scalar extension keeps this finite-type fiber nonzero. A maximal ideal and the Nullstellensatz give a Qbar point of that same fiber. Thus transcendental lower coefficients create no additional geometric parameter values. QED.

Membership of a parameter in this set means existence of SOME pure phase over Qbar with that parameter. It is not a sufficient test for purity of a specified phase having that value. It is also not a rational-point decision algorithm for its coefficient fiber.

For comparison, the family f=T+t has mu_1=T and mu_2=T^2+1/3. Theorem1 can stop at H=T, which leaves the false candidate 0. The complete valid linear-phase prefix gives I_B=(1), since 3*mu_2-3*T*mu_1=1, and Theorem4 returns the empty set. In the family f=T*t the complete moment ideal is(T^2), which is nonreduced; the exact parameter set is still cut out by its squarefree generator T. These examples distinguish necessary-prefix elimination, full pure recognition, and reduction of the parameter ideal.

## 6. The concrete quintic outputs

Apply the algorithms at d=5 to the EXACT supported families in PRIME_FACE_PARAMETER_RIGIDITY.md.

(I) f=P(y,t)+B*x*t^2, with j+r<=5 and j+2r<=7 in P; put A=[y*t^3]P and C=[y^3*t^2]P, and require BC!=0.

(II) f=P(y,t)+x*t^2*Q(t), degQ<=2, with j+r<=4 and 3j+2r<=9 in P; put A=[y*t^3]P, C=[y^3]P and B=[t^2]Q, and require BC!=0.

In both families the actual angular and overall-scalar action permits B=C=1 over an algebraic extension; the retained normalized parameter is lambda=A and the original invariant is T=A^2/(BC)=lambda^2. These operations preserve the prescribed supports and all moment vanishing. The PFP finite-image conclusions apply to the full supported families, including arbitrary lower coefficients.

Thus Theorem4 computes, separately for(I) and(II), a squarefree rational h_lambda whose roots are precisely the normalized pure parameter values. It also computes an EXACT polynomial for the original invariant T: adjoin T-lambda^2 and eliminate, or equivalently take the squarefree monic part of

    Res_lambda(h_lambda(lambda),lambda^2-T).             (9)

If h_lambda=1 use 1. The zero set of(9) is the image of the finite root set under squaring. Every original pure phase normalizes to one such phase after scalar extension; every normalized phase is already an original phase. Hence this is the exact geometric invariant-parameter image, not merely an unspecified finite superset. PFP separately excludes T=0,-2, so the computed squarefree polynomial is coprime to T*(T+2).

The same conclusion holds for the PFP even D=0 subfamily

    f=P(y,t)+x*t^2*(q0+q2*t^2),
    r even, j+r<=5 in P,
    D=[y^5]P=0, C=[y^3*t^2]P, B=q2, BC!=0.              (10)

Normalize B=C=1 and retain lambda=[y*t^4]P. The resulting exact squarefree polynomial for lambda^2 again excludes 0,-2. This is only the specified D=0 chart, not the entire even space.

For the full even family with BD!=0, PFP uses B=D=1 and retains

    alpha=A/(B*rho^2), beta=C/(B*rho^4), rho^6=D/B,

where A=[y*t^4]P, C=[y^3*t^2]P. Lemmas2--3 compute a nonzero elimination ideal J in Q[alpha,beta] defining the EXACT Zariski closure of the normalized pure-parameter image. Properness follows from PFP. Its generators, or any selected nonzero generator G, are effectively computable. The closure may be empty or zero-dimensional. No equality between this closure and the actual two-parameter image, and no realizability of every point of V(J), is claimed.

Theorem1 additionally computes a nonzero necessary polynomial separately for EACH FIXED tail threshold in any of these charts. This assertion uses PFP's fixed-tail properness directly. It does not apply the pure cutoff B(5) to a tail beginning at m0>1, compute an exact fixed-tail set, or produce one finite exceptional set for a union of unbounded thresholds.

## 7. Singular and exceptional charts

There is no specialization of one generic differential recurrence in these algorithms. All coefficient relations, including singular points, are retained in the polynomial moment ideals. The all-field prefix(U) already handles the singular-origin issue by its actual indicial test and generic-component algorithm. A generic telescoper whose coefficients vanish upon specialization is not substituted for(U).

Required nonvanishing conditions such as BC!=0 or BD!=0 are legitimate because they define the actual PFP chart. Other pivots are not automatic units. If a triangular elimination divides by pivot polynomials d_i(q), its generic chart must adjoin w*product_i d_i(q)-1, and all pivot-zero branches must be retained separately with their compatibility equations. Alternatively use the undivided moment ideal throughout. A rational-function-field computation in Q(q) alone does not cover algebraic q specializations.

Nonreduced moment components cause no difficulty in Theorem1, which uses a power, or in Lemma2, which asserts only radical equality. Exact equation(6) may retain nilpotence; squarefree reduction in one parameter removes it only at the final set level. No ordinary-ideal equality or nilpotence exponent bound has been inferred.

## 8. What this does not compute

The earlier global recognition theorem already supplies the terminating algorithm for B(d). This note supplies no new proof of that original effectiveness result. The bridge here makes the PFP exceptional-parameter conclusions EFFECTIVELY extractable, including an exact finite geometric PURE set in the one-parameter cases. Neither B(5), the output moment index, h_lambda, a polynomial for lambda^2, nor the two-parameter ideal J has been run to completion or displayed here. No practical complexity bound is given.

D-finiteness or bounds on differential order and coefficient degrees alone would not suffice: (z*D_z-N)z^N=0 permits an arbitrarily delayed first coefficient at fixed order and polynomial degree. The actual global result uses regular rational telescoping, a coefficient-specific singular indicial test, and a terminating coefficient-space radical algorithm. Likewise a generic recurrence or a computed finite prefix without the complete cutoff premise is not an all-order purity certificate.

No ordinary-prime existence claim, universal marker-safety criterion, rational-point decision or global quintic nonexistence conclusion follows. The canonical PFP proof's literal statement that it supplies no explicit polynomial/list/cutoff remains true.

## 9. Exact source comparison

The mathematical inputs reused here are GLOBAL_SPHERE_PURITY_RECOGNITION.md, Theorem ComputableUniformSpherePrefix and Sections3--4; PRIME_FACE_PARAMETER_RIGIDITY.md, Theorem1 and Corollaries2--4; and the narrow existing effective-Groebner and Nullstellensatz inputs named in Section1. The proofs in Sections2--6 explain the additional elimination bridge.

The retained Bostan--Lairez--Salvy primary is arXiv1301.4313v2. Proposition4 and Theorem12/Section7 concern REGULAR certificates for closed pole-free rational periods, including an affine singular-denominator deformation. They underlie the existing global theorem and are not being reimported as an unverified generic ODE. The actual official Singular4.4.1 eliminate documentation describes ideal contraction to a polynomial subring, and Stacks Theorem10.34.1 describes the finite-type Nullstellensatz. These are classical algorithmic and geometric tools, not new moment assertions. No priority claim, new computational checker, external human review or formal verification is asserted.
