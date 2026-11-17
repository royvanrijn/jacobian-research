# Sharp raw13 through degree seven for distinct nonzero diagonal weights (KCUB38)

Subsequent [KCUB40](EQUAL_NEGATIVE_FULL_DEGREE7.md) closes full111 through
degree seven by a changed curve/Euler-Wronskian proof. The current
degree<=7 diagonal-weight raw<=12 frontier is empty:only112 can collide,
with exactly13 terms. All larger frontier and full111-unresolved descriptions
below are historical stages;original timeout records are unchanged.

Subsequent [KCUB39](ZERO_WEIGHT_CLASSIFICATION.md) closes every one-zero-weight
three-variable class in all degrees. The current degree<=7 diagonal-weight
raw<=12 frontier contains only full111. All fourteen/twelve/ten-type
frontiers and unresolved zero-weight descriptions below are historical
stages, not current next steps; their solver records remain unchanged.

In three variables over a characteristic-zero field, consider Keller maps
with pairwise distinct nonzero source/output diagonal G_m weights in the
coordinates being counted. The minimum ordinary degree of a map with an
algebraic-closure collision is seven. Among maps of degree at most seven,
the minimum number of raw nonlinear monomial occurrences is thirteen.

The existing weight(1,-1,-2) representative attains both bounds. Every
collision in this restricted degree range has that weight type, up to
permutation and overall weight sign, and [KCUB15](WEIGHTED_COMPLETE_DEGREE7_CLASSIFICATION.md)
classifies its normalized map as a diagonal transform of the established
13-term representative. There is no new positive witness here.

The new work closes both remaining mixed types(3,-1,-2),(4,-1,-2) from
[KCUB37](WEIGHTED_124_125_EXCLUSION.md). The full diagonal-weight search for
a degree<=7 raw<=12 collision now has exactly ten retained types:
full(1,-1,-1), and nine(0,p,-q),gcd(p,q)=1,1<=p<=q,p+q<=7.
Those ten are unresolved, not constructed. Arbitrary source-frame changes,
nondiagonal actions in the counted coordinates, unweighted maps, and
higher-degree sparsity are outside this sharp raw-count result. Unrestricted support and degree minima, and dimension minima for other
special forms, remain open.

## Minimum dimension for this additional form

Allow dimension to vary while requiring a nontrivial diagonal G_m action
with distinct nonzero weights in the displayed coordinates. The minimum
collision dimension is three, with no degree cutoff. In dimensions one
and two, every active native relation lies in w-perp, whose dimension is
at most one. KCUB36's all-degree rank-one lemma therefore makes every
normalized map an automorphism. The explicitly rechecked three-variable
control attains dimension three. This is a bound for the additional form,
not an unrestricted dimension bound or a claim about other special forms.

## Complete source models and the finite cover

Distinct weights make the linear part diagonal; output scaling normalizes
it to identity and det JF=1 without changing degree or occurrences.
For p=3 or4 use the nonzero-a cover

    a=s^p, b=u/s, c=v/s^2, F=(s^p*f,s^-1*g,s^-2*h).

The full native models have15/14 nonlinear atoms, respectively, and26
Jacobian coefficient equations each. All atoms are retained. Exact
source differentiation gives

    J=f*{g,h}+g*{f,h}/p+2h*{g,f}/p=1,
    {P,Q}=P_u*Q_v-P_v*Q_u.                              (1)

A quotient zero-weight monomial u^i*v^j descends to
 a^((i+2j)/p)*b^i*c^j. The checker verifies every exponent is an integer
and compares the entire original-source determinant, including a=0.
The cover is not asserted to give polynomial invariant coordinates.

## Weight(4,-1,-2): complete triangular classification

Put T=a*c^2,V=a*b^2*c,U=a*b^4, with V^2=T*U. Complete weight modules are
F=(a*f,b*g,c*h+b^2*k(U)); f,g,h have constant one and may contain
T,V,U,T^2, while k=L+M*U. At b=0, exactly

    J=g0(T)*[h0(T)*(f0(T)+T*f0'(T))+2T*f0(T)*h0'(T)]=1.

Polynomial units and the nonzero leading scalar1+deg f0+2deg h0 force
f0=g0=h0=1. This removes all six T-only coefficients, retaining

    f=1+A*V+B*U, g=1+C*V+D*U, h=1+E*V+G*U, k=L+M*U.

The actual source determinant lies in K[U,V]. These two monomials are
algebraically independent, so formal coefficient comparison is valid.
The coefficients [U^0*V^j],j=1,2,3, of J-1 are

    2A+3C+2E, 4AC+3AE+4CE, 5ACE.

For any nonempty actual support of A,C,E, its size-th equation is the
product of active variables times
1+1_(A active)+2*1_(C active)+1_(E active). All seven cases have a nonzero
positive integer scalar, forcing A=C=E=0 over every characteristic-zero
field. Then

    J=(1+G*U)*(1+(2B+5D)*U+6BD*U^2)=1.

Units give G=0; the two remaining coefficients force B=D=0. Thus the exact
complete normalized class is

    F=(a,b,c+L*b^2+M*a*b^6).

Its polynomial inverse subtracts the shear, with arbitrary L,M. No external
theorem is needed for this weight type.

## Weight(3,-1,-2): full model and the L=0 boundary

Write its complete quotient model as

    f=1+mu*u*v+P*u^3+T*v^3+Q*u^2*v^2,
    g=u+L*v^2+M*u^2*v+N*u^4+A*u*v^3+B*u^3*v^2,
    h=v+S*u^2+H*u*v^2+J*u^3*v+K*u^5+W*v^4+X*u^2*v^3. (2)

If L=0, the second source output is divisible by b. At b=0 put Z=a^2*c^3.
The determinant factors as

    (1+A*Z)*[(1+W*Z)*(1+3T*Z)+3W*Z*(1+T*Z)]=1.

Equivalently use f0=1+T*Z,h0=1+W*Z and the leading scalar
1+2deg f0+3deg h0. Polynomial units give A=T=W=0.
Now put U=a*b^3,V=a*b*c. The full remaining source shapes are

    F1=a*(1+mu*V+P*U+Q*V^2),
    F2=b*(1+M*V+N*U+B*V^2),
    F3=c*(1+H*V+J*U+X*V^2)+b^2*(S+K*U).

The source determinant lies in K[U,V], and the coefficient of U^0 is

    d/dV[V*(1+mu*V+Q*V^2)*(1+M*V+B*V^2)*(1+H*V+X*V^2)].

It equals one. Integration at zero and polynomial units force all three
factors to one. This is formal coefficient comparison in algebraically
independent U,V, not a claim that every point with U=0,V!=0 is realized
by source coordinates. After removing those six coefficients,

    J=(1+J*U)*(1+(2P+4N)*U+5PN*U^2)=1

forces P=N=J=0. Precisely the shears(a,b,c+S*b^2+K*a*b^5) survive on L=0.

## Complete normalization for L!=0

A diagonal b conjugacy makes L=1 without roots. On this slice the residual
source conjugacy(a,b,c)->(a/t,t*b,t*c), with corresponding inverse output
scaling, preserves L and scales(mu,T) by t and(P,Q) by t^2. Therefore every
nonconstant first quotient belongs, over an algebraic closure, to one of
four complete strata:

| stratum | fixed first coefficients |
|---|---|
| mu1 | mu=1; T,P,Q arbitrary |
| mu0T1 | mu=0,T=1; P,Q arbitrary |
| muT0P1 | mu=T=0,P=1; Q arbitrary |
| muTP0Q1 | mu=T=P=0,Q=1 |

The last two use a square root of a nonzero field coefficient only to choose
a diagonal frame over the algebraic closure. No zero coefficient is
inverted silently; all earlier zero strata are separate. A contradiction
over the algebraic closure excludes the original field point too.

The three strata with mu=0 have exact full-ring rational units after six
constant-pivot eliminations each. The mu1 reduced model retains seven
variables and twenty coefficient equations; its changed ten-second
Groebner/lift query timed out. That timeout remains unchanged. The next
boundary step, not a larger cap on the same ideal, resolves the entire mu1
stratum.

## The complete mu1 boundary cover

In the native polynomial ring, the b=0 coefficients in powers of Z=a^2*c^3
are exactly

    3T+A+4W-2L*H-mu*L,
    4A*W+3T*A+6T*W-3T*L*H,
    6T*A*W.                                             (3)

On mu=L=1 eliminate H=(3T+A+4W-1)/2 using the first equation. Twice the
second becomes

    A*(8W+3T)-3T*(3T-1)=0, T*A*W=0.                    (4)

If T=0, (4) gives A*W=0. If T!=0, it gives A*W=0 again; for W=0 the
second equation gives A=3T-1, and for A=0 it gives T=1/3. Thus the four
following parameter charts cover every field point, including all overlaps:

    T=0,W=0; T=0,A=0; W=0,A=3T-1; A=0,T=1/3.

No saturation is needed: the third chart is allowed to include T=0 too.
All source coefficients not fixed by these necessary equations remain free.
Each chart has a full-ring rational unit after complete constant-pivot
elimination from the original native Jacobian equations.

## Exact constant pivots and seven multiplied units

Every reduction step selects an original native coefficient equation whose
derivative with respect to the eliminated variable is a nonzero rational
constant. Solving that equation gives a unique polynomial expression with
only rational constant denominators. The checker independently reconstructs
the equation at that step, verifies the pivot, verifies the solved expression,
substitutes it into every original equation and reconstructs every remaining
generator. No parameter denominator or coefficient-zero root is discarded.

| accepted case | free variables | constant pivots | generators | multiplier monomials |
|---|---:|---:|---:|---:|
| mu0T1 | 6 | 6 | 20 | 127 |
| muT0P1 | 5 | 6 | 18 | 275 |
| muTP0Q1 | 4 | 6 | 17 | 2 |
| mu1,T0W0 | 5 | 6 | 16 | 161 |
| mu1,T0A0 | 6 | 5 | 17 | 7 |
| mu1,W0,A=3T-1 | 5 | 6 | 18 | 386 |
| mu1,A0,T=1/3 | 4 | 7 | 17 | 138 |

All seven identities have sum_i multiplier_i*generator_i=1, over Q, with
1,096 multiplier monomials total and maximum multiplier degree seven.
The independent checker multiplies every term. It also verifies the six
constant pivots of the timed-out mu1 reduction as retained model provenance:
48 pivot steps overall,42 in accepted unit models. Primary Groebner results,
finite coefficient grids and generic parameter-field units are not premises.

Every nonlinear first quotient is now excluded. Thus every surviving
weight(3,-1,-2) Keller map has first source output a. The explicit imported
[GLDC10](../quartic-three-audit/DEGREE7_COMPONENT_EXCLUSION.md) excludes
collisions in this degree-seven linear-component sector. We need no new
plane classification or explicit inverse formula for its general tail.
The L=0 and weight412 subfamilies above have direct inverse formulas.

## Matching bound and current ten-type frontier

The complete [KCUB36 census](WEIGHT_RESONANCE_CENSUS7.md), together with
[KCUB37](WEIGHTED_124_125_EXCLUSION.md) and the new exclusions, leaves only
weight112 among distinct nonzero diagonal-weight collision maps through
degree seven. KCUB15 gives its exact13 count and degree-seven minimum.
For matching attainment, the checker reconstructs the existing representative
without trusting a manifest. Put v=1+2ab:

    F1=a-3a^2*b-a^3*c,
    F2=3v^2*a*c+(9v^2-6v-2)*b,
    F3=v^3*c+(6v^2+2v)*b^2.

Its full symbolic determinant is one, its linear part is identity, its
nonlinear row counts are(2,5,6), and its ordinary degree is seven. The
independent replay verifies the distinct points(0,1,-8),(-1,1,4) both map
to(0,1,0), and compares every nonlinear coefficient with the retained
raw13 input. Thus this is a matching existing control, not a new witness.

All smaller degrees are covered by the degree-seven exclusions, so the
minimum ordinary collision degree is seven even without an upper cutoff
on degree. The raw13 term minimum retains its explicit degree<=7 cutoff;
higher-degree sparsity is not settled. Both statements require the displayed
nonzero distinct diagonal-weight coordinates. General source-frame changes
can change counts, including for the known112 map, and are not excluded.

For the broader diagonal-weight raw<=12 objective, precisely ten weight types
remain: full111 and nine opposite-sign zero-weight types. Their existence
and their own sharp minima are unresolved. No unrestricted full profile is
newly excluded;62 excluded/12 remaining and all global positive bounds are
unchanged.

## Reproduction and trust boundaries

Run

    python verify_weighted3412.py

The [receipt](weighted-312-412/validation.json) pins both full sources,
the structural model, all eight reduced inputs, all seven unit certificates,
the old native timeout checkpoint, the new reduced timeout, the imported
proof scopes and the directly rechecked13-term upper control. See the
[artifact record](weighted-312-412/README.md). The original four native
five-second queries and the changed mu1 ten-second query remain terminal
timeout artifacts, not successful solver runs. No unchanged ideal was
rerun or given a larger cap.

Universal boundary/unit/degree arguments, normalization and field case
coverage are written mathematics, not formalized or externally independently
reviewed. GLDC10,KCUB15,KCUB36,KCUB37 retain their explicit theorem/import
boundaries; this replay does not reprove their universal arguments. No
novelty priority, unrestricted term optimum, new counterexample mechanism,
unrestricted degree/support/dimension optimum, or dimension optimum for
other special forms, is claimed.
