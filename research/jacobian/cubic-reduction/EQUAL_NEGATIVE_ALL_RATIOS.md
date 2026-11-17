# Equal-negative weights: all ratios in the strict degree-seven census

KCUB32,15 November2026. Let p,q be coprime positive integers. For a
three-variable characteristic-zero Keller map with the same source/output
weights(p,-q,-q) and ordinary degree<=7:

- If(p,q)!=(1,1), there is no algebraic-closure collision, with no flag or
  support restriction. After tangent-block normalization, its first output
  is the first source coordinate.
- If(p,q)=(1,1), import the existing KCUB23 exclusion only for the
  affine-invariant flag strictly below(7,6,4). Its full degree7 class,
  including two independent degree7 outputs, remains unresolved.

Consequently every equal-negative ratio is excluded from the strict
three-variable degree7 census. A common nonzero integer multiple of the
weights makes no difference, so primitive normalization loses no class.
The new unequal-ratio proof covers the full degree7 classes, including
profiles outside the strict census. This is a weighted subcase theorem:
it closes no additional unrestricted whole profile. The partition remains
62 excluded/12 remaining, and all positive simplicity bounds remain
unchanged. No smaller counterexample or global optimum is produced.

The new arithmetic consists of30 generic Jacobian identities over ten
complete weighted models and14 independent unit identities for145
coefficient generators. The units have376 multiplier monomials, maximum
multiplier degree13. Every coefficient and product is reconstructed with
standard-library rational arithmetic. The weight-lattice, field/power,
normalization and imported plane arguments below are written, not
proof-assistant formalized or externally reviewed.

## Tangent normalization and the degree lattice

The weights distinguish the positive one-dimensional linear block and the
negative two-dimensional block. The Keller condition makes both blocks
invertible. A constant output normalization preserves weights, maximum
degree, collision existence and the invariant flag, and gives derivative I.
All outputs vanish at0 because the three weights are nonzero.

A monomial in the positive output has exponents

```
a^(1+q*j) * binary_form_of_degree(p*j),
```

and a negative output has

```
a^(q*j) * binary_form_of_degree(p*j+1).
```

Here j>=0, by positivity and gcd(p,q)=1. In both cases its ordinary degree
is1+(p+q)*j. Thus p+q>=7 permits only the normalized linear map through
degree7. The remaining unequal primitive ratios are exactly

```
p+q=3: (1,2),(2,1), two nonlinear layers of degrees4,7;
p+q=4: (1,3),(3,1), one layer of degree5;
p+q=5: (1,4),(2,3),(3,2),(4,1), one layer of degree6;
p+q=6: (1,5),(5,1), one layer of degree7.
```

The independent checker enumerates every compatible nonlinear monomial
through degree7 for these ten ratios,147 atoms in total, and reconstructs
the entire coefficient maps. The infinite p+q>=7 conclusion follows from
the written degree formula, not a bounded weight scan.

## Automatic collision if the positive output is nonlinear

Every normalized positive output is a times a factor. In the one- or
two-layer range it is

```
a*(1+a^q*Q(b,c)+a^(2*q)*T(b,c)),
```

where deg Q=p and deg T=2p, and absent layers are zero. If Q or T is
nonzero,1+Q+T is a nonconstant polynomial and has a zero over the algebraic
closure. Set a=1 at that zero. Its image has first coordinate0.
The plane a=0 maps as(0,b,c), so it supplies a distinct point with the same
image. Hence any Keller completion of a nonlinear positive output would
automatically be a counterexample.

We use [GLDC5](../small-component-gate/PROOF.md) in its linear-component
case: a degree<=7 three-variable Keller map with a linear output
combination has no collision. Its plane degree is within the explicitly
imported Moh100 bound. This is an existing written/imported premise, not
reproved by the new arithmetic. No three-variable quartic-invertibility
import is needed for the new unequal-ratio proof.

## One nonlinear layer

Write A=a^q and let

```
F=(a+a*A*Q, b+A*R, c+A*S),
deg Q=p, deg R=deg S=p+1.
```

Let {f,g}=f_b*g_c-f_c*g_b, D=R_b+S_c and d={R,S}.
Euler's identity gives adj(J(R,S))*(R,S)^t=d*(b,c)^t/(p+1).
The top coefficient is therefore

```
[A^3]det JF = (p+q+1)/(p+1) * Q*d.
```

If Q!=0, the Keller condition forces d=0. Homogeneous binary forms of the
same degree with zero bracket are proportional: dehomogenization reduces
their ratio to a rational function with derivative zero in characteristic
zero. Some nonzero output combination of F2,F3 is then linear. GLDC5
contradicts the automatic collision. If Q=0, F1=a is already linear and
GLDC5 excludes a collision directly. This treats all eight one-layer
ratios. The generic first/top identities are independently checked for
every one of their complete coefficient models.

## Two layers: complete leading reduction

Only ratios(2,1) and(1,2) remain. Write

```
F1=a+a*A*Q+a*A^2*T
F2=b+A*R+A^2*U
F3=c+A*S+A^2*V,
```

where deg(Q,T,R,S,U,V)=(p,2p,p+1,p+1,2p+1,2p+1).
The independent checker differentiates the full28-atom model for(2,1)
and19-atom model for(1,2). It reconstructs these exact identities:

```
[A]det JF = R_b+S_c+(q+1)*Q
[A^6]det JF = (2*p+2*q+1)/(2*p+1) * T*{U,V}.
```

When T=0, the next top coefficient is

```
[A^5]det JF = (2*p+q+1)/(2*p+1) * Q*{U,V}.
```

If the positive output is nonlinear, at least one of Q,T is nonzero.
Thus either U=V=0, or U,V are proportional homogeneous binary forms.
In the latter case a coupled linear source/output conjugacy in(b,c)
puts the nonzero high-degree tail in F2 and gives V=0, while preserving
tangent identity and all degree bounds. This covers both the strict
single-leading-row and the nonstrict two-leading-row cases.

When U=V=0 and T!=0,

```
[A^4]det JF = (p+2*q+1)/(p+1) * T*{R,S}.
```

Hence R,S are proportional, giving a linear output combination and the
GLDC5 contradiction. With T=0 this is the one-layer argument. If Q=T=0,
F1=a is linear. All no-high-tail boundaries are therefore covered.

We may assume V=0, U!=0. If S=0, F3=c is linear and GLDC5 applies.
The remaining cases have S,U nonzero and either T=0,Q!=0 or T!=0.

## The correct homogeneous power relations

For T=0 the top coefficient simplifies, using homogeneous Euler identities,
to

```
[A^4]det JF = ((p+q+1)*Q*{U,S}-q*U*{Q,S})/(p+1).
```

For T!=0 it is

```
[A^5]det JF = ((p+2*q+1)*T*{U,S}-q*U*{T,S})/(p+1).
```

All these identities are checked on the full generic coefficient models.
The first gives {Q^q/U^(p+q+1),S}=0. The second gives
{T^q/U^(p+2*q+1),S}=0, in the rational function field.

For completeness, if H,S are nonzero homogeneous rational functions,
S has degree e>0 and H has degree -n*e, then {H,S}=0 and Euler imply

```
grad H = -n*(H/S)*grad S,
```

so grad(H*S^n)=0. In characteristic zero this product is a nonzero field
constant. Applying this observation gives the exact power relations

```
T=0:    Q^q * S^(2*p+q+1) = constant * U^(p+q+1)
T!=0:   T^q * S^(2*p+2*q+1) = constant * U^(p+2*q+1).
```

The degree exponents are essential. In particular the quartic-first
(1,2) relation is Q^2*S^5=constant*U^4.

Unique factorization of binary forms over the algebraic closure gives
precisely the following possibilities, with all displayed scalar
coefficients nonzero and L,M nonzero linear forms:

| Ratio | Case | Forced forms |
| --- | --- | --- |
|(2,1)|T=0|Q=alpha*L^2, S=beta*L*M^2, U=gamma*L^2*M^3|
|(2,1)|T!=0|T=alpha*L*M^3, S=beta*L^2*M, U=gamma*L^3*M^2|
|(1,2)|T=0|Q=alpha*L, S=beta*L^2, U=gamma*L^3|
|(1,2)|T!=0|T=alpha*L^2, S=beta*L^2, U=gamma*L^3|

In the mixed T!=0 cases Q remains an arbitrary form of degree p.
The checker independently enumerates all possible line-factor valuations
and degree partitions in these four power equations. Coincident lines
are included. No squarefree or generic-root assumption is made.

## Every source normalization stratum

Once the high tail is in the b direction, coupled upper-triangular linear
changes in(b,c) preserve that direction and the identity linear part.
For an ordered pair L,M they have exactly five needed representatives:

```
(L,M)=(b,b),(b,c),(b,b+c),(c,c),(c,b).
```

Indeed a line not proportional to c can be sent to b. Its stabilizer is
diagonal up to rescaling the line, putting the other line in b,c or b+c.
If the first line is proportional to c, the other is either also c or can
be sent to b by a triangular shear. All scalar factors are retained as
arbitrary nonzero parameters. For a single line there are two charts,
L=b and L=c. This covers every zero-coefficient and coincident-line case.

The first Jacobian coefficient forces

```
R_b=-S_c-(q+1)*Q.
```

Thus each chart retains exactly

```
R=integral_b(-S_c-(q+1)*Q)+r*c^(p+1),
```

with r arbitrary. Only division by positive integers is used in this
integration, not division by a parameter. The scalar product alpha*beta*gamma
is saturated with a fresh variable sigma, imposing
sigma*alpha*beta*gamma-1=0. This imposes no extra condition on any of the
nonzero-form strata and excludes no zero stratum that was not already
handled separately.

There are five quartic-first and five mixed charts for(2,1), and two of
each for(1,2): fourteen charts in total. The remaining full determinant
coefficient equations and saturation have an explicit polynomial unit
identity in every chart. These are reconstructed independently from the
original three-variable map, not merely read from primary Groebner output.

| Ratio/case | Charts | Unit multiplier monomials |
| --- | ---: | ---: |
|(2,1),T=0|5|16|
|(2,1),T!=0|5|307|
|(1,2),T=0|2|7|
|(1,2),T!=0|2|46|
|Total|14|376|

All145 generators,376 multiplier terms and fourteen sums equal to1 are
checked over the rationals. The maximum multiplier degree is13. Hence no
Keller completion exists in any remaining nonlinear-first case. Together
with the handled linear/no-high-tail boundaries, the new full unequal-ratio
exclusion follows. Once the first output is linear, GLDC5 excludes any
algebraic-closure collision.

## Strict census and the remaining equal-ratio exception

The [KCUB23 theorem](../equal-negative-completion/PROOF.md) separately
excludes the primitive equal ratio(1,1) in the strict flag range below764.
Combining it with the new full unequal-ratio exclusion closes that strict
range for every positive equal-negative ratio.

The full weight(1,-1,-1) degree7 class is not closed. Its two independent
degree7 outputs can put the flag in77*, outside the strict target. The
older chartII timeout and unreplayed broader chartI certificate remain
exactly as recorded; neither is used or silently resolved here. A new
method would be required before revisiting those unchanged experiments.

No additional whole unrestricted profile is excluded. The twelve remaining
profiles, global3D degree[5,7], raw13, collision-normalized[10,48],
identity-linear-part cubic12/36, free-linear-part cubic31/17 and
cubic-linear89 bounds remain unchanged. A strict low-profile candidate
must leave equal-negative weights, but a broader sparsity candidate in the
unclosed weight11177* sector is not ruled out by this theorem.

## Independent replay, correction history and assurance

Run

```
.venv/bin/python research/jacobian/cubic-reduction/verify_equal_negative_allratios.py
```

The checker independently enumerates the weighted monomials, reconstructs
all ten generic maps and30 structural identities, checks ten valid
non-affine triangular maps with both inverse compositions, enumerates the
four valuation problems, rebuilds the fourteen normal charts and all145
generators, and multiplies the376-term unit certificates exactly. It uses
standard-library Fraction arithmetic and the pinned sparse polynomial
helper; SymPy, Singular and their elimination/factorization results are not
trusted. The complete written weight, field and Borel-cover reasoning is
not formalized. GLDC5/Moh100 and the KCUB23 strict-only theorem remain
explicit imported premises with their original trust boundaries.

- [Structural producer](equal_negative_structural.py) and [models](equal-negative-smallweights-v2/structural-models.json)
- [Corrected normal-chart producer](equal_negative_smallweights_v2.py)
- [Primary exact certificate producer](certify_equal_negative_smallweights_v2.py)
- [Independent checker](verify_equal_negative_allratios.py) and [receipt](equal-negative-smallweights-v2/validation.json)

An initial211-only screen is retained as discovery history. The next
first-pass17-chart list included three nonparallel quartic-first(1,2)
charts that are not on the correct leading power locus. Their unit ideals
did not constitute the normalization proof. The
[historical first-pass note](equal-negative-smallweights/README.md) preserves
all models, scripts, logs and certificates. The corrected v2 models use
the actual power relation and assert that the top coefficient vanishes
on every normal template. The independent valuation enumeration and full
Jacobian reconstruction verify the correction. No theorem is inferred
from the discarded list or from a bounded point search. There was no
mathematical timeout in the new complete proof.
