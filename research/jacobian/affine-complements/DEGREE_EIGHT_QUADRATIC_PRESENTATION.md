# An exact quadratic presentation of the rooted degree-eight necessary core

14 November 2026. PF2RULED29 gives a global characteristic-zero quadratic
presentation of the unsolved degree-eight necessary problem. It has38
variables and39 equations of degree at most2 over each marked coefficient
field. Its27 auxiliary coefficients are uniquely recoverable by monic
polynomial division. Eliminating them recovers exactly the previous
necessary system, with the nonzero-kappa condition preserved.

This is an equivalent presentation, not an exclusion. Native finite-field
saturation separately returned the unit basis for one large-prime residue
in each marked sector. Those observations identify a concrete exact-
certificate target; they are not promoted to characteristic zero.

## Scope and notation

Keep the actual-global-root hypothesis and the monic degree-eight boundary
perturbation of [PF2RULED26](SQUARE_BOUNDARY_DEGREE_GAP.md). The global
four-parameter core and polynomial root-constant recovery of
[PF2RULED28](DEGREE_EIGHT_GLOBAL_CONSTANT_RECOVERY.md) remain valid.

In K=Q[y]/(P8), use the complete four-parameter family of actual S in V3
with symbol H and the cusp equation. Its jets S0,S1,S2 are affine linear
in a0,...,a3. The marked constants are fixed:

    h=v-s, lambda=3h/4+delta, r=S0'/h.

All S satisfy S0'(s)=0. The complete normalized boundary difference is

    R=h^8+c7*h^7+...+c2*h^2, L=R'/h.

Both divisions by h are exact. The six c parameters and the four a
parameters are the original ten high variables. Neither individual
critical-root multiplicities nor later pivot coefficients are restricted.
The original nonzero leading coefficient of R is restored by reversing
the scalar normalization from PF2RULED26; kappa here is its normalized
Jacobian constant, and nonvanishing is equivalent.

The critical polynomial is the actual normal Jacobian

    C=J1(S,T_R)
     =lambda^2*r*L'
      +[2h*(S2+v*lambda*r)-lambda^2*r']*L.

Its coefficients are bilinear in the root and perturbation parameters.
The checker also verifies this expression against the full J1 formula.

## Thirty-nine quadratic equations

Introduce three univariate polynomials with auxiliary coefficients:

    B=b0+...+b8*v^8,
    M=m0+...+m5*v^5,
    A=n0+...+n10*v^10.

Use the following identities and nonzero equation:

    C-2*S0'*B-M=0,                                    (1)
    2*M*S0+2*h*L*A-kappa=0,                           (2)
    kappa*invkappa-1=0.                               (3)

Identity(1) has20 coefficient equations in degrees0,...,19. Identity(2)
has18 in degrees0,...,17. Together with(3), there are39 equations.
The variables are

    4 root + 6 perturbation + 9 quotient + 6 remainder
    + 11 normal control + kappa + invkappa = 38.

Every equation has total variable degree at most2 over K: S and L are
affine linear, C is bilinear, and all products in(1)-(3) use two affine
or independent auxiliary factors. Coefficients in K are constants for
this degree assertion; an extra field-generator variable is not counted
as a quadratic unknown.

The [exact retained presentations](degree8_quadratic_presentations.json.gz)
contain every polynomial and its variable order in both coefficient
fields. The checker constructs the complete S family from its actual
Cox sections and exact cusp inverse, then verifies all39 degrees and
coefficient counts. It does not specialize the S or c parameters.

## Global equivalence, not a generic elimination

The leading coefficient of S0' is a fixed nonzero field scalar. For the
present normalization it is -12, so the leading coefficient of2*S0'
is -24. The nine high equations of(1) solve b8,...,b0 by a triangular
matrix with diagonal24. The six low coefficients then solve m0,...,m5
with diagonal -1. The remaining five rows are exactly the primary
high-remainder equations in degrees6,...,10. The checker compares each
of those five polynomials to the retained PF2RULED26 input.

Equivalently, divide C by2*S0'. Its quotient has degree<=8 and the
remainder has degree<=10. Identity(1) holds precisely when that remainder
has degree<=5; B is the quotient and M is the remainder. No coefficient
polynomial in a or c has been inverted.

Evaluation of(2) at v=s gives the unique value

    kappa=2*M(s)*S0(s).

The leading scalar of2*h*L is16. Its eleven high coefficient equations
solve n10,...,n0 by a triangular matrix with diagonal16. The remaining
rows are exactly the divisibility condition

    kappa-2*M*S0 divisible by h*L.

Multiplication by the constant2 does not change divisibility. The
numerator has degree<=17 and the divisor degree7, so every admissible
quotient A has degree<=10, matching the full normal-control bound in
PF2RULED26. No allowed normal coefficient is omitted.

Thus all27 auxiliary coefficients B,M,A,kappa are unique polynomial
functions of the original ten high variables when the necessary
conditions hold. Monic division makes these recovery maps polynomial
on the entire parameter space, including all rank-drop and critical-root
collision cases. The inverse variable in(3) is the same original
nonzero condition, not a new genericity assumption. Its comparison with
PF2RULED28's t*M(s)*Zstar-1 uses only the field unit2.

These observations prove a global isomorphism of the normalized
necessary presentations after eliminating the27 auxiliary variables.
Combining it with PF2RULED27/28 recovers exactly the same four-parameter
core and its saturation. This equivalence does not establish actual
section lifts or a constant source Jacobian on the whole surface.

## Separate modular evidence

The retained native F4 inputs are reductions of the complete primary
necessary core, not sampled values of the four core parameters. They
use the following prime and marked residue values:

| Sector | Prime | y residue | Recorded saturated basis |
| --- | --- | --- | --- |
| (2,7) | 1073741827 | 743393613 | [1] |
| (4,5) | 1073741827 | 1037248271 | [1] |

The [sector27 input](sources/degree8-native-saturation/sector27.ms) and
[output](sources/degree8-native-saturation/sector27.out), and
[sector45 input](sources/degree8-native-saturation/sector45.ms) and
[output](sources/degree8-native-saturation/sector45.out), are preserved.
The checker verifies that both y residues satisfy P8 modulo the prime,
reconstructs the full primary D,E,M(s) system, and matches every native
input polynomial. The outputs record the unit basis. An optional bounded
replay reruns the native solver and compares the result.

The solver handles the nonzero condition by native saturation of the
last polynomial M(s)*Z; it does not require an additional inverse
variable. One first-sector unsaturated native input also returned[1];
its [input](sources/degree8-native-saturation/sector27-unsaturated.ms)
and [output](sources/degree8-native-saturation/sector27-unsaturated.out)
are retained as further discovery evidence. No second-sector unsaturated
or characteristic-zero conclusion is inferred from it.

These computations suggest a possible exact unit certificate, but a
finite-field unit basis alone does not prove characteristic-zero
emptiness. No exact K-coefficient identity1 in the necessary ideal has
been extracted here. The finite-field reductions cover two residues,
whereas the characteristic-zero presentation equivalence covers all
sixteen complex embeddings. Those are different assurance scopes.

Short bounded expanded and quadratic Singular attempts, and a native
quadratic attempt, did not yield results within their caps. The preserved
[expanded input](sources/degree8-native-saturation/expanded_singular_attempt.py.txt)
and [quadratic input](sources/degree8-native-saturation/quadratic_singular_attempt.py.txt)
explain those approaches. Observation timeouts or interrupted calculations
are not treated as proof of emptiness. There is no live solver being
polled by this document.

## The remaining exact certificate gate

The next gate is to produce and replay an exact characteristic-zero
identity over K. One sufficient target is1 in the ideal(1),(2),(3).
If1 already belongs to the coefficient ideal(1),(2) before saturation,
that stronger identity also suffices; the first-sector unsaturated
observation is a reason to test this, not evidence that it is achieved.
A valid saturated certificate may instead prove a power of kappa lies
in the unsaturated ideal. All coefficients and hypotheses must be
retained and both marked fields covered.

Until such an identity or a valid global lift is obtained, the rooted
degree-eight nonzero-kappa core is open. The actual-global-root
hypothesis, general PF2RULED25 determinant locus and other rank-drop
strata, later normal coefficients throughu^9, contact cases, other
height-nine pairs, higher heights, all-height pole-four, field degree4,
surface finiteness and actual compatible embedding remain unresolved.
No whole marked candidate or quartic packet is excluded. The first
open pairs remain(5,9),(6,9),(7,9),(8,9); the full goal stays active.

## Replay and proof boundary

Run from the repository root:

```sh
sage -python research/jacobian/affine-complements/verify_degree8_quadratic_presentation.py
sage -python research/jacobian/affine-complements/verify_degree8_quadratic_presentation.py --replay-modular
.venv/bin/python research/jacobian/affine-complements/audit_pole_integration.py
```

The default [checker](verify_degree8_quadratic_presentation.py) verifies
the exact characteristic-zero equations, monic pivots, old-input
coefficient identities and recorded modular inputs. The optional command
performs two bounded native msolve runs; it is not a characteristic-zero
search. The [receipt](degree8_quadratic_receipt.json) records whether
that optional replay was actually performed, the binary hash and tool
versions. The [observations](degree8_native_saturation_observations.json)
state explicitly that no characteristic-zero exclusion is proved.

Global monic-elimination equivalence is a written proof supported by the
exact identities. No formal verification, independent whole-proof
implementation, external review, new imported theorem or novelty
priority is claimed.
