# A smooth algebraic cubic whose first nonzero Gaussian moment is fourteen

Current quantitative refinement, 15 November 2026: the [exact Euler representative proof](EULER_REPRESENTATIVES_GLOBAL_297.md) gives global pure297, odd296, conic/full-W-linear258 and all-homogeneous258. Ordinary-node complements are297/271/244, with homogeneous nonisotropic-node101; cusp255/tangent284 are unchanged. Current intervals are global[20,297], odd[20,296] and all-homogeneous[14,258]. Earlier353/336 bounds below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I297 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [joint-degree actual-seed proof](JOINT_DEGREE_GLOBAL_353.md) gives global pure353, odd352, full W-linear336, rank-one140, mixed-parent255 and cusp255/tangent284. Current least-prefix intervals are global[20,353], odd[20,352] and all-homogeneous[14,336]. Earlier590/888 bounds and rank-one425/mixed-parent320 estimates below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I353 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [balanced logarithmic seed proof](BALANCED_LOGARITHMIC_GLOBAL_590.md) gives global pure590, all-homogeneous336, conic336 and full W-linear425. Current intervals are global/odd[20,590] and all-homogeneous[14,336]; smooth[20,94] and homogeneous-smooth[14,56] are unchanged. Earlier888/1855 bounds and frame-plus888 statements below retain their historical scopes. Leastness,20/14 sufficiency, expanded I590 cofactors and independent review remain open.

Current global refinement,15 November2026: the [rescaled singular seed proof](RESCALED_SINGULAR_GLOBAL_888.md)
now proves888 sufficient as a global pure-only prefix, without an extra
frame condition. Current intervals are global[20,888], odd[20,888] and
all-homogeneous[14,888]. Earlier1855 bounds and reservations about888
alone below are historical. The new proof remains internal and
independently unreviewed; leastness and expanded global cofactors are open.

Current cubic refinement,15 November2026: the [graded seed proof](GRADED_SEED_SMOOTH_CUBIC_PREFIXES.md)
gives sufficient smooth-leading94 and homogeneous-smooth56 gates, with
current intervals[20,94]/[14,56]. Earlier188 bounds remain valid coarser
evidence; the analytic input and independent-review boundary are unchanged.

15 November 2026. This supplies a concrete algebraic witness for the
homogeneous lower bound in the [invariant argument](ODD_CUBIC_CLASSIFICATION_AND_PREFIX_LOWER_BOUND.md).
It improves the homogeneous-smooth lower bound from six to fourteen:

    14 <= least homogeneous-smooth ternary cubic pure prefix <= 56.

The subsequent [odd twentieth-moment witness](ODD_CUBIC_TWENTIETH_MOMENT_WITNESS.md)
together with the subsequent graded upper gives interval[20,94] when
lower terms are allowed. This
earlier homogeneous witness and its first nonzero order remain unchanged.

The upper bound is the retained [smooth-leading theorem](COERCIVE_GRADIENT_MOMENTS.md),
including arbitrary lower terms. This note proves the new lower bound
directly, without using the retained analytic classification of pure-null
cubics. It does not prove fourteen sufficient on any whole cubic class.
The global interval[20,1855] is unchanged; the odd19-zero witness has now
been constructed in the subsequent note linked above.

The [certificate](homogeneous_prefix_witness_certificate.json) specifies
six algebraic coefficients by rational polynomial equations and a rational
isolating polydisc. The separate [verifier](verify_homogeneous_prefix_witness.py)
and [accepted receipt](homogeneous_prefix_witness_verification.json) check
the inequalities in the proof below with outward-rounded complex ball
arithmetic. These are certified intervals, not a tolerance-based decision
that small residuals are zero. Independent replay, external review and
formal verification remain outstanding. The [source and failure audit](homogeneous_prefix_witness_source_audit.json)
binds inputs, library documentation and the preserved initial parser failure.

## The polynomial and its exact algebraic coefficient definition

Let Z=(X+iY)/sqrt(2), W=(X-iY)/sqrt(2), T=U for independent standard
real Gaussians X,Y,U. Thus E(ZW)=E(T^2)=1 and all other pairings vanish.
Set

    P=Z^3+ZWT+b W^3+c Z^2 T+d W^2 T
        +e(Z^2W-2ZT^2)+f(ZW^2-2WT^2)+h(T^3-3ZWT). (1)

For even n define the rational polynomial F_n(b,c,d,e,f,h) by

    F_n=E(P^n)/(3n+1)!!
       =sum_r r!*(3n-2r-1)!!
          *[Z^r W^r T^(3n-2r)]P^n /(3n+1)!!.     (2)

The sum is over0<=r<=3n/2, with (-1)!!=1. Formula(2), or ordinary sparse
polynomial multiplication followed by Wick contraction, specifies every
coefficient exactly. No quadrature or absolute moment occurs here.

The six coefficients in(1) are the unique root of

    F_2=F_4=F_6=F_8=F_10=F_12=0                 (3)

in the closed complex sup-norm polydisc of radius10^-40 centered at
the six exact rational complex numbers in the certificate's `center`.
Those decimal strings are interpreted as exact rationals. For orientation
only, the center starts as follows; these shortened decimals do not define
the isolating polydisc:

| Coefficient | Approximate value |
|---|---|
| b | -0.038004435478 - 0.115229219649 i |
| c | 1.234195555043 + 0.070003860703 i |
| d | -0.224811608459 + 0.138546678913 i |
| e | -0.464550029180 + 0.377806306621 i |
| f | -0.104632474480 - 0.099158038695 i |
| h | 0.135180075506 - 0.324522313681 i |

Existence and uniqueness in this region are proved next. The Jacobian is
invertible there, so the root is an isolated simple solution of a rational
polynomial system. Its coordinates are algebraic: the Jacobian-nonzero
part of the zero set is zero dimensional and has finitely many points;
elimination over Q gives an algebraic equation for each coordinate. Thus
(1) specifies a polynomial over Qbar also in the original Gaussian variables.
No minimal polynomials or radical expressions for its coefficients are claimed.

## A contraction certificate proves the moment equations exactly

Let x0 be the rational center, r=10^-40, F the six-component map in(3),
and B the exact rational complex six-by-six matrix stored as
`preconditioner`. Let D={x: max_i |x_i-x0_i|<=r}.

The verifier evaluates F(x0) and encloses J_F(x) for every x in D. It
uses the slightly larger rectangular product with real and imaginary
radius r per coefficient, so its enclosure includes all of D. It proves

    eta=||B F(x0)||_infinity < 10^-60,
    q=sup_(x in D) ||I-B J_F(x)||_infinity < 10^-20. (4)

The matrix norm is the maximum absolute row sum, matching the complex
sup norm. In particular q<1/2 and eta+q*r<r. The verifier retains exact
rational upper bounds for eta and q, not just their printed decimals.

For G(x)=x-BF(x), integration of the polynomial derivative along the
straight segment joining any two points of the convex polydisc gives

    ||G(x)-G(y)||_infinity <= q*||x-y||_infinity.

Also ||G(x0)-x0||<=eta, hence G maps D strictly into itself. Iterating
G from x0 gives a Cauchy sequence: successive distances shrink by q,
so their geometric-series tail tends to zero. Completeness of D gives
a fixed point, and the same inequality proves uniqueness.

The bound ||I-BJ_F(x0)||<1 makes BJ_F(x0) invertible by its convergent
geometric inverse series. Therefore B is invertible. The fixed-point
equation BF(x)=0 is consequently F(x)=0. The bound likewise makes
J_F invertible at the root. This proves all six equations(3) exactly.

The same rigorous evaluation encloses the next moment throughout D:

    2*10^-11 < Re F_14 < 4*10^-11,
    -9*10^-11 < Im F_14 < -8*10^-11.             (5)

More closely F_14 is approximately
2.8720448772595*10^-11 - 8.2996500971797*10^-11 i.
Since E(P^14)=43!!*F_14, the fourteenth moment is nonzero. The cubic
is odd, so every odd moment vanishes by Gaussian sign symmetry.
Equations(3) and(5) therefore prove that its first nonzero pure moment
is exactly fourteen, with moments1 through13 all zero.

## Smoothness is certified without a discriminant formula

Form the15-by18 matrix of the linear map

    (quadratic polynomials)^3 -> quartic polynomials,
    (A,B,C) -> A*P_Z+B*P_W+C*P_T.                (6)

The code lists each monomial basis explicitly. It encloses the determinant
of the fifteen columns selected in the certificate, uniformly on D. The
accepted enclosure is approximately

    -47.3541558855274 + 34.5617621326157 i,

with real and imaginary error less than10^-30. Its absolute lower bound
is positive. Hence(6) is onto at the root. A nonzero common zero of the
three partial derivatives would make every quartic vanish at that point,
including Z^4,W^4,T^4, which is impossible. Thus the projective plane
cubic P=0 is smooth. This is a direct rank certificate, not a genericity
assumption or a numerical discriminant check.

As an auxiliary check, (2*partial_Z partial_W+partial_T^2)P=2T
identically in the six parameters. The initial search used this slice
because the retained pure-null classification makes an anisotropic
Laplacian a useful obstruction. Acceptance of the final witness instead
uses the direct contraction and nonzero-moment proof above.

## Reproduction, search boundary and remaining work

Run the acceptance gate from the repository root:

```sh
.venv/bin/python research/moments/vanishing/verify_homogeneous_prefix_witness.py
```

The bounded [numerical exploration](homogeneous_prefix_witness_exploration.json)
found three approximate candidates in nine starts. The
[builder](build_homogeneous_prefix_witness.py) refined the first candidate
and proposed rational data. Neither step is the proof. The verifier
recomputes every acceptance inequality from the stored rational data,
and its integer Wick regressions independently match symbolic powers.
The first checker failed before any inequality because `fmpq` rejected
decimal strings; its input and failure receipt are preserved. Exact
`Fraction` parsing fixes that issue without changing the certificate.

The arithmetic semantics are documented in the primary
[Acb](https://python-flint.readthedocs.io/en/latest/acb.html),
[Arb](https://python-flint.readthedocs.io/en/latest/arb.html), and
[complex matrix](https://python-flint.readthedocs.io/en/latest/acb_mat.html)
references, version0.9.0 matching the installed wrapper. No theorem on
Gaussian moment vanishing is imported from these software references.

The old sixth-moment root-of-unity example remains useful evidence for
its all-dimension construction; it is superseded as the strongest smooth
ternary lower-bound witness. The homogeneous thirteen-zero witness and
its first nonzero order are now concrete. The odd nineteen-zero witness is also now constructed, with first
nonzero order20 and smooth leading form. The least global/homogeneous/
smooth thresholds, expanded global radical
cofactors, independent review, exceptional coefficient schemes,
term-count minima and lower-power GVC remain separate unfinished work.
