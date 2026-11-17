# An algebraic odd cubic whose first nonzero Gaussian moment is twenty

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

15 November 2026. This constructs the odd nineteen-zero witness left open
by the [invariant-ring argument](ODD_CUBIC_CLASSIFICATION_AND_PREFIX_LOWER_BOUND.md).
The cubic leading form is smooth, so the current least-prefix intervals are:

| Class | Current interval |
|---|---:|
| All ternary polynomials of degree at most three | [20,1855] |
| Odd ternary polynomials of degree at most three | [20,1854] |
| Smooth-leading ternary cubics, with arbitrary lower terms | **[20,94]** |
| Homogeneous ternary cubics | [14,1854] |

The [global1855](GLOBAL_TERNARY_CUBIC_PURE_PREFIX.md) and
[smooth94](GRADED_SEED_SMOOTH_CUBIC_PREFIXES.md) upper proofs are separate;
the original coarser smooth188 proof remains valid.
The new polynomial is not homogeneous. It does not increase the
homogeneous lower bound from the retained [fourteenth-moment witness](SMOOTH_CUBIC_FOURTEENTH_MOMENT_WITNESS.md).
Neither twenty-moment sufficiency nor any least threshold is proved.

The [certificate](odd_prefix_witness_certificate.json) gives exact
rational equations and an isolating region for nine algebraic coefficients.
The separate [verifier](verify_odd_prefix_witness.py) and
[accepted receipt](odd_prefix_witness_verification.json) prove the required
interval inequalities. The lower-bound construction is independent of
the retained analytic all-order classifications and invariant-ring
existence argument. It uses Wick's rule, a written contraction proof and
rigorous Acb arithmetic. Independent replay, external review and formal
verification remain outstanding. The [source audit](odd_prefix_witness_source_audit.json)
binds the calculation, previous frontier, software semantics and literature scope.

## A nine-coefficient polynomial with exact algebraic specification

In circular Gaussian coordinates Z=(X+iY)/sqrt(2), W=(X-iY)/sqrt(2),
T=U, where X,Y,U are independent standard real Gaussians, put

    P=Z^3+ZWT+b W^3+c Z^2 T+d W^2 T
       +e(Z^2W-2ZT^2)+f(ZW^2-2WT^2)+h(T^3-3ZWT)
       +uZ+vW+wT.                               (1)

Thus E(ZW)=E(T^2)=1 and all other pairings vanish. Define

    F_n=E(P^n)/(3n+1)!!
       =sum_(i,k>=0) i!*(2k-1)!!
          *[Z^i W^i T^(2k)]P^n /(3n+1)!!.        (2)

Only finitely many coefficients occur. Formula(2) is a rational polynomial
in the nine coefficients for each n. Unlike the homogeneous formula, its
sum includes *all* total degrees occurring in P^n; one must not restrict
to degree3n. The denominator is just a fixed positive integer scaling of
each equation, not a radial factorization of this inhomogeneous moment.

The coefficient vector (b,c,d,e,f,h,u,v,w) is the unique root of

    F_2=F_4=F_6=F_8=F_10=F_12=F_14=F_16=F_18=0  (3)

in the closed complex sup-norm polydisc of radius10^-40 about the nine
exact rational complex entries of `center` in the certificate. Every
decimal string there is an exact rational number. The following shortened
values are only for inspection, not for defining the isolating region:

| Coefficient | Approximate value |
|---|---|
| b | -0.120897625083 - 0.047806024445 i |
| c | 0.091103334551 - 0.104476634561 i |
| d | -0.054137022998 - 0.004050541471 i |
| e | -0.055555276058 - 0.000730098624 i |
| f | -0.014357166540 + 0.007206493210 i |
| h | 0.090450167193 - 0.010375049211 i |
| u | -1.554936812272 - 0.080273226879 i |
| v | -0.195925591048 - 0.242592255013 i |
| w | -0.404286624930 - 0.217942295767 i |

## The contraction inequalities and nonzero twentieth moment

Write x0 for that rational center, r=10^-40, F for the nine-component
polynomial map in(3), and B for the exact rational complex nine-by-nine
`preconditioner` in the certificate. For D={x:max_i|x_i-x0_i|<=r}, the
verifier proves

    eta=||B F(x0)||_infinity < 10^-65,
    q=sup_(x in D)||I-B J_F(x)||_infinity < 10^-15. (4)

The matrix norm is the maximum absolute row sum. The calculation encloses
J_F on the larger rectangular product of real and imaginary radius r,
so it includes every point of the complex polydisc. The receipt retains
exact rational upper bounds. In particular eta+q*r<r and q<1/2.

For G(x)=x-BF(x), integrate its derivative along a straight segment in D.
Equation(4) gives ||G(x)-G(y)||<=q||x-y|| and ||G(x0)-x0||<=eta.
Thus G maps D strictly into itself and contracts distances. Iteration
from x0 is Cauchy, since successive increments have a geometric bound.
Completeness gives a unique fixed point in D. Also BJ_F(x0) is invertible
by the geometric-series inverse of I-(I-BJ_F(x0)), so B is invertible.
Consequently the fixed point solves F=0 exactly. The same bound makes
J_F invertible at the root.

The root is therefore an isolated simple point of a rational polynomial
system. Its coefficients are algebraic: the Jacobian-nonzero zero locus
is zero dimensional with finitely many points, and elimination gives an
algebraic equation for each coordinate. This is an exact isolated-root
representation over Qbar, not a claim to display minimal polynomials or
radical expressions.

The verifier simultaneously proves, throughout D,

    2*10^-21 < Re F_20 < 3*10^-21,
    1*10^-21 < Im F_20 < 2*10^-21.               (5)

The tighter accepted enclosure starts

    F_20 = 2.6650234345855777515385842*10^-21
             +1.20476545551986511161340628*10^-21 i,

with real and imaginary enclosure errors less than10^-46. Hence
M20=61!!*F_20 is nonzero. Equation(1) is odd, so Gaussian sign symmetry
annuls every odd moment. Equation(3) therefore proves M1=...=M19=0,
and the first nonzero pure moment is exactly twenty.

## The cubic leading form is smooth

Let F3 denote only the cubic part of(1). Form the15-by18 matrix

    (quadratics)^3 -> quartics,
    (A,B,C) -> A*(F3)_Z+B*(F3)_W+C*(F3)_T.        (6)

The verifier removes the linear terms before forming this matrix. A
fifteen-column determinant selected in the certificate is enclosed,
uniformly on D, near

    0.00363565941890542427 + 0.00092274986146764662 i,

with positive absolute lower bound. Thus(6) is onto. If the three partials
had a common nonzero projective zero, every quartic would vanish there,
including Z^4,W^4,T^4, a contradiction. This proves smoothness directly.
It is not inferred from an approximate discriminant or a generic chart.

The slice also has (2*partial_Z partial_W+partial_T^2)F3=2T identically.
This motivated the search, but neither that classification-based exclusion
nor the smooth188 analytic proof is needed to accept the explicit lower
witness. They enter only the separate class discussion and upper bounds.

## Reproduction and literature boundary

Run the acceptance gate from the repository root:

```sh
.venv/bin/python research/jacobian-followup/moment-vanishing-audit/verify_odd_prefix_witness.py
```

The bounded [exploration](odd_prefix_witness_exploration.json) lifted the
retained homogeneous coefficient center by adding three linear parameters.
It found a numerical proposal on its first start. The
[builder](build_odd_prefix_witness.py) refined this proposal at384-bit
precision and stored rational data. The acceptance gate uses only those
rational data, not numerical residual tolerances or the search outcome.
Independent integer regressions compare the mixed-degree Wick evaluation
with direct symbolic powers through order six.

The interval semantics reuse the audited
[Acb documentation](https://python-flint.readthedocs.io/en/latest/acb.html)
and retained0.9.0 docstring snapshot. The installed core is FLINT3.6.0.

The current literature pass also read Bedratyuk,
[2610.13229v1,9 October2026](https://arxiv.org/html/2610.13229v1),
including its correspondence and nonfinite-generation arguments.
It concerns projective differential image moments in an infinite formal
moment algebra. Its SL(3) cubic invariants do not supply orthogonal
Gaussian pure-prefix certificates. No theorem from it is imported here;
its scope does not contradict finite generation in our fixed coefficient space.

Both lower-bound coefficient-witness obligations, homogeneous13-zero and
odd19-zero, are now complete as isolated algebraic representations. Their
first nonzero orders are14 and20. Least thresholds, twenty-moment
sufficiency, expanded global radical cofactors, independent review,
exceptional coefficient schemes, global term-count minima and lower-power
GVC remain unfinished. No full coefficient-ring radical search was run.
