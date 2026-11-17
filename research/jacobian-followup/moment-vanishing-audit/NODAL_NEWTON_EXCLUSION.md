# Homogeneous ternary cubics: the nodal Newton exclusion

Current quantitative refinement, 15 November 2026: the [exact Euler representative proof](EULER_REPRESENTATIVES_GLOBAL_297.md) gives global pure297, odd296, conic/full-W-linear258 and all-homogeneous258. Ordinary-node complements are297/271/244, with homogeneous nonisotropic-node101; cusp255/tangent284 are unchanged. Current intervals are global[20,297], odd[20,296] and all-homogeneous[14,258]. Earlier353/336 bounds below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I297 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [joint-degree actual-seed proof](JOINT_DEGREE_GLOBAL_353.md) gives global pure353, odd352, full W-linear336, rank-one140, mixed-parent255 and cusp255/tangent284. Current least-prefix intervals are global[20,353], odd[20,352] and all-homogeneous[14,336]. Earlier590/888 bounds and rank-one425/mixed-parent320 estimates below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I353 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [balanced logarithmic seed proof](BALANCED_LOGARITHMIC_GLOBAL_590.md) gives global pure590, all-homogeneous336, conic336 and full W-linear425. Current intervals are global/odd[20,590] and all-homogeneous[14,336]; smooth[20,94] and homogeneous-smooth[14,56] are unchanged. Earlier888/1855 bounds and frame-plus888 statements below retain their historical scopes. Leastness,20/14 sufficiency, expanded I590 cofactors and independent review remain open.

Current global refinement,15 November2026: the [rescaled singular seed proof](RESCALED_SINGULAR_GLOBAL_888.md)
now proves888 sufficient as a global pure-only prefix, without an extra
frame condition. Current intervals are global[20,888], odd[20,888] and
all-homogeneous[14,888]. Earlier1855 bounds and reservations about888
alone below are historical. The new proof remains internal and
independently unreviewed; leastness and expanded global cofactors are open.

Current lower-bound refinement,15 November2026: the [odd-invariant proof](ODD_CUBIC_CLASSIFICATION_AND_PREFIX_LOWER_BOUND.md)
raises the least global pure-prefix interval to[20,1855] and proves
homogeneous lower14, by algebraic existence without displayed coefficients.
Earlier lower10 statements remain valid explicit-witness evidence.

14 November 2026. Internal written proof importing Sabbah's affine Newton
filtration theorem. The checker verifies the algebra, not that theorem or
the analytic argument. No external review, formal verification or priority
is claimed.

## Result and complete homogeneous classification

Every homogeneous irreducible nodal ternary cubic has a nonzero positive
Gaussian pure moment. Together with the preceding smooth, cuspidal and
reducible results, this closes GMC for **every homogeneous ternary cubic**.
It does not close GMC for arbitrary nonhomogeneous ternary cubics.

More precisely, for a nonzero homogeneous cubic P in three standard Gaussian
coordinates, these conditions are equivalent:

1. Every positive pure moment E(P^m) vanishes.
2. P admits a factorization P=(l^T*x)*(x^T*A*x), with A symmetric,
   l nonzero, l^T*l=0, and A*l proportional to l.

For every such P and every polynomial R,

    E(R*P^m)=0 whenever m>deg R.

The zero polynomial satisfies the same statement with A=0. This mixed
cutoff is sharp by the [reducible theorem](REDUCIBLE_CUBIC_CLASSIFICATION.md).
The existing sharp pure gate10 applies to reducibles and cusps; sufficient188
applies to smooth cubics. This note supplies no explicit uniform numerical
pure gate for the nodal class or for all homogeneous ternary cubics.

## The imported theorem and its growth consequence

Use the affine-polynomial case of [Sabbah, *Hypergeometric periods for a
tame polynomial*, Section12, Lemma12.2, printed pages32--33](https://perso.pages.math.cnrs.fr/users/claude.sabbah/articles/sabbah_bbases.pdf).
For a convenient polynomial g nondegenerate with respect to its Newton
polyhedron at infinity, its exponential Gauss--Manin module G is regular
singular at tau=0. If the outer facet forms L_sigma have value1 on their
facets, the Newton degree of u*dx is

    delta*(u)=max_sigma,max_(nu in supp u) L_sigma(nu+1).

Lemma12.2 identifies N_alpha G with the increasing V_alpha G. In
particular, a polynomial top form of Newton degree<=alpha lies in V_alpha G.
Section1, equations(1.1)--(1.2) and Lemma1.3 fix the convention:
tau*partial_tau+beta is nilpotent on gr_beta^V G.

Consequently, for a top form omega in V_alpha G, every horizontal
functional h satisfies, on a closed smaller sector at tau=0,

    |h(omega)|=O(|tau|^(-alpha)*(1+|log tau|)^K)       (1)

for some finite K. To check the sign and the inequality, use the convergent
regular-singular decomposition (1.2): a graded piece of index beta evaluates
as tau^(-beta) times a polynomial in log tau; multiplication by tau lowers
the index by one. V_alpha contains only indices<=alpha. Analytic coefficients
are bounded on a smaller sector. There is no stationary-phase or Morse
assumption here. We use polynomial g on affine three-space, not the Laurent
torus variant whose convenience hypothesis is different.

## Formal moments supply the functional that must be excluded

Let Q(x)=x^T*x/2 and let P be homogeneous cubic. Over K=C(t), form the
top twisted de Rham quotient

    H=K[x]/sum_i (partial_i+t*P_i-x_i)K[x],
    nabla_t[u]=[partial_t u+P*u].                     (2)

Gaussian Wick integration by parts defines a horizontal K-linear map

    J:H -> C((t)),  J(u)=E(u*exp(tP)).                (3)

Laurent expansion makes rational coefficients legitimate. No convergence
of the original exponential series is assumed. If all positive pure
moments vanish, J([1])=1 and J(nabla_t^j[1])=0 for every j>=1. On the
differential submodule S generated by [1], these values are rational;
thus J|S is a nonzero morphism of differential modules to the trivial module.

At an ordinary nonzero t-point, extend this functional on the fiber of S
to the fiber of H and solve the dual linear ODE. Its horizontal continuation
restricts to J|S by uniqueness. It therefore evaluates [1] as the exact
constant1. This step only needs a finite-dimensional connection, supplied
below by the Newton-tame module after a rational change of variables. It
does not require the seed to generate H, a splitting over K, or convergence
of (3). A rapid-decay period representation is optional, not an extra
assumption in this proof.

## Normal form and the Gaussian metric

Let the singular point be the T-axis. An ordinary node has a nondegenerate
binary tangent cone. A linear choice of u,v makes this cone uv, so

    P=Tuv+a*u^3+b*u^2*v+c*u*v^2+d*v^3.

Replace T by T-bu-cv to obtain

    P=Tuv+a*u^3+d*v^3,  a*d!=0.                      (4)

If either a or d were zero the cubic would factor. Every transformation
here is invertible and linear; retain the transformed Gaussian quadratic

    Q(z)=z^T*M*z/2,  M nonsingular, z=(u,v,T).

The node is nonisotropic precisely when M_TT!=0: its original vector is
the image of the T-axis, and its squared Gaussian length is M_TT. The shear
in (4) does not change that axis. After this fixed coordinate change, the
operators in (2) have Gaussian term -Q_i in place of -x_i. The constant
volume determinant only changes the nonzero normalization of the functional;
the parameter-dependent determinant below must be retained. Put f=P-Q and g=-f.

## All Newton faces, including the boundaries

For a nonisotropic node the outer vertices of the Newton polyhedron of f
are

    A=(3,0,0), B=(0,3,0), C=(1,1,1), D=(0,0,2),

with the origin as the remaining vertex. Their coefficients are respectively
a,d,1,-M_TT/2, all nonzero. The three outer facets and normalized forms are

| Facet | L_sigma(u exponent,v exponent,T exponent) |
|---|---|
| ABC | (1/3,1/3,1/3) |
| ACD | (1/3,1/6,1/2) |
| BCD | (1/6,1/3,1/2) |

Every other quadratic monomial is strictly below all outer facets. In
particular, the shear in (4) only introduces interior quadratic terms.
The polynomial meets all coordinate axes and is convenient. Its outer
facet polynomials, up to nonzero signs, are

    a*u^3+d*v^3+Tuv,
    a*u^3+Tuv-(M_TT/2)*T^2,
    d*v^3+Tuv-(M_TT/2)*T^2.

Their respective partial derivatives in T,v,u are the monomials uv,Tu,Tv,
which never vanish on the algebraic torus. Each proper outer face also
has a partial derivative consisting of a single nonzero monomial; the
[checker](verify_nodal_newton.py) enumerates all13 faces and records the
witnesses. Hence g is Newton-nondegenerate at infinity for every allowed
metric and coefficients. No genericity or omitted determinant chart is
being imposed. Absorbing b,c before this check matters: leaving arbitrary
mixed binary cubic coefficients would invalidate this face argument.

All three facet weights sum to1. Therefore the constant volume form

    omega=du*dv*dT

has Newton degree1, so omega belongs to V_1 G. Equation(1) bounds every
horizontal evaluation of omega by O(|tau|^(-1)*(1+|log tau|)^K).

## The incompatible exponent

Make the exact change z=t*x and lambda=t^(-2). Homogeneity gives

    t*P(x)-Q(x)=lambda*(P(z)-Q(z))=lambda*f(z),
    dx=t^(-3)*dz.                                   (5)

Thus (2) is isomorphic to the pullback of the exponential module for
exp(lambda*f)=exp(-lambda*g). The top-form seed maps to t^(-3)*omega.
This is an isomorphism of connections, including the Jacobian, not only
an equality of exponents. For example its seed compatibility is the
twisted-divergence identity

    mapped nabla_t[1]-nabla_t[t^(-3)*omega]
      =t^(-4)*sum_i (partial_i+lambda*f_i)(z_i)*omega=0.

If all pure moments vanished, the functional above would give

    1=t^(-3)*h(omega), hence h(omega)=t^3=lambda^(-3/2). (6)

Choose a branch on a small sector. Horizontal continuation reaches lambda=0
on that sector; the Newton module has no singularities at nonzero lambda.
But lambda^(-3/2) violates (1) with alpha=1, for every finite K. This proves
pure nonvanishing for every nonisotropic-node irreducible homogeneous cubic.

The strict inequality is essential: a nonsingular Gaussian quadratic has
volume Newton degree3/2, and its lambda^(-3/2) period is perfectly compatible
with (1). If the node is isotropic, M_TT=0 and the T-axis intercept disappears;
the present Newton argument is not applied to that boundary.

## Closing the remaining homogeneous cases

For an isotropic node, the [Gaussian frame criterion](RANK_ONE_WT2_CLASSIFICATION.md)
puts P=W*H2(Z,T)+G3(Z,T). The ordinary-node tangent cone H2 is nondegenerate.
The existing [nondegenerate conic exclusion](QUADRATIC_W_LOG_PERIODS.md)
therefore gives a nonzero pure moment on this entire boundary. Together
with (6), this exhausts irreducible nodal cubics.

An irreducible singular plane cubic over C has either an ordinary node or
a cusp: its singular multiplicity is two (multiplicity three would leave
a binary cubic, which factors); a rank-two tangent cone gives the node,
and a rank-one tangent cone with nonzero cubic on its tangent line gives
the cusp. If that latter coefficient vanished there would be a line factor.
The [cuspidal theorem](CUSPIDAL_CUBIC_SHARP_THRESHOLD.md) excludes the full
cuspidal class. [Smooth cubics](COERCIVE_GRADIENT_MOMENTS.md) are already
excluded, and [reducibles](REDUCIBLE_CUBIC_CLASSIFICATION.md) have the stated
isotropic-eigenvector factor classification and sharp mixed cutoff. This
proves the complete homogeneous assertion without an unclassified orbit.

## Reproduction and surviving obligations

Run `.venv/bin/python research/jacobian-followup/moment-vanishing-audit/verify_nodal_newton.py`.
The [receipt](nodal_newton_verification.json) checks the universal normal
form, all facets and13 outer-face witnesses, all six quadratic supports,
the exact exponent and connection rescalings, the isotropic negative
control and Gaussian exponent calibration. The [source audit](nodal_newton_source_audit.json)
binds the cached Sabbah PDF and specifies the sections actually read.
The imported theorem and regular-singular argument remain written mathematics.

The subsequent [Gaussian-phase critical profiles](NODAL_GAUSSIAN_CRITICAL_PROFILES.md)
now account for all seven critical points in the irreducible nonhomogeneous
nodal case. Every escape has a nonzero action, while the center-line
boundary has bounded critical points. The subsequent
[center-line cut-decay theorem](NODAL_CENTER_LINE_DECAY.md) now proves
normalized center-line period decay and excludes the zero-critical-line
and closed-half-plane critical-value classes, including Gaussian means.
The later [dominance theorem](NODAL_MORSE_DOMINANCE.md) closes the
strict triangle and the finite-critical Morse/distinct-value locus with
nonzero quartic gamma. Only finite degeneracy/value-collision/resonance
strata remain in this irreducible nodal chart.

The previous [nodal conic reduction](NODAL_CUBIC_REDUCTION.md) remains useful
for nonhomogeneous cubics. This result does not establish its general conic
period decay estimate. The equality (5) fails with arbitrary lower-degree
perturbations: after rescaling, their contributions depend on lambda in
additional powers. They cannot be treated as coefficients of the fixed f
used here. General nonhomogeneous singular-leading cubics, a numerical
uniform nodal pure gate, full coefficient schemes in dimensions3,4, exact
generic quaternary null torsion and lower-power GVC remain open.
