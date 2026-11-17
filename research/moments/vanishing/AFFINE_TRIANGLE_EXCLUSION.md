# Three independent leading lines: all affine perturbations are pure-excluded

Current quantitative status,15 November2026: the [global pure-prefix proof](GLOBAL_TERNARY_CUBIC_PURE_PREFIX.md)
now gives a sufficient1855 pure-only gate and least-prefix interval[10,1855].
Earlier absence-of-global-upper or unquantified-node statements below are
historical. Their earlier class/certificate scopes remain unchanged; the
new proof remains internal and independently unreviewed.

14 November 2026. Internal written proof with exact auxiliary algebra.
No independent specialist replay, external review or formal verification.

## Theorem and remaining boundary

**G3TRIANGLEAFFINE20261114.** If a ternary cubic's leading homogeneous
form is a product of three linearly independent linear forms, some
positive pure Gaussian moment is nonzero, for arbitrary lower terms.
The full transformed nondegenerate affine Gaussian is retained.

This closes the projective triangle class for GMC. It is stronger than
the earlier centered-product theorem because the lower terms need not
factor with the leading lines. Three concurrent distinct leading lines
are a different class and are not included. The remaining possible
ternary cubic GMC failures have a line-smooth-conic leading curve,
three concurrent leading lines, or a repeated leading line. All
homogeneous cubics were already safe. No numerical uniform affine
triangle prefix is asserted; coefficient schemes and lower-power GVC
remain separate.

The subsequent [concurrent-line theorem](CONCURRENT_LINE_EXCLUSION.md)
closes that distinct-line class as well. Remaining possible cubic GMC
failures have line-smooth-conic or repeated-line leading curves; the
triangle theorem below retains its independent-factor hypothesis.

## 1. Reduction and rank five

The leading triangle has three ordinary nodes. If a node is isotropic
and F2 vanishes there, the retained
[W-linear frame criterion](RANK_ONE_WT2_CLASSIFICATION.md) puts P in the
nondegenerate quadratic-W conic class: ordinary-node tangent rank is2.
The [nondegenerate conic theorem](QUADRATIC_W_LOG_PERIODS.md) pure-excludes
it. If F2 is nonzero at every node, the
[nodal coercive-gradient criterion](COERCIVE_GRADIENT_MOMENTS.md) applies.

Otherwise choose a nonisotropic node with F2=0. Linear coordinates put
F3=Tuv with that node on the T-axis. Translations of u,v remove uT,vT;
a T translation removes the quadratic uv term. Thus

    P=T*(uv+h)+b*u^2+c*v^2+eu+fv+k,
    Q=(A*u^2+B*v^2+m*T^2+2*Cuv+2*DuT+2*EvT)/2
        +l_u*u+l_v*v+l_T*T+q0.                    (1)

Keep det Hess Q!=0 and m!=0. The other two zero-F2 nodes cannot be
isotropic, since that alternative was already excluded. Hence

    b=0 implies A!=0; c=0 implies B!=0.           (2)

No assumption is imposed on A when b!=0 or on B when c!=0.

For large |t|, S_t=tP-Q has nonzero quadratic value at all three
leading nodes: -m/2,bt-A/2,ct-B/2. Its complex gradient is linearly
coercive. The three projective quadratic gradient equations have
exactly the three leading nodes at infinity, each of length1, by the
nonzero local Jacobian. Bezout leaves affine length8-3=5. The retained
tame comparison gives top twisted rank5. Exceptional finite t-values
are irrelevant over C(t) and are avoided by taking |t| large.

## 2. Exhaustive critical-set split

For h!=0, set v=-h/u. Finite P-critical points are the critical points
in C* of

    R(u)=b*u^2+eu+k-fh/u+c*h^2/u^2.               (3)

The normal coordinates uv+h and a T shear give the exact local type
alpha+s*w+y^(mu+1) at every isolated point. If R is nonconstant, let
S be its nonzero exponent support among2,1,-1,-2. The finite Milnor
rank is max S-min S, including zero when S is a singleton. Indeed uR'
is a Laurent polynomial with exactly those extreme nonzero coefficients.
Its nonzero-root multiplicities sum to that width. It is at most4.
If R is constant, b=c=e=f=0, and the critical set is the smooth conic
uv=-h,T=0 with value k.

For h=0, on v=0,u!=0 an isolated critical point exists exactly when
b*e!=0; it is Morse, with u=-e/(2b),T=-f/u. The other axis is symmetric.
The critical set is finite exactly when both (b,e) and (c,f) are nonzero
pairs and at least one of e,f is nonzero. Its rank is

    F=indicator(b*e!=0)+indicator(c*f!=0).         (4)

All non-finite cases are exactly the following two families, after a
permutation of the leading coordinates:

    critical conic: P=w*(xy+d*w+a)+k, a!=0;
    critical axes: P=Tuv+b*u^2+c*v^2+k.           (5)

For example b=e=0,h=0,f!=0 gives the first family with w=v,a=f,d=c;
the symmetric alternative is analogous. If e=f=h=0, the second family
applies. The [checker](verify_affine_triangle_profiles.py) exhausts all32
zero/nonzero support patterns in h,b,c,e,f:23 finite,5 conic and4 axes.
Coefficient degeneracies within each finite support are retained through
their full Milnor multiplicities, not sampled away.

## 3. All escapes on the finite-critical locus

Only nodes with F2=0 can carry escapes as t tends to infinity. Nodes
with F2!=0 already have a fixed-P linear gradient bound near their
direction. At a zero-F2 node use local height coordinates

    P=Z*(XY+H)+G(X,Y), deg G<=2,

with nonzero Gaussian Z^2 coefficient M. Write E1,F1 for the linear
coefficients of X,Y and C2 for the Y^2 coefficient after swapping X,Y
if needed. The finite-critical possibilities are exactly:

| Conditions | Escape count | Leading action after k*t |
|---|---:|---|
| H!=0 | 1 | H^2*t^2/(2M) |
| H=0,E1*F1!=0 | 3 | -3*M*K^2*t^(2/3)/2, K^3=E1*F1/M |
| H=F1=0,E1*C2!=0 | 4 | -M*K^2*t^(1/2), K^4=-2*C2*E1^2/M |

All leading coefficients are nonzero. These are the explicit formal
implicit/Gaussian rescalings of the
[nodal critical-profile note](NODAL_GAUSSIAN_CRITICAL_PROFILES.md), with
the omitted pure cubic coefficients now zero. They were not needed
in these three leading models. If the last quadratic term is also
absent, or both linear terms vanish, a critical conic/axis occurs,
already assigned to (5), not an omitted finite-critical escape.
The physical Gaussian fluctuation along the node axis is O(1);
all transverse fluctuations tend to zero. Prefactors are respectively
t^-2,t^-4/3,t^-5/4. No seed-cyclicity or action-distinctness assumption
is made.

For explicit rank exhaustion, when h!=0 the T-node supplies1. If b=0,
the u-node supplies1 when e!=0,3 when e=0,f!=0,4 when e=f=0,c!=0.
If b!=0 it supplies0. At the v-node interchange b/e with c/f.
Together with the width in (3) the sum is always5. When h=0 on the
finite locus, the T-node supplies3 if e*f!=0, otherwise4. The u-node
supplies indicator(b=0), the v-node indicator(c=0). With (4) the sum
is again5. Thus the displayed formal implicit branches and finite
Milnor characters exhaust all five dimensions.

The finite characters have scalar form exp(alpha*t)*t^(-1-j/(mu+1))
times power series in1/t,1<=j<=mu, by the explicit power-phase
recurrences in the [full formal-frame proof](NODAL_FORMAL_FRAME.md).
Every finite zero-exponential seed value has positive valuation.
Every escape has nonzero exponential factor, including resonant
opposite fourth roots.

## 4. Confluent numerator interpolation handles different escape scales

Here is a full-frame argument that avoids normalizing ordinary Lagrange
polynomials independently across very different node scales. Choose a
constant generic linear form L that separates finite points, has nonzero
derivative along each finite univariate direction, and separates every
leading escape coefficient of the same scale. Write finite nodes L_p,
with multiplicities mu_p, and escaping nodes L_i(t). Define

    B(X)=product_p (X-L_p)^mu_p * product_i (X-L_i(t)).

Its degree is5. Use the five polynomial amplitudes

    finite p,j: B(L)/(L-L_p)^(mu_p-j), 0<=j<mu_p;
    escape i:   B(L)/(L-L_i(t)).                   (6)

Every quotient in (6) is a polynomial, of degree at most4. They are
linearly independent by uniqueness of the partial fractions with
denominator B, hence a basis of degree<5 polynomials in L.

To check their formal period matrix, for a finite point p divide its
normalized character columns by B_p(L_p), B_p=B/(X-L_p)^mu_p. The
finite local jet matrix is triangular with nonzero diagonal, since
B_p(X)/B_p(L_p) has bounded Taylor coefficients near p, L is an invertible
univariate coordinate, and the volume/Gaussian amplitude is a unit.
At another finite point the rows have a zero of the prescribed Milnor
order, so the leading cross jets vanish.

For an escape column i divide by B'(L_i). A different simple escape row
j evaluated at L=L_i+delta has the exact factorization

    [B(L)/(L-L_j)]/B'(L_i)
      =delta/(L_i-L_j)
           *product_(other nodes l) (1+delta/(L_i-L_l)), (7)

with the finite-node factors raised to their multiplicities. The
factor for j is omitted. All separations from an escape tend to
infinity, and delta has bounded formal Gaussian moments. Thus the
escape block has a nonzero Gaussian diagonal and vanishing cross
terms, regardless of the relative rates t,t^1/3,t^1/4. A finite row
at escape i instead carries delta/(L_i-L_p)^(mu_p-j), so it too
vanishes after this column normalization. An escape row at finite p
contains (L-L_p)^mu_p divided by the large L_p-L_i and bounded relative
factors. Its leading finite jets vanish as well.

Consequently the full five-by-five period matrix is invertible over
a Puiseux Laurent field. Column normalizations here test its determinant;
they do not redefine the horizontal functionals or their seed values.
The original complete formal frame consists of the finite power
characters and the escaping Gaussian characters. After ramification
it is semisimple. The
[positive formal seed criterion](CUSPIDAL_AFFINE_EXCLUSION.md), Section2,
excludes a rational quotient taking the seed to1: all its zero-exponential
seed values have strictly positive valuation. This proves the entire
finite-critical case, including coincident values and all resonance.

## 5. Critical conics give a fixed Laurent block of rank four

Take the first family in (5). Conditions (2) at the x/y nodes give
nonzero Gaussian diagonal coefficients A_x,A_y. On xy=-a,w=0 the
Gaussian restriction is a Laurent phase

    psi(x)=-A_x*x^2/2-A_y*a^2/(2*x^2)
              -l_x*x+a*l_y/x+constant.            (8)

Use s=xy+a, y=(s-a)/x, s=rS,w=rW,r=t^-1/2. Exactly,

    S_t=k*t+SW+dW^2+psi(x)+O(r),
    dx dy dw=r^2*(dx/x) dS dW.                   (9)

The perturbation is a polynomial in r with Laurent-polynomial x
coefficients. The normal Gaussian pair has determinant -1 for every d.
The fixed Laurent quotient has rank4. For log measure dx/x its recurrence
on x^j is

    j*x^j-A_x*x^(j+2)-l_x*x^(j+1)
         -a*l_y*x^(j-1)+A_y*a^2*x^(j-2)=0.       (10)

Four consecutive moments, for example powers0,1,2,3, determine every
other moment: use A_x!=0 toward positive powers and A_y*a^2!=0 toward
negative powers. They are four free initial values, so these powers
form a basis. Tensoring with the Gaussian pair keeps rank4.
The coefficientwise fixed-phase conjugacy ell(exp(F_r-F0)*a) gives
four formal characters, all with scalar seed exp(k*t)*t^-1*C[[r]].

If d!=0 the fifth character is the isolated Morse point
(x,y,w)=(0,0,-a/(2d)), with value k-a^2/(4d) and Hessian determinant
-a^2/(2d). Its seed prefactor is t^-3/2. Polynomial cutoffs in w
separate this point from the conic: eighth powers of1-w/w_p and w/w_p
make cross terms vanish. Together with1,x,x^2,x^3 on the conic they
prove a full rank5 frame.

If d=0 there is instead one escape w~(a/Q_ww)*t, whose action has
nonzero quadratic coefficient a^2/(2Q_ww). Q_ww!=0 follows from the
already excluded isotropic zero-F2 alternative. Use the same cutoffs
with the formal escaping w_i(t) in place of w_p. Their cross terms
again vanish. The conic rank4 and escape rank1 exhaust5.

Fixed Laurent characters are formally trivial after the common
exp(k*t) and the volume power are removed; no Morse hypothesis on
psi is needed. Their zero-exponential seed values and any finite
Morse value have positive valuation. The formal seed criterion proves
the critical-conic exclusion.

## 6. Critical axes use uniform normalized cut decay

For the second family in (5), all P-critical values equal k. There
is a uniform complex gradient bound for S_t outside one fixed ball,
for every sufficiently large |t|. Here is the missing uniformity check.
If ||x||->infinity and grad S_t=o(||x||), its limiting direction is
one of the three leading nodes. Near the T-node, divide the u,v
equations by tT. Their coefficient matrix tends uniformly to the
invertible off-diagonal unit matrix, so u,v=O(1/t). The T equation
then has leading term -mT, a contradiction. Near the u-node the v,T
equations similarly give v,T=O(1/t). The u equation then has leading
(2bt-A)u, whose coefficient has a uniform nonzero lower bound for
large |t| by (2). The v-node is symmetric. The same estimates exclude
bounded-t subsequences in that large-|t| region. This proves

    ||grad S_t(x)||>=c0*||x|| for ||x||>=R0,
    c0,R0 independent of large t.                 (11)

All phase-critical branches are bounded. For epsilon=1/t, stationarity
of P-epsilon Q gives derivative -Q along each critical branch. Therefore
their actions are k*t+O(1), including all collisions.

Apply the full mechanism of the
[uniform cut/current-limit proof](NODAL_CENTER_LINE_DECAY.md), now with
rank5 rather than7: (11) gives radius-squared carriers <=C0+C1*s;
compatible semialgebraic cut pairs give a flat basis with fixed finite
multiplicities; uniform mass and boundary-mass bounds control each
finite cut band; the normalized tails are uniformly exponentially
small. Every fixed-band limit lies in the complex fiber P=k, where
the holomorphic top volume vanishes. Hence every normalized cut period
tends to zero. These are exactly the hypotheses used in that proof;
no isolated critical point or isolated critical line is required here.

Choose a t-ray with Re(k*t)<=0. All five actual scalar basis periods
then tend to zero. The retained rational seed/actual-period bridge
would represent all pure zeros by a nonzero constant scalar period,
which is impossible. This handles intersecting critical axes and their
possible logarithmic formal behavior without a semisimplicity assertion.

## Evidence and scope

The [receipt](affine_triangle_profiles_verification.json) binds exact
checks of the32-pattern partition, all23 finite rank counts, infinity
Jacobians, local escape actions/Hessians, the12 confluent numerator
partition cases, factorization (7), conic identities and Laurent
recurrence. The [source audit](affine_triangle_source_audit.json) separates
those checks from the written formal completeness and uniform analytic
argument. The first checker stopped at a presentation-sensitive symbolic
x^j simplification; its exact input and terminal failure are retained.
The corrected universal coefficient identity and integer-power checks pass.

Every affine triangle-leading cubic is now pure-excluded. Concurrent
leading lines, line-smooth-conic curves and repeated-line degenerations
still need their own complete analysis; neither a triangle coordinate
change nor the centered homogeneous theorem covers them automatically.
