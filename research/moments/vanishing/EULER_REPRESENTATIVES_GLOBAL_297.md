# Exact Euler representatives: global297 and homogeneous258

15 November 2026. Internal written refinement of the
[joint-degree353 proof](JOINT_DEGREE_GLOBAL_353.md), using the same original
Gaussian Wick functional, twisted quotient and analytic exclusions.

| Coefficient scope | Sufficient pure gate or exclusion | Least-prefix interval |
|---|---:|---:|
| Every complex ternary polynomial of total degree at most three | **297** | **[20,297]** |
| Odd ternary cubics | **296** | **[20,296]** |
| All homogeneous ternary cubics | **258** | **[14,258]** |
| All sixteen W-linear cubic coefficients | **258** | Not computed |
| Nondegenerate quadratic W coefficient, cubic G | **258**, an exclusion | Not computed |

The global test uses only the initial pure moments. No separate Gaussian
frame condition is imposed. Pure nullity then gives the retained Gaussian
W-linear frame and sharp mixed cutoff m>3*deg Q. The smooth-leading and
homogeneous-smooth intervals stay[20,94] and[14,56]. Both isolated algebraic
lower witnesses with first nonzero moments20 and14 remain complete.
No leastness or20/14 sufficiency is proved.

The [checker](verify_euler_seed_bounds.py),
[certificate](euler_seed_bounds_certificate.json),
[receipt](euler_seed_bounds_verification.json) and
[source audit](euler_seed_bounds_source_audit.json) bind exact auxiliary algebra.
Full coefficient coverage, rank completeness, seed/cofactor degree inductions
and analytic nonzero zeroth coefficients remain written mathematics.
Independent replay, external review, formal verification and expanded global
radical cofactors remain outstanding.

## 1. Change the representative of the connection, keeping the same connection

For the ordinary-node binary phase use the retained normalization

    I=constant*t^-2*D^-1/2*integral exp(W/delta) dx dy,
    delta=t^2*D, deg_t D<=1, D!=0 in C[t],
    W=x^2*y^2/2+A*x^3+B*y^3+lower amplitude terms,
    deg_(t,x,y) W<=4.

The previous principal coefficient opens and full node/frame case split
remain in force. In particular every member with any isotropic zero-F2
ordinary node has already been assigned to the conic frame branch.
The current calculation does not replace that full case split by a check
at one chosen node.

Put E=x*partial_x+y*partial_y and G_i=delta*partial_i+W_i. The retained
normalized source is

    C=S/delta^2,
    S=delta*W_t-delta_prime*W-2*t^3*D^2-t^4*D*D_prime/2.

For every amplitude f,

    G_x(x*f)+G_y(y*f)=delta*(E(f)+2f)+E(W)*f.

Therefore adding delta_prime/(4*delta^2) times this exact twisted divergence
to the old connection does not change its class. A new representative is

    nabla_E(f)=f_t+(delta_prime/(4*delta))*E(f)+(S_E/delta^2)*f,
    S_E=delta*W_t-delta_prime*(W-E(W)/4)-t^3*D^2.       (1)

This is equality of operators on the same quotient, not a new period or a
change of the normalization. In particular the volume divergence2f is
required. The affine Gaussian means, mixed terms and lower coefficients
stay arbitrary.

The only amplitude-four term of W is x^2*y^2/2, whose coefficient is
constant in t. It is killed both by W_t and W-E(W)/4. Thus

    amplitude degree S_E<=3,
    joint t/x/y degree S_E<=6.                         (2)

It is not sufficient to replace S by S_E alone: the Euler derivative in
(1) must act on every later seed.

## 2. Actual ordinary-node seeds and smaller complement gates

Write the actual derivative seed as R_j/delta^(2j), with R_0=1. Formula(1)
gives the complete recurrence

    R_(j+1)=delta^2*R_j_prime-2j*delta*delta_prime*R_j
              +(delta*delta_prime/4)*E(R_j)+S_E*R_j.  (3)

Its difference from the old recurrence applied to the same R_j is the
exact divergence with witnesses delta_prime*x*R_j/4 and
delta_prime*y*R_j/4. This proves by induction that the original moments,
not merely a new family of normalized periods, are being differentiated.
From(2)--(3),

    joint degree R_j<=6j, amplitude degree R_j<=3j.

Use the unchanged complete principal projections and bases. Their whole
cleared corrections have joint costs k=2,3,4 for nodal, transverse and
triangle types. Every non-basis correction strictly lowers ordinary
amplitude degree; equal-degree remainders are already basis terms.
The required depths are now3j-2,3j-2,3j-1. Therefore

    H_j=delta^(2j)*A0^r_j,
    deg_t N_(i,j)<=b_j-d_i,
    b_j=12j-4 (nodal),15j-6 (transverse),18j-4 (triangle).

The d_i are the ordinary basis grades. Every higher scalar coefficient
contains H_j, with valuation at least4j. The retained minimal-seed
cofactor argument, including its analytic nonzero zeroth coefficient,
gives the complete complement exclusions

| Type | Rank | Zeroth degree bound | Gate |
|---|---:|---|---:|
| Irreducible nodal cubic | 7 | 12*28-4*7-11 | **297** |
| Transverse line-smooth-conic | 6 | 15*21-6*6-8 | **271** |
| Three independent lines | 5 | 18*15-4*5-6 | **244** |

All smaller seed ranks and row subsets containing the unit row give smaller
bounds. These are complement bounds until section5 adds back the frames.
No generic seed cyclicity or Morse specialization is assumed.

For a homogeneous irreducible nodal cubic with nonisotropic node, the
retained linear normalization makes D=m constant and W independent of t.
The Euler cancellation now gives joint degree S_E<=4, instead of5.
Constant principal divisors preserve joint degree in the reductions, so
deg_t N_(i,j)<=4j-d_i and the complement gate improves129 to

    4*28-11=101.                                     (4)

Zero Gaussian means are essential to this homogeneous specialization.

## 3. A conic endpoint is a cube of a low-degree polynomial

Retain the balanced logarithmic conic phase and its original seed:

    I=constant*t^-1*integral exp(W/t^2) dy/y,
    W=sum_(k=-3)^3 A_k(t)*y^k, deg_t A_k<=3-k,
    S=t*W_t-2W-t^2.

The old proof covers all actual endpoint orders0<=p,q<=3, including
regular ends and the trapped g=0 formal component. When q<3 its existing
all-rank bounds have maximum145. Only the q=3 cases need improvement.

In every q=3 conic chart,

    A_-3=c*K(t)^3, c!=0 constant in t, deg_t K<=2.    (5)

Here K is nonzero as a polynomial, and it may vanish at t=0. To prove(5),
in the f!=0 chart recall eta=f*(rho-sigma) and

    U=((rho*c-b)*t+1)/eta,
    V=((b-sigma*c)*t-1)/eta,
    K=U*V-a*t^2/f.

The negative Laurent coefficients of the two rescaled coordinates are
-f*K/eta and rho*f*K/eta. Only the homogeneous degree-three part of G
can supply y^-3; evaluating it at this constant direction proves(5).
The Gaussian term has order at worst-2. For f=0,e!=0, the negative
coefficient of Zhat is zero and that of That is

    K=(t*z0-t^2*A(z0))/e, z0=-c/e.

Thus A_-3=[T^3]G*K^3. The trapped component has no negative coordinate
order and cannot occur in q=3. If a candidate K vanishes identically, the
actual negative order is below3 and the old endpoint cases apply.

Set

    B=t*K_prime-2*K/3,
    kappa=B/(t*K), h=t^3*K.

For the logarithmic measure the exact operator is

    Dlog(f)=y*f_y+(y*W_y/t^2)*f.

There is no additional volume divergence constant in this formula. Adding
kappa*Dlog to the old connection gives the same quotient connection, with

    nabla_K(f)=f_t+kappa*y*f_y+(S_K/h)*f,
    S_K=K*S+B*y*W_y.                                (6)

Using(5), its y^-3 coefficient cancels exactly. Hence S_K has Laurent
support[-2,p] and coefficient degrees deg_t(S_K)_k<=5-k. Clearing K
rather than A_-3=constant*K^3 is the denominator saving. Cancelling the
endpoint with an unrestricted degree-six leading coefficient would lose
that saving.

## 4. Keep the Euler derivative and separate endpoint denominators

The actual logarithmic seeds are R_j/h^j, with

    R_0=1,
    R_(j+1)=h*R_j_prime-j*h_prime*R_j
               +t^2*B*y*(R_j)_y+S_K*R_j.           (7)

Induction gives support[-2j,pj] and deg_t(R_j)_k<=5j-k. The Euler derivative
in(7) is essential to retaining the original normalized moment series.
No new singular t=0 Gaussian integral is defined.

For p>0 the unchanged full logarithmic basis is y^-3,...,y^(p-1).
High and low reductions do not feed outside each other's already retained
basis boundary. With the separate endpoint divisors, sufficient depths are

    high_j=p*(j-1)+1,
    low_j=max(0,2j-3).

For p=0 the basis is1,y^-1,y^-2, high_j=0 and
low_j=max(0,2j-2). The unit must not be omitted at this regular positive end.
Using Aplus=A_p for p>0 and Aminus=A_-3,

    H_j=h^j*Aplus^high_j*Aminus^low_j,
    deg_t N_(k,j)<=b_j-k,
    b_j=5j+(3-p)*high_j+6*low_j.

The polynomial H_j has valuation at least3j, even if K(0)=0. Thus there
is no extra derivative-order loss in the scalar zeroth-cofactor comparison.
The largest bound over every actual seed rank and row subset containing0 is
at full quotient rank:

| Positive order p, with q=3 | 0 | 1 | 2 | 3 |
|---|---:|---:|---:|---:|
| Sufficient exclusion | 69 | 130 | 201 | **258** |

For p=3, the seed degree bounds are5,16,33,50,67,84, summing to255;
the row exponents -3,-2,-1,0,1,2 sum to-3, yielding258. Together with
all retained q<3 bounds<=145, this gives conic exclusion258 for cubic G.
It is not a numeric bound for arbitrary-degree G.

The retained mixed-parent255 and rank-one140 bounds imply full f=0 gate258
and full sixteen-coefficient W-linear gate258. Their sharp mixed cutoff
3*deg Q is unchanged. For all homogeneous cubics, combine smooth56,
cusp10, reducible10, nonisotropic-node101 and isotropic-node conic258.
This gives all-homogeneous258; parity does not reduce an even258 to257.

## 5. Complete global cover and remaining obligations

| Leading type, arbitrary lower terms | Full-class sufficient gate/exclusion |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible nodal cubic | max(297,258)=297 |
| Irreducible cuspidal cubic | 255 |
| Transverse line-smooth-conic | max(271,258)=271 |
| Tangent line-smooth-conic | 284 |
| Three independent lines | max(244,258)=258 |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 99 |

The maximum297 proves the pure-only global test. Odd polynomials have all
odd pure moments zero automatically, so their first296 conditions include
every needed condition through297. No unrestricted296 claim follows.
In the rational twenty-coefficient ring,
radical(I297)=radical(I_infinity), without ordinary ideal equality or
expanded radical cofactors. At most ceil(log2(297-19))=9 exact monotone
search rounds against I297 would find the least global K; none has been run.

The checker supplies symbolic source/divergence identities, six full-Gaussian
ordinary-node connection models with three actual seeds each, the homogeneous
source cancellation, both conic factor charts and eight Laurent connection
models with three seeds each. Cases include D(0)=0, K(0)=0 with simple and
double zero, K constant, every q=3 positive endpoint and every retained q<3
case. The initial rational-constant denominator failure is preserved; clearing
polynomial t-denominators permits rational scalar coefficients.

Selected reading of [Crew--Feldbrugge2610.22767v2, Appendix B](https://arxiv.org/html/2610.22767v2)
rechecked the divergence correction and Euler relation, including the volume
weights. The appendix assumes tame Morse phases with distinct critical
values. Its generic good-basis and rank conclusions are not used to cover
our degenerate strata. The two representative identities above are checked
directly, and the prior complete analytic exclusions remain the inputs.

Ordinary-node297 and tangent284 are now the largest class gates. Sharper
source representatives or reductions may improve them, but exact leastness,
20/14 sufficiency, expanded global radical certificates, independent review,
exceptional coefficient schemes, term-count minima and lower-power GVC remain
open. The qualitative cubic case cover and both lower-witness constructions
stay closed internally.
