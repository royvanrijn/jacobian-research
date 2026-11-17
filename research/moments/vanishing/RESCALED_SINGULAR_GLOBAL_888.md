# Rescaled singular seeds give a global888 pure-moment gate

Current quantitative refinement, 15 November 2026: the [exact Euler representative proof](EULER_REPRESENTATIVES_GLOBAL_297.md) gives global pure297, odd296, conic/full-W-linear258 and all-homogeneous258. Ordinary-node complements are297/271/244, with homogeneous nonisotropic-node101; cusp255/tangent284 are unchanged. Current intervals are global[20,297], odd[20,296] and all-homogeneous[14,258]. Earlier353/336 bounds below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I297 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [joint-degree actual-seed proof](JOINT_DEGREE_GLOBAL_353.md) gives global pure353, odd352, full W-linear336, rank-one140, mixed-parent255 and cusp255/tangent284. Current least-prefix intervals are global[20,353], odd[20,352] and all-homogeneous[14,336]. Earlier590/888 bounds and rank-one425/mixed-parent320 estimates below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I353 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [balanced logarithmic seed proof](BALANCED_LOGARITHMIC_GLOBAL_590.md) gives global pure590, all-homogeneous336, conic336 and full W-linear425. Current intervals are global/odd[20,590] and all-homogeneous[14,336]; smooth[20,94] and homogeneous-smooth[14,56] are unchanged. Earlier888/1855 bounds and frame-plus888 statements below retain their historical scopes. Leastness,20/14 sufficiency, expanded I590 cofactors and independent review remain open.

15 November 2026. Internal written refinement of the retained complete
ternary cubic case cover. The new conclusion is

    deg P<=3, P in C[X,Y,T],
    E(P^j)=0 for j=1,...,888
      iff every positive pure moment vanishes.    (1)

No separate Gaussian-frame hypothesis is needed. The previously retained
frame-plus888 theorem now becomes a pure-only test on the whole coefficient
space. The least global prefix is in **[20,888]**, improving[20,1855].
The [odd20 witness](ODD_CUBIC_TWENTIETH_MOMENT_WITNESS.md) supplies the lower
bound. The odd-only and all-homogeneous intervals become[20,888] and[14,888];
the smooth intervals remain[20,94] and homogeneous-smooth[14,56].

The first888 zero moments also force the Gaussian W-linear frame and the
retained sharp mixed cutoff m>3*deg Q. These conclusions use the retained
[complete cubic theorem](TERNARY_CUBIC_GMC_CLOSURE.md). Optimality of888,
twenty-moment sufficiency, ordinary ideal generation and expanded global
radical cofactors are not claimed.

The [checker](verify_rescaled_singular_prefixes.py),
[certificate](rescaled_singular_prefixes_certificate.json),
[receipt](rescaled_singular_prefixes_verification.json), and
[source audit](rescaled_singular_prefixes_source_audit.json) bind the auxiliary
identities and ten model calculations. They do not replace the written
all-coefficient rank, seed and analytic arguments. Independent specialist
replay, external review and formal verification remain outstanding.

## 1. Rescale the binary variables and retain the volume Jacobian

The [ordinary-node normal forms](GLOBAL_TERNARY_CUBIC_PURE_PREFIX.md)
complete the quadratic T variable and remove the mixed binary cubic terms.
First remove every member having any isotropic ordinary node with F2=0;
that branch uses the retained888 gate. On its complement the chosen D
and all required axis coefficients are nonzero by the retained full node
case split. Nonzero D at just one selected node is not a substitute for
that split. The remaining normal forms give a generic-t normalized seed

    I(t)=constant*D^(-1/2)*integral exp(V(t,u,v)/delta) du dv,
    D=m-2alpha*t !=0, delta=t^2*D.                (2)

Set x=t*u,y=t*v and W(t,x,y)=V(t,x/t,y/t). The exact transformed seed is

    I(t)=constant*t^(-2)*D^(-1/2)
              *integral exp(W(t,x,y)/delta) dx dy. (3)

The rescaled phase is a polynomial satisfying total degree

    deg_(t,x,y) W<=4.                            (4)

To see this before the ordinary-node translation, write the completed
phase numerator as

    U=D*(t*(G3+g2)-q2)+(t*(H2+ell)-J)^2/2.

Then t^2*U(t,x/t,y/t) is polynomial of total degree<=4: its square term
is the square of H2(x,y)+t*ell1+t^2*ell0-J1-t*J0. Its other terms have
the same bound. The ordinary-node translation becomes a constant x,y
translation, preserving(4). All affine Gaussian means, mixed terms and
lower coefficients are retained. The same argument applies directly to
the cusp/tangent phases in the [weighted proof](SQUAREFREE_SINGULAR_FINITE_PREFIXES.md).

The ordinary-node phase has the same nine monomials as before, but their
coefficient degree bounds are now

| Monomial | x^2*y^2 | x^3,y^3 | x^2,xy,y^2 | x,y | 1 |
|---|---:|---:|---:|---:|---:|
| Maximum t-degree | 0 | 1 | 2 | 3 | 4 |

In particular C=1/2, A=aD, B=dD in the nodal case. The transverse and
triangle axis coefficients are nonzero by the same retained node/frame
case split. Their principal divisors have degrees2,3,4, instead of12.

Define G_x=delta*partial_x+W_x, G_y=delta*partial_y+W_y. The normalized
connection source is S/delta^2, where

    S=W_t*delta-W*delta'
          -2*t^3*D^2-t^4*D*D'/2,
    deg_(t,x,y) S<=6.                            (5)

The term -2*t^3*D^2 is the logarithmic derivative -2/t of the volume
Jacobian in(3). Replacing the prefactor by delta^(-1/2) would be wrong.
Equivalently, the new source and the pulled-back old source differ by
the exact twisted image -(G_x(x)+G_y(y))/(t*delta). The checker verifies
this identity with the full symbolic ordinary-node Gaussian data.

All transformations are invertible over C(t). The retained original
formal Wick functional pulls back to the same seed I(t); it is not
redefined by evaluating a singular eliminated Gaussian at t=0. Cases
D(0)=0 are included. The final scalar equation acts on the original
I=1+sum_(j>=1) E(P^j)*t^j/j!.

## 2. Ordinary nodes: parameter degree plus amplitude degree

Use the four actual divergence operators and complete principal projections
from the ordinary-node proof, replacing V by W. They retain the full bases

    nodal:      1,x,y,x^2,xy,y^2,y^3, rank7;
    transverse: 1,x,y,x^2,xy,y^2,     rank6;
    triangle:   1,x,y,xy,y^2,         rank5.

The original written rank argument proves independence, and the explicit
strictly decreasing divergence reduction proves spanning. No ordinary
critical-ideal quotient is substituted for the twisted vector quotient.

Give t,x,y weight one. The phase gradients have total degree<=3 and
delta has degree<=3. For a monomial of ordinary amplitude degree M, one
complete principal projection, cleared by A0, produces coefficients and
remainders of total degree<=M+k, where

    A0=CAB, k=2 (nodal);
    A0=CAF, k=3 (transverse);
    A0=CEF, k=4 (triangle).                     (6)

Here each principal polynomial is homogeneous in the assigned upper
parameter degrees: C has degree0, A,B degree1 and E,F degree2.
The listed cleared principal cofactors have total degree<=M+k-3.
Multiplication by a phase gradient restores at most three; a divergence
term delta*partial restores at most two. This proves the bound for every
projection family and every amplitude degree. Same-degree remainders are
already basis elements. Every other correction reduces ordinary amplitude
degree by at least one.

For the unreduced actual seed derivatives write

    nabla^j[1]=[R_j/delta^(2j)], R_0=1,
    R_(j+1)=delta^2*R_j'-2j*delta*delta'*R_j+S*R_j. (7)

Induction using(5) gives total degree<=6j and ordinary x,y degree<=4j
for R_j. The required reduction depths for j>=1 are

    r_j=4j-2 (nodal/transverse), 4j-1 (triangle).  (8)

These use minimum reducible degrees3,3,2 and the strict degree decrease.
After padding by unused powers of A0, reduce to N_j/H_j with

    H_j=delta^(2j)*A0^r_j,
    deg_t N_(i,j)<=b_j-d_i, deg_t H_j<=b_j,
    b_j=14j-4,18j-6,22j-4, respectively.          (9)

The d_i are ordinary basis degrees. Uniform bound(9) comes from
6j+k*r_j and preserves every original parameter specialization for which
the selected principal divisor is nonzero.

## 3. Cusp and tangent: retain weighted reduction depths

The rescaled weighted top phases are now

    cusp:    x^4/2+D*y^3,    weights(3,4), rank6;
    tangent: x^4/2+D*x*y^2,  weights(2,3), rank5.

The retained weighted bases and exact principal divisions remain valid.
Their common divisors A0=a*b and a*b^2, with a=1/2,b=D, have maximum
t-degrees1 and2. A cleared principal cofactor has degree at most1 or2;
the phase-gradient coefficients have t-degree at most3, as does delta.
Thus each whole correction costs at most4 or5 numerator degrees.
In the tangent division of x^3, the already-normal y^2 term costs at
most3 degrees, which is also covered. No same-weight non-basis reduction
chain is discarded.

The lower weighted gaps are still2 and1. The seed R_j has weighted
amplitude degree at most12j or8j and parameter degree at most6j. Hence

    cusp:    r_j=6j-3, b_j=6j+4r_j=30j-12;
    tangent: r_j=8j-4, b_j=6j+5r_j=46j-20.       (10)

Reducing gives N_j/H_j with the same H_j as in(9), and all numerator
degrees at most b_j. Its denominator degree is no larger. The weighted
basis/rank proof and full analytic exclusions are retained unchanged.

## 4. Scalar cofactors and the initial-order comparison

Let s be the actual minimal seed rank. Set H_0=1 and N_0 equal to the
constant basis column. Select a nonzero s-row minor of the independent
columns N_0,...,N_(s-1); the rows include the constant row. Polynomial
cofactors yield a scalar equation

    sum_(j=0)^s a_j(t)*I^(j)=0,
    a_j=H_j*c_j.

The retained no-constant-normalized-period theorem makes a_0 nonzero.
This uses the actual minimal seed equation, not an assumption of full
seed cyclicity or of Morse critical points.

For ordinary nodes, keeping the selected basis grades gives

    deg a_0 <= sum_(j=1)^s b_j-sum_selected d_i.  (11)

For cusp and tangent, we use the coarser bound sum b_j. Every higher
coefficient contains H_j, and

    ord_(t=0) H_j >=4j,                          (12)

since delta=t^2D and A0 is polynomial. Thus if I=1+O(t^(K+1)), every
j>=1 error term has order at least K+1-j+4j>=K+1. The j=0 error has
the same lower bound. Taking K equal to the upper degree bound for a_0
contradicts its being a nonzero polynomial. There is no additional '+s'
derivative loss.

The worst rank is the full rank in each row: adding a column increases
the sum of b_j by more than any ordinary basis grade. Their resulting
D!=0 complement exclusions are

| Type | Rank | Degree sum used | Nonzero pure moment by |
|---|---:|---|---:|
| Irreducible node | 7 | 14*28-4*7-11 | **353** |
| Transverse line-conic | 6 | 18*21-6*6-8 | **334** |
| Independent triangle | 5 | 22*15-4*5-6 | **304** |
| Irreducible cusp | 6 | 30*21-12*6 | **558** |
| Tangent line-conic | 5 | 46*15-20*5 | **590** |

For ordinary nodes, the omitted isotropic zero-F2 branch has the retained
nondegenerate quadratic W coefficient and sufficient pure exclusion888.
For cusp and tangent, the D=0 branches use the retained425/320 frame
gates. These branches must be added back;353/334/304 are not full-class
ordinary-node gates. Consequently the full class gates are888 for each
ordinary-node type,558 for cusp and590 for tangent.

## 5. The complete global cover

| Leading type | Sufficient full-class pure gate/exclusion |
|---|---:|
| Degree at most two | 7 |
| Smooth cubic | 94 |
| Irreducible nodal cubic | 888 |
| Irreducible cuspidal cubic | 558 |
| Transverse line-smooth-conic | 888 |
| Tangent line-smooth-conic | 590 |
| Three independent lines | 888 |
| Three distinct concurrent lines | 104 |
| Repeated leading line(s) | 99 |

The last two bounds are retained from their complete canonical proofs.
The maximum is888. This proves(1), with arbitrary lower terms and no
extra frame test. The global pure-null members then have the Gaussian
W-linear frame and sharp mixed coefficient3 by the retained full theorem.

In the rational twenty-coefficient ring, radical(I888)=radical(I_infinity)
follows from(1) and Nullstellensatz. This is not ordinary ideal equality
or an expanded radical-membership replay. Exact radical tests against
I888 give at most ceil(log2(888-19))=10 binary-search rounds to find the
least K in[20,888]; no such full coefficient-ring round has been run.
All-homogeneous and odd subclasses inherit the even upper888. Parity
does not reduce an even gate to887.

## Evidence and remaining work

The checker verifies the full ordinary-node coordinate/gauge identity,
135 symbolic ordinary-node divergence identities and63 symbolic weighted
divergence identities. Its ten Gaussian-normal-form models include D(0)=0,
nondegenerate full Gaussian metrics and the retained zero axis-quadratic
charts. They supply58 complete connection-column and20 direct seed-column
integration-by-parts certificates. Both initial and strengthened successful
checker versions are preserved. No moments near888 are expanded.

The basis/rank completeness, every leading/frame case, the nonzero scalar
zeroth coefficient and the final all-order conclusion remain written
mathematics using retained analytic exclusions. The new rescaling and
degree arguments import no new external theorem. Current metadata checks
of Crew--Feldbrugge and Helmer--Tellander confirm their previously audited
versions; their generic formulas are not substituted for the degenerate
analytic arguments. Numerical leastness,20/14 sufficiency, independent
review, expanded radical cofactors, exceptional coefficient schemes,
global term-count minima and lower-power GVC remain unfinished.
