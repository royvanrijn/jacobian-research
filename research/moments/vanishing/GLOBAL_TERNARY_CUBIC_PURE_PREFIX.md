# A numerical global pure-moment prefix for ternary cubics

Current quantitative refinement, 15 November 2026: the [exact Euler representative proof](EULER_REPRESENTATIVES_GLOBAL_297.md) gives global pure297, odd296, conic/full-W-linear258 and all-homogeneous258. Ordinary-node complements are297/271/244, with homogeneous nonisotropic-node101; cusp255/tangent284 are unchanged. Current intervals are global[20,297], odd[20,296] and all-homogeneous[14,258]. Earlier353/336 bounds below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I297 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [joint-degree actual-seed proof](JOINT_DEGREE_GLOBAL_353.md) gives global pure353, odd352, full W-linear336, rank-one140, mixed-parent255 and cusp255/tangent284. Current least-prefix intervals are global[20,353], odd[20,352] and all-homogeneous[14,336]. Earlier590/888 bounds and rank-one425/mixed-parent320 estimates below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I353 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [balanced logarithmic seed proof](BALANCED_LOGARITHMIC_GLOBAL_590.md) gives global pure590, all-homogeneous336, conic336 and full W-linear425. Current intervals are global/odd[20,590] and all-homogeneous[14,336]; smooth[20,94] and homogeneous-smooth[14,56] are unchanged. Earlier888/1855 bounds and frame-plus888 statements below retain their historical scopes. Leastness,20/14 sufficiency, expanded I590 cofactors and independent review remain open.

Current global refinement,15 November2026: the [rescaled singular seed proof](RESCALED_SINGULAR_GLOBAL_888.md)
now proves888 sufficient as a global pure-only prefix, without an extra
frame condition. Current intervals are global[20,888], odd[20,888] and
all-homogeneous[14,888]. Earlier1855 bounds and reservations about888
alone below are historical. The new proof remains internal and
independently unreviewed; leastness and expanded global cofactors are open.

15 November 2026. Internal written proof with exact auxiliary divergence
identities and model connection certificates. Analytic exclusions are
retained from the earlier canonical class proofs. No independent specialist
replay, external review or formal verification. The full twenty-coefficient
moment ideal is not expanded, and no moment1855 computation is performed.

The subsequent [graded smooth seed refinement](GRADED_SEED_SMOOTH_CUBIC_PREFIXES.md)
reduces the smooth class bound188 to94 (homogeneous-smooth56). The
original cover below remains a valid coarser proof; its global maximum1855
and least-prefix interval[20,1855] are unchanged.

## Results

For complex ternary cubics with arbitrary lower terms:

| Claim | Leading curve | Some positive pure moment is nonzero by |
|---|---|---:|
| **G3NODALFINITE20261115** | Irreducible ordinary nodal cubic | 1855 |
| **G3TRANSVERSEFINITE20261115** | A line transverse to a smooth conic | 1140 |
| **G3TRIANGLEFINITE20261115** | Three linearly independent lines | 995 |

**G3PUREGLOBAL20261115.** For every complex polynomial P of total
degree<=3 in three independent standard real Gaussian coordinates,

    E(P^j)=0 for j=1,...,1855
      iff E(P^j)=0 for every j>=1.               (1)

Every such P has a Gaussian W-linear frame and satisfies
E(QP^m)=0 whenever m>3*deg Q. The mixed coefficient3 is sharp by the
retained cubic shear. The sufficient pure prefix1855 is not claimed
optimal. The least global pure prefix K satisfies

    20<=K<=1855.                                 (2)

The subsequent [odd-invariant proof](ODD_CUBIC_CLASSIFICATION_AND_PREFIX_LOWER_BOUND.md)
raises the lower bound from the displayed cusp10 to20 by algebraic
existence of19 initial zero moments without all-order nullity. The subsequent
[odd isolated-root witness](ODD_CUBIC_TWENTIETH_MOMENT_WITNESS.md) now
supplies coefficients by an exact algebraic isolating representation and
proves the twentieth moment nonzero, with smooth leading form. The upper
proof below is unchanged. Thus a numerical global upper bound is proved internally. The
optimal numerical value in **OP-GMC3-PURE-PREFIX** remains open. Equation(1)
uses only pure moments; the separate frame condition in the older
frame-plus888 certificate is not required at the larger prefix1855.
Neither888 nor any smaller unproved number is silently promoted to a
global pure-only gate.

## 1. Ordinary-node normal forms and the frame alternative

If a leading ordinary node v is isotropic and F2(v)=0, the exact
Gaussian frame criterion gives a W-linear frame with a nondegenerate
binary quadratic W coefficient: the tangent cone of an ordinary node
has rank2. The retained [conic finite theorem](QUADRATIC_W_LOG_PERIODS.md)
pure-excludes that member by888. Remove this already quantified branch.
Then every zero-F2 leading node is nonisotropic.

Choose an ordinary node as T-axis. Normalize its tangent cone to uv.
A linear T shear absorbs the mixed binary cubic terms; translations
of u,v remove uT,vT; a constant T translation removes the quadratic uv
term. None of these operations discards the T^2 coefficient. They give

    P=alpha*T^2+T*(uv+h)+a*u^3+d*v^3
                  +b*u^2+c*v^2+eu+fv+k,
    Q=m*T^2/2+T*(ju*u+jv*v+j0)+q2(u,v),         (3)

with the full nondegenerate affine Gaussian Q retained. The three
leading types are a*d!=0 (irreducible nodal), a!=0,d=0 (transverse
line-conic, after swapping), and a=d=0 (independent triangle). They
have respectively1,2,3 ordinary nodes. The residual line-conic is
u*(Tv+a*u^2), with smooth quadratic cofactor and two transverse nodes.

Put D=m-2alpha*t. It is nonzero over C[t]: D=0 would mean that the
chosen node is isotropic with F2=0, already removed. Completing T gives
prefactor constant*D^-1/2 and phase U/D, where

    U=D*(tG-q2)+(t*(uv+h)-ju*u-jv*v-j0)^2/2.

Now translate old_u=u+jv/t, old_v=v+ju/t. The volume Jacobian is1,
the mixed cubic coefficients cancel, and the phase is V/delta with

    delta=t^2*D, deg_t V<=4, deg_t delta<=3,
    V=C*u^2*v^2+A*u^3+E*u^2+B*v^3+F*v^2
                         +J*u*v+L*u+N*v+K0.     (4)

Here the capital letters in(4) are coefficient polynomials in t, each
of degree<=4; they are unrelated to the earlier full Gaussian matrix
entries. The [accepted exact checkpoint](ordinary_node_phase_shift_checkpoint.json)
and its script verify

    C=t^4/2, A=a*t^3*D, B=d*t^3*D,
    E=t^2*D*(tb+3a*jv-q_uu/2),
    F=t^2*D*(tc+3d*ju-q_vv/2).                  (5)

Consequently the necessary principal coefficients are all nonzero:

| Type | Nonzero coefficients |
|---|---|
| Nodal | C,A,B |
| Transverse | C,A,F; B=0 |
| Triangle | C,E,F; A=B=0 |

For the latter two rows, if a corresponding original binary quadratic
coefficient is zero, the original Q value at that other leading node
is nonzero, since the isotropic zero-F2 alternative was removed. Thus
tc-q_vv/2 and, for triangles, tb-q_uu/2 cannot vanish identically.
This is a full parameter cover, not a generic diagonal-norm assumption.

All generic-t coordinate changes are invertible and preserve the unit
volume seed up to a nonzero constant. The rational translation contributes
no determinant factor. Its parameter derivative differs from the
coordinatewise connection only by vertical exact derivatives. The
prefactor remains D^-1/2, not delta^-1/2. D(0)=0 is allowed: the formal
Wick functional is inherited from the original full Gaussian before
elimination; no singular conditional average at t=0 is used.

## 2. Full ranks7,6,5, including the alpha!=0 charts

For generic t, the quadratic homogeneous part of S_t=tP-Q is nonzero
at every leading node. Indeed its value there is tF2(v)-Q2(v); if
F2(v)=0, nonisotropy makes Q2(v)!=0. Avoiding finitely many t zeros
therefore suffices. The retained ordinary-node gradient criterion makes
S_t complex linearly gradient-coercive at each such t.

The three homogenized quadratic critical equations have only the
ordinary leading nodes at infinity. In a chart around one node, the
linear transverse block is the rank2 tangent Hessian times t, and the
remaining linear coefficient is the nonzero quadratic value of S_t
along the node. Hence its local three-equation Jacobian is invertible,
so each infinity node has intersection length1. A positive-dimensional
affine critical component would have a projective closure meeting
infinity, contradicting the isolation of these points. The quadrics
thus form a proper complete intersection of degree8. Removing the1/2/3
infinity points leaves affine critical length7/6/5.

The retained coercive/tame comparison identifies these lengths, including
all finite Milnor multiplicities, with the full top twisted de Rham ranks.
Exact quadratic elimination and the rational translations preserve them.
This extends the retained alpha=0 count to every alpha chart here; it
does not infer rank from a sampled critical ideal.

## 3. Divergence operators, rather than an ordinary Jacobian ideal

Write G_u=delta*partial_u+V_u and G_v=delta*partial_v+V_v. Their images
define the twisted vector quotient. It is not a ring quotient, so a
polynomial relation cannot simply be multiplied without derivative terms.
Use the following actual image operators on every amplitude z(u,v):

    R1(z)=G_u(z), R2(z)=G_v(z),
    R3(z)=G_u(u*z)-G_v(v*z),
    R4(z)=G_u(C*u*v*z)+G_v((-C*v^2-3A*u/2-epsilon*E)*z), (6)

where epsilon=0 for nodal and1 for transverse/triangle. The polynomial
part of R3 (its two divergence terms cancel) is

    r3=3A*u^3-3B*v^3+2E*u^2-2F*v^2+L*u-N*v.

The polynomial part of R4 is

    C*v*r3-(3A*u/2+epsilon*E)*V_v-delta*C*v.     (7)

The last term is essential: it is the divergence of the coefficients
in(6). It is not a commutative S-polynomial correction.

The principal ordinary-degree polynomials p_i and full bases are:

| Type | p1,p2,p3,p4 | Basis |
|---|---|---|
| Nodal | 2Cuv^2,2Cu^2v,3Au^3-3Bv^3,-3CBv^4 | 1,u,v,u^2,uv,v^2,v^3 |
| Transverse | 2Cuv^2,2Cu^2v,3Au^3,-2CFv^3 | 1,u,v,u^2,uv,v^2 |
| Triangle | 2Cuv^2,2Cu^2v,2Eu^2-2Fv^2,-2CFv^3 | 1,u,v,uv,v^2 |

The other polynomial terms of each R_i have strictly lower degree than
p_i. Their differential coefficients have degree<=2 in u,v.

For clarity, these are complete principal projections. Mixed monomials
divisible by uv^2 or u^2v use p1 or p2 (the latter is unnecessary in
the triangle row). In the nodal row,

    u^3=(B/A)*v^3+p3/(3A),
    u^i=u^(i-3)*p3/(3A)+B*u^(i-4)*v*p1/(2AC), i>=4,

and pure v^j,j>=4 use -v^(j-4)*p4/(3CB). In the transverse row pure
u^i,i>=3 and v^j,j>=3 use p3/(3A) and -p4/(2CF) with the appropriate
monomial multiplier. In the triangle row,

    u^2=(F/E)*v^2+p3/(2E),
    u^2*v=v*p3/(2E)-p4/(2EC),
    u^i*v^j=u^(i-2)*v^j*p3/(2E)
             +F*u^(i-3)*v^j*p1/(2EC), i>=3,j<=1.

The remaining mixed monomials use p1 and pure v^j,j>=3 use p4.
Same-degree remainders are already basis terms; no unbounded sequence
of equal-degree corrections is hidden here.

For a principal projection f=NF+sum h_i*p_i, replace p_i*h_i by the
negative lower polynomial and differential parts of R_i(h_i). Every
non-basis correction has ordinary degree at most deg f-1. Induction
proves spanning by the listed7/6/5 monomials. The full rank count in
Section2 then proves independence. Thus these are full quotient bases,
without assuming a commutative Groebner basis is a twisted basis.

## 4. Uniform parameter degrees and numerical prefixes

Use common principal divisors

    A0=C*A*B (nodal), C*A*F (transverse), C*E*F (triangle).

Each has t-degree<=12. After clearing by A0, the projections above have
polynomial coefficients in t. The combined original G_u/G_v cofactors
have t-degree<=8. Multiplying by V_u,V_v raises this to at most12;
their derivative/divergence terms multiplied by delta have degree<=11.
The cleared basis remainder has degree<=12. Thus every whole correction
level costs at most one A0 denominator and12 numerator degrees.

The normalized connection source is

    S/delta^2=(V_t*delta-V*delta')/delta^2-D'/(2D).

Its numerator S is polynomial of t-degree<=6 and ordinary u,v degree<=4,
since delta^2*D'/D=t^4*D*D'. A connection column starts at degree at
most7 (nodal) or6 (other rows). The minimum reducible degrees are3,3,2.
Strict degree decrease therefore gives correction depths r=5,4,5.
Padding basis terms by the unused powers of A0 gives a connection N/H:

    H=delta^2*A0^r,
    deg_t H, deg_t N_ij <= B0=6+12r.

The B0 values are66,54,66. Minimal-seed polynomial cofactors, exactly
as in the retained finite-prefix proofs, yield scalar coefficient degree
at most B0*s*(s+1)/2 at seed rank s<=7/6/5. The zeroth coefficient is
nonzero: otherwise its minimal seed equation would admit the constant
solution, giving a rational trivial seed quotient. The retained complete
class exclusions forbid that quotient (or its actual constant-period
realization on the nonsemisimple critical-axis/line branches). No cyclic
seed in the full ambient module is assumed.

The formal initial-order comparison then excludes a zero prefix of
B0*s*(s+1)/2+s. At the full ranks these are

    nodal:66*7*8/2+7=1855;
    transverse:54*6*7/2+6=1140;
    triangle:66*5*6/2+5=995.

The isotropic zero-F2 ordinary-node branch was already excluded by888,
below all three numbers. This proves the three class statements on
their entire coefficient spaces. The numbers are sufficient, not optimal.

## 5. The global cover and the optimality interval

The exhaustive projective leading-cubic cover now has numerical bounds:

| Class | Sufficient pure-zero gate, or pure exclusion |
|---|---:|
| Degree<=2 | 7 |
| Smooth leading cubic | 188 |
| Irreducible nodal | 1855 |
| Irreducible cusp | 720 |
| Transverse line-smooth-conic | 1140 |
| Tangent line-smooth-conic | 935 |
| Independent three lines | 995 |
| Concurrent three distinct lines | 104 |
| Repeated leading lines | 99 |

The retained canonical class proofs and the new ordinary-node calculation
cover every member, with arbitrary lower terms. The largest bound is1855,
so the first1855 zero moments force all-order pure nullity. The retained
internal full cubic theorem then gives a W-linear frame and sharp mixed
cutoff3*deg Q. The converse in(1) is immediate. The retained sharp cusp
Z^2*T+W^3, whose first9 moments vanish and tenth does not, proves the
earlier lower bound10. The subsequent odd-invariant argument above raises
it to20 without displayed coefficients. Nothing here determines the optimal K.

In the rational twenty-coefficient ring, put I_j=(M1,...,M_j). Equation(1)
and Nullstellensatz give radical(I1855)=radical(I_infinity), without any
claim of ideal equality. The least K can now be found by a bounded exact
search: test all M_(K+1),...,M1855 in radical(I_K), using inverse-variable
unit-ideal certificates. Monotonicity and the known false value19/true
value1855 permit at most ceil(log2(1855-19))=11 binary-search rounds.
Each round may contain many large exact polynomial calculations. No
such global round, expanded radical cofactor certificate or numerical
optimality computation has been run in this follow-up.

The same1855 bound applies uniformly to algebraically specified constrained
cubic families. Their smaller least thresholds remain computable by the
retained constructible-family algorithm. This does not determine moment
ideal generation, nilpotence, primary decompositions, term-count minima
or lower-power GVC.

## Evidence and assurance

The [checker](verify_ordinary_node_finite_prefixes.py),
[receipt](ordinary_node_finite_prefixes_verification.json),
[six model certificates](ordinary_node_connection_models.json),
[accepted first-run record](ordinary_node_finite_prefix_first_run.json)
and [source audit](ordinary_node_finite_source_audit.json) bind92 universal
principal/divergence tests and all36 connection-column identities. Three
models have D(0)=0; zero original axis-quadratic coefficients are retained.
All model Gaussian Hessians are nondegenerate. The earlier normalization
checkpoint and unverified plan remain preserved as history.

The auxiliary checker does not establish universal rank/tameness, the
complete leading-case cover, the analytic nonzero zeroth coefficient or
the all-order global conclusion. Those are the written arguments and
retained analytic dependencies stated here. Independent specialist replay,
external review, formal verification and explicit full-coefficient radical
cofactors remain outstanding. The broader research goal is not complete.
