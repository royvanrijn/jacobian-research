# Square deformation closes the marked quadratic ruling-six class

PF2RULED188. 17 November 2026. Actual complex polynomial plane Keller maps,
within the complete marked quadratic obstruction model.

There is NO noninvertible actual Keller pair in the retained complete
quadratic source model with a target coordinate of ruling degree SIX.
This excludes the ENTIRE marked quadratic ruling-six obstruction class,
for arbitrary leading degree D, arbitrary source-coefficient degrees,
arbitrary companion degree and any number of pole supports.

The proof closes the last perfect-square first quartic left by PF187.
A deformation before the parameter cap contradicts two constant moment
relations; deformation at the cap gives a quadratic first parameter
response at every possible actual pole.

This is a closure of the stated model. It does NOT prove that every
global JC2/HC4 candidate enters this model, cover higher ruling degrees,
establish target trace perfection, or prove HC4/JC2.

## 1. The original actual model and the remaining first square

Retain the original source, fields and full normalizations from PF162:

    L0=C(x,t), E=C(u,w), F=C(P,Q), u=t^2,
    w^2=N+M*z, N=u*Sigma^2, M=u*A,
    w=t*(Sigma+A*x), z=2Sigma*x+A*x^2,
    dP wedge dQ=du wedge dw/M.                           (1)

Here V=Spec(C[u,w,z]/(w^2-N-M*z)) is the complete source image
used throughout. The coefficients are complex polynomials in u, with the retained
central marking ord_0 M>=2 and Sigma(0)!=0. The target pair is ACTUAL
and noninvertible, not an abstract embedding of a comparison algebra.
Every original full-normalization prime is retained.

Suppose P has ruling degree six, meaning deg_w P*=6 after substituting
z=(w^2-N)/M. This is not the generic map degree. The complete-end
argument PF162/PF164 makes every affine point, M-root end and
u-infinity end holomorphic for the actual coframe dQ.
The simple-escape proof PF165 excludes the B/C pole faces.
Thus there must be at least one ACTUAL A pole support c off M:
a nonconstant rational Q on the complete generic P curves has poles,
and all those poles must be among the retained A supports.

Use the unrestricted regular normal form and balanced remainder:

    P=H^2+R, H=H3 in O(V),
    R=d(u)z+e(u)w+f(u), d,e in C[u] nonzero.              (2)

No source function H or R is put into F or assumed integral over
the target. Fix any actual support c, put v=u-c and Y=v^L*w.
The previous results give

    L>=3, M(c)!=0, ord_c h=2L, ord_c K=L,
    v^L*H(v,v^-L*Y) mod v=C*Y*(Y-gamma)^2,
    C,gamma in C*,
    t=L+ord_c e, L<=t<=2L-2.                             (3)

The entire equation and its monic separated depressed quartic are

    F0(v,Y)-p*v^(2L)=F2*G,
    G=X^4+2J*X^2+q*X+J^2+Delta,
    X=Y-m_G(v,p),                                      (4)

where F2 is a unit near the quartic cluster, all scaled coefficients
are regular, and p first occurs at order 2L.

PF187 proves that every first nonsquare initial, including repeated
initial roots and the first-scale cap edge, is impossible. Therefore
the ONLY remaining first region satisfies

    ell=ord_v J is finite, ell<L,
    t=ord_v q, 3ell<2t,
    delta=ord_v Delta>2ell.                             (5)

These orders are first taken over C(p). The universal parameter
precision theorem PF182 gives delta<=2L. Thus no order or degree
cutoff is added here.

## 2. The full constant companion relations

PF165 gives a globally reduced actual companion degree 3k,
k>=3 odd. The full actual Jacobian and noncube Kummer equation of
PF173/PF176 yield complex constants c_n, c_k!=0, with

    Q-sum_(n=0)^k c_n*P^(n/2)=O(w^-3).                  (6)

The branch of sqrt(P) is selected by its polynomial part H.
Define coefficient residues for every positive odd n,

    m_n=[w^-1]P^(n/2),
    l_n=[w^-2]P^(n/2)=[w^-1]w*P^(n/2).                 (7)

The polynomial Q and every even half-power have zero coefficients.
Consequently the SAME constants satisfy

    sum_(1<=n<=k,n odd) c_n*m_n=0,
    sum_(1<=n<=k,n odd) c_n*l_n=0.                       (8)

For any local p-independent Laurent center b(v), the centered coefficient

    lhat_n=l_n-b(v)*m_n

has the same zero relation. Multiplication by a scalar Puiseux power
of v also preserves these local identities. No such scalar or center
is asserted to be a target function.

## 3. The exact remaining square deformation and the parameter split

Put

    rho=ell/2, 0<rho<L/2,
    m=min(delta,t+rho).

By (5) and the parameter cap,

    4rho<m<=2L,
    nu=m-4rho>0.                                       (9)

We distinguish m<2L and m=2L BEFORE evaluating any cap coefficient
at p=0. All coefficients below order 2L are independent of p.

If m<2L, its first nonzero deformation is strictly before the
parameter term, so its leading coefficient is unchanged at p=0.
This remains true even if a later coefficient at order 2L vanishes
there. No generic constant target translation is required.

At p=0 use X=v^rho*Z and the exact quartic midpoint center
b(v)=v^-L*m_G(v,0). The full quartic becomes

    v^(-4rho)*G(v,v^rho*Z,0)
       =(Z^2+j(v))^2
              +v^nu*(A0*Z+B0+higher positive-order terms),
    j(0)=j0!=0, (A0,B0)!=(0,0).                         (10)

The remainder is polynomial in Z of degree at most one.
The complete square j(v)=v^-ell*J(v,0) is retained.
Specifically the q term first occurs at order t-3rho,
and the Delta term at order delta-4rho; their minimum is nu.
If both orders agree, either leading coefficient may cancel at one
square center, but A0 and B0 cannot both vanish.

After a finite Puiseux substitution all rescaled coefficients
are convergent power series. The entire source polynomial has

    P(v,w)=v^(4rho-2L)*U(v,Z)
                 *[(Z^2+j(v))^2+v^nu*T(v,Z)],
    w=b(v)+v^(rho-L)*Z,
    U(0,Z)=tau^2 in C*, T(0,Z)=A0*Z+B0.                 (11)

The other separated factor U is an analytic unit on a fixed disc
containing all four escape roots, and has a holomorphic square root
there. This retains its full variation, not a frozen quotient.

## 4. Before the cap, the entire moment vector grows strictly

The source polynomial specializes to a nonconstant quadratic at c.
Its two finite roots stay inside a fixed circle; the other four
escape. Both clusters have EVEN root counts, so the selected sqrt(P)
is single-valued outside the two discs. Cauchy's theorem splits each
coefficient-residue integral into its finite and escape integrals,
as in PF184.

For each fixed positive odd n,

    finite contribution to m_n is O(1),
    finite contribution to lhat_n is O(v^-L).             (12)

Both finite roots are retained in these estimates.

Write n=2r-1, r>=1. On a fixed escape circle expand the whole factor

    U(v,Z)^(n/2)*(Z^2+j(v))^n
      *sum_j binom(n/2,j)
            [v^nu*T(v,Z)/(Z^2+j(v))^2]^j.               (13)

Every term j<r is holomorphic on the ENTIRE escape disc:
its quadratic-base exponent n-2j is nonnegative, and the square
root of U is holomorphic. Its integral is therefore EXACTLY zero.
Multiplying these integrands by either 1 or Z preserves
holomorphy, so this applies to BOTH vector components.
The first possible term is j=r. All terms j>r are uniformly of
strictly higher v order on the contour.

Choose z0 in C* with z0^2=-j0, and put

    U0=B0+A0*z0, V0=B0-A0*z0,

    M_r=(U0^r-V0^r)/(2z0),
    N_r=(U0^r+V0^r)/2.                                  (14)

These are respectively the coefficient residues of
(A0Z+B0)^r/(Z^2+j0) and Z*(A0Z+B0)^r/(Z^2+j0).
The pair (M_r,N_r) is NONZERO for EVERY positive integer r:
otherwise U0^r=V0^r=0 and hence A0=B0=0.
This includes one-sign residual cancellation and all root-of-unity
ratios between U0 and V0.

The exact escape asymptotic of the ALIGNED vector is

    (m_n, v^(L-rho)*lhat_n)
       =tau^n*binom(n/2,r)
          *v^(-p_r)*(M_r,N_r)+o(v^(-p_r)),
    p_r=r*(2L-m)+rho.                                   (15)

Indeed dw contributes v^(rho-L), P^(n/2) contributes
v^((4rho-2L)*n/2), and the first residue term contributes
v^(nu*r). Their sum is -rho-r*(2L-m).
The second coefficient has one extra factor v^(rho-L),
removed by the displayed alignment.

The finite contributions in (12) have aligned upper pole bounds
zero and rho. If m<2L,

    p_r>rho and p_(r+1)>p_r.                             (16)

Thus they cannot change the nonzero leading vector (15).
The half-integer binomial coefficient is also nonzero.
The exact WHOLE vector has strictly growing pole order with
the odd index.

The top nonzero constant c_k in the two zero identities (8)
therefore multiplies a vector with strictly larger pole than every
lower odd character. Its nonzero leading vector cannot cancel.
This excludes the ENTIRE case m<2L, for arbitrary slope,
coefficients and companion index.

## 5. At the cap, the whole first response is quadratic

It remains m=2L. Use p=0 only to choose a p-independent formal
center, not to assert a nonzero specialized cap coefficient.
At that specialization J still has exact order ell, while

    ord Delta(v,0)>=2L, t+rho>=2L.

Let x_sign=+sqrt(-J(v,0)) or its negative, and let Yc be this
center added to the p=0 quartic midpoint. Both are
C-coefficient Puiseux series. The entire p-independent polynomial
F0 satisfies at Yc

    ord F0>=2L,
    ord F0_Y>=t,
    ord F0_YY=ell EXACT,
    ord F0_YYY>=rho,
    all higher Y derivatives have order>=0.             (17)

To check (17), use the full p=0 factorization F0=F2*G.
At X=x_sign,

    G=Delta+q*x_sign has order at least 2L,
    G_X=q has order t,
    G_XX=-8J has order ell,
    G_XXX=24*x_sign has order rho.

The other factor is a regular unit. Since ell<L<=t,
the second derivative has a unique leading product term
and cannot cancel. This proves the bounds for the WHOLE
source polynomial.

Put R0=L-rho. The Taylor substitution

    Y=Yc+v^R0*Z

has exact quadratic order ell+2R0=2L. Every other nonconstant
Taylor term is strictly later:

    linear: t+R0>=3L-ell>2L,
    cubic: rho+3R0=3L-ell>2L,
    degree at least four: 4R0>2L.                       (18)

Thus the entire P=p equation, after division by v^(2L),
has initial

    a2*Z^2+c0-p, a2 in C*, c0 in C.                    (19)

The p=0 cap coefficient is allowed to vanish or to cancel at
one center; only its lower bound in (17) was used.
The generic term -p in (19) cannot vanish.

At each of the two centers (19) has TWO distinct roots over
an algebraic closure of C(p), with nonzero derivative.
They lift to four Puiseux roots of the full equation.
Their displacements are smaller than the center separation,
because L-rho>rho. They exhaust the four escaping sheets;
none is omitted or counted as a different normalized point
without its ramification multiplicity.

At any ACTUAL pole among these sheets, the first coefficient
outside C satisfies a genuine QUADRATIC equation
a2*lambda^2+c0=p. Replacing the infinite formal center by
its finite C-coefficient truncation changes lambda only by
a constant and does not change degree two.
PF159 requires the full first response at every actual
coframe pole to be AFFINE. This is a contradiction.

There must be an actual pole among these four sheets, because
c was chosen as an actual A support and the other two roots
are finite and holomorphic. One may also check directly that
t<=2L-2 and t+rho>=2L force rho>=2: the cap packets then
have mu=L+rho and are all actual Q poles. The argument
does not need that extra observation.

This excludes m=2L, including every possible leading cap
cancellation and all ramification.

## 6. Closure and exact mathematical scope

PF187 excluded all nonsquare first initials and the first-scale
edge rho=L/2. Its only remaining region was (5).
Sections 4 and 5 now exclude BOTH possibilities in (9).
Hence no actual A pole support exists.

The complete-end and simple-escape results require at least
one such support for every actual noninvertible ruling-six
pair in the retained marked quadratic model. This final
contradiction proves the theorem.

The proof covers every leading degree D, every slope L,
every number of supports and all coefficient and companion
degrees. There is no surviving D33 or larger-D frontier
inside THIS ruling-six class. Conditional numerical bounds
proved in earlier notes remain valid intermediate implications;
they are not the remaining obstruction after this closure.

The theorem does not assert that arbitrary polynomial maps
with a degree-six expression fit this model, nor that all
global JC2/HC4 cases enter it. Higher ruling degrees and
unverified global model/import coverage, target derivative
integrality, trace perfection and HC4/JC2 remain unresolved.

## 7. Exact evidence and replay

The checker verifies the unrestricted residue pair, exact scaling
and strict growth identities, whole sextic moments with two finite
simple roots for integer and half-integer first scales, one-sign
cancellation and pure linear/constant residual controls, the exact
generic-order split and all four whole-polynomial cap packets.
Finite controls supplement the written proof; they do not prove
unbounded statements by a coefficient prefix.

Replay:

    .venv/bin/python research/jacobian/affine-complements/verify_ruling_six_square_deformation_closure.py

Essential retained sources include
[PF159](POLYNOMIAL_MATE_RESPONSES_AND_QUADRATIC_RULING_FIVE_ESCAPE_COVERAGE.md),
[PF162](RULING_SIX_FINITE_END_COVERAGE_AND_QUADRATIC_BOUNDARY_ALGEBRA.md),
[PF164](RULING_SIX_AFFINE_POLE_BUDGET_AND_SOURCE_INFINITY.md),
[PF165](RULING_SIX_SIMPLE_ESCAPE_GRAPHS_AND_QUARTIC_FACE_REDUCTION.md),
[PF171](RULING_SIX_REGULAR_CUBIC_APPROXIMATE_ROOT_AND_QUADRATIC_REMAINDER.md),
[PF173](RULING_SIX_LEADING_PRIMITIVE_AND_POLYNOMIAL_ODE.md),
[PF176](RULING_SIX_COMPANION_RESONANCE_AND_LINEAR_REMAINDER_CLOSURE.md),
[PF178](RULING_SIX_BALANCED_QUADRATIC_REMAINDER_AND_LEADING_RESONANCE.md),
[PF182](RULING_SIX_PARAMETER_PRECISION_AND_SLOPE6_DEGREE30_EXCLUSION.md),
[PF184](RULING_SIX_TWO_MOMENT_RESONANCE_AND_DEGREE30_CUBIC_REDUCTION.md),
and
[PF187](RULING_SIX_UNIVERSAL_THREE_MOMENT_FIRST_QUARTIC_REDUCTION.md).
All actual-source and original-normalization hypotheses are retained.
Internal derivation and replay are not formal verification or
external human review.
