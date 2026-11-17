# Graded seed bounds: smooth cubic pure prefixes94 and56

Current quantitative refinement, 15 November 2026: the [exact Euler representative proof](EULER_REPRESENTATIVES_GLOBAL_297.md) gives global pure297, odd296, conic/full-W-linear258 and all-homogeneous258. Ordinary-node complements are297/271/244, with homogeneous nonisotropic-node101; cusp255/tangent284 are unchanged. Current intervals are global[20,297], odd[20,296] and all-homogeneous[14,258]. Earlier353/336 bounds below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I297 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [joint-degree actual-seed proof](JOINT_DEGREE_GLOBAL_353.md) gives global pure353, odd352, full W-linear336, rank-one140, mixed-parent255 and cusp255/tangent284. Current least-prefix intervals are global[20,353], odd[20,352] and all-homogeneous[14,336]. Earlier590/888 bounds and rank-one425/mixed-parent320 estimates below are historical valid coarser evidence. Leastness,20/14 sufficiency, expanded I353 cofactors and independent review remain open.

Current quantitative refinement, 15 November 2026: the [balanced logarithmic seed proof](BALANCED_LOGARITHMIC_GLOBAL_590.md) gives global pure590, all-homogeneous336, conic336 and full W-linear425. Current intervals are global/odd[20,590] and all-homogeneous[14,336]; smooth[20,94] and homogeneous-smooth[14,56] are unchanged. Earlier888/1855 bounds and frame-plus888 statements below retain their historical scopes. Leastness,20/14 sufficiency, expanded I590 cofactors and independent review remain open.

Current global refinement,15 November2026: the [rescaled singular seed proof](RESCALED_SINGULAR_GLOBAL_888.md)
now proves888 sufficient as a global pure-only prefix, without an extra
frame condition. Current intervals are global[20,888], odd[20,888] and
all-homogeneous[14,888]. Earlier1855 bounds and reservations about888
alone below are historical. The new proof remains internal and
independently unreviewed; leastness and expanded global cofactors are open.

15 November 2026. Internal written refinement of the retained
[smooth-leading finite-prefix proof](COERCIVE_GRADIENT_MOMENTS.md).
The analytic no-constant-period input is unchanged. The improvement comes
from reducing the actual seed derivatives P^j, keeping the degrees of
the quotient basis and the different orders of the scalar coefficients.

| Class of complex ternary cubics | Sufficient pure exclusion prefix | Current least-prefix interval |
|---|---:|---:|
| Smooth cubic leading form, arbitrary lower terms | **94** | **[20,94]** |
| Smooth homogeneous cubic | **56** | **[14,56]** |

Some positive pure Gaussian moment is nonzero by the indicated order.
These are class exclusions, not global gates for singular-leading cubics.
The [odd20 witness](ODD_CUBIC_TWENTIETH_MOMENT_WITNESS.md) and
[homogeneous14 witness](SMOOTH_CUBIC_FOURTEENTH_MOMENT_WITNESS.md)
supply the lower bounds. Optimality is not asserted. The global interval
remains[20,1855], and all-homogeneous[14,1854]: in particular irreducible
nodal homogeneous cubics are not silently assigned the smooth56 bound.

There is also an all-dimension cubic refinement. For N>=1 and mu=2^N,

    smooth cubic: K <= mu*(3mu+3-N)/2-2;
    smooth homogeneous cubic: K <= mu*(mu-1).      (1)

The old general-degree formula remains a valid coarser bound in degree
three and retains its stated scope in degrees greater than three.

The [checker](verify_smooth_cubic_seed_bounds.py),
[model certificate](smooth_cubic_seed_bounds_certificate.json), and
[receipt](smooth_cubic_seed_bounds_verification.json) support the algebra
below. Universal filtered reduction, the seed argument and the retained
analytic exclusion remain written mathematics. No independent specialist
replay, external review or formal verification is claimed. The
[source audit](smooth_cubic_seed_bounds_source_audit.json) keeps this boundary
and the immutable former numerical frontier.

## The graded twisted quotient

Put z=1/t and let P=F3+F2+F1+F0 with smooth homogeneous leading F3.
The regular sequence partial_i F3 has quotient Hilbert series(1+w)^N.
Choose any homogeneous representatives b_i of that Jacobian quotient,
including b_0=1. Write d_i=deg b_i. They lift to a basis of the twisted
vector quotient

    H=C(t)[x]/sum_i (partial_i+t*partial_i P-x_i)C(t)[x].

This is a vector quotient, not an ordinary quotient ring. Its dimension
is mu=2^N. For N=3 its basis grades are

    0,1,1,1,2,2,2,3.                             (2)

All coefficient-dependent basis choices are made for the actual smooth
form; no fixed Hesse chart or minor is imposed on the whole class.

If a homogeneous input component has degree M, use a homogeneous leading
division by partial_i F3. For a quotient term (partial_i F3)*h_i, with
deg h_i=M-2, the exact twisted relation is

    [(partial_i F3)*h_i]
      =-[(partial_i F2+partial_i F1)*h_i]
        +z*[x_i*h_i-partial_i h_i].              (3)

The four resulting degree changes are:

| Term | New x-degree | New z-power |
|---|---:|---:|
| -(partial_i F2)*h_i | M-1 | 0 |
| -(partial_i F1)*h_i | M-2 | 0 |
| z*x_i*h_i | M-1 | 1 |
| -z*partial_i h_i | M-3 | 1 |

Every replacement strictly lowers x-degree, so reduction terminates.
Starting from an input of degree at most M, the coefficient of b_i is
a polynomial in z of degree at most M-d_i: each acquired z lowers the
x-degree by at least one.

The constant row has a stronger bound, degree at most M-2 for M>=2.
A path that creates a constant must end either with a derivative term
(degree loss three for one z), or with the linear-part gradient
(degree loss two without a z). Neither a linear term x_i*h_i nor the
quadratic-part gradient can be a nonzero constant. Homogeneous normal
projection preserves degree. A constant already present in the input
has z-degree zero. Thus the final constant-producing step loses two
more units of x-degree than the number of z factors gained, proving
the bound. Inputs below degree two require no reduction and satisfy
the same bound whenever the overall degree bound M>=2 is used.

## Actual seed columns rather than the full connection bound

For nabla[u]=[partial_t u+P*u], the seed derivatives are exactly

    v_j=nabla^j[1]=[P^j].                         (4)

The representative P^j itself has no t dependence. Consequently, for j>=1,

    v_(i,j) in C[z],
    deg_z v_(i,j)<=3j-d_i,
    deg_z v_(0,j)<=3j-2.                         (5)

These bounds replace the older common-denominator estimate obtained by
applying the whole connection matrix repeatedly. They preserve exceptional
seed ranks and do not require the seed to generate all of H.

Let s<=mu be the actual minimal seed rank: v_0,...,v_(s-1) are independent,
and v_s is dependent. Choose a nonzero s-row minor of the independent
columns and use the same rows for all cofactors c_j(z), j=0,...,s.
The row set includes the constant row because v_0 is supported only there.
Let S be the sum of the selected nonconstant row grades, and J=s(s+1)/2.
Then

    deg c_0 <= A=3J-S-2,
    deg c_j <= 3(J-j)-S=A+2-3j,  j>=1.           (6)

In the second line expansion along v_0 deletes its constant row. Since
c_s is nonzero, its upper degree bound forces A>=3s-2, so the powers of t
used below have nonnegative exponent. The
relation sum_j c_j(z)*v_j=0 gives a scalar equation for I(t)=E(exp(tP)).
Its zeroth coefficient c_0 is nonzero. Otherwise the minimal seed equation
would have the constant solution, contradicting the retained smooth-leading
no-constant-period proof. This is where the analytic input enters; the
degree estimates alone do not establish c_0!=0.

Multiply the equation by t^A and put a_j(t)=t^A*c_j(1/t). These are
polynomials: for j>=1 the bound in(6) is at most A-1. Moreover

    deg a_0<=A,
    ord_(t=0) a_j >=3j-2,  j>=1.                 (7)

If the first K positive moments vanished, I=1+O(t^(K+1)). For j>=1,

    ord(a_j*partial_t^j(I-1))
        >=K+1-j+3j-2=K-1+2j>=K+1.              (8)

The j=0 error has order at least K+1 as well. Hence a_0=O(t^(K+1)).
Taking K=A contradicts the nonzero polynomial degree bound in(7).
There is no additional '+s' loss: the higher derivatives come with
precisely the extra vanishing orders retained in(7).

For ternary cubics, the largest possible A at each seed rank is

| s | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| A | 1 | 6 | 14 | 25 | 38 | 54 | 73 | **94** |

The bound uses the smallest available sum S of s-1 grades from(2).
Thus94 suffices uniformly, including rank drops.

In dimension N the grades have binomial multiplicities, total grade sum
N*mu/2 and constant grade zero. The maximum of A occurs at s=mu. Indeed
the sorted kth nonconstant grade is at most k (there are N degree-one
classes, and every grade is at most N); adjoining the next row increases
3s(s+1)/2-S-2 strictly. At full rank this is

    A=3mu(mu+1)/2-N*mu/2-2,

proving the first formula in(1). This is an all-N written argument;
the finite dimension table in the checker is only a regression.

## Homogeneous cubics gain a lower z-order

If P=F3, only the last two terms in(3) remain. Each replacement gains one
z and lowers x-degree by one or three. Therefore every nonzero monomial
z^r in v_(i,j) satisfies

    r <= 3j-d_i <= 3r,
    r == 3j-d_i modulo2.                         (9)

For the constant row this gives r>=j. For every nonconstant row,

    r>=j+2-d_i.                                 (10)

For d_i=1, the lower inequality in(9) gives r>=j and parity raises it
to j+1. For d_i=2 it gives r>=j. For d_i>=3, the inequality
j-d_i/3>=j+2-d_i already suffices. Rows that vanish identically cause
no problem.

It follows that the zeroth cofactor has

    ord_(z=0) c_0 >= L=J+2(s-1)-S.               (11)

Together with(6), this sharpens only the degree of a_0:

    deg_t a_0 <= A-L=s(s-1).                    (12)

The other a_j still satisfy(7), so exactly the same error-order argument
now takes K=s(s-1)<=mu(mu-1). This proves the second formula in(1),
including the ternary56 gate. If a nominal very small seed rank made
the nonzero a_0 impossible, that rank simply cannot occur in this class;
no extra constant solution is introduced.

## Exact evidence and boundaries

The new checker retains the three original Hesse model hashes and adds
an asymmetric homogeneous model with full seed rank eight. It verifies
32 complete twisted-divergence connection columns, constructs all seed
columns through j=8 via the exact derivative recurrence, checks(5) and
the homogeneous support restriction(9), and independently reduces P,P^2,P^3.
It checks all128 choices of ternary nonconstant rows and the arithmetic
leading to94 and56. The certificate stores every model connection and
seed column, not just hashes or rank observations. A rank-eight evaluation
proves full seed rank in three models; the symmetric model's rank-four
evaluation is not asserted to prove its full rational seed rank.

These models do not prove the universal filtered reduction or analytic
no-constant-period theorem. The former is the all-coefficient argument
above; the latter remains in the earlier canonical proof with its named
period-comparison and degenerate stationary-phase sources. No new analytic
theorem is imported and no high-order global moment ideal is expanded.

The sufficient94/56 bounds improve their smooth classes only. They do
not prove twenty/fourteen sufficient, settle singular-leading numerical
optimality, close coefficient-scheme or GVC questions, or provide the
expanded global radical cofactors. Both previously constructed lower-bound
witnesses remain unchanged and need not be recomputed.
