# Exact even-moment lattices and a sixth-moment trace test

2026-11-16. HC4MYLZ0/1. ACTUAL complex polynomial PLANE Keller maps,
ALL DEGREES and ALL positive integer moment orders.

HC4 and the plane Jacobian conjecture remain open.
The remaining unproved implication is: ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, or full marked-energy integrality; OPEN.
The higher-even-moment spaces are finite projective lattices, not unknown
analytic completions. Their exact divisors separate integrability from
multiplication closure and trace perfection. The sixth-moment lattice has
an integral trace form whose defect retains EVERY ramification index.
Its automatic perfection is NOT proved. No moment subproblem replaces
the full HC4/JC2 objective.

## 1. Setup and the exact GLOBAL power argument: LZ0

Use B=C[P,Q] subset R=C[x,y], full finite S, codifferent Lambda subset R,
Omega=dP wedge dQ and

    Ram=sum_i a_i E_i, a_i=e_i-1>=0.

Keep LW0's ONE bounded target ball V meeting the image of EVERY omitted
prime, U=K^-1(V), and its global theorem R intersect L2(U)=Lambda, including
finite cyclic quotient points. Then for EVERY integer m>=1,

    T_m={r in R:integral_U |r|^(2m)<infinity}
       ={r in R:r^m in Lambda}
       =O_Z(floor(Ram/m)).                         (1)

Indeed |r|^(2m)=|r^m|^2, so the first equality with a power locus is exact.
The height-one canonical bound is m*v_E(r)>=-a_E, equivalently
v_E(r)>=-floor(a_E/m). On source divisors r is regular. Conversely that
divisorial fractional module is contained in Lambda subset R, and its
mth powers lie in Lambda by normality/reflexivity. This proves (1) at
ALL finite points. In particular there is NO missing codimension-two
integrability correction for these EVEN INTEGER exponents: it is absorbed
by the global canonical L2 statement applied to the polynomial r^m.

Every T_m is a finite reflexive S-module and B-PROJECTIVE, by LY0's surface
depth argument. Thus a locus initially defined by powers and an integral
is in fact an embedded linear finite module. No bounded-degree assumption,
normal-crossing boundary, smooth normalization or root extraction of r
is required. The inherited canonical L2 theorem retains its stated
quotient-singularity and whole-plane hypotheses.

## 2. Multiplication and the EXACT terminal lattice

Let M=1+max_i a_i, with max(empty)=0. This is the largest finite
ramification index, with M=1 in the unramified case. Then

    Lambda=T_1 contains T_2 contains T_3 contains ... contains S,
    T_m=S iff m>=M.                               (2)

If m<M, some floor(a_i/m)>0 and the coherent divisorial module attains
that pole at its generic prime, so T_m properly contains S. In particular
the last jump is T_(M-1) != T_M=S when M>=2. Since M<=d=[L:F], this
recovers LW0's sufficient cutoff L^(2d); no new field-degree bound is claimed.

The complete multiplicative statements include

    products of m elements of T_m lie in Lambda,
    T_(2m)*T_(2m) subset T_m,
    (T_m^m)^(reflexive S-hull)=O_Z(m*floor(Ram/m)). (3)

They follow either from valuations or from Holder's inequality and (1).
The last module equals Lambda exactly when m divides EVERY coefficient
a_i. This statement is about a REFLEXIVE HULL. Equality of ordinary
products with that hull at singular points is NOT asserted.

Moreover

    T_m is closed under multiplication iff T_m=S iff m>=M. (4)

The nontrivial implication uses the determinant trick: T_m is a finite
B-module containing 1, so if it is an algebra, its elements are integral
over B and lie in the full integral closure S. Conversely S is an algebra.
Thus high-moment multiplication closure is EXPECTED even under putative
failure: it eventually gives the ramified normalization, not the inverse.
Only T_1=Lambda being closed would force S=Lambda and inversion.

## 3. Trace duality of EVERY even-moment lattice

For m>=2, put b_i=floor(a_i/m). Its target trace dual is

    T_m^#=O_Z(Ram-floor(Ram/m)),
    T_m subset T_m^# subset Lambda.                (5)

The source-DVR exponents are -b_i and -(a_i-b_i), as in LY0; both are
reflexive surface modules, so the equality is GLOBAL. Trace(T_m*T_m)
is ALWAYS regular because 2b_i<=a_i. The finite projective trace Gram
map therefore has the perfect residual linking module

    Q_m=T_m^#/T_m.                                (6)

This quotient is Cohen--Macaulay of pure curve support when nonzero, with
no independent point-supported piece. At a target prime D its length,
or equivalently its trace-discriminant valuation, is EXACTLY

    sum_(E over D) f_E*(a_E-2*floor(a_E/m)),         (7)

where f_E is the residue degree. One way to verify (7) is to compare the
inclusion S subset T_m at the target DVR: the lattice index is
sum f_E*b_E, and a Gram determinant changes by its SQUARE. The original
finite discriminant exponent is sum f_E*a_E. No global basis or
orthogonal-bundle trivialization is used. Locally the trace Gram matrix
is regular; its perfection is equivalent to Q_m=0.

For m=2, Q_2 is precisely LY0's residual even-index form: its coefficient
a-2floor(a/2) is zero for even a and one for odd a. The polynomial L4
space T_2 is therefore AUTOMATICALLY the lower half-ramification lattice.
In an all-odd ramification case it is already trace-self-dual and all its
elements are L4, even under a putative failure. This does NOT give the
missing LW estimate that EVERY element of the LARGER Lambda is L4.

## 4. LZ1: a UNIFORM sixth-moment trace test sees every positive index

For m>=3 and a>0,

    a-2*floor(a/m) >= a/3 > 0.                    (8)

Consequently for EVERY integer m>=3 the finite regular trace Gram map
of T_m is perfect if and only if Ram=0. In the ACTUAL plane Keller setup,

    trace form on the polynomial L6 lattice T_3 is PERFECT
       iff the full Keller inverse exists.        (9)

The same is true for each higher even moment T_m, m>=3. The forward
implication retains KN0's finite-cover purity/inversion argument; absence
of ramification cannot be substituted for absence of every omitted
unramified prime without that argument. Conversely inversion has T_m=S=B=R.

The support of Q_3 is EXACTLY the finite branch-discriminant support,
including odd-index primes that Q_2 can lose. It is NOT claimed to equal
every component of the nonproperness curve: images losing only unramified
sheets are a different support. Although Q_3 keeps every ramified component,
its Witt class can still be zero; no metabolic-to-zero shortcut is made.

This is a field-degree-independent REGULAR trace-matrix criterion. It
does not prove its Gram determinant a unit, any target trace tensor
regularity beyond the stated quadratic integrality, or the original
full codifferent pairing. The obstruction has been represented by a
different finite projective lattice, not eliminated.

## 5. Exact finite-cover controls: all indices and all moment orders

Use the connected constant-Jacobian Laurent controls from LW, without
replaying its registered checker:

    P=u^e, Q=v/u^(e-1), Jac=e,
    S=C[u,Q], Lambda=u^(1-e)S, source=Gm times A1.

For EVERY e>=2 and m>=1,

    T_m=u^-b S, b=floor((e-1)/m).                 (10)

On the target unit bidisk, the monomial source norm is

    integral |u|^(-2m q) dV_(u,v)
        =pi^2/(e-m*q) if m*q<e,

with divergence at equality or above. For integer q this is exactly
q<=floor((e-1)/m), verifying the endpoint strictness in (1). At m>=2,
the B-basis tau_j=u^(j-b), 0<=j<e, has REGULAR trace Gram entries and

    det Gram(T_m)=(-1)^((e-1)*(e-2)/2)*e^e*
                      P^(e-1-2b).                (11)

Thus the odd-index m=2 constant-Gram controls from LY become nonperfect
again at m=3; every even index is already detected. For instance index
three has b=1 and unit Gram at m=2, but b=0 and discriminant exponent two
at m=3. These are local/finite-cover tests, NOT a reopened whole-plane
generic-degree-three or degree-four realization queue.

No fixed high finite moment of ONE element suffices in all degrees:
choose e=m+1 and r=u^-1. Then r is in T_m but is not integral over B.
Its field degree grows with m. This does not refute the cutoff using a
fixed field degree or the full-lattice perfection test (9).

All these controls have nonconstant source units, a denominator in Q and
no retained sheet at P0. They are not polynomial maps on the whole plane
and do not supply an HC4/JC2 counterexample.

## 6. Why the EVEN-INTEGER proof does not justify every real exponent

Here is a NEW actual finite normal polynomial-cover control, using the
pulled-back target volume. On the WHOLE smooth (u,v) plane set

    P=u, Q=v^5/5-u^9*v,
    F=v^4-u^9, Omega=dP wedge dQ=F du wedge dv.

The cover is finite of degree five since v satisfies
v^5-5P^9*v-5Q=0. Its finite normalization is S=C[u,v], and its
codifferent is (1/F)S, by the monogenic derivative formula (the scalar
5 is a unit). Restrict the source to F!=0 and let r=1/F. Intrinsic
integrability with respect to |Omega|^2 has r L2 locally, since
r*Omega=du wedge dv. At EVERY generic smooth point of F0, r passes
the p=3 test: |r|^3|Omega|^2=|F|^-1 is integrable in a transverse
complex coordinate.

It nevertheless FAILS p=3 at the origin. On the anisotropic regions

    rho^4 <= |u| <= 2*rho^4, |v| <= rho^9/2,

F is bounded above and below by nonzero constants times rho^36. The
Euclidean volume is a nonzero constant times rho^26. Thus the integral
of |F|^-1 on each region is at least C*rho^-10, unbounded as rho->0.
These regions avoid F0 and map into every small target neighborhood of
the origin. This proof uses scaling only; no log-canonical-threshold
formula or resolution calculation is imported.

The control is NOT an actual whole-plane Keller map: its Jacobian F
vanishes, the chosen etale source removes a singular curve, and that
curve is not smooth in the smooth normalization as the retained actual
boundary theorem requires. Its measure here is the PULLBACK target
volume, not an unweighted constant Euclidean source volume. It only
refutes extending generic-divisor Lp tests to all finite normal covers
without point analysis. No counterexample to the actual even-integer
formula or to a stronger future theorem for actual boundary pairs is
asserted. Arbitrary non-even p at actual quotient points remains outside
the scope proved here.

## 7. Attribution, verification, and the still-open multiplication step

LW0 supplies the global canonical L2 theorem. LY0 supplies trace duality
of reflexive divisorial surface modules and projectivity over the target.
KN0 supplies the actual finite ramification class and final inverse
criterion. The proof of (1) is the direct polynomial-power argument;
it newly resolves global higher EVEN moments without importing a general
Lp theorem for singular pairs. Literature on complex singularity exponents
and log canonical thresholds was consulted for context only. No such
theorem, new external classification, or primary PDF cache is adopted.

New exact controls check the all-parameter norm, floor and trace-discriminant
identities, the uniform m>=3 distinction and the finite-cover scaling
example. No registered mathematical checker is replayed. The universal
analytic/reflexive-module proofs are written, not certified by a finite
test list. Independent replay, formal verification and external review
remain absent.

The actual missing implication is unchanged: force Lambda=S, equivalently
regularity of Tr(Lambda Lambda) or the full marked-energy condition.
The exact moment filtration and sixth-moment matrix retain this obstruction
globally; neither supplies automatic trace perfection or multiplication
closure of the full Lambda. Full HC4/JC2 remains OPEN, including all branch
indices, and the original goal stays active.
