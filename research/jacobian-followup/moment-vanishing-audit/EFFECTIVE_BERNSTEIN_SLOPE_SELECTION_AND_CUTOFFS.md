# Effective slope selection and cutoffs from actual Bernstein certificates

16 November 2026. Internal written theorem and exact implementation.
Independent replay, formal verification and external review remain false.
No novelty is claimed. See the [producer](bernstein_slope_producer.py),
[controls](verify_bernstein_slope_producer.py),
[receipt](bernstein_slope_producer_verification.json) and
[source audit](bernstein_slope_producer_source_audit.json).

This closes the finite extraction step left by the
[rational first-survival theorem](BOUNDED_BERNSTEIN_STEPS_AND_RATIONAL_SURVIVAL_SLOPES.md).
From SUPPLIED VALID actual Bernstein certificates, one finite set of
original-power evaluations determines the actual rational slope or
eventual annihilation, decides qualitative all-marker safety, and
produces numerical affine cutoffs. On a supplied all-positive-pure
input it produces a numerical linear cutoff when safe. The extraction
works in every dimension and degree. It is not a general generator or
verifier for the underlying infinite-formal Bernstein certificates.

## 1. Exact certificate contract

Let K be characteristic zero, P a nonzero polynomial of positive
degree D, lambda a nonzero zero-constant formal symbol, and

    T=lambda(partial), tau_k=min{j>=0:T^kP^j!=0},
    v_k=lambda(xi)^k*1 in the retained ACTUAL module.

Supply a nonzero p in Q[X,Y] and integers k0>=1,J>=0 certifying:

1. For every k>=k0 with finite tau_k, p(k,tau_k)=0.
2. For every k>=k0, v_k!=0 implies v_(k+1)!=0, and then
   0<=tau_(k+1)-tau_k<=J.

The prior proof produces these data from its two Bernstein equations,
finite localization-torsion exponents, and an index past the roots of
the ordinary scalar Bernstein factor. In particular item2 includes
NONVANISHING PROPAGATION, not just a step bound conditional on both
orders being finite. It prevents a finite query from preceding an
unknown later death of the sequence.

For arbitrary K a rational component polynomial exists by the earlier
Q-linear-span argument. An effective producer must supply that
component and its validity index; the current code does not compute
a Q-basis inside a general abstract field or verify infinitely many
formal coefficient identities. Merely declaring the two properties
true in a file is not a mathematical certificate.

The postprocessor consumes their finite SUMMARY (p,k0,J), conditional
on those proven properties. It then performs exact rational algebra.
Original-power tests can be supplied by any exact effective-field
backend. The included evaluator implements Q coefficients and a
complete finite symbol jet. No purity premise is needed until the
linear-cutoff or critical-marker output is requested.

## 2. An effective barrier for each possible slope

The rational theorem says that, in the all-finite alternative, the
actual slope rho=a/b>0 in lowest terms is the slope of a rational
AFFINE LINE FACTOR of p and rho<=J. Factor p over Q and collect the
finitely many distinct such slopes in (0,J]. For each one form

    q(X,Z)=p(X,(aX+Z)/b)=sum_(i=0)^d X^i q_i(Z),
    f(Z)=q_d(Z)!=0, e=degree(f), alpha=leading coefficient(f).

All coefficients here are rational; clearing denominators is optional.
Let R be a rational Cauchy bound for the moduli of all complex roots
of f. If f is constant take R=0; otherwise take

    R=1+max_(j<e)|coefficient_(Z^j)(f)/alpha|.

The elementary triangle inequality proves the root bound. Put

    B=ceil(R)+1, L=a+bJ, A=B+L+1,
    E=sum_(i<d,j>=0)|coefficient_(X^i Z^j)(q)|*A^j,
    K_rho=max(1,floor(E/|alpha|)+1).                       (1)

For real |Z|>=B, factor f over C. Each root has modulus<=R,
so |f(Z)|>=|alpha|*(B-R)^e>=|alpha|. For real X>=1 and
|Z|<=A, the lower-X-degree terms satisfy

    |X^(-d)q(X,Z)-f(Z)|<=E/X.

Consequently for EVERY real X>=K_rho the two closed collars

    B<=Z<=B+L+1,  -B-L-1<=Z<=-B                          (2)

contain NO real root of q(X,Z). This is a certified inequality from
explicit coefficient norms, not a sampled root plot.

For the actual selected slope, z_k=b*tau_k-a*k is an integer and
|z_(k+1)-z_k|<=L. The rational theorem's finite parallel-line
conclusion places every sufficiently late z_k at a constant line
offset c with q(X,c) identically zero. Such c is a root of f, hence
|c|<=R<B. The collars (2) and their width L+1 prevent a bounded step
from crossing in EITHER direction between |z|<B and |z|>B+L+1.
No root lies inside either collar. Since the actual sequence is
eventually inside, it was inside for EVERY k>=max(k0,K_rho):

    |b*tau_k-a*k|<B.                                    (3)

This backward use of an eventual result is legitimate because the
root-free barriers hold uniformly at ALL indices after the computed
start. It does not assume that the exact affine line selected at an
early index will be selected forever. Bounded nonlinear branches
can still occur at early indices within the strip.

## 3. One computable original-order query selects the alternative

Choose a common integer K>=k0 and every K_rho. Increase it until:

* For every two distinct candidate slopes rho_i,rho_j,
  |rho_i-rho_j|*K>B_i/b_i+B_j/b_j.
* The coefficient of the highest Y power in p(K,Y) is nonzero.

The first is a finite set of rational inequalities. For the second,
choose K beyond a Cauchy bound for the roots of that nonzero
coefficient polynomial in X. Thus p(K,Y) is a nonzero univariate
rational polynomial. Its NONNEGATIVE INTEGER roots form a finite,
effectively computable list S. There may be no candidate slopes or
no roots; neither is silently replaced by an observed numerical root.

For every j in S, evaluate the ORIGINAL polynomial T^K(P^j).
Only the symbol jet through order D*j is needed. If all evaluations
are zero, tau_K=infinity: a finite first survivor would have belonged
to S and given a nonzero evaluation. Hence all later orders are
infinite, and this is an actual annihilation certificate for every
phase power at operator index K.

Otherwise let j0 be the least successful j. It equals tau_K, since
any smaller first survivor would again have been in S. Item2 of the
contract now prevents every future order from becoming infinite.
The rational theorem applies to this full finite tail, and (3)
places tau_K in the strip for its actual slope. The pairwise strip
separation makes that slope UNIQUE. Select it by

    |b*j0-a*K|<B.                                       (4)

Thus the alternative and actual rational slope are determined by
one finite collection of original-power tests. The largest required
formal-symbol jet is D*max(S), or0 when S is empty. This is not a
finite jet bound independent of the supplied certificates.

## 4. Numerical cutoffs and critical bad-marker candidates

If tau_K is infinite, every actual height h_m is at most K-1.
The formal commutator lemma gives

    m>deg(Q)+K-1 => T^m(QP^m)=0.                         (5)

Thus C=1,E=K-1 is an affine cutoff. Under supplied all-order purity,
C_linear=K is a sufficient linear coefficient, using purity for
constant markers and deg(Q)>=1 for the others.

For a selected slope rho=a/b>1, set

    H=max(K-1,ceil(B/a)),
    C=ceil(a/(a-b)), E=ceil(a*H/(a-b)).                   (6)

For h_m>=K, nonvanishing gives tau_(h_m)<=m, and (3) implies
h_m<=m/rho+B/a. For h_m<K the bound h_m<=m/rho+H is immediate.
Hence it holds globally. If m>C*deg(Q)+E, then

    (1-1/rho)m>deg(Q)+H,
    m-deg(Q)>h_m.

The same commutator lemma kills the output. Under all-positive
purity, C_linear=C+E suffices. These are numerical sufficient
bounds, not least bounds. The supplied certificate may give very
large K and therefore a large intercept.

If rho<=1, the retained rational theorem certifies failure of
qualitative all-marker safety. Supplied all-positive purity excludes
rho<1. At rho=1 collect the POSITIVE INTEGER offsets d of ALL
parallel line factors Y-X-d of p. At least one listed fixed marker
Q=P^d has infinitely many nonzero outputs. This is a finite bad-marker
candidate family, not an implemented selection of which member is
bad. The finite parallel-line theorem does not give a periodic
selection rule. The producer explicitly reports this limit.

## 5. Implemented interface and reproducibility

The implementation has three operations, also importable as functions:

    python bernstein_slope_producer.py plan input.json plan.json
    python bernstein_slope_producer.py evaluate input.json evaluations.json
    python bernstein_slope_producer.py resolve input.json result.json

The [checked example bundle](bernstein_slope_producer_example.json) can
serve as input to all three commands. It selects slope3/2 for
P=x^2,T=partial_x^3 and emits sufficient affine and pure-linear
cutoffs. The checker exercises the complete CLI as well as the API.

The plan input has polynomial_terms, validity_start and step_bound.
Polynomial terms are [[[i,j],"rational coefficient"],...], representing
X^iY^j. The plan includes every strip, its root bound and explicit
collar inequality data, the computed operator index K, the complete
integer-root candidate list, and the required jet order per phase
degree. Each operation recomputes the plan to detect inconsistent
derived data.

The rational evaluator input contains the plan, phase_terms,
symbol_jet_terms and complete_through. Both term lists use exponent
vectors of the original spatial dimension and rational coefficients.
Omitted coefficients THROUGH that order are asserted to be zero;
coefficients beyond it are not needed. Sparse multiplication and
literal repeated differentiation compute each original T^K(P^j)
exactly, including cancellations. The receipt stores nonzero results,
term counts and output hashes. It does not substitute the roots of
p for the original-power evaluations.

Resolve takes the plan and those evaluations; an optional
all_positive_pure premise enables the linear cutoff or critical
bad-marker-family output. The code does not infer that premise
from the finite query. Empty or insufficient jets, missing or
extra evaluations, altered plans, and selections outside every
unique strip are rejected. The implementation supports arbitrary
dimension and degrees over Q in principle, with no complexity claim.
The mathematical postprocessor works over any effective field with
an exact original-power evaluator and the supplied rational data.

Exact controls include genuine unary safe/unsafe inputs, an inactive
direction with immediate eventual annihilation, the retained genuine
critical unsafe DVG seed, multiple candidate slopes, extra nonlinear
branches with the same slope, and a bounded nonlinear branch. They
exercise the implementation; the universal proof is Sections1-4.

For the DVG control, P=(t+c)(ad+bt), R=partial_a partial_d-
partial_b partial_c and T=partial_t R. If j<k, R^kP^j=0 by
degree in (a,b); at j=k the retained identity is
R^kP^k=(k!)^2*c^k. At j=k+1, direct binomial expansion gives

    R^kP^(k+1)=((k+1)!)^2*ad*c^k*(t+c)
       +(k+1)(k!)^2*b*t*(c^(k+1)-(-1)^(k+1)t^(k+1)).

After partial_t^k the coefficient of b*t^2 is
(-1)^k*(k+1)(k!)^2*(k+2)!/2, nonzero. Thus tau_k=k+1
for EVERY k>=1, not by extrapolating the finite control.

The checker also retains a deliberately FALSE certificate p=Y-X-1
for P=x,T=partial_x^2. It passes the single original query at K=1
but predicts critical slope1 instead of the actual slope2. This
confirms an essential boundary: finite original-power queries do not
verify the supplied all-index certificate contract. Results are
conditional on its actual proof, and the program labels them so.

## 6. Bare coefficient queries cannot supply the missing classification

There is no algorithm which, for EVERY nonzero zero-constant formal
symbol given only by coefficient queries, halts after finitely many
queries and correctly distinguishes eventual infinite survival orders
from the finite rational-slope alternative. This remains true on
binary, degree-one, ALL-POSITIVE-PURE SAFE inputs.

Indeed take P=x and first run the alleged algorithm on lambda=zeta,
the inactive y-direction. All tau_k with k>=1 are infinite. A
halting execution queries only finitely many exponent vectors.
Choose an integer M>=2 whose vector (M,0) was not queried, and
instead supply lambda_M=zeta+xi^M. Its replies to every query in
that execution are identical, including adaptive queries. Yet
T_M restricted to K[P] is partial_x^M, so tau_k=M*k is finite
with slope M. The two symbols are both nonzero, zero-constant
and polynomial, and both inputs are pure and safe. The algorithm
must give the same answer on them, a contradiction.

This is a representation/oracle obstruction. A finite polynomial
expression, a known support bound, or actual Bernstein certificates
provide more information than that bare oracle. It does not prove
that safety itself is undecidable from coefficient access, or forbid
a common sufficient cutoff for this family: m>2deg(Q) follows from
h_m<=m/2 for lambda_M and h_m=0 for the inactive symbol. No such
broader negative claim is made.

## 7. Ownership and remaining general obligation

This thread owns the GENERAL GAUSSIAN and FORMAL CONSTANT-COEFFICIENT
programme, under the user's focused ownership amendment. Retain the
broader objective: read relevant research, adopt useful methods, update
knowledge/open questions, and pursue general classifications, full
theorems and valid counterexamples. The prior field, all-ring and
rational-slope closures stay proved at their exact scopes, with
their existing assurance limits and countercontrols preserved.

The next implementation obligation is to construct effective ACTUAL
Bernstein certificates, their finite torsion bounds and valid
specialization indices for general effective input classes, rather
than merely postprocess supplied data. For arbitrary infinite formal
symbols, neither an effective presentation/oracle nor a certificate
generator is automatic. This theorem does not remove those missing
inputs or implement a least cutoff. General ordinary-ideal algorithms,
complete Gaussian safe normal forms and other retained structural
questions remain part of the wider programme.

The separate Mathieu thread owns direct proof/counterexample work on
the full exceptional SL2/2O and SL2/2I MZ cases. Record concrete
transfers when available, but do not replace this general programme
with parallel direct MZ architecture searches. No Mathieu conclusion
follows from these formal certificates alone. No external outreach,
broad compute campaign, all-ring slope promotion, periodicity theorem
or independent proof assurance is claimed.
