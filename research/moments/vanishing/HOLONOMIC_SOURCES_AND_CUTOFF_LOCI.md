# Holonomic sources, linear operator cutoffs and closed cutoff loci

15 November 2026. Internal all-degree deductions; independent review remains
outstanding. The general algorithms are **written, not implemented or
executed**. The [checker](verify_holonomic_source_cutoffs.py) verifies explicit
source recurrences, cone filters, nilpotent controls and examples. These
finite checks do not certify the universal arguments.

## 1. What this adds

The qualitative-to-linear implication extends from Gaussian moments and
full SIC to **every constant-coefficient operator**. If Lambda,P satisfy

    Lambda^m(P^m)=0 for all m>=1,
    for every Q, Lambda^m(QP^m)=0 for all sufficiently large m,

then there is a rational C(Lambda,P)>=0 with

    Lambda^m(QP^m)=0 whenever m>C(Lambda,P) deg Q.             (1)

The original differential order and variable count are retained. This is a
conditional bound for safe pairs, not a positive GVC theorem.

In every finite rational Gaussian, GVC or full SIC coefficient family,
the locus satisfying any **specified integer coefficient N>=0** in (1)
is an effectively computable closed set. More precisely, its complete
ordinary coefficient ideal is computable, with nilpotents retained.
This needs no specialization-valid universal b-function chart.

The whole pure-null safe locus S is an ascending countable union of these
closed sets Z_N. Over C, the following are equivalent:

    S is constructible;
    S is closed;
    one common finite linear coefficient works for every member of S;
    S=Z_N for some N.                                       (2)

This proves equivalence, not any of these assertions for the unrestricted
Gaussian or operator coefficient space. On a supplied rational family
whose pure-null members are all safe, a common coefficient is computable
from its generic components. Its least sufficient **integer** coefficient
is computable too. All general computations remain unexecuted.

The [finite-family decision proof](ALGEBRAIC_SPECIALIZATION_AND_FINITE_FAMILY_DECISION.md)
and its algebraic-witness construction are retained. No dense SIC(2)
bidegree-(3,3) answer, least Laplacian power, new cubic prefix or
operator compression is obtained here.

## 2. A general Taylor-source recurrence

Let K be a characteristic-zero field, and suppose

    F(t,s)=sum_(alpha in N^q) F_alpha(t) s^alpha

is a formal series with a supplied holonomic polynomial Weyl-module
presentation annihilating it; here s=(s_1,...,s_q) and q>=1. The coefficients
F_alpha belong to K[[t]]. The seed of interest has F_0=1.
Effective claims use an exact effective presentation of K and computable
coefficients of F.

Write E_s=sum_i s_i partial_(s_i). There is a nonzero b(z) in K[z] and a
finite normal-order operator R with

    (b(E_s)+R)F=0,
    R=sum_l c_l t^(k_l) partial_t^(j_l) s^(u_l) partial_s^(v_l),
    r_l=|u_l|-|v_l|>=1.                                    (3)

Here is the classical input and the sign convention. The algebraic Fourier
automorphism preserves the Bernstein total-degree filtration, hence module
dimension and holonomicity. In Fourier variables y use
the integration b-function identity

    b(theta_y)+R_y,     theta_y=-sum partial_(y_i)y_i,
    wt(y_i)=1, wt(partial_(y_i))=-1, wt(t)=wt(partial_t)=0,

where R_y has weight at most -1. Transform back by

    y_i -> partial_(s_i),       partial_(y_i) -> -s_i.

Then theta_y -> E_s. Every transformed lower-weight term, even after
normal ordering, has |u|-|v|>=1. Commutators preserve that difference.
This proves (3) using the already retained integration b-function theorem
and its lift. In particular, the Fourier transform applies to the
**holonomic module**, not to an asserted convergent Fourier integral.

The imported b-function statement is
[Oaku 2018, theorem 5.4 and algorithm 5.6](https://www.lab.twcu.ac.jp/oaku/MSJSI2015oakuLecture.pdf),
already audited in the
[Gaussian linear-bound proof](LINEAR_MIXED_CUTOFF_FROM_B_FUNCTION.md).
It holds for a cyclic element of a holonomic module; no special
Gaussian exponential form is required at this stage. A finite
presentation permits effective Weyl-module computation of a lift.
This general search has not been run here.

For d=|alpha|, extracting s^alpha in (3) gives

    b(d) F_alpha =
      -sum_(l:alpha>=u_l) c_l
          (alpha-u_l+v_l)!/(alpha-u_l)!
          t^(k_l) partial_t^(j_l) F_(alpha-u_l+v_l).          (4)

Factorials are products over coordinates. Every right-side index has
total degree d-r_l<d. The condition alpha>=u_l is coordinatewise;
terms failing it are zero. This is an ordinary Taylor coefficient
formula. If exponential source coefficients are used, the factorial
normalizations change scalar factors but not time degree or vanishing.

Since every u_l is nonzero, the degree-zero extraction gives b(0)F_0=0.
Thus b(0)=0 when F_0=1. Let B be the largest nonnegative integer root of b.
The finitely many coefficients |alpha|<=B determine all others through
(4). All F_alpha are polynomial in t iff these finitely many are
polynomial. The retained scalar endpoint-indicial test decides their
polynomiality over an effective K.

Put

    C_rec=max(0,max_l (k_l-j_l)/r_l).

When all base series are polynomial and F_0=1, take

    C=max(C_rec,
          max_(1<=|alpha|<=B,F_alpha!=0) deg_t F_alpha/|alpha|), (5)

omitting an empty maximum. Induction in (4) proves

    deg_t F_alpha <= C |alpha|                              (6)

for every alpha. Indeed a nonzero right-side term has degree at most
C(d-r_l)+k_l-j_l<=Cd. A finite linear combination of monomial multipliers
has the same degree bound.

More generally, if all base series are polynomial but F_0 is a polynomial
of positive degree, the same argument gives an affine bound A+C|alpha|;
the zero-offset bound specifically uses F_0=1. All constants obtained
are sufficient rational numbers, with no optimal-real-slope claim.

Thus any holonomic joint source series with constant pure seed and
polynomial time coefficients has a linear time/source-degree support
bound. Pointwise separate holonomicity without a joint presentation
is not substituted for this hypothesis.

## 3. Joint sources for the three moment formulations

For Gaussian moments use

    F(t,s)=E exp(tP(X)+s dot X).

Each coefficient is a finite polynomial Wick calculation. Algebraic
holonomic integration of the exponential connection gives a joint
holonomic presentation in (t,s). The coefficient of s^alpha is
E(X^alpha exp(tP))/alpha!. Pure nullity makes F(t,0)=1.

For full SIC, with original variables (xi,U), retain independent output x:

    F(t,r,s,x)=E exp(t f(Z,x+W)+r dot Z+s dot(x+W)).          (7)

The coefficients are Ecal(xi^a U^b f^m)(x)/(a!b!m!).
Work over K(x), using the
[nonseparable translation](ALGEBRAIC_SPECIALIZATION_AND_FINITE_FAMILY_DECISION.md).
The joint series has 2n source variables and is holonomic. It detects
original polynomial identities; x is never set to zero.

For Lambda=lambda(partial_x), set

    F(t,s,x)=sum_(m,alpha)
      Lambda^m(x^alpha P(x)^m)
          t^m s^alpha / ((m!)^2 alpha!).                    (8)

Its [Gaussian--Laurent representation](CONSTANT_COEFFICIENT_OPERATOR_DECISION.md)
is

    CT_v E exp(v lambda(Z)+(t/v)P(x+W)+s dot(x+W)).           (9)

Retain localization at v and logarithmic measure dv/v, so the v-connection
contains -1/v. Adding the polynomial source phase leaves a holonomic
rational connection in all variables. Full algebraic localization and
integration therefore give a joint holonomic presentation in (t,s)
over K(x).

Formal coefficient finiteness still holds. At fixed time degree m and
source multi-index alpha, the Laurent residue selects exactly the
m-th lambda power. For any fixed Laurent-polynomial amplitude, its
v-shift changes that index by a fixed amount; all polynomial contractions
are finite. Gaussian integration by parts and the Laurent total-derivative
identity remain valid coefficientwise with sources present. Consequently
the period map is jointly D_(t,s)-linear. No unrestricted substitution
into an undefined multivariate Laurent completion is used.

For any of these models, every original monomial multiplier has an
eventually zero sequence iff every F_alpha(t) is polynomial. Nonzero
rational factorial normalizations do not affect that assertion or its
degree. The pure hypothesis is exactly F_0=1.
Applying section 2 proves (1) in every constant-coefficient order,
dimension and degree, and recovers the Gaussian and full SIC bounds.
It also gives a second finite test using **original monomial** multiplier
degrees up to B, without needing to transfer individual Laurent module
generators back to original multipliers.

This does not make unsafe pairs safe. Retained Dvorsky and SIC
counterexamples must have a nonpolynomial base coefficient in every
valid source b-function test. No variable-coefficient extension is asserted.

## 4. Prescribing a cutoff is one holonomic zero problem

Fix an integer N>=0. For any series F(t,s) as above, form the rational
power series

    R_N(t,s)=t / ((1-t) product_i(1-s_i t^N))
            =sum_(m>N|alpha|) t^m s^alpha.                  (10)

This formula includes N=0: the product is product_i(1-s_i).
Define the coefficientwise Hadamard product

    B_N=F odot R_N.                                         (11)

Its coefficients are exactly the moments forbidden by the prescribed
bound m>N deg Q. The strict inequality is essential: equality m=N|alpha|
is allowed. At alpha=0 the filter includes every positive pure index,
so pure nullity is part of B_N=0 automatically.

Products and partial diagonals of formal D-finite series are D-finite.
Explicitly F odot G is obtained by pairing and diagonalizing the
variables of F(y)G(z). Partial diagonals are needed to retain all
independent time/source coordinates; a single complete diagonal
would test only a thin subsequence. Therefore B_N has a computable
holonomic/D-finite presentation.

The diagonal closure input is classical. We read the January 2026
revision of
[Chen--Chyzak--Ma--Zhu, 2506.22011v2](https://arxiv.org/html/2506.22011v2),
sections 2--3, especially theorems 3.1/3.2 and the structured-operator
construction and proof in 3.25--3.27. It repairs an earlier proof route.
We use closure and its finite operator constructions, with the already
imported effective holonomic operations. Its quantitative bounds are
not an end-to-end complexity bound for our moment-ideal algorithms.

To reduce the multivariate zero test to the retained scalar one, choose
independent parameters a_0,...,a_q over K and substitute

    t=a_0 u,       s_i=a_i u.                               (12)

Each coefficient of u is a finite sum. Its distinct monomials
a_0^m a^alpha record the original coefficients without cancellation.
Thus B_N=0 iff the scalar series after (12) is zero.

This substitution is effective and preserves D-finiteness. Choose a
finite rational derivative-space presentation for B_N. Every nonzero
polynomial denominator remains a nonzero rational function after (12),
because the independent a variables distinguish its original monomials.
The chain derivative is sum_i a_i partial_i. A finite-dimensional
linear-dependence computation supplies a nonzero scalar differential
annihilator over K(a)[u]. The origin-indicial zero test is valid even
when u=0 is singular. Setting all a_i=1 instead would allow cancellation
and is not an adequate zero test.

Hence a specified integer N can be tested pointwise. For a safe input
section 2 provides a finite upper bound; testing N=0,1,... up to its
ceiling computes the least sufficient integer coefficient.

## 5. The complete ordinary cutoff ideal is computable

Let R=Q[c]/J_0 be a specified affine parameter algebra for one of the
three models. Original moments/output coefficients are polynomial in c;
output x remains independent where applicable. Let I_N be the ideal
generated by all original coefficients with

    m>N|alpha|.

Then Z_N=V(I_N) is the closed prescribed-cutoff locus, including pure
nullity. We now justify **effective generation of I_N**, not just its
abstract Noetherian finite generation.

Given a proposed ordinary ideal J containing J_0, use its full primary
decomposition, including embedded components, as in the
[primary-thickening construction](PRIMARY_THICKENING_MOMENT_IDEAL_ALGORITHM.md).
Each primary quotient injects into a finite possibly nonreduced algebra
A over a normalization fraction field K. Extend by independent output
x as needed. Polynomial moment coefficients vanish in that algebra
exactly when they belong to the primary ideal.

The joint source connection is still finite-rank and flat over K(t,s,X)
with coefficients given by commuting multiplication matrices of A.
Holonomic localization/integration applies to this finite-rank module.
Apply the scalar rational Hadamard filter (10) to each of its finitely
many K-components and then (12). The resulting finitely many scalar
holonomic zero tests decide whether all forbidden coefficients vanish
in A. No radical replacement or reduction to closed points occurs.

Combining the primary components decides I_N subset J. If false,
coefficient search, known to terminate after a negative test, returns
one explicit forbidden moment/output coefficient outside J. Adjoin it
and repeat. Every failure strictly enlarges the ordinary ideal.
ACC terminates, and the terminal ideal is exactly J_0+I_N.

This proves full ordinary ideal computability for each specified N,
including nonreduced parameter constraints. It is a written algorithm;
no general Hadamard/primary implementation has been run. It does not
give a uniform numerical generating-prefix bound for all N.

As N increases, I_(N+1) subset I_N, so the closed sets Z_N increase.
This is a **descending** chain of ideals. ACC does not force it to
terminate. Even eventual stabilization of the closed sets would concern
their radicals and would not assert stabilization of scheme thickness.

## 6. Uniform bounds on any universally safe rational family

Suppose a supplied finite rational constructible family has only safe
pure-null members. Its pure-null coefficient locus is computable by the
previous ordinary-ideal algorithm, with inverse-coordinate charts for
inequations. Take its minimal prime components p_i.

The generic all-multiplier tests all pass, by the
[finite-family decision theorem](ALGEBRAIC_SPECIALIZATION_AND_FINITE_FAMILY_DECISION.md).
At each generic field, the source construction in sections 2--3 computes
a rational coefficient C_i. For every original monomial X^alpha and
every m>C_i|alpha|, its moment/output coefficient lies in p_i.
Those polynomial identities vanish at **every** complex point of the
component, including every denominator exception of the generic
b-function or annihilator certificate. There is no need to specialize
the certificate itself.

The finite maximum C=max_i C_i works for the entire pure-null family.
An empty pure locus may use C=0. For finitely many constructible charts
take their maximum. This proves a computable uniform bound on every
supplied universally safe rational family, even if generic differential
module certificates degenerate at some points.

The finite-family decision algorithm first decides whether that
hypothesis holds. If it does, the new procedure returns C.
For a closed chart of that pure-null locus with constraint ideal J, N works at
every complex point exactly when I_N subset radical(J); ordinary
scheme identities would instead require I_N subset J. Sections 4--5
make both tests effective. Testing integers up to ceil(C) therefore
finds the least sufficient integer coefficient for the given family.

This is not a uniformity theorem for the safe subset of a family
containing both safe and failing points. That subset may not be
presented by finitely many equations/inequations. Its geometry is
the separate issue addressed next.

## 7. Constructibility, closedness and uniformity are the same remaining question

Work over C in an affine finite-dimensional coefficient family X,
with pure-null safe locus S. At every safe point the source theorem
gives a finite rational coefficient. Taking its ceiling and using
the converse implication from a linear bound to eventuality yields

    S=union_(N>=0) Z_N,      Z_N subset Z_(N+1),              (13)

where each Z_N is closed relative to X.

We use the following elementary uncountability fact: a nonempty open
subset of an irreducible complex variety cannot be covered by
countably many proper closed subsets. Here is a field proof.
Collect coefficients of defining equations for the variety, the open
complement and the countable closed sets into one countable algebraically
closed subfield k of C. A generic point over k can be realized as a
complex point: choose finitely many algebraically independent elements
of C over k for a transcendence basis of k(Y), and extend the resulting
finite algebraic field embedding into algebraically closed C.
This point avoids every proper closed k-defined set and belongs to
the chosen nonempty open. Dimension zero is immediate.

Now suppose S is constructible. On every irreducible component Y of
its closure, S contains a nonempty relative open subset. That subset
is covered by the Z_N. The preceding fact forces some Z_N to contain
all of Y. Thus every component of the closure lies in S, and S is
closed. There are finitely many components, so one common maximum N
satisfies S=Z_N. Conversely such an equality implies a common coefficient
and closedness, and closedness implies constructibility.
A common real/rational finite coefficient can be enlarged to an integer,
so it gives the same condition. This proves (2).

The result is relative on a locally closed chart as well. It uses
uncountability of C, not a countability assertion over the algebraic
closure of Q. It asserts equality of point sets, not ordinary ideals.

Thus the following programmes are equivalent for a fixed coefficient
family: proving constructibility of safety, proving closedness of safety,
or proving a uniform linear coefficient on all safe points.
The earlier sufficient b-function-chart hypotheses remain a valid route,
but are no longer the only description of the issue. Computing each
closed Z_N does not decide whether their increasing union stabilizes.
Neither finite-family failure-existence decision nor Noetherian ACC
supplies such a stopping rule.

## 8. Controls and limitations

For Lambda=partial_x and P=c constant, the joint source series is

    F(t,s,x)=exp(xs) sum_(m>=0) (cst)^m/(m!)^2.

Its exact source identity is

    [E_s^2-2x s^2 partial_s-x s+x^2 s^2-cst]F=0.            (14)

Every remainder term has a positive source-degree drop, and (5) gives
C=1. This is sharp for c!=0: Q=x^d has a nonzero m=d contraction.
For c=0 every time series is constant and coefficient zero suffices.
The checker verifies the recurrence on exact polynomial source
coefficients and also tests a nonminimal b with a non-even sign control.

The general holonomic family

    H_a(t,s)=1+s*((1+t)^a-1)

has pure seed H_a(t,0)=1. All source coefficients are polynomials in t
exactly when a is a nonnegative integer. Its prescribed-N ordinary ideal is

    I_N=(a(a-1)...(a-N)),                                   (15)

because the first forbidden coefficient is binom(a,N+1) and every later
one is a polynomial multiple. Hence Z_N={0,1,...,N} increases forever,
its union is nonconstructible, and no common coefficient exists.
This is a holonomic source example, **not** a Gaussian, GVC or SIC
moment model. It shows that joint holonomicity and a constant pure seed
alone do not prove global parameter uniformity or constructibility.

The checker also retains a nonzero nilpotent forbidden coefficient in
Q[e]/(e^2) and a cancellation example showing why independent parameters
in (12) are essential. It verifies the rational filter including its
strict boundary and N=0. These are finite transcription controls, not
a universal algorithm implementation.

No unrestricted positive vanishing statement follows. The remaining
Gaussian/GVC/SIC global-safe-locus question is sharper: establish
stabilization S=Z_N, or give a genuine moment-family obstruction.
Full-image kernel/pole/irregular geometry, normalized higher-rank phase
classification, valid operator compression, concrete lower-power
frontiers and independent proof review remain separate obligations.

[Source and assurance audit](holonomic_source_cutoff_source_audit.json).

