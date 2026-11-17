# The full ordinary Gaussian moment ideal is computable in every fixed degree

15 November 2026. Internal constructive extension of the
[generic-point radical algorithm](GENERIC_POINT_FINITE_MOMENT_ALGORITHM.md).
Classical effective primary decomposition and algebraic holonomic integration
are imported explicitly. No novelty, independent review or formal verification
is claimed. The general algorithm below is written, not implemented or run.
The [small replay](verify_primary_moment_ideal_algorithm.py) checks one complete
linear-family ideal construction and controls that detect loss of nilpotents
and embedded components. No full cubic coefficient calculation is undertaken.

## 1. The ordinary ideal, including infinitesimal structure

Let R=Q[c_1,...,c_s], let P(c,X) be any polynomial in n Gaussian variables,
and put

    M_m(c)=E(P(c,X)^m),       I_infinity=(M_m:m>=1) in R.

**Theorem.** There is a terminating exact algorithm computing I_infinity,
a finite K with I_infinity=(M_1,...,M_K), and the least such K.
The result also holds in a specified finitely presented quotient R/H,
including a nonreduced quotient.

In particular it applies to universal polynomials of every fixed dimension
n and degree bound D. After obtaining the finite ideal, standard effective
commutative algebra computes its radical, associated primes and a rational
primary decomposition. These are computability statements, not computed
descriptions of all those objects for arbitrary n,D. A decomposition after
extension to C can be obtained by effective algebraic-number factorization
if needed; the core algorithm stays over rational polynomial rings and
rational function fields.

This goes beyond knowing the radical. The previous radical proof deliberately
did not retain nilpotents. Its scope remains correct. The present proof adds
a separate oracle over primary thickenings, then an ordinary ascending-ideal
argument with an exact stopping test. It does not infer scheme structure
from a set-theoretic result.

## 2. An injective finite-algebra representation of each primary component

Given a proper ideal I in R, compute a full rational primary decomposition

    I=intersection_(i=1)^r Q_i,       sqrt(Q_i)=p_i.             (1)

Keep embedded components. For each Q=Q_i choose an effective Noether
normalization B=Q[y_1,...,y_d] of R/p, with the y_j represented by polynomials
in R and algebraically independent modulo p. Then

    R/Q is finite over B,
    A=Frac(B) tensor_B (R/Q) is finite dimensional over F=Frac(B),
    R/Q -> A is injective.                                    (2)

Here is why both finiteness and injectivity hold. The algebra R/p is finite
over B. Each coordinate of R/p satisfies a monic B-polynomial. In R/Q the
same polynomial has value in the nilpotent ideal p/Q, so a power of that
monic polynomial vanishes. Thus every coordinate is integral over B.
There are finitely many coordinates, hence R/Q is a finite B-module.

No nonzero element of B lies in p. Since Q is p-primary, multiplication by
such an element on R/Q is injective. Localizing at B minus zero therefore
does not kill an element of R/Q. This proves (2), including its nilpotent
directions. A Groebner computation over F gives a finite vector-space basis
of A and all multiplication matrices. A is a commutative finite algebra,
not a field; it may have nonzero nilpotents.

Rational primary decomposition and Noether normalization are standard
effective characteristic-zero commutative-algebra inputs. The
[GTZ reference](https://doi.org/10.1016/S0747-7171(88)80040-3)
is the same classical foundation used by the earlier construction.
Its publisher body was not freshly proof-audited. The installed Singular
primary-decomposition library documents rational characteristic-zero input;
that producer implementation is not executed here. All decompositions in
(1) occur in R, not in an algebraic extension or quotient implementation.

## 3. A holonomic zero test over a finite algebra retains nilpotents

Let F be a rational function field over Q, let A be any specified
finite-dimensional commutative F-algebra, and let P in A[X]. Define

    G(t)=E exp(tP) in A[[t]],       Y=G-1.

There is an effective test for Y=0, returning a nonzero positive moment
in A when Y is nonzero.

To prove this, set O=F[t,X] and view V=A tensor_F O as a free O-module of
finite rank. Give it the polynomial connection

    nabla_t = partial_t + multiplication_by(P),
    nabla_i = partial_(X_i) + multiplication_by(t P_i-X_i).    (3)

These operators commute: the coefficient algebra is commutative, and the
remaining derivative commutators cancel. Thus V is a finite-rank integrable
connection. Its order filtration concentrated in degree zero has associated
graded module supported on the zero section. It is therefore holonomic.
This direct filtration argument does not require A to be reduced.

For the classical connection criterion see
[Milicic, chapter V, section 2, theorem 2.1 and corollary 2.2](https://www.math.utah.edu/~milicic/Eprints/dmodules.pdf).
Holonomicity of direct images is in chapter I, section 13, lemma 13.4 and
theorem 13.5. Equivalently use the algebraic integration construction in
[Oaku, section 3](https://www.lab.twcu.ac.jp/oaku/Oaku2013.pdf)
and its underlying module algorithms. These facts concern algebraic
modules; no convergent integral or analytic boundary hypothesis is used.

Integration in X now supplies a nonzero scalar operator
L=sum_j p_j(t) partial_t^j with p_j in F[t], and polynomial amplitudes
U_i in A[t,X] satisfying

    sum_j p_j(t) P^j
      =sum_i [partial_(X_i) U_i+(t P_i-X_i)U_i] in A[t,X].     (4)

The class of 1 in the integrated module has a nonzero scalar annihilator:
a holonomic one-variable module has finite rank after localizing in t,
and t-torsion also has a nonzero polynomial annihilator. Clearing t
denominators yields a polynomial relation. The quotient description by
total X derivatives supplies (4).

One can find (4) without a general module implementation: enumerate finite
degree boxes for p_j and for the coordinates of U_i in a fixed F-basis of
A, and solve the resulting linear equations over F. Normalize an enumerated
nonzero scalar coefficient of L to one. Existence guarantees termination.
Only scalar coefficients in F are allowed in L. An alleged nonzero
A-coefficient operator can be a zero divisor and is not enough.

Gaussian integration by parts, coefficient by coefficient, proves LG=0.
The rational function field can be extended faithfully to an algebraically
closed characteristic-zero field for the holonomicity theorem. Existence
of each finite linear identity then descends by exact linear algebra.
Construct a nonzero scalar annihilator M of Y exactly as before:

    M=L if p_0=0;
    M=(p_0 partial_t-p_0') composed L otherwise.               (5)

Write M=sum_j q_j(t) partial_t^j and set

    mu=min_j(ord_0 q_j-j),
    B(z)=sum_(ord_0 q_j-j=mu) lc_0(q_j) (z)_j in F[z].         (6)

The falling factorials make B nonzero. If a nonzero A-valued series Y has
first nonzero coefficient y_N, the first possible coefficient of MY is
B(N)y_N. Since B(N) is a scalar in the field F, B(N)!=0 would act invertibly
on the F-vector space A. Hence B(N)=0. This is precisely where a scalar
annihilator matters: a nonzero element of A need not act injectively.

The possible N are the finitely many nonnegative integer roots of B.
Compute them by clearing F-denominators and taking the gcd of the rational
z-polynomials in its parameter-monomial coefficients. Check Y through the
largest such root, in every vector-space coordinate of A. If there are no
roots, Y must be zero. This decides Y=0 effectively, even at an irregular
t=0 singularity, without treating a nilpotent as zero.

## 4. Deciding whether every moment belongs to a supplied ordinary ideal

Apply section 3 to P reduced in each finite algebra A_i from (2).
By injectivity and (1),

    all M_m in I
      iff all M_m vanish in every A_i
      iff G_i(t)=1 for every i.                              (7)

Each test terminates. A failed test returns an actual moment M_k with
nonzero image in some A_i, which proves M_k not in I. A successful collection
of tests proves ordinary membership for every positive moment. Handle
I=R immediately.

The same oracle works for finitely many marked Gaussian series
E(Q_j(c,X) exp(tP)) minus their constant terms: take their finitely many
seeds in the integrated connection and use common scalar annihilators.
This decides all-order zero identities for a specified finite set of
amplitudes. It does not decide whether every polynomial multiplier has an
eventually zero tail.

## 5. A genuine stopping algorithm for the full moment ideal

Start with I=(0) in R, or H if working in R/H. Repeat:

1. Use (7) to decide whether every positive moment belongs to I.
2. If yes, return I.
3. Otherwise obtain M_k not in I and replace I by I+(M_k).

Every unsuccessful step strictly increases I. The ascending chain condition
for R therefore forces termination. This is more than the bare assertion
that some finite prefix generates the ideal: the oracle decides the
otherwise unavailable stopping condition.

Every generator added is a moment. At termination all moments belong to I,
so I=I_infinity (or H+I_infinity in R). Keep all added indices and let K be
their maximum, with K=0 if no moment was added. Then I_K=I_infinity.
The finite target is now known. Ordinary Groebner membership tests

    I_infinity subset (M_1,...,M_k)

are monotone in k and give the least ordinary prefix by finite binary
search on 0<=k<=K. The other inclusion holds automatically. Radical
membership similarly recovers the least set-theoretic prefix if desired.
Their least values need not agree.

This supplies a computable ordinary ideal-prefix function K_ideal(n,D),
not an explicit numerical bound. Nilpotence, embedded associated components
and full scheme information can in principle be recovered after the finite
ideal is obtained. Particular unsolved coefficient schemes and their
numeric invariants remain uncomputed. No equality of two consecutive
prefix ideals is used as a stopping criterion.

## 6. Exact replay and two false shortcuts

For P=a+bX in one Gaussian variable,

    G=exp(a t+b^2 t^2/2),
    I_infinity=(a,b^2),       sqrt(I_infinity)=(a,b).

The replay starts at zero, detects M_1=a, then on (a) detects M_2=b^2.
At I=(a,b^2), retain A=Q[b]/(b^2), not its reduced quotient. The operator
partial_t has certificate U=-b in (4), since P=bX and b^2=0.
The constant term of G-1 is zero; the indicial polynomial z then proves
the whole A-valued series is zero. This certifies ordinary ideal equality,
with least ordinary prefix two. It is a tiny full loop, not a general
primary-decomposition/normalization implementation.

**Nilpotent control.** Modulo (b^3), the same P=bX has

    G-1=b^2 t^2/2 !=0.

The scalar operator partial_t^3 is valid since P^3=0. Its indicial
polynomial z(z-1)(z-2) forces a bound of two and detects the nonzero second
moment. Passing to the radical (b) would incorrectly accept this ideal.

**Embedded control.** In R=Q[a,b],

    I=(a^2,ab)=(a) intersection (a^2,b).

For the constant Gaussian polynomial P=a, all moments vanish modulo the
minimal prime (a). But in the embedded component's algebra Q[a]/(a^2),
G-1=a t is nonzero. The scalar operator partial_t^2 and its bound one
detect the missing first moment. Checking only minimal components is
therefore wrong for ordinary ideal membership.

Localizing the whole nonprimary R/I by inverting b also kills a, even
though a is nonzero in R/I. Decomposing into all primary components before
localization prevents exactly this loss. By contrast, the primary
Q=(a^3) in Q[a,b] injects into Q(b)[a]/(a^3); the replay retains its
three-dimensional multiplication table and the nonzero square.

The controls are exact identities and explicitly rejected shortcuts.
They are not failed mathematical claims or prior failed checker runs.

There is also a scalar-coefficient control: in Q[b]/(b^2), the nonzero
operator b*partial_t kills the nonzero series b*t. Its putative indicial
coefficient at index one is b, also nonzero. Multiplying by b is not
injective. The field-valued scalar requirement in sections 3--4 cannot
be dropped.

## Recent primary-decomposition source boundary

[Miller, arXiv:2609.25409v1, 21 September 2026](https://arxiv.org/html/2609.25409v1),
sections 1--3 including theorem 3.1 and algorithm 3.4, gives a nonrecursive
exponent-matrix decomposition for monomial ideals. The nilpotent and
embedded controls here are monomial examples, but general Gaussian moment
ideals need not be monomial. The paper does not replace the general
primary-decomposition input or prove a Gaussian stopping theorem. Its
Macaulay2 appendix was viewed but not executed; no theorem from it is
imported into the argument above.

## 7. What is closed, and what remains

General theoretical computability of ordinary pure-moment ideals, their
least generating prefixes, radicals and primary decompositions is now
closed by a written constructive argument, subject to independent review
of this synthesis and the stated classical algebraic inputs. This
strengthens the previous radical-only computability result without
editing its hash-bound proof.

It does not classify the geometry of every resulting ideal, run those
algorithms at full cubic size, compute any new cubic prefix or thickness,
or solve GMC safety for all multipliers. General higher-rank rational
trivial quotients, exceptional cancellation normal forms and valid
operator transfers remain structural research targets. The separate
analytic logarithmic-curve and cubic review obligations are unchanged.
[Source audit](primary_moment_ideal_source_audit.json).

