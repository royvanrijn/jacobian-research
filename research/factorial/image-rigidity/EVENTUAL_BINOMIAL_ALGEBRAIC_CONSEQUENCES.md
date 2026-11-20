# Algebraic consequences of the fixed-support eventual binomial theorem

This note uses Theorem EveryFixedOrdinaryBinomialEventualStrongFactorial in
ALL_ORDINARY_BINOMIAL_EVENTUAL_STRONG_FACTORIAL.md for each fixed support.
The complete-intersection and independent-adjunction mechanisms are classical;
their support range is extended by that eventual input. No all-order assertion
is inferred.
Fix integers a>=0,d>=1 and let N_0>=2 be a coefficient-uniform cutoff for
that support. Write Lambda_R(X^j)=j! over a commutative Q-algebra R and

    mu_m(A,B)=Lambda_R[(A X^a+B X^(a+d))^m]
      =sum_(j=0)^m binom(m,j)*(am+dj)!*A^(m-j)*B^j.

The premise says that mu_(N-1),mu_N have no common nonzero complex
coefficient pair for every N>=N_0. It follows from the actual consecutive
moment theorem, and includes the two coefficient axes.

## 1. The ordinary coefficient scheme is a sharp complete intersection

**Theorem EventualBinomialCoefficientCompleteIntersection.**
Fix ONE terminal N>=N_0. In S=Q[A,B] put

    I_N=(mu_(N-1),mu_N), J=(A,B).

The generators are nonzero homogeneous forms of degrees N-1,N. A common
nonconstant factor would have a projective zero over C and would violate
the premise. Hence their gcd is one. Since S is a UFD, multiplication by
mu_N on S/(mu_(N-1)) is injective: divisibility of h*mu_N by mu_(N-1)
forces divisibility of h. They therefore form a regular sequence.

The two graded multiplication exact sequences give

    Hilb_(S/I_N)(t)=(1-t^(N-1))*(1-t^N)/(1-t)^2
                  =(1+...+t^(N-2))*(1+...+t^(N-1)).     (A1)

Thus the quotient is an Artinian complete intersection of length N(N-1).
Its positive-degree ideal is nilpotent and its quotient by that ideal is Q,
so it is local with residue field Q. Its top nonzero homogeneous degree is 2N-3, of dimension one.
In particular the ORDINARY, not merely radical, containments are

    J^(2N-2) subset I_N,  J^(2N-3) not subset I_N.       (A2)

**Corollary SharpLargePairCoefficientNilpotence.**
Let A,B now be elements of any commutative Q-algebra R and suppose their
two actual moments at N-1,N vanish. Specialization of (A2) gives

    (A,B)^(2N-2)=0,
    [A X^a+B X^(a+d)]^(2N-2)=0 in R[X].                (A3)

Every original polynomial marked moment consequently vanishes at every
exponent at least 2N-2: the marked polynomial itself is zero before the
functional is applied. This statement allows nilpotent coefficients;
it does not force A=B=0 over a nonreduced ring.

The exponent is universally sharp for THIS terminal pair. In the
universal quotient R_N=S/I_N, the ideal J^(2N-3) is nonzero. For each r
the coefficients of [A X^a+B X^(a+d)]^r are

    binom(r,j)*A^(r-j)*B^j, 0<=j<=r.

Their exponents ar+dj are distinct, and their nonzero rational binomial
coefficients are units. Thus these coefficients generate exactly J^r.
The universal phase has nonzero power 2N-3 and zero power 2N-2.

This is the same regular-sequence/Hilbert-series mechanism already
written in PRIME_POWER_TWO_TERM_SFC.tex Section6 (with exact certificate
context in Section7). Its extended scope
here is EVERY ordinary support at every sufficiently large terminal N;
the new analytic theorem supplies precisely that larger input class.

## 2. A genuine existential reverse over coefficient rings

**Theorem ExistentialLargePairNilpotenceEquivalence.**
For a fixed commutative Q-algebra R and fixed A,B the following are equivalent:

1. The coefficient ideal (A,B) is nilpotent.
2. The phase A X^a+B X^(a+d) is nilpotent in R[X].
3. All sufficiently high factorial moments vanish.
4. There exists N>=N_0 for which the two moments N-1,N vanish.

The coefficient calculation above proves 1 iff2; 2 implies3 and 3
implies4. Statement4 implies1 by (A3). No evaluation of N_0 is required
for this equivalence. Equivalently one may replace4 by arbitrarily
large zero pairs. The same fixed-pair assertion is NOT an iff: for
R=Q[e]/(e^(N+1)), A=e,B=0, the coefficient ideal is nilpotent but the
(N-1)st moment is (a(N-1))!*e^(N-1), which is nonzero. Nilpotence need
not have occurred by a prescribed terminal order.

## 3. Independent linear adjunction preserves the terminal cutoff

**Theorem IndependentLinearAdjunctionEventualWindow.**
Let K be a characteristic-zero field, f=A X^a+B X^(a+d) nonzero, and
let Y_1,...,Y_t be new independent exponential variables. Put

    g=f+lambda_1 Y_1+...+lambda_t Y_t, lambda_j in K^*.

For any phase h let h_m=Lambda_K(h^m)/m!, with h_0=1. Independence and
Lambda(Y_j^q)=q! give the formal power-series identity

    sum_(m>=0)g_m z^m
       =[sum_(m>=0)f_m z^m]*product_(j=1)^t(1-lambda_j z)^(-1).

No analytic convergence assertion is involved. If e_r is the rth
elementary symmetric function of the lambda_j, coefficient extraction
gives the exact finite recurrence

    f_m=sum_(r=0)^min(t,m)(-1)^r e_r*g_(m-r).            (A4)

If a block of t+2 consecutive positive moments of g ending at terminal
N>=max(N_0,t+2) vanished, (A4) at m=N-1,N would force both corresponding
moments of f to vanish, contradicting the premise. Hence every such
block contains a nonzero moment, with the SAME support cutoff N_0 apart
from the necessary positive-index requirement. The conclusion is
uniform over all field coefficients, including the lambda_j. For a
genuine two-term f it is the support-count window t+2. If f has only one
nonzero coefficient, its window is one at all orders. The resulting window
is t+1 for every terminal N>=t+1; no eventual support cutoff is needed in
that monomial case.

This transfer is classical: PRIME_POWER_TWO_TERM_SFC.tex Section4
reproduces Rocks's independent-linear-adjunction construction, cited
there to Theorem3.12 of the pinned2015 author preprint. Only its new
eventual input class is extended here; no new adjunction mechanism is
claimed.

**Corollary IndependentMonomialFactors.**
Finally multiply any such phase by a monomial in further independent
variables, M=product Z_i^(c_i), with c_i>=0. Then

    Lambda[(M g)^m]=[product(c_i*m)!]*Lambda[g^m].       (A5)

The displayed rational scalar is nonzero, so exactly the same zero
windows and terminal cutoff persist. These constructions do not imply
all-order SFC, a uniform cutoff over supports, or any statement for
arbitrary dependent additions or arbitrary polynomial supports.


## 4. Proof and source boundary

The field-valued eventual theorem is used only to exclude a nonzero complex
coefficient pair and hence a common homogeneous factor over Q. The ordinary
ideal containments then follow in Q[A,B] and specialize to every commutative
Q-algebra. No field-only statement is applied directly to nilpotent coefficients.

The complete-intersection claim, its sharp exponent, the existential ring
equivalence, and the adjunction identities have written proofs above.
Finite exact controls check representative actual moment ideals, explicit
ordinary membership witnesses and a surviving top-degree functional,
nonreduced-ring examples, independent multivariate moments and scope
counterexamples. They support this proof and its source integrity, not a
formal verification or a finite-prefix inference of the universal claims.

The primary source for independent linear adjunction is Brady Rocks's
2015 author preprint, Theorem3.12 and Corollary3.13 (printed12--13),
retained as sources/rocks-2015-author-preprint.pdf. Its normalized
one-variable recurrence is reproduced in the retained prime-power
note and generalized to t independent variables in(A4). The CI/Hilbert
mechanism is already written in that retained note's Section6.
No historical-priority claim, external review or formal verification is made.
