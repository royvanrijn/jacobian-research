# Support-uniform coefficient nilpotence and the Gamma1 radial Gaussian model

Use the support-uniform ordinary-binomial eventual theorem in
SUPPORT_UNIFORM_EVENTUAL_ORDINARY_BINOMIAL_SYNTHESIS.md. It supplies one
unevaluated integer M_*>=2 such that the actual consecutive moment
polynomials are coprime for every support a>=0,d>=1 and every terminal
M>=M_*. This note changes no canonical source or registry.

All coefficient rings below are unital COMMUTATIVE Q-algebras.
Write S=Q[U,V], J=(U,V), and for a fixed ordinary support put

    mu_m(U,V)=sum_(j=0)^m binom(m,j)(am+dj)! U^(m-j)V^j.

These are the literal factorial moments of f(Z)=U Z^a+V Z^(a+d).
The support parameters a,d are integers, a>=0,d>=1.

## 1. The genuine reverse is pointwise in the terminal order

**Theorem FixedPairUniversalCoefficientEquivalence.**
Fix ONE M>=2 and ONE support(a,d). Put D_M=2M-2.
The following predicates are equivalent:

1. mu_(M-1),mu_M have no common nonzero coefficient vector over C.
2. For every commutative Q-algebra R and A,B in R, vanishing of those
   two literal moments implies (A,B)^D_M=0.
3. Under the same universal quantifiers, vanishing of those moments
   implies [A Z^a+B Z^(a+d)]^D_M=0 in R[Z].
4. Over every reduced commutative Q-algebra, vanishing of those moments
   forces A=B=0.

These are equivalences of UNIVERSAL predicates. They do not say that,
for a fixed ring and coefficient pair, nilpotence is equivalent to the
vanishing of the prescribed two moments.

Here is the proof, including the exact exponent. Both moment forms are
nonzero homogeneous forms of degrees M-1 and M. Predicate1 implies their
gcd over Q is1: a positive-degree common homogeneous factor would have
a projective complex zero. Since Q[U,V] is a UFD, divisibility of
h*mu_M by mu_(M-1) forces divisibility of h. Hence the two forms are
a regular sequence. Its two graded multiplication exact sequences give

    Hilb_(S/I_M)(t)=(1-t^(M-1))(1-t^M)/(1-t)^2
                  =(1+...+t^(M-2))(1+...+t^(M-1)),
    I_M=(mu_(M-1),mu_M).                               (C1)

The quotient is Artinian, of Q-dimension M(M-1). Its top nonzero
degree is 2M-3 and has dimension1. Thus the ORDINARY ideal statements are

    J^(2M-2) subset I_M,   J^(2M-3) not subset I_M.       (C2)

Every coefficient evaluation S->R with the two moments zero factors
through S/I_M, proving1=>2. The coefficients of f^r are

    binom(r,j) A^(r-j)B^j at the distinct powers ar+dj.

Their nonzero rational binomial multipliers are units in R. They
therefore generate exactly (A,B)^r, proving2 iff3.
Predicate2 implies4 because a reduced ring has no nonzero nilpotents.
Finally apply either2,3 or4 to R=C. A nonzero coefficient vector
generates the unit ideal in C and a nonzero polynomial in C[Z];
neither can be nilpotent. This proves the remaining reverses.

## 2. Sharpness, uniformity and their quantifiers

**Corollary SupportUniformEventualCoefficientNilpotence.**
There is ONE unevaluated M_* such that for EVERY M>=M_*, EVERY ordinary
support(a,d), EVERY commutative Q-algebra R and A,B in R,

    mu_(M-1)(A,B)=mu_M(A,B)=0
      => (A,B)^(2M-2)=0
      => [A Z^a+B Z^(a+d)]^(2M-2)=0.                  (C3)

The same M_* works independently of the ring and coefficients.
Conversely, the universal implication in(C3) for EVERY M>=M_*
implies the support-uniform eventual field theorem by R=C and rational Bezout.
This is the genuine eventual equivalence.

For each chosen M>=M_* and each support the exponent2M-2 is sharp
as a universal ring bound for THIS pair. Take

    R_M=Q[U,V]/I_M,  A=class(U), B=class(V).

The two moments vanish by definition, while(C1)--(C2) give

    (A,B)^(2M-3) !=0,    (A,B)^(2M-2)=0.               (C4)

The distinct-exponent coefficient calculation also gives f^(2M-3)!=0
and f^(2M-2)=0. The ring has length M(M-1); it is local with residue
field Q because every positive-degree element is nilpotent.

One MAY choose the single order pair M_*-1,M_* and obtain the one
support-independent exponent2M_*-2 for phases vanishing at that pair.
One MAY NOT infer the full eventual theorem from a universal assertion
about that one pair alone: its converse supplies only that fixed-pair
coprimality predicate. The full formulation quantifies over EVERY
M>=M_*.

There is likewise no terminal-independent finite nilpotence exponent
for phases known to vanish at an unspecified arbitrarily late pair.
The quotients R_M in(C4), for arbitrarily large M, have unbounded
nilpotence index2M-2. Even the weaker assertion that all sufficiently
high moments vanish permits arbitrary coefficient nilpotence index:
take A=epsilon,B=0 in Q[epsilon]/(epsilon^L).

For a prescribed pair, coefficient nilpotence is not a converse to its
vanishing. For example in Q[epsilon]/(epsilon^(2M-2)), A=epsilon,B=0
has (A,B)^(2M-2)=0 but its (M-1)st moment is

    (a(M-1))! epsilon^(M-1) !=0.

Over an arbitrary Q-algebra, a NONNILPOTENT ordinary-binomial phase
therefore cannot have a zero pair at any M>=M_*: otherwise(C3) would
make it nilpotent. In a reduced coefficient ring every nonzero phase
is nonnilpotent, so the same cutoff works there for every nonzero phase.

## 3. The exact two-dimensional Gaussian substitution

Let X,Y be independent REAL standard normal variables, each of variance1,
and let G be their normalized joint moment functional. Define

    T=(X^2+Y^2)/2.

The precise identity is

    G(T^k)=k!   for EVERY integer k>=0.                (C5)

One exact derivation is the double-factorial convolution. Odd Gaussian
moments vanish, G(X^(2j))=(2j-1)!!, and independence gives

    G(T^k)=2^(-k) sum_(j=0)^k binom(k,j)
                      (2j-1)!!(2(k-j)-1)!!.

For formal z,

    sum_(k>=0)G(T^k)z^k/k!
       =[sum_(j>=0)binom(2j,j)(z/4)^j]^2
       =[(1-z)^(-1/2)]^2=(1-z)^(-1).

Coefficient comparison proves(C5). The convention(-1)!!=1 includes k=0.
Equivalently polar coordinates give that T has the Gamma1 exponential
law. No other dimension or Gamma shape is asserted.

The formula is a rational polynomial-moment identity, so it extends
coefficientwise to every commutative Q-algebra, without interpreting
that ring as a probability space. For every polynomial F in R[Z],

    G_R[F((X^2+Y^2)/2)]=Lambda_R[F(Z)],
    Lambda_R(Z^k)=k!.                                 (C6)

In particular the actual radial-binomial phase

    H(X,Y)=T^a(A+B T^d)

satisfies G_R(H^m)=mu_m(A,B) at the SAME order m. Consequently all
fixed-pair equivalences and the cutoff in(C3) transfer exactly to these
radial phases. This does not assert a theorem for arbitrary bivariate
Gaussian polynomials.

## 4. Coefficient ideals and sharpness survive the radial substitution

The coefficient ideal of H is exactly(A,B). Every coefficient is a
rational multiple of A or B; conversely its homogeneous degrees2a and
2(a+d) differ, and the coefficients of X^(2a),X^(2(a+d)) are respectively
2^(-a)A and2^(-(a+d))B. Both multipliers are units.

More generally the substitution R[Z]->R[X,Y], Z->T, is injective even
when R is nonreduced. If sum c_j T^j=0, the homogeneous degree2j part
has X^(2j) coefficient2^(-j)c_j, hence c_j=0.
Thus f^r=0 iff H^r=0, and the sharp universal quotients(C4) also have

    H^(2M-3)!=0,   H^(2M-2)=0.                         (C7)

If the two radial Gaussian moments at a terminal M>=M_* vanish,
H^(2M-2)=0 before a functional is applied. Every polynomial-marked
Gaussian moment G_R[Q(X,Y)H^r] consequently vanishes for r>=2M-2.
This marked conclusion uses literal nilpotence; it is not an assertion
about arbitrary nonnilpotent Gaussian phases.

## 5. Source and scope boundary

The CI/Hilbert-series and independent-variable algebra are already
written in EVENTUAL_BINOMIAL_ALGEBRAIC_CONSEQUENCES.md, Section1,
and its credited retained prime-power note. Here the new field theorem
makes the input cutoff uniform over ALL ordinary supports. The mechanism
and exponent formula are classical; no priority claim is made.

The proof is restricted to commutative Q-algebras. For instance over F2,
A=1,B=0,a=1 has its second and third factorial moments zero while(A,B)
is the unit ideal. No positive-characteristic extension follows.

The Gaussian statement is exactly variance1 in dimension2 with T divided
by2. In dimension3 the analogous variable has mean3/2, already different
from Gamma1. Dropping the divisor2 in dimension2 multiplies the kth
moment by2^k, so it changes the literal formula even though nonzero scalar
factors preserve which moments vanish.

The shared cutoff and any containing finite exceptional box remain
unevaluated. No all-order ordinary-binomial or arbitrary-Gaussian
Strong Factorial claim is obtained. Finite exact algebra, Macaulay ranks
and Gaussian controls support this written proof rather than supplying
a formal universal analytic verification.
