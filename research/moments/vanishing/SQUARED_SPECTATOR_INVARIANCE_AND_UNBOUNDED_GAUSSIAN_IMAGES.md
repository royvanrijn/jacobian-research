# Squared spectators preserve ordinary moments and produce unbounded actual rank

16 November 2026. Internal all-degree deductions from independence, the
retained rational-image theorem and elementary differential-operator
filtrations. No novelty, independent, formal or external assurance is
claimed. This gives structural transfers and an unbounded-rank obstruction,
not a classification of every Gaussian phase.

The [checker](verify_squared_spectator_geometry.py),
[receipt](squared_spectator_geometry_verification.json) and
[source audit](squared_spectator_geometry_source_audit.json) retain exact
operator, source and connection transcriptions. General statements have
the written proofs below; their scope is not inferred from finite checks.

## 1. An exact all-ring transfer preserves the entire ordinary ideal profile

Let A be ANY commutative Q-algebra, B in A[X], and let Y be independent
standard real Gaussian variables. Fix a nonzero R in QQ[Y] and put

    C(Y)=R(Y)^2,       P(X,Y)=C(Y) B(X).

For every d,m>=0 define the ordinary A-ideal

    V_(d,m)(B)=(E[X^alpha B^m]: |alpha|<=d),
    V_(d,m)(P)=(E[X^alpha Y^beta P^m]: |alpha|+|beta|<=d).

Then EXACTLY

    V_(d,m)(P)=V_(d,m)(B), for EVERY d,m.                 (1)

Proof. Independence gives

    E[X^alpha Y^beta P^m]
      =E[X^alpha B^m] E[Y^beta C^m].                     (2)

The second factor is rational, so the left generators lie in V_(d,m)(B).
For beta=0 this rational number is E[R^(2m)]>0, including1 at m=0.
A nonzero real polynomial is nonzero on a nonempty open set of positive
Gaussian measure, so its even moments are strictly positive. They are
nonzero RATIONAL numbers and hence UNITS in every Q-algebra A. Taking
beta=0 gives the reverse inclusion. This proof needs neither a reduced
target nor a Noetherian target.

Consequences hold for arbitrary polynomial parameter families, including
their ordinary base ideals and every specialization:

* The full ordinary pure-moment ideal is unchanged, even at each prefix.
* Every nonnegative real prescribed-cutoff ideal
  (E[X^alpha B^m]:m>c|alpha|+e), c,e>=0, is unchanged.
  In particular all linear I_N and affine K_N schemes, with their complete
  nilpotent and embedded structure, are unchanged.
* Every exact last-nonzero-moment function at bounded marker degree is
  unchanged whenever it is finite. Every least cutoff of the specified
  shape, including a least real coefficient when it exists, is unchanged.
* Pure nullity and all-original-marker eventual safety are each equivalent
  for P and B over every Q-algebra. Marker coefficients in A cause no
  difficulty because they are finite linear combinations of monomials.
* The ordinary ideal of scalars making every actual marked period
  polynomial is unchanged as well, by applying the same rational units
  after multiplying every moment by that scalar.

For the cutoff assertion, a forbidden P generator has m>c(|alpha|+|beta|)+e,
which implies m>c|alpha|+e since c>=0; the reverse again uses beta=0.
For qualitative safety use (2) in the forward direction and its beta=0
unit in the reverse. No general actual-image equivalence is asserted.

There is a useful nonsquare variant. If C in QQ[Y] is any nonzero real
polynomial, safety is still equivalent. Its even moments are units. For
even m use beta-marker1, and for odd m use the FIXED polynomial marker C,
so E[C^(m+1)] is a unit. A sufficient affine bound c*degQ+e for CB gives
c*degQ+e+c*degC for B; a bound for B works unchanged for CB. Thus the
infimum asymptotic slope allowing an arbitrary finite intercept is the
same. Pure-null equivalence and exact ideal equality use the stronger
all-positive-moments condition, automatic for C=R^2.

The real/rational condition matters for reflection. An independent
holomorphic isotropic variable Z has E[Z^m]=0, and every fixed marked
Z-period is polynomial. Multiplying ANY independent B by Z is safe by
charge/independence, including unsafe B. Over a general Q-algebra the
nonzero rational scalar UNIT, rather than positivity inside that algebra,
is the reason the square transfer works.

## 2. How a finite rational pole packet transforms

Work now over an algebraically closed characteristic-zero field K; analytic
growth assertions below use K=C. Let B be pure-null and suppose its FULL
actual Gaussian period image has generic rank one. The sealed
[rank-one image theorem](RANK_ONE_IMAGES_AND_NILPOTENT_CUTOFFS.md) gives

    M_B=K[t, product_(lambda in S)(1-lambda*t)^(-1)],

where S is a finite set of DISTINCT nonzero scalars. If S is empty B is
safe. Fix ANY integer ell>=1, add one independent real Gaussian Y, and set

    P_ell=Y^(2ell) B,
    a_m=(2ell*m-1)!!=(2ell*m)!/[2^(ell*m)(ell*m)!],
    h_lambda(t)=sum_(m>=0) a_m lambda^m t^m.               (3)

Every h_lambda belongs to the actual image M_(P_ell). Indeed the element
1/(1-lambda*t) of M_B is represented by some polynomial amplitude A(t,X).
Substitution A(tY^(2ell),X) remains a polynomial amplitude and its actual
period is h_lambda. This uses the actual map, not just an annihilator.

Conversely every marked B period is a polynomial plus principal parts at
these poles. For a pole order p>=1 and an even Y-marker Y^(2k), direct
coefficient calculation gives the operator identity ON ACTUAL SERIES

    E[Y^(2k)(1-lambda*tY^(2ell))^(-p)]
      = [ product_(j=1..p-1)(theta+j)/(p-1)!
          * product_(j=0..k-1)(2ell*theta+2j+1)] h_lambda,
    theta=t D_t.                                        (4)

The empty products equal1. An odd Y-marker gives0. Polynomial parts
give polynomials in t. The factorial ratio proving (4) is
E[Y^(2ell*m+2k)]/a_m=product_(j=0..k-1)(2ell*m+2j+1).

Let

    b_ell(z)=product_(j=0..ell-1)(2ell*z+2j+1),
    A_lambda=1-lambda*t*b_ell(theta).

The exact coefficient recurrence a_(m+1)=b_ell(m)a_m proves

    A_lambda h_lambda=1,
    (D_t composed A_lambda)h_lambda=0.                   (5)

Hence K[t] lies in D_t h_lambda for every lambda, and (4) proves the FULL
global polynomial-Weyl image identity

    M_(P_ell)=sum_(lambda in S) D_t h_lambda, S nonempty. (6)

For S empty the image is K[t], by (1) and the actual-image criterion.
Equation (6) retains the polynomial lattice as a concrete D_t-span; it
does not substitute its meromorphic connection for that lattice.
Neither an actual image-kernel computation nor the minimality of a CAS
annihilator is needed to prove (3)-(6).

## 3. A self-contained weighted proof of exact rank

Put E=K((t))<theta>, with theta*a=a*theta+t*a'. For a finite operator
L=sum a_i(t)theta^i set

    w(L)=min_i(ell*ord_t(a_i)-i).

Commuting theta through a coefficient produces a derivative term of no
smaller t-valuation but one LESS theta power. Thus commutators have
strictly higher w-weight. Leading symbols multiply in the commutative
domain K[T,T^-1,xi], and

    in_w(A_lambda)=1-lambda*(2ell)^ell*T*xi^ell.          (7)

This polynomial is PRIME: its quotient ring is K[xi,xi^-1], with
T=(lambda*(2ell)^ell)^(-1)xi^(-ell). If A_lambda=L1 L2,
the product of the initial symbols is (7). Its xi-degree is ell, the
SAME as the total differential order. Since initial xi-degree never
exceeds actual order, equality holds for both factors. Primality forces
one initial symbol to be a unit and have xi-degree0, so that factor has
order0. Therefore A_lambda is irreducible over E and over K(t)<D_t>.

For lambda!=mu, the two symbols (7) are relatively prime. A homomorphism
between the simple cyclic modules for A_lambda and A_mu is represented
by R of order<ell satisfying A_lambda R=Q A_mu. Taking initial symbols
would force in(A_mu) to divide in(R), impossible by xi-degree. Thus these
irreducible differential modules are pairwise nonisomorphic, even formally.
This is the elementary Newton-weight argument, not an appeal to an
unaudited general factorization algorithm.

The scalar h_lambda is not rational. Over C its coefficients have radius
of convergence0, since log(a_m)/(m log m)->ell. Rational formal power
series at0 have positive radius. Over any characteristic-zero extension
the same nonrationality follows by descent of a finite constant-coefficient
linear recurrence to QQ, after replacing t by lambda*t.

Let V=K(t)M_(P_ell). Its quotient by K(t)1 is a quotient of the direct
sum of the simple, pairwise nonisomorphic A_lambda modules. Every summand
has NONZERO image because h_lambda is nonrational. A quotient of this
semisimple sum therefore contains each summand exactly once. Consequently

    rank_actual(M_(P_ell))=1+ell*|S|,
    V/K(t)1 = direct_sum_(lambda in S) K(t)<D_t>/(A_lambda).
                                                               (8)

An exact basis is 1 and theta^j h_lambda, 0<=j<ell. If
b_ell(z)=sum_(j=0..ell) b_j z^j, then the complete generic connection is

    theta(1)=0,
    theta(theta^j h_lambda)=theta^(j+1)h_lambda, j<ell-1,
    theta^ell h_lambda
      =(h_lambda-1)/(lambda*t*b_ell_lead)
          -sum_(j<ell)(b_j/b_ell_lead)theta^j h_lambda.   (9)

Thus the exact actual generic rank and connection, not merely an ambient
rank bound, are known in every ell and every finite pole packet.

## 4. Irregularity, essential unit and the unique formal period

At t=0 the Newton polygon of A_lambda has one positive edge from (0,0)
to (ell,1). Its slope is1/ell with horizontal length ell; its Newton
irregularity is1. The quotient in (8) therefore has positive slope1/ell
with multiplicity ell|S| and total irregularity |S|. The trivial unit
adds a slope0 block. At infinity (9) has bounded theta-matrix, so the
generic connection is regular singular there. On the punctured affine
line (9) has no poles, and (6) restricted there is exactly this connection.

The slope convention and classical interpretation are checked against
Joris van der Hoeven, *Efficient accelero-summation of holonomic functions*,
author preprint9December2006, complete physical pages17,21-22:
[author copy](https://www.texmacs.org/joris/reshol/reshol.pdf).
The published article is J. Symbolic Computation42(4),389-428(2007),
[DOI10.1016/j.jsc.2006.12.005](https://doi.org/10.1016/j.jsc.2006.12.005).
These pages were visually read because this PDF's text encoding is corrupt.
No accelero-summation algorithm, Stokes matrix or analytic boundary theorem
is imported. The irreducibility proof above is given separately in full.
Stirling's formula directly gives exact Gevrey order ell for h_lambda.

The generic extension by the unit is NONSPLIT. For any nonzero Laurent
series f, A_lambda f has the SAME initial t-term as f, because its other
term is divisible by t and theta preserves order. Hence A_lambda f=0
has no nonzero solution in K((t)). A splitting would embed one simple
A_lambda summand into the actual scalar series field, a contradiction.
The same observation shows K(t)1 is essential in V. Clearing denominators
and differentiating a nonzero polynomial then shows K[t] is essential
in the GLOBAL actual image M as well. It is its unique simple submodule.

This does not assert that the global quotient M/K[t] is semisimple or
that its complete t=0 lattice has been classified. Equation (8) classifies
the GENERIC quotient; equation (6) is the exact global span. Additional
point-lattice assertions would require their own proof.

There is nevertheless exactly ONE normalized formal period map. Since
A_lambda=1-t*(polynomial in theta) is bijective on K[[t]], the equation
A_lambda g=c has the unique solution c*h_lambda. Any D_t-map
M->K[[t]] sends1 to a constant c and every h_lambda to c*h_lambda.
Thus Hom_D(M,K[[t]]) has dimension1, even when actual rank is arbitrarily
large. When S is nonempty Hom_D(M,K[t])=0. A unique normalized formal
period or unit seed does not force actual rank one or convergence.

## 5. Four variables already realize unbounded actual rank

Use circular covariance E[ZW]=1 and independent standard real S,Y. The
retained quartic pure-null phase is

    B=W(1+Z)-(1+3Z/2+Z^2/2)S^2,
    F_Z(t)=t/(1-t),      M_B=K[t,(1-t)^-1].

These are the SEALED
[full-source and image identities](FINITE_MULTIPLIER_MIXED_DECISION.md),
not new bounded moment conjectures. Put

    P_ell=Y^(2ell) B, ell>=1.

It is pure-null in exactly FOUR Gaussian variables, of degree2ell+4.
Its full actual image is D_t h_1, of exact generic rank ell+1; the
normalized quotient is irreducible of rank ell and slope1/ell. The unit
extension is nonsplit. The LINEAR marker Z has period h_1-1 and infinitely
many nonzero coefficients, so the least failing marker degree is1.
Actual ranks and Gevrey orders are therefore unbounded with fixed
Gaussian dimension4. This is not a new minimal-dimension/degree failure:
the earlier ternary quartic remains the smaller retained counterexample.

Every finite pole packet can also be realized: take independent copies
of lambda*B for distinct nonzero lambda and ADD them. Their unit periods
are1; product markers give rational functions with just that pole union,
and each single-pole image is retained by marking only its own factor.
The full image is precisely the corresponding localization. Multiplying
that sum by Y^(2ell) realizes (8) for every finite S. No computation at
successively larger degrees is needed for either assertion.

## 6. Identical ordinary schemes and sharp nilpotent geometry

Use the already retained one-parameter quartic family

    B_e=W(1+eZ)-e(1+3eZ/2+e^2Z^2/2)S^2,
    P_(e,ell)=Y^(2ell)B_e.

Equation (1) says that the ENTIRE ordinary marked-moment ideal profile
over QQ[e] is independent of ell. Yet for e!=0 the actual generic rank
is ell+1, while at e=0 it is1. Thus even the complete collection of
ordinary vanishing/cutoff schemes, with all nilpotents, does not determine
or bound actual generic rank. The coefficient family is genuinely
polynomial; no pole specialization of a generic operator is assumed.

The sharp real affine geometry is explicit at every nilpotent depth.
For c,b>=0 let J_(c,b) be the full ordinary forbidden ideal for
m>c*degQ+b. The circular scaling identity for B_e gives every nonzero
monomial moment the e-exponent m+beta_W-alpha_Z. For a nonconstant marker
of total degree delta this exponent is >=m-delta. Consequently

    J_(c,b)=(1),                      0<=c<1,
    J_(c,b)=(e^floor(c+b)),           c>=1.               (10)

For c<1 choose arbitrarily large m=delta>b/(1-c) and marker Z^m.
Its moment is m! times a_m, a nonzero rational constant. For c>=1,
m>c*delta+b implies the e-exponent is >c+b-1 and hence >=floor(c+b).
The marker Z at m=floor(c+b)+1 attains this exponent and has value
m!*a_m*e^(m-1). This proves both inclusions, for EVERY real c,b, not a
finite fit. Pure constant markers contribute0. Extra Y-markers do not
alter the ideal by (1).

Over QQ[e]/(e^L), L>=1, the EXACT feasible region is therefore

    c>=1, b>=0, c+b>=L.                                  (11)

The least zero-intercept real/integer coefficient is L; the least common
affine integer N in m>N(degQ+1) is ceil(L/2); the asymptotic slope with
arbitrary finite intercept is1. In QQ[e], I_N=(e^N) for N>=1 as before,
and the full affine K_N=(e^(2N)). Every one of these ordinary schemes
and sharp nilpotent cutoffs is the SAME for all ell, despite unbounded
actual ranks and Gevrey orders. This extends the retained all-N proof;
it does not reopen its hash-bound artifact.

## 7. Structural boundary and credit

This supplies a valid all-ring product transfer and a general conversion
from finite rational pole defects to arbitrary-rank irregular actual
quotients. A classification cannot assume all pure failures are rational,
regular, rank one, or bounded in actual rank at fixed dimension. Nor can
the ordinary safe-ideal topology or all cutoff schemes alone recover
the connection. Essential unit and unique formal period persist in this
unbounded family.

The remaining task is still exhaustive safe normal forms/actual quotient
geometry beyond these factor architectures, including regular nonrational
possibilities and the global lattice at an irregular point. This paper
does not classify every Gaussian image or implement a general actual
image-kernel producer. The previous ordinary/least-integer gates stay closed;
both capped jobs stay terminal UNKNOWN without restarts.

The coordinator's two-O(2)-block restriction is adopted at its EXISTING
registry scopes FIRRADIAL2RINGS20261116/FIRRADIAL2NORM20261116. The
[complete internal Long/EIT audit](../../factorial/image-rigidity/LONG_EIT_AND_FC2_AUDIT.tex)
and exact Gamma bridge close phases f((X1^2+X2^2)/2,(X3^2+X4^2)/2): over
fields only0 is eventually pure or safe; over arbitrary Q-algebras the
coefficient ideal must be nilpotent and normalization is its ordinary
torsion. Long/EIT credit and the unrefereed-source assurance limits remain.
This does not cover arbitrary Gamma shapes or more blocks. Our family is
outside that radial sector: (Z partial_Z-W partial_W)P_ell is nonzero,
already through its -Y^(2ell)W term. No source restriction is being bypassed.

Gaussian independence, Wick contraction, the retained Long-type quartic,
Euler recurrences, Newton weights and Gevrey terminology are classical
ingredients; no priority claim is made for the combination. Sealed Oaku
operators, torsion witnesses and specialization thresholds are unchanged.
Sabbah/Maisonobe retain convergent rational-slope credit. Arbitrary
nonconvergent formal extension and all-ring descent remain separately
scoped, with priority unassessed. The broad goal remains active.
