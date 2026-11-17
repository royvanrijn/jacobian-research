# A terminating finite-moment algorithm in every fixed dimension and degree

15 November 2026. Internal constructive theorem using classical algebraic
holonomic integration and effective commutative algebra. The proof does not
assume GMC and does not use the independently unreviewed contour exclusion.
No novelty claim, independent replay, external review or formal verification
is made. The [checker](verify_generic_point_moment_algorithm.py) verifies
the indicial zero test, formal telescoping identities and an exact
positive-dimensional linear-family descent. It is not a general implementation
of the full coefficient-variety algorithm, and no full cubic coefficient
campaign is run.

## 1. The general question is computable, not merely Noetherian

Fix integers n>=1,D>=0. Write the universal degree-at-most-D polynomial

    P_c(X)=sum_(|a|<=D) c_a X^a,
    R=Q[c_a],           M_m(c)=E(P_c(X)^m).

**Theorem.** There is a terminating exact algorithm which computes:

* the radical ideal J in R defining every complex coefficient point at which
  all positive pure Gaussian moments vanish;
* a finite index K with radical(M_1,...,M_K)=J; and
* the least such K, using finite exact radical-membership tests afterwards.

This works in every fixed dimension and degree, even where GMC is false.
The same construction works for a polynomial parameter family on an
explicit closed algebraic coefficient variety over Q. It does not provide
a degree-independent K, an efficient complexity bound, ordinary moment-ideal
generation, mixed GMC classification, or numerical values beyond the
small replay stated below.

The [earlier quantifier map](ALL_DEGREE_STRUCTURE_AND_BARRIERS.md) correctly
observed that Noetherianity alone supplies no effective stopping criterion.
The additional ingredient here is an effective identity test at each generic
point, followed by dimension descent. No equality between two consecutive
prefix ideals is used as a stopping test.

## 2. Formal Gaussian integration gives an effective generic identity oracle

Let p be a prime ideal of R. Work in the computable field

    K=Frac(R/p),

with the generic polynomial P obtained by reducing its coefficients. Set

    F(t)=E exp(tP)=sum_(m>=0) M_m*t^m/m! in K[[t]],
    Psi(t,X)=tP(X)-sum_i X_i^2/2.

The exponential connection e^Psi is an algebraic holonomic module. The
algebraic integration module in the X variables is holonomic over the
t-line. Thus its cyclic class has a nonzero annihilator

    L=sum_(j=0)^r p_j(t)*partial_t^j,       p_j in K[t].      (1)

This is the imported classical input: [Oaku, section 3, definition 4,
Theorem 1 and Algorithm 3](https://www.lab.twcu.ac.jp/oaku/Oaku2013.pdf),
with [Oaku--Takayama's algebraic functor algorithms](https://www.lab.twcu.ac.jp/~oaku/OakuTakayama2001.pdf).
We use the algebraic holonomicity/integration statement, not a claim that
the original Gaussian exponential is analytically integrable. The
[author preprint](sources/oaku_integration_author_preprint.pdf) is retained.

Here is a direct way to obtain a certificate without implementing a whole
integration package. Enumerate finite degree boxes for nonzero p_j(t) and
polynomials A_i(t,X), with coefficients in K, and solve the finite linear
identity

    sum_j p_j(t)*P(X)^j
       =sum_i [partial_(X_i) A_i+(t*P_(X_i)-X_i)*A_i].       (2)

The integration ideal guarantees that some such identity exists. To see
the polynomial-certificate form, write an integration-ideal operator as
a sum of total X derivatives plus annihilators of e^Psi and apply it to
e^Psi. Divide by e^Psi. All derivatives have polynomial amplitudes; clearing
t denominators gives (2). The nonzero-operator condition can be enforced
by normalizing one nonzero coefficient of L while enumerating that choice.
Field arithmetic and exact linear algebra suffice for this search.

Coefficientwise Gaussian integration by parts makes the right side of (2)
zero after multiplication by exp(tP). Therefore LF=0 formally. This
requires neither convergence nor a vanishing analytic boundary integral.
It also does not assume isolated critical points: general holonomic
integration, not a finite Jacobi-ring assumption, supplies existence.

The field K can be represented by fractions modulo a rational Groebner
basis of p. Equality is ideal reduction of the numerator, and a nonzero
fraction is inverted by exchanging numerator and denominator. The
algebraic holonomicity statement descends to this field; equivalently,
extend it to a complex coefficient field and use faithful field extension.
The finite linear identities descend as well.

## 3. An exact finite zero test at a possibly irregular singular point

We need to decide Y=F-1, not just whether F satisfies an equation. If the
zeroth coefficient p_0 of L is zero, L already annihilates Y. Otherwise

    A=p_0*partial_t-p_0',         M=A composed with L          (3)

annihilates both F and 1, and hence Y. M is nonzero because the differential
operator ring over a field has no zero divisors. Clear t denominators if
necessary, and write M=sum_j q_j(t)*partial_t^j, with q_j in K[t].

For its nonzero coefficients put

    mu=min_j(ord_0 q_j-j),
    B(z)=sum_(ord q_j-j=mu) lc_0(q_j)*(z)_j,                  (4)

where (z)_j=z(z-1)...(z-j+1). The polynomial B is nonzero: the distinct
falling factorials in this sum are linearly independent over K.
If Y has first nonzero coefficient at index N>=0, the coefficient of
t^(N+mu) in MY is that coefficient times B(N). Later coefficients of Y
cannot contribute at this power. Thus

    Y!=0 implies B(N)=0.                                   (5)

Find all nonnegative integer roots of B in K and let Nmax be their largest.
If there are no such roots, the only formal power-series solution is zero.
Otherwise compute Y through Nmax. It is zero exactly when these finitely
many coefficients are all zero. This is a formal zero test at a singular
or irregular point; it is not an ordinary-point initial-value assumption.

Integer-root detection is effective over K. Clear the parameter denominators
of B, reduce each z coefficient modulo p, and write its normal form as a
finite sum of standard parameter monomials times polynomials in Q[z].
B(N)=0 is equivalent to all those rational polynomials vanishing at N.
Compute their gcd and its nonnegative integer roots over Q. No infinite
search through possible integer values is required.

If Y is nonzero, this oracle returns a first nonzero M_k modulo p. If it
is zero, it certifies that every M_m lies in p. The latter is a polynomial
coefficient identity on the whole variety V(p), including points at which
a denominator in a generic certificate vanishes. We do not specialize that
certificate at such a point; its zero coefficient identities already lie
in p.

## 4. Generic-point dimension descent terminates

Use effective rational prime decomposition and Groebner operations; these
are the standard commutative-algebra algorithms, e.g. Gianni--Trager--Zacharias,
[*Groebner bases and primary decomposition of polynomial ideals*](https://doi.org/10.1016/S0747-7171(88)80040-3).
Only minimal primes are required; no primary thickness is computed. Run
these decompositions in the original rational polynomial ring, not in an
algebraic extension or quotient-ring implementation that has different
termination requirements. The generic-field step uses exact linear algebra.

Start at p=(0), or at the rational prime components of the specified
coefficient variety. At each prime node:

1. Apply the exact generic identity oracle of sections 2--3.
2. If F=1 in K[[t]], accept V(p) as entirely pure-null and stop this branch.
3. Otherwise let M_k be the returned nonzero polynomial modulo p. Every
   pure-null point of V(p) lies in V(p+(M_k)). Compute its rational minimal
   primes, and recurse on them. If this variety is empty, reject the branch.

The cut is proper because M_k is nonzero in R/p. Each child therefore has
strictly smaller dimension. Every decomposition is finite and a branch has
length at most the number of coefficient parameters plus one. The finitely
branching tree with this bounded depth is finite. Each individual oracle
terminates, so the complete algorithm terminates.

All coefficients may still be arbitrary complex numbers. Rational prime
components need not remain irreducible over C: their complex zero sets
still cover the correct cut. A generic-zero certificate says M_m is in
p for every m and hence vanishes on every complex point of that component.
This is enough; geometric primality is not assumed.

Let J be the intersection of the accepted prime ideals. It defines exactly
the complex pure-null locus. Also collect every cut index into a finite
set S. A point annihilating all M_k for k in S must survive every cut along
a containing branch, so it reaches an accepted leaf. Conversely every
pure-null point annihilates these moments. Consequently

    radical(M_k:k in S)=J,
    radical(M_1,...,M_K)=J,          K=max(S).                (6)

Allow K=0 if the starting coefficient variety is already entirely null.
For a constrained initial variety, include its ideal in (6). All equalities
can be read over Q and then extended to C; the equality of zero loci and
faithful extension contract the radical equality back to Q.

## 5. Least prefix and exactly what has become computable

A finite replay record can retain each node's prime generators, its
polynomial telescoping certificate, its indicial polynomial and rational
integer-root gcd, the initial moments through Nmax, and its exact cut
decomposition. Accepted leaves carry their zero-jet certificate. These
finite objects check the generic identity and every branch transition;
no unrecorded analytic-period judgment is required by this algorithm.

Once (6) supplies a finite target J and an upper index K, determine the
least prefix by exact monotone tests

    J subset radical(M_1,...,M_k).

The reverse inclusion is already known. Test finitely many generators of
J by the usual Rabinowitsch Groebner membership construction. Binary search
on 0<=k<=K terminates and returns the least prefix. This avoids guessing
that a moment ideal has stabilized after one unchanged step.

Thus K(n,D) is computable as a function of n,D, without a known numeric
formula or GMC in that dimension. Its least value is computable too. This
does not make a uniform degree-free K possible: the retained growing-degree
families still have arbitrarily delayed first moments. Nor does (6) give
ordinary equality of the full moment ideals or their embedded components.

This result closes the general *theoretical computability* question.
It does not compute the open least cubic value in[20,94], expand I94
cofactors, certify a classification of mixed failures, or authorize a large
coefficient campaign without a structural reason. Those remain separate.

## 6. A complete small descent and a specialization failure control

For the universal two-variable linear family

    P=a+bX+cY,        v=b^2+c^2,
    F=exp(a*t+v*t^2/2),

the certified operator is L=partial_t-a-vt. A polynomial certificate for
(2) is A_X=-b,A_Y=-c. At the root, the indicial bound for F-1 is one and
the first nonzero moment is a. On the cut a=0, recomputing gives bound two
and the first nonzero moment is v. On (a,v), the operator is partial_t,
the zeroth coefficient of F-1 is zero, and the whole component is accepted.
Thus J=(a,b^2+c^2), with least prefix two. This includes nonzero isotropic
linear forms over C and is not merely the zero-polynomial leaf.

The root's bound one must not be reused on a=0. For a=0,b=1,c=0, the first
moment vanishes but the second is one. The generic operator's leading
valuation changes on the child, so the child needs its own field-level
certificate and indicial bound. The checker retains this exact mutation
control rather than silently using a generic bound on exceptional points.

## 7. The recent reduction paper and the assurance boundary

[Krason--Milewski, arXiv:2609.38221v1](https://arxiv.org/html/2609.38221v1),
sections 2 and 5, was read for finite-band reduction and its resonant indices.
Its nonzero leading recurrence coefficient and retained boundary terms
matter. The paper treats indefinite integral reductions; it is not cited
as proving this Gaussian null-locus algorithm. Its Mathematica package
was not executed. The general holonomic integration input is the classical
one cited in section 2.

The all-dimension termination proof, generic specialization logic and
minimal-prefix argument are written mathematics. The small replay checks
their algebraic components and one full tree; it does not execute the
general prime-variety/certificate enumeration algorithm. Independent review
of the new synthesis remains outstanding. Previous contour and cubic
assurance is not upgraded. [Source audit](generic_point_moment_source_audit.json).
