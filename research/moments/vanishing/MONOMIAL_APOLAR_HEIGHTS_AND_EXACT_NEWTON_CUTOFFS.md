# Monomial apolar heights and exact Newton cutoffs

16 November 2026. Internal written proofs and small exact transcription
controls; no novelty, independent replay, external review or formal
verification is claimed. See the [checker](verify_monomial_newton.py),
[receipt](monomial_newton_verification.json) and
[source audit](monomial_newton_source_audit.json).

The finite-symbol monomial-phase GVC theorem is classical:
van den Essen--Willems--Zhao,
[0903.1478v2](https://arxiv.org/html/0903.1478v2), Theorem5.1 and
Corollary5.3. We credit it explicitly. Here the symbol may be an ARBITRARY
formal power series. We prove an exact optimal multiplier cutoff and a
rational, bounded-error height law for the actual rectangular apolar
quotients. These do not prove the gap implication for general phases.

## 1. Finite Newton data without truncating a formal symbol

Let K be a characteristic-zero field, let a be in N^n, P=x^a, and
T=lambda(partial), where lambda is a nonzero element of K[[Z]]. Put

    N(lambda)=conv(Supp(lambda))+R^n_(>=0).

Dickson's lemma gives a FINITE set M of coordinatewise minimal support
exponents, and

    N(lambda)=conv(M)+R^n_(>=0).                            (1)

Thus this is a closed rational polyhedron, even for a transcendental formal
symbol. This is a mathematical finite description, NOT an algorithm to
recognize the full support of an arbitrary series from a finite jet.

We will need the following elementary compact-face observation. Every
coordinatewise minimal point b of N(lambda) belongs to a compact face
exposed by a strictly positive rational weight w. Indeed, all facet
inequalities of this upward-closed polyhedron have nonnegative inward
normals. Sum the normals of the facets active at b. Every coordinate of
the sum is positive: otherwise a sufficiently small move b-epsilon e_i
would preserve the active inequalities and, by their strict slack, all
inactive inequalities, contradicting minimality. This sum exposes the
minimal face containing b; its positive coordinates make that face compact.
For rational b the normals and weight can be chosen rational.

Let f be lambda's terms on this compact face. It is a FINITE polynomial.
For every k with kb integral,

    [Z^(kb)]lambda^k=[Z^(kb)]f^k.                          (2)

Every factor contributing to the indicated coefficient must attain the
minimum w-weight. In particular, later formal terms cannot change (2).

If b is rational, choose N>=1 with beta=Nb integral. The Laurent polynomial

    L=Z^(-beta) f^N

has0 in its Newton polytope. The Duistermaat--van der Kallen theorem
(VWZ Theorem4.5) implies CT(L^j) is nonzero for infinitely many j. Otherwise
choose J beyond all nonzero indices and apply the theorem to L^J: all its
positive power constant terms would vanish although its polytope contains0.
Consequently (2) is nonzero at k=Nj for infinitely many j.

The theorem is used ONLY for this finite Laurent polynomial. It is never
applied directly to a Laurent formal series. Its characteristic-zero field
version follows by embedding the finitely generated field of its finitely
many coefficients into C and reflecting coefficient nonvanishing.

## 2. Full formal purity and the exact least real cutoff

Define the rational number

    kappa(lambda,a)=min_(b in N(lambda)) sum_i max(b_i-a_i,0). (3)

The minimum is attained at a rational point: it is a rational linear
programme with auxiliary variables u_i>=b_i-a_i,u_i>=0, and its bounded
sublevel sets in b are compact. Choose a minimizer that also minimizes
sum b_i. It is coordinatewise minimal, so section1 applies to it.

**GVCMONOMIALNEWTON20261116.** The following are equivalent:

* T^m(P^m)=0 for EVERY positive m;
* T^m(P^m)=0 for all sufficiently large m;
* N(lambda) has empty intersection with the box [0,a];
* kappa(lambda,a)>0.

When these hold, the EXACT least nonnegative REAL constant in

    T^m(Q P^m)=0 whenever m>C deg Q                         (4)

for all polynomial Q and all positive m is

    C_min=1/kappa(lambda,a).                               (5)

Here a_i may be0, so unused variables and the constant monomial phase are
included. Multiplication of P by a nonzero scalar, affine translation, and
an invertible linear coordinate change give the corresponding result for
products of powers of independent affine linear forms. Degrees are
preserved by these changes. The zero symbol is a separate trivial case:
all positive outputs vanish, with C_min=0.

Proof. For any polynomial monomial x^v,

    lambda(partial)^k x^v !=0
       iff some gamma in Supp(lambda^k) satisfies gamma<=v. (6)

Distinct gamma give distinct output monomials, with nonzero factorial
coefficients. This is an output-polynomial assertion, not evaluation at0.
Each coefficient in (6) uses a finite jet, so it holds for formal lambda.

If kappa>0, every gamma in Supp(lambda^m) lies in m N(lambda), and hence

    sum_i max(gamma_i-ma_i,0)>=m kappa.                     (7)

If its derivative survives on x^(ma+v), then the left side is at most |v|.
For a general marker Q a nonzero output requires at least one surviving
monomial of Q. Thus deg Q>=m kappa, proving (4)-(5) as an upper bound and
all positive pure outputs zero.

If kappa=0, choose the rational minimal point b<=a just described.
Section1 supplies infinitely many m=Nj with a nonzero coefficient at mb.
By (6), T^m P^m is nonzero at all those indices. This proves both necessity
and the equivalence with EVENTUAL pure nullity.

For sharpness when kappa>0, use instead a rational minimizing b in(3).
At infinitely many m=Nj its coefficient at mb is nonzero. Take the ORIGINAL
polynomial marker

    Q_m=x^(m(b-a)_+),    deg Q_m=m kappa.                   (8)

Its exponents are integral. The derivative at gamma=mb survives on Q_mP^m,
and distinct derivatives cannot cancel its output monomial. Therefore every
C<1/kappa fails (4). Markers in this sharpness argument depend on m, as is
appropriate to an optimal COMMON degree bound. They are NOT one fixed
marker witnessing individual-eventuality failure. This completes the proof.

For a support-separated pair, the same sufficient cutoff proof (7) works
over EVERY commutative Q-algebra if its symbol support lies in the specified
polyhedron. Necessity, optimality and the height assertions below are FIELD
claims. No full ordinary nilpotent coefficient-family gate is inferred.

## 3. A bounded-error rational law for the actual apolar heights

Now assume a_i>=1 for every i, lambda(0)=0 and lambda!=0. No purity
hypothesis is required. Set

    h_m=max{k>=0:lambda(partial)^k P^m!=0},
    c=min{t>=0:t a in N(lambda)}
     =min_(b in N(lambda)) max_i b_i/a_i.                   (9)

Then 0<c<infinity is rational. The actual apolar ideal is

    Ann(P^m)=(Z_1^(ma_1+1),...,Z_n^(ma_n+1)).               (10)

**GVCMONOMIALHEIGHT20261116.** There is a finite constant B depending on
the specified pair such that, for EVERY m>=1,

    m/c-B <= h_m <= floor(m/c).                            (11)

In particular h_m=m/c+O(1), with a RATIONAL slope. Its bounded error also
implies a bounded quasi-additivity error for this actual family, but not
for general polynomial phases. No numerical B is asserted or computed.

The upper bound follows from(6): gamma<=ma with gamma in kN(lambda)
implies (m/k)a in N(lambda), hence k<=m/c.

For the lower bound choose a rational b<=ca in N(lambda), minimizing sum
b_i among such points. It is coordinatewise minimal. Use its compact face
f, N,beta and finite Laurent L from section1. We need a bounded-gap version
of the infinitely many nonzero constant terms there.

The series A(s)=sum_(j>=0) CT(L^j)s^j is D-finite. For completeness this
uses a PROPER ordinary rational diagonal. Choose an integer q>=1 such that
H(Z)=(product Z_i)^q L(Z) is polynomial. The complete diagonal, in t and
all Z_i, of the ordinary rational power series

    1/(1-t^q H(Z))

is exactly A(u^q): only equal exponents qj contribute, with coefficient
[Z_1^(qj)...Z_n^(qj)]H^j=CT(L^j). Classical diagonal closure therefore
gives D-finiteness of A(u^q), and algebraic substitution u=s^(1/q) and
finite-field descent give D-finiteness of A(s). Only closure, not any
quantitative bound, is used from Chen--Chyzak--Ma--Zhu,
[2506.22011v2](https://arxiv.org/html/2506.22011v2), Theorems3.1-3.2.
The retained algebraic-symbol note already audits this classical input.

A univariate D-finite series has a nonzero polynomial-coefficient recurrence
for its coefficients. Beyond the finitely many integer zeros of its leading
coefficient, a block of r consecutive zeros propagates forever, where r is
the recurrence width. Section1 excludes an eventually zero sequence. Thus
CT(L^j) has bounded consecutive zero gaps after a finite initial segment.

Choose such a nonzero index j at most m/(Nc), within a bounded distance
below that number. At k=Nj, equation(2) gives a nonzero coefficient at
kb<=ma, so h_m>=Nj>=m/c-B. Increasing B handles the finite initial m.
This proves(11), without assuming a recurrence for the full formal source.

When lambda is pure, section2 gives c>1, so its actual deficit satisfies
b_m=m-h_m=(1-1/c)m+O(1), with strictly positive slope. Thus an unbounded
but sublinear deficit cannot occur in this entire class. If some a_i=0,
restrict lambda to the active dual coordinates before applying(9)-(11).
If that restriction is0, h_m=0 identically. The exact mixed cutoff(5)
still uses the FULL symbol, since markers can involve passive variables.

## 4. A concise ternary, genuinely nonalgebraic control

Let phi(Z)=sum_(k>=3) Z^(k!), P=xyz, and

    lambda=X^2+Y^4-phi(Z).                                 (12)

Its Newton polyhedron has compact intercept face (2,0,0),(0,4,0),(0,0,6).
The weight(6,3,2) has symbol minimum12 and phase weight11. In(3),
b=(7/6,1,1) gives kappa<=1/6. The normalized weight(1,1/2,1/3)
gives the reverse inequality, so kappa=1/6 and C_min=6.
For every m=12l, the coefficient at (14l,12l,12l) is

    (12l)!/((7l)!(3l)!(2l)!) !=0,

coming only from the compact face. The marker is Q=x^(2l), attaining the
boundary m=6 deg Q. Moreover,

    h_m=floor(m/2)+floor(m/4)+floor(m/6),  slope=11/12.      (13)

Coordinate orders bound the three allocations from above; their three
lowest terms attain all three bounds simultaneously with a unique nonzero
coefficient. The infinitely many later coefficients of phi cannot alter it.

This principal formal ideal contains NO nonzero algebraic formal series.
Here is a boundary proof, rather than an inference from nonalgebraicity of
lambda itself. If an algebraic nonzero A belonged to(lambda), its algebraic
equation has a nonzero polynomial constant coefficient B(X,Y,Z), and
lambda divides B. Divide B by the monic polynomial X^2+Y^4-U, with U an
independent indeterminate. The remainder r_0(Y,Z,U)+X r_1(Y,Z,U) is not
identically0: this monic relation has zero contraction to K[X,Y,Z]. Formal
Weierstrass division and specialization U=phi(Z) imply
r_i(Y,Z,phi(Z))=0. A nonzero coefficient of a power of Y in a nonzero r_i
is then a nonzero polynomial over K[Z] annihilating phi, a contradiction.
The lacunary phi is transcendental by the retained gap-recurrence argument
in [FORMAL_SYMBOL_SOURCE_RIGIDITY.md](FORMAL_SYMBOL_SOURCE_RIGIDITY.md):
a P-recursive sequence cannot have the arbitrarily long internal zero
blocks followed by its infinitely many coefficients1.

Thus(12) is outside unit-times-algebraic principal symbols. Also
J_P=intersection_j Ann(P^j)=0, by(10); the phase-power quotient cannot
replace its symbol by a different representative. The full raw source is
nonholonomic by retained GVCJOINTSOURCERIGID20261116. The finite compact-face
period used above remains D-finite. This distinguishes their roles.
Weighted degree already proves safety of this example; it is NOT a new
unsafe seed or a claim of novelty for a safe family.

## 5. Source boundaries and the remaining general gap

VWZ Example3.5 has the NONPOLYNOMIAL phase x+exp(y), with a finite
operator. Example4.11 has a Laurent formal series y^(-1)(1+x^(-1)exp(y)).
Neither justifies applying the finite Laurent theorem directly to an
infinite symbol, and neither contradicts the compact-face argument(2).
The source's generated HTML date differs from the arXiv version history;
v2 is dated16January2011 and the journal publication is2015.

The proof here succeeds because(10) is an ACTUAL monomial complete
intersection at every m, so output monomials cannot cancel in(6). A general
polynomial P has different annihilators, and cancellations between its
different monomials invalidate that step. Neither a monomial degeneration
of Ap(P) nor the symmetric-power subquotient comparison identifies all the
embedded Ann(P^m) with(10). No general bounded-error height law, bounded
quasi-subadditivity, finite Rees presentation, or formal eventuality-to-linear
theorem is claimed. Both full exceptional MZ cases also remain open.

The checker only checks the displayed finite identities and exact LP primal/
dual certificates, a sparse constant-term recurrence, unit and passive-variable
boundaries, and guard refusals. No increasing rank, degree or prime campaign,
general recurrence producer, or recognition of arbitrary formal Newton support
has been implemented. Closure-first work and deferred outreach persist.
