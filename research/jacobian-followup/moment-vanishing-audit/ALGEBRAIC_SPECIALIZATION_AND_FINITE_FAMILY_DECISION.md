# Algebraic witnesses and universal decision in finite coefficient families

15 November 2026. Internal constructive deductions, with independent review
outstanding. No novelty claim is made. The general holonomic, primary-ideal,
normalization and family-decision algorithms below are **written, not
implemented or executed**. The [small checker](verify_algebraic_specialization.py)
checks the contraction translation, finite-dimensional specialization
mechanism and explicit controls. It does not certify the general proofs.

## 1. Scope and conclusion

For each rational polynomial coefficient family of bounded size, one can
theoretically decide whether it contains a pure-null counterexample to
any of the following eventual-vanishing properties:

* Gaussian moments: E(P^m)=0 for m>=1, but E(QP^m) is nonzero for infinitely
  many m for some fixed polynomial Q.
* Constant-coefficient operators: Lambda^m(P^m)=0 as polynomials for m>=1,
  but Lambda^m(QP^m) is nonzero for infinitely many m.
* The full image/SIC contraction: Ecal(f^m)=0 as output polynomials for
  m>=1, but Ecal(g f^m) is nonzero for infinitely many m.

The family can have finitely many polynomial equations and inequations
over Q. Thus this covers fixed dimension and degree, fixed operator order,
homogeneous or prescribed-support families, rank conditions expressed by
minors, and finite unions of such charts. If a counterexample exists over
C, one exists with all input coefficients algebraic, with a monomial
multiplier of the same degree as a selected bad monomial. In particular,
the construction can return an exact number-field witness.

This is a decision theorem for **existence of a failure in an entire
finite family**, not a computation of that answer for an outstanding large
family. It supplies neither a classification of all safe points nor a
uniform mixed cutoff. It does not transfer a counterexample to a different
operator order or a smaller number of variables. The dense SIC(2)
bidegree-(3,3) question and the least ordinary-Laplacian power in the
retained three-variable frontier remain numerically unsettled here.

The coefficient families have a finite exact rational presentation. A
general oracle for unrepresented complex coefficients is not assumed.
All series and integrals below are formal and coefficientwise.

## 2. The full, possibly nonseparable SIC contraction is Gaussian

For multi-indices a,b define the usual polynomial-valued contraction

    Ecal(xi^a U^b)(x) = partial_U^a U^b evaluated at U=x.

Use independent circular Gaussian pairs (Z_i,W_i), with

    E(Z^a W^b) = a! if a=b, and 0 otherwise.

Taylor expansion, followed by this pairing, gives for every polynomial h

    Ecal(h)(x) = E h(Z,x+W).                                  (1)

This does not require h, f or g to factor in the xi and U variables.
In particular, with P_x(Z,W)=f(Z,x+W),

    Ecal(g f^m)(x) = E[g(Z,x+W) P_x(Z,W)^m].                  (2)

Keep x independent and work over K(x). Equality to zero over that field is
equivalent to the original output polynomial being identically zero.
Evaluation at x=0 would lose information: for f=U, all positive Gaussian
moments at x=0 vanish, whereas Ecal(f^m)=x^m.

The original all-g eventual property is equivalent to the full Gaussian
all-amplitude eventual property for P_x. One implication follows from (2).
For the other, every monomial Gaussian amplitude has the expansion

    Z^a W^b =
      sum_(j<=b) binom(b,j) (-x)^(b-j) Z^a (x+W)^j.            (3)

Thus it is a finite K(x)-linear combination of translated original
monomial multipliers xi^a U^j. Extension from K to K(x), or to C, adds only
finite linear combinations for each fixed polynomial, so cannot create
a bad multiplier if all monomials have polynomial moment series.

The [finite Gaussian multiplier algorithm](FINITE_MULTIPLIER_MIXED_DECISION.md)
therefore decides the original all-g property, produces a finite
original-multiplier degree bound, and finds the least failure degree.
If the Gaussian test bound is b, it suffices to test all original monomials
of total degree <=b: (3) never increases total degree. No separability
hypothesis is introduced. Conversely a translated original multiplier
also has Gaussian amplitude degree at most its original total degree.

The [primary finite-algebra oracle](PRIMARY_THICKENING_MOMENT_IDEAL_ALGORITHM.md)
also computes the full ordinary ideal of all output coefficients of
Ecal(f^m), and the least finite prefix generating it, in any fixed
polynomial parameter family. In a primary component extend its
normalization field by independent x; the injection A[x] -> A tensor K(x)
preserves all nilpotents. A failed identity supplies an output coefficient
outside the proposed ideal. Strict ordinary-ideal ascent terminates.

Finally suppose f is pure-null and passes the all-g test. Apply the
[Gaussian linear-cutoff theorem](LINEAR_MIXED_CUTOFF_FROM_B_FUNCTION.md)
to P_x over K(x). It gives a rational C(f)>=0 such that

    Ecal(g f^m)=0 whenever m>C(f) deg(g).                      (4)

Translation does not increase amplitude degree, so this is an original
SIC multiplier bound. Over an effective field C(f) is theoretically
computable. This proves the qualitative-to-linear implication for the full
SIC contraction in every degree and pair dimension. It does not assert
the corresponding linear statement for the Laurent-encoded arbitrary
operator problem, nor a uniform C over coefficient families.

[Zhao, 0907.3990v2](https://arxiv.org/html/0907.3990v2), section 4.1,
equations (4.3)--(4.4), gives the polynomial contraction formula; section
4.2 relates the image formulation to Mathieu subspaces. The old conjecture
status is superseded by the retained counterexamples. We use the
polynomial map, not an unrestricted substitution into a completed
power-series ring: each coefficient of our time series is a polynomial,
so (1) applies separately to it.

## 3. A bounded-space algebraic specialization lemma

Let B=Q[a_1,...,a_d]. Suppose S is any collection, even infinite, of
nonzero polynomials in B with degree in each a_i at most M. There is an
explicit algebraic point that annihilates none of S.

Choose an integer b>max(M,1), put N=b^d, choose alpha satisfying

    alpha^N=2,

and set

    a_i = alpha^(b^(i-1)).                                   (5)

The polynomial T^N-2 is Eisenstein at 2 and hence irreducible over Q.
Distinct monomials a^e with 0<=e_i<=M map to distinct powers alpha^k
with 0<=k<N, by base-b expansion. These powers are Q-linearly independent.
Consequently evaluation (5) is injective on the entire bounded polynomial
vector space, and every member of S remains nonzero simultaneously.

The same assertion holds for polynomials with additional independent
variables x or z: apply it coefficientwise in those variables. It is not
an embedding of Q(a) into a number field. A polynomial outside the stated
degree box can vanish; for example a_2-a_1^b vanishes under (5).
For d=0 the assertion is vacuous.

This mechanism handles infinitely many conditions only because they all
belong to one explicitly bounded-dimensional rational polynomial space.

## 4. Preserving nonpolynomiality of a holonomic series

Let A=Q[c]/p be an affine domain, K=Frac(A), and let

    F(t)=sum_(m>=0) f_m(x)t^m,       f_m(x) in A[x],

be an exactly computable series with a nonzero scalar annihilator
L in K(x)[t]<partial_t>. Assume F is not a polynomial over K(x).
We show how to find an algebraic point of Spec A where F remains
nonpolynomial. Any finite list of nonzero elements of A can also be
required to remain nonzero.

Write L=sum p_j(t) partial_t^j and put

    h=max_j(deg p_j-j),
    C(z)=sum_(deg p_j-j=h) lc(p_j)(z)_j.                      (6)

The falling-factorial polynomials make C nonzero. Every nonzero polynomial
solution of L has degree a nonnegative integer root of C, as proved in
the [scalar polynomiality test](FINITE_MULTIPLIER_MIXED_DECISION.md).
Compute those roots, and let B_*=max roots or -1 if there are none.
As F is nonpolynomial, coefficient search finds an index n>B_* and
a nonzero coefficient q in A of the output polynomial f_n(x).
This search terminates; alternatively the finite origin-indicial
nonpolynomiality certificate identifies a suitable coefficient.

### 4.1 A finite normalization model

Compute Noether normalization B=Q[a_1,...,a_d] subset A with A finite over
B. Its fraction field K is finite separable over Q(a).
Choose an integral primitive element y in A. Such an element can be
found by testing rational linear combinations of integral algebra
generators: the primitive-element exceptional choices are a finite
union of proper rational linear conditions, since the extension is
separable. The [primitive element theorem](https://stacks.math.columbia.edu/tag/030N)
supplies the field statement.

The monic minimal polynomial f(a,Y) of y belongs to B[Y], because B is
integrally closed and y is integral. After inverting one computable
nonzero d(a), all original coordinates lie in B_d[y], and

    A_d = B_d[Y]/(f).                                       (7)

Indeed express each generator in the fraction-field basis
1,y,...,y^(ell-1), clear its finitely many B denominators, and use
monic division by f. Enlarge d to include the discriminant of f.
Thus every algebraic specialization with d!=0 gives a reduced
degree-ell algebra, and any root of the specialized f gives a point
of Spec A on this open model.

Use multiplication matrices in this basis to compute the field norm

    G(a,x,z)=Norm_(K(x)/Q(a,x)) C(z).                        (8)

The norm is the determinant of the multiplication map, as in
[Stacks, section 9.20](https://stacks.math.columbia.edu/tag/0BIE).
It is nonzero. For every integer r,

    C(r)!=0 in K(x) iff G(a,x,r)!=0 in Q(a,x).               (9)

Clear one common denominator of G to obtain a polynomial in a,x,z.
Its degree in the a variables is bounded independently of r.

All rational coefficients of L, the multiplication matrices and (7) use
only finitely many additional denominators. Include their nonzero norms,
the norm of q, and the norms of the requested open conditions.
For any denominator involving x, retain a nonzero numerator coefficient
in x of its norm; nonvanishing of that coefficient makes it a nonzero
rational function in x after specialization. Also retain one nonzero
coefficient ensuring the highest weight h of L survives. These are
finitely many nonzero rational functions in a; clear denominators and
include their numerators and denominators.

There is now a computable M bounding the individual a-degrees of:
all these finitely many polynomials, and all coefficient polynomials
in a of the cleared G(a,x,r), for every integer r. Taking coefficients
after substituting r cannot increase those a-degrees.

### 4.2 One algebraic specialization works at every integer

Apply (5) with this M. Every nonzero G(a,x,r) remains nonzero as a
polynomial/rational function of x, simultaneously for all integers r.
All finite model conditions also survive. Pick any algebraic root beta
of f(a_special,Y). Equation (7) gives exact algebraic c coordinates.

Specialization of multiplication matrices computes multiplication by
the specialized C(r) in the finite product algebra of roots of f.
Its determinant is nonzero whenever generic C(r) was nonzero.
Thus C_special(r) is nonzero on **every** root beta in those cases.
There are no new nonnegative integer roots of the infinity indicial
polynomial; its root set is contained in the generic root set.

The specialized L is nonzero, has the same highest weight h, and
annihilates the specialized F. To see the last assertion, clear its
finite denominators first. Every coefficient of LF=0 is then an
identity in A[x], hence vanishes at every point where the chosen
denominators remain nonzero. This is an all-coefficient identity,
not a finite-jet inference.

If F_special were nonzero polynomial, (6) would force its degree <=B_*.
But q_special!=0 supplies a nonzero coefficient of index n>B_*.
This also excludes the zero polynomial. Therefore F_special is
nonpolynomial. This proves the claimed effective specialization.

When d=0, K is already a number field. Any embedding with the original
exact presentation preserves the identities and nonzero elements,
so no Kronecker step is needed.

The common mistake would be to preserve just one nonzero coefficient
of F. In the family (1+t)^a, integer a>=0 gives a polynomial with
arbitrarily high degree. Preserving a finite early jet does not preserve
nonpolynomiality. Our construction instead prevents every new integer
infinity root at once. This family is a holonomic control, not a
Gaussian moment example or a Gaussian nonconstructibility theorem.

## 5. Decide a finite family by its pure-null components

Here is a terminating construction, conditional only on the classical
effective algebraic and holonomic algorithms already imported in the
retained pointwise/ordinary-ideal proofs.

1. Express the input as a finite union of rational closed affine charts.
   Replace an inequation q(c)!=0 by an auxiliary coordinate v with
   v*q(c)-1=0. A finite Boolean combination of equations/inequations can
   be converted into finitely many such charts.

2. Compute the full ordinary pure-output coefficient ideal in each chart,
   using the [Gaussian primary algorithm](PRIMARY_THICKENING_MOMENT_IDEAL_ALGORITHM.md),
   the [operator extension](CONSTANT_COEFFICIENT_OPERATOR_DECISION.md),
   or section 2 for SIC. This terminates by the exact ideal-membership
   oracle and strict ascending chain condition, retaining primary
   nilpotents during this computation. Take its finitely many minimal
   primes p_i. If there are none, there is no pure-null member.

3. For each generic field K_i=Frac(Q[c,v]/p_i), run the appropriate
   pointwise **all-multiplier** eventual test. The output variable x
   remains independent for polynomial-valued contractions. Each test
   terminates and, on failure, gives an original monomial multiplier
   with a nonpolynomial scalar series. If a computed candidate multiplier
   is not monomial, test its finite monomial expansion to obtain one.

4. If every generic component passes, return universal safety for the
   input family. For any fixed monomial multiplier Q on a component,
   its generic moment series is a polynomial. Hence every sufficiently
   high moment/output coefficient, an original rational polynomial
   parameter expression, belongs to p_i. It therefore vanishes at
   **every** point of that component, including denominator exceptions
   of the generic holonomic certificate. Finite linear combinations
   give the statement for arbitrary complex multipliers.
   This uses the tail coefficient identities, not direct specialization
   of a generic annihilator or generator list.

5. If one generic component fails, return failure. Its finitely
   generated field embeds in C, already establishing a complex
   counterexample. For an exact algebraic witness apply section 4
   to the nonpolynomial series of the selected monomial. All pure
   coefficient identities hold in p_i, so remain zero on specialization.
   The inverse-chart equations preserve the original inequations.
   The selected bad series remains nonpolynomial. Exact algebraic
   arithmetic gives a finite number-field representation and a scalar
   indicial/coefficient certificate for its infinite nonzero tail.

In the operator case the series is
sum Lambda^m(QP^m)t^m/(m!)^2; in the other cases it is the exponential
moment series with denominator m!. These nonzero rational normalizations
do not change which coefficients vanish. All output coefficients lie
in the original affine chart algebra, as required in section 4.

The construction is finite: ordinary-ideal ascent terminates; primary
decomposition and minimal primes are effective; there are finitely many
generic all-multiplier tests; and each has a finite scalar certificate.
The potentially large algebraic witness degree in (5) is computable,
with no efficiency claim.

One need not classify the safe exceptions inside a generically failing
component to decide whether that component contains **some** failure.
Conversely a generically safe component is safe everywhere by the
coefficient argument above. These one-way facts are sufficient.
They do not imply constructibility of the pointwise safe locus.

## 6. Consequences and exact boundaries

There is a terminating theoretical truth test for GMC in every fixed
n,D family, for constant-coefficient GVC in every fixed n,operator-order,
polynomial-degree family, and for full nonseparable SIC in every fixed
pair dimension and degree/bidegree family. A fixed symbol such as Delta^k
is simply a polynomial constraint; the original differential operator is
unchanged throughout.

Any complex counterexample to one of these three formulations can be
replaced by an algebraic-coefficient counterexample with the same finite
degree/order bounds and the same selected bad monomial. For this existence
statement take p to be the ideal of rational algebraic relations among
the finitely many original coefficients. The given complex point embeds
its generic field in C. Pure identities vanish in p; a monomial is bad
because otherwise every finite linear combination would be eventually
zero. Section 4 produces the algebraic point. This is an existence
argument even when p is not available effectively from black-box complex
inputs; the finite-family algorithm computes its own relevant primes.

Finite polynomial inequations can preserve exact degree, a specified
nonzero coefficient, rank-minor data or an open chart. Equalities
preserve homogeneous/support restrictions. Thus minimal discrete
dimension/degree/order metrics for the existence of a counterexample
agree over C and over the algebraic closure of Q whenever the metric
is expressed by such finite charts. This is not a rational-coefficient
witness theorem.

The known three-variable Delta^6 failure supplies a finite upper bound
for the least power in that particular homogeneous frontier. Testing
k=1,...,6 by section 5 would in principle compute the minimum. No such
family computation is executed here, and k<=5 remains open. The same
distinction applies to dense SIC(2) bidegree (3,3): a termination proof
does not choose its answer. Retained independently replayed SIC2C4,
Dvorsky's order-three/five-variable example, the one-pair completion
obstruction and their precise nonseparability/order boundaries are
unchanged.

No single algorithm deciding an unrestricted all-dimension/all-degree
positive conjecture follows without a finite upper bound on a potential
counterexample. Searching degrees for a failure need not terminate if
none exists. Uniform Gaussian b-function charts and resonant-image
degree bounds, pointwise coefficient-locus geometry, normalized
higher-rank phase classification and independent proof review remain
separate open obligations.

[Source and assurance audit](algebraic_specialization_source_audit.json).

