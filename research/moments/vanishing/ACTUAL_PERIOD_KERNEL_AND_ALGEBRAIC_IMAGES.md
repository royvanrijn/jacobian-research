# Computing the actual period kernel and testing algebraic images

15 November 2026. Internal constructive synthesis; independent review
outstanding. General minimization, connection localization, Weyl-kernel
and algebraicity algorithms are **written, not implemented or executed**.
The [checker](verify_actual_period_kernel.py) verifies small exact cyclic,
minimization, torsion, lattice and algebraic-cover controls. It does not
certify the general constructions.

## 1. The actual image can be recovered, including its polynomial lattice

Let K be an exactly presented effective finitely generated characteristic-
zero field. Let D=K[t]<partial_t>, let H be a finitely presented holonomic
D-module, and let

    lambda:H -> K[[t]]

be a D-linear formal period map with computable coefficients on the
specified generators. There is a terminating theoretical construction of:

* the full polynomial-Weyl submodule ker(lambda);
* a finite presentation of the actual image M=H/ker(lambda), with the
  images of the original generators;
* its exact generic rank, not just a source-rank upper bound;
* when a specified class maps to 1, the canonical obstruction module

      C=M/D[1]=M/K[t].                                      (1)

In the pure Gaussian, full SIC and constant-coefficient GVC settings,
C=0 is exactly the original all-multiplier eventual property by the
retained full-amplitude transfer theorems. If the actual rank is d,
C has generic rank d-1. At d=1 its only obstructions are the finite
point quotients from the
[rank-one classification](RANK_ONE_IMAGES_AND_NILPOTENT_CUTOFFS.md).

This is pointwise computability, not a closed-form classification of
higher-rank Gaussian phases. The actual kernel must be recomputed at
exceptional coefficient points. No constructibility or uniformity theorem
for its parameter variation is asserted.

For operator/SIC polynomial outputs retain independent x and use
K(x) as coefficient field. The construction below applies to this
effective finitely generated field too. It does not replace a polynomial
identity by evaluation at one output point.

## 2. Effective scalar minimization over the required coefficient fields

The classical ingredients are finite right-factor degree bounds and
scalar holonomic zero testing. The general-field extension is explained
here rather than inferred from the rational-coefficient implementations.

Let f in K[[t]] be computable and let nonzero L in K[t]<partial_t>
annihilate it. First test f=0 by the retained origin-indicial oracle.
For nonzero f, its monic least-order annihilator L_min over K(t)
is a right factor of L by Ore division.

### 2.1 A computable faithful complex embedding

A finite rational prime presentation of K gives a finite separable
extension of Q(a_1,...,a_d) by Noether normalization and a primitive
element. A faithful computable embedding K -> C can be constructed.

Enumerate the nonzero rational polynomials in the a variables. Choose
nested closed rational boxes with nonempty interiors, diameter tending
effectively to zero, such that the n-th polynomial is bounded away from
zero on the n-th box. Inside the preceding box choose a rational point
where that polynomial is nonzero and a smaller box whose nonvanishing
is certified by rational polynomial bounds. Nonzero polynomials cannot
vanish on an open box. The computable limit tuple is algebraically
independent; each earlier nonvanishing condition persists on all later
boxes. This embeds Q(a) faithfully into computable reals.

The primitive minimal polynomial remains separable at this tuple.
Select a computable complex root by certified root approximation.
The nonzero field norms remain nonzero at the algebraically independent
tuple, so the resulting map on K is injective. For d=0 this is ordinary
number-field embedding. This construction need not match a user's
particular complex embedding: the requested operator identities and
degrees are invariant under the abstract field isomorphism.

Only finitely many coefficient evaluations and strict upper bounds are
needed below. A nonzero computable denominator can be bounded away from
zero by increasing precision; no equality decision on arbitrary
computable complex numbers is being assumed. Exact field equality
continues to use the given rational prime presentation.

### 2.2 A finite right-factor coefficient bound

The retained author version of
[Bostan--Rivoal--Salvy, Explicit degree bounds for right factors](https://rivoal.perso.math.cnrs.fr/articles/factorization.pdf),
dated 12 June 2019, theorem 1, gives over C, for a monic order-r right
factor of an order-m, polynomial-degree-q operator,

    deg_rat(right factor) <=
      r^2(q+1)E + r(q+1)N + rq
        + r^2(r-1)(3q+1)/2.                                (2)

E bounds the moduli of generalized local exponents and N bounds the
largest Newton slope, including infinity. This is the explicitly
retained preprint formula; later sharper variants are not substituted.
One may take N=m+q. Formal local algorithms compute the finitely many
generalized exponents by finite algebraic operations. Under the
computable embedding above, a strict integer upper bound E is obtainable
by root bounds and certified approximation. Apparent singularities may
be included for a larger safe bound.

Thus (2) supplies a computable integer B_r for every possible order.
It bounds individual numerator/denominator degrees of the monic
rational coefficients. Clearing their at most r denominators gives
a polynomial-coefficient right factor of degree at most (r+1)B_r.
The extra factor is deliberate; a rational coefficient bound alone
is not silently treated as a common-denominator bound.

The complex theorem bounds all right factors after the chosen faithful
embedding, hence bounds every factor over the original K. No theorem
about arithmetic heights of transcendental parameters is needed.
The general embedding and formal-exponent computations are not executed
by our small checker.

### 2.3 Exact finite linear algebra finds the least order

For each order r from 1 through ord(L)-1, search the finite vector space

    A=sum_(j<=r,k<=C_r) a_jk t^k partial_t^j,
    C_r=(r+1)B_r,        a_jk in K.                          (3)

The conditions A(f)=0 are linear in the a_jk and can be decided exactly.
For example compute a common scalar annihilator of the finitely many
series t^k f^(j); an origin-indicial bound then makes one finite jet
sufficient for every linear combination. Alternatively extend the jet,
compute its nullspace, and test every nullspace basis operator by the
existing scalar zero oracle. Stop when all basis operators annihilate
f. The descending finite-dimensional nullspaces eventually equal the
true relation space; testing the whole basis, not mere stabilization of
two successive jets, makes this stopping condition certified.

At the first r admitting a nonzero order-r operator, normalize it to
monic form. It is L_min. Negative lower-order tests are conclusive because
(2) bounds every possible minimal right factor. If no smaller order
occurs, normalize L. Check exact right divisibility as an output
certificate. Scalar origin tests permit singular or irregular origins
and delayed first nonzero terms.

[Bostan--Rivoal--Salvy, Mathematics of Computation93(2024),1427--1472](https://mathexp.eu/bostan/publications/BoRiSa24.pdf),
algorithm 1 and theorem 2.2, implement the corresponding minimization
principle over Q and discuss number fields. Our K-extension uses (2)
and the faithful computable embedding explicitly; it is not a claim
that their implementation supports arbitrary function fields.
Their E-function special-value applications are not needed here.

## 3. A cyclic vector with an exact finite construction

Localize H to the finite-rank differential module

    V=K(t) tensor_(K[t]) H.

All polynomial t-torsion of H is killed by lambda, because K[[t]]
is torsion-free. Compute a basis of V and its rational connection matrix
A(t), with nabla v=v'+A v for coefficient columns.

Here is an elementary effective cyclic-vector construction. Let r=dim V
and choose an ordinary rational point a of A. Such a point can be found
by exact tests against the finitely many denominators. Set w_0=e_1.
Inductively, after w_0,...,w_(j-1) have been chosen, let

    v_<j=sum_(i<j) w_i(t-a)^i/i!,
    w_j=e_(j+1)-(nabla^j v_<j)(a),      1<=j<r.

Then

    v=sum_(j<r) w_j(t-a)^j/j!

satisfies (nabla^j v)(a)=e_(j+1) for 0<=j<r. The determinant of its
first r connection iterates equals one at a, so v is cyclic over K(t).
The triangular jet argument proves this directly: the new degree-j
term affects nabla^j at a by exactly w_j and affects lower iterates by
zero.

Lift v to H after multiplying by a common nonzero denominator if needed.
Multiplication by a rational scalar preserves cyclicity by a triangular
change of the derivative columns. Call the lifted element h and write
v for its image in V. Its period

    f=lambda(h)

is computable and formal at zero. The cyclic relation for V gives a
scalar annihilator L of order r. If r=0 the whole map is zero.
If f=0, cyclicity makes the map zero on V and hence on H.
Otherwise compute L_min by section 2. Then

    actual generic rank = ord(L_min).                       (4)

Indeed the original generator images are rational combinations of
f,f',...,f^(r-1); their span has dimension exactly ord(L_min).

## 4. Contract the generic kernel back to H

Let R=K(t)<partial_t>. The cyclic identification is V=R/RL.
Its generic period kernel is

    (R L_min)/(R L).                                        (5)

For each original H-generator h_i, compute A_i in R with h_i=A_i v
in V. Right-divide A_i by L_min; its remainder of order less than
d=ord(L_min) gives the coordinates of its actual period in the basis

    f,f',...,f^(d-1).

Minimality says these d functions are independent over K(t).
Thus the companion connection W=R/RL_min injects into K((t)) by
sending its cyclic generator to f.

Choose one polynomial p(t) containing every denominator of the companion
matrix and of these finitely many generator coordinates. The meromorphic
connection

    W_p=K[t,p^(-1)]^d

has a computable holonomic polynomial-Weyl presentation by the retained
[effective localization framework](CONSTANT_COEFFICIENT_OPERATOR_DECISION.md).
Present the cleared connection relations and actually localize at p;
clearing denominators without localization can leave spurious
boundary modules.

The coordinate assignment defines an explicit D-linear map

    psi:H -> W_p.

Its finite relations can be checked by exact rational-vector reduction.
Standard Weyl-module syzygies/elimination compute its kernel from the
two finite presentations. Because W_p injects into K((t)) via the
independent derivatives of f,

    ker(psi)=ker(lambda).                                   (6)

Thus the finite quotient H/ker(psi), with its specified generator map,
is the full actual image over K[t]<partial_t>. This automatically
retains the original t-torsion in the kernel, and retains any finite
point quotients or unsaturated lattice information in the image.

**Do not replace the image by W_p.** The latter is only an ambient
meromorphic connection and can be larger. For example a safe image
K[t] with cyclic period 1+t may use p=1+t in this construction, but
does not acquire a pole at -1. Conversely the image K[t,(1-t)^(-1)]
has the same generic connection as K[t] and its pole quotient must
not disappear. Equation (6) contracts back to H to distinguish them.

If the specified pure seed maps to 1, its cyclic submodule in the
actual image is exactly K[t]. Taking the finite quotient (1) and
testing whether it is zero gives a canonical module certificate of
safety. This recovers the earlier finite-multiplier decision through
the actual kernel. In the operator case use the full localized Laurent
module and the earlier original-amplitude transfer; no rationality
transfer through Laurent shifts is assumed.

## 5. Deciding whether an actual image is algebraic

A fresh source is
[Bostan--Salvy--Singer,2504.16697v2,17September2026](https://arxiv.org/html/2504.16697v2).
Its complete algorithm 1/theorem 7 decides algebraicity of a specified
Q-coefficient D-finite series from its operator and initial data.
Minimization is essential; its faster sufficient transcendence test
can return FAIL and is not a complete algebraicity test.

For a number field K, write f=sum beta_j f_j in a rational basis of K,
with f_j in Q[[t]]. A common left multiple of the finitely many conjugate
operators gives a rational-coefficient annihilator, which kills every
f_j by basis independence. Moreover f is algebraic over K(t) iff all
f_j are algebraic over Q(t). The forward implication follows by taking
conjugates and trace-dual basis projections; the reverse follows by
finite sums of algebraic functions. Thus the complete rational test
extends effectively to number-field series.

For a number-field scalar period map, test its finitely many full
polynomial-Weyl generator images. Every image is algebraic iff all
generators are: algebraic functions are closed under differentiation,
rational scalar multiplication and finite sums. If they pass, searching
for algebraic relations, with the retained exact zero oracle, terminates
and supplies a finite algebraic function-field cover containing them.
After inverting finitely many branch/pole factors the actual connection
embeds in that finite-cover connection. It is regular singular with
finite geometric monodromy. This does not make it rank one or safe.
If a generator fails, an actual transcendental image series is detected.

No extension of the complete algebraicity algorithm to arbitrary
transcendental coefficient fields is imported here.

### 5.1 Bounded polynomial outputs cover balanced operator frontiers

The number-field test also applies when each full generator period
F_i(t,x) has a known finite polynomial x-degree bound B_i.
After clearing the finitely many rational x-coefficients in the amplitude,
evaluate on the finite rational grid {0,...,B_i} in each output coordinate.
Compute fresh scalar period annihilators at those points rather than
blindly specializing a possibly exceptional generic certificate.
Tensor Vandermonde interpolation recovers every x-coefficient series.

A polynomial in independent x with coefficients in an extension of K(t)
is algebraic over K(t,x) iff all its coefficient elements are algebraic
over K(t). One proof extends every K(t)-derivation of the coefficient
field by delta(x)=0: algebraicity forces delta(F)=0, hence delta kills
each coefficient; in characteristic zero this characterizes algebraic
elements. Therefore the finite grid tests are equivalent to algebraicity
of the entire generator, not just evidence at sampled output points.

Two all-degree cases have such bounds.

* For full SIC with every monomial of f having equal total xi and U
  degrees, Ecal(g f^m) has fixed output degree deg_U(g)-deg_xi(g)
  for each monomial g. Reverse Gaussian translation preserves this
  bound after the outside x factors are included. Balanced
  bidegree-(d,d) families in every pair dimension are covered.

* For GVC let r_0 be the smallest differential order of a nonzero
  monomial of lambda and suppose deg(P)<=r_0. For a full Laurent atom
  v^ell Z^a W^b, its m-th coefficient uses Lambda^(m-ell) on P^m
  with the translated monomial markers. Its output degree is at most

      |b|-|a|+ell*r_0+m*(deg(P)-r_0)
        <= |b|-|a|+ell*r_0.                                (7)

  Negative derivative indices give zero. A time shift changes m to
  a nonnegative shifted index without invalidating the same bound.
  Polynomial x-coefficients of an amplitude add their finite degree.
  Thus every computed full generator has an explicit finite bound.
  In particular homogeneous degree-r P with homogeneous order-r
  Lambda, including Delta^k with degree-2k P and the retained order-
  three Dvorsky setting, are covered without changing their operators.

The finite tests decide image algebraicity in these bounded-output
settings over number-field input coefficients. They do not compute
a dense SIC33 answer, a new least Delta power, or a general answer for
unbounded-output transcendental-parameter algebraicity.

## 6. Controls and remaining mathematical work

The checker constructs cyclic vectors with determinant one at the
chosen ordinary point, and minimizes explicit series inside supplied
finite coefficient bounds. It includes a Kummer polynomial whose
minimal first-order factor has coefficient degree larger than the
input equation's polynomial degree. Equal initial constants for 1 and
exp(t) do not identify the same solution. A nonminimal equation with
logarithmic solutions must not be used to declare the selected constant
period transcendental.

Separate controls retain a killed point module and distinguish the safe
polynomial lattice from its larger meromorphic ambient connection.
Algebraic 1+sqrt(1-t), logarithmic and exponential models show that
higher rank, reducibility, algebraicity and irregularity are different
questions. These are module controls, not new Gaussian counterexamples.

The full actual kernel, image presentation and exact generic rank are
now theoretically computable for every exactly presented input of the
three retained moment formulations. Image algebraicity has the explicit
number-field/bounded-output scope above. No large original module is
computed here. Geometric classification of nonrational image families,
normalized constant-b higher-rank trivial quotients, exceptional
coefficient behavior, global safe-locus stabilization, valid operator
compression, concrete lower-power frontiers and independent review
remain open.

[Source and assurance audit](actual_period_kernel_source_audit.json).

