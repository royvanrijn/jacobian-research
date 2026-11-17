# Critical normalization and arithmetic rigidity of rank-one Gaussian images

17 November 2026. Internal all-dimension/all-degree written proofs with
finite exact controls. Independent replay, formal verification and external
review remain false; priority unassessed. Older proof packets are unchanged.

The new bridge is a common arithmetic model for ALL original marked periods
after removing the selected critical action. It implies that a Gaussian
unit period with rational logarithmic derivative is a rational exponential
times a Kummer algebraic function. Consequently an actual rank-one Gaussian
image cannot be both globally regular and transcendental. This is a field
theorem, not an assertion about higher rank or general coefficient rings.

## 1. Uniform denominators without convergence

First use standard Gaussian covariance I. If P,Q have integer coefficients
and every monomial of P has spatial degree at least two, then

    Phi_Q(2t) = E[Q exp(2tP)] belongs to Z[[t]].         (1)

The same factor 2 works for EVERY original polynomial amplitude Q, in any
dimension and phase degree. Here and below coefficient integrality does
NOT imply complex convergence: P=X^4 is a divergent control.

Here is a termwise proof, also valid with indeterminate coefficient
variables. Choose one amplitude monomial X^gamma, phase monomials
X^alpha_a, and their occupation numbers k_a. Put m=sum k_a and
beta=gamma+sum k_a alpha_a. The contribution before integer coefficients
is E[X^beta]/product k_a!. Odd beta coordinates give zero. Otherwise, for
an odd prime p and q=p^j, Legendre's formula gives

    v_p(E[X^beta])
      = sum_j sum_i (floor(beta_i/q)-floor(beta_i/(2q)))
      = sum_j sum_i ceil(floor(beta_i/q)/2).

Since beta_i >= q sum_a alpha_ai floor(k_a/q), the sum over i for each
j is at least (1/2)sum_a |alpha_a| floor(k_a/q), hence at least
sum_a floor(k_a/q). This cancels v_p(product k_a!) at EVERY odd p.
At 2 the nonzero Gaussian moments are odd integers, and
sum v_2(k_a!) <= m; the factor 2^m clears the remaining denominator.
This proves (1) coefficientwise, without an analytic argument.

For P in D^-1 Z[X] and Q in D_Q^-1 Z[X] the conclusion is
D_Q Phi_Q(2D t) integral under the same spatial-degree premise.
One may replace Z by a polynomial coefficient ring throughout.

## 2. Exact critical normalization for every phase and every marker

For arbitrary P over a Q-algebra, let the unique branch g(0)=0 solve

    g=t grad P(g),
    B=I-t Hess P(g), L=B^(-1/2),
    psi=t P(g)-g^T g/2.                               (2)

Both the implicit recursion and the normalized binomial matrix root are
t-adically well defined over EVERY Q-algebra. Over fields g and psi are
algebraic over K(t): the implicit branch is an etale, relative-dimension
zero point of the polynomial critical equations. Matrix L is algebraic
too, for instance from its polynomial root equation and its invertible
linearization 2I at the origin. No infinite Morse coordinate change is
being assumed.

Write

    R_g(Y)=P(g+Y)-P(g)-grad P(g).Y-Y^T Hess P(g)Y/2.

This is a FINITE polynomial in Y of spatial order at least three. Formal
Gaussian translation followed by quadratic normalization gives

    Phi_Q = exp(psi) det(L)
       * E[Q(g+L Y) exp(t R_g(L Y))].                  (3)

For clarity, translation is the formal identity
E[F(X)]=exp(-g^T g/2) E[F(Y+g) exp(-g.Y)]; the critical
equation cancels the linear term in the exponent. The remaining quadratic
part is absorbed because L is symmetric and L^T B L=I. These operations
are coefficientwise finite: g=O(t), B-I=O(t), and the remaining interaction
has a factor t. Thus (3) retains EVERY original Q and remains valid over
nonreduced and non-Noetherian Q-algebras. Its proof is algebraic, even
when the original Gaussian series diverges.

If P,Q are integral, the sharper useful statement is

    H_Q=exp(-psi) Phi_Q,  H_Q(4t) belongs to Z[[t]].    (4)

Indeed g has integral coefficients by its recursive equation.
binom(-1/2,j)=(-1)^j binom(2j,j)/4^j implies L(4t) is
integral, as are det L(4t) and Q(g(4t)+L(4t)Y). The Taylor residual
has integral coefficients before substituting L: its degree>=3
coefficients are integer binomial expansions. Every monomial coefficient
of 4t R_(g(4t))(L(4t)Y) is divisible by 4t. At time order n an occupation
of m<=n such vertices has denominator dividing 2^m by Section 1, and
the factor 4^m clears it. This proves (4), including all amplitudes with
one phase time scale. For rational data, D_Q H_Q(4D t) is integral if
DP and D_Q Q are integral. This is an exponential denominator bound, not
a G-function or convergence claim.

Equivalently, for a phase with coefficients in any finitely generated
Z-domain A, all normalized marked coefficients belong to A[1/2]. The
explicit 4^n bound is stronger when a time denominator bound is wanted.
Amplitudes with extra coefficients use the same enlarged finite ring.

For a fixed nonsingular rational symmetric covariance C, first pass to a
finite algebraic constant extension with C=U U^T and set X=U Y. The
arithmetic conclusion needed below is preserved: choose ONE finite
arithmetic ring containing U, U^-1 and the coefficients. In original
coordinates the selected branch and action are

    g=t C grad P(g), psi=tP(g)-g^T C^-1 g/2.           (5)

They are invariant under this change. The explicit factors 2 and 4 in
(1),(4) are stated for standard covariance; no such numerical bound is
silently asserted for an arbitrary rational C.

## 3. A reduction lemma forcing the critical action to be rational

Let K be a characteristic-zero field, r in K(t), and psi in K[[t]]
algebraic over K(t). Suppose a unit h in K[[t]] satisfies

    h'=(r-psi')h,                                    (6)

and ALL coefficients of both h and the selected algebraic branch psi lie
in a single finitely generated Z-domain inside K. Then psi' belongs to
K(t). The critical recursion in Section 2 supplies this premise explicitly. The same conclusion holds after
adjoining finitely many algebraic constants. Here is the spreading
argument, including the hypotheses needed for reduction.

Reduce to a finitely generated characteristic-zero subfield containing
the finite algebraic presentation and the common arithmetic ring. Put
F=K(t,psi), or use a larger finite algebraic field containing g. This
field is regular over K: K is relatively algebraically closed in K((t)),
and characteristic zero gives separability. A smooth K-curve with this
function field therefore has geometrically integral generic fiber.
Choose a finite basis 1,e_1,...,e_(d-1) over K(t), write psi' in that
basis, and spread its finite presentation, denominators and discriminant
over a finitely generated Z-domain A. Enlarge A to include the common
coefficient ring of h and the branch coefficients of g when using (2).
For (2) those branch coefficients already lie in the ring of phase
coefficients; no infinitely generated ring is introduced.

After one finite localization of A, geometric integrality of the fibers
follows from Stacks [0559](https://stacks.math.columbia.edu/tag/0559) and
[0578](https://stacks.math.columbia.edu/tag/0578). A monic primitive
polynomial after inverting a polynomial in t preserves the finite degree,
and its nonzero discriminant preserves separability. The selected branch
then embeds each good fiber field into kappa(v)((t)). This last embedding
is essential: integrality at a single unrelated formal place would not
establish an identity on the fiber field. The finite model can include
the critical equations and the chosen etale branch. More explicitly,
choose the primitive minimal polynomial F(t,z) of psi, spread its
coefficients and preserve its z-degree and separability. Regularity over
K makes F geometrically irreducible. The integral branch gives a map
F_bar(t,psi_bar)=0 into kappa(v)((t)); irreducibility and preservation
of the degree make that map an EMBEDDING of the fiber function field.
Thus no rationality assertion is lifted from a reducible or unrelated
specialized branch. All integrality here is part of the stated premise.

For a good fiber of characteristic p, D^p is zero on kappa(v)((t)),
and also on the finite separable extension F_v/kappa(v)(t): D^p is a
derivation zero on the base, and a separable algebraic extension has no
nonzero base derivations. Since h mod v is a unit, the rank-one
Jacobson formula for (6) gives

    0=(r_bar-psi_bar')^p+D^(p-1)(r_bar-psi_bar'),
    (psi_bar')^p=r_bar^p+D^(p-1)r_bar in kappa(v)(t).   (7)

Separability forces psi_bar' itself into the base field. If a nonbase
basis coordinate of psi' were nonzero in characteristic zero, invert
one of its nonzero numerator coefficients in A. It remains nonzero
in every good fiber, contradicting (7). Such a positive-characteristic
fiber exists after every finite localization: the dominant finite-type
image in Spec Z contains an open set, by
[Chevalley, Stacks054K](https://stacks.math.columbia.edu/tag/054K).
This proves psi' rational. It does not invoke the higher-rank p-curvature
conjecture or assume quasi-unipotence of an exponential direct image.

Finally an algebraic primitive of a rational function is rational up to
a constant. One quick proof is trace: if [K(t,psi):K(t)]=d, then
(Tr psi/d)'=psi'; their difference has derivative zero and hence is in
K because it lies in K((t)). Therefore psi is rational.

The rank-one formula used in (7) is classical Jacobson; see
[Bostan--Caruso--Roques2304.05061v2, Theorem3.12 and its full proof](https://arxiv.org/html/2304.05061v2).
We use the convention h'=a h, so the vanishing expression is
a^p+D^(p-1)a. No sign convention from a displayed D+a equation is
copied without adjustment.

## 4. Arithmetic rank-one lemma over arbitrary characteristic-zero fields

Suppose h(0)=1, h'=a(t)h, a in K(t), and the coefficients of h lie in
one finitely generated Z-domain. Then

    h^N=q(t), N>=1, q in K(t), q(0)=1.                (8)

For number fields this is the classical Honda/Kronecker argument. At
every prime outside ONE finite set the ordinary series h reduces to a
unit, so a_bar^p+D^(p-1)a_bar=0. Comparing partial fractions shows that a
has no polynomial part, has only simple poles, and every residue c obeys
c_bar^p=c_bar at all primes above almost every rational p. After adjoining
the finitely many poles, the same assertion holds in a number field
containing all the residues. Kronecker then forces each residue into Q.
See Bostan--Caruso--Roques, Proposition3.1, full proof of Theorem3.6,
Theorem3.12 and the splitting formulation Theorem3.25. This proof applies
to number-field coefficients as well as Q: reductions are required at
ALL primes above each sufficiently large rational prime.

Here is the additional constant-field descent, rather than an unsupported
appeal to a theorem stated only over number fields. Spread all finite
rational data and the common arithmetic coefficient ring over a finitely
generated Q-domain B. Every closed characteristic-zero specialization
has number-field residue and inherits the ENTIRE formal equation and a
finite arithmetic coefficient ring. A nonzero polynomial part or a pole
of order>1 persists on a nonempty open of such specializations, so these
cannot occur. For a simple pole, adjoin that pole and its residue c to
the finite model and invert the finitely many quantities preserving its
simple-pole presentation. Every closed specialization now sends c to Q.
If c were transcendental over Q, Spec B -> A1_Q given by c would be
dominant, hence its constructible image would contain a nonempty open.
That open contains a nonrational algebraic closed point, whose nonempty
fiber has a closed number-field point: contradiction. If c were algebraic
but nonrational, every specialization of the contained field Q(c) would
keep it nonrational, also a contradiction. Thus ALL residues lie in Q.

Choose N clearing their denominators. The integer-residue product over
the poles is Galois invariant and has logarithmic derivative N a, so
it is q in K(t); normalize q(0)=1. The formal equation makes h^N/q
constant1, proving (8). The residue at infinity is rational as well.
This is an arithmetic FIELD statement. No residue fields, rational
functions, or Kummer classification are asserted over arbitrary rings.

The recent [Fuernsinn--Pannier2510.00892v3,10November2026](https://arxiv.org/html/2510.00892v3)
provides an effective first-order p-curvature/Kronecker development.
Its introduction and Section2 are checked for scope and residue-resultant
context only. Its Hermite--Pade prime bounds are neither needed nor
implemented here. Honda, Jacobson and Kronecker retain classical credit.

## 5. Consequences for actual Gaussian images

Let P be ANY polynomial over ANY characteristic-zero field, with a fixed
nonsingular rational covariance. Write f=Phi_1. Even without an actual
rank premise, assume ONLY r=f'/f in K(t). Apply Sections 2--4 to
h=exp(-psi)f. Its coefficients have the common arithmetic model (4),
and (6) holds. It follows that

    psi=R in K(t), R(0)=0,
    f=exp(R) a, a^N=q in K(t), a(0)=q(0)=1.           (9)

The finite constant extension used to whiten C causes no descent gap:
psi already has coefficients in K; the rational function and the
normalized identity a^N=q descend, or take a finite norm/multiple N.

In particular (9) holds for EVERY actual rank-one Gaussian image:
D-stability implies f'/f rational, and every original marked period is
a rational multiple of f. Combining with the sealed
[rank-one global lattice theorem](STEIN_ACTION_AND_RANK_ONE_LATTICES.md),

    M=(exp(R) a/p) K[t,d^-1].                         (10)

Here p,d and the original marker maps remain essential; (9) does not
replace them. The rank-one connection is globally regular on P1 iff
R=0: dlog a has only simple poles with rational residues, whereas the
derivative of any nonconstant rational R has a higher-order pole as a
one-form somewhere on P1. Hence a globally regular rank-one Gaussian
image has ONLY algebraic periods and finite cyclic generic geometric
monodromy. The abstract exp(arctan t) rank-one model from the previous
frontier is now EXCLUDED as an actual Gaussian rank-one image.

Other general consequences are useful:

* A nonrational selected critical action obstructs a rational logarithmic
  derivative of the unit period, hence obstructs actual rank one.
* If P has zero constant and linear parts, the selected branch is g=0,
  so R=0. Every rank-one such phase has algebraic marked periods.
* This does not infer actual rank one from a scalar annihilator or jets.
  It constrains any actual rank-one image AFTER rank is established.
* Finite arithmetic models alone do not force algebraicity in higher
  rank. Elliptic/Kloosterman/Airy phenomena and general regular
  transcendental realizability in rank>=2 remain separate questions.

For P=X+X^3 the critical equation is g=t(1+3g^2). Parameterizing by g,

    t=g/(1+3g^2), psi=g^2(1-g^2)/(2(1+3g^2)).

The other root is 1/(3g); these two action values are distinct. The
selected action is therefore not in K(t), so the actual image has
rank>=2. This is an exact algebraic obstruction, not a moment-prefix
inference. P=X^2/2 instead has g=psi=0 and algebraic seed (1-t)^(-1/2).
For a nonzero linear phase psi=t^2 P_1^2/2 removes the exponential;
without this normalization the common arithmetic model would fail.

## 6. What remains open and what is retained

This closes the rank-one regular-transcendental Gaussian realization
question by exclusion, without a weight-frame premise. It does NOT
prove field safety=>purity: a nonconstant polynomial unit period is
compatible with (9) with R=0. The remaining obstruction is still the
finite nonzero Stein summand of an ORIGINAL nonpure safe Gaussian,
outside the closed filtered class. Exhaustive safe PHASE structure with
original marker filtration and higher-rank regular transcendental
realization/obstruction remain open.

The [constructor](gaussian_saddle_arithmetic.py),
[checker](verify_gaussian_saddle_arithmetic.py),
[certificate](gaussian_saddle_arithmetic_certificate.json),
[receipt](gaussian_saddle_arithmetic_verification.json) and
[source audit](gaussian_saddle_arithmetic_source_audit.json) retain focused
exact identities, denominator controls, matrix root/translation controls,
the critical conjugate obstruction, Kummer residue certificates and
p-curvature controls. They test finite algebra; the universal arithmetic
spreading arguments are written mathematics, not machine formalization.

All Oaku operators, torsion witnesses and specialization thresholds remain
sealed. Sabbah/Maisonobe retain convergent rational-slope credit; arbitrary
nonconvergent formal and all-ring extensions stay separately scoped with
priority unassessed. No terminal UNKNOWN job, broad search campaign,
delegation, outreach or assurance promotion is involved. Broad goal ACTIVE.
