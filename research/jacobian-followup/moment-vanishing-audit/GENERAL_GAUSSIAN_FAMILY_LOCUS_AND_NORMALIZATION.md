# General Gaussian families, ambient restriction and exact normalization

16 November 2026. Internal written results and a general implemented family
algorithm. Independent, formal and external assurance remain false. The
registry is authoritative. This closes the previously missing executable
**reduced family-locus** gate and the normalization-ideal decision for
finitely presented rational coefficient algebras. General numerical cutoff
and ordinary cutoff-ideal producers, least bounds and safe normal forms
remain separate obligations.

The [Python producer](gaussian_family_producer.py),
[Macaulay2 core](gaussian_family_core.m2),
[checker](verify_gaussian_family_producer.py),
[receipt](gaussian_family_producer_verification.json),
[source audit](gaussian_family_producer_source_audit.json) and retained
`gaussian_family_examples/` bind the implementation. It computes complete
families, with no fixed dimension, phase-degree or parameter-degree limit.
A cap gives UNKNOWN. Its holonomic computations are trusted CAS executions,
unlike the separate [pointwise finite IBP/jet certificates](WHOLE_GAUSSIAN_SOURCE_HOM_AND_EFFECTIVE_CUTOFFS.md).

## 1. Full mixed safety is a closed field locus, without a pure hypothesis

Fix n,d and let R=Q[c_beta:|beta|<=d], q=binomial(n+d,d), with universal
phase P_c. Write

    mu_(alpha,m)=E[X^alpha P_c^m],
    K_N=(mu_(alpha,m):m>N(|alpha|+1)),  N>=0.               (1)

These are ordinary coefficient ideals in R, including the constant markers.
Their zero sets increase with N. The retained pointwise affine-cutoff theorem
gives

    S_mixed=union_N V(K_N),                               (2)

where S_mixed means eventual vanishing for EVERY fixed original multiplier,
without an all-positive pure-null premise. The retained joint holonomic
[fiberwise-Hom theorem](FIBERWISE_POLYNOMIAL_HOM_AND_GAUSSIAN_FIELD_UNIFORMITY.md)
makes S_mixed constructible: the exact ordinary Gaussian source fiber has
polynomial Hom dimension0 or1, and safety is dimension1.

For completeness, a constructible subset of a finite-type complex variety
which is an increasing countable union of closed subsets is already one
of those closed subsets. Its closure has finitely many irreducible
components. On each component the constructible subset contains a dense
open. An irreducible variety over an uncountable field cannot have a
nonempty open covered by countably many proper closed subsets. Therefore
one member contains that entire component. Take the maximum of the
finitely many indices. This proves

    S_mixed=V(K_(N0)),  J_mixed=sqrt(K_N) for N>=N0.        (3)

Consequently there is a common sufficient affine cutoff C*deg(Q)+E on
the ENTIRE mixed-safe field locus. For example C=E=N0 works. This is
not the claim that all pure-null phases are safe. The protected Gaussian
counterexamples remain unsafe points.

Everything is defined over Q. A phase over any characteristic-zero field
has finitely many coefficients, generating a finitely generated field
which embeds in C; all original moment identities are preserved and
reflected. Thus (3) holds over all such fields. The same affine bound
holds on reduced Q-algebras whose prime residue phases belong to (3).
The ideal J_mixed represents this reduced safety functor. Unlike the
older pure-safe ideal, it does not include the requirement that the first
positive moments be zero.

## 2. All coefficient rings and the full affine cutoff topology

The earlier [nilpotent-lifting proof](NILPOTENT_LIFTING_AND_FORMAL_SAFE_COMPLETION.md)
extends to (3) with its intercept retained. Indeed for |gamma|=k<=m,

    partial_c^gamma mu_(alpha,m)
      =(m)_k E[X^(alpha+sum gamma_beta beta) P_c^(m-k)].    (4)

Its new marker has degree at most |alpha|+dk. On V(J_mixed), (4) is zero
when

    m>C|alpha|+E+(Cd+1)k.                                 (5)

Apply the same explicitly retained Zariski--Nagata radical differential-
power theorem and Ein--Lazarsfeld--Smith uniform symbolic containment in
the REGULAR universal coefficient polynomial ring over C, then contract
the resulting ordinary ideal power to R. The ambient dimension q gives

    m>C|alpha|+E+(Cd+1)(qL-1)
       implies mu_(alpha,m) in J_mixed^L.                 (6)

These imported theorems are the same ones, with the same regular-ambient
hypotheses, as in the earlier proof; no symbolic theorem is applied in a
singular parameter quotient. For a nonzero proper radical ideal, q can
be replaced by the maximum codimension of its components. The zero and
unit ideal cases have their direct interpretations.

Let A be ANY commutative Q-algebra and let I=J_mixed(P)A be the ideal of
evaluated universal safe generators. If I^L=0, (6) gives the sufficient
affine bound

    E_A[Q P^m]=0 for m>C deg(Q)+E+(Cd+1)(qL-1).            (7)

Conversely safety descends to every prime residue field, which forces
every generator of I into every prime ideal of A. I is a finitely
generated nil ideal, hence nilpotent, even when A is not Noetherian and
its whole nilradical has no bounded exponent. We obtain the full
unrestricted equivalences

    P mixed-safe over A
      iff P mixed-safe over A_red
      iff EVERY prime residue phase is mixed-safe
      iff I is nilpotent.                                (8)

Purity is not assumed in (7)-(8). If it is additionally true, the affine
intercept can be absorbed in a linear bound for positive-degree markers;
constant markers are handled by the actual pure identities. The older
pure-null theorem and its ordinary ideal boundary remain intact.

The full affine ideals (1) are cofinal with J_mixed powers. Explicitly,

    N>=ceil max(C,E+(Cd+1)(qL-1)) implies K_N subset J_mixed^L,
    for each N>=N0 some J_mixed^s(N) is contained in K_N.  (9)

The second containment is radical equality and Noetherianity of R. Thus
the full mixed-safe functor on arbitrary discrete Q-algebras is the
formal completion of Spec R along V(J_mixed): continuous maps to A are
exactly maps R->A on which J_mixed has nilpotent image. Here no pure
quotient is needed, because (1) allows an intercept. In the constant
phase control P=e, K_N=(e^(N+1)), which is indeed cofinal with (e)^L;
the earlier zero-intercept ideals I_N=(e) had a different boundary.

These statements descend by polynomial substitution to every finite
rational parameter algebra B/I0. Its reduced field-safe ideal is
sqrt(J_mixed B+I0). Cofinality after substitution uses (9) for one
direction and radical equality/Noetherianity of B/I0 for the other.
There is no ordinary stabilization claim and no bound independent of
nilpotent thickness on all nonreduced targets.

## 3. The EXACT polynomial-Hom normalization ideal is coefficient torsion

For the exact Gaussian source H_P over A, the retained completion theorem
identifies Hom_Dt(H_P,A[[t]]) with A Phi_P. Its polynomial-Hom normalization
ideal is

    J_P={a in A:a Phi_P(H_P) subset A[t]}.

For a fixed a, the displayed condition is precisely safety of the phase
over A/Ann_A(a): each coefficient is zero after multiplication by a iff
its image is zero in that quotient. By (8), this is equivalent to some
power of I being contained in Ann_A(a). Therefore

    J_P=Gamma_I(A)=union_L Ann_A(I^L)=H_I^0(A).            (10)

This holds for EVERY polynomial phase and EVERY Q-algebra, without a
pure assumption. Safety remains 1 in J_P; nonzero J_P alone is insufficient.
For a finitely presented rational A, (10) is an ordinary ideal saturation
and is computable once the family safe ideal is known. It avoids relative
polynomial-Hom base change entirely.

For P=aX over Q[a,b]/(ab), I=(a) and J_P=(b). Over
Q[a,b]/(a^2 b^2), J_P=(b^2), showing why ordinary thickness must be kept.
Over Q[e]/e^3 with P=eX, J_P=A: the actual unit period is
1+e^2 t^2/2, safe but not pure. The field/radical criterion alone must not
be reported as equality to a reduced scheme on these nonreduced targets.

## 4. Ambient local cohomology handles arbitrary exceptional strata

Let B=A^r and let M be the EXACT jointly holonomic Gaussian source on
A1_t times B, with differentiated parameter directions. For an irreducible
closed subset S of codimension c, form the ambient module

    M^[S]=H^c_(A1 times S)(M).                            (11)

On the smooth locus of S, the closed embedding i satisfies the classical
identity

    R Gamma_(A1 times S)(M)=i_+ i^dagger M,
    i^dagger=Li^*[-c].

Closed direct image is exact. Taking H^c therefore identifies (11), there,
with i_+ H^0(Li^* M), the pushforward of the ORDINARY restricted source.
Its ordinary parameter fibers are exactly H_(P_a), by associativity of
ordinary tensor products. This is Example6.62, printed64-65/PDF70-71,
in [Mustață's notes](https://websites.umich.edu/mmustata/DmodulesNotes.pdf),
with the shift on printed63/PDF69. It works for characteristic embeddings;
one must take degree c, not degree0, and must reuse the ORIGINAL M at
every node. It avoids constructing rings of differential operators on
nonlinear S. Singular points are put into the next exceptional set.

This ambient local-cohomology strategy is classical:
[Walther's stratification paper](https://www.math.purdue.edu/~walther/research/reg.pdf),
Section2/Algorithms2.1-2.2, takes arbitrary HOLONOMIC input for its
algebraic connection strata and ranks. Regularity is needed for its
subsequent identification with the holomorphic solution complex. We use
the former algebraic mechanism and the separately proved Gaussian
polynomial-Hom formula; we import neither a regularity assumption nor
an irregular analytic-solution interpretation. Walther also explicitly
warns that the graded module carries information lost by its annihilator.
The conventions and H^c selection here are fixed by Mustață's shifted
inverse image and the actual codimension controls, rather than copying
Walther's restriction indices without translation.

Compute, still in ambient Weyl algebras,

    N^[S]=pi_+ D(M^[S]).                                 (12)

Duality commutes with the proper CLOSED embedding, and pi composed with
that embedding is the closed pushforward of the time projection on S.
Thus the cohomology of (12) is the closed pushforward of the corresponding
cohomology for the ordinary S-source. We never commute duality with the
nonproper time projection.

On a dense good open of S, the polynomial-Hom dimension equals the
connection rank of H^(-1) in (12), by the retained fiberwise-Hom theorem.
Inside the ambient cotangent bundle that rank is the characteristic-cycle
multiplicity of the conormal closure of S. Retain the associated graded
MODULE, not just its annihilator: doubling a connection doubles this
multiplicity while leaving its annihilator unchanged. HTT Definition2.2.2
and Example2.2.4, printed60/PDF69, give the length and rank interpretation.

The implementation isolates this conormal component by saturating the
graded relation module away from every other characteristic component.
Its generic multiplicity is

    degree(isolated graded module)/degree(conormal prime).

The equality follows by a finite prime filtration and additivity of
top-dimensional degree. Lower-dimensional embedded pieces do not affect
it. This is checked on noncyclic modules and nonlinear conormals.

The companion [restriction experiments](https://www.math.purdue.edu/~walther/research/dexp.pdf),
Example1.1, records x^2 D_x+1: its origin restriction is zero although its
convergent analytic Ext has a nonzero class, killed formally by a divergent
factorial series. It is a warning against precisely the analytic
identification we do not use. These papers supply classical credit and
an implementation audit, not a new Gaussian safety theorem by themselves.

## 5. The finite family algorithm

Given P(c,X) with rational polynomial coefficients and an arbitrary
parameter ideal I0 in Q[c], construct the exact joint exponential ideal

    (D_Xi+(H X)_i-t P_i, D_t-P, D_cj-t P_cj),             (13)

where H is a FIXED rational symmetric nondegenerate Gaussian precision
matrix; the default is identity. Full top spatial integration gives M,
retaining every parameter derivative. Begin with the minimal primes of I0.
At a node S use (11)-(12) and the conormal multiplicity above.

Here is an explicit computable proper exceptional set. Characteristic
components are conormal closures to their base projections Z. For the
source on A1 times S, ignore Z=A1 times S. If Z does not dominate S,
remove the closure of its parameter image. Otherwise Z is generically
finite over S; remove the image of

    Z intersect V(partial_t f: f among generators of I_Z).

Over smooth S this includes the singular and vertical-tangent points of
Z, and characteristic zero makes its image proper. For each direct-image
cohomology module, remove all characteristic base projections properly
contained in S. Finally remove Sing(S). Polynomial elimination, Jacobian
minors and prime decomposition implement these operations.

On the resulting dense open the generic polynomial-Hom dimension is0 or1.
If it is1, that open is safe; closedness from Section1 makes the ENTIRE
closure S safe, including its degenerations. If it is0, all remaining safe
points lie in the computed proper exceptional set. Recurse on its minimal
primes. Every step lowers dimension, so the tree is finite. Intersect the
ideals of all safe closures to obtain the radical field-safe ideal J_f.
Overlaps are harmless: this is a labeled finite closed cover, not a claim
of disjoint strata. Repeated identical prime closures are processed once.

Finally set A=Q[c]/I0, keeping I0 ORDINARY, and compute

    Gamma_(J_f A)(A)=(0:_A (J_f A)^infinity).             (14)

This is (10) for the presented family, since nilpotence of the pullback
universal safe ideal is equivalent to nilpotence of its finite-type radical.
The phase over this coefficient algebra is safe iff (14) contains1.
The producer also finds a finite L with Gamma*I^L=0 and checks
Gamma:I=Gamma. These two finite ideal conditions certify that the
returned ideal is exactly I-torsion, rather than a guessed subideal.
The producer returns both J_f and the lift of the ordinary normalization
ideal to Q[c]. Algebraic number fields, quotient rings with nilpotents,
and finite presentations of localizations are included through I0.

## 6. Execution and limits

The input contains `dimension`, `parameter_count`, sparse `phase_terms`
whose powers list the spatial variables followed by the parameters,
optional sparse `parameter_ideal`, and optional rational `precision_matrix`.
The function `universal_family(n,d)` constructs the ENTIRE degree<=d
phase coefficient space. No support or all-positive-purity premise is used.

    .venv/bin/python gaussian_family_producer.py input.json result.json
    .venv/bin/python gaussian_family_producer.py input.json result.json --timeout 0

Macaulay2 1.22, Dmodules1.4.1.1 and BernsteinSato1.0 perform the joint
operations. Local cohomology explicitly uses the Oaku--Takayama graph
strategy, not a localization requiring an unproved saturation hypothesis.
Quotient presentations are restored before operations requiring them.
The code returns UNKNOWN on its default30second cap; errors fail closed.
Uncapped termination is the classical holonomic algorithms plus the
strict dimension descent above, not an inference from examples.

Small controls include nonlinear and singular safe closures, codimension2
specialization, a safe split cubic family, an algebraic coefficient field,
and ordinary zero-divisor thickness. A more complicated cusp degeneration
exhausted its30second cap and is retained UNKNOWN; it was not restarted.
This is not a practical complexity claim or an exhaustive computation of
large universal safe schemes. The outputs rely on the CAS functor and
commutative-algebra implementations; independent universal assurance is
not asserted.

The producer computes the reduced family locus and the exact polynomial-
Hom normalization ideal of a presented coefficient algebra. It does NOT
compute the numerical universal field constants C,E, all individual
ordinary cutoff ideals K_N, their least bounds, or general geometric safe
normal forms. Existing pure-safe uniformity, formal nonconvergent/all-ring
closures and every protected counterexample remain unchanged. Classical
algorithmic credit is explicit; Sabbah/Maisonobe retain the convergent
rational-slope credit, and extension priority remains unassessed.

## 7. The next numerical-bound bridge needs only generic SAFE components

Let S be any irreducible safe component found above, and let K=Frac(Q[c]/S).
An actual pointwise certificate for its generic phase giving a numerical
bound C deg(Q)+E gives the SAME bound on EVERY field point of S, including
its singular and denominator-zero points. Indeed each forbidden original
moment is a regular coefficient polynomial in Q[c]/S. Its image is zero
in the fraction field by the generic certificate, so it was already zero
in that domain. No denominators are inverted at the exceptional points.
Take maxima of C and E over the finitely many safe components.

Thus a numerical family producer need not stratify all possible integer
roots of centrally parametrized scalar operators. Extend the pointwise
IBP/jet backend to these finitely presented generic coefficient fields,
retain its original identities modulo the prime, and extend the actual
moment bound by this regular-coefficient argument. This is an additional
general implementation route, not an executed generic-field backend.
The Euler family t D_t-a still forbids transferring arbitrary generic
scalar solution-degree bounds to all its parameter fibers; it is the
ACTUAL Gaussian moment identities on a SAFE component that extend here.
