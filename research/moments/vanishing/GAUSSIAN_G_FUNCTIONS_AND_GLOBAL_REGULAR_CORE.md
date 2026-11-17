# Gaussian G-functions, critical induction and the globally regular core

17 November 2026. Internal all-dimension/all-degree application of classical
arithmetic differential equations. All independent replay, formal verification
and external review flags remain false; priority unassessed. Earlier Oaku,
formal/all-ring, Stein and critical-arithmetic packets remain unchanged.

The [critical normalization](CRITICAL_NORMALIZATION_AND_RANK_ONE_ARITHMETIC.md)
now gives a higher-rank bridge. Over a number field, EVERY conjugate-convergent
normalized original period is a G-function. If the selected critical action
is nonzero, the actual generic Gaussian image has NO nonzero globally regular
submodule, in ANY rank and over ANY characteristic-zero field. When the action
is zero, its globally regular core over a number field is exactly its
conjugate-convergent part. This does not construct a regular transcendental
Gaussian period, or prove that none exists in rank>=2.

## 1. Arithmetic plus convergence, with every hypothesis retained

Let K be a number field, P in K[X], and C a fixed nonsingular rational
symmetric covariance. Use the selected branch

    g=t C grad P(g), g(0)=0,
    psi=tP(g)-g^T C^-1 g/2,
    H_Q=exp(-psi) Phi_Q.                              (1)

The critical-normalization theorem gives algebraic psi and, after a finite constant
extension and ONE phase scale, exponentially bounded common denominators
for EVERY original polynomial amplitude Q. For standard covariance and
integral data the precise bound is H_Q(4t) integral. A rational time
combination of finitely many H_Q also has exponential denominators after
removing a finite Laurent principal part and changing the finite scale:
Taylor coefficients of any fixed rational function regular at zero have
such a bound by geometric expansion of its denominator.

Each H_Q is D-finite over K(t). Indeed Phi_Q is D-finite by Gaussian
holonomicity; on K(t,g), exp(-psi) satisfies a rank-one differential equation,
and finite restriction of scalars makes it D-finite over K(t). Products of
D-finite functions are D-finite. This argument does NOT claim that the full
actual image is an algebraic-exponential connection in general.

Suppose Phi_Q converges at the origin for EVERY embedding K into C.
The algebraic branch psi converges at every such embedding, hence so does
H_Q. Cauchy's estimate, separately at finitely many embeddings, bounds all
conjugates of its n-th coefficient by A B^n. Combined with the denominator
bound, these are precisely the STRICT G-function hypotheses. No analytic
convergence has been inferred from integrality alone. A finite Laurent
shift handles generic actual periods with a pole at zero in the same way.

Therefore

    conjugate convergence of an original Gaussian period
      => its critical-normalized period is a G-function.          (2)

For K=Q there is just one complex coefficient embedding. For a general
number field convergence at one embedding is not substituted for the stated
all-embedding hypothesis. Coefficients in an arbitrary Q-algebra are not
called G-functions.

## 2. Classical imports and the density convention

The published primary source
[Andre, Annals151(2000),705-740,Section3](https://annals.math.princeton.edu/2000/151-2/p05)
states the Chudnovsky minimal-operator theorem3.2, the Andre arithmetic-radius
comparison3.3, Katz regularity/rational exponents3.4, and differential-module
closure3.6. The retained PDF pages14-16 (printed718-720) were read, with
pages14-15 visually checked. The primary author is restating some earlier
deep inputs through cited proofs; those earlier proofs are NOT independently
reproved or claimed fully verified here.

[Lepetit, NWEJM7(2021),81-149](https://nwejm.univ-lille.fr/index.php/nwejm/article/view/23)
Definition2 and Theorem1, printed82-83, gives the strict definition and a
stronger broad-sense theorem. We need only the strict case. Its proof chain
through Sections2-3, especially the final synthesis printed110, and its
nilpotence convention printed103-107 were checked for exact scope. A separate
strict-case exposition is [Lepetit2109.10239v1](https://arxiv.org/html/2109.10239v1),
Theorem5 and Section3.5; this is a thesis, NOT the published paper with a
different title. No claim of independently redoing every Pade estimate is made.

The classical conclusion for a NONZERO G-function is that its MINIMAL
operator is globally regular singular with rational local exponents.
Its connection has quasi-unipotent local monodromy and nilpotent p-curvature
on a density-one set of primes. In the source's terminology this is
"global nilpotence": the exceptional set has Dirichlet density zero.
We explicitly retain that convention. We do NOT upgrade it to a cofinite
prime assertion, vanishing p-curvature, a full horizontal basis of original
periods modulo p, or an unconditional geometric-origin theorem.

An arbitrary annihilator may contain extra irregular factors. We apply the
theorem to the ACTUAL function's minimal cyclic module, then use finite sums
and subquotients. No assertion about an ambient annihilator's entire solution
space is substituted for the actual image. This is classical
Andre--Chudnovsky--Katz credit, not a novelty claim for G-function theory.

## 3. The convergent marked image embeds in a critical exponential induction

Let F=K(t), L=K(t,g) in K((t)), and V=F Phi(K[t,X]) be the FULL actual
generic image. It is finite dimensional and D-stable. Suppose EVERY original
marked period converges at every coefficient embedding; it suffices to have
this for a faithful finite D-generating set of ORIGINAL periods.

Put

    N=L-span{ exp(-psi) h : h in V } inside K((t)).       (3)

It is finite dimensional over L, and D-stable by
D(exp(-psi)h)=exp(-psi)(Dh-psi' h). Choose a finite L-basis from normalized
original periods. Each is a G-function by (2). Its differential L-span is a
quotient of the pullback of its actual G-function cyclic module. Therefore
N is a quotient of a finite direct sum of such pullbacks. In particular
N is regular singular, quasi-unipotent, and nilpotent in the density-one
arithmetic sense on the smooth projective critical curve.

There is a FAITHFUL F-differential embedding

    V -> Res_(L/F)( exp(psi) tensor N ), h |-> h,         (4)
    dim_F V <= [L:F] dim_L N.

The map is the original scalar-function inclusion; hence injective. We do
NOT assume that L tensor_F V -> L V is injective: algebraic relations can
make that map lose rank. This distinction preserves the actual image while
avoiding a false tensor-product rank calculation. All original amplitude
maps survive (3)-(4). No full image-kernel computation is needed.

Equation(4) is a general structural constraint on the CONVERGENT class.
It does not say that N is algebraic, Liouvillian or geometrically realized.
Nor is it an algorithm recognizing all-marker convergence from a finite
moment prefix. Regularity and actual Gaussian realization are additional hypotheses; no unrestricted finite algebraic-exponential induction theorem is asserted.

## 4. Every nonzero globally regular actual period forces psi=0

First work over a number field. Let 0!=h in V have a globally regular
minimal cyclic connection. After multiplication by a rational function,
h is a polynomial-time combination of original periods, regular at zero.
Such a multiplication preserves global regularity. Frobenius regular-singular
theory implies convergence of that formal power series at every embedding.
Thus H=exp(-psi)h is a G-function by Section1, and its cyclic connection is
globally regular by the classical theorem.

Pass to L. Let W_h=L<D>h and W_H=L<D>H be the ACTUAL scalar-function
spans. These are quotients of the pullbacks of the two regular connections,
so both are regular. They have the same positive dimension r and

    W_H = exp(-psi) tensor W_h.

At determinant level their connection forms differ by -r dpsi. Regular
rank-one connection forms have at most simple poles after a meromorphic
gauge; a meromorphic gauge changes a connection form by a logarithmic
derivative, also with at most simple poles. At any pole of the meromorphic
function psi, dpsi has a pole of order>=2. This cannot occur in the difference
of those determinant connections. Thus psi has no pole on its projective
critical curve; it is constant. Since psi(0)=0, psi=0.

This argument works in every rank and requires just ONE nonzero globally
regular actual scalar period. It is not a claim that a globally regular
quotient of an arbitrary larger source provides such an embedded period.

Here is descent to ANY characteristic-zero field, without guessing that
actual ranks survive specialization. Retain one exact polynomial-coefficient
Fuchsian operator annihilating the actual nonzero h. A minimal globally regular
operator is such an operator; a nonminimal verified Fuchsian annihilator also
suffices. Write it as sum_(i=0)^r p_i(t) D^i. After factoring p_r into
squarefree pairwise coprime factors q_j with multiplicities m_j, the exact
Fuchs conditions are

    q_j^max(m_j-r+i,0) divides p_i for every i<r,
    deg p_i <= deg p_r-r+i, i<r.                      (5)

These are finite polynomial identities and degree bounds. Spread their
coefficients, quotient witnesses and factors over a finitely generated
Q-domain B. Invert leading coefficients of p_r and the q_j, their
discriminants and pairwise resultants, and all coefficient denominators.
Then every number-field specialization is a Fuchsian operator of the same
order, annihilating the specialized ORIGINAL period. The
[Fuchsian utility](gaussian_gfunction_regular_core.py) records these finite
quotients and a sufficient specialization guard over QQ(parameter fields).
Its output alone does not verify that a supplied operator annihilates a
Gaussian period.

If psi were nonzero, invert one of its nonzero formal coefficients and
one nonzero coefficient of h in B. These are finite algebraic data in the
critical recursion and original moment formula. A closed point of Spec B
has number-field residue. It preserves both nonzero series and the Fuchsian
annihilation identity. The number-field argument gives psi_specialized=0,
contradicting the retained nonzero coefficient. Therefore psi=0 before
specialization. No actual-image rank equality in a fiber was used.

We have proved, over EVERY characteristic-zero FIELD and in EVERY
dimension/degree,

    psi != 0 => V has NO nonzero globally regular submodule.        (6)

In particular global regularity of the unit period, global regularity of
the full actual image, and safety each imply zero selected critical action.
The necessary condition psi=0 is not sufficient for safety or regularity.
The general field safety=>purity question remains open: polynomial seeds
and their finite nonzero Stein summands are compatible with psi=0.

## 5. Exact regular-core description over number fields

Let V_reg be the largest globally regular differential submodule of V.
Let V_conv consist of those actual Laurent series that converge at zero
at EVERY coefficient embedding. It is an F-vector space stable under D.
The denominator theorem applies to each of its members after clearing its
finite rational time denominators, even if convergence results from an exact
cancellation between divergent marked periods.

If psi=0, (2) says every h in V_conv is, up to a finite Laurent shift, a
G-function. Hence its actual cyclic module is globally regular, so
V_conv is contained in V_reg. Conversely every member of V_reg converges
at every embedding by Frobenius theory. Together with (6),

    V_reg = V_conv if psi=0;
    V_reg = 0      if psi!=0.                          (7)

In the first case the regular core is quasi-unipotent and has density-one
nilpotent p-curvature, by choosing finitely many convergent generators and
taking finite sums/subquotients of their G-modules. In particular,

    full V globally regular
      iff psi=0 AND EVERY original marked period is
          convergent at EVERY coefficient embedding.              (8)

Equation(8) is an all-degree structural criterion, not an implemented
convergence oracle or a claim that one convergent unit controls all markers.
Rational core and regular core are different notions. No transcendental
regular Gaussian realization is inferred from the existence of abstract
G-functions such as an elliptic hypergeometric period.

## 6. Verification and remaining questions

The [checker](verify_gaussian_gfunction_regular_core.py),
[certificate](gaussian_gfunction_regular_core_certificate.json) and
[receipt](gaussian_gfunction_regular_core_verification.json) retain exact
Fuchsian divisibility, infinity bounds, specialization guards and indicial
controls. They do not prove the classical ACK theorem, decide convergence
from finite jets or compute an actual-image kernel.

Convergence at zero is insufficient for global regularity: exp(t^2/2) is
convergent there but globally irregular. Zero critical action and arithmetic
integrality alone are insufficient as well. The theorem requires every
original marked period and its stated arithmetic hypotheses.

The higher-rank question is to construct a genuinely transcendental convergent
marked G-part at zero critical action with every original marker controlled,
or prove an additional Gaussian obstruction. The rank-one original-vacuum
problem remains separate. Independent mathematical review and formal
verification are outstanding.
