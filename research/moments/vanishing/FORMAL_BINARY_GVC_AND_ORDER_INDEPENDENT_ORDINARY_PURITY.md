# Formal binary GVC and order-independent ordinary purity

16 November 2026. Internal written proofs using retained internal Hall,
shifted-ray, cutoff and coefficient-ideal theorems. No novelty, independent
replay, external review or formal verification is claimed. See the
[checker](verify_formal_binary_gvc.py),
[receipt](formal_binary_gvc_verification.json) and
[source audit](formal_binary_gvc_source_audit.json).

This closes GVC for EVERY formal constant-coefficient symbol in TWO variables,
with polynomial phases and polynomial markers of arbitrary degree. The sharp
field coefficient remains

    C_D=1+floor((D+1)^2/4).

It also closes the full ordinary purity/cutoff problem in this entire binary
class over EVERY commutative Q-algebra: for each phase-degree bound D there
is a finite pure prefix, a finite symbol jet and one common linear mixed
cutoff, independent of operator order and target nilpotent thickness.
The sufficient ring constants are theoretically computable, not numerically
computed here, and the sharp field coefficient is not asserted over rings.

The general formal-symbol eventuality-to-linearity question remains OPEN in
spatial dimension at least3 and phase degree at least3. The full exceptional
Haar MZ cases also remain open; this theorem is not an exceptional transfer.

## 1. Formal lower Newton envelopes are finite

Let K have characteristic zero and lambda in K[[xi,eta]]. Its formal
constant-coefficient operator acts locally finitely on every polynomial.
For any positive weight, every lowest-weight face of its support is finite.
More strongly the ENTIRE lower envelope is controlled by a fixed finite set.

Here is the two-variable Dickson argument explicitly. In a nonempty support
S subset N^2 choose its least first coordinate a0 and the least b0 with
(a0,b0) in S. Any other coordinatewise minimal point must have second
coordinate less than b0. For each such second coordinate there is at most
one minimal point, its least first coordinate. Thus there are finitely many
minimal points M. Every point of S dominates one of M, by descending inside
its finite coordinate rectangle. Consequently, for ALL positive weights w,

    min_(alpha in S) w.alpha=min_(alpha in M) w.alpha,       (1)

and the full equality face is already in M. This remains true after any
linear coordinate change, since a formal series still has support in N^2.
No infinite support is replaced by a finite Taylor jet in this step.

## 2. The complete binary field theorem, including formal symbols

**GVCBINARYFORMAL20261116.** Let P,Q be polynomials in K[x,y], P nonzero of
degree d, and Lambda=lambda(partial) for ANY formal lambda. If

    Lambda^m(P^m)=0 for every m>=1,

then

    Lambda^m(QP^m)=0 for m>C_d degQ.                        (2)

For each d>=0 the coefficient C_d is optimal, even among finite symbols.
Pure nullity is also equivalent, for nonzero lambda and P, to existence
of a linear K-frame and a positive rational weight strictly separating
the FULL formal symbol support above the polynomial support.

We verify each extension step of the retained finite proof
[GVC2ENV/GVC2SC](../../foundations/extended-geometry/BINARY_GVC_FINITE_CERTIFICATE.md).
If lambda has nonzero constant coefficient, the top ordinary-degree part of
Lambda^m(P^m) forces P=0, so this case is excluded for nonzero P.
If lambda=0 everything is zero. Otherwise let r be its least positive order.
For r>d ordinary order gives cutoff1. For d>=r, the top ordinary-degree
part gives lambda_r(partial)^m(P_d^m)=0. These are FINITE homogeneous forms.
The retained Hall argument applies without change, including d=r: it gives
a unique K-rational direction of multiplicity e with minimum xi-exponent e
in lambda_r and maximum x-exponent t<e in P_d. Galois descent of the frame
is the same finite-form argument as in GVC2SC section2.

For s>=1 form L(s)=min(s alpha_x+alpha_y) over the FULL formal symbol
and U(s)=max(s beta_x+beta_y) over P. By(1) both envelopes have only finitely
many linear pieces. Just right of1 their active x-exponents are e>t.
At every positive rational slope, highest weighted output in the pure
identities gives A_s(partial)^m(B_s^m)=0 for the complete extremal faces.
They are finite polynomials by(1). Every output expression itself uses only
finitely many derivatives, so no infinite coefficient sum or output-weight
limit is being interchanged.

The retained shifted-ray theorem GVC2UFT therefore forbids overlap of their
x-intervals while U-L>0. Their right-left order cannot reverse before a
common threshold. The last-envelope-piece contradiction in GVC2ENV remains
valid because its ONLY support-finiteness use is finite lower and upper
envelopes, now established by(1). At the first common threshold, the
zero-output shifted-ray theorem makes the equality faces disjoint and
ordered in the same direction. The right derivative of U-L is strictly
negative. Moving slightly right gives a strict positive-weight gap. For
d=r, that negative derivative is already t-e at s=1. All faces, breakpoints
and directions used in this argument are finite and algebraic over K.

Apply the retained optimal integer-support lemma
[GVC2OC, section7](../../foundations/extended-geometry/BINARY_GVC_FINITE_CERTIFICATE.md)
to the finite minimal set M and suppP. It supplies a separating weight w
with max(w_1,w_2)/gap<=C_d. Equation(1) extends that SAME gap to every formal
symbol monomial. Weighted derivative counting proves(2), including constants
Q by the pure identities. Conversely a strict gap directly kills every pure
output. Constant phases have the elementary cutoff1=C_0. The retained sharp
finite family remains inside the formal class, proving optimality.

This proof does not assume finite formal-symbol support or an unproved
uniform finite pure prefix. It uses retained written finite Hall/shifted-ray
theorems with their original internal assurance boundaries intact.

## 3. A general bootstrap from formal field cutoffs to finite ordinary seeds

**GVCFORMALCBOOT20261116.** Fix n,D>=1 and integer C>=1. Assume over EVERY
characteristic-zero field that EVERY pure pair with an arbitrary FORMAL
symbol and phase degree<=D has the common cutoff C. Then finite integers
K,N exist, depending only on n,D,C and finite coefficient certificates,
such that over EVERY Q-algebra and every formal symbol, the first K pure
outputs characterize full purity and pure inputs satisfy m>N degQ.
The full ordinary pure ideal is a finite-jet cylinder and equals the
ordinary full N-cutoff ideal. This requires the FORMAL field hypothesis,
not merely a conclusion for finite symbols of unbounded order.

Put

    R=D(CD+1), lambda_0=lambda_<=R.                         (3)

The retained [high-tail theorem](HIGH_ORDER_OPERATOR_TAILS_AND_FINITE_SAFETY_CERTIFICATES.md)
shows over fields that

    lambda pure iff lambda_0 pure.                         (4)

For the forward implication, use the assumed cutoff C for the ORIGINAL
formal pure lambda, then REMOVE its tail, of order>=R+1>D(CD+1).
All positive pure outputs are unchanged. For the reverse implication the
finite core is a formal pure pair too, so has cutoff C; ADD the same tail.
The same exact pure-output equality applies. Mixed outputs need not agree.

Here a finite-order core purity gate cannot simply be applied to the first
tests of lambda: higher symbol coefficients can occur in those tests.
We supply the missing field-valued-point argument instead.

Let B_R be the finite rational polynomial ring in the phase coefficients
through D and symbol coefficients through R. Let J be the radical of its
full pure-output ideal. It is theoretically computable by retained
[GVCPUREIDEAL20261115](CONSTANT_COEFFICIENT_OPERATOR_DECISION.md).
Let B_infty adjoin all higher symbol coefficients, and let I_infty be the
ordinary ideal of ALL output coefficients of Lambda^m(P^m), m>=1.
Then

    sqrt(I_infty)=J B_infty.                                (5)

To prove this without an infinite-ring Nullstellensatz, take ANY prime ideal
of B_infty, quotient and pass to its fraction field. Its infinitely many
symbol coefficients define a legitimate formal symbol, and every individual
output is finite. Equation(4) says this field point kills I_infty iff it
kills J. Thus the two ideals have the same prime sets; J B_infty is radical
as a polynomial extension of a radical ideal. This proves(5) in an arbitrary
coefficient ring with countably many variables, using only prime detection.

Because J has finitely many generators, some ordinary power J^a lies in
I_infty. Its finitely many generators each have a finite ideal expression.
Choose k0 large enough to contain all those pure equations and Dk0>=R.
Thus in the finite ring B_(Dk0),

    J^a subset I_<=k0.                                     (6)

One can find such k0,a by enumerating them and exact ideal-membership tests
in finite polynomial rings. Existence is proved by(5), so this enumeration
terminates in principle. No bounded search is substituted for that proof.
The first k0 equations now force J=0 at EVERY field point, so(4) gives full
formal purity and cutoff C. This supplies the two field inputs of the
retained [GVCFORMALSEED20261115](FINITE_SEEDS_AND_FORMAL_OPERATOR_UNIFORMITY.md).

For clarity its sufficient ring constants can be written as follows. In
B_(Dk0) set A0=I_<=k0, J0=sqrt(A0), choose nu with J0^nu subset A0, and
let h bound component codimension. By(4)-(6), J0 is the extension of J, so
one may take h=binomial(n+D,D)+binomial(n+R,R). Put

    E=(CD+1)(h nu-1), K=max(k0,E), N=C+E.                    (7)

The retained direct coefficient-power/ELS argument proves the FULL ORDINARY
identities, not just radical equalities,

    I_pure,infty=I_<=K B_infty,
    I_N,infty=I_pure,infty.                                 (8)

The seed through K uses only symbol jet DK. All later coefficients are free
even on the nonreduced purity scheme. Specialization gives the assertion over
every Q-algebra, including non-Noetherian and nonreduced ones. The first k0
tests alone give the affine eventual bound m>C degQ+E; impose the first K
for all positive pure outputs as well. No infinite-ring ACC is assumed.

## 4. The entire formal binary class over every coefficient algebra

**GVCBINARYALLRING20261116.** For every D>=1 take

    C=C_D=1+floor((D+1)^2/4), R_D=D(C_D D+1), n=2.

Section2 supplies the formal field hypothesis of section3. Consequently
(7)-(8) give a finite ordinary pure prefix K_D and sufficient integer cutoff
N_D for EVERY formal binary symbol over EVERY commutative Q-algebra and
EVERY polynomial phase of degree<=D. They do not depend on operator order
or target nilpotent thickness. Full purity depends on only symbol jet D K_D.
Across fields it already depends on jet R_D by(4), but this smaller field jet
is NOT asserted for the full ordinary scheme.

No numerical general k0,nu,K_D,N_D, least ordinary prefix or full cutoff ideal
is computed here. Their theoretical finite-ring producer is described above.
The optimal field coefficient C_D is NOT transported unchanged to nonreduced
rings. General formal phases/multipliers and variable-coefficient operators
are outside the theorem. The general formal eventuality-to-linearity question
is now restricted to n>=3 and phase degree>=3; binary inputs are safe even
without separately assuming individual-marker eventuality.

## 5. Controls, sources and remaining work

The checker supplies the arithmetic thresholds C_D,R_D, finite lower-envelope
transcriptions, selected genuinely formal geometric-tail actions on polynomial
inputs, and the retained sharp-family endpoint coefficient formula. These are
exact boundary controls, not finite proofs of the Hall/shifted-ray inputs or
the universal ideal containments. No new degree/support/prime campaign,
large radical computation or ideal solver is executed.

The fresh primary-source check identifies the current paper record as
[Hexagon2610.00580v2](https://hexagonmath.org/2610.00580v2), not an arXiv ID;
its history marks v2 current. Its stated proof-review limitation is retained.
The mathematical dependence here is on the repository's exact retained
Hall-envelope and optimal-cutoff proofs, not a new assurance claim for that
paper. Long's current GitHub torus source has the same pinned Git BLOB
719ca23022821868320fc0a5b9824c051105ff7e; no new torus theorem was imported.

The next major tasks remain the full exceptional zero-value logarithmic
extension problem and general formal n>=3, phase degree>=3 eventuality versus
linearity. No statement that existing pure UNSAFE higher-dimensional seeds
contradict that conditional implication is made: they fail its individual
eventuality premise. Closure-first priority remains, with external review
and outreach deferred and all assurance flags truthful.
