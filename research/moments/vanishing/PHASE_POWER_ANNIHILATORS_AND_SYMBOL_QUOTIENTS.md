# Phase-power annihilators and symbol quotients

16 November 2026. Internal written proofs, in every dimension and phase
degree. No novelty, independent replay, formal verification or external
review is claimed. See the [checker](verify_phase_power_quotients.py),
[receipt](phase_power_quotients_verification.json) and
[source audit](phase_power_quotients_source_audit.json).

The full formal-symbol eventuality-to-linearity implication remains OPEN.
The reduction here preserves that entire question while enlarging a closed
class: algebraicity may be tested modulo the annihilator of ALL phase powers,
rather than only up to a unit in the ambient formal symbol ring. No procedure
for computing that annihilator or recognizing such representatives is supplied.

## 1. The phase-power annihilator

Let A be any commutative Q-algebra, P in A[x_1,...,x_n], and D>=1 a bound
for its degree. Formal constant-coefficient operators act locally finitely
on polynomials. Define the ideal

    J_P = intersection_(j>=0) Ann(P^j) in A[[Z_1,...,Z_n]],
    Ann(f) = {r : r(partial)f=0}.                           (1)

Each Ann(f) is an ideal, since constant-coefficient operators commute.
The term j=0 includes P^0=1 and forces r(0)=0. The ideal J_P is a common
kernel in the ACTUAL power-apolar algebras; it is not Ann(P) alone.

**GVCPHASEQUOT20261116.** Suppose lambda-mu belongs to J_P, and put
T=lambda(partial), B=mu(partial), R=T-B. Then, for EVERY j,k>=0,

    T^k P^j = B^k P^j.                                    (2)

Indeed R P^j=0, and in the commuting binomial expansion every term
containing R kills P^j. Thus every pure prefix and every unmarked shifted
output is literally identical. Over any Q-algebra, all-marker eventual
safety for (T,P) is EQUIVALENT to that for (B,P), without a purity premise:
use (2) and the full formal shifted-power theorem
[GVCFORMALCOMM20261116](FORMAL_SHIFTED_POWERS_AND_APOLAR_NILPOTENCE_GAPS.md).
All shifts and their separate eventuality thresholds must be retained.

There is also a quantitative statement over EVERY Q-algebra. A common
linear cutoff C>=0 for B means B^m(QP^m)=0 whenever m>C degQ, m>=1.
It includes pure nullity by Q=1. Then T has cutoff

    C' = CD+1.                                            (3)

To prove it, let delta=degQ and m>(CD+1)delta. Then k=m-delta>=1 and
k>CD delta. Apply the base cutoff at index k with marker P^delta, of
degree at most D delta, to obtain B^(m-delta)P^m=0. Equation (2) gives
T^(m-delta)P^m=0. Formal ad-nilpotence on coefficients of degree delta
then gives T^m(QP^m)=0. This argument includes delta=0 and the zero marker
is automatic. Interchanging B and T proves equivalence of linear-cutoff
existence. It does NOT assert equality of least cutoffs or marked outputs.

Over a characteristic-zero field, for nonzero P on the pure locus, (2)
preserves each exact height h_j=max{k:T^kP^j!=0}, hence every deficit
b_j=j-h_j, running maximum, asymptotic slope and infimum deficit fraction.
There is no asymptotic or bounded-sample approximation in this equality.

## 2. Ordinary ideals before taking radicals

One can impose (1) without presupposing reduced coefficients. In the
universal polynomial coefficient ring for a degree<=D phase and two formal
symbols mu,r, let H be the ordinary ideal of ALL coefficients of
r(partial)P^j, j>=0, and set lambda=mu+r. Work in the arbitrary Q-algebra
obtained by quotienting by H. The preceding proof applies there.

Consequently the two pure-prefix ideals are equal modulo H, at EVERY
prefix length. Their full shifted-cutoff ideals E_N are equal modulo H
for every N, with E_N as in GVCEXCESSIDEALS20261116. For the ordinary
marked-cutoff ideals I_C, the literal containments are

    I_(CD+1)(T) subset H+I_C(B),
    I_(CD+1)(B) subset H+I_C(T).                            (4)

For (4), quotient by H+I_C(B), apply (3) over that Q-algebra, and reverse
roles for the other containment. This retains all nilpotents. It does not
claim H is finite-jet, finitely generated in the infinite coefficient ring,
or decidable by finitely many tests. No infinite-ring ACC is invoked.

## 3. Algebraic representatives enlarge the closed conditional class

**GVCPHASERESALG20261116.** Over any characteristic-zero field, suppose
there are a formal unit U and an algebraic formal symbol mu with

    lambda-U mu in J_P.                                   (5)

For a pure pair (lambda,P), individual eventuality for every polynomial
marker implies one common linear cutoff. In fact section1 transfers pure
nullity and qualitative safety to (U mu,P); exact unit invariance transfers
them to (mu,P). Apply
[GVCALGFORMALLINEAR20261116](ALGEBRAIC_FORMAL_SYMBOLS_AND_GLOBAL_LINEARITY.md),
then transfer a resulting cutoff back by unit invariance and (3).
If C is a cutoff for (mu,P), CD+1 suffices for (lambda,P).

The zero representative is allowed: lambda in J_P has cutoff1 directly,
since T P^m=0 and the commutator implication kills every marker at
m>degQ. This assertion is valid over every Q-algebra, independently of the
field-only algebraic-source theorem. For nonzero P over a field on the pure
locus, (5) also preserves the entire height sequence of mu: first use (2),
then invertibility of U(partial)^k.

This is conditional linearity, not full positive GVC for algebraic residues.
The general residual OP can therefore be represented on a concise phase
and requires a class in K[[Z]]/J_P that is NOT associated to the image of an
algebraic formal symbol. This is stronger than merely requiring the
ambient principal ideal (lambda) to have no algebraic generator. It is
still the same unresolved actual-apolar-gap problem, not a decision of it.

## 4. Arbitrary passive derivative tails, including binary phase closure

Write x=(t,y) and suppose P=P(y) is independent of ALL t variables.
Every symbol in the ideal (Z_t) belongs to J_P. For ANY formal lambda put
mu(Z_y)=lambda(0,Z_y). There is no order, weight or algebraicity condition
on lambda-mu. Then sections1-2 apply. Pure-prefix ideals are now literally
equal extensions from the active symbol ring, with no extra H necessary.
All coefficients involving at least one passive derivative are free, even
over nonreduced bases. Bad active markers retain their exact nonzero outputs.

**GVCPHASEBINARY20261116.** In any number of ambient variables, every
pure formal pair whose phase depends on at most TWO linear coordinates is
safe over a characteristic-zero field. For degP<=D, a sufficient common
field cutoff is

    D C_D+1,  C_D=1+floor((D+1)^2/4),                       (6)

by the full formal binary theorem and (3). Padding by dummy variables
covers unary and constant phases. No sharpness of (6) is asserted.

More strongly, over EVERY Q-algebra with a prescribed splitting, for which
P lies in two polynomial coordinates, the retained binary all-ring
constants K_D,N_D give the SAME ordinary pure gate K_D, using only the
ACTIVE symbol jet through order D K_D, and common mixed cutoff

    D N_D+1.                                              (7)

Indeed the full pure ideal is exactly the active binary pure ideal
extended by freely adjoining every passive-tail coefficient. In its
universal pure quotient the active cutoff N_D holds, so (3) gives (7).
Conversely all marked-cutoff ideals contain the pure ideal. Thus in that
universal ring I_(D N_D+1)(lambda)=I_pure, with all nilpotents retained.
These constants are independent of ambient dimension and passive-tail
order. General numerical K_D,N_D remain uncomputed, as in the retained
binary theorem. Over rings we require an actual splitting (or invertible
linear frame); a rank condition on a module need not supply one.

Over a field let W={v:D_vP=0}. Choosing a basis of W as passive coordinate
directions realizes P in n-dim(W) essential linear variables. Therefore a
minimal-dimensional counterexample to the remaining formal-linearity
implication may be taken concise, with essential dimension>=3. This
concerns PHASE variables: it differs from the already retained theorem
for operators differentiating in at most two active directions.

## 5. Sharp loss, a concise example, and necessary boundaries

The loss CD+1 in (3) is optimal as a general real cutoff bound. For any
integers r>D>=1 take P=y^D, B=partial_y^r and
T=partial_y^r+partial_t. The base sharp real cutoff is C=1/(r-D), from
degree counting and the marker y^((r-D)m) at index m. A term on the marker
t^a y^b can survive T^m only if, for some l<=a,

    r(m-l)<=Dm+b,  hence (r-D)m<=ra+b<=r(a+b).

Thus the full cutoff r/(r-D)=CD+1 suffices. At m=r L, marker
t^((r-D)L) has EXACTLY one surviving allocation, and the nonzero output is

    binomial(r L,(r-D)L) ((r-D)L)! (D r L)!.

Equality is attained for every L>=1, proving sharpness. Even though
T^kP^j=B^kP^j for every j,k, their least marked cutoffs differ.

The new algebraic-residue condition is strictly broader than the retained
ambient algebraic-principal-divisor condition. Put
phi(Y)=sum_(k>=2) Y^(k!), lambda=Z^r+X-phi(Y), P=z^D, r>D.
Then lambda-Z^r belongs to J_P. The ambient principal ideal (lambda) has
NO nonzero algebraic formal element: the invertible polynomial formal
coordinate change X'=X+Z^r takes it to (X'-phi(Y)), to which the retained
proof in GVCQUALNONHOL20261116 applies, including extra variables.
Algebraicity is preserved by this polynomial change and its inverse.
Nevertheless the pair is pure safe with the sharp cutoff r/(r-D): every
passive factor differentiates x,y at least once, so the preceding degree
allocation bound applies, and the same x-power marker realizes equality.
The full raw Taylor source remains nonholonomic by source rigidity; (5)
does NOT make the original raw source holonomic.

Finally J_P need not be generated by passive LINEAR directions. In four
variables take

    P=a+2tb+t^2 c,
    R=partial_a partial_c-(1/4)partial_b^2.

For every j, R(P^j)=j(j-1)(t^2-(1/4)(2t)^2)P^(j-2)=0
(j<2 separately). No nonzero constant vector annihilates P: coefficients
of b,c in D_vP first force v_t=0, then 1,t,t^2 force v_a=v_b=v_c=0.
Thus P is concise, while J_P contains this nonzero quadratic symbol.
This example illustrates the quotient, not a counterexample to safety.

Annihilating P alone is insufficient. For P=y, r=Y^2 annihilates P but
not P^2. Taking B=0,T=partial_y^2 gives T^0P^2=P^2 but T P^2=2,
while B P^2=0. Likewise exact unmarked equality must not be promoted to
equality of marked outputs: for P=y and T=partial_t,B=0, T(tP)=y.

## 6. Verification and remaining gate

The checker uses a handful of exact shifted-output controls, a small
dual-number coefficient identity, the sharp r=3,D=2 boundary, and the
concise quadratic-annihilator example. These are transcription controls,
not proofs of the universal statements. All-order claims follow from
commuting binomial expansion, the retained formal commutator theorem and
the stated retained binary/algebraic theorems. No new external theorem or
fresh exhaustive literature coverage is claimed.

The decisive open gate remains: for actual power-apolar algebras of a
concise polynomial phase and a residual symbol class as above, does
unbounded deficit imply a positive deficit fraction? No arbitrary sequence,
finite zero window, quotient dimension count or Noetherianity claim resolves
it. The full exceptional Haar MZ cases and their actual complement-logarithm
obstruction remain open. Closure-first and deferred external review persist.
