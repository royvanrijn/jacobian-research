# Prime-index congruences and drift reduction in every operator order

16 November 2026. Internal written proofs, without novelty, independent
replay, formal verification or external review claims. See the
[checker](verify_prime_drift.py), [receipt](prime_drift_verification.json)
and [source audit](prime_drift_source_audit.json).

This proves drift annihilation in EVERY dimension and phase degree for
every finite symbol, every algebraic formal symbol, and a larger class
of formal symbols with one finitely generated arithmetic coefficient
ring. It then proves full positive GVC in THREE variables for those
symbols with nonzero linear part, even after multiplication by an arbitrary
formal unit. No individual-marker eventuality premise is needed in that
class. The full arbitrary-formal apolar-gap question stays OPEN.

## 1. A prime-index identity with its integrality hypothesis explicit

Let p be a rational prime, R a p-torsion-free Z_(p)-algebra, P in R[x],
and lambda=sum_(alpha!=0)c_alpha Z^alpha a formal symbol with ZERO
constant term. Write T=lambda(partial) and

    b_i=c_(e_i),   D_b=sum_i b_i partial_i.

Only finitely many symbol coefficients enter T^p P^p. All of them are
assumed in R. Then T^pP^p is divisible by p! in R[x], and

    (1/p!) T^p P^p = (D_b P)^p modulo pR[x].                (1)

Here the division is in R, using p!=p*(p-1)! and the fact that the
second factor is a unit in Z_(p). No reduction of an expression containing
an uncontrolled 1/p coefficient is permitted.

To prove (1), truncate lambda far enough for this output. The commuting
multinomial expansion gives

    T^p=sum_alpha c_alpha^p partial^(p alpha)+p S,           (2)

where S is a finite constant-coefficient operator with coefficients in R
and NO constant term. Every positive derivative of P^p is in pR[x]:
modulo p, Frobenius makes every first derivative of P^p zero, hence every
higher derivative zero too. Thus p S(P^p) is in p^2 R[x].

For |alpha|>=2, partial^(p alpha) on any polynomial is divisible by
product_i (p alpha_i)!, which has p-adic valuation at least |alpha|>=2.
These diagonal terms also vanish modulo p after division by p!.
For alpha=e_i, the Hasse derivative identity gives

    (partial_i^p/p!) P^p
      =[u^p] P(x+u e_i)^p
      =(partial_i P)^p modulo p.                          (3)

Combining the surviving terms uses Frobenius once more:
sum_i b_i^p(partial_iP)^p=(D_bP)^p. This proves (1), including p=2.
It uses neither an asymptotic estimate nor an assumption p>degP.

If degP<=D, D>=1, it suffices to retain lambda through order

    r_p=(D-1)p+1.                                         (4)

Indeed a factor of order>r_p in T^p, together with the other p-1 factors
of positive order, differentiates more than Dp times. This explains
exactly how far the p-integrality hypothesis must reach for a formal
symbol at prime index p.

## 2. Drift annihilation over characteristic-zero fields

**GVCPRIMEDRIFT20261116.** Let K be any characteristic-zero field and
assume ALL coefficients of lambda and P lie in ONE finitely generated
Z-subalgebra R_0 of K. If

    T^p P^p=0 for infinitely many rational primes p,         (5)

then D_bP=0. In particular all-order pure nullity implies this identity.
If D_bP!=0, every sufficiently large PRIME index has a nonzero pure
output. This conclusion is not a claim about every sufficiently large
integer index.

Here is a specialization proof that retains the infinite symbol. Suppose
a nonzero coefficient g of D_bP survives in K. The algebra
R_0[g^-1] tensor Q is a nonzero finitely generated Q-algebra. A maximal
ideal specializes it to a number field L, by Zariski's lemma. Every
coefficient of the ENTIRE symbol specializes, because it lies in R_0;
every output specializes coefficientwise by a finite calculation. The
image of g is nonzero. The images of a finite set of generators of
R_0[g^-1] are integral away from finitely many rational primes. Choose
a p in the infinite set (5) outside these primes, and reduce at a prime
of L over p. The specialized image of g is a unit there. Equation (1)
would force the nonzero reduced polynomial (D_bP)^p to vanish over the
residue field, contradiction.

This covers every POLYNOMIAL symbol over K: finitely many coefficients,
together with those of P, generate such an R_0. It also includes arbitrary
integer-coefficient formal symbols, including lacunary transcendental ones.
Having infinitely many symbol coefficients is not itself an obstruction;
having no common finite arithmetic model is the possible obstruction.

## 3. Why all algebraic formal symbols qualify

We use one precisely scoped classical input: the Denef--Lipshitz étale
rational presentation, stated and proved as Theorem4.1 in Yurkevich,
[2010.14386v1](https://arxiv.org/html/2010.14386v1#S4.SS1).
For algebraic lambda over K, it gives

    lambda(Z)=A(Z,h(Z))/B(Z,h(Z)),                         (6)

with A,B polynomial, B(0,0)!=0, and h(0)=0 satisfying a polynomial
equation S(Z,h)=0 with a=S_h(0,0)!=0. Adjoin the finitely many coefficients
of A,B,S,P and inverses of a,B(0,0) to Z, inside K, obtaining R_0.

The formal implicit coefficient recursion puts EVERY coefficient of h
in R_0: at positive multi-degree alpha its new coefficient appears as
a h_alpha, while all remaining terms use already computed lower degrees.
There is no division by alpha! or by an increasing list of primes.
Expanding B^-1 similarly uses only its constant inverse, so (6) places
every coefficient of lambda in the same R_0. Section2 applies.

The imported input is this presentation, not any new GVC result. The
retained source's final displayed line in the proof of Theorem1.1 writes
W(x,h(x))=h(x); equation(6) there shows the last h should be f. We use
Theorem4.1 and the explicit coefficient recursion above, so that harmless
terminal typo is not a hypothesis or an identity in this argument.

More generally lambda=U mu, where U is ANY formal unit and mu has the
common arithmetic model above, has the same all-order pure tests as mu
by exact unit invariance. For nonzero P on the pure locus, mu(0)=0:
a nonzero scalar constant would make mu(partial) invertible on each
bounded-degree polynomial space. Hence

    lambda_1=U(0) mu_1,   D_(lambda_1)P=0.                  (7)

Thus arbitrary unit factors do not remove the drift conclusion. No
arithmetic model is required for U itself. The zero phase is automatic.

## 4. A full three-variable positive theorem and the general dimension drop

**GVCSMOOTHTERNARY20261116.** Let lambda=U mu be as in section3 over a
characteristic-zero field, with lambda_1!=0, and suppose (lambda,P) is
pure in ALL positive indices. In any ambient dimension n, choose a linear
frame in which D_(lambda_1)=partial_t. Then (7) forces

    P=P(y_1,...,y_(n-1)).                                 (8)

Putting nu(Z_y)=lambda(0,Z_y), the retained
[phase-power quotient theorem](PHASE_POWER_ANNIHILATORS_AND_SYMBOL_QUOTIENTS.md)
gives exact equality of all unmarked shifted outputs with (nu,P), and
equivalence of all-marker qualitative safety and linear-cutoff existence.
Every active bad marker retains its exact output. A cutoff C for the
active pair transfers to CD+1 in the full pair. Since (8) removes the
whole linear symbol, nu has order at least2, or is zero. No formal heat
conjugation or transformed multiplication rule is assumed.

For n<=3 the phase uses at most TWO essential linear coordinates, so
the full formal binary theorem proves ALL-marker safety, with sufficient
common cutoff

    D C_D+1,   C_D=1+floor((D+1)^2/4),                     (9)

for degP<=D. Thus full GVC is positive, in every phase degree, for ALL
three-variable finite symbols with nonzero linear part, ALL algebraic
formal symbols with nonzero linear part, and ALL unit multiples of
symbols with a common finite arithmetic model and nonzero linear part.
This is stronger than a conditional eventuality-to-linearity theorem.
The bound (9) is sufficient, with no sharpness claim here.

This gives the COMPLETE spatial dimension classification of the universal
GVC property within the stated NONZERO-DRIFT symbol class: it holds for
n<=3 and fails for every n>=4. For failure add one unused variable s to
the retained [GVC3HC](../../foundations/extended-geometry/THREE_VARIABLE_HOMOGENEOUS_GVC_COUNTEREXAMPLE.md)
phase P_0 of degree12, with its fixed marker Q=x^2, and set

    T=partial_s+(4 partial_x partial_y+partial_t^2)^6.       (10)

Both P_0 and Q are independent of s, so ALL its pure and marked outputs
are exactly those of the retained ternary seed. It is pure, while Q
fails at EVERY positive index. Its symbol is finite with integer
coefficients and has the nonzero linear part partial_s, so it lies even
in the polynomial subclass. Padding gives every n>4. The minimum ambient
dimension of failure in this RESTRICTED nonzero-drift class is therefore4.
This is a direct all-order dimension theorem, not a search or a change to
the unrestricted GVC dimension3 endpoint. The new failure still does not
satisfy individual eventuality and does not refute the conditional gap OP.

For arbitrary n, a concise nonzero pure phase in this arithmetic class
necessarily has ZERO linear symbol. This is a structural reduction, not
a proof of positive GVC in every dimension. Existing higher-order pure
unsafe seeds remain valid and are not changed. No minimum dimension,
degree or operator-order frontier in the UNRESTRICTED problem is changed.

## 5. Exact boundary for unrestricted rational formal symbols

For lambda in Q[[Z]], let P in Q[x] have degree<=D and assume
D_bP!=0. Outside finitely many rational primes determined by P and b,
some coefficient of D_bP remains nonzero modulo p. If T^pP^p=0 at such
a prime, (1)-(4) imply that the symbol jet through (D-1)p+1 must contain
a coefficient with NEGATIVE p-adic valuation:

    some |alpha|<=(D-1)p+1 has v_p(c_alpha)<0.              (11)

For a pure pair with nonpassive drift this is required at EVERY
sufficiently large prime. In particular infinitely many good prime jets
already force drift annihilation for rational symbols, even when the
entire coefficient set lacks a common finite arithmetic model.

This is a necessary denominator obstruction, not a realization of a
counterexample. Arbitrary formal coefficients can have new denominators
at every prime, so (1) does NOT prove drift annihilation for all formal
symbols. Nor does raw holonomicity follow from the arithmetic model:
integral lacunary series can be nonholonomic. General formal drift and
the full actual-apolar-gap implication remain OPEN outside the stated
class. A symbol presented as a unit times a model symbol must first be
stripped as in section3; (11) concerns the displayed rational symbol.

## 6. The field identity fails over nonreduced rings

**GVCNILDRIFTBOUNDARY20261116.** Over A=Q[epsilon]/(epsilon^2), take

    T=partial_x+partial_y^2,
    P=z+epsilon*(x-y^2/2).                                (12)

Then TP=0 and (partial_yP)^2=epsilon^2*y^2=0. The order-two product rule
gives, for EVERY j>=1,

    T(P^j)=j P^(j-1)TP+j(j-1)P^(j-2)(partial_yP)^2=0.      (13)

Thus all positive pure outputs vanish, and P is not nilpotent (its
reduction is z). The linear drift coefficient is a UNIT, but

    partial_xP=epsilon!=0.                                (14)

Literal drift annihilation is therefore false over general Q-algebras,
even for order2 and a nonnilpotent phase. This pair is SAFE with cutoff1,
because T annihilates every phase power and the commutator theorem applies.
It is not a counterexample to any field theorem or to the retained ordinary
all-ring binary theorem.

There is a precise ordinary-ideal consequence for finite operator order.
Fix n,D,r and lambda(0)=0 in the finite universal Q-coefficient ring.
Let J be the ideal of coefficients of D_bP and I the full ordinary pure
ideal. Section2, applied over fraction fields of all prime quotients,
proves J subset sqrt(I). Finite generation then gives some integers N,M
with

    J^N subset I_(pure,<=M).                              (15)

No numerical N,M or order-independent formal gate is computed. Example
(12) proves J need not be contained in I itself. One must retain this
nilpotent thickening rather than replace ordinary purity by a reduced
drift equation over rings. This differs from the exact all-ring passive
reduction, which starts from an ACTUAL phase splitting.

## 7. Verification and the remaining major questions

The checker verifies a few prime congruences with mixed-order symbols,
the Hasse/Frobenius coefficient identity, rational p-integral controls,
an exact example showing how a 1/p coefficient defeats an unqualified
prime reduction, a small algebraic implicit recursion, and (12)-(14)
modulo epsilon^2. They are transcription and boundary controls. The
universal congruence, arithmetic specialization, algebraic coefficient
model and dimension reduction are the written proofs above.

The general formal implication remains OPEN for residual concise phases
and symbol classes. In the new arithmetic class concise phases force
order>=2; in dimension3 nonzero drift closes the full positive problem.
No increasing prime campaign, general symbol-model recognizer or arbitrary
formal solver is implemented. Both full exceptional Haar MZ decisions
remain OPEN at their actual complement-logarithm obstruction. The earlier
radial covering obstruction and actual four-Gaussian radial bridge are
preserved. Closure-first and deferred external review/outreach persist.
