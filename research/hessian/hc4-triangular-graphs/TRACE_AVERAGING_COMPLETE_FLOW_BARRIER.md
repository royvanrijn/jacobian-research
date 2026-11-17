# Trace averaging does not preserve complete additive flows

Date: 2026-11-16. Claims HC4MYLF0 and HC4MYLF1.

FULL THEOREM TARGET: prove HC4/JC2 or construct an actual counterexample.
CURRENT LAST UNSOLVED IMPLICATION: an ACTUAL polynomial plane source forces
Tr(Lambda_L*Lambda_L) subset B, equivalently regular inverse Gram entries;
the alternative full marked-energy integrality implication is also OPEN.
Neither is proved here. External review/outreach remains DEFERRED.

LE0/1 reduces actual complete-flow inversion to ONE component derivative
and divergence descending to the target field. KX1 also constructs source
additive actions on the full finite normalization when an integral source
coordinate is supplied. It is tempting to combine these results by tracing
the source velocities and treating their averages as a complete target flow.
This note proves precisely what tracing does preserve and breaks the extra
completeness inference in EVERY ODD DEGREE, by a written universal family.

## 1. LF0: trace gives a derivation and a divergence formula

Let B=C[P,Q] subset S be a finite normal domain extension with
F=Frac(B), L=Frac(S)=C(x,y), n=[L:F]. Suppose P,Q are algebraically
independent in this source rational plane, and put
j=Jac(P,Q), a nonzero rational function. A polynomial derivation delta
of C[x,y] which preserves S defines

    epsilon(b)=(1/n)*Tr_(L/F)(delta(b)), b in B.          (1)

The normalized trace is B-linear. Thus epsilon is a DERIVATION of B:
apply delta to b*c, then move each base factor outside the trace.
Trace sends S to B, since those traces are integral over normal B and
belong to F. In the source-plane models below S=C[x,y] itself.

The exact rational divergence identity is

    div_Y(epsilon)
      =(1/n)*Tr(div_(x,y)(delta)+delta(j)/j).            (2)

Indeed the Lie derivative of dP wedge dQ=j*dx wedge dy gives
partial_P(delta P)+partial_Q(delta Q)=div(delta)+delta(j)/j.
Trace commutes with the uniquely extended target derivations, by their
uniqueness on a finite separable extension and the conjugate-sum formula.
Taking trace proves (2). If delta preserves the pulled-back target volume
j*dx wedge dy, its averaged target divergence is ZERO. A polynomial
divergence-zero plane field is Hamiltonian: integrate the closed polynomial
one-form epsilon(P)*dQ-epsilon(Q)*dP explicitly. These conclusions require
no action descent and assert NO completeness or local finiteness.

For an ACTUAL plane Keller map j is constant. Thus every normalization-
preserving source LND produces a polynomial divergence-zero target field
by (1). This is a correct general construction. But (1) is a LINEAR
projection of velocities, not an equivariant quotient of their flows.
The following family proves that even preservation of the finite order,
pulled-back volume and the entire ramification divisor cannot by themselves
justify completeness. It also shows why this construction does not supply
LE1's actual field membership of the UNAVERAGED component derivative.

## 2. LF1: a counterfamily in every odd degree

For ANY integer k>=1 put n=2k+1. Use independent source coordinates x,y,
target coordinates P,Q, and define

    P=x^2+y,
    q_k(P,x)=sum_(j=0)^k (-1)^j*binomial(k,j)
                         *P^(k-j)*x^(2j+1)/(2j+1),
    Q=q_k(x^2+y,x),
    pi_k=(P,Q):A2_(x,y)->A2_(P,Q).                      (3)

Equivalently q_k is the polynomial antiderivative from0 to x of
(P-t^2)^k with P held constant. In coordinates (P,x),

    partial_x q_k=(P-x^2)^k,
    Jac_(x,y)(P,Q)=-y^k.                               (4)

The coefficient of x^n in q_k(P,x) is c_k=(-1)^k/n, a nonzero CONSTANT.
Hence c_k^-1*(q_k(P,T)-Q) is monic of degree n in T over B=C[P,Q].
It gives the exact free presentation

    S=C[x,y]=B[T]/(c_k^-1*(q_k(P,T)-Q)), y=P-T^2.        (5)

Indeed the relation is linear in Q before taking the quotient, leaving
the polynomial ring C[P,T]. It is a domain, and the displayed monic
basis proves field degree EXACTLY n. The map pi_k is FINITE and S is its
full normal finite normalization. Its exact ramification divisor is
k*E, where E={y=0}. Its entire etale locus is

    X={y!=0}=A1_x times Gm_y,                           (6)

which is NOT THE WHOLE AFFINE PLANE. In particular (3) is NOT Keller:
the Jacobian -y^k is nonconstant. This remains explicit throughout.

Take the source field

    delta_k=y^(k+1)*partial_x.                          (7)

It is locally nilpotent on the FULL S with polynomial flow
(x,y) maps to (x+t*y^(k+1),y). It fixes the entire E POINTWISE and
preserves the ramification divisor. On X it is a free complete additive
action. Its source divergence is zero and delta_k(-y^k)=0, so it also
preserves the EXACT pulled-back target volume

    pi_k^*(dP wedge dQ)=-y^k*dx wedge dy.

That volume has the global regular primitive P*dQ. Contraction by delta_k
has the global polynomial Hamiltonian

    i_delta pi_k^*(dP wedge dQ)
      =-y^(2k+1)*dy=d[-y^(2k+2)/(2k+2)].                (8)

Thus lack of exactness, volume preservation, normalization stability or
completeness of the SOURCE is not the reason tracing fails.

## 3. Universal averaged field and its finite-time pole

Let epsilon_k be the normalized trace field (1) for delta_k. There is a
rational number gamma_k, whose explicit value is NOT needed, such that

    epsilon_k(P)=-2*kappa_k*P*Q,
    epsilon_k(Q)=kappa_k*Q^2+gamma_k*P^n,
    kappa_k=n^2/(n-2)!=0.                              (9)

Here is a universal proof without expanding the trace at arbitrary k.
Give x weight1 and y weight2, so P has weight2 and Q weightn. The finite
free presentation is homogeneous, so trace preserves weights. The source
field delta_k has weight n. Thus the only possible target monomials in
epsilon_k(P), of weight n+2, are P*Q; the only possibilities in
epsilon_k(Q), of weight2n, are Q^2 and P^n. This uses n odd and n>=3.
No polynomial coefficient of an unspecified degree is being guessed.

We determine the Q^2 coefficient by setting P=0. There

    Q=c_k*x^n, y=-x^2,
    (partial_x Q)|_(fixed y,P=0)
       =(-1)^(k+1)*x^(n-1)/(n-2).                     (10)

For (10), the partial derivative with P held fixed is (-1)^k*x^(n-1).
The chain-rule correction is 2x times partial_P q_k. At P=0 only the
j=k-1 term of that partial contributes, namely
k*(-1)^(k-1)*x^(n-2)/(n-2). Their sum is (10).
Therefore on P=0,

    delta_k(Q)=x^(2n)/(n-2)=n^2*Q^2/(n-2).              (11)

This equality holds on every generic point of that fibre, so its
normalized trace is the same number. Base change of finite FREE trace
justifies this specialization. Equation(2) makes the target divergence
zero, so the coefficient of P*Q must be -2*kappa_k. This proves (9).

The target line P=0 is invariant and its exact equation is

    Q'=kappa_k*Q^2,
    Q(t)=Q0/(1-kappa_k*Q0*t).                           (12)

For Q0!=0 it has a finite-time pole. Hence epsilon_k is NOT
HOLOMORPHICALLY COMPLETE and NOT LOCALLY FINITE. This holds for EVERY
k>=1, with no parameter or coefficient search. All epsilon_k are
nevertheless polynomial Hamiltonian fields, with Hamiltonian

    -kappa_k*P*Q^2-gamma_k*P^(n+1)/(n+1).

We intentionally leave gamma_k unevaluated in the general theorem.
It is irrelevant to the invariant line and to the obstruction.

## 4. Field and image constraints do not repair the averaging step

This family also has a target-prime unramified sheet and the same
no-Galois-intermediate field property used in LE0. Let

    beta_k=integral_0^1(1-t^2)^k dt
          =2^(2k)*(k!)^2/(2k+1)! !=0.

The image of E is the irreducible cusp

    Q^2=beta_k^2*P^n, parametrized by (P,Q)=(x^2,beta_k*x^n).

For P!=0 the two critical points of q_k(P,x) are +-sqrt(P), each
of derivative multiplicity k. Their two critical values are nonzero
and opposite. Over a generic point of this irreducible branch curve,
exactly ONE critical root is present, of multiplicity k+1, and there are
k OTHER SIMPLE ROOTS. Thus divisorial inertia has cycle type(k+1)*1^k
and fixes an unramified sheet. All other target prime divisors are
unramified. The same fixed-sheet/purity/simple-connectedness argument
written in retained KMFI1 makes the sheet stabilizer normally generate
the Galois group. Therefore this finite extension has NO NONTRIVIAL
INTERMEDIATE FIELD GALOIS OVER F. No particular full monodromy group is
needed or asserted in this proof.

Also P is an actual source polynomial coordinate in S=C[P,x], so C(P)
is relatively algebraically closed in L. Both constants regularity and
Galois-intermediate exclusion from LE0 are therefore present. They do
not make the TRACE of a source additive field complete.

The etale source X in (6) maps onto the target minus the SINGLE ORIGIN.
For P!=0 a fibre can lose at most its one critical root, leaving k simple
roots. For P=0,Q!=0 all n roots are nonzero and etale. Over(0,0) the
only source point is x=y=0, which was removed. This even retains the
finite target-image complement property. The critical missing hypothesis
for a plane Keller counterexample is still the WHOLE etale source A2;
X is Gm times A1, has nonconstant units and Euler characteristic0.
No degree-three, quartic, HC4 or JC2 counterexample is claimed.

## 5. Exact cubic bracket, flow and boundary failures

At k=1, this is the standard universal cubic normalization, in polynomial
coordinates already present in KU0:

    P=x^2+y, Q=P*x-x^3/3=xy+(2/3)x^3,
    x^3-3*P*x+3*Q=0.

For delta=y^2*partial_x, exact trace gives

    epsilon=-18*P*Q*partial_P+(11*P^3+9*Q^2)*partial_Q.  (13)

The source field eta=partial_x COMMUTES with delta. Both are complete
source LNDs preserving the finite normalization, ramification support
and pulled-back target volume. Their averaged fields are

    T(eta)=3*P*partial_Q,
    [T(eta),T(delta)]
        =-54*P^2*partial_P+108*P*Q*partial_Q !=0.        (14)

So trace averaging is NOT A LIE ALGEBRA HOMOMORPHISM, even on this
commuting pair of source LNDs. A derived Lie bracket or commuting-flow
argument cannot silently be pushed through (1).

Nor do averaged time operators T_t(b)=Tr(exp(t*delta)b)/3 form the
flow generated by epsilon. Already on P,

    Tr(delta^2 P)/3-epsilon^2(P)
       =220*P^4-90*P*Q^2 !=0.                          (15)

The second time derivatives at zero differ. Finite polynomial source
flows and regular polynomial averaged velocities do not give a target
group law.

Finally take the branch equation B0=9*Q^2-4*P^3. Although delta fixes
the ENTIRE ramification curve pointwise,

    epsilon(B0) mod(B0)=486*P^3*Q !=0.                  (16)

Thus the averaged target field is NOT LOGARITHMIC even for that whole
branch curve. Contributions of the other retained unramified sheets
matter. Preserving or fixing a source boundary does not prove that the
averaged target velocity preserves its image.

## 6. Consequence for the HC4 route and assurance

The previous goal turn is PROGRESS: LE0/1 changed the mathematical authority
and supplied one-component descent and full pencil-group constraints.
This turn breaks a proposed way of making the required descending field
automatic. For an actual source-normalization LND, regular trace does give
a polynomial Hamiltonian target field; it does NOT by the finite-cover,
volume and boundary hypotheses above supply completeness, a group law,
Lie compatibility or whole nonproperness preservation.

One may still seek a theorem which uses the actual WHOLE polynomial etale
source A2 to supply those conclusions. Nothing here disproves such a
Keller-specific theorem. The actual LE1 unaveraged component landing,
uniform trace-dual regularity and full marked-energy integrality remain
OPEN. The family is a control of an invalid general inference, not a
replacement for the full HC4/JC2 objective.

The basic field-trace normalization, conjugate formula and multiplication-
matrix definition were rechecked in the primary Stacks Project Section9.20,
tag0BIE, including Lemmas9.20.2-3 and their proofs:

    https://stacks.math.columbia.edu/tag/0BIE

They are retained standard foundations, not a new major theorem import.
The polynomial integral/cusp models and basic trace methods are not claimed
new constructions or priority results. The universal averaging barrier is
proved above. The checker tests k=1,2,3 as exact transcription controls
and the parameter-independent formulas separately; it does NOT infer the
all-k theorem from those cases. No higher-degree coefficient campaign or
external completeness classification is used. Old registered mathematical
checkers and source bytes remain untouched. Independent replay, formal
verification, external review and original goal completion remain false.
