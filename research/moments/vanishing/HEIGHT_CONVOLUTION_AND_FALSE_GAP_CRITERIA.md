# Height convolution and false finite-gap criteria for GVC

16 November 2026. Internal written all-order proofs. No independent replay,
formal verification, external review or novelty is claimed. The
[checker](verify_gvc_false_gaps.py),
[receipt](gvc_false_gaps_verification.json) and
[source audit](gvc_false_gaps_source_audit.json) separate exact controls from
the universal arguments.

There are explicit pure GVC failures in EIGHT variables with arbitrarily
long DOUBLING blocks of zero first-shift outputs, while the marker P^2
is nonzero at EVERY positive index. This disproves a generic GVC extension
of the symmetric-inverse doubling-gap test. It is not a counterexample to
formal eventuality-to-linearity: these pairs fail individual eventuality.
The current minimum dimensions, orders, degrees and Laplacian frontiers
are unchanged. No large expansion or search is used.

## 1. Exact height calculus on independent variable blocks

Work over a characteristic-zero FIELD K. Let T_i be formal constant-
coefficient operators with zero constant symbol coefficient, and let P_i
be nonzero polynomials on disjoint sets of variables. Write

    h_i(m)=max{k>=0:T_i^k P_i^m!=0}.

Zero operators are allowed. All heights are finite by degree lowering.
For P=product_i P_i, the following identities are exact:

    h_(sum T_i,P)(m)=sum_i h_i(m),
    h_(product T_i,P)(m)=min_i h_i(m),
    h_(T^e,P)(m)=floor(h_(T,P)(m)/e), e>=1.               (1)

**GVCHEIGHTCONV20261116.** For the sum formula, expand by the multinomial
theorem. At power H=sum h_i(m), the ONLY potentially nonzero allocation
has exponent h_i(m) in every block. Its multinomial scalar and its product
of nonzero polynomials are nonzero over K. Every higher power vanishes.
Lower powers cannot vanish, since applying more copies would then kill
the nonzero H-th output. The product formula factors the action blockwise;
the power formula follows from nilpotence itself.

These are polynomial-output identities and nonvanishing statements; no
evaluation at output zero is used. General nonreduced rings are excluded
from the height assertion: products of nonzero coefficients can vanish.
Pointwise sums of heights do not justify adding their limsup slopes without
checking synchronization, and no such assertion is made here.

## 2. Pure safe pairs already disprove subadditivity of nilpotence indices

For n>=2 and integer r>n, let

    T=sum_(i=1)^n partial_(x_i)^r,      P=product_(i=1)^n x_i.

The expansion of T^k P^m has a nonzero term precisely when k copies can
be distributed with at most floor(m/r) in each coordinate. Distinct
allocations give distinct monomials, and their nonzero coefficients are
positive integers. Thus

    h_m=n floor(m/r),       a_m=h_m+1.                    (2)

The pair is pure and safe, with sufficient real cutoff 1/(r-n), by total
degree: rm>degQ+nm implies T^m(QP^m)=0. This coefficient is sharp. At
m=r*l take Q=x_1^((r-n)m); the allocation ((r-n+1)l,l,...,l) gives a
nonzero constant at the equality boundary. The least integer cutoff is1.

**GVCNILSUBADDFAIL20261116.** At indices1 and r-1,

    a_r-a_1-a_(r-1)=n-1.                                 (3)

In particular n=2,r=3 gives a_3=3>a_1+a_2=2, despite full purity and
safety. A proposed dimension-free quasi-subadditivity error is also
impossible, since (3) is unbounded as n grows. This does NOT disprove a
bound depending on the actual phase/symbol or on dimension/degree; such
a bound has not been established. A bounded-error inequality of the right
direction would need a separate structural proof before use in the general
apolar-gap question. This family is a proof-method control, not a new
programme of higher-order safe examples.

## 3. The retained failure seed has an exact height

Use [DVG1](../../foundations/extended-geometry/DVORSKY_GVC5_COUNTEREXAMPLE.md)
on variables a,b,c,d,t:

    F=(t+c)(ad+bt),    B=partial_t(partial_a partial_d-partial_b partial_c).

Its retained all-order theorem says B^mF^m=0 for every m>=1, while
B^m(cF^m) is nonzero for every m>=2. The latter marker has degree1.
The [formal commutator lemma](FORMAL_SHIFTED_POWERS_AND_APOLAR_NILPOTENCE_GAPS.md)
therefore implies B^(m-1)F^m!=0 for m>=2: if it were zero, the marked
output would vanish. At m=1, BF=0 but F!=0. Consequently

    h_(B,F)(m)=m-1 for EVERY m>=1.                         (4)

This uses the retained all-order seed proof, not extrapolation from its
small checker. No new independently reviewed Dvorsky proof is claimed.

## 4. A fixed eight-variable family defeats arbitrary doubling gaps

Fix ANY integer M>=1 and set r=2M+3. On three new variables u,v,w put

    G_r=u v w^(r-2),
    A_r=partial_u^r+partial_v^r+partial_w^r,
    P_r=F G_r,         Lambda_r=(B+A_r)^2.                 (5)

This is an explicit finite constant-coefficient operator and a homogeneous
phase in eight variables. The phase has degree r+3; the operator has
maximum order2r and minimum order6. For r>=5 its terms have orders
6,r+3,2r, so it is not an ordinary Laplacian or an operator of one
homogeneous order. No smaller known endpoint is improved.

Write m=qr+s with 0<=s<r. Direct monomial allocation gives

    h_(A_r,G_r)(m)=2 floor(m/r)+floor((r-2)m/r)
                  =m-ceil(2s/r).                         (6)

Use (1) and (4). The height for Lambda_r on P_r^m is

    h_m=floor((2m-1-ceil(2s/r))/2)
       =m-1  if 0<=s<=M+1,
        m-2  if M+2<=s<=2M+2.                            (7)

All coefficients and exponents in (5)-(7) are explicit. In the second
case m>=M+2>=3, so no negative height occurs. In particular h_m<=m-1
at every m, proving full PURE vanishing.

**GVCFALSEDOUBLINGGAP20261116.** For the pair (5),

    Lambda_r^k P_r^(k+1)=0 for EVERY M+1<=k<=2M,           (8)
    Lambda_r^k P_r^(k+2)!=0 for EVERY k>=1.                (9)

For (8), m=k+1 lies in [M+2,2M+1], so (7) gives h_m=m-2<k.
For (9), put m=k+2. Both cases of (7) give h_m>=m-2=k, so that
output is nonzero. Thus Q=P_r^2 is one fixed marker bad at EVERY index.
There are also infinitely many nonzero first-shift outputs: at k=jr-1,
m=jr has h_m=m-1=k. The first-shift support is exactly periodic:
it vanishes iff (k+1) mod r is in [M+2,2M+2].

The same formula also handles the exact inverse-component indexing.
Put L=M+1>=2. For EVERY L+1<=j<=2L, one has

    Lambda_r^(j-1)P_r^j=0.                               (10)

Thus neither convention for shifting the index rescues the generic
doubling-gap test. Equation (10) has precisely the component-index window
used in the symmetric inverse theorem, without asserting that these
general-operator expressions are genuine inverse components.

The deficit b_m=m-h_m is always1 or2. It does not tend to infinity.
Hence this is NOT an eventual-but-nonlinear pair and does not refute
OP-GVC-FORMAL-LINEAR. Its role is an all-order counterexample to the
proposed finite doubling-gap certificate for general constant operators.
Even all degree<=1 markers vanish at every time m in [M+2,2M+2], by
the commutator lemma and b_m=2, while the pair remains unsafe.

## 5. Why this does not contradict the genuine symmetric gap theorem

Wright's [math/0511214v2](https://arxiv.org/html/math/0511214v2), Section2.3,
Theorems2.5-2.6, concerns the homogeneous components of a symmetric formal
inverse. A quadratic recurrence for those components propagates a doubling
block of zeros. The precise hypotheses and recurrence were read in the
primary HTML, together with the version history. Arbitrary GVC shifted
powers do not automatically satisfy that recurrence. Equations (5)-(9)
disprove its unrestricted transplantation; they do not challenge Wright's
theorem or its valid Laplacian/Hessian-nilpotent applications.

No theorem from that paper is imported to prove the construction. The new
proof uses retained DVG1, the formal commutator lemma, and (1). The
uniform-Mathieu/co-integral literature discovered in this turn supplies no
verified general linear-gap theorem here; search snippets are not proof
dependencies. No conclusion is drawn from unrelated notions of nilpotence.

## 6. Verification boundary and the next decisive question

The checker verifies small independent-block convolution/power identities,
the sparse auxiliary height formula by direct differentiation, low seed
heights, the exact subadditivity failure and sharp safe boundary, and
integer instances of the all-M periodic formula and doubling window. It
does not expand P_r^m or Lambda_r^m in eight variables at large parameters.
The all-order conclusions follow from (1),(4),(6), not a numerical prefix.

The actual formal apolar-gap implication remains OPEN. Neither index
subadditivity nor an unqualified symmetric doubling-gap rule is available
as a shortcut. A bounded-error structural theorem, if proposed, must be
proved for the actual formal-symbol/polynomial-phase family; (3) sets a
necessary dimensional allowance. No arbitrary abstract deficit sequence
has been realized. Both full exceptional Haar MZ cases remain open at
actual complement logarithms. All protected frontiers and radial transfer
remain unchanged; closure-first and deferred review/outreach persist.
