# Codifferent preimages under target derivations

17 November 2026. Finite normal complex surface covers, with an additional
uniform conclusion for an actual polynomial plane Keller source.

Let B=C[P,Q], F=C(P,Q), L/F finite, and S the full finite normal
integral closure of B in L. The ORDINARY trace-dual is

    Lambda={a in L : Tr_(L/F)(a*S) subset B}.

Write D_P,D_Q for the unique extensions to L of the coordinate
derivations of F. Then

    S={f in L : D_P(f) and D_Q(f) belong to Lambda}.       (1)

For an ACTUAL plane Keller inclusion B subset R=C[x,y], one generic
constant target derivation D works simultaneously for EVERY source
polynomial:

    S={f in R : D(f) belongs to Lambda}.                  (2)

The exceptional target directions are finite and depend only on the
images of the omitted divisors. They do not depend on f. The exact
multiplier ideal is computed in Section4.

Consequently, for generic target coordinate P* and source coordinates
(x,y), one has

    P*_x in Lambda  iff  y in S.                          (3)

In a source-finite frame, this single codifferent membership forces
the Keller inverse. Automatic membership is NOT proved. Lambda here
is the ordinary trace-dual, not an unidentified square-trace or
higher-moment lattice.

## 1. The complete height-one codifferent test

For each prime divisor E of Spec S let C be its target image and e_E
its normal ramification index. Characteristic zero gives

    a in Lambda  iff  v_E(a)>=1-e_E for EVERY E.           (4)

Here is the global justification. Normality of B tests regular target
traces at all height-one target primes. Localizing there gives a finite
semilocal Dedekind algebra. Its trace-dual localizations are the inverse
differents at all source DVRs. In characteristic zero the residue
extensions are separable and the tame different exponent is e_E-1.
Thus the local trace-dual is u^(1-e_E)S_E. The localization tests for
Tr(aS) then give exactly (4); there is no additional finite-point
condition. This is also the canonical coefficient identity used in
[the full-normalization boundary proof](CANONICAL_BOUNDARY_AND_STABLE_COMPLETE_INTERSECTION_OBSTRUCTION.md).

The tame exponent can be checked after an unramified residue extension
and a unit-root change, where t=u^e. The basis 1,u,...,u^(e-1) and
Tr(u^j)=0 unless e divides j identify the inverse trace lattice as
u^(1-e)S_E. The residue extension has perfect unramified trace pairing
and introduces no normal valuation. These operations preserve (4).

We need one differential calculation at those same primes. Say a target
derivation D is TRANSVERSE to C when D(H_C) is nonzero modulo H_C, for
a reduced local equation H_C of C. Then for q>0,

    v_E(f)=-q  implies  v_E(Df)=-q-e_E.                  (5)

For any target derivation regular on B, one always has

    v_E(f)>=0  implies  v_E(Df)>=1-e_E.                 (6)

At a general smooth point, away from the finite singular and residue
ramification loci, choose target normal and tangential coordinates
(t,z) and source coordinates (u,z) with t=u^e. This is the tame local
normal form, after the harmless unit-root change. A regular target
derivation has D(z) regular and D(u)=D(t)/(e*u^(e-1)).
Transversality makes D(t) a unit along E.

If f=u^-q*a(u,z), with a a unit, the leading term of Df is

    -(q/e)*a*D(t)*u^(-q-e).

Terms differentiating a are strictly later in u order. The nonzero
integer q and the unit D(t) prevent cancellation, proving (5).
For regular f the same expression for D(u) proves (6).
The proof uses generic divisor points to compute valuations, so neither
smooth branch crossings nor a surface classification is required.

The local ramification input is classical. [Stacks, Lemmas58.31.5--6](https://stacks.math.columbia.edu/tag/0EYH)
give the tame form t=u^e near a regular branch divisor after etale base
change. [Lemma49.12.6](https://stacks.math.columbia.edu/tag/0BWJ) identifies
the different with the Jacobian, and [Lemma49.8.2](https://stacks.math.columbia.edu/tag/0BW2)
identifies the trace-dual with the relative dualizing module. Applied at
generic divisor points in characteristic zero, these give the inverse
different u^(1-e) used in (4). The derivative-order calculation and the
uniform generic-direction conclusions below are deductions here; the
references do not assert automatic codifferent membership of a source
velocity.

## 2. Two target derivatives recover the full finite order

If f belongs to S, (6) and (4) show D_P f,D_Q f in Lambda.

Conversely, suppose f has a pole at some prime E. At the smooth generic
point of its image curve C, at least one of D_P,D_Q is transverse:
the two coordinate derivatives cannot both annihilate a local normal
parameter modulo that parameter. Formula (5) gives for that derivative

    v_E(D_i f)=-q-e_E<1-e_E,

contradicting (4). Thus both derivative conditions forbid every
height-one pole. Normality of S proves f in S. This establishes (1).

In particular the automatic statement is

    D_P(S)+D_Q(S) subset Lambda,                         (7)

rather than stability of S or Lambda under the target derivations.

## 3. One generic target derivative works on the actual whole source

Now let B subset R=C[x,y] be an actual polynomial plane Keller
inclusion. Keep its FULL finite normalization S. Its canonical source
map X=Spec R->Z=Spec S is an open immersion. There are finitely many
omitted prime divisors E_1,...,E_s; write Gamma_i for their image
curves in the target. No curve or finite point is discarded from Z.

Choose a nonzero constant combination

    D=alpha*D_P+beta*D_Q

transverse to every distinct Gamma_i. Such choices form the complement
of a finite set of projective directions. At the generic point of one
curve its tangent space has dimension one, so at most one projective
constant direction can be tangent there. Taking the union over the
finitely many images proves the assertion. This uses no theorem
excluding affine-line nonproperness components.

The target derivations preserve the ACTUAL source ring. If
J=Jac_(x,y)(P,Q) in C*, their expressions are

    D_P=J^-1*(Q_y*partial_x-Q_x*partial_y),
    D_Q=J^-1*(-P_y*partial_x+P_x*partial_y).              (8)

Thus Df is a source polynomial whenever f is. Also Lambda subset R:
at every retained source divisor the ramification index is one, so
(4) gives regularity there; normality of R gives source polynomiality.

For f in R, every possible pole on Z lies on an omitted E_i.
If f is not in S, formula (5) at such a pole contradicts Df in Lambda.
Conversely f in S gives Df in Lambda by (6). This proves (2),
simultaneously for all f in R.

The source condition in (2) is essential. Even for the identity cover
S=B=C[P,Q], D=D_Q and f=1/P give Df=0 in Lambda while f is not in S.
The general full-field statement therefore needs both derivatives as
in (1); it cannot be replaced by one fixed derivative on all of L.

There is also an exact analogue for the entire canonical etale
completion U=Spec A, A=B[Lambda]:

    A={f in R : Df in A}.                               (9)

Target derivations are regular on U, so they preserve A. If f in R
has a pole at a boundary divisor of U, that prime is unramified and
(5) gives a pole of order q+1 for Df. Normality of A excludes it.
This proves the reverse implication. No statement about surjectivity
of D on any quotient or a classification of its kernel follows.

## 4. The exact codifferent-clearing ideal

For a fixed f in R, let q_E=-v_E(f)>0 at its polar omitted primes.
For each distinct target image Gamma of such primes, choose its
reduced target equation H_Gamma. Define

    n_Gamma=max_(E over Gamma, q_E>0) ceil((q_E+1)/e_E),
    J_f=product_Gamma H_Gamma^(n_Gamma).                 (10)

An empty product is one. For the generic D chosen in Section3,

    {b in B : b*Df in Lambda}=J_f*B.                    (11)

Indeed at a polar E, (4)--(5) require

    e_E*v_Gamma(b)-q_E-e_E >= 1-e_E,

which is exactly v_Gamma(b)>=ceil((q_E+1)/e_E).
At primes where f is regular, (6) already places Df in Lambda and
multiplication by b causes no problem. Intersect all these tests.
Factoriality of B gives (11), including coincident final image curves.
The complete test (4) shows that no finite-point condition is hidden.

If s_1,...,s_N is ANY finite B-generating family for S, then (11) is
also the exact common denominator ideal of the finitely many traces

    Tr(s_i*Df), 1<=i<=N.                                (12)

This follows directly from the DEFINITION of the ordinary codifferent.
No global free basis or unmarked power-trace criterion is assumed.

For a strictly split local control t=u^e, f=u^-q, and j the remainder
of q modulo e, one sees the pole explicitly:

    Tr(u^j*D_t(f))=-q*t^(-floor(q/e)-1).                 (13)

Multiplication by t^n removes it exactly when
n>=floor(q/e)+1=ceil((q+1)/e). Formula (13) is a split local control,
not a replacement for residue-field trace duality in the general proof.

## 5. The trace defect has a derivative realization in the source

For actual Keller R and the same generic D, differentiation induces

    Lambda/S -> R/Lambda,
    [lambda] |-> [D(lambda)].                           (14)

It is well-defined by (7), and is B-linear: D(b*lambda)-b*D(lambda)
equals D(b)*lambda, which belongs to Lambda. Its kernel is zero by
(2). Thus EVERY nonzero ordinary trace-defect class differentiates to
a source polynomial lying outside the ordinary codifferent.

The image is a finite B-module isomorphic to Lambda/S. In particular
it retains the reduced-branch support proved in
[the trace-defect pairing theorem](REDUCED_BRANCH_PAIRING_AND_SIXTH_MOMENT_DEFECT.md).
The reduced branch equation H kills it directly: H*Lambda subset S,
so

    H*D(lambda)=D(H*lambda)-D(H)*lambda in Lambda.

This is not stability of Lambda under D. Rather,

    D(Lambda) subset Lambda  iff  Lambda=S.             (15)

In the actual plane setting the latter is the retained inverse
criterion. No differential IC/source-layer theorem is used in
(14), and no vanishing of its source module is asserted.

## 6. One first velocity detects the transverse source coordinate

Choose a target linear coordinate P* nonconstant on every omitted
image Gamma_i, and complete it to target coordinates (P*,Q*) with
source Jacobian one after constant rescaling. This excludes only
finitely many target directions. Then D_(Q*) is transverse to every
Gamma_i. In any fixed polynomial source coordinates (x,y),

    D_(Q*)(y)=P*_x.

Apply (2) to y. This proves (3), in EVERY source frame. The statement
does not require a generic source direction until the next implication.

Assume the source frame is SOURCE-FINITE, meaning R is finite over
B[y]. A sufficient condition is that some target polynomial has a
nonzero constant leading coefficient as a polynomial in x over C[y].
This holds in every generic constant source direction by its highest
homogeneous part.

If P*_x in Lambda, (3) gives y in S. Since x is integral over B[y],
it is integral over S; normality then gives x in S and R=S.
The original Keller map is now finite etale over A2 and hence has a
polynomial inverse. This is also immediate from the retained
full-normalization inverse criterion. The principal theorem needs no
Chau component theorem.

Under the separately proved ordinary L2 theorem in
[the Bergman trace-dual note](BERGMAN_TRACE_DUAL_AND_MULTIPLIER_DOMAIN_RIGIDITY.md),
Lambda is exactly the polynomial source L2 space for its specified
bounded target ball seeing the entire boundary. Only under that exact
identification may (3) be restated as a one-velocity L2 test.
No identification with a sixth-moment lattice is inferred here.

## 7. Exact support and remaining obligation

The checker verifies tame derivative orders with unit perturbations,
the two-derivative and source-polynomial scope barriers, inverse
Jacobian/coframe identities, exact split marked traces, clearing
exponents, reduced-branch derivative classes, and finite source
integrality arithmetic. Local branched and Laurent controls are
identified as controls, not actual noninvertible whole-plane Keller maps.

The proofs establish the unrestricted identities above. Their input is
ordinary trace duality, tame codimension-one differentiation, normal
height-one regularity and the actual source open immersion. The
analytic restatement additionally needs the separately proved LW0
theorem and its finite-quotient input. Neither it nor the algebraic
argument proves automatic codifferent membership of a first velocity,
full codifferent multiplication, T3 trace perfection, HC4 or JC2.

Replay:

    .venv/bin/python research/hessian/hc4-triangular-graphs/verify_codifferent_derivation_preimages.py

## 8. One source-finite velocity detects every omitted image

**Corollary SourceFiniteOneVelocityFullBoundaryTrace.**
Let \(B\subset R=\mathbb C[x,y]\) be an actual plane Keller inclusion
in a source-finite polynomial frame, so \(R\) is finite over \(B[y]\).
Use any polynomial target coordinate pair \(P^*,Q^*\) of \(B\).
At every omitted prime \(E\) of the full finite normalization one has
\(q_E=-v_E(y)>0\). For each distinct image \(\Gamma\) of an omitted
prime, let \(H_\Gamma\) be its reduced target equation and put
\[
 n_\Gamma=\max_{E\to\Gamma}\left(\left\lfloor q_E/e_E\right\rfloor+1\right).
\]
For any finite \(B\)-generating list \(s_1,\ldots,s_m\) of the
full normalization, the traces
\(\operatorname{Tr}(s_i P_x^*)\) have common denominator ideal
\(\prod_\Gamma H_\Gamma^{n_\Gamma}B\). A single constant-linear
combination of the markers can attain the maximal pole order
\(n_\Gamma\) simultaneously along every such curve.

If \(v_E(y)\ge0\), then \(B[y]\) lies in the source valuation ring
of \(E\). By source finiteness and integral closedness that ring
contains \(x\) as well. Thus \(x,y\) give a map of its valuation
spectrum to the original source plane. Composing with the source
open immersion into the full finite normalization agrees with
the original valuation map on the generic point. Separatedness
makes the maps agree on the closed point, contradicting that
\(E\) is omitted. Therefore \(q_E>0\) at every omitted prime.

An omitted image is a component of the nonproperness curve of
the actual map. Apply
[Chau's no-line component theorem](../../jacobian/affine-complements/sources/global-boundary-topology/chau-2003-two-remarks.pdf)
after the polynomial target change to \(P^*,Q^*\).
Every omitted image projects nonconstantly to the \(P^*\)-line,
since a vertical irreducible curve would be that affine line.
Its equation therefore has positive \(Q^*\)-degree, and its
\(Q^*\)-derivative is nonzero modulo the equation by the strict
degree drop in characteristic zero. This proves transversality
of \(D_{Q^*}\) at every omitted image for the selected frame.
The exact relation is \(P_x^*=J^*D_{Q^*}(y)\), with \(J^*\ne0\)
constant. The clearing calculation of Section 4 now gives the
asserted ideal, supported on the entire nonproperness curve,
including wholly unramified omissions.

At each \(\Gamma\), at least one trace in the generating list
has pole order \(n_\Gamma\), because otherwise its factor in the
common denominator could be lowered. The coefficient of pole
order \(n_\Gamma\) depends \(\mathbb C\)-linearly on the marker
coefficients and is a nonzero linear map to the residue field
of \(\Gamma\). Avoiding its kernel for the finitely many curves
gives one marker with all the asserted orders. Concretely, if
there are \(b\) distinct omitted images, among any
\(1+b(m-1)\) distinct integer choices for \(c\) some marker
\(\sum_{i=1}^m c^{i-1}s_i\) works. At each curve its residue
polynomial is nonzero of degree at most \(m-1\).

This determines the affected traces throughout the actual full
normalization. Eliminating the omitted primes, or proving that
the marked poles cannot arise for an actual Keller map, remains
open.
