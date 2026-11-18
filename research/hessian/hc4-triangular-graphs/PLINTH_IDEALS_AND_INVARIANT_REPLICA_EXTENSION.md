# Plinth ideals and invariant replica extension

18 November 2026. Actual complex polynomial plane Keller maps, with the
full finite normalization and its entire canonical etale completion
kept distinct.

Let B=C[P,Q] subset R=C[x,y], F=Frac(B), L=Frac(R), with nonzero
constant source Jacobian. Let S be the FULL finite normalization of
B in L, Lambda its ORDINARY trace-dual, and

    A=B[Lambda], U=Spec A.

Use the proved common source coordinate y and its translation
delta=partial_x on A. There are two exact assertions:

1. The plinth ideal delta(A) intersect C[y] and the ideal of all
   positive-x leading coefficients of A are the SAME explicit
   principal ideal g(y)C[y]. If U is proper, g has at least double
   zeros at every added fibre.
2. A nonzero invariant replica q(y)delta extends to S, or stabilizes
   Lambda, if and only if y already belongs to S.

The first assertion concerns the canonical etale ring A. The second
concerns the full finite ring S, including all ramified primes outside
U. Neither is an automatic landing theorem. In particular the common
frame is never source-finite when U is proper.

## 1. Common-coordinate boundary data

The [exact-volume common-pencil theorem](EXACT_VOLUME_BOUNDARY_MULTIPLICITIES_AND_COMMON_SOURCE_PENCIL.md)
and [whole root charts](COMMON_COORDINATE_FORCES_SMOOTH_BOUNDARY_AND_WHOLE_ROOT_CHARTS.md)
supply polynomial source coordinates (x,y) with

    y in A, delta(A) subset A, ker(delta|A)=C[y].

Write D_1,...,D_s for the added UNRAMIFIED boundary curves of U minus X.
They are whole smooth affine lines. For each i,

    y|D_i=lambda_i,
    v_Di(y-lambda_i)=r_i>=1,
    v_Di(x)=-p_i, p_i>=r_i+1.                            (1)

The values p_i are computed in this common complementary coordinate x.
To check this last point, any individually adapted coordinate x_i
has the form a_i*x+h_i(y), with a_i nonzero constant. The second
summand is regular at D_i and cannot cancel the negative value of x_i.
Thus the individual whole-chart pole order is also p_i in (1).

No omitted ramified prime is included among these D_i. Such primes
remain in Z=Spec S outside U. Normality of A tests regularity of a
source polynomial using the D_i alone.

For each distinct lambda among the lambda_i put

    c_lambda=max_(i:lambda_i=lambda) ceil(p_i/r_i),
    g(y)=product_lambda (y-lambda)^(c_lambda).            (2)

If there is no added boundary, g=1 and A=R. Otherwise every
c_lambda>=2 by the strict exact-volume inequality in (1).

## 2. Exact plinth and leading-coefficient ideals

Define the plinth ideal in the explicit invariant ring by

    pl(delta)=delta(A) intersect C[y].

Then

    pl(delta)=g(y)C[y].                                  (3)

Indeed if delta(F)=q(y), then source polynomiality gives

    F=q(y)*x+h(y), h in C[y].

At D_i, h is regular. A pole of qx cannot cancel against it, so F
is regular there exactly when

    r_i*ord_(y=lambda_i)(q)>=p_i.

Conversely these inequalities make qx regular at every D_i, and it
is already polynomial on X. Normality puts qx in A, with delta(qx)=q.
Grouping the inequalities by lambda proves (3), including every
closed point of U. An arbitrary scalar normalization of g gives the
same ideal.

For every F in A of positive x-degree M,

    delta^M(F)=M!*LC_x(F) in pl(delta).

The preceding derivative delta^(M-1)(F) still lies in A. Therefore

    LC_x(F) belongs to g(y)C[y].                         (4)

Conversely each q in gC[y] occurs as the leading coefficient of
the degree-one element qx in A. Thus the set of nonzero positive-x
leading coefficients, together with zero, is exactly gC[y].

In particular every positive-x leading coefficient of an element
of S or B is divisible by g, since B subset S subset A. Equality
of their own leading-coefficient ideals with gC[y] is NOT asserted.
The plinth generator is not silently placed in B, S, or Lambda.

Localizing gives

    A[g^-1]=C[x,y,g^-1].                                 (5)

For g*x and y are in A, so x is obtained after that localization.
This is the precise invariant denominator appearing in the local
slice. When U is proper, the equation delta(F)=1 has NO solution
in A although delta is locally nilpotent and has a complete source
translation flow. Solving delta(F)=q(y) is possible exactly when
g divides q. Finite Taylor subtraction therefore does not supply
unconditional integration inside A.

## 3. The proper common frame is never source-finite

Suppose R were finite over B[y]. Because y belongs to A,

    B[y] subset A subset R.

A finite B[y]-generating family for R also generates it over A.
Thus R is finite over A. Lying-over makes Spec R->Spec A surjective.
But this is the original open immersion X->U. A surjective open
immersion is an isomorphism, so A=R and U=X.

Consequently

    U proper  implies  R is NOT finite over B[y].        (6)

This is stronger than failure of a chosen sufficient monicity test.
It explains why the common-pencil translation cannot simply be
inserted into a theorem requiring a generic source-finite frame.
The exact leading-coefficient ideal (4) also shows directly that,
in the proper case, no positive-x target polynomial has nonzero
constant leading coefficient.

When U=X, g=1 and every source translation exists on A=R. This says
nothing about whether S=R: ramified primes outside U can remain.

## 4. A replica extends to the full finite ring exactly when y is integral

For the GIVEN source coordinates, the following are equivalent:

    (a) y belongs to S;
    (b) some q in C[y], q!=0, has q(y)*delta(S) subset S;
    (c) some q in C[y], q!=0, has q(y)*delta(Lambda) subset Lambda. (7)

Here the q in (b) and (c) can be chosen the same. Every successful
replica is locally nilpotent on S, since it is locally nilpotent on
the source ring R and its coefficient is killed by delta.

First assume (b). Choose a generic constant target linear coordinate
P* nonconstant on every omitted image. The exact pole calculation of
[codifferent preimages](CODIFFERENT_PREIMAGES_UNDER_TARGET_DERIVATIONS.md)
says that if y has a pole of order a>0 at an omitted prime E of
index e, then

    v_E(delta(P*))=-a-e.

For a nonzero polynomial q(y) of degree d,

    v_E(q(y))=-d*a,

because its highest power is uniquely dominant at a pole of y.
Therefore

    v_E(q(y)*delta(P*))=-(d+1)*a-e<0.                    (8)

But P* belongs to B subset S, contradicting (b). Thus y has no
boundary pole. It is already regular on X, so normality gives y in S.
For the common coordinate all possible poles in this argument are
at RAMIFIED primes outside U, since y is regular on U.

Conversely assume y in S. Localize at the nonzero elements of C[y].
At least one of P,Q has positive x-degree over C(y), and its leading
coefficient is a unit in that field. Thus x is integral over
C(y)[P,Q], and normality of S tensor_(C[y]) C(y) gives

    S tensor_(C[y]) C(y)=C(y)[x].                         (9)

A denominator h(y)!=0 can consequently be chosen with

    S[h^-1]=R[h^-1].

For example clear a denominator in an expression for x in (9).
Choose finite C-algebra generators s_j of S. Their source derivatives
delta(s_j) belong to R[h^-1]=S[h^-1]. One power h^N clears these
finitely many denominators. The replica h(y)^N*delta then preserves S.
This proves (a) implies (b), with no source-finiteness assumption
on the original frame.

This sufficiency step is the classical extension of a translation from
a principal cylinder. In the proof of
[Kishimoto--Prokhorov--Zaidenberg, *Group actions on affine cones*,
arXiv:0905.4647v2, Proposition3.5, printed p16](https://arxiv.org/pdf/0905.4647v2),
a translation on D(f)=Z times A1 is multiplied by f^N to preserve the
finitely generated ambient coordinate ring. Here (9) and denominator
clearing identify D(h) in Spec(S) with
Spec(C[y,h^-1]) times A1_x, and h is invariant under partial_x, so the
displayed construction is exactly that proof. Its statement and both
directions of the proof were inspected. The converse pole obstruction
and the equivalence with preservation of the ordinary codifferent in(7)
remain separate parts of the present argument; the cylinder theorem
supplies no automatic integrality of y.

To relate the trace-dual, recall the exact multiplier algebra

    (Lambda:Lambda)=S.                                   (10)

One inclusion is S-stability of Lambda. If multiplication by f
preserves Lambda, the determinant trick on this finite faithful
B-module makes f integral over B; full integral closedness gives
f in S. No freeness of a chosen global basis is required.

If D=q(y)delta preserves Lambda, then for s in S and lambda in Lambda,

    (D s)*lambda=D(s*lambda)-s*D(lambda) in Lambda.

Equation (10) gives D(s) in S, proving (c) implies (b).

If D preserves S, it is a regular derivation at every height-one
point of the normal Z. Moreover it preserves the SOURCE volume:

    Lie_D(dx wedge dy)=0,

since delta(q)=0. The ordinary codifferent coefficients identify
regular reflexive canonical forms as Lambda*(dx wedge dy).
Lie derivative of regular canonical forms remains regular at every
smooth height-one point, and

    Lie_D(lambda*Omega)=D(lambda)*Omega.

The complete height-one test for Lambda gives D(Lambda) subset Lambda.
This proves (b) implies (c). It uses only ordinary finite trace/canonical
duality, not the unverified differential source-layer theorem.

## 5. What finite Taylor and the canonical module do not force

The implication in (7) is an extension criterion, not an existence
theorem for q. Formula (8) shows why invariant polynomial multipliers
cannot repair a pole of y at a ramified prime: every positive degree
makes the pole WORSE. Target-polynomial multipliers may clear source
derivations into S, but the retained target multiplier theorem does
not make those rescalings locally finite source flows.

The [ramified fold control](RAMIFIED_FOLD_VALUATIONS_AND_COMMON_COORDINATE_POLES.md)
contains the explicit smooth affine surface with marked source

    a=x, b=xy, c=xy^2-1, d=xy^3-y.

Its whole noncontracting plane open has boundary E=A1. On the other
chart y=u^-1, x=u^2+u^3*t, one has

    Omega=dx wedge dy=u*du wedge dt,
    theta=x*dy=-(1+u*t)*du.

Both Omega and theta are globally regular, with dtheta=Omega, while
y has a pole at E. The linked source proves the affine presentation,
chart covers, smoothness and complete marked source; they are not
reproved here.

Its canonical coefficient module relative to this volume is

    M=O(E)=S0+S0*y,   S0=O(surface),   S0[M]=C[x,y].

This follows from the two charts: away from E the generator 1
suffices, and at E the function y=u^-1 generates the allowed simple
pole. Nevertheless the source LND delta=partial_x has

    delta(b)=y in M minus S0,
    delta(c)=y^2 outside M,
    delta^2(b)=delta^2(c)=0.

Thus a finite Taylor chain can end in a nonintegral invariant or
leave a finite canonical lattice. This control has NO supplied
finite Keller plane target; M is a canonical coefficient module,
NOT claimed to be an ordinary trace-dual for such a target.
It does not refute an argument genuinely using that extra trace
compatibility. It blocks deriving the landing from the canonical
completion, exact volume, source noncontraction and local nilpotence
alone.

## 6. Evidence and remaining scope

The checker verifies the valuation inequalities giving the exact
plinth exponents, leading-coefficient extraction, source-frame
localizations, invariant multiplier pole growth, regular coordinate
denominator clearing, the canonical-module commutator, and the
linked affine control's Taylor obstruction. The unrestricted proofs
are written above; finite controls do not prove a valuation or surface
classification.

The input for (1)--(6) is the proved actual common-pencil and
whole-component chart theorems. The extension criterion (7) itself
holds in any given polynomial source frame. It does not need common
coordinates, Chau's component theorem or a source LND classification.

No source velocity is proved automatically integral or automatically
in the ordinary codifferent. The source-finite one-velocity inverse
criterion still needs its own hypotheses. Full trace perfection,
automatic cubic trace regularity, HC4 and JC2 remain open.

Replay:

    .venv/bin/python research/hessian/hc4-triangular-graphs/verify_plinth_replica_extension.py
