# A whole-plane chart does not remove the complete-ruling obstruction

16 November 2026. SOURCE/ROUTE AUDIT ONLY. No new HC4 theorem is
registered, no external theorem is newly adopted, and no mathematical
checker is executed. HC4MYML0 and every retained checker remain frozen.

The proposed next strengthening of ML0 was to use the ACTUAL whole-plane
open to rule out a complete ruling's fibre without multiplicity-one
components. The classical higher-power Danielewski atlas already blocks
that geometric argument: its entire bad fibre can lie OUTSIDE a whole
plane chart, while its nowhere-zero volume is exact. This does not supply
the missing étale target map, a finite Keller normalization, or a
counterexample to HC4. The full original finite trace problem stays OPEN.

## 1. Retained evidence and source-reading limits

UADANCHART1 already gives, in every exponent and degree, the classical
Danielewski plane atlas, its transition functions, and the precise
compatible Keller-pair problem. UADANVOL1 already gives exact volume for
every higher exponent and excludes exponent1 by a nonzero period. The
surfaces, exact primitives and atlas are PRIOR work, not newly discovered
examples. Their checkers are NOT rerun.

The retained primary Basel2009 programme PDF, page2, contains Dubouloz's
talk abstract *Surfaces de Danielewski et Conjecture Jacobienne*. It
explicitly states that naive differential obstructions do not immediately
decide whether x^2*z=y^2-1 admits an étale morphism to A2. This is a
TALK ABSTRACT, not the full talk or a full proof. The two-page cached
programme was inspected; no slides or missing artifact are reconstructed.
The exact primary programme URL is
[Basel programme](https://dubouloz.perso.math.cnrs.fr/GdTs/GdTAutos/Prog-Basel-09-11.pdf).

PF2WRIGHT1 already identifies Wright's marked index-three ring, its source
plane open, omitted prime and volume zero. We read its complete canonical
note, not its historical open-case lists as current instructions. The
author-uploaded primary of Wright's 1997 paper was screened for its opening
setup and bundle terminology; its full exclusion theorem/proof was NOT
reproduced or newly adopted. The retained marked calculation is sufficient
for Section4 below. Published generic-degree-four PLANE exclusion remains
in force, and no quartic realization queue is reopened.

## 2. The existing atlas gives a complete ruling in every exponent

Let p(z) be square-free, of degree d>=2, and choose a root alpha. For
ANY integer m>=1 consider the known smooth affine surface

    W_m={u^m*v=p(z)},
    D_beta={u=0,z=beta}, beta a root of p.

Write p(z)=(z-alpha)*p_alpha(z). The retained plane chart is

    U_alpha=W_m minus union_(beta!=alpha) D_beta = A2_(u,w),
    z=alpha+u^m*w,       v=w*p_alpha(alpha+u^m*w).        (1)

On u!=0, w=(z-alpha)/u^m; near D_alpha it is v/p_alpha(z).
These regular descriptions agree on their overlap. They were already
proved as the classical atlas in UADANCHART1.

Define a map rho_alpha:W_m->P1, on the dense open u!=0, by

    rho_alpha=[u^m:z-alpha].                            (2)

The rational expression extends EVERYWHERE. On U_alpha it is [1:w].
Near D_beta for beta!=alpha, z-alpha is a unit, and it is

    [mu:1],       mu=u^m/(z-alpha).

Those charts cover W_m. Thus rho_alpha is a genuine COMPLETE-BASE
affine-line ruling, not a pencil with unresolved base points. For every
finite value w=c, its entire fibre is the A1 in the plane chart (1).
Its sole degenerate fibre is

    rho_alpha^-1(infinity)=sum_(beta!=alpha) m*D_beta.   (3)

Indeed u is a local parameter at each smooth D_beta and z-alpha is a
unit there, so mu vanishes to order m. No extra component or point can
occur, since the complementary finite-value chart is exactly U_alpha.

The morphism is AFFINE. Its finite-value chart has preimage U_alpha=A2.
The other base chart has preimage the complement of one fibre in the
smooth affine surface W_m. That fibre is an effective Cartier divisor;
its complement immersion is affine, locally a principal-open immersion.
Hence its total preimage is affine. This checks affineness over both
standard affine base charts, without inferring it from the word 'ruling'.

For m>=2, the known volume

    Omega_m=du wedge dz/u^m

is nowhere zero and EXACT, by retained UADANVOL1, and pulls back in (1)
to du wedge dw. Therefore ALL higher exponents give actual smooth affine
exact-volume surfaces with a whole-plane open and a complete ruling
whose fibre (3) has no multiplicity-one component. In fact that WHOLE
fibre lies OUTSIDE this plane open. Arbitrarily many components and
arbitrarily high common multiplicity occur. For m=1 the fibre is reduced,
but the volume is NOT exact, consistent with ML0.

These are classical surfaces. The explicit reading of their existing
atlas as (2) only records the hypothesis boundary relevant here. It
does not establish novelty, a rigid surface, or an actual Keller target.
These surfaces have additive actions, so they do NOT refute a separate
conjecture that the genuine canonical completion always admits one.

## 3. The missing condition remains the compatible étale TARGET

None of (1)--(3) provides a map e:W_m->A2 which is étale everywhere.
UADANCHART1 already writes that map as a system of polynomial Keller
pairs H_alpha on the ENTIRE plane atlas, satisfying

    H_alpha(u,w)=H_beta(u,w+(alpha-beta)/u^m).           (4)

The same target pair must glue across EVERY root component. This is
not accomplished by an exact primitive, a complete ruling, or the local
description of its nonreduced fibre. A hypothetical such e would restrict
to a nonautomorphic plane Keller map on EVERY atlas chart, by the retained
separated étale section argument. No compatible pair in (4) is constructed.

Consequently these controls do not identify a full finite normalization S,
its trace-dual Lambda, or its entire maximal étale open over an actual
polynomial target B. In particular one must not label W_m a verified
canonical completion simply because it passes the geometric hypotheses.
The all-degree trace-pairing requirement is still the original one.

## 4. The FINITE normalization is different from its maximal étale open

There is a second, independent scope hazard. In an actual finite Keller
normalization N=Spec S, the pulled-back source volume has divisor

    div(omega)=sum_E (e_E-1)*E

on omitted ramified primes. It is exact as d(P*dQ), but has ZEROS there.
ML0 requires the two-form to be NOWHERE ZERO. That hypothesis holds on
U=Spec B[Lambda], the maximal étale open, not automatically on all of N.
An application of ML0 to an affine-line bundle model of the FINITE N
therefore needs a new argument removing those zeros. No such argument
is present in the retained Wright reduction.

The known marked ring makes the distinction concrete:

    A=C[Y,XY,X^2*Y,X^3*Y+X] subset C[X,Y].               (5)

PF2WRIGHT1 already proves its full kernel, source plane open and single
omitted A1. Near that omitted line, its actual regular coordinates are

    X=u^-1,       Y=u^3*v-u^2,
    dX wedge dY=-u*du wedge dv.                         (6)

This form is regular and has a SIMPLE ZERO on u=0. Moreover the source
primitive itself extends there:

    -Y*dX=(u*v-1)*du.                                  (7)

Its exterior derivative is exactly the form in (6). These local expressions
agree with the polynomial source form and hence give a regular global
primitive on the known marked surface. No nowhere-zero assertion follows.
The model is not a constructed finite Keller target, and its former
quartic-realization wording is historical, superseded by the ledger's
published full-plane degree-four exclusion.

For perspective, the same zero is visible in the standard bundle
transition calculation. On a supplied A1-bundle over P1, suppose its
two charts have transition

    s=1/t,       w=t^n*z+b(t),       n>=2.

The source-chart form transforms as

    dt wedge dz=-s^(n-2)*ds wedge dw.                  (8)

It is regular at the other chart, but has order n-2 along its zero fibre.
Chart primitives -z*dt and s^(n-2)*w*ds differ by

    -b(t)*t^-n*dt.                                    (9)

The residue condition for exactness kills the coefficient of t^(n-1)
in b. For n=2 that is the whole cotangent-torsor class used by ML0.
For n>=3 it need not kill the whole bundle class. In (6), the transition
is precisely n=3, b=t: the residue is zero but (8) has a simple zero.
This calculation does not assert affineness or finite Keller realizability
for arbitrary formal transition parameters. The retained ring (5), rather
than such an inference, supplies the actual marked affine control.

## 5. Consequence for the research direction

The strengthening tested here is closed negatively: whole-plane inclusion
and exact nonvanishing volume do NOT rule out complete rulings with their
all-nonreduced fibre wholly outside the plane. The classical atlas gives
the controls in every exponent and degree. This does not settle whether
the genuine completion is rigid, and it does not construct an étale target.

Do not continue a multiplicity enumeration or repeat the Danielewski atlas
as a new construction. A useful further step must use the ACTUAL target
pair, the full finite normalization, or its trace multiplication. Neither
a ruling alone nor abstract surface recognition supplies that step.

FULL THEOREM TARGET: prove or break full HC4/JC2.
CURRENT LAST UNSOLVED IMPLICATION remains

    ACTUAL polynomial plane source
       ==> Tr(Lambda_L*Lambda_L) subset B.               OPEN

Full marked-energy integrality remains OPEN. No new theorem/import/cache,
mathematical checker invocation, registered replay or assurance upgrade.
Independent replay/formal verification/external review remain false.
Review and outreach remain DEFERRED.
