# Exact volume obstructs reduced component selections in every complete ruling

HC4MYML0. 16 November 2026. All degrees and all complete-type affine-line
rulings of the full canonical étale completion; no ruling is assumed to exist.

Let X be a smooth complex AFFINE surface with a regular nowhere-zero
EXACT two-form omega. Then X contains NO open subset which is a Zariski
locally trivial A1-bundle over P1. Consequently, for EVERY affine
A1-fibration pi:X->P1, at least one scheme-theoretic fibre has NO
MULTIPLICITY-ONE component. This is stronger than merely requiring one
nonreduced component. It does NOT say that the gcd of that fibre's component
multiplicities exceeds1, and it does not rule out the fibration altogether.

The first assertion has a direct two-chart proof below. The passage from
fibre components to a locally trivial bundle uses explicitly credited
published affine-line-fibration structure. For an actual plane Keller
map it applies to the ENTIRE canonical étale completion U=Spec B[Lambda].
It does not prove an additive action on U, trace regularity, or HC4/JC2.

## 1. Local-structure input and reading boundary

We adopt one standard affine-line-fibration input package, as stated in
Adrien Dubouloz,
[A1-cylinders over smooth A1-fibered surfaces, arXiv:2002.12838v1](https://arxiv.org/html/2002.12838v1),
Section1, Proposition3 and the paragraph following Definition4:

* For an affine A1-fibration on a smooth complex surface over a smooth
  curve, each degenerate fibre is a disjoint union of affine lines with
  positive scheme multiplicities. The source credits Miyanishi, Lemma1.4.2.
* If every scheme-theoretic fibre is irreducible and reduced, the fibration
  is a ZARISKI locally trivial A1-bundle. Proposition3 cites its earlier
  published curve-base result [19,Lemma1.3].

The complete relevant statements, Definitions1/4 and Examples2/5, the
entire Theorem9 proof SKETCH and the Section1.3 torsor discussion were read.
The source explicitly says that some results there are given without full
detailed proofs. Its cited local-structure foundations are NOT independently
reproduced here; this is an adopted standard input, not a consequence of
the checker. We use neither cancellation nor the general algebraic-space
classification to prove the two-chart obstruction. The primary PDF is
retained with its byte hash. No whole-paper or independent proof replay
is claimed.

The stronger nonaffine-total-space theorem of Kambayashi--Wright was also
screened in its author-uploaded primary text; the publisher PDF was blocked.
It is NOT a second adopted input. Our selected subfibration below is an
AFFINE morphism, so the explicit Proposition3 above suffices.

Retained KU0/1 and LA0 provide the full completion, and LM0/1 provide the
exact-volume/additive-rigidity boundary. UADANVOL1's period exclusion of
ordinary Danielewski surfaces is PRIOR work, not a discovery made here.
The cotangent-torsor calculation below is a standard geometric mechanism;
no priority claim is made for it. Its use here supplies a uniform condition
on all complete-type rulings of the actual canonical finite-algebra model.

## 2. An open A1-bundle over P1 is incompatible with exact affine volume

Suppose V subset X is open and rho:V->P1 is a Zariski locally trivial
A1-bundle. Use the standard cover A1_t and A1_s with s=1/t. On either
chart the bundle is trivial: its affine-linear transition first defines
a line bundle, trivial on A1, and its translation class is in H1(A1,O)=0.
Thus the chart inverse images are A2, with coordinates (t,z) and (s,w).

The restriction of omega is nowhere zero on each A2. Its coefficient
relative to dt wedge dz, or ds wedge dw, is a unit of a polynomial ring,
hence a nonzero CONSTANT. Rescale fibre coordinates so that

    omega=dt wedge dz=ds wedge dw.

On the overlap, the affine-line-bundle transition has the form

    s=1/t,       w=a(t)*z+b(t),
    a in C[t,t^-1]*, b in C[t,t^-1].

The equality of two-forms forces

    a(t)=-t^2,
    w=-t^2*z+b(t).                                      (1)

Translations z'=z+f(t), w'=w+g(s) change b by

    b'(t)=b(t)+t^2*f(t)+g(1/t).                         (2)

All Laurent monomials of exponent>=2 can be removed by f, and all of
exponent<=0 by g. The only remaining coefficient is the exponent1 term:

    w=-t^2*z+c*t,       c in C.                         (3)

This is the usual one-dimensional cotangent-torsor class, derived here
without a surface classification.

Define regular chart primitives

    theta_0=-z*dt,       theta_infty=-w*ds.

Their differentials are omega. By (3), on the overlap

    theta_infty-theta_0=c*dt/t.                         (4)

Since omega=d theta on X, it is exact on V. On each of its A2 charts,
theta-theta_0 and theta-theta_infty are closed polynomial one-forms,
hence polynomial differentials by the algebraic Poincaré lemma. An
elementary proof of that lemma integrates each positive homogeneous
component with the Euler contraction; division by its positive degree
is valid over C. Therefore (4) must be the differential of a regular
function on the overlap Gm times A1.

But c*dt/t is exact in C[t,t^-1,z] ONLY if c=0. Indeed, a primitive
would have zero z derivative and so belong to C[t,t^-1]. The derivative
of a Laurent polynomial has ZERO coefficient of t^-1: a monomial t^n
differentiates to n*t^(n-1), and the sole possible n=0 contributes0.
This covers EVERY exponent, not a finite Laurent window. Thus c=0.

For c=0, (3) is a line bundle transition and its two zero sections glue
to a copy of P1 in V. This is a nonconstant complete curve in the affine
surface X, impossible: every regular function on P1 is constant, so
every morphism P1->X=Spec O(X) is constant. This proves the assertion.

Both sides of the contradiction matter. If c!=0 the volume is NOT exact;
if c=0 the total space has a complete curve and cannot be open in an
affine surface. We do not replace this argument by a false assertion
that every affine-line bundle over P1 is a line bundle.

## 3. Choose a multiplicity-one component in each fibre: a global contradiction

Let pi:X->P1 be an AFFINE A1-fibration. The adopted local structure says
that its finitely many degenerate fibres have disjoint A1 components.
Suppose EVERY fibre has a component of multiplicity1. Choose ONE such
component in each degenerate fibre and delete ALL other components of
those fibres. Let V be the resulting open subset of X. No generic fibre
is changed, and no base point is lost. Every scheme-theoretic fibre of

    rho=pi|V:V->P1

is now ONE REDUCED A1. In particular this is not just normalization of a
nonreduced fibre, nor a selection which omits an exceptional base value.

We check the hypotheses needed for the adopted bundle theorem. The deleted
finite union of curves is an effective Cartier divisor because X is smooth.
Its complement immersion V->X is affine: locally on X the divisor is one
equation and its complement is a principal affine open. Hence rho is an
affine finite-type morphism, being a composite of affine morphisms. It is
surjective. It is flat over the smooth curve P1: on a local base DVR its
source coordinate modules are torsion-free, since V is integral and rho
is dominant. The generic fibre is unchanged A1 over C(P1), and each closed
fibre is the chosen reduced A1. Proposition3 therefore makes rho a ZARISKI
locally trivial A1-bundle. Section2 gives a contradiction.

We have proved the all-fibre selection obstruction

    for EVERY complete ruling pi:X->P1,
    SOME fibre pi^-1(b)=sum_i m_i*C_i has ALL m_i>=2.    (5)

In particular an everywhere reduced complete-type ruling is impossible.
A mixed fibre with both a multiplicity1 and a higher-multiplicity component
does not by itself meet (5). Conversely (5) allows multiplicities2 and3
in one fibre: their gcd is1. There is NO conclusion that a whole fibre
is a multiple of another divisor, or that pi1(X) is nontrivial. The component
selection is the condition used, not an orbifold gcd calculation.

## 4. Application to the actual canonical étale completion

Let K:A2->A2 be any actual complex plane Keller map, with the retained
full normalization S and trace-dual Lambda over B=C[P,Q]. KU0/1 gives

    U=Spec B[Lambda],

the ENTIRE maximal étale open of Spec S, including every omitted unramified
prime. LA0/LM0 give its smooth affineness and the actual plane open.
Its target functions give the exact nowhere-zero form

    omega=dP wedge dQ=d(P*dQ).

Thus Section2 applies to U in EVERY generic sheet degree, and EVERY
complete-type affine-line ruling of U obeys (5). Its complete base curve,
if any, is P1: the rational surface U contains a dense plane, and a
nonconstant generic source line maps to the base, forcing that smooth
complete curve to be rational. This can also be read from the inclusion
of its one-variable function field in C(x,y).

This is a constraint on the actual finite-algebra completion, not on an
arbitrarily substituted rational surface. It does NOT prove that a ruling
exists, that every ruling has affine base, or that U admits a nonzero LND.
In particular LM0's ADDITIVELY RIGID branch remains open. LM1 says that
noncyclic generic deck groups require that branch; the new result only
forces any complete-type ruling there to have a fibre as in (5).
The generic deck group is not the full monodromy group.

Even U=A2, or absence of all added UNRAMIFIED primes, does not by itself
remove omitted RAMIFIED primes of the finite normalization. No inference
from (5) to Tr(Lambda*Lambda) subset B is made. PF128/129's finite
intermediate restrictions remain separate; no intermediate is assumed
to be a plane and no solvable-monodromy inversion theorem is used.

## 5. Sharp hypothesis controls and remaining goal

* The nonaffine cotangent bundle T*P1 has its exact canonical two-form and
  its zero section P1. It shows why affineness of the ambient surface is
  necessary in Section2. It is not an HC4 source.
* The smooth affine quadric uv=z^2-1 is the standard NONTRIVIAL A1-bundle
  over P1 in Dubouloz Example2. Its nowhere-zero volume is NOT exact, as
  already proved by the retained UADANVOL1 period calculation. It shows
  why trivial canonical bundle alone cannot replace exactness.
* Retained higher-power Danielewski surfaces x^m*y=p(z), m>=2 with p
  square-free, have exact nonvanishing volume and an AFFINE-base ruling
  whose fibres are reduced. They are not excluded by (5), whose base
  must be COMPLETE. No étale target pair is supplied by these controls.
* The integer profile(2,3) meets the no-multiplicity-one requirement but
  has gcd1. This is an arithmetic boundary check, not an asserted realized
  plane Keller fibre or a counterexample.

The checker transcribes the transition determinant, primitive difference,
Laurent coboundary normalization and characteristic-zero residue identity,
plus provenance and hypothesis controls. The all-exponent and all-fibre
arguments are written proofs; a finite list of checks does not prove them.
Earlier registered mathematical checkers are NOT replayed. The adopted
local-structure package has its stated source/reproof limits. No independent
replay, formal verification or external review is claimed.

FULL THEOREM TARGET remains proof or counterexample for full HC4/JC2.
CURRENT LAST UNSOLVED IMPLICATION remains

    ACTUAL polynomial plane source
       ==> Tr(Lambda_L*Lambda_L) subset B.                OPEN

Full marked-energy integrality remains OPEN. The rigid completion and
the finite trace-multiplication gap are not replaced by a new component
classification queue. Review and outreach remain DEFERRED.
