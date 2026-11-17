# Quotient points of the full Keller normalization and canonical L2 extension

Let K:X=A2_C -> Y=A2_C be an actual polynomial Keller map, in any
polynomial or generic degree. Let S be the full finite normalization of
B=C[Y] in C(X), and Z=Spec(S). Every singular point of Z is a finite
quotient surface singularity. Consequently every regular reflexive
canonical form on Z is locally square integrable, including at all
omitted finite points. This supplies the point-singularity input to the
global L2 trace-dual theorem. It does not prove multiplication closure,
trace perfection, or JC2.

## 1. The exact primary theorem and its application

The primary theorem is Flenner--Zaidenberg,
[*Rational curves and rational singularities*, math/0109221v1](https://arxiv.org/pdf/math/0109221v1),
Section2.1, Theorem2.1 and its proof (PDF pp13--14; published Math. Z.
244 (2003), 549--575, DOI 10.1007/s00209-003-0497-z). Their result says
that a singular point of a normal affine surface is a quotient point
if general surface points lie on polynomial curves avoiding that point.
It requires no C*-action, chosen fibration, or assumed rationality of the
singularity. The authors prove it using Lemma2.2, Proposition2.3 and the
classical cylinder theorem attributed to Miyanishi in Theorem2.4.
Those arguments and their affineness discussion were read. The theorem
is adopted with that classical dependency, rather than claimed as an
independent reproduction of the original cylinder theorem.

The canonical lift j:X->Z is an open immersion by the normalization and
Zariski--Main argument in KN0 Section1. The surface Z is normal, finite
over affine Y and therefore affine; j(X) is dense. All singular points
lie in Z minus j(X), because X is smooth.

Fix any such singular point p. A general point q lies in j(X). Choose
a source affine line L=A1 through j^-1(q). Its composite with j is a
nonconstant polynomial morphism A1->Z whose image avoids p. We check
that its image is a CLOSED polynomial curve in Z: affineness makes this
point stronger than mere closedness inside the source open.

Let A be the image of S in C[t] under restriction to L. Some a(t) in A
is nonconstant, since j restricts injectively to L. If d=deg a, the
monic polynomial obtained by dividing a(T)-a(t) by its nonzero leading
coefficient makes t integral over C[a(t)]. Thus C[t] is finite over
C[a(t)], and the generators 1,t,...,t^(d-1) also generate it as an
A-module. Therefore Spec C[t]->Spec A is finite and surjective. Its
scheme-theoretic image Spec A is closed in Z, so the actual line image
has no additional boundary limit points in Z. In particular it still
avoids p. The family of source lines gives the exact polynomial curves
required through general points. Flenner--Zaidenberg Theorem2.1 now
shows that p is a finite quotient point. Since p was arbitrary, the
conclusion holds at every singularity of the FULL normalization.

The cylinder X=A1 times A1 would also meet the hypotheses of the
stronger Miyanishi statement restated as Flenner--Zaidenberg Theorem2.4,
which has a cyclic conclusion. Cyclicity is unnecessary here. We use
the primary authors quotient theorem and preserve that weaker exact
scope, without asserting that the original 1981 proof was inspected.

## 2. Why arbitrary finite quotient points suffice for LW0

Let alpha be a regular reflexive canonical form near a quotient point
p of Z. Use a sufficiently small analytic presentation of the germ as
C2/G with finite small G, so the quotient cover is etale off the origin.
Its pullback alpha is an ordinary holomorphic top form on a punctured
smooth cover chart. Hartogs extends its coefficient across the origin.
The extended coefficient is bounded on a smaller chart, so the canonical
L2 integral upstairs is finite there. Outside the measure-zero branch
point, change of variables identifies the downstairs canonical integral
with 1/|G| times the upstairs integral. Thus alpha is locally L2 at p.
This argument uses no cyclic representation or positivity of a bilinear
trace form.

At smooth points alpha is already holomorphic. A finite normalization
has compact preimage of any compact target set, so finitely many smaller
charts prove the global canonical L2 bound used in LW0 for a bounded
target ball. The generic-divisor calculation there still detects exactly
the codifferent pole orders. Therefore quotient points supply its missing
converse at ALL finite singular points; arbitrary higher Lp thresholds
are not inferred directly from this local top-form argument. The even
integer exponents use the already written polynomial-power reduction.

## 3. Keep the differential source layer separate

The former KQ0/PF2BATT1 citation bundled three claims: the plane open
immersion, its quotient-point boundary geometry, and a semisimple
curve-supported differential quotient. KN0 proves the open immersion;
Section1 above supplies the quotient-point claim using the adopted
primary theorem. Neither fact proves the differential IC-socle and
source-layer assertions invoked by KU2 in Section4 of
CANONICAL_TRACE_DUAL_ETALE_RECONSTRUCTION. That differential input remains
unverified unless its own complete proof is supplied.

This point theorem is enough for LW0 and its analytic even-moment
successors. Those results should not import the entire stronger KQ0
package as a premise. The precise primary bytes and assurance limits are
retained in [sources/quotient-point-context.json](sources/quotient-point-context.json).
No new checker, formal verification, independent external review or
novelty claim accompanies this application of the classical theorem.
