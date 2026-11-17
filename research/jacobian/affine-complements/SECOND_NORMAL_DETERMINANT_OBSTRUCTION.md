# A determinant obstruction at the second normal coefficient

14 November 2026. PF2RULED25 excludes a nonempty Zariski open set of
first coordinates in every one of the sixteen marked noncontact
pole(6,9) cusp-lift families. It allows the Jacobian constant to vary:
this is not merely the failure of a chosen normalization kappa=1.

Every surviving first coordinate must satisfy an explicit nonzero
polynomial Delta of degree at most69 in its31 affine parameters.
This necessary equation also covers the cases where the coefficient
matrix drops rank. On a specified rank69 chart the full second-normal
test reduces to24 rational compatibility equations in P alone; kappa
is uniquely determined, and Q, if it exists, has exactly17 free
directions. The nonzero-kappa solution locus remains unresolved.
No entire marked candidate, pole pair or quartic packet is excluded.

## The full fixed-first-coordinate system

Keep [PF2RULED24's actual global first-cusp-jet spaces](FULL_NEXT_NORMAL_LIFTS.md)
and its exact marked coefficient fields. For either K=Q[y]/(P8), write

    P(a)=P*+sum_(i=1)^31 a_i*U_i,
    Q(b)=Q*+sum_(j=1)^85 b_j*V_j.

The references and full kernel bases use exactly the inherited cusp
inverse circuits, not a sampled subfamily. Symbols H^2,H^3 and the
complete equation L(F)=h*F1-lambda*F0'=0 are fixed. In particular
P0 has degree24 for every a; the lower V5 terms cannot change it.
Target translations have not been divided out.

Put J1=[u]Jac_uv and J2=[u^2]Jac_uv. The existing exact formulas are

    J1(P,Q)=2*(P2*Q0'-P0'*Q2)+P1*Q1'-P1'*Q1,

    J2(P,Q)=3*(P3*Q0'-P0'*Q3)
            +2*(P2*Q1'-P1'*Q2)+P1*Q2'-P2'*Q1.

Both are linear in Q for fixed P. J1 has degree<=48 throughout the
whole first-cusp-jet family. Once J1 is constant, PF2RULED24's intrinsic
source-Jacobian bound gives deg J2<=43. Therefore the complete system

    J1(P,Q)=kappa, J2(P,Q)=0                             (1)

has49+44=93 coefficient rows. The unlisted J2 coefficients above43
are forced to vanish after the first49 rows are imposed; no extra
assumption on them is made.

Let B(P) be the93-by85 matrix of the two displayed operators on V_j,
let e have value1 only in the row for [v^0]J1, and let r(P) be minus
the operator applied to Q*. The system is

    M(P)*(b,kappa)=r(P), M(P)=[B(P),-e].                (2)

Each column of B and r is affine linear in a. Only the kappa column
is constant. Requiring kappa!=0 remains part of the actual Keller
problem; allowing every kappa in (2) makes the exclusion below stronger.

## A universal rank bound, before evaluating any parameters

There are15 independent sections of V8 whose coefficients F0,...,F3
all vanish. The complete rational Cox jet matrix has precisely this
kernel, as replayed by the checker. Every one is in ker L and is
annihilated by both operators for every P.

There are two more kernel directions: the functions1 and P. The
filtration inclusion V6 subset V8 puts P in the same Q-direction
space, and both have zero Jacobian with P. They are independent
modulo the15 hidden sections, because their restrictions1 and P0
are independent polynomials, with deg P0=24. Thus universally

    rank B(P)<=85-17=68,
    rank M(P)<=69.                                    (3)

This is a global bound on the entire31-dimensional P-space, including
its rank-drop strata. It is not an upper bound inferred from modular
or numerical samples.

## Nonzero minors give a global necessary polynomial

Evaluate at a=(1,...,1), the exact circuit point already used in
PF2RULED24. The new residue matrices at prime1009 give

| Root multiplicities | rank B | rank M | rank [M,r] | rank69 minor D | augmented minor Delta |
| --- | --- | --- | --- | --- | --- |
| (2,7) | 68 | 69 | 70 | 280 | 145 |
| (4,5) | 68 | 69 | 70 | 326 | 885 |

The last two columns are nonzero residues modulo1009. All divisions
in the inherited cusp circuits have residue units. The same local-ring
argument as PF2RULED24 consequently proves that these minors, evaluated
at this exact characteristic-zero circuit point, are nonzero in K.
The two irreducible fields then cover every complex embedding, hence
all sixteen marked candidates.

The certificate specifies69 columns I and69 rows J giving

    N(P)=M(P)_(J,I), D(P)=det N(P).

It also specifies one additional row t and defines the polynomial

    Delta(a)=det [M(P(a)),r(P(a))]_(J union {t}, I union {r}). (4)

Here the order is exactly the retained row/column order. Formula (4),
the complete Cox matrices and the inherited inverse circuits specify
an exact polynomial over K without expanding its potentially large
coefficients. A numerical interpolation or an unpinned symbolic
elimination is not used.

The minor includes the kappa column, whose only nonzero entry is -1.
Expanding along it gives deg Delta<=69. The nonzero residues145/885
show that Delta is not the zero polynomial. If (2) has any solution,
its augmented rank is at most69 by (3), so

    Delta(a)=0                                         (5)

is necessary for every solution, with any kappa. This implication
requires no division by D. It therefore retains every rank-drop
solution as well. Since Delta is nonzero at a=(1,...,1) in every
embedding, its complement is a nonempty excluded Zariski open set
in each of the sixteen complex P-families.

The proof has two distinct inputs: exact residue-unit nonvanishing
and the written universal rank bound (3). A nonzero augmented minor
alone, without a valid global rank bound, would not justify (5).

## Complete local elimination of Q and kappa

On the chart D(P)!=0, (3) gives rank M=69 and rank B=68. In particular
the constant column e is independent of the Q columns. Solve the69
pivot rows with the remaining coefficient controls set to zero:

    x_I=N(P)^-1*r_J(P).

The selected columns span the full column space of M. Hence (2) is
equivalent on this chart to the24 remaining-row equations

    r_l(P)-M_(l,I)(P)*N(P)^-1*r_J(P)=0, l not in J.    (6)

These are24 rational conditions on31 P parameters, not assertions
of independence or a dimension formula for their zero locus. Their
denominators are powers of D; clearing them requires retaining D!=0.
In the retained order, the extra row t is row69, corresponding to
[v^20]J2. Its residual is nonzero at the circuit witness and satisfies

    Delta=D*(that residual).

If (6) is consistent, kappa is uniquely determined. Indeed a difference
of solutions with nonzero change in kappa would place e in the column
space of B, contradicting rank B=68 and rank M=69. The homogeneous
kernel in the86 unknowns has dimension17 and is exactly

    (15 sections invisible through u^3) plus span{1,P},

with zero kappa component. Thus the full Q fiber, if nonempty, is one
particular section plus those17 directions. The two visible directions
are precisely a translation and a linear target shear. The hidden15
must return at later normal coefficients.

The chart test still needs kappa(P)!=0 and the actual singleton-image
conditions. Its complement D=0 is not excluded: cover the appropriate
rank strata with other minors and their compatibility tests. The global
necessary polynomial (5) applies to all of them, but is not a sufficient
test or a replacement for those strata.

## A reversible target normalization leaves30 P parameters

The derivative matrices and Delta are unchanged by translating P. Set
P0(s)=0. Its leading coefficient in degree24 is1 throughout the marked
family, so a shear Q<-Q-a*P with a=[v^24]Q0 removes that Q coefficient
without division. Then translate Q to impose Q0(s)=0. These three target
automorphisms preserve the two leading symbols, L, every Jacobian
coefficient, the target algebra, field degree and surface finiteness.
They do not set or change kappa.

Thus one may solve (6) in30 normalized P parameters, still with24
rational compatibility equations and the nonzero-kappa requirement.
If a Q fiber exists on the rank69 chart, the two visible directions
span{1,P} are uniquely fixed by the Q normalizations. Precisely the15
hidden directions remain. This is a reversible normalization of the
full problem, not deletion of valid target presentations or an assertion
of independence of the24 equations. Rank-drop strata still require
their own complete compatibility analysis.

## Why zero-kappa points must not count as realizations

There are exact points satisfying all relaxed equations with kappa=0.
Use the same T in V3 and the surjective map L:V2->K[v]_(<=7) of
PF2RULED22 to construct

    S=T+R2, L(S)=0,
    P=S^2, Q=S^3.

A recorded8-by8 residue-unit minor supplies the exact inverse circuit
for R2. The lower differences P-T^2 and Q-T^3 lie in V5 and V8,
respectively, so this pair belongs to the complete affine families.
The universal identity d(S^2) wedge d(S^3)=0 proves all its Jacobian
coefficients vanish in characteristic zero, not just after reduction.

Consequently Delta vanishes at these power points, whereas it is
nonzero at a=(1,...,1). The necessary polynomial is nonconstant and
its zero locus is nonempty. The relaxed full determinantal compatibility
locus is also nonempty. None of this establishes a point with
kappa!=0: these maps are algebraically dependent and are not Keller
or finite degree-four realizations. Saturating by the nonzero-kappa
condition remains essential in solving the actual problem.

## Reproduction, remaining scope and assurance

Run from the repository root:

```sh
sage -python research/jacobian/affine-complements/verify_second_normal_obstruction.py
.venv/bin/python research/jacobian/affine-complements/audit_pole_integration.py
```

The [checker](verify_second_normal_obstruction.py),
[determinant circuits and exact certificates](second_normal_obstruction_certificates.json)
and [receipt](second_normal_obstruction_receipt.json) retain the complete
rational hidden kernel, both marked matrix constructions, all pivot
choices, the nonzero residues, the24 residual equations and the
zero-kappa power controls. Universal rank bounds, characteristic-zero
unit lifting and completeness of the local elimination are written
proofs using the already proved section spaces. Shared exact utilities
are used; no formal verification, independent whole-proof implementation,
external review, new imported theorem or novelty priority is claimed.

The nonzero-kappa second-normal locus, every rank-drop case, the later
coefficients u^3 throughu^9, contact cases, other height-nine pairs,
higher heights, all-height pole-four coordinates, field degree4,
global surface finiteness and actual compatible embedding remain open.
The first open pairs remain (5,9),(6,9),(7,9),(8,9). The full original
goal remains active.
