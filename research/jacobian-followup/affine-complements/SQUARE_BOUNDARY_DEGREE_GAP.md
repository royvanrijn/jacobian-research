# A boundary degree gap around the square/cube locus

14 November 2026. PF2RULED26 excludes a family inside the unresolved
nonzero-kappa problem of [PF2RULED25](SECOND_NORMAL_DETERMINANT_OBSTRUCTION.md).
At every one of the sixteen marked noncontact pole(6,9) candidates,
suppose the first coordinate has boundary restriction S0^2+c for an
actual pole-three section S with symbol H and the complete cusp identity.
If the companion differs from S0^3 on the boundary by a polynomial of
degree at most7, its first-normal Jacobian cannot be a nonzero constant.

Every allowed lower normal correction is retained. This is not only an
exclusion of the exact powers S^2,S^3, nor a generic-rank argument. It
covers the full four-parameter family of such actual S, including its
special parameter values. The existence of this global pole-three root
is an explicit hypothesis, not a conclusion for every first coordinate
on the determinant zero locus. No whole marked candidate is excluded.

## Precise hypotheses in the actual ring

Use the specified Wright inclusion and omitted-prime chart of
[PF2WRIGHT1](WRIGHT_MARKED_MODEL.md). Write F=F0+uF1+u^2F2+... and put

    h=v-s, lambda=3h/4+delta, delta!=0,
    L(F)=h*F1-lambda*F0'.

Keep one of the sixteen fixed common nonics H. Let P in V6 and Q in
V9 have symbols H^2,H^3 and satisfy the complete equations L(P)=L(Q)=0.
Assume there is S in V3 with symbol H and L(S)=0 such that

    P0=S0^2+c, R(v)=Q0-S0^3, deg R<=7.                 (1)

Then J1=[u]Jac_uv(P,Q) cannot be any nonzero constant. Equivalently,
any pair passing the first-normal test under the root hypothesis in
(1) must have deg(Q0-S0^3)>=8.

The theorem requires only the first-normal Jacobian equation, not
surface finiteness, field degree4 or any higher Jacobian coefficient.
It applies to actual compatible finite quartic candidates as a necessary
boundary test, but does not assume or prove their existence.

The root condition is exhaustive over actual S: after fixing the
symbol H, the complete V2 correction space has dimension12 and its
cusp operator has rank8. Thus all S are one exact reference plus its
four-dimensional kernel. The checker constructs that entire family
in each irreducible marked degree-eight field, without selecting a
numerical embedding.

## Retaining the normal corrections

Translate P to remove c and translate Q to make R(s)=0. These target
translations change neither the source Jacobian nor the target algebra.
Evaluation of L at s gives S0'(s)=Q0'(s)=0, so R'(s)=0. Put

    r=S0'/h, L_R=R'/h.

Because P and S^2 have the same symbol, F=P-S^2 lies in V5. Its
coefficients F0,F1 vanish. The complete first-two-jet kernel of
[PF2RULED22](FULL_CUSP_FIRST_JET_LIFTS.md) therefore gives

    A=F2, deg A<=10.                                  (2)

No higher coefficient of F is set to zero.

For any polynomial R with R'(s)=0, the actual global function

    T_R=R(D)-A1*lambda(D)*R'(D)/(D-s), A1=XY,

has jets

    (T_R)0=R, (T_R)1=lambda*L_R, (T_R)2=-v*lambda*L_R.

Here A1=-u+u^2v exactly. These formulas involve a polynomial quotient;
there is no untracked rational pole at s. For degree7, T_R can have
pole9. We do not infer that it is in V8. Subtract it only to write

    G=Q-S^3-T_R, G0=G1=0, G2=B.

The proof permits B to be any polynomial; it imposes no extra degree
bound on this normal correction. Thus every actual allowed Q correction
is retained, even when the chosen auxiliary T_R has higher pole than
Q-S^3. All later coefficients of G are also unrestricted by this test.

## The exact first-normal equation

For any two cusp lifts, write

    J1(U,V)=2*(U2*V0'-U0'*V2)+U1*V1'-U1'*V1.

Let C_R=J1(S,T_R) and absorb the freely occurring term by setting
Btilde=B-(3/2)*S0*A. Direct use of the product rule gives

    J1(P,Q)=2*S0*C_R+2*h*A*L_R-4*h*S0*r*Btilde.        (3)

This identity is independent of the discarded higher normal orders.
The proof can also use the linear expression

    C_R=lambda^2*r*L_R'
        +[2h*(S2+v*lambda*r)-lambda^2*r']*L_R.

Suppose (3) is a constant kappa. If R is constant, its reduction modulo
S0 immediately gives kappa=0. If deg R<=2, then L_R is constant and

    kappa=2h*A*L_R modulo S0.

The right side has degree at most11, whereas S0 has degree12.
The congruence is an equality of polynomials. Evaluation at s forces
kappa=0. A nonconstant linear R is already impossible because R'(s)=0.
This treats every lower-degree case, without a division by a perturbation
coefficient.

Now let d=deg R>=3. Divide R by its nonzero highest coefficient to
normalize that coefficient to1. This is only a division in equation
(3): divide kappa and Btilde by the same scalar. S, P and their leading
symbols are unchanged. Write R=h^d+sum_(j=2)^(d-1)c_j*h^j.

Reduction of (3) modulo S0 gives

    2h*A*L_R=kappa+c(v)*S0, deg c<=d-3.               (4)

Substitution into (3) gives

    C_R-2*S0'*Btilde=-c(v)/2.

Divide C_R by S0' and write C_R=S0'*U+M, deg M<=10.
Since d<=8 in the computations retained below, a necessary condition is

    deg M<=d-3, c=-2M,
    kappa=2*M(s)*S0(s),
    kappa-2*M*S0 is divisible by h*L_R.                (5)

The division by S0' has fixed nonzero leading coefficient, because
S0 has the fixed degree12 monic marked leading term up to sign. The
division by h*L_R also has fixed nonzero leading coefficient d. These
are polynomial operations on the entire S-parameter space. Neither
critical-root discriminants nor variable rank minors are inverted.

Conditions (5) are used only as necessary conditions. No unsupported
global realization is inferred from a solution of this reduced system.

## Exact certificates close degrees3 through7

For each of the two fields K=Q[y]/(P8), the producer uses four variables
for the complete S family and d-2 variables for the normalized boundary
perturbation. It computes C_R and its complete remainder M. The high
coefficient ideal is

    I_d=([v^j]M : j=d-2,...,10).

Every calculation is in characteristic zero. The full exclusion is:

| d | High-remainder stage | Further divisibility stage |
| --- | --- | --- |
| 3,4 | exact unit ideal | not needed |
| 5,6,7 | six linear Groebner generators | exact unit ideal after (5) |

For d=5,6,7, reduce S0,L_R,M modulo the exact high-remainder ideal.
The remaining test is polynomial divisibility of
2*M(s)*S0(s)-2*M*S0 by h*L_R. Adding its remainder coefficients to the
high ideal gives1. This conclusion does not use a nonzero-kappa
saturation: for these degrees no constant first-normal Jacobian is
possible in this family, even the zero constant.

The [certificates](square_boundary_gap_certificates.json) retain explicit
Bezout coefficients whose exact sums equal1 for all ten field/degree
cases. They also retain identities in both directions between the high
generators and their Groebner basis before the second-stage reduction.
The checker replays all these polynomial equalities over K. It does not
rely on a reported Groebner flag or a successful reduction modulo a prime.
Irreducibility and the exact identities transport the result to all
eight complex embeddings of each field, covering the sixteen marked
points.

This is a complete finite algebraic obstruction for the stated boundary
family. Section coverage, the normal-Jacobian reduction and why this
family is only part of the global determinant locus are written proofs.
There is no formal verification, independent whole-proof implementation,
external review, new imported theorem or novelty priority claim.

## Degree eight is retained as an unresolved checkpoint

At d=8 the high-remainder ideal is not the unit ideal. In both fields
it has a dimension5 locus and nonlinear Groebner generators. The exact
high basis, its membership identities, the reduced divisibility equations
and the condition invk*kappa-1 are retained in the compressed
[degree-eight checkpoint](square_boundary_degree8_checkpoint.json.gz).
These are necessary equations; their solvability and actual global
realization have not been proved.

One direct elimination of the saturated degree-eight equations in the
(2,7) sector was observed live at105 seconds and then interrupted. Its
process returned KeyboardInterrupt. It supplied no theorem and is not
running. The [exploratory input](sources/square_boundary_degree8_exploration.py.txt)
is preserved in the checkpoint and source tree. The default checker
rebuilds the checkpoint but never restarts that full elimination.

Do not infer emptiness from the stopped attempt, or nonzero-kappa
solutions from the nonempty high-remainder locus. The general P locus
of PF2RULED25, the cases lacking an actual global pole-three root, all
rank-drop strata, later coefficients u^3 throughu^9, contact cases,
other height-nine pairs, higher heights, all-height pole-four, field
degree4, surface finiteness and actual compatible embedding remain open.
The first open pairs are still (5,9),(6,9),(7,9),(8,9). The original
full goal remains active.

## Replay

Run from the repository root:

```sh
sage -python research/jacobian-followup/affine-complements/verify_square_boundary_gap.py
.venv/bin/python research/jacobian-followup/affine-complements/audit_pole_integration.py
```

The [checker](verify_square_boundary_gap.py) and
[receipt](square_boundary_gap_receipt.json) pin all proof inputs, the exact
certificates and the compressed checkpoint. Python and Sage versions
are recorded. The replay proves the degree-seven exclusion and only
constructs the degree-eight input; it does not launch its unsolved search.
