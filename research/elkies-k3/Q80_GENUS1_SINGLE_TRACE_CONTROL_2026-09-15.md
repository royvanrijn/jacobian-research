# Exact rank18 on an arithmetic progression of genus-one Q80 covers

The height-eight pencil with trace **P1** (zero-based source index) has
exact function-field rank18 at parameter0 and at every integer parameter

```
lambda = 175361103335727332945 * n,    n in Z.
```

These quadratic covers of the original parameter line have genus1. None
has two new independent directions, at any section height. Rational points
on their covering curves are not certified or needed for this exclusion.
This is a scoped control for the [open correlated-gain objective](CORRELATED_QUADRATIC_GAINS_2026-09-12.md),
not its solution or an exclusion of the entire pencil.

## Literal construction

Use the actual17-section [direct11952 model](R17_NORM12_ORBIT11952_DIRECT_FIBRATION_2026-09-03.md).
Write P1=(Nx/h²,Ny/h³), with h monic of degree2, and choose

```
M0 = -Ny/Nx mod h²,       degree M0 < 4,
M_lambda = M0 + lambda*h²,
q_lambda = (M_lambda^4 - 6*M_lambda²*Nx - 8*M_lambda*Ny
            - 3*Nx² - 4*A*h^4)/h^6.
```

The regular chord theorem gives a polynomial quartic q_lambda. Its explicit
five coefficient polynomials in lambda are recorded in the
[complete pencil note](Q80_COMPLETE_GENUS_ONE_BRANCH_INJECTIVITY_2026-09-13.md#projective-branch-maps).
The [input](../artifacts/generated-results/elkies-k3-q80-genus1-single-trace-control-v1/input.json)
freezes the generic source hash and literal h,Nx,Ny,M0,q0. Parameter0 was
chosen before testing branch arithmetic. No exceptional specialized curve
or point, search-bank target, or reconstruction of missing artifacts enters.

On C_lambda: w²=q_lambda, the new section is

```
x = (M_lambda²-Nx)/(2h²) + (h/2)*w,
y = (M_lambda*(M_lambda²-3Nx)-2Ny)/(2h³) + (M_lambda/2)*w.
```

At lambda0 its polynomial coefficient degrees are4,2,6,3. The checker
verifies both coefficients of the curve identity and verifies that its
branch point doubles to P1 by exact arithmetic modulo q0.

## Branch certificate and exact rank

The quartic q0 reduces irreducibly modulo5. The checker proves this with
Frobenius remainders, rather than trusting a factorization label. At the
following simple roots of q0, all17 parent sections reduce on smooth fibres:

| prime | roots used |
|---|---|
|53|23|
|73|33|
|101|0|
|109|62|
|113|106|
|127|14,60,77,119|
|131|35|
|137|86|
|139|71|

Complete enumeration of each finite elliptic group and its subgroup of
doubles gives a combined matrix of rank16 over F2. Its kernel is the line
spanned by e1. Hence, at the single closed quartic branch alpha,

```
ker(M/2M -> E_alpha(Q(alpha))/2E_alpha(Q(alpha))) <= span(e1).
```

At prime23, t=0 is a simple branch root and the smooth fibre has order33.
Good reduction injects branch-field2-torsion into this group, so there is
none. By the retained [unramified Kummer and trace-norm theorem](Q80_GOOD_BRANCH_CODE_AND_TRACE_NORMS_2026-09-14.md),
the actual anti-invariant branch image is zero and the gain is at most1.
The displayed section differs from its conjugate. All base-changed fibres
are irreducible and chi=4; the height formula excludes nonzero torsion.
Thus its anti-invariant difference supplies one direction, proving exact
rank18 and equality of the halving kernel with span(e1).

## Why the entire progression has the same rank

All five coefficient polynomials of q_lambda are integral at

```
5,23,53,73,101,109,113,127,131,137,139.
```

The checker verifies this coefficientwise through valuation bounds. For
lambda divisible by their product N above, q_lambda reduces to q0 at all
these primes. Irreducibility modulo5 ensures that q_lambda is an irreducible
squarefree degree4 polynomial over Q. Its branch field has every simple
residue root listed above by Hensel lifting. The same parent-section reductions
therefore give the same rank16 matrix and the same no-2-torsion witness.
The prime23 witness also prevents q_lambda from dividing the parent
discriminant, so all branch fibres are smooth. The regular chord construction
still supplies a non-invariant section. The rank argument consequently
applies to every lambda in N*Z, without a parameter or section-height cutoff.
This is an infinite set of pencil parameters; pairwise isomorphism classes
of the underlying abstract genus-one curves are not being counted.

The [checker](scripts/verify_q80_genus1_single_trace_control.py) and
[result](../artifacts/generated-results/elkies-k3-q80-genus1-single-trace-control-v1/result.json)
retain the exact identities, residue matrices, group orders and congruence
modulus. Replay takes below0.1 seconds, under20 CPU seconds and1GiB:

```
.venv/bin/python research/elkies-k3/scripts/verify_q80_genus1_single_trace_control.py
```

The parent rank/saturation, global unramified theorem, regular pencil theorem
and height argument are inherited or written mathematics, not formal
verification. Only the finite arithmetic and indicated polynomial identities
are replayed here. A two-gain cover outside this progression remains possible;
a larger local kernel would only be a necessary gate, not global solubility.
